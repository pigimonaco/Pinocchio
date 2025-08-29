/**
 * Handles the import and export of the initial density.
 *
 * The initial density is the mass overdensity δ = (ρ - <ρ>)/<ρ> at redshift 0.
 * Reading the initial density bypasses what otherwise happens in `GenIC.c`.
 *
 * With the functions provided in this module, the density can be written to
 * or read from a file. It is a headerless binary file containing a series of
 * 64-bit floating-point numbers representing the mass overdensity at each grid
 * point. It uses C-like index order, with the z-index varying fastest (memory
 * stride of one number, i.e. neighboring points on the z-axis are stored next
 * to each other in memory) and the x-index varying slowest (memory stride of
 * one y–z plane).
 *
 * In Python, the file could be read like so:
 * ```python
 * δ = numpy.fromfile(file_name, dtype='float64').reshape( (nx, ny, nz) )
 * ```
 *
 * The file uses the same memory order as `rvector_fft`, the internal input
 * buffer for the forward Fourier transform. The file must be C-contiguous so
 * that we can efficiently read and write it in parallel across MPI tasks.
 * The internal memory order is reflected by the `COORD_TO_INDEX()` macro,
 * which elsewhere in the code base handles the transformation from (x, y, z)
 * coordinates to those linear C-like memory indices. We don't use that macro
 * in the code below, but do rely on that same memory order.
 */

#include "pinocchio.h"
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>


// We have no use for tagging MPI messages, so always use the same tag.
int const MPI_TAG_0 = 0;


/**
 * Creates a custom MPI data type for the distributed array `rvector_fft`.
 *
 * Pass in the `grid` parameters, then find `array_type` and `chunk_size` as
 * output.
 *
 * The chunk size is the number of grid points of `rvector_fft` that this very
 * MPI task stores. In the unlikely event that there are more MPI tasks than
 * grid points along the x-axis, it may be zero. In that case we return
 * `MPI_DOUBLE` for `array_type`. It acts as a dummy, since going forward we
 * will have to make the same MPI calls in each task so the job as a whole
 * would not hang.
 *
 * Make sure you free the array type after use by calling
 * `free_rvector_array_type(&array_type, chunk_size)`.
*/
void create_rvector_array_type(
    grid_data const grid,
    MPI_Datatype*   array_type,
    int*            chunk_size)
{
    // Shorten variables for clarity.
    ptrdiff_t const nx = grid.GSglobal[_x_];
    ptrdiff_t const ny = grid.GSglobal[_y_];
    ptrdiff_t const nz = grid.GSglobal[_z_];
    ptrdiff_t const dx = grid.GSlocal[_x_];
    ptrdiff_t const dy = grid.GSlocal[_y_];
    ptrdiff_t const dz = grid.GSlocal[_z_];
    ptrdiff_t const x0 = grid.GSstart[_x_];
    ptrdiff_t const y0 = grid.GSstart[_y_];
    ptrdiff_t const z0 = grid.GSstart[_z_];

    // We implicitly cast from `ptrdiff_t` to `int` for the array sizes before
    // passing them to `MPI_Type_creates_subarray()` as `mpicc` would otherwise
    // complain.
    int const global_shape[3] = {nx, ny, nz};
    int const local_shape[3]  = {dx, dy, dz};
    int const start_point[3]  = {x0, y0, z0};

    // Create array type unless this very tasks will store no data.
    *chunk_size = dx*dy*dz;
    if (*chunk_size == 0)
        *array_type = MPI_DOUBLE;
    else {
        MPI_Type_create_subarray(
            3, global_shape, local_shape, start_point,
            MPI_ORDER_C, MPI_DOUBLE, array_type
        );
        MPI_Type_commit(array_type);
    }
}


/**
 * Frees the custom MPI data type for `rvector_fft` after use.
 */
void free_rvector_array_type(MPI_Datatype* array_type, int chunk_size)
{
    if (chunk_size > 0)
        MPI_Type_free(array_type);
}


/**
 * Reads the initial mass overdensity from a raw binary file.
 *
 * The file must be a headerless binary file of 64-bit floating-point numbers
 * representing the initial mass overdensity array as a function of `(x, y, z)`
 * in row-major (C-like) index order. After reading the file, we
 * Fourier-transform the data and store it in `kdensity`, which would typically
 * point to the global array that Pinocchio uses as the basis for the
 * simulation. It is assumed that memory has already been allocated for that
 * array.
 *
 * The `grid_index` is typically 0, referrring to the first of (in principle)
 * multiple `grids`.
 *
 * Returns 0 on success, a non-zero error code otherwise.
 */
int read_initial_density(
    double**         kdensity,
    char*      const file_name,
    int        const grid_index,
    grid_data* const grids)
{
    // Get MPI rank (task number) of this process and total number of tasks.
    int rank, tasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);

    // The root task (a.k.a. main task) is the one at MPI rank 0.
    int const root = 0;

    // Read in parallel from the input file.
    MPI_File file;
    MPI_Datatype array_type;
    int chunk_size;
    create_rvector_array_type(grids[grid_index], &array_type, &chunk_size);
    int error = MPI_File_open(
        MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file
    );
    if (error != MPI_SUCCESS && rank == root) {
        printf("Error: Could not open file \"%s\".\n", file_name);
        return 1;
    }
    MPI_File_set_view(
        file, 0, MPI_DOUBLE, array_type, "native", MPI_INFO_NULL
    );
    MPI_File_read_all(
        file, rvector_fft[grid_index], chunk_size,
        MPI_DOUBLE, MPI_STATUS_IGNORE
    );
    MPI_File_close(&file);
    free_rvector_array_type(&array_type, chunk_size);

    // Fourier-transform the density.
    // The global input buffer for the forward FFT is `rvector_fft`, where we
    // have put the data. Results are found in the output buffer `cvector_fft`,
    // which we then copy over to `kdensity`.
    forward_transform(grid_index);
    write_from_cvector(grid_index, kdensity[grid_index]);

    return 0;
}


/**
 * Writes the initial mass overdensity to a raw binary file.
 *
 * We expect the Fourier components `kdensity` of the density as input and
 * reverse-transform it to the physical domain before writing it out to a
 * file named `file_name` (which may include a relative or absolute path). It
 * is a headerless binary file of 64-bit floating-point numbers representing
 * the `density` array as a function of `(x, y, z)` in row-major (C-like)
 * index order.
 *
 * The `grid_index` is typically 0, referrring to the first of (in principle)
 * multiple `grids`.
 *
 * Returns 0 on success, a non-zero error code otherwise.
 */
int write_initial_density(
    double**   const kdensity,
    char*      const file_name,
    int        const grid_index,
    grid_data* const grids)
{
    // Get MPI rank (task number) of this process and total number of tasks.
    int rank, tasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);

    // The root task (a.k.a. main task) is the one at MPI rank 0.
    int const root = 0;

    // Reverse Fourier transform of the density.
    // The global input buffer for the reverse FFT is `cvector_fft`, so we
    // copy into it, then find the result in the output buffer `rvector_fft`.
    write_in_cvector(grid_index, kdensity[grid_index]);
    reverse_transform(grid_index);

    // Write in parallel to output file.
    MPI_File file;
    MPI_Datatype array_type;
    int chunk_size;
    create_rvector_array_type(grids[grid_index], &array_type, &chunk_size);
    int error = MPI_File_open(
        MPI_COMM_WORLD, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL, &file
    );
    if (error != MPI_SUCCESS && rank == root) {
        printf("Error: Could not open file \"%s\".\n", file_name);
        return 1;
    }
    MPI_File_set_view(
        file, 0, MPI_DOUBLE, array_type, "native", MPI_INFO_NULL
    );
    MPI_File_write_all(
        file, rvector_fft[grid_index], chunk_size,
        MPI_DOUBLE, MPI_STATUS_IGNORE
    );
    MPI_File_close(&file);
    free_rvector_array_type(&array_type, chunk_size);

    return 0;
}
