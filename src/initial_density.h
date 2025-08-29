#ifndef INITIAL_DENSITY_H
#define INITIAL_DENSITY_H

int read_initial_density(
    double**         kdensity,
    char*      const file_name,
    int        const grid_index,
    grid_data* const grids);

int write_initial_density(
    double**   const kdensity,
    char*      const file_name,
    int        const grid_index,
    grid_data* const grid);

#endif
