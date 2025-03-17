"""
Validation Script for PINOCCHIO Simulation
------------------------------------------
This script compiles and runs a test PINOCCHIO simulation, then analyzes the Halo Mass Function (HMF). 

It performs the following steps:

- Compiles the PINOCCHIO code using the Makefile with selected compiler options.
  (See below and the PINOCCHIO documentation for available flags)
- Runs the simulation using MPI-OpenMP and an input parameter file.
- Compares the output HMF with either:
    1) Theoretical fit (Watson Fit) 
    2) Another mass function from a previous run.
- Saves validation and simulation logs for further analysis.


OUTPUT_DIR Structure:
The script will automatically create an output directory 'OUTPUT_DIR' as specified by the user. Inside this directory, it will:

- Copy the executable from 'MAKEFILE_DIR' to 'OUTPUT_DIR'
- Copy the parameter file from 'PARAM_SOURCE' to 'OUTPUT_DIR'
- Copy the outputs files from 'OUTPUTS_SOURCE' to 'OUTPUT_DIR'
- Modify the parameter file with the necessary settings ('RunFlag', 'BoxSize', 'GridSize').
- Perform HMF validation only at z = 0.

### Usage
Run the script using:
```bash
python HMF_validation.py
```

Dependencies:
- Python: numpy, matplotlib
- Ensure all dependencies to correctly run PINOCCHIO are installed and aligned with `SYSTYPE` in the Makefile.
"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import shutil


# =========================== #
#   COMPILER and EXEC name    #
# =========================== #
EXEC_NAME    = "pinocchio.x"
COMPILER_CC  = "mpicc"   
COMPILER_CPP = "mpicc++" 

# =========================== #
#          USER SETTINGS      #
# =========================== #
RUN_FLAG   = "test"      # Simulation tag
BOX_SIZE   = "128"       # Box size in Mpc/h
GRID_SIZE  = "128"       # Grid size (resolution)

NUM_PROCESSES = 1  # Number of MPI processes
NUM_THREADS   = 1  # Number of OpenMP threads

# =========================== #
#        FILE LOCATIONS       #
# =========================== #
MAKEFILE_DIR   = "../src"                                        # Directory containing the Makefile
OUTPUT_DIR     = "../HMF_Validation"                             # Directory containing the final validation output
LOG_FILE       = os.path.join(OUTPUT_DIR, "VALIDATION_log.txt")  # Name of the validation log

EXEC_SOURCE    = os.path.join(MAKEFILE_DIR, EXEC_NAME)           # Executable file source
EXEC_DEST      = os.path.join(OUTPUT_DIR, EXEC_NAME)             # Executbale file destionation
PARAM_SOURCE   = "../example/parameter_file"                     # Parameter file source
PARAM_DEST     = os.path.join(OUTPUT_DIR, "parameter_file")      # Parameter file destination
OUTPUTS_SOURCE = "../example/outputs"                            # Outputs file source
OUTPUTS_DEST   = os.path.join(OUTPUT_DIR, "outputs")             # Outputs file destination


# =========================== #
#    COMPILATION OPTIONS      #
# =========================== #
COMPILATION_OPTIONS = {
    "DEBUG"        : "NO",
    "OMP"          : "YES"
}

# =========================== #
#     PREPOCESSOR FLAGS       #
# =========================== #
PREPROCESSOR_FLAGS = {

    # Displacement LPT order # 
    "-DTWO_LPT"   : "YES",
    "-DTHREE_LPT" : "YES",

    # PLC reconstruction
    "-DPLC": "NO",

    # Dynamics of triaxial collapse #
    "-DELL_CLASSIC"  : "YES",
    "-DELL_SNG"      : "NO",
    "-DTABULATED_CT" : "NO",

    # Building groups and fragmentation #
    "-DCLASSIC_FRAGMENTATION": "NO",

    # Output #
    "-DSNAPSHOT"                  : "NO",
    "-DLIGHT_OUTPUT"              : "NO",
    "-DLONGIDS"                   : "NO",
    "-DDOUBLE_PRECISION_PRODUCTS" : "NO",

    # Beyond LambdaCDM models #

    # These are for neutrino cosmology
    "-DSCALE_DEPENDENT"         : "NO",
    "-DREAD_PK_TABLE"           : "NO",
    "-DONLY_MATTER_POWER"       : "NO",
    "-DRECOMPUTE_DISPLACEMENTS" : "NO",

    # Add also these for f(R) gravity #
    "-DMOD_GRAV_FR": "NO",
    "-DFR0=1.e-8"  : "NO",
    
    # Other options # 
    "-DWHITENOISE" : "NO",
    "-DNORADIATION": "YES",

    # - This option impacts CPU code only
    # - It is used by default using the GPU with OMP
    # "-DCUSTOM_INTERPOLATION": "NO"
}

# =========================== #
#       HELPER FUNCTIONS      #
# =========================== #

def copy_files(files_to_copy):
    """
    Copies multiple files from source to destination.
    
    Args:
        files_to_copy (dict): A dictionary where keys are source file paths and values are destination paths.
    """
    for src, dest in files_to_copy.items():
        if not os.path.exists(src):
            print(f"Error: File not found at {src}")
            exit(1)

        shutil.copy(src, dest)
        print(f"Copied {src} -> {dest}")

def update_param_file():
    """ Updates RunFlag, BoxSize, and GridSize in the parameter file. """
    
    if not os.path.isfile(PARAM_DEST):
        print(f"Error: Parameter file not found at {PARAM_DEST}")
        exit(1)

    updated_lines = []
    with open(PARAM_DEST, "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) > 1:
                key, value = parts[0], parts[1]
                if key == "RunFlag":
                    updated_lines.append(f"RunFlag                {RUN_FLAG}\n")
                elif key == "BoxSize":
                    updated_lines.append(f"BoxSize                {BOX_SIZE}\n")
                elif key == "GridSize":
                    updated_lines.append(f"GridSize               {GRID_SIZE}\n")
                else:
                    updated_lines.append(line)
            else:
                updated_lines.append(line)  # Keep empty lines unchanged

    with open(PARAM_DEST, "w") as file:
        file.writelines(updated_lines)

    print(f"Updated parameter file at {PARAM_DEST}")

def parse_parameters(param_file):
    """ Reads parameter file to extract RunFlag, BoxSize, and GridSize. """
    run_flag, box_size, grid_size = "example", "128", "128"
    
    with open(param_file, "r") as file:
        for line in file:
            parts = line.split()
            if len(parts) > 1:
                if "RunFlag" in parts[0]:
                    run_flag = parts[1]
                elif "BoxSize" in parts[0]:
                    box_size = parts[1]
                elif "GridSize" in parts[0]:
                    grid_size = parts[1]
    
    return run_flag, box_size, grid_size

def compile_code():
    """ Compiles the PINOCCHIO code. """
    os.makedirs(OUTPUT_DIR, exist_ok=True)  # Ensure output dir exists
    print("Compiling the code...")

    # Extract active compilation options
    compilation_str = " ".join(f"{key}={value}" for key, value in COMPILATION_OPTIONS.items())
    
    # Extract active preprocessor flags
    preprocessor_str = " ".join(flag for flag, enabled in PREPROCESSOR_FLAGS.items() if enabled == "YES")

    # Run `make clean`
    clean_result = subprocess.run(["make", "clean"], cwd=MAKEFILE_DIR, capture_output=True, text=True)
    if clean_result.returncode != 0:
        print("Error during 'make clean':", clean_result.stderr)
        exit(1)

    # Run `make` with dynamically selected options
    build_command = [
        "make",
        f"EXEC={EXEC_NAME}",
        f"CC={COMPILER_CC}",
        compilation_str,
        f"OPTIONS={preprocessor_str}",
        "-j16"  # Parallel compilation
    ]
    result = subprocess.run(build_command, cwd=MAKEFILE_DIR, capture_output=True, text=True)

    # Save compilation log
    with open(LOG_FILE, "w") as log:
        log.write("=== Compilation Log ===\n")
        log.write(" ".join(build_command) + "\n")
        log.write(result.stdout)
        log.write("\n")

    if result.returncode == 0:
        print("Compilation successful.")
    else:
        print("Error during compilation:", result.stderr)
        exit(1)

def run_simulation():
    """ Runs the PINOCCHIO simulation with MPI support and parameter file. """

    # Creating source-destination dictionary
    files_to_copy = {
        EXEC_SOURCE:    EXEC_DEST,
        PARAM_SOURCE:   PARAM_DEST,
        OUTPUTS_SOURCE: OUTPUTS_DEST
    }

    # Copy needed file in the OUTPUT_DIR (including the executable)
    copy_files(files_to_copy)
    
    # Update parameter file
    update_param_file()

    # Extract run tags from parameter file
    run_flag, box_size, grid_size = parse_parameters(PARAM_DEST)

    print("Running PINOCCHIO simulation...")
    
    # Set environment variable OMP_NUM_THREADS
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(NUM_THREADS)

    # Run the simulation
    run_log_file = os.path.join(OUTPUT_DIR, "log_RUN.txt")

    command = f"mpirun -np {NUM_PROCESSES} {EXEC_DEST} {PARAM_DEST} | tee {run_log_file}"
    result = subprocess.run(command, cwd=OUTPUT_DIR, shell=True, capture_output=True, env=env)

    # Save execution logs
    with open(LOG_FILE, "a") as log:
        log.write("\n=== Simulation INFO ===\n")
        log.write(f"Using RunFlag: {run_flag}, BoxSize: {box_size}, GridSize: {grid_size}")
        log.write(f"\nSimulation log saved to: {run_log_file}\n")
    
    if result.returncode == 0:
        print("Simulation completed successfully.")
    else:
        print("Error in simulation:", result.stderr)
        exit(1)

def plot_hmf(compare_with_fit=True, previous_run_flag=None):
    """ 
    Plots the Halo Mass Function (HMF) and allows the user to compare it either with:
    - The Watson Fit (default)
    - Another mass function from a previous run
    
    Args:
        compare_with_fit (bool): If True, compares with Watson Fit; otherwise, compares with another MF.
        previous_run_flag (str, optional): RunFlag of a previous simulation to use for comparison.
    """
    run_flag, _, _ = parse_parameters(PARAM_DEST)
    hmf_file = os.path.join(OUTPUT_DIR, f"pinocchio.0.0000.{run_flag}.mf.out")

    if not os.path.exists(hmf_file):
        print(f"Error: HMF output file {hmf_file} not found!")
        return None
    
    # Load data for the current run
    m, nm, fit = np.loadtxt(hmf_file, unpack=True, usecols=(0, 1, 5))

    if compare_with_fit:
        comparison_label = "Watson Fit"
        comparison_data  = fit
        HMF_plot_name    = "MHF_Validation_with_Watson_fit.png"
    else:
        while True:
            prev_hmf_file = os.path.join(OUTPUT_DIR, f"pinocchio.0.0000.{previous_run_flag}.mf.out")

            if os.path.exists(prev_hmf_file):
                break  # File exists, proceed
            else:
                print(f"Error: The file 'pinocchio.0.0000.{previous_run_flag}.mf.out' is not present. Please check and enter a correct RunFlag.")
                previous_run_flag = input("Enter the previous RunFlag for comparison: ").strip()
        _ , prev_nm, _ = np.loadtxt(prev_hmf_file, unpack=True, usecols=(0, 1, 5))
        comparison_label = f"Reference Run ({previous_run_flag})"
        comparison_data = prev_nm
        HMF_plot_name    = "MHF_Validation_with_Reference_Run.png"
        
    # Compute residuals
    residuals = (nm - comparison_data) / comparison_data

    # Mask out residuals that are NaN or exactly -1
    valid_residuals = residuals[~np.isnan(residuals) & (residuals != -1)]

    # Compute the average of the valid residuals
    avg_residual = np.mean(np.abs(valid_residuals))

    # Create figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

    # First subplot: Mass function
    axs[0].plot(m, m * nm, label=f'PINOCCHIO ({run_flag})', ls='-', lw=2, c='crimson')
    axs[0].plot(m, m * comparison_data, label=comparison_label, ls='--', lw=2, c='dodgerblue')
    axs[0].set_ylabel(r'M n(M) (Mpc$^{-3}$)', fontsize=14)
    axs[0].set_xlabel(r'M (M$_\odot$)', fontsize=14)
    axs[0].set_title('Mass function at z=0')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].legend(frameon=True)

    # Second subplot: Residuals
    axs[1].plot(m, residuals, linestyle='-', c='coral', label='Residuals')
    axs[1].axhline(0, color='black', linestyle='-', lw=1)  # Reference line at 0
    axs[1].set_ylabel(r'(PINOCCHIO - Comparison)/Comparison', fontsize=14)
    axs[1].set_xlabel(r'M (M$_\odot$)', fontsize=14)
    axs[1].set_title('Residuals at z=0')
    axs[1].set_xscale('log')

    # Ensure the output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save plot in the output directory
    plot_path = os.path.join(OUTPUT_DIR, f"{HMF_plot_name}")
    plt.savefig(plot_path)
    print(f"Plot saved to {plot_path}")

    # Append results to log
    with open(LOG_FILE, "a") as log:
        log.write("\n=== Analysis Results ===\n")
        log.write(f"HMF Average Residual: {avg_residual:.5e}\n")
        log.write(f"Plot saved to: {plot_path}\n")
          
        log.write("\nCompilation Options:\n")
        for key, value in COMPILATION_OPTIONS.items():
            log.write(f"  {key}: {value}\n")

        log.write("\nPreprocessor Flags:\n")
        for flag, enabled in PREPROCESSOR_FLAGS.items():
            log.write(f"  {flag}: {enabled}\n")


    return avg_residual

# Main execution
def main():

    compile_code()
    run_simulation()

    # Ask the user whether to compare with the Watson Fit or a previous RunFlag
    compare_with_fit  = input("Compare the new HMF with Watson Fit? (yes/no): ").strip().lower() == 'yes'
    previous_run_flag = None

    if not compare_with_fit:
        while True:  # Keep asking until they enter something
            previous_run_flag = input("Enter the previous RunFlag for comparison: ").strip()
            if previous_run_flag:
                break
            print("Error: You must enter a valid RunFlag.")

    avg_residual = plot_hmf(compare_with_fit, previous_run_flag)

    # Print final summary
    print("\n=== Summary ===")
    print("Compilation Options:")
    for key, value in COMPILATION_OPTIONS.items():
        print(f"  {key}: {value}")

    print("\nPreprocessor Flags:")
    for flag, enabled in PREPROCESSOR_FLAGS.items():
        print(f"  {flag}: {enabled}")

    print(f"VALIDATION Log saved to: {LOG_FILE}")
    print(f"HMF Average Residual: {avg_residual:.5e}")
    print("Validation completed successfully.")


if __name__ == "__main__":
    main()