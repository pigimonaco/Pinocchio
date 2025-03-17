"""
Validation Script for PINOCCHIO Simulation
------------------------------------------
This script compiles and runs the PINOCCHIO simulation, then analyzes the output mass function (HMF). 

It performs the following steps:

- Compiles the PINOCCHIO code using the Makefile with selected compiler options.
  (See below and the PINOCCHIO documentation for available flags)
- Runs the simulation using MPI-OpenMP and an input parameter file.
- Compares the output mass function with a theoretical fit (or a reference HMF) and plots the results.
- Saves validation and simulation logs for further analysis.

### Usage
Run the script using:
```bash
python HMF_validation.py
```

Dependencies:
- Python: numpy, matplotlib
- System: MPI (e.g., OpenMPI, MPICH), a compatible compiler (e.g., nvc, gcc, mpicc)
- Ensure all dependencies are installed and aligned with `SYSTYPE` in the Makefile.

"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import shutil

# === Simulation Configuration ===
NUM_PROCESSES = 1  # Number of MPI processes
NUM_THREADS   = 1  # Number of OpenMP threads

# Paths
EXECUTABLE    = "pinocchio.x"
MAKEFILE_DIR  = "../src"  # Directory containing the Makefile
OUTPUT_DIR    = "../HMF_Validation"
HMF_FILE      = os.path.join(OUTPUT_DIR, "pinocchio.0.0000.test.mf.out")
LOG_FILE      = os.path.join(OUTPUT_DIR, "VALIDATION_log.txt")

# === Select Compiler ===
COMPILER_CC  = "nvc"   
COMPILER_CPP = "nvc++" 

# === Compilation Options (Set YES/NO) === #
COMPILATION_OPTIONS = {
    "DEBUG"        : "NO",
    "OMP"          : "YES",
    "GPU"          : "NO",
    "FULL_GPU"     : "YES",
    "SCOREP_NO"    : "NO",
    "NSIGHT"       : "NO",
    "ENERGY_CPU"   : "NO",
    "ENERGY_GPU"   : "NO"
}

# === Preprocessor Flags (Set YES/NO) === #
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

# Compile the code using Makefile
def compile_code():
    
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
        f"COMPILER_CC={COMPILER_CC}",
        f"COMPILER_CPP={COMPILER_CPP}",
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

# Run PINOCCHIO simulation with MPI support and parameter file
def run_simulation():
    
     # Ensure executable is in the output directory
    executable_source = os.path.join(MAKEFILE_DIR, EXECUTABLE)
    executable_dest   = os.path.join(OUTPUT_DIR, EXECUTABLE)
    shutil.copy(executable_source, executable_dest)

    if not os.path.isfile(executable_source):
        print(f"Error: Executable '{EXECUTABLE}' not found in {MAKEFILE_DIR}")
        return
    
    print(f"Executable copied to {OUTPUT_DIR}")
    
    # Define source and destination paths for parameter file and outputs file
    param_file_source = os.path.join("..", "example", "parameter_file")
    param_file_dest   = os.path.join(OUTPUT_DIR, "parameter_file")

    z_output_source = os.path.join("..", "example", "outputs")
    z_output_dest   = os.path.join(OUTPUT_DIR, "outputs")
    
    # Check and copy parameter file
    if not os.path.isfile(param_file_source):
        print(f"Error: Parameter file not found at {param_file_source}")
        return
    
    shutil.copy(param_file_source, param_file_dest)
    print(f"Parameter file copied to {OUTPUT_DIR}")

    # Check and copy outputs file
    if not os.path.isfile(z_output_source):
        print(f"Error: Outputs file not found at {z_output_source}")
        return
    
    shutil.copy(z_output_source, z_output_dest)
    print(f"Outputs file copied to {OUTPUT_DIR}")

    print("Running PINOCCHIO simulation...")
    
    # Set environment variable OMP_NUM_THREADS
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(NUM_THREADS)

    # Run the simulation
    run_log_file = os.path.join(OUTPUT_DIR, "log_RUN.txt")

    command = f"mpirun -np {NUM_PROCESSES} {executable_dest} {param_file_dest} | tee {run_log_file}"
    result = subprocess.run(command, cwd=OUTPUT_DIR, shell=True, capture_output=True, env=env)

    # Save execution logs
    with open(LOG_FILE, "a") as log:
        log.write("\n=== Simulation Log ===\n")
        log.write(f"Simulation log saved to: {run_log_file}\n")
    
    if result.returncode == 0:
        print("Simulation completed successfully.")
    else:
        print("Error in simulation:", result.stderr)
        exit(1)

# Plot HMF
def plot_hmf():
    hmf_file = HMF_FILE

    if not os.path.exists(hmf_file):
        print(f"Error: HMF output file {hmf_file} not found!")
        return None
    
    # Load data
    m, nm, fit = np.loadtxt(hmf_file, unpack=True, usecols=(0, 1, 5))

    # Compute residuals
    residuals    = (nm - fit) / fit
    avg_residual = np.mean(np.abs(residuals))

    # Create figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

    # First subplot: Mass function
    axs[0].plot(m, m * nm, label='PINOCCHIO', ls='-', lw=2, c='crimson')
    axs[0].plot(m, m * fit, label='Watson Fit', ls='--', lw=2, c='dodgerblue')
    axs[0].set_ylabel(r'M n(M) (Mpc$^{-3}$)', fontsize=14)
    axs[0].set_xlabel(r'M (M$_\odot$)', fontsize=14)
    axs[0].set_title('Mass function at z=0')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].legend(frameon=True)

    # Second subplot: Residuals
    axs[1].plot(m, residuals, linestyle='-', c='coral', label='Residuals')
    axs[1].axhline(0, color='black', linestyle='-', lw=1)  # Reference line at 0
    axs[1].set_ylabel(r'(PINOCCHIO - Watson)/Watson', fontsize=14)
    axs[1].set_xlabel(r'M (M$_\odot$)', fontsize=14)
    axs[1].set_title('Residuals at z=0')
    axs[1].set_xscale('log')

    # Ensure the output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save plot in the output directory
    plot_path = os.path.join(OUTPUT_DIR, 'HMF.png')
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
    avg_residual = plot_hmf()

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