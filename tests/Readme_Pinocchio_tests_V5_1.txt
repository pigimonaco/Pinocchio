** PINOCCHIO usage for ICs generation **

1) MakeFile options ON:

* OPTIONS += -DTWO_LPT     # For 2LPT ICs
* OPTIONS += -DTHREE_LPT   # For 3LPT ICs 

NOTE: If both are compiled, the code will use 3LPT. To use 2LPT initial conditions, comment out OPTIONS += -DTHREE_LPT

* OPTIONS += -DSNAPSHOT    

2) Runtime example:

mpirun -np #task executable_name paramenter_file 3   # When run with the argument 3 after the parameter file, Pinocchio will operate in a special mode designed exclusively for generating the LPT initial conditions for an N-body


** READ_PK_TABLE e SCALE DEPENDENT GROWTH RATE **

The behavior in this case is well described in the DOCUMENTATION (Section 8) on Git (branch fivedotone). There are minor details to be added to correctly run with the PKs table and the transfer function:

1) MakeFile options to enable reading the PKs Table:

* OPTIONS += -DREAD_PK_TABLE 
* OPTIONS += -DSCALE_DEPENDENT

2) ParameterFile setup:

* FileWithInputSpectrum CAMBTable

3) Makefile options to enable reading the TransferFunction:

* OPTIONS += -DONLY_MATTER_POWER       

NOTE_1 for TransferFunction usage: 
The input transfer function file must contain 7 columns in the following order:

|k/h| |delta_cdm| |delta_phot| |delta_nu_1 (massless neutrinos)| |delta_nu_2 (massive neutrinos)| |delta_baryon| |delta_tot|

NOTE_2 for the CAMBRedshiftsFile:
The input redshift file for which the PKs are available must contain two columns in the following order:

|000| |highest_z|
|001| |second highest_z| 
  .      .
  .      .
  .      .
|***| |lowest_z|


** MODIFIED GRAVITY RUN ** 

1) MakeFile options to enable MODIFIED GRAVITY RUN:    # NOTE_1: This part should be used with extreme caution as it is still under development	

* OPTIONS += -DSCALE_DEPENDENT
* OPTIONS += -DMOD_GRAV_FR
* OPTIONS += -DFR0=1.e-8
* OPTIONS += -DTABULATED_CT  # The modified gravity part only work with collapse time calculation from a table (interpolation)

We offer three different types of interpolation options:

* OPTIONS += -DTRILINEAR:       Slowest but most precise.
* OPTIONS += -DBILINEAR_SPLINE: A good compromise between speed and precision.
* OPTIONS += -DALL_SPLINE:      Fastest but least precise

2) ParameterFile setup:

CTtablefile none (if not provided)  # The name of the precalculated table of collapse times must be specified if available


** GENERAL NOTE for running with one external PK file ** 
When using an external PK file, set FileWithInputSpectrum != no in the ParameterFile and specify the file name (or its path). 
For this type of run, there is no need to compile with the READ_PK_TABLE option in the Makefile.
 



