MODULE temp_module

LOGICAL, SAVE, PUBLIC :: vegparmnew    = .FALSE. ! using new format input file (BP dec 2007)
LOGICAL, SAVE, PUBLIC :: spinup        = .FALSE. ! model spinup to soil state equilibrium?
LOGICAL, SAVE, PUBLIC :: spincasa      = .FALSE. ! TRUE: CASA-CNP Will spin mloop times, FALSE: no spin up
LOGICAL, SAVE, PUBLIC :: CASAONLY      = .FALSE. ! ONLY Run CASA-CNP
LOGICAL, SAVE, PUBLIC :: l_casacnp     = .FALSE. ! using CASA-CNP with CABLE
LOGICAL, SAVE, PUBLIC :: l_landuse     = .FALSE. ! using CASA-CNP with CABLE
LOGICAL, SAVE, PUBLIC :: l_laiFeedbk   = .FALSE. ! using prognostic LAI
LOGICAL, SAVE, PUBLIC :: l_vcmaxFeedbk = .FALSE. ! using prognostic Vcmax

REAL, SAVE, PUBLIC :: delsoilM ! allowed variation in soil moisture for spin up
REAL, SAVE, PUBLIC :: delsoilT ! allowed variation in soil temperature for spin up
REAL, SAVE, PUBLIC :: delgwM = 1e-4

INTEGER, SAVE, PUBLIC :: LALLOC = 0            ! alloc coeff passed to spincasa

END MODULE temp_module
