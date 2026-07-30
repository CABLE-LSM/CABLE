MODULE casa_files_type_mod

IMPLICIT NONE

TYPE casa_files_type

  CHARACTER(LEN=99) :: cnpbiome     = '' ! biome-specific BGC parameters
  CHARACTER(LEN=99) :: cnppoint     = '' ! file for point-specific BGC inputs
  CHARACTER(LEN=99) :: cnpepool     = '' ! file for end-of-run pool sizes
  CHARACTER(LEN=99) :: cnpipool     = '' ! file for inital pool sizes
  CHARACTER(LEN=99) :: cnpmetin     = '' ! met file for spin up
  CHARACTER(LEN=99) :: cnpmetout    = '' ! met file for spin up
  CHARACTER(LEN=99) :: cnpspin      = '' ! input file for spin up
  CHARACTER(LEN=99) :: dump_cnpspin = '' ! dump file for spinning casa-cnp
  CHARACTER(LEN=99) :: phen         = '' ! leaf phenology datafile
  CHARACTER(LEN=99) :: cnpflux      = '' ! modelled mean yearly CNP fluxes
  CHARACTER(LEN=99) :: c2cdumppath  = '' ! cable2casa dump for casa spinup
  LOGICAL           :: l_ndep       = .FALSE.
  CHARACTER(LEN=99) :: ndep         = '' ! N deposition input file

END TYPE casa_files_type

TYPE(casa_files_type) :: casafile

END MODULE casa_files_type_mod

