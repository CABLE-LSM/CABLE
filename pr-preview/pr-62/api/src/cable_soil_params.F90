MODULE cable_soil_params_mod

USE grid_constants_mod_cbl, ONLY: nsoil_max   ! # of soil types [9]
IMPLICIT NONE 

TYPE soilin_type

  CHARACTER(LEN=70) ::  desc(nsoil_max) ! decriptns of soil type 

   REAL, DIMENSION(nsoil_max) ::                                        &
      silt,    & !
      clay,    & !
      sand,    & !
      swilt,   & !
      sfc,     & !
      ssat,    & !
      bch,     & !
      hyds,    & !
      sucs,    & !
      rhosoil, & !
      css,     & !
      c3         !

END TYPE soilin_type

TYPE(soilin_type), SAVE  :: soilin

CHARACTER(LEN=70), DIMENSION(nsoil_max) ::  soil_desc    

CONTAINS

subroutine cable_soil_params()

! Gets parameter values for each vegetation type and soil type.
USE cable_def_types_mod, ONLY : mstype

implicit none

integer :: ERROR
integer, parameter :: namelist_unit=711178
integer :: j
CHARACTER(LEN=*), parameter :: iomessage='something wrong with your soil params file' 
CHARACTER(LEN=*), parameter :: nml_dir='./' 
CHARACTER(LEN=*), PARAMETER :: routinename='cable_soil_params'

NAMELIST / cable_soilparm / soilin

mstype = nsoil_max 
 
!SOIL parameters are assigned as TYPE soilin% but later mapped to soil%

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
write (6,*) "Reading CABLE_SOILPARM namelist..."

OPEN( namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'cable_soilparm.nml'),          &
      STATUS='old', POSITION='rewind', ACTION='read', IOSTAT  = ERROR )
IF ( ERROR /= 0 ) write (6,*) "Error opening  CABLE_SOILPARM namelist..."

READ(namelist_unit, NML = cable_soilparm, IOSTAT = ERROR )
IF ( ERROR /= 0 ) write (6,*) "Error reading  CABLE_SOILPARM namelist..."
                 
CLOSE(namelist_unit, IOSTAT = ERROR)

soil_desc = soilin%desc

End subroutine cable_soil_params

END MODULE cable_soil_params_mod

