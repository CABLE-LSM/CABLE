MODULE cable_soil_params_mod

USE grid_constants_mod_cbl, ONLY: nsoil_max   ! # of soil types [9]

IMPLICIT NONE 

CHARACTER(LEN=*), parameter :: nml_dir='./' 
CHARACTER(LEN=*), PARAMETER :: filename="soil.nml"

TYPE soilin_type

  CHARACTER(LEN=70) ::  desc(nsoil_max) ! decriptns of soil type 

  REAL, DIMENSION(nsoil_max) ::                                                &
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
    c3         ! IDK where this came from?

END TYPE soilin_type

TYPE(soilin_type) :: soilin

CHARACTER(LEN=70), DIMENSION(nsoil_max) ::  soil_desc    

CONTAINS

SUBROUTINE cable_read_soil_params()

IMPLICIT NONE

INTEGER :: error
INTEGER :: j

INTEGER, PARAMETER          :: namelist_unit=711178
CHARACTER(LEN=*), parameter :: iomessage='error with your soil params file'
CHARACTER(LEN=*), PARAMETER :: routinename='cable_soil_params'

NAMELIST / cable_soilparm / soilin

!SOIL parameters are assigned as TYPE soilin% but later mapped to soil%

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
WRITE(6,*) "Reading CABLE_SOILPARM namelist ", TRIM(nml_dir) // '/' // filename

OPEN( namelist_unit, FILE=(TRIM(nml_dir) // '/' // filename ),                 &
      STATUS='old', POSITION='rewind', ACTION='read', IOSTAT  = ERROR )

IF ( ERROR /= 0 ) THEN
  WRITE (6,*) "Error opening  CABLE_SOILPARM namelist..."
  STOP
ENDIF

READ(namelist_unit, NML = cable_soilparm, IOSTAT = ERROR )
                 
IF ( ERROR /= 0 ) THEN
  WRITE (6,*) "Error reading  CABLE_SOILPARM namelist..."
  STOP
ENDIF

CLOSE(namelist_unit, IOSTAT = ERROR)

soil_desc = soilin%desc

RETURN
END SUBROUTINE cable_read_soil_params

END MODULE cable_soil_params_mod

