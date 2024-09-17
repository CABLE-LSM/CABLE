MODULE cable_model_env_mod

USE cable_model_env_opts_mod, ONLY: icycle

IMPLICIT NONE

NAMELIST / cable_model_environment /                                           &
  icycle 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CABLE_MODEL_ENV_MOD'

CONTAINS

SUBROUTINE set_derived_variables_cable_model_env

USE cable_model_env_opts_mod, ONLY: l_casacnp 

IMPLICIT NONE

IF( icycle == 1 ) THEN
  l_casacnp = .TRUE.
  WRITE(6,*) "CABLE-CASA-CNP configured for Carbon cycle"
ELSEIF( icycle == 2 ) THEN
  l_casacnp = .TRUE.
  WRITE(6,*) "CABLE-CASA-CNP configured for Carbon, Nitrogen cycle"
ELSEIF( icycle == 3 ) THEN
  l_casacnp = .TRUE.
  WRITE(6,*) "CABLE-CASA-CNP configured for Carbon, Nitrogen, Phosphorus cycle"
ELSE
  l_casacnp = .FALSE.
  WRITE(6,*) "CABLE-CASA-CNP configured for NO biogeochemical cycle"
ENDIF

RETURN

END SUBROUTINE set_derived_variables_cable_model_env

SUBROUTINE read_nml_cable_model_env(unitnumber)

! Description:
!  Read the cable_model_env namelist

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*),   PARAMETER :: RoutineName='READ_NML_CABLE_MODEL_ENV'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 1

TYPE my_namelist
  !!SEQUENCE
  INTEGER :: icycle
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = cable_model_environment, IOSTAT = errorstatus,    &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist cable_model_environment",iomessage)

  my_nml % icycle = icycle
END IF
  
CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  icycle = my_nml % icycle

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_cable_model_env

END MODULE cable_model_env_mod

