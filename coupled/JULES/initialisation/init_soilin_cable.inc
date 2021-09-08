#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_soilin_cable(nml_dir)

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     rmdi

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE cable_other_constants_mod, ONLY: n_soiltypes

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

USE cable_params_mod, ONLY: soilin

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the non-vegetation parameters
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists
! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

CHARACTER(LEN=*), PARAMETER :: routinename='INIT_SOILPARM_CABLE'

!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm) ::                                                     &
      silt_io(n_soiltypes),                                                   &
      clay_io(n_soiltypes),                                                   &
      sand_io(n_soiltypes),                                                   &
      swilt_io(n_soiltypes),                                                  &
      sfc_io(n_soiltypes),                                                    &
      ssat_io(n_soiltypes),                                                   &
      bch_io(n_soiltypes),                                                    &
      hyds_io(n_soiltypes),                                                   &
      sucs_io(n_soiltypes),                                                   &
      rhosoil_io(n_soiltypes),                                                &
      css_io(n_soiltypes)

!-----------------------------------------------------------------------------
! Namelist definition
!-----------------------------------------------------------------------------
NAMELIST / cable_soilparm/ silt_io, clay_io, sand_io, swilt_io,               &
          sfc_io, ssat_io, bch_io, hyds_io, sucs_io, rhosoil_io, css_io

!-----------------------------------------------------------------------------
! Initialise namelist values before reading them
!-----------------------------------------------------------------------------
silt_io(:n_soiltypes)       = rmdi
clay_io(:n_soiltypes)    = rmdi
sand_io(:n_soiltypes)    = rmdi
swilt_io(:n_soiltypes)   = rmdi
sfc_io(:n_soiltypes)     = rmdi
ssat_io(:n_soiltypes)    = rmdi
bch_io(:n_soiltypes)     = rmdi
hyds_io(:n_soiltypes)    = rmdi
sucs_io(:n_soiltypes)    = rmdi
rhosoil_io(:n_soiltypes) = rmdi
css_io(:n_soiltypes)     = rmdi

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info(routinename, "Reading CABLE_SOILPARM namelist...")

! Open the CABLE soil parameters namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' //                             &
               'cable_soilparm.nml'),                                         &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(routinename,                                                 &
                 "Error opening namelist file cable_soilparm.nml " //         &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

READ(namelist_unit, NML = cable_soilparm, IOSTAT = error,                     &
                 IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(routinename,                                                 &
                 "Error reading namelist CABLE_SOILPARM " //                  &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(routinename,                                                 &
                 "Error closing namelist file cable_soilparm.nml " //         &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")


!-----------------------------------------------------------------------------
! Process the namelist values
!-----------------------------------------------------------------------------
! Copy values from dedicated I/O arrays into the soil parameter data type
soilin%silt(1:n_soiltypes)    = silt_io(1:n_soiltypes)
soilin%clay(1:n_soiltypes)    = clay_io(1:n_soiltypes)
soilin%sand(1:n_soiltypes)    = sand_io(1:n_soiltypes)
soilin%swilt(1:n_soiltypes)   = swilt_io(1:n_soiltypes)
soilin%sfc(1:n_soiltypes)     = sfc_io(1:n_soiltypes)
soilin%ssat(1:n_soiltypes)    = ssat_io(1:n_soiltypes)
soilin%bch(1:n_soiltypes)     = bch_io(1:n_soiltypes)
soilin%hyds(1:n_soiltypes)    = hyds_io(1:n_soiltypes)
soilin%sucs(1:n_soiltypes)    = sucs_io(1:n_soiltypes)
soilin%rhosoil(1:n_soiltypes) = rhosoil_io(1:n_soiltypes)
soilin%css(1:n_soiltypes)     = css_io(1:n_soiltypes)

!-----------------------------------------------------------------------------
! Check that all variables were present in the namelist.
! The namelist variables were initialised to rmdi.
!-----------------------------------------------------------------------------
error = 0
IF ( ANY( ABS( soilin%silt(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for silt")
END IF
IF ( ANY( ABS( soilin%clay(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for clay")
END IF
IF ( ANY( ABS( soilin%sand(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for sand")
END IF
IF ( ANY( ABS( soilin%swilt(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for swilt")
END IF
IF ( ANY( ABS( soilin%sfc(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for sfc")
END IF
IF ( ANY( ABS( soilin%ssat(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for ssat")
END IF
IF ( ANY( ABS( soilin%bch(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for bch")
END IF
IF ( ANY( ABS( soilin%hyds(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for hyds")
END IF
IF ( ANY( ABS( soilin%sucs(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for sucs")
END IF
IF ( ANY( ABS( soilin%rhosoil(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for rhosoil")
END IF
IF ( ANY( ABS( soilin%css(1:n_soiltypes) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for css")
END IF

IF ( error /= 0 )                                                             &
  CALL log_fatal(routinename,                                                 &
                 "Variable(s) missing from namelist - see earlier " //        &
                 "error message(s)")

RETURN

END SUBROUTINE init_soilin_cable

#endif
