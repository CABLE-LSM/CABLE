!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Error management for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: Developed by Gab Abramowitz and Harvey Davies
!
!
! ==============================================================================


MODULE cable_abort_module

  IMPLICIT NONE

CONTAINS

  !==============================================================================
  !
  ! Name: abort
  !
  ! Purpose: Prints an error message and stops the code
  !
  ! CALLed from: get_default_inits
  !              get_restart_data
  !              get_default_lai
  !              open_met_file
  !              get_met_data
  !              load_parameters
  !              open_output_file
  !              write_output
  !              read_gridinfo
  !              countpatch
  !              get_type_parameters
  !              readpar_i
  !              readpar_r
  !              readpar_rd
  !              readpar_r2
  !              readpar_r2d
  !              define_output_variable_r1
  !              define_output_variable_r2
  !              define_output_parameter_r1
  !              define_output_parameter_r2
  !              write_output_variable_r1
  !              write_output_variable_r2
  !              write_output_parameter_r1
  !              write_output_parameter_r1d
  !              write_output_parameter_r2
  !              write_output_parameter_r2d
  !
  !==============================================================================


  SUBROUTINE abort( message )

    ! Input arguments
    CHARACTER(LEN=*), INTENT(IN) :: message

    WRITE(*, *) message
    STOP 1

  END SUBROUTINE abort

  !==============================================================================
  !
  ! Name: nc_abort
  !
  ! Purpose: For NETCDF errors. Prints an error message then stops the code
  !
  ! CALLed from: get_restart_data
  !              extrarestart
  !              get_default_lai
  !              open_met_file
  !              get_met_data
  !              close_met_file
  !              load_parameters
  !              open_output_file
  !              write_output
  !              close_output_file
  !              create_restart
  !              read_gridinfo
  !              readpar_i
  !              readpar_r
  !              readpar_rd
  !              readpar_r2
  !              readpar_r2d
  !              define_output_variable_r1
  !              define_output_variable_r2
  !              define_output_parameter_r1
  !              define_output_parameter_r2
  !              write_output_variable_r1
  !              write_output_variable_r2
  !              write_output_parameter_r1
  !              write_output_parameter_r1d
  !              write_output_parameter_r2
  !              write_output_parameter_r2d
  !
  ! MODULEs used: netcdf
  !
  !==============================================================================


  SUBROUTINE nc_abort( ok, message )

    USE netcdf

    ! Input arguments
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER, INTENT(IN) :: ok

    WRITE(*,*) message ! error from subroutine
    WRITE(*,*) NF90_STRERROR(ok) ! netcdf error details

    STOP

  END SUBROUTINE nc_abort

  !==============================================================================
  !
  ! Name: range_abort
  !
  ! Purpose: Prints an error message and localisation information then stops the
  !          code
  !
  ! CALLed from: write_output_variable_r1
  !              write_output_variable_r2
  !
  ! MODULEs used: cable_def_types_mod
  !               cable_IO_vars_module
  !
  !==============================================================================



  SUBROUTINE range_abort(message,ktau,met,value,var_range,                       &
       i,xx,yy)

    USE cable_def_types_mod, ONLY: met_type
    USE cable_IO_vars_module, ONLY: latitude,longitude,landpt,lat_all,lon_all

    ! Input arguments
    CHARACTER(LEN=*), INTENT(IN) :: message

    INTEGER, INTENT(IN) ::                                                      &
         ktau, & ! time step
         i       ! landpt number of erroneous grid square

    INTEGER,INTENT(IN),OPTIONAL ::                                              &
         xx, & ! coordinates of erroneous grid square
         yy    ! coordinates of erroneous grid square


    TYPE(met_type),INTENT(IN) :: met  ! met data

    REAL,INTENT(IN) :: value ! value deemed to be out of range

    REAL,DIMENSION(2),INTENT(IN) :: var_range ! appropriate var range

    WRITE(*,*) "in SUBR range_abort"
    WRITE(*,*) message ! error from subroutine

    IF( PRESENT(yy) ) THEN ! i.e. using rectangular land/sea grid

       WRITE(*,*) 'Site lat, lon:',lat_all(xx,yy),lon_all(xx,yy)
       WRITE(*,*) 'Output timestep',ktau,                                       &
            ', or ', met%hod(landpt(i)%cstart),' hod, ',                  &
            INT( met%doy( landpt(i)%cstart) ),'doy, ',                   &
            INT( met%year(landpt(i)%cstart) )

    ELSE ! i.e. using compressed land only grid

       WRITE(*,*) 'Site lat, lon:', latitude(i), longitude(i)
       WRITE(*,*) 'Output timestep', ktau,                                       &
            ', or ', met%hod( landpt(i)%cstart ), ' hod, ',               &
            INT( met%doy( landpt(i)%cstart) ),'doy, ',                    &
            INT( met%year( landpt(i)%cstart) )

    END IF

    WRITE(*,*) 'Specified acceptable range (checks.f90):', var_range(1),        &
         'to',var_range(2)

    WRITE(*,*) 'Value:',value

    STOP

  END SUBROUTINE range_abort

  !==============================================================================
END MODULE cable_abort_module
