MODULE cable_driver_init_mod
  !! Module for CABLE offline driver initialisation.
  USE cable_namelist_util, ONLY : get_namelist_file_name
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cable_driver_init

CONTAINS

  SUBROUTINE cable_driver_init()
    !! Model initialisation routine for the CABLE offline driver.

    !check to see if first argument passed to cable is
    !the name of the namelist file
    !if not use cable.nml
    CALL get_namelist_file_name()

#ifdef __MPI__
    ! MPI specific initialisation
#else
    ! Serial specific initialisation
#endif

  END SUBROUTINE cable_driver_init

END MODULE cable_driver_init_mod
