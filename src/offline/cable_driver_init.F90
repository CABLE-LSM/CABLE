MODULE cable_driver_init_mod
  !! Module for CABLE offline driver initialisation.
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cable_driver_init

CONTAINS

  SUBROUTINE cable_driver_init()
    !! Model initialisation routine for the CABLE offline driver.

#ifdef __MPI__
    ! MPI specific initialisation
#else
    ! Serial specific initialisation
#endif

  END SUBROUTINE cable_driver_init

END MODULE cable_driver_init_mod
