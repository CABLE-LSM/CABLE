MODULE def_cable_grid_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: def_cable_grid

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='def_cable_grid_mod'

CONTAINS

SUBROUTINE def_cable_grid(mp)

USE init_active_tile_mask_mod,        ONLY: init_active_tile_mask_cbl
USE ancil_info,           ONLY: surft_pts

IMPLICIT NONE

INTEGER :: mp
!-----------------------------------------------------------------------------
! Description:
!   Allocates the CABLE model arrays using sizes determined during
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Determine the number of active tiles
mp = SUM(surft_pts)

CALL init_active_tile_mask_cbl()

END SUBROUTINE def_cable_grid

END MODULE def_cable_grid_mod
