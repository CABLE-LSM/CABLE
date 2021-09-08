!jhan:Althoughthis isonly calling subrs and not using data - still needs revision to only call once and use L_tile_pts from here
module init_active_tile_mask_mod
  public L_tile_pts
  public init_active_tile_mask_cbl
!H! remove SAVE attr later 
  !mask TRUE where tile fraction is greater than zero
  LOGICAL, allocatable,SAVE :: L_tile_pts(:,:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_active_tile_mask_cbl()

USE cable_types_mod,          ONLY: l_tile_pts_types => l_tile_pts
USE ancil_info,               ONLY: frac_surft, land_pts
USE jules_surface_types_mod,  ONLY: ntype

IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!   Initialises the JULES/CABLE grid array, which aligns JULES grid points
!   with CABLE land points
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

INTEGER :: i, j

! Determine active tiles map
IF ( .NOT. ALLOCATED(l_tile_pts)) ALLOCATE( l_tile_pts(land_pts, ntype) )

l_tile_pts(:,:) = .FALSE.

DO j = 1, ntype
  DO i = 1, land_pts
    IF ( frac_surft(i,j)  >   0.0 ) THEN
      l_tile_pts(i,j) = .TRUE.
    END IF
  END DO
END DO 

l_tile_pts_types = l_tile_pts

RETURN

END SUBROUTINE init_active_tile_mask_cbl

End module init_active_tile_mask_mod
