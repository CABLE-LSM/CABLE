!******************************COPYRIGHT********************************************
! (c) CSIRO 2022.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT********************************************

MODULE init_active_tile_mask_mod

!------------------------------------------------------------------------------
! Description:
!   Initialises the JULES/CABLE grid array, which aligns JULES grid points
!   with CABLE land points
!
! This MODULE is USEd by:
!      cable_land_albedo_mod.F90
!
! This MODULE contains 1 public Subroutine:
!      init_active_tile_mask_cbl
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!------------------------------------------------------------------------------

IMPLICIT NONE
PUBLIC :: init_active_tile_mask_cbl
PRIVATE

CONTAINS

SUBROUTINE init_active_tile_mask_cbl(l_tile_pts, land_pts, nsurft, tile_frac )

! Description:
!   Nothing further to add to module description.

IMPLICIT NONE

LOGICAL, INTENT(OUT), ALLOCATABLE :: L_tile_pts(:,:)
INTEGER, INTENT(IN) :: land_pts, nsurft
REAL,    INTENT(IN) :: tile_frac(land_pts,nsurft)        !fraction of each surf type

!Local vars:
INTEGER :: i, j

! Determine active tiles map
IF ( .NOT. ALLOCATED(l_tile_pts)) ALLOCATE( l_tile_pts(land_pts, nsurft) )

l_tile_pts(:,:) = .FALSE.

DO j = 1, nsurft
  DO i = 1, land_pts
    IF ( tile_frac(i,j)  >   0.0 ) THEN
      l_tile_pts(i,j) = .TRUE.
    END IF
  END DO
END DO

RETURN

END SUBROUTINE init_active_tile_mask_cbl

END MODULE init_active_tile_mask_mod
