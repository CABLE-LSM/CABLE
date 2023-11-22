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
MODULE cbl_LAI_canopy_height_mod

!-----------------------------------------------------------------------------
! Description:
!   Restricts the range of canopy height and LAI inherited from JULES/UM
!   spatial maps
!
! This MODULE is USEd by:
!      cable_land_albedo_mod_cbl.F90
!
! This MODULE contains 1 public Subroutine:
!      limit_HGT_LAI
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE
PUBLIC :: limit_HGT_LAI
PRIVATE

CONTAINS

SUBROUTINE limit_HGT_LAI( LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, ntiles,      &
                          npft, tile_pts, tile_index, tile_frac,L_tile_pts,    &
                          LAI_pft, HGT_pft, CLAI_thresh )

! Description:
!   Nothing further to add to module description.

IMPLICIT NONE
INTEGER, INTENT(IN) :: mp
REAL, INTENT(OUT) :: LAI_pft_cbl(mp)
REAL, INTENT(OUT) :: HGT_pft_cbl(mp)
INTEGER, INTENT(IN) :: land_pts, ntiles, npft
REAL, INTENT(IN) :: LAI_pft(land_pts, ntiles)
REAL, INTENT(IN) :: HGT_pft(land_pts, ntiles)
REAL, INTENT(IN):: tile_frac(land_pts,ntiles)
REAL, INTENT(IN) :: Clai_thresh                 !The minimum LAI below which a "cell" is considred NOT vegetated
INTEGER, INTENT(IN) :: tile_pts(ntiles)
INTEGER, INTENT(IN):: tile_index(land_pts,ntiles)
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts,ntiles)

!local vars
REAL :: HGT_pft_temp(land_pts,ntiles) ! needed to filter spatail map
REAL :: LAI_pft_temp(land_pts,ntiles) ! needed to filter spatail map
INTEGER :: i,j, n

! init everywhere, even where tile_frac=0
LAI_pft_temp = 0.0
HGT_pft_temp = 0.0

DO n=1,ntiles
  DO j=1,tile_pts(n)

    i = tile_index(j,n)  ! It must be landpt index

    IF ( tile_frac(i,n)  >   0.0 ) THEN
      ! LAI set either just below threshold OR from INput field
      LAI_pft_temp(i,n) = MAX(CLAI_thresh*.99,LAI_pft(i,n))
      IF (n>npft)  LAI_pft_temp(i,n) = 0.0 !to match offline Loobos
       ! hard-wired vegetation type numbers need to be removed
      IF (n < 5 ) THEN ! set trees min. height
        HGT_pft_temp(i,n) = MAX(1.0,HGT_pft(i,n))
      ELSE  ! set shrubs/grass min. height
        HGT_pft_temp(i,n) = MAX(0.1, HGT_pft(i,n))
      END IF

    END IF

  END DO
END DO

! pack filtered JULE/UM maps to CABLE variables
LAI_pft_cbl  = PACK(LAI_pft_temp, l_tile_pts)
HGT_pft_cbl  = PACK(HGT_pft_temp, l_tile_pts)

END SUBROUTINE limit_HGT_LAI

END MODULE cbl_LAI_canopy_height_mod

