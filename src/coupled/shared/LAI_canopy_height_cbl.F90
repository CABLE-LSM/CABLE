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

SUBROUTINE Limit_HGT_LAI( HGT_pft_temp, LAI_pft_cbl, HGT_pft_cbl, mp, land_pts,&
                          ntiles, npft, tile_pts, tile_index, tile_frac,       &
                          L_tile_pts, LAI_pft, HGT_pft, CLAI_thresh )

!*## Purpose
!
! This SUBROUTINE checks that input values of leaf area and canopy height
! lie between set values (defined by PFT) and overwrites if necessary.
! These limits are needed to ensure that the canopy roughness properties of
! the land point are well behaved (see [[ruff_resist]]).
!
! The mapping from grid cells (land_pts,ntiles) to the (mp) structure
! is also undertaken.
!
! **This SUBROUTINE is active when CABLE is run within a coupled model
! (i.e. ACCESS)** 
!
!## Method
!
! The cell-level input (land_pts,ntile) variables for leaf area `LAI_pft` and
! canopy height `HGT_pft` are checked for
!
! - whether the land_pt has a non-zero fraction of that tile (if zero
!   fraction then LAI and canopy height are set to zero)
! - whether the LAI lies above a minimum value `CLAI_thresh`
! - whether the canopy height lies above a minimum value (which is PFT dependent)
! - whether non-vegetated tiles have a non-zero canopy height.
!
! Outputs are `LAI_pft_cbl` and `HGT_pft_cbl`
!
! **WARNINGS**
!
! - INTENT statements need to be added to the argument lists
! - hardwired indexing is used throughout.  This SUBROUTINE assumes that
!     1. non-vegetated tiles take indexes 14 and onwards
!     2. tall vegetation occupies tile indexes 1-4
!     3. other vegetation (shrubs/grasses/crops/wetlands) occupy tile indexes
!     5-13
! - the limits on canopy height are hardwired
!     1. canopy height for tall vegetation is limited to be greater than 1m
!     2. canopy height for other vegetation is limited to be greater than 0.1m

USE cable_surface_types_mod, ONLY: shrub_cable, aust_xeric, wetland
                                          

IMPLICIT NONE

INTEGER, INTENT(IN)  :: land_pts, ntiles, npft        !! # land points here 
                                                      !! max # of tiles, PFTs
!! CABLE 1-D vector of length "mp" active tiles
INTEGER, INTENT(IN)  :: mp                            !! total number of tiles 

!! range limited LAI/ canopy height 
REAL,    INTENT(OUT) :: HGT_pft_temp(land_pts,ntiles) !! UM dimensions (m)

!! PACKed to  CABLE 1-D vector of length mp
REAL,    INTENT(OUT) :: LAI_pft_cbl(mp)               !! veg%vlai 
REAL,    INTENT(OUT) :: HGT_pft_cbl(mp)               !! veg%hc

!! From UM spatial maps via restart or ancillary. 
REAL,    INTENT(IN)  :: tile_frac(land_pts,ntiles)    !! frac. veg/non-veg type
INTEGER, INTENT(IN)  :: tile_pts(ntiles)              !! # tiles per PFT type
INTEGER, INTENT(IN)  :: tile_index(land_pts,ntiles)   !! land_pt index per tile
!! updates monthly
REAL,    INTENT(IN)  :: LAI_pft(land_pts, npft)       !! LAI (m\(^2\)m\(^{-2}\))
REAL,    INTENT(IN)  :: HGT_pft(land_pts, npft)       !! canopy height (m)         
!! logical mask. TRUE if tilefrac > 0 0 (OR some threshold)
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts,ntiles)

!! scalar constant threshold for "cell" to be considred vegetated
REAL,    INTENT(IN)  :: Clai_thresh                   ! minimum LAI threshold 

!local vars
INTEGER :: i,j, n
REAL :: LAI_pft_temp(land_pts,ntiles) ! needed to filter spatail map
REAL :: hgt_max

!Retain init where tile_frac=0
LAI_pft_temp(:,:) = 0.0
HGT_pft_temp(:,:) = 0.0

! surface types ordered as vegetated = 1:13. non vegetated = 14:17 
DO n=1,ntiles
  DO j=1,tile_pts(n)

    i = tile_index(j,n)  ! landpt index

    IF( tile_frac(i,n) .GT. 0.0 ) THEN

      ! Check canopy height & LAI > lower limit for vegetated PFTs 
      IF ( n <= npft ) THEN

        ! Max height for all trees
        hgt_max = 1.0 
        ! Max height for shrubs, grasses, crops, tundra and wetlands
        IF ( n >= shrub_cable .AND. n <=  wetland ) THEN 
           hgt_max = 0.1
        END IF
        
        ! 99% factor allows LAI_pft_temp to drop below CLAI_thresh, and 
        ! subsequently trigger appropriate science in dryLeaf, photosynthesis
        ! Coherent treatment of "vegetated" status/science would be preferable 
        LAI_pft_temp(i,n) = MAX( 0.99*CLAI_thresh, LAI_pft(i,n) )
        HGT_pft_temp(i,n) = MAX( hgt_max, HGT_pft(i,n) )

      ENDIF

    ENDIF

  ENDDO
ENDDO

!surface_type = PACK(surface_type_temp, um1%L_TILE_PTS)
LAI_pft_cbl  = PACK(LAI_pft_temp, l_tile_pts)
HGT_pft_cbl  = PACK(HGT_pft_temp, l_tile_pts)

END SUBROUTINE limit_HGT_LAI

END MODULE cbl_LAI_canopy_height_mod

