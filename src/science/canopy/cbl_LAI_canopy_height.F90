module cbl_LAI_canopy_height_mod
  !* This MODULE contains the SUBROUTINE [[limit_HGT_LAI]].  This routine
  ! places limits on the values of the leaf area and canopy height to ensure
  ! realism.

contains

subroutine limit_HGT_LAI( HGT_pft_temp, LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, ntiles, &
                          tile_pts, tile_index, tile_frac,L_tile_pts, &
                          LAI_pft, HGT_pft, CLAI_thresh )
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

implicit none
real :: Clai_thresh                 !! The minimum LAI below which a tile is considered to be NOT vegetated (m\(^2\)m\(^{-2}\))
integer :: land_pts, ntiles         !! number of land points, max number of tiles that can be carried in any land point (-)
integer :: tile_pts(ntiles)         !! number of land_pts in the array that include an active tile of a particular tile type (-) 
integer:: tile_index(land_pts,ntiles) !! specifies the tile types within each land point (-)
real:: tile_frac(land_pts,ntiles)   !! fraction of cell assigned to each tile (-)
logical :: L_tile_pts(land_pts,ntiles) !! switch to denote whether a tile is active (-)
integer :: mp                       !! total number of tiles carried by land vector (-)
real :: LAI_pft(land_pts, ntiles)   !! IN leaf area (m\(^2\)m\(^{-2}\))
real :: HGT_pft(land_pts, ntiles)   !! IN canopy height (m)

!return vars
real :: HGT_pft_temp(land_pts,ntiles) !! OUT adjusted canopy height in (land_pts,ntile) form (m)
real :: LAI_pft_cbl(mp)             !! OUT adjusted leaf area in (mp) form (m\(^2\)m\(^{-2}\))
real :: HGT_pft_cbl(mp)             !! OUT adjusted canopy height in (mp) form (m)

!local vars
real :: LAI_pft_temp(land_pts,ntiles)
integer :: i,j, N 

!Retain init where tile_frac=0
LAI_pft_temp = 0. 
HGT_pft_temp = 0.

DO N=1,NTILES
  DO J=1,TILE_PTS(N)
      
    i = TILE_INDEX(j,N)  ! It must be landpt index

    IF( TILE_FRAC(i,N) .gt. 0.0 ) THEN
      
       ! hard-wired vegetation type numbers need to be removed
       IF(N < 5 ) THEN ! trees 
        LAI_pft_temp(i,N) = max(CLAI_thresh,LAI_pft(i,N)) 
          HGT_pft_temp(i,N) = max(1.,HGT_pft(i,N)) 
      ELSE IF(N > 4 .AND. N < 14 ) THEN  ! shrubs/grass
        LAI_pft_temp(i,N) = max(CLAI_thresh,LAI_pft(i,N)) 
          HGT_pft_temp(i,N) = max(0.1, HGT_pft(i,N)) 
      ELSE IF(N > 13 ) THEN ! non-vegetated
        LAI_pft_temp(i,N) = 0.0
        HGT_pft_temp(i,N) = 0.0
       ENDIF

    ENDIF

   ENDDO
ENDDO
  
  !surface_type = PACK(surface_type_temp, um1%L_TILE_PTS)
  LAI_pft_cbl  = PACK(LAI_pft_temp, L_TILE_PTS)
  HGT_pft_cbl  = PACK(HGT_pft_temp, L_TILE_PTS)

End subroutine limit_HGT_LAI

!==============================================================================
 
End module cbl_LAI_canopy_height_mod

