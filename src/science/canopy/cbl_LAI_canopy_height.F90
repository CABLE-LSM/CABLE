module cbl_LAI_canopy_height_mod

contains

subroutine limit_HGT_LAI( HGT_pft_temp, LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, ntiles, &
                          tile_pts, tile_index, tile_frac,L_tile_pts, &
                          LAI_pft, HGT_pft, CLAI_thresh )

implicit none
real :: Clai_thresh                 !The minimum LAI below which a "cell" is considred NOT vegetated
integer :: land_pts, ntiles
integer :: tile_pts(ntiles) 
integer:: tile_index(land_pts,ntiles) 
real:: tile_frac(land_pts,ntiles) 
logical :: L_tile_pts(land_pts,ntiles) 
integer :: mp
real :: LAI_pft(land_pts, ntiles)
real :: HGT_pft(land_pts, ntiles)

!return vars
real :: HGT_pft_temp(land_pts,ntiles)
real :: LAI_pft_cbl(mp)
real :: HGT_pft_cbl(mp)

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
      
      LAI_pft_temp(i,N) = max(CLAI_thresh*.99,LAI_pft(i,N)) 
      if(N>13)  LAI_pft_temp(i,N) = 0.0 !to match offline Loobos
       ! hard-wired vegetation type numbers need to be removed
       IF(N < 5 ) THEN ! trees 
          HGT_pft_temp(i,N) = max(1.,HGT_pft(i,N)) 
       ELSE ! shrubs/grass
          HGT_pft_temp(i,N) = max(0.1, HGT_pft(i,N)) 
       ENDIF

    ENDIF

   ENDDO
ENDDO
  
  !surface_type = PACK(surface_type_temp, um1%L_TILE_PTS)
  LAI_pft_cbl  = PACK(LAI_pft_temp, L_TILE_PTS)
  HGT_pft_cbl  = PACK(HGT_pft_temp, L_TILE_PTS)

End subroutine limit_HGT_LAI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
End module cbl_LAI_canopy_height_mod

