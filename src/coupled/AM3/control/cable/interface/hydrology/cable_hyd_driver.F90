MODULE cable_hyd_driv_mod
  
CONTAINS

SUBROUTINE cable_hyd_driver( land_pts, ntiles, tile_frac, lying_snow, SNOW_TILE, SURF_ROFF,&
                             SUB_SURF_ROFF, TOT_TFALL )

  implicit none

  !___ re-decl input args

  integer :: land_pts, ntiles
  

  REAL, INTENT(OUT), DIMENSION(LAND_PTS,NTILES) ::                    &
    SNOW_TILE   ! IN Lying snow on tiles (kg/m2)        

  REAL, INTENT(OUT), DIMENSION(LAND_PTS) ::                               &
    LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
    SUB_SURF_ROFF, & !
    SURF_ROFF,     & !
    TOT_TFALL        !

  !___ local vars


  REAL :: miss =0. 
 
  
return

END SUBROUTINE cable_hyd_driver


End module cable_hyd_driv_mod



