MODULE cbl_hruff_mod

  IMPLICIT NONE

  PUBLIC HgtAboveSnow

CONTAINS
!Computes height above surface given that therre is snow present
!variable formerly known as rough%hruff
subroutine HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft, &
                         SnowDepth, SnowDensity )
  implicit none

  !re-decl input args  
  integer  :: mp
  !result to return
  real :: HeightAboveSnow(mp) 
  
  !re-decl input args  
  real :: z0surf_min
  real :: HGT_pft(mp) 
  real :: SnowDepth(mp) !snow depth (liquid water ) 
  real :: SnowDensity(mp) 

 
  !local_vars: 
  real, parameter  :: fmin = 10. ! [meters]?
  real, parameter  :: SnowDensity_min = 100. ! [meters]?
  
  real :: SnowDensity_eff(mp) 
  real :: HgtAboveSnow_min(mp)  
  real :: HgtAboveSnow_comp(mp)  

  ! restrict Effective snow density to be .GE. "a" minimum
  SnowDensity_eff= max( SnowDensity_min, SnowDensity ) 
  
  ! min. allowed canopy height (fixed @ 10* min. surface roughness)
  HgtAboveSnow_min =  fmin * z0surf_min 
  
  !Canopy Hgt above snow given computed snow depth & PFT height 
  HgtAboveSnow_comp =   HGT_pft - ( 1.2 * SnowDepth / SnowDensity_eff )

  ! Finally Set Effective canopy height above snow level 
  HeightAboveSnow = max( HgtAboveSnow_min,  HgtAboveSnow_comp )      

  return

END subroutine HgtAboveSnow

END MODULE cbl_hruff_mod
