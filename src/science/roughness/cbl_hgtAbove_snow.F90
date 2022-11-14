MODULE cbl_hruff_mod

  IMPLICIT NONE

  PUBLIC HgtAboveSnow

CONTAINS

subroutine HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft, &
     SnowDepth, SnowDensity )

  !* Subroutine computes height of canopy above ground/snow surface when there
  !  is snow present
  !
  !  inputs:
  !
  !  * mp - number of land points
  !  * z0surf_min - minimum roughness length of surface (m)
  !  * HGT_pft(mp) - height of canopy with no snow (m)
  !  * SnowDepth(mp) - amount of snow (mm m\(^{-2}\) liquid water)
  !  * SnowDensity(mp) - density of snow (in kg m\(^{-3}\))
  !
  !  outputs:
  !
  !  * HeightAboveSnow (mp) - effective height of canopy (m)
  !
  !  The output variable was formerly known as rough%hruff
  !

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
  real, parameter  :: fmin = 10.             
  real, parameter  :: SnowDensity_min = 100. 
  !* two internal parameters: fmin = 10 (scalar multiplier)
  !      and SnowDensity_min = 100.0 (in kg/m3)
  
  real :: SnowDensity_eff(mp) 
  real :: HgtAboveSnow_min(mp)  
  real :: HgtAboveSnow_comp(mp)  

  SnowDensity_eff= max( SnowDensity_min, SnowDensity )
  ! restricts the Effective snow density to be .GE. "a" minimum
  
  HgtAboveSnow_min =  fmin * z0surf_min
  ! evaluates min. allowed canopy height (fixed @ 10* min. surface roughness)
   
  HgtAboveSnow_comp =   HGT_pft - ( 1.2 * SnowDepth / SnowDensity_eff )
  !* The height of the canopy above snow level is simply the canopy height
  !  without snow minus the true depth of snow (in m).
  !  The depth of snow is evaluated from the amount of snow (in mm m\(^{-2}\
  !  of water) and the snow density.
  !
  ! \( h_{c,abovesnow} = \max[10 z_{0,min}, h_{c,nosnow} - 1.2 d_{snow}/\rho_{snow}] \)
  !
  !  A minimum value for the effective canopy height of
  !  '10* the minimum roughness length of the surface' is applied
  !  (a numerics requirement).
  !  Converting snow amount to snow depth requires the snow density;
  !  the snow density has a minimum value of 100 kg/m3 enforced.
  !
  !  **the multipler 1.2 needs to be followed up**

  HeightAboveSnow = max( HgtAboveSnow_min,  HgtAboveSnow_comp )
  ! Finally return the Effective canopy height above snow 

  return

END subroutine HgtAboveSnow

END MODULE cbl_hruff_mod
