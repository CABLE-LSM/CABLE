MODULE cbl_hruff_mod

  IMPLICIT NONE

  PUBLIC HgtAboveSnow

CONTAINS

subroutine HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft, &
     SnowDepth, SnowDensity )

  !* Computes the height of the canopy above ground/snow surface when there
  !  is snow present
  !
  ! The height of canopy above snow level is simply the difference between
  ! the height of the canopy above the ground and the depth of snow.
  ! \[ HgtAboveSnow\_comp = HGT\_pft - SnowHeight \]
  ! The snow height is obtained from the snow amount and the snow density:
  ! \[ SnowHeight = 1.2 * \frac{SnowDepth}{SnowDensity} \]
  ! **the multipler 1.2 needs to be followed up**
  !
  ! Additionally, we ensure the snow density is larger than a minimum and the
  ! canopy height above snow is higher than 10 times the minimum surface
  ! roughness.

  implicit none

  !re-decl input args  
  integer  :: mp !! number of land points (-)
  !result to return
  real :: HeightAboveSnow(mp)  
    !! Effective height of canopy (output) (m) (also stored in rough%hruff)
  
  !re-decl input args  
  real :: z0surf_min !! minimum roughness length of surface (m)
  real :: HGT_pft(mp)  !! height of canopy with no snow (m)
  real :: SnowDepth(mp) !! amount of snow (\(mm/m^2\) liquid water)
  real :: SnowDensity(mp) !! density of snow (\(kg/m^3\))

  !local_vars: 
  real, parameter  :: fmin = 10.    ! scalar multiplier         
  real, parameter  :: SnowDensity_min = 100. ! min. snow density (kg/m3)
  
  real :: SnowDensity_eff(mp) 
  real :: HgtAboveSnow_min(mp)  
  real :: HgtAboveSnow_comp(mp)  

  !|### Structure
  ! * restricts the Effective snow density to be >= "a" set minimum
  SnowDensity_eff= max( SnowDensity_min, SnowDensity )
  
  !! * evaluates the mininimum allowed canopy height 
  !!  (fixed at 10 * min. surface roughness)
  HgtAboveSnow_min =  fmin * z0surf_min

  !! * calculates the height of the canopy above snow 
  HgtAboveSnow_comp =   HGT_pft - ( 1.2 * SnowDepth / SnowDensity_eff )

  !* * returns the effective canopy height above snow after comparison to
  ! the minimum allowed canopy height.
  HeightAboveSnow = max( HgtAboveSnow_min,  HgtAboveSnow_comp )

  return

END subroutine HgtAboveSnow

END MODULE cbl_hruff_mod
