MODULE cbl_hruff_mod

  !* The procedure in this module evaluates the canopy height
  !  given the effect of any snow present.  

  IMPLICIT NONE

  PUBLIC HgtAboveSnow

CONTAINS


subroutine HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft, &
     SnowDepth, SnowDensity )

  !* This subroutine computes the height of the canopy above ground/snow
  !  surface when there is snow present. 
  
  implicit none

  !re-decl input args  
  integer  :: mp !! Number of land tiles (-)
  !result to return
  real :: HeightAboveSnow(mp) 
    !! Output. Effective height of canopy, known as `rough%hruff` elsewhere. (m)
  
  !re-decl input args  
  real :: z0surf_min !! minimum roughness length of the surface (m)
  real :: HGT_pft(mp) !! height of the canopy without snow (m)
  real :: SnowDepth(mp) !! amount of snow (mm m\(^{-2}\) liquid water) 
  real :: SnowDensity(mp) !! density of snow (kg m\(^{-3}\))

  !local_vars: 
  real, parameter  :: fmin = 10.  
    ! Multiplier scalar to fix the minimum allowed canopy height. (-)           
  real, parameter  :: SnowDensity_min = 100. 
    ! Minimum allowed snow density (kg m\(^{-3}\))
  
  real :: SnowDensity_eff(mp) 
  real :: HgtAboveSnow_min(mp)  
  real :: HgtAboveSnow_comp(mp)  

  ! restricts the Effective snow density to be >= a set minimum
  SnowDensity_eff= max( SnowDensity_min, SnowDensity )
  
  ! evaluates mininimum allowed canopy height for numerical stability 
  ! (fixed at 10 * mininimum surface roughness)
  HgtAboveSnow_min =  fmin * z0surf_min
   
  HgtAboveSnow_comp =   HGT_pft - ( 1.2 * SnowDepth / SnowDensity_eff )
  !* The height of the canopy above snow level is simply the canopy height
  ! without snow minus the true depth of snow (in m).
  ! The depth of snow is evaluated from the amount of snow (in mm m\(^{-2}\)
  ! of liquid water) and the snow density.
  !
  ! \[ h_{c,abovesnow} = \max[10*z_{0,min}, h_{c,nosnow} - 1.2*d_{snow}/\rho_{snow}] \]
  !
  !  A minimum value for the effective canopy height of
  !  '10 * the minimum roughness length of the surface' is applied
  !  (a numerics requirement).
  !
  !  The snow density has a minimum value of 100 kg m\(^{-3}\) enforced when
  !  used to convert snow amount to snow depth.
  !
  !  ** Warning: The multiplier 1.2 needs to be followed up**

  ! Finally return the Effective canopy height above snow 
  HeightAboveSnow = max( HgtAboveSnow_min,  HgtAboveSnow_comp )

  return

END subroutine HgtAboveSnow

END MODULE cbl_hruff_mod
