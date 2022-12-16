!******************************************************************************
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) can be found
! at https://github.com/CABLE-LSM/CABLE/blob/main/
!
!******************************************************************************

MODULE hruff_eff_LAI_mod_cbl

!-----------------------------------------------------------------------------
! Description:
!   Computes the height above the surface given that there is snow present
!   and the effective LAI of the canopy given the effect of snow
!
! This MODULE is USEd in:
!     cable_land_albedo_mod_cbl.F90 (JULES)
!
! This MODULE contains 2 public Subroutines:
!     HgtAboveSnow,
!     LAI_eff
!
! Module specific documentation: https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow: https://trac.nci.org.au/trac/cable/wiki/TBC
!-----------------------------------------------------------------------------

!* The first procedure in this module evaluates the canopy height
!  given the effect of any snow present.  

!* The secone procedure in this module computes the effective LAI of a canopy
!  given the effect of any snow present

IMPLICIT NONE

PUBLIC :: HgtAboveSnow
PUBLIC :: LAI_eff

CONTAINS

SUBROUTINE HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft,             &
                         SnowDepth, SnowDensity )

!* This subroutine computes the height of the canopy above ground/snow
!  surface when there is snow present. 

IMPLICIT NONE

INTEGER, INTENT(IN)   :: mp    !! Number of tiles (-)     

REAL, INTENT(OUT) :: HeightAboveSnow(mp) 
    !! Output. Effective height of canopy, known as `rough%hruff` elsewhere. (m)

REAL, INTENT(IN) :: z0surf_min      !! minimum roughness length of surface (m)
REAL, INTENT(IN) :: HGT_pft(mp)     !! height of the canopy without snow (m)      
REAL, INTENT(IN) :: SnowDepth(mp)   !! snow amount (mm m\(^{-2}\) liquid water)
REAL, INTENT(IN) :: SnowDensity(mp) !! density of snow (kg m\(^{-3}\))           

!local_vars:
REAL, PARAMETER  :: fmin = 10.0 ! [meters]?
    ! Multiplier scalar to fix the minimum allowed canopy height. (-)           
REAL, PARAMETER  :: SnowDensity_min = 100.0 ! min. snow density
    ! Minimum allowed snow density (kg m\(^{-3}\))

REAL :: SnowDensity_eff(mp)         ! Effective snow density range restricted
REAL :: HgtAboveSnow_min(mp)        ! min. canopy height
REAL :: HgtAboveSnow_comp(mp)       ! computed canopy height above snow

! restricts the Effective snow density to be >= a set minimum
SnowDensity_eff= MAX( SnowDensity_min, SnowDensity )

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

! Finally Set Effective canopy height above snow level
HeightAboveSnow = MAX( HgtAboveSnow_min,  HgtAboveSnow_comp )

RETURN

END SUBROUTINE HgtAboveSnow

SUBROUTINE LAI_eff( mp, lai_pft, Hgt_PFT, HgtAboveSnow,                        &
                    reducedLAIdue2snow )

!* This subroutine computes the leaf-area index of a canopy when there is snow
! present. The formulae assume that the leaf area is distributed uniformly
! in the vertical within the canopy.

IMPLICIT NONE

INTEGER, INTENT(IN)   :: mp           ! CABLE VECTOR LENGTH
! return result - considered LAI seen given snow coverage
REAL, INTENT(OUT) :: reducedLAIdue2snow(mp)
REAL, INTENT(IN) :: lai_pft(mp)       ! LAI
REAL, INTENT(IN) :: HGT_pft(mp)       ! canopy height
REAL, INTENT(IN) :: HgtAboveSnow(mp)  ! computed canopy height above snow

!local_vars:
REAL :: Hgt_eff(mp)
REAL :: FracOfCanopyAboveSnow(mp)

!Fraction Of Canopy Above Snow
FracOfCanopyAboveSnow = HgtAboveSnow/ MAX( 0.01, Hgt_PFT)

! LAI decreases due to snow:
reducedLAIdue2snow = lai_pft * FracOfCanopyAboveSnow
!* Leaf area, LAI, accounting for the presence of snow, is given by
!
! \[ LAI_{snow} = LAI_{nosnow} \frac{h_{c,snow}}{\max(0.01, h_{c,nosnow})} \]
!
! where \(h_{c,snow}\) is evaluated in subroutine [[HgtAboveSnow]].
! The LAI is decreased proportionally to the canopy height with/without snow.
!
! The effective canopy height used to evaluate the
! proportional decrease in LAI takes a minimum value of 0.01m.
! This condition prevents numerical issues when the canopy height
! (without snow) is small.
!
! **Warning: To follow up**
! 
! Can HgtAboveSnow be greater than MAX(0.01, Hgt_PFT)?
! If yes then there is an issue to resolve.
! HgtAboveSnow can sit at its minimum value of 10 \(z_{0,min}\) so if
! \(10 z_{0,min} > 0.01 > h_{c,nosnow} \) then LAI can increase if
! there is snow!  
!

RETURN
END SUBROUTINE LAI_eff

END MODULE hruff_eff_LAI_mod_cbl
