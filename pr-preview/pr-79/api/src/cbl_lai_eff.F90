MODULE cbl_LAI_eff_mod

   !* The procedure in this module computes the effective LAI of a canopy
   !  given the effect of any snow present

   IMPLICIT NONE

   PUBLIC LAI_eff 

CONTAINS


SUBROUTINE LAI_eff( mp, LAI_PFT, Hgt_PFT, HgtAboveSnow,  &
                    reducedLAIdue2snow ) 

  !* This subroutine computes the leaf-area index of a canopy when there is snow
  ! present. The formulae assume that the leaf area is distributed uniformly
  ! in the vertical within the canopy.

  !re-decl input args  
  integer  :: mp !! Number of land points (-)
  real :: LAI_PFT(mp) !! Leaf area index without snow (m\(^2\) m\(^{-2})\)
  real :: Hgt_PFT(mp) !! height of canopy without snow (m)
  real :: HgtAboveSnow(mp) !! height of canopy above the snow surface (m)
  real :: reducedLAIdue2snow(mp) 
   !! Output. Modified leaf area index for snow (m\(^2\) m\(^{-2}\))
 
  !local_vars: 
  real :: Hgt_eff(mp)
  real :: FracOfCanopyAboveSnow(mp)

  !Fraction Of Canopy Above Snow
  FracOfCanopyAboveSnow = HgtAboveSnow/ MAX( 0.01, Hgt_PFT)
  
  ! LAI decreases due to snow:
  reducedLAIdue2snow = LAI_PFT * FracOfCanopyAboveSnow
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
  

END SUBROUTINE LAI_eff

End MODULE cbl_LAI_eff_mod
