MODULE cbl_LAI_eff_mod

   IMPLICIT NONE

   PUBLIC LAI_eff 

CONTAINS

!Computes Effective LAI of exposed canopy given effect of snow present
!variable formerly known as canopy%vlaiw
SUBROUTINE LAI_eff( mp, LAI_PFT, Hgt_PFT, HgtAboveSnow,  &
                    reducedLAIdue2snow ) 

  !* Subroutine computes the leaf-area index of a canopy when there is snow
  ! present.  This subroutine assumes that leaf area is distributed uniformly
  ! in the vertical.
  !
  ! inputs:
  !
  ! * mp - number of land points</li>
  ! * LAI_PFT (mp) - leaf area with no snow (in m2/m2)
  ! * HGT_pft(mp) - height of canopy with no snow (in m)
  ! * HgtAboveSnow (mp) - height of canopy above snow surface (in m) 
  !
  ! outputs:
  !
  ! * reducedLAIdue2snow (m2/m2) - leaf area with snow (in m2/m2)
  !
  ! The output variable was formerly known as rough%vlaiw
  !

  !re-decl input args  
  integer  :: mp
  real :: LAI_PFT(mp)
  real :: Hgt_PFT(mp)
  real :: HgtAboveSnow(mp) 
  real :: reducedLAIdue2snow(mp) 
 
  !local_vars: 
  real :: Hgt_eff(mp)
  real :: FracOfCanopyAboveSnow(mp)

  !Fraction Of Canopy Above Snow
  FracOfCanopyAboveSnow = HgtAboveSnow/ MAX( 0.01, Hgt_PFT)
  
  ! LAI decreases due to snow:
  reducedLAIdue2snow = LAI_PFT * FracOfCanopyAboveSnow
  !* Leaf area, accounting for the presence of snow, is given by
  !
  ! \( LAI_{snow} = LAI_{nosnow} h_{c,snow} / \max[0.01, h_{c,nosnow} ] \)
  !
  ! where \(h_{c,snow}\) is evaluated in cbl_HgtAboveSnow.F90.
  ! The LAI is decreased proportionally to the canopy height with/without snow.
  !
  ! The effective canopy height used to evaluate the
  ! proportional decrease in LAI takes a minimum value of 0.01m.
  ! This conditions prevents numerical issues when the canopy height
  ! (without snow) is small.
  !
  ! **To follow up** Can HgtAboveSnow be greater than MAX(0.01, Hgt_PFT)?
  ! HgtAboveSnow takes a minimum value of 10 \(z_{0,min}\) so if
  ! 10 \(z_{0,min} > 0.01 > h_{c,nosnow} \) then LAI is increased when
  ! there is snow!  
  !
  

END SUBROUTINE LAI_eff

End MODULE cbl_LAI_eff_mod
