MODULE cbl_LAI_eff_mod

   IMPLICIT NONE

   PUBLIC LAI_eff 

CONTAINS

!Computes Effective LAI of exposed canopy given effect of snow present
!variable formerly known as canopy%vlaiw
SUBROUTINE LAI_eff( mp, LAI_PFT, Hgt_PFT, HgtAboveSnow,  &
                    reducedLAIdue2snow ) 

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

END SUBROUTINE LAI_eff

End MODULE cbl_LAI_eff_mod
