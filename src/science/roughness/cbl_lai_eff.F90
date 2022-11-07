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
  ! inputs:  <ul> <li> mp - number of land points</li>
  !           <li> LAI_PFT (mp) - leaf area with no snow (in m2/m2) </li>
  !           <li> HGT_pft(mp) - height of canopy with no snow (in m)</li>
  !           <li> HgtAboveSnow (mp) - height of canopy above snow surface (in m)</li></ul>
  !
  !  outputs: <ul> <li> reducedLAIdue2snow (m2/m2) - leaf area above snow (in m2/m2) </li></ul>
  !
  !           output variable formerly known as rough%hruff
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
  !* leaf area acconting for the presence of snow is given by
  !
  ! \( LAI_{snow} = LAI_{nosnow} * h_{c,snow} / \max[0.01, h_{c,nosnow} \)
  !
  ! where \(h_{c,snow}\) is evaluated in cbl_HgtAboveSnow.F90.
  !  The LAI is decreased proportionally to the canopy height with/without snow.
  !
  ! The effective canopy height takes a minimum value of 0.01m.
  !
  ! **There is potential for misbehaviour - if  \( h_{c,snow} < 0.01 \)
  ! then LAI could be increased (unphysical).  \(h_{c,snow} \ge 10 z{0,min}\)
  ! only which need not be greater than 0.01m.

END SUBROUTINE LAI_eff

End MODULE cbl_LAI_eff_mod
