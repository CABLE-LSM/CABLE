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
MODULE cbl_masks_mod

!-----------------------------------------------------------------------------
! Description:
!   Computes various masks for CABLE:
!     - veg_mask: TRUE where LAI > LAI threshold value
!     - sunlit_mask: TRUE where the cosine of zenith angle is higher than a
!                    specified tolerance.
!     - sunlit_veg_mask: TRUE where veg_mask and sunlit_mask are TRUE.
!
! This MODULE is USEd in:
!     cable_land_albedo_mod_cbl.F90 (JULES),
!     cable_cbm.F90 (ESM1.5),
!     cable_rad_driver.F90 (ESM1.5),
!     cbl_model_driver_offline.F90 (CABLE)
!
! This MODULE contains 3 public Subroutine:
!     fveg_mask,
!     fsunlit_mask,
!     fsunlit_veg_mask
!
! Module specific documentation: https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow: https://trac.nci.org.au/trac/cable/wiki/TBC
!-----------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC fveg_mask
PUBLIC fsunlit_mask
PUBLIC fsunlit_veg_mask

CONTAINS

SUBROUTINE fveg_mask( veg_mask, mp, lai_thresh, reducedLAIdue2snow )
! Description:
!   Computes mask where LAI is higher than a given threshold.

IMPLICIT NONE
LOGICAL, INTENT(OUT) :: veg_mask(:)

INTEGER, INTENT(IN) :: mp
REAL,    INTENT(IN) :: lai_thresh
REAL,    INTENT(IN) :: reducedLAIdue2snow(mp)
!local vars
INTEGER :: i

veg_mask(:) = .FALSE.
! Define vegetation mask:
DO i=1, mp
  IF ( reducedLAIdue2snow(i) > LAI_thresh ) veg_mask(i) = .TRUE.
END DO

RETURN
END SUBROUTINE fveg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE fsunlit_mask( sunlit_mask, mp, coszen_tols, coszen )
! Description:
!   Computes mask where the sun is above the horizon within a given tolerance.

IMPLICIT NONE
LOGICAL, INTENT(OUT) :: sunlit_mask(:)

INTEGER, INTENT(IN) :: mp
REAL,    INTENT(IN) :: coszen_tols
REAL,    INTENT(IN) :: coszen(mp)
!local vars
INTEGER :: i

sunlit_mask = .FALSE.
! Define sunlit mask:
DO i=1, mp
  IF ( coszen(i) > coszen_tols ) sunlit_mask(i) = .TRUE.
END DO

RETURN
END SUBROUTINE fsunlit_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE fsunlit_veg_mask(  sunlit_veg_mask, veg_mask, sunlit_mask, mp )
! Description:
!   Computes mask where the sun is above the horizon within a given tolerance
!   and the LAI is above a given threshold. This is the union of the
!   masks from fsunlit_mask() and fveg_mask().

IMPLICIT NONE
LOGICAL, INTENT(OUT) :: sunlit_veg_mask(:)

INTEGER, INTENT(IN) :: mp
LOGICAL, INTENT(IN) :: veg_mask(mp)
LOGICAL, INTENT(IN) :: sunlit_mask(mp)
!local vars
INTEGER :: i

sunlit_veg_mask = .FALSE.
! Define sunlit AND vegetation mask:
DO i=1, mp
  IF ( veg_mask(i) .AND.  sunlit_mask(i) ) sunlit_veg_mask(i) = .TRUE.
END DO

RETURN
END SUBROUTINE fsunlit_veg_mask

END MODULE cbl_masks_mod
