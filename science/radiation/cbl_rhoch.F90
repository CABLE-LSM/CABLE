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

MODULE cbl_rhoch_module

!-----------------------------------------------------------------------------
! Description:
!   Computes the reflection of black horizontal leaves
!   from Goudriaan and van Laar, 1994.
!
! This MODULE is USEd in:
!     cbl_init_radiation.F90 (JULES, CABLE, ESM1.5)
!
! This MODULE contains 1 public Subroutine:
!     calc_rhoch
!
! Module specific documentation: https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow: https://trac.nci.org.au/trac/cable/wiki/TBC
!-----------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC :: calc_rhoch
PRIVATE

CONTAINS

! this subroutine called from _init_radiation on cable_albedo.F90 pathway, explicit and implict
SUBROUTINE calc_rhoch( c1,rhoch, mp, nrb, taul, refl )
! Description:
!   Nothing further to add to the module description.

IMPLICIT NONE

!model dimensions
INTEGER, INTENT(IN) :: mp                     ! number of "tiles"
INTEGER, INTENT(IN) :: nrb                    ! # of rad. bands VIS,NIR(,LW)

REAL, INTENT(OUT) :: c1(mp,nrb)
REAL, INTENT(OUT) :: rhoch(mp,nrb)            ! REFLection black horiz leaves
REAL, INTENT(IN) :: taul(mp,nrb)              ! Leaf Transmisivity
REAL, INTENT(IN) :: refl(mp,nrb)              ! Leaf Reflectivity

c1(:,1) = SQRT(1.0 - taul(:,1) - refl(:,1))
c1(:,2) = SQRT(1.0 - taul(:,2) - refl(:,2))
c1(:,3) = 1.0

! Canopy C%REFLection black horiz leaves
! (eq. 6.19 in Goudriaan and van Laar, 1994):
rhoch = (1.0 - c1) / (1.0 + c1)

RETURN
END SUBROUTINE calc_rhoch

END MODULE cbl_rhoch_module
