MODULE cbl_rhoch_module

  IMPLICIT NONE

  PUBLIC calc_rhoch
  PRIVATE

CONTAINS

! this subroutine called from _init_radiation on cable_albedo.F90 pathway, explicit and implict
SUBROUTINE calc_rhoch( c1,rhoch, mp, nrb, taul, refl )

integer :: mp
integer :: nrb
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: taul(mp,nrb)
REAL :: refl(mp,nrb)

c1(:,1) = SQRT(1. - taul(:,1) - refl(:,1))
c1(:,2) = SQRT(1. - taul(:,2) - refl(:,2))
c1(:,3) = 1.

! Canopy C%REFLection black horiz leaves
! (eq. 6.19 in Goudriaan and van Laar, 1994):
rhoch = (1.0 - c1) / (1.0 + c1)

END SUBROUTINE calc_rhoch

END MODULE cbl_rhoch_module
