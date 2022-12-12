MODULE cbl_qsat_module

PUBLIC qsatfjh
PUBLIC qsatfjh2

CONTAINS

SUBROUTINE qsatfjh(mp, var, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, tair,pmb)

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
REAL, INTENT(IN) :: CRMH2O, CRMAIR
REAL, INTENT(IN) :: CTETENA, CTETENB, CTETENC
REAL, INTENT(IN) :: tair(mp)            ! air temperature (C)
REAL, INTENT(IN) :: pmb(mp)             ! pressure PMB (mb)

REAL, INTENT(OUT) :: var(mp)         ! result; sat sp humidity

!local vars
INTEGER :: j

DO j=1,mp

  var(j) = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair(j) /(CTETENC+tair(j))))    &
          / pmb(j)
ENDDO

RETURN
END SUBROUTINE qsatfjh

!==============================================================================

SUBROUTINE qsatfjh2( var, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, tair,pmb)
IMPLICIT NONE

REAL, INTENT(IN) :: CRMH2O, CRMAIR
REAL, INTENT(IN) :: CTETENA, CTETENB, CTETENC
REAL, INTENT(IN) :: tair            ! air temperature (C)
REAL, INTENT(IN) :: pmb             ! pressure PMB (mb)

REAL, INTENT(OUT) :: var             ! result; sat sp humidity

var = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair/(CTETENC+tair))) / pmb

RETURN
END SUBROUTINE qsatfjh2

END MODULE cbl_qsat_module
