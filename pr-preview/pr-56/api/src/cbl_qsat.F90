MODULE cbl_qsat_module

PUBLIC qsatfjh
PUBLIC qsatfjh2

CONTAINS

SUBROUTINE qsatfjh(mp, var, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, tair,pmb)

IMPLICIT  NONE   
integer :: mp 
REAL, INTENT(OUT) :: var(mp)         ! result; sat sp humidity
REAL :: CRMH2o, Crmair, CTETENA, CTETENB, CTETENC
REAL, INTENT(IN) ::                                          &
  tair(mp),                        & ! air temperature (C)
  pmb(mp)                            ! pressure PMB (mb)

INTEGER :: j

DO j=1,mp

  var(j) = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair(j)/(CTETENC+tair(j))))  &
          / pmb(j)
ENDDO

END SUBROUTINE qsatfjh



SUBROUTINE qsatfjh2( var, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, tair,pmb)

REAL :: CRMH2o, Crmair, CTETENA, CTETENB, CTETENC
  REAL, INTENT(IN) ::                                                         &
       tair,         & ! air temperature (C)
       pmb             ! pressure PMB (mb)

  REAL, INTENT(OUT) ::                                                        &
       var             ! result; sat sp humidity

  var = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair/(CTETENC+tair))) / pmb

END SUBROUTINE qsatfjh2

END MODULE cbl_qsat_module
