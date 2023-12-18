MODULE cbl_qsat_module
  !* This MODULE contains two SUBROUTINEs that calculate
  ! the specific humidity at saturation as a function of air pressure
  ! and temperature.
  !
  ! The two SUBROUTINEs differ only in that [[qsatfjh]] operates on an array
  ! whereas [[qsatfjh2]] operates on single element REALs.
  !
  ! **Warning:** [[qsatfjh2]] is redundant and the code should be changed to only use [[qsatfjh]].

PUBLIC qsatfjh
PUBLIC qsatfjh2

CONTAINS

SUBROUTINE qsatfjh(mp, var, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, tair,pmb)
    !*## Purpose
    !
    ! This SUBROUTINE evaluates the specific humidity (water vapour mixing ratio)
    ! at saturation (in kgkg\(^{-1}\)), at a given air temperature (in \(^{\circ}\)C)
    ! and air pressure (in hPa) for an array of size mp.
    !
    !## Method
    !
    ! The Teten's formula for specific humidity at saturation is used,
    ! based on the values for the mass of a mole of water,`CRMH2o`, 
    ! mass of a mole of dry air `Crmair`, and the Teten constants.
    !
    !## Reference
    !
    ! [Murray F. W., 1967](https://doi.org/10.1175/1520-0450(1967)006%3C0203:OTCOSV%3E2.0.CO;2)

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp               !! size of array of land points (-)
REAL, INTENT(IN) :: CRMH2O, CRMAIR
REAL, INTENT(IN) :: CTETENA, CTETENB, CTETENC
REAL, INTENT(IN) :: tair(mp)            !! air temperature (\(^{\circ}\)C)
REAL, INTENT(IN) :: pmb(mp)             ! surface air pressure (hPa)

REAL, INTENT(OUT) :: var(mp)            !! specific humidity at saturation (kgkg\(^{-1}\))

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
  !*## Purpose
  !
  ! This SUBROUTINE evaluates the specific humidity (water vapour mixing ratio)
  ! at saturation (in kgkg\(^{-1}\)), at a given air temperature (in \(^{\circ}\)C)
  ! and air pressure (in hPa).
  !
  !## Method
  !
  ! The Teten's formula for specific humidity at saturation is used,
  ! based on the values for the mass of a mole of water,`CRMH2o`, 
  ! mass of a mole of dry air `Crmair`, and the Teten constants.
  !
  !## Reference
  !
  ! [Murray F. W., 1967](https://doi.org/10.1175/1520-0450(1967)006%3C0203:OTCOSV%3E2.0.CO;2)
IMPLICIT NONE

REAL, INTENT(IN) :: CRMH2O, CRMAIR
REAL, INTENT(IN) :: CTETENA, CTETENB, CTETENC
REAL, INTENT(IN) :: tair            !! air temperature (\(^{\circ}\)C),
REAL, INTENT(IN) :: pmb             !! surface air pressure (hPa)

REAL, INTENT(OUT) :: var            !! specific humidity at saturation (kgkg\(^{-1}\))

var = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair/(CTETENC+tair))) / pmb

RETURN
END SUBROUTINE qsatfjh2

END MODULE cbl_qsat_module
