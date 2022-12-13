MODULE cbl_qsat_module
  !* This MODULE contains two SUBROUTINEs that calculate
  ! the specific humdity at saturation as a function of air pressure
  ! and temperature.
  !
  ! The two SUBROUTINEs differ only in that [[qsatfjh]] operates on an array
  ! whereas as [[qsatjfh2]] operates on single element REALs.

PUBLIC qsatfjh
PUBLIC qsatfjh2

CONTAINS

  SUBROUTINE qsatfjh(mp, var, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, tair,pmb)
    !*## Purpose
    !
    ! This SUBROUTINE evaluates the specific humdity (water vapour mixing ratio)
    ! at saturation (in kgkg\(^{-1}\)), at air temperature (in \(^{\circ}\)C)
    ! and air pressure (in hPa) across an array of size mp.
    !
    !## Method
    !
    ! The standard formula for specific humdity at saturation is used,
    ! see e.g. Garratt (1995), based on the values for the mass of a mole
    ! of water,`CRMH2o`, mass of a mole of dry air `Crmair`,
    ! and the Teten constants.

IMPLICIT  NONE   
integer :: mp                        !! size of array of land points
REAL, INTENT(OUT) :: var(mp)         !! sat specific humidity (kgkg\(^{-1}\))
REAL :: CRMH2o, Crmair, CTETENA, CTETENB, CTETENC
REAL, INTENT(IN) ::                                          &
  tair(mp),                        & !! air temperature (\(^{\circ}\)C)
  pmb(mp)                            !! surface air pressure (hPa)

INTEGER :: j

DO j=1,mp

  var(j) = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair(j)/(CTETENC+tair(j))))  &
          / pmb(j)
ENDDO

END SUBROUTINE qsatfjh



SUBROUTINE qsatfjh2( var, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, tair,pmb)
  !*## Purpose
  !
  ! This SUBROUTINE evaluates the specific humdity (water vapour mixing ratio)
  ! at saturation (in kgkg\(^{-1}\)), at air temperature (in \(^{\circ}\)C) and
  ! and air pressure (in hPa).
  !
  !## Method
  !
  ! The standard formula for specific humdity at saturation is used,
  ! see e.g. Garratt (1995), based on the values for the mass of a mole
  ! of water,`CRMH2o`, mass of a mole of dry air `Crmair`,
  ! and the Teten constants.
  
REAL :: CRMH2o, Crmair, CTETENA, CTETENB, CTETENC
  REAL, INTENT(IN) ::                                                         &
       tair,         & !! air temperature (\(^{\circ}\)C),
       pmb             !! surface air pressure (hPa)

  REAL, INTENT(OUT) ::                                                        &
       var             !! specific humidity at saturation (kgkg\(^{-1}\))

  var = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair/(CTETENC+tair))) / pmb

END SUBROUTINE qsatfjh2

END MODULE cbl_qsat_module
