
MODULE cbl_sinbet_mod

  USE cable_common_module, ONLY : is_leapyear

  PUBLIC sinbet

contains

ELEMENTAL FUNCTION sinbet(year,doy,xslat,hod,jules_version) RESULT(z)

USE cable_math_constants_mod, ONLY : PI, PI180
! calculate sin(bet), bet = elevation angle of sun
! calculations according to goudriaan & van laar 1994 p30
INTEGER, INTENT(IN) ::                                                      &
  year       ! calendar year

REAL, INTENT(IN) ::                                                         &
  doy,     & ! day of year
  xslat,   & ! latitude (degrees north)
  hod        ! hour of day

LOGICAL, INTENT(IN) :: jules_version

REAL ::                                                                     &
  loy,     & ! length of year
  sindec,  & ! sine of maximum declination
  z          ! result

IF (jules_version) THEN
   loy = 365.0
   IF (is_leapyear(year)) loy = 366.0
   sindec = -SIN( 23.4 * PI180 * COS( 2. * PI * ( doy + 10.0 ) / loy ) ) 
ELSE
   sindec = -SIN( 23.45 * PI180 ) * COS( 2. * PI * ( doy + 10.0 ) / 365.0 )
END IF
   
z = MAX( SIN( PI180 * xslat ) * sindec                                 &
     + COS( PI180 * xslat ) * SQRT( 1. - sindec * sindec )              &
     * COS( PI * ( hod - 12.0 ) / 12.0 ), 1e-8 )

END FUNCTION sinbet

End MODULE cbl_sinbet_mod
