
MODULE cbl_sinbet_mod

PUBLIC sinbet

contains

ELEMENTAL FUNCTION sinbet(doy,xslat,hod) RESULT(z)

USE cable_math_constants_mod, ONLY : PI, PI180
! calculate sin(bet), bet = elevation angle of sun
! calculations according to goudriaan & van laar 1994 p30
REAL, INTENT(IN) ::                                                         &
  doy,     & ! day of year
  xslat,   & ! latitude (degrees north)
  hod        ! hour of day

REAL ::                                                                     &
  sindec,  & ! sine of maximum declination
  z          ! result

sindec = -SIN( 23.45 * PI180 ) * COS( 2. * PI * ( doy + 10.0 ) / 365.0 )

z = MAX( SIN( PI180 * xslat ) * sindec                                 &
     + COS( PI180 * xslat ) * SQRT( 1. - sindec * sindec )              &
     * COS( PI * ( hod - 12.0 ) / 12.0 ), 1e-8 )

END FUNCTION sinbet

End MODULE cbl_sinbet_mod
