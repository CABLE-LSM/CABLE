MODULE cbl_pot_evap_snow_module

IMPLICIT NONE

PUBLIC Penman_Monteith
PUBLIC Humidity_deficit_method
PRIVATE

CONTAINS

FUNCTION Penman_Monteith( mp, Ctfrz, CRMH2o, Crmair, CTETENA, CTETENB,         &
                          CTETENC, veg_clitt, cable_user_litter,               &
                          air_dsatdk, air_psyc, air_rho, air_rlam,             & 
                          met_tvair, met_pmb, met_qvair,                       &
                          ground_H_flux, canopy_fns, canopy_DvLitt,            & 
                          ssnow_rtsoil, ssnow_isflag )        RESULT(ssnowpotev)
USE cbl_qsat_module, ONLY : qsatfjh 
IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
!REAL, INTENT(OUT)   :: ssnowpotev(mp)     ! returned result of function
REAL   :: ssnowpotev(mp)     ! returned result of function

REAL, INTENT(IN) :: Ctfrz 
REAl :: CRMH2o
REAl :: Crmair
REAl :: CTETENA, CTETENB, CTETENC
REAL, INTENT(IN) :: veg_clitt(mp)
LOGICAL, INTENT(IN) :: cable_user_litter
REAL, INTENT(IN) :: air_dsatdk(mp)
REAL, INTENT(IN) :: air_psyc(mp)
REAL, INTENT(IN) :: air_rho(mp)
REAL, INTENT(IN) :: air_rlam(mp)
REAL, INTENT(IN) :: met_tvair(mp) 
REAL, INTENT(IN) :: met_qvair(mp) 
REAL, INTENT(IN) :: met_pmb(mp)
REAL, INTENT(IN) :: ground_H_flux(mp)
REAL, INTENT(IN) :: canopy_fns(mp)
REAL, INTENT(IN) :: canopy_DVLitt(mp)
REAL, INTENT(IN) :: ssnow_rtsoil(mp)
INTEGER, INTENT(IN) :: ssnow_isflag(mp)

!local vars
REAL :: sss(mp)          ! var for Penman-Monteith soil evap  
REAL :: cc1(mp)          ! var for Penman-Monteith soil evap                                                    &
REAL :: cc2(mp)          ! var for Penman-Monteith soil evap
REAL :: qsatfvar(mp)
INTEGER :: j

! Penman-Monteith formula
sss=air_dsatdk
cc1=sss/(sss+air_psyc )
cc2=air_psyc /(sss+air_psyc )

CALL qsatfjh( mp, qsatfvar, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC,         &
              met_tvair-CTfrz,met_pmb)

!INH 10-1-2017 - this P-M implementation is incorrect over snow.
!variable ssnowpotev is actually the latent heat flux associated with
!potential evaporation.
!Needs to be addressed/simplified at a later date - involves changes
!to HDM method and latent_heat_flux() and elsewhere

IF (cable_user_litter) THEN
  ! vh_js !
  ssnowpotev = cc1 * (canopy_fns - ground_H_flux) + &
               cc2 * air_rho * air_rlam*(qsatfvar - met_qvair)/ &
               (ssnow_rtsoil+ REAL((1-ssnow_isflag))*veg_clitt*0.003/canopy_DvLitt)
ELSE
  ssnowpotev = cc1 * (canopy_fns - ground_H_flux) + &
               cc2 * air_rho * air_rlam*(qsatfvar  - met_qvair)/ssnow_rtsoil
ENDIF

RETURN
END FUNCTION Penman_Monteith



FUNCTION Humidity_deficit_method( mp, Ctfrz, veg_clitt,cable_user_or_evap,     &
                                 cable_user_gw_model, cable_user_litter,       &
                                 air_rho,air_rlam,           & 
                                 dq,dqu,qstss,   & 
                                 canopy_DvLitt,      &
                                 ssnow_isflag, ssnow_satfrac, ssnow_rtsoil,    &
                                 ssnow_rtevap_sat, ssnow_rtevap_unsat,      & 
                                 ssnow_snowd, ssnow_tgg &
                                 ) RESULT(ssnowpotev)
IMPLICIT NONE

INTEGER :: mp
REAL ::  ssnowpotev(mp)

REAL:: Ctfrz 
REAL, INTENT(IN) :: veg_clitt(mp)
LOGICAL :: cable_user_or_evap, cable_user_gw_model, cable_user_litter
REAL :: air_rho(mp)    !
REAL :: air_rlam(mp)    !
REAL :: dq(mp)       ! sat spec hum diff.
REAL :: dqu(mp)      ! sat spec hum diff.
REAL :: qstss(mp)    !dummy var for compilation
REAL :: canopy_DvLitt(mp)    ! 
INTEGER, INTENT(IN) :: ssnow_isflag(mp)
REAL :: ssnow_snowd(mp)    ! 
REAL :: ssnow_tgg(mp)    ! 
REAL :: ssnow_satfrac(mp)    !
REAL :: ssnow_rtsoil(mp)    !
REAL :: ssnow_rtevap_sat(mp)    !
REAL :: ssnow_rtevap_unsat(mp)    !

!local vars
INTEGER :: j
REAL, DIMENSION(mp) :: q_air

q_air = qstss - dq

DO j=1,mp

  IF( ssnow_snowd(j)>1.0 .OR. ssnow_tgg(j) .EQ. Ctfrz ) THEN
    dq(j) = MAX( -0.1e-3, dq(j))
    dqu(j) = MAX( -0.1e-3, dqu(j))
  END IF

  IF (dq(j) .LE. 0.0 .AND. dqu(j) .LT. dq(j)) THEN
    dqu(j) = dq(j)
  END IF

  IF (dq(j) .GE. 0.0 .AND. dqu(j) .LT. 0.0) THEN
    dqu(j) = 0.0
  ENDIF
ENDDO

IF (cable_user_or_evap .or. cable_user_gw_model) then
  write(6,*) "GW or ORevepis not an option right now"
  !H!        IF (cable_user_or_evap) THEN
  !H!          do j=1,mp
  !H!       
  !H!             if (veg_iveg(j) .lt. 16 .and. ssnow_snowd(j) .lt. 1e-7) THEN
  !H!       
  !H!                if (dq(j) .le. 0.0) THEN
  !H!                   ssnow_rtevap_sat(j) = min(rtevap_max,canopy_sublayer_dz(j)/rt_Dff)
  !H!                end if
  !H!       
  !H!                if (dqu(j) .le. 0.0) THEN
  !H!                   ssnow_rtevap_unsat(j) = min(rtevap_max,canopy_sublayer_dz(j)/rt_Dff)
  !H!                end if
  !H!       
  !H!             end if
  !H!
  !H!          end do
  !H!
  !H!        END IF

  ssnowpotev = air_rho * air_rlam * ( &
               REAL(ssnow_satfrac) * dq /(ssnow_rtsoil + REAL(ssnow_rtevap_sat)) + &
               (1.0 - REAL(ssnow_satfrac))* dqu/( &
               ssnow_rtsoil + REAL(ssnow_rtevap_unsat)) )

 ELSEIF (cable_user_litter) THEN
         ! vh_js !
  ssnowpotev = air_rho * air_rlam * dq /( ssnow_rtsoil +                       &
                          REAL((1-ssnow_isflag))* veg_clitt*0.003/canopy_DvLitt)
 ELSE
  ssnowpotev = air_rho * air_rlam * dq / ssnow_rtsoil
 ENDIF

RETURN
END FUNCTION Humidity_deficit_method

END MODULE cbl_pot_evap_snow_module
