MODULE cable_psm

  USE cable_def_types_mod, ONLY : r_2,ms,mp,air_type,met_type,soil_snow_type,&
       canopy_type,soil_parameter_type,veg_parameter_type,&
       roughness_type
  USE cable_common_module, ONLY : cable_user

  IMPLICIT NONE


  REAL(r_2), PARAMETER ::rt_Dff=2.5e-5, & !diffusivity in air
       lm=1.73e-5, &       !converts units
       c2 = 2.0,&                  !params
       litter_thermal_diff=2.7e-5  !param based on vh thermal diffusivity

  REAL(r_2), PARAMETER :: rtevap_max = 10000.0
  REAL(r_2), DIMENSION(0:8), PARAMETER :: gamma_pre = &
       (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
       771.32342877765313, -176.61502916214059, 12.507343278686905, &
       -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)
  INTEGER, PARAMETER                      :: c_gamma = 7
  REAL(r_2),PARAMETER :: pi_r_2=3.14159


  PUBLIC  or_soil_evap_resistance,update_or_soil_resis,rtevap_max,rt_Dff

CONTAINS

  RECURSIVE FUNCTION my_gamma(a) RESULT(g)


    REAL(r_2), INTENT(in) :: a
    REAL(r_2) :: g

    ! these precomputed values are taken by the sample code in Wikipedia,
    ! and the sample itself takes them from the GNU Scientific Library
    !real(r_2), dimension(0:8), parameter :: gamma_pre = &
    !     (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
    !     771.32342877765313, -176.61502916214059, 12.507343278686905, &
    !     -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)

    REAL(r_2) :: t, w, x
    INTEGER :: i

    x = a

    IF ( x < 0.5 ) THEN
       g = (pi_r_2) / ( SIN((pi_r_2)*x) * my_gamma(1.0-x) )
    ELSE
       x = x - 1.0
       t = gamma_pre(0)
       DO i=1, c_gamma+2
          t = t + gamma_pre(i-1)/(x+REAL(i,r_2))
       END DO
       w = x + REAL(c_gamma,r_2) + 0.5
       g = SQRT(2.0*(pi_r_2)) * w**(x+0.5) * EXP(-w) * t
    END IF
  END FUNCTION my_gamma

  SUBROUTINE or_soil_evap_resistance(soil,air,met,canopy,ssnow,veg,rough)

    TYPE (air_type), INTENT(INOUT)       :: air
    TYPE (met_type), INTENT(INOUT)       :: met
    TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT) :: veg
    TYPE (roughness_type), INTENT(INOUT) :: rough



    REAL(r_2), DIMENSION(mp) :: sublayer_dz, eddy_shape,eddy_mod,soil_moisture_mod, &
         soil_moisture_mod_sat, wb_liq, &
         pore_size,pore_radius, rel_s,hk_zero,hk_zero_sat,time_scale  !note pore_size in m

    REAL(r_2), DIMENSION(mp) :: litter_dz

    INTEGER, DIMENSION(mp) :: int_eddy_shape

    !

    INTEGER :: i,j,k

    IF (cable_user%litter) THEN
       litter_dz(:) = veg%clitt*0.003
    ELSE
       litter_dz(:) = 0.0
    ENDIF

    pore_radius(:) = 0.148  / (1000.0*9.81*ABS(soil%sucs_vec(:,1))/1000.0)  !should replace 0.148 with surface tension, unit coversion, and angle
    pore_size(:) = pore_radius(:)*SQRT((pi_r_2))

    !scale ustar according to the exponential wind profile, assuming we are a mm from the surface
    eddy_shape = 0.3*met%ua/ MAX(1.0e-4,MAX(1.0e-3,canopy%us)*&
         EXP(-rough%coexp*(1.0-canopy%sublayer_dz/MAX(1e-2,rough%hruff))))
    int_eddy_shape = FLOOR(eddy_shape)
    eddy_mod(:) = 0.0
    DO i=1,mp
       IF (veg%iveg(i) .LT. 16) THEN
          eddy_mod(i) = 2.2*SQRT(112.0*(pi_r_2)) / (2.0**(eddy_shape(i)+1.0) * SQRT(eddy_shape(i)+1.0))

          IF (int_eddy_shape(i) .GT. 0) THEN
             eddy_mod(i) = eddy_mod(i) / my_gamma(eddy_shape(i)+1.0) * (2.0*eddy_shape(i)+1.0)
             DO k=1,int_eddy_shape(i)
                eddy_mod(i) = eddy_mod(i) * (2.0*(eddy_shape(i) - k) + 1.0)
             END DO
          END IF
          canopy%sublayer_dz(i) = MIN(0.05, MAX(eddy_mod(i) * air%visc(i) / MAX(1.0e-3,canopy%us(i))*&
               EXP(-rough%coexp(i)*(1.0-canopy%sublayer_dz(i)/MAX(1e-2,rough%hruff(i)))),1e-7)  )

       ELSE

          canopy%sublayer_dz(i) = 0.0

       END IF

    END DO

    DO i=1,mp
       IF (veg%iveg(i) .LT. 16) THEN

          wb_liq(i) = REAL(MAX(0.0001,MIN((pi_r_2)/4.0, &
               (ssnow%wb(i,1)-ssnow%wbice(i,1) - ssnow%satfrac(i)*soil%ssat_vec(i,1)) / &
               MAX((1._r_2 - ssnow%satfrac(i)),1e-5) ) ) )

          rel_s(i) = REAL( MAX(wb_liq(i)-soil%watr(i,1),0._r_2)/(soil%ssat_vec(i,1)-soil%watr(i,1)) )
          hk_zero(i) = MAX(0.001*soil%hyds_vec(i,1)*(MIN(MAX(rel_s(i),0.001_r_2),1._r_2)**(2._r_2*soil%bch_vec(i,1)+3._r_2) ),1e-8)
          hk_zero_sat(i) = MAX(0.001*soil%hyds_vec(i,1),1e-8)

          soil_moisture_mod(i)     = 1.0/(pi_r_2)/SQRT(wb_liq(i))* ( SQRT((pi_r_2)/(4.0*wb_liq(i)))-1.0)
          soil_moisture_mod_sat(i) = 1.0/(pi_r_2)/SQRT(soil%ssat_vec(i,1))* ( SQRT((pi_r_2)/(4.0*soil%ssat_vec(i,1)))-1.0)

          IF (ssnow%isflag(i) .EQ. 0 .AND. (ssnow%snowd(i) .LT. 1.0e-4) ) THEN

             canopy%sublayer_dz(i) = canopy%sublayer_dz(i) + litter_dz(i)

          ELSE

             canopy%sublayer_dz(i) = 0.0
             hk_zero(i)     = 1.0e15
             hk_zero_sat(i) = 1.0e15
             soil_moisture_mod(i)     = 0.0
             soil_moisture_mod_sat(i) = 0.0

          END IF

          IF (canopy%sublayer_dz(i) .GE. 1.0e-7) THEN
             ssnow%rtevap_unsat(i) = MIN(rtevap_max, &
                  rough%z0soil(i)/canopy%sublayer_dz(i) * (lm/ (4.0*hk_zero(i)) +&
                  (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod(i)) / rt_Dff))
             ssnow%rtevap_sat(i)  = MIN(rtevap_max, &
                  rough%z0soil(i)/canopy%sublayer_dz(i) * (lm/ (4.0*hk_zero_sat(i)) + &
                  (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod_sat(i)) / rt_Dff))

             ssnow%rt_qh_sublayer(i) = canopy%sublayer_dz(i) / litter_thermal_diff

          ELSE
             ssnow%rtevap_unsat(i) = MIN(rtevap_max, &
                  lm/ (4.0*hk_zero(i)) + (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod(i)) / rt_Dff)
             ssnow%rtevap_sat(i)  = MIN(rtevap_max, &
                  lm/ (4.0*hk_zero_sat(i)) + (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod_sat(i)) / rt_Dff)

             ssnow%rt_qh_sublayer(i) = 0.0
          END IF

       ELSE
          !no additional evap resistane over lakes
          ssnow%rtevap_unsat(i) = 0.0
          ssnow%rt_qh_sublayer(i) = 0.0
          ssnow%satfrac(i) = 0.5
          IF (veg%iveg(i) .EQ. 16 .AND. met%tk(i) .LT. 268.15 ) &
               ssnow%rtevap_sat(i) = 0.41*ssnow%rtsoil(i)

       END IF

    END DO


  END SUBROUTINE or_soil_evap_resistance


  SUBROUTINE update_or_soil_resis(ssnow,canopy,veg,dq,dqu)

    TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (veg_parameter_type), INTENT(INOUT) :: veg
    REAL, DIMENSION(mp), INTENT(IN) :: dq,&
         dqu

    INTEGER :: i


    DO i=1,mp

       IF (veg%iveg(i) .LT. 16 .AND. ssnow%snowd(i) .LT. 1e-7) THEN

          IF (dq(i) .LE. 0.0) THEN
             ssnow%rtevap_sat(i) = MIN(rtevap_max,canopy%sublayer_dz(i)/rt_Dff)
          END IF

          IF (dqu(i) .LE. 0.0) THEN
             ssnow%rtevap_unsat(i) = MIN(rtevap_max,canopy%sublayer_dz(i)/rt_Dff)
          END IF

       END IF
    END DO






  END SUBROUTINE update_or_soil_resis

END MODULE cable_psm
