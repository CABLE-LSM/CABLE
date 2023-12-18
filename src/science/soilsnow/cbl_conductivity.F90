
MODULE total_soil_conductivity_mod

USE cbl_ssnow_data_mod

PUBLIC  total_soil_conductivity

CONTAINS

! soil thermal conductivity (incl water/ice)
FUNCTION total_soil_conductivity(ssnow,soil)

    REAL(r_2), DIMENSION(mp,ms) ::  total_soil_conductivity

    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL(r_2) :: exp_arg
    REAL(r_2) :: dels_r2
    REAL(r_2) :: Ko,Ktmp
    REAL(r_2), DIMENSION(mp,ms) :: Ke,quartz,Sr,Ksat,liq_frac
    REAL      :: tfreeze
    INTEGER :: k,j,i

    total_soil_conductivity(:,:) = soil%cnsd_vec(:,:)

    DO k = 1, ms
       DO j = 1, mp
          IF (soil%isoilm(j) .EQ. 9) THEN
             total_soil_conductivity(j,k) = snow_ccnsw
          ELSE
             quartz(j,k) = MAX(0.0,MIN(0.8,soil%sand_vec(j,k)*0.92))
             IF (quartz(j,k) .GT. 0.2) THEN
                Ko = 2.0
             ELSE
                Ko = 3.0
             END IF

             Ktmp      = ( (7.7**(quartz(j,k))) * &
                  (Ko**(1.0-quartz(j,k))) ) **(1.0-soil%ssat_vec(j,k))

             IF (ssnow%wb(j,k) .GE. 1.0e-15) THEN
                liq_frac(j,k) = MIN(1._r_2, MAX(0._r_2, ssnow%wbliq(j,k) / ssnow%wb(j,k)))
             ELSE
                liq_frac(j,k) = 0.0
             END IF

             Ksat(j,k) =  Ktmp * &
                  (2.2 ** (soil%ssat_vec(j,k)*(1.0-liq_frac(j,k) ) ) )*&
                  (0.57**(liq_frac(j,k)))

             Sr(j,k) = MIN( 0.9999 , &
                  MAX(0., ssnow%wb(j,k)-soil%watr(j,k))/(soil%ssat_vec(j,k)-soil%watr(j,k)) )

             !frozen or not?
             IF (Sr(j,k) .GE. 0.05) THEN
                Ke(j,k) = 0.7*LOG10(Sr(j,k)) + 1.0
             ELSE
                Ke(j,k) = 0.0
             END IF

             IF ((ssnow%wbice(j,k) .GT. 0.0) .OR. &
                  (ssnow%tgg(j,k) .LT. CTFRZ) .OR. &
                  (ssnow%isflag(j) .NE. 0) .OR.     &
                  (ssnow%snowd(j) .GE. 0.1) )   THEN

                Ke(j,k) = Sr(j,k)

             END IF

             total_soil_conductivity(j,k) = Ke(j,k)*Ksat(j,k) + &
                  (1.0-Ke(j,k))*soil%cnsd_vec(j,k)

             total_soil_conductivity(j,k) = MIN(Ksat(j,k), MAX(soil%cnsd_vec(j,k),&
                  total_soil_conductivity(j,k) ) )


          ENDIF

       END DO

    END DO

  END FUNCTION total_soil_conductivity

END MODULE total_soil_conductivity_mod


