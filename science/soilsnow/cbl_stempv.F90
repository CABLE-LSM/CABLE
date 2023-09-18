MODULE stempv_mod

USE cbl_ssnow_data_mod

PUBLIC  stempv

CONTAINS 

! calculates temperatures of the soil
! tgg - new soil/snow temperature
! ga - heat flux from the atmosphere (ground heat flux)
! ccnsw - soil thermal conductivity, including water/ice
SUBROUTINE stempv(dels, canopy, ssnow, soil,heat_cap_lower_limit)

USE cable_common_module, ONLY: cable_user
USE trimb_mod,                   ONLY: trimb
USE total_soil_conductivity_mod, ONLY: total_soil_conductivity 
USE old_soil_conductivity_mod,   ONLY: old_soil_conductivity 
IMPLICIT NONE
    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(canopy_type),    INTENT(INOUT) :: canopy
    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL, DIMENSION(mp) ::                                                      &
         coefa, coefb,  & !
         sgamm            !

    REAL(r_2), DIMENSION(mp) ::                                                 &
         dtg,     & !
         ew,      & !
         xx,      & !
         wblfsp     !

    REAL(r_2), DIMENSION(mp,ms) ::                                              &
         ccnsw  ! soil thermal conductivity (incl water/ice)

    REAL(r_2), DIMENSION(mp, -2:ms) ::                                          &
         at, bt, ct, rhs !

    REAL(r_2), DIMENSION(mp,-2:ms+1) :: coeff

    REAL(r_2), DIMENSION(mp,ms+3)    :: tmp_mat ! temp. matrix for tggsn & tgg

    INTEGER :: j,k
    REAL :: exp_arg
    LOGICAL :: direct2min = .FALSE.
    REAL :: heat_cap_lower_limit(mp,ms)    ! best to declare INTENT here - rk4417 - phase2

    at = 0.0
    bt = 1.0
    ct = 0.0
    coeff = 0.0

    ssnow%otgg(:,:) = ssnow%tgg  ! FEEDBACK (OK to insert this line as per MMY code?) --rk4417 
    
    IF (cable_user%soil_thermal_fix) THEN
       ccnsw = total_soil_conductivity(ssnow,soil)
    ELSE
       ccnsw = old_soil_conductivity(ssnow,soil)
    ENDIF

    xx = 0.

    WHERE(ssnow%isflag == 0)
       xx = MAX( 0., ssnow%snowd / ssnow%ssdnn )
       ccnsw(:,1) = ( ccnsw(:,1) - 0.2 ) * ( soil%zse(1) / ( soil%zse(1) + xx ) &
            ) + 0.2
    END WHERE

    DO k = 3, ms

       WHERE (ssnow%isflag == 0)
          coeff(:,k) = 2.0 / ( soil%zse(k-1) / ccnsw(:,k-1) + soil%zse(k) /     &
               ccnsw(:,k) )
       END WHERE
    END DO

    k = 1
    WHERE( ssnow%isflag == 0 )
       coeff(:,2) = 2.0 / ( ( soil%zse(1) + xx ) / ccnsw(:,1) + soil%zse(2) /   &
            ccnsw(:,2) )
       coefa = 0.0
       coefb = REAL( coeff(:,2) )

       wblfsp = ssnow%wblf(:,k)

       xx = heat_cap_lower_limit(:,k)

       ssnow%gammzz(:,k) = MAX( heat_cap_lower_limit(:,k), &
            ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil   &
            + soil%ssat * ( wblfsp * Ccswat * Cdensity_liq +            &
            ssnow%wbfice(:,k) * Ccsice * Cdensity_liq * 0.9 ) )     &
            * soil%zse(k)

       ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + Ccgsnow * ssnow%snowd

       dtg = dels / ssnow%gammzz(:,k)

       at(:,k) = - dtg * coeff(:,k)
       ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
       bt(:,k) = 1.0 - at(:,k) - ct(:,k)

    END WHERE

    DO k = 2, ms

       WHERE( ssnow%isflag == 0 )

          wblfsp = ssnow%wblf(:,k)
          xx = soil%css * soil%rhosoil

          ssnow%gammzz(:,k) = MAX( REAL(heat_cap_lower_limit(:,k)), &
               ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil   &
               + soil%ssat * ( wblfsp * Ccswat * Cdensity_liq +            &
               ssnow%wbfice(:,k) * Ccsice * Cdensity_liq * 0.9 ) )     &
               * soil%zse(k)

          dtg = dels / ssnow%gammzz(:,k)
          at(:,k) = - dtg * coeff(:,k)
          ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
          bt(:,k) = 1.0 - at(:,k) - ct(:,k)

       END WHERE

    END DO

    WHERE( ssnow%isflag == 0 )
       bt(:,1) = bt(:,1) - canopy%dgdtg * dels / ssnow%gammzz(:,1)
       ssnow%tgg(:,1) = ssnow%tgg(:,1) + ( canopy%ga - ssnow%tgg(:,1)           &
            * REAL( canopy%dgdtg ) ) * dels / REAL( ssnow%gammzz(:,1) )
    END WHERE

    coeff(:,1-3) = 0.0  ! coeff(:,-2)

    ! 3-layer snow points done here
    WHERE( ssnow%isflag /= 0 )

       ssnow%sconds(:,1) = MAX( 0.2, MIN( 2.876e-6 * ssnow%ssdn(:,1)**2         &
            + 0.074, max_sconds ) )
       ssnow%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,2)**2 &
            & + 0.074, max_sconds) )
       ssnow%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,3)**2 &
            & + 0.074, max_sconds) )
       coeff(:,-1) = 2.0 / (ssnow%sdepth(:,1) / ssnow%sconds(:,1) &
            & + ssnow%sdepth(:,2) / ssnow%sconds(:,2) )
       coeff(:,0) = 2.0 / (ssnow%sdepth(:,2) / ssnow%sconds(:,2) &
            & + ssnow%sdepth(:,3) / ssnow%sconds(:,3) )
       coeff(:,1) = 2.0 / (ssnow%sdepth(:,3) / ssnow%sconds(:,3) &
            & + soil%zse(1) / ccnsw (:,1) )
    END WHERE

    DO k = 2, ms

       WHERE( ssnow%isflag /= 0 )                                               &
            coeff(:,k) = 2.0 / ( soil%zse(k-1) / ccnsw(:,k-1) + soil%zse(k) /     &
            ccnsw(:,k) )

    END DO

    WHERE( ssnow%isflag /= 0 )
       coefa = REAL( coeff (:,-1) )
       coefb = REAL( coeff (:,1) )
    END WHERE

    DO k = 1, 3

       WHERE( ssnow%isflag /= 0 )
          sgamm = ssnow%ssdn(:,k) * Ccgsnow * ssnow%sdepth(:,k)
          dtg = dels / sgamm
          at(:,k-3) = - dtg * coeff(:,k-3)
          ct(:,k-3) = - dtg * coeff(:,k-2)
          bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
       END WHERE

    END DO

    DO k = 1, ms

       WHERE( ssnow%isflag /= 0 )
          wblfsp = ssnow%wblf(:,k)
          xx = soil%css * soil%rhosoil

          ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css *             &
               soil%rhosoil + soil%ssat * ( wblfsp * Ccswat *     &
               Cdensity_liq + ssnow%wbfice(:,k) * Ccsice * Cdensity_liq *     &
               0.9) , &
               heat_cap_lower_limit(:,k) ) * soil%zse(k)

          dtg = dels / ssnow%gammzz(:,k)
          at(:,k) = - dtg * coeff(:,k)
          ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
          bt(:,k) = 1.0 - at(:,k) - ct(:,k)

       END WHERE

    END DO

    WHERE( ssnow%isflag /= 0 )
       sgamm = ssnow%ssdn(:,1) * Ccgsnow * ssnow%sdepth(:,1)

       bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm

       ssnow%tggsn(:,1) = ssnow%tggsn(:,1) + ( canopy%ga - ssnow%tggsn(:,1 )    &
            * REAL( canopy%dgdtg ) ) * dels / sgamm

       rhs(:,1-3) = ssnow%tggsn(:,1)
    END WHERE


    !     note in the following that tgg and tggsn are processed together
    tmp_mat(:,1:3) = REAL(ssnow%tggsn,r_2)
    tmp_mat(:,4:(ms+3)) = REAL(ssnow%tgg,r_2)

    CALL trimb( at, bt, ct, tmp_mat, ms + 3 )

    ssnow%tggsn = REAL( tmp_mat(:,1:3) )
    ssnow%tgg   = REAL( tmp_mat(:,4:(ms+3)) )
    canopy%sghflux = coefa * ( ssnow%tggsn(:,1) - ssnow%tggsn(:,2) )
    canopy%ghflux = coefb * ( ssnow%tgg(:,1) - ssnow%tgg(:,2) ) ! +ve downwards

END SUBROUTINE stempv

END MODULE stempv_mod
