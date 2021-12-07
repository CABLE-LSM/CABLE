MODULE snowdensity_mod

USE cbl_ssnow_data_mod

PUBLIC  snowdensity

CONTAINS

SUBROUTINE snowdensity (dels, ssnow, soil)

    REAL, INTENT(IN) :: dels   ! integration time step (s)

    TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL, DIMENSION(mp) :: ssnow_tgg_min1
   REAL, DIMENSION(mp,3) :: ssnow_tgg_min

    ssnow_tgg_min1 = MIN( CTFRZ, ssnow%tgg(:,1) )

    WHERE( ssnow%snowd > 0.1 .AND. ssnow%isflag == 0 )

       ssnow%ssdn(:,1) = MIN( max_ssdn, MAX( 120.0, ssnow%ssdn(:,1) + dels      &
            * ssnow%ssdn(:,1) * 3.1e-6 * EXP( -0.03 * ( 273.15 -   &
            ssnow_tgg_min1 ) - MERGE( 0.046, 0.0,                  &
            ssnow%ssdn(:,1) >= 150.0 ) * ( ssnow%ssdn(:,1) - 150.0)&
            ) ) )

       ssnow%ssdn(:,1) = MIN(max_ssdn,ssnow%ssdn(:,1) + dels * 9.806      &
            & * ssnow%ssdn(:,1) * 0.75 * ssnow%snowd                             &
            & / (3.0e7 * EXP(0.021 * ssnow%ssdn(:,1) + 0.081                     &
            & * (273.15 - MIN(CTFRZ, ssnow%tgg(:,1) ) ) ) ) )

       ! permanent ice: fix hard-wired number in next version
       WHERE( soil%isoilm /= 9 ) ssnow%ssdn(:,1) = MIN( 450.0, ssnow%ssdn(:,1) )

       ssnow%sconds(:,1) = MAX( 0.2, MIN( 2.876e-6 * ssnow%ssdn(:,1)**2         &
            + 0.074, max_sconds ) )

       ssnow%sconds(:,2) = ssnow%sconds(:,1)
       ssnow%sconds(:,3) = ssnow%sconds(:,1)

       ssnow%ssdnn = ssnow%ssdn(:,1)

       ssnow%ssdn(:,2) = ssnow%ssdn(:,1)
       ssnow%ssdn(:,3) = ssnow%ssdn(:,1)

    END WHERE


    WHERE (ssnow%isflag == 1)

       ssnow%ssdn(:,1) = ssnow%ssdn(:,1) + dels * ssnow%ssdn(:,1) * 3.1e-6      &
            * EXP( -0.03 * (273.15 - MIN(CTFRZ, ssnow%tggsn(:,1)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,1) >= 150.0)                      &
            * (ssnow%ssdn(:,1) - 150.0) )

       ssnow%ssdn(:,2) = ssnow%ssdn(:,2) + dels * ssnow%ssdn(:,2) * 3.1e-6      &
            * EXP( -0.03 * (273.15 - MIN(CTFRZ, ssnow%tggsn(:,2)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,2) >= 150.0)                      &
            * (ssnow%ssdn(:,2) - 150.0) )

       ssnow%ssdn(:,3) = ssnow%ssdn(:,3) + dels * ssnow%ssdn(:,3) * 3.1e-6      &
            * EXP( -0.03 * (273.15 - MIN(CTFRZ, ssnow%tggsn(:,3)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,3) >= 150.0)                      &
            * (ssnow%ssdn(:,3) - 150.0) )

       ssnow%ssdn(:,1) = ssnow%ssdn(:,1) + dels * 9.806 * ssnow%ssdn(:,1)       &
            * ssnow%t_snwlr*ssnow%ssdn(:,1)                                    &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,1) + 0.081                      &
            * (273.15 - MIN(CTFRZ, ssnow%tggsn(:,1)))))

       ssnow%ssdn(:,2) = ssnow%ssdn(:,2) + dels * 9.806 * ssnow%ssdn(:,2)       &
            * (ssnow%t_snwlr * ssnow%ssdn(:,1) + 0.5 * ssnow%smass(:,2) )      &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,2) + 0.081                      &
            * (273.15 - MIN(CTFRZ, ssnow%tggsn(:,2)))))

       ssnow%ssdn(:,3) = ssnow%ssdn(:,3) + dels * 9.806 * ssnow%ssdn(:,3)       &
            * (ssnow%t_snwlr*ssnow%ssdn(:,1) + ssnow%smass(:,2)                &
            + 0.5*ssnow%smass(:,3))                                            &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,3) + 0.081                      &
            * (273.15 - MIN(CTFRZ, ssnow%tggsn(:,3)))))

       ssnow%sdepth(:,1) =  ssnow%smass(:,1) / ssnow%ssdn(:,1)
       ssnow%sdepth(:,2) =  ssnow%smass(:,2) / ssnow%ssdn(:,2)
       ssnow%sdepth(:,3) =  ssnow%smass(:,3) / ssnow%ssdn(:,3)

       ssnow%ssdnn = (ssnow%ssdn(:,1) * ssnow%smass(:,1) + ssnow%ssdn(:,2)      &
            * ssnow%smass(:,2) + ssnow%ssdn(:,3) * ssnow%smass(:,3) )          &
            / ssnow%snowd

       ssnow%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,1) ** 2         &
            & + 0.074, max_sconds) )
       ssnow%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,2) ** 2 &
            & + 0.074, max_sconds) )
       ssnow%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,3) ** 2 &
            & + 0.074, max_sconds) )
    END WHERE

END SUBROUTINE snowdensity

END MODULE snowdensity_mod
