MODULE snowl_adjust_mod

USE cbl_ssnow_data_mod

PUBLIC  snowl_adjust

CONTAINS

SUBROUTINE snowl_adjust(dels, ssnow, canopy )

IMPLICIT NONE
    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE(canopy_type), INTENT(INOUT)    :: canopy

    INTEGER :: k

    REAL(r_2), DIMENSION(mp) ::                                                 &
         excd,    & !
         excm,    & !
         frac,    & !
         xfrac     !

    REAL, DIMENSION(mp) :: osm

    INTEGER :: api ! active patch counter


    ! adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    WHERE( ssnow%isflag > 0 )

       WHERE( ssnow%sdepth(:,1) > ssnow%t_snwlr )

          excd = ssnow%sdepth(:,1) - ssnow%t_snwlr
          excm = excd * ssnow%ssdn(:,1)
          ssnow%sdepth(:,1) = ssnow%sdepth(:,1) - REAL(excd)
          osm = ssnow%smass(:,1)
          ssnow%smass(:,1) = ssnow%smass(:,1) - REAL(excm)

          osm = ssnow%smass(:,2)
          ssnow%smass(:,2) = MAX( 0.01, ssnow%smass(:,2) + REAL(excm) )

          ssnow%ssdn(:,2) = REAL( MAX( 120.0_r_2, MIN( REAL( max_ssdn, r_2 ),   &
               ssnow%ssdn(:,2) * osm / ssnow%smass(:,2) +          &
               ssnow%ssdn(:,1) * excm / ssnow%smass(:,2) ) ) )

          ssnow%sdepth(:,2) =  ssnow%smass(:,2) / ssnow%ssdn(:,2)

          ssnow%tggsn(:,2) = REAL( ssnow%tggsn(:,2) * osm / ssnow%smass(:,2)   &
               + ssnow%tggsn(:,1) * excm / ssnow%smass(:,2) )

          ! following line changed to fix -ve sdepth (EK 21Dec2007)
          ssnow%smass(:,3) = MAX( 0.01, ssnow%snowd - ssnow%smass(:,1)         &
               - ssnow%smass(:,2) )

       ELSEWHERE ! ssnow%sdepth(:,1) < ssnow%t_snwlr

          ! 1st layer
          excd = ssnow%t_snwlr - ssnow%sdepth(:,1)
          excm = excd * ssnow%ssdn(:,2)
          osm = ssnow%smass(:,1)
          ssnow%smass(:,1) = ssnow%smass(:,1) + REAL(excm)
          ssnow%sdepth(:,1) = ssnow%t_snwlr
          ssnow%ssdn(:,1) = REAL( MAX( 120.0_r_2, MIN( REAL( max_ssdn,r_2 ),    &
               ssnow%ssdn(:,1) * osm / ssnow%smass(:,1)            &
               + ssnow%ssdn(:,2) * excm / ssnow%smass(:,1) ) ) )

          ssnow%tggsn(:,1) = REAL( ssnow%tggsn(:,1) * osm / ssnow%smass(:,1)   &
               + ssnow%tggsn(:,2) * excm / ssnow%smass(:,1) )

          ! 2nd layer
          ssnow%smass(:,2) = MAX( 0.01, ssnow%smass(:,2) - REAL(excm) )
          ssnow%sdepth(:,2) = ssnow%smass(:,2) / ssnow%ssdn(:,2)

          ! following line changed to fix -ve sdepth (EK 21Dec2007)
          ssnow%smass(:,3) = MAX( 0.01, ssnow%snowd - ssnow%smass(:,1)          &
               - ssnow%smass(:,2) )

       END WHERE

    END WHERE

    DO  api=1,mp

       IF( ssnow%isflag(api).GT.0 ) THEN

          frac(api) = ssnow%smass(api,2) / MAX( 0.02, ssnow%smass(api,3) )
          ! if frac > 0.6 or frac < 0.74 do nothing
          ! HOW TO translate this to xfrac
          xfrac(api) = 2.0/3.0/ frac(api)

          IF( xfrac(api) > 1.0 ) THEN

             excm(api) = (xfrac(api) - 1.0) * ssnow%smass(api,2)
             osm(api) = ssnow%smass(api,2)

             ! changed 0.02 to 0.01 to fix -ve sdepth (EK 21Dec2007)
             ssnow%smass(api,2) = MAX( 0.01, ssnow%smass(api,2) +               &
                  REAL( excm(api) ) )

             ssnow%tggsn(api,2) = ssnow%tggsn(api,2) * osm(api) /               &
                  ssnow%smass(api,2) +  ssnow%tggsn(api,3)      &
                  * REAL( excm(api) )/ ssnow%smass(api,2)

             ssnow%ssdn(api,2) = MAX( 120.0, MIN( max_ssdn, ssnow%ssdn(api,2) * &
                  osm(api) / ssnow%smass(api,2) +                &
                  ssnow%ssdn(api,3) * REAL( excm(api) )          &
                  / ssnow%smass(api,2) ) )

             ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
             ssnow%smass(api,3) = MAX( 0.01, ssnow%snowd(api) -                 &
                  ssnow%smass(api,1) - ssnow%smass(api,2) )

             ssnow%sdepth(api,3) = MAX( 0.02, ssnow%smass(api,3) /              &
                  ssnow%ssdn(api,3) )

          ELSE! xfrac < 1

             excm(api) = ( 1 - xfrac(api) ) * ssnow%smass(api,2)
             ssnow%smass(api,2) = MAX(0.01, ssnow%smass(api,2) - REAL(excm(api)))
             ssnow%sdepth(api,2) = MAX(0.02, ssnow%smass(api,2) /                &
                  ssnow%ssdn(api,2) )

             osm(api) = ssnow%smass(api,3)
             ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
             ssnow%smass(api,3) = MAX(0.01, &
                  ssnow%snowd(api) - ssnow%smass(api,1) -        &
                  ssnow%smass(api,2) )


             ssnow%tggsn(api,3) = ssnow%tggsn(api,3) * osm(api) /                &
                  ssnow%smass(api,3) +  ssnow%tggsn(api,2) *     &
                  REAL( excm(api) ) / ssnow%smass(api,3)
             ssnow%ssdn(api,3) = MAX(120.0, MIN( max_ssdn, ssnow%ssdn(api, 3 )*  &
                  osm(api) / ssnow%smass(api,3) +                 &
                  ssnow%ssdn(api,2) * REAL( excm(api) )           &
                  / ssnow%smass(api,3) ) )
             ssnow%sdepth(api,3) = ssnow%smass(api,3) /  ssnow%ssdn(api,3)

          END IF

          ssnow%isflag(api) = 1

          ssnow%ssdnn(api) = ( ssnow%ssdn(api,1) * ssnow%sdepth(api,1) +         &
               ssnow%ssdn(api,2) * ssnow%sdepth(api,2) +           &
               ssnow%ssdn(api,3) * ssnow%sdepth(api,3) )           &
               / ( ssnow%sdepth(api,1) + ssnow%sdepth(api,2)       &
               + ssnow%sdepth(api,3) )

       END IF

    END DO

END SUBROUTINE snowl_adjust

END MODULE snowl_adjust_mod

