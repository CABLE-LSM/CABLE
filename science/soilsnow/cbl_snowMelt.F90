MODULE snow_melting_mod

USE cbl_ssnow_data_mod

PUBLIC  snow_melting

CONTAINS

SUBROUTINE snow_melting (dels, snowmlt, ssnow, soil )

    USE cable_common_module

    REAL, INTENT(IN) :: dels   ! integration time step (s)

    REAL, DIMENSION(mp), INTENT(OUT) :: snowmlt ! snow melt

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)   :: ssnow  ! soil+snow variables

    INTEGER                 :: k,j

    REAL, DIMENSION(mp) ::                                                      &
         osm,     & !
         sgamm,   & !
         snowflx    !

    REAL, DIMENSION(mp,0:3) :: smelt1

    snowmlt= 0.0
    smelt1 = 0.0

    DO j=1,mp

       IF( ssnow%snowd(j) > 0.0 .AND. ssnow%isflag(j) == 0                      &
            .AND. ssnow%tgg(j,1) >= CTFRZ ) THEN

          ! snow covered land
          ! following done in sflux  via  ga= ... +cls*egg + ...
          ! ** land,snow,melting
          snowflx(j) = REAL((ssnow%tgg(j,1) - CTFRZ) * ssnow%gammzz(j,1))

          ! prevent snow depth going negative
          snowmlt(j) = MIN(snowflx(j) / CHLF, ssnow%snowd(j) )

          ssnow%dtmlt(j,1) = ssnow%dtmlt(j,1) + snowmlt(j) * CHLF              &
               / ssnow%gammzz(j,1)

          ssnow%snowd(j) = ssnow%snowd(j) - snowmlt(j)
          ssnow%tgg(j,1) = REAL( ssnow%tgg(j,1) - snowmlt(j) *                  &
               CHLF / ssnow%gammzz(j,1) )
       ENDIF

    END DO

    smelt1(:,0) = 0.0

    DO k = 1, 3

       !where there is snow
       WHERE( ssnow%snowd > 0.0 .AND. ssnow%isflag > 0 )

          sgamm = ssnow%ssdn(:,k) * Ccgsnow * ssnow%sdepth(:,k)

          ! snow melt refreezing
          snowflx = smelt1(:,k-1) * CHLF / dels

          ssnow%tggsn(:,k) = ssnow%tggsn(:,k) + ( snowflx * dels +              &
               smelt1(:,k-1)*Ccswat *( CTFRZ-ssnow%tggsn(:,k) ) ) &
               / ( sgamm + Ccswat*smelt1(:,k-1) )

          ! increase density due to snowmelt
          osm = ssnow%smass(:,k)
          ssnow%smass(:,k) = ssnow%smass(:,k) + smelt1(:,k-1)
          ssnow%ssdn(:,k) = MAX( 120.0, MIN( ssnow%ssdn(:,k) * osm /            &
               ssnow%smass(:,k) + Cdensity_liq * ( 1.0 - osm /           &
               ssnow%smass(:,k)), max_ssdn ) )

          ! permanent ice: fix hard-wired number in next version
          WHERE( soil%isoilm /= 9 )                                             &
               ssnow%ssdn(:,k) = MIN( 450.0, ssnow%ssdn(:,k) )

          ssnow%sdepth(:,k) = ssnow%smass(:,k) / ssnow%ssdn(:,k)

          sgamm = ssnow%smass(:,k) * Ccgsnow

          smelt1(:,k-1) = 0.0
          smelt1(:,k) = 0.0

          ! snow melting
          WHERE (ssnow%tggsn(:,k) > CTFRZ)

             snowflx = ( ssnow%tggsn(:,k) - CTFRZ ) * sgamm

             smelt1(:,k) = MIN( snowflx / CHLF, 0.6 * ssnow%smass(:,k) )

             ssnow%dtmlt(:,k) = ssnow%dtmlt(:,k) + smelt1(:,k) * CHLF / sgamm

             osm = ssnow%smass(:,k)

             ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)

             ssnow%tggsn(:,k) = ssnow%tggsn(:,k) - smelt1(:,k) * CHLF / sgamm

             ssnow%sdepth(:,k) = ssnow%smass(:,k) / ssnow%ssdn(:,k)

          END WHERE
          ! END snow melting

       END WHERE
       ! END where there is snow

    END DO

    WHERE( ssnow%snowd > 0.0 .AND. ssnow%isflag > 0 )
       snowmlt = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
       ssnow%snowd = ssnow%snowd - snowmlt
    END WHERE

  END SUBROUTINE snow_melting

END MODULE snow_melting_mod

