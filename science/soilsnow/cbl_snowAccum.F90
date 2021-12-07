MODULE snow_accum_mod

USE cbl_ssnow_data_mod

PUBLIC  snow_accum

CONTAINS 

  SUBROUTINE snow_accum ( dels,  canopy, met, ssnow, soil )

    USE cable_common_module

    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(canopy_type), INTENT(INOUT)         :: canopy ! vegetation variables
    TYPE(met_type), INTENT(INOUT)            :: met   ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters

    REAL, DIMENSION(mp) ::                                                      &
         osm,     & !
         sgamm,   & !
         snowmlt, & !
         xxx        !

    INTEGER             :: i,j,k

    DO i=1,mp

       IF(canopy%precis(i) > 0.0 .AND. ssnow%isflag(i) == 0) THEN

          ! accumulate solid part
          ssnow%snowd(i) = MAX( ssnow%snowd(i) + met%precip_sn(i), 0.0 )

          canopy%precis(i) = canopy%precis(i) - met%precip_sn(i)

          ssnow%ssdn(i,1) = MAX( 120.0, ssnow%ssdn(i,1)                            &
               * ssnow%osnowd(i) / MAX( 0.01, ssnow%snowd(i) )              &
               + 120.0 * met%precip_sn(i) / MAX( 0.01, ssnow%snowd(i) ) )

          ssnow%ssdnn(i) = ssnow%ssdn(i,1)

          IF( canopy%precis(i) > 0.0 .AND. ssnow%tgg(i,1) < CTFRZ ) THEN

             ssnow%snowd(i) = MAX(ssnow%snowd(i) + canopy%precis(i), 0.0)

             ssnow%tgg(i,1) = ssnow%tgg(i,1) + canopy%precis(i) * CHLF               &
                  / ( REAL( ssnow%gammzz(i,1) ) + Ccswat *canopy%precis(i) )
             ! change density due to water being added
             ssnow%ssdn(i,1) = MIN( max_ssdn, MAX( 120.0, ssnow%ssdn(i,1)          &
                  * ssnow%osnowd(i) / MAX( 0.01, ssnow%snowd(i) ) + Cdensity_liq  &
                  * canopy%precis(i) / MAX( 0.01, ssnow%snowd(i) )  ) )

             ! permanent ice: fix hard-wired number in next version
             IF( soil%isoilm(i) /= 9 )                                             &
                  ssnow%ssdn(i,1) = MIN( 450.0, ssnow%ssdn(i,1) )

             canopy%precis(i) = 0.0
             ssnow%ssdnn(i) = ssnow%ssdn(i,1)

          END IF

       END IF ! (canopy%precis > 0. .and. ssnow%isflag == 0)

       IF(canopy%precis(i) > 0.0 .AND.  ssnow%isflag(i) > 0) THEN

          ! add solid precip
          ssnow%snowd(i) = MAX( ssnow%snowd(i) + met%precip_sn(i), 0.0 )

          canopy%precis(i) = canopy%precis(i) - met%precip_sn(i)  ! remaining liquid precip

          ! update top snow layer with fresh snow
          osm(i) = ssnow%smass(i,1)
          ssnow%smass(i,1) = ssnow%smass(i,1) + met%precip_sn(i)
          ssnow%ssdn(i,1) = MAX( 120.0,ssnow%ssdn(i,1) * osm(i) / ssnow%smass(i,1)    &
               + 120.0 * met%precip_sn(i) / ssnow%smass(i,1) )

          ssnow%sdepth(i,1) = MAX( 0.02, ssnow%smass(i,1) / ssnow%ssdn(i,1) )

          ! add liquid precip
          IF( canopy%precis(i) > 0.0 ) THEN

             ssnow%snowd(i) = MAX( ssnow%snowd(i) + canopy%precis(i), 0.0 )
             sgamm(i) = ssnow%ssdn(i,1) * Ccgsnow * ssnow%sdepth(i,1)
             osm(i) = ssnow%smass(i,1)

             ssnow%tggsn(i,1) = ssnow%tggsn(i,1) + canopy%precis(i) * CHLF           &
                  * osm(i) / (sgamm(i) * ssnow%osnowd(i) )
             ssnow%smass(i,1) = ssnow%smass(i,1) + canopy%precis(i)                   &
                  * osm(i)/ssnow%osnowd(i)

             ssnow%ssdn(i,1) = MAX( 120.0, MIN( ssnow%ssdn(i,1) * osm(i) /            &
                  ssnow%smass(i,1) +  Cdensity_liq *                        &
                  ( 1.0 - osm(i) / ssnow%smass(i,1) ), max_ssdn ) )

             ! permanent ice: fix hard-wired number in next version
             IF( soil%isoilm(i) /= 9 ) &
                  ssnow%ssdn(i,1) = MIN( 450.0, ssnow%ssdn(i,1) )

             ssnow%sdepth(i,1) = ssnow%smass(i,1)/ssnow%ssdn(i,1)

             !layer 2
             sgamm(i) = ssnow%ssdn(i,2) * Ccgsnow * ssnow%sdepth(i,2)
             osm(i) = ssnow%smass(i,2)
             ssnow%tggsn(i,2) = ssnow%tggsn(i,2) + canopy%precis(i) * CHLF           &
                  * osm(i) / ( sgamm(i) * ssnow%osnowd(i) )
             ssnow%smass(i,2) = ssnow%smass(i,2) + canopy%precis(i)                   &
                  * osm(i) / ssnow%osnowd(i)
             ssnow%ssdn(i,2) = MAX( 120.0, MIN( ssnow%ssdn(i,2) * osm(i) /            &
                  ssnow%smass(i,2) + Cdensity_liq *                         &
                  ( 1.0 - osm(i) / ssnow%smass(i,2) ), max_ssdn ) )

             ! permanent ice: fix hard-wired number in next version
             IF( soil%isoilm(i) /= 9 ) &
                  ssnow%ssdn(i,2) = MIN( 450.0, ssnow%ssdn(i,2) )

             ssnow%sdepth(i,2) = ssnow%smass(i,2) / ssnow%ssdn(i,2)

             !layer 3
             sgamm(i) = ssnow%ssdn(i,3) * Ccgsnow * ssnow%sdepth(i,3)
             osm(i) = ssnow%smass(i,3)
             ssnow%tggsn(i,3) = ssnow%tggsn(i,3) + canopy%precis(i) * CHLF           &
                  * osm(i) / ( sgamm(i) * ssnow%osnowd(i) )
             ssnow%smass(i,3) = ssnow%smass(i,3) + canopy%precis(i)                   &
                  * osm(i) / ssnow%osnowd(i)
             ssnow%ssdn(i,3) = MAX( 120.0, MIN( ssnow%ssdn(i,3) * osm(i) /             &
                  ssnow%smass(i,3) + Cdensity_liq *                          &
                  ( 1.0 - osm(i) / ssnow%smass(i,3) ), max_ssdn ) )

             ! permanent ice: fix hard-wired number in next version
             IF( soil%isoilm(i) /= 9 ) &
                  ssnow%ssdn(i,3) = MIN(450.0,ssnow%ssdn(i,3))

             ssnow%sdepth(i,3) = ssnow%smass(i,3) / ssnow%ssdn(i,3)

             canopy%precis(i) = 0.0

          ENDIF

       ENDIF

    ENDDO

    ! 'fess' is for soil evap and 'fes' is for soil evap plus soil puddle evap
    canopy%segg = canopy%fess / CHL
    canopy%segg = ( canopy%fess + canopy%fes_cor ) / CHL

    ! Initialise snow evaporation:
    ssnow%evapsn = 0
    DO i=1,mp
       ! Snow evaporation and dew on snow
       ! NB the conditions on when %fes applies to %segg or %evapsn MUST(!!)
       ! match those used to set %cls in the latent_heat_flux calculations
       ! for moisture conservation purposes
       ! Ticket 137 - using %cls as the trigger not %snowd
       IF( ssnow%cls(i) == 1.1335 ) THEN
          !WHERE( ssnow%snowd > 0.1 )

          ssnow%evapsn(i) = dels * ( canopy%fess(i) + canopy%fes_cor(i) ) / ( CHL + CHLF )
          xxx(i) = ssnow%evapsn(i)

          IF( ssnow%isflag(i) == 0 .AND. canopy%fess(i) + canopy%fes_cor(i).GT. 0.0 )    &
               ssnow%evapsn(i) = MIN( ssnow%snowd(i), xxx(i) )

          IF( ssnow%isflag(i)  > 0 .AND. canopy%fess(i) + canopy%fes_cor(i) .GT. 0.0 )   &
               ssnow%evapsn(i) = MIN( 0.9 * ssnow%smass(i,1), xxx(i) )

          ssnow%snowd(i) = ssnow%snowd(i) - ssnow%evapsn(i)

          IF( ssnow%isflag(i) > 0 ) THEN
             ssnow%smass(i,1) = ssnow%smass(i,1)  - ssnow%evapsn(i)
             ssnow%sdepth(i,1) = MAX( 0.02, ssnow%smass(i,1) / ssnow%ssdn(i,1) )
          ENDIF

          canopy%segg(i) = 0.0

          !INH: cls package
          !we still need to conserve moisture/energy when evapsn is limited
          !this is a key point of moisture non-conservation

       ENDIF

    ENDDO
  END SUBROUTINE snow_accum

END MODULE snow_accum_mod
