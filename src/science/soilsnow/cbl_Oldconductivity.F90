MODULE old_soil_conductivity_mod

USE cbl_ssnow_data_mod

CONTAINS

  FUNCTION old_soil_conductivity(ssnow, soil)
    TYPE(soil_snow_type), INTENT(IN) :: ssnow
    TYPE(soil_parameter_type), INTENT(IN) :: soil

    REAL(r_2), DIMENSION(mp,ms) ::                                              &
         old_soil_conductivity  ! soil thermal conductivity (incl water/ice)

    REAL, DIMENSION(mp) ::                                                 &
         dtg,     & !
         ew       !

    INTEGER :: j,k
    REAL :: exp_arg
    LOGICAL :: direct2min = .FALSE.

    DO k = 1, ms

       DO j = 1, mp

          IF( soil%isoilm(j) == 9 ) THEN
             ! permanent ice: fix hard-wired number in next version
             old_soil_conductivity(j,k) = snow_ccnsw
          ELSE
             ew(j) = ssnow%wblf(j,k) * soil%ssat(j)
             exp_arg = ( ew(j) * LOG( 60.0 ) ) + ( ssnow%wbfice(j,k)            &
                  * soil%ssat(j) * LOG( 250.0 ) )

             IF( exp_arg > 30 ) direct2min = .TRUE.

             IF( direct2min) THEN

                old_soil_conductivity(j,k) = 1.5 * MAX( 1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 *      &
                     soil%ssat(j) /                                     &
                     MIN( ew(j), 0.5_r_2 * soil%ssat(j) ) ) ) )

             ELSE

                old_soil_conductivity(j,k) = MIN( soil%cnsd(j) * EXP( exp_arg ), 1.5_r_2 )      &
                     * MAX( 1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 *          &
                     soil%ssat(j) /                                     &
                     MIN( ew(j), 0.5_r_2 * soil%ssat(j) ) ) ) )

             ENDIF

             direct2min = .FALSE.

          ENDIF

       END DO

    END DO

  END FUNCTION old_soil_conductivity

END MODULE old_soil_conductivity_mod
