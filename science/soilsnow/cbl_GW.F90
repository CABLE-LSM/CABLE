MODULE GWstempv_mod

USE cbl_ssnow_data_mod

PUBLIC  GWstempv

CONTAINS

  SUBROUTINE GWstempv(dels, canopy, ssnow, soil)

    !*## Purpose
    ! updates soil temp and ground heat flux

    USE cable_common_module,         ONLY: cable_user
    USE total_soil_conductivity_mod, ONLY: total_soil_conductivity
    USE old_soil_conductivity_mod,   ONLY: old_soil_conductivity
    USE trimb_mod,                   ONLY : trimb

    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(canopy_type),    INTENT(INOUT) :: canopy
    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

!    REAL, DIMENSION(mp) ::    & ! rk4417 - phase2
    REAL(r_2), DIMENSION(mp) ::                                                      & 
         coefa, coefb,  & !
         sgamm            !

!    REAL, DIMENSION(mp) ::    & ! rk4417 - phase2
    REAL(r_2), DIMENSION(mp) ::                                                 &
         dtg,     & !
!         ew,      & !  not used   ! rk4417 - phase2
         xx  !,     & !
!         wblfsp     ! not used   ! rk4417 - phase2

    REAL(r_2), DIMENSION(mp,ms) ::                                              &
         ccnsw,&  ! soil thermal conductivity (incl water/ice)
         gammzz_snow

    REAL(r_2), DIMENSION(mp, -2:ms) ::                                          &
         at, bt, ct, rhs !

    REAL(r_2), DIMENSION(mp,-2:ms+1) :: coeff

    REAL(r_2), DIMENSION(mp,ms+3)    :: tmp_mat ! temp. matrix for tggsn & tgg

    INTEGER :: j,k,i
    REAL(r_2) :: exp_arg,dels_r2

    LOGICAL :: direct2min = .FALSE.
    
    dels_r2 = real(dels,r_2)      ! rk4417 - phase2
    
    at = 0._r_2             ! MMY@23Apr2023 are these taken from CABLE-GW? 
    bt = 1._r_2             ! accept them gammzz_snow is used the eq later
    ct = 0._r_2
    coeff = 0._r_2

   ssnow%otgg(:,:) = ssnow%tgg(:,:) ! MMY??? ssnow%otgg has gotten value in SUBROUTINE soil_snow_gw before call snow_processes_soil_thermal

   gammzz_snow(:,:) = 0._r_2

   k=1
   do i=1,mp
      if (ssnow%isflag(i) .ne. 0) then
         gammzz_snow(i,k) = real(Ccgsnow * ssnow%snowd(i),r_2)
      end if
   end do        

    IF (cable_user%soil_thermal_fix) THEN
       ccnsw = total_soil_conductivity(ssnow,soil)
    ELSE
       ccnsw = old_soil_conductivity(ssnow,soil)
    END IF

    xx(:) = 0.

!    WHERE(ssnow%isflag == 0)                            ! rk4417 - phase2
!       xx(:) = MAX( 0., ssnow%snowd / ssnow%ssdnn )
!       ccnsw(:,1) = ( ccnsw(:,1) - 0.2 ) * ( soil%zse(1) / ( soil%zse(1) + xx(:) ) &
!            ) + 0.2
!    END WHERE

    WHERE(ssnow%isflag == 0)
       xx(:) = MAX( 0._r_2, real(ssnow%snowd / ssnow%ssdnn,r_2) )
       ccnsw(:,1) = ( ccnsw(:,1) - 0.2_r_2 ) * ( soil%zse_vec(:,1) / ( soil%zse_vec(:,1) + xx(:) ) &
            ) + 0.2_r_2
    END WHERE
    
!    DO k = 3, ms                            ! rk4417 - phase2
!       WHERE (ssnow%isflag == 0)
!          coeff(:,k) = 2.0 / ( soil%zse(k-1) / ccnsw(:,k-1) + soil%zse(k) /     &
!               ccnsw(:,k) )
!       END WHERE
!    END DO

    DO k = 3, ms
       WHERE (ssnow%isflag == 0)
          coeff(:,k) = 2.0 / ( soil%zse_vec(:,k-1) / ccnsw(:,k-1) + soil%zse_vec(:,k) /     &
               ccnsw(:,k) )
       END WHERE
    END DO

!    k = 1                                           ! rk4417 - phase2
!    WHERE( ssnow%isflag == 0 )
!       coeff(:,2) = 2.0 / ( ( soil%zse(1) + xx(:) ) / ccnsw(:,1) + soil%zse(2) /   &
!            ccnsw(:,2) )
!       coefa = 0.0
!       coefb = REAL( coeff(:,2) )
!
!       wblfsp = ssnow%wblf(:,k)
!
!       ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)), &
!            ( 1.0 - soil%ssat_vec(:,k) ) * &
!            soil%css_vec(:,k) * soil%rhosoil_vec(:,k)   &
!            + soil%ssat_vec(:,k) * ( wblfsp * Ccs_rho_wat +            &
!            ssnow%wbfice(:,k) * Ccs_rho_ice ) )     &
!            * soil%zse_vec(:,k)
!
!       ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + Ccgsnow * ssnow%snowd
!
!       dtg = dels / ssnow%gammzz(:,k)
!
!       at(:,k) = - dtg * coeff(:,k)
!       ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
!       bt(:,k) = 1.0 - at(:,k) - ct(:,k)
!    END WHERE

    k = 1
    WHERE( ssnow%isflag == 0 )

       coeff(:,2) = 2._r_2 / ( ( soil%zse_vec(:,1) + xx(:) ) / ccnsw(:,1) + soil%zse_vec(:,2) /   &
            ccnsw(:,2) )
       coefa = 0._r_2
       coefb = coeff(:,2)
       
       ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)), &
            ( 1.0 - soil%ssat_vec(:,k) ) * &
            soil%css_vec(:,k) * soil%rhosoil_vec(:,k)   &
            + ssnow%wbliq(:,k)*real(Ccswat*Cdensity_liq,r_2)           &
            !+ ssnow%wbice(:,k)*real(C%csice*C%density_liq*0.9,r_2) )      & ! MMY
            + ssnow%wbice(:,k)*real(Ccsice*Cdensity_ice,r_2) )      & ! MMY
            * soil%zse_vec(:,k) + gammzz_snow(:,k)
       
       dtg = dels_r2 / ssnow%gammzz(:,k)
       
       at(:,k) = - dtg * coeff(:,k)
       ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
       bt(:,k) = 1.0 - at(:,k) - ct(:,k)
    END WHERE

!    DO k = 2, ms                         ! rk4417 - phase2
!
!       WHERE( ssnow%isflag == 0 )
!
!          wblfsp = ssnow%wblf(:,k)
!
!          ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)), &
!               ( 1.0 - soil%ssat_vec(:,k) ) * &
!               soil%css_vec(:,k) * soil%rhosoil_vec(:,k)   &
!               + soil%ssat_vec(:,k) * ( wblfsp * Ccs_rho_wat +            &
!               ssnow%wbfice(:,k) * Ccs_rho_ice ) )     &
!               * soil%zse_vec(:,k)
!
!          dtg = dels / ssnow%gammzz(:,k)
!          at(:,k) = - dtg * coeff(:,k)
!          ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
!          bt(:,k) = 1.0 - at(:,k) - ct(:,k)
!
!       END WHERE
!
!    END DO

    DO k = 2, ms

       WHERE( ssnow%isflag == 0 )
          
          ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)), &
               ( 1.0 - soil%ssat_vec(:,k) ) * &
               soil%css_vec(:,k) * soil%rhosoil_vec(:,k)   &
               + ssnow%wbliq(:,k)*real(Ccswat*Cdensity_liq,r_2)           &
               !+ ssnow%wbice(:,k)*real(C%csice*C%density_liq*0.9,r_2) )      & ! MMY
               + ssnow%wbice(:,k)*real(Ccsice*Cdensity_ice,r_2) )      & ! MMY
               * soil%zse_vec(:,k) + gammzz_snow(:,k)
          
          dtg = dels_r2 / ssnow%gammzz(:,k)
          
          at(:,k) = - dtg * coeff(:,k)
          ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
          bt(:,k) = 1.0 - at(:,k) - ct(:,k)
          
       END WHERE

    END DO

    
!    WHERE( ssnow%isflag == 0 )                            ! rk4417 - phase2
!       bt(:,1) = bt(:,1) - canopy%dgdtg * dels / ssnow%gammzz(:,1)
!       ssnow%tgg(:,1) = ssnow%tgg(:,1) + ( canopy%ga - ssnow%tgg(:,1)           &
!            * REAL( canopy%dgdtg ) ) * dels / REAL( ssnow%gammzz(:,1) )
!    END WHERE

    WHERE( ssnow%isflag == 0 )
       bt(:,1) = bt(:,1) - canopy%dgdtg * dels_r2 / ssnow%gammzz(:,1)
       ssnow%tgg(:,1) = ssnow%tgg(:,1) + real(( real(canopy%ga,r_2) - real(ssnow%tgg(:,1),r_2)           &
            * REAL( canopy%dgdtg ) ) * dels_r2 /  ssnow%gammzz(:,1) )
    END WHERE
    
    coeff(:,1-3) = 0.0  ! coeff(:,-2)

!    ! 3-layer snow points done here                     ! rk4417 - phase2
!    WHERE( ssnow%isflag /= 0 )
!
!       ssnow%sconds(:,1) = MAX( 0.2, MIN( 2.876e-6 * ssnow%ssdn(:,1)**2         &
!            + 0.074, max_sconds ) )
!       ssnow%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,2)**2 &
!            & + 0.074, max_sconds) )
!       ssnow%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,3)**2 &
!            & + 0.074, max_sconds) )
!       coeff(:,-1) = 2.0 / (ssnow%sdepth(:,1) / ssnow%sconds(:,1) &
!            & + ssnow%sdepth(:,2) / ssnow%sconds(:,2) )
!       coeff(:,0) = 2.0 / (ssnow%sdepth(:,2) / ssnow%sconds(:,2) &
!            & + ssnow%sdepth(:,3) / ssnow%sconds(:,3) )
!       coeff(:,1) = 2.0 / (ssnow%sdepth(:,3) / ssnow%sconds(:,3) &
!            & + soil%zse(1) / ccnsw (:,1) )
!    END WHERE

    ! 3-layer snow points done here
    WHERE( ssnow%isflag /= 0 )

       ssnow%sconds(:,1) = MAX( 0.2, MIN( 2.876e-6 * ssnow%ssdn(:,1)**2         &
            + 0.074, max_sconds ) )
       ssnow%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,2)**2 &
            & + 0.074, max_sconds) )
       ssnow%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,3)**2 &
            & + 0.074, max_sconds) )
       coeff(:,-1) = 2._r_2 / (real(ssnow%sdepth(:,1) / ssnow%sconds(:,1),r_2) &
            & + real(ssnow%sdepth(:,2) / ssnow%sconds(:,2),r_2) )
       coeff(:,0) = 2._r_2 / (real(ssnow%sdepth(:,2) / ssnow%sconds(:,2),r_2) &
            & + real(ssnow%sdepth(:,3) / ssnow%sconds(:,3),r_2) )
       coeff(:,1) = 2._r_2 / (real(ssnow%sdepth(:,3) / ssnow%sconds(:,3),r_2) &
            & + soil%zse_vec(:,1) / ccnsw (:,1) )
    END WHERE

!    DO k = 2, ms                            ! rk4417 - phase2
!       WHERE( ssnow%isflag /= 0 )                                               &
!            coeff(:,k) = 2.0 / ( soil%zse(k-1) / ccnsw(:,k-1) + soil%zse(k) /     &
!            ccnsw(:,k) )
!    END DO

    DO k = 2, ms
       WHERE( ssnow%isflag /= 0 )                                               &
            coeff(:,k) = 2._r_2 / ( soil%zse_vec(:,k-1) / ccnsw(:,k-1) + soil%zse_vec(:,k) /     &
            ccnsw(:,k) )
    END DO
    
!    WHERE( ssnow%isflag /= 0 )                 ! rk4417 - phase2
!       coefa = REAL( coeff (:,-1) )
!       coefb = REAL( coeff (:,1) )
!    END WHERE

    WHERE( ssnow%isflag /= 0 )
       coefa = coeff (:,-1)         
       coefb = coeff (:,1)
    END WHERE
    
!    DO k = 1, 3                                     ! rk4417 - phase2
!       WHERE( ssnow%isflag /= 0 )
!          sgamm = ssnow%ssdn(:,k) * Ccgsnow * ssnow%sdepth(:,k)
!          dtg = dels / sgamm
!          at(:,k-3) = - dtg * coeff(:,k-3)
!          ct(:,k-3) = - dtg * coeff(:,k-2)
!          bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
!       END WHERE
!    END DO

    DO k = 1, 3
       WHERE( ssnow%isflag /= 0 )
          sgamm = real(ssnow%ssdn(:,k) * Ccgsnow * ssnow%sdepth(:,k),r_2)
          dtg = dels_r2 / sgamm
          at(:,k-3) = - dtg * coeff(:,k-3)
          ct(:,k-3) = - dtg * coeff(:,k-2)
          bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
       END WHERE
    END DO
    
!    DO k = 1, ms                                 ! rk4417 - phase2
!       WHERE( ssnow%isflag /= 0 )
!          wblfsp = ssnow%wblf(:,k)
!
!          ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)),&
!               ( 1.0 - soil%ssat_vec(:,k) ) * soil%css_vec(:,k) *             &
!               soil%rhosoil_vec(:,k) + soil%ssat_vec(:,k) * ( wblfsp * Ccs_rho_wat +&
!               ssnow%wbfice(:,k) * Ccs_rho_ice)) * &
!               soil%zse_vec(:,k)
!
!          dtg = dels / ssnow%gammzz(:,k)
!          at(:,k) = - dtg * coeff(:,k)
!          ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
!          bt(:,k) = 1.0 - at(:,k) - ct(:,k)
!       END WHERE
!    END DO

    DO k = 1, ms
       WHERE( ssnow%isflag /= 0 )

          ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)), &
               ( 1.0 - soil%ssat_vec(:,k) ) * &
               soil%css_vec(:,k) * soil%rhosoil_vec(:,k)   &
               + ssnow%wbliq(:,k)*real(Ccswat*Cdensity_liq,r_2)           &
               !+ ssnow%wbice(:,k)*real(C%csice*C%density_liq*0.9,r_2) )      & ! MMY
               + ssnow%wbice(:,k)*real(Ccsice*Cdensity_ice,r_2) )      & ! MMY
               * soil%zse_vec(:,k) + gammzz_snow(:,k)
          
          dtg = dels_r2 / ssnow%gammzz(:,k)
          at(:,k) = - dtg * coeff(:,k)
          ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
          bt(:,k) = 1._r_2 - at(:,k) - ct(:,k)
       END WHERE
    END DO
    
!    WHERE( ssnow%isflag /= 0 )                                ! rk4417 - phase2
!       sgamm = ssnow%ssdn(:,1) * Ccgsnow * ssnow%sdepth(:,1)
!
!       bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm
!
!       ssnow%tggsn(:,1) = ssnow%tggsn(:,1) + ( canopy%ga - ssnow%tggsn(:,1 )    &
!            * REAL( canopy%dgdtg ) ) * dels / sgamm
!
!       rhs(:,1-3) = ssnow%tggsn(:,1)
!    END WHERE

    WHERE( ssnow%isflag /= 0 )
       sgamm = real(ssnow%ssdn(:,1) * Ccgsnow * ssnow%sdepth(:,1),r_2)
       
       bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels_r2 / sgamm
       
       ssnow%tggsn(:,1) = ssnow%tggsn(:,1) +real( ( real(canopy%ga,r_2) - real(ssnow%tggsn(:,1),r_2)    &
            * (canopy%dgdtg) * dels_r2) / sgamm )
       
       rhs(:,1-3) = ssnow%tggsn(:,1)
    END WHERE

    !     note in the following that tgg and tggsn are processed together
    tmp_mat(:,1:3) = REAL(ssnow%tggsn,r_2)
    tmp_mat(:,4:(ms+3)) = REAL(ssnow%tgg,r_2)

    CALL trimb( at, bt, ct, tmp_mat, ms + 3 )

    ssnow%tggsn = REAL( tmp_mat(:,1:3) )
    ssnow%tgg   = REAL( tmp_mat(:,4:(ms+3)) )
!    canopy%sghflux = coefa * ( ssnow%tggsn(:,1) - ssnow%tggsn(:,2) )                 ! rk4417 - phase2
!    canopy%ghflux = coefb * ( ssnow%tgg(:,1) - ssnow%tgg(:,2) ) ! +ve downwards
    canopy%sghflux = real(coefa) * ( ssnow%tggsn(:,1) - ssnow%tggsn(:,2) )
    canopy%ghflux = real(coefb) * ( ssnow%tgg(:,1) - ssnow%tgg(:,2) ) ! +ve downwards
  END SUBROUTINE GWstempv

END MODULE GWstempv_mod
