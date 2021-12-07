MODULE hydraulic_redistribution_mod

USE cbl_ssnow_data_mod

PUBLIC  hydraulic_redistribution

CONTAINS

  !+++++++++++++++++++  Hydraulic Redistribution Section  ++++++++++++++++++++++
  ! Science from Ryel et al. Oecologia, 2002; Lee et al., 2005, PNAS
  ! Code by LiLH 16 Feb, 2011
  ! Fixed problem of negative wb in global run by BP Mar 2011
  SUBROUTINE hydraulic_redistribution(dels, soil, ssnow, canopy, veg, met)

    USE cable_common_module, ONLY : wiltParam, satuParam

    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(soil_parameter_type), INTENT(IN) :: soil
    TYPE(canopy_type),         INTENT(IN) :: canopy
    TYPE(veg_parameter_type),  INTENT(IN) :: veg

    TYPE(soil_snow_type),   INTENT(INOUT) :: ssnow
    TYPE(met_type),         INTENT(INOUT) :: met

    REAL, PARAMETER ::                                                         &
         thetas=0.45,         & ! from Belk et al., 2007, WRR
         thetar=0.20 ,        & ! from Belk et al., 2007, WRR
         n_hr = 3.22,         & ! --
         wpsy50 = -1.0,       & ! MPa
         n_VG = 2.06,         & ! -- 2.06
         m_VG = 1.0-1.0/n_VG, & ! --
         alpha_VG = 0.00423,  & ! cm^{-1} Note: 1cmH2O=100Pa
         CRT = 125.0            ! cm MPa^-1 h^-1, default value (0.097)
    ! from Ryel et al., 2002
    REAL, DIMENSION(mp) ::                                                      &
         frootX,      & ! --
         Dtran,       & ! Swith for hr
         available,   &
         accommodate, &
         totalmoist,  &
         totalice,    &
         total2,      &
         zsetot,      &
         temp

    REAL, DIMENSION(mp,ms)::                                                    &
         S_VG, & ! --
         wpsy, & ! MPa
         C_hr    ! --

    REAL, DIMENSION(mp,ms,ms) ::                                                &
         hr_term,    & ! cm/hour
         hr_perTime    !

    INTEGER :: j, k

    zsetot = SUM(soil%zse)
    totalmoist(:) = 0.0
    totalice(:) = 0.0
    DO k=1, ms
       totalmoist(:) = totalmoist(:) + ssnow%wb(:,k)*soil%zse(k)/zsetot
       totalice(:) = totalice(:) + ssnow%wbice(:,k)*soil%zse(k)/zsetot
    ENDDO

    Dtran=0.0
    WHERE( canopy%fevc < 10.0 .AND.  totalice  < 1.e-2 )  Dtran=1.0

    DO k=1, ms
       S_VG(:,k) = MIN( 1.0, MAX( 1.0E-4, REAL(ssnow%wb(:,k)) - soil%swilt )          &
            / ( soil%ssat - soil%swilt ) )
       ! VG model, convert from cm to Pa by (*100), to MPa (*1.0E-6)
       wpsy(:,k) = -1.0 / alpha_VG * ( S_VG(:,k)**(-1.0/m_VG) - 1.0 )**(1/n_VG) &
            * 100 * 1.0E-6

       C_hr(:,k) = 1./(1+(wpsy(:,k)/wpsy50)**n_hr)
    ENDDO

    temp(:)        = 0.0
    hr_term(:,:,:) = 0.0    ! unit: cm h^{-1}
    hr_perTime(:,:,:) = 0.0

    ! setting hr_term=0 for top layer, follows Lee et al., 2005, PNAS
    DO k = ms, 3, -1

       DO j = k-1, 2, -1

          temp(:)        = 0.0
          available(:)   = 0.0
          accommodate(:) = 0.0
          frootX= MAX(0.01,MAX( veg%froot(:,k),veg%froot(:,j)))
          hr_term(:,k,j) = CRT*(wpsy(:,j)-wpsy(:,k))*MAX(C_hr(:,k),C_hr(:,j)) &
               *(veg%froot(:,k)*veg%froot(:,j))/(1-frootX) * Dtran
          hr_perTime(:,k,j) = hr_term(:,k,j)*1.0E-2/3600.0*dels ! m per timestep
          hr_perTime(:,j,k) = -1.0 * hr_perTime(:,k,j)
          hr_perTime(:,k,j) = hr_perTime(:,k,j)/soil%zse(k)
          hr_perTime(:,j,k) = hr_perTime(:,j,k)/soil%zse(j)

          ! Overwrite to give zero redistribution for all types except
          ! evergreen broadleaf (2) and c4 grass (7)
          ! NB: Hard-wired numbers should be removed in future version
          WHERE( .NOT.(veg%iveg == 2 .OR. veg%iveg == 7 ) )
             hr_perTime(:,k,j) = 0.0
             hr_perTime(:,j,k) = 0.0
          ENDWHERE

          WHERE( hr_perTime(:,k,j) < 0.0 )

             available(:)   = MAX( 0.0_r_2, ssnow%wb(:,k) -                         &
                  ( soil%swilt(:) + ( soil%sfc(:) - soil%swilt(:) )  &
                  / 3. ) )
             accommodate(:) = MAX( 0.0_r_2, soil%ssat(:) - ssnow%wb(:,j) )

             temp(:) = MAX( hr_perTime(:,k,j),                                  &
                  -1.0 * wiltParam * available(:),                     &
                  -1.0 * satuParam * accommodate(:) * soil%zse(j) /    &
                  soil%zse(k) )

             hr_perTime(:,k,j) = temp(:)
             hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)

          ELSEWHERE (hr_perTime(:,j,k) < 0.0)

             available(:)   = MAX( 0.0_r_2, ssnow%wb(:,j) -                          &
                  ( soil%swilt(:) + ( soil%sfc(:) - soil%swilt(:) )  &
                  / 3. ) )

             accommodate(:) = MAX( 0.0_r_2, soil%ssat(:) - ssnow%wb(:,k) )

             temp(:) = MAX( hr_perTime(:,j,k),                                   &
                  - 1.0 * wiltParam * available(:),                     &
                  -1.0 * satuParam * accommodate(:) * soil%zse(k) /     &
                  soil%zse(j) )

             hr_perTime(:,j,k) = temp(:)
             hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)

          ENDWHERE

          ssnow%wb(:,k) = ssnow%wb(:,k) + hr_perTime(:,k,j)
          ssnow%wb(:,j) = ssnow%wb(:,j) + hr_perTime(:,j,k)

       ENDDO

    ENDDO

    WHERE( met%tk < CTFRZ + 5.  ) Dtran=0.0

    DO k=1, ms
       S_VG(:,k) = MIN( 1.0, MAX( 1.0E-4, REAL(ssnow%wb(:,k)) - soil%swilt )          &
            / ( soil%ssat - soil%swilt ) )

       ! VG model, convert from cm to Pa by (*100), to MPa (*1.0E-6)
       wpsy(:,k) = -1.0 / alpha_VG * ( S_VG(:,k)**(-1.0/m_VG) -1.0 )**(1/n_VG)  &
            * 100 * 1.0E-6

       C_hr(:,k) = 1./(1+(wpsy(:,k)/wpsy50)**n_hr)
    ENDDO
    hr_term(:,:,:) = 0.0    ! unit: cm h^{-1}
    hr_perTime(:,:,:) = 0.0

    DO k = 1,ms-2

       DO j = k+1,ms-1

          temp(:)        = 0.0
          available(:)   = 0.0
          accommodate(:) = 0.0
          frootX= MAX(0.01,MAX( veg%froot(:,k),veg%froot(:,j)))
          hr_term(:,k,j) = CRT*(wpsy(:,j)-wpsy(:,k))*MAX(C_hr(:,k),C_hr(:,j)) &
               *(MAX(0.01,veg%froot(:,k))*MAX(0.01,veg%froot(:,j)))/(1-frootX)*Dtran
          hr_perTime(:,k,j) = hr_term(:,k,j)*1.0E-2/3600.0*dels ! m per timestep
          hr_perTime(:,j,k) = -1.0 * hr_perTime(:,k,j)
          hr_perTime(:,k,j) = hr_perTime(:,k,j)/soil%zse(k)
          hr_perTime(:,j,k) = hr_perTime(:,j,k)/soil%zse(j)

          ! Overwrite to give zero redistribution for all types except
          ! evergreen broadleaf (2) and c4 grass (7)
          ! NB: Hard-wired numbers should be removed in future version
          WHERE( .NOT.( veg%iveg == 2 .OR. veg%iveg == 7 ) )
             hr_perTime(:,k,j) = 0.0
             hr_perTime(:,j,k) = 0.0
          ENDWHERE

          WHERE( hr_perTime(:,k,j) < 0.0 )

             available(:)   = MAX( 0.0_r_2, ssnow%wb(:,k) - soil%sfc(:) )
             accommodate(:) = MAX( 0.0_r_2, soil%ssat(:) - ssnow%wb(:,j) )

             temp(:) = MAX(hr_perTime(:,k,j),                                   &
                  -1.0 * wiltParam*available(:),                       &
                  -1.0 * satuParam * accommodate(:) * soil%zse(j) /    &
                  soil%zse(k) )

             hr_perTime(:,k,j) = temp(:)
             hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)

          ELSEWHERE (hr_perTime(:,j,k) < 0.0)

             available(:)   = MAX( 0.0_r_2, ssnow%wb(:,j)- soil%sfc(:) )
             accommodate(:) = MAX( 0.0_r_2, soil%ssat(:)-ssnow%wb(:,k) )

             temp(:) = MAX(hr_perTime(:,j,k),                                   &
                  -1.0 * wiltParam*available(:),                            &
                  -1.0 * satuParam * accommodate(:) * soil%zse(k) /         &
                  soil%zse(j) )

             hr_perTime(:,j,k) = temp(:)
             hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)

          ENDWHERE

          ssnow%wb(:,k) = ssnow%wb(:,k) + hr_perTime(:,k,j)
          ssnow%wb(:,j) = ssnow%wb(:,j) + hr_perTime(:,j,k)
       ENDDO
    ENDDO

  END SUBROUTINE hydraulic_redistribution

END MODULE hydraulic_redistribution_mod

