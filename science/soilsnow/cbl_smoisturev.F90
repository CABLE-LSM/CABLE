MODULE smoisturev_mod

USE cbl_ssnow_data_mod

PUBLIC  smoisturev

CONTAINS

!      Solves implicit soil moisture equation
!      Science development by Eva Kowalczyk and John McGregor, CMAR
SUBROUTINE smoisturev (dels,ssnow,soil,veg)
USE trimb_mod,                    ONLY: trimb
USE cable_common_module
IMPLICIT NONE

    REAL, INTENT(IN) :: dels    ! time step size (s)

    TYPE(soil_snow_type),      INTENT(INOUT) ::                                &
         ssnow ! soil and snow variables

    TYPE(soil_parameter_type), INTENT(INOUT) ::                                &
         soil  ! soil parameters

    TYPE(veg_parameter_type), INTENT(IN)  :: veg

    ! nmeth selects the solution method
    ! Values as follows:
    !  -1 for simple implicit D  - PREFERRED METHOD
    !   1 for fully implicit solution
    !   2 for simpler implicit
    !   3 for simple implicit D, explicit K
    !   4 for simple implicit D, implicit K
    !   0 for simple implicit D, new jlm TVD K
    INTEGER, PARAMETER ::                                                       &
         nmeth = -1

    REAL, DIMENSION(mp) ::                                                      &
         totwba,     & ! diagnostic
         totwbb,     & !
         totwbc,     & !
         totwblb,    & !
         totwblc,    & !
         wbficemx,   & !
         wblfmn,     & !
         wblfmx        !

    REAL(r_2), DIMENSION(mp) ::                                                 &
         fact,       & !
         fact2,      & !
         fluxhi,     & !
         fluxlo,     & !
         hydss,      & ! hydraulic conductivity adjusted for ice
         phi,        & !
         pwb,        & !
         rat,        & !
         speed_k,    & !
         ssatcurr_k, & !
         wbh_k,      & !
         wbl_k,      & !
         wbl_kp,     & !
         wh,         & !
         z3_k,       & !
         pwb_wbh       !


    REAL, DIMENSION(mp,ms+1) :: z1mult

    ! change dimension of at,bt,ct from 3*ms to ms (BP Jun2010)
    REAL(r_2), DIMENSION(mp,ms) ::                                              &
         at,      & ! coef "A" in finite diff eq
         bt,      & ! coef "B" in finite diff eq
         ct,      & ! coef "C" in finite diff eq
         ssatcurr   !

    REAL(r_2), DIMENSION(mp,ms+1) ::                                            &
         wbh,     & !
         z1,      & !
         z2,      & !
         z3         !

    REAL(r_2), DIMENSION(mp,0:ms) ::                                            &
         fluxh,      & !
         delt,       & !
         dtt           !

    LOGICAL :: is_open     ! Is file open?
    INTEGER ::                                                                  &
         u,    & ! I/O unit
         k

    at = 0.0
    bt = 1.0
    ct = 0.0
    z1mult(:,1) = 0.0       ! corresponds to 2b+3
    z1mult(:,ms+1) = 0.0    ! corresponds to 2b+3
    z1(:,1) = 0.0           ! i.e. K(.5),    value at surface
    z1(:,ms+1) = 0.0        ! i.e. K(ms+.5), value at bottom
    z2(:,1) = 0.0           ! i.e. K(.5),    value at surface
    z2(:,ms+1) = 0.0        ! i.e. K(ms+.5), value at bottom
    z3(:,1) = 0.0           ! i.e. K(.5),    value at surface
    z3(:,ms+1) = 0.0        ! i.e. K(ms+.5), value at bottom

    ! nmeth: equation solution technique
    IF (nmeth <= 0) THEN

       ! jlm split TVD version
       ! all land points
       delt(:,0) = 0.0
       fluxh(:,0) = 0.0
       fluxh(:,ms) = 0.0

       DO k = 1, ms-1

          ! Calculate amount of liquid soil water:
          IF (.NOT. cable_user%l_new_runoff_speed) THEN
             wbl_k = MAX( 0.01_r_2, ssnow%wb(:,k) - ssnow%wbice(:,k) )
             wbl_kp = MAX( 0.01_r_2, ssnow%wb(:,k+1) - ssnow%wbice(:,k+1) )
          ELSE
             wbl_k = MAX( 0.001_r_2, ssnow%wb(:,k) - ssnow%wbice(:,k) )
             wbl_kp = MAX( 0.001_r_2, ssnow%wb(:,k+1) - ssnow%wbice(:,k+1) )
          ENDIF

          ! Calculate difference in liq soil water b/w consecutive layers:
          delt(:,k) = wbl_kp - wbl_k

          ! especially to allow for isolated frozen layers, use min speed
          wh = MIN( wbl_k, wbl_kp )
          WHERE( ssnow%wbice(:,k) > 0.05 .OR. ssnow%wbice(:,k+1) > 0.01 )       &
               wh = 0.9*wbl_k + 0.1*wbl_kp

          ! with 50% wbice, reduce hyds by 1.e-5
          ! Calculate hyd conductivity adjusted for ice:
          hydss = soil%hyds

          speed_k = hydss * (wh / soil%ssat )**( soil%i2bp3 - 1 )

          ! update wb by TVD method
          rat = delt(:,k - 1) / ( delt(:,k) + SIGN( REAL( 1.0e-20, r_2 ),       &
               delt(:,k) ) )

          phi = MAX( 0.0_r_2, MIN( 1.0_r_2, 2.0_r_2 * rat ),                    &
               MIN( 2.0_r_2, rat ) ) ! 0 for -ve rat
          fluxhi = wh
          fluxlo = wbl_k

          ! scale speed to grid lengths per dels & limit speed for stability
          ! 1. OK too for stability
          speed_k = MIN( speed_k, REAL( 0.5 * soil%zse(k) / dels , r_2 ) )
          fluxh(:,k) = speed_k * ( fluxlo + phi * ( fluxhi - fluxlo ) )

       END DO

       ! calculate drainage (this code replaces the code in the surfb)
       k = ms

       IF (.NOT. cable_user%l_new_runoff_speed) THEN

          WHERE( ssnow%wb(:,ms) > soil%sfc(:) )

             wbl_k = MAX( 0.001_r_2, ssnow%wb(:,ms) - ssnow%wbice(:,ms) )
             wbl_kp = MAX( 0.001_r_2, soil%ssat(:) - ssnow%wbice(:,ms) )

             wh = MIN( wbl_k, wbl_kp )

             WHERE( ssnow%wbice(:,ms) .GT. 0.05 ) wh = 0.9 * wbl_k + 0.1 * wbl_kp

             ! Calculate hyd conductivity adjusted for ice:
             hydss = soil%hyds

             speed_k = hydss * ( wh / soil%ssat )**( soil%i2bp3 - 1 )
             speed_k =  0.5 * speed_k / ( 1. - MIN( 0.5_r_2, 10. * ssnow%wbice(:,ms) ) )
             fluxlo = wbl_k

             ! scale speed to grid lengths per dt & limit speed for stability
             speed_k = MIN( 0.5 * speed_k, 0.5_r_2 * soil%zse(ms) / dels )
             fluxh(:,ms) = MAX( 0.0_r_2, speed_k * fluxlo )

          END WHERE

       ELSE

          WHERE( ssnow%wb(:,ms) > soil%sfc(:) )

             wbl_k = MAX( 0.001_r_2, ssnow%wb(:,ms) - ssnow%wbice(:,ms) )
             wbl_kp = MAX( 0.001_r_2, soil%ssat(:) - ssnow%wbice(:,ms) )

             wh = MIN( wbl_k, wbl_kp )

             WHERE( ssnow%wbice(:,ms) .GT. 0.05 ) wh = 0.9 * wbl_k + 0.1 * wbl_kp

             ! Calculate hyd conductivity adjusted for ice:
             hydss = soil%hyds

             speed_k = hydss * ( wh / soil%ssat )**( soil%i2bp3 - 1 )
             speed_k =  speed_k / ( 1. - MIN( 0.5_r_2, 10. * ssnow%wbice(:,ms) ) )
             fluxlo = wbl_k

             ! scale speed to grid lengths per dt & limit speed for stability
             speed_k = MIN( speed_k, REAL(0.5 * soil%zse(ms) / dels, r_2) )
             fluxh(:,ms) = MAX( 0.0_r_2, speed_k * fluxlo )

          END WHERE

       ENDIF

       ! update wb by TVD method
       DO k = ms, 1, -1

          IF(  nmeth == -1 ) THEN ! each new wb constrained by ssat
             fluxh(:,k-1) = MIN( fluxh (:,k-1), ( soil%ssat - ssnow%wb(:,k) )   &
                  * soil%zse(k) / dels + fluxh(:,k) )
          END IF

          ! fluxh (:,ms) is drainage
          ssnow%wb(:,k) = ssnow%wb(:,k) + dels * ( fluxh(:,k-1) - fluxh(:,k) )  &
               / soil%zse(k)

          ! re-calculate wblf
          ssatcurr_k = soil%ssat - ssnow%wbice(:,k)
          dtt(:,k) = dels / ( soil%zse(k) * ssatcurr_k )

          ! this defn of wblf has different meaning from previous one in surfbv
          ! N.B. are imposing wbice<wb, so wblf <1
          ssnow%wblf(:,k) = ( ssnow%wb(:,k) - ssnow%wbice(:,k) ) / ssatcurr_k

       END DO

       ssnow%rnof2 = dels * REAL( fluxh(:,ms) ) * 1000.0

       ! wbh_k represents wblf(k-.5)
       DO k = 2, ms

          ssatcurr_k = REAL( soil%ssat, r_2 ) - ssnow%wbice(:,k)
          wbh_k = ( soil%zse(k) * ssnow%wblf(:,k-1) + soil%zse(k-1)             &
               * ssnow%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
          ! i.e. wbh**(bch+1)
          fact = wbh_k**( soil%ibp2 - 1 )

          ! with 50% wbice, reduce hbsh by 1.e-5
          pwb_wbh = ( soil%hsbh * ( 1. - MIN( 2. * MIN (0.1_r_2, MAX(           &
               ssnow%wbice(:,k-1) / MAX( 0.01_r_2, ssnow%wb(:,k-1) ),       &
               ssnow%wbice(:,k)   / MAX( 0.01_r_2, ssnow%wb(:,k) ) ) )      &
               , 0.1_r_2 ) ) )                                              &
               * MAX( soil%pwb_min, wbh_k * fact )

          ! moisture diffusivity (D) is  wbh*pwb; hsbh includes b
          ! i.e. D(k-.5)/soil%zshh(k)
          z3_k = pwb_wbh / soil%zshh (k)

          ! where dtt=dels/(soil%zse(k)*ssatcurr_k)
          at (:,k) = - dtt(:,k) * z3_k
          ct (:,k-1) = - dtt(:,k-1) * z3_k

       END DO

       bt = 1. - at - ct
       ssnow%wblf(:,1) = ssnow%wblf(:,1) + dtt(:,1) * ssnow%fwtop1 / Cdensity_liq
       ssnow%wblf(:,2) = ssnow%wblf(:,2) + dtt(:,2) * ssnow%fwtop2 / Cdensity_liq
       ssnow%wblf(:,3) = ssnow%wblf(:,3) + dtt(:,3) * ssnow%fwtop3 / Cdensity_liq

    END IF
    ! END: IF (nmeth <= 0) THEN

    IF ( nmeth > 0 ) THEN

       wbficemx = 0.0

       DO k = 1, ms

          ssatcurr(:,k) = REAL(soil%ssat,r_2) - ssnow%wbice(:,k)

          ! this defn of wblf has different meaning from previous one in surfbv
          ! N.B. are imposing wbice<wb, so wblf <1
          ssnow%wblf(:,k) = ( ssnow%wb(:,k) - ssnow%wbice(:,k) ) / ssatcurr(:,k)

          ssnow%wbfice(:,k) = REAL( ssnow%wbice(:,k) ) / soil%ssat

          wbficemx = MAX( wbficemx, REAL(ssnow%wbfice(:,k)) )
          dtt(:,k) = dels / ( soil%zse(k) * ssatcurr(:,k) )

       END DO

       IF( nmeth == 1 ) THEN ! full implicit method

          DO k = 2, ms

             wbh(:,k) = ( soil%zse(k) * ssnow%wblf(:,k-1) + soil%zse(k-1)       &
                  * ssnow%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )

             fact = wbh(:,k)**( soil%ibp2 - 1 ) ! i.e. wbh**(bch+1)
             fact2 = fact * fact
             pwb = soil%hsbh * fact

             ! moisture diffusivity (D) is  wbh*pwb
             ! other term (K) is wbh*soil%hyds*fact2
             z1(:,k) = wbh(:,k) * ( (soil%i2bp3 - 1 ) * soil%hyds * fact2       &
                  - soil%ibp2 * pwb *                                      &
                  ( ssnow%wblf(:,k) - ssnow%wblf(:,k-1 ) ) / soil%zshh (k) )

             z2(:,k) = - soil%i2bp3 * soil%hyds * fact2 + soil%ibp2 * pwb       &
                  * ( ssnow%wblf(:,k) - ssnow%wblf(:,k-1) ) / soil%zshh (k)

             z3(:,k) = pwb * wbh(:,k) / soil%zshh (k)

             at(:,k) = dtt(:,k) * (z2(:,k) * 0.5 * soil%zse(k) / soil%zshh (k)  &
                  - z3(:,k) )

          END DO

          DO k = 1, ms - 1

             ct(:,k) = dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k)             &
                  / soil%zshh (k+1) - z3(:,k+1) )

             bt(:,k) = 1.0 + dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k+1)     &
                  / soil%zshh (k+1) + z2(:,k) * 0.5 * soil%zse( MAX( k-1,  &
                  1 ) ) / soil%zshh (k) + z3(:,k+1) + z3(:,k) )

          END DO

          bt(:,ms) = 1.0 + dtt(:,ms) * ( z2(:,ms) * 0.5 * soil%zse(ms)          &
               / soil%zshh (ms) + z3(:,ms) )

          DO k = 1, ms
             ssnow%wblf(:,k) = ssnow%wblf(:,k) + dtt(:,k) *                     &
                  ( z1(:,k+1) - z1(:,k) )
          END DO

       END IF ! (nmeth == 1)

       IF (nmeth >= 2) THEN ! part implicit method

          DO k = 2, ms
             z1mult(:,k) = soil%i2bp3 ! corresponds to 2b+3
          END DO

          DO k = 2, ms ! wbh(k) represents wblf(k-.5)
             wbh(:,k) = ( soil%zse(k) * ssnow%wblf(:,k-1) + soil%zse(k-1)       &
                  * ssnow%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )

             fact = wbh(:,k)**( soil%ibp2 - 1 ) ! i.e. wbh**(bch+1)

             IF (nmeth == 2) pwb_wbh = soil%hsbh * wbh(:,k) * fact

             IF (nmeth >= 3)                                                    &
                  pwb_wbh = soil%hsbh * MAX( soil%pwb_min,wbh(:,k) * fact)

             ! moisture diffusivity (D) is  wbh*pwb
             ! other term (K) is wbh*soil%hyds*fact2
             z1(:,k) = soil%hyds * fact2 !  i.e. K(k-.5)/wbh(:,k)
             z3(:,k) = pwb_wbh / soil%zshh(k) !  i.e. D(k-.5)/soil%zshh(k)
             at(:,k) = - dtt(:,k) * z3(:,k)
             ct(:,k-1) = - dtt(:,k-1) * z3(:,k)

          END DO

          bt = 1. - at - ct

          IF (nmeth == 4) THEN ! for simple implicit D, implicit K
             bt(:,1) = bt(:,1) + dtt(:,1) * z1mult(:,1+1) &
                  * z1(:,1+1) * soil%zse(1+1) / (soil%zse(1) + soil%zse(1+1) )
             DO k = 2, ms
                at(:,k)   = at(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k)           &
                     * soil%zse(k) / (soil%zse(k) + soil%zse(k-1) )

                ct(:,k-1) = ct(:,k-1) + dtt(:,k-1) * z1mult(:,k) * z1(:,k)       &
                     * soil%zse(k-1) / ( soil%zse(k) + soil%zse(k-1) )

                bt(:,k)   = bt(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k)           &
                     * soil%zse(k-1) / ( soil%zse(k) + soil%zse(k-1) )    &
                     + dtt(:,k) * z1mult(:,k+1) * z1(:,k+1)               &
                     * soil%zse(k+1) / ( soil%zse(k) + soil%zse(k+1) )

             END DO

          END IF ! (nmeth == 4)

          DO k = 2, ms
             ! i.e. now K(k-.5)
             z1(:,k) = wbh(:,k) * z1(:,k)
          END DO

          ! the following top & bottom b.c.'s will preserve a uniform column
          !     z1(1) =z1(2)   ! simple dk/dz=0
          !     z1(ms+1)=z1(ms) ! simple dk/dz=0
          ! N.B. z1 are here +ve
          z1(:,1) = MIN( z1(:,2), z1(:,ms) )
          z1(:,ms + 1) = z1(:,1)

          ! no gravit. term if too much ice 11/12/00
          DO k = 1, ms

             IF (nmeth == 4) THEN

                WHERE( wbficemx < 0.75 )                                        &
                     ssnow%wblf(:,k) = ssnow%wblf(:,k) + dtt(:,k)                 &
                     * ( ( z1mult(:,k+1) - 1.0 ) * z1(:,k+1)    &
                     - (z1mult(:,k) - 1.0) * z1(:,k) )

             ELSE

                WHERE( wbficemx < 0.75 )                                        &
                     ssnow%wblf(:,k) = ssnow%wblf(:,k) + dtt(:,k)                 &
                     * ( z1(:,k) - z1(:,k+1) )

             END IF

          END DO

       END IF

       IF (nmeth == 3) THEN
          ! artificial fix applied here for safety (explicit nmeth only)
          DO k = 1, ms
             ssnow%wblf(:,k) = MAX( 0.0_r_2, MIN( ssnow%wblf(:,k), 1.0_r_2 ) )
          END DO
       END IF

       ssnow%wblf(:,1) = ssnow%wblf(:,1) + dtt(:,1) * ssnow%fwtop1 / Cdensity_liq
       ssnow%wblf(:,2) = ssnow%wblf(:,2) + dtt(:,2) * ssnow%fwtop2 / Cdensity_liq
       ssnow%wblf(:,3) = ssnow%wblf(:,3) + dtt(:,3) * ssnow%fwtop3 / Cdensity_liq

    END IF  ! IF (nmeth > 0)

    CALL trimb(at, bt, ct, ssnow%wblf, ms)

    DO k = 1, ms
       ssatcurr(:,k) = soil%ssat - ssnow%wbice(:,k)
       ssnow%wb(:,k) = ssnow%wblf(:,k) * ssatcurr(:,k) + ssnow%wbice(:,k)
       ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k), frozen_limit * ssnow%wb(:,k) )
    END DO

END SUBROUTINE smoisturev

END MODULE smoisturev_mod
