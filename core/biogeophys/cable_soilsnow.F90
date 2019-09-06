!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.
! ==============================================================================
!
! Purpose: All routines for calculating soil temperature and moisture
!          and snow calculations
!
! Contact: Eva.Kowalczyk@csiro.au
!
! History: v2.0 Tighter water budget
!          v2.0 Hydraulic redistribution subroutine (with namelist switch).
!               NB Currently hard-wired to veg types 2 and 7
!                  (usually evergreen broadleaf and c4 grass)
!          v2.0 ssoil variable renamed ssnow
!          Aug 2017 - applied changes for cls and rev_corr packages
!
! ==============================================================================

MODULE cable_soil_snow_module

  USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
       veg_parameter_type, canopy_type, met_type,        &
       balances_type, r_2, ms, mp, bgc_pool_type

  USE cable_data_module, ONLY : issnow_type, point2constants

  USE cable_common_module, ONLY: cable_user,snow_ccnsw,snmin,&
       max_ssdn,max_sconds,frozen_limit,&
       max_glacier_snowd

  IMPLICIT NONE

  PRIVATE

  TYPE ( issnow_type ), SAVE :: C

  ! This module contains the following subroutines:
  PUBLIC soil_snow ! must be available outside this module
  PUBLIC snow_processes_soil_thermal,trimb

CONTAINS


  SUBROUTINE trimb (a, b, c, rhs, kmax)

    INTEGER, INTENT(IN)                  :: kmax ! no. of discrete layers

    REAL(r_2), DIMENSION(:,:), INTENT(IN) ::                                    &
         a,    & ! coef "A" in finite diff eq
         b,    & ! coef "B" in finite diff eq
         c       ! coef "C" in finite diff eq

    REAL(r_2), DIMENSION(:,:), INTENT(INOUT)  :: rhs ! right hand side of eq

    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) ::                                &
         e, temp, g

    INTEGER :: k   ! do lloop counter


    e(:,1) = c(:,1) / b(:,1)
    DO k = 2, kmax - 1
       temp(:,k) = 1. / ( b(:,k) - a(:,k) * e(:,k-1) )
       e(:,k) = c(:,k) * temp(:,k)
    END DO

    g(:,1) = rhs(:,1) / b(:,1)
    DO k = 2, kmax - 1
       g(:,k) = ( rhs(:,k) - a(:,k) * g(:,k-1) ) * temp(:,k)
    END DO

    ! do back substitution to give answer now
    rhs(:,kmax) = ( rhs(:,kmax) - a(:,kmax) * g(:,kmax-1) )                     &
         / ( b(:,kmax) - a(:,kmax) * e(:,kmax-1) )

    DO k = kmax - 1, 1, - 1
       rhs(:,k) = g(:,k) - e(:,k) * rhs(:,k + 1)
    END DO

  END SUBROUTINE trimb

  ! -----------------------------------------------------------------------------

  !      Solves implicit soil moisture equation
  !      Science development by Eva Kowalczyk and John McGregor, CMAR
  SUBROUTINE smoisturev (dels,ssnow,soil,veg)

    USE cable_common_module

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
       ssnow%wblf(:,1) = ssnow%wblf(:,1) + dtt(:,1) * ssnow%fwtop1 / C%density_liq
       ssnow%wblf(:,2) = ssnow%wblf(:,2) + dtt(:,2) * ssnow%fwtop2 / C%density_liq
       ssnow%wblf(:,3) = ssnow%wblf(:,3) + dtt(:,3) * ssnow%fwtop3 / C%density_liq

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

       ssnow%wblf(:,1) = ssnow%wblf(:,1) + dtt(:,1) * ssnow%fwtop1 / C%density_liq
       ssnow%wblf(:,2) = ssnow%wblf(:,2) + dtt(:,2) * ssnow%fwtop2 / C%density_liq
       ssnow%wblf(:,3) = ssnow%wblf(:,3) + dtt(:,3) * ssnow%fwtop3 / C%density_liq

    END IF  ! IF (nmeth > 0)

    CALL trimb(at, bt, ct, ssnow%wblf, ms)

    DO k = 1, ms
       ssatcurr(:,k) = soil%ssat - ssnow%wbice(:,k)
       ssnow%wb(:,k) = ssnow%wblf(:,k) * ssatcurr(:,k) + ssnow%wbice(:,k)
       ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k), frozen_limit * ssnow%wb(:,k) )
    END DO

  END SUBROUTINE smoisturev

  ! -----------------------------------------------------------------------------

  SUBROUTINE snowdensity (dels, ssnow, soil)

    REAL, INTENT(IN) :: dels   ! integration time step (s)

    TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    INTEGER, DIMENSION(mp,3) :: ssnow_isflag_ssdn
    REAL, DIMENSION(mp) :: ssnow_tgg_min1
    REAL, DIMENSION(mp,3) :: dels_ssdn, ssnow_tgg_min

    ssnow_isflag_ssdn = SPREAD( ssnow%isflag,2,mp)

    dels_ssdn = SPREAD( SPREAD( dels, 1, mp ), 2,  mp )
    ssnow_tgg_min1 = MIN( C%TFRZ, ssnow%tgg(:,1) )

    WHERE( ssnow%snowd > 0.1 .AND. ssnow%isflag == 0 )

       ssnow%ssdn(:,1) = MIN( max_ssdn, MAX( 120.0, ssnow%ssdn(:,1) + dels      &
            * ssnow%ssdn(:,1) * 3.1e-6 * EXP( -0.03 * ( 273.15 -   &
            ssnow_tgg_min1 ) - MERGE( 0.046, 0.0,                  &
            ssnow%ssdn(:,1) >= 150.0 ) * ( ssnow%ssdn(:,1) - 150.0)&
            ) ) )

       ssnow%ssdn(:,1) = MIN(max_ssdn,ssnow%ssdn(:,1) + dels * 9.806      &
            & * ssnow%ssdn(:,1) * 0.75 * ssnow%snowd                             &
            & / (3.0e7 * EXP(0.021 * ssnow%ssdn(:,1) + 0.081                     &
            & * (273.15 - MIN(C%TFRZ, ssnow%tgg(:,1) ) ) ) ) )

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
            * EXP( -0.03 * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,1)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,1) >= 150.0)                      &
            * (ssnow%ssdn(:,1) - 150.0) )

       ssnow%ssdn(:,2) = ssnow%ssdn(:,2) + dels * ssnow%ssdn(:,2) * 3.1e-6      &
            * EXP( -0.03 * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,2)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,2) >= 150.0)                      &
            * (ssnow%ssdn(:,2) - 150.0) )

       ssnow%ssdn(:,3) = ssnow%ssdn(:,3) + dels * ssnow%ssdn(:,3) * 3.1e-6      &
            * EXP( -0.03 * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,3)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,3) >= 150.0)                      &
            * (ssnow%ssdn(:,3) - 150.0) )

       ssnow%ssdn(:,1) = ssnow%ssdn(:,1) + dels * 9.806 * ssnow%ssdn(:,1)       &
            * ssnow%t_snwlr*ssnow%ssdn(:,1)                                    &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,1) + 0.081                      &
            * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,1)))))

       ssnow%ssdn(:,2) = ssnow%ssdn(:,2) + dels * 9.806 * ssnow%ssdn(:,2)       &
            * (ssnow%t_snwlr * ssnow%ssdn(:,1) + 0.5 * ssnow%smass(:,2) )      &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,2) + 0.081                      &
            * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,2)))))

       ssnow%ssdn(:,3) = ssnow%ssdn(:,3) + dels * 9.806 * ssnow%ssdn(:,3)       &
            * (ssnow%t_snwlr*ssnow%ssdn(:,1) + ssnow%smass(:,2)                &
            + 0.5*ssnow%smass(:,3))                                            &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,3) + 0.081                      &
            * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,3)))))

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

  ! -----------------------------------------------------------------------------

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
            .AND. ssnow%tgg(j,1) >= C%TFRZ ) THEN

          ! snow covered land
          ! following done in sflux  via  ga= ... +cls*egg + ...
          ! ** land,snow,melting
          snowflx(j) = REAL((ssnow%tgg(j,1) - C%TFRZ) * ssnow%gammzz(j,1))

          ! prevent snow depth going negative
          snowmlt(j) = MIN(snowflx(j) / C%HLF, ssnow%snowd(j) )

          ssnow%dtmlt(j,1) = ssnow%dtmlt(j,1) + snowmlt(j) * C%HLF              &
               / ssnow%gammzz(j,1)

          ssnow%snowd(j) = ssnow%snowd(j) - snowmlt(j)
          ssnow%tgg(j,1) = REAL( ssnow%tgg(j,1) - snowmlt(j) *                  &
               C%HLF / ssnow%gammzz(j,1) )
       ENDIF

    END DO

    smelt1(:,0) = 0.0

    DO k = 1, 3

       !where there is snow
       WHERE( ssnow%snowd > 0.0 .AND. ssnow%isflag > 0 )

          sgamm = ssnow%ssdn(:,k) * C%cgsnow * ssnow%sdepth(:,k)

          ! snow melt refreezing
          snowflx = smelt1(:,k-1) * C%HLF / dels

          ssnow%tggsn(:,k) = ssnow%tggsn(:,k) + ( snowflx * dels +              &
               smelt1(:,k-1)*C%cswat *( C%TFRZ-ssnow%tggsn(:,k) ) ) &
               / ( sgamm + C%cswat*smelt1(:,k-1) )

          ! increase density due to snowmelt
          osm = ssnow%smass(:,k)
          ssnow%smass(:,k) = ssnow%smass(:,k) + smelt1(:,k-1)
          ssnow%ssdn(:,k) = MAX( 120.0, MIN( ssnow%ssdn(:,k) * osm /            &
               ssnow%smass(:,k) + C%density_liq * ( 1.0 - osm /           &
               ssnow%smass(:,k)), max_ssdn ) )

          ! permanent ice: fix hard-wired number in next version
          WHERE( soil%isoilm /= 9 )                                             &
               ssnow%ssdn(:,k) = MIN( 450.0, ssnow%ssdn(:,k) )

          ssnow%sdepth(:,k) = ssnow%smass(:,k) / ssnow%ssdn(:,k)

          sgamm = ssnow%smass(:,k) * C%cgsnow

          smelt1(:,k-1) = 0.0
          smelt1(:,k) = 0.0

          ! snow melting
          WHERE (ssnow%tggsn(:,k) > C%TFRZ)

             snowflx = ( ssnow%tggsn(:,k) - C%TFRZ ) * sgamm

             smelt1(:,k) = MIN( snowflx / C%HLF, 0.6 * ssnow%smass(:,k) )

             ssnow%dtmlt(:,k) = ssnow%dtmlt(:,k) + smelt1(:,k) * C%HLF / sgamm

             osm = ssnow%smass(:,k)

             ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)

             ssnow%tggsn(:,k) = ssnow%tggsn(:,k) - smelt1(:,k) * C%HLF / sgamm

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

  ! -----------------------------------------------------------------------------

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

          IF( canopy%precis(i) > 0.0 .AND. ssnow%tgg(i,1) < C%TFRZ ) THEN

             ssnow%snowd(i) = MAX(ssnow%snowd(i) + canopy%precis(i), 0.0)

             ssnow%tgg(i,1) = ssnow%tgg(i,1) + canopy%precis(i) * C%HLF               &
                  / ( REAL( ssnow%gammzz(i,1) ) + C%cswat *canopy%precis(i) )
             ! change density due to water being added
             ssnow%ssdn(i,1) = MIN( max_ssdn, MAX( 120.0, ssnow%ssdn(i,1)          &
                  * ssnow%osnowd(i) / MAX( 0.01, ssnow%snowd(i) ) + C%density_liq  &
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
             sgamm(i) = ssnow%ssdn(i,1) * C%cgsnow * ssnow%sdepth(i,1)
             osm(i) = ssnow%smass(i,1)

             ssnow%tggsn(i,1) = ssnow%tggsn(i,1) + canopy%precis(i) * C%HLF           &
                  * osm(i) / (sgamm(i) * ssnow%osnowd(i) )
             ssnow%smass(i,1) = ssnow%smass(i,1) + canopy%precis(i)                   &
                  * osm(i)/ssnow%osnowd(i)

             ssnow%ssdn(i,1) = MAX( 120.0, MIN( ssnow%ssdn(i,1) * osm(i) /            &
                  ssnow%smass(i,1) +  C%density_liq *                        &
                  ( 1.0 - osm(i) / ssnow%smass(i,1) ), max_ssdn ) )

             ! permanent ice: fix hard-wired number in next version
             IF( soil%isoilm(i) /= 9 ) &
                  ssnow%ssdn(i,1) = MIN( 450.0, ssnow%ssdn(i,1) )

             ssnow%sdepth(i,1) = ssnow%smass(i,1)/ssnow%ssdn(i,1)

             !layer 2
             sgamm(i) = ssnow%ssdn(i,2) * C%cgsnow * ssnow%sdepth(i,2)
             osm(i) = ssnow%smass(i,2)
             ssnow%tggsn(i,2) = ssnow%tggsn(i,2) + canopy%precis(i) * C%HLF           &
                  * osm(i) / ( sgamm(i) * ssnow%osnowd(i) )
             ssnow%smass(i,2) = ssnow%smass(i,2) + canopy%precis(i)                   &
                  * osm(i) / ssnow%osnowd(i)
             ssnow%ssdn(i,2) = MAX( 120.0, MIN( ssnow%ssdn(i,2) * osm(i) /            &
                  ssnow%smass(i,2) + C%density_liq *                         &
                  ( 1.0 - osm(i) / ssnow%smass(i,2) ), max_ssdn ) )

             ! permanent ice: fix hard-wired number in next version
             IF( soil%isoilm(i) /= 9 ) &
                  ssnow%ssdn(i,2) = MIN( 450.0, ssnow%ssdn(i,2) )

             ssnow%sdepth(i,2) = ssnow%smass(i,2) / ssnow%ssdn(i,2)

             !layer 3
             sgamm(i) = ssnow%ssdn(i,3) * C%cgsnow * ssnow%sdepth(i,3)
             osm(i) = ssnow%smass(i,3)
             ssnow%tggsn(i,3) = ssnow%tggsn(i,3) + canopy%precis(i) * C%HLF           &
                  * osm(i) / ( sgamm(i) * ssnow%osnowd(i) )
             ssnow%smass(i,3) = ssnow%smass(i,3) + canopy%precis(i)                   &
                  * osm(i) / ssnow%osnowd(i)
             ssnow%ssdn(i,3) = MAX( 120.0, MIN( ssnow%ssdn(i,3) * osm(i) /             &
                  ssnow%smass(i,3) + C%density_liq *                          &
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
    canopy%segg = canopy%fess / C%HL
    canopy%segg = ( canopy%fess + canopy%fes_cor ) / C%HL

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

          ssnow%evapsn(i) = dels * ( canopy%fess(i) + canopy%fes_cor(i) ) / ( C%HL + C%HLF )
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

  ! -----------------------------------------------------------------------------

  SUBROUTINE surfbv (dels, met, ssnow, soil, veg, canopy )

    USE cable_common_module

    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(canopy_type), INTENT(IN)       :: canopy

    TYPE(met_type),       INTENT(INOUT) :: met    ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow  ! soil+snow variables

    TYPE(veg_parameter_type),  INTENT(IN)     :: veg
    TYPE(soil_parameter_type), INTENT(INOUT)  :: soil  ! soil parameters

    !jhan:cable.nml
    INTEGER            :: nglacier  ! 0 original, 1 off, 2 new Eva

    REAL, DIMENSION(mp) ::                                                      &
         rnof5,      & !
         sfact,      & !
         sgamm,      & !
         smasstot,   & !
         talb,       & ! snow albedo
         tmp           ! temporary value
    REAL(r_2), DIMENSION(mp) :: xxx

    REAL, DIMENSION(mp,0:3) :: smelt1

    REAL :: wb_lake_T, rnof2_T, ratio
    INTEGER :: k,j

    IF( cable_runtime%UM ) THEN
       nglacier = 0
    ELSE
       nglacier = 2
    ENDIF

    CALL smoisturev( dels, ssnow, soil, veg )

    DO k = 1, ms
       xxx = REAL( soil%ssat,r_2 )
       ssnow%rnof1 = ssnow%rnof1 + REAL( MAX( ssnow%wb(:,k) - xxx, 0.0_r_2 )  &
            * 1000.0 )  * soil%zse(k)
       ssnow%wb(:,k) = MAX( 0.01_r_2, MIN( ssnow%wb(:,k), xxx ) )
    END DO
    ! for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
    ! account for soil/ice cracking
    ! fracm = MIN(0.2, 1. - MIN(1., ssnow%wb(:,ms) / soil%sfc ) )
    ! ssnow%wb(:,ms) = ssnow%wb(:,ms) &
    !                  + fracm*ssnow%rnof1/(1000.0*soil%zse(ms))
    ! ssnow%rnof1 = (1. - fracm) * ssnow%rnof1

    ! Scaling  runoff to kg/m^2/s to match rest of the model
    !jhan:replace nested wheres

    !---  glacier formation
    rnof5= 0.

    IF (nglacier == 2) THEN

       smelt1=0.
       WHERE( ssnow%snowd > max_glacier_snowd )

          rnof5 = MIN( 0.1, ssnow%snowd - max_glacier_snowd )

          !---- change local tg to account for energy - clearly not best method
          WHERE( ssnow%isflag == 0 )
             smasstot = 0.0
             ssnow%tgg(:,1) = ssnow%tgg(:,1) - rnof5 * C%HLF                    &
                  / REAL( ssnow%gammzz(:,1) )
             ssnow%snowd = ssnow%snowd - rnof5
          ELSEWHERE
             smasstot = ssnow%smass(:,1) + ssnow%smass(:,2) + ssnow%smass(:,3)
          END WHERE

       END WHERE

       DO k = 1, 3

          WHERE( ssnow%snowd > max_glacier_snowd  .AND.  ssnow%isflag > 0 )
             sgamm = ssnow%ssdn(:,k) * C%cgsnow * ssnow%sdepth(:,k)
             smelt1(:,k) = MIN( rnof5 * ssnow%smass(:,k) / smasstot,            &
                  0.2 * ssnow%smass(:,k) )
             ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)
             ssnow%snowd = ssnow%snowd - smelt1(:,k)
          END WHERE

       END DO

       WHERE( ssnow%isflag > 0 ) rnof5 = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)

    END IF

    !  Rescale drainage to remove water added to lakes (wb_lake)
    ssnow%sinfil = 0.0
    WHERE( veg%iveg == 16 )
       ssnow%sinfil  = MIN( ssnow%rnof1, ssnow%wb_lake ) ! water that can be extracted from the rnof1
       ssnow%rnof1   = MAX( 0.0, ssnow%rnof1 - ssnow%sinfil )
       ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
       ssnow%sinfil  = MIN( ssnow%rnof2, ssnow%wb_lake ) ! water that can be extracted from the rnof2
       ssnow%rnof2   = MAX( 0.0, ssnow%rnof2 - ssnow%sinfil )
       ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
       xxx = MAX(0.0_r_2, (ssnow%wb(:,ms) - REAL(soil%sfc(:),r_2))*soil%zse(ms)*1000.0)
       ssnow%sinfil  = MIN( REAL(xxx), ssnow%wb_lake )
       ssnow%wb(:,ms) = ssnow%wb(:,ms) - REAL(ssnow%sinfil / (soil%zse(ms)*1000.0), r_2)
       ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
       xxx = MAX(0.0_r_2, (ssnow%wb(:,ms) - 0.5*(soil%sfc + soil%swilt))*soil%zse(ms)*1000.0)
       ssnow%sinfil  = MIN( REAL(xxx), ssnow%wb_lake )
       ssnow%wb(:,ms) = ssnow%wb(:,ms) - ssnow%sinfil / (soil%zse(ms)*1000.0)
       ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
    ENDWHERE

    !wb_lake_T = sum( ssnow%wb_lake )
    !rnof2_T = sum( ssnow%rnof2 )
    !ratio = min( 1., wb_lake_T/max(rnof2_T,1.))
    !ssnow%rnof2 = ssnow%rnof2 - ratio*ssnow%rnof2
    !ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ratio*ssnow%rnof2)

    !  Rescale drainage to remove water added to lakes (wb_lake)
    !wb_lake_T = 0.0
    !rnof2_T = 0.
    !DO j=1,mp
    !   IF( ssnow%wb_lake(j) >  0.0 ) wb_lake_T = wb_lake_T + ssnow%wb_lake(j)
    !   rnof2_T = rnof2_T + ssnow%rnof2(j)
    !END DO
    !ratio = min( 1., wb_lake_T/max(rnof2_T,1.))
    !ssnow%rnof2 = ssnow%rnof2 - ratio*ssnow%rnof2
    !ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ratio*ssnow%rnof2)

    ssnow%rnof1 = ssnow%rnof1 / dels + rnof5/dels
    ssnow%rnof2 = ssnow%rnof2 / dels
    ssnow%runoff = ssnow%rnof1 + ssnow%rnof2

  END SUBROUTINE surfbv

  ! -----------------------------------------------------------------------------

  ! calculates temperatures of the soil
  ! tgg - new soil/snow temperature
  ! ga - heat flux from the atmosphere (ground heat flux)
  ! ccnsw - soil thermal conductivity, including water/ice
  SUBROUTINE stempv(dels, canopy, ssnow, soil)
    USE cable_common_module, ONLY: cable_user
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

    at = 0.0
    bt = 1.0
    ct = 0.0
    coeff = 0.0

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

       xx = soil%heat_cap_lower_limit(:,1)

       ssnow%gammzz(:,k) = MAX( (soil%heat_cap_lower_limit(:,1)), &
            ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil   &
            + soil%ssat * ( wblfsp * C%cswat * C%density_liq +            &
            ssnow%wbfice(:,k) * C%csice * C%density_liq * 0.9 ) )     &
            * soil%zse(k)

       ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + C%cgsnow * ssnow%snowd

       dtg = dels / ssnow%gammzz(:,k)

       at(:,k) = - dtg * coeff(:,k)
       ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
       bt(:,k) = 1.0 - at(:,k) - ct(:,k)

    END WHERE

    DO k = 2, ms

       WHERE( ssnow%isflag == 0 )

          wblfsp = ssnow%wblf(:,k)
          xx = soil%css * soil%rhosoil

          ssnow%gammzz(:,k) = MAX( REAL(soil%heat_cap_lower_limit(:,1)), &
               ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil   &
               + soil%ssat * ( wblfsp * C%cswat * C%density_liq +            &
               ssnow%wbfice(:,k) * C%csice * C%density_liq * 0.9 ) )     &
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
          sgamm = ssnow%ssdn(:,k) * C%cgsnow * ssnow%sdepth(:,k)
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
               soil%rhosoil + soil%ssat * ( wblfsp * C%cswat *     &
               C%density_liq + ssnow%wbfice(:,k) * C%csice * C%density_liq *     &
               0.9) , &
               (soil%heat_cap_lower_limit(:,k)) ) * soil%zse(k)

          dtg = dels / ssnow%gammzz(:,k)
          at(:,k) = - dtg * coeff(:,k)
          ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
          bt(:,k) = 1.0 - at(:,k) - ct(:,k)

       END WHERE

    END DO

    WHERE( ssnow%isflag /= 0 )
       sgamm = ssnow%ssdn(:,1) * C%cgsnow * ssnow%sdepth(:,1)

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

  ! -----------------------------------------------------------------------------

  SUBROUTINE snowcheck(dels, ssnow, soil, met )

    USE cable_common_module

    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE(met_type),       INTENT(INOUT) :: met ! all met forcing

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters

    INTEGER :: k,j


    DO j=1,mp

       IF( ssnow%snowd(j) <= 0.0 ) THEN

          ssnow%isflag(j) = 0
          ssnow%ssdn(j,:) = 120.0
          ssnow%ssdnn(j) = 120.0
          ssnow%tggsn(j,:) = C%TFRZ
          ssnow%sdepth(j,1) = ssnow%snowd(j) / ssnow%ssdn(j,1)

          ssnow%sdepth(j,2) = 0.
          ssnow%sdepth(j,3) = 0.

          ssnow%smass(j,1) = ssnow%snowd(j)
          ssnow%smass(j,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
          ssnow%smass(j,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007

          ! in loop: IF( ssnow%snowd(j) <= 0.0 ) THEN
       ELSEIF( ssnow%snowd(j) < snmin * ssnow%ssdnn(j) ) THEN

          IF( ssnow%isflag(j) == 1 ) THEN
             ssnow%ssdn(j,1) = ssnow%ssdnn(j)
             ssnow%tgg(j,1) = ssnow%tggsn(j,1)
          ENDIF

          ssnow%isflag(j) = 0
          ssnow%ssdnn(j) = MIN( 400.0, MAX( 120.0, ssnow%ssdn(j,1) ) )

          ssnow%tggsn(j,:) = MIN( C%TFRZ,ssnow%tgg(j,1) )

          ssnow%sdepth(j,1) = ssnow%snowd(j) / ssnow%ssdn(j,1)
          ssnow%sdepth(j,2) = 0.0
          ssnow%sdepth(j,3) = 0.0

          ssnow%smass(j,1) = ssnow%snowd(j)
          ssnow%smass(j,2) = 0.0
          ssnow%smass(j,3) = 0.0

          ssnow%ssdn(j,:) = ssnow%ssdnn(j)

          IF( .NOT.cable_user%CABLE_RUNTIME_COUPLED ) THEN
             IF( soil%isoilm(j) == 9 .AND. ktau_gl <= 2 )                       &
                                ! permanent ice: fixed hard-wired number in next version
                  ssnow%ssdnn(j) = 700.0
          ENDIF


       ELSE ! in loop: IF( ssnow%snowd(j) <= 0.0 ) THEN
          ! sufficient snow now for 3 layer snowpack

          IF( ssnow%isflag(j) == 0 ) THEN

             ssnow%tggsn(j,:) = MIN( C%TFRZ, ssnow%tgg(j,1) )

             ssnow%ssdn(j,2) = ssnow%ssdn(j,1)
             ssnow%ssdn(j,3) = ssnow%ssdn(j,1)

             IF( .NOT. cable_user%cable_runtime_coupled) THEN
                IF( soil%isoilm(j) == 9 .AND. ktau_gl <= 2 ) THEN
                   ! permanent ice: fix hard-wired number in next version
                   ssnow%ssdn(j,1)  = 450.0
                   ssnow%ssdn(j,2)  = 580.0
                   ssnow%ssdn(j,3)  = 600.0
                ENDIF
             ENDIF

             ssnow%sdepth(j,1) = ssnow%t_snwlr(j)

             ssnow%smass(j,1)  =  ssnow%t_snwlr(j) * ssnow%ssdn(j,1)

             ssnow%smass(j,2)  = ( ssnow%snowd(j) - ssnow%smass(j,1) ) * 0.4
             ssnow%smass(j,3)  = ( ssnow%snowd(j) - ssnow%smass(j,1) ) * 0.6

             ssnow%sdepth(j,2) = ssnow%smass(j,2) / ssnow%ssdn(j,2)
             ssnow%sdepth(j,3) = ssnow%smass(j,3) / ssnow%ssdn(j,3)

             ssnow%ssdnn(j) = ( ssnow%ssdn(j,1) * ssnow%smass(j,1) +            &
                  ssnow%ssdn(j,2) * ssnow%smass(j,2) +             &
                  ssnow%ssdn(j,3) * ssnow%smass(j,3) )             &
                  / ssnow%snowd(j)

          ENDIF

          ssnow%isflag(j) = 1

       ENDIF ! END: IF( ssnow%snowd(j) <= 0.0 ) THEN


    ENDDO ! END: DO j=1,mp

  END SUBROUTINE snowcheck

  ! -----------------------------------------------------------------------------

  SUBROUTINE snowl_adjust(dels, ssnow, canopy )

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

  ! -----------------------------------------------------------------------------

  SUBROUTINE soilfreeze(dels, soil, ssnow)
    USE cable_common_module
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    REAL(r_2), DIMENSION(mp)           :: sicefreeze
    REAL(r_2), DIMENSION(mp)           :: sicemelt
    REAL, DIMENSION(mp)           :: xx
    INTEGER k

    xx = 0.
    DO k = 1, ms

       WHERE (ssnow%tgg(:,k) < C%TFRZ &
            & .AND. frozen_limit * ssnow%wb(:,k) - ssnow%wbice(:,k) > .001)

          sicefreeze = MIN( MAX( 0.0_r_2, ( frozen_limit * ssnow%wb(:,k) -      &
               ssnow%wbice(:,k) ) ) * soil%zse(k) * 1000.0,             &
               ( C%TFRZ - ssnow%tgg(:,k) ) * ssnow%gammzz(:,k) / C%HLF )
          ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k) + sicefreeze / (soil%zse(k)  &
               * 1000.0), frozen_limit * ssnow%wb(:,k) )
          xx = soil%css * soil%rhosoil
          ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)),           &
               REAL((1.0 - soil%ssat) * soil%css * soil%rhosoil ,r_2)            &
               + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(C%cswat * C%density_liq,r_2)   &
               + ssnow%wbice(:,k) * REAL(C%csice * C%density_liq * 0.9,r_2))* &
               REAL( soil%zse(k),r_2 )

          WHERE (k == 1 .AND. ssnow%isflag == 0)
             ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + C%cgsnow * ssnow%snowd
          END WHERE
          ssnow%tgg(:,k) = ssnow%tgg(:,k) + REAL(sicefreeze)                    &
               * C%HLF / REAL(ssnow%gammzz(:,k) )

       ELSEWHERE( ssnow%tgg(:,k) > C%TFRZ .AND. ssnow%wbice(:,k) > 0. )

          sicemelt = MIN( ssnow%wbice(:,k) * soil%zse(k) * 1000.0,              &
               ( ssnow%tgg(:,k) - C%TFRZ ) * ssnow%gammzz(:,k) / C%HLF )

          ssnow%wbice(:,k) = MAX( 0.0_r_2, ssnow%wbice(:,k) - sicemelt          &
               / (soil%zse(k) * 1000.0) )
          xx = soil%css * soil%rhosoil
          ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)),       &
               REAL((1.0-soil%ssat) * soil%css * soil%rhosoil,r_2)             &
               + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(C%cswat*C%density_liq,r_2)   &
               + ssnow%wbice(:,k) * REAL(C%csice * C%density_liq * 0.9,r_2))            &
               * REAL(soil%zse(k),r_2)
          WHERE (k == 1 .AND. ssnow%isflag == 0)
             ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + C%cgsnow * ssnow%snowd
          END WHERE
          ssnow%tgg(:,k) = ssnow%tgg(:,k) - REAL(sicemelt)                     &
               * C%HLF / REAL(ssnow%gammzz(:,k))

       END WHERE

    END DO

  END SUBROUTINE soilfreeze

  ! -----------------------------------------------------------------------------

  SUBROUTINE remove_trans(dels, soil, ssnow, canopy, veg)

    USE cable_common_module, ONLY : redistrb, cable_user

    ! Removes transpiration water from soil.
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    REAL(r_2), DIMENSION(mp,0:ms) :: diff
    REAL(r_2), DIMENSION(mp)      :: xx,xxd,evap_cur
    REAL(r_2), DIMENSION(mp)      :: needed, available, difference
    INTEGER k, i

    IF (cable_user%FWSOIL_SWITCH == 'hydraulics') THEN
       ! This follows the default extraction logic, but instead of weighting
       ! by froot, we are weighting by the frac uptake we calculated when we
       ! were weighting the soil water potential.
       !
       ! Martin De Kauwe, 22/02/19

       needed = 0._r_2
       difference = 0._r_2
       available = 0._r_2

       DO k = 1, ms
          WHERE (canopy%fevc > 0.0)

             ! Calculate the amount of water we wish to extract from each
             ! layer, kg/m2
             needed = canopy%fevc * dels / C%HL * &
                         ssnow%fraction_uptake(:,k)

             ! Calculate the amount of water available in the layer
             available = MAX(0.0, ssnow%wb(:,k) - soil%swilt) * &
                            (soil%zse(k) * C%density_liq)

             difference = available - needed

             ! Calculate new layer water balance
             WHERE (difference < 0.0)
                ! We don't have sufficent water to supply demand, extract only
                ! the remaining SW in the layer
                ssnow%wb(:,k) = ssnow%wb(:,k) - available / &
                                   (soil%zse(k) * C%density_liq)
             ELSEWHERE
                ! We have sufficent water to supply demand, extract needed SW
                ! from the layer
                ssnow%wb(:,k) = ssnow%wb(:,k) - needed / &
                                   (soil%zse(k) * C%density_liq)
             ENDWHERE

          END WHERE   !fvec > 0
       END DO   !ms

    ELSE IF (cable_user%FWSOIL_switch.NE.'Haverd2013') THEN
       xx = 0.; xxd = 0.; diff(:,:) = 0.
       DO k = 1,ms

          ! Removing transpiration from soil:
          WHERE (canopy%fevc > 0.0 )     ! convert to mm/dels

             ! Calculate the amount (perhaps moisture/ice limited)
             ! which can be removed:
             xx = canopy%fevc * dels / C%HL * veg%froot(:,k) + diff(:,k-1)   ! kg/m2
             diff(:,k) = MAX( 0.0_r_2, ssnow%wb(:,k) - soil%swilt) &      ! m3/m3
                  * soil%zse(k)*1000.0
             xxd = xx - diff(:,k)

             WHERE ( xxd .GT. 0.0 )
                ssnow%wb(:,k) = ssnow%wb(:,k) - diff(:,k) / (soil%zse(k)*1000.0)
                diff(:,k) = xxd
             ELSEWHERE
                ssnow%wb(:,k) = ssnow%wb(:,k) - xx / (soil%zse(k)*1000.0)
                diff(:,k) = 0.0
             ENDWHERE

          END WHERE

       END DO

    ELSE
       WHERE (canopy%fevc .LT. 0.0_r_2)
          canopy%fevw = canopy%fevw+canopy%fevc
          canopy%fevc = 0.0_r_2
       END WHERE
       DO k = 1,ms
          ssnow%wb(:,k) = ssnow%wb(:,k) - ssnow%evapfbl(:,k)/(soil%zse(k)*1000.0)

          !  write(59,*) k,  ssnow%wb(:,k),  ssnow%evapfbl(:,k)/(soil%zse(k)*1000.0)
          !  write(59,*)
       ENDDO


    ENDIF

  END SUBROUTINE remove_trans

  ! -----------------------------------------------------------------------------

  ! Inputs:
  !        dt_in - time step in sec
  !        ktau_in - time step no.
  !        ga      - ground heat flux W/m^2
  !        dgdtg   -
  !        condxpr - total precip reaching the ground (liquid and solid)
  !        scondxpr - precip (solid only)
  !        fev   - transpiration (W/m2)
  !        fes   - soil evaporation (W/m2)
  !        isoil - soil type
  !        ivegt - vegetation type
  ! Output
  !        ssnow
  SUBROUTINE soil_snow(dels, soil, ssnow, canopy, met, bal, veg, bgc)
    USE cable_common_module
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    TYPE (balances_type), INTENT(INOUT)      :: bal
    TYPE (bgc_pool_type), INTENT(IN)         :: bgc

    INTEGER             :: k, i
    REAL, DIMENSION(mp) :: snowmlt
    REAL, DIMENSION(mp) :: totwet
    REAL, DIMENSION(mp) :: weting
    REAL, DIMENSION(mp) :: xx, tgg_old, tggsn_old
    REAL(r_2), DIMENSION(mp) :: xxx,deltat,sinfil1,sinfil2,sinfil3
    REAL                :: zsetot
    INTEGER, SAVE :: ktau =0

    REAL, DIMENSION(ms) :: root_length

    CALL point2constants( C )

    ktau = ktau +1

    !jhan - make switchable
    ! appropriate for ACCESS1.0
    !max_glacier_snowd = 50000.0
    ! appropriate for ACCESS1.3
    !max_glacier_snowd = 1100.0

    zsetot = SUM(soil%zse)
    ssnow%tggav = 0.
    DO k = 1, ms
       ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
    END DO


    IF( cable_runtime%offline .OR. cable_runtime%mk3l ) THEN
       ssnow%t_snwlr = 0.05
    ENDIF

    ssnow%fwtop1 = 0.0
    ssnow%fwtop2 = 0.0
    ssnow%fwtop3 = 0.0
    ssnow%runoff = 0.0 ! initialise total runoff
    ssnow%rnof1 = 0.0 ! initialise surface runoff
    ssnow%rnof2 = 0.0 ! initialise deep drainage
    ssnow%smelt = 0.0 ! initialise snowmelt
    ssnow%dtmlt = 0.0
    ssnow%osnowd = ssnow%snowd

    IF (cable_user%soil_thermal_fix) THEN
       soil%heat_cap_lower_limit(:,:) = 0.01  !never allow /0
    ELSE
       DO k=1,ms
          soil%heat_cap_lower_limit(:,k) = soil%css(:) * soil%rhosoil(:)
       END DO
    END IF

    IF( .NOT.cable_user%cable_runtime_coupled ) THEN

       IF( ktau_gl <= 1 ) THEN

          IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
          ! after discussion with BP
          ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
          ssnow%wbtot = 0.0
          DO k = 1, ms
             ssnow%wb(:,k)  = MIN( soil%ssat, MAX( REAL(ssnow%wb(:,k)), soil%swilt ) )
          END DO

          ssnow%wb(:,ms-2)  = MIN( soil%ssat, MAX( REAL(ssnow%wb(:,ms-2)),           &
               0.5 * ( soil%sfc + soil%swilt ) ) )
          ssnow%wb(:,ms-1)  = MIN( soil%ssat, MAX( REAL(ssnow%wb(:,ms-1)),           &
               0.8 * soil%sfc ) )
          ssnow%wb(:,ms)    = MIN( soil%ssat, MAX( REAL(ssnow%wb(:,ms)), soil%sfc ) )

          DO k = 1, ms

             WHERE( ssnow%tgg(:,k) <= C%TFRZ .AND. ssnow%wbice(:,k) <= 0.01 )   &
                  ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)

             WHERE( ssnow%tgg(:,k) < C%TFRZ)                                    &
                  ssnow%wbice(:,k) = frozen_limit * ssnow%wb(:,k)

          END DO

          WHERE (soil%isoilm == 9)
             ! permanent ice: fix hard-wired number in next version
             ssnow%snowd = max_glacier_snowd
             ssnow%osnowd = max_glacier_snowd
             ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0
             ssnow%wb(:,1) = 0.95 * soil%ssat
             ssnow%wb(:,2) = 0.95 * soil%ssat
             ssnow%wb(:,3) = 0.95 * soil%ssat
             ssnow%wb(:,4) = 0.95 * soil%ssat
             ssnow%wb(:,5) = 0.95 * soil%ssat
             ssnow%wb(:,6) = 0.95 * soil%ssat
             ssnow%wbice(:,1) = 0.90 * ssnow%wb(:,1)
             ssnow%wbice(:,2) = 0.90 * ssnow%wb(:,2)
             ssnow%wbice(:,3) = 0.90 * ssnow%wb(:,3)
             ssnow%wbice(:,4) = 0.90 * ssnow%wb(:,4)
             ssnow%wbice(:,5) = 0.90 * ssnow%wb(:,5)
             ssnow%wbice(:,6) = 0.90 * ssnow%wb(:,6)
          ENDWHERE

          xx=REAL(soil%heat_cap_lower_limit(:,1))

          ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
               & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * C%cswat * C%density_liq &
               & + ssnow%wbice(:,1) * C%csice * C%density_liq * .9, xx ) * soil%zse(1)
       END IF
    ENDIF  ! if(.NOT.cable_runtime_coupled)


    IF (ktau <= 1)       THEN

       xx=soil%heat_cap_lower_limit(:,1)
       ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil      &
            & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * C%cswat * C%density_liq           &
            & + ssnow%wbice(:,1) * C%csice * C%density_liq * .9, xx ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * C%cgsnow * ssnow%snowd

    END IF

    ssnow%wbliq = ssnow%wb - ssnow%wbice


    DO k = 1, ms ! for stempv

       ! Set liquid soil water fraction (fraction of saturation value):
       ssnow%wblf(:,k) = MAX( 0.01_r_2, (ssnow%wb(:,k) - ssnow%wbice(:,k)) )    &
            & / REAL(soil%ssat,r_2)

       ! Set ice soil water fraction (fraction of saturation value):
       ssnow%wbfice(:,k) = REAL(ssnow%wbice(:,k)) / soil%ssat
    END DO

    CALL snowcheck (dels, ssnow, soil, met )

    CALL snowdensity (dels, ssnow, soil)

    CALL snow_accum (dels, canopy, met, ssnow, soil )

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add snow melt to global snow melt variable:
    ssnow%smelt = snowmlt

    ! Adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    CALL snowl_adjust(dels, ssnow, canopy )

    CALL stempv(dels, canopy, ssnow, soil)

    ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add new snow melt to global snow melt variable:
    ssnow%smelt = ssnow%smelt + snowmlt

    ! PH: mgk576, 13/10/17, added two funcs
    IF (cable_user%FWSOIL_SWITCH == 'hydraulics') THEN
       DO i = 1, mp

          CALL calc_soil_root_resistance(ssnow, soil, veg, bgc, root_length, i)
          CALL calc_swp(ssnow, soil, i)
          CALL calc_weighted_swp_and_frac_uptake(ssnow, soil, canopy, &
                                                 root_length, i)

       END DO
    END IF

    CALL remove_trans(dels, soil, ssnow, canopy, veg)

    CALL  soilfreeze(dels, soil, ssnow)

    totwet = canopy%precis + ssnow%smelt

    ! total available liquid including puddle
    weting = totwet + MAX(0._r_2,ssnow%pudsto - canopy%fesp/C%HL*dels)
    xxx=soil%ssat - ssnow%wb(:,1)

    sinfil1 = MIN( 0.95*xxx*soil%zse(1)*C%density_liq, weting) !soil capacity
    xxx=soil%ssat - ssnow%wb(:,2)
    sinfil2 = MIN( 0.95*xxx*soil%zse(2)*C%density_liq, weting - REAL(sinfil1)) !soil capacity
    xxx=soil%ssat - ssnow%wb(:,3)
    sinfil3 = MIN( 0.95*xxx*soil%zse(3)*C%density_liq,weting-REAL(sinfil1)-REAL(sinfil2))

    ! net water flux to the soil
    ssnow%fwtop1 = sinfil1 / dels - canopy%segg
    ssnow%fwtop2 = sinfil2 / dels
    ssnow%fwtop3 = sinfil3 / dels

    ! Puddle for the next time step
    ssnow%pudsto = MAX( 0._r_2, weting - sinfil1 - sinfil2 - sinfil3 )
    ssnow%rnof1 = MAX(0.,ssnow%pudsto - ssnow%pudsmx)
    ssnow%pudsto = ssnow%pudsto - ssnow%rnof1

    CALL surfbv(dels, met, ssnow, soil, veg, canopy )


    ! correction required for energy balance in online simulations
    IF( cable_runtime%um ) THEN

       !cls package - rewritten for flexibility
       canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
       !canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%dfe_ddq * ssnow%ddq_dtg)
       canopy%fes_cor = ssnow%dtmlt(:,1)*ssnow%dfe_dtg

       canopy%fhs = canopy%fhs+canopy%fhs_cor
       canopy%fes = canopy%fes+canopy%fes_cor

       !REV_CORR associated changes to other energy balance terms
       !NB canopy%fns changed not rad%flws as the correction term needs to
       !pass through the canopy in entirety, not be partially absorbed
       IF (cable_user%L_REV_CORR) THEN
          canopy%fns_cor = ssnow%dtmlt(:,1)*ssnow%dfn_dtg
          canopy%ga_cor = ssnow%dtmlt(:,1)*canopy%dgdtg

          canopy%fns = canopy%fns + canopy%fns_cor
          canopy%ga = canopy%ga + canopy%ga_cor

          canopy%fess = canopy%fess + canopy%fes_cor
       ENDIF
    ENDIF

    ! redistrb (set in cable.nml) by default==.FALSE.
    IF( redistrb )                                                              &
         CALL hydraulic_redistribution( dels, soil, ssnow, canopy, veg, met )

    ssnow%smelt = ssnow%smelt/dels

    ! Set weighted soil/snow surface temperature
    ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

    ssnow%wbliq = ssnow%wb - ssnow%wbice

    ssnow%wbtot = 0.0
    DO k = 1, ms
       ssnow%wbtot = ssnow%wbtot + REAL(ssnow%wb(:,k)*1000.0*soil%zse(k),r_2)
    END DO

  END SUBROUTINE soil_snow

  ! -----------------------------------------------------------------------------

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
          !WHERE( .NOT.(veg%iveg == 2 .OR. veg%iveg == 7 ) )
          WHERE( .NOT.(veg%iveg == 2 .OR. veg%iveg == 7 .OR. veg%iveg >= 18) )
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

    WHERE( met%tk < C%TFRZ + 5.  ) Dtran=0.0

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
          !WHERE( .NOT.( veg%iveg == 2 .OR. veg%iveg == 7 ) )
          WHERE( .NOT.( veg%iveg == 2 .OR. veg%iveg == 7 .OR. veg%iveg >= 18) )
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
                  (ssnow%tgg(j,k) .LT. C%TFRZ) .OR. &
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

  SUBROUTINE snow_processes_soil_thermal(dels,ssnow,soil,veg,canopy,met,bal,snowmlt)
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    TYPE (balances_type), INTENT(INOUT)      :: bal
    REAL, DIMENSION(:),  INTENT(INOUT)       :: snowmlt

    INTEGER             :: k,i

    CALL point2constants( C )

    CALL snowcheck (dels, ssnow, soil, met )

    CALL snowdensity (dels, ssnow, soil)

    CALL snow_accum (dels, canopy, met, ssnow, soil )

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add snow melt to global snow melt variable:
    ssnow%smelt(:) = snowmlt(:)
    ! Adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    CALL snowl_adjust(dels, ssnow, canopy )

    CALL GWstempv(dels, canopy, ssnow, soil)

    !do the soil and snow melting, freezing prior to water movement
    DO i=1,mp
       ssnow%tss(i) =  (1-ssnow%isflag(i))*ssnow%tgg(i,1) + ssnow%isflag(i)*ssnow%tggsn(i,1)
    END DO

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add new snow melt to global snow melt variable:
    ssnow%smelt(:) = ssnow%smelt(:) + snowmlt(:)

  END SUBROUTINE snow_processes_soil_thermal

  SUBROUTINE GWstempv(dels, canopy, ssnow, soil)
    USE cable_common_module, ONLY: cable_user
    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(canopy_type),    INTENT(INOUT) :: canopy
    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL, DIMENSION(mp) ::                                                      &
         coefa, coefb,  & !
         sgamm            !

    REAL, DIMENSION(mp) ::                                                 &
         dtg,     & !
         ew,      & !
         xx,     & !
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

    at = 0.0
    bt = 1.0
    ct = 0.0
    coeff = 0.0

    IF (cable_user%soil_thermal_fix) THEN

       ccnsw = total_soil_conductivity(ssnow,soil)

    ELSE

       ccnsw = old_soil_conductivity(ssnow,soil)
    END IF

    xx(:) = 0.

    WHERE(ssnow%isflag == 0)
       xx(:) = MAX( 0., ssnow%snowd / ssnow%ssdnn )
       ccnsw(:,1) = ( ccnsw(:,1) - 0.2 ) * ( soil%zse(1) / ( soil%zse(1) + xx(:) ) &
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
       coeff(:,2) = 2.0 / ( ( soil%zse(1) + xx(:) ) / ccnsw(:,1) + soil%zse(2) /   &
            ccnsw(:,2) )
       coefa = 0.0
       coefb = REAL( coeff(:,2) )

       wblfsp = ssnow%wblf(:,k)

       ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)), &
            ( 1.0 - soil%ssat_vec(:,k) ) * &
            soil%css_vec(:,k) * soil%rhosoil_vec(:,k)   &
            + soil%ssat_vec(:,k) * ( wblfsp * C%cs_rho_wat +            &
            ssnow%wbfice(:,k) * C%cs_rho_ice ) )     &
            * soil%zse_vec(:,k)

       ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + C%cgsnow * ssnow%snowd

       dtg = dels / ssnow%gammzz(:,k)

       at(:,k) = - dtg * coeff(:,k)
       ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
       bt(:,k) = 1.0 - at(:,k) - ct(:,k)

    END WHERE

    DO k = 2, ms

       WHERE( ssnow%isflag == 0 )

          wblfsp = ssnow%wblf(:,k)

          ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)), &
               ( 1.0 - soil%ssat_vec(:,k) ) * &
               soil%css_vec(:,k) * soil%rhosoil_vec(:,k)   &
               + soil%ssat_vec(:,k) * ( wblfsp * C%cs_rho_wat +            &
               ssnow%wbfice(:,k) * C%cs_rho_ice ) )     &
               * soil%zse_vec(:,k)

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
          sgamm = ssnow%ssdn(:,k) * C%cgsnow * ssnow%sdepth(:,k)
          dtg = dels / sgamm
          at(:,k-3) = - dtg * coeff(:,k-3)
          ct(:,k-3) = - dtg * coeff(:,k-2)
          bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
       END WHERE

    END DO

    DO k = 1, ms

       WHERE( ssnow%isflag /= 0 )
          wblfsp = ssnow%wblf(:,k)

          ssnow%gammzz(:,k) = MAX((soil%heat_cap_lower_limit(:,k)),&
               ( 1.0 - soil%ssat_vec(:,k) ) * soil%css_vec(:,k) *             &
               soil%rhosoil_vec(:,k) + soil%ssat_vec(:,k) * ( wblfsp * C%cs_rho_wat +&
               ssnow%wbfice(:,k) * C%cs_rho_ice)) * &
               soil%zse_vec(:,k)

          dtg = dels / ssnow%gammzz(:,k)
          at(:,k) = - dtg * coeff(:,k)
          ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
          bt(:,k) = 1.0 - at(:,k) - ct(:,k)

       END WHERE

    END DO

    WHERE( ssnow%isflag /= 0 )
       sgamm = ssnow%ssdn(:,1) * C%cgsnow * ssnow%sdepth(:,1)

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

  END SUBROUTINE GWstempv

  ! ----------------------------------------------------------------------------
  SUBROUTINE calc_soil_root_resistance(ssnow, soil, veg, bgc, root_length, i)
     ! Calculate root & soil hydraulic resistance following SPA approach
     ! (Williams et al.)
     !
     ! Root hydraulic resistance declines linearly with increasing root
     ! biomass according to root resistivity (400) [MPA s m2 mmol-1].
     !
     ! Soil hydraulic resistance depends on soil conductivity, root length,
     ! depth of layer and distance between roots.
     !
     ! In units conversion, useful to recall that:
     ! m s-1 = m3 m-1 m-1 s-1
     ! m3 (amount of water) m-1 (per unit length) m-1 (per unit hydraulic head,
     !                                                 measured in meters) s-1
     !
     ! References:
     ! -----------
     ! * Duursma, R. A. 2008. Predicting the decline in daily maximum
     !   transpiration rate of two pine stands during drought based on
     !   constant minimum leaf water potential and plant hydraulic conductance.
     !   Tree Physiology, 28, 265–276.
     ! * Gardner, W.R. 1964. Relation of root distribution to water uptake
     !   and availability. Agron. J. 56:41–45.
     ! * Newman, E.I. 1969. Resistance to water flow in soil and plant. I.
     !   Soil resistance in relation to amounts of root: theoretical
     !   estimates. J. Appl. Ecol. 6:1–12.
     ! * Williams, M. et al. 1996. Modeling the soil–plant–atmosphere continuum
     !   in a Quercus–Acer stand at Harvard Forest: the regulation of stomatal
     !   conductance by light, nitrogen and soil/plant hydraulic properties.
     !   Plant Cell Environ. 19:911–927.
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_def_types_mod
     USE cable_common_module

     IMPLICIT NONE

     TYPE (soil_snow_type), INTENT(INOUT)        :: ssnow
     TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
     TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
     TYPE (bgc_pool_type),  INTENT(IN)           :: bgc

     ! All from Williams et al. 2001, Tree phys
     REAL, PARAMETER :: pi = 3.1415927
     REAL, PARAMETER :: root_radius = 0.0005                 ! m
     REAL, PARAMETER :: root_xsec_area = pi * root_radius**2 ! m2
     REAL, PARAMETER :: root_density = 0.5e6                ! g biomass m-3 root
     REAL, PARAMETER :: root_resistivity = 25.             ! MPa s g mmol-1, Bonan
     REAL, PARAMETER :: root_k = 100.0

     ! unit conv
     REAL, PARAMETER :: head = 0.009807             ! head of pressure  (MPa/m)
     REAL, PARAMETER :: MM_TO_M = 0.001
     REAL, PARAMETER :: KPA_2_MPa = 0.001
     REAL, PARAMETER :: M_HEAD_TO_MPa = 9.8 * KPA_2_MPa
     REAL, PARAMETER :: G_WATER_TO_MOLE = 1.0 / 18.01528
     REAL, PARAMETER :: CUBIC_M_WATER_2_GRAMS = 1E6
     REAL, PARAMETER :: MOL_2_MMOL = 1000.0
     REAL, PARAMETER :: TINY_NUMBER = 1E-35
     REAL, PARAMETER :: HUGE_NUMBER = 1E35
     REAL, PARAMETER :: BIG_NUMBER = 1E9

     REAL, DIMENSION(ms) :: depth
     REAL                :: root_mass, rs, Ksoil, root_biomass, root_depth
     REAL                :: soil_resist, rsum, conv


     REAL, DIMENSION(:), INTENT(INOUT) :: root_length
     ! ratio Dry matter mass to g(C)
     REAL, PARAMETER                   :: gC2DM = 1./0.49

     INTEGER, INTENT(IN) :: i
     INTEGER :: j
     INTEGER, PARAMETER :: ROOT_INDEX = 3

     REAL               :: ht, stem_biomass
     REAL, PARAMETER    :: Kbiometric = 50.0 ! cst in height-diameter relationship
     REAL, PARAMETER    :: WD = 300.0 ! Wood density kgC/m3
     INTEGER, PARAMETER :: STEM_INDEX = 2

     stem_biomass = bgc%cplant(i,STEM_INDEX) * gC2DM
     ht = (Kbiometric**(3.0/4.0))*(4.*stem_biomass/(WD*PI))**(1.0/4.0)
     !print*, ht

     ! convert from gC to g biomass, i.e. twice the C content
     root_biomass = bgc%cplant(i,ROOT_INDEX) * gC2DM

     ! Always provide a minimum root biomass
     root_biomass = MAX(5., root_biomass)

     root_length = 0.0
     DO j = 1, ms ! Loop over 6 soil layers

        ! Root biomass density (g biomass m-3 soil)
        ! Divide root mass up by the frac roots in the layer (g m-3)
        ! plant carbon is g C m-2
        root_mass = root_biomass * veg%froot(i,j)

        ! Root length density (m root m-3 soil)
        root_length(j) = root_mass / (root_density * root_xsec_area)

     END DO

     ! Store each layers resistance, used in LWP calculatons
     rsum = 0.0

     DO j = 1, ms ! Loop over 6 soil layers

        ! Soil Hydraulic conductivity (m s-1), Campbell 1974
        Ksoil = soil%hyds(i) * (ssnow%wb(i,j) / &
                  soil%ssat(i))**(2.0 * soil%bch(i) + 3.0)

        ! converts from m s-1 to m2 s-1 MPa-1
        Ksoil = Ksoil / head

        ! prevent floating point error
        IF (Ksoil < TINY_NUMBER) THEN
           ssnow%soilR(i,j) = HUGE_NUMBER
        ELSE
           ! Conductance of the soil-to-root pathway can be estimated
           ! assuming that the root system consists of one long root that
           ! has access to a surrounding cylinder of soil
           ! (Gardner 1960, Newman 1969)
           rs = SQRT(1.0 / (root_length(j) * pi))

           ! Soil-to-root resistance (MPa s m2 mmol-1 H2O)
           soil_resist = LOG(rs / root_radius) / &
                              (2.0 * pi * root_length(j) * soil%zse(j) * Ksoil)

           ! convert from MPa s m2 m-3 to MPa s m2 mmol-1
           soil_resist = soil_resist * 1E-6 * 18. * 0.001

           ! MPa s m2 mmol-1 H2O
           ! root_resistance is commented out : don't use root-component of
           ! resistance (is part of plant resistance)
           ssnow%soilR(i,j) = soil_resist !+ root_resist
        END IF

        IF (ssnow%soilR(i,j) .GT. 0.0) THEN
           ! Need to combine resistances in parallel, but we only want the
           ! soil term as the root component is part of the plant resistance
           rsum = rsum + ( 1.0 / ssnow%soilR(i,j) )
        ENDIF

     END DO
     ssnow%tot_bg_resist(i) = 1.0 / rsum

  END SUBROUTINE calc_soil_root_resistance
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  SUBROUTINE calc_swp(ssnow, soil, i)
     ! Calculate the soil water potential.
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_def_types_mod
     USE cable_common_module

     IMPLICIT NONE

     TYPE (soil_snow_type), INTENT(INOUT)        :: ssnow
     TYPE (soil_parameter_type), INTENT(INOUT)   :: soil

     INTEGER             :: j
     INTEGER, INTENT(IN) :: i
     REAL                :: psi_sat, t_over_t_sat, cond_per_layer
     REAL, PARAMETER     :: sucmin = -1E5 ! minimum soil pressure head [m]

     REAL, PARAMETER :: KPA_2_MPa = 0.001
     REAL, PARAMETER :: M_HEAD_TO_MPa = 9.8 * KPA_2_MPa

     ssnow%psi_soil(i,:) = 0.0 ! MPa

     ! Soil matric potential at saturation (m of head to MPa: 9.81 * KPA_2_MPA)
     psi_sat = soil%sucs(i) * M_HEAD_TO_MPa

     DO j = 1, ms ! Loop over 6 soil layers

        ! Below the wilting point (-1.5 MPa) the water potential drops to
       ! silly value. This really only an issue for the v.top
       ! two layers and has negligble impact on the weighted psi_soil which is
       ! what is used anyway
       t_over_t_sat = MAX(1.0e-9, MIN(1.0, ssnow%wb(i,j) / soil%ssat(i)))
       ssnow%psi_soil(i,j) = psi_sat * t_over_t_sat**(-soil%bch(i))
       ssnow%psi_soil(i,j) = MAX(MIN(ssnow%psi_soil(i,j), soil%sucs(i)), sucmin)

    END DO

  END SUBROUTINE calc_swp
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  SUBROUTINE calc_weighted_swp_and_frac_uptake(ssnow, soil, canopy, &
                                               root_length, i)
     !
     ! Determine weighted SWP given the the maximum rate of water supply from
     ! each rooted soil layer and hydraulic resistance of each layer. We are
     ! also calculating a weighting fraction for water extraction. This is
     ! achieved by roughly estimating the maximum rate of water supply from each
     ! rooted soil layer, using SWP and hydraulic resistance of each layer.
     ! Actual water from each layer is determined using the estimated value as a
     ! weighted factor.
     !
     ! Martin De Kauwe, 4th March, 2019

     USE cable_def_types_mod
     USE cable_common_module

     IMPLICIT NONE

     TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow
     TYPE (soil_parameter_type), INTENT(INOUT) :: soil
     TYPE(canopy_type), INTENT(INOUT)          :: canopy ! vegetation variables

     REAL, PARAMETER :: MM_TO_M = 0.001
     REAL, PARAMETER :: KPA_2_MPa = 0.001
     REAL, PARAMETER :: M_HEAD_TO_MPa = 9.8 * KPA_2_MPa

     ! The minimum root water potential (MPa), used in determining fractional
     ! water uptake in soil layers
     REAL, PARAMETER :: min_root_wp = -3

     REAL, DIMENSION(ms)            :: swp, est_evap
     REAL, DIMENSION(:), INTENT(IN) :: root_length
     REAL                           :: total_est_evap, swp_diff

     INTEGER, INTENT(IN) :: i
     INTEGER             :: j

     ! SPA method to figure out relative water uptake.
     LOGICAL :: SPA_relative_uptake
     SPA_relative_uptake = .TRUE.

     total_est_evap = 0.0
     est_evap = 0.0
     ssnow%weighted_psi_soil(i) = 0.0
     ssnow%fraction_uptake(i,:) = 0.0

     ! Estimate max transpiration from gradient-gravity / soil resistance
     DO j = 1, ms ! Loop over 6 soil layers

        IF (ssnow%soilR(i,j) .GT. 0.0) THEN
           est_evap(j) = MAX(0.0, &
                        (ssnow%psi_soil(i,j) - min_root_wp) / ssnow%soilR(i,j))
        ELSE
           est_evap(j) = 0.0 ! when no roots present
        ENDIF
        ! NEED TO ADD SOMETHING IF THE SOIL IS FROZEN, what is ice in CABLE?
        !IF ( iceprop(i) .gt. 0. ) THEN
        !  est_evap(i) = 0.0
        !ENDIF

        ! Soil water potential weighted by layer Emax (from SPA)
        ssnow%weighted_psi_soil(i) = ssnow%weighted_psi_soil(i) + &
                                     ssnow%psi_soil(i,j) * est_evap(j)
     END DO
     total_est_evap = SUM(est_evap)

     ! calculate the weighted psi_soil
     IF (total_est_evap > 0.0) THEN
        ! Soil water potential is weighted by total_est_evap.
        ssnow%weighted_psi_soil(i) = ssnow%weighted_psi_soil(i) / total_est_evap
     ELSE
        ssnow%weighted_psi_soil(i) = 0.0
        DO j = 1, ms ! Loop over 6 soil layers
           ssnow%weighted_psi_soil(i) = ssnow%weighted_psi_soil(i) + &
                                          ssnow%psi_soil(i,j) * soil%zse(j)
        END DO
        ssnow%weighted_psi_soil(i) = ssnow%weighted_psi_soil(i) / SUM(soil%zse)
     END IF

     ! SPA method to figure out relative water uptake.
     ! Fraction uptake in each layer by Emax in each layer
     IF (SPA_relative_uptake) THEN

        IF (total_est_evap > 0.0) THEN
           DO j = 1, ms ! Loop over 6 soil layers

              ! fraction of water taken from layer, I've lower bounded frac
              ! uptake because when soilR is set to a huge number
              ! (see calc_soil_root_resistance), then frac_uptake will be so
              ! small you end up with numerical issues.
              ssnow%fraction_uptake(i,j) = MAX(1E-09, &
                                             est_evap(j) / total_est_evap)

              IF ((ssnow%fraction_uptake(i,j) > 1.0) .OR. &
                  (ssnow%fraction_uptake(i,j) < 0.0)) THEN
                 PRINT *, 'Problem with the uptake fraction (either >1 or 0<)'
                 STOP
              END IF

           END DO
        ELSE
           ! No water was evaporated
           ssnow%fraction_uptake(i,:) = 1.0 / FLOAT(ms)
        END IF

     ! Use Taylor-Keppler root water uptake distribution.
     ELSE
        ! Taylor and Keppler: relative water uptake is
        ! proportional to root length density and Psi difference.
        ! See : Taylor, H.M. and B. Keppler. 1975. Water uptake by cotton root
        ! systems: an examination of assumptions in the single root model.
        ! Soil Science. 120:57-67.
        DO j = 1, ms ! Loop over 6 soil layers

           IF (total_est_evap .GT. 0.) THEN
              swp_diff = MAX(0., (ssnow%psi_soil(i,j) - min_root_wp))
              ssnow%fraction_uptake(i,j) = root_length(j) * swp_diff
           ELSE
              ! no water uptake possible
              ssnow%fraction_uptake(i,j) = 0.0
           END IF
        END DO

        IF (SUM(ssnow%fraction_uptake) .GT. 0) THEN
           ! Make sure that it sums to 1.
           ssnow%fraction_uptake = ssnow%fraction_uptake / &
                                      SUM(ssnow%fraction_uptake)
        ELSE
           ssnow%fraction_uptake = 0.0
        ENDIF

     ENDIF

  END SUBROUTINE calc_weighted_swp_and_frac_uptake
  ! ----------------------------------------------------------------------------

END MODULE cable_soil_snow_module
