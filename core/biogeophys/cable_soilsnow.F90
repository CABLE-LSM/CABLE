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
!
! ==============================================================================

MODULE cable_soil_snow_module

   USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             r_2, ms, mp
   USE cable_data_module, ONLY : issnow_type, point2constants

   IMPLICIT NONE

   PRIVATE

   TYPE ( issnow_type ) :: C

   REAL, PARAMETER :: &
      cgsnow = 2090.0,     & ! specific heat capacity for snow
      csice = 2.100e3,     & ! specific heat capacity for ice
      cswat = 4.218e3,     & ! specific heat capacity for water
      rhowat = 1000.0,     & ! density of water
      snmin = 1.,          & ! for 3-layer;
      max_ssdn = 750.0,    & !
      max_sconds = 2.51
      
   REAL(r_2), PARAMETER :: frozen_limit = 0.85_r_2  ! EAK Feb2011 (could be 0.95)

   REAL :: cp    ! specific heat capacity for air

   !jhan:make parameter
   REAL :: max_glacier_snowd

   ! This module contains the following subroutines:
   PUBLIC :: soil_snow ! must be available outside this module
   PRIVATE :: snowdensity, snow_melting, snowcheck, snowl_adjust
   PRIVATE :: trimb, smoisturev, snow_accum, stempv
   PRIVATE :: remove_trans, soilfreeze_serial ! , soilfreeze

CONTAINS

! SUBROUTINE trimb
!
!      this routine solves the system
!          a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kmax-1
!          with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
!          and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kmax
!
!        the Thomas algorithm is used for solving sets of linear equation
!        rhs initially contains rhs; leaves with answer (jlm)
!        n.b. this one does not assume b = 1-a-c
!
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

! SUBROUTINE smoisturev (fwtop,dels,ssnow,soil)
!      Solves implicit soil moisture equation
!      Science development by Eva Kowalczyk and John McGregor, CMAR
!
SUBROUTINE smoisturev(dels, ssnow, soil)

   USE cable_common_module

   REAL, INTENT(IN) :: dels    ! time step size (s)
   TYPE(soil_snow_type),      INTENT(INOUT) ::                                &
      ssnow ! soil and snow variables
   TYPE(soil_parameter_type), INTENT(INOUT) ::                                &
      soil  ! soil parameters

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
      wbficemx

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

   INTEGER :: k

   at = 0.0
   bt = 1.0
   ct = 0.0
   z1mult(:,1) = 0.0       ! corresponds to 2b+3
   z1mult(:,ms+1) = 0.0    ! corresponds to 2b+3
   z1(:,1) = 0.0           ! i.e. K(.5),    value at surface
   z1(:,ms+1) = 0.0        ! i.e. K(ms+.5), value at bottom

   ! nmeth: equation solution technique
   IF (nmeth <= 0) THEN

      ! jlm split TVD version
      ! all land points
      delt(:,0) = 0.0
      fluxh(:,0) = 0.0
      fluxh(:,ms) = 0.0

      DO k = 1, ms-1

         ! Calculate amount of liquid soil water:
         IF (.not. cable_user%l_new_runoff_speed) THEN
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

      IF (.not. cable_user%l_new_runoff_speed) then

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
         speed_k = MIN( speed_k, real(0.5 * soil%zse(ms) / dels, r_2) )
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
      ssnow%wblf(:,1) = ssnow%wblf(:,1) + dtt(:,1) * ssnow%fwtop1 / rhowat
      ssnow%wblf(:,2) = ssnow%wblf(:,2) + dtt(:,2) * ssnow%fwtop2 / rhowat
      ssnow%wblf(:,3) = ssnow%wblf(:,3) + dtt(:,3) * ssnow%fwtop3 / rhowat

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

         wbficemx = MAX( wbficemx, real(ssnow%wbfice(:,k)) )
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

      ssnow%wblf(:,1) = ssnow%wblf(:,1) + dtt(:,1) * ssnow%fwtop1 / rhowat
      ssnow%wblf(:,2) = ssnow%wblf(:,2) + dtt(:,2) * ssnow%fwtop2 / rhowat
      ssnow%wblf(:,3) = ssnow%wblf(:,3) + dtt(:,3) * ssnow%fwtop3 / rhowat

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
   REAL, DIMENSION(mp,3) :: dels_ssdn !, ssnow_tgg_min

   !MC? ssnow_isflag_ssdn = SPREAD( ssnow%isflag,2,mp)
   ssnow_isflag_ssdn = SPREAD(ssnow%isflag,2,3)

   !MC? dels_ssdn = SPREAD( SPREAD( dels, 1, mp ), 2,  mp )
   dels_ssdn = SPREAD(SPREAD(dels, 1, mp), 2,  3)
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

      IF( ssnow%snowd(j) > 0.0 .AND. ssnow%isflag(j) == 0 &
          .AND. ssnow%tgg(j,1) >= C%TFRZ ) THEN

         ! snow covered land
         ! following done in sflux  via  ga= ... +cls*egg + ...
         ! ** land,snow,melting
         snowflx(j) = REAL((ssnow%tgg(j,1) - C%TFRZ) * ssnow%gammzz(j,1))

         ! prevent snow depth going negative
         snowmlt(j) = MIN(snowflx(j) / C%HLF, ssnow%snowd(j) )

         ssnow%dtmlt(j,1) = ssnow%dtmlt(j,1) + snowmlt(j) * C%HLF &
                            / real(ssnow%gammzz(j,1))

         ssnow%snowd(j) = ssnow%snowd(j) - snowmlt(j)
         ssnow%tgg(j,1) = ssnow%tgg(j,1) - snowmlt(j) * &
                          C%HLF / real(ssnow%gammzz(j,1))
      ENDIF

   END DO

   smelt1(:,0) = 0.0

   DO k = 1, 3

      !where there is snow
      WHERE( ssnow%snowd > 0.0 .AND. ssnow%isflag > 0 )

         sgamm = ssnow%ssdn(:,k) * cgsnow * ssnow%sdepth(:,k)

         ! snow melt refreezing
         snowflx = smelt1(:,k-1) * C%HLF / dels

         ssnow%tggsn(:,k) = ssnow%tggsn(:,k) + ( snowflx * dels +              &
                            smelt1(:,k-1)*cswat *( C%TFRZ-ssnow%tggsn(:,k) ) ) &
                            / ( sgamm + cswat*smelt1(:,k-1) )

         ! increase density due to snowmelt
         osm = ssnow%smass(:,k)
         ssnow%smass(:,k) = ssnow%smass(:,k) + smelt1(:,k-1)
         ssnow%ssdn(:,k) = MAX( 120.0, MIN( ssnow%ssdn(:,k) * osm /            &
                           ssnow%smass(:,k) + rhowat * ( 1.0 - osm /           &
                           ssnow%smass(:,k)), max_ssdn ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,k) = MIN( 450.0, ssnow%ssdn(:,k) )

         ssnow%sdepth(:,k) = ssnow%smass(:,k) / ssnow%ssdn(:,k)

         sgamm = ssnow%smass(:,k) * cgsnow

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
      xxx        !

   WHERE (canopy%precis > 0.0 .and. ssnow%isflag == 0)
      ! accumulate solid part
      ssnow%snowd = MAX( ssnow%snowd + met%precip_sn, 0.0 )

      canopy%precis = canopy%precis - met%precip_sn

      ssnow%ssdn(:,1) = MAX( 120.0, ssnow%ssdn(:,1)                            &
                        * ssnow%osnowd / MAX( 0.01, ssnow%snowd )              &
                        + 120.0 * met%precip_sn / MAX( 0.01, ssnow%snowd ) )

      ssnow%ssdnn = ssnow%ssdn(:,1)

      WHERE( canopy%precis > 0.0 .AND. ssnow%tgg(:,1) < C%TFRZ )

         ssnow%snowd = MAX(ssnow%snowd + canopy%precis, 0.0)

         ssnow%tgg(:,1) = ssnow%tgg(:,1) + canopy%precis * C%HLF               &
                          / ( REAL( ssnow%gammzz(:,1) ) + cswat *canopy%precis )
         ! change density due to water being added
         ssnow%ssdn(:,1) = MIN( max_ssdn, MAX( 120.0, ssnow%ssdn(:,1)          &
                           * ssnow%osnowd / MAX( 0.01, ssnow%snowd ) + rhowat  &
                           * canopy%precis / MAX( 0.01, ssnow%snowd )  ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,1) = MIN( 450.0, ssnow%ssdn(:,1) )

         canopy%precis = 0.0
         ssnow%ssdnn = ssnow%ssdn(:,1)

      END WHERE

   END WHERE ! (canopy%precis > 0. .and. ssnow%isflag == 0)

   WHERE (canopy%precis > 0.0 .and.  ssnow%isflag > 0)

      ! add solid precip
      ssnow%snowd = MAX( ssnow%snowd + met%precip_sn, 0.0 )

      canopy%precis = canopy%precis - met%precip_sn  ! remaining liquid precip

      ! update top snow layer with fresh snow
      osm = ssnow%smass(:,1)
      ssnow%smass(:,1) = ssnow%smass(:,1) + met%precip_sn
      ssnow%ssdn(:,1) = MAX( 120.0,ssnow%ssdn(:,1) * osm / ssnow%smass(:,1)    &
                        + 120.0 * met%precip_sn / ssnow%smass(:,1) )

      ssnow%sdepth(:,1) = MAX( 0.02, ssnow%smass(:,1) / ssnow%ssdn(:,1) )

      ! add liquid precip
      WHERE( canopy%precis > 0.0 )

         ssnow%snowd = MAX( ssnow%snowd + canopy%precis, 0.0 )
         sgamm = ssnow%ssdn(:,1) * cgsnow * ssnow%sdepth(:,1)
         osm = ssnow%smass(:,1)

         ssnow%tggsn(:,1) = ssnow%tggsn(:,1) + canopy%precis * C%HLF           &
                            * osm / (sgamm * ssnow%osnowd )
         ssnow%smass(:,1) = ssnow%smass(:,1) + canopy%precis                   &
                            * osm/ssnow%osnowd

         ssnow%ssdn(:,1) = MAX( 120.0, MIN( ssnow%ssdn(:,1) * osm /            &
                           ssnow%smass(:,1) +  rhowat *                        &
                           ( 1.0 - osm / ssnow%smass(:,1) ), max_ssdn ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,1) = MIN( 450.0, ssnow%ssdn(:,1) )

         ssnow%sdepth(:,1) = ssnow%smass(:,1)/ssnow%ssdn(:,1)

         !layer 2
         sgamm = ssnow%ssdn(:,2) * cgsnow * ssnow%sdepth(:,2)
         osm = ssnow%smass(:,2)
         ssnow%tggsn(:,2) = ssnow%tggsn(:,2) + canopy%precis * C%HLF           &
                            * osm / ( sgamm * ssnow%osnowd )
         ssnow%smass(:,2) = ssnow%smass(:,2) + canopy%precis                   &
                            * osm / ssnow%osnowd
         ssnow%ssdn(:,2) = MAX( 120.0, MIN( ssnow%ssdn(:,2) * osm /            &
                           ssnow%smass(:,2) + rhowat *                         &
                           ( 1.0 - osm / ssnow%smass(:,2) ), max_ssdn ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,2) = MIN( 450.0, ssnow%ssdn(:,2) )

         ssnow%sdepth(:,2) = ssnow%smass(:,2) / ssnow%ssdn(:,2)

         !layer 3
         sgamm = ssnow%ssdn(:,3) * cgsnow * ssnow%sdepth(:,3)
         osm = ssnow%smass(:,3)
         ssnow%tggsn(:,3) = ssnow%tggsn(:,3) + canopy%precis * C%HLF           &
                            * osm / ( sgamm * ssnow%osnowd )
         ssnow%smass(:,3) = ssnow%smass(:,3) + canopy%precis                   &
                            * osm / ssnow%osnowd
        ssnow%ssdn(:,3) = MAX( 120.0, MIN( ssnow%ssdn(:,3) * osm /             &
                          ssnow%smass(:,3) + rhowat *                          &
                          ( 1.0 - osm / ssnow%smass(:,3) ), max_ssdn ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,3) = MIN(450.0,ssnow%ssdn(:,3))

         ssnow%sdepth(:,3) = ssnow%smass(:,3) / ssnow%ssdn(:,3)

         canopy%precis = 0.0

      END WHERE

   END WHERE


   ! 'fess' is for soil evap and 'fes' is for soil evap plus soil puddle evap
   ! canopy%segg = real(canopy%fess) / C%HL
   canopy%segg = real(canopy%fess + canopy%fes_cor) / C%HL

   ! Initialise snow evaporation:
   ssnow%evapsn = 0

   ! Snow evaporation and dew on snow
   WHERE (ssnow%snowd > 0.1)

      ssnow%evapsn = dels * real(canopy%fess + canopy%fes_cor) / ( C%HL + C%HLF )
      xxx = ssnow%evapsn

      WHERE ((ssnow%isflag == 0) .AND. ((canopy%fess+canopy%fes_cor).GT. 0.0_r_2) ) &
         ssnow%evapsn = MIN( ssnow%snowd, xxx )

      WHERE ((ssnow%isflag  > 0) .AND. ((canopy%fess+canopy%fes_cor) .GT. 0.0_r_2) ) &
         ssnow%evapsn = MIN( 0.9 * ssnow%smass(:,1), xxx )

      ssnow%snowd = ssnow%snowd - ssnow%evapsn

      WHERE( ssnow%isflag > 0 )
         ssnow%smass(:,1) = ssnow%smass(:,1)  - ssnow%evapsn
         ssnow%sdepth(:,1) = MAX( 0.02, ssnow%smass(:,1) / ssnow%ssdn(:,1) )
      END WHERE

      canopy%segg = 0.0

   END WHERE

END SUBROUTINE snow_accum

! -----------------------------------------------------------------------------

SUBROUTINE surfbv(dels, ssnow, soil, veg)

   USE cable_common_module

   REAL, INTENT(IN) :: dels ! integration time step (s)
   TYPE(soil_snow_type), INTENT(INOUT) :: ssnow  ! soil+snow variables
   TYPE(soil_parameter_type), INTENT(INOUT)  :: soil  ! soil parameters
   TYPE(veg_parameter_type),  INTENT(IN)     :: veg

   !jhan:cable.nml
   INTEGER, PARAMETER      :: nglacier = 2 ! 0 original, 1 off, 2 new Eva

   REAL, DIMENSION(mp) ::                                                      &
      rnof5,      & !
      sgamm,      & !
      smasstot
   REAL(r_2), DIMENSION(mp) :: xxx

   REAL, DIMENSION(mp,0:3) :: smelt1

   INTEGER :: k

   CALL smoisturev(dels, ssnow, soil)

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
            sgamm = ssnow%ssdn(:,k) * cgsnow * ssnow%sdepth(:,k)
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
      xxx = MAX(0.0_r_2, (ssnow%wb(:,ms) - real(soil%sfc(:),r_2))*soil%zse(ms)*1000.0)
      ssnow%sinfil  = MIN( real(xxx), ssnow%wb_lake )
      ssnow%wb(:,ms) = ssnow%wb(:,ms) - real(ssnow%sinfil / (soil%zse(ms)*1000.0), r_2)
      ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
      xxx = MAX(0.0_r_2, (ssnow%wb(:,ms) - 0.5*(soil%sfc + soil%swilt))*soil%zse(ms)*1000.0)
      ssnow%sinfil  = MIN( real(xxx), ssnow%wb_lake )
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

   INTEGER   :: j,k
   REAL      :: snow_ccnsw
   REAL(r_2) :: exp_arg
   LOGICAL   :: direct2min = .FALSE.


   at         = 0.0_r_2
   bt         = 1.0_r_2
   ct         = 0.0_r_2
   coeff      = 0.0_r_2
   snow_ccnsw = 2.0

   DO k = 1, ms

      DO j = 1, mp

         IF ( soil%isoilm(j) == 9 ) THEN
            ! permanent ice: fix hard-wired number in next version
            ccnsw(j,k) = snow_ccnsw
         ELSE
            ew(j) = ssnow%wblf(j,k) * real(soil%ssat(j),r_2)
            exp_arg = (ew(j) * LOG(60.0_r_2)) + (ssnow%wbfice(j,k) &
                      * real(soil%ssat(j),r_2) * LOG(250.0_r_2))

            IF ( exp_arg > 30.0_r_2 ) direct2min = .TRUE.

            IF ( direct2min) THEN

               ccnsw(j,k) = 1.5_r_2 * MAX( 1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5_r_2 * &
                            real(soil%ssat(j),r_2) / &
                            MIN( ew(j), 0.5_r_2 * real(soil%ssat(j),r_2) ) ) ) )

            ELSE

               ccnsw(j,k) = MIN( soil%cnsd(j) * EXP(exp_arg), 1.5_r_2 ) &
                            * MAX( 1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5_r_2 * &
                            real(soil%ssat(j),r_2) / &
                            MIN( ew(j), 0.5_r_2 * real(soil%ssat(j),r_2) ) ) ) )

            ENDIF

            direct2min = .FALSE.

         ENDIF

      END DO

   END DO

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

      xx = soil%css * soil%rhosoil

      ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil   &
                          + soil%ssat * ( wblfsp * cswat * rhowat +            &
                          ssnow%wbfice(:,k) * csice * rhowat * 0.9 ), xx )     &
                          * soil%zse(k)

      ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + cgsnow * ssnow%snowd

      dtg = dels / ssnow%gammzz(:,k)

      at(:,k) = - dtg * coeff(:,k)
      ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
      bt(:,k) = 1.0 - at(:,k) - ct(:,k)

   END WHERE

   DO k = 2, ms

      WHERE( ssnow%isflag == 0 )

         wblfsp = ssnow%wblf(:,k)
         xx = soil%css * soil%rhosoil

         ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil&
                             + soil%ssat * ( wblfsp * cswat * rhowat +         &
                             ssnow%wbfice(:,k) * csice * rhowat * 0.9 ), xx )  &
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
         sgamm = ssnow%ssdn(:,k) * cgsnow * ssnow%sdepth(:,k)
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
                             soil%rhosoil + soil%ssat * ( wblfsp * cswat *     &
                             rhowat + ssnow%wbfice(:,k) * csice * rhowat *     &
                             0.9) , xx ) * soil%zse(k)

         dtg = dels / ssnow%gammzz(:,k)
         at(:,k) = - dtg * coeff(:,k)
         ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
         bt(:,k) = 1.0 - at(:,k) - ct(:,k)

      END WHERE

   END DO

   WHERE( ssnow%isflag /= 0 )
      sgamm = ssnow%ssdn(:,1) * cgsnow * ssnow%sdepth(:,1)

      bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm

      ssnow%tggsn(:,1) = ssnow%tggsn(:,1) + ( canopy%ga - ssnow%tggsn(:,1 )    &
                         * REAL( canopy%dgdtg ) ) * dels / sgamm

      rhs(:,1-3) = ssnow%tggsn(:,1)
   END WHERE


   !     note in the following that tgg and tggsn are processed together
   tmp_mat(:,1:3)      = REAL(ssnow%tggsn,r_2)
   tmp_mat(:,4:(ms+3)) = REAL(ssnow%tgg,r_2)

   CALL trimb( at, bt, ct, tmp_mat, ms + 3 )

   ssnow%tggsn = REAL( tmp_mat(:,1:3) )
   ssnow%tgg   = REAL( tmp_mat(:,4:(ms+3)) )
   canopy%sghflux = coefa * ( ssnow%tggsn(:,1) - ssnow%tggsn(:,2) )
   canopy%ghflux  = coefb * ( ssnow%tgg(:,1) - ssnow%tgg(:,2) ) ! +ve downwards

END SUBROUTINE stempv

! -----------------------------------------------------------------------------

SUBROUTINE snowcheck(ssnow, soil)

   USE cable_common_module

   TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters

   INTEGER :: j

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

SUBROUTINE snowl_adjust(ssnow)

   TYPE(soil_snow_type), INTENT(INOUT) :: ssnow

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

! not used
! SUBROUTINE soilfreeze(dels, soil, ssnow)
!    USE cable_common_module
!    REAL, INTENT(IN)                    :: dels ! integration time step (s)
!    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
!    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
!    REAL(r_2), DIMENSION(mp)           :: sicefreeze
!    REAL(r_2), DIMENSION(mp)           :: sicemelt
!    REAL, DIMENSION(mp)           :: xx
!    INTEGER k

!    xx = 0.
!    DO k = 1, ms

!       WHERE (ssnow%tgg(:,k) < C%TFRZ &
!           & .AND. frozen_limit * ssnow%wb(:,k) - ssnow%wbice(:,k) > .001)

!          sicefreeze = MIN( MAX( 0.0_r_2, ( frozen_limit * ssnow%wb(:,k) -      &
!                       ssnow%wbice(:,k) ) ) * soil%zse(k) * 1000.0,             &
!                       ( C%TFRZ - ssnow%tgg(:,k) ) * ssnow%gammzz(:,k) / C%HLF )
!          ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k) + sicefreeze / (soil%zse(k)  &
!                             * 1000.0), frozen_limit * ssnow%wb(:,k) )
!          xx = soil%css * soil%rhosoil
!          ssnow%gammzz(:,k) = MAX(                                              &
!              REAL((1.0 - soil%ssat) * soil%css * soil%rhosoil ,r_2)            &
!              + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(cswat * rhowat,r_2)   &
!              + ssnow%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2),              &
!              REAL(xx,r_2)) * REAL( soil%zse(k),r_2 )

!          WHERE (k == 1 .AND. ssnow%isflag == 0)
!             ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + cgsnow * ssnow%snowd
!          END WHERE
!          ssnow%tgg(:,k) = ssnow%tgg(:,k) + REAL(sicefreeze)                    &
!                           * C%HLF / REAL(ssnow%gammzz(:,k) )

!       ELSEWHERE( ssnow%tgg(:,k) > C%TFRZ .AND. ssnow%wbice(:,k) > 0. )

!          sicemelt = MIN( ssnow%wbice(:,k) * soil%zse(k) * 1000.0,              &
!                     ( ssnow%tgg(:,k) - C%TFRZ ) * ssnow%gammzz(:,k) / C%HLF )

!          ssnow%wbice(:,k) = MAX( 0.0_r_2, ssnow%wbice(:,k) - sicemelt          &
!                             / (soil%zse(k) * 1000.0) )
!          xx = soil%css * soil%rhosoil
!          ssnow%gammzz(:,k) = MAX(                                              &
!               REAL((1.0-soil%ssat) * soil%css * soil%rhosoil,r_2)             &
!               + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(cswat*rhowat,r_2)   &
!               + ssnow%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2),            &
!               REAL(xx,r_2) ) * REAL(soil%zse(k),r_2)
!          WHERE (k == 1 .AND. ssnow%isflag == 0)
!             ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + cgsnow * ssnow%snowd
!          END WHERE
!          ssnow%tgg(:,k) = ssnow%tgg(:,k) - REAL(sicemelt)                     &
!                           * C%HLF / REAL(ssnow%gammzz(:,k))

!       END WHERE

!    END DO

! END SUBROUTINE soilfreeze

! same as soilfreeze but with do loop to catch underflows
SUBROUTINE soilfreeze_serial(soil, ssnow)
  
  USE cable_common_module
  
   TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil

   REAL(r_2), DIMENSION(mp) :: sicefreeze
   REAL(r_2), DIMENSION(mp) :: sicemelt
   REAL,      DIMENSION(mp) :: xx
   real(r_2), dimension(ms) :: zse
   INTEGER :: k, i

   zse = real(soil%zse,r_2)
   
   xx = 0.
   DO k = 1, ms
      do i=1, mp
         if (ssnow%tgg(i,k) < C%TFRZ &
              .AND. frozen_limit * ssnow%wb(i,k) - ssnow%wbice(i,k) > 0.001_r_2) then
            sicefreeze(i) = MIN( MAX( 0.0_r_2, ( frozen_limit * ssnow%wb(i,k) - &
                 ssnow%wbice(i,k) ) ) * zse(k) * 1000.0_r_2, &
                 ( C%TFRZ - ssnow%tgg(i,k) ) * ssnow%gammzz(i,k) / C%HLF )
            ssnow%wbice(i,k) = MIN( ssnow%wbice(i,k) + sicefreeze(i) / (zse(k) &
                 * 1000.0_r_2), frozen_limit * ssnow%wb(i,k) )
            xx(i) = soil%css(i) * soil%rhosoil(i)
            ssnow%gammzz(i,k) = MAX( &
                 REAL((1.0 - soil%ssat(i)) * soil%css(i) * soil%rhosoil(i), r_2) &
                 + (ssnow%wb(i,k) - ssnow%wbice(i,k)) * REAL(cswat * rhowat, r_2) &
                 + ssnow%wbice(i,k) * REAL(csice * rhowat * 0.9, r_2), &
                 REAL(xx(i),r_2)) * zse(k)
            if (k == 1 .AND. ssnow%isflag(i) == 0) ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + real(cgsnow * ssnow%snowd(i),r_2)
            if (sicefreeze(i)/ssnow%gammzz(i,k) > real(tiny(1.0),r_2)) &
                 ssnow%tgg(i,k) = ssnow%tgg(i,k) + REAL(sicefreeze(i)/ssnow%gammzz(i,k)) * C%HLF
         elseif (ssnow%tgg(i,k) > C%TFRZ .AND. ssnow%wbice(i,k) > 0._r_2 ) then
            sicemelt(i) = MIN( ssnow%wbice(i,k) * zse(k) * 1000.0_r_2, &
                 real(ssnow%tgg(i,k) - C%TFRZ,r_2) * ssnow%gammzz(i,k) / real(C%HLF,r_2) )
            ssnow%wbice(i,k) = MAX( 0.0_r_2, ssnow%wbice(i,k) - sicemelt(i) &
                 / (zse(k) * 1000.0_r_2) )
            xx(i) = soil%css(i) * soil%rhosoil(i)
            ssnow%gammzz(i,k) = MAX( &
                 REAL((1.0-soil%ssat(i)) * soil%css(i) * soil%rhosoil(i),r_2) &
                 + (ssnow%wb(i,k) - ssnow%wbice(i,k)) * REAL(cswat*rhowat,r_2)   &
                 + ssnow%wbice(i,k) * REAL(csice * rhowat * 0.9,r_2), &
                 REAL(xx(i),r_2) ) * zse(k)
            if (k == 1 .AND. ssnow%isflag(i) == 0) ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + cgsnow * ssnow%snowd(i)
            if (sicemelt(i)/ssnow%gammzz(i,k) > real(tiny(1.0),r_2)) &
                 ssnow%tgg(i,k) = ssnow%tgg(i,k) - REAL(sicemelt(i) / ssnow%gammzz(i,k)) * C%HLF
         endif
      enddo
   END DO

END SUBROUTINE soilfreeze_serial

! -----------------------------------------------------------------------------

SUBROUTINE remove_trans(dels, soil, ssnow, canopy, veg)

   USE cable_common_module, ONLY : cable_user

   ! Removes transpiration water from soil.
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(canopy_type), INTENT(INOUT)         :: canopy
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
   REAL(r_2), DIMENSION(mp,0:ms) :: diff
   REAL(r_2), DIMENSION(mp)      :: xx, xxd
   INTEGER :: k

  IF (cable_user%FWSOIL_switch.ne.'Haverd2013') THEN
     xx  = 0.0_r_2
     xxd = 0.0_r_2
     diff(:,:) = 0.0_r_2
     DO k = 1, ms

        ! Removing transpiration from soil:
        WHERE (canopy%fevc > 0.0_r_2)     ! convert to mm/dels
           ! Calculate the amount (perhaps moisture/ice limited)
           ! which can be removed:
           xx = canopy%fevc * real(dels / C%HL * veg%froot(:,k),r_2) + diff(:,k-1)   ! kg/m2
           diff(:,k) = MAX( 0.0_r_2, ssnow%wb(:,k) - real(soil%swilt,r_2)) &      ! m3/m3
                * real(soil%zse(k),r_2)*1000.0_r_2
           xxd = xx - diff(:,k)
           WHERE ( xxd .GT. 0.0_r_2 )
              ssnow%wb(:,k) = ssnow%wb(:,k) - diff(:,k) / real(soil%zse(k)*1000.0,r_2)
              diff(:,k) = xxd
           ELSEWHERE
              ssnow%wb(:,k) = ssnow%wb(:,k) - xx / real(soil%zse(k)*1000.0,r_2)
              diff(:,k) = 0.0_r_2
           ENDWHERE

        END WHERE

     END DO

  ELSE
     WHERE (canopy%fevc .lt. 0.0_r_2)
        canopy%fevw = canopy%fevw + real(canopy%fevc)
        canopy%fevc = 0.0_r_2
     END WHERE
     DO k = 1,ms
        ssnow%wb(:,k) = ssnow%wb(:,k) - real(ssnow%evapfbl(:,k)/(soil%zse(k)*1000.0),r_2)
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
SUBROUTINE soil_snow(dels, soil, ssnow, canopy, met, veg)
   USE cable_common_module
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(canopy_type), INTENT(INOUT)         :: canopy
   TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
   TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
   INTEGER             :: k
   REAL, DIMENSION(mp) :: snowmlt
   REAL, DIMENSION(mp) :: totwet
   REAL, DIMENSION(mp) :: weting
   REAL, DIMENSION(mp) :: xxx
   REAL(r_2), DIMENSION(mp) :: xx
   REAL, DIMENSION(mp) :: sinfil1, sinfil2, sinfil3
   REAL                :: zsetot
   INTEGER, SAVE :: ktau =0

   CALL point2constants( C )
   cp = C%CAPP

   ktau = ktau +1

   !jhan - make switchable
   ! appropriate for ACCESS1.0
   !max_glacier_snowd = 50000.0
   ! appropriate for ACCESS1.3
   max_glacier_snowd = 1100.0

   zsetot = sum(soil%zse)
   ssnow%tggav = 0.
   DO k = 1, ms
      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
   END DO

   IF( cable_runtime%offline .or. cable_runtime%mk3l ) THEN
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


   IF( .NOT.cable_user%cable_runtime_coupled ) THEN

      IF( ktau_gl <= 1 ) THEN

         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
                                                  ! after discussion with BP
         ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
         ssnow%wbtot = 0.0
         DO k = 1, ms
            ssnow%wb(:,k)  = MIN( soil%ssat, MAX( real(ssnow%wb(:,k)), soil%swilt ) )
         END DO

         ssnow%wb(:,ms-2)  = MIN( soil%ssat, MAX( real(ssnow%wb(:,ms-2)),           &
                             0.5 * ( soil%sfc + soil%swilt ) ) )
         ssnow%wb(:,ms-1)  = MIN( soil%ssat, MAX( real(ssnow%wb(:,ms-1)),           &
                             0.8 * soil%sfc ) )
         ssnow%wb(:,ms)    = MIN( soil%ssat, MAX( real(ssnow%wb(:,ms)), soil%sfc ) )

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

         xx = soil%css * soil%rhosoil

         ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
              & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * rhowat &
              & + ssnow%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1)
      END IF
   ENDIF  ! if(.NOT.cable_runtime_coupled)

   xx=soil%css * soil%rhosoil
   IF (ktau <= 1)                                                              &
     ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil      &
            & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * rhowat           &
            & + ssnow%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * cgsnow * ssnow%snowd

   DO k = 1, ms ! for stempv

      ! Set liquid soil water fraction (fraction of saturation value):
      ssnow%wblf(:,k) = MAX( 0.01_r_2, (ssnow%wb(:,k) - ssnow%wbice(:,k)) )    &
           & / REAL(soil%ssat,r_2)

      ! Set ice soil water fraction (fraction of saturation value):
      ssnow%wbfice(:,k) = ssnow%wbice(:,k) / real(soil%ssat, r_2)
   END DO

   CALL snowcheck(ssnow, soil)
   
   CALL snowdensity(dels, ssnow, soil)

   CALL snow_accum(dels, canopy, met, ssnow, soil)

   CALL snow_melting(dels, snowmlt, ssnow, soil)

   ! Add snow melt to global snow melt variable:
   ssnow%smelt = snowmlt

   ! Adjust levels in the snowpack due to snow accumulation/melting,
   ! snow aging etc...
   CALL snowl_adjust(ssnow)

   CALL stempv(dels, canopy, ssnow, soil)

   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   CALL snow_melting (dels, snowmlt, ssnow, soil )

   ! Add new snow melt to global snow melt variable:
   ssnow%smelt = ssnow%smelt + snowmlt

   CALL remove_trans(dels, soil, ssnow, canopy, veg)

   CALL soilfreeze_serial(soil, ssnow)

   totwet = canopy%precis + ssnow%smelt

   ! total available liquid including puddle
   weting = totwet + max(0.0, ssnow%pudsto - real(canopy%fesp)/C%HL*dels)
   xxx    = soil%ssat - real(ssnow%wb(:,1))
   sinfil1 = MIN(0.95*xxx*soil%zse(1)*rhowat, weting) !soil capacity
   xxx     = soil%ssat - real(ssnow%wb(:,2))
   sinfil2 = MIN(0.95*xxx*soil%zse(2)*rhowat, weting - sinfil1) !soil capacity
   xxx     = soil%ssat - real(ssnow%wb(:,3))
   sinfil3 = MIN(0.95*xxx*soil%zse(3)*rhowat, weting - sinfil1 - sinfil2)

   ! net water flux to the soil
   ssnow%fwtop1 = sinfil1 / dels - canopy%segg
   ssnow%fwtop2 = sinfil2 / dels
   ssnow%fwtop3 = sinfil3 / dels

   ! Puddle for the next time step
   ssnow%pudsto = max(0., weting - sinfil1 - sinfil2 - sinfil3)
   ssnow%rnof1  = max(0., ssnow%pudsto - ssnow%pudsmx)
   ssnow%pudsto = ssnow%pudsto - ssnow%rnof1

   CALL surfbv(dels, ssnow, soil, veg)

   ! correction required for energy balance in online simulations
   IF( cable_runtime%um ) THEN
      canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
      canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%dfe_ddq * ssnow%ddq_dtg)

      canopy%fhs = canopy%fhs+canopy%fhs_cor
      canopy%fes = canopy%fes+canopy%fes_cor
   ENDIF

   ! redistrb (set in cable.nml) by default==.FALSE.
   IF( redistrb )                                                              &
      CALL hydraulic_redistribution( dels, soil, ssnow, canopy, veg, met )

   ssnow%smelt = ssnow%smelt/dels

   ! Set weighted soil/snow surface temperature
   ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

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

   zsetot = sum(soil%zse)
   totalmoist(:) = 0.0
   totalice(:)   = 0.0
   DO k=1, ms
     totalmoist(:) = totalmoist(:) + real(ssnow%wb(:,k))    * soil%zse(k) / zsetot
     totalice(:)   = totalice(:)   + real(ssnow%wbice(:,k)) * soil%zse(k) / zsetot
   ENDDO

   Dtran=0.0
   WHERE( canopy%fevc < 10.0 .and.  totalice  < 1.e-2 )  Dtran=1.0

   DO k=1, ms
      S_VG(:,k) = MIN( 1.0, MAX( 1.0E-4, real(ssnow%wb(:,k)) - soil%swilt )          &
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
         frootX= max(0.01,max( veg%froot(:,k),veg%froot(:,j)))
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

            available(:)   = MAX( 0.0, real(ssnow%wb(:,k)) -                   &
                            ( soil%swilt(:) + ( soil%sfc(:) - soil%swilt(:) )  &
                             / 3. ) )
            accommodate(:) = MAX( 0.0, soil%ssat(:) - real(ssnow%wb(:,j)) )

            temp(:) = MAX( hr_perTime(:,k,j),                                  &
                          -1.0 * wiltParam * available(:),                     &
                          -1.0 * satuParam * accommodate(:) * soil%zse(j) /    &
                          soil%zse(k) )

            hr_perTime(:,k,j) = temp(:)
            hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)

         ELSEWHERE (hr_perTime(:,j,k) < 0.0)

           available(:)   = MAX( 0.0, real(ssnow%wb(:,j)) -                    &
                            ( soil%swilt(:) + ( soil%sfc(:) - soil%swilt(:) )  &
                            / 3. ) )

           accommodate(:) = MAX( 0.0, soil%ssat(:) - real(ssnow%wb(:,k)) )

           temp(:) = MAX( hr_perTime(:,j,k),                                   &
                         - 1.0 * wiltParam * available(:),                     &
                         -1.0 * satuParam * accommodate(:) * soil%zse(k) /     &
                         soil%zse(j) )

           hr_perTime(:,j,k) = temp(:)
           hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)

         ENDWHERE

         ssnow%wb(:,k) = ssnow%wb(:,k) + real(hr_perTime(:,k,j),r_2)
         ssnow%wb(:,j) = ssnow%wb(:,j) + real(hr_perTime(:,j,k),r_2)

      ENDDO

   ENDDO

   WHERE( met%tk < C%TFRZ + 5.  ) Dtran=0.0

   DO k=1, ms
      S_VG(:,k) = MIN( 1.0, MAX( 1.0E-4, real(ssnow%wb(:,k)) - soil%swilt )          &
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
         frootX= max(0.01,max( veg%froot(:,k),veg%froot(:,j)))
         hr_term(:,k,j) = CRT*(wpsy(:,j)-wpsy(:,k))*MAX(C_hr(:,k),C_hr(:,j)) &
            *(max(0.01,veg%froot(:,k))*max(0.01,veg%froot(:,j)))/(1-frootX)*Dtran
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

            available(:)   = MAX( 0.0, real(ssnow%wb(:,k)) - soil%sfc(:) )
            accommodate(:) = MAX( 0.0, soil%ssat(:) - real(ssnow%wb(:,j)) )

            temp(:) = MAX(hr_perTime(:,k,j),                                   &
                          -1.0 * wiltParam*available(:),                       &
                          -1.0 * satuParam * accommodate(:) * soil%zse(j) /    &
                          soil%zse(k) )

            hr_perTime(:,k,j) = temp(:)
            hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)

         ELSEWHERE (hr_perTime(:,j,k) < 0.0)

            available(:)   = MAX( 0.0, real(ssnow%wb(:,j)) - soil%sfc(:) )
            accommodate(:) = MAX( 0.0, soil%ssat(:) - real(ssnow%wb(:,k)) )

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

END MODULE cable_soil_snow_module
