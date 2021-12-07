!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from
! http://www.accessimulator.org.au/cable
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
!          Mark Decker - used ssnow as base for ssgw.  Could be part of same module
!          Aug 2017 - applied changes for cls and rev_corr packages
!                   - partner changes applied in soilsnow module
!
! ==============================================================================

MODULE cable_gw_hydro_module

  USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
       veg_parameter_type, canopy_type, met_type,        &
       balances_type, r_2, ms, mp

  USE cable_common_module, ONLY : gw_params,cable_user,&
       cable_runtime,&
       max_glacier_snowd

!distribute these per sbr
USE cable_phys_constants_mod, ONLY : CTFRZ => TFRZ
USE cable_phys_constants_mod, ONLY : CHL => HL
USE cable_phys_constants_mod, ONLY : CHLF => HLF
USE cable_phys_constants_mod, ONLY : CHLS => HLS
USE cable_phys_constants_mod, ONLY : Cdensity_liq => density_liq
USE cable_phys_constants_mod, ONLY : Cdensity_ice => density_ice
USE cable_phys_constants_mod, ONLY : Ccgsnow => cgsnow
USE cable_phys_constants_mod, ONLY : Ccswat => cswat
USE cable_phys_constants_mod, ONLY : Ccs_rho_wat => cs_rho_wat
USE cable_phys_constants_mod, ONLY : Ccs_rho_ice => cs_rho_ice
USE cable_math_constants_mod, ONLY : CPI => PI

  IMPLICIT NONE

  PRIVATE

  !mrd561 GW params
  REAL(r_2), SAVE :: smp_cor = 8.0
  REAL(r_2), PARAMETER :: sucmin       = -1.0e8      ! minimum soil pressure head [mm]
  REAL(r_2), PARAMETER :: volwatmin    = 1e-4        !min soil water [mm]
  REAL(r_2), PARAMETER :: wtd_uncert   = 0.1         ! uncertaintiy in wtd calcultations [mm]
  REAL(r_2), PARAMETER :: wtd_max      = 1000000.0   ! maximum wtd [mm]
  REAL(r_2), PARAMETER :: wtd_min      = 10.0        ! minimum wtd [mm]
  REAL(r_2), PARAMETER :: close_to_one = 0.9999
  REAL(r_2), PARAMETER :: Sy_deep = 0.1
  REAL(r_2), PARAMETER :: dz_deep = 50000.0

  INTEGER, PARAMETER :: wtd_iter_max = 20 ! maximum number of iterations to find the water table depth

  ! ! This module contains the following subroutines that
  !are called from other modules
  PUBLIC :: soil_snow_gw,calc_srf_wet_fraction,sli_hydrology,&
       pore_space_relative_humidity,set_unsed_gw_vars

CONTAINS


  ! -----------------------------------------------------------------------------

  SUBROUTINE GWsoilfreeze(dels, soil, ssnow,tgg_old)
    !NOTE: this is only included because gw_model uses parameters XXX_vec
    !these are r_2.  this breaks bitwise compatibility with trunk
    !if acceptable this routine does the same thing but with r_2 soil params
    ! if max_ice_frac always set to frozen_limit and tgg_tmp is always CTFRZ
IMPLICIT NONE

    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil
    REAL, INTENT(INOUT), DIMENSION(mp,ms) :: tgg_old
    REAL(r_2), DIMENSION(mp)           :: sicefreeze
    REAL(r_2), DIMENSION(mp)           :: sicemelt
    REAL(r_2), DIMENSION(mp,ms)        :: wbice_delta,avail_por
    REAL(r_2), DIMENSION(mp)           :: ice_mass,liq_mass,tot_mass
    INTEGER :: i,j,k
    REAL, DIMENSION(mp,ms) :: tgg_tmp
    REAL(r_2),DIMENSION(mp,ms) :: xx,max_ice_frac,iceF,den_css  !Decker and Zeng 2009

    max_ice_frac(:,:) = 0.0
    DO k=1,ms
       DO i=1,mp
          IF  (ssnow%tgg(i,k) .LT. CTFRZ .AND. soil%ssat_vec(i,k) .GT. 1.0e-8) THEN
             max_ice_frac(i,k) = (1. - EXP(-2.*(ssnow%wb(i,k)/soil%ssat_vec(i,k))**4.0 *&
                  (ssnow%tgg(i,k)-CTFRZ)))/EXP(1. - ssnow%wb(i,k)/soil%ssat_vec(i,k))
             max_ice_frac(i,k) = MAX(0.4,max_ice_frac(i,k))*ssnow%wb(i,k)

             wbice_delta(i,k) = MAX(0.,max_ice_frac(i,k) - ssnow%wbice(i,k))

             avail_por(i,k) = soil%ssat_vec(i,k) - ssnow%wbliq(i,k)+&
                  Cdensity_ice/Cdensity_liq*wbice_delta(i,k) - &
                  ssnow%wbice(i,k)

             wbice_delta(i,k) = MIN(wbice_delta(i,k),avail_por(i,k))

             max_ice_frac(i,k) = ssnow%wbice(i,k) + wbice_delta(i,k)

          ENDIF
       END DO
    END DO


    tgg_tmp(:,:) = tgg_old(:,:)
    DO k=1,ms
       DO i=1,mp
          IF (soil%isoilm(i) .EQ. 9) THEN
             tgg_tmp(i,k) = CTFRZ
          ELSE
             IF  (ssnow%tgg(i,k) .LE. CTFRZ) THEN
                IF (tgg_old(i,k) .GT. CTFRZ) THEN
                   tgg_tmp(i,k) = CTFRZ
                END IF
             ELSE
                IF (tgg_old(i,k) .LE. CTFRZ) THEN
                   tgg_tmp(i,k) = CTFRZ
                END IF
             END IF
          END IF
       END DO
    END DO



    !allow more freezing for permenant glacier ice regions
    DO i=1,mp
       IF (soil%isoilm(i) .EQ. 9) max_ice_frac(i,:) = 0.85_r_2*ssnow%wb(i,:)
    END DO

    DO k = 1, ms
       DO i=1,mp

          ice_mass(i) = ssnow%wbice(i,k)*soil%zse_vec(i,k)*REAL(Cdensity_ice,r_2)
          liq_mass(i) = ssnow%wbliq(i,k)*soil%zse_vec(i,k)*REAL(Cdensity_liq,r_2)
          tot_mass(i) = liq_mass(i) + ice_mass(i)

          IF (ssnow%tgg(i,k) .LE. CTFRZ .AND. &
               ssnow%tgg(i,k) .LT. tgg_tmp(i,k) .AND. &
               max_ice_frac(i,k) - ssnow%wbice(i,k) > .001) THEN

             sicefreeze(i) = MIN( MAX( 0.0_r_2, ( max_ice_frac(i,k)  -      &
                  ssnow%wbice(i,k) ) ) * soil%zse_vec(i,k) * Cdensity_ice,             &
                  ( tgg_tmp(i,k) - ssnow%tgg(i,k) ) * ssnow%gammzz(i,k) / Chlf )

             ssnow%wbice(i,k) = MIN( ssnow%wbice(i,k) +&
                  sicefreeze(i)/soil%zse_vec(i,k)/Cdensity_ice,&
                  max_ice_frac(i,k) )
             ssnow%gammzz(i,k) = MAX(soil%heat_cap_lower_limit(i,k),  &
                  (1.0-soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k) &
                  + (ssnow%wb(i,k)-ssnow%wbice(i,k)) * REAL(Ccs_rho_wat,r_2)   &
                  + ssnow%wbice(i,k) * REAL(Ccs_rho_ice,r_2)&
                  )*soil%zse_vec(i,k)

             IF (k .EQ. 1 .AND. ssnow%isflag(i) .EQ. 0) THEN
                ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + Ccgsnow * ssnow%snowd(i)
             END IF

             ssnow%tgg(i,k) = ssnow%tgg(i,k) + REAL(sicefreeze(i))                    &
                  * Chlf / REAL(ssnow%gammzz(i,k) )

          ELSEIF ( ssnow%tgg(i,k) .GT. tgg_tmp(i,k) .AND. &
               ssnow%wbice(i,k) .GT.  max_ice_frac(i,k) ) THEN

             sicemelt(i) = MIN( ssnow%wbice(i,k) * soil%zse_vec(i,k) * Cdensity_ice,              &
                  ( ssnow%tgg(i,k) - tgg_tmp(i,k) ) * ssnow%gammzz(i,k) / Chlf )

             ssnow%wbice(i,k) = MAX( 0.0_r_2, ssnow%wbice(i,k) - sicemelt(i)          &
                  / (soil%zse_vec(i,k) * Cdensity_ice) )

             ssnow%gammzz(i,k) =MAX(soil%heat_cap_lower_limit(i,k), &
                  (1.0-soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k)&
                  + (ssnow%wb(i,k)-ssnow%wbice(i,k)) * REAL(Ccs_rho_wat,r_2)   &
                  + ssnow%wbice(i,k) * REAL(Ccs_rho_ice,r_2)&
                  )*soil%zse_vec(i,k)

             IF (k .EQ. 1 .AND. ssnow%isflag(i) .EQ. 0) THEN
                ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + Ccgsnow * ssnow%snowd(i)
             END IF
             ssnow%tgg(i,k) = ssnow%tgg(i,k) - REAL(sicemelt(i))                     &
                  * Chlf / REAL(ssnow%gammzz(i,k))
             !if for some reason end up here
          ELSEIF( tgg_tmp(i,k) .GE. CTFRZ .AND. &
               ssnow%tgg(i,k) > tgg_tmp(i,k) .AND. ssnow%wbice(i,k) > 0. ) THEN

             sicemelt(i) = MIN( ssnow%wbice(i,k) * soil%zse_vec(i,k) * Cdensity_ice,              &
                  ( ssnow%tgg(i,k) - tgg_tmp(i,k) ) * ssnow%gammzz(i,k) / Chlf )

             ssnow%wbice(i,k) = MAX( 0.0_r_2, ssnow%wbice(i,k) - sicemelt(i)          &
                  / (soil%zse_vec(i,k) * Cdensity_ice) )

             ssnow%gammzz(i,k) =MAX(soil%heat_cap_lower_limit(i,k), &
                  (1.0-soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k)&
                  + (ssnow%wb(i,k)-ssnow%wbice(i,k)) * REAL(Ccs_rho_wat,r_2)   &
                  + ssnow%wbice(i,k) * REAL(Ccs_rho_ice,r_2)&
                  )*soil%zse_vec(i,k)

             IF (k .EQ. 1 .AND. ssnow%isflag(i) .EQ. 0) THEN
                ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + Ccgsnow * ssnow%snowd(i)
             END IF
             ssnow%tgg(i,k) = ssnow%tgg(i,k) - REAL(sicemelt(i))                     &
                  * Chlf / REAL(ssnow%gammzz(i,k))

          END IF
          !update the liq and ice volume and mass
          ice_mass(i)   = ssnow%wbice(i,k)*soil%zse_vec(i,k)*REAL(Cdensity_ice,r_2)
          liq_mass(i)   = tot_mass(i) - ice_mass(i)
          ssnow%wbliq(i,k) = liq_mass(i) / soil%zse_vec(i,k)/REAL(Cdensity_liq,r_2)
          ssnow%wbice(i,k) = ice_mass(i) / soil%zse_vec(i,k)/REAL(Cdensity_ice,r_2)
          ssnow%wb(i,k)    = ssnow%wbliq(i,k) + ssnow%wbice(i,k)

       END DO
    END DO

  END SUBROUTINE GWsoilfreeze
  ! -----------------------------------------------------------------------------
  !
  !! -----------------------------------------------------------------------------
  !
  SUBROUTINE remove_transGW(dels, soil, ssnow, canopy, veg)
    !NOTE: this is only included because gw_model uses parameters XXX_vec
    !these are r_2.  this breaks bitwise compatibility with trunk
    !if acceptable this routine does the same thing but with r_2 soil params

    ! Removes transpiration water from soil.
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    REAL(r_2), DIMENSION(mp,0:ms+1) :: diff
    REAL(r_2), DIMENSION(mp)      :: xx,xxd,evap_cur
    REAL(r_2), DIMENSION(mp,ms) :: zse_mp_mm
    INTEGER :: k,i

    DO k=1,ms
       DO i=1,mp
          zse_mp_mm(i,k)  = REAL(soil%zse_vec(i,k)*Cdensity_liq,r_2)
       END DO
    END DO

    IF (cable_user%FWSOIL_switch.NE.'Haverd2013') THEN

       xx(:) = 0._r_2
       xxd(:) = 0._r_2
       diff(:,:) = 0._r_2

       DO k = 1,ms

          DO i=1,mp

             IF (canopy%fevc(i) .GT. 0._r_2) THEN

                xx(i) = canopy%fevc(i) * dels / Chl * veg%froot(i,k) + diff(i,k-1)
                diff(i,k) = MAX(0._r_2,ssnow%wbliq(i,k)-soil%swilt_vec(i,k)) &
                     * zse_mp_mm(i,k)
                xxd(i) = xx(i) - diff(i,k)

                IF (xxd(i) .GT. 0._r_2) THEN
                   ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - diff(i,k)/zse_mp_mm(i,k)
                   diff(i,k) = xxd(i)
                ELSE
                   ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - xx(i)/zse_mp_mm(i,k)
                   diff(i,k) = 0._r_2
                END IF


             END IF  !fvec > 0

          END DO  !mp
       END DO     !ms

    ELSE

       WHERE (canopy%fevc .LT. 0.0_r_2)
          canopy%fevw = canopy%fevw+canopy%fevc
          canopy%fevc = 0.0_r_2
       END WHERE
       DO k = 1,ms
          ssnow%wbliq(:,k) = ssnow%wbliq(:,k) - ssnow%evapfbl(:,k)/zse_mp_mm(:,k)
       ENDDO

    ENDIF

    DO k=1,ms
       DO i=1,mp
          ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*zse_mp_mm(i,k)!mass
          ssnow%wmtot(i,k) = ssnow%wmliq(i,k) + ssnow%wmice(i,k)  !mass
          ssnow%wb(i,k)    = ssnow%wbliq(i,k) + ssnow%wbice(i,k)  !volume
       END DO
    END DO


  END SUBROUTINE remove_transGW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!MD GW code from here on!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-------------------------------------------------------------------------
  SUBROUTINE ovrlndflx (dels, ssnow, soil,veg, canopy,sli_call )
    USE cable_common_module, ONLY : gw_params,cable_user

    IMPLICIT NONE
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil ! soil parameters
    TYPE(veg_parameter_type) , INTENT(INOUT)    :: veg  ! veg parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy
    LOGICAL, INTENT(IN)                      :: sli_call
    INTEGER                                  :: nglacier ! 0 original, 1 off, 2 new Eva
    INTEGER                                  :: k, i, j
    REAL, DIMENSION(mp)                :: rnof5
    REAL, DIMENSION(mp)                :: sgamm
    REAL, DIMENSION(mp)                :: smasstot
    REAL, DIMENSION(mp,0:3)            :: smelt1                   !snow melt
    REAL(r_2), DIMENSION(mp)           :: icef,efpor               !tmp vars, fraction of ice in gridcell
    REAL(r_2)                          :: tmpa,tmpb,qinmax         !tmp vars, maximum infiltration [mm/s]
    REAL(r_2), DIMENSION(mp)           :: satfrac_liqice,S       !saturated fraction of cell, wtd in m
    REAL(r_2)                          :: liqmass,icemass,totmass  !liquid mass,ice mass, total mass [mm]
    REAL(r_2)                          :: fice
    REAL(r_2)                          :: dzmm,slopeSTDmm

    !For now assume there is no puddle
    dzmm = 1000._r_2 * soil%zse(1)

    IF (sli_call) THEN
       DO i=1,mp
          IF (canopy%through(i) .GE. canopy%through_sn(i)) THEN
             ssnow%fwtop(i)  = MAX((canopy%through(i)-canopy%through_sn(i))/dels , 0.)             ! liq precip rate (m s-1)
          ELSE
             ssnow%fwtop(i) = MAX(canopy%through(i), 0.)
          END IF
       END DO
    END IF
    !amount of ice in surface layer
    DO i = 1,mp
       efpor(i) = MAX(0.001_r_2, soil%ssat_vec(i,1) - ssnow%wbice(i,1))
       icemass  = ssnow%wbice(i,1) * dzmm
       liqmass  = (ssnow%wb(i,1)-ssnow%wbice(i,1)) * dzmm
       totmass  = MAX(liqmass+icemass,REAL(1e-2,r_2))
       icef(i)     = MAX(0._r_2,MIN(1._r_2,gw_params%IceBeta*icemass / totmass))
    END DO

    !sat fraction assuming topo controlled subgrid soil moisture distribution
    !called from cable_canopy for srf wet fraction alrady
    !call saturated_fraction(ssnow,soil,veg)

    !srf frozen fraction.  should be based on topography
    DO i = 1,mp
       fice = (EXP(-gw_params%IceAlpha*(1._r_2-icef(i)))-EXP(-gw_params%IceAlpha))/(1._r_2-EXP(-gw_params%IceAlpha))
       fice  = MIN(MAX(fice,0._r_2),1._r_2)
       satfrac_liqice(i)   = MAX(0.,MIN(0.95,fice + (1._r_2-fice)*ssnow%satfrac(i) ) )
    END DO

    DO i=1,mp
       tmpa = ssnow%wbliq(i,1) / efpor(i)
       tmpb = MAX( (tmpa-satfrac_liqice(i))/MAX(0.01_r_2,(1._r_2-satfrac_liqice(i))), 0._r_2)
       tmpa = -2._r_2*soil%bch_vec(i,1)*soil%sucs_vec(i,1)/dzmm
       qinmax = (1._r_2 + tmpa*(tmpb-1._r_2))*soil%hyds_vec(i,1)*EXP(-gw_params%hkrz*(0.5*dzmm/1000.0_r_2-gw_params%zdepth))

       ssnow%rnof1(i) = satfrac_liqice(i) * ssnow%fwtop(i) + &
            (1._r_2-satfrac_liqice(i))*MAX((ssnow%fwtop(i)-qinmax) , 0._r_2)

       ssnow%fwtop(i) = ssnow%fwtop(i) - ssnow%rnof1(i)

    END DO  !mp

    !add back to the lakes to keep saturated instead of drying
    DO i=1,mp
       IF (veg%iveg(i) .EQ. 16) THEN
          ssnow%fwtop(i) = ssnow%fwtop(i) + ssnow%rnof1(i)
          ssnow%rnof1(i) = 0._r_2
       END IF
    END DO

    !---  glacier formation
    rnof5= 0.

    IF (sli_call .OR. cable_runtime%UM) THEN
       nglacier = 0
    ELSE
       nglacier = 2
    END IF

    IF (nglacier == 2) THEN
       smelt1=0.
       WHERE( ssnow%snowd > max_glacier_snowd )

          rnof5 = MIN( 0.1, ssnow%snowd - max_glacier_snowd )

          !---- change local tg to account for energy - clearly not best method
          WHERE( ssnow%isflag == 0 )
             smasstot = 0.0
             ssnow%tgg(:,1) = ssnow%tgg(:,1) - rnof5 * Chlf                    &
                  / REAL( ssnow%gammzz(:,1) )
             ssnow%snowd = ssnow%snowd - rnof5
          ELSEWHERE
             smasstot = ssnow%smass(:,1) + ssnow%smass(:,2) + ssnow%smass(:,3)
          END WHERE

       END WHERE

       DO k = 1, 3

          WHERE( ssnow%snowd > max_glacier_snowd  .AND.  ssnow%isflag > 0 )
             sgamm = ssnow%ssdn(:,k) * Ccgsnow * ssnow%sdepth(:,k)
             smelt1(:,k) = MIN( rnof5 * ssnow%smass(:,k) / smasstot,            &
                  0.2 * ssnow%smass(:,k) )
             ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)
             ssnow%snowd = ssnow%snowd - smelt1(:,k)
          END WHERE

       END DO

       WHERE( ssnow%isflag > 0 ) rnof5 = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)

       ssnow%rnof1 = ssnow%rnof1 + rnof5/dels   !include this runoff in suface runoff term

    END IF

  END SUBROUTINE ovrlndflx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE simple_wtd(ssnow, soil, veg)
    !This was only for testing purposes
    IMPLICIT NONE
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg

    REAL(r_2), DIMENSION(mp)            :: fz, wmean,ztot
    REAL(r_2), DIMENSION(mp,ms)         :: stot
    INTEGER                             :: k,i

    DO i=1,mp
       wmean(i) = 0._r_2
       fz(i)    = 5._r_2
       ztot(i)  = 0._r_2
       stot(i,:) = (ssnow%wb(i,:)-soil%watr(i,:)) / (soil%ssat_vec(i,:)-soil%watr(i,:))
    END DO
    DO k  = 1, ms
       DO i=1,mp
          wmean(i) = wmean(i) + stot(i,k)*soil%zse(k)*1000._r_2
          ztot(i)  = ztot(i) + soil%zse(k)*1000._r_2
       END DO
    END DO

    DO i=1,mp
       wmean(i) = wmean(i) + ssnow%GWwb(i)/soil%GWssat_vec(i) * soil%GWdz(i)*1000._r_2
       ztot(i)  = ztot(i) + soil%GWdz(i)*1000._r_2

       ssnow%wtd(i) = MIN(200000._r_2, fz(i) * (ztot(i) - wmean(i)))
    END DO

  END SUBROUTINE simple_wtd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !----------------------------------------------------------------------
  ! SUBROUTINE iterative_wtd
  !
  ! Iteratively calcs the water table depth by equating the mass of water in the
  ! soil column to the mass of a hydrostatic column inegrated from the surface to the
  ! water table depth
  !
  SUBROUTINE iterative_wtd (ssnow, soil, veg, include_aquifer)
    IMPLICIT NONE
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    LOGICAL, INTENT(IN)                       :: include_aquifer  !use GWwb or only wb to find wtd?


    !Local vars
    REAL(r_2), DIMENSION(mp,ms)   :: dzmm_mp,tmp_def
    REAL(r_2), DIMENSION(0:ms)    :: zimm
    REAL(r_2), DIMENSION(ms)      :: zmm
    REAL(r_2), DIMENSION(mp)      :: GWzimm,temp
    REAL(r_2), DIMENSION(mp)      :: def,defc,total_depth_column

    REAL(r_2)                     :: deffunc,tempa,tempb,derv,calc,tmpc
    REAL(r_2), DIMENSION(mp)      :: invB,Nsucs_vec  !inverse of C&H B,Nsucs_vec
    INTEGER :: k,i,wttd,jlp

    !make code cleaner define these here
    invB     = 1._r_2/soil%bch_vec(:,ms)                                !1 over C&H B
    Nsucs_vec  = soil%sucs_vec(:,ms)                                !psi_saturated mm
    dzmm_mp  = REAL(SPREAD((soil%zse(:)) * 1000.0,1,mp),r_2)    !layer thickness mm
    zimm(0)  = 0.0_r_2                                          !depth of layer interfaces mm

    !total depth of soil column
    DO k=1,ms
       zimm(k) = zimm(k-1) + soil%zse(k)*1000._r_2
    END DO

    def(:) = 0._r_2

    IF (include_aquifer) THEN  !do we include the aquifer in the calculation of wtd?

       DO i=1,mp
          total_depth_column(i) = zimm(ms) + soil%GWdz(i)*1000._r_2
          def(i) = def(i) + MAX(0._r_2,soil%GWssat_vec(i)-ssnow%GWwb(i))*soil%GWdz(i)*1000._r_2
       END DO

    END IF

    !comute the total mass away from full saturation
    DO k=1,ms
       DO i=1,mp

          def(i) = def(i) +                                                           &
               MAX(0._r_2,(soil%ssat_vec(i,k)-(ssnow%wbliq(i,k)+ssnow%wbice(i,k)))*dzmm_mp(i,k))
       END DO  !mp
    END DO  !ms

    !find the deficit if the water table is at the bottom of the soil column
    DO i=1,mp
       defc(i) = (soil%ssat_vec(i,ms))*(total_depth_column(i)+Nsucs_vec(i)/(1._r_2-invB(i))*            &
            (1._r_2-((Nsucs_vec(i)+total_depth_column(i))/Nsucs_vec(i))**(1._r_2-invB(i))))
       defc(i) = MAX(0.1_r_2,defc(i))

       !initial guess at wtd
       ssnow%wtd(:) = total_depth_column(:)*def(:)/defc(:)
    END DO


    !use newtons method to solve for wtd, note this assumes homogenous column but
    !that is ok
    DO i=1,mp
       IF ((soil%isoilm(i) .NE. 9) .AND. (veg%iveg(i) .NE. 16)) THEN

          IF (defc(i) > def(i)) THEN                 !iterate tfor wtd

             jlp=0

             mainloop: DO

                tempa   = 1.0_r_2
                tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-invB(i))
                derv    = (soil%ssat_vec(i,ms))*(tempa-tempb) + &
                     soil%ssat_vec(i,ms)

                IF (ABS(derv) .LT. REAL(1e-8,r_2)) derv = SIGN(REAL(1e-8,r_2),derv)

                tempa   = 1.0_r_2
                tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-invB(i))
                deffunc = (soil%ssat_vec(i,ms))*(ssnow%wtd(i) +&
                     Nsucs_vec(i)/(1-invB(i))* &
                     (tempa-tempb)) - def(i)
                calc    = ssnow%wtd(i) - deffunc/derv

                IF ((ABS(calc-ssnow%wtd(i))) .LE. wtd_uncert) THEN

                   ssnow%wtd(i) = calc
                   EXIT mainloop

                ELSEIF (jlp .GE. wtd_iter_max) THEN

                   EXIT mainloop

                ELSE

                   jlp=jlp+1
                   ssnow%wtd(i) = calc

                END IF

             END DO mainloop  !defc .gt. def

          ELSEIF (defc(i) .LT. def(i)) THEN

             jlp=0

             mainloop2: DO

                tmpc     = Nsucs_vec(i)+ssnow%wtd(i)-total_depth_column(i)
                tempa    = (ABS(tmpc/Nsucs_vec(i)))**(-invB(i))
                tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-invB(i))
                derv     = (soil%ssat_vec(i,ms))*(tempa-tempb)
                IF (ABS(derv) .LT. REAL(1e-8,r_2)) derv = SIGN(REAL(1e-8,r_2),derv)

                tempa    = (ABS((Nsucs_vec(i)+ssnow%wtd(i)-total_depth_column(i))/Nsucs_vec(i)))**(1._r_2-invB(i))
                tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-invB(i))
                deffunc  = (soil%ssat_vec(i,ms))*(total_depth_column(i) +&
                     Nsucs_vec(i)/(1._r_2-invB(i))*(tempa-tempb))-def(i)
                calc     = ssnow%wtd(i) - deffunc/derv

                IF ((ABS(calc-ssnow%wtd(i))) .LE. wtd_uncert) THEN

                   ssnow%wtd(i) = calc
                   EXIT mainloop2

                ELSEIF (jlp==wtd_iter_max) THEN

                   EXIT mainloop2

                ELSE

                   jlp=jlp+1
                   ssnow%wtd(i) = calc

                END IF

             END DO mainloop2  !defc .lt. def

          ELSE  !water table depth is exactly on bottom boundary

             ssnow%wtd(i) = total_depth_column(i)

          ENDIF

       ENDIF  !check veg and soils

    END DO   !mp loop

    !limit wtd to be within a psecified range
    DO i=1,mp
       IF (veg%iveg(i) .GE. 16) ssnow%wtd(i) = wtd_min
       ssnow%wtd(i) = MIN(wtd_max,MAX(wtd_min,ssnow%wtd(i) ) )
    END DO


  END SUBROUTINE iterative_wtd

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoistgw (fwtop,dt,ktau,ssnow,soil,prin)
  ! solves the modified richards equation (Zeng and Decker 2009) to find
  ! vertical mocement of soil water.  Bottom boundary condition is determined
  ! using a single layer groundwater module
  !
  SUBROUTINE smoistgw (dels,ktau,ssnow,soil,veg,canopy)
    USE cable_common_module
USE trimb_mod,                       ONLY : trimb

    IMPLICIT NONE

    REAL, INTENT(IN)                          :: dels  ! time step size (s)
    INTEGER, INTENT(IN)                       :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    TYPE(canopy_type), INTENT(INOUT)          :: canopy ! vegetation variables

    !Local variables.
    REAL(r_2), DIMENSION(mp,ms+1)       :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: rt

    INTEGER                             :: k,kk,i
    REAL(r_2), DIMENSION(mp,ms)         :: eff_por,old_wb  !effective porosity (mm3/mm3),wb(mm3/mm3),mass (mm) of eff_por
    REAL(r_2), DIMENSION(mp,ms)         :: msliq,msice             !mass of the soil liquid and ice water
    REAL(r_2), DIMENSION(mp)            :: den
    REAL(r_2), DIMENSION(mp)            :: dne
    REAL(r_2), DIMENSION(mp)            :: num
    REAL(r_2), DIMENSION(mp)            :: qin
    REAL(r_2), DIMENSION(mp)            :: qout
    REAL(r_2), DIMENSION(mp)            :: dqidw0
    REAL(r_2), DIMENSION(mp)            :: dqidw1
    REAL(r_2), DIMENSION(mp)            :: dqodw0
    REAL(r_2), DIMENSION(mp)            :: dqodw1,dqodw2
    REAL(r_2), DIMENSION(mp)            :: s1,s2,tmpi,temp0,voleq1,tempi
    REAL(r_2), DIMENSION(ms)            :: dzmm
    REAL(r_2), DIMENSION(mp,ms)         :: dzmm_mp
    REAL(r_2), DIMENSION(0:ms+1)        :: zimm
    REAL(r_2), DIMENSION(ms)            :: zmm
    REAL(r_2), DIMENSION(mp)            :: GWzimm,xs,zaq,s_mid,GWdzmm
    REAL(r_2), DIMENSION(mp)            :: xs1,GWmsliq!xsi    !mass (mm) of liquid over/under saturation, mass of aquifer water
    REAL(r_2)                           :: xsi
    REAL(r_2), DIMENSION(mp,ms+1)       :: del_wb
    !MD DEBUG VARS
    INTEGER :: imp,ims,k_drain

    !make code cleaner define these here
    dzmm(:) = 1000.0_r_2 * REAL(soil%zse(:),r_2)
    DO i=1,mp
       dzmm_mp(i,:) = dzmm(:)
    END DO

    zimm(0) = 0._r_2
    DO k=1,ms
       zimm(k) = zimm(k-1) + dzmm(k)
    END DO
    zmm(1:ms)  = zimm(0:(ms-1)) + 0.5_r_2*dzmm(1:ms)

    DO i=1,mp
       GWdzmm(i) = REAL(soil%GWdz(i),r_2)*1000._r_2
       GWzimm(i) = zimm(ms)+GWdzmm(i)
       zaq(i)    = zimm(ms) + 0.5_r_2*GWdzmm(i)
    END DO

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! preset to allow for non-land & snow points in trimb
    DO k=1,ms
       DO i=1,mp
          old_wb(i,k) = ssnow%wb(i,k)
          rt(i,k) = 0._r_2
          at(i,k) = 0._r_2
          bt(i,k) = 0._r_2
          ct(i,k) = 0._r_2
       END DO
    END DO

    !equilibrium water content
    CALL calc_equilibrium_water_content(ssnow,soil)

    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    CALL subsurface_drainage(ssnow,soil,veg,dzmm)

    k = 1     !top soil layer
    DO i=1,mp
       qin(i)     = ssnow%sinfil(i)
       den(i)     = (zmm(k+1)-zmm(k))
       dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
       num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
       dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
       dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
       rt(i,k) =  qin(i) - qout(i)
       at(i,k) =  0._r_2
       bt(i,k) =  dzmm(k)/dels + dqodw1(i)
       ct(i,k) =  dqodw2(i)
    END DO
    DO k = 2, ms - 1     !middle soil layers
       DO i=1,mp
          den(i)     = (zmm(k) - zmm(k-1))
          dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
          num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
          qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
          dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          den(i)     = (zmm(k+1)-zmm(k))
          dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
          num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
          qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
          dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
          dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
          rt(i,k) =  qin(i) - qout(i)
          at(i,k) = -dqidw0(i)
          bt(i,k) =  dzmm(k)/dels - dqidw1(i) + dqodw1(i)
          ct(i,k) =  dqodw2(i)
       END DO
    END DO

    k = ms   !Bottom soil layer
    DO i=1,mp
       den(i)     = (zmm(k) - zmm(k-1))
       dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
       num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
       qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
       dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       den(i)     = zaq(i) - zmm(k)
       dne(i)     = (ssnow%GWzq(i)-ssnow%zq(i,k))
       num(i)     =  (ssnow%GWsmp(i)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = 0._r_2
       dqodw1(i)  = 0._r_2
       dqodw2(i)  = 0._r_2
       rt(i,k) =  qin(i) - qout(i)
       at(i,k) = -dqidw0(i)
       bt(i,k) =  dzmm(k)/dels - dqidw1(i) + dqodw1(i)
       ct(i,k) =  dqodw2(i)
    END DO

    CALL aquifer_recharge(dels,ssnow,soil,veg,zaq,zmm,dzmm)

    CALL trimb(at,bt,ct,rt,ms)                       !use the defulat cable tridiag solution

    DO k=1,ms
       DO i=1,mp
          ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + rt(i,k) - ssnow%qhlev(i,k)*dels/dzmm(k)   !volutermic liquid
       END DO
    END DO

    DO i=1,mp
       ssnow%wbliq(i,ms) = ssnow%wbliq(i,ms) - ssnow%Qrecharge(i)*dels/dzmm(ms)
    END DO
    DO i=1,mp
       ssnow%GWwb(i) = ssnow%GWwb(i)  +  (ssnow%Qrecharge(i)-ssnow%qhlev(i,ms+1))*dels/GWdzmm(i)
    END DO

    !determine the available pore space
    !volumetric
    DO k=1,ms
       DO i=1,mp
          eff_por(i,k)  = MAX(0._r_2, soil%ssat_vec(i,k) - ssnow%wbice(i,k) )
       END DO
    END DO

    DO i=1,mp
       xsi = 0._r_2

       IF (ssnow%GWwb(i) .GT. soil%GWssat_vec(i)) THEN
          xsi = (ssnow%GWwb(i) - soil%GWssat_vec(i))*GWdzmm(i)
          ssnow%GWwb(i) = soil%GWssat_vec(i)
       END IF

       DO k=1,ms
          IF (ssnow%wbliq(i,k) .GT. eff_por(i,k)) THEN
             xsi = xsi + (ssnow%wbliq(i,k) - eff_por(i,k))*dzmm(k)
             ssnow%wbliq(i,k) = eff_por(i,k)
          END IF
       END DO

       DO k = ms,1,-1  !loop from bottom to top adding extra water to each layer
          IF (xsi .GT. 0._r_2) THEN
             IF (xsi .LT. (eff_por(i,k)-ssnow%wbliq(i,k))*dzmm(k)) THEN
                ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + xsi/dzmm(k)
                xsi = 0._r_2
             ELSE
                xsi = xsi - (eff_por(i,k) - ssnow%wbliq(i,k))*dzmm(k)
                ssnow%wbliq(i,k) = eff_por(i,k)
             END IF
          END IF
       END DO  !ms loop

       IF (xsi .GT. 0._r_2) THEN
          ssnow%qhz(i) = ssnow%qhz(i) + xsi/dels
          xsi = 0._r_2
       END IF

       DO k = 1,ms
          xsi = 0._r_2             !should be a single float (array not needed)
          IF (ssnow%wbliq(i,k) .LT. volwatmin) THEN
             xsi = (volwatmin - ssnow%wbliq(i,k))*dzmm(k)  !in mm
             ssnow%wbliq(i,k) = volwatmin
             IF (k .LT. ms) THEN
                ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) - xsi/dzmm(k+1)
             ELSE
                ssnow%GWwb(i) = ssnow%GWwb(i) - xsi / GWdzmm(i)
             END IF
          END IF
       END DO  !ms loop

       IF ( (ssnow%GWwb(i) .LT. volwatmin) .AND. (soil%isoilm(i) .NE. 9) ) THEN
          xsi = (volwatmin - ssnow%GWwb(i)) / GWdzmm(i)  !mm
          ssnow%GWwb(i) = volwatmin
          ssnow%qhz(i) = ssnow%qhz(i) - xsi / dels
       END IF

    END DO

    DO k=1,ms
       DO i=1,mp

          !update mass variables
          ssnow%wmliq(i,k)      = ssnow%wbliq(i,k) * &
               soil%zse_vec(i,k)*Cdensity_liq
          ssnow%wmice(i,k)      = ssnow%wbice(i,k) * &
               soil%zse_vec(i,k)*Cdensity_ice
          ssnow%wb(i,k)         = ssnow%wbliq(i,k) + ssnow%wbice(i,k)
          ssnow%wmtot(i,k)      = ssnow%wmliq(i,k) + ssnow%wmice(i,k)
       END DO
    END DO

    DO i=1,mp
       ssnow%rnof2(i)        = ssnow%qhz(i)               !rnof2 is default cable deep runoff var
    END DO  !mp loop


  END SUBROUTINE smoistgw



  ! Inputs:
  !	 dt_in - time step in sec
  !	 ktau_in - time step no.
  !	 ga	 - ground heat flux W/m^2
  !	 dgdtg	 -
  !	 condxpr - total precip reaching the ground (liquid and solid)
  !	 scondxpr - precip (solid only)
  !	 fev   - transpiration (W/m2)
  !	 fes   - soil evaporation (W/m2)
  !	 isoil - soil type
  !	 ivegt - vegetation type
  ! Output
  !	 ssnow
  SUBROUTINE soil_snow_gw(dels, soil, ssnow, canopy, met, bal, veg)
    USE cable_IO_vars_module, ONLY: wlogn

    USE cable_common_module
USE snow_processes_soil_thermal_mod, ONLY : snow_processes_soil_thermal
    REAL                     , INTENT(IN)     :: dels ! integration time step (s)
    TYPE(soil_parameter_type), INTENT(INOUT)  :: soil
    TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
    TYPE(canopy_type)        , INTENT(INOUT)  :: canopy
    TYPE(veg_parameter_type) , INTENT(INOUT)  :: veg
    TYPE(met_type)           , INTENT(INOUT)  :: met ! all met forcing
    TYPE (balances_type)     , INTENT(INOUT)  :: bal

    INTEGER             :: k,i
    REAL, DIMENSION(mp) :: snowmlt
    REAL, DIMENSION(mp) :: GWwb_ic
    REAL, DIMENSION(mp,ms) :: tgg_old
    REAL, DIMENSION(mp) :: tggsn_old,wbtot_ic,del_wbtot
    REAL(r_2), DIMENSION(mp) :: xx
    REAL                :: zsetot
    INTEGER, SAVE :: ktau =0
    REAL(r_2) :: wb_lake_T, rnof2_T
    LOGICAL :: use_sli
    LOGICAL, SAVE :: first_gw_hydro_call = .TRUE.

    use_sli = .FALSE.

    ktau = ktau +1

    zsetot = SUM(soil%zse)
    ssnow%tggav = 0.

    DO k = 1, ms

       ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot

    END DO


    IF( cable_runtime%offline .OR. cable_runtime%mk3l ) ssnow%t_snwlr = 0.05_r_2

    DO i=1,mp
       ssnow%fwtop1(i) = 0.0
       ssnow%fwtop2(i) = 0.0
       ssnow%fwtop3(i) = 0.0
       ssnow%runoff(i) = 0.0 ! initialise total runoff
       ssnow%rnof1(i) = 0.0 ! initialise surface runoff
       ssnow%rnof2(i) = 0.0 ! initialise deep drainage
       ssnow%smelt(i) = 0.0 ! initialise snowmelt
       ssnow%dtmlt(i,:) = 0.0
       ssnow%osnowd(i) = ssnow%snowd(i)
       ! Scaling  runoff to kg/m^2/s (mm/s) to match rest of the model
       ssnow%sinfil(i) = 0.0
       ssnow%qhz(i) = 0.0
    END DO

    IF (cable_user%soil_thermal_fix) THEN
       soil%heat_cap_lower_limit(:,:) = 0.01  !never allow /0
    ELSE
       soil%heat_cap_lower_limit(:,:) = soil%css_vec(:,:) * soil%rhosoil_vec(:,:)
    END IF

    IF( (.NOT.cable_user%cable_runtime_coupled ) .AND. (first_gw_hydro_call)) THEN

       IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
       ! after discussion with BP
       ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
       ssnow%wbtot = 0.0
       ssnow%wb(:,:)  = MIN( soil%ssat_vec(:,:), MAX ( ssnow%wb(:,:), 0.5*soil%swilt_vec(:,:) ) )

       DO k = 1, ms

          WHERE( ssnow%tgg(:,k) <= CTFRZ .AND. ssnow%wbice(:,k) <= 0.001*ssnow%wb(:,k) )   &
               ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)

          !WHERE( ssnow%tgg(:,k) < CTFRZ)                                    &
          !   ssnow%wbice(:,k) = 0.8 * ssnow%wb(:,k)

       END DO

       WHERE ( soil%isoilm .EQ. 9)! .and. ssnow%snowd .le. 0.1*max_glacier_snowd)

          ! permanent ice: fix hard-wired number in next version
          ssnow%snowd = max_glacier_snowd
          ssnow%osnowd = max_glacier_snowd
          ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0

       END WHERE

       WHERE ( SPREAD(soil%isoilm,2,ms) .EQ. 9 )

          ssnow%wb    = 0.95 * soil%ssat_vec
          ssnow%wbice = 0.95 * ssnow%wb

       END WHERE

    END IF


    tgg_old = ssnow%tgg

    !Start with wb and wbice.  Need wbliq, wmliq,wmice,wmtot
    !find the mass of ice and liq from the prognostic volumetric values
    DO k=1,ms
       DO i=1,mp
          ssnow%wbliq(i,k) = ssnow%wb(i,k) - ssnow%wbice(i,k)                     !liquid volume
          ssnow%wmice(i,k) = ssnow%wbice(i,k)*REAL(Cdensity_ice*soil%zse(k),r_2) !ice mass
          ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*REAL(Cdensity_liq*soil%zse(k),r_2) !liquid mass
          ssnow%wmtot(i,k) = ssnow%wmice(i,k) + ssnow%wmliq(i,k)                  !liq+ice mass
          ssnow%wblf(i,k)   = MAX(ssnow%wbliq(i,k)/soil%ssat_vec(i,k),0.01_r_2)
          ssnow%wbfice(i,k) = MAX(ssnow%wbice(i,k)/soil%ssat_vec(i,k),0._r_2)

       END DO
    END DO


    IF( first_gw_hydro_call ) THEN

       DO i=1,mp
          ssnow%gammzz(i,1) = MAX(soil%heat_cap_lower_limit(i,1),&
               (1.0-soil%ssat_vec(i,1))*&
               soil%css_vec(i,1) * soil%rhosoil_vec(i,1)  &
               & + ssnow%wbliq(i,1) * Ccs_rho_wat           &
               & + ssnow%wbice(i,1) * Ccs_rho_ice ) * soil%zse(1) +   &
               & (1. - ssnow%isflag(i)) * Ccgsnow * ssnow%snowd(i)

       END DO

    ENDIF  ! if(.NOT.cable_runtime_coupled) and first_gw_hydro_call


    DO i=1,mp
       !initial water in the soil column
       wbtot_ic(i)  = SUM(ssnow%wbliq(i,:)*Cdensity_liq*soil%zse(:),1) + &
            SUM(ssnow%wbice(i,:)*Cdensity_ice*soil%zse(:),1) + &
            ssnow%GWwb(i)*soil%GWdz(i)*Cdensity_liq

       GWwb_ic(i) = ssnow%GWwb(i)

    END DO

    !improve hiding, call single soilsnow subroutine to do all the
    !snow processes and thermal soil calculations

    CALL snow_processes_soil_thermal(dels,ssnow,soil,veg,canopy,met,bal,snowmlt)

    !leave here for now, could move into soilsnow as well
    CALL remove_transGW(dels, soil, ssnow, canopy, veg)        !transpiration loss per soil layer

    CALL  GWsoilfreeze(dels, soil, ssnow,tgg_old)

    ssnow%fwtop = canopy%precis/dels + ssnow%smelt/dels   !water from canopy and snowmelt [mm/s]

    CALL iterative_wtd (ssnow, soil, veg, .TRUE. )

    CALL ovrlndflx (dels, ssnow, soil, veg, canopy,use_sli )         !surface runoff, incorporate ssnow%pudsto?

    ssnow%sinfil = ssnow%fwtop - canopy%segg  !canopy%fes/Chl               !remove soil evap from throughfall

    CALL smoistgw (dels,ktau,ssnow,soil,veg,canopy)               !vertical soil moisture movement.

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

    DO i=1,mp
       ssnow%pudsto(i) = 0.0  !no puddle
       ssnow%smelt(i)  = ssnow%smelt(i)/dels    !change units to mm/s.  cable_driver then reverts back to mm
       ssnow%runoff(i) = (ssnow%rnof1(i) + ssnow%rnof2(i))!*dels  !cable_driver converts from mm/s to mm
       !rnof1 and rnof2 are already in mm/s
       ! Set weighted soil/snow surface temperature
       ssnow%tss(i) =  (1-ssnow%isflag(i))*ssnow%tgg(i,1) + ssnow%isflag(i)*ssnow%tggsn(i,1)

       !total water mass at the end of the soilsnow_GW routine
       ssnow%wbtot(i)  = SUM(ssnow%wbliq(i,:)*Cdensity_liq*soil%zse(:),dim=1) + &
            SUM(ssnow%wbice(i,:)*Cdensity_ice*soil%zse(:),dim=1) + &
            ssnow%GWwb(i)*soil%GWdz(i)*Cdensity_liq

       !for debug water balance.  del_wbtot = fluxes = infiltration [though-evap] - trans - qhorz drainage
       del_wbtot(i)   = dels * (ssnow%sinfil(i) - ssnow%rnof2(i) - canopy%fevc(i) / Chl)
       !set below to keep track of water imbalance within the GW module explicitly.  also must change cable_checks
       !ssnow%wbtot(i) = ssnow%wbtot(i)-wbtot_ic(i)

    END DO

  END SUBROUTINE soil_snow_gw

  SUBROUTINE calc_equilibrium_water_content(ssnow,soil)
    !find layer mean soil moisture and potential at equilibrium with wtd

    IMPLICIT NONE

    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    !local variables
    REAL(r_2), DIMENSION(mp)    :: zaq      !node depth of the aquifer
    REAL(r_2), DIMENSION(ms)    :: dzmm     !layer thickness for single tile
    REAL(r_2), DIMENSION(mp)    :: GWdzmm   !aquifer thickness at each tile
    REAL(r_2), DIMENSION(mp)    :: GWzimm   !aquifer layer interface depth
    REAL(r_2), DIMENSION(0:ms)  :: zimm     !layer interface depth in mm
    REAL(r_2), DIMENSION(ms)    :: zmm      !node depths in mm
    REAL(r_2)                   :: tempi, temp0,voleq1,wbrat
    REAL(r_2), DIMENSION(mp,ms+1) :: ice_correction

    INTEGER :: k,i

    IF (gw_params%ssgw_ice_switch) THEN
       smp_cor = 8.0
    ELSE
       smp_cor = 0.0
    END IF

    !make code cleaner define these here
    dzmm    = 1000.0_r_2 * REAL(soil%zse(:),r_2)
    zimm(0) = 0._r_2

    DO k=1,ms
       zimm(k) = zimm(k-1) + dzmm(k)
       zmm(k)  = zimm(k-1) + 0.5_r_2*dzmm(k)
    END DO

    DO i=1,mp
       GWdzmm(i) = REAL(soil%GWdz(i),r_2)*1000._r_2
       GWzimm(i) = zimm(ms)+GWdzmm(i)
       zaq(i)    = zimm(ms) + 0.5_r_2*GWdzmm(i)
    END DO

    IF (.NOT.gw_params%ssgw_ice_switch) THEN
       ice_correction(:,:) = 1._r_2

    ELSE

       DO k=1,ms
          DO i=1,mp
             ice_correction(i,k)    = 1._r_2 + smp_cor * ssnow%wbice(i,k)
             ice_correction(i,k)    =  ice_correction(i,k)**(2.0/soil%bch_vec(i,k))
          END DO
       END DO
       DO i=1,mp
          ice_correction(i,ms+1)    = 1._r_2 + smp_cor * ssnow%wbice(i,ms)
          ice_correction(i,ms+1)    =  ice_correction(i,ms+1)**(2.0/soil%GWbch_vec(i))
       END DO

    END IF


    !equilibrium water content
    DO k=1,ms
       DO i=1,mp

          IF ((ssnow%wtd(i) .LE. zimm(k-1))) THEN         !fully saturated

             ssnow%wbeq(i,k) = soil%ssat_vec(i,k)

          ELSEIF ((ssnow%wtd(i) .LE. zimm(k)) .AND. &
               (ssnow%wtd(i) .GT. zimm(k-1))) THEN

             tempi = 1._r_2
             temp0 = &
                  (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(k-1))/&
                  soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
             voleq1 = -soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
                  (1._r_2-1._r_2/soil%bch_vec(i,k))/&
                  (ssnow%wtd(i)-zimm(k-1))*(tempi-temp0)
             ssnow%wbeq(i,k) = (voleq1*(ssnow%wtd(i)-zimm(k-1)) +&
                  (soil%ssat_vec(i,k)-soil%watr(i,k))&
                  *(zimm(k)-ssnow%wtd(i)))/(zimm(k)-zimm(k-1))&
                  + soil%watr(i,k)
             ssnow%wbeq(i,k) = ssnow%wbeq(i,k)*ice_correction(i,k)
          ELSE

             tempi = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(k))/&
                  soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
             temp0 = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(k-1))/&
                  soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
             ssnow%wbeq(i,k) = -soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
                  (1._r_2-1._r_2/soil%bch_vec(i,k))/&
                  (zimm(k)-zimm(k-1))*(tempi-temp0)+soil%watr(i,k)
             ssnow%wbeq(i,k) = ssnow%wbeq(i,k)*ice_correction(i,k)
          END IF

          ssnow%wbeq(i,k) = MIN(MAX(ssnow%wbeq(i,k),soil%watr(i,k)),soil%ssat_vec(i,k))

          ssnow%wbeq(i,k) = ssnow%wbeq(i,k)*ice_correction(i,k)

          wbrat = MIN(MAX((&
               ssnow%wbeq(i,k) - soil%watr(i,k))/(soil%ssat_vec(i,k)-soil%watr(i,k)),&
               0.001_r_2),1._r_2)

          ssnow%zq(i,k) = MAX(&
               -soil%sucs_vec(i,k)*(wbrat**(-soil%bch_vec(i,k))),sucmin)

          IF (gw_params%ssgw_ice_switch) THEN
             ssnow%zq(i,k) = MAX(&
                  -soil%sucs_vec(i,k)*(1._r_2+smp_cor*ssnow%wbice(i,k))*(wbrat**(-soil%bch_vec(i,k))),sucmin)
          END IF


       END DO  !mp
    END DO  !ms

    DO i=1,mp
       !Aquifer Equilibrium water content
       IF (ssnow%wtd(i) .LE. zimm(ms)) THEN      !fully saturated

          ssnow%GWwbeq(i) = soil%GWssat_vec(i)-soil%GWwatr(i)

       ELSEIF ((ssnow%wtd(i) .GT. GWzimm(i)))   THEN     !fully unsaturated

          tempi = &
               (((soil%GWsucs_vec(i)+ssnow%wtd(i)-GWzimm(i))/&
               soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
          temp0 = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(ms))/&
               soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
          ssnow%GWwbeq(i) = -soil%GWsucs_vec(i)*soil%GWssat_vec(i)/&
               (1._r_2-1._r_2/soil%GWbch_vec(i))/&
               (GWzimm(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)

       ELSE

          tempi  = 1._r_2
          temp0  = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(ms))/&
               soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
          voleq1 = -soil%GWsucs_vec(i)*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
               (1._r_2-1._r_2/soil%GWbch_vec(i))/&
               (ssnow%wtd(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)
          ssnow%GWwbeq(i) = (voleq1*(ssnow%wtd(i)-zimm(ms)) + &
               (soil%GWssat_vec(i)-soil%GWwatr(i))*&
               (GWzimm(i)-ssnow%wtd(i)))/(GWzimm(i)-zimm(ms)) + soil%GWwatr(i)

       END IF

       ssnow%GWwbeq(i) = MIN(MAX(ssnow%GWwbeq(i),soil%GWwatr(i)),soil%GWssat_vec(i))

       ssnow%GWzq(i) = -soil%GWsucs_vec(i)*(MAX((ssnow%GWwbeq(i)-soil%GWwatr(i))/     &
            (soil%GWssat_vec(i)-soil%GWwatr(i)),0.001_r_2))**(-soil%GWbch_vec(i))
       ssnow%GWzq(i) = MAX(sucmin, ssnow%GWzq(i))

    END DO

  END SUBROUTINE calc_equilibrium_water_content

  SUBROUTINE calc_srf_wet_fraction(ssnow,soil,met,veg)

    IMPLICIT NONE
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters
    TYPE (met_type), INTENT(IN)       :: met
    TYPE (veg_parameter_type), INTENT(IN)    :: veg

    !local variables
    REAL(r_2), DIMENSION(mp)           :: icef,satfrac_liqice,S
    REAL(r_2)                          :: fice,xx
    REAL(r_2)                          :: dzmm_one,liqmass,icemass,totmass
    INTEGER                            :: i,j,k
    REAL(r_2)                          :: wb_unsat,wb_lin,funcval
    REAL(r_2)                          :: derv,slopeSTDmm,func_step
    REAL(r_2)                          :: wb_evap_threshold

    IF (cable_user%or_evap) THEN

       CALL saturated_fraction(ssnow,soil,veg)

       ssnow%wetfac(:) = 1.0

       DO i=1,mp
          IF( ssnow%snowd(i) > 0.1) ssnow%wetfac(i) = 0.9

          IF ( veg%iveg(i) == 16 .AND. met%tk(i) >= CTFRZ + 5. )   &
               ssnow%wetfac(i) = 1.0 ! lakes: hard-wired number to be removed

          IF( veg%iveg(i) == 16 .AND. met%tk(i) < CTFRZ + 5. )   &
               ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed
       END DO

    ELSEIF (cable_user%gw_model) THEN

       CALL saturated_fraction(ssnow,soil,veg)


       DO i = 1,mp
          dzmm_one  = 1000._r_2 * REAL(soil%zse_vec(i,1),r_2)
          icemass  = ssnow%wbice(i,1) * dzmm_one
          liqmass  = (ssnow%wb(i,1)-ssnow%wbice(i,1)) * dzmm_one
          totmass  = MAX(liqmass+icemass,REAL(1e-2,r_2))
          icef(i)     = MAX(0._r_2,MIN(1._r_2, gw_params%IceBeta*icemass / totmass))
       END DO


       !srf frozen fraction.  should be based on topography
       DO i = 1,mp
          fice = (EXP(-gw_params%IceAlpha*(1._r_2-icef(i)))-&
               EXP(-gw_params%IceAlpha))/&
               (1._r_2-EXP(-gw_params%IceAlpha))
          fice = MIN(1._r_2,MAX(0._r_2,fice))

          satfrac_liqice(i) = fice + (1._r_2-fice)*ssnow%satfrac(i)

          wb_unsat = ((ssnow%wb(i,1)-ssnow%wbice(i,1)) -&
               ssnow%satfrac(i)*soil%ssat_vec(i,1))/(1.-ssnow%satfrac(i))
          wb_unsat = MIN(soil%ssat_vec(i,1),MAX(0.,wb_unsat))

          wb_evap_threshold = MIN( MAX( &
               gw_params%SoilEvapAlpha*soil%sfc_vec(i,1), &
               soil%swilt_vec(i,1) ), soil%ssat_vec(i,1) )

          !Sakguchi and Zeng 2009
          IF (wb_unsat .GE. wb_evap_threshold) THEN
             xx = 1.
          ELSE
             xx = 0.25 * (1._r_2 - COS(Cpi*wb_unsat/(wb_evap_threshold)))**2.0
          END IF

          ssnow%wetfac(i) = MAX(0.0,MIN(1.0,satfrac_liqice(i) +&
               (1. - satfrac_liqice(i))*xx ) )

       END DO

    ELSE  !Default formulation

       !call saturated_fraction(ssnow,soil,veg)
       ssnow%satfrac(:) = 1.0e-8
       ssnow%rh_srf(:)  = 1.0

       ssnow%wetfac = MAX( 1.e-6, MIN( 1.0,&
            ( REAL (ssnow%wb(:,1) ) - soil%swilt/ 2.0 )                  &
            / ( soil%sfc - soil%swilt/2.0 ) ) )

       DO i=1,mp

          IF( ssnow%wbice(i,1) > 0. )&
               ssnow%wetfac(i) = ssnow%wetfac(i) * &
               REAL(MAX( 0.5_r_2, 1._r_2 - MIN( 0.2_r_2, &
               ( ssnow%wbice(i,1) / ssnow%wb(i,1) )**2 ) ) )

          IF( ssnow%snowd(i) > 0.1) ssnow%wetfac(i) = 0.9

          IF ( veg%iveg(i) == 16 .AND. met%tk(i) >= Ctfrz + 5. )   &
               ssnow%wetfac(i) = 1.0 ! lakes: hard-wired number to be removed

          IF( veg%iveg(i) == 16 .AND. met%tk(i) < Ctfrz + 5. )   &
               ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed

       ENDDO
       ! owetfac introduced to reduce sharp changes in dry regions,
       ! especially in offline runs in which there may be discrepancies b/n
       ! timing of precip and temperature change (EAK apr2009)
       ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)

    ENDIF  !or_evap, gw_model, or default wetfac parameterization

  END SUBROUTINE calc_srf_wet_fraction

  !SUBROUTINE calc_soil_hydraulic_props(ssnow,soil,veg)
  !   USE cable_common_module
  !   TYPE(soil_parameter_type), INTENT(INOUT)     :: soil
  !   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
  !   TYPE(veg_parameter_type) , INTENT(INOUT)     :: veg
  !
  !   INTEGER :: i,k,kk
  !
  !   REAL(r_2), DIMENSION(mp) :: s1, &  !temporary variables for calculating hydraulic properties
  !                               s2, &
  !                               s_mid, &
  !                               liq_ratio, &
  !                               Dliq_ratio_Dz
  !
  !   REAL(r_2), DIMENSION(0:ms) :: zimm  !depths at interface between layers
  !   REAL(r_2), dimension(mp,ms) ::wb_temp
  !
  !    !soil matric potential, hydraulic conductivity, and
  !         derivatives of each with respect to water (calculated using total (not liquid))
  !
  !    do k=1,ms
  !       do i=1,mp
  !          ssnow%icefrac(i,k) = ssnow%wbice(i,k)/(max(ssnow%wb(i,k),0.01_r_2))
  !          ssnow%fracice(i,k) = (exp(-gw_params%IceAlpha*(1._r_2-ssnow%icefrac(i,k)))&
  !                               -exp(-gw_params%IceAlpha))/(1._r_2-exp(-gw_params%IceAlpha))
  !       end do
  !    end do
  !
  !    ssnow%fracice(:,:) = max( min( ssnow%fracice, 1._r_2), 0._r_2)
  !
  !    if (gw_params%ssgw_ice_switch) then
  !       wb_temp = ssnow%wbliq
  !    else
  !       wb_temp = ssnow%wb
  !    end if
  !
  !    do k=1,ms-1
  !       kk=k+1
  !       do i=1,mp
  !
  !          s1(i) = 0.5_r_2*(max(wb_temp(i,k)-soil%watr(i,k),0.) + &
  !                           max(wb_temp(i,kk)-soil%watr(i,kk),0.)) / &
  !                         (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + &
  !                         (soil%ssat_vec(i,kk)-soil%watr(i,kk))))
  !
  !          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
  !          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)
  !
  !          ssnow%hk(i,k)    =  s1(i)*s2(i)
  !          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*s2(i)*&
  !                            0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))
  !          if (.not.gw_params%ssgw_ice_switch) then
  !             ssnow%hk(i,k)    = ssnow%hk(i,k)*&
  !                               (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))
  !             ssnow%dhkdw(i,k) = ssnow%dhkdw(i,k)*&
  !                               (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))
  !          end if
  !       end do
  !    end do
  !
  !    liq_ratio(:) = 1.0
  !
  !    k = ms
  !       do i=1,mp
  !
  !          if (gw_params%ssgw_ice_switch) then
  !             liq_ratio(i) =min(1.,max(0.,wb_temp(i,k)/max(ssnow%wb(i,k),1e-6) ) )
  !          end if
  !
  !          s1(i) = 0.5_r_2*(max(wb_temp(i,k)-soil%watr(i,k),0.) + &
  !                           max(liq_ratio(i)*ssnow%GWwb(i)-soil%GWwatr(i),0.)) / &
  !                  (0.5_r_2*(soil%ssat_vec(i,k)-soil%watr(i,k) +&
  !                            soil%GWssat_vec(i)-soil%GWwatr(i)))
  !
  !          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
  !          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)
  !
  !          ssnow%hk(i,k)    = s1(i)*s2(i)
  !          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*&
  !                             s2(i)*0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))
  !          if (.not.gw_params%ssgw_ice_switch) then
  !             ssnow%hk(i,k)    =  (1.-ssnow%fracice(i,k))*ssnow%hk(i,k)
  !             ssnow%dhkdw(i,k) =  (1.-ssnow%fracice(i,k))*ssnow%dhkdw(i,k)
  !          endif
  !       end do
  !
  !    do k=1,ms
  !       do i=1,mp
  !          s_mid(i) = (ssnow%wb(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
  !              (soil%ssat_vec(i,k)-soil%watr(i,k))
  !
  !          s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)
  !
  !          ssnow%smp(i,k) = -soil%sucs_vec(i,k)*s_mid(i)**(-soil%bch_vec(i,k))
  !
  !          ssnow%smp(i,k) = max(min(ssnow%smp(i,k),-soil%sucs_vec(i,k)),sucmin)
  !
  !          ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/&
  !                    (max(s_mid(i)*(soil%ssat_vec(i,k)-soil%watr(i,k)),0.001_r_2))
  !       end do
  !    end do
  !
  !    do i=1,mp
  !       !Aquifer properties
  !       s_mid(i) = (ssnow%GWwb(i)*liq_ratio(i)-soil%GWwatr(i))/&
  !                    (soil%GWssat_vec(i)-soil%GWwatr(i))
  !       s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)
  !       s2(i)    = soil%GWhyds_vec(i)*s_mid(i)**(2._r_2*soil%GWbch_vec(i)+2._r_2)
  !
  !       ssnow%GWhk(i)     =s_mid(i)*s2(i)
  !
  !       ssnow%GWdhkdw(i)  =  (2._r_2*soil%GWbch_vec(i)+3._r_2)*&
  !                           s2(i)*0.5_r_2/(soil%GWssat_vec(i)-soil%GWwatr(i))
  !
  !       if (.not.gw_params%ssgw_ice_switch) then
  !          ssnow%GWhk(i)    = (1.-ssnow%fracice(i,ms)) * ssnow%GWhk(i)
  !          ssnow%GWdhkdw(i) = (1.-ssnow%fracice(i,ms)) * ssnow%GWdhkdw(i)
  !       endif
  !
  !       s_mid(i) = (ssnow%GWwb(i)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
  !       s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)
  !
  !       ssnow%GWsmp(i)    = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))
  !       ssnow%GWsmp(i)    = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)
  !       ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/&
  !                            (s_mid(i)*(soil%GWssat_vec(i)-soil%GWwatr(i)))
  !    end do
  !
  !END SUBROUTINE calc_soil_hydraulic_props

  SUBROUTINE calc_soil_hydraulic_props(ssnow,soil,veg)
    USE cable_common_module
    TYPE(soil_parameter_type), INTENT(INOUT)     :: soil
    TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
    TYPE(veg_parameter_type) , INTENT(INOUT)     :: veg

    INTEGER :: i,k,kk

    REAL(r_2), DIMENSION(mp) :: s1, &  !temporary variables for calculating hydraulic properties
         s2, &
         s_mid, &
         liq_ratio, &
         Dliq_ratio_Dz

    REAL(r_2), DIMENSION(0:ms) :: zimm  !depths at interface between layers
    REAL(r_2), DIMENSION(mp,ms+1) ::wb_temp
    REAL(r_2), DIMENSION(mp,ms+1) :: hk_ice_factor

    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water (calculated using total (not liquid))

    DO k=1,ms
       DO i=1,mp
          ssnow%icefrac(i,k) = ssnow%wbice(i,k)/(MAX(ssnow%wb(i,k),0.01_r_2))
          ssnow%fracice(i,k) = (EXP(-gw_params%IceAlpha*(1._r_2-ssnow%icefrac(i,k)))&
               -EXP(-gw_params%IceAlpha))/(1._r_2-EXP(-gw_params%IceAlpha))
       END DO
    END DO

    ssnow%fracice(:,:) = MAX( MIN( ssnow%fracice, 1._r_2), 0._r_2)

    IF (gw_params%ssgw_ice_switch) THEN
       wb_temp(:,1:ms) =  ssnow%wbliq(:,:)
       wb_temp(:,ms+1) = ssnow%GWwb(:)
       smp_cor = 8.0
       DO k=1,ms
          kk = MIN(k+1,ms)
          DO i=1,mp
             IF (soil%isoilm(i) .EQ. 9) THEN
                hk_ice_factor(i,k) = 10.0**(-gw_params%ice_impedence)
             ELSE
                hk_ice_factor(i,k) = 10.0**(-gw_params%ice_impedence* &
                     ( 0.5*(ssnow%wbice(i,k)/MAX(1.0e-8,ssnow%wb(i,k)) + &
                     ssnow%wbice(i,kk)/MAX(1.0e-8,ssnow%wb(i,kk))) ) &
                     )
             END IF
          END DO
       END DO
       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
    ELSE
       wb_temp(:,1:ms) = ssnow%wb(:,:)
       wb_temp(:,ms+1) = ssnow%GWwb(:)
       smp_cor = 0.0
       DO k=1,ms
          kk = MIN(k+1,ms)
          DO i=1,mp
             hk_ice_factor(i,k) = (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))
          END DO
       END DO
       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
    END IF

    liq_ratio(:) = 1.0
    k = ms
    DO i=1,mp
       liq_ratio(i) =MIN(1.,MAX(0.,wb_temp(i,k)/MAX(ssnow%wb(i,k),1e-6) ) )
    END DO
    !aquifer ice
    wb_temp(:,ms+1) = liq_ratio(:) * wb_temp(:,ms+1)

    !potential from soil water rention function
    !defined as layer average
    DO k=1,ms
       DO i=1,mp
          s_mid(i) = (wb_temp(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
               (soil%ssat_vec(i,k)-soil%watr(i,k))

          s_mid(i) = MIN(MAX(s_mid(i),0.001_r_2),1._r_2)

          ssnow%smp(i,k) = -soil%sucs_vec(i,k)*s_mid(i)**(-soil%bch_vec(i,k))*&
               ((1._r_2 + smp_cor*ssnow%wbice(i,k))**2.0)

          ssnow%smp(i,k) = MAX(MIN(ssnow%smp(i,k),-soil%sucs_vec(i,k)),sucmin)

          ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/&
               (MAX(s_mid(i)*(soil%ssat_vec(i,k)-soil%watr(i,k)),0.001_r_2)) *&
               ((1._r_2 + smp_cor*ssnow%wbice(i,k))**2.0)
       END DO
    END DO

    !Aquifer potential
    DO i=1,mp
       s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/&
            (soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = MIN(MAX(s_mid(i),0.001_r_2),1._r_2)
       s2(i)    = soil%GWhyds_vec(i)*s_mid(i)**(2._r_2*soil%GWbch_vec(i)+2._r_2)

       !ssnow%GWhk(i)     =s_mid(i)*s2(i) * hk_ice_factor(i,ms+1)
       !ssnow%GWdhkdw(i)  =  (2._r_2*soil%GWbch_vec(i)+3._r_2)*&
       !                    s2(i)*0.5_r_2/(soil%GWssat_vec(i)-soil%GWwatr(i)) *&
       !                    hk_ice_factor(i,ms+1)

       ssnow%GWhk(i)     = soil%GWhyds_vec(i) * hk_ice_factor(i,ms+1)*&
            EXP(-ssnow%wtd(i)/1000._r_2/&
            (1.0/(120*(soil%drain_dens(i)+1.0e-3))))
       !d(h)*Sy=dW
       ssnow%GWdhkdw(i)  = ssnow%GWhk(i)/(soil%GWssat_vec(i)-soil%GWwatr(i))*&
            (0.001/(120*(soil%drain_dens(i)+1.0e-3)))

       s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = MIN(MAX(s_mid(i),0.001_r_2),1._r_2)

       ssnow%GWsmp(i)    = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))
       ssnow%GWsmp(i)    = MAX(MIN(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)
       ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/&
            (s_mid(i)*(soil%GWssat_vec(i)-soil%GWwatr(i)))
    END DO

    !hydraulic conductivity
    !Interfacial so uses layer i and i+1
    DO k=1,ms
       kk=MIN(ms+1,k+1)
       DO i=1,mp

          IF (k .LT. ms) THEN
             s1(i) = 0.5_r_2*((wb_temp(i,k)-soil%watr(i,k)) + &
                  (wb_temp(i,kk)-soil%watr(i,kk))) / &
                  (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + &
                  (soil%ssat_vec(i,kk)-soil%watr(i,kk))))
          ELSE
             s1(i) = 0.5_r_2*((wb_temp(i,k)-soil%watr(i,k)) + &
                  (wb_temp(i,kk)-soil%GWwatr(i))) / &
                  (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + &
                  (soil%GWssat_vec(i)-soil%GWwatr(i))))
          END IF
          s1(i) = MIN(MAX(s1(i),0.01_r_2),1._r_2)
          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)

          ssnow%hk(i,k)    =  s1(i)*s2(i)*hk_ice_factor(i,k)
          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*s2(i)*&
               0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))*&
               hk_ice_factor(i,k)
       END DO
    END DO

    !conductivity between soil column and aquifer
    k = ms
    kk = ms+1
    DO i=1,mp

       s1(i) = 0.5_r_2*(MAX(wb_temp(i,k )-soil%watr(i,k),0.) + &
            MAX(wb_temp(i,kk)-soil%GWwatr(i),0.)) / &
            (0.5_r_2*(soil%ssat_vec(i,k)-soil%watr(i,k) +&
            soil%GWssat_vec(i)-soil%GWwatr(i)))

       s1(i) = MIN(MAX(s1(i),0.01_r_2),1._r_2)
       s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)

       ssnow%hk(i,k)    = s1(i)*s2(i)*hk_ice_factor(i,k)
       ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*&
            hk_ice_factor(i,k)*&
            s2(i)*0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))
    END DO


  END SUBROUTINE calc_soil_hydraulic_props


  SUBROUTINE aquifer_recharge(dt,ssnow,soil,veg,zaq,zmm,dzmm)
    USE cable_common_module

    IMPLICIT NONE
    REAL, INTENT(in) :: dt
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    REAL(r_2), DIMENSION(:), INTENT(in)       :: zaq
    REAL(r_2), DIMENSION(:), INTENT(in)       :: zmm,dzmm

    INTEGER :: i
    !Doing the recharge outside of the soln of Richards Equation makes it easier to track total recharge amount.
    !Add to ssnow at some point
    DO i=1,mp
       IF ((ssnow%wtd(i) .LE. SUM(dzmm,dim=1)) .OR. &
            (veg%iveg(i) .GE. 16) .OR. &
            (soil%isoilm(i) .EQ. 9))  THEN

          ssnow%Qrecharge(i) = 0._r_2
       ELSE
          ssnow%Qrecharge(i) = -0.5*(ssnow%hk(i,ms)*ssnow%GWhk(i))*&
               ((-ssnow%smp(i,ms)) -&
               (-ssnow%zq(i,ms))) / &
               (ssnow%wtd(i) - &
               (zmm(ms)-0.5*soil%zse_vec(i,ms)*1000.0))
       END IF
    END DO


  END SUBROUTINE aquifer_recharge

  SUBROUTINE subsurface_drainage(ssnow,soil,veg,dzmm)
    USE cable_common_module

    IMPLICIT NONE

    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    REAL(r_2), DIMENSION(:), INTENT(in)       :: dzmm
    REAL(r_2), DIMENSION(mp)                  :: sm_tot,&!sum var
         ice_factor_tot  !avg ice factor
    REAL(r_2), DIMENSION(mp,ms)               :: ice_factor  !ice limitation on
    !subsurface drainage
    INTEGER, DIMENSION(mp)                    :: k_drain
    INTEGER :: i,k

    REAL(r_2), DIMENSION(17) :: Efold_mod

    Efold_mod(:) = 1.0
    !Efold_mod(1:4) = (/0.2,0.2,0.2,0.2/)
    !Efold_mod(9) = 0.25
    DO i=1,mp

       !Note: future revision will have interaction with river here. nned to
       !work on router and add river type cells
       ssnow%qhz(i)  = MIN(MAX(soil%slope(i),0.000001),0.9)*&
            gw_params%MaxHorzDrainRate* &
            EXP(-ssnow%wtd(i)/1000._r_2/&
            (1.0/(60.0*gw_params%EfoldHorzDrainRate*&
            (soil%drain_dens(i)+1.0e-3))*Efold_mod(veg%iveg(i))))


       IF (gw_params%subsurface_sat_drainage) THEN
          !drain from sat layers
          k_drain(i) = ms+1
          DO k=ms,2,-1
             IF (ssnow%wtd(i) .LE. SUM(dzmm(1:k),dim=1)) THEN
                k_drain(i) = k + 1
             END IF
          END DO
          k_drain(i) = MAX(k_drain(i),3)
       ELSE
          k_drain(i) = 2
       END IF

    END DO

    IF (gw_params%ssgw_ice_switch) THEN
       DO k=1,ms
          DO i=1,mp
             ice_factor(i,k) = (10.0**(-gw_params%ice_impedence*ssnow%wbice(i,k)/&
                  (ssnow%wb(i,k)+1.0e-12)))
          END DO
       END DO
    ELSE
       DO k=1,ms
          DO i=1,mp
             ice_factor(i,k) = (1._r_2-ssnow%fracice(i,k))
          END DO
       END DO

    END IF
    DO i=1,mp

       ice_factor(i,:)  = 0._r_2
       ice_factor_tot(i)= 0._r_2
       ssnow%qhlev(i,:) = 0._r_2
       sm_tot(i)        = 0._r_2

    END DO

    DO i=1,mp

       IF (gw_params%subsurface_sat_drainage) THEN
          sm_tot(i) = MAX((ssnow%GWwb(i) - soil%watr(i,ms)),0._r_2)*&
               ice_factor(i,ms)

          ice_factor_tot(i) = (SUM(ice_factor(i,k_drain(i):ms)*&
               soil%zse_vec(i,k_drain(i):ms),dim=1)+&
               ice_factor(i,ms)*soil%GWdz(i))/&
               (SUM(soil%zse_vec(i,:),dim=1)+&
               soil%GWdz(i))

          DO k=k_drain(i),ms
             sm_tot(i) = sm_tot(i) + MAX(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)*&
                  ice_factor(i,k)
          END DO

          ssnow%qhz(i) = ssnow%qhz(i)*ice_factor_tot(i)  !reduced due to ice

          IF (sm_tot(i) .GE. 1.0e-12) THEN
             DO k=k_drain(i),ms
                ssnow%qhlev(i,k) = ssnow%qhz(i)*ice_factor(i,k)*&
                     MAX(0.0,ssnow%wbliq(i,k)-soil%watr(i,k))/sm_tot(i)
             END DO
             ssnow%qhlev(i,ms+1) = MAX((ssnow%GWwb(i) - soil%watr(i,ms)),0.0)*&
                  ice_factor(i,ms)*ssnow%qhz(i)/sm_tot(i)
          ENDIF

       ELSE  !second option
          IF (k_drain(i) .LE. ms) THEN
             DO k=k_drain(i),ms
                sm_tot(i) = sm_tot(i) + MAX(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)*&
                     ice_factor(i,k)
             END DO
             ice_factor_tot(i) = (SUM(ice_factor(i,k_drain(i):ms)*&
                  soil%zse_vec(i,k_drain(i):ms),dim=1))/&
                  (SUM(soil%zse_vec(i,:),dim=1))

             ssnow%qhz(i) = ssnow%qhz(i)*ice_factor_tot(i)  !reduced due to ice
             IF (sm_tot(i) .GE. 1.0e-12) THEN
                DO k=k_drain(i),ms
                   ssnow%qhlev(i,k) = (ssnow%qhz(i)*ice_factor(i,k)/sm_tot(i))*&
                        MAX(0.0,ssnow%wbliq(i,k)-soil%watr(i,k))
                END DO
             ENDIF

          ELSE
             ice_factor_tot(i) = ice_factor(i,ms)
             ssnow%qhz(i)        = ssnow%qhz(i)*ice_factor_tot(i)
             ssnow%qhlev(i,ms+1) = ssnow%qhz(i)*MAX(ssnow%GWwb(i)-soil%watr(i,ms),0.0)

          END IF

       END IF
       !incase every layer is frozen very dry
       ssnow%qhz(i) = SUM(ssnow%qhlev(i,:),dim=1)

       !Keep "lakes" saturated forcing qhz = 0.  runoff only from lakes
       !overflowing
       IF (soil%isoilm(i) .EQ. 9 .OR. veg%iveg(i) .GE. 16) THEN
          ssnow%qhz(i) = 0._r_2
          ssnow%qhlev(i,:) = 0._r_2
       END IF

    END DO


  END SUBROUTINE subsurface_drainage


  SUBROUTINE saturated_fraction(ssnow,soil,veg)
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE(veg_parameter_type) , INTENT(IN)    :: veg  ! veg parameters

    REAL(r_2), DIMENSION(mp) :: S
    REAL(r_2) :: slopeSTDmm
    INTEGER :: i,k

    !if !gw_model and !or_evap:
    !in cable_um_init_subrs.F90 or cable_parameters:
    !  ssat_vec(i,:) = ssat
    !  zse_vec(i,:)  = zse
    !  UM: slope_std read in
    ! offline: slope_std read in or set to const
    !  all gw_params set by default in cable_common
    ! doesn do anything but cannot hurt


    S(:) = 0._r_2
    DO k=1,gw_params%level_for_satfrac
       S(:) = S(:) + MAX(0.01,MIN(1.0, &
            (ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/&
            MAX(0.001,soil%ssat_vec(:,k)-soil%watr(:,k)) ) )*soil%zse_vec(:,k)
    END DO
    S(:) = S(:)/SUM(soil%zse(1:gw_params%level_for_satfrac),dim=1)
    !srf frozen fraction.  should be based on topography
    DO i = 1,mp
       !Saturated fraction
       IF (gw_params%MaxSatFraction .GT. 1e-7 .AND. veg%iveg(i) .LT. 16) THEN
          slopeSTDmm = SQRT(MIN(MAX(&
               gw_params%MaxSatFraction*soil%slope_std(i),&
               1e-5),10000._r_2)) ! ensure some variability
          ssnow%satfrac(i)    = MAX(0._r_2,MIN(0.99_r_2,&
                                !note UM wants std03, and erf is not included then
               1._r_2 - my_erf( slopeSTDmm / SQRT(2.0* S(i)) ) ) )
       ELSEIF (veg%iveg(i) .LT. 16) THEN
          ssnow%satfrac(i) = 0._r_2
       ELSE
          ssnow%satfrac(i) = 0.975
       END IF
    END DO


  END SUBROUTINE saturated_fraction

  SUBROUTINE pore_space_relative_humidity(ssnow,soil,veg)
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE(veg_parameter_type), INTENT(INOUT)      :: veg

    REAL(r_2), DIMENSION(mp) :: unsat_wb,unsat_smp

    ! Need a matching array of ones to use in Mark's call to the intrinsic
    ! sign func below
    ! mgk, 24/07/2018
    REAL(r_2), DIMENSION(mp,ms) :: ones

    INTEGER :: i

    !if gw_model = true
    !cable_um_init_subrs.F90 or cable_parameters:
    ! ssat(i) = ssat_vec(i,1)
    !if gw_model = false
    !cable_um_init_subrs.F90 or cable_parameters:
    !ssat_vec(i,:) = ssat(i)
    !so ssat_vec can be used although soilsnow uses ssat

    DO i=1,mp
       IF (veg%iveg(i) .LT. 16 .AND. soil%isoilm(i) .NE. 9 .AND. &
            ssnow%snowd(i) .LE. 1e-8 ) THEN

          unsat_wb(i) =  (ssnow%wb(i,1) - soil%ssat_vec(i,1)*&
               MIN(0.95,MAX(0.0,ssnow%satfrac(i))))/(1.0 - MIN(0.95,MAX(0.0,ssnow%satfrac(i))))

          unsat_wb(i) = MAX(soil%watr(i,1)+1e-2, MIN(soil%ssat_vec(i,1), unsat_wb(i) ) )

          !unsat_smp(i) = sign(soil%sucs_vec(i,1),-1.0) * min(1.0, &
          !                         (max(0.001, (unsat_wb(i)-soil%watr(i,1))/(soil%ssat_vec(i,1)-&
          !                         soil%watr(i,1)) ) )** (-soil%bch_vec(i,1)) )

          ! mgk, 24/07/2018 - fix to compile
          unsat_smp(i) = SIGN(soil%sucs_vec(i,1),ones(i,1)) * MIN(1.0, &
               (MAX(0.001, (unsat_wb(i)-soil%watr(i,1))/(soil%ssat_vec(i,1)-&
               soil%watr(i,1)) ) )** (-soil%bch_vec(i,1)) )

          unsat_smp(i) = MAX(-1.0e7,unsat_smp(i) )/1000._r_2 !m

          ssnow%rh_srf(i) = MAX(0.,MIN(1., &
               EXP(9.81*unsat_smp(i)/(ssnow%tgg(i,1)*461.4)) ) )

       ELSE

          ssnow%rh_srf(i) = 1.0

       END IF
    END DO


  END SUBROUTINE pore_space_relative_humidity

  SUBROUTINE sli_hydrology(dels,ssnow,soil,veg,canopy)
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil ! soil parameters
    TYPE(veg_parameter_type) , INTENT(INOUT)    :: veg  ! veg parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy

    LOGICAL, SAVE :: sli_call = .TRUE.

    REAL(r_2), DIMENSION(ms) :: dzmm
    REAL(r_2), DIMENSION(mp) :: zmm
    REAL(r_2), DIMENSION(mp) :: zaq

    CALL iterative_wtd (ssnow, soil, veg, cable_user%test_new_gw)

    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    CALL  ovrlndflx (dels, ssnow, soil, veg,canopy,sli_call )

    dzmm = REAL(soil%zse(:),r_2)*1000._r_2

    CALL subsurface_drainage(ssnow,soil,veg,dzmm)

    zmm(:) = 1000._r_2*(SUM(REAL(soil%zse,r_2),dim=1))
    zaq(:) = zmm(:) + 0.5_r_2*soil%GWdz(:)*1000._r_2

    CALL aquifer_recharge(dels,ssnow,soil,veg,zaq,zmm,zmm)




  END SUBROUTINE sli_hydrology


  SUBROUTINE set_unsed_gw_vars(ssnow,soil,canopy)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil ! soil parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy

    ssnow%qhlev = 0.
    ssnow%Qrecharge = 0.
    ssnow%fwtop = 0.
    ssnow%wtd = 0.
    ssnow%satfrac = 1.0
    ssnow%qhz = 0.
    ssnow%wbliq = ssnow%wb - ssnow%wbice
    canopy%sublayer_dz = 0.0
    ssnow%rtevap_sat = 0.0
    ssnow%rtevap_unsat = 0.0

    ssnow%GWwb = 0.9*soil%ssat


  END SUBROUTINE set_unsed_gw_vars

  REAL(r_2) FUNCTION my_erf(x)

    IMPLICIT NONE

    REAL(r_2), INTENT(in) :: x
    REAL(r_2)             :: tmp_val, tmp


    tmp = 1.0 / ( 1.0 + 0.5 * ABS(x) )

    tmp_val =  tmp * EXP(-ABS(x) * ABS(x) - 1.26551223 + tmp *     &
         ( 1.00002368 + tmp * ( 0.37409196 + tmp *          &
         ( 0.09678418 + tmp * (-0.18628806 + tmp *              &
         ( 0.27886807 + tmp * (-1.13520398 + tmp *          &
         ( 1.48851587 + tmp * (-0.82215223 + tmp * 0.17087277 )))))))))

    IF ( x.LT.0.0 ) tmp_val = 2.0 - tmp_val

    my_erf = 1.0 - tmp_val

  END FUNCTION my_erf

END MODULE cable_gw_hydro_module
