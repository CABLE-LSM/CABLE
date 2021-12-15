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
                                   max_glacier_snowd,ktau_gl,psi_c,psi_o

   USE cable_soil_snow_module, ONLY : trimb, snow_processes_soil_thermal

   USE cable_data_module, only: C=>PHYS!issnow_type,point2constants

   IMPLICIT NONE

   PRIVATE

   !TYPE(issnow_type), SAVE :: C

   !mrd561 GW params
   REAL(r_2), PARAMETER :: sucmin       = -1.0e8      ! minimum soil pressure head [mm]
   REAL(r_2), PARAMETER :: volwatmin    = 1e-4        !min soil water [mm]
   REAL(r_2), PARAMETER :: wtd_uncert   = 0.1         ! uncertaintiy in wtd calcultations [mm]
   REAL(r_2), PARAMETER :: wtd_max      = 1000000.0   ! maximum wtd [mm]
   REAL(r_2), PARAMETER :: wtd_min      = 10.0        ! minimum wtd [mm]
   REAL(r_2), PARAMETER :: close_to_one = 0.9999
   REAL(r_2), PARAMETER :: Sy_deep = 0.1
   REAL(r_2), PARAMETER :: dz_deep = 50000.0
   REAL(r_2), PARAMETER :: m2mm = 1000.0
   REAL(r_2), PARAMETER :: mm2m = 0.001
   REAL(r_2), SAVE      :: den_rat=0.921
   REAL(r_2), PARAMETER :: zero=0.0

  INTEGER, PARAMETER :: wtd_iter_max = 20 ! maximum number of iterations to find the water table depth
   abstract interface
      subroutine swc_smp_func(soil,ssnow)
        Import :: soil_parameter_type,soil_snow_type
        type(soil_parameter_type) :: soil
         type(soil_snow_type) :: ssnow
      end subroutine
   end interface

  procedure(swc_smp_func), pointer, save :: swc_smp_dsmpdw=>null()



  ! ! This module contains the following subroutines that
  !are called from other modules
   PUBLIC :: soil_snow_gw,calc_srf_wet_fraction,sli_hydrology,&
             pore_space_relative_humidity,set_unsed_gw_vars,den_rat,set_den_rat

CONTAINS


! -----------------------------------------------------------------------------

SUBROUTINE GWsoilfreeze(dels, soil, ssnow)
  !NOTE: this is only included because gw_model uses parameters XXX_vec
  !these are r_2.  this breaks bitwise compatibility with trunk
  !if acceptable this routine does the same thing but with r_2 soil params
  ! if max_ice_frac always set to frozen_limit and tgg_tmp is always C%TFRZ

   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(INOUT)    :: soil
   REAL     , DIMENSION(mp,ms)        :: tgg_tmp !calc wb ice at approx tgg
   REAL(r_2), DIMENSION(mp,ms)        :: wbice_delta,avail_por,delta_ice_vol
   REAL(r_2), DIMENSION(mp)           :: ice_mass,liq_mass,tot_mass
   INTEGER :: i,j,k
   REAL(r_2) :: func,funcderv,Aconst,Dconst,t_zero,t_one,dtmp
   REAL, DIMENSION(mp,ms) :: gammzz_snow
   REAL(r_2),DIMENSION(mp,ms) :: xx,max_ice_frac,den_css  !Decker and Zeng 2009 ! MMY iceF is useless
   REAL(r_2) :: delta_wbliq,tmp_var

   !call point2constants( C )

   max_ice_frac(:,:) = 0.0
   delta_ice_vol(:,:) = 0.0
   tgg_tmp(:,:) = ssnow%tgg(:,:)
   gammzz_snow(:,:) = 0._r_2

   k=1
   do i=1,mp
      if (ssnow%isflag(i) .eq. 0 .and. soil%isoilm(i) .ne. 9) then
           gammzz_snow(i,k) = real(C%cgsnow,r_2) * real(ssnow%snowd(i),r_2)
      end if
   end do

   do k=1,ms
   do i=1,mp

      if ( (ssnow%tgg(i,k) .lt. C%TFRZ) .and. &
          (ssnow%tgg(i,k) .lt. ssnow%otgg(i,k)) ) then ! MMY freezing

            ssnow%otgg(i,k) = min(ssnow%otgg(i,k),C%TFRZ)

            ! ___________________ MMY not work at all ______________________
            !if (ssnow%wb(i,k) .gt. 1.0e-4) then
            !   iceF(i,k) = 1._r_2 - max(0.1_r_2,min(0.9_r_2,ssnow%wbliq(i,k)/max(ssnow%wb(i,k),1.0e-8)))
            !else
            !   iceF(i,k) = 0._r_2
            !end if
            ! __________________________________________________________________

            tgg_tmp(i,k) = ssnow%tgg(i,k)

            Aconst = -2.0*( (0.2+ssnow%wb(i,k)/soil%ssat_vec(i,k))**4.0 )
            Dconst = exp(1.  - min(1.0,0.2+ssnow%wb(i,k)/soil%ssat_vec(i,k)))

            max_ice_frac(i,k) = (1._r_2 - exp(2._r_2*(min(1.,ssnow%wb(i,k)/soil%ssat_vec(i,k))**4.0) *&
                          real(tgg_tmp(i,k)-C%TFRZ,r_2)))/exp(1._r_2- min(1.0,ssnow%wb(i,k)/soil%ssat_vec(i,k)))

            if (soil%isoilm(i) .eq. 9) max_ice_frac(i,k) = 0.85_r_2

            ! MMY maximum ice can get under current wb
            max_ice_frac(i,k) = min(0.9_r_2,max_ice_frac(i,k))

            !delta_ice_vol(i,k) = max(0._r_2, ssnow%wb(i,k)*max_ice_frac(i,k) - ssnow%wbice(i,k)) ! MMY
            ! MMY ssnow%wb(i,k)*max_ice_frac(i,k) here should be ice volume, thus divided by den_rat
            delta_ice_vol(i,k) = max(0._r_2, ssnow%wb(i,k)*max_ice_frac(i,k)/den_rat - ssnow%wbice(i,k)) ! MMY

            ! MMY wbliq can be lower than watr due to condensation, but to avoid large negative water potential
            !     set watr as wbliq lower boundary
            !check amount of water we have
            delta_ice_vol(i,k) = min((ssnow%wbliq(i,k)-soil%watr(i,k))/den_rat, &
                                 max(0._r_2, delta_ice_vol(i,k) ) )

            delta_ice_vol(i,k) = min(delta_ice_vol(i,k),max(0._r_2,&
                                   real(ssnow%otgg(i,k)-ssnow%tgg(i,k),r_2)*ssnow%gammzz(i,k)/&
                                     (soil%zse_vec(i,k)*real(C%HLF*C%density_ice,r_2)) ) )

      elseif ((ssnow%tgg(i,k) .gt. C%TFRZ) .and. &
              (ssnow%tgg(i,k) .gt. ssnow%otgg(i,k)) .and. ssnow%wbice(i,k) .gt.  0.0) then ! MMY melting

              ssnow%otgg(i,k) = C%TFRZ

             delta_ice_vol(i,k) = ssnow%wbice(i,k)

             delta_ice_vol(i,k) = min(delta_ice_vol(i,k), max(0._r_2,&
                                  real(ssnow%tgg(i,k)-ssnow%otgg(i,k),r_2) * ssnow%gammzz(i,k) /&
                                     (soil%zse_vec(i,k)*real(C%HLF*C%density_ice,r_2)) ) )

      endif
   end do
   end do

   DO k = 1, ms
      DO i=1,mp

          if ( (ssnow%tgg(i,k) .lt. C%TFRZ) .and. &
               (delta_ice_vol(i,k) .gt. 0.0)  ) then

            ssnow%wbice(i,k) = ssnow%wbice(i,k) + delta_ice_vol(i,k)
            ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - delta_ice_vol(i,k)*den_rat

            ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),  &
                              (1.0-soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k) &
                            + ssnow%wbliq(i,k) * REAL(C%cswat*C%density_liq,r_2)   &
                            + ssnow%wbice(i,k) * REAL(C%csice*C%density_ice,r_2)&
                            )*soil%zse_vec(i,k)  + gammzz_snow(i,k)

            ! MMY freezing releases heat
            ssnow%tgg(i,k) = ssnow%tgg(i,k) + real( delta_ice_vol(i,k)*soil%zse_vec(i,k)  *&
                                       real(C%hlf*C%density_ice,r_2) / ssnow%gammzz(i,k) )


            ssnow%wb(i,k) = ssnow%wbliq(i,k) + den_rat*ssnow%wbice(i,k)

            ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*soil%zse_vec(i,k)*C%density_liq
            ssnow%wmice(i,k) = ssnow%wbice(i,k)*soil%zse_vec(i,k)*C%density_ice
            ssnow%wmtot(i,k) = ssnow%wb(i,k)   *soil%zse_vec(i,k)*C%density_liq

          elseif ((ssnow%tgg(i,k) .gt. C%TFRZ) .and. &
              delta_ice_vol(i,k) .gt. 0.0 ) then ! MMY
            ! ssnow%wbice(i,k) .gt. 0.0 ) then ! MMY It's not right.
            ! MMY When previous tgg > current tgg > 0, wbice won't melt

            ssnow%wbice(i,k) = ssnow%wbice(i,k)  -  delta_ice_vol(i,k)
            ssnow%wbliq(i,k) = ssnow%wbliq(i,k)  +  delta_ice_vol(i,k)*den_rat

            ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),  &
                              (1.0-soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k) &
                            + ssnow%wbliq(i,k) * REAL(C%cswat*C%density_liq,r_2)   &
                            + ssnow%wbice(i,k) * REAL(C%csice*C%density_ice,r_2)&
                            )*soil%zse_vec(i,k)  + gammzz_snow(i,k)

            ! MMY melting absorbs heat
            ssnow%tgg(i,k) = ssnow%tgg(i,k) - real( delta_ice_vol(i,k)*soil%zse_vec(i,k)  *&
                                          real(C%hlf*C%density_ice,r_2) / ssnow%gammzz(i,k) )

            ssnow%wb(i,k) = ssnow%wbliq(i,k) + den_rat*ssnow%wbice(i,k)

            ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*soil%zse_vec(i,k)*C%density_liq
            ssnow%wmice(i,k) = ssnow%wbice(i,k)*soil%zse_vec(i,k)*C%density_ice
            ssnow%wmtot(i,k) = ssnow%wb(i,k)   *soil%zse_vec(i,k)*C%density_liq
          END IF


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

   do k=1,ms
      do i=1,mp
         zse_mp_mm(i,k)  = real(soil%zse_vec(i,k)*C%density_liq,r_2)
      end do
   end do

   IF (cable_user%FWSOIL_switch.ne.'Haverd2013') THEN

      xx(:) = 0._r_2
      xxd(:) = 0._r_2
      diff(:,:) = 0._r_2

      DO k = 1,ms

         DO i=1,mp

            if (canopy%fevc(i) .gt. 0._r_2) then

               xx(i) = canopy%fevc(i) * dels / C%hl * veg%froot(i,k) + diff(i,k-1)
               diff(i,k) = max(0._r_2,ssnow%wbliq(i,k)-soil%swilt_vec(i,k)) &
                          * zse_mp_mm(i,k)
               xxd(i) = xx(i) - diff(i,k)

               if (xxd(i) .gt. 0._r_2) then
                  ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - diff(i,k)/zse_mp_mm(i,k)
                  diff(i,k) = xxd(i)
               else
                  ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - xx(i)/zse_mp_mm(i,k)
                  diff(i,k) = 0._r_2
               end if


             end if  !fvec > 0

         END DO  !mp
      END DO     !ms

   ELSE

     WHERE (canopy%fevc .lt. 0.0_r_2)
        canopy%fevw = canopy%fevw+canopy%fevc
        canopy%fevc = 0.0_r_2
     END WHERE
     DO k = 1,ms
        ssnow%wbliq(:,k) = ssnow%wbliq(:,k) - ssnow%evapfbl(:,k)/(soil%zse_vec(:,k)*m2mm)
     ENDDO

  ENDIF

  do k=1,ms
     do i=1,mp
        ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*zse_mp_mm(i,k)!mass
        ssnow%wmtot(i,k) = ssnow%wmliq(i,k) + ssnow%wmice(i,k)  !mass
        ssnow%wb(i,k)    = ssnow%wbliq(i,k) + den_rat * ssnow%wbice(i,k)  !volume
     end do
  end do


END SUBROUTINE remove_transGW
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
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
   dzmm = m2mm * soil%zse(1)

   if (sli_call) then
      do i=1,mp
         if (canopy%through(i) .ge. canopy%through_sn(i)) then
           ssnow%fwtop(i)  = max((canopy%through(i)-canopy%through_sn(i))/dels , 0.)             ! liq precip rate (m s-1)
         else
           ssnow%fwtop(i) = max(canopy%through(i), 0.)
         end if
      end do
   end if
   !amount of ice in surface layer
   do i = 1,mp
      efpor(i) = max(0.01_r_2, soil%ssat_vec(i,1) - den_rat*ssnow%wbice(i,1))
      icemass  = ssnow%wmice(i,1)
      liqmass  = ssnow%wmliq(i,1)
      totmass  = max(liqmass+icemass,real(1e-2,r_2))
      icef(i)     = max(0._r_2,min(1._r_2,icemass / totmass))
   end do

   !sat fraction assuming topo controlled subgrid soil moisture distribution
   !called from cable_canopy for srf wet fraction alrady
   call saturated_fraction(ssnow,soil,veg)

   !srf frozen fraction.  should be based on topography
   do i = 1,mp
      fice = (exp(-gw_params%IceAlpha*(1._r_2-icef(i)))-exp(-gw_params%IceAlpha))/(1._r_2-exp(-gw_params%IceAlpha))
      fice  = min(max(fice,0._r_2),1._r_2)
      satfrac_liqice(i)   = max(0.,min(0.95,fice + (1._r_2-fice)*ssnow%satfrac(i) ) )
   end do

   do i=1,mp
      tmpa = ssnow%wbliq(i,1) / efpor(i)
      tmpb = max( (tmpa-satfrac_liqice(i))/max(0.01_r_2,(1._r_2-satfrac_liqice(i))), 0._r_2)
      tmpa = -2._r_2*soil%bch_vec(i,1)*soil%sucs_vec(i,1)/soil%zse_vec(i,1)/m2mm
      qinmax = (1._r_2 + tmpa*(tmpb-1._r_2))*soil%hyds_vec(i,1)

      ssnow%rnof1(i) = satfrac_liqice(i) * ssnow%fwtop(i) + &
                         (1._r_2-satfrac_liqice(i))*max((ssnow%fwtop(i)-qinmax) , 0._r_2)

      ssnow%fwtop(i) = ssnow%fwtop(i) - ssnow%rnof1(i)

   end do  !mp

  !add back to the lakes to keep saturated instead of drying
  do i=1,mp
     if (veg%iveg(i) .eq. 16) then
        ssnow%fwtop(i) = ssnow%fwtop(i) + ssnow%rnof1(i)
        ssnow%rnof1(i) = 0._r_2
     end if
  end do

   !---  glacier formation
   rnof5= 0.

   if (sli_call .or. cable_runtime%UM .or. cable_user%gw_model) then
      nglacier = 0
   else
     nglacier = 2
   end if

   IF (nglacier == 2) THEN
      smelt1=0.
      WHERE( ssnow%snowd > max_glacier_snowd )

         rnof5 = MIN( 0.1, ssnow%snowd - max_glacier_snowd )

         !---- change local tg to account for energy - clearly not best method
         WHERE( ssnow%isflag == 0 )
            smasstot = 0.0
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - rnof5 * C%hlf                    &
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

  do i=1,mp
     wmean(i) = 0._r_2
     fz(i)    = 5._r_2
     ztot(i)  = 0._r_2
     stot(i,:) = (ssnow%wb(i,:)-soil%watr(i,:)) / (soil%ssat_vec(i,:)-soil%watr(i,:))
  end do
  do k  = 1, ms
     do i=1,mp
        wmean(i) = wmean(i) + stot(i,k)*soil%zse(k)*m2mm
        ztot(i)  = ztot(i) + soil%zse(k)*m2mm
     end do
  end do

  do i=1,mp
     wmean(i) = wmean(i) + ssnow%GWwb(i)/soil%GWssat_vec(i) * soil%GWdz(i)*m2mm
     ztot(i)  = ztot(i) + soil%GWdz(i)*m2mm

     ssnow%wtd(i) = min(200000._r_2, fz(i) * (ztot(i) - wmean(i)))
  end do

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
  REAL(r_2), DIMENSION(mp,ms)   :: tmp_def
  REAL(r_2), DIMENSION(mp)      :: temp
  REAL(r_2), DIMENSION(mp)      :: def,defc,total_depth_column
  REAL(r_2)                     :: deffunc,tempa,tempb,derv,calc,tmpc
  REAL(r_2), DIMENSION(mp)      :: lam,Nsucs_vec  !inverse of C&H B,Nsucs_vec
  INTEGER :: k,i,wttd,jlp

  !make code cleaner define these here
  lam(:)        = 1._r_2/soil%bch_vec(:,ms)                                !1 over C&H B
  Nsucs_vec(:)  = abs(soil%sucs_vec(:,ms))                                !psi_saturated mm

  if (include_aquifer) then  !do we include the aquifer in the calculation of wtd?

     do i=1,mp
        total_depth_column(i) = soil%GWdz(i)*m2mm
        def(i) = max(0._r_2,soil%GWssat_vec(i)-ssnow%GWwb(i))*soil%GWdz(i)*m2mm
     end do

  else
     def(:) = 0._r_2
     total_depth_column(:) = 0._r_2
  end if

  !total depth of soil column
  do k=1,ms
     do i=1,mp
         total_depth_column(i) = total_depth_column(i) + soil%zse_vec(i,k)*m2mm
     end do
  end do

  !comute the total mass away from full saturation
  do k=1,ms
     do i=1,mp

       def(i) = def(i) +                                                           &
                max(0._r_2,(soil%ssat_vec(i,k)-(ssnow%wbliq(i,k)+den_rat*ssnow%wbice(i,k)))*soil%zse_vec(i,k)*m2mm)
      end do  !mp
  end do  !ms

  !find the deficit if the water table is at the bottom of the soil column
  do i=1,mp
     defc(i) = (soil%ssat_vec(i,ms))*(total_depth_column(i)+Nsucs_vec(i)/(1._r_2-lam(i))*            &
             (1._r_2-((Nsucs_vec(i)+total_depth_column(i))/Nsucs_vec(i))**(1._r_2-lam(i))))
     defc(i) = max(0.1_r_2,defc(i))
  end do

  !initial guess at wtd
  ssnow%wtd(:) = total_depth_column(:)*def(:)/defc(:)

 !use newtons method to solve for wtd, note this assumes homogenous column but
 !that is ok
  do i=1,mp
    if ((soil%isoilm(i) .ne. 9) .and. (veg%iveg(i) .ne. 16)) then

      if (defc(i) > def(i)) then                 !iterate tfor wtd

        jlp=0

        mainloop: DO

          tempa   = 1.0_r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-lam(i))
          derv    = (soil%ssat_vec(i,ms))*(tempa-tempb) + &
                                       soil%ssat_vec(i,ms)

          if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

          tempa   = 1.0_r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-lam(i))
          deffunc = (soil%ssat_vec(i,ms))*(ssnow%wtd(i) +&
                           Nsucs_vec(i)/(1-lam(i))* &
                     (tempa-tempb)) - def(i)
          calc    = ssnow%wtd(i) - deffunc/derv

          IF ((abs(calc-ssnow%wtd(i))) .le. wtd_uncert) THEN

            ssnow%wtd(i) = calc
            EXIT mainloop

          ELSEIF (jlp .ge. wtd_iter_max) THEN

            EXIT mainloop

          ELSE

            jlp=jlp+1
            ssnow%wtd(i) = calc

          END IF

        END DO mainloop  !defc .gt. def

      elseif (defc(i) .lt. def(i)) then

        jlp=0

        mainloop2: DO

          tmpc     = Nsucs_vec(i)+ssnow%wtd(i)-total_depth_column(i)
          tempa    = (abs(tmpc/Nsucs_vec(i)))**(-lam(i))
          tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-lam(i))
          derv     = (soil%ssat_vec(i,ms))*(tempa-tempb)
          if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

          tempa    = (abs((Nsucs_vec(i)+ssnow%wtd(i)-total_depth_column(i))/Nsucs_vec(i)))**(1._r_2-lam(i))
          tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-lam(i))
          deffunc  = (soil%ssat_vec(i,ms))*(total_depth_column(i) +&
                     Nsucs_vec(i)/(1._r_2-lam(i))*(tempa-tempb))-def(i)
          calc     = ssnow%wtd(i) - deffunc/derv

          IF ((abs(calc-ssnow%wtd(i))) .le. wtd_uncert) THEN

            ssnow%wtd(i) = calc
            EXIT mainloop2

          ELSEIF (jlp==wtd_iter_max) THEN

            EXIT mainloop2

          ELSE

            jlp=jlp+1
            ssnow%wtd(i) = calc

          END IF

        END DO mainloop2  !defc .lt. def

      else  !water table depth is exactly on bottom boundary

        ssnow%wtd(i) = total_depth_column(i)

      endif

    endif  !check veg and soils

  end do   !mp loop

  !limit wtd to be within a psecified range
  do i=1,mp
     if (veg%iveg(i) .ge. 16) ssnow%wtd(i) = wtd_min
     ssnow%wtd(i) = min(wtd_max,max(wtd_min,ssnow%wtd(i) ) )
  end do


  END SUBROUTINE iterative_wtd

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoistgw (fwtop,dt,ktau,ssnow,soil,prin)
  ! solves the modified richards equation (Zeng and Decker 2009) to find
  ! vertical mocement of soil water.  Bottom boundary condition is determined
  ! using a single layer groundwater module
  !
  SUBROUTINE smoistgw (dels,ktau,ssnow,soil,veg,canopy)
  USE cable_common_module

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
    REAL(r_2), DIMENSION(mp,0:ms+1)        :: zimm
    REAL(r_2), DIMENSION(mp,ms)            :: zmm
    REAL(r_2), DIMENSION(mp)            :: GWzimm,xs,zaq,s_mid,GWdzmm
    REAL(r_2), DIMENSION(mp)            :: xs1,GWmsliq!xsi    !mass (mm) of liquid over/under saturation, mass of aquifer water
    REAL(r_2)                           :: xsi
    REAL(r_2), DIMENSION(mp,ms+1)       :: del_wb
    !MD DEBUG VARS
    INTEGER :: imp,ims,k_drain

    zimm(:,0) = 0._r_2
    do k=1,ms
       zimm(:,k) = zimm(:,k-1) + soil%zse_vec(:,k)*m2mm
    end do
    zmm(:,1:ms)  = zimm(:,1:ms) - 0.5*soil%zse_vec(:,1:ms)*m2mm

    do i=1,mp
       GWzimm(i) = zimm(i,ms)+soil%GWdz(i)*m2mm
       zaq(i)    = zimm(i,ms) + 0.5_r_2*soil%GWdz(i)*m2mm
    end do

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! preset to allow for non-land & snow points in trimb
    do k=1,ms
       do i=1,mp
          old_wb(i,k) = ssnow%wb(i,k)
          rt(i,k) = 0._r_2
          at(i,k) = 0._r_2
          bt(i,k) = 0._r_2
          ct(i,k) = 0._r_2
       end do
    end do

    !soil hydraulic state/props
    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    !equilibrium water content
    CALL calc_equilibrium_water_content(ssnow,soil)

    CALL subsurface_drainage(ssnow,soil,veg)

    k = 1     !top soil layer
    do i=1,mp
       qin(i)     = ssnow%sinfil(i)
       den(i)     = 0.5*(soil%zse_vec(i,k)+soil%zse_vec(i,k+1) )*m2mm
       dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
       num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
       dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
       dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
       rt(i,k) =  qin(i) - qout(i)
       at(i,k) =  0._r_2
       bt(i,k) =  m2mm*soil%zse_vec(i,k)/dels + dqodw1(i)
       ct(i,k) =  dqodw2(i)
    end do
    do k = 2, ms - 1     !middle soil layers
       do i=1,mp
          den(i)     = 0.5*(soil%zse_vec(i,k)+soil%zse_vec(i,k-1) )*m2mm
          dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
          num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
          qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
          dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          den(i)     = 0.5*(soil%zse_vec(i,k)+soil%zse_vec(i,k+1) )*m2mm
          dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
          num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
          qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
          dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
          dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
          rt(i,k) =  qin(i) - qout(i)
          at(i,k) = -dqidw0(i)
          bt(i,k) =  m2mm*soil%zse_vec(i,k)/dels - dqidw1(i) + dqodw1(i)
          ! MMY ??? why the sign od bt is different to CLM5's
          ct(i,k) =  dqodw2(i)
       end do
    end do

    k = ms   !Bottom soil layer
    do i=1,mp
       den(i)     = 0.5*(soil%zse_vec(i,k)+soil%zse_vec(i,k-1) )*m2mm
       dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
       num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
       qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
       dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       den(i)     = zaq(i) - zmm(i,k)
       dne(i)     = (ssnow%GWzq(i)-ssnow%zq(i,k))
       num(i)     =  (ssnow%GWsmp(i)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = 0._r_2
       dqodw1(i)  = 0._r_2
       dqodw2(i)  = 0._r_2
       rt(i,k) =  qin(i) - qout(i)
       at(i,k) = -dqidw0(i)
       bt(i,k) =  m2mm*soil%zse_vec(i,k)/dels - dqidw1(i) + dqodw1(i)
       ct(i,k) =  dqodw2(i)
    end do

    CALL aquifer_recharge(dels,ssnow,soil,veg)

    CALL trimb(at,bt,ct,rt,ms)                       !use the defulat cable tridiag solution

    do k=1,ms
       do i=1,mp
          ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + rt(i,k) - ssnow%qhlev(i,k)*dels/(m2mm*soil%zse_vec(i,k))   !volutermic liquid
       end do
    end do

    do i=1,mp
       ssnow%wbliq(i,ms) = ssnow%wbliq(i,ms) - ssnow%Qrecharge(i)*dels/(m2mm*soil%zse_vec(i,ms))
    end do
    do i=1,mp
       ssnow%GWwb(i) = ssnow%GWwb(i)  +  (ssnow%Qrecharge(i)-ssnow%qhlev(i,ms+1))*dels/(m2mm*soil%GWdz(i))
    end do

    !determine the available pore space
    !volumetric
    do k=1,ms
       do i=1,mp
          eff_por(i,k)  = max(0._r_2, soil%ssat_vec(i,k) - den_rat * ssnow%wbice(i,k) )
       end do
    end do

    do i=1,mp
       xsi = 0._r_2

       if (ssnow%GWwb(i) .gt. soil%GWssat_vec(i)) then
          xsi = (ssnow%GWwb(i) - soil%GWssat_vec(i))*m2mm*soil%GWdz(i)
          ssnow%GWwb(i) = soil%GWssat_vec(i)
       end if

       do k=1,ms
          if (ssnow%wbliq(i,k) .gt. eff_por(i,k)) then
             xsi = xsi + (ssnow%wbliq(i,k) - eff_por(i,k))*(m2mm*soil%zse_vec(i,k))
             ssnow%wbliq(i,k) = eff_por(i,k)
           end if
       end do

       do k = ms,1,-1  !loop from bottom to top adding extra water to each layer
          if (xsi .gt. 0._r_2) then
             if (xsi .lt. (eff_por(i,k)-ssnow%wbliq(i,k))*(m2mm*soil%zse_vec(i,k))) then
                ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + xsi/(m2mm*soil%zse_vec(i,k))
                xsi = 0._r_2
             else
                xsi = xsi - (eff_por(i,k) - ssnow%wbliq(i,k))*(m2mm*soil%zse_vec(i,k))
                ssnow%wbliq(i,k) = eff_por(i,k)
             end if
          end if
       end do  !ms loop

       if (xsi .gt. 0._r_2) then
          ssnow%qhz(i) = ssnow%qhz(i) + xsi/dels
          xsi = 0._r_2
       end if

       ! __________________________________ MMY ________________________________
       ! do k = 1,ms
       !    xsi = 0._r_2             !should be a single float (array not needed)
       !    if (ssnow%wbliq(i,k) .lt. volwatmin) then
       !       xsi = (volwatmin - ssnow%wbliq(i,k))*(m2mm*soil%zse_vec(i,k))  !in mm
       !       ssnow%wbliq(i,k) = volwatmin
       !       if (k .lt. ms) then
       !          ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) - xsi/(m2mm*soil%zse_vec(i,k+1))
       !       else
       !          ssnow%GWwb(i) = ssnow%GWwb(i) - xsi / (m2mm*soil%GWdz(i))
       !       end if
       !    end if
       ! end do  !ms loop
       !
       ! if ( (ssnow%GWwb(i) .lt. volwatmin) .and. (soil%isoilm(i) .ne. 9) ) then
       !    xsi = (volwatmin - ssnow%GWwb(i)) / (m2mm*soil%GWdz(i))  !mm
       !    ssnow%GWwb(i) = volwatmin
       !    ssnow%qhz(i) = ssnow%qhz(i) - xsi / dels
       ! end if

       ! MMY use watr to replace volwatmin to avoid questionable smp or hk (when wbliq<watr)
       do k = 1,ms
          xsi = 0._r_2             !should be a single float (array not needed)
          if (ssnow%wbliq(i,k) .lt. soil%watr(i,k)) then         ! MMY
             xsi = (soil%watr(i,k) - ssnow%wbliq(i,k))*(m2mm*soil%zse_vec(i,k)) ! MMY
             ssnow%wbliq(i,k) = soil%watr(i,k)   ! MMY
             if (k .lt. ms) then
                ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) - xsi/(m2mm*soil%zse_vec(i,k+1))
             else
                ssnow%GWwb(i) = ssnow%GWwb(i) - xsi / (m2mm*soil%GWdz(i))
             end if
          end if
       end do  !ms loop

       if ( (ssnow%GWwb(i) .lt. soil%GWwatr(i)) .and. (soil%isoilm(i) .ne. 9) ) then ! MMY
          xsi = (soil%GWwatr(i) - ssnow%GWwb(i)) / (m2mm*soil%GWdz(i))  !mm ! MMY
          ssnow%GWwb(i) = soil%GWwatr(i) ! MMY
          ssnow%qhz(i) = ssnow%qhz(i) - xsi / dels
          if (ssnow%qhz(i) .lt. 0.) print *, " MMY ===> Soil is too dry, found in SUBROUTINE smoistgw"
       end if
       ! _______________________________________________________________________
   end do

   do k=1,ms
      do i=1,mp

       !update mass variables
         ssnow%wmliq(i,k)      = ssnow%wbliq(i,k) * &
                                         soil%zse_vec(i,k)*real(C%density_liq,r_2)
         ssnow%wmice(i,k)      = ssnow%wbice(i,k) * &
                                         soil%zse_vec(i,k)*real(C%density_ice,r_2)
         ssnow%wb(i,k)         = ssnow%wbliq(i,k) + den_rat * ssnow%wbice(i,k)
         ssnow%wmtot(i,k)      = ssnow%wmliq(i,k) + ssnow%wmice(i,k)
      end do
   end do

   do i=1,mp
       ssnow%rnof2(i)        = ssnow%qhz(i)               !rnof2 is default cable deep runoff var
   end do  !mp loop


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
   !USE cable_IO_vars_module, ONLY: wlogn

   USE cable_common_module
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
   REAL, DIMENSION(mp) :: tggsn_old,wbtot_ic,del_wbtot
   REAL(r_2), DIMENSION(mp) :: xx
   REAL(r_2), DIMENSION(mp,ms) :: gammzz_snow
   REAL                :: zsetot
   INTEGER, SAVE :: ktau =0
   REAL(r_2) :: wb_lake_T, rnof2_T
   LOGICAL :: use_sli
   LOGICAL, SAVE :: first_gw_hydro_call = .true.

   use_sli = .false.

   ktau = ktau +1

   zsetot = sum(soil%zse)
   ssnow%tggav = 0.

   !CALL point2constants( C )

   DO k = 1, ms
      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
   END DO

   if (first_gw_hydro_call) then
      call set_den_rat()
      ! ____________________________ MMY _____________________________
      ! check wb wbice when first_gw_hydro_call = true
      DO k = 1, ms
          ssnow%wb(:,k)    = MIN(soil%ssat_vec(:,k), MAX(real(ssnow%wb(:,k)), soil%watr(:,k)))
          ssnow%wbice(:,k) = MIN(real(ssnow%wb(:,k))/den_rat, real(ssnow%wbice(:,k)))
      END DO
      ! ______________________________________________________________
   end if

   IF( cable_runtime%offline .or. cable_runtime%mk3l ) ssnow%t_snwlr = 0.05_r_2

   ssnow%otgg = ssnow%tgg
   do i=1,mp
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
   end do

   IF (cable_user%soil_thermal_fix) THEN
      soil%heat_cap_lower_limit(:,:) = 0._r_2  !allow /0 to show bugs
   ELSE
      soil%heat_cap_lower_limit(:,:) = soil%css_vec(:,:) * soil%rhosoil_vec(:,:)
   END IF

   IF( (.NOT.cable_user%cable_runtime_coupled ) .and. (first_gw_hydro_call)) THEN

         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
   END IF

   gammzz_snow(:,:) = 0._r_2
   do i=1,mp
      gammzz_snow(i,1) = real((1. - ssnow%isflag(i)) * C%cgsnow * ssnow%snowd(i),r_2)
   end do

   !Start with wb and wbice.  Need wbliq, wmliq,wmice,wmtot
   !find the mass of ice and liq from the prognostic volumetric values
   do k=1,ms
      do i=1,mp
         ssnow%wbliq(i,k) = ssnow%wb(i,k) - den_rat*ssnow%wbice(i,k)  ! MMY      !liquid volume
         ssnow%wmice(i,k) = ssnow%wbice(i,k)*real(C%density_ice,r_2)*soil%zse_vec(i,k)!ice mass
         ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*real(C%density_liq,r_2)*soil%zse_vec(i,k) !liquid mass
         ssnow%wmtot(i,k) = ssnow%wmice(i,k) + ssnow%wmliq(i,k)                  !liq+ice mass

         ssnow%wblf(i,k)   = max(ssnow%wbliq(i,k)/soil%ssat_vec(i,k),0.01_r_2)
         ssnow%wbfice(i,k) = max(ssnow%wbice(i,k)/soil%ssat_vec(i,k),0._r_2)
         !do not pass with mpi
         ssnow%sucs_hys(i,k) = ssnow%hys_fac(i,k)*soil%sucs_vec(i,k)
         ssnow%smp_hys(i,k)  = min(ssnow%smp_hys(i,k),-ssnow%sucs_hys(i,k))
      end do
   end do

   IF( first_gw_hydro_call ) THEN
     do k=1,ms
       do i=1,mp
            ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),&
                               (1.0-soil%ssat_vec(i,k))*&
                               soil%css_vec(i,k) * soil%rhosoil_vec(i,k)  &
                & + ssnow%wbliq(i,k) * real(C%cswat*C%density_liq,r_2)  &
                & + ssnow%wbice(i,k) * real(C%csice*C%density_ice,r_2) )* &
                 soil%zse_vec(i,k) +   gammzz_snow(i,k)

       end do
     end do

   ENDIF  ! if(.NOT.cable_runtime_coupled) and first_gw_hydro_call

   ssnow%wbliq_old = ssnow%wbliq

   do i=1,mp
      !initial water in the soil column
      wbtot_ic(i) = sum(ssnow%wmtot(i,:),dim=1)+ssnow%GWwb(i)*soil%GWdz(i)*C%density_liq

      GWwb_ic(i) = ssnow%GWwb(i)

   end do

   !improve hiding, call single soilsnow subroutine to do all the
   !snow processes and thermal soil calculations

   CALL snow_processes_soil_thermal(dels,ssnow,soil,veg,canopy,met,bal)

   !leave here for now, could move into soilsnow as well
   CALL remove_transGW(dels, soil, ssnow, canopy, veg)        !transpiration loss per soil layer

   CALL  GWsoilfreeze(dels, soil, ssnow)

   ssnow%fwtop = canopy%precis/dels + ssnow%smelt/dels   !water from canopy and snowmelt [mm/s]

   CALL iterative_wtd (ssnow, soil, veg, .true. )

   CALL ovrlndflx (dels, ssnow, soil, veg, canopy,use_sli )         !surface runoff, incorporate ssnow%pudsto?

   ssnow%sinfil = ssnow%fwtop - canopy%segg  !canopy%fes/C%hl               !remove soil evap from throughfall

   CALL smoistgw (dels,ktau,ssnow,soil,veg,canopy)               !vertical soil moisture movement.

   IF (gw_params%BC_hysteresis)  &
             CALL swc_hyst_direction(soil,ssnow,veg)

   !call swc_smp_dsmpdw(soil,ssnow)
   if (gw_params%BC_hysteresis) then
      call brook_corey_hysteresis_swc_smp(soil,ssnow)
   elseif (gw_params%HC_SWC) then
      call hutson_cass_swc_smp(soil,ssnow)
   else
      call brook_corey_swc_smp(soil,ssnow)
   end if

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

   do i=1,mp
      ssnow%pudsto(i) = 0.0  !no puddle
      ssnow%smelt(i)  = ssnow%smelt(i)/dels    !change units to mm/s.  cable_driver then reverts back to mm
      ssnow%runoff(i) = (ssnow%rnof1(i) + ssnow%rnof2(i))!*dels  !cable_driver converts from mm/s to mm
                                                        !rnof1 and rnof2 are already in mm/s
      ! Set weighted soil/snow surface temperature
      ssnow%tss(i) =  (1-ssnow%isflag(i))*ssnow%tgg(i,1) + ssnow%isflag(i)*ssnow%tggsn(i,1)

      !total water mass at the end of the soilsnow_GW routine
      ssnow%wbtot(i)  = sum(ssnow%wbliq(i,:)*C%density_liq*soil%zse(:),dim=1) + &
                     sum(ssnow%wbice(i,:)*C%density_ice*soil%zse(:),dim=1) + &
                     ssnow%GWwb(i)*soil%GWdz(i)*C%density_liq

      !for debug water balance.  del_wbtot = fluxes = infiltration [though-evap] - trans - qhorz drainage
      del_wbtot(i)   = dels * (ssnow%sinfil(i) - ssnow%rnof2(i) - canopy%fevc(i) / C%hl)
      !set below to keep track of water imbalance within the GW module explicitly.  also must change cable_checks
      !ssnow%wbtot(i) = ssnow%wbtot(i)-wbtot_ic(i)

   end do

   first_gw_hydro_call=.false.

END SUBROUTINE soil_snow_gw

SUBROUTINE calc_equilibrium_water_content(ssnow,soil)
  !find layer mean soil moisture and potential at equilibrium with wtd

  IMPLICIT NONE

    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    !local variables
    REAL(r_2), dimension(mp)    :: zaq      !node depth of the aquifer
    REAL(r_2), dimension(ms)    :: dzmm     !layer thickness for single tile
    REAL(r_2), dimension(mp)    :: GWdzmm   !aquifer thickness at each tile
    REAL(r_2), dimension(mp)    :: GWzimm   !aquifer layer interface depth
    REAL(r_2), dimension(mp,0:ms)  :: zimm     !layer interface depth in mm
    REAL(r_2), dimension(mp,ms)    :: zmm      !node depths in mm
    REAL(r_2)                   :: tempi, temp0,voleq1,wbrat,&
                                   zi_smpc,tmp_const,voleq2
    real(r_2), dimension(mp,ms+1) :: ztop,zbot

    INTEGER :: k,i

    !make code cleaner define these here
    zimm(:,:) = 0._r_2
    do k=1,ms
       zimm(:,k) = zimm(:,k-1) + m2mm*soil%zse_vec(:,k)
       zmm(:,k)  = zimm(:,k-1) + 0.5_r_2*m2mm*soil%zse_vec(:,k)
    end do

    do i=1,mp
       GWzimm(i) = zimm(i,ms)+m2mm*soil%GWdz(i)
       zaq(i)    = zimm(i,ms) + 0.5_r_2*m2mm*soil%GWdz(i)
    end do

    soil%wbc_GW(:) = 0.0
    soil%wbc_vec(:,:) = 0.0
    soil%smpc_vec(:,:) = 1.0e+30
    soil%smpc_GW(:) = 1.0e+30

    if (.not.gw_params%BC_hysteresis) then !not needed but make code clear
       ssnow%watr_hys = soil%watr(:,:)
       ssnow%ssat_hys = soil%ssat_vec(:,:)
       ssnow%sucs_hys = soil%sucs_vec

       if (gw_params%HC_SWC) then
         soil%wbc_vec(:,:) = 2.0*soil%bch_vec(:,:)*soil%ssat_vec(:,:)/&
                                 (1.0+2.0*soil%bch_vec(:,:))
         soil%smpc_vec(:,:) = -soil%sucs_vec(:,:) * &
                               (2.0*soil%bch_vec(:,:)/&
                               ((1.0+2.0*soil%bch_vec(:,:))))**(-soil%bch_vec(:,:))
         soil%wbc_GW(:) = 2.0*soil%GWbch_vec(:)*soil%GWssat_vec(:)/&
                                 (1.0+2.0*soil%GWbch_vec(:))
         soil%smpc_GW(:) = -soil%GWsucs_vec(:) * &
                               (2.0*soil%GWbch_vec(:)/&
                               ((1.0+2.0*soil%GWbch_vec(:))))**(-soil%GWbch_vec(:))
       end if
  end if

    !!equilibrium water content
    do k=1,ms
       do i=1,mp

          if ((ssnow%wtd(i) .le. zimm(i,k-1))) then         !fully saturated

             ssnow%wbeq(i,k) = ssnow%ssat_hys(i,k)

          elseif ((ssnow%wtd(i) .le. zimm(i,k)) .and. &
                  (ssnow%wtd(i) .gt. zimm(i,k-1))) then

             tempi = 1._r_2
             temp0 = &
                   (((ssnow%sucs_hys(i,k)+ssnow%wtd(i)-zimm(i,k-1))/&
                      ssnow%sucs_hys(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
             voleq1 = -ssnow%sucs_hys(i,k)*(ssnow%ssat_hys(i,k)-ssnow%watr_hys(i,k))/&
                       (1._r_2-1._r_2/soil%bch_vec(i,k))/&
                       (ssnow%wtd(i)-zimm(i,k-1))*(tempi-temp0) + ssnow%watr_hys(i,k)
             ssnow%wbeq(i,k) = (voleq1*(ssnow%wtd(i)-zimm(i,k-1)) +&
                               (ssnow%ssat_hys(i,k))&
                               *(zimm(i,k)-ssnow%wtd(i)))/(zimm(i,k)-zimm(i,k-1))
          else

             tempi = (((ssnow%sucs_hys(i,k)+ssnow%wtd(i)-zimm(i,k))/&
                      ssnow%sucs_hys(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
             temp0 = (((ssnow%sucs_hys(i,k)+ssnow%wtd(i)-zimm(i,k-1))/&
                     ssnow%sucs_hys(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
             ssnow%wbeq(i,k) = -ssnow%sucs_hys(i,k)*(ssnow%ssat_hys(i,k)-ssnow%watr_hys(i,k))/&
                               (1._r_2-1._r_2/soil%bch_vec(i,k))/&
                               (zimm(i,k)-zimm(i,k-1))*(tempi-temp0)+ssnow%watr_hys(i,k)
          end if

         if (.not.gw_params%HC_SWC) then
             ssnow%wbeq(i,k) = min(max(ssnow%wbeq(i,k),ssnow%watr_hys(i,k)),ssnow%ssat_hys(i,k))

             wbrat = min(max((&
                     ssnow%wbeq(i,k) - ssnow%watr_hys(i,k))/&
                             (ssnow%ssat_hys(i,k)-ssnow%watr_hys(i,k)),&
                          0.01_r_2),1._r_2)
              ssnow%zq(i,k) = max(&
                            -ssnow%sucs_hys(i,k)*(wbrat**(-soil%bch_vec(i,k))),sucmin)

          else
             if (ssnow%wbeq(i,k) .lt. soil%wbc_vec(i,k)) then
               wbrat = min(max((&
                       ssnow%wbeq(i,k) - soil%watr(i,k))/(soil%ssat_vec(i,k)-soil%watr(i,k)),&
                               0.01_r_2),1._r_2)

               ssnow%zq(i,k) = max(&
                                 -soil%sucs_vec(i,k)*(wbrat**(-soil%bch_vec(i,k))),sucmin)
             else
                ssnow%zq(i,k) = -soil%sucs_vec(i,k)* &
                             sqrt(1._r_2 - ssnow%wbeq(i,k)/soil%ssat_vec(i,k))/&
                                  (sqrt(1._r_2-soil%wbc_vec(i,k)/soil%ssat_vec(i,k))*&
                                 ( (soil%wbc_vec(i,k)/soil%ssat_vec(i,k))**soil%bch_vec(i,k) ))

             end if

          end if
       end do  !mp
    end do  !ms
    do i=1,mp
    !Aquifer Equilibrium water content
       if (ssnow%wtd(i) .le. zimm(i,ms)) then      !fully saturated

          ssnow%GWwbeq(i) = soil%GWssat_vec(i)

       elseif ((ssnow%wtd(i) .gt. GWzimm(i)))   then     !fully unsaturated

          tempi = &
                (((soil%GWsucs_vec(i)+ssnow%wtd(i)-GWzimm(i))/&
                  soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
          temp0 = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(i,ms))/&
                     soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
          ssnow%GWwbeq(i) = -soil%GWsucs_vec(i)*soil%GWssat_vec(i)/&
                          (1._r_2-1._r_2/soil%GWbch_vec(i))/&
                           (GWzimm(i)-zimm(i,ms))*(tempi-temp0) + soil%GWwatr(i)

       else

             tempi = 1._r_2
             temp0 = &
                   (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(i,ms))/&
                      soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
             voleq1 = -soil%GWsucs_vec(i)*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
                       (1._r_2-1._r_2/soil%GWbch_vec(i))/&
                       (ssnow%wtd(i)-zimm(i,ms))*(tempi-temp0) + soil%GWwatr(i)
             ssnow%GWwbeq(i) = (voleq1*(ssnow%wtd(i)-zimm(i,ms)) +&
                               (soil%GWssat_vec(i))&
                               *(GWzimm(i)-ssnow%wtd(i)))/(GWzimm(i)-zimm(i,ms))

       end if

       if (.not.gw_params%HC_SWC) then

          ssnow%GWwbeq(i) = min(max(ssnow%GWwbeq(i),soil%GWwatr(i)),soil%GWssat_vec(i))

          ssnow%GWzq(i) = -soil%GWsucs_vec(i)*(max((ssnow%GWwbeq(i)-soil%GWwatr(i))/     &
                       (soil%GWssat_vec(i)-soil%GWwatr(i)),0.01_r_2))**(-soil%GWbch_vec(i))
          ssnow%GWzq(i) = max(sucmin, ssnow%GWzq(i))

       else
          if (ssnow%GWwbeq(i) .lt. soil%wbc_GW(i)) then
            wbrat = min(max((&
                    ssnow%wbeq(i,k) - soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i)),&
                            0.01_r_2),1._r_2)

            ssnow%GWzq(i) = max(&
                              -soil%GWsucs_vec(i)*(wbrat**(-soil%GWbch_vec(i))),sucmin)
          else
             ssnow%gwzq(i) = -soil%gwsucs_vec(i)* &
                          sqrt(1._r_2 - ssnow%gwwbeq(i)/soil%gwssat_vec(i))/&
                               (sqrt(1._r_2-soil%wbc_gw(i)/soil%gwssat_vec(i))*&
                              ( (soil%wbc_gw(i)/soil%GWssat_vec(i))**soil%GWbch_vec(i) ))

          end if
        end if

    end do
    !equilibrium water content
!   else! (gw_params%HC_SWC) then
!
!   do k=1,ms
!      do i=1,mp
!         ztop(i,k) = -soil%sucs_vec(i,k)-ssnow%wtd(i) + zimm(i,k-1)
!         zbot(i,k) = -soil%sucs_vec(i,k)-ssnow%wtd(i) + zimm(i,k)
!      end do
!   end do
!   do i=1,mp
!      ztop(i,ms+1) = -soil%GWsucs_vec(i)-ssnow%wtd(i) + zimm(i,ms)
!      zbot(i,ms+1) = -soil%GWsucs_vec(i)-ssnow%wtd(i) + GWzimm(i)
!   end do
!
!    do k=1,ms
!       do i=1,mp
!
!          if ((ssnow%wtd(i) .le. zimm(i,k-1))) then         !fully saturated
!
!             ssnow%wbeq(i,k) = soil%ssat_vec(i,k)
!
!          elseif ((ssnow%wtd(i) .le. zimm(i,k)) .and. &
!                  (ssnow%wtd(i) .gt. zimm(i,k-1))) then
!
!
!             if (zbot(i,k) .lt. soil%smpc_vec(i,k) .and. &
!                 ztop(i,k) .lt. soil%smpc_vec(i,k) ) then
!
!             tempi = 1._r_2
!             temp0 = &
!                   (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(i,k-1))/&
!                      soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
!             voleq1 = -soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
!                       (1._r_2-1._r_2/soil%bch_vec(i,k))/&
!                       (ssnow%wtd(i)-zimm(i,k-1))*(tempi-temp0) + soil%watr(i,k)
!             ssnow%wbeq(i,k) = (voleq1*(ssnow%wtd(i)-zimm(i,k-1)) +&
!                               (soil%ssat_vec(i,k))&
!                               *(zimm(i,k)-ssnow%wtd(i)))/(zimm(i,k)-zimm(i,k-1))
!             elseif (zbot(i,k) .ge. soil%smpc_vec(i,k) .and. &
!                     ztop(i,k) .lt. soil%smpc_vec(i,k) ) then
!
!             zi_smpc = ssnow%wtd(i)-(soil%sucs_vec(i,k)+soil%smpc_vec(i,k))
!
!             tmp_const = (1.0/3.0)*soil%ssat_vec(i,k)*(1.0-soil%wbc_vec(i,k)/soil%ssat_vec(i,k))/&
!                                    (soil%sucs_vec(i,k)*soil%sucs_vec(i,k)*&
!                                     (soil%wbc_vec(i,k)/soil%ssat_vec(i,k))**(-soil%bch_vec(i,k)))
!
!             voleq2 = soil%ssat_vec(i,k)*(ssnow%wtd(i)-zi_smpc) + &
!                      tmp_const*( (-soil%sucs_vec(i,k))**3.0 - &
!                         (-soil%sucs_vec(i,k)+ssnow%wtd(i)-zi_smpc)**3.0)
!             tempi = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zi_smpc)/&
!                      soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
!             temp0 = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(i,k-1))/&
!                     soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
!             voleq1 =-soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
!                               (1._r_2-1._r_2/soil%bch_vec(i,k))/&
!                               (zi_smpc-zimm(i,k-1))*(tempi-temp0) + soil%watr(i,k)
!             ssnow%wbeq(i,k) =( (zimm(i,k)-ssnow%wtd(i))*soil%ssat_vec(i,k)+&
!                                (ssnow%wtd(i)-zi_smpc)*voleq2 + &
!                                (zi_smpc-zimm(i,k-1))*voleq1 ) / &
!                                 (zimm(i,k)-zimm(i,k-1))
!             else
!
!             tmp_const = (1.0/3.0)*soil%ssat_vec(i,k)*(1.0-soil%wbc_vec(i,k)/soil%ssat_vec(i,k))/&
!                                    (soil%sucs_vec(i,k)*soil%sucs_vec(i,k)*&
!                                     (soil%wbc_vec(i,k)/soil%ssat_vec(i,k))**(-soil%bch_vec(i,k)))
!
!             zi_smpc = zimm(i,k-1)
!             voleq2 = soil%ssat_vec(i,k)*(ssnow%wtd(i)-zi_smpc) + &
!                      tmp_const*( (-soil%sucs_vec(i,k))**3.0 - &
!                         (-soil%sucs_vec(i,k)+ssnow%wtd(i)-zi_smpc)**3.0)
!
!             ssnow%wbeq(i,k) =( (zimm(i,k)-ssnow%wtd(i))*soil%ssat_vec(i,k)+&
!                                 (ssnow%wtd(i)-zi_smpc)*voleq2) / (zimm(i,k)-zimm(i,k-1))
!
!             end if
!
!          else
!             if (zbot(i,k) .lt. soil%smpc_vec(i,k) .and.  &
!                 ztop(i,k) .lt. soil%smpc_vec(i,k) ) then
!
!             tempi = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(i,k))/&
!                      soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
!             temp0 = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(i,k-1))/&
!                     soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
!             ssnow%wbeq(i,k) = -soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
!                               (1._r_2-1._r_2/soil%bch_vec(i,k))/&
!                               (zimm(i,k)-zimm(i,k-1))*(tempi-temp0)+soil%watr(i,k)
!
!             elseif (zbot(i,k) .lt. soil%smpc_vec(i,k) .and.&
!                     ztop(i,k) .ge. soil%smpc_vec(i,k) ) then
!               !watr is zero for HC SWC
!             zi_smpc = zimm(i,k) - (soil%sucs_vec(i,k)+soil%smpc_vec(i,k))
!
!             tmp_const = (1.0/3.0)*soil%ssat_vec(i,k)*(1.0-soil%wbc_vec(i,k)/soil%ssat_vec(i,k))/&
!                                    (soil%sucs_vec(i,k)*soil%sucs_vec(i,k)*&
!                                     (soil%wbc_vec(i,k)/soil%ssat_vec(i,k))**(-soil%bch_vec(i,k)))
!
!             voleq2 = soil%ssat_vec(i,k)*(zimm(i,k) - zi_smpc) + &
!                      tmp_const*( (-soil%sucs_vec(i,k))**3.0 - &
!                         (-soil%sucs_vec(i,k)+ssnow%wtd(i)-zi_smpc)**3.0)
!             tempi = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zi_smpc)/&
!                      soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
!             temp0 = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(i,k-1))/&
!                     soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
!             voleq1 =-soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
!                               (1._r_2-1._r_2/soil%bch_vec(i,k))/&
!                               (zi_smpc-zimm(i,k-1))*(tempi-temp0) + soil%watr(i,k)
!             ssnow%wbeq(i,k) =( &
!                                (zimm(i,k)-zi_smpc)*voleq2 + &
!                                (zi_smpc-zimm(i,k-1))*voleq1 ) / &
!                                 (zimm(i,k)-zimm(i,k-1))
!             else
!             tmp_const = (1.0/3.0)*soil%ssat_vec(i,k)*(1.0-soil%wbc_vec(i,k)/soil%ssat_vec(i,k))/&
!                                    (soil%sucs_vec(i,k)*soil%sucs_vec(i,k)*&
!                                     (soil%wbc_vec(i,k)/soil%ssat_vec(i,k))**(-soil%bch_vec(i,k)))
!
!             ssnow%wbeq(i,k)  = (soil%ssat_vec(i,k)*(zimm(i,k) - zimm(i,k-1)) + &
!                               tmp_const*( (min(-sucmin,zbot(i,k)))**3.0 - &
!                                   (min(-sucmin,ztop(i,k)))**3.0) )/ (zimm(i,k)-zimm(i,k-1))
!
!
!             end if
!
!          end if
!
!          ssnow%wbeq(i,k) = min(max(ssnow%wbeq(i,k),soil%watr(i,k)),soil%ssat_vec(i,k))
!
!          if (ssnow%wbeq(i,k) .lt. soil%wbc_vec(i,k)) then
!            wbrat = min(max((&
!                    ssnow%wbeq(i,k) - soil%watr(i,k))/(soil%ssat_vec(i,k)-soil%watr(i,k)),&
!                            0.01_r_2),1._r_2)
!
!            ssnow%zq(i,k) = max(&
!                              -soil%sucs_vec(i,k)*(wbrat**(-soil%bch_vec(i,k))),sucmin)
!          else
!             ssnow%zq(i,k) = -soil%sucs_vec(i,k)* &
!                          sqrt(1._r_2 - ssnow%wbeq(i,k)/soil%ssat_vec(i,k))/&
!                               (sqrt(1._r_2-soil%wbc_vec(i,k)/soil%ssat_vec(i,k))*&
!                              ( (soil%wbc_vec(i,k)/soil%ssat_vec(i,k))**soil%bch_vec(i,k) ))
!
!          end if
!
!
!       end do  !mp
!    end do  !ms
!
!    do i=1,mp
!    !Aquifer Equilibrium water content
!       if (ssnow%wtd(i) .le. zimm(i,ms)) then      !fully saturated
!
!          ssnow%GWwbeq(i) = soil%GWssat_vec(i)
!
!       elseif ((ssnow%wtd(i) .gt. GWzimm(i)))   then     !fully unsaturated
!         if (zbot(i,ms+1) .lt. soil%smpc_GW(i) .and.  &
!             ztop(i,ms+1) .lt. soil%smpc_GW(i) ) then
!
!          tempi = &
!                (((soil%GWsucs_vec(i)+ssnow%wtd(i)-GWzimm(i))/&
!                  soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
!          temp0 = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(i,ms))/&
!                     soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
!          ssnow%GWwbeq(i) = -soil%GWsucs_vec(i)*soil%GWssat_vec(i)/&
!                          (1._r_2-1._r_2/soil%GWbch_vec(i))/&
!                           (GWzimm(i)-zimm(i,ms))*(tempi-temp0) + soil%GWwatr(i)
!
!         elseif (zbot(i,ms+1) .lt. soil%smpc_GW(i) .and.  &
!                ztop(i,ms+1) .ge. soil%smpc_GW(i) ) then
!               !watr is zero for HC SWC
!             zi_smpc = GWzimm(i) - (soil%GWsucs_vec(i)+soil%smpc_GW(i))
!
!             tmp_const = (1.0/3.0)*soil%GWssat_vec(i)*(1.0-soil%wbc_GW(i)/soil%GWssat_vec(i))/&
!                                    (soil%GWsucs_vec(i)*soil%GWsucs_vec(i)*&
!                                     (soil%wbc_GW(i)/soil%GWssat_vec(i))**(-soil%GWbch_vec(i)))
!
!             voleq2 = soil%GWssat_vec(i)*(GWzimm(i) - zi_smpc) + &
!                      tmp_const*( (-soil%GWsucs_vec(i))**3.0 - &
!                         (-soil%GWsucs_vec(i)+ssnow%wtd(i)-zi_smpc)**3.0)
!             tempi = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zi_smpc)/&
!                      soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
!             temp0 = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(i,ms))/&
!                     soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
!             voleq1 =-soil%GWsucs_vec(i)*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
!                               (1._r_2-1._r_2/soil%GWbch_vec(i))/&
!                               (zi_smpc-zimm(i,ms))*(tempi-temp0) + soil%GWwatr(i)
!             ssnow%GWwbeq(i) =( &
!                                (GWzimm(i)-zi_smpc)*voleq2 + &
!                                (zi_smpc-zimm(i,ms))*voleq1 ) / &
!                                 (zimm(i,k)-zimm(i,k-1))
!             else
!             tmp_const = (1.0/3.0)*soil%GWssat_vec(i)*(1.0-soil%wbc_GW(i)/soil%GWssat_vec(i))/&
!                                    (soil%GWsucs_vec(i)*soil%GWsucs_vec(i)*&
!                                     (soil%wbc_GW(i)/soil%GWssat_vec(i))**(-soil%GWbch_vec(i)))
!
!             ssnow%wbeq(i,k)  = (soil%GWssat_vec(i)*(GWzimm(i) - zimm(i,ms)) + &
!                               tmp_const*( (zbot(i,k))**3.0 - &
!                                   (ztop(i,k))**3.0) )/ (GWzimm(i)-zimm(i,ms))
!             end if
!       else
!             if (zbot(i,k) .lt. soil%smpc_GW(i) .and. &
!                 ztop(i,k) .lt. soil%smpc_GW(i) ) then
!
!             tempi = 1._r_2
!             temp0 = &
!                   (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(i,ms))/&
!                      soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
!             voleq1 = -soil%GWsucs_vec(i)*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
!                       (1._r_2-1._r_2/soil%GWbch_vec(i))/&
!                       (ssnow%wtd(i)-zimm(i,ms))*(tempi-temp0) + soil%GWwatr(i)
!             ssnow%GWwbeq(i) = (voleq1*(ssnow%wtd(i)-zimm(i,ms)) +&
!                               (soil%GWssat_vec(i))&
!                               *(GWzimm(i)-ssnow%wtd(i)))/(GWzimm(i)-zimm(i,ms))
!             elseif (zbot(i,k) .ge. soil%smpc_GW(i) .and. &
!                     ztop(i,k) .lt. soil%smpc_GW(i) ) then
!
!             zi_smpc = ssnow%wtd(i)-(soil%GWsucs_vec(i)+soil%smpc_GW(i))
!
!             tmp_const = (1.0/3.0)*soil%GWssat_vec(i)*(1.0-soil%wbc_GW(i)/soil%GWssat_vec(i))/&
!                                    (soil%GWsucs_vec(i)*soil%GWsucs_vec(i)*&
!                                     (soil%wbc_GW(i)/soil%GWssat_vec(i))**(-soil%GWbch_vec(i)))
!
!             voleq2 = soil%GWssat_vec(i)*(ssnow%wtd(i)-zi_smpc) + &
!                      tmp_const*( (-soil%GWsucs_vec(i))**3.0 - &
!                         (-soil%GWsucs_vec(i)+ssnow%wtd(i)-zi_smpc)**3.0)
!             tempi = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zi_smpc)/&
!                      soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
!             temp0 = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(i,ms))/&
!                     soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
!             voleq1 =-soil%GWsucs_vec(i)*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
!                               (1._r_2-1._r_2/soil%GWbch_vec(i))/&
!                               (zi_smpc-zimm(i,ms))*(tempi-temp0) + soil%GWwatr(i)
!             ssnow%GWwbeq(i) =( (GWzimm(i)-ssnow%wtd(i))*soil%GWssat_vec(i)+&
!                                (ssnow%wtd(i)-zi_smpc)*voleq2 + &
!                                (zi_smpc-zimm(i,ms))*voleq1 ) / &
!                                 (GWzimm(i)-zimm(i,ms))
!             else
!
!             tmp_const = (1.0/3.0)*soil%GWssat_vec(i)*(1.0-soil%wbc_GW(i)/soil%GWssat_vec(i))/&
!                                    (soil%GWsucs_vec(i)*soil%GWsucs_vec(i)*&
!                                     (soil%wbc_GW(i)/soil%GWssat_vec(i))**(-soil%GWbch_vec(i)))
!
!             zi_smpc = zimm(i,ms)
!             voleq2 = soil%GWssat_vec(i)*(ssnow%wtd(i)-zi_smpc) + &
!                      tmp_const*( (-soil%GWsucs_vec(i))**3.0 - &
!                         (-soil%GWsucs_vec(i)+ssnow%wtd(i)-zi_smpc)**3.0)
!
!             ssnow%GWwbeq(i) =( (GWzimm(i)-ssnow%wtd(i))*soil%GWssat_vec(i)+&
!                                 (ssnow%wtd(i)-zi_smpc)*voleq2) / (GWzimm(i)-zimm(i,ms))
!
!             end if
!
!
!       end if
!
!       ssnow%GWwbeq(i) = min(max(ssnow%GWwbeq(i),soil%GWwatr(i)),soil%GWssat_vec(i))
!
!       if (ssnow%GWwbeq(i) .lt. soil%wbc_GW(i)) then
!          ssnow%GWzq(i) = -soil%GWsucs_vec(i)*(max((ssnow%GWwbeq(i)-soil%GWwatr(i))/     &
!                       (soil%GWssat_vec(i)-soil%GWwatr(i)),0.001_r_2))**(-soil%GWbch_vec(i))
!       else
!          ssnow%GWzq(i) = -soil%GWsucs_vec(i)* &
!                           sqrt(1._r_2 - ssnow%GWwbeq(i)/soil%GWssat_vec(i))/&
!                            (sqrt(1._r_2-soil%wbc_GW(i)/soil%GWssat_vec(i))*&
!                           ( (soil%wbc_GW(i)/soil%GWssat_vec(i))**soil%GWbch_vec(i) ))
!       end if
!       ssnow%GWzq(i) = max(sucmin, ssnow%GWzq(i))
!
!    end do
!
!    end if

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

    !CALL point2constants( C )


    IF (cable_user%or_evap) THEN

       call saturated_fraction(ssnow,soil,veg)

       ssnow%wetfac(:) = 1.0

      do i=1,mp
         IF( ssnow%snowd(i) > 0.1) ssnow%wetfac(i) = 0.9

         IF ( veg%iveg(i) == 16 .and. met%tk(i) >= C%TFRZ + 5. )   &
                 ssnow%wetfac(i) = 1.0 ! lakes: hard-wired number to be removed

         IF( veg%iveg(i) == 16 .and. met%tk(i) < C%TFRZ + 5. )   &
                 ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed
      end do

    ELSEIF (cable_user%gw_model) THEN

       call saturated_fraction(ssnow,soil,veg)

      do i = 1,mp
        ! MMY it doesn't match to icef in SUBROUTINE ovrlndflx. I made them consistent
         dzmm_one  = m2mm * real(soil%zse_vec(i,1),r_2)
         icemass  = ssnow%wmice(i,1) ! MMY ssnow%wbice(i,1) * dzmm_one
         liqmass  = ssnow%wmliq(i,1) ! MMY (ssnow%wb(i,1)-den_rat*ssnow%wbice(i,1)) * dzmm_one
         totmass  = max(liqmass+icemass,real(1e-2,r_2))
         icef(i)  = max(0._r_2,min(1._r_2,icemass / totmass)) ! MMY
               ! MMY max(0._r_2,min(1._r_2, gw_params%IceBeta*icemass / totmass))
     end do

      !srf frozen fraction.  should be based on topography
      do i = 1,mp
         fice = (exp(-gw_params%IceAlpha*(1._r_2-icef(i)))-&
                 exp(-gw_params%IceAlpha))/&
                 (1._r_2-exp(-gw_params%IceAlpha))
         fice = min(1._r_2,max(0._r_2,fice))

         satfrac_liqice(i) = max(0.,min(0.95,fice + (1._r_2-fice)*ssnow%satfrac(i) ) ) ! MMY
                             ! MMY fice + (1._r_2-fice)*ssnow%satfrac(i)

         wb_unsat = ((ssnow%wb(i,1)-ssnow%wbice(i,1)*den_rat) -& ! MMY add *den_rat
                     ssnow%satfrac(i)*soil%ssat_vec(i,1))/(1.-ssnow%satfrac(i))
         wb_unsat = min(soil%ssat_vec(i,1),max(0.,wb_unsat))

         wb_evap_threshold = min( max( &
                             gw_params%SoilEvapAlpha*soil%sfc_vec(i,1), &
                             soil%swilt_vec(i,1) ), soil%ssat_vec(i,1) )

         !Sakguchi and Zeng 2009
         if (wb_unsat .ge. wb_evap_threshold) then
            xx = 1.
         else
            xx = 0.25 * (1._r_2 - cos(3.14159_r_2*wb_unsat/(wb_evap_threshold)))**2.0
         end if

         ssnow%wetfac(i) = max(0.0,min(1.0,satfrac_liqice(i) +&
                                        (1. - satfrac_liqice(i))*xx ) )

      end do

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
                                real(MAX( 0.5_r_2, 1._r_2 - MIN( 0.2_r_2, &
                                ( ssnow%wbice(i,1) * den_rat & ! MMY not sure whether needed but add * den_rat
                                / ssnow%wb(i,1) )**2 ) ) )

          IF( ssnow%snowd(i) > 0.1) ssnow%wetfac(i) = 0.9

          IF ( veg%iveg(i) == 16 .and. met%tk(i) >= C%tfrz + 5. )   &
               ssnow%wetfac(i) = 1.0 ! lakes: hard-wired number to be removed

          IF( veg%iveg(i) == 16 .and. met%tk(i) < C%tfrz + 5. )   &
               ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed

       ENDDO
       ! owetfac introduced to reduce sharp changes in dry regions,
       ! especially in offline runs in which there may be discrepancies b/n
       ! timing of precip and temperature change (EAK apr2009)
       ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)

   ENDIF  !or_evap, gw_model, or default wetfac parameterization

END SUBROUTINE calc_srf_wet_fraction

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
   REAL(r_2), dimension(mp,ms+1) ::wb_temp
   REAL(r_2), DIMENSION(mp,ms+1) :: hk_ice_factor

   if (gw_params%BC_hysteresis) then
      !swc_smp_dsmpdw => brook_corey_hysteresis_swc_smp
      ssnow%sucs_hys(:,:) = ssnow%hys_fac(:,:)*soil%sucs_vec(:,:)
   elseif (gw_params%HC_SWC) then
      !swc_smp_dsmpdw => hutson_cass_swc_smp
      ssnow%sucs_hys(:,:) = soil%sucs_vec(:,:)
   else
      !swc_smp_dsmpdw => brook_corey_swc_smp
      ssnow%sucs_hys(:,:) = soil%sucs_vec(:,:)
   end if

   ! ___ MMY COMMENT: for Hutson Cass SWC potential ___
   do k=1,ms
      do i=1,mp
         soil%wbc_vec(i,k) = 2.0*soil%bch_vec(i,k)*soil%ssat_vec(i,k)/&
                                 (1.0+2.0*soil%bch_vec(i,k))
         soil%smpc_vec(i,k) = -soil%sucs_vec(i,k) * &
                               (2.0*soil%bch_vec(i,k)/&
                               ((1.0+2.0*soil%bch_vec(i,k))))**(-soil%bch_vec(i,k))
      end do
   end do
   do i=1,mp
      soil%wbc_GW(i) = 2.0*soil%GWbch_vec(i)*soil%GWssat_vec(i)/&
                                 (1.0+2.0*soil%GWbch_vec(i))
      soil%smpc_GW(i) = -soil%GWsucs_vec(i) * &
                               (2.0*soil%GWbch_vec(i)/&
                               ((1.0+2.0*soil%GWbch_vec(i))))**(-soil%GWbch_vec(i))
   end do
  ! ___________________________________________________
    do k=1,ms
       do i=1,mp
          ! _____ MMY it doesn't match to icef in SUBROUTINE ovrlndflx ______
          !ssnow%icefrac(i,k) = ssnow%wbice(i,k)/(max(ssnow%wb(i,k),0.01_r_2)) ! MMY
          ! ________________ MMY to unify the calculation ____________________
          ssnow%icefrac(i,k) = max(0._r_2, min(1._r_2, &
                               ssnow%wmice(i,k) / max(ssnow%wmtot(i,k), 0.01_r_2)) )
          ! __________________________________________________________________
          ssnow%fracice(i,k) = (exp(-gw_params%IceAlpha*(1._r_2-ssnow%icefrac(i,k)))&
                               -exp(-gw_params%IceAlpha))/(1._r_2-exp(-gw_params%IceAlpha))
       end do
    end do

    ssnow%fracice(:,:) = max( min( ssnow%fracice, 1._r_2), 0._r_2)

    if (gw_params%ssgw_ice_switch) then
       wb_temp(:,1:ms) =  ssnow%wbliq(:,:)
       wb_temp(:,ms+1) = ssnow%GWwb(:)
       do k=1,ms
          kk = min(k+1,ms)
          do i=1,mp
             if (soil%isoilm(i) .eq. 9) then
                !hk_ice_factor(i,k) = (1.0-soil%ssat_vec(i,k))**(gw_params%ice_impedence)
                hk_ice_factor(i,k) = 10.0**(gw_params%ice_impedence)
             else
                !hk_ice_factor(i,k) = sqrt((1.0-ssnow%wbice(i,k))**(gw_params%ice_impedence) *&
                !                     (1.0-ssnow%wbice(i,kk))**(gw_params%ice_impedence))
                hk_ice_factor(i,k) = min(10.0**(gw_params%ice_impedence*ssnow%wbice(i,k)/(soil%ssat_vec(i,k)-soil%watr(i,k))) ,&
                                          10.0**(gw_params%ice_impedence*ssnow%wbice(i,kk)/(soil%ssat_vec(i,kk)-soil%watr(i,kk))))
             end if
          end do
       end do
       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
    else
       wb_temp(:,1:ms) = ssnow%wb(:,:)
       wb_temp(:,ms+1) = ssnow%GWwb(:)
       do k=1,ms
          kk = min(k+1,ms)
          do i=1,mp
             hk_ice_factor(i,k) = (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))
          end do
       end do
       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
    end if

    liq_ratio(:) = 1.0
    k = ms
    do i=1,mp
       liq_ratio(i) =min(1.,max(0.,wb_temp(i,k)/max(ssnow%wb(i,k),1e-6) ) )
    end do

    !aquifer ice ! MMY should be aquifer liq
    wb_temp(:,ms+1) = liq_ratio(:) * wb_temp(:,ms+1)
    ! MMY for ssgw_ice_switch true, wb_temp(:,ms+1) = ssnow%GWwb(:) * wb_temp(:,ms)/max(ssnow%wb(:,ms))
    ! MMY for ssgw_ice_switch false, wb_temp(:,ms+1) = ssnow%GWwb(:)

    !soil potential (head) calculations can use brooks-corey
    !or hutson-cass SWC

    !call swc_smp_dsmpdw(soil,ssnow)
    if (gw_params%BC_hysteresis) then
       call brook_corey_hysteresis_swc_smp(soil,ssnow)
    elseif (gw_params%HC_SWC) then
       call hutson_cass_swc_smp(soil,ssnow)
    else
       call brook_corey_swc_smp(soil,ssnow)
    end if


    !hydraulic conductivity
    !Interfacial so uses layer i and i+1
    do k=1,ms-1
       ! kk=min(ms+1,k+1) ! MMY
       kk= k+1 ! MMY
       do i=1,mp
           ! ____________________ MMY ______________________
           s1(i) = 0.5_r_2*((wb_temp(i,k)-soil%watr(i,k)) + &
                            (wb_temp(i,kk)-soil%watr(i,kk))) / &
                          (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + &
                          (soil%ssat_vec(i,kk)-soil%watr(i,kk))))
           s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
           s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)
          ! _________________________________________________
          ssnow%hk(i,k)    =  s1(i)*s2(i)*hk_ice_factor(i,k)
          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*s2(i)*&
                            0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))*&
                            hk_ice_factor(i,k)
          ! MMY Note that the dhkdw equation doesn't exactly follow the finite difference
       end do
    end do

    !conductivity between soil column and aquifer
    k = ms
    kk = ms+1
       do i=1,mp

          s1(i) = 0.5_r_2*(max(wb_temp(i,k )-soil%watr(i,k),0.) + &
                           max(wb_temp(i,kk)-soil%GWwatr(i),0.)) / &
                  (0.5_r_2*(soil%ssat_vec(i,k)-soil%watr(i,k) +&
                            soil%GWssat_vec(i)-soil%GWwatr(i)))

          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)

          ssnow%hk(i,k)    = s1(i)*s2(i)*hk_ice_factor(i,k)
          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*&
                             hk_ice_factor(i,k)*&
                             s2(i)*0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))

           !Aquifer

          s2(i) = soil%GWhyds_vec(i)*s1(i)**(2._r_2*soil%GWbch_vec(i)+2._r_2)
          ssnow%GWhk(i)     =s1(i)*s2(i) * hk_ice_factor(i,ms+1)
          ssnow%GWdhkdw(i)  =  (2._r_2*soil%GWbch_vec(i)+3._r_2)*&
                              s2(i)*0.5_r_2/(soil%GWssat_vec(i)-soil%GWwatr(i)) *&
                              hk_ice_factor(i,ms+1)
       end do


END SUBROUTINE calc_soil_hydraulic_props


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
!   REAL(r_2), dimension(mp,ms+1) ::wb_temp
!   REAL(r_2), DIMENSION(mp,ms+1) :: hk_ice_factor
!
!
!    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water (calculated using total (not liquid))
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
!       wb_temp(:,1:ms) =  ssnow%wbliq(:,:)
!       wb_temp(:,ms+1) = ssnow%GWwb(:)
!       do k=1,ms
!          kk = min(k+1,ms)
!          do i=1,mp
!             if (soil%isoilm(i) .eq. 9) then
!                hk_ice_factor(i,k) = (1.0-soil%ssat_vec(i,k))**(gw_params%ice_impedence)
!             else
!                hk_ice_factor(i,k) = sqrt((1.0-ssnow%wbice(i,k))**(gw_params%ice_impedence) *&
!                                     (1.0-ssnow%wbice(i,kk))**(gw_params%ice_impedence))
!             end if
!          end do
!       end do
!       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
!    else
!       wb_temp(:,1:ms) = ssnow%wb(:,:)
!       wb_temp(:,ms+1) = ssnow%GWwb(:)
!       do k=1,ms
!          kk = min(k+1,ms)
!          do i=1,mp
!             hk_ice_factor(i,k) = (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))
!          end do
!       end do
!       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
!    end if
!
!    liq_ratio(:) = 1.0
!    k = ms
!    do i=1,mp
!       liq_ratio(i) =min(1.,max(0.,wb_temp(i,k)/max(ssnow%wb(i,k),1e-6) ) )
!    end do
!    !aquifer ice
!    wb_temp(:,ms+1) = liq_ratio(:) * wb_temp(:,ms+1)
!
!    !potential from soil water rention function
!    !defined as layer average
!    do k=1,ms
!       do i=1,mp
!          s_mid(i) = (wb_temp(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
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
!    !Aquifer potential
!    do i=1,mp
!       s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/&
!                    (soil%GWssat_vec(i)-soil%GWwatr(i))
!       s_mid(i) = min(max(s_mid(i),0.01_r_2),1._r_2)
!       s2(i)    = soil%GWhyds_vec(i)*s_mid(i)**(2._r_2*soil%GWbch_vec(i)+2._r_2)
!
!       ssnow%GWhk(i)     =s_mid(i)*s2(i) * hk_ice_factor(i,ms+1)
!       ssnow%GWdhkdw(i)  =  (2._r_2*soil%GWbch_vec(i)+3._r_2)*&
!                           s2(i)*0.5_r_2/(soil%GWssat_vec(i)-soil%GWwatr(i)) *&
!                           hk_ice_factor(i,ms+1)
!
!
!       s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
!       s_mid(i) = min(max(s_mid(i),0.01_r_2),1._r_2)
!
!       ssnow%GWsmp(i)    = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))
!       ssnow%GWsmp(i)    = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)
!       ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/&
!                            (s_mid(i)*(soil%GWssat_vec(i)-soil%GWwatr(i)))
!    end do
!
!    !hydraulic conductivity
!    !Interfacial so uses layer i and i+1
!    do k=1,ms
!       kk=min(ms+1,k+1)
!       do i=1,mp
!
!          if (k .lt. ms) then
!          s1(i) = 0.5_r_2*((wb_temp(i,k)-soil%watr(i,k)) + &
!                           (wb_temp(i,kk)-soil%watr(i,kk))) / &
!                         (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + &
!                         (soil%ssat_vec(i,kk)-soil%watr(i,kk))))
!          else
!          s1(i) = 0.5_r_2*((wb_temp(i,k)-soil%watr(i,k)) + &
!                           (wb_temp(i,kk)-soil%GWwatr(i))) / &
!                         (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + &
!                         (soil%GWssat_vec(i)-soil%GWwatr(i))))
!          end if
!          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
!          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)
!
!          ssnow%hk(i,k)    =  s1(i)*s2(i)*hk_ice_factor(i,k)
!          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*s2(i)*&
!                            0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))*&
!                            hk_ice_factor(i,k)
!       end do
!    end do
!
!    !conductivity between soil column and aquifer
!    k = ms
!    kk = ms+1
!       do i=1,mp
!
!          s1(i) = 0.5_r_2*(max(wb_temp(i,k )-soil%watr(i,k),0.) + &
!                           max(wb_temp(i,kk)-soil%GWwatr(i),0.)) / &
!                  (0.5_r_2*(soil%ssat_vec(i,k)-soil%watr(i,k) +&
!                            soil%GWssat_vec(i)-soil%GWwatr(i)))
!
!          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
!          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)
!
!          ssnow%hk(i,k)    = s1(i)*s2(i)*hk_ice_factor(i,k)
!          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*&
!                             hk_ice_factor(i,k)*&
!                             s2(i)*0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))
!       end do
!
!
!END SUBROUTINE calc_soil_hydraulic_props


  SUBROUTINE aquifer_recharge(dt,ssnow,soil,veg)
  USE cable_common_module

  IMPLICIT NONE
    real, intent(in) :: dt
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    REAL(r_2), dimension(mp)   :: zaq
    REAL(r_2), dimension(mp,ms) :: zmm,csum_dzmm

    integer :: i,k

    csum_dzmm(:,1) = m2mm*soil%zse_vec(:,1)
    do k=2,ms
       csum_dzmm(:,k) = csum_dzmm(:,k-1) + m2mm*soil%zse_vec(:,k)
    end do
    !Doing the recharge outside of the soln of Richards Equation makes it easier to track total recharge amount.
    !Add to ssnow at some point
    select case (gw_params%aquifer_recharge_function)
       case(0)
          ssnow%Qrecharge(:) = 0._r_2
       case(1)
          do i=1,mp
             if ((ssnow%wtd(i) .le. csum_dzmm(i,ms-1)) .or. &
                 (veg%iveg(i) .ge. 16) .or. &
                 (soil%isoilm(i) .eq. 9))  then

                ssnow%Qrecharge(i) = 0._r_2
             else
                ssnow%Qrecharge(i) = -0.5*(ssnow%hk(i,ms)+ssnow%GWhk(i))*&
                                     ((-ssnow%smp(i,ms)) -&
                                      (-ssnow%zq(i,ms))) / &
                                      (ssnow%wtd(i) - &
                                       (csum_dzmm(i,ms)-0.5*soil%zse_vec(i,ms)*m2mm))
             end if
          end do

       case default
          do i=1,mp
             if ((ssnow%wtd(i) .le. csum_dzmm(i,ms)) .or. &
                 (veg%iveg(i) .ge. 16) .or. &
                 (soil%isoilm(i) .eq. 9))  then

                ssnow%Qrecharge(i) = 0._r_2
             else
                ssnow%Qrecharge(i) = -(ssnow%hk(i,ms)+ssnow%GWhk(i))*&
                                     ((ssnow%GWsmp(i)-ssnow%smp(i,ms)) -&
                                      (ssnow%GWzq(i)-ssnow%zq(i,ms))) / &
                                     (m2mm*(soil%GWdz(i)+soil%zse_vec(i,ms)))
             end if
          end do
    end select

  END SUBROUTINE aquifer_recharge

  SUBROUTINE subsurface_drainage(ssnow,soil,veg)
  USE cable_common_module

  IMPLICIT NONE

    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    REAL(r_2), dimension(mp)                  :: sm_tot,&!sum var
                                                 ice_factor_tot  !avg ice factor
    REAL(r_2), dimension(mp,ms)               :: ice_factor,&  !ice limitation on
                                                 dzmm !subsurface drainage
    INTEGER, dimension(mp)                    :: k_drain
    integer :: i,k

    real(r_2), dimension(17) :: Efold_mod

    Efold_mod(:) = 1.0
    !Efold_mod(1:4) = (/0.2,0.2,0.2,0.2/)
    !Efold_mod(9) = 0.25

    dzmm(:,:) = m2mm*soil%zse_vec(:,:)

    do i=1,mp

       !Note: future revision will have interaction with river here. nned to
       !work on router and add river type cells
       ssnow%qhz(i)  = soil%slope(i)*soil%qhz_max(i)*&
                        exp( -mm2m*ssnow%wtd(i)/ soil%qhz_efold(i) )
       ! MMY have no idea where the equation comes from, doubtful

       !drain from sat layers
       k_drain(i) = ms+1
       do k=ms,2,-1
          if (ssnow%wtd(i) .le. sum(dzmm(i,1:k),dim=1)) then
             k_drain(i) = k + 1
          end if
       end do
       k_drain(i) = max(k_drain(i),3)

   end do

   if (gw_params%ssgw_ice_switch) then
      do k=1,ms
         do i=1,mp
            ice_factor(i,k) = (1.0-ssnow%wbice(i,k))**gw_params%ice_impedence
         end do
      end do
   else
      do k=1,ms
         do i=1,mp
            ice_factor(i,k) = (1._r_2-ssnow%fracice(i,k))
         end do
      end do

   end if

   ssnow%qhlev(:,:) = 0._r_2
   sm_tot(:)        = 1._r_2

   do i=1,mp

       if (gw_params%subsurface_sat_drainage) then
          ! MMY it doesn't make sense since sm_tot == 0.
          ! sm_tot(i) = sum(ssnow%qhlev(i,k_drain(i):ms+1),dim=1) ! MMY
          ! ______ MMY find this calculation from version : cable-2.2.3-pore-scale-model _______________
          sm_tot(i) = max(ssnow%GWwb(i)-soil%GWwatr(i),0._r_2)
          do k=k_drain(i),ms
            sm_tot(i) = sm_tot(i) + max(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)!*dzmm(k)
          end do
          ! ____________________________________________________________________________________________
          if (sm_tot(i) .lt. 1.0e-8) then
             sm_tot(i) = 1._r_2
          end if
       end if
   end do


   do i=1,mp

       ssnow%qhlev(i,ms+1) = max((ssnow%GWwb(i) - soil%GWwatr(i)),0._r_2)*&
                       ice_factor(i,ms)*ssnow%qhz(i)/sm_tot(i)

       do k=k_drain(i),ms
          ssnow%qhlev(i,k) = max(ssnow%wbliq(i,k)-ssnow%watr_hys(i,k),0._r_2)*&
                                   ice_factor(i,k)*ssnow%qhz(i)/sm_tot(i)
       end do

       !incase every layer is frozen very dry
       ssnow%qhz(i) = sum(ssnow%qhlev(i,:),dim=1)

       !Keep "lakes" saturated forcing qhz = 0.  runoff only from lakes
       !overflowing
       if (soil%isoilm(i) .eq. 9 .or. veg%iveg(i) .ge. 16) then
          ssnow%qhz(i) = 0._r_2
          ssnow%qhlev(i,:) = 0._r_2
       end if

    end do


  END SUBROUTINE subsurface_drainage


  SUBROUTINE saturated_fraction(ssnow,soil,veg)
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE(veg_parameter_type) , INTENT(IN)    :: veg  ! veg parameters

    REAL(r_2), DIMENSION(mp) :: S
    REAL(r_2) :: slopeSTDmm
    INTEGER :: i,k

     S(:) = 0._r_2
     do k=1,gw_params%level_for_satfrac
       S(:) = S(:) + max(0.01,min(1.0, &
              (ssnow%wb(:,k)-den_rat*ssnow%wbice(:,k)-ssnow%watr_hys(:,k))/&
               max(0.001,ssnow%ssat_hys(:,k)-soil%watr(:,k)-den_rat*ssnow%wbice(:,k)) ) )*soil%zse_vec(:,k)
     end do
     S(:) = S(:)/sum(soil%zse(1:gw_params%level_for_satfrac),dim=1)
     !srf frozen fraction.  should be based on topography
      do i = 1,mp
         !Saturated fraction
          if (gw_params%MaxSatFraction .gt. 1e-7 .and. veg%iveg(i) .lt. 16) then
             slopeSTDmm = sqrt(min(max(&
                           gw_params%MaxSatFraction*soil%slope_std(i),&
                           1e-5),10000.0)) ! ensure some variability
             ssnow%satfrac(i)    = max(0._r_2,min(0.99_r_2,&
                 !note UM wants std03, and erf is not included then
                                   1._r_2 - my_erf( slopeSTDmm * sqrt(2.0* S(i)) ) ) ) ! MMY change / to *@24 Sep 2021
          elseif (veg%iveg(i) .lt. 16) then
             ssnow%satfrac(i) = 0._r_2
          else
             ssnow%satfrac(i) = 0.975
          end if
      end do


  END SUBROUTINE saturated_fraction

  SUBROUTINE pore_space_relative_humidity(ssnow,soil,veg)
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE(veg_parameter_type), INTENT(INOUT)      :: veg

    REAL(r_2), DIMENSION(mp) :: unsat_wb,unsat_smp
    INTEGER :: i

    ! Need a matching array of ones to use in Mark's call to the intrinsic
    ! sign func below
    ! mgk, 24/07/2018
    REAL(r_2), DIMENSION(mp,ms) :: minus_ones

    !CALL point2constants( C )

    !if gw_model = true
                   !cable_um_init_subrs.F90 or cable_parameters:
                           ! ssat(i) = ssat_vec(i,1)
    !if gw_model = false
                   !cable_um_init_subrs.F90 or cable_parameters:
                            !ssat_vec(i,:) = ssat(i)
                            !so ssat_vec can be used although soilsnow uses ssat
    minus_ones = -1.

    do i=1,mp
       if (veg%iveg(i) .lt. 16 .and. soil%isoilm(i) .ne. 9 .and. &
           ssnow%snowd(i) .le. 1e-8 ) then

          unsat_wb(i) =  (ssnow%wb(i,1) - ssnow%ssat_hys(i,1)*&
                      min(0.95,max(0.0,ssnow%satfrac(i))))/(1.0 - min(0.95,max(0.0,ssnow%satfrac(i))))
          ! MMY ??? I'm not sure whether we should take wbice*den_rat into consideration here,
          !     ??? so I didn't modify.
          unsat_wb(i) = max(ssnow%watr_hys(i,1)+1e-2, min(ssnow%ssat_hys(i,1), unsat_wb(i) ) )

          !unsat_smp(i) = sign(soil%sucs_vec(i,1),-1.0) * min(1.0, &
          !                         (max(0.001, (unsat_wb(i)-soil%watr(i,1))/(soil%ssat_vec(i,1)-&
          !                         soil%watr(i,1)) ) )** (-soil%bch_vec(i,1)) )

          ! mgk, 24/07/2018 - fix to compile
          unsat_smp(i) = SIGN(soil%sucs_vec(i,1),minus_ones(i,1)) * MIN(1.0, &
               (MAX(0.001, (unsat_wb(i)-soil%watr(i,1))/(soil%ssat_vec(i,1)-&
               soil%watr(i,1)) ) )** (-soil%bch_vec(i,1)) )

          unsat_smp(i) = MAX(-1.0e7,unsat_smp(i) )/1000._r_2 !m

          ssnow%rh_srf(i) = max(0.,min(1., &
                         exp(9.81*unsat_smp(i)/(ssnow%tgg(i,1)*461.4)) ) )

       else

          ssnow%rh_srf(i) = 1.0

       end if
    end do


  END SUBROUTINE pore_space_relative_humidity

  SUBROUTINE sli_hydrology(dels,ssnow,soil,veg,canopy)
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil ! soil parameters
    TYPE(veg_parameter_type) , INTENT(INOUT)    :: veg  ! veg parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy

    LOGICAL, SAVE :: sli_call = .true.

    REAL(r_2), DIMENSION(ms) :: dzmm
    REAL(r_2), DIMENSION(mp) :: zmm
    REAL(r_2), DIMENSION(mp) :: zaq

    !CALL point2constants( C )

    ssnow%sucs_hys = soil%sucs_vec
    ssnow%ssat_hys = soil%ssat_vec
    ssnow%watr_hys = soil%watr

    call iterative_wtd (ssnow, soil, veg, cable_user%test_new_gw)

    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    call  ovrlndflx (dels, ssnow, soil, veg,canopy,sli_call )

    CALL subsurface_drainage(ssnow,soil,veg)

    call aquifer_recharge(dels,ssnow,soil,veg)




  END SUBROUTINE sli_hydrology


  SUBROUTINE set_unsed_gw_vars(ssnow,soil,canopy)
    ! MMY Actually, this subrountine is called by nowhere
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

  real(r_2) function my_erf(x)

  implicit none

  real(r_2), intent(in) :: x
  real(r_2)             :: tmp_val, tmp


   tmp = 1.0 / ( 1.0 + 0.5 * abs(x) )

   tmp_val =  tmp * exp(-abs(x) * abs(x) - 1.26551223 + tmp *     &
             ( 1.00002368 + tmp * ( 0.37409196 + tmp *          &
         ( 0.09678418 + tmp * (-0.18628806 + tmp *              &
                     ( 0.27886807 + tmp * (-1.13520398 + tmp *          &
         ( 1.48851587 + tmp * (-0.82215223 + tmp * 0.17087277 )))))))))

  if ( x.lt.0.0 ) tmp_val = 2.0 - tmp_val

  my_erf = 1.0 - tmp_val

  end function my_erf


   subroutine brook_corey_swc_smp(soil,ssnow)
      type(soil_parameter_type), intent(inout) :: soil
      type(soil_snow_type),      intent(inout) :: ssnow

      real(r_2), dimension(mp,ms+1) :: wb_temp
      real(r_2), dimension(mp) :: s_mid
      integer :: i,k

      if (gw_params%ssgw_ice_switch) then
         wb_temp(:,1:ms) = ssnow%wbliq(:,:)
         wb_temp(:,ms+1) = ssnow%GWwb(:)
         do i=1,mp
          if (ssnow%wb(i,ms) .gt. 1.0e-8) then
            wb_temp(i,ms+1) = ssnow%GWwb(i) * ssnow%wbliq(i,ms)/ssnow%wb(i,ms)
          end if
         end do
      else
         wb_temp(:,1:ms) = ssnow%wb(:,:)
         wb_temp(:,ms+1) = ssnow%GWwb(:)
      end if
     !potential from soil water rention function
     !defined as layer average
     do k=1,ms
        do i=1,mp
           s_mid(i) = (wb_temp(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
               (soil%ssat_vec(i,k)-soil%watr(i,k))

           s_mid(i) = min(max(s_mid(i),0.01_r_2),1._r_2)

           ssnow%smp(i,k) = -soil%sucs_vec(i,k)*s_mid(i)**(-soil%bch_vec(i,k))

           ssnow%smp(i,k) = max(min(ssnow%smp(i,k),-soil%sucs_vec(i,k)),sucmin)

           ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/&
                             (max(s_mid(i)*(wb_temp(i,k)-soil%watr(i,k)),0.01_r_2)) !MMY
                     ! MMY BUG in (max(s_mid(i)*(soil%ssat_vec(i,k)-soil%watr(i,k)),0.01_r_2))
        end do
     end do
     ! _____________________ MMY BUG it forgot to calculate the aquifer __________________
          !Aquifer potential
          do i=1,mp
             s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
             s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)

             ssnow%GWsmp(i)    = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))
             ssnow%GWsmp(i)    = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)
             ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/&
                                  (max(s_mid(i)*(wb_temp(i,ms+1)-soil%GWwatr(i)),0.01_r_2)) !MMY
          end do
     ! ____________________________________________________________________________________
  end subroutine brook_corey_swc_smp

   subroutine hutson_cass_swc_smp(soil,ssnow)
      type(soil_parameter_type), intent(inout) :: soil
      type(soil_snow_type),      intent(inout) :: ssnow

      real(r_2), dimension(mp,ms+1) :: wb_temp
      real(r_2), dimension(mp) :: s_mid
      integer :: i,k

   do k=1,ms
      do i=1,mp
         soil%wbc_vec(i,k) = 2.0*soil%bch_vec(i,k)*soil%ssat_vec(i,k)/&
                                 (1.0+2.0*soil%bch_vec(i,k))
         soil%smpc_vec(i,k) = -soil%sucs_vec(i,k) * &
                               (2.0*soil%bch_vec(i,k)/&
                               ((1.0+2.0*soil%bch_vec(i,k))))**(-soil%bch_vec(i,k))
      end do
   end do

   do i=1,mp
      soil%wbc_GW(i) = 2.0*soil%GWbch_vec(i)*soil%GWssat_vec(i)/&
                                 (1.0+2.0*soil%GWbch_vec(i))
      soil%smpc_GW(i) = -soil%GWsucs_vec(i) * &
                               (2.0*soil%GWbch_vec(i)/&
                               ((1.0+2.0*soil%GWbch_vec(i))))**(-soil%GWbch_vec(i))
   end do
    !potential from soil water rention function
    !defined as layer average
    do k=1,ms
       do i=1,mp
          s_mid(i) = (wb_temp(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
              (soil%ssat_vec(i,k)-soil%watr(i,k))

          s_mid(i) = min(max(s_mid(i),0.01_r_2),0.9999999_r_2)

          if (wb_temp(i,k) .lt. soil%wbc_vec(i,k)) then
             ssnow%smp(i,k) = -soil%sucs_vec(i,k)*s_mid(i)**(-soil%bch_vec(i,k))

             ssnow%smp(i,k) = max(min(ssnow%smp(i,k),-soil%sucs_vec(i,k)),sucmin)

             ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/&
                               (max(s_mid(i)*(wb_temp(i,k)-soil%watr(i,k)),0.001_r_2)) ! MMY
                               ! MMY BUG in (max(s_mid(i)*(soil%ssat_vec(i,k)-soil%watr(i,k)),0.001_r_2))
          else
            !!! MMY BUG: in Hutson & Cass 1987 using 'wb/ssat' as independent factor
            !!! MMY      but here Mark uses effective saturation (wb-watr)/(ssat-watr)
            !!! MMY      in wet soil and wb/ssat in wet soil. This leads to a iscontinuity
            !!! MMY      in the water retention curve (smp-wb relationship), and causes
            !!! MMY      abs(smp) in wet soil is larger than in dry soil. I suggest to
            !!! MMY      use wb/ssat in HC_SWC
             ssnow%smp(i,k) = -soil%sucs_vec(i,k)* sqrt(1._r_2 - wb_temp(i,k)/soil%ssat_vec(i,k))/&
                               (sqrt(1._r_2-soil%wbc_vec(i,k)/soil%ssat_vec(i,k))*&
                              ( (soil%wbc_vec(i,k)/soil%ssat_vec(i,k))**soil%bch_vec(i,k) ))

             ssnow%dsmpdw(i,k) = -ssnow%smp(i,k)/(soil%ssat_vec(i,k)*&
                                (1._r_2 - wb_temp(i,k)/soil%ssat_vec(i,k)) )
          end if
       end do
    end do

    !Aquifer potential
    do i=1,mp
       s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.01_r_2),1._r_2)
       if (wb_temp(i,ms+1) .lt. soil%wbc_GW(i)) then
          ssnow%GWsmp(i) = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))

          ssnow%GWsmp(i) = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)

          ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/&
                    (max(s_mid(i)*(soil%GWssat_vec(i)-soil%GWwatr(i)),0.001_r_2))
       else
          ssnow%GWsmp(i) = -soil%GWsucs_vec(i)* sqrt(1._r_2 - wb_temp(i,ms+1)/soil%GWssat_vec(i))/&
                            (sqrt(1._r_2-soil%wbc_GW(i)/soil%GWssat_vec(i))*&
                           ( (soil%wbc_GW(i)/soil%GWssat_vec(i))**soil%GWbch_vec(i) ))

          ssnow%GWsmp(i) = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)

          ssnow%GWdsmpdw(i) = -ssnow%GWsmp(i)/(soil%GWssat_vec(i)*&
                             (1._r_2 - wb_temp(i,ms+1)/soil%GWssat_vec(i)) )
       end if
    end do
   end subroutine hutson_cass_swc_smp

   subroutine brook_corey_hysteresis_swc_smp(soil,ssnow)
   type(soil_parameter_type), intent(inout) :: soil
   type(soil_snow_type),      intent(inout) :: ssnow

   real(r_2), dimension(mp,ms) :: delta_wbliq
   real(r_2), dimension(mp) :: s_mid
   integer :: i,k,klev
   real(r_2), dimension(mp,ms+1) :: wb_temp
   real(r_2) :: tmp_smp

      if (gw_params%ssgw_ice_switch) then
         wb_temp(:,1:ms) = ssnow%wbliq(:,:)
         wb_temp(:,ms+1) = ssnow%GWwb(:) * ssnow%wbliq(:,ms)/soil%ssat_vec(:,ms)
      else
         wb_temp(:,1:ms) = ssnow%wb(:,:)
         wb_temp(:,ms+1) = ssnow%GWwb(:)
      end if

     !ensure sucs_hys is set
     ssnow%sucs_hys = ssnow%hys_fac * soil%sucs_vec
!    do k=1,ms
!       do i=1,mp
!          !on dry or wet?
!          IF (ssnow%hys_fac(i,k) .gt. 0.75) then
!             !on drying
!             ssnow%watr_hys(i,k) = soil%watr(i,k)
!             ssnow%hys_fac(i,k)      = 1.0  !drying
!             ssnow%sucs_hys(i,k) = ssnow%hys_fac(i,k)*soil%sucs_vec(i,k)
!             if (ssnow%smp_hys(i,k) .lt. -ssnow%sucs_hys(i,k)) then
!                tmp_smp = (min(0.99,-ssnow%smp_hys(i,k)/ssnow%sucs_hys(i,k)))**(-soil%bch_vec(i,k))
!                ssnow%ssat_hys(i,k) = (ssnow%wb_hys(i,k)  - ssnow%watr_hys(i,k)*(1.0-tmp_smp))/&
!                                       tmp_smp
!             else
!                ssnow%ssat_hys(i,k) = ssnow%wb_hys(i,k)
!             end if
!             ssnow%ssat_hys(i,k) = max(0.001+ssnow%watr_hys(i,k),ssnow%ssat_hys(i,k))
!
!          ELSE  !wetting
!             ssnow%ssat_hys(i,k) = gw_params%ssat_wet_factor*soil%ssat_vec(i,k)
!             ssnow%hys_fac(i,k)        = 0.5  !wetting
!             ssnow%sucs_hys(i,k) = ssnow%hys_fac(i,k)*soil%sucs_vec(i,k)
!
!             if (ssnow%smp_hys(i,k) .lt. -ssnow%sucs_hys(i,k) ) then
!                tmp_smp = (min(0.99,-ssnow%smp_hys(i,k)/ssnow%sucs_hys(i,k)))**(-soil%bch_vec(i,k))
!                ssnow%watr_hys(i,k) = (ssnow%wb_hys(i,k) - ssnow%ssat_hys(i,k) *tmp_smp)/&
!                                      (1._r_2 - tmp_smp)
!             else
!                ssnow%watr_hys(i,k) = soil%watr(i,k)
!             end if
!
!             ssnow%watr_hys(i,k) = max(0._r_2,min(ssnow%ssat_hys(i,k)-0.001, ssnow%watr_hys(i,k) ) )
!          END IF
!       END DO
!    END DO

    do k=1,ms
       do i=1,mp
          s_mid(i) = (wb_temp(i,k)        - ssnow%watr_hys(i,k)) / &
                     (ssnow%ssat_hys(i,k) - ssnow%watr_hys(i,k))

          s_mid(i) = min(max(s_mid(i),0.01_r_2),1._r_2)
          ! MMY ??? why 0.01_r_2 here but 0.001_r_2 in aquifer
          ssnow%smp(i,k) = -ssnow%sucs_hys(i,k)*s_mid(i)**(-soil%bch_vec(i,k))

          ssnow%smp(i,k) = max(min(ssnow%smp(i,k),-ssnow%sucs_hys(i,k)),sucmin)

          ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/s_mid(i)
       end do
    end do

    !Aquifer potential
    do i=1,mp
       s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)
       ! MMY ???
       ssnow%GWsmp(i)    = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))
       ssnow%GWsmp(i)    = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)
       ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/&
                        (max(s_mid(i)*(wb_temp(i,ms+1)-soil%GWwatr(i)),0.01_r_2)) !MMY
    end do

  end subroutine brook_corey_hysteresis_swc_smp

subroutine swc_hyst_direction(soil,ssnow,veg)
   type(soil_parameter_type), intent(inout) :: soil
   type(soil_snow_type),      intent(inout) :: ssnow
   TYPE(veg_parameter_type) , INTENT(INOUT)    :: veg  ! veg parameters

   real(r_2), dimension(mp,ms) :: delta_wbliq,psi_tmp
   integer :: i,k,klev
   real(r_2) :: tmp_smp

   delta_wbliq = ssnow%wbliq - ssnow%wbliq_old
   !soil hydraulic state/props  so smp out matches the wb out
   !call swc_smp_dsmpdw(soil,ssnow)
   if (gw_params%BC_hysteresis) then
      call brook_corey_hysteresis_swc_smp(soil,ssnow)
   elseif (gw_params%HC_SWC) then
      call hutson_cass_swc_smp(soil,ssnow)
   else
      call brook_corey_swc_smp(soil,ssnow)
   end if

     !switch drying/wetting curve
    do k=1,ms
       do i=1,mp
          if (delta_wbliq(i,k) .gt. 0.0 .and. &
              ssnow%hys_fac(i,k) .gt. 0.75.and.&
              ssnow%wb(i,k).lt.gw_params%ssat_wet_factor * soil%ssat_vec(i,k)) then  !avoid testing .eq. 1.0
           !layer was drying now it is wetting!
              ssnow%smp_hys(i,k) = ssnow%smp(i,k)
              ssnow%wb_hys(i,k)  = ssnow%wb(i,k)
              ssnow%ssat_hys(i,k)= gw_params%ssat_wet_factor * soil%ssat_vec(i,k)
              ssnow%hys_fac(i,k) = 0.5
              ssnow%sucs_hys(i,k) = ssnow%hys_fac(i,k)*soil%sucs_vec(i,k)
          elseif (delta_wbliq(i,k) .le. 0.0 .and. &
                  ssnow%hys_fac(i,k) .lt. 0.75 .and.&
                  ssnow%wb(i,k) .gt. soil%watr(i,k)) then
             !swtiched wetting to drying
             ssnow%smp_hys(i,k) = ssnow%smp(i,k)
             ssnow%wb_hys(i,k)  = ssnow%wb(i,k)
             ssnow%hys_fac(i,k) = 1.0
             ssnow%sucs_hys(i,k) = ssnow%hys_fac(i,k)*soil%sucs_vec(i,k)
          end if

          if (delta_wbliq(i,k) .gt. 0.0) then

             ssnow%ssat_hys(i,k)= gw_params%ssat_wet_factor * soil%ssat_vec(i,k)
             if (ssnow%smp_hys(i,k) .lt. -ssnow%sucs_hys(i,k) ) then

                tmp_smp = (min(0.99,-ssnow%smp_hys(i,k)/ssnow%sucs_hys(i,k)))**(-soil%bch_vec(i,k))
                ssnow%watr_hys(i,k) = (ssnow%wb_hys(i,k) - soil%ssat_vec(i,k) * tmp_smp)/&
                                     (1.0 - tmp_smp)
             else
                ssnow%watr_hys(i,k) = ssnow%wb_hys(i,k)
             end if

             ssnow%watr_hys(i,k) = max(0._r_2,min(ssnow%ssat_hys(i,k)-0.001, ssnow%watr_hys(i,k) ) )

          else !drying ! MMY add a space
             ssnow%watr_hys(i,k)= soil%watr(i,k)  !this is how we tell if drying
                                               !watr_hys .eq. watr
             if (ssnow%smp_hys(i,k) .lt. -ssnow%sucs_hys(i,k)) then
                tmp_smp = (min(0.99,-ssnow%smp_hys(i,k)/ssnow%sucs_hys(i,k)))**(-soil%bch_vec(i,k))
                ssnow%ssat_hys(i,k) = (ssnow%wb_hys(i,k)  - soil%watr(i,k)*(1.0-tmp_smp))/&
                                       tmp_smp
             else
                ssnow%ssat_hys(i,k) = ssnow%wb_hys(i,k)
             end if

             ssnow%ssat_hys(i,k) = max(ssnow%watr_hys(i,k)+0.001,ssnow%ssat_hys(i,k))
          end if

       end do
    end do

  if (cable_user%gw_model .and. gw_params%bc_hysteresis) then
      do klev=1,ms
         do i=1,mp
            if (soil%isoilm(i) .ne. 9 .and. veg%iveg(i) .le. 16) then

               psi_tmp(i,klev) = abs(psi_c(veg%iveg(i)))

               soil%swilt_vec(i,klev) = (ssnow%ssat_hys(i,klev)-ssnow%watr_hys(i,klev)) * &
                                        (psi_tmp(i,klev)/soil%sucs_vec(i,klev))&
                                         **(-1.0/soil%bch_vec(i,klev))+&
                                        ssnow%watr_hys(i,klev)
               soil%sfc_vec(i,klev) = (gw_params%sfc_vec_hk/soil%hyds_vec(i,klev))&
                                       **(1.0/(2.0*soil%bch_vec(i,klev)+3.0)) *&
                                       (ssnow%ssat_hys(i,klev)-ssnow%watr_hys(i,klev)) + ssnow%watr_hys(i,klev)
            end if
        end do
     end do

   end if

    end subroutine swc_hyst_direction

!`!taken from
!`!/g/data1/w35/mrd561/CABLE/CMIP6-GM2_dev_testing_tiles/core/biogeophys
!`! to test
!`SUBROUTINE GWsoilfreeze(dels, soil, ssnow)
!`  !NOTE: this is only included because gw_model uses parameters XXX_vec
!`  !these are r_2.  this breaks bitwise compatibility with trunk
!`  !if acceptable this routine does the same thing but with r_2 soil params
!`  ! if max_ice_frac always set to frozen_limit and tgg_tmp is always C%TFRZ
!`
!`   REAL, INTENT(IN)                    :: dels ! integration time step (s)
!`   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
!`   TYPE(soil_parameter_type), INTENT(INOUT)    :: soil
!`   REAL     , DIMENSION(mp,ms) :: tgg_old,tgg_new,tgg_tmp !tgg_old is previous point of when crosses freezing
!`   REAL(r_2), DIMENSION(mp,ms)           :: sicefreeze
!`   REAL(r_2), DIMENSION(mp,ms)           :: sicemelt
!`   REAL(r_2), DIMENSION(mp,ms)        :: wbice_delta,avail_por,delta_ice_vol
!`   REAL(r_2), DIMENSION(mp)           :: ice_mass,liq_mass,tot_mass
!`   INTEGER :: i,j,k
!`   REAL(r_2) :: func,funcderv,Aconst,Dconst,t_zero,t_one,dtmp
!`   REAL, DIMENSION(mp,ms) :: gammzz_snow
!`   REAL(r_2),DIMENSION(mp,ms) :: xx,max_ice_frac,iceF,den_css  !Decker and Zeng 20.9
!`   REAL(r_2) :: delta_wbliq,tmp_var
!`
!`   max_ice_frac(:,:) = zero
!`   delta_ice_vol(:,:) = zero
!`   tgg_old(:,:) = ssnow%otgg(:,:)
!`   tgg_new(:,:) = ssnow%tgg(:,:)
!`   tgg_tmp(:,:) = tgg_old(:,:)
!`
!`   gammzz_snow(:,:) = zero
!`   k=1
!`   do i=1,mp
!`      if (ssnow%isflag(i) .eq. zero.and. soil%isoilm(i) .ne. 9) then
!`           gammzz_snow(i,k) = real(C%cgsnow,r_2) * real(ssnow%snowd(i),r_2)
!`      end if
!`   end do
!`
!`   sicefreeze(:,:) = -1._r_2
!`   sicemelt(:,:) = -1._r_2
!`
!`   do k=1,ms
!`   do i=1,mp
!`
!`      ssnow%wmice(i,k)  = ssnow%wbice(i,k)*soil%zse_vec(i,k)*real(C%density_ice,r_2)
!`      ssnow%wmliq(i,k)  = ssnow%wbliq(i,k)*soil%zse_vec(i,k)*real(C%density_liq,r_2)
!`      ssnow%wmtot(i,k)  = ssnow%wmice(i,k) + ssnow%wmliq(i,k)
!`
!`      if  ((ssnow%tgg(i,k) .lt. C%TFRZ)  .and. &
!`           (ssnow%tgg(i,k) .lt. ssnow%otgg(i,k))) then
!`
!`            ssnow%otgg(i,k) = min(ssnow%otgg(i,k),C%TFRZ)
!`
!`            if (ssnow%wb(i,k) .gt. (soil%watr(i,k)+1.0e-6) ) then
!`               iceF(i,k) = max(zero,min(0.95_r_2,max(zero,ssnow%wbliq(i,k)-soil%watr(i,k))/&
!`                                 (ssnow%wb(i,k)-soil%watr(i,k))  )  )
!`            else
!`              iceF(i,k) = 0._r_2
!`            end if
!`            tgg_tmp(i,k) = iceF(i,k)*ssnow%otgg(i,k) + &
!`                           (1._r_2 - iceF(i,k))*ssnow%tgg(i,k)
!`
!`            Aconst =max(0.1,min(0.99_r_2,(ssnow%wb(i,k)-soil%watr(i,k))/(soil%ssat_vec(i,k)-soil%watr(i,k)) ) )
!`            Dconst = exp(1. - Aconst  )
!`
!`            if (tgg_tmp(i,k) .lt. C%TFRZ) then
!`            max_ice_frac(i,k) = max(zero,ssnow%wb(i,k)*&
!`                                 (1._r_2 - exp(2._r_2*(tgg_tmp(i,k)-C%TFRZ)*Aconst*Aconst)) / Dconst)
!`
!`            else
!`                max_ice_frac(i,k) = zero
!`            end if
!`
!`            max_ice_frac(i,k) = min(ssnow%wb(i,k)*0.9, max_ice_frac(i,k) )
!`
!`            delta_ice_vol(i,k) = max(zero, max_ice_frac(i,k) - ssnow%wbice(i,k) )
!`
!`            !check amount of water we have
!`            delta_ice_vol(i,k) = min(ssnow%wbliq(i,k)*real(C%density_liq/C%density_ice,r_2), delta_ice_vol(i,k) )
!`
!`            sicefreeze(i,k) = min(delta_ice_vol(i,k)*soil%zse_vec(i,k)*c%density_ice, &
!`                                     max(zero,(ssnow%otgg(i,k)-ssnow%tgg(i,k))*ssnow%gammzz(i,k)/C%HLF) )
!`
!`      elseif ((ssnow%tgg(i,k) .gt. C%TFRZ) .and. &
!`              (ssnow%tgg(i,k) .gt. ssnow%otgg(i,k)) .and. &
!`              ssnow%wbice(i,k) .gt. zero ) then
!`
!`            ssnow%otgg(i,k) = min(ssnow%otgg(i,k),C%TFRZ)
!`
!`            tgg_tmp(i,k) = C%TFRZ
!`
!`            delta_ice_vol(i,k) = ssnow%wbice(i,k)
!`
!`            sicemelt(i,k) = min(delta_ice_vol(i,k)*soil%zse_vec(i,k)*c%density_ice, &
!`                                     max(zero,(ssnow%tgg(i,k)-ssnow%otgg(i,k))*ssnow%gammzz(i,k)/C%HLF))
!`
!`      endif
!`   end do
!`   end do
!`
!`   DO k = 1, ms
!`   DO i=1,mp
!`
!`
!`      if (sicefreeze(i,k) .gt. zero .and. ssnow%tgg(i,k).lt.c%tfrz .and. &
!`               ssnow%tgg(i,k).lt.ssnow%otgg(i,k) ) then
!`
!`
!`         ssnow%wbice(i,k) = ssnow%wbice(i,k) +&
!`                                  sicefreeze(i,k)/soil%zse_vec(i,k)/real(C%density_ice,r_2)
!`
!`         delta_wbliq = max(zero,min(ssnow%wbliq(i,k),delta_ice_vol(i,k)*real(C%density_ice/C%density_liq,r_2)))
!`
!`         ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),  &
!`                              (1.0- soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k) &
!`                            + (ssnow%wbliq(i,k) - delta_wbliq) * REAL(C%cswat*C%density_liq,r_2)   &
!`                            + ssnow%wbice(i,k) * REAL(C%csice*C%density_ice,r_2)&
!`                            )*soil%zse_vec(i,k) + gammzz_snow(i,k)
!`
!`         ssnow%tgg(i,k) = ssnow%tgg(i,k) + real(sicefreeze(i,k) )&
!`                             * C%hlf / real(ssnow%gammzz(i,k) )
!`
!`         ssnow%wmice(i,k) = ssnow%wbice(i,k)*soil%zse_vec(i,k)*real(C%density_ice,r_2)
!`         ssnow%wmliq(i,k) = ssnow%wmtot(i,k) - ssnow%wmice(i,k)
!`
!`         ssnow%wbliq(i,k) = ssnow%wmliq(i,k) / (soil%zse_vec(i,k)*m2mm)
!`         ssnow%wb(i,k)    = ssnow%wbliq(i,k) + den_rat * ssnow%wbice(i,k)
!`
!`         elseif (sicemelt(i,k) .gt. zero .and. ssnow%tgg(i,k).gt.c%tfrz &
!`                .and. ssnow%wbice(i,k) .gt. 0._r_2) then
!`
!`         ssnow%wbice(i,k) = ssnow%wbice(i,k)  -  delta_ice_vol(i,k)
!`
!`         delta_wbliq = delta_ice_vol(i,k)*real(C%density_ice/C%density_liq,r_2)
!`
!`         ssnow%gammzz(i,k) =max(soil%heat_cap_lower_limit(i,k), &
!`                           (1.0 - soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k)&
!`              + (ssnow%wbliq(i,k)+delta_wbliq) * REAL(C%cswat*C%density_liq,r_2)   &
!`              + ssnow%wbice(i,k) * REAL(C%csice*C%density_ice,r_2)&
!`                             )*soil%zse_vec(i,k)  + gammzz_snow(i,k)
!`
!`         ssnow%tgg(i,k) = ssnow%tgg(i,k) - real(sicemelt(i,k) )&
!`                             * C%hlf / real(ssnow%gammzz(i,k) )
!`
!`         ssnow%wmice(i,k) = ssnow%wbice(i,k)*soil%zse_vec(i,k)*real(C%density_ice,r_2)
!`         ssnow%wmliq(i,k) = ssnow%wmtot(i,k) - ssnow%wmice(i,k)
!`
!`         ssnow%wbliq(i,k) = ssnow%wmliq(i,k) / (soil%zse_vec(i,k)*m2mm)
!`         ssnow%wb(i,k)    = ssnow%wbliq(i,k) + den_rat * ssnow%wbice(i,k)
!`
!`      END IF
!`
!`
!`   END DO
!`   END DO
!`
!`END SUBROUTINE GWsoilfreeze

  subroutine set_den_rat()
      den_rat = real(C%density_ice/C%density_liq,r_2)
  end subroutine set_den_rat

END MODULE cable_gw_hydro_module
