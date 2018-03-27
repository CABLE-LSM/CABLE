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
                                   max_glacier_snowd,ktau_gl

   USE cable_soil_snow_module, ONLY : trimb, snow_processes_soil_thermal

   USE cable_data_module, only: C=>PHYS!issnow_type,point2constants


   IMPLICIT NONE

   PRIVATE

   !mrd561 GW params
   REAL(r_2), PARAMETER :: sucmin       = -1.0e12      ! minimum soil pressure head [mm]
   REAL(r_2), PARAMETER :: volwatmin    = 1e-4        !min soil water [mm]      
   REAL(r_2), PARAMETER :: wtd_uncert   = 0.1         ! uncertaintiy in wtd calcultations [mm]
   REAL(r_2), PARAMETER :: wtd_max      = 1000000.0   ! maximum wtd [mm]
   REAL(r_2), PARAMETER :: wtd_min      = 10.0        ! minimum wtd [mm]
   REAL(r_2), PARAMETER :: close_to_one = 0.9999
   REAL(r_2), PARAMETER :: m_to_mm = 1000.0
   REAL(r_2), PARAMETER :: mm_to_m = 0.001
   REAL(r_2), SAVE :: r2pi=3.14159,&
                      r2_density_ice=921.0,&
                      r2_density_liq=1000.0,&
                      den_rat=0.921,&
                      r2_cgsnow = 2090.0,&     
                      r2_cs_rho_ice = 1.9341e6,&  
                      r2_cs_rho_wat = 4.218e6,&   
                      r2_csice = 2.100e3,&     
                      r2_cswat = 4.218e3  

  INTEGER, PARAMETER :: wtd_iter_max = 20 ! maximum number of iterations to find the water table depth                    

  ! ! This module contains the following subroutines that
  !are called from other modules
   PUBLIC :: soil_snow_gw,calc_srf_wet_fraction,sli_hydrology,& 
             pore_space_relative_humidity,set_unsed_gw_vars,&
             set_r2_gw_params,set_wbliq

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
   REAL,      DIMENSION(mp,ms)        :: tgg_old
   REAL(r_2), DIMENSION(mp)           :: sicefreeze
   REAL(r_2), DIMENSION(mp)           :: sicemelt
   REAL(r_2), DIMENSION(mp,ms)        :: wbice_delta,avail_por
   REAL(r_2), DIMENSION(mp)           :: ice_mass,liq_mass,tot_mass
   INTEGER :: i,j,k
   INTEGER, DIMENSION(mp,ms) :: phase_change_needed
   REAL, DIMENSION(mp,ms) :: tgg_tmp,t_guess,tgg_new
   REAL ,DIMENSION(mp,ms) :: avail_energy
   REAL(r_2),DIMENSION(mp,ms) :: satF,xx,max_ice_frac,iceF,den_css  !Decker and Zeng 2009
   real(r_2) :: wgt

   phase_change_needed(:,:) = 0
   tgg_tmp(:,:) = ssnow%otgg(:,:)


   do k=1,ms
   do i=1,mp
      if (ssnow%tgg(i,k) .lt. ssnow%otgg(i,k) .and. ssnow%tgg(i,k) .le. C%TFRZ) then
         phase_change_needed(i,k) = 2
         if (ssnow%otgg(i,k) .gt. C%TFRZ) then
            tgg_tmp(i,k) = C%TFRZ
         else
            tgg_tmp(i,k) =  ssnow%otgg(i,k)
         end if
         !negative
      elseif ((ssnow%tgg(i,k) .gt. ssnow%otgg(i,k)) .and. &
              (ssnow%tgg(i,k) .gt. gw_params%Tmlt)) then
         if (ssnow%otgg(i,k) .lt. gw_params%Tmlt) then
            tgg_tmp(i,k) = gw_params%Tmlt
         else 
            tgg_tmp(i,k) = min(ssnow%otgg(i,k), C%TFRZ)
         end if
         phase_change_needed(i,k) = -2     
      end if

   end do
   end do
 
   max_ice_frac(:,:) = 0.0

   do k=1,ms
   do i=1,mp
   
      if  (phase_change_needed(i,k) .gt. 1) then
         !when freezing weight toward needing to freeze ice when little ice
         if (ssnow%wb(i,k) .gt. 1.0e-4) then
            wgt = 1.0 - max(0.05,min(0.95,ssnow%wbliq(i,k)/ssnow%wb(i,k)))
         else
            wgt = 1.0
         end if
         t_guess(i,k) = wgt*ssnow%tgg(i,k)+ (1.0-wgt)*tgg_tmp(i,k)
      elseif (phase_change_needed(i,k) .lt. -1) then
         if (ssnow%wb(i,k) .gt. 1.0e-4) then
            wgt = max(0.05,min(0.95,ssnow%wbliq(i,k)/ssnow%wb(i,k)))
         else
            wgt = 0.0
         end if
         t_guess(i,k) = wgt*ssnow%tgg(i,k)+ (1.0-wgt)*tgg_tmp(i,k)
      else
         t_guess(i,k) =ssnow%tgg(i,k)
      endif
   end do
  end do


   max_ice_frac(:,:) = 0._r_2
   do k=1,ms
   do i=1,mp
      if (t_guess(i,k) .lt. C%TFRZ)  then
          wgt = min(0.9,max(0.01,ssnow%wbice(i,k)/max(1.0e-4,ssnow%wmtot(i,k)) ))

          max_ice_frac(i,k) = ssnow%wb(i,k) * &
                              min(0.9,max(0.0, (1.0 - exp(-2.0*(wgt**4.0)*(C%TFRZ-t_guess(i,k))))/&
                              exp(1.0 - wgt) ) )

!          max_ice_frac(i,k) = min((ssnow%wb(i,k)-soil%watr(i,k)), &
!                              (ssnow%wb(i,k) - soil%watr(i,k)-&
!                                (soil%ssat_vec(i,k)-soil%watr(i,k)) *&
!                              ( abs((1000.0*C%hlf/C%grav/soil%sucs_vec(i,k)*&
!                              (t_guess(i,k) - C%TFRZ/C%TFRZ)))**(-1.0/soil%bch_vec(i,k)))))
!
      end if
   end do
   end do 


   !allow more freezing for permenant glacier ice regions
   do i=1,mp
      if (soil%isoilm(i) .eq. 9) max_ice_frac(i,:) = 0.85_r_2*ssnow%wb(i,:)*den_rat
   end do

   ssnow%wmliq  = ssnow%wbliq * soil%zse_vec * r2_density_liq
   ssnow%wmice  = ssnow%wbice * soil%zse_vec * r2_density_ice
   ssnow%wmtot  = ssnow%wmice + ssnow%wmliq

   DO k = 1, ms
   DO i=1,mp
      if (phase_change_needed(i,k) .gt. 1 .and. &
           max_ice_frac(i,k) - ssnow%wbice(i,k)  .gt. 1.0e-6) then
         
         sicefreeze(i) = MIN( MAX( 0.0_r_2, (max_ice_frac(i,k)-ssnow%wbice(i,k)))*&  
                         soil%zse_vec(i,k) * r2_density_ice,             &
                      real( tgg_tmp(i,k) - ssnow%tgg(i,k),r_2 ) &
                      * ssnow%gammzz(i,k) / real(C%hlf,r_2) )

         sicefreeze(i) = min(sicefreeze(i),ssnow%wbliq(i,k)*soil%zse_vec(i,k) * r2_density_liq )

         ssnow%wbice(i,k) = ssnow%wbice(i,k) +&
                                  sicefreeze(i)/(soil%zse_vec(i,k)*r2_density_ice)

         ssnow%wbliq(i,k) = ssnow%wbliq(i,k) -  sicefreeze(i)/&
                                    (soil%zse_vec(i,k)*r2_density_liq)

         ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),  &
                              (1.0-soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k) &
                            + ssnow%wbliq(i,k) * r2_cs_rho_wat   &
                            + ssnow%wbice(i,k) * r2_cs_rho_ice   &
                            )*soil%zse_vec(i,k) 

        if (k .eq. 1 .and. ssnow%isflag(i) .eq. 0) then
           ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + r2_cgsnow * real(ssnow%snowd(i),r_2)
        end if

         ssnow%tgg(i,k) = ssnow%tgg(i,k) + REAL(sicefreeze(i))                    &
                             * C%hlf / REAL(ssnow%gammzz(i,k) )

         ssnow%wb(i,k) = ssnow%wbliq(i,k) + den_rat*ssnow%wbice(i,k)

         ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*soil%zse_vec(i,k)*r2_density_liq
         ssnow%wmice(i,k) = ssnow%wbice(i,k)*soil%zse_vec(i,k)*r2_density_ice
         ssnow%wmtot(i,k) = ssnow%wmliq(i,k) + ssnow%wmice(i,k)

     elseif (phase_change_needed(i,k) .lt. -1 .and. &
              ssnow%wbice(i,k) .gt.  max_ice_frac(i,k) ) then

         sicemelt(i) = MAX(0._r_2, MIN((ssnow%wbice(i,k)-max_ice_frac(i,k))*&
                       soil%zse_vec(i,k) * r2_density_ice, &
                       real( ssnow%tgg(i,k) - tgg_tmp(i,k),r_2 ) &
                        * ssnow%gammzz(i,k) / real(C%hlf,r_2) ))
         
         ssnow%wbice(i,k) = MAX( 0.0_r_2, ssnow%wbice(i,k) - sicemelt(i)          &
                            / (soil%zse_vec(i,k) * r2_density_ice) )

         ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + sicemelt(i)/&
                                    (soil%zse_vec(i,k)*r2_density_liq)

         ssnow%gammzz(i,k) =max(soil%heat_cap_lower_limit(i,k), &
                           (1.0-soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k)&
              + ssnow%wbliq(i,k) * r2_cs_rho_wat     &
              + ssnow%wbice(i,k) * r2_cs_rho_ice     &
                             )*soil%zse_vec(i,k)            

         if (k .eq. 1 .and. ssnow%isflag(i) .eq. 0) then
           ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + r2_cgsnow * real(ssnow%snowd(i),r_2)
         end if
         ssnow%tgg(i,k) = ssnow%tgg(i,k) - REAL(sicemelt(i))                     &
                          * C%hlf / REAL(ssnow%gammzz(i,k))

         ssnow%wb(i,k) = ssnow%wbliq(i,k) + den_rat*ssnow%wbice(i,k)

         ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*soil%zse_vec(i,k)*r2_density_liq
         ssnow%wmice(i,k) = ssnow%wbice(i,k)*soil%zse_vec(i,k)*r2_density_ice
         ssnow%wmtot(i,k) = ssnow%wmliq(i,k) + ssnow%wmice(i,k)

      !if for some reason end up here     
!      ELSEIF( tgg_tmp(i,k) .ge. C%TFRZ .and. &
!              ssnow%tgg(i,k) > tgg_tmp(i,k) .AND. ssnow%wbice(i,k) > 0. ) then
!         
!         sicemelt(i) = MIN( ssnow%wbice(i,k) * soil%zse_vec(i,k) * r2_density_ice,              &
!                    ( ssnow%tgg(i,k) - tgg_tmp(i,k) ) * ssnow%gammzz(i,k) / C%hlf )
!         
!         ssnow%wbice(i,k) = MAX( 0.0_r_2, ssnow%wbice(i,k) - sicemelt(i)          &
!                            / (soil%zse_vec(i,k) * r2_density_ice) )
!
!         ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + sicemelt(i)/&
!                                    (soil%zse_vec(i,k)*r2_density_liq)
!
!         ssnow%gammzz(i,k) =max(soil%heat_cap_lower_limit(i,k), &
!                           (1.0-soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k)&
!              + ssnow%wbliq(i,k) * r2_cs_rho_wat   &
!              + ssnow%wbice(i,k) * r2_cs_rho_ice&
!                             )*soil%zse_vec(i,k)            
!
!         if (k .eq. 1 .and. ssnow%isflag(i) .eq. 0) then
!           ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + C%cgsnow * ssnow%snowd(i)
!         end if
!         ssnow%tgg(i,k) = ssnow%tgg(i,k) - REAL(sicemelt(i))                     &
!                          * C%hlf / REAL(ssnow%gammzz(i,k))
!
!
!         ssnow%wb(i,k) = ssnow%wbliq(i,k) + den_rat*ssnow%wbice(i,k)
!         ssnow%wmliq(i,k) = ssnow%wb(i,k)*soil%zse_vec(i,k)*r2_density_liq
!         ssnow%wmice(i,k) = ssnow%wb(i,k)*soil%zse_vec(i,k)*r2_density_ice
!         ssnow%wmtot(i,k) = ssnow%wmliq(i,k) + ssnow%wmice(i,k)
!       
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
         zse_mp_mm(i,k)  = real(soil%zse_vec(i,k)*r2_density_liq,r_2)
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
        ssnow%wbliq(:,k) = ssnow%wbliq(:,k) - ssnow%evapfbl(:,k)/zse_mp_mm(:,k)
     ENDDO

  ENDIF

  do k=1,ms
     do i=1,mp
        ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*zse_mp_mm(i,k)!mass
        ssnow%wmtot(i,k) = ssnow%wmliq(i,k) + ssnow%wmice(i,k)  !mass
        ssnow%wb(i,k)    = ssnow%wbliq(i,k) + den_rat*ssnow%wbice(i,k)  !volume
     end do
  end do


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
      efpor(i) = max(0.001_r_2, soil%ssat_vec(i,1) - den_rat*ssnow%wbice(i,1))
      icemass  = ssnow%wbice(i,1) * dzmm *den_rat
      liqmass  = (ssnow%wb(i,1)-den_rat*ssnow%wbice(i,1)) * dzmm
      totmass  = max(liqmass+icemass,real(1e-2,r_2))
      icef(i)     = max(0._r_2,min(1._r_2,icemass / totmass))
   end do

   !sat fraction assuming topo controlled subgrid soil moisture distribution
   !called from cable_canopy for srf wet fraction alrady
   !call saturated_fraction(ssnow,soil,veg)

   !srf frozen fraction.  should be based on topography
   do i = 1,mp
      fice = (exp(-gw_params%IceAlpha*(1._r_2-icef(i)))-exp(-gw_params%IceAlpha))/(1._r_2-exp(-gw_params%IceAlpha))
      fice  = min(max(fice,0._r_2),1._r_2)
      satfrac_liqice(i)   = max(0.,min(0.95,fice + (1._r_2-fice)*ssnow%satfrac(i) ) )
   end do

   do i=1,mp
      tmpa = ssnow%wbliq(i,1) / efpor(i)
      tmpb = max( (tmpa-satfrac_liqice(i))/max(0.01_r_2,(1._r_2-satfrac_liqice(i))), 0._r_2)
      tmpa = -2._r_2*soil%bch_vec(i,1)*soil%sucs_vec(i,1)/dzmm
      qinmax = (1._r_2 + tmpa*(tmpb-1._r_2))*soil%hyds_vec(i,1)*exp(-gw_params%hkrz*(0.5*dzmm/1000.0_r_2-gw_params%zdepth))

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

   if (sli_call .or. cable_runtime%UM) then
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
        wmean(i) = wmean(i) + stot(i,k)*soil%zse(k)*1000._r_2
        ztot(i)  = ztot(i) + soil%zse(k)*1000._r_2
     end do
  end do

  do i=1,mp
     wmean(i) = wmean(i) + ssnow%GWwb(i)/soil%GWssat_vec(i) * soil%GWdz(i)*1000._r_2
     ztot(i)  = ztot(i) + soil%GWdz(i)*1000._r_2

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
        total_depth_column(i) = soil%GWdz(i)*m_to_mm
        def(i) = max(0._r_2,soil%GWssat_vec(i)-ssnow%GWwb(i))*soil%GWdz(i)*r2_density_liq
     end do   

  else
     def(:) = 0._r_2
     total_depth_column(:) = 0._r_2
  end if

  !total depth of soil column
  do k=1,ms
     do i=1,mp
         total_depth_column(i) = total_depth_column(i) + soil%zse_vec(i,k)*m_to_mm
     end do
  end do
  
  !comute the total mass away from full saturation
  do k=1,ms
     do i=1,mp

       def(i) = def(i) +                                                           &
                max(0._r_2,(soil%ssat_vec(i,k)-(ssnow%wb(i,k)))*soil%zse_vec(i,k)*m_to_mm)
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
    if (veg%iveg(i) .eq. 16) then      
       ssnow%wtd(i) = 0.0

    else

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
     ssnow%wtd(i) = min(wtd_max, max(wtd_min,ssnow%wtd(i) ) )
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
    REAL(r_2), DIMENSION(ms)            :: dzmm
    REAL(r_2), DIMENSION(mp,ms)         :: dzmm_mp,dqhdw,macro
    REAL(r_2), DIMENSION(0:ms+1)        :: zimm
    REAL(r_2), DIMENSION(ms)            :: zmm
    REAL(r_2), DIMENSION(mp)            :: GWzimm,xs,zaq,s_mid,GWdzmm
    REAL(r_2), DIMENSION(mp)            :: xs1,GWmsliq!xsi    !mass (mm) of liquid over/under saturation, mass of aquifer water
    REAL(r_2)                           :: xsi
    REAL(r_2), DIMENSION(mp,ms+1)       :: del_wb
    !MD DEBUG VARS
    INTEGER :: imp,ims,k_drain

    dqhdw(:,:) = 0._r_2
    !make code cleaner define these here
    dzmm(:) = 1000.0_r_2 * real(soil%zse(:),r_2)
    do i=1,mp
       dzmm_mp(i,:) = dzmm(:)
    end do

    zimm(0) = 0._r_2
    zimm(1:ms) = zimm(0:(ms-1)) + dzmm(1:ms)
    zmm(1:ms)  = zimm(0:(ms-1)) + 0.5_r_2*dzmm(1:ms)

    do i=1,mp
       GWdzmm(i) = real(soil%GWdz(i),r_2)*1000._r_2
       GWzimm(i) = zimm(ms)+GWdzmm(i)
       zaq(i)    = zimm(ms) + 0.5_r_2*GWdzmm(i)
    end do

    macro(:,:) = real(gw_params%macro_down,r_2)
    
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
    
    !equilibrium water content
    CALL calc_equilibrium_water_content(ssnow,soil)

    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    CALL subsurface_drainage(ssnow,soil,veg,dzmm)

    CALL aquifer_recharge(dels,ssnow,soil,veg,zaq,zmm,dzmm)

    k = 1     !top soil layer
    do i=1,mp
       qin(i)     = ssnow%sinfil(i)
       den(i)     = (zmm(k+1)-zmm(k))
       dne(i)     = macro(i,k)*(ssnow%zq(i,k+1)-ssnow%zq(i,k))
       num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
       dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   +  num(i)*ssnow%dhkdw(i,k))/den(i) 
       dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
       rt(i,k) =  qin(i) - qout(i) - ssnow%qhlev(i,k)*max(0.,ssnow%wbliq(i,k)-soil%watr(i,k))
       at(i,k) =  0._r_2
       bt(i,k) =  dzmm(k)/dels + dqodw1(i) !-ssnow%qhlev(i,k)
       ct(i,k) =  dqodw2(i)   !+ssnow%qhlev(i,k+1)
    end do
    do k = 2, ms - 1     !middle soil layers
       do i=1,mp
          den(i)     = (zmm(k) - zmm(k-1))
          dne(i)     = macro(i,k-1)*(ssnow%zq(i,k)-ssnow%zq(i,k-1))
          num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
          qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
          dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          den(i)     = (zmm(k+1)-zmm(k))
          dne(i)     = macro(i,k)*(ssnow%zq(i,k+1)-ssnow%zq(i,k))
          num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
          qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
          dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i) 
          dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
          rt(i,k) =  qin(i) - qout(i) - ssnow%qhlev(i,k)*max(0.,ssnow%wbliq(i,k)-soil%watr(i,k))
          at(i,k) = -dqidw0(i)
          bt(i,k) =  dzmm(k)/dels - dqidw1(i) + dqodw1(i)! -ssnow%qhlev(i,k)
          ct(i,k) =  dqodw2(i)  !+ssnow%qhlev(i,k+1)
       end do
    end do
       
    k = ms   !Bottom soil layer
    do i=1,mp
       den(i)     = (zmm(k) - zmm(k-1))
       dne(i)     = macro(i,k-1)*(ssnow%zq(i,k)-ssnow%zq(i,k-1))
       num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
       qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
       dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       den(i)     = zaq(i) - zmm(k)
       dne(i)     = macro(i,k)*(ssnow%GWzq(i)-ssnow%zq(i,k))
       num(i)     =  (ssnow%GWsmp(i)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = 0._r_2
       dqodw1(i)  = 0._r_2
       dqodw2(i)  = 0._r_2
       rt(i,k) =  qin(i) - qout(i) - &
                   (ssnow%qhlev(i,k)*max(0.,ssnow%wbliq(i,k)-soil%watr(i,k))+ssnow%Qrecharge(i))
       at(i,k) = -dqidw0(i)
       bt(i,k) =  dzmm(k)/dels - dqidw1(i) + dqodw1(i)! -ssnow%qhlev(i,ms)
       ct(i,k) =  dqodw2(i) !+0.0
    end do
       
    CALL trimb(at,bt,ct,rt,ms)          !use the defulat cable tridiag solution

    !leave here to facilitate including how qhlev changes as function of wbliq
    ssnow%qhz(:) = ssnow%qhlev(:,ms+1)
    do k=1,ms
       do i=1,mp
           ssnow%qhz(i) = ssnow%qhz(i) + ssnow%qhlev(i,k)*(&
                                        max(0.,ssnow%wbliq(i,k)-soil%watr(i,k)))
       end do
    end do

    do k=1,ms
       do i=1,mp
          ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + rt(i,k)
       end do
    end do

    do i=1,mp
       !mreoved because Qrecharge added to solution above
       !snow%wbliq(i,ms) = snow%wbliq(i,ms) -&
       !                    ssnow%Qrecharge(i)*dels/(soil%zse_vec(i,ms)*r2_density_liq)
       ssnow%GWwb(i) = ssnow%GWwb(i)  +  (ssnow%Qrecharge(i)-ssnow%qhlev(i,ms+1))*dels/GWdzmm(i)
    end do

    !determine the available pore space
    !volumetric
    do k=1,ms
       do i=1,mp
          eff_por(i,k)  = max(0._r_2, soil%ssat_vec(i,k) - den_rat*ssnow%wbice(i,k) )
       end do
    end do

    do i=1,mp
       xsi = 0._r_2

       if (ssnow%GWwb(i) .gt. soil%GWssat_vec(i)) then
          xsi = (ssnow%GWwb(i) - soil%GWssat_vec(i))*GWdzmm(i)
          ssnow%GWwb(i) = soil%GWssat_vec(i)
       end if

       do k=1,ms
          if (ssnow%wbliq(i,k) .gt. eff_por(i,k)) then
             xsi = xsi + (ssnow%wbliq(i,k) - eff_por(i,k))*dzmm(k)
             ssnow%wbliq(i,k) = eff_por(i,k)
           end if
       end do

       if (xsi .gt. 0._r_2) then
          if (xsi .lt. (soil%GWssat_vec(i)-ssnow%GWwb(i))*GWdzmm(i)) then
             ssnow%GWwb(i) = ssnow%GWwb(i) + xsi/GWdzmm(i)
             xsi = 0._r_2
          else
             xsi = xsi - max(0._r_2,(soil%GWssat_vec(i)-ssnow%GWwb(i))*GWdzmm(i))
             ssnow%GWwb(i) = soil%GWssat_vec(i)
          end if
       end if
 
       do k = ms,1,-1  !loop from bottom to top adding extra water to each layer
          if (xsi .gt. 0._r_2) then
             if (xsi .lt. (eff_por(i,k)-ssnow%wbliq(i,k))*dzmm(k)) then
                ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + xsi/dzmm(k)
                xsi = 0._r_2
             else
                xsi = xsi - (eff_por(i,k) - ssnow%wbliq(i,k))*dzmm(k)
                ssnow%wbliq(i,k) = eff_por(i,k)
             end if
          end if
       end do  !ms loop

       if (xsi .gt. 0._r_2) then
          ssnow%qhz(i) = ssnow%qhz(i) + xsi/dels
       end if

       xsi = 0._r_2

       do k = 1,ms
          xsi = 0._r_2             !should be a single float (array not needed)
          if (ssnow%wbliq(i,k) .lt. volwatmin) then
             xsi = xsi + (volwatmin - ssnow%wbliq(i,k))*dzmm(k)
             ssnow%wbliq(i,k) = volwatmin
          end if

          if (k .lt. ms) then
             ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) - xsi/dzmm(k+1)
          else
             ssnow%GWwb(i) = ssnow%GWwb(i) - xsi / GWdzmm(i)
          end if

       end do  !ms loop
 
       if ( (ssnow%GWwb(i) .lt. (soil%GWwatr(i)+volwatmin)) ) then
          xsi = ((soil%GWwatr(i)+volwatmin) - ssnow%GWwb(i)) * GWdzmm(i)  !mm
          ssnow%GWwb(i) = soil%GWwatr(i)+volwatmin
          ssnow%qhz(i) = ssnow%qhz(i) - xsi / dels
       end if

   end do

   do k=1,ms
      do i=1,mp
         
       !update mass variables
         ssnow%wmliq(i,k)      = ssnow%wbliq(i,k) * &
                                         soil%zse_vec(i,k)*r2_density_liq
         ssnow%wmice(i,k)      = ssnow%wbice(i,k) * &
                                         soil%zse_vec(i,k)*r2_density_ice
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
   USE cable_IO_vars_module, ONLY: wlogn

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
   REAL, DIMENSION(mp,ms) :: tgg_old
   REAL, DIMENSION(mp) :: tggsn_old,wbtot_ic,del_wbtot
   REAL(r_2), DIMENSION(mp) :: xx
   real(r_2), dimension(mp,ms) :: gammzz_snow
   REAL                :: zsetot
   INTEGER, SAVE :: ktau =0 
   REAL(r_2) :: wb_lake_T, rnof2_T
   LOGICAL :: use_sli
   LOGICAL, SAVE :: first_gw_hydro_call = .true.

   use_sli = .false. 

   if (first_gw_hydro_call)   &
                     call set_r2_gw_params(C%density_ice,C%density_liq,&
                                 C%csice,C%cswat,C%Cgsnow)

   ktau = ktau +1 
   
   zsetot = sum(soil%zse) 
   ssnow%tggav = 0.

   DO k = 1, ms

      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot

   END DO


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
      soil%heat_cap_lower_limit(:,:) = 0.01  !never allow /0
   ELSE
      soil%heat_cap_lower_limit(:,:) = soil%css_vec(:,:) * soil%rhosoil_vec(:,:)
   END IF

!   IF( (.NOT.cable_user%cable_runtime_coupled ) .and. (first_gw_hydro_call)) THEN
!   
!         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
!!                                                  ! after discussion with BP
!!         ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
!!         ssnow%wbtot = 0.0
!!         ssnow%wb(:,:)  = MIN( soil%ssat_vec(:,:), MAX ( ssnow%wb(:,:), 0.5*soil%swilt_vec(:,:) ) )   
!!
!            
!            WHERE( ssnow%tgg(:,:) <= C%TFRZ)  &
!               ssnow%wbice(:,:) = ssnow%wb(:,:)*(1.0-exp(2.0*((ssnow%wb(:,:)/soil%ssat_vec(:,:))**2.0)*&
!                                  (ssnow%tgg(:,:)-C%TFRZ)))/exp(1.0-ssnow%wb(:,:)/soil%ssat_vec(:,:))
!!
!         WHERE ( soil%isoilm .eq. 9 .and. ssnow%snowd .le. 0.1*max_glacier_snowd) 
!!
!            ! permanent ice: fix hard-wired number in next version
!            ssnow%snowd = max_glacier_snowd
!            ssnow%osnowd = max_glacier_snowd
!            ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0
!
!         END WHERE
!!
!!         WHERE ( spread(soil%isoilm,2,ms) .eq. 9 )
!!
!!              ssnow%wb    = 0.95 * soil%ssat_vec
!!              ssnow%wbice = 0.95 * ssnow%wb
!!
!!         END WHERE
!!
!   END IF

   gammzz_snow(:,:) = 0._r_2
   do i=1,mp
      gammzz_snow(i,1) = real((1. - ssnow%isflag(i)) * C%cgsnow * ssnow%snowd(i),r_2)
   end do


   !Start with wb and wbice.  Need wbliq, wmliq,wmice,wmtot
   !find the mass of ice and liq from the prognostic volumetric values
   do k=1,ms
      do i=1,mp
         ssnow%wbliq(i,k) = ssnow%wb(i,k) - den_rat*ssnow%wbice(i,k)                !liquid volume
         ssnow%wmice(i,k) = ssnow%wbice(i,k)*r2_density_ice*soil%zse_vec(i,k) !ice mass
         ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*r2_density_liq*soil%zse_vec(i,k) !;liq mass
         ssnow%wmtot(i,k) = ssnow%wmice(i,k) + ssnow%wmliq(i,k)                  !liq+ice mass

         ssnow%wblf(i,k)   = max(ssnow%wbliq(i,k)/soil%ssat_vec(i,k),0.01_r_2)
         ssnow%wbfice(i,k) = max(ssnow%wbice(i,k)/soil%ssat_vec(i,k),0._r_2)

      end do
   end do

   IF( first_gw_hydro_call ) THEN
      do k=1,ms
        do i=1,mp
             ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),&
                                (1.0-soil%ssat_vec(i,k))*&
                                soil%css_vec(i,k) * soil%rhosoil_vec(i,k)*  &
                              soil%zse_vec(i,k) &
                  + ssnow%wmliq(i,k) * r2_cswat &
                  + ssnow%wmice(i,k) * r2_csice) &
                  + gammzz_snow(i,k)
      
        end do
      end do
   ENDIF  ! if(.NOT.cable_runtime_coupled) and first_gw_hydro_call

   !initial water in the soil column
   wbtot_ic(:)  = sum(ssnow%wmliq(:,:),2) + &
                  sum(ssnow%wmice(:,:),2) + &
                  ssnow%GWwb(:)*soil%GWdz(:)*r2_density_liq
              
   GWwb_ic(:) = ssnow%GWwb(:)


   gammzz_snow(:,:) = 0._r_2
   do i=1,mp
      gammzz_snow(i,1) = real((1. - ssnow%isflag(i)) * C%cgsnow * ssnow%snowd(i),r_2)
   end do

   do k=1,ms
     do i=1,mp
          ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),&
                             (1.0-soil%ssat_vec(i,k))*&
                             soil%css_vec(i,k) * soil%rhosoil_vec(i,k)*  &
                           soil%zse_vec(i,k) &
               + ssnow%wmliq(i,k) * r2_cswat &
               + ssnow%wmice(i,k) * r2_csice) &
               + gammzz_snow(i,k)
   
     end do
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

   ssnow%pudsto(:) = 0.0  !no puddle
   ssnow%smelt(:)  = ssnow%smelt(:)/dels    !change units to mm/s.  cable_driver then reverts back to mm
   ssnow%runoff(:) = (ssnow%rnof1(:) + ssnow%rnof2(:))!*dels  !cable_driver converts from mm/s to mm
                                                     !rnof1 and rnof2 are already in mm/s
   ! Set weighted soil/snow surface temperature
   ssnow%tss(:) =  (1-ssnow%isflag(:))*ssnow%tgg(:,1) + ssnow%isflag(:)*ssnow%tggsn(:,1)

   !total water mass at the end of the soilsnow_GW routine
   ssnow%wbtot(:)  = sum(ssnow%wmtot,dim=2) + ssnow%GWwb(:)*soil%GWdz(:)*r2_density_liq
                 

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
    REAL(r_2), dimension(0:ms)  :: zimm     !layer interface depth in mm  
    REAL(r_2), dimension(ms)    :: zmm      !node depths in mm
    REAL(r_2)                   :: tempi, temp0,voleq1,wbrat
    REAL(r_2), dimension(mp,ms+1) :: ice_correction

    INTEGER :: k,i

    !make code cleaner define these here
    dzmm    = 1000.0_r_2 * real(soil%zse(:),r_2)
    zimm(0) = 0._r_2

    do k=1,ms
       zimm(k) = zimm(k-1) + dzmm(k)
       zmm(k)  = zimm(k-1) + 0.5_r_2*dzmm(k)
    end do 

    do i=1,mp
       GWdzmm(i) = real(soil%GWdz(i),r_2)*1000._r_2
       GWzimm(i) = zimm(ms)+GWdzmm(i)
       zaq(i)    = zimm(ms) + 0.5_r_2*GWdzmm(i)
    end do

    !equilibrium water content
    do k=1,ms
       do i=1,mp

          if ((ssnow%wtd(i) .le. zimm(k-1))) then         !fully saturated

             ssnow%wbeq(i,k) = soil%ssat_vec(i,k)
             
          elseif ((ssnow%wtd(i) .le. zimm(k)) .and. &
                  (ssnow%wtd(i) .gt. zimm(k-1))) then

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
          else

             tempi = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(k))/&
                      soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
             temp0 = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(k-1))/&
                     soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))   
             ssnow%wbeq(i,k) = -soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
                               (1._r_2-1._r_2/soil%bch_vec(i,k))/&
                               (zimm(k)-zimm(k-1))*(tempi-temp0)+soil%watr(i,k)
          end if
          
          ssnow%wbeq(i,k) = min(max(ssnow%wbeq(i,k),soil%watr(i,k)),soil%ssat_vec(i,k))

          wbrat = min(max((&
                  ssnow%wbeq(i,k) - soil%watr(i,k))/(soil%ssat_vec(i,k)-soil%watr(i,k)),&
                          0.01_r_2),1._r_2)

          ssnow%zq(i,k) = max(&
                            -soil%sucs_vec(i,k)*(wbrat**(-soil%bch_vec(i,k))),sucmin)

       end do  !mp
    end do  !ms
 
    do i=1,mp
    !Aquifer Equilibrium water content
       if (ssnow%wtd(i) .le. zimm(ms)) then      !fully saturated

          ssnow%GWwbeq(i) = soil%GWssat_vec(i)

       elseif ((ssnow%wtd(i) .gt. GWzimm(i)))   then     !fully unsaturated

          tempi = &
                (((abs(soil%GWsucs_vec(i))+ssnow%wtd(i)-GWzimm(i))/&
                  abs(soil%GWsucs_vec(i))))**(1._r_2-1._r_2/soil%GWbch_vec(i))
          temp0 = (((abs(soil%GWsucs_vec(i))+ssnow%wtd(i)-zimm(ms))/&
                     abs(soil%GWsucs_vec(i))))**(1._r_2-1._r_2/soil%GWbch_vec(i))   
          ssnow%GWwbeq(i) = -abs(soil%GWsucs_vec(i))*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
                          (1._r_2-1._r_2/soil%GWbch_vec(i))/&
                           (GWzimm(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)   

       else           

          tempi  = 1._r_2
          temp0  = (((abs(soil%GWsucs_vec(i))+ssnow%wtd(i)-zimm(ms))/&
                     abs(soil%GWsucs_vec(i))))**(1._r_2-1._r_2/soil%GWbch_vec(i))               
          voleq1 = -abs(soil%GWsucs_vec(i))*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
                    (1._r_2-1._r_2/soil%GWbch_vec(i))/&
                    (ssnow%wtd(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)
          ssnow%GWwbeq(i) = (voleq1*(ssnow%wtd(i)-zimm(ms)) + &
                            (soil%GWssat_vec(i)-soil%GWwatr(i))*&
                         (GWzimm(i)-ssnow%wtd(i)))/(GWzimm(i)-zimm(ms)) + soil%GWwatr(i)

       end if



       wbrat = min(1._r_2,max(0.01_r_2, (ssnow%GWwbeq(i)-soil%GWwatr(i))/&
                                         (soil%GWssat_vec(i)-soil%GWwatr(i))))

       ssnow%GWzq(i) = -abs(soil%GWsucs_vec(i))*((wbrat)**(-soil%GWbch_vec(i)))
       ssnow%GWzq(i) = max(sucmin, ssnow%GWzq(i))
       
    end do

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
          dzmm_one  = 1000._r_2 * real(soil%zse_vec(i,1),r_2)
          icemass  = ssnow%wbice(i,1) * dzmm_one
          liqmass  = (ssnow%wb(i,1)-den_rat*ssnow%wbice(i,1)) * dzmm_one
          totmass  = max(liqmass+icemass,real(1e-2,r_2))
          icef(i)     = max(0._r_2,min(1._r_2, icemass / totmass))
      end do
   
   
      !srf frozen fraction.  should be based on topography
      do i = 1,mp
         fice = (exp(-gw_params%IceAlpha*(1._r_2-icef(i)))-&
                 exp(-gw_params%IceAlpha))/&
                 (1._r_2-exp(-gw_params%IceAlpha))
         fice = min(1._r_2,max(0._r_2,fice))
   
         satfrac_liqice(i) = fice + (1._r_2-fice)*ssnow%satfrac(i)
   
         wb_unsat = ((ssnow%wb(i,1)-den_rat*ssnow%wbice(i,1)) -&
                     ssnow%satfrac(i)*soil%ssat_vec(i,1))/(1.-ssnow%satfrac(i))
         wb_unsat = min(soil%ssat_vec(i,1),max(0.,wb_unsat))
   
         wb_evap_threshold = min( max( &
                             gw_params%SoilEvapAlpha*soil%sfc_vec(i,1), &
                             soil%swilt_vec(i,1) ), soil%ssat_vec(i,1) )
   
         !Sakguchi and Zeng 2009
         if (wb_unsat .ge. wb_evap_threshold) then
            xx = 1.
         else
            xx = 0.25 * (1._r_2 - cos(r2pi*wb_unsat/(wb_evap_threshold)))**2.0
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
                                 ( ssnow%wbice(i,1) / ssnow%wb(i,1) )**2 ) ) )
   
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
   REAL(r_2), DIMENSION(mp,1:ms+1) :: hk_ice_factor
   REAL(r_2), DIMENSION(mp,ms) :: smp_ice_min
   REAL(r_2), DIMENSION(0:ms) :: zimm  !depths at interface between layers
   REAL(r_2), pointer, dimension(:,:) ::wb_temp 

    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water (calculated using total (not liquid))

    do k=1,ms
       do i=1,mp
          ssnow%icefrac(i,k) = ssnow%wbice(i,k)/(max(ssnow%wb(i,k),0.01_r_2))
          ssnow%fracice(i,k) = (exp(-gw_params%IceBeta*(1._r_2-ssnow%icefrac(i,k)))&
                               -exp(-gw_params%IceBeta))/(1._r_2-exp(-gw_params%IceBeta))
       end do
    end do

    ssnow%fracice(:,:) = max( min( ssnow%fracice, 1._r_2), 0._r_2)
    smp_ice_min(:,:) = sucmin
    if (gw_params%ssgw_ice_switch) then
       do k=1,ms
          wb_temp => ssnow%wbliq(:,:)
          kk = min(k+1,ms)
          do i=1,mp
             hk_ice_factor(i,k) = 0.5*( (1.0 - ssnow%wbice(i,k))**gw_params%ice_impedence +&
                                        (1.0 - ssnow%wbice(i,kk))**gw_params%ice_impedence)
          end do
       end do
       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
    else
       wb_temp => ssnow%wb(:,:)
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
    !aquifer ice 

    !potential from soil water rention function
    !defined as layer average
    do k=1,ms 
       do i=1,mp
          s_mid(i) = (wb_temp(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
              (soil%ssat_vec(i,k)-soil%watr(i,k))

          s_mid(i) = min(max(s_mid(i),0.01_r_2),1._r_2)

          ssnow%smp(i,k) = -soil%sucs_vec(i,k)*s_mid(i)**(-soil%bch_vec(i,k))

          ssnow%smp(i,k) = max(ssnow%smp(i,k),sucmin)

          ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/s_mid(i)

       end do   
    end do

    !Aquifer potential
    do i=1,mp

       s_mid(i) = (ssnow%GWwb(i)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.01_r_2),1._r_2)

       ssnow%GWsmp(i)    = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))
       ssnow%GWsmp(i)    = max(ssnow%GWsmp(i),sucmin)
       ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/s_mid(i)

    end do

    !bottom layer and aquifer conductivity
    do i=1,mp
       s_mid(i) = 0.5*((liq_ratio(i)*ssnow%GWwb(i)-soil%GWwatr(i))/&
                    (soil%GWssat_vec(i)-soil%GWwatr(i)) + &
                   (wb_temp(i,ms)-soil%watr(i,ms))/&
                    (soil%ssat_vec(i,ms)-soil%watr(i,ms)))

       s_mid(i) = min(max(s_mid(i),0.01_r_2),1._r_2)
       s2(i)    = hk_ice_factor(i,ms+1)*0.5*&
                  (soil%hyds_vec(i,ms)+ soil%GWhyds_vec(i))*&
                    s_mid(i)**(2._r_2*soil%bch_vec(i,ms)+2._r_2)

       ssnow%GWhk(i)     =s_mid(i)*s2(i)
       ssnow%GWdhkdw(i)  =  (2._r_2*soil%bch_vec(i,ms)+3._r_2)*&
                           s2(i)/(soil%GWssat_vec(i)+soil%ssat_vec(i,ms)&
                                  -soil%GWwatr(i)-soil%watr(i,ms))

       s2(i)    = hk_ice_factor(i,ms)*&
                  soil%hyds_vec(i,ms)*&
                    s_mid(i)**(2._r_2*soil%bch_vec(i,ms)+2._r_2)
      
       ssnow%hk(i,ms)     = s_mid(i)*s2(i)
       ssnow%dhkdw(i,ms)  =  (2._r_2*soil%bch_vec(i,ms)+3._r_2)*&
                           s2(i)/(soil%GWssat_vec(i)+soil%ssat_vec(i,ms)&
                                  -soil%GWwatr(i)-soil%watr(i,ms))

    end do
    !hydraulic conductivity
    !Interfacial so uses layer i and i+1
    do k=1,ms-1
       kk=k+1
       do i=1,mp

          s1(i) = 0.5_r_2*((wb_temp(i,k)-soil%watr(i,k)) + &
                           (wb_temp(i,kk)-soil%watr(i,kk))) / &
                         (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + &
                         (soil%ssat_vec(i,kk)-soil%watr(i,kk))))

          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
          s2(i) = hk_ice_factor(i,k)*soil%hyds_vec(i,k)*&
                     s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)

          ssnow%hk(i,k)    =  s1(i)*s2(i)
          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*s2(i)/&
                            (soil%ssat_vec(i,k)-soil%watr(i,k)+ &
                             soil%ssat_vec(i,kk)-soil%watr(i,kk)) 
       end do
    end do

END SUBROUTINE calc_soil_hydraulic_props


  SUBROUTINE aquifer_recharge(dt,ssnow,soil,veg,zaq,zmm,dzmm)
  USE cable_common_module

  IMPLICIT NONE
    real, intent(in) :: dt 
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    REAL(r_2), dimension(:), intent(in)       :: zaq
    REAL(r_2), dimension(:), intent(in)       :: zmm,dzmm

    integer :: i    

    !Doing the recharge outside of the soln of Richards Equation makes it easier to track total recharge amount.
    !Add to ssnow at some point 
    do i=1,mp
       if ((ssnow%wtd(i) .le. (sum(dzmm(1:ms-1),dim=1)+0.75*dzmm(ms)) ) .or. &
           (veg%iveg(i) .eq. 16) ) then

          ssnow%Qrecharge(i) = 0._r_2
       else
          ssnow%Qrecharge(i) = -ssnow%hk(i,ms)*&
                               ((ssnow%GWsmp(i)-ssnow%smp(i,ms)) -&
                            real(gw_params%GWmacro_down,r_2)*&
                               (ssnow%GWzq(i)-ssnow%zq(i,ms))) / &
                                 (500._r_2*(soil%GWdz(i)+soil%zse_vec(i,ms)) )
       end if
    end do


  END SUBROUTINE aquifer_recharge

  SUBROUTINE subsurface_drainage(ssnow,soil,veg,dzmm)
  USE cable_common_module

  IMPLICIT NONE
  
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    REAL(r_2), dimension(:), intent(in)       :: dzmm
    REAL(r_2), dimension(mp)                  :: sm_tot,drain_dens_fac
    INTEGER, dimension(mp)                    :: k_drain
    integer :: i,k

    real(r_2), dimension(17) :: Efold_mod


    do i=1,mp
       drain_dens_fac(i) = min(4.0,max(1.0,0.001_r_2/soil%drain_dens(i)))
    end do
    Efold_mod(:) = 1.0
    !Efold_mod(1:4) = (/0.2,0.2,0.2,0.2/)
    !Efold_mod(9) = 0.25
    do i=1,mp

       !Note: future revision will have interaction with river here. nned to
       !work on router and add river type cells
       ssnow%qhz(i)  = min(max(soil%slope(i),0.00001),0.1)*&
                       gw_params%MaxHorzDrainRate* &
                        exp(-ssnow%wtd(i)/(1000._r_2*&
                       (gw_params%EfoldHorzDrainRate+drain_dens_fac(i))))

       !drain from sat layers
       k_drain(i) = ms+1
       do k=ms,2,-1
          if (ssnow%wtd(i) .le. sum(dzmm(1:k),dim=1)) then
             k_drain(i) = k + 1
          end if
       end do
       k_drain(i) = max(k_drain(i),3)

   end do

   do i=1,mp
       ssnow%qhlev(i,:) = 0._r_2
       sm_tot(i) = 0._r_2
       ssnow%qhlev(i,:) = 0._r_2

       sm_tot(i) = max((ssnow%GWwb(i) - soil%watr(i,ms))*&
                       (1._r_2-ssnow%fracice(i,ms)), 0.)

       do k=k_drain(i),ms
          sm_tot(i) = sm_tot(i) +max(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)
       end do

       if (sm_tot(i) .ge. 1.0e-12) then
           ssnow%qhz(i) = ssnow%qhz(i) / sm_tot(i)
           do k=k_drain(i),ms
              ssnow%qhlev(i,k) = ssnow%qhz(i)
           end do
           ssnow%qhlev(i,ms+1) = max((ssnow%GWwb(i) - soil%watr(i,ms))*&
                                 (1._r_2-ssnow%fracice(i,ms)), 0.)*ssnow%qhz(i)
       endif

       !incase every layer is frozen very dry
       !ssnow%qhz(i) = sum(ssnow%qhlev(i,:),dim=1)

       !Keep "lakes" saturated forcing qhz = 0.  runoff only from lakes
       !overflowing
       if (veg%iveg(i) .eq. 16) then
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

     !if !gw_model and !or_evap:
                   !in cable_um_init_subrs.F90 or cable_parameters:
                   !  ssat_vec(i,:) = ssat
                   !  zse_vec(i,:)  = zse
                   !  UM: slope_std read in
                   ! offline: slope_std read in or set to const
                   !  all gw_params set by default in cable_common
                   ! doesn do anything but cannot hurt


     S(:) = 0._r_2
     do k=1,gw_params%level_for_satfrac
       S(:) = S(:) + max(0.01,min(1.0, &
              (ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/&
               max(0.001,soil%ssat_vec(:,k)-soil%watr(:,k)) ) )*soil%zse_vec(:,k)
     end do
     S(:) = S(:)/sum(soil%zse(1:gw_params%level_for_satfrac),dim=1)
     !srf frozen fraction.  should be based on topography
      do i = 1,mp
         !Saturated fraction
          if (gw_params%MaxSatFraction .gt. 1e-7 .and. veg%iveg(i) .lt. 16) then 
             slopeSTDmm = sqrt(min(max(&
                           gw_params%MaxSatFraction*soil%slope_std(i),&
                           1e-5),10000._r_2)) ! ensure some variability
             ssnow%satfrac(i)    = max(0._r_2,min(0.99_r_2,&
                 !note UM wants std03, and erf is not included then
                                   1._r_2 - my_erf( slopeSTDmm / sqrt(2.0* S(i)) ) ) )  
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

    !if gw_model = true 
                   !cable_um_init_subrs.F90 or cable_parameters:
                           ! ssat(i) = ssat_vec(i,1)
    !if gw_model = false 
                   !cable_um_init_subrs.F90 or cable_parameters:
                            !ssat_vec(i,:) = ssat(i)
                            !so ssat_vec can be used although soilsnow uses ssat

    do i=1,mp
       if (veg%iveg(i) .lt. 16 .and. soil%isoilm(i) .ne. 9 .and. &
           ssnow%snowd(i) .le. 1e-8 ) then

          unsat_wb(i) =  (ssnow%wb(i,1) - soil%ssat_vec(i,1)*&
                      min(0.95,max(0.0,ssnow%satfrac(i))))/(1.0 - min(0.95,max(0.0,ssnow%satfrac(i))))

          unsat_wb(i) = max(soil%watr(i,1)+1e-2, min(soil%ssat_vec(i,1), unsat_wb(i) ) )

          unsat_smp(i) = sign(soil%sucs_vec(i,1),-1.0) * min(1.0, &
                         (max(0.001, (unsat_wb(i)-soil%watr(i,1))/(soil%ssat_vec(i,1)-&
                         soil%watr(i,1)) ) )** (-soil%bch_vec(i,1)) )

          unsat_smp(i) = max(-1.0e7,unsat_smp(i) )/1000._r_2 !m

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

    call iterative_wtd (ssnow, soil, veg, cable_user%test_new_gw)
 
    CALL calc_soil_hydraulic_props(ssnow,soil,veg)
 
    call  ovrlndflx (dels, ssnow, soil, veg,canopy,sli_call )
 
    dzmm = real(soil%zse(:),r_2)*1000._r_2
 
    CALL subsurface_drainage(ssnow,soil,veg,dzmm)
 
    zmm(:) = 1000._r_2*(sum(real(soil%zse,r_2),dim=1))
    zaq(:) = zmm(:) + 0.5_r_2*soil%GWdz(:)*1000._r_2
 
    call aquifer_recharge(dels,ssnow,soil,veg,zaq,zmm,zmm)

  END SUBROUTINE sli_hydrology


  SUBROUTINE set_unsed_gw_vars(ssnow,soil,canopy)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil ! soil parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy
    
    !CALL point2constants( C )

    CALL set_r2_gw_params(C%density_ice,C%density_liq,&
                              C%csice,C%cswat,C%cgsnow)

       ssnow%qhlev = 0.
       ssnow%Qrecharge = 0.
       ssnow%fwtop = 0.
       ssnow%wtd = 0.
       ssnow%satfrac = 1.0
       ssnow%qhz = 0.
       ssnow%wbliq = ssnow%wb - den_rat*ssnow%wbice
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


  subroutine set_r2_gw_params(C_density_ice,C_density_liq,&
                              C_csice,C_cswat,C_cgsnow)
     use cable_common_module, only: gw_params
     use cable_def_types_mod, only: r_2

     real, intent(in) :: C_density_ice,&
                         C_density_liq,&
                         C_csice,&
                         C_cswat,&
                         C_cgsnow


     !density of ice / density of liq (r_2 type)
     r2_density_ice = real(C_density_ice,r_2)
     r2_density_liq = real(C_density_liq,r_2)

     den_rat       = r2_density_ice/r2_density_liq
     r2_cgsnow     = real(C_cgsnow,r_2)
     r2_csice      = real(C_csice,r_2)
     r2_cswat      = real(C_cswat,r_2)

     r2_cs_rho_ice = r2_csice*r2_density_ice
     r2_cs_rho_wat = r2_cswat*r2_density_liq
     r2pi          = 3.14159_r_2

  end subroutine set_r2_gw_params

  function set_wbliq(wb,wbice)
      real(r_2), dimension(mp,ms),intent(in) :: wb,wbice

      real(r_2), dimension(mp,ms) :: set_wbliq

      if (cable_user%gw_model) then
          set_wbliq(:,:) = wb(:,:) - wbice(:,:)*(real(C%density_ice,r_2)/&
                                                real(C%density_liq,r_2))
      else
          set_wbliq(:,:) = wb(:,:) - wbice(:,:)

      end if

  end function


END MODULE cable_gw_hydro_module
