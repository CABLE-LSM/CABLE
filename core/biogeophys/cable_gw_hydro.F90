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
! Purpose: large scale hydrological paraeterizations of subgrid scale generation
! of surface and subsurface runoff due to subgrid variations in soil moisture
! and an unconfied aquifer model

! Contact: m.decker@unsw.edu.au
!
!
! ==============================================================================

MODULE cable_gw_hydro_module

   USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             r_2, ms, mp

   USE cable_data_module, ONLY : C=>PHYS

   USE cable_common_module, ONLY : cable_user, cable_runtime,gw_params,  &
                                   cgsnow,csice,cswat,nglacier,snmin,    &
                                   max_ssdn,max_sconds,max_glacier_snowd,&
                                   knode_gl

   USE cable_soil_snow_module, only :  snowdensity, snow_melting, snowcheck, &
                                      snowl_adjust,snow_accum, stempv,trimb, &
                                      glacier_adjustment

   IMPLICIT NONE

   REAL(r_2), PARAMETER :: sucmin  = -1.0e8, &! minimum soil pressure head [mm]
                      volwatmin    = 1e-4,         &! min soil water [mm]      
                      wtd_uncert   = 0.1,          &! uncertaintiy in wtd calcultations [mm]
                      wtd_max      = 1000000.0,    &! maximum wtd [mm]
                      wtd_min      = 100.0          ! minimum wtd [mm]

   INTEGER, PARAMETER :: wtd_iter_max = 20 ! maximum number of iterations to find the water table depth                    

  ! ! This module contains the following subroutines:
   PUBLIC soil_snow_gw,calc_srf_wet_fraction,sli_hydrology ! must be available outside this module
   ! This module contains the following subroutines:
   PRIVATE subsurface_drainage,iterative_wtd,ovrlndflx,aquifer_recharge, calc_soil_hydraulic_props
   PRIVATE calc_equilibrium_water_content,GWsoilfreeze, remove_transGW
   PRIVATE smoistgw,my_erf

CONTAINS


SUBROUTINE GWsoilfreeze(dels, soil, ssnow)
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(IN)    :: soil
   REAL(r_2), DIMENSION(mp)           :: sicefreeze
   REAL(r_2), DIMENSION(mp)           :: sicemelt
   REAL, DIMENSION(mp)                :: xx, ice_mass,liq_mass,tot_mass
   INTEGER :: i,j,k
   REAL(r_2),DIMENSION(mp,ms) :: var_frz_lim,iceF  !Decker and Zeng 2009

   do k=1,ms
   do i=1,mp
      if (soil%isoilm(i) .ne. 9) then
         if  (ssnow%tgg(i,k) .le. C%TFRZ) then
            var_frz_lim(i,k) = (1. - exp(-2.*(ssnow%wb(i,k)/soil%ssat_vec(i,k))**4.0 *&
                             (ssnow%tgg(i,k)-C%TFRZ)))/exp(1. - ssnow%wb(i,k)/soil%ssat_vec(i,k))
            var_frz_lim(i,k) = max(0.5,var_frz_lim(i,k))
         else
            var_frz_lim(i,k) = 0.
         endif
      else
         var_frz_lim(i,k) = 0.95
      end if
   end do
   end do

   xx = 0.
   DO k = 1, ms
   DO i=1,mp

      ice_mass = ssnow%wbice(i,k)*real(soil%zse(k)*C%density_ice,r_2)
      liq_mass = ssnow%wbliq(i,k)*real(soil%zse(k)*C%density_liq,r_2)
      tot_mass = liq_mass + ice_mass
      
      if (ssnow%tgg(i,k) .lt. C%TFRZ .and. &
           var_frz_lim(i,k) * ssnow%wb(i,k) - ssnow%wbice(i,k) > .001) then
         
         sicefreeze(i) = MIN( MAX( 0.0_r_2, ( var_frz_lim(i,k) * ssnow%wb(i,k) -      &
                      ssnow%wbice(i,k) ) ) * soil%zse(k) * C%density_ice,             &
                      ( C%TFRZ - ssnow%tgg(i,k) ) * ssnow%gammzz(i,k) / C%HLF )
         ssnow%wbice(i,k) = MIN( ssnow%wbice(i,k) + sicefreeze(i) / (soil%zse(k)  &
                            * 1000.0), var_frz_lim(i,k) * ssnow%wb(i,k) )
         xx(i) = soil%css(i) * soil%rhosoil(i)
         ssnow%gammzz(i,k) = MAX(                                              &
             REAL((1.0 - soil%ssat_vec(i,k)) * soil%css(i) * soil%rhosoil(i) ,r_2)            &
             + (ssnow%wb(i,k) - ssnow%wbice(i,k)) * REAL(cswat * C%density_liq,r_2)   &
             + ssnow%wbice(i,k) * REAL(csice * C%density_ice,r_2),              &
             REAL(xx(i),r_2)) * REAL( soil%zse(k),r_2 )

        if (k .eq. 1 .and. ssnow%isflag(i) .eq. 0) then
           ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + cgsnow * ssnow%snowd(i)
        end if

         ssnow%tgg(i,k) = ssnow%tgg(i,k) + REAL(sicefreeze(i))                    &
                             * C%HLF / REAL(ssnow%gammzz(i,k) )

     
      ELSEIF( ssnow%tgg(i,k) > C%TFRZ .AND. ssnow%wbice(i,k) > 0. ) then
         
         sicemelt(i) = MIN( ssnow%wbice(i,k) * soil%zse(k) * C%density_ice,              &
                    ( ssnow%tgg(i,k) - C%TFRZ ) * ssnow%gammzz(i,k) / C%HLF )
         
         ssnow%wbice(i,k) = MAX( 0.0_r_2, ssnow%wbice(i,k) - sicemelt(i)          &
                            / (soil%zse(k) * C%density_ice) )
         xx(i) = soil%css(i) * soil%rhosoil(i)
         ssnow%gammzz(i,k) = MAX(                                              &
              REAL((1.0-soil%ssat_vec(i,k)) * soil%css(i) * soil%rhosoil(i),r_2)             &
              + (ssnow%wb(i,k) - ssnow%wbice(i,k)) * REAL(cswat*C%density_liq,r_2)   &
              + ssnow%wbice(i,k) * REAL(csice * C%density_ice,r_2),            &
              REAL(xx(i),r_2) ) * REAL(soil%zse(k),r_2)
         if (k .eq. 1 .and. ssnow%isflag(i) .eq. 0) then
           ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + cgsnow * ssnow%snowd(i)
         end if
         ssnow%tgg(i,k) = ssnow%tgg(i,k) - REAL(sicemelt(i))                     &
                          * C%HLF / REAL(ssnow%gammzz(i,k))
       
      END IF
      !update the liq and ice volume and mass
      ice_mass(i)   = ssnow%wbice(i,k)*real(soil%zse(k)*C%density_ice,r_2)
      liq_mass(i)   = tot_mass(i) - ice_mass(i)
      ssnow%wbliq(i,k) = liq_mass(i) / real(soil%zse(k)*C%density_liq,r_2)
      ssnow%wbice(i,k) = ice_mass(i) / real(soil%zse(k)*C%density_ice,r_2)
      ssnow%wb(i,k)    = ssnow%wbliq(i,k) + ssnow%wbice(i,k)
    
   END DO
   END DO

END SUBROUTINE GWsoilfreeze
!
!!-----------------------------------------------------------------------------------------
!
!! -----------------------------------------------------------------------------
!
SUBROUTINE remove_transGW(dels, soil, ssnow, canopy, veg)
   !note same as remove_trans from soilsnow, except for
   !uses wbliq instead of wb, swilt_vec, and C%density_liq
   !swilt_vec, prevents it from being bitwise reproducable to previous version

   ! Removes transpiration water from soil.
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(canopy_type), INTENT(INOUT)         :: canopy
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(IN)    :: soil
   TYPE(veg_parameter_type), INTENT(IN)  :: veg
   REAL(r_2), DIMENSION(mp,0:ms+1) :: diff 
   REAL(r_2), DIMENSION(mp)      :: xx,xxd,evap_cur
   REAL(r_2), DIMENSION(mp,ms) :: zse_mp_mm
   INTEGER :: k,i

   do k=1,ms
      do i=1,mp
         zse_mp_mm(i,k)  = real(soil%zse(k)*C%density_liq,r_2)
      end do
   end do

   IF (cable_user%FWSOIL_switch.eq.'Haverd2013') THEN

     WHERE (canopy%fevc .lt. 0.0_r_2)
        canopy%fevw = canopy%fevw+canopy%fevc
        canopy%fevc = 0.0_r_2
     END WHERE
     DO k = 1,ms 
        ssnow%wbliq(:,k) = ssnow%wbliq(:,k) - ssnow%evapfbl(:,k)/zse_mp_mm(:,k)
     ENDDO

   ELSE
 
      xx(:) = 0._r_2
      xxd(:) = 0._r_2
      diff(:,:) = 0._r_2
   
      DO k = 1,ms  
   
         DO i=1,mp
   
            if (canopy%fevc(i) .gt. 0._r_2) then
   
               xx(i) = canopy%fevc(i) * dels / C%HL * veg%froot(i,k) + diff(i,k-1)
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

  ENDIF

  ssnow%wb = ssnow%wb + ssnow%wbice


END SUBROUTINE remove_transGW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!MD GW code from here on!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  !-------------------------------------------------------------------------
  SUBROUTINE ovrlndflx (dels, ssnow, soil,veg, canopy,sli_call )

  IMPLICIT NONE
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters
    TYPE(veg_parameter_type) , INTENT(IN)    :: veg  ! veg parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy
    LOGICAL, INTENT(IN)                      :: sli_call
    INTEGER                                  :: k, i, j
    REAL(r_2), DIMENSION(mp)           :: icef,efpor               !tmp vars, fraction of ice in gridcell
    REAL(r_2)                          :: tmpa,tmpb,qinmax         !tmp vars, maximum infiltration [mm/s]
    REAL(r_2), DIMENSION(mp)           :: satfrac_liqice,S       !saturated fraction of cell, wtd in m
    REAL(r_2)                          :: liqmass,icemass,totmass  !liquid mass,ice mass, total mass [mm]
    REAL(r_2)                          :: fice
    REAL(r_2)                          :: dzmm,slopeSTDmm

   !For now assume there is no puddle
   dzmm = real( C%density_liq * soil%zse(1), r_2)

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
      efpor(i) = max(0.001_r_2, soil%ssat_vec(i,1) - ssnow%wbice(i,1))
      icemass  = ssnow%wbice(i,1) * dzmm * C%density_ice/C%density_liq
      liqmass  = (ssnow%wb(i,1)-ssnow%wbice(i,1)) * dzmm
      totmass  = max(liqmass+icemass,real(1e-2,r_2))
      icef(i)     = max(0._r_2,min(1._r_2,gw_params%IceBeta*icemass / totmass))
   end do

   !sat fraction assuming topo controlled subgrid soil moisture distribution
   call saturated_fraction(ssnow,soil,veg)

   !srf frozen fraction.  should be based on topography
   do i = 1,mp
      fice = (exp(-gw_params%IceAlpha*(1._r_2-icef(i)))-exp(-gw_params%IceAlpha))/(1._r_2-exp(-gw_params%IceAlpha))
      fice  = min(max(fice,0._r_2),1._r_2)
      if (soil%isoilm(i) .ne. 9) then
         satfrac_liqice(i)   = max(0.,min(0.95,fice + (1._r_2-fice)*ssnow%satfrac(i) ) )
      else
         satfrac_liqice(i)   = max(0.,min(0.95,fice ) )
      end if

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

  if (nglacier .eq. 2) then
     call glacier_adjustment(dels,ssnow)
  end if

           
  END SUBROUTINE ovrlndflx


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
  TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
  TYPE (veg_parameter_type), INTENT(IN)     :: veg
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
  dzmm_mp  = real(spread((soil%zse(:)) * 1000.0,1,mp),r_2)    !layer thickness mm
  zimm(0)  = 0.0_r_2                                          !depth of layer interfaces mm

  !total depth of soil column
  do k=1,ms
    zimm(k) = zimm(k-1) + soil%zse(k)*1000._r_2
  end do

  def(:) = 0._r_2

  if (include_aquifer) then  !do we include the aquifer in the calculation of wtd?

     do i=1,mp
        total_depth_column(i) = zimm(ms) + soil%GWdz(i)*1000._r_2
        def(i) = def(i) + max(0._r_2,soil%GWssat_vec(i)-ssnow%GWwb(i))*soil%GWdz(i)*1000._r_2
     end do   

  end if

  !comute the total mass away from full saturation
  do k=1,ms
     do i=1,mp

       def(i) = def(i) +                                                           &
                max(0._r_2,(soil%ssat_vec(i,k)-(ssnow%wbliq(i,k)+&
                   C%density_ice/C%density_liq*ssnow%wbice(i,k)))*dzmm_mp(i,k))
      end do  !mp
  end do  !ms

  !find the deficit if the water table is at the bottom of the soil column
  do i=1,mp
     defc(i) = (soil%ssat_vec(i,ms))*(total_depth_column(i)+Nsucs_vec(i)/(1._r_2-invB(i))*            &
             (1._r_2-((Nsucs_vec(i)+total_depth_column(i))/Nsucs_vec(i))**(1._r_2-invB(i)))) 
     defc(i) = max(0.1_r_2,defc(i)) 

     !initial guess at wtd
     ssnow%wtd(:) = total_depth_column(:)*def(:)/defc(:)
  end do


 !use newtons method to solve for wtd, note this assumes homogenous column but
 !that is ok 
  do i=1,mp
    if ((soil%isoilm(i) .ne. 9) .and. (veg%iveg(i) .lt. 16)) then      

      if (defc(i) > def(i)) then                 !iterate tfor wtd

        jlp=0

        mainloop: DO

          tempa   = 1.0_r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-invB(i))
          derv    = (soil%ssat_vec(i,ms))*(tempa-tempb) + &
                                       soil%ssat_vec(i,ms)

          if (abs(derv) .lt. real(1e-12,r_2)) derv = sign(real(1e-12,r_2),derv)

          tempa   = 1.0_r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-invB(i))
          deffunc = (soil%ssat_vec(i,ms))*(ssnow%wtd(i) +&
                           Nsucs_vec(i)/(1-invB(i))* &
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
          tempa    = (abs(tmpc/Nsucs_vec(i)))**(-invB(i))
          tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-invB(i))
          derv     = (soil%ssat_vec(i,ms))*(tempa-tempb)
          if (abs(derv) .lt. real(1e-12,r_2)) derv = sign(real(1e-12,r_2),derv)

          tempa    = (abs((Nsucs_vec(i)+ssnow%wtd(i)-total_depth_column(i))/Nsucs_vec(i)))**(1._r_2-invB(i))
          tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-invB(i))
          deffunc  = (soil%ssat_vec(i,ms))*(total_depth_column(i) +&
                     Nsucs_vec(i)/(1._r_2-invB(i))*(tempa-tempb))-def(i)
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
     ssnow%wtd(i) = min(wtd_max,max(wtd_min,ssnow%wtd(i) ) )
  end do


  END SUBROUTINE iterative_wtd

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoistgw (fwtop,dt,ssnow,soil,prin)
  ! solves the modified richards equation (Zeng and Decker 2009) to find
  ! vertical mocement of soil water.  Bottom boundary condition is determined
  ! using a single layer groundwater module
  !
  SUBROUTINE smoistgw (dels,ssnow,soil,veg,canopy)

  IMPLICIT NONE
  
    REAL, INTENT(IN)                          :: dels  ! time step size (s)
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)     :: veg
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
    real(r_2) :: needed_space
    !MD DEBUG VARS
    INTEGER :: imp,ims,k_drain

    !make code cleaner define these here
    dzmm(:) = 1000.0_r_2 * real(soil%zse(:),r_2)
    do k=1,ms
       do i=1,mp
          dzmm_mp(i,k) = dzmm(k)
       end do
    end do

    zimm(0) = 0._r_2
    do k=1,ms
       zimm(k) = zimm(k-1) + dzmm(k)
    end do

    do k=1,ms
       zmm(k) = zimm(k-1) + 0.5_r_2*dzmm(k)
    end do

    do i=1,mp
       GWdzmm(i) = real(soil%GWdz(i),r_2)*1000._r_2
       GWzimm(i) = zimm(ms)+GWdzmm(i)
       zaq(i)    = zimm(ms) + 0.5_r_2*GWdzmm(i)
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
    
    !equilibrium water content
    CALL calc_equilibrium_water_content(ssnow,soil)

    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    CALL subsurface_drainage(ssnow,soil,veg,dzmm)

    qin(:) = 0._r_2
    den(:) = 0._r_2
    dne(:) = 0._r_2
    qout(:) = 0._r_2
    dqodw1(:) = 0._r_2
    dqodw2(:) = 0._r_2
    dqidw0(:) = 0._r_2
    dqidw1(:) = 0._r_2

    k = 1     !top soil layer
    do i=1,mp
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
    end do

    do k = 2, ms - 1     !middle soil layers
       do i=1,mp
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
       end do
    end do
       
    k = ms   !Bottom soil layer
    do i=1,mp
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
    end do
       
    CALL aquifer_recharge(dels,ssnow,soil,veg,zaq,zmm,dzmm)

    CALL trimb(at,bt,ct,rt,ms)                       !use the defulat cable tridiag solution

    do k=1,ms
       do i=1,mp
          ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + rt(i,k) - ssnow%qhlev(i,k)*dels/dzmm(k)   !volutermic liquid
       end do
    end do

    do i=1,mp
       ssnow%wbliq(i,ms) = ssnow%wbliq(i,ms) - ssnow%Qrecharge(i)*dels/dzmm(ms)
    end do
    do i=1,mp
       ssnow%GWwb(i) = ssnow%GWwb(i)  +  (ssnow%Qrecharge(i)-ssnow%qhlev(i,ms+1))*dels/GWdzmm(i)
    end do

    !determine the available pore space
    !volumetric
    do k=1,ms
       do i=1,mp
          eff_por(i,k)  = soil%ssat_vec(i,k) - ssnow%wbice(i,k)
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
          xsi = 0._r_2
       end if

       do k = 1,ms
          xsi = 0._r_2             !should be a single float (array not needed)
          if (ssnow%wbliq(i,k) .lt. volwatmin) then
             xsi = (volwatmin - ssnow%wbliq(i,k))*dzmm(k)  !in mm
             ssnow%wbliq(i,k) = volwatmin
             if (k .lt. ms) then
                ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) - xsi/dzmm(k+1)
             else
                ssnow%GWwb(i) = ssnow%GWwb(i) - xsi / GWdzmm(i)
             end if
          end if
       end do  !ms loop
 
       if ( (ssnow%GWwb(i) .lt. volwatmin) .and. (soil%isoilm(i) .ne. 9) ) then
          xsi = (volwatmin - ssnow%GWwb(i)) / GWdzmm(i)  !mm
          ssnow%GWwb(i) = volwatmin
          ssnow%qhz(i) = ssnow%qhz(i) - xsi / dels
       end if

   end do

   do k=1,ms
      do i=1,mp
        ssnow%wb(i,k)  = ssnow%wbliq(i,k) + ssnow%wbice(i,k)
      end do
   end do

   do i=1,mp
       ssnow%rnof2(i)        = ssnow%qhz(i)               !rnof2 is default cable deep runoff var 
   end do  !mp loop
       

 END SUBROUTINE smoistgw



! Inputs:
!	 dt_in - time step in sec
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
SUBROUTINE soil_snow_gw(dels, soil, ssnow, canopy, met, veg)

   REAL                     , INTENT(IN)     :: dels ! integration time step (s)
   TYPE(soil_parameter_type), INTENT(INOUT)  :: soil
   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
   TYPE(canopy_type)        , INTENT(INOUT)  :: canopy
   TYPE(veg_parameter_type) , INTENT(INOUT)  :: veg
   TYPE(met_type)           , INTENT(INOUT)  :: met ! all met forcing

   INTEGER             :: k,i
   REAL, DIMENSION(mp) :: snowmlt
   REAL, DIMENSION(mp) :: GWwb_ic
   REAL, DIMENSION(mp) :: tgg_old, tggsn_old,wbtot_ic,del_wbtot
   REAL(r_2), DIMENSION(mp) :: xx,xxx
   REAL                :: zsetot
   LOGICAL :: use_sli
   LOGICAL, SAVE :: first_gw_hydro_call = .true.
   !character(len=4) :: cnode
   !character(len=4) :: c_step
   !integer :: level_num
   integer, save :: call_num = 0
  
   use_sli = .false. 

   call_num = call_num + 1

   !write(cnode,FMT='(I4.4)') int(knode_gl,4)
   !write(c_step,FMT='(I4.4)') int(call_num,4)
   !level_num = 1031 + knode_gl 
   !open (level_num,file="/g/data1/w35/mrd561/my_ouput/681/gw_hydro_node-"//cnode//"_step"//c_step,status="replace")

   zsetot = sum(soil%zse) 
   ssnow%tggav = 0.

   DO k = 1, ms

      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot

   END DO


   IF( cable_runtime%offline .or. cable_runtime%mk3l ) ssnow%t_snwlr = 0.05_r_2

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
   end do

   IF( (.NOT.cable_user%cable_runtime_coupled ) .and. (call_num <= 1)) THEN
   
         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
                                                  ! after discussion with BP
         ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
         ssnow%wbtot = 0.0
         ssnow%wb(:,:)  = MIN( soil%ssat_vec(:,:), MAX ( ssnow%wb(:,:), soil%swilt_vec(:,:) ) )   

         DO k = 1, ms
            
            WHERE( ssnow%tgg(:,k) <= C%TFRZ .AND. ssnow%wbice(:,k) <= 0.01 )   &
               ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)

            WHERE( ssnow%tgg(:,k) < C%TFRZ)                                    &
               ssnow%wbice(:,k) = 0.8 * ssnow%wb(:,k)
            
         END DO

         WHERE ( soil%isoilm .eq. 9 ) 

            ! permanent ice: fix hard-wired number in next version
            ssnow%snowd = max_glacier_snowd
            ssnow%osnowd = max_glacier_snowd
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0

         END WHERE

         WHERE ( spread(soil%isoilm,2,ms) .eq. 9 )

              ssnow%wb    = 0.95 * soil%ssat_vec
              ssnow%wbice = 0.95 * ssnow%wb

         END WHERE
         
         xx=soil%css * soil%rhosoil

         ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat_vec(:,1)) * soil%css * soil%rhosoil &
              & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * C%density_liq &
              & + ssnow%wbice(:,1) * csice * C%density_ice, xx ) * soil%zse(1)

   ENDIF  ! if(.NOT.cable_runtime_coupled) and first_gw_hydro_call

   do k=1,ms
      do i=1,mp
         ssnow%wbliq(i,k) = ssnow%wb(i,k) - ssnow%wbice(i,k)                     !liquid volume
         ssnow%wmice(i,k) = ssnow%wbice(i,k)*real(C%density_ice*soil%zse(k),r_2) !ice mass
         ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*real(C%density_liq*soil%zse(k),r_2) !liquid mass
         ssnow%wmtot(i,k) = ssnow%wmice(i,k) + ssnow%wmliq(i,k)                  !liq+ice mass
         ssnow%wblf(i,k)   = max(ssnow%wbliq(i,k)/soil%ssat_vec(i,k),0.01_r_2)
         ssnow%wbfice(i,k) = max(ssnow%wbice(i,k)/soil%ssat_vec(i,k),0._r_2)

      end do
   end do

   xx = soil%css * soil%rhosoil
   IF (call_num <= 0)                   &
     ssnow%gammzz(:,1) = MAX( &
                        (1.0 - soil%ssat_vec(:,1)) * soil%css(:) * soil%rhosoil(:)  &
         & + ssnow%wbliq(:,1) * cswat * C%density_liq           &
         & + ssnow%wbice(:,1) * csice * C%density_ice, xx(:) ) * soil%zse(1) +   &
         & (1. - ssnow%isflag(:)) * cgsnow * ssnow%snowd(:)

   !initial water in the soil column
   do i=1,mp
      wbtot_ic(i)  = sum(ssnow%wbliq(i,:)*C%density_liq*soil%zse(:),1) + &
                     sum(ssnow%wbice(i,:)*C%density_ice*soil%zse(:),1) + &
                     ssnow%GWwb(i)*soil%GWdz(i)*C%density_liq
   end do
               
   GWwb_ic(:) = ssnow%GWwb(:)

   !write(level_num,*) 'Before calling anything in ge hydro '
   !call write_stuff()

   CALL snowcheck (dels, ssnow, soil, met )

   CALL snowdensity (dels, ssnow, soil)

   CALL snow_accum (dels, canopy, met, ssnow, soil )

   CALL snow_melting (dels, snowmlt, ssnow, soil )

   !write(level_num,*) 'after snow melt '
   !call write_stuff()

   ! Add snow melt to global snow melt variable:
   ssnow%smelt(:) = snowmlt(:)
   ! Adjust levels in the snowpack due to snow accumulation/melting,
   ! snow aging etc...

   CALL snowl_adjust(dels, ssnow, canopy )

   CALL stempv(dels, canopy, ssnow, soil)

   !write(level_num,*) 'after stempv '
   !call write_stuff()

   !do the soil and snow melting, freezing prior to water movement
   ssnow%tss(:) =  (1-ssnow%isflag(:))*ssnow%tgg(:,1) + ssnow%isflag(:)*ssnow%tggsn(:,1)

   CALL snow_melting (dels, snowmlt, ssnow, soil )
   
   ! Add new snow melt to global snow melt variable: 
   ssnow%smelt(:) = ssnow%smelt(:) + snowmlt(:)

   !write(level_num,*) 'before removetrans '
   !call write_stuff()

   CALL remove_transGW(dels, soil, ssnow, canopy, veg)        !transpiration loss per soil layer

   CALL  GWsoilfreeze(dels, soil, ssnow)

   !write(level_num,*) 'after freeze '
   !call write_stuff()
   ssnow%fwtop = canopy%precis/dels + ssnow%smelt/dels   !water from canopy and snowmelt [mm/s]   

   CALL iterative_wtd (ssnow, soil, veg, .true. )  
   !write(level_num,*) 'after wtd '
   !call write_stuff()

   CALL ovrlndflx (dels, ssnow, soil, veg, canopy,use_sli )         !surface runoff, incorporate ssnow%pudsto?

   ssnow%sinfil = ssnow%fwtop - canopy%segg             !remove soil evap from throughfall

   !write(level_num,*) 'after overland '
   !call write_stuff()

   CALL smoistgw (dels,ssnow,soil,veg,canopy)               !vertical soil moisture movement. 
   !write(level_num,*) 'after smoistgw '
   !call write_stuff()
 
   ! correction required for energy balance in online simulations 
   IF( cable_runtime%um) THEN

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
 
   ssnow%smelt  = ssnow%smelt/dels    !change units to mm/s.  cable_driver then reverts back to mm
   ssnow%runoff = (ssnow%rnof1 + ssnow%rnof2)!*dels  !cable_driver converts from mm/s to mm
                                                     !rnof1 and rnof2 are already in mm/s
   ! Set weighted soil/snow surface temperature
   !write(level_num,*) 'max tgg ',maxval(ssnow%tgg,dim=1)
   !write(level_num,*) 'min tgg ',minval(ssnow%tgg,dim=1)
   !write(level_num,*) 'max tggsn ',maxval(ssnow%tggsn,dim=1)
   !write(level_num,*) 'min tggsn ',minval(ssnow%tggsn,dim=1)
   !write(level_num,*) 'max isflag ',maxval(ssnow%isflag,dim=1)
   !write(level_num,*) 'min isflag ',minval(ssnow%isflag,dim=1)

   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   !total water mass at the end of the soilsnow_GW routine
   do i=1,mp
      ssnow%wbtot(i)  = sum(ssnow%wbliq(i,:)*C%density_liq*soil%zse(:),1) + &
                        sum(ssnow%wbice(i,:)*C%density_ice*soil%zse(:),1) + &
                        ssnow%GWwb(i)*soil%GWdz(i)*C%density_liq
   end do
                 
   !for debug water balance.  del_wbtot = fluxes = infiltration [though-evap] - trans - qhorz drainage
   del_wbtot   = dels * (ssnow%sinfil - ssnow%rnof2 - canopy%fevc / C%HL)

   first_gw_hydro_call = .false.

    !do k=1,ms
    !do i=1,mp
    !   write(6,*) knode_gl,' 3108 ssnow%wb ',k,'-',i,'',ssnow%wb(i,k)
    !end do
    !end do
   !write(level_num,*) 'at end '
   !call write_stuff()

   !call flush(level_num)
   !close(level_num)
!below is simply to aid in debugging
!contains
!
!    subroutine write_more_stuff()
!
!    write(level_num,*) 'ssnow%tgg '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%tgg(i,1:6)
!    end do
!    write(level_num,*) 'ssnow%wb '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%wb(i,1:6)
!    end do
!   
!    write(level_num,*) 'ssnow%wbice '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%wbice(i,1:6)
!    end do
!    write(level_num,*) 'ssnow%wbliq '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%wbliq(i,1:6)
!    end do
!
!    write(level_num,*) 'ssnow%GWwb '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%GWwb(i)
!    end do
!    write(level_num,*) 'ssnow%qhz '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%GWwb(i)
!    end do
!    write(level_num,*) 'ssnow%qhlev '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%qhlev(i,1:6)
!    end do
!    write(level_num,*) 'ssnow%qrechrage '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%qrecharge(i)
!    end do
!    write(level_num,*) 'ssnow%snowd '
!    do i=1,mp
!       write(level_num,*) i,' ',ssnow%snowd(i)
!    end do
!
!    end subroutine write_more_stuff
!
!    subroutine write_stuff()
!
!    real(r_2), dimension(ms) :: tmp_lev_vals
!
!    do k=1,ms
!       tmp_lev_vals(k) = maxval(ssnow%tgg(:,k),dim=1)
!    end do
!    write(level_num,*) 'ssnow%tgg max ',tmp_lev_vals(:)
!    do k=1,ms
!       tmp_lev_vals(k) = minval(ssnow%tgg(:,k),dim=1)
!    end do
!    write(level_num,*) 'ssnow%tgg min ',tmp_lev_vals(:)
!
!    do k=1,ms
!       tmp_lev_vals(k) = maxval(ssnow%wb(:,k),dim=1)
!    end do
!    write(level_num,*) 'max ssnow%wb ',tmp_lev_vals(:)
!    do k=1,ms
!       tmp_lev_vals(k) = minval(ssnow%wb(:,k),dim=1)
!    end do
!    write(level_num,*) 'min ssnow%wb ',tmp_lev_vals(:)
!   
!    write(level_num,*) 'max ssnow%wtd ',maxval(ssnow%wtd,dim=1)
!    write(level_num,*) 'max canopy%fes ',maxval(canopy%fes,dim=1)
!    write(level_num,*) 'max canopy%fess ',maxval(canopy%fess,dim=1)
!    write(level_num,*) 'max canopy%fev ',maxval(canopy%fev,dim=1)
!    write(level_num,*) 'max canopy%fns ',maxval(canopy%fns,dim=1)
!    write(level_num,*) 'max canopy%tscrn ',maxval(canopy%tscrn,dim=1)
!    write(level_num,*) 'max canopy%fe ',maxval(canopy%fe,dim=1)
!    write(level_num,*) 'max canopy%fh ',maxval(canopy%fh,dim=1)
!    write(level_num,*) 'max canopy%us ',maxval(canopy%us,dim=1)
!    write(level_num,*) 'max met%tvair ',maxval(met%tvair,dim=1)
!    write(level_num,*) 'max met%qvair ',maxval(met%qvair,dim=1)
!
!    write(level_num,*) 'min ssnow%wtd ',minval(ssnow%wtd,dim=1)
!    write(level_num,*) 'min canopy%fes ',minval(canopy%fes,dim=1)
!    write(level_num,*) 'min canopy%fess ',minval(canopy%fess,dim=1)
!    write(level_num,*) 'min canopy%fev ',minval(canopy%fev,dim=1)
!    write(level_num,*) 'min canopy%fns ',minval(canopy%fns,dim=1)
!    write(level_num,*) 'min canopy%tscrn ',minval(canopy%tscrn,dim=1)
!    write(level_num,*) 'min canopy%fe ',minval(canopy%fe,dim=1)
!    write(level_num,*) 'min canopy%fh ',minval(canopy%fh,dim=1)
!    write(level_num,*) 'min canopy%us ',minval(canopy%us,dim=1)
!    write(level_num,*) 'min met%tvair ',minval(met%tvair,dim=1)
!    write(level_num,*) 'min met%qvair ',minval(met%qvair,dim=1)
!
!    end subroutine write_stuff

END SUBROUTINE soil_snow_gw


SUBROUTINE calc_equilibrium_water_content(ssnow,soil)
  !find layer mean soil moisture and potential at equilibrium with wtd

  IMPLICIT NONE
  
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    !local variables
    REAL(r_2), dimension(mp)    :: zaq      !node depth of the aquifer
    REAL(r_2), dimension(ms)    :: dzmm     !layer thickness for single tile
    REAL(r_2), dimension(mp)    :: GWdzmm   !aquifer thickness at each tile
    REAL(r_2), dimension(mp)    :: GWzimm   !aquifer layer interface depth
    REAL(r_2), dimension(0:ms)  :: zimm     !layer interface depth in mm  
    REAL(r_2), dimension(ms)    :: zmm      !node depths in mm
    REAL(r_2)                   :: tempi, temp0,voleq1,wbrat

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
                          0.001_r_2),1._r_2)

          ssnow%zq(i,k) = max(sucmin,min(-soil%sucs_vec(i,k),&
                            -soil%sucs_vec(i,k)*(wbrat**(-soil%bch_vec(i,k))) ))
          
          
       end do  !mp
    end do  !ms
 
    do i=1,mp
    !Aquifer Equilibrium water content
       if (ssnow%wtd(i) .le. zimm(ms)) then      !fully saturated

          ssnow%GWwbeq(i) = soil%GWssat_vec(i)-soil%GWwatr(i)

       elseif ((ssnow%wtd(i) .gt. GWzimm(i)))   then     !fully unsaturated

          tempi = &
                (((soil%GWsucs_vec(i)+ssnow%wtd(i)-GWzimm(i))/&
                  soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
          temp0 = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(ms))/&
                     soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))   
          ssnow%GWwbeq(i) = -soil%GWsucs_vec(i)*soil%GWssat_vec(i)/&
                          (1._r_2-1._r_2/soil%GWbch_vec(i))/&
                           (GWzimm(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)   

       else           

          tempi  = 1._r_2
          temp0  = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(ms))/&
                     soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))               
          voleq1 = -soil%GWsucs_vec(i)*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
                    (1._r_2-1._r_2/soil%GWbch_vec(i))/&
                    (ssnow%wtd(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)
          ssnow%GWwbeq(i) = (voleq1*(ssnow%wtd(i)-zimm(ms)) + &
                            (soil%GWssat_vec(i)-soil%GWwatr(i))*&
                         (GWzimm(i)-ssnow%wtd(i)))/(GWzimm(i)-zimm(ms)) + soil%GWwatr(i)

       end if

       ssnow%GWwbeq(i) = min(max(ssnow%GWwbeq(i),soil%GWwatr(i)),soil%GWssat_vec(i))

       ssnow%GWzq(i) = -soil%GWsucs_vec(i)*(max((ssnow%GWwbeq(i)-soil%GWwatr(i))/     &
                    (soil%GWssat_vec(i)-soil%GWwatr(i)),0.001_r_2))**(-soil%GWbch_vec(i))
       ssnow%GWzq(i) = max(sucmin, ssnow%GWzq(i))
       
    end do

END SUBROUTINE calc_equilibrium_water_content

SUBROUTINE calc_srf_wet_fraction(ssnow,soil,met,veg)
  USE cable_data_module, only : math,phys

  IMPLICIT NONE
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)    :: veg
    TYPE (met_type), INTENT(IN)   :: met

    !local variables
    REAL(r_2), DIMENSION(mp)           :: icef,satfrac_liqice,S
    REAL(r_2)                          :: fice,xx
    REAL(r_2)                          :: dzmm_one,liqmass,icemass,totmass
    INTEGER                            :: i,j,k
    REAL(r_2)                          :: wb_unsat,wb_lin,funcval
    REAL(r_2)                          :: derv,slopeSTDmm,func_step
    REAL(r_2)                          :: wb_evap_threshold
    real, pointer :: tfrozen,ice_density,liq_density,pi

    tfrozen => PHYS%tfrz
    ice_density => PHYS%density_ice
    liq_density => PHYS%density_liq
    pi  => math%pi_c

    call saturated_fraction(ssnow,soil,veg)

    IF (cable_user%or_evap) THEN
       do i=1,mp
          IF( veg%iveg(i) == 16 .and. met%tk(i) < tfrozen + 5. )   THEN
             ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed
          ELSE
             ssnow%wetfac = 1.0
          ENDIF
       END DO

    ELSEIF (cable_user%gw_model) THEN

       dzmm_one  = real(soil%zse(1),r_2)
   
       do i = 1,mp
          icemass  = ssnow%wbice(i,1) * dzmm_one * ice_density
          liqmass  = (ssnow%wb(i,1)-ssnow%wbice(i,1)) * dzmm_one*liq_density
          totmass  = max(liqmass+icemass,real(1e-2,r_2))
          icef(i)     = max(0._r_2,min(1._r_2, gw_params%IceBeta*icemass / totmass))
      end do
   
      !srf frozen fraction.  should be based on topography
      do i = 1,mp
         fice = (exp(-gw_params%IceAlpha*(1._r_2-icef(i)))-&
                 exp(-gw_params%IceAlpha))/&
                 (1._r_2-exp(-gw_params%IceAlpha))
         fice = min(1._r_2,max(0._r_2,fice))
   
         satfrac_liqice(i) = fice + (1._r_2-fice)*ssnow%satfrac(i)
   
         wb_unsat = ((ssnow%wb(i,1)-ssnow%wbice(i,1)) -&
                     ssnow%satfrac(i)*soil%ssat_vec(i,1))/(1.-ssnow%satfrac(i))
         wb_unsat = min(soil%ssat_vec(i,1),max(0.,wb_unsat))
   
         wb_evap_threshold = min( max( &
                             gw_params%SoilEvapAlpha*soil%sfc_vec(i,1), &
                             soil%swilt_vec(i,1) ), soil%ssat_vec(i,1) )
   
         !Sakguchi and Zeng 2009
         if (wb_unsat .ge. wb_evap_threshold) then
            xx = 1.
         else
            xx = 0.25 * (1._r_2 - cos(real(pi,r_2)*wb_unsat/(wb_evap_threshold)))**2.0
         end if
   
         ssnow%wetfac(i) = max(0.0,min(1.0,satfrac_liqice(i) +&
                                        (1. - satfrac_liqice(i))*xx ) )
         IF ( veg%iveg(i) == 16 .and. met%tk(i) >= tfrozen + 5. )   &
              ssnow%wetfac(i) = 1.0 ! lakes: hard-wired number to be removed
   
         IF( veg%iveg(i) == 16 .and. met%tk(i) < tfrozen + 5. )   &
              ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed

      end do

   ELSE  !Default formulation

       ssnow%wetfac = MAX( 1.e-6, MIN( 1.0,&
            ( REAL (ssnow%wb(:,1) ) - soil%swilt/ 2.0 )                  &
            / ( soil%sfc - soil%swilt/2.0 ) ) )
   
       DO i=1,mp
   
          IF( ssnow%wbice(i,1) > 0. )&
               ssnow%wetfac(i) = ssnow%wetfac(i) * &
                                real(MAX( 0.5_r_2, 1._r_2 - MIN( 0.2_r_2, &
                                 ( ssnow%wbice(i,1) / ssnow%wb(i,1) )**2 ) ) )
   
          IF( ssnow%snowd(i) > 0.1) ssnow%wetfac(i) = 0.9
   
          IF ( veg%iveg(i) == 16 .and. met%tk(i) >= tfrozen + 5. )   &
               ssnow%wetfac(i) = 1.0 ! lakes: hard-wired number to be removed
   
          IF( veg%iveg(i) == 16 .and. met%tk(i) < tfrozen + 5. )   &
               ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed
   
       ENDDO
       ! owetfac introduced to reduce sharp changes in dry regions,
       ! especially in offline runs in which there may be discrepancies b/n
       ! timing of precip and temperature change (EAK apr2009)
       ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)

   ENDIF  !or_evap, gw_model, or default wetfac parameterization

END SUBROUTINE calc_srf_wet_fraction

SUBROUTINE calc_soil_hydraulic_props(ssnow,soil,veg)
   TYPE(soil_parameter_type), INTENT(IN)     :: soil 
   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
   TYPE(veg_parameter_type) , INTENT(IN)     :: veg

   INTEGER :: i,k,kk

   REAL(r_2), DIMENSION(mp) :: s1, &  !temporary variables for calculating hydraulic properties
                               s2, &
                               s_mid

   REAL(r_2), DIMENSION(0:ms) :: zimm  !depths at interface between layers
    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water (calculated using total (not liquid))

    do k=1,ms
       do i=1,mp
          ssnow%icefrac(i,k) = ssnow%wbice(i,k)/(max(ssnow%wb(i,k),0.01_r_2))
          ssnow%fracice(i,k) = (exp(-gw_params%IceAlpha*(1._r_2-ssnow%icefrac(i,k)))&
                               -exp(-gw_params%IceAlpha))/(1._r_2-exp(-gw_params%IceAlpha))
       end do
    end do

    ssnow%fracice(:,:) = max( min( ssnow%fracice, 1._r_2), 0._r_2)

    do k=1,ms-1
       kk=k+1
       do i=1,mp

          s1(i) = 0.5_r_2*(max(ssnow%wb(i,k)-soil%watr(i,k),0.) + &
                           max(ssnow%wb(i,kk)-soil%watr(i,kk),0.)) / &
                         (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + &
                         (soil%ssat_vec(i,kk)-soil%watr(i,kk))))

          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)

          ssnow%hk(i,k)    =  (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))*s1(i)*s2(i)
          ssnow%dhkdw(i,k) = (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))* &
                             (2._r_2*soil%bch_vec(i,k)+3._r_2)*s2(i)*0.5_r_2/&
                              (soil%ssat_vec(i,k)-soil%watr(i,k))

       end do
    end do

    k = ms 
       do i=1,mp

          s1(i) = 0.5_r_2*(max(ssnow%wb(i,k)-soil%watr(i,k),0.) + &
                           max(ssnow%GWwb(i)-soil%GWwatr(i),0.)) / &
                  (0.5_r_2*(soil%ssat_vec(i,k)-soil%watr(i,k) +&
                            soil%GWssat_vec(i)-soil%GWwatr(i)))

          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)
          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)

          ssnow%hk(i,k)    = (1._r_2 - ssnow%fracice(i,k))*s1(i)*s2(i)
          ssnow%dhkdw(i,k) = (1._r_2 - ssnow%fracice(i,k))*&
                             (2._r_2*soil%bch_vec(i,k)+3._r_2)*&
                             s2(i)*0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))
       end do
 
    do k=1,ms 
       do i=1,mp
          s_mid(i) = (ssnow%wb(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
              (soil%ssat_vec(i,k)-soil%watr(i,k))

          s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)

          ssnow%smp(i,k) = -soil%sucs_vec(i,k)*s_mid(i)**(-soil%bch_vec(i,k))

          ssnow%smp(i,k) = max(min(ssnow%smp(i,k),-soil%sucs_vec(i,k)),sucmin)

          ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/&
                    (max(s_mid(i)*(soil%ssat_vec(i,k)-soil%watr(i,k)),0.001_r_2))       
       end do   
    end do

    do i=1,mp
       !Aquifer properties
       s_mid(i) = (ssnow%GWwb(i)-soil%GWwatr(i))/&
                    (soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)
       s2(i)    =&
            soil%GWhyds_vec(i)*s_mid(i)**(2._r_2*soil%GWbch_vec(i)+2._r_2)

       ssnow%GWhk(i)     =s_mid(i)*s2(i)

       ssnow%GWdhkdw(i)  =  (2._r_2*soil%GWbch_vec(i)+3._r_2)*&
                           s2(i)*0.5_r_2/(soil%GWssat_vec(i)-soil%GWwatr(i))

          ssnow%GWhk(i)    = (1.-ssnow%fracice(i,ms)) * ssnow%GWhk(i)
          ssnow%GWdhkdw(i) = (1.-ssnow%fracice(i,ms)) * ssnow%GWdhkdw(i)

       s_mid(i) = (ssnow%GWwb(i)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)

       ssnow%GWsmp(i)    = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))
       ssnow%GWsmp(i)    = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)
       ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/&
                            (s_mid(i)*(soil%GWssat_vec(i)-soil%GWwatr(i)))
    end do

END SUBROUTINE calc_soil_hydraulic_props

  SUBROUTINE aquifer_recharge(dt,ssnow,soil,veg,zaq,zmm,dzmm)

  IMPLICIT NONE
    real, intent(in) :: dt 
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)     :: veg
    REAL(r_2), dimension(:), intent(in)       :: zaq
    REAL(r_2), dimension(:), intent(in)       :: zmm,dzmm

    integer :: i    

    !Doing the recharge outside of the soln of Richards Equation makes it easier to track total recharge amount.
    !Add to ssnow at some point 
    do i=1,mp
       if ((ssnow%wtd(i) .le. sum(dzmm,dim=1)) .or. &
           (veg%iveg(i) .ge. 16) .or. &
           (soil%isoilm(i) .eq. 9))  then

          ssnow%Qrecharge(i) = 0._r_2
       else
          ssnow%Qrecharge(i) = -ssnow%hk(i,ms)*&
                               ((ssnow%GWsmp(i)-ssnow%smp(i,ms)) -&
                                (ssnow%GWzq(i)-ssnow%zq(i,ms))) / &
                                (zaq(i) - zmm(ms))
       end if
    end do


  END SUBROUTINE aquifer_recharge

  SUBROUTINE subsurface_drainage(ssnow,soil,veg,dzmm)

  IMPLICIT NONE
  
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)     :: veg
    REAL(r_2), dimension(:), intent(in)       :: dzmm
    REAL(r_2), dimension(mp)                  :: sm_tot
    INTEGER, dimension(mp)                    :: k_drain
    integer :: i,k

    real(r_2), dimension(17) :: Efold_mod

    Efold_mod(:) = 1.0
    !Efold_mod(1:4) = (/0.2,0.2,0.2,0.2/)
    !Efold_mod(9) = 0.25
    do i=1,mp

       !Note: future revision will have interaction with river here. nned to
       !work on router and add river type cells
       ssnow%qhz(i)  = min(max(soil%slope(i),0.00001),0.1)*&
                       gw_params%MaxHorzDrainRate* &
                        exp(-ssnow%wtd(i)/(1000._r_2*&
                       (gw_params%EfoldHorzDrainRate*Efold_mod(veg%iveg(i)))))

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

       if (gw_params%subsurface_sat_drainage) then
          sm_tot(i) = max((ssnow%GWwb(i) - soil%watr(i,ms))*&
                          (1._r_2-ssnow%fracice(i,ms)), 0.)

          do k=k_drain(i),ms
             sm_tot(i) = sm_tot(i) +max(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)
          end do

          if (sm_tot(i) .ge. 1.0e-12) then
              do k=k_drain(i),ms
                 ssnow%qhlev(i,k) = ssnow%qhz(i)*max(&
                                   ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)/sm_tot(i)
              end do
              ssnow%qhlev(i,ms+1) = max((ssnow%GWwb(i) - soil%watr(i,ms))*&
                                    (1._r_2-ssnow%fracice(i,ms)), 0.)*ssnow%qhz(i)/sm_tot(i)
          endif

       else  !second option
          if (k_drain(i) .le. ms) then
             sm_tot(i) = 0.
             do k=k_drain(i),ms
                sm_tot(i) = sm_tot(i) + max(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)
             end do

             if (sm_tot(i) .ge. 1.0e-12) then
                 do k=k_drain(i),ms
                    ssnow%qhlev(i,k) = ssnow%qhz(i)*max(&
                                        ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)/sm_tot(i)
                 end do
             endif

          else

             ssnow%qhlev(i,ms+1) = ssnow%qhz(i)*(1._r_2-ssnow%fracice(i,ms))

          end if

       end if
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
    TYPE (veg_parameter_type), INTENT(IN)     :: veg


    REAL(r_2), DIMENSION(mp) :: S
    REAL(r_2) :: slopeSTDmm
    INTEGER :: i,k

     S(:) = 0._r_2
     do k=1,ms
       S(:) = S(:) + max(0.001,min(1.0, &
              (ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/&
               max(0.001,soil%ssat_vec(:,k)-soil%watr(:,k)) ) )*soil%zse(k)
     end do
     S(:) = S(:)/sum(soil%zse(1:ms),dim=1)
     !srf frozen fraction.  should be based on topography
      do i = 1,mp
         !Saturated fraction
          if (gw_params%MaxSatFraction .gt. 1e-7) then 
             if (veg%iveg(i) .lt. 16 ) then
               slopeSTDmm = sqrt(min(max(&
                             gw_params%MaxSatFraction*soil%slope_std(i),&
                             1e-5),10000._r_2)) ! ensure some variability
               ssnow%satfrac(i)    = max(0._r_2,min(0.95_r_2,&
                   !note UM wants std03, and erf is not included then
                                     1._r_2 - my_erf( slopeSTDmm / sqrt(2.0* S(i)) ) ) )  
             else
                ssnow%satfrac(i) = 0.5
             end if

             ssnow%satfrac(i) = min(0.995,max(0.005,ssnow%satfrac(i) ) )
          else
             ssnow%satfrac(i) = 0._r_2
          end if
      end do


  END SUBROUTINE saturated_fraction

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

   !this function is only here because the um by default compiles with
   !-std03 and the 2003 standard does not include or allow
   !the erf function
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


END MODULE cable_gw_hydro_module

