!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Calls CABLE routines including define_air, surface_albedo,
!          define_canopy, soilsnow, carbon
!          Note that cbm is called once per timestep in the offline case but
!          twice per timestep in the ACCESS case. Not all parts of cbm
!          are executed in each of the ACCESS calls.
!
! Called from: cable_driver for offline version
!              cable_explicit_driver, cable_implicit_driver for ACCESS
!
! Contact: Yingping.Wang@csiro.au
!
! History: Calling sequence changes for ACCESS compared to v1.4b
!
!          REV_CORR package of fixes for the sensitivity/correction terms 
!
! ==============================================================================

!#define NO_CASA_YET 1

MODULE cable_cbm_module


   IMPLICIT NONE

   PRIVATE
   PUBLIC cbm

CONTAINS

   SUBROUTINE cbm( ktau,dels, air, bgc, canopy, met,                                &
                   bal, rad, rough, soil,                                      &
                   ssnow, sum_flux, veg, climate )

   USE cable_common_module
   USE cable_carbon_module

   USE cable_def_types_mod
   USE cable_roughness_module, only : ruff_resist
   USE cable_air_module, only : define_air
#ifndef NO_CASA_YET
   USE casadimension,     only : icycle ! used in casa_cnp
#endif
   !mrd561
   USE cable_gw_hydro_module, only : sli_hydrology,&
                                     soil_snow_gw
   USE cable_canopy_module, only : define_canopy
   USE sli_main_mod, only : sli_main
                                   
!subrs:
USE cbl_albedo_mod,             ONLY: albedo
USE cbl_init_radiation_module,  ONLY: init_radiation
USE cbl_soil_snow_main_module, ONLY : soil_snow
USE cbl_masks_mod, ONLY: fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
USE snow_aging_mod, ONLY : snow_aging 

!jhan:pass these !data
USE cable_other_constants_mod, ONLY: Ccoszen_tols => coszen_tols
USE cable_other_constants_mod,  ONLY : Crad_thresh => rad_thresh
USE cable_other_constants_mod, ONLY: clai_thresh => lai_thresh
USE cable_other_constants_mod, ONLY: cgauss_w => gauss_w
USE cable_math_constants_mod,  ONLY: cpi => pi
USE cable_math_constants_mod,  ONLY: cpi180 => pi180
USE cable_phys_constants_mod,  ONLY: cEMLEAF=> EMLEAF
USE cable_phys_constants_mod,  ONLY: cEMSOIL=> EMSOIL
USE cable_phys_constants_mod,  ONLY: cSBOLTZ=> SBOLTZ
USE grid_constants_mod_cbl, ONLY : ICE_SoilType, lakes_cable

   !ptrs to local constants
   ! CABLE model variables
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (bgc_pool_type),  INTENT(INOUT) :: bgc
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (balances_type),  INTENT(INOUT) :: bal
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (sum_flux_type),  INTENT(INOUT) :: sum_flux
   TYPE (climate_type), INTENT(IN)      :: climate

   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg

   REAL, INTENT(IN)               :: dels ! time setp size (s)
   INTEGER, INTENT(IN) :: ktau
   INTEGER :: k,kk,j
   logical, save :: first_call = .true.
#ifdef NO_CASA_YET
   INTEGER :: ICYCLE
   ICYCLE = 0
#endif
LOGICAL :: veg_mask(mp), sunlit_mask(mp), sunlit_veg_mask(mp)

character(len=*), parameter :: subr_name = "cbm"
LOGICAL :: cbl_standalone= .true.
LOGICAL :: jls_standalone= .false.
LOGICAL :: jls_radiation= .false.

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)


cable_user%soil_struc="default"

cable_runtime%um_radiation = .FALSE.

IF( cable_runtime%um_explicit ) THEN
  CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )
ENDIF
   ! Height adjustment not used in ACCESS CM2. See CABLE ticket 197
   ! met%tk = met%tk + C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
CALL define_air (met, air)

call fveg_mask( veg_mask, mp, Clai_thresh, canopy%vlaiw )
call fsunlit_mask( sunlit_mask, mp, CRAD_THRESH,( met%fsd(:,1)+met%fsd(:,2) ) )
call fsunlit_veg_mask( sunlit_veg_mask, veg_mask, sunlit_mask, mp )

CALL init_radiation( rad%extkb, rad%extkd,                                     &
                     !ExtCoeff_beam, ExtCoeff_dif,
                     rad%extkbm, rad%extkdm, Rad%Fbeam,                        &
                     !EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,
                     c1, rhoch, xk,                                            &
                     mp,nrb,                                                   &
                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,         &
                     cbl_standalone, jls_standalone, jls_radiation,            &
                     subr_name,                                                &
                     veg_mask,                                                 &
                     veg%Xfang, veg%taul, veg%refl,                            &
                     !VegXfang, VegTaul, VegRefl
                     met%coszen, int(met%DoY), met%fsd,                        &
                     !coszen, metDoY, SW_down,
                     canopy%vlaiw                                              &
                   ) !reducedLAIdue2snow 

!Ticket 331 refactored albedo code for JAC
CALL snow_aging(ssnow%snage,mp,dels,ssnow%snowd,ssnow%osnowd,ssnow%tggsn(:,1),&
         ssnow%tgg(:,1),ssnow%isflag,veg%iveg,soil%isoilm) 

IF( cable_runtime%um_explicit ) THEN

  call Albedo( ssnow%AlbSoilsn, soil%AlbSoil,                                &
               !AlbSnow, AlbSoil,              
               mp, nrb,                                                      &
               ICE_SoilType, lakes_cable,                                    &
               jls_radiation,                                                &
               veg_mask,                                                     &  
               Ccoszen_tols, CGAUSS_W,                                       & 
               veg%iveg, soil%isoilm, veg%refl, veg%taul,                    & 
               !surface_type, VegRefl, VegTaul,
               met%coszen, canopy%vlaiw,                                     &
               !coszen, reducedLAIdue2snow,
               ssnow%snowd, ssnow%ssdnn, ssnow%tgg(:,1), ssnow%snage,        & 
               !SnowDepth, SnowDensity, SoilTemp, SnowAge, 
               xk, c1, rhoch,                                                & 
               rad%fbeam, rad%albedo,                                        &
               !RadFbeam, RadAlbedo,
               rad%extkb, rad%extkd,                                         & 
               !ExtCoeff_beam, ExtCoeff_dif,
               rad%extkbm, rad%extkdm,                                       & 
               !EffExtCoeff_beam, EffExtCoeff_dif,                
               rad%rhocbm, rad%rhocdf,                                       &
               !CanopyRefl_beam,CanopyRefl_dif,
               rad%cexpkbm, rad%cexpkdm,                                     & 
               !CanopyTransmit_beam, CanopyTransmit_dif, 
               rad%reffbm, rad%reffdf                                        &
             ) !EffSurfRefl_beam, EffSurfRefldif_
  
ENDIF

! Calculate canopy variables:

!! vh_js !!
!CABLE_LSM:check
IF( first_call ) then
  ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1) 
  ssnow%otss = ssnow%tss
  first_call = .false.
endif
ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
ssnow%otss = ssnow%tss

CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, sunlit_veg_mask,  canopy%vlaiw)
! RML moved out of following IF after discussion with Eva
ssnow%owetfac = ssnow%wetfac

IF( cable_runtime%um_implicit ) THEN
   IF (cable_user%gw_model) then
      CALL soil_snow_gw(dels, soil, ssnow, canopy, met, bal,veg)
   ELSE
       CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
    ENDIF
 ENDIF

ssnow%deltss = ssnow%tss-ssnow%otss

! correction required for energy balance in online simulations
! REV_CORR - multiple changes to address %cls bugs and revised correction
! terms.  Also - do not apply correction terms if using SLI
! SSEB package will move these calculations to within soilsnow
IF( cable_user%SOIL_STRUC=='default') THEN

   canopy%fhs = canopy%fhs + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
   canopy%fhs_cor = canopy%fhs_cor + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
   canopy%fh = canopy%fhv + canopy%fhs

   !canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) *                    &
   !          ( ssnow%dfe_ddq * ssnow%ddq_dtg )
   !          !( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )
   !
   !Ticket 137 - remove double couting of %cls
   !canopy%fes_cor = canopy%fes_cor + ( ssnow%tss-ssnow%otss ) *            &
   !                      ( ssnow%dfe_ddq * ssnow%ddq_dtg )
   !               ( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )

   !INH rewritten in terms of %dfe_dtg - NB factor %cls above was a bug
   canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) * ssnow%dfe_dtg

   !INH NB factor %cls in %fes_cor above was a bug - see Ticket #135 #137
   canopy%fes_cor = canopy%fes_cor + (ssnow%tss-ssnow%otss) * ssnow%dfe_dtg
   !canopy%fes_cor = canopy%fes_cor + ssnow%cls*(ssnow%tss-ssnow%otss) & 
   !       * ssnow%dfe_dtg
   
   IF (cable_user%L_REV_CORR) THEN
      !INH need to add on corrections to all terms in the soil energy balance
      canopy%fns_cor = canopy%fns_cor + (ssnow%tss-ssnow%otss)*ssnow%dfn_dtg

      !NB %fns_cor also added onto out%Rnet and out%LWnet in cable_output and
      !cable_checks as the correction term needs to pass through the 
      !canopy in entirity not be partially absorbed and %fns not used there 
      !(as would be the case if rad%flws were changed)
      canopy%fns = canopy%fns + ( ssnow%tss-ssnow%otss )*ssnow%dfn_dtg

      canopy%ga_cor = canopy%ga_cor + ( ssnow%tss-ssnow%otss )*canopy%dgdtg
      canopy%ga = canopy%ga + ( ssnow%tss-ssnow%otss )*canopy%dgdtg

      !assign all the correction to %fes to %fess - none to %fesp
      canopy%fess = canopy%fess + ( ssnow%tss-ssnow%otss ) * ssnow%dfe_dtg

   ENDIF
ENDIF
   

   ! need to adjust fe after soilsnow
   canopy%fev  = canopy%fevc + canopy%fevw

   ! Calculate total latent heat flux:
   canopy%fe = canopy%fev + canopy%fes

   ! Calculate net radiation absorbed by soil + veg
   canopy%rnet = canopy%fns + canopy%fnv

   ! Calculate radiative/skin temperature:
   if (cable_runtime%um) then
       !Jan 2018: UM assumes a single emissivity for the surface in the radiation scheme
       !To accommodate this a single value of is 1. is assumed in ACCESS
       ! any leaf/soil emissivity /=1 must be incorporated into rad%trad.  
       ! check that emissivities (pft and nvg) set = 1 within the UM i/o configuration
       ! CM2 - further adapted to pass the correction term onto %trad correctly
       rad%trad = ( ( 1.-rad%transd ) * Cemleaf * canopy%tv**4 +                      &
              rad%transd * Cemsoil * ssnow%otss**4 + canopy%fns_cor/CSBOLTZ )**0.25
   else       
   rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
              rad%transd * ssnow%tss**4 )**0.25
   endif

   ! rml 17/1/11 move all plant resp and soil resp calculations here
   ! from canopy. in UM only call on implicit step.
   ! put old and new soil resp calculations into soilcarb subroutine
   ! make new plantcarb subroutine
   IF (.not.cable_runtime%um_explicit .AND. icycle == 0) THEN

      !calculate canopy%frp
      CALL plantcarb(veg,bgc,met,canopy)

      !calculate canopy%frs
      CALL soilcarb(soil, ssnow, veg, bgc, met, canopy)

      CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

      canopy%fnpp = -1.0* canopy%fpn - canopy%frp
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

   ENDIF


END SUBROUTINE cbm

END MODULE cable_cbm_module


