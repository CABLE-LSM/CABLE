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

   USE cable_canopy_module
   USE cable_albedo_module

   IMPLICIT NONE

   PRIVATE
   PUBLIC cbm

CONTAINS

   SUBROUTINE cbm( ktau,dels, air, bgc, canopy, met,                                &
                   bal, rad, rough, soil,                                      &
                   ssnow, sum_flux, veg, climate )

   USE cable_common_module
   USE cable_carbon_module
   USE cable_soil_snow_module, only : soil_snow
   USE cable_def_types_mod
   USE cable_roughness_module
   USE cable_radiation_module
   USE cable_air_module
#ifndef NO_CASA_YET
   USE casadimension,     only : icycle ! used in casa_cnp
#endif
   USE cable_data_module, ONLY : icbm_type, point2constants
   !mrd561
   USE cable_gw_hydro_module

   !ptrs to local constants
   TYPE( icbm_type ) :: C
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
   TYPE (climate_type), INTENT(IN)    :: climate

   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg

   REAL, INTENT(IN)               :: dels ! time setp size (s)
   INTEGER, INTENT(IN) :: ktau
   INTEGER :: k,kk,j

#ifdef NO_CASA_YET
   INTEGER :: ICYCLE
   ICYCLE = 0
#endif

   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)

   IF( cable_runtime%um ) THEN

      cable_runtime%um_radiation = .FALSE.

      IF( cable_runtime%um_explicit ) THEN
         CALL ruff_resist(veg, rough, ssnow, canopy)
         met%tk = met%tk + C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
      ENDIF

      CALL define_air (met, air)

   ELSE
      call ruff_resist(veg, rough, ssnow, canopy)
   ENDIF
   
   CALL init_radiation(met,rad,veg, canopy) ! need to be called at every dt

   IF( cable_runtime%um ) THEN

      IF( cable_runtime%um_explicit ) THEN
         CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
      ENDIF

   ELSE
      CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
   ENDIf

   !! vh_js !!
   ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
   ssnow%otss = ssnow%tss

   ! Calculate canopy variables:
   IF (.not.(cable_user%or_evap .or. cable_user%gw_model .or. &
            (cable_user%SOIL_STRUC=='sli' .and. cable_user%test_new_gw) ) ) THEN
      call set_unsed_gw_vars(ssnow,soil,canopy)
   END IF

   CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate)
   
   !ssnow%otss_0 = ssnow%otss
   !ssnow%otss = ssnow%tss

   ! RML moved out of following IF after discussion with Eva
   ssnow%owetfac = ssnow%wetfac

   IF( cable_runtime%um ) THEN

     IF( cable_runtime%um_implicit ) THEN
        IF (cable_user%gw_model) then
           if (mp .ge. 1) CALL soil_snow_gw(dels, soil, ssnow, canopy, met, bal,veg)
        ELSE
            CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
         ENDIF
      ENDIF

   ELSE
      IF(cable_user%SOIL_STRUC=='default') THEN
        IF (cable_user%gw_model) then
           if (mp .ge. 1) CALL soil_snow_gw(dels, soil, ssnow, canopy, met, bal,veg)
        ELSE
            CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
         ENDIF
      ELSEIF (cable_user%SOIL_STRUC=='sli') THEN

         IF (cable_user%test_new_gw) &
            CALL sli_hydrology(dels,ssnow,soil,veg,canopy)

         CALL sli_main(ktau,dels,veg,soil,ssnow,met,canopy,air,rad,0)
      ENDIF
   ENDIF
   

   ssnow%deltss = ssnow%tss-ssnow%otss
   ! correction required for energy balance in online simulations
   ! REV_CORR - multiple changes to address %cls bugs and revised correction
   ! terms.  Also - do not apply correction terms if using SLI
   ! SSEB package will move these calculations to within soilsnow
   IF( cable_runtime%um .AND. cable_user%SOIL_STRUC=='default') THEN

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
   rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
              rad%transd * ssnow%tss**4 )**0.25

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


