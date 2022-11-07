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
!
! ==============================================================================

!#define NO_CASA_YET 1

MODULE cable_cbm_module

   USE cable_canopy_module
   USE cable_albedo_module

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: cbm

CONTAINS

   SUBROUTINE cbm(ktau, dels, air, bgc, canopy, met, &
                  bal, rad, rough, soil, &
                  ssnow, veg, climate )

   USE cable_common_module
   USE cable_carbon_module
   USE cable_soil_snow_module
   USE cable_def_types_mod
   USE cable_roughness_module
   USE cable_radiation_module
   USE cable_air_module
#ifndef NO_CASA_YET
   USE casadimension,     only: icycle ! used in casa_cnp
#endif
   USE cable_data_module, ONLY: icbm_type, point2constants
   use cable_sli_main,    only: sli_main

   ! CABLE model variables
   INTEGER,                   INTENT(IN)    :: ktau
   REAL,                      INTENT(IN)    :: dels ! time setp size (s)
   TYPE(air_type),            INTENT(INOUT) :: air
   TYPE(bgc_pool_type),       INTENT(INOUT) :: bgc
   TYPE(canopy_type),         INTENT(INOUT) :: canopy
   TYPE(met_type),            INTENT(INOUT) :: met
   TYPE(balances_type),       INTENT(INOUT) :: bal
   TYPE(radiation_type),      INTENT(INOUT) :: rad
   TYPE(roughness_type),      INTENT(INOUT) :: rough
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow
   TYPE(veg_parameter_type),  INTENT(INOUT) :: veg
   TYPE(climate_type),        INTENT(IN)    :: climate

   ! ptrs to local constants
   TYPE(icbm_type) :: C

#ifdef NO_CASA_YET
   INTEGER :: ICYCLE
   ICYCLE = 0
#endif

   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)

   IF (cable_runtime%um) THEN
      cable_runtime%um_radiation = .FALSE.

      IF (cable_runtime%um_explicit) THEN
         CALL ruff_resist(veg, rough, ssnow, canopy)
         met%tk = met%tk + C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
      ENDIF

      CALL define_air(met, air)
   ELSE
      call ruff_resist(veg, rough, ssnow, canopy)
   ENDIF

   CALL init_radiation(met, rad, veg, canopy) ! need to be called at every dt

   IF (cable_runtime%um) THEN
      IF (cable_runtime%um_explicit) THEN
         CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
      ENDIF
   ELSE
      CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
   ENDIf

   !! vh_js !!

   ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
   ssnow%otss   = ssnow%tss
   ssnow%owetfac = ssnow%wetfac ! MC should also be before canopy

   ! Calculate canopy variables
   CALL define_canopy(bal, rad, rough, air, met, dels, ssnow, soil, veg, canopy, climate)

   ! write(*,*) 'hod, TVeg: ', met%hod(1), canopy%fevc(1), canopy%fwsoil(1)
   ! if (met%hod(1).gt.12.0) stop

   !ssnow%otss_0 = ssnow%otss
   !ssnow%otss = ssnow%tss

   ! RML moved out of following IF after discussion with Eva
   ! ssnow%owetfac = ssnow%wetfac

   IF (cable_runtime%um) THEN
      IF (cable_runtime%um_implicit) THEN
         CALL soil_snow(dels, soil, ssnow, canopy, met, veg)
      ENDIF
   ELSE
      IF (cable_user%SOIL_STRUC=='default') THEN
         call soil_snow(dels, soil, ssnow, canopy, met, veg)
      ELSEIF (cable_user%SOIL_STRUC=='sli') THEN
         ! print*, 'SLIMAIN01 ', ktau, dels
         ! call print_cbm_var(veg)
         ! call print_cbm_var(soil)
         ! call print_cbm_var(ssnow)
         ! call print_cbm_var(met)
         ! call print_cbm_var(canopy)
         ! call print_cbm_var(air)
         ! call print_cbm_var(rad)
         CALL sli_main(ktau, dels, veg, soil, ssnow, met, canopy, air, rad, 0)
         ! print*, 'SLIMAIN02 ', ktau, dels
         ! call print_cbm_var(veg)
         ! call print_cbm_var(soil)
         ! call print_cbm_var(ssnow)
         ! call print_cbm_var(met)
         ! call print_cbm_var(canopy)
         ! call print_cbm_var(air)
         ! call print_cbm_var(rad)
      ENDIF
   ENDIF

   ssnow%deltss = ssnow%tss-ssnow%otss
   ! correction required for energy balance in online simulations
   IF (cable_runtime%um) THEN

      canopy%fhs = canopy%fhs + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg

      canopy%fhs_cor = canopy%fhs_cor + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg

      canopy%fh = canopy%fhv + canopy%fhs

      canopy%fes = canopy%fes + real(ssnow%tss-ssnow%otss, r_2) * &
                real(ssnow%dfe_ddq * ssnow%ddq_dtg, r_2)
                !( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )
      canopy%fes_cor = canopy%fes_cor + real(ssnow%tss-ssnow%otss, r_2) * &
                    real(ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg, r_2)

   ENDIF

   ! need to adjust fe after soilsnow
   canopy%fev  = real(canopy%fevc) + canopy%fevw

   ! Calculate total latent heat flux:
   canopy%fe = canopy%fev + real(canopy%fes)

   ! Calculate net radiation absorbed by soil + veg
   canopy%rnet = canopy%fns + canopy%fnv

   ! Calculate radiative/skin temperature:
   rad%trad = ( (1. - rad%transd) * canopy%tv**4 + &
        rad%transd * ssnow%tss**4 )**0.25

   ! rml 17/1/11 move all plant resp and soil resp calculations here
   ! from canopy. in UM only call on implicit step.
   ! put old and new soil resp calculations into soilcarb subroutine
   ! make new plantcarb subroutine
   IF ((.not. cable_runtime%um_explicit) .AND. (icycle == 0)) THEN

      !calculate canopy%frp
      CALL plantcarb(veg, bgc, met, canopy)

      !calculate canopy%frs
      CALL soilcarb(soil, ssnow, veg, bgc, canopy)

      CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

      canopy%fnpp = -1.0 * canopy%fpn - canopy%frp
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

   ENDIF

END SUBROUTINE cbm

END MODULE cable_cbm_module
