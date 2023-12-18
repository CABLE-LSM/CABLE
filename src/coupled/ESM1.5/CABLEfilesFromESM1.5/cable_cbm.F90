MODULE cable_cbm_module
   
   USE cable_canopy_module
  
   IMPLICIT NONE
  
   PRIVATE
   PUBLIC cbm 

CONTAINS

   SUBROUTINE cbm( dels, air, bgc, canopy, met,                                &
                   bal, rad, rough, soil,                                      &
                   ssnow, sum_flux, veg, climate )
    
USE cbl_init_radiation_module, ONLY: init_radiation
USE cbl_albedo_mod, ONLY: albedo
USE cbl_masks_mod, ONLY: fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
USE cbl_soil_snow_main_module, ONLY: soil_snow
USE snow_aging_mod, ONLY : snow_aging 

!jhan:pass these !data
USE cable_other_constants_mod, ONLY: Ccoszen_tols => coszen_tols
USE cable_other_constants_mod,  ONLY : Crad_thresh => rad_thresh
USE cable_other_constants_mod, ONLY: clai_thresh => lai_thresh
USE cable_other_constants_mod, ONLY: cgauss_w => gauss_w
USE cable_math_constants_mod,  ONLY: cpi => pi
USE cable_math_constants_mod,  ONLY: cpi180 => pi180
USE grid_constants_mod_cbl, ONLY : ICE_SoilType, lakes_cable


   USE cable_common_module
   USE cable_carbon_module
   USE cable_def_types_mod
   USE cable_roughness_module
   USE cable_air_module
#ifndef NO_CASA_YET
   USE casadimension,     only : icycle ! used in casa_cnp
#endif
   USE cable_data_module, ONLY : icbm_type, point2constants 
   
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
   TYPE(climate_type),    INTENT(INOUT) :: climate
    
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg  

   REAL, INTENT(IN)               :: dels ! time setp size (s)
    
!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)

LOGICAL :: veg_mask(mp), sunlit_mask(mp), sunlit_veg_mask(mp)
CHARACTER(LEN=*), PARAMETER :: subr_name = "cbl_model_driver"
LOGICAL :: jls_standalone= .FALSE.
LOGICAL :: jls_radiation= .FALSE.
LOGICAL :: cbl_standalone = .FALSE.    


   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)    

      cable_runtime%um_radiation = .FALSE.
      
      IF( cable_runtime%um_explicit ) THEN
        CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )
         met%tk = met%tk + C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
      ENDIF
      
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
 
IF( cable_runtime%um_explicit ) THEN
   
 !Ticket 331 refactored albedo code for JAC
 CALL snow_aging(ssnow%snage,mp,dels,ssnow%snowd,ssnow%osnowd,ssnow%tggsn(:,1),&
         ssnow%tgg(:,1),ssnow%isflag,veg%iveg,soil%isoilm) 
         
call Albedo( ssnow%AlbSoilsn, soil%AlbSoil,                                &
             mp, nrb, ICE_SoilType, lakes_cable, jls_radiation, veg_mask,       &
             Ccoszen_tols, cgauss_w,                                        & 
             veg%iveg, soil%isoilm, veg%refl, veg%taul,                     & 
             !surface_type, VegRefl, VegTaul,
             met%coszen, canopy%vlaiw,                                      &
             !coszen, reducedLAIdue2snow,
             ssnow%snowd, ssnow%ssdnn, ssnow%tgg(:,1), ssnow%snage,         &
             !SnowDepth, SnowDensity, SoilTemp, SnowAge,  
             xk, c1, rhoch,                                                 & 
             rad%fbeam, rad%albedo,                                         &
             !RadFbeam, RadAlbedo,
             rad%extkb, rad%extkd,                                          & 
             !ExtCoef_beamf, ExtCoeff,
             rad%extkbm, rad%extkdm,                                        & 
             !EffExtCoeff_beam, EffExtCoeff_dif,                
             rad%rhocbm, rad%rhocdf,                                        &
             !CanopyRefl_beam,CanopyRefl_dif,
             rad%cexpkbm, rad%cexpkdm,                                      & 
             !CanopyTransmit_beam, CanopyTransmit_dif, 
             rad%reffbm, rad%reffdf                                        &
           ) !EffSurfRefl_beam, EffSurfRefl_dif


         
      ENDIF
   
   ! Calculate canopy variables:
   CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy, climate, sunlit_veg_mask, canopy%vlaiw )

   ssnow%otss_0 = ssnow%otss
   ssnow%otss = ssnow%tss
   ! RML moved out of following IF after discussion with Eva
   ssnow%owetfac = ssnow%wetfac

   IF( cable_runtime%um ) THEN
      
     IF( cable_runtime%um_implicit ) THEN
         CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
      ENDIF

   ELSE
      call soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
   ENDIF

   ssnow%deltss = ssnow%tss-ssnow%otss
   ! correction required for energy balance in online simulations
   IF( cable_runtime%um ) THEN
   
      canopy%fhs = canopy%fhs + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
      
      canopy%fhs_cor = canopy%fhs_cor + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
      
      canopy%fh = canopy%fhv + canopy%fhs

   canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) *                        &
                ( ssnow%dfe_ddq * ssnow%ddq_dtg )
                !( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )
   
   canopy%fes_cor = canopy%fes_cor + ( ssnow%tss-ssnow%otss ) *                &
                    ( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )

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

END SUBROUTINE cbm

END MODULE cable_cbm_module


