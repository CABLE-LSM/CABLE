MODULE cable_cbm_module
   
   USE cable_canopy_module
  
   IMPLICIT NONE
  
   PRIVATE
   PUBLIC cbm 

CONTAINS

   SUBROUTINE cbm( dels, air, bgc, canopy, met,                                &
                   bal, rad, rough, soil,                                      &
                   ssnow, sum_flux, veg )
    
   USE cable_common_module
   USE cable_carbon_module
USE cbl_soil_snow_main_module,  ONLY: soil_snow
!restrict with ONLY syntax
   USE cable_def_types_mod
   USE cable_roughness_module
   USE cbl_init_radiation_module, ONLY: init_radiation
   USE cable_air_module
!CBL3 
USE cbl_albedo_mod, ONLY: albedo
USE cbl_masks_mod, ONLY: fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
USE cbl_masks_mod, ONLY: veg_mask,  sunlit_mask,  sunlit_veg_mask
!jhan:pass these !data
USE cable_other_constants_mod, ONLY: Ccoszen_tols => coszen_tols
USE cable_other_constants_mod,  ONLY : Crad_thresh => rad_thresh
USE cable_other_constants_mod, ONLY: clai_thresh => lai_thresh
USE cable_other_constants_mod, ONLY: cgauss_w => gauss_w
USE cable_other_constants_mod, ONLY : cmax_kLAI => max_kLAI
USE cable_math_constants_mod,  ONLY: cpi => pi
USE cable_math_constants_mod,  ONLY: cpi180 => pi180
USE cable_climate_type_mod, ONLY : climate_cbl

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
    
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg  

   REAL, INTENT(IN)               :: dels ! time setp size (s)
    
   INTEGER :: k,kk,j  

CHARACTER(LEN=*), PARAMETER :: subr_name = "cbl_model_driver"
LOGICAL :: jls_standalone = .FALSE.
LOGICAL :: jls_radiation  = .FALSE.
LOGICAL :: cbl_standalone = .FALSE.    

   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)    

      cable_runtime%um_radiation = .FALSE.
      
      IF( cable_runtime%um_explicit ) THEN
!d1!      met%DoY = met%DoY + 1.
        CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )
         met%tk = met%tk + C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
      ENDIF
      
      CALL define_air (met, air)
   
call fveg_mask( veg_mask, mp, Clai_thresh, canopy%vlaiw )
call fsunlit_mask( sunlit_mask, mp, CRAD_THRESH,( met%fsd(:,1)+met%fsd(:,2) ) )
call fsunlit_veg_mask( sunlit_veg_mask, mp )

!Ticket 334 - as ESM1.5 ensure that it uses the original kbeam at low sun angles
cable_user%use_new_beam_coef = .FALSE.
CALL init_radiation( &
                     rad%extkb, rad%extkd,                                     &
                     !ExtCoeff_beam, ExtCoeff_dif,                             &
                     rad%extkbm, rad%extkdm, Rad%Fbeam,                        &
                     !EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,             &
                     c1, rhoch, xk,                                            &
                     mp,nrb,                                                   &
                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,         &
                     cbl_standalone, jls_standalone, jls_radiation,            &
                     cable_user%use_new_beam_coef, subr_name,                  &
                     veg_mask, sunlit_mask, sunlit_veg_mask,                   &
                     veg%Xfang, veg%taul, veg%refl,                            &
                     !VegXfang, VegTaul, VegRefl                               &
                     met%coszen, int(met%DoY), met%fsd,                        &
                     !coszen, metDoY, SW_down,                                 &
                     canopy%vlaiw  ) !reducedLAIdue2snow 
 
IF( cable_runtime%um_explicit ) THEN
   
 !Ticket 331 refactored albedo code for JAC
 CALL snow_aging(ssnow%snage,mp,dels,ssnow%snowd,ssnow%osnowd,ssnow%tggsn(:,1),&
         ssnow%tgg(:,1),ssnow%isflag,veg%iveg,soil%isoilm) 

 !Ticket 334 as ESM force to use old scheme
 cable_user%limit_all_exp = .FALSE.
 CALL Albedo( ssnow%AlbSoilsn, soil%AlbSoil,                                 &
             !AlbSnow, AlbSoil,              
             mp, nrb,                                                       &
             jls_radiation, cable_user%limit_all_exp, Cmax_kLAI,            &
             veg_mask, sunlit_mask, sunlit_veg_mask,                        &  
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
             rad%extkd, rad%extkb,                                          & 
             !ExtCoeff_dif, ExtCoeff_beam,
             rad%extkdm, rad%extkbm,                                        & 
             !EffExtCoeff_dif, EffExtCoeff_beam,                
             rad%rhocdf, rad%rhocbm,                                        &
             !CanopyRefl_dif,CanopyRefl_beam,
             rad%cexpkdm, rad%cexpkbm,                                      & 
             !CanopyTransmit_dif, CanopyTransmit_beam, 
             rad%reffdf, rad%reffbm                                        &
           ) !EffSurfRefl_dif, EffSurfRefl_beam 

      ENDIF
   
   ! Calculate canopy variables:
   CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy, climate_cbl, sunlit_veg_mask, canopy%vlaiw, cable_user%limit_all_exp, Cmax_kLAI )

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


