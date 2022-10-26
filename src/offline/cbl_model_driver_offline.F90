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

MODULE cable_cbm_module

IMPLICIT NONE

PRIVATE
PUBLIC cbm

CONTAINS

SUBROUTINE cbm( ktau,dels, air, bgc, canopy, met,                                &
       bal, rad, rough, soil,                                      &
       ssnow, sum_flux, veg, climate, xk, c1, rhoch )

    USE cable_common_module
    USE cable_carbon_module
    USE cbl_soil_snow_main_module, ONLY : soil_snow
    USE cable_def_types_mod
    USE cable_roughness_module, ONLY : ruff_resist
    USE cbl_init_radiation_module, ONLY : init_radiation
    USE cable_air_module, ONLY : define_air
    USE casadimension,     ONLY : icycle ! used in casa_cnp
! physical constants
USE cable_phys_constants_mod, ONLY : CGRAV  => GRAV
USE cable_phys_constants_mod, ONLY : CCAPP   => CAPP
USE cable_phys_constants_mod, ONLY : CEMLEAF => EMLEAF
USE cable_phys_constants_mod, ONLY : CEMSOIL => EMSOIL
USE cable_phys_constants_mod, ONLY : CSBOLTZ => SBOLTZ
    !mrd561
    USE cable_gw_hydro_module, ONLY : sli_hydrology,&
         soil_snow_gw
    USE cable_canopy_module, ONLY : define_canopy
    USE cbl_albedo_mod, ONLY : albedo
    USE sli_main_mod, ONLY : sli_main
    USE snow_aging_mod,               ONLY: snow_aging
    
!data !jhan:pass these
USE cable_other_constants_mod, ONLY : CLAI_THRESH => lai_thresh
USE cable_other_constants_mod,  ONLY : Crad_thresh => rad_thresh
USE cable_other_constants_mod,  ONLY : Ccoszen_tols => coszen_tols
USE cable_other_constants_mod, ONLY : CGAUSS_W => gauss_w
USE cable_math_constants_mod, ONLY : CPI => pi
USE cable_math_constants_mod, ONLY : CPI180 => pi180
use cbl_masks_mod, ONLY :  fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
use cbl_masks_mod, ONLY :  veg_mask,  sunlit_mask,  sunlit_veg_mask

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
character(len=*), parameter :: subr_name = "cbm"
LOGICAL :: cbl_standalone= .true.
LOGICAL :: jls_standalone= .false.
LOGICAL :: jls_radiation= .false.

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)

!iFor testing
cable_user%soil_struc="default"

CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )

!jhan: this call to define air may be redundant
CALL define_air (met, air)

call fveg_mask( veg_mask, mp, Clai_thresh, canopy%vlaiw )
!call fsunlit_mask( sunlit_mask, mp, Ccoszen_tols, met%coszen )
call fsunlit_mask( sunlit_mask, mp, CRAD_THRESH,( met%fsd(:,1)+met%fsd(:,2) ) )
call fsunlit_veg_mask( sunlit_veg_mask, mp )

CALL init_radiation( rad%extkb, rad%extkd,                                     &
                     !ExtCoeff_beam, ExtCoeff_dif,
                     rad%extkbm, rad%extkdm, Rad%Fbeam,                        &
                     !EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,
                     c1, rhoch, xk,                                            &
                     mp,nrb,                                                   &
                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,         &
                     cbl_standalone, jls_standalone, jls_radiation,            &
                     subr_name,                                                &
                     veg_mask, sunlit_mask, sunlit_veg_mask,                   &
                     veg%Xfang, veg%taul, veg%refl,                            &
                     !VegXfang, VegTaul, VegRefl
                     met%coszen, int(met%DoY), met%fsd,                        &
                     !coszen, metDoY, SW_down,
                     canopy%vlaiw                                              &
                   ) !reducedLAIdue2snow 

!Ticket 331 refactored albedo code for JAC
CALL snow_aging(ssnow%snage,mp,dels,ssnow%snowd,ssnow%osnowd,ssnow%tggsn(:,1),&
         ssnow%tgg(:,1),ssnow%isflag,veg%iveg,soil%isoilm) 

call Albedo( ssnow%AlbSoilsn, soil%AlbSoil,                                &
             !AlbSnow, AlbSoil,              
             mp, nrb,                                                      &
             jls_radiation,                                                &
             veg_mask, sunlit_mask, sunlit_veg_mask,                       &  
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
             rad%extkd, rad%extkb,                                         & 
             !ExtCoeff_dif, ExtCoeff_beam,
             rad%extkdm, rad%extkbm,                                       & 
             !EffExtCoeff_dif, EffExtCoeff_beam,                
             rad%rhocdf, rad%rhocbm,                                       &
             !CanopyRefl_dif,CanopyRefl_beam,
             rad%cexpkdm, rad%cexpkbm,                                     & 
             !CanopyTransmit_dif, CanopyTransmit_beam, 
             rad%reffdf, rad%reffbm                                        &
           ) !EffSurfRefl_dif, EffSurfRefl_beam 

ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
ssnow%otss = ssnow%tss

CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, sunlit_veg_mask,  canopy%vlaiw)
    
ssnow%owetfac = ssnow%wetfac

CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)

ssnow%deltss = ssnow%tss-ssnow%otss

    ! need to adjust fe after soilsnow
    canopy%fev  = canopy%fevc + canopy%fevw

    ! Calculate total latent heat flux:
    canopy%fe = canopy%fev + canopy%fes

    ! Calculate net radiation absorbed by soil + veg
    canopy%rnet = canopy%fns + canopy%fnv

    ! Calculate radiative/skin temperature:
    rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
            rad%transd * ssnow%tss**4 )**0.25
    IF (icycle == 0) THEN

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
