MODULE cbl_model_driver_mod

IMPLICIT NONE

PRIVATE
PUBLIC cbl_model_driver

CONTAINS

SUBROUTINE cbl_model_driver( explicit_path, mp,nrb, land_pts, npft, ktau,dels, air,          &
                    bgc, canopy, met,                                         &
                    bal, rad, rough, soil,                                    &
                    ssnow, sum_flux, veg, z0surf_min,                         &
                    LAI_pft, HGT_pft, RmetDoY, reducedLAIdue2snow )

!subrs
USE cbl_albedo_mod, ONLY: albedo
USE cbl_masks_mod, ONLY: fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
USE cbl_masks_mod, ONLY: veg_mask,  sunlit_mask,  sunlit_veg_mask

!data
USE cable_other_constants_mod,  ONLY: Ccoszen_tols => coszen_tols

!jhan:pass these
USE cable_other_constants_mod, ONLY: clai_thresh => lai_thresh

USE cable_other_constants_mod, ONLY: cgauss_w => gauss_w
USE cable_math_constants_mod, ONLY: cpi => pi
USE cable_math_constants_mod, ONLY: cpi180 => pi180

USE cable_common_module
USE cable_carbon_module
USE cable_air_type_mod,       ONLY: air_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_params_mod,         ONLY: soil_parameter_type
USE cable_def_types_mod,      ONLY: climate_type
!icycle switch. [0]=None,1=C,2=CN,3=CNP
USE casadimension,            ONLY : icycle

!subrs
USE cbl_soil_snow_main_module,  ONLY: soil_snow
USE cable_roughness_module, ONLY: ruff_resist
USE cbl_init_radiation_module, ONLY: init_radiation

USE cable_air_module, ONLY: define_air
   
USE cable_phys_constants_mod, ONLY : csboltz => sboltz
USE cable_phys_constants_mod, ONLY : cemleaf => emleaf
USE cable_phys_constants_mod, ONLY : cemsoil => emsoil
   
USE cable_canopy_module, ONLY: define_canopy
USE cable_canopy_module_explicit, ONLY: define_canopy_explicit
                                   
LOGICAL :: explicit_path
INTEGER :: mp
INTEGER :: nrb
INTEGER :: land_pts, npft
REAL :: AlbSoil(mp,nrb)
REAL :: MetTk(mp)
REAL :: SnowDepth(mp) 
REAL :: SnowDensity(mp)
REAL :: SnowODepth(mp)
REAL :: SoilTemp(mp)
REAL :: SnowAge(mp)
INTEGER:: SnowFlag_3L(mp)
INTEGER:: surface_type(mp)

REAL :: z0surf_min
REAL  :: hgt_pft(mp)
REAL  :: lai_pft(mp) 
REAL :: RmetDoY(mp)          !Day of the Year [formerly met%doy]
INTEGER :: metDoY(mp)          !Day of the Year [formerly met%doy]
  
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
TYPE (climate_type) :: climate

TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg

REAL, INTENT(IN)               :: dels ! time setp size (s)
INTEGER, INTENT(IN) :: ktau
LOGICAL, SAVE :: first_call = .TRUE.

CHARACTER(LEN=*), PARAMETER :: subr_name = "cbl_model_driver"
LOGICAL :: jls_standalone= .TRUE.
LOGICAL :: jls_radiation= .FALSE.
LOGICAL :: cbl_standalone = .FALSE.    

!make local to rad_driver and also again in cbl_model_driver
!CABLE variables to keep for all CABLE pathways across the timestep 
REAL :: reducedLAIdue2snow(mp)

!co-efficients usoughout init_radiation ` called from _albedo as well - need to come from toplevel
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)
REAL :: CanopyRefl_dif(mp,nrb)
REAL :: CanopyRefl_beam(mp,nrb)

!is this even nencessary anymore?
IF ( first_call )  met%ca = 3.999999e-4 
      
metDoy = INT(RmetDoy)
!iFor testing
!ICYCLE = 0
cable_user%soil_struc="default"

CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )
!H!CALL ruff_resist( veg, rough, ssnow, canopy, LAI_pft, HGT_pft, reducedLAIdue2snow )
reducedLAIdue2snow = canopy%vlaiw
!jhan: this call to define air may be redundant
CALL define_air (met, air)

call fveg_mask( veg_mask, mp, Clai_thresh, canopy%vlaiw )
call fsunlit_mask( sunlit_mask, mp, Ccoszen_tols, met%coszen )
call fsunlit_veg_mask( sunlit_veg_mask, mp )

CALL init_radiation( rad%extkb, rad%extkd,                                    &
                     !ExtCoeff_beam, ExtCoeff_dif,
                     rad%extkbm, rad%extkdm, Rad%Fbeam,                       &
                     !EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,
                     c1, rhoch, xk,                                           &
                     mp,nrb,                                                  &
                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,        &
                     cbl_standalone, jls_standalone, jls_radiation,           &
                     subr_name,                                               &
                     veg_mask, sunlit_mask, sunlit_veg_mask,                  &
                     veg%Xfang, veg%taul, veg%refl,                           &
                     !VegXfang, VegTaul, VegRefl
                     met%coszen, INT(met%DoY), met%fsd,                       &
                     !coszen, metDoY, SW_down,
                     canopy%vlaiw                                              &
                   ) !reducedLAIdue2snow 
 
!!IF ( cable_runtime%um_explicit )                                              &
CALL Albedo( ssnow%AlbSoilsn, soil%AlbSoil,                                 &
             !AlbSnow, AlbSoil,              
             mp, nrb,                                                       &
             jls_radiation,                                                 &
             veg_mask, sunlit_mask, sunlit_veg_mask,                        &  
             Ccoszen_tols, cgauss_w,                                        & 
             veg%iveg, soil%isoilm,veg%refl, veg%taul,                      & 
             !surface_type, VegRefl, VegTaul,
             met%tk, met%coszen, canopy%vlaiw,                              &
             !metTk, coszen, reducedLAIdue2snow,
             ssnow%snowd, ssnow%osnowd, ssnow%isflag,                       & 
             !SnowDepth, SnowODepth, SnowFlag_3L, 
             ssnow%ssdnn, ssnow%tgg(:,1), ssnow%tggsn(:,1), ssnow%snage,                      & 
             !SnowDensity, SoilTemp, SnowAge, 
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

!CABLE_LSM:check
!IF ( explicit_path ) THEN
IF ( ktau_gl == 1) THEN
  ssnow%tss=(1 - ssnow%isflag) * ssnow%tgg(:,1) + ssnow%isflag * ssnow%tggsn(:,1) 
  ssnow%otss = ssnow%tss
  first_call = .FALSE.
ENDIF
ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
ssnow%otss = ssnow%tss

IF ( explicit_path ) THEN
  CALL define_canopy_explicit(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, sunlit_veg_mask,  canopy%vlaiw)
ELSE
  CALL define_canopy         ( bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, &
  sunlit_veg_mask, reducedLAIdue2snow )
END IF

ssnow%owetfac = ssnow%wetfac

IF ( cable_runtime%um_implicit ) &
  CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)

ssnow%deltss = ssnow%tss - ssnow%otss

! need to adjust fe after soilsnow
canopy%fev  = canopy%fevc + canopy%fevw

! Calculate total latent heat flux:
canopy%fe = canopy%fev + canopy%fes

! Calculate net radiation absorbed by soil + veg
canopy%rnet = canopy%fns + canopy%fnv

! Calculate radiative/skin temperature:
rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 + rad%transd * ssnow%tss**4 )**0.25
! Calculate radiative/skin temperature:
       rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
            rad%transd * ssnow%tss**4 )**0.25

IF (.NOT.cable_runtime%um_explicit .AND. icycle == 0) THEN

  !calculate canopy%frp
  CALL plantcarb(veg,bgc,met,canopy)

  !calculate canopy%frs
  CALL soilcarb(soil, ssnow, veg, bgc, met, canopy)

  CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

  canopy%fnpp = -1.0* canopy%fpn - canopy%frp
  canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

ENDIF

END SUBROUTINE cbl_model_driver

END MODULE cbl_model_driver_mod


