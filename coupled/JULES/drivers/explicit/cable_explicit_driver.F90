MODULE cable_explicit_driv_mod
  
CONTAINS

SUBROUTINE cable_explicit_driver(                                             &
!corresponding name (if differs) of varaible on other side of call/subroutine shown in "[]" 

!Variables to be calculated and returned by CABLE
!------------------------------------------------------------------------------
ftl_tile, & !surface sensible heat flux  [W/m2]? (up=+ve?) -"FTL_tile" in CABLE
fqw_tile, & !surface moisture flux flux  [kg/m^2/s/m2](up=+ve?) units could be changed?
tstar_tile,         & !surface temperature [K] per tile
u_s,                & ! land point surface friction velocity [m/s]
u_s_std_tile,       & ! per tile surface friction velocity [m/s] canopy%us
cd_tile,                                                                      &
ch_tile,                                                                      &
radnet_tile,        & !Net radiation at surface [W/m2]
fraca,& !Fraction of surface moisture flux with only aerodynamic resistance for
                      !snow-free land tiles.
resfs,& !Combined soil, stomatal and aerodynamic resistance factor for fraction
                      !(1-FRACA) of snow-free land tiles.
resft, & !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1 for snow.
Z0H_tile,           & ! Tile roughness lengths for heat and moisture (m).
Z0M_tile,           & ! OUT Tile roughness lengths for momentum.
RECIP_L_MO_tile,   & ! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
EPOT_tile,          & ! Potential evaporation from surface, per tile
!------------------------------------------------------------------------------

!This is an "outlier" and possibly misleading. CABLE does not actually calculate radiation:
!Generally we speak of 4-band radiation. This is actually only 2-bands VIS/NIR in the SW.
!We further split each of these into Direct Beam and Diffuse components. 
!Offline CABLE splits the total SW forcing into VIS/NIR using a Spitter() fuction
!We include this variable here to connect back to JULES toplevel routines because:
!Online the UM radiation scheme DOES compute surf_down_sw using a more sophisticated model 
!than that which we use Offline, however not until AFTER the surface albedos have been 
!calculated which IS what is done AND here and technically does not require knowledge of 
!the downward SW. JULES aggregates this SW and threads this to the LSM as it is called 
!explicitly. 
!------------------------------------------------------------------------------
surf_down_sw,     & ! ShortWave radiation per rad band (row_length,rows,4) 
!------------------------------------------------------------------------------

!Mostly model dimensions and associated
!------------------------------------------------------------------------------
row_length,         & !grid cell x
rows,               & !grid cell y
land_pts,           & !grid cell land points on the x,y grid
ntiles,             & !grid cell number of surface types [] 
sm_levels,          & !grid cell number of soil levels 
npft,               & !grid cell number of PFTs 
tile_pts,           & !Number of land points per PFT [] 
tile_index,         & !Index of land point in (land_pts) array[] 
land_index, & !Index of land points in (x,y) array - see  corresponding *decs.inc
timestep_width,     & !bin width in seconds of timestep
endstep,            & !last timestep of experiment
timestep_number,                                                              &
doy,                                                                          &
mp,                                                                           &
nrb,                                                                          &
!------------------------------------------------------------------------------

!Surface descriptions generally parametrized
!------------------------------------------------------------------------------
Fland,              & !fraction of land per land point (could be coastal) 
tile_frac,         & !fraction of each surface type per land point [frac_surft]
L_tile_pts,       & !Logical mask TRUE where tile frac > 0. used to PACK/UNPACK
LAI_pft_um,             & !Leaf area index. [LAI_pft/LAI_pft_um in radiation]
canht_pft_um,           & !Canopy height [canht_pft/HGT_pft_um in radiation]
albsoil,    & !(albsoil)Snow-free, bare soil albedo [albsoil_soilt(:,1) in um ]
z0surf_min,                                                                   &
dzsoil,             & !soil thicknesses in each layer  
bexp,                                                                         &
hcon,                                                                         &
satcon,                                                                       &
sathh,                                                                        &
smvcst,                                                                       &
smvcwt,                                                                       &
smvccl,                                                                       &
!------------------------------------------------------------------------------

!Variables passed from JULES/UM
!------------------------------------------------------------------------------
latitude,           & !latitude
longitude,          & !longitude
cosine_zenith_angle,& ! cosine_zenith_angle [cosz_ij]
sin_theta_latitude,                                                           &
sthu,                                                                         &
lw_down,                                                                      &
ls_rain,                                                                      &
ls_snow,                                                                      &
tl_1,                                                                         &
qw_1,                                                                         &
vshr_land,                                                                    &
pstar,                                                                        &
z1_tq,                                                                        &
z1_uv,                                                                        &
snow_tile,          & !snow depth equivalent (in water?) [snow_surft]
                      !This is the total snow depth per tile. CABLE also has depth per layer
canopy_tile,                                                                  &
co2_mmr,                                                                      &
!------------------------------------------------------------------------------

!CABLE prognostics
!------------------------------------------------------------------------------
iThreeLayerSnowFlag, & ! flag indicating whether enough snow to treat as 3 layers
                      ! [real(ThreeLayerSnowFlag_CABLE)]
OneLyrSnowDensity, & ! density considering snow as 1 layer [OneLyrSnowDensity_CABLE
SnowAge,                                                                      &
SnowDensity,                                                                  &
SnowDepth,                                                                    &
SnowMass,                                                                     &
SnowTemp,                                                                     &
SoilMoisture,                                                                 &
FrozenSoilFrac,                                                               &
SoilTemp,                                                                     &
!The simplest way to replace old CABLE types with these new ones issimply to 
!del the "_cbl" tag from the suffix 
air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
soilin, soil_cbl )

!subrs called 
USE cable_um_init_mod, ONLY: interface_UM_data
USE cbl_model_driver_mod, ONLY: cbl_model_driver
!data
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
USE cable_params_mod,         ONLY: soilin_type
USE cable_params_mod,         ONLY: soil_parameter_type
USE cable_phys_constants_mod, ONLY : cdensity_liq => density_liq

USE cable_runtime_opts_mod ,ONLY: cable_user
USE cable_runtime_opts_mod ,ONLY: satuparam
USE cable_runtime_opts_mod ,ONLY: wiltparam

  !processor number, timestep number / width, endstep
USE cable_common_module, ONLY: knode_gl, ktau_gl, kwidth_gl, kend_gl
USE cable_common_module, ONLY: cable_runtime
!--- vars common to CABLE declared 
   
USE cable_def_types_mod, ONLY: ms

!USE casavariable
!USE casa_types_mod
!USE cable_climate_mod
!made this a module so can build in the UM re:dependencies etc
!did the same dor sli_main. promote everything to modules
!USE casa_cable
   
IMPLICIT NONE
REAL :: z0surf_min
  
  !___ re-decl input args
   
  !___IN: UM dimensions, array indexes, flags
INTEGER ::                                                                    &
  row_length, rows, & ! UM grid resolution
  land_pts,         & ! # of land points being processed
  ntiles,           & ! # of tiles 
  npft,             & ! # of plant functional types
  sm_levels           ! # of soil layers 

INTEGER :: mp, nrb
! index of land points being processed
INTEGER, DIMENSION(land_pts) :: land_index 

! # of land points on each tile
INTEGER,  DIMENSION(ntiles) :: tile_pts 
  
LOGICAL,DIMENSION(land_pts, ntiles) ::                                        &
   L_tile_pts  ! true IF vegetation (tile) fraction is greater than 0
  
INTEGER,  DIMENSION(land_pts, ntiles) ::                                      &
  tile_index ,& ! index of tile points being processed
  iThreeLayerSnowFlag   ! 3 layer snow flag

!___UM parameters: water density, soil layer thicknesses 
REAL,  DIMENSION(sm_levels) :: dzsoil

!___UM soil/snow/radiation/met vars
REAL,  DIMENSION(land_pts) ::                                                 &
  bexp,    & ! => parameter b in Campbell equation 
  hcon,    & ! Soil thermal conductivity (W/m/K).
  satcon,  & ! hydraulic conductivity @ saturation [mm/s]
  sathh,                                                                      &
  smvcst,                                                                     &
  smvcwt,                                                                     &
  smvccl,                                                                     &
  albsoil,                                                                    &
  fland 
  
REAL,  DIMENSION(land_pts) ::                                                 &
  slope_avg,                                                                  &
  slope_std,                                                                  &
    dz_gw,sy_gw,perm_gw,drain_gw

  
REAL,  DIMENSION(row_length,rows) ::                                          &
  cosine_zenith_angle,                                                        &
  latitude,                                                                   &
  longitude,                                                                  &
  sw_down,    & ! NOT SW forcing. surf_down_sw IS 
  lw_down,                                                                    &
  ls_rain,                                                                    &
  ls_snow,                                                                    &
  tl_1,                                                                       &
  qw_1,                                                                       &
  vshr_land,                                                                  &
  pstar,                                                                      &
  z1_tq,                                                                      &
  z1_uv
  
REAL,  DIMENSION(land_pts, ntiles) ::                                         &
!REAL, DIMENSION(row_length, rows, 4) ::                         &
  surf_down_sw 
   
REAL,  DIMENSION(land_pts, ntiles) ::                                         &
  snow_tile,                                                                  &
  tile_frac,                                                                  &
  OneLyrSnowDensity,                                                          &
  SnowAge,                                                                    &
  canopy_tile,                                                                &
  visc_sublayer_depth 

REAL, DIMENSION(land_pts, npft) ::                                            &
  canht_pft_um, lai_pft_um 
   
REAL,  DIMENSION(land_pts, ntiles,3) ::                                       &
  snow_cond,                                                                  &
  SnowDensity,                                                                &
  SnowDepth,                                                                  &
  SnowMass,                                                                   &
  SnowTemp
   
REAL, DIMENSION(land_pts, sm_levels) ::                                       &
  sthu 
   
REAL, DIMENSION(land_pts, ntiles, sm_levels) ::                               &
  sthu_tile,                                                                  &
  FrozenSoilFrac,                                                             &
  SoilMoisture,                                                               &
  SoilTemp
   
REAL :: co2_mmr
   
REAL :: sin_theta_latitude(row_length,rows) 
    
!___return fluxes AND miscelaneous 
REAL, DIMENSION(land_pts,ntiles) ::                                           &
  ftl_tile,   &  ! Surface FTL for land tiles     
  fqw_tile,   &  ! Surface FQW for land tiles     
  tstar_tile,                                                                 &
  z0h_tile,                                                                   &
  z0m_tile,                                                                   &
  cd_tile,    &     ! Drag coefficient
  ch_tile,    &     ! Transfer coefficient for heat & moisture
  u_s_std_tile, &      ! Surface friction velocity
  radnet_tile,   & ! Surface net radiation
  resfs,         & ! Combined soil, stomatal & aerodynamic resistance
                   ! factor for fraction (1-FRACA) of snow-free land tiles
  resft,         & ! Total resistance factor.
                   ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                   ! 1 for snow.    
  fraca,         & ! Fraction of surface moisture
  recip_l_mo_tile,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
  epot_tile,                                                                  &
  smgw_tile

REAL, DIMENSION(row_length,rows)  ::                                          &
  u_s               ! Surface friction velocity (m/s)
   
! end step of experiment, this step, step width, processor num
INTEGER :: endstep, timestep_number, mype
REAL ::  timestep_width     
   
TYPE(air_type),       INTENT(INOUT)  :: air_cbl
TYPE(met_type),       INTENT(INOUT)  :: met_cbl
TYPE(radiation_type),       INTENT(INOUT)  :: rad_cbl
TYPE(roughness_type),     INTENT(INOUT)  :: rough_cbl
TYPE(canopy_type),    INTENT(INOUT)  :: canopy_cbl
TYPE(soil_snow_type),     INTENT(INOUT)  :: ssnow_cbl
TYPE(bgc_pool_type),       INTENT(INOUT)  :: bgc_cbl
TYPE(balances_type),       INTENT(INOUT)  :: bal_cbl
TYPE(sum_flux_type),  INTENT(INOUT)  :: sum_flux_cbl
TYPE(veg_parameter_type),   INTENT(INOUT) :: veg_cbl
TYPE(soilin_type),  INTENT(INOUT) ::  soilin  
TYPE(soil_parameter_type),  INTENT(INOUT) ::  soil_cbl

  !CASA vars 
REAL, DIMENSION(land_pts,ntiles,10) ::                                        &
  cpool_tile,    & ! Carbon Pools
  npool_tile       ! Nitrogen Pools

REAL, DIMENSION(land_pts,ntiles,12) ::                                        &
  ppool_tile       ! Phosphorus Pools

REAL, DIMENSION(land_pts) ::                                                  &
  soil_order,    & ! Soil Order (1 to 12)
  nidep,         & ! Nitrogen Deposition
  nifix,         & ! Nitrogen Fixation
  pwea,          & ! Phosphorus from Weathering
  pdust            ! Phosphorus from Dust

REAL, DIMENSION(land_pts,ntiles) ::                                           &
  glai,         &  ! Leaf Area Index for Prognostics LAI
  phenphase,    &  ! Phenology Phase for Casa-CNP
  npp_ft_acc,                                                                 &
  resp_w_ft_acc

INTEGER :: doy
  
  !___ local vars
CHARACTER(LEN=200), PARAMETER ::                                              &
    runtime_vars_file = 'cable.nml'! path/name namelist def runtime vars

!___ 1st call in RUN (!=ktau_gl -see below) 
LOGICAL, SAVE :: first_cable_call = .TRUE.

REAL  :: rho_water, rho_ice
  
! std template args 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_explicit_driver"
LOGICAL, PARAMETER :: explicit_path = .TRUE.

rho_water = cdensity_liq
rho_ice   = cdensity_liq
!---------------------------------------------------------------------!
!--- initialize CABLE using UM forcings etc. these args are passed ---!
!--- down from UM.                                                 ---! 
!---------------------------------------------------------------------!
!H!pass cable% types
!HACK
dz_gw (:) = 20.0
slope_avg(:) = 0.02
slope_std(:) = .005 
drain_gw(:) = 0.0008
perm_gw(:) = 3.0e-6
npp_ft_acc(:,:)  = 0.0 
resp_w_ft_acc(:,:)  = 0.0 

CALL interface_UM_data( mp, row_length, rows, land_pts, ntiles, npft,         &
                        sm_levels, ktau_gl, latitude, longitude,              &
                        land_index, tile_frac, tile_pts, tile_index,          &
                        bexp, hcon, satcon, sathh, smvcst, smvcwt,            &
                        smvccl, albsoil, slope_avg,slope_std,dz_gw,           &
                        perm_gw,drain_gw,snow_tile, OneLyrSnowDensity ,       &
                        SnowAge, iThreeLayerSnowFlag, SnowDensity, snow_cond, &
                        SnowDepth, SnowTemp, SnowMass, sw_down,               &
                        lw_down, cosine_zenith_angle, surf_down_sw, ls_rain,  &
                        ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,         &
                        z1_uv, rho_water, rho_ice,L_tile_pts,                 &
                        visc_sublayer_depth, canopy_tile, Fland,              &
                        co2_mmr,                                              &
                        sthu_tile, SoilMoisture, smgw_tile,FrozenSoilFrac,    &
                        sthu, SoilTemp, canht_pft_um, lai_pft_um,                     &
                        sin_theta_latitude, dzsoil,                           &
                        cpool_tile, npool_tile, ppool_tile, soil_order,       &
                        nidep, nifix, pwea, pdust, glai, phenphase,           &
                        npp_ft_acc,resp_w_ft_acc,                             &
air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
soilin, soil_cbl )

met_cbl%doy = REAL(doy)
canopy_cbl%oldcansto = canopy_cbl%cansto
rad_cbl%otrad = rad_cbl%trad
  !---------------------------------------------------------------------!
  !--- cbm "mainly" controls the calling of model components         ---!  
  !---------------------------------------------------------------------!
CALL cbl_model_driver( explicit_path, mp, nrb, land_pts, npft, ktau_gl,       &
          timestep_width,      &
          air_cbl ,                                                           &
          bgc_cbl,                                                            &
          canopy_cbl ,                                                        &
          met_cbl ,                                                           &
          bal_cbl ,                                                           &
          rad_cbl,                                                            &
          rough_cbl,                                                          &
          soil_cbl ,                                                          &
          ssnow_cbl ,                                                         &
          sum_flux_cbl ,                                                      &
          veg_cbl ,                                                           &
          z0surf_min,                                                         &
          !H!shouuld already work from here LAI_pft, HGT_pft )
          !H!veg_cbl %vlai, veg_cbl %hc, met_cbl %Doy, canopy_cbl %vlaiw )
          veg_cbl %vlai, veg_cbl %hc, met_cbl %Doy, canopy_cbl %vlaiw )

RETURN

END SUBROUTINE cable_explicit_driver

END MODULE cable_explicit_driv_mod
