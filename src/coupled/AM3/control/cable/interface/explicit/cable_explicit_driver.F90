MODULE cable_explicit_driv_mod
  
CONTAINS

SUBROUTINE cable_explicit_driver(                                              & 
      ! IN: UM/JULES/CABLE model/grid parameters, fields, mappings
      mype, row_length, rows, land_pts, nsurft, npft, sm_levels, dzsoil,       &
      timestep, timestep_number, mp, nrb, land_index, surft_pts, surft_index,  &
      l_tile_pts, latitude, longitude, cos_zenith_angle, Fland, tile_frac,     &
      
      ! IN: soil parameters !1 is only allowable index in UM
      bexp, hcon, satcon, sathh, smvcst, smvcwt, smvccl, albsoil,              &
      
      ! IN: SW forcing: manipulated for CABLE 
      sw_down_VIS, sw_down_NIR, beamFrac_VIS, beamFrac_NIR, beamFrac_TOT,      &
      
      ! IN: Met forcing: 
      lw_down, ls_rain, ls_snow,                                               &
      tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy_tile,                 &
      ! This an outlier IN here. INOUT @ implicit. (was)OUT at extras
      ! I think we are dealing with it OK now but confusion could be removed  
      snow_tile,                                                               &
 
      ! IN: canopy height, LAI seasonally presecribed, potentially prognostic 
      ! IN: CO2 mass mixing ratio  
      canht_pft, lai_pft, CO2_MMR,                                             &
      
      ! IN: carries vegin/soilin - potentially redundant  
      pars,                                                                    &
      
      ! IN: tiled soil/snow prognostics - IN here. INOUT @ implicit  
      progs_soiltemp, progs_soilmoisture, progs_FrozenSoilFrac,                &
      progs_ThreeLayerSnowFlag, progs_SnowDepth, progs_SnowMass,               &
      progs_SnowTemp, progs_SnowDensity, progs_snowage, progs_snowosurft,      &
      progs_OneLyrSnowDensity,                                                 &  
      
      ! INOUT: CABLE TYPEs roughly grouped fields per module 
      rad, met, rough, canopy, veg, soil, ssnow, bal, air, bgc, sum_flux,      &
      
      !OUT: currently being passed back to UM in veg%iveg, soil%isoilm      
      SurfaceType, SoilType,                                                   &
      !OUT: currently being passed back to UM in veg%hc, veg%vlai 
      HGT_pft_cbl, LAI_pft_cbl,                                                &
      
      !IN: currently being passed from prev radiation call through work%
      ! jhan:quirky, snow (in turn reduced LAI due to snow) can evolve through a 
      ! constant rad dt. However reducedLAIdue2snow used ubiquitously as trigger 
      ! Further, snow does NOT evolve in explicit AND reducedLAIdue2snow absent 
      ! in implicit
      reducedLAIdue2snow,                                                      & 
      
          !GW
          !visc_sublayer_depth, smgw_tile, slope_avg, slope_std,
          !dz_gw, perm_gw, drain_gw,                           
          !casa progs
          !CPOOL_TILE, NPOOL_TILE, PPOOL_TILE, SOIL_ORDER, NIDEP,
          !NIFIX, PWEA, PDUST, GLAI, PHENPHASE,
      
      !IN: if not passed a dangling argument would ensue 
      npp_pft_acc, resp_w_pft_acc )

! subrs called 
USE cbl_um_init_mod,   ONLY: init_data
USE cbl_um_update_mod, ONLY: update_data
USE cable_cbm_module,  ONLY: cbm_expl

! data
USE grid_constants_mod_cbl,   ONLY: ICE_SoilType, nsl, nsnl
USE cable_phys_constants_mod, ONLY: density_liq, density_ice, tfrz
USE cable_surface_types_mod,  ONLY: ICE_SurfaceType => ICE_cable


USE params_io_mod_cbl, ONLY: params_io_data_type
USE params_io_mod_cbl, ONLY: params_io_type


USE cable_def_types_mod, ONLY : climate_type
USE cable_def_types_mod, ONLY : met_type, radiation_type, veg_parameter_type,  &
                                soil_parameter_type, roughness_type,           &
                                canopy_type, soil_snow_type, balances_type,    &
                                air_type, bgc_pool_type, sum_flux_type

!--- processor number, timestep number, timestep width !ultimately get rid of these - pass %runtime through parent
USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, cable_runtime, cable_user, redistrb, satuparam,wiltparam
!block!USE casavariable
!block!USE casa_types_mod

IMPLICIT NONE
  
INTEGER, INTENT(IN) :: mype
INTEGER, INTENT(IN) :: timestep_number     
REAL,    INTENT(IN) :: timestep
INTEGER, INTENT(IN) :: row_length           ! # columns in spatial grid
INTEGER, INTENT(IN) :: rows                 ! # rows in spatial grid
INTEGER, INTENT(IN) :: land_pts             ! # land points
INTEGER, INTENT(IN) :: nsurft               ! # tiles 
INTEGER, INTENT(IN) :: npft                 ! # plant functional types
INTEGER, INTENT(IN) :: sm_levels            ! # soil layers 
REAL,    INTENT(IN) :: dzsoil(sm_levels)    ! soil layer thicknesses 
INTEGER, INTENT(IN) :: mp                   ! # active land points
INTEGER, INTENT(IN) :: nrb                  ! # radiation bands
INTEGER, INTENT(IN) :: surft_pts(nsurft)    ! # land points on each tile
INTEGER, INTENT(IN) :: land_index(land_pts) ! index of land points
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! index of tile points
LOGICAL, INTENT(IN) :: l_tile_pts(land_pts, nsurft)

REAL, INTENT(IN) :: canht_pft(land_pts, npft) 
REAL, INTENT(IN) :: lai_pft(land_pts, npft)
REAL, INTENT(IN) :: fland(land_pts)

REAL, INTENT(IN) :: co2_mmr
REAL, INTENT(IN) :: tile_frac(land_pts, nsurft)
REAL, INTENT(IN) :: cos_zenith_angle(row_length,rows)
   
REAL, INTENT(IN) :: latitude(row_length,rows)
REAL, INTENT(IN) :: longitude(row_length,rows)

REAL, INTENT(IN) :: bexp (land_pts, sm_levels)
        ! => parameter b in Campbell equation 
REAL, INTENT(IN) :: satcon(land_pts, sm_levels)
       ! hydraulic conductivity @ saturation [mm/s]
REAL, INTENT(IN) :: sathh(land_pts, sm_levels)
REAL, INTENT(IN) :: smvcst(land_pts, sm_levels)
REAL, INTENT(IN) :: smvcwt(land_pts, sm_levels)
REAL, INTENT(IN) :: smvccl(land_pts, sm_levels)

REAL, INTENT(IN) :: hcon(land_pts)      ! Soil thermal conductivity (W/m/K).
REAL, INTENT(IN) :: albsoil(land_pts)
REAL, INTENT(IN) :: reducedLAIdue2snow(mp)

TYPE(params_io_data_type), INTENT(IN)       :: pars
TYPE(met_type),            INTENT(OUT) :: met
TYPE(radiation_type),      INTENT(OUT) :: rad
TYPE(roughness_type),      INTENT(OUT) :: rough
TYPE(soil_snow_type),      INTENT(OUT) :: ssnow 
TYPE(balances_type),       INTENT(OUT) :: bal 
TYPE(canopy_type),         INTENT(OUT) :: canopy
TYPE(air_type),            INTENT(OUT) :: air
TYPE(bgc_pool_type),       INTENT(OUT) :: bgc
TYPE(sum_flux_type),       INTENT(OUT) :: sum_flux
TYPE(veg_parameter_type),  INTENT(OUT) :: veg        ! vegetation parameters
TYPE(soil_parameter_type), INTENT(OUT) :: soil      ! soil parameters

! "forcing"
REAL, INTENT(IN) :: lw_down(row_length,rows)
REAL, INTENT(IN) :: ls_rain(row_length,rows)
REAL, INTENT(IN) :: ls_snow(row_length,rows)
REAL, INTENT(IN) :: tl_1(row_length,rows)
REAL, INTENT(IN) :: qw_1(row_length,rows)
REAL, INTENT(IN) :: vshr_land(row_length,rows)
REAL, INTENT(IN) :: pstar(row_length,rows)
REAL, INTENT(IN) :: z1_tq(row_length,rows)
REAL, INTENT(IN) :: z1_uv(row_length,rows)

REAL, INTENT(IN) :: sw_down_VIS(row_length,rows)
REAL, INTENT(IN) :: sw_down_NIR(row_length,rows)
REAL, INTENT(IN) :: beamFrac_VIS(row_length,rows)
REAL, INTENT(IN) :: beamFrac_NIR(row_length,rows)
REAL, INTENT(IN) :: beamFrac_TOT(row_length,rows)
  
! prognostics   
REAL, INTENT(IN) :: canopy_tile(land_pts, nsurft)
REAL, INTENT(IN) :: snow_tile(land_pts, nsurft)
REAL, INTENT(IN) :: progs_soiltemp(land_pts, nsurft, sm_levels)
REAL, INTENT(IN) :: progs_soilmoisture(land_pts, nsurft, sm_levels)
REAL, INTENT(IN) :: progs_FrozenSoilFrac(land_pts, nsurft, sm_levels)
REAL, INTENT(IN) :: progs_SnowDepth(land_pts, nsurft,nsnl)
REAL, INTENT(IN) :: progs_SnowTemp(land_pts, nsurft,nsnl)
REAL, INTENT(IN) :: progs_SnowMass(land_pts, nsurft,nsnl)
REAL, INTENT(IN) :: progs_SnowDensity(land_pts, nsurft,nsnl)
REAL, INTENT(IN) :: progs_snowosurft(land_pts, nsurft)
REAL, INTENT(IN) :: progs_OneLyrSnowDensity(land_pts, nsurft)
REAL, INTENT(IN) :: progs_snowage(land_pts, nsurft)
INTEGER, INTENT(IN) :: progs_ThreeLayerSnowFlag(land_pts, nsurft)

INTEGER, INTENT(IN) :: SurfaceType(mp) !CABLE surface tile PFT/nveg
INTEGER, INTENT(IN) :: SoilType(mp)    !CABLE soil type per tile
REAL, INTENT(OUT)   :: HGT_pft_cbl(mp)
REAL, INTENT(OUT)   :: LAI_pft_cbl(mp)
    
!!jh:8/23 this was a problem BUT was never an issue because it is defined before it is used in soilsnow()    
!!snow_cond,     &
!GW progs:
!!  REAL,  DIMENSION(land_pts) :: & 
!!    slope_avg,&
!!    slope_std,&
!!    dz_gw,sy_gw,perm_gw,drain_gw
!!REAL ::  smgw_tile(land_pts,nsurft)
!!REAL,  DIMENSION(land_pts, nsurft) ::                         &
!!   !visc_sublayer_depth 
!GW progs: End
  
!CASA progs:
!!REAL, DIMENSION(land_pts,nsurft,10) :: &
!!  CPOOL_TILE,    & ! Carbon Pools
!!  NPOOL_TILE       ! Nitrogen Pools

!!REAL, DIMENSION(land_pts,nsurft,12) :: &
!!  PPOOL_TILE       ! Phosphorus Pools

!!REAL, DIMENSION(land_pts) :: &
!!  SOIL_ORDER,    & ! Soil Order (1 to 12)
!!  NIDEP,         & ! Nitrogen Deposition
!!  NIFIX,         & ! Nitrogen Fixation
!!  PWEA,          & ! Phosphorus from Weathering
!!  PDUST            ! Phosphorus from Dust

!!  GLAI,         &  ! Leaf Area Index for Prognostics LAI
!!  PHENPHASE,    &  ! Phenology Phase for Casa-CNP
REAL, INTENT(IN) :: npp_pft_acc(land_pts,npft)
REAL, INTENT(IN) :: resp_w_pft_acc(land_pts,npft) 
!CASA progs: End
  
!___ local vars
!jhan: this can be moved and USEd - needed  to pass arg
TYPE (climate_type) :: climate     ! climate variables
REAL  :: rho_water, rho_ice
  
LOGICAL, SAVE :: first_call = .TRUE.
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_explicit_driver"

rho_water = density_liq
rho_ice   = density_liq

IF (cable_user%GW_model) then
  rho_ice = density_ice
ENDIF      

IF(first_call) THEN  
  !--- fill CABLE fields from UM ancillaries/fields, CABLE prognostics
  CALL init_data( row_length, rows, land_pts, nsurft, npft, sm_levels,         &
                  nsnl, dzsoil, mp, nrb, CO2_MMR, tfrz, ICE_SurfaceType,       &
                  ICE_SoilType, land_index, surft_pts, surft_index, tile_frac, &
                  L_tile_pts, albsoil, bexp, hcon, satcon, sathh, smvcst,      &
                  smvcwt, smvccl, pars, tl_1, snow_tile, progs_soiltemp,       &
                  progs_soilmoisture, progs_FrozenSoilFrac,                    &
                  progs_OneLyrSnowDensity, progs_snowage,                      &
                  progs_ThreeLayerSnowFlag, progs_SnowDensity, progs_SnowDepth,&
                  progs_SnowTemp, progs_SnowMass, rad%trad, met%tk, veg, soil, &
                  canopy, ssnow, bgc, sum_flux, SurfaceType, SoilType,         &
                  npp_pft_acc,resp_w_pft_acc )

  !CALL init_data_sci( nsl, nsnl, soil%zse, mp, tfrz, ICE_SoilType, rad%trad,   &
  !                    met%tk, veg, soil, canopy, ssnow )
                      
  first_call = .FALSE.

ENDIF     

!--- update CABLE fields from UM forcings and equivalent fields at tis timestep 

CALL update_data( row_length, rows, land_pts, nsurft, npft, sm_levels,         &
                  nsnl, timestep, timestep_number, mp, nrb, CO2_MMR, canht_pft,&
                  lai_pft, land_index, surft_pts, surft_index, tile_frac,      &
                  L_tile_pts, cos_zenith_angle, latitude, longitude,           &
                  sw_down_VIS, sw_down_NIR, beamFrac_VIS, beamFrac_NIR,        &
                  beamFrac_TOT, lw_down, ls_rain, ls_snow, tl_1, qw_1,         &
                  vshr_land, pstar, z1_tq, z1_uv, canopy_tile, rad, met,       &
                  veg, soil, rough, canopy, ssnow, HGT_pft_cbl, LAI_pft_cbl,   &
                  reducedLAIdue2snow )

!CALL update_data_sci( mp, rad, met, veg, soil, canopy, ssnow, &
!                      canopy%vlaiw )
 
!---------------------------------------------------------------------!
!--- Feedback prognostic vcmax and daily LAI from casaCNP to CABLE ---!
!---------------------------------------------------------------------!
!block!IF(l_vcmaxFeedbk) call casa_feedback(ktau_gl,veg,casabiome,casapool,casamet)
!block!IF(l_laiFeedbk) veg%vlai(:) = casamet%glai(:)

!---------------------------------------------------------------------!
!--- cbm "mainly" controls the calling of model components         ---!  
!---------------------------------------------------------------------!
CALL cbm_expl( mp, nrb, timestep_number, timestep, air, bgc, canopy, met, bal, &
          rad, rough, soil, ssnow, sum_flux, veg, climate )

canopy%cansto =  canopy%oldcansto

RETURN
END SUBROUTINE cable_explicit_driver

END MODULE cable_explicit_driv_mod
