MODULE hydrol_mod_cbl
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HYDROL_MOD'

CONTAINS

SUBROUTINE hydrol_cbl (                                                            &
     land_pts, soil_pts, lice_pts, sm_levels, npft, nsurft, dim_cs1,           &
     asteps_since_triffid, timestep,                                           &
     l_inland, l_pdm, l_soil_sat_down, l_top, stf_sub_surf_roff,               &
     soil_index, lice_index,                                                   &
     surft_pts, surft_index, nsnow_surft,                                      &
     bexp_soilt, hcap_soilt, hcon_soilt,                                       &
     satcon_soilt, sathh_soilt,                                                &
     smvcst_soilt, smvcwt_soilt,                                               &
     catch_surft, frac_irr_soilt,  infil_surft,                                &
     non_lake_frac, slope_gb, tile_frac,                                       &
     con_rainfrac_land,  con_rain_land, ls_rainfrac_land, ls_rain_land,        &
     surf_ht_flux_ld,                                                          &
     ecan_surft, ext_irr_soilt,  ext_soilt,                                    &
     snowdepth_surft, melt_surft, snow_melt,                                   &
     snow_soil_htf,                                                            &
     a_fsat_soilt, c_fsat_soilt, a_fwet_soilt,  c_fwet_soilt,                  &
     fexp_soilt, ti_mean_soilt,                                                &
     npp_soilt, inlandout_atm_gb,                                              &
     canopy_surft, smcl_soilt, sthf_soilt,                                     &
     sthu_soilt,  sthu_irr_soilt, tsoil_deep_gb,                               &
     t_soil_soilt, t_soil_soilt_acc, tsurf_elev_surft,                         &
     fsat_soilt, fwetl_soilt, sthzw_soilt, zw_soilt,                           &
     cs_pool_soilt,resp_s_soilt,                                               &
     cs_ch4_soilt, fch4_wetl_acc_soilt,                                        &
     substr_ch4, mic_ch4, mic_act_ch4, acclim_ch4,                             &
     n_inorg_avail_pft, n_inorg_soilt_lyrs,                                    &
     n_leach_gb_acc,                                                           &
     canopy_gb, smc_soilt,                                                     &
     drain_soilt, dun_roff_soilt, sub_surf_roff_gb,                            &
     surf_roff_gb, tot_tfall_gb, tot_tfall_surft,                              &
     w_flux_soilt, qbase_soilt, qbase_l_soilt, qbase_zw_soilt,                 &
     fch4_wetl_soilt, fch4_wetl_cs_soilt,                                      &
     fch4_wetl_npp_soilt, fch4_wetl_resps_soilt,                               &
     n_leach_soilt,                                                            &
     progs_snow_surft, snow_mass_gb, work_cbl ) 

!Use in relevant subroutines
USE ancil_info,               ONLY: dim_cslayer, nsoilt
USE calc_baseflow_jules_mod,  ONLY: calc_baseflow_jules
USE calc_zw_inund_mod,        ONLY: calc_zw_inund
USE ch4_wetl_mod,             ONLY: ch4_wetl
USE elev_htc_mod,             ONLY: elev_htc
USE ereport_mod,              ONLY: ereport
USE ice_htc_mod,              ONLY: ice_htc
USE n_leach_mod,              ONLY: n_leach
USE soil_hyd_wt_mod,          ONLY: soil_hyd_wt
USE soil_hyd_update_mod,      ONLY: soil_hyd_update
USE soil_hyd_mod,             ONLY: soil_hyd
USE soil_htc_mod,             ONLY: soil_htc
USE soilmc_mod,               ONLY: soilmc
USE soilt_mod,                ONLY: soilt
USE surf_hyd_mod,             ONLY: surf_hyd
USE um_types,                 ONLY: real_jlslsm
USE water_constants_mod,      ONLY: rho_water  ! Density of pure water (kg/m3).

USE jules_hydrology_mod,      ONLY: l_wetland_unfrozen, ti_max, zw_max
USE jules_vegetation_mod,     ONLY: l_nitrogen
USE jules_irrig_mod,          ONLY: l_irrig_dmd
USE jules_surface_mod,        ONLY: l_elev_land_ice

USE jules_soil_biogeochem_mod, ONLY:                                           &
  ! imported scalar parameters
  soil_model_ecosse, soil_model_rothc, soil_model_1pool,                       &
  ! imported scalar variables
  soil_bgc_model, l_ch4_tlayered, l_layeredc, dim_ch4layer

USE jules_soil_mod, ONLY:                                                      &
  dzsoil,                                                                      &
     ! Thicknesses of the soil layers (m).
  dzsoil_elev, ns_deep


#if !defined(UM_JULES)
USE update_mod,     ONLY: l_imogen
USE model_time_mod, ONLY: timestep_len
#endif

!CABLE_LSM: Make avail CABLE hydrology interface subr to CALL
USE cable_hyd_main_mod, ONLY: cable_hyd_main
USE work_vars_mod_cbl,  ONLY: work_vars_type      ! and some kept thru timestep

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! The order of arguments and declarations should follow UMDP3.
! Within that, where possible variables should be grouped so that related
! variables (e.g. from the same area of science) are together.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Scalar arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Number of gridpoints.
  soil_pts,                                                                    &
    ! Number of soil points.
  lice_pts,                                                                    &
    ! Number of land ice points.
  sm_levels,                                                                   &
    ! Number of soil moisture levels.
  npft,                                                                        &
    ! Number of plant functional types
  nsurft,                                                                      &
    ! Number of tiles
  dim_cs1,                                                                     &
    ! Number of soil carbon pools
  asteps_since_triffid
    ! Number of atmospheric timesteps since last call to TRIFFID.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep
    ! Model timestep (s).

LOGICAL , INTENT(IN) ::                                                        &
  l_inland,                                                                    &
    ! True if re-routing inland basin flow to soil moisture.
  l_pdm,                                                                       &
    ! Flag for PDM hydrology.
  l_soil_sat_down,                                                             &
    ! Switch controlling direction of movement of
    ! soil moisture in excess of saturation.
  l_top,                                                                       &
    ! Flag for TOPMODEL-based hydrology.
  stf_sub_surf_roff
    ! Stash flag for sub-surface runoff.

TYPE(work_vars_type), INTENT(IN) :: work_cbl
!-----------------------------------------------------------------------------
! Array arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  soil_index(land_pts),                                                        &
    ! Array of soil points.
  lice_index(land_pts),                                                        &
    ! Array of land ice points.
  surft_pts(nsurft),                                                           &
    ! Number of tile points.
  surft_index(land_pts,nsurft),                                                &
    ! Index of tile points.
  nsnow_surft(land_pts,nsurft)
    ! Number of snow layers

! Soil characteristics.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  bexp_soilt(land_pts,nsoilt,sm_levels),                                       &
    ! Brooks & Corey exponent.
  hcap_soilt(land_pts,nsoilt,sm_levels),                                       &
    ! Soil heat capacity (J/K/m3).
  hcon_soilt(land_pts,nsoilt,0:sm_levels),                                     &
    ! Soil thermal conductivity (W/m/K).
  satcon_soilt(land_pts,nsoilt,0:sm_levels),                                   &
    ! Saturated hydraulic conductivity (kg/m2/s).
  sathh_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Saturated soil water pressure (m).
  smvcst_soilt(land_pts,nsoilt,sm_levels),                                     &
    ! Volumetric soil moisture concentration at saturation (m3 H2O/m3 soil).
  smvcwt_soilt(land_pts,nsoilt,sm_levels)
    ! Volumetric soil moisture concentration below which
    ! stomata close (m3 H2O/m3 soil).

! Surface characteristics.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  catch_surft(land_pts,nsurft),                                                &
    ! Canopy/surface capacity of land tiles (kg/m2).
  frac_irr_soilt(land_pts, nsoilt),                                            &
    ! Irrigation fraction for each soil tile.
  infil_surft(land_pts,nsurft),                                                &
    ! Maximum surface infiltration (kg m-2 s-1).
  non_lake_frac(land_pts),                                                     &
    ! Sum of fractions of surface tiles linked to soil (i.e. not using FLake)
  slope_gb(land_pts),                                                          &
    ! Terrain slope.
  tile_frac(land_pts,nsurft)
    ! Tile fractions.

! Precipitation variables and heat fluxes.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  con_rainfrac_land(land_pts),                                                 &
    ! Convective rain fraction
  con_rain_land(land_pts),                                                     &
    ! Convective rain (kg/m2/s).
  ls_rain_land(land_pts),                                                      &
    ! Large-scale rain (kg/m2/s).
  ls_rainfrac_land(land_pts),                                                  &
    ! large scale rain fraction
  surf_ht_flux_ld(land_pts)
    ! Net downward surface heat flux (W/m2).

! Evaporation variables.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ecan_surft(land_pts,nsurft),                                                 &
    ! Canopy evaporation from land tiles (kg/m2/s).
  ext_irr_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Extraction of water from each soil layer over irrigation (kg m-2 s-1).
  ext_soilt(land_pts,nsoilt,sm_levels)
    ! Extraction of water from each soil layer (kg/m2/s).

! Snow and melt
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  snowdepth_surft(land_pts,nsurft),                                            &
    ! Snow depth (on ground) (m).
  melt_surft(land_pts,nsurft),                                                 &
    ! Snowmelt on tiles (kg/m2/s).
  !snow_melt(land_pts),                                                         &
  !  ! Snowmelt (kg/m2/s).
  !  ! for CABLE this needs to be an INOUT var
  snow_soil_htf(land_pts,nsurft)
    ! Tiled snowpack-> soil heat flux.
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                       &
  snow_melt(land_pts)
    ! Snowmelt (kg/m2/s)

! TOPMODEL variables.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  a_fsat_soilt(land_pts,nsoilt),                                               &
    ! Fitting parameter for Fsat in LSH model.
  c_fsat_soilt(land_pts,nsoilt),                                               &
    ! Fitting parameter for Fsat in LSH model.
  a_fwet_soilt(land_pts,nsoilt),                                               &
    ! Fitting parameter for Fwet in LSH model.
  c_fwet_soilt(land_pts,nsoilt),                                               &
    ! Fitting parameter for Fwet in LSH model.
  fexp_soilt(land_pts,nsoilt),                                                 &
    ! Decay factor in Sat. Conductivity in deep LSH/TOPMODEL layer.
  ti_mean_soilt(land_pts,nsoilt)
    ! Mean topographic index.

! Variables related to methane calculation.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  npp_soilt(land_pts,nsoilt)
    ! Gridbox mean net primary productivity (kg C/m2/s).

! Other variables.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  inlandout_atm_gb(land_pts)
    ! IN TRIP INLAND BASIN OUTFLOW FOR LAND POINTS ONLY,kg/m2/s=mm.
!-----------------------------------------------------------------------------
! Array arguments with INTENT(IN OUT):
!-----------------------------------------------------------------------------

! Canopy water.
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  canopy_surft(land_pts,nsurft)
    ! Canopy water content for land tiles (kg/m2).

! Soil water and temperature.
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  smcl_soilt(land_pts,nsoilt,sm_levels),                                       &
    ! Soil moisture content of each
!                          !       layer (kg/m2).
  sthf_soilt(land_pts,nsoilt,sm_levels),                                       &
    ! Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
  sthu_soilt(land_pts,nsoilt,sm_levels),                                       &
    ! Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
  sthu_irr_soilt(land_pts, nsoilt, sm_levels),                                 &
    !  Unfrozen soil wetness over irrigation.
  tsoil_deep_gb(land_pts,ns_deep),                                             &
    ! Deep soil temperature (K).
  t_soil_soilt(land_pts,nsoilt,sm_levels),                                     &
    ! Sub-surface temperatures (K).
  t_soil_soilt_acc(land_pts,nsoilt,sm_levels),                                 &
    ! Sub-surface temperature on layers and soil tiles
    ! accumulated over TRIFFID timestep (K).
  tsurf_elev_surft(land_pts,nsurft)
    ! Tiled sub-surface temperatures (K).

! TOPMODEL variables.
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  fsat_soilt(land_pts,nsoilt),                                                 &
     ! Surface saturation fraction.
  fwetl_soilt(land_pts,nsoilt),                                                &
     ! Wetland fraction.
  sthzw_soilt(land_pts,nsoilt),                                                &
     ! Soil moisture fraction in deep LSH/TOPMODEL layer.
  zw_soilt(land_pts,nsoilt)
     ! Water table depth (m).

! Soil carbon variables.
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1),                          &
    ! Soil carbon (kg C/m2).
!                          !   For RothC (dim_cs1=4), the pools
!                          !   are DPM, RPM, biomass and humus.
  resp_s_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)
    ! Soil respiration in pools (kg C/m2/s).

! Methane variables.
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  cs_ch4_soilt(land_pts,nsoilt),                                               &
    ! Soil carbon used in CH4 wetlands if TRIFFID is switched off (kg C/m2).
  fch4_wetl_acc_soilt(land_pts,nsoilt),                                        &
    ! Accumulated scaled wetland methane flux (kg C/m2).
  substr_ch4(land_pts,dim_ch4layer),                                           &
    ! Dissolved substrate that methaogens consume (kg C/m2)
  mic_ch4(land_pts,dim_ch4layer),                                              &
    ! Methanogenic biomass (kg C/m2)
  mic_act_ch4(land_pts,dim_ch4layer),                                          &
    ! Activity level of methanogenic biomass (fraction)
  acclim_ch4(land_pts,dim_ch4layer)
    ! Acclimation factor for microbial trait adaptation

! Nitrogen variables.
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  n_inorg_avail_pft(land_pts,npft,dim_cslayer),                                &
    ! Available inorganic N for PFTs (kg N m-2).
  n_inorg_soilt_lyrs(land_pts,nsoilt,dim_cslayer),                             &
    ! Inorganic N pool on soil levels (kg N m-2).
  n_leach_gb_acc(land_pts)
    ! Accumulated leached nitrogen diagnostic on TRIFFID timesteps
    ! (kg m-2 (360days)-1).

!-----------------------------------------------------------------------------
! Array arguments with INTENT(OUT):
!-----------------------------------------------------------------------------
! Canopy water and soil moisture.
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  canopy_gb(land_pts),                                                         &
    ! Gridbox canopy water content (kg/m2).
  smc_soilt(land_pts,nsoilt)
    ! Available soil moisture in a layer at the surface (kg/m2)

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                       &
  tot_tfall_surft(land_pts,nsurft),                                            &
    ! Surface-tiled contributions to throughfall.
  dun_roff_soilt(land_pts,nsoilt)
    ! Dunne part of sfc runoff (kg/m2/s).

! Water fluxes.
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  drain_soilt(land_pts,nsoilt),                                                &
    ! Drainage out of sm_levels'th level (kg/m2/s).
  sub_surf_roff_gb(land_pts),                                                  &
    ! Sub-surface runoff (kg/m2/s).
  surf_roff_gb(land_pts),                                                      &
    ! Surface runoff (kg/m2/s).
  tot_tfall_gb(land_pts),                                                      &
    ! Total throughfall (kg/m2/s).
  w_flux_soilt(land_pts,nsoilt,0:sm_levels)
    ! Fluxes of water between layers (kg/m2/s).

! TOPMODEL variables.
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  qbase_soilt(land_pts,nsoilt),                                                &
    ! Base flow (kg/m2/s).
  qbase_l_soilt(land_pts,nsoilt,sm_levels+1),                                  &
    ! Base flow from each level (kg/m2/s).
  qbase_zw_soilt(land_pts,nsoilt)
    ! Base flow from deep LSH/TOPMODEL layer (kg/m2/s).

! Methane and nitrigen variables.
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  fch4_wetl_soilt(land_pts,nsoilt),                                            &
    ! Scaled wetland methane flux (default substrate) for use in
    ! atmos chemistry model (10^-9 kg C/m2/s).
  fch4_wetl_cs_soilt(land_pts,nsoilt),                                         &
    ! Scaled methane flux (soil carbon substrate) (kg C/m2/s).
  fch4_wetl_npp_soilt(land_pts,nsoilt),                                        &
    ! Scaled methane flux (npp substrate) (kg C/m2/s).
  fch4_wetl_resps_soilt(land_pts,nsoilt),                                      &
    ! Scaled methane flux (soil respiration substrate) (kg C/m2/s).
  n_leach_soilt(land_pts,nsoilt)
    ! Leached N (kg m-2 s-1).

!CABLE
REAL, INTENT(OUT) :: progs_snow_surft(land_pts,nsurft)
REAL, INTENT(OUT) :: snow_mass_gb(land_pts)        ! OUT Gridbox snowmass (kg/m2)     
!-----------------------------------------------------------------------------
! Local parameters:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER :: to_kg_conversion = 1.0e-9
    ! multiplier for converting to kgC for wetland CH4 and IMOGEN

!-----------------------------------------------------------------------------
! Local scalar variables:
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  i, j,                                                                        &
  n,                                                                           &
    ! Counter for soil level.
  m,                                                                           &
    ! Counter for soil tile.
  errorstatus

!-----------------------------------------------------------------------------
! Local array variables:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  dsmc_dt_soilt(land_pts,nsoilt),                                              &
    ! Rate of change of soil moisture due to water falling onto the
    ! surface after surface runoff (kg/m2/s).
  ksz_soilt(land_pts,nsoilt,0:sm_levels),                                      &
    ! Saturated hydraulic conductivity in layer (kg/m2/s).
  qbase_unfr_soilt(land_pts,nsoilt),                                           &
    ! Base flow in unfrozen soil (kg/m2/s).
  qbase_l_unfr_soilt(land_pts,nsoilt,sm_levels+1),                             &
    ! As qbase_l but for unfrozen soil (kg/m2/s).
  top_crit_soilt(land_pts,nsoilt),                                             &
    ! Critical TI when ZW <=0.0
  dumtop_crit_soilt(land_pts,nsoilt),                                          &
    ! Dummy for top_crit_soilt
  dumsthf_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Dummy Frozen soil moisture content of each layer as a fraction of
    ! saturation (always set to 0).
  zdepth(0:sm_levels),                                                         &
    ! Lower soil layer boundary depth (m).
  tsoil_d_soilt(land_pts,nsoilt),                                              &
    ! Soil temperature in the top metre
  zw_inund_soilt(land_pts,nsoilt),                                             &
    ! Water table depth used
  wutot_soilt(land_pts,nsoilt),                                                &
    ! Ratio of unfrozen to total soil moisture at ZW.
  surf_roff_inc_soilt(land_pts,nsoilt),                                        &
    ! Increment to tiled surface runoff (kg m-2 s-1).
  surf_roff_soilt(land_pts,nsoilt),                                            &
    ! Soil-tiled contributions to surface runoff (kg m-2 s-1).
  sub_surf_roff_soilt(land_pts,nsoilt),                                        &
    ! Soil-tiled contributions to subsurface runoff (kg m-2 s-1).
    smcl_old_soilt(land_pts,nsoilt,sm_levels)
    ! Retain oriignal soil moisture for leaching code

! Variables required for irrigation code
REAL(KIND=real_jlslsm) ::                                                      &
  w_flux_irr_soilt(land_pts,nsoilt,0:sm_levels),                               &
    ! The fluxes of water between layers in irrigated fraction (kg/m2/s).
  w_flux_nir_soilt(land_pts,nsoilt,0:sm_levels),                               &
    ! The fluxes of water between layers in non-irrigated fraction (kg/m2/s).
  smcl_irr_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! Total soil moisture contents of each layer in irrigated
    ! fraction (kg/m2).
  smcl_nir_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! Total soil moisture contents of each layer in non-irrigated
    ! fraction (kg/m2).
  sthu_nir_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation in irrigated fraction.
  ext_nir_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Extraction of water from each soil layer in non-irrigated fraction
    ! (kg/m2/s).
  smclsat_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! The saturation moisture content of each layer (kg/m2).
  smclzw_soilt(land_pts,nsoilt),                                               &
    ! moisture content in deep layer(kg/m2).
  smclsatzw_soilt(land_pts,nsoilt)
    ! moisture content in deep layer (kg/m2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HYDROL_cable'

! End of header--------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if !defined(UM_JULES)
 timestep =  REAL(timestep_len)
#endif

canopy_gb = canopy_gb
!CM2.1!!CABLE_LSM: call CABLE hydrology
call cable_hyd_main( land_pts, nsurft, tile_frac, timestep,                    &
                     snow_mass_gb, progs_snow_surft, snow_melt,                &
                     surf_roff_gb, sub_surf_roff_gb,                           &
                     tot_tfall_gb, work_cbl%snow_surft, melt_surft,            &
                     work_cbl%lying_snow, work_cbl%surf_roff,                  &
                     work_cbl%sub_surf_roff, work_cbl%tot_tfall )
                    

! Calculate soil carbon for use in the wetland CH4 scheme only
! (only used if single-pool C model is used):
IF ( soil_bgc_model == soil_model_1pool ) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,m,n)              &
!$OMP SHARED(soil_pts, nsoilt, soil_index, cs_ch4_soilt, dim_cslayer,          &
!$OMP        cs_pool_soilt)
  DO m = 1, nsoilt
    DO j = 1, soil_pts
      i = soil_index(j)
      cs_ch4_soilt(i,m) = 0.0
      DO n = 1,dim_cslayer
        cs_ch4_soilt(i,m) = cs_ch4_soilt(i,m) + cs_pool_soilt(i,m,n,1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,m,n)                                    &
!$OMP SHARED(sm_levels, nsoilt, smcl_old_soilt, smcl_soilt, w_flux_soilt,      &
!$OMP        w_flux_irr_soilt, land_pts)
! smcl_old_soilt calculated for leaching code
DO n = 1, sm_levels
  DO m = 1, nsoilt
!$OMP DO SCHEDULE(STATIC)
    DO i = 1, land_pts
      smcl_old_soilt(i,m,n) = smcl_soilt(i,m,n)
    END DO
!$OMP END DO NOWAIT
  END DO
END DO

! Initialise w_flux variables that are used in irrigation code
DO n = 0, sm_levels
  DO m = 1, nsoilt
!$OMP DO SCHEDULE(STATIC)
    DO i = 1, land_pts
      w_flux_soilt(i,m,n)     = 0.0
      w_flux_irr_soilt(i,m,n) = 0.0
      ! to prevent random values reported over areas
      ! that are not included as soil points (i.e. ice points)
    END DO
!$OMP END DO NOWAIT
  END DO
END DO
!$OMP END PARALLEL
!CABLE_LSM: call CABLE hydrology
!!call cable_hyd_main( npnts, nsurft, lying_snow, snow_surft, surf_roff,                     &
 !!                    sub_surf_roff, tot_tfall )

!CM2 vn starts aboout here
!-----------------------------------------------------------------------------
! Set up variables required for LSH scheme:
!-----------------------------------------------------------------------------
! set level zero depth to 0
zdepth(0) = 0.0

DO n = 1, sm_levels
  zdepth(n) = zdepth(n-1) + dzsoil(n)
END DO
!jhan:fudge
dsmc_dt_soilt(:,:) = 0.0
!-----------------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------------
!!CABLE_LSM: omit for CABLE: BUT initialize dsmc_dt for fpe0, otherwise?  
!!CALL surf_hyd (land_pts, nsurft, sm_levels, soil_pts, timestep, l_top, l_pdm,  &
!!               surft_pts, surft_index, soil_index,                             &
!!               catch_surft, ecan_surft, tile_frac, non_lake_frac, infil_surft, &
!!               con_rain_land, ls_rain_land, con_rainfrac_land,                 &
!!               ls_rainfrac_land, melt_surft, slope_gb, snow_melt,              &
!!               fsat_soilt, smvcst_soilt,                                       &
!!               sthf_soilt, sthu_soilt,                                         &
!!               canopy_surft, canopy_gb, dsmc_dt_soilt,                         &
!!               surf_roff_gb, tot_tfall_gb, tot_tfall_surft,                    &
!!               dun_roff_soilt, surf_roff_soilt)
!
!-----------------------------------------------------------------------------
! Specify the reduction of hydraulic conductivity with depth.
! Initialiase base flow to zero.
!-----------------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, m, n)                               &
!$OMP SHARED(sm_levels, nsoilt, soil_pts, soil_index, ksz_soilt, satcon_soilt, &
!$OMP        smclsat_soilt, qbase_l_soilt, qbase_l_unfr_soilt, dumsthf_soilt,  &
!$OMP        dzsoil, smvcst_soilt, land_pts, qbase_soilt,                      &
!$OMP        qbase_zw_soilt, wutot_soilt, drain_soilt, qbase_unfr_soilt,       &
!$OMP        zw_inund_soilt, surf_roff_inc_soilt)
DO n = 0, sm_levels
  DO m = 1, nsoilt
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, soil_pts
      i = soil_index(j)
      ksz_soilt(i,m,n) = satcon_soilt(i,m,n)
    END DO
!$OMP END DO NOWAIT
  END DO
END DO

DO n = 1, sm_levels
  DO m = 1, nsoilt
!$OMP DO SCHEDULE(STATIC)
    DO j = 1,soil_pts
      i = soil_index(j)
      smclsat_soilt(i,m,n)      = rho_water * dzsoil(n) * smvcst_soilt(i,m,n)
      qbase_l_soilt(i,m,n)      = 0.0
      qbase_l_unfr_soilt(i,m,n) = 0.0
      dumsthf_soilt(i,m,n)      = 0.0
    END DO
!$OMP END DO NOWAIT
  END DO
END DO

DO m = 1, nsoilt
!$OMP DO SCHEDULE(STATIC)
  DO i = 1,land_pts
    qbase_soilt(i,m)      = 0.0
    qbase_zw_soilt(i,m)   = 0.0
    wutot_soilt(i,m)      = 0.0
    drain_soilt(i,m)      = 0.0
    qbase_unfr_soilt(i,m) = 0.0
    zw_inund_soilt(i,m)   = 0.0
    ! Initialise runoff increment.
    surf_roff_inc_soilt(i,m) = 0.0
  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL

!!CABLE_LSM: omit for CABLE
!IF (l_top) THEN
!  IF (soil_pts /= 0) THEN
!    DO m = 1, nsoilt
!      CALL calc_baseflow_jules(                                                &
!        land_pts, sm_levels, soil_pts, soil_index,                             &
!        bexp_soilt(:,m,:), fexp_soilt(:,m), sthf_soilt(:,m,:),                 &
!        ti_mean_soilt(:,m), zdepth, zw_soilt(:,m), ksz_soilt(:,m,:),           &
!        qbase_soilt(:,m), qbase_l_soilt(:,m,:), dumtop_crit_soilt(:,m) )
!    END DO
!
!    IF (l_wetland_unfrozen) THEN
!      DO m = 1,nsoilt
!        CALL calc_zw_inund(land_pts, sm_levels, soil_pts, soil_index, zdepth,  &
!          bexp_soilt(:,m,:), sathh_soilt(:,m,:), smclsat_soilt(:,m,:),         &
!          smcl_soilt(:,m,:), sthu_soilt(:,m,:), sthzw_soilt(:,m),              &
!          zw_soilt(:,m), zw_inund_soilt(:,m), wutot_soilt(:,m))
!
!        ! Now call again to get the unfrozen equivalents to calculate fsat and
!        ! fwet:
!        CALL calc_baseflow_jules(                                              &
!          land_pts, sm_levels, soil_pts, soil_index,                           &
!          bexp_soilt(:,m,:), fexp_soilt(:,m), dumsthf_soilt(:,m,:),            &
!          ti_mean_soilt(:,m), zdepth, zw_inund_soilt(:,m), ksz_soilt(:,m,:),   &
!          qbase_unfr_soilt(:,m), qbase_l_unfr_soilt(:,m,:),                    &
!          top_crit_soilt(:,m) )
!      END DO
!    ELSE
!      DO m = 1, nsoilt
!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i)                    &
!!$OMP SHARED(m, land_pts, top_crit_soilt, dumtop_crit_soilt)
!        DO i = 1,land_pts
!          top_crit_soilt(i,m) = dumtop_crit_soilt(i,m)
!        END DO
!!$OMP END PARALLEL DO
!      END DO
!    END IF
!
!  END IF
!END IF  !  l_top
!!CABLE_LSM:End

IF (l_inland) THEN
  DO i = 1,land_pts

    ! Add inland basin outflow to change in soil moisture store.
    ! Note for soil tiling- this is only used by the riv_intctl_1a, which is
    ! not compatible with nsoilt > 1.
    dsmc_dt_soilt(i,1) = dsmc_dt_soilt(i,1) + inlandout_atm_gb(i)
  END DO
END IF

!-----------------------------------------------------------------------------
! Update the layer soil moisture contents and calculate the
! gravitational drainage.
!-----------------------------------------------------------------------------
!!!CABLE_LSM: omit for CABLE 
!!!CM2IF (soil_pts /= 0) THEN
!!
!!  !CM2IF (l_irrig_dmd) THEN
!!
!!  !CM2  !-------------------------------------------------------------------------
!!  !CM2  ! If l_irrig_dmd = TRUE, call soil_hyd separately for irrigated and
!!  !CM2  ! non-irrigated fraction
!!  !CM2  ! afterwards, call soil_hyd ONLY to update water table/drainage with
!!  !CM2  ! gridbox total w_flux_soilt, smcl_soilt
!!  !CM2  !-------------------------------------------------------------------------
!!
!!  !CM2  !--------------------------------------------------------------------------
!!  !CM2  ! Split tile values into values for irrigated and non-irrigated fractions.
!!  !CM2  !--------------------------------------------------------------------------
!!  !CM2  DO n = 1,sm_levels
!!  !CM2    DO m = 1, nsoilt
!!  !CM2      DO j = 1,soil_pts
!!  !CM2        i = soil_index(j)
!!
!!  !CM2        ! hadrd - gridbox sthu_soilt is assumed to be combination of
!!  !CM2        ! sthu_soilt of non-irrigated fraction and sthu_irr_soilt, i.e.
!!  !CM2        ! sthu_soilt = frac_irr_soilt*sthu_irr_soilt + (1-frac_irr_soilt)
!!  !CM2        !              *sthu_nir_soilt
!!  !CM2        IF ( frac_irr_soilt(i,m) < 1.0 ) THEN
!!  !CM2          sthu_nir_soilt(i,m,n) = (sthu_soilt(i,m,n) - frac_irr_soilt(i,m)   &
!!  !CM2                                  * sthu_irr_soilt(i,m,n))                   &
!!  !CM2                                  / (1.0 - frac_irr_soilt(i,m))
!!  !CM2          ext_nir_soilt(i,m,n)  = (ext_soilt(i,m,n)  - frac_irr_soilt(i,m)   &
!!  !CM2                                  * ext_irr_soilt(i,m,n))                    &
!!  !CM2                                  / (1.0 - frac_irr_soilt(i,m))
!!  !CM2        ELSE
!!  !CM2          sthu_nir_soilt(i,m,n) = sthu_soilt(i,m,n)
!!  !CM2          ext_nir_soilt(i,m,n)  = ext_soilt(i,m,n)
!!  !CM2        END IF
!!
!!  !CM2        smcl_irr_soilt(i,m,n) = smcl_soilt(i,m,n)                            &
!!  !CM2                                + (sthu_irr_soilt(i,m,n)                     &
!!  !CM2                                   - sthu_soilt(i,m,n))                      &
!!  !CM2                                * smclsat_soilt(i,m,n)
!!  !CM2        smcl_nir_soilt(i,m,n) = smcl_soilt(i,m,n)                            &
!!  !CM2                                + (sthu_nir_soilt(i,m,n)                     &
!!  !CM2                                   - sthu_soilt(i,m,n))                      &
!!  !CM2                                * smclsat_soilt(i,m,n)
!!  !CM2      END DO
!!  !CM2    END DO
!!  !CM2  END DO
!!
!!  !CM2  DO m = 1, nsoilt
!!  !CM2    ! First call soil_hyd for non-irrigated fraction.
!!  !CM2    ! Note that the values of smclsat_soilt, smclzw_soilt and smclsatzw_soilt
!!  !CM2    ! calculated here and in the next call are later replaced by the
!!  !CM2    ! recalculated values from soil_hyd_update.
!!  !CM2    CALL soil_hyd (                                                          &
!!  !CM2      land_pts, sm_levels, soil_pts, timestep, l_top, l_soil_sat_down,       &
!!  !CM2      soil_index, bexp_soilt(:,m,:), dzsoil,                                 &
!!  !CM2      ext_nir_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),            &
!!  !CM2      sathh_soilt(:,m,:), sthzw_soilt(:,m), smvcst_soilt(:,m,:),             &
!!  !CM2      qbase_l_soilt(:,m,:), zdepth,                                          &
!!  !CM2      smcl_nir_soilt(:,m,:), sthu_nir_soilt(:,m,:), smclsat_soilt(:,m,:),    &
!!  !CM2      w_flux_nir_soilt(:,m,:), smclzw_soilt(:,m), smclsatzw_soilt(:,m))
!!
!!  !CM2    ! Next call soil_hyd for irrigated fraction.
!!  !CM2    CALL soil_hyd (                                                          &
!!  !CM2      land_pts, sm_levels, soil_pts, timestep, l_top, l_soil_sat_down,       &
!!  !CM2      soil_index, bexp_soilt(:,m,:), dzsoil,                                 &
!!  !CM2      ext_irr_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),            &
!!  !CM2      sathh_soilt(:,m,:), sthzw_soilt(:,m), smvcst_soilt(:,m,:),             &
!!  !CM2      qbase_l_soilt(:,m,:), zdepth,                                          &
!!  !CM2      smcl_irr_soilt(:,m,:), sthu_irr_soilt(:,m,:), smclsat_soilt(:,m,:),    &
!!  !CM2      w_flux_irr_soilt(:,m,:), smclzw_soilt(:,m), smclsatzw_soilt(:,m))
!!  !CM2  END DO
!!
!!  !CM2  !--------------------------------------------------------------------------
!!  !CM2  ! Calculate values for the whole soil tile by combining irrigated and
!!  !CM2  ! non-irrigated values.
!!  !CM2  !--------------------------------------------------------------------------
!!  !CM2  DO m = 1,nsoilt
!!  !CM2    DO n = 0,sm_levels
!!  !CM2      DO j = 1,soil_pts
!!  !CM2        i = soil_index(j)
!!
!!  !CM2        ! Ensure sensible values if irrigation fraction is very small by
!!  !CM2        ! using values from the non-irrigated fraction.
!!  !CM2        IF ( frac_irr_soilt(i,m) <= EPSILON(1.0) ) THEN
!!  !CM2          w_flux_soilt(i,m,n)     = w_flux_nir_soilt(i,m,n)
!!  !CM2          IF ( n >= 1 ) THEN
!!  !CM2            sthu_irr_soilt(i,m,n) = sthu_nir_soilt(i,m,n)
!!  !CM2            smcl_soilt(i,m,n)     = smcl_nir_soilt(i,m,n)
!!  !CM2            sthu_soilt(i,m,n)     = sthu_nir_soilt(i,m,n)
!!  !CM2          END IF
!!  !CM2        ELSE
!!  !CM2          w_flux_soilt(i,m,n) = frac_irr_soilt(i,m)                          &
!!  !CM2                                * w_flux_irr_soilt(i,m,n)                    &
!!  !CM2                                + ( 1.0 - frac_irr_soilt(i,m) )              &
!!  !CM2                                * w_flux_nir_soilt(i,m,n)
!!  !CM2          IF ( n >= 1 ) THEN
!!  !CM2            smcl_soilt(i,m,n) = frac_irr_soilt(i,m) * smcl_irr_soilt(i,m,n)  &
!!  !CM2                                + ( 1.0 - frac_irr_soilt(i,m) )              &
!!  !CM2                                * smcl_nir_soilt(i,m,n)
!!  !CM2            sthu_soilt(i,m,n) = frac_irr_soilt(i,m) * sthu_irr_soilt(i,m,n)  &
!!  !CM2                                + ( 1.0 - frac_irr_soilt(i,m) )              &
!!  !CM2                                * sthu_nir_soilt(i,m,n)
!!  !CM2          END IF
!!  !CM2        END IF
!!  !CM2      END DO  !  soil points
!!  !CM2    END DO  !  layers
!!  !CM2  END DO  !  tiles
!!
!!  !CM2  !--------------------------------------------------------------------------
!!  !CM2  ! Recalculate values for the whole soil tile that were over-written by the
!!  !CM2  ! separate calls to soil_hyd above.
!!  !CM2  !--------------------------------------------------------------------------
!!  !CM2  DO m = 1, nsoilt
!!  !CM2    CALL soil_hyd_update(land_pts, sm_levels, soil_pts, soil_index, dzsoil,  &
!!  !CM2                         smvcst_soilt(:,m,:), zdepth, smclzw_soilt(:,m),     &
!!  !CM2                         sthzw_soilt(:,m), smclsat_soilt(:,m,:),             &
!!  !CM2                         smclsatzw_soilt(:,m))
!!  !CM2  END DO
!!
!!  !CM2ELSE
!!    ! .NOT. l_irrig_dmd
!!
!!    !CM2DO m = 1, nsoilt
!!    !CM2  CALL soil_hyd (                                                          &
!!    !CM2    land_pts, sm_levels, soil_pts, timestep, l_top, l_soil_sat_down,       &
!!    !CM2    soil_index, bexp_soilt(:,m,:), dzsoil,                                 &
!!    !CM2    ext_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),                &
!!    !CM2    sathh_soilt(:,m,:), sthzw_soilt(:,m), smvcst_soilt(:,m,:),             &
!!    !CM2    qbase_l_soilt(:,m,:), zdepth,                                          &
!!    !CM2    smcl_soilt(:,m,:), sthu_soilt(:,m,:), smclsat_soilt(:,m,:),            &
!!    !CM2    w_flux_soilt(:,m,:), smclzw_soilt(:,m), smclsatzw_soilt(:,m))
!!    !CM2END DO
!!
!!  !CM2END IF  !  l_irrig_dmd

  !----------------------------------------------------------------------------
  ! Update further stores and fluxes.
  !----------------------------------------------------------------------------
  !!IF (nsoilt == 1) THEN
  !!  !To maintain bit-comparability, we need to call with the _gb version of
  !!  !surface runoff.
  !!  m = 1
  !!  CALL soil_hyd_wt (                                                         &
  !!    land_pts, sm_levels, soil_pts, timestep, stf_sub_surf_roff,              &
  !!    l_top,  soil_index,                                                      &
  !!    bexp_soilt(:,m,:), dsmc_dt_soilt(:,m), sathh_soilt(:,m,:),               &
  !!    smclsat_soilt(:,m,:), smclsatzw_soilt(:,m), smvcst_soilt(:,m,:),         &
  !!    w_flux_soilt(:,m,:), smcl_soilt(:,m,:),                                  &
  !!    surf_roff_gb,                                                            &
  !!    qbase_l_soilt(:,m,:), smclzw_soilt(:,m), zw_soilt(:,m),                  &
  !!    sub_surf_roff_soilt(:,m), drain_soilt(:,m), qbase_soilt(:,m),            &
  !!    sthzw_soilt(:,m), surf_roff_inc_soilt(:,m) )

  !!  ! For a single soil tile, simply copy across to the output variable.
  !!  sub_surf_roff_gb(:) = sub_surf_roff_soilt(:,m)
  !!  ! Update the tiled surface runoff.
  !!  surf_roff_soilt(:,m) = surf_roff_soilt(:,m) + surf_roff_inc_soilt(:,m)
  !!ELSE

  !!  ! Initialise output variable
  !!  DO i = 1, land_pts
  !!    sub_surf_roff_gb(i) = 0.0
  !!  END DO
  !!  DO m = 1, nsoilt
  !!    CALL soil_hyd_wt (                                                       &
  !!      land_pts, sm_levels, soil_pts, timestep, stf_sub_surf_roff,            &
  !!      l_top, soil_index,                                                     &
  !!      bexp_soilt(:,m,:), dsmc_dt_soilt(:,m), sathh_soilt(:,m,:),             &
  !!      smclsat_soilt(:,m,:), smclsatzw_soilt(:,m), smvcst_soilt(:,m,:),       &
  !!      w_flux_soilt(:,m,:), smcl_soilt(:,m,:),                                &
  !!      surf_roff_soilt(:,m),                                                  &
  !!      qbase_l_soilt(:,m,:), smclzw_soilt(:,m), zw_soilt(:,m),                &
  !!      sub_surf_roff_soilt(:,m), drain_soilt(:,m), qbase_soilt(:,m),          &
  !!      sthzw_soilt(:,m), surf_roff_inc_soilt(:,m) )

  !!    ! For multiple soil tiles, add up the contributions, allowing for frac
  !!    sub_surf_roff_gb(:) = sub_surf_roff_gb(:)                                &
  !!                          + ( tile_frac(:,m) * sub_surf_roff_soilt(:,m))
  !!    surf_roff_gb(:) = surf_roff_gb(:) + tile_frac(:,m)                       &
  !!                                        * surf_roff_inc_soilt(:,m)

  !!  END DO
  !!END IF !nsoilt == 1

  !---------------------------------------------------------------------------
  ! Calculate surface saturation and wetland fractions:
  !---------------------------------------------------------------------------
  !CM2IF (l_top) THEN
  !CM2  DO m = 1, nsoilt

  !CM2    DO i = 1,land_pts
  !CM2      fsat_soilt(i,m)  = 0.0
  !CM2      fwetl_soilt(i,m) = 0.0

  !CM2      ! Zero soil porosity over land ice:
  !CM2      IF (smvcst_soilt(i,m,sm_levels) <= 0.0) THEN
  !CM2        zw_soilt(i,m) = zw_max
  !CM2      END IF
  !CM2    END DO

  !CM2    DO j = 1,soil_pts
  !CM2      i = soil_index(j)
  !CM2      qbase_zw_soilt(i,m) = qbase_l_soilt(i,m,sm_levels+1)

  !CM2      !Now use fit for fsat_soilt and fwet:
  !CM2      IF (l_wetland_unfrozen) THEN
  !CM2        fsat_soilt(i,m)  = wutot_soilt(i,m) * a_fsat_soilt(i,m)              &
  !CM2                           * EXP(-c_fsat_soilt(i,m) * top_crit_soilt(i,m))
  !CM2        fwetl_soilt(i,m) = wutot_soilt(i,m) * a_fwet_soilt(i,m)              &
  !CM2                           * EXP(-c_fwet_soilt(i,m) * top_crit_soilt(i,m))
  !CM2      ELSE
  !CM2        fsat_soilt(i,m)  = a_fsat_soilt(i,m)                                 &
  !CM2                           * EXP(-c_fsat_soilt(i,m) * top_crit_soilt(i,m))
  !CM2        fwetl_soilt(i,m) = a_fwet_soilt(i,m)                                 &
  !CM2                           * EXP(-c_fwet_soilt(i,m) * top_crit_soilt(i,m))
  !CM2      END IF

  !CM2      IF (top_crit_soilt(i,m) >= ti_max) THEN
  !CM2        fsat_soilt(i,m)  = 0.0
  !CM2        fwetl_soilt(i,m) = 0.0
  !CM2      END IF

  !CM2    END DO
  !CM2  END DO
  !CM2END IF  !  l_top

!CM2ELSE ! soil pts

  !---------------------------------------------------------------------------
  ! If required by STASH flag and there are no soil points,
  ! set sub-surface runoff to zero.
  !---------------------------------------------------------------------------
  !CM2IF (stf_sub_surf_roff) THEN
  !CM2  DO i = 1,land_pts
  !CM2    sub_surf_roff_gb(i) = 0.0
  !CM2  END DO
  !CM2END IF

!CM2END IF  !  soil_pts

!-----------------------------------------------------------------------------
! Update the soil temperatures and the frozen moisture fractions
!-----------------------------------------------------------------------------

!=============================================================================
! *NOTICE REGARDING SOIL TILING**
!
!The following section facilitates the use of soil tiling. As implemented,
!there are two soil tiling options:
!
!nsoilt == 1
!Operate as with a single soil tile, functionally identical to JULES upto
! at least vn4.7 (Oct 2016)
! This means that a soilt variable being passed 'up' to the surface is
! broadcast to the surft variable (with weighting by frac if requred)
!
!nsoilt > 1
!Operate with nsoilt = nsurft, with a direct mapping between them
! This means that a soilt variable being passed 'up' to the surface is simply
! copied into the surft variable
!
! This will need to be refactored for other tiling approaches. This note
! will be replicated elsewhere in the code as required
!
!These comments apply until **END NOTICE REGARDING SOIL TILING**
!=============================================================================
!CM2IF (soil_pts /= 0) THEN
!CM2  ! When using soil tiling, we can use the surface tiled version of
!CM2  ! surf_ht_flux_ld, snow_soil_htf. The _ld version is a gridbox mean
!CM2  ! calculated in snow and passed through.
!CM2  IF (nsoilt == 1) THEN
!CM2    m = 1
!CM2    CALL soil_htc (                                                            &
!CM2      land_pts, sm_levels, nsurft, soil_pts, timestep, soil_index,             &
!CM2      surft_pts, surft_index, nsnow_surft,                                     &
!CM2      bexp_soilt(:,m,:), dzsoil, tile_frac,                                    &
!CM2      non_lake_frac, hcap_soilt(:,m,:), hcon_soilt(:,m,:),                     &
!CM2      sathh_soilt(:,m,:), smcl_soilt(:,m,:), snowdepth_surft,                  &
!CM2      surf_ht_flux_ld, smvcst_soilt(:,m,:), w_flux_soilt(:,m,:),               &
!CM2      sthf_soilt(:,m,:), sthu_soilt(:,m,:), sthu_irr_soilt(:,m,:),             &
!CM2      t_soil_soilt(:,m,:), tsoil_deep_gb )
!CM2  ELSE
!CM2    ! Surface and soil tiles map directly on to each other.
!CM2    DO m = 1, nsoilt
!CM2      n = m
!CM2      CALL soil_htc (                                                          &
!CM2        land_pts, sm_levels, nsurft, soil_pts, timestep, soil_index,           &
!CM2        surft_pts, surft_index, nsnow_surft,                                   &
!CM2        bexp_soilt(:,m,:), dzsoil, tile_frac,                                  &
!CM2        non_lake_frac, hcap_soilt(:,m,:), hcon_soilt(:,m,:),                   &
!CM2        sathh_soilt(:,m,:), smcl_soilt(:,m,:), snowdepth_surft,                &
!CM2        snow_soil_htf(:,n), smvcst_soilt(:,m,:), w_flux_soilt(:,m,:),          &
!CM2        sthf_soilt(:,m,:), sthu_soilt(:,m,:),  sthu_irr_soilt(:,m,:),          &
!CM2        t_soil_soilt(:,m,:), tsoil_deep_gb )
!CM2    END DO
!CM2  END IF
!CM2END IF
!CM2!=============================================================================
!CM2! *END NOTICE REGARDING SOIL TILING**
!CM2!=============================================================================
!CM2
!CM2!-----------------------------------------------------------------------------
!CM2! Update the sub-surface temperatures for land ice.
!CM2!-----------------------------------------------------------------------------
!CM2IF (lice_pts /= 0) THEN
!CM2  IF ( .NOT. l_elev_land_ice) THEN
!CM2    DO m = 1, nsoilt
!CM2      CALL ice_htc (lice_pts, land_pts, sm_levels, timestep, lice_index,       &
!CM2                    dzsoil, surf_ht_flux_ld, t_soil_soilt(:,m,:))
!CM2    END DO
!CM2  ELSE
!CM2    CALL elev_htc (lice_pts, land_pts, nsurft, dzsoil_elev, timestep,          &
!CM2                   lice_index, snow_soil_htf, tsurf_elev_surft)
!CM2  END IF
!CM2END IF
!CABLE_LSM:End
!jhan:up to here
!-----------------------------------------------------------------------------
! Diagnose the available soil moisture in a layer at the surface.
!-----------------------------------------------------------------------------
DO m = 1, nsoilt
  CALL soilmc ( land_pts,sm_levels,soil_pts,soil_index,                        &
                dzsoil,sthu_soilt(:,m,:),smvcst_soilt(:,m,:),                  &
                smvcwt_soilt(:,m,:),smc_soilt(:,m) )
END DO

!-----------------------------------------------------------------------------
! Calculate mean soil temperature and scaled CH4 flux:
!-----------------------------------------------------------------------------
DO m = 1, nsoilt
  DO i = 1,land_pts
    fch4_wetl_soilt(i,m)       = 0.0
    fch4_wetl_cs_soilt(i,m)    = 0.0
    fch4_wetl_npp_soilt(i,m)   = 0.0
    fch4_wetl_resps_soilt(i,m) = 0.0
  END DO
END DO
!!!CABLE_LSM: omit for CABLE 
!!  IF ( l_top .AND. soil_pts /= 0 ) THEN
!!    SELECT CASE ( soil_bgc_model )
!!    CASE ( soil_model_1pool, soil_model_rothc )
!!      IF ( l_ch4_tlayered ) THEN
!!        ! This variable is not used with layered CH4 calc.
!!        tsoil_d_soilt(:,m) = 0.0
!!      ELSE
!!        CALL soilt(land_pts, sm_levels, soil_pts, soil_index,                  &
!!                   dzsoil, t_soil_soilt(:,m,:), tsoil_d_soilt(:,m))
!!      END IF
!!      CALL ch4_wetl(land_pts, sm_levels, soil_pts, dim_cs1, timestep,          &
!!                    l_ch4_tlayered, soil_index,                                &
!!                    t_soil_soilt(:,m,:), tsoil_d_soilt(:,m), cs_ch4_soilt(:,m),&
!!                    resp_s_soilt(:,m,:,:), npp_soilt(:,m), fwetl_soilt(:,m),   &
!!                    sthu_soilt(:,m,:), bexp_soilt(:,m,:),                      &
!!                    fch4_wetl_soilt(:,m), fch4_wetl_cs_soilt(:,m),             &
!!                    fch4_wetl_npp_soilt(:,m), fch4_wetl_resps_soilt(:,m),      &
!!                    substr_ch4, mic_ch4, mic_act_ch4, acclim_ch4,              &
!!                    cs_pool_soilt(:,m,:,:))
!!#if !defined(UM_JULES)
!!      IF (l_imogen) THEN
!!        DO i = 1,land_pts
!!          ! fch4_wetl_acc_soilt in kg/m2 and fch4_wetl_soilt in 10e9kg/m2/s
!!          fch4_wetl_acc_soilt(i,m) = fch4_wetl_acc_soilt(i,m) +                &
!!                (fch4_wetl_soilt(i,m) * to_kg_conversion * timestep_len)
!!        END DO
!!      END IF
!!#endif
!!    END SELECT
!!  END IF
!!END DO  !  tiles
!!
!!IF ( (soil_bgc_model == soil_model_rothc .AND. l_layeredc )                    &
!!     .OR. soil_bgc_model == soil_model_ecosse ) THEN
!!  !-----------------------------------------------------------------------
!!  !accumulate soil temperature for layered soil carbon and nitrogen
!!  !-----------------------------------------------------------------------
!!  DO m = 1, nsoilt
!!    IF (asteps_since_triffid == 1) THEN
!!      t_soil_soilt_acc(:,m,:) = 0.0
!!    END IF
!!    DO j = 1,soil_pts
!!      i = soil_index(j)
!!      t_soil_soilt_acc(i,m,:) = t_soil_soilt_acc(i,m,:) + t_soil_soilt(i,m,:)
!!    END DO
!!  END DO
!!END IF
!CABLE_LSM: omit for CABLE 

!-----------------------------------------------------------------------
! Calculate Nitrogen Leaching
!-----------------------------------------------------------------------
IF (soil_bgc_model == soil_model_rothc .AND. l_nitrogen) THEN
  IF (nsoilt > 1) THEN
    errorstatus = 101
    CALL ereport("check hydrol_jls", errorstatus,                              &
                 "nsoilt>1 and l_nitrogen - n_leach not currently coded")
  END IF
  CALL n_leach(land_pts, asteps_since_triffid, dim_cslayer, timestep,          &
               w_flux_soilt, qbase_l_soilt, sub_surf_roff_gb, smcl_soilt,      &
               smcl_old_soilt, n_inorg_avail_pft, n_inorg_soilt_lyrs,          &
               n_leach_gb_acc, n_leach_soilt)
END IF

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE hydrol_cbl
END MODULE hydrol_mod_cbl
