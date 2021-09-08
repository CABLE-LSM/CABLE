! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYDROL--------------------------------------------------------

! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with MOSES II land surface
!     scheme.

! Documentation : UM Documentation Paper 25

MODULE hydrol_cbl_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HYDROL_MOD'

CONTAINS

SUBROUTINE hydrol_cbl (                                                           &
    lice_pts, lice_index, soil_pts, soil_index,  nsnow_surft,                 &
    land_pts, sm_levels, bexp_soilt, catch_surft, con_rain_land,              &
    ecan_surft, ext_soilt, hcap_soilt, hcon_soilt, ls_rain_land,              &
    con_rainfrac_land,  ls_rainfrac_land,                                     &
    satcon_soilt, sathh_soilt, snowdepth_surft,  snow_soil_htf,               &
    surf_ht_flux_ld, timestep,                                                &
    smvcst_soilt, smvcwt_soilt, canopy_surft,                                 &
    stf_sub_surf_roff, smcl_soilt, sthf_soilt, sthu_soilt,                    &
    t_soil_soilt, tsurf_elev_surft, canopy_gb, smc_soilt, snow_melt,          &
    sub_surf_roff_gb, surf_roff_gb, tot_tfall_gb,                             &
    ! add new inland basin variable
    inlandout_atm_gb, l_inland,                                               &
    ! Additional variables for MOSES II
    nsurft, surft_pts, surft_index,                                           &
    infil_surft,  melt_surft, tile_frac,                                      &
    ! Additional variables required for large-scale hydrology:
    l_top, l_pdm, fexp_soilt, ti_mean_soilt, cs_ch4_soilt, cs_pool_soilt,     &
    dun_roff_soilt, drain_soilt, fsat_soilt, fwetl_soilt, qbase_soilt,        &
    qbase_l_soilt, qbase_zw_soilt, w_flux_soilt,                              &
    zw_soilt, sthzw_soilt, a_fsat_soilt, c_fsat_soilt, a_fwet_soilt,          &
    c_fwet_soilt,                                                             &
    resp_s_soilt, npp_soilt, fch4_wetl_soilt,                                 &
    fch4_wetl_cs_soilt, fch4_wetl_npp_soilt, fch4_wetl_resps_soilt,           &
    dim_cs1, l_soil_sat_down, asteps_since_triffid,                           &
    !CABLE_LSM:additional existing vars
    snow_surft, lying_snow,                                                   &
    air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                         &
    ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                       &
    soil_cbl )

!Use in relevant subroutines
USE calc_zw_inund_mod,        ONLY: calc_zw_inund
USE ch4_wetl_mod,             ONLY: ch4_wetl
USE surf_hyd_mod,             ONLY: surf_hyd
USE soilt_mod,                ONLY: soilt
USE soilmc_mod,               ONLY: soilmc
USE soil_hyd_wt_mod,          ONLY: soil_hyd_wt
USE soil_hyd_update_mod,      ONLY: soil_hyd_update
USE soil_hyd_mod,             ONLY: soil_hyd
USE soil_htc_mod,             ONLY: soil_htc
USE ice_htc_mod,              ONLY: ice_htc
USE calc_baseflow_jules_mod,  ONLY: calc_baseflow_jules
USE n_leach_mod,              ONLY: n_leach
USE prognostics,              ONLY: t_soil_soilt_acc

!Use in relevant variables
USE ancil_info,               ONLY:                                           &
  nsoilt

USE jules_hydrology_mod, ONLY:                                                &
  l_wetland_unfrozen, ti_max, zw_max

USE ancil_info, ONLY: dim_cslayer

USE elev_htc_mod, ONLY: elev_htc

USE jules_soil_biogeochem_mod, ONLY:                                          &
  ! imported scalar parameters
  soil_model_ecosse, soil_model_rothc, soil_model_1pool,                      &
  ! imported scalar variables
  soil_bgc_model, l_ch4_tlayered, l_layeredc

USE jules_soil_mod, ONLY:                                                     &
  dzsoil,                                                                     &
     ! Thicknesses of the soil layers (m).
  dzsoil_elev

USE crop_vars_mod, ONLY: sthu_irr_soilt, frac_irr_soilt, ext_irr_soilt

USE jules_vegetation_mod, ONLY: l_irrig_dmd, l_nitrogen

USE jules_surface_mod, ONLY: l_elev_land_ice

USE ereport_mod, ONLY: ereport

USE water_constants_mod, ONLY:                                                &
  rho_water  ! Density of pure water (kg/m3).

#if !defined(UM_JULES)
USE top_pdm,        ONLY: fch4_wetl_acc_soilt
USE update_mod,     ONLY: l_imogen
USE model_time_mod, ONLY: timestep_len
#endif

!CABLE_LSM: Make avail CABLE hydrology interface subr to CALL
USE cable_hyd_main_mod, ONLY : cable_hyd_main
USE cable_air_type_mod,       ONLY : air_type
USE cable_met_type_mod,       ONLY : met_type
USE cable_radiation_type_mod, ONLY : radiation_type
USE cable_roughness_type_mod, ONLY : roughness_type
USE cable_canopy_type_mod,    ONLY : canopy_type
USE cable_soil_snow_type_mod, ONLY : soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY : bgc_pool_type
USE cable_balances_type_mod,  ONLY : balances_type
USE cable_sum_flux_type_mod,  ONLY : sum_flux_type
USE cable_params_mod,         ONLY : veg_parameter_type
USE cable_params_mod,         ONLY : soil_parameter_type


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  lice_pts,                                                                   &
    ! Number of land ice points.
  land_pts,                                                                   &
    ! Number of gridpoints.
  sm_levels,                                                                  &
    ! Number of soil moisture levels.
  soil_pts,                                                                   &
    ! Number of soil points.
  nsurft,                                                                     &
    ! Number of tiles
  dim_cs1
    ! Number of soil carbon pools

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  timestep
    ! Model timestep (s).

LOGICAL , INTENT(IN) ::                                                       &
  stf_sub_surf_roff,                                                          &
    ! Stash flag for sub-surface runoff.
  l_top,                                                                      &
    ! Flag for TOPMODEL-based hydrology.
  l_pdm,                                                                      &
    ! Flag for PDM hydrology.
  l_soil_sat_down
    ! Switch controlling direction of movement of
    ! soil moisture in excess of saturation.

!-----------------------------------------------------------------------------
! Array arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  lice_index(land_pts),                                                       &
    ! Array of land ice points.
  soil_index(land_pts),                                                       &
    ! Array of soil points.
  nsnow_surft(land_pts,nsurft)
    ! Number of snow layers

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  bexp_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Clapp-Hornberger exponent.
  catch_surft(land_pts,nsurft),                                               &
    ! Canopy/surface capacity of land tiles (kg/m2).
  con_rain_land(land_pts),                                                    &
    ! Convective rain (kg/m2/s).
  ecan_surft(land_pts,nsurft),                                                &
    ! Canopy evaporation from land tiles (kg/m2/s).
  ext_soilt(land_pts,nsoilt,sm_levels),                                       &
    ! Extraction of water from each soil layer (kg/m2/s).
  hcap_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Soil heat capacity (J/K/m3).
  hcon_soilt(land_pts,nsoilt,0:sm_levels),                                    &
    ! Soil thermal conductivity (W/m/K).
  ls_rain_land(land_pts),                                                     &
    ! Large-scale rain (kg/m2/s).
  con_rainfrac_land(land_pts),                                                &
    ! Convective rain fraction
  ls_rainfrac_land(land_pts),                                                 &
    ! large scale rain fraction
  satcon_soilt(land_pts,nsoilt,0:sm_levels),                                  &
    ! Saturated hydraulic conductivity (kg/m2/s).
  sathh_soilt(land_pts,nsoilt,sm_levels),                                     &
    ! Saturated soil water pressure (m).
  snow_melt(land_pts),                                                        &
    ! Snowmelt (kg/m2/s).
  snowdepth_surft(land_pts,nsurft),                                           &
    ! Snow depth (on ground) (m)
  snow_soil_htf(land_pts,nsurft),                                             &
    ! Tiled snowpack-> soil heat flux.
  surf_ht_flux_ld(land_pts),                                                  &
    ! Net downward surface heat flux (W/m2).
  smvcst_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric soil moisture concentration at saturation (m3 H2O/m3 soil).
  smvcwt_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric soil moisture concentration below which
    ! stomata close (m3 H2O/m3 soil).
  fexp_soilt(land_pts,nsoilt),                                                &
    ! Decay factor in Sat. Conductivity in deep LSH/TOPMODEL layer.
  ti_mean_soilt(land_pts,nsoilt),                                             &
    ! Mean topographic index.
  npp_soilt(land_pts,nsoilt),                                                 &
    ! Gridbox mean net primary productivity (kg C/m2/s).
  a_fsat_soilt(land_pts,nsoilt),                                              &
    ! Fitting parameter for Fsat in LSH model.
  c_fsat_soilt(land_pts,nsoilt),                                              &
    ! Fitting parameter for Fsat in LSH model.
  a_fwet_soilt(land_pts,nsoilt),                                              &
    ! Fitting parameter for Fwet in LSH model.
  c_fwet_soilt(land_pts,nsoilt)
    ! Fitting parameter for Fwet in LSH model.

!-----------------------------------------------------------------------------
! Array arguments with INTENT(INOUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  cs_ch4_soilt(land_pts,nsoilt),                                              &
    ! Soil carbon used in CH4 wetlands if TRIFFID is switched off (kg C/m2).
  canopy_surft(land_pts,nsurft),                                              &
    ! Canopy water content for land tiles (kg/m2).
  smcl_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Soil moisture content of each
!                          !       layer (kg/m2).
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
  t_soil_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Sub-surface temperatures (K).
  tsurf_elev_surft(land_pts,nsurft),                                          &
    ! Tiled sub-surface temperatures (K).
  cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1),                         &
    ! Soil carbon (kg C/m2).
!                          !   For RothC (dim_cs1=4), the pools
!                          !   are DPM, RPM, biomass and humus.
  resp_s_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1),                          &
    ! Soil respiration in pools (kg C/m2/s).
  fsat_soilt(land_pts,nsoilt),                                                &
     ! Surface saturation fraction.
  fwetl_soilt(land_pts,nsoilt),                                               &
     ! Wetland fraction.
  zw_soilt(land_pts,nsoilt),                                                  &
     ! Water table depth (m).
  sthzw_soilt(land_pts,nsoilt)
     ! Soil moisture fraction in deep LSH/TOPMODEL layer.

!-----------------------------------------------------------------------------
! Array arguments with INTENT(OUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  canopy_gb(land_pts),                                                        &
    ! Gridbox canopy water content (kg/m2).
  smc_soilt(land_pts,nsoilt),                                                 &
    ! Available soil moisture in a layer at the surface (kg/m2)
  sub_surf_roff_gb(land_pts),                                                 &
    ! Sub-surface runoff (kg/m2/s).
  surf_roff_gb(land_pts),                                                     &
    ! Surface runoff (kg/m2/s).
  tot_tfall_gb(land_pts),                                                     &
    ! Total throughfall (kg/m2/s).
  dun_roff_soilt(land_pts,nsoilt),                                            &
    ! Dunne part of sfc runoff (kg/m2/s).
  qbase_soilt(land_pts,nsoilt),                                               &
    ! Base flow (kg/m2/s).
  qbase_l_soilt(land_pts,nsoilt,sm_levels+1),                                 &
    ! Base flow from each level (kg/m2/s).
  qbase_zw_soilt(land_pts,nsoilt),                                            &
    ! Base flow from deep LSH/TOPMODEL layer (kg/m2/s).
  w_flux_soilt(land_pts,nsoilt,0:sm_levels),                                  &
    ! Fluxes of water between layers (kg/m2/s).
  drain_soilt(land_pts,nsoilt),                                               &
    ! Drainage out of sm_levels'th level (kg/m2/s).
  fch4_wetl_soilt(land_pts,nsoilt),                                           &
    ! Scaled wetland methane flux (default substrate) for use in
    ! atmos chemistry model (10^-9 kg C/m2/s).
  fch4_wetl_cs_soilt(land_pts,nsoilt),                                        &
    ! Scaled methane flux (soil carbon substrate) (kg C/m2/s).
  fch4_wetl_npp_soilt(land_pts,nsoilt),                                       &
    ! Scaled methane flux (npp substrate) (kg C/m2/s).
  fch4_wetl_resps_soilt(land_pts,nsoilt)
    ! Scaled methane flux (soil respiration substrate) (kg C/m2/s).

! Additional variables for MOSES II
INTEGER, INTENT(IN) ::                                                        &
  surft_pts(nsurft),                                                          &
    ! Number of tile points.
  surft_index(land_pts,nsurft)
    ! Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  infil_surft(land_pts,nsurft),                                               &
    ! Maximum surface infiltration
  melt_surft(land_pts,nsurft),                                                &
    ! Snowmelt on tiles (kg/m2/s).
  tile_frac(land_pts,nsurft),                                                 &
    ! Tile fractions.
! Declare variable for inland basin outflow
  inlandout_atm_gb(land_pts)
    ! IN TRIP INLAND BASIN OUTFLOW FOR LAND POINTS ONLY,kg/m2/s=mm.

LOGICAL, INTENT(IN) ::                                                        &
  l_inland
    ! True if re-routing inland basin flow to soil moisture.

!-----------------------------------------------------------------------------
! Local scalars:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i, j,                                                                       &
  n,                                                                          &
    ! Counter for soil level.
  m,                                                                          &
    ! Counter for soil tile.
  errorstatus

!-----------------------------------------------------------------------------
! Local arrays:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  dsmc_dt_soilt(land_pts,nsoilt),                                             &
    ! Rate of change of soil moisture due to water falling onto the
    ! surface after surface runoff (kg/m2/s).
  ksz_soilt(land_pts,nsoilt,0:sm_levels),                                     &
    ! Saturated hydraulic conductivity in layer (kg/m2/s).
  qbase_unfr_soilt(land_pts,nsoilt),                                          &
    ! Base flow in unfrozen soil (kg/m2/s).
  qbase_l_unfr_soilt(land_pts,nsoilt,sm_levels+1),                            &
    ! As qbase_l but for unfrozen soil (kg/m2/s).
  top_crit_soilt(land_pts,nsoilt),                                            &
    ! Critical TI when ZW <=0.0
  dumtop_crit_soilt(land_pts,nsoilt),                                         &
    ! Dummy for top_crit_soilt
  dumsthf_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! Dummy Frozen soil moisture content of each layer as a fraction of
    ! saturation (always set to 0).
  zdepth(0:sm_levels),                                                        &
    ! Lower soil layer boundary depth (m).
  tsoil_d_soilt(land_pts,nsoilt),                                             &
    ! Soil temperature in the top metre
  zw_inund_soilt(land_pts,nsoilt),                                            &
    ! Water table depth used
  wutot_soilt(land_pts,nsoilt),                                               &
    ! Ratio of unfrozen to total soil moisture at ZW.
  surf_roff_inc_soilt(land_pts,nsoilt),                                       &
    ! Increment to tiled surface runoff (kg m-2 s-1).
  surf_roff_soilt(land_pts,nsoilt),                                           &
    ! Soil-tiled contributions to surface runoff (kg m-2 s-1).
  sub_surf_roff_soilt(land_pts,nsoilt)
    ! Soil-tiled contributions to subsurface runoff (kg m-2 s-1).

REAL(KIND=real_jlslsm), PARAMETER :: to_kg_conversion = 1.0e-9
    ! multiplier for converting to kgC for wetland CH4 and IMOGEN

! Variables required for irrigation code
REAL(KIND=real_jlslsm) ::                                                     &
  w_flux_irr_soilt(land_pts,nsoilt,0:sm_levels),                              &
    ! The fluxes of water between layers in irrigated fraction (kg/m2/s).
  w_flux_nir_soilt(land_pts,nsoilt,0:sm_levels),                              &
    ! The fluxes of water between layers in non-irrigated fraction (kg/m2/s).
  smcl_irr_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Total soil moisture contents of each layer in irrigated
    ! fraction (kg/m2).
  smcl_nir_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Total soil moisture contents of each layer in non-irrigated
    ! fraction (kg/m2).
  sthu_nir_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation in irrigated fraction.
  ext_nir_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! Extraction of water from each soil layer in non-irrigated fraction
    ! (kg/m2/s).
  smclsat_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! The saturation moisture content of each layer (kg/m2).
  smclzw_soilt(land_pts,nsoilt),                                              &
    ! moisture content in deep layer(kg/m2).
  smclsatzw_soilt(land_pts,nsoilt)
    ! moisture content in deep layer (kg/m2).

!CABLE_LSM:
real ::  lying_snow(land_pts)
real ::  snow_surft(land_pts,nsurft)
TYPE(air_type),       INTENT(inout) :: air_cbl
TYPE(met_type),       INTENT(inout) :: met_cbl
TYPE(radiation_type), INTENT(inout) :: rad_cbl
TYPE(roughness_type), INTENT(inout) :: rough_cbl
TYPE(canopy_type),    INTENT(inout) :: canopy_cbl
TYPE(soil_snow_type), INTENT(inout) :: ssnow_cbl
TYPE(bgc_pool_type),  INTENT(inout) :: bgc_cbl
TYPE(balances_type),  INTENT(inout) :: bal_cbl
TYPE(sum_flux_type),  INTENT(inout) :: sum_flux_cbl
TYPE(veg_parameter_type),  INTENT(inout) :: veg_cbl
TYPE(soil_parameter_type), INTENT(inout) :: soil_cbl

LOGICAL :: l_triffid
INTEGER :: asteps_since_triffid

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HYDROL'

! End of header---------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate soil carbon for use in the wetland CH4 scheme only
! (only used if single-pool C model is used):
IF ( soil_bgc_model == soil_model_1pool ) THEN
  DO m = 1,nsoilt
    DO j = 1,soil_pts
      i = soil_index(j)
      cs_ch4_soilt(i,m) = 0.0
      DO n = 1,dim_cslayer
        cs_ch4_soilt(i,m) = cs_ch4_soilt(i,m) + cs_pool_soilt(i,m,n,1)
      END DO
    END DO
  END DO
END IF

!CABLE_LSM: call CABLE hydrology
call cable_hyd_main( land_pts, nsurft, lying_snow, snow_surft, surf_roff_gb,   &
                     sub_surf_roff_gb, tot_tfall_gb, air_cbl, met_cbl, rad_cbl,& 
                     rough_cbl, canopy_cbl, ssnow_cbl, bgc_cbl, bal_cbl,       &
                     sum_flux_cbl, veg_cbl, soil_cbl )
! Initialise w_flux variables that are used in irrigation code
w_flux_soilt(:,:,:)     = 0.0
w_flux_irr_soilt(:,:,:) = 0.0 ! to prevent random values reported over areas
                    ! that are not included as soil points (i.e. ice points)

!-----------------------------------------------------------------------------
! Set up variables required for LSH scheme:
!-----------------------------------------------------------------------------
zdepth(:) = 0.0

DO n = 1,sm_levels
  zdepth(n) = zdepth(n-1) + dzsoil(n)
END DO

! Initialise runoff increment.
surf_roff_inc_soilt(:,:) = 0.0

!-----------------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------------
!CABLE_LSM: omit for CABLE: BUT initialize dsmc_dt for fpe0, otherwise?  
!C!CALL surf_hyd (land_pts, nsurft, surft_pts, surft_index,                      &
!C!               catch_surft, ecan_surft, tile_frac, infil_surft, con_rain_land,&
!C!               ls_rain_land,  con_rainfrac_land,  ls_rainfrac_land,           &
!C!               melt_surft, snow_melt, timestep,                               &
!C!               canopy_surft, canopy_gb, dsmc_dt_soilt,                        &
!C!               l_top, l_pdm, sm_levels, soil_pts, soil_index,                 &
!C!               surf_roff_gb, tot_tfall_gb,                                    &
!C!               dun_roff_soilt, fsat_soilt, smvcst_soilt, sthu_soilt,          &
!C!               sthf_soilt, surf_roff_soilt)

!-----------------------------------------------------------------------------
! Specify the reduction of hydraulic conductivity with depth.
! Initialiase base flow to zero.
!-----------------------------------------------------------------------------
DO m = 1, nsoilt
  DO n = 0,sm_levels
    !CDIR NODEP
    DO j = 1,soil_pts
      i = soil_index(j)
      ksz_soilt(i,m,n) = satcon_soilt(i,m,n)
    END DO
  END DO
END DO

DO m = 1, nsoilt
  DO n = 1,sm_levels
    !CDIR NODEP
    DO j = 1,soil_pts
      smclsat_soilt(soil_index(j),m,n)      = rho_water * dzsoil(n) *         &
                                              smvcst_soilt(soil_index(j),m,n)
      qbase_l_soilt(soil_index(j),m,n)      = 0.0
      qbase_l_unfr_soilt(soil_index(j),m,n) = 0.0
      dumsthf_soilt(soil_index(j),m,n)      = 0.0
    END DO
  END DO
END DO

DO m = 1, nsoilt
  DO i = 1,land_pts
    qbase_soilt(i,m)      = 0.0
    qbase_zw_soilt(i,m)   = 0.0
    wutot_soilt(i,m)      = 0.0
    drain_soilt(i,m)      = 0.0
    qbase_unfr_soilt(i,m) = 0.0
    zw_inund_soilt(i,m)   = 0.0
  END DO
END DO

!CABLE_LSM: omit for CABLE
!C!IF (l_top) THEN
!C!  IF (soil_pts /= 0) THEN
!C!    DO m = 1, nsoilt
!C!      CALL calc_baseflow_jules(                                               &
!C!        soil_pts, soil_index, land_pts, sm_levels,                            &
!C!        zdepth, ksz_soilt(:,m,:),                                             &
!C!        bexp_soilt(:,m,:), fexp_soilt(:,m), ti_mean_soilt(:,m), zw_soilt(:,m),&
!C!        sthf_soilt(:,m,:),                                                    &
!C!        dumtop_crit_soilt(:,m), qbase_soilt(:,m), qbase_l_soilt(:,m,:))
!C!    END DO
!C!
!C!    IF (l_wetland_unfrozen) THEN
!C!      DO m = 1,nsoilt
!C!        CALL calc_zw_inund(land_pts, sm_levels, soil_pts, soil_index, zdepth, &
!C!          bexp_soilt(:,m,:), sathh_soilt(:,m,:), smclsat_soilt(:,m,:),        &
!C!          smcl_soilt(:,m,:), sthu_soilt(:,m,:), sthzw_soilt(:,m),             &
!C!          zw_soilt(:,m), zw_inund_soilt(:,m), wutot_soilt(:,m))
!C!
!C!        ! Now call again to get the unfrozen equivalents to calculate fsat and
!C!        ! fwet:
!C!        CALL calc_baseflow_jules(                                             &
!C!          soil_pts, soil_index, land_pts, sm_levels,                          &
!C!          zdepth, ksz_soilt(:,m,:),                                           &
!C!          bexp_soilt(:,m,:), fexp_soilt(:,m), ti_mean_soilt(:,m),             &
!C!          zw_inund_soilt(:,m), dumsthf_soilt(:,m,:),                          &
!C!          top_crit_soilt(:,m), qbase_unfr_soilt(:,m),                         &
!C!          qbase_l_unfr_soilt(:,m,:))
!C!      END DO
!C!    ELSE
!C!      top_crit_soilt(:,:) = dumtop_crit_soilt(:,:)
!C!    END IF
!C!
!C!  END IF
!C!END IF  !  l_top
!CABLE_LSM:End

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
!CABLE_LSM: omit for CABLE 
!C!IF (soil_pts /= 0) THEN
!C!
!C!  IF (l_irrig_dmd) THEN
!C!
!C!    !-------------------------------------------------------------------------
!C!    ! If l_irrig_dmd = TRUE, call soil_hyd separately for irrigated and
!C!    ! non-irrigated fraction
!C!    ! afterwards, call soil_hyd ONLY to update water table/drainage with
!C!    ! gridbox total w_flux_soilt, smcl_soilt
!C!    !-------------------------------------------------------------------------
!C!
!C!    ! Split into irrigated / non-irrigated fraction.
!C!    DO n = 1,sm_levels
!C!      DO m = 1, nsoilt
!C!        DO j = 1,soil_pts
!C!          i = soil_index(j)
!C!          smclsat_soilt(i,m,n) = rho_water * dzsoil(n) * smvcst_soilt(i,m,n)
!C!
!C!          ! hadrd - gridbox sthu_soilt is assumed to be combination of
!C!          ! sthu_soilt of non-irrigated fraction and sthu_irr_soilt, i.e.
!C!          ! sthu_soilt = frac_irr_soilt*sthu_irr_soilt + (1-frac_irr_soilt)
!C!          !              *sthu_nir_soilt
!C!          IF ( frac_irr_soilt(i,m) < 1.0 ) THEN
!C!            sthu_nir_soilt(i,m,n) = (sthu_soilt(i,m,n) - frac_irr_soilt(i,m)  &
!C!                                    * sthu_irr_soilt(i,m,n))                  &
!C!                                    / (1.0 - frac_irr_soilt(i,m))
!C!            ext_nir_soilt(i,m,n)  = (ext_soilt(i,m,n)  - frac_irr_soilt(i,m)  &
!C!                                    * ext_irr_soilt(i,m,n))                   &
!C!                                    / (1.0 - frac_irr_soilt(i,m))
!C!          ELSE
!C!            sthu_nir_soilt(i,m,n) = sthu_soilt(i,m,n)
!C!            ext_nir_soilt(i,m,n)  = ext_soilt(i,m,n)
!C!          END IF
!C!
!C!          smcl_irr_soilt(i,m,n) = smcl_soilt(i,m,n)                           &
!C!                                  + (sthu_irr_soilt(i,m,n)                    &
!C!                                     - sthu_soilt(i,m,n))                     &
!C!                                  * smclsat_soilt(i,m,n)
!C!          smcl_nir_soilt(i,m,n) = smcl_soilt(i,m,n)                           &
!C!                                  + (sthu_nir_soilt(i,m,n)                    &
!C!                                     - sthu_soilt(i,m,n))                     &
!C!                                  * smclsat_soilt(i,m,n)
!C!        END DO
!C!      END DO
!C!    END DO
!C!
!C!    DO m = 1, nsoilt
!C!      ! First call soil_hyd for non-irrigated fraction.
!C!      CALL soil_hyd (                                                         &
!C!        land_pts, sm_levels, soil_pts, soil_index, bexp_soilt(:,m,:), dzsoil, &
!C!        ext_nir_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),           &
!C!        sathh_soilt(:,m,:), timestep, smvcst_soilt(:,m,:),                    &
!C!        smcl_nir_soilt(:,m,:), sthu_nir_soilt(:,m,:), w_flux_nir_soilt(:,m,:),&
!C!        sthzw_soilt(:,m), zdepth, qbase_l_soilt(:,m,:), l_top,                &
!C!        l_soil_sat_down, smclzw_soilt(:,m), smclsatzw_soilt(:,m),             &
!C!        smclsat_soilt(:,m,:))
!C!
!C!      ! Next call soil_hyd for irrigated fraction.
!C!      CALL soil_hyd (                                                         &
!C!        land_pts, sm_levels, soil_pts, soil_index, bexp_soilt(:,m,:), dzsoil, &
!C!        ext_irr_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),           &
!C!        sathh_soilt(:,m,:), timestep, smvcst_soilt(:,m,:),                    &
!C!        smcl_irr_soilt(:,m,:), sthu_irr_soilt(:,m,:), w_flux_irr_soilt(:,m,:),&
!C!        sthzw_soilt(:,m), zdepth, qbase_l_soilt(:,m,:), l_top,                &
!C!        l_soil_sat_down, smclzw_soilt(:,m), smclsatzw_soilt(:,m),             &
!C!        smclsat_soilt(:,m,:))
!C!    END DO
!C!
!C!    ! Re-calculate total grid box soil moisture.
!C!    ! hadrd - perhaps this could be done more efficiently with WHERE
!C!    DO m = 1,nsoilt
!C!      DO n = 0,sm_levels
!C!        DO j = 1,soil_pts
!C!          i = soil_index(j)
!C!
!C!          ! Ensure sensible values if irrigation fraction is very small:
!C!          IF ( frac_irr_soilt(i,m) <= EPSILON(1.0) ) THEN
!C!            w_flux_irr_soilt(i,m,n) = 0.0
!C!            w_flux_soilt(i,m,n)     = w_flux_nir_soilt(i,m,n)
!C!            IF ( n >= 1 ) THEN
!C!              sthu_irr_soilt(i,m,n) = sthu_nir_soilt(i,m,n)
!C!              smcl_soilt(i,m,n)     = smcl_nir_soilt(i,m,n)
!C!              sthu_soilt(i,m,n)     = sthu_nir_soilt(i,m,n)
!C!            END IF
!C!          ELSE
!C!            w_flux_soilt(i,m,n) = frac_irr_soilt(i,m)                         &
!C!                                  * w_flux_irr_soilt(i,m,n)                   &
!C!                                  + ( 1.0 - frac_irr_soilt(i,m) )             &
!C!                                  * w_flux_nir_soilt(i,m,n)
!C!            IF ( n >= 1 ) THEN
!C!              smcl_soilt(i,m,n) = frac_irr_soilt(i,m) * smcl_irr_soilt(i,m,n) &
!C!                                  + ( 1.0 - frac_irr_soilt(i,m) )             &
!C!                                  * smcl_nir_soilt(i,m,n)
!C!              sthu_soilt(i,m,n) = frac_irr_soilt(i,m) * sthu_irr_soilt(i,m,n) &
!C!                                  + ( 1.0 - frac_irr_soilt(i,m) )             &
!C!                                  * sthu_nir_soilt(i,m,n)
!C!            END IF
!C!          END IF
!C!        END DO  !  soil points
!C!      END DO  !  layers
!C!    END DO  !  tiles
!C!
!C!    DO m = 1, nsoilt
!C!      CALL soil_hyd_update(land_pts, sm_levels, soil_pts, soil_index, dzsoil, &
!C!                           smvcst_soilt(:,m,:), zdepth, smclzw_soilt(:,m),    &
!C!                           sthzw_soilt(:,m), smclsat_soilt(:,m,:),            &
!C!                           smclsatzw_soilt(:,m))
!C!    END DO
!C!
!C!  ELSE
!C!    ! .NOT. l_irrig_dmd
!C!
!C!    DO m = 1, nsoilt
!C!      CALL soil_hyd (                                                         &
!C!        land_pts, sm_levels, soil_pts, soil_index, bexp_soilt(:,m,:), dzsoil, &
!C!        ext_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),               &
!C!        sathh_soilt(:,m,:), timestep, smvcst_soilt(:,m,:),                    &
!C!        smcl_soilt(:,m,:), sthu_soilt(:,m,:), w_flux_soilt(:,m,:),            &
!C!        sthzw_soilt(:,m), zdepth, qbase_l_soilt(:,m,:), l_top,                &
!C!        l_soil_sat_down, smclzw_soilt(:,m), smclsatzw_soilt(:,m),             &
!C!        smclsat_soilt(:,m,:))
!C!    END DO
!C!
!C!  END IF  !  l_irrig_dmd
!C!
!C!  IF (nsoilt == 1) THEN
!C!    !To maintain bit-comparability, we need to call with the _gb version of
!C!    !surface runoff.
!C!    m = 1
!C!    CALL soil_hyd_wt (                                                        &
!C!      land_pts, sm_levels, soil_pts, soil_index,                              &
!C!      bexp_soilt(:,m,:), dsmc_dt_soilt(:,m), sathh_soilt(:,m,:), timestep,    &
!C!      smvcst_soilt(:,m,:), sub_surf_roff_soilt(:,m), smcl_soilt(:,m,:),       &
!C!      surf_roff_gb,                                                           &
!C!      w_flux_soilt(:,m,:), stf_sub_surf_roff, zw_soilt(:,m),                  &
!C!      sthzw_soilt(:,m), qbase_soilt(:,m), qbase_l_soilt(:,m,:),               &
!C!      drain_soilt(:,m), l_top, smclzw_soilt(:,m), smclsatzw_soilt(:,m),       &
!C!      smclsat_soilt(:,m,:), surf_roff_inc_soilt(:,m) )
!C!
!C!    ! For a single soil tile, simply copy across to the output variable.
!C!    sub_surf_roff_gb(:) = sub_surf_roff_soilt(:,m)
!C!    ! Update the tiled surface runoff.
!C!    surf_roff_soilt(:,m) = surf_roff_soilt(:,m) + surf_roff_inc_soilt(:,m)
!C!  ELSE
!C!
!C!    ! Initialise output variable
!C!    sub_surf_roff_gb(:) = 0.0
!C!    DO m = 1, nsoilt
!C!      CALL soil_hyd_wt (                                                      &
!C!        land_pts, sm_levels, soil_pts, soil_index,                            &
!C!        bexp_soilt(:,m,:), dsmc_dt_soilt(:,m), sathh_soilt(:,m,:), timestep,  &
!C!        smvcst_soilt(:,m,:), sub_surf_roff_soilt(:,m), smcl_soilt(:,m,:),     &
!C!        surf_roff_soilt(:,m),                                                 &
!C!        w_flux_soilt(:,m,:), stf_sub_surf_roff, zw_soilt(:,m),                &
!C!        sthzw_soilt(:,m), qbase_soilt(:,m), qbase_l_soilt(:,m,:),             &
!C!        drain_soilt(:,m), l_top, smclzw_soilt(:,m), smclsatzw_soilt(:,m),     &
!C!        smclsat_soilt(:,m,:), surf_roff_inc_soilt(:,m) )
!C!
!C!      ! For multiple soil tiles, add up the contributions, allowing for frac
!C!      sub_surf_roff_gb(:) = sub_surf_roff_gb(:)                               &
!C!                            + ( tile_frac(:,m) * sub_surf_roff_soilt(:,m))
!C!      surf_roff_gb(:) = surf_roff_gb(:) + tile_frac(:,m)                      &
!C!                                          * surf_roff_inc_soilt(:,m)
!C!
!C!    END DO
!C!  END IF !nsoilt == 1
!C!
!C!  !---------------------------------------------------------------------------
!C!  ! Calculate surface saturation and wetland fractions:
!C!  !---------------------------------------------------------------------------
!C!  IF (l_top) THEN
!C!    DO m = 1, nsoilt
!C!
!C!      DO i = 1,land_pts
!C!        fsat_soilt(i,m)  = 0.0
!C!        fwetl_soilt(i,m) = 0.0
!C!
!C!        ! Zero soil porosity over land ice:
!C!        IF (smvcst_soilt(i,m,sm_levels) <= 0.0) THEN
!C!          zw_soilt(i,m) = zw_max
!C!        END IF
!C!      END DO
!C!
!C!      DO j = 1,soil_pts
!C!        i = soil_index(j)
!C!        qbase_zw_soilt(i,m) = qbase_l_soilt(i,m,sm_levels+1)
!C!
!C!        !Now use fit for fsat_soilt and fwet:
!C!        IF (l_wetland_unfrozen) THEN
!C!          fsat_soilt(i,m)  = wutot_soilt(i,m) * a_fsat_soilt(i,m)             &
!C!                             * EXP(-c_fsat_soilt(i,m) * top_crit_soilt(i,m))
!C!          fwetl_soilt(i,m) = wutot_soilt(i,m) * a_fwet_soilt(i,m)             &
!C!                             * EXP(-c_fwet_soilt(i,m) * top_crit_soilt(i,m))
!C!        ELSE
!C!          fsat_soilt(i,m)  = a_fsat_soilt(i,m)                                &
!C!                             * EXP(-c_fsat_soilt(i,m) * top_crit_soilt(i,m))
!C!          fwetl_soilt(i,m) = a_fwet_soilt(i,m)                                &
!C!                             * EXP(-c_fwet_soilt(i,m) * top_crit_soilt(i,m))
!C!        END IF
!C!
!C!        IF (top_crit_soilt(i,m) >= ti_max) THEN
!C!          fsat_soilt(i,m)  = 0.0
!C!          fwetl_soilt(i,m) = 0.0
!C!        END IF
!C!
!C!      END DO
!C!    END DO
!C!  END IF  !  l_top
!C!
!C!ELSE ! soil pts
!C!
!C!  !---------------------------------------------------------------------------
!C!  ! If required by STASH flag and there are no soil points,
!C!  ! set sub-surface runoff to zero.
!C!  !---------------------------------------------------------------------------
!C!  IF (stf_sub_surf_roff) THEN
!C!    DO i = 1,land_pts
!C!      sub_surf_roff_gb(i) = 0.0
!C!    END DO
!C!  END IF
!C!
!C!END IF  !  soil_pts
!C!
!C!!-----------------------------------------------------------------------------
!C!! Update the soil temperatures and the frozen moisture fractions
!C!!-----------------------------------------------------------------------------
!C!
!C!!=============================================================================
!C!! *NOTICE REGARDING SOIL TILING**
!C!!
!C!!The following section facilitates the use of soil tiling. As implemented,
!C!!there are two soil tiling options:
!C!!
!C!!nsoilt == 1
!C!!Operate as with a single soil tile, functionally identical to JULES upto
!C!! at least vn4.7 (Oct 2016)
!C!! This means that a soilt variable being passed 'up' to the surface is
!C!! broadcast to the surft variable (with weighting by frac if requred)
!C!!
!C!!nsoilt > 1
!C!!Operate with nsoilt = nsurft, with a direct mapping between them
!C!! This means that a soilt variable being passed 'up' to the surface is simply
!C!! copied into the surft variable
!C!!
!C!! This will need to be refactored for other tiling approaches. This note
!C!! will be replicated elsewhere in the code as required
!C!!
!C!!These comments apply until **END NOTICE REGARDING SOIL TILING**
!C!!=============================================================================
!C!IF (soil_pts /= 0) THEN
!C!  ! When using soil tiling, we can use the surface tiled version of
!C!  ! surf_ht_flux_ld, snow_soil_htf. The _ld version is a gridbox mean
!C!  ! calculated in snow and passed through.
!C!  IF (nsoilt == 1) THEN
!C!    m = 1
!C!    CALL soil_htc (                                                           &
!C!      land_pts, sm_levels, nsurft, soil_pts, soil_index,                      &
!C!      surft_pts, surft_index, nsnow_surft,                                    &
!C!      bexp_soilt(:,m,:), dzsoil, tile_frac, hcap_soilt(:,m,:),                &
!C!      hcon_soilt(:,m,:),                                                      &
!C!      sathh_soilt(:,m,:), surf_ht_flux_ld, timestep, smvcst_soilt(:,m,:),     &
!C!      w_flux_soilt(:,m,:), sthu_irr_soilt(:,m,:), smcl_soilt(:,m,:),          &
!C!      snowdepth_surft, sthu_soilt(:,m,:), sthf_soilt(:,m,:),                  &
!C!      t_soil_soilt(:,m,:))
!C!  ELSE
!C!    ! Surface and soil tiles map directly on to each other.
!C!    DO m = 1, nsoilt
!C!      n = m
!C!      CALL soil_htc (                                                         &
!C!        land_pts, sm_levels, nsurft, soil_pts, soil_index,                    &
!C!        surft_pts, surft_index, nsnow_surft,                                  &
!C!        bexp_soilt(:,m,:), dzsoil, tile_frac, hcap_soilt(:,m,:),              &
!C!        hcon_soilt(:,m,:),                                                    &
!C!        sathh_soilt(:,m,:), snow_soil_htf(:,n), timestep, smvcst_soilt(:,m,:),&
!C!        w_flux_soilt(:,m,:), sthu_irr_soilt(:,m,:), smcl_soilt(:,m,:),        &
!C!        snowdepth_surft, sthu_soilt(:,m,:), sthf_soilt(:,m,:),                &
!C!        t_soil_soilt(:,m,:))
!C!    END DO
!C!  END IF
!C!END IF
!C!!=============================================================================
!C!! *END NOTICE REGARDING SOIL TILING**
!C!!=============================================================================
!C!
!C!!-----------------------------------------------------------------------------
!C!! Update the sub-surface temperatures for land ice.
!C!!-----------------------------------------------------------------------------
!C!IF (lice_pts /= 0) THEN
!C!  IF ( .NOT. l_elev_land_ice) THEN
!C!    DO m = 1, nsoilt
!C!      CALL ice_htc (land_pts, sm_levels, lice_pts, lice_index, dzsoil,        &
!C!                    surf_ht_flux_ld, timestep,                                &
!C!                    t_soil_soilt(:,m,:))
!C!    END DO
!C!  ELSE
!C!    CALL elev_htc (land_pts, lice_pts, lice_index, nsurft,                    &
!C!                   dzsoil_elev, snow_soil_htf, timestep,                      &
!C!                   tsurf_elev_surft)
!C!  END IF
!C!END IF
!CABLE_LSM:End

!-----------------------------------------------------------------------------
! Diagnose the available soil moisture in a layer at the surface.
!-----------------------------------------------------------------------------
DO m = 1, nsoilt
  CALL soilmc ( land_pts,sm_levels,soil_pts,soil_index,                       &
                dzsoil,sthu_soilt(:,m,:),smvcst_soilt(:,m,:),                 &
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
  !CABLE_LSM: omit for CABLE 
  IF ( l_top .AND. soil_pts /= 0 ) THEN
   SELECT CASE ( soil_bgc_model )
   CASE ( soil_model_1pool, soil_model_rothc )
  !C!    IF ( l_ch4_tlayered ) THEN
  !C!      ! This variable is not used with layered CH4 calc.
  !C!      tsoil_d_soilt(:,m) = 0.0
  !C!    ELSE
  !C!      CALL soilt(land_pts, sm_levels, soil_pts, soil_index,                 &
  !C!                 dzsoil, t_soil_soilt(:,m,:), tsoil_d_soilt(:,m))
  !C!    END IF
  !C!    CALL ch4_wetl(land_pts, soil_pts, dim_cs1, soil_index, sm_levels,       &
  !C!                  tsoil_d_soilt(:,m), cs_ch4_soilt(:,m),                    &
  !C!                  t_soil_soilt(:,m,:), cs_pool_soilt(:,m,:,:),              &
  !C!                  resp_s_soilt(:,m,:,:), npp_soilt(:,m), fwetl_soilt(:,m),  &
  !C!                  fch4_wetl_soilt(:,m), fch4_wetl_cs_soilt(:,m),            &
  !C!                  fch4_wetl_npp_soilt(:,m), fch4_wetl_resps_soilt(:,m),     &
  !C!                  timestep, l_ch4_tlayered)
#if !defined(UM_JULES)
      IF (l_imogen) THEN
        DO i = 1,land_pts
          ! fch4_wetl_acc_soilt in kg/m2 and fch4_wetl_soilt in 10e9kg/m2/s
          fch4_wetl_acc_soilt(i,m) = fch4_wetl_acc_soilt(i,m) +               &
                (fch4_wetl_soilt(i,m) * to_kg_conversion * timestep_len)
        END DO
      END IF
#endif
    END SELECT
  END IF
END DO  !  tiles

IF ( (soil_bgc_model == soil_model_rothc .AND. l_layeredc )                   &
     .OR. soil_bgc_model == soil_model_ecosse ) THEN
  !-----------------------------------------------------------------------
  !accumulate soil temperature for layered soil carbon and nitrogen
  !-----------------------------------------------------------------------
  DO m = 1, nsoilt
    IF (asteps_since_triffid == 1) THEN
      t_soil_soilt_acc(:,m,:) = 0.0
    END IF
    DO j = 1,soil_pts
      i = soil_index(j)
      t_soil_soilt_acc(i,m,:) = t_soil_soilt_acc(i,m,:) + t_soil_soilt(i,m,:)
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate Nitrogen Leaching
!-----------------------------------------------------------------------
IF (soil_bgc_model == soil_model_rothc .AND. l_nitrogen) THEN
  IF (nsoilt > 1) THEN
    errorstatus = 101
    CALL ereport("check hydrol_jls", errorstatus,                             &
                 "nsoilt>1 and l_nitrogen - n_leach not currently coded")
  END IF
  CALL n_leach(land_pts, timestep, smcl_soilt, w_flux_soilt,                  &
               sub_surf_roff_gb, qbase_l_soilt, asteps_since_triffid)
END IF

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE hydrol_cbl
END MODULE hydrol_cbl_mod
