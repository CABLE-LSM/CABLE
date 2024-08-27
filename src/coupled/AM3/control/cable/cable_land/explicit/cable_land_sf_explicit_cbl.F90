! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE cable_land_sf_explicit_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
                  ModuleName='CABLE_LAND_SF_EXPLICIT_MOD'

CONTAINS
!  SUBROUTINE CABLE_LAND_SF_EXPLICIT ---------------------------------
!
!  Purpose: Calculate explicit surface fluxes of heat, moisture and
!           momentum over land. Also calculates surface exchange
!           coefficients required for implicit update of surface
!           fluxes and surface information required by the
!           explicit boundary layer routine
!
!
!  Documentation: UMDP 24.
!
!  !CM3#55 - this routine is derived from jules_land_sf_explicit and
!  and modified for the purposes of linking CABLE and JULES together
!  in the surf_couple_explicit section - work commencing 6/12/2023
!
!  All edits will be accompanied by comments starting with !CM3#55-x
!  referring to the git issue - where possible the included number x will
!  refer to sub-issue as per the git discussion
!
!---------------------------------------------------------------------
!    Arguments :-
SUBROUTINE cable_land_sf_explicit (                                            &
! IN date-related values
 curr_day_number,                                                              &
! IN values defining field dimensions and subset to be processed :
 land_pts,                                                                     &
! IN  parameters for iterative SISL scheme
 numcycles, cycleno,                                                           &
! IN parameters required from boundary-layer scheme :
 bq_1,bt_1,z1_uv,z1_uv_top,z1_tq,z1_tq_top,qw_1,tl_1,                          &
! IN soil/vegetation/land surface data :
 land_index,nsurft,sm_levels,canopy,catch,catch_snow,hcon_soilt,               &
 ho2r2_orog, flandg,                                                           &
 snow_surft,sil_orog_land,smvccl_soilt,smvcst_soilt,smvcwt_soilt,sthf_soilt,   &
 sthu_soilt,z0_surft,z0h_surft_bare, z0m_soil_in,                              &
! IN input data from the wave model
 charnock_w,                                                                   &
! IN everything not covered so far :
 pstar,lw_down,sw_surft,zh,ddmfx,                                              &
 co2_mmr,co2_3d,l_co2_interactive,l_phenol,                                    &
 asteps_since_triffid,cs_pool_soilt,veg_state,frac,canht_pft,                  &
 photosynth_act_rad,lai_pft,                                                   &
 l_mr_physics,t_soil_soilt,tsurf_elev_surft,tstar_surft,z_land,                &
 albsoil_soilt,cos_zenith_angle,l_aero_classic,l_dust,l_dust_diag,             &
 clay_soilt,o3, l_emis_surft_set,                                              &
! INOUT diagnostics
 sf_diag,                                                                      &
! INOUT data :
 emis_surft,gs,gc_corr,g_leaf_acc,npp_pft_acc,resp_w_pft_acc,resp_s_acc_soilt, &
 rhostar,fqw_1,ftl_1,t1_sd,q1_sd,vshr,vshr_land,                               &
! OUT Diagnostic not requiring STASH flags :
 ftl_surft,                                                                    &
! OUT variables for message passing
 rhokm_land, cdr10m,                                                           &
! OUT data required for mineral dust scheme
 u_s_std_surft,                                                                &
! OUT data required elsewhere in boundary layer or surface code
 alpha1,ashtf_prime_surft,fqw_surft,epot_surft,fraca,                          &
 resfs,resft,rhokh_surft,dtstar_surft,z0h_surft, z0m_surft,                    &
 chr1p5m,smc_soilt,hcons_soilt,gpp,npp,resp_p,g_leaf,gpp_pft,npp_pft,          &
 resp_p_pft,resp_s_soilt,resp_s_tot_soilt,resp_l_pft,resp_r_pft,               &
 resp_w_pft,n_leaf,n_root,n_stem,lai_bal,gc_surft,canhc_surft,wt_ext_surft,    &
 flake,surft_index,surft_pts,tile_frac,fsmc_pft,emis_soil,                     &
! OUT required for classic aerosols
 cd_land,rib_surft,ch_surft_classic,cd_std_classic,                            &
! OUT required for sea and sea-ice calculations
 l_cdr10m_snow,                                                                &
 !New arguments replacing USE statements
 !Fluxes (IN)
 t_home_gb, t_growth_gb,                                                       &
 !urban_param (IN)
 emisr_gb, emisw_gb, hwr_gb,                                                   &
 !jules_mod (IN OUT)
 albobs_scaling_surft,                                                         &
 !jules_chemvars_mod (OUT)
 isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                          &
 methanol_gb, methanol_pft, acetone_gb, acetone_pft,                           &
 !trif_vars_mod (OUT)
 fapar_diag_pft, apar_diag_pft, apar_diag_gb, gpp_gb_acc, gpp_pft_acc,         &
 !crop_vars_mod (IN)
 rootc_cpft, sthu_irr_soilt, frac_irr_soilt, frac_irr_surft, dvi_cpft,         &
 !crop_vars_mod (IN OUT)
 resfs_irr_surft,                                                              &
 !crop_vars_mod (OUT)
 gs_irr_surft, smc_irr_soilt, wt_ext_irr_surft, gc_irr_surft,                  &
 !p_s_parms (IN)
 bexp_soilt, sathh_soilt, v_close_pft, v_open_pft,                             &
 !urban_param (IN)
 wrr_gb,                                                                       &
 !Fluxes (IN OUT)
 anthrop_heat_surft,                                                           &
 !prognostics (IN)
 nsnow_surft, sice_surft, sliq_surft, snowdepth_surft,                         &
                        tsnow_surft, ds_surft,                                 &
 !c_elevate (OUT)
 surf_hgt_surft, lw_down_elevcorr_surft,                                       &
 !jules_mod (OUT)
 snowdep_surft,                                                                &
 !urban_param (IN)
 hgt_gb, disp_gb,                                                              &
 !lake_mod (IN)
 lake_t_ice_gb,lake_t_mxl_gb, lake_h_ice_gb,lake_depth_gb,                     &
 g_dt_gb, non_lake_frac,                                                       &
 !lake_mod (OUT)
 nusselt_gb, ts1_lake_gb, hcon_lake,                                           &
 !ancil_info
 l_lice_point, l_soil_point,                                                   &
 !jules_surface_types (IN)
 diff_frac,                                                                    &
 !chemvars (OUT)
 flux_o3_pft, fo3_pft,                                                         &
 !CABLE_LSM:CM3
 progs_cbl, work_cbl, pars_io_cbl, progs_cnp,                                  &
 mype, timestep_number, satcon_soilt,                                          &
 latitude, longitude, u_s, ls_rain, ls_snow )

USE ancil_info,              ONLY: dim_cslayer, l_lice_surft, nsoilt, rad_nband
USE atm_fields_bounds_mod,      ONLY: pdims_s, pdims, tdims
USE bl_option_mod,              ONLY: l_quick_ap2
USE c_elevate,                  ONLY: l_elev_absolute_height
USE c_z0h_z0m,                  ONLY: z0h_z0m, z0h_z0m_classic
USE calc_air_dens_mod,          ONLY: calc_air_dens
USE can_drag_mod,               ONLY: can_drag_z0, can_drag_phi_m_h
USE csigma,                     ONLY: sbcon
USE dust_param,                 ONLY: z0_soil
USE elevate_mod,                ONLY: elevate
USE fcdch_mod,                  ONLY: fcdch
USE gen_anthrop_heat_mod,       ONLY: generate_anthropogenic_heat
USE heat_con_mod,               ONLY: heat_con
USE physiol_mod,                ONLY: physiol
USE planet_constants_mod,       ONLY: cp, vkman, r, c_virtual,epsil=>repsilon
USE qsat_mod,                   ONLY: qsat, qsat_mix
USE sf_diags_mod,               ONLY: strnewsfdiag
USE sf_flux_mod_cbl,            ONLY: sf_flux_cbl
USE sf_orog_mod,                ONLY: sf_orog
USE sf_resist_mod,              ONLY: sf_resist
USE sf_rib_mod,                 ONLY: sf_rib
USE sfl_int_mod,                ONLY: sfl_int
USE snowtherm_mod,              ONLY: snowtherm
USE solinc_data,                ONLY: sky, l_skyview
USE stdev1_mod,                 ONLY: stdev1
USE stochastic_physics_run_mod, ONLY: l_rp2, i_rp_scheme, i_rp2b, z0hm_pft_rp
USE theta_field_sizes,          ONLY: t_i_length,t_j_length
USE tilepts_mod,                ONLY: tilepts
USE timestep_mod,               ONLY: timestep
USE urban_param_mod,            ONLY: z0m_mat
USE urbanz0_mod,                ONLY: urbanz0
USE veg_param,                  ONLY: secs_per_360days
USE veg3_field_mod,             ONLY: veg_state_type
USE water_constants_mod,        ONLY: lc, rho_ice, tm

USE jules_soil_biogeochem_mod, ONLY:                                           &
! imported scalar parameters
   soil_model_rothc,                                                           &
! imported scalar variables (IN)
   soil_bgc_model

USE jules_soil_mod, ONLY: dzsoil, dzsoil_elev, hcice, hcwat, hcondeep

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype,                         &
                                   urban_canyon, urban_roof, soil, lake, ncpft

#if defined(UM_JULES)
USE atm_step_local, ONLY:  dim_cs1, co2_dim_len,co2_dim_row
#else
USE ancil_info, ONLY:  dim_cs1, co2_dim_len, co2_dim_row
#endif

USE jules_snow_mod, ONLY: cansnowtile                                          &
                          ,rho_snow_const                                      &
                          ,snow_hcon                                           &
                          ,l_snowdep_surf                                      &
                          ,l_snow_nocan_hc                                     &
                          ,nsmax                                               &
                          ,unload_rate_u

USE jules_surface_mod, ONLY: l_aggregate, formdrag, l_anthrop_heat_src,        &
                             i_aggregate_opt, cor_mo_iter,                     &
                             use_correct_ustar, iscrntdiag,                    &
                             l_flake_model,l_elev_lw_down,                     &
                             l_mo_buoyancy_calc, effective_z0,                 &
                             IP_ScrnDecpl2, IP_ScrnDecpl3,                     &
                             l_vary_z0m_soil, l_elev_land_ice, ls

USE jules_vegetation_mod, ONLY: can_model, can_rad_mod, ilayers, l_triffid,    &
                                l_vegdrag_surft

USE jules_irrig_mod, ONLY: l_irrig_dmd

USE jules_sea_seaice_mod, ONLY: l_ctile, charnock, ip_ss_solid

USE jules_urban_mod, ONLY: l_moruses_rough_surft, l_moruses_storage

USE jules_science_fixes_mod, ONLY: ctile_orog_fix, correct_sea_adjust_land,    &
                                   l_fix_wind_snow, l_accurate_rho,            &
                                   l_fix_moruses_roof_rad_coupling

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE progs_cbl_vars_mod,      ONLY: progs_cbl_vars_type ! CABLE intro-ed progs
USE work_vars_mod_cbl,       ONLY: work_vars_type      ! some kept thru timestep
USE params_io_mod_cbl,       ONLY: params_io_data_type
USE progs_cnp_vars_mod,      ONLY: progs_cnp_vars_type ! CASA-CNP intro-ed progs
USE cable_explicit_main_mod, ONLY: cable_explicit_main

IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs :-
!-----------------------------------------------------------------------
! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.
INTEGER, INTENT(IN) ::                                                         &
 curr_day_number,                                                              &
            ! IN current day of year
 land_pts,                                                                     &
            ! IN No of land points being processed.
 numcycles,                                                                    &
            ! Number of cycles (iterations) for iterative SISL.
 cycleno
            ! Iteration no

! Defining vertical grid of model atmosphere.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! IN A buoyancy parameter
                             !    (beta q tilde).
,bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! IN A buoyancy parameter
                             !    (beta T tilde).
,z1_uv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN Height of lowest uv level (m).
,z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN Height of lowest tq level (m).
                             !    Note, if the grid used is
                             !    staggered in the vertical,
                             !    Z1_UV and Z1_TQ can be
                             !    different.
,qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! IN Total water content
,tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! IN Ice/liquid water temperature

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
! Charnock's coefficient from wave model

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
                    z1_uv_top(tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end)
                             ! Height of top of lowest uv-layer
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
                    z1_tq_top(tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end)
                             ! Height of top of lowest Tq-layer

! (c) Soil/vegetation/land surface parameters (mostly constant).
LOGICAL, INTENT(IN) ::                                                         &
 l_co2_interactive                                                             &
                             ! IN Switch for 3D CO2 field
,l_phenol
                             ! IN Indicates whether phenology
                             !    in use

INTEGER, INTENT(IN) ::                                                         &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
                             !    point in ROW_LENGTH,ROWS is the
                             !    land point.

INTEGER, INTENT(IN) ::                                                         &
 sm_levels                                                                     &
                             ! IN No. of soil moisture levels
,nsurft                                                                        &
                             ! IN No. of land-surface tiles
,asteps_since_triffid
                             ! IN Number of atmospheric
                             !    timesteps since last call
                             !    to TRIFFID.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 canopy(land_pts,nsurft)                                                       &
                             ! IN Surface/canopy water for
                             !    snow-free land tiles (kg/m2)
,catch(land_pts,nsurft)                                                        &
                             ! IN Surface/canopy water capacity
                             !    of snow-free land tiles (kg/m2).
,catch_snow(land_pts,nsurft)                                                   &
                             ! IN Snow interception capacity of
                             !    tiles (kg/m2).
,hcon_soilt(land_pts,nsoilt)                                                   &
                             ! IN Soil thermal conductivity
                             !    (W/m/K).
,snow_surft(land_pts,nsurft)                                                   &
                             ! IN Lying snow on tiles (kg/m2)
,smvccl_soilt(land_pts,nsoilt,sm_levels)                                       &
                             ! IN Critical volumetric SMC
                             !    (cubic m per cubic m of soil).
,smvcst_soilt(land_pts,nsoilt,sm_levels)                                       &
                             ! IN Volumetric saturation point
                             !    (m3/m3 of soil).
,smvcwt_soilt(land_pts,nsoilt,sm_levels)                                       &
                             ! IN Volumetric wilting point
                             !    (cubic m per cubic m of soil).
,sthf_soilt(land_pts,nsoilt,sm_levels)                                         &
                             ! IN Frozen soil moisture content of
                             !    each layer as a fraction of
                             !    saturation.
,sthu_soilt(land_pts,nsoilt,sm_levels)                                         &
                             ! IN Unfrozen soil moisture content
                             !    of each layer as a fraction of
                             !    saturation.
,z0_surft(land_pts,nsurft)                                                     &
                             ! IN Tile roughness lengths (m).
,z0h_surft_bare(land_pts,nsurft)                                               &
                             ! IN Tile thermal roughness lengths
                             !    without snow cover(m).
,z0m_soil_in(land_pts)                                                         &
                             ! IN bare soil momentum z0 (m).
,sil_orog_land(land_pts)                                                       &
                             ! IN Silhouette area of unresolved
                             !    orography per unit horizontal
                             !    area on land points only.
,ho2r2_orog(land_pts)                                                          &
                             ! IN Standard Deviation of orography.
                             !    equivilent to peak to trough
                             !    height of unresolved orography
,flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! IN Land fraction on all tiles.
                             !    divided by 2SQRT(2) on land
                             !    points only (m)

! (f) Atmospheric + any other data not covered so far, incl control.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)                    &
                             ! IN Surface pressure (Pascals).
,lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! IN Surface downward LW radiation
                             !    (W/m2).
,sw_surft(land_pts,nsurft)                                                     &
                             ! IN Surface net SW radiation on
                             !    land tiles (W/m2).
,zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                       &
                             ! IN Height above surface of top of
                             !    boundary layer (metres).
,ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN Convective downdraught
                             !    mass-flux at cloud base
,co2_mmr                                                                       &
                             ! IN CO2 Mass Mixing Ratio
,co2_3d(co2_dim_len,co2_dim_row)                                               &
                             ! IN 3D CO2 field if required.
,cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)                            &
                             ! IN Soil carbon (kg C/m2).
,frac(land_pts,ntype)                                                          &
                             ! IN Fractions of surface types.
,canht_pft(land_pts,npft)                                                      &
                             ! IN Canopy height (m)
,photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! IN Net downward shortwave radiation
                             !    in band 1 (w/m2).
,lai_pft(land_pts,npft)                                                        &
                             ! IN Leaf area index
,t_soil_soilt(land_pts,nsoilt,sm_levels)                                       &
                             ! IN Soil temperatures (K).
,tsurf_elev_surft(land_pts,nsurft)                                             &
                             ! IN Tiled ice sub-surface temperature (K)
,z_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN Land height (m).
,albsoil_soilt(land_pts,nsoilt)                                                &
                             ! IN Soil albedo.
, cos_zenith_angle(tdims%i_start:tdims%i_end,                                  &
                   tdims%j_start:tdims%j_end)                                  &
                             ! IN Cosine of the zenith angle
,clay_soilt(land_pts,nsoilt,dim_cslayer)                                       &
                             ! IN Soil clay fraction.
,o3(land_pts)
                             ! IN Surface ozone concentration (ppb).

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                       &
tstar_surft(land_pts,nsurft)
                             ! IN Surface tile temperatures

LOGICAL, INTENT(IN) ::                                                         &
 l_aero_classic                                                                &
                             ! IN switch for using CLASSIC aerosol
                             !    scheme
,l_dust                                                                        &
                             ! IN switch for mineral dust
,l_dust_diag                                                                   &
                             ! IN Switch for diagnostic mineral dust
                             !    lifting
,l_mr_physics                                                                  &
                             ! IN Switch for when mixing ratios are used
,l_emis_surft_set(nsurft)
                             ! IN Switch for varying grey surface emissivity

!-----------------------------------------------------------------------
!  In/outs :-
!-----------------------------------------------------------------------
!Diagnostics
TYPE (strnewsfdiag), INTENT(IN OUT) :: sf_diag

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 emis_surft(land_pts,nsurft)                                                   &
                              ! INOUT Emissivity for land tiles
,gs(land_pts)                                                                  &
                              ! INOUT "Stomatal" conductance to
                              !        evaporation (m/s).
,g_leaf_acc(land_pts,npft)                                                     &
                              ! INOUT Accumulated G_LEAF
,npp_pft_acc(land_pts,npft)                                                    &
                              ! INOUT Accumulated NPP_pft
,resp_w_pft_acc(land_pts,npft)                                                 &
                              ! INOUT Accum RESP_W_pft
,resp_s_acc_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)                         &
                              ! INOUT Accumulated RESP_S
,rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                              ! INOUT Surface air density
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                              ! INOUT Moisture flux between layers
                              !       (kg per square metre per sec).
                              !       FQW(,1) is total water flux
                              !       from surface, 'E'.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                              ! INOUT FTL(,K) contains net turbulent
                              !       sensible heat flux into layer K
                              !       from below; so FTL(,1) is the
                              !       surface sensible heat, H.(W/m2)
,t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                              ! OUT Standard deviation of turbulent
                              !     fluctuations of layer 1 temp;
                              !     used in initiating convection.
,q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                              ! OUT Standard deviation of turbulent
                              !     flucs of layer 1 humidity;
                              !     used in initiating convection.
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                              ! OUT Magnitude of surface-to-lowest
                              !     atm level wind shear (m per s).
,vshr_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                              ! OUT Magnitude of surface-to-lowest
                              !     atm level wind shear (m per s).

!-----------------------------------------------------------------------
!  Outputs :-
!-----------------------------------------------------------------------
!-1 Diagnostic (or effectively so - includes coupled model requisites):-
INTEGER, INTENT(OUT) ::                                                        &
 surft_index(land_pts,ntype)                                                   &
                              ! OUT Index of tile points
,surft_pts(ntype)             ! OUT Number of tile points

!  (a) Calculated anyway (use STASH space from higher level) :-

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 ftl_surft(land_pts,nsurft)                                                    &
                             ! OUT Surface FTL for land tiles
,u_s_std_surft(land_pts,nsurft)                                                &
                             ! OUT Surface friction velocity
                             !     (standard value)
                             !     for mineral dust
,emis_soil(land_pts)                                                           &
                             ! OUT Emissivity of underlying soil
,rhokm_land(pdims_s%i_start:pdims_s%i_end,                                     &
            pdims_s%j_start:pdims_s%j_end),                                    &
 cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

!-2 Genuinely output, needed by other atmospheric routines :-

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 alpha1(land_pts,nsurft)                                                       &
                             ! OUT Mean gradient of saturated
                             !     specific humidity with respect
                             !     to temperature between the
                             !     bottom model layer and tile
                             !     surfaces
,ashtf_prime_surft(land_pts,nsurft)                                            &
                             ! OUT Coefficient to calculate
                             !     surface heat flux into land
                             !     tiles.
,fqw_surft(land_pts,nsurft)                                                    &
                             ! OUT Surface FQW for land tiles
,epot_surft(land_pts,nsurft)                                                   &
                             ! OUT Local EPOT for land tiles.
,fraca(land_pts,nsurft)                                                        &
                             ! OUT Fraction of surface moisture
                             !     flux with only aerodynamic
                             !     resistance for snow-free land
                             !     tiles.
,resfs(land_pts,nsurft)                                                        &
                             ! OUT Combined soil, stomatal
                             !     and aerodynamic resistance
                             !     factor for fraction (1-FRACA)
                             !     of snow-free land tiles.
,resft(land_pts,nsurft)                                                        &
                             ! OUT Total resistance factor.
                             !     FRACA+(1-FRACA)*RESFS for
                             !     snow-free land, 1 for snow.
,rhokh_surft(land_pts,nsurft)                                                  &
                             ! OUT Surface exchange coefficients
                             !     for land tiles
,dtstar_surft(land_pts,nsurft)                                                 &
                             ! OUT Change in TSTAR over timestep
                             !     for land tiles
,z0h_surft(land_pts,nsurft)                                                    &
                             ! OUT Tile roughness lengths for heat
                             !     and moisture (m).
,z0m_surft(land_pts,nsurft)                                                    &
                             ! OUT Tile roughness lengths for
                             !     momentum.
,chr1p5m(land_pts,nsurft)                                                      &
                             ! OUT Ratio of coefffs for
                             !     calculation of 1.5m temp for
                             !     land tiles.
,smc_soilt(land_pts,nsoilt)                                                    &
                             ! OUT Available moisture in the
                             !     soil profile (mm).
,hcons_soilt(land_pts,nsoilt)                                                  &
                             ! OUT Soil thermal conductivity
                             !     including water and ice
,gpp(land_pts)                                                                 &
                             ! OUT Gross primary productivity
                             !     (kg C/m2/s).
,npp(land_pts)                                                                 &
                             ! OUT Net primary productivity
                             !     (kg C/m2/s).
,resp_p(land_pts)                                                              &
                             ! OUT Plant respiration (kg C/m2/s).
,g_leaf(land_pts,npft)                                                         &
                             ! OUT Leaf turnover rate (/360days).
,gpp_pft(land_pts,npft)                                                        &
                             ! OUT Gross primary productivity
                             !     on PFTs (kg C/m2/s).
,npp_pft(land_pts,npft)                                                        &
                             ! OUT Net primary productivity
                             !     (kg C/m2/s).
,resp_p_pft(land_pts,npft)                                                     &
                             ! OUT Plant respiration on PFTs
                             !     (kg C/m2/s).
,resp_s_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)                             &
                          ! OUT Soil respiration (kg C/m2/s).
,resp_s_tot_soilt(land_pts,nsoilt)                                             &
                             ! OUT Total soil respiration
                             ! (kg C/m2/s).
,resp_l_pft(land_pts,npft)                                                     &
                             ! OUT Leaf maintenance respiration
                             !     (kg C/m2/s).
,resp_r_pft(land_pts,npft)                                                     &
                             ! OUT Root maintenance respiration
                             !     (kg C/m2/s).
,resp_w_pft(land_pts,npft)                                                     &
                             ! OUT Wood maintenance respiration
                             !     (kg C/m2/s).
,n_leaf(land_pts,npft)                                                         &
                             ! OUT Leaf N content scaled by LAI
                             !     (kg N/m2).
,n_root(land_pts,npft)                                                         &
                             ! OUT Root N content scaled by LAI_bal
                             !     (kg N/m2).
,n_stem(land_pts,npft)                                                         &
                             ! OUT Stem N content scaled by LAI_bal
                             !     (kg N/m2).
,lai_bal(land_pts,npft)                                                        &
                             ! OUT LAI_bal
,gc_surft(land_pts,nsurft)                                                     &
                             ! OUT "Stomatal" conductance to
                             !      evaporation for land tiles
                             !      (m/s).
,canhc_surft(land_pts,nsurft)                                                  &
                             ! OUT Areal heat capacity of canopy
                             !    for land tiles (J/K/m2).
,wt_ext_surft(land_pts,sm_levels,nsurft)                                       &
                             ! OUT Fraction of evapotranspiration
                             !     which is extracted from each
                             !     soil layer by each tile.
,flake(land_pts,nsurft)                                                        &
                             ! OUT Lake fraction.
,tile_frac(land_pts,nsurft)                                                    &
                             ! OUT Tile fractions including
                             !     snow cover in the ice tile.
,fsmc_pft(land_pts,npft)                                                       &
                             ! OUT Moisture availability factor.
,gc_corr(land_pts,npft)
                             ! OUT "Stomatal" conductance
                             !     without bare soil evaporation

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 cd_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! OUT Bulk transfer coefficient for
                             !     momentum over land.
,rib_surft(land_pts,nsurft)                                                    &
                             ! OUT RIB for land tiles.
,ch_surft_classic(land_pts,nsurft)                                             &
                             ! OUT Bulk transfer coefficient for
                             !     heat for aerosol deposition.
,cd_std_classic(land_pts,nsurft)
                             ! OUT Bulk transfer coefficient for
                             !     momentum for aerosol deposition.
LOGICAL, INTENT(OUT) ::                                                        &
 l_cdr10m_snow
                             ! OUT Flag indicating if cdr10m
                             !     (an interpolation coefficient) is
                             !     to be calculated for use with
                             !     snow unloading.

!New arguments replacing USE statements
!urban_param
REAL(KIND=real_jlslsm), INTENT(IN) :: emisr_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: emisw_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: hwr_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: wrr_gb(land_pts)

!jules_mod
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
          albobs_scaling_surft(land_pts,ntype,rad_nband)

!p_s_parms (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: bexp_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: sathh_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: v_close_pft(land_pts,sm_levels,npft)
REAL(KIND=real_jlslsm), INTENT(IN) :: v_open_pft(land_pts,sm_levels,npft)

!crop_vars_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: rootc_cpft(land_pts,ncpft)
REAL(KIND=real_jlslsm), INTENT(IN) :: sthu_irr_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: frac_irr_soilt(land_pts,nsoilt)
REAL(KIND=real_jlslsm), INTENT(IN) :: frac_irr_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) :: dvi_cpft(land_pts,ncpft)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: resfs_irr_surft(land_pts,nsurft)

!Fluxes (IN OUT)
REAL(KIND=real_jlslsm), INTENT(IN OUT)    :: anthrop_heat_surft(land_pts,nsurft)

!Fluxes
REAL(KIND=real_jlslsm), INTENT(IN) :: t_home_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: t_growth_gb(land_pts)

!veg_state
TYPE(veg_state_type),   INTENT(IN OUT) :: veg_state

!jules_chemvars_mod
REAL(KIND=real_jlslsm), INTENT(OUT) :: isoprene_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: terpene_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: methanol_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: acetone_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: isoprene_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: terpene_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: methanol_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: acetone_pft(land_pts,npft)

!trif_vars_mod
REAL(KIND=real_jlslsm), INTENT(OUT) :: fapar_diag_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: apar_diag_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: apar_diag_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gpp_gb_acc(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gpp_pft_acc(land_pts,npft)

!crop_vars_mod (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gs_irr_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: smc_irr_soilt(land_pts,nsoilt)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
        wt_ext_irr_surft(land_pts,sm_levels,nsurft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gc_irr_surft(land_pts,nsurft)

!prognostics (IN)
INTEGER, INTENT(IN) :: nsnow_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) :: sice_surft(land_pts,nsurft,nsmax),       &
                                      sliq_surft(land_pts,nsurft,nsmax),       &
                                      snowdepth_surft(land_pts,nsurft),        &
                                      tsnow_surft(land_pts,nsurft,nsmax),      &
                                      ds_surft(land_pts,nsurft,nsmax)

!c_elevate (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: surf_hgt_surft(land_pts,nsurft),        &
                                       lw_down_elevcorr_surft(land_pts,nsurft)

!jules_mod (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: snowdep_surft(land_pts,nsurft)

!urban_param (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: hgt_gb(land_pts),                        &
                                      disp_gb(land_pts)

!lake_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: lake_t_ice_gb(land_pts),                 &
                                      lake_t_mxl_gb(land_pts),                 &
                                      lake_h_ice_gb(land_pts),                 &
                                      lake_depth_gb(land_pts),                 &
                                      g_dt_gb(land_pts),                       &
                                      non_lake_frac(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: ts1_lake_gb(land_pts),                  &
                                      nusselt_gb(land_pts),                    &
                                      hcon_lake(land_pts)

!ancil_info (IN)
LOGICAL, INTENT(IN) :: l_lice_point(land_pts)
LOGICAL, INTENT(IN) :: l_soil_point(land_pts)

!JULES surface_types_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: diff_frac(t_i_length * t_j_length)

!chemvars (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: flux_o3_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: fo3_pft(land_pts,npft)


!-----------------------------------------------------------------------
! LOCAL variables
!-----------------------------------------------------------------------
!  Workspace :-
REAL(KIND=real_jlslsm) :: work_clay         ! working variable

REAL(KIND=real_jlslsm) ::                                                      &
 vfrac_surft(land_pts,nsurft)                                                  &
                             ! Fractional canopy coverage for
                             ! land tiles.
,radnet_surft(land_pts,nsurft)                                                 &
                             ! Surface net radiation on tiles
,csnow(land_pts,nsmax)                                                         &
                             ! Areal heat capacity of snow (J/K/m2)
,ksnow(land_pts,nsmax)                                                         &
                             ! Thermal conductivity of snow (W/m/K)
,hcons_snow(land_pts,nsurft)                                                   &
                             ! Snow thermal conductivity
,resp_frac(land_pts,dim_cslayer)                                               &
                             ! respired fraction of RESP_S
,gc_stom_surft(land_pts,nsurft)
                             ! canopy conductance

REAL(KIND=real_jlslsm) ::                                                      &
 lh0                         ! Latent heat for snow free surface
                             !   =LS for sea-ice, =LC otherwise


!  Workspace for sea-ice and marginal ice zone
REAL(KIND=real_jlslsm) ::                                                      &
 ch_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Bulk transfer coefficient for
                             !      het and moisture over land.

!  Workspace for land tiles
REAL(KIND=real_jlslsm) ::                                                      &
 cd_std(land_pts,nsurft)                                                       &
                             ! Local drag coefficient for calc
                             ! of interpolation coefficient
,cd_surft(land_pts,nsurft)                                                     &
                             ! Drag coefficient
,ch_surft(land_pts,nsurft)                                                     &
                             ! Transfer coefficient for heat and
                             ! moisture
,chn(land_pts,nsurft)                                                          &
                             ! Neutral value of CH.
,dq(land_pts)                                                                  &
                             ! Sp humidity difference between
                             ! surface and lowest atmospheric lev
,epdt(land_pts)                                                                &
                             ! "Potential" Evaporation * Timestep
,pstar_land(land_pts)                                                          &
                             ! Surface pressure for land points.
,qstar_surft(land_pts,nsurft)                                                  &
                             !Surface saturated sp humidity.
,rhokh_can(land_pts,nsurft)                                                    &
                             ! Exchange coefficient for canopy
                             ! air to surface
,rhokm_1_surft(land_pts,nsurft)                                                &
                             ! Momentum exchange coefficient.
,tsurf(land_pts,nsurft)                                                        &
                             ! Surface layer temp (snow or soil) (
,lake_ice_mid_temp(land_pts)                                                   &
                             ! Median temperature of the lake ice (K)
,dzsurf(land_pts,nsurft)                                                       &
                             ! Surface layer thickness
                             ! (snow or soil) (m)
,canhc_surf(land_pts,nsurft)                                                   &
                             ! Surface layer thickness
                             ! (snow or soil) (m)
,hcons_surf(land_pts,nsurft)                                                   &
                             ! Thermal conductivity
                             ! (snow or soil)  (W/m/K)
,wind_profile_factor(land_pts,nsurft)                                          &
                             ! For transforming effective surface
                             ! transfer coefficients to those
                             ! excluding form drag.
,z0m_eff_surft(land_pts,nsurft)                                                &
                             ! Effective momentum roughness length
,db_surft(land_pts,nsurft)                                                     &
                             ! Buoyancy difference for surface
                             ! tile
,v_s_surft(land_pts,nsurft)                                                    &
                             ! Surface layer scaling velocity
                             ! for tiles (m/s).
,v_s_std(land_pts,nsurft)                                                      &
                             ! Surface layer scaling velocity
                             ! for tiles excluding orographic
                             ! form drag (m/s).
,u_s_iter_surft(land_pts,nsurft)                                               &
                             ! Scaling velocity from middle of
                             ! MO scheme - picked up in error by
                             ! dust code!
,recip_l_mo_surft(land_pts,nsurft)                                             &
                             ! Reciprocal of the Monin-Obukhov
                             ! length for tiles (m^-1).
,z0m_soil(land_pts,nsurft)                                                     &
                             ! Bare soil momentum roughness length
                             ! for use in 1 tile dust scheme
,z0h_soil(land_pts,nsurft)                                                     &
                             ! Bare soil roughness length for heat
                             ! for use in 1 tile dust scheme
,wind_profile_fac_soil(land_pts,nsurft)                                        &
                             ! Equivalent of wind_profile_factor
                             ! for use in 1 tile dust scheme
,cd_surft_soil(land_pts,nsurft), ch_surft_soil(land_pts,nsurft)                &
,cd_std_soil(land_pts,nsurft), v_s_surft_soil(land_pts,nsurft)                 &
,recip_l_mo_surft_soil(land_pts,nsurft)                                        &
                             ! Dummy output variables from extra
                             ! call to fcdch needed for
                             ! 1 tile dust scheme
,v_s_std_soil(land_pts,nsurft)                                                 &
                             ! Bare soil surface layer scaling
                             ! velocity for tiles excluding
                             ! orographic form drag (m/s)
                             ! for use in 1 tile dust scheme
,u_s_iter_soil(land_pts,nsurft)                                                &
                             ! Bare soil scaling velocity from
                             ! middle of MO scheme for use in
                             ! 1 tile dust scheme - picked up in
                             ! error by dust code
,z0h_surft_classic(land_pts,nsurft)                                            &
                             ! z0h to be used in calculation for
                             ! CLASSIC aerosol deposition
,cd_surft_classic(land_pts,nsurft)                                             &
,v_s_surft_classic(land_pts,nsurft)                                            &
,recip_l_mo_surft_classic(land_pts,nsurft)                                     &
,v_s_std_classic(land_pts,nsurft)                                              &
,u_s_iter_classic(land_pts,nsurft)                                             &
                             ! Dummy output variables from extra
                             ! call to fcdch needed for aerosol
                             ! deposition with different z0h
,t_elev(land_pts,nsurft)                                                       &
                             ! Temperature at elevated height (k)
,q_elev(land_pts,nsurft)                                                       &
                             ! Specific humidity at elevated
                             !     height (kg per kg air)
,qs1_elev(land_pts,nsurft)                                                     &
                             ! Saturated specific humidity at elev
                             !     height (kg per kg air)
,scaling_urban(land_pts,nsurft)                                                &
                             ! MORUSES: ground heat flux scaling;
                             ! canyon tile only coupled to soil
,zdt_surft(land_pts,nsurft)                                                    &
                             ! Difference between the canopy height and
                             ! displacement height (m)
,phi_m(land_pts)                                                               &
                             ! Monin-Obukhov stability function for momentum
                             ! integrated to the model's lowest wind level.
,phi_h(land_pts)                                                               &
                             ! Monin-Obukhov stability function for scalars
                             ! integrated to the model's lowest temperature
                             ! and humidity level.
,rhostar_mom(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! Surface air density for momentum
,ashtf_surft(land_pts,nsurft)                                                  &
                             ! Coefficient to calculate surface
                             ! heat flux into soil (W/m2/K).
,lw_down_surftsum(land_pts)                                                    &
                             ! Gridbox sum of elevation corrections to
                             ! downward longwave radiation
,lw_down_surftabs(land_pts)
                             ! Gridbox sum of absolution changes to downward
                             ! longwave radiation from elevation corrections

! dummy arrays required for sea and se-ice to create universal
! routines for all surfaces
REAL(KIND=real_jlslsm) ::                                                      &
 array_zero(t_i_length * t_j_length)                                           &
                                ! Array of zeros
,zdt_dummy(t_i_length * t_j_length)
                                ! Dummy array for zdt

!Gridbox mean values calculated from soil tiled versions for FLAKE
REAL(KIND=real_jlslsm) ::                                                      &
  hcons_mean_soil(land_pts),                                                   &
  tsoil_mean_soil(land_pts)


!  Local scalars :-

INTEGER ::                                                                     &
 i,j                                                                           &
             ! Loop counter (horizontal field index).
,k                                                                             &
             ! Loop counter (tile field index).
,l                                                                             &
             ! Loop counter (land point field index).
,n                                                                             &
             ! Loop counter (tile index).
,nn                                                                            &
             ! Loop counter (soil carbon layers)
,m                                                                             &
             ! Index for soil tile
,n_veg
             ! Actual or dummy pointer to array
             ! defined only on PFTs

REAL(KIND=real_jlslsm) ::                                                      &
 ds_ratio                                                                      &
             ! 2 * snowdepth / depth of top soil layer.
,d_t                                                                           &
             ! Temporary in calculation of alpha1.
,zetam                                                                         &
             ! Temporary in calculation of CHN.
,zetah                                                                         &
             ! Temporary in calculation of CHN.
,zeta1                                                                         &
             ! Work space
,z0                                                                            &
             ! yet more workspace

             ! Temporary variables for adjustment of downwelling
             ! logwave to elevation tiles and correction back to
             ! conserve gridbox mean
,t_rad

LOGICAL ::                                                                     &
 l_vegdrag_active_here                                                         &
             ! Logical to indicate whether the vegetative drag scheme
             ! is active on the current surface tile in cases where
             ! l_vegdrag_surft itself may not be applicable
,l_shallow_lake_depth(land_pts)
             ! Logical to indicate unsuitable shallow lake depth

REAL(KIND=real_jlslsm) :: sea_point

INTEGER :: n_diag

INTEGER :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CABLE_LAND_SF_EXPLICIT'

!CABLE_LSM:CM2 
!CABLE TYPES containing field data (IN OUT)
TYPE(progs_cbl_vars_type), INTENT(IN OUT)  :: progs_cbl
TYPE(work_vars_type),      INTENT(IN OUT)  :: work_cbl
TYPE(params_io_data_type), INTENT(IN OUT)  :: pars_io_cbl
TYPE(progs_cnp_vars_type), INTENT(IN OUT)  :: progs_cnp

INTEGER :: mype, timestep_number
REAL :: satcon_soilt(land_pts, sm_levels) !0:sm_levels
REAL :: latitude(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL :: longitude(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL :: u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL :: ls_rain(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL :: ls_snow(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

  
logical, save :: first_call = .true.
!CABLE_LSM: End

!CM3#55 As of 6/12/2023 no attempt is being made to tidy up the workspace


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
 
!CABLE_LSM:initia;lize intent(OUT) - HOWEVER, most of these won't even be req^d in cable_land_expl AND others will get a value from
!CABLE 

IF( first_call ) THEN
  dtstar_surft(:,:)           =  0.0 
  resp_l_pft(:,:)             =  0.0                        
  resp_r_pft(:,:)             =  0.0                                 
  resp_w_pft(:,:)             =  0.0                                                                         
  n_leaf(:,:)                 =  0.0                           
  n_root(:,:)                 =  0.0                           
  n_stem(:,:)                 =  0.0                           
  lai_bal(:,:)                =  0.0                            
  wt_ext_surft(:,:,:)         =  0.0                                
  tile_frac(:,:)              =  0.0    
  fsmc_pft(:,:)               =  0.0                             
  gc_corr(:,:)                =  0.0                            
  isoprene_gb(:)              =  0.0    
  terpene_gb(:)               =  0.0    
  methanol_gb(:)              =  0.0    
  acetone_gb(:)               =  0.0    
  isoprene_pft(:,:)           =  0.0    
  terpene_pft(:,:)            =  0.0    
  methanol_pft(:,:)           =  0.0    
  acetone_pft(:,:)            =  0.0    
  fapar_diag_pft(:,:)         =  0.0    
  apar_diag_pft(:,:)          =  0.0    
  apar_diag_gb(:)             =  0.0    
  gpp_gb_acc(:)               =  0.0    
  gpp_pft_acc(:,:)            =  0.0    
  gs_irr_surft(:,:)           =  0.0    
  smc_irr_soilt(:,:)          =  0.0    
  wt_ext_irr_surft(:,:,:)     =  0.0    
  gc_irr_surft(:,:)           =  0.0    
  flux_o3_pft(:,:)            =  0.0    
  fo3_pft(:,:)                =  0.0    
  gpp(:)                      =  0.0    
  npp(:)                      =  0.0    
  smc_soilt(:,:)              =  0.0    
  gpp_pft(:,:)                =  0.0    
  resp_p(:)                   =  0.0    
  g_leaf(:,:)                 =  0.0    
  first_call = .FALSE. 
END IF

!-----------------------------------------------------------------------
!  0. Initialisations
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, l, n)                               &
!$OMP SHARED(t_i_length, t_j_length, array_zero, tdims, rhokm_land, cd_land,   &
!$OMP ch_land, zdt_dummy, land_pts, lake_ice_mid_temp, l_shallow_lake_depth,   &
!$OMP zdt_surft, db_surft, rhokh_can, nsurft,ftl_surft,fqw_surft,              &
!$OMP rib_surft,z0m_surft,u_s_std_surft,chr1p5m,resfs,alpha1,radnet_surft,     &
!$OMP scaling_urban,snowdep_surft,snowdepth_surft,sf_diag)
!$OMP DO SCHEDULE(STATIC)
DO i = 1,t_i_length * t_j_length
  array_zero(i)     = 0.0
  zdt_dummy(i)      = 0.0
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    rhokm_land(i,j) = 0.0
    cd_land(i,j)    = 0.0
    ch_land(i,j)    = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO l = 1,land_pts
  lake_ice_mid_temp(l)    = 0.0
  l_shallow_lake_depth(l) = .FALSE.
END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
!  1. Initialise FTL_SURFT and RIB_SURFT on all tiles at all points,
!     to allow STASH to process these as diagnostics.
!-----------------------------------------------------------------------
DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    ! MORUSES Initialise urban roughness array
    ftl_surft(l,n)     = 0.0
    fqw_surft(l,n)     = 0.0
    rib_surft(l,n)     = 0.0
    z0m_surft(l,n)     = 0.0    !jh:can't 0 this AND cancel roughness
    u_s_std_surft(l,n) = 0.0
    chr1p5m(l,n)       = 0.0
    resfs(l,n)         = 0.0
    alpha1(l,n)        = 0.0
    !CN2ish!radnet_surft(l,n)  = 0.0
    scaling_urban(l,n) = 1.0
    zdt_surft(l,n)     = 0.0
    db_surft(l,n)      = 0.0
    rhokh_can(l,n)     = 0.0
    ! Equivalent snowdepth for surface calculations.
    snowdep_surft(l,n) = snowdepth_surft(l,n)

    IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
      sf_diag%resfs_stom(l,n) = 0.0
    END IF
  END DO
!$OMP END DO NOWAIT
END DO
IF (sf_diag%l_tau_surft) THEN
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,land_pts
      sf_diag%tau_surft(l,n) = 0.0
    END DO
!$OMP END DO NOWAIT
  END DO
END IF
!$OMP END PARALLEL

!-----------------------------------------------------------------------
! Call TILEPTS to calculate surft_pts and surft_index for surface types
!-----------------------------------------------------------------------
CALL tilepts(land_pts, frac, surft_pts, surft_index, l_lice_point)

!CABLE_LSM: following CM2 - dodge OMP here
DO n = 1,nsurft
  DO l = 1,land_pts
    canhc_surft(l,n) = 0.0
    vfrac_surft(l,n) = 0.0
    IF( frac(l,n) == 0.0 ) THEN
      cd_surft(l,n) = 0.0
      ch_surft(l,n) = 0.0
      z0h_surft(l,n) = 0.0
      z0m_eff_surft(l,n) = 0.0
      !not In 7.1!rhokpm(l,n) = 0.0 !doesnt go anywhere
      radnet_surft(l,n) = 0.0
      fraca(l,n) = 0.0
      resfs(l,n) = 0.0
      resft(l,n) = 0.0
    END IF
  END DO
END DO
!CABLE_LSM: End 

!CM3#55-1
!!-----------------------------------------------------------------------
!!   Generate the anthropogenic heat for surface calculations
!!-----------------------------------------------------------------------
!IF ( l_anthrop_heat_src .AND. .NOT. l_aggregate ) THEN
!  CALL generate_anthropogenic_heat( curr_day_number, land_pts, frac,           &
!                                  surft_pts, surft_index,                      &
!                                  !New arguments replacing USE statements
!                                 !urban_param (IN)
!                                  wrr_gb,                                      &
!                                  !Fluxes (IN OUT)
!                                  anthrop_heat_surft)
!END IF

!CABLE_LSM:CM2! Even if this is left in, crashes immediately in call to alb_pft..
!CM2! Even if this is left in, crashes immediately in call to alb_pft..
!CM2!!-----------------------------------------------------------------------
!CM2!! Call physiology routine to calculate surface conductances and carbon
!CM2!! fluxes.
!CM2!!-----------------------------------------------------------------------
!CM2!CALL physiol (                                                                 &
!CM2!  land_pts,land_index,                                                         &
!CM2!  sm_levels,nsurft,surft_pts,surft_index,                                      &
!CM2!  dim_cs1,                                                                     &
!CM2!  co2_mmr,co2_3d,co2_dim_len, co2_dim_row,l_co2_interactive,                   &
!CM2!  can_model,cs_pool_soilt,veg_state,frac,canht_pft,photosynth_act_rad,         &
!CM2!  lai_pft,pstar,qw_1,sthu_soilt,sthf_soilt,t_soil_soilt,tstar_surft,           &
!CM2!  smvccl_soilt,smvcst_soilt,smvcwt_soilt,vshr,z0_surft,z1_uv,o3,               &
!CM2!  canhc_surft,vfrac_surft,emis_surft,l_emis_surft_set,emis_soil,flake,         &
!CM2!  g_leaf,gs,gc_surft,gc_stom_surft,gc_corr,gpp,gpp_pft,npp,npp_pft,            &
!CM2!  resp_p,resp_p_pft,resp_s_soilt,resp_l_pft,                                   &
!CM2!  resp_r_pft,resp_w_pft,n_leaf,                                                &
!CM2!  n_root,n_stem,lai_bal,                                                       &
!CM2!  smc_soilt,wt_ext_surft,fsmc_pft,                                             &
!CM2!  albsoil_soilt,cos_zenith_angle,                                              &
!CM2!  can_rad_mod,ilayers,flux_o3_pft,fo3_pft,sf_diag,asteps_since_triffid,        &
!CM2!  non_lake_frac,                                                               &
!CM2!  !New arguments replacing USE statements
!CM2!  !Fluxes (IN)
!CM2!  t_home_gb,t_growth_gb,                                                       &
!CM2!  !urban_param (IN)
!CM2!  emisr_gb, emisw_gb, hwr_gb,                                                  &
!CM2!  !jules_mod (IN OUT)
!CM2!  albobs_scaling_surft,                                                        &
!CM2!  !jules_chemvars_mod (OUT)
!CM2!  isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                         &
!CM2!  methanol_gb, methanol_pft, acetone_gb, acetone_pft,                          &
!CM2!  !trif_vars_mod (OUT)
!CM2!  fapar_diag_pft, apar_diag_pft, apar_diag_gb, gpp_gb_acc, gpp_pft_acc,        &
!CM2!  !crop_vars_mod (IN)
!CM2!  rootc_cpft, sthu_irr_soilt, frac_irr_soilt, frac_irr_surft, dvi_cpft,        &
!CM2!  !crop_vars_mod (OUT)
!CM2!  gs_irr_surft, smc_irr_soilt, wt_ext_irr_surft, gc_irr_surft,                 &
!CM2!  !p_s_parms (IN)
!CM2!  bexp_soilt, sathh_soilt, v_close_pft, v_open_pft,                            &
!CM2!  !ancil_info
!CM2!  l_soil_point,                                                                &
!CM2!  !jules_surface_types (IN)
!CM2!  diff_frac)

!CM3#55-2
!! Update gc_surft for canopy snow if using the canopy snow scheme
!IF ( .NOT. l_aggregate .AND. can_model == 4) THEN
!  DO n = 1,npft
!    IF ( cansnowtile(n) ) THEN
!!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)       &
!!$OMP          SHARED(surft_pts, surft_index, land_index, t_i_length,          &
!!$OMP                 snow_surft, gc_surft, catch_snow, tstar_surft,           &
!!$OMP                 vshr_land, n)   SCHEDULE(STATIC)
!      DO k = 1,surft_pts(n)
!        l = surft_index(k,n)
!        IF (snow_surft(l,n) >  0.0) THEN
!          j = (land_index(l) - 1) / t_i_length + 1
!          i = land_index(l) - (j-1) * t_i_length
!          gc_surft(l,n) = 0.06 * snow_surft(l,n)**0.6 * catch_snow(l,n)**0.4   &
!                       * 2.06e-5 * (tm / tstar_surft(l,n))**1.75               &
!                       * (1.79+3 * SQRT(vshr_land(i,j)))                       &
!                       / (2 * rho_ice * 5.0e-4**2)
!        END IF
!      END DO
!!$OMP END PARALLEL DO
!    END IF
!  END DO
!END IF

!CM3#55-3 - note that at some point CASA variables need to make their way
!           into these TRIFFID vars for output purposes
!!----------------------------------------------------------------------
!! If TRIFFID is being used apply any correction to the land-atmosphere
!! fluxes on the first timestep after the last TRIFFID call. Such a
!! correction will typically be associated with a total depletion of
!! carbon or with maintanence of the seed fraction. The corrections
!! are stored in the accumulation variables after the call to TRIFFID.
!! The correction is added to the instantaneous land-atmosphere fluxes
!! (so that the atmospheric carbon budget is corrected) but is not
!! included in the accumulation variables which drive TRIFFID, since
!! this has already been dealt with during the last TRIFFID call.
!!----------------------------------------------------------------------
!IF (l_triffid .AND. (asteps_since_triffid == 1)                                &
!    .AND. ( cycleno == numcycles .OR. l_quick_ap2) ) THEN
!!jhan
!  DO n = 1,nnpft
!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                    &
!!$OMP SHARED(land_pts, npp_pft, npp_pft_acc, timestep, resp_p_pft, n)
!    DO l = 1,land_pts
!      npp_pft(l,n) = npp_pft(l,n) + npp_pft_acc(l,n) / timestep
!      resp_p_pft(l,n) = resp_p_pft(l,n) - npp_pft_acc(l,n) / timestep
!      npp_pft_acc(l,n)=-npp_pft_acc(l,n)
!    END DO
!!$OMP END PARALLEL DO
!  END DO
!
!  ! Here we have assumed that RothC must be used with TRIFFID, and is called
!  ! on the same timestep.
!  IF ( soil_bgc_model == soil_model_rothc ) THEN
!    DO n = 1,dim_cs1
!      DO nn = 1,dim_cslayer
!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                    &
!!$OMP SHARED(land_pts, resp_s_soilt, resp_s_acc_soilt, timestep, n, nn)
!        DO l = 1,land_pts
!          !soil tiling is not compatible with triffid. OK to hard-code soil
!          !tile index to 1 here
!          resp_s_soilt(l,1,nn,n) = resp_s_soilt(l,1,nn,n)                      &
!                                   + (resp_s_acc_soilt(l,1,nn,n) / timestep)
!          resp_s_acc_soilt(l,1,nn,n) = -resp_s_acc_soilt(l,1,nn,n)
!        END DO
!!$OMP END PARALLEL DO
!      END DO
!    END DO
!  END IF
!
!END IF

!CM3#55-4 - note that at some point we need to ensure that CABLE values
!           into these JULES vars for output purposes
!!----------------------------------------------------------------------
!! Increment accumulation of leaf turnover rate.
!! This is required for leaf phenology and/or TRIFFID, either of
!! which can be enabled independently of the other.
!!----------------------------------------------------------------------
!IF ( cycleno == numcycles .OR. l_quick_ap2 ) THEN
!
!  IF (l_phenol .AND. .NOT. l_triffid) THEN
!    DO n = 1,nnpft
!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                    &
!!$OMP SHARED(n, land_pts, g_leaf_acc, g_leaf, timestep)
!      DO l = 1,land_pts
!        g_leaf_acc(l,n) = g_leaf_acc(l,n) +                                    &
!                          g_leaf(l,n) * ( timestep / secs_per_360days )
!      END DO
!!$OMP END PARALLEL DO
!    END DO
!  END IF

!!jhan
!  !----------------------------------------------------------------------
!  ! Increment accumulation prognostics for TRIFFID
!  !----------------------------------------------------------------------
!  IF (l_triffid) THEN
!    DO n = 1,nnpft
!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                    &
!!$OMP SHARED(n, land_pts, g_leaf_acc, g_leaf, timestep,                        &
!!$OMP        npp_pft_acc, npp_pft, resp_w_pft_acc, resp_w_pft)
!      DO l = 1,land_pts
!        g_leaf_acc(l,n) = g_leaf_acc(l,n) +                                    &
!                          g_leaf(l,n) * ( timestep / secs_per_360days )
!        npp_pft_acc(l,n) = npp_pft_acc(l,n) + npp_pft(l,n) * timestep
!        resp_w_pft_acc(l,n) = resp_w_pft_acc(l,n)                              &
!                              + resp_w_pft(l,n) * timestep
!      END DO
!!$OMP END PARALLEL DO
!    END DO
!  END IF
!
!IF ( soil_bgc_model == soil_model_rothc ) THEN
!    DO n = 1,dim_cs1
!      DO nn = 1,dim_cslayer
!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                    &
!!$OMP SHARED(n, nn, land_pts,resp_s_soilt, resp_s_acc_soilt, timestep)
!        DO l = 1,land_pts
!          !soil tiling is not compatible with triffid. OK to hard-code soil
!          !tile index to 1 here
!          ! calculated before resp_frac applied
!          resp_s_acc_soilt(l,1,nn,n) = resp_s_acc_soilt(l,1,nn,n)              &
!                                       + resp_s_soilt(l,1,nn,n) * timestep
!        END DO
!!$OMP END PARALLEL DO
!      END DO
!    END DO
!  END IF
!
!END IF ! CycleNo == NumCycles

!CM3#55-5
!!-----------------------------------------------------------------------
!! calculate CO2:(BIO+HUM) ratio, dependent on soil clay content, and
!! sum soil respiration components
!! (RESP_FRAC here then contains the fraction of soil respiration which
!! is respired to the atmos. the rest is re-partitioned into BIO+HUM)
!
!! resp_s_acc_soilt contains the full amount, and this is carried forward to
!! VEG_CTL for use in updating soil carbon pools. RESP_S_TOT calculated
!! here is passed to BL_TRMIX as the fraction which is respired as CO2
!! to the atmosphere. RESP_S_TOT, and RESP_S are also passed out for
!! storage in diagnostics 3293, and 3467-470.
!
!!-----------------------------------------------------------------------
!IF ( soil_bgc_model == soil_model_rothc ) THEN
!  DO j = 1, nsoilt
!    DO i = 1, land_pts
!      resp_s_tot_soilt(i,j)=0.0
!    END DO
!  END DO
!  !soil tiling is not compatible with triffid. OK to hard-code soil tile
!  !index to 1 here by setting m = 1
!  m = 1
!  DO nn = 1,dim_cslayer
!    DO i = 1,land_pts
!      work_clay = EXP(-0.0786 * 100.0 * clay_soilt(i,m,nn))
!      resp_frac(i,nn) = (3.0895+2.672 * work_clay) /                           &
!                        (4.0895+2.672 * work_clay)
!      resp_s_soilt(i,m,nn,1)  = resp_s_soilt(i,m,nn,1) * resp_frac(i,nn)
!      resp_s_soilt(i,m,nn,2)  = resp_s_soilt(i,m,nn,2) * resp_frac(i,nn)
!      resp_s_soilt(i,m,nn,3)  = resp_s_soilt(i,m,nn,3) * resp_frac(i,nn)
!      resp_s_soilt(i,m,nn,4)  = resp_s_soilt(i,m,nn,4) * resp_frac(i,nn)
!      resp_s_tot_soilt(i,m)   = resp_s_tot_soilt(i,m)                          &
!                                + resp_s_soilt(i,m,nn,1)                       &
!                                + resp_s_soilt(i,m,nn,2)                       &
!                                + resp_s_soilt(i,m,nn,3)                       &
!                                + resp_s_soilt(i,m,nn,4)
!    END DO  !  layers
!  END DO  !  points
!END IF

!CM3#55-6 - CABLE always runs with tiles
!!-----------------------------------------------------------------------
!! Reset surft_pts and surft_index and set tile fractions to 1 if aggregate
!! tiles are used (L_AGGREGATE=.T.).
!! Otherwise, set tile fractions to surface type fractions.
!!-----------------------------------------------------------------------
!IF (l_aggregate) THEN
!  surft_pts(1) = land_pts
!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                    &
!!$OMP SHARED(land_pts, tile_frac, surft_index)
!  DO l = 1,land_pts
!    tile_frac(l,1) = 1.0
!    surft_index(l,1) = l
!  END DO
!!$OMP END PARALLEL DO
!ELSE
  DO n = 1,ntype
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                    &
!$OMP SHARED(land_pts, tile_frac, frac, n)
    DO l = 1, land_pts
      tile_frac(l,n) = frac(l,n)
    END DO
!$OMP END PARALLEL DO
  END DO
!END IF

IF (land_pts >  0) THEN    ! Omit if no land points

   !CM3#55-7 - CABLE needs to provide values for hcons_soilt and hcons_snow
   !           retain for now but ideally we'd use CABLE science not JULES science 
  !-----------------------------------------------------------------------
  ! Calculate the thermal conductivity of the top soil layer.
  !-----------------------------------------------------------------------
  DO m = 1, nsoilt
    CALL heat_con (land_pts,hcon_soilt,sthu_soilt(:,m,1),                      &
                   sthf_soilt(:,m,1),smvcst_soilt(:,m,1),hcons_soilt(:,m))
  END DO

  ! Thermal conductvity of top snow layer if nsmax > 0
  !CM2!IF (nsmax > 0) THEN
  !CM2!  DO n = 1,nsurft
  !CM2!    CALL snowtherm(land_pts,surft_pts(n),nsnow_surft(:,n),                   &
  !CM2!                   surft_index(:,n),ds_surft(:,n,:),sice_surft(:,n,:),       &
  !CM2!                   sliq_surft(:,n,:),csnow,ksnow)
  !CM2!    DO l = 1,land_pts
  !CM2!      hcons_snow(l,n) = ksnow(l,1)
  !CM2!    END DO
  !CM2!  END DO
  !CM2!END IF

END IF                     ! End test on land points

!CM3#55-8 - correctly removed already
!-----------------------------------------------------------------------
! Calculate net radiation on land tiles
!-----------------------------------------------------------------------
!CM2!$OMP PARALLEL                                                                 &
!CM2!$OMP DEFAULT(NONE)                                                            &
!CM2!$OMP PRIVATE(l,k,j,i,n)                                                       &
!CM2!$OMP SHARED(surft_pts,surft_index,land_index,pdims,radnet_surft,sw_surft,     &
!CM2!$OMP        emis_surft,sky,lw_down,tstar_surft,nsurft,l_skyview)
!CABLE_LSM:CM2
!CM2IF (l_skyview) THEN
!CM2  DO n = 1,nsurft
!CM2!$OMP DO SCHEDULE(STATIC)
!CM2    DO k = 1,surft_pts(n)
!CM2      l = surft_index(k,n)
!CM2      j=(land_index(l) - 1) / pdims%i_end + 1
!CM2      i = land_index(l) - (j-1) * pdims%i_end
!CM2      radnet_surft(l,n) = sw_surft(l,n) + emis_surft(l,n) *                    &
!CM2        sky(i,j) * ( lw_down(i,j) - sbcon * tstar_surft(l,n)**4 )
!CM2    END DO
!CM2!$OMP END DO NOWAIT
!CM2  END DO
!CM2ELSE
!CM2  DO n = 1,nsurft
!CM2!$OMP DO SCHEDULE(STATIC)
!CM2    DO k = 1,surft_pts(n)
!CM2      l = surft_index(k,n)
!CM2      j=(land_index(l) - 1) / pdims%i_end + 1
!CM2      i = land_index(l) - (j-1) * pdims%i_end
!CM2      radnet_surft(l,n) = sw_surft(l,n) + emis_surft(l,n) *                    &
!CM2                 ( lw_down(i,j) - sbcon * tstar_surft(l,n)**4 )
!CM2    END DO
!CM2!$OMP END DO NOWAIT
!CM2  END DO
!CM2END IF
!CM2!!$OMP END PARALLEL


!-----------------------------------------------------------------------
! 4.  Surface turbulent exchange coefficients and "explicit" fluxes
!     (P243a, routine SF_EXCH).
!     Wind mixing "power" and some values required for other, later,
!     diagnostic calculations, are also evaluated if requested.
!-----------------------------------------------------------------------

!CM3#55-9
!IF ( l_rp2 .AND. i_rp_scheme == i_rp2b) THEN
!  DO n = 1,npft
!    z0h_z0m(n) = z0hm_pft_rp(n)
!  END DO
!END IF

!CM3#55-10 - note that snowdep_surft needs to take CABLE value - check

!$OMP PARALLEL PRIVATE(i,j,l,n) DEFAULT(NONE) IF(land_pts>1)                   &
!$OMP SHARED(land_pts,nsurft,land_index, flandg,z_land,t_i_length,             &
!$OMP can_model,cansnowtile,l_snowdep_surf,snow_surft,rho_snow_const,          &
!$OMP snowdep_surft,surf_hgt_surft,ctile_orog_fix,l_ctile)
!DO n = 1,nsurft
!  IF ( (can_model == 4) .AND. cansnowtile(n) .AND. l_snowdep_surf ) THEN
!!$OMP DO SCHEDULE(STATIC)
!    DO l = 1,land_pts
!      snowdep_surft(l,n) = snow_surft(l,n) / rho_snow_const
!    END DO
!!$OMP END DO NOWAIT
!  END IF
!END DO

!CM3#55-11
IF (ctile_orog_fix == correct_sea_adjust_land .AND. l_ctile) THEN
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    IF (flandg(i,j) > 0.0 .AND. flandg(i,j) < 1.0) THEN
      ! calculate height of orography relative to grid-box mean
      ! limit this to 1000m. z_land already limited to be > 0
      ! in ni_bl_ctl
      surf_hgt_surft(l,:) = MIN(z_land(i,j) * (1.0 / flandg(i,j) - 1.0),1000.0)
    END IF
  END DO
!$OMP END DO NOWAIT
END IF
!$OMP END PARALLEL

!CM3#55-12
! Calculate temperature and specific humidity for elevation bands
CALL elevate(                                                                  &
 land_pts,nsurft,surft_pts,land_index,surft_index,                             &
 tl_1,qw_1,pstar,surf_hgt_surft,l_elev_absolute_height,z_land,                 &
 t_elev,q_elev)


!-----------------------------------------------------------------------
!  2.  Calculate QSAT values required later.
!-----------------------------------------------------------------------
!CM3#55-13

!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(i, j, l)              &
!$OMP SHARED(land_index, land_pts, pstar, pstar_land, t_i_length)              &
!$OMP SCHEDULE(STATIC)
DO l = 1,land_pts
  j=(land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length
  pstar_land(l) = pstar(i,j)
END DO
!$OMP END PARALLEL DO

IF (l_mr_physics) THEN
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(n)                      &
!$OMP SHARED(nsurft,qstar_surft,tstar_surft,pstar_land,land_pts,               &
!$OMP qs1_elev,t_elev)  SCHEDULE(STATIC)
  DO n = 1,nsurft
    CALL qsat_mix(qstar_surft(:,n),tstar_surft(:,n),pstar_land,land_pts)
    CALL qsat_mix(qs1_elev(:,n),t_elev(:,n),pstar_land,land_pts)
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(n)                      &
!$OMP SHARED(nsurft,qstar_surft,tstar_surft,pstar_land,land_pts,               &
!$OMP qs1_elev,t_elev)  SCHEDULE(STATIC)
  DO n = 1,nsurft
    CALL qsat(qstar_surft(:,n),tstar_surft(:,n),pstar_land,land_pts)
    CALL qsat(qs1_elev(:,n),t_elev(:,n),pstar_land,land_pts)
  END DO
!$OMP END PARALLEL DO
END IF

!CM3#55-14
!-----------------------------------------------------------------------
!!  Calculate gradient of saturated specific humidity for use in
!!  calculation of surface fluxes
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,l,n,d_t, i, j)                          &
!$OMP SHARED(nsurft,surft_pts,surft_index,tstar_surft,t_elev,alpha1,           &
!$OMP        qstar_surft,qs1_elev,epsil,c_virtual,r, tdims,rhostar_mom,rhostar)
DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    d_t = tstar_surft(l,n) - t_elev(l,n)
    IF (d_t > 0.05 .OR. d_t < -0.05) THEN
      alpha1(l,n) = (qstar_surft(l,n) - qs1_elev(l,n)) / d_t
    ELSE IF (t_elev(l,n) > tm) THEN
      alpha1(l,n) = epsil * lc * qs1_elev(l,n) *                               &
                  (1.0 + c_virtual * qs1_elev(l,n)) /                          &
                  ( r * t_elev(l,n) * t_elev(l,n))
    ELSE
      alpha1(l,n) = epsil * ls * qs1_elev(l,n) *                               &
                  (1.0 + c_virtual * qs1_elev(l,n)) /                          &
                  ( r * t_elev(l,n) * t_elev(l,n))
    END IF
  END DO
!$OMP END DO NOWAIT
END DO

!-----------------------------------------------------------------------
! If requested, improve accuracy of air density, rhostar
! On input rhostar = pstar/R*Tstar
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    ! original approximation for surface air density
    rhostar_mom(i,j) = rhostar(i,j)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (l_accurate_rho) THEN
  ! More accurate expressions for surface air density.
  ! Use bottom level vapour as a better approximation over land
  ! than qsat(tstar)
  CALL calc_air_dens(l_mr_physics,qw_1,rhostar,rhostar_mom)

END IF ! l_accurate_rho

!CM3#55-16
!! Initialise scaling_urban to 1.0 so that it only affects urban tiles when
!! MORUSES used with no aggregation.
!IF ( .NOT. l_aggregate .AND. l_moruses_storage ) THEN
!  n = urban_canyon
!!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(land_pts,   &
!!$OMP             n, scaling_urban, tile_frac, urban_roof) SCHEDULE(STATIC)
!  DO l = 1, land_pts
!    IF ( tile_frac(l,n) > 0.0 ) THEN
!      scaling_urban(l,n) =                                                     &
!         ( tile_frac(l,n) + tile_frac(l,urban_roof) ) /                        &
!         tile_frac(l,n)
!    END IF
!  END DO
!!$OMP END PARALLEL DO
!END IF

!CM3#55-17
!! Calculate average layer temperature and conductivity for lakes.
!! This is a fudge - a layer with average properties won't
!! really behave like a stack of layers with different properties.
!!
!IF (     l_flake_model                                                         &
!    .AND. ( .NOT. l_aggregate)) THEN
!
!  !==============================================================================
!  ! *NOTICE REGARDING SOIL TILING**
!  !
!  !The following section facilitates the use of soil tiling. As implemented,
!  !there are two soil tiling options:
!  !
!  !nsoilt == 1
!  !Operate as with a single soil tile, functionally identical to JULES upto
!  ! at least vn4.7 (Oct 2016)
!  ! This means that a soilt variable being passed 'up' to the surface is
!  ! broadcast to the surft variable (with weighting by frac if requred)
!  !
!  !nsoilt > 1
!  !Operate with nsoilt = nsurft, with a direct mapping between them
!  ! This means that a soilt variable being passed 'up' to the surface is simply
!  ! copied into the surft variable
!  !
!  ! This will need to be refactored for other tiling approaches. This note
!  ! will be replicated elsewhere in the code as required
!  !
!  !These comments apply until **END NOTICE REGARDING SOIL TILING**
!  !==============================================================================
!
!!$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, m)                                     &
!!$OMP SHARED(l_flake_model, l_aggregate, land_pts, nsoilt, hcon_lake,          &
!!$OMP        dzsoil, hcons_mean_soil, hcons_soilt, t_soil_soilt, tile_frac,    &
!!$OMP        lake_h_ice_gb, snow_hcon, lake_depth_gb, nusselt_gb,              &
!!$OMP        g_dt_gb, ts1_lake_gb, lake_t_mxl_gb, lake_t_ice_gb,               &
!!$OMP        tsoil_mean_soil, lake_ice_mid_temp, l_shallow_lake_depth)
!
!  ! Initialise thermal variables to zero and then set to mean soil values
!!$OMP DO SCHEDULE(STATIC)
!  DO l = 1, land_pts
!    hcons_mean_soil(l) = 0.0
!    tsoil_mean_soil(l) = 0.0
!    hcon_lake(l)       = 0.0
!    ts1_lake_gb(l)     = 0.0
!  END DO
!!$OMP END DO
!  IF (nsoilt == 1) THEN
!    !Just 1 soil tile
!    m = 1
!!$OMP DO SCHEDULE(STATIC)
!    DO l = 1, land_pts
!      hcons_mean_soil(l) = hcons_soilt(l,m)
!      tsoil_mean_soil(l) = t_soil_soilt(l,m,1)
!    END DO
!!$OMP END DO
!  ELSE
!    !Surface tiles map directly on to soil tiles
!!$OMP DO SCHEDULE(STATIC)
!    DO m = 1,nsoilt
!      DO l = 1, land_pts
!        hcons_mean_soil(l) = hcons_mean_soil(l)                                &
!                            + (tile_frac(l,m) * hcons_soilt(l,m))
!        tsoil_mean_soil(l) = tsoil_mean_soil(l)                                &
!                            + (tile_frac(l,m) * t_soil_soilt(l,m,1))
!      END DO
!    END DO
!!$OMP END DO
!  END IF
!
!  !==============================================================================
!  ! *END NOTICE REGARDING SOIL TILING**
!  !==============================================================================
!
!  IF (dzsoil(1) <= 0.0) THEN
!    !
!    ! catch-all for sillies - just use soil value
!    !
!!$OMP DO SCHEDULE(STATIC)
!    DO l = 1, land_pts
!      hcon_lake(l)   = hcons_mean_soil(l)
!      ts1_lake_gb(l) = tsoil_mean_soil(l)
!    END DO
!!$OMP END DO NOWAIT
!
!  ELSE
!
!!$OMP DO SCHEDULE(STATIC)
!    DO l = 1, land_pts
!
!      IF ( dzsoil(1) <= lake_h_ice_gb(l) ) THEN
!
!        ! Near surface layer is entirely within the ice. Set conductivity
!        ! to ice value and use piecewise linear interpolation to find
!        ! temperature at midpoint of near surface layer.
!        hcon_lake(l)   = hcice
!        ts1_lake_gb(l) = lake_t_ice_gb(l) + (dzsoil(1) / 2.0) *                &
!                         (lake_t_mxl_gb(l) - lake_t_ice_gb(l))                 &
!                         / lake_h_ice_gb(l)
!
!      ELSE IF (      (dzsoil(1)  >  lake_h_ice_gb(l))                          &
!               .AND. (dzsoil(1) <= (lake_h_ice_gb(l)                           &
!                               +lake_depth_gb( l)))) THEN
!
!        ! Set conductivity as weighted average of ice and water based upon
!        ! realtive thicknesses and the Nusselt number
!        nusselt_gb(l) = g_dt_gb(l) * (dzsoil(1) - lake_h_ice_gb(l) )           &
!                     / ( 2.0 * hcwat )
!        nusselt_gb(l) = MAX( nusselt_gb(l), 1.0 )
!        hcon_lake(l)  = ( hcice * lake_h_ice_gb(l)                             &
!                         +hcwat * (dzsoil(1) - lake_h_ice_gb(l)))              &
!                                * nusselt_gb(l)                                &
!                       / dzsoil(1)
!
!        ! Find temperature mid-way through the ice layer
!        ! (temperature varies linearly across ice layer)
!        lake_ice_mid_temp(l) = lake_t_ice_gb(l) + 0.5 * (lake_t_mxl_gb(l)      &
!                            - lake_t_ice_gb(l))
!
!        ! Calculate temperature as weighted average of temperature in ice layer
!        ! and mixed layer temp for near surface water beneath the ice.
!        ts1_lake_gb(l) = (lake_ice_mid_temp(l) * lake_h_ice_gb(l)              &
!                          + lake_t_mxl_gb(l) * (dzsoil(1) - lake_h_ice_gb(l))) &
!                         / dzsoil(1)
!
!      ELSE
!
!        ! Lake depth is less than first soil layer thickness so set logical
!        ! to write our warning message.
!        l_shallow_lake_depth(l) = .TRUE.
!
!        ! Set conductivity as weighted average of ice and water and soil
!        ! based upon realtive thicknesses and the Nusselt number
!        ! Use soil value for the temperature
!        nusselt_gb(l)        = g_dt_gb(l) * lake_depth_gb(l)                   &
!                               / ( 2.0 * hcwat )
!        nusselt_gb(l)        = MAX( nusselt_gb(l), 1.0 )
!        hcon_lake(l)         = ( hcice * lake_h_ice_gb(l)                      &
!                               + hcwat * lake_depth_gb(l) * nusselt_gb(l)      &
!                               + hcons_mean_soil(l) * (dzsoil(1)               &
!                                 - lake_h_ice_gb(l)  - lake_depth_gb(l)))      &
!                               / dzsoil(1)
!        ts1_lake_gb(l)       = tsoil_mean_soil(l)
!
!      END IF
!    END DO
!!$OMP END DO NOWAIT
!  END IF
!!$OMP END PARALLEL
!
!  ! Write out warning message for negative top soil layer thinkness
!  ! or lake depth less than first soil layer thinkness
!  IF (dzsoil(1) <= 0.0) THEN
!    errcode = -1
!    WRITE(cmessage, '(A,F16.4)') 'Negative value of dzsoil = ', dzsoil(1)
!    CALL ereport(routinename, errcode, cmessage)
!  END IF
!
!  IF ( ANY(l_shallow_lake_depth(1:land_pts)) ) THEN
!    errcode = -1
!    WRITE(cmessage, '(A,F16.4)') 'Unusual value of lake depth found '//        &
!                                 'when computing hcon_lake. Found '  //        &
!                                 'lake depth comparable to first ' //          &
!                                 'soil layer: dzsoil= ', dzsoil(1)
!    CALL ereport(routinename, errcode, cmessage)
!  END IF
!
!END IF

!CM3#55-18
DO n = 1,nsurft
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                    &
!$OMP SHARED(n, land_pts, lw_down_elevcorr_surft)
  DO l = 1,land_pts
    lw_down_elevcorr_surft(l,n) = 0.0
  END DO
!$OMP END PARALLEL DO
END DO

!CM3#55-19 - moving the necessary parts to after the call to CABLE
!IF( first_call ) THEN
!  emis_surft = 1.0 
!  emis_soil  = 1.0 
!  hcons_snow = hcons_soilt(1,1)
!  first_call = .FALSE.
!END IF
!hcons_surf = hcons_soilt(1,1)
!hcons_snow = hcons_soilt(1,1)

!CM3#55-20
IF (l_elev_lw_down) THEN

  ! for tiles at different elevations, adjust downwelling longwave
  ! according to the anomaly in surface temperature that has been
  ! calculated with LW ~ T^4
  DO l = 1,land_pts
    lw_down_surftsum(l) = 0.0
    lw_down_surftabs(l) = 0.0
  END DO

  DO n = 1,nsurft
    DO l = 1,land_pts
      j = (land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      IF (lw_down(i,j) > 0.0) THEN

        ! Adjust radiative temperature and net longwave
        t_rad = 0.0
        IF (t_elev(l,n) > 0.0 ) THEN
          t_rad = (lw_down(i,j) / sbcon)**(1.0 / 4.0)
          t_rad = t_rad + t_elev(l,n) - tl_1(i,j)
          lw_down_elevcorr_surft(l,n) =  sbcon * (t_rad**4) - lw_down(i,j)
        END IF

        ! Keep track of total adjustments to longwave
        lw_down_surftsum(l) = lw_down_surftsum(l) +                            &
                              lw_down_elevcorr_surft(l,n) * tile_frac(l,n)
        lw_down_surftabs(l) = lw_down_surftabs(l) +                            &
                          ABS(lw_down_elevcorr_surft(l,n)) * tile_frac(l,n)
      END IF
    END DO
  END DO

  DO l = 1,land_pts
    j = (land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    IF (lw_down(i,j) > 0.0) THEN
      IF (lw_down_surftabs(l) > EPSILON(0.0)) THEN
        ! correct each adjustment to preserve the gridbox mean.
        ! size of correction in proportion to the size of adjustment
        ! so that unadjusted tiles remain unaffected
        DO n = 1,nsurft
          lw_down_elevcorr_surft(l,n) = lw_down_elevcorr_surft(l,n)            &
                                        - lw_down_surftsum(l)                  &
                                        * ABS(lw_down_elevcorr_surft(l,n))     &
                                        / lw_down_surftabs(l)
        END DO
      END IF
    END IF
  END DO

  !CM3#55-20 - these lines are correctly removed
!CM2!  ! Adjust the net radiation
!CM2!  DO n = 1,nsurft
!CM2!    DO l = 1,land_pts
!CM2!      j = (land_index(l) - 1) / t_i_length + 1
!CM2!      i = land_index(l) - (j-1) * t_i_length
!CM2!      IF (lw_down(i,j) > 0.0) THEN
!CM2!        IF (l_skyview) THEN
!CM2!          radnet_surft(l,n) = radnet_surft(l,n)                                &
!CM2!                              + sky(i,j) * emis_surft(l,n)                     &
!CM2!                                * lw_down_elevcorr_surft(l,n)
!CM2!        ELSE
!CM2!          radnet_surft(l,n) = radnet_surft(l,n)                                &
!CM2!                              + emis_surft(l,n) * lw_down_elevcorr_surft(l,n)
!CM2!        END IF
!CM2!      END IF
!CM2!    END DO
!CM2!  END DO
END IF

!=============================================================================
!CM3#55 - at this point all of the forcing for CABLE has been created -
!         so call CABLE and then map remaining values back to JULES variables
!         where needed.

!CABLE_LSM:CM2.5 Oddly tstar_surft is declared as INTENT(IN) above!!
! This is important may need to be followed up as the INTENTs are likely inherited
! - the simplest work around would be to not UNPACK %trad to tstar_surft
!   in the explicit section.
!
! However we do NEED to provide a value for dtstar_surft - this is missing
!
! Also we may wish to force cable with lw_down_elevcorr_surft not lw_down
CALL cable_explicit_main(                                                      &
      ! IN: UM/JULES model/grid parameters, fields, mappings
      mype, timestep, timestep_number, tdims%i_end, tdims%j_end, land_pts,     &
      nsurft, npft, sm_levels, dzsoil, land_index, surft_pts, surft_index,     &
      cos_zenith_angle, latitude, longitude, Flandg, tile_frac,                &
      
      ! IN: soil parameters !1 is only allowable index in UM
      bexp_soilt(:,1,:), hcon_soilt(:,1), satcon_soilt(:,:),                   &
      sathh_soilt(:,1,:), smvcst_soilt(:,1,:), smvcwt_soilt(:,1,:),            &
      smvccl_soilt(:,1,:), albsoil_soilt(:,1),                                 & 
      
      ! IN: Met forcing: 
      lw_down, sw_surft, ls_rain, ls_snow,                                     &
      tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy,                      & 
      ! This an outlier IN here. INOUT @ implicit. (was)OUT at extras
      ! I think we are dealing with it OK now but confusion could be removed  
      snow_surft,                                                              &
      
      ! IN: canopy height, LAI seasonally presecribed, potentially prognostic 
      ! IN: CO2 mass mixing ratio  
      canht_pft, lai_pft, co2_mmr,                                             & 
      
      ! TYPEs passed from top_level to maintain scope, access to UM STASH 
      ! IN: tiled soil/snow prognostics - IN here. INOUT @ implicit  
      ! INOUT: Carries fields needed by CABLE b/n pathways (rad, explicit etc) 
      !        Currently carrying CABLE TYPEs (canopy%, rad% etc).   
      ! IN: pars carries vegin/soilin - potentially redundant 
      progs_cbl, work_cbl, pars_io_cbl,                                        &
      
      ! OUT: UM fields UNPACKed from CABLE (@ explicit) 
      ftl_surft, fqw_surft, tstar_surft, dtstar_surft, u_s, u_s_std_surft,     &
      cd_surft, ch_surft, radnet_surft, fraca, resfs, resft, z0h_surft,        &
      z0m_surft, recip_l_MO_surft, epot_surft, npp_pft_acc, resp_w_pft_acc )
      
! hard-wire CABLE emissivities into respective JULES vars - CM3#55-19
! hcons_surf and hcons_snow should not need filling again
DO l=1,land_pts
   emis_soil(l) = 1.0
END DO

DO n=1,nsurft
   DO l=1,land_pts
      emis_surft(l,n) = 1.0
      flake(l,n) = 0.0  !CM2-era overwrite - CM#55-24
   END DO
END DO

!pass CABLE non-topography impacted values for drag & fr velo into JULES variables
DO n=1,nsurft
   DO l=1,land_pts
      cd_std(l,n) = cd_surft(l,n)
      v_s_surft(l,n) = u_s_std_surft(l,n)  !estimate of fr vel with topography
      v_s_std(l,n) = u_s_std_surft(l,n)    !value of fraction velocity no topography
      z0m_eff_surft(l,n) = z0m_surft(l,n)  !default topograhic roughness (?needed)
   END DO
END DO


!CM3#55-22 CM3#55-23
!recip_l_mo_surft = 1.0 !?
!v_s_std = 1.0e-6 !?
!CM3#55-24 - noting the link to resfs
!flake =0.0 !CABLE_LSM:CM2clobbered flake

!CM3#55 - from here on goes a set of mappings from CABLE to JULES and/or from
!         JULES to JULES given a CABLE variable has been given earlier.
!
! variables that need thought are resft, hcons_soilt, hcons_snow, cd_surft, ch_surft,
! other parts of classic aerosol and dust code, ashtf_surft and ashtf_prime_surft
!
! we are looking to avoid any calls to fcdch (the science is complicated and does not
! follow CABLE science) and sf_flux
!CABLE_LSM: End
!==============================================================================

!==============================================================================
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
!==============================================================================

!CM3#55-21 - hoping to be able to set ashtf_surft = 1.0 and ahstf_prime_surft=1.0
! and bypass all of this
!DO n = 1,nsurft
!
!  !Set the current soil tile (see notice above)
!  IF (nsoilt == 1) THEN
!    !There is only 1 soil tile
!    m = 1
!  ELSE ! nsoilt == nsurft
!    !Soil tiles map directly on to surface tiles
!    m = n
!  END IF !nsoilt
!
!  ! Set up soil and surface thermal properties
!!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE)                               &
!!!$OMP PRIVATE(i, j, l, ds_ratio)                                               &
!!$OMP SHARED(land_pts, land_index, t_i_length, tsurf, t_soil_soilt, t_elev,    &
!!$OMP        tl_1,                                                             &
!!$OMP        dzsurf, dzsoil, hcons_surf, hcons_soilt, canhc_surf, canhc_surft, &
!!$OMP        nsmax, nsnow_surft, tsnow_surft, ds_surft, hcons_snow,            &
!!$OMP        cansnowtile, snowdepth_surft, snow_hcon,                          &
!!$OMP        l_snow_nocan_hc, l_flake_model, l_aggregate, lake,                &
!!$OMP        n, m, ts1_lake_gb, hcon_lake, l_elev_land_ice, l_lice_point,      &
!!$OMP        tsurf_elev_surft, dzsoil_elev, l_lice_surft, hcondeep,            &
!!$OMP        l_moruses_storage, urban_roof, l_fix_moruses_roof_rad_coupling,   &
!!$OMP        vfrac_surft, ashtf_surft, scaling_urban, l_soil_point)            &
!!$OMP SCHEDULE(STATIC)
!  DO l = 1,land_pts
!    j = (land_index(l) - 1) / t_i_length + 1
!    i = land_index(l) - (j-1) * t_i_length
!    hcons_surf(l,n) = hcons_soilt(l,m)
!    IF (l_elev_land_ice .AND. l_lice_point(l)) THEN
!
!      ! Land ice
!      tsurf(l,n)        = tsurf_elev_surft(l,n)
!      dzsurf(l,n)       = dzsoil_elev
!      canhc_surf(l,n)   = 0.0
!      IF (l_lice_surft(n)) THEN
!        hcons_surf(l,n) = snow_hcon
!      ELSE
!        hcons_surf(l,n) = hcondeep
!      END IF
!    ELSE
!
!      ! Soil
!      tsurf(l,n)      = t_soil_soilt(l,m,1) + t_elev(l,n) - tl_1(i,j)
!      dzsurf(l,n)     = dzsoil(1)
!      hcons_surf(l,n) = hcons_soilt(l,m)
!      canhc_surf(l,n) = canhc_surft(l,n)
!    END IF
!    IF (     (l_flake_model   )                                                &
!        .AND. ( .NOT. l_aggregate)                                             &
!        .AND. (n == lake  )) THEN
!
!      ! Lake
!      tsurf(l,n)      = ts1_lake_gb(l) + t_elev(l,n) - tl_1(i,j)
!      hcons_surf(l,n) = hcon_lake(l)
!      dzsurf(l,n)     = dzsoil(1)
!    END IF
!    IF ((nsmax > 0) .AND. (nsnow_surft(l,n) > 0)) THEN
!
!      ! Snow
!      tsurf(l,n) = tsnow_surft(l,n,1)
!      ! change the effective surface layer thickness for snow
!      dzsurf(l,n)     = ds_surft(l,n,1)
!      hcons_surf(l,n) = hcons_snow(l,n)
!      IF ( ( .NOT. cansnowtile(n)) .AND. l_snow_nocan_hc .AND.                 &
!           (nsmax > 0) .AND. (nsnow_surft(l,n) > 0) ) canhc_surf(l,n) = 0.0
!    END IF
!
!    ! MORUSES: Uncouple the roof for perfect insulation. hcons should only be zero
!    ! to change the "conductive" coupling to "uncoupled" otherwise it is radiatively
!    ! coupled.
!    IF ( .NOT. l_aggregate .AND. l_moruses_storage .AND. n == urban_roof ) THEN
!      IF ( l_fix_moruses_roof_rad_coupling ) THEN
!        IF (vfrac_surft(l,n) == 0.0) THEN
!          hcons_surf(l,n) = 0.0
!        END IF
!      ELSE
!        hcons_surf(l,n)   = 0.0
!      END IF
!    END IF
!
!    ! Set up surface soil condictivity
!    ashtf_surft(l,n) = 2.0 *hcons_surf(l,n) / MAX( dzsurf(l,n), dzsoil(1) )
!    ! Except when n == urban_canyon when MORUSES is used
!    ! scaling_urban(l) = 1.0
!    ashtf_surft(l,n) = ashtf_surft(l,n) * scaling_urban(l,n)
!
!    ! Adjust surface soil condictivity for snow
!    IF (snowdepth_surft(l,n) > 0.0 .AND. l_soil_point(l)                       &
!        .AND. nsnow_surft(l,n) == 0 ) THEN
!      IF ( l_moruses_storage .AND. n == urban_roof ) THEN
!        ! This required as HCONS(L) = 0 in this case.
!        ashtf_surft(l,n)   =  0.0
!      ELSE
!        ds_ratio           = 2.0 * snowdepth_surft(l,n) / dzsurf(l,n)
!        IF (ds_ratio <= 1.0) THEN
!          ashtf_surft(l,n) =  ashtf_surft(l,n) /                               &
!                      (1.0 + ds_ratio * (hcons_surf(l,n) / snow_hcon - 1.0))
!        ELSE
!          ashtf_surft(l,n) =  ashtf_surft(l,n) *                               &
!                                    snow_hcon / hcons_surf(l,n)
!        END IF
!      END IF
!    END IF
!
!  END DO
!!$OMP END PARALLEL DO
!END DO
!==============================================================================
! *END NOTICE REGARDING SOIL TILING**
!==============================================================================


!-----------------------------------------------------------------------
!  3. Calculation of transfer coefficients and surface layer stability
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  3.1 Calculate neutral roughness lengths
!-----------------------------------------------------------------------

! Land tiles
!jh:could be able to dodge all these roughness calcs - test re-instating them
!jh!alla CM2
!jh! Z0_SURFT contains the appropriate value for land-ice points, but has to
!jh! be modified for snow-cover on non-land-ice points. Thermal roughness
!jh! lengths are set to be proportional to the tiled roughness length in
!jh! the case of multiple tiles, but if tiled properties have already
!jh! been aggregated, an adjustment for snow cover is required. In this case
!jh! a ratio of 0.1 between the thermal and momentum roughness lengths over
!jh! snow is assumed; to do otherwise would require reaggregation. In the
!jh! case of multiple tiles, the assignment is delayed until the urban
!jh! options have been considered.
!jh!
!jh!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n, z0, zeta1)     &
!jh!$OMP          SHARED(nsurft, surft_pts, surft_index, snow_surft,              &
!jh!$OMP                 l_soil_point, z0_surft, snowdep_surft, l_aggregate,      &
!jh!$OMP                 i_aggregate_opt, z0h_surft_bare, z0m_surft,              &
!jh!$OMP                 l_moruses_rough_surft, z0h_z0m, z0h_surft,               &
!jh!$OMP                 z0h_surft_classic,                                       &
!jh!$OMP                 z0h_z0m_classic)   SCHEDULE(STATIC)
!jh
!jh!CM3#55-22 - only need to fill a value for classic stuff here 
!jh!DO n = 1,nsurft
!jh  !! MORUSES parameterises z0m and z0h differently and are not affected by snow
!jh  !! in the same way so are dealt with independently
!jh  !IF ( .NOT. l_moruses_rough_surft(n) ) THEN
!jh  !  DO k = 1,surft_pts(n)
!jh  !    l = surft_index(k,n)
!jh      !IF ( snow_surft(l,n) > 0.0 .AND. l_soil_point(l) ) THEN
!jh      !  z0 = z0_surft(l,n) - 0.1 * snowdep_surft(l,n)
!jh      !  zeta1 = MIN( 5.0e-4 , z0_surft(l,n)  )
!jh      !  z0m_surft(l,n) = MAX( zeta1 , z0 )
!jh      !  !     Set z0h_surft explicitly if this option is selected,
!jh      !  !     otherwise, it will be set for the first tile below.
!jh      !  IF (l_aggregate .AND. i_aggregate_opt == 1) THEN
!jh      !    z0 = z0h_surft_bare(l,n) - 0.1 * 0.1 * snowdep_surft(l,n)
!jh      !    zeta1 = MIN( 5.0e-5 , z0h_surft_bare(l,n)  )
!jh      !    z0h_surft(l,n) = MAX( zeta1 , z0 )
!jh      !  END IF
!jh      !ELSE
!jh      !  z0m_surft(l,n) = z0_surft(l,n)
!jh      !  !     Set z0h_surft explicitly if this option is selected,
!jh      !  !     otherwise, it will be set for the first tile below.
!jh      !  IF (l_aggregate .AND. i_aggregate_opt == 1)                            &
!jh      !     z0h_surft(l,n) = z0h_surft_bare(l,n)
!jh      !END IF
!jh
!jh      !!   Set the thermal roughness length if aggregation is not being
!jh      !!   carried out, or if the original scheme is being used.
!jh      !!   It must be done here for consistency with the urban options.
!jh      !IF ( ( .NOT. l_aggregate) .OR.                                           &
!jh      !   (l_aggregate .AND. i_aggregate_opt == 0) )                            &
!jh      !   z0h_surft(l,n) = z0h_z0m(n) * z0m_surft(l,n)
!jh
!jh      ! Also set additional roughness length for use in CLASSIC aerosol
!jh      ! deposition scheme - CM3 retain - now done later
!jh      !z0h_surft_classic(l,n) = z0h_z0m_classic(n) * z0m_surft(l,n)
!jh   ! END DO
!jh  !END IF
!jh!END DO
!jh!$OMP END PARALLEL DO
!jh
!jh!CM3#55-22
!jh!! MORUSES does not yet contain a parametrisation for snow, which should not
!jh!! affect the behaviour of bluff bodies. It will instead affect the material
!jh!! roughness of the road & roof and not the walls, which are essentially
!jh!! snow-free. Snow could be added to the material roughness length for
!jh!! momentum before passing to urbanz0, which calculates roughness length for heat
!jh!! Two calls are required; one for urban_canyon and one for urban_roof.
!jh!!
!jh!! If l_aggregate MORUSES roughness lengths are set in the usual way above with
!jh!! z0_surft = ztm in sparm
!jh!IF ( .NOT. l_aggregate ) THEN
!jh!  IF ( ANY( l_moruses_rough_surft(1:nsurft) ) ) THEN
!jh!    n = urban_canyon
!jh!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,l,j,i)              &
!jh!!$OMP SHARED(surft_pts,n,surft_index,land_index,t_i_length,z1_uv,z1_tq,hgt_gb, &
!jh!!$OMP hwr_gb,disp_gb,z0m_surft,z0_surft,z0h_surft,urban_roof,                  &
!jh!!$OMP z0h_surft_classic)
!jh!    DO k = 1,surft_pts(n)
!jh!      l = surft_index(k,n)
!jh!      z0m_surft(l,n)          = z0_surft(l,n)
!jh!      z0m_surft(l,urban_roof) = z0_surft(l,urban_roof)
!jh!      j = ( land_index(l) - 1 ) / t_i_length + 1
!jh!      i = land_index(l) - ( j - 1 ) * t_i_length
!jh!      CALL urbanz0(                                                            &
!jh!         n, z1_uv(i,j), z1_tq(i,j), hgt_gb(l), hwr_gb(l), disp_gb(l),          &
!jh!         z0m_mat, z0m_surft(l,n), z0h_surft(l,n) )
!jh!      CALL urbanz0(                                                            &
!jh!         urban_roof, z1_uv(i,j), z1_tq(i,j), hgt_gb(l), hwr_gb(l), disp_gb(l), &
!jh!         z0m_mat, z0m_surft(l,urban_roof), z0h_surft(l,urban_roof) )
!jh!      ! Make CLASSIC aerosol roughness length for urban tiles consistent
!jh!      ! with those for heat and momentum
!jh!      z0h_surft_classic(l,n)          = z0h_surft(l,n)
!jh!      z0h_surft_classic(l,urban_roof) = z0h_surft(l,urban_roof)
!jh!    END DO
!jh!!$OMP END PARALLEL DO
!jh!  END IF
!jh!END IF
!jh
!jh!CM3#55-22
!jh!! Calculate roughness length affected by roughness sublayer in neutral
!jh!! conditions.
!jh!DO n = 1,nsurft
!jh!  IF (l_vegdrag_surft(n)) THEN
!jh!    CALL can_drag_z0(                                                          &
!jh!      land_pts, surft_pts(n), surft_index(:,n),                                &
!jh!      array_zero, canht_pft(:,n), lai_pft(:,n),                                &
!jh!      z0m_surft(:,n), z0h_surft(:,n), zdt_surft(:,n))
!jh!  END IF
!jh!END DO

! Calculate orographic effective parameter for neutral conditions
! if using orographic roughness scheme
IF (formdrag == effective_z0) THEN
  DO n = 1,nsurft
    CALL sf_orog (                                                             &
     land_pts,surft_pts(n),land_index,surft_index(:,n),                        &
     ho2r2_orog,rib_surft(:,n),sil_orog_land,z0m_surft(:,n),z1_uv,             &
     wind_profile_factor(:,n),z0m_eff_surft(:,n)                               &
     )
  END DO
ELSE
  wind_profile_factor(:,:) = 1.0
  z0m_eff_surft(:,:)       = z0m_surft(:,:)
END IF


!CM3#55 - CABLE now called earlier
!!CABLE_LSM:CM2.5 Oddly tstar_surft is declared as INTENT(IN) above!!
!CALL cable_explicit_main(                                                      &
!          mype, timestep, timestep_number,                                     &
!          tdims%i_end,tdims%j_end, land_pts, nsurft, npft, sm_levels, dzsoil,  &
!          land_index, surft_pts, surft_index, canht_pft, lai_pft, Flandg,      &
!          co2_mmr, tile_frac, cos_zenith_angle, latitude, longitude,           &
!          bexp_soilt(:,1,:), hcon_soilt(:,1), satcon_soilt(:,:),               &
!          sathh_soilt(:,1,:), smvcst_soilt(:,1,:), smvcwt_soilt(:,1,:),        &
!          smvccl_soilt(:,1,:), sthu_soilt(:,1,:), albsoil_soilt(:,1),          & 
!          lw_down, sw_surft,                                                   &
!          ls_rain, ls_snow,                                                    &
!          !not fully checked. iup to here all IN ! unpacked ? 
!          tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy, snow_surft,      &
!          progs_cbl, work_cbl, pars_io_cbl,                                    &
!          ftl_surft, fqw_surft,                                                &
!          tstar_surft, u_s, u_s_std_surft, cd_surft, ch_surft,                 &
!          radnet_surft, fraca, resfs, resft, z0h_surft, z0m_surft,             &
!          recip_l_MO_surft, epot_surft, &
!          npp_pft_acc, resp_w_pft_acc )
!
!CD_STD = CD_surft
!V_S_surft = U_S_STD_surft
!V_S_STD = U_S_STD_surft
!z0m_eff_surft(:,:)       = z0m_surft(:,:)
!
!recip_l_mo_surft = 1.0 !?
!v_s_std = 1.0e-6 !?
!flake =0.0 !CABLE_LSM:CM2clobbered flake
!!CABLE_LSM: End

!-----------------------------------------------------------------------
! Calculate RESFT with neutral CH and EPDT = 0 for use in calculation
! of Richardson number. RESFT=1 for snow and land-ice.
!-----------------------------------------------------------------------

!CM3#55-25 - noting that ideally we'd move the evaluation of resft to earlier 
!DO n = 1,nsurft
  !IF (l_vegdrag_surft(n)) THEN
  !  CALL can_drag_phi_m_h(                                                     &
  !    land_pts, surft_pts(n), surft_index(:,n), land_index,                    &
  !    array_zero, z1_uv, z1_tq, canht_pft(:,n), lai_pft(:,n),                  &
  !    z0m_surft(:,n), z0h_surft(:,n),                                          &
  !    phi_m, phi_h)
  !
!!$OMP PARALLEL DO IF(surft_pts(n)>1) DEFAULT(NONE)                             &
!!$OMP PRIVATE(i, j, k, l)                                                      &
!!$OMP SHARED(surft_pts, surft_index, t_i_length, land_index,                   &
!!$OMP        phi_m, phi_h, chn, wind_profile_factor, dq, qw_1,                 &
!!$OMP        qstar_surft, epdt, n, resft, resfs, flake, fraca) SCHEDULE(STATIC)
!    DO k = 1,surft_pts(n)
!      l = surft_index(k,n)
!      j=(land_index(l) - 1) / t_i_length + 1
!      i = land_index(l) - (j-1) * t_i_length
!      chn(l,n) = vkman**2 / (phi_m(l) * phi_h(l)) * wind_profile_factor(l,n)
!      dq(l) = qw_1(i,j) - qstar_surft(l,n)
!      epdt(l) = 0.0
!      !CABLE_LSM:CM2 
!      resft(l,n) =  MIN(1.,flake(l,n) + (1. - flake(l,n)) *        &
!                      ( fraca(l,n) + (1. - fraca(l,n))*resfs(l,n) )) 
!    END DO
!!$OMP END PARALLEL DO

!  ELSE ! l_vegdrag_surft = F

!!$OMP PARALLEL DO IF(surft_pts(n)>1) DEFAULT(NONE)                             &
!!$OMP PRIVATE(i, j, k, l, zetah, zetam)                                        &
!!$OMP SHARED(surft_pts, surft_index, t_i_length, land_index, z1_uv, resft,     &
!!$OMP        z0m_surft, z1_tq, z0h_surft, chn, wind_profile_factor, dq, qw_1,  &
!!$OMP        qstar_surft, epdt,resfs,fraca, flake,n)    SCHEDULE(STATIC)
    !DO k = 1,surft_pts(n)
    !  l = surft_index(k,n)
    !  j=(land_index(l) - 1) / t_i_length + 1
    !  i = land_index(l) - (j-1) * t_i_length
    !  !zetam = LOG ( (z1_uv(i,j) + z0m_surft(l,n)) / z0m_surft(l,n) )
    !  !zetah = LOG ( (z1_tq(i,j) + z0m_surft(l,n)) / z0h_surft(l,n) )
    !  !chn(l,n) = (vkman / zetah) * (vkman / zetam) *                           &
    !  !  wind_profile_factor(l,n)
    !  !dq(l) = qw_1(i,j) - qstar_surft(l,n)
    !  !epdt(l) = 0.0
    !  !CABLE_LSM: from 8.5
    !  !CM3#55-25 - now done later
    !  !resft(l,n) =  MIN(1.,flake(l,n) + (1. - flake(l,n)) *        &
    !  !              ( fraca(l,n) + (1. - fraca(l,n))*resfs(l,n) )) 
   ! END DO
!!$OMP END PARALLEL DO

!  END IF

!  ! We should only attempt to access sf_diag%resfs_stom(:,n) if it has
!  ! been fully allocated.
!  IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
!    n_diag = n
!  ELSE
!    n_diag = 1
!  END IF

!CABLE_LSM:CM2
!CM2!  CALL sf_resist (                                                             &
!CM2!   land_pts,surft_pts(n),land_index,surft_index(:,n),cansnowtile(n),           &
!CM2!   canopy(:,n),catch(:,n),chn(:,n),dq,epdt,flake(:,n),gc_surft(:,n),           &
!CM2!   gc_stom_surft(:,n),snowdep_surft(:,n),snow_surft(:,n),vshr_land,            &
!CM2!   fraca(:,n),resfs(:,n),resft(:,n),                                           &
!CM2!   sf_diag%resfs_stom(:,n_diag),sf_diag%l_et_stom,sf_diag%l_et_stom_surft)

!END DO

!CM3#55-26
!-----------------------------------------------------------------------
!  3.2 Calculate bulk Richardson number for the lowest model level.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  CALL sf_rib (                                                                &
   land_pts,surft_pts(n),land_index,surft_index(:,n),                          &
   bq_1,bt_1,qstar_surft(:,n),q_elev(:,n),resft(:,n),t_elev(:,n),              &
   tstar_surft(:,n),vshr_land,z0h_surft(:,n),z0m_surft(:,n),zdt_surft(:,n),    &
   z1_tq,z1_uv,l_vegdrag_surft(n),                                             &
   rib_surft(:,n),db_surft(:,n)                                                &
   )
END DO

!-----------------------------------------------------------------------
!  3.3 Calculate stability corrected effective roughness length.
!  Stability correction only applies to land points.
!-----------------------------------------------------------------------

IF (formdrag == effective_z0) THEN
  DO n = 1,nsurft
    CALL sf_orog (                                                             &
     land_pts,surft_pts(n),land_index,surft_index(:,n),                        &
     ho2r2_orog,rib_surft(:,n),sil_orog_land,z0m_surft(:,n),z1_uv,             &
     wind_profile_factor(:,n),z0m_eff_surft(:,n)                               &
     )
  END DO
END IF

!CM3#53 - once the roughnesses have been fully evaluated including topograhpic effects
!then cd_surft, driction velocities etc. needs to be revised here
! also do ashtf_surft and ashtf_prime_surft here for readability (could go ealier)
!
!what appears to be needed is
! 1. cd_surft needs to be adjusted to reflect topographic effects.
! 2. cd_std remains unimpacted by topography (may not be needed)
! 3. u_s_std_surft needs to be unimpacted by topgraphic effects
! 4. _classic variables need to use the classic values where
! z0h is replaced by z0h_classic and z0m is replaced by z0m_eff
!
! the equations in fcdch() are that
! v_s = v_s_std / windfactor
! cd_surft = cd_std * v_s/v_s_std * windfactor
! where v_s_std and cd_std are those provided by CABLE (tiled, no topography)
!
! and u_s_std_surft = v_s_std
!
! since we don't know whether z0h from CABLE = z0h_classic we will
! need to recalculate the classic outputs but chose to approximate this as
! cd_std_classic = cd_std = cd_surft (from CABLE) but
! ch_surft_classic = ch_surft * log(z1_tq/z0h_surft) / log(z1_tq/z0h_classic)

! first, evaluations that are needed regardless of configuration
! - these could go earlier
DO n=1,nsurft
   DO l=1,land_pts
      !fill ashtf_surft and ashtf_prime with default value
      !this needs to be non-zero to function with im_sf_pt2
      ashtf_surft(l,n) = 1.0
      ashtf_prime_surft(l,n) = 1.0

      !fill in resft
      resft(l,n) =  MIN(1.,flake(l,n) + (1. - flake(l,n)) *        &
                    ( fraca(l,n) + (1. - fraca(l,n))*resfs(l,n) )) 
   END DO
END DO

!if topographic form drag enabled revise the tile drag coefficient
IF (formdrag == effective_z0) THEN
   DO n=1,nsurft
      DO l=1,land_pts
         v_s_surft(l,n) = v_s_std(l,n) / wind_profile_factor(l,n)
         cd_surft(l,n) = cd_std(l,n) * v_s_surft(l,n) / v_s_std(l,n) * &
              wind_profile_factor(l,n)
      END DO
   END DO
END IF

!if classic aerosol scheme used evaluate the tile aerosol deposition coefficients 
IF (l_aero_classic) THEN
   DO n=1,nsurft
      DO k = 1,surft_pts(n)
         l = surft_index(k,n)
         j=(land_index(l) - 1) / t_i_length + 1
         i = land_index(l) - (j-1) * t_i_length
         z0h_surft_classic(l,n) = z0h_z0m_classic(n) * z0m_surft(l,n)

         !using zetam and zetah as temporary vars
         zetam = LOG ( (z1_tq(i,j) + z0m_surft(l,n)) / z0h_surft(l,n) )
         zetah = LOG ( (z1_tq(i,j) + z0m_surft(l,n)) / z0h_surft_classic(l,n) )
      
         cd_std_classic(l,n) = cd_std(l,n)
         ch_surft_classic(l,n) = ch_surft(l,n) * zetam / zetah
      END DO
   END DO
END IF

         
!-----------------------------------------------------------------------
!  3.4 Calculate CD, CH via routine FCDCH.
!-----------------------------------------------------------------------

!CABLE_LSM:CM2
!CM2!! Land tiles
!CM2!DO n = 1,nsurft
!CM2!  IF (l_vegdrag_surft(n)) THEN
!CM2!    n_veg = n
!CM2!    z0m_eff_surft(:,n) = z0m_surft(:,n)
!CM2!  ELSE
!CM2!    n_veg = 1
!CM2!  END IF
!CM2!  CALL fcdch (                                                                 &
!CM2!    cor_mo_iter,land_pts,surft_pts(n),                                         &
!CM2!    surft_index(:,n),land_index,                                               &
!CM2!    db_surft(:,n),vshr_land,                                                   &
!CM2!    z0m_eff_surft(:,n),z0h_surft(:,n),zdt_surft(:,n),zh,                       &
!CM2!    z1_uv,z1_uv_top,z1_tq,z1_tq_top,wind_profile_factor(:,n),                  &
!CM2!    ddmfx,ip_ss_solid,charnock,                                                &
!CM2!    charnock_w,                                                                &
!CM2!    l_vegdrag_surft(n),canht_pft(:,n_veg),lai_pft(:,n_veg),                    &
!CM2!    nsnow_surft(:,n),n,l_mo_buoyancy_calc,cansnowtile(n),l_soil_point,         &
!CM2!    canopy(:,n),catch(:,n),flake(:,n),gc_surft(:,n),                           &
!CM2!    snowdep_surft(:,n),snow_surft(:,n),canhc_surf(:,n),                        &
!CM2!    dzsurf(:,n),qstar_surft(:,n),q_elev(:,n),radnet_surft(:,n),                &
!CM2!    snowdepth_surft(:,n),timestep,t_elev(:,n),tsurf(:,n),tstar_surft(:,n),     &
!CM2!    vfrac_surft(:,n),emis_surft(:,n),emis_soil,anthrop_heat_surft(:,n),        &
!CM2!    scaling_urban(:,n),alpha1(:,n),hcons_surf(:,n),ashtf_surft(:,n),           &
!CM2!    rhostar,bq_1,bt_1,                                                         &
!CM2!    cd_surft(:,n),ch_surft(:,n),cd_std(:,n),                                   &
!CM2!    v_s_surft(:,n),v_s_std(:,n),recip_l_mo_surft(:,n),                         &
!CM2!    u_s_iter_surft(:,n)                                                        &
!CM2!    )
!CM2!END DO
!CM2!
!CM2!! As roughness length have been changed by vegetation drag effect, effective
!CM2!! roughness length should be updated.
!CM2!DO n = 1,nsurft
!CM2!  IF (l_vegdrag_surft(n)) THEN
!CM2!    z0m_surft(:,n) = z0m_eff_surft(:,n)
!CM2!    IF (formdrag == effective_z0) THEN
!CM2!      CALL sf_orog (                                                           &
!CM2!       land_pts,surft_pts(n),land_index,surft_index(:,n),                      &
!CM2!       ho2r2_orog,rib_surft(:,n),sil_orog_land,z0m_surft(:,n),z1_uv,           &
!CM2!       wind_profile_factor(:,n),z0m_eff_surft(:,n)                             &
!CM2!       )
!CM2!    END IF
!CM2!  END IF
!CM2!END DO
!CM2!
!CM2!!$OMP PARALLEL IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n)                   &
!CM2!!$OMP          SHARED(cor_mo_iter, nsurft, surft_pts, surft_index,             &
!CM2!!$OMP                 u_s_iter_surft, u_s_std_surft, v_s_std)
!CM2!IF ( cor_mo_iter >= use_correct_ustar ) THEN
!CM2!  !       Use correct "standard" ustar
!CM2!!$OMP DO SCHEDULE(STATIC)
!CM2!  DO n = 1,nsurft
!CM2!    DO k = 1,surft_pts(n)
!CM2!      l = surft_index(k,n)
!CM2!      u_s_std_surft(l,n) = v_s_std(l,n)
!CM2!    END DO
!CM2!  END DO
!CM2!!$OMP END DO
!CM2!ELSE
!CM2!  !       Use ustar from mid-iteration
!CM2!!$OMP DO SCHEDULE(STATIC)
!CM2!  DO n = 1,nsurft
!CM2!    DO k = 1,surft_pts(n)
!CM2!      l = surft_index(k,n)
!CM2!      u_s_std_surft(l,n) = u_s_iter_surft(l,n)
!CM2!    END DO
!CM2!  END DO
!CM2!!$OMP END DO
!CM2!END IF
!CM2!!$OMP END PARALLEL

!-----------------------------------------------------------------------
!  3.5 Recalculate friction velocity for the dust scheme using the
!      bare soil roughness length if using only 1 aggregated tile
!-----------------------------------------------------------------------

!CM3#55-27
!IF ((l_dust .OR. l_dust_diag) .AND. l_aggregate) THEN
!
!  n = 1
!  ! Calculate z0m and z0h for bare soil
!  IF (l_vary_z0m_soil) THEN
!!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(k, l)             &
!!$OMP             SHARED(n, soil, surft_pts, surft_index, z0m_soil_in,         &
!!$OMP                   z0h_soil,z0h_z0m, z0m_soil, wind_profile_fac_soil)     &
!!$OMP             SCHEDULE(STATIC)
!    DO k = 1, surft_pts(n)
!      l = surft_index(k,n)
!      z0m_soil(l,n) = z0m_soil_in(l)
!      z0h_soil(l,n) = z0h_z0m(soil) * z0m_soil(l,n)
!      !     Set wind profile factor to 1 as not using orog term in z0m
!      wind_profile_fac_soil(l,n) = 1.0
!    END DO
!!$OMP END PARALLEL DO
!  ELSE
!!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(k, l)             &
!!$OMP             SHARED(n, soil, surft_pts, surft_index, z0_soil, z0h_soil,   &
!!$OMP                   z0h_z0m, z0m_soil, wind_profile_fac_soil)              &
!!$OMP             SCHEDULE(STATIC)
!    DO k = 1, surft_pts(n)
!      l = surft_index(k,n)
!      z0m_soil(l,n) = z0_soil
!      z0h_soil(l,n) = z0h_z0m(soil) * z0m_soil(l,n)
!      !     Set wind profile factor to 1 as not using orog term in z0m
!      wind_profile_fac_soil(l,n) = 1.0
!    END DO
!!$OMP END PARALLEL DO
!  END IF

!  ! Call fcdch again to calculate dust friction velocity on bare soil.
!  ! The canopy drag scheme is not available on the aggregated tile
!  ! and is disabled.
!  l_vegdrag_active_here = .FALSE.
!  CALL fcdch (                                                                 &
!    cor_mo_iter,land_pts,surft_pts(n),                                         &
!    surft_index(:,n),land_index,                                               &
!    db_surft(:,n),vshr_land,                                                   &
!    z0m_soil(:,n),z0h_soil(:,n),zdt_dummy,zh,                                  &
!    z1_uv,z1_uv_top,z1_tq,z1_tq_top,wind_profile_fac_soil(:,n),                &
!    ddmfx,ip_ss_solid,charnock,                                                &
!    charnock_w,                                                                &
!    l_vegdrag_active_here,array_zero,array_zero,                               &
!    nsnow_surft(:,n),n,.FALSE.,cansnowtile(n),l_soil_point,                    &
!    canopy(:,n),catch(:,n),flake(:,n),gc_surft(:,n),                           &
!    snowdep_surft(:,n),snow_surft(:,n),canhc_surf(:,n),                        &
!    dzsurf(:,n),qstar_surft(:,n),q_elev(:,n),radnet_surft(:,n),                &
!    snowdepth_surft(:,n),timestep,t_elev(:,n),tsurf(:,n),tstar_surft(:,n),     &
!    vfrac_surft(:,n),emis_surft(:,n),emis_soil,anthrop_heat_surft(:,n),        &
!    scaling_urban(:,n),alpha1(:,n),hcons_surf(:,n),ashtf_surft(:,n),           &
!    rhostar,bq_1,bt_1,                                                         &
!  ! Following tiled outputs (except v_s_std_soil and u_s_iter_soil)
!  ! are dummy variables not needed from this call
!    cd_surft_soil(:,n),ch_surft_soil(:,n),cd_std_soil(:,n),                    &
!    v_s_surft_soil(:,n),v_s_std_soil(:,n),recip_l_mo_surft_soil(:,n),          &
!    u_s_iter_soil(:,n)                                                         &
!    )

!!$OMP PARALLEL IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n)                   &
!!$OMP          SHARED(cor_mo_iter, nsurft, surft_index, surft_pts,             &
!!$OMP                 u_s_iter_soil, u_s_std_surft, v_s_std_soil)
!  IF ( cor_mo_iter >= use_correct_ustar ) THEN
!    !       Use correct "standard" ustar
!!$OMP DO SCHEDULE(STATIC)
!    DO n = 1,nsurft
!      DO k = 1,surft_pts(n)
!        l = surft_index(k,n)
!        u_s_std_surft(l,n) = v_s_std_soil(l,n)
!      END DO
!    END DO
!!$OMP END DO
!  ELSE
!    !       Use ustar from mid-iteration
!!$OMP DO SCHEDULE(STATIC)
!    DO n = 1,nsurft
!      DO k = 1,surft_pts(n)
!        l = surft_index(k,n)
!        u_s_std_surft(l,n) = u_s_iter_soil(l,n)
!      END DO
!    END DO
!!$OMP END DO
!  END IF
!!$OMP END PARALLEL
!
!END IF

!-----------------------------------------------------------------------
!  3.6 Recalculate cd, ch etc. using z0h=z0h_classic. The parameters
!       calculated using this additional roughness length are for
!       CLASSIC aerosol deposition only.
!-----------------------------------------------------------------------
!CM3#55-28 - we've implemented CABLE estimates for this earlier
! remove to facilitate simplification of ashtf evaluations
!IF (l_aero_classic) THEN
!  ! The canopy drag scheme is not supported for this aerosol scheme and
!  ! is turned off.
!  l_vegdrag_active_here = .FALSE.
!  ! Land tiles
!  DO n = 1,nsurft
!    CALL fcdch (                                                               &
!    ! Input variables identical to main call except using different z0h
!      cor_mo_iter,land_pts,surft_pts(n),                                       &
!      surft_index(:,n),land_index,                                             &
!      db_surft(:,n),vshr_land,                                                 &
!      z0m_eff_surft(:,n),z0h_surft_classic(:,n),zdt_dummy,zh,                  &
!      z1_uv,z1_uv_top,z1_tq,z1_tq_top,                                         &
!      wind_profile_factor(:,n),                                                &
!      ddmfx,ip_ss_solid,charnock,                                              &
!      charnock_w,                                                              &
!      l_vegdrag_active_here,array_zero,array_zero,                             &
!      nsnow_surft(:,n),n,.FALSE.,cansnowtile(n),l_soil_point,                  &
!      canopy(:,n),catch(:,n),flake(:,n),gc_surft(:,n),                         &
!      snowdep_surft(:,n),snow_surft(:,n),canhc_surf(:,n),                      &
!      dzsurf(:,n),qstar_surft(:,n),q_elev(:,n),radnet_surft(:,n),              &
!      snowdepth_surft(:,n),timestep,t_elev(:,n),tsurf(:,n),tstar_surft(:,n),   &
!      vfrac_surft(:,n),emis_surft(:,n),emis_soil,anthrop_heat_surft(:,n),      &
!      scaling_urban(:,n),alpha1(:,n),hcons_surf(:,n),ashtf_surft(:,n),         &
!      rhostar,bq_1,bt_1,                                                       &
!    ! Following tiled outputs (except cd_std_classic and ch_surft_classic)
!    ! are dummy variables not needed from this call
!      cd_surft_classic(:,n),ch_surft_classic(:,n),                             &
!      cd_std_classic(:,n),v_s_surft_classic(:,n),                              &
!      v_s_std_classic(:,n),recip_l_mo_surft_classic(:,n),                      &
!      u_s_iter_classic(:,n)                                                    &
!      )
!  END DO
!END IF

!-----------------------------------------------------------------------
! Calculate gridbox-means of transfer coefficients.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    cd_land(i,j) = cd_land(i,j) + tile_frac(l,n) * cd_surft(l,n)
    ch_land(i,j) = ch_land(i,j) + tile_frac(l,n) * ch_surft(l,n)
  END DO
END DO

! aerodynamic resistance diagnostic
IF (sf_diag%l_ra) THEN
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    sf_diag%ra(l) = 1.0 / (ch_land(i,j) * vshr_land(i,j))
  END DO
END IF

!-----------------------------------------------------------------------
!  4.3 Calculate the surface exchange coefficients RHOK(*) and
!       resistances for use in CLASSIC aerosol scheme
!       (Note that CD_STD, CH and VSHR should never = 0)
!     RHOSTAR * CD * VSHR stored for diagnostic output before
!     horizontal interpolation.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    rhokm_1_surft(l,n) = rhostar_mom(i,j) * cd_surft(l,n) * vshr_land(i,j)
    !                                                         ! P243.124
    rhokm_land(i,j) = rhokm_land(i,j) +                                        &
           tile_frac(l,n) * rhokm_1_surft(l,n)
    rhokh_surft(l,n) = rhostar_mom(i,j) * ch_surft(l,n) * vshr_land(i,j)
    !                                                         ! P243.125
  END DO
END DO

!-----------------------------------------------------------------------
!  Calculate local and gridbox-average surface fluxes of heat and
!  moisture.
!-----------------------------------------------------------------------

! Adjust ASHTF for sens. heat flux to ground beneath coniferous canopy
!CM3#55-missed in first pass - remove as can_model==4 not relevant to CABLE
!IF ( .NOT. l_aggregate .AND. can_model == 4) THEN
!  DO n = 1,npft
!    IF ( cansnowtile(n) ) THEN
!!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE( k, l, j, i)      &
!!$OMP SHARED(n, cd_surft, land_index, rhokh_can, rhostar, surft_pts,           &
!!$OMP        surft_index, t_i_length, vshr_land, cp) SCHEDULE(STATIC)
!      DO k = 1,surft_pts(n)
!        l = surft_index(k,n)
!        j=(land_index(l) - 1) / t_i_length + 1
!        i = land_index(l) - (j-1) * t_i_length
!        rhokh_can(l,n) = rhostar(i,j) * cp /                                   &
!                     (43.0 / (SQRT(cd_surft(l,n)) * vshr_land(i,j)))
!      END DO
!!$OMP END PARALLEL DO
!    END IF
!  END DO
!END IF

!-----------------------------------------------------------------------
!  4.1 Recalculate RESFT using "true" CH and EPDT for land tiles
!-----------------------------------------------------------------------

lh0 = lc
sea_point = 0.0

DO n = 1,nsurft
!!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)       &
!!$OMP SHARED(ch_surft, land_index, qstar_surft, qw_1, dq, epdt,n,              &
!!$OMP rhostar, surft_pts, surft_index, timestep, t_i_length,                   &
!!$OMP vshr_land) SCHEDULE(STATIC)

   !CM3#55-29 - we do not need ashtf_surft or ashtf_prime_surft to take a true value
   ! only non-zero is needed - see above - remove sf_flux_cbl
  !! Calcualte humidity gradient and rate of change of potential evaporation
  !! with time
  !DO k = 1,surft_pts(n)
  !  l = surft_index(k,n)
  !  j = (land_index(l) - 1) / t_i_length + 1
  !  i = land_index(l) - (j-1) * t_i_length
  !  dq(l)   = qw_1(i,j) - qstar_surft(l,n)
  !  epdt(l) = - rhostar(i,j) * ch_surft(l,n) * vshr_land(i,j)                  &
  !    *dq(l) * timestep
  !END DO
!!$OMP END PARALLEL DO

!  ! We should only attempt to access sf_diag%resfs_stom(:,n) if it has
!  ! been fully allocated.
!  IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
!    n_diag = n
!  ELSE
!    n_diag = 1
!  END IF

!CABLE_LSM:CM2{ This is all being done above BUT to ashtf_surft - in CM2 we used
!elements of sf_flux to update ashtf_prime and others, however it seems the
!others were unnecessary
!CM2!  CALL sf_resist (                                                             &
!CM2!   land_pts,surft_pts(n),land_index,surft_index(:,n),cansnowtile(n),           &
!CM2!   canopy(:,n),catch(:,n),ch_surft(:,n),dq,epdt,flake(:,n),gc_surft(:,n),      &
!CM2!   gc_stom_surft(:,n),snowdep_surft(:,n),snow_surft(:,n),vshr_land,            &
!CM2!   fraca(:,n),resfs(:,n),resft(:,n),                                           &
!CM2!   sf_diag%resfs_stom(:,n_diag),sf_diag%l_et_stom,sf_diag%l_et_stom_surft)

!  CALL sf_flux_cbl (                                                               &
!   land_pts,surft_pts(n),                                                      &
!   land_index,surft_index(:,n),                                                &
!   nsnow_surft(:,n),n,canhc_surf(:,n),dzsurf(:,n),hcons_surf(:,n),             &
!   ashtf_surft(:,n),qstar_surft(:,n),q_elev(:,n),                              &
!   radnet_surft(:,n),resft(:,n),rhokh_surft(:,n),l_soil_point,                 &
!   snowdepth_surft(:,n),timestep,t_elev(:,n),tsurf(:,n),                       &
!   tstar_surft(:,n),vfrac_surft(:,n),rhokh_can(:,n),z0h_surft(:,n),            &
!   z0m_eff_surft(:,n),zdt_surft(:,n),z1_tq,lh0,emis_surft(:,n),emis_soil,      &
!   1.0,anthrop_heat_surft(:,n),scaling_urban(:,n),l_vegdrag_surft(n),          &
!   alpha1(:,n),ashtf_prime_surft(:,n),fqw_surft(:,n),                          &
!   epot_surft(:,n),ftl_surft(:,n),dtstar_surft(:,n),sea_point                  &
! )

  ! update gridbox means and diagnostics
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)       &
!$OMP SHARED(n, surft_pts, surft_index, land_index, t_i_length,                &
!$OMP        l_irrig_dmd, l_aggregate, can_model, cansnowtile,                 &
!$OMP        snow_surft, resfs_irr_surft, resfs, gc_irr_surft,                 &
!$OMP        ftl_1, flandg, tile_frac, ftl_surft, fqw_1, fqw_surft,            &
!$OMP        ch_surft, vshr_land) SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j = (land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length

    ! Calculate gridbox mean fluxes of heat and moisture
    ftl_1(i,j) = ftl_1(i,j) + flandg(i,j) * tile_frac(l,n) * ftl_surft(l,n)
    fqw_1(i,j) = fqw_1(i,j) + flandg(i,j) * tile_frac(l,n) * fqw_surft(l,n)

    !CM3#55-29 - may need the OMP directives changed
    !! Calculate surface resistance term for irrigated surfaces
    !IF (l_irrig_dmd) THEN
    !  IF ( .NOT. l_aggregate .AND. can_model == 4 .AND. cansnowtile(n) .AND.   &
    !         snow_surft(l,n) > 0.0) THEN
    !    resfs_irr_surft(l,n) = resfs(l,n)
    !  ELSE
    !    resfs_irr_surft(l,n) = gc_irr_surft(l,n) /                             &
    !                ( gc_irr_surft(l,n) + ch_surft(l,n) * vshr_land(i,j) )
    !  END IF
    !END IF
  END DO
!$OMP END PARALLEL DO

END DO

! Calculate surface momentum flux
IF (sf_diag%l_tau_surft) THEN
  DO n = 1,nsurft
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)       &
!$OMP SHARED(n, surft_pts, surft_index, land_index, t_i_length,                &
!$OMP        rhokm_1_surft, vshr_land, sf_diag) SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j = (land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      sf_diag%tau_surft(l,n) = rhokm_1_surft(l,n) * vshr_land(i,j)
    END DO
!$OMP END PARALLEL DO
  END DO
END IF

! Set surface stress on tiles diagnostic
IF (sf_diag%l_tau_1) THEN
  DO n = 1,nsurft
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)       &
!$OMP SHARED(n, surft_pts, surft_index, land_index, t_i_length,                &
!$OMP        flandg, tile_frac, rhokm_1_surft, vshr_land,                      &
!$OMP        sf_diag) SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      sf_diag%tau_1(i,j) = sf_diag%tau_1(i,j) +                                &
                           flandg(i,j) * tile_frac(l,n) *                      &
                           rhokm_1_surft(l,n) * vshr_land(i,j)
    END DO
!$OMP END PARALLEL DO
  END DO
END IF


!-----------------------------------------------------------------------
!  4.4   Calculate the standard deviations of layer 1 turbulent
!        fluctuations of temperature and humidity using approximate
!        formulae from first order closure.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  CALL stdev1 (                                                                &
   land_pts,surft_pts(n),land_index,surft_index(:,n),flandg,                   &
   bq_1,bt_1,fqw_surft(:,n),ftl_surft(:,n),rhokm_1_surft(:,n),                 &
   rhostar,vshr_land,z0m_surft(:,n),z1_tq,tile_frac(:,n),                      &
   q1_sd,t1_sd                                                                 &
   )
END DO

!-----------------------------------------------------------------------
! Call SFL_INT to calculate CDR10M and CHR1P5M - interpolation coeffs
! used to calculate screen temperature, humidity and 10m winds.
!-----------------------------------------------------------------------

! Set flag if cdr10m is to be calculated for use with snow unloading.
! If l_fix_wind_snow=F this calculation will depend on the su10 and sv10
! switches in sf_diag.
l_cdr10m_snow = .FALSE.
!CM3#55-31
!IF ( l_fix_wind_snow ) THEN
!  DO n = 1,nsurft
!    IF ( canSnowTile(n) .AND. unload_rate_u(n) /= 0.0 ) l_cdr10m_snow = .TRUE.
!  END DO
!END IF

IF (sf_diag%suv10m_n) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j)                  &
!$OMP SHARED(pdims,sf_diag)
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      sf_diag%cdr10m_n(i,j) = 0.0
      sf_diag%cd10m_n(i,j)  = 0.0
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Land tiles
IF (sf_diag%su10 .OR. sf_diag%sv10 .OR. sf_diag%sq1p5 .OR.                     &
    sf_diag%st1p5 .OR. sf_diag%suv10m_n .OR.                                   &
    l_cdr10m_snow .OR.                                                         &
    (IScrnTDiag == IP_ScrnDecpl2) .OR.                                         &
    (IScrnTDiag == IP_ScrnDecpl3) ) THEN
  DO n = 1,nsurft

     !CM2!!CABLE_LSM:land_pts=0 crashes
     !CM3#55 - need to look at this in more detail - not convinced this is correct
    CALL sfl_int (                                                             &
     land_pts,surft_pts(n),l_cdr10m_snow,surft_index(:,n),land_index,flandg,   &
     vshr_land,cd_std(:,n),cd_surft(:,n),ch_surft(:,n),                        &
     tile_frac(:,n),                                                           &
     z0m_eff_surft(:,n),z0m_surft(:,n),z0h_surft(:,n),                         &
     recip_l_mo_surft(:,n),                                                    &
     v_s_surft(:,n),v_s_std(:,n),                                              &
     z1_uv,z1_tq,db_surft(:,n),                                                &
     sf_diag,                                                                  &
     cdr10m,sf_diag%cdr10m_n,sf_diag%cd10m_n,chr1p5m(:,n)                      &
     )
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE cable_land_sf_explicit
END MODULE cable_land_sf_explicit_mod
