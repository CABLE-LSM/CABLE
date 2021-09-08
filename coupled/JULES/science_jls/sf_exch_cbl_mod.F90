! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE sf_exch_cbl_mod
 
USE sf_aero_mod, ONLY: sf_aero
USE fcdch_mod, ONLY: fcdch
USE stdev1_mod, ONLY: stdev1
USE sf_rib_mod, ONLY: sf_rib
USE sf_resist_mod, ONLY: sf_resist
 
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_EXCH_MOD'

CONTAINS

!   SUBROUTINE SF_EXCH------------------------------------------------

!  Purpose: Calculate coefficients of turbulent exchange between
!           the surface and the lowest atmospheric layer, and
!           "explicit" fluxes between the surface and this layer.

!  Suitable for Single Column use.

!  Documentation: UM Documentation Paper No 24, section P243.
!                 See especially sub-section (ix).

!---------------------------------------------------------------------

! Arguments :-

SUBROUTINE sf_exch_cbl (                                                          &
 land_pts,nsurft,land_index,                                                  &
 surft_index,surft_pts,flandg,                                                &
 nice, nice_use,                                                              &
 nsnow,ds,hcons_snow,hconsdz_sice,                                            &
 bq_1,bt_1,canhc_surft,canht_pft,lai_pft,canopy,catch,dzsoil,dzsoil_elev,     &
 flake,gc,gc_stom_surft,hcons_soilt,                                          &
 can_model,catch_snow, lq_mix_bl,                                             &
 ho2r2_orog,ice_fract_cat,snowdepth,snow_surft,pstar,qw_1,                    &
 radnet_sea,radnet_sice,radnet_surft,sil_orog,tile_frac,timestep,             &
 surf_hgt,l_elev_absolute_height,emis_surft,emis_soil,                        &
 tl_1,ti,ti_cat,t_soil_soilt,                                                 &
 tsurf_elev_surft,                                                            &
 tsnow,                                                                       &
 tstar_surft,tstar_sea,tstar_sice_cat,z_land,                                 &
 l_ctile,seasalinityfactor,                                                   &
 !IN input data from the wave model
 charnock_w,                                                                  &
 tstar,l_aggregate,l_spec_z0,z0m_scm,z0h_scm,                                 &
 l_aero_classic,l_dust,l_dust_diag,                                           &
 vfrac_surft,vshr_land,vshr_ssi,zh,ddmfx,                                     &
 z0_surft, z0h_surft_bare, z0m_soil_in, z1_uv, z1_uv_top, z1_tq,              &
 z1_tq_top, sf_diag,formdrag,fd_stab_dep,                                     &
 orog_drag_param,z0msea,                                                      &
 lw_down,lw_down_elevcorr_surft,                                              &
 alpha1,alpha1_sea,alpha1_sice,ashtf_prime,ashtf_prime_sea,ashtf_prime_surft, &
 recip_l_mo_sea,cdr10m,                                                       &
 chr1p5m,chr1p5m_ssi_mean,fqw_1,fqw_surft,epot_surft,                         &
 fqw_ice,                                                                     &
 ftl_1,ftl_surft,ftl_ice,fraca,h_blend_orog,charnock,                         &
 rhostar,resfs,resft,rib,                                                     &
 fb_surf,u_s,q1_sd,t1_sd,z0hssi,z0h_surft,                                    &
 z0mssi,z0m_surft,z0m_eff,rho_aresist,aresist,resist_b,                       &
 rho_aresist_surft,aresist_surft,resist_b_surft,                              &
 r_b_dust,cd_std_dust,u_s_std_surft,                                          &
 rhokh_1,rhokh_1_sice_ncats,rhokh_1_sea,rhokm_1,rhokm_land,                   &
 rhokm_ssi,                                                                   &
 dtstar_surft,dtstar_sea,dtstar_sice,rhokh_gb,anthrop_heat,                   &
!CABLE_LSM:{pass additional existing & CABLE state vars
 cycleno, numcycles, sm_levels, sw_surft, co2_mmr, sthu, fland,               &
 curr_day_number,                                                             &
 air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl, ssnow_cbl, bgc_cbl,        & 
 bal_cbl, sum_flux_cbl, veg_cbl,  soilin, soil_cbl )
!CABLE_LSM:}

USE sf_flux_mod, ONLY: sf_flux
!CABLE_LSM:we use a bunch of modules in 10.6. Limit these here - PASS instead  
!Make availabe explicit CALL to CABLE
USE cable_explicit_main_mod, ONLY : cable_explicit_main
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
USE cable_params_mod,         ONLY : soilin_type

USE pftparm, ONLY: emis_pft
USE nvegparm, ONLY: emis_nvg
!CABLE_LSM:}
USE sfl_int_mod, ONLY: sfl_int
USE ice_formdrag_lupkes_mod, ONLY: ice_formdrag_lupkes
USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY: t_i_length,t_j_length

USE ancil_info, ONLY: l_soil_point, l_lice_point, l_lice_surft

USE c_z0h_z0m, ONLY: z0h_z0m, z0h_z0m_classic
USE planet_constants_mod, ONLY: cp, g, r, vkman
USE blend_h, ONLY: lb
USE water_constants_mod, ONLY: lc, rho_ice, rhosea, tm
USE csigma
USE dust_param, ONLY: ndiv, z0_soil
USE c_kappai, ONLY: kappa_seasurf,dzsea

USE bl_option_mod, ONLY: max_stress_grad, on

USE jules_surface_types_mod, ONLY:                                            &
   npft, urban_canyon, urban_roof, soil, lake

USE jules_snow_mod, ONLY: cansnowtile                                         &
                          ,rho_snow_const                                     &
                          ,snow_hcon                                          &
                          ,l_snowdep_surf                                     &
                          ,l_snow_nocan_hc                                    &
                          ,nsmax                                              &
                          ,unload_rate_cnst                                   &
                          ,unload_rate_u

USE jules_mod,   ONLY: snowdep_surft

USE jules_surface_mod, ONLY: i_modiscopt, i_aggregate_opt,                    &
                              cor_mo_iter, use_correct_ustar,                 &
                              iscrntdiag, l_flake_model,                      &
                              l_elev_lw_down,                                 &
                              effective_z0, h_blend_min, ls,                  &
                              IP_ScrnDecpl2, IP_ScrnDecpl3,                   &
                              l_vary_z0m_soil, l_elev_land_ice

USE jules_vegetation_mod, ONLY: l_irrig_dmd, l_vegdrag_surft

USE jules_sea_seaice_mod, ONLY: iseasurfalg, ip_ss_solid,                     &
                                 ip_ss_fixed, ip_ss_surf_div,                 &
                                 ip_ss_surf_div_coupled,                      &
                                 z0miz, z0sice, z0h_z0m_miz,                  &
                                 z0sice, z0h_z0m_sice, z0hsea,                &
                                 emis_sea, emis_sice, hcap_sea,               &
                                 l_iceformdrag_lupkes, beta_evap,             &
                                 ip_hwdrag_limited, ip_hwdrag_reduced_v1,     &
                                 i_high_wind_drag, cdn_max_sea, cdn_hw_sea,   &
                                 u_cdn_max, u_cdn_hw, z_10m

USE jules_internal, ONLY: unload_backgrnd_pft

USE urban_param, ONLY: hgt_gb, hwr_gb, ztm_gb, disp_gb, z0m_mat

USE switches_urban, ONLY: l_moruses_rough, l_moruses_storage

USE sf_diags_mod, ONLY: strnewsfdiag

USE ancil_info, ONLY: ssi_pts, sea_pts, ssi_index, sea_index,                 &
                       fssi_ij, sea_frac, sice_pts_ncat,                      &
                       sice_index_ncat, sice_frac_ncat,                       &
                       l_soil_point, l_lice_point, nsoilt

USE lake_mod, ONLY: lake_t_ice_gb                                             &
                    ,lake_t_snow_gb                                           &
                    ,lake_t_mxl_gb                                            &
                    ,lake_t_mean_gb                                           &
                    ,lake_h_snow_gb                                           &
                    ,lake_h_ice_gb                                            &
                    ,lake_h_mxl_gb                                            &
                    ,lake_depth_gb                                            &
                    ,ts1_lake_gb                                              &
                    ,nusselt_gb                                               &
                    ,g_dt_gb

USE jules_soil_mod, ONLY: hcice, hcwat, hcondeep

USE jules_science_fixes_mod, ONLY: l_fix_ctile_orog, l_fix_wind_snow

USE crop_vars_mod, ONLY: gc_irr_surft, resfs_irr_surft

USE solinc_data, ONLY: sky, l_skyview

USE elevate_mod, ONLY: elevate
USE urbanz0_mod, ONLY: urbanz0
USE sf_orog_mod, ONLY: sf_orog
USE sf_orog_gb_mod, ONLY: sf_orog_gb

USE stochastic_physics_run_mod, ONLY: l_rp2, i_rp_scheme, i_rp2b, z0hm_pft_rp

USE can_drag_mod, ONLY: can_drag_z0, can_drag_phi_m_h

USE qsat_mod, ONLY: qsat_new => qsat,                                         &
                    qsat_mix_new => qsat_mix,                                 &
                    l_new_qsat_jules

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE


INTEGER ::                                                                    &
 land_pts                                                                     &
                       ! IN No of land points being processed.
,nsurft                                                                       &
                       ! IN Number of land tiles per land point.
,land_index(land_pts)                                                         &
                       ! IN Index of land points.
,surft_index(land_pts,nsurft)                                                 &
                       ! IN Index of tile points.
,surft_pts(nsurft)                                                            &
                       ! IN Number of tile points.
,nsnow(land_pts,nsurft)                                                       &
                       ! IN Number of snow layers
,can_model                                                                    &
                       ! IN Switch for thermal vegetation canopy.
,nice                                                                         &
                       ! IN total number of sea ice categories
,nice_use                                                                     &
                       ! IN Number of sea ice categories used
                       !     fully in surface exchange
,formdrag                                                                     &
                       ! IN Switch for orographic form drag
,fd_stab_dep                                                                  &
                       ! IN Switch to implement stability
                       !    dependence of orog form drag
,qsat_tile
                       ! WORK Index of tile being passed into qsat_mix


REAL(KIND=real_jlslsm) ::                                                     &
 bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! IN A buoyancy parameter for lowest atm
                       !    level ("beta-q twiddle").
,bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! IN A buoyancy parameter for lowest atm
                       !    level ("beta-T twiddle").
,canhc_surft(land_pts,nsurft)                                                 &
                       ! IN Areal heat capacity of canopy for
                       !    land tiles (J/K/m2).
,canht_pft(land_pts,npft)                                                     &
                       ! IN Canopy height (m)
,lai_pft(land_pts,npft)                                                       &
                       ! IN Leaf area index
,canopy(land_pts,nsurft)                                                      &
                       ! IN Surface water for land tiles
                       !    (kg/m2).
,catch(land_pts,nsurft)                                                       &
                       ! IN Surface capacity (max. surface water)
                       !    of land tiles (kg/m2).
,catch_snow(land_pts,nsurft)                                                  &
                       ! IN Snow interception capacity of NLT
                       !    tile (kg/m2).
,ds(land_pts,nsurft,nsmax)                                                    &
                       ! IN Snow layer thicknesses (m)
,dzsoil                                                                       &
                       ! IN Soil or land-ice surface layer
                       !    thickness (m).
,dzsoil_elev                                                                  &
                       ! IN Tiled land-ice surface layer
,flake(land_pts,nsurft)                                                       &
                       ! IN Lake fraction.
,gc(land_pts,nsurft)                                                          &
                       ! IN "Stomatal" conductance to evaporation
                       !    for land tiles (m/s).
,gc_stom_surft(land_pts,nsurft)                                               &
                       ! IN canopy conductance to evaporation
,hcons_soilt(land_pts,nsoilt)                                                 &
                       ! IN Soil thermal conductivity including
                       !    effects of water and ice (W/m/K).
,hcons_snow(land_pts,nsurft)                                                  &
                       ! IN Snow thermal conductivity (W/m/K)
,hconsdz_sice(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end,nice_use)                             &
                       ! IN 2 * Sea ice thermal conductivity divided
                       !  by surface layer thickness  (in coupled mode,
                       !  this is passed in from the sea ice model) (W/m2/K)
,ho2r2_orog(land_pts)                                                         &
                       ! IN Peak to trough height of unresolved
                       !    orography divided by 2SQRT(2) (m).
,orog_drag_param                                                              &
                       ! IN Drag coefficient for orographic
                       !    form drag
,ice_fract_cat(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end,nice_use)                            &
                       ! IN Fraction of gridbox which is sea-ice.
                       ! (per category if nice_use>1)
,flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                       ! IN Land fraction on all tiles.
,snowdepth(land_pts,nsurft)                                                   &
                       ! IN Depth of lying snow (m)
,snow_surft(land_pts,nsurft)                                                  &
                       ! IN Lying snow on land tiles (kg/m2).
,pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Surface pressure (Pascals).
,qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! IN Total water content of lowest
                       !    atmospheric layer (kg per kg air).
,radnet_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                       ! IN Open sea net surface radiation (W/m2)
,radnet_sice(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                       ! IN Sea-ice net surface radiation (W/m2)
,radnet_surft(land_pts,nsurft)                                                &
                       ! IN Land tile net surface radiation (W/m2)
,anthrop_heat(land_pts,nsurft)                                                &
                       ! IN Anthropogenic Urban heat source (W/m2)
,sil_orog(land_pts)                                                           &
                       ! IN Silhouette area of unresolved
                       !    orography per unit horizontal area
,tile_frac(land_pts,nsurft)                                                   &
                       ! IN Tile fractions.
,timestep                                                                     &
                       ! IN Timestep in seconds for EPDT calc.
,tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! IN Liquid/frozen water temperature for
                       !    lowest atmospheric layer (K).
,ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                       ! IN Temperature of sea-ice surface layer
                       !    (K)
,ti_cat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)             &
,t_soil_soilt(land_pts,nsoilt)                                                &
                       ! IN Temperature of top soil or land-ice
                       !    layer (K)
,tsurf_elev_surft(land_pts,nsurft)                                            &
                       ! IN Tiled ice sub-surface temperature (K)
,tsnow(land_pts,nsurft,nsmax)                                                 &
                       !  IN Snow layer temperatures (K)
,tstar_surft(land_pts,nsurft)                                                 &
                       ! IN Tile surface temperatures (K).
,tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                       ! IN Open sea surface temperature (K).
,tstar_sice_cat(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end,nice_use)                           &
                       ! IN Sea-ice surface temperature (K).
,z_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                       ! IN Land height (m).
,tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Gridbox Mean Surface Temperature (K)
,vfrac_surft(land_pts,nsurft)                                                 &
                       ! IN Fractional canopy coverage for
                       !    land tiles.
,vshr_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                       ! IN Magnitude of land sfc-to-lowest-level
                       !    wind shear
,vshr_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                       ! IN Mag. of mean sea sfc-to-lowest-level
                       !    wind shear
,zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                       ! IN Height above surface of top of
                       !    boundary layer (metres).
,ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Convective downdraught
                       !    mass-flux at cloud base
,z0_surft(land_pts,nsurft)                                                    &
                       ! IN Tile roughness lengths (m).
,z0h_surft_bare(land_pts,nsurft)                                              &
                       ! IN Tile thermal roughness lengths (m)
                       !    without snow cover
,z0m_soil_in(land_pts)                                                        &
                       ! IN bare soil momentum z0 (m)
,z1_uv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Height of lowest uv level (m).
,z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Height of lowest tq level (m).
                       !    Note, if the grid used is staggered in
                       !    the vertical, Z1_UV and Z1_TQ can be
                       !    different.
,surf_hgt(land_pts,nsurft)                                                    &
                       ! IN Height of elevated tile above
                       !    mean gridbox surface (m)
,emis_surft(land_pts,nsurft)                                                  &
                       ! IN Emissivity for land tiles
,emis_soil(land_pts)                                                          &
                       ! IN Emissivity of underlying soil
,charnock              ! Charnock parameter for sea surface

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                      ! Charnock's coefficient from wave model
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                       ! Height of top of lowest uv-layer
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                       ! Height of top of lowest Tq-layer

REAL(KIND=real_jlslsm) ::                                                     &
 z0h_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                         ! IN Namelist input z0h (if >0)
                         !    (if <=0 use Z0HSEA)
                         !    Used in SCM Configurations
,z0m_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                         ! IN Namelist input z0m (if >0)
                         !    (if <=0 use standard Z0MSEA)
                         !    Used in SCM Configurations

LOGICAL ::                                                                    &
 lq_mix_bl                                                                    &
                       ! IN TRUE if mixing ratios used in
                       !    boundary layer code
,l_aero_classic                                                               &
                       ! IN switch for using CLASSIC aerosol scheme
,l_dust                                                                       &
                       ! IN switch for prognostic mineral dust
,l_dust_diag                                                                  &
                       ! IN switch for diagnostic mineral dust
                       !    lifting
,l_aggregate                                                                  &
                       ! IN Logical to set aggregate surface schem
,l_ctile                                                                      &
                       ! IN switch for coastal tiling
,l_spec_z0                                                                    &
                       ! IN T if using prescribed
                       ! sea surface roughness lengths
,l_elev_absolute_height(nsurft)
                       ! IN switch for tile heights (relative
                       ! to the gridbox mean or absolute)

REAL(KIND=real_jlslsm), INTENT(IN) :: seasalinityfactor                       &
!       Factor allowing for the effect of the salinity of
!       sea water on the evaporative flux.
,lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

!  Modified (INOUT) variables.
!Diagnostics
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

REAL(KIND=real_jlslsm) ::                                                     &
 z0msea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                       ! INOUT Sea-surface roughness length for
                       !       momentum (m).  F617.
,lw_down_elevcorr_surft(land_pts,nsurft)


!  Output variables.

REAL(KIND=real_jlslsm) ::                                                     &
 alpha1(land_pts,nsurft)                                                      &
                       ! OUT Gradients of saturated specific
                       !     humidity with respect to temperature
                       !     between the bottom model layer and
                       !     tile surface
,alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                       ! OUT ALPHA1 for sea.
,alpha1_sice(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                       ! OUT ALPHA1 for sea-ice.
,ashtf_prime(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                       ! OUT Adjusted SEB coefficient for sea-ice
,ashtf_prime_sea(tdims%i_start:tdims%i_end,                                   &
                 tdims%j_start:tdims%j_end)                                   &
                       ! OUT Adjusted SEB coefficient for sea
,ashtf_prime_surft(land_pts,nsurft)                                           &
                       ! OUT Adjusted SEB coefficient for land
                       !     points
,cd_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                       ! OUT Bulk transfer coefficient for
                       !      momentum over sea mean.
,recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
                       ! OUT Reciprocal of the Monin-Obukhov
                       !     length for sea/ice points (m^-1).
,ch_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                       ! OUT Bulk transfer coefficient for heat
                       !    and/or moisture over sea mean.
,cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)          &
                       ! OUT Reqd for calculation of 10m wind
                       !     (u & v).
                       !     NBB: This is output on the UV-grid,
                       !     but with the first and last rows set
                       !     to a "missing data indicator".
                       !     Sea-ice leads ignored.
,chr1p5m(land_pts,nsurft)                                                     &
                       ! OUT Reqd for calculation of 1.5m temp for
                       !     land tiles.
,chr1p5m_ssi_mean(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                       ! OUT CHR1P5M for sea and sea-ice
                       !     (leads ignored).
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! OUT "Explicit" surface flux of QW (i.e.
                       !     evaporation), on P-grid (kg/m2/s).
                       !     for whole grid-box
,fqw_surft(land_pts,nsurft)                                                   &
                       ! OUT Local FQW_1 for land tiles.
,fqw_ice(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end,nice_use)                                  &
                       ! OUT GBM FQW_1 for sea-ice.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! OUT "Explicit" surface flux of TL = H/CP.
                       !     (sensible heat / CP). grid-box mean
,ftl_surft(land_pts,nsurft)                                                   &
                       ! OUT Local FTL_1 for land tiles.
,ftl_ice(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end,nice_use)                                  &
                       ! OUT GBM FTL_1 for sea-ice.
,fraca(land_pts,nsurft)                                                       &
                       ! OUT Fraction of surface moisture flux
                       !     with only aerodynamic resistance
                       !     for land tiles.
,h_blend_orog(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)            &
                       ! OUT Blending height for orographic
                       !     roughness
,rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                       ! OUT Surface air density
,resfs(land_pts,nsurft)                                                       &
                       ! OUT Combined soil, stomatal and
                       !     aerodynamic resistance factor for
                       !     fraction 1-FRACA of land tiles
,resft(land_pts,nsurft)                                                       &
                       ! OUT Total resistance factor
                       !     FRACA+(1-FRACA)*RESFS for snow-free
                       !     tiles, 1 for snow and land-ice.
,rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                       ! OUT Mean bulk Richardson number for
                       !     lowest layer
,fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                       ! OUT Surface flux buoyancy over
                       !     density (m^2/s^3)
,u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                       ! OUT Surface friction velocity (m/s)
,q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! OUT Standard deviation of turbulent
                       !     fluctuations of surface layer
                       !     specific humidity (kg/kg).
,t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! OUT Standard deviation of turbulent
                       !     fluctuations of surface layer
                       !     temperature (K).
,z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                       ! OUT Roughness length for heat and
                       !     moisture over sea/sea ice (m)
,z0h_surft(land_pts,nsurft)                                                   &
                       ! OUT Tile roughness lengths for heat
                       !     and moisture
,z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                       ! OUT Roughness length for momentum
                       !     over sea/sea ice (m)
,z0m_surft(land_pts,nsurft)                                                   &
                       ! OUT Tile roughness lengths for momentum
,z0m_eff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                       ! OUT Effective roughness length for
                       !     momentum
,rho_aresist(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)             &
                       ! OUT RHOSTAR*CD_STD*VSHR
                       !     for CLASSIC aerosol scheme
,aresist(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                       ! OUT 1/(CD_STD*VSHR)
                       !     for CLASSIC aerosol scheme
,resist_b(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                       ! OUT (1/CH-1/CD_STD)/VSHR
                       !     for CLASSIC aerosol scheme
,rho_aresist_surft(land_pts,nsurft)                                           &
                       ! OUT RHOSTAR*CD_STD*VSHR on land tiles
                       !     for CLASSIC aerosol scheme
,aresist_surft(land_pts,nsurft)                                               &
                       ! OUT 1/(CD_STD*VSHR) on land tiles
                       !     for CLASSIC aerosol scheme
,resist_b_surft(land_pts,nsurft)                                              &
                       ! OUT (1/CH-1/CD_STD)/VSHR on land tiles
                       !     for CLASSIC aerosol scheme
,r_b_dust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,ndiv)           &
                       ! OUT surf layer res for mineral dust
,cd_std_dust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)             &
                       ! OUT Bulk transfer coef. for
                       !     momentum, excluding orographic effects
                       !     for mineral dust
,u_s_std_surft(land_pts,nsurft)
                       ! OUT Surface layer scaling velocity
                       !     for tiles excluding orographic
                       !     form drag (m/s) for mineral dust

! Surface exchange coefficients;passed to subroutine IMPL_CAL
REAL(KIND=real_jlslsm) ::                                                     &
 rhokh_1(land_pts,nsurft)                                                     &
!                            ! OUT Surface exchange coefficient for land
!                            !     tiles.
,rhokh_1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)            &
!                            ! OUT Surface exchange coefficient for
!                            !     sea-ice.
,rhokh_1_sice_ncats(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
                                                                nice_use)     &
!                            ! OUT Surface exchange coefficient for
!                            !     sea-ice.
,rhokh_1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)             &
!                            ! OUT Surface exchange coefficient for
!                            !     sea.
,rhokm_1(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end)         &
!                            ! OUT For momentum. NB: This is output on
!                            !     UV-grid, but with the first and last
!                            !     rows set to "missing data indicator".
,rhokm_land(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end)      &
!                            ! OUT For land momentum. NB: This is output
!                            !     on UV-grid, but with the first and
!                            !      last rows set to "missing data".
,rhokm_ssi(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end)       &
!                            ! OUT For mean sea mom. NB: This is output
!                            !     on UV-grid, but with the first and
!                            !     last rows set to "missing data".
,dtstar_surft(land_pts,nsurft)                                                &
!                            ! OUT Change in TSTAR over timestep for
!                            !     land tiles
,dtstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
!                            ! OUT Change is TSTAR over timestep for
!                            !     open sea
,dtstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)    &
!                            ! OUT Change is TSTAR over timestep for
!                            !     sea-ice
,rhokh_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
!                            ! OUT Grid-box surface exchange coefficient
,epot_surft(land_pts,nsurft)
                             ! OUT EPOT for land tiles.

!   Define local storage.

!   (a) Workspace.

REAL(KIND=real_jlslsm) ::                                                     &
 qs1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! Sat. specific humidity
!                                  ! qsat(TL_1,PSTAR)
,rhokm_ssi_nohalo(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
!                                  ! like RHOKM_SSI, but with no halo
,lh0                         ! Latent heat for snow free surface
!                                  !   =LS for sea-ice, =LC otherwise

!  Workspace for sea and sea-ice leads
REAL(KIND=real_jlslsm) ::                                                     &
 cd_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! Drag coefficient
,ch_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! Transfer coefficient for heat and
!                                  ! moisture
,qstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! Surface saturated sp humidity
,rib_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! Bulk Richardson number
,z0h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! Roughness length for heat and
!                                  ! moisture transport
,db_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! Buoyancy difference for sea points
,v_s_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! Surface layer scaling velocity
,hcons_sea(t_i_length * t_j_length)                                           &
                             ! Heat conductivity into sea
,canhc_sea(t_i_length * t_j_length)                                           &
                             ! Heat capacity for sea
,u_s_std_sea(t_i_length * t_j_length)                                         &
                             ! Surface friction velocity for sea
!                                  ! (dummy variable for sea)
,v_s_std_sea(t_i_length * t_j_length)                                         &
                             ! Surface layer scaling velocity
!                                  ! for sea excluding orographic
!                                  ! form drag (m/s).
!                                  ! (dummy variable for sea)
,cd_std_sea(t_i_length * t_j_length)
!                                  ! Local drag coefficient for calc
!                                  ! of interpolation coefficient
!                                  ! (dummy variable for sea)
!                                  ! for sea points (m/s).

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 rhokm_1_sice(:,:)                                                            &
!                                  ! Surface momentum exchange coefficient
!                                  !     for sea-ice.
,rhokm_1_sice_ncats(:,:,:)                                                    &
!                                  ! Surface momentum exchange coefficient
!                                  !      for sea-ice.
,rhokm_1_sea(:,:)                                                             &
!                                  ! Surface momentum exchange coefficient
!                                  !     for sea.
,tau_sea(:,:)                                                                 &
                                   ! GBM tau_1 for sea.
,tau_ice(:,:,:)
                                   ! GBM tau_1 for sea-ice.


!  Workspace for sea-ice and marginal ice zone
REAL(KIND=real_jlslsm) ::                                                     &
 cd_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
                                                        nice_use)             &
                             ! Drag coefficient
,cd_ice_mean(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)             &
                      ! Gridbox aggregate sea ice drag coefficient
,cd_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! Bulk transfer coefficient for
!                                  !      momentum over land.
,ch_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! Bulk transfer coefficient for
!                                  !      het and moisture over land.
,cd_miz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! Drag coefficient
,ch_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
                                                        nice_use)             &
!                             ! Transfer coefficient for heat and
!                             ! moisture
,chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
                                                      nice_use)               &
                       ! CHR1P5M for sea ice, on categories
                       !     (leads ignored).
,chr1p5m_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)             &
                       ! CHR1P5M for sea
                       !     (leads ignored).
,ch_ice_mean(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)             &
!                            ! Gridbox aggregate transfer
!                            ! coefficient for heat and moisture
,ch_miz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! Transfer coefficient for heat and
!                                  ! moisture
,qstar_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! Surface saturated sp humidity
,qstar_ice_cat(tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nice_use)                           &
                       ! Sea-ice surface sat sp humidity on cats
,rib_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
                                                        nice_use)             &
                             ! Bulk Richardson number
,rib_surft(land_pts,nsurft)                                                   &
                       ! RIB for land tiles.
,z0m_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
                                                        nice_use)             &
                             ! Momentum Roughness length.
,z0h_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
                                                        nice_use)             &
                             ! Thermal Roughness length.
,z0m_miz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! Momentum Roughness length.
,z0h_miz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! Thermal Roughness length.
,db_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
                                                        nice_use)             &
                             ! Buoyancy difference for sea ice
,db_ice_mean(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)             &
                             ! Gridbox mean of db_ice
,v_s_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
                                                        nice_use)             &
                             ! Surface layer scaling velocity
!                                  ! for sea ice (m/s).
,v_s_miz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! Surface layer scaling velocity
!                                  ! for marginal sea ice (m/s).
,recip_l_mo_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
                                                           nice_use)          &
!                                  ! Reciprocal of the Monin-Obukhov
!                                  ! length for sea ice (m^-1).
,recip_l_mo_miz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
!                                  ! Reciprocal of the Monin-Obukhov
!                                  ! length for marginal sea ice (m^-1).
,u_s_std_ice(t_i_length * t_j_length,nice_use)                                &
                             ! Surface friction velocity for sea-i
!                                  ! (dummy variable for sea-ice)
,u_s_std_miz(t_i_length * t_j_length)                                         &
                             ! Surface friction velocity for
!                                  ! marginal sea-ice
!                                  ! (dummy variable for marginal sea-ic
,v_s_std_ice(t_i_length * t_j_length,nice_use)                                &
                             ! Surface layer scaling velocity
!                                  ! for sea-ice excluding orographic
!                                  ! form drag (m/s).
!                                  ! (dummy variable for sea-ice)
,v_s_std_miz(t_i_length * t_j_length)                                         &
                             ! Surface layer scaling velocity
!                                  ! for marginal sea-ice excluding
!                                  ! orographic form drag (m/s).
!                                  ! (dummy variable for marginal sea-ic
,cd_std_ice(t_i_length * t_j_length,nice_use)                                 &
                             ! Local drag coefficient for calc
!                                  ! of interpolation coefficient
!                                  ! (dummy variable for sea-ice)
,cd_std_miz(t_i_length * t_j_length)                                          &
                             ! Local drag coefficient for calc
!                                  ! of interpolation coefficient
!                                  ! (dummy variable for marginal sea-ic
,epot_sea(t_i_length * t_j_length)                                            &
                             ! Potential evaporation from
                             ! sea surface
                             ! (dummy variable for sea surface)
,epot_ice(t_i_length * t_j_length)                                            &
                             ! Potential evaporation from sea-ice
                             ! (dummy variable for sea-ice)
! The next fields are temporary and used while the option is still
! available for sea ice exchange coeffs not to be calculated per
! category.  This functionality may be removed in future model versions.
,tstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! Sea ice sfc T averaged over all cats
,ice_fract(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! Sea ice fraction summed over all cats
,sice_frac_use(ssi_pts)
                             ! Sea ice fractions

! The next fields are temporary and used while the option is still
! available for sea ice exchange coeffs not to be calculated per
! category.  This functionality may be removed in future model versions.
INTEGER::                                                                     &
 sice_pts_use,                                                                &
 sice_index_use(ssi_pts)

REAL(KIND=real_jlslsm) ::                                                     &
 z1_tq_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
!                            ! Height of lowest model level
!                            ! relative to sea.
 z1_tq_ctile(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                            ! Height of lowest model level
!                            ! relative to sea dependent on coastal
!                            ! tiling option.
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 z1_tq_top_sea(:,:)                                                           &
                             ! Top of lowest Tq layer over sea
,z1_tq_top_ctile(:,:)                                                         &
                             ! Top of lowest Tq layer over sea
                             ! dependent on coastal tiling option
,chr10m_sice(:,:,:)                                                           &
                       ! CHR10M for sea ice, on categories
                       !     (leads ignored).
,chr10m_sea(:,:)
                       ! CHR10M for sea
                       !     (leads ignored).

!  Workspace for land tiles
REAL(KIND=real_jlslsm) ::                                                     &
 e_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! Evaporation from sea times leads
                       !     fraction (kg/m2/s). Zero over land.
,h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! OUT Surface sensible heat flux over sea
                       !     times leads fraction (W/m2).
                       !     Zero over land.
,cd_std(land_pts,nsurft)                                                      &
                             ! Local drag coefficient for calc
!                                  ! of interpolation coefficient
,cd_surft(land_pts,nsurft)                                                    &
                             ! Drag coefficient
,ch_surft(land_pts,nsurft)                                                    &
                             ! Transfer coefficient for heat and
!                                  ! moisture
,chn(land_pts,nsurft)                                                         &
                             ! Neutral value of CH.
,dq(land_pts)                                                                 &
                             ! Sp humidity difference between
!                                  ! surface and lowest atmospheric lev
,epdt(land_pts)                                                               &
                             ! "Potential" Evaporation * Timestep
,fz0(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! Aggregation function for Z0.
,fz0h(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! Aggregation function for Z0H.
,pstar_land(land_pts)                                                         &
                             ! Surface pressure for land points.
,qstar_surft(land_pts,nsurft)                                                 &
                             !Surface saturated sp humidity.
,rhokh_can(land_pts,nsurft)                                                   &
                             ! Exchange coefficient for canopy
!                                  ! air to surface
,rhokm_1_surft(land_pts,nsurft)                                               &
!                                  ! Momentum exchange coefficient.
,tsurf(land_pts)                                                              &
                             ! Surface layer temp (snow or soil) (
,dzsurf(land_pts)                                                             &
                             ! Surface layer thickness
!                                  ! (snow or soil) (m)
,canhc_surf(land_pts)                                                         &
                             ! Surface layer thickness
!                                  ! (snow or soil) (m)
,dzssi(t_i_length * t_j_length)                                               &
                             ! Surface layer thickness
!                                  ! (sea) (m)
,dzdummy(t_i_length * t_j_length)                                             &
                             ! Dummy field for surface layer thickness
!                                  ! (sea-ice) (m)
,hcons_surf(land_pts)                                                         &
                             ! Thermal conductivity
!                                  ! (snow or soil)  (W/m/K)
,wind_profile_factor(land_pts,nsurft)                                         &
!                                  ! For transforming effective surface
!                                  ! transfer coefficients to those
!                                  ! excluding form drag.
,z0h_eff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                       ! Effective roughness length for heat
,z0m_gb_land(land_pts)                                                        &
                             ! GBM momentum land roughness length
,z0h_gb_land(land_pts)                                                        &
                             ! GBM land roughness length for heat
,z0m_eff_surft(land_pts,nsurft)                                               &
!                                  ! Effective momentum roughness length
,db_surft(land_pts,nsurft)                                                    &
                             ! Buoyancy difference for surface
!                                  ! tile
,v_s_surft(land_pts,nsurft)                                                   &
                             ! Surface layer scaling velocity
!                                  ! for tiles (m/s).
,v_s_std(land_pts,nsurft)                                                     &
!                                  ! Surface layer scaling velocity
!                                  ! for tiles excluding orographic
!                                  ! form drag (m/s).
,u_s_iter_surft(land_pts,nsurft)                                              &
!                                  ! Scaling velocity from middle of
!                                  ! MO scheme - picked up in error by
!                                  ! dust code!
,recip_l_mo_surft(land_pts,nsurft)                                            &
!                                  ! Reciprocal of the Monin-Obukhov
!                                  ! length for tiles (m^-1).
,z0m_soil(land_pts,nsurft)                                                    &
                             ! Bare soil momentum roughness length
                             ! for use in 1 tile dust scheme
,z0h_soil(land_pts,nsurft)                                                    &
                             ! Bare soil roughness length for heat
                             ! for use in 1 tile dust scheme
,wind_profile_fac_soil(land_pts,nsurft)                                       &
                             ! Equivalent of wind_profile_factor
                             ! for use in 1 tile dust scheme
,cd_surft_soil(land_pts,nsurft), ch_surft_soil(land_pts,nsurft)               &
,cd_std_soil(land_pts,nsurft), v_s_surft_soil(land_pts,nsurft)                &
,recip_l_mo_surft_soil(land_pts,nsurft)                                       &
                             ! Dummy output variables from extra
                             ! call to fcdch needed for
                             ! 1 tile dust scheme
,v_s_std_soil(land_pts,nsurft)                                                &
                             ! Bare soil surface layer scaling
                             ! velocity for tiles excluding
                             ! orographic form drag (m/s)
                             ! for use in 1 tile dust scheme
,u_s_iter_soil(land_pts,nsurft)                                               &
                             ! Bare soil scaling velocity from
                             ! middle of MO scheme for use in
                             ! 1 tile dust scheme - picked up in
                             ! error by dust code
,z0h_surft_classic(land_pts,nsurft)                                           &
                             ! z0h to be used in calculation for
                             ! CLASSIC aerosol deposition
,ch_surft_classic(land_pts,nsurft)                                            &
,cd_std_classic(land_pts,nsurft)                                              &
                             ! Additional variables from extra
                             ! call to fcdch needed for aerosol
                             ! deposition with different z0h
,cd_surft_classic(land_pts,nsurft)                                            &
,v_s_surft_classic(land_pts,nsurft)                                           &
,recip_l_mo_surft_classic(land_pts,nsurft)                                    &
,v_s_std_classic(land_pts,nsurft)                                             &
,u_s_iter_classic(land_pts,nsurft)                                            &
                             ! Dummy output variables from extra
                             ! call to fcdch needed for aerosol
                             ! deposition with different z0h
,t_elev(land_pts,nsurft)                                                      &
                             ! Temperature at elevated height (k)
,q_elev(land_pts,nsurft)                                                      &
                             ! Specific humidity at elevated
!                                  !     height (kg per kg air)
,qs1_elev(land_pts,nsurft)                                                    &
                             ! Saturated specific humidity at elev
!                                  !     height (kg per kg air)
,scaling_urban(land_pts,nsurft)                                               &
                             ! MORUSES: ground heat flux scaling;
                             ! canyon tile only coupled to soil
,hcon_lake(land_pts)                                                          &
!                            ! "average" thermal conductivity of
!                            ! lake-ice, lake and soil sandwich (W/m/K).
,zdt_surft(land_pts,nsurft)                                                   &
                             ! Difference between the canopy height and
                             ! displacement height (m)
,phi_m(land_pts)                                                              &
                             ! Monin-Obukhov stability function for momentum
                             ! integrated to the model's lowest wind level.
,phi_h(land_pts)
                             ! Monin-Obukhov stability function for scalars
                             ! integrated to the model's lowest temperature
                             ! and humidity level.

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 tau_surft(:,:)
                             ! Local tau_1 for land tiles.

! dummy arrays required for sea and se-ice to create universal
! routines for all surfaces
REAL(KIND=real_jlslsm) ::                                                     &
 array_zero(t_i_length * t_j_length)                                          &
                                ! Array of zeros
,array_one(t_i_length * t_j_length)                                           &
                                ! Array of ones

,array_emis(t_i_length * t_j_length)                                          &
                                ! Emissivity expanded to an array
,zdt_dummy(t_i_length * t_j_length)
                                ! Dummy array for zdt

!Gridbox mean values calculated from soil tiled versions for FLAKE
REAL(KIND=real_jlslsm) ::                                                     &
  hcons_mean_soil(land_pts),                                                  &
  tsoil_mean_soil(land_pts)

LOGICAL ::                                                                    &
 array_false(t_i_length * t_j_length)
                                ! Array of .FALSE.
INTEGER ::                                                                    &
 array_zero_int(t_i_length * t_j_length)    ! Array of zeros


!   (b) Scalars.

INTEGER ::                                                                    &
 i,j                                                                          &
             ! Loop counter (horizontal field index).
,k                                                                            &
             ! Loop counter (tile field index).
,l                                                                            &
             ! Loop counter (land point field index).
,n                                                                            &
             ! Loop counter (tile index).
,m                                                                            &
              ! Index for soil tile
,jits                                                                         &
             ! Counter for iteration for Z0H
,first_counter                                                                &
             ! Used for setting sea ice fields
,n_veg
             ! Actual or dummy pointer to array
             ! defined only on PFTs
REAL(KIND=real_jlslsm) ::                                                     &
 tau                                                                          &
             ! Magnitude of surface wind stress over sea.
,zetam                                                                        &
             ! Temporary in calculation of CHN.
,zetah                                                                        &
             ! Temporary in calculation of CHN.
,zeta1                                                                        &
             ! Work space
,z0                                                                           &
             ! yet more workspace
,ustr_l                                                                       &
             ! Low-wind estimate of friction velocity
,ustr_n                                                                       &
             ! Neutral estimate of friction velocity
,tol_ustr_n                                                                   &
             ! Tolerance for USTR_N
,tol_ustr_l                                                                   &
             ! Tolerance for USTR_L (see below)
,bl_stress_grad                                                               &
             ! Stress gradient across boundary layer
,ws10                                                                         &
             ! 10-m wind speed, used in calculation of unloading of snow from
             ! vegetation (m s-1)

             ! Temporary variables for adjustment of downwelling
             ! logwave to elevation tiles and correction back to
             ! conserve gridbox mean
,t_rad                                                                        &
,lw_down_surftsum                                                             &
,lw_down_surftabs                                                             &
,z0msea_max                                                                   &
             ! Sea roughness at maximum neutral drag
,u10n                                                                         &
             ! Neutral wind speed at 10 m.
,cdn_lim_loc
             ! Local limiting value of the neutral drag coefficient


LOGICAL :: l_cdr10m_snow
             ! Flag indicating if cdr10m (an interpolation coefficient) is
             ! to be calculated for use with snow unloading.
LOGICAL :: l_vegdrag_active_here
             ! Logical to indicate whether the vegetative drag scheme
             ! is active on the current surface tile in cases where 
             ! l_vegdrag_surft itself may not be applicable
LOGICAL, PARAMETER :: l_vegdrag_ssi = .FALSE.
             ! Logical to indicate that the canopy drag scheme cannot
             ! be applied at sea or sea-ice points

REAL(KIND=real_jlslsm) :: sea_point

INTEGER :: n_diag

INTEGER :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_EXCH'

!CABLE_LSM: Dec. vars passed from sf_expl for CABLE
integer :: mype, timestep_number, endstep,  cycleno, numcycles,              &
  sm_levels         ! # of soil layers 

integer :: curr_day_number

REAL :: fland(land_pts)

TYPE(air_type),       INTENT(inout)  :: air_cbl
TYPE(met_type),       INTENT(inout)  :: met_cbl
TYPE(radiation_type), INTENT(inout)  :: rad_cbl
TYPE(roughness_type), INTENT(inout)  :: rough_cbl
TYPE(canopy_type),    INTENT(inout)  :: canopy_cbl
TYPE(soil_snow_type), INTENT(inout)  :: ssnow_cbl
TYPE(bgc_pool_type),  INTENT(inout)  :: bgc_cbl
TYPE(balances_type),  INTENT(inout)  :: bal_cbl
TYPE(sum_flux_type),  INTENT(inout)  :: sum_flux_cbl
TYPE(veg_parameter_type),  INTENT(inout) :: veg_cbl
TYPE(soilin_type),         INTENT(inout) :: soilin  
TYPE(soil_parameter_type), INTENT(inout) :: soil_cbl

REAL,  DIMENSION( tdims%i_end,tdims%j_end ) ::                             &
  true_latitude,   &
  true_longitude

!___UM soil/snow/radiation/met vars
REAL,  DIMENSION(land_pts) :: & 
  bexp_gb,    & ! => parameter b in Campbell equation 
  hcon_gb,    & ! 
  satcon_gb,  & ! hydraulic conductivity @ saturation [mm/s]
  sathh_gb,   &
  smvcst_gb,  &
  smvcwt_gb,  &
  smvccl_gb,  &
  soil_alb

REAL :: sw_surft(land_pts,nsurft)

REAL,  DIMENSION( tdims%i_end,tdims%j_end ) ::                             &
  ls_rain_cable,    &
  ls_snow_cable

REAL,  DIMENSION( tdims%i_end,tdims%j_end ) ::                             &
  cosz_gb, &
  sin_theta_latitude

REAL,  DIMENSION(sm_levels) :: ardzsoil

REAL :: co2_mmr

REAL, DIMENSION(land_pts, sm_levels) ::                         &
  sthu 

!Declare these locally. ACCESS1.3 uses ashtf, ashtf_tile
!and local decs of these vars in SEB factor calc. 10.6 does it this way 
Real, parameter :: fEpsilon = 0.62198
Real ::                                                                        &
  fc_virtual,                                                               &
  fD_T,                                                                        &
  fDS_RATIO,                                                                   &
  fLH
real :: rhokpm(land_pts,nsurft)    
logical, save :: first_call = .true.

!CABLE_LSM: End 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!CABLE_LSM:
fc_virtual =  1. / fepsilon - 1.

!-----------------------------------------------------------------------
!  0. Initialisations
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(I,j) &
!$OMP SHARED(t_i_length,t_j_length,array_zero,array_one,array_false,         &
!$OMP array_zero_int,tdims,rhokm_land,rhostar,pstar,tstar,r,cd_land,         &
!$OMP ch_land,ftl_1,fqw_1,q1_sd,t1_sd,cdr10m,zdt_dummy)
!$OMP DO SCHEDULE(STATIC)
DO i = 1,t_i_length * t_j_length
  array_zero(i)     = 0.0
  array_one(i)      = 1.0
  array_false(i)    = .FALSE.
  array_zero_int(i) = 0
  zdt_dummy(i)      = 0.0
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    rhokm_land(i,j) = 0.0
    rhostar(   i,j) = pstar(i,j) / ( r * tstar(i,j) )
    !                        ... surface air density from ideal gas equation
    cd_land(i,j) = 0.0
    ch_land(i,j) = 0.0
    ftl_1(i,j) = 0.0
    fqw_1(i,j) = 0.0
    q1_sd(i,j) = 0.0
    t1_sd(i,j) = 0.0
    cdr10m(i,j) = 0.0
  END DO
END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL 


IF (sf_diag%l_tau_1) THEN
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      sf_diag%tau_1(i,j) = 0.0
    END DO
  END DO
END IF

IF (sf_diag%l_tau_surft .OR. sf_diag%l_tau_1) THEN
  ALLOCATE(rhokm_1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end))
  ALLOCATE(rhokm_1_sice_ncats(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use))
  ALLOCATE(rhokm_1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end))
  ALLOCATE(tau_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end))
  ALLOCATE(tau_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use))
  ALLOCATE(tau_surft(land_pts,nsurft))
ELSE
  ALLOCATE(rhokm_1_sice(1,1))
  ALLOCATE(rhokm_1_sice_ncats(1,1,nice_use))
  ALLOCATE(rhokm_1_sea(1,1))
  ALLOCATE(tau_sea(1,1))
  ALLOCATE(tau_ice(1,1,nice_use))
  ALLOCATE(tau_surft(1,nsurft))
END IF



IF ( l_rp2 .AND. i_rp_scheme == i_rp2b) THEN
  DO n = 1,npft
    z0h_z0m(n) = z0hm_pft_rp(n)
  END DO
END IF
!-----------------------------------------------------------------------
!  1. Initialise FTL_SURFT and RIB_SURFT on all tiles at all points,
!     to allow STASH to process these as diagnostics.
!-----------------------------------------------------------------------
!$OMP PARALLEL PRIVATE(i,j,l,n) DEFAULT(NONE)                                 &
!$OMP SHARED(land_pts,nsurft,ftl_surft,fqw_surft,tau_surft,                   &
!$OMP rib_surft,z0m_surft,u_s_std_surft,chr1p5m,resfs,snowdep_surft,          &
!$OMP snowdepth,sf_diag)                                                      &
!$OMP SHARED(land_index, flandg,surf_hgt,z_land,t_i_length,can_model,         &
!$OMP cansnowtile,l_snowdep_surf,snow_surft,rho_snow_const,                   &
!$OMP   l_fix_ctile_orog, l_ctile, canhc_surft, vfrac_surft, tile_frac,       &
!$OMP   cd_surft, ch_surft, z0h_surft, z0m_eff_surft, rhokpm,                 &
!$OMP   radnet_surft, fraca, resft)    
DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    ! MORUSES Initialise urban roughness array
    ftl_surft(l,n) = 0.0
    fqw_surft(l,n) = 0.0
    rib_surft(l,n) = 0.0
    z0m_surft(l,n) = 0.0
    u_s_std_surft(l,n) = 0.0
    chr1p5m(l,n) = 0.0
    resfs(l,n) = 0.0
    ! Equivalent snowdepth for surface calculations.
    snowdep_surft(l,n) = snowdepth(l,n)
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
      tau_surft(l,n) = 0.0
      !CABLE_LSM:
      canhc_surft(l,n) = 0.0
      vfrac_surft(l,n) = 0.0
      if (tile_frac(l,n).eq.0.0) then
        cd_surft(l,n) = 0.0
        ch_surft(l,n) = 0.0
        z0h_surft(l,n) = 0.0
        z0m_eff_surft(l,n) = 0.0
        rhokpm(l,n) = 0.0
        radnet_surft(l,n) = 0.0
        fraca(l,n) = 0.0
        resfs(l,n) = 0.0
        resft(l,n) = 0.0
      end if
    END DO
!$OMP END DO NOWAIT
  END DO
END IF


DO n = 1,nsurft
  IF ( (can_model == 4) .AND. cansnowtile(n) .AND. l_snowdep_surf ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,land_pts
      snowdep_surft(l,n) = snow_surft(l,n) / rho_snow_const
    END DO
!$OMP END DO NOWAIT
  END IF
END DO
!-----------------------------------------------------------------------
!  2.  Calculate QSAT values required later.
!-----------------------------------------------------------------------


! Calculate temeprature and specific humidity for elevation bands
IF (l_fix_ctile_orog .AND. l_ctile) THEN
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    IF (flandg(i,j) > 0.0 .AND. flandg(i,j) < 1.0) THEN
      ! calculate height of orography relative to grid-box mean
      ! limit this to 1000m. z_land already limited to be > 0
      ! in ni_bl_ctl
      surf_hgt(l,:) = MIN(z_land(i,j) * (1.0 / flandg(i,j) - 1.0),1000.0)
    END IF
  END DO
!$OMP END DO NOWAIT
END IF
!$OMP END PARALLEL


CALL elevate(                                                                 &
 land_pts,nsurft,surft_pts,land_index,surft_index,                            &
 tl_1,qw_1,pstar,surf_hgt,l_elev_absolute_height,z_land,                      &
 t_elev,q_elev)

! DEPENDS ON: qsat_mix
qsat_tile = 0

IF (l_new_qsat_jules) THEN
  IF (lq_mix_bl) THEN
    CALL qsat_mix_new(qs1,tl_1,pstar,t_i_length,t_j_length)
  ELSE
    CALL qsat_new(qs1,tl_1,pstar,t_i_length,t_j_length)
  END IF
ELSE
  CALL qsat_mix(qs1,tl_1,pstar,t_i_length * t_j_length,lq_mix_bl)
END IF

!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(i, j, l)             &
!$OMP SHARED(land_index, land_pts, pstar, pstar_land, t_i_length)             &
!$OMP SCHEDULE(STATIC)
DO l = 1,land_pts
  j=(land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length
  pstar_land(l) = pstar(i,j)
END DO
!$OMP END PARALLEL DO

IF (l_new_qsat_jules) THEN
  IF (lq_mix_bl) THEN
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(n)                     &
!$OMP SHARED(nsurft,qstar_surft,tstar_surft,pstar_land,land_pts,              &
!$OMP qs1_elev,t_elev)  SCHEDULE(STATIC)
    DO n = 1,nsurft
      CALL qsat_mix_new(qstar_surft(:,n),tstar_surft(:,n),pstar_land,land_pts)
      CALL qsat_mix_new(qs1_elev(:,n),t_elev(:,n),pstar_land,land_pts)
    END DO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(n)                     &
!$OMP SHARED(nsurft,qstar_surft,tstar_surft,pstar_land,land_pts,              &
!$OMP qs1_elev,t_elev)  SCHEDULE(STATIC)
    DO n = 1,nsurft
      CALL qsat_new(qstar_surft(:,n),tstar_surft(:,n),pstar_land,land_pts)
      CALL qsat_new(qs1_elev(:,n),t_elev(:,n),pstar_land,land_pts)
    END DO
!$OMP END PARALLEL DO
  END IF
ELSE
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(n)                     &
!$OMP SHARED(nsurft,qstar_surft,tstar_surft,pstar_land,land_pts,              &
!$OMP lq_mix_bl,qs1_elev,t_elev)  SCHEDULE(STATIC)
  DO n = 1,nsurft
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(qstar_surft(:,n),tstar_surft(:,n),                          &
                  pstar_land,land_pts,lq_mix_bl)
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(qs1_elev(:,n),t_elev(:,n),                                  &
                  pstar_land,land_pts,lq_mix_bl)
  END DO
!$OMP END PARALLEL DO
END IF

!-----------------------------------------------------------------------
!  3. Calculation of transfer coefficients and surface layer stability
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  3.1 Calculate neutral roughness lengths
!-----------------------------------------------------------------------

! Land tiles
! Z0_SURFT contains the appropriate value for land-ice points, but has to
! be modified for snow-cover on non-land-ice points. Thermal roughness
! lengths are set to be proportional to the tiled roughness length in
! the case of multiple tiles, but if tiled properties have already
! been aggregated, an adjustment for snow cover is required. In this case
! a ratio of 0.1 between the thermal and momentum roughness lengths over
! snow is assumed; to do otherwise would require reaggregation. In the
! case of multiple tiles, the assignment is delayed until the urban
! options have been considered.
!
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n, z0, zeta1)   &
!$OMP          SHARED(nsurft, surft_pts, surft_index, snow_surft,            &
!$OMP                 l_soil_point, z0_surft, snowdep_surft, l_aggregate,    &
!$OMP                 i_aggregate_opt, z0h_surft_bare, z0m_surft,            &
!$OMP                 l_moruses_rough, urban_canyon, urban_roof,              &
!$OMP                 ztm_gb, z0h_z0m, db_surft, z0h_surft,                  &
!$OMP                 wind_profile_factor, z0m_eff_surft, z0h_surft_classic, &
!$OMP                 z0h_z0m_classic)   SCHEDULE(STATIC)
DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    IF ( snow_surft(l,n) > 0.0 .AND. l_soil_point(l) ) THEN
      z0 = z0_surft(l,n) - 0.1 * snowdep_surft(l,n)
      zeta1 = MIN( 5.0e-4 , z0_surft(l,n)  )
      z0m_surft(l,n) = MAX( zeta1 , z0 )
      !     Set z0h_surft explicitly if this option is selected,
      !     otherwise, it will be set for the first tile below.
      IF (l_aggregate .AND. i_aggregate_opt == 1) THEN
        z0 = z0h_surft_bare(l,n) - 0.1 * 0.1 * snowdep_surft(l,n)
        zeta1 = MIN( 5.0e-5 , z0h_surft_bare(l,n)  )
        z0h_surft(l,n) = MAX( zeta1 , z0 )
      END IF
    ELSE
      z0m_surft(l,n) = z0_surft(l,n)
      !     Set z0h_surft explicitly if this option is selected,
      !     otherwise, it will be set for the first tile below.
      IF (l_aggregate .AND. i_aggregate_opt == 1)                             &
        z0h_surft(l,n) = z0h_surft_bare(l,n)
    END IF

    ! MORUSES: z0m_surft (z0_surft previously set in sparm) reset with ztm_gb
    ! (calculated in init_urban) to remove effects of snow. MORUSES does not
    ! yet contain a parametrisation for snow although snow should not affect
    ! the behaviour of bluff bodies rather it affects the  material roughness
    ! of the road & roof
    ! Note: The test on smvcst (not replaced with l_soil_point) is used to alter
    !       just the roof tile and not the
    !       land-ice tiles as the tile is shared at the moment. This will not
    !       be required when flexible tiles are introduced.

    IF ( .NOT. l_aggregate                                                    &
       .AND. l_moruses_rough                                                  &
       .AND. l_soil_point(l)                                                  &
       .AND. ( n == urban_canyon .OR. n == urban_roof ) ) THEN
      ! Snow could be added here to the material roughness length for
      ! momentum before passing to urbanz0. Only the road & roof should be
      ! affected as the walls will be essentially snow-free
      z0m_surft(l,n) = ztm_gb(l)
    END IF

    !
    !   Set the thermal roughness length if aggregation is not being
    !   carried out, or if the original scheme is being used.
    !   It must be done here for consistency with the urban options.
    IF ( ( .NOT. l_aggregate) .OR.                                            &
         (l_aggregate .AND. i_aggregate_opt == 0) )                           &
      z0h_surft(l,n) = z0h_z0m(n) * z0m_surft(l,n)
    db_surft(l,n) = 0.0
    wind_profile_factor(l,n) = 1.0
    z0m_eff_surft(l,n) = z0m_surft(l,n)

    ! Also set additional roughness length for use in CLASSIC aerosol
    ! deposition scheme
    z0h_surft_classic(l,n) = z0h_z0m_classic(n) * z0m_surft(l,n)
  END DO
END DO
!$OMP END PARALLEL DO

! MORUSES: urbanz0 calculates r.l. for heat and updates z0h_surft, which is
! previously set above. Two calls are required; one for urban_canyon and one
! for urban_roof. It is done this way to remove the need for do loop over tile
! type and to avoid over-writing any land-ice tiles. If l_aggregate roughness
! lengths are already calculated: z0m in sparm & z0h above.
IF ( .NOT. l_aggregate ) THEN
  IF ( l_moruses_rough ) THEN
    n = urban_canyon
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,l,j,i) &
!$OMP SHARED(surft_pts,n,surft_index,land_index,t_i_length,z1_uv,z1_tq,hgt_gb, &
!$OMP hwr_gb,disp_gb,ztm_gb,z0h_surft,urban_roof,                             &
!$OMP z0h_surft_classic)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j = ( land_index(l) - 1 ) / t_i_length + 1
      i = land_index(l) - ( j - 1 ) * t_i_length
      CALL urbanz0(                                                           &
         n, z1_uv(i,j), z1_tq(i,j), hgt_gb(l), hwr_gb(l), disp_gb(l),         &
         z0m_mat, ztm_gb(l),                                                  &
         z0h_surft(l,n) )
      CALL urbanz0(                                                           &
         urban_roof, z1_uv(i,j), z1_tq(i,j), hgt_gb(l), hwr_gb(l), disp_gb(l),&
         z0m_mat, ztm_gb(l),                                                  &
         z0h_surft(l,urban_roof) )
      ! Make CLASSIC aerosol roughness length for urban tiles consistent
      ! with those for heat and momentum
      z0h_surft_classic(l,n) = z0h_surft(l,n)
      z0h_surft_classic(l,urban_roof) = z0h_surft(l,urban_roof)
    END DO
!$OMP END PARALLEL DO

  END IF

END IF


! Calculate roughness length affected by roughness sublayer in neutral
! conditions.
zdt_surft(:,:) = 0.0
DO n = 1,nsurft
  IF (l_vegdrag_surft(n)) THEN
    CALL can_drag_z0(                                                         &
      land_pts, surft_pts(n), surft_index(:,n),                               &
      array_zero, canht_pft(:,n), lai_pft(:,n),                               &
      z0m_surft(:,n), z0h_surft(:,n), zdt_surft(:,n))
  END IF
END DO

! Calculate orographic effective parameter for neutral conditions
! if using orographic roughness scheme
IF (formdrag ==  effective_z0) THEN
  DO n = 1,nsurft
    CALL sf_orog (                                                            &
     land_pts,surft_pts(n),land_index,surft_index(:,n),                       &
     fd_stab_dep,orog_drag_param,                                             &
     ho2r2_orog,rib_surft(:,n),sil_orog,z0m_surft(:,n),z1_uv,                 &
     wind_profile_factor(:,n),z0m_eff_surft(:,n)                              &
     )
  END DO
END IF

!CABLE_LSM: CABLE explicit call
  CALL cable_explicit_main(                                                    &
            !mype, timestep, timestep_number, endstep, cycleno, numcycles,      &
            timestep, cycleno, numcycles,      &
            !tdims%i_end,tdims%j_end, land_pts, nsurft, npft,                   &
            land_pts, nsurft, npft,                   &
            sm_levels,                                                         &
            !true_latitude, true_longitude,                                     &
            land_index, tile_frac, surft_pts, surft_index,                     &
            !bexp_gb, hcon_gb, satcon_gb, sathh_gb,                             &
            !smvcst_gb, smvcwt_gb, smvccl_gb, soil_alb,                         & 
            snow_surft, lw_down,                                               &
            !cosz_gb, 
            sw_surft, & 
            !ls_rain_cable, ls_snow_cable,         &
            tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy, Fland,         &
            !passed
            co2_mmr, sthu, canht_pft, lai_pft,                                 &
            sin_theta_latitude, &!ardzsoil,                                       &
            ftl_surft, fqw_surft,                                              &
            tstar_surft, u_s, u_s_std_surft, cd_surft, ch_surft,               &
            radnet_surft, fraca, resfs, resft, z0h_surft, z0m_surft,           &
            recip_l_MO_surft, epot_surft, curr_day_number,                     &
            air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                  &
            ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                &
            soilin, soil_cbl )

  CD_STD = CD_surft
  V_S_surft = U_S_STD_surft
  V_S_STD = U_S_STD_surft
!CABLE_LSM: End
  z0m_eff_surft = z0m_surft

!-----------------------------------------------------------------------
! Calculate RESFT with neutral CH and EPDT = 0 for use in calculation
! of Richardson number. RESFT=1 for snow and land-ice.
!-----------------------------------------------------------------------

!CABLE_LSM:clobbered flake 
flake =0.0
!CABLE_LSM:mod OMP directive below
DO n = 1,nsurft
  IF (l_vegdrag_surft(n)) THEN
    CALL can_drag_phi_m_h(                                                    &
      land_pts, surft_pts(n), surft_index(:,n), land_index,                   &
      array_zero, z1_uv, z1_tq, canht_pft(:,n), lai_pft(:,n),                 &
      z0m_surft(:,n), z0h_surft(:,n),                                         &
      phi_m, phi_h)

!$OMP PARALLEL DO IF(surft_pts(n)>1) DEFAULT(NONE)                           &
!$OMP PRIVATE(i, j, k, l)                                                    &
!$OMP SHARED(surft_pts, surft_index, t_i_length, land_index,                 &
!$OMP        phi_m, phi_h, chn, wind_profile_factor, dq, qw_1,               &
!$OMP        qstar_surft, epdt, n, resft, resfs, flake, fraca) SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      chn(l,n) = vkman**2 / (phi_m(l) * phi_h(l)) * wind_profile_factor(l,n)
      dq(l) = qw_1(i,j) - qstar_surft(l,n)
      epdt(l) = 0.0
      !CABLE_LSM: from 8.5 - necessary on both sides of IF condition?
      RESFT(L,N) =  min(1.,FLAKE(L,N) + (1. - FLAKE(L,N)) *        &
                      ( FRACA(L,N) + (1. - FRACA(L,N))*RESFS(L,N) )) 
    END DO
!$OMP END PARALLEL DO

  ELSE ! l_vegdrag_surft = F

!$OMP PARALLEL DO IF(surft_pts(n)>1) DEFAULT(NONE)                            &
!$OMP PRIVATE(i, j, k, l, zetah, zetam)                                       &
!$OMP SHARED(surft_pts, surft_index, t_i_length, land_index, z1_uv,           &
!$OMP       z0m_surft, z1_tq, z0h_surft, chn, wind_profile_factor, dq, qw_1,  &
!$OMP        qstar_surft, epdt, n, resft, flake, fraca, resfs) SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      zetam = LOG ( (z1_uv(i,j) + z0m_surft(l,n)) / z0m_surft(l,n) )
      zetah = LOG ( (z1_tq(i,j) + z0m_surft(l,n)) / z0h_surft(l,n) )
      chn(l,n) = (vkman / zetah) * (vkman / zetam) *                          &
        wind_profile_factor(l,n)
      dq(l) = qw_1(i,j) - qstar_surft(l,n)
      epdt(l) = 0.0
      !CABLE_LSM: from 8.5 - necessary on both sides of IF condition?
      RESFT(L,N) =  min(1.,FLAKE(L,N) + (1. - FLAKE(L,N)) *        &
                      ( FRACA(L,N) + (1. - FRACA(L,N))*RESFS(L,N) )) 
    END DO
!$OMP END PARALLEL DO
  END IF

  ! We should only attempt to access sf_diag%resfs_stom(:,n) if it has
  ! been fully allocated.
  IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
    n_diag = n
  ELSE
    n_diag = 1
  END IF

!CABLE_LSM:{ 
!C!  CALL sf_resist (                                                            &
!C!   land_pts,surft_pts(n),land_index,surft_index(:,n),                         &
!C!   canopy(:,n),catch(:,n),chn(:,n),dq,epdt,flake(:,n),gc(:,n),                &
!C!   gc_stom_surft(:,n),snowdep_surft(:,n),vshr_land,fraca(:,n),resfs(:,n),     &
!C!   resft(:,n),gc_irr_surft(:,n),resfs_irr_surft(:,n),                         &
!C!   sf_diag%resfs_stom(:,n_diag),sf_diag%l_et_stom, sf_diag%l_et_stom_surft)  

END DO

! RESFT < 1 for snow on canopy if canopy snow model used
! N.B. chn is calculated in the preceding loop over tiles and
! used in the next loop over npft. This works only if the
! first npft tiles are the vegetated ones.
IF ( .NOT. l_aggregate .AND. can_model == 4) THEN
  DO n = 1,npft
    IF ( cansnowtile(n) ) THEN
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)     &
!$OMP          SHARED(surft_pts, surft_index, land_index, t_i_length,        &
!$OMP                 snow_surft, gc, catch_snow, tstar_surft, vshr_land,    &
!$OMP                 fraca, resfs, chn, resft, l_irrig_dmd,                 &
!$OMP                 resfs_irr_surft, n)   SCHEDULE(STATIC)
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j-1) * t_i_length
        IF (snow_surft(l,n) >  0.0) THEN
          gc(l,n) = 0.06 * snow_surft(l,n)**0.6 * catch_snow(l,n)**0.4        &
                       * 2.06e-5 * (tm / tstar_surft(l,n))**1.75              &
                       * (1.79+3 * SQRT(vshr_land(i,j)))                      &
                       / (2 * rho_ice * 5.0e-4**2)
          fraca(l,n) = 0.0
          resfs(l,n) = gc(l,n) / (gc(l,n) + chn(l,n) * vshr_land(i,j))
          resft(l,n) = resfs(l,n)
          IF (l_irrig_dmd) THEN
            resfs_irr_surft(l,n) = resfs(l,n)
            ! print*,'changing resfs cos of canopy snow?'
          END IF
        END IF
      END DO
!$OMP END PARALLEL DO
    END IF
  END DO
END IF

!-----------------------------------------------------------------------
!  3.2 Calculate bulk Richardson number for the lowest model level.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  CALL sf_rib (                                                               &
   land_pts,surft_pts(n),land_index,surft_index(:,n),                         &
   bq_1,bt_1,qstar_surft(:,n),q_elev(:,n),resft(:,n),t_elev(:,n),             &
   tstar_surft(:,n),vshr_land,z0h_surft(:,n),z0m_surft(:,n),zdt_surft(:,n),   &
   z1_tq,z1_uv,l_vegdrag_surft(n),                                            &
   rib_surft(:,n),db_surft(:,n)                                               &
   )
END DO

!-----------------------------------------------------------------------
!  3.3 Calculate stability corrected effective roughness length.
!  Stability correction only applies to land points.
!-----------------------------------------------------------------------

IF (formdrag ==  effective_z0) THEN
  DO n = 1,nsurft
    CALL sf_orog (                                                            &
     land_pts,surft_pts(n),land_index,surft_index(:,n),                       &
     fd_stab_dep,orog_drag_param,                                             &
     ho2r2_orog,rib_surft(:,n),sil_orog,z0m_surft(:,n),z1_uv,                 &
     wind_profile_factor(:,n),z0m_eff_surft(:,n)                              &
     )
  END DO
END IF

!-----------------------------------------------------------------------
!  3.4 Calculate CD, CH via routine FCDCH.
!-----------------------------------------------------------------------

!CABLE_LSM:{follows 10.6, however std JLS calls slightly modified here CHECK 
!C!! Land tiles
!C!DO n = 1,nsurft
!C!  IF (l_vegdrag_surft(n)) THEN
!C!    n_veg = n
!C!    z0m_eff_surft(:,n) = z0m_surft(:,n)
!C!  ELSE
!C!    n_veg = 1
!C!  END IF
!C!  CALL fcdch (                                                                &
!C!   cor_mo_iter,land_pts,surft_pts(n),                                         &
!C!   surft_index(:,n),land_index,                                               &
!C!   db_surft(:,n),vshr_land,                                                   &
!C!   z0m_eff_surft(:,n),z0h_surft(:,n),zdt_surft(:,n),zh,                       &
!C!   z1_uv,z1_uv_top,z1_tq,z1_tq_top,wind_profile_factor(:,n),                  &
!C!   ddmfx,ip_ss_solid,charnock,                                                &
!C!   charnock_w,                                                                &
!C!   l_vegdrag_surft(n),canht_pft(:,n_veg),lai_pft(:,n_veg),                    &
!C!   cd_surft(:,n),ch_surft(:,n),cd_std(:,n),                                   &
!C!   v_s_surft(:,n),v_s_std(:,n),recip_l_mo_surft(:,n),                         &
!C!   u_s_iter_surft(:,n)                                                        &
!C!   )
!C!END DO
!C!
!C!! As roughness length have been changed by vegetation drag effect, effective
!C!! roughness length should be updated.
!C!DO n = 1,nsurft
!C!  IF (l_vegdrag_surft(n)) THEN
!C!    z0m_surft(:,n) = z0m_eff_surft(:,n)
!C!    IF (formdrag == effective_z0) THEN
!C!      CALL sf_orog (                                                          &
!C!       land_pts,surft_pts(n),land_index,surft_index(:,n),                     &
!C!       fd_stab_dep,orog_drag_param,                                           &
!C!       ho2r2_orog,rib_surft(:,n),sil_orog,z0m_surft(:,n),z1_uv,               &
!C!       wind_profile_factor(:,n),z0m_eff_surft(:,n)                            &
!C!       )
!C!    END IF
!C!  END IF
!C!END DO
!C!
!C!!$OMP PARALLEL IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n)                   &
!C!!$OMP          SHARED(cor_mo_iter, nsurft, surft_pts, surft_index,             &
!C!!$OMP                 u_s_iter_surft, u_s_std_surft, v_s_std)
!C!IF ( cor_mo_iter >= use_correct_ustar ) THEN
!C!  !       Use correct "standard" ustar
!C!!$OMP DO SCHEDULE(STATIC)
!C!  DO n = 1,nsurft
!C!    DO k = 1,surft_pts(n)
!C!      l = surft_index(k,n)
!C!      u_s_std_surft(l,n) = v_s_std(l,n)
!C!    END DO
!C!  END DO
!C!!$OMP END DO
!C!ELSE
!C!  !       Use ustar from mid-iteration
!C!!$OMP DO SCHEDULE(STATIC)
!C!  DO n = 1,nsurft
!C!    DO k = 1,surft_pts(n)
!C!      l = surft_index(k,n)
!C!      u_s_std_surft(l,n) = u_s_iter_surft(l,n)
!C!    END DO
!C!  END DO
!C!!$OMP END DO
!C!END IF
!C!!$OMP END PARALLEL
!CABLE_LSM:} 
!-----------------------------------------------------------------------
!  3.5 Recalculate friction velocity for the dust scheme using the
!      bare soil roughness length if using only 1 aggregated tile
!-----------------------------------------------------------------------

IF ((l_dust .OR. l_dust_diag) .AND. l_aggregate) THEN

  n = 1
  ! Calculate z0m and z0h for bare soil
  IF (l_vary_z0m_soil) THEN
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(k, l)            &
!$OMP             SHARED(n, soil, surft_pts, surft_index, z0m_soil_in,        &
!$OMP                   z0h_soil,z0h_z0m, z0m_soil, wind_profile_fac_soil)    &
!$OMP             SCHEDULE(STATIC)
    DO k = 1, surft_pts(n)
      l = surft_index(k,n)
      z0m_soil(l,n) = z0m_soil_in(l)
      z0h_soil(l,n) = z0h_z0m(soil) * z0m_soil(l,n)
      !     Set wind profile factor to 1 as not using orog term in z0m
      wind_profile_fac_soil(l,n) = 1.0
    END DO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(k, l)            &
!$OMP             SHARED(n, soil, surft_pts, surft_index, z0_soil, z0h_soil,  &
!$OMP                   z0h_z0m, z0m_soil, wind_profile_fac_soil)             &
!$OMP             SCHEDULE(STATIC)
    DO k = 1, surft_pts(n)
      l = surft_index(k,n)
      z0m_soil(l,n) = z0_soil
      z0h_soil(l,n) = z0h_z0m(soil) * z0m_soil(l,n)
      !     Set wind profile factor to 1 as not using orog term in z0m
      wind_profile_fac_soil(l,n) = 1.0
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Call fcdch again to calculate dust friction velocity on bare soil.
  ! The canopy drag scheme is not available on the aggregated tile
  ! and is disabled.
  l_vegdrag_active_here = .FALSE.
  CALL fcdch (                                                                &
   cor_mo_iter,land_pts,surft_pts(n),                                         &
   surft_index(:,n),land_index,                                               &
   db_surft(:,n),vshr_land,                                                   &
   z0m_soil(:,n),z0h_soil(:,n),zdt_dummy,zh,                                  &
   z1_uv,z1_uv_top,z1_tq,z1_tq_top,wind_profile_fac_soil(:,n),                &
   ddmfx,ip_ss_solid,charnock,                                                &
   charnock_w,                                                                &
   l_vegdrag_active_here,array_zero,array_zero,                               &
  !  Following tiled outputs (except v_s_std_soil and u_s_iter_soil)
  !  are dummy variables not needed from this call
     cd_surft_soil(:,n),ch_surft_soil(:,n),cd_std_soil(:,n),                  &
     v_s_surft_soil(:,n),v_s_std_soil(:,n),recip_l_mo_surft_soil(:,n),        &
     u_s_iter_soil(:,n)                                                       &
     )

!$OMP PARALLEL IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n)                   &
!$OMP          SHARED(cor_mo_iter, nsurft, surft_index, surft_pts,             &
!$OMP                 u_s_iter_soil, u_s_std_surft, v_s_std_soil)
  IF ( cor_mo_iter >= use_correct_ustar ) THEN
    !       Use correct "standard" ustar
!$OMP DO SCHEDULE(STATIC)
    DO n = 1,nsurft
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        u_s_std_surft(l,n) = v_s_std_soil(l,n)
      END DO
    END DO
!$OMP END DO
  ELSE
    !       Use ustar from mid-iteration
!$OMP DO SCHEDULE(STATIC)
    DO n = 1,nsurft
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        u_s_std_surft(l,n) = u_s_iter_soil(l,n)
      END DO
    END DO
!$OMP END DO
  END IF
!$OMP END PARALLEL

END IF

!-----------------------------------------------------------------------
!  3.6 Recalculate cd, ch etc. using z0h=z0h_classic. The parameters
!       calculated using this additional roughness length are for
!       CLASSIC aerosol deposition only.
!-----------------------------------------------------------------------

IF (l_aero_classic) THEN
  ! The canopy drag scheme is not supported for this aerosol scheme and
  ! is turned off.
  l_vegdrag_active_here = .FALSE.
  ! Land tiles
  DO n = 1,nsurft
    CALL fcdch (                                                              &
    !  Input variables identical to main call except using different z0h
         cor_mo_iter,land_pts,surft_pts(n),                                   &
         surft_index(:,n),land_index,                                         &
         db_surft(:,n),vshr_land,                                             &
         z0m_eff_surft(:,n),z0h_surft_classic(:,n),zdt_dummy,zh,              &
         z1_uv,z1_uv_top,z1_tq,z1_tq_top,                                     &
         wind_profile_factor(:,n),                                            &
         ddmfx,ip_ss_solid,charnock,                                          &
         charnock_w,                                                          &
         l_vegdrag_active_here,array_zero,array_zero,                         &
    !  Following tiled outputs (except cd_std_classic and ch_surft_classic)
    !  are dummy variables not needed from this call
         cd_surft_classic(:,n),ch_surft_classic(:,n),                         &
         cd_std_classic(:,n),v_s_surft_classic(:,n),                          &
         v_s_std_classic(:,n),recip_l_mo_surft_classic(:,n),                  &
         u_s_iter_classic(:,n)                                                &
         )
  END DO
END IF

!-----------------------------------------------------------------------
!  4.1 Recalculate RESFT using "true" CH and EPDT for land tiles
!-----------------------------------------------------------------------

DO n = 1,nsurft
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l) &
!$OMP SHARED(ch_surft, land_index, qstar_surft, qw_1, dq, epdt,n,        &
!$OMP rhostar, surft_pts, surft_index, timestep, t_i_length,             &
!$OMP vshr_land) SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    dq(l) = qw_1(i,j) - qstar_surft(l,n)
    epdt(l) = - rhostar(i,j) * ch_surft(l,n) * vshr_land(i,j)                 &
      *dq(l) * timestep
  END DO
!$OMP END PARALLEL DO

  ! We should only attempt to access sf_diag%resfs_stom(:,n) if it has
  ! been fully allocated.
  IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
    n_diag = n
  ELSE
    n_diag = 1
  END IF

!CABLE_LSM:{ 
!C!  CALL sf_resist (                                                            &
!C!   land_pts,surft_pts(n),land_index,surft_index(:,n),                         &
!C!   canopy(:,n),catch(:,n),ch_surft(:,n),dq,epdt,flake(:,n),                   &
!C!   gc(:,n),gc_stom_surft(:,n),snowdep_surft(:,n),vshr_land,fraca(:,n),        &
!C!   resfs(:,n),resft(:,n),                                                     &
!C!   gc_irr_surft(:,n),resfs_irr_surft(:,n),sf_diag%resfs_stom(:,n_diag),       &
!C!   sf_diag%l_et_stom,sf_diag%l_et_stom_surft)   

END DO

!C!IF ( .NOT. l_aggregate .AND. can_model == 4) THEN
!C!  DO n = 1,npft
!C!    IF ( cansnowtile(n) ) THEN
!C!!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)       &
!C!!$OMP             SHARED(ch_surft, fraca, gc, land_index, l_irrig_dmd, n,      &
!C!!$OMP                    snow_surft, resfs, resfs_irr_surft, resft, surft_pts, &
!C!!$OMP                    surft_index, t_i_length, vshr_land) SCHEDULE(STATIC)
!C!      DO k = 1,surft_pts(n)
!C!        l = surft_index(k,n)
!C!        j=(land_index(l) - 1) / t_i_length + 1
!C!        i = land_index(l) - (j-1) * t_i_length
!C!        IF (snow_surft(l,n) >  0.0) THEN
!C!          fraca(l,n) = 0.0
!C!          resfs(l,n) = gc(l,n) /                                              &
!C!                  (gc(l,n) + ch_surft(l,n) * vshr_land(i,j))
!C!          resft(l,n) = resfs(l,n)
!C!          IF (l_irrig_dmd) THEN
!C!            resfs_irr_surft(l,n) = resfs(l,n)
!C!          END IF
!C!        END IF
!C!      END DO
!C!!$OMP END PARALLEL DO
!C!    END IF
!C!  END DO
!C!END IF
!CABLE_LSM:} 

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
    rhokm_1_surft(l,n) = rhostar(i,j) * cd_surft(l,n) * vshr_land(i,j)
    !                                                         ! P243.124
    rhokm_land(i,j) = rhokm_land(i,j) +                                       &
           tile_frac(l,n) * rhokm_1_surft(l,n)
    rhokh_1(l,n) = rhostar(i,j) * ch_surft(l,n) * vshr_land(i,j)
    !                                                         ! P243.125
  END DO
END DO

!-----------------------------------------------------------------------
!  Calculate local and gridbox-average surface fluxes of heat and
!  moisture.
!-----------------------------------------------------------------------

! Adjust ASHTF for sens. heat flux to ground beneath coniferous canopy
rhokh_can(:,:) = 0.0

IF ( .NOT. l_aggregate .AND. can_model == 4) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n, k, l, j, i) SHARED(cansnowtile,     &
!$OMP             cd_surft, land_index, npft, rhokh_can, rhostar, surft_pts,   &
!$OMP             surft_index, t_i_length, vshr_land, cp) SCHEDULE(STATIC)
  DO n = 1,npft
    IF ( cansnowtile(n) ) THEN
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j-1) * t_i_length
        rhokh_can(l,n) = rhostar(i,j) * cp /                                  &
                     (43.0 / (SQRT(cd_surft(l,n)) * vshr_land(i,j)))
      END DO
    END IF
  END DO
!$OMP END PARALLEL DO
END IF

! Initialise scaling_urban to 1.0 so that it only affects urban tiles when
! MORUSES used with no aggregation.
scaling_urban(:,:) = 1.0
IF ( .NOT. l_aggregate .AND. l_moruses_storage ) THEN
  n = urban_canyon
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(land_pts,   &
!$OMP             n, scaling_urban, tile_frac, urban_roof) SCHEDULE(STATIC)
  DO l = 1, land_pts
    IF ( tile_frac(l,n) > 0.0 ) THEN
      scaling_urban(l,n) =                                                    &
         ( tile_frac(l,n) + tile_frac(l,urban_roof) ) /                       &
         tile_frac(l,n)
    END IF
  END DO
!$OMP END PARALLEL DO
END IF

! Calculate average layer temperature and conductivity for lakes.
! This is a fudge - a layer with average properties won't
! really behave like a stack of layers with different properties.
!
IF (     l_flake_model                                                        &
    .AND. ( .NOT. l_aggregate)) THEN

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


    !Calculate GBM values for soilt tiled variables
  IF (nsoilt == 1) THEN
    !Just 1 soil tile
    m = 1
    hcons_mean_soil(:) = hcons_soilt(:,m)
    tsoil_mean_soil(:) = t_soil_soilt(:,m)
  ELSE
    !Surface tiles map directly on to soil tiles
    hcons_mean_soil(:) = 0.0
    tsoil_mean_soil(:) = 0.0
    DO m = 1,nsoilt
      hcons_mean_soil(:) = hcons_mean_soil(:)                                 &
                          + (tile_frac(:,m) * hcons_soilt(:,m))
      tsoil_mean_soil(:) = tsoil_mean_soil(:)                                 &
                          + (tile_frac(:,m) * t_soil_soilt(:,m))
    END DO
  END IF


  !==============================================================================
  ! *END NOTICE REGARDING SOIL TILING**
  !==============================================================================

!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE)                              &
!$OMP PRIVATE(l, errcode, cmessage)                                           &
!$OMP SHARED(land_pts, hcon_lake, dzsoil, hcons_mean_soil, lake_h_snow_gb,    &
!$OMP        lake_h_ice_gb, snow_hcon, lake_depth_gb, nusselt_gb,             &
!$OMP        g_dt_gb, ts1_lake_gb, lake_t_snow_gb, lake_t_mxl_gb,             &
!$OMP        lake_t_ice_gb, lake_h_mxl_gb, lake_t_mean_gb, tsoil_mean_soil)   &
!$OMP SCHEDULE(STATIC)
  DO l = 1, land_pts

    !
    ! first calculate HCON_LAKE
    !---------------------------
    hcon_lake(l) = 0.0

    IF (dzsoil <= 0.0) THEN
      !
      ! catch-all for sillies - just use land value
      !
      errcode = -1
      WRITE(cmessage, '(A,F16.4)') 'Negative value of dzsoil = ', dzsoil
      CALL ereport(routinename, errcode, cmessage)
      hcon_lake(l) = hcons_mean_soil(l)

    ELSE IF (     (dzsoil  > 0.0)                                             &
             .AND. (dzsoil <= lake_h_snow_gb(l)) ) THEN

      hcon_lake(l) = snow_hcon

    ELSE IF (     (dzsoil  >  lake_h_snow_gb(l))                              &
             .AND. (dzsoil <= (lake_h_snow_gb(l)                              &
                             +lake_h_ice_gb( l)))) THEN

      hcon_lake(l) = ( snow_hcon * lake_h_snow_gb(l)                          &
                      +hcice * (dzsoil - lake_h_snow_gb(l)))                  &
                    / dzsoil

    ELSE IF (     (dzsoil  > (lake_h_snow_gb(l)                               &
                             +lake_h_ice_gb( l)))                             &
             .AND. (dzsoil <= (lake_h_snow_gb(l)                              &
                             +lake_h_ice_gb( l)                               &
                             +lake_depth_gb( l)))) THEN

      nusselt_gb(l) = g_dt_gb(l) * (dzsoil - lake_h_snow_gb(l)                &
                                     - lake_h_ice_gb(l) )                     &
                   / ( 2.0 * hcwat )
      nusselt_gb(l) = MAX( nusselt_gb(l), 1.0 )
      hcon_lake(l) = ( snow_hcon * lake_h_snow_gb(l)                          &
                      +hcice * lake_h_ice_gb(l)                               &
                      +hcwat * (dzsoil - lake_h_snow_gb(l)                    &
                                       - lake_h_ice_gb( l)))                  &
                             * nusselt_gb(l)                                  &
                    / dzsoil

    ELSE

      errcode = -1
      WRITE(cmessage, '(A,F16.4)') 'Unusual value of dzsoil found when '//    &
                                   'computing hcon_lake. Found depth '  //    &
                                   'of first soil layer comparable to ' //    &
                                   'lake depth: dzsoil= ', dzsoil
      CALL ereport(routinename, errcode, cmessage)
      nusselt_gb(l) = g_dt_gb(l) * lake_depth_gb(l)                           &
                   / ( 2.0 * hcwat )
      nusselt_gb(l) = MAX( nusselt_gb(l), 1.0 )
      hcon_lake(l) = ( snow_hcon * lake_h_snow_gb(l)                          &
                      +hcice * lake_h_ice_gb(l)                               &
                      +hcwat * lake_depth_gb(l)                               &
                             * nusselt_gb(l)                                  &
                      +hcons_mean_soil(l) * (dzsoil - lake_h_snow_gb(l)       &
                                             - lake_h_ice_gb( l)              &
                                             - lake_depth_gb( l)))            &
                    / dzsoil

    END IF
    !
    ! second calculate TS1_LAKE
    !----------------------------
    ts1_lake_gb(l) = 0.0

    IF (dzsoil <= 0.0) THEN
      !
      ! catch-all for sillies - just use land value
      !

      ts1_lake_gb(l) = tsoil_mean_soil(l)

    ELSE IF (     ((dzsoil / 2.0)  > 0.0)                                     &
             .AND. ((dzsoil / 2.0) <= lake_h_snow_gb(l)) ) THEN

      ts1_lake_gb(l) = lake_t_snow_gb(l) +                                    &
                    (lake_t_ice_gb(l) - lake_t_snow_gb(l))                    &
                   *(dzsoil / 2.0)                                            &
                   / lake_h_snow_gb(l)

    ELSE IF (     ((dzsoil / 2.0)  >  lake_h_snow_gb(l))                      &
             .AND. ((dzsoil / 2.0) <= (lake_h_snow_gb(l)                      &
                                   +lake_h_ice_gb( l)))) THEN

      ts1_lake_gb(l) = lake_t_ice_gb(l) +                                     &
                    (tm - lake_t_ice_gb(l))                                   &
                   *((dzsoil / 2.0) - lake_h_snow_gb(l))                      &
                   / lake_h_ice_gb(l)

    ELSE IF (     ((dzsoil / 2.0)  > (lake_h_snow_gb(l)                       &
                                   +lake_h_ice_gb( l)))                       &
             .AND. ((dzsoil / 2.0) <= (lake_h_snow_gb(l)                      &
                                   +lake_h_ice_gb( l)                         &
                                   +lake_h_mxl_gb( l)))) THEN

      ts1_lake_gb(l) = lake_t_mxl_gb(l)

    ELSE IF (     ((dzsoil / 2.0)  > (lake_h_snow_gb(l)                       &
                                   +lake_h_snow_gb(l)                         &
                                   +lake_h_mxl_gb( l)))                       &
             .AND. ((dzsoil / 2.0) <= (lake_h_snow_gb(l)                      &
                                   +lake_h_ice_gb( l)                         &
                                   +lake_depth_gb( l)))) THEN

      ts1_lake_gb(l) = lake_t_mean_gb(l)

    ELSE
      errcode = -1
      WRITE(cmessage, '(A,F16.4)') 'Unusual value of dzsoil found when '//    &
                                   'computing ts1_lake_gb. Found depth '//    &
                                   'of first soil layer comparable to ' //    &
                                   'lake depth: dzsoil= ', dzsoil
      CALL ereport(routinename, errcode, cmessage)
      ts1_lake_gb(l) = tsoil_mean_soil(l)

    END IF
  END DO
!$OMP END PARALLEL DO
END IF

!CABLE_LSM:{ move from here to start?
if(first_call) emis_surft = 1.0
!CABLE_LSM}:

lw_down_elevcorr_surft(:,:) = 0.0

IF (l_elev_lw_down) THEN
  ! for tiles at different elevations, adjust downwelling longwave
  ! according to the anomaly in surface temperature that has been
  ! calculated with LW ~ T^4
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length

    IF (lw_down(i,j) > 0.0) THEN

      lw_down_surftsum = 0.0
      lw_down_surftabs = 0.0

      DO n = 1,nsurft
        t_rad = 0.0
        IF (t_elev(l,n) > 0.0 ) THEN
          t_rad = (lw_down(i,j) / sbcon)**(1.0 / 4.0)
          t_rad = t_rad + t_elev(l,n) - tl_1(i,j)
          lw_down_elevcorr_surft(l,n) =  sbcon * (t_rad**4) - lw_down(i,j)
        END IF
        lw_down_surftsum = lw_down_surftsum +                                 &
                           lw_down_elevcorr_surft(l,n) * tile_frac(l,n)
        lw_down_surftabs = lw_down_surftabs +                                 &
                           ABS(lw_down_elevcorr_surft(l,n)) * tile_frac(l,n)
      END DO

      IF (lw_down_surftabs > EPSILON(0.0)) THEN
        ! correct each adjustment to preserve the gridbox mean.
        ! size of correction in proportion to the size of adjustment
        ! so that unadjusted tiles remain unaffected
        DO n = 1,nsurft
          lw_down_elevcorr_surft(l,n) = lw_down_elevcorr_surft(l,n)           &
                                        - lw_down_surftsum                    &
                                          * ABS(lw_down_elevcorr_surft(l,n))  &
                                            / lw_down_surftabs
        END DO
      END IF

!CABLE_LSM:{ 
!C!      DO n = 1,nsurft
!C!        IF (l_skyview) THEN
!C!          radnet_surft(l,n) = radnet_surft(l,n)                               &
!C!                              + sky(i,j) * emis_surft(l,n)                    &
!C!                                * lw_down_elevcorr_surft(l,n)
!C!        ELSE
!C!          radnet_surft(l,n) = radnet_surft(l,n)                               &
!C!                              + emis_surft(l,n) * lw_down_elevcorr_surft(l,n)
!C!        END IF
!C!      END DO
!CABLE_LSM:} 

    END IF

  END DO
END IF

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

DO n = 1,nsurft
  lh0 = lc

  !Set the current soil tile (see notice above)
  IF (nsoilt == 1) THEN
    !There is only 1 soil tile
    m = 1
  ELSE ! nsoilt == nsurft
    !Soil tiles map directly on to surface tiles
    m = n
  END IF !nsoilt

!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(i, j, l, t_rad)      &
!$OMP SHARED(land_pts, land_index, t_i_length, tsurf, t_soil_soilt, t_elev,   &
!$OMP        tl_1,                                                            &
!$OMP        dzsurf, dzsoil, hcons_surf, hcons_soilt, canhc_surf, canhc_surft,&
!$OMP        nsmax, nsnow, tsnow, ds, hcons_snow, cansnowtile,                &
!$OMP        snowdepth, snow_hcon,                                            &
!$OMP        l_snow_nocan_hc, l_flake_model, l_aggregate, lake,               &
!$OMP        n, m, ts1_lake_gb, hcon_lake, l_elev_land_ice, l_lice_point,     &
!$OMP        tsurf_elev_surft, dzsoil_elev, l_lice_surft, hcondeep)           &
!$OMP SCHEDULE(STATIC)
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    IF (l_elev_land_ice .AND. l_lice_point(l)) THEN
      tsurf(l) = tsurf_elev_surft(l,n)
      dzsurf(l) = dzsoil_elev
      canhc_surf(l) = 0.0
      IF (l_lice_surft(n)) THEN
        hcons_surf(l) = snow_hcon
      ELSE
        hcons_surf(l) = hcondeep
      END IF
    ELSE
      tsurf(l) = t_soil_soilt(l,m) + t_elev(l,n) - tl_1(i,j)
      dzsurf(l) = dzsoil
      hcons_surf(l) = hcons_soilt(l,m)
      canhc_surf(l) = canhc_surft(l,n)
    END IF
    IF ((nsmax > 0) .AND. (nsnow(l,n) > 0)) THEN
      tsurf(l) = tsnow(l,n,1)
      ! change the effective surface layer thickness for snow
      dzsurf(l) = ds(l,n,1)
      hcons_surf(l) = hcons_snow(l,n)
      IF ( ( .NOT. cansnowtile(n)) .AND. l_snow_nocan_hc .AND.                &
           (nsmax > 0) .AND. (nsnow(l,n) > 0) ) canhc_surf(l) = 0.0
    END IF

    IF (     (l_flake_model   )                                               &
        .AND. ( .NOT. l_aggregate)                                            &
        .AND. (n == lake       )                                              &
        ! single-layer snow scheme OR negligible snow
        .AND. ( (nsmax == 0) .OR. (snowdepth(l,n) < EPSILON(1.0) ) )          &
       ) THEN
      tsurf(l) = ts1_lake_gb(l)
      hcons_surf(l) = hcon_lake(l)
      dzsurf(l) = dzsoil
    END IF

  END DO
!$OMP END PARALLEL DO

  !==============================================================================
  ! *END NOTICE REGARDING SOIL TILING**
  !==============================================================================

  sea_point = 0.0
!CABLE_LSM:{
!C!  CALL sf_flux (                                                              &
!C!   land_pts,surft_pts(n),flandg,                                              &
!C!   land_index,surft_index(:,n),                                               &
!C!   nsnow(:,n),n,canhc_surf(:),dzsurf,hcons_surf,                              &
!C!   qs1_elev(:,n),qstar_surft(:,n),q_elev(:,n),                                &
!C!   radnet_surft(:,n),resft(:,n),rhokh_1(:,n),l_soil_point,                    &
!C!   snowdepth(:,n),tile_frac(:,n),timestep,t_elev(:,n),tsurf,                  &
!C!   tstar_surft(:,n),vfrac_surft(:,n),rhokh_can(:,n),z0h_surft(:,n),           &
!C!   z0m_eff_surft(:,n),zdt_surft(:,n),z1_tq,lh0,emis_surft(:,n),emis_soil,     &
!C!   1.0,anthrop_heat(:,n),scaling_urban(:,n),l_vegdrag_surft(n),               &
!C!   fqw_1,ftl_1,                                                               &
!C!   alpha1(:,n),ashtf_prime_surft(:,n),fqw_surft(:,n),                         &
!C!   epot_surft(:,n),ftl_surft(:,n),                                            &
!C!   rhokm_1_surft(:,n),vshr_land,tau_surft(:,n),                               &
!C!   dtstar_surft(:,n),sea_point, sf_diag                                       &
!C! )
!CABLE_LSM:}
END DO

! Set surface stress on tiles diagnostic
IF (sf_diag%l_tau_surft) THEN
  DO n = 1,nsurft
    DO l = 1,land_pts
      sf_diag%tau_surft(l,n) = tau_surft(l,n)
    END DO
  END DO
END IF

! tau_surft no longer required
DEALLOCATE(tau_surft)

!CABLE_LSM:{ Elemnts of sf_flux for SEB calc, carried from ACCESS1.3 - review

!jhan:added from if(nsoilt) block as was inside flake section and we use
!hcons_mean_soil here
hcons_mean_soil(:) = hcons_soilt(:,1)

DO n = 1,nsurft
   DO K=1,surft_PTS(N)
     L = surft_INDEX(K,N)
     J=(LAND_INDEX(L)-1)/t_i_length+ 1
     I = LAND_INDEX(L) - (J-1)*t_i_length
     fD_T = TSTAR_surft(L,N) - TL_1(I,J)
     IF (fD_T  >   0.05 .OR. fD_T  <   -0.05) THEN
       ALPHA1(L,N) = (QSTAR_surft(L,N) - QS1(I,J)) / fD_T
     ELSEIF (TL_1(I,J)  >   TM) THEN
       ALPHA1(L,N) = fepsilon*LC*QS1(I,J)*                              &
         (1. + fC_VIRTUAL*QS1(I,J)) / ( R*TL_1(I,J)*TL_1(I,J))
     ELSE
       ALPHA1(L,N) = FePsilon*LS*QS1(I,J)*                              &
         (1. + fC_VIRTUAL*QS1(I,J)) / ( R*TL_1(I,J)*TL_1(I,J))
     ENDIF
    ! This needs to be look at:
     ASHTF_prime_surft(L,N) = 2.0 * hcons_mean_soil(L) / DZSOIL
     IF (SNOW_surft(L,N) >  0.0 .AND. SMVCST_gb(L) /= 0.) THEN
       fDS_RATIO = 2.0 * SNOW_surft(L,N) / (RHO_SNOW_const * DZSOIL)
       IF (fDS_RATIO <= 1.0) THEN
         ASHTF_prime_surft(L,N) =  ASHTF_prime_surft(L,N) /                         &
                     (1. + fDS_RATIO*(hcons_mean_soil(L)/SNOW_HCON - 1.))
       ELSE
         ASHTF_prime_surft(L,N) =  ASHTF_prime_surft(L,N)*SNOW_HCON / hcons_mean_soil(L)
       ENDIF
     ENDIF
     !jhan: LH may-be lh0 here 
     fLH = LC
     IF (SNOW_surft(L,N)  >   0.) fLH = LS
     RHOKPM(L,N) = RHOKH_1(L,N) / ( ASHTF_prime_surft(L,N)  +                &
                   RHOKH_1(L,N)*(fLH*ALPHA1(L,N)*RESFT(L,N) + CP))
     !jhan:ACCESS1.3 left this in but never used. only in unused calc of epot
     !RHOKPM_POT(L,N)=RHOKH_1(L,N) / ( ASHTF_prime_TILE(L,N)  +              &
     !                 RHOKH_1(L,N)*(cable% tmp% LH*ALPHA1(L,N) + CP) )
     FTL_1(I,J) = FTL_1(I,J)+FLAND(L)*tile_FRAC(L,N)*FTL_surft(L,N)
     FQW_1(I,J) = FQW_1(I,J)+FLAND(L)*tile_FRAC(L,N)*FQW_surft(L,N)
  
   ENDDO
!CABLE_LSM:End 

END DO



!-----------------------------------------------------------------------
!  4.4   Calculate the standard deviations of layer 1 turbulent
!        fluctuations of temperature and humidity using approximate
!        formulae from first order closure.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  CALL stdev1 (                                                               &
   land_pts,surft_pts(n),land_index,surft_index(:,n),flandg,                  &
   bq_1,bt_1,fqw_surft(:,n),ftl_surft(:,n),rhokm_1_surft(:,n),                &
   rhostar,vshr_land,z0m_surft(:,n),z1_tq,tile_frac(:,n),                     &
   q1_sd,t1_sd                                                                &
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
IF ( l_fix_wind_snow ) THEN
  DO n = 1,nsurft
    IF ( canSnowTile(n) .AND. unload_rate_u(n) /= 0.0 ) l_cdr10m_snow = .TRUE.
  END DO
END IF

IF (sf_diag%suv10m_n) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j) &
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
IF (sf_diag%su10 .OR. sf_diag%sv10 .OR. sf_diag%sq1p5 .OR.                    &
    sf_diag%st1p5 .OR. sf_diag%suv10m_n .OR.                                  &
    l_cdr10m_snow .OR.                                                        &
    (IScrnTDiag == IP_ScrnDecpl2) .OR.                                        &
    (IScrnTDiag == IP_ScrnDecpl3) ) THEN
  DO n = 1,nsurft
    CALL sfl_int (                                                            &
     land_pts,surft_pts(n),l_cdr10m_snow,surft_index(:,n),land_index,flandg,  &
     vshr_land,cd_std(:,n),cd_surft(:,n),ch_surft(:,n),                       &
     tile_frac(:,n),                                                          &
     z0m_eff_surft(:,n),z0m_surft(:,n),z0h_surft(:,n),                        &
     recip_l_mo_surft(:,n),                                                   &
     v_s_surft(:,n),v_s_std(:,n),                                             &
     z1_uv,z1_tq,db_surft(:,n),                                               &
     sf_diag,                                                                 &
     cdr10m,sf_diag%cdr10m_n,sf_diag%cd10m_n,chr1p5m(:,n)                     &
     )
  END DO

END IF

!-----------------------------------------------------------------------
! Sea and sea-ice surface calculations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Calculate height of lowest model level relative to sea.
!-----------------------------------------------------------------------

IF (i_modiscopt == on) THEN
  ALLOCATE(z1_tq_top_sea(t_i_length,t_j_length))
  ALLOCATE(z1_tq_top_ctile(t_i_length,t_j_length))
ELSE
  ALLOCATE(z1_tq_top_sea(1,1))
  ALLOCATE(z1_tq_top_ctile(1,1))
END IF

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,z1_tq_sea,z1_tq,z1_tq_ctile,flandg,z_land,l_ctile,         &
!$OMP        l_fix_ctile_orog,z1_tq_top_sea,z1_tq_top,z1_tq_top_ctile,        &
!$OMP        i_modiscopt)

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    z1_tq_sea(i,j)   = z1_tq(i,j)
    z1_tq_ctile(i,j) = z1_tq(i,j)
  END DO
END DO
!$OMP END DO NOWAIT

IF (l_ctile) THEN !ij
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF ((flandg(i,j) > 0.0) .AND. (flandg(i,j) < 1.0)) THEN
        z1_tq_sea(i,j) = z1_tq(i,j) + z_land(i,j)
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (l_ctile .AND. .NOT. l_fix_ctile_orog) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF ((flandg(i,j) > 0.0) .AND. (flandg(i,j) < 1.0)) THEN
        z1_tq_ctile(i,j) = z1_tq(i,j) + z_land(i,j)
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
END IF


IF (i_modiscopt == on) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      z1_tq_top_sea(i,j) = z1_tq_top(i,j)
      z1_tq_top_ctile(i,j) = z1_tq_top(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (l_ctile) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        IF ((flandg(i,j) > 0.0) .AND. (flandg(i,j) < 1.0)) THEN
          z1_tq_top_sea(i,j) = z1_tq_top(i,j) + z_land(i,j)
        END IF
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_ctile .AND. .NOT. l_fix_ctile_orog) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        IF ((flandg(i,j) > 0.0) .AND. (flandg(i,j) < 1.0)) THEN
          z1_tq_top_ctile(i,j) = z1_tq_top(i,j) + z_land(i,j)
        END IF
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
END IF

!$OMP END PARALLEL

IF (nice_use >  1) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,n) &
!$OMP SHARED(tdims,ice_fract,tstar_sice,nice_use,ice_fract_cat,tstar_sice_cat)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      ice_fract(i,j) = 0.0
      tstar_sice(i,j) = 0.0
      DO n = 1, nice_use
        IF (ice_fract_cat(i,j,n) >  0.0) THEN
          ice_fract(i,j)  = ice_fract(i,j) + ice_fract_cat(i,j,n)
          tstar_sice(i,j) = tstar_sice(i,j) +                                 &
                            tstar_sice_cat(i,j,n) * ice_fract_cat(i,j,n)

        END IF
      END DO
      IF (ice_fract(i,j) >  0.0)                                              &
         tstar_sice(i,j) = tstar_sice(i,j) / ice_fract(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO
  sice_pts_use = 0
  sice_index_use(:) = 0
  sice_frac_use(:) = 0.0
  DO l = 1, ssi_pts
    first_counter = 0
    DO n = 1, nice_use
      IF (sice_frac_ncat(l,n) >  0.0) THEN
        IF (first_counter == 0) THEN
          sice_pts_use = sice_pts_use + 1
          sice_index_use(sice_pts_use) = l
          first_counter = 1
        END IF
        sice_frac_use(l) = sice_frac_use(l) + sice_frac_ncat(l,n)
      END IF
    END DO
  END DO
ELSE
  ice_fract(:,:) = ice_fract_cat(:,:,1)
  tstar_sice(:,:) = tstar_sice_cat(:,:,1)
  sice_pts_use = sice_pts_ncat(1)
  sice_index_use(:) = sice_index_ncat(:,1)
  sice_frac_use(:) = sice_frac_ncat(:,1)
END IF

! DEPENDS ON: qsat_mix
qsat_tile = 0

IF (l_new_qsat_jules) THEN
  IF (lq_mix_bl) THEN
    CALL qsat_mix_new(qstar_sea,tstar_sea,pstar,t_i_length,t_j_length)
    CALL qsat_mix_new(qstar_ice,tstar_sice,pstar,t_i_length,t_j_length)
  ELSE
    CALL qsat_new(qstar_sea,tstar_sea,pstar,t_i_length,t_j_length)
    CALL qsat_new(qstar_ice,tstar_sice,pstar,t_i_length,t_j_length)
  END IF
ELSE
  CALL qsat_mix(qstar_sea,tstar_sea,pstar,t_i_length * t_j_length,lq_mix_bl)
  CALL qsat_mix(qstar_ice,tstar_sice,pstar,t_i_length * t_j_length,lq_mix_bl)
END IF

!
! Also calculate surface specific saturated humidity on ice by category, as
! it's needed in the call to sf_flux.
IF (l_new_qsat_jules) THEN
  IF (lq_mix_bl) THEN
    DO n = 1, nice_use
      CALL qsat_mix_new(qstar_ice_cat(:,:,n),tstar_sice_cat(:,:,n),pstar,     &
                    t_i_length,t_j_length)
    END DO
  ELSE
    DO n = 1, nice_use
      CALL qsat_new(qstar_ice_cat(:,:,n),tstar_sice_cat(:,:,n),pstar,         &
                    t_i_length,t_j_length)
    END DO
  END IF
ELSE
  DO n = 1, nice_use
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(qstar_ice_cat(:,:,n),tstar_sice_cat(:,:,n),pstar,           &
                  t_i_length * t_j_length,lq_mix_bl)
  END DO
END IF

! Sea, sea-ice leads, sea-ice and marginal ice zone
hcons_sea(:)    = kappa_seasurf
canhc_sea(:)    = hcap_sea


!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,z0h_sea,z0m_miz,z0miz,z0h_z0m_miz,z0m_ice,z0h_ice,         &
!$OMP        z0h_z0m_sice, z0hsea,rib_sea,rib_ice,db_sea,db_ice,rhokm_ssi,    &
!$OMP        rhokm_ssi_nohalo,cd_ssi,ch_ssi,ftl_ice,fqw_ice,h_sea,e_sea,      &
!$OMP        z0h_miz,z0sice)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    z0h_sea(i,j)          = z0hsea
    !   Potentially remove these settings if using form drag
    z0m_miz(i,j)          = z0miz
    z0h_miz(i,j)          = z0h_z0m_miz * z0m_miz(i,j)
    z0m_ice(i,j,:)        = z0sice
    z0h_ice(i,j,:)        = z0h_z0m_sice * z0m_ice(i,j,:)
    rib_sea(i,j)          = 0.0
    rib_ice(i,j,:)        = 0.0
    db_sea(i,j)           = 0.0
    db_ice(i,j,:)         = 0.0
    rhokm_ssi(i,j)        = 0.0
    rhokm_ssi_nohalo(i,j) = 0.0
    cd_ssi(i,j)           = 0.0
    ch_ssi(i,j)           = 0.0
    ftl_ice(i,j,:)        = 0.0
    fqw_ice(i,j,:)        = 0.0
    h_sea(i,j)            = 0.0
    e_sea(i,j)            = 0.0
  END DO
END DO
!$OMP END PARALLEL DO


SELECT CASE (iseasurfalg)
  !
CASE (ip_ss_surf_div)

  !       Composite formulation for thermal roughness lengths,
  !       incoporating the smooth aerodynamic limit for low
  !       wind speeds and a value based on surface divergence
  !       theory for higher wind speeds.

  !       The friction velocity is diagnosed in the surface
  !       transfer scheme, using z0m from the previous time-step.
  !       z0[T,q] is also required but depends on u_* in this
  !       scheme. For practical purposes, it is sufficient to
  !       infer it from u_* determined from z0m, but in general
  !       a given value of z0m corresponds to two values of
  !       u_*, so we need to know whether we are on the low or
  !       high wind speed side of the minimum value of z0m.
  !       If we are on the high side, z0[T,q] will be inversely
  !       proportional to z0m, but on the low side it may follow
  !       this relationship, or be aerodynamically smooth. In
  !       the smooth case we iterate u_* from the expression for
  !       the roughness length and then take the maximum of the
  !       smooth and high-wind expressions for z0[T,q]. An
  !       iteration for the low-wind value of u_*, USTR_L is
  !       carried out. This will converge to the correct limit only
  !       on the low-wind side of the minimum, and the standard
  !       criterion that the gradient of a fixed-point iteration
  !       should be less than 1 in modulus gives a more precise
  !       boundary, TOL_USTR_L. For consistency with earlier versions
  !       of the modset, hard-wired values are retained for the
  !       operational value of Charnock's parameter. An additional
  !       check is required, since z0m can be large at low or at
  !       high wind-speeds. This is less precise and a fixed
  !       value of 0.07 is used to test USTR_N, which was determined
  !       by inspection of a graph of z0m against u_*: it is
  !       unlikely that this will need to be changed unless
  !       Charnock's constant is greatly altered.

  IF (charnock == 0.018) THEN
    tol_ustr_l = 0.055
    tol_ustr_n = 0.07
  ELSE
    tol_ustr_l = 0.75 * (1.54e-6 * g / (2.0 * charnock))**0.33333
    tol_ustr_n = 0.07
  END IF

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, jits, ustr_l, ustr_n)          &
!$OMP SHARED( tdims, vshr_ssi, z1_uv, z0msea, tol_ustr_n, tol_ustr_l,        &
!$OMP         charnock, g, z0h_sea)                                          &
!$OMP             SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      !           We need to infer u_* from Z0M.
      IF (vshr_ssi(i,j)  > 0.0) THEN
        !             Compute u_* using neutral stratification.
        ustr_n = vkman * vshr_ssi(i,j) /                                      &
          LOG(1.0 + z1_uv(i,j) / z0msea(i,j) )
        !             Compute u_* using low wind approximation.
        ustr_l = 1.54e-06 /  z0msea(i,j) - 1.0e-05
        !             Since Z0M could be large for low and high u_*, we use
        !             USTR_N as an additional check on the regime.
        IF ( (ustr_n < tol_ustr_n) .AND.                                      &
             (ustr_l < tol_ustr_l) ) THEN
          !               Iterate u_* for low winds.
          DO jits = 1, 5
            ustr_l = 1.54e-06 / (z0msea(i,j) - (charnock / g) * ustr_l**2)    &
              -1.0e-05
          END DO
          !               Take the maximum of the smooth and high-wind values.
          !               A lower limit is imposed on the friction velocity to
          !               allow for the exceptional case of very low winds: the
          !               value of 10^-5 is the same as the limit for the momentum
          !               roughness length.
          z0h_sea(i,j) = MAX( 2.52e-6 / (ustr_l+1.0e-05),                     &
            2.56e-9 / z0msea(i,j) )
        ELSE
          !               Take the high-wind value, but limit it to the molecular
          !               mean free path (we should not hit this limit
          !               in practice).
          z0h_sea(i,j) = MAX( 2.56e-9 / z0msea(i,j), 7.0e-08 )
        END IF
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO

CASE (ip_ss_surf_div_coupled)

  ! tol_ustr_n = 0.07 - written as a constant below for optimization
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, jits, ustr_l, ustr_n, tol_ustr_l)    &
!$OMP SHARED( tdims, vshr_ssi, z1_uv, z0msea, charnock_w, g, z0h_sea)         &
!$OMP             SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      tol_ustr_l = 0.75 * (1.54e-6 * g / (2.0 * charnock_w(i,j)))**0.33333

      !           We need to infer u_* from Z0M.
      IF (vshr_ssi(i,j)  > 0.0) THEN
        !             Compute u_* using neutral stratification.
        ustr_n = vkman * vshr_ssi(i,j) /                                      &
          LOG(1.0 + z1_uv(i,j) / z0msea(i,j) )
        !             Compute u_* using low wind approximation.
        ustr_l = 1.54e-06 /  z0msea(i,j) - 1.0e-05
        !             Since Z0M could be large for low and high u_*, we use
        !             USTR_N as an additional check on the regime.
        IF ( (ustr_n < 0.07) .AND.                                            &
             (ustr_l < tol_ustr_l) ) THEN
          !               Iterate u_* for low winds.
          DO jits = 1, 5
            ustr_l = 1.54e-06 / (z0msea(i,j) - (charnock_w(i,j) / g) * ustr_l**2) &
              -1.0e-05
          END DO
          !               Take the maximum of the smooth and high-wind values.
          !               A lower limit is imposed on the friction velocity to
          !               allow for the exceptional case of very low winds: the
          !               value of 10^-5 is the same as the limit for the momentum
          !               roughness length.
          z0h_sea(i,j) = MAX( 2.52e-6 / (ustr_l+1.0e-05),                     &
            2.56e-9 / z0msea(i,j) )
        ELSE
          !               Take the high-wind value, but limit it to the molecular
          !               mean free path (we should not hit this limit
          !               in practice).
          z0h_sea(i,j) = MAX( 2.56e-9 / z0msea(i,j), 7.0e-08 )
        END IF
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO

CASE (ip_ss_fixed)

  !       Use a fixed thermal roughness length.
  z0h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) = z0hsea

CASE DEFAULT
  !     Do not alter the roughness lengths at this point.
  !
END SELECT


IF ( l_spec_z0 ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,z0h_scm,z0h_sea)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF (z0h_scm(i,j)  >   0.0) THEN

        ! Set Z0H from SCM namelist
        ! (if specified) for sea points
        z0h_sea(i,j) = z0h_scm(i,j)

      END IF  ! Z0H_SCM
    END DO  ! I
  END DO  ! J
!$OMP END PARALLEL DO
END IF

!-----------------------------------------------------------------------
!  3.2 Calculate bulk Richardson number for the lowest model level.
!-----------------------------------------------------------------------
!
! Sea, sea-ice and sea-ice leads
CALL sf_rib (                                                                 &
 ssi_pts,sea_pts,ssi_index,sea_index,                                         &
 bq_1,bt_1,qstar_sea,qw_1,array_one,tl_1,                                     &
 tstar_sea,vshr_ssi,z0h_sea,z0msea,array_zero,                                &
 z1_tq_sea,z1_uv,l_vegdrag_ssi,                                               &
 rib_sea,db_sea                                                               &
 )
!
IF (nice_use >  1) THEN
  DO n = 1, nice_use
    CALL sf_rib (                                                             &
     ssi_pts,sice_pts_ncat(n),ssi_index,sice_index_ncat(:,n),                 &
     bq_1,bt_1,qstar_ice_cat(:,:,n),qw_1,array_one,tl_1,                      &
     tstar_sice_cat(:,:,n),vshr_ssi,z0h_ice(:,:,n),z0m_ice(:,:,n),array_zero, &
     z1_tq_sea,z1_uv,l_vegdrag_ssi,                                           &
     rib_ice(:,:,n),db_ice(:,:,n)                                             &
      )
  END DO
ELSE
  CALL sf_rib (                                                               &
    ssi_pts,sice_pts_use,ssi_index,sice_index_use,                            &
    bq_1,bt_1,qstar_ice_cat(:,:,1),qw_1,array_one,tl_1,                       &
    tstar_sice_cat(:,:,1),vshr_ssi,z0h_ice(:,:,1),z0m_ice(:,:,1),array_zero,  &
    z1_tq_sea,z1_uv,l_vegdrag_ssi,                                            &
    rib_ice(:,:,1),db_ice(:,:,1)                                              &
     )

END IF
!
! Calculate gridbox mean db_ice
db_ice_mean(:,:) = 0.0
IF (nice_use >  1) THEN
  DO n = 1, nice_use
    DO k = 1,sice_pts_ncat(n)
      l = sice_index_ncat(k,n)
      j=(ssi_index(l) - 1) / t_i_length + 1
      i = ssi_index(l) - (j-1) * t_i_length
      IF (ice_fract(i,j) >  0.0) THEN
        db_ice_mean(i,j) = db_ice_mean(i,j)                                   &
            + (sice_frac_ncat(l,n) * db_ice(i,j,n)                            &
                         / ice_fract(i,j))
      END IF
    END DO   ! K
  END DO   ! N
ELSE
  db_ice_mean(:,:) = db_ice(:,:,1)
END IF
!
!-----------------------------------------------------------------------
!  3.4 Calculate CD, CH via routine FCDCH.
!-----------------------------------------------------------------------
!
! Initialisations:
cd_ice(:,:,:) = 0.0
ch_ice(:,:,:) = 0.0
cd_std_ice(:,:) = 0.0
v_s_ice(:,:,:) =  0.0
v_s_std_ice(:,:) = 0.0
recip_l_mo_ice(:,:,:) = 0.0
u_s_std_ice(:,:) = 0.0
!
CALL fcdch (                                                                  &
   cor_mo_iter,ssi_pts,sea_pts,                                               &
   sea_index,ssi_index,                                                       &
   db_sea,vshr_ssi,                                                           &
   z0msea,z0h_sea,zdt_dummy,zh,                                               &
   z1_uv,z1_uv_top,z1_tq_ctile,z1_tq_top_ctile,array_one,                     &
   ddmfx,iseasurfalg,charnock,                                                &
   charnock_w,                                                                &
   l_vegdrag_ssi,array_zero,array_zero,                                       &
   cd_sea,ch_sea,cd_std_sea,                                                  &
   v_s_sea,v_s_std_sea,recip_l_mo_sea,                                        &
   u_s_std_sea                                                                &
   )

!
IF (nice_use >  1) THEN
  DO n = 1, nice_use
    CALL fcdch (                                                              &
     cor_mo_iter,ssi_pts,sice_pts_ncat(n),                                    &
     sice_index_ncat(:,n),ssi_index,                                          &
     db_ice(:,:,n),vshr_ssi,                                                  &
     z0m_ice(:,:,n),z0h_ice(:,:,n),zdt_dummy,zh,                              &
     z1_uv,z1_uv_top,z1_tq_ctile,z1_tq_top_ctile,array_one,                   &
     ddmfx,ip_ss_solid,charnock,                                              &
     charnock_w,                                                              &
     l_vegdrag_ssi,array_zero,array_zero,                                     &
     cd_ice(:,:,n),ch_ice(:,:,n),cd_std_ice(:,n),                             &
     v_s_ice(:,:,n),v_s_std_ice(:,n),recip_l_mo_ice(:,:,n),                   &
     u_s_std_ice(:,n)                                                         &
     )
  END DO
ELSE
  CALL fcdch (                                                                &
   cor_mo_iter,ssi_pts,sice_pts_use,                                          &
   sice_index_use,ssi_index,                                                  &
   db_ice(:,:,1),vshr_ssi,                                                    &
   z0m_ice(:,:,1),z0h_ice(:,:,1),zdt_dummy,zh,                                &
   z1_uv,z1_uv_top,z1_tq_ctile,z1_tq_top_ctile,array_one,                     &
   ddmfx,ip_ss_solid,charnock,                                                &
   charnock_w,                                                                &
   l_vegdrag_ssi,array_zero,array_zero,                                       &
   cd_ice(:,:,1),ch_ice(:,:,1),cd_std_ice(:,1),                               &
   v_s_ice(:,:,1),v_s_std_ice(:,1),recip_l_mo_ice(:,:,1),                     &
   u_s_std_ice(:,1)                                                           &
   )
END IF
!
! If using form drag, put the appropriate coefficients into the arrays
! for marginal ice. Scalar transfer is assumed to take place only through
! the interfacial route.
! for form drag.
IF (l_iceformdrag_lupkes) THEN
  CALL ice_formdrag_lupkes (                                                  &
   flandg, ice_fract,                                                         &
   z0m_ice(:,:,1), z0msea, cd_ice(:,:,1), cd_sea,                             &
   z1_tq_sea, z1_tq_top_sea,                                                  &
   cd_miz                                                                     &
 )
ELSE
  !
  CALL fcdch (                                                                &
     cor_mo_iter,ssi_pts,sice_pts_use,                                        &
     sice_index_use,ssi_index,                                                &
     db_ice_mean,vshr_ssi,                                                    &
     z0m_miz,z0h_miz,zdt_dummy,zh,                                            &
     z1_uv,z1_uv_top,z1_tq_ctile,z1_tq_top_ctile,array_one,                   &
     ddmfx,ip_ss_solid,charnock,                                              &
     charnock_w,                                                              &
     l_vegdrag_ssi,array_zero,array_zero,                                     &
     cd_miz,ch_miz,cd_std_miz,                                                &
     v_s_miz,v_s_std_miz,recip_l_mo_miz,                                      &
     u_s_std_miz                                                              &
     )
END IF
!
! z1_tq_top_sea is no longer required.
DEALLOCATE(z1_tq_top_sea)
DEALLOCATE(z1_tq_top_ctile)
!
! Calculate gridbox mean CD and CH for sea ice:
cd_ice_mean(:,:) = 0.0
ch_ice_mean(:,:) = 0.0
IF (nice_use >  1) THEN
  DO n = 1,nice_use
    DO k = 1,sice_pts_ncat(n)
      l = sice_index_ncat(k,n)
      j=(ssi_index(l) - 1) / t_i_length + 1
      i = ssi_index(l) - (j-1) * t_i_length
      IF (ice_fract(i,j) >  0.0) THEN
        cd_ice_mean(i,j) = cd_ice_mean(i,j)                                   &
          + (sice_frac_ncat(l,n) * cd_ice(i,j,n) / ice_fract(i,j))
        ch_ice_mean(i,j) = ch_ice_mean(i,j)                                   &
          + (sice_frac_ncat(l,n) * ch_ice(i,j,n) / ice_fract(i,j))
      END IF
    END DO
  END DO
ELSE
  cd_ice_mean(:,:) = cd_ice(:,:,1)
  ch_ice_mean(:,:) = ch_ice(:,:,1)
END IF
!
! Sea and sea-ice
rhokh_1_sice(:,:) = 0.0
rhokh_1_sice_ncats(:,:,:) = 0.0
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j)                                    &
!$OMP SHARED(tdims, flandg, ice_fract, cd_ssi, cd_miz, cd_sea, ch_ssi, ch_miz,&
!$OMP        ch_sea, cd_ice_mean, ch_ice_mean, nice_use, rhokm_ssi, rhostar,  &
!$OMP        vshr_ssi, rhokm_ssi_nohalo, rhokh_1_sea, rhokh_1_sice, sf_diag,  &
!$OMP        rhokh_1_sice_ncats, ice_fract_cat, ch_ice, l_iceformdrag_lupkes, &
!$OMP        rhokm_1_sea,rhokm_1_sice_ncats,rhokm_1_sice,cd_ice)

IF (l_iceformdrag_lupkes) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF (flandg(i,j) <  1.0 ) THEN
        cd_ssi(i,j) = (1.0 - ice_fract(i,j)) * cd_sea(i,j) +                  &
                      ice_fract(i,j) * ( cd_ice_mean(i,j) +                   &
                                         cd_miz(i,j) )
        !         No contribution to thermal roughness from form drag.
        ch_ssi(i,j) = (1.0 - ice_fract(i,j)) * ch_sea(i,j) +                  &
                      ice_fract(i,j) * ch_ice_mean(i,j)
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF (flandg(i,j) <  1.0 ) THEN
        IF ( ice_fract(i,j) <  0.7 ) THEN
          cd_ssi(i,j) = ( ice_fract(i,j) * cd_miz(i,j) +                      &
                (0.7 - ice_fract(i,j)) * cd_sea(i,j) ) / 0.7  ! P2430.5
          ch_ssi(i,j) = ( ice_fract(i,j) * ch_miz(i,j) +                      &
                (0.7 - ice_fract(i,j)) * ch_sea(i,j) ) / 0.7  ! P2430.4
        ELSE
          cd_ssi(i,j) = ( (1.0 - ice_fract(i,j)) * cd_miz(i,j) +              &
               (ice_fract(i,j) - 0.7) * cd_ice_mean(i,j) ) / 0.3  ! P2430.7
          ch_ssi(i,j) = ( (1.0 - ice_fract(i,j)) * ch_miz(i,j) +              &
               (ice_fract(i,j) - 0.7) * ch_ice_mean(i,j) ) / 0.3  ! P2430.7
        END IF
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

! Sea and sea ice surface drag coefficients
IF (sf_diag%l_cd_ssi) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      sf_diag%cd_ssi(i,j) = cd_ssi(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (sf_diag%l_ch_ssi) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      sf_diag%ch_ssi(i,j) = ch_ssi(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

! Sea and sea-ice
IF (sf_diag%l_tau_surft .OR. sf_diag%l_tau_1) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF ( flandg(i,j) <  1.0 ) THEN
        IF (nice_use > 1) THEN
          !     ! Include effect of leads.
          rhokm_1_sea(i,j) = rhostar(i,j) * cd_ssi(i,j) * vshr_ssi(i,j)
          DO n = 1,nice_use
            rhokm_1_sice_ncats(i,j,n) = rhostar(i,j) * cd_ice(i,j,n) * vshr_ssi(i,j)
            IF (ice_fract(i,j) /= 0.0) THEN
              rhokm_1_sice(i,j) = rhokm_1_sice(i,j) +                         &
                  (rhokm_1_sice_ncats(i,j,n) * ice_fract_cat(i,j,n)           &
                           / ice_fract(i,j))
            END IF
          END DO
        ELSE
          !     ! Original scheme without leads.
          rhokm_1_sice_ncats(i,j,1) = rhostar(i,j) * cd_ssi(i,j) * vshr_ssi(i,j)
          rhokm_1_sice(i,j) = rhokm_1_sice_ncats(i,j,1)
        END IF
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    IF ( flandg(i,j) <  1.0 ) THEN
      IF (nice_use > 1) THEN
        !     ! Include effect of leads.
        rhokh_1_sea(i,j) = rhostar(i,j) * ch_ssi(i,j) * vshr_ssi(i,j)
        DO n = 1,nice_use
          rhokh_1_sice_ncats(i,j,n) = rhostar(i,j) * ch_ice(i,j,n) * vshr_ssi(i,j)
          IF (ice_fract(i,j) /= 0.0) THEN
            rhokh_1_sice(i,j) = rhokh_1_sice(i,j) +                           &
                (rhokh_1_sice_ncats(i,j,n) * ice_fract_cat(i,j,n)             &
                         / ice_fract(i,j))
          END IF
        END DO
        rhokm_ssi(i,j) = ((1.0 - ice_fract(i,j)) * rhostar(i,j)               &
                       *cd_ssi(i,j) * vshr_ssi(i,j))                          &
                      + (ice_fract(i,j) * rhostar(i,j)                        &
                       *cd_ice_mean(i,j) * vshr_ssi(i,j))
      ELSE
        !     ! Original scheme without leads.
        rhokh_1_sice_ncats(i,j,1) = rhostar(i,j) * ch_ssi(i,j) * vshr_ssi(i,j)
        rhokh_1_sice(i,j) = rhokh_1_sice_ncats(i,j,1)
        rhokm_ssi(i,j) = rhostar(i,j) * cd_ssi(i,j) * vshr_ssi(i,j)
      END IF
      !                                                          ! P243.124
      rhokm_ssi_nohalo(i,j) = rhokm_ssi(i,j)
      !
      !                                                           ! P243.125
    ELSE
      rhokm_ssi(       i,j) = 0.0
      rhokm_ssi_nohalo(i,j) = 0.0
      rhokh_1_sea(    i,j)    = 0.0
      rhokh_1_sice_ncats(    i,j,:) = 0.0
      rhokh_1_sice(    i,j) = 0.0
    END IF
  END DO
END DO
!$OMP END DO


!$OMP END PARALLEL

! Sea
lh0 = lc
dzssi(:) = dzsea
array_emis(:) = emis_sea


IF (nice_use > 1) THEN
  ! Include effect of leads
  sea_point = 1.0
  CALL sf_flux (                                                              &
     ssi_pts,sea_pts,fssi_ij,                                                 &
     ssi_index,sea_index,                                                     &
     array_zero_int,0,canhc_sea,dzssi,hcons_sea,                              &
     qs1,qstar_sea,qw_1,radnet_sea,                                           &
     array_one * beta_evap,rhokh_1_sea,array_false,array_zero,                &
     sea_frac,timestep,tl_1,tstar_sea,tstar_sea,                              &
     array_zero,array_zero,z0h_sea,z0msea,array_zero,z1_tq_sea,lh0,           &
     array_emis,array_one,                                                    &
     seasalinityfactor,array_zero,array_one,l_vegdrag_ssi,                    &
     fqw_1,ftl_1,                                                             &
     alpha1_sea,ashtf_prime_sea,e_sea,epot_sea,h_sea,                         &
     rhokm_1_sea,vshr_ssi,tau_sea,                                            &
     dtstar_sea,sea_point, sf_diag                                            &
     )

ELSE
  ! Do not include leads
  sea_point = 1.0
  CALL sf_flux (                                                              &
       ssi_pts,sea_pts,fssi_ij,                                               &
       ssi_index,sea_index,                                                   &
       array_zero_int,0,canhc_sea,dzssi,hcons_sea,                            &
       qs1,qstar_sea,qw_1,radnet_sea,                                         &
       array_one * beta_evap,rhokh_1_sice,array_false,array_zero,             &
       sea_frac,timestep,tl_1,tstar_sea,tstar_sea,                            &
       array_zero,array_zero,z0h_sea,z0msea,array_zero,z1_tq_sea,lh0,         &
       array_emis,array_one,                                                  &
       seasalinityfactor,array_zero,array_one,l_vegdrag_ssi,                  &
       fqw_1,ftl_1,                                                           &
       alpha1_sea,ashtf_prime_sea,e_sea,epot_sea,h_sea,                       &
       rhokm_1_sice,vshr_ssi,tau_sea,                                         &
       dtstar_sea,sea_point, sf_diag                                          &
       )
END IF

! rhokm_1_sea and tau_sea no longer required
DEALLOCATE(rhokm_1_sea)
DEALLOCATE(rhokm_1_sice)
DEALLOCATE(tau_sea)


! Sea ice

lh0 = ls

dzdummy(:) = 2.0       ! as hconsdz_sice equals 2*kappa/dz
! So in sf_flux, ashtf = 2 * hconsdz_sice/dzdummy, this will then give
!        ashtf = hconsdz_sice as required
array_emis(:) = emis_sice

IF (nice_use >  1) THEN
  DO n = 1, nice_use
    sea_point = 0.0
    CALL sf_flux (                                                            &
     ssi_pts,sice_pts_ncat(n),fssi_ij,                                        &
     ssi_index,sice_index_ncat(:,n),                                          &
     array_zero_int,0,array_zero,dzdummy,hconsdz_sice(:,:,n),                 &
     qs1,qstar_ice_cat(:,:,n),qw_1,radnet_sice(:,:,n),array_one,              &
     rhokh_1_sice_ncats(:,:,n),array_false,array_zero,                        &
     sice_frac_ncat(:,n),timestep,tl_1,ti_cat(:,:,n),                         &
     tstar_sice_cat(:,:,n),                                                   &
     array_zero,array_zero,z0h_ice(:,:,n),z0m_ice(:,:,n),array_zero,          &
     z1_tq_sea,lh0,                                                           &
     array_emis,array_one,                                                    &
     1.0,array_zero,array_one,l_vegdrag_ssi,                                  &
     fqw_1,ftl_1,                                                             &
     alpha1_sice(:,:,n),ashtf_prime(:,:,n),fqw_ice(:,:,n),epot_ice,           &
     ftl_ice(:,:,n),                                                          &
     rhokm_1_sice_ncats(:,:,n),vshr_ssi,tau_ice(:,:,n),                       &
     dtstar_sice(:,:,n),sea_point, sf_diag                                    &
     )
  END DO
ELSE
  sea_point = 0.0
  CALL sf_flux (                                                              &
   ssi_pts,sice_pts_ncat(1),fssi_ij,                                          &
   ssi_index,sice_index_ncat(:,1),                                            &
   array_zero_int,0,array_zero,dzdummy,hconsdz_sice(:,:,1),                   &
   qs1,qstar_ice,qw_1,radnet_sice(:,:,1),array_one,                           &
   rhokh_1_sice_ncats(:,:,1),array_false,array_zero,                          &
   sice_frac_ncat(:,1),timestep,tl_1,ti,tstar_sice_cat(:,:,1),                &
   array_zero,array_zero,z0h_ice(:,:,1),z0m_ice(:,:,1),array_zero,            &
   z1_tq_sea,lh0,                                                             &
   array_emis,array_one,                                                      &
   1.0,array_zero,array_one,l_vegdrag_ssi,                                    &
   fqw_1,ftl_1,                                                               &
   alpha1_sice(:,:,1),ashtf_prime(:,:,1),fqw_ice(:,:,1),epot_ice,             &
   ftl_ice(:,:,1),                                                            &
   rhokm_1_sice_ncats(:,:,1),vshr_ssi,tau_ice(:,:,1),                         &
   dtstar_sice(:,:,1),sea_point, sf_diag                                      &
   )
END IF

! rhokm_1_sice, rhokm_1_sive_ncats and tau_ice no longer required
DEALLOCATE(rhokm_1_sice_ncats)
DEALLOCATE(tau_ice)

CALL stdev1 (                                                                 &
   ssi_pts,sea_pts,ssi_index,sea_index,fssi_ij,                               &
   bq_1,bt_1,e_sea,h_sea,rhokm_ssi_nohalo,                                    &
   rhostar,vshr_ssi,z0msea,z1_tq_ctile,sea_frac,                              &
   q1_sd,t1_sd                                                                &
   )
!
IF (nice_use >  1) THEN
  DO n = 1,nice_use
    CALL stdev1 (                                                             &
     ssi_pts,sice_pts_ncat(n),ssi_index,sice_index_ncat(:,n),fssi_ij,         &
     bq_1,bt_1,fqw_ice(:,:,n),ftl_ice(:,:,n),rhokm_ssi_nohalo,                &
     rhostar,vshr_ssi,z0m_ice(:,:,n),z1_tq_ctile,                             &
     sice_frac_ncat(:,n),                                                     &
     q1_sd,t1_sd                                                              &
     )
  END DO
ELSE
  CALL stdev1 (                                                               &
   ssi_pts,sice_pts_use,ssi_index,sice_index_use,fssi_ij,                     &
   bq_1,bt_1,fqw_ice(:,:,1),ftl_ice(:,:,1),rhokm_ssi_nohalo,                  &
   rhostar,vshr_ssi,z0m_ice(:,:,1),z1_tq_ctile,                               &
   sice_frac_ncat(:,1),                                                       &
   q1_sd,t1_sd                                                                &
   )

END IF
!
!-----------------------------------------------------------------------
!  4.6 For sea points, calculate the wind mixing energy flux and the
!      sea-surface roughness length on the P-grid, using time-level n
!      quantities.
!-----------------------------------------------------------------------

SELECT CASE (iseasurfalg)
CASE (ip_ss_fixed, ip_ss_surf_div)
  SELECT CASE (i_high_wind_drag)
  CASE (ip_hwdrag_limited)
    ! Limit the neutral drag coefficient at 10m by calculating the
    ! equivalent roughness length and capping the roughness length.
    ! The equivalent wind speed will depend on the value of Charnock's
    ! coefficient.
    z0msea_max = z_10m / ( EXP(vkman / SQRT(cdn_max_sea) ) - 1.0)
  CASE (ip_hwdrag_reduced_v1)
    ! Calculate the maximum roughness length and the high-wind limit
    ! to limit the drag coefficient.
    z0msea_max = z_10m / ( EXP(vkman / SQRT(cdn_max_sea)) - 1.0)
  END SELECT
END SELECT

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, tau, u10n, cdn_lim_loc)            &
!$OMP SHARED(tdims, flandg, rhokm_ssi, vshr_ssi, ice_fract, rhostar, cd_sea,  &
!$OMP        sf_diag, iseasurfalg, charnock, g, z0msea, l_spec_z0, z0m_scm,   &
!$OMP        z0h_scm, z0h_sea, z0hsea,                                        &
!$OMP        i_high_wind_drag, cdn_max_sea, cdn_hw_sea,                       &
!$OMP        u_cdn_max, u_cdn_hw, z0msea_max, charnock_w)

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end

    IF (flandg(i,j) <  1.0) THEN
      tau = rhokm_ssi(i,j) * vshr_ssi(i,j)             ! P243.130
      IF (ice_fract(i,j) >  0.0)                                              &
        tau = rhostar(i,j) * cd_sea(i,j)                                      &
          * vshr_ssi(i,j) * vshr_ssi(i,j)

      IF (sf_diag%sfme)                                                       &
        sf_diag%fme(i,j) = (1.0 - ice_fract(i,j)) * tau * SQRT(tau / rhosea)
      !                                                            ! P243.96
      !     Recalculate the momentum roughness length under the older
      !     non-interactive treatments.
      SELECT CASE (iseasurfalg)
        !
      CASE (ip_ss_fixed, ip_ss_surf_div)

        ! Limit Z0MSEA to 0.154m for TAU very small
        IF ( rhostar(i,j) > EPSILON(0.0) ) THEN
          z0msea(i,j) = 1.54e-6 / (SQRT(tau / rhostar(i,j)) + 1.0e-5)         &
                      +  (charnock / g) * (tau / rhostar(i,j))
          z0msea(i,j) = MAX ( z0hsea , z0msea(i,j) )
          !                                       ... P243.B6 (Charnock formula)
          !                    TAU/RHOSTAR is "mod VS squared", see eqn P243.131
        ELSE
          z0msea(i,j) = z0hsea
        END IF
        !
      CASE (ip_ss_surf_div_coupled)

        ! Limit Z0MSEA to 0.154m for TAU very small
        IF ( rhostar(i,j) > EPSILON(0.0) ) THEN
          z0msea(i,j) = 1.54e-6 / (SQRT(tau / rhostar(i,j)) + 1.0e-5)         &
                      +  (charnock_w(i,j) / g) * (tau / rhostar(i,j))
          z0msea(i,j) = MAX ( z0hsea , z0msea(i,j) )
          !                                       ... P243.B6 (Charnock formula)
          !                    TAU/RHOSTAR is "mod VS squared", see eqn P243.131
        ELSE
          z0msea(i,j) = z0hsea
        END IF
        !
      CASE DEFAULT
        !
        !       The momentum roughness length has already been calculated
        !       within the iteration for the Obukhov length.
        !
      END SELECT
      SELECT CASE (i_high_wind_drag)
      CASE (ip_hwdrag_limited)
        z0msea(i,j) = MIN(z0msea(i,j), z0msea_max)
      CASE (ip_hwdrag_reduced_v1)
        !         Calculate 10-m neutral wind based on current stress
        u10n = (SQRT(tau / rhostar(i,j)) / vkman) *                           &
               LOG (1.0 + z_10m / z0msea(i,j))
        !         Determine a limiting value of cd
        IF (u10n <= u_cdn_max) THEN
          cdn_lim_loc = cdn_max_sea
        ELSE IF ( (u10n > u_cdn_max) .AND. (u10n < u_cdn_hw) ) THEN
          cdn_lim_loc = cdn_max_sea - (cdn_max_sea - cdn_hw_sea) *            &
            (u10n - u_cdn_max) / (u_cdn_hw - u_cdn_max)
        ELSE
          cdn_lim_loc = cdn_hw_sea
        END IF
        !         Reset the roughness length consistently, leaving aside very
        !         light winds.
        IF (u10n > 1.0)                                                       &
          z0msea(i,j) = MIN(z0msea(i,j),                                      &
            z_10m / ( EXP(vkman / SQRT(cdn_lim_loc) ) - 1.0))
      END SELECT

    END IF

  END DO
END DO
!$OMP END DO


IF ( l_spec_z0 ) THEN
  ! Check for prescribed surface roughness lengths specified in SCM
  ! NAMELIST.  If specified in the &INPROF then they will be used
  ! instead of Model calculated values
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF ( z0m_scm(i,j) > 0.0 ) THEN
        ! Set z0m from SCM namelist for sea points
        z0msea(i,j)  = z0m_scm(i,j)
      END IF
      IF ( z0h_scm(i,j) > 0.0 ) THEN
        ! Set z0h from SCM namelist for sea points
        z0h_sea(i,j) = z0h_scm(i,j)
      END IF
    END DO
  END DO
!$OMP END DO
END IF
!$OMP END PARALLEL

! The following lines have been commented out since unless accompanied by
! similar changes for cd_std_sea and v_s_std_sea, which do not feature in
! JULES they will cause an inconsistency in the diagnosis of 10-m winds.
! The question of how winds and temperatures should be diagnosed over
! sea ice is being reconsidered.
!--!-----------------------------------------------------------------------
!--! If sea ice is present then set RECIP_L_MO_SEA to its value over
!--! the ice, a long-standing choice for the screen diagnostics.
!--! Note that RECIP_L_MO_SEA is also used in BDY_EXPL2 to diagnose
!--! shear-dominated boundary layer types.  To preserve bit
!--! reproducibility when screen diagnostics are switched on,
!--! this change has been moved outside the if-test on stash logicals
!--!-----------------------------------------------------------------------
!--!DO j=1,rows
!--!  DO i=1,row_length
!--!    IF (flandg(i,j) <  1.0 .AND. ice_fract(i,j) >  0.0 ) THEN
!--!      recip_l_mo_sea(i,j) = recip_l_mo_ice(i,j)
!--!    END IF
!--!  END DO
!--!END DO

IF (sf_diag%l_t10m .OR. sf_diag%l_q10m) THEN
  ALLOCATE(chr10m_sice(tdims%i_start:tdims%i_end,                             &
                       tdims%j_start:tdims%j_end,nice_use))
  ALLOCATE(chr10m_sea(tdims%i_start:tdims%i_end,                              &
                      tdims%j_start:tdims%j_end))
ELSE
  ALLOCATE(chr10m_sice(1,1,1))
  ALLOCATE(chr10m_sea(1,1))
END IF

! Sea and sea-ice (leads ignored)
! The unloading rate for plant canopies is tested since we are
! using the GBM drag coefficient at coastal points for simplicity in the
! calculation of unloading.
IF (sf_diag%su10 .OR. sf_diag%sv10 .OR. sf_diag%sq1p5 .OR.                    &
    sf_diag%st1p5 .OR. sf_diag%suv10m_n .OR.                                  &
    sf_diag%l_t10m .OR. sf_diag%l_q10m .OR.                                   &
    l_cdr10m_snow .OR.                                                        &
    (IScrnTDiag == IP_ScrnDecpl2) .OR.                                        &
    (IScrnTDiag == IP_ScrnDecpl3)) THEN
  chr1p5m_sea(:,:) = 0.0
  CALL sfl_int (                                                              &
     ssi_pts,sea_pts,l_cdr10m_snow,sea_index,ssi_index,fssi_ij,               &
     vshr_ssi,cd_std_sea,cd_sea,ch_sea,                                       &
     sea_frac,                                                                &
     z0msea,z0msea,z0h_sea,                                                   &
     recip_l_mo_sea,                                                          &
     v_s_sea,v_s_std_sea,                                                     &
     z1_uv,z1_tq_ctile,db_sea,                                                &
     sf_diag,                                                                 &
     cdr10m,sf_diag%cdr10m_n,sf_diag%cd10m_n,chr1p5m_sea,chr10m_sea           &
     )
  !
  chr1p5m_sice(:,:,:) = 0.0        ! Initialise
  DO n = 1,nice_use
    CALL sfl_int (                                                            &
       ssi_pts,sice_pts_ncat(n),l_cdr10m_snow,sice_index_ncat(:,n),ssi_index, &
       fssi_ij,vshr_ssi,cd_std_ice(:,n),cd_ice(:,:,n),ch_ice(:,:,n),          &
       sice_frac_ncat(:,n),                                                   &
       z0m_ice(:,:,n),z0m_ice(:,:,n),z0h_ice(:,:,n),                          &
       recip_l_mo_ice(:,:,n),                                                 &
       v_s_ice(:,:,n),v_s_std_ice(:,n),                                       &
       z1_uv,z1_tq_ctile,db_ice(:,:,n),                                       &
       sf_diag,                                                               &
       cdr10m,sf_diag%cdr10m_n,sf_diag%cd10m_n,chr1p5m_sice(:,:,n),           &
       chr10m_sice(:,:,n)                                                     &
       )
  END DO
  !
  IF (sf_diag%l_t10m .OR. sf_diag%l_q10m) THEN
    IF (nice_use > 1) THEN
      ! Take account of leads, and use multiple thickness categories,
      ! if these are present in CICE.
      sf_diag%chr10m(:,:) = ((1.0 - ice_fract(:,:)) * chr10m_sea(:,:))
      DO n = 1,nice_use
        DO k = 1,sice_pts_ncat(n)
          l = sice_index_ncat(k,n)
          j=(ssi_index(l) - 1) / t_i_length + 1
          i = ssi_index(l) - (j-1) * t_i_length
          sf_diag%chr10m(i,j) = sf_diag%chr10m(i,j)                           &
                    + (sice_frac_ncat(l,n) * chr10m_sice(i,j,n))
        END DO
      END DO
    ELSE
      ! Do not include leads or thickness categories.

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,ice_fract,sf_diag,chr10m_sice,chr10m_sea)
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          IF (ice_fract(i,j) >  0.0) THEN
            sf_diag%chr10m(i,j) = chr10m_sice(i,j,1)
          ELSE
            sf_diag%chr10m(i,j) = chr10m_sea(i,j)
          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF
  END IF
  DEALLOCATE(chr10m_sice)
  DEALLOCATE(chr10m_sea)

  IF (nice_use > 1) THEN
    ! Take account of leads, and use multiple thickness categories,
    ! if these are present in CICE.
    chr1p5m_ssi_mean(:,:) = ((1.0 - ice_fract(:,:)) * chr1p5m_sea(:,:))
    DO n = 1,nice_use
      DO k = 1,sice_pts_ncat(n)
        l = sice_index_ncat(k,n)
        j=(ssi_index(l) - 1) / t_i_length + 1
        i = ssi_index(l) - (j-1) * t_i_length
        chr1p5m_ssi_mean(i,j) = chr1p5m_ssi_mean(i,j)                         &
                    + (sice_frac_ncat(l,n) * chr1p5m_sice(i,j,n))
      END DO
    END DO
  ELSE
    ! Do not include leads or thickness categories.
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,ice_fract,chr1p5m_ssi_mean,chr1p5m_sice,chr1p5m_sea)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        IF (ice_fract(i,j) >  0.0) THEN
          chr1p5m_ssi_mean(i,j) = chr1p5m_sice(i,j,1)
        ELSE
          chr1p5m_ssi_mean(i,j) = chr1p5m_sea(i,j)
        END IF
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF


END IF
!
!-----------------------------------------------------------------------
! GBM diagnstic calculations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculate effective roughness lengths, orographic blending heights
! and gridbox-average Richardson numbers.
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j)                                   &
!$OMP SHARED(tdims, fz0, fz0h, rib, h_blend_orog, flandg,                    &
!$OMP        ice_fract, rib_sea, z0msea, z0h_sea)

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    fz0(i,j) = 0.0
    fz0h(i,j) = 0.0
    rib(i,j) = 0.0
    h_blend_orog(i,j) = h_blend_min
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    rib(i,j) = rib(i,j) +                                                     &
               (1.0 - flandg(i,j)) * (1.0 - ice_fract(i,j)) * rib_sea(i,j)
    fz0(i,j) = fz0(i,j) + (1.0 - flandg(i,j)) * (1.0 - ice_fract(i,j))        &
                    / (LOG(lb / z0msea(i,j))**2)
    fz0h(i,j) = fz0h(i,j) +                                                   &
                    (1.0 - flandg(i,j)) * (1.0 - ice_fract(i,j)) /            &
                    (LOG(lb / z0msea(i,j)) * LOG(lb / z0h_sea(i,j)))
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
!
IF (nice_use >  1) THEN
  DO n = 1,nice_use
    DO k = 1,sice_pts_ncat(n)
      l = sice_index_ncat(k,n)
      j=(ssi_index(l) - 1) / t_i_length + 1
      i = ssi_index(l) - (j-1) * t_i_length
      rib(i,j) = rib(i,j) +                                                   &
             (1.0 - flandg(i,j)) * sice_frac_ncat(l,n) * rib_ice(i,j,n)
      fz0(i,j) = fz0(i,j) + (1.0 - flandg(i,j))                               &
                    * sice_frac_ncat(l,n) /                                   &
                       (LOG(lb / z0m_ice(i,j,n))**2)
      fz0h(i,j) = fz0h(i,j) + (1.0 - flandg(i,j))                             &
                                *sice_frac_ncat(l,n) /                        &
                 (LOG(lb / z0m_ice(i,j,n)) * LOG(lb / z0h_ice(i,j,n)))
    END DO
  END DO
ELSE
  rib(:,:) = rib(:,:) +                                                       &
        (1.0 - flandg(:,:)) * ice_fract(:,:) * rib_ice(:,:,1)
  fz0(:,:) = fz0(:,:) + (1.0 - flandg(:,:))                                   &
                    * ice_fract(:,:) /                                        &
                       (LOG(lb / z0m_ice(:,:,1))**2)
  fz0h(:,:) = fz0h(:,:) + (1.0 - flandg(:,:))                                 &
                                *ice_fract(:,:) /                             &
                 (LOG(lb / z0m_ice(:,:,1)) * LOG(lb / z0h_ice(:,:,1)))
END IF
!
DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    rib(i,j) = rib(i,j) + flandg(i,j) * tile_frac(l,n) * rib_surft(l,n)
    fz0(i,j) = fz0(i,j)                                                       &
      + flandg(i,j) * tile_frac(l,n) / (LOG(lb / z0m_surft(l,n))**2)
    fz0h(i,j) = fz0h(i,j) + flandg(i,j) * tile_frac(l,n) /                    &
                    (LOG(lb / z0m_surft(l,n)) * LOG(lb / z0h_surft(l,n)))
  END DO
END DO

IF (sf_diag%l_rib_surft) THEN
  sf_diag%rib_surft = rib_surft
END IF

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(j,i,l)                                   &
!$OMP SHARED(tdims,fz0,fz0h,z0m_eff,z0h_eff,sf_diag,land_pts,land_index,      &
!$OMP        t_i_length,z0h_gb_land,z0m_gb_land)
!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    z0m_eff(i,j) = lb * EXP( - SQRT(1.0 / fz0(i,j)) )
    z0h_eff(i,j) = lb * EXP( - SQRT(fz0(i,j)) / fz0h(i,j))
  END DO
END DO
!$OMP END DO

IF (sf_diag%l_z0m_gb) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      sf_diag%z0m_gb(i,j) = z0m_eff(i,j)
    END DO
  END DO
!$OMP END DO
END IF

!$OMP DO SCHEDULE(STATIC)
DO l = 1,land_pts
  j=(land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length
  z0m_gb_land(l) = z0m_eff(i,j)
  z0h_gb_land(l) = z0h_eff(i,j)
END DO
!$OMP END DO
!$OMP END PARALLEL


IF (formdrag ==  effective_z0) THEN
  CALL sf_orog_gb(                                                            &
   land_pts,land_index,                                                       &
   fd_stab_dep,orog_drag_param,                                               &
   ho2r2_orog,rib,sil_orog,z0m_gb_land,z1_uv,                                 &
   h_blend_orog,z0m_eff,sf_diag,z0h_gb_land,z0h_eff                           &
   )
END IF

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(j,i)                                     &
!$OMP SHARED(sf_diag,tdims,z0h_eff,rhokm_1,flandg,rhokm_land,rhokm_ssi)
IF (sf_diag%l_z0h_eff_gb) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      sf_diag%z0h_eff_gb(i,j) = z0h_eff(i,j)
    END DO
  END DO
!$OMP END DO
END IF

!-----------------------------------------------------------------------
! Set grid-box surface exchange coefficients
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    rhokm_1(i,j) = flandg(i,j) * rhokm_land(i,j) +                            &
                       (1.0 - flandg(i,j)) * rhokm_ssi(i,j)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (nice_use > 1) THEN
  ! Take account of leads, and use multiple thickness categories,
  ! if these are present in CICE.
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j) &
!$OMP SHARED(tdims,flandg,rhokh_gb,ice_fract,rhokh_1_sea)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF ( flandg(i,j) <  1.0) THEN
        rhokh_gb(i,j) = (1.0 - flandg(i,j)) * (1.0 - ice_fract(i,j))          &
                         * rhokh_1_sea(i,j)
      ELSE
        rhokh_gb(i,j) = 0.0
      END IF
    END DO
  END DO
!$OMP  END PARALLEL DO
  DO n = 1,nice_use
    DO k = 1,sice_pts_ncat(n)
      l = sice_index_ncat(k,n)
      j=(ssi_index(l) - 1) / t_i_length + 1
      i = ssi_index(l) - (j-1) * t_i_length
      rhokh_gb(i,j) = rhokh_gb(i,j) +                                         &
              (sice_frac_ncat(l,n) * rhokh_1_sice_ncats(i,j,n))
    END DO
  END DO
ELSE
  ! Do not include leads or thickness categories.
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j) &
!$OMP SHARED(tdims,flandg,rhokh_gb,rhokh_1_sice)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF ( flandg(i,j) <  1.0) THEN
        rhokh_gb(i,j) = (1.0 - flandg(i,j)) * rhokh_1_sice(i,j)
      ELSE
        rhokh_gb(i,j) = 0.0
      END IF
    END DO
  END DO
!$OMP  END PARALLEL DO
END IF


DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    rhokh_gb(i,j) = rhokh_gb(i,j)                                             &
      + flandg(i,j) * tile_frac(l,n) * rhokh_1(l,n)
  END DO
END DO


!-----------------------------------------------------------------------
! Calculate scaling parameters required for non-local BL scheme
!-----------------------------------------------------------------------

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(bl_stress_grad, i, j)                 &
!$OMP SHARED(tdims, u_s, flandg, cd_land, vshr_land, cd_ssi, vshr_ssi, zh,    &
!$OMP        fb_surf, g, bt_1, ftl_1, bq_1, fqw_1, rhostar)                   &
!$OMP SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    u_s(i,j) = SQRT( flandg(i,j) * cd_land(i,j) * vshr_land(i,j) *            &
                     vshr_land(i,j) + (1.0 - flandg(i,j)) * cd_ssi(i,j) *     &
                     vshr_ssi(i,j) * vshr_ssi(i,j) )
    !   !------------------------------------------------------
    !   ! Limit the explicitly calculated surface stress,
    !   ! used to scale the non-local parametrizations,
    !   ! such that the implied stress gradient across the BL
    !   ! is less than Max_Stress_Grad.
    !   !------------------------------------------------------
    bl_stress_grad = u_s(i,j) * u_s(i,j) / zh(i,j)
    IF (bl_stress_grad > max_stress_grad) THEN
      u_s(i,j) = SQRT(zh(i,j) * max_stress_grad)
    END IF

    fb_surf(i,j) = g * ( bt_1(i,j) * ftl_1(i,j) +                             &
                     bq_1(i,j) * fqw_1(i,j) ) / rhostar(i,j)
  END DO
END DO
!$OMP END PARALLEL DO


! Calculate parameters required for CLASSIC aerosol
! Note that the values for cd_std and ch_surft are now those
! set using the "aerosol roughness length"
CALL sf_aero (                                                                &
 land_pts,nsurft,land_index,surft_index,surft_pts,                            &
 l_aero_classic,l_dust,l_dust_diag,                                           &
 flandg,tile_frac,pstar,rhostar,tstar,vshr_land,vshr_ssi,                     &
 cd_ssi,ch_ssi,cd_std_classic,ch_surft_classic,                               &
 rho_aresist,aresist,resist_b,rho_aresist_surft,aresist_surft,                &
 resist_b_surft,r_b_dust,cd_std_dust                                          &
 )

! set these roughness lengths which otherwise are unspecified
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,z0mssi,z0msea,z0hssi,z0h_sea)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    z0mssi(i,j) = z0msea(i,j)
    z0hssi(i,j) = z0h_sea(i,j)
  END DO
END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
! Atmospherically determined variables for the snow scheme.
! (These calculations are placed here to use the full GBM wind speed
! with coastal tiling and for possible later convenience when integrated
! with sea-ice.)
!-----------------------------------------------------------------------
!
DO n = 1, nsurft
  IF (cansnowtile(n)) THEN
    unload_backgrnd_pft(:,n) = unload_rate_cnst(n)
    IF ( unload_rate_u(n) /= 0.0 ) THEN
      DO k = 1, surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j-1) * t_i_length
        !       Use GBM value of CD for simplicity. The complexity of
        !       distinguishing coastal points does not seem justified.
        ws10 = cdr10m(i,j) * vshr_land(i,j)
        unload_backgrnd_pft(l,n) = unload_backgrnd_pft(l,n) +                 &
          unload_rate_u(n) * ws10
      END DO
    END IF
  END IF
END DO

first_call = .false.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_exch_cbl
END MODULE sf_exch_cbl_mod
