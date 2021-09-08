! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE sf_impl2_cbl_mod

USE sf_melt_mod, ONLY: sf_melt
USE screen_tq_cbl_mod, ONLY: screen_tq_cbl
!C!USE sf_evap_mod, ONLY: sf_evap
!USE im_sf_pt2_mod, ONLY: im_sf_pt2
USE im_sf_pt2_cbl_mod, ONLY: im_sf_pt2_cbl
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_IMPL2_MOD'

CONTAINS
!  SUBROUTINE SF_IMPL2-----------------------------------------------
!
!  Purpose: Calculate implicit correction to surface fluxes of heat,
!           moisture and momentum to be used by the unconditionally
!           stable and non-oscillatory BL numerical solver.  Also
!           calculates screen level temperature and humidity as well
!           as 10 m winds.
!
!--------------------------------------------------------------------
!    Arguments :-
SUBROUTINE sf_impl2_cbl (                                                     &
! IN values defining field dimensions and subset to be processed :
 land_pts,land_index,nice,nice_use,nsurft,surft_index,surft_pts,              &
 sm_levels,canhc_surft,canopy,flake,smc_soilt,tile_frac,wt_ext_surft,         &
 fland,flandg,lq_mix_bl,                                                      &
! IN sea/sea-ice data :
 ice_fract,ice_fract_ncat,k_sice,u_0,v_0,                                     &
! IN everything not covered so far :
 pstar,lw_down,sw_surft,                                                      &
 t_soil_soilt,qw_1,tl_1,u_1,v_1,rhokm_u_1,rhokm_v_1,r_gamma,                  &
 gamma1,gamma2,alpha1,alpha1_sea,alpha1_sice,                                 &
 ashtf_prime,ashtf_prime_sea,ashtf_prime_surft,                               &
 dtrdz_charney_grid_1,du_1,dv_1,                                              &
 fraca,resfs,resft,rhokh,rhokh_surft,rhokh_sice,rhokh_sea,z1,                 &
 z0hssi,z0mssi,z0h_surft,z0m_surft,cdr10m_u,cdr10m_v,                         &
 chr1p5m,chr1p5m_sice,ctctq1,                                                 &
 dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,                         &
 l_correct,flandg_u,flandg_v,                                                 &
 emis_surft,ti,snow_surft,                                                    &
! IN variables used to calculate cooling at the screen level
 l_co2_interactive, co2_mmr, co2_3d,rho1, f3_at_p, ustargbm,                  &
! INOUT data :
 epot_surft,fqw_ice,ftl_ice,dtstar_surft,dtstar_sea,dtstar_sice,              &
 tstar_sice_cat,tstar_ssi,tstar_surft,tstar_sea,radnet_sice,fqw_surft,        &
 fqw_1,ftl_1,ftl_surft,olr,taux_land,taux_ssi,tauy_land,tauy_ssi,             &
 TScrnDcl_SSI,TScrnDcl_SURFT,tStbTrans,sf_diag,                               &
 taux_land_star,tauy_land_star,taux_ssi_star, tauy_ssi_star,                  &
! OUT Diagnostic not requiring STASH flags :
 ecan,ei_surft,esoil_surft,sea_ice_htf,surf_ht_flux,                          &
 surf_ht_flux_land,surf_ht_flux_sice,surf_htf_surft,                          &
! OUT data required elsewhere in UM system :
 tstar,tstar_land,tstar_sice,le_surft,radnet_surft,e_sea,h_sea,               &
 taux_1,tauy_1,ecan_surft,ei,ei_sice,esoil_soilt,ext_soilt,                   &
 snowmelt, melt_surft,rhokh_mix,error,                                        &
!CABLE_LSM: see ACCESS-CM2 version for extra vars here
 air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                            &
 ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                          &
 soil_cbl )

USE csigma,                   ONLY:                                           &
 sbcon

USE planet_constants_mod,     ONLY: cp

USE atm_fields_bounds_mod,    ONLY:                                           &
  tdims, udims, vdims, pdims, tdims_s, udims_s, vdims_s, pdims_s

USE theta_field_sizes,        ONLY:                                           &
  t_i_length, t_j_length

USE jules_surface_mod,        ONLY:                                           &
  l_aggregate, l_flake_model, ls, ip_scrndecpl2, ip_scrndecpl3

USE jules_radiation_mod,      ONLY:                                           &
  l_dolr_land_black

USE jules_sea_seaice_mod,     ONLY:                                           &
  l_tstar_sice_new, emis_sea, emis_sice, l_use_dtstar_sea

USE jules_snow_mod,           ONLY:                                           &
  nsmax, rho_snow_const, cansnowtile, l_snowdep_surf, l_snow_nocan_hc

USE jules_vegetation_mod,     ONLY:                                           &
  can_model

USE jules_surface_types_mod,  ONLY:                                           &
  lake

USE lake_mod,                 ONLY:                                           &
  surf_ht_flux_lake_ij, lake_h_ice_gb
USE sf_diags_mod, ONLY: strnewsfdiag
USE timestep_mod,             ONLY:                                           &
  timestep

USE fluxes,                   ONLY:                                           &
  anthrop_heat_surft, surf_ht_store_surft, sw_sicat

USE ancil_info,               ONLY:                                           &
  ssi_pts, sea_pts, sice_pts, sice_pts_ncat, ssi_index, sea_index,            &
  sice_index,sice_index_ncat, fssi_ij, sea_frac, sice_frac, sice_frac_ncat,   &
  nsoilt

USE jules_surface_mod,        ONLY:                                           &
  iscrntdiag, l_neg_tstar

USE c_elevate,                ONLY: lw_down_elevcorr_surft

USE prognostics,              ONLY:                                           &
  nsnow_surft, snowdepth_surft

USE jules_mod,                ONLY:                                           &
  snowdep_surft

USE water_constants_mod,      ONLY:                                           &
  lc, lf, rho_ice, tfs

USE solinc_data,              ONLY:                                           &
  sky, l_skyview

USE parkind1,                 ONLY:                                           &
  jprb, jpim
!CABLE_LSM:Make avail call to CABLE implicit version
USE cable_implicit_main_mod,  ONLY: cable_implicit_main
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


USE yomhook,                  ONLY:                                           &
  lhook, dr_hook

USE jules_print_mgr,          ONLY:                                           &
  jules_message, jules_print

IMPLICIT NONE
!--------------------------------------------------------------------
!  Inputs :-
! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.
!--------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
 land_pts    ! IN No of land points

! (c) Soil/vegetation/land surface parameters (mostly constant).
INTEGER, INTENT(IN) ::                                                        &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
!                                  !    point in ROW_LENGTH,ROWS is the
!                                  !    Ith land point.

INTEGER, INTENT(IN) ::                                                        &
 sm_levels                                                                    &
                             ! IN No. of soil moisture levels
,nsurft                                                                       &
                             ! IN No. of land tiles
,surft_index(land_pts,nsurft)                                                 &
                             ! IN Index of tile points
,surft_pts(nsurft)                                                            &
                             ! IN Number of tile points
,nice                                                                         &
                             ! IN Number of sea ice categories
,nice_use                    ! IN Number of sea ice categories used
                             !    fully in surface exchange

REAL, INTENT(IN) ::                                                           &
 canhc_surft(land_pts,nsurft)                                                 &
                             ! IN Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
,canopy(land_pts,nsurft)                                                      &
                             ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
,flake(land_pts,nsurft)                                                       &
                             ! IN Lake fraction.
,smc_soilt(land_pts,nsoilt)                                                   &
                             ! IN Available soil moisture (kg/m2).
,tile_frac(land_pts,nsurft)                                                   &
                             ! IN Tile fractions including
!                                  ! snow cover in the ice tile.
,wt_ext_surft(land_pts,sm_levels,nsurft)                                      &
!                                  ! IN Fraction of evapotranspiration
!                                  !    extracted from each soil layer
!                                  !    by each tile.
,fland(land_pts)                                                              &
                             ! IN Land fraction on land pts.
,flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
!                                  ! IN Land fraction on all pts.
,flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)                &
                             ! IN Land fraction on U grid.
,flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)                &
                             ! IN Land fraction on V grid.
,emis_surft(land_pts,nsurft)                                                  &
                             ! IN Emissivity for land tiles
,ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)                 &
                             ! IN   Sea-ice surface layer
                             !       temperature (K).
,snow_surft(land_pts,nsurft)
                             ! IN Lying snow on tiles (kg/m2)

! (d) Sea/sea-ice data.
REAL, INTENT(IN) ::                                                           &
 ice_fract(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! IN Fraction of gridbox covered by
!                            !     sea-ice (decimal fraction).
,ice_fract_ncat(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end,nice)                               &
                             ! IN Fraction of gridbox
!                            !  covered by sea-ice on categories.
,k_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)             &
                             ! IN sea ice effective conductivity in
!                             !     sfc layer on categories (W/m2/k)
,u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end)                     &
                             ! IN W'ly component of surface
!                                  !    current (m/s).
,v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                             ! IN S'ly component of surface
!                                  !    current (m/s).

! (f) Atmospheric + any other data not covered so far, incl control.

REAL, INTENT(IN) ::                                                           &
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)                   &
                             ! IN Surface pressure (Pascals).
,lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! IN Surface downward LW radiation
!                                  !    (W/m2).
,sw_surft(land_pts,nsurft)                                                    &
                             ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
,t_soil_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Soil temperatures (K).
,qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN Total water content
,tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN Ice/liquid water temperature
,u_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end)             &
                             ! IN W'ly wind component (m/s)
,v_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end)             &
                             ! IN S'ly wind component (m/s)
,rhokm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)               &
                             ! IN Exchange coefficients for
!                                  !    momentum (on U-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
,rhokm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                             ! IN Exchange coefficients for
!                                  !    momentum (on V-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")

REAL, INTENT(IN) :: r_gamma        ! IN implicit weight in level 1

REAL, INTENT(IN) ::                                                           &
 gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)                  &
                             ! weights for new BL solver
,gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL, INTENT(IN) ::                                                           &
 alpha1(land_pts,nsurft)                                                      &
                             ! IN Mean gradient of saturated
!                                  !    specific humidity with respect
!                                  !    to temperature between the
!                                  !    bottom model layer and tile
!                                  !    surfaces
,alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! IN ALPHA1 for sea.
,alpha1_sice(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                             ! IN ALPHA1 for sea-ice.
,ashtf_prime(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                             ! IN Adjusted SEB coefficient for
                             !    sea-ice
,ashtf_prime_sea(tdims%i_start:tdims%i_end,                                   &
                 tdims%j_start:tdims%j_end)                                   &
                             ! IN Adjusted SEB coefficient for
                             !    sea
,ashtf_prime_surft(land_pts,nsurft)                                           &
                             ! IN Adjusted SEB coefficient for
                             !    land tiles.
,dtrdz_charney_grid_1(pdims%i_start:pdims%i_end,                              &
                      pdims%j_start:pdims%j_end)                              &
!                                  ! IN -g.dt/dp for model layers.
,fraca(land_pts,nsurft)                                                       &
                             ! IN Fraction of surface moisture
!                                  !    flux with only aerodynamic
!                                  !    resistance for snow-free land
!                                  !    tiles.
,resfs(land_pts,nsurft)                                                       &
                             ! IN Combined soil, stomatal
!                                  !    and aerodynamic resistance
!                                  !    factor for fraction (1-FRACA) of
!                                  !    snow-free land tiles.
,resft(land_pts,nsurft)                                                       &
                             ! IN Total resistance factor.
!                                  !    FRACA+(1-FRACA)*RESFS for
!                                  !    snow-free land, 1 for snow.
,rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN Grid-box surface exchange
!                                  !     coefficients
!                                  !    (not used for JULES)
,rhokh_surft(land_pts,nsurft)                                                 &
                             ! IN Surface exchange coefficients
!                                  !    for land tiles
,rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
                                                        nice_use)             &
!                             ! IN Surface exchange coefficients
!                                  !    for sea-ice
,rhokh_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                             ! IN Surface exchange coefficients
!                                  !    for sea

REAL, INTENT(IN) ::                                                           &
z1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                       &
                            ! IN Height of lowest level (i.e.
!                                  !    height of middle of lowest
!                                  !    layer).
,z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
,z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! IN Roughness lengths over sea (m)
,z0h_surft(land_pts,nsurft)                                                   &
                             ! IN Tile roughness lengths for heat
!                                  !    and moisture (m).
,z0m_surft(land_pts,nsurft)                                                   &
                             ! IN Tile roughness lengths for
!                                  !    momentum.
,cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)                &
                             ! IN Ratio of CD's reqd for
!                                  !    calculation of 10 m wind. On
!                                  !    U-grid; comments as per RHOKM.
,cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)                &
                             ! IN Ratio of CD's reqd for
!                                  !    calculation of 10 m wind. On
!                                  !    V-grid; comments as per RHOKM.
,chr1p5m(land_pts,nsurft)                                                     &
                             ! IN Ratio of coefffs for calculation
!                                  !    of 1.5m temp for land tiles.
,chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)            &
!                                  ! IN CHR1P5M for sea and sea-ice
!                                  !    (leads ignored).
,cq_cm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)               &
                             ! IN Coefficient in U tri-diagonal
!                                  !    implicit matrix
,cq_cm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)               &
                             ! IN Coefficient in V tri-diagonal
!                                  !    implicit matrix
,du_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end)            &
                             ! IN Level 1 increment to u wind
!                                  !    field
,dv_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end)            &
                             ! IN Level 1 increment to v wind
!                                  !    field
,ctctq1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)                  &
,dqw1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)                  &
,dtl1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)                  &
,du_star1(udims_s%i_start:udims_s%i_end,                                      &
          udims_s%j_start:udims_s%j_end)                                      &
,dv_star1(vdims_s%i_start:vdims_s%i_end,                                      &
          vdims_s%j_start:vdims_s%j_end)
!                                  ! IN Additional arrays needed by the
!                                  !    uncond stable BL numerical solver
! IN Additional variables for screen-level diagnostics
LOGICAL, INTENT(IN) :: l_co2_interactive
                                ! Flag for interactive 3-D CO2
REAL, INTENT(IN)    :: co2_mmr
                                ! Initial or fixed mass mixing ratio
                                ! of CO2
REAL, INTENT(IN)    ::                                                        &
  co2_3d(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end)
  ! 3-D field of CO2
REAL, INTENT(IN)    ::                                                        &
  rho1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
  ! Density on lowest level
REAL, INTENT(IN)    ::                                                        &
  f3_at_p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
  ! Coriolis parameter
REAL, INTENT(IN)    ::                                                        &
  uStarGBM(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
  ! GBM surface friction velocity

LOGICAL, INTENT(IN) ::                                                        &
 l_correct                   ! IN flag used by the new BL solver

LOGICAL, INTENT(IN) ::                                                        &
 lq_mix_bl                   ! TRUE if mixing ratios used in
                             ! boundary layer code

!--------------------------------------------------------------------
!  In/outs :-
!--------------------------------------------------------------------
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag
REAL, INTENT(INOUT) ::                                                        &
 epot_surft(land_pts,nsurft)                                                  &
                             ! INOUT surface tile potential
!                                  !    evaporation
,fqw_ice(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end,nice_use)                                  &
                             ! INOUT Surface FQW for sea-ice
,ftl_ice(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end,nice_use)                                  &
                             ! INOUT Surface FTL for sea-ice
,dtstar_surft(land_pts,nsurft)                                                &
                             ! INOUT Change in TSTAR over timestep
!                                  !     for land tiles
,dtstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! INOUT Change is TSTAR over timestep
!                                  !     for open sea
,dtstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)    &
                             ! INOUT Change is TSTAR over timestep
!                                  !     for sea-ice
,tstar_sice_cat(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end,nice_use)                           &
                             ! INOUT   Sea-ice sfc temperature (K).
,tstar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! INOUT Sea mean sfc temperature (K).
,tstar_surft(land_pts,nsurft)                                                 &
                             ! INOUT Surface tile temperatures
,tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! IN    Open sea sfc temperature (K).
,radnet_sice(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                             ! INOUT Surface net radiation on
!                                  !       sea-ice (W/m2)
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! INOUT Moisture flux between layers
!                                  !       (kg per square metre per sec)
!                                  !       FQW(,1) is total water flux
!                                  !       from surface, 'E'.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! INOUT FTL(,K) contains net
!                                  !       turbulent sensible heat flux
!                                  !       into layer K from below; so
!                                  !       FTL(,1) is the surface
!                                  !       sensible heat, H.(W/m2)
,ftl_surft(land_pts,nsurft)                                                   &
                             ! INOUT Surface FTL for land tiles
,fqw_surft(land_pts,nsurft)                                                   &
                             ! INOUT Surface FQW for land tiles
,olr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! IN    TOA - surface upward LW on
!                                  !       last radiation timestep
!                                  ! OUT   Corrected TOA outward LW
,taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end)               &
                             ! INOUT W'ly component of surface
!                                  !       wind stress over land
!                                  !       (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
,taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end)                &
                             ! INOUT W'ly component of surface
!                                  !       wind stress over mean sea
!                                  !       (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
,tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)               &
                             ! INOUT S'ly component of land sfc
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX
,tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX

REAL, INTENT(INOUT) ::                                                        &
  TScrnDcl_SSI(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                            !    Decoupled screen-level temperature
                            !    over sea or sea-ice
REAL, INTENT(INOUT) ::                                                        &
  TScrnDcl_SURFT(land_pts,nsurft)
                            !    Decoupled screen-level temperature
                            !    over land tiles
REAL, INTENT(INOUT) ::                                                        &
  tStbTrans(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                            !    Time since the transition to stable
                            !    conditions
REAL, INTENT(INOUT) ::                                                        &
 taux_land_star(udims%i_start:udims%i_end, udims%j_start:udims%j_end)         &
,tauy_land_star(vdims%i_start:vdims%i_end, vdims%j_start:vdims%j_end)         &
,taux_ssi_star( udims%i_start:udims%i_end, udims%j_start:udims%j_end)         &
,tauy_ssi_star( vdims%i_start:vdims%i_end, vdims%j_start:vdims%j_end)
  
!--------------------------------------------------------------------
!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!--------------------------------------------------------------------
REAL, INTENT(OUT) ::                                                          &
 ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! OUT Gridbox mean evaporation from
!                                  !     canopy/surface store (kg/m2/s).
!                                  !     Zero over sea.
,esoil_surft(land_pts,nsurft)                                                 &
                             ! OUT ESOIL for snow-free land tiles
,sea_ice_htf(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice)                                  &
                             ! OUT Heat flux through sea-ice
!                                  !     (W/m2, positive downwards).
!                                  !     (Not used for JULES)
,surf_ht_flux(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end)                                      &
!                                  ! OUT Net downward heat flux at
!                                  !     surface over land and sea-ice
!                                  !     fraction of gridbox (W/m2).
,surf_ht_flux_land(tdims%i_start:tdims%i_end,                                 &
                   tdims%j_start:tdims%j_end)                                 &
!                                  ! OUT Net downward heat flux at
!                                  !     surface over land
!                                  !     fraction of gridbox (W/m2).
,surf_ht_flux_sice(tdims%i_start:tdims%i_end,                                 &
                   tdims%j_start:tdims%j_end,nice)                            &
!                                  ! OUT Net category downward heat flux at
!                                  !     surface over sea-ice
!                                  !     fraction of gridbox (W/m2).
,surf_htf_surft(land_pts,nsurft)
!                                  ! OUT Net downward surface heat flux
!                                  !     on tiles (W/m2)

!-2 Genuinely output, needed by other atmospheric routines :-

REAL, INTENT(OUT) ::                                                          &
 tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! OUT   GBM surface temperature (K).
,tstar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! OUT   Land mean sfc temperature (K)
,tstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! OUT Ice area mean sea ice surface temperature
,le_surft(land_pts,nsurft)                                                    &
                             ! OUT Surface latent heat flux for
!                                  !     land tiles
,radnet_surft(land_pts,nsurft)                                                &
                             ! OUT Surface net radiation on
!                                  !       land tiles (W/m2)
,e_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! OUT Evaporation from sea times
!                                  !       leads fraction. Zero over
!                                  !       land. (kg per square metre
!                                  !       per sec).
,h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! OUT Surface sensible heat flux
!                                  !       over sea times leads fraction
!                                  !       (W/m2)
,taux_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)                  &
                             ! OUT   W'ly component of surface
!                                  !       wind stress (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
,tauy_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)                  &
                             ! OUT   S'ly component of surface
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX
,ei(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                             ! OUT Sublimation from lying snow or
!                                  !     sea-ice (kg/m2/s).
,ei_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)        &
                             ! OUT Sublimation from sea-ice
!                                  !     (kg/m2/s).
,ei_surft(land_pts,nsurft)                                                    &
                             ! OUT EI for land tiles.
,ecan_surft(land_pts,nsurft)                                                  &
                             ! OUT ECAN for snow-free land tiles
,esoil_soilt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nsoilt)      &
                             ! OUT Surface evapotranspiration
!                                  !     from soil moisture store
!                                  !     (kg/m2/s).
,ext_soilt(land_pts,nsoilt,sm_levels)                                         &
                             ! OUT Extraction of water from each
!                                  !     soil layer (kg/m2/s).
,snowmelt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                             ! OUT Snowmelt (kg/m2/s).
,melt_surft(land_pts,nsurft)                                                  &
                             ! OUT Snowmelt on land tiles (kg/m2/s
,rhokh_mix(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! OUT Exchange coeffs for moisture.
                             !     (Not used for JULES)

INTEGER, INTENT(OUT) ::                                                       &
 error          ! OUT 0 - AOK;
!                     !     1 to 7  - bad grid definition detected;
!--------------------------------------------------------------------
!  Workspace :-
!--------------------------------------------------------------------
REAL ::                                                                       &
 elake_surft(land_pts,nsurft)                                                 &
                             ! Lake evaporation.
,melt_ice_surft(land_pts,nsurft)                                              &
                             ! Ice melt on FLake lake tile (kg/m2/s)
,lake_ice_mass(land_pts)                                                      &
                             ! areal density equivalent to
                             ! lake ice of a given depth (kg/m2)
,qim_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! Implicit value of first model level
!                                  ! humidity
,tim_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! Implicit value of first model level
!                                  ! temperature
,tstar_rad4(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! Effective surface radiative
!                                  ! temperature for land and sea-ice
,tstar_surft_old(land_pts,nsurft)                                             &
!                                  ! Tile surface temperatures at
!                                  ! beginning of timestep.
,tstar_sice_cat_old(tdims%i_start:tdims%i_end,                                &
           tdims%j_start:tdims%j_end,nice_use)                                &
                             ! Sea ice surface T at beginning of timestep
,tstar_sea_old(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end)                                     &
                             ! Sea surface T at beginning of timestep
,sice_melt(tdims%i_start:tdims%i_end,                                         &
           tdims%j_start:tdims%j_end,nice)                                    &
                             !Melt at surface sea-ice category
,non_lake_frac(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end)                                     &
                             ! total tile fraction for surface types
                             ! other than inland water

,dftl_sice_ncat(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end)                                    &
                             ! Increment for ftl_ice from sea-ice
                             ! melt calculated for each category,
                             ! but un-weighted by category fractions
,dfqw_sice_ncat(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end)                                    &
                             ! Increment for fqw_ice from sea-ice
                             ! melt calculated for each category,
                             ! but un-weighted by category fractions
,dei_sice_ncat(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end)                                     &
                             ! Increment for ei_sice from sea-ice
                             ! melt calculated for each category,
                             ! but un-weighted by category fractions
,ice_fract_cat_use(tdims%i_start:tdims%i_end,                                 &
                   tdims%j_start:tdims%j_end,nice_use)                        &
                             ! Sea ice category fractions
                             ! If nice_use=1, this is the total ice
                             ! fraction
,ei_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! OUT Sublimation from lying snow
!                                  !     (kg/m2/s).

REAL, ALLOCATABLE :: tstar_ssi_old(:,:)
                                   ! Sea and sea-ice surface temperature
                                   ! at beginning of timestep --
                                   ! Required only for decoupled diagnosis,
                                   ! so allocatable, and local since it is
                                   ! used only on the predictor step

REAL, ALLOCATABLE :: tstar_sic(:,:,:)
                             !Ice category surface temperature
                             ! Only used if nice_use EQ 1

! dummy arrays required for sea and se-ice to create universal
! routines for all surfaces
REAL ::                                                                       &
 array_one(t_i_length * t_j_length)                                           &
                             ! Array of ones
,array_one_e_six(t_i_length * t_j_length)
                             ! Array of 1.0e6

REAL ::                                                                       &
 surf_ht_flux_sice_sm(tdims%i_start:tdims%i_end,                              &
                      tdims%j_start:tdims%j_end)                              &
                             ! Sea area mean seaice surface heat flux
,ei_sice_sm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Sea area mean sea ice sublimation

REAL ::                                                                       &
 canhc_surf(land_pts)
                             ! Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).

!  Local scalars :-

INTEGER ::                                                                    &
 i,j                                                                          &
            ! LOCAL Loop counter (horizontal field index).
,k                                                                            &
            ! LOCAL Tile pointer
,l                                                                            &
            ! LOCAL Land pointer
,n                                                                            &
            ! LOCAL Loop counter (tile index).
,m
            ! Loop counter for soil tiles

LOGICAL ::                                                                    &
 l_sice_new_code   ! Controls the sea ice temperature calculation
                   ! (See comments in sea ice section below)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!CABLE_LSM: see ACCESS-CM2 version for extra vars here
integer, parameter :: mm = 1
CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_IMPL2_cbl'

TYPE(air_type),       INTENT(inout)  :: air_cbl
TYPE(met_type),       INTENT(inout)  :: met_cbl
TYPE(radiation_type),       INTENT(inout)  :: rad_cbl
TYPE(roughness_type),     INTENT(inout)  :: rough_cbl
TYPE(canopy_type),    INTENT(inout)  :: canopy_cbl
TYPE(soil_snow_type),     INTENT(inout)  :: ssnow_cbl
TYPE(bgc_pool_type),       INTENT(inout)  :: bgc_cbl
TYPE(balances_type),       INTENT(inout)  :: bal_cbl
TYPE(sum_flux_type),  INTENT(inout)  :: sum_flux_cbl
TYPE(veg_parameter_type),   INTENT(inout) :: veg_cbl
TYPE(soil_parameter_type),  INTENT(inout) ::  soil_cbl

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error = 0

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,n,j,i)                                                        &
!$OMP SHARED(t_i_length,t_j_length,array_one,array_one_e_six,nice_use,        &
!$OMP        tdims,ice_fract_cat_use,ice_fract_ncat,ice_fract)

!$OMP DO SCHEDULE(STATIC)
DO l = 1, t_i_length * t_j_length
  array_one(l)       = 1.0
  array_one_e_six(l) = 1.0e6
END DO
!$OMP END DO NOWAIT

! Set up sea ice field depending on nice_use
IF (nice_use > 1) THEN
  ! Use all categories fully in surface exchange
  DO n = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        ice_fract_cat_use(i,j,n) = ice_fract_ncat(i,j,n)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
ELSE  ! nice_use=1
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      ice_fract_cat_use(i,j,1) = ice_fract(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL
!CABLE_LSM: jhan!!!
CALL im_sf_pt2_cbl (                                                              &
!CALL im_sf_pt2(                                                              &
 land_pts,land_index,nsurft,surft_index,surft_pts                             &
,flandg,tile_frac,snow_surft,nice_use,ice_fract,ice_fract_cat_use             &
,r_gamma,gamma1,gamma2,alpha1,alpha1_sea,alpha1_sice                          &
,ashtf_prime,ashtf_prime_sea,ashtf_prime_surft                                &
,resft,dtstar_surft,dtstar_sea,dtstar_sice                                    &
,rhokm_u_1,rhokm_v_1,rhokh_surft,rhokh_sice,rhokh_sea                         &
,ctctq1,dqw1_1,dtl1_1                                                         &
,cq_cm_u_1,cq_cm_v_1,du_1,dv_1,du_star1,dv_star1                              &
,flandg_u,flandg_v                                                            &
,fqw_1,ftl_1                                                                  &
,taux_1,taux_land,taux_land_star,taux_ssi,taux_ssi_star,tauy_1                &
,tauy_land,tauy_land_star,tauy_ssi,tauy_ssi_star                              &
,fqw_surft,epot_surft,ftl_surft,fqw_ice,ftl_ice,e_sea,h_sea                   &
,l_correct                                                                    &
)


!-----------------------------------------------------------------------

! Calculate surface scalar fluxes, temperatures only at the 1st call
! of the subroutine (first stage of the new BL solver) using standard
! MOSES2 physics and equations. These are the final values for this
! timestep and there is no need to repeat the calculation.
!-----------------------------------------------------------------------

tstar_surft_old  = 0.0
canhc_surf = 0.0
IF ( .NOT. l_correct ) THEN

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,n,j,i,k)                                                      &
!$OMP SHARED(tdims,ftl_1,h_sea,nice_use,ftl_ice,nsurft,surft_pts,surft_index, &
!$OMP        ftl_surft, radnet_surft, le_surft, nsoilt, land_pts,             &
!$OMP        t_soil_soilt, MELT_surft, soil_cbl, veg_cbl, sum_flux_cbl,       &
!$OMP        bal_cbl, bgc_cbl, ssnow_cbl, rough_cbl, rad_cbl, met_cbl,        &
!$OMP        air_cbl, sf_diag, ei_surft, esoil_surft, ecan_surft,             &
!$OMP        surf_htf_surft, surf_ht_flux_land, fqw_surft, fqw_1, ctctq1,     &
!$OMP        dqw1_1, dtl1_1, qw_1, tl_1, fland, sm_levels, canopy_cbl,        &
!$OMP        tstar_surft_old,tstar_surft,dtstar_surft,cp,error)

  !-----------------------------------------------------------------------
  ! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
  !-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      ftl_1(i,j) = ftl_1(i,j) * cp
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      h_sea(i,j) = cp * h_sea(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT

  DO n = 1,nice_use
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        ftl_ice(i,j,n) = cp * ftl_ice(i,j,n)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      ftl_surft(l,n) = cp * ftl_surft(l,n)
    END DO
!$OMP END DO NOWAIT
  END DO
!CABLE_LSM:_CM2 inits to avoid packing problems. copied from below location 
  ! initialise diagnostics to 0 to avoid packing problems
  DO n = 1, nsurft
!$OMP DO SCHEDULE(STATIC)
    DO l = 1, land_pts
      radnet_surft(l,n) = 0.0
      le_surft(l,n) = 0.0
    END DO
!$OMP END DO NOWAIT
  END DO

  !-----------------------------------------------------------------------
  ! Land surface calculations
  !-----------------------------------------------------------------------

 !CABLE_LSM:
  DO N=1,Nsurft
    DO L=1,LAND_PTS
      MELT_surft(L,N) = 0.
    ENDDO
  ENDDO

 call cable_implicit_main( & 
        1, &!C!cycleno, (in standalone control.F() declares as param=1) &
        !row_length, rows, dimensions & land fraction
        tdims%i_end, tdims%j_end, land_pts, nsurft, sm_levels, Fland,          &
        !forcing: temp and humidity & Jan 2018: Ticket #132 needs ctctq1
        tl_1, qw_1, dtl1_1,  dqw1_1, ctctq1,                                   & 
        !returned fluxes etc
        ftl_1, ftl_surft, fqw_1, fqw_surft,  tstar_surft,                      &
        surf_ht_flux_land, surf_htf_surft, ecan_surft, esoil_surft,            &
        ei_surft, radnet_surft,  sf_diag%t1p5m_surft, sf_diag% q1p5m_surft,    &
        melt_surft, dtstar_surft,                                              &  
        ! CABLE state vars  
        air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl, ssnow_cbl, bgc_cbl,  &
        bal_cbl, sum_flux_cbl, veg_cbl, soil_cbl )
!CABLE_LSM:End

  !-----------------------------------------------------------------------
  ! Optional error check : test for negative top soil layer temperature
  !-----------------------------------------------------------------------
  IF (l_neg_tstar) THEN
    DO m = 1,nsoilt
!$OMP DO SCHEDULE(STATIC)
      DO l = 1,land_pts
        IF (t_soil_soilt(l,m,1) < 0) THEN
          error = 1
          WRITE(jules_message,*) '*** ERROR DETECTED BY ROUTINE SF_IMPL2 ***'
          CALL jules_print('sf_impl2_jls',jules_message)
          WRITE(jules_message,*) 'NEGATIVE TEMPERATURE IN TOP SOIL LAYER AT '
          CALL jules_print('sf_impl2_jls',jules_message)
          WRITE(jules_message,*) 'LAND POINT ',l
          CALL jules_print('sf_impl2_jls',jules_message)
        END IF
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF

  !-----------------------------------------------------------------------
  !   Diagnose the land surface temperature
  !-----------------------------------------------------------------------

  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      tstar_surft_old(l,n) = tstar_surft(l,n)
!CABLE_LSM: tstar_surft(l,n) = tstar_surft_old(l,n) + dtstar_surft(l,n)
    END DO
!$OMP END DO NOWAIT
  END DO

!$OMP END PARALLEL

  !-----------------------------------------------------------------------
  ! 7.  Surface evaporation components and updating of surface
  !     temperature (P245, routine SF_EVAP).
  !-----------------------------------------------------------------------
!CABLE_LSM: 
!C!  CALL sf_evap (                                                              &
!C!    land_pts,nsurft,                                                          &
!C!    land_index,surft_index,surft_pts,sm_levels,fland,                         &
!C!    ashtf_prime_surft,canopy,dtrdz_charney_grid_1,flake,fraca,                &
!C!    snow_surft,resfs,resft,rhokh_surft,tile_frac,smc_soilt,wt_ext_surft,      &
!C!    timestep,r_gamma,fqw_1,fqw_surft,ftl_1,ftl_surft,tstar_surft,             &
!C!    ecan,ecan_surft,elake_surft,esoil_soilt,esoil_surft,ei_surft,ext_soilt,   &
!C!    sf_diag)

!CABLE_LSM:-CM2 hassection in here for elake, ecan and esoil - possibly missed
!CABLE_LSM: sf_evap stuff 
    DO N=1,Nsurft
      DO L=1,LAND_PTS
        ELAKE_surft(L,N) = 0.
      ENDDO
    ENDDO

    DO j=tdims%j_start,tdims%j_end
     DO i=tdims%i_start,tdims%i_end
      ECAN(I,J) = 0.
      ESOIL_soilt(I,J,:) = 0.
     ENDDO
    ENDDO
!write(6,*) "WARNING: Hard-wiring of mm forces ssumption of 1 soil_type. sf_impl2_cbl"

    DO N=1,Nsurft
     DO K=1,surft_PTS(N)
       L = surft_INDEX(K,N)
       j=(land_index(l)-1)/tdims%i_end + 1
       i = land_index(l) - (j-1)*tdims%i_end
       ECAN(I,J) = ECAN(I,J) + TILE_FRAC(L,N)*ECAN_surft(L,N)
       ESOIL_soilt(I,J,mm) = ESOIL_soilt(I,J,mm) + TILE_FRAC(L,N)*ESOIL_surft(L,N)          
     ENDDO
    ENDDO   
    EXT_soilt = 0. ! MRD
!CABLE_LSM:End



  !-----------------------------------------------------------------------
  !     Surface melting of sea-ice and snow on land tiles.
  !-----------------------------------------------------------------------

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,n,j,i)                                                        &
!$OMP SHARED(tdims,ei_land,snowmelt,nsurft,land_pts,melt_ice_surft)

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      ei_land(i,j)  = 0.0
      snowmelt(i,j) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

  ! Lake initialisation
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,land_pts
      melt_ice_surft(l,n) = 0.0
    END DO
!$OMP END DO NOWAIT
  END DO

!$OMP END PARALLEL

!CABLE_LSM: -CM2 dodges 
  DO n = 1,nsurft
!C!    CALL sf_melt (                                                            &
!C!      land_pts,land_index,                                                    &
!C!      surft_index(:,n),surft_pts(n),flandg,                                   &
!C!      alpha1(:,n),ashtf_prime_surft(:,n),dtrdz_charney_grid_1,                &
!C!      resft(:,n),rhokh_surft(:,n),tile_frac(:,n),timestep,r_gamma,            &
!C!      ei_surft(:,n),fqw_1,ftl_1,fqw_surft(:,n),ftl_surft(:,n),                &
!C!      tstar_surft(:,n),snow_surft(:,n),snowdep_surft(:,n),                    &
!C!      melt_surft(:,n)                                                         &
!C!      )
!C!
!C!    !-----------------------------------------------------------------------
!C!    ! thermodynamic, flux contribution of melting ice on the FLake lake tile
!C!    !-----------------------------------------------------------------------
!C!    IF (     (l_flake_model   )                                               &
!C!        .AND. ( .NOT. l_aggregate)                                            &
!C!        .AND. (n == lake       ) ) THEN
!C!
!C!      ! lake_h_ice_gb is only initialised if FLake is on.
!C!
!C!!$OMP PARALLEL DO                                                             &
!C!!$OMP SCHEDULE(STATIC)                                                        &
!C!!$OMP DEFAULT(NONE)                                                           &
!C!!$OMP PRIVATE(l)                                                              &
!C!!$OMP SHARED(land_pts,lake_ice_mass,lake_h_ice_gb)
!C!      DO l = 1, land_pts
!C!        lake_ice_mass(l) = lake_h_ice_gb(l) * rho_ice
!C!      END DO
!C!!$OMP END PARALLEL DO
!C!
!C!      CALL sf_melt (                                                          &
!C!        land_pts,land_index,                                                  &
!C!        surft_index(:,n),surft_pts(n),flandg,                                 &
!C!        alpha1(:,n),ashtf_prime_surft(:,n),dtrdz_charney_grid_1,              &
!C!        resft(:,n),rhokh_surft(:,n),tile_frac(:,n),timestep,r_gamma,          &
!C!        ei_surft(:,n),fqw_1,ftl_1,fqw_surft(:,n),ftl_surft(:,n),              &
!C!        tstar_surft(:,n),lake_ice_mass,lake_ice_mass / rho_snow_const,        &
!C!        melt_ice_surft(:,n)                                                   &
!C!          )
!C!    END IF

    !-----------------------------------------------------------------------
    !  Increment snow by sublimation and melt
    !-----------------------------------------------------------------------

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,l,j,i)                                                        &
!$OMP SHARED(surft_pts,surft_index,land_index,t_i_length,ei_land,tile_frac,   &
!$OMP        ei_surft,snowmelt,melt_surft,n)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      ei_land(i,j) = ei_land(i,j) + tile_frac(l,n) * ei_surft(l,n)
      snowmelt(i,j) = snowmelt(i,j) +                                         &
                      tile_frac(l,n) * melt_surft(l,n)
    END DO
!$OMP END PARALLEL DO

  END DO

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(l,n,j,i,k)

  IF (sf_diag%smlt) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        sf_diag%snomlt_surf_htf(i,j) = lf * snowmelt(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      surf_ht_flux_land(i,j) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

!CABLE_LSM:ACCESS-1.3 had stuff here - see comment in -CM2 though 
  IF (     (l_flake_model   )                                                 &
      .AND. ( .NOT. l_aggregate) ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        surf_ht_flux_lake_ij(i,j) = 0.0
        ! initialise the non-lake fraction to one, not zero,
        ! in case there should ever be more than one lake tile, see below
        non_lake_frac(    i,j) = 1.0
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    tstar_land(i,j) = 0.0
  END DO
!$OMP END DO NOWAIT

!CABLE_LSM: 
!CABLE_LSM:this section copied to above location 
!C!  ! initialise diagnostics to 0 to avoid packing problems
!C!  DO n = 1, nsurft
!C!!$OMP DO SCHEDULE(STATIC)
!C!    DO l = 1, land_pts
!C!      radnet_surft(l,n) = 0.0
!C!      le_surft(l,n) = 0.0
!C!    END DO
!C!!$OMP END DO NOWAIT
!C!  END DO

  IF (sf_diag%l_lw_surft) THEN
    DO n = 1, nsurft
!$OMP DO SCHEDULE(STATIC)
      DO l = 1, land_pts
        sf_diag%lw_up_surft(l,n) = 0.0
        sf_diag%lw_down_surft(l,n) = 0.0
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF

!$OMP BARRIER

  IF (l_skyview) THEN
!CABLE_LSM: 
!C!    DO n = 1,nsurft
!C!!$OMP DO SCHEDULE(STATIC)
!C!      DO k = 1,surft_pts(n)
!C!        l = surft_index(k,n)
!C!        j=(land_index(l) - 1) / tdims%i_end + 1
!C!        i = land_index(l) - (j-1) * tdims%i_end
!C!        radnet_surft(l,n) = sw_surft(l,n) +   emis_surft(l,n) *               &
!C!          sky(i,j) * ( lw_down(i,j) + lw_down_elevcorr_surft(l,n)             &
!C!                                  - sbcon * tstar_surft(l,n)**4 )
!C!      END DO
!C!!$OMP END DO
!C!    END DO
    IF (sf_diag%l_lw_surft) THEN
      DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
        DO k = 1,surft_pts(n)
          l = surft_index(k,n)
          j=(land_index(l) - 1) / tdims%i_end + 1
          i = land_index(l) - (j-1) * tdims%i_end
          sf_diag%lw_up_surft(l,n)   = emis_surft(l,n) * sky(i,j) *           &
                                       sbcon * tstar_surft(l,n)**4
          sf_diag%lw_down_surft(l,n) = emis_surft(l,n) * sky(i,j) *           &
                                       (lw_down(i,j) +                        &
                                       lw_down_elevcorr_surft(l,n))
        END DO
!$OMP END DO
      END DO
    END IF
  ELSE
!CABLE_LSM: 
!C!    DO n = 1,nsurft
!C!!$OMP DO SCHEDULE(STATIC)
!C!      DO k = 1,surft_pts(n)
!C!        l = surft_index(k,n)
!C!        j=(land_index(l) - 1) / tdims%i_end + 1
!C!        i = land_index(l) - (j-1) * tdims%i_end
!C!        radnet_surft(l,n) = sw_surft(l,n) +   emis_surft(l,n) *               &
!C!                   ( lw_down(i,j) + lw_down_elevcorr_surft(l,n)               &
!C!                                  - sbcon * tstar_surft(l,n)**4 )
!C!      END DO
!C!!$OMP END DO
!C!    END DO
    IF (sf_diag%l_lw_surft) THEN
      DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
        DO k = 1,surft_pts(n)
          l = surft_index(k,n)
          j=(land_index(l) - 1) / tdims%i_end + 1
          i = land_index(l) - (j-1) * tdims%i_end
          sf_diag%lw_up_surft(l,n)   = emis_surft(l,n) * sbcon *              &
                                       tstar_surft(l,n)**4
          sf_diag%lw_down_surft(l,n) = emis_surft(l,n) * (lw_down(i,j) +      &
                                       lw_down_elevcorr_surft(l,n))
        END DO
!$OMP END DO
      END DO
    END IF
  END IF

  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      canhc_surf(l) = canhc_surft(l,n)
      IF ( ( .NOT. cansnowtile(n)) .AND. l_snow_nocan_hc .AND.                &
           (nsmax > 0) .AND. (nsnow_surft(l,n) > 0) ) canhc_surf(l) = 0.0
      le_surft(l,n) = lc * ecan_surft(l,n) + lc * esoil_surft(l,n) +          &
                     lc * elake_surft(l,n) + ls * ei_surft(l,n)
      surf_ht_store_surft(l,n) = (canhc_surf(l) / timestep) *                 &
                           (tstar_surft(l,n) - tstar_surft_old(l,n))
!CABLE_LSM: from -CM2 - replace w CABLE field 
!C!      surf_htf_surft(l,n) = radnet_surft(l,n) + anthrop_heat_surft(l,n) -     &
!C!                          ftl_surft(l,n) -                                    &
!C!                          le_surft(l,n) -                                     &
!C!                          lf * (melt_surft(l,n) + melt_ice_surft(l,n)) -      &
!C!                          surf_ht_store_surft(l,n)
      ! separate out the lake heat flux for FLake
      ! and replace the snow-melt (NSMAX=0 only) and ice-melt heat fluxes
      ! so Flake can do its melting
      IF (     (l_flake_model   )                                             &
          .AND. ( .NOT. l_aggregate)                                          &
          .AND. (n == lake       ) ) THEN
        IF (nsmax == 0) THEN
          surf_ht_flux_lake_ij(i,j) = surf_htf_surft(l,n)                     &
                        + lf * (melt_surft(l,n) + melt_ice_surft(l,n))
          non_lake_frac(    i,j) = non_lake_frac(i,j) - tile_frac(l,n)
        ELSE
          surf_ht_flux_lake_ij(i,j) = surf_htf_surft(l,n)                     &
                        + lf * melt_ice_surft(l,n)
          non_lake_frac(    i,j) = non_lake_frac(i,j) - tile_frac(l,n)
        END IF
      ELSE
        surf_ht_flux_land(i,j) = surf_ht_flux_land(i,j)                       &
                          + tile_frac(l,n) * surf_htf_surft(l,n)
      END IF
      tstar_land(i,j) = tstar_land(i,j)                                       &
                 + tile_frac(l,n) * tstar_surft(l,n)
    END DO
!$OMP END DO
  END DO

!CABLE_LSM:from -CM2
!C!  ! normalise the non-lake surface heat flux
!C!  IF (     (l_flake_model   )                                                 &
!C!      .AND. ( .NOT. l_aggregate) ) THEN
!C!!$OMP DO SCHEDULE(STATIC)
!C!    DO j = tdims%j_start,tdims%j_end
!C!      DO i = tdims%i_start,tdims%i_end
!C!        ! be careful about gridboxes that are all lake
!C!        IF (non_lake_frac(i,j) > EPSILON(0.0)) THEN
!C!          surf_ht_flux_land(i,j) =   surf_ht_flux_land(i,j)                   &
!C!                                   / non_lake_frac(i,j)
!C!        END IF
!C!      END DO
!C!    END DO
!C!!$OMP END DO
!C!  END IF


  !-----------------------------------------------------------------------
  ! Optional error check : test for negative surface temperature
  !-----------------------------------------------------------------------
  IF (l_neg_tstar) THEN
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,land_pts
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      IF (tstar_land(i,j) < 0) THEN
        error = 1
        WRITE(jules_message,*) '*** ERROR DETECTED BY ROUTINE SF_IMPL2 ***'
        CALL jules_print('sf_impl2_jls',jules_message)
        WRITE(jules_message,*) 'NEGATIVE SURFACE TEMPERATURE AT LAND POINT ',l
        CALL jules_print('sf_impl2_jls',jules_message)
      END IF
    END DO
!$OMP END DO
  END IF

!$OMP END PARALLEL

  !-----------------------------------------------------------------------
  ! Sea and sea-ice surface calculations
  !-----------------------------------------------------------------------
  ! Allow tstar_sea to be modified by dtstar_sea
  IF (l_use_dtstar_sea) THEN
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        IF ( flandg(i,j) < 1.0 ) THEN
          tstar_sea_old(i,j) = tstar_sea(i,j)
          tstar_sea(i,j) = tstar_sea_old(i,j) + dtstar_sea(i,j)
        END IF
      END DO
    END DO
  END IF

  ! Set control logical (l_sice_new_code) to determine how to calculate the
  ! sea ice surface temperature.  l_sice_new_code = T is the method that is
  ! compatible with using the sea ice categories fully in the radiation and
  ! surface exchange code (so nice_use=nice).  l_sice_new_code = F is the
  ! old method that only uses the categories in the implicit surface exchange
  ! code (so nice_use = 1).

  ! NOTE that nice_use=nice=1 can use either scheme and the choice is
  ! determined by the logical l_tstar_sice_new set in switches.F90.  The
  ! difference between the 2 methods is small but causes differences in the
  ! results that are larger than bit level, hence the need to control which
  ! method is used.

  IF (nice_use == 1 .AND. nice == 1) THEN
    ! If nice_use=nice=1 : Choice method determined by logical
    ! l_tstar_sice_new
    l_sice_new_code = l_tstar_sice_new

  ELSE IF (nice_use /= nice) THEN
    ! Old calculation must be used
    l_sice_new_code = .FALSE.

  ELSE IF (nice > 1) THEN
    ! New calculation must be used
    l_sice_new_code = .TRUE.
  END IF

  IF ( .NOT. l_sice_new_code) THEN

    ALLOCATE(tstar_sic(tdims%i_start:tdims%i_end,                             &
                       tdims%j_start:tdims%j_end,nice))

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(n,j,i)                                                          &
!$OMP SHARED(tstar_sic,nice,tdims)
    DO n = 1, nice
!$OMP DO SCHEDULE(STATIC)
      DO j = tdims%j_start,tdims%j_end
        DO i= tdims%i_start,tdims%i_end
          tstar_sic(i,j,n) = 0.0
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
!$OMP END PARALLEL
  END IF

  !-----------------------------------------------------------------------
  ! Store old surface temperature for sea and sea-ice if using the
  ! decoupled diagnostic.
  !-----------------------------------------------------------------------
  IF ((IScrnTDiag == IP_ScrnDecpl2) .OR. (IScrnTDiag == IP_ScrnDecpl3)) THEN
    ALLOCATE(tstar_ssi_old(tdims%i_start:tdims%i_end,                         &
                           tdims%j_start:tdims%j_end))

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,tstar_ssi_old,tstar_ssi)
    DO j = tdims%j_start,tdims%j_end
      DO i= tdims%i_start,tdims%i_end
        tstar_ssi_old(i,j) = tstar_ssi(i,j)
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE
    ALLOCATE(tstar_ssi_old(1,1))
  END IF

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(n,j,i)                                                          &
!$OMP SHARED(surf_ht_flux_sice_sm,sice_melt,sea_ice_htf,ei_sice,fqw_ice,      &
!$OMP        tdims,nice,nice_use,sf_diag)

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i= tdims%i_start,tdims%i_end
      surf_ht_flux_sice_sm(i,j) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

  DO n = 1, nice
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i= tdims%i_start,tdims%i_end
        sice_melt(i,j,n)          = 0.0
        sea_ice_htf(i,j,n)        = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

  DO n = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i= tdims%i_start,tdims%i_end
        ei_sice(i,j,n)            = fqw_ice(i,j,n)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

  IF (sf_diag%simlt) THEN
    DO n = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
      DO j = tdims%j_start,tdims%j_end
        DO i= tdims%i_start,tdims%i_end
          sf_diag%sice_mlt_htf(i,j,n) = 0.0
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF

!$OMP END PARALLEL

  !-----------------------------------------------------------------------
  ! Diagnose the surface temperature for points with sea-ice
  ! Note that k_sice = 2.0*thermal conductivity/surface layer thickness
  !-----------------------------------------------------------------------

  IF (l_sice_new_code) THEN

    ! Update tstar_sice_cat using dtstar_sice

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(n,j,i)                                                          &
!$OMP SHARED(nice_use,tdims,flandg,ice_fract_cat_use,tstar_sice_cat_old,      &
!$OMP        tstar_sice_cat,dtstar_sice)
    DO n = 1,nice_use
!$OMP DO SCHEDULE(STATIC)
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          IF ( flandg(i,j) < 1.0 .AND. ice_fract_cat_use(i,j,n) > 0 ) THEN
            tstar_sice_cat_old(i,j,n) = tstar_sice_cat(i,j,n)
            tstar_sice_cat(i,j,n) = tstar_sice_cat_old(i,j,n) +               &
                                    dtstar_sice(i,j,n)
          END IF
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
!$OMP END PARALLEL

  ELSE  ! Use old code

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,j,i,n)                                                        &
!$OMP SHARED(tdims,flandg,ice_fract,surf_ht_flux_sice_sm,radnet_sice,         &
!$OMP        tstar_sice_cat,dtstar_sice,ftl_ice,fqw_ice,ice_fract_ncat,       &
!$OMP        tstar_sic,ti,k_sice,emis_sice,nice)
    ! Update tstar_sic using the surface heat equation
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        IF ( flandg(i,j) < 1.0 .AND. ice_fract(i,j) > 0.0 ) THEN
          surf_ht_flux_sice_sm(i,j) = radnet_sice(i,j,1) -                    &
             4.0 * emis_sice * sbcon * (tstar_sice_cat(i,j,1)**3.0) *         &
                 dtstar_sice(i,j,1) - ftl_ice(i,j,1) - ls * fqw_ice(i,j,1)
          DO n = 1,nice
            IF ( ice_fract_ncat(i,j,n) > 0.0 ) THEN
              tstar_sic(i,j,n) = ti(i,j,n) +                                  &
                       surf_ht_flux_sice_sm(i,j) / k_sice(i,j,n)
            END IF
          END DO
        END IF
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF (l_sice_new_code) THEN

    DO n = 1, nice_use
      CALL sf_melt (                                                          &
        ssi_pts,ssi_index,                                                    &
        sice_index_ncat(:,n),sice_pts_ncat(n),fssi_ij,                        &
        alpha1_sice(:,:,n),ashtf_prime(:,:,n),dtrdz_charney_grid_1,           &
        array_one,rhokh_sice(:,:,n),sice_frac_ncat(:,n),timestep,             &
        r_gamma,                                                              &
        ei_sice(:,:,n),fqw_1,ftl_1,fqw_ice(:,:,n),ftl_ice(:,:,n),             &
        tstar_sice_cat(:,:,n),array_one_e_six,                                &
        array_one_e_six / rho_snow_const,                                     &
        sice_melt(:,:,n)                                                      &
        )

      IF (sf_diag%simlt) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,l)                                                        &
!$OMP SHARED(sice_pts_ncat,n,sice_index_ncat,ssi_index,t_i_length,sf_diag,    &
!$OMP        sice_melt)
        DO k = 1,sice_pts_ncat(n)
          l = sice_index_ncat(k,n)
          j=(ssi_index(l) - 1) / t_i_length + 1
          i = ssi_index(l) - (j-1) * t_i_length
          sf_diag%sice_mlt_htf(i,j,n) = lf * sice_melt(i,j,n)
        END DO
!$OMP END PARALLEL DO
      END IF

    END DO

  ELSE ! Use old code

    DO n = 1,nice

      ! Since sea-ice categories are not actually tiled for their surface
      ! fluxes here, then the increment to ftl_ice, fqw_ice and
      ! ei_sice are not correctly weighted in sf_melt. Hence need to keep
      ! the increments and update ftl_ice, fqw_ice and ei_sice with
      ! weighted contributions below

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,dftl_sice_ncat,dfqw_sice_ncat,dei_sice_ncat)
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          dftl_sice_ncat(i,j) = 0.0
          dfqw_sice_ncat(i,j) = 0.0
          dei_sice_ncat(i,j)  = 0.0
        END DO
      END DO
!$OMP END PARALLEL DO

      CALL sf_melt (                                                          &
        ssi_pts,ssi_index,                                                    &
        sice_index_ncat(:,n),sice_pts_ncat(n),fssi_ij,                        &
        alpha1_sice(:,:,1),ashtf_prime(:,:,1),dtrdz_charney_grid_1,           &
        array_one,rhokh_sice(:,:,1),sice_frac_ncat(:,n),timestep,             &
        r_gamma,                                                              &
        dei_sice_ncat,fqw_1,ftl_1,dfqw_sice_ncat,dftl_sice_ncat,              &
        tstar_sic(:,:,n),array_one_e_six,                                     &
        array_one_e_six / rho_snow_const,                                     &
        sice_melt(:,:,n)                                                      &
        )

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,l,j,i)                                                        &
!$OMP SHARED(sice_pts_ncat,n,sice_index_ncat,ssi_index,t_i_length,sf_diag,    &
!$OMP        sice_melt,ftl_ice,sice_frac_ncat,sice_frac,dftl_sice_ncat,       &
!$OMP        fqw_ice,dfqw_sice_ncat,ei_sice,dei_sice_ncat)
      DO k = 1,sice_pts_ncat(n)
        l = sice_index_ncat(k,n)
        j=(ssi_index(l) - 1) / t_i_length + 1
        i = ssi_index(l) - (j-1) * t_i_length
        IF (sf_diag%simlt) sf_diag%sice_mlt_htf(i,j,n) = lf * sice_melt(i,j,n)
        ! Add weighted increments to ftl_ice, fqw_ice and ei_sice
        ftl_ice(i,j,1) = ftl_ice(i,j,1)                                       &
           +( sice_frac_ncat(l,n) / sice_frac(l) ) * dftl_sice_ncat(i,j)
        fqw_ice(i,j,1) = fqw_ice(i,j,1)                                       &
           +( sice_frac_ncat(l,n) / sice_frac(l) ) * dfqw_sice_ncat(i,j)
        ei_sice(i,j,1) = ei_sice(i,j,1)                                       &
           +( sice_frac_ncat(l,n) / sice_frac(l) ) * dei_sice_ncat(i,j)
      END DO
!$OMP END PARALLEL DO
    END DO

  END IF

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(i,j,k,l,n)

  !-----------------------------------------------------------------------
  !     Gridbox-mean surface temperature and net surface heat fluxes
  !-----------------------------------------------------------------------

  ! Assign SW flux on categories to radnet_sice to avoid cumbersome
  ! indexing in the following loops.

  DO n = 1,nice_use
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        radnet_sice(i,j,n) = 0.0
      END DO
    END DO
!$OMP END DO
  END DO

  IF (nice_use > 1) THEN
    !   In this case, nice_use = nice.
    !   Radiative fluxes are on all categories in use, so use the
    !   full arrays.
    DO n = 1,nice_use
!$OMP DO SCHEDULE(STATIC)
      DO k = 1,sice_pts_ncat(n)
        l = sice_index_ncat(k,n)
        j=(ssi_index(l) - 1) / t_i_length + 1
        i = ssi_index(l) - (j-1) * t_i_length
        radnet_sice(i,j,n) = sw_sicat(l,n)
      END DO
!$OMP END DO
    END DO
  ELSE
    !   In this case n_ice_use must be 1, so indexing is over all sea-ice
    !   points.
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,sice_pts
      l = sice_index(k)
      j=(ssi_index(l) - 1) / t_i_length + 1
      i = ssi_index(l) - (j-1) * t_i_length
      radnet_sice(i,j,1) = sw_sicat(l,1)
    END DO
!$OMP END DO
  END IF


  IF (l_sice_new_code) THEN

!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        surf_ht_flux_sice_sm(i,j) = 0.0
        DO n = 1,nice_use
          surf_ht_flux_sice(i,j,n)= 0.0
        END DO
        IF ( flandg(i,j) < 1.0 .AND. ice_fract(i,j) > 0.0 ) THEN
          tstar_sice(i,j)= 0.0
          DO n = 1,nice_use
            IF ( ice_fract_ncat(i,j,n) > 0.0 ) THEN
              tstar_sice(i,j) = tstar_sice(i,j) +                             &
                                (ice_fract_ncat(i,j,n) /                      &
                               ice_fract(i,j)) * tstar_sice_cat(i,j,n)
            ELSE
              tstar_sice_cat(i,j,n) = tfs   ! copy setting in sice_htf
            END IF
          END DO
          tstar_ssi(i,j) = (1.0 - ice_fract(i,j)) * tstar_sea(i,j) +          &
                            ice_fract(i,j) * tstar_sice(i,j)


          DO n = 1,nice_use
            IF ( ice_fract_cat_use(i,j,n) > 0.0 ) THEN

              radnet_sice(i,j,n) = radnet_sice(i,j,n) + emis_sice *           &
                ( lw_down(i,j) - sbcon * tstar_sice_cat(i,j,n)**4 )

              surf_ht_flux_sice(i,j,n) = radnet_sice(i,j,n) -                 &
                               ftl_ice(i,j,n) - ls * fqw_ice(i,j,n) -         &
                               lf * sice_melt(i,j,n)
              surf_ht_flux_sice_sm(i,j) = surf_ht_flux_sice_sm(i,j) +         &
                           (ice_fract_cat_use(i,j,n) / ice_fract(i,j))*       &
                                surf_ht_flux_sice(i,j,n)
            END IF
          END DO

        ELSE IF (flandg(i,j) < 1.0) THEN  ! non-icy ocean point
          tstar_sice_cat(i,j,:) = tfs
          tstar_ssi(i,j) = tstar_sea(i,j)
        END IF
        IF (nice == 1) THEN
          tstar_sice(i,j) = tstar_sice_cat(i,j,1)  ! Ensure these are the
                                                   ! same at all points
        END IF
      END DO
    END DO
!$OMP END DO NOWAIT

  ELSE   ! Use old code

!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        surf_ht_flux_sice_sm(i,j) = 0.0
        DO n = 1,nice
          surf_ht_flux_sice(i,j,n) = 0.0
        END DO
        IF ( flandg(i,j) < 1.0 .AND. ice_fract(i,j) > 0.0 ) THEN
          tstar_sice(i,j) = 0.0
          DO n = 1,nice
            IF ( ice_fract_ncat(i,j,n) > 0.0 ) THEN
              tstar_sice(i,j) = tstar_sice(i,j) +                             &
                                (ice_fract_ncat(i,j,n) /                      &
                                 ice_fract(i,j)) * tstar_sic(i,j,n)
            END IF
          END DO
          tstar_ssi(i,j) = (1.0 - ice_fract(i,j)) * tstar_sea(i,j) +          &
                            ice_fract(i,j) * tstar_sice(i,j)

          radnet_sice(i,j,1) = radnet_sice(i,j,1) + emis_sice *               &
            ( lw_down(i,j) - sbcon * tstar_sice(i,j)**4 )

          DO n = 1,nice
            IF ( ice_fract_ncat(i,j,n) > 0.0 ) THEN
              surf_ht_flux_sice(i,j,n) = radnet_sice(i,j,1) -                 &
                               4.0 * emis_sice * sbcon *                      &
                               tstar_sice(i,j)**3 *                           &
                               (tstar_sic(i,j,n) - tstar_sice(i,j)) -         &
                               ftl_ice(i,j,1) - ls * fqw_ice(i,j,1) -         &
                               lf * sice_melt(i,j,n)
              surf_ht_flux_sice_sm(i,j) = surf_ht_flux_sice_sm(i,j) +         &
                           (ice_fract_ncat(i,j,n) / ice_fract(i,j))*          &
                            surf_ht_flux_sice(i,j,n)

            END IF
          END DO

        ELSE IF (flandg(i,j) < 1.0) THEN  ! non-icy ocean point
          tstar_sice(i,j) = tfs
          tstar_ssi(i,j) = tstar_sea(i,j)
        END IF

        tstar_sice_cat(i,j,1) = tstar_sice(i,j)  ! Ensure these are the
                                                 ! same at all points
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  ! Convert sea and sea-ice fluxes to be fraction of grid-box
  ! (as required by sea and sea-ice modellers)

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      h_sea(i,j)=(1.0 - ice_fract(i,j)) * h_sea(i,j)
      e_sea(i,j)=(1.0 - ice_fract(i,j)) * e_sea(i,j)
      surf_ht_flux_sice_sm(i,j) = ice_fract(i,j) * surf_ht_flux_sice_sm(i,j)
      ei_sice_sm(i,j) = 0.0
      DO n = 1,nice_use
        ei_sice(i,j,n) = ice_fract_cat_use(i,j,n) * ei_sice(i,j,n)
        ei_sice_sm(i,j)= ei_sice_sm(i,j) + ei_sice(i,j,n)
        ftl_ice(i,j,n) = ice_fract_cat_use(i,j,n) * ftl_ice(i,j,n)
        fqw_ice(i,j,n) = ice_fract_cat_use(i,j,n) * fqw_ice(i,j,n)
        radnet_sice(i,j,n) = ice_fract_cat_use(i,j,n) * radnet_sice(i,j,n)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (sf_diag%l_ftl_ice_sm) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        sf_diag%ftl_ice_sm(i,j) = SUM(ftl_ice(i,j,:))
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (sf_diag%l_tstar_sice_weighted) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        sf_diag%tstar_sice_weighted(i,j) = ice_fract(i,j) * tstar_sice(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (sf_diag%l_tstar_sice_weighted_cat) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        DO n = 1,nice_use
          sf_diag%tstar_sice_weighted_cat(i,j,n) =                            &
             ice_fract_cat_use(i,j,n) * tstar_sice_cat(i,j,n)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  ! Compute sea ice time fraction variable, required for CMIP6.
  IF (sf_diag%l_ice_present_cat) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        DO n = 1,nice_use
          IF (ice_fract_cat_use(i,j,n) > 0.0) THEN
            sf_diag%ice_present_cat(i,j,n) = 1.0
          ELSE
            sf_diag%ice_present_cat(i,j,n) = 0.0
          END IF
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (sf_diag%l_ice_present) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        IF (ice_fract(i,j) > 0.0) THEN
          sf_diag%ice_present(i,j) = 1.0
        ELSE
          sf_diag%ice_present(i,j) = 0.0
        END IF
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  ! Calculate surface upward LW over sea ice, weighted by ice fraction.  This
  ! is required for CMIP6.
  IF (sf_diag%l_lw_up_sice_weighted_cat) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        DO n = 1,nice_use
          IF ((ice_fract_cat_use(i,j,n) > 0.0) .AND. (flandg(i,j) < 1.0)) THEN
            sf_diag%lw_up_sice_weighted_cat(i,j,n) =                          &
                                          ice_fract_cat_use(i,j,n)            &
                                          * (emis_sice                        &
                                          * sbcon * tstar_sice_cat(i,j,n)**4  &
                                          + (1.0 - emis_sice) * lw_down(i,j))
          ELSE
            sf_diag%lw_up_sice_weighted_cat(i,j,n) = 0.0
          END IF ! ice_fract_cat_use
        END DO   ! n
      END DO    ! i
    END DO    ! j
!$OMP END DO NOWAIT
  END IF

  IF (sf_diag%l_lw_up_sice_weighted) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        IF ((ice_fract(i,j) > 0.0) .AND. (flandg(i,j) < 1.0)) THEN
          sf_diag%lw_up_sice_weighted(i,j) = ice_fract(i,j)                   &
                                             * (emis_sice * sbcon             &
                                                * tstar_sice(i,j)**4          &
                                                + (1.0 - emis_sice)           &
                                                * lw_down(i,j))
        ELSE
          sf_diag%lw_up_sice_weighted(i,j) = 0.0
        END IF    ! flandg, ice_fract
      END DO    ! i
    END DO    ! j
!$OMP END DO NOWAIT
  END IF

  !-----------------------------------------------------------------------
  ! GBM diagnostic calculations
  !-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      qim_1(i,j) = qw_1(i,j) + dqw1_1(i,j) - ctctq1(i,j) * fqw_1(i,j)
      tim_1(i,j) = tl_1(i,j) + dtl1_1(i,j) - ctctq1(i,j) * ftl_1(i,j) / cp
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      tstar(i,j) = flandg(i,j) * tstar_land(i,j)                              &
                   + (1.0 - flandg(i,j)) * tstar_ssi(i,j)
      ei(i,j) = flandg(i,j) * ei_land(i,j)                                    &
                + (1.0 - flandg(i,j)) * ei_sice_sm(i,j)
      surf_ht_flux(i,j) = flandg(i,j) * surf_ht_flux_land(i,j)                &
                          + (1.0 - flandg(i,j)) * surf_ht_flux_sice_sm(i,j)
      rhokh_mix(i,j) = rhokh(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT

  ! TOA outward LW radiation after boundary layer

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      tstar_rad4(i,j) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      !     The contribution from the sea is removed in UM imp_solver to keep it
      !     consistent with the UM radiation scheme (see olr comment there)
      tstar_rad4(i,j) = tstar_rad4(i,j) + (1.0 - flandg(i,j))*                &
                        (1.0 - ice_fract(i,j)) * emis_sea *                   &
                        tstar_sea(i,j)**4
      DO n = 1,nice_use
        tstar_rad4(i,j) = tstar_rad4(i,j) + (1.0 - flandg(i,j))*              &
                   ice_fract_cat_use(i,j,n) * emis_sice *                     &
                   tstar_sice_cat(i,j,n)**4

      END DO
    END DO
  END DO
!$OMP END DO

  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      !     For historical compatibility, the addjustment of the OLR
      !     may be made with or without the surface emissivity.
      IF (l_dolr_land_black) THEN
        tstar_rad4(i,j) = tstar_rad4(i,j) + flandg(i,j) *                     &
                          tile_frac(l,n) * tstar_surft(l,n)**4
      ELSE
        tstar_rad4(i,j) = tstar_rad4(i,j) + flandg(i,j) *                     &
                          tile_frac(l,n) * emis_surft(l,n) *                  &
                          tstar_surft(l,n)**4
      END IF
    END DO
!$OMP END DO
  END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      olr(i,j) = olr(i,j) + sbcon * tstar_rad4(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

!CABLE_LSM: cable vn - make sure modified
  !-----------------------------------------------------------------------
  !     Specific humidity and temperature at 1.5 metres.
  !-----------------------------------------------------------------------
  CALL screen_tq_cbl (                                                         &
    land_pts,nsurft,                                                          &
    land_index,surft_index,surft_pts,flandg,                                  &
    sf_diag,chr1p5m,chr1p5m_sice,pstar,qim_1,resft,                           &
    tile_frac,tim_1,tstar_ssi,tstar_surft,                                    &
    z0hssi,z0h_surft,z0mssi,z0m_surft,z1,                                     &
    timestep,tstar_ssi_old,tstar_surft_old,                                   &
    l_co2_interactive, co2_mmr, co2_3d,                                       &
    f3_at_p, uStarGBM, rho1,                                                  &
    TScrnDcl_SSI,TScrnDcl_SURFT,tStbTrans,                                    &
    lq_mix_bl                                                                 &
    )

  ! Release space allocated for the transitional diagnostic.
  DEALLOCATE(tstar_ssi_old)
  IF ( .NOT. l_sice_new_code) DEALLOCATE(tstar_sic)

  !-----------------------------------------------------------------------
  ! 9.  Calculate surface latent heat flux.
  !-----------------------------------------------------------------------

  IF (sf_diag%slh) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,sf_diag,fqw_1,flandg,ei_land,ei_sice_sm)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        sf_diag%latent_heat(i,j) = lc * fqw_1(i,j)                            &
                                   + lf * (flandg(i,j) * ei_land(i,j) +       &
                                   (1.0 - flandg(i,j)) * ei_sice_sm(i,j))
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF


  !-----------------------------------------------------------------------
  ! Rescale FTL_1 as it should be used to update the botom row of the
  ! discrete equation handled by the new BL solver at the next (2nd)
  ! stage of the scheme.
  !-----------------------------------------------------------------------
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(ftl_1,tdims,cp)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      ftl_1(i,j) = ftl_1(i,j) / cp
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE ! L_correct = true: 2nd stage of the scheme

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(j,i)

  !-----------------------------------------------------------------------
  ! Rescale to Watts/m^2 as this is the final call to the imp BL solver
  ! and FTL_1 will be used by stash diagnostics
  !-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      ftl_1(i,j) = cp * ftl_1(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT
  !-----------------------------------------------------------------------
  !  U_V will be updated at 2nd stage of the scheme as the equations
  !  providing the implicit surface stresses have been modified
  !  consistently with the new scheme.
  !-----------------------------------------------------------------------
  ! U component of 10m wind
  IF (sf_diag%su10) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = udims%j_start,udims%j_end
      DO i = udims%i_start,udims%i_end
        sf_diag%u10m(i,j) = (u_1(i,j) + du_star1(i,j) + (du_1(i,j) -          &
                     cq_cm_u_1(i,j) * taux_1(i,j)) -                          &
                     u_0(i,j)) * cdr10m_u(i,j) + u_0(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  ! V component of 10m wind
  IF (sf_diag%sv10) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = vdims%j_start,vdims%j_end
      DO i = vdims%i_start,vdims%i_end
        sf_diag%v10m(i,j) = (v_1(i,j) + dv_star1(i,j) + (dv_1(i,j) -          &
                     cq_cm_v_1(i,j) * tauy_1(i,j)) -                          &
                     v_0(i,j)) * cdr10m_v(i,j) + v_0(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  ! Similar calculations for the neutral winds.
  IF (sf_diag%suv10m_n) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = udims%j_start,udims%j_end
      DO i = udims%i_start,udims%i_end
        sf_diag%u10m_n(i,j) = (u_1(i,j) + du_star1(i,j) + (du_1(i,j) -        &
                     cq_cm_u_1(i,j) * taux_1(i,j)) -                          &
                     u_0(i,j)) * sf_diag%cdr10m_n_u(i,j) + u_0(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO j = vdims%j_start,vdims%j_end
      DO i = vdims%i_start,vdims%i_end
        sf_diag%v10m_n(i,j) = (v_1(i,j) + dv_star1(i,j) + (dv_1(i,j) -        &
                     cq_cm_v_1(i,j) * tauy_1(i,j)) -                          &
                     v_0(i,j)) * sf_diag%cdr10m_n_v(i,j) + v_0(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  ! Correct surface stress diagnostics
!$OMP DO SCHEDULE(STATIC)
  DO j = udims%j_start,udims%j_end
    DO i = udims%i_start,udims%i_end
      taux_land(i,j) = taux_land(i,j) + taux_land_star(i,j)
      taux_ssi(i,j)  = taux_ssi(i,j)  + taux_ssi_star(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO j = vdims%j_start,vdims%j_end
    DO i = vdims%i_start,vdims%i_end
      tauy_land(i,j) = tauy_land(i,j) + tauy_land_star(i,j)
      tauy_ssi(i,j)  = tauy_ssi(i,j)  + tauy_ssi_star(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

END IF ! IF .NOT. L_correct

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_impl2_cbl
END MODULE sf_impl2_cbl_mod
