! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE sf_expl_l_cbl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_EXPL_L_cbl_MOD'

CONTAINS
!  SUBROUTINE SF_EXPL------------------------------------------------
!
!  Purpose: Calculate explicit surface fluxes of heat, moisture and
!           momentum. Also calculates surface exchange coefficients
!           required for implicit update of surface fluxes and surface
!           information required by the explicit boundary layer
!           routine
!
!
!  Documentation: UMDP 24.
!
!---------------------------------------------------------------------
!    Arguments :-
SUBROUTINE sf_expl_l_cbl (                                                        &
! IN date-related values
 curr_day_number,                                                             &
! IN values defining field dimensions and subset to be processed :
 land_pts, nice, nice_use,                                                    &
! IN  parameters for iterative SISL scheme
 numcycles, cycleno,                                                          &
! IN parameters required from boundary-layer scheme :
 bq_1,bt_1,z1_uv_ij,z1_uv_top,z1_tq,z1_tq_top,qw_1,tl_1,                      &
! IN soil/vegetation/land surface data :
 land_index,nsurft,sm_levels,canopy,catch,catch_snow,hcon_soilt,              &
 ho2r2_orog, flandg,                                                          &
 snow_surft,sil_orog_land,smvccl_soilt,smvcst_soilt,smvcwt_soilt,sthf_soilt,  &
 sthu_soilt,z0_surft,z0h_surft_bare, z0m_soil,                                &
! IN sea/sea-ice data :
 ice_fract_cat,k_sice,                                                        &
! IN input data from the wave model
 charnock_w,                                                                  &
! IN everything not covered so far :
 pstar,lw_down,sw_surft,zh,ddmfx,                                             &
 co2_mmr,co2_3d,l_co2_interactive,l_phenol,                                   &
 asteps_since_triffid,cs_pool_soilt,frac,canht_pft,photosynth_act_rad,lai_pft,&
 lq_mix_bl,t_soil_soilt,tsurf_elev_surft,ti,ti_cat,tstar,tstar_sea,           &
 tstar_sice_cat,tstar_surft,z_land,albsoil_soilt,cos_zenith_angle,            &
 l_aero_classic,l_dust,l_dust_diag,clay_soilt,o3,                             &
! IN idealised and SCM things
 l_spec_z0, z0m_scm, z0h_scm,                                                 &
! IN variables for message passing
 u_1_px, v_1_px, u_0_px, v_0_px,                                              &
! INOUT diagnostics
 sf_diag,                                                                     &
! INOUT data :
 z0msea,gs,g_leaf_acc,npp_pft_acc,resp_w_pft_acc,resp_s_acc_soilt,            &
! OUT Diagnostic not requiring STASH flags :
 recip_l_mo_sea,fqw_1,ftl_1,ftl_surft,                                        &
 radnet_sice,rhokm_1,rib,                                                     &
! OUT variables for message passing
 flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m,                            &
! OUT data required for tracer mixing :
 rho_aresist,aresist,resist_b,                                                &
 rho_aresist_surft,aresist_surft,resist_b_surft,                              &
! OUT data required for mineral dust scheme
 r_b_dust,cd_std_dust,u_s_std_surft,                                          &
! OUT data required elsewhere in UM system :
 fb_surf,u_s,t1_sd,q1_sd,                                                     &
! OUT data required elsewhere in boundary layer or surface code
 alpha1,alpha1_sea,alpha1_sice,ashtf_prime,ashtf_prime_sea,ashtf_prime_surft, &
 fqw_surft,epot_surft,fqw_ice,ftl_ice,fraca,rhostar,resfs,resft,              &
 rhokh,rhokh_surft,rhokh_sice,rhokh_sea,                                      &
 dtstar_surft,dtstar_sea,dtstar_sice,                                         &
 h_blend_orog,z0hssi,z0h_surft, z0mssi,z0m_surft,                             &
 z0m_eff,chr1p5m,chr1p5m_sice,smc_soilt,hcons_soilt,vshr,vshr_land,vshr_ssi,  &
 gpp,npp,resp_p,g_leaf,gpp_pft,npp_pft,                                       &
 resp_p_pft,resp_s_soilt,resp_s_tot_soilt,resp_l_pft,resp_r_pft,resp_w_pft,   &
 n_leaf,n_root,n_stem,lai_bal,gc_surft,canhc_surft,wt_ext_surft,flake,        &
 surft_index,surft_pts,tile_frac,fsmc_pft,emis_surft,emis_soil,               &
 !CABLE_LSM:{pass additional existing & CABLE state vars
 fland,                                                                       &
 air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,  ssnow_cbl, bgc_cbl,       & 
 bal_cbl, sum_flux_cbl, veg_cbl, soilin, soil_cbl )

USE sf_exch_cbl_mod,           ONLY: sf_exch_cbl
USE snowtherm_mod,         ONLY: snowtherm
USE gen_anthrop_heat_mod,  ONLY: generate_anthropogenic_heat
USE tilepts_mod,           ONLY: tilepts
USE heat_con_mod,          ONLY: heat_con
USE physiol_mod,           ONLY: physiol

!Use in relevant variables
USE theta_field_sizes, ONLY: t_i_length

USE dust_param, ONLY: ndiv
USE csigma
USE missing_data_mod

USE jules_soil_biogeochem_mod, ONLY:                                          &
! imported scalar parameters
   soil_model_rothc,                                                          &
! imported scalar variables (IN)
   soil_bgc_model

USE jules_soil_mod, ONLY: dzsoil, dzsoil_elev

USE jules_snow_mod, ONLY: nsmax

USE jules_sea_seaice_mod, ONLY: emis_sice, emis_sea

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype
USE fluxes, ONLY: anthrop_heat_surft, sw_sicat, sw_sea

USE jules_sea_seaice_mod, ONLY: l_ctile, buddy_sea, charnock, SeaSalinityFactor

USE jules_vegetation_mod, ONLY: can_model, can_rad_mod, ilayers, l_triffid

USE veg_param, ONLY: secs_per_360days

USE jules_surface_mod, ONLY: l_aggregate, formdrag, orog_drag_param,          &
                             fd_stab_dep, l_anthrop_heat_src
USE sf_diags_mod, ONLY: strnewsfdiag
USE timestep_mod, ONLY: timestep

#if defined(UM_JULES)
USE atm_step_local, ONLY: land_pts_trif, npft_trif,dim_cs1, dim_cs2,          &
     co2_dim_len,co2_dim_row
#else
USE ancil_info, ONLY: land_pts_trif, npft_trif,dim_cs1, dim_cs2,              &
     co2_dim_len,co2_dim_row
#endif

USE ancil_info, ONLY: sice_pts_ncat, sice_index_ncat,                         &
     ssi_index, sea_pts, sea_index, nsoilt, dim_cslayer

USE prognostics, ONLY: nsnow_surft, sice_surft, sliq_surft, snowdepth_surft,  &
                       tsnow_surft, ds_surft

USE c_elevate, ONLY: surf_hgt_surft, lw_down_elevcorr_surft,                  &
                      l_elev_absolute_height

USE bl_option_mod, ONLY: on, l_quick_ap2
USE solinc_data, ONLY: sky, l_skyview

USE atm_fields_bounds_mod, ONLY:                                              &
   pdims_s, pdims, tdims

USE ozone_vars, ONLY: flux_o3_pft, fo3_pft

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
USE cable_params_mod,         ONLY : soilin_type
USE cable_params_mod,         ONLY : soil_parameter_type

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs :-
!-----------------------------------------------------------------------
! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.
INTEGER, INTENT(IN) ::                                                        &
 curr_day_number,                                                             &
            ! IN current day of year
 land_pts,                                                                    &
            ! IN No of land points being processed.
 numcycles,                                                                   &
            ! Number of cycles (iterations) for iterative SISL.
 cycleno,                                                                     &
            ! Iteration no
 nice,                                                                        &
            ! Total number of sea ice categories
 nice_use   ! No. of sea ice categories used fully in surface calculations

! Defining vertical grid of model atmosphere.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN A buoyancy parameter
!                                  !    (beta q tilde).
,bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN A buoyancy parameter
!                                  !    (beta T tilde).
,z1_uv_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                             ! IN Height of lowest uv level (m).
,z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN Height of lowest tq level (m).
!                                  !    Note, if the grid used is
!                                  !    staggered in the vertical,
!                                  !    Z1_UV and Z1_TQ can be
!                                  !    different.
,qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN Total water content
,tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! IN Ice/liquid water temperature

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 u_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),         &
 v_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),         &
 u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),         &
 v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
! Charnock's coefficient from wave model

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
                    z1_uv_top(tdims%i_start:tdims%i_end,                      &
                              tdims%j_start:tdims%j_end)
                             ! Height of top of lowest uv-layer
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
                    z1_tq_top(tdims%i_start:tdims%i_end,                      &
                              tdims%j_start:tdims%j_end)
                             ! Height of top of lowest Tq-layer

! (c) Soil/vegetation/land surface parameters (mostly constant).
LOGICAL, INTENT(IN) ::                                                        &
 l_co2_interactive                                                            &
                             ! IN Switch for 3D CO2 field
,l_phenol                                                                     &
                             ! IN Indicates whether phenology
!                                  !    in use
,l_spec_z0
                             ! IN T if using prescribed
!                                  !    sea surface roughness lengths

INTEGER, INTENT(IN) ::                                                        &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
                                  !    point in ROW_LENGTH,ROWS is the
                                  !    land point.

INTEGER, INTENT(IN) ::                                                        &
 sm_levels                                                                    &
                             ! IN No. of soil moisture levels
,nsurft                                                                       &
                             ! IN No. of land-surface tiles
,asteps_since_triffid
                             ! IN Number of atmospheric
                                  !    timesteps since last call
                                  !    to TRIFFID.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 canopy(land_pts,nsurft)                                                      &
                             ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
,catch(land_pts,nsurft)                                                       &
                             ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
,catch_snow(land_pts,nsurft)                                                  &
                             ! IN Snow interception capacity of
!                                  !    tiles (kg/m2).
,hcon_soilt(land_pts,nsoilt)                                                  &
                             ! IN Soil thermal conductivity
!                                  !    (W/m/K).
,snow_surft(land_pts,nsurft)                                                  &
                             ! IN Lying snow on tiles (kg/m2)
,smvccl_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Critical volumetric SMC
!                                  !    (cubic m per cubic m of soil).
,smvcst_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Volumetric saturation point
!                                  !    (m3/m3 of soil).
,smvcwt_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Volumetric wilting point
!                                  !    (cubic m per cubic m of soil).
,sthf_soilt(land_pts,nsoilt,sm_levels)                                        &
                             ! IN Frozen soil moisture content of
!                                  !    each layer as a fraction of
!                                  !    saturation.
,sthu_soilt(land_pts,nsoilt,sm_levels)                                        &
                             ! IN Unfrozen soil moisture content
!                                  !    of each layer as a fraction of
!                                  !    saturation.
,z0_surft(land_pts,nsurft)                                                    &
                             ! IN Tile roughness lengths (m).
,z0h_surft_bare(land_pts,nsurft)                                              &
                             ! IN Tile thermal roughness lengths
                             ! without snow cover(m).
,z0m_soil(land_pts)                                                           &
                            ! IN bare soil momentum z0 (m).
,sil_orog_land(land_pts)                                                      &
                             ! IN Silhouette area of unresolved
!                                  !    orography per unit horizontal
!                                  !    area on land points only.
,ho2r2_orog(land_pts)                                                         &
                             ! IN Standard Deviation of orography.
!                                  !    equivilent to peak to trough
!                                  !    height of unresolved orography
,flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
!                                  ! IN Land fraction on all tiles.
!                                  !    divided by 2SQRT(2) on land
!                                  !    points only (m)
! (d) Sea/sea-ice data.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 ice_fract_cat(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end,nice_use)                            &
                       ! IN Fraction of gridbox covered by
!                      ! category sea-ice (decimal fraction).
                       ! If nice_use=1, this is the sum of the categories
,k_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)
                       ! IN 2 * Sea ice thermal conductivity divided
                       !  by surface layer thickness  (in coupled mode,
                       !  this is passed in from the sea ice model) (W/m2/K)

! (f) Atmospheric + any other data not covered so far, incl control.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)                   &
                             ! IN Surface pressure (Pascals).
,lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! IN Surface downward LW radiation
!                                  !    (W/m2).
,sw_surft(land_pts,nsurft)                                                    &
                             ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
,zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                             ! IN Height above surface of top of
!                            !    boundary layer (metres).
,ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
!                            ! IN Convective downdraught
!                            !    mass-flux at cloud base
,co2_mmr                                                                      &
                             ! IN CO2 Mass Mixing Ratio
,co2_3d(co2_dim_len,co2_dim_row)                                              &
!                                  ! IN 3D CO2 field if required.
,cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)                           &
                        ! IN Soil carbon (kg C/m2).
,frac(land_pts,ntype)                                                         &
                             ! IN Fractions of surface types.
,canht_pft(land_pts,npft)                                                     &
                             ! IN Canopy height (m)
,photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
!                                  ! IN Net downward shortwave radiation
!                                  !    in band 1 (w/m2).
,lai_pft(land_pts,npft)                                                       &
                             ! IN Leaf area index
,t_soil_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Soil temperatures (K).
,tsurf_elev_surft(land_pts,nsurft)                                            &
                             ! Tiled ice sub-surface temperature (K)
,ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                             ! IN Sea-ice surface layer
                             !    temperature (K).
,ti_cat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)             &
,tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! IN Open sea surface temperature (K)
,tstar_sice_cat(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end,nice_use)                           &
                             ! IN Sea-ice surface temperature (K).
,tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN GBM surface temperature (K).
,tstar_surft(land_pts,nsurft)                                                 &
                             ! IN Surface tile temperatures
,z_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! IN Land height (m).
,albsoil_soilt(land_pts,nsoilt)                                               &
!                                  ! Soil albedo.
, cos_zenith_angle(tdims%i_start:tdims%i_end,                                 &
                   tdims%j_start:tdims%j_end)                                 &
!                                  ! Cosine of the zenith angle
,clay_soilt(land_pts,nsoilt,dim_cslayer)                                      &
   ! Soil clay fraction.
,o3(land_pts)                                                                 &
                            ! IN Surface ozone concentration (ppb).
,z0m_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! IN Fixed Sea-surface roughness
!                                  !    length for momentum (m).(SCM)
,z0h_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! IN Fixed Sea-surface roughness
!                                  !    length for heat (m). (SCM)

LOGICAL, INTENT(IN) ::                                                        &
 l_aero_classic                                                               &
                             ! IN switch for using CLASSIC aerosol
                             !    scheme
,l_dust                                                                       &
                             ! IN switch for mineral dust
,l_dust_diag                                                                  &
                             ! IN Switch for diagnostic mineral dust
                             !    lifting
,lq_mix_bl
!-----------------------------------------------------------------------
!  In/outs :-
!-----------------------------------------------------------------------
!Diagnostics
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
 z0msea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! INOUT Sea-surface roughness
!                                  !       length for momentum (m).
,gs(land_pts)                                                                 &
                             ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
,g_leaf_acc(land_pts,npft)                                                    &
                             ! INOUT Accumulated G_LEAF
,npp_pft_acc(land_pts_trif,npft_trif)                                         &
!                                  ! INOUT Accumulated NPP_pft
,resp_w_pft_acc(land_pts_trif,npft_trif)                                      &
!                                  ! INOUT Accum RESP_W_pft
,resp_s_acc_soilt(land_pts_trif,nsoilt,dim_cslayer,dim_cs1)                   &
                                   ! INOUT Accumulated RESP_S
,n_leaf(land_pts,npft)                                                        &
                            ! OUT Leaf N content scaled by LAI
!                                 !     (kg N/m2).
,n_root(land_pts,npft)                                                        &
                            ! OUT Root N content scaled by LAI_bal
!                                 !     (kg N/m2).
,n_stem(land_pts,npft)                                                        &
                            ! OUT Stem N content scaled by LAI_bal
!                                 !     (kg N/m2).
,lai_bal(land_pts,npft)                                                       &
                            ! OUT LAI_bal
,resp_l_pft(land_pts,npft)                                                    &
                            ! OUT Leaf maintenance respiration
!                                 !     (kg C/m2/s).
,resp_r_pft(land_pts,npft)
                            ! OUT Root maintenance respiration
!                                 !     (kg C/m2/s).

!-----------------------------------------------------------------------
!  Outputs :-
!-----------------------------------------------------------------------
!-1 Diagnostic (or effectively so - includes coupled model requisites):-
INTEGER, INTENT(OUT) ::                                                       &
 surft_index(land_pts,ntype)                                                  &
                             ! OUT Index of tile points
,surft_pts(ntype)             ! OUT Number of tile points

!  (a) Calculated anyway (use STASH space from higher level) :-

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
!                                  ! OUT Reciprocal of the surface
!                                  !     Obukhov  length at sea
!                                  !     points. (m-1).
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! OUT Moisture flux between layers
!                                  !     (kg per square metre per sec).
!                                  !     FQW(,1) is total water flux
!                                  !     from surface, 'E'.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! OUT FTL(,K) contains net turbulent
!                                  !     sensible heat flux into layer K
!                                  !     from below; so FTL(,1) is the
!                                  !     surface sensible heat, H.(W/m2)
,ftl_surft(land_pts,nsurft)                                                   &
                             ! OUT Surface FTL for land tiles
,radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
             nice_use)                                                        &
                             ! OUT Surface net radiation on
!                                  !     sea-ice (W/m2)
,rhokm_1(pdims_s%i_start:pdims_s%i_end,                                       &
         pdims_s%j_start:pdims_s%j_end)                                       &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum on P-grid
,rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! OUT Mean bulk Richardson number for
!                                  !     lowest layer.
,rho_aresist(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)             &
                                   ! OUT RHOSTAR*CD_STD*VSHR
!                                  ! for CLASSIC aerosol scheme
,aresist(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                                   ! OUT 1/(CD_STD*VSHR)
!                                  ! for CLASSIC aerosol scheme
,resist_b(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                                   ! OUT (1/CH-1/(CD_STD)/VSHR
!                                  ! for CLASSIC aerosol scheme
,rho_aresist_surft(land_pts,nsurft)                                           &
!                                  ! OUT RHOSTAR*CD_STD*VSHR on land
!                                  ! tiles for CLASSIC aerosol scheme
,aresist_surft(land_pts,nsurft)                                               &
!                                  ! OUT 1/(CD_STD*VSHR) on land tiles
!                                  ! for CLASSIC aerosol scheme
,resist_b_surft(land_pts,nsurft)                                              &
!                                  ! OUT (1/CH-1/CD_STD)/VSHR on land
!                                  ! tiles for CLASSIC aerosol scheme
, r_b_dust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,ndiv)          &
                                   ! OUT surf layer res for dust
, cd_std_dust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)            &
                                   ! OUT Bulk transfer coef. for
!                                  ! momentum, excluding orographic
!                                  ! effects for mineral dust
, u_s_std_surft(land_pts,nsurft)                                              &
                                   ! OUT Surface friction velocity
!                                  ! (standard value)
!                                  ! for mineral dust
,emis_surft(land_pts,nsurft)                                                  &
                             ! OUT Emissivity for land tiles
,emis_soil(land_pts)
                             ! OUT Emissivity of underlying soil

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 flandfac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
 fseafac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
 rhokm_land(pdims_s%i_start:pdims_s%i_end,                                    &
            pdims_s%j_start:pdims_s%j_end),                                   &
 rhokm_ssi(pdims_s%i_start:pdims_s%i_end,                                     &
           pdims_s%j_start:pdims_s%j_end),                                    &
 cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

!-2 Genuinely output, needed by other atmospheric routines :-
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! OUT Surface flux buoyancy over
!                                  !     density (m^2/s^3)
,u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! OUT Surface friction velocity (m/s)
,t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! OUT Standard deviation of turbulent
!                                  !     fluctuations of layer 1 temp;
!                                  !     used in initiating convection.
,q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! OUT Standard deviation of turbulent
!                                  !     flucs of layer 1 humidity;
!                                  !     used in initiating convection.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 alpha1(land_pts,nsurft)                                                      &
                             ! OUT Mean gradient of saturated
!                                  !     specific humidity with respect
!                                  !     to temperature between the
!                                  !     bottom model layer and tile
!                                  !     surfaces
,alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! OUT ALPHA1 for sea.
,alpha1_sice(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                             ! OUT ALPHA1 for sea-ice.
,ashtf_prime(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                             ! OUT Coefficient to calculate
!                                  !     surface heat flux into sea-ice.
,ashtf_prime_sea(tdims%i_start:tdims%i_end,                                   &
                 tdims%j_start:tdims%j_end)                                   &
                             ! OUT Coefficient to calculate
!                                  !     surface heat flux into sea.
,ashtf_prime_surft(land_pts,nsurft)                                           &
                             ! OUT Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
,fqw_surft(land_pts,nsurft)                                                   &
                             ! OUT Surface FQW for land tiles
,epot_surft(land_pts,nsurft)                                                  &
                             ! OUT Local EPOT for land tiles.
,fqw_ice(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end,nice_use)                                  &
                             ! OUT Surface FQW for sea-ice
,ftl_ice(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end,nice_use)                                  &
                             ! OUT Surface FTL for sea-ice
,fraca(land_pts,nsurft)                                                       &
                             ! OUT Fraction of surface moisture
!                                  !     flux with only aerodynamic
!                                  !     resistance for snow-free land
!                                  !     tiles.
,rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! OUT Surface air density
,resfs(land_pts,nsurft)                                                       &
                             ! OUT Combined soil, stomatal
!                                  !     and aerodynamic resistance
!                                  !     factor for fraction (1-FRACA)
!                                  !     of snow-free land tiles.
,resft(land_pts,nsurft)                                                       &
                             ! OUT Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
,rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! OUT Grid-box surface exchange
!                                  !     coefficients
,rhokh_surft(land_pts,nsurft)                                                 &
                             ! OUT Surface exchange coefficients
!                                  !     for land tiles
,rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
                                                        nice_use)             &
                             ! OUT Surface exchange coefficients
!                                  !     for sea-ice
,rhokh_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! OUT Surface exchange coefficients
!                                  !     for sea
,dtstar_surft(land_pts,nsurft)                                                &
                             ! OUT Change in TSTAR over timestep
!                                  !     for land tiles
,dtstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! OUT Change is TSTAR over timestep
!                                  !     for open sea
,dtstar_sice(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,nice_use)                              &
                             ! OUT Change is TSTAR over timestep
!                                  !     for sea-ice
,h_blend_orog(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)            &
!                                  ! OUT Blending height used as part of
!                                  !     effective roughness scheme
,z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! OUT Roughness length for heat and
!                                  !     moisture over sea (m).
,z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! OUT Roughness length for momentum
!                                  !     over sea (m).
,z0h_surft(land_pts,nsurft)                                                   &
                             ! OUT Tile roughness lengths for heat
!                                  !     and moisture (m).
,z0m_surft(land_pts,nsurft)                                                   &
                             ! OUT Tile roughness lengths for
!                                  !     momentum.
,z0m_eff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! OUT Effective grid-box roughness
!                                  !     length for momentum
,chr1p5m(land_pts,nsurft)                                                     &
                             ! OUT Ratio of coefffs for
!                                  !     calculation of 1.5m temp for
!                                  !     land tiles.
,chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)            &
!                                  ! OUT CHR1P5M for sea and sea-ice
!                                  !     (leads ignored).
,smc_soilt(land_pts,nsoilt)                                                   &
                             ! OUT Available moisture in the
!                                  !     soil profile (mm).
,hcons_soilt(land_pts,nsoilt)                                                 &
                             ! OUT Soil thermal conductivity
!                                  !     including water and ice
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
,vshr_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
,vshr_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                             ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
,gpp(land_pts)                                                                &
                             ! OUT Gross primary productivity
!                                  !     (kg C/m2/s).
,npp(land_pts)                                                                &
                             ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
,resp_p(land_pts)                                                             &
                             ! OUT Plant respiration (kg C/m2/s).
,g_leaf(land_pts,nsurft)                                                        &
                             ! OUT Leaf turnover rate (/360days).
,gpp_pft(land_pts,nsurft)                                                       &
                             ! OUT Gross primary productivity
!                                  !     on PFTs (kg C/m2/s).
,npp_pft(land_pts,nsurft)                                                       &
                             ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
,resp_p_pft(land_pts,nsurft)                                                    &
                             ! OUT Plant respiration on PFTs
!                                  !     (kg C/m2/s).
,resp_s_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)                            &
                          ! OUT Soil respiration (kg C/m2/s).
,resp_s_tot_soilt(dim_cs2,nsoilt)                                             &
                            ! OUT Total soil respiration
                            ! (kg C/m2/s).
,resp_w_pft(land_pts,npft)                                                    &
                             ! OUT Wood maintenance respiration
!                                  !     (kg C/m2/s).
,gc_surft(land_pts,nsurft)                                                    &
                             ! OUT "Stomatal" conductance to
!                                  !      evaporation for land tiles
!                                  !      (m/s).
,canhc_surft(land_pts,nsurft)                                                 &
                             ! IN Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
,wt_ext_surft(land_pts,sm_levels,nsurft)                                      &
!                                  ! IN Fraction of evapotranspiration
!                                  !    which is extracted from each
!                                  !    soil layer by each tile.
,flake(land_pts,nsurft)                                                       &
                             ! IN Lake fraction.
,tile_frac(land_pts,nsurft)                                                   &
                             ! OUT Tile fractions including
!                                  !     snow cover in the ice tile.
,fsmc_pft(land_pts,npft)         ! OUT Moisture availability factor.
!-----------------------------------------------------------------------
! LOCAL variables
!-----------------------------------------------------------------------
!  Workspace :-
REAL(KIND=real_jlslsm) :: work_clay         ! working variable

REAL(KIND=real_jlslsm) ::                                                     &
 vfrac_surft(land_pts,nsurft)                                                 &
                          ! Fractional canopy coverage for
                          ! land tiles.
,radnet_surft(land_pts,nsurft)                                                &
                          ! Surface net radiation on tiles
,csnow(land_pts,nsmax)                                                        &
                          ! Areal heat capacity of snow (J/K/m2)
,ksnow(land_pts,nsmax)                                                        &
                          ! Thermal conductivity of snow (W/m/K)
,hcons_snow(land_pts,nsurft)                                                  &
                          ! Snow thermal conductivity
,resp_frac(dim_cs2,dim_cslayer)                                               &
                          ! respired fraction of RESP_S
,gc_stom_surft(land_pts,nsurft)                                               &
                             ! canopy conductance
,radnet_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Surface net radiation on open sea (W/m2)

!  Local scalars :-

INTEGER ::                                                                    &
 i,j,k,l,n,nn,m                                                               &
            ! LOCAL Loop counter (horizontal field index).
,IS,js                                                                        &
            ! Loop counter for coastal point stencil
,COUNT      ! Counter for average wind speed

REAL(KIND=real_jlslsm) ::                                                     &
 ushear                                                                       &
              ! U-component of surface-to-lowest-level wind shear.
,vshear                                                                       &
              ! V-component of surface-to-lowest-level wind shear.
,vshr2        ! Square of magnitude of surface-to-lowest-level
!                   ! wind shear.

REAL(KIND=real_jlslsm) :: seawind  ! average wind speed adjacent to coast
REAL(KIND=real_jlslsm) :: fseamax  ! Maximum factor to apply to coast wind speed

     ! Minimum factor allowed to convert coastal wind speed to land part
REAL(KIND=real_jlslsm) :: flandmin
PARAMETER(flandmin = 0.2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_EXPL_L_cbl'
!CABLE_LSM: Passed from surf_couple_
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

!H!REAL,  DIMENSION( tdims%i_end,tdims%j_end ) ::                                 &
!H!  true_latitude,   &
!H!  true_longitude
!H!!___UM soil/snow/radiation/met vars
!H!REAL,  DIMENSION(land_pts) :: & 
!H!  bexp_gb,    & ! => parameter b in Campbell equation 
!H!  hcon_gb,    & ! this is already passed as hcon? 
!H!  satcon_gb,  & ! hydraulic conductivity @ saturation [mm/s]
!H!  sathh_gb,   &
!H!  soil_alb
!H!REAL,  DIMENSION( tdims%i_end,tdims%j_end ) ::                                 &
!H!  ls_rain_cable,    &
!H!  ls_snow_cable
!H!REAL,  DIMENSION( tdims%i_end,tdims%j_end ) ::                                 &
!H!  cosz_gb, &
!H!  sin_theta_latitude
!CABLE_LSM: End

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Call TILEPTS to calculate surft_pts and surft_index for surface types
!-----------------------------------------------------------------------
CALL tilepts(land_pts,frac,surft_pts,surft_index)


!-----------------------------------------------------------------------
!   Generate the anthropogenic heat for surface calculations
!-----------------------------------------------------------------------
CALL generate_anthropogenic_heat( curr_day_number, land_pts, frac,            &
                                  l_anthrop_heat_src)


!-----------------------------------------------------------------------
! Calculate wind shear between level 1 and the surface
!-----------------------------------------------------------------------

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i,ushear,vshear,vshr2)                                        &
!$OMP SHARED(tdims,u_1_px,u_0_px,v_1_px,v_0_px,vshr_ssi,flandg,vshr_land,vshr)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    IF (flandg(i,j) <  1.0) THEN
      ushear = u_1_px(i,j) - u_0_px(i,j)
      vshear = v_1_px(i,j) - v_0_px(i,j)
      vshr2 = MAX (1.0e-6 , ushear * ushear + vshear * vshear)
      vshr_ssi(i,j) = SQRT(vshr2)
    ELSE
      vshr_ssi(i,j) = 0.0
    END IF

    IF (flandg(i,j) >  0.0) THEN
      vshr2 = MAX (1.0e-6 , u_1_px(i,j) * u_1_px(i,j)                         &
        + v_1_px(i,j) * v_1_px(i,j))
      vshr_land(i,j) = SQRT(vshr2)
    ELSE
      vshr_land(i,j) = 0.0
    END IF

    vshr(i,j)= flandg(i,j) * vshr_land(i,j)                                   &
      + (1.0 - flandg(i,j)) * vshr_ssi(i,j)
  END DO
END DO
!$OMP END PARALLEL DO

#if !defined(SCMA)

IF (l_ctile .AND. buddy_sea == on) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i,COUNT,is,js,ushear,vshear,vshr2,seawind,fseamax)            &
!$OMP SHARED(tdims,fseafac,flandfac,flandg,u_1_px,u_0_px,v_1_px,v_0_px,       &
!$OMP        vshr,vshr_ssi,vshr_land)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      fseafac(i,j)  = 1.0
      flandfac(i,j) = 1.0

      IF ( flandg(i,j) > 0.01 .AND. flandg(i,j) < 0.99 ) THEN
        !           !-----------------------------------------------------
        !           ! Calculate average windspeed over adjacent sea points
        !           !-----------------------------------------------------
        seawind = 0.0
        COUNT = 0
        DO IS = i - 1,i+1
          DO js = j - 1,j+1
            IF ( flandg(IS,js) < 0.001 ) THEN
              !               ! ie. this is basically a sea point
              ushear = u_1_px(IS,js) - u_0_px(IS,js)
              vshear = v_1_px(IS,js) - v_0_px(IS,js)
              vshr2 = MAX (1.0e-10 , ushear * ushear + vshear * vshear)
              seawind = seawind + SQRT( vshr2 )
              COUNT = COUNT + 1
            END IF
          END DO
        END DO
        !           !-----------------------------------------------------
        !           ! Calculate multiplicative factor, FSEAFAC, to convert
        !           ! from the GBM VSHR to an appropriate marine VSHR
        !           !-----------------------------------------------------
        IF (COUNT > 0) THEN
          seawind = seawind / REAL(COUNT)
          !             ! Restrict FSEAFAC so FLANDFAC>FLANDMIN
          fseamax = MIN( 1.0 / flandmin,                                      &
                        (1.0 - flandmin * flandg(i,j)) / (1.0 - flandg(i,j)) )
          !             ! First limit is to keep fseamax sensible as FLANDG -> 1
          !             ! Second limit is to keep fland > flandmin, remembering
          !             !   that the we want FLANDG-weighted sum of factors =1
          !             !   to preserve the gridbox mean VSHR
          fseafac(i,j) = MAX(1.0,                                             &
                         MIN( fseamax, seawind / vshr(i,j) ))
        END IF

        vshr_ssi(i,j) = vshr(i,j) * fseafac(i,j)

        flandfac(i,j) = ( 1.0 - fseafac(i,j) * (1.0 - flandg(i,j)) )          &
                      / flandg(i,j)
        vshr_land(i,j) = vshr(i,j) * flandfac(i,j)


        vshr(i,j)= flandg(i,j) * vshr_land(i,j)                               &
          + (1.0 - flandg(i,j)) * vshr_ssi(i,j)
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO
END IF  ! test on buddy_sea switch
#endif


!-----------------------------------------------------------------------
! Call physiology routine to calculate surface conductances and carbon
! fluxes.
!-----------------------------------------------------------------------
!CABLE_LSM:{ comment out
!C!CALL physiol (                                                                &
!C!  land_pts,land_index,                                                        &
!C!  sm_levels,nsurft,surft_pts,surft_index,                                     &
!C!  dim_cs1,dim_cs2,                                                            &
!C!  co2_mmr,co2_3d,co2_dim_len, co2_dim_row,l_co2_interactive,                  &
!C!  can_model,cs_pool_soilt,frac,canht_pft,photosynth_act_rad,                  &
!C!  lai_pft,pstar,qw_1,sthu_soilt,sthf_soilt,t_soil_soilt,tstar_surft,          &
!C!  smvccl_soilt,smvcst_soilt,smvcwt_soilt,vshr,z0_surft,z1_uv_ij,o3,           &
!C!  canhc_surft,vfrac_surft,emis_surft,emis_soil,flake,                         &
!C!  g_leaf,gs,gc_surft,gc_stom_surft,gpp,gpp_pft,npp,npp_pft,                   &
!C!  resp_p,resp_p_pft,resp_s_soilt,resp_l_pft,                                  &
!C!  resp_r_pft,resp_w_pft,n_leaf,                                               &
!C!  n_root,n_stem,lai_bal,                                                      &
!C!  smc_soilt,wt_ext_surft,fsmc_pft,                                            &
!C!  albsoil_soilt,cos_zenith_angle,                                             &
!C!  can_rad_mod,ilayers,flux_o3_pft,fo3_pft,sf_diag,asteps_since_triffid)
!CABLE_LSM:}
!----------------------------------------------------------------------
! If TRIFFID is being used apply any correction to the land-atmosphere
! fluxes on the first timestep after the last TRIFFID call. Such a
! correction will typically be associated with a total depletion of
! carbon or with maintanence of the seed fraction. The corrections
! are stored in the accumulation variables after the call to TRIFFID.
! The correction is added to the instantaneous land-atmosphere fluxes
! (so that the atmospheric carbon budget is corrected) but is not
! included in the accumulation variables which drive TRIFFID, since
! this has already been dealt with during the last TRIFFID call.
!----------------------------------------------------------------------
IF (l_triffid .AND. (asteps_since_triffid == 1)                               &
    .AND. ( cycleno == numcycles .OR. l_quick_ap2) ) THEN

  DO n = 1,nnpft
    DO l = 1,land_pts
      npp_pft(l,n) = npp_pft(l,n) + npp_pft_acc(l,n) / timestep
      resp_p_pft(l,n) = resp_p_pft(l,n) - npp_pft_acc(l,n) / timestep
      npp_pft_acc(l,n)=-npp_pft_acc(l,n)
    END DO
  END DO

  ! Here we have assumed that RothC must be used with TRIFFID, and is called
  ! on the same timestep.
  IF ( soil_bgc_model == soil_model_rothc ) THEN
    DO n = 1,dim_cs1
      DO l = 1,land_pts
        DO nn = 1,dim_cslayer
          !soil tiling is not compatible with triffid. OK to hard-code soil
          !tile index to 1 here
          resp_s_soilt(l,1,nn,n) = resp_s_soilt(l,1,nn,n)                     &
                                   + (resp_s_acc_soilt(l,1,nn,n) / timestep)
          resp_s_acc_soilt(l,1,nn,n) = -resp_s_acc_soilt(l,1,nn,n)
        END DO
      END DO
    END DO
  END IF

END IF

!----------------------------------------------------------------------
! Increment accumulation of leaf turnover rate.
! This is required for leaf phenology and/or TRIFFID, either of
! which can be enabled independently of the other.
!----------------------------------------------------------------------
IF ( cycleno == numcycles .OR. l_quick_ap2 ) THEN

  IF (l_phenol .OR. l_triffid) THEN
    DO n = 1,nnpft
      DO l = 1,land_pts
        g_leaf_acc(l,n) = g_leaf_acc(l,n) +                                   &
                          g_leaf(l,n) * ( timestep / secs_per_360days )
      END DO
    END DO
  END IF

  !----------------------------------------------------------------------
  ! Increment accumulation prognostics for TRIFFID
  !----------------------------------------------------------------------
  IF (l_triffid) THEN
    DO n = 1,nnpft
      DO l = 1,land_pts
        npp_pft_acc(l,n) = npp_pft_acc(l,n) + npp_pft(l,n) * timestep
        resp_w_pft_acc(l,n) = resp_w_pft_acc(l,n)                             &
                              + resp_w_pft(l,n) * timestep
      END DO
    END DO
  END IF

  IF ( soil_bgc_model == soil_model_rothc ) THEN
    DO n = 1,dim_cs1
      DO l = 1,land_pts
        DO nn = 1,dim_cslayer
          !soil tiling is not compatible with triffid. OK to hard-code soil
          !tile index to 1 here
          resp_s_acc_soilt(l,1,nn,n) = resp_s_acc_soilt(l,1,nn,n)             &
                                       + resp_s_soilt(l,1,nn,n) * timestep
        END DO
      END DO
    END DO
  END IF

END IF ! CycleNo == NumCycles

!-----------------------------------------------------------------------
! calculate CO2:(BIO+HUM) ratio, dependent on soil clay content, and
! sum soil respiration components
! (RESP_FRAC here then contains the fraction of soil respiration which
! is respired to the atmos. the rest is re-partitioned into BIO+HUM)

! resp_s_acc_soilt contains the full amount, and this is carried forward to
! VEG_CTL for use in updating soil carbon pools. RESP_S_TOT calculated
! here is passed to BL_TRMIX as the fraction which is respired as CO2
! to the atmosphere. RESP_S_TOT, and RESP_S are also passed out for
! storage in diagnostics 3293, and 3467-470.

!-----------------------------------------------------------------------
IF ( soil_bgc_model == soil_model_rothc ) THEN
  DO i = 1,land_pts
    !soil tiling is not compatible with triffid. OK to hard-code soil tile
    !index to 1 here by setting m = 1
    m = 1
    resp_s_tot_soilt(i,m) = 0.0
    DO nn = 1,dim_cslayer
      work_clay = EXP(-0.0786 * 100.0 * clay_soilt(i,m,nn))
      resp_frac(i,nn) = (3.0895+2.672 * work_clay) /                          &
                        (4.0895+2.672 * work_clay)
      resp_s_soilt(i,m,nn,1)  = resp_s_soilt(i,m,nn,1) * resp_frac(i,nn)
      resp_s_soilt(i,m,nn,2)  = resp_s_soilt(i,m,nn,2) * resp_frac(i,nn)
      resp_s_soilt(i,m,nn,3)  = resp_s_soilt(i,m,nn,3) * resp_frac(i,nn)
      resp_s_soilt(i,m,nn,4)  = resp_s_soilt(i,m,nn,4) * resp_frac(i,nn)
      resp_s_tot_soilt(i,m)   = resp_s_tot_soilt(i,m)                         &
                                + resp_s_soilt(i,m,nn,1)                      &
                                + resp_s_soilt(i,m,nn,2)                      &
                                + resp_s_soilt(i,m,nn,3)                      &
                                + resp_s_soilt(i,m,nn,4)
    END DO  !  layers
  END DO  !  points
END IF
!-----------------------------------------------------------------------
! Reset surft_pts and surft_index and set tile fractions to 1 if aggregate
! tiles are used (L_AGGREGATE=.T.).
! Otherwise, set tile fractions to surface type fractions.
!-----------------------------------------------------------------------
IF (l_aggregate) THEN
  surft_pts(1) = land_pts
  DO l = 1,land_pts
    tile_frac(l,1) = 1.0
    surft_index(l,1) = l
  END DO
ELSE
  DO n = 1,ntype
    DO l = 1, land_pts
      tile_frac(l,n) = frac(l,n)
    END DO
  END DO
END IF

IF (land_pts >  0) THEN    ! Omit if no land points

  !-----------------------------------------------------------------------
  ! Calculate the thermal conductivity of the top soil layer.
  !-----------------------------------------------------------------------
  DO m = 1, nsoilt
    CALL heat_con (land_pts,hcon_soilt,sthu_soilt(:,m,1),sthf_soilt(:,m,1),   &
                   smvcst_soilt(:,m,1),hcons_soilt(:,m))
  END DO

!CABLE_LSM:{ comment out - BUT hcons_snow needs definite value for sf_flux_cable
!C!  ! Thermal conductvity of top snow layer if nsmax > 0
!C!  IF (nsmax > 0) THEN
!C!    DO n = 1,nsurft
!C!      CALL snowtherm(land_pts,surft_pts(n),nsnow_surft(:,n),                  &
!C!                     surft_index(:,n),ds_surft(:,n,:),sice_surft(:,n,:),      &
!C!                     sliq_surft(:,n,:),csnow,ksnow)
!C!      DO l = 1,land_pts
!C!        hcons_snow(l,n) = ksnow(l,1)
!C!      END DO
!C!    END DO
!C!  END IF
!CABLE_LSM:}

END IF                     ! End test on land points

!-----------------------------------------------------------------------
! Calculate net radiation on land tiles, sea and sea-ice
!-----------------------------------------------------------------------

!CABLE_LSM: dont overwrite radnet_tile calculated by cable
!C!radnet_surft(:,:)  = 0.0
radnet_sea(:,:)    = 0.0
radnet_sice(:,:,:) = 0.0

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,k,j,i,n)                                                      &
!$OMP SHARED(surft_pts,surft_index,land_index,pdims,radnet_surft,sw_surft,    &
!$OMP        emis_surft,sky,lw_down,tstar_surft,sice_pts_ncat,sice_index_ncat,&
!$OMP        ssi_index,t_i_length,radnet_sice,sw_sicat,emis_sice,nsurft,      &
!$OMP        tstar_sice_cat,flandg,ice_fract_cat,l_skyview,                   &
!$OMP        nice_use, sea_pts, sea_index, radnet_sea, sw_sea, emis_sea,      &
!$OMP        tstar_sea)
IF (l_skyview) THEN
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / pdims%i_end + 1
      i = land_index(l) - (j-1) * pdims%i_end
!CABLE_LSM: dont overwrite radnet_tile calculated by cable
!C!      radnet_surft(l,n) = sw_surft(l,n) + emis_surft(l,n) *                   &
!C!        sky(i,j) * ( lw_down(i,j) - sbcon * tstar_surft(l,n)**4 )
    END DO
!$OMP END DO NOWAIT
  END DO
ELSE
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / pdims%i_end + 1
      i = land_index(l) - (j-1) * pdims%i_end
!CABLE_LSM: dont overwrite radnet_tile calculated by cable
!C!      radnet_surft(l,n) = sw_surft(l,n) + emis_surft(l,n) *                   &
!C!                 ( lw_down(i,j) - sbcon * tstar_surft(l,n)**4 )
    END DO
!$OMP END DO NOWAIT
  END DO
END IF

!$OMP DO SCHEDULE(STATIC)
DO k = 1, sea_pts
  l = sea_index(k)
  j = (ssi_index(l) - 1) / t_i_length + 1
  i = ssi_index(l) - (j-1) * t_i_length
  radnet_sea(i,j) = sw_sea(l) + emis_sea *                                    &
             (lw_down(i,j) - sbcon * tstar_sea(i,j)**4)
END DO
!$OMP END DO NOWAIT

DO n = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, sice_pts_ncat(n)
    l = sice_index_ncat(k,n)
    j=(ssi_index(l) - 1) / t_i_length + 1
    i = ssi_index(l) - (j-1) * t_i_length
    radnet_sice(i,j,n) = sw_sicat(l,n) + emis_sice *                          &
               ( lw_down(i,j) - sbcon * tstar_sice_cat(i,j,n)**4 )
  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL

IF (sf_diag%l_radnet_sea) THEN
  sf_diag%radnet_sea(:,:) = radnet_sea(:,:)
END IF

!-----------------------------------------------------------------------
! 4.  Surface turbulent exchange coefficients and "explicit" fluxes
!     (P243a, routine SF_EXCH).
!     Wind mixing "power" and some values required for other, later,
!     diagnostic calculations, are also evaluated if requested.
!-----------------------------------------------------------------------
!CABLE_LSM: call CABLE version of sf_exch
CALL sf_exch_cbl (                                                                &
 land_pts,nsurft,land_index,                                                  &
 surft_index,surft_pts,                                                       &
 flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
 nice, nice_use,                                                              &
 nsnow_surft,ds_surft,hcons_snow,k_sice,                                      &
 bq_1,bt_1,canhc_surft,canht_pft,lai_pft,canopy,catch,dzsoil(1),dzsoil_elev,  &
 flake,gc_surft,gc_stom_surft,hcons_soilt,                                    &
 can_model,catch_snow,lq_mix_bl,                                              &
 ho2r2_orog,ice_fract_cat,snowdepth_surft,snow_surft,pstar,qw_1,              &
 radnet_sea,radnet_sice,                                                      &
 radnet_surft,sil_orog_land,tile_frac,timestep,                               &
 surf_hgt_surft,l_elev_absolute_height,emis_surft,emis_soil,                  &
 tl_1,ti,ti_cat,t_soil_soilt(:,:,1),                                          &
 tsurf_elev_surft,                                                            &
 tsnow_surft,                                                                 &
 tstar_surft,tstar_sea,tstar_sice_cat,z_land,                                 &
 l_ctile,seasalinityfactor,                                                   &
 !IN input data from the wave model
 charnock_w,                                                                  &
 tstar,l_aggregate,l_spec_z0,z0m_scm,z0h_scm,                                 &
 l_aero_classic,l_dust,l_dust_diag,                                           &
 vfrac_surft,vshr_land,vshr_ssi,zh,ddmfx,                                     &
 z0_surft,z0h_surft_bare,z0m_soil,z1_uv_ij,z1_uv_top,z1_tq,z1_tq_top,         &
 sf_diag,formdrag,fd_stab_dep,                                                &
 orog_drag_param,z0msea,                                                      &
 lw_down,lw_down_elevcorr_surft,                                              &
 alpha1,alpha1_sea,alpha1_sice,ashtf_prime,ashtf_prime_sea,ashtf_prime_surft, &
 recip_l_mo_sea,cdr10m, chr1p5m,                                              &
 chr1p5m_sice,fqw_1,fqw_surft,epot_surft,fqw_ice,                             &
 ftl_1,ftl_surft,ftl_ice,fraca,h_blend_orog,charnock,                         &
 rhostar,resfs,resft,rib,                                                     &
 fb_surf,u_s,q1_sd,t1_sd,z0hssi,z0h_surft,                                    &
 z0mssi,z0m_surft,z0m_eff,rho_aresist,aresist,resist_b,                       &
 rho_aresist_surft,aresist_surft,resist_b_surft,                              &
 r_b_dust,cd_std_dust,u_s_std_surft,                                          &
 rhokh_surft,rhokh_sice,rhokh_sea,rhokm_1,rhokm_land,rhokm_ssi,               &
 dtstar_surft,dtstar_sea,dtstar_sice,rhokh,anthrop_heat_surft,                &
 !CABLE_LSM:{pass additional existing & CABLE state vars
 cycleno, numcycles, sm_levels,                                                & 
 sw_surft, co2_mmr, sthu_soilt(:,1,:), fland, curr_day_number,                 &
 air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
 ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
 soilin, soil_cbl )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sf_expl_l_cbl
END MODULE sf_expl_l_cbl_mod
