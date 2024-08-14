MODULE cable_land_sf_implicit_mod

USE sf_melt_mod, ONLY: sf_melt
USE screen_tq_mod, ONLY: screen_tq
USE sf_evap_mod, ONLY: sf_evap
USE sice_htf_mod, ONLY: sice_htf

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
                  ModuleName='CABLE_LAND_SF_IMPLICIT_MOD'

CONTAINS
!  SUBROUTINE CABLE_LAND_SF_IMPLICIT --------------------------------
!
!  Purpose: Calculate implicit correction for land point to surface
!           fluxes of heat,moisture and momentum, to be used by
!           the unconditionally stable and non-oscillatory BL
!           numerical solver.
!
!  edits associated with creation of CM3 denoted by CM3#56  
!--------------------------------------------------------------------
!    Arguments :-
SUBROUTINE cable_land_sf_implicit (                                            &
! IN values defining field dimensions and subset to be processed :
 land_pts,land_index,nsurft,surft_index,surft_pts,sm_levels,                   &
 canhc_surft,canopy,flake,smc_soilt,tile_frac,wt_ext_surft,fland,flandg,       &
! IN everything not covered so far :
 lw_down,sw_surft,t_soil_soilt,r_gamma,alpha1,ashtf_prime_surft,               &
 dtrdz_charney_grid_1,fraca,resfs,resft,rhokh_surft,                           &
 emis_surft,snow_surft,dtstar_surft,                                           &
! INOUT data :
 tstar_surft,fqw_surft,fqw_1,ftl_1,ftl_surft,sf_diag,                          &
! OUT Diagnostic not requiring STASH flags :
 ecan,ei_surft,esoil_surft,surf_ht_flux_land,ei_land,surf_htf_surft,           &
! OUT data required elsewhere in UM system :
 tstar_land,le_surft,radnet_surft,ecan_surft,esoil_soilt,                      &
 ext_soilt,melt_surft,tstar_surft_old,ERROR,                                   &
 !New arguments replacing USE statements
 ! lake_mod (IN)
 lake_h_ice_gb,                                                                &
 ! lake_mod (OUT)
 surf_ht_flux_lake_ij, non_lake_frac,                                          &
 ! fluxes (IN)
 anthrop_heat_surft,                                                           &
 ! fluxes (OUT)
 surf_ht_store_surft,                                                          &
 ! c_elevate (IN)
 lw_down_elevcorr_surft,                                                       &
 ! prognostics (IN)
 nsnow_surft,                                                                  &
 ! jules_mod (IN)
 snowdep_surft,                                                                &
 ! JULES Types containing field data (IN OUT)
 crop_vars,                                                                    &
!!!CABLE_LSM:
 progs_cbl, work_cbl, progs_cnp,                                               &
 !CM3:  
 cycleno, numcycles, npft, dim_cs1, smvcst,                           &
 ls_rain, ls_snow, con_rain, con_snow,                                         &
 tl_1, dtl1_1, qw_1, dqw1_1, ctctq1,                                           &
 canopy_gb, canopy_surft, smcl_soilt, sthf_soilt, sthu_soilt,                  &
 gs, gs_surft,                                                                 &
 npp_gb, npp_PFT, npp_acc_pft, gpp_gb, gpp_pft,                                &
 resp_s, resp_s_tot, &   !resp_s_tile, !Kathy intro-ed as diag         &
 resp_p, resp_p_pft, g_leaf_pft )
!!!CABLE_LSM: End 

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type

USE csigma,                   ONLY: sbcon

USE planet_constants_mod,     ONLY: cp

USE atm_fields_bounds_mod,    ONLY: tdims, pdims

USE theta_field_sizes,        ONLY: t_i_length, t_j_length

USE jules_surface_mod,        ONLY: l_aggregate, l_flake_model, ls

USE jules_snow_mod,           ONLY:                                            &
  nsmax, rho_snow_const, cansnowtile, l_snow_nocan_hc

USE jules_surface_types_mod,  ONLY: lake

USE sf_diags_mod,             ONLY: strnewsfdiag

#if defined(UM_JULES)
USE timestep_mod, ONLY: timestep, timestep_number
USE submodel_mod,   ONLY: atmos_im 
USE UM_parcore,       ONLY: mype
#else
USE model_time_mod, ONLY: timestep_number => timestep
USE model_time_mod, ONLY: timestep        => timestep_len  
USE model_grid_mod, ONLY: latitude, longitude
#endif

USE ancil_info,               ONLY: nsoilt

USE jules_surface_mod,        ONLY: l_neg_tstar

USE water_constants_mod,      ONLY: lc, lf, rho_ice

USE solinc_data,              ONLY: sky, l_skyview

USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook

USE jules_print_mgr,          ONLY: jules_message, jules_print

! In general CABLE utilizes a required subset of tbe JULES types, however;
USE progs_cbl_vars_mod, ONLY: progs_cbl_vars_type ! CABLE requires extra progs
USE work_vars_mod_cbl,  ONLY: work_vars_type      ! and some kept thru timestep
USE progs_cnp_vars_mod, ONLY: progs_cnp_vars_type ! CASA requires extra progs
!CABLE_LSM:Make avail call to CABLE implicit version
USE cable_implicit_main_mod,  ONLY: cable_implicit_main

IMPLICIT NONE
!--------------------------------------------------------------------
!  Inputs :-
! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.
!--------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
 land_pts    ! IN No of land points

! (c) Soil/vegetation/land surface parameters (mostly constant).
INTEGER, INTENT(IN) ::                                                         &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
                             !    point in ROW_LENGTH,ROWS is the
                             !    Ith land point.

INTEGER, INTENT(IN) ::                                                         &
 sm_levels                                                                     &
                             ! IN No. of soil moisture levels
,nsurft                                                                        &
                             ! IN No. of land tiles
,surft_index(land_pts,nsurft)                                                  &
                             ! IN Index of tile points
,surft_pts(nsurft)
                             ! IN Number of tile points

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 canhc_surft(land_pts,nsurft)                                                  &
                             ! IN Areal heat capacity of canopy
                             !    for land tiles (J/K/m2).
,canopy(land_pts,nsurft)                                                       &
                             ! IN Surface/canopy water for
                             !    snow-free land tiles (kg/m2)
,flake(land_pts,nsurft)                                                        &
                             ! IN Lake fraction.
,smc_soilt(land_pts,nsoilt)                                                    &
                             ! IN Available soil moisture (kg/m2).
,tile_frac(land_pts,nsurft)                                                    &
                             ! IN Tile fractions including
                             !    snow cover in the ice tile.
,wt_ext_surft(land_pts,sm_levels,nsurft)                                       &
                             ! IN Fraction of evapotranspiration
                             !    extracted from each soil layer
                             !    by each tile.
,fland(land_pts)                                                               &
                             ! IN Land fraction on land pts.
,flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN Land fraction on all pts.
,emis_surft(land_pts,nsurft)
                             ! IN Emissivity for land tiles
                             ! IN Lying snow on tiles (kg/m2)


! (f) Atmospheric + any other data not covered so far, incl control.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! IN Surface downward LW radiation
                                  !    (W/m2).
,sw_surft(land_pts,nsurft)
                             ! IN Surface net SW radiation on
                                  !    land tiles (W/m2).

REAL(KIND=real_jlslsm), INTENT(IN) :: r_gamma
                             ! IN implicit weight in level 1

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 alpha1(land_pts,nsurft)                                                       &
                             ! IN Mean gradient of saturated
                             !    specific humidity with respect
                             !    to temperature between the
                             !    bottom model layer and tile
                             !    surfaces
,ashtf_prime_surft(land_pts,nsurft)                                            &
                             ! IN Adjusted SEB coefficient for
                             !    land tiles.
,dtrdz_charney_grid_1(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end)                               &
                             ! IN -g.dt/dp for model layers.
,fraca(land_pts,nsurft)                                                        &
                             ! IN Fraction of surface moisture
                             !    flux with only aerodynamic
                             !    resistance for snow-free land
                             !    tiles.
,resfs(land_pts,nsurft)                                                        &
                             ! IN Combined soil, stomatal
                             !    and aerodynamic resistance
                             !    factor for fraction (1-FRACA) of
                             !    snow-free land tiles.
,resft(land_pts,nsurft)                                                        &
                             ! IN Total resistance factor.
                             !    FRACA+(1-FRACA)*RESFS for
                             !    snow-free land, 1 for snow.
,rhokh_surft(land_pts,nsurft)
                             ! IN Surface exchange coefficients
                             !    for land tiles


REAL(KIND=real_jlslsm), INTENT(OUT) :: t_soil_soilt(land_pts,nsoilt,sm_levels)
                             ! Soil temperatures (K).
REAL(KIND=real_jlslsm), INTENT(OUT) :: snow_surft(land_pts,nsurft)
!!!CABLE_LSM:
INTEGER, INTENT(IN) :: npft
INTEGER, INTENT(IN) :: dim_cs1
INTEGER, INTENT(IN) :: cycleno, numcycles
REAL, INTENT(IN) :: smvcst(land_pts,sm_levels)
                    ! IN Volumetric saturation point
!forcing
REAL, INTENT(IN) :: ls_rain(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: ls_snow(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: con_rain(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: con_snow(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    ! temperature
REAL, INTENT(IN) :: qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    ! humidity 
REAL, INTENT(IN) :: dtl1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  ! increment to temp
REAL, INTENT(IN) :: dqw1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  ! increment to humidity 
REAL, INTENT(IN) :: ctctq1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  ! temp/humidity increment
!re-declared as IN OUT as CABLE OUTs dtstar 
REAL, INTENT(IN OUT) ::  dtstar_surft(land_pts,nsurft)
                             ! Change in TSTAR over timestep
!!!CABLE_LSM: End 
 
!--------------------------------------------------------------------
!  In/outs :-
!--------------------------------------------------------------------
TYPE (strnewsfdiag), INTENT(IN OUT) :: sf_diag
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 tstar_surft(land_pts,nsurft)                                                  &
                             ! INOUT Surface tile temperatures
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! INOUT Moisture flux between layers
                             !       (kg per square metre per sec)
                             !       FQW(,1) is total water flux
                             !       from surface, 'E'.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! INOUT FTL(,K) contains net
                             !       turbulent sensible heat flux
                             !       into layer K from below; so
                             !       FTL(,1) is the surface
                             !       sensible heat, H.(W/m2)
,ftl_surft(land_pts,nsurft)                                                    &
                             ! INOUT Surface FTL for land tiles
,fqw_surft(land_pts,nsurft)
                             ! INOUT Surface FQW for land tiles

!!!CABLE_LSM:CM3
!CM3#56 Potentially problematically - while these additional variables are available
! some of the CNP vars would not be available in the implicit sectio.  This implies
! either partner changes through the UM boundary-layer scheme, USE of data modules
! and/or the work_cbl TYPE to pass stuff around in the longer term.
!
! As of 12/12/2023 - INH thinks it's okay
REAL(KIND=real_jlslsm), INTENT(INOUT) :: smcl_soilt(land_pts,nsoilt,sm_levels)
                                         ! INOUT Soil Moisture
REAL(KIND=real_jlslsm), INTENT(INOUT) :: sthf_soilt(land_pts,nsoilt,sm_levels)
                                         ! INOUT Soil Frozen fraction 
REAL(KIND=real_jlslsm), INTENT(INOUT) :: sthu_soilt(land_pts,nsoilt,sm_levels)
                                         ! INOUT Soil Unfrozen 
REAL(KIND=real_jlslsm), INTENT(INOUT) ::ext_soilt(land_pts,nsoilt,sm_levels)
                             ! OUT Extraction of water from each
                             !     soil layer (kg/m2/s).
REAL, INTENT(OUT) :: canopy_gb(land_pts)
REAL, INTENT(OUT) :: canopy_surft(land_pts, nsurft)
REAL, INTENT(OUT) :: gs( land_pts )
REAL, INTENT(OUT) :: gs_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: npp_acc_pft(land_pts,nsurft)
REAL, INTENT(OUT) :: npp_gb(land_pts)
REAL, INTENT(OUT) :: npp_pft(land_pts,nsurft)
REAL, INTENT(OUT) :: gpp_gb(land_pts)
REAL, INTENT(OUT) :: gpp_pft(land_pts,nsurft)
REAL, INTENT(OUT) :: resp_s(land_pts, dim_cs1) ! Soil respiration (kg C/m2/s)
REAL, INTENT(OUT) :: resp_s_tot(land_pts)       ! Total soil resp'n (kg C/m2/s)
REAL, INTENT(OUT) :: resp_p(land_pts)
REAL, INTENT(OUT) :: resp_p_pft(land_pts,nsurft)
REAL, INTENT(OUT) :: g_leaf_pft(land_pts,npft)
!CABLE TYPES containing field data (IN OUT)
TYPE(progs_cbl_vars_type), INTENT(IN OUT) :: progs_cbl
TYPE(work_vars_type),      INTENT(IN OUT) :: work_cbl
TYPE(progs_cnp_vars_type), INTENT(IN OUT) :: progs_cnp
!!!CABLE_LSM: End 

!--------------------------------------------------------------------
!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!--------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                             ! OUT Gridbox mean evaporation from
                             !     canopy/surface store (kg/m2/s).
                             !     Zero over sea.
,esoil_surft(land_pts,nsurft)                                                  &
                             ! OUT ESOIL for snow-free land tiles
,surf_ht_flux_land(tdims%i_start:tdims%i_end,                                  &
                   tdims%j_start:tdims%j_end)                                  &
                             ! OUT Net downward heat flux at
                             !     surface over land
                             !     fraction of gridbox (W/m2).
,ei_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! OUT Sublimation from lying snow
                             !     (kg/m2/s).
,surf_htf_surft(land_pts,nsurft)
                             ! OUT Net downward surface heat flux
                             !     on tiles (W/m2)

!-2 Genuinely output, needed by other atmospheric routines :-

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 tstar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! OUT   Land mean sfc temperature (K)
,le_surft(land_pts,nsurft)                                                     &
                             ! OUT Surface latent heat flux for
                             !     land tiles
,radnet_surft(land_pts,nsurft)                                                 &
                             ! OUT Surface net radiation on
                             !     land tiles (W/m2)
,ei_surft(land_pts,nsurft)                                                     &
                             ! OUT EI for land tiles.
,ecan_surft(land_pts,nsurft)                                                   &
                             ! OUT ECAN for snow-free land tiles
,esoil_soilt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nsoilt)       &
                             ! OUT Surface evapotranspiration
                             !     from soil moisture store
                             !     (kg/m2/s).
,melt_surft(land_pts,nsurft)                                                   &
                             ! OUT Snowmelt on land tiles (kg/m2/s
,tstar_surft_old(land_pts,nsurft)                                              &
                             ! OUT Tile surface temperatures at
                             !     beginning of timestep.
,non_lake_frac(land_pts)
                             ! OUT total tile fraction for surface types
                             ! other than inland water

INTEGER, INTENT(OUT) ::                                                        &
 ERROR                       ! OUT 0 - AOK;
                             !     1 to 7  - bad grid definition detected;

!New arguments replacing USE statements
! lake_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: lake_h_ice_gb(land_pts)
! lake_mod (OUT)
REAL(KIND=real_jlslsm)              :: surf_ht_flux_lake_ij(t_i_length,t_j_length)
! fluxes (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: anthrop_heat_surft(land_pts,nsurft)
! fluxes (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: surf_ht_store_surft(land_pts,nsurft)
! c_elevate (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: lw_down_elevcorr_surft(land_pts,nsurft)
! prognostics (IN)
INTEGER, INTENT(IN) :: nsnow_surft(land_pts,nsurft)
! jules_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: snowdep_surft(land_pts,nsurft)

!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars

!--------------------------------------------------------------------
!  Workspace :-
!--------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 elake_surft(land_pts,nsurft)                                                  &
                             ! Lake evaporation.
,melt_ice_surft(land_pts,nsurft)                                               &
                             ! Ice melt on FLake lake tile (kg/m2/s)
,lake_ice_mass(land_pts)                                                       &
                             ! areal density equivalent to
                             ! lake ice of a given depth (kg/m2)
,snowmelt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Snowmelt (kg/m2/s).

REAL(KIND=real_jlslsm) ::                                                      &
 canhc_surf(land_pts)
                             ! Areal heat capacity of canopy
                             ! for land tiles (J/K/m2).

!  Local scalars :-

INTEGER ::                                                                     &
 i,j                                                                           &
            ! LOCAL Loop counter (horizontal field index).
,k                                                                             &
            ! LOCAL Tile pointer
,l                                                                             &
            ! LOCAL Land pointer
,n                                                                             &
            ! LOCAL Loop counter (tile index).
,m
            ! Loop counter for soil tiles

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CABLE_LAND_SF_IMPLICIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------

! Calculate surface scalar fluxes, temperatures only at the 1st call
! of the subroutine (first stage of the new BL solver) using standard
! MOSES2 physics and equations. These are the final values for this
! timestep and there is no need to repeat the calculation.
!-----------------------------------------------------------------------

ERROR = 0

!-----------------------------------------------------------------------
! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
!-----------------------------------------------------------------------

DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    ftl_surft(l,n) = cp * ftl_surft(l,n)
  END DO
END DO

!-----------------------------------------------------------------------
! Land surface calculations
!-----------------------------------------------------------------------
!CABLE_LSM:likely unecessary initialization of radnet_surft
! initialise diagnostics to 0 to avoid packing problems
DO n = 1, nsurft
  DO l = 1, land_pts
    radnet_surft(l,n) = 0.0
    le_surft(l,n) = 0.0
  END DO
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
!CABLE_LSM: End

!CM3#56 - need to take a copy of tstart_surft prior to CABLE
! not needed to evaluate surf_ht_store_surft BUT it is an output needed by the UM
DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    tstar_surft_old(l,n) = tstar_surft(l,n)
  END DO
!$OMP END DO NOWAIT
END DO

!CM3#56 - use of REAL() in the call
!       - ideally map %fe to the output le_surft in the unpack
!       - check whether dtstar_surft is needed
!           - this impacts on the CABLE redeclaration so cannot be a long-term fix
!       - at some point we will need to sort out the CASA variables
call cable_implicit_main( tdims%i_end, tdims%j_end, land_pts, nsurft,  npft,   &
                          sm_levels, dim_cs1, cycleno, numcycles,     &
                          timestep, timestep_number, land_index,       &
                          surft_pts, surft_index,   &
                          Fland, tile_frac, smvcst, &
                          ls_rain, ls_snow, con_rain, con_snow,                &
                          tl_1, dtl1_1, qw_1, dqw1_1, ctctq1,                  &
                          canopy_gb, canopy_surft, t_soil_soilt(:,1,:),        &
                          smcl_soilt(:,1,:),  sthf_soilt(:,1,:),               &
                          sthu_soilt(:,1,:), snow_surft,                       & 
                          !returned fluxes etc
                          ftl_1, ftl_surft, fqw_1, fqw_surft, le_surft,        &
                          tstar_surft, dtstar_surft, surf_ht_flux_land, surf_htf_surft,      &
                          ecan_surft, esoil_surft, ei_surft, radnet_surft,     &
                          gs, gs_surft, sf_diag% t1p5m_surft,                  &
                          sf_diag% q1p5m_surft, melt_surft,                    &
                          npp_gb, npp_pft, npp_acc_pft, gpp_gb, gpp_pft,       &
                          resp_s, resp_s_tot, resp_p, resp_p_pft, g_leaf_pft,  &
                          progs_cbl, work_cbl )                                                          

!!CABLE_LSM:End

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,n,j,i,k)                                                       &
!$OMP SHARED(tdims,nsurft,surft_pts,surft_index,                               &
!$OMP        ftl_surft,nsoilt,land_pts,t_soil_soilt,                           &
!$OMP        tstar_surft_old,tstar_surft,dtstar_surft,cp,error,tile_frac,      &
!$OMP        non_lake_frac,lake,l_flake_model,l_aggregate)


!-----------------------------------------------------------------------
! Optional error check : test for negative top soil layer temperature
!-----------------------------------------------------------------------
IF (l_neg_tstar) THEN
  DO m = 1,nsoilt
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,land_pts
      IF (t_soil_soilt(l,m,1) < 0) THEN
        ERROR = 1
        WRITE(jules_message,*)                                                 &
              '*** ERROR DETECTED BY ROUTINE JULES_LAND_SF_IMPLICIT ***'
        CALL jules_print('jules_land_sf_implicit_jls',jules_message)
        WRITE(jules_message,*) 'NEGATIVE TEMPERATURE IN TOP SOIL LAYER AT '
        CALL jules_print('jules_land_sf_implicit_jls',jules_message)
        WRITE(jules_message,*) 'LAND POINT ',l
        CALL jules_print('jules_land_sf_implicit_jls',jules_message)
      END IF
    END DO
!$OMP END DO NOWAIT
  END DO
END IF

!-----------------------------------------------------------------------
!   Diagnose the land surface temperature
!-----------------------------------------------------------------------

!CM3#56 - if needed this evaluation of tstart_surft_old should be before cable.
!DO n = 1,nsurft
!!$OMP DO SCHEDULE(STATIC)
!  DO k = 1,surft_pts(n)
!    l = surft_index(k,n)
!    tstar_surft_old(l,n) = tstar_surft(l,n)
!!CABLE_LSM: use CABLE's tstar_tile
!!C!      tstar_surft(l,n) = tstar_surft_old(l,n) + dtstar_surft(l,n)
!  END DO
!!$OMP END DO NOWAIT
!END DO

!-----------------------------------------------------------------------
!   Calculate non_lake_frac
!-----------------------------------------------------------------------
!CM3#56 - this needs some thought and follow through as to where non_lake_frac is used
! - this appears to be linked solely to use of FLAKE so should remain as 1.0
!$OMP DO SCHEDULE(STATIC)
DO l = 1,land_pts
  ! initialise the non-lake fraction to one, not zero,
  ! in case there should ever be more than one lake tile, see below
  non_lake_frac(l) = 1.0
END DO
!$OMP END DO NOWAIT

!IF ( ( l_flake_model ) .AND. ( .NOT. l_aggregate) ) THEN
!!$OMP DO SCHEDULE(STATIC)
!  DO l = 1,land_pts
!    ! Remove FLake tile fraction.
!    non_lake_frac(l) = non_lake_frac(l) - tile_frac(l,lake)
!  END DO
!!$OMP END DO NOWAIT
!END IF

!$OMP END PARALLEL

!-----------------------------------------------------------------------
! 7.  Surface evaporation components and updating of surface
!     temperature (P245, routine SF_EVAP).
!-----------------------------------------------------------------------
!CABLE_LSM:
!CM2!CALL sf_evap (                                                                 &
!CM2!  land_pts,nsurft,                                                             &
!CM2!  land_index,surft_index,surft_pts,sm_levels,fland,                            &
!CM2!  ashtf_prime_surft,canopy,dtrdz_charney_grid_1,flake,fraca,                   &
!CM2!  snow_surft,resfs,resft,rhokh_surft,tile_frac,smc_soilt,wt_ext_surft,         &
!CM2!  timestep,r_gamma,fqw_1,fqw_surft,ftl_1,ftl_surft,tstar_surft,                &
!CM2!  ecan,ecan_surft,elake_surft,esoil_soilt,esoil_surft,ei_surft,ext_soilt,      &
!CM2!  sf_diag, non_lake_frac,                                                      &
!CM2!  ! crop_vars_mod (IN)
!CM2!  crop_vars%frac_irr_soilt, crop_vars%frac_irr_surft,                          &
!CM2!  crop_vars%wt_ext_irr_surft, crop_vars%resfs_irr_surft,                       &
!CM2!  ! crop_vars_mod (IN OUT)
!CM2!  crop_vars%smc_irr_soilt,                                                     &
!CM2!  ! crop_vars_mod (OUT)
!CM2!  crop_vars%ext_irr_soilt)

DO N=1,Nsurft
  DO L=1,LAND_PTS
    ELAKE_surft(L,N) = 0.
  ENDDO
ENDDO

DO j=tdims%j_start,tdims%j_end
 DO i=tdims%i_start,tdims%i_end
  ecan(i,j) = 0.
  esoil_soilt(i,j,1) = 0.
 ENDDO
ENDDO

DO N=1,Nsurft
 DO K=1,surft_PTS(N)
   L = surft_INDEX(K,N)
   j=(land_index(l)-1)/tdims%i_end + 1
   i = land_index(l) - (j-1)*tdims%i_end
   ecan(i,j) = ecan(i,j) + tile_frac(l,n)*ecan_surft(l,n)
   esoil_soilt(i,j,1) = esoil_soilt(i,j,1) + tile_frac(l,n)*esoil_surft(l,n)
 ENDDO
ENDDO   
!CABLE_LSM:End

!-----------------------------------------------------------------------
!     Surface melting of sea-ice and snow on land tiles.
!-----------------------------------------------------------------------

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,n,j,i)                                                         &
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

!CABLE_LSM:
DO n = 1,nsurft
!CM2!  CALL sf_melt (                                                               &
!CM2!    land_pts,land_index,                                                       &
!CM2!    surft_index(:,n),surft_pts(n),flandg,                                      &
!CM2!    alpha1(:,n),ashtf_prime_surft(:,n),dtrdz_charney_grid_1,                   &
!CM2!    resft(:,n),rhokh_surft(:,n),tile_frac(:,n),timestep,r_gamma,               &
!CM2!    ei_surft(:,n),fqw_1,ftl_1,fqw_surft(:,n),ftl_surft(:,n),                   &
!CM2!    tstar_surft(:,n),snow_surft(:,n),snowdep_surft(:,n),                       &
!CM2!    melt_surft(:,n)                                                            &
!CM2!    )
!CM2!
!CM2!  !-----------------------------------------------------------------------
!CM2!  ! thermodynamic, flux contribution of melting ice on the FLake lake tile
!CM2!  !-----------------------------------------------------------------------
!CM2!  IF (     (l_flake_model   )                                                  &
!CM2!    .AND. ( .NOT. l_aggregate)                                                 &
!CM2!    .AND. (n == lake       ) ) THEN
!CM2!
!CM2!    ! lake_h_ice_gb is only initialised if FLake is on.
!CM2!
!CM2!!$OMP PARALLEL DO                                                              &
!CM2!!$OMP SCHEDULE(STATIC)                                                         &
!CM2!!$OMP DEFAULT(NONE)                                                            &
!CM2!!$OMP PRIVATE(l)                                                               &
!CM2!!$OMP SHARED(land_pts,lake_ice_mass,lake_h_ice_gb)
!CM2!    DO l = 1, land_pts
!CM2!      lake_ice_mass(l) = lake_h_ice_gb(l) * rho_ice
!CM2!    END DO
!CM2!!$OMP END PARALLEL DO
!CM2!
!CM2!    CALL sf_melt (                                                             &
!CM2!      land_pts,land_index,                                                     &
!CM2!      surft_index(:,n),surft_pts(n),flandg,                                    &
!CM2!      alpha1(:,n),ashtf_prime_surft(:,n),dtrdz_charney_grid_1,                 &
!CM2!      resft(:,n),rhokh_surft(:,n),tile_frac(:,n),timestep,r_gamma,             &
!CM2!      ei_surft(:,n),fqw_1,ftl_1,fqw_surft(:,n),ftl_surft(:,n),                 &
!CM2!      tstar_surft(:,n),lake_ice_mass,lake_ice_mass / rho_snow_const,           &
!CM2!      melt_ice_surft(:,n)                                                      &
!CM2!        )
!CM2!  END IF

  !-----------------------------------------------------------------------
  !  Increment snow by sublimation and melt
  !-----------------------------------------------------------------------

!$OMP PARALLEL DO                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(k,l,j,i)                                                         &
!$OMP SHARED(surft_pts,surft_index,land_index,t_i_length,ei_land,tile_frac,    &
!$OMP        ei_surft,snowmelt,melt_surft,n)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    ei_land(i,j) = ei_land(i,j) + tile_frac(l,n) * ei_surft(l,n)
    snowmelt(i,j) = snowmelt(i,j) +                                            &
                    tile_frac(l,n) * melt_surft(l,n)
  END DO
!$OMP END PARALLEL DO

END DO

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(SHARED)                                                          &
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

!CM3#56 - remove FLAKE model
!IF (     (l_flake_model   )                                                    &
!    .AND. ( .NOT. l_aggregate) ) THEN
!!$OMP DO SCHEDULE(STATIC)
!  DO j = tdims%j_start,tdims%j_end
!    DO i = tdims%i_start,tdims%i_end
!      surf_ht_flux_lake_ij(i,j) = 0.0
!    END DO
!  END DO
!!$OMP END DO NOWAIT
!END IF

!$OMP DO SCHEDULE(STATIC)
DO l = 1,land_pts
  j=(land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length
  tstar_land(i,j) = 0.0
END DO
!$OMP END DO NOWAIT

!CABLE_LSM:
! initialise diagnostics to 0 to avoid packing problems
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

!CM3#56 - remove option for skyview
!IF (l_skyview) THEN
!CABLE_LSM:
!CM2!  DO n = 1,nsurft
!CM2!!$OMP DO SCHEDULE(STATIC)
!CM2!    DO k = 1,surft_pts(n)
!CM2!      l = surft_index(k,n)
!CM2!      j=(land_index(l) - 1) / tdims%i_end + 1
!CM2!      i = land_index(l) - (j-1) * tdims%i_end
!CM2!      radnet_surft(l,n) = sw_surft(l,n) +   emis_surft(l,n) *                  &
!CM2!        sky(i,j) * ( lw_down(i,j) + lw_down_elevcorr_surft(l,n)                &
!CM2!                                - sbcon * tstar_surft(l,n)**4 )
!CM2!    END DO
!CM2!!$OMP END DO
!CM2!  END DO
!  IF (sf_diag%l_lw_surft) THEN
!    DO n = 1,nsurft
!!$OMP DO SCHEDULE(STATIC)
!      DO k = 1,surft_pts(n)
!        l = surft_index(k,n)
!        j=(land_index(l) - 1) / tdims%i_end + 1
!        i = land_index(l) - (j-1) * tdims%i_end
!        sf_diag%lw_up_surft(l,n)   = emis_surft(l,n) * sky(i,j) *              &
!                                     sbcon * tstar_surft(l,n)**4               &
!                                   + (1.0 - emis_surft(l,n)) *                 &
!                                     sky(i,j) * (lw_down(i,j) +                &
!                                     lw_down_elevcorr_surft(l,n))
!        sf_diag%lw_down_surft(l,n) = sky(i,j) * (lw_down(i,j) +                &
!                                     lw_down_elevcorr_surft(l,n))
!      END DO
!!$OMP END DO
!    END DO
!  END IF
!ELSE
!CABLE_LSM:
!!  DO n = 1,nsurft
!!!$OMP DO SCHEDULE(STATIC)
!!    DO k = 1,surft_pts(n)
!!      l = surft_index(k,n)
!!      j=(land_index(l) - 1) / tdims%i_end + 1
!!      i = land_index(l) - (j-1) * tdims%i_end
!!      radnet_surft(l,n) = sw_surft(l,n) +   emis_surft(l,n) *                  &
!!                 ( lw_down(i,j) + lw_down_elevcorr_surft(l,n)                  &
!!                                - sbcon * tstar_surft(l,n)**4 )
!!    END DO
!!!$OMP END DO
!!  END DO
  IF (sf_diag%l_lw_surft) THEN
    DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / tdims%i_end + 1
        i = land_index(l) - (j-1) * tdims%i_end
        sf_diag%lw_up_surft(l,n)   = emis_surft(l,n) * sbcon *                 &
                                     tstar_surft(l,n)**4                       &
                                   + (1.0 - emis_surft(l,n)) *                 &
                                     (lw_down(i,j) +                           &
                                     lw_down_elevcorr_surft(l,n))
        sf_diag%lw_down_surft(l,n) = lw_down(i,j) +                            &
                                     lw_down_elevcorr_surft(l,n)
      END DO
!$OMP END DO
    END DO
  END IF
!END IF !CM3#56

!CM3#56 - remove surf_ht_store_surft and look to fill inside CABLE if needed
!       - remove FLAKE code
!       - remove evaluation of any tiled flux and look to fill from CABLE
DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    !canhc_surf(l) = canhc_surft(l,n)
    !IF ( ( .NOT. cansnowtile(n)) .AND. l_snow_nocan_hc .AND.                   &
    !     (nsmax > 0) .AND. (nsnow_surft(l,n) > 0) ) canhc_surf(l) = 0.0

    !CM3#56 - set to zero for now - should be coming from CABLE anyway
    surf_ht_store_surft(l,n) = 0.0
!    surf_ht_store_surft(l,n) = (canhc_surf(l) / timestep) *                    &
!                         (tstar_surft(l,n) - tstar_surft_old(l,n))
!!!CABLE_LSM: Replaced with CABLE field 
!!    surf_htf_surft(l,n) = radnet_surft(l,n) + anthrop_heat_surft(l,n) -        &
!!                        ftl_surft(l,n) -                                       &
!!                        le_surft(l,n) -                                        &
!!                        lf * (melt_surft(l,n) + melt_ice_surft(l,n)) -         &
!!                        surf_ht_store_surft(l,n)
    ! separate out the lake heat flux for FLake
    ! and replace the snow-melt (NSMAX=0 only) and ice-melt heat fluxes
    ! so Flake can do its melting
    !IF (     (l_flake_model   )                                                &
    !    .AND. ( .NOT. l_aggregate)                                             &
    !    .AND. (n == lake       ) ) THEN
    !  IF (nsmax == 0) THEN
    !    surf_ht_flux_lake_ij(i,j) = surf_htf_surft(l,n)                        &
    !                  + lf * (melt_surft(l,n) + melt_ice_surft(l,n))
    !  ELSE
    !    surf_ht_flux_lake_ij(i,j) = surf_htf_surft(l,n)                        &
    !                  + lf * melt_ice_surft(l,n)
    !  END IF
    !ELSE
      surf_ht_flux_land(i,j) = surf_ht_flux_land(i,j)                          &
                        + tile_frac(l,n) * surf_htf_surft(l,n)
    !END IF
    tstar_land(i,j) = tstar_land(i,j)                                          &
               + tile_frac(l,n) * tstar_surft(l,n)
  END DO
!$OMP END DO
END DO

  ! normalise the non-lake surface heat flux
!CABLE_LSM: 
!CM2!IF ( l_flake_model .AND. ( .NOT. l_aggregate) ) THEN
!CM2!!$OMP DO SCHEDULE(STATIC)
!CM2!  DO l = 1,land_pts
!CM2!    j=(land_index(l) - 1) / t_i_length + 1
!CM2!    i = land_index(l) - (j-1) * t_i_length
!CM2!    ! be careful about gridboxes that are all lake
!CM2!    IF (non_lake_frac(l) > EPSILON(0.0)) THEN
!CM2!      surf_ht_flux_land(i,j) = surf_ht_flux_land(i,j) / non_lake_frac(l)
!CM2!    END IF
!CM2!  END DO
!CM2!!$OMP END DO
!CM2!END IF

IF (sf_diag%l_lh_land) THEN
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    sf_diag%lh_land(l) = SUM((tile_frac(l,:) * le_surft(l,:)))
  END DO
!$OMP END DO
END IF

!-----------------------------------------------------------------------
! Optional error check : test for negative surface temperature
!-----------------------------------------------------------------------
!CM3#56 revise to refer to CABLE not JULES
IF (l_neg_tstar) THEN
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    IF (tstar_land(i,j) < 0) THEN
      ERROR = 1
      WRITE(jules_message,*)                                                   &
           '*** ERROR DETECTED BY ROUTINE CABLE_LAND_SF_IMPLICIT ***'
      CALL jules_print('cable_land_sf_implicit_cbl',jules_message)
      WRITE(jules_message,*) 'NEGATIVE SURFACE TEMPERATURE AT LAND POINT ',l
      CALL jules_print('cable_land_sf_implicit_cbl',jules_message)
    END IF
  END DO
!$OMP END DO
END IF

!$OMP END PARALLEL


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE cable_land_sf_implicit
END MODULE cable_land_sf_implicit_mod
