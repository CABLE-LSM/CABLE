MODULE cable_implicit_main_mod
  
CONTAINS

SUBROUTINE Cable_implicit_main( row_length, rows, land_pts, nsurft, npft,      &
                                sm_levels, dim_cs1, cycleno, numcycles,        &
                                timestep, timestep_number, land_index,         &
                                surft_pts, surft_index,                        & 
                                Fland, tile_frac, smvcst,                      &
                                ls_rain, ls_snow, conv_rain, conv_snow,        &
                                tl_1, dtl1_1, qw_1, dqw1_1, ctctq1,            &
                                canopy_gb, canopy_surft, t_soil,               & 
                                smcl, sthf, sthu, snow_surft,                  & 
                                ftl_1, ftl_surft, fqw_1, fqw_surft, le_surft,  &
                                tstar_surft, dtstar_surft,                     &
                                surf_ht_flux_land, surf_htf_surft,             &
                                ecan_surft, esoil_surft, ei_surft,             &
                                radnet_surft, gs, gs_surft,                    &
                                t1p5m_surft, q1p5m_surft, melt_surft,          &
                                NPP_gb, NPP_pft, NPP_acc_pft, GPP_gb, GPP_pft, &
                                resp_s, resp_s_tot, resp_p, resp_p_pft,        &
                                g_leaf_pft, &!RESP_S_TILE, !Kathy-ed as diag 
                                progs, work)                                                          
  
!subrs called 
USE cable_implicit_driv_mod,   ONLY: cable_implicit_driver
USE cable_implicit_unpack_mod, ONLY: implicit_unpack
USE prognostic_bank_mod_cbl,   ONLY: cable_reinstate_prognostics
USE prognostic_bank_mod_cbl,   ONLY: cable_store_prognostics
USE init_active_tile_mask_mod,  ONLY: init_active_tile_mask_cbl

! data: TYPE definitions of passed asarguments
USE progs_cbl_vars_mod,      ONLY: progs_cbl_vars_type ! CABLE's extra progs
USE work_vars_mod_cbl,       ONLY: work_vars_type      ! CABLE's types etc 
                                                       ! kept thru timestep
USE prognostic_bank_mod_cbl, ONLY: ProgBank 

! data: Scalars
USE cable_common_module,    ONLY: knode_gl
USE grid_constants_mod_cbl, ONLY: nrb, nsnl, mp
USE grid_constants_mod_cbl, ONLY: nsCs          ! # soil carbon stores
USE grid_constants_mod_cbl, ONLY: nvCs          ! # vegetation carbon stores 

USE cable_common_module, ONLY : cable_runtime !jhan:have to sort this out for JAC
!jhan: Leave for reference
!# if defined(UM_JULES)
!  ! CABLE prognostics declared at top_level
!  USE atm_fields_real_mod, ONLY : soil_temp_cable, soil_moist_cable,           &
!                                  soil_froz_frac_cable, snow_dpth_cable,       & 
!                                  snow_mass_cable, snow_temp_cable,            &
!                                  snow_rho_cable, 
!                                  snow_age_cable, snow_flg_cable,              &
!                                  aquifer_moist_cable
!  ! CASA prognostics declared at top_level
!  USE atm_fields_real_mod, ONLY : C_pool_casa, N_pool_casa, P_pool_casa,       &
!                                  SOIL_ORDER_casa, N_DEP_casa, N_FIX_casa,     &
!                                  P_DUST_casa, P_weath_casa, LAI_casa,         &
!                                  PHENPHASE_casa, RSP_W_PFT_ACC
!  !UM: time info 
!  USE model_time_mod, ONLY:    target_end_stepim, i_day

!  USE atmos_physics2_alloc_mod, ONLY : resp_s_tile

# if defined(UM_JULES)
USE model_time_mod, ONLY: i_day_number
#else
USE model_time_mod, ONLY: timesteps_in_day
#endif

IMPLICIT NONE

!___ re-decl input args (INTENT?)
INTEGER :: row_length, rows, land_pts, nsurft, npft, sm_levels
INTEGER :: dim_cs1
INTEGER :: cycleno, numcycles
REAL    :: timestep
INTEGER :: timestep_number
INTEGER :: surft_pts(nsurft) ! # land points per tile
INTEGER :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER :: land_index(land_pts)  ! tangled cell index of land_pt

REAL :: tile_frac(land_pts, nsurft)
REAL :: Fland(land_pts)
REAL :: smvcst(land_pts,sm_levels)
                       ! IN Volumetric saturation point
      
REAL :: ls_rain(row_length,rows)   !forcing%rain precip: large scale
REAL :: ls_snow(row_length,rows)   !forcing%snow precip: large scale
REAL :: conv_rain(row_length,rows) !forcing%rain precip: convective 
REAL :: conv_snow(row_length,rows) !forcing%snow precip: convective 
  
!prog: canopy water stores
REAL :: canopy_gb(land_pts)      !prog: aggregate over tiles 
REAL :: canopy_surft(land_pts, nsurft) !prog: per tile

REAL :: tl_1(row_length,rows)    !forcing: temp
REAL :: qw_1(row_length,rows)    !forcing: humidity 
REAL :: dtl1_1(row_length,rows)  !forcing: increment to temp
REAL :: dqw1_1(row_length,rows)  !forcing: increment to humidity 
!Ticket #132 needs ctctq1
REAL :: ctctq1(row_length,rows)  !forcing: temp/humidity increment

REAL :: ftl_1(row_length,rows)
REAL :: fqw_1(row_length,rows)

REAL, INTENT(OUT) :: ftl_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: fqw_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: le_surft ( land_pts, nsurft )       ! latent heat flux
    
!prog: UM soil quantities,  per land_pt
REAL, INTENT(OUT) :: t_soil(land_pts,sm_levels)
REAL, INTENT(OUT) :: smcl(land_pts,sm_levels)
REAL, INTENT(OUT) :: sthf(land_pts,sm_levels)
REAL, INTENT(OUT) :: sthu(land_pts,sm_levels)
REAL, INTENT(OUT) :: snow_surft( land_pts, nsurft ) ! Lying snow [kg/m2]

REAL, INTENT(OUT) :: tstar_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: dtstar_surft( land_pts, nsurft )

REAL, INTENT(OUT) :: surf_htf_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: surf_ht_flux_land(row_length,rows)
REAL, INTENT(OUT) :: ecan_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: esoil_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: ei_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: radnet_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: gs_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: t1p5m_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: q1p5m_surft( land_pts, nsurft )
REAL, INTENT(OUT) :: melt_surft( land_pts, nsurft )

REAL, INTENT(OUT) :: gs( land_pts )  ! conductance for use in dust scheme

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

TYPE(progs_cbl_vars_type), INTENT(IN OUT)  :: progs
TYPE(work_vars_type), INTENT(IN OUT)       :: work

!___ local vars
!Instantiate prognostic storage (bank)
TYPE (ProgBank) :: pb

# if !defined(UM_JULES)
INTEGER :: i_day_number
#endif
REAL, ALLOCATABLE :: dtrad(:)
LOGICAL, ALLOCATABLE :: l_tile_pts(:,:)
LOGICAL, SAVE :: first_call = .true.
LOGICAL, SAVE :: zero_points_warning= .true.
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_implicit_main"

IF( land_pts  ==  0 ) THEN
  IF( zero_points_warning ) THEN
    WRITE(6,*) "Reached CABLE ", subr_name,                                      &
               " even though zero land_points on processor ", knode_gl 
  END IF
  zero_points_warning = .FALSE. 
  RETURN
END IF

cable_runtime%um          = .TRUE.
cable_runtime%um_implicit = .TRUE.

! this done every call (maybe we hould pass this through work%)
!------------------------------------------------------------------------------
! Determine the number of active tiles
mp = SUM(surft_pts)

IF( .NOT. ALLOCATED(l_tile_pts) ) ALLOCATE( l_tile_pts(land_pts, nsurft) ) 

! Define mapping mask. i.e. l_tile_pts =TRUE (active) , where tile_frac > 0
CALL init_active_tile_mask_cbl(l_tile_pts, land_pts, nsurft, tile_frac )
!------------------------------------------------------------------------------

IF( .NOT. ALLOCATED(dtrad) ) ALLOCATE( dtrad(mp) ) 

# if !defined(UM_JULES)
i_day_number = FLOOR( REAL(timestep_number) / REAL(timesteps_in_day) )
#endif

!Due to ENDGAME, CABLE(any LSM) is called twice on implicit step.   
CALL cable_store_prognostics( pb, mp, sm_levels, nsCs,nvCs, work%ssnow, &
                              work%canopy, work%bgc )

!jhan:check - in CM2 but left out here? do this closer to evaluation
IF (first_call ) THEN
  T1P5M_surft = 0.0 
  Q1P5M_surft = 0.0 
ENDIF

CALL cable_implicit_driver( cycleno, numcycles, i_day_number,                   &
                            timestep, timestep_number,                          & 
                            row_length, rows, land_pts, nsurft, npft,          &
                            sm_levels, dim_cs1, mp, nrb,                       &
                            land_index, surft_pts, surft_index, l_tile_pts,    &
                            ls_rain, conv_rain, ls_snow, conv_snow,            &
                            tl_1, qw_1, ftl_1, fqw_1, dtl1_1, dqw1_1, ctctq1,  &
                            work%rad, work%met, work%rough, work%canopy,       &
                            work%veg, work%soil, work%ssnow, work%bal,         &
                            work%air, work%bgc, work%sum_flux )
                            !CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,               &
                            !GLAI, PHENPHASE)
                          
IF (cycleno .NE. numcycles ) THEN
  CALL cable_reinstate_prognostics( pb, work%ssnow, work%canopy, work%bgc )
ENDIF

CALL implicit_unpack( cycleno, row_length, rows, land_pts, nsurft, npft,       &
                      sm_levels, dim_cs1, timestep, mp, nsnl, land_index,      &
                      surft_pts, surft_index, tile_frac, l_tile_pts, smvcst,   &
                      t_soil, smcl, sthf, sthu, snow_surft,                    &
                      FTL_1, FTL_surft, FQW_1, FQW_surft,  le_surft,            &
                      TSTAR_surft, SURF_HT_FLUX_LAND, ECAN_surft,              &
                      ESOIL_surft, EI_surft, RADNET_surft,                     &
                      canopy_surft, gs, gs_surft, t1p5m_surft, q1p5m_surft,    &
                      canopy_gb, fland, melt_surft,                            &
                      NPP_gb, npp_pft, GPP_gb, GPP_pft, resp_s, resp_s_tot,    &
                      RESP_P, resp_p_pft, g_leaf_pft,                          &
                      NPP_acc_pft, surf_htf_surft,                             &
                      dtstar_surft, progs, work,                               & 
                      work%rad, work%met, work%rough, work%canopy,             &
                      work%veg, work%soil, work%ssnow, work%bal,               &
                      work%air, work%bgc, work%sum_flux )
 
first_call = .FALSE.        
cable_runtime%um_implicit = .FALSE.

RETURN

END SUBROUTINE cable_implicit_main
  
END MODULE cable_implicit_main_mod


