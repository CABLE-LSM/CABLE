MODULE cable_explicit_main_mod

CONTAINS

SUBROUTINE cable_explicit_main(                                                &
      ! IN: UM/JULES model/grid parameters, fields, mappings
      mype, timestep_len, timestep_number, row_length, rows, land_pts,         &
      nsurft, npft, sm_levels, dzsoil, land_index, surft_pts, surft_index,     &
      cosine_zenith_angle, latitude, longitude, Fland, tile_frac,              &
      
      ! IN: soil parameters !1 is only allowable index in UM
      bexp, hcon, satcon, sathh, smvcst, smvcwt, smvccl, albsoil,              &
      
      ! IN: Met forcing: 
      lw_down, sw_surft, ls_rain, ls_snow,                                     &
      tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy_tile,                 &
      ! This an outlier IN here. INOUT @ implicit. (was)OUT at extras
      ! I think we are dealing with it OK now but confusion could be removed  
      snow_tile,                                                               & 
      
      ! IN: canopy height, LAI seasonally presecribed, potentially prognostic 
      ! IN: CO2 mass mixing ratio  
      canht_pft, lai_pft, CO2_MMR,                                             & 
      
      ! TYPEs passed from top_level to maintain scope, access to UM STASH 
      ! IN: tiled soil/snow prognostics - IN here. INOUT @ implicit  
      ! INOUT: Carries fields needed by CABLE b/n pathways (rad, explicit etc) 
      !        Currently carrying CABLE TYPEs (canopy%, rad% etc).   
      ! IN: pars carries vegin/soilin - potentially redundant 
      progs, work, pars,                                                       & 
      
      ! OUT: UM fields UNPACKed from CABLE (@ explicit) 
      ftl_tile, fqw_tile, tstar_tile, dtstar_surft,                            &
      u_s, u_s_std_tile, cd_tile, ch_tile,                                     &
      radnet_tile, fraca, resfs, resft, z0h_tile, z0m_tile,                    &
      recip_l_mo_tile, epot_tile, npp_pft_acc, resp_w_pft_acc )
  
! subrs 
USE cable_explicit_driv_mod,    ONLY: cable_explicit_driver
USE cable_expl_unpack_mod,      ONLY: cable_expl_unpack
USE init_active_tile_mask_mod,  ONLY: init_active_tile_mask_cbl

! data: TYPE definitions of passed asarguments
USE progs_cbl_vars_mod, ONLY: progs_cbl_vars_type ! CABLE requires extra progs
USE work_vars_mod_cbl,  ONLY: work_vars_type      ! and some kept thru timestep
USE params_io_mod_cbl,  ONLY: params_io_data_type
USE params_io_mod_cbl,  ONLY: params_io_type

! data: Scalars
USE grid_constants_mod_cbl, ONLY: nrb, nrs, mp
USE cable_common_module,    ONLY: knode_gl, ktau_gl, kwidth_gl, cable_runtime, &
                                  cable_user, redistrb, satuparam, wiltparam,  &
                                  l_casacnp_cd => l_casacnp
USE cable_model_env_opts_mod, ONLY: icycle  ! 0=No CASA- [1=C,2=CN,3=CNP]
USE cable_model_env_opts_mod, ONLY: l_casacnp 
USE casadimension,            ONLY: icycle_cd =>  icycle   
!Leave for reference
!!  USE atm_fields_real_mod, ONLY : soil_temp_cable, soil_moist_cable, etc       &
!!                                  C_pool_casa, N_pool_casa, P_pool_casa,       &
!!                                  SOIL_ORDER_casa, N_DEP_casa, N_FIX_casa,     &
!!                                  P_DUST_casa, P_weath_casa, LAI_casa,         &
!!                                  PHENPHASE_casa, NPP_PFT_ACC, RSP_W_PFT_ACC,  &
!!                                  aquifer_moist_cable,aquifer_thickness_cable, &
!!                                  slope_avg_cable,slope_std_cable,&
!!                                  visc_sublayer_depth,aquifer_perm_cable,&
!!                                  aquifer_draindens_cable

IMPLICIT NONE
 
INTEGER, INTENT(IN) :: mype              ! # processor   
REAL,    INTENT(IN) :: timestep_len      ! # seconds  (cucurrently 1200)
INTEGER, INTENT(IN) :: timestep_number   ! # timestep (cucurrently 3 per hr) 
INTEGER, INTENT(IN) :: row_length        ! # columns in spatial grid
INTEGER, INTENT(IN) :: rows              ! # rows in spatial grid
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
INTEGER, INTENT(IN) :: npft              ! # plant functional types
INTEGER, INTENT(IN) :: sm_levels         ! # soil layers 
REAL,    INTENT(IN) :: dzsoil(sm_levels) ! soil layer thicknesses 

INTEGER, INTENT(IN) :: land_index(land_pts)          ! land point indices 
                                                     ! recipe back to (i,j) cell  
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! tile points indices
                                                     ! recipe back to land_index

REAL, INTENT(IN) :: canht_pft(land_pts, npft)   ! canopy height (seasonal)
REAL, INTENT(IN) :: lai_pft(land_pts, npft)     ! LAI           (seasonal)
REAL, INTENT(IN) :: fland(land_pts)             ! land fraction (<1 for coastal) 
REAL, INTENT(IN) :: co2_mmr                     ! prescribed MMR
REAL, INTENT(IN) :: tile_frac(land_pts, nsurft) ! tile fraction

REAL, INTENT(IN) :: cosine_zenith_angle(row_length,rows)
REAL, INTENT(IN) :: latitude (row_length,rows)
REAL, INTENT(IN) :: longitude(row_length,rows)

! soil parameters
REAL, INTENT(IN) :: bexp (land_pts, sm_levels)  ! parameter b in Campbell eqn 
REAL, INTENT(IN) :: sathh(land_pts, sm_levels)  !
REAL, INTENT(IN) :: smvcst(land_pts, sm_levels) ! 
REAL, INTENT(IN) :: smvcwt(land_pts, sm_levels) !
REAL, INTENT(IN) :: smvccl(land_pts, sm_levels) !
REAL, INTENT(IN) :: hcon(land_pts)              ! Soil thermal conductivity (W/m/K).
REAL, INTENT(IN) :: albsoil(land_pts)           ! bare soil albedo
REAL, INTENT(IN) :: satcon(land_pts, sm_levels) ! hydraulic conductivity
                                                !  @ saturation [mm/s]

! "forcing"
REAL, INTENT(IN) :: lw_down(row_length,rows)
REAL, INTENT(IN) :: ls_rain(row_length,rows)
REAL, INTENT(IN) :: ls_snow(row_length,rows)
REAL, INTENT(IN) :: sw_surft(land_pts, nsurft)
REAL, INTENT(IN) :: tl_1(row_length,rows)
REAL, INTENT(IN) :: qw_1(row_length,rows)
REAL, INTENT(IN) :: vshr_land(row_length,rows)
REAL, INTENT(IN) :: pstar(row_length,rows)
REAL, INTENT(IN) :: z1_tq(row_length,rows)
REAL, INTENT(IN) :: z1_uv(row_length,rows)

! prognostics   
REAL, INTENT(IN) :: canopy_tile(land_pts, nsurft)
REAL, INTENT(IN) :: snow_tile(land_pts, nsurft) ! Lying snow [kg/m2]

! TYPEs passed from top_level to maintain scope, access to UM STASH 
! IN: tiled soil/snow prognostics - IN here. INOUT @ implicit  
! INOUT: Carries fields needed by CABLE b/n pathways (rad, explicit etc) 
!        Currently carrying CABLE TYPEs (canopy%, rad% etc).   
! IN: pars carries vegin/soilin - potentially redundant 
TYPE(progs_cbl_vars_type), INTENT(IN)    :: progs
TYPE(params_io_data_type), INTENT(IN)    :: pars
TYPE(work_vars_type),      INTENT(INOUT) :: work

! OUT: UM fields UNPACKed from CABLE (@ explicit) 
REAL, INTENT(OUT) :: ftl_tile(land_pts,nsurft)    ! surface FTL for land tiles     
REAL, INTENT(OUT) :: fqw_tile(land_pts,nsurft)    ! surface FQW for land tiles     
REAL, INTENT(OUT) :: tstar_tile(land_pts,nsurft)  ! radiative surf. temperature
REAL, INTENT(OUT) :: dtstar_surft(land_pts,nsurft) ! 
REAL, INTENT(OUT) :: u_s(row_length,rows)         ! friction velocity (m/s)
REAL, INTENT(OUT) :: u_s_std_tile(land_pts,nsurft)!
REAL, INTENT(OUT) :: cd_tile(land_pts,nsurft)     ! Drag coefficient
REAL, INTENT(OUT) :: ch_tile(land_pts,nsurft)     ! Transfer coefficient
REAL, INTENT(OUT) :: radnet_tile(land_pts,nsurft) ! Surface net radiation
REAL, INTENT(OUT) :: z0h_tile(land_pts,nsurft)    ! roughness
REAL, INTENT(OUT) :: z0m_tile(land_pts,nsurft)    ! roughness
REAL, INTENT(OUT) :: epot_tile(land_pts,nsurft)   !
REAL, INTENT(OUT) :: recip_l_mo_tile(land_pts,nsurft) ! Reciprocal:Monin-Obukhov
                                                      ! length for tiles (m^-1)
REAL, INTENT(OUT) :: fraca(land_pts,nsurft)       ! Fraction - surface moisture
REAL, INTENT(OUT) :: RESFS(land_pts,nsurft)
                      ! Combined soil, stomatal & aerodynamic resistance
                      ! factor for fraction (1-FRACA) of snow-free land tiles
REAL, INTENT(OUT) :: RESFT(land_pts,nsurft)
                      ! Total resistance factor.
                      ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                      ! 1 for snow.    
REAL, INTENT(IN) :: npp_pft_acc(land_pts,npft)
REAL, INTENT(IN) :: resp_w_pft_acc (land_pts,npft)
     
!___ local vars, may be passed as args downstream
LOGICAL :: cbl_standalone = .FALSE. !needs to be set from namelist  
LOGICAL :: jls_standalone = .FALSE. !needs to be set from namelist 
LOGICAL :: jls_radiation  = .FALSE. !needs to be set from n amelist 

INTEGER :: isnow_flg_cable(land_pts, nsurft)
REAL    :: radians_degrees
REAL    :: latitude_deg(row_length,rows)
REAL    :: longitude_deg(row_length,rows)
REAL    :: sw_down_ij(row_length,rows,nrs)
REAL    :: sw_down_TOT(row_length,rows)
REAL    :: sw_down_DIR(row_length,rows)
REAL    :: sw_down_VIS(row_length,rows)
REAL    :: sw_down_NIR(row_length,rows)
REAL    :: beamFrac_VIS(row_length,rows)
REAL    :: beamFrac_NIR(row_length,rows)
REAL    :: beamFrac_TOT(row_length,rows)

LOGICAL, ALLOCATABLE :: l_tile_pts(:,:)
    
INTEGER       :: i,j,k,l,n
LOGICAL, SAVE :: zero_points_warning = .TRUE.
 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_explicit_main"
    
IF( land_pts  ==  0 ) THEN
  IF( zero_points_warning ) THEN
    WRITE(6,*) "Reached CABLE ", subr_name,                                      &
               " even though zero land_points on processor ", mype 
  END IF
  zero_points_warning = .FALSE.
  RETURN
END IF

!--- Set up some cable-globals -------------------------------------------------
cable_runtime%um          = .TRUE.
cable_runtime%um_explicit = .TRUE.

! this done every call (maybe we hould pass this through work%)
!------------------------------------------------------------------------------
! Determine the number of active tiles
mp = SUM(surft_pts)

IF( .NOT. ALLOCATED(l_tile_pts)  ) ALLOCATE( l_tile_pts(land_pts, nsurft) ) 

! Define mapping mask. i.e. l_tile_pts =TRUE (active) , where tile_frac > 0
CALL init_active_tile_mask_cbl( l_tile_pts, land_pts, nsurft, tile_frac )
!-------------------------------------------------------------------------------

!!extracted from ap/um/rose-app.conf au-aa809@2729
![namelist:cable]
cable_user%diag_soil_resp='ON'
cable_user%fwsoil_switch='Haverd2013'
cable_user%gs_switch='medlyn'
cable_user%gw_model=.false.
cable_user%l_rev_corr=.true.
cable_user%l_revised_coupling=.true.
cable_user%or_evap=.false.
!cable_user%soil_thermal_fix=.true.
cable_user%soil_thermal_fix=.false.! fudge - worked to dt=4
cable_user%ssnow_potev='HDM'
redistrb=.false.
satuparam=0.8
wiltparam=0.5
! set icycle/lcasacnp seen thru-out model from namelist read version  
icycle_cd    = icycle   
l_casacnp_cd = l_casacnp

! initialize processor number, timestep len
knode_gl  = mype 
kwidth_gl = INT(timestep_len)
ktau_gl   = timestep_number

!--- Convert lat/long to degrees
radians_degrees = 180.0 / ( 4.0*atan(1.0) ) ! 180 / PI
latitude_deg  = latitude * radians_degrees
longitude_deg  = longitude * radians_degrees

isnow_flg_cable = INT(progs%ThreeLayerSnowFlag_CABLE)

!--- Fix SW for CABLE  ----------------------------------------------------------------------------
sw_down_ij(:,:,:)  = 0.0
sw_down_TOT(:,:)   = 0.0
sw_down_DIR(:,:)   = 0.0
sw_down_VIS(:,:)   = 0.0
sw_down_NIR(:,:)   = 0.0
beamFrac_VIS(:,:)  = 0.0
beamFrac_NIR(:,:)  = 0.0
beamFrac_TOT(:,:)  = 0.0

IF(jls_standalone) THEN

  DO n = 1, nsurft
    ! loop over number of points per tile
    DO k = 1, surft_pts(n)
      l = surft_index(k, n)
      j = (land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      sw_down_VIS(i, j) = sw_surft(l,n) / 2.0
    END DO
  END DO
  sw_down_NIR(:,:) = sw_down_VIS(:,:)

ELSE
  
  ! in all cases zenith angle needs to be applied
  DO n = 1, nrs
    sw_down_ij(:,:,n) = work%sw_down_ij(:,:,n) * cosine_zenith_angle(:,:)
  END DO
  
  ! SUM over ALL components of sw_down_ij
  sw_down_TOT(:,:) = sw_down_ij(:,:,1) + sw_down_ij(:,:,2) +                        &
                sw_down_ij(:,:,3) + sw_down_ij(:,:,4)  
  
  ! SUM DIRect components of sw_down_ij(in VIS & NIR )
  sw_down_DIR(:,:) = sw_down_ij(:,:,1) + sw_down_ij(:,:,3)
  
  ! SUM VIS components of sw_down_ij(incl DIR & DIF )
  sw_down_VIS(:,:) = sw_down_ij(:,:,1) + sw_down_ij(:,:,2)  
  
  ! SUM NIR components of sw_down_ij(incl DIR & DIF )
  sw_down_NIR(:,:) = sw_down_ij(:,:,3) + sw_down_ij(:,:,4)

  ! beam(DIR) fraction in VISible spectrum 
  beamFrac_VIS(:,:) = sw_down_ij(:,:,1) / MAX( 0.1, sw_down_VIS(:,:) )  
  
  ! beam(DIR) fraction in NIR spectrum 
  beamFrac_NIR(:,:) = sw_down_ij(:,:,3) / MAX( 0.1, sw_down_NIR(:,:) )  
                                 
  ! beam(DIR) fraction for all solar 
  beamFrac_TOT(:,:) = sw_down_DIR(:,:) / MAX( 0.1, sw_down_TOT(:,:) )  
 
ENDIF

!----------------------------------------------------------------------------
!--- CALL _driver to run specific and necessary components of CABLE with IN -
!--- args PACKED to force CABLE
!-------------------------------------------------------------------------------
CALL cable_explicit_driver(                                                    & 
  ! IN: UM/JULES/CABLE model/grid parameters, fields, mappings
  mype, row_length, rows, land_pts, nsurft, npft, sm_levels, dzsoil,           &
  timestep_len, timestep_number, mp, nrb, land_index, surft_pts, surft_index,  &
  l_tile_pts, latitude_deg, longitude_deg, cosine_zenith_angle, Fland,         &
  tile_frac,                                                                   &

  ! IN: soil parameters !1 is only allowable index in UM
  bexp, hcon, satcon, sathh, smvcst, smvcwt, smvccl, albsoil,                  &
  
  ! IN: SW forcing: manipulated for CABLE 
  sw_down_VIS, sw_down_NIR, beamFrac_VIS, beamFrac_NIR, beamFrac_TOT,          &
      
  ! IN: Met forcing: 
  lw_down, ls_rain, ls_snow,                                                   &
  tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy_tile,                     &
  ! This an outlier IN here. INOUT @ implicit. (was)OUT at extras
  ! I think we are dealing with it OK now but confusion could be removed  
  snow_tile,                                                                   &
 
  ! IN: canopy height, LAI seasonally presecribed, potentially prognostic 
  ! IN: CO2 mass mixing ratio  
  canht_pft, lai_pft, CO2_MMR,                                                 &

  ! IN: carries vegin/soilin - potentially redundant  
  pars,                                                                        &
      
  ! IN: tiled soil/snow prognostics - IN here. INOUT @ implicit  
  progs%soiltemp_CABLE, progs%soilmoisture_CABLE, progs%FrozenSoilFrac_CABLE,  &
  isnow_flg_cable, progs%SnowDepth_CABLE, progs%SnowMass_CABLE,                &
  progs%SnowTemp_CABLE, progs%SnowDensity_CABLE, progs%snowage_CABLE,          &
  progs%snowosurft, progs%OneLyrSnowDensity_CABLE,                             &
      
  ! INOUT: CABLE TYPEs roughly grouped fields per module 
  work%rad, work%met, work%rough, work%canopy, work%veg, work%soil,            &
  work%ssnow, work%bal, work%air, work%bgc, work%sum_flux,                     &
  
  !IN: persistent veg%iveg, soil%isoilm are initialized on first rad/alb call     
  work%veg%iveg, work%soil%isoilm,                                             &
  !OUT: currently being passed back to UM in veg%hc, veg%vlai 
  work%veg%hc, work%veg%vlai,                                                  &
  
  !IN: currently being passed from prev radiation call through work%
  ! jhan:quirky, snow (in turn reduced LAI due to snow) can evolve through a 
  ! constant rad dt. However reducedLAIdue2snow used ubiquitously as trigger 
  ! Further, snow does NOT evolve in explicit AND reducedLAIdue2snow absent 
  ! in implicit
  work%reducedLAIdue2snow,                                                     & 
      
          !GW
          !visc_sublayer_depth, smgw_tile, slope_avg, slope_std,
          !dz_gw, perm_gw, drain_gw,                           
          !casa progs
          !CPOOL_TILE, NPOOL_TILE, PPOOL_TILE, SOIL_ORDER, NIDEP,
          !NIFIX, PWEA, PDUST, GLAI, PHENPHASE,
      
  !IN: if not passed a dangling argument would ensue 
  npp_pft_acc, resp_w_pft_acc )
    
!----------------------------------------------------------------------------
!--- CALL _unpack to unpack variables from CABLE back to UM format to return
!----------------------------------------------------------------------------
call cable_expl_unpack( & 
  ! IN: UM/JULES/CABLE model/grid parameters, fields, mappings
  row_length, rows, land_pts, nsurft, npft, mp, land_index, surft_pts,         &
  surft_index, l_tile_pts, fland, tile_frac, latitude, longitude,              &
  
  !OUT: UM fields to be updated
  ftl_tile, fqw_tile, tstar_tile, dtstar_surft , u_s, u_s_std_tile, cd_tile,   &
  ch_tile, radnet_tile, fraca, resfs, resft, z0h_tile, z0m_tile,               &
  recip_l_mo_tile, epot_tile,                                                  &
 
  !IN: UM fields to be updated FROM these CABLE fields
  work%canopy%fh, work%canopy%fes, work%canopy%fev, work%canopy%us,            &
  work%canopy%cdtq,work%canopy%fwet, work%canopy%wetfac_cs,                    & 
  work%canopy%rnet, work%canopy%zetar, work%canopy%epot, work%rad%trad,        &
  work%rad%otrad, work%rad%transd, work%rough%z0m, work%rough%zref_tq,         &
  
  !IN: UM fields used in derivation of fields to be updated
  work%ssnow%snowd, work%ssnow%cls, work%air%rlam, work%air%rho, work%met%ua )

cable_runtime%um_explicit = .FALSE.

RETURN

END SUBROUTINE cable_explicit_main
  
END MODULE cable_explicit_main_mod

