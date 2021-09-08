
MODULE cable_implicit_driv_mod
  
CONTAINS

SUBROUTINE cable_implicit_driver( i_day_number, cycleno, &! num_cycles 
                          row_length, rows, land_pts, ntiles, npft,           &
                          sm_levels, dim_cs1, dim_cs2, Fland,                 &
                          ls_rain, con_rain, ls_snow, conv_snow,              &
                          dtl_1,dqw_1, ctctq1, tsoil, tsoil_tile, smcl,       &
                          smcl_tile, &!SMGW_TILE, 
                          timestep_width, smvcst, sthf,                       &
                          sthf_tile, sthu, snow_tile, snow_rho1l, isnow_flg3l,&
                          snow_depth3l, snow_mass3l, snow_rho3l, snow_tmp3l,  &
                          ftl_1, ftl_tile, fqw_1, fqw_tile, tstar_tile,       &
                          surf_ht_flux_land, ecan_tile, esoil_tile, ei_tile,  &
                          radnet_tile, snow_age, canopy_tile, gs, gs_tile,    &
                          t1p5m_tile, q1p5m_tile, canopy_gb, melt_tile,       &
                          !NPP, NPP_FT, GPP, GPP_FT, RESP_S,                    &
                          !RESP_S_TOT,  RESP_P, RESP_P_FT,          &
                          !G_LEAF, 
                          tl_1, qw_1, surf_htf_tile,                          &
                          !CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,                  &
                          !GLAI, PHENPHASE, NPP_FT_ACC, RESP_W_FT_ACC, 
                          dtrad,                                              &
air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
soil_cbl)
  !subrs called 
USE cbl_model_driver_mod, ONLY: cbl_model_driver
USE cable_um_init_subrs_mod, ONLY: um2cable_rr
!H!  USE casa_cable, only : bgcdriver, sumcflux
  
!data
USE cable_other_constants_mod, ONLY: z0surf_min
USE cable_phys_constants_mod, ONLY : ccapp => capp
USE cable_def_types_mod, ONLY: mp, msn, ncs,ncp, nrb
USE cable_um_tech_mod,   ONLY: um1, conv_rain_prevstep, conv_snow_prevstep
USE cable_common_module, ONLY: cable_runtime, cable_user, l_casacnp,          &
                                l_vcmaxFeedbk, knode_gl, ktau_gl, kend_gl

USE cable_air_type_mod,       ONLY: air_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_params_mod,         ONLY: soil_parameter_type

IMPLICIT NONE
!___ re-decl input args
TYPE (air_type),       INTENT(INOUT) :: air_cbl
TYPE (bgc_pool_type),  INTENT(INOUT) :: bgc_cbl
TYPE (canopy_type),    INTENT(INOUT) :: canopy_cbl
TYPE (met_type),       INTENT(INOUT) :: met_cbl
TYPE (balances_type),  INTENT(INOUT) :: bal_cbl
TYPE (radiation_type), INTENT(INOUT) :: rad_cbl
TYPE (roughness_type), INTENT(INOUT) :: rough_cbl
TYPE (soil_snow_type), INTENT(INOUT) :: ssnow_cbl
TYPE (sum_flux_type),  INTENT(INOUT) :: sum_flux_cbl
!H!TYPE (climate_type) :: climate

TYPE (soil_parameter_type), INTENT(INOUT)   :: soil_cbl
TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg_cbl

!necessary as arg checking is enforce in modular structure that now present 
! - HOWEVER *NB*  this POP is not initialized anywhere
  
INTEGER :: cycleno
INTEGER :: row_length,rows, land_pts, ntiles, npft, sm_levels
INTEGER :: dim_cs1, dim_cs2 
  
REAL,  DIMENSION(land_pts) :: fland       ! IN Land fraction on land tiles
   
REAL, DIMENSION(row_length,rows) ::                                           &
  ls_rain,  & ! IN Large scale rain
  ls_snow,  & ! IN Large scale snow
  con_rain, & ! IN Convective rain
  conv_snow,& ! IN Convective snow
  tl_1,     & !
  qw_1,     & !
  dtl_1,    & ! IN Level 1 increment to T field 
  dqw_1,    & ! IN Level 1 increment to q field 
  ctctq1,   & ! IN information needed for increment to T an q field   
  surf_ht_flux_land, & ! Net DW heat flux at surface (W/m2).
  fqw_1,   & ! Moisture flux between layers. (kg/m^2/sec).
             !--- FQW(,1) is total water flux from surface, 'E'.
  ftl_1      !  FTL(,K) =net turbulent sensible heat flux into layer K
             !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)

REAL :: timestep_width

REAL, DIMENSION(land_pts) ::                                                  &
  gs,      &  ! OUT "Stomatal" conductance to
  smvcst

REAL, DIMENSION(land_pts,ntiles) ::                                           &
  surf_htf_tile,                                                              &
  ftl_tile, fqw_tile, & ! Surface FTL, FQL for land tiles
  le_tile, melt_tile, & ! latent heat flux, melting
  gs_tile,            & ! tiled stomatatal conductance
  radnet_tile,        & ! Surface net radiation on tiles (W/m2)
  tot_alb,     & ! total albedo
  ei_tile,     & ! OUT EI for land tiles.
  ecan_tile,   & ! OUT ECAN for snow-free land tiles
  esoil_tile     ! evapotranspiration from soil moisture store (kg/m2/s) 

REAL, DIMENSION(land_pts,sm_levels) ::                                        &
  smcl,       & ! UM aggregated soil moisture
  sthf,       & ! UM aggregated soil frozen fraction
  sthu,       & ! UM aggregated soil unfrozen fraction
  tsoil         ! UM aggregated soil temp.

! (tiled) soil prognostics: as above 
REAL, DIMENSION(land_pts,ntiles,sm_levels) ::                                 &
  smcl_tile, & !
  sthf_tile, & !
  sthu_tile, & !
  tsoil_tile

REAL, DIMENSION(land_pts,ntiles) ::                                           &
  smgw_tile

!___flag for 3 layer snow pack
INTEGER :: isnow_flg3l(land_pts,ntiles)
   
!___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
REAL, DIMENSION(land_pts,ntiles,3) ::                                         &
  snow_depth3l,  & ! 
  snow_mass3l,   & !
  snow_rho3l,    & !
  snow_tmp3l,    & !
  snow_cond        !
  
REAL, DIMENSION(land_pts,ntiles) ::                                           &
  resp_p_ft,                                                                  &
  g_leaf,                                                                     &
  npp_ft,                                                                     &
  gpp_ft      

REAL :: dtrad(mp)          !change in Trad over the time step

REAL, DIMENSION(land_pts) ::                                                  &
  snow_grd,    & !
  canopy_gb,   & !
  resp_p,      & !
  npp,         & !
  gpp
      
REAL, DIMENSION( land_pts,ntiles ) ::                                         &
  snow_tile,                                                                  &
  snow_rho1l,    &  ! Mean snow density
  snow_age,                                                                   &
  canopy_tile,                                                                &
  t1p5m_tile,                                                                 &
  q1p5m_tile,                                                                 &
  tstar_tile,                                                                 &
  resp_s_tile,                                                                &
  transp_tile

REAL ::                                                                       &
  resp_s(land_pts,dim_cs1),                                                   &
  resp_s_tot(dim_cs2)    

REAL, DIMENSION(land_pts,ntiles,10) ::                                        &
  cpool_tile,                                                                 &
   npool_tile     
REAL, DIMENSION(land_pts,ntiles,12) ::                                        &
  ppool_tile
REAL, DIMENSION(land_pts,ntiles) ::                                           &
  glai,                                                                       &
  phenphase

REAL, DIMENSION(land_pts,ntiles) ::                                           &
  npp_ft_acc,                                                                 &
  resp_w_ft_acc

INTEGER ::                                                                    &
  ktauday,      & ! day counter for CASA-CNP
  i_day_number, & ! day of year (1:365) counter for CASA-CNP
  idoy            ! day of year (1:365) counter for CASA-CNP
INTEGER, SAVE ::                                                              &
  kstart = 1

REAL, DIMENSION(mp) ::                                                        &
  dtlc,                                                                       &
  dqwc

!Ticket 132 - need ctctq1, incoming values of ftl_1 and fqw_1 on tiles
REAL, DIMENSION(mp) ::                                                        &
  ctctq1c,      &  ! UM boundary layer coefficient
  ftl1c,        &  ! grid box averaged FTL
  fqw1c            ! gird box averaged FQW
  
REAL, DIMENSION(land_pts) ::                                                  &
  lying_snow,    & ! OUT Gridbox snowmass (kg/m2)        
  sub_surf_roff, & !
  surf_roff,     & !
  tot_tfall        !

!___ local vars
  
!This is a quick fix. These can be organised through namelists
LOGICAL :: spinup = .FALSE., spinconv = .FALSE.,                              &
           dump_read = .FALSE., dump_write = .FALSE.
INTEGER :: loy = 365, lalloc = 0
  
!___ 1st call in RUN (!=ktau_gl -see below) 
LOGICAL, SAVE :: first_cable_call = .TRUE.
REAL, ALLOCATABLE:: fwork(:,:,:)

!Prog Bank copies all prognostic and other variables whose
!values need to be retain from UM timestep to UM timestep.
!NOTE that canopy%cansto is a prognostic variable but is handled
!differently through the canopy%oldcansto variable
TYPE ProgBank
     
  REAL, DIMENSION(:,:), ALLOCATABLE ::                                        &
    tsoil, smcl, sthf,                                                        &
    snow_depth,                                                               &
    snow_mass, snow_tmp, snow_rho ,                                           &
    sli_s, SLI_Tsoil, SLI_sconds, SLI_snowliq,                                &
    cplant, csoil !carbon variables
    
  REAL, DIMENSION(:), ALLOCATABLE ::                                          &
    snow_rho1l, snow_age, snow_tile,                                          &
    ocanopy,                                                                  &
    fes_cor,fhs_cor, osnowd,owetfac,otss,GWwb,tss0,                           &
    puddle, rtsoil, wblake, GWaq,                                             &
    SLI_h0, SLI_Tsurf
    
  INTEGER, DIMENSION(:), ALLOCATABLE ::                                       &
    snow_flg3l, SLI_nsnow
    
END TYPE ProgBank

INTEGER, PARAMETER :: cpb  = 2 !grab from UM, hardwired to ENDGAME norm 10.6
  
!Instantiate, NB:using dims=cpb 
TYPE (ProgBank), DIMENSION(cpb), SAVE :: pb
   
INTEGER :: ipb

! std template args 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_implicit_driver"
LOGICAL, PARAMETER :: explicit_path = .FALSE.
  
!-------- Unique subroutine body -----------
        
!Due to ENDGAME, CABLE(any LSM) is called twice on implicit step.   
CALL cable_store_prognostics()

IF (ipb == cpb) CALL cable_reinstate_prognostics()

dtlc = 0.0 ; dqwc = 0.0

!--- All these subrs do is pack a CABLE var with a UM var.
!-------------------------------------------------------------------
!--- UM met forcing vars needed by CABLE which have UM dimensions
!---(rowlength,rows)[_rr], which is no good to CABLE. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!--- if the land point is/has an active tile
!--- generic format:
!--- um2cable_rr( UM var, default value for snow tile, CABLE var, mask )
!--- where mask tells um2cable_rr whether or not to use default value 
!--- for snow tile 
!-------------------------------------------------------------------
  
CALL um2cable_rr( (ls_rain + con_rain) * timestep_width, met_cbl%precip)

CALL um2cable_rr( (ls_snow + conv_snow) * timestep_width, met_cbl%precip_sn)
      
CALL um2cable_rr( tl_1, met_cbl%tk)
      
CALL um2cable_rr( qw_1, met_cbl%qv)
      
CALL um2cable_rr( dtl_1, dtlc)
      
CALL um2cable_rr( dqw_1, dqwc)
      
!--- conv_rain(snow)_prevstep are added to precip. in explicit call
CALL um2cable_rr( (con_rain) * timestep_width, conv_rain_prevstep)
  
CALL um2cable_rr( (CONV_snow) * timestep_width, conv_snow_prevstep)

!#define UM_JULES - coupling to UM *NOT* JULES standalone
!Ticket #132 implementation --------------------------------------------
!dtlc, dqwc found on tiles above - these are the corrected dtlc and dqwc 
IF (cable_user%l_revised_coupling) THEN
  CALL um2cable_rr( ctctq1, ctctq1c)
  CALL um2cable_rr( ftl_1, ftl1c)
  CALL um2cable_rr( fqw_1, fqw1c) 
  dtlc = dtlc - ctctq1c * ftl1c / ccapp  !NB FTL_1 is in W/m2 hence / CAPP
  dqwc = dqwc - ctctq1c * fqw1c
END IF
!-----------------------------------------------------------------------

met_cbl%precip   =  met_cbl%precip + met_cbl%precip_sn
met_cbl%tk = met_cbl%tk
!hacks to match offline !@ Loobos
!met_cbl%tk = met_cbl%tk + dtlc
met_cbl%qv = met_cbl%qv
!hacks to match offline !@ Loobos
!met_cbl%qv = met_cbl%qv + dqwc
met_cbl%tvair = met_cbl%tk
met_cbl%tvrad = met_cbl%tk

canopy_cbl%cansto = canopy_cbl%oldcansto

CALL cbl_model_driver( explicit_path, mp, nrb, land_pts, npft, ktau_gl,       &
          timestep_width,      &
          air_cbl, bgc_cbl, canopy_cbl, met_cbl, bal_cbl, &
          rad_cbl, rough_cbl, soil_cbl, ssnow_cbl, sum_flux_cbl, veg_cbl, z0surf_min, &
          !H!shouuld already work from here LAI_pft, HGT_pft )
          veg_cbl%vlai, veg_cbl%hc, met_cbl%doy, canopy_cbl%vlaiw )

    ! Integrate wb_lake over the river timestep.
    ! Used to scale river flow within ACCESS
    ! Zeroed each river step in subroutine cable_lakesriver and on restarts.
    !  ssnow_wb_lake in kg/m^2
    !C!if (ipb == cpb) THEN
    !C!  ssnow%totwblake = ssnow%totwblake + ssnow%wb_lake/river_step
    !C!end if
 
!Jun 2018 - change in Trad over time step
dtrad = rad_cbl%trad - rad_cbl%otrad

! Lestevens - temporary ?
ktauday = INT(24.0 * 3600.0 / timestep_width)
idoy = i_day_number
  
!Jan 2018: Only call carbon cycle prognostics updates on the last call to 
!cable_implicit per atmospheric time step
IF (ipb == cpb) THEN
  !Call CASA-CNP
  !H!if (l_casacnp) & 
  !H!  CALL bgcdriver(ktau_gl,kstart,kend_gl,timestep,met,ssnow,canopy,veg,soil, &
  !H!                 climate,casabiome,casapool,casaflux,casamet,casabal,phen, &
  !H!                 pop, spinConv,spinup, ktauday, idoy,loy, dump_read,   &
  !H!                 dump_write, LALLOC)
  !H! have to comment out as we dont havecasaflux yet
      !H!CALL sumcflux(ktau_gl,kstart,kend_gl,TIMESTEP,bgc,canopy,soil,ssnow,      &
      !H!              sum_flux,veg,met,casaflux,l_vcmaxFeedbk)
END IF

! Only call carbon cycle prognostics updates on the last call to 
! cable_implicit per atmospheric time step
! Call CASA-CNP collect pools
!H!if (ipb==cpb .AND. l_casacnp) & 
!H!  CALL casa_poolout_unpk(casapool,casaflux,casamet,casabal,phen,  &
!H!                        CPOOL_TILE,NPOOL_TILE,PPOOL_TILE, &
!H!                        GLAI,PHENPHASE)

RETURN

CONTAINS

SUBROUTINE cable_store_prognostics()
IMPLICIT NONE

!cpb = cable% um% numcycles
IF ( .NOT. ALLOCATED(pb(1) %tsoil) ) THEN
 
  DO ipb = 1, cpb 
      
    ALLOCATE( pb(ipb)%tsoil(mp,sm_levels) )
    ALLOCATE( pb(ipb)%smcl(mp,sm_levels) )
    ALLOCATE( pb(ipb)%sthf(mp,sm_levels) )
    ALLOCATE( pb(ipb)%snow_depth(mp,3) )
    ALLOCATE( pb(ipb)%snow_mass(mp,3) )
    ALLOCATE( pb(ipb)%snow_tmp(mp,3) )
    ALLOCATE( pb(ipb)%snow_rho(mp,3) )
    ALLOCATE( pb(ipb)%snow_rho1l(mp) )
    ALLOCATE( pb(ipb)%snow_age(mp) )
    ALLOCATE( pb(ipb)%snow_flg3l(mp) )
    ALLOCATE( pb(ipb)%snow_tile(mp) )
    ALLOCATE( pb(ipb)%ocanopy(mp) )
    !Jan 2018 new PB variables
    ALLOCATE( pb(ipb)%fes_cor(mp) )
    ALLOCATE( pb(ipb)%puddle(mp) )
    ALLOCATE( pb(ipb)%owetfac(mp) )
    ALLOCATE( pb(ipb)%rtsoil(mp) )
    ALLOCATE( pb(ipb)%wblake(mp) )
    !carbon variables - may not be needed unless CASA
    ALLOCATE( pb(ipb)%cplant(mp,ncp) )
    ALLOCATE( pb(ipb)%csoil(mp,ncs) )
    !GW model variables no need to test always has a value
    !so do not introduce issues with restarting a GW run 
    ALLOCATE (pb(ipb)%GWaq(mp) )
    !otss - July 2018
    ALLOCATE(pb(ipb)%tss0(mp) )
      
    !SLI variables - Jhan please check the second dimension
    IF (cable_user%soil_struc == 'sli') THEN
      ALLOCATE(pb(ipb)%SLI_nsnow(mp) )
      ALLOCATE(pb(ipb)%sli_s(mp,sm_levels) )
      ALLOCATE(pb(ipb)%SLI_Tsoil(mp,sm_levels) )
      ALLOCATE(pb(ipb)%SLI_sconds(mp,3) )
      ALLOCATE(pb(ipb)%SLI_h0(mp) )
      ALLOCATE(pb(ipb)%SLI_Tsurf(mp) )
      ALLOCATE(pb(ipb)%SLI_snowliq(mp,3) )
    END IF
      
  END DO
  
END IF !.NOT. allocated   

ipb = cycleno

pb(ipb)%tsoil     = ssnow_cbl%tgg
pb(ipb)%smcl      = ssnow_cbl%wb
pb(ipb)%sthf      = ssnow_cbl%wbice
pb(ipb)%snow_depth= ssnow_cbl%sdepth
pb(ipb)%snow_mass = ssnow_cbl%smass
pb(ipb)%snow_tmp  = ssnow_cbl%tggsn
pb(ipb)%snow_rho  = ssnow_cbl%ssdn
pb(ipb)%snow_rho1l= ssnow_cbl%ssdnn
pb(ipb)%snow_age  = ssnow_cbl%snage
pb(ipb)%snow_flg3l= ssnow_cbl%isflag
pb(ipb)%snow_tile = ssnow_cbl%snowd
pb(ipb)%ocanopy   = canopy_cbl%oldcansto
!Jan 2018 new PB variables
pb(ipb)%fes_cor   = canopy_cbl%fes_cor
pb(ipb)%puddle    = ssnow_cbl%pudsto
pb(ipb)%rtsoil    = ssnow_cbl%rtsoil !?needed
pb(ipb)%owetfac   = ssnow_cbl%owetfac
pb(ipb)%wblake    = ssnow_cbl%wb_lake 
!carbon variables - may not be needed unless CASA
pb(ipb)%cplant = bgc_cbl%cplant
pb(ipb)%csoil = bgc_cbl%csoil
!GW model variables
pb(ipb)%GWaq   = ssnow_cbl%GWwb
!otss added July 2018
pb(ipb)%tss0 = ssnow_cbl%tss
  
!SLI variables
IF (cable_user%soil_struc == 'sli') THEN
  pb(ipb)%SLI_nsnow   = ssnow_cbl%nsnow
  pb(ipb)%sli_s       = ssnow_cbl%s
  pb(ipb)%Tsoil       = ssnow_cbl%Tsoil
  pb(ipb)%SLI_sconds  = ssnow_cbl%sconds
  pb(ipb)%SLI_h0      = ssnow_cbl%h0
  pb(ipb)%SLI_Tsurf   = ssnow_cbl%Tsurface
  pb(ipb)%SLI_snowliq = ssnow_cbl%snowliq
END IF


END SUBROUTINE cable_store_prognostics

SUBROUTINE cable_reinstate_prognostics()
IMPLICIT NONE

ssnow_cbl%tgg     = pb(1)%tsoil
ssnow_cbl%wb      = pb(1)%smcl
ssnow_cbl%wbice   = pb(1)%sthf
ssnow_cbl%sdepth  = pb(1)%snow_depth
ssnow_cbl%smass   = pb(1)%snow_mass
ssnow_cbl%tggsn   = pb(1)%snow_tmp
ssnow_cbl%ssdn    = pb(1)%snow_rho
ssnow_cbl%ssdnn   = pb(1)%snow_rho1l
ssnow_cbl%snage   = pb(1)%snow_age
ssnow_cbl%isflag  = pb(1)%snow_flg3l
ssnow_cbl%snowd   = pb(1)%snow_tile
canopy_cbl%oldcansto = pb(1)%ocanopy
!Jan 2018 new PB variables 
canopy_cbl%fes_cor = pb(1)%fes_cor
ssnow_cbl%pudsto  = pb(1)%puddle
ssnow_cbl%rtsoil  = pb(1)%rtsoil  !?needed
ssnow_cbl%owetfac = pb(1)%owetfac
ssnow_cbl%wb_lake = pb(1)%wblake  
!carbon variables - may not be needed unless CASA
bgc_cbl%cplant = pb(1)%cplant
bgc_cbl%csoil = pb(1)%csoil
!GW model variables
ssnow_cbl%GWwb = pb(1)%GWaq
!%tss added July 2018
ssnow_cbl%tss = pb(1)%tss0
  
!SLI variables
IF (cable_user%soil_struc == 'sli') THEN
  ssnow_cbl%nsnow    = pb(1)%SLI_nsnow
  ssnow_cbl%s        = pb(1)%sli_s
  ssnow_cbl%Tsoil    = pb(1)%SLI_Tsoil
  ssnow_cbl%sconds   = pb(1)%SLI_sconds
  ssnow_cbl%h0       = pb(1)%SLI_h0
  ssnow_cbl%Tsurface = pb(1)%SLI_Tsurf
  ssnow_cbl%snowliq  = pb(1)%SLI_snowliq
END IF

END SUBROUTINE cable_reinstate_prognostics
 
END SUBROUTINE cable_implicit_driver

END MODULE cable_implicit_driv_mod
