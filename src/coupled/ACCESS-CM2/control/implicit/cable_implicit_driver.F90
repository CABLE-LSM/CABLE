!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
! 
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located 
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Updates CABLE variables (as altered by first pass through boundary 
!          layer and convection scheme), calls cbm, passes CABLE variables back 
!          to UM. 'Implicit' is the second call to cbm in each UM timestep.
!
! Called from: UM codecable_implicit_main 
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
! ==============================================================================

module cable_implicit_driv_mod
  
contains

subroutine cable_implicit_driver( i_day_number, cycleno, &! num_cycles 
                          row_length, rows, land_pts, ntiles, npft,            &
                          sm_levels, dim_cs1, dim_cs2, Fland,                  &
                          LS_RAIN, CON_RAIN, LS_SNOW, CONV_SNOW,               &
                          DTL_1,DQW_1, ctctq1, TSOIL, TSOIL_TILE, SMCL,        &
                          SMCL_TILE, SMGW_TILE, timestep, SMVCST, STHF,        &
                          STHF_TILE, STHU, snow_tile, SNOW_RHO1L, ISNOW_FLG3L, &
                          SNOW_DEPTH3L, SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L,   &
                          FTL_1, FTL_TILE, FQW_1, FQW_TILE, TSTAR_TILE,        &
                          SURF_HT_FLUX_LAND, ECAN_TILE, ESOIL_TILE, EI_TILE,   &
                          RADNET_TILE, SNOW_AGE, CANOPY_TILE, GS, GS_TILE,     &
                          T1P5M_TILE, Q1P5M_TILE, CANOPY_GB, MELT_TILE,        &
                          NPP, NPP_FT, GPP, GPP_FT, RESP_S,                    &
                          RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT,          &
                          G_LEAF, TL_1, QW_1, SURF_HTF_TILE,                   &
                          CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,                  &
                          GLAI, PHENPHASE, NPP_FT_ACC, RESP_W_FT_ACC, DTRAD )

  !subrs called 
  USE cable_um_init_subrs_mod, ONLY : um2cable_rr
  USE cable_cbm_module,    ONLY : cbm
  
  USE cable_def_types_mod, ONLY : mp, msn, ncs,ncp, climate_type
USE cable_phys_constants_mod,  ONLY: TFRZ, CAPP 
  USE cable_um_tech_mod,   ONLY : um1, conv_rain_prevstep, conv_snow_prevstep,&
                                 air, bgc, canopy, met, bal, rad, rough, soil,&
                                 ssnow, sum_flux, veg, basic_diag
  USE cable_common_module, ONLY : cable_runtime, cable_user, l_casacnp,       &
                                  l_vcmaxFeedbk, knode_gl, ktau_gl, kend_gl
  
  USE casavariable
  USE phenvariable
  USE casa_types_mod
  USE casa_um_inout_mod
  use POP_TYPES, only : pop_type
  USE river_inputs_mod,   ONLY: river_step
  

  implicit none

  !___ re-decl input args
  TYPE (climate_type)  :: climate     ! climate variables
  !necessary as arg checking is enforce in modular structure that now present 
  ! - HOWEVER *NB*  this POP is not initialized anywhere
  TYPE(POP_TYPE) :: POP 
  
  integer :: cycleno
  integer :: row_length,rows, land_pts, ntiles, npft, sm_levels
  integer :: dim_cs1, dim_cs2 
  
  REAL,  DIMENSION(land_pts) :: fland       ! IN Land fraction on land tiles
   
  REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                 &
    LS_RAIN,  & ! IN Large scale rain
    LS_SNOW,  & ! IN Large scale snow
    CON_RAIN, & ! IN Convective rain
    CONV_SNOW,& ! IN Convective snow
    TL_1,     & !
    QW_1,     & !
    DTL_1,    & ! IN Level 1 increment to T field 
    DQW_1,    & ! IN Level 1 increment to q field 
    ctctq1,   & ! IN information needed for increment to T an q field   
    SURF_HT_FLUX_LAND, & ! Net DW heat flux at surface (W/m2).
    FQW_1,   & ! Moisture flux between layers. (kg/m^2/sec).
               !--- FQW(,1) is total water flux from surface, 'E'.
    FTL_1      !  FTL(,K) =net turbulent sensible heat flux into layer K
               !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)

  REAL :: timestep

  REAL, DIMENSION(land_pts) ::                                            &
    GS,      &  ! OUT "Stomatal" conductance to
    SMVCST

  REAL, DIMENSION(land_pts,ntiles) ::                              &
    SURF_HTF_TILE, &
    FTL_TILE, FQW_TILE, & ! Surface FTL, FQL for land tiles
    LE_TILE, MELT_TILE, & ! latent heat flux, melting
    GS_TILE,            & ! tiled stomatatal conductance
    RADNET_TILE,        & ! Surface net radiation on tiles (W/m2)
    TOT_ALB,     & ! total albedo
    EI_TILE,     & ! OUT EI for land tiles.
    ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
    ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

  REAL, dimension(land_pts,sm_levels) ::                           &
    SMCL,       & ! UM aggregated soil moisture
    STHF,       & ! UM aggregated soil frozen fraction
    STHU,       & ! UM aggregated soil unfrozen fraction
    TSOIL         ! UM aggregated soil temp.

  ! (tiled) soil prognostics: as above 
  REAL, dimension(land_pts,ntiles,sm_levels) ::                &
    SMCL_TILE, & !
    STHF_TILE, & !
    STHU_TILE, & !
    TSOIL_TILE

  REAL, dimension(land_pts,ntiles) ::                &
    SMGW_TILE

  !___flag for 3 layer snow pack
  INTEGER :: ISNOW_FLG3L(LAND_PTS,NTILES)
   
  !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
  REAL, dimension(land_pts,ntiles,3) :: &
    SNOW_DEPTH3L,  & ! 
    SNOW_MASS3L,   & !
    SNOW_RHO3L,    & !
    SNOW_TMP3L,    & !
    SNOW_COND        !
  
  REAL, DIMENSION(land_pts,ntiles) ::                              &
    RESP_P_FT,     &
    G_LEAF,        &
    NPP_FT,     &
    GPP_FT      

  REAL, DIMENSION(land_pts) ::                                         &
    SNOW_GRD,    & !
    CANOPY_GB,   & !
    RESP_P,      & !
    NPP,         & !
    GPP

  REAL, DIMENSION(mp) ::  DTRAD          !change in Trad over the time step
      
  REAL, DIMENSION( land_pts,ntiles ) ::                               &
    SNOW_TILE,     &
    SNOW_RHO1L,    &  ! Mean snow density
    SNOW_AGE,    &
    CANOPY_TILE,   &
    T1P5M_TILE,    &
    Q1P5M_TILE,    &
    TSTAR_TILE,    &
    RESP_S_TILE,   & 
    TRANSP_TILE

  REAL ::                                                                     &
    RESP_S(LAND_PTS,DIM_CS1),     &
    RESP_S_TOT(DIM_CS2)    

  REAL, DIMENSION(LAND_PTS,NTILES,10) ::                              &
    CPOOL_TILE, &
     NPOOL_TILE     
  REAL, DIMENSION(LAND_PTS,NTILES,12) ::                              &
    PPOOL_TILE
  REAL, DIMENSION(LAND_PTS,NTILES) ::                                 &
    GLAI, &
    PHENPHASE

  REAL, DIMENSION(LAND_PTS,NTILES) ::                                 &
    NPP_FT_ACC, &
    RESP_W_FT_ACC

  INTEGER ::     &
    ktauday,      & ! day counter for CASA-CNP
    i_day_number, & ! day of year (1:365) counter for CASA-CNP
    idoy            ! day of year (1:365) counter for CASA-CNP
  INTEGER, SAVE :: &
    kstart = 1

  REAL, DIMENSION(mp) ::                                                      & 
    dtlc, & 
    dqwc

  !Ticket 132 - need ctctq1, incoming values of ftl_1 and fqw_1 on tiles
  REAL, DIMENSION(mp) ::                                                      & 
    ctctq1c,      &  ! UM boundary layer coefficient
    ftl1c,        &  ! grid box averaged FTL
    fqw1c            ! gird box averaged FQW
  
  REAL, DIMENSION(LAND_PTS) ::                               &
    LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
    SUB_SURF_ROFF, & !
    SURF_ROFF,     & !
    TOT_TFALL        !

  !___ local vars
  

  !This is a quick fix. These can be organised through namelists
  logical :: spinup=.false., spinconv=.false.,                   &
             dump_read=.false., dump_write=.false.
  integer :: loy=365, lalloc=0
  
  !___ 1st call in RUN (!=ktau_gl -see below) 
  LOGICAL, SAVE :: first_cable_call = .TRUE.
  REAL, ALLOCATABLE:: fwork(:,:,:)

  !Prog Bank copies all prognostic and other variables whose
  !values need to be retain from UM timestep to UM timestep.
  !NOTE that canopy%cansto is a prognostic variable but is handled
  !differently through the canopy%oldcansto variable
  type ProgBank
     
    real, dimension(:,:), allocatable ::                                       &
      TSOIL, SMCL, STHF,                                                       & 
      SNOW_DEPTH,                                                              &
      SNOW_MASS, SNOW_TMP, SNOW_RHO ,                                          &
      SLI_S, SLI_Tsoil, SLI_sconds, SLI_snowliq,                               &
      cplant, csoil !carbon variables
    
    real, dimension(:), allocatable ::                                         &
      SNOW_RHO1L, SNOW_AGE, SNOW_TILE,                                         &
      OCANOPY,                                                                 &
      fes_cor,fhs_cor, osnowd,owetfac,otss,GWwb,tss0,                        &
      puddle, rtsoil, wblake, GWaq,                                            &
      SLI_h0, SLI_Tsurf
    
    integer, dimension(:), allocatable ::                                      &
      SNOW_FLG3L, SLI_nsnow
    
  End type ProgBank

  integer, parameter :: cpb =2 !grab from UM, hardwired to ENDGAME norm 10.6
  
  !Instantiate, NB:using dims=cpb 
  type (ProgBank), dimension(cpb), save :: PB
   
  integer :: ipb

  ! std template args 
  character(len=*), parameter :: subr_name = "cable_implicit_driver"

  !-------- Unique subroutine body -----------
        
  !Due to ENDGAME, CABLE(any LSM) is called twice on implicit step.   
  call cable_store_prognostics()

  if (ipb == cpb) call cable_reinstate_prognostics()

      dtlc = 0. ; dqwc = 0.

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
  ! Leave here as example
  CALL um2cable_rr( (LS_RAIN+CON_RAIN)*TIMESTEP, met%precip)

  CALL um2cable_rr( (LS_SNOW+CONV_SNOW)*TIMESTEP, met%precip_sn)
      
  CALL um2cable_rr( TL_1, met%tk)
      
  CALL um2cable_rr( QW_1, met%qv)
      
  CALL um2cable_rr( dtl_1, dtlc)
      
  CALL um2cable_rr( dqw_1, dqwc)
      
  !--- conv_rain(snow)_prevstep are added to precip. in explicit call
  CALL um2cable_rr( (CON_RAIN)*TIMESTEP, conv_rain_prevstep)
  
  CALL um2cable_rr( (CONV_snow)*TIMESTEP, conv_snow_prevstep)

  !Ticket #132 implementation --------------------------------------------
  !dtlc, dqwc found on tiles above - these are the corrected dtlc and dqwc 
  if (cable_user%l_revised_coupling) then
    CALL um2cable_rr( ctctq1, ctctq1c)
    CALL um2cable_rr( FTL_1, ftl1c)
    CALL um2cable_rr( FQW_1, fqw1c) 
    dtlc = dtlc - ctctq1c*ftl1c/CAPP  !NB FTL_1 is in W/m2 hence / CAPP
    dqwc = dqwc - ctctq1c*fqw1c
  endif
  !-----------------------------------------------------------------------

  met%precip   =  met%precip + met%precip_sn
  met%tk = met%tk + dtlc
  met%qv = met%qv + dqwc
  met%tvair = met%tk
  met%tvrad = met%tk

  canopy%cansto = canopy%oldcansto

  CALL cbm( ktau_gl,timestep, air, bgc, canopy, met, bal,                             &
            rad, rough, soil, ssnow, sum_flux, veg, climate )

      ! Integrate wb_lake over the river timestep.
      ! Used to scale river flow within ACCESS
      ! Zeroed each river step in subroutine cable_lakesriver and on restarts.
      !  ssnow_wb_lake in kg/m^2
      if (ipb == cpb) THEN
        ssnow%totwblake = ssnow%totwblake + ssnow%wb_lake/river_step
      end if
 
  !Jun 2018 - change in Trad over time step
  DTRAD = rad%trad - rad%otrad

  ! Lestevens - temporary ?
  ktauday = int(24.0*3600.0/TIMESTEP)
  idoy=i_day_number
  
  !Jan 2018: Only call carbon cycle prognostics updates on the last call to 
  !cable_implicit per atmospheric time step
  if (ipb==cpb) then
    !Call CASA-CNP
    !CM2!if (l_casacnp) & 
    !CM2!  CALL bgcdriver(ktau_gl,kstart,kend_gl,timestep,met,ssnow,canopy,veg,soil, &
    !CM2!                 climate,casabiome,casapool,casaflux,casamet,casabal,phen, &
    !CM2!                 pop, spinConv,spinup, ktauday, idoy,loy, dump_read,   &
    !CM2!                 dump_write, LALLOC)

    !CM2!CALL sumcflux(ktau_gl,kstart,kend_gl,TIMESTEP,bgc,canopy,soil,ssnow,      &
    !CM2!              sum_flux,veg,met,casaflux,l_vcmaxFeedbk)
  endif

  ! Only call carbon cycle prognostics updates on the last call to 
  ! cable_implicit per atmospheric time step
  ! Call CASA-CNP collect pools
  if (ipb==cpb .AND. l_casacnp) & 
    CALL casa_poolout_unpk(casapool,casaflux,casamet,casabal,phen,  &
                          CPOOL_TILE,NPOOL_TILE,PPOOL_TILE, &
                          GLAI,PHENPHASE)

  !-------- End Unique subroutine body -----------

!End Testing puroses:

return

CONTAINS

subroutine cable_store_prognostics()
  implicit none

  !cpb = cable% um% numcycles
  if (.NOT. allocated(PB(1) %tsoil) ) then
 
    do ipb = 1, cpb 
      
      allocate( PB(ipb)%TSOIL(mp,sm_levels) )
      allocate( PB(ipb)%SMCL(mp,sm_levels) )
      allocate( PB(ipb)%STHF(mp,sm_levels) )
      allocate( PB(ipb)%snow_depth(mp,3) )
      allocate( PB(ipb)%snow_mass(mp,3) )
      allocate( PB(ipb)%snow_tmp(mp,3) )
      allocate( PB(ipb)%snow_rho(mp,3) )
      allocate( PB(ipb)%snow_rho1l(mp) )
      allocate( PB(ipb)%snow_age(mp) )
      allocate( PB(ipb)%snow_flg3l(mp) )
      allocate( PB(ipb)%snow_tile(mp) )
      allocate( PB(ipb)%ocanopy(mp) )
      !Jan 2018 new PB variables
      allocate( PB(ipb)%fes_cor(mp) )
      allocate( PB(ipb)%puddle(mp) )
      allocate( PB(ipb)%owetfac(mp) )
      allocate( PB(ipb)%rtsoil(mp) )
      allocate( PB(ipb)%wblake(mp) )
      !carbon variables - may not be needed unless CASA
      allocate( PB(ipb)%cplant(mp,ncp) )
      allocate( PB(ipb)%csoil(mp,ncs) )
      !GW model variables no need to test always has a value
      !so do not introduce issues with restarting a GW run 
      allocate (PB(ipb)%GWaq(mp) )
      !otss - July 2018
      allocate(PB(ipb)%tss0(mp) )
      
      !SLI variables - Jhan please check the second dimension
      if (cable_user%soil_struc=='sli') then
        allocate(PB(ipb)%SLI_nsnow(mp) )
        allocate(PB(ipb)%SLI_S(mp,sm_levels) )
        allocate(PB(ipb)%SLI_Tsoil(mp,sm_levels) )
        allocate(PB(ipb)%SLI_sconds(mp,3) )
        allocate(PB(ipb)%SLI_h0(mp) )
        allocate(PB(ipb)%SLI_Tsurf(mp) )
        allocate(PB(ipb)%SLI_snowliq(mp,3) )
      endif
      
    enddo
  
  endif !.NOT. allocated   

  ipb = cycleno

  PB(ipb)%tsoil     = ssnow%tgg
  PB(ipb)%smcl      = ssnow%wb
  PB(ipb)%sthf      = ssnow%wbice
  PB(ipb)%snow_depth= ssnow%sdepth
  PB(ipb)%snow_mass = ssnow%smass
  PB(ipb)%snow_tmp  = ssnow%tggsn
  PB(ipb)%snow_rho  = ssnow%ssdn
  PB(ipb)%snow_rho1l= ssnow%ssdnn
  PB(ipb)%snow_age  = ssnow%snage
  PB(ipb)%snow_flg3l= ssnow%isflag
  PB(ipb)%snow_tile = ssnow%snowd
  PB(ipb)%ocanopy   = canopy%oldcansto
  !Jan 2018 new PB variables
  PB(ipb)%fes_cor   = canopy%fes_cor
  PB(ipb)%puddle    = ssnow%pudsto
  PB(ipb)%rtsoil    = ssnow%rtsoil !?needed
  PB(ipb)%owetfac   = ssnow%owetfac
  PB(ipb)%wblake    = ssnow%wb_lake 
  !carbon variables - may not be needed unless CASA
  PB(ipb)%cplant = bgc%cplant
  PB(ipb)%csoil = bgc%csoil
  !GW model variables
  PB(ipb)%GWaq   = ssnow%GWwb
  !otss added July 2018
  PB(ipb)%tss0 = ssnow%tss
  
  !SLI variables
  if (cable_user%soil_struc=='sli') then
    PB(ipb)%SLI_nsnow   = ssnow%nsnow
    PB(ipb)%SLI_S       = ssnow%S
    PB(ipb)%Tsoil       = ssnow%Tsoil
    PB(ipb)%SLI_sconds  = ssnow%sconds
    PB(ipb)%SLI_h0      = ssnow%h0
    PB(ipb)%SLI_Tsurf   = ssnow%Tsurface
    PB(ipb)%SLI_snowliq = ssnow%snowliq
  endif


End SUBROUTINE cable_store_prognostics

SUBROUTINE cable_reinstate_prognostics()
  implicit none

  ssnow%tgg     = PB(1)%tsoil
  ssnow%wb      = PB(1)%smcl
  ssnow%wbice   = PB(1)%sthf
  ssnow%sdepth  = PB(1)%snow_depth
  ssnow%smass   = PB(1)%snow_mass
  ssnow%tggsn   = PB(1)%snow_tmp
  ssnow%ssdn    = PB(1)%snow_rho
  ssnow%ssdnn   = PB(1)%snow_rho1l
  ssnow%snage   = PB(1)%snow_age
  ssnow%isflag  = PB(1)%snow_flg3l
  ssnow%snowd   = PB(1)%snow_tile
  canopy%oldcansto = PB(1)%ocanopy
  !Jan 2018 new PB variables 
  canopy%fes_cor = PB(1)%fes_cor
  ssnow%pudsto  = PB(1)%puddle
  ssnow%rtsoil  = PB(1)%rtsoil  !?needed
  ssnow%owetfac = PB(1)%owetfac
  ssnow%wb_lake = PB(1)%wblake  
  !carbon variables - may not be needed unless CASA
  bgc%cplant = PB(1)%cplant
  bgc%csoil = PB(1)%csoil
  !GW model variables
  ssnow%GWwb = PB(1)%GWaq
  !%tss added July 2018
  ssnow%tss = PB(1)%tss0
  
  !SLI variables
  if (cable_user%soil_struc=='sli') then
    ssnow%nsnow    = PB(1)%SLI_nsnow
    ssnow%S        = PB(1)%SLI_S
    ssnow%Tsoil    = PB(1)%SLI_Tsoil
    ssnow%sconds   = PB(1)%SLI_sconds
    ssnow%h0       = PB(1)%SLI_h0
    ssnow%Tsurface = PB(1)%SLI_Tsurf
    ssnow%snowliq  = PB(1)%SLI_snowliq
  endif

End SUBROUTINE cable_reinstate_prognostics
 
END SUBROUTINE cable_implicit_driver


        

End module cable_implicit_driv_mod
