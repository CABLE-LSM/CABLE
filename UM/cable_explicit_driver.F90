!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Passes UM variables to CABLE, calls cbm, passes CABLE variables 
!          back to UM. 'Explicit' is the first of two routines that call cbm at 
!          different parts of the UM timestep.
!
! Called from: UM code sf_exch
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
!
! ==============================================================================


SUBROUTINE cable_explicit_driver( row_length, rows, land_pts, ntiles,npft,     &
                                  sm_levels, timestep, latitude, longitude,    &
                                  land_index, tile_frac,  tile_pts, tile_index,&
                                  bexp, hcon, satcon, sathh, smvcst,           &
                                  smvcwt,  smvccl, albsoil, snow_tile,         &
                                  snow_rho1l, snage_tile, isnow_flg3l,         &
                                  snow_rho3l, snow_cond, snow_depth3l,         &
                                  snow_tmp3l, snow_mass3l, sw_down, lw_down,   &
                                  cos_zenith_angle, surf_down_sw, ls_rain,     &
                                  ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,&
                                  z1_uv, rho_water, L_tile_pts, canopy_tile,   &
                                  Fland,                                       &
! rml 2/7/13 pass 3d co2 through to cable if required
                   CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,   &
                                  sthu_tile, smcl_tile,                        &
                                  sthf_tile, sthu, tsoil_tile, canht_ft,       &
                                  lai_ft, sin_theta_latitude, dzsoil,          &
                                  LAND_MASK, FTL_TILE_CAB, FTL_CAB, FTL_TILE,  &
                                  FQW_TILE, LE_TILE_CAB, LE_CAB, TSTAR_TILE,   &
                                  TSTAR_TILE_CAB, TSTAR_CAB, U_S, U_S_STD_TILE,&
                                  U_S_CAB, CH_CAB, CD_CAB, CD_TILE, CH_TILE,   &
                                  RADNET_TILE, FRACA, RESFS, RESFT, Z0H_TILE,  &
                                  Z0M_TILE, RECIP_L_MO_TILE, EPOT_TILE,        &
                                  CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,          &
                                  SOIL_ORDER, NIDEP, NIFIX, PWEA, PDUST,       &
                                  GLAI, PHENPHASE, PREV_YR_SFRAC,              &
                                  WOOD_HVEST_C,WOOD_HVEST_n,WOOD_HVEST_p,      &
                                  WOOD_FLUX_C,WOOD_FLUX_n,WOOD_FLUX_p,         &
                                  WRESP_C,WRESP_n,WRESP_p,THINNING,            &
                                  NPP_FT_ACC, RESP_W_FT_ACC, RESP_S_ACC,       &
                                  iday, endstep, timestep_number, mype )    
   
   !--- reads runtime and user switches and reports
   USE cable_um_tech_mod, ONLY : cable_um_runtime_vars, air, bgc, canopy,      &
                                 met, bal, rad, rough, soil, ssnow, sum_flux, veg 
   
   !--- vars common to CABLE declared 
   USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl,         &
                                   knode_gl, kwidth_gl, kend_gl,               &
                                   report_version_no,                          & 
                                   l_vcmaxFeedbk, l_laiFeedbk,l_luc
   
   !--- subr to (manage)interface UM data to CABLE
   USE cable_um_init_mod, ONLY : interface_UM_data
   
   !--- subr to call CABLE model
   USE cable_cbm_module, ONLY : cbm

   USE cable_def_types_mod, ONLY : mp

   !--- include subr called to write data for testing purposes 
   USE cable_diag_module
   USE casa_um_inout_mod
   USE casavariable
   USE casa_types_mod

  USE feedback_mod

   IMPLICIT NONE
 
 
  
   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM sf_exch() --------------------------------------------
   !-------------------------------------------------------------------------- 
   
   !___IN: UM dimensions, array indexes, flags
   INTEGER, INTENT(IN) ::                                                      & 
      row_length, rows, & ! UM grid resolution
      land_pts,         & ! # of land points being processed
      ntiles,           & ! # of tiles 
      npft,             & ! # of plant functional types
      sm_levels           ! # of soil layers 

   ! index of land points being processed
   INTEGER, INTENT(IN), DIMENSION(land_pts) :: land_index 

   ! # of land points on each tile
   INTEGER, INTENT(IN), DIMENSION(ntiles) :: tile_pts 
   
   INTEGER, INTENT(INOUT), DIMENSION(land_pts, ntiles) ::                         & 
      tile_index ,& ! index of tile points being processed
      isnow_flg3l   ! 3 layer snow flag

   !--- TRUE if land, F elsewhere.
   !jhan:rm land_mask
   LOGICAL,DIMENSION(row_length,rows) :: land_mask   

   !___UM parameters: water density, soil layer thicknesses 
   REAL, INTENT(IN) :: rho_water 
   REAL, INTENT(IN), DIMENSION(sm_levels) :: dzsoil

   !___UM soil/snow/radiation/met vars
   REAL, INTENT(IN), DIMENSION(land_pts) :: & 
      bexp,    & ! => parameter b in Campbell equation 
      hcon,    & ! Soil thermal conductivity (W/m/K).
      satcon,  & ! hydraulic conductivity @ saturation [mm/s]
      sathh,   &
      smvcst,  &
      smvcwt,  &
      smvccl,  &
      albsoil, &
      fland 
   
   REAL, INTENT(INOUT), DIMENSION(row_length,rows) :: &
      sw_down,          & 
      cos_zenith_angle
   
   REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             &
      latitude,   &
      longitude,  &
      lw_down,    &
      ls_rain,    &
      ls_snow,    &
      tl_1,       &
      qw_1,       &  
      vshr_land,  &
      pstar,      &
      z1_tq,      &
      z1_uv

   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles) ::                         &
      snow_tile

   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles) ::                            &
      tile_frac,  &    
      snow_rho1l, &
      snage_tile
   
   REAL, INTENT(IN), DIMENSION(row_length, rows, 4) ::                         &
      surf_down_sw 
   
   REAL, INTENT(IN), DIMENSION(land_pts, npft) ::                              &
      canht_ft, lai_ft 
   
   REAL, INTENT(INOUT),DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile
   
   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles,3) ::                       &
      snow_cond
   
   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles,3) ::                          &
      snow_rho3l,    &
      snow_depth3l,  &
      snow_mass3l,   &
      snow_tmp3l
   
   REAL, INTENT(IN), DIMENSION(land_pts, sm_levels) ::                         &
      sthu 
   
   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles, sm_levels) :: & 
      sthu_tile, &
      sthf_tile, &
      smcl_tile, &
      tsoil_tile
   
   REAL, INTENT(IN) :: co2_mmr
! rml 2/7/13 Extra atmospheric co2 variables
   LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
   INTEGER, INTENT(IN) ::                              &
      CO2_DIM_LEN                                      &
     ,CO2_DIM_ROW
   REAL, INTENT(IN) :: CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)  ! co2 mass mixing ratio

   !___true IF vegetation (tile) fraction is greater than 0
   LOGICAL, INTENT(INOUT), DIMENSION(land_pts, ntiles) :: L_tile_pts
  
   REAL :: sin_theta_latitude(row_length,rows) 
     
   !___return fluxes
   REAL, INTENT(OUT), DIMENSION(land_pts) ::   &
      FTL_CAB, &
      LE_CAB

   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
      FTL_TILE_CAB, &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      FQW_TILE,   &  ! Surface FQW for land tiles     
      LE_TILE_CAB

   !___return temp and roughness
   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
      TSTAR_TILE_CAB,   &
      TSTAR_TILE,       & 
      Z0H_TILE,         &
      Z0M_TILE

   REAL, INTENT(OUT), DIMENSION(land_pts) ::  TSTAR_CAB

   !___return friction velocities/drags/ etc
   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
      CD_TILE,    &     ! Drag coefficient
      CH_TILE,    &     ! Transfer coefficient for heat & moisture
      U_S_STD_TILE      ! Surface friction velocity

   REAL, INTENT(OUT), DIMENSION(row_length,rows)  :: &
      U_S               ! Surface friction velocity (m/s)
   
   REAL, INTENT(OUT), DIMENSION(land_pts) ::                  &
      CH_CAB,  &  ! Turbulent surface exchange
      CD_CAB,  &  ! Turbulent surface exchange
      U_S_CAB     ! Surface friction velocity (m/s)

   ! end step of experiment, this step, step width, processor num
   INTEGER, INTENT(IN) :: endstep, timestep_number, mype, iday
   REAL, INTENT(IN) ::  timestep     
   
   INTEGER:: itimestep
   INTEGER:: mtau
    
   !___return miscelaneous 
   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
      RADNET_TILE,   & ! Surface net radiation
      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                       ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,         & ! Total resistance factor.
                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                       ! 1 for snow.    
      FRACA,         & ! Fraction of surface moisture
      RECIP_L_MO_TILE,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
      EPOT_TILE

! Lestevens Sept2012 - CASA-CNP Pools
   REAL, INTENT(INOUT), DIMENSION(land_pts,ntiles,10) :: &
      CPOOL_TILE,    & ! Carbon Pools
      NPOOL_TILE       ! Nitrogen Pools

   REAL, INTENT(INOUT), DIMENSION(land_pts,ntiles,12) :: &
      PPOOL_TILE       ! Phosphorus Pools

   REAL, INTENT(INOUT), DIMENSION(land_pts) :: &
      SOIL_ORDER,    & ! Soil Order (1 to 12)
      NIDEP,         & ! Nitrogen Deposition
      NIFIX,         & ! Nitrogen Fixation
      PWEA,          & ! Phosphorus from Weathering
      PDUST            ! Phosphorus from Dust

   REAL, INTENT(INOUT), DIMENSION(land_pts,ntiles) :: &
      GLAI,         &    ! Leaf Area Index for Prognostics LAI
      PHENPHASE,    &    ! Phenology Phase for Casa-CNP
      PREV_YR_SFRAC,&    ! user_anc1, previous years surface fractions
      WOOD_FLUX_C,  &
      WOOD_FLUX_N,  &
      WOOD_FLUX_P,  &
      THINNING

   REAL, INTENT(INOUT), DIMENSION(land_pts,ntiles,3) :: &
      WOOD_HVEST_C,&
      WOOD_HVEST_N,&
      WOOD_HVEST_P,&
      WRESP_C,&
      WRESP_N,&
      WRESP_P
                                  
   REAL, INTENT(INOUT), DIMENSION(land_pts,ntiles) :: &
      NPP_FT_ACC,     &
      RESP_W_FT_ACC,  &
      RESP_S_ACC  
     
   !-------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM sf_exch() ----------------------------------------
   !-------------------------------------------------------------------------- 
   


   !___ declare local vars 
   
   !___ location of namelist file defining runtime vars
   CHARACTER(LEN=200), PARAMETER ::                                            & 
      runtime_vars_file = 'cable.nml'


   !___ 1st call in RUN (!=ktau_gl -see below) 
   LOGICAL, SAVE :: first_cable_call = .TRUE.
 

   

   !--- initialize cable_runtime% switches 
   IF(first_cable_call) THEN
      cable_runtime%um = .TRUE.
      write(6,*) ""
      write(6,*) "CABLE_log"
      CALL report_version_no(6) ! wriite revision number to stdout(6)
   ENDIF
   
write(6,*) "jhan:ESM1.5 test SB,BW 2"
      
   !--- basic info from global model passed to cable_common_module 
   !--- vars so don't need to be passed around, just USE _module
   ktau_gl = timestep_number     !timestep of EXPERIMENT not necesarily 
                                 !the same as timestep of particular RUN
   knode_gl = mype               !which processor am i on?
   itimestep = INT(timestep)     !realize for 'call cbm' pass
   kwidth_gl = itimestep         !width of timestep (secs)
   kend_gl = endstep             !timestep of EXPERIMENT not necesarily 

   !--- internal FLAGS def. specific call of CABLE from UM
   !--- from cable_common_module
   cable_runtime%um_explicit = .TRUE.

   
   !--- user FLAGS, variables etc def. in cable.nml is read on 
   !--- first time step of each run. these variables are read at 
   !--- runtime and for the most part do not require a model rebuild.
   IF(first_cable_call) THEN
      CALL cable_um_runtime_vars(runtime_vars_file) 
      first_cable_call = .FALSE.
   ENDIF      
   

  mtau = mod(ktau_gl,int(24.*3600./timestep))
  if (l_luc .and. iday==1 .and. mtau==1) then
   ! resdistr(frac,in,out) - Lestevens 10oct17
   call redistr_luc(PREV_YR_SFRAC,tsoil_tile   ,tsoil_tile)   ! ssnow%tgg
   call redistr_luc(PREV_YR_SFRAC,smcl_tile    ,smcl_tile)    ! ssnow%wb
   call redistr_luc(PREV_YR_SFRAC,sthf_tile    ,sthf_tile)    ! ssnow%wbice
   call redistr_luc(PREV_YR_SFRAC,snow_depth3L ,snow_depth3l) ! ssnow%sdepth
   call redistr_luc(PREV_YR_SFRAC,snow_mass3L  ,snow_mass3l)  ! ssnow%smass
   call redistr_luc(PREV_YR_SFRAC,snow_tmp3L   ,snow_tmp3l)   ! ssnow%tggsn
   call redistr_luc(PREV_YR_SFRAC,snow_rho3L   ,snow_rho3l)   ! ssnow%ssdn
   call redistr_luc(PREV_YR_SFRAC,snow_rho1l   ,snow_rho1l)   ! ssnow%ssdnn
   call redistr_luc(PREV_YR_SFRAC,snage_tile   ,snage_tile)   ! ssnow%snage
   !call redistr_luc_i(PREV_YR_SFRAC,isnow_flg3l  ,isnow_flg3l)  ! ssnow%isflag
   call redistr_luc(PREV_YR_SFRAC,snow_tile    ,snow_tile)    ! ssnow%snowd
   call redistr_luc(PREV_YR_SFRAC,snow_cond    ,snow_cond)    ! scond
   call redistr_luc(PREV_YR_SFRAC,canopy_tile  ,canopy_tile)  ! canopy%oldcansto
   call redistr_luc(PREV_YR_SFRAC,npp_ft_acc   ,npp_ft_acc)   ! frs
   call redistr_luc(PREV_YR_SFRAC,resp_w_ft_acc,resp_w_ft_acc)! frp
   call redistr_luc(PREV_YR_SFRAC,resp_s_acc,resp_s_acc)! npp
  endif

   !---------------------------------------------------------------------!
   !--- initialize CABLE using UM forcings etc. these args are passed ---!
   !--- down from UM.                                                 ---! 
   !---------------------------------------------------------------------!
   CALL interface_UM_data( row_length, rows, land_pts, ntiles, npft,           & 
                           sm_levels, itimestep, latitude, longitude,          &
                           land_index, tile_frac, tile_pts, tile_index,        &
                           bexp, hcon, satcon, sathh, smvcst, smvcwt,          &
                           smvccl, albsoil, snow_tile, snow_rho1l,             &
                           snage_tile, isnow_flg3l, snow_rho3l, snow_cond,     &
                           snow_depth3l, snow_tmp3l, snow_mass3l, sw_down,     &
                           lw_down, cos_zenith_angle, surf_down_sw, ls_rain,   &
                           ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,       &
                           z1_uv, rho_water, L_tile_pts, canopy_tile, Fland,   &
! rml 2/7/13 pass 3d co2 through to cable if required
                   CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,   &
                           sthu_tile, smcl_tile, sthf_tile,                    &
                           sthu, tsoil_tile, canht_ft, lai_ft,                 &
                           sin_theta_latitude, dzsoil,                         &
                           CPOOL_TILE, NPOOL_TILE, PPOOL_TILE, SOIL_ORDER,     &
                           NIDEP, NIFIX, PWEA, PDUST, GLAI, PHENPHASE,         &
                           WOOD_HVEST_C,WOOD_HVEST_N,WOOD_HVEST_P,             &
                           WOOD_FLUX_C,WOOD_FLUX_N,WOOD_FLUX_P,                &
                           WRESP_C,WRESP_N,WRESP_P, THINNING,                  &
                           PREV_YR_SFRAC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,  &
                           iday )

   !---------------------------------------------------------------------!
   !--- Feedback prognostic vcmax and daily LAI from casaCNP to CABLE ---!
   !---------------------------------------------------------------------!
   IF(l_vcmaxFeedbk) call casa_feedback(ktau_gl,veg,casabiome,casapool,casamet)
   IF(l_laiFeedbk) veg%vlai(:) = casamet%glai(:)

   canopy%oldcansto=canopy%cansto


   !---------------------------------------------------------------------!
   !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
   !--- req'd by Mk3L  --------------------------------------------------!
   !---------------------------------------------------------------------!
   CALL cbm( timestep, air, bgc, canopy, met, bal,                             &
             rad, rough, soil, ssnow, sum_flux, veg )

! output CO2_MMR value used in CABLE (passed from UM)
  if ( (knode_gl.eq.1) .and. (ktau_gl.eq.1) ) then
        write(6,*) 'CO2_MMR in CABLE: ',  met%ca(1)*44./28.966
  end if



   !---------------------------------------------------------------------!
   !--- pass land-surface quantities calc'd by CABLE in explicit call ---!
   !--- back to UM.                                                   ---!
   !---------------------------------------------------------------------!
   call cable_expl_unpack( FTL_TILE_CAB, FTL_CAB, FTL_TILE, FQW_TILE,          &
                           LE_TILE_CAB, LE_CAB, TSTAR_TILE, TSTAR_TILE_CAB,    &
                           TSTAR_CAB, U_S, U_S_STD_TILE, U_S_CAB, CH_CAB,      &
                           CD_CAB, CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
                           rad%transd, rough%z0m, rough%zref_tq )


   ! dump bitwise reproducible testing data
   IF( cable_user%RUN_DIAG_LEVEL == 'zero')                                    &
      call cable_diag( 1, "FLUXES", mp, kend_gl, ktau_gl, knode_gl,            &
                          "FLUXES", canopy%fe + canopy%fh )
                

   cable_runtime%um_explicit = .FALSE.


END SUBROUTINE cable_explicit_driver




!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM.                                                   ---!
!---------------------------------------------------------------------!
SUBROUTINE cable_expl_unpack( FTL_TILE_CAB, FTL_CAB, FTL_TILE, FQW_TILE,       &
                           LE_TILE_CAB, LE_CAB, TSTAR_TILE, TSTAR_TILE_CAB,    &
                           TSTAR_CAB, U_S, U_S_STD_TILE, U_S_CAB, CH_CAB,      &
                           CD_CAB, CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
                           ssnow_snowd, ssnow_cls, air_rlam, air_rho,          &
                           canopy_fe, canopy_fh, canopy_us, canopy_cdtq,       &
                           canopy_fwet, canopy_wetfac_cs, canopy_rnet,         &
                           canopy_zetar, canopy_epot, met_ua, rad_trad,        &
                           rad_transd, rough_z0m, rough_zref_tq )

   USE cable_def_types_mod, ONLY : mp, NITER 
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1
   USE cable_common_module, ONLY : cable_runtime, cable_user, &
                                   ktau_gl, knode_gl 
   IMPLICIT NONE         


   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 


   !___ UM variables to recieve unpacked CABLE vars

   !___return fluxes
   REAL, INTENT(OUT), DIMENSION(um1%land_pts) ::   &
      FTL_CAB, &
      LE_CAB
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      FTL_TILE_CAB, &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      FQW_TILE,   &  ! Surface FQW for land tiles     
      LE_TILE_CAB

   !___return temp and roughness
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      TSTAR_TILE_CAB, TSTAR_TILE,  Z0H_TILE, Z0M_TILE
   REAL, INTENT(OUT), DIMENSION(um1%land_pts) ::                  &
      TSTAR_CAB

   !___return friction velocities/drags/ etc
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      CD_TILE,    &     ! Drag coefficient
      CH_TILE,    &     ! Transfer coefficient for heat & moisture
      U_S_STD_TILE      ! Surface friction velocity
   REAL, INTENT(OUT), DIMENSION(um1%row_length,um1%rows)  :: &
      U_S               ! Surface friction velocity (m/s)
   REAL, INTENT(OUT), DIMENSION(um1%land_pts) ::                  &
      CH_CAB,  &  ! Turbulent surface exchange
      CD_CAB,  &  ! Turbulent surface exchange
      U_S_CAB     ! Surface friction velocity (m/s)

   !___return miscelaneous 
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      RADNET_TILE,   & ! Surface net radiation
      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                       ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,         & ! Total resistance factor.
                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,
                       ! 1 for snow.    
      FRACA,         & ! Fraction of surface moisture
      RECIP_L_MO_TILE,&! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
      EPOT_TILE
   
   LOGICAL,DIMENSION(um1%land_pts,um1%ntiles) :: l_tile_pts

   !___UM vars used but NOT returned 
   REAL, INTENT(IN), DIMENSION(um1%land_pts) ::   &
      FLAND(um1%land_pts)              ! IN Land fraction on land tiles.




   !___ decs of intent(in) CABLE variables to be unpacked

   ! snow depth (liquid water), factor for latent heat
   REAL, INTENT(IN), DIMENSION(mp) :: ssnow_snowd, ssnow_cls
   
   ! surface wind speed (m/s)
   REAL, INTENT(IN), DIMENSION(mp) :: met_ua 
   
   ! latent heat for water (j/kg), dry air density (kg m-3)
   REAL, INTENT(IN), DIMENSION(mp) :: air_rlam, air_rho 
   
   ! frac SW diffuse transmitted thru canopy, rad. temp. (soil and veg)
   REAL, INTENT(IN), DIMENSION(mp) :: rad_trad,rad_transd 
   
   ! total latent heat (W/m2), total sensible heat (W/m2)
   REAL, INTENT(IN), DIMENSION(mp) :: canopy_fe, canopy_fh  
   
   ! fraction of canopy wet
   REAL, INTENT(IN), DIMENSION(mp) :: canopy_fwet, canopy_wetfac_cs
   
   ! friction velocity, drag coefficient for momentum
   REAL, INTENT(IN), DIMENSION(mp) :: canopy_us, canopy_cdtq
   
   ! net rad. absorbed by surface (W/m2), total potential evaporation 
   REAL, INTENT(IN), DIMENSION(mp) :: canopy_rnet, canopy_epot        
   
   ! stability correction
   REAL, INTENT(IN), DIMENSION(mp,niter) :: canopy_zetar
   
   ! roughness length, Reference height for met forcing
   REAL, INTENT(IN), DIMENSION(mp) :: rough_z0m, rough_zref_tq 
 
   !-------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 


   
        
   !___vars in local calc. of latent heat fluxes
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
      FQW_TILE_CAB,  &
      LE_TILE

   !___vars in local calc of Surface friction velocities
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
      CD_CAB_TILE,   &  
      CH_CAB_TILE,   &  ! (bulk transfer) coeff. for momentum
      U_S_TILE
   REAL, DIMENSION(mp)  :: &
      CDCAB,CHCAB

   !___local miscelaneous
   REAL, DIMENSION(mp)  :: &
   THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
   INTEGER :: i,j,k,N,L
   REAL :: miss = 0.0
   LOGICAL, SAVE :: first_cable_call = .true.
   REAL, POINTER :: CAPP 
   
      CAPP => PHYS%CAPP
      
      !___return fluxes
      FTL_TILE_CAB = UNPACK(canopy_fh,  um1%l_tile_pts, miss)
      FTL_CAB = SUM(um1%TILE_FRAC * FTL_TILE_CAB,2)
      FQW_TILE_CAB = UNPACK(canopy_fe,  um1%l_tile_pts, miss)
      LE_TILE_CAB = UNPACK(canopy_fe,  um1%l_tile_pts, miss)
      LE_CAB = SUM(um1%TILE_FRAC * LE_TILE_CAB,2)
      fe_dlh = canopy_fe/(air_rlam*ssnow_cls)
      FTL_TILE = UNPACK(canopy_fh,  um1%l_tile_pts, miss)
      FTL_TILE = FTL_TILE / capp
      FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)

      !___return temp and roughness
      TSTAR_TILE_CAB = UNPACK(rad_trad, um1%l_tile_pts, miss)
      TSTAR_CAB = SUM(um1%TILE_FRAC * TSTAR_TILE_CAB,2)
      TSTAR_TILE = UNPACK(rad_trad,  um1%l_tile_pts, miss)
      Z0M_TILE = UNPACK(rough_z0m,  um1%l_tile_pts, miss)
      Z0H_TILE = 0.1*Z0M_TILE
      
      !___return friction velocities/drags/ etc
      U_S_TILE  =  UNPACK(canopy_us, um1%l_tile_pts, miss)
      U_S_CAB  = SUM(um1%TILE_FRAC *  U_S_TILE,2)
      CDCAB = canopy_us**2/met_ua**2   ! met%ua is always above umin = 0.1m/s
      ! for Cable CD*
      CD_CAB_TILE =  UNPACK(CDCAB,um1%l_tile_pts, miss)
      CD_CAB= SUM(um1%TILE_FRAC * CD_CAB_TILE,2)
      ! for Cable CH*
      CH_CAB_TILE =  UNPACK(canopy_cdtq,um1%l_tile_pts, miss)
      CH_CAB= SUM(um1%TILE_FRAC * CH_CAB_TILE,2)

      U_S_STD_TILE=U_S_TILE
      CD_TILE = CD_CAB_TILE
      CH_TILE = CH_CAB_TILE

      U_S = 0.
      DO N=1,um1%ntiles
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
            I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
            U_S(I,J) = U_S(I,J)+um1%TILE_FRAC(L,N)*U_S_TILE(L,N)
         ENDDO
      ENDDO




      !___return miscelaneous 
      fraca_cab = canopy_fwet * (1.-rad_transd)
      WHERE( ssnow_snowd > 1.0 ) fraca_cab = 1.0
      rfsfs_cab = MIN( 1., MAX( 0.01, canopy_wetfac_cs - fraca_cab ) /         &
                  MAX( 0.01,1. - fraca_cab ) )
      FRACA = UNPACK( fraca_cab, um1%l_tile_pts, miss )
      RESFT = UNPACK( canopy_wetfac_cs,um1%l_tile_pts, miss )
      RESFS = UNPACK( rfsfs_cab , um1%l_tile_pts, miss )

      RADNET_TILE = UNPACK( canopy_rnet , um1%l_tile_pts, miss )
      THETAST = ABS( canopy_fh ) / ( air_rho * capp*canopy_us )
      RECIPLMOTILE =  canopy_zetar(:,niter) / rough_zref_tq
      RECIP_L_MO_TILE = UNPACK( RECIPLMOTILE, um1%l_tile_pts, miss )
      EPOT_TILE = UNPACK( canopy_epot, um1%l_tile_pts, miss )
      

      IF(first_cable_call) THEN 
         l_tile_pts = um1%l_tile_pts
         first_cable_call = .FALSE.
      ENDIF

   
END SUBROUTINE cable_expl_unpack
    
!============================================================================
!============================================================================
!============================================================================

