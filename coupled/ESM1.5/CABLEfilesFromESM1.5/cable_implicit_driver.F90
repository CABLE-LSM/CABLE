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
! Purpose: Updates CABLE variables (as altered by first pass through boundary 
!          layer and convection scheme), calls cbm, passes CABLE variables back 
!          to UM. 'Implicit' is the second call to cbm in each UM timestep.
!
! Called from: UM code sf_impl
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
! ==============================================================================


   !USE cable_data_module,   ONLY : PHYS
   !REAL, POINTER :: TFRZ
   !   TFRZ => PHYS%TFRZ
subroutine cable_implicit_driver( LS_RAIN, CON_RAIN, LS_SNOW, CONV_SNOW,       &
                                  DTL_1,DQW_1, TSOIL, TSOIL_TILE, SMCL,        &
                                  SMCL_TILE, timestep, SMVCST,STHF, STHF_TILE, &
                                  STHU,STHU_TILE, snow_tile, SNOW_RHO1L,       &
                                  ISNOW_FLG3L, SNOW_DEPTH3L, SNOW_MASS3L,      &
                                  SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,           &
                                  FTL_TILE_CAB, FTL_CAB, LE_TILE_CAB,          &
                                  LE_CAB, FTL_1, FTL_TILE, FQW_1, FQW_TILE,    &
                                  TSTAR_TILE, TSTAR_TILE_CAB, TSTAR_CAB,       &
                                  SMCL_CAB, TSOIL_CAB, SURF_HTF_CAB,           &
                                  SURF_HT_FLUX_LAND, ECAN_TILE, ESOIL_TILE,    &
                                  EI_TILE, RADNET_TILE, TOT_ALB, SNAGE_TILE,   &
                                  CANOPY_TILE, GS, GS_TILE,T1P5M_TILE, Q1P5M_TILE,     &
                                  CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,        &
                                  DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,   &
                                  RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT,  &
                                  G_LEAF, TRANSP_TILE, CPOOL_TILE, NPOOL_TILE, &
                                  PPOOL_TILE, GLAI, PHENPHASE,  &
                                  !WOOD_HVEST_C,WOOD_HVEST_N,WOOD_HVEST_P,      &
                                  !WOOD_FLUX_C,WOOD_FLUX_N,WOOD_FLUX_P,         &
                                  !WRESP_C,WRESP_N,WRESP_P,THINNING,            &
                                  NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,          &
                                  FNSNET,FNLEACH,FNUP,FNLOSS,FNDEP,FNFIX,idoy )

   USE cable_def_types_mod, ONLY : mp
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1, conv_rain_prevstep, conv_snow_prevstep,&
                                  air, bgc, canopy, met, bal, rad, rough,      &
                                  ssnow, sum_flux
    USE cable_params_mod, ONLY : veg => veg_cbl 
    USE cable_params_mod, ONLY : soil => soil_cbl 
   USE cable_common_module, ONLY : cable_runtime, cable_user, l_casacnp,       &
                                   l_vcmaxFeedbk, knode_gl, ktau_gl, kend_gl
   USE cable_um_init_subrs_mod, ONLY : um2cable_rr
   USE cable_cbm_module,    ONLY : cbm

   USE casavariable
   USE phenvariable
   USE casa_types_mod
   USE bgcdriver_mod, ONLY : bgcdriver
   USE sumcflux_mod, ONLY : sumcflux
   USE casa_um_inout_mod
USE cbl_rhoch_ESM1pt5_module, ONLY : rhoch =>rhoch_gl,  &
                                     c1 => c1_gl, &
                                     xk => xk_gl

   IMPLICIT NONE
        
   REAL, DIMENSION(um1%ROW_LENGTH,um1%ROWS) ::                                 &
      LS_RAIN,  & ! IN Large scale rain
      LS_SNOW,  & ! IN Large scale snow
      CON_RAIN, & ! IN Convective rain
      CONV_SNOW,& ! IN Convective snow
      DTL_1,    & ! IN Level 1 increment to T field 
      DQW_1       ! IN Level 1 increment to q field 

   REAL :: timestep

   INTEGER ::                                                                  &
      DIM_CS1, DIM_CS2 

   REAL, DIMENSION(um1%land_pts) ::                                            &
      GS,      &  ! OUT "Stomatal" conductance to
      SMVCST,  &  ! IN Volumetric saturation point
      FLAND       ! IN Land fraction on land tiles
   
   REAL, DIMENSION(um1%ROW_LENGTH,um1%ROWS) ::                                 &
      !--- Net downward heat flux at surface over land.
      !--- fraction of gridbox (W/m2).
      SURF_HT_FLUX_LAND,                                                       &   
      !--- Moisture flux between layers. (kg/m^2/sec).
      !--- FQW(,1) is total water flux from surface, 'E'.
      FQW_1,                                                                   &  
      !--- FTL(,K) =net turbulent sensible heat flux into layer K
      !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
      FTL_1         

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                              &
      !___Surface FTL, FQL for land tiles
      FTL_TILE, FQW_TILE, FQW_TILE_CAB,                                     &
      
      !___(tiled) latent heat flux, melting, stomatatal conductance
     LE_TILE, MELT_TILE, GS_TILE,                                           &
     
     !___ INOUT Surface net radiation on tiles (W/m2)
     RADNET_TILE, &
     TOT_ALB,     & ! total albedo
     EI_TILE,     & ! OUT EI for land tiles.
     ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
     ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
   !___ runoff ??
   REAL, dimension(um1%land_pts,um1%sm_levels) ::                           &
      SMCL,       & ! 
      STHF,       & !
      STHU,       & !
      SMCL_CAB,   & !
      TSOIL_CAB,  & !
      TSOIL,      & !
      SURF_CAB_ROFF !       

   !___(tiled) soil prognostics: as above 
   REAL, dimension(um1%land_pts,um1%ntiles,um1%sm_levels) ::                &
      SMCL_TILE, & !
      STHU_TILE, & !
      TSOIL_TILE,& !
      STHF_TILE    !

   !___flag for 3 layer snow pack
   INTEGER :: ISNOW_FLG3L(um1%LAND_PTS,um1%NTILES)
   
   !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
   REAL, dimension(um1%land_pts,um1%ntiles,3) :: &
      SNOW_DEPTH3L,  & ! 
      SNOW_MASS3L,   & !
      SNOW_RHO3L,    & !
      SNOW_TMP3L,    & !
      SNOW_COND        !

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                              &
      FRS_TILE,   & ! Local
      NEE_TILE,   & ! Local
      NPP_TILE,   & ! Local
      GPP_TILE,   & ! Local
      GLEAF_TILE, & ! Local, kdcorbin, 10/10
      FRP_TILE,   &
      NPP_FT,     &
      NPP_FT_old, &
      GPP_FT,     &
      GPP_FT_old       

   REAL, DIMENSION(um1%land_pts) ::                                         &
      SNOW_GRD,    & !
      CANOPY_GB,   & !
      FTL_CAB,     & !   
      LE_CAB,      & !
      TSTAR_CAB,   & !  
      SURF_HTF_CAB,& !
      RESP_P,      & !
      NPP,         & !
      GPP            !
      
   REAL, DIMENSION( um1%land_pts,um1%ntiles ) ::                               &
      SNOW_TILE,     &
      SNOW_RHO1L,    &  ! Mean snow density
      SNAGE_TILE,    &
      CANOPY_TILE,   &
      FTL_TILE_CAB,  & 
      LE_TILE_CAB,   &
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      TSTAR_TILE_CAB,&
      TSTAR_TILE,    &
      SURF_HTF_T_CAB,& 
      RESP_S_TILE,   & 
      RESP_P_FT,     &
      RESP_P_FT_old, &
      G_LEAF,        &
      TRANSP_TILE

   REAL ::                                                                     &
      RESP_S(um1%LAND_PTS,DIM_CS1),     &
      RESP_S_old(um1%LAND_PTS,DIM_CS1), &
      RESP_S_TOT(DIM_CS2)    

   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES,10) ::                              &
      CPOOL_TILE, &
      NPOOL_TILE     
   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES,12) ::                              &
      PPOOL_TILE
   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
      GLAI,       &
   !INTEGER, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                              &
      PHENPHASE!,  &
      !WOOD_FLUX_C,&
      !WOOD_FLUX_N,&
      !WOOD_FLUX_P,&
      !THINNING

   !REAL, DIMENSION(um1%LAND_PTS,um1%NTILES,3) ::                                 &
   !   WOOD_HVEST_C,&
   !   WOOD_HVEST_N,&
   !   WOOD_HVEST_P!,&
      !WRESP_C,&
      !WRESP_N,&
      !WRESP_P

   ! Lestevens 23apr13
   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
      NPP_FT_ACC, &
      RESP_W_FT_ACC, &
      RESP_S_ACC

   ! TZ 010219 nitrogen fluxes 
   REAL, DIMENSION( um1%land_pts,um1%ntiles ) ::                               &
      FNSNET,&
      FNLEACH,&
      FNUP,&
      FNLOSS,&
      FNDEP,& 
      FNFIX 

   INTEGER ::     &
      ktauday,    &  ! day counter for CASA-CNP
      idoy           ! day of year (1:365) counter for CASA-CNP
   INTEGER, SAVE :: &
      kstart = 1

   REAL, DIMENSION(mp) ::                                                      & 
      dtlc, & 
      dqwc

   REAL, POINTER :: TFRZ
   
      TFRZ => PHYS%TFRZ
   
      ! FLAGS def. specific call to CABLE from UM
      cable_runtime%um_explicit = .FALSE.
      cable_runtime%um_implicit = .TRUE.
   
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
      CALL um2cable_rr( (LS_RAIN+CON_RAIN)*um1%TIMESTEP, met%precip)
      CALL um2cable_rr( (LS_SNOW+CONV_SNOW)*um1%TIMESTEP, met%precip_sn)
      CALL um2cable_rr( dtl_1, dtlc)
      CALL um2cable_rr( dqw_1, dqwc)
      
      !--- conv_rain(snow)_prevstep are added to precip. in explicit call
      CALL um2cable_rr( (CON_RAIN)*um1%TIMESTEP, conv_rain_prevstep)
      CALL um2cable_rr( (CONV_snow)*um1%TIMESTEP, conv_snow_prevstep)
      
      met%precip   =  met%precip + met%precip_sn
      met%tk = met%tk + dtlc
      met%qv = met%qv + dqwc
      met%tvair = met%tk
      met%tvrad = met%tk
      met%doy = idoy + 1
 
      canopy%cansto = canopy%oldcansto

   CALL cbm( timestep, air, bgc, canopy, met, bal,                             &
             rad, rough, soil, ssnow, sum_flux, veg, xk, c1, rhoch )

      ! Lestevens - temporary ?
      ktauday = int(24.0*3600.0/TIMESTEP)
      IF(idoy==0) idoy =365
  
! Lestevens Sept2012 - Call CASA-CNP
      if (l_casacnp) then
      CALL bgcdriver(ktau_gl,kstart,kend_gl,TIMESTEP,met,ssnow,canopy,veg,soil, &
                     casabiome,casapool,casaflux,casamet,casabal,phen,          &
                     .FALSE., .FALSE., ktauday, idoy, .FALSE., .FALSE. )
!                     spinConv, spinup, ktauday, idoy, cable_user%casa_dump_read,&
!                     cable_user%casa_dump_write )
      endif

      CALL sumcflux(ktau_gl,kstart,kend_gl,TIMESTEP,bgc,canopy,soil,ssnow,      &
                    sum_flux,veg,met,casaflux,l_vcmaxFeedbk)

      CALL implicit_unpack( TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                            FTL_TILE_CAB, FTL_CAB, LE_TILE_CAB, LE_CAB,        &
                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            TSTAR_TILE_CAB, TSTAR_CAB, SMCL_CAB, TSOIL_CAB,    &
                            SURF_HTF_CAB, SURF_HT_FLUX_LAND, ECAN_TILE,        &
                            ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
                            SNAGE_TILE, CANOPY_TILE, GS, GS_TILE, T1P5M_TILE,           &
                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,  &
                            DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF,& 
                            TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC, RESP_S_ACC,&
                            FNSNET,FNLEACH,FNUP,FNLOSS,FNDEP,FNFIX)

! Lestevens Sept2012 - Call CASA-CNP
      if (l_casacnp) then
      !if (l_casacnp .and. ktau_gl==kend_gl ) then
        if (knode_gl==0 .and. ktau_gl==kend_gl) then
        !if (knode_gl==0) then
         print *, '  '; print *, 'CASA_log:'
         print *, '  Calling CasaCNP - Poolout '
         print *, '  l_casacnp = ',l_casacnp
         print *, '  ktau_gl, kend_gl = ',ktau_gl,kend_gl
         print *, 'End CASA_log:'; print *, '  '
        endif
       CALL casa_poolout_unpk(casapool,casaflux,casamet,casabal,phen,  &
                              CPOOL_TILE,NPOOL_TILE,PPOOL_TILE, &
                              GLAI,PHENPHASE)
      endif
       
      cable_runtime%um_implicit = .FALSE.
  
END SUBROUTINE cable_implicit_driver


!========================================================================= 
!========================================================================= 
!========================================================================= 
        
SUBROUTINE implicit_unpack( TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                            FTL_TILE_CAB, FTL_CAB, LE_TILE_CAB, LE_CAB,        &
                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            TSTAR_TILE_CAB, TSTAR_CAB, SMCL_CAB, TSOIL_CAB,    &
                            SURF_HTF_CAB, SURF_HT_FLUX_LAND, ECAN_TILE,        &
                            ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
                            SNAGE_TILE, CANOPY_TILE, GS, GS_TILE, T1P5M_TILE,           &
                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,  &
                            DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF,&
                            TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC, RESP_S_ACC,&
                            FNSNET,FNLEACH,FNUP,FNLOSS,FNDEP,FNFIX)
 
   USE cable_def_types_mod, ONLY : mp
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1 ,canopy, rad, ssnow, air
    USE cable_params_mod, ONLY : soil => soil_cbl 
   USE cable_common_module, ONLY : cable_runtime, cable_user
   USE casa_types_mod
   IMPLICIT NONE
 
   !jhan:these need to be cleaned out to what is actualllly passed
   INTEGER :: DIM_CS1 ,DIM_CS2 

   REAL, DIMENSION(um1%land_pts) ::                                            &
      GS,         &  ! OUT "Stomatal" conductance to
      SMVCST,     &  ! IN Volumetric saturation point
      FLAND          ! IN Land fraction on land tiles
   
   real, dimension(um1%ROW_LENGTH,um1%ROWS) ::                                 &
      !--- Net downward heat flux at surface over land.
      !--- fraction of gridbox (W/m2).
      SURF_HT_FLUX_LAND,           &
      !--- Moisture flux between layers. (kg/m^2/sec).
      !--- FQW(,1) is total water flux from surface, 'E'.
      FQW_1,       &  
      !--- FTL(,K) =net turbulent sensible heat flux into layer K
      !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
      FTL_1         

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      !___Surface FTL, FQL for land tiles
      FTL_TILE, FQW_TILE, FQW_TILE_CAB,   &  
      !___(tiled) latent heat flux, melting, stomatatal conductance
      LE_TILE, MELT_TILE, GS_TILE,     &  
      RADNET_TILE, & ! INOUT Surface net radiation on tiles (W/m2)
      TOT_ALB,     & ! total albedo
      EI_TILE,     & ! OUT EI for land tiles.
      ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
      ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
   !___ runoff ??
   REAL, DIMENSION(um1%land_pts,um1%sm_levels) ::                              &
      SMCL,       & !
      STHF,       &
      STHU,       &
      SMCL_CAB,   &
      TSOIL_CAB,  &
      TSOIL,      &
      SURF_CAB_ROFF !      

   !___(tiled) soil prognostics: as above 
   REAL, DIMENSION(um1%land_pts,um1%ntiles,um1%sm_levels) ::                   &
      SMCL_TILE,  & 
      STHU_TILE,  &
      TSOIL_TILE, &
      STHF_TILE  

   !___flag for 3 layer snow pack
   INTEGER :: ISNOW_FLG3L(um1%LAND_PTS,um1%NTILES)
   
   !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
   REAL, DIMENSION(um1%land_pts,um1%ntiles,3) ::                               &
      SNOW_DEPTH3L,  &
      SNOW_MASS3L,   &
      SNOW_RHO3L,    &
      SNOW_TMP3L,    &
      SNOW_COND 

   REAL, dimension(um1%land_pts,um1%ntiles) ::                                 &
      FRS_TILE,   & ! Local
      NEE_TILE,   & ! Local
      NPP_TILE,   & ! Local
      GPP_TILE,   & ! Local
      GLEAF_TILE, & ! Local, kdcorbin, 10/10
      FRP_TILE,   &
      NPP_FT,     &
      NPP_FT_old, &
      GPP_FT,     &
      GPP_FT_old       

   REAL, dimension(um1%land_pts) ::                                            &
      SNOW_GRD,   &  
      CANOPY_GB,  &
      FTL_CAB,    &
      LE_CAB,     &
      TSTAR_CAB,  &
      SURF_HTF_CAB,&
      RESP_P,     & 
      NPP,        & 
      GPP

   ! Lestevens 23apr13
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      NPP_FT_ACC,    & ! sresp for CASA-CNP
      RESP_W_FT_ACC, & ! presp for CASA-CNP
      RESP_S_ACC       ! NPP for CASA-CNP   

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      SNOW_TILE,     & !
      SNOW_RHO1L,    & ! Mean snow density
      SNAGE_TILE,    & !
      CANOPY_TILE,   & !
      FTL_TILE_CAB,  & !
      LE_TILE_CAB,   & !
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      TSTAR_TILE_CAB,&
      TSTAR_TILE,    &
      SURF_HTF_T_CAB,& 
      RESP_S_TILE,   & 
      RESP_P_FT,     &
      RESP_P_FT_old, &
      G_LEAF,        &
      TRANSP_TILE

   ! TZ 010219 nitrogen fluxes
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      FNSNET,&
      FNLEACH,&
      FNUP,&
      FNLOSS,&
      FNDEP,&
      FNFIX

   REAL ::                                                                     &
      RESP_S(um1%LAND_PTS,DIM_CS1),    & !
      RESP_S_old(um1%LAND_PTS,DIM_CS1),& !
      RESP_S_TOT(DIM_CS2)                !
  
   REAL, DIMENSION(mp) ::                                                                     &
      fe_dlh,    & !
      fes_dlh,   & !
      fev_dlh      !

   !--- Local vars
   INTEGER :: i,j,l,k,n

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
         !--- Local buffer surface FTL, FQL @ prev dt
         FTL_TILE_old, FQW_TILE_old

   INTEGER:: i_miss = 0
   REAL :: miss = 0.0
   
   REAL, POINTER :: TFRZ
   
   
      TFRZ => PHYS%TFRZ
  
      !--- set UM vars to zero
      TSOIL_CAB = 0.; SMCL_CAB = 0.; TSOIL_TILE = 0.; 
      SMCL_TILE = 0.; STHF_TILE = 0.; STHU_TILE = 0.

      DO j = 1,um1%SM_LEVELS
         TSOIL_TILE(:,:,j)= UNPACK(ssnow%tgg(:,j), um1%L_TILE_PTS, miss)
         TSOIL_CAB(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)
         SMCL_TILE(:,:,j)= UNPACK(REAL(ssnow%wb(:,j)), um1%L_TILE_PTS, miss)
         SMCL_TILE(:,:,j)=SMCL_TILE(:,:,j)*soil%zse(j)*um1%RHO_WATER
         SMCL_CAB(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
         STHF_TILE(:,:,j)= UNPACK(REAL(ssnow%wbice(:,j)), um1%L_TILE_PTS, miss)
         SMCL(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
         TSOIL(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)
         
         DO N=1,um1%NTILES
            DO K=1,um1%TILE_PTS(N)
               I = um1%TILE_INDEX(K,N)
               IF ( SMVCST(I) > 0. ) THEN ! Exclude permanent ice - mrd
                  STHF_TILE(I,N,J)= STHF_TILE(I,N,J)/SMVCST(I)
                  STHU_TILE(I,N,J)= MAX( 0., SMCL_TILE(I,N,J) -                &
                                    STHF_TILE(I,N,J) * SMVCST(I) * soil%zse(J) &
                                    * um1%RHO_WATER ) / ( soil%zse(J) *        &
                                    um1%RHO_WATER * SMVCST(I) )
               ENDIF
            ENDDO
         ENDDO

         STHF(:,J) = SUM(um1%TILE_FRAC * STHF_TILE(:,:,J),2)
         STHU(:,J) = SUM(um1%TILE_FRAC * STHU_TILE(:,:,J),2)
      ENDDO

      !--- unpack snow vars 
      SNOW_RHO1L  = UNPACK(ssnow%ssdnn, um1%L_TILE_PTS, miss)
      ISNOW_FLG3L = UNPACK(ssnow%isflag, um1%L_TILE_PTS, i_miss)
      MELT_TILE   = UNPACK(ssnow%smelt, um1%L_TILE_PTS, miss)
      SNOW_TILE= UNPACK(ssnow%snowd, um1%L_TILE_PTS, miss)
      SNOW_GRD=  SUM(um1%TILE_FRAC * SNOW_TILE,2)  ! gridbox snow mass & snow below canopy 

      !--- unpack layered snow vars 
      do k = 1,3
        SNOW_TMP3L(:,:,k) = UNPACK(ssnow%tggsn(:,k), um1%L_TILE_PTS, miss)
        SNOW_MASS3L(:,:,k)= UNPACK(ssnow%smass(:,k), um1%L_TILE_PTS, miss)
        SNOW_RHO3L(:,:,k) = UNPACK(ssnow%ssdn(:,k), um1%L_TILE_PTS, miss)
        SNOW_COND(:,:,k)  = UNPACK(ssnow%sconds(:,k),um1%L_TILE_PTS,miss)
        SNOW_DEPTH3L(:,:,k)  = UNPACK(ssnow%sdepth(:,k),um1%L_TILE_PTS,miss)
      enddo

      
      canopy%gswx_T = canopy%gswx_T/air%cmolar
      GS_TILE = UNPACK(canopy%gswx_T,um1%L_TILE_PTS,miss)
      GS =  SUM(um1%TILE_FRAC * GS_TILE,2)

      !___return fluxes
      FTL_TILE_CAB = UNPACK(canopy%fh,  um1%l_tile_pts, miss)
      FTL_CAB = SUM(um1%TILE_FRAC * FTL_TILE_CAB,2)
      FQW_TILE_CAB = UNPACK(canopy%fe,  um1%l_tile_pts, miss)
      LE_TILE_CAB = FQW_TILE_CAB 
      LE_CAB = SUM(um1%TILE_FRAC * LE_TILE_CAB,2)
!      fe_dlh = canopy%fe/(air%rlam*ssnow%cls)
      fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
      fev_dlh = canopy%fev/air%rlam
      fe_dlh =  fev_dlh + fes_dlh

      !---preserve fluxes from the previous time step for the coastal grids
      FTL_TILE_old = FTL_TILE
      FQW_TILE_old = FQW_TILE
      
      !---update fluxes 
      FTL_TILE = FTL_TILE_CAB 
      FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)


      !___return temp and roughness
      TSTAR_TILE_CAB = UNPACK(rad%trad, um1%l_tile_pts, miss)
      TSTAR_CAB = SUM(um1%TILE_FRAC * TSTAR_TILE_CAB,2)
      TSTAR_TILE = TSTAR_TILE_CAB 

      !___return miscelaneous 
      RADNET_TILE = unpack( canopy%rnet , um1%l_tile_pts, miss)

      SURF_HTF_T_CAB = UNPACK(canopy%ga,um1%L_TILE_PTS,miss)
      SURF_HTF_CAB = SUM(um1%TILE_FRAC * SURF_HTF_T_CAB,2)

     TOT_ALB=UNPACK(rad%albedo_T,um1%L_TILE_PTS, miss) 
     ESOIL_TILE = UNPACK(fes_dlh, um1%L_TILE_PTS, miss)
     ECAN_TILE = UNPACK(fev_dlh,  um1%L_TILE_PTS, miss)
     EI_TILE = 0.
     SNAGE_TILE = UNPACK(ssnow%snage, um1%L_TILE_PTS, miss) 
     TRANSP_TILE = UNPACK(canopy%fevc, um1%L_TILE_PTS, miss) 

     !unpack screen level (1.5m) variables
     !Convert back to K 
     t1p5m_tile     = UNPACK(canopy%tscrn+tfrz, um1%L_TILE_PTS, miss)
     q1p5m_tile     = UNPACK(canopy%qscrn, um1%L_TILE_PTS, miss)
     CANOPY_TILE    = UNPACK(canopy%cansto, um1%L_TILE_PTS, miss)
     CANOPY_GB      = SUM(um1%TILE_FRAC * CANOPY_TILE,2)

     ! TZ 010219 nitrogen fluxes
    FNSNET        = UNPACK(casaflux%Nsnet, um1%L_TILE_PTS, miss)
    FNLEACH       = UNPACK(casaflux%Nminleach, um1%L_TILE_PTS, miss)
    FNUP          = UNPACK(casaflux%Nminuptake, um1%L_TILE_PTS, miss)
    FNLOSS        = UNPACK(casaflux%Nminloss, um1%L_TILE_PTS, miss)
    FNDEP         = UNPACK(casaflux%Nmindep, um1%L_TILE_PTS, miss)
    FNFIX         = UNPACK(casaflux%Nminfix, um1%L_TILE_PTS, miss)

     ! Lestevens - Passing CO2 from CABLE to bl_trmix_dd.F90
!     FRS_TILE       = UNPACK(canopy%frs, um1%L_TILE_PTS, miss)
!     NEE_TILE       = UNPACK(canopy%fnee, um1%L_TILE_PTS, miss)
!     NPP_TILE       = UNPACK(canopy%fnpp, um1%L_TILE_PTS, miss)
!     GLEAF_TILE     = UNPACK(canopy%frday,um1%L_TILE_PTS, miss)

! TZ: output casa fluxes instead of canopy fluxes
     FRS_TILE    = UNPACK((casaflux%crsoil)/86400.0, um1%L_TILE_PTS, miss)
     NPP_TILE    = UNPACK((casaflux%cnpp)/86400.0, um1%L_TILE_PTS, miss)
     NEE_TILE    = UNPACK((casaflux%crsoil-casaflux%cnpp)/86400.0, um1%L_TILE_PTS, miss)
     GLEAF_TILE  = UNPACK((casaflux%crmplant(:,1))/86400.0, um1%L_TILE_PTS, miss)
     FRP_TILE    = UNPACK((casaflux%crmplant(:,2)+casaflux%crmplant(:,3)+casaflux%crgplant)/86400.0, um1%L_TILE_PTS, miss)
     IF( cable_user%leaf_respiration == 'on' .OR. cable_user%leaf_respiration == 'ON') THEN
        GPP_TILE = UNPACK((casaflux%cnpp+casaflux%crmplant(:,2)+casaflux%crmplant(:,3)+casaflux%crgplant)/86400.0, um1%L_TILE_PTS, miss)
     ELSE  
        GPP_TILE = UNPACK((casaflux%cnpp+casaflux%crmplant(:,1)+casaflux%crmplant(:,2)+casaflux%crmplant(:,3)+casaflux%crgplant)/86400.0, um1%L_TILE_PTS, miss)
     ENDIF

!      IF( cable_user%leaf_respiration == 'on' .OR.                             &
!           cable_user%leaf_respiration == 'ON') THEN
!         GPP_TILE = UNPACK(canopy%fnpp+canopy%frp, um1%L_TILE_PTS, miss)
!      ELSE 
!         GPP_TILE = UNPACK(canopy%fnpp+canopy%frp+canopy%frday,  &
!                            um1%L_TILE_PTS, miss)
!      ENDIF

!     FRP_TILE       = UNPACK(canopy%frp, um1%L_TILE_PTS, miss)
     NPP_FT_old     = NPP_FT
     GPP_FT_old     = GPP_FT
     RESP_P_FT_old  = RESP_P_FT
     RESP_S_old     = RESP_S

     !initialse full land grids and retain coastal grid fluxes
      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           J=(um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
           I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
           IF( FLAND(L) == 1.0) THEN 
             FTL_1(I,J) =  0.0
             FQW_1(I,J) =  0.0
           ELSE
             !retain sea/ice contribution and remove land contribution
             FTL_1(I,J) = FTL_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *         &
                          FTL_TILE_old(L,N)
             FQW_1(I,J) = FQW_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *         &
                          FQW_TILE_old(L,N)
           ENDIF
           SURF_HT_FLUX_LAND(I,J) = 0.
         ENDDO
     ENDDO

      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           J=(um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
           I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
           FTL_1(I,J) = FTL_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FTL_TILE(L,N)
           FQW_1(I,J) = FQW_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FQW_TILE(L,N)
           SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                   &
                                    FLAND(L)*um1%TILE_FRAC(L,N) *              &
                                    SURF_HTF_T_CAB(L,N)
         ENDDO
      ENDDO

      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            IF( FLAND(L) == 1.0) THEN
               NPP(L)=0.; NPP_FT(L,N)=0.; GPP(L)=0.; GPP_FT(L,N)=0.
               RESP_P(L)=0.; RESP_P_FT(L,N)=0.; RESP_S(L,:)=0.; G_LEAF(L,N)=0.   
            ELSE
               ! For coastal points: currently no contribution
               NPP(L)=NPP(L)-FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT_old(L,N)
               GPP(L)=GPP(L)-FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT_old(L,N)
               RESP_P(L)=RESP_P(L)-FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT_old(L,N)
               !--- loop for soil respiration
               DO I=1,DIM_CS1
                  RESP_S(L,I)=RESP_S(L,I)-FLAND(L)*RESP_S_old(L,I)
               ENDDO
               RESP_S_TOT(L)=sum(RESP_S(L,:))
            ENDIF
         ENDDO
      ENDDO

     RESP_S_TILE=FRS_TILE*1.e-3
     ! Lestevens 23apr13 - possible miss match ntiles<-->npft
     DO N=1,um1%NTILES
        DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           !---convert units to kg C m-2 s-1
           NPP_FT_ACC(L,N)    = FRS_TILE(L,N)*1.e-3
           RESP_W_FT_ACC(L,N) = FRP_TILE(L,N)*1.e-3
           RESP_S_ACC(L,N) = NPP_TILE(L,N)*1.e-3
        ENDDO
     ENDDO

      DO N=1,um1%NTILES 
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            !add leaf respiration to output
            G_LEAF(L,N)=GLEAF_TILE(L,N)*1.e-3
            NPP_FT(L,N)=NPP_TILE(L,N)*1.e-3
            NPP(L)=NPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT(L,N)
            GPP_FT(L,N)=GPP_TILE(L,N)*1.e-3
            GPP(L)=GPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT(L,N)

            !loop for soil resp. - all UM levels = single CABLE output 
            DO I=1,DIM_CS1
               RESP_S(L,I) = RESP_S(L,I) + &
                             FLAND(L)*um1%TILE_FRAC(L,N)*FRS_TILE(L,N)*1.e-3
            ENDDO

            RESP_S_TOT(L)=sum(RESP_S(L,:))
            RESP_P_FT(L,N)=FRP_TILE(L,N)*1.e-3
            RESP_P(L)=RESP_P(L)+FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT(L,N)
         ENDDO
      ENDDO

END SUBROUTINE Implicit_unpack


