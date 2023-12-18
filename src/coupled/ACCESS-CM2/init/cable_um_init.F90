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
! Purpose: Initialize and update CABLE variables from UM forcing, calls to 
!          memory allocation and initialization subroutines
!
! Called from: cable_explicit_driver
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: No significant change from v1.8
!
!
! ==============================================================================

MODULE cable_um_init_mod
   IMPLICIT NONE   
   PUBLIC :: interface_UM_data
   PRIVATE  

CONTAINS

SUBROUTINE interface_UM_data( row_length, rows, land_pts, ntiles,              &
                              npft, sm_levels, itimestep, latitude, longitude, &
                              land_index, tile_frac, tile_pts, tile_index,     &
                              bexp, hcon, satcon, sathh, smvcst, smvcwt,       &
                              smvccl, albsoil, slope_avg,slope_std,dz_gw,           &
                              perm_gw,drain_gw,snow_tile, snow_rho1l,    &
                              snow_age, isnow_flg3l, snow_rho3l, snow_cond,  &
                              snow_depth3l, snow_tmp3l, snow_mass3l, sw_down,  &
                              lw_down, cos_zenith_angle, surf_down_sw, ls_rain,&
                              ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,    &
                              z1_uv, rho_water, rho_ice,L_tile_pts,            &
                              visc_sublayer_depth, canopy_tile, Fland,         &
                              ! rml 2/7/13 pass 3d co2 through to cable if required
                              CO2_MMR,&
                              !r935 CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,   &
                              sthu_tile, smcl_tile, smgw_tile,sthf_tile, sthu, &
                              tsoil_tile, canht_ft, lai_ft, sin_theta_latitude,&
                              dzsoil, CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,      &
                              SOIL_ORDER, NIDEP, NIFIX, PWEA, PDUST, GLAI,     &
                              PHENPHASE, NPP_FT_ACC, RESP_W_FT_ACC )

   USE cable_um_init_subrs_mod          ! where most subrs called from here reside
  USE casa_um_inout_mod
   
   USE cable_um_tech_mod,   ONLY :                                             &
      alloc_um_interface_types,  & ! mem. allocation subr (um1, kblum%) 
      um1,                       & ! um1% type UM basics 4 convenience
      kblum_veg                    ! kblum_veg% reset UM veg vars 4 CABLE use

   USE cable_common_module, ONLY :                                             &
      cable_user,          & ! cable_user% type inherits user definition
                             ! via namelist (cable.nml) 
      l_casacnp,           & !
      knode_gl               !

   USE cable_def_types_mod, ONLY : mp ! number of points CABLE works on
  
  use cable_pft_params_mod, ONLY : cable_pft_params ! subr
  use cable_soil_params_mod, ONLY : cable_soil_params !subr
  
   !USE casa_um_inout_mod


   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 
   !___UM dimensions, array indexes, flags
   INTEGER ::                                                      &
      row_length, rows, & ! UM resolution
      land_pts,         & ! number of land_pts
      ntiles,           & ! number of tiles
      npft,             & ! number of Plant Functional Types
      sm_levels           ! number of soil layers

   INTEGER, DIMENSION(land_pts) ::                                 &
      land_index  ! index of land point 
   
   INTEGER, DIMENSION(ntiles) ::                                   &
      tile_pts    ! number of land_pts per tile type
  
   INTEGER, DIMENSION(land_pts, ntiles) ::                         &
      tile_index, &  ! index of tile 
      isnow_flg3l    ! flag for 3-layer snow 

   !___UM parameters 
   INTEGER :: itimestep
   REAL :: rho_water ,rho_ice
   REAL, DIMENSION(sm_levels) ::                                   &
      dzsoil

   !___UM soil/snow/radiation/met vars
   REAL, DIMENSION(land_pts) ::                                    &
      bexp, &     !
      hcon, &     !  
      satcon, &   !
      sathh,  &   !
      smvcst, &   !
      smvcwt, &   !
      smvccl, &   !
      albsoil,&   !
      fland       !
       
   REAL, DIMENSION(land_pts) ::                                    &
      slope_avg,&
      slope_std,&
      dz_gw,&
      sy_gw,&
      perm_gw,&
      drain_gw

   REAL, DIMENSION(row_length,rows) ::                          &
      sw_down, &        !
      cos_zenith_angle  !

   REAL, DIMENSION(row_length,rows) ::                             &
      latitude, longitude, &
      lw_down, &  !
      ls_rain, &  !
      ls_snow, &  !   
      tl_1,    &  !
      qw_1,    &  !
      vshr_land,& !
      pstar,   &  !
      z1_tq,   &  !
      z1_uv       !

   REAL, DIMENSION(land_pts, ntiles) ::                         &
      snow_tile   !
   
   REAL, DIMENSION(land_pts, ntiles) ::                            & 
      tile_frac, &   !   
      snow_rho1l,&   !
      snow_age     !

   REAL, DIMENSION(row_length, rows, 4) ::                         &
      surf_down_sw 

   REAL, DIMENSION(land_pts, npft) ::                              &
      canht_ft, & !
      lai_ft      !

   REAL,DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile !

   REAL,DIMENSION(land_pts, ntiles) ::                             &
      visc_sublayer_depth !
   REAL, DIMENSION(land_pts, ntiles,3) ::                       &
      snow_cond   !

   REAL, DIMENSION(land_pts, ntiles,3) ::                          &
      snow_rho3l, &     !
      snow_depth3l, &   ! 
      snow_mass3l,  &   ! 
      snow_tmp3l

   REAL, DIMENSION(land_pts, sm_levels) ::                         &
      sthu  !

   REAL, DIMENSION(land_pts, ntiles, sm_levels) ::                 &
      sthu_tile, &   !
      sthf_tile, &   !
      smcl_tile, &   !
      tsoil_tile     !

   REAL, DIMENSION(land_pts, ntiles) ::                 &
      smgw_tile

   REAL :: co2_mmr
! ~r935 rml 2/7/13 Extra atmospheric co2 variables
   LOGICAL  :: L_CO2_INTERACTIVE
   INTEGER ::                              &
      CO2_DIM_LEN                                     &
     ,CO2_DIM_ROW 
   REAL :: CO2_3D(1,1)  ! co2 mass mixing ratio

   LOGICAL,DIMENSION(land_pts, ntiles) ::                       &
      L_tile_pts  ! true IF vegetation (tile) fraction is greater than 0
  
   REAL, DIMENSION(row_length,rows) ::                             & 
      sin_theta_latitude

! Les 28 Spet 2012 - CASA-CNP Pools
   REAL, DIMENSION(land_pts,ntiles,10) :: &
      CPOOL_TILE, &
      NPOOL_TILE

   REAL, DIMENSION(land_pts,ntiles,12) :: &
      PPOOL_TILE

   REAL, DIMENSION(land_pts) :: &
      SOIL_ORDER, &
      NIDEP,      &
      NIFIX,      &
      PWEA,       &
      PDUST

   REAL, DIMENSION(land_pts,ntiles) :: &
      GLAI, &
   !INTEGER, DIMENSION(land_pts,ntiles) :: &
      PHENPHASE

   REAL, DIMENSION(land_pts,ntiles) :: &
      NPP_FT_ACC,   &
      RESP_W_FT_ACC

   !------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM cable_explicit_driver() -------------------------
   !------------------------------------------------------------------------- 

   !___ local decs
   !___defs 1st call to CABLE in this run. necesary in UM & coupled.
   LOGICAL, SAVE :: first_call = .TRUE. 
   INTEGER :: i,j

   INTEGER, DIMENSION(land_pts, ntiles) ::                                     &
     tile_index_mp   !tile # on each land point, ntile for packing to CABLE

   !--- logn, vegparmnew can be set thru cable.nml
   INTEGER :: logn=6       ! 6=write to std out
   LOGICAL :: vegparmnew=.true.   ! true=read std veg params false=CASA file 


      !---------------------------------------------------------------------!
      !--- code to create type um1% conaining UM basic vars describing    --! 
      !--- dimensions etc, which are used frequently throughout interfaces. !
      !--- and then assign every time explicit is called (shouldn't need to)!
      !---------------------------------------------------------------------!
      IF(first_call) THEN
         CALL alloc_um_interface_types( row_length, rows, land_pts,            &
                                        ntiles, sm_levels )
      ENDIF       

      CALL assign_um_basics_to_um1( row_length, rows, land_pts, ntiles,        &
                                    npft, sm_levels, itimestep, latitude,      &
                                    longitude, land_index, tile_frac,          &
                                    tile_pts, tile_index, l_tile_pts,          &
                                    rho_water,rho_ice  )

      !---------------------------------------------------------------------!
      !--- CABLE vars are initialized/updated from passed UM vars       ----!
      !---------------------------------------------------------------------!
      
      !---def. vector length for cable(mp) & logical l_tile_pts
      !--- IF the tile is "active"
      IF ( first_call ) THEN
      
         um1%L_TILE_PTS = .FALSE.
         mp = SUM(um1%TILE_PTS)
         
         CALL alloc_cable_types()

         DO i=1,land_pts
            DO j=1,ntiles
               
               IF( um1%TILE_FRAC(i,j) .GT. 0.0 ) THEN 
                  um1%L_TILE_PTS(i,j) = .TRUE.
                  tile_index_mp(i,j) = j 
               ENDIF
            
            ENDDO
         ENDDO
      
      ENDIF

      !jhan: turn this off until implementation finalised
      !--- initialize latitude/longitude & mapping IF required
      if ( first_call ) & 
         call initialize_maps(latitude,longitude, tile_index_mp)

      !--- read in soil (and veg) parameters 
      IF(first_call ) then
        call cable_pft_params()
        call cable_soil_params()
      endif

      !--- initialize veg   
      CALL initialize_veg( canht_ft, lai_ft, dzsoil ) 
 
      !--- initialize soil
      CALL initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt,      &
                            smvccl, albsoil, tsoil_tile, sthu, sthu_tile,   &
                            dzsoil, slope_avg, slope_std, dz_gw,            &
                            perm_gw,drain_gw ) 
      !--- initialize roughness
      CALL initialize_roughness( z1_tq, z1_uv, kblum_veg%htveg ) 
       
      !--- initialize soilsnow
      CALL initialize_soilsnow( smvcst, tsoil_tile, sthf_tile, smcl_tile,   &
                                smgw_tile, snow_tile, snow_rho1l, snow_age, &
                                isnow_flg3l, snow_rho3l, snow_cond,         &
                                snow_depth3l, snow_mass3l, snow_tmp3l,      &
                                fland, sin_theta_latitude ) 

      !--- initialize canopy   
      CALL initialize_canopy(CANOPY_TILE,visc_sublayer_depth)
 
      !--- initialize radiation & met forcing
      CALL initialize_radiation( sw_down, lw_down, cos_zenith_angle,        &
                                 surf_down_sw, sin_theta_latitude, ls_rain, &
                                 ls_snow, tl_1, qw_1, vshr_land, pstar,     &
! rml 2/7/13 pass 3d co2 through to cable if required
                   CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE )   

      IF( first_call ) THEN
         CALL init_bgc_vars() 
         CALL init_sumflux_zero() 

      !--- initialize respiration for CASA-CNP
      !--- Lestevens 23apr13
      !IF (l_casacnp) THEN ?
         CALL init_respiration(NPP_FT_ACC,RESP_W_FT_ACC)

      ! Lestevens 28 Sept 2012 - Initialize CASA-CNP here
         if (l_casacnp) then
           if (knode_gl==0) then
             print *, '  '; print *, 'CASA_log:'
             print *, '  Calling CasaCNP - Initialise '
             print *, '  l_casacnp = ',l_casacnp
             print *, 'End CASA_log:'; print *, '  '
           endif
           call init_casacnp(sin_theta_latitude,cpool_tile,npool_tile,&
                             ppool_tile,soil_order,nidep,nifix,pwea,pdust,&
                             GLAI,PHENPHASE)
         endif

         !CM2!CALL dealloc_vegin_soilin()
         first_call = .FALSE. 
      ENDIF      

 return
END SUBROUTINE interface_UM_data
                                   
!============================================================================
!============================================================================
!============================================================================

SUBROUTINE assign_um_basics_to_um1( row_length, rows, land_pts, ntiles,     &
                                    npft, sm_levels, timestep, latitude,    &
                                    longitude, land_index, tile_frac,       &
                                    tile_pts, tile_index, l_tile_pts,       &
                                    rho_water,rho_ice  )
   USE cable_um_tech_mod,   ONLY : um1
   USE cable_common_module, ONLY : cable_user

   INTEGER, INTENT(IN) :: row_length, rows, land_pts, ntiles, npft, sm_levels
   INTEGER, INTENT(IN) :: timestep 
   REAL, INTENT(IN), DIMENSION(row_length,rows) :: latitude, longitude 
   REAL,INTENT(IN):: rho_water,rho_ice
   INTEGER, INTENT(IN), DIMENSION(land_pts)  :: land_index 
   INTEGER, INTENT(IN), DIMENSION(ntiles)  :: tile_pts 
   INTEGER, INTENT(IN), DIMENSION(land_pts, ntiles)  :: tile_index
   REAL, INTENT(IN), DIMENSION(land_pts, ntiles)  :: tile_frac 
   LOGICAL, INTENT(IN), DIMENSION(land_pts,ntiles)  :: l_tile_pts 
     
      um1%row_length = row_length
      um1%rows = rows
      um1%land_pts = land_pts
      um1%ntiles = ntiles   
      um1%npft = npft   
      um1%sm_levels = sm_levels
      um1%timestep = timestep
      um1%latitude = latitude
      um1%longitude = longitude
      um1%land_index = land_index
      um1%tile_frac = tile_frac
      um1%tile_pts = tile_pts
      um1%tile_index = tile_index
      um1%rho_water= rho_water 
      um1%rho_ice  = rho_ice

END SUBROUTINE assign_um_basics_to_um1

!============================================================================
!============================================================================
!============================================================================

END MODULE cable_um_init_mod
