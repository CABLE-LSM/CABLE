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
                              smvccl, albsoil, snow_tile, snow_rho1l,          &
                              snage_tile, isnow_flg3l, snow_rho3l, snow_cond,  &
                              snow_depth3l, snow_tmp3l, snow_mass3l, sw_down,  &
                              lw_down, cos_zenith_angle, surf_down_sw, ls_rain,&
                              ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,    &
                              z1_uv, rho_water, L_tile_pts, canopy_tile, Fland,&
! rml 2/7/13 pass 3d co2 through to cable if required
                   CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,   &
                              sthu_tile, smcl_tile, sthf_tile, sthu,           &
                              tsoil_tile, canht_ft, lai_ft, sin_theta_latitude,&
                              dzsoil, CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,      &
                              SOIL_ORDER, NIDEP, NIFIX, PWEA, PDUST, GLAI,     &
                              PHENPHASE,WOOD_HVEST_C,WOOD_HVEST_N,WOOD_HVEST_P,& 
                              WOOD_FLUX_C,WOOD_FLUX_N,WOOD_FLUX_P,             & 
                              WRESP_C,WRESP_N,WRESP_P,THINNING,                & 
                              PREV_YR_SFRAC, NPP_FT_ACC, RESP_W_FT_ACC,        &
                              RESP_S_ACC, iday )

   USE cable_um_init_subrs_mod          ! where most subrs called from here reside
   
   USE cable_um_tech_mod,   ONLY :                                             &
      alloc_um_interface_types,  & ! mem. allocation subr (um1, kblum%) 
      dealloc_vegin_soilin,      & ! mem. allocation subr (vegin%,soilin%)
      um1,soil,                  & ! um1% type UM basics 4 convenience
      kblum_veg                    ! kblum_veg% reset UM veg vars 4 CABLE use

   USE cable_common_module, ONLY :                                             &
      cable_user,          & ! cable_user% type inherits user definition
                             ! via namelist (cable.nml) 
      get_type_parameters, & ! veg and soil parameters READ subroutine  
                             !
      l_casacnp,           & !
      knode_gl               !

   USE cable_def_types_mod, ONLY : mp, mland ! number of points CABLE works on

   USE casa_um_inout_mod


   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 
   !___UM dimensions, array indexes, flags
   INTEGER, INTENT(IN) ::                                                      &
      row_length, rows, & ! UM resolution
      land_pts,         & ! number of land_pts
      ntiles,           & ! number of tiles
      npft,             & ! number of Plant Functional Types
      sm_levels           ! number of soil layers

   INTEGER, INTENT(IN), DIMENSION(land_pts) ::                                 &
      land_index  ! index of land point 
   
   INTEGER, INTENT(IN), DIMENSION(ntiles) ::                                   &
      tile_pts    ! number of land_pts per tile type
  
   INTEGER, INTENT(IN), DIMENSION(land_pts, ntiles) ::                         &
      tile_index, &  ! index of tile 
      isnow_flg3l    ! flag for 3-layer snow 

   !___UM parameters 
   INTEGER, INTENT(IN) :: itimestep
   REAL, INTENT(IN) :: rho_water 
   REAL, INTENT(IN), DIMENSION(sm_levels) ::                                   &
      dzsoil

   !___UM soil/snow/radiation/met vars
   REAL, INTENT(IN), DIMENSION(land_pts) ::                                    &
      bexp, &     !
      hcon, &     !  
      satcon, &   !
      sathh,  &   !
      smvcst, &   !
      smvcwt, &   !
      smvccl, &   !
      albsoil,&   !
      fland       !
       
   REAL, INTENT(INOUT), DIMENSION(row_length,rows) ::                          &
      sw_down, &        !
      cos_zenith_angle  !

   REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             &
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

   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles) ::                         &
      snow_tile   !
   
   REAL, INTENT(IN), DIMENSION(land_pts, ntiles) ::                            & 
      tile_frac, &   !   
      snow_rho1l,&   !
      snage_tile     !

   REAL, INTENT(IN), DIMENSION(row_length, rows, 4) ::                         &
      surf_down_sw 

   REAL, INTENT(IN), DIMENSION(land_pts, npft) ::                              &
      canht_ft, & !
      lai_ft      !

   REAL, INTENT(IN),DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile !

   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles,3) ::                       &
      snow_cond   !

   REAL, INTENT(IN), DIMENSION(land_pts, ntiles,3) ::                          &
      snow_rho3l, &     !
      snow_depth3l, &   ! 
      snow_mass3l,  &   ! 
      snow_tmp3l

   REAL, INTENT(IN), DIMENSION(land_pts, sm_levels) ::                         &
      sthu  !

   REAL, INTENT(IN), DIMENSION(land_pts, ntiles, sm_levels) ::                 &
      sthu_tile, &   !
      sthf_tile, &   !
      smcl_tile, &   !
      tsoil_tile     !

   REAL, INTENT(IN) :: co2_mmr
! rml 2/7/13 Extra atmospheric co2 variables
   LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
   INTEGER, INTENT(IN) ::                              &
      CO2_DIM_LEN                                      &
     ,CO2_DIM_ROW &
     , iday
   REAL, INTENT(IN) :: CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)  ! co2 mass mixing ratio

   LOGICAL, INTENT(INOUT),DIMENSION(land_pts, ntiles) ::                       &
      L_tile_pts  ! true IF vegetation (tile) fraction is greater than 0
  
   REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             & 
      sin_theta_latitude

! Les 28 Spet 2012 - CASA-CNP Pools
   REAL, INTENT(INOUT), DIMENSION(land_pts,ntiles,10) :: &
      CPOOL_TILE, &
      NPOOL_TILE

   REAL, INTENT(INOUT), DIMENSION(land_pts,ntiles,12) :: &
      PPOOL_TILE

   REAL, INTENT(INOUT), DIMENSION(land_pts) :: &
      SOIL_ORDER, &
      NIDEP,      &
      NIFIX,      &
      PWEA,       &
      PDUST

   REAL, INTENT(INOUT), DIMENSION(land_pts,ntiles) :: &
      GLAI,         &
   !INTEGER, INTENT(INOUT), DIMENSION(land_pts,ntiles) :: &
      PHENPHASE,    &
      PREV_YR_SFRAC,&
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
      NPP_FT_ACC,   &
      RESP_W_FT_ACC, &
      RESP_S_ACC

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
                                    rho_water  )


      !print *, 'DEBUG tile_frac', tile_frac(6805,:), um1%tile_frac(6805,:)



      !---------------------------------------------------------------------!
      !--- CABLE vars are initialized/updated from passed UM vars       ----!
      !---------------------------------------------------------------------!
      
      !---def. vector length for cable(mp) & logical l_tile_pts
      !--- IF the tile is "active"
      IF ( first_call ) THEN
      
         um1%L_TILE_PTS = .FALSE.
         mp = SUM(um1%TILE_PTS)
         mland = LAND_PTS
         
         CALL alloc_cable_types()
         
         DO i=1,land_pts
            DO j=1,ntiles
               
               IF( um1%TILE_FRAC(i,j) .GT. 0.0 ) THEN 
                     um1%L_TILE_PTS(i,j) = .TRUE.
                  !jhan:can set veg%iveg from  here ?
                  tile_index_mp(i,j) = j 
               ENDIF
            
            ENDDO
         ENDDO
      
      ENDIF
         
      !jhan: turn this off until implementation finalised
      !--- initialize latitude/longitude & mapping IF required
      !if ( first_call ) & 
      !   call initialize_maps(latitude,longitude, tile_index_mp)



      !--- read in soil (and veg) parameters 
      IF(first_call)                                                        & 
         CALL  get_type_parameters(logn,vegparmnew)

      !--- initialize veg   
      CALL initialize_veg( canht_ft, lai_ft ) 
 
      !--- initialize soil
      CALL initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt,      &
                            smvccl, albsoil, tsoil_tile, sthu, sthu_tile,   &
                            dzsoil ) 
        
      !--- initialize roughness
      CALL initialize_roughness( z1_tq, z1_uv, kblum_veg%htveg ) 
       
      !--- initialize soilsnow
      CALL initialize_soilsnow( smvcst, tsoil_tile, sthf_tile, smcl_tile,   &
                                snow_tile, snow_rho1l, snage_tile,          &
                                isnow_flg3l, snow_rho3l, snow_cond,         &
                                snow_depth3l, snow_mass3l, snow_tmp3l,      &
                                fland, sin_theta_latitude ) 

      !--- initialize canopy   
      CALL initialize_canopy(CANOPY_TILE)
 
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
         CALL init_respiration(NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC)
      ENDIF      

      ! Lestevens 28 Sept 2012 - Initialize CASA-CNP here
         if (l_casacnp) then
           if (knode_gl==0) then
             print *, '  '; print *, 'CASA_log:'
             print *, '  Calling CasaCNP - Initialise '
             print *, '  l_casacnp = ',l_casacnp
             print *, 'End CASA_log:'; print *, '  '
           endif
      IF( first_call ) THEN
           call init_casacnp(sin_theta_latitude,cpool_tile,npool_tile,&
                             ppool_tile,soil_order,nidep,nifix,pwea,pdust,&
                             wood_hvest_c,wood_hvest_n,wood_hvest_p, &
                             wood_flux_c,wood_flux_n,wood_flux_p, &
                             wresp_c,wresp_n,wresp_p,thinning, &
                             GLAI,PHENPHASE,PREV_YR_SFRAC,iday)
      ENDIF      
      ! Lestevens 4 sept 2018 - init Ndep every tstep for ancil updates
           call casa_ndep_pk(nidep)

         endif

      IF( first_call ) THEN
         CALL dealloc_vegin_soilin()
         first_call = .FALSE. 
      ENDIF      
      
END SUBROUTINE interface_UM_data
                                   
!============================================================================
!============================================================================
!============================================================================

SUBROUTINE assign_um_basics_to_um1( row_length, rows, land_pts, ntiles,     &
                                    npft, sm_levels, timestep, latitude,    &
                                    longitude, land_index, tile_frac,       &
                                    tile_pts, tile_index, l_tile_pts,       &
                                    rho_water  )
   USE cable_um_tech_mod,   ONLY : um1
   USE cable_common_module, ONLY : cable_user

   INTEGER, INTENT(IN) :: row_length, rows, land_pts, ntiles, npft, sm_levels
   INTEGER, INTENT(IN) :: timestep 
   REAL, INTENT(IN), DIMENSION(row_length,rows) :: latitude, longitude 
   REAL,INTENT(IN):: rho_water
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

END SUBROUTINE assign_um_basics_to_um1

!============================================================================
!============================================================================
!============================================================================

END MODULE cable_um_init_mod
