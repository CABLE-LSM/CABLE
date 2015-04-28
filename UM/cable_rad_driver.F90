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
! Purpose: Pass CABLE albedo to UM variable for UM radiation scheme
!
! Called from: UM code glue_rad_ctl
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
!
! ==============================================================================

SUBROUTINE cable_rad_driver(                                                   &
                             ! IN atmospheric forcing
                             surf_down_sw, cos_zenith_angle,                   &
                             ! IN soil/vegetation/land surface data :
                             SNOW_TILE, SNOW_TMP3L, SNOW_RHO1L, TSOIL_TILE,    &
                             ISNOW_FLG3L, ALBSOIL,                             &
                             ! OUT
                             LAND_ALBEDO_CABLE, ALB_TILE, LAND_ALB_CABLE ) 

   USE cable_def_types_mod, ONLY : mp
   USE cable_albedo_module, ONLY : surface_albedo
   USE cable_um_tech_mod,   ONLY : kblum_rad, um1, soil, ssnow, rad, veg,      &
                                   met, canopy
   USE cable_um_init_subrs_mod, ONLY : update_kblum_radiation,  um2cable_met_rad,  &
                                   um2cable_lp 
   USE cable_common_module, ONLY : cable_runtime, cable_user
   
   IMPLICIT NONE                     

   INTEGER, DIMENSION(um1%LAND_PTS,um1%NTILES) :: isnow_flg3l    
   
   REAL :: ALBSOIL(um1%LAND_PTS)          !          &! IN soil albedo 
   
   REAL, DIMENSION(um1%row_length,um1%rows) ::                                 &
      LAND_ALB_CABLE,      & ! Land albedo calculated by Cable
      SW_DOWN,             & ! Surface downward SW radiation (W/m2).
      cos_zenith_angle

   REAL, DIMENSION(um1%row_length,um1%rows,4) ::                               &
      LAND_ALBEDO_CABLE, & ! Land albedo calculated by Cable NIR/VIS/Beam/Diffuse
      surf_down_sw         ! IN Surface downward SW radiation

   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
      LAND_ALB_CABLE_TILE,  & ! Land albedo calculated by Cable
      SNOW_TILE,            & ! IN Lying snow on tiles (kg/m2) 
      SNOW_RHO1L             ! snow cover in the ice tile.
   
   REAL, DIMENSION( um1%LAND_PTS, um1%NTILES, 4 ) ::                           &
      ALB_TILE    ! Land albedo calculated by Cable
      
   REAL, DIMENSION( um1%LAND_PTS, um1%NTILES, um1%SM_LEVELS ) ::               &
      TSOIL_TILE  ! Mean snow density  (or 1 layer)

   REAL, DIMENSION( um1%LAND_PTS, um1%NTILES, 3 ) ::                           &
      SNOW_TMP3L  ! Snow temperature (3 layer)
                                                                   
   INTEGER :: i,J,N,K,L
   REAL :: miss = 0.0
   LOGICAL :: skip =.TRUE. 
   
   REAL :: rad_vis(mp), rad_nir(mp), met_fsd_tot_rel(mp), rad_albedo_tot(mp) 

      !jhan:check that these are reset after call done
      cable_runtime%um_radiation= .TRUE.
      
      !     **** surf_down_sw is from the previous time step  ****
      !--- re-set UM rad. forcings to suit CABLE. also called in explicit call to 
      !--- CABLE from subr cable_um_expl_update() 
      CALL update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
   
      !--- set met. and rad. forcings to CABLE. also called in explicit call to 
      !--- CABLE from subr update_explicit_vars() 
      !--- subr.  um2cable_met_rad_alb() USES CABLE types met%, rad%, soil%
      !--- and kblum%rad. calculated in  update_kblum_radiation() above 
      CALL um2cable_met_rad( cos_zenith_angle)
      
      !--- CABLE soil albedo forcing
      CALL um2cable_lp( albsoil, albsoil, soil%albsoil(:,1), soil%isoilm, skip )
      !-------------------------------------------------------------------
  
      !At present only single value is used for each land point 
      soil%albsoil(:,2) = REAL(0) 
      soil%albsoil(:,3) = REAL(0) 

      ! soil%albsoil should be set to geograpically explicit data for 
      ! snow free soil albedo in VIS and NIR  
      ! get the latest surface condtions
      ssnow%snowd  =     PACK( SNOW_TILE, um1%L_TILE_PTS )
      ssnow%ssdnn  =     PACK( SNOW_RHO1L, um1%L_TILE_PTS )
      ssnow%isflag =     PACK( ISNOW_FLG3L, um1%L_TILE_PTS )
      ssnow%tggsn(:,1) = PACK( SNOW_TMP3L(:,:,1), um1%L_TILE_PTS )
      ssnow%tgg(:,1) =   PACK( TSOIL_TILE(:,:,1), um1%L_TILE_PTS )


      CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)

      ! only for land points, at present do not have a method for treating 
      ! mixed land/sea or land/seaice points as yet.
      ALB_TILE(:,:,1) = UNPACK(rad%reffbm(:,1),um1%L_TILE_PTS, miss)
      ALB_TILE(:,:,2) = UNPACK(rad%reffdf(:,1),um1%L_TILE_PTS, miss)
      ALB_TILE(:,:,3) = UNPACK(rad%reffbm(:,2),um1%L_TILE_PTS, miss)
      ALB_TILE(:,:,4) = UNPACK(rad%reffdf(:,2),um1%L_TILE_PTS, miss)

      rad_vis = ( (1.0-rad%fbeam(:,1) ) * rad%reffdf(:,1) + rad%fbeam(:,1) *   &
                rad%reffbm(:,1) ) 
      rad_nir = ( (1.0-rad%fbeam(:,2) ) * rad%reffdf(:,2) + rad%fbeam(:,2) *   &
                rad%reffbm(:,2)) 

      met_fsd_tot_rel = met%fsd(:,1) / MAX( 0.1, ( met%fsd(:,1)+met%fsd(:,2) ) )
      rad_albedo_tot = met_fsd_tot_rel  * rad_vis                              &
                       + ( 1.- met_fsd_tot_rel ) * rad_nir

      LAND_ALBEDO_CABLE =0.
      LAND_ALB_CABLE =0.
      LAND_ALB_CABLE_TILE = UNPACK( rad_albedo_tot, um1%L_TILE_PTS, miss )

      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            J=(um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
            I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
            
            ! direct beam visible
            LAND_ALBEDO_CABLE(I,J,1) = LAND_ALBEDO_CABLE(I,J,1) +              &
                                       um1%TILE_FRAC(L,N)*ALB_TILE(L,N,1)

            ! diffuse beam visible
            LAND_ALBEDO_CABLE(I,J,2) = LAND_ALBEDO_CABLE(I,J,2) +              &
                                       um1%TILE_FRAC(L,N)*ALB_TILE(L,N,2)

            ! direct beam nearinfrared 
            LAND_ALBEDO_CABLE(I,J,3) = LAND_ALBEDO_CABLE(I,J,3) +              &
                                       um1%TILE_FRAC(L,N)*ALB_TILE(L,N,3)

            ! diffuse beam nearinfrared
            LAND_ALBEDO_CABLE(I,J,4) = LAND_ALBEDO_CABLE(I,J,4) +              &
                                       um1%TILE_FRAC(L,N)*ALB_TILE(L,N,4)
            LAND_ALB_CABLE(I,J) = LAND_ALB_CABLE(I,J) +                        &
                                um1%TILE_FRAC(L,N)*LAND_ALB_CABLE_TILE(L,N)
         ENDDO
      ENDDO

      cable_runtime%um_radiation= .FALSE.

END SUBROUTINE cable_rad_driver

     

