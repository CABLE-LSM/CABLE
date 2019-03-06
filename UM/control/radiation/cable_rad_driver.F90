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

module cable_rad_driv_mod
  
contains

subroutine cable_rad_driver( cycleno, row_length, rows, land_pts, ntiles,      &
                             sm_levels, tile_frac, snow_tmp3L, snow_rho1L,     &
                             Tsoil_tile, isnow_flg3L, surf_down_sw,            &
                             cos_zenith_angle, snow_tile, albsoil,             &
                             land_albedo_cable, alb_tile, land_alb_cable )
 
 !subrs
  USE cable_albedo_module, ONLY : surface_albedo
  USE cable_um_init_subrs_mod, ONLY : update_kblum_radiation,                  &
      um2cable_met_rad, um2cable_lp 
  USE cable_radiation_module, only : init_radiation

  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir

  USE cable_diag_module  

  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl, knode_gl

  !data  
  USE cable_def_types_mod, ONLY : mp
  USE cable_um_tech_mod,   ONLY : kblum_rad, um1, soil, ssnow, rad, veg,      &
                                  met, canopy, basic_diag
  USE cable_decs_mod, ONLY : L_tile_pts, sw_down_diag

  implicit none

  !___ re-decl input args
  integer :: mype, timestep_number, cycleno, numcycles
  integer :: row_length, rows, land_pts, ntiles, sm_levels ! grid
  real :: surf_down_sw(row_length,rows,4) ! 4-band ShortWave forcing
  real :: tile_frac(land_pts,ntiles)      ! surface type fraction 
  real :: cos_zenith_angle(row_length,rows)        ! cosine_zenith_angle          
  real :: snow_tile(land_pts,ntiles)     ! formerly snow tile        
  real :: albsoil(land_pts)              ! Snow-free, bare soil albedo: 
  real :: land_albedo_cable(row_length,rows,4) 
  real :: alb_tile(Land_pts,ntiles,4)
  real :: land_alb_cable(row_length,rows)         ! Mean land_albedo
  integer :: isnow_flg3L(land_pts,ntiles)      ! surface type fraction 
 
  REAL, DIMENSION(LAND_PTS,NTILES) ::                                 &
     SNOW_RHO1L             ! snow cover in the ice tile.
  
  REAL, DIMENSION( LAND_PTS, NTILES, SM_LEVELS ) ::               &
     TSOIL_TILE  ! Mean snow density  (or 1 layer)

  REAL, DIMENSION( LAND_PTS, NTILES, 3 ) ::                           &
     SNOW_TMP3L  ! Snow temperature (3 layer)
                                                                  
  !___ local vars
  INTEGER :: i,J,N,K,L,M
  REAL :: miss = 0.0
  LOGICAL :: skip =.TRUE. 
   
  real :: land_alb_cable_tile(land_pts,ntiles)         ! Mean land_albedo
  REAL :: rad_vis(mp), rad_nir(mp), met_fsd_tot_rel(mp), rad_albedo_tot(mp) 

  ! std template args 
  character(len=*), parameter :: subr_name = "cable_rad_driver"
  
# include "../../../core/utils/diag/cable_fprint.txt"
  
  !-------- Unique subroutine body -----------
    
  !     **** surf_down_sw is from the previous time step  ****
  !--- re-set UM rad. forcings to suit CABLE. also called in explicit call to 
  !--- CABLE from subr cable_um_expl_update() 
  CALL update_kblum_radiation( sw_down_diag, cos_zenith_angle, surf_down_sw )
  
  !--- set met. and rad. forcings to CABLE. also called in explicit call to 
  !--- CABLE from subr update_explicit_vars() 
  !--- subr.  um2cable_met_rad_alb() USES CABLE types met%, rad%, soil%
  !--- and kblum%rad. calculated in  update_kblum_radiation() above 
  CALL um2cable_met_rad( cos_zenith_angle)
  
  CALL init_radiation(met, rad, veg, canopy)

  !--- CABLE soil albedo forcing
  CALL um2cable_lp( albsoil, albsoil, soil%albsoil(:,1), soil%isoilm, skip )
  !-------------------------------------------------------------------
  
  !At present only single value is used for each land point 
  soil%albsoil(:,2) = REAL(0) 
  soil%albsoil(:,3) = REAL(0) 

  ! soil%albsoil should be set to geograpically explicit data for 
  ! snow free soil albedo in VIS and NIR  
  ! get the latest surface condtions
  ssnow%snowd  =     PACK( SNOW_TILE, L_TILE_PTS )
  ssnow%ssdnn  =     PACK( SNOW_RHO1L, L_TILE_PTS )
  ssnow%isflag =     PACK( ISNOW_FLG3L, L_TILE_PTS )
  ssnow%tggsn(:,1) = PACK( SNOW_TMP3L(:,:,1), L_TILE_PTS )
  ssnow%tgg(:,1) =   PACK( TSOIL_TILE(:,:,1), L_TILE_PTS )


  CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)

  ! only for land points, at present do not have a method for treating 
  ! mixed land/sea or land/seaice points as yet.
  ALB_TILE(:,:,1) = UNPACK(rad%reffbm(:,1),L_TILE_PTS, miss)
  ALB_TILE(:,:,2) = UNPACK(rad%reffdf(:,1),L_TILE_PTS, miss)
  ALB_TILE(:,:,3) = UNPACK(rad%reffbm(:,2),L_TILE_PTS, miss)
  ALB_TILE(:,:,4) = UNPACK(rad%reffdf(:,2),L_TILE_PTS, miss)

  rad_vis = ( (1.0-rad%fbeam(:,1) ) * rad%reffdf(:,1) + rad%fbeam(:,1) *   &
            rad%reffbm(:,1) ) 
  rad_nir = ( (1.0-rad%fbeam(:,2) ) * rad%reffdf(:,2) + rad%fbeam(:,2) *   &
            rad%reffbm(:,2)) 

  met_fsd_tot_rel = met%fsd(:,1) / MAX( 0.1, ( met%fsd(:,1)+met%fsd(:,2) ) )
  rad_albedo_tot = met_fsd_tot_rel  * rad_vis                              &
                   + ( 1.- met_fsd_tot_rel ) * rad_nir
  LAND_ALB_CABLE=0
  LAND_ALBEDO_CABLE=0

  LAND_ALB_CABLE_TILE = UNPACK( rad_albedo_tot, L_TILE_PTS, miss )

  DO N=1,NTILES
    DO K=1,um1%TILE_PTS(N)
      L = um1%TILE_INDEX(K,N)
      J=(um1%LAND_INDEX(L)-1)/ROW_LENGTH + 1
      I = um1%LAND_INDEX(L) - (J-1)*ROW_LENGTH
      
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

    !-------- End Unique subroutine body -----------

  fprintf_dir=trim(fprintf_dir_root)//trim(unique_subdir)//"/"
  if(L_cable_fprint) then 
    !basics to std output stream
    if (knode_gl == 0 .and. ktau_gl == 1)  call cable_fprintf(subr_name, .true.) 
    !more detailed output
    vname=trim(subr_name//'_')
    call cable_fprintf( cDiag00, vname, knode_gl, ktau_gl, .true. )
  endif

  if(L_cable_Pyfprint) then 
    !vname='tscrn'; dimx=land_pts
    !call cable_Pyfprintf( cDiag2, vname, t1p5m, dimx, .true.)
    !vname='tstar_tile'; dimx=size(tstar_tile,1); dimy=size(tstar_tile,2)
    !call cable_Pyfprintf( cDiag4, vname, tstar_tile, dimx, dimy, .true.)
  endif

return
 
END SUBROUTINE cable_rad_driver
     
End module cable_rad_driv_mod

