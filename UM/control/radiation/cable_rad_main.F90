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
! Purpose:
!
! Called from: JULES: surf_couple_radiation, tile_albedo_cable
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_rad_main_mod
  
contains

SUBROUTINE cable_rad_main( mype, timestep_number, cycleno, numcycles,          &
                           row_length, rows, land_pts, ntiles, sm_levels,      &
                           tile_frac, surf_down_sw, cosine_zenith_angle,       &
                           snow_tile, soil_alb, land_albedo, alb_surft,        &
                           land_alb )
  !subrs
  USE cable_rad_driv_mod, ONLY : cable_rad_driver
  
  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir
  
  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl, knode_gl

  !data
  USE atm_fields_real_mod, ONLY : snow_temp_cable, snow_avg_rho_cable, soil_temp_cable,&      
                                snow_flg_cable 
  USE cable_decs_mod, ONLY : sw_down_diag

  implicit none
  !--- IN ARGS FROM surf_couple_radiation() ------------------------------------
  
  integer :: mype, timestep_number, cycleno, numcycles
  integer :: row_length, rows, land_pts, ntiles, sm_levels ! grid
  real :: surf_down_sw(row_length,rows,4) ! 4-band ShortWave forcing
  real :: tile_frac(land_pts,ntiles)      ! surface type fraction 
  real :: cosine_zenith_angle(row_length,rows)        ! cosine_zenith_angle          
  real :: snow_tile(land_pts,ntiles)     ! formerly snow tile        
  real :: soil_alb(land_pts)              ! Snow-free, bare soil albedo: 
  real :: land_albedo(row_length,rows,4) 
  real :: alb_surft(Land_pts,ntiles,4)
  real :: land_alb(row_length,rows)         ! Mean land_albedo

  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_rad_main"
  logical, save :: first_call = .true.
  integer,  DIMENSION(land_pts, ntiles) :: isnow_flg_cable
  
# include "../../../core/utils/diag/cable_fprint.txt"
  
  !-------- Unique subroutine body -----------
  
  !--- initialize cable_runtime% switches 
  cable_runtime%um =          .TRUE.
  cable_runtime%um_radiation= .TRUE.
  
  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  isnow_flg_cable = int(snow_flg_cable)
  if( .NOT. allocated(sw_down_diag) )allocate( sw_down_diag(row_length, rows, 4) ) 
  if(mype==0) then
    write (6, *) "CABLE_LSM:Subr: ", subr_name,  "calling _rad_driver: "
  endif
     
  call cable_rad_driver( cycleno, row_length, rows, land_pts, ntiles,          &
                         sm_levels, tile_frac, snow_temp_cable,                &
                         snow_avg_rho_cable, soil_temp_cable,                  &
                         isnow_flg_cable, surf_down_sw, cosine_zenith_angle,   &
                         snow_tile, soil_alb, land_albedo, alb_surft, land_alb )
  
  
  snow_flg_cable = real(isnow_flg_cable)
  
  first_call = .false.        

  cable_runtime%um_radiation= .FALSE.
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

End subroutine cable_rad_main
  
End module cable_rad_main_mod











































