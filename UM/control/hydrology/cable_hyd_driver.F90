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
! Purpose: Converts selected CABLE hydrology variables to UM variables for 
!          UM hydrology code
!
! Called from: UM code hydrol
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Written for CABLE v1.8. No change for CABLE v2.0
!          In future could be combined with standard unpacking of cable 
!          variables at end of implicit call
!
!          2018 WB_LAKE extraction removed and placed in a separate routine
!          as needed for interation with ACCESS rivers' scheme
!
! ==============================================================================
module cable_hyd_driv_mod
  
contains

SUBROUTINE cable_hyd_driver( land_pts, ntiles, lying_snow, SNOW_TILE, SURF_ROFF,&
                             SUB_SURF_ROFF, TOT_TFALL )

  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime
  USE cable_data_module, ONLY : cable
  !from old version
  !USE cable_common_module!, only : cable_runtime, cable_user
  USE cable_data_module,   ONLY : PHYS, OTHER
  USE cable_um_tech_mod, only : um1, ssnow, canopy, veg
  USE cable_decs_mod, ONLY : L_tile_pts

  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir
  USE cable_def_types_mod, ONLY : mp !only need for fprint here
  
  implicit none

  !___ re-decl input args

  integer :: land_pts, ntiles
  
  REAL, INTENT(OUT), DIMENSION(LAND_PTS,NTILES) ::                    &
    SNOW_TILE   ! IN Lying snow on tiles (kg/m2)        

  REAL, INTENT(OUT), DIMENSION(LAND_PTS) ::                               &
    LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
    SUB_SURF_ROFF, & !
    SURF_ROFF,     & !
    TOT_TFALL        !

  !___ local vars

  REAL, DIMENSION(LAND_PTS,NTILES) ::                                 &
    SURF_CAB_ROFF,    &
    TOT_TFALL_TILE                

  REAL :: miss =0. 
  REAL, POINTER :: TFRZ
 
  ! std template args 
  character(len=*), parameter :: subr_name = "cable_explicit_main"

# include "../../../core/utils/diag/cable_fprint.txt"
  
  !-------- Unique subroutine body -----------
 
  TFRZ => PHYS%TFRZ
   
  SNOW_TILE= UNPACK(ssnow%snowd, L_TILE_PTS, miss) 

  LYING_SNOW = SUM(um1%TILE_FRAC * SNOW_TILE,2) !gridbox snow mass

  SURF_CAB_ROFF  = UNPACK(ssnow%rnof1, L_TILE_PTS, miss)
  SURF_ROFF      = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)

  SURF_CAB_ROFF  = UNPACK(ssnow%rnof2, L_TILE_PTS, miss)
  SUB_SURF_ROFF  = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)

  ! %through is /dels in UM app. for STASH output  
  canopy%through = canopy%through / kwidth_gl
  TOT_TFALL_TILE = UNPACK(canopy%through, L_TILE_PTS, miss)
  TOT_TFALL      = SUM(um1%TILE_FRAC * TOT_TFALL_TILE,2)

  !-------- End Unique subroutine body -----------
  
  fprintf_dir=trim(fprintf_dir_root)//trim(unique_subdir)//"/"
  if(L_cable_fprint) then 
    !basics to std output stream
    if (knode_gl == 0 .and. ktau_gl == 1)  call cable_fprintf(subr_name, .true.) 
    !more detailed output
    vname=trim(subr_name//'_')
    call cable_fprintf( cDiag00, vname, knode_gl, ktau_gl, .true. )
  endif

  if(L_cable_Pyfprint .and. ktau_gl == 1) then 
    !vname='latitude'; dimx=mp
    !call cable_Pyfprintf( cDiag1, vname, cable%lat, dimx, .true.)
  endif

return

END SUBROUTINE cable_hyd_driver

SUBROUTINE cable_lakesrivers(TOT_WB_LAKE)
    
    USE cable_um_tech_mod, ONLY : um1, ssnow
    USE cable_decs_mod, ONLY : L_tile_pts
    
    IMPLICIT NONE
    
    !routine extracts daily integrated ssnow%totwblake - water added to keep 
    !lake tiles saturated - and grid cell averages (over land fraction)
    !for use in river flow scaling routines
    
    REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS) :: TOT_WB_LAKE
    
    !working variables
    REAL :: miss = 0.
    REAL, DIMENSION(um1%LAND_PTS, um1%ntiles) :: TOT_WB_LAKE_TILE
    
    TOT_WB_LAKE_TILE = UNPACK(ssnow%totwblake, um1%L_TILE_PTS, miss)
    TOT_WB_LAKE = SUM(um1%TILE_FRAC * TOT_WB_LAKE_TILE,2)
      
    !zero the current integration
    ssnow%totwblake = 0.

END SUBROUTINE cable_lakesrivers

End module cable_hyd_driv_mod



