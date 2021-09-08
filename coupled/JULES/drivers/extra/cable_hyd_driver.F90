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
MODULE cable_hyd_driv_mod
  
CONTAINS

SUBROUTINE cable_hyd_driver( land_pts, ntiles, L_tile_pts, lying_snow, snow_tile, surf_roff, &
                             sub_surf_roff, tot_tfall, ssnow, canopy, veg )

  !processor number, timestep number / width, endstep
USE cable_common_module, ONLY: knode_gl, ktau_gl, kwidth_gl, kend_gl
USE cable_common_module, ONLY: cable_runtime
USE cable_um_tech_mod, ONLY: um1 

USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_params_mod,         ONLY: veg_parameter_type

USE cable_common_module, ONLY: cable_runtime, cable_user, l_casacnp,          &
                                l_vcmaxFeedbk, knode_gl, ktau_gl, kend_gl
  
IMPLICIT NONE
TYPE (canopy_type),    INTENT(INOUT) :: canopy
TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg

!___ re-decl input args
INTEGER :: land_pts, ntiles
LOGICAL :: L_tile_pts(land_pts,ntiles)     ! packing mask 
  
REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) ::                              &
  snow_tile   ! IN Lying snow on tiles (kg/m2)        

REAL, INTENT(OUT), DIMENSION(land_pts) ::                                     &
  lying_snow,    & ! OUT Gridbox snowmass (kg/m2)        
  sub_surf_roff, & !
  surf_roff,     & !
  tot_tfall        !

!___ local vars

REAL, DIMENSION(land_pts,ntiles) ::                                           &
  surf_cab_roff,                                                              &
  tot_tfall_tile                

REAL :: miss  = 0.0 
 
! std template args 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_explicit_main"

!-------- Unique subroutine body -----------
 
snow_tile= UNPACK(ssnow%snowd, l_tile_pts, miss) 

lying_snow = SUM(um1%tile_frac * snow_tile,2) !gridbox snow mass

surf_cab_roff  = UNPACK(ssnow%rnof1, l_tile_pts, miss)
surf_roff      = SUM(um1%tile_frac * surf_cab_roff,2)

surf_cab_roff  = UNPACK(ssnow%rnof2, l_tile_pts, miss)
sub_surf_roff  = SUM(um1%tile_frac * surf_cab_roff,2)

! %through is /dels in UM app. for STASH output  
canopy%through = canopy%through / kwidth_gl
tot_tfall_tile = UNPACK(canopy%through, l_tile_pts, miss)
tot_tfall      = SUM(um1%tile_frac * tot_tfall_tile,2)

!-------- End Unique subroutine body -----------
  
RETURN

END SUBROUTINE cable_hyd_driver

!H!SUBROUTINE cable_lakesrivers(TOT_WB_LAKE)
!H!    
!H!    USE cable_um_tech_mod, ONLY : um1, ssnow
!H!    !jhan : thisversion of L_tile_pts is currently inctive
!H!    USE cbl_masks_mod, ONLY : L_tile_pts
!H!    
!H!    IMPLICIT NONE
!H!    
!H!    !routine extracts daily integrated ssnow%totwblake - water added to keep 
!H!    !lake tiles saturated - and grid cell averages (over land fraction)
!H!    !for use in river flow scaling routines
!H!    
!H!    REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS) :: TOT_WB_LAKE
!H!    
!H!    !working variables
!H!    REAL :: miss = 0.
!H!    REAL, DIMENSION(um1%LAND_PTS, um1%ntiles) :: TOT_WB_LAKE_TILE
!H!    
!H!    TOT_WB_LAKE_TILE = UNPACK(ssnow%totwblake, um1%L_TILE_PTS, miss)
!H!    TOT_WB_LAKE = SUM(um1%TILE_FRAC * TOT_WB_LAKE_TILE,2)
!H!      
!H!    !zero the current integration
!H!    ssnow%totwblake = 0.
!H!
!H!END SUBROUTINE cable_lakesrivers

END MODULE cable_hyd_driv_mod



