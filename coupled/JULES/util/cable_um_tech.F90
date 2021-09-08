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
! Purpose:  Module declares and allocates um1% type. Dimensions etc from parent model
!           which are used frequently. This legacy and should be removed in favour of
!           using vars directly passedfrom the parent model
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Rewrite of code from v1.8 (ACCESS1.3)
!          soil_snow_type now ssnow (instead of ssoil)
! ==============================================================================

MODULE cable_um_tech_mod
   
USE cable_def_types_mod
IMPLICIT NONE
  
TYPE um_dimensions 
  INTEGER :: row_length, rows, land_pts, ntiles, npft, sm_levels, timestep 
  REAL    :: rho_water,rho_ice
  INTEGER, ALLOCATABLE :: tile_pts(:), land_index(:)
  INTEGER, ALLOCATABLE :: tile_index(:,:)
  REAL, ALLOCATABLE :: tile_frac(:,:)
  REAL, ALLOCATABLE :: latitude(:,:), longitude(:,:)
  LOGICAL, ALLOCATABLE :: l_tile_pts(:,:) 
END TYPE um_dimensions 

TYPE(um_dimensions)    :: um1

REAL,ALLOCATABLE, DIMENSION(:) :: conv_rain_prevstep, conv_snow_prevstep 

CONTAINS

SUBROUTINE alloc_um_interface_types( row_length, rows, land_pts, ntiles,      &
                                     sm_levels )
USE cable_common_module, ONLY: cable_runtime, cable_user
      
INTEGER,INTENT(IN) :: row_length, rows, land_pts, ntiles, sm_levels   

IF ( .NOT. ALLOCATED ( um1%land_index ) )  ALLOCATE( um1%land_index(land_pts) )
IF ( .NOT. ALLOCATED ( um1%tile_pts ) )    ALLOCATE( um1%tile_pts(ntiles) )
IF ( .NOT. ALLOCATED ( um1%tile_frac ) )   ALLOCATE( um1%tile_frac(land_pts, ntiles) )
IF ( .NOT. ALLOCATED ( um1%tile_index ) )  ALLOCATE( um1%tile_index(land_pts, ntiles) )
IF ( .NOT. ALLOCATED ( um1%latitude ) )    ALLOCATE( um1%latitude(row_length, rows) )
IF ( .NOT. ALLOCATED ( um1%longitude ) )   ALLOCATE( um1%longitude(row_length, rows) )
IF ( .NOT. ALLOCATED ( um1%l_tile_pts ) )  ALLOCATE( um1%l_tile_pts(land_pts, ntiles) ) 
         
END SUBROUTINE alloc_um_interface_types 

END MODULE cable_um_tech_mod
