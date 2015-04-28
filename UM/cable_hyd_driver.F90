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
!
! ==============================================================================

SUBROUTINE cable_hyd_driver( SNOW_TILE, LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF,  &
                             TOT_TFALL, WB_LAKE )

   USE cable_data_module,   ONLY : PHYS, OTHER
   USE cable_common_module!, only : cable_runtime, cable_user
   USE cable_um_tech_mod, only : um1, ssnow, canopy, veg
   IMPLICIT NONE

   REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS,um1%NTILES) ::                    &
      SNOW_TILE   ! IN Lying snow on tiles (kg/m2)        

   REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS) ::                               &
      LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
      SUB_SURF_ROFF, & !
      SURF_ROFF,     & !
      TOT_TFALL        !

   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
      SURF_CAB_ROFF,    &
      TOT_TFALL_TILE                

   ! Lestevens 25sep13 - water balance fix for lakes
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      WB_LAKE         ! unpack CABLE wb_lake

   REAL :: miss =0. 
   REAL, POINTER :: TFRZ
      
      TFRZ => PHYS%TFRZ
   
      SNOW_TILE= UNPACK(ssnow%snowd, um1%L_TILE_PTS, miss) 
      LYING_SNOW = SUM(um1%TILE_FRAC * SNOW_TILE,2) !gridbox snow mass

      SURF_CAB_ROFF  = UNPACK(ssnow%rnof1, um1%L_TILE_PTS, miss)
      SURF_ROFF      = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)
      
      SURF_CAB_ROFF  = UNPACK(ssnow%rnof2, um1%L_TILE_PTS, miss)
      SUB_SURF_ROFF  = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)

      ! %through is /dels in UM app. for STASH output  
      canopy%through = canopy%through / kwidth_gl
      TOT_TFALL_TILE = UNPACK(canopy%through, um1%L_TILE_PTS, miss)
      TOT_TFALL      = SUM(um1%TILE_FRAC * TOT_TFALL_TILE,2)

      ! Lest 25sep13 - wb_lake fix
      WB_LAKE        = UNPACK(ssnow%wb_lake, um1%L_TILE_PTS, miss)
      
END SUBROUTINE cable_hyd_driver
      




