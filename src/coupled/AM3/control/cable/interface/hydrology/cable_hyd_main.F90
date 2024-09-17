MODULE Cable_hyd_main_mod
  
CONTAINS

SUBROUTINE cable_hyd_main(land_pts, nsurft, tile_frac, timestep_len,                    &
                     lying_snow, snow_surft,                                            &
                     snow_melt_gb, surf_roff, sub_surf_roff,                            &
                     tot_tfall, work_snow_surft,                                        &
                     melt_surft, work_lying_snow, work_surf_roff,                       &
                     work_sub_surf_roff, work_tot_tfall)

! This routine passes some of CABLE's hydrology variables to JULES for purpose
! of output.
                    
USE cable_common_module, ONLY : cable_runtime
USE cable_common_module, ONLY : knode_gl

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nsurft
REAL, INTENT(IN) :: timestep_len
REAL, INTENT(IN) :: tile_frac(land_pts, nsurft)
REAL, INTENT(IN) :: work_snow_surft(land_pts,nsurft)
REAL, INTENT(IN) :: melt_surft(land_pts,nsurft)      !  tile snow melt (from implicit)
REAL, INTENT(IN) :: work_lying_snow(land_pts)        !  Gridbox snowmass (kg/m2)     
REAL, INTENT(IN) :: work_sub_surf_roff(land_pts)     !  Sub-surface runoff (kg/m2/s).
REAL, INTENT(IN) :: work_surf_roff(land_pts)         !  Surface runoff (kg/m2/s).
REAL, INTENT(IN) :: work_tot_tfall(land_pts)         !  Total throughfall (kg/m2/s).

REAL, INTENT(OUT) :: snow_surft(land_pts,nsurft)
REAL, INTENT(OUT) :: snow_melt_gb(land_pts)      ! OUT gridbox snow melt (kg/m2/s)
REAL, INTENT(OUT) :: lying_snow(land_pts)        ! OUT Gridbox snowmass (kg/m2)     
REAL, INTENT(OUT) :: sub_surf_roff(land_pts)     ! OUT Sub-surface runoff (kg/m2/s).
REAL, INTENT(OUT) :: surf_roff(land_pts)         ! OUT Surface runoff (kg/m2/s).
REAL, INTENT(OUT) :: tot_tfall(land_pts)         ! OUT Total throughfall (kg/m2/s).
   
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_hyd_main"
LOGICAL, SAVE :: zero_points_warning = .true.
 
IF( land_pts  ==  0 ) THEN
  IF( zero_points_warning ) THEN
    WRITE(6,*) "Reached CABLE ", subr_name,                                      &
               " even though zero land_points on processor ", knode_gl 
  END IF
  zero_points_warning = .FALSE. 
  RETURN
END IF

!--- initialize cable_runtime% switches 
cable_runtime%um =          .TRUE.
cable_runtime%um_hydrology =.TRUE.

!jhan:This is all very dubious - mixing up depth and mass
!snow_surft     = work_snow_surft !this is overwriting with zero 
lying_snow     = work_lying_snow 

surf_roff      = work_surf_roff      
sub_surf_roff  = work_sub_surf_roff  
tot_tfall      = work_tot_tfall 

!CM3 new variables
! melt_surft unpacked in cable_implicit from ssnow%smelt
snow_melt_gb    = SUM(tile_frac * melt_surft,2)

cable_runtime%um_hydrology =.FALSE.
  
RETURN

END SUBROUTINE cable_hyd_main

!!jhan:we'll have to organize ho we will get the right thigs here 
! For CM2 we have to pass mp & %totwblake and TILE_FRAC down through argument list from *extras*
SUBROUTINE cable_lakesrivers(land_pts, nsurft, mp, totwblake, TILE_FRAC, &
                    TOT_WB_LAKE, L_tile_pts)
       
    IMPLICIT NONE
    
    !routine extracts daily integrated ssnow%totwblake - water added to keep 
    !lake tiles saturated - and grid cell averages (over land fraction)
    !for use in river flow scaling routines

    !This routine is called from riv_intcl-riv_ic1a
    
    INTEGER, INTENT(IN) :: land_pts, nsurft, mp
    LOGICAL :: l_tile_pts(land_pts, nsurft)
    REAL, INTENT(INOUT), DIMENSION(mp) :: totwblake
    REAL, INTENT(IN) :: TILE_FRAC(land_pts, nsurft)
    REAL, INTENT(OUT), DIMENSION(land_pts) :: TOT_WB_LAKE
    
    !working variables
    REAL :: miss = 0.
    REAL, DIMENSION(LAND_PTS, nsurft) :: TOT_WB_LAKE_TILE
    
    !CM2 era code
    !TOT_WB_LAKE_TILE = UNPACK(ssnow%totwblake, L_TILE_PTS, miss)

    !updated for CM3
    TOT_WB_LAKE_TILE = UNPACK(totwblake, L_TILE_PTS, miss)
    TOT_WB_LAKE = SUM(TILE_FRAC * TOT_WB_LAKE_TILE,2)
      
    !zero the current integration
    totwblake = 0.

END SUBROUTINE cable_lakesrivers

 
End module cable_hyd_main_mod











































