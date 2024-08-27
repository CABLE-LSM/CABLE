MODULE cbl_um_update_radiation_mod
   
IMPLICIT NONE

CONTAINS

SUBROUTINE update_radiation( mp, row_length, rows, timestep, land_pts, nsurft, &
                             surft_pts, surft_index, land_index, L_tile_pts,   &
                             latitude, longitude, cos_zenith_angle,            & 
                             sw_down_VIS, sw_down_NIR, beamFrac_VIS,           &
                             beamFrac_NIR, beamFrac_TOT,  lw_down, rad, met )
                                                               
USE cable_pack_mod,             ONLY: cable_pack_rr
USE cable_other_constants_mod,  ONLY: RAD_THRESH 
USE cable_def_types_mod,        ONLY: radiation_type, met_type
IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: row_length        ! # columns in spatial grid
INTEGER, INTENT(IN) :: rows              ! # rows in spatial grid
REAL,    INTENT(IN) :: timestep 
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft            ! # tiles 

INTEGER, INTENT(IN) :: surft_pts(nsurft) ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: land_index(land_pts)  ! tangled cell index of land_pt
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)  ! TRUE if active tile
REAL,    INTENT(IN)    :: cos_zenith_angle(row_length,rows)
REAL,    INTENT(IN)    :: latitude(row_length,rows)
REAL,    INTENT(IN)    :: longitude(row_length,rows)

! radiation streams 
REAL,    INTENT(IN)    :: sw_down_VIS(row_length,rows)
REAL,    INTENT(IN)    :: sw_down_NIR(row_length,rows)
REAL,    INTENT(IN)    :: beamFrac_VIS(row_length,rows)
REAL,    INTENT(IN)    :: beamFrac_NIR(row_length,rows)
REAL,    INTENT(IN)    :: beamFrac_TOT(row_length,rows)
REAL,    INTENT(IN)    :: lw_down(row_length,rows)
  
TYPE(radiation_type), INTENT(OUT) :: rad
TYPE(met_type),       INTENT(OUT) :: met
   
!local vars
INTEGER :: i
     
rad%otrad = rad%trad

CALL cable_pack_rr( met%coszen, cos_zenith_angle, mp, l_tile_pts, row_length,  &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

met%coszen =    max(met%coszen,1e-8)  ! is this really required now

CALL cable_pack_rr( met%fsd(:,1), sw_down_VIS, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( met%fsd(:,2), sw_down_NIR, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( rad%fbeam(:,1), beamFrac_VIS, mp, l_tile_pts, row_length,  &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( rad%fbeam(:,2), beamFrac_NIR, mp, l_tile_pts, row_length,  &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( rad%fbeam(:,3), beamFrac_TOT, mp, l_tile_pts, row_length,  &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )
DO i=1,mp
  IF( met%coszen(i) < RAD_THRESH ) THEN 
    rad%fbeam(i,1) = REAL(0) 
    rad%fbeam(i,2) = REAL(0) 
    rad%fbeam(i,3) = REAL(0) 
  END IF
END DO

CALL cable_pack_rr( met%fld, lw_down, mp, l_tile_pts, row_length,              &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( rad%latitude, latitude, mp, l_tile_pts, row_length,        &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( rad%longitude, longitude, mp, l_tile_pts, row_length,      &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

RETURN
END SUBROUTINE update_radiation

END MODULE cbl_um_update_radiation_mod
