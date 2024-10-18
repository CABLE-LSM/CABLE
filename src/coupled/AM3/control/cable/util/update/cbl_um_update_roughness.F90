MODULE cbl_um_update_roughness_mod
   
IMPLICIT NONE

CONTAINS
     
SUBROUTINE update_roughness( row_length, rows, mp, land_pts, nsurft, npft,     &
                             lai_thresh, surft_pts, surft_index, land_index,   &
                             L_tile_pts, z1_tq, z1_uv, HGT_pft_um, HGT_pft_cbl,&
                             rough, veg )

USE cable_pack_mod,       ONLY: cable_pack_rr
USE cable_def_types_mod,  ONLY: roughness_type, veg_parameter_type

IMPLICIT NONE
                                                               
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
INTEGER, INTENT(IN) :: npft 
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
REAL,    INTENT(IN) :: lai_thresh 
INTEGER, INTENT(IN) :: surft_pts(nsurft) ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: land_index(land_pts)  ! tangled cell index of land_pt
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)  ! TRUE if active tile
REAL,    INTENT(IN) :: z1_tq(row_length, rows)
REAL,    INTENT(IN) :: z1_uv(row_length, rows)
REAL,    INTENT(IN) :: HGT_pft_um(land_pts,nsurft)
REAL,    INTENT(IN) :: HGT_pft_cbl(mp)

TYPE(roughness_type),     INTENT(OUT) :: rough
TYPE(veg_parameter_type), INTENT(IN)  :: veg        ! vegetation parameters

!local vars 
INTEGER :: i,j,k,L,n
REAL    :: jhruff(land_pts,nsurft)
REAL    :: jhwork(land_pts,nsurft)

!--- CABLE roughness type forcings
CALL cable_pack_rr( rough%za_tq, z1_tq, mp, l_tile_pts, row_length,            &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( rough%za_uv, z1_uv, mp, l_tile_pts, row_length,            &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

!Veg height changes seasonally in MOSES hence no updates here due to snow
jhwork = 0.
DO n=1,nsurft
  DO k=1,surft_pts(n)
    l = surft_index(k,n)
    jhwork(l,n) = MAX( .01, HGT_pft_um(L,N) )
  ENDDO
ENDDO

jhruff= 0.01 
DO l=1,land_pts
  DO n=1,nsurft     
    IF( jhruff(L,N) .lt. jhwork(l,n)) THEN
      jhruff(L,:) =  jhwork(l,n)
    END IF
  ENDDO
ENDDO

! CM2 set hruff from veg%hc - will require review with POP implementation       
rough%hruff= MAX( 0.01, HGT_pft_cbl ) 
rough%hruff_grmx = PACK(jhruff, l_tile_pts) 

RETURN
END SUBROUTINE update_roughness

END MODULE cbl_um_update_roughness_mod
