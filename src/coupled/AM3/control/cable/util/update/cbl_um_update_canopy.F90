MODULE cbl_um_update_canopy_mod
   
IMPLICIT NONE

PUBLIC update_canopy

CONTAINS
     
SUBROUTINE update_canopy( mp, land_pts, nsurft, L_tile_pts, canopy_tile,       &
                          reducedLAIdue2snow, canopy )

USE cable_def_types_mod, ONLY: canopy_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: land_pts                      ! # land points 
INTEGER, INTENT(IN) :: nsurft                        ! # tiles 
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)
REAL,    INTENT(IN) :: canopy_tile(land_pts, nsurft)
REAL,    INTENT(IN) :: reducedLAIdue2snow(mp)
TYPE(canopy_type), INTENT(OUT) :: canopy
   
!---set canopy storage (already in dim(land_pts,ntiles) ) 
canopy%cansto    = PACK(CANOPY_TILE, l_tile_pts)
canopy%oldcansto = canopy%cansto
canopy%vlaiw     = reducedLAIdue2snow

RETURN
END SUBROUTINE update_canopy

END MODULE cbl_um_update_canopy_mod




