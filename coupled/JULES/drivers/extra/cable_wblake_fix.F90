MODULE cable_wblake_mod
IMPLICIT NONE

REAL  :: wblake_ratio               ! ratio of wblake/subroff

REAL, ALLOCATABLE, SAVE ::                                                    &
  WBLAKE_cable(:,:), TOT_WBLAKE_cable(:), TOT_SUBRUN_cable(:)

CONTAINS

SUBROUTINE cable_wblake_fix_alloc( land_points, ntiles )
INTEGER :: land_points, ntiles

IF ( .NOT. ALLOCATED( wblake_cable ) )                                        &
  ALLOCATE ( WBLAKE_cable(land_points,ntiles) )

END SUBROUTINE cable_wblake_fix_alloc

END MODULE cable_wblake_mod
