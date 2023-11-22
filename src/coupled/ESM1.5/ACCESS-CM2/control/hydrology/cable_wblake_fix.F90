module cable_wblake_mod
  implicit none

  real  :: wblake_ratio               ! ratio of wblake/subroff

  real, allocatable, save ::                                                   &
    WBLAKE_cable(:,:), TOT_WBLAKE_cable(:), TOT_SUBRUN_cable(:)

contains

subroutine cable_wblake_fix_alloc( land_points, ntiles )
  integer :: land_points, ntiles

  if(.NOT. allocated( wblake_cable ) ) & 
    allocate ( WBLAKE_cable(land_points,ntiles) )

End subroutine cable_wblake_fix_alloc

End module cable_wblake_mod
