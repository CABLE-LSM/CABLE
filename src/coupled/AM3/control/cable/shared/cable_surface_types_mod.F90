MODULE cable_surface_types_mod

  IMPLICIT NONE

  PUBLIC

  ! cable_surface_type (nml) Index
  INTEGER, PARAMETER :: evergreen_needleleaf = 1
  INTEGER, PARAMETER :: evergreen_broadleaf  = 2
  INTEGER, PARAMETER :: deciduous_needleleaf = 3
  INTEGER, PARAMETER :: deciduous_broadleaf  = 4
  INTEGER, PARAMETER :: shrub_cable          = 5
  INTEGER, PARAMETER :: c3_grassland         = 6
  INTEGER, PARAMETER :: c4_grassland         = 7
  INTEGER, PARAMETER :: tundra               = 8
  INTEGER, PARAMETER :: c3_cropland          = 9 
  INTEGER, PARAMETER :: c4_cropland          = 10 
  INTEGER, PARAMETER :: wetland              = 11 
  INTEGER, PARAMETER :: aust_mesic           = 12
  INTEGER, PARAMETER :: aust_xeric           = 13
  INTEGER, PARAMETER :: barren_cable         = 14
  INTEGER, PARAMETER :: urban_cable          = 15
  INTEGER, PARAMETER :: lakes_cable          = 16 
  INTEGER, PARAMETER :: ice_cable            = 17 

CONTAINS

  SUBROUTINE link_JULES_surface_ids_to_CABLE(surface_type_ids, nnpft, ncpft,&
      npft, nnvg, ntype, urban, lake, soil, ice)
    ! JULES permits surface IDs to be assigned in the namelist, and needs to
    ! link these IDs to I/O machinery (via the surface_type_ids) array and to
    ! the science via the integer parameters. The parameters specify how many
    ! surfaces are PFTs (they always occupy the first N slots), how many
    ! non-PFT surfaces there are to round out the total count, which IDs are
    ! urban, lake, soil or ice as they get special treatment.
    INTEGER, INTENT(INOUT) :: surface_type_ids
    INTEGER, INTENT(OUT) :: nnpft, ncpft, npft, nnvp, ntype, urban, lake,&
      soil, ice

    surface_type_ids(evergreen_needleleaf) = 1
    surface_type_ids(evergreen_broadleaf ) = 2
    surface_type_ids(deciduous_needleleaf) = 3
    surface_type_ids(deciduous_broadleaf ) = 4
    surface_type_ids(shrub_cable         ) = 5
    surface_type_ids(c3_grassland        ) = 6
    surface_type_ids(c4_grassland        ) = 7
    surface_type_ids(tundra              ) = 8
    surface_type_ids(c3_cropland         ) = 9
    surface_type_ids(c4_cropland         ) = 10
    surface_type_ids(wetland             ) = 11
    surface_type_ids(aust_mesic          ) = 12
    surface_type_ids(aust_xeric          ) = 13
    surface_type_ids(barren_cable        ) = 14
    surface_type_ids(urban_cable         ) = 15
    surface_type_ids(lakes_cable         ) = 16
    surface_type_ids(ice_cable           ) = 17

    npft = 13
    nnvg = 4
    urban = urban_cable
    lake = lake_cable
    soil = barren_cable
    ice = ice_cable

    ntype = npft - nnvg
    nnpft = npft - ncpft

  END SUBROUTINE link_JULES_surface_ids_to_CABLE

END MODULE cable_surface_types_mod
