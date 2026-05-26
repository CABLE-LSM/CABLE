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

  SUBROUTINE set_JULES_surface_types(nnpft, npft, nnvg, ntype, urban, lake,&
      soil, ice)
    INTEGER, INTENT(OUT) :: nnpft, npft, nnvp, ntype, urban, lake,&
      soil, ice

    npft = 13
    nnvg = 4
    urban = urban_cable
    lake = lake_cable
    soil = barren_cable
    ice = ice_cable

    ntype = npft - nnvg
    nnpft = npft - ncpft

  END SUBROUTINE set_JULES_surface_types

END MODULE cable_surface_types_mod
