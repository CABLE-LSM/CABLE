!#define UM_CBL YES
!******************************************************************************
! This source code is part of the Community Atmosphere Biosphere Land Exchange
! (CABLE) model. This work is licensed under the CSIRO Open Source Software
! License Agreement (variation of the BSD / MIT License).You may not use this
! this file except in compliance with this License. A copy of the License is
! available at https://trac.nci.org.au/trac/cable/wiki/license.
!******************************************************************************

MODULE cable_surface_types_mod 

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
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
INTEGER, PARAMETER :: c3_cropland          = 10 
INTEGER, PARAMETER :: wetland              = 11 
INTEGER, PARAMETER :: empty1               = 12
INTEGER, PARAMETER :: empty2               = 13
INTEGER, PARAMETER :: barren_cable         = 14
INTEGER, PARAMETER :: urban_cable          = 15
INTEGER, PARAMETER :: lakes_cable          = 16 
INTEGER, PARAMETER :: ice_cable            = 17 

END MODULE cable_surface_types_mod

