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
! Purpose: defines parameters, variables and derived types, allocation and
!          deallocation of these derived types
!
! Contact: Bernard.Pak@csiro.au
!
! History: Brings together define_dimensions and define_types from v1.4b
!          rs20 now in veg% instead of soil%
!          fes split into fess and fesp (though fes still defined)
!
! Jan 2016: Now includes climate% for use in climate variables required for
! prognostic phenology and potential veg type
! ==============================================================================

MODULE cable_def_types_mod

! Scattered throughout CABLE (science code) there are USEs of these vars which 
! are/were declared/defined in multiple places. We are migrating to only having 
! them defined in one place. Bringing them into the header here satisfies legacy
USE grid_constants_mod_cbl, ONLY: nrb, nrs, mp, swb
USE grid_constants_mod_cbl, ONLY: mstype => nsoil_max ! # of soil types [9]
USE grid_constants_mod_cbl, ONLY: mvtype => ntype_max ! # of PFT  types [17]
USE grid_constants_mod_cbl, ONLY: ms     => nsl       ! # soil layers !sm_levels in JULES IO                        
USE grid_constants_mod_cbl, ONLY: msn    => nsnl      ! # snow layers                              
USE grid_constants_mod_cbl, ONLY: ncp    => nvCs      ! # vegetation carbon s
USE grid_constants_mod_cbl, ONLY: ncs    => nsCs      ! # soil carbon stores
USE grid_constants_mod_cbl, ONLY: ntype_max           ! Max # tiles ! compile time constant
USE grid_constants_mod_cbl, ONLY: mf                  ! # leaves (sunlit, shaded)
USE grid_constants_mod_cbl, ONLY: niter               ! # iterations for za/L
USE cable_other_constants_mod, ONLY: r_2              ! currently DOUBLE precision was 
                                                      ! SELECTED_REAL_KIND(12, 50)

! Scattered throughout CABLE (science code) there are USEs of these vars. We are  
! ascertaining (per "module") how we can better describe this data structure. As 
! a first step they are now declared/defined in independent, per module, files. 
! Bringing them into the header here satisfies legacy
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_air_type_mod,       ONLY: air_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_soil_type_mod,      ONLY: soil_parameter_type => soil_type
USE cable_veg_type_mod,       ONLY: veg_parameter_type  => veg_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_climate_type_mod,   ONLY: climate_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type

IMPLICIT NONE

PUBLIC

! CABLE special KINDs for representing INTEGER/REAL values with at least 
! 10-digit precision. NA in UM/JULES anyway as -i8 -r8 compile flags overrride

INTEGER, PARAMETER :: i_d  = KIND(9)    ! this is useless but needs to be def
   
INTEGER :: mland                        ! # land grid cells where is this used? 
INTEGER, PARAMETER :: n_ktherm = 3      ! where is this used? remove? local? 

END MODULE cable_def_types_mod
