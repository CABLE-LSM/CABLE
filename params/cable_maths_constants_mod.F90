!******************************************************************************
! This source code is part of the Community Atmosphere Biosphere Land Exchange
! (CABLE) model. This work is licensed under the CSIRO Open Source Software
! License Agreement (variation of the BSD / MIT License).You may not use this
! this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) can be found
! at https://github.com/CABLE-LSM/CABLE/blob/main/
!******************************************************************************
MODULE cable_math_constants_mod

!------------------------------------------------------------------------------
! Description:
!   Mathematical constants used in CABLE.
!
! This MODULE is USEd in:
!     surf_couple_radiation_mod.F90 (JULES),
!     cable_cbm.F90 (ESM1.5),
!     cable_rad_driver.F90 (ESM1.5),
!     cbl_model_driver_offline.F90 (CABLE),
!     cable_canopy.F90 (CABLE),
!     cable_gw_hydro.F90 (CABLE),
!     cbl_sinbet.F90 (CABLE)
!
! Module specific documentation:https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow:https://trac.nci.org.au/trac/cable/wiki/TBC
!------------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC

REAL, PARAMETER :: pi    = 3.1415927
REAL, PARAMETER :: pi180 = pi / 180.0 ! radians / degree

END MODULE cable_math_constants_mod
