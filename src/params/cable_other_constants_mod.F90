!******************************************************************************
! This source code is part of the Community Atmosphere Biosphere Land Exchange
! (CABLE) model. This work is licensed under the CSIRO Open Source Software
! License Agreement (variation of the BSD / MIT License).You may not use this
! this file except in compliance with this License. A copy of the License is
! available at https://trac.nci.org.au/trac/cable/wiki/license.
!******************************************************************************
MODULE cable_other_constants_mod

!-----------------------------------------------------------------------------
! Description:
!   Other CABLE constants
!
! This MODULE is USEd throughout CABLE
!
! Module specific documentation:https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow:https://trac.nci.org.au/trac/cable/wiki/TBC
!-----------------------------------------------------------------------------
!CABLE science not yet in JAC uses msn to describe number of snow layers
USE grid_constants_mod_cbl, ONLY: nrb, nsl, nsCs, nvCs, msn => nsnl
USE grid_constants_mod_cbl, ONLY: n_soiltypes => nsoil_max ! # of soil types [9]
USE grid_constants_mod_cbl, ONLY: niter                    ! # iterations za/L
USE grid_constants_mod_cbl, ONLY: mf                       !sunlit/shaded leaves
USE grid_constants_mod_cbl, ONLY: swb  ! 2 shortwave bands (VIS,NIR)

IMPLICIT NONE

PUBLIC

REAL, PARAMETER :: gauss_w(nrb)=[0.308,0.514,0.178 ] ! Gaussian integ. weights

REAL, PARAMETER :: rad_thresh = 0.001 ! min. zenithal angle for downward SW
REAL, PARAMETER :: lai_thresh = 0.001 ! min. LAI to be considered as vegetated

INTEGER, PARAMETER :: r_2  = KIND(1.d0) ! SELECTED_REAL_KIND(12, 50)

REAL, PARAMETER ::                                                             &
  max_snow_depth = 50000.0,  & ! maximum depth of lying snow on tiles (kg/m2)
  init_snow_rho1l = 140.0      ! Initial value for snow mean density

! minimum (cosine)zenith angle of sun signalling sunrise
REAL, PARAMETER :: coszen_tols = 1.0e-4

REAL, PARAMETER :: z0surf_min = 1.0e-7 ! min. roughness of bare soil surface
!H!REAL, PARAMETER :: z0snow_min = 1.e-7 ! min. roughness of bare snow surface

REAL, PARAMETER :: wilt_limitfactor = 2.0 ! Used in lower limit of soil moisture

END MODULE cable_other_constants_mod
