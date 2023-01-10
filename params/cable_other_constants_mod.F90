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

USE grid_constants_mod_cbl, ONLY : nrb, nsl, nsCs, nvCs
USE grid_constants_mod_cbl, ONLY : msn =>  nsnl

IMPLICIT NONE

PUBLIC

REAL, PARAMETER :: gauss_w(nrb)=[0.308,0.514,0.178 ] ! Gaussian integ. weights

REAL, PARAMETER :: rad_thresh = 0.001 ! min. zenithal angle for downward SW 
REAL, PARAMETER :: lai_thresh = 0.001 ! min. LAI to be considered as vegetated

INTEGER, PARAMETER ::                                                          &
  swb = 2,           & ! 2 shortwave bands (initial division - visible /
                       ! near infrared)
  n_sw_bands = 4,    & ! total number of shortwave radiation bands
                       ! (above divided into direct / diffuse)
  mf = 2,            & ! types of leaves (sunlit / shaded)
  r_2  = SELECTED_REAL_KIND(12, 50), &!this will be removed
                       ! double precision real dimension
  niter = 4,         & ! number of iterations for za/L
  n_assim_rates = 3, & ! Rubisco, RuBP and Sink-limited rates of photosynthesis
  n_soiltypes = 9      ! number of soil types

REAL, PARAMETER ::                                                             &
  max_snow_depth = 50000.0,  & ! maximum depth of lying snow on tiles (kg/m2)
  init_snow_rho1l = 140.0      ! Initial value for snow mean density

! minimum (cosine)zenith angle of sun signalling sunrise
REAL, PARAMETER :: coszen_tols = 1.0e-4

REAL, PARAMETER :: z0surf_min = 1.0e-7 ! min. roughness of bare soil surface
!H!REAL, PARAMETER :: z0snow_min = 1.e-7 ! min. roughness of bare snow surface

END MODULE cable_other_constants_mod
