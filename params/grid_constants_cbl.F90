!******************************************************************************
! This source code is part of the Community Atmosphere Biosphere Land Exchange
! (CABLE) model. This work is licensed under the CSIRO Open Source Software
! License Agreement (variation of the BSD / MIT License).You may not use this
! this file except in compliance with this License. A copy of the License is
! available at https://trac.nci.org.au/trac/cable/wiki/license.
!******************************************************************************

MODULE grid_constants_mod_cbl

!-----------------------------------------------------------------------------
! Description:
!   CABLE constants that define the "grid" ie. dimensions of arrays/vectors
!   (replace constants USEd from elsewhere AS they are implemented)
!   Potentially these dims could be read/set from CABLE nmls. However,
!   NB. Simply changing these dims arbitrarily would render the model useless
!
! This MODULE is USEd throughout CABLE.
!
! Module specific documentation:https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow:https://trac.nci.org.au/trac/cable/wiki/TBC
!******************************************************************************

IMPLICIT NONE

PUBLIC

! Although it seems these model "dimensions" may indeed be be a configuration
! variable that can/should be set thru a namelist, CABLE has been developed
! (even at an algorithm level) with the implicit assumption of the values
! defined here and is unlikely to yield meaningful results (if it runs at all)
! should these be arbitrarily changed at runtime

! # of land-cover/soil types. Types of land-cover for veg+nonveg i.e.(13+4)
!-----------------------------------------------------------------------------
! Req'd to be defined at compile time to read in pars. strictly speaking these
! only need to be greater than ntiles, nsoils (below). However, there is no
! point in allocating useless space here
INTEGER, PARAMETER :: ntype_max = 17 ! Max # tiles ! compile time constant
INTEGER, PARAMETER :: nsoil_max = 9  ! Max # soils ! req'd to read in pars

INTEGER, PARAMETER :: nsl       = 6  ! # soil layers !sm_levels in JULES IO
INTEGER, PARAMETER :: nsnl      = 3  ! # snow layers
INTEGER, PARAMETER :: nrb       = 3  ! # rad bands VISual/NIR + Legacy incl LW
INTEGER, PARAMETER :: nrs       = 4  ! # streams (VIS+NIR)*(Direct+Diffuse)=4
INTEGER, PARAMETER :: nsCs      = 2  ! # soil carbon stores
INTEGER, PARAMETER :: nvCs      = 3  ! # vegetation carbon stores
INTEGER, PARAMETER :: ICE_SoilType = 9 ! SoilType Index (soilparm_cable.nml JAC)
INTEGER, PARAMETER :: lakes_cable  = 16! SoilType Index (soilparm_cable.nml JAC)

INTEGER, PARAMETER :: mf = 2          ! # leaves (sunlit, shaded)
INTEGER, PARAMETER :: niter = 4       ! number of iterations for za/L

! Strictly NOT a constant. # of active tiles, length of CABLE working vectors
INTEGER :: mp

END MODULE grid_constants_mod_cbl
