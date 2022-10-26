MODULE grid_constants_mod_cbl

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   CABLE constants that define the "grid" ie. dimensions of arrays/vectors
!   (replace constants USEd from elsewhere AS they are implemented)
!   Potentially these dims could be read/set from CABLE nmls. However,
!   NB. Simply changing these dims arbitrarily would render the model useless
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE Science
!-----------------------------------------------------------------------------

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
INTEGER, PARAMETER :: ntype_max = 17 ! Max # tiles
INTEGER, PARAMETER :: nsoil_max = 9  ! Max # soils

INTEGER, PARAMETER :: ntiles    = 17 ! # tiles !npft+nveg in JULES IO
INTEGER, PARAMETER :: nsoils    = 9  ! # soils
!-----------------------------------------------------------------------------

INTEGER, PARAMETER :: nsl       = 6  ! # soil layers !sm_levels in JULES IO
INTEGER, PARAMETER :: nsnl      = 3  ! # snow layers

INTEGER, PARAMETER :: nrb       = 3  ! # rad bands VISual/NIR + Legacy incl LW
INTEGER, PARAMETER :: nsCs      = 2  ! # soil carbon stores
INTEGER, PARAMETER :: nvCs      = 3  ! # vegetation carbon stores

END MODULE grid_constants_mod_cbl
