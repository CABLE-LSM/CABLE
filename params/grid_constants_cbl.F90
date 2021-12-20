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

INTEGER, PARAMETER :: ntype_max = 17 ! Max # tiles ! compile time constant
INTEGER, PARAMETER :: nsoil_max = 9  ! Max # soils ! req'd to read in pars
INTEGER, PARAMETER :: nsl       = 6  ! # soil layers
INTEGER, PARAMETER :: nsnl      = 3  ! # snow layers
INTEGER, PARAMETER :: nrb       = 3  ! # rad bands VISual/NIR + Legacy incl LW
INTEGER, PARAMETER :: nsCs      = 2  ! # soil carbon stores
INTEGER, PARAMETER :: nvCs      = 3  ! # vegetation carbon stores

END MODULE grid_constants_mod_cbl
