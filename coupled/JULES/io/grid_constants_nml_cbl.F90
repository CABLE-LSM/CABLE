MODULE grid_constants_cbl_mod

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   Other CABLE constants
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in 
!-----------------------------------------------------------------------------

! number of radiation "bands" normally in use
INTEGER, PARAMETER ::                                                          &
  n_rad_components = 2,  &  ! number of components  DIRect(beam)/DIFfuse         
  n_rad_sw_bands   = 2,  &  ! number of spectral bands VISual/NIR
  !H!trb = n_rad_components * n_rad_sw_bands,  & ! total number of "bands" 
  trb = 3                   ! total number of "bands" 

INTEGER, PARAMETER :: tsl = 3     ! # snow layers

INTEGER, parameter :: nsoiltypes=9    ! # total no of soil types

INTEGER, PARAMETER :: mf = 2          ! # leaves (sunlit, shaded)
INTEGER, PARAMETER :: niter = 4       ! number of iterations for za/L

INTEGER, PARAMETER :: ncp = 3         ! # vegetation carbon stores
INTEGER, PARAMETER :: ncs = 2         ! # soil carbon stores


END MODULE grid_constants_cbl_mod
