MODULE grid_constants_cbl_mod

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   CABLE constants
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE Science
!-----------------------------------------------------------------------------

! This file will need some reorganisation. Some of those values will need
! to be moved to namelists. We will also need to make sure CABLE science
! code will be using the values from this file and not create duplicates.

! number of radiation "bands" normally in use
INTEGER, PARAMETER ::                                                         &
  n_rad_components = 2,  &  ! number of components  DIRect(beam)/DIFfuse         
  n_rad_sw_bands   = 2,  &  ! number of spectral bands VISual/NIR
  !H!trb = n_rad_components * n_rad_sw_bands,  & ! total number of "bands" 
  trb = 3                   ! total number of "bands" 

INTEGER, PARAMETER :: tsl = 3     ! # snow layers
INTEGER, PARAMETER :: nsl = 6     ! # soil layers

INTEGER, PARAMETER :: nsoiltypes = 9    ! # total no of soil types

INTEGER, PARAMETER :: mf = 2          ! # leaves (sunlit, shaded)
INTEGER, PARAMETER :: niter = 4       ! number of iterations for za/L

INTEGER, PARAMETER :: ncp = 3         ! # vegetation carbon stores
INTEGER, PARAMETER :: ncs = 2         ! # soil carbon stores

!jh!elevate these to namelist status
INTEGER, PARAMETER :: mvtype = 17       ! total # vegetation types,   from input

!jh!elevate these to namelist status
INTEGER, PARAMETER ::                                                         &
       i_d  = KIND(9),                                                        &
#ifdef UM_BUILD 
       r_2  = KIND(1.0),&!SELECTED_REAL_KIND(12, 50), &
#else       
       r_2  = KIND(1.0d0),&!SELECTED_REAL_KIND(12, 50), &
#endif
       n_tiles = 17     ! ccc: not sure what that is # possible no of different

INTEGER, PARAMETER :: n_ktherm = 3  ! ccc not sure what that is

INTEGER, PARAMETER :: ICE_soiltype_cbl = 9 !indice of ice-soil_typw

INTEGER, PARAMETER :: n_assim_rates = 3 ! Rubisco,
                                        ! RuBP and 
                                        ! Sink-limited rates of photosynthesis

REAL, PARAMETER ::                                                            &
  max_snow_depth = 50000.0,  & ! maximum depth of lying snow on tiles (kg/m2)
  init_snow_rho1l = 140.0      ! Initial value for snow mean density

REAL, PARAMETER :: gauss_w(trb)=(/0.308,0.514,0.178 /) ! Gaussian integ. weights
REAL, PARAMETER :: rad_thresh = 0.001
                        ! minimum zenithal angle for downward SW radiation
REAL, PARAMETER :: lai_thresh = 0.001
                        ! threshold for minimum significant LAI

REAL, PARAMETER :: coszen_tols = 1.0e-4
                        !minimum (cosine)zenith angle of sun signalling sunrise 
!minimum roughness of bare soil surface
REAL, PARAMETER :: z0surf_min = 1.0e-7
!minimum roughness of bare snow surface
!H!real, parameter :: z0snow_min = 1.e-7

END MODULE grid_constants_cbl_mod
