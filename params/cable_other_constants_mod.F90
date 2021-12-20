MODULE cable_other_constants_mod

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   Other CABLE constants
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

INTEGER, PARAMETER ::                                                          &
  swb = 2,           & ! 2 shortwave bands (initial division - visible /
                       !   near infrared)
  n_sw_bands = 4,    & ! total number of shortwave radiation bands
                       !   (above divided into direct / diffuse)
  nrb = 3,           & ! number of radiation bands normally in use
  mf = 2,            & ! types of leaves (sunlit / shaded)
  msn = 3,           & ! max # snow layers
  r_2  = SELECTED_REAL_KIND(12, 50), &!this will be removed
                       ! double precision real dimension
  niter = 4,         & ! number of iterations for za/L
  nsl = 6,           & ! number of soil layers
  nscs = 2,          & ! number of soil carbon stores
  nvcs = 3,          & ! number of vegetation carbon stores
  n_assim_rates = 3, & ! Rubisco, RuBP and Sink-limited rates of photosynthesis
  n_soiltypes = 9      ! number of soil types

REAL, PARAMETER ::                                                             &
  max_snow_depth = 50000.0,  & ! maximum depth of lying snow on tiles (kg/m2)
  init_snow_rho1l = 140.0      ! Initial value for snow mean density
! Gaussian integ. weights
!REAL, PARAMETER :: gauss_w(nrb)=(/0.308,0.514,0.178 /) ! F90 
REAL, PARAMETER :: gauss_w(nrb)=[0.308,0.514,0.178 ]    ! F03
REAL, PARAMETER :: rad_thresh = 0.001
                        ! minimum zenithal angle for downward SW radiation
REAL, PARAMETER :: lai_thresh = 0.001
                        ! threshold for minimum significant LAI

! minimum (cosine)zenith angle of sun signalling sunrise 
REAL, PARAMETER :: coszen_tols = 1.0e-4

REAL, PARAMETER :: z0surf_min = 1.e-7 ! min. roughness of bare soil surface
!H!REAL, PARAMETER :: z0snow_min = 1.e-7 ! min. roughness of bare snow surface

END MODULE cable_other_constants_mod
