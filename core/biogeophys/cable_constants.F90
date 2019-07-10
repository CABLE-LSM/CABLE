! cable_constants.f90
!
! Source file containing constants for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach, Matthias Cuntz
!
! bugs to bernard.pak@csiro.au
!
! This file contains modules:
!  math_constants
!  physical_constants
!  other_constants
!  photosynthetic_constants
!  spatial_heterogeneity

MODULE math_constants

  USE cable_def_types_mod, ONLY : i_d, r_2

  IMPLICIT NONE

  REAL, PARAMETER :: pi     = 3.141592653589793238462643383279502884197
  REAL, PARAMETER :: pi180  = pi / 180.0 ! radians / degree
  REAL, PARAMETER :: two_pi = 2.0 * pi
  REAL(r_2), PARAMETER :: pi_r_2 = 3.141592653589793238462643383279502884197

END MODULE math_constants

!=========================================================================

MODULE physical_constants

  USE cable_def_types_mod, ONLY : i_d

  IMPLICIT NONE

  !mrd561
  !constants from/for soilsnow and gw_hydro
  REAL,    PARAMETER :: hl = 2.5014e6  ! air spec. heat (J/kg/K)
  REAL,    PARAMETER :: hlf = 0.334e6  ! latent heat of fusion
  REAL,    PARAMETER :: hls = 2.8350e6  ! latent heat of sublimation (J/kg)

  REAL,    PARAMETER :: capp   = 1004.64  ! air spec. heat capacity (J/kg/K)
  REAL,    PARAMETER :: dheat  = 21.5E-6  ! molecular diffusivity for heat
  REAL,    PARAMETER :: grav   = 9.80     ! gravity acceleration (m/s2)
  REAL,    PARAMETER :: rgas   = 8.3143   ! universal gas const  (J/mol/K)
  REAL,    PARAMETER :: rmair  = 0.02897  ! molecular wt: dry air (kg/mol)
  REAL,    PARAMETER :: rmh2o  = 0.018016 ! molecular wt: water       (kg/mol)
  REAL,    PARAMETER :: sboltz = 5.67e-8  ! Stefan-Boltz. constant (W/m2/K4)
  REAL,    PARAMETER :: tfrz   = 273.16   ! Temp (K) corresp. to 0 C
  ! Teten coefficients
  REAL,    PARAMETER :: tetena = 6.106    ! ??? refs?
  REAL,    PARAMETER :: tetenb = 17.27
  REAL,    PARAMETER :: tetenc = 237.3
  !mrd561 the parameters for sat above ice
  REAL,    PARAMETER :: tetena_ice = 6.1078  ! ??? refs?
  REAL,    PARAMETER :: tetenb_ice = 21.875 
  REAL,    PARAMETER :: tetenc_ice = 265.5 
  ! Aerodynamic parameters, diffusivities, water density:
  REAL,    PARAMETER :: vonk   = 0.40     ! von Karman constant
  REAL,    PARAMETER :: a33    = 1.25     ! inertial sublayer sw/us
  REAL,    PARAMETER :: csw    = 0.50     ! canopy sw decay (Weil theory)
  REAL,    PARAMETER :: ctl    = 0.40     ! Wagga wheat (RDD 1992, Challenges)
  REAL,    PARAMETER :: apol   = 0.70     ! Polhausen coeff: single-sided plate
  REAL,    PARAMETER :: prandt = 0.71     ! Prandtl number: visc/diffh
  REAL,    PARAMETER :: schmid = 0.60     ! Schmidt number: visc/diffw
  REAL,    PARAMETER :: diffwc = 1.60     ! diffw/diffc = H2O/CO2 diffusivity
  REAL,    PARAMETER :: rhow   = 1000.0   ! liquid water density   [kg/m3]
  REAL,    PARAMETER :: emleaf = 1.0      ! leaf emissivity
  REAL,    PARAMETER :: emsoil = 1.0      ! soil emissivity
  REAL,    PARAMETER :: cr     = 0.3      ! element drag coefficient
  REAL,    PARAMETER :: cs     = 0.003    ! substrate drag coefficient
  REAL,    PARAMETER :: beta   = cr/cs    ! ratio cr/cs
  REAL,    PARAMETER :: ccd    = 15.0     ! constant in d/h equation
  REAL,    PARAMETER :: ccw    = 2.0      ! ccw=(zw-d)/(h-d)
  REAL,    PARAMETER :: usuhm  = 0.3      ! (max of us/uh)
  ! Turbulence  parameters:
  INTEGER(i_d), PARAMETER :: niter  = 10       ! number of iterations for za/L
  REAL,    PARAMETER :: zetmul = 0.4      ! if niter=2, final zeta=zetmul*zetar(2)
  REAL,    PARAMETER :: zeta0  = 0.0      ! initial value of za/L
  REAL,    PARAMETER :: zetneg = -10.0    ! negative limit on za/L when niter>=3
  REAL,    PARAMETER :: zetpos = 0.5      ! positive limit on za/L when niter>=3
  REAL,    PARAMETER :: zdlin  = 1.0      ! height frac of d below which TL linear
  REAL,    PARAMETER :: umin   = 0.1

END MODULE physical_constants

!=========================================================================

MODULE other_constants

  USE cable_def_types_mod, ONLY : i_d, nrb,r_2

  IMPLICIT NONE

  REAL,    PARAMETER, DIMENSION(nrb) :: gauss_w =(/0.308,0.514,0.178/)   ! Gaussian integ. weights
  ! values in refl and taul are slightly modified since Oct07 and Mar08 (YP)
  ! leaf reflectance
  REAL,    PARAMETER, DIMENSION(nrb) :: refl    = (/ 0.1, 0.425, 0.02 /) ! mar08
  ! leaf transmittance
  REAL,    PARAMETER, DIMENSION(nrb) :: taul    = (/ 0.1, 0.425, 0.02 /) ! mar08
  INTEGER(i_d), PARAMETER                 :: istemp  = 4                      ! soil temp:     1,2,3,4 = FR,kf,mrr,mrrkf
  INTEGER(i_d), PARAMETER                 :: ismois  = 2                      ! soil moist:  1,2,3     = MP84,NP89,Richards
  INTEGER(i_d), PARAMETER                 :: isinf   = 2                      ! soil infilt: 1,2       = MP84, FC96
  INTEGER(i_d), PARAMETER                 :: isevap  = 2                      ! soil evap: 1,2,3 = alfa,beta,threshold
  INTEGER(i_d), PARAMETER                 :: itherm  = 1                      ! VW or KGK algorithm for hconds,rkapps
  INTEGER(i_d), PARAMETER                 :: irktem  = 5                      ! RK steps in soil temp schemes
  INTEGER(i_d), PARAMETER                 :: irkmoi  = 5                      ! RK steps in soil moisture schemes
  ! soil water  parameters:
  REAL,    PARAMETER                 :: etarct  = 0.7                    ! rel soil moisture for finding zst1,zst2
  REAL,    PARAMETER                 :: dbde    = 1.3333                 ! d(beta)/d(etar): 4/3, 8 for D78,KGK91
  INTEGER(i_d), PARAMETER                 :: istsw   = 1                      !
  INTEGER(i_d), PARAMETER                 :: iresp   = 0                      ! unscaled (iresp=0) or scaled (iresp=1) respiration

END MODULE other_constants


!=========================================================================

MODULE photosynthetic_constants

  USE cable_def_types_mod, ONLY : i_d

  IMPLICIT NONE

  INTEGER(i_d), PARAMETER :: maxiter         = 20       ! max # interations for leaf temperature
  ! a1c3 is defined inside cable_canopy.f90
  REAL,    PARAMETER :: a1c4_default    = 4.0
  REAL,    PARAMETER :: a1c3_default    = 9.0
  REAL,    PARAMETER :: d0c3_hawkesbury = 5000.0   ! Empirical coef for vpd sensitivity of stomata
  REAL,    PARAMETER :: d0c3_default    = 1500.0   ! Empirical coef for vpd sensitivity of stomata
  REAL,    PARAMETER :: d0c4_default    = 1500.0   ! Empirical coef for vpd sensitivity of stomata
  REAL,    PARAMETER :: alpha3          = 0.2
  REAL,    PARAMETER :: alpha4          = 0.05
  REAL,    PARAMETER :: cfrd3           = 0.015
  REAL,    PARAMETER :: cfrd4           = 0.025
  REAL,    PARAMETER :: conkc0          = 302.0E-6 !mol mol^-1
  REAL,    PARAMETER :: conko0          = 256.0E-3 !mol mol^-1
  REAL,    PARAMETER :: convx3          = 0.01
  REAL,    PARAMETER :: convx4          = 0.8
  REAL,    PARAMETER :: ekc             = 59430.0  !J mol^-1
  REAL,    PARAMETER :: eko             = 36000.0  !J mol^-1
  REAL,    PARAMETER :: gam0            = 28.0E-6  !mol mol^-1 @ 20C = 36.9 @ 25C
  REAL,    PARAMETER :: gam1            = 0.0509
  REAL,    PARAMETER :: gam2            = 0.0010
  REAL,    PARAMETER :: gsw03           = 0.01
  REAL,    PARAMETER :: gsw04           = 0.04
  REAL,    PARAMETER :: rgbwc           = 1.32
  REAL,    PARAMETER :: rgswc           = 1.57
  REAL,    PARAMETER :: tmaxj           = 45.0
  REAL,    PARAMETER :: tmaxv           = 45.0
  REAL,    PARAMETER :: tminj           = -5.0
  REAL,    PARAMETER :: tminv           = -5.0
  REAL,    PARAMETER :: toptj           = 20.0
  REAL,    PARAMETER :: toptv           = 20.0
  REAL,    PARAMETER :: trefk           = 298.2    !reference temperature K

END MODULE photosynthetic_constants
