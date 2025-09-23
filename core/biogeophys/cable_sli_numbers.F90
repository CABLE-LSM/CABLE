MODULE sli_numbers

  USE cable_def_types_mod,  ONLY: r2, i4

  IMPLICIT NONE

  ! define some numbers
  REAL(r2), PARAMETER :: zero      = 0.0_r2
  REAL(r2), PARAMETER :: half      = 0.5_r2
  REAL(r2), PARAMETER :: one       = 1.0_r2
  REAL(r2), PARAMETER :: two       = 2.0_r2
  REAL(r2), PARAMETER :: four      = 4.0_r2
  REAL(r2), PARAMETER :: thousand  = 1000._r2
  REAL(r2), PARAMETER :: e1        = 1.e-1_r2
  REAL(r2), PARAMETER :: e2        = 1.e-2_r2
  REAL(r2), PARAMETER :: e3        = 1.e-3_r2
  REAL(r2), PARAMETER :: e4        = 1.e-4_r2
  REAL(r2), PARAMETER :: e5        = 1.e-5_r2
  REAL(r2), PARAMETER :: e6        = 1.e-6_r2
  REAL(r2), PARAMETER :: e7        = 1.e-7_r2

  ! define some constants
  REAL(r2), PARAMETER :: pi        = 3.1415927_r2
  REAL(r2), PARAMETER :: Tzero     = 273.16_r2          ! Celcius -> Kelvin
  REAL(r2), PARAMETER :: gravity   = 9.8086_r2          ! gravitation constant [m/s2]
  REAL(r2), PARAMETER :: Mw        = 0.018016_r2        ! weight of 1 mol of water [kg]
  REAL(r2), PARAMETER :: rmair     = 0.02897_r2         ! molecular wt: dry air (kg/mol)
  REAL(r2), PARAMETER :: Mw18      = 0.018_r2           ! weight of 1 mol of water [kg] (main isotopologue only)
  REAL(r2), PARAMETER :: cpa       = 1004.64_r2         ! specific heat capacity of dry air at 0-40 degC [J/kgK]
  REAL(r2), PARAMETER :: esata     = 6.106_r2*100.0_r2 ! constants for saturated vapour pressure calculation
  REAL(r2), PARAMETER :: esatb     = 17.27_r2           ! %
  REAL(r2), PARAMETER :: esatc     = 237.3_r2           ! %

  REAL(r2), PARAMETER :: rlambda   = 2.442e6_r2  ! latent heat of condensation at 25 degC [J/kg]
  REAL(r2), PARAMETER :: lambdaf   = 335000._r2  ! latent heat of fusion (J kg-1)
  REAL(r2), PARAMETER :: lambdas   = 2835000._r2 ! latent heat of sublimation (J kg-1)
  REAL(r2), PARAMETER :: Dva       = 2.17e-5_r2  ! vapour diffusivity of water in air at 0 degC [m2/s]
  REAL(r2), PARAMETER :: rhow      = 1000.0_r2   ! denisty of water [kg/m3]
  REAL(r2), PARAMETER :: rhoi      = 920._r2     ! density of ice (kg m-3)

  REAL(r2), PARAMETER :: rhoa      = 1.184_r2    ! denisty of dry air at std (25 degC) [kg/m3]
  REAL(r2), PARAMETER :: rhocp     = 1189.8_r2   ! cpa*rhoa at std (25 degC) [J/m3K]

  REAL(r2), PARAMETER :: esata_ice = 611.2_r2   ! constants for saturated vapour pressure calculation over ice (WMO, 2008)
  REAL(r2), PARAMETER :: esatb_ice = 22.46_r2   ! %
  REAL(r2), PARAMETER :: esatc_ice = 272.62_r2  ! %
  REAL(r2), PARAMETER :: csice     = 2.100e3_r2 ! specific heat capacity for ice
  REAL(r2), PARAMETER :: cswat     = 4.218e3_r2 ! specific heat capacity for water
  REAL(r2), PARAMETER :: rgas      = 8.3143_r2  ! universal gas const  (J/mol/K)
  REAL(r2), PARAMETER :: kw        = 0.58_r2    ! dito

  ! numerical limits
  REAL(r2), PARAMETER :: dSfac        = 5.25_r2
  !REAL(r2), PARAMETER :: dSfac        = 1.25_r2
  REAL(r2), PARAMETER :: dpmaxr       = 0.5_r2
  REAL(r2), PARAMETER :: h0min        = -5.e-3_r2
  ! Changed from 0.005 to 0.0005 because of results at Sodankyla in 2000-2003 for ESMSnowMIP
  REAL(r2), PARAMETER :: snmin        = 0.0005_r2 ! depth of snowpack (m) without dedicated snow layer(s)
  REAL(r2), PARAMETER :: fsnowliq_max = 0.03_r2  ! max fraction of snow water in liquid phase
  INTEGER(i4), PARAMETER :: nsnow_max = 1     ! maximum number of dedicated snow layers (1 or 2)

  REAL(r2), PARAMETER :: dh0max    = 0.0001_r2
  REAL(r2), PARAMETER :: SLmax     = 1.01_r2
  REAL(r2), PARAMETER :: SLmin     = 0.001_r2
  REAL(r2), PARAMETER :: Smax      = 1.05_r2
  REAL(r2), PARAMETER :: h0max     = 0.005_r2
  REAL(r2), PARAMETER :: qprecmax  = 1.0e10_r2
  !REAL(r2), PARAMETER :: dSmax     = 0.5_r2
  !REAL(r2), PARAMETER :: dSmaxr    = 0.5_r2
  REAL(r2), PARAMETER :: dSmax     = 0.1_r2
  !REAL(r2), PARAMETER :: dSmaxr    = 0.1_r2
  REAL(r2), PARAMETER :: dSmaxr    = 0.4_r2
  REAL(r2), PARAMETER :: dtmax     = 86400._r2
  REAL(r2), PARAMETER :: dtmin     = 0.01_r2
  REAL(r2), PARAMETER :: dsmmax    = 1.0_r2
  REAL(r2), PARAMETER :: dTsoilmax = 30.0_r2
  REAL(r2), PARAMETER :: dTLmax    = 30.0_r2
  !REAL(r2), PARAMETER :: dTsoilmax = 1.0_r2
  !REAL(r2), PARAMETER :: dTLmax    = 1.0_r2
  REAL(r2), PARAMETER :: tol_dthetaldT = 1.e-12_r2
  INTEGER(i4), PARAMETER :: nsteps_ice_max = 20
  INTEGER(i4), PARAMETER :: nsteps_max = 200
  !MC - ToDo: Identify why smaller time-steps are needed for isotopes.
  ! With isotopes  REAL(r2), PARAMETER :: dSmax=0.1_r2, dSmaxr=0.1_r2, dtmax=0.05_r2*3600.,
  !                                        dtmin=0.0_r2, dsmmax=1.0_r2
  ! With isotopes   REAL(r2), PARAMETER ::  dTsoilmax=1.0_r2, dTLmax=1.0_r2
  REAL(r2), PARAMETER :: gf        = 1.0_r2    ! gravity factor for flow direction (usually 1 for straight down).
  REAL(r2), PARAMETER :: hmin      = -1.0e6_r2 ! minimum matric head h (used by MF).
  REAL(r2), PARAMETER :: csol      = 0.0_r2    ! solute concentration (mol kg-1 soil water)
  REAL(r2), PARAMETER :: rhmin     = 0.05_r2   ! minimum relative humidity in soil and litter

  ! boundary conditions
  REAL(r2), PARAMETER :: hbot  = 0.0_r2
  CHARACTER(LEN=20)    :: botbc = "free drainage"
  ! CHARACTER(LEN=20)    :: botbc = "zero flux"
  ! CHARACTER(LEN=20)    :: botbc = "aquifer"
  ! CHARACTER(LEN=20)    :: botbc = "constant head"
  ! CHARACTER(LEN=20)    :: botbc = "seepage"

  ! Special setups for sli stand-alone, such as 1-8: testcases of Haverd & Cuntz (2010);
  !  0: normal run
  ! 11: Mizoguchi (1990) / Hansson et al. (2004) lab experiment of freezing unsaturated soil; etc.
  ! 16: Loetschental
  INTEGER(i4) :: experiment = 0

  ! Steeper freezing curve factor: 1=normal, >1=steeper (e.g. 1.5-2.0)
  REAL(r2), PARAMETER :: freezefac = 1.0_r2

  ! Topmodel approach
  ! 0: normal ponding and free drainage
  ! 1: topmodel surface runoff
  ! 2: topmodel deep drainage
  ! 3: topmodel surface runoff and deep drainage
  INTEGER(i4), PARAMETER :: topmodel = 0
  REAL(r2),    PARAMETER :: alpha    = 0.1_r2 ! anistropy param for lateral flow (topmodel)
  REAL(r2),    PARAMETER :: fsat_max = 2.0_r2 ! exponent for vertical profile of Ksat (topmodel)

  ! Thermal conductivity of soil
  ! 0: Campbell (1985)
  ! 1: Van de Griend and O'Neill (1986)
  INTEGER(i4) :: ithermalcond = 0

  ! define types
  TYPE vars_met
     REAL(r2) :: Ta, rha, rbw, rbh, rrc, Rn, Da, cva, civa, phiva, Rnsw
  END TYPE vars_met

  TYPE vars
     INTEGER(i4) :: isat
     REAL(r2)    :: h, phi, phiS, K, KS, Dv, cvsat, rh, phiv, phivS, kH
     REAL(r2)    :: kE, kth, csoil, eta_th, hS, rhS, sl, cv, cvsatT, cvS, kv
     INTEGER(i4) :: iice
     REAL(r2)    :: thetai, thetal, phiT, KT, lambdav, lambdaf, Sliq
     REAL(r2)    :: he, phie, Ksat ! air-entry potential values (different to core sp params for frozen soil)
     REAL(r2)    :: dthetaldT
     REAL(r2)    :: Tfrz ! freezing point (deg C) depends on csol and S
     REAL(r2)    :: csoileff
     REAL(r2)    :: zsat ! height of saturated/unsaturated boundary (relative to bottom of layer)
     REAL(r2)    :: macropore_factor
  END TYPE vars

  TYPE vars_snow
     REAL(r2), DIMENSION(nsnow_max):: depth, hsnow, hliq, dens, tsn, kH, kE, kth, &
          Dv, cv, sl, melt, &
          Jsensible, Jlatent,  deltaJlatent, deltaJsensible, fsnowliq_max
     REAL(r2) ::  wcol, Qadv_snow, Qadv_rain, totdepth, J, &
          Qadv_melt, Qadv_vap, Qcond_net, &
          Qadv_transfer, Qmelt, Qtransfer, FluxDivergence, deltaJ, &
          Qvap, MoistureFluxDivergence, Qprec, Qevap, deltawcol
     INTEGER(i4) :: nsnow, nsnow_last
  END TYPE vars_snow

  TYPE vars_aquifer
     INTEGER(i4) :: isat
     REAL(r2)    :: zdelta, zsoil, zzero, K, Wa, discharge, f, Rsmax, Sy
  END TYPE vars_aquifer

  TYPE params
     REAL(r2) :: the, thre, he, lam, Ke, eta, thr
     REAL(r2) :: KSe, phie, phiSe, rho, thw, thfc, kd, css, clay, tortuosity
     INTEGER(i4) :: ishorizon ! horizon number with different soil properties
     REAL(r2) :: zeta
     REAL(r2) :: fsatmax
     REAL(r2) :: lambc        ! original lam for storage
     REAL(r2) :: LambdaS      ! thermal inertia of saturation for van de Griend & O'Neill (1986) thermal conductivity
  END TYPE params

  TYPE rapointer
     REAL(r2), DIMENSION(:), POINTER :: p => null()
  END TYPE rapointer

  TYPE solve_type ! for energy and moisture balance in rh0_sol, etc.
     INTEGER(i4) :: k
     REAL(r2)    :: T1, Ta, cva, Rnet, hr1, hra, Dv, gv, gh, Dh, dz, phie, he, K1, eta,lambda, Ks, lambdav
  END TYPE solve_type

END MODULE sli_numbers
