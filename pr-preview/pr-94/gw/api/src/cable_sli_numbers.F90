MODULE sli_numbers

  USE cable_def_types_mod,  ONLY: r_2, i_d

  IMPLICIT NONE

  ! define some numbers
  REAL(r_2), PARAMETER :: zero      = 0.0_r_2
  REAL(r_2), PARAMETER :: half      = 0.5_r_2
  REAL(r_2), PARAMETER :: one       = 1.0_r_2
  REAL(r_2), PARAMETER :: two       = 2.0_r_2
  REAL(r_2), PARAMETER :: four      = 4.0_r_2
  REAL(r_2), PARAMETER :: thousand  = 1000._r_2
  REAL(r_2), PARAMETER :: e1        = 1.e-1_r_2
  REAL(r_2), PARAMETER :: e2        = 1.e-2_r_2
  REAL(r_2), PARAMETER :: e3        = 1.e-3_r_2
  REAL(r_2), PARAMETER :: e4        = 1.e-4_r_2
  REAL(r_2), PARAMETER :: e5        = 1.e-5_r_2
  REAL(r_2), PARAMETER :: e6        = 1.e-6_r_2
  REAL(r_2), PARAMETER :: e7        = 1.e-7_r_2

  ! define some constants
  REAL(r_2), PARAMETER :: pi        = 3.1415927_r_2
  REAL(r_2), PARAMETER :: Tzero     = 273.16_r_2          ! Celcius -> Kelvin
  REAL(r_2), PARAMETER :: gravity   = 9.8086_r_2          ! gravitation constant [m/s2]
  REAL(r_2), PARAMETER :: Mw        = 0.018016_r_2        ! weight of 1 mol of water [kg]
  REAL(r_2), PARAMETER :: rmair     = 0.02897_r_2         ! molecular wt: dry air (kg/mol)
  REAL(r_2), PARAMETER :: Mw18      = 0.018_r_2           ! weight of 1 mol of water [kg] (main isotopologue only)
  REAL(r_2), PARAMETER :: cpa       = 1004.64_r_2         ! specific heat capacity of dry air at 0-40 degC [J/kgK]
  REAL(r_2), PARAMETER :: esata     = 6.106_r_2*100.0_r_2 ! constants for saturated vapour pressure calculation
  REAL(r_2), PARAMETER :: esatb     = 17.27_r_2           ! %
  REAL(r_2), PARAMETER :: esatc     = 237.3_r_2           ! %

  REAL(r_2), PARAMETER :: rlambda   = 2.442e6_r_2  ! latent heat of condensation at 25 degC [J/kg]
  REAL(r_2), PARAMETER :: lambdaf   = 335000._r_2  ! latent heat of fusion (J kg-1)
  REAL(r_2), PARAMETER :: lambdas   = 2835000._r_2 ! latent heat of sublimation (J kg-1)
  REAL(r_2), PARAMETER :: Dva       = 2.17e-5_r_2  ! vapour diffusivity of water in air at 0 degC [m2/s]
  REAL(r_2), PARAMETER :: rhow      = 1000.0_r_2   ! denisty of water [kg/m3]
  REAL(r_2), PARAMETER :: rhoi      = 920._r_2     ! density of ice (kg m-3)

  REAL(r_2), PARAMETER :: rhoa      = 1.184_r_2    ! denisty of dry air at std (25 degC) [kg/m3]
  REAL(r_2), PARAMETER :: rhocp     = 1189.8_r_2   ! cpa*rhoa at std (25 degC) [J/m3K]

  REAL(r_2), PARAMETER :: esata_ice = 611.2_r_2   ! constants for saturated vapour pressure calculation over ice (WMO, 2008)
  REAL(r_2), PARAMETER :: esatb_ice = 22.46_r_2   ! %
  REAL(r_2), PARAMETER :: esatc_ice = 272.62_r_2  ! %
  REAL(r_2), PARAMETER :: csice     = 2.100e3_r_2 ! specific heat capacity for ice
  REAL(r_2), PARAMETER :: cswat     = 4.218e3_r_2 ! specific heat capacity for water
  REAL(r_2), PARAMETER :: rgas      = 8.3143_r_2  ! universal gas const  (J/mol/K)
  REAL(r_2), PARAMETER :: kw        = 0.58_r_2    ! dito

  ! numerical limits
  REAL(r_2), PARAMETER :: dSfac        = 5.25_r_2
  !REAL(r_2), PARAMETER :: dSfac        = 1.25_r_2
  REAL(r_2), PARAMETER :: dpmaxr       = 0.5_r_2
  REAL(r_2), PARAMETER :: h0min        = -5.e-3_r_2
  REAL(r_2), PARAMETER :: snmin        = 0.005_r_2 ! depth of snowpack (m) without dedicated snow layer(s)
  REAL(r_2), PARAMETER :: fsnowliq_max = 0.03_r_2  ! max fraction of snow water in liquid phase
  INTEGER(i_d), PARAMETER :: nsnow_max = 1     ! maximum number of dedicated snow layers (1 or 2)

  REAL(r_2), PARAMETER :: dh0max    = 0.0001_r_2
  REAL(r_2), PARAMETER :: SLmax     = 1.01_r_2
  REAL(r_2), PARAMETER :: SLmin     = 0.001_r_2
  REAL(r_2), PARAMETER :: Smax      = 1.05_r_2
  REAL(r_2), PARAMETER :: h0max     = 0.005_r_2
  REAL(r_2), PARAMETER :: qprecmax  = 1.0e10_r_2
  !REAL(r_2), PARAMETER :: dSmax     = 0.5_r_2
  !REAL(r_2), PARAMETER :: dSmaxr    = 0.5_r_2
  REAL(r_2), PARAMETER :: dSmax     = 0.1_r_2
  !REAL(r_2), PARAMETER :: dSmaxr    = 0.1_r_2
  REAL(r_2), PARAMETER :: dSmaxr    = 0.4_r_2
  REAL(r_2), PARAMETER :: dtmax     = 86400._r_2
  REAL(r_2), PARAMETER :: dtmin     = 0.01_r_2
  REAL(r_2), PARAMETER :: dsmmax    = 1.0_r_2
  REAL(r_2), PARAMETER :: dTsoilmax = 30.0_r_2
  REAL(r_2), PARAMETER :: dTLmax    = 30.0_r_2
  !REAL(r_2), PARAMETER :: dTsoilmax = 1.0_r_2
  !REAL(r_2), PARAMETER :: dTLmax    = 1.0_r_2
  REAL(r_2), PARAMETER :: tol_dthetaldT = 1.e-12_r_2
  INTEGER(i_d), PARAMETER :: nsteps_ice_max = 20
  INTEGER(i_d), PARAMETER :: nsteps_max = 200
  !MC-ToDo! Identify why smaller time-steps are needed for isotopes
  ! With isotopes  REAL(r_2), PARAMETER :: dSmax=0.1_r_2, dSmaxr=0.1_r_2, dtmax=0.05_r_2*3600.,
  !                                        dtmin=0.0_r_2, dsmmax=1.0_r_2
  ! With isotopes   REAL(r_2), PARAMETER ::  dTsoilmax=1.0_r_2, dTLmax=1.0_r_2
  REAL(r_2), PARAMETER :: gf        = 1.0_r_2    ! gravity factor for flow direction (usually 1 for straight down).
  REAL(r_2), PARAMETER :: hmin      = -1.0e6_r_2 ! minimum matric head h (used by MF).
  REAL(r_2), PARAMETER :: csol      = 0.0_r_2    ! solute concentration (mol kg-1 soil water)
  REAL(r_2), PARAMETER :: rhmin     = 0.05_r_2   ! minimum relative humidity in soil and litter

  ! boundary conditions
  REAL(r_2), PARAMETER :: hbot  = 0.0_r_2
  CHARACTER(LEN=20)    :: botbc = "free drainage"
  ! CHARACTER(LEN=20)    :: botbc = "zero flux"
  ! CHARACTER(LEN=20)    :: botbc = "aquifer"
  ! CHARACTER(LEN=20)    :: botbc = "constant head"
  ! CHARACTER(LEN=20)    :: botbc = "seepage"

  ! Special setups for sli stand-alone, such as 1-8: testcases of Haverd & Cuntz (2010);
  !  0: normal run
  ! 11: Mizoguchi (1990) / Hansson et al. (2004) lab experiment of freezing unsaturated soil; etc.
  ! 16: Loetschental
  INTEGER(i_d) :: experiment = 0

  ! Steeper freezing curve factor: 1=normal, >1=steeper (e.g. 1.5-2.0)
  REAL(r_2), PARAMETER :: freezefac = 1.0_r_2

  ! Topmodel approach
  ! 0: normal ponding and free drainage
  ! 1: topmodel surface runoff
  ! 2: topmodel deep drainage
  ! 3: topmodel surface runoff and deep drainage
  INTEGER(i_d), PARAMETER :: topmodel = 0
  REAL(r_2),    PARAMETER :: alpha    = 0.1_r_2 ! anistropy param for lateral flow (topmodel)
  REAL(r_2),    PARAMETER :: fsat_max = 2.0_r_2 ! exponent for vertical profile of Ksat (topmodel)

  ! Thermal conductivity of soil
  ! 0: Campbell (1985)
  ! 1: Van de Griend and O'Neill (1986)
  INTEGER(i_d) :: ithermalcond = 0

  ! define types
  TYPE vars_met
     REAL(r_2) :: Ta, rha, rbw, rbh, rrc, Rn, Da, cva, civa, phiva, Rnsw
  END TYPE vars_met

  TYPE vars
     INTEGER(i_d) :: isat
     REAL(r_2)    :: h, phi, phiS, K, KS, Dv, cvsat, rh, phiv, phivS, kH
     REAL(r_2)    :: kE, kth, csoil, eta_th, hS, rhS, sl, cv, cvsatT, cvS, kv
     INTEGER(i_d) :: iice
     REAL(r_2)    :: thetai, thetal, phiT, KT, lambdav, lambdaf, Sliq
     REAL(r_2)    :: he, phie, Ksat ! air-entry potential values (different to core sp params for frozen soil)
     REAL(r_2)    :: dthetaldT
     REAL(r_2)    :: Tfrz ! freezing point (deg C) depends on csol and S
     REAL(r_2)    :: csoileff
     REAL(r_2)    :: zsat ! height of saturated/unsaturated boundary (relative to bottom of layer)
     REAL(r_2)    :: macropore_factor
  END TYPE vars

  TYPE vars_snow
     REAL(r_2), DIMENSION(nsnow_max):: depth, hsnow, hliq, dens, tsn, kH, kE, kth, &
          Dv, cv, sl, melt, &
          Jsensible, Jlatent,  deltaJlatent, deltaJsensible, fsnowliq_max
     REAL(r_2) ::  wcol, Qadv_snow, Qadv_rain, totdepth, J, &
          Qadv_melt, Qadv_vap, Qcond_net, &
          Qadv_transfer, Qmelt, Qtransfer,FluxDivergence, deltaJ, &
          Qvap, MoistureFluxDivergence, Qprec, Qevap, deltawcol
     INTEGER(i_d) :: nsnow, nsnow_last
  END TYPE vars_snow

  TYPE vars_aquifer
     INTEGER(i_d) :: isat
     REAL(r_2)    :: zdelta, zsoil, zzero, K, Wa, discharge, f, Rsmax, Sy
  END TYPE vars_aquifer

  TYPE params
     REAL(r_2) :: the, thre, he, lam, Ke, eta, thr
     REAL(r_2) :: KSe, phie, phiSe, rho, thw, thfc, kd, css, clay, tortuosity
     INTEGER(i_d) :: ishorizon ! horizon number with different soil properties
     REAL(r_2) :: zeta
     REAL(r_2) :: fsatmax
     REAL(r_2) :: lambc        ! original lam for storage
     REAL(r_2) :: LambdaS      ! thermal inertia of saturation for van de Griend & O'Neill (1986) thermal conductivity
  END TYPE params

  TYPE rapointer
     REAL(r_2), DIMENSION(:), POINTER :: p
  END TYPE rapointer

  TYPE solve_type ! for energy and moisture balance in rh0_sol, etc.
     INTEGER(i_d) :: k
     REAL(r_2)    :: T1, Ta, cva, Rnet, hr1, hra, Dv, gv, gh, Dh, dz, phie, he, K1, eta,lambda, Ks, lambdav
  END TYPE solve_type

END MODULE sli_numbers
