!******************************************************************************
! This source code is part of the Community Atmosphere Biosphere Land Exchange
! (CABLE) model. This work is licensed under the CSIRO Open Source Software
! License Agreement (variation of the BSD / MIT License).You may not use this
! this file except in compliance with this License. A copy of the License is
! available at https://trac.nci.org.au/trac/cable/wiki/license.
!******************************************************************************
MODULE cable_phys_constants_mod

!-----------------------------------------------------------------------------
! Description:
!   CABLE physical constants
!
! This MODULE is USEd throughout CABLE.
!
! Module specific documentation:https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow:https://trac.nci.org.au/trac/cable/wiki/TBC
!******************************************************************************

IMPLICIT NONE

PUBLIC

REAL, PARAMETER :: tfrz   = 273.16        ! Temp (K) corresp. to 0 C
REAL, PARAMETER :: sboltz = 5.67e-8       ! Stefan-Boltz. const (W/m2/K4)
REAL, PARAMETER :: emsoil = 1.0           ! soil emissivity
REAL, PARAMETER :: emleaf = 1.0           ! leaf emissivity
REAL, PARAMETER :: capp   = 1004.64    ! air spec. heat (J/kg/K)
REAL, PARAMETER :: hl = 2.5014e6       ! latent heat of vaporization (J/s/m2)
!Below are constants used in CABLE model which are not as yet used in JAC-6.2
REAL, PARAMETER :: hlf = 0.334e6          ! latent heat of fusion
REAL, PARAMETER :: hls = 2.8350e6         ! latent heatOFsublimation (J/kg)
REAL, PARAMETER :: dheat  = 21.5e-6       ! molecular diffusivity for heat
REAL, PARAMETER :: grav   = 9.8086        ! gravity acceleration (m/s2)
REAL, PARAMETER :: rgas   = 8.3143        ! universal gas const  (J/mol/K)
REAL, PARAMETER :: rmair  = 0.02897       ! molecular wt: dry air (kg/mol)
REAL, PARAMETER :: rmh2o  = 0.018016      ! molecular wt: water (kg/mol)
REAL, PARAMETER :: cgsnow = 2090.0        ! specific heat capacity for snow
REAL, PARAMETER :: cs_rho_ice = 1.9341e6  !heat capacity * density ice
REAL, PARAMETER :: cs_rho_wat = 4.218e6   ! heat capacity * density  water
REAL, PARAMETER :: csice = 2.100e3        ! specific heat capacity for ice
REAL, PARAMETER :: cswat = 4.218e3        ! specific heat capacity water
REAL, PARAMETER :: density_liq = 1000.0   ! density of liquid water
REAL, PARAMETER :: density_ice = 921.0    ! denisty of ice

! Teten coefficients
REAL, PARAMETER :: tetena = 6.106         ! Magnus Tetans (Murray 1967)
REAL, PARAMETER :: tetenb = 17.27
REAL, PARAMETER :: tetenc = 237.3
! mrd561 the parameters for sat above ice
REAL, PARAMETER :: tetena_ice = 6.1078    ! ??? refs?
REAL, PARAMETER :: tetenb_ice = 21.875
REAL, PARAMETER :: tetenc_ice = 265.5

! Aerodynamic parameters, diffusivities, water density:
REAL, PARAMETER :: vonk   = 0.40          ! von Karman constant
REAL, PARAMETER :: a33    = 1.25          ! inertial sublayer sw/us
REAL, PARAMETER :: csw    = 0.50          ! canopy sw decay (Weil theory)
REAL, PARAMETER :: ctl    = 0.40          ! Wagga wheat (RDD 1992-Challenges)
REAL, PARAMETER :: apol   = 0.70          ! Polhausen coeff: one-sided plate
REAL, PARAMETER :: prandt = 0.71          ! Prandtl number: visc/diffh
REAL, PARAMETER :: schmid = 0.60          ! Schmidt number: visc/diffw
REAL, PARAMETER :: diffwc = 1.60          ! diffw/diffc = H2O/CO2 diffusivity
REAL, PARAMETER :: rhow   = 1000.0        ! liquid water density   [kg/m3]
REAL, PARAMETER :: crd = 0.3              ! element drag coefficient
REAL, PARAMETER :: csd = 0.003            ! substrate drag coefficient

!jhan:hardwire for now. note beta2 = crd/csd
REAL, PARAMETER ::  beta2 = 0.3/0.003     ! ratio cr/cs
REAL, PARAMETER ::  ccd   = 15.0          ! constant in d/h equation
REAL, PARAMETER ::  ccw_c = 2.0           ! ccw=(zw-d)/(h-d)
REAL, PARAMETER ::  usuhm = 0.3           ! (max of us/uh)

! Turbulence parameters:
REAL, PARAMETER :: zetmul = 0.4     ! if niter=2, final zeta=zetmul*zetar(2)
                                    ! niter=4 ATM see cable_define_types.F90
REAL, PARAMETER :: zeta0  = 0.0     ! initial value of za/L
REAL, PARAMETER :: zetneg = -15.0   ! negative limit on za/L when niter>=3
REAL, PARAMETER :: zetpos = 1.0     ! positive limit on za/L when niter>=3
REAL, PARAMETER :: zdlin  = 1.0     ! height frac of d below which TL linear
REAL, PARAMETER :: umin   = 0.1     ! guarantees convergence, was 0.01

!model parameter shared across subroutines -> cable_phys_constants
REAL, PARAMETER :: snow_depth_thresh = 1.0

END MODULE cable_phys_constants_mod
