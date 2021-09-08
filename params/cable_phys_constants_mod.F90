MODULE cable_phys_constants_mod

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   CABLE physical constants
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in 
!-----------------------------------------------------------------------------

  REAL, PARAMETER :: tfrz   = 273.16     ! Temp (K) corresp. to 0 C
  REAL, PARAMETER :: sboltz = 5.67e-8    ! Stefan-Boltz. constant (W/m2/K4)
  REAL, PARAMETER :: emsoil = 1.0        ! soil emissivity
  REAL, PARAMETER :: emleaf = 1.0        ! leaf emissivity
  REAL, PARAMETER :: capp   = 1004.64    ! air spec. heat (J/kg/K)
  REAL, PARAMETER :: hl = 2.5014e6       ! latent heat of vaporization (J/s/m2)

      real ::                                                                  &
      hlf = 0.334e6, & ! latent heat of fusion
      hls = 2.8350e6, & ! latent heat of sublimation (J/kg)
      !hl = 2.5104e6, & ! air spec. heat (J/kg/K)
      !hlf = 0.335e6, & ! latent heat of fusion
      dheat  = 21.5E-6, & ! molecular diffusivity for heat
      !grav   = 9.80, & ! gravity acceleration (m/s2)
      grav   = 9.8086, & ! gravity acceleration (m/s2)
      rgas   = 8.3143, & ! universal gas const  (J/mol/K)
      rmair  = 0.02897, & ! molecular wt: dry air (kg/mol)
      rmh2o  = 0.018016, & ! molecular wt: water        (kg/mol)
      cgsnow = 2090.0,&      ! specific heat capacity for snow
      cs_rho_ice = 1.9341e6,&    !heat capacity * density ice
      cs_rho_wat = 4.218e6,&    ! heat capacity * density  water
      csice = 2.100e3,&      ! specific heat capacity for ice
      cswat = 4.218e3,&      ! specific heat capacity for water
      density_liq = 1000.0,  &  !density of liquid water
      density_ice = 921.0,&     !denisty of ice

      ! Teten coefficients
      tetena = 6.106, & ! ??? refs?  mrd561 Magnus Tetans (Murray 1967)
      tetenb = 17.27, &
      tetenc = 237.3, &
      !mrd561 the parameters for sat above ice
      tetena_ice = 6.1078, & ! ??? refs?
      tetenb_ice = 21.875, &
      tetenc_ice = 265.5, &

      ! Aerodynamic parameters, diffusivities, water density:
      vonk   = 0.40, & ! von Karman constant
      a33    = 1.25, & ! inertial sublayer sw/us
      csw    = 0.50, & ! canopy sw decay (Weil theory)
      ctl    = 0.40, & ! Wagga wheat (RDD 1992, Challenges)
      apol   = 0.70, & ! Polhausen coeff: single-sided plate
      prandt = 0.71, & ! Prandtl number: visc/diffh
      schmid = 0.60, & ! Schmidt number: visc/diffw
      diffwc = 1.60, & ! diffw/diffc = H2O/CO2 diffusivity
      rhow   = 1000.0, & ! liquid water density   [kg/m3]
      crd = 0.3,    & ! element drag coefficient
      csd = 0.003,  & ! substrate drag coefficient

      !jhan:hardwire for now. note beta2 = crd/csd
      beta2 = 0.3/0.003, & ! ratio cr/cs
      ccd   = 15.0,  & ! constant in d/h equation
      ccw_c = 2.0,   & ! ccw=(zw-d)/(h-d)
      usuhm = 0.3,   & ! (max of us/uh)

      ! Turbulence parameters:
      zetmul = 0.4,  & ! if niter=2, final zeta=zetmul*zetar(2)
                       ! NB> niter currently=4 see cable_define_types.F90
      zeta0  = 0.0,  & ! initial value of za/L
      zetneg = -15.0, & ! negative limit on za/L when niter>=3
      zetpos = 1.0,  & ! positive limit on za/L when niter>=3
      zdlin  = 1.0,  & ! height frac of d below which TL linear
      !revised upwards from 0.01 to guarantee convergence, unnecessary now ?
      umin   = 0.1

END MODULE cable_phys_constants_mod
