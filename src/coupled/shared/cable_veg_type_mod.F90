MODULE cable_veg_type_mod

IMPLICIT NONE

PUBLIC :: veg_parameter_type

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Vegetation parameters used:
TYPE veg_parameter_type
  
  !jhan: surely %meth doesn't need to be 'mp' here
  INTEGER, DIMENSION(:), POINTER ::                                        &
    meth, & ! method for calculation of canopy fluxes and temp.
    iveg, & ! vegetation type
    iLU     ! land use type
  
  REAL, DIMENSION(:), POINTER ::                                           &
    canst1,  & ! max intercepted water by canopy (mm/LAI)
    dleaf,   & ! chararacteristc legnth of leaf (m)
    ejmax,   & ! max pot. electron transp rate top leaf(mol/m2/s)
    frac4,   & ! fraction of c4 plants
    hc,      & ! roughness height of canopy (veg - snow)
    vlai,    & ! leaf area index
    xalbnir, &
    rp20,    & ! plant respiration coefficient at 20 C
    rpcoef,  & ! temperature coef nonleaf plant respiration (1/C)
    rs20,    & ! soil respiration at 20 C [mol m-2 s-1]
    shelrb,  & ! sheltering factor (dimensionless)
    vegcf,   & ! kdcorbin, 08/10
    tminvj,  & ! min temperature of the start of photosynthesis
    toptvj,  & ! opt temperature of the start of photosynthesis
    tmaxvj,  & ! max temperature of the start of photosynthesis
    vbeta,   & !
    vcmax,   & ! max RuBP carboxylation rate top leaf (mol/m2/s)
    xfang,   & ! leaf angle PARAMETER
    extkn,   & ! extinction coef for vertical
    vlaimax, & ! extinction coef for vertical
    wai,     & ! wood area index (stem+branches+twigs)
    a1gs,    & ! a1 parameter in stomatal conductance model
    d0gs,    & ! d0 in stomatal conductance model
    alpha,   & ! initial slope of J-Q response curve
    convex,  & ! convexity of J-Q response curve
    cfrd,    & ! ratio of day respiration to vcmax
    gswmin,  & ! minimal stomatal conductance
    conkc0,  & ! Michaelis-menton constant for carboxylase
    conko0,  & ! Michaelis-menton constant for oxygenase
    ekc,     & ! activation energy for caroxylagse
    eko,     & ! acvtivation enegery for oxygenase
    g0,      & ! Belinda's stomatal model intercept, Ticket #56.
    g1         ! Belinda's stomatal model slope, Ticket #56.

  LOGICAL, DIMENSION(:), POINTER ::                                        &
    deciduous ! flag used for phenology fix

  REAL, DIMENSION(:,:), POINTER ::                                         &
    refl,    &
    taul,    &
    froot      ! fraction of root in each soil layer

  ! Additional  veg parameters:
  REAL, DIMENSION(:), POINTER :: rootbeta ! parameter for estimating vertical 
                                          ! root mass distribution (froot)
  REAL, DIMENSION(:), POINTER :: gamma    ! parameter: root efficiency function
                                          ! (Lai and Katul 2000)
  REAL, DIMENSION(:), POINTER :: ZR       ! maximum rooting depth (cm)
  REAL, DIMENSION(:), POINTER :: F10      ! fraction of roots in top 10 cm
  REAL, DIMENSION(:), POINTER :: clitt    !

  ! Additional POP veg param
  INTEGER, DIMENSION(:,:), POINTER ::  disturbance_interval
  REAL, DIMENSION(:,:),    POINTER ::  disturbance_intensity

  REAL,  POINTER ::                                                            &
    soil(:,:),                                                                 &
    csoil(:,:),                                                                &
    cplant(:,:),                                                               &
    ratecs(:,:),                                                               &
    ratecp(:,:)

END TYPE veg_parameter_type

END MODULE cable_veg_type_mod
