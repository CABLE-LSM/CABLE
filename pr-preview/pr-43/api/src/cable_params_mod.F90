MODULE cable_params_mod
!jhan: This is currently hard-wired as USEd module was NA in ESM1.5 and NOT ideal to be 
! USEing any data, especially JULES da
!USE max_dimensions,             ONLY: ntype_max ! defined PARAMETER @ compile time

!H! Elevate these to namelist definable
USE cable_other_constants_mod,  ONLY: nsl
USE cable_other_constants_mod,  ONLY: nrb       ! # radiation "bANDS" 
                                                !dir/dif components in bands VIS/NIR
USE cable_other_constants_mod,  ONLY: nscs      ! number of soil carbon stores
USE cable_other_constants_mod,  ONLY: nvcs      ! number of vegetation carbon stores
USE grid_constants_cbl_mod, ONLY : mstype => nsoiltypes  ! # of soil types [9]

IMPLICIT NONE

PUBLIC :: veg_parameter_type
PUBLIC :: vegin_type
PUBLIC :: soil_parameter_type
PUBLIC :: soilin_type

!jhan: eventuaLLY REMOVE THESE instances as PUBLIC
PUBLIC :: veg_cbl
PUBLIC :: vegin
PUBLIC :: soil_cbl
PUBLIC :: soilin

INTEGER, PARAMETER :: ntype_max = 17 
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
          meth,    & ! method for calculation of canopy fluxes and temp.
          iveg , &      ! vegetation type
          iLU ! land use type
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
          conkc0,  &  ! Michaelis-menton constant for carboxylase
          conko0,  &  ! Michaelis-menton constant for oxygenase
          ekc,     &  ! activation energy for caroxylagse
          eko,     &  ! acvtivation enegery for oxygenase
          g0,      & ! Belinda's stomatal model intercept, Ticket #56.
          g1         ! Belinda's stomatal model slope, Ticket #56.

     LOGICAL, DIMENSION(:), POINTER ::                                        &
          deciduous ! flag used for phenology fix

     REAL, DIMENSION(:,:), POINTER ::                                         &
          refl,    &
          taul,    &
          froot      ! fraction of root in each soil layer

     ! Additional  veg parameters:
     REAL, DIMENSION(:), POINTER :: rootbeta ! parameter for estimating vertical root mass distribution (froot)
     REAL, DIMENSION(:), POINTER :: gamma    ! parameter in root efficiency function (Lai and Katul 2000)
     REAL, DIMENSION(:), POINTER :: ZR       ! maximum rooting depth (cm)
     REAL, DIMENSION(:), POINTER :: F10      ! fraction of roots in top 10 cm

     REAL, DIMENSION(:), POINTER :: clitt     !

     ! Additional POP veg param
     INTEGER, DIMENSION(:,:), POINTER ::  disturbance_interval
     REAL, DIMENSION(:,:), POINTER ::  disturbance_intensity

  REAL,  POINTER ::                                                                     &
       soil(:,:),                                                  &
       csoil(:,:),                                                  &
       cplant(:,:),                                                 &
       ratecs(:,:),                                                            &
       ratecp(:,:)

  END TYPE veg_parameter_type


! Vegetation parameters I/O:
TYPE vegin_type
  REAL ::                                                                     &
       canst1(ntype_max),                                                      &
       dleaf(ntype_max),                                                       &
       length(ntype_max),                                                      &
       width(ntype_max),                                                       &
       vcmax(ntype_max),                                                       &
       ejmax(ntype_max),                                                       &
       hc(ntype_max),                                                          &
       xfang(ntype_max),                                                       &
       rp20(ntype_max),                                                        &
       rpcoef(ntype_max),                                                      &
       rs20(ntype_max),                                                        &
       wai(ntype_max),                                                         &
       rootbeta(ntype_max),                                                    &
       shelrb(ntype_max),                                                      &
       vegcf(ntype_max),                                                       &
       frac4(ntype_max),                                                       &
       xalbnir(ntype_max),                                                     &
       extkn(ntype_max),                                                       &
       tminvj(ntype_max),                                                      &
       tmaxvj(ntype_max),                                                      &
       vbeta(ntype_max),                                                       &
       a1gs(ntype_max),                                                        &
       d0gs(ntype_max),                                                        &
       alpha(ntype_max),                                                       &
       convex(ntype_max),                                                      &
       cfrd(ntype_max),                                                        &
       gswmin(ntype_max),                                                      &
       conkc0(ntype_max),                                                      &
       conko0(ntype_max),                                                      &
       ekc(ntype_max),                                                         &
       eko(ntype_max),                                                         &
       g0(ntype_max),                                                          &
       g1(ntype_max),                                                          &
       zr(ntype_max),                                                          &
       clitt(ntype_max),                                                       &
       froot(nsl,ntype_max),                                                   &
       csoil(nscs,ntype_max),                                                  &
       ratecs(nscs,ntype_max),                                                 &
       cplant(nvcs,ntype_max),                                                 &
       ratecp(nvcs,ntype_max),                                                 &
       refl(nrb,ntype_max),                                                    &
       taul(nrb,ntype_max)
END TYPE vegin_type

! Soil parameters used:
  TYPE soil_parameter_type

     INTEGER, DIMENSION(:), POINTER ::                                        &
          isoilm     ! integer soil type

     REAL, DIMENSION(:), POINTER ::                                           &
          bch,     & ! parameter b in Campbell equation
          c3,      & ! c3 drainage coeff (fraction)
          clay,    & ! fraction of soil which is clay
          css,     & ! soil specific heat capacity [kJ/kg/K]
          hsbh,    & ! difsat * etasat (=hyds*abs(sucs)*bch)
          hyds,    & ! hydraulic conductivity @ saturation [m/s], Ksat
          i2bp3,   & ! par. one in K vis suction (=nint(bch)+2)
          ibp2,    & ! par. two in K vis suction (fn of pbch)
          rhosoil, & ! soil density [kg/m3]
          sand,    & ! fraction of soil which is sand
          sfc,     & ! vol H2O @ field capacity
          silt,    & ! fraction of soil which is silt
          ssat,    & ! vol H2O @ saturation
          sucs,    & ! suction at saturation (m)
          swilt,   & ! vol H2O @ wilting
          zse,     & ! thickness of each soil layer (1=top) in m
          zshh,    & ! distance between consecutive layer midpoints (m)
                                ! vars intro for Ticket #27
          soilcol, & ! keep color for all patches/tiles
          albsoilf   ! soil reflectance

     REAL, DIMENSION(:,:), POINTER :: &
          zse_vec,css_vec,cnsd_vec

     REAL, DIMENSION(:), POINTER ::                                      &
          cnsd,    & ! thermal conductivity of dry soil [W/m/K]
          pwb_min    ! working variable (swilt/ssat)**ibp2

     REAL, DIMENSION(:,:), POINTER ::                                         &
          albsoil    ! soil reflectance (2nd dim. BP 21Oct2009)
     !mrd561
     !MD parameters for GW module that vary with soil layer
     REAL, DIMENSION(:,:), POINTER ::                                    &
          sucs_vec, & !psi at saturation in [mm]
          hyds_vec,  & !saturated hydraulic conductivity  [mm/s]
          bch_vec, & !C and H B [none]
          clay_vec,  & !fraction of soil that is clay [frac]
          sand_vec,  & !fraction of soil that is sand [frac]
          silt_vec,  & !fraction of soil that is silt [frac]
          org_vec,   & !fration of soil made of organic soils [frac]
          rhosoil_vec,& !soil density  [kg/m3]
          ssat_vec, & !volumetric water content at saturation [mm3/mm3]
          watr,   & !residual water content of the soil [mm3/mm3]
          sfc_vec, & !field capcacity (hk = 1 mm/day)
          swilt_vec     ! wilting point (hk = 0.02 mm/day)

     REAL, DIMENSION(:), POINTER ::                                      &
          drain_dens,&!  drainage density ( mean dist to rivers/streams )
          elev, &  !elevation above sea level
          elev_std, &  !elevation above sea level
          slope,  &  !mean slope of grid cell
          slope_std  !stddev of grid cell slope

     !MD parameters for GW module for the aquifer
     REAL, DIMENSION(:), POINTER ::                                       &
          GWsucs_vec,  &  !head in the aquifer [mm]
          GWhyds_vec,   &  !saturated hydraulic conductivity of the aquifer [mm/s]
          GWbch_vec,  & !clapp and horn b of the aquifer   [none]
          GWssat_vec,  & !saturated water content of the aquifer [mm3/mm3]
          GWwatr,    & !residual water content of the aquifer [mm3/mm3]
          GWz,       & !node depth of the aquifer    [m]
          GWdz,      & !thickness of the aquifer   [m]
          GWrhosoil_vec    !density of the aquifer substrate [kg/m3]

     ! Additional SLI parameters
     INTEGER,   DIMENSION(:),   POINTER :: nhorizons ! number of soil horizons
     INTEGER,   DIMENSION(:,:), POINTER :: ishorizon ! horizon number 1:nhorizons
     REAL, DIMENSION(:),   POINTER :: clitt     ! litter (tC/ha)
     REAL, DIMENSION(:),   POINTER :: zeta      ! macropore parameter
     REAL, DIMENSION(:),   POINTER :: fsatmax   ! variably saturated area parameter

  END TYPE soil_parameter_type



! Soil parameters I/O:
TYPE soilin_type
  REAL ::                                                                     &
       silt(mstype),                                                        &
       clay(mstype),                                                        &
       sand(mstype),                                                        &
       swilt(mstype),                                                       &
       sfc(mstype),                                                         &
       ssat(mstype),                                                        &
       bch(mstype),                                                         &
       hyds(mstype),                                                        &
       sucs(mstype),                                                        &
       rhosoil(mstype),                                                     &
       css(mstype)
END TYPE soilin_type

!Instantiate types
TYPE(vegin_type) :: vegin !read from namelist
TYPE(soilin_type) :: soilin !read from namelist

TYPE(soil_parameter_type) :: soil_cbl !used in CABLE
TYPE(veg_parameter_type) :: veg_cbl  !used in CABLE

END MODULE cable_params_mod
