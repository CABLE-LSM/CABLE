MODULE cable_soil_type_mod

IMPLICIT NONE

PUBLIC :: soil_parameter_type

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

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

   REAL, DIMENSION(:,:), POINTER ::                                            &
    heat_cap_lower_limit, &
    zse_vec,              &
    css_vec,              &
    cnsd_vec

   REAL, DIMENSION(:), POINTER ::                                              &
     cnsd,    & ! thermal conductivity of dry soil [W/m/K]
     pwb_min    ! working variable (swilt/ssat)**ibp2

   REAL, DIMENSION(:,:), POINTER ::                                            &
     albsoil    ! soil reflectance (2nd dim. BP 21Oct2009)
   
   ! parameters for GW module that vary with soil layer
   REAL, DIMENSION(:,:), POINTER ::                                            &
     sucs_vec,   & !psi at saturation in [mm]
     hyds_vec,   & !saturated hydraulic conductivity  [mm/s]
     bch_vec,    & !C and H B [none]
     clay_vec,   & !fraction of soil that is clay [frac]
     sand_vec,   & !fraction of soil that is sand [frac]
     silt_vec,   & !fraction of soil that is silt [frac]
     org_vec,    & !fration of soil made of organic soils [frac]
     rhosoil_vec,& !soil density  [kg/m3]
     ssat_vec,   & !volumetric water content at saturation [mm3/mm3]
     watr,       & !residual water content of the soil [mm3/mm3]
     sfc_vec,    & !field capcacity (hk = 1 mm/day)
     swilt_vec     ! wilting point (hk = 0.02 mm/day)

   REAL, DIMENSION(:), POINTER ::                                              &
     drain_dens, & !  drainage density ( mean dist to rivers/streams )
     elev,       & ! elevation above sea level
     elev_std,   & ! elevation above sea level
     slope,      & ! mean slope of grid cell
     slope_std     ! stddev of grid cell slope

   !MD parameters for GW module for the aquifer
   REAL, DIMENSION(:), POINTER ::                                       &
     GWsucs_vec, & ! head in the aquifer [mm]
     GWhyds_vec, & ! saturated hydraulic conductivity of aquifer [mm/s]
     GWbch_vec,  & ! clapp and horn b of the aquifer   [none]
     GWssat_vec, & ! saturated water content of the aquifer [mm3/mm3]
     GWwatr,     & ! residual water content of the aquifer [mm3/mm3]
     GWz,        & ! node depth of the aquifer    [m]
     GWdz,       & ! thickness of the aquifer   [m]
     GWrhosoil_vec !density of the aquifer substrate [kg/m3]

   ! Additional SLI parameters
   INTEGER, POINTER :: nhorizons(:)   ! number of soil horizons
   INTEGER, POINTER :: ishorizon(:,:) ! horizon number 1:nhorizons
   REAL,    POINTER :: clitt(:)       ! litter (tC/ha)
   REAL,    POINTER :: zeta(:)        ! macropore parameter
   REAL,    POINTER :: fsatmax(:)     ! variably saturated area parameter

END TYPE soil_parameter_type

END MODULE cable_soil_type_mod

