!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Reads vegetation and soil parameter files, fills vegin, soilin
!          NB. Most soil parameters overwritten by spatially explicit datasets
!          input as ancillary file (for ACCESS) or surface data file (for offline)
!          Module enables accessibility of variables throughout CABLE
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: v2.0 vegin%dleaf now calculated from leaf length and width
!          Parameter files were read elsewhere in v1.8 (init_subrs)
!
! ==============================================================================

MODULE cable_common_module

USE cable_runtime_opts_mod ,ONLY : cable_user
USE cable_runtime_opts_mod ,ONLY : satuparam
USE cable_runtime_opts_mod ,ONLY : wiltparam

  IMPLICIT NONE

  !---allows reference to "gl"obal timestep in run (from atm_step)
  !---total number of timesteps, and processing node
  INTEGER, SAVE :: ktau_gl, kend_gl, knode_gl, kwidth_gl

  LOGICAL :: L_fudge = .FALSE.

  INTEGER, SAVE :: CurYear  ! current year of multiannual run

  ! set from environment variable $HOME
  CHARACTER(LEN=200) ::                                                       &
       myhome

  ! switch to calc sil albedo using soil colour - Ticket #27
  LOGICAL :: calcsoilalbedo = .FALSE.
  !---Lestevens Sept2012
  !---CASACNP switches and cycle index
  LOGICAL, SAVE :: l_casacnp,l_laiFeedbk,l_vcmaxFeedbk
   LOGICAL :: l_luc = .FALSE.
   LOGICAL :: l_thinforest = .FALSE.
   LOGICAL :: l_landuse = .FALSE.

  !---CABLE runtime switches def in this type
  TYPE kbl_internal_switches
     LOGICAL :: um = .FALSE., um_explicit = .FALSE., um_implicit = .FALSE.,   &
          um_radiation = .FALSE., um_hydrology = .FALSE., esm15 = .FALSE.
     LOGICAL :: offline = .FALSE., mk3l = .FALSE.
  END TYPE kbl_internal_switches

  ! instantiate internal switches
  TYPE(kbl_internal_switches), SAVE :: cable_runtime

  ! external files read/written by CABLE
  TYPE filenames_type

     CHARACTER(LEN=500) ::                                                        &
          met,        & ! name of file for CABLE input
          path='./',       & ! path for output and restart files for CABLE and CASA
          out,        & ! name of file for CABLE output
          log,        & ! name of file for execution log
          restart_in = ' ', & ! name of restart file to read
          restart_out,& ! name of restart file to read
          LAI,        & ! name of file for default LAI
          TYPE,       & ! file for default veg/soil type
          veg,        & ! file for vegetation parameters
          soil,       & ! name of file for soil parameters
          soilcolor,  & ! file for soil color(soilcolor_global_1x1.nc)
          inits,      & ! name of file for initialisations
          soilIGBP,   & ! name of file for IGBP soil map
          gw_elev,    & !name of file for gw/elevation data
          fxpft,      & !filename for PFT fraction and transition,wood harvest, secondary harvest
          fxluh2cable,& !filename for mapping 12 luc states into 17 CABLE PFT
          gridnew       !filename for updated gridinfo file                       


  END TYPE filenames_type

   TYPE(filenames_type) :: filename

  ! hydraulic_redistribution switch _soilsnow module
  LOGICAL :: redistrb = .FALSE.  

  TYPE organic_soil_params
     !Below are the soil properties for fully organic soil

     REAL ::    &
          hyds_vec_organic  = 1.0e-4,&
          sucs_vec_organic = 10.3,   &
          clappb_organic = 2.91,     &
          ssat_vec_organic = 0.9,    &
          watr_organic   = 0.1,     &
          sfc_vec_hk      = 1.157407e-06, &
          swilt_vec_hk      = 2.31481481e-8

  END TYPE organic_soil_params

  TYPE gw_parameters_type

     REAL ::                   &
          MaxHorzDrainRate=2e-4,  & !anisintropy * q_max [qsub]
          EfoldHorzDrainRate=2.0, & !e fold rate of q_horz
          MaxSatFraction=2500.0,     & !parameter controll max sat fraction
          hkrz=0.5,               & !hyds_vec variation with z
          zdepth=1.5,             & !level where hyds_vec(z) = hyds_vec(no z)
          frozen_frac=0.05,       & !ice fraction to determine first non-frozen layer for qsub
          SoilEvapAlpha = 1.0,    & !modify field capacity dependence of soil evap limit
          IceAlpha=3.0,           &
          IceBeta=1.0

     REAL :: ice_impedence=5.0

     TYPE(organic_soil_params) :: org

     INTEGER :: level_for_satfrac = 6
     LOGICAL :: ssgw_ice_switch = .FALSE.

     LOGICAL :: subsurface_sat_drainage = .TRUE.

  END TYPE gw_parameters_type

  TYPE(gw_parameters_type), SAVE :: gw_params

  REAL, SAVE ::        &!should be able to change parameters!
       max_glacier_snowd=1100.0,&
       snow_ccnsw = 2.0, &
                                !jh!an:clobber - effectively force single layer snow
                                !snmin = 100.0,      & ! for 1-layer;
       snmin = 1.,          & ! for 3-layer;
       max_ssdn = 750.0,    & !
       max_sconds = 2.51,   & !
       frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95)

contains 

  ELEMENTAL FUNCTION IS_LEAPYEAR( YYYY )
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: YYYY
    LOGICAL :: IS_LEAPYEAR

    IS_LEAPYEAR = .FALSE.
    IF ( ( ( MOD( YYYY,  4 ) .EQ. 0 .AND. MOD( YYYY, 100 ) .NE. 0 ) .OR. &
         MOD( YYYY,400 ) .EQ. 0 ) ) IS_LEAPYEAR = .TRUE.

  END FUNCTION IS_LEAPYEAR

  FUNCTION LEAP_DAY( YYYY )
    IMPLICIT NONE
    INTEGER :: YYYY, LEAP_DAY

    IF ( IS_LEAPYEAR ( YYYY ) ) THEN
       LEAP_DAY = 1
    ELSE
       LEAP_DAY = 0
    END IF
  END FUNCTION LEAP_DAY



END MODULE cable_common_module
