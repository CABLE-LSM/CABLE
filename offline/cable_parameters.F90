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
! Purpose:       This module file reads default parameter sets and basic
!                initialisations for CABLE. Parameters values are chosen based
!                on a global map of vegetation and soil types, currently based
!                on a 1x1-degree grid for offline case and host-model grid for
!                online case. Default initialisations are obtained from monthly
!                climatology in GSWP and Mk3L runs for offline and online
!                respectively.
!
! Contact: Bernard.Pak@csiro.au
!
! History: Changes since v1.4b for global offline (GSWP) cases, read in new
!          inpu files
!          Two subroutines moved to cable_common (reading veg and soil parameter
!          files)
!          Addition of code for CASA-CNP
!
!
! ==============================================================================
! CALLed from:   cable_input.F90
!
! MODULEs used:  cable_abort_module
!                cable_common_module
!                cable_def_types_mod
!                casadimension
!                casaparm
!                cable_IO_vars_module
!                phenvariable
!                physical_constants
!                netcdf

! CALLs:         get_default_params
!                read_gridinfo
!                spatialSoil
!                NSflip
!                countPatch
!                write_default_params
!                write_cnp_params
!                derived_parameters
!                check_parameter_values
!                report_parameters
!

MODULE cable_param_module

  USE cable_def_types_mod
  USE casadimension, ONLY: icycle
  USE casavariable
  USE phenvariable
  USE cable_abort_module
  USE cable_IO_vars_module
  USE cable_common_module, ONLY: cable_user, hide, &
                                 gw_params,init_veg_from_vegin,filename,&
                                 report_parameters_to_log,psi_c,psi_o
  USE cable_pft_params_mod
  USE cable_soil_params_mod
  USE CABLE_LUC_EXPT, ONLY: LUC_EXPT, LUC_EXPT_TYPE, LUC_EXPT_SET_TILES
  IMPLICIT NONE
  PRIVATE
  PUBLIC get_default_params, write_default_params, derived_parameters,         &
         check_parameter_values, report_parameters, parID_type,                &
         write_cnp_params,GWspatialParameters
  INTEGER :: patches_in_parfile=4 ! # patches in default global parameter
                                       ! file . which one?

  CHARACTER(LEN=4)  :: classification

  ! Variables below are temporary - for file read-in:
  INTEGER, DIMENSION(:, :, :),    ALLOCATABLE :: inVeg
  REAL,    DIMENSION(:, :, :),    ALLOCATABLE :: inPFrac
  INTEGER, DIMENSION(:, :),       ALLOCATABLE :: inSoil
  REAL,    DIMENSION(:, :, :, :), ALLOCATABLE :: inWB
  REAL,    DIMENSION(:, :, :, :), ALLOCATABLE :: inTGG
  REAL,    DIMENSION(:),          ALLOCATABLE :: inLon
  REAL,    DIMENSION(:),          ALLOCATABLE :: inLat
  REAL,    DIMENSION(:, :, :, :), ALLOCATABLE :: inALB
  REAL,    DIMENSION(:, :, :, :), ALLOCATABLE :: inSND
  REAL,    DIMENSION(:, :, :),    ALLOCATABLE :: inLAI
  REAL,    DIMENSION(:, :),       ALLOCATABLE :: inArea
  INTEGER, DIMENSION(:, :),       ALLOCATABLE :: inSorder
  REAL,    DIMENSION(:, :),       ALLOCATABLE :: inNdep
  REAL,    DIMENSION(:, :),       ALLOCATABLE :: inNfix
  REAL,    DIMENSION(:, :),       ALLOCATABLE :: inPwea
  REAL,    DIMENSION(:, :),       ALLOCATABLE :: inPdust

  ! Temporary values for reading IGBP soil map Q.Zhang @ 12/20/2010
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: inswilt
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: insfc
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: inssat
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: inbch
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: inhyds
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: insucs
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: inrhosoil
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: incss
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: incnsd
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: inclay
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: insilt
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: insand
  REAL,    DIMENSION(:, :),     ALLOCATABLE :: inORG

  ! vars intro for Ticket #27
  INTEGER, DIMENSION(:, :),     ALLOCATABLE :: inSoilColor


  INTERFACE get_gw_data
     MODULE PROCEDURE get_gw_2d_var_constdef
     MODULE PROCEDURE get_gw_3d_var_constdef
     MODULE PROCEDURE get_gw_4d_var_constdef
     MODULE PROCEDURE get_gw_2d_var
     MODULE PROCEDURE get_gw_3d_var
     MODULE PROCEDURE get_gw_4d_var
  END INTERFACE

CONTAINS

  SUBROUTINE get_default_params(logn, vegparmnew, LUC_EXPT,soil,ssnow)
    use cable_common_module, only : filename,             &
                                    calcsoilalbedo,cable_user
  ! Load parameters for each veg type and each soil type. (get_type_parameters)
  ! Also read in initial information for each grid point. (read_gridinfo)
  ! Count to obtain 'landpt', 'max_vegpatches' and 'mp'. (countPatch)
  !
  ! New input structure using netcdf and introduced 'month' to initialize
  ! soil profiles with the correct monthly average values (BP apr2010)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logn     ! log file unit number
    LOGICAL,      INTENT(IN) :: vegparmnew ! new format input file (BP dec2007)
    TYPE (LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (soil_snow_type),      INTENT(INOUT) :: ssnow

    ! local variables
    INTEGER :: npatch
    INTEGER :: nlon
    INTEGER :: nlat

    ! Get parameter values for all default veg and soil types:
    !CALL get_type_parameters(logn, vegparmnew, classification)
    CALL cable_pft_params()
    CALL cable_soil_params()
    WRITE(logn,*) ' Reading grid info from ', TRIM(filename%type)
    WRITE(logn,*) ' And assigning C4 fraction according to veg classification.'
    WRITE(logn,*)

    CALL read_gridinfo(nlon,nlat,npatch)

! Overwrite veg type and inital patch frac with land-use info
    IF (CABLE_USER%POPLUC) then
       CALL get_land_index(nlon, nlat)
       CALL LUC_EXPT_SET_TILES(inVeg, inPfrac, LUC_EXPT)
    ENDIF


    IF (soilparmnew) THEN
      PRINT *,      'Use spatially-specific soil properties; ', nlon, nlat
      WRITE(logn,*) 'Use spatially-specific soil properties; ', nlon, nlat
      CALL spatialSoil(nlon, nlat, logn)
    ENDIF

    ! include prescribed soil colour in determining albedo - Ticket #27
    IF (calcsoilalbedo) THEN
       CALL read_soilcolor(logn)
    END IF

    IF (cable_user%force_npatches_as .gt. 0) npatch=cable_user%force_npatches_as

    ! count to obtain 'landpt', 'max_vegpatches' and 'mp'
    CALL countPatch(nlon, nlat, npatch)

  END SUBROUTINE get_default_params
  !=============================================================================
  SUBROUTINE read_gridinfo(nlon, nlat, npatch)
  ! Reads in veg type, patch fraction, soil type, soil moisture and temperature
  ! profiles; also grid area and nutrients
  !
  ! Input variables:
  !   filename%type  - via cable_IO_vars_module
  !   classification - via cable_param_module
  ! Output variables:
  !   nlon           - # longitudes in input data set
  !   nlat           - # latitudes  in input data set
  !   npatch         - # patches in each grid from input data set
  !   inVeg          - via cable_param_module
  !   inPFrac        - via cable_param_module
  !   inSoil         - via cable_param_module
  !   inWB           - via cable_param_module
  !   inTGG          - via cable_param_module
  !   inLon          - via cable_param_module
  !   inLat          - via cable_param_module
  !   inALB          - via cable_param_module
  !   inSND          - via cable_param_module
  !   inLAI          - via cable_param_module
  !
  ! New input structure using netcdf and introduced 'month' to initialize
  ! soil profiles with the correct monthly average values (BP apr2010)

    USE netcdf
    use cable_common_module, only : filename

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: nlon
    INTEGER, INTENT(OUT) :: nlat
    INTEGER, INTENT(OUT) :: npatch

    ! local variables
    INTEGER :: ncid, ok
    INTEGER :: xID, yID, pID, sID, tID, bID
    INTEGER :: varID
    INTEGER :: nslayer, ntime, nband
    INTEGER :: ii, jj, kk,pp
    INTEGER, DIMENSION(:, :),     ALLOCATABLE :: idummy
    REAL,    DIMENSION(:, :),     ALLOCATABLE :: rdummy
    REAL,    DIMENSION(:, :, :),  ALLOCATABLE :: r3dum, r3dum2

    ok = NF90_OPEN(filename%type, 0, ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error opening grid info file.')

    ok = NF90_INQ_DIMID(ncid, 'longitude', xID)
    IF (ok /= NF90_NOERR) ok = NF90_INQ_DIMID(ncid, 'x', xID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring x dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, xID, LEN=nlon)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting x dimension.')
    ok = NF90_INQ_DIMID(ncid, 'latitude', yID)
    IF (ok /= NF90_NOERR) ok = NF90_INQ_DIMID(ncid, 'y', yID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring y dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, yID, LEN=nlat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting y dimension.')
    ok = NF90_INQ_DIMID(ncid, 'patch', pID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring patch dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, pID, LEN=npatch)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting patch dimension.')
    ok = NF90_INQ_DIMID(ncid, 'soil', sID)
    ok = NF90_INQUIRE_DIMENSION(ncid, sID, LEN=nslayer)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting soil dimension.')
    ok = NF90_INQ_DIMID(ncid, 'time', tID)
    ok = NF90_INQUIRE_DIMENSION(ncid, tID, LEN=ntime)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting time dimension.')
    ok = NF90_INQ_DIMID(ncid, 'rad', bID)
    ok = NF90_INQUIRE_DIMENSION(ncid, bID, LEN=nband)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting rad dimension.')

    !pass npatch to the def_types module
    !enables a future where npatch can be set to any value by reading in
    !tile info from gridinfo like file

    ! check dimensions of soil-layers and time
     !! vh_js !!
      IF ( (nslayer /= ms) .OR. (ntime /= 12)) THEN
         PRINT *, 'Variable dimensions do not match:'
         PRINT *, 'nslayer and ms = ', nslayer, ms
         PRINT *, 'ntime not equal 12 months: ', ntime
         IF (ntime /=12) THEN
            CALL abort('Variable dimensions do not match (read_gridinfo)')
         ELSE
            PRINT*, 'warning: soil layers below nslayer will be initialsed with moisture'
            PRINT*,    'and temperature of lowest layer in grid_info'
         ENDIF
      END IF


    ALLOCATE( inLon(nlon), inLat(nlat) )
    ALLOCATE( inVeg(nlon, nlat, npatch) )
    ALLOCATE( inPFrac(nlon, nlat, npatch) )
    ALLOCATE( inSoil(nlon, nlat) )
    ALLOCATE( idummy(nlon, nlat) )
    ALLOCATE( rdummy(nlon, nlat) )
    ALLOCATE(  inWB(nlon, nlat, nslayer,ntime) )
    ALLOCATE( inTGG(nlon, nlat, nslayer,ntime) )
    ALLOCATE( inALB(nlon, nlat, npatch,nband) )
    ALLOCATE( inSND(nlon, nlat, npatch,ntime) )
    ALLOCATE( inLAI(nlon, nlat, ntime) )
    ALLOCATE( r3dum(nlon, nlat, nband) )
    ALLOCATE( r3dum2(nlon, nlat, ntime) )

    ok = NF90_INQ_VARID(ncid, 'longitude', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error finding variable longitude.')
    ok = NF90_GET_VAR(ncid, varID, inLon)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error reading variable longitude.')

    !ensure this longitude is -180->180
    !as for GSWP3 it is 0-360
    WHERE (inLON > 180.0)
       inLON = inLON - 360.0
    ENDWHERE

    ok = NF90_INQ_VARID(ncid, 'latitude', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable latitude.')
    ok = NF90_GET_VAR(ncid, varID, inLat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable latitude.')

    ok = NF90_INQ_VARID(ncid, 'iveg', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable iveg.')
   ok = NF90_GET_VAR(ncid, varID, idummy)
    ok = NF90_GET_VAR(ncid, varID, inVeg)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable iveg.')
    inVeg(:, :, 1) = idummy(:,:) ! npatch=1 in 1x1 degree input

    ok = NF90_INQ_VARID(ncid, 'patchfrac', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error finding variable patchfrac.')
    ok = NF90_GET_VAR(ncid, varID, inPFrac)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error reading variable patchfrac.')
!CLN    inPFrac(:, :, 1) = rdummy(:, :)

    ok = NF90_INQ_VARID(ncid, 'isoil', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable isoil.')
    ok = NF90_GET_VAR(ncid, varID, inSoil)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable isoil.')

    ok = NF90_INQ_VARID(ncid, 'SoilMoist', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error finding variable SoilMoist.')
    ok = NF90_GET_VAR(ncid, varID, inWB)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error reading variable SoilMoist.')

    ok = NF90_INQ_VARID(ncid, 'SoilTemp', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable SoilTemp.')
    ok = NF90_GET_VAR(ncid, varID, inTGG)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable SoilTemp.')

    ok = NF90_INQ_VARID(ncid, 'Albedo', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable Albedo.')
    ok = NF90_GET_VAR(ncid, varID, r3dum)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable Albedo.')
 !   DO kk = 1, nband
 !     inALB(:,:,1,kk) = r3dum(:,:,kk)
 !   ENDDO
! vh!
    DO kk = 1, nband
       DO pp = 1,npatch
          inALB(:,:,pp,kk) = r3dum(:,:,kk)
       ENDDO
    ENDDO

    ok = NF90_INQ_VARID(ncid, 'SnowDepth', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error finding variable SnowDepth.')
    ok = NF90_GET_VAR(ncid,varID,r3dum2)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error reading variable SnowDepth.')
   ! DO kk = 1, ntime
   !   inSND(:, :, 1, kk) = r3dum2(:, :, kk)
   ! ENDDO

    DO kk = 1, ntime
       DO pp = 1,npatch
          inSND(:, :, pp, kk) = r3dum2(:, :, kk)
       ENDDO
    ENDDO

    ok = NF90_INQ_VARID(ncid, 'LAI', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable LAI.')
    ok = NF90_GET_VAR(ncid,varID,inLAI)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable LAI.')

    IF (icycle > 0) THEN
      ! casaCNP parameters
      ALLOCATE( inArea(nlon, nlat) )
      ALLOCATE( inSorder(nlon, nlat) )
      ALLOCATE( inNdep(nlon, nlat) )
      ALLOCATE( inNfix(nlon, nlat) )
      ALLOCATE( inPwea(nlon, nlat) )
      ALLOCATE( inPdust(nlon, nlat) )

      ok = NF90_INQ_VARID(ncid, 'area', varID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding area.')
      ok = NF90_GET_VAR(ncid, varID, inArea)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading area.')

      ok = NF90_INQ_VARID(ncid, 'SoilOrder', varID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding SoilOrder.')
      ok = NF90_GET_VAR(ncid, varID, inSorder)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading SoilOrder.')

      ok = NF90_INQ_VARID(ncid, 'Ndep', varID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding Ndep.')
      ok = NF90_GET_VAR(ncid, varID, inNdep)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading Ndep.')

      ok = NF90_INQ_VARID(ncid, 'Nfix', varID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding Nfix.')
      ok = NF90_GET_VAR(ncid, varID, inNfix)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading Nfix.')

      ok = NF90_INQ_VARID(ncid, 'Pwea', varID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding Pwea.')
      ok = NF90_GET_VAR(ncid, varID, inPwea)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading Pwea.')

      ok = NF90_INQ_VARID(ncid, 'Pdust', varID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding Pdust.')
      ok = NF90_GET_VAR(ncid, varID, inPdust)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading Pdust.')

      ! change units from g/m2/yr to g/m2/day
      inNdep  = inNdep  / 365.0
      inNfix  = inNfix  / 365.0
      inPwea  = inPwea  / 365.0
      inPdust = inPdust / 365.0

    ENDIF

    ok = NF90_CLOSE(ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error closing grid info file.')

  END SUBROUTINE read_gridinfo
  !============================================================================
  SUBROUTINE spatialSoil(nlon, nlat, logn)
  ! Read in spatially-specific soil properties including snow-free albedo
  ! plus soil texture; all these from UM ancilliary file
  !
  ! Input variables:
  !   nlon,nlat         - # longitudes and latitudes in the previous input file
  !   filename%soilIGBP - via cable_IO_vars_module
  ! Output variables:
  !   inswilt   - via cable_param_module
  !   insfc     - via cable_param_module
  !   inssat    - via cable_param_module
  !   inbch     - via cable_param_module
  !   inhyds    - via cable_param_module
  !   insucs    - via cable_param_module
  !   inrhosoil - via cable_param_module
  !   incss     - via cable_param_module
  !   incnsd    - via cable_param_module
  !   inclay    - via cable_param_module
  !   insilt    - via cable_param_module
  !   insand    - via cable_param_module
  !   inALB     - via cable_param_module

    USE netcdf
    use cable_common_module, only : filename

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlon
    INTEGER, INTENT(IN) :: nlat
    INTEGER, INTENT(IN) :: logn ! log file unit number

    ! local variables
    INTEGER :: ncid, ok, ii, jj, kk, ok2, ncid_elev
    INTEGER :: xID, yID, fieldID
    INTEGER :: xlon, xlat,ok_gw_elev
    REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: indummy
    REAL, DIMENSION(:,:),     ALLOCATABLE :: sfact, dummy2
    REAL, DIMENSION(:,:),     ALLOCATABLE :: in2alb

    ok = NF90_OPEN(filename%type, 0, ncid)

    ALLOCATE(    in2alb(nlon, nlat) ) ! local
    ALLOCATE(    dummy2(nlon, nlat) ) ! local
    ALLOCATE(     sfact(nlon, nlat) ) ! local
    ALLOCATE(   inswilt(nlon, nlat) )
    ALLOCATE(     insfc(nlon, nlat) )
    ALLOCATE(    inssat(nlon, nlat) )
    ALLOCATE(     inbch(nlon, nlat) )
    ALLOCATE(    inhyds(nlon, nlat) )
    ALLOCATE(    insucs(nlon, nlat) )
    ALLOCATE( inrhosoil(nlon, nlat) )
    ALLOCATE(     incss(nlon, nlat) )
    ALLOCATE(    incnsd(nlon, nlat) )
    ALLOCATE(    inclay(nlon, nlat) )
    ALLOCATE(    insilt(nlon, nlat) )
    ALLOCATE(    insand(nlon, nlat) )
    ALLOCATE(       inORG(nlon, nlat) )


    ! 1
    ok = NF90_INQ_VARID(ncid, 'swilt', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable swilt.')
    ok = NF90_GET_VAR(ncid, fieldID, inswilt)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable swilt.')
    ! 2
    ok = NF90_INQ_VARID(ncid, 'sfc', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable sfc.')
    ok = NF90_GET_VAR(ncid, fieldID, insfc)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable sfc.')
    ! 3
    ok = NF90_INQ_VARID(ncid, 'ssat', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable ssat.')
    ok = NF90_GET_VAR(ncid, fieldID, inssat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable ssat.')
    ! 4
    ok = NF90_INQ_VARID(ncid, 'bch', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable bch.')
    ok = NF90_GET_VAR(ncid, fieldID, inbch)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable bch.')
    ! 5
    ok = NF90_INQ_VARID(ncid,'hyds',fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable hyds.')
    ok = NF90_GET_VAR(ncid, fieldID, inhyds)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable hyds.')
    ! 6
    ok = NF90_INQ_VARID(ncid, 'sucs', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable sucs.')
    ok = NF90_GET_VAR(ncid, fieldID, insucs)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable sucs.')
    ! 7
    ok = NF90_INQ_VARID(ncid, 'rhosoil', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable rhosoil.')
    ok = NF90_GET_VAR(ncid,fieldID,inrhosoil)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable rhosoil.')
    ! 8
    ok = NF90_INQ_VARID(ncid, 'cnsd', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable cnsd.')
    ok = NF90_GET_VAR(ncid, fieldID, incnsd)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable cnsd.')
    ! 9
    ok = NF90_INQ_VARID(ncid, 'css', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable css.')
    ok = NF90_GET_VAR(ncid, fieldID, incss)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable css.')
    ! 10
    ok = NF90_INQ_VARID(ncid, 'clay', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable clay.')
    ok = NF90_GET_VAR(ncid, fieldID, inclay)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable clay.')
    ! 11
    ok = NF90_INQ_VARID(ncid, 'silt', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable silt.')
    ok = NF90_GET_VAR(ncid, fieldID, insilt)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable silt.')
    ! 12
    ok = NF90_INQ_VARID(ncid, 'sand', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable sand.')
    ok = NF90_GET_VAR(ncid, fieldID, insand)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable sand.')
    ! 13 UM albedo
    ok = NF90_INQ_VARID(ncid, 'albedo2', fieldID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable UM albedo')
    ok = NF90_GET_VAR(ncid, fieldID, in2alb)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable UM albedo')

    ok = NF90_INQ_VARID(ncid, 'organic', fieldID)
    ok2 = ok
    IF (ok .eq. NF90_NOERR) then
      ok2 = NF90_GET_VAR(ncid, fieldID, inORG)
      write(logn,*) 'READ FORG FROM THE DATA FILE, yeidling '
      write(logn,*) 'A maximum value of ',maxval(inORG),' and min val of',minval(inORG)
    end if
    IF ((ok2 .ne. NF90_NOERR) .or. (ok .ne. NF90_NOERR)) then
      inORG(:,:) = 0.0
      write(logn,*) 'COULD NOT READ FORG FROM THR SRF FILE setting to 0.0'
    END IF



! Use this code if need to process original UM file soil fields into CABLE
! offline format
!    ! 1
!    ok = NF90_INQ_VARID(ncid,'field329',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable swilt.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable swilt.')
!    inswilt(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,inswilt)
!    ! 2
!    ok = NF90_INQ_VARID(ncid,'field330',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable sfc.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable sfc.')
!    insfc(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,insfc)
!    ! 3
!    ok = NF90_INQ_VARID(ncid,'field332',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable ssat.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable ssat.')
!    inssat(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,inssat)
!    ! 4
!    ok = NF90_INQ_VARID(ncid,'field1381',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable bch.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable bch.')
!    inbch(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,inbch)
!    ! 5
!    ok = NF90_INQ_VARID(ncid,'field333',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable hyds.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable hyds.')
!    inhyds(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,inhyds)
!    ! 6
!    ok = NF90_INQ_VARID(ncid,'field342',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable sucs.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable sucs.')
!    insucs(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,insucs)
!    ! 7
!    ok = NF90_INQ_VARID(ncid,'field2011',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable rhosoil.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable rhosoil.')
!    inrhosoil(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,inrhosoil)
!    ! 8
!    ok = NF90_INQ_VARID(ncid,'field335',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable css.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable css.')
!    incss(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,incss)
!    ! 9
!    ok = NF90_INQ_VARID(ncid,'field336',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable cnsd.')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable cnsd.')
!    incnsd(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,incnsd)
!    ! 10 albedo
!    ok = NF90_INQ_VARID(ncid,'field1395',fieldID)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable albedo')
!    ok = NF90_GET_VAR(ncid,fieldID,indummy)
!    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable albedo')
!    in2alb(:,:) = indummy(:,:,1,1)
!    CALL NSflip(nlon,nlat,in2alb)

    ok = NF90_CLOSE(ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error closing IGBP soil map.')

! Code if using UM soil file
! unit change and glacial-point check were done in preprocessing
!    ! change unit to m/s
!    inhyds = inhyds * 1.0E-3
!    ! Assign values to glacial points which are zeroes
!    WHERE(inswilt==0.)     inswilt = 0.216
!    WHERE( insfc ==0.)       insfc = 0.301
!    WHERE(inssat ==0.)      inssat = 0.479
!    WHERE( inbch ==0.)       inbch = 7.1
!    WHERE(inhyds ==0.)      inhyds = 1.E-3
!    WHERE(insucs ==0.)      insucs = 0.153
!    WHERE(inrhosoil==0.) inrhosoil = 1455
!    WHERE(incnsd ==0.)      incnsd = 0.272
!    WHERE(incss > 630000.0)
!      incss = incss / inrhosoil   ! normal points need unit conversion
!    ELSEWHERE
!      incss = 2100.0   ! glacial points
!    ENDWHERE

    ! Calculate albedo for radiation bands and overwrite previous
    ! initialization
    PRINT *,      'When choosing spatially-specific soil properties,'
    PRINT *,      'snow-free albedo is also overwritten by this data set.'
    WRITE(logn, *) 'When choosing spatially-specific soil properties,'
    WRITE(logn, *) 'snow-free albedo is also overwritten by this data set.'
    sfact = 0.68
    WHERE (in2alb <= 0.14)
       sfact = 0.5
    ELSEWHERE (in2alb > 0.14)
       sfact = 0.62
    END WHERE
    WHERE (in2alb > 1.0e19)   ! ocean points
      in2alb = -1.0
    END WHERE
    dummy2(:, :) = 2.0 * in2alb(:, :) / (1.0 + sfact(:, :))
    inALB(:, :, 1, 2) = dummy2(:, :)
    inALB(:, :, 1, 1) = sfact(:, :) * dummy2(:, :)


    DEALLOCATE(in2alb, sfact, dummy2)
!    DEALLOCATE(in2alb,sfact,dummy2,indummy)

  END SUBROUTINE spatialSoil
  !=============================================================================
  !subr to read soil color for albed o calc - Ticket #27
  SUBROUTINE read_soilcolor(logn)
  ! Read soil color
  !
  ! Input variables:
  !   filename%soilcolor  - via cable_IO_vars_module
  ! Output variables:
  !   soilcol    - via cable_param_module
  !
  ! New input structure using netcdf

    USE netcdf
    USE cable_common_module, ONLY : filename, calcsoilalbedo
    ! USE cable_IO_vars_module, ONLY : soilcol

    IMPLICIT NONE
    ! INTEGER, DIMENSION(:), INTENT(INOUT) :: soilcol
    ! TYPE (soil_parameter_type), INTENT(OUT) :: soil
    INTEGER, INTENT(IN) ::  logn ! log file unit number

    ! local variables
    ! INTEGER, DIMENSION(:, :),     ALLOCATABLE :: inSoilColor
    INTEGER :: ncid, ok
    INTEGER :: nlon
    INTEGER :: nlat
    INTEGER :: xID, yID
    INTEGER :: varID
    INTEGER :: r, e

    REAL,    DIMENSION(:),          ALLOCATABLE :: inLonSoilCol
    REAL,    DIMENSION(:),          ALLOCATABLE :: inLatSoilCol

    ok = NF90_OPEN(filename%soilcolor, 0, ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error opening soil color file.')

    ok = NF90_INQ_DIMID(ncid, 'longitude', xID)
    IF (ok /= NF90_NOERR) ok = NF90_INQ_DIMID(ncid, 'x', xID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring x dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, xID, LEN=nlon)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting x dimension.')
    ok = NF90_INQ_DIMID(ncid, 'latitude', yID)
    IF (ok /= NF90_NOERR) ok = NF90_INQ_DIMID(ncid, 'y', yID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring y dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, yID, LEN=nlat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting y dimension.')


    ALLOCATE( inLonSoilCol(nlon), inLatSoilCol(nlat) )
    ALLOCATE( inSoilColor(nlon, nlat) )
    ! ALLOCATE( soilcol(mp) )

    ok = NF90_INQ_VARID(ncid, 'longitude', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error finding variable longitude.')
    ok = NF90_GET_VAR(ncid, varID, inLonSoilCol)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error reading variable longitude.')

    DO r = 1, nlon
      IF ( inLonSoilCol(r) /= inLon(r) ) CALL nc_abort(ok,                     &
                                               'Wrong resolution in longitude.')
    END DO

    ok = NF90_INQ_VARID(ncid, 'latitude', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable latitude.')
    ok = NF90_GET_VAR(ncid, varID, inLatSoilCol)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable latitude.')

    DO r = 1, nlat
      IF ( inLatSoilCol(r) /= inLat(r) ) CALL nc_abort(ok,                     &
                                               'Wrong resolution in latitude.')
    END DO

    ok = NF90_INQ_VARID(ncid, 'SOIL_COLOR', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable soil color.')
    ok = NF90_GET_VAR(ncid, varID, inSoilColor)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable soil color.')

    ok = NF90_CLOSE(ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error closing soil color file.')

  END SUBROUTINE read_soilcolor
  !=============================================================================
  SUBROUTINE NSflip(nlon, nlat, invar)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlon
    INTEGER, INTENT(IN) :: nlat
    REAL, INTENT(INOUT) :: invar(nlon,nlat)

    ! local variables
    INTEGER :: ii, jj
    REAL    :: rdummy(nlon, nlat)

    DO jj = 1, nlat
    DO ii = 1, nlon
      rdummy(ii, jj) = invar(ii, nlat - jj + 1)
    ENDDO
    ENDDO
    invar(:, :) = rdummy(:, :)

  END SUBROUTINE NSflip
  !=============================================================================
 SUBROUTINE get_land_index(nlon, nlat)
  !
  ! fill the index variable 'landpt%ilat, landpt%ilon'
  !
  ! Input variables:
  !   nlon           - # longitudes in input data set
  !   nlat           - # latitudes  in input data set
  !   npatch         - # patches in each grid from input data set
  !   inLon          - via cable_param_module
  !   inLat          - via cable_param_module
  !   longitude      - via cable_IO_vars_module, dim(mland), not patches
  !   latitude       - via cable_IO_vars_module, dim(mland), not patches
  !   nmetpatches    - via cable_IO_vars_module
  !   vegtype_metfile - via cable_IO_vars_module, dim(mland,nmetpatches)
  !   soiltype_metfile- via cable_IO_vars_module, dim(mland,nmetpatches)
  ! Output variables:
  !   max_vegpatches - via cable_IO_vars_module
  !   landpt%type    - via cable_IO_vars_module (%nap,cstart,cend,ilon,ilat)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlon, nlat

    ! local variables
    REAL :: lon2, distance, newLength
    INTEGER :: ii, jj, kk, tt, ncount

    ! range of longitudes from input file (inLon) should be -180 to 180,
    ! and longitude(:) has already been converted to -180 to 180 for CABLE.
    landpt(:)%ilon = -999
    landpt(:)%ilat = -999
    ncount = 0
    DO kk = 1, mland
       distance = 5300.0 ! initialise, units are degrees
       DO jj = 1, nlat
          DO ii = 1, nlon
             IF (inVeg(ii,jj, 1) > 0) THEN
                newLength = SQRT((inLon(ii) - longitude(kk))**2                      &
                     + (inLat(jj) -  latitude(kk))**2)
                IF (newLength < distance) THEN
                   distance = newLength
                   landpt(kk)%ilon = ii
                   landpt(kk)%ilat = jj
                END IF
             END IF
          END DO
       END DO

       IF (landpt(kk)%ilon < -900 .OR. landpt(kk)%ilat < -900) THEN
          PRINT *, 'Land point ', kk, ' cannot find the nearest grid!'
          PRINT *, 'lon, lat = ', longitude(kk), latitude(kk)
          PRINT *, 'inLon range:', MINVAL(inLon), MAXVAL(inLon)
          PRINT *, 'inLat range:', MINVAL(inLat), MAXVAL(inLat)
          STOP
       END IF

    END DO


  END SUBROUTINE get_land_index


 !=============================================================================


  SUBROUTINE countPatch(nlon, nlat, npatch)
  ! count the total number of active patches and
  ! fill the index variable 'landpt'
  !
  ! Input variables:
  !   nlon           - # longitudes in input data set
  !   nlat           - # latitudes  in input data set
  !   npatch         - # patches in each grid from input data set
  !   inLon          - via cable_param_module
  !   inLat          - via cable_param_module
  !   longitude      - via cable_IO_vars_module, dim(mland), not patches
  !   latitude       - via cable_IO_vars_module, dim(mland), not patches
  !   nmetpatches    - via cable_IO_vars_module
  !   vegtype_metfile - via cable_IO_vars_module, dim(mland,nmetpatches)
  !   soiltype_metfile- via cable_IO_vars_module, dim(mland,nmetpatches)
  ! Output variables:
  !   max_vegpatches - via cable_IO_vars_module
  !   landpt%type    - via cable_IO_vars_module (%nap,cstart,cend,ilon,ilat)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlon, nlat, npatch

    ! local variables
    REAL :: lon2, distance, newLength
    INTEGER :: ii, jj, kk, tt, ncount

    !set mpatch, the publically available number of patches
    mpatch = npatch

    ! range of longitudes from input file (inLon) should be -180 to 180,
    ! and longitude(:) has already been converted to -180 to 180 for CABLE.
    landpt(:)%ilon = -999
    landpt(:)%ilat = -999
    ncount = 0
    DO kk = 1, mland
       distance = 5300.0 ! initialise, units are degrees
       DO jj = 1, nlat
          DO ii = 1, nlon
             IF (inVeg(ii,jj, 1) > 0) THEN
                newLength = SQRT((inLon(ii) - longitude(kk))**2                      &
                     + (inLat(jj) -  latitude(kk))**2)
                IF (newLength < distance) THEN
                   distance = newLength
                   landpt(kk)%ilon = ii
                   landpt(kk)%ilat = jj
                END IF
             END IF
          END DO
       END DO

       IF (longitude(kk)-inLon(landpt(kk)%ilon) /= 0.) print *, "longitude: ", longitude(kk), &   ! MMY
                                                                "inLon: ", inLon(landpt(kk)%ilon) ! MMY
       IF (latitude(kk)-inLat(landpt(kk)%ilat) /= 0.) print *, "latitude: ", latitude(kk),    &   ! MMY
                                                               "inLat: ", inLat(landpt(kk)%ilat)  ! MMY

       IF (landpt(kk)%ilon < -900 .OR. landpt(kk)%ilat < -900) THEN
          PRINT *, 'Land point ', kk, ' cannot find the nearest grid!'
          PRINT *, 'lon, lat = ', longitude(kk), latitude(kk)
          PRINT *, 'inLon range:', MINVAL(inLon), MAXVAL(inLon)
          PRINT *, 'inLat range:', MINVAL(inLat), MAXVAL(inLat)
          STOP
       END IF

       landpt(kk)%nap = 0
       landpt(kk)%cstart = ncount + 1
       IF (ASSOCIATED(vegtype_metfile)) THEN
          DO tt = 1, nmetpatches
             IF (vegtype_metfile(kk,tt) > 0) ncount = ncount + 1
             landpt(kk)%nap = landpt(kk)%nap + 1
          END DO
          landpt(kk)%cend = ncount
          IF (landpt(kk)%cend < landpt(kk)%cstart) THEN
             PRINT *, 'Land point ', kk, ' does not have veg type!'
             PRINT *, 'landpt%cstart, cend = ', landpt(kk)%cstart, landpt(kk)%cend
             PRINT *, 'vegtype_metfile = ', vegtype_metfile(kk,:)
             STOP
          END IF
          ! CLN added for npatches
       ELSE IF ( npatch .GT. 1 .AND. cable_user%force_npatches_as .le. 0) THEN
          DO tt = 1, npatch
             IF (inVeg(landpt(kk)%ilon,landpt(kk)%ilat,tt) > 0) THEN
                landpt(kk)%nap = landpt(kk)%nap + 1
             ENDIF
          END DO
          ncount = ncount + landpt(kk)%nap
          landpt(kk)%cend = ncount
          IF (landpt(kk)%cend < landpt(kk)%cstart) THEN
             PRINT *, 'Land point ', kk, ' does not have veg type!'
             PRINT *, 'landpt%cstart, cend = ', landpt(kk)%cstart, landpt(kk)%cend
             PRINT *, 'vegtype_metfile = ', vegtype_metfile(kk,:)
             STOP
        END IF
          ! CLN added for npatches
       ELSE IF ( npatch .GT. 1 .AND. cable_user%force_npatches_as .eq. npatch) THEN
          !force the number of tiles per cell
          !will read in the tile info from the GWspatial subroutine
          landpt(kk)%nap = npatch
          ncount = ncount + landpt(kk)%nap
          landpt(kk)%cend = ncount
          IF (landpt(kk)%cend < landpt(kk)%cstart) THEN
             PRINT *, 'Land point ', kk, ' does not have veg type!'
             PRINT *, 'landpt%cstart, cend = ', landpt(kk)%cstart, landpt(kk)%cend
             PRINT *, 'vegtype_metfile = ', vegtype_metfile(kk,:)
             STOP
        END IF
     ELSE
        ! assume nmetpatches to be 1
        IF (nmetpatches == 1) THEN
           ncount = ncount + 1
           landpt(kk)%nap = 1
           landpt(kk)%cend = ncount
        ELSE
          PRINT *, 'nmetpatches = ', nmetpatches, '. Should be 1.'
          PRINT *, 'If soil patches exist, add new code.'
          STOP
        END IF
     END IF
  END DO
  ! CLN IF (ncount > mland * nmetpatches) THEN
  IF (ncount > mland * nmetpatches .AND. npatch == 1) THEN
     PRINT *, ncount, ' should not be greater than mland*nmetpatches.'
     PRINT *, 'mland, nmetpatches = ', mland, nmetpatches
     STOP
  END IF
  DEALLOCATE(inLon, inLat)

    ! Set the maximum number of active patches to that read from met file:
  max_vegpatches = MAXVAL(landpt(:)%nap)
  mpatch         = max_vegpatches
  !CLN    IF (max_vegpatches /= nmetpatches) THEN
    IF (max_vegpatches /= nmetpatches .and. npatch == 1) THEN
      PRINT *, 'Error! Met file claiming to have more active patches than'
      PRINT *, 'it really has. Check met file.'
      STOP
    END IF
    IF (npatch < nmetpatches) THEN
      PRINT *, 'Warning! Met file data have more patches than the global file.'
      PRINT *, 'Remember to check final veg type and patch fractions.'
    END IF

    ! Write to total # patches - used to allocate all of CABLE's variables:
    mp = ncount
    PRINT *, 'Total number of patches (countPatch): ', ncount

  END SUBROUTINE countPatch
  !=============================================================================
  SUBROUTINE write_default_params(met,  air,    ssnow, veg, bgc,               &
                                  soil, canopy, rough, rad, logn,              &
                                  vegparmnew, month, TFRZ, LUC_EXPT)
  ! Initialize many canopy_type, soil_snow_type, soil_parameter_type and
  ! roughness_type variables;
  ! Calculate 'froot' from 'rootbeta' parameter;
  ! Assign values from input file to their proper variables in soil_snow_type,
  ! soil_parameter_type, veg_parameter_type and patch_type;
  ! Prescribe parameters for each point based on its veg/soil type.
  !
  ! New input structure using netcdf and introduced 'month' to initialize
  ! soil profiles with the correct monthly average values (BP apr2010)
  !
  ! Input variables:
  !   longitude      - via cable_IO_vars_module, dim(mland), not patches
  !   latitude       - via cable_IO_vars_module, dim(mland), not patches
  !   nmetpatches    - via cable_IO_vars_module
  !   vegtype_metfile - via cable_IO_vars_module, dim(mland,nmetpatches)
  !   soiltype_metfile- via cable_IO_vars_module, dim(mland,nmetpatches)
  ! Output variables:
  !   max_vegpatches - via cable_IO_vars_module
  !   landpt(mp)%type- via cable_IO_vars_module (%nap,cstart,cend,ilon,ilat)
  !   patch(mp)%type - via cable_IO_vars_module (%frac,longitude,latitude)

    USE cable_common_module, only : vegin, soilin, &
                  calcsoilalbedo,cable_user,init_veg_from_vegin

    IMPLICIT NONE
    INTEGER,               INTENT(IN)    :: logn  ! log file unit number
    INTEGER,               INTENT(IN)    :: month ! month of year
    LOGICAL,                    INTENT(IN)    :: vegparmnew ! new format input
    REAL,                       INTENT(IN)    :: TFRZ
    TYPE (met_type),            INTENT(INOUT) :: met
    TYPE (air_type),            INTENT(INOUT) :: air
    TYPE (soil_snow_type),      INTENT(INOUT) :: ssnow
    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg
    TYPE (bgc_pool_type),       INTENT(INOUT) :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (canopy_type),         INTENT(INOUT)   :: canopy
    TYPE (roughness_type),      INTENT(INOUT)   :: rough
    TYPE (radiation_type),      INTENT(INOUT)   :: rad
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT

    INTEGER,dimension(:), ALLOCATABLE :: ALLVEG
    INTEGER :: e,f,h,i,klev  ! do loop counter
    INTEGER :: is     ! YP oct07
    INTEGER :: ir     ! BP sep2010
    REAL :: totdepth  ! YP oct07
    REAL :: tmp       ! BP sep2010
    REAL, ALLOCATABLE, DIMENSION(:,:) :: tmp_mp_ms

!    The following is for the alternate method to calculate froot by Zeng 2001
!    REAL :: term1(17), term2(17)                ! (BP may2010)
!    REAL :: aa(17), bb(17)   ! new parameter values for IGBP types
!    DATA aa /6.706,7.344,7.066,5.990,4.453,6.326,7.718,7.604,8.235,10.740,
!    10.740,5.558,5.558,5.558,4.372,4.372,4.372/
!    DATA bb /2.175,1.303,1.953,1.955,1.631,1.567,1.262,2.300,1.627,2.608,
!    2.608,2.614,2.614,2.614,0.978,0.978,0.978/
!    (BP may2010)

    ! *******************************************************************
    ! Site independent initialisations (all gridcells):
    canopy%cansto  = 0.0 ! canopy water storage (mm or kg/m2)
    canopy%sghflux = 0.0
    canopy%ghflux  = 0.0
    ssnow%ssdn   = 120.0 ! snow density per layer (kg/m3)
    ssnow%ssdnn  = 120.0 ! overall snow density (kg/m3)
    ssnow%tggsn  = tfrz  ! snow temperature per layer (K)
    ssnow%isflag = 0     ! snow layer scheme flag (0 = no/little snow, 1=snow)
    ssnow%snowd  = 0.0   ! snow liquid water equivalent depth (mm or kg/m2)
    ssnow%osnowd = 0.0   ! snow depth prev timestep (mm or kg/m2)
    ssnow%sdepth = 0.0   ! snow depth for each snow layer (BP jul2010)
    ssnow%snage  = 0.0   ! snow age
    ssnow%wbice  = 0.0   ! soil ice
    ssnow%thetai = 0.0   ! soil ice
    ssnow%smass  = 0.0   ! snow mass per layer (kg/m^2)
    ssnow%runoff = 0.0   ! runoff total = subsurface + surface runoff
    ssnow%rnof1  = 0.0   ! surface runoff (mm/timestepsize)
    ssnow%rnof2  = 0.0   ! deep drainage (mm/timestepsize)
    ssnow%rtsoil = 100.0 ! turbulent resistance for soil
    ssnow%wb_lake = 0.0
    canopy%ga     = 0.0  ! ground heat flux (W/m2)
    canopy%dgdtg  = 0.0  ! derivative of ground heat flux wrt soil temp
    canopy%fev    = 0.0  ! latent heat flux from vegetation (W/m2)
    canopy%fes    = 0.0  ! latent heat flux from soil (W/m2)
    canopy%fhs    = 0.0  ! sensible heat flux from soil (W/m2)
    !! vh_js !!
    canopy%us = 0.1 ! friction velocity (needed in roughness before first call to canopy: should in be retart?
    canopy%fh    = 0.0  ! sensible heat flux
    canopy%fe    = 0.0  ! sensible heat flux
    !mrd
    ssnow%qrecharge = 0.0
    ssnow%wtd = 1.0
    canopy%sublayer_dz = 0.01  !could go into restart to ensure starting/stopping runs gives identical results
                                !however the impact is negligible

   !IF(hide%Ticket49Bug2) THEN
    canopy%ofes    = 0.0  ! latent heat flux from soil (W/m2)
    canopy%fevc     = 0.0 !vh!
    canopy%fevw     = 0.0 !vh!
    canopy%fns      = 0.0
    canopy%fnv     = 0.0
    canopy%fhv     = 0.0
    canopy%fwsoil = 1.0 ! vh -should be calculated from soil moisture or
                        ! be in restart file

    ssnow%kth = 0.3  ! vh ! should be calculated from soil moisture or be in restart file
    ssnow%sconds(:,:) = 0.06_r_2    ! vh snow thermal cond (W m-2 K-1),
                                    ! should be in restart file

    ! parameters that are not spatially dependent

     select case(ms)

     case(6)
        soil%zse = (/.022, .058, .154, .409, 1.085, 2.872/) ! layer thicknessnov03
     case(12)
        soil%zse = (/.022,  0.0500,    0.1300 ,   0.3250 ,   0.3250 ,   0.3000,&
             0.3000,    0.3000 ,   0.3000,    0.3000,    0.7500,  1.50 /)
     case(13)
        soil%zse = (/.02,  0.0500,  0.06,  0.1300 ,   0.300 ,   0.300 ,   0.3000,&
             0.3000,    0.3000 ,   0.3000,    0.3000,    0.7500,  1.50 /)
     end select

    soil%zse_vec = real(spread(soil%zse,1,mp),r_2)

    rough%za_uv = 40.0 ! lowest atm. model layer/reference height
    rough%za_tq = 40.0

    veg%meth = 1 ! canopy turbulence parameterisation method: 0 or 1

    ! calculate vegin%froot from using rootbeta and soil depth
    ! (Jackson et al. 1996, Oceologica, 108:389-411)
    totdepth = 0.0
    DO is = 1, ms
       totdepth = totdepth + soil%zse(is) * 100.0  ! unit in centimetres
       vegin%froot(is, :) = MIN(1.0, 1.0-vegin%rootbeta(:)**totdepth)
    END DO
    DO is = ms, 2, -1
       vegin%froot(is, :) = vegin%froot(is, :)-vegin%froot(is-1, :)
    END DO

    ALLOCATE(defaultLAI(mp, 12))

    DO e = 1, mland ! over all land grid points

      ! Write to CABLE variables from temp variables saved in
      ! get_default_params
      veg%iveg(landpt(e)%cstart:landpt(e)%cend) =                              &
                          inVeg(landpt(e)%ilon, landpt(e)%ilat, 1:landpt(e)%nap)
      patch(landpt(e)%cstart:landpt(e)%cend)%frac =                            &
                        inPFrac(landpt(e)%ilon, landpt(e)%ilat, 1:landpt(e)%nap)

!write(*,*) 'iveg', e,  veg%iveg(landpt(e)%cstart:landpt(e)%cend)
!write(*,*) 'patchfrac', e,  patch(landpt(e)%cstart:landpt(e)%cend)%frac

      ! set land use (1 = primary; 2 = secondary, 3 = open)
      if (cable_user%popluc) then
         veg%iLU(landpt(e)%cstart:landpt(e)%cend)= 1
         if (landpt(e)%nap.gt.1) then
            veg%iLU(landpt(e)%cstart+1) = 2
            veg%iLU(landpt(e)%cend) = 3
         endif
      endif
      ! Check that patch fractions total to 1
      tmp = 0
      IF (landpt(e)%cstart == landpt(e)%cend) THEN
        patch(landpt(e)%cstart)%frac = 1.0
      ELSE
        DO is = landpt(e)%cstart, landpt(e)%cend
          tmp = tmp + patch(is)%frac
        END DO
        IF (ABS(1.0 - tmp) > 0.001) THEN
          IF ((1.0 - tmp) < -0.001 .OR. (1.0 - tmp) > 0.5) THEN
            PRINT *, 'Investigate the discrepancy in patch fractions:'
            PRINT *, 'patch%frac = ',                                          &
                                     patch(landpt(e)%cstart:landpt(e)%cend)%frac
            PRINT *, 'landpoint # ', e
            PRINT *, 'veg types = ', veg%iveg(landpt(e)%cstart:landpt(e)%cend)
            STOP
          END IF
          patch(landpt(e)%cstart)%frac = patch(landpt(e)%cstart)%frac + 1.0    &
                                         - tmp
        END IF
      END IF

      patch(landpt(e)%cstart:landpt(e)%cend)%longitude = longitude(e)
      patch(landpt(e)%cstart:landpt(e)%cend)%latitude  = latitude(e)
      soil%isoilm(landpt(e)%cstart:landpt(e)%cend) =                           &
                                          inSoil(landpt(e)%ilon, landpt(e)%ilat)
      ! Set initial soil temperature and moisture according to starting month
!! vh_js !!
    DO is = 1, ms
       ! Work around set everything above last input layer to the last input layer
       ssnow%tgg(landpt(e)%cstart:landpt(e)%cend, is) =                       &
              inTGG(landpt(e)%ilon,landpt(e)%ilat, min(is,size(inTGG,3)), month)
       ssnow%wb(landpt(e)%cstart:landpt(e)%cend, is) =                        &
              inWB(landpt(e)%ilon, landpt(e)%ilat, min(is,size(inTGG,3)), month)
    END DO

      ! Set initial snow depth and snow-free soil albedo


      DO is = 1, landpt(e)%cend - landpt(e)%cstart + 1  ! each patch
        DO ir = 1, nrb
           IF (CABLE_USER%POPLUC) then !vh! use same soilalbedo for all land-use tiles
              ssnow%albsoilsn(landpt(e)%cstart + is - 1, ir)                       &
                   = inALB(landpt(e)%ilon, landpt(e)%ilat, 1, ir) ! various rad band
           ELSE
                                 ! each band
              ssnow%albsoilsn(landpt(e)%cstart + is - 1, ir)                       &
                   = inALB(landpt(e)%ilon, landpt(e)%ilat, is, ir) ! various rad band
           ENDIF
        END DO
        ! total depth, change from m to mm
        ssnow%snowd(landpt(e)%cstart + is - 1)                                 &
                     = inSND(landpt(e)%ilon, landpt(e)%ilat, is, month) * 1000.0
      END DO

      ! Set default LAI values
      DO is = 1, 12
        defaultLAI(landpt(e)%cstart:landpt(e)%cend,is) =                       &
                                        inLAI(landpt(e)%ilon,landpt(e)%ilat,is)
      END DO

      ! Set IGBP soil texture values, Q.Zhang @ 12/20/2010.
      IF (soilparmnew) THEN

      soil%swilt(landpt(e)%cstart:landpt(e)%cend) =                            &
                                         inswilt(landpt(e)%ilon, landpt(e)%ilat)
      soil%sfc(landpt(e)%cstart:landpt(e)%cend) =                              &
                                           insfc(landpt(e)%ilon, landpt(e)%ilat)
      soil%ssat(landpt(e)%cstart:landpt(e)%cend) =                             &
                                          inssat(landpt(e)%ilon, landpt(e)%ilat)
      soil%bch(landpt(e)%cstart:landpt(e)%cend) =                              &
                                           inbch(landpt(e)%ilon, landpt(e)%ilat)
      soil%hyds(landpt(e)%cstart:landpt(e)%cend) =                             &
                                          inhyds(landpt(e)%ilon, landpt(e)%ilat)
      soil%sucs(landpt(e)%cstart:landpt(e)%cend) =                             &
                                     -1.* abs(insucs(landpt(e)%ilon, landpt(e)%ilat)) !ensure negative
      soil%rhosoil(landpt(e)%cstart:landpt(e)%cend) =                          &
                                       inrhosoil(landpt(e)%ilon, landpt(e)%ilat)
      soil%css(landpt(e)%cstart:landpt(e)%cend) =                              &
                                           incss(landpt(e)%ilon, landpt(e)%ilat)
      soil%cnsd(landpt(e)%cstart:landpt(e)%cend) =                             &
                                          incnsd(landpt(e)%ilon, landpt(e)%ilat)
      soil%sand(landpt(e)%cstart:landpt(e)%cend) =                             &
                                          insand(landpt(e)%ilon, landpt(e)%ilat)
      soil%clay(landpt(e)%cstart:landpt(e)%cend) =                             &
                                          inclay(landpt(e)%ilon, landpt(e)%ilat)
      soil%silt(landpt(e)%cstart:landpt(e)%cend) =                             &
                                          insilt(landpt(e)%ilon, landpt(e)%ilat)
      ENDIF

      ! vars intro for Ticket #27
      IF (calcsoilalbedo) THEN
        soil%soilcol(landpt(e)%cstart:landpt(e)%cend) =                        &
                                     inSoilColor(landpt(e)%ilon, landpt(e)%ilat)
      END IF

! offline only below
       ! If user defined veg types are present in the met file then use them.
       ! This means that if met file just has veg type and no other parameters,
       ! the other veg parameters will be chosen as a function of this type:
       ! N.B. for offline run only
       IF(ASSOCIATED(vegtype_metfile)) THEN ! i.e. iveg found in the met file
          ! Overwrite iveg for those patches available in met file,
          ! which are currently set to def values above:
          veg%iveg(landpt(e)%cstart:landpt(e)%cstart + nmetpatches - 1) =      &
                                                           vegtype_metfile(e, :)
          ! In case gridinfo file provides more patches than met file(BP may08)
          DO f = nmetpatches+1, landpt(e)%nap
             IF (patch(landpt(e)%cstart + f - 1)%frac > 0.0) THEN
                patch(landpt(e)%cstart)%frac = patch(landpt(e)%cstart)%frac    &
                                          + patch(landpt(e)%cstart + f - 1)%frac
                patch(landpt(e)%cstart + f - 1)%frac = 0.0
             END IF
          END DO
       END IF
       ! Similarly, if user defined soil types are present then use them:
       IF(ASSOCIATED(soiltype_metfile)) THEN ! i.e. isoil found in the met file
          soil%isoilm(landpt(e)%cstart:landpt(e)%cstart + nmetpatches - 1) =   &
                                                          soiltype_metfile(e, :)
       END IF
! offline only above
       !call veg% init that is common
       CALL init_veg_from_vegin(landpt(e)%cstart, landpt(e)%cend, veg)

       ! Prescribe parameters for current gridcell based on veg/soil type (which
       ! may have loaded from default value file or met file):
       DO h = landpt(e)%cstart, landpt(e)%cend ! over each patch in current grid
          bgc%cplant(h,:) = vegin%cplant(:, veg%iveg(h))
          bgc%csoil(h,:)  = vegin%csoil(:, veg%iveg(h))
          bgc%ratecp(:)   = vegin%ratecp(:, veg%iveg(h))
          bgc%ratecs(:)   = vegin%ratecs(:, veg%iveg(h))

          IF (.NOT. soilparmnew) THEN   ! Q,Zhang @ 12/20/2010
            soil%swilt(h)   =  soilin%swilt(soil%isoilm(h))
            soil%sfc(h)     =  soilin%sfc(soil%isoilm(h))
            soil%ssat(h)    =  soilin%ssat(soil%isoilm(h))
            soil%bch(h)     =  soilin%bch(soil%isoilm(h))
            soil%hyds(h)    =  soilin%hyds(soil%isoilm(h))
            soil%sucs(h)    =  soilin%sucs(soil%isoilm(h))
            soil%rhosoil(h) =  soilin%rhosoil(soil%isoilm(h))
            soil%css(h)     =  soilin%css(soil%isoilm(h))

            soil%silt(h)    =  soilin%silt(soil%isoilm(h))
            soil%clay(h)    =  soilin%clay(soil%isoilm(h))
            soil%sand(h)    =  soilin%sand(soil%isoilm(h))

             !MDeck
            IF (.not.cable_user%GW_MODEL) THEN
                    WRITE(logn,*) 'DO NOT USE GW_MODEL and soilparmen=.false.!'
                    WRITE(logn,*) 'IF YOU REALLY WANT THIS, CHANGE THE CODE!'
                    STOP
            END IF

            do klev=1,ms
               soil%clay_vec(h,klev) = real(soilin%clay(soil%isoilm(h)),r_2)
               soil%sand_vec(h,klev) = real(soilin%sand(soil%isoilm(h)),r_2)
               soil%silt_vec(h,klev) = real(soilin%silt(soil%isoilm(h)),r_2)
               soil%rhosoil_vec(h,klev) = real(soilin%rhosoil(soil%isoilm(h)),r_2)
               soil%watr(h,klev)    = 0.01
            end do
            soil%GWsucs_vec(h)  = real(abs(soilin%sucs(soil%isoilm(h)))*1000.0,r_2)
            soil%GWhyds_vec(h)   = real(soilin%hyds(soil%isoilm(h))*1000.0,r_2)
            soil%GWbch_vec(h)  = real(soilin%bch(soil%isoilm(h)),r_2)
            soil%GWssat_vec(h)  = real(soilin%ssat(soil%isoilm(h)),r_2)
            soil%GWwatr(h)    = 0.01

          END IF
          rad%latitude(h) = latitude(e)
            !IF(hide%Ticket49Bug4) &
          rad%longitude(h) = longitude(e)
          !jhan:is this done online? YES
          veg%ejmax(h) = 2.0 * veg%vcmax(h)
       END DO ! over each veg patch in land point
    END DO ! over all land points
    soil%albsoil = ssnow%albsoilsn

    !should not be needed
    !soil%elev = max(0.0,min(1.0e10,soil%elev))
    !soil%slope = max(0.000001,min(0.9,soil%slope))
    !soil%slope_std = max(0.0005,min(0.01,soil%slope_std))
    !soil%GWdz = max(7.0,min(150.0,soil%GWdz))
    !soil%drain_dens = max(1.0e-8,min(0.01,soil%drain_dens))

    ! check tgg and alb
    IF(ANY(ssnow%tgg > 350.0) .OR. ANY(ssnow%tgg < 180.0))                     &
            ! PRINT *, "sand are", soil%sand ! MMY
            ! PRINT *, "tgg are", ssnow%tgg ! MMY
            CALL abort('Soil temps nuts')
    IF(ANY(ssnow%albsoilsn > 1.0) .OR. ANY(ssnow%albsoilsn < 0.0))             &
           CALL abort('Albedo nuts')

    if (cable_user%GSWP3) then
       rough%za_uv = 2.0 + veg%hc ! lowest atm. model layer/reference height
       rough%za_tq = 2.0 + veg%hc
    end if
    ! Deallocate temporary variables:
    IF (soilparmnew) DEALLOCATE(inswilt, insfc, inssat, inbch, inhyds,         &
                       insucs, inrhosoil, incss, incnsd) ! Q,Zhang @ 12/20/2010
    IF (calcsoilalbedo) DEALLOCATE(inSoilColor) ! vars intro for Ticket #27
    DEALLOCATE(inVeg, inPFrac, inSoil, inWB, inTGG)
    DEALLOCATE(inLAI, inSND, inALB)
!    DEALLOCATE(soiltemp_temp,soilmoist_temp,patchfrac_temp,isoilm_temp,&
!         frac4_temp,iveg_temp)
!    IF(ASSOCIATED(vegtype_metfile)) DEALLOCATE(vegtype_metfile)
!    IF(ASSOCIATED(soiltype_metfile)) DEALLOCATE(soiltype_metfile)
    DEALLOCATE(soilin%silt, soilin%clay, soilin%sand, soilin%swilt,            &
               soilin%sfc, soilin%ssat, soilin%bch, soilin%hyds, soilin%sucs,  &
               soilin%rhosoil, soilin%css, vegin%canst1, vegin%dleaf,          &
               vegin%vcmax, vegin%ejmax, vegin%hc, vegin%xfang, vegin%rp20,    &
               vegin%rpcoef, vegin%rs20, vegin%shelrb, vegin%frac4,            &
               vegin%wai, vegin%vegcf, vegin%extkn, vegin%tminvj,              &
               vegin%tmaxvj, vegin%vbeta,vegin%clitt, vegin%zr, vegin%rootbeta, vegin%froot,         &
               vegin%cplant, vegin%csoil, vegin%ratecp, vegin%ratecs,          &
               vegin%xalbnir, vegin%length, vegin%width,                       &
               vegin%g0, vegin%g1,                                             &
               vegin%a1gs, vegin%d0gs, vegin%alpha, vegin%convex, vegin%cfrd,  &
               vegin%gswmin, vegin%conkc0,vegin%conko0,vegin%ekc,vegin%eko   )
    !         vegf_temp,urbanf_temp,lakef_temp,icef_temp, &

    if (allocated(inORG     )) deallocate(inORG)

    ! BP dec 2007
    !    IF (.NOT. vegparmnew) THEN
    veg%deciduous = .FALSE.
    IF (mvtype == 13) THEN
      WHERE (veg%iveg == 2 .OR. veg%iveg == 5) veg%deciduous = .TRUE.
    ELSE IF (mvtype == 15 .or. mvtype == 16 .or. mvtype == 17) THEN
      WHERE (veg%iveg == 3 .OR. veg%iveg == 4) veg%deciduous = .TRUE.
    ELSE
      STOP 'Warning. Check number of vegetation types.'
    END IF
    !    END IF

    ! Only the following snow inits are necessary,
    ! soilsnow will update other variables.
    WHERE(ssnow%snowd(:) > 0.0) ! in cm
      ssnow%ssdnn(:)    = 120.0 ! overall snow density (kg/m3)
      ssnow%ssdn(:, 1)  = 120.0 ! snow density per layer (kg/m3)
      ssnow%ssdn(:, 2)  = 120.0 ! snow density per layer (kg/m3)
      ssnow%ssdn(:, 3)  = 120.0 ! snow density per layer (kg/m3)
      ssnow%snage(:)    = 0.0   ! snow age (fresh)
      ssnow%isflag(:)   = 0
    ELSEWHERE
      ssnow%ssdnn(:)    = 140.0 ! overall snow density (kg/m3)
      ssnow%osnowd(:)   = 0.0   ! snow depth prev timestep (mm or kg/m2)
      ssnow%snage(:)    = 0.0   ! snow age
      ssnow%isflag(:)   = 0     ! snow layer scheme flag
                               ! (0 = no/little snow, 1=snow)
      ssnow%tggsn(:, 1) = 273.1 ! snow temperature per layer (K)
      ssnow%tggsn(:, 2) = 273.1 ! snow temperature per layer (K)
      ssnow%tggsn(:, 3) = 273.1 ! snow temperature per layer (K)
      ssnow%ssdn(:, 1)  = 140.0 ! snow density per layer (kg/m3)
      ssnow%ssdn(:, 2)  = 140.0 ! snow density per layer (kg/m3)
      ssnow%ssdn(:, 3)  = 140.0 ! snow density per layer (kg/m3)
      ssnow%smass(:, 1) = 0.0   ! snow mass per layer (kg/m^2)
      ssnow%smass(:, 2) = 0.0   ! snow mass per layer (kg/m^2)
      ssnow%smass(:, 3) = 0.0   ! snow mass per layer (kg/m^2)
    ENDWHERE
    ! Soil ice:

   WHERE(ssnow%tgg(:, :) < 273.16)
     ssnow%wbice(:,:) = ssnow%wb(:, :) * 0.8
   ELSEWHERE
     ssnow%wbice(:, :) = 0.0
   END WHERE

   ssnow%Qrecharge = 0.0
   ssnow%rtevap_sat = 0.0
   ssnow%rtevap_unsat = 0.0
   ssnow%satfrac = 0.5
   ssnow%wbliq = ssnow%wb - ssnow%wbice
   ssnow%GWwb = soil%GWssat_vec

   ssnow%wb_hys = -1.0e+36
   ssnow%hys_fac = 1.0
   ssnow%watr_hys = soil%watr
   ssnow%ssat_hys = soil%ssat_vec
   ssnow%smp_hys  = -1.0e+36

   !IF(hide%Ticket49Bug5) THEN

!! vh_js !! neeed to remove this if to enable the code below

      ! SLI specific initialisations:
    !  IF(cable_user%SOIL_STRUC=='sli') THEN
         ssnow%h0(:) = 0.0
         ssnow%S(:,:) = ssnow%wb(:,:)/SPREAD(soil%ssat,2,ms)
         ssnow%snowliq(:,:) = 0.0
         ssnow%Tsurface = 25.0
         ssnow%nsnow = 0
         ssnow%Tsoil = ssnow%tgg - 273.15
         ssnow%thetai = 0.0
         soil%zeta = 0.0
         soil%fsatmax = 0.0
   !   END IF

      IF(cable_user%SOIL_STRUC=='sli') THEN
         soil%nhorizons = 1 ! use 1 soil horizon globally
        ! veg%clitt = 5.0 ! (tC / ha)
         veg%F10 = 0.85
         veg%ZR = 5.0
      END IF

      IF(cable_user%SOIL_STRUC=='sli'.or.cable_user%FWSOIL_SWITCH=='Haverd2013') THEN
         veg%gamma = 3.e-2
         !veg%clitt = 5.0 ! (tC / ha)
      ENDIF
!! vh_js !!
      IF(cable_user%CALL_POP) THEN
         veg%disturbance_interval = 100
         veg%disturbance_intensity = 0.
      ENDIF

   !set vectorized versions as same as defaut for now
!   soil%swilt_vec(:,:)  = real(spread(soil%swilt(:),2,ms),r_2)
!   soil%sfc_vec(:,:)  = real(spread(soil%sfc(:),2,ms),r_2)
!   soil%sucs_vec(:,:)  = real(spread(soil%sucs(:),2,ms),r_2)
!   soil%bch_vec(:,:)  = real(spread(soil%bch(:),2,ms),r_2)
!   soil%ssat_vec(:,:)  = real(spread(soil%ssat(:),2,ms),r_2)
!   soil%hyds_vec(:,:)  = real(spread(soil%hyds(:),2,ms),r_2)
!   soil%css_vec(:,:)  = spread(soil%css(:),2,ms)
!   soil%cnsd_vec(:,:)  = spread(soil%cnsd(:),2,ms)
!   IF (.not.cable_user%gw_model) then
!   soil%rhosoil_vec(:,:) = spread(soil%rhosoil,2,ms)
!   ELSE
!      soil%sucs_vec(:,:)  = abs(soil%sucs_vec(:,:) )
!   ENDIF
!
  END SUBROUTINE write_default_params
  !=============================================================================
  SUBROUTINE write_cnp_params(veg, casaflux, casamet)
  ! Input variables:
  !   landpt(mp)%type- via cable_IO_vars_module (%cstart,cend,ilon,ilat)
  !   patch(mp)%type - via cable_IO_vars_module (%frac)
  !   inSorder       - via cable_param_module
  !   inArea         - via cable_param_module
  !   inNdep         - via cable_param_module
  !   inNfix         - via cable_param_module
  !   inPdust        - via cable_param_module
  !   inPwea         - via cable_param_module

    USE casaparm, ONLY: cropland, croplnd2
    IMPLICIT NONE
    TYPE (veg_parameter_type),  INTENT(IN)    :: veg
    TYPE (casa_flux),           INTENT(INOUT) :: casaflux
    TYPE (casa_met),            INTENT(INOUT) :: casamet

    ! local variables
    INTEGER :: ee, hh

    DO ee=1, mland ! over all land grid points
      casamet%isorder(landpt(ee)%cstart:landpt(ee)%cend) =                     &
                                       inSorder(landpt(ee)%ilon,landpt(ee)%ilat)
      DO hh = landpt(ee)%cstart, landpt(ee)%cend  ! each patch in current grid
        casamet%lon(hh) = patch(hh)%longitude
        casamet%lat(hh) = patch(hh)%latitude
        casamet%areacell(hh) = patch(hh)%frac                                  &
                               * inArea(landpt(ee)%ilon, landpt(ee)%ilat)
        casaflux%Nmindep(hh) = patch(hh)%frac                                  &
                               * inNdep(landpt(ee)%ilon, landpt(ee)%ilat)
        casaflux%Nminfix(hh) = patch(hh)%frac                                  &
                               * inNfix(landpt(ee)%ilon, landpt(ee)%ilat)
        casaflux%Pdep(hh)    = patch(hh)%frac                                  &
                               * inPdust(landpt(ee)%ilon, landpt(ee)%ilat)
        casaflux%Pwea(hh)    = patch(hh)%frac                                  &
                               * inPwea(landpt(ee)%ilon, landpt(ee)%ilat)
        !! vh !! fluxes shouldn't be weighted by patch frac.
     !   IF (CABLE_USER%POPLUC) then
           casaflux%Nmindep(hh) =  inNdep(landpt(ee)%ilon, landpt(ee)%ilat)
           casaflux%Nminfix(hh) = max( inNfix(landpt(ee)%ilon, landpt(ee)%ilat), &
                8.0e-4)
!vh ! minimum fixation rate of 3 kg N ha-1y-1 (8e-4 g N m-2 d-1)
! Cleveland, Cory C., et al. "Global patterns of terrestrial biological nitrogen (N2) &
!fixation in natural ecosystems." Global biogeochemical cycles 13.2 (1999): 623-645.
           casaflux%Pdep(hh)    = inPdust(landpt(ee)%ilon, landpt(ee)%ilat)
           casaflux%Pwea(hh)    = inPwea(landpt(ee)%ilon, landpt(ee)%ilat)
      !  ENDIF

        ! fertilizer addition is included here
        IF (veg%iveg(hh) == cropland .OR. veg%iveg(hh) == croplnd2) then
          ! P fertilizer =13 Mt P globally in 1994
          casaflux%Pdep(hh)    = casaflux%Pdep(hh)                             &
                                 + patch(hh)%frac * 0.7 / 365.0
          casaflux%Nmindep(hh) = casaflux%Nmindep(hh)                          &
                                 + patch(hh)%frac * 4.0 / 365.0
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(inSorder, inArea, inNdep, inNfix, inPwea, inPdust)

!write(668,*) 'in write_cnp_params: Ndep, Nfix', casaflux%Nmindep(1), casaflux%Nminfix(1)

  END SUBROUTINE write_cnp_params
  !============================================================================
  SUBROUTINE derived_parameters(soil, sum_flux, bal, ssnow, veg, rough)
    ! Gives values to parameters that are derived from other parameters.
    TYPE (soil_snow_type),      INTENT(INOUT)    :: ssnow
    TYPE (veg_parameter_type),  INTENT(IN)    :: veg
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (sum_flux_type),       INTENT(INOUT) :: sum_flux
    TYPE (balances_type),       INTENT(INOUT) :: bal
    TYPE (roughness_type),      INTENT(INOUT) :: rough

    INTEGER :: j,i,klev,k ! do loop counter
    REAL(r_2)    :: temp(mp)
    REAL    :: tmp2(mp)

    REAL(r_2), dimension(mp,ms) :: perc_frac
    !REAL(r_2), DIMENSION(17)    :: psi_o,psi_c
    REAL(r_2), DIMENSION(mp,ms) :: psi_tmp
    REAL(r_2), DIMENSION(mp,ms) :: soil_depth,rhosoil_temp
    REAL(r_2), DIMENSION(:,:), ALLOCATABLE :: ssat_bounded,rho_soil_bulk

    where(veg%iveg .eq. 17) soil%isoilm = 9

    soil_depth(:,1) = soil%zse_vec(:,1)
    do klev=2,ms
       soil_depth(:,klev) = soil_depth(:,klev-1) + soil%zse_vec(:,klev)
    end do

    !psi_o(1:3)  = -66000._r_2
    !psi_o(4)    = -35000._r_2
    !psi_o(5)    = -83000._r_2
    !psi_o(6:17) = -74000._r_2
    !psi_c(1:3)  = -2550000._r_2
    !psi_c(4)    = -2240000._r_2
    !psi_c(5)    = -4280000._r_2
    !psi_c(6:17) = -2750000._r_2
    ! Construct derived parameters and zero initialisations,
    ! regardless of where parameters and other initialisations
    ! have loaded from:
    soil%zshh(1) = 0.5 * soil%zse(1) ! distance between consecutive layer
                                     ! midpoints:
    soil%zshh(ms + 1) = 0.5 * soil%zse(ms)
    soil%zshh(2:ms)   = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))


    IF (cable_user%GW_MODEL) then

       soil%qhz_max(:) = real(gw_params%MaxHorzDrainRate,r_2)  !enable distributed values
       soil%hkrz(:)    = real(gw_params%hkrz,r_2)
       soil%zdepth(:)  = real(gw_params%zdepth,r_2)

       !placeholder
       soil%srf_frac_ma(:) = 0._r_2
       soil%edepth_ma(:)   = 0._r_2

       soil%qhz_efold(:) = real(gw_params%EfoldHorzDrainScale,r_2)*soil%drain_dens(:)&
                                   +   real(gw_params%EfoldHorzDrainRate,r_2)

       !include organin impact.  fraction of grid cell where percolation through
       !organic macropores dominates
       soil%org_vec = max(0._r_2,soil%org_vec)
       soil%org_vec = min(1._r_2,soil%org_vec)

       DO klev=1,ms
          do i=1,mp
             if (abs(soil%sand_vec(i,klev) + soil%clay_vec(i,klev) +&
                soil%silt_vec(i,klev)-1.0) .gt. 0.1) then
                soil%sand_vec(i,klev) = 0.4
                soil%clay_vec(i,klev) = 0.2
                soil%silt_vec(i,klev) = 0.4
             endif
          END DO
       END DO

       !use single var (uni) or two var (multi) regression from COsby 1984
       if (gw_params%cosby_univariate .or. &
           gw_params%cosby_multivariate .or. &
           gw_params%HC_SWC ) THEN

           rhosoil_temp(:,:) = (1.0-soil%org_vec(:,:))* 2.7*(1.0-soil%ssat_vec(:,:))+&
                          soil%org_vec(:,:)*gw_params%org%ssat_vec

          DO klev=1,ms;do i=1,mp
            if (gw_params%cosby_univariate) then
                soil%hyds_vec(i,klev) = 0.0070556*10.0**(-0.884 + 1.53*soil%sand_vec(i,klev))* &
                                         exp(-soil%hkrz(i)*(soil_depth(i,klev)-soil%zdepth(i)))
                soil%sucs_vec(i,klev) = 10.0 * 10.0**(1.88 -1.31*soil%sand_vec(i,klev))
                soil%bch_vec(i,klev) = 2.91 + 15.9*soil%clay_vec(i,klev)
                soil%ssat_vec(i,klev) = min(0.489,max(0.1, 0.489 - 0.126*soil%sand_vec(i,klev) ) )
                !forgot source but not from cosby and not for BC characteristic function
                soil%watr(i,klev) = 0.02 + 0.018*soil%clay_vec(i,klev) !forgot
              !!!2 parameters
                 soil%wbc_vec(i,klev) = 0.0
                 soil%smpc_vec(i,klev) = 0.0
            elseif (gw_params%cosby_multivariate) then
                soil%hyds_vec(i,klev) = 0.00706*(10.0**(-0.60 + 1.26*soil%sand_vec(i,klev) + &
                                                        -0.64*soil%clay_vec(i,klev) ) )*&
                                         exp(-soil%hkrz(i)*(soil_depth(i,klev)-soil%zdepth(i)))
                soil%sucs_vec(i,klev) = 10.0 * 10.0**(1.54 - 0.95*soil%sand_vec(i,klev) + &
                                                          0.63*soil%silt_vec(i,klev) )
                soil%bch_vec(i,klev) = 3.1 + 15.4*soil%clay_vec(i,klev) -  &
                                               0.3*soil%sand_vec(i,klev)
                soil%ssat_vec(i,klev) = 0.505 - 0.142*soil%sand_vec(i,klev) - &
                                                0.037*soil%clay_vec(i,klev)
                !forgot source but not from cosby and not for BC characteristic function
                soil%watr(i,klev) = 0.02 + 0.018*soil%clay_vec(i,klev)
                 soil%wbc_vec(i,klev) = 0.0
                 soil%smpc_vec(i,klev) = 0.0

              elseif (gw_params%HC_SWC) THEN
                !Hutson-Cass SWC : seperate dry/wet
                !avoid discont in derv at smp=sucs
                !pedo transfer from T. Mayr, N.J. JarÕisr Geoderma 91
                ! 1999
                soil%sucs_vec(i,klev) = 10.0 * 10.0** ( -4.98403 +&
                                         5.0922*soil%sand_vec(i,klev) +            &
                                         15.751*soil%silt_vec(i,klev) +            &
                                         0.124090*rhosoil_temp(i,klev) -     &
                                         16.4000*soil%org_vec(i,klev) -            &
                                         21.76*(soil%silt_vec(i,klev)**2.0) +      &
                                         14.382*(soil%silt_vec(i,klev)**3.0) +     &
                                         8.0407*(soil%clay_vec(i,klev)**2.0) +     &
                                         44.06*(soil%org_vec(i,klev)**2.0) )

                soil%bch_vec(i,klev) = 10.0**(1.0 / (-0.84669 - &
                                            0.4680*soil%sand_vec(i,klev) &
                                            +.9246*soil%silt_vec(i,klev)  &
                                            -0.4543*rhosoil_temp(i,klev) &
                                            -0.04979*soil%org_vec(i,klev) &
                                            +3.2947*(soil%sand_vec(i,klev)**2.0) &
                                            -1.689*(soil%sand_vec(i,klev)**3.0) &
                                            +11.2*(soil%org_vec(i,klev)**3.0) ))

                soil%ssat_vec(i,klev) = 0.234597 +    &
                                       0.466142*soil%sand_vec(i,klev) + &
                                       0.88163*soil%silt_vec(i,klev) + &
                                       0.643386*soil%clay_vec(i,klev) - &
                                       0.3028160*rhosoil_temp(i,klev) + &
                                       0.179762*(soil%sand_vec(i,klev)**2.0) - &
                                       0.03134631*(soil%silt_vec(i,klev)**2.0)

                soil%hyds_vec(i,klev) = 0.00706*(10.0**(&
                                             -0.60 + 1.26*soil%sand_vec(i,klev)+&
                                             -0.64*soil%clay_vec(i,klev) ) )*&
                                         exp(-soil%hkrz(i)*&
                                 (soil_depth(i,klev)-soil%zdepth(i)))

                soil%watr(i,klev) = 0.0
                if (klev .eq. 1) soil%GWwatr(i) = 0.0

                else

                 soil%hyds_vec(i,klev) = soil%hyds_vec(i,klev) *&
                                           exp(-soil%hkrz(i)*(soil_depth(i,klev)-soil%zdepth(i)))

              end if

          end do; end do

          if (.not.gw_params%HC_SWC) then
             DO klev=1,ms  !0-23.3 cm, data really is to 30cm
                do i=1,mp
                   soil%hyds_vec(i,klev)  = (1.-soil%org_vec(i,klev))*soil%hyds_vec(i,klev) + &
                                                soil%org_vec(i,klev)*gw_params%org%hyds_vec*&
                                                 exp(-soil%hkrz(i)*(soil_depth(i,klev)-soil%zdepth(i)))
                   soil%sucs_vec(i,klev) = (1.-soil%org_vec(i,klev))*soil%sucs_vec(i,klev) + &
                                               soil%org_vec(i,klev)*gw_params%org%sucs_vec
                   soil%bch_vec(i,klev) = (1.-soil%org_vec(i,klev))*soil%bch_vec(i,klev) +&
                                              soil%org_vec(i,klev)*gw_params%org%bch_vec
                   soil%ssat_vec(i,klev) = (1.-soil%org_vec(i,klev))*soil%ssat_vec(i,klev) + &
                                               soil%org_vec(i,klev)*gw_params%org%ssat_vec
                   soil%watr(i,klev)   = (1.-soil%org_vec(i,klev))*soil%watr(i,klev) + &
                                             soil%org_vec(i,klev)*gw_params%org%watr
                END DO
             END DO
          end if
          !!vegetation dependent field capacity (point plants get stressed) and
          do klev=1,ms
             do i=1,mp
                if (soil%isoilm(i) .ne. 9 .and. veg%iveg(i) .le. 16) then

                   psi_tmp(i,klev) = abs(psi_c(veg%iveg(i)))

                   soil%swilt_vec(i,klev) = (ssnow%ssat_hys(i,klev)-ssnow%watr_hys(i,klev)) * &
                                            (psi_tmp(i,klev)/soil%sucs_vec(i,klev))&
                                             **(-1.0/soil%bch_vec(i,klev))+&
                                            ssnow%watr_hys(i,klev)
                   soil%sfc_vec(i,klev) = (gw_params%sfc_vec_hk/soil%hyds_vec(i,klev))&
                                           **(1.0/(2.0*soil%bch_vec(i,klev)+3.0)) *&
                                           (ssnow%ssat_hys(i,klev)-ssnow%watr_hys(i,klev)) + ssnow%watr_hys(i,klev)

                   !soil%swilt_vec(i,klev) = (soil%ssat_vec(i,klev)-soil%watr(i,klev)) * &
                   !                         (psi_tmp(i,klev)/soil%sucs_vec(i,klev))&
                   !                          **(-1.0/soil%bch_vec(i,klev))+&
                   !                         soil%watr(i,klev)
                   !soil%sfc_vec(i,klev) = (gw_params%sfc_vec_hk/soil%hyds_vec(i,klev))&
                   !                        **(1.0/(2.0*soil%bch_vec(i,klev)+3.0)) *&
                   !                        (soil%ssat_vec(i,klev)-soil%watr(i,klev)) + soil%watr(i,klev)

                   soil%swilt_vec(i,klev) = min(0.95*soil%sfc_vec(i,klev),soil%swilt_vec(i,klev))

                else

                   soil%swilt_vec(i,klev) = soil%swilt(i)
                   soil%sfc_vec(i,klev) = soil%sfc(i)

                end if
             end do
          end do
       ! _________________ MMY remove hyds change with depth ________________
       !ELSE

       !   DO klev=1,ms
       !       soil%hyds_vec(:,klev) = soil%hyds_vec(:,klev)*exp(-soil%hkrz(:)*(soil_depth(:,klev)-soil%zdepth(:)))
       !   END DO
       ! ___________________________________________________________________
       END IF  !use either uni or multi cosby transfer func

       !set the non-vectored values to srf value
       soil%sfc(:) = real(soil%sfc_vec(:,3))
       soil%swilt(:) = real(soil%swilt_vec(:,3))
       soil%css(:)   = soil%css_vec(:,3)
       soil%cnsd(:)   = soil%cnsd_vec(:,3)
       soil%rhosoil(:) = soil%rhosoil_vec(:,3)

       !convert the units back to what default uses and GW only uses the
       !vectored versions
       soil%hyds = real(soil%hyds_vec(:,1))/1000.0
       soil%sucs = real(soil%sucs_vec(:,1))/1000.0
       soil%ssat = real(soil%ssat_vec(:,1))
       soil%bch  = real(soil%bch_vec(:,1))
       soil%sand = real(soil%sand_vec(:,1))
       soil%clay = real(soil%clay_vec(:,1))
       soil%silt = real(soil%silt_vec(:,1))

       do i=1,mp
          soil%slope(i) = min(0.9,max(1e-5,soil%slope(i)))
          soil%slope_std(i) = min(0.9,max(1e-5,soil%slope_std(i)))
       end do

       if ((gw_params%MaxSatFraction .lt. -9999.9) .and. (mp .eq. 1)) soil%slope(:) = 0.01

    ELSE  !not gw model

      !These are not used when gw_model == false
      soil%watr = 0._r_2
      soil%GWwatr = 0._r_2

    END IF



    IF (cable_user%soil_thermal_fix) then
       IF ((.not.(gw_params%cosby_univariate.or.gw_params%cosby_multivariate)) .or. .not.cable_user%gw_model) THEN
          WRITE(logn,*) 'OVER WRITING CSS_VEC, RHOSOIL_VEC, CNSD_VEC read from gw_elev input'
          WRITE(logn,*) 'Forcing values consistant with soil_thermal_fix = true '
          WRITE(*,*)    'OVER WRITING CSS_VEC, RHOSOIL_VEC, CNSD_VEC read from gw_elev input'
          WRITE(*,*)    'Forcing values consistant with soil_thermal_fix = true '
       END IF

       if (allocated(ssat_bounded)) deallocate(ssat_bounded)
       if (allocated(rho_soil_bulk)) deallocate(rho_soil_bulk)

       allocate(ssat_bounded(size(soil%ssat_vec,dim=1),&
                            size(soil%ssat_vec,dim=2) ) )

       ssat_bounded(:,:) = min( 0.8, max(0.1, &
                                         soil%ssat_vec(:,:) ) )

       allocate(rho_soil_bulk(size(soil%rhosoil_vec,dim=1),&
                             size(soil%rhosoil_vec,dim=2) ) )

       rho_soil_bulk(:,:) = min(2500.0, max(500.0 , &
                             (2700.0*(1.0 - ssat_bounded(:,:)) ) ) )

       do klev=1,ms
          do i=1,mp


            ! if (soil%isoilm(i) .ne. 9 .and. veg%iveg(i) .lt. 16) then

                !soil%rhosoil_vec(i,klev) = 2700.0

                soil%cnsd_vec(i,klev) = ( (0.135*(1.0-ssat_bounded(i,klev))) +&
                                    (64.7/soil%rhosoil_vec(i,klev)) ) / &
                                  (1.0 - 0.947*(1.0-ssat_bounded(i,klev)))

                soil%rhosoil_vec(i,klev) = soil%rhosoil_vec(i,klev)/(1.0-soil%ssat_vec(i,klev))
                !took avg of results from A New Perspective on Soil Thermal Properties Ochsner, Horton,Tucheng
                !Soil Sci Soc America 2001
                !to find what silt (1.0-sand-clay) is !simply regress to his means !in J/kg/K
                 soil%css_vec(i,klev) =  max(910.6479*soil%silt_vec(i,klev) +&
                                          916.4438 * soil%clay_vec(i,klev) +&
                                          740.7491*soil%sand_vec(i,klev), 800.0)

           !  end if

          end do
       end do

       k=1
       do i=1,mp
          !if (soil%isoilm(i) .ne. 9) then
             soil%rhosoil(i) = soil%rhosoil_vec(i,1)
             soil%cnsd(i)    = soil%cnsd_vec(i,1)
             soil%css(i)     = soil%css_vec(i,1)
          !end if
       end do

       IF (cable_user%gw_model) then  !organic correction?
          do klev=1,ms
             do i=1,mp
                soil%css_vec(i,klev) = (1.-soil%org_vec(i,klev))*soil%css_vec(i,klev) + &
                                          real(soil%org_vec(i,klev)*gw_params%org%css_vec)
                soil%cnsd_vec(i,klev) = (1.-soil%org_vec(i,klev))*soil%cnsd_vec(i,klev) + &
                                          soil%org_vec(i,klev)*gw_params%org%cnsd_vec
              end do
           end do
      END IF
      if (allocated(ssat_bounded)) deallocate(ssat_bounded)
      if (allocated(rho_soil_bulk)) deallocate(rho_soil_bulk)

    ELSEIF ( .NOT. soilparmnew) THEN  ! Q,Zhang @ 12/20/2010
      soil%cnsd  = soil%sand * 0.3 + soil%clay * 0.25                          &
                   + soil%silt * 0.265 ! set dry soil thermal conductivity
                                       ! [W/m/K]
      soil%cnsd_vec = spread(soil%cnsd,2,ms)
    END IF

    soil%hsbh   = soil%hyds*ABS(soil%sucs) * soil%bch ! difsat*etasat
    soil%ibp2   = NINT(soil%bch) + 2
    ! Ticket #66
    where( soil%ssat > 0.) & ! Avoid divide by
    soil%pwb_min = (soil%swilt/soil%ssat)**soil%ibp2
    soil%i2bp3  = 2 * NINT(soil%bch) + 3
    rough%hruff = max(0.01, veg%hc - 1.2 * ssnow%snowd/max(ssnow%ssdnn, 100.))
    rough%hruff_grmx = rough%hruff
    ! owetfac introduced by EAK apr2009
    ssnow%owetfac = MAX(0.0, MIN(1.0,                                          &
                   (REAL(ssnow%wb(:, 1)) - soil%swilt) /                  &
                   (soil%sfc - soil%swilt)))
    temp(:) = 0.0
    tmp2(:) = 0.0
    WHERE ( ssnow%wbice(:, 1) > 0. ) ! Prevents divide by zero at glaciated
                                     ! points where wb and wbice=0.
      temp(:) = ssnow%wbice(:, 1) / ssnow%wb(:, 1)
      tmp2(:) = REAL(temp(:))
      ssnow%owetfac = ssnow%owetfac * (1.0 - tmp2(:)) ** 2

    END WHERE
    ssnow%pudsto = 0.0
    ssnow%pudsmx = 0.0

    ! Initialise sum flux variables:
    sum_flux%sumpn  = 0.0
    sum_flux%sumrp  = 0.0
    sum_flux%sumrpw = 0.0
    sum_flux%sumrpr = 0.0
    sum_flux%sumrs  = 0.0
    sum_flux%sumrd  = 0.0
    sum_flux%dsumpn = 0.0
    sum_flux%dsumrp = 0.0
    sum_flux%dsumrd = 0.0
    ! Initialise conservation variables:
    bal%precip_tot = 0.0
    bal%rnoff_tot  = 0.0
    bal%evap_tot   = 0.0
    bal%wbal_tot   = 0.0
    bal%ebal_tot   = 0.0
    bal%ebal_tot_cncheck = 0.0
    bal%drybal = 0.0
    bal%wetbal = 0.0
    bal%wbtot0 = 0.0
    bal%RadbalSum = 0.0
    DO j=1, ms
       bal%wbtot0 = bal%wbtot0 + REAL(ssnow%wb(:, j)) * soil%zse(j)       &
                    * 1000.0
    END DO
    bal%osnowd0 = ssnow%osnowd

  !! vh_js !! comment out hide% condition
   ! IF (hide%Ticket49Bug6) THEN

   IF(cable_user%SOIL_STRUC=='sli') THEN
      ! Only 1 horizon by default !
      soil%nhorizons = 1
      soil%ishorizon = 1
   END IF
   ! END IF

  where(soil%ssat_vec .gt. 0.0)
     ssnow%wblf = max(0.01_r_2,ssnow%wbliq/soil%ssat_vec)
     ssnow%wbfice = ssnow%wbice / soil%ssat_vec
  elsewhere
    ssnow%wblf =0.01
    ssnow%wbfice =0.99
  endwhere

  do k=1,ms
     do i=1,mp
        if (ssnow%wb_hys(i,k) .lt. 0._r_2) then
           ssnow%wb_hys(i,k) = ssnow%wb(i,k)
        end if
        ssnow%wb_hys(i,k)  = max(soil%watr(i,k) ,min(soil%ssat_vec(i,k), ssnow%wb_hys(i,k)))

       if (ssnow%smp_hys(i,k) .lt. -1.0e+30_r_2) then  !set to missing, calc
          ssnow%smp_hys(i,k) = -soil%sucs_vec(i,k)*  &
                               ( (ssnow%wb_hys(i,k)-ssnow%watr_hys(i,k))/&
                                 (ssnow%ssat_hys(i,k)-ssnow%watr_hys(i,k)) )**&
                                (-1._r_2/soil%bch_vec(i,k) )
       end if
       ssnow%smp_hys(i,k) = max(-1.0e10,min(-soil%sucs_vec(i,k),ssnow%smp_hys(i,k) ))
    end do
  end do

  if (cable_user%gw_model .and. gw_params%bc_hysteresis) then
      do klev=1,ms
         do i=1,mp
            if (soil%isoilm(i) .ne. 9 .and. veg%iveg(i) .le. 16) then

               psi_tmp(i,klev) = abs(psi_c(veg%iveg(i)))

               soil%swilt_vec(i,klev) = (ssnow%ssat_hys(i,klev)-ssnow%watr_hys(i,klev)) * &
                                        (psi_tmp(i,klev)/soil%sucs_vec(i,klev))&
                                         **(-1.0/soil%bch_vec(i,klev))+&
                                        ssnow%watr_hys(i,klev)
               soil%sfc_vec(i,klev) = (gw_params%sfc_vec_hk/soil%hyds_vec(i,klev))&
                                       **(1.0/(2.0*soil%bch_vec(i,klev)+3.0)) *&
                                       (ssnow%ssat_hys(i,klev)-ssnow%watr_hys(i,klev)) + ssnow%watr_hys(i,klev)
            end if
        end do
     end do

   end if
  END SUBROUTINE derived_parameters
  !============================================================================
  SUBROUTINE check_parameter_values(soil, veg, ssnow)
    ! Checks for basic inconsistencies in parameter values
    ! INH - changed soil & veg to INOUT.
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameter data
    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg   ! vegetation parameter
                                                       ! data
    TYPE (soil_snow_type),      INTENT(INOUT) :: ssnow ! soil and snow
                                                       ! variables
    INTEGER :: i, j ! do loop counter

    DO i = 1, mland
       ! Check all veg types make sense:
       IF(ANY(veg%iveg(landpt(i)%cstart:(landpt(i)%cstart + landpt(i)%nap      &
          - 1)) < 1 ) .OR. ANY(veg%iveg(landpt(i)%cstart:(landpt(i)%cstart +   &
          landpt(i)%nap - 1)) > mvtype)) THEN
          WRITE(*, *) 'SUBROUTINE load_parameters:'
          WRITE(*, *) 'Land point number:', i
          WRITE(*, *) 'Veg types:', veg%iveg(landpt(i)%cstart:                 &
                      (landpt(i)%cstart + landpt(i)%nap - 1))
          CALL abort('Unknown vegetation type! Aborting.')
       END IF
       ! Check all soil types make sense:
       IF(ANY(soil%isoilm(landpt(i)%cstart:(landpt(i)%cstart + landpt(i)%nap   &
          - 1)) < 1 ) .OR. ANY(soil%isoilm(landpt(i)%cstart:(landpt(i)%cstart  &
          + landpt(i)%nap - 1)) > mstype) .and. .not.soilparmnew) THEN
          WRITE(*,*) 'SUBROUTINE load_parameters: soil < 1 or > ',mstype
          DO j=landpt(i)%cstart,landpt(i)%cend
             IF (soil%isoilm(j) .lt. 1 .or. soil%isoilm(j) .gt. mstype) THEN
                write(*,*) 'AT ',i,j,' indices '
                write(*,*) '~~~~~~~~~~~~~~lat=>',&
                            patch(landpt(i)%cstart:landpt(i)%cend)%latitude
                write(*,*) '~~~~~~~~~~~~~~lon=>',&
                            patch(landpt(i)%cstart:landpt(i)%cend)%longitude
                write(*,*) '~~~~~~~~~~~~~~soil%isoilm',soil%isoilm(j)
             END IF
          ENDDO
          WRITE(*,*) 'Land point number:',i
          CALL abort('Unknown soil type! Aborting.')
       END IF
       ! Check patch fractions sum to 1 in each grid cell:
       IF((SUM(patch(landpt(i)%cstart:landpt(i)%cend)%frac) - 1.0)             &
          > 1.0E-6) THEN
          WRITE(*,*) 'SUBROUTINE load_parameters:'
          WRITE(*,*) 'At land point number', i
          WRITE(*,*) 'And patch numbers:  ', landpt(i)%cstart, landpt(i)%cend
          WRITE(*,*) 'patchfrac values are: ',                                 &
                      patch(landpt(i)%cstart:landpt(i)%cend)%frac
          WRITE(*,*) 'veg types are:        ',                                 &
                      veg%iveg(landpt(i)%cstart:landpt(i)%cend)
          WRITE(*,*) 'patch longitudes are: ',                                 &
                      patch(landpt(i)%cstart:landpt(i)%cend)%longitude
          WRITE(*,*) 'patch latitudes are:  ',                                 &
                      patch(landpt(i)%cstart:landpt(i)%cend)%latitude
          CALL abort ('Sum of fractional coverage of vegetation patches /= 1!')
       END IF
       !      ! Check sum of surface type fractions is 1:
       !      IF(landpt(i)%veg%frac + landpt(i)%urban%frac +                   &
       !         landpt(i)%lake%frac + landpt(i)%ice%frac /= 1) THEN
       !        WRITE(*,*) 'SUBROUTINE load_parameters:'
       !        WRITE(*,*) 'At land point number', i
       !        CALL abort ('Sum of fractional coverage of surface types /= 1!')
       !      END IF
    END DO
    ! Check sand+soil+clay fractions sum to 1:
    DO i = 1, mland
       DO j = 1, landpt(i)%nap
 ! vh changed limits from 1.0000001, 0.999999 to 1.01 and 0.99 for compatibility with gridinfo
          IF((soil%sand(landpt(i)%cstart + j - 1)                              &
              + soil%silt(landpt(i)%cstart + j - 1)                            &
              + soil%clay(landpt(i)%cstart + j - 1)) > 1.01 .OR.          &
             (soil%sand(landpt(i)%cstart + j - 1)                              &
              + soil%silt(landpt(i)%cstart + j - 1)                            &
              + soil%clay(landpt(i)%cstart + j - 1)) < 0.99) THEN
             WRITE(*,*) 'SUBROUTINE load_parameters:'
             WRITE(*,*) 'At land point number:', i
             WRITE(*,*) '        patch number:', j
             WRITE(*,*) 'Clay fraction is ',soil%clay(landpt(i)%cstart + j - 1)
             WRITE(*,*) 'Sand fraction is ',soil%sand(landpt(i)%cstart + j - 1)
             WRITE(*,*) 'Silt fraction is ',soil%silt(landpt(i)%cstart + j - 1)
             WRITE(*,*) 'Silt fraction is ',soil%isoilm(landpt(i)%cstart + j - 1)
             WRITE(*,*) 'SUM:',soil%sand(landpt(i)%cstart + j - 1)             &
                               + soil%silt(landpt(i)%cstart + j - 1)           &
                               + soil%clay(landpt(i)%cstart + j - 1)
             !mrd561 error where fraction was slightly off summing to 1.  Fix rather than abort.
             soil%silt(landpt(i)%cstart + j - 1) = 1.0 -                       &
                                         soil%clay(landpt(i)%cstart + j - 1) - &
                                         soil%sand(landpt(i)%cstart + j - 1)
          END IF
       END DO
    END DO
    ! Check that fraction of roots in each layer sum to 1:
    DO i = 1, mland
       DO j = 1, landpt(i)%nap
          IF(ABS(1 - SUM(veg%froot((landpt(i)%cstart + j - 1), :)))            &
             > 0.00001) THEN
             WRITE(*,*) 'SUBROUTINE load_parameters:'
             WRITE(*,*) 'At land point number:', i, 'patch:', j
             WRITE(*,*) 'Froot:',veg%froot((landpt(i)%cstart + j - 1), :)
             veg%froot((landpt(i)%cstart+j-1),ms) = veg%froot((landpt(i)%cstart+j-1),ms) + &
                      (1. - SUM(veg%froot((landpt(i)%cstart + j - 1), :)))
          END IF
       END DO
    END DO

    ! Check that wilting pt < field capacity < saturation value:
    IF(ANY(soil%swilt > soil%sfc) .OR. ANY(soil%sfc > soil%ssat)) THEN
       DO i = 1, mland
          DO j = 1, landpt(i)%nap
             IF(soil%swilt(landpt(i)%cstart + j - 1) >                         &
                soil%sfc(landpt(i)%cstart + j - 1)                             &
                .OR. soil%sfc(landpt(i)%cstart + j - 1) >                      &
                soil%ssat(landpt(i)%cstart + j - 1)) THEN
                WRITE(*, *) 'SUBROUTINE load_parameters:'
                WRITE(*, *) 'At land point number', i, 'patch:', j
                CALL abort ('Wilting pt < field capacity < saturation '//      &
                            'violated!')
             END IF
          END DO
       END DO
    END IF
    ! Ensure soil moisture values are reasonable (possible restart precision
    ! issue):
    !actually if denliq .ne. denice than ssnow%wb > ssat_vec is possible due to
    !the expansion during freezing
    !mrd561 Left using a real to set wb since that is what trunk does
    WHERE(ssnow%wb  > real(soil%ssat_vec)) ! Can only happen due to i/o issues
       ssnow%wb = 0.9999 * real(soil%ssat_vec)
    END WHERE

  END SUBROUTINE check_parameter_values
!===============================================================================
SUBROUTINE report_parameters(logn, soil, veg, bgc, rough,                    &
                               ssnow, canopy, casamet, casapool, casaflux,     &
                               phen, vegparmnew, verbose )

   USE cable_pft_params_mod, ONLY : veg_desc
   USE cable_soil_params_mod, ONLY : soil_desc
   IMPLICIT NONE
   INTEGER,      INTENT(IN)  :: logn        ! log file unit number
   LOGICAL,      INTENT(IN)  :: vegparmnew  ! are we using the new format?
   LOGICAL,      INTENT(IN)  :: verbose     ! write all parameter details to
                                            ! log file?
   TYPE (soil_parameter_type), INTENT(IN)  :: soil
   TYPE (veg_parameter_type),  INTENT(IN)  :: veg
   TYPE (bgc_pool_type),       INTENT(IN)  :: bgc
   TYPE (roughness_type),      INTENT(IN)  :: rough
   TYPE (soil_snow_type),      INTENT(IN)  :: ssnow
   TYPE (canopy_type),         INTENT(IN)  :: canopy
   TYPE (casa_met)     , INTENT(IN)  :: casamet
   TYPE (casa_pool)    , INTENT(IN)  :: casapool
   TYPE (casa_flux)    , INTENT(IN)  :: casaflux
   TYPE (phen_variable), INTENT(IN)  :: phen

   INTEGER :: e, f, g ! do loop counter
   CHARACTER(LEN=15) :: patchfmtr  ! patch format specifier for real numbers
   CHARACTER(LEN=13) :: patchfmti  ! patch format specifier for integer
                                   ! numbers
   CHARACTER(LEN=15) :: patchfmte  ! patch format specifier for expon. numbers
   CHARACTER(LEN=15) :: patchfmte2 ! patch format specifier for expon. numbers

   ! Get vegetation/soil type descriptions in case they haven't yet been
   ! loaded (i.e. if restart file + met file contains all parameter/init/LAI
   ! info). This will not overwrite any parameter values.
   ! CALL get_type_parameters(filename_veg, filename_soil, logn, vegparmnew)

   ! Only report parameters for active vegetation patches:
   IF (report_parameters_to_log) THEN
   DO e = 1, mland
      WRITE(logn, *) '==================================================='//  &
                     '======'
      WRITE(logn, '(A36, I8, 1X, A1)') ' CABLE setup details for land'//      &
            ' point ',e,':'
      WRITE(logn, *) '==================================================='//  &
                     '======'
      !      WRITE(logn,'(A21)') ' Surface type ratios:'
      !      WRITE(logn,*) '---------------------------------------------'//  &
      !                    '------------'
      !      ! Write surface type ratios to log file:
      !      WRITE(logn,'(A30,I3,A1)') '                   vegetated: ',&
      !           INT(landpt(e)%veg%frac*100.0),'%'
      !      WRITE(logn,'(A30,I3,A1)') '                       urban: ',&
      !           INT(landpt(e)%urban%frac*100.0),'%'
      !      WRITE(logn,'(A30,I3,A1)') '                       lakes: ',&
      !           INT(landpt(e)%lake%frac*100.0),'%'
      !      WRITE(logn,'(A30,I3,A1)') '                    land ice: ',&
      !           INT(landpt(e)%ice%frac*100.0),'%'
      !      ! Report patch details to log file:
      !      WRITE(logn,*) '---------------------------------------------'//  &
      !                    '------------'
      WRITE(logn, '(A43)') ' Proportions of each active veg/soil patch:'
      WRITE(logn, *) '---------------------------------------------------'//  &
                     '------'
      DO g = 1, landpt(e)%nap
         WRITE(logn, '(A7, I2, A3, F6.2, A11, I3, 1X, A30)') ' patch ',       &
               g,':  ', patch(landpt(e)%cstart + g - 1)%frac * 100.0,         &
               '% veg type ', veg%iveg(landpt(e)%cstart + g - 1),             &
               TRIM(veg_desc(veg%iveg(landpt(e)%cstart + g - 1)))
         WRITE(logn,'(18X, A11, I3, 1X, A45)') '  soil type',                 &
               soil%isoilm(landpt(e)%cstart + g - 1),                         &
               TRIM(soil_desc(soil%isoilm(landpt(e)%cstart + g - 1)))
      END DO


      IF(verbose) THEN
         ! Set up format specifier for writing active patch details below:
         WRITE(patchfmtr,'(A8, I1, A6)') '(4X,A50,', landpt(e)%nap, 'F12.4)'
         WRITE(patchfmti,'(A8, I1, A4)') '(4X,A50,', landpt(e)%nap, 'I12)'
         WRITE(patchfmte,'(A8, I1, A6)') '(4X,A50,', landpt(e)%nap, 'E12.4)'
         WRITE(patchfmte2,'(A8, I1, A6)') '(4X,A50,', landpt(e)%nap, 'E12.4)'
         ! Write parameter set details to log file:
         WRITE(logn, *) '------------------------------------------------'//  &
                        '---------'
         WRITE(logn, '(A36, I8, 1X, A2)') ' CABLE parameter values (land '//  &
               'point ', e, '):'
         WRITE(logn, *) '------------------------------------------------'//  &
                        '---------'
         WRITE(logn,'(4X, A50, 2F10.4)') 'reference height (m): ',            &
                               ! AJA MODIFIED
                               ! rough%za(e*max_vegpatches)
               rough%za_uv(landpt(e)%cend - landpt(e)%cstart + 1),            &
               rough%za_tq(landpt(e)%cend - landpt(e)%cstart + 1)
         WRITE(logn, *) ' Vegetation parameters: '
         WRITE(logn, patchfmti) 'Veg type for each active (>0% gridcell) '//  &
               'patch: ', veg%iveg(landpt(e)%cstart:landpt(e)%cend)
         WRITE(logn, patchfmtr) 'Vegetation height (m): ',                    &
               veg%hc(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'Fraction of roots in layer 1 (-): ',         &
               veg%froot(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 1)
         WRITE(logn, patchfmtr) 'Fraction of roots in layer 2 (-): ',         &
               veg%froot(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 2)
         WRITE(logn, patchfmtr) 'Fraction of roots in layer 3 (-): ',         &
               veg%froot(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 3)
         WRITE(logn, patchfmtr) 'Fraction of roots in layer 4 (-): ',         &
               veg%froot(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 4)
         WRITE(logn, patchfmtr) 'Fraction of roots in layer 5 (-): ',         &
               veg%froot(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 5)
         WRITE(logn, patchfmtr) 'Fraction of roots in layer 6 (-): ',         &
               veg%froot(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 6)
         WRITE(logn, patchfmtr) 'Fraction of plants which are C4 (-): ',      &
               veg%frac4(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmtr) 'Maximum canopy water storage (mm/LAI): ',    &
               veg%canst1(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap  &
               - 1))
         WRITE(logn, patchfmte)                                               &
               'Max pot elec transport rate top leaf (mol/m2/s): ',           &
               veg%ejmax(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmte)                                               &
               'Max RuBP carboxylation rate top leaf (mol/m^2/s): ',          &
               veg%vcmax(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmtr) 'Plant respiration coeff @ 20 C '//           &
               '(mol/m^2/s): ', veg%rp20(landpt(e)%cstart:(landpt(e)%cstart   &
               + landpt(e)%nap - 1))
         WRITE(logn, patchfmtr)                                               &
               'Temperature coef nonleaf plant respiration (1/C): ',          &
               veg%rpcoef(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap  &
               - 1))
         WRITE(logn, patchfmtr) 'Sheltering factor (-): ',                    &
               veg%shelrb(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap  &
               - 1))
         WRITE(logn, patchfmtr) 'Chararacteristic legnth of leaf (m): ',      &
               veg%dleaf(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmtr) 'Leaf angle parameter (-): ',                 &
               veg%xfang(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmtr)                                               &
               'Min temperature for start of photosynthesis (C): ',           &
               veg%tminvj(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap  &
               - 1))
         WRITE(logn, patchfmtr)                                               &
               'Max temperature for start of photosynthesis (C): ',           &
               veg%tmaxvj(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap  &
               - 1))
         WRITE(logn, patchfmtr) 'Stomatal sensitivity to soil water: ',       &
               veg%vbeta(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmtr) 'Modifier for surface albedo in near IR '//   &
               'band: ', veg%xalbnir(landpt(e)%cstart:(landpt(e)%cstart +     &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'a1 parameter in leaf stomatal model  ',      &
               veg%a1gs(landpt(e)%cstart:(landpt(e)%cstart +                  &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'd0 parameter in leaf stomatal model  ',      &
               veg%d0gs(landpt(e)%cstart:(landpt(e)%cstart +                  &
               landpt(e)%nap - 1))
         IF (icycle == 0) THEN
           WRITE(logn,'(4X, A50, F12.4)')                                     &
                 'Plant carbon rate constant pool 1 (1/year): ', bgc%ratecp(1)
           WRITE(logn,'(4X, A50, F12.4)')                                     &
                 'Plant carbon rate constant pool 2 (1/year): ', bgc%ratecp(2)
           WRITE(logn,'(4X, A50, F12.4)')                                     &
                 'Plant carbon rate constant pool 3 (1/year): ', bgc%ratecp(3)
         ENDIF
         WRITE(logn, *) '------------------------------------------------'//  &
                        '---------'
         WRITE(logn, *) ' Soil parameters: '
         WRITE(logn, patchfmti)        'Soil type for each active (>0%) '//   &
               'patch: ', soil%isoilm(landpt(e)%cstart:landpt(e)%cend)
         WRITE(logn, patchfmtr) 'Fraction of soil which is sand (-): ',       &
               soil%sand(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmtr) 'Fraction of soil which is silt (-): ',       &
               soil%silt(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmtr) 'Fraction of soil which is clay (-): ',       &
               soil%clay(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn, patchfmtr)                                               &
               'Volumetric soil moisture at saturation (m^3/m^3): ',          &
               soil%ssat(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1))
         WRITE(logn,patchfmtr)                                                &
               'Vol. soil moisture at field capacity (m^3/m^3): ',            &
               soil%sfc(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap    &
               - 1))
         WRITE(logn, patchfmtr) 'Vol. soil moisture at wilting point '//      &
               '(m^3/m^3): ', soil%swilt(landpt(e)%cstart:(landpt(e)%cstart   &
               + landpt(e)%nap - 1))
         WRITE(logn,patchfmtr) 'Soil respiration coeff @ 20C (mol/m^2/s): ',  &
               veg%rs20(landpt(e)%cstart:(landpt(e)%cstart+landpt(e)%nap-1))
!              soil%rs20(landpt(e)%cstart:(landpt(e)%cstart+landpt(e)%nap-1))
         WRITE(logn, patchfmtr) 'Suction at saturation (m): ',                &
               soil%sucs(landpt(e)%cstart:(landpt(e)%cstart+landpt(e)%nap-1))
         WRITE(logn, patchfmtr) 'Soil density (kg/m^3): ',                    &
               soil%rhosoil(landpt(e)%cstart:(landpt(e)%cstart +              &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'Soil specific heat capacity (kJ/kg/K): ',    &
               soil%css(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap    &
               - 1))
         WRITE(logn, patchfmtr) 'Parameter b in Campbell equation: ',         &
               soil%bch(landpt(e)%cstart:(landpt(e)%cstart+landpt(e)%nap      &
               - 1))
         WRITE(logn, patchfmte2) 'Hydraulic conductivity @ saturation '//     &
               '(m/s): ', soil%hyds(landpt(e)%cstart:(landpt(e)%cstart +      &
               landpt(e)%nap - 1))
         IF (icycle == 0) THEN
           WRITE(logn,'(4X, A50, F12.4)')                                     &
              'Soil carbon rate constant pool 1 (1/year): ', bgc%ratecs(1)
            WRITE(logn,'(4X, A50, F12.4)')                                    &
              'Soil carbon rate constant pool 2 (1/year): ', bgc%ratecs(2)
         ENDIF

         WRITE(logn, patchfmtr) 'Bare soil albedo, vis (-): ',                &
               soil%albsoil(landpt(e)%cstart:(landpt(e)%cstart +              &
               landpt(e)%nap - 1), 1)
         WRITE(logn, patchfmtr) 'Bare soil albedo, nir (-): ',                &
               soil%albsoil(landpt(e)%cstart:(landpt(e)%cstart +              &
               landpt(e)%nap - 1), 2)
         WRITE(logn, *) '------------------------------------------------'//  &
                        '---------'
         WRITE(logn,'(A35, I8, 1X, A2)') ' CABLE initialisations (land '//    &
               'point ', e, '):'
         WRITE(logn, *) '------------------------------------------------'//  &
                        '---------'
         WRITE(logn, *) ' Soil-specific initialisations, per patch: -----'//  &
                        '---------'
         WRITE(logn, patchfmtr) 'Soil moisture, layer 1: ',                   &
               ssnow%wb(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap    &
               - 1), 1)
         WRITE(logn, patchfmtr) 'Soil moisture, layer 2: ',                   &
               ssnow%wb(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap    &
               - 1), 2)
         WRITE(logn, patchfmtr) 'Soil moisture, layer 3: ',                   &
               ssnow%wb(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap    &
               - 1), 3)
         WRITE(logn, patchfmtr) 'Soil moisture, layer 4: ',                   &
               ssnow%wb(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap    &
               - 1), 4)
         WRITE(logn, patchfmtr) 'Soil moisture, layer 5: ',                   &
               ssnow%wb(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap    &
               - 1), 5)
         WRITE(logn, patchfmtr) 'Soil moisture, layer 6: ',                   &
               ssnow%wb(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap    &
               - 1), 6)
         DO f = landpt(e)%cstart, (landpt(e)%cstart + landpt(e)%nap - 1)
            IF(ANY(ssnow%wb(f, :) < soil%swilt(f)))                           &
                 WRITE(logn, '(3X, A6, I2, A47)')                             &
                 'PATCH ',f - landpt(e)%cstart + 1,                           &
                 ' SOIL MOISTURE INITIALISED BELOW WILTING POINT!'
            IF(ANY(ssnow%wb(f,:)>soil%ssat(f)))                               &
                 WRITE(logn,'(3X, A6, I2, A50)')                              &
                 'PATCH ',f - landpt(e)%cstart + 1,                           &
                 ' SOIL MOISTURE INITIALISED ABOVE SATURATION VALUE!'
         END DO
         WRITE(logn, patchfmtr) 'Soil temperature, layer 1: ',                &
               ssnow%tgg(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 1)
         WRITE(logn, patchfmtr) 'Soil temperature, layer 2: ',                &
               ssnow%tgg(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 2)
         WRITE(logn, patchfmtr) 'Soil temperature, layer 3: ',                &
               ssnow%tgg(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 3)
         WRITE(logn, patchfmtr) 'Soil temperature, layer 4: ',                &
               ssnow%tgg(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 4)
         WRITE(logn, patchfmtr) 'Soil temperature, layer 5: ',                &
               ssnow%tgg(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 5)
         WRITE(logn, patchfmtr) 'Soil temperature, layer 6: ',                &
               ssnow%tgg(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 6)
         IF (icycle == 0) THEN
            WRITE(logn, patchfmtr) 'Soil carbon pool size (g C/m2), pool 1: ',&
               bgc%csoil(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 1)
            WRITE(logn, patchfmtr) 'Soil carbon pool size (g C/m2), pool 2: ',&
               bgc%csoil(landpt(e)%cstart:(landpt(e)%cstart + landpt(e)%nap   &
               - 1), 2)
         ENDIF
         WRITE(logn, patchfmtr) 'Volumetric soil ice, layer 1: ',             &
               ssnow%wbice(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 1)
         WRITE(logn, patchfmtr) 'Volumetric soil ice, layer 2: ',             &
               ssnow%wbice(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 2)
         WRITE(logn, patchfmtr) 'Volumetric soil ice, layer 3: ',             &
               ssnow%wbice(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 3)
         WRITE(logn, patchfmtr) 'Volumetric soil ice, layer 4: ',             &
               ssnow%wbice(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 4)
         WRITE(logn, patchfmtr) 'Volumetric soil ice, layer 5: ',             &
               ssnow%wbice(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 5)
         WRITE(logn, patchfmtr) 'Volumetric soil ice, layer 6: ',             &
               ssnow%wbice(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 6)
         WRITE(logn, patchfmtr) 'Turbulent resistance for soil: ',            &
               ssnow%rtsoil(landpt(e)%cstart:(landpt(e)%cstart +              &
               landpt(e)%nap - 1))
         WRITE(logn, *) ' Snow-specific initialisations, per patch: '//       &
                        '--------------'
         WRITE(logn, patchfmtr) 'Snow liquid water equivalent depth (mm): ',  &
               ssnow%snowd(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr)                                               &
               'Snow liq. water equiv. depth previous tstep (mm): ',          &
               ssnow%osnowd(landpt(e)%cstart:(landpt(e)%cstart +              &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'Overall snow density (kg/m^3): ',            &
               ssnow%ssdnn(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'Snow age (-): ',                             &
               ssnow%snage(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'Snow temperature (K), layer 1: ',            &
               ssnow%tggsn(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 1)
         WRITE(logn, patchfmtr) 'Snow temperature (K), layer 2: ',            &
               ssnow%tggsn(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 2)
         WRITE(logn, patchfmtr) 'Snow temperature (K), layer 3: ',            &
               ssnow%tggsn(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 3)
         WRITE(logn, patchfmtr) 'Snow density (kg/m^3), layer 1: ',           &
               ssnow%ssdn(landpt(e)%cstart:(landpt(e)%cstart +                &
               landpt(e)%nap - 1), 1)
         WRITE(logn, patchfmtr) 'Snow density (kg/m^3), layer 2: ',           &
               ssnow%ssdn(landpt(e)%cstart:(landpt(e)%cstart +                &
               landpt(e)%nap - 1), 2)
         WRITE(logn, patchfmtr) 'Snow density (kg/m^3), layer 3: ',           &
               ssnow%ssdn(landpt(e)%cstart:(landpt(e)%cstart +                &
               landpt(e)%nap - 1), 3)
         WRITE(logn, patchfmtr) 'Snow mass (kg/m^2), layer 1: ',              &
               ssnow%smass(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 1)
         WRITE(logn, patchfmtr) 'Snow mass (kg/m^2), layer 2: ',              &
               ssnow%smass(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 2)
         WRITE(logn, patchfmtr) 'Snow mass (kg/m^2), layer 3: ',              &
               ssnow%smass(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1), 3)
         WRITE(logn, patchfmti) 'Snow layer scheme flag: ',                   &
               ssnow%isflag(landpt(e)%cstart:(landpt(e)%cstart +              &
               landpt(e)%nap - 1))
         WRITE(logn, *) ' Vegetation-specific initialisations, per patch:'//  &
                        ' --------'
         WRITE(logn, patchfmtr) 'Canopy surface water storage (mm): ',        &
               canopy%cansto(landpt(e)%cstart:(landpt(e)%cstart +             &
               landpt(e)%nap - 1))
         WRITE(logn,*) '                 default monthly Leaf area index: '
         DO f = 1, 12
           WRITE(logn,patchfmtr) ' ', defaultLAI(e,f)
         ENDDO
         IF (exists%LAI) THEN
           WRITE(logn,*) 'Check LAI in output for values provided in met file.'
         ELSE
           WRITE(logn,*) 'These default values are used as no LAI in met file.'
         ENDIF
         IF (icycle == 0) THEN
            WRITE(logn, patchfmtr)                                            &
               'Plant carbon pool size (g C/m2), pool 1: ',                   &
                bgc%cplant(landpt(e)%cstart:(landpt(e)%cstart +               &
                landpt(e)%nap - 1), 1)
            WRITE(logn, patchfmtr)                                            &
               'Plant carbon pool size (g C/m2), pool 2: ',                   &
               bgc%cplant(landpt(e)%cstart:(landpt(e)%cstart +                &
               landpt(e)%nap - 1), 2)
         ENDIF
         WRITE(logn, *) ' Other initialisations, per patch: '//               &
                        '----------------------'
         WRITE(logn, patchfmtr) 'Soil+snow albedo (-), visible: ',            &
               ssnow%albsoilsn(landpt(e)%cstart:(landpt(e)%cstart +           &
               landpt(e)%nap - 1), 1)
         WRITE(logn, patchfmtr) 'Soil+snow albedo (-), near infrared: ',      &
               ssnow%albsoilsn(landpt(e)%cstart:(landpt(e)%cstart +           &
               landpt(e)%nap - 1), 2)
         WRITE(logn, patchfmtr) 'Soil+snow albedo (-), thermal: ',            &
               ssnow%albsoilsn(landpt(e)%cstart:(landpt(e)%cstart +           &
               landpt(e)%nap - 1), 3)
         WRITE(logn, patchfmtr) 'Runoff total (mm/time step): ',              &
               ssnow%runoff(landpt(e)%cstart:(landpt(e)%cstart +              &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'Surface runoff (mm/time step): ',            &
               ssnow%rnof1(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1))
         WRITE(logn, patchfmtr) 'Deep drainage runoff (mm/time step): ',      &
               ssnow%rnof2(landpt(e)%cstart:(landpt(e)%cstart +               &
               landpt(e)%nap - 1))
         WRITE(logn, *) '================================================'//  &
                        '========='
         WRITE(logn, *) '================================================'//  &
                        '========='
         WRITE(logn, *)

!jhan:reviewformatting.spacing of whole module
         IF (icycle >= 1) THEN

            WRITE(logn,*) 'CASA-CNP initialisations, per patch:'

            WRITE(logn,patchfmti)                                              &
               '  veg class (0=noveg,1=grassy,2=shrub,3=woody): ',             &
               casamet%iveg2( landpt(e)%cstart :                               &
               (landpt(e)%cstart + landpt(e)%nap-1) )

            WRITE(logn,patchfmti)                                              &
               '                                    soil order: ',             &
               casamet%isorder( landpt(e)%cstart :                             &
               (landpt(e)%cstart + landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               ' patch area (10^9 m^2): ',                                     &
               casamet%areacell(landpt(e)%cstart :                             &
               (landpt(e)%cstart + landpt(e)%nap-1) )*1.0e-9

            WRITE(logn,patchfmtr)                                              &
               '          Nitrogen deposition   (g N/m^2/year): ',             &
               casaflux%Nmindep(landpt(e)%cstart :                             &
               (landpt(e)%cstart + landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               '          Nitrogen fixation     (g N/m^2/year): ',             &
               casaflux%Nminfix(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               '          Phosphorus weathering (g P/m^2/year): ',             &
               casaflux%Pwea(landpt(e)%cstart :                                &
               (landpt(e)%cstart+landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               '          Phosphorus in dust    (g P/m^2/year): ',             &
               casaflux%Pdep(landpt(e)%cstart :                                &
               (landpt(e)%cstart+landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               '  Leaf area index in CASA-CNP: ',                              &
               casamet%glai(landpt(e)%cstart :                                 &
               (landpt(e)%cstart+landpt(e)%nap-1) )

            WRITE(logn,patchfmti)                                              &
               '  Phenological phase: ',                                       &
               phen%phase(landpt(e)%cstart :                                   &
               (landpt(e)%cstart+landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               '  Carbon pools (g C/m^2)       - labile: ',                    &
               casapool%clabile(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               '               - plant - Leaf: ',                              &
               casapool%cplant(landpt(e)%cstart :                              &
               (landpt(e)%cstart+landpt(e)%nap-1), 1 )

            WRITE(logn,patchfmtr)                                              &
               '               - plant - Wood: ',                              &
               casapool%cplant(landpt(e)%cstart :                              &
               (landpt(e)%cstart+landpt(e)%nap-1), 2)

            WRITE(logn,patchfmtr)                                              &
               '               - plant - Root: ',                              &
               casapool%cplant(landpt(e)%cstart :                              &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

            WRITE(logn,patchfmtr)                                              &
               '               - litter - MTB: ',                              &
               casapool%clitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 1 )

            WRITE(logn,patchfmtr)                                              &
               '               - litter - STR: ',                              &
               casapool%clitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 2 )

            WRITE(logn,patchfmtr)                                              &
               '               - litter - CWD: ',                              &
               casapool%clitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

            WRITE(logn,patchfmtr)                                              &
               '               - soil - micro: ',                              &
               casapool%csoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 1 )

            WRITE(logn,patchfmtr)                                              &
               '               - soil -  slow: ',                              &
               casapool%csoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 2 )

            WRITE(logn,patchfmtr)                                              &
               '               - soil - passv: ',                              &
               casapool%csoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

         ENDIF

         IF (icycle >= 2) THEN

            WRITE(logn,patchfmtr) '  Nitrogen pools (g N/m^2) - plant - Leaf: ',&
               casapool%nplant(landpt(e)%cstart :                              &
               (landpt(e)%cstart+landpt(e)%nap-1), 1  )

            WRITE(logn,patchfmtr) '               - plant - Wood: ',           &
               casapool%nplant(landpt(e)%cstart:                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 2 )

            WRITE(logn,patchfmtr) '               - plant - Root: ',           &
               casapool%nplant(landpt(e)%cstart :                              &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

            WRITE(logn,patchfmtr) '               - litter - MTB: ',           &
               casapool%nlitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 1 )

            WRITE(logn,patchfmtr) '               - litter - STR: ',           &
               casapool%nlitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 2 )

            WRITE(logn,patchfmtr) '               - litter - CWD: ',           &
               casapool%nlitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

            WRITE(logn,patchfmtr) '               - soil - micro: ',           &
               casapool%nsoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 1 )

            WRITE(logn,patchfmtr) '               - soil -  slow: ',           &
               casapool%nsoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 2 )

            WRITE(logn,patchfmtr) '               - soil - passv: ',           &
               casapool%nsoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

            WRITE(logn,patchfmtr) '  Mineral nitrogen (inorganic): ',          &
               casapool%nsoilmin(landpt(e)%cstart :                            &
               (landpt(e)%cstart+landpt(e)%nap-1) )

         ENDIF

         IF (icycle == 3) THEN

            WRITE(logn,patchfmtr)                                              &
               '  Phosphorus pools (g P/m^2) - plant - Leaf: ',                &
               casapool%pplant(landpt(e)%cstart :                              &
               (landpt(e)%cstart+landpt(e)%nap-1), 1 )

            WRITE(logn,patchfmtr)                                              &
               '                   - plant - Wood: ',                          &
               casapool%pplant(landpt(e)%cstart :                              &
               (landpt(e)%cstart+landpt(e)%nap-1), 2 )

            WRITE(logn,patchfmtr)                                              &
               '                   - plant - Root: ',                          &
               casapool%pplant(landpt(e)%cstart :                              &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

            WRITE(logn,patchfmtr)                                              &
               '                   - litter - MTB: ',                          &
               casapool%plitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 1 )

            WRITE(logn,patchfmtr)                                              &
               '                   - litter - STR: ',                          &
               casapool%plitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 2 )

            WRITE(logn,patchfmtr)                                              &
               '                   - litter - CWD: ',                          &
               casapool%plitter(landpt(e)%cstart :                             &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

            WRITE(logn,patchfmtr)                                              &
               '                   - soil - micro: ',                          &
               casapool%psoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 1 )

            WRITE(logn,patchfmtr)                                              &
               '                   - soil -  slow: ',                          &
               casapool%psoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 2 )

            WRITE(logn,patchfmtr)                                              &
               '                   - soil - passv: ',                          &
               casapool%psoil(landpt(e)%cstart :                               &
               (landpt(e)%cstart+landpt(e)%nap-1), 3 )

            WRITE(logn,patchfmtr)                                              &
               '  Mineral phosphorus -  Labile: ',                             &
               casapool%psoillab(landpt(e)%cstart :                            &
               (landpt(e)%cstart+landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               '  - Adsorbed: ',                                               &
               casapool%psoilsorb(landpt(e)%cstart :                           &
               (landpt(e)%cstart+landpt(e)%nap-1) )

            WRITE(logn,patchfmtr)                                              &
               '  - Occluded: ',                                               &
               casapool%psoilocc(landpt(e)%cstart :                            &
               (landpt(e)%cstart+landpt(e)%nap-1) )

         ENDIF

      WRITE(logn,*) '========================================================='
      WRITE(logn,*) '========================================================='
      WRITE(logn,*)

      END IF ! if verbose
   END DO

   ENDIF  !report_parameters

END SUBROUTINE report_parameters

  !============================================================================
  SUBROUTINE GWspatialParameters(logn,soil,ssnow)
  ! Read in spatially-specific groundwater parameters

    USE netcdf
    use cable_common_module, only : filename,cable_user,gw_params

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logn ! log file unit number
    TYPE (soil_parameter_type), INTENT(INOUT)  :: soil
    TYPE (soil_snow_type),      INTENT(INOUT)  :: ssnow

    ! local variables
    INTEGER :: ncid_elev, ok, ok2, file_status
    INTEGER :: lon_id,lat_id,tile_id,layer_id
    INTEGER :: klev,ii, jj, kk,e,i,j,k, n
    INTEGER :: nlon,nlat,npatch,nhorz
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: inGW4dtmp
    REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: inGW3dtmp
    REAL, ALLOCATABLE, DIMENSION(:,:)     :: inGWtmp
    !MD Aquifer properties

    !open gw_elev file
    IF (cable_user%gw_model) then

       file_status = NF90_OPEN(trim(filename%type),NF90_NOWRITE,ncid_elev) ! MMY gw_elev -> type

       IF( NF90_INQ_DIMID(ncid_elev,'longitude',lon_id)  .eq. nf90_noerr) THEN
          IF (NF90_INQUIRE_DIMENSION(ncid_elev,lon_id,LEN=nlon) .ne. nf90_noerr) nlon = xdimsize
       ELSE
          nlon = -1
       ENDIF

       IF( NF90_INQ_DIMID(ncid_elev,'latitude',lat_id)  .eq. nf90_noerr) THEN
          IF (NF90_INQUIRE_DIMENSION(ncid_elev,lat_id,LEN=nlat) .ne. nf90_noerr) nlat = ydimsize
       ELSE
          nlat = -1
       ENDIF

       !Have not used tiling yet, will be the future :-)
       IF(  NF90_INQ_DIMID(ncid_elev,'num_tiles',tile_id) .eq. nf90_noerr) THEN
         IF ( NF90_INQUIRE_DIMENSION(ncid_elev,tile_id,LEN=npatch) .ne. nf90_noerr) npatch=1
       ELSE
          npatch=-1
       ENDIF

       IF(  NF90_INQ_DIMID(ncid_elev,'nhorz',layer_id) .eq. nf90_noerr) THEN
         IF ( NF90_INQUIRE_DIMENSION(ncid_elev,layer_id,LEN=nhorz) .ne. nf90_noerr) nhorz=1
       ELSE
          nhorz=ms
       ENDIF

       IF ((nlon*nlat*npatch*nhorz .ne. xdimsize*ydimsize*ms*mpatch) .and. &
           xdimsize*ydimsize .ne. 1) THEN
          WRITE(logn,*) 'Errors reading the dimensions from '//filename%gw_elev

          IF (nlon .lt. 0 .or. nlon .ne. xdimsize) &
             WRITE(logn,*) 'nlon: found ',nlon,' need to have ',xdimsize,' to match forcing.'
          IF (nlat .lt. 0 .or. nlat .ne. ydimsize) &
             WRITE(logn,*) 'nlat: found ',nlat,' need to have ',ydimsize,' to match forcing.'
          IF (nhorz .lt. 0 .or. nhorz .ne. ms) &
             WRITE(logn,*) 'nhorz: found ',nhorz,' need to have ',ms,' to match forcing.'
          IF (npatch .lt. 0) &
             WRITE(logn,*) 'npatch: found ',npatch,' need to have ',mpatch,' to match forcing.'

             WRITE(logn,*) 'Setting dims to forcing file values, CHECK THE OUTPUT!'
             !amu561: code should probably stop here instead of re-writing
             !dimensions?
             nlon=xdimsize; nlat=ydimsize; npatch=mpatch; nhorz=ms
       ENDIF

       allocate(inGWtmp(nlon,nlat))
       allocate(inGW3dtmp(nlon,nlat,ms))
       allocate(inGW4dtmp(nlon,nlat,ms,mpatch))

    !1
    soil%elev(:) = get_gw_data(ncid_elev,file_status,'elevation',50.0,nlon,nlat)
    !2
    soil%elev_std(:) = get_gw_data(ncid_elev,file_status,'elevation_std',25.0,nlon,nlat)
    !3
    soil%slope(:) = get_gw_data(ncid_elev,file_status,'slope',0.08,nlon,nlat)
    !4
    soil%slope_std(:) = get_gw_data(ncid_elev,file_status,'slope_std',0.03,nlon,nlat)
    !5
    soil%GWdz(:) = get_gw_data(ncid_elev,file_status,'dtb',25.0,nlon,nlat)
    where(soil%GWdz(:) .lt. 5.0)  soil%GWdz(:) = 5.0
    where(soil%GWdz(:) .gt. 50.0) soil%GWdz(:) = 50.0
    !6
    soil%drain_dens(:) = get_gw_data(ncid_elev,file_status,'drainage_density',0.0008,nlon,nlat)
    where( soil%drain_dens(:) .lt. 1.0e-6) soil%drain_dens(:)=1.0e-6
    where( soil%drain_dens(:) .gt. 0.02  ) soil%drain_dens(:)=0.02
    !7
    soil%GWhyds_vec(:) = get_gw_data(ncid_elev,file_status,'permeability',1.0e-6,nlon,nlat)
    soil%GWhyds_vec(:) = 1000._r_2 * soil%GWhyds_vec(:)
    !8
    soil%GWssat_vec(:) = get_gw_data(ncid_elev,file_status,'Sy',inssat(:,:),nlon,nlat)
    soil%GWssat_vec(:) = max(0.23_r_2,soil%GWssat_vec(:))

    soil%GWwatr(:) = 0.0

    soil%ssat_vec(:,:) = get_gw_data(ncid_elev,file_status,'ssat_vec',inssat(:,:),nlon,nlat,ms)
    inGWtmp(:,:) = 0.66*inssat(:,:)

    soil%sfc_vec(:,:) = get_gw_data(ncid_elev,file_status,'sfc_vec',inGWtmp(:,:),nlon,nlat,ms)
    inGWtmp(:,:) = 0.15*inssat(:,:)

    soil%swilt_vec(:,:) = get_gw_data(ncid_elev,file_status,'swilt_vec',inGWtmp(:,:),nlon,nlat,ms)

    inGWtmp(:,:) = 0.01*inssat(:,:)
    soil%watr(:,:) = get_gw_data(ncid_elev,file_status,'watr',inGWtmp(:,:),nlon,nlat,ms)

    inGWtmp(:,:) = 1000.0*inhyds(:,:)
    soil%hyds_vec(:,:) = get_gw_data(ncid_elev,file_status,'hyds_vec',inGWtmp(:,:),nlon,nlat,ms)

    inGWtmp(:,:) = abs(insucs(:,:))
    soil%sucs_vec(:,:) = get_gw_data(ncid_elev,file_status,'sucs_vec',inGWtmp(:,:),nlon,nlat,ms)
    soil%sucs_vec(:,:) =  1000._r_2*abs(soil%sucs_vec(:,:)  )

    !add last laery to aquifer
    !should have zero head at top of aquifer, however this can be well below
    !the soil column, so reading below we can treat top part of whenm dry
    !as unsat flow
    soil%GWsucs_vec(:) = soil%sucs_vec(:,ms)

    soil%bch_vec(:,:) = get_gw_data(ncid_elev,file_status,'bch_vec',inbch(:,:),nlon,nlat,ms)
    soil%GWbch_vec(:) = soil%bch_vec(:,ms)

    soil%rhosoil_vec(:,:) = get_gw_data(ncid_elev,file_status,'rhosoil_vec',inrhosoil(:,:),nlon,nlat,ms)

    soil%css_vec(:,:) = get_gw_data(ncid_elev,file_status,'css_vec',incss(:,:),nlon,nlat,ms)

    soil%cnsd_vec(:,:) = get_gw_data(ncid_elev,file_status,'cnsd_vec',incnsd(:,:),nlon,nlat,ms)

    if (file_status .eq. nf90_noerr) &
         file_status = nf90_close(ncid_elev)

    !set the default IC for hysteresis state
    ssnow%smp_hys(:,:) = -soil%sucs_vec(:,:)
    ssnow%hys_fac(:,:) = 1.0
    ssnow%watr_hys(:,:) = soil%watr(:,:)
    ssnow%ssat_hys(:,:) = soil%ssat_vec(:,:)

    ELSE  !gw_model=false

    DO e=1,mland
        soil%elev(landpt(e)%cstart:landpt(e)%cend) = 50.0
       soil%elev_std(landpt(e)%cstart:landpt(e)%cend) =25
       soil%slope(landpt(e)%cstart:landpt(e)%cend) =0.08
       soil%slope_std(landpt(e)%cstart:landpt(e)%cend) =0.03
       soil%GWdz(landpt(e)%cstart:landpt(e)%cend) =25.0
       soil%drain_dens(landpt(e)%cstart:landpt(e)%cend) =0.008
       soil%GWhyds_vec(landpt(e)%cstart:landpt(e)%cend) =1.0e-6
       soil%GWssat_vec(landpt(e)%cstart:landpt(e)%cend) =0.23
       soil%GWwatr(landpt(e)%cstart:landpt(e)%cend) =0.0
       DO klev=1,ms
          !layered_in_soils(landpt(e)%cstart:landpt(e)%cend,klev,1) =&
          soil%ssat_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             inssat(landpt(e)%ilon,landpt(e)%ilat)


          soil%sfc_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             insfc(landpt(e)%ilon,landpt(e)%ilat)


          soil%swilt_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             inswilt(landpt(e)%ilon,landpt(e)%ilat)


          soil%hyds_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             inhyds(landpt(e)%ilon,landpt(e)%ilat)


          soil%sucs_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             insucs(landpt(e)%ilon,landpt(e)%ilat)


          soil%bch_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             inbch(landpt(e)%ilon,landpt(e)%ilat)


          soil%rhosoil_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             inrhosoil(landpt(e)%ilon,landpt(e)%ilat)


          soil%cnsd_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             incnsd(landpt(e)%ilon,landpt(e)%ilat)


          soil%css_vec(landpt(e)%cstart:landpt(e)%cend,klev) =& !1
             incss(landpt(e)%ilon,landpt(e)%ilat)

          end do
       end do
   end if

    DO e=1,mland
      ii=landpt(e)%cstart
      jj=landpt(e)%cend
      DO klev=1,ms
         soil%sand_vec(ii:jj,klev) = insand(landpt(e)%ilon,landpt(e)%ilat)
         soil%clay_vec(ii:jj,klev) = inclay(landpt(e)%ilon,landpt(e)%ilat)
         soil%silt_vec(ii:jj,klev) = insilt(landpt(e)%ilon,landpt(e)%ilat)
         soil%org_vec(ii:jj,klev)  = inorg(landpt(e)%ilon,landpt(e)%ilat)
         !
         IF (klev .gt. 3) soil%org_vec(ii,klev) =  0.0
      END DO
    ENDDO

    !removed gw_soils for now
  
    !___________________ MMY GW_vec = ms6 ____________________
     soil%GWhyds_vec(:) = soil%hyds_vec(:,ms) 
     soil%GWssat_vec(:) = soil%ssat_vec(:,ms) 
     soil%GWsucs_vec(:) = soil%sucs_vec(:,ms)
     soil%GWbch_vec(:)  = soil%bch_vec(:,ms) 
     soil%GWwatr(:)     = soil%watr(:,ms)    
    !_________________________________________________________


  END SUBROUTINE GWspatialParameters

   function get_gw_2d_var_constdef(ncfile_id,try_it,varname,default_const,nlon,nlat) result(data_vec)
         use netcdf
         use cable_def_types_mod, only : r_2,mp
         real, intent(in) :: default_const
         integer, intent(in) :: try_it,nlon,nlat
         character(len=*), intent(in) :: varname
         integer, intent(in)               :: ncfile_id
         real(r_2), dimension(mp) :: data_vec

         real, dimension(nlon,nlat) :: GW2d_data
         real, dimension(nlon,nlat) :: default_data
         integer :: i,j,k,n
         integer :: varid
         integer :: varinq_status,varget_status

         default_data(:,:) = default_const

         if (.not.cable_user%gw_model) then

            GW2d_data(:,:) = default_data(:,:)

         else

         if (try_it .ne. nf90_noerr) then
            varinq_status = try_it
            varget_status = try_it
         else

            varinq_status = nf90_inq_varid(ncfile_id,trim(varname),varid)
            varget_status = varinq_status

            if (varinq_status .eq. nf90_noerr) then
               varget_status = nf90_get_var(ncfile_id,varid,GW2d_data)
            endif

         end if

         if (try_it.ne. nf90_noerr .or. varget_status .ne. nf90_noerr) then
            GW2d_data(:,:) = default_data(:,:)
         end if
         end if

         do i=1,mland
            data_vec(landpt(i)%cstart:landpt(i)%cend) = &
                   real(GW2d_data(landpt(i)%ilon,landpt(i)%ilat),r_2)
         end do

         return

   end function


   function get_gw_2d_var(ncfile_id,try_it,varname,default_data,nlon,nlat) result(data_vec)
         use netcdf
         use cable_def_types_mod, only : r_2,mp
         real, dimension(:,:), intent(in) :: default_data
         integer, intent(in) :: try_it,nlon,nlat
         character(len=*), intent(in) :: varname
         integer, intent(in)               :: ncfile_id
         real(r_2), dimension(mp) :: data_vec

         real, dimension(nlon,nlat) :: GW2d_data
         integer :: i,j,k,n
         integer :: varid
         integer :: varinq_status,varget_status

         if (.not.cable_user%gw_model) then

            GW2d_data(:,:) = default_data(:,:)

         else

         if (try_it .ne. nf90_noerr) then
            varinq_status = try_it
            varget_status = try_it
         else

            varinq_status = nf90_inq_varid(ncfile_id,trim(varname),varid)
            varget_status = varinq_status

            if (varinq_status .eq. nf90_noerr) then
               varget_status = nf90_get_var(ncfile_id,varid,GW2d_data)
            endif

         end if

         if (try_it.ne. nf90_noerr .or. varget_status .ne. nf90_noerr) then
            GW2d_data(:,:) = default_data(:,:)
         end if
         end if

         do i=1,mland
            data_vec(landpt(i)%cstart:landpt(i)%cend) = &
                   real(GW2d_data(landpt(i)%ilon,landpt(i)%ilat),r_2)
         end do

         return

   end function



   function get_gw_3d_var_constdef(ncfile_id,try_it,varname,default_const,nlon,nlat,ms) result(data_vec)
         use netcdf
         use cable_def_types_mod, only : r_2,mp
         real, intent(in) :: default_const
         integer, intent(in) :: try_it,nlon,nlat,ms
         character(len=*), intent(in) :: varname
         integer, intent(in)               :: ncfile_id
         real(r_2), dimension(mp,ms) :: data_vec

         real, dimension(nlon,nlat,ms) :: GW3d_data
         integer :: i,j,k,n
         integer :: varid
         integer :: varinq_status,varget_status
         real, dimension(nlon,nlat) :: default_data

          default_data(:,:) = default_const

         if (.not. cable_user%gw_model) then
            do k=1,ms
               GW3d_data(:,:,k) = default_data(:,:)
            end do

         else
         if (try_it .ne. nf90_noerr) then
            varinq_status = try_it
            varget_status = try_it
         else

            varinq_status = nf90_inq_varid(ncfile_id,trim(varname),varid)
            varget_status = varinq_status

            if (varinq_status .eq. nf90_noerr) then
               varget_status = nf90_get_var(ncfile_id,varid,GW3d_data)
            endif

         end if

         if (try_it.ne. nf90_noerr .or. varget_status .ne. nf90_noerr) then
            do k=1,ms
               GW3d_data(:,:,k) = default_data(:,:)
            end do
         end if


         end if

         do i=1,mland
            do k=1,ms
               data_vec(landpt(i)%cstart:landpt(i)%cend,k) = &
                      real(GW3d_data(landpt(i)%ilon,landpt(i)%ilat,k),r_2)
            end do
         end do

         return

   end function



   function get_gw_3d_var(ncfile_id,try_it,varname,default_data,nlon,nlat,ms) result(data_vec)
         use netcdf
         use cable_def_types_mod, only : r_2,mp
         real, dimension(:,:), intent(in) :: default_data
         integer, intent(in) :: try_it,nlon,nlat,ms
         character(len=*), intent(in) :: varname
         integer, intent(in)               :: ncfile_id
         real(r_2), dimension(mp,ms) :: data_vec

         real, dimension(nlon,nlat,ms) :: GW3d_data
         integer :: i,j,k,n
         integer :: varid
         integer :: varinq_status,varget_status

         if (.not. cable_user%gw_model) then
            do k=1,ms
               GW3d_data(:,:,k) = default_data(:,:)
            end do

         else
         if (try_it .ne. nf90_noerr) then
            varinq_status = try_it
            varget_status = try_it
         else

            varinq_status = nf90_inq_varid(ncfile_id,trim(varname),varid)
            varget_status = varinq_status

            if (varinq_status .eq. nf90_noerr) then
               varget_status = nf90_get_var(ncfile_id,varid,GW3d_data)
            endif

         end if

         if (try_it.ne. nf90_noerr .or. varget_status .ne. nf90_noerr) then
            do k=1,ms
               GW3d_data(:,:,k) = default_data(:,:)
            end do
         end if


         end if

         do i=1,mland
            do k=1,ms
               data_vec(landpt(i)%cstart:landpt(i)%cend,k) = &
                      real(GW3d_data(landpt(i)%ilon,landpt(i)%ilat,k),r_2)
            end do
         end do

         return

   end function



   function get_gw_4d_var_constdef(ncfile_id,try_it,varname,default_const,nlon,nlat,ms,npatch) result(data_vec)
         use netcdf
         use cable_def_types_mod, only : r_2,mp
         real, intent(in) :: default_const
         integer, intent(in) :: try_it,nlon,nlat,ms,npatch
         character(len=*),  intent(in) :: varname
         integer, intent(in)               :: ncfile_id
         real(r_2), dimension(mp,ms,mpatch) :: data_vec

         integer :: i,j,k,n
         integer :: varid
         integer :: varinq_status,varget_status
         real, dimension(nlon,nlat,ms,mpatch) :: GW4d_data
         real, dimension(nlon,nlat,ms)        :: GW3d_data
         real, dimension(nlon,nlat) :: default_data

         default_data(:,:) = default_const

         if (.not.cable_user%gw_model) then
            do k=1,ms
               do n=1,npatch
                  GW4d_data(:,:,k,n) = default_data(:,:)
               end do
            end do
         else
         if (try_it .ne. nf90_noerr) then
            varinq_status = try_it
            varget_status = try_it
         else

            varinq_status = nf90_inq_varid(ncfile_id,trim(varname),varid)
            varget_status = varinq_status

            if (varinq_status .eq. nf90_noerr) then
               if (mpatch .gt. 1) then
                  varget_status = nf90_get_var(ncfile_id,varid,GW4d_data)
               else
                  varget_status = nf90_get_var(ncfile_id,varid,GW3d_data)
                  GW4d_data(:,:,:,1) = GW3d_data(:,:,:)
               end if
            endif

         end if

         if (try_it.ne. nf90_noerr .or. varget_status .ne. nf90_noerr) then
            do k=1,ms
               do n=1,mpatch
                  GW4d_data(:,:,k,n) = default_data(:,:)
               end do
            end do
         end if


         end if

         do i=1,mland
            do k=1,ms
               do j=1,mpatch
                  data_vec(landpt(i)%cstart:landpt(i)%cend,k,j) = &
                         real(GW4d_data(landpt(i)%ilon,landpt(i)%ilat,k,j),r_2)
               end do
            end do
         end do

         return

   end function




   function get_gw_4d_var(ncfile_id,try_it,varname,default_data,nlon,nlat,ms,npatch) result(data_vec)
         use netcdf
         use cable_def_types_mod, only : r_2,mp
         real, dimension(:,:), intent(in) :: default_data
         integer, intent(in) :: try_it,nlon,nlat,ms,npatch
         character(len=*),  intent(in) :: varname
         integer, intent(in)               :: ncfile_id
         real(r_2), dimension(mp,ms,mpatch) :: data_vec

         integer :: i,j,k,n
         integer :: varid
         integer :: varinq_status,varget_status
         real, dimension(nlon,nlat,ms,mpatch) :: GW4d_data
         real, dimension(nlon,nlat,ms)        :: GW3d_data


         if (.not.cable_user%gw_model) then
            do k=1,ms
               do n=1,npatch
                  GW4d_data(:,:,k,n) = default_data(:,:)
               end do
            end do
         else
         if (try_it .ne. nf90_noerr) then
            varinq_status = try_it
            varget_status = try_it
         else

            varinq_status = nf90_inq_varid(ncfile_id,trim(varname),varid)
            varget_status = varinq_status

            if (varinq_status .eq. nf90_noerr) then
               if (mpatch .gt. 1) then
                  varget_status = nf90_get_var(ncfile_id,varid,GW4d_data)
               else
                  varget_status = nf90_get_var(ncfile_id,varid,GW3d_data)
                  GW4d_data(:,:,:,1) = GW3d_data(:,:,:)
               end if

            endif

         end if

         if (try_it.ne. nf90_noerr .or. varget_status .ne. nf90_noerr) then
            do k=1,ms
               do n=1,mpatch
                  GW4d_data(:,:,k,n) = default_data(:,:)
               end do
            end do
         end if


         end if

         do i=1,mland
            do k=1,ms
               do j=1,mpatch
                  data_vec(landpt(i)%cstart:landpt(i)%cend,k,j) = &
                         real(GW4d_data(landpt(i)%ilon,landpt(i)%ilat,k,j),r_2)
               end do
            end do
         end do

         return

   end function


END MODULE cable_param_module
