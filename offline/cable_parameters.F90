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
!          input files
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
  IMPLICIT NONE
  PRIVATE
  PUBLIC get_default_params, write_default_params, derived_parameters,         &
         check_parameter_values, report_parameters, parID_type,                &
         write_cnp_params
  INTEGER :: patches_in_parfile=4 ! # patches in default global parameter
                                       ! file

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

  ! vars intro for Ticket #27 
  INTEGER, DIMENSION(:, :),     ALLOCATABLE :: inSoilColor

CONTAINS

  SUBROUTINE get_default_params(logn, vegparmnew)
    use cable_common_module, only : get_type_parameters, filename,             &
                                    calcsoilalbedo
  ! Load parameters for each veg type and each soil type. (get_type_parameters)
  ! Also read in initial information for each grid point. (read_gridinfo)
  ! Count to obtain 'landpt', 'max_vegpatches' and 'mp'. (countPatch)
  !
  ! New input structure using netcdf and introduced 'month' to initialize
  ! soil profiles with the correct monthly average values (BP apr2010)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logn     ! log file unit number
    LOGICAL,      INTENT(IN) :: vegparmnew ! new format input file (BP dec2007)

    ! local variables
    INTEGER :: npatch
    INTEGER :: nlon
    INTEGER :: nlat

    ! Get parameter values for all default veg and soil types:
    CALL get_type_parameters(logn, vegparmnew, classification)

    WRITE(logn,*) ' Reading grid info from ', TRIM(filename%type)
    WRITE(logn,*) ' And assigning C4 fraction according to veg classification.'
    WRITE(logn,*) 
    CALL read_gridinfo(nlon,nlat,npatch)

    IF (soilparmnew) THEN
      PRINT *,      'Use spatially-specific soil properties; ', nlon, nlat
      WRITE(logn,*) 'Use spatially-specific soil properties; ', nlon, nlat
      CALL spatialSoil(nlon, nlat, logn)
    ENDIF

    ! include prescribed soil colour in determining albedo - Ticket #27
    IF (calcsoilalbedo) THEN
       CALL read_soilcolor(logn)
    END IF

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
    INTEGER :: ii, jj, kk
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

    ! check dimensions of soil-layers and time
    IF (nslayer /= ms .OR. ntime /= 12) THEN
      PRINT *, 'Variable dimensions do not match:'
      PRINT *, 'nslayer and ms = ', nslayer, ms
      PRINT *, 'ntime not equal 12 months: ', ntime
      CALL abort('Variable dimensions do not match (read_gridinfo)')
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

    ok = NF90_INQ_VARID(ncid, 'latitude', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable latitude.')
    ok = NF90_GET_VAR(ncid, varID, inLat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable latitude.')

    ok = NF90_INQ_VARID(ncid, 'iveg', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error finding variable iveg.')
    ok = NF90_GET_VAR(ncid, varID, idummy)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error reading variable iveg.')
    inVeg(:, :, 1) = idummy(:,:) ! npatch=1 in 1x1 degree input

    ok = NF90_INQ_VARID(ncid, 'patchfrac', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error finding variable patchfrac.')
    ok = NF90_GET_VAR(ncid, varID, inPFrac)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error reading variable patchfrac.')
    inPFrac(:, :, 1) = rdummy(:, :)

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
    DO kk = 1, nband
      inALB(:,:,1,kk) = r3dum(:,:,kk)
    ENDDO

    ok = NF90_INQ_VARID(ncid, 'SnowDepth', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error finding variable SnowDepth.')
    ok = NF90_GET_VAR(ncid,varID,r3dum2)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                                        'Error reading variable SnowDepth.')
    DO kk = 1, ntime
      inSND(:, :, 1, kk) = r3dum2(:, :, kk)
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
    INTEGER :: ncid, ok, ii, jj
    INTEGER :: xID, yID, fieldID
    INTEGER :: xlon, xlat
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
    ELSEWHERE (in2alb > 0.14 .and. in2alb <= 0.20)
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
  integer :: logn
  ! Read soil color
  !
  ! Input variables:
  !   filename%soilcolor  - via cable_IO_vars_module
  ! Output variables:
  !   soilcol    - via cable_param_module
  !
  ! New input structure using netcdf

  End SUBROUTINE read_soilcolor
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

    ! range of longitudes from input file (inLon) should be -180 to 180,
    ! and longitude(:) has already been converted to -180 to 180 for CABLE.
    landpt(:)%ilon = -999
    landpt(:)%ilat = -999
    ncount = 0
    DO kk = 1, mland
      distance = 3.0 ! initialise, units are degrees
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
    IF (ncount > mland * nmetpatches) THEN
      PRINT *, ncount, ' should not be greater than mland*nmetpatches.'
      PRINT *, 'mland, nmetpatches = ', mland, nmetpatches
      STOP
    END IF
    DEALLOCATE(inLon, inLat)

    ! Set the maximum number of active patches to that read from met file:
    max_vegpatches = MAXVAL(landpt(:)%nap)
    IF (max_vegpatches /= nmetpatches) THEN
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
                                  vegparmnew, month, TFRZ)
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

    USE cable_common_module, only : vegin, soilin, calcsoilalbedo
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
    TYPE (canopy_type),         INTENT(OUT)   :: canopy
    TYPE (roughness_type),      INTENT(OUT)   :: rough
    TYPE (radiation_type),      INTENT(OUT)   :: rad

    INTEGER :: e,f,h  ! do loop counter
    INTEGER :: is     ! YP oct07
    INTEGER :: ir     ! BP sep2010
    REAL :: totdepth  ! YP oct07
    REAL :: tmp       ! BP sep2010

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
    ssnow%smass  = 0.0   ! snow mass per layer (kg/m^2)
    ssnow%runoff = 0.0   ! runoff total = subsurface + surface runoff
    ssnow%rnof1  = 0.0   ! surface runoff (mm/timestepsize)
    ssnow%rnof2  = 0.0   ! deep drainage (mm/timestepsize)
    ssnow%rtsoil = 100.0 ! turbulent resistance for soil
    canopy%ga     = 0.0  ! ground heat flux (W/m2)
    canopy%dgdtg  = 0.0  ! derivative of ground heat flux wrt soil temp
    canopy%fev    = 0.0  ! latent heat flux from vegetation (W/m2)
    canopy%fes    = 0.0  ! latent heat flux from soil (W/m2)
    canopy%fhs    = 0.0  ! sensible heat flux from soil (W/m2)

    ! *******************************************************************
    ! parameters that are not spatially dependent
    soil%zse = (/.022, .058, .154, .409, 1.085, 2.872/) ! layer thickness nov03

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
      DO is = 1, ms
        ssnow%tgg(landpt(e)%cstart:landpt(e)%cend, is) =                       &
                                 inTGG(landpt(e)%ilon,landpt(e)%ilat, is, month)
        ssnow%wb(landpt(e)%cstart:landpt(e)%cend, is) =                        &
                                 inWB(landpt(e)%ilon, landpt(e)%ilat, is, month)
      END DO

      ! Set initial snow depth and snow-free soil albedo
      DO is = 1, landpt(e)%cend - landpt(e)%cstart + 1  ! each patch
        DO ir = 1, nrb                                  ! each band
          ssnow%albsoilsn(landpt(e)%cstart + is - 1, ir)                       &
              = inALB(landpt(e)%ilon, landpt(e)%ilat, is, ir) ! various rad band
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
                                     -1.* insucs(landpt(e)%ilon, landpt(e)%ilat)
      soil%rhosoil(landpt(e)%cstart:landpt(e)%cend) =                          &
                                       inrhosoil(landpt(e)%ilon, landpt(e)%ilat)
      soil%css(landpt(e)%cstart:landpt(e)%cend) =                              &
                                           incss(landpt(e)%ilon, landpt(e)%ilat)
      soil%cnsd(landpt(e)%cstart:landpt(e)%cend) =                             &
                                          incnsd(landpt(e)%ilon, landpt(e)%ilat)

      ENDIF


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

       ! Prescribe parameters for current gridcell based on veg/soil type (which
       ! may have loaded from default value file or met file):
       DO h = landpt(e)%cstart, landpt(e)%cend ! over each patch in current grid
          veg%frac4(h)    = vegin%frac4(veg%iveg(h))
          veg%taul(h,1)    = vegin%taul(1,veg%iveg(h))
          veg%taul(h,2)    = vegin%taul(2,veg%iveg(h))
          veg%refl(h,1)    = vegin%refl(1,veg%iveg(h))
          veg%refl(h,2)    = vegin%refl(2,veg%iveg(h))
          veg%canst1(h)   = vegin%canst1(veg%iveg(h))
          veg%dleaf(h)    = vegin%dleaf(veg%iveg(h))
          veg%vcmax(h)    = vegin%vcmax(veg%iveg(h))
          veg%ejmax(h)    = vegin%ejmax(veg%iveg(h))
          veg%hc(h)       = vegin%hc(veg%iveg(h))
          veg%xfang(h)    = vegin%xfang(veg%iveg(h))
          veg%vbeta(h)    = vegin%vbeta(veg%iveg(h))
          veg%xalbnir(h)  = vegin%xalbnir(veg%iveg(h))
          veg%rp20(h)     = vegin%rp20(veg%iveg(h))
          veg%rpcoef(h)   = vegin%rpcoef(veg%iveg(h))
          veg%rs20(h)     = vegin%rs20(veg%iveg(h))
          veg%shelrb(h)   = vegin%shelrb(veg%iveg(h))
          veg%wai(h)      = vegin%wai(veg%iveg(h))
          veg%vegcf(h)    = vegin%vegcf(veg%iveg(h))
          veg%extkn(h)    = vegin%extkn(veg%iveg(h))
          veg%tminvj(h)   = vegin%tminvj(veg%iveg(h))
          veg%tmaxvj(h)   = vegin%tmaxvj(veg%iveg(h))
          bgc%cplant(h,:) = vegin%cplant(:, veg%iveg(h))
          bgc%csoil(h,:)  = vegin%csoil(:, veg%iveg(h))
          bgc%ratecp(:)   = vegin%ratecp(:, veg%iveg(h))
          bgc%ratecs(:)   = vegin%ratecs(:, veg%iveg(h))
          veg%froot(h,:)  = vegin%froot(:, veg%iveg(h))
          soil%silt(h)    =  soilin%silt(soil%isoilm(h))
          soil%clay(h)    =  soilin%clay(soil%isoilm(h))
          soil%sand(h)    =  soilin%sand(soil%isoilm(h))
          IF (.NOT. soilparmnew) THEN   ! Q,Zhang @ 12/20/2010
            soil%swilt(h)   =  soilin%swilt(soil%isoilm(h))
            soil%sfc(h)     =  soilin%sfc(soil%isoilm(h))
            soil%ssat(h)    =  soilin%ssat(soil%isoilm(h))
            soil%bch(h)     =  soilin%bch(soil%isoilm(h))
            soil%hyds(h)    =  soilin%hyds(soil%isoilm(h))
            soil%sucs(h)    =  soilin%sucs(soil%isoilm(h))
            soil%rhosoil(h) =  soilin%rhosoil(soil%isoilm(h))
            soil%css(h)     =  soilin%css(soil%isoilm(h))
          END IF
          rad%latitude(h) = latitude(e)
          veg%ejmax(h) = 2.0 * veg%vcmax(h)
       END DO ! over each veg patch in land point
    END DO ! over all land points
    soil%albsoil = ssnow%albsoilsn

    ! check tgg and alb
    IF(ANY(ssnow%tgg > 350.0) .OR. ANY(ssnow%tgg < 180.0))                     &
           CALL abort('Soil temps nuts')
    IF(ANY(ssnow%albsoilsn > 1.0) .OR. ANY(ssnow%albsoilsn < 0.0))             &
           CALL abort('Albedo nuts')

    WRITE(logn, *)

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
               vegin%tmaxvj, vegin%vbeta, vegin%rootbeta, vegin%froot,         &
               vegin%cplant, vegin%csoil, vegin%ratecp, vegin%ratecs,          &
               vegin%xalbnir, vegin%length, vegin%width                        &
               )
    !         vegf_temp,urbanf_temp,lakef_temp,icef_temp, &

    ! if using old format veg_parm input file, need to define veg%deciduous
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
    WHERE(ssnow%tgg(:, :) < 273.15)
      ssnow%wbice(:,:) = ssnow%wb(:, :) * 0.8
    ELSEWHERE
      ssnow%wbice(:, :) = 0.0
    END WHERE

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

  END SUBROUTINE write_cnp_params
  !============================================================================
  SUBROUTINE derived_parameters(soil, sum_flux, bal, ssnow, veg, rough)
    ! Gives values to parameters that are derived from other parameters.
    TYPE (soil_snow_type),      INTENT(IN)    :: ssnow
    TYPE (veg_parameter_type),  INTENT(IN)    :: veg
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (sum_flux_type),       INTENT(INOUT) :: sum_flux
    TYPE (balances_type),       INTENT(INOUT) :: bal
    TYPE (roughness_type),      INTENT(INOUT) :: rough

    INTEGER :: j ! do loop counter
    REAL(r_2)    :: temp(mp)
    REAL    :: tmp2(mp)

    ! Construct derived parameters and zero initialisations,
    ! regardless of where parameters and other initialisations 
    ! have loaded from:
    soil%zshh(1) = 0.5 * soil%zse(1) ! distance between consecutive layer
                                     ! midpoints:
    soil%zshh(ms + 1) = 0.5 * soil%zse(ms)
    soil%zshh(2:ms)   = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
    IF ( .NOT. soilparmnew) THEN  ! Q,Zhang @ 12/20/2010
      soil%cnsd  = soil%sand * 0.3 + soil%clay * 0.25                          &
                   + soil%silt * 0.265 ! set dry soil thermal conductivity
                                       ! [W/m/K]
    END IF
    soil%hsbh   = soil%hyds*ABS(soil%sucs) * soil%bch ! difsat*etasat
    soil%ibp2   = NINT(soil%bch) + 2
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
!      ssnow%owetfac = ssnow%owetfac &
!                    * (1.0 - REAL(ssnow%wbice(:,1)/ssnow%wb(:,1)))**2
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
    DO j=1, ms
       bal%wbtot0 = bal%wbtot0 + REAL(ssnow%wb(:, j)) * soil%zse(j)       &
                    * 1000.0
    END DO
    bal%osnowd0 = ssnow%osnowd

  END SUBROUTINE derived_parameters
  !============================================================================
  SUBROUTINE check_parameter_values(soil, veg, ssnow)
    ! Checks for basic inconsistencies in parameter values
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameter data
    TYPE (veg_parameter_type),  INTENT(IN)    :: veg   ! vegetation parameter
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
          + landpt(i)%nap - 1)) > mstype)) THEN
          WRITE(*,*) 'SUBROUTINE load_parameters:'
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
          IF((soil%sand(landpt(i)%cstart + j - 1)                              &
              + soil%silt(landpt(i)%cstart + j - 1)                            &
              + soil%clay(landpt(i)%cstart + j - 1)) > 1.0000001 .OR.          &
             (soil%sand(landpt(i)%cstart + j - 1)                              &
              + soil%silt(landpt(i)%cstart + j - 1)                            &
              + soil%clay(landpt(i)%cstart + j - 1)) < 0.9999999) THEN
             WRITE(*,*) 'SUBROUTINE load_parameters:'
             WRITE(*,*) 'At land point number:', i
             WRITE(*,*) '        patch number:', j
             WRITE(*,*) 'Clay fraction is ',soil%clay(landpt(i)%cstart + j - 1)
             WRITE(*,*) 'Sand fraction is ',soil%sand(landpt(i)%cstart + j - 1)
             WRITE(*,*) 'Silt fraction is ',soil%silt(landpt(i)%cstart + j - 1)
             WRITE(*,*) 'SUM:',soil%sand(landpt(i)%cstart + j - 1)             &
                               + soil%silt(landpt(i)%cstart + j - 1)           &
                               + soil%clay(landpt(i)%cstart + j - 1)
             CALL abort ('clay+sand+silt fraction does not sum to 1!')
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
             CALL abort ('Sum of fraction of roots in each soil layer /= 1!')
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
    DO i = 1, ms
       WHERE(ssnow%wb(:, i) > soil%ssat) ! Can only happen due to i/o issues
          ssnow%wb(:, i) = 0.9999 * soil%ssat
       END WHERE    
    END DO

  END SUBROUTINE check_parameter_values
!===============================================================================
SUBROUTINE report_parameters(logn, soil, veg, bgc, rough,                    &
                               ssnow, canopy, casamet, casapool, casaflux,     &
                               phen, vegparmnew, verbose )
   USE cable_common_module, ONLY : veg_desc, soil_desc 
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

END SUBROUTINE report_parameters


END MODULE cable_param_module


