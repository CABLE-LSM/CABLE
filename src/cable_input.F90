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
! Purpose: Input module for CABLE offline version
!
! Contact: Bernard.Pak@csiro.au
!
! History: Developed by Gab Abramowitz
!          Rewritten for v2.0 for new input files (1x1 deg instead of CCAM ~2x2 deg)
!          LAI and casa-cnp nutrient inputs included in 1x1 deg file
!
! ==============================================================================
!
! MODULEs used: cable_abort_module
!               cable_def_types_mod
!               cable_IO_vars_module
!               cable_read_module
!               netcdf
!               casadimension
!               casavariable
!               phenvariable
!               cable_param_module
!               cable_checks_module
!               cable_radiation_module
!               cable_init_module
!
!==============================================================================

MODULE cable_input_module
  ! Note that any precision changes from r_1 to REAL(4) enable running with -r8
  !
  USE cable_abort_module,      ONLY: abort, nc_abort
  USE cable_def_types_mod
  USE casadimension,     ONLY: icycle
  USE casavariable
  USE casaparm, ONLY: forest, shrub
  USE phenvariable
  !! vh_js !!
  USE POP_Types,               ONLY: POP_TYPE
  USE POPLUC_Types,               ONLY: POPLUC_TYPE
  USE cable_param_module
  USE cable_checks_module,     ONLY: ranges, rh_sh
  USE cbl_sinbet_mod,  ONLY: sinbet
  USE cable_IO_vars_module
  USE cable_read_module,       ONLY: readpar
  USE cable_init_module
  USE netcdf ! link must be made in cd to netcdf-x.x.x/src/f90/netcdf.mod
  USE cable_common_module, ONLY : filename, cable_user, CurYear, is_leapyear
  USE casa_ncdf_module, ONLY: HANDLE_ERR
  USE casa_inout_module, ONLY: casa_readbiome, casa_readphen, casa_init

  IMPLICIT NONE

  PRIVATE
  PUBLIC get_default_lai, open_met_file, close_met_file,load_parameters,      &
       allocate_cable_vars, get_met_data, &
       ncid_met,        &
       ncid_rain,       &
       ncid_snow,       &
       ncid_lw,         &
       ncid_sw,         &
       ncid_ps,         &
       ncid_qa,         &
       ncid_ta,         &
       ncid_wd,         &
       ncid_mask

  INTEGER                      ::                                        &
       ncid_met,        & ! met data netcdf file ID
       ncid_rain,       & ! following are netcdf file IDs for gswp run
       ncid_snow,       &
       ncid_lw,         &
       ncid_sw,         &
       ncid_ps,         &
       ncid_qa,         &
       ncid_ta,         &
       ncid_wd,         &
       ncid_mask,       &
       ok                 ! netcdf error status
  ! - see ALMA compress by gathering
  INTEGER,POINTER,DIMENSION(:) :: landGrid ! for ALMA compressed variables
  REAL,POINTER,DIMENSION(:)    ::                                        &
       elevation,       & ! site/grid cell elevation
       avPrecip           ! site/grid cell average precip
  TYPE met_varID_type
     INTEGER                   ::                                        &
          SWdown,       &
          LWdown,       &
          Wind,         &
          Wind_E,       &
          PSurf,        &
          Tair,         &
          Qair,         &
          Rainf,        &
          Snowf,        &
          CO2air,       &
          Elev,         &
          LAI,          &
          avPrecip,     &
          iveg,         &
          isoil,        &
          patchfrac
  END TYPE met_varID_type
  TYPE(met_varID_type)              :: id ! netcdf variable IDs for input met variables
  TYPE met_units_type
     CHARACTER(LEN=20)              ::                                        &
          SWdown,       &
          LWdown,       &
          Wind,         &
          Wind_E,       &
          PSurf,        &
          Tair,         &
          Qair,         &
          Rainf,        &
          Snowf,        &
          CO2air,       &
          Elev,         &
          avPrecip
  END TYPE met_units_type
  TYPE(met_units_type)              :: metunits ! units for meteorological variables
  TYPE convert_units_type
     REAL                      ::                                        &
          PSurf,        &
          Tair,         &
          Qair,         &
          Rainf,        &
          CO2air,       &
          Elev
  END TYPE convert_units_type
  TYPE(convert_units_type)          :: convert ! units change factors for met variables

  !$OMP THREADPRIVATE(ok,exists)

CONTAINS

  !==============================================================================
  !
  ! Name: get_default_lai
  !
  ! Purpose: Reads all monthly LAI from default gridded netcdf file
  !
  ! CALLed from: load_parameters
  !
  ! CALLs: nc_abort
  !        abort
  !
  ! Input file: [LAI].nc
  !
  !==============================================================================

  SUBROUTINE get_default_lai

    ! Input variables:
    !   filename%LAI      - via cable_IO_vars_module
    !   landpt            - via cable_IO_vars_module (%nap,cstart,cend,ilon,ilat)
    !   exists%laiPatch   - via cable_IO_vars_module
    ! Output variables:
    !   defaultLAI(mp,12) - via cable_IO_vars_module

    ! Local variables
    INTEGER                              ::                                &
         ncid,                                   &
         xID,                                    &
         yID,                                    &
         pID,                                    &
         tID,                                    &
         laiID,                                  &
         nlon,                                   &
         nlat,                                   &
         nLaiPatches,                            &
         ntime,                                  &
         e, tt,i                                     ! do loop counter
    REAL, DIMENSION(:,:,:),  ALLOCATABLE :: inLai3D
    REAL, DIMENSION(:,:,:,:),ALLOCATABLE :: inLai4D

    ! Allocate default LAI variable: changed mland to mp (BP apr2010)
    ALLOCATE(defaultLAI(mp,12))   ! mp = mp_global

    WRITE(logn,*) ' Loading LAI from default file ', TRIM(filename%LAI)
    ! Open netcdf file
    ok = NF90_OPEN(filename%LAI,0,ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error opening default LAI file.')

    ok = NF90_INQ_DIMID(ncid,'x',xID)
    IF (ok /= NF90_NOERR) THEN   ! added to read some input files (BP mar2011)
       ok = NF90_INQ_DIMID(ncid,'longitude',xID)
       IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error inquiring x dimension.')
    END IF

    ok = NF90_INQUIRE_DIMENSION(ncid,xID,LEN=nlon)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error getting x dimension.')
    ok = NF90_INQ_DIMID(ncid,'y',yID)
    IF (ok /= NF90_NOERR) THEN   ! added to read some input files (BP mar2011)
       ok = NF90_INQ_DIMID(ncid,'latitude',yID)
       IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error inquiring y dimension.')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid,yID,LEN=nlat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error getting y dimension.')
    ok = NF90_INQ_DIMID(ncid,'patch',pID)
    IF(ok/=NF90_NOERR) THEN ! if failed
       exists%laiPatch = .FALSE.
       WRITE(logn,*) ' **ALL patches will be given the same LAI**'
    ELSE
       exists%laiPatch = .TRUE.
       ok = NF90_INQUIRE_DIMENSION(ncid,pID,len=nLaiPatches)
       IF(ANY(landpt(:)%nap>nLaiPatches)) THEN
          WRITE(logn,*) ' **Some patches will be given the same LAI**'
       END IF
       IF(nLaiPatches>max_vegpatches) THEN ! input file can have more info
          WRITE(*,*) ' Note that LAI input file has ', nLaiPatches, ' patches.'
          WRITE(*,*) ' while the model has max at ', max_vegpatches, ' patches.'
       END IF
    END IF
    ok = NF90_INQ_DIMID(ncid,'time',tID)
    IF (ok /= NF90_NOERR) THEN   ! added to read some input files (BP mar2011)
       ok = NF90_INQ_DIMID(ncid,'month',tID)
       IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error inquiring t dimension.')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid,tID,LEN=ntime)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error getting time dimension.')
    IF (ntime /= 12) CALL abort('Time dimension not 12 months.')

    ok = NF90_INQ_VARID(ncid,'LAI',laiID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding LAI variable.')

    ! Read LAI values:
    IF (exists%laiPatch) THEN
       ALLOCATE(inLai4D(nlon,nlat,nLaiPatches,ntime))
       ok = NF90_GET_VAR(ncid,laiID,inLai4D)
       IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading 4D LAI variable.')
       DO e = 1, mland  ! over all land grid points
          DO tt = 1, ntime
             defaultLAI(landpt(e)%cstart:landpt(e)%cend,tt) = &
                  inLai4D(landpt(e)%ilon,landpt(e)%ilat,1:landpt(e)%nap,tt)
          END DO
       END DO
    ELSE
       ALLOCATE(inLai3D(nlon,nlat,ntime))
       ok = NF90_GET_VAR(ncid,laiID,inLai3D)
       IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading 3D LAI variable.')
       DO e = 1, mland
          DO tt = 1, ntime
             defaultLAI(landpt(e)%cstart:landpt(e)%cend,tt) = &
                  inLai3D(landpt(e)%ilon,landpt(e)%ilat,tt)
          END DO
       END DO
    END IF

    ! Close netcdf file
    ok = NF90_CLOSE(ncid)


  END SUBROUTINE get_default_lai
  !==============================================================================
  !
  ! Name: open_met_file
  !
  ! Purpose: Opens netcdf file containing meteorological (LSM input) data
  !          and determines:
  !   1. Spatial details - number of sites/grid cells, latitudes, longitudes
  !   2. Timing details - time step size, number of timesteps, starting date,
  !      and whether time coordinate is local or GMT
  !   3. Checks availability, including units issues, of all required
  !      meteorological input variables. Also checks whether or not LAI is
  !      present, and fetches prescribed veg and soil type if present.
  !
  !
  ! CALLed from: cable_offline_driver
  !
  ! CALLs: abort
  !        nc_abort
  !        date_and_time
  !
  ! Input file: [SiteName].nc
  !             [GSWP_Snowf].nc
  !             [GSWP_LWdown].nc
  !             [GSWP_SWdown].nc
  !             [GSWP_PSurf].nc
  !             [GSWP_Qair].nc
  !             [GSWP_Tair].nc
  !             [GSWP_wind].nc
  !             [GSWP_Rainf].nc
  !
  !==============================================================================

  SUBROUTINE open_met_file(dels,koffset,kend,spinup, TFRZ)

    USE CABLE_COMMON_MODULE, ONLY : IS_LEAPYEAR
    USE casa_ncdf_module, ONLY: HANDLE_ERR, YMDHMS2DOYSOD, DOYSOD2YMDHMS
    USE CABLE_METUTILS_MODULE, ONLY: possible_varnames, find_metvarid

    IMPLICIT NONE
    ! Input arguments
    REAL, INTENT(OUT) :: dels   ! time step size
    REAL, INTENT(IN) :: TFRZ
    INTEGER, INTENT(INOUT)      :: koffset ! offset between met file and desired period
    INTEGER, INTENT(OUT)        :: kend   ! number of time steps in simulation
    LOGICAL, INTENT(IN)              :: spinup ! will a model spinup be performed?

    ! Local variables
    INTEGER                     ::                                         &
         timevarID,              & ! time variable ID number
         xdimID,                 & ! x dimension ID numbers
         ydimID,                 & ! y dimension ID numbers
         patchdimID,             & ! patch dimension ID
         monthlydimID,           & ! month dimension ID for LAI info
         maskID,                 & ! mask variable ID
         landID,                 & ! land variable ID
         landdimID,              & ! land dimension ID
         latitudeID,             & ! lat variable IDs
         longitudeID,            & ! lon variable IDs
         edoy,                   & ! end time day-of-year
         eyear,                  & ! end time year
         jump_days,              & ! days made by first "time" entry
         sdoytmp,                & ! used to determine start time hour-of-day
         mland_ctr,              & ! counter for number of land points read from file
         mland_fromfile,         & ! number of land points in file
         lai_dims,               & ! number of dims of LAI var if in met file
         iveg_dims,              & ! number of dims of iveg var if in met file
         isoil_dims,             & ! number of dims of isoil var if in met file
         tsmin,tsdoy,tsyear,     & ! temporary variables
         x,y,i,j,                & ! do loop counters
         tempmonth,                              &
         ssod, &
         nsod, &
         LOY, &
         iday,&
         imin,&
         isec,&
         ishod, &
         dnsec  = 0,&
         ntstp
    INTEGER,DIMENSION(1)        ::                                         &
         timedimID,              & ! time dimension ID number
         data1i                    ! temp variable for netcdf reading
    INTEGER,DIMENSION(4)        :: laidimids ! for checking lai variable
    INTEGER,DIMENSION(1,1)      :: data2i ! temp variable for netcdf reading
    INTEGER,POINTER,DIMENSION(:)     ::land_xtmp,land_ytmp ! temp indicies
    REAL,POINTER, DIMENSION(:)  :: lat_temp, lon_temp ! lat and lon
    REAL                        ::                                         &
         tshod,                  & ! temporary variable
         ehod,                   & ! end time hour-of-day
         precipTot,              & ! used for spinup adj
         avPrecipInMet             ! used for spinup adj
    CHARACTER(LEN=10)                :: todaydate, nowtime ! used to timestamp log file
    REAL(4),DIMENSION(1)             :: data1 ! temp variable for netcdf reading
    REAL(4),DIMENSION(1,1)           :: data2 ! temp variable for netcdf reading
    REAL(4), DIMENSION(:),     ALLOCATABLE :: temparray1  ! temp read in variable
    REAL(4), DIMENSION(:,:),   ALLOCATABLE :: &
         tempPrecip2,            & ! used for spinup adj
         temparray2                ! temp read in variable
    REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: tempPrecip3 ! used for spinup adj
    LOGICAL                          ::                                         &
         all_met,LAT1D,LON1D     ! ALL required met in met file (no synthesis)?

    ! Initialise parameter loading switch - will be set to TRUE when
    ! parameters are loaded:
    exists%parameters = .FALSE. ! initialise
    ! Initialise initialisation loading switch - will be set to TRUE when
    ! initialisation data are loaded:
    exists%initial = .FALSE. ! initialise

    LAT1D = .FALSE.
    LON1D = .FALSE.

    ! Write filename to log file:
    WRITE(logn,*) '============================================================'
    WRITE(logn,*) 'Log file for offline CABLE run:'
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate=todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime=nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    WRITE(logn,*) TRIM(nowtime),' ',TRIM(todaydate)
    WRITE(logn,*) '============================================================'

    ! Open netcdf file:
    IF (ncciy > 0) THEN

       IF( globalMetfile%l_gpcc ) THEN
       WRITE(logn,*) 'Opening met data file: ', TRIM(globalMetfile%rainf), ' and 7 more'
       ELSE
         WRITE(logn,*) 'Opening met data file: ', TRIM(gswpfile%rainf), ' and 7 more'
       ENDIF
 
       IF( globalMetfile%l_gpcc ) THEN
       ok = NF90_OPEN(globalMetfile%rainf,0,ncid_rain)
       ELSE
         ok = NF90_OPEN(gswpfile%rainf,0,ncid_rain)
       ENDIF
       IF (ok /= NF90_NOERR) THEN
          PRINT*,'rainf'
          CALL handle_err( ok )
       ENDIF
       IF(.NOT. globalMetfile%l_gpcc) THEN
         ok = NF90_OPEN(gswpfile%snowf,0,ncid_snow)
         IF (ok /= NF90_NOERR) THEN
            PRINT*,'snow'
            CALL handle_err( ok )
         ENDIF
       ENDIF
       
       IF( globalMetfile%l_gpcc ) THEN
       ok = NF90_OPEN(globalMetfile%LWdown,0,ncid_lw)
       ELSE
         ok = NF90_OPEN(gswpfile     %LWdown,0,ncid_lw)
       ENDIF
       IF (ok /= NF90_NOERR) THEN
          PRINT*,'lw'
          CALL handle_err( ok )
       ENDIF
       
       IF( globalMetfile%l_gpcc ) THEN
       ok = NF90_OPEN(globalMetfile%SWdown,0,ncid_sw)
       ELSE
         ok = NF90_OPEN(gswpfile%SWdown,0,ncid_sw)
       ENDIF
       IF (ok /= NF90_NOERR) THEN
          PRINT*,'sw'
          CALL handle_err( ok )
       ENDIF
       
       IF( globalMetfile%l_gpcc ) THEN
       ok = NF90_OPEN(globalMetfile%PSurf,0,ncid_ps)
       ELSE
         ok = NF90_OPEN(gswpfile%PSurf,0,ncid_ps)
       ENDIF
       IF (ok /= NF90_NOERR) THEN
          PRINT*,'ps'
          CALL handle_err( ok )
       ENDIF
       
       IF( globalMetfile%l_gpcc ) THEN
       ok = NF90_OPEN(globalMetfile%Qair,0,ncid_qa)
       ELSE
         ok = NF90_OPEN(gswpfile%Qair,0,ncid_qa)
       ENDIF
       IF (ok /= NF90_NOERR) THEN
          PRINT*,'qa'
          CALL handle_err( ok )
       ENDIF
       
       IF( globalMetfile%l_gpcc ) THEN
       ok = NF90_OPEN(globalMetfile%Tair,0,ncid_ta)
       ELSE
         ok = NF90_OPEN(gswpfile%Tair,0,ncid_ta)
       ENDIF
       IF (ok /= NF90_NOERR) THEN
          PRINT*,'ta'
          CALL handle_err( ok )
       ENDIF
       
       IF( globalMetfile%l_gpcc ) THEN
       ok = NF90_OPEN(globalMetfile%wind,0,ncid_wd)
       ELSE
         ok = NF90_OPEN(gswpfile%wind,0,ncid_wd)
       ENDIF
       IF (ok /= NF90_NOERR) THEN
          PRINT*,'wind',ncid_wd
          CALL handle_err( ok )
       ENDIF
       IF (cable_user%GSWP3) THEN
          ok = NF90_OPEN(gswpfile%mask,0,ncid_mask)
          IF (ok .NE. NF90_NOERR) THEN
             CALL nc_abort(ok, "Error opening GSWP3 mask file")
          END IF
          LAT1D = .TRUE.   !GSWP3 forcing has 1d lat/lon variables
          LON1D = .TRUE.
       ELSE
          ncid_mask = ncid_rain
       END IF
       ncid_met = ncid_rain
    ELSE
       WRITE(logn,*) 'Opening met data file: ', TRIM(filename%met)
       ok = NF90_OPEN(filename%met,0,ncid_met) ! open met data file
       IF (ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error opening netcdf met forcing file '//TRIM(filename%met)// &
            ' (SUBROUTINE open_met_file)')
    ENDIF

    !!=====================VV Determine spatial details VV=================
    ! Determine number of sites/gridcells.
    ! Find size of 'x' or 'lon' dimension:
    ok = NF90_INQ_DIMID(ncid_met,'x', xdimID)
    IF(ok/=NF90_NOERR) THEN ! if failed
       ! Try 'lon' instead of x
       ok = NF90_INQ_DIMID(ncid_met,'lon', xdimID)
       IF(ok/=NF90_NOERR) CALL nc_abort &
            (ok,'Error finding x dimension in '&
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid_met,xdimID,len=xdimsize)
    IF(ok/=NF90_NOERR) CALL nc_abort &
         (ok,'Error determining size of x dimension in ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Find size of 'y' dimension:
    ok = NF90_INQ_DIMID(ncid_met,'y', ydimID)
    IF(ok/=NF90_NOERR) THEN ! if failed
       ! Try 'lat' instead of y
       ok = NF90_INQ_DIMID(ncid_met,'lat', ydimID)
       IF(ok/=NF90_NOERR) CALL nc_abort &
            (ok,'Error finding y dimension in ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid_met,ydimID,len=ydimsize)
    IF(ok/=NF90_NOERR) CALL nc_abort &
         (ok,'Error determining size of y dimension in ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Determine number of gridcells in netcdf file:
    ngridcells = xdimsize*ydimsize
    WRITE(logn,'(A28,I7)') 'Total number of gridcells: ', ngridcells

    ! Get all latitude and longitude values.
    ! Find latitude variable (try 'latitude' and 'nav_lat'(ALMA)):
    CALL find_metvarid(ncid_met, possible_varnames%LatNames, latitudeID, ok)

    IF (ok /= NF90_NOERR) CALL nc_abort &
      (ok,'Error finding latitude variable in ' &
      //TRIM(filename%met)//' (SUBROUTINE open_met_file)')

    ! Allocate space for lat_all variable and its temp counterpart:
    ALLOCATE(lat_all(xdimsize,ydimsize))
    ALLOCATE(temparray2(xdimsize,ydimsize))
    !MDeck allow for 1d lat called 'lat'
    ! Get latitude values for entire region:
    IF (.NOT.LAT1D) THEN
       ok= NF90_GET_VAR(ncid_met,latitudeID,temparray2)
    ELSE
       IF (ALLOCATED(temparray1)) DEALLOCATE(temparray1)
       ALLOCATE(temparray1(ydimsize))
       ok= NF90_GET_VAR(ncid_met,latitudeID,temparray1)
       temparray2 = SPREAD(temparray1,1,xdimsize)
       DEALLOCATE(temparray1)
    END IF
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error reading latitude variable in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Needed since r_1 will be double precision with -r8:
    lat_all = REAL(temparray2)

    ! Find longitude variable (try 'longitude' and 'nav_lon'(ALMA)):
    CALL find_metvarid(ncid_met, possible_varnames%LonNames, longitudeID, ok)

    IF(ok /= NF90_NOERR) CALL nc_abort &
      (ok,'Error finding longitude variable in ' &
      //TRIM(filename%met)//' (SUBROUTINE open_met_file)')

      ! Allocate space for lon_all variable:
    ALLOCATE(lon_all(xdimsize,ydimsize))
    ! Get longitude values for entire region:
    !MDeck allow for 1d lon
    IF (.NOT.LON1D) THEN
       ok= NF90_GET_VAR(ncid_met,longitudeID,temparray2)
    ELSE
       ALLOCATE(temparray1(xdimsize))
       ok= NF90_GET_VAR(ncid_met,longitudeID,temparray1)
       temparray2 = SPREAD(temparray1,2,ydimsize)
       DEALLOCATE(temparray1)
    END IF
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error reading longitude variable in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Needed since r_1 will be double precision with -r8:
    lon_all = REAL(temparray2)
    DEALLOCATE(temparray2)

    ! Check for "mask" variable or "land" variable to tell grid type
    ! (and allow neither if only one gridpoint). "mask" is a 2D variable
    ! with dims x,y and "land" is a 1D variable.
    CALL find_metvarid(ncid_mask, possible_varnames%MaskNames, maskID, ok)
    IF(ok /= NF90_NOERR) THEN ! if error, i.e. no "mask" variable:
       ! Check for "land" variable:
       ok = NF90_INQ_VARID(ncid_met, 'land', landID)
       IF(ok /= NF90_NOERR) THEN ! ie no "land" or "mask"
          IF(ngridcells==1) THEN
             ! Allow no explicit grid system if only one gridpoint
             ALLOCATE(mask(xdimsize,ydimsize)) ! Allocate "mask" variable
             metGrid='mask' ! Use mask system, one gridpoint.
             mask = 1
             ALLOCATE(latitude(1),longitude(1))
             latitude = lat_all(1,1)
             longitude = lon_all(1,1)
             mland_fromfile=1
             ALLOCATE(land_x(mland_fromfile),land_y(mland_fromfile))
             land_x = 1
             land_y = 1
          ELSE
             ! Call abort if more than one gridcell and no
             ! recognised grid system:
             CALL nc_abort &
                  (ok,'Error finding grid system ("mask" or "land") variable in ' &
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          END IF
       ELSE ! i.e. "land" variable exists
          metGrid='land'
          ! Check size of "land" dimension:
          ok = NF90_INQ_DIMID(ncid_met,'land', landdimID)
          IF(ok/=NF90_NOERR) CALL nc_abort &
               (ok,'Error finding land dimension in ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ok = NF90_INQUIRE_DIMENSION(ncid_met,landdimID,len=mland_fromfile)
          IF(ok/=NF90_NOERR) CALL nc_abort &
               (ok,'Error determining size of land dimension in ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ! Allocate landGrid variable and its temporary counterpart:
          ALLOCATE(landGrid(mland_fromfile))
          ALLOCATE(temparray1(mland_fromfile))
          ! Get values of "land" variable from file:
          ok= NF90_GET_VAR(ncid_met,landID,temparray1)
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading "land" variable in ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ! Needed since r_1 will be double precision with -r8:
          landGrid = REAL(temparray1)
          DEALLOCATE(temparray1)
          ! Allocate latitude and longitude variables:
          ALLOCATE(latitude(mland_fromfile),longitude(mland_fromfile))
          ! Write to indicies of points in all-grid which are land
          ALLOCATE(land_x(mland_fromfile),land_y(mland_fromfile))
          ! Allocate "mask" variable:
          ALLOCATE(mask(xdimsize,ydimsize))
          ! Initialise all gridpoints as sea:
          mask = 0
          DO j=1, mland_fromfile ! over all land points
             ! Find x and y coords of current land point
             y = INT((landGrid(j)-1)/xdimsize)
             x = landGrid(j) - y * xdimsize
             y=y+1
             ! Write lat and lon to land-only lat/lon vars:
             latitude(j) = lat_all(x,y)
             longitude(j) = lon_all(x,y)
             ! Write to mask variable:
             mask(x,y)=1
             ! Save indicies:
             land_x(j) = x
             land_y(j) = y
          END DO
       END IF ! does "land" variable exist
    ELSE ! i.e. "mask" variable exists
       ! Allocate "mask" variable:
       ALLOCATE(mask(xdimsize,ydimsize))
       metGrid='mask' ! Use mask system
       ! Get mask values from file:
       ok= NF90_GET_VAR(ncid_mask,maskID,mask)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading "mask" variable in ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       !gswp3 uses 1 for sea and 0 for land, make is opposite
       IF (cable_user%gswp3) mask = 1-mask

       ! Allocate space for extracting land lat/lon values:
       ALLOCATE(lat_temp(ngridcells),lon_temp(ngridcells))
       ! Allocate space for extracting index of mask which is land
       ALLOCATE(land_xtmp(ngridcells),land_ytmp(ngridcells))
       ! Cycle through all gridsquares:
       mland_ctr = 0 ! initialise
       DO y=1,ydimsize
          DO x=1,xdimsize
             IF(mask(x,y)==1) THEN ! If land
                mland_ctr = mland_ctr + 1
                ! Store lat and lon for land points
                lat_temp(mland_ctr) = lat_all(x,y)
                lon_temp(mland_ctr) = lon_all(x,y)
                ! Store indicies of points in mask which are land
                land_xtmp(mland_ctr) = x
                land_ytmp(mland_ctr) = y
             END IF
          END DO
       END DO
       ! Record number of land points
       mland_fromfile = mland_ctr
       ! Allocate latitude and longitude variables:
       ALLOCATE(latitude(mland_fromfile),longitude(mland_fromfile))
       ! Write to latitude and longitude variables:
       latitude = lat_temp(1:mland_fromfile)
       longitude = lon_temp(1:mland_fromfile)
       ! Write to indicies of points in mask which are land
       ALLOCATE(land_x(mland_fromfile),land_y(mland_fromfile))
       land_x = land_xtmp(1:mland_fromfile)
       land_y = land_ytmp(1:mland_fromfile)
       ! Clear lon_temp, lat_temp,land_xtmp,land_ytmp
       DEALLOCATE(lat_temp,lon_temp,land_xtmp,land_ytmp)
    END IF ! "mask" variable or no "mask" variable

    ! Set global mland value (number of land points), used to allocate
    ! all of CABLE's arrays:
    mland = mland_fromfile

    ! Write number of land points to log file:
    WRITE(logn,'(24X,I7,A29)') mland_fromfile, ' of which are land grid cells'

    ! Check if veg/soil patch dimension exists (could have
    ! parameters with patch dimension)
    ok = NF90_INQ_DIMID(ncid_met,'patch', patchdimID)
    IF(ok/=NF90_NOERR) THEN ! if failed
       exists%patch = .FALSE.
       nmetpatches = 1       ! initialised so that old met files without patch
       ! data can still be run correctly (BP apr08)
    ELSE ! met file does have patch dimension
       exists%patch = .TRUE.
       ok = NF90_INQUIRE_DIMENSION(ncid_met,patchdimID,len=nmetpatches)
    END IF

    ! Check if monthly dimension exists for LAI info
    ok = NF90_INQ_DIMID(ncid_met,'monthly', monthlydimID)
    IF(ok==NF90_NOERR) THEN ! if found
       ok = NF90_INQUIRE_DIMENSION(ncid_met,monthlydimID,len=tempmonth)
       IF(tempmonth/=12) CALL abort ('Number of months in met file /= 12.')
    END IF

    ! Set longitudes to be [-180,180]:
    WHERE(longitude>180.0)
       longitude = longitude - 360.0
    END WHERE
    ! Check ranges for latitude and longitude:
    IF(ANY(longitude>180.0).OR.ANY(longitude<-180.0)) &
         CALL abort('Longitudes read from '//TRIM(filename%met)// &
         ' are not [-180,180] or [0,360]! Please set.')
    IF(ANY(latitude>90.0).OR.ANY(latitude<-90.0)) &
         CALL abort('Latitudes read from '//TRIM(filename%met)// &
         ' are not [-90,90]! Please set.')

    !!=================^^ End spatial details ^^========================

    !!=========VV Determine simulation timing details VV================
    ! Inquire 'time' variable's ID:
    CALL find_metvarid(ncid_met, possible_varnames%TimeNames, timevarID, ok)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding time variable in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get ID for dimension upon which time depends:
    ok = NF90_INQUIRE_VARIABLE(ncid_met,timevarID,dimids=timedimID)
    IF(ok/=NF90_NOERR) CALL nc_abort &
         (ok,'Error determining "time" dimension dimension in ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Determine number of time steps:
    ok = NF90_INQUIRE_DIMENSION(ncid_met,timedimID(1),len=kend)
    IF(ok/=NF90_NOERR) CALL nc_abort &
         (ok,'Error determining number of timesteps in ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Allocate time variable:
    ALLOCATE(timevar(kend))
    ! Fetch 'time' variable:
    ok= NF90_GET_VAR(ncid_met,timevarID,timevar)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error reading time variable in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF (cable_user%gswp3) THEN         !Hack the GSWP3 time units to make from start of year
       timevar(:) = (timevar(:)-timevar(1))*3600.0 + 1.5*3600.0  !convert hours to seconds
    END IF
    ! Set time step size:
    dels = REAL(timevar(2) - timevar(1))
    WRITE(logn,'(1X,A29,I8,A3,F10.3,A5)') 'Number of time steps in run: ',&
         kend,' = ', REAL(kend)/(3600/dels*24),' days'


    ! CLN READJUST kend referring to Set START & END
    ! if kend > # days in selected episode


    !********* gswp input file has bug in timevar **************
    IF (ncciy > 0) THEN
       PRINT *, 'original timevar(kend) = ', timevar(kend)
       DO i = 1, kend - 1
          timevar(i+1) = timevar(i) + dels
       ENDDO
       PRINT *, 'New      timevar(kend) = ', timevar(kend)
       !! hacking (BP feb2011)
       !      kend = 16   ! 2 days for 1986
       !!      kend = 480  ! 2 months for 1986
       !      PRINT *, 'Hacked   timevar(kend) = ', timevar(kend)
       !      PRINT *, 'Hacked kend = ', kend
       !! end hacking
    END IF
    !********* done bug fixing for timevar in gswp input file **

    ! Write time step size to log file:
    WRITE(logn,'(1X,A17,F8.1,1X,A7)') 'Time step size:  ', dels, 'seconds'
    ! Get units for 'time' variable:
    ok = NF90_GET_ATT(ncid_met,timevarID,'units',timeunits)
    IF (.NOT.cable_user%GSWP3) THEN
       ok = NF90_GET_ATT(ncid_met,timevarID,'units',timeunits)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding time variable units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ELSE
       !Hack the GSWP3 time units to make from start of year
       WRITE(*,*) 'writing timeunits'
       WRITE (timeunits, "('seconds since ',I4.4,'-01-01 00:00:00')") ncciy
       WRITE(*,*) 'wrote time units'
    END IF
    !MDeck
    !write(*,*) timeunits
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding time variable units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')

    !****** PALS met file has timevar(1)=0 while timeunits from 00:30:00 ******
    !!CLN CRITICAL! From my point of view, the information in the file is correct...
    !!CLN WHY DO the input files all have bugs???
    IF (timevar(1) == 0.0) THEN
       READ(timeunits(29:30),*) tsmin
       IF (tsmin*60.0 >= dels) THEN
          tsmin = tsmin - INT(dels / 60)
          timevar = timevar + dels
          WRITE(timeunits(29:30),'(i2.2)') tsmin
       ENDIF
    ENDIF
    !****** done bug fixing for timevar in PALS met file **********************

    !********* gswp input file has bug in timeunits ************
    IF (ncciy > 0) WRITE(timeunits(26:27),'(i2.2)') 0
    !********* done bug fixing for timeunits in gwsp file ******
    WRITE(logn,*) 'Time variable units: ', timeunits
    ! Get coordinate field:
    ok = NF90_GET_ATT(ncid_met,timevarID,'coordinate',time_coord)
    ! If error getting coordinate field (i.e. it doesn't exist):
    IF(ok /= NF90_NOERR) THEN
       ! Assume default time coordinate:
       IF(mland_fromfile==1.AND.(TRIM(cable_user%MetType) .NE. 'gswp')) THEN ! If single site, this is local time
          time_coord = 'LOC' ! 12am is 12am local time, at site/gridcell
       ELSE ! If multiple/global/regional, use GMT
          time_coord = 'GMT' ! 12am is GMT time, local time set by longitude
       END IF
    ELSE IF((ok==NF90_NOERR.AND.time_coord=='LOC'.AND.mland_fromfile>1)) THEN
       ! Else if local time is selected for regional simulation, abort:
       CALL abort('"time" variable must be GMT for multiple site simulation!' &
            //' Check "coordinate" field in time variable.' &
            //' (SUBROUTINE open_met_file)')
    ELSE IF(time_coord/='LOC'.AND.time_coord/='GMT') THEN
       CALL abort('Meaningless time coordinate in met data file!' &
            // ' (SUBROUTINE open_met_file)')
    END IF

    ! Use internal files to convert "time" variable units (giving the run's
    ! start time) from character to integer; calculate starting hour-of-day,
    ! day-of-year, year:
    IF (.NOT.cable_user%GSWP3) THEN
       READ(timeunits(15:18),*) syear
       READ(timeunits(20:21),*) smoy ! integer month
       READ(timeunits(23:24),*) sdoytmp ! integer day of that month
       READ(timeunits(26:27),*) shod  ! starting hour of day
    ELSE
       syear=ncciy
       smoy=1
       sdoytmp=1
       shod=0
    END IF

    ! if site data, shift start time to middle of timestep
    ! only do this if not already at middle of timestep
    !! vh_js !!
    IF ((TRIM(cable_user%MetType).EQ.'' .OR. &
         TRIM(cable_user%MetType).EQ.'site').AND. MOD(shod*3600, dels)==0 .AND. &
         (shod.GT.dels/3600./2.) ) THEN
       shod = shod - dels/3600./2.
    ELSEIF (TRIM(cable_user%MetType).EQ.''.OR. &
         (TRIM(cable_user%MetType).EQ.'site') .AND. MOD(shod*3600, dels)==0 .AND. &
         (shod.LT.dels/3600./2.) ) THEN
       shod = shod + dels/3600./2.
    ENDIF

    ! Decide day-of-year for non-leap year:
    CALL YMDHMS2DOYSOD( syear, smoy, sdoytmp, INT(shod), 0, 0, sdoy, ssod )
    ! Number of days between start position and 1st timestep:
    sdoy = sdoy + INT((timevar(1)/3600.0 + shod)/24.0)
    nsod = MOD(INT((timevar(1) + shod*3600)),86400)

    DO
       LOY = 365
       IF ( IS_LEAPYEAR( syear ) ) LOY = 366
       IF ( sdoy .GT. LOY ) THEN
          sdoy  = sdoy - LOY
          syear = syear + 1
       ELSE
          EXIT
       END IF
    END DO


    CALL DOYSOD2YMDHMS( syear, sdoy, nsod, smoy, iday, ishod, imin, isec )
    shod = REAL(ishod) + REAL(imin)/60. + REAL(isec)/3600.
    ! Cycle through days to find leap year inclusive starting date:
    ! Now all start time variables established, report to log file:
    WRITE(logn,'(1X,A12,F5.2,A14,I3,A14,I4,2X,A3,1X,A4)') &
         'Run begins: ',shod,' hour-of-day, ',sdoy, ' day-of-year, ',&
         syear, time_coord, 'time'
    ! Determine ending time of run...
    IF(leaps) THEN ! If we're using leap year timing...
       eyear = syear
       edoy  = sdoy + INT(((timevar(kend)-timevar(1))/3600.0 + shod)/24.0)
       ehod  = MOD(((timevar(kend)-timevar(1)/3600.) + shod),24._r_2)

       DO
          LOY = 365
          IF ( IS_LEAPYEAR( eyear ) ) LOY = 366
          IF ( edoy .GT. LOY ) THEN
             edoy  = edoy - LOY
             eyear = eyear + 1
          ELSE
             EXIT
          END IF
       END DO


    ELSE ! if not using leap year timing
       ! Update shod, sdoy, syear for first "time" value:
       ehod = MOD(REAL((timevar(kend)-timevar(1))/3600.0 + shod),24.0)
       edoy = MOD(INT(((timevar(kend)-timevar(1))/3600.0 + shod)/24.0) &
            + sdoy, 365)
       eyear = INT(REAL(INT(((timevar(kend)-timevar(1)) &
            /3600.0+shod)/24.0)+sdoy)/365.0)+syear
       !       ehod = MOD(REAL((timevar(kend)-timevar(1)+dels)/3600.0 + shod),24.0)
       !       edoy = MOD(INT(((timevar(kend)-timevar(1)+dels)/3600.0 + shod)/24.0) &
       !            + sdoy, 365)
       !       eyear = INT(REAL(INT(((timevar(kend)-timevar(1)+dels) &
       !            /3600.0+shod)/24.0)+sdoy)/365.0)+syear
    END IF
    ! IF A CERTAIN PERIOD IS DESIRED AND WE ARE NOT RUNNING ON GSWP DATA or special site
    ! RECALCULATE STARTING AND ENDING INDICES
    IF ( CABLE_USER%YEARSTART .GT. 0 .AND. .NOT. ncciy.GT.0 .AND. &
         TRIM(cable_user%MetType) .NE. "site") THEN
       IF ( syear.GT.CABLE_USER%YEARSTART .OR. eyear.LE.CABLE_USER%YEAREND .OR. &
            ( syear.EQ.CABLE_USER%YEARSTART .AND. sdoy.GT.1 ) ) THEN
          WRITE(*,*) "Chosen period doesn't match dataset period!"
          WRITE(*,*) "Chosen period: ",CABLE_USER%YEARSTART,1,CABLE_USER%YEAREND,365
          WRITE(*,*) "Data   period: ",syear,sdoy, eyear,edoy
          WRITE(*,*) "For using the metfile's time set CABLE_USER%YEARSTART = 0 !"
          STOP
       ENDIF

       ! Find real kstart!
       dnsec = 0
       DO y = syear, CABLE_USER%YEARSTART-1
          LOY = 365
          IF ( IS_LEAPYEAR( y ) ) LOY = 366
          IF ( y .EQ. syear ) THEN
             dnsec = ( LOY - sdoy ) * 86400 + (24 - shod) * 3600
          ELSE
             dnsec = dnsec + LOY * 86400
          ENDIF
       END DO
       koffset = INT(REAL(dnsec)/REAL(dels)) - 1
       ! Find real kend
       kend = 0
       DO y = CABLE_USER%YEARSTART, CABLE_USER%YEAREND
          LOY = 365
          IF ( IS_LEAPYEAR( y ) ) LOY = 366
          kend = kend + INT( REAL(LOY) * 86400./REAL(dels) )
       END DO
    ENDIF

    ! Report finishing time to log file:
    WRITE(logn,'(1X,A12,F5.2,A14,I3,A14,I4,2X,A3,1X,A4)') 'Run ends:   ',&
         ehod,' hour-of-day, ',edoy, &
         ' day-of-year, ', eyear, time_coord, 'time'
    !!===================^^ End timing details ^^==========================

    !!===================VV Look for met variables VV======================
    all_met = .TRUE. ! initialise
    ! Look for SWdown (essential):- - - - - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_sw

   ! option was added by Chris Lu to allow for different variable names between GPCC and GSWP forcings
   ! added by ypwang 30/oct/2012 
    CALL find_metvarid(ncid_met, possible_varnames%SWdownNames, id%SWdown, ok)

    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding SWdown in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get SWdown units and check okay:
    ok = NF90_GET_ATT(ncid_met,id%SWdown,'units',metunits%SWdown)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding SWdown units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    !! vh_js !! fixed bug in logic
    IF(.NOT.(metunits%SWdown(1:4)/='W/m2'.OR.metunits%SWdown(1:5) &
         /='W/m^2'.OR.metunits%SWdown(1:5)/='Wm^-2' &
         .OR.metunits%SWdown(1:4)/='Wm-2'.OR.metunits%SWdown(1:5) /= 'W m-2')) THEN
       WRITE(*,*) metunits%SWdown
       CALL abort('Unknown units for SWdown'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Look for Tair (essential):- - - - - - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_ta
    CALL find_metvarid(ncid_met, possible_varnames%TairNames, id%Tair, ok)

    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Tair in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get Tair units and check okay:
    ok = NF90_GET_ATT(ncid_met,id%Tair,'units',metunits%Tair)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Tair units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Tair(1:1)=='C'.OR.metunits%Tair(1:1)=='c') THEN
       ! Change from celsius to kelvin:
       convert%Tair = tfrz
       WRITE(logn,*) 'Temperature will be converted from C to K'
    ELSE IF(metunits%Tair(1:1)=='K'.OR.metunits%Tair(1:1)=='k') THEN
       ! Units are correct
       convert%Tair = 0.0
    ELSE
       WRITE(*,*) metunits%Tair
       CALL abort('Unknown units for Tair'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Look for Qair (essential):- - - - - - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_qa
    Call find_metvarid(ncid_met, possible_varnames%QairNames, id%Qair, ok)

    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Qair in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get Qair units:
    ok = NF90_GET_ATT(ncid_met,id%Qair,'units',metunits%Qair)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Qair units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Qair(1:1)=='%'.OR.metunits%Qair(1:1)=='-') THEN
       ! Change from relative humidity to specific humidity:
       convert%Qair = -999.0
       WRITE(logn,*) 'Humidity will be converted from relative to specific'
    ELSE IF(metunits%Qair(1:3)=='g/g'.OR.metunits%Qair(1:5)=='kg/kg' &
         .OR.metunits%Qair(1:3)=='G/G'.OR.metunits%Qair(1:5)=='KG/KG'.OR.metunits%Qair(1:7)=='kg kg-1') THEN
       ! Units are correct
       convert%Qair=1.0
    ELSE
       WRITE(*,*) metunits%Qair
       CALL abort('Unknown units for Qair'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF

    ! Look for Rainf (essential):- - - - - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_rain
    CALL find_metvarid(ncid_met, possible_varnames%RainNames, id%Rainf, ok)

    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Rainf in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get Rainf units:
    ok = NF90_GET_ATT(ncid_met,id%Rainf,'units',metunits%Rainf)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Rainf units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Rainf(1:8)=='kg/m^2/s'.OR.metunits%Rainf(1:6)=='kg/m2s'.OR.metunits%Rainf(1:10)== &
         'kgm^-2s^-1'.OR.metunits%Rainf(1:4)=='mm/s'.OR. &
         metunits%Rainf(1:6)=='mms^-1'.OR. &
         metunits%Rainf(1:7)=='kg/m^2s'.OR.metunits%Rainf(1:10)=='kg m-2 s-1'.OR.metunits%Wind(1:5)/='m s-1') THEN
       ! Change from mm/s to mm/time step:
       convert%Rainf = dels
    ELSE IF(metunits%Rainf(1:4)=='mm/h'.OR.metunits%Rainf(1:6)== &
         'mmh^-1') THEN
       ! Change from mm/h to mm/time step:
       convert%Rainf = dels/3600.0
    ELSE
       WRITE(*,*) metunits%Rainf
       CALL abort('Unknown units for Rainf'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Multiply acceptable Rainf ranges by time step size:
    !ranges%Rainf = ranges%Rainf*dels ! range therefore depends on dels ! vh ! why has this been commented out?
    ranges%Rainf = ranges%Rainf*dels ! range therefore depends on dels
    ! Look for Wind (essential):- - - - - - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_wd
    CALL find_metvarid(ncid_met, possible_varnames%WindNames, id%Wind, ok)

    IF(ok /= NF90_NOERR) THEN
       ! Look for vector wind:
       ok = NF90_INQ_VARID(ncid_met,'Wind_N',id%Wind)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding Wind_N in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       ok = NF90_INQ_VARID(ncid_met,'Wind_E',id%Wind_E)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding Wind_E in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       exists%Wind = .FALSE. ! Use vector wind when reading met
    ELSE
       exists%Wind = .TRUE. ! 'Wind' variable exists
    END IF

    ! The following does not work with vector winds. Do we want to keep
    ! vector winds?
    ! Get Wind units:
    ok = NF90_GET_ATT(ncid_met,id%Wind,'units',metunits%Wind)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Wind units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF (metunits%Wind(1:3)/='m/s'.AND.metunits%Wind(1:2)/='ms'.AND.metunits%Wind(1:5)/='m s-1') THEN
       WRITE(*,*) metunits%Wind
       CALL abort('Unknown units for Wind'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Now "optional" variables:
    ! Look for LWdown (can be synthesised):- - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_lw
    CALL find_metvarid(ncid_met, possible_varnames%LWdownNames, id%LWdown, ok)

    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%LWdown = .TRUE. ! LWdown is present in met file
       ! Get LWdown units and check okay:
       ok = NF90_GET_ATT(ncid_met,id%LWdown,'units',metunits%LWdown)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding LWdown units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       !! vh_js !! fixed bug in logic
!!$       IF(metunits%LWdown(1:4)/='W/m2'.AND.metunits%LWdown(1:5) &
!!$            /='W/m^2'.AND.metunits%LWdown(1:5)/='Wm^-2' &
!!$            .AND.metunits%LWdown(1:4)/='Wm-2') THEN
       IF(.NOT.(metunits%LWdown(1:4)/='W/m2'.OR.metunits%LWdown(1:5) &
            /='W/m^2'.OR.metunits%LWdown(1:5)/='Wm^-2' &
            .OR.metunits%LWdown(1:4)/='Wm-2'.OR.metunits%SWdown(1:5) /= 'W m-2')) THEN

          WRITE(*,*) metunits%LWdown
          CALL abort('Unknown units for LWdown'// &
               ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE
       exists%LWdown = .FALSE. ! LWdown is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Note this in log file:
       WRITE(logn,*) 'LWdown not present in met file; ', &
            'values will be synthesised based on air temperature.'
    END IF
    ! Look for PSurf (can be synthesised):- - - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_ps
    CALL find_metvarid(ncid_met, possible_varnames%PSurfNames, id%PSurf, ok)

    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%PSurf = .TRUE. ! PSurf is present in met file
       ! Get PSurf units and check:
       ok = NF90_GET_ATT(ncid_met,id%PSurf,'units',metunits%PSurf)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding PSurf units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       IF(metunits%PSurf(1:2)=='Pa'.OR.metunits%PSurf(1:2)=='pa'.OR. &
            metunits%PSurf(1:2)=='PA' ) THEN
          ! Change from pa to mbar (cable uses mbar):
          convert%PSurf = 0.01
          WRITE(logn,*) 'Pressure will be converted from Pa to mb'
       ELSE IF(metunits%PSurf(1:2)=='KP'.OR.metunits%PSurf(1:2)=='kP' &
            .OR.metunits%PSurf(1:2)=='Kp'.OR.metunits%PSurf(1:2)=='kp') THEN
          ! convert from kPa to mb
          convert%PSurf = 10.0
          WRITE(logn,*) 'Pressure will be converted from kPa to mb'
       ELSE IF(metunits%PSurf(1:2)=='MB'.OR.metunits%PSurf(1:2)=='mB' &
            .OR.metunits%PSurf(1:2)=='Mb'.OR.metunits%PSurf(1:2)=='mb') THEN
          ! Units are correct
          convert%PSurf = 1.0
       ELSE
          WRITE(*,*) metunits%PSurf
          CALL abort('Unknown units for PSurf'// &
               ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE ! If PSurf not present
       exists%PSurf = .FALSE. ! PSurf is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Look for "elevation" variable to approximate pressure based
       ! on elevation and temperature:
       CALL find_metvarid(ncid_met, possible_varnames%ElevNames, id%Elev, ok)
       IF(ok == NF90_NOERR) THEN ! elevation present
          ! Get elevation units:
          ok = NF90_GET_ATT(ncid_met,id%Elev,'units',metunits%Elev)
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error finding elevation units in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ! Units should be metres or feet:
          IF(metunits%Elev(1:1)=='m'.OR.metunits%Elev(1:1)=='M') THEN
             ! This is the expected unit - metres
             convert%Elev = 1.0
          ELSE IF(metunits%Elev(1:1)=='f'.OR.metunits%Elev(1:1)=='F') THEN
             ! Convert from feet to metres:
             convert%Elev = 0.3048
          ELSE
             CALL abort('Unknown units for Elevation'// &
                  ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
          END IF
          ! Allocate space for elevation variable:
          ALLOCATE(elevation(mland))
          ! Get site elevations:
          IF(metGrid=='mask') THEN
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%Elev,data2, &
                     start=(/land_x(i),land_y(i)/),count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading elevation in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                elevation(i)=REAL(data2(1,1))*convert%Elev
             END DO
          ELSE IF(metGrid=='land') THEN
             ! Collect data from land only grid in netcdf file:
             ok= NF90_GET_VAR(ncid_met,id%Elev,data1)
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading elevation in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             elevation = REAL(data1) * convert%Elev
          END IF
       ELSE ! If both PSurf and elevation aren't present, abort:
          CALL abort &
               ('Error finding PSurf or Elevation in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       END IF
       ! Note static pressure based on elevation in log file:
       WRITE(logn,*) 'PSurf not present in met file; values will be ', &
            'synthesised based on elevation and temperature.'
    END IF
    ! Look for CO2air (can be assumed to be static):- - - - - - - - - - -
    CALL find_metvarid(ncid_met, possible_varnames%CO2Names, id%CO2air, ok)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%CO2air = .TRUE. ! CO2air is present in met file
       ! Get CO2air units:
       ok = NF90_GET_ATT(ncid_met,id%CO2air,'units',metunits%CO2air)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding CO2air units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       IF(metunits%CO2air(1:3)/='ppm') THEN
          WRITE(*,*) metunits%CO2air
          CALL abort('Unknown units for CO2air'// &
               ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE ! CO2 not present
       exists%CO2air = .FALSE. ! CO2air is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Note this in log file:
       WRITE(logn,'(A33,A24,I4,A5)') ' CO2air not present in met file; ', &
            'values will be fixed at ',INT(fixedCO2),' ppmv'
    END IF
    ! Look for Snowf (could be part of Rainf variable):- - - - - - - - - -
    !IF (ncciy > 0 .AND. (globalMetfile%l_gpcc .OR. globalMetfile%l_gswp .OR. globalMetfile%l_ncar))Then
    IF (ncciy > 0 .AND. ( globalMetfile%l_gswp .OR. globalMetfile%l_ncar))Then
       ncid_met = ncid_snow
    END IF

    CALL find_metvarid(ncid_met, possible_varnames%SnowNames, id%Snowf, ok)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%Snowf = .TRUE. ! Snowf is present in met file
       ! Get Snowf units:
       ok = NF90_GET_ATT(ncid_met,id%Snowf,'units',metunits%Snowf)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding Snowf units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       ! Make sure Snowf units are the same as Rainf units:
       IF(metunits%Rainf/=metunits%Snowf) CALL abort &
            ('Please ensure Rainf and Snowf units are the same'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    ELSE
       exists%Snowf = .FALSE. ! Snowf is not present in met file
       !  all_met=.FALSE. not required; Snowf assumed to be in Rainf
       ! Note this in log file:
       WRITE(logn,*) 'Snowf not present in met file; ', &
            'Assumed to be contained in Rainf variable'
    END IF
    ! Look for LAI - - - - - - - - - - - - - - - - - - - - - - - - -
    CALL find_metvarid(ncid_met, possible_varnames%LAINames, id%LAI, ok)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%LAI = .TRUE. ! LAI is present in met file
       ! LAI will be read in which ever land grid is used
       ! Check dimension of LAI variable:
       ok=NF90_INQUIRE_VARIABLE(ncid_met,id%LAI, &
            ndims=lai_dims,dimids=laidimids)
       ! If any of LAI's dimensions are the time dimension
       IF(ANY(laidimids==timedimID(1))) THEN
          exists%LAI_T = .TRUE. ! i.e. time varying LAI
          WRITE(logn,*) 'LAI found in met file - time dependent;'
       ELSE
          exists%LAI_T = .FALSE. ! i.e. not time varying LAI
       END IF
       IF(ANY(laidimids==monthlydimID)) THEN
          exists%LAI_M = .TRUE. ! i.e. time varying LAI, but monthly only
          WRITE(logn,*) 'LAI found in met file - monthly values;'
       ELSE
          exists%LAI_M = .FALSE.
       END IF
       IF(ANY(laidimids==patchdimID)) THEN
          exists%LAI_P = .TRUE. ! i.e. patch varying LAI
          WRITE(logn,*) 'LAI found in met file - patch-specific values'
       ELSE
          exists%LAI_P = .FALSE. ! i.e. not patch varying LAI
       END IF
    ELSE
       exists%LAI = .FALSE. ! LAI is not present in met file
       ! Report to log file
       WRITE(logn,*) 'LAI not present in met file; ', &
            'Will use MODIS coarse grid monthly LAI'
    END IF
    ! If a spinup is to be performed:
    IF(spinup) THEN
       ! Look for avPrecip variable (time invariant - used for spinup):
       CALL find_metvarid(ncid_met, possible_varnames%APrecipNames, id%avPrecip, ok)
       IF(ok == NF90_NOERR) THEN ! If inquiry is okay and avPrecip exists
          ! Report to log file than modified spinup will be used:
          WRITE(logn,*) 'Spinup will use modified precip - avPrecip variable found'
          WRITE(logn,*) '  precip will be rescaled to match these values during spinup:'
          WRITE(*,*) 'Spinup will use modified precip - avPrecip variable found'
          WRITE(*,*) '  precip will be rescaled to match these values during spinup'
          ! Spinup will modify precip values:
          exists%avPrecip = .TRUE.
          ! Get avPrecip units:
          ok = NF90_GET_ATT(ncid_met,id%avPrecip,'units',metunits%avPrecip)
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error finding avPrecip units in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          IF(metunits%avPrecip(1:2)/='mm') CALL abort( &
               'Unknown avPrecip units in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ! Allocate space for avPrecip variable:
          ALLOCATE(avPrecip(mland))
          ! Get avPrecip from met file:
          IF(metGrid=='mask') THEN
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%avPrecip,data2, &
                     start=(/land_x(i),land_y(i)/),count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading avPrecip in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                avPrecip(i)=REAL(data2(1,1))
             END DO
          ELSE IF(metGrid=='land') THEN
             ! Allocate single preciaion temporary variable:
             ALLOCATE(temparray1(mland))
             ! Collect data from land only grid in netcdf file:
             ok= NF90_GET_VAR(ncid_met,id%avPrecip,temparray1)
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading avPrecip in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             ! Needed since r_1 will be double precision with -r8:
             avPrecip = REAL(temparray1)
             DEALLOCATE(temparray1)
          END IF
          ! Now find average precip from met data, and create rescaling
          ! factor for spinup:
          ALLOCATE(PrecipScale(mland))
          DO i = 1, mland
             IF(metGrid=='mask') THEN
                ! Allocate space for temporary precip variable:
                ALLOCATE(tempPrecip3(1,1,kend))
                ! Get all data for this grid cell:
                ok= NF90_GET_VAR(ncid_met,id%Rainf,tempPrecip3, &
                     start=(/land_x(i),land_y(i),1+koffset/),count=(/1,1,kend/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading Rainf in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Store total Rainf for this grid cell:
                PrecipTot = REAL(SUM(SUM(SUM(tempPrecip3,3),2))) &
                     * convert%Rainf
                ! Get snowfall data for this grid cell:
                IF(exists%Snowf) THEN
                   ok= NF90_GET_VAR(ncid_met,id%Snowf,tempPrecip3, &
                        start=(/land_x(i),land_y(i),1+koffset/),count=(/1,1,kend/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading Snowf in met data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                   ! Add total Snowf to this grid cell total:
                   PrecipTot = PrecipTot + &
                        (REAL(SUM(SUM(SUM(tempPrecip3,3),2))) &
                        * convert%Rainf)
                END IF
                DEALLOCATE(tempPrecip3)
             ELSE IF(metGrid=='land') THEN
                ! Allocate space for temporary precip variable:
                ALLOCATE(tempPrecip2(1,kend))
                ! Get rainfall data for this land grid cell:
                ok= NF90_GET_VAR(ncid_met,id%Rainf,tempPrecip2, &
                     start=(/i,1+koffset/),count=(/1,kend/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading Rainf in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Store total Rainf for this land grid cell:
                PrecipTot = REAL(SUM(SUM(tempPrecip2,2)))*convert%Rainf
                IF(exists%Snowf) THEN
                   ok= NF90_GET_VAR(ncid_met,id%Snowf,tempPrecip2, &
                        start=(/i,1+koffset/),count=(/1,kend/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading Snowf in met data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                   ! Add total Snowf to this land grid cell total:
                   PrecipTot = PrecipTot + (REAL(SUM(SUM(tempPrecip2,2))) &
                        * convert%Rainf)
                END IF
                DEALLOCATE(tempPrecip2)
             END IF
             ! Create rescaling factor for this grid cell to ensure spinup
             ! rainfall/snowfall is closer to average rainfall:
             ! First calculate annual average precip in met data:
             avPrecipInMet = PrecipTot/REAL(kend) * 3600.0/dels * 365 * 24
             PrecipScale(i) = avPrecipInMet/avPrecip(i)
             WRITE(logn,*) '  Site number:',i
             WRITE(logn,*) '  average precip quoted in avPrecip variable:', &
                  avPrecip(i)
             WRITE(logn,*) '  average precip in met data:',avPrecipInMet
          END DO ! over each land grid cell
          DEALLOCATE(avPrecip)
       ELSE ! avPrecip doesn't exist in met file
          ! Spinup will not modify precip values:
          exists%avPrecip = .FALSE.
          WRITE(logn,*) 'Spinup will repeat entire data set until states converge'
          WRITE(logn,*) '  (see below for convergence criteria);'
          WRITE(*,*) 'Spinup will repeat entire data set until states converge:'
       END IF
    END IF  ! if a spinup is to be performed

    ! Look for veg type - - - - - - - - - - - - - - - - -:
    CALL find_metvarid(ncid_met, possible_varnames%IVegNames, id%iveg, ok)
    IF(ok == NF90_NOERR) THEN ! If 'iveg' exists in the met file
       ! Note existence of at least one model parameter in the met file:
       exists%parameters = .TRUE.
       ! Allocate space for user-defined veg type variable:
       ALLOCATE(vegtype_metfile(mland,nmetpatches))
       IF(exists%patch)  ALLOCATE(vegpatch_metfile(mland,nmetpatches))

       ! Check dimension of veg type:
       ok=NF90_INQUIRE_VARIABLE(ncid_met,id%iveg,ndims=iveg_dims)
       IF(metGrid=='mask') THEN ! i.e. at least two spatial dimensions
          IF(iveg_dims==2) THEN ! no patch specific iveg information, just x,y
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%iveg,data2i, & ! get iveg data
                     start=(/land_x(i),land_y(i)/),count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading iveg in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Set all veg patches in grid cell to be this single type
                vegtype_metfile(i,:)=data2i(1,1)
             END DO
          ELSE IF(iveg_dims==3) THEN ! i.e. patch specific iveg information
             ! Patch-specific iveg variable MUST be accompanied by
             ! patchfrac variable with the same dimensions. So,
             ! Make sure that the patchfrac variable exists:
             CALL find_metvarid(ncid_met, possible_varnames%PFracNames, id%patchfrac, ok)
             IF(ok /= NF90_NOERR) CALL nc_abort & ! check read ok
                  (ok,'Patch-specific vegetation type (iveg) must be accompanied '// &
                  'by a patchfrac variable - this was not found in met data file '&
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             DO i = 1, mland
                ! Then, get the patch specific iveg data:
                ok= NF90_GET_VAR(ncid_met,id%iveg,vegtype_metfile(i,:), &
                     start=(/land_x(i),land_y(i),1/),count=(/1,1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort & ! check read ok
                     (ok,'Error reading iveg in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')

              IF(exists%patch) then
       
                !Anna: also read patch fractions
                ok= NF90_GET_VAR(ncid_met,id%patchfrac,vegpatch_metfile(i,:), &
                     start=(/land_x(i),land_y(i),1/),count=(/1,1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort & ! check read ok
                     (ok,'Error reading patchfrac in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
              END IF
             END DO
          END IF
       ELSE IF(metGrid=='land') THEN
          ! Collect data from land only grid in netcdf file:
          IF(iveg_dims==1) THEN ! i.e. no patch specific iveg information
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%iveg,data1i, &
                     start=(/i/),count=(/1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading iveg in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Set all veg patches in grid cell to be this single type
                vegtype_metfile(i,:) = data1i(1)
             END DO
          ELSE IF(iveg_dims==2) THEN ! i.e. patch specific iveg information
             ! Patch-specific iveg variable MUST be accompanied by
             ! patchfrac variable with same dimensions. So,
             ! Make sure that the patchfrac variable exists:
            CALL find_metvarid(ncid_met, possible_varnames%PFracNames, id%patchfrac, ok)
            IF(ok /= NF90_NOERR) CALL nc_abort & ! check read ok
                  (ok,'Patch-specific vegetation type (iveg) must be accompanied'// &
                  'by a patchfrac variable - this was not found in met data file '&
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')

              IF(exists%patch) then
                !Anna: also read patch fractions
                ok= NF90_GET_VAR(ncid_met,id%patchfrac,vegpatch_metfile(i,:), &
                     start=(/land_x(i),land_y(i),1/),count=(/1,1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort & ! check read ok
                     (ok,'Error reading patchfrac in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
              END IF

             DO i = 1, mland
                ! Then, get the patch specific iveg data:
                ok= NF90_GET_VAR(ncid_met, id%iveg, &
                     vegtype_metfile(i,:),&
                     start=(/i,1/), count=(/1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading iveg in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             END DO
          END IF
       END IF
    ELSE
       NULLIFY(vegtype_metfile)
       IF(exists%patch) NULLIFY(vegpatch_metfile)
       
    END IF

    ! Look for soil type:
    CALL find_metvarid(ncid_met, possible_varnames%ISoilNames, id%isoil, ok)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       ! Note existence of at least one model parameter in the met file:
       exists%parameters = .TRUE.
       ! Check dimension of soil type:
       ok=NF90_INQUIRE_VARIABLE(ncid_met,id%isoil,ndims=isoil_dims)
       ! Allocate space for user-defined soil type variable:
       ALLOCATE(soiltype_metfile(mland,nmetpatches))
       ! Get soil type from met file:
       IF(metGrid=='mask') THEN
          IF(isoil_dims==2) THEN ! i.e. no patch specific isoil information
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%isoil,data2i, &
                     start=(/land_x(i),land_y(i)/),count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading isoil in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Set all soil patches in grid cell to be this single type
                soiltype_metfile(i,:)=data2i(1,1)
             END DO
          ELSE IF(isoil_dims==3) THEN ! i.e. patch specific isoil information
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%isoil, &
                     soiltype_metfile(i,:), &
                     start=(/land_x(i),land_y(i),1/),count=(/1,1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading isoil in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             END DO
          END IF
       ELSE IF(metGrid=='land') THEN
          IF(isoil_dims==1) THEN ! i.e. no patch specific isoil information
             ! Collect data from land only grid in netcdf file:
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%isoil,data1i, &
                     start=(/i/),count=(/1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading isoil in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Set all veg patches in grid cell to be this single type
                soiltype_metfile(i,:) = data1i(1)
             END DO
          ELSE IF(isoil_dims==2) THEN ! i.e. patch specific isoil information
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met, id%isoil, &
                     soiltype_metfile(i,:), &
                     start=(/i,1/), count=(/1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading isoil in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             END DO
          END IF
       END IF
    ELSE
       NULLIFY(soiltype_metfile)
    END IF

    ! Report finding met variables to log file:
    IF(all_met) THEN
       WRITE(logn,*) 'Found all met variables in met file.'
    ELSE
       WRITE(logn,*) 'Found all ESSENTIAL met variables in met file,', &
            ' some synthesised (as above).'
    END IF

    !!=================^^ End met variables search^^=======================
  END SUBROUTINE open_met_file
  !==============================================================================
  !
  ! Name: get_met_data
  !
  ! Purpose: Fetches meteorological forcing data from the netcdf met forcing file
  !          for a single time step, including LAI if it exists.
  !          Note that currently met forcing is duplicated for every vegetated
  !          patch in a single gridcell.
  !
  ! CALLed from: cable_offline_driver
  !
  ! MODULEs used: cable_common_module
  !
  ! CALLs: abort
  !        nc_abort
  !        rh_sh
  !        sinbet
  !
  ! Input file: [SiteName].nc
  !
  !==============================================================================

  SUBROUTINE get_met_data(spinup,spinConv,met,soil,rad,                          &
       veg,kend,dels, TFRZ, ktau, kstart )
    ! Precision changes from REAL(4) to r_1 enable running with -r8


    ! Input arguments
    LOGICAL, INTENT(IN)                    ::                                   &
         spinup,         & ! are we performing a spinup?
         spinConv          ! has model spinup converged?
    TYPE(met_type),             INTENT(INOUT) :: met     ! meteorological data
    TYPE (soil_parameter_type),INTENT(IN)  :: soil
    TYPE (radiation_type),INTENT(IN)       :: rad
    TYPE(veg_parameter_type),INTENT(INOUT) :: veg ! LAI retrieved from file
    INTEGER, INTENT(IN)               :: ktau, &  ! timestep in loop including spinup
         kend, & ! total number of timesteps in run
         kstart  ! starting timestep
    REAL,INTENT(IN)                   :: dels ! time step size
    REAL, INTENT(IN) :: TFRZ

    ! Local variables
    REAL(KIND=4),DIMENSION(1,1,1)          :: data3 ! temp variable for netcdf reading
    REAL(KIND=4),DIMENSION(1,1,1,1)        :: data4 !  " " "
    REAL(KIND=4),DIMENSION(1,1)            :: data2 ! " "
    REAL(KIND=4),DIMENSION(1)              :: data1 ! " "
    INTEGER                           :: i,j ! do loop counter
    INTEGER                           :: ndims ! tempvar for number of dims in variable
    REAL(KIND=4),ALLOCATABLE,DIMENSION(:)       :: tmpDat1
    REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)     :: tmpDat2, tmpDat2x
    REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:)   :: tmpDat3, tmpDat3x
    REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:,:) :: tmpDat4, tmpDat4x

    DO i=1,mland ! over all land points/grid cells
       ! First set timing variables:
       ! All timing details below are initially written to the first patch
       ! of each gridcell, then dumped to all patches for the gridcell.
       IF(ktau==kstart) THEN ! initialise...
          SELECT CASE(time_coord)
          CASE('LOC')! i.e. use local time by default
             ! hour-of-day = starting hod
             met%hod(landpt(i)%cstart) = shod
             met%doy(landpt(i)%cstart) = sdoy
             met%moy(landpt(i)%cstart) = smoy
             met%year(landpt(i)%cstart) = syear
          CASE('GMT')! use GMT
             ! hour-of-day = starting hod + offset from GMT time:
             met%hod(landpt(i)%cstart) = shod + (longitude(i)/180.0)*12.0
             ! Note above that all met%* vars have dim mp,
             ! while longitude and latitude have dimension mland.
             met%doy(landpt(i)%cstart) = sdoy
             met%moy(landpt(i)%cstart) = smoy
             met%year(landpt(i)%cstart) = syear
          CASE DEFAULT
             CALL abort('Unknown time coordinate! ' &
                  //' (SUBROUTINE get_met_data)')
          END SELECT
       ELSE
          ! increment hour-of-day by time step size:
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) + dels/3600.0
       END IF
       !
       IF(met%hod(landpt(i)%cstart)<0.0) THEN ! may be -ve since longitude
          ! has range [-180,180]
          ! Reduce day-of-year by one and ammend hour-of-day:
          met%doy(landpt(i)%cstart) = met%doy(landpt(i)%cstart) - 1
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) + 24.0
          ! If a leap year AND we're using leap year timing:
          IF (is_leapyear(met%year(landpt(i)%cstart))) THEN
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(0) ! ie Dec previous year
                met%moy(landpt(i)%cstart) = 12
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) - 1
                met%doy(landpt(i)%cstart) = 365 ! prev year not leap year as this is
             CASE(31) ! Jan
                met%moy(landpt(i)%cstart) = 1
             CASE(60) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(91) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(121)
                met%moy(landpt(i)%cstart) = 4
             CASE(152)
                met%moy(landpt(i)%cstart) = 5
             CASE(182)
                met%moy(landpt(i)%cstart) = 6
             CASE(213)
                met%moy(landpt(i)%cstart) = 7
             CASE(244)
                met%moy(landpt(i)%cstart) = 8
             CASE(274)
                met%moy(landpt(i)%cstart) = 9
             CASE(305)
                met%moy(landpt(i)%cstart) = 10
             CASE(335)
                met%moy(landpt(i)%cstart) = 11
             END SELECT
          ELSE ! not a leap year or not using leap year timing
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(0) ! ie Dec previous year
                met%moy(landpt(i)%cstart) = 12
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) - 1
                ! If previous year is a leap year
                IF (is_leapyear(met%year(landpt(i)%cstart))) THEN
                   met%doy(landpt(i)%cstart) = 366
                ELSE
                   met%doy(landpt(i)%cstart) = 365
                END IF
             CASE(31) ! Jan
                met%moy(landpt(i)%cstart) = 1
             CASE(59) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(90)
                met%moy(landpt(i)%cstart) = 3
             CASE(120)
                met%moy(landpt(i)%cstart) = 4
             CASE(151)
                met%moy(landpt(i)%cstart) = 5
             CASE(181)
                met%moy(landpt(i)%cstart) = 6
             CASE(212)
                met%moy(landpt(i)%cstart) = 7
             CASE(243)
                met%moy(landpt(i)%cstart) = 8
             CASE(273)
                met%moy(landpt(i)%cstart) = 9
             CASE(304)
                met%moy(landpt(i)%cstart) = 10
             CASE(334)
                met%moy(landpt(i)%cstart) = 11
             END SELECT
          END IF ! if leap year or not
       ELSE IF(met%hod(landpt(i)%cstart)>=24.0) THEN
          ! increment or GMT adj has shifted day
          ! Adjust day-of-year and hour-of-day:
          met%doy(landpt(i)%cstart) = met%doy(landpt(i)%cstart) + 1
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) - 24.0
          ! If a leap year AND we're using leap year timing:
          !! vh_js !! use is_leapyear function here instead of multiple conditions
          IF (is_leapyear(met%year(landpt(i)%cstart))) THEN
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(32) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(61) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(92)
                met%moy(landpt(i)%cstart) = 4
             CASE(122)
                met%moy(landpt(i)%cstart) = 5
             CASE(153)
                met%moy(landpt(i)%cstart) = 6
             CASE(183)
                met%moy(landpt(i)%cstart) = 7
             CASE(214)
                met%moy(landpt(i)%cstart) = 8
             CASE(245)
                met%moy(landpt(i)%cstart) = 9
             CASE(275)
                met%moy(landpt(i)%cstart) = 10
             CASE(306)
                met%moy(landpt(i)%cstart) = 11
             CASE(336)
                met%moy(landpt(i)%cstart) = 12
             CASE(367)! end of year; increment
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) + 1
                met%moy(landpt(i)%cstart) = 1
                met%doy(landpt(i)%cstart) = 1
             END SELECT
             ! ELSE IF not leap year and Dec 31st, increment year
          ELSE
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(32) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(60) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(91)
                met%moy(landpt(i)%cstart) = 4
             CASE(121)
                met%moy(landpt(i)%cstart) = 5
             CASE(152)
                met%moy(landpt(i)%cstart) = 6
             CASE(182)
                met%moy(landpt(i)%cstart) = 7
             CASE(213)
                met%moy(landpt(i)%cstart) = 8
             CASE(244)
                met%moy(landpt(i)%cstart) = 9
             CASE(274)
                met%moy(landpt(i)%cstart) = 10
             CASE(305)
                met%moy(landpt(i)%cstart) = 11
             CASE(335)
                met%moy(landpt(i)%cstart) = 12
             CASE(366)! end of year; increment
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) + 1
                met%moy(landpt(i)%cstart) = 1
                met%doy(landpt(i)%cstart) = 1
             END SELECT
          END IF ! if leap year or not
       END IF ! if increment has pushed hod to a different day
       ! Now copy these values to all veg/soil patches in the current grid cell:
       met%hod(landpt(i)%cstart:landpt(i)%cend) = met%hod(landpt(i)%cstart)
       met%doy(landpt(i)%cstart:landpt(i)%cend) = met%doy(landpt(i)%cstart)
       met%moy(landpt(i)%cstart:landpt(i)%cend) = met%moy(landpt(i)%cstart)
       met%year(landpt(i)%cstart:landpt(i)%cend) = met%year(landpt(i)%cstart)
    ENDDO

    IF(metGrid=='mask') THEN
      ! N.B. not for GSWP runs, therefore only one met file here.
      ! Also, xdimsize and ydimsize are passed from io_variables.

       ALLOCATE(tmpDat2(xdimsize,ydimsize))
       ALLOCATE(tmpDat3(xdimsize,ydimsize,1))
       ALLOCATE(tmpDat4(xdimsize,ydimsize,1,1))
       ALLOCATE(tmpDat3x(xdimsize,ydimsize,nmetpatches))
       ALLOCATE(tmpDat4x(xdimsize,ydimsize,nmetpatches,1))

       ! Get SWdown data for mask grid:
       IF (cable_user%GSWP3) ncid_met=ncid_sw ! since GSWP3 multiple met files
       ok= NF90_GET_VAR(ncid_met,id%SWdown,tmpDat3, &
            start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading SWdown in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
       ! Assign value to met data variable (no units change required):
       !jhan:quick fix, use (1/2) dimension 1 here arbitrarily
       DO i=1,mland ! over all land points/grid cells
          met%fsd(landpt(i)%cstart:landpt(i)%cend,1) = &
               0.5 * REAL(tmpDat3(land_x(i),land_y(i),1))
          met%fsd(landpt(i)%cstart:landpt(i)%cend,2) = &
               0.5 * REAL(tmpDat3(land_x(i),land_y(i),1))
       ENDDO

       ! Get Tair data for mask grid:- - - - - - - - - - - - - - - - - -
       IF(cable_user%GSWP3) ncid_met = ncid_ta ! since GSWP3 multiple met files
       ! Find number of dimensions of Tair:
       ok = NF90_INQUIRE_VARIABLE(ncid_met,id%Tair,ndims=ndims)
       IF(ndims==3) THEN ! 3D var, either on grid or new ALMA format single site
         ok= NF90_GET_VAR(ncid_met,id%Tair,tmpDat3, &
              start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
         IF(ok /= NF90_NOERR) CALL nc_abort & ! check for error
              (ok,'Error reading Tair in met data file ' &
              //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
         ! Assign value to met data variable with units change:
         DO i=1,mland ! over all land points/grid cells
            met%tk(landpt(i)%cstart:landpt(i)%cend) = &
                 REAL(tmpDat3(land_x(i),land_y(i),1)) + convert%Tair
         ENDDO
       ELSE ! i.e. ndims==4, the older ALMA format for Tair
         ok= NF90_GET_VAR(ncid_met,id%Tair,tmpDat4, &
              start=(/1,1,1,ktau/),count=(/xdimsize,ydimsize,1,1/))
         IF(ok /= NF90_NOERR) CALL nc_abort &
              (ok,'Error reading Tair in met data file HERE' &
              //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
         ! Assign value to met data variable with units change:
         DO i=1,mland ! over all land points/grid cells
            met%tk(landpt(i)%cstart:landpt(i)%cend) = &
                 REAL(tmpDat4(land_x(i),land_y(i),1,1)) + convert%Tair
         ENDDO
       END IF

       ! Get PSurf data for mask grid:- - - - - - - - - - - - - - - - - -
       IF (cable_user%GSWP3) ncid_met = ncid_ps ! since GSWP3 multiple met files
       IF(exists%PSurf) THEN ! IF PSurf is in met file:
         ! Find number of dimensions of PSurf:
         ok = NF90_INQUIRE_VARIABLE(ncid_met,id%PSurf,ndims=ndims)
         IF(ndims==3) THEN ! 3D var, either grid or new ALMA format single site
           ok= NF90_GET_VAR(ncid_met,id%PSurf,tmpDat3, &
                start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
           IF(ok /= NF90_NOERR) CALL nc_abort &
                (ok,'Error reading PSurf in met data file ' &
                //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
           DO i=1,mland ! over all land points/grid cells
              met%pmb(landpt(i)%cstart:landpt(i)%cend) = &
                   REAL(tmpDat3(land_x(i),land_y(i),1)) * convert%PSurf
           ENDDO
         ELSE ! i.e. ndims==4, the older ALMA format for PSurf
           ok= NF90_GET_VAR(ncid_met,id%PSurf,tmpDat4, &
                start=(/1,1,1,ktau/),count=(/xdimsize,ydimsize,1/))
           IF(ok /= NF90_NOERR) CALL nc_abort &
                (ok,'Error reading PSurf in met data file ' &
                //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
           DO i=1,mland ! over all land points/grid cells
              met%pmb(landpt(i)%cstart:landpt(i)%cend) = &
                   REAL(tmpDat4(land_x(i),land_y(i),1,1)) * convert%PSurf
           ENDDO
         END IF
       ELSE ! PSurf must be fixed as a function of site elevation and T:
         DO i=1,mland ! over all land points/grid cells
            met%pmb(landpt(i)%cstart:landpt(i)%cend)=1013.25* &
                 (met%tk(landpt(i)%cstart)/(met%tk(landpt(i)%cstart) + 0.0065* &
                 elevation(i)))**(9.80665/287.04/0.0065)
         ENDDO
       END IF

       ! Get Qair data for mask grid: - - - - - - - - - - - - - - - - - -
       IF(cable_user%GSWP3) ncid_met = ncid_qa ! since GSWP3 multiple met files
       ! Find number of dimensions of Qair:
       ok = NF90_INQUIRE_VARIABLE(ncid_met,id%Qair,ndims=ndims)
       IF(ndims==3) THEN ! 3D var, either grid or new ALMA format single site
         ok= NF90_GET_VAR(ncid_met,id%Qair,tmpDat3, & ! read 3D Qair var
              start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
         IF(ok /= NF90_NOERR) CALL nc_abort & ! check for error
              (ok,'Error reading Qair in met data file ' &
              //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
         IF(convert%Qair==-999.0) THEN
           ! Convert relative value using only first veg/soil patch values
           ! (identical)
           DO i=1,mland ! over all land points/grid cells
             CALL rh_sh(REAL(tmpDat3(land_x(i),land_y(i),1)), &
                  met%tk(landpt(i)%cstart), &
                  met%pmb(landpt(i)%cstart),met%qv(landpt(i)%cstart))
                  met%qv(landpt(i)%cstart:landpt(i)%cend) = met%qv(landpt(i)%cstart)
           ENDDO
         ELSE
           DO i=1,mland ! over all land points/grid cells
              met%qv(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(tmpDat3(land_x(i),land_y(i),1))
           ENDDO
         END IF
       ELSE   ! i.e. ndims==4, the older ALMA format for Qair
         ok= NF90_GET_VAR(ncid_met,id%Qair,tmpDat4, &
            start=(/1,1,1,ktau/),count=(/xdimsize,ydimsize,1,1/))
         IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading Qair in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
         IF(convert%Qair==-999.0) THEN
           ! Convert relative value using only first veg/soil patch values
           ! (identical)
           DO i=1,mland ! over all land points/grid cells
              CALL rh_sh(REAL(tmpDat4(land_x(i),land_y(i),1,1)), &
                met%tk(landpt(i)%cstart), &
                met%pmb(landpt(i)%cstart),met%qv(landpt(i)%cstart))
                met%qv(landpt(i)%cstart:landpt(i)%cend) = met%qv(landpt(i)%cstart)
           ENDDO
         ELSE
           DO i=1,mland ! over all land points/grid cells
              met%qv(landpt(i)%cstart:landpt(i)%cend) = &
                    REAL(tmpDat4(land_x(i),land_y(i),1,1))
           ENDDO
         END IF
       END IF

       ! Get Wind data for mask grid: - - - - - - - - - - - - - - - - - -
       IF(cable_user%GSWP3) ncid_met = ncid_wd ! since GSWP3 multiple met files
       IF(exists%Wind) THEN ! Scalar Wind
         ! Find number of dimensions of Wind:
         ok = NF90_INQUIRE_VARIABLE(ncid_met,id%Wind,ndims=ndims)
         IF(ndims==3) THEN ! 3D var, either grid or new ALMA format single site
           ok= NF90_GET_VAR(ncid_met,id%Wind,tmpDat3, &
                start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
           IF(ok /= NF90_NOERR) CALL nc_abort &
                (ok,'Error reading Wind in met data file ' &
                //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
           ! Assign value to met data variable (no units change required):
           DO i=1,mland ! over all land points/grid cells
              met%ua(landpt(i)%cstart:landpt(i)%cend) = &
                   REAL(tmpDat3(land_x(i),land_y(i),1))
           ENDDO
         ELSE ! i.e. ndims==4, the older ALMA format for Wind
           ok= NF90_GET_VAR(ncid_met,id%Wind,tmpDat4, &
                start=(/1,1,1,ktau/),count=(/xdimsize,ydimsize,1,1/))
           IF(ok /= NF90_NOERR) CALL nc_abort &
                (ok,'Error reading Wind in met data file ' &
                //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
           ! Assign value to met data variable (no units change required):
           DO i=1,mland ! over all land points/grid cells
              met%ua(landpt(i)%cstart:landpt(i)%cend) = &
                   REAL(tmpDat4(land_x(i),land_y(i),1,1))
           ENDDO
         END IF ! 3 or 4D for 'Wind' variable
       ELSE ! Vector wind
         ! Find number of dimensions of Wind_N:
         ok = NF90_INQUIRE_VARIABLE(ncid_met,id%Wind,ndims=ndims)
         IF(ndims==3) THEN ! 3D var, either grid or new ALMA format single site
           ok= NF90_GET_VAR(ncid_met,id%Wind,tmpDat3, &
                start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
           IF(ok /= NF90_NOERR) CALL nc_abort &
                (ok,'Error reading Wind_N in met data file ' &
                //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
           ! Write part of wind variable to met%ua:
           DO i=1,mland ! over all land points/grid cells
              met%ua(landpt(i)%cstart) = REAL(tmpDat3(land_x(i),land_y(i),1))
           ENDDO
           ! Then fetch 3D Wind_E, and combine:
           ok= NF90_GET_VAR(ncid_met,id%Wind_E,tmpDat3, &
                start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
           IF(ok /= NF90_NOERR) CALL nc_abort &
                (ok,'Error reading Wind_E in met data file ' &
                //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
           ! Write final scalar Wind value:
           DO i=1,mland ! over all land points/grid cells
              met%ua(landpt(i)%cstart:landpt(i)%cend) = &
                   SQRT(met%ua(landpt(i)%cstart)**2 + &
                   REAL(tmpDat3(land_x(i),land_y(i),1))**2)
           ENDDO
         ELSE ! i.e. ndims==4, the older ALMA format for Wind_N and _E
           ! Get 4D Wind_N:
           ok= NF90_GET_VAR(ncid_met,id%Wind,tmpDat4, &
                start=(/1,1,1,ktau/),count=(/xdimsize,ydimsize,1,1/))
           IF(ok /= NF90_NOERR) CALL nc_abort &
                (ok,'Error reading Wind_N in met data file ' &
                //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
           ! Write part of wind variable to met%ua:
           DO i=1,mland ! over all land points/grid cells
              met%ua(landpt(i)%cstart) = REAL(tmpDat4(land_x(i),land_y(i),1,1))
           ENDDO
           ! Then fetch 4D Wind_E, and combine:
           ok= NF90_GET_VAR(ncid_met,id%Wind_E,tmpDat4, &
                start=(/1,1,1,ktau/),count=(/xdimsize,ydimsize,1,1/))
           IF(ok /= NF90_NOERR) CALL nc_abort &
                (ok,'Error reading Wind_E in met data file ' &
                //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
           ! Write final scalar Wind value:
           DO i=1,mland ! over all land points/grid cells
              met%ua(landpt(i)%cstart:landpt(i)%cend) = &
                   SQRT(met%ua(landpt(i)%cstart)**2 + &
                   REAL(tmpDat4(land_x(i),land_y(i),1,1))**2)
           ENDDO
         END IF ! 3 or 4D for 'Wind_N' and 'Wind_E' variables
       END IF ! scalar or vector wind - 'Wind' or 'Wind_N'/'Wind_E'

       ! Get Rainf and Snowf data for mask grid:- - - - - - - - - - - - -
       IF (cable_user%GSWP3) ncid_met = ncid_rain
       ok= NF90_GET_VAR(ncid_met,id%Rainf,tmpDat3, &
            start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading Rainf in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
       DO i=1,mland ! over all land points/grid cells
          met%precip(landpt(i)%cstart:landpt(i)%cend) = &
               REAL(tmpDat3(land_x(i),land_y(i),1)) ! store Rainf
       ENDDO
       IF(exists%Snowf) THEN
          IF (cable_user%GSWP3) ncid_met = ncid_snow
          ok= NF90_GET_VAR(ncid_met,id%Snowf,tmpDat3, &
               start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Snowf in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          ! store Snowf value (EK nov2007)
          DO i=1,mland ! over all land points/grid cells
             met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(tmpDat3(land_x(i),land_y(i),1))
          ENDDO
       ELSE
          met%precip_sn(:) = 0.0
       END IF
       ! combine Rainf and Snowf data
       met%precip(:) = met%precip(:) + met%precip_sn(:)
       ! Convert units:
       met%precip(:) = met%precip(:) * convert%Rainf
       met%precip_sn(:) = met%precip_sn(:) * convert%Rainf
       ! If we're performing a spinup, the spinup hasn't converged,
       ! and an avPrecip variable has been found, modify precip to
       ! ensure reasonable equilibration:
       IF(spinup.AND.(.NOT.spinConv).AND.exists%avPrecip) THEN
          ! Rescale precip to average rainfall for this site:
          DO i=1,mland ! over all land points/grid cells
             met%precip(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip(landpt(i)%cstart) / PrecipScale(i)
             ! Added for snow (EK nov2007)
             met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip_sn(landpt(i)%cstart) / PrecipScale(i)
          ENDDO
       END IF

       ! Get LWdown data for mask grid: - - - - - - - - - - - - - - - - -
       IF (cable_user%GSWP3) ncid_met = ncid_lw
       IF(exists%LWdown) THEN ! If LWdown exists in met file
          ok= NF90_GET_VAR(ncid_met,id%LWdown,tmpDat3, &
               start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading LWdown in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          DO i=1,mland ! over all land points/grid cells
             met%fld(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(tmpDat3(land_x(i),land_y(i),1))
          ENDDO
       ELSE ! Synthesise LWdown based on temperature
          ! Use Swinbank formula:
          met%fld(:) = 0.0000094*0.0000000567*(met%tk(:)**6.0)
       END IF

       ! Get CO2air data for mask grid:- - - - - - - - - - - - - - - - - -
       IF(exists%CO2air) THEN ! If CO2air exists in met file
          ok= NF90_GET_VAR(ncid_met,id%CO2air,tmpDat4, &
               start=(/1,1,1,ktau/),count=(/xdimsize,ydimsize,1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading CO2air in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          DO i=1,mland ! over all land points/grid cells
             met%ca(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(tmpDat4(land_x(i),land_y(i),1,1))/1000000.0
          ENDDO
       ELSE
          ! Fix CO2 air concentration:
          met%ca(:) = fixedCO2 /1000000.0
       END IF

       ! Get LAI, if it's present, for mask grid:- - - - - - - - - - - - -
       IF(exists%LAI) THEN ! If LAI exists in met file
          IF(exists%LAI_T) THEN ! i.e. time dependent LAI
             IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat4x, &
                     start=(/1,1,1,ktau/),count=(/xdimsize,ydimsize,nmetpatches,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met1 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   DO j=1,nmetpatches
                      veg%vlai(landpt(i)%cstart+j-1) = &
                           REAL(tmpDat4x(land_x(i),land_y(i),j,1))
                   ENDDO
                ENDDO
             ELSE ! i.e. patch independent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat3, &
                     start=(/1,1,ktau/),count=(/xdimsize,ydimsize,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met2 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = &
                        REAL(tmpDat3(land_x(i),land_y(i),1))
                ENDDO
             END IF
          ELSEIF(exists%LAI_M) THEN ! i.e. monthly LAI (BP apr08)
             IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat4x, &
                     start=(/1,1,1,met%moy/),count=(/xdimsize,ydimsize,nmetpatches,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met3 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   DO j=1,nmetpatches
                      veg%vlai(landpt(i)%cstart+j-1) = &
                           REAL(tmpDat4x(land_x(i),land_y(i),j,1))
                   ENDDO
                ENDDO
             ELSE ! i.e. patch independent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat3, &
                     start=(/1,1,met%moy/),count=(/xdimsize,ydimsize,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met4 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = &
                        REAL(tmpDat3(land_x(i),land_y(i),1))
                ENDDO
             END IF
          ELSE ! i.e. time independent LAI
             IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat3x, &
                     start=(/1,1,1/),count=(/xdimsize,ydimsize,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met5 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   DO j=1,nmetpatches
                      veg%vlai(landpt(i)%cstart+j-1) = &
                           REAL(tmpDat3x(land_x(i),land_y(i),j))
                   ENDDO
                ENDDO
             ELSE ! i.e. patch independent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat2, &
                     start=(/1,1/),count=(/xdimsize,ydimsize/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met6 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = &
                        REAL(tmpDat2(land_x(i),land_y(i)))
                ENDDO
             END IF
          END IF
       ELSE
          ! If not in met file, use default LAI value:
          DO i=1,mland ! over all land points/grid cells

             veg%vlai(landpt(i)%cstart:landpt(i)%cend) =  &
                  defaultLAI(landpt(i)%cstart:landpt(i)%cend,met%moy(landpt(i)%cstart))

          ENDDO
       END IF


       DEALLOCATE(tmpDat2,tmpDat3,tmpDat4,tmpDat3x,tmpDat4x)

    ELSE IF(metGrid=='land') THEN

       ! Collect data from land only grid in netcdf file:
       ALLOCATE(tmpDat1(mland))
       ALLOCATE(tmpDat2(mland,1))
       ALLOCATE(tmpDat2x(mland,nmetpatches))
       ALLOCATE(tmpDat3(mland,nmetpatches,1))

       ! Get SWdown data for land-only grid: - - - - - - - - - - - - -
       IF (ncciy > 0) ncid_met = ncid_sw
       ok= NF90_GET_VAR(ncid_met,id%SWdown,tmpDat2, &
            start=(/1,ktau/),count=(/mland,1/))
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading SWdown in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
       ! Assign value to met data variable (no units change required):
       DO i=1,mland ! over all land points/grid cells
          met%fsd(landpt(i)%cstart:landpt(i)%cend,1) = 0.5*REAL(tmpDat2(i,1))
          met%fsd(landpt(i)%cstart:landpt(i)%cend,2) = 0.5*REAL(tmpDat2(i,1))
       ENDDO

       ! Get Tair data for land-only grid:- - - - - - - - - - - - - - -
       IF (ncciy > 0) ncid_met = ncid_ta
       ok= NF90_GET_VAR(ncid_met,id%Tair,tmpDat2, &
            start=(/1,ktau/),count=(/mland,1/))
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading Tair in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
       ! Assign value to met data variable with units change:
       DO i=1,mland ! over all land points/grid cells
          met%tk(landpt(i)%cstart:landpt(i)%cend) = &
               REAL(tmpDat2(i,1)) + convert%Tair
       ENDDO

       ! Get PSurf data for land-only grid:- -- - - - - - - - - - - - - -
       IF (ncciy > 0) ncid_met = ncid_ps
       IF(exists%PSurf) THEN ! IF PSurf is in met file:
          IF ((ncciy == 1986) .AND. (ktau == 2184)) THEN
             !hzz to fix the problem of ps data on time step 2184
             ok= NF90_GET_VAR(ncid_met,id%PSurf,tmpDat2, &
                  start=(/1,2176/),count=(/mland,1/)) ! fixing bug in GSWP ps data
          ELSE
             ok= NF90_GET_VAR(ncid_met,id%PSurf,tmpDat2, &
                  start=(/1,ktau/),count=(/mland,1/))
          ENDIF
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading PSurf in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          DO i=1,mland ! over all land points/grid cells
             met%pmb(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(tmpDat2(i,1)) * convert%PSurf
          ENDDO
       ELSE ! PSurf must be fixed as a function of site elevation and T:
          DO i=1,mland ! over all land points/grid cells
             met%pmb(landpt(i)%cstart:landpt(i)%cend) = 1013.25 &
                  *(met%tk(landpt(i)%cstart)/(met%tk(landpt(i)%cstart) &
                  + 0.0065*elevation(i)))**(9.80665/287.04/0.0065)
          ENDDO
       END IF

       ! Get Qair data for land-only grid:- - - - - - - - - - - - - - - -
       IF (ncciy > 0) ncid_met = ncid_qa
       ok= NF90_GET_VAR(ncid_met,id%Qair,tmpDat2, &
            start=(/1,ktau/),count=(/mland,1/))
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading Qair in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
       IF(convert%Qair==-999.0) THEN
          DO i=1,mland ! over all land points/grid cells
             CALL rh_sh(REAL(tmpDat2(i,1)), met%tk(landpt(i)%cstart), &
                  met%pmb(landpt(i)%cstart),met%qv(landpt(i)%cstart))
             met%qv(landpt(i)%cstart:landpt(i)%cend)=met%qv(landpt(i)%cstart)
          ENDDO
       ELSE
          DO i=1,mland ! over all land points/grid cells
             met%qv(landpt(i)%cstart:landpt(i)%cend) = REAL(tmpDat2(i,1))
          ENDDO
       END IF

       ! Get Wind data for land-only grid: - - - - - - - - - - - - - - - -
       IF (ncciy > 0) ncid_met = ncid_wd
       IF(exists%Wind) THEN ! Scalar Wind
          ok= NF90_GET_VAR(ncid_met,id%Wind,tmpDat2, &
               start=(/1,ktau/),count=(/mland,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Wind in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable (no units change required):
          DO i=1,mland ! over all land points/grid cells
             met%ua(landpt(i)%cstart:landpt(i)%cend) = REAL(tmpDat2(i,1))
          ENDDO
       ELSE ! Vector wind
          ! Get Wind_N:
          ok= NF90_GET_VAR(ncid_met,id%Wind,tmpDat2, &
               start=(/1,ktau/),count=(/mland,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Wind_N in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          ! write part of the wind variable
          met%ua(landpt(:)%cstart) = REAL(tmpDat2(:,1))
          ok= NF90_GET_VAR(ncid_met,id%Wind_E,tmpDat2, &
               start=(/1,ktau/),count=(/mland,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Wind_E in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          ! Write final scalar Wind value:
          DO i=1,mland ! over all land points/grid cells
             met%ua(landpt(i)%cstart:landpt(i)%cend) = &
                  SQRT(met%ua(landpt(i)%cstart)**2 + REAL(tmpDat2(i,1))**2)
          ENDDO
       END IF

       ! Get Rainf and Snowf data for land-only grid: - - - - - - - - - - -
       IF (ncciy > 0) ncid_met = ncid_rain
       ok= NF90_GET_VAR(ncid_met,id%Rainf,tmpDat2, &
            start=(/1,ktau/),count=(/mland,1/))
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading Rainf in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
       DO i=1,mland ! over all land points/grid cells
          met%precip(landpt(i)%cstart:landpt(i)%cend) = REAL(tmpDat2(i,1))
       ENDDO

       IF (ncciy > 0) ncid_met = ncid_snow
       IF(exists%Snowf) THEN
          ok= NF90_GET_VAR(ncid_met,id%Snowf,tmpDat2, &
               start=(/1,ktau/),count=(/mland,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Snowf in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          DO i=1,mland ! over all land points/grid cells
             met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(tmpDat2(i,1))
          ENDDO
       ELSE
          met%precip_sn(:) = 0.0
       END IF

       ! combine Rainf and Snowf data
       met%precip(:) = met%precip(:) + met%precip_sn(:)
       ! Convert units:
       met%precip(:) = met%precip(:) * convert%Rainf
       met%precip_sn(:) = met%precip_sn(:) * convert%Rainf
       ! If we're performing a spinup, the spinup hasn't converged,
       ! and an avPrecip variable has been found, modify precip to
       ! ensure reasonable equilibration:
       IF(spinup.AND.(.NOT.spinConv).AND.exists%avPrecip) THEN
          ! Rescale precip to average rainfall for this site:
          DO i=1,mland ! over all land points/grid cells
             met%precip(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip(landpt(i)%cstart) / PrecipScale(i)
             met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip_sn(landpt(i)%cstart) / PrecipScale(i)
          ENDDO
       END IF

       ! Get LWdown data for land-only grid: - - - - - - - - - - - - - -
       IF (ncciy > 0) ncid_met = ncid_lw
       IF(exists%LWdown) THEN ! If LWdown exists in met file
          ok= NF90_GET_VAR(ncid_met,id%LWdown,tmpDat2, &
               start=(/1,ktau/),count=(/mland,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading LWdown in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          DO i=1,mland ! over all land points/grid cells
             met%fld(landpt(i)%cstart:landpt(i)%cend) = REAL(tmpDat2(i,1))
          ENDDO
       ELSE ! Synthesise LWdown based on temperature
          ! Use Swinbank formula:
          met%fld(:) = 0.0000094*0.0000000567*(met%tk(:)**6.0)
       END IF

       ! Get CO2air data for land-only grid:- - - - - - - - - - - - - -
       IF(exists%CO2air) THEN ! If CO2air exists in met file
          ok= NF90_GET_VAR(ncid_met,id%CO2air,tmpDat2, &
               start=(/1,ktau/),count=(/mland,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading CO2air in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          DO i=1,mland ! over all land points/grid cells
             met%ca(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(tmpDat2(i,1)) / 1000000.0
          ENDDO
       ELSE
          met%ca(:) = fixedCO2 /1000000.0
       END IF

       ! Get LAI data, if it exists, for land-only grid:- - - - - - - - -
       IF(exists%LAI) THEN ! If LAI exists in met file
          IF(exists%LAI_T) THEN ! i.e. time dependent LAI
             IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat3, &
                     start=(/1,1,ktau/),count=(/mland,nmetpatches,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met7 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   IF ( (landpt(i)%cend - landpt(i)%cstart + 1) < nmetpatches) THEN
                      PRINT *, 'not enough patches at land point ', i
                      STOP
                   END IF
                   DO j=1, nmetpatches
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(tmpDat3(i,j,1))
                   ENDDO
                ENDDO
             ELSE ! i.e. patch independent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat2, &
                     start=(/1,ktau/),count=(/mland,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met8 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = &
                        REAL(tmpDat2(i,1))
                ENDDO
             END IF
          ELSEIF(exists%LAI_M) THEN ! i.e. monthly LAI (BP apr08)
             IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat3, &
                     start=(/1,1,met%moy/),count=(/mland,nmetpatches,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met9 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   DO j=1, nmetpatches
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(tmpDat3(i,j,1))
                   ENDDO
                ENDDO
             ELSE ! i.e. patch independent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat2, &
                     start=(/1,met%moy/),count=(/mland,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met10 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = &
                        REAL(tmpDat2(i,1))
                ENDDO
             END IF
          ELSE ! LAI time independent
             IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat2x, &
                     start=(/1,1/),count=(/mland,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met11 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   DO j=1, nmetpatches
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(tmpDat2x(i,j))
                   ENDDO
                ENDDO
             ELSE ! i.e. patch independent LAI
                ok= NF90_GET_VAR(ncid_met,id%LAI,tmpDat1, &
                     start=(/1/),count=(/mland/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading LAI in met12 data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                DO i=1,mland ! over all land points/grid cells
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = REAL(tmpDat1(i))
                ENDDO
             END IF
          END IF
       ELSE
          ! If not in met file, use default LAI value:
          DO i=1,mland ! over all land points/grid cells
             !! vh_js !! corrected indices of defaultLAI
             veg%vlai(landpt(i)%cstart:landpt(i)%cend) =  &
                  defaultLAI(landpt(i)%cstart:landpt(i)%cend,met%moy(landpt(i)%cstart))



          ENDDO

       END IF
       DEALLOCATE(tmpDat1, tmpDat2, tmpDat3, tmpDat2x)

    ELSE
       CALL abort('Unrecognised grid type')
    END IF ! grid type

    IF ((.NOT. exists%Snowf) .OR. ALL(met%precip_sn == 0.0)) THEN ! honour snowf input
       DO i=1,mland ! over all land points/grid cells
          ! Set solid precip based on temp
          met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = 0.0 ! (EK nov2007)
          IF( met%tk(landpt(i)%cstart) <= tfrz ) &
               met%precip_sn(landpt(i)%cstart:landpt(i)%cend) &
               = met%precip(landpt(i)%cstart) ! (EK nov2007)
       END DO ! 1, mland over all land grid points
    ENDIF

    ! Set cosine of zenith angle (provided by GCM when online):
    met%coszen = sinbet(met%doy, rad%latitude, met%hod)
    ! initialise within canopy air temp
    met%tvair = met%tk
    met%tvrad = met%tk
    IF(check%ranges) THEN
       ! Check ranges are okay:
       !jhan:quick fix, use dimension 1 here arbitrarily
       IF(ANY(met%fsd(:,1)<ranges%SWdown(1)).OR.ANY(met%fsd(:,1)>ranges%SWdown(2))) &
            CALL abort('SWdown out of specified ranges!')
       IF(ANY(met%fsd(:,2)<ranges%SWdown(1)).OR.ANY(met%fsd(:,2)>ranges%SWdown(2))) &
            CALL abort('SWdown out of specified ranges!')
       IF(ANY(met%fld<ranges%LWdown(1)).OR.ANY(met%fld>ranges%LWdown(2))) &
            CALL abort('LWdown out of specified ranges!')
       IF(ANY(met%qv<ranges%Qair(1)).OR.ANY(met%qv>ranges%Qair(2))) &
            CALL abort('Qair out of specified ranges!')
       IF(ANY(met%precip<ranges%Rainf(1)).OR.ANY(met%precip>ranges%Rainf(2))) THEN
          CALL abort('Rainf out of specified ranges!')
       ENDIF
       IF(ANY(met%ua<ranges%Wind(1)).OR.ANY(met%ua>ranges%Wind(2))) &
            CALL abort('Wind out of specified ranges!')
       IF(ANY(met%tk<ranges%Tair(1)).OR.ANY(met%tk>ranges%Tair(2))) &
            CALL abort('Tair out of specified ranges!')
       IF(ANY(met%pmb<ranges%PSurf(1)).OR.ANY(met%pmb>ranges%PSurf(2))) THEN
          WRITE(*,*) "min, max Psurf", MINVAL(met%pmb), MAXVAL(met%pmb),ranges%Psurf(1), ranges%Psurf(2)
          CALL abort('PSurf out of specified ranges!')
       ENDIF
    END IF

  END SUBROUTINE get_met_data
  !==============================================================================
  !
  ! Name: close_met_file
  !
  ! Purpose: Close the file with the meteorological data
  !
  ! CALLed from: cable_offline_driver
  !
  ! CALLs: nc_abort
  !
  ! Input file: [SiteName].nc
  !
  !==============================================================================

  SUBROUTINE close_met_file

    ok=NF90_CLOSE(ncid_met)
    IF(ok /= NF90_NOERR) CALL nc_abort (ok,'Error closing met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE close_met_file)')
    ! Clear lat_all and lon_all variables

  END SUBROUTINE close_met_file

  !==============================================================================
  !
  ! Name: load_parameters
  !
  ! Purpose: Checks where parameters and initialisations should be loaded from.
  !          If they can be found in either the met file or restart file, they
  !          will load from there, with the met file taking precedence. Otherwise,
  !          they'll be chosen from a coarse global grid of veg and soil types,
  !          based on the lat/lon coordinates.
  !
  ! CALLed from: cable_offline_driver
  !
  ! CALLs: get_default_params
  !        allocate_cable_vars
  !        alloc_casavariable
  !        alloc_phenvariable
  !        write_default_params
  !        write_cnp_params
  !        casa_readbiome
  !        casa_readphen
  !        casa_init
  !        abort
  !        get_restart_data
  !        get_parameters_met
  !        derived_parameters
  !        check_parameter_values
  !        report_parameters
  !
  ! Input file: [restart].nc
  !
  !==============================================================================

  SUBROUTINE load_parameters(met,air,ssnow,veg,climate,bgc,soil,canopy,rough,rad,        &
       sum_flux,bal,logn,vegparmnew,casabiome,casapool,    &
       casaflux,sum_casapool, sum_casaflux,casamet,casabal,phen,POP,spinup,EMSOIL, &
       TFRZ, LUC_EXPT, POPLUC)
    ! Input variables not listed:
    !   filename%type  - via cable_IO_vars_module
    !   exists%type    - via cable_IO_vars_module
    !   smoy           - via cable_IO_vars_module
    ! Output variables not listed:
    !   (determined here or from sub get_default_params <- countPatch)
    !   landpt%type    - via cable_IO_vars_module (nap,cstart,cend,ilon,ilat)
    !   max_vegpatches - via cable_IO_vars_module
    !! vh_js !!
    USE POPmodule, ONLY: POP_INIT
    USE POPLUC_module, ONLY: POPLUC_INIT
    USE CABLE_LUC_EXPT, ONLY: LUC_EXPT_TYPE

    IMPLICIT NONE

    ! Input arguments
    TYPE (met_type), INTENT(INOUT)          :: met
    TYPE (air_type), INTENT(INOUT)          :: air
    TYPE (soil_snow_type), INTENT(OUT)      :: ssnow
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg
    TYPE (climate_type), INTENT(INOUT)          :: climate
    TYPE (bgc_pool_type), INTENT(OUT)       :: bgc
    TYPE (soil_parameter_type), INTENT(OUT) :: soil
    TYPE (canopy_type), INTENT(OUT)         :: canopy
    TYPE (roughness_type), INTENT(OUT)      :: rough
    TYPE (radiation_type),INTENT(OUT)       :: rad
    TYPE (sum_flux_type), INTENT(OUT)       :: sum_flux
    TYPE (balances_type), INTENT(OUT)       :: bal
    TYPE (casa_biome)  , INTENT(OUT)        :: casabiome
    TYPE (casa_pool)   , INTENT(OUT)        :: casapool
    TYPE (casa_flux)   , INTENT(OUT)        :: casaflux
    TYPE (casa_pool)   , INTENT(OUT)        :: sum_casapool
    TYPE (casa_flux)   , INTENT(OUT)        :: sum_casaflux
    TYPE (casa_met)    , INTENT(OUT)        :: casamet
    TYPE (casa_balance), INTENT(OUT)        :: casabal
    TYPE(phen_variable), INTENT(OUT)        :: phen
    TYPE( POP_TYPE ), INTENT(INOUT)         :: POP
    TYPE( POPLUC_TYPE ), INTENT(INOUT)         :: POPLUC
    TYPE (LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT
    INTEGER,INTENT(IN)                      :: logn     ! log file unit number
    LOGICAL,INTENT(IN)                      :: &
         vegparmnew, &  ! are we using the new format?
                                !! vh_js !!
         spinup         ! for POP (initialize pop)
    REAL, INTENT(IN) :: TFRZ, EMSOIL

    ! Local variables
    REAL,POINTER,DIMENSION(:)          :: pfractmp ! temp store of patch fraction
    LOGICAL                                 :: completeSet ! was a complete parameter set found?
    LOGICAL                            :: EXRST = .FALSE. ! does a RunIden restart file exist?
    INTEGER                            ::                                  &
         mp_restart,        & ! total number of patches in restart file
         mpID,              &
         napID,             &
         i , j                   ! do loop variables
    !! vh_js !!
    ! CHARACTER :: frst_in*100, CYEAR*4
    CHARACTER :: frst_in*200, CYEAR*4

    INTEGER   :: IOS
    CHARACTER :: TACC*20
    INTEGER,DIMENSION(:), ALLOCATABLE :: ALLVEG
    !! vh_js !!
    INTEGER :: mp_POP
    INTEGER, DIMENSION(:), ALLOCATABLE :: Iwood

    ! Allocate spatial heterogeneity variables:
    ALLOCATE(landpt(mland))

    WRITE(logn,*) '-------------------------------------------------------'
    WRITE(logn,*) 'Looking for parameters and initial states....'
    WRITE(logn,*) ' Loading initialisations from default grid.'

    ! Parameter values and some grid info are read in.
    ! They will be overwritten by values from the restart file, if present.
    ! Those variables found in the met file will again overwrite existing ones.

    CALL get_default_params(logn,vegparmnew,LUC_EXPT)
    CALL allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssnow, &
         sum_flux,veg,mp)
    WRITE(logn,*) ' CABLE variables allocated with ', mp, ' patch(es).'

    IF (icycle > 0 .OR. CABLE_USER%CASA_DUMP_WRITE ) &
         CALL alloc_casavariable(casabiome,casapool,casaflux, &
         casamet,casabal,mp)
    !mpdiff
    CALL alloc_sum_casavariable(sum_casapool,sum_casaflux,mp)
    IF (icycle > 0) THEN
       CALL alloc_phenvariable(phen,mp)
    ENDIF

    ! Write parameter values to CABLE's parameter variables:
    CALL write_default_params(met,air,ssnow,veg,bgc,soil,canopy,rough, &
         rad,logn,vegparmnew,smoy, TFRZ, LUC_EXPT)



    ! Zero out lai where there is no vegetation acc. to veg. index
    WHERE ( veg%iveg(:) .GE. 14 ) veg%vlai = 0.

    IF (icycle > 0) THEN
       CALL write_cnp_params(veg,casaflux,casamet)
       CALL casa_readbiome(veg,soil,casabiome,casapool,casaflux, &
            casamet,phen)
       IF (cable_user%PHENOLOGY_SWITCH.EQ.'MODIS') CALL casa_readphen(veg,casamet,phen)

       CALL casa_init(casabiome,casamet,casaflux,casapool,casabal,veg,phen)
       !! vh_js !!
       IF ( CABLE_USER%CALL_POP ) THEN
          ! evaluate mp_POP and POP_array
          mp_POP = COUNT(casamet%iveg2==forest)+COUNT(casamet%iveg2==shrub)

          ALLOCATE(Iwood(mp_POP))
          j = 1
          DO i=1,mp
             IF (casamet%iveg2(i)==forest .OR. casamet%iveg2(i)==shrub) THEN
                Iwood(j) = i
                j = j+1
             ENDIF
          ENDDO

          CALL POP_init( POP, veg%disturbance_interval(Iwood,:), mp_POP, Iwood )
          IF ( .NOT. (spinup .OR. CABLE_USER%POP_fromZero )) &
               CALL POP_IO( POP, casamet, cable_user%YearStart, "READ_rst " , .TRUE.)

       ENDIF

       IF (CABLE_USER%POPLUC) THEN

          ! initialise POPLUC structure and params
          !zero biomass in secondary forest tiles
          ! read POP_LUC restart file here
          ! set POP%LU here for secondary tiles if cable_user%POPLUC_RunType is not 'static'
          CALL POPLUC_init(POPLUC,LUC_EXPT, casapool, casaflux, casabiome, veg, POP, mland)

       ENDIF

    ENDIF

    ! removed get_default_inits and get_default_lai as they are already done
    ! in write_default_params
    !    ! Load default initialisations from Mk3L climatology:
    !    CALL get_default_inits(met,soil,ssnow,canopy,logn)
    !
    !    ! load default LAI values from global data:
    !    CALL get_default_lai

    ! Look for explicit restart file (which will have parameters):
    IF ( TRIM(filename%restart_in) .EQ. '' ) filename%restart_in = './'
    frst_in = filename%restart_in
    ok = NF90_OPEN(TRIM(frst_in),NF90_NOWRITE,ncid_rin)
    IF ( ok == NF90_NOERR ) EXRST = .TRUE.

    ! If not an explicit rstfile, search for RunIden_YEAR...nc
    ! use (filename%restart_in) as path
    IF ( .NOT. EXRST .AND. CABLE_USER%YEARSTART .GT. 0 ) THEN
       WRITE( CYEAR,FMT="(I4)" ) CurYear
       frst_in = TRIM(filename%restart_in)//'/'//TRIM(cable_user%RunIden)//&
            '_'//CYEAR//'_cable_rst.nc'
       INQUIRE( FILE=TRIM( frst_in ), EXIST=EXRST )
    ENDIF

    IF ( EXRST ) THEN
       ok = NF90_OPEN(TRIM(frst_in),NF90_NOWRITE,ncid_rin) ! open restart file
       IF (ok /= NF90_NOERR) CALL HANDLE_ERR(ok)
       ! Any restart file exists, parameters and init will be loaded from it.
       WRITE(logn,*) ' Overwriting initialisations with values in ', &
            'restart file: ', TRIM(frst_in)
       WRITE(*,*)    ' Overwriting initialisations with values in ', &
            'restart file: ', TRIM(frst_in)

       ! Check total number of patches in restart file:
       ok = NF90_INQ_DIMID(ncid_rin,'mp',mpID)
       IF(ok /= NF90_NOERR) THEN
          ok = NF90_INQ_DIMID(ncid_rin,'mp_patch',mpID)
          IF(ok /= NF90_NOERR)  CALL nc_abort &
               (ok,'Error finding mp or mp_patch dimension in restart file ' &
               //TRIM(frst_in)//' (SUBROUTINE load_parameters) ' &
               //'Recommend running without restart file.')
       END IF
       ok = NF90_INQUIRE_DIMENSION(ncid_rin,mpID,len=mp_restart)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding total number of patches in restart file ' &
            //TRIM(frst_in)//' (SUBROUTINE load_parameters) ' &
            //'Recommend running without restart file.')
       ! Check that mp_restart = mp from default/met values
       IF(mp_restart /= mp) CALL abort('Number of patches in '// &
            'restart file '//TRIM(frst_in)//' does not equal '// &
            'to number in default/met file settings. (SUB load_parameters) ' &
            //'Recommend running without restart file.')

       ! Load initialisations and parameters from restart file:
       CALL get_restart_data(logn,ssnow,canopy,rough,bgc,bal,veg, &
            soil,rad,vegparmnew, EMSOIL )

    ELSE
       ! With no restart file, use default parameters already loaded
       WRITE(logn,*) ' Could neither find restart file ', TRIM(filename%restart_in)
       WRITE(logn,*) ' nor ', TRIM(frst_in)
       WRITE(logn,*) ' Pre-loaded default initialisations are used.'
       WRITE(*,*)    ' Could neither find restart file ', TRIM(filename%restart_in)
       WRITE(*,*)    ' nor ', TRIM(frst_in)
       WRITE(*,*)    ' Pre-loaded default initialisations are used.'
    END IF ! if restart file exists

    ! Overwrite default values by those available in met file:
    CALL get_parameters_met(soil,veg,bgc,rough,completeSet)

    ! Results of looking for parameters in the met file:
    WRITE(logn,*)
    IF(exists%parameters.AND.completeSet) THEN
       ! All pars were found in met file:
       WRITE(logn,*) ' Loaded all parameters from met input file: ', &
            TRIM(filename%met)
    ELSE IF(exists%parameters.AND..NOT.completeSet) THEN
       ! Only some pars were found in met file:
       WRITE(logn,*) ' Loaded some parameters from met input file: ', &
            TRIM(filename%met), & ! write to log file
            ' the rest are default values'
       WRITE(*,*)    ' Loaded some parameters from met input file: ', &
            TRIM(filename%met), & ! write to screen
            ' the rest are default values - check log file'
    ELSE
       ! No parameters were found in met file:
       WRITE(logn,*) ' Met file has no params; all parameters remain as default.'
    END IF
    WRITE(logn,*)

    ! Construct derived parameters and zero initialisations, regardless
    ! of where parameters and other initialisations have loaded from:
    CALL derived_parameters(soil,sum_flux,bal,ssnow,veg,rough)

    ! Check for basic inconsistencies in parameter values:
    CALL check_parameter_values(soil,veg,ssnow)

    ! Write per-site parameter values to log file if requested:
    CALL report_parameters(logn,soil,veg,bgc,rough,ssnow,canopy, &
         casamet,casapool,casaflux,phen,vegparmnew,verbose)


  END SUBROUTINE load_parameters


  !==============================================================================
  !
  ! Name: get_parameters_met
  !
  ! Purpose: This subroutine looks for parameters in the met file, and loads
  !          those that are found.
  !
  ! CALLed from: load_parameters
  !              old_load_parameters
  !
  ! CALLs: readpar
  !
  ! Input file: [SiteName].nc
  !
  !==============================================================================

  SUBROUTINE get_parameters_met(soil,veg,bgc,rough,completeSet)

    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE (bgc_pool_type), INTENT(INOUT)       :: bgc
    TYPE (roughness_type), INTENT(INOUT)      :: rough
    LOGICAL, INTENT(OUT)                      :: completeSet ! were all pars found?

    ! Local variables
    INTEGER                              :: parID ! parameter's netcdf ID

    ! removed the following section because already in IGBP types (BP apr08)
    !    ! First, if user defined surface type ratios are present in the
    !    ! met file then use them:
    !    IF(ASSOCIATED(vegfrac_user)) THEN
    !       DO i=1,mland
    !          ! Overwrite landpt(i)%*%frac, which will be set either by restart
    !          ! or default values:
    !          landpt(i)%veg%frac = vegfrac_user(i)
    !          landpt(i)%urban%frac = urbanfrac_user(i)
    !          landpt(i)%lake%frac = lakefrac_user(i)
    !          landpt(i)%ice%frac = icefrac_user(i)
    !       END DO
    !    END IF

    completeSet=.TRUE. ! initialise (assume all param will load from met file)

    ! Get parameter values:
    ! Arguments: netcdf file ID; parameter name; complete set check;
    !   parameter value; filename for error messages; number of veg/soil patches
    !   in met file; switch to indicate size of dimensions of the parameter.
    ! ! Use 'defd' for single dim double precision.
    ! veg and soil types already obtained in sub open_met_file
    !    CALL readpar(ncid_met,'iveg',completeSet,veg%iveg,filename%met, &
    !         nmetpatches,'def')
    CALL readpar(ncid_met,'patchfrac',completeSet,patch(:)%frac,filename%met,   &
         nmetpatches,'def')
    !    CALL readpar(ncid_met,'isoil',completeSet,soil%isoilm,filename%met, &
    !         nmetpatches,'def')
    CALL readpar(ncid_met,'clay',completeSet,soil%clay,filename%met,            &
         nmetpatches,'def')
    CALL readpar(ncid_met,'sand',completeSet,soil%sand,filename%met,            &
         nmetpatches,'def')
    CALL readpar(ncid_met,'silt',completeSet,soil%silt,filename%met,            &
         nmetpatches,'def')
    CALL readpar(ncid_met,'ssat',completeSet,soil%ssat,filename%met,            &
         nmetpatches,'def')
    CALL readpar(ncid_met,'sfc',completeSet,soil%sfc,filename%met,              &
         nmetpatches,'def')
    CALL readpar(ncid_met,'swilt',completeSet,soil%swilt,filename%met,          &
         nmetpatches,'def')
    CALL readpar(ncid_met,'bch',completeSet,soil%bch,filename%met,              &
         nmetpatches,'def')
    CALL readpar(ncid_met,'hyds',completeSet,soil%hyds,filename%met,            &
         nmetpatches,'def')
    CALL readpar(ncid_met,'sucs',completeSet,soil%sucs,filename%met,            &
         nmetpatches,'def')
    CALL readpar(ncid_met,'css',completeSet,soil%css,filename%met,              &
         nmetpatches,'def')
    CALL readpar(ncid_met,'rhosoil',completeSet,soil%rhosoil,filename%met,      &
         nmetpatches,'def')
    CALL readpar(ncid_met,'rs20',completeSet,veg%rs20,filename%met,             &
         nmetpatches,'def')
    CALL readpar(ncid_met,'albsoil',completeSet,soil%albsoil,filename%met,      &
         nmetpatches,'nrb')
    CALL readpar(ncid_met,'froot',completeSet,veg%froot,filename%met,           &
         nmetpatches,'ms')
    CALL readpar(ncid_met,'hc',completeSet,veg%hc,filename%met,                 &
         nmetpatches,'def')
    CALL readpar(ncid_met,'canst1',completeSet,veg%canst1,filename%met,         &
         nmetpatches,'def')
    CALL readpar(ncid_met,'dleaf',completeSet,veg%dleaf,filename%met,           &
         nmetpatches,'def')
    CALL readpar(ncid_met,'frac4',completeSet,veg%frac4,filename%met,           &
         nmetpatches,'def')
    CALL readpar(ncid_met,'ejmax',completeSet,veg%ejmax,filename%met,           &
         nmetpatches,'def')
    CALL readpar(ncid_met,'vcmax',completeSet,veg%vcmax,filename%met,           &
         nmetpatches,'def')
    CALL readpar(ncid_met,'rp20',completeSet,veg%rp20,filename%met,             &
         nmetpatches,'def')
    CALL readpar(ncid_met,'rpcoef',completeSet,veg%rpcoef,filename%met,         &
         nmetpatches,'def')
    CALL readpar(ncid_met,'shelrb',completeSet,veg%shelrb,filename%met,         &
         nmetpatches,'def')
    CALL readpar(ncid_met,'xfang',completeSet,veg%xfang,filename%met,           &
         nmetpatches,'def')
    CALL readpar(ncid_met,'wai',completeSet,veg%wai,filename%met,               &
         nmetpatches,'def')
    CALL readpar(ncid_met,'vegcf',completeSet,veg%vegcf,filename%met,           &
         nmetpatches,'def')
    CALL readpar(ncid_met,'extkn',completeSet,veg%extkn,filename%met,           &
         nmetpatches,'def')
    CALL readpar(ncid_met,'tminvj',completeSet,veg%tminvj,filename%met,         &
         nmetpatches,'def')
    CALL readpar(ncid_met,'tmaxvj',completeSet,veg%tmaxvj,filename%met,         &
         nmetpatches,'def')
    CALL readpar(ncid_met,'vbeta',completeSet,veg%vbeta,filename%met,           &
         nmetpatches,'def')
    CALL readpar(ncid_met,'xalbnir',completeSet,veg%xalbnir,filename%met,       &
         nmetpatches,'def')
    CALL readpar(ncid_met,'meth',completeSet,veg%meth,filename%met,             &
         nmetpatches,'def')
    CALL readpar(ncid_met,'g0',completeSet,veg%g0,filename%met,            &
         nmetpatches,'def') ! Ticket #56
    CALL readpar(ncid_met,'g1',completeSet,veg%g1,filename%met,             &
         nmetpatches,'def') ! Ticket #56
    ok = NF90_INQ_VARID(ncid_met,'za',parID)
    IF(ok == NF90_NOERR) THEN ! if it does exist
       CALL readpar(ncid_met,'za',completeSet,rough%za_uv,filename%met,         &
            nmetpatches,'def')
       CALL readpar(ncid_met,'za',completeSet,rough%za_tq,filename%met,         &
            nmetpatches,'def')
    ELSE
       CALL readpar(ncid_met,'za_uv',completeSet,rough%za_uv,filename%met,      &
            nmetpatches,'def')
       CALL readpar(ncid_met,'za_tq',completeSet,rough%za_tq,filename%met,      &
            nmetpatches,'def')
    ENDIF
    CALL readpar(ncid_met,'zse',completeSet,soil%zse,filename%met,              &
         nmetpatches,'ms')
    CALL readpar(ncid_met,'ratecp',completeSet,bgc%ratecp,filename%met,         &
         nmetpatches,'ncp')
    CALL readpar(ncid_met,'ratecs',completeSet,bgc%ratecs,filename%met,         &
         nmetpatches,'ncs')

  END SUBROUTINE get_parameters_met

  !==============================================================================
  !
  ! Name: allocate_cable_vars
  !
  ! Purpose: Allocate CABLE's main variables.
  !
  ! CALLed from: load_parameters
  !              old_load_parameters
  !
  ! CALLs: alloc_cbm_var
  !
  !==============================================================================

  SUBROUTINE allocate_cable_vars(air,bgc,canopy,met,bal,                         &
       rad,rough,soil,ssnow,sum_flux,                  &
       veg,arraysize)
    TYPE (met_type), INTENT(INOUT)            :: met
    TYPE (air_type), INTENT(INOUT)            :: air
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE (bgc_pool_type), INTENT(INOUT)       :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (canopy_type), INTENT(INOUT)         :: canopy
    TYPE (roughness_type), INTENT(INOUT)      :: rough
    TYPE (radiation_type),INTENT(INOUT)       :: rad
    TYPE (sum_flux_type), INTENT(INOUT)       :: sum_flux
    TYPE (balances_type), INTENT(INOUT)       :: bal
    INTEGER, INTENT(IN)                       :: arraysize

    CALL alloc_cbm_var(air, arraysize)
    CALL alloc_cbm_var(bgc, arraysize)
    CALL alloc_cbm_var(canopy, arraysize)
    CALL alloc_cbm_var(met, arraysize)
    CALL alloc_cbm_var(bal, arraysize)
    CALL alloc_cbm_var(rad, arraysize)
    CALL alloc_cbm_var(rough, arraysize)
    CALL alloc_cbm_var(soil, arraysize)
    CALL alloc_cbm_var(ssnow, arraysize)
    CALL alloc_cbm_var(sum_flux, arraysize)
    CALL alloc_cbm_var(veg, arraysize)


    ! Allocate patch fraction variable:
    ALLOCATE(patch(arraysize))

  END SUBROUTINE allocate_cable_vars

END MODULE cable_input_module
!==============================================================================
