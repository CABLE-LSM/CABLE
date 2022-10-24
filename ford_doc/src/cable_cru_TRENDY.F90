MODULE CABLE_CRU

  USE netcdf                         ! Access to netcdf routines
       
  USE casa_ncdf_module, ONLY: HANDLE_ERR, GET_UNIT                       ! Finds an unused unit number for file opensYMDHMS2DOYSOD, DOYSOD2YMDHMS

  USE cable_IO_vars_module, ONLY: &  ! Selected cable_iovars.F90 variables:
       logn,            &             ! Log file unit number
       land_x, land_y,  &             ! Col (x) & row (y) indices of each land point in land mask (dimension mland)
       exists                         ! Only for exists%Snowf, which we will set to .FALSE. because there is no snow
  ! in CRU-NCEP. Setting this ensures snow will be determined in CABLE from temperature.

  IMPLICIT NONE

  ! Define a type for CRU-NCEP information, and the subtype METVALS

  TYPE CRU_MET_TYPE
     REAL, DIMENSION(:), ALLOCATABLE :: METVALS  ! Define a spatial vector of meteorology for one timestep
  END TYPE CRU_MET_TYPE

  TYPE CRU_TYPE
     INTEGER  :: mland                    ! Number of land cells
     INTEGER  :: NMET                     ! Number of met variable types (rain, lwdn etc) NOT INCLUDING prevTmax and nextTmin
     INTEGER  :: xdimsize, ydimsize       ! Landmask grid size dimensions (x=cols, y=rows)
     INTEGER  :: tdimsize                 ! Time dimension of metfiles (met data timesteps per annual file)
     INTEGER  :: CYEAR                    ! Current run year, same as CurYear, not necessarily the same as MetYear
     INTEGER  :: MetStart                 ! First year of met
     INTEGER  :: MetEnd                   ! Last year of met
     INTEGER  :: CTSTEP                   ! Current met data timestep (1 to tdimsize, i.e. 365 for CRU-NCEP annual daily files)
     INTEGER  :: DTsecs                   ! Model timestep in seconds, converted from namelist value in hours
     INTEGER  :: ktau                     ! Current model timestep, reset at the start of a new year of met
     INTEGER, DIMENSION(9) :: F_ID, V_ID  ! NetCDF object id's for files and variables (NetCDF bookkeeping stuff)
     REAL, DIMENSION(:) ,ALLOCATABLE :: AVG_LWDN ! Avg of one day's diurnal cycle of lwdn calculated by Swinbank. AVG_LWDN
     ! is used to rescale the diurnal cycle to match the day's CRUNCEP lwdn. (dim=mland)
     REAL, DIMENSION(:) ,ALLOCATABLE :: CO2VALS  ! Global annual CO2 values (dim is the number of years of data, or 1 if time-invariant)
     LOGICAL  :: DirectRead     ! Flag to do with reading small numbers of points efficiently. Set true for small numbers of points
     LOGICAL  :: LeapYears      ! Flag for whether leaps years occur, required by CABLE. Always false for CRUNCEP (no Feb 29th)
     LOGICAL,DIMENSION(:,:),ALLOCATABLE :: LandMask ! Logical landmask, true for land, false for non-land
     !
     CHARACTER(len=30)  :: Run            ! Where run type is      : "S0_TRENDY", "S1_TRENDY", "S2_TRENDY"
     CHARACTER(len=15)  :: CO2            ! CO2 takes value        : "static1860", "1860_1900", "1901_2015"
     CHARACTER(len=15)  :: Ndep            ! Ndep takes value        : "static1860", "1860_1900", "1901_2015"
     CHARACTER(len=15)  :: Forcing        ! Met Forcing takes value: "spinup",        "spinup", "1901_2015"
     !
     CHARACTER(len=200) :: BasePath       ! Full path for the location of data used for CRU runs "/x/y"
     CHARACTER(len=200) :: MetPath        ! Full path for the location of the met files "/x/y"
     CHARACTER(len=200) :: LandMaskFile   ! Land mask filename, without path
     CHARACTER(len=30) ,DIMENSION(9) :: VAR_NAME  ! Netcdf variable 'Name' for each type of met (dim=# of met vars). Note: Name, not 'Title'
     CHARACTER(len=200),DIMENSION(9) :: MetFile   ! Met file names incl metpath, constructed in CRU_GET_FILENAME (dim=# of met vars)
     TYPE(CRU_MET_TYPE), DIMENSION(11) :: MET  ! Met data vectors (METVALS) for one timestep, dim=# of met vars + 2 for prev Tmax and next Tmin
     REAL,  DIMENSION(:), ALLOCATABLE :: NdepVALS
     INTEGER :: NdepF_ID, NdepV_ID
     INTEGER  :: Ndep_CTSTEP   ! counter for Ndep in input file
  END TYPE CRU_TYPE

  TYPE (CRU_TYPE):: CRU  ! Define the variable CRU, of type CRU_TYPE

  ! Define local parameter names representing the position of each met var within variable MET.
  ! prevTmax and nextTmin are special cases of Tmax and Tmin that do not count as extra met variables per se.
  INTEGER, PRIVATE, PARAMETER :: &
       rain     =  1, &
       lwdn     =  2, &
       swdn     =  3, &
       pres     =  4, &
       qair     =  5, &
       tmax     =  6, &
       tmin     =  7, &
       uwind    =  8, &
       vwind    =  9, &
                                !
       prevTmax = 10, &
       nextTmin = 11

  ! Error status of various operations (mostly netcdf-related). Typically 0 means ok, > 0 means unexpected condition.
  INTEGER, PRIVATE :: ErrStatus

  REAL, PRIVATE, PARAMETER :: SecDay = 86400. ! Number of seconds in a day

  ! Filename prefix expected in the names of met files. Used by CRU_GET_FILENAME to construct met file names.
  CHARACTER(len=6), DIMENSION(9), PARAMETER, PRIVATE :: &
       PREF = (/ "rain  ", "lwdown", "swdown", "press ", "qair  ", "tmax  ", "tmin  ", "uwind ", "vwind " /)

CONTAINS

  !**************************************************************************************************

  SUBROUTINE CRU_INIT( CRU )

    ! Initialise the contents of the CRU defined type collection, from the CRU namelist file
    ! and by obtaining dimensions from the landmask

    !**************************************************************************************************

    USE cable_IO_vars_module, ONLY: &
         latitude, longitude, & ! (R) Lat and long of landcells only (?)
         nmetpatches,         & ! (I) Size of patch dimension in met file, if it exists
         mask,                & ! (I) Land/sea mask (1,0)
         metGrid,             & ! (C4) Either 'land' or 'mask' for whether the data are packed or not (?)
         sdoy, smoy, syear,   & ! (I) Start time day of year, month, year
         shod,                & ! (R) Start time hour of day
         xdimsize, ydimsize,  & ! (I) Size of grid dimensions
         lat_all, lon_all       ! (R) Grids with the lat or lon of each cell (i.e. repetition along rows/cols), for CABLE.

    USE cable_def_types_mod,  ONLY: mland  ! (I) Number of land cells

    IMPLICIT NONE

    TYPE (CRU_TYPE):: CRU

    INTEGER              :: ErrStatus  ! Error status returned by nc routines (zero=ok, non-zero=error)
    INTEGER              :: nmlunit    ! Unit number for reading namelist file
    INTEGER              :: FID        ! NetCDF id for the landmask file
    INTEGER              :: latID, lonID, timID  ! NetCDF ids for dimensions in the landmask file
    INTEGER              :: landID     ! NetCDF id for the landmask variable in the landmask file
    INTEGER              :: tdimsize   ! Time dimension in the met file, not used here apparently (??)
    INTEGER              :: landcnt    ! Manually incremented counter for the number of land cells
    INTEGER              :: xcol, yrow ! Column and row position in the data file grids
    INTEGER              :: imetvar    ! loop counter through met variables

    ! Temporary local names for CRU% variables as they are read from the namelist file.
    ! Note that CRU%CO2 and CRU%Forcing are assigned based on the value of Run, not read as options from the namelist file.
    LOGICAL              :: DirectRead = .FALSE.
    CHARACTER(len=30)    :: Run
    CHARACTER(len=200)   :: BasePath
    CHARACTER(len=200)   :: MetPath
    CHARACTER(len=200)   :: LandMaskFile
    REAL                 :: DThrs   ! CABLE timestep (hrs), converted immediately to integer seconds for CRU%DTsecs
    REAL,DIMENSION(:)     ,ALLOCATABLE :: CRU_lats, CRU_lons  ! Lat/long values for each grid rows/cols from landmask.
    INTEGER,DIMENSION(:,:),ALLOCATABLE :: landmask

    ! Flag for errors
    LOGICAL              :: ERR = .FALSE.

    NAMELIST /CRUNML/ BasePath, MetPath, LandMaskFile, Run, DThrs, DirectRead

    ! Read CRU namelist settings
    CALL GET_UNIT(nmlunit)  ! CABLE routine finds spare unit number
    OPEN (nmlunit,FILE="cru.nml",STATUS='OLD',ACTION='READ')
    READ (nmlunit,NML=CRUNML)
    CLOSE(nmlunit)

    ! Assign namelist settings to corresponding CRU defined-type elements
    CRU%BasePath     = BasePath
    CRU%MetPath      = MetPath
    CRU%LandMaskFile = LandMaskFile
    CRU%Run          = Run
    CRU%DTsecs       = INT(DThrs * 3600.)  ! in seconds
    CRU%DirectRead   = DirectRead

    ! Assign Forcing and CO2 labels based only on the value of CRU%Run
    SELECT CASE (TRIM(CRU%Run))
    CASE( "S0_TRENDY" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep     = "static1860"
       WRITE(*   ,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
       WRITE(logn,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
    CASE( "S0_TRENDY_CO2" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "1901_2015"
       CRU%Ndep     = "static1860"
       WRITE(*   ,*)"Run = 'S0_CO2': Therefore Forcing = 'S0_CO2', CO2 = '1860_2015'"
       WRITE(logn,*)"Run = 'S0_CO2': Therefore Forcing = 'S0_CO2', CO2 = '1860_2015'"
    CASE( "S0_TRENDY_Ndep" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep     = "1901_2015"
       WRITE(*   ,*)"Run = 'S0_Ndep': Therefore Forcing = 'S0_Ndep', Ndep = '1860_2015'"
       WRITE(logn,*)"Run = 'S0_Ndep': Therefore Forcing = 'S0_Ndep', Ndep = '1860_2015'"
    CASE( "S0_TRENDY_Precip" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep     = "static1860"
       WRITE(*   ,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
       WRITE(logn,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
    CASE( "S0_TRENDY_Temp" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep     = "static1860"
       WRITE(*   ,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
       WRITE(logn,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
    CASE( "S0_TRENDY_Temp_Precip" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep     = "static1860"
       WRITE(*   ,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
       WRITE(logn,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
    CASE( "S0_TRENDY_CO2_Temp" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "1901_2015"
       CRU%Ndep     = "static1860"
       WRITE(*   ,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
       WRITE(logn,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
    CASE( "S0_TRENDY_CO2_Precip" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "1901_2015"
       CRU%Ndep     = "static1860"
       WRITE(*   ,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
       WRITE(logn,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
    CASE( "S0_TRENDY_CO2_Temp_Precip" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "1901_2015"
       CRU%Ndep     = "static1860"
       WRITE(*   ,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
       WRITE(logn,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
    CASE( "S1_TRENDY" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "1860_1900"
       CRU%Ndep     = "1860_1900"
       WRITE(*   ,*)"Run = 'S1_TRENDY': Therefore Forcing = 'spinup', CO2 = '1860_1900'"
       WRITE(logn,*)"Run = 'S1_TRENDY': Therefore Forcing = 'spinup', CO2 = '1860_1900'"
    CASE( "S2_TRENDY" )
       CRU%Forcing = "1901_2015"
       CRU%CO2     = "1901_2015"
       CRU%Ndep     = "1901_2015"
       WRITE(*   ,*)"Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
       WRITE(logn,*)"Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
    CASE( "S2_TRENDY_precip" )
       CRU%Forcing = "1901_2015"
       CRU%CO2     = "1901_2015"
       CRU%Ndep     = "1901_2015"
       WRITE(*   ,*)"Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
       WRITE(logn,*)"Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
    CASE( "S2_TRENDY_precip0" )
       CRU%Forcing = "1901_2015"
       CRU%CO2     = "1901_2015"
       CRU%Ndep     = "1901_2015"
       WRITE(*   ,*)"Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
       WRITE(logn,*)"Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
    CASE  default
       WRITE(*   ,*)"Wrong CRU%Run: ",CRU%Run
       WRITE(*   ,*)"Use: S0_TRENDY, S1_TRENDY, or S2_TRENDY!"
       WRITE(logn,*)"Wrong CRU%Run: ",CRU%Run
       WRITE(logn,*)"Use: S0_TRENDY, S1_TRENDY, or S2_TRENDY!"
       ERR = .TRUE.
    END SELECT

    ! Print settings
    WRITE(*   ,*)"========================================= CRU ============"
    WRITE(*   ,*)"CRU settings chosen:"
    WRITE(*   ,*)" BasePath: ",TRIM(CRU%BasePath)
    WRITE(*   ,*)" LandMask: ",TRIM(CRU%LandMaskFile)
    WRITE(*   ,*)" Run                : ",TRIM(CRU%Run)
    WRITE(*   ,*)" Forcing (assigned) : ",TRIM(CRU%Forcing)
    WRITE(*   ,*)" CO2     (assigned) : ",TRIM(CRU%CO2)
    WRITE(*   ,*)" Ndep     (assigned) : ",TRIM(CRU%Ndep)
    WRITE(*   ,*)" DT(secs): ",CRU%DTsecs
    WRITE(logn,*)"========================================= CRU ============"
    WRITE(logn,*)"CRU settings chosen:"
    WRITE(logn,*)" BasePath: ",TRIM(CRU%BasePath)
    WRITE(logn,*)" LandMask: ",TRIM(CRU%LandMaskFile)
    WRITE(logn,*)" Run                : ",TRIM(CRU%Run)
    WRITE(logn,*)" Forcing (assigned) : ",TRIM(CRU%Forcing)
    WRITE(logn,*)" CO2     (assigned) : ",TRIM(CRU%CO2)
    WRITE(logn,*)" Ndep     (assigned) : ",TRIM(CRU%Ndep)
    WRITE(logn,*)" DT(secs): ",CRU%DTsecs

    ! Error trap for bad namelist.
    IF ( ERR ) THEN
       WRITE(logn,*)"Invalid settings in CRU_INIT"
       STOP "Invalid settings in CRU_INIT"
    ENDIF

    ! If this is a S0_TRENDY run look for met data in the spinup directory instead.
    !IF (TRIM(CRU%Run) .EQ. "S0_TRENDY") THEN
    !   CRU%MetPath = TRIM(CRU%MetPath)//"/spinup_data"
    !ENDIF

    ! Set variable names to their NetCDF 'Names' (i.e. not their 'Titles')
    CRU%NMET = 9
    CRU%VAR_NAME(rain)  = "Total_Precipitation"
    CRU%VAR_NAME(lwdn)  = "Incoming_Long_Wave_Radiation"
    CRU%VAR_NAME(swdn)  = "Incoming_Short_Wave_Radiation"
    CRU%VAR_NAME(pres)  = "Pression"
    CRU%VAR_NAME(qair)  = "Air_Specific_Humidity"
    ! CRU%VAR_NAME(tmax)  = "maximum_6h_air_temperature"
    ! CRU%VAR_NAME(tmin)  = "minimum_6h_air_temperature"
    CRU%VAR_NAME(tmax)  = "maximum_air_temperature"
    CRU%VAR_NAME(tmin)  = "minimum_air_temperature"
    CRU%VAR_NAME(uwind) = "U_wind_component"
    CRU%VAR_NAME(vwind) = "V_wind_component"

    WRITE(*   ,*)"========================================= CRU ============"
    WRITE(logn,*)"========================================= CRU ============"

    ! Now read landmask file
    ! Landmask file into init! Get LAt, LON etc. from there
    ! LMFILE = TRIM(CRU%LandMaskFile)
    WRITE(*   ,*) 'Opening CRU landmask file: ',TRIM(LandMaskFile)
    WRITE(logn,*) 'Opening CRU landmask file: ',TRIM(LandMaskFile)

    ! Open the land mask file
    ErrStatus = NF90_OPEN(TRIM(LandMaskFile), NF90_NOWRITE, FID)
    CALL HANDLE_ERR(ErrStatus, "Opening CRU Land-mask file"//TRIM(LandMaskFile))

    ! Latitude: Get the dimension ID, find the size of the dimension, assign it to CRU.
    ErrStatus = NF90_INQ_DIMID(FID,'latitude',latID)
    ErrStatus = NF90_INQUIRE_DIMENSION(FID,latID,len=ydimsize)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'lat'"//TRIM(LandMaskFile))
    CRU%ydimsize = ydimsize

    ! Collect the latitudes into CRU_lats
    ALLOCATE( CRU_lats ( ydimsize ) )
    ErrStatus = NF90_INQ_VARID(FID,'latitude',latID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'latitudes'"//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,latID,CRU_lats)
    CALL HANDLE_ERR(ErrStatus, "Reading 'latitudes'"//TRIM(LandMaskFile))

    ! Longitude: Get the dimension ID, find the size of the dimension, assign it to CRU.
    ErrStatus = NF90_INQ_DIMID(FID,'longitude',lonID)
    ErrStatus = NF90_INQUIRE_DIMENSION(FID,lonID,len=xdimsize)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'lon'"//TRIM(LandMaskFile))
    CRU%xdimsize = xdimsize

    ! Collect the longitudes into CRU_lons
    ALLOCATE( CRU_lons ( xdimsize ) )
    ErrStatus = NF90_INQ_VARID(FID,'longitude',lonID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'longitudes'"//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,lonID,CRU_lons)
    CALL HANDLE_ERR(ErrStatus, "Reading 'longitudes'"//TRIM(LandMaskFile))

    ! Allocate the landmask arrays for...
    ALLOCATE( CRU%landmask ( xdimsize, ydimsize) )  ! Passing out to other CRU routines (logical)
    ALLOCATE( landmask ( xdimsize, ydimsize) )      ! Local use in this routine (integer)
    ALLOCATE ( mask( xdimsize, ydimsize) )          ! Use by CABLE

    ! Check that the land mask variable is called "land" in the land mask file,
    ! and read it into local variable landmask
    ErrStatus = NF90_INQ_VARID(FID,'land',landID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'land' "//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,landID,landmask)
    CALL HANDLE_ERR(ErrStatus, "Reading 'land' "//TRIM(LandMaskFile))

    ! Convert the integer landmask into the logical CRU%landmask
    WHERE ( landmask .GT. 0 )
       CRU%landmask = .TRUE.
       mask           = 1
    ELSEWHERE
       CRU%landmask = .FALSE.
       mask           = 0
    END WHERE

    ! Count the number of land cells -> mland
    CRU%mland = COUNT(CRU%landmask)

    ! Allocate CABLE land-only vectors for lat/long and row/col values/indices.
    ALLOCATE( latitude(CRU%mland), longitude(CRU%mland) )
    ALLOCATE( land_y  (CRU%mland), land_x   (CRU%mland) )

    ! Allocate vectors for each of the different met quantities, including extra
    ! prev/next temperatures for the Cesarracio temperature calculations in the
    ! weather generator.
    DO imetvar = 1, CRU%NMET
       ALLOCATE( CRU%MET(imetvar)%METVALS(CRU%mland) )
    END DO
    ALLOCATE( CRU%MET(prevTmax)%METVALS(CRU%mland) )
    ALLOCATE( CRU%MET(nextTmin)%METVALS(CRU%mland) )
    ! allocate array for Nitrogen deposition input data
    ALLOCATE( CRU%NdepVALS(CRU%mland) )

    ! Copy the col/row and lat/long positions of each land cell into the corresponding
    ! land only CABLE vectors. Q: We know mland at this point. Why not use landcnt to confirm
    ! the correct value of mland?
    landcnt = 1
    DO yrow = 1, ydimsize
       DO xcol = 1, xdimsize
          IF ( .NOT. CRU%landmask(xcol,yrow) ) CYCLE   ! Go to next iteration if not a land cell

          !          WRITE(6,FMT='(A15,I5,2(1X,F8.2),2(1x,I3))')"i, lo,la, xcol,yrow",landcnt,CRU_lons(xcol),CRU_lats(yrow),xcol, yrow

          ! C
          land_x(landcnt)    = xcol
          land_y(landcnt)    = yrow
          longitude(landcnt) = CRU_lons(xcol)
          latitude(landcnt)  = CRU_lats(yrow)
          landcnt = landcnt + 1
       END DO
    END DO

    ! Set global CABLE variables
    metGrid     = "mask"
    ALLOCATE( mask(xdimsize, ydimsize) )
    mask        = landmask
    mland       = CRU%mland
    nmetpatches = 1
    ALLOCATE( lat_all(xdimsize, ydimsize), lon_all(xdimsize, ydimsize) )
    DO xcol = 1, xdimsize
       lat_all(xcol,:) = CRU_lats
    END DO
    DO yrow = 1, ydimsize
       lon_all(:,yrow) = CRU_lons
    END DO

    ! CABLE TIME-UNITS needed by load-parameters (only on CABLE_init)
    shod        = 0.
    sdoy        = 1
    smoy        = 1
    syear       = CRU%CYEAR

    ! Used to rescale the diurnal cycle from Swinbank calculation to match CRU-NCEP provided value.
    ALLOCATE( CRU%AVG_LWDN(mland) )

    DEALLOCATE ( landmask, CRU_lats, CRU_lons )

    ErrStatus = NF90_CLOSE(FID)
    CALL HANDLE_ERR(ErrStatus, "Closing mask-file"//TRIM(LandMaskFile))

  END SUBROUTINE CRU_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CRU_GET_FILENAME ( CRU, cyear, par, FN )

    ! Build the filename FN: One annual file of daily met for one met quantity.

    IMPLICIT NONE

    TYPE(CRU_TYPE), INTENT(INOUT) :: CRU  ! Information about CRU
    INTEGER,              INTENT(IN)    :: cyear  ! Current year as an integer
    INTEGER,              INTENT(IN)    :: par    ! Index (1-9) of which met quantity will be sought
    CHARACTER(LEN=200),   INTENT(OUT)   :: FN     ! Met filename (outgoing)
    INTEGER   :: i, idx
    CHARACTER(4)   :: cy  ! Character representation of cyear
    CHARACTER(200) :: mp  ! Local repr of met path

    ! Create a character version of the year for building that part of the filename.
    WRITE(cy,FMT='(I4)')cyear

    ! Initialise the filename with the met path
    mp = CRU%MetPath
    FN = TRIM(mp)

    ! Build the rest of the filename according to the value of par, which references 11 possible
    ! types of met through the parameter names rain, lwdn, etc.
!!$    SELECT CASE ( par )
!!$    CASE(rain) ; FN = TRIM(FN)//"/rain/cruncep2015_1_rain_"//cy//".daymean.nc"
!!$    CASE(lwdn) ; FN = TRIM(FN)//"/lwdown/cruncep2015_1_lwdown_"//cy//".daymean.nc"
!!$    CASE(swdn) ; FN = TRIM(FN)//"/swdown/cruncep2015_1_swdown_"//cy//".daymean.nc"
!!$    CASE(pres) ; FN = TRIM(FN)//"/press/cruncep2015_1_press_"//cy//".daymean.nc"
!!$    CASE(qair) ; FN = TRIM(FN)//"/qair/cruncep2015_1_qair_"//cy//".daymean.nc"
!!$    CASE(tmax,PrevTmax) ; FN = TRIM(FN)//"/tmax/cruncep2015_1_tair_"//cy//".daymax.nc"
!!$    CASE(tmin,NextTmin) ; FN = TRIM(FN)//"/tmin/cruncep2015_1_tair_"//cy//".daymin.nc"
!!$    CASE(uwind) ; FN = TRIM(FN)//"/uwind/cruncep2015_1_uwind_"//cy//".daymean.nc"
!!$    CASE(vwind) ; FN = TRIM(FN)//"/vwind/cruncep2015_1_vwind_"//cy//".daymean.nc"
!!$    END SELECT


    SELECT CASE ( par )
    CASE(rain) ; FN = TRIM(FN)//"/rain/cruncepV8_rain_"//cy//".daytot.nc"
    CASE(lwdn) ; FN = TRIM(FN)//"/lwdown/cruncepV8_lwdown_"//cy//".daymean.nc"
    CASE(swdn) ; FN = TRIM(FN)//"/swdown/cruncepV8_swdown_"//cy//".daymean.nc"
    CASE(pres) ; FN = TRIM(FN)//"/press/cruncepV8_press_"//cy//".daymean.nc"
    CASE(qair) ; FN = TRIM(FN)//"/qair/cruncepV8_qair_"//cy//".daymean.nc"
    CASE(tmax,PrevTmax) ; FN = TRIM(FN)//"/tmax/cruncepV8_tmax_"//cy//".daymax.nc"
    CASE(tmin,NextTmin) ; FN = TRIM(FN)//"/tmin/cruncepV8_tmin_"//cy//".daymin.nc"
    CASE(uwind) ; FN = TRIM(FN)//"/uwind/cruncepV8_uwind_"//cy//".daymean.nc"
    CASE(vwind) ; FN = TRIM(FN)//"/vwind/cruncepV8_vwind_"//cy//".daymean.nc"
    END SELECT

  END SUBROUTINE CRU_GET_FILENAME

  !**************************************************************************************************

  SUBROUTINE GET_CRU_CO2( CRU, CO2air )

    ! Get CO2 values for use with a CRU-NCEP run. Assign a static 1860 value if specified otherwise
    ! on the first call read all the annual values from a file into the CRU%CO2VALS array. On the first
    ! and subsequent

    IMPLICIT NONE

    TYPE(CRU_TYPE) :: CRU           ! All the info needed for CRU met runs
    REAL, INTENT(OUT)    :: CO2air  ! A single annual value of CO2air in ppm for the current year.

    INTEGER              :: i, iunit, iyear, IOS = 0
    CHARACTER            :: CO2FILE*200
    LOGICAL,        SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine

    ! For S0_TRENDY, use only static 1860 CO2 value and return immediately
    IF ( TRIM(CRU%CO2) .EQ. "static1860") THEN
       CO2air = 286.42   ! CO2 in ppm for 1860
       RETURN

       ! If not S0_TRENDY, varying CO2 values will be used...
    ELSE

       ! On the first call, allocate the CRU%CO2VALS array to store the entire history of annual CO2
       ! values, open the (ascii) CO2 file and read the values into the array.
       IF (CALL1) THEN
          ALLOCATE( CRU%CO2VALS( 1750:2016 ) )
          CO2FILE = TRIM(CRU%BasePath)//"/co2/1750_2015_globalCO2_time_series.csv"
          CALL GET_UNIT(iunit)
          OPEN (iunit, FILE=TRIM(CO2FILE), STATUS="OLD", ACTION="READ")
          DO WHILE( IOS .EQ. 0 )
             READ(iunit, FMT=*, IOSTAT=IOS) iyear, CRU%CO2VALS(iyear)
          END DO
          CLOSE(iunit)

          CALL1 = .FALSE.

       END IF

       ! In all varying CO2 cases, return the element of the array for the current year
       ! as a single CO2 value.
       !
       CO2air = CRU%CO2VALS( CRU%CYEAR )

    END IF

  END SUBROUTINE GET_CRU_CO2

  !**************************************************************************************************

  SUBROUTINE GET_CRU_Ndep( CRU )

    ! Get Ndep values for use with a CRU-NCEP run. Assign a static 1860 value if specified otherwise
    ! on the first call read all the annual values from a file into the CRU%CO2VALS array. On the first
    ! and subsequent

    IMPLICIT NONE

    TYPE(CRU_TYPE), INTENT(INOUT) :: CRU           ! All the info needed for CRU met runs
    REAL    :: tmparr(720,360)        ! Temporary array for reading one day of met before
    ! packing into CRU%NdepVALS(k)

    INTEGER              :: i, iunit, iyear, IOS = 0, k, t
    INTEGER :: xds, yds        ! Ndep file dimensions of long (x), lat (y)

    LOGICAL,        SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine
    CHARACTER(200) :: NdepFILE

    ! Abbreviate dimensions for readability.
    xds = CRU%xdimsize
    yds = CRU%ydimsize

    ! For S0_TRENDY, use only static 1860 CO2 value and return immediately



    ! On the first call, allocate the CRU%CO2VALS array to store the entire history of annual CO2
    ! values, open the (ascii) CO2 file and read the values into the array.
    IF (CALL1) THEN

       NdepFILE = TRIM(CRU%BasePath)// &
            "/ndep/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual.nc"

       ! Open the NDep and access the variables by their name and variable id.
       WRITE(*   ,*) 'Opening ndep data file: ', NdepFILE
       WRITE(logn,*) 'Opening ndep data file: ', NdepFILE


       ErrStatus = NF90_OPEN(TRIM(NdepFILE), NF90_NOWRITE, CRU%NdepF_ID)
       CALL HANDLE_ERR(ErrStatus, "Opening CRU file "//NdepFILE )
       ErrStatus = NF90_INQ_VARID(CRU%NdepF_ID,'N_deposition', CRU%NdepV_ID)
       CALL HANDLE_ERR(ErrStatus, "Inquiring CRU var "//"N_deposition"// &
            " in "//NdepFILE )

       ! Set internal counter
       CRU%Ndep_CTSTEP = 1

       IF ( TRIM(CRU%Ndep) .EQ. "static1860") THEN
          ! read Ndep at year 1860 (noting that file starts at 1850)
          CRU%Ndep_CTSTEP = 11
          t =  CRU%Ndep_CTSTEP
          ErrStatus = NF90_GET_VAR(CRU%NdepF_ID, CRU%NdepV_ID, tmparr, &
               start=(/1,1,t/),count=(/xds,yds,1/) )
          CALL HANDLE_ERR(ErrStatus, "Reading from "//NdepFILE )
          DO k = 1, CRU%mland
             CRU%NdepVALS(k) = tmparr( land_x(k), land_y(k) )
          END DO


       END IF
       CALL1 = .FALSE.
    END IF

    IF ( TRIM(CRU%Ndep) .NE. "static1860") THEN

       ! read Ndep at current year (noting that file starts at 1850 and ends in 2015)
       CRU%Ndep_CTSTEP = MIN(CRU%CYEAR, 2015) - 1850 + 1
       t =  CRU%Ndep_CTSTEP
       ErrStatus = NF90_GET_VAR(CRU%NdepF_ID, CRU%NdepV_ID, tmparr, &
            start=(/1,1,t/),count=(/xds,yds,1/) )
       CALL HANDLE_ERR(ErrStatus, "Reading from "//NdepFILE )
       DO k = 1, CRU%mland
          CRU%NdepVALS(k) = tmparr( land_x(k), land_y(k) )
       END DO

    END IF





  END SUBROUTINE GET_CRU_Ndep

  !**************************************************************************************************

  SUBROUTINE OPEN_CRU_MET( CRU )

    ! Opens each of the met files required for one year. This is where the distinction is made between
    ! the nominal run year (CYEAR) and the year of met required (MetYear), which is different for
    ! S0_TRENDY and S1_TRENDY than for a standard run (S2_TRENDY).

    USE cable_IO_vars_module, ONLY: timeunits ! (Char33) Name of time units read from nc file

    IMPLICIT NONE

    TYPE( CRU_TYPE ), INTENT(INOUT) :: CRU ! All CRU-NCEP related quantities and flags

    INTEGER             :: iVar            ! Loop counter through met variables
    INTEGER             :: tID             ! Numerical variable identifier returned by NetCDF routines,
    ! in this case for time. Needed to retrieve the time units.
    INTEGER             :: MetYear         ! Year of met to access. Equals CYEAR for normal runs, but
    ! must be calculated for S0_TRENDY and initialisation runs.
    INTEGER, SAVE       :: RunStartYear    ! The value of CRU%CYEAR on the first call, also equals syear.
    ! Allows the calculation of MetYear during S0_TRENDY and init runs.
    LOGICAL, SAVE       :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine

    ! Keep the initial value of CYEAR for calculation of different MetYear if required.
    !IF (CALL1) RunStartYear = 1710 ! edit vh !
    IF (CALL1) RunStartYear = 1691 ! edit vh !
    DO iVar = 1, CRU%NMET  ! For each met variable

       ! For S0_TRENDY and initialisation, calculate the required met year for repeatedly cycling through the
       ! 30 years of 1901-1930 spinup meteorology. For normal runs 1901-2015, MetYear = CYEAR.
!!$    IF ( TRIM(CRU%Run) .EQ. 'S0_TRENDY' .OR.  ( TRIM(CRU%Run) .EQ. 'S1_TRENDY' )) THEN
!!$      MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
!!$    ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY' ) THEN
!!$      MetYear = CRU%CYEAR
!!$    ENDIF
       IF ( TRIM(CRU%Run) .EQ. 'S0_TRENDY' .OR.  ( TRIM(CRU%Run) .EQ. 'S1_TRENDY' ) &
            .OR.  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2') &
            .OR.  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Ndep' )) THEN
          MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
       ELSEIF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Precip' .OR. &
            TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Precip'.OR. &
            TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp_Precip'.OR. &
            TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp_Precip' ) THEN
          IF (iVar.EQ.1) THEN
             MetYear = CRU%CYEAR
          ELSE
             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
          ENDIF

          IF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp_Precip'.OR. &
               TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp_Precip' ) THEN
             IF (iVar.EQ.6 .OR. iVar.EQ.7) THEN
                MetYear = CRU%CYEAR
             ELSE
                MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
             ENDIF
          ENDIF

       ELSEIF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp' .OR. &
            TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp' ) THEN

          IF (iVar.EQ.6 .OR. iVar.EQ.7) THEN
             MetYear = CRU%CYEAR
          ELSE
             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
          ENDIF
       ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY' ) THEN
          MetYear = CRU%CYEAR
       ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY_precip0' ) THEN
          IF (iVar.EQ.1) THEN
             ! special for baseline precip
             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
          ELSE
             MetYear = CRU%CYEAR
          ENDIF
       ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY_precip' ) THEN
          IF (iVar.NE.1) THEN
             ! special for baseline non-precip
             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
          ELSE
             MetYear = CRU%CYEAR
          ENDIF
       ENDIF

       CALL CRU_GET_FILENAME( CRU, MetYear, iVar, CRU%MetFile(iVar) ) ! Call routine to build the filenames.

       ! Open the new met files and access the variables by their name and variable id.
       WRITE(*   ,*) 'Opening met data file: ', CRU%MetFile(iVar)
       WRITE(logn,*) 'Opening met data file: ', CRU%MetFile(iVar)

       ErrStatus = NF90_OPEN(TRIM(CRU%MetFile(iVar)), NF90_NOWRITE, CRU%F_ID(iVar))
       CALL HANDLE_ERR(ErrStatus, "Opening CRU file "//CRU%MetFile(iVar) )
       ErrStatus = NF90_INQ_VARID(CRU%F_ID(iVar),TRIM(CRU%VAR_NAME(iVar)), CRU%V_ID(iVar))
       CALL HANDLE_ERR(ErrStatus, "Inquiring CRU var "//TRIM(CRU%VAR_NAME(iVar))// &
            " in "//CRU%MetFile(iVar) )
    END DO

    ! Set internal counter
    CRU%CTSTEP = 1

    CALL1 = .FALSE. ! No longer the first call (saved).

  END SUBROUTINE OPEN_CRU_MET

  !**************************************************************************************************

  SUBROUTINE CRU_GET_DAILY_MET( CRU, LastDayOfYear, LastYearOfMet )

    IMPLICIT NONE

    TYPE(CRU_TYPE) :: CRU
    LOGICAL, INTENT(IN)  :: LastDayOfYear, LastYearOfMet
    REAL    :: tmparr(720,360)        ! Temporary array for reading one day of met before
    ! packing into CRU%MET(iVar)%METVALS(k)
    REAL    :: tmp, stmp(365)
    INTEGER :: iVar, ii, k, x, y, realk
    INTEGER :: t, tplus1              ! The current and next timestep
    INTEGER :: fid, vid, tid          ! Netcdf id's for file, variable, and time
    INTEGER :: xds, yds, tds          ! Metfile dimensions of long (x), lat (y), and time (t)
    INTEGER :: MetYear                ! Year of meteorology currently in use
    INTEGER :: NextMetYear            ! Next met year: Where to look for the nextTmin on Dec 31st
    CHARACTER(LEN=200) :: filename

    INTEGER, SAVE       :: RunStartYear    ! The value of CRU%CYEAR on the first call, also equals syear.
    ! Allows the calculation of MetYear during S0_TRENDY and init runs.
    LOGICAL, SAVE :: CALL1 = .TRUE.   ! A *local* variable recording the first call of this routine


    ! If first call...
    ! Keep the initial value of CYEAR for calculation of different MetYear if required.
    IF (CALL1) THEN
       !RunStartYear = CRU%CYEAR
       RunStartYear = 1691
       ! If this is not the first call, capture the existing Tmax from the previous day as the
       ! 'previous Tmax' before reading another one. Move the existing next day's Tmin into the current
       ! Tmin before reading a new 'next Tmin'. These previous and next values are required for the
       ! Cesaraccio et al. algorithm used by the weather generator to interpolate daily into subdiurnal
       ! temperatures.
    ELSE
       CRU%MET(prevTmax)%METVALS(:) = CRU%MET(  Tmax  )%METVALS(:)
       CRU%MET(  Tmin  )%METVALS(:) = CRU%MET(NextTmin)%METVALS(:)
    ENDIF

    ! For S0_TRENDY and initialisation, calculate the year of meteorology as mod 50, so we repeatedly cycle
    ! through the 50 years of 1951-2000 spinup meteorology. For normal runs 1901-2015, MetYear = CYEAR.
    ! Stop with error for anything else.


!!$  IF ( TRIM(CRU%Run) .EQ. 'S0_TRENDY' .OR.  ( TRIM(CRU%Run) .EQ. 'S1_TRENDY' )) THEN
!!$    MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
!!$  ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY' ) THEN
!!$    MetYear = CRU%CYEAR
!!$  ELSE
!!$    STOP 'Error in cable_cru.F90: CRU%Run not S0_TRENDY, S1_TRENDY, or 1901-2015'
!!$  ENDIF

    !print *, "runstartyear, metyear", runstartyear, metyear

    ! Abbreviate dimensions for readability.
    xds = CRU%xdimsize
    yds = CRU%ydimsize

    ! Loop through all 9 met variables (not including prevTmax and nextTmin, which are addressed
    ! separately as special cases of Tmax and Tmin)

    !print *,  "CRU%CTSTEP, LastDayOfYear, LastYearOfMet", CRU%CTSTEP, LastDayOfYear, LastYearOfMet
    DO iVar = 1, CRU%NMET

       IF ( TRIM(CRU%Run) .EQ. 'S0_TRENDY' .OR.  ( TRIM(CRU%Run) .EQ. 'S1_TRENDY' ) &
            .OR.  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2') &
            .OR.  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Ndep' )) THEN
          MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
       ELSEIF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Precip' &
            .OR.  TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Precip'.OR. &
            TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp_Precip'.OR. &
            TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp_Precip' ) THEN
          IF (iVar.EQ.1) THEN
             MetYear = CRU%CYEAR
          ELSE
             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
          ENDIF
          IF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp_Precip'.OR. &
               TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp_Precip' ) THEN
             IF (iVar.EQ.6 .OR. iVar.EQ.7) THEN
                MetYear = CRU%CYEAR
             ELSE
                MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
             ENDIF
          ENDIF

       ELSEIF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp' &
            .OR.  TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp') THEN
          IF (iVar.EQ.6 .OR. iVar.EQ.7) THEN
             MetYear = CRU%CYEAR
          ELSE
             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
          ENDIF
       ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY' ) THEN
          MetYear = CRU%CYEAR
       ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY_precip0' ) THEN
          IF (iVar.EQ.1) THEN
             ! special for baseline precip
             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
          ELSE
             MetYear = CRU%CYEAR
          ENDIF
       ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY_precip' ) THEN
          IF (iVar.NE.1) THEN
             ! special for baseline non-precip
             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
          ELSE
             MetYear = CRU%CYEAR
          ENDIF
       ENDIF


       SELECT CASE (iVar)

          !--------------------------------------------------------------------------------------------------

       CASE (Tmin) ! When Tmin comes up we will jump ahead to deal with nextTmin, since we
          ! are assigning Tmin from the previous value of nextTmin

          ! Bump the variable index ii to be the value for nextTmin and the time index ahead by 1 day.
          ii = nextTmin
          t  = CRU%CTSTEP
          tplus1 = CRU%CTSTEP + 1

          ! The first call exception for Tmin is to read both the first (t) and second (t+1) timesteps.
          ! The first goes into Tmin, the second into nextTmin.
          IF ( CALL1 ) THEN

             ! In the case of a small number of points, it is more efficient to pull the data for each land point
             ! directly from the file (DirectRead) rather than reading the whole array and extracting the land
             ! points from that.
             IF ( CRU%DirectRead ) THEN

                ! For each land cell, read the first day value of Tmin into tmp, and copy it to the METVALS vector.

                DO k = 1, CRU%mland

                   ! t (first value) into iVar (tmin)
                   ErrStatus = NF90_GET_VAR( CRU%F_ID(iVar), CRU%V_ID(iVar), tmp, &
                        start=(/land_x(k),land_y(k),t/) )
                   CALL HANDLE_ERR(ErrStatus, "Reading direct from "//CRU%MetFile(iVar) )
                   CRU%MET(iVar)%METVALS(k) = tmp

                   ! tplus1 (second value) into ii (nextTmin)
                   ErrStatus = NF90_GET_VAR( CRU%F_ID(iVar), CRU%V_ID(iVar), tmp, &
                        start=(/land_x(k),land_y(k),tplus1/) )
                   CALL HANDLE_ERR(ErrStatus, "Reading direct from "//CRU%MetFile(iVar) )
                   CRU%MET(ii)%METVALS(k) = tmp

                END DO

                ! Not DirectRead: Read the whole spatial grid of first day Tmin into tmparr, and copy it to the
                ! METVALS vector.
             ELSE

                ! t (first value) into iVar (tmin)
                ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
                     start=(/1,1,t/),count=(/xds,yds,1/) )
                CALL HANDLE_ERR(ErrStatus, "Reading from "//CRU%MetFile(iVar) )

                DO k = 1, CRU%mland
                   CRU%MET(iVar)%METVALS(k) = tmparr( land_x(k), land_y(k) )
                END DO

                ! tplus1 (second value) into ii (nextTmin)
                ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
                     start=(/1,1,tplus1/),count=(/xds,yds,1/) )
                CALL HANDLE_ERR(ErrStatus, "Reading from "//CRU%MetFile(iVar) )

                DO k = 1, CRU%mland
                   CRU%MET(ii)%METVALS(k) = tmparr( land_x(k), land_y(k) )
                END DO

             ENDIF  ! End of if DirectRead

          END IF   ! End of if CALL1

          ! If this is the last day of the year we will need to get the nextTmin
          ! value from the next year's met file.
          IF ( LastDayOfYear ) THEN

             IF ( LastYearOfMet ) THEN  ! If there is no more met to open...

                CONTINUE ! Do nothing. The current value of nextTmin will just be reused.

             ELSE  ! There is another year of Tmin data available

                t = 1   ! Time index is set to the first day of the next year

                ! Add one to the calculation of MetYear
!!$          IF ( TRIM(CRU%Run) .EQ. 'S0_TRENDY' .OR.  ( TRIM(CRU%Run) .EQ. 'S1_TRENDY' )) THEN
!!$            NextMetYear = 1901 + MOD(CRU%CYEAR + 1 - RunStartYear,30)
!!$          ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY' ) THEN
!!$            NextMetYear = CRU%CYEAR + 1
!!$          ENDIF
                IF ( TRIM(CRU%Run) .EQ. 'S0_TRENDY' .OR.  ( TRIM(CRU%Run) .EQ. 'S1_TRENDY' ) &
                     .OR.  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2') &
                     .OR.  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Ndep' )) THEN
                   NextMetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
                ELSEIF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Precip'&
                     .OR.  TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Precip' .OR. &
                     TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp_Precip'.OR. &
                     TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp_Precip' ) THEN
                   IF (iVar.EQ.1) THEN
                      NextMetYear = CRU%CYEAR
                   ELSE
                      NextMetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
                   ENDIF
                   IF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp_Precip'.OR. &
                        TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp_Precip' ) THEN
                      IF (iVar.EQ.6 .OR. iVar.EQ.7) THEN
                         NextMetYear = CRU%CYEAR
                      ELSE
                         NextMetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
                      ENDIF
                   ENDIF


                ELSEIF  ( TRIM(CRU%Run) .EQ. 'S0_TRENDY_Temp' &
                     .OR.  TRIM(CRU%Run) .EQ. 'S0_TRENDY_CO2_Temp' ) THEN
                   IF (iVar.EQ.6 .OR. iVar.EQ.7) THEN
                      NextMetYear = CRU%CYEAR
                   ELSE
                      NextMetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
                   ENDIF
                ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY' ) THEN
                   NextMetYear = CRU%CYEAR
                ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY_precip0' ) THEN
                   IF (iVar.EQ.1) THEN
                      ! special for baseline precip
                      NextMetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
                   ELSE
                      NextMetYear = CRU%CYEAR
                   ENDIF
                ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY_precip' ) THEN
                   IF (iVar.NE.1) THEN
                      ! special for baseline non-precip
                      NextMetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
                   ELSE
                      NextMetYear = CRU%CYEAR
                   ENDIF
                ENDIF

                ! Open the Tmin file for next year in preparation for reading Jan 1st...
                CALL CRU_GET_FILENAME( CRU, NextMetYear, Tmin, filename )
                ErrStatus = NF90_OPEN(TRIM(filename), NF90_NOWRITE, fid)
                CALL HANDLE_ERR(ErrStatus, "Opening CRU file "//filename )
                ErrStatus = NF90_INQ_VARID(fid,TRIM(CRU%VAR_NAME(iVar)), vid)
                CALL HANDLE_ERR(ErrStatus, "Inquiring CRU var "//filename )


                ! If DirectRead is specified (small domain) pull the points directly from the file.
                IF ( CRU%DirectRead ) THEN

                   DO k = 1, CRU%mland
                      ErrStatus = NF90_GET_VAR(fid, vid, CRU%MET(ii)%METVALS(k), &  ! Values to ii (nextTmin)
                           start=(/land_x(k), land_y(k),t/) )
                      CALL HANDLE_ERR(ErrStatus, "Reading directly from "//filename )
                   END DO

                   ! Else not DirectRead: Read full grid into temp array and extract domain from that.
                ELSE

                   ErrStatus = NF90_GET_VAR(fid, vid, tmparr, &
                        start=(/1,1,t/),count=(/xds,yds,1/) )
                   CALL HANDLE_ERR(ErrStatus, "Reading from "//filename )
                   DO k = 1, CRU%mland
                      CRU%MET(ii)%METVALS(k) = tmparr( land_x(k), land_y(k) )   ! Values to ii (nextTmin)
                   END DO

                ENDIF  ! End of If DirectRead

                ! Close the next year's Tmin met file
                ErrStatus = NF90_CLOSE(fid)
                CALL HANDLE_ERR(ErrStatus, "Closing CRU file "//filename)

             ENDIF ! End of If LastYearOfMet

          END IF  ! End of If LastDayOfYear

          ! If this is not the first call or last day of the year (i.e. somewhere in the middle)
          ! Just read the next timestep (tplus1) of Tmin data (iVar) into the nextTmin (ii) variable.

          IF ((.NOT. CALL1) .AND. (.NOT. LastDayOfYear)) THEN
             IF ( CRU%DirectRead ) THEN

                DO k = 1, CRU%mland
                   ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), CRU%MET(ii)%METVALS(k), &
                        start=(/land_x(k),land_y(k),tplus1/) )
                   CALL HANDLE_ERR(ErrStatus, "Reading directly from "//CRU%MetFile(iVar) )
                END DO

             ELSE

                ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
                     start=(/1,1,tplus1/),count=(/xds,yds,1/) )
                CALL HANDLE_ERR(ErrStatus, "Reading from "//CRU%MetFile(iVar) )
                DO k = 1, CRU%mland
                   CRU%MET(ii)%METVALS(k) = tmparr( land_x(k), land_y(k) )
                END DO
             ENDIF
          END IF

          !--------------------------------------------------------------------------------------------------

       CASE (TMax)

          ! For Tmax, we need the previous day's Tmax for the Cesaraccio et al subdiurnal
          ! temperature interpolation algorithm.

          IF ( CALL1 ) THEN   ! On the first day of the run...

             ! Assign the met variable index to be prevTmax. iVar will still refer to Tmax.
             ii = prevTmax

             ! If the very first MetYear is more than 1901, there is a previous year of Tmax data available.
             ! such as during spinups, when the starting year is 1951.
             IF ( MetYear .GT. 1901 ) THEN

                ! Open the previous year's Tmax file (MetYear-1)
                CALL CRU_GET_FILENAME( CRU, MetYear-1, iVar, filename )
                ErrStatus = NF90_OPEN(TRIM(filename), NF90_NOWRITE, fid)
                CALL HANDLE_ERR(ErrStatus, "Opening CRU file "//filename )

                ! Obtain the size of the time dimension (tds), which locates the data for Dec 31st in the file.
                ErrStatus = NF90_INQ_DIMID(fid,'time',tid)
                ErrStatus = NF90_INQUIRE_DIMENSION(fid,tid,len=tds)
                CALL HANDLE_ERR(ErrStatus, "Inquiring 'time'"//TRIM(filename))

                ! Obtain the variable ID of the Tmax variable.
                ErrStatus = NF90_INQ_VARID(fid,TRIM(CRU%VAR_NAME(iVar)), vid)
                CALL HANDLE_ERR(ErrStatus, "Inquiring CRU var "//filename )

                ! If DirectRead is specified (small domain) pull the data for Dec 31st (tds) directly from the file.
                ! Otherwise, read the whole Dec 31st grid into tmparr and pack the land points from there.
                IF ( CRU%DirectRead ) THEN
                   !print *,"about to read last value of MetYear-1 ", MetYear - 1, trim(filename)
                   DO k = 1, CRU%mland
                      ErrStatus = NF90_GET_VAR(fid, vid, CRU%MET(ii)%METVALS(k), &
                           start=(/land_x(k), land_y(k),tds/) )
                      CALL HANDLE_ERR(ErrStatus, "Reading from "//filename )
                   END DO

                ELSE

                   ErrStatus = NF90_GET_VAR(fid, vid, tmparr, &
                        start=(/1,1,tds/),count=(/xds,yds,1/) )
                   CALL HANDLE_ERR(ErrStatus, "Reading from "//filename )

                   DO k = 1, CRU%mland
                      CRU%MET(ii)%METVALS(k) = tmparr( land_x(k), land_y(k) )
                   END DO

                ENDIF  ! End of If DirectRead

                ErrStatus = NF90_CLOSE(fid)
                CALL HANDLE_ERR(ErrStatus, "Closing CRU file "//filename)

                ! Having read the last day of the previous year's Tmax into prevTmax, and closed the file,
                ! now read the first day of Tmax from the current year's file (which is already open).

                t  = CRU%CTSTEP  ! CRU%CTSTEP here should be 1 (?)

                ! Read the current points directly into the vector for small domains, or read the whole grid into tmparr
                ! and extract them from there. Read the current timestep value into the iVar variable (Tmax).
                IF ( CRU%DirectRead ) THEN

                   DO k = 1, CRU%mland
                      ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), CRU%MET(iVar)%METVALS(k), &
                           start=(/land_x(k),land_y(k),t/) )
                      CALL HANDLE_ERR(ErrStatus, "Reading directly from "//CRU%MetFile(iVar) )
                   END DO

                ELSE

                   ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
                        start=(/1,1,t/),count=(/xds,yds,1/) )
                   CALL HANDLE_ERR(ErrStatus, "Reading from "//CRU%MetFile(iVar) )
                   DO k = 1, CRU%mland
                      CRU%MET(iVar)%METVALS(k) = tmparr( land_x(k), land_y(k) )
                   END DO
                ENDIF

             ELSE ! MetYear = 1901

                ! Where MetYear is 1901, there is no previous Tmax value available, so we read the
                ! Jan 1 value into Tmax, then assign it to prevTmax as well

                ! Open the 1901 Tmax file (MetYear)
                CALL CRU_GET_FILENAME( CRU, 1901, iVar, filename )
                ErrStatus = NF90_OPEN(TRIM(filename), NF90_NOWRITE, fid)
                CALL HANDLE_ERR(ErrStatus, "Opening CRU file "//filename )

                t  = CRU%CTSTEP  ! CRU%CTSTEP here should be 1 (?)

                ! Read the current points directly into the vector for small domains, or read the whole grid into tmparr
                ! and extract them from there. Read the current timestep value into the iVar variable (Tmax).
                IF ( CRU%DirectRead ) THEN

                   DO k = 1, CRU%mland
                      ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), CRU%MET(iVar)%METVALS(k), &
                           start=(/land_x(k),land_y(k),t/) )
                      CALL HANDLE_ERR(ErrStatus, "Reading directly from "//CRU%MetFile(iVar) )
                   END DO

                ELSE

                   ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
                        start=(/1,1,t/),count=(/xds,yds,1/) )
                   CALL HANDLE_ERR(ErrStatus, "Reading from "//CRU%MetFile(iVar) )
                   DO k = 1, CRU%mland
                      CRU%MET(iVar)%METVALS(k) = tmparr( land_x(k), land_y(k) )
                   END DO
                ENDIF

                ! Now set the value for prevTmax (ii) equal to the current value for Tmax.
                CRU%MET(ii)%METVALS(:) = CRU%MET(iVar)%METVALS(:)

             ENDIF  ! End of If MetYear > 1901

          ELSE  ! Not CALL1

             ! This is not the first call, so the only thing remaining to do for Tmax is read the new value.
             ! (Current Tmax was assigned to prevTmax at the top of the routine).

             ! Set the variable index to iVar (Tmax). The timestep is the current timestep.
             t  = CRU%CTSTEP

             ! Read the current points directly into the vector for small domains, or read the whole grid into tmparr
             ! and extract them from there.
             IF ( CRU%DirectRead ) THEN

                DO k = 1, CRU%mland
                   ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), CRU%MET(iVar)%METVALS(k), &
                        start=(/land_x(k),land_y(k),t/) )
                   CALL HANDLE_ERR(ErrStatus, "Reading directly from "//CRU%MetFile(iVar) )
                END DO

             ELSE

                ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
                     start=(/1,1,t/),count=(/xds,yds,1/) )
                CALL HANDLE_ERR(ErrStatus, "Reading from "//CRU%MetFile(iVar) )
                DO k = 1, CRU%mland
                   CRU%MET(iVar)%METVALS(k) = tmparr( land_x(k), land_y(k) )
                END DO
             ENDIF

          END IF   ! End of If CALL1

          !--------------------------------------------------------------------------------------------------

       CASE DEFAULT   ! All variables other than Tmin or Tmax...

          ! iVar is not Tmin or Tmax so the variable index and timestep index are unchanged.
          ii = iVar
          t  = CRU%CTSTEP

          ! Standard read of the current variable, for the current timestep:
          ! Directly read the current points into the met vector (more efficient for small domains),
          ! or read the whole grid into tmparr and extract them from there.
          IF ( CRU%DirectRead ) THEN

             DO k = 1, CRU%mland
                ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), CRU%MET(ii)%METVALS(k), &
                     start=(/land_x(k),land_y(k),t/) )
                CALL HANDLE_ERR(ErrStatus, "Reading directly from "//CRU%MetFile(iVar) )
             END DO

          ELSE

             ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
                  start=(/1,1,t/),count=(/xds,yds,1/) )
             CALL HANDLE_ERR(ErrStatus, "Reading from "//CRU%MetFile(iVar) )
             DO k = 1, CRU%mland
                CRU%MET(ii)%METVALS(k) = tmparr( land_x(k), land_y(k) )
             END DO
          ENDIF

       END SELECT

    END DO  ! End of loop through all met variables

    ! Convert pressure Pa -> hPa

    CRU%MET(pres)%METVALS(:) = CRU%MET(pres)%METVALS(:) / 100.

    !print *, 'CRU%CTSTEP ', CRU%CTSTEP

    ! Increment the internal timestep counter
    CRU%CTSTEP = CRU%CTSTEP + 1

    ! CALL1 can only happen once!
    IF (CALL1) CALL1 = .FALSE.

    !  print *, 'CRU%MET(tmin)%METVALS(:)'
    !  print *, CRU%MET(tmin)%METVALS(:)
    !  print *, 'CRU%MET(nexttmin)%METVALS(:)'
    !  print *, CRU%MET(nexttmin)%METVALS(:)
    !  print *, 'CRU%MET(prevtmax)%METVALS(:)'
    !  print *, CRU%MET(prevtmax)%METVALS(:)
    !  print *, 'CRU%MET(tmax)%METVALS(:)'
    !  print *, CRU%MET(tmax)%METVALS(:)
    !  print *, 'CRU%MET(pres)%METVALS(:)'
    !  print *, CRU%MET(pres)%METVALS(:)
    !  print *, 'CRU%MET(rain)%METVALS(:)'
    !  print *, CRU%MET(rain)%METVALS(:)

  END SUBROUTINE CRU_GET_DAILY_MET

  !**************************************************************************************************

  SUBROUTINE CRU_GET_SUBDIURNAL_MET(CRU, MET, CurYear, ktau, kend, LastYearOfMet )

    ! Obtain one day of CRU-NCEP meteorology, subdiurnalise it using a weather
    ! and return the result to the CABLE driver.

    USE cable_def_types_mod,   ONLY: MET_TYPE
    USE cable_IO_vars_module,  ONLY: LANDPT, latitude
  USE casa_ncdf_module, ONLY: DOYSOD2YMDHMS
    USE cable_weathergenerator,ONLY: WEATHER_GENERATOR_TYPE, WGEN_INIT, &
         WGEN_DAILY_CONSTANTS, WGEN_SUBDIURNAL_MET
    USE cable_checks_module,   ONLY: rh_sh

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: CurYear, ktau, kend
    LOGICAL, INTENT(IN)  :: LastYearOfMet

    TYPE(CRU_TYPE) :: CRU

    ! Define MET the CABLE version, different from the MET defined and used within the CRU variable.
    ! The structure of MET_TYPE is defined in cable_define_types.F90
    TYPE(MET_TYPE) :: MET

    ! Local variables
    LOGICAL   :: newday, LastDayOfYear  ! Flags for occurence of a new day (0 hrs) and the last day of the year.
    INTEGER   :: iland                  ! Loop counter through 'land' cells (cells in the spatial domain)
    INTEGER   :: itimestep              ! Loop counter through subdiurnal timesteps in a day
    INTEGER   :: imetvar                ! Loop counter through met variables
    INTEGER   :: dM, dD                 ! Met date as year, month, and day returned from DOYSOD2YMDHMS
    INTEGER   :: is, ie                 ! Starting and ending vegetation type per land cell
    REAL      :: dt                     ! Timestep in seconds
    REAL      :: CO2air                 ! CO2 concentration in ppm
    REAL      :: etime
    CHARACTER :: LandMaskFile*200       ! Name of the land mask file

    TYPE(WEATHER_GENERATOR_TYPE), SAVE :: WG
    LOGICAL,                      SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine

    ! Purely for readability...
    dt = CRU%DTsecs

    ! On first step read and check CRU settings and read land-mask
    IF ( CALL1 ) CALL WGEN_INIT( WG, CRU%mland, latitude, dt )

    ! Pass time-step information to CRU
    CRU%CYEAR = CurYear


    CRU%ktau  = ktau     ! ktau is the current timestep in the year.

!!!!  this only works with CANBERRA cable_driver, as ktau    !!!!
!!!!  restarts on Jan 1                                      !!!!
    ! Based on the ktau timestep, calculate date and time information (the same for the whole spatial dimension.)
    met%hod (:) = REAL(MOD( (ktau-1) * NINT(dt), INT(SecDay)) ) / 3600.  ! Hour of the day
    met%doy (:) = INT(REAL(ktau-1) * dt / SecDay ) + 1                   ! Day of Year = days since 0 hr 1st Jan
    met%year(:) = CurYear                                                ! Current year

    ! Using the day-of-year and seconds-of-day calculate the month and day-of-month, using the time information
    ! for the first land cell only (because they will be the same across the domain).
    !                            In      In         In        Out       Out       Optional Out
    ! SUBROUTINE DOYSOD2YMDHMS( Year, Yearday, SecondsOfDay, Month, DayOfMonth, [Hour, Min, Sec])

    CALL DOYSOD2YMDHMS(CurYear, INT(met%doy(1)), INT(met%hod(1)) * 3600, dM, dD)

    met%moy (:) = dM     ! Record the month

    ! It's a new day if the hour of the day is zero.
    newday = ( met%hod(landpt(1)%cstart).EQ. 0 )

    ! Beginning-of-year accounting
    IF (ktau .EQ. 1) THEN  ! ktau is always reset to 1 at the start of the year.

       ! Read a new annual CO2 value and convert it from ppm to mol/mol
       CALL GET_CRU_CO2( CRU, CO2air )
       met%ca(:) = CO2air / 1.e+6  !

       CALL GET_CRU_Ndep( CRU )
       DO iland = 1, CRU%mland
          met%Ndep(landpt(iland)%cstart:landpt(iland)%cend) = &
               CRU%NdepVALS(iland)*86400000.  ! kg/m2/s > g/m2/d (1000.*3600.*24.)
       END DO

       ! Open a new annual CRU-NCEP met file.
       CALL OPEN_CRU_MET( CRU )

    ENDIF

    ! %%%%%% PRB to add his own comments from here down for this routine.
    ! Now get the Met data for this day

    IF ( newday ) THEN

       !print *, CRU%CTSTEP, ktau, kend
       !   CALL CPU_TIME(etime)
       !   PRINT *, 'b4 daily ', etime, ' seconds needed '

       LastDayOfYear = (ktau .EQ. kend-((SecDay/dt)-1))
       CALL CRU_GET_DAILY_MET( CRU, LastDayOfYear, LastYearOfMet )
       !   CALL CPU_TIME(etime)
       !   PRINT *, 'after daily ', etime, ' seconds needed '
       !   STOP

       ! Air pressure assumed to be constant over day
       DO iland = 1, CRU%mland
          met%pmb(landpt(iland)%cstart:landpt(iland)%cend) = CRU%MET(pres)%METVALS(iland) !CLN interpolation??
       END DO

       ! Convert wind from u and v components to wind speed by Pythagorean Theorem
       WG%WindDay        = SQRT( (CRU%MET(uwind)%METVALS * CRU%MET(uwind)%METVALS) +  &
            (CRU%MET(vwind)%METVALS * CRU%MET(vwind)%METVALS) )
       ! Convert all temperatures from K to C
       WG%TempMinDay     = CRU%MET(  Tmin  )%METVALS - 273.15
       WG%TempMaxDay     = CRU%MET(  Tmax  )%METVALS - 273.15
       WG%TempMinDayNext = CRU%MET(NextTmin)%METVALS - 273.15
       WG%TempMaxDayPrev = CRU%MET(PrevTmax)%METVALS - 273.15
       ! Convert solar radiation from J /m2/s to MJ/m2/d
       WG%SolarMJDay     = CRU%MET(  swdn  )%METVALS * 1.e-6 * SecDay ! ->[MJ/m2/d]
       ! Convert precip from mm to m/day
       WG%PrecipDay      = CRU%MET(  rain  )%METVALS  / 1000. ! ->[m/d]
       WG%SnowDay        = 0.0

       CALL WGEN_DAILY_CONSTANTS( WG, CRU%mland, INT(met%doy(1))+1 )

       ! To get the diurnal cycle for lwdn get whole day and scale with
       ! LWDN from file later

       CRU%AVG_LWDN(:) = 0.
       DO itimestep = 1, NINT(SecDay/dt)
          CALL WGEN_SUBDIURNAL_MET( WG, CRU%mland, itimestep-1 )
          CRU%AVG_LWDN = CRU%AVG_LWDN + WG%PhiLD
       END DO
       CRU%AVG_LWDN = CRU%AVG_LWDN / (SecDay/dt)
    END IF ! End of If newday

    ! Decision has been made, that first tstep of the day is at 0:01 am

    CALL WGEN_SUBDIURNAL_MET( WG, CRU%mland, NINT(met%hod(1)*3600./dt) )

    ! assign to cable variables

    ! strangely, met% is not save over Wait all in MPI...!
    !  met%pmb = CRU%MET(pres)%METVALS !CLN interpolation??

    ! Assign weather-generated data, or daily values, as required to CABLE variables.
    !
    DO iland = 1, CRU%mland
       is = landpt(iland)%cstart
       ie = landpt(iland)%cend

       met%precip    (is:ie)   = WG%Precip(iland)


       ! Cable's swdown is split into two components, visible and nir, which
       ! get half of the CRU-NCEP swdown each.
       met%fsd       (is:ie,1) = WG%PhiSD(iland) * 0.5    ! Visible
       met%fsd       (is:ie,2) = WG%PhiSD(iland) * 0.5  ! NIR

       ! Convert C to K for cable's tk
       met%tk        (is:ie)   = WG%Temp(iland) + 273.15

       met%ua        (is:ie)   = WG%Wind(iland)
       met%coszen    (is:ie)   = WG%coszen(iland)

       ! For longwave down, scale the diurnal series returned by the weather generator (WG%PhiLD(iland))
       ! to the daily value from CRU-NCEP.
       met%fld       (is:ie)   = CRU%MET(lwdn)%METVALS(iland) * WG%PhiLD(iland) / CRU%AVG_LWDN(iland)

       ! Specific humidity (qair g/g) was not sent to the weather generator. Here we assign the
       ! daily value to the whole diurnal cycle
       met%qv        (is:ie)   = CRU%MET(qair)%METVALS(iland)

       ! calculate snowfall based on total precip and air T
       !(ref Jin et al. Table II, Hyd Proc, 1999)

!!$    if (WG%Temp(iland) > 2.5) then
!!$       met%precip_sn(is:ie) = 0.0
!!$    elseif ((WG%Temp(iland) <= 2.5) .and. (WG%Temp(iland) > 2.0)) then
!!$       met%precip_sn(is:ie) = 0.6* met%precip(is:ie)
!!$    elseif ((WG%Temp(iland) <= 2.0) .and. (WG%Temp(iland) > 0.0)) then
!!$       met%precip_sn(is:ie) = (1.0 - (54.62 - 0.2 *(WG%Temp(iland) + 273.15)))* met%precip(is:ie) ! this facr can be > 1 !!!
!!$    elseif (WG%Temp(iland) <= 0.0) then
!!$       met%precip_sn(is:ie) = met%precip(is:ie)
!!$    endif

       IF (WG%Temp(iland) <= 0.0) THEN
          met%precip_sn(is:ie) = met%precip(is:ie)
       ELSE
          met%precip_sn(is:ie) = 0.0
       ENDIF
       !met%precip(is:ie) = met%precip(is:ie) -  met%precip_sn(is:ie)



    END DO

    ! initialise within canopy air temp
    met%tvair     = met%tk
    met%tvrad     = met%tk

    ! If this is the end of the year or the end of the met, close the current met files.
    !print *, "ktau, kend, LastYearOfMet as close test:", ktau, kend
    IF (ktau .EQ. kend) THEN
       DO imetvar=1, CRU%NMET
          !print *, 'Close CRU%MetFile(imetvar)', CRU%MetFile(imetvar)
          ErrStatus = NF90_CLOSE(CRU%F_ID(imetvar))
          CALL HANDLE_ERR(ErrStatus, "Closing CRU file"//CRU%MetFile(imetvar))
       END DO
    END IF

    ! CALL1 is over...
    CALL1 = .FALSE.

  END SUBROUTINE CRU_GET_SUBDIURNAL_MET

END MODULE CABLE_CRU
