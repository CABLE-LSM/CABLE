module cable_cru

  use netcdf, only: &               ! Access to netcdf routines
       nf90_open, nf90_nowrite, nf90_inq_dimid, nf90_inquire_dimension, &
       nf90_inq_varid, nf90_get_var, nf90_close
  use cable_common_module, only: &  ! Selected cable_common.f90 routines:
       handle_err,  &               ! Print error status info returned by netcdf file operations
       get_unit, &                  ! Finds an unused unit number for file opens
       cable_user                   ! used for diffuse radiation switch
  use cable_IO_vars_module, only: & ! Selected cable_iovars.F90 variables:
       logn,            &           ! Log file unit number
       land_x, land_y,  &           ! Col (x) & row (y) indices of each land point
                                    ! in land mask (dimension mland)
       exists                       ! Only for exists%Snowf, which we will set to .FALSE.
                                    ! because there is no snow in CRU-NCEP. Setting this ensures
                                    ! snow will be determined in CABLE from temperature.
  implicit none

  integer, parameter :: sp = kind(1.0)

  ! Define a type for CRU-NCEP information, and the subtype METVALS

  type cru_met_type
     ! define a spatial vector of meteorology for one timestep
     real, dimension(:), allocatable :: metvals
  end type cru_met_type

  type cru_type
     integer :: mland              ! Number of land cells
     integer :: NMET               ! Number of met variable types (rain, lwdn etc) NOT INCLUDING prevTmax and nextTmin
     integer :: xdimsize, ydimsize ! Landmask grid size dimensions (x=cols, y=rows)
     integer :: tdimsize           ! Time dimension of metfiles (met data timesteps per annual file)
     integer :: CYEAR              ! Current Run Year, Same As Curyear, Not Necessarily The Same As Metyear
     integer :: Metstart           ! First Year Of Met
     integer :: Metend             ! Last Year Of Met
     integer :: Ctstep             ! Current Met Data Timestep (1 To Tdimsize, I.E. 365 For Cru-Ncep Annual Daily Files)
     integer :: Dtsecs             ! Model Timestep In Seconds, Converted From Namelist Value In Hours
     integer :: Ktau               ! Current model timestep, reset at the start of a new year of met
     integer :: metrecyc=20        ! number of years for the met recycling
     integer, dimension(10) :: f_id, v_id ! NetCDF object id's for files and variables (NetCDF bookkeeping stuff)
     ! Avg of one day's diurnal cycle of lwdn calculated by Swinbank. AVG_LWDN
     ! is used to rescale the diurnal cycle to match the day's CRUNCEP lwdn. (dim=mland)
     real, dimension(:), allocatable :: avg_lwdn
     ! Global annual CO2 values (dim is the number of years of data, or 1 if time-invariant)
     real, dimension(:), allocatable :: co2vals
     logical :: DirectRead   ! Flag to do with reading small numbers of points efficiently. Set true for small numbers of points
     logical :: LeapYears    ! Flag for whether leaps years occur, required by CABLE. Always false for CRUNCEP (no Feb 29th)
     logical :: ReadDiffFrac ! Is fraction of diffuse radiation read in from met file (TRUE) or calculated in rad routine (FALSE)?
     logical, dimension(:,:), allocatable :: LandMask ! Logical landmask, true for land, false for non-land
     !
     character(len=30)  :: run            ! Where run type is      : "S0_TRENDY", "S1_TRENDY", "S2_TRENDY"
     character(len=15)  :: CO2            ! CO2 takes value        : "static1860", "1860_1900", "1901_2015"
     character(len=15)  :: ndep           ! Ndep takes value        : "static1860", "1860_1900", "1901_2015"
     character(len=15)  :: forcing        ! Met Forcing takes value: "spinup",        "spinup", "1901_2015"
     character(len=200) :: BasePath       ! Full path for the location of data used for CRU runs "/x/y"
     character(len=200) :: MetPath        ! Full path for the location of the met files "/x/y"
     character(len=50)  :: MetVersion     ! Met Forcing Version (currently CRUJRA_YEAR and VERIFY_2021)
     character(len=200) :: LandMaskFile   ! Land mask filename, without path
     ! Netcdf variable 'Name' for each type of met (dim=# of met vars). Note: Name, not 'Title'
     character(len=30),  dimension(10)  :: var_name
     ! Met file names incl metpath, constructed in CRU_GET_FILENAME (dim=# of met vars)
     character(len=200), dimension(10)  :: MetFile
     ! Met data vectors (METVALS) for one timestep, dim=# of met vars + 2 for prev Tmax and next Tmin
     type(cru_met_type), dimension(12) :: met
     real, dimension(:), allocatable :: ndepvals
     integer :: ndepf_id, ndepv_id
     integer :: ndep_ctstep   ! counter for Ndep in input file
  end type cru_type

  ! TYPE(CRU_TYPE) :: CRU  ! Define the variable CRU, of type CRU_TYPE

  ! Define local parameter names representing the position of each met var within variable MET.
  ! prevTmax and nextTmin are special cases of Tmax and Tmin that do not count as extra met variables per se.
  integer, private, parameter :: &
       rain     =  1, &
       lwdn     =  2, &
       swdn     =  3, &
       pres     =  4, &
       qair     =  5, &
       tmax     =  6, &
       tmin     =  7, &
       uwind    =  8, &
       vwind    =  9, &
       fdiff    = 10, &
       prevTmax = 11, &
       nextTmin = 12

  ! Error status of various operations (mostly netcdf-related). Typically 0 means ok, > 0 means unexpected condition.
  integer, private :: ErrStatus

  real, private, parameter :: SecDay = 86400. ! Number of seconds in a day

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------

  subroutine CRU_INIT(CRU)

    ! Initialise the contents of the CRU defined type collection, from the CRU namelist file
    ! and by obtaining dimensions from the landmask

    use cable_IO_vars_module, only: &
         latitude, longitude, & ! (R) Lat and long of landcells only (?)
         nmetpatches,         & ! (I) Size of patch dimension in met file, if it exists
         mask,                & ! (I) Land/sea mask (1,0)
         metGrid,             & ! (C4) Either 'land' or 'mask' for whether the data are packed or not (?)
         sdoy, smoy, syear,   & ! (I) Start time day of year, month, year
         shod,                & ! (R) Start time hour of day
         xdimsize, ydimsize,  & ! (I) Size of grid dimensions
         lat_all, lon_all       ! (R) Grids with the lat or lon of each cell
                                ! (i.e. repetition along rows/cols), for CABLE.
    use cable_def_types_mod,  only: mland  ! (I) Number of land cells
#ifdef __MPI__
    use mpi,                  only: MPI_Abort
#endif

    implicit none

    type(CRU_TYPE), intent(inout) :: CRU

    integer :: ErrStatus  ! Error status returned by nc routines (zero=ok, non-zero=error)
    integer :: nmlunit    ! Unit number for reading namelist file
    integer :: FID        ! NetCDF id for the landmask file
    integer :: latID, lonID  ! NetCDF ids for dimensions in the landmask file
    integer :: landID     ! NetCDF id for the landmask variable in the landmask file
    integer :: landcnt    ! Manually incremented counter for the number of land cells
    integer :: xcol, yrow ! Column and row position in the data file grids
    integer :: imetvar    ! loop counter through met variables

    ! Temporary local names for CRU% variables as they are read from the
    ! namelist file. Note that CRU%CO2 and CRU%Forcing are assigned based on the
    ! value of Run, not read as options from the namelist file.
    logical            :: DirectRead = .false.
    logical            :: ReadDiffFrac = .false.
    character(len=30)  :: Run
    character(len=200) :: BasePath
    character(len=200) :: MetPath
    character(len=50)  :: MetVersion
    character(len=200) :: LandMaskFile
    ! CABLE timestep (hrs), converted immediately to integer seconds for CRU%DTsecs
    real               :: DThrs
    ! Lat/long values for each grid rows/cols from landmask.
    real,    dimension(:),   allocatable :: CRU_lats, CRU_lons
    integer, dimension(:,:), allocatable :: landmask

    ! Flag for errors
    logical :: ERR = .false.
#ifdef __MPI__
    integer :: ierr
#endif

    namelist /crunml/ BasePath, MetPath, MetVersion, ReadDiffFrac, LandMaskFile, Run, DThrs, DirectRead

    ! Read CRU namelist settings
    call get_unit(nmlunit)  ! CABLE routine finds spare unit number
    open(nmlunit, file="cru.nml", status='old', action='read')
    read(nmlunit, nml=crunml)
    close(nmlunit)

    ! Assign namelist settings to corresponding CRU defined-type elements
    CRU%BasePath     = BasePath
    CRU%MetPath      = MetPath
    CRU%MetVersion   = trim(MetVersion)
    CRU%ReadDiffFrac = ReadDiffFrac
    CRU%LandMaskFile = trim(LandMaskFile)
    CRU%Run          = Run
    CRU%DTsecs       = int(DThrs * 3600.)  ! in seconds
    CRU%DirectRead   = DirectRead

    ! diffuse fraction not available for all Metversions
    if ((CRU%ReadDiffFrac == .true.) .and. &
      ((trim(CRU%MetVersion) /= "CRUJRA_2022") .and. (trim(CRU%MetVersion) /= "CRUJRA_2023"))) then
       write(*,'(a)') "Diffuse Fraction only available for CRUJRA_2022 and CRUJRA_2023!"
       write(logn,*)  "Diffuse Fraction only available for CRUJRA_2022 and CRUJRA_2023!"
    endif

    ! Assign Forcing and CO2 labels based only on the value of CRU%Run
    select case(trim(CRU%Run))
    case("S0_TRENDY")
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep    = "static1860"
       write(*,'(a)') "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
       write(logn,*)  "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
    case("S0_TRENDY_CO2")
       CRU%Forcing = "spinup"
       CRU%CO2     = "1901_2015"
       CRU%Ndep    = "static1860"
       write(*,'(a)') "Run = 'S0_CO2': Therefore Forcing = 'S0_CO2', CO2 = '1860_2015'"
       write(logn,*)  "Run = 'S0_CO2': Therefore Forcing = 'S0_CO2', CO2 = '1860_2015'"
    case("S0_TRENDY_Ndep")
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep    = "1901_2015"
       write(*,'(a)') "Run = 'S0_Ndep': Therefore Forcing = 'S0_Ndep', Ndep = '1860_2015'"
       write(logn,*)  "Run = 'S0_Ndep': Therefore Forcing = 'S0_Ndep', Ndep = '1860_2015'"
    case("S0_TRENDY_Precip")
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep    = "static1860"
       write(*,'(a)') "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
       write(logn,*)  "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
    case("S0_TRENDY_Temp")
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep    = "static1860"
       write(*,'(a)') "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
       write(logn,*)  "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
    case("S0_TRENDY_Temp_Precip")
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1860"
       CRU%Ndep    = "static1860"
       write(*,'(a)') "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
       write(logn,*)  "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1860'"
    case("S0_TRENDY_CO2_Temp")
       CRU%Forcing = "spinup"
       CRU%CO2     = "1901_2015"
       CRU%Ndep    = "static1860"
       write(*,'(a)') "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
       write(logn,*)  "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
    case("S0_TRENDY_CO2_Precip")
       CRU%Forcing = "spinup"
       CRU%CO2     = "1901_2015"
       CRU%Ndep    = "static1860"
       write(*,'(a)') "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
       write(logn,*)  "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
    case("S0_TRENDY_CO2_Temp_Precip")
       CRU%Forcing = "spinup"
       CRU%CO2     = "1901_2015"
       CRU%Ndep    = "static1860"
       write(*,'(a)') "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
       write(logn,*)  "Run = 'spinup': Therefore Forcing = 'spinup', CO2 = '1860_2015'"
    case("S1_TRENDY")
       CRU%Forcing = "spinup"
       CRU%CO2     = "1860_1900"
       CRU%Ndep    = "1860_1900"
       write(*,'(a)') "Run = 'S1_TRENDY': Therefore Forcing = 'spinup', CO2 = '1860_1900'"
       write(logn,*)  "Run = 'S1_TRENDY': Therefore Forcing = 'spinup', CO2 = '1860_1900'"
    case("S2_TRENDY")
       CRU%Forcing = "1901_2015"
       CRU%CO2     = "1901_2015"
       CRU%Ndep    = "1901_2015"
       write(*,'(a)') "Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
       write(logn,*)  "Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
    case("S2_TRENDY_precip")
       CRU%Forcing = "1901_2015"
       CRU%CO2     = "1901_2015"
       CRU%Ndep    = "1901_2015"
       write(*,'(a)') "Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
       write(logn,*)  "Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
    case("S2_TRENDY_precip0")
       CRU%Forcing = "1901_2015"
       CRU%CO2     = "1901_2015"
       CRU%Ndep    = "1901_2015"
       write(*,'(a)') "Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
       write(logn,*)  "Run = 'S2_TRENDY': Therefore Forcing = 'spinup', CO2 = '1901_2015'"
    case default
       write(*,'(a)') "Wrong CRU%Run: "//trim(CRU%Run)
       write(*,'(a)') "Use: S0_TRENDY, S1_TRENDY, or S2_TRENDY!"
       write(logn,*)  "Wrong CRU%Run: ",trim(CRU%Run)
       write(logn,*)  "Use: S0_TRENDY, S1_TRENDY, or S2_TRENDY!"
       ERR = .true.
    end select

    ! Print settings
    write(*,'(a)') "========================================= CRU ============"
    write(*,'(a)') "CRU settings chosen:"
    write(*,'(a)') "  BasePath: "//trim(CRU%BasePath)
    write(*,'(a)') "  LandMask: "//trim(CRU%LandMaskFile)
    write(*,'(a)') "  Run               : "//trim(CRU%Run)
    write(*,'(a)') "  Forcing (assigned): "//trim(CRU%Forcing)
    write(*,'(a)') "  CO2     (assigned): "//trim(CRU%CO2)
    write(*,'(a)') "  Ndep    (assigned): "//trim(CRU%Ndep)
    write(*,'(a,i5)') "  DT(secs): ", CRU%DTsecs
    write(logn,*) "========================================= CRU ============"
    write(logn,*) "CRU settings chosen:"
    write(logn,*) " BasePath: ",trim(CRU%BasePath)
    write(logn,*) " LandMask: ",trim(CRU%LandMaskFile)
    write(logn,*) " Run               : ",trim(CRU%Run)
    write(logn,*) " Forcing (assigned): ",trim(CRU%Forcing)
    write(logn,*) " CO2     (assigned): ",trim(CRU%CO2)
    write(logn,*) " Ndep    (assigned): ",trim(CRU%Ndep)
    write(logn,*) " DT(secs): ",CRU%DTsecs
    if (CRU%ReadDiffFrac) then
       write(*,'(a)') "  Reading diffuse fraction from met file"
       write(logn,*) " Reading diffuse fraction from met file"
    endif

    ! Error trap for bad namelist.
    if (ERR) then
       write(logn,*) "Invalid settings in CRU_INIT"
       write(*,*) "Invalid settings in CRU_INIT"
#ifdef __MPI__
       call MPI_Abort(0, 5, ierr) ! Do not know comm nor rank here
#else
       stop 5
#endif
    endif

    ! ! If this is a S0_TRENDY run look for met data in the spinup directory instead.
    ! IF (TRIM(CRU%Run) .EQ. "S0_TRENDY") THEN
    !   CRU%MetPath = TRIM(CRU%MetPath)//"/spinup_data"
    ! ENDIF

    ! ! Set variable names to their NetCDF 'Names' (i.e. not their 'Titles')
    ! CRU%NMET = 9
    ! CRU%VAR_NAME(rain)  = "Total_Precipitation"
    ! CRU%VAR_NAME(lwdn)  = "Incoming_Long_Wave_Radiation"
    ! CRU%VAR_NAME(swdn)  = "Incoming_Short_Wave_Radiation"
    ! CRU%VAR_NAME(pres)  = "Pression"
    ! CRU%VAR_NAME(qair)  = "Air_Specific_Humidity"
    ! ! CRU%VAR_NAME(tmax)  = "maximum_6h_air_temperature"
    ! ! CRU%VAR_NAME(tmin)  = "minimum_6h_air_temperature"
    ! CRU%VAR_NAME(tmax)  = "maximum_air_temperature"
    ! CRU%VAR_NAME(tmin)  = "minimum_air_temperature"
    ! CRU%VAR_NAME(uwind) = "U_wind_component"
    ! CRU%VAR_NAME(vwind) = "V_wind_component"

    ! Default
    CRU%NMET = 9
    CRU%VAR_NAME(rain)  = "pre"
    CRU%VAR_NAME(lwdn)  = "dlwrf"
    CRU%VAR_NAME(swdn)  = "dswrf"
    CRU%VAR_NAME(pres)  = "pres"
    CRU%VAR_NAME(qair)  = "spfh"
    CRU%VAR_NAME(tmax)  = "tmax"
    CRU%VAR_NAME(tmin)  = "tmin"
    CRU%VAR_NAME(uwind) = "ugrd"
    CRU%VAR_NAME(vwind) = "vgrd"

    CRU%Metstart = 1901

    if (trim(CRU%MetVersion) == "CRUJRA_2021") then
       CRU%VAR_NAME(swdn) = "tswrf"
    else if (trim(CRU%MetVersion) == "CRUJRA_2022" .or. trim(CRU%MetVersion) == "CRUJRA_2023") then
       CRU%VAR_NAME(swdn) = "tswrf"
       if (CRU%ReadDiffFrac) then
          CRU%NMET = 10
          CRU%VAR_NAME(fdiff) = "fd"
       endif
    else if (trim(CRU%MetVersion) == "VERIFY_2021") then
       CRU%VAR_NAME(rain)  = "Precipalign"
       CRU%VAR_NAME(lwdn)  = "LWdownnoalign"
       CRU%VAR_NAME(swdn)  = "SWdownalign"
       CRU%VAR_NAME(pres)  = "Psurfnoalign"
       CRU%VAR_NAME(qair)  = "Qairnoalign"
       CRU%VAR_NAME(tmax)  = "Tmaxalign"
       CRU%VAR_NAME(tmin)  = "Tminalign"
       CRU%VAR_NAME(uwind) = "Wind_Enoalign"
       CRU%VAR_NAME(vwind) = "Wind_Nnoalign"
    end if

    write(*,'(a)') "========================================= CRU ============"
    write(logn,*)  "========================================= CRU ============"

    ! Now read landmask file
    ! Landmask file into init! Get LAt, LON etc. from there
    ! LMFILE = TRIM(CRU%LandMaskFile)
    write(*,'(a)') 'Opening CRU landmask file: ' // trim(LandMaskFile)
    write(logn,*)  'Opening CRU landmask file: ', trim(LandMaskFile)

    ! Open the land mask file
    ErrStatus = NF90_OPEN(trim(LandMaskFile), NF90_NOWRITE, FID)
    call HANDLE_ERR(ErrStatus, "Opening CRU Land-mask file" // trim(LandMaskFile))

    ! Latitude: Get the dimension ID, find the size of the dimension, assign it to CRU.
    ErrStatus = NF90_INQ_DIMID(FID, 'latitude', latID)
    ErrStatus = NF90_INQUIRE_DIMENSION(FID, latID, len=ydimsize)
    call HANDLE_ERR(ErrStatus, "Inquiring 'lat'" // trim(LandMaskFile))
    CRU%ydimsize = ydimsize

    ! Collect the latitudes into CRU_lats
    allocate(CRU_lats(ydimsize))
    ErrStatus = NF90_INQ_VARID(FID, 'latitude', latID)
    call HANDLE_ERR(ErrStatus, "Inquiring 'latitudes'" // trim(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID, latID, CRU_lats)
    call HANDLE_ERR(ErrStatus, "Reading 'latitudes'" // trim(LandMaskFile))

    ! Longitude: Get the dimension ID, find the size of the dimension, assign it to CRU.
    ErrStatus = NF90_INQ_DIMID(FID, 'longitude', lonID)
    ErrStatus = NF90_INQUIRE_DIMENSION(FID, lonID, len=xdimsize)
    call HANDLE_ERR(ErrStatus, "Inquiring 'lon'" // trim(LandMaskFile))
    CRU%xdimsize = xdimsize

    ! Collect the longitudes into CRU_lons
    allocate(CRU_lons(xdimsize))
    ErrStatus = NF90_INQ_VARID(FID, 'longitude', lonID)
    call HANDLE_ERR(ErrStatus, "Inquiring 'longitudes'" // trim(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID, lonID, CRU_lons)
    call HANDLE_ERR(ErrStatus, "Reading 'longitudes'" // trim(LandMaskFile))
    where (CRU_lons > 180.0)
       CRU_lons = CRU_lons - 360.0
    end where

    ! Allocate the landmask arrays for...
    allocate(CRU%landmask(xdimsize, ydimsize))  ! Passing out to other CRU routines (logical)
    allocate(landmask(xdimsize, ydimsize))      ! Local use in this routine (integer)
    allocate(mask(xdimsize, ydimsize))          ! Use by CABLE

    ! Check that the land mask variable is called "land" in the land mask file,
    ! and read it into local variable landmask
    ErrStatus = NF90_INQ_VARID(FID, 'land', landID)
    call HANDLE_ERR(ErrStatus, "Inquiring 'land' " // trim(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID, landID, landmask)
    call HANDLE_ERR(ErrStatus, "Reading 'land' " // trim(LandMaskFile))

    ! Convert the integer landmask into the logical CRU%landmask
    where (landmask > 0)
       CRU%landmask = .true.
       mask         = 1
    elsewhere
       CRU%landmask = .false.
       mask         = 0
    end where

    ! Count the number of land cells -> mland
    CRU%mland = count(CRU%landmask)

    ! Allocate CABLE land-only vectors for lat/long and row/col values/indices.
    allocate(latitude(CRU%mland), longitude(CRU%mland))
    allocate(land_y(CRU%mland), land_x(CRU%mland))

    ! Allocate vectors for each of the different met quantities, including extra
    ! prev/next temperatures for the Cesarracio temperature calculations in the
    ! weather generator.
    do imetvar=1, CRU%NMET
       allocate(CRU%MET(imetvar)%METVALS(CRU%mland))
    end do
    allocate(CRU%MET(prevTmax)%METVALS(CRU%mland))
    allocate(CRU%MET(nextTmin)%METVALS(CRU%mland))
    ! allocate array for Nitrogen deposition input data
    allocate(CRU%NdepVALS(CRU%mland))

    ! Copy the col/row and lat/long positions of each land cell into the corresponding
    ! land only CABLE vectors. Q: We know mland at this point. Why not use landcnt to confirm
    ! the correct value of mland?
    landcnt=1
    do yrow=1, ydimsize
       do xcol=1, xdimsize
          if (.not. CRU%landmask(xcol, yrow)) cycle   ! Go to next iteration if not a land cell
          ! WRITE(6,FMT='(A15,I5,2(1X,F8.2),2(1x,I3))') "i, lo,la, xcol,yrow", landcnt, &
          !    CRU_lons(xcol),CRU_lats(yrow),xcol, yrow
          land_x(landcnt)    = xcol
          land_y(landcnt)    = yrow
          longitude(landcnt) = CRU_lons(xcol)
          latitude(landcnt)  = CRU_lats(yrow)
          landcnt = landcnt + 1
       end do
    end do

    ! Set global CABLE variables
    metGrid     = "mask"
    allocate(mask(xdimsize, ydimsize))
    mask        = landmask
    mland       = CRU%mland
    nmetpatches = 1
    allocate(lat_all(xdimsize, ydimsize), lon_all(xdimsize, ydimsize))
    do xcol=1, xdimsize
       lat_all(xcol,:) = CRU_lats
    end do
    do yrow=1, ydimsize
       lon_all(:,yrow) = CRU_lons
    end do

    ! CABLE TIME-UNITS needed by load-parameters (only on CABLE_init)
    shod  = 0.
    sdoy  = 1
    smoy  = 1
    syear = CRU%CYEAR

    ! Used to rescale the diurnal cycle from Swinbank calculation to match CRU-NCEP provided value.
    allocate(CRU%AVG_LWDN(mland))

    deallocate(landmask, CRU_lats, CRU_lons)

    ErrStatus = NF90_CLOSE(FID)
    FID = -1
    call HANDLE_ERR(ErrStatus, "Closing mask-file" // trim(LandMaskFile))

    ! set units to -1
    CRU%f_id = -1
    CRU%Ndepf_id = -1

  end subroutine CRU_INIT

  ! ------------------------------------------------------------------

  subroutine CRU_GET_FILENAME(CRU, cyear, par, fn)

    ! Build the filename FN: One annual file of daily met for one met quantity.

    implicit none

    type(cru_type),     intent(in)  :: CRU    ! Information about CRU
    integer,            intent(in)  :: cyear  ! Current year as an integer
    integer,            intent(in)  :: par    ! Index (1-9) of which met quantity will be sought
    character(len=200), intent(out) :: fn     ! Met filename (outgoing)

    character(len=4)   :: cy     ! Character representation of cyear
    character(len=200) :: metp   ! Local repr of met path
    character(len=50)  :: cruver ! cru version as in filename

    ! Create a character version of the year for building that part of the filename.
    write(cy, fmt='(i4)') cyear

    ! Initialise the filename with the met path
    metp = trim(CRU%MetPath)
    fn   = trim(metp)

    select case(trim(CRU%MetVersion))
    case("CRUJRA_2021")
       cruver="crujra.v2.2"
    case("CRUJRA_2022")
       cruver="crujra.v2.3"
    case("CRUJRA_2023")
       cruver="crujra.v2.4"
    case("VERIFY_2021")
       cruver="cru_verify"
    end select

    ! Build the rest of the filename according to the value of par, which references 11 possible
    ! types of met through the parameter names rain, lwdn, etc.

    if (trim(CRU%MetVersion) == "VERIFY_2021") then
       select case(par)
       case(rain)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_Precipalign.nc"
       case(lwdn)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_LWdownnoalign.nc"
       case(swdn)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_SWdownalign.nc"
       case(pres)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_Psurfnoalign.nc"
       case(qair)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_Qairnoalign.nc"
       case(tmax, PrevTmax)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_Tmaxalign.nc"
       case(tmin, NextTmin)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_Tminalign.nc"
       case(uwind)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_Wind_Enoalign.nc"
       case(vwind)
          FN = trim(FN)//"/"//trim(cruver)//"_"//cy//"_daily_Wind_Nnoalign.nc"
       end select
    else
       select case(par)
       case(rain)
          fn = trim(fn)//"/pre/"//trim(cruver)//".5d.pre."//cy//".365d.noc.daytot.1deg.nc"
       case(lwdn)
          fn = trim(fn)//"/dlwrf/"//trim(cruver)//".5d.dlwrf."//cy//".365d.noc.daymean.1deg.nc"
       case(swdn)
          if (trim(CRU%MetVersion) == "CRUJRA_2021") then
             fn = trim(fn)//"/tswrf/tswrf_v10_"//cy//".daymean.1deg.nc"
          else if (trim(CRU%MetVersion) == "CRUJRA_2022") then
             fn = trim(fn)//"/tswrf/tswrf_v11_"//cy//".daymean.1deg.nc"
          else if (trim(CRU%MetVersion) == "CRUJRA_2023") then
             fn = trim(fn)//"/tswrf/tswrf_v12_"//cy//".daymean.1deg.nc"
          else
             fn = trim(fn)//"/dswrf/"//trim(cruver)//".5d.dswrf."//cy//".365d.noc.daymean.1deg.nc"
          endif
       case(pres)
          fn = trim(fn)//"/pres/"//trim(cruver)//".5d.pres."//cy//".365d.noc.daymean.1deg.nc"
       case(qair)
          fn = trim(fn)//"/spfh/"//trim(cruver)//".5d.spfh."//cy//".365d.noc.daymean.1deg.nc"
       case(tmax, PrevTmax)
          fn = trim(fn)//"/tmax/"//trim(cruver)//".5d.tmax."//cy//".365d.noc.daymax.1deg.nc"
       case(tmin, NextTmin)
          fn = trim(fn)//"/tmin/"//trim(cruver)//".5d.tmin."//cy//".365d.noc.daymin.1deg.nc"
       case(uwind)
          fn = trim(fn)//"/ugrd/"//trim(cruver)//".5d.ugrd."//cy//".365d.noc.daymean.1deg.nc"
       case(vwind)
          fn = trim(fn)//"/vgrd/"//trim(cruver)//".5d.vgrd."//cy//".365d.noc.daymean.1deg.nc"
       case(fdiff)
          if (trim(CRU%MetVersion) == "CRUJRA_2022") then
             fn = trim(fn)//"/fd/fd_v11_"//cy//".daymean.1deg.nc"
          else if (trim(CRU%MetVersion) == "CRUJRA_2023") then
             fn = trim(fn)//"/fd/fd_v12_"//cy//".daymean.1deg.nc"
          end if
       end select
    endif

  end subroutine CRU_GET_FILENAME

  ! ------------------------------------------------------------------

  function cru_get_metyear(cru, cyear, runstartyear, ivar)
    ! Get the year of the meteorology input file for a given model year

    implicit none

    type(cru_type), intent(in) :: cru    ! CRU type
    integer,        intent(in) :: cyear  ! current model year
    integer,        intent(in) :: runstartyear  ! start year of whole model chain
    integer,        intent(in) :: ivar   ! years for special runs depend on variable
    integer                    :: cru_get_metyear

    if ( (trim(cru%run) == 'S0_TRENDY') &
         .or. (trim(cru%run) == 'S1_TRENDY') &
         .or. (trim(cru%run) == 'S0_TRENDY_CO2') &
         .or. (trim(cru%run) == 'S0_TRENDY_Ndep') ) then
       cru_get_metyear = 1901 + mod(cyear-runstartyear, cru%metrecyc)
    else if ( (trim(cru%run) == 'S0_TRENDY_Precip') &
         .or. (trim(cru%run) == 'S0_TRENDY_CO2_Precip') &
         .or. (trim(cru%run) == 'S0_TRENDY_CO2_Temp_Precip') &
         .or. (trim(cru%run) == 'S0_TRENDY_Temp_Precip') ) then
       if (ivar == 1) then
          cru_get_metyear = cyear
       else
          cru_get_metyear = 1901 + mod(cyear-runstartyear, cru%metrecyc)
       endif
       if  ( (trim(cru%run) == 'S0_TRENDY_CO2_Temp_Precip') &
            .or. (trim(cru%run) == 'S0_TRENDY_Temp_Precip') ) then
          if ((ivar == 6) .or. (ivar == 7)) then
             cru_get_metyear = cyear
          else
             cru_get_metyear = 1901 + mod(cyear-runstartyear, cru%metrecyc)
          endif
       endif
    else if ( (trim(cru%run) == 'S0_TRENDY_Temp') &
         .or. (trim(cru%run) == 'S0_TRENDY_CO2_Temp') ) then
       if ((ivar == 6) .or. (ivar == 7)) then
          cru_get_metyear = cyear
       else
          cru_get_metyear = 1901 + mod(cyear-runstartyear, cru%metrecyc)
       endif
    else if (trim(cru%run) == 'S2_TRENDY') then
       cru_get_metyear = cyear
    else if (trim(cru%run) == 'S2_TRENDY_precip0') then
       if (ivar == 1) then
          ! special for baseline precip
          cru_get_metyear = 1901 + mod(cyear-runstartyear, cru%metrecyc)
       else
          cru_get_metyear = cyear
       endif
    else if (trim(cru%run) == 'S2_TRENDY_precip') then
       if (ivar /= 1) then
          ! special for baseline non-precip
          cru_get_metyear = 1901 + mod(cyear-runstartyear, cru%metrecyc)
       else
          cru_get_metyear = cyear
       endif
    endif

  end function cru_get_metyear

  ! ------------------------------------------------------------------

  subroutine cru_read_metvals(cru, fid, vid, ii, t, land_x, land_y, filename)
    ! read the points land_x/land_y of input variable from input file

    implicit none

    type(cru_type),        intent(inout) :: cru  ! CRU type
    integer,               intent(in)    :: fid  ! file id
    integer,               intent(in)    :: vid  ! variable id
    integer,               intent(in)    :: ii   ! index of cru%met
    integer,               intent(in)    :: t    ! time steps
    integer, dimension(:), intent(in)    :: land_x  ! x indices
    integer, dimension(:), intent(in)    :: land_y  ! y indices
    character(len=*),      intent(in)    :: filename  ! file name

    integer :: k
    integer :: errstatus
    integer :: xds, yds
    real, dimension(:, :), allocatable :: tmparr

    if (cru%DirectRead) then
       ! read the current points directly into the vector for small domains
       do k=1, cru%mland
          errstatus = nf90_get_var(fid, vid, cru%met(ii)%metvals(k), &
               start=(/land_x(k), land_y(k), t/) )
          call handle_err(errstatus, "Reading directly from " // trim(filename))
       end do
    else
       ! or read the whole grid into tmparr and extract them from there
       xds = CRU%xdimsize
       yds = CRU%ydimsize
       allocate(tmparr(xds, yds))
       errstatus = nf90_get_var(fid, vid, tmparr, &
            start=(/1, 1, t/), count=(/xds, yds, 1/) )
       call handle_err(errstatus, "Reading from " // trim(filename))
       do k=1, cru%mland
          cru%met(ii)%metvals(k) = tmparr(land_x(k), land_y(k))
       end do
       deallocate(tmparr)
    endif

  end subroutine cru_read_metvals

  ! ------------------------------------------------------------------

  subroutine GET_CRU_CO2(CRU, CO2air)

    ! Get CO2 values for use with a CRU-NCEP run. Assign a static 1860 value if specified otherwise
    ! on the first call read all the annual values from a file into the CRU%CO2VALS array. On the first
    ! and subsequent

    implicit none

    type(CRU_TYPE), intent(inout) :: CRU    ! All the info needed for CRU met runs
    real,           intent(OUT)   :: CO2air ! A single annual value of CO2air in ppm for the current year.

    integer              :: iunit, iyear, IOS = 0
    character(len=200)   :: CO2FILE
    logical,        save :: CALL1 = .true.  ! A *local* variable recording the first call of this routine

    ! For S0_TRENDY, use only static 1860 CO2 value and return immediately
    if (trim(CRU%CO2) == "static1860") then
       !CO2air = 286.42   ! CO2 in ppm for 1860
       CO2air = 276.59   ! CO2 in ppm for 1700
       return
    else if (trim(CRU%CO2) == "static2011") then
       CO2air = 389.78   ! CO2 in ppm for 2011
       return
    else ! If not S0_TRENDY, varying CO2 values will be used...
       ! On the first call, allocate the CRU%CO2VALS array to store the entire history of annual CO2
       ! values, open the (ascii) CO2 file and read the values into the array.
       if (CALL1) then
          select case(trim(CRU%MetVersion))
          case("CRUJRA_2021", "VERIFY_2021")
             allocate(CRU%CO2VALS(1700:2020))
             CO2FILE = trim(CRU%BasePath)//"/co2/global_co2_ann_1700_2020.txt"
          case("CRUJRA_2022")
             allocate(CRU%CO2VALS(1700:2021))
             CO2FILE = trim(CRU%BasePath)//"/co2/global_co2_ann_1700_2021.txt"
          case("CRUJRA_2023")
             allocate(CRU%CO2VALS(1700:2022))
             CO2FILE = trim(CRU%BasePath)//"/co2/global_co2_ann_1700_2022.txt"
          end select

          call GET_UNIT(iunit)
          open(iunit, FILE=trim(CO2FILE), STATUS="OLD", ACTION="READ")
          do while (IOS == 0)
             read(iunit, FMT=*, IOSTAT=IOS) iyear, CRU%CO2VALS(iyear)
          end do
          close(iunit)

          CALL1 = .false.
       end if
       ! In all varying CO2 cases, return the element of the array for the current year
       ! as a single CO2 value.
       CO2air = CRU%CO2VALS(CRU%CYEAR)
    end if

  end subroutine GET_CRU_CO2

  ! ------------------------------------------------------------------

  subroutine GET_CRU_Ndep(CRU)

    ! Get Ndep values for use with a CRU-NCEP run. Assign a static 1860 value if specified otherwise
    ! on the first call read all the annual values from a file into the CRU%CO2VALS array. On the first
    ! and subsequent

    implicit none

    type(CRU_TYPE), intent(INOUT) :: CRU ! All the info needed for CRU met runs

    real, allocatable :: tmparr(:,:)
    integer :: k, t
    integer :: xds, yds ! Ndep file dimensions of long (x), lat (y)

    logical, save  :: CALL1 = .true. ! A *local* variable recording the first call of this routine
    character(400) :: NdepFILE

    ! Abbreviate dimensions for readability.
    xds = CRU%xdimsize
    yds = CRU%ydimsize
    allocate(tmparr(xds, yds))
    ! For S0_TRENDY, use only static 1860 CO2 value and return immediately

    ! On the first call, allocate the CRU%CO2VALS array to store the entire history of annual CO2
    ! values, open the (ascii) CO2 file and read the values into the array.
    if (CALL1) then
       if (trim(CRU%MetVersion) == "VERIFY_2021") then
          NdepFILE = trim(CRU%BasePath) // "/ndep/" // &
               "NOy_plus_NHx_dry_plus_wet_deposition_1850_2099_annual.0.125deg_Europe.nc"
       else
          NdepFILE = trim(CRU%BasePath) // "/ndep/" // &
               "NOy_plus_NHx_dry_plus_wet_deposition_1850_2099_annual.1deg.nc"
       end if

       ! Open the NDep and access the variables by their name and variable id.
       write(*,'(a)') 'Opening ndep data file: ' // trim(NdepFILE)
       write(logn,*)  'Opening ndep data file: ', trim(NdepFILE)

       ErrStatus = NF90_OPEN(trim(NdepFILE), NF90_NOWRITE, CRU%NdepF_ID)
       call HANDLE_ERR(ErrStatus, "Opening CRU file " // trim(NdepFILE))
       ErrStatus = NF90_INQ_VARID(CRU%NdepF_ID,'N_deposition', CRU%NdepV_ID)
       call HANDLE_ERR(ErrStatus, "Inquiring CRU var " // "N_deposition" // " in " // trim(NdepFILE))

       ! Set internal counter
       CRU%Ndep_CTSTEP = 1

       if ((trim(CRU%Ndep) == "static1860") &
            .or. (trim(CRU%Ndep) == "static2011") &
            .or. (CRU%CYEAR <= 1850)) then
          ! read Ndep at year 1850 (file starts at 1850)
          ! prior to TRENDYv10: year 1860
          CRU%Ndep_CTSTEP = 1
          t =  CRU%Ndep_CTSTEP
          ErrStatus = NF90_GET_VAR(CRU%NdepF_ID, CRU%NdepV_ID, tmparr, &
               start=(/1, 1, t/), count=(/xds, yds, 1/) )
          call HANDLE_ERR(ErrStatus, "Reading from " // trim(NdepFILE))
          do k = 1, CRU%mland
             CRU%NdepVALS(k) = tmparr(land_x(k), land_y(k))
          end do
       end if
       CALL1 = .false.
    end if

    if ((trim(CRU%Ndep) /= "static1860") &
         .and. (trim(CRU%Ndep) /= "static2011") &
         .and. (CRU%CYEAR > 1850)) then
       ! read Ndep at current year (noting that file starts at 1850 and ends in 2099)
       CRU%Ndep_CTSTEP = min(CRU%CYEAR, 2099) - 1850 + 1
       t =  CRU%Ndep_CTSTEP
       ErrStatus = NF90_GET_VAR(CRU%NdepF_ID, CRU%NdepV_ID, tmparr, &
            start=(/1, 1, t/), count=(/xds, yds, 1/))
       call HANDLE_ERR(ErrStatus, "Reading from " // trim(NdepFILE))
       do k = 1, CRU%mland
          CRU%NdepVALS(k) = tmparr(land_x(k), land_y(k))
       end do
    end if

  end subroutine GET_CRU_Ndep

  ! ------------------------------------------------------------------

  subroutine OPEN_CRU_MET(CRU)

    ! Opens each of the met files required for one year. This is where the distinction is made between
    ! the nominal run year (CYEAR) and the year of met required (MetYear), which is different for
    ! S0_TRENDY and S1_TRENDY than for a standard run (S2_TRENDY).

    implicit none

    type(CRU_TYPE), intent(INOUT) :: CRU ! All CRU-NCEP related quantities and flags

    integer       :: iVar            ! Loop counter through met variables
    integer       :: MetYear         ! Year of met to access. Equals CYEAR for normal runs, but
    ! must be calculated for S0_TRENDY and initialisation runs.
    ! integer, save :: RunStartYear    ! The value of CRU%CYEAR on the first call, also equals syear.
    integer :: RunStartYear    ! The value of CRU%CYEAR on the first call, also equals syear.
    ! ! Allows the calculation of MetYear during S0_TRENDY and init runs.
    ! logical, save :: CALL1 = .true.  ! A *local* variable recording the first call of this routine

    ! Keep the initial value of CYEAR for calculation of different MetYear if required.
    ! IF (CALL1) RunStartYear = 1710 ! edit vh !
    ! IF (CALL1) RunStartYear = 1691 ! edit vh !
    ! if (CALL1) RunStartYear = 1501 ! edit jk !
    RunStartYear = 1501 ! edit mc !

    do iVar = 1, CRU%NMET  ! For each met variable

       ! For S0_TRENDY and initialisation, calculate the required met year for repeatedly cycling through the
       ! 30 years of 1901-1930 spinup meteorology. For normal runs 1901-2015, MetYear = CYEAR.
       ! ! IF ( TRIM(CRU%Run) .EQ. 'S0_TRENDY' .OR.  ( TRIM(CRU%Run) .EQ. 'S1_TRENDY' )) THEN
       ! !   MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
       ! ! ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY' ) THEN
       ! !   MetYear = CRU%CYEAR
       ! ! ENDIF
       ! JK: according to v9 protocol, met should be recycled from 1901-1920 (i.e. 20 years)
       MetYear = cru_get_metyear(CRU, CRU%CYEAR, RunStartYear, iVar)
       call CRU_GET_FILENAME(CRU, MetYear, iVar, CRU%MetFile(iVar))

       ! Open the new met files and access the variables by their name and variable id.
       write(*,'(a)') 'Opening met data file: ' // trim(CRU%MetFile(iVar))
       write(logn,*)  'Opening met data file: ', CRU%MetFile(iVar)

       ErrStatus = NF90_OPEN(trim(CRU%MetFile(iVar)), NF90_NOWRITE, CRU%F_ID(iVar))
       call HANDLE_ERR(ErrStatus, "Opening CRU file " // trim(CRU%MetFile(iVar)))
       ErrStatus = NF90_INQ_VARID(CRU%F_ID(iVar),trim(CRU%VAR_NAME(iVar)), CRU%V_ID(iVar))
       call HANDLE_ERR(ErrStatus, "Inquiring CRU var " // trim(CRU%VAR_NAME(iVar)) // &
            " in " // trim(CRU%MetFile(iVar)))
    end do

    ! Set internal counter
    CRU%CTSTEP = 1

    ! CALL1 = .false. ! No longer the first call (saved).

  end subroutine OPEN_CRU_MET

  ! ------------------------------------------------------------------

  subroutine CRU_GET_DAILY_MET(CRU, LastDayOfYear)

    implicit none

    type(CRU_TYPE), intent(inout) :: CRU
    logical,        intent(IN)    :: LastDayOfYear

    integer :: iVar, ii
    integer :: t, tplus1              ! The current and next timestep
    integer :: fid, vid, tid          ! Netcdf id's for file, variable, and time
    integer :: xds, yds, tds          ! Metfile dimensions of long (x), lat (y), and time (t)
    integer :: MetYear                ! Year of meteorology currently in use
    integer :: NextMetYear            ! Next met year: Where to look for the nextTmin on Dec 31st
    integer :: PrevMetYear            ! Previous met year: Where to look for the prevTmax on Jan 1st
    character(LEN=200) :: filename

    ! integer, save     :: RunStartYear   ! The value of CRU%CYEAR on the first call, also equals syear.
    integer :: RunStartYear   ! The value of CRU%CYEAR on the first call, also equals syear.
    ! Allows the calculation of MetYear during S0_TRENDY and init runs.
    logical, save     :: CALL1 = .true. ! A *local* variable recording the first call of this routine
    real, allocatable :: tmparr(:,:)    ! packing into CRU%MET(iVar)%METVALS(k)
    logical :: existfile

    ! RunStartYear = CRU%CYEAR
    ! RunStartYear = 1691
    RunStartYear = 1501

    ! Abbreviate dimensions for readability.
    xds = CRU%xdimsize
    yds = CRU%ydimsize
    allocate(tmparr(xds, yds))

    ! Loop through all 9 met variables (not including prevTmax and nextTmin,
    ! which are addressed separately as special cases of Tmax and Tmin)

    do iVar=1, CRU%NMET
       PrevMetYear = cru_get_metyear(CRU, CRU%CYEAR-1, RunStartYear, iVar)
       MetYear     = cru_get_metyear(CRU, CRU%CYEAR,   RunStartYear, iVar)
       NextMetYear = cru_get_metyear(CRU, CRU%CYEAR+1, RunStartYear, iVar)
       ! NextMetYear = cru_get_metyear(CRU, CRU%CYEAR, RunStartYear, iVar)

       select case(iVar)
       case(Tmin)
          ! When Tmin comes up we will jump ahead to deal with nextTmin, since
          ! we are assigning Tmin from the previous value of nextTmin

          ! Bump the variable index ii to be the value for nextTmin and the time
          ! index ahead by 1 day.
          ii = nextTmin
          t  = CRU%CTSTEP
          tplus1 = CRU%CTSTEP + 1

          if (CALL1) then
             ! The first call reads the first timestep (t) into Tmin and
             ! the second timestep (t+1) into nextTmin.
             call cru_read_metvals(CRU, CRU%F_ID(iVar), CRU%V_ID(iVar), iVar, &
                  t, land_x(:), land_y(:), CRU%MetFile(iVar))
             call cru_read_metvals(CRU, CRU%F_ID(iVar), CRU%V_ID(iVar), ii, &
                  tplus1, land_x(:), land_y(:), CRU%MetFile(iVar))
          else if (LastDayOfYear) then
             ! old nextTmin is current Tmin
             cru%met(iVar)%metvals(:) = cru%met(ii)%metvals(:)
             ! If it is the last day of the year, open next year's metfile if
             ! possible, otherwise reuse current nextTmin
             t = 1 ! Time index is set to the first day of the next year
             call CRU_GET_FILENAME(CRU, NextMetYear, Tmin, filename)
             inquire(file=trim(filename), exist=existfile)
             if (existfile) then
                ! NextMetyear file exists
                ErrStatus = NF90_OPEN(trim(filename), NF90_NOWRITE, fid)
                call HANDLE_ERR(ErrStatus, "Opening CRU file " // trim(filename))
                ErrStatus = NF90_INQ_VARID(fid,trim(CRU%VAR_NAME(iVar)), vid)
                call HANDLE_ERR(ErrStatus, "Inquiring CRU var " // trim(filename))
                call cru_read_metvals(CRU, fid, vid, ii, &
                     t, land_x(:), land_y(:), filename)
                ! Close the next year's Tmin met file
                ErrStatus = NF90_CLOSE(fid)
                fid = -1
                call HANDLE_ERR(ErrStatus, "Closing CRU file "//trim(filename))
             ! else
                ! There is no more met file -> use old vale of nextTmin
                ! cru%met(ii)%metvals(:) = cru%met(ii)%metvals(:)
             endif
          else
             ! Common case:
             ! old nextTmin is current Tmin
             ! new nextTmin from met file
             cru%met(iVar)%metvals(:) = cru%met(ii)%metvals(:)
             call cru_read_metvals(CRU, CRU%F_ID(iVar), CRU%V_ID(iVar), ii, &
                  tplus1, land_x(:), land_y(:), CRU%MetFile(iVar))
          end if   ! End of if CALL1

       case(TMax)
          ! For Tmax, we need the previous day's Tmax for the Cesaraccio et al.
          ! subdiurnal temperature interpolation algorithm.
          ii = prevTmax
          t  = CRU%CTSTEP

          ! Common case:
          ! old Tmax is current prevTmax
          ! new Tmax comes from met file
          cru%met(ii)%metvals(:) = cru%met(iVar)%metvals(:)
          call cru_read_metvals(CRU, CRU%F_ID(iVar), CRU%V_ID(iVar), iVar, &
               t, land_x(:), land_y(:), CRU%MetFile(iVar))
          if (CALL1) then
             ! The first call reads the last timestep (tds) of previous year's
             ! met file into prevTmax if possible, otherwise use current Tmax
             call CRU_GET_FILENAME(CRU, PrevMetYear, iVar, filename)
             inquire(file=trim(filename), exist=existfile)
             if (existfile) then
                ! previous year met file exists
                ErrStatus = NF90_OPEN(trim(filename), NF90_NOWRITE, fid)
                call HANDLE_ERR(ErrStatus, "Opening CRU file "//trim(filename))
                ! obtain the size of the time dimension (tds), which locates the
                ! data for Dec 31st in the file
                ErrStatus = NF90_INQ_DIMID(fid, 'time', tid)
                ErrStatus = NF90_INQUIRE_DIMENSION(fid, tid, len=tds)
                call HANDLE_ERR(ErrStatus, "Inquiring 'time'" // trim(filename))
                ! obtain the variable ID of the Tmax variable
                ErrStatus = NF90_INQ_VARID(fid, trim(CRU%VAR_NAME(iVar)), vid)
                call HANDLE_ERR(ErrStatus, "Inquiring CRU var " // trim(filename))
                ! read Tmax
                call cru_read_metvals(CRU, fid, vid, ii, &
                     tds, land_x(:), land_y(:), filename)
                ! close previous met file
                ErrStatus = NF90_CLOSE(fid)
                fid = -1
                call HANDLE_ERR(ErrStatus, "Closing CRU file "//trim(filename))
             else
                ! There is no previous met file -> use current Tmax
                cru%met(ii)%metvals(:) = cru%met(iVar)%metvals(:)
             endif  ! End of If MetYear > CRU%Metstart
          endif

       case default  ! All variables other than Tmin and Tmax
          ii = iVar
          t  = CRU%CTSTEP
          ! Standard read of the current variable, for the current timestep
          call cru_read_metvals(CRU, CRU%F_ID(iVar), CRU%V_ID(iVar), ii, &
               t, land_x(:), land_y(:), CRU%MetFile(iVar))

       end select

    end do  ! do iVar=1, CRU%NMET -> end of loop through all met variables
    
    ! Convert pressure Pa -> hPa
    cru%met(pres)%metvals(:) = cru%met(pres)%metvals(:) / 100.

    ! Increment the internal timestep counter
    CRU%CTSTEP = CRU%CTSTEP + 1

    ! if (CALL1) CALL1 = .false.
    CALL1 = .false.

  end subroutine CRU_GET_DAILY_MET

  ! ------------------------------------------------------------------

  subroutine CRU_GET_SUBDIURNAL_MET(CRU, MET, CurYear, ktau, kend)

    ! Obtain one day of CRU-NCEP meteorology, subdiurnalise it using a weather
    ! and return the result to the CABLE driver.

    use cable_def_types_mod,   only: MET_TYPE, r_2
    use cable_IO_vars_module,  only: LANDPT, latitude
    use cable_common_module,   only: DOYSOD2YMDHMS
    use cable_weathergenerator,only: WEATHER_GENERATOR_TYPE, WGEN_INIT, &
         WGEN_DAILY_CONSTANTS, WGEN_SUBDIURNAL_MET
    use mo_utils,              only: eq

    implicit none

    type(CRU_TYPE), intent(inout) :: CRU
    integer,        intent(IN)    :: CurYear, ktau, kend

    ! Define MET the CABLE version, different from the MET defined and used
    ! within the CRU variable.
    ! The structure of MET_TYPE is defined in cable_define_types.F90
    type(MET_TYPE) :: MET

    ! Local variables
    logical   :: newday, LastDayOfYear  ! Flags for occurence of a new day (0 hrs) and the last day of the year.
    integer   :: iland                  ! Loop counter through 'land' cells (cells in the spatial domain)
    integer   :: itimestep              ! Loop counter through subdiurnal timesteps in a day
    integer   :: imetvar                ! Loop counter through met variables
    integer   :: dM, dD                 ! Met date as year, month, and day returned from DOYSOD2YMDHMS
    integer   :: is, ie                 ! Starting and ending vegetation type per land cell
    real      :: dt                     ! Timestep in seconds
    real      :: CO2air                 ! CO2 concentration in ppm
    type(WEATHER_GENERATOR_TYPE), save :: WG
    logical,                      save :: CALL1 = .true.  ! A *local* variable recording the first call of this routine

    ! Purely for readability...
    dt = CRU%DTsecs

    ! On first step read and check CRU settings and read land-mask
    if (CALL1) then
       if (CRU%ReadDiffFrac) then
          cable_user%calc_fdiff = .false.
       else
          cable_user%calc_fdiff = .true.
       endif
       
       call WGEN_INIT(WG, CRU%mland, latitude, dt)
    endif
       
    ! Pass time-step information to CRU
    CRU%CYEAR = CurYear
    CRU%ktau  = ktau     ! ktau is the current timestep in the year.

    ! this only works with CANBERRA cable_driver, as ktau ! restarts on Jan 1 !
    ! Based on the ktau timestep, calculate date and time information (the same
    ! for the whole spatial dimension.)
    met%hod (:) = real(mod((ktau-1) * nint(dt), int(SecDay))) / 3600.  ! Hour of the day
    met%doy (:) = int(real(ktau-1) * dt / SecDay ) + 1                 ! Day of Year = days since 0 hr 1st Jan
    met%year(:) = CurYear                                              ! Current year

    ! Using the day-of-year and seconds-of-day calculate the month and
    ! day-of-month, using the time information for the first land cell only
    ! (because they will be the same across the domain).
    !                            In      In         In        Out       Out       Optional Out
    ! SUBROUTINE DOYSOD2YMDHMS( Year, Yearday, SecondsOfDay, Month, DayOfMonth, [Hour, Min, Sec])
    call DOYSOD2YMDHMS(CurYear, int(met%doy(1)), int(met%hod(1)) * 3600, dM, dD)
    met%moy (:) = dM     ! Record the month

    ! It's a new day if the hour of the day is zero.
    newday = eq(met%hod(landpt(1)%cstart), 0.0)

    ! Beginning-of-year accounting
    if (ktau == 1) then  ! ktau is always reset to 1 at the start of the year.
       ! Read a new annual CO2 value and convert it from ppm to mol/mol
       call GET_CRU_CO2(CRU, CO2air)
       met%ca(:) = CO2air / 1.0e6

       call GET_CRU_Ndep(CRU)
       do iland = 1, CRU%mland
          met%Ndep(landpt(iland)%cstart:landpt(iland)%cend) = &
               CRU%NdepVALS(iland)*86400000.  ! kg/m2/s > g/m2/d (1000.*3600.*24.)
       end do
       ! Open a new annual CRU-NCEP met file.
       call OPEN_CRU_MET(CRU)
    endif

    ! %%%%%% PRB to add his own comments from here down for this routine.
    ! Now get the Met data for this day
    if (newday) then
       !print *, CRU%CTSTEP, ktau, kend
       !   CALL CPU_TIME(etime)
       !   PRINT *, 'b4 daily ', etime, ' seconds needed '
       LastDayOfYear = ktau == (kend-(nint(SecDay/dt)-1))

       call CRU_GET_DAILY_MET(CRU, LastDayOfYear)
       !   CALL CPU_TIME(etime)
       !   PRINT *, 'after daily ', etime, ' seconds needed '
       !   STOP

       ! Air pressure assumed to be constant over day
       do iland = 1, CRU%mland
          met%pmb(landpt(iland)%cstart:landpt(iland)%cend) = CRU%MET(pres)%METVALS(iland) !CLN interpolation??
       end do

       ! Convert wind from u and v components to wind speed by Pythagorean Theorem
       WG%WindDay        = real(sqrt((CRU%MET(uwind)%METVALS * CRU%MET(uwind)%METVALS) +  &
            (CRU%MET(vwind)%METVALS * CRU%MET(vwind)%METVALS)), r_2)
       ! Convert all temperatures from K to C
       WG%TempMinDay     = real(CRU%MET(  Tmin  )%METVALS - 273.15, r_2)
       WG%TempMaxDay     = real(CRU%MET(  Tmax  )%METVALS - 273.15, r_2)
       WG%TempMinDayNext = real(CRU%MET(NextTmin)%METVALS - 273.15, r_2)
       WG%TempMaxDayPrev = real(CRU%MET(PrevTmax)%METVALS - 273.15, r_2)
       ! Convert solar radiation from J /m2/s to MJ/m2/d
       WG%SolarMJDay     = real(CRU%MET(swdn)%METVALS * 1.e-6 * SecDay, r_2) ! ->[MJ/m2/d]
       ! Convert precip from mm to m/day
       WG%PrecipDay      = real(max(CRU%MET(rain)%METVALS  / 1000., 0.0), r_2) ! ->[m/d]
       !WG%PrecipDay      = max(CRU%MET(rain)%METVALS  / 1000., 0.0)/2.0 ! ->[m/d] ! test vh !
       WG%SnowDay        = 0.0_r_2
       WG%VapPmbDay = real(esatf(real(WG%TempMinDay, sp)), r_2)
       call WGEN_DAILY_CONSTANTS(WG, CRU%mland, int(met%doy(1))+1)

       ! To get the diurnal cycle for lwdn get whole day and scale with
       ! LWDN from file later

       CRU%AVG_LWDN(:) = 0.
       do itimestep = 1, nint(SecDay/dt)
          call WGEN_SUBDIURNAL_MET(WG, CRU%mland, itimestep-1)
          CRU%AVG_LWDN = CRU%AVG_LWDN + real(WG%PhiLD)
       end do
       CRU%AVG_LWDN = CRU%AVG_LWDN / (SecDay/dt)
    end if ! End of If newday

    ! Decision has been made, that first tstep of the day is at 0:01 am
    call WGEN_SUBDIURNAL_MET(WG, CRU%mland, nint(met%hod(1)*3600./dt))

    ! assign to cable variables

    ! strangely, met% is not save over Wait all in MPI...!
    !  met%pmb = CRU%MET(pres)%METVALS !CLN interpolation??

    ! Assign weather-generated data, or daily values, as required to CABLE
    ! variables.
    do iland = 1, CRU%mland
       is = landpt(iland)%cstart
       ie = landpt(iland)%cend

       met%precip(is:ie) = real(WG%Precip(iland)) ! test vh !

       ! Cable's swdown is split into two components, visible and nir, which
       ! get half of the CRU-NCEP swdown each.
       met%fsd(is:ie,1) = real(WG%PhiSD(iland) * 0.5_r_2)  ! Visible
       met%fsd(is:ie,2) = real(WG%PhiSD(iland) * 0.5_r_2)  ! NIR

       if (CRU%ReadDiffFrac) then ! read from met forcing, otherwise calculated in cable_radiation.F90
          met%fdiff(is:ie) = CRU%MET(fdiff)%METVALS(iland)
       endif

       ! Convert C to K for cable's tk
       met%tk(is:ie)     = real(WG%Temp(iland) + 273.15_r_2)

       met%ua(is:ie)     = real(WG%Wind(iland))
       met%coszen(is:ie) = real(WG%coszen(iland))

       ! For longwave down, scale the diurnal series returned by the weather
       ! generator (WG%PhiLD(iland)) to the daily value from CRU-NCEP.
       met%fld(is:ie) = CRU%MET(lwdn)%METVALS(iland) * real(WG%PhiLD(iland)) / CRU%AVG_LWDN(iland)

       ! Specific humidity (qair g/g) was not sent to the weather generator.
       ! Here we assign the daily value to the whole diurnal cycle
       met%qv(is:ie) = CRU%MET(qair)%METVALS(iland)

       ! rel humidity (%)
       met%rhum(is:ie)  = real(WG%VapPmb(iland))/esatf(real(WG%Temp(iland),sp)) * 100.0
       met%u10(is:ie)   = met%ua(is:ie)
       ! initialise within canopy air temp
       met%tvair(is:ie) = met%tk(is:ie)
       met%tvrad(is:ie) = met%tk(is:ie)
       met%pdep = 0.0

       ! calculate snowfall based on total precip and air T
       ! (ref Jin et al. Table II, Hyd Proc, 1999)
       !    if (WG%Temp(iland) > 2.5) then
       !       met%precip_sn(is:ie) = 0.0
       !    elseif ((WG%Temp(iland) <= 2.5) .and. (WG%Temp(iland) > 2.0)) then
       !       met%precip_sn(is:ie) = 0.6* met%precip(is:ie)
       !    elseif ((WG%Temp(iland) <= 2.0) .and. (WG%Temp(iland) > 0.0)) then
       !       ! this factor can be > 1 !!!
       !       met%precip_sn(is:ie) = (1.0 - (54.62 - 0.2 *(WG%Temp(iland) + 273.15)))* met%precip(is:ie)
       !    elseif (WG%Temp(iland) <= 0.0) then
       !       met%precip_sn(is:ie) = met%precip(is:ie)
       !    endif

       if (WG%Temp(iland) <= 0.0_r_2) then
          met%precip_sn(is:ie) = met%precip(is:ie)
       else
          met%precip_sn(is:ie) = 0.0
       endif
       ! met%precip(is:ie) = met%precip(is:ie) - met%precip_sn(is:ie)
    end do

    ! initialise within canopy air temp
    met%tvair = met%tk
    met%tvrad = met%tk

    ! If this is the end of the year or the end of the met, close the current met files.
    if (ktau == kend) then
       do imetvar=1, CRU%NMET
          ErrStatus = NF90_CLOSE(CRU%F_ID(imetvar))
          CRU%F_ID(imetvar) = -1
          call HANDLE_ERR(ErrStatus, "Closing CRU file"//trim(CRU%MetFile(imetvar)))
       end do
    end if

    ! CALL1 is over...
    CALL1 = .false.

  contains

    elemental function Esatf(TC)
      ! ------------------------------------------------------------------------------
      ! At temperature TC [deg C], return saturation water vapour pressure Esatf [mb]
      ! from Teten formula.
      ! MRR, xx/1987
      ! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
      !                 just like intrinsic functions.
      ! MRR, 12-mar-02: Convert Qsatf (specific humidity routine) to Esatf
      ! ------------------------------------------------------------------------------
      implicit none

      real(sp), intent(in) :: TC          ! temp [deg C]
      real(sp)             :: Esatf       ! saturation vapour pressure [mb]

      real(sp) :: TCtmp                   ! local
      real(sp),parameter:: A = 6.106      ! Teten coefficients
      real(sp),parameter:: B = 17.27      ! Teten coefficients
      real(sp),parameter:: C = 237.3      ! Teten coefficients

      TCtmp = TC                            ! preserve TC
      if (TCtmp > 100.0) TCtmp = 100.0   ! constrain TC to (-40.0, 100.0)
      if (TCtmp < -40.0) TCtmp = -40.0

      Esatf = A*exp(B*TCtmp/(C+TCtmp))    ! sat vapour pressure (mb)

    end function Esatf

  end subroutine CRU_GET_SUBDIURNAL_MET

  ! ------------------------------------------------------------------

  subroutine cru_close(CRU)

    implicit none

    type(cru_type), intent(inout) :: CRU

    integer :: i

    write(*,*) 'Closing CRU files.'
    do i=1, CRU%nmet
       if (CRU%f_id(i) > -1) then
          errstatus = nf90_close(CRU%f_id(i))
          call handle_err(errstatus, "Closing CRU met file " // trim(CRU%MetFile(i)))
       end if
    end do

    if (CRU%Ndepf_id > -1) then
       errstatus = nf90_close(CRU%Ndepf_id)
       call handle_err(errstatus, "Closing CRU Ndep file.")
    end if

  end subroutine cru_close

end module CABLE_CRU
