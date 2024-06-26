MODULE CABLE_CRU

!*## Purpose
!
! The master module to handle the meteorological forcing when using the CRU
! dataset.

USE NetCDF,                 ONLY: NF90_OPEN, NF90_NOWRITE, NF90_INQ_DIMID,&
                                  NF90_INQUIRE_DIMENSION, NF90_INQ_VARID,&
                                  NF90_GET_VAR, NF90_CLOSE, NF90_NOERR
USE cable_abort_module,     ONLY: cable_abort, nc_abort
USE cable_common_module,    ONLY: handle_err, get_unit, is_leapyear,&
                                cable_user, doysod2ymdhms
USE cable_io_vars_module,   ONLY: logn, land_x, land_y, exists, nMetPatches,&
                                latitude, longitude, lat_all, lon_all, landPt,&
                                xDimSize, yDimSize, sdoy, smoy, syear, shod,&
                                mask, metGrid
USE cable_input_module,     ONLY: SPATIO_TEMPORAL_DATASET, find_variable_ID,&
                                prepare_spatiotemporal_dataset, read_metvals,&
                                open_at_first_file
USE cable_def_types_mod,    ONLY: MET_TYPE, r_2, mLand
USE cable_weathergenerator, ONLY: WEATHER_GENERATOR_TYPE, WGEN_INIT,&
                                WGEN_DAILY_CONSTANTS, WGEN_SUBDIURNAL_MET
USE mo_utils,               ONLY: eq
#IFDEF __MPI__
USE mpi,                    ONLY: MPI_Abort
#ENDIF

IMPLICIT NONE

TYPE CRU_MET_TYPE
  REAL, DIMENSION(:), ALLOCATABLE :: MetVals
END TYPE CRU_MET_TYPE

TYPE CRU_TYPE
  INTEGER   :: mLand              ! Number of land cells in the run
  INTEGER   :: nMet               ! Number of Meteorological variables
  INTEGER   :: xDimSize, yDimSize ! Size of longitude and latitude dimensions
  INTEGER   :: cYear              ! CABLE main year
  INTEGER   :: CTStep             ! Day of year I think?
  INTEGER   :: DtSecs             ! Size of the timestep in seconds
  INTEGER   :: MetRecyc           ! Period of the met forcing recycling
  INTEGER   :: Ktau

  INTEGER, DIMENSION(10)    :: FileID, VarId    ! File and variable IDs for
  INTEGER                   :: NDepFId, NDepVId ! reading netCDF

  REAL, DIMENSION(:), ALLOCATABLE ::  avg_lwdn  ! Average longwave down rad,
                                                ! for Swinbank method

  REAL, DIMENSION(:), ALLOCATABLE ::  CO2Vals   ! Atmospheric carbon values
  REAL, DIMENSION(:), ALLOCATABLE ::  NDepVals  ! Nitrogen deposition values

  LOGICAL   :: DirectRead         ! Whether to read in entire met array
  LOGICAL   :: ReadDiffFrac       ! Read diff fraction or calculate it
  LOGICAL   :: LeapYears

  LOGICAL, DIMENSION(10)    :: isRecycled
  LOGICAL, DIMENSION(:,:), ALLOCATABLE  :: LandMask ! The logical landmask

  CHARACTER(LEN=16)   :: CO2Method   ! Method for choosing atmospheric CO2
  CHARACTER(LEN=16)   :: NDepMethod     ! Method for choosing N deposition

  TYPE(CRU_MET_TYPE), DIMENSION(12)   :: Met

  ! The spatio-temporal datasets we use to generate the forcing data.
  TYPE(SPATIO_TEMPORAL_DATASET), DIMENSION(10)  :: MetDatasets
  TYPE(SPATIO_TEMPORAL_DATASET)                 :: NDepDataset
END TYPE CRU_TYPE

! Set some private parameters
REAL, PRIVATE, PARAMETER  :: SecDay = 86400.
INTEGER, PRIVATE, PARAMETER  :: rain = 1, lwdn = 2, swdn = 3, pres = 4, qair = 5,&
                                tmax = 6, tmin = 7, uwind = 8, vwind = 9, fdiff = 10,&
                                prevTmax = 11, nextTmin = 12
INTEGER, PRIVATE, PARAMETER  :: sp = kind(1.0)
INTEGER, PRIVATE             :: ErrStatus

CONTAINS

SUBROUTINE CRU_INIT(CRU)
  !*## Purpose
  !
  ! Initialise the program to handle met forcing from the CRU dataset.
  !
  !## Method
  !
  ! Reads the cru.nml namelist and configures the passed ```CRU``` derived type
  ! with information about the way the met forcing is handled, the location of
  ! the datasets on the filesystem, and reads in the landmask to globally shared
  ! variables.

  ! The CRU derived type which is used to control the met forcing.
  TYPE(CRU_TYPE), INTENT(OUT)    :: CRU

  ! We want one filename for each variable, stored in a predefined index.
  CHARACTER(LEN=256), dimension(13) :: InputFiles

  ! Landmask to spatially filter the data
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: LandMask

  ! Iterator variable for the variables
  INTEGER   :: VarIndx

  ! Checker for the NetCDF io
  INTEGER   :: ok
  LOGICAL :: IsOpen

  ! Start with the things we want from the namelist. The namelist must set
  ! the filenames to read from, the method of choosing atmospheric carbon and
  ! nitrogen deposition, and the timestep.
  CALL read_MET_namelist_cbl(InputFiles, CRU)

  ! Read the landmask and allocate appropriate memory for the array variables
  CALL read_landmask(InputFiles(13), CRU)

  ! Build the spatio-temporal datasets for each necessary datatype
  BuildKeys: DO VarIndx = 1, CRU%nMet
    CALL prepare_spatiotemporal_dataset(InputFiles(VarIndx),&
      CRU%MetDatasets(VarIndx))
  END DO BuildKeys

  ! Set the possible variable names for the main Met variables
  CALL read_variable_names(CRU%MetDatasets)

  ! Open the datasets at the first file so we don't need CALL1 behaviour later
  InitialiseDatasets: DO VarIndx = 1, CRU%nMet
    CALL open_at_first_file(CRU%MetDatasets(VarIndx))
  END DO InitialiseDatasets

  ! Set up the carbon reader
  CALL prepare_temporal_dataset(InputFiles(11), CRU%CO2Vals)

  ! Set up the nitrogen deposition reader
  CALL prepare_spatiotemporal_dataset(InputFiles(12), CRU%NDepDataset)
  ! For now, set the file index to 1 since we know its only one file
  CRU%NDepDataset%CurrentFileIndx = 1

  ok = NF90_OPEN(CRU%NDepDataset%Filenames(1), NF90_NOWRITE,&
    CRU%NDepFID)
  CALL handle_err(ok, "Opening NDep file")

  ok = NF90_INQ_VARID(CRU%NDepFID, "N_deposition",&
    CRU%NDepVID)
  CALL handle_err(ok, "Finding NDep variable")
END SUBROUTINE CRU_INIT


SUBROUTINE CRU_GET_SUBDIURNAL_MET(CRU, MET, CurYear, ktau, kend)

  ! Obtain one day of CRU-NCEP meteorology, subdiurnalise it using a weather
  ! generator and return the result to the CABLE driver.

  IMPLICIT NONE

  TYPE(CRU_TYPE), INTENT(inout) :: CRU
  INTEGER,        INTENT(IN)    :: CurYear, ktau, kend

  ! Define MET the CABLE version, different from the MET defined and used
  ! within the CRU variable.
  ! The structure of MET_TYPE is defined in cable_define_types.F90
  type(MET_TYPE) :: MET

  ! Local variables
  logical   :: newday, LastDayOfYear  ! Flags for occurence of a new day (0 hrs) and the last day of the year.
  INTEGER   :: iland                  ! Loop counter through 'land' cells (cells in the spatial domain)
  INTEGER   :: itimestep              ! Loop counter through subdiurnal timesteps in a day
  INTEGER   :: imetvar                ! Loop counter through met variables
  INTEGER   :: dM, dD                 ! Met date as year, month, and day returned from DOYSOD2YMDHMS
  INTEGER   :: is, ie                 ! Starting and ending vegetation type per land cell
  REAL      :: dt                     ! Timestep in seconds
 ! Store the CO2Air as an array
  REAL, DIMENSION(:), ALLOCATABLE    :: CO2air                 ! CO2 concentration in ppm
  type(WEATHER_GENERATOR_TYPE), save :: WG
  logical,                      save :: CALL1 = .true.  ! A *local* variable recording the first call of this routine
  INTEGER :: VarIter

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

  ! Allocate memory for the CO2 values
  ALLOCATE(CO2Air(SIZE(met%ca)))

  ! Beginning-of-year accounting
  if (ktau == 1) then  ! ktau is always reset to 1 at the start of the year.
     CRU%CTStep = 1
     ! Read a new annual CO2 value and convert it from ppm to mol/mol
     call GET_CRU_CO2(CRU, CO2air)
     met%ca(:) = CO2air(:) / 1.0e6

     call GET_CRU_Ndep(CRU)
     do iland = 1, CRU%mland
        met%Ndep(landpt(iland)%cstart:landpt(iland)%cend) = &
             CRU%NdepVALS(iland)*86400000.  ! kg/m2/s > g/m2/d (1000.*3600.*24.)
     end do
     ! Open a new annual CRU-NCEP met file.
     !call OPEN_CRU_MET(CRU)
  endif

  ! %%%%%% PRB to add his own comments from here down for this routine.
  ! Now get the Met data for this day
  if (newday) then
     !print *, CRU%CTSTEP, ktau, kend
     !   CALL CPU_TIME(etime)
     !   PRINT *, 'b4 daily ', etime, ' seconds needed '
     LastDayOfYear = ktau == (kend-(nint(SecDay/dt)-1))

     call CRU_GET_DAILY_MET(CRU, LastDayOfYear)
     ! Scale presuure to hPa
     CRU%Met(pres)%MetVals(:) = CRU%Met(pres)%MetVals(:) / 100.

     CRU%CTStep = CRU%CTStep + 1
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
  !if (ktau == kend) then
  !   do imetvar=1, CRU%NMET
  !      ErrStatus = NF90_CLOSE(CRU%F_ID(imetvar))
  !      CRU%F_ID(imetvar) = -1
  !      call HANDLE_ERR(ErrStatus, "Closing CRU file"//trim(CRU%MetFile(imetvar)))
  !   end do
  !end if

  ! CALL1 is over...
  CALL1 = .false.

end subroutine CRU_GET_SUBDIURNAL_MET

SUBROUTINE cru_get_daily_met(CRU, LastDayOfYear)
  TYPE(CRU_TYPE), INTENT(INOUT)   :: CRU
  LOGICAL, INTENT(IN)             :: LastDayOfYear

  ! The year of met forcing we use depends on our choice of configuration.
  ! Sometimes, recycle through a subset of data, and others we use the sim year.
  ! So we want to store both.
  INTEGER   :: RecycledYear, MetYear

  ! We want to handle the nextTmin and prevTmax things separately so we don't
  ! mess with the data stored in CRU, initialise variables to handle this
  INTEGER   :: DummyYear, DummyDay, DaysInYear

  ! Define iteration variable
  INTEGER   :: VarIndx

  ! Determine the recycled and sim year
  RecycledYear = 1901 + MOD(CRU%cYear - 1501, CRU%metRecyc)

  ! Iterate through the base variables
  IterateVariables: DO VarIndx = 1, CRU%nMet
    MetYear = CRU%cYear
    ! If the variable is time recycled
    IF (CRU%isRecycled(VarIndx)) THEN
      MetYear = RecycledYear
    END IF

    ! CRU%CTStep is the current day, which we use to index the netCDF arrays
    CALL read_metvals(CRU%MetDatasets(VarIndx), CRU%Met(VarIndx)%MetVals,&
      land_x, land_y, MetYear, CRU%CTStep, CRU%LeapYears, CRU%xDimSize,&
      CRU%yDimSize, CRU%DirectRead)
  END DO IterateVariables

  ! Now the variables with special handling- nextTmin and prevTmax
  ! The easiest way to do this is to simply change the day by 1, check if we
  ! need to change the year, and go through the same read process. This is a
  ! little inefficient as we'll be messing with the Tmin and Tmax dataset
  ! structs and in select instances (at the end of the era of a particular file)
  ! opening and closing the io stream unnecessarily, but I think it's a minor
  ! evil.

  ! Address prevTmax first
  ! Assume it's a regular day of the year
  DummyYear = CRU%cYear
  DummyDay = CRU%CTStep - 1

  ! Special handling at the first day of the year
  IF (CRU%CTSTEP == 1) THEN
    ! Go back to previous year
    DummyYear = CRU%cYear - 1
    DummyDay = 365
    IF ((CRU%LeapYears) .AND. (is_leapyear(DummyYear))) THEN
      DummyDay = 366
    END IF
  END IF

  ! Was the Tmax recycled?
  IF (CRU%isRecycled(Tmax)) THEN
    DummyYear = 1901 + MOD(DummyYear - 1501, CRU%metRecyc)
  END IF

  ! Now we just need to call cru_read_metvals with the Tmax Dataset reader and
  ! the prevTmax array to write to
  CALL read_metvals(CRU%MetDatasets(Tmax), CRU%Met(prevTmax)%MetVals, land_x,&
    land_y, DummyYear, DummyDay, CRU%LeapYears, CRU%xDimSize, CRU%yDimSize,&
    CRU%DirectRead)

  ! Now do nextTmin
  ! Assume it's a regular day of the year
  DummyDay = CRU%CTStep + 1
  DummyYear = CRU%cYear

  ! Special handling at the last day of the year
  DaysInYear = 365
  IF ((CRU%LeapYears) .AND. (is_leapyear(DummyYear))) THEN
    DaysInYear = 366
  END IF

  IF (CRU%CTStep == DaysInYear) THEN
    ! We're at the end of a year
    DummyDay = 1
    DummyYear = CRU%cYear + 1
  END IF

  ! Was the Tmin recycled?
  IF (CRU%isRecycled(Tmin)) THEN
    DummyYear = 1901 + MOD(DummyYear - 1501, CRU%metRecyc)
  END IF

  CALL read_metvals(CRU%MetDatasets(Tmin), CRU%Met(nextTmin)%MetVals, land_x,&
    land_y, DummyYear, DummyDay, CRU%LeapYears, CRU%xDimSize, CRU%yDimSize,&
    CRU%DirectRead)

END SUBROUTINE cru_get_daily_met

SUBROUTINE read_MET_namelist_cbl(InputFiles, CRU)
  !*## Purpose
  !
  ! Set the metadata for the met forcing.
  !
  !## Method
  !
  ! Reads the namelist at ```cru.nml```, extracts a series of filename
  ! templates used to build the ```SPATIO_TEMPORAL_DATASET```(s) used to manage
  ! the input met data and sets metadata about the met forcing in the ```CRU```
  ! derived type.

  ! Internally, we're going to pick up the individually separated inputs and
  ! pack them into a more convenient data structure. We don't expect the user
  ! to know the order in which we store the MET inputs internally, so we need
  ! to access them by name.
  CHARACTER(LEN=256), DIMENSION(13), INTENT(OUT) :: InputFiles

  TYPE(CRU_TYPE), INTENT(OUT)  :: CRU    ! The master CRU structure

  ! Now we assign variables to temporarily store the files in the namelist
  ! In an ideal world, we'd have the users write directory to the InputFiles
  ! string array, but we can't expect everyone to know the order the files are
  ! expected in. So we first read them from a recognisable name, then pass them
  ! to the array.
  CHARACTER(LEN=256)  :: rainFile, lwdnFile, swdnFile, presFile, qairFile,&
                         TmaxFile, TminFile, uwindFile, vwindFile, fdiffFile,&
                         CO2File, NDepFile, landmaskFile
  LOGICAL             :: rainRecycle, lwdnRecycle, swdnRecycle, presRecycle,&
                         qairRecycle, TmaxRecycle, TminRecycle, uwindRecycle,&
                         vwindRecycle, fdiffRecycle
  CHARACTER(LEN=16)   :: CO2Method, NDepMethod
  INTEGER             :: MetRecyc
  REAL                :: DtHrs
  LOGICAL             :: ReadDiffFrac, LeapYears, DirectRead

  ! Need a unit to handle the io
  INTEGER             :: nmlUnit

  ! Set the initial values for the filenames, as there are many instances where
  ! not all are required. We will set them all initially to "None", which we
  ! can easily check against later.
  rainFile = "None"
  lwdnFile = "None"
  swdnFile = "None"
  presFile = "None"
  qairFile = "None"
  TmaxFile = "None"
  TminFile = "None"
  uwindFile = "None"
  vwindFile = "None"
  fdiffFile = "None"
  CO2File = "None"
  NDepFile = "None"
  landmaskFile = "None"

  ! For the recycling booleans
  rainRecycle = .FALSE.
  lwdnRecycle = .FALSE.
  swdnRecycle = .FALSE.
  presRecycle = .FALSE.
  qairRecycle = .FALSE.
  TmaxRecycle = .FALSE.
  TminRecycle = .FALSE.
  uwindRecycle = .FALSE.
  vwindRecycle = .FALSE.
  fdiffRecycle = .FALSE.

  ! Defaults for the other inputs
  CO2Method = "Yearly"
  NDepMethod = "Yearly"
  MetRecyc = 20
  ReadDiffFrac = .TRUE.
  DirectRead = .FALSE.

  ! Set up and read the namelist
  NAMELIST /crunml/ rainFile, lwdnFile, swdnFile, presFile, qairFile,&
                    TmaxFile, TminFile, uwindFile, vwindFile, fdiffFile,&
                    CO2File, NDepFile, landmaskFile,&
                    rainRecycle, lwdnRecycle, swdnRecycle, presRecycle,&
                    qairRecycle, TmaxRecycle, TminRecycle, uwindRecycle,&
                    vwindRecycle, fdiffRecycle,&
                    ReadDiffFrac, CO2Method, NDepMethod, MetRecyc, LeapYears,&
                    DtHrs, DirectRead

  ! Get a temporary unique ID and use it to read the namelist
  CALL get_unit(nmlUnit)
  OPEN(nmlUnit, file = "cru.nml", status = 'old', action = 'read')
  READ(nmlUnit, nml = crunml)
  CLOSE(nmlUnit)

  ! Now pack the filepaths into the data structure we want to transport around
  InputFiles(rain) = rainFile
  InputFiles(lwdn) = lwdnFile
  InputFiles(swdn) = swdnFile
  InputFiles(pres) = presFile
  InputFiles(qair) = qairFile
  InputFiles(Tmax) = TmaxFile
  InputFiles(Tmin) = TminFile
  InputFiles(uwind) = uwindFile
  InputFiles(vwind) = vwindFile
  InputFiles(fdiff) = fdiffFile
  InputFiles(11) = CO2File
  InputFiles(12) = NDepFile
  InputFiles(13) = landmaskFile

  ! Set the recycling booleans in the CRU struct
  CRU%IsRecycled(rain) = rainRecycle
  CRU%IsRecycled(lwdn) = lwdnRecycle
  CRU%IsRecycled(swdn) = swdnRecycle
  CRU%IsRecycled(pres) = presRecycle
  CRU%IsRecycled(qair) = qairRecycle
  CRU%IsRecycled(Tmax) = TmaxRecycle
  CRU%IsRecycled(Tmin) = TminRecycle
  CRU%IsRecycled(uwind) = uwindRecycle
  CRU%IsRecycled(vwind) = vwindRecycle
  CRU%IsRecycled(fdiff) = fdiffRecycle

  ! Bind the remaining config variables to the CRU structure
  CRU%CO2Method = CO2Method
  CRU%NDepMethod = NDepMethod
  ! Convert the hourly timestep to seconds
  CRU%DtSecs = int(DtHrs * 3600.)
  CRU%MetRecyc = MetRecyc
  CRU%ReadDiffFrac = ReadDiffFrac
  CRU%LeapYears = LeapYears
  CRU%DirectRead = DirectRead

  ! Set the number of met variables, based on ReadDiffFrac
  IF (CRU%ReadDiffFrac) THEN
    CRU%nMet = 10
  ELSE
    CRU%nMet = 9
  END IF

END SUBROUTINE read_MET_namelist_cbl

!------------------------------------------------------------------------------!

SUBROUTINE read_variable_names(STDatasets)
  !*## Purpose
  !
  ! Prepare the possible NetCDF names for the met forcing variables.
  !
  !## Method
  !
  ! Read the ```met_names.nml``` namelist to build a series of arrays
  ! containing the possible NetCDF names for each variable. The namelist
  ! contains, for each variable, the number of possible names to ```ALLOCATE```
  ! string arrays, and then the string arrays themselves.

  TYPE(SPATIO_TEMPORAL_DATASET), DIMENSION(10)  :: STDatasets

  ! We need a unit to reference the namelist file, arrays to store the set of
  ! names for each variable, and integers to store how many names there are for
  ! each variable. We set the size of the name arrays to 10 for now, because it
  ! seems sufficient. If it ever ends up requiring more, it's probably an
  ! indication we should switch to option 1.
  INTEGER   :: nmlUnit
  CHARACTER(LEN=32), DIMENSION(10)    :: rainNames, lwdnNames, swdnNames,&
                                         presNames, qairNames, TmaxNames,&
                                         TminNames, uWindNames, vWindNames,&
                                         fdiffNames

  INTEGER   :: rainNo, lwdnNo, swdnNo, presNo, qairNo, TmaxNo, TminNo, uwindNo,&
               vwindNo, fdiffNo

  ! Set up and read the namelist
  NAMELIST /MetNames/ rainNames, lwdnNames, swdnNames, presNames, qairNames,&
                    TmaxNames, TminNames, uwindNames, vwindNames, fdiffNames,&
                    rainNo, lwdnNo, swdnNo, presNo, qairNo, TmaxNo, TminNo,&
                    uwindNo, vwindNo, fdiffNo

  CALL get_unit(nmlUnit)
  OPEN(nmlUnit, FILE = "met_names.nml", STATUS = 'old', ACTION = 'read')
  READ(nmlUnit, NML = MetNames)
  CLOSE(nmlUnit)

  ! Allocate memory for the names in the SpatioTemporalDataset struct
  ALLOCATE(STDatasets(1)%VarNames(rainNo), STDatasets(2)%VarNames(lwdnNo),&
    STDatasets(3)%VarNames(swdnNo), STDatasets(4)%VarNames(presNo),&
    STDatasets(5)%VarNames(qairNo), STDatasets(6)%VarNames(TmaxNo),&
    STDatasets(7)%VarNames(TminNo), STDatasets(8)%VarNames(uwindNo),&
    STDatasets(9)%VarNames(vwindNo), STDatasets(10)%VarNames(fdiffNo))

  STDatasets(1)%VarNames = rainNames(1:rainNo)
  STDatasets(2)%VarNames = lwdnNames(1:lwdnNo)
  STDatasets(3)%VarNames = swdnNames(1:swdnNo)
  STDatasets(4)%VarNames = presNames(1:presNo)
  STDatasets(5)%VarNames = qairNames(1:qairNo)
  STDatasets(6)%VarNames = TmaxNames(1:TmaxNo)
  STDatasets(7)%VarNames = TminNames(1:TminNo)
  STDatasets(8)%VarNames = uwindNames(1:uwindNo)
  STDatasets(9)%VarNames = vwindNames(1:vwindNo)
  STDatasets(10)%VarNames = fdiffNames(1:fdiffNo)
END SUBROUTINE read_variable_names

!------------------------------------------------------------------------------!

SUBROUTINE read_landmask(LandmaskFile, CRU)
  !*## Purpose
  !
  ! Reads a landmask and allocates the data arrays in the ```CRU``` derived
  ! type, sets up indexing arrays for reading from met forcing datasets and
  ! global landmask variables.
  !
  !*## Method
  !
  ! Uses the NetCDF dataset at ```LandmaskFile``` to determine which cells are
  ! computed in this run of CABLE.

  ! Read the landmask file to determine which land cells we're looking at.
  ! Define the input variables
  CHARACTER(len=256), INTENT(IN)  :: LandmaskFile
  TYPE(cru_type), INTENT(INOUT)   :: CRU

  ! Assign an integer to store the file ID, and to record the error status
  INTEGER       :: FileID, ok

  ! And some IDs to store the latitude, longitude and mask variable IDs
  INTEGER       :: LatID, LonID, LandID

  ! And then the actual values along the axes
  REAL, DIMENSION(:), ALLOCATABLE :: Latitudes, Longitudes

  ! Need a counter and iterators when iterating through the points in the land mask
  INTEGER       :: MaskCounter, LatIndx, LonIndx, VarIndx

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: LandMask

  ! Now we attempt to read the netCDF file
  ok = NF90_OPEN(TRIM(LandmaskFile), nf90_nowrite, FileID)
  CALL handle_err(ok, "Error opening landmask file at "//TRIM(LandmaskFile))

  ! Inquire about the latitude and longitude dimensions
  ok = NF90_INQ_DIMID(FileID, 'latitude', LatID)
  CALL handle_err(ok, "Error finding latitude DIMID.")

  ok = NF90_INQUIRE_DIMENSION(FileID, LatID, LEN = yDimSize)
  CALL handle_err(ok, "Error inquiring latitude dimension.")

  ok = NF90_INQ_DIMID(FileID, 'longitude', LonID)
  CALL handle_err(ok, "Error finding longitude DIMID.")

  ok = NF90_INQUIRE_DIMENSION(FileID, LonID, LEN = xDimSize)
  CALL handle_err(ok, "Error inquiring longitude dimension.")

  ! And read the dimensions into arrays
  ALLOCATE(Latitudes(yDimSize))
  ok = NF90_INQ_VARID(FileID, 'latitude', LatID)
  call handle_err(ok, "Error finding latitude VARID.")
  ok = NF90_GET_VAR(FileID, LatID, Latitudes)
  call handle_err(ok, "Error reading latitude variable.")

  ALLOCATE(Longitudes(xDimSize))
  ok = NF90_INQ_VARID(FileID, 'longitude', LonID)
  call handle_err(ok, "Error finding longitude VARID.")
  ok = NF90_GET_VAR(FileID, LonID, Longitudes)
  call handle_err(ok, "Error reading longitude variable.")

  ! Now read the mask, this mask is an integer (0/1) version
  ALLOCATE(mask(xDimSize, yDimSize))
  ALLOCATE(LandMask(xDimSize, yDimSize))
  ALLOCATE(CRU%LandMask(xDimSize, yDimSize))
  ok = NF90_INQ_VARID(FileID, 'land', LandID)
  CALL handle_err(ok, "Error finding land VARID.")
  ok = NF90_GET_VAR(FileID, LandID, LandMask)
  CALL handle_err(ok, "Error reading land variable.")

  ! Finally close the mask file
  ok = NF90_CLOSE(FileID)
  CALL handle_err(ok, "Error closing landmask file.")

  ! For some reason we also want a logical version
  !CRU%LandMask = to_logical(mask)
  WHERE (LandMask > 0)
    CRU%LandMask = .TRUE.
    mask = 1
  ELSEWHERE
    CRU%LandMask = .FALSE.
    mask = 0
  END WHERE

  ! Set the total number of land cells by counting the TRUEs in the landmask.
  CRU%mLand = COUNT(CRU%LandMask)

  ! Now we know how many elements in our specific region, allocate global (!!!)
  ! arrays for latitude and longitude value/indices.
  ALLOCATE(Latitude(CRU%mLand), Longitude(CRU%mLand))
  ALLOCATE(Land_y(CRU%mLand), Land_x(CRU%mLand))

  ! Now pull out the coordinates in our current mask
  MaskCounter = 1
  LoopThroughLats: DO LatIndx = 1, yDimSize
    LoopThroughLons: DO LonIndx = 1, xDimSize
      IF (CRU%LandMask(LonIndx, LatIndx)) THEN
        Land_x(MaskCounter) = LonIndx
        Land_y(MaskCounter) = LatIndx
        Longitude(MaskCounter) = Longitudes(LonIndx)
        Latitude(MaskCounter) = Latitudes(LatIndx)
        MaskCounter = MaskCounter + 1
      END IF
    END DO LoopThroughLons
  END DO LoopThroughLats

  ! Allocate memory for the actual meteorological forcing data
  AllocateMemory: DO VarIndx = 1, SIZE(CRU%Met)
    ALLOCATE(CRU%Met(VarIndx)%MetVals(CRU%mLand))
  END DO AllocateMemory

  ! And do the same for nitrogen deposition and atmospheric CO2
  ALLOCATE(CRU%NDepVals(CRU%mLand))
  ALLOCATE(CRU%CO2Vals(CRU%mLand))

  ! Allocate global data- this practice must be purged in future.
  metGrid = "mask"
  mLand = CRU%mLand
  nMetPatches = 1

  ALLOCATE(Lat_all(xDimSize, yDimSize), Lon_all(xDimSize, yDimSize))
  FillLatitudes: DO LonIndx = 1, xDimSize
    Lat_all(LonIndx, :) = Latitudes
  END DO FillLatitudes

  FillLongitudes: DO LatIndx = 1, yDimSize
    Lon_all(:, LatIndx) = Longitudes
  END DO FillLongitudes

  ! Some CABLE time units?
  shod = 0.
  sdoy = 1
  smoy = 1
  syear = CRU%CYEAR

  ! Allocate memory for the average long-wave radiation down for some reason?
  ALLOCATE(CRU%AVG_LWDN(mLand))
END SUBROUTINE read_landmask

!!------------------------------------------------------------------------------!

SUBROUTINE prepare_temporal_dataset(FileName, TargetArray)
  !*## Purpose
  !
  ! Read a temporally varying dataset that exists as a text file.
  !
  !## Method
  !
  ! Read the ```FileName``` text file which contains the yearly CO2 data in
  ! ```YEAR VALUE``` pairs into ```TargetArray```.

  ! If we have a dataset that's only varying in time, then we should be able to
  ! just the whole thing into memory immediately.

  CHARACTER(LEN=256), INTENT(IN)    :: FileName
  REAL, DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: TargetArray

  ! A unique file ID
  INTEGER   :: FileID

  ! A placeholder for the line, that we inspect to check for comment lines
  CHARACTER(LEN=256)    :: LineInFile

  ! Placeholders for the Key in the line, and years which form the start and end
  ! indexes of the array
  INTEGER               :: Time, StartTime, EndTime
  REAL                  :: DummyValue

  ! The status holder, and the file and header counters
  INTEGER :: ios, LineCounter, HeaderCounter, iter

  ! We're going to assume that we only have a single file.
  ! We need to iterate through the file twice, first to inspect the number of
  ! entries in the file to allocate the correct amount of memory, second to
  ! actually write the data to the array.
  CALL get_unit(FileID)
  OPEN(FileID, FILE = FileName, STATUS = "old", ACTION = "read", IOSTAT = ios)
  IF (ios < 0) THEN
    WRITE(*,*) "Open of temporal dataset file failed with status:", ios
    CALL EXIT(5)
  END IF

  LineCounter = 0
  HeaderCounter = 0

  DetermineSize: DO
    ! Read line by line, checking for a header line
    READ(FileID, '(A)', IOSTAT = ios) LineInFile

    IF (ios < 0) THEN
      ! Read completed successfully and hit EOF

      ! Take the current time value and set it as the end index of the array
      EndTime = Time

      EXIT DetermineSize
    END IF

    ! Check for comment characters- assume they're all at the top of the file
    IF ((LineInFile(1:1) == "#") .OR. (LineInFile(1:1) == "!")) THEN
      HeaderCounter = HeaderCounter + 1
    ELSE
      ! Otherwise, its a line of useful data
      LineCounter = LineCounter + 1

      ! Read the line of data as a (time, value) pair. We don't want to use the
      ! value this time around because we haven't yet allocated the array to
      ! store it.
      READ(LineInFile, FMT = '(I4, F)') Time, DummyValue
    END IF

    ! If the LineCounter is 1, we just read the first line of data, so get the
    ! starting index for our array
    IF (LineCounter == 1) THEN
      StartTime = Time
    END IF
  END DO DetermineSize

  ! Rewind the pointer to the start of the file
  REWIND(FileID)

  ! Now allocate the appropriate memory
  ALLOCATE(TargetArray(StartTime:EndTime))

  ! First step over the header lines
  SkipHeader: DO iter = 1, HeaderCounter
    READ(FileID, FMT = '(A)', IOSTAT = ios) LineInFile
  END DO SkipHeader

  ! Now the useful contents of the file into the array
  ReadValues: DO iter = StartTime, EndTime
    READ(FileID, FMT = '(A)', IOSTAT = ios) LineInFile
    READ(LineInFile, FMT = '(I4, F)', IOSTAT = ios) Time, TargetArray(iter)
  END DO ReadValues
  CLOSE(FileID)
END SUBROUTINE prepare_temporal_dataset

!------------------------------------------------------------------------------!

SUBROUTINE get_cru_co2(CRU, CO2Air)
  !*## Purpose
  !
  ! Read the CO2 value from a specific year.
  !
  !## Method
  !
  ! Get the CO2 values from the already-defined CO2 data array and fill in the
  ! ```CO2Air``` array with this data. The year is dependent on the
  ! ```CRU%CO2Method```.

  ! Get CO2 values from the MET data.
  TYPE(CRU_TYPE), INTENT(IN)                  :: CRU
  REAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: CO2Air

  ! We treat the CO2 as an array always, even in the instances where we're using
  ! a single global value. We pass to the MET struct as an array anyway, so we
  ! don't really get any benefit by treating it as a single value.

  ! We may want to take a constant CO2 value from a particular year, so bind the
  ! year we use here to another variable
  INTEGER         :: CO2Year

  ! What method are we using for CO2?
  IF (CRU%CO2Method == "Yearly") THEN
    ! We use the same year as the simulation year
    CO2Year = CRU%cYear
    
    ! Make sure it's in the bounds of the array
    CO2Year = MAX(LBOUND(CRU%CO2Vals, 1), CO2Year)
    CO2Year = MIN(UBOUND(CRU%CO2Vals, 1), CO2Year)

    CO2Air(:) = CRU%CO2Vals(CO2Year)
  ELSE
    ! The user has specified the year
    READ(CRU%CO2Method, '(I)') CO2Year
    CO2Air(:) = CRU%CO2Vals(CO2Year)
  END IF
END SUBROUTINE get_cru_co2

!------------------------------------------------------------------------------!

SUBROUTINE get_cru_ndep(CRU)
  !*## Purpose
  !
  ! Read the NDep value from a specific year.
  !
  !## Method
  !
  ! Get the NDep values using the ```SPATIO_TEMPORAL_DATASET``` associated with
  ! NDep and fill in the NDep data array.

  ! Get the nitrogen deposition values from the MET data.
  TYPE(CRU_TYPE), INTENT(INOUT)             :: CRU

  ! The nitrogen deposition is a spatio_temporal_dataset netCDF dataset.
  ! The data is unfortunately arranged, in that the frequency is yearly, but the
  ! time index is monthly with month fractions i.e. 5.922..., 17.922... which
  ! makes it very inconvenient to index. For now, we have to rely on the prior
  ! knowledge that the dataset begins at 1850.
  INTEGER         :: TimeIndex

  ! Temporary array to store the full set of spatial data
  REAL, DIMENSION(:, :), ALLOCATABLE        :: TmpArray

  ! Iterator variable
  INTEGER  :: GridCell

  ! Status checker
  INTEGER  :: ok

  ! Which method are we using for Nitrogen deposition?
  IF (CRU%NDepMethod == "Yearly") THEN
    ! We're using the simulation year for the deposition data
    TimeIndex = CRU%cYear - 1849
  ELSE
    ! We're using a specific year of data
    READ(CRU%NDepMethod, '(I4)') TimeIndex
    TimeIndex = TimeIndex - 1849
  END IF

  ! If the year is before the time period of the data, get the earliest
  ! year.
  TimeIndex = MAX(1, TimeIndex)

  ! And finally grab the data. Since we're using landmask, we have to extract
  ! the full set of data first and then read the unmasked points into the NDep
  ! values.
  IF (CRU%DirectRead) THEN

    ApplyLandmaskDirect: DO GridCell = 1, CRU%mLand
      ok = NF90_GET_VAR(CRU%NDepFID, CRU%NDepVID, CRU%NDepVals(GridCell),&
      START = (/land_x(GridCell), land_y(GridCell), TimeIndex/))
      CALL handle_err(ok, "Reading from NDep")
    END DO ApplyLandMaskDirect
  ELSE

    ALLOCATE(TmpArray(CRU%xDimSize, CRU%yDimSize))

    ok = NF90_GET_VAR(CRU%NDepFID, CRU%NDepVID, TmpArray,&
      START = (/1, 1, TimeIndex/), COUNT = (/CRU%xDimSize, CRU%yDimSize, 1/))
    CALL handle_err(ok, "Reading from NDep")

    ApplyLandmaskIndirect: DO GridCell = 1, CRU%mLand
      CRU%NDepVals(GridCell) = TmpArray(land_x(GridCell), land_y(GridCell))
    END DO ApplyLandmaskIndirect
  END IF


END SUBROUTINE get_cru_ndep

!!------------------------------------------------------------------------------!

ELEMENTAL FUNCTION Esatf(TC)
! ------------------------------------------------------------------------------
! At temperature TC [deg C], return saturation water vapour pressure Esatf [mb]
! from Teten formula.
! MRR, xx/1987
! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
!                 just like intrinsic functions.
! MRR, 12-mar-02: Convert Qsatf (specific humidity routine) to Esatf
! ------------------------------------------------------------------------------
IMPLICIT NONE

REAL(sp), INTENT(IN) :: TC          ! temp [deg C]
REAL(sp)             :: Esatf       ! saturation vapour pressure [mb]

REAL(sp) :: TCtmp                   ! local
REAL(sp),PARAMETER:: A = 6.106      ! Teten coefficients
REAL(sp),PARAMETER:: B = 17.27      ! Teten coefficients
REAL(sp),PARAMETER:: C = 237.3      ! Teten coefficients

TCtmp = TC                            ! preserve TC
IF (TCtmp > 100.0) TCtmp = 100.0   ! constrain TC to (-40.0, 100.0)
IF (TCtmp < -40.0) TCtmp = -40.0

Esatf = A*exp(B*TCtmp/(C+TCtmp))    ! sat vapour pressure (mb)

END FUNCTION Esatf

SUBROUTINE cru_close(CRU)
  !*## Purpose
  !
  ! Close the NetCDF file at the end of the run.
  !
  !## Method
  !
  ! Use the ```CurrentFileID``` to determine the currently open NetCDF file in
  ! the dataset and then close the file attached to that ID.

  ! Close the open NetCDF files
  TYPE(CRU_TYPE), INTENT(INOUT) :: CRU

  ! Iterator
  INTEGER  :: iter, ok

  CloseFiles: DO iter = 1, SIZE(CRU%MetDatasets)
    ok = NF90_CLOSE(CRU%MetDatasets(iter)%CurrentFileID)
  END DO CloseFiles
END SUBROUTINE cru_close

END MODULE CABLE_CRU
