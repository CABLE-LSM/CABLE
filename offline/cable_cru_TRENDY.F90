! Author: Lachlan Whyborn
! Last Modified: Tue 14 May 2024 14:17:11

MODULE CABLE_CRU

! Import from other modules. We should absolutely not be importing data from
! other modules, like those from cable_io_vars_module, but this behaviour
! permeates the code and will take months of work to refactor out and will be
! addressed when CABLE receives a fundamental rewrite.
! For some reason, some bits of data are imported here, and others are imported
! within subroutines. We've combined into a single import at the start of the
! module.

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

!TYPE SPATIO_TEMPORAL_DATASET
  !! When we actually want a piece of data, we pass this along with the year
  !! and day to a read routine. The routine checks the passed year against
  !! StartYear(CurrentFileIndx) and EndYear(CurrentFileIndx), and decides whether
  !! to either a) continue with the current file if it falls within that range,
  !!        or b) move to a new file if it falls outside that range.
  !! By storing the ID associated with the currently open file, we only need to
  !! open a new IO stream a couple of times per simulation, rather than every
  !! time we start a new year or a new time step.

  !! This data type is used to store effectively a reference to spatio-temporal
  !! datasets which span multiple files, which is common for external forcings.
  !! The idea is to store the list of files in the dataset, and their start and
  !! end years respectively. We also need to know the internal variable names in
  !! the netCDF dataset so we can call it correctly. 
  !! This can possibly permit series with different sampling intervals as well.
  !CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: FileNames
  !INTEGER, DIMENSION(:), ALLOCATABLE            :: StartYear, EndYear
  !CHARACTER(LEN=16), DIMENSION(:), ALLOCATABLE  :: VarNames

  !! As well as the metadata about the dataset contents, we want some data that
  !! is updated to inform which file is currently open and it's period. If open a
  !! file relevant to a particular year, and move onto the next year that's still
  !! contained in the current file's period, we don't want to go through open
  !! another file, we just want to continue using the current file.
  !! To do this, we need to know the current FileID, and which index it is in
  !! our list of files so we can access it's period.
  !INTEGER :: CurrentFileID, CurrentVarId, CurrentFileIndx = 0
!END TYPE SPATIO_TEMPORAL_DATASET

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
  ! A rewrite of the CRU_INIT subroutine in cable_cru_TRENDY.f90

  ! The main goal of this routine is to mutate the CRU data structure to prep
  ! it for the experiment.
  TYPE(CRU_TYPE), INTENT(OUT)    :: CRU

  ! Start with the things we want from the namelist. The namelist must set
  ! the filenames to read from, the method of choosing atmospheric carbon and
  ! nitrogen deposition, and the timestep.
  ! We want one filename for each variable, stored in a list predefined index.
  CHARACTER(LEN=256), dimension(13) :: InputFiles

  ! To determine the data we want from the Met data, we need to look at the
  ! landmask
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: LandMask
 
  ! Iterator variable for the variables
  INTEGER   :: VarIndx

  ! Checker for the NetCDF io
  INTEGER   :: ok
  LOGICAL :: IsOpen

  ! Read the namelist and pull out the input files and configuration variables
  CALL read_MET_namelist_cbl(InputFiles, CRU)

  ! Read the landmask and allocate appropriate memory for the array variables
  CALL read_landmask(InputFiles(13), CRU)

  ! We know that the first 9 Met variables (indices 1-9) are always going to be
  ! time series, so always build their keys.
  BuildKeys: DO VarIndx = 1, 9
    CALL prepare_spatiotemporal_dataset(InputFiles(VarIndx), CRU%MetDatasets(VarIndx))
  END DO BuildKeys

  ! Our building of NDep, Carbon and FDiff keys depend on our choices of method.
  IF (CRU%ReadDiffFrac) THEN
    CALL prepare_spatiotemporal_dataset(InputFiles(10), CRU%MetDatasets(10))
  END IF

  ! Set the possible variable names for the main Met variables
  CALL read_variable_names(CRU%MetDatasets)

  ! Open the datasets at the first file so we don't need CALL1 behaviour later
  InitialiseDatasets: DO VarIndx = 1, 9
    CALL open_at_first_file(CRU%MetDatasets(VarIndx))
  END DO InitialiseDatasets

  IF (CRU%ReadDiffFrac) THEN
    CALL open_at_first_file(CRU%MetDatasets(10))
  END IF

  ! Set up the carbon reader
  IF (CRU%CO2Method == "Yearly") THEN
    ! If it's year by year, or single year atmospheric CO2, we want to read the
    ! temporal dataset.
    CALL prepare_temporal_dataset(InputFiles(11), CRU%CO2Vals)
  ELSEIF (CRU%CO2Method == "Spatial") THEN
    ! This is a placeholder, not currently implemented
  ELSE
    CALL prepare_temporal_dataset(InputFiles(11), CRU%CO2Vals)
  END IF

  ! Set up the nitrogen deposition reader
  IF (CRU%NDepMethod == "Yearly") THEN
    CALL prepare_temporal_dataset(InputFiles(12), CRU%NDepVals)
  ELSEIF (CRU%NDepMethod == "Spatial") THEN
    CALL prepare_spatiotemporal_dataset(InputFiles(12), CRU%NDepDataset)
    ! For now, set the file index to 1 since we know its only one file
    ! We don't want this to be a long term solution, just a temporary practical
    ! solution to get the LUC and Demography running with ACCESS Met forcing.
    CRU%NDepDataset%CurrentFileIndx = 1

    ok = NF90_OPEN(CRU%NDepDataset%Filenames(1), NF90_NOWRITE,&
      CRU%NDepFID)
    CALL handle_err(ok, "Reading NDep file")

    ok =NF90_INQ_VARID(CRU%NDepFID, "N_deposition",&
      CRU%NDepVID)
    CALL handle_err(ok, "Reading NDep variable")
  ELSE
    CALL prepare_spatiotemporal_dataset(InputFiles(12), CRU%NDepDataset)
    ! For now, set the file index to 1 since we know its only one file
    CRU%NDepDataset%CurrentFileIndx = 1
    
    ok = NF90_OPEN(CRU%NDepDataset%Filenames(1), NF90_NOWRITE,&
      CRU%NDepFID)
    CALL handle_err(ok, "Reading NDep file")

    ok = NF90_INQ_VARID(CRU%NDepFID, "N_deposition",&
      CRU%NDepVID)
    CALL handle_err(ok, "Reading NDep variable")
  END IF
END SUBROUTINE CRU_INIT


SUBROUTINE CRU_GET_SUBDIURNAL_MET(CRU, MET, CurYear, ktau, kend)

  ! Obtain one day of CRU-NCEP meteorology, subdiurnalise it using a weather
  ! and return the result to the CABLE driver.

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
     !       met%precip_sn(is:ie) = (1.0 - (54.62 - 0.2 *(WG%Temp(iland) + 273.15)))* met%precip(is:ie)trailer made it look like a Godzilla story told from the personal perspective of Bryan Cranston as an 'everyman' type character. But then
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

  ! One of the tricky things here is the requirement to have the previous and
  ! next days Tmax and Tmin. This is handled by the dataset reader, if we try to
  ! read outside the range of our data, we just use the first/last record.

  ! The year of met forcing we use depends on our choice of configuration.
  ! Sometimes, recycle through a subset of data, and others we use the sim year.
  ! So we want to store both.
  INTEGER   :: RecycledYear, MetYear

  ! We want to handle the nextTmin and prevTmax things separately so we don't
  ! mess with the data stored in CRU
  INTEGER   :: DummyYear, DummyDay

  ! Define iteration variable
  INTEGER   :: VarIndx 

  ! Determine the recycled and sim year
  RecycledYear = 1901 + MOD(CRU%cYear - 1501, CRU%metRecyc)

  ! Iterate through the base variables
  IterateVariables: DO VarIndx = 1, 9
    ! If the variable is time recycled
    IF (CRU%isRecycled(VarIndx)) THEN
      MetYear = RecycledYear
    ELSE
      MetYear = CRU%cYear
    END IF

    ! CRU%CTStep is the current day, which we use to index the netCDF arrays
    CALL read_metvals(CRU%MetDatasets(VarIndx), CRU%Met(VarIndx)%MetVals,&
      land_x, land_y, MetYear, CRU%CTStep, CRU%LeapYears)
  END DO IterateVariables

  ! Only read the DiffFrac sometimes, with the same process as the rest of the
  ! met variables
  IF (CRU%ReadDiffFrac) THEN
    IF (CRU%isRecycled(10)) THEN
      MetYear = RecycledYear
    ELSE
      MetYear = CRU%cYear
    END IF

    CALL read_metvals(CRU%MetDatasets(10), CRU%Met(10)%MetVals, land_x, land_y,&
      MetYear, CRU%CTStep, CRU%LeapYears)
  END IF

  ! Now the variables with special handling- nextTmin and prevTmax
  ! The easiest way to do this is to simply change the day by 1, check if we
  ! need to change the year, and go through the same read process. This is a
  ! little inefficient as we'll be messing with the Tmin and Tmax dataset
  ! structs and in select instances (at the end of the era of a particular file)
  ! opening and closing the io stream unnecessarily, but I think it's a minor
  ! evil.

  ! Address prevTmax first
  ! Check what the day is, and whether we need to change the year
  IF (CRU%CTSTEP == 1) THEN
    ! Go back to previous year
    DummyYear = CRU%cYear - 1
    IF (CRU%LeapYears) THEN
      ! Check what the day should be
      IF (is_leapyear(DummyYear)) THEN
        DummyDay = 366
      ELSE
        DummyDay = 365
      END IF
    ELSE
      ! No leapyears
      DummyDay = 365
    END IF
  ELSE
    ! Not the first day of the year, treat normally
    DummyYear = CRU%cYear
    DummyDay = CRU%CTStep - 1
  END IF

  ! Was the Tmax recycled?
  IF (CRU%isRecycled(6)) THEN
    MetYear = 1901 + MOD(DummyYear - 1501, CRU%metRecyc)
  ELSE
    MetYear = CRU%cYear
  END IF

  ! Now we just need to call cru_read_metvals with the Tmax Dataset reader and
  ! the prevTmax array to write to
  CALL read_metvals(CRU%MetDatasets(6), CRU%Met(11)%MetVals, land_x, land_y,&
    MetYear, DummyDay, CRU%LeapYears)

  ! Now do nextTmin
  ! Check what the day is, and whether we need to change the year
  ! This changes whether its a leapyear or not
  IF ((CRU%LeapYears) .AND. (CRU%CTStep == 366)) THEN
    ! We're at the last day of a leapyear
    DummyDay = 1
    DummyYear = CRU%cYear + 1
  ELSEIF (CRU%CTStep == 365) THEN
    ! We're at the end of a normal year
    DummyDay = 1
    DummyYear = CRU%cYear + 1
  ELSE
    DummyDay = CRU%CTStep + 1
    DummyYear = CRU%cYear
  END IF

  ! Was the Tmin recycled?
  IF (CRU%isRecycled(7)) THEN
    MetYear = 1901 + MOD(DummyYear - 1501, CRU%metRecyc)
  ELSE
    MetYear = CRU%cYear
  END IF

  CALL read_metvals(CRU%MetDatasets(7), CRU%Met(12)%MetVals, land_x, land_y,&
    MetYear, DummyDay, CRU%LeapYears)
END SUBROUTINE cru_get_daily_met

SUBROUTINE read_MET_namelist_cbl(InputFiles, CRU)
  ! # read_MET_namelist_cbl(InputFiles, CRU)
  ! Read the supplied cru.nml namelist and extract the desired input file names
  ! and a series of configuration options regarding the usage of the MET data.
  ! ## Namelist inputs
  ! The entries expected in the namelist are:
  ! ```
  ! &MetConfig
  ! rainFile = filename_template
  ! lwdnFile = ...
  ! swdnFile = ...
  ! presFile = ...
  ! qairFile = ...
  ! TmaxFile = ...
  ! TminFile = ...
  ! uwindFile = ...
  ! vwindFile = ...
  ! fdiffFile = ...
  ! carbonFile = ...
  ! NDepFile = ...
  ! CO2Method = "yearly", "spatial", <year>
  ! NDepMethod = "yearly", "spatial", <year>
  ! ReadDiffFrac = T/F
  ! LeapYears = T/F
  ! MetRecyc = <NYears>
  ! ! And then variables for the recycling of MET variables.
  ! rainRecyc = T/F
  ! lwdnRecyc = T/F
  ! ! and so on.
  ! ```

  ! It's not clear why this recycling is necessary, but it's included in the
  ! TRENDY configuration, so keep it for imply having the repconsistency.

  ! The filename template is "path/to/file/prefix<startdate>_<enddate>suffix".
  ! For example, there may be a series of files describing the rain data at
  ! ```"data/rain/cru_rain_data_18500101_18991231_av.nc"```,
  ! ```"data/rain/cru_rain_data_19000101_19491231_av.nc"``` and so on.
  ! The entry for rainFile would be 
  ! ```"data/rain/cru_rain_data_<startdate>_<enddate>_av.nc".
  !
  ! The CO2Method and NDepMethod are "yearly", with globally constant values
  ! for each year, "spatial" with spatially varying values for each year, or
  ! a specific year which fixes the value at that specific year.

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
  LOGICAL             :: ReadDiffFrac, LeapYears

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
  CO2Method = "yearly"
  NDepMethod = "yearly"
  MetRecyc = 20
  ReadDiffFrac = .TRUE.

  ! Set up and read the namelist
  NAMELIST /MetConfig/ rainFile, lwdnFile, swdnFile, presFile, qairFile,&
                    TmaxFile, TminFile, uwindFile, vwindFile, fdiffFile,&
                    CO2File, NDepFile, landmaskFile,&
                    rainRecycle, lwdnRecycle, swdnRecycle, presRecycle,&
                    qairRecycle, TmaxRecycle, TminRecycle, uwindRecycle,&
                    vwindRecycle, fdiffRecycle,&
                    ReadDiffFrac, CO2Method, NDepMethod, MetRecyc, LeapYears,&
                    DtHrs

  ! Get a temporary unique ID and use it to read the namelist
  CALL get_unit(nmlUnit)
  OPEN(nmlUnit, file = "cru.nml", status = 'old', action = 'read')
  READ(nmlUnit, nml = MetConfig)
  CLOSE(nmlUnit)

  ! Now pack the filepaths into the data structure we want to transport around
  InputFiles(1) = rainFile
  InputFiles(2) = lwdnFile
  InputFiles(3) = swdnFile
  InputFiles(4) = presFile
  InputFiles(5) = qairFile
  InputFiles(6) = TmaxFile
  InputFiles(7) = TminFile
  InputFiles(8) = uwindFile
  InputFiles(9) = vwindFile
  InputFiles(10) = fdiffFile
  InputFiles(11) = CO2File
  InputFiles(12) = NDepFile
  InputFiles(13) = landmaskFile

  ! Set the recycling booleans in the CRU struct
  CRU%IsRecycled(1) = rainRecycle
  CRU%IsRecycled(2) = lwdnRecycle
  CRU%IsRecycled(3) = swdnRecycle
  CRU%IsRecycled(4) = presRecycle
  CRU%IsRecycled(5) = qairRecycle
  CRU%IsRecycled(6) = TmaxRecycle
  CRU%IsRecycled(7) = TminRecycle
  CRU%IsRecycled(8) = uwindRecycle
  CRU%IsRecycled(9) = vwindRecycle
  CRU%IsRecycled(10) = fdiffRecycle

  ! Bind the remaining config variables to the CRU structure
  CRU%CO2Method = CO2Method
  CRU%NDepMethod = NDepMethod
  ! Convert the hourly timestep to seconds
  CRU%DtSecs = int(DtHrs * 3600.)
  CRU%MetRecyc = MetRecyc
  CRU%ReadDiffFrac = ReadDiffFrac
  CRU%LeapYears = LeapYears
END SUBROUTINE read_MET_namelist_cbl

!------------------------------------------------------------------------------!

SUBROUTINE read_variable_names(STDatasets)
  ! Read in the possible variable names from a namelist, and assign them to the
  ! SpatioTemporalDataset structs.
  ! While the CF convention exists to dictate what specific variables should be
  ! named in a dataset, this convention is yet to hold much sway and the
  ! variables in a NetCDF dataset may be named any number of things.
  ! The current solution is to have an internally maintained namelist which
  ! contains the possible NetCDF names for a given variable, and hope that there
  ! is some overlap between names across macro datasets.
  ! Having an internal namelist was chosen because it strikes a balance between
  ! putting too much onus on the user, and having such options set opaquely deep
  ! within the source code. If a user comes with a new macro dataset, with a new
  ! set of variable names, they should not have to dig deep into the source code
  ! to make an addition to the list of possible NetCDF names for a variable.
  ! There are two alternatives to the current approach:
  ! 1. Have the user simply state the variable names used in the NetCDF data.
  !   This may be a better long term option, if it turns out there is minimal
  !   overlap in the namings.
  ! 2. Have a list of possible names defined explicitly in the code. This is not
  !   desirable due to the reasons outlined above, in that it makes the process
  !   opaque and difficult to adapt for a new user.

  ! The namelist contains both a list of possible names for an internal variable
  ! and the number of names in this list. This allows us to first allocate 
  ! memory for the list of possible names, then to write to those arrays.

  ! The more I think about this the more I think going with option 1 is the
  ! better long term option, but let's run with this for now.

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
  ! Read the landmask file to determine which land cells we're looking at.
  ! Define the input variables
  CHARACTER(len=256), INTENT(IN)  :: LandmaskFile
  TYPE(cru_type), INTENT(INOUT)   :: CRU

  ! Assign an integer to store the file ID, and to record the error status
  INTEGER       :: FileID, ErrStatus

  ! And some IDs to store the latitude, longitude and mask variable IDs
  INTEGER       :: LatID, LonID, LandID

  ! Need the lengths of the axes
  !INTEGER       :: xDimSize, yDimSize

  ! And then the actual values along the axes
  REAL, DIMENSION(:), ALLOCATABLE :: Latitudes, Longitudes

  ! Need a counter and iterators when iterating through the points in the land mask
  INTEGER       :: MaskCounter, LatIndx, LonIndx, VarIndx

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: LandMask

  ! Now we attempt to read the netCDF file
  ErrStatus = NF90_OPEN(TRIM(LandmaskFile), nf90_nowrite, FileID)
  IF (ErrStatus /= NF90_NOERR) THEN
    CALL NC_ABORT(ErrStatus, "Error reading the landmask file.")
  END IF

  ! Inquire about the latitude and longitude dimensions
  ErrStatus = nf90_inq_dimid(FileID, 'latitude', LatID)
  IF (ErrStatus /= NF90_NOERR) THEN
    CALL NC_ABORT(ErrStatus, "Error reading the latitude ID.")
  END IF
  ErrStatus = nf90_inquire_dimension(FileID, LatID, LEN = yDimSize)
  IF (ErrStatus /= NF90_NOERR) THEN
    CALL NC_ABORT(ErrStatus, "Error reading the latitude length.")
  END IF

  ErrStatus = nf90_inq_dimid(FileID, 'longitude', LonID)
  IF (ErrStatus /= NF90_NOERR) THEN
    CALL NC_ABORT(ErrStatus, "Error reading the longitude ID.")
  END IF
  ErrStatus = nf90_inquire_dimension(FileID, LonID, LEN = xDimSize)
  IF (ErrStatus /= NF90_NOERR) THEN
    CALL NC_ABORT(ErrStatus, "Error reading the longitude length.")
  END IF

  ! And read the dimensions into arrays
  ALLOCATE(Latitudes(yDimSize))
  ErrStatus = nf90_inq_varid(FileID, 'latitude', LatID)
  ErrStatus = nf90_get_var(FileID, LatID, Latitudes)

  ALLOCATE(Longitudes(xDimSize))
  ErrStatus = nf90_inq_varid(FileID, 'longitude', LonID)
  ErrStatus = nf90_get_var(FileID, LonID, Longitudes)

  ! Now read the mask, this mask is an integer (0/1) version
  ALLOCATE(mask(xDimSize, yDimSize))
  ALLOCATE(LandMask(xDimSize, yDimSize))
  ALLOCATE(CRU%LandMask(xDimSize, yDimSize))
  ErrStatus = nf90_inq_varid(FileID, 'land', LandID)
  ErrStatus = nf90_get_var(FileID, LandID, LandMask)

  ! Finally close the mask file
  ErrStatus = nf90_close(FileID)

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

!------------------------------------------------------------------------------!

!SUBROUTINE prepare_spatiotemporal_dataset(FileTemplate, Dataset)
  !! Prepare the met netCDF files for reading. We follow a process similar to
  !! that of JULES to handle time series data split over multiple files, see
  !! 6.31.2.2 at https://jules-lsm.github.io/latest/namelists/drive.nml.html.

  !! The user specifies a template matching the file names for each variable,
  !! containing the start and end date of the data in the file. For example,
  !! rain data files may be named "pre/precip_1850101_18991231.nc",
  !! "pre/precip_19000101_19491231.nc" and so on. Then the user would specify
  !! rainFile = "pre/precip_<startdate>_<enddate>.nc" in the namelist.

  !! We then use the ls command to find all files matching this pattern and 
  !! write the output to a text file. We then inspect the filenames to
  !! determine the time period for each filej.

  !! We then create a data structure which effectively holds metadata about
  !! the possibly multi-file dataset.
  !! The data structure contains the list of files, their periods, as well as
  !! a reference to the current open file.
  !! Specify the intent of the arguments
  !CHARACTER(len=256), INTENT(IN)  :: FileTemplate
  !TYPE(SPATIO_TEMPORAL_DATASET), INTENT(OUT)  :: Dataset

  !! We need to have a definition of the substrings we're going to replace.
  !CHARACTER(len=11) :: StartDate = "<startdate>"
  !CHARACTER(len=9)  :: EndDate = "<enddate>"

  !! For clarity in the looping, and so we don't mutate the original inputs,
  !! we extract the filename from the array.
  !CHARACTER(len=256):: CurrentFile

  !! We write the outputs from the `ls` command to a specific reserved filename
  !CHARACTER(len=32) :: ReservedOutputFile = "__FileNameWithNoClashes__.txt"

  !! We read the start and end years to strings before writing them to the key.
  !CHARACTER(len=4)  :: StartYear, EndYear

  !! Initialise integers to store the status of the command.
  !INTEGER           :: ExStat, CStat

  !! We're going to want to remember where "<startdate>" and "<enddate>" occur
  !! in the filename templates, so we can retrieve the time periods of each
  !! file later.
  !INTEGER           :: IndxStartDate = 0, IndxEndDate = 0

  !! INTEGERs for the file ID that reads the output from the 'ls' command.
  !INTEGER           :: InputUnit

  !! And status integers for the IO
  !INTEGER           :: ios

  !! Need to compute the number of files in the dataset
  !INTEGER           :: FileCounter

  !! Iterators
  !INTEGER           :: CharIndx, FileIndx

  !! Get a unique ID here.
  !CALL get_unit(InputUnit)

  !! First thing we have to do is determine the location of the <startdate>
  !! and <enddate> strings in the supplied filenames. To do that, we compare
  !! substrings of the relevant length, with the substring moving along the
  !! length of the filename string.
  !! We iterate separately for <startdate> and <enddate> so that we make
  !! absolutely sure we never read past the end of the allocated string.
  !! We stop the iteration when we either find <startdate>/<enddate> or we 
  !! reach the character 11/9 elements from the end of the string (11/9 are 
  !! the lengths of the substrings we're comparing against).

  !! Let's keep a record of the original template, so make a copy to mutate
  !CurrentFile = FileTemplate

  !! Iterate through the string for <startdate>
  !FindStart: DO CharIndx = 1, LEN(CurrentFile) - 11
    !! Check against <startdate>, if found set the index
    !IF (CurrentFile(CharIndx:CharIndx+10) == StartDate) THEN
      !IndxStartDate = CharIndx
      !EXIT FindStart
    !END IF
  !END DO FindStart

  !! Iterate through for <enddate>
  !FindEnd: DO CharIndx = 1, LEN(CurrentFile) - 9
    !! Check against <enddate>, if found set the index
    !IF (CurrentFile(CharIndx:CharIndx+8) == EndDate) THEN
      !IndxEndDate = CharIndx
      !EXIT FindEnd
    !END IF
  !END DO FindEnd

  !! Check that both occurrences were found- this triggers if either remain 0
  !IF ((IndxStartDate == 0) .OR. (IndxEndDate == 0)) THEN
    !write(*,*) "File for "//CurrentFile//" does not match the template."
    !CALL EXIT(5)
  !END IF

  !! Now we go and replace the <startdate> and <enddate> strings with "*" so
  !! that we can pass the filename to the unix ls command. We also want to
  !! shift the remaining section of the string to the right of <startdate>/
  !! <enddate> 10/8 characters to the left to fill in the new gap.

  !! Do <enddate> first so we don't mess up the indexing
  !CurrentFile(IndxEndDate:IndxEndDate) = "*"
  !CurrentFile(IndxEndDate+1:) = CurrentFile(IndxEndDate+9:)

  !CurrentFile(IndxStartDate:IndxStartDate) = "*"
  !CurrentFile(IndxStartDate+1:) = CurrentFile(IndxStartDate+11:)

  !! Now we've processed the supplied string, pass it to the ls unix command.
  !! We need to write the output to a temporary file and read it back in to
  !! get the date ranges for each file.
  !CALL execute_command_line("ls -1 "//TRIM(CurrentFile)//" >&
     !__FileNameWithNoClashes__.txt", EXITSTAT = ExStat, CMDSTAT = CStat)

  !! Now we want to read this temporary file back in and inspect the contents
  !! to determine the date ranges for each file. We then write back to a
  !! reserved name file, that now contains the original filenames as well as
  !! the year ranges for each file as described in the opening preamble to
  !! this subroutine.

  !OPEN(InputUnit, FILE = "__FileNameWithNoClashes__.txt", IOSTAT = ios)

  !! Now iterate through this temporary file and check how many files are in the
  !! dataset, so we can allocate memory in the SpatioTemporalDataset
  !FileCounter = 0
  !CountFiles: DO
    !! Read in the line and check the status. We write the line to the variable,
    !! but it's not actually necessary at the moment.
    !READ(InputUnit, '(A)', IOSTAT = ios) CurrentFile

    !IF (ios == -1) THEN
      !! Read completed successfully
      !EXIT CountFiles
    !ELSEIF (ios /= 0) THEN
      !! Read failed for some other reason
      !CALL EXIT(5)
    !END If
    !! Otherwise, we read a line, so increase the counter
    !FileCounter = FileCounter + 1
  !END DO CountFiles

  !! Now allocate memory for the filenames and periods
  !ALLOCATE(Dataset%Filenames(FileCounter), Dataset%StartYear(FileCounter),&
          !Dataset%EndYear(FileCounter))

  !! Now rewind the file and read it again
  !REWIND(InputUnit)

  !BuildDataset: DO FileIndx = 1, FileCounter
    !! Read in the line, this time we do need to store the line
    !READ(InputUnit, '(A)', IOSTAT = ios) CurrentFile

    !! Place the filename in the dataset structure
    !Dataset%Filenames(FileIndx) = CurrentFile

    !! Reading the start year is easy- just the first 4 characters starting
    !! from the recorded IndxStartDate.
    !READ(CurrentFile(IndxStartDate:IndxStartDate+3), *)&
      !Dataset%StartYear(FileIndx)

    !! The end year is not as simple as taking the 4 characters at
    !! IndxEndDate, because the size of the <startdate> string is 3
    !! characters longer than the YYYYMMDD date specifier. Therefore, shift
    !! the index back 3 units and then take the next 4 characters
    !READ(CurrentFile(IndxEndDate-3:IndxEndDate), *) Dataset%EndYear(FileIndx)

  !END DO BuildDataset

  !! Close the file
  !CLOSE(InputUnit)
  
  !! Finished the work (we assign the VarNames later). Now delete the temporary
  !! file we used to store the `ls` output.
  !CALL execute_command_line("rm __FileNameWithNoClashes__.txt")
!END SUBROUTINE prepare_spatiotemporal_dataset

!!------------------------------------------------------------------------------!

!SUBROUTINE open_at_first_file(Dataset)
  !! Open the spatiotemporal dataset at the first file as an initialisation
  !TYPE(SPATIO_TEMPORAL_DATASET), INTENT(INOUT) :: Dataset

  !INTEGER :: ok

  !! Set the current FileIndx to 1, and open that NCdataset
  !Dataset%CurrentFileIndx = 1

  !ok = NF90_OPEN(Dataset%FileNames(Dataset%CurrentFileIndx), NF90_NOWRITE,&
    !Dataset%CurrentFileID)
  !IF (ok /= NF90_NOERR) THEN
    !CALL handle_err(ok, "Error opening "//Dataset%FileNames(Dataset%CurrentFileIndx))  
  !END IF

  !CALL find_variable_ID(Dataset)
!END SUBROUTINE open_at_first_file

!!------------------------------------------------------------------------------!

!SUBROUTINE find_variable_ID(Dataset)
  !! Find the desired variable name in the dataset
  !TYPE(SPATIO_TEMPORAL_DATASET), INTENT(INOUT) :: Dataset

  !INTEGER :: ok, VarNameIter

  !FindVar: DO VarNameIter = 1, SIZE(Dataset%VarNames)
    !ok = NF90_INQ_VARID(Dataset%CurrentFileID, Dataset%VarNames(VarNameIter),&
      !Dataset%CurrentVarID)
    !IF (ok == NF90_NOERR) THEN
      !EXIT FindVar
    !END IF
  !END DO FindVar

!END SUBROUTINE find_variable_ID

!!------------------------------------------------------------------------------!

SUBROUTINE prepare_temporal_dataset(FileName, TargetArray)
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
    READ(FileID, FMT = '(I4, F)', IOSTAT = ios) Time, TargetArray(iter)
  END DO ReadValues
  CLOSE(FileID)
END SUBROUTINE prepare_temporal_dataset

!------------------------------------------------------------------------------!

SUBROUTINE get_cru_co2(CRU, CO2Air)
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
  IF (CRU%CO2Method == "yearly") THEN
    ! We use the same year as the simulation year
    CO2Year = CRU%cYear
    CO2Air(:) = CRU%CO2Vals(CO2Year)
  ELSE
    ! The user has specified the year
    READ(CRU%CO2Method, '(I)') CO2Year
    CO2Air(:) = CRU%CO2Vals(CO2Year)
  END IF
END SUBROUTINE get_cru_co2

!------------------------------------------------------------------------------!

SUBROUTINE get_cru_ndep(CRU)
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
  IF (CRU%NDepMethod == "spatial") THEN
    ! We're using the simulation year for the deposition data
    TimeIndex = CRU%cYear - 1849
  ELSE
    ! We're using a specific year of data
    READ(CRU%NDepMethod, '(I4)') TimeIndex
    TimeIndex = TimeIndex - 1849
  END IF

  ! And finally grab the data. Since we're using landmask, we have to extract
  ! the full set of data first and then read the unmasked points into the NDep
  ! values.
  ALLOCATE(TmpArray(CRU%xDimSize, CRU%yDimSize))

  ApplyLandmask: DO GridCell = 1, CRU%mLand
    ok = NF90_GET_VAR(CRU%NDepFID, CRU%NDepVID, CRU%NDepVals(GridCell),&
    START = (/land_x(GridCell), land_y(GridCell), TimeIndex/))
    IF (ok /= NF90_NOERR) THEN
      CALL handle_err(ok, "Reading from NDep")
    END IF
  END DO ApplyLandMask

END SUBROUTINE get_cru_ndep

!!------------------------------------------------------------------------------!

!SUBROUTINE cru_read_metvals(STD, MetData, Year, DayOfYear,&
    !LeapYears)
  !! Get the data from a day
  !TYPE(SPATIO_TEMPORAL_DATASET), INTENT(INOUT)  :: STD
  !TYPE(CRU_MET_TYPE), INTENT(INOUT)             :: MetData
  !INTEGER, INTENT(IN)                           :: DayOfYear, Year
  !LOGICAL, INTENT(IN)                           :: LeapYears

  !! We'll need to compute the record index to grab
  !INTEGER     :: YearIndex, TimeIndex

  !! Iterators
  !INTEGER     :: FileIndx, VarNameIndx, YearIter, LandCell

  !! Status checker
  !INTEGER  :: ok

  !CHARACTER(LEN=256) :: CheckFilename
  !LOGICAL :: IsOpen

  !TimeIndex = DayOfYear
  !YearIndex = Year

  !! We've already opened a file, check whether the current year is open
  !IF ((Year >= STD%StartYear(STD%CurrentFileIndx)) .AND.&
    !(Year <= STD%EndYear(STD%CurrentFileIndx))) THEN
    !! In this instance, we don't need to do anything
    !CONTINUE
  !ELSE
    !ok = NF90_CLOSE(STD%CurrentFileID)

    !! First, check if we're outside the range of our data
    !IF (YEAR < STD%StartYear(1)) THEN
      !! Before the first year
      !STD%CurrentFileIndx = 1
      !YearIndex = STD%StartYear(1)
      !TimeIndex = 1
    !ELSE IF (YEAR > STD%EndYear(SIZE(STD%EndYear))) THEN
      !! After the last year
      !STD%CurrentFileIndx = SIZE(STD%EndYear)
      !YearIndex = STD%EndYear(STD%CurrentFileIndx)
      !IF ((LeapYears) .AND. (is_leapyear(YearIndex))) THEN
        !TimeIndex = 366
      !ELSE
        !TimeIndex = 365
      !END IF
    !ELSE
      !! Normal operation, we're in the year of our data
      !FindFile: DO FileIndx = 1, SIZE(STD%FileNames)
        !IF ((Year >= STD%StartYear(FileIndx)) .AND. &
          !(Year <= STD%EndYear(FileIndx))) THEN
          !STD%CurrentFileIndx = FileIndx
          !EXIT FindFile
        !END IF
      !END DO FindFile
    !END IF

    !! Open the new file
    !ok = NF90_OPEN(STD%FileNames(STD%CurrentFileIndx), NF90_NOWRITE,&
      !STD%CurrentFileID)

    !CALL find_variable_id(STD)
  !END IF

  !! Now read the desired time step

  !! Due to leapyears, we can't just do add 365 * number of years from startyear.
  !IF (LeapYears) THEN
    !CountDays: DO YearIter = STD%StartYear(STD%CurrentFileIndx), YearIndex
      !IF (is_leapyear(YearIndex)) THEN
        !TimeIndex = TimeIndex + 366
      !ELSE
        !TimeIndex = TimeIndex + 365
      !END IF
    !END DO CountDays
  !ELSE
    !TimeIndex = TimeIndex + 365 * (YearIndex - STD%StartYear(STD%CurrentFileIndx))
  !END IF
  !! Now we have the index, we can grab the data

  !! Read from the global array to the masked array
  !ApplyMask: DO LandCell = 1, mLand
    !ok = NF90_GET_VAR(STD%CurrentFileID,&
      !STD%CurrentVarID, MetData%MetVals(LandCell), START = (/land_x(LandCell), land_y(LandCell), TimeIndex/))
      !IF (ok /= NF90_NOERR) THEN
        !CALL handle_err(ok, "Failed reading "//STD%FileNames(STD%CurrentFileIndx))
      !END IF
  !END DO ApplyMask

!END SUBROUTINE cru_read_metvals

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
  ! Close the open NetCDF files
  TYPE(CRU_TYPE), INTENT(INOUT) :: CRU

  ! Iterator
  INTEGER  :: iter, ok

  CloseFiles: DO iter = 1, SIZE(CRU%MetDatasets)
    ok = NF90_CLOSE(CRU%MetDatasets(iter)%CurrentFileID)
  END DO CloseFiles
END SUBROUTINE cru_close

END MODULE CABLE_CRU
