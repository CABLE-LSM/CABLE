! Author: Lachlan Whyborn
! Last Modified: Mon 08 Apr 2024 01:36:00 PM AEST

MODULE CRU

  IMPLICIT NONE

TYPE SPATIO_TEMPORAL_DATASET
  ! This data type is used to store effectively a reference to spatio-temporal
  ! datasets which span multiple files, which is common for external forcings.
  ! The idea is to store the list of files in the dataset, and their start and
  ! end years respectively. We also need to know the internal variable names in
  ! the netCDF dataset so we can call it correctly. 
  ! This can possibly permit series with different sampling intervals as well.
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: FileNames
  INTEGER, DIMENSION(:), ALLOCATABLE            :: StartYear, EndYear
  CHARACTER(LEN=16), DIMENSION(:), ALLOCATABLE  :: VarNames
  ! As well as the metadata about the dataset contents, we want some data that
  ! is updated to inform which file is currently open and it's period. If open a
  ! file relevant to a particular year, and move onto the next year that's still
  ! contained in the current file's period, we don't want to go through open
  ! another file, we just want to continue using the current file.
  ! To do this, we need to know the current FileID, and which index it is in
  ! our list of files so we can access it's period.
  INTEGER CurrentFileID, CurrentFileIndx

END TYPE SPATIO_TEMPORAL_DATASET

TYPE CRU_TYPE
  INTEGER   :: mLand              ! Number of land cells in the run
  INTEGER   :: nMet               ! Number of Meteorological variables
  INTEGER   :: xDimSize, yDimSize ! Size of longitude and latitude dimensions
  INTEGER   :: cYear              ! CABLE main year
  INTEGER   :: DtSecs             ! Size of the timestep in seconds
  INTEGER   :: MetRecyc           ! Period of the met forcing recycling
  
  INTEGER, DIMENSION(10)    :: FileID, VarId    ! File and variable IDs for
  INTEGER                   :: NDepFId, NDepVId ! reading netCDF

  REAL, DIMENSION(:), ALLOCATABLE ::  avg_lwdn  ! Average longwave down rad,
                                                ! for Swinbank method

  REAL, DIMENSION(:), ALLOCATABLE ::  CO2Vals   ! Atmospheric carbon values
  REAL, DIMENSION(:), ALLOCATABLE ::  NDepVals  ! Nitrogen deposition values

  LOGICAL   :: DirectRead         ! Whether to read in entire met array
  LOGICAL   :: ReadDiffFrac       ! Read diff fraction or calculate it

  LOGICAL, DIMENSION(:,:), ALLOCATABLE  :: LandMask ! The logical landmask

  CHARACTER(LEN=16)   :: CarbonMethod   ! Method for choosing atmospheric CO2
  CHARACTER(LEN=16)   :: NDepMethod     ! Method for choosing N deposition

  CHARACTER(LEN=16), DIMENSION(10)    :: VarNames

  TYPE(CRU_MET_TYPE), DIMENSION(12)   :: Met

  TYPE(SPATIO_TEMPORAL_DATASET), DIMENSION(10)  :: MetDatasets
  TYPE(SPATIO_TEMPORAL_DATASET)                 :: NDepDataset
END TYPE CRU_TYPE

SUBROUTINE CRU_INIT(CRU)

  ! A rewrite of the CRU_INIT subroutine in cable_cru_TRENDY.f90

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
  ! The main goal of this routine is to mutate the CRU data structure to prep
  ! it for the experiment.
  TYPE(CRU_TYPE), INTENT(OUT)    :: CRU

  ! Start with the things we want from the namelist. The namelist must set
  ! the filenames to read from, the method of choosing atmospheric carbon and
  ! nitrogen deposition, and the timestep.
  ! We want one filename for each variable, stored in a list predefined index.
  CHARACTER(LEN=256), dimension(12) :: InputFiles
  CHARACTER(LEN=16)                 :: CarbonMethod, NDepMethod
  REAL                              :: DtSecs

  ! To determine the data we want from the Met data, we need to look at the
  ! landmask
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: LandMask
 
  ! Read the namelist and pull out the input files and configuration variables
  CALL read_MET_namelist_cbl(InputFiles, CRU)

  ! Read the landmask and allocate appropriate memory for the array variables
  CALL read_landmask(InputFiles(13), CRU)

  ! We know that the first 9 Met variables (indices 1-9) are always going to be
  ! time series, so always build their keys.
  BuildKeys: DO VarIndx = 1, 9
    prepare_spatialtemporal_dataset(InputFiles(VarIndx), CRU%VarNames(VarIndx),&
          CRU%MetDatasets(VarIndx))
  END DO BuildKeys

  ! Our building of NDep, Carbon and FDiff keys depend on our choices of method.
  IF (CRU%ReadDiffFrac) THEN
    prepare_spatiotemporal_dataset(InputFiles(10), VarNames(10),&
          CRU%MetDatasets(VarIndx))
  END IF

  IF (CRU%CarbonMethod == "Yearly") THEN
    prepare_temporal_dataset(InputFiles(11), CRU%CO2Vals)
  ELSEIF (CRU%CarbonMethod == "Spatial") THEN
    prepare_spatiotemporal_dataset(InputFiles(11), "CO2", CRU%CO2Vals)
  END IF

  IF (CRU%NDepMethod == "Yearly") THEN
    prepare_temporal_dataset(InputFiles(12), CRU%NDepVals)
  ELSEIF (CRU%NDepMethod == "Spatial") THEN
    prepare_spatiatemporal_dataset(InputFiles(12), "N_deposition",&
          CRU%NDepDataset)
  END IF
END SUBROUTINE CRU_INIT

SUBROUTINE read_MET_namelist_cbl(InputFiles, CRU)
  ! # read_MET_namelist_cbl(InputFiles, CarbonMethod, NDepMethod, DtSecs)
  ! Read the supplied cru.nml namelist and extract the desired input file names,
  ! method for reading the Carbon and Nitrogen deposition and the MET forcing
  ! timestep.
  ! ## Namelist inputs
  ! The entries expected in the namelist are:
  ! ```
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
  ! CarbonMethod = "yearly", "spatial", <year>
  ! NDepMethod = "yearly", "spatial", <year>
  ! ReadDiffFrac = T/F
  ! DtHrs = <real>
  ! ```
  ! The filename template is "path/to/file/prefix<startdate>_<enddate>suffix".
  ! For example, there may be a series of files describing the rain data at
  ! ```"data/rain/cru_rain_data_18500101_18991231_av.nc"```,
  ! ```"data/rain/cru_rain_data_19000101_19491231_av.nc"``` and so on.
  ! The entry for rainFile would be 
  ! ```"data/rain/cru_rain_data_<startdate>_<enddate>_av.nc".
  !
  ! The CarbonMethod and NDepMethod are "yearly", with globally constant values
  ! for each year, "spatial" with spatially varying values for each year, or
  ! a specific year which fixes the value at that specific year.

  ! Internally, we're going to pick up the individually separated inputs and
  ! pack them into a more convenient data structure.
  CHARACTER(LEN=256), DIMENSION(13), INTENT(OUT) :: InputFiles

  TYPE(CRU_TYPE), INTENT(INOUT) :: CRU    ! The master CRU structure

  ! Now we assign variables to temporarily store the files in the namelist
  ! In an ideal world, we'd have the users write directory to the InputFiles
  ! string array, but we can't expect everyone to know the order the files are
  ! expected in. So we first read them from a recognisable name, then pass them
  ! to the array.
  CHARACTER(LEN=256)  :: rainFile, lwdnFile, swdnFile, presFile, qairFile,&
                         TmaxFile, TminFile, uwindFile, vwindFile, fdiffFile,&
                         carbonFile, NDepFile, landmaskFile
  CHARACTER(LEN=16)   :: CarbonMethod, NDepMethod
  REAL                :: DtHrs
  INTEGER             :: MetRecyc
  LOGICAL             :: ReadDiffFrac, RecycleMet

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
  carbonFile = "None"
  NDepFile = "None"
  landmaskFile = "None"

  ! Defaults for the other inputs
  CarbonMethod = "Yearly"
  NDepMethod = "Yearly"
  DtHrs = 3
  MetRecyc = 20
  ReadDiffFrac = .TRUE.
  RecycleMet = .TRUE.

  ! Set up and read the namelist
  NAMELIST /metnml/ rainFile, lwdnFile, swdnFile, presFile, qairFile, TmaxFile,&
                    TminFile, uwindFile, vwindFile, fdiffFile, carbonFile,&
                    NDepFile, landmaskFile, CarbonMethod, NDepMethod, DtHrs

  ! Get a temporary unique ID and use it to read the namelist
  CALL get_unit(nmlUnit)
  OPEN(nmlUnit, file = "met.nml", status = 'old', action = 'read')
  READ(nmlUnit, nml = metnml)
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
  InputFiles(9) = vinwdFile
  InputFiles(10) = fdiffFile
  InputFiles(11) = CarbonFile
  InputFiles(12) = NDepFile
  InputFiles(13) = landmaskFile

  ! Bind the remaining config variables to the CRU structure
  CRU%CarbonMethod = CarbonMethod
  CRU%NDepMethod = NDepMethod
  ! Convert the hourly timestep to seconds
  CRU%DtSecs = DtHrs * 3600
  CRU%MetRecyc = 20
  CRU%ReadDiffFrac = ReadDiffFrac
  CRU%RecycleMet = RecycleMet
END SUBROUTINE read_MET_namelist_cbl

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
  INTEGER       :: xDimSize, yDimSize

  ! And then the actual values along the axes
  REAL, DIMENSION(:), ALLOCATABLE :: Latitudes, Longitudes

  ! Need a counter when iterating through the points in the land mask
  INTEGER       :: MaskCounter

  ! Now we attempt to read the netCDF file
  ErrStatus = nf90_open(TRIM(LandmaskFile), nf90_nowrite, FileID)

  ! Inquire about the latitude and longitude dimensions
  ErrStatus = nf90_inq_dimid(FileID, 'latitude', LatID)
  ErrStatus = nf90_inqure_dimension(FileID, LatID, len = yDimSize)

  ErrStatus = nf90_inq_dimid(FileID, 'longitude', LonID)
  ErrStatus = nf90_inqure_dimension(FileID, LonID, len = xDimSize)

  ! And read the dimensions into arrays
  ALLOCATE(Latitudes(yDimSize)
  ErrStatus = nf90_inq_varid(FileID, 'latitude', LatID)
  ErrStatus = nf90_get_var(FileID, LatID, Latitudes)

  ALLOCATE(Longitudes(xDimSize))
  ErrStatus = nf90_inq_varid(FileID, 'longitude', LonID)
  ErrStatus = nf90_get_var(FileID, LonID, Longitudes)

  ! Now read the mask, this mask is an integer (0/1) version
  ALLOCATE(mask(xDimSize, yDimSize))
  ErrStatus = nf90_inq_varid(FileID, 'land', LandID)
  ErrStatus = nf90_get_var(FileID, LandID, mask)

  ! Finally close the mask file
  ErrStatus = nf90_close(FileID)

  ! For some reason we also want a logical version
  ALLOCATE(CRU%LandMask(xDimSize, yDimSize))
  CRU%LandMask = to_logical(mask)

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
        Land_x = LonIndx
        Land_y = LatIndx
        Longitude(MaskCounter) = Longitudes(LonIndx)
        Latitude(MaskCounter) = Latitudes(LatIndx)
        MaskCounter = MaskCounter + 1
      END IF
    END DO LoopThoughLons
  END DO LoopThroughLats

  ! Allocate memory for the actual meteorological forcing data
  AllocateMemory: DO VarIndx = 1, SIZE(CRU%Met)
    ALLOCATE(CRU%Met(VarIndx)%MetVals(CRU%mLand))
  END DO AllocateMemory

  ! And do the same for nitrogen deposition
  ALLOCATE(CRU%NDepVals(CRU%mLand))

  ! Allocate global data- this practice must be purged in future.
  metGrid = "mask"
  mLand = CRU%mLand
  nMetPatches = 1

  ALLOCATE(Lat_all(xDimSize, yDimSize), Lon_all(xDimSize, yDimSize))
  FillLatitudes: DO LonIndx = 1, xDimSize
    Lat_all(LonIndx, :) = Latitudes
  END DO FillLatitudes

  FillLongitudes DO LatIndx = 1, yDimSize
    Lon_all(LatIndx, :) = Longitudes
  END DO FillLongitudes

  ! Some CABLE time units?
  shod = 0
  sdoy = 1
  smoy = 1
  syear = CRU%CYEAR

  ! Allocate memory for the average long-wave radiation down for some reason?
  ALLOCATE(CRU%AVG_LWDN(mLand))
END SUBROUTINE read_landmask

ELEMENTAL PURE FUNCTION to_logical(IntegerIn)
  ! Convert an integer (0/1) input to a logical
  ! We put the onus on the user to ensure the inputs are valid
  INTEGER, INTENT(IN) :: IntegerIn

  IF IntegerIn == 0 THEN
    to_logical = .FALSE.
  ELSE
    to_logical = .TRUE.
  END IF
END FUNCTION integer_to_logical

SUBROUTINE prepare_temporal_dataset(FileName, TargetArray)
  ! If we have a dataset that's only varying in time, then we should be able to
  ! just the whole thing into memory immediately.

  CHARACTER(LEN=256), INTENT(IN)    :: FileName
  REAL, DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: TargetArray

  ! A unique file ID
  INTEGER   :: FileID

  ! A placeholder for the line, that we inspect to check for comment lines
  CHARACTER(LEN=256)    :: LineInFile

  ! Placeholders for the year key in the line. We also want to store the start
  ! and end years so that we can set the custom array indexing for the data
  INTEGER    :: Time, StartTime, EndTime

  ! Need a dummy variable to store the values on the first read of the file
  REAL       :: DummyValue

  ! The status holder, and the file and header counters
  INTEGER :: ios, LineCounter, HeaderCounter

  ! We're going to assume that we only have a single file.
  ! We need to iterate through the file twice, first to inspect the number of
  ! We need to iterate through the file twice, first to inspect the number of
  ! entries in the file to allocate the correct amount of memory, second to
  ! actually write the data to the array.
  CALL get_unit(FileID)
  OPEN(FileID, FileName, STATUS = "old", ACTION = "read")

  LineCounter = 0
  DetermineFileSize: DO
    ! Read line by line, checking for a header line
    READ(FileID, '(A)', IOSTAT = ios) LineInFile
    
    IF (ios == -1) THEN
      ! Read completed successfully and hit EOF
      ! Assign Time to EndTime
      EndTime = Time

      EXIT DetermineSize
    END IF

    ! Check for comment characters- assume they're all at the top of the file
    IF (LineInFile(1) == '#') .OR. (LineInFile(1) == '!') THEN
      HeaderCounter = HeaderCounter + 1
    ELSE
      ! Otherwise, its a line of useful data
      LineCounter = LineCounter + 1
      
      ! To correctly index the array, we need to pull out the time interval
      ! from the line
      READ(LineInFile, '(I) (F)') Time, DummyValue
    END IF

    IF (LineCounter == 1) THEN
      ! Assign Time to StartTime
      StartTime = Time
    END IF

  END DO DetermineFileSize

  ! Rewind the cursor to the start of the file
  REWIND(FileID)

  ! Now allocate the appropriate memory
  ALLOCATE(TemporaryArray(StartTime:EndTime))

  ! First step over the header lines
  SkipHeader: DO iter = 1, HeaderCounter
    READ(FileID, '(A)', IOSTAT = ios) LineInFile
  END DO SkipHeader

  ! Now the useful contents of the file into the array
  ReadValues: DO iter = StartTime, EndTime
    READ(FileID, '(I4) (F)', IOSTAT = ios) Time, TargetArray(iter)
  END DO ReadValues
END SUBROUTINE prepare_temporal_dataset

SUBROUTINE prepare_spatiotemporal_dataset(FileTemplate, VariableName, Dataset)
  ! Prepare the met netCDF files for reading. We follow a process similar to
  ! that of JULES to handle time series data split over multiple files, see
  ! 6.31.2.2 at https://jules-lsm.github.io/latest/namelists/drive.nml.html.

  ! The user specifies a template matching the file names for each variable,
  ! containing the start and end date of the data in the file. For example,
  ! rain data files may be named "pre/precip_1850101_18991231.nc",
  ! "pre/precip_19000101_19491231.nc" and so on. Then the user would specify
  ! rainFile = "pre/precip_<startdate>_<enddate>.nc" in the namelist.

  ! We then use the ls command to find all files matching this pattern and 
  ! write the output to a text file. We then inspect the filenames to
  ! determine the time period for each file.

  !! Then add to them to a file so what we end up with is a file that contains:
  !! "<filename_1>.nc" startYear endYear
  !! "<filename_2>.nc" startYear endYear
  !! So for the example above, it would look like
  !! "pre/precip_18500101_18991231.nc" 1850 1899
  !! "pre/precip_19000101_19491231.nc" 1900 1949

  !! We can then read these filelists with a standardised reader.
  !! This approach should be able to be extended to any time series data, but
  !! we're just applying it to the Met data for now.

  ! Alternatively, we can create a data structure which behaves like this file.
  ! The data structure contains the list of files, their periods, as well as
  ! a reference to the current open file.
  ! Specify the intent of the arguments
  CHARACTER(len=256), INTENT(IN)  :: FileTemplate
  CHARACTER(len=8)  , INTENT(IN)  :: VariableName
  TYPE(SPATIO_TEMPORAL_DATASET), INTENT(OUT)  :: Dataset

  ! We need to have a definition of the substrings we're going to replace.
  CHARACTER(len=11) :: StartDate = "<startdate>"
  CHARACTER(len=9)  :: EndDate = "<enddate>"

  ! For clarity in the looping, and so we don't mutate the original inputs,
  ! we extract the filename from the array.
  CHARACTER(len=256):: CurrentFile

  !! We write the outputs to a specific reserved filename
  !CHARACTER(len=32) :: ReservedOutputFile

  ! We read the start and end years to strings before writing them to the key.
  CHARACTER(len=4)  :: StartYear, EndYear

  ! Initialise integers to store the status of the command.
  INTEGER           :: ExStat, CStat
  ! Initialise integers to store the status of the command.
  INTEGER           :: ExStat, CStat

  ! We're going to want to remember where "<startdate>" and "<enddate>" occur
  ! in the filename templates, so we can retrieve the time periods of each
  ! file later.
  INTEGER           :: IndxStartDate = 0, IndxEndDate = 0

  ! INTEGERs for the file IDs that read and write the file metadata.
  INTEGER           :: InputUnit!, OutputUnit

  ! And status integers for those files
  INTEGER           :: iosIn!, iosOut

  ! Need to compute the number of files in the dataset
  INTEGER           :: FileCounter

  ! Get a unique ID here.
  CALL get_unit(InputUnit)
  !CALL get_unit(OutputUnit)

  ! First thing we have to do is determine the location of the <startdate>
  ! and <enddate> strings in the supplied filenames. To do that, we compare
  ! substrings of the relevant length, with the substring moving along the
  ! length of the filename string.
  ! We iterate separately for <startdate> and <enddate> so that we make
  ! absolutely sure we never read past the end of the allocated string.
  ! We stop the iteration when we either find <startdate>/<enddate> or we 
  ! reach the character 11/9 elements from the end of the string (11/9 are 
  ! the lengths of the substrings we're comparing against).

  ! Build the reserved filename for the output
  ReservedOutputFile = VariableName//"Key.txt"

  ! Let's keep a record of the original template, so make a copy to mutate
  CurrentFile = FileTemplate

  ! Iterate through the string for <startdate>
  FindStart: DO CharIndx = 1, LEN(CurrentFile) - 11
    ! Check against <startdate>, if found set the index
    IF (CurrentFile(CharIndx:CharIndx+10) == StartDate) THEN
      IndxStartDate = CharIndx
      EXIT FindStart
    END IF
  END DO FindStart

  ! Iterate through for <enddate>
  FindEnd: DO CharIndx = 1, LEN(CurrentFile) - 9
    ! Check against <enddate>, if found set the index
    IF (CurrentFile(CharIndx:CharIndx+8) == EndDate) THEN
      IndxEndDate = CharIndx
      EXIT FindEnd
    END IF

  END DO FindEnd

  ! Check that both occurrences were found- this triggers if either remain 0
  IF ((.NOT. IndxStartDate) .OR. (.NOT. IndxEndDate)) THEN
    write(*,*) "Supplied file for "//CurrentFile//" is invalid."
    CALL EXIT(5)
  END IF

  ! Now we go and replace the <startdate> and <enddate> strings with "*" so
  ! that we can pass the filename to the unix ls command. We also want to
  ! shift the remaining section of the string to the right of <startdate>/
  ! <enddate> 10/8 characters to the left to fill in the new gap.

  ! Do <enddate> first so we don't mess up the indexing
  CurrentFile(IndxEndDate:IndxEndDate) = "*"
  CurrentFile(IndxEndDate+1:) = CurrentFile(IndxEndDate+9:)

  CurrentFile(IndxStartDate:IndxStartDate) = "*"
  CurrentFile(IndxStartDate+1:) = CurrentFile(IndxStartDate+11:)

  ! Now we've processed the supplied string, pass it to the ls unix command.
  ! We need to write the output to a temporary file and read it back in to
  ! get the date ranges for each file.
  CALL execute_command_line("ls -1 "//TRIM(CurrentFile)//" >&
     __tmpFileWithNoClashes__.txt", EXITSTAT = ExStat, CMDSTAT = CStat)

  ! Now we want to read this temporary file back in and inspect the contents
  ! to determine the date ranges for each file. We then write back to a
  ! reserved name file, that now contains the original filenames as well as
  ! the year ranges for each file as described in the opening preamble to
  ! this subroutine.

  OPEN(InputUnit, "__tmpFileWithNoClashes__.txt", ios = ios1)

  ! Now iterate through this temporary file and check how many files are in the
  ! dataset
  FileCounter = 0
  CountFiles: DO
    ! Read in the line and check the status
    READ(InputUnit, '(A)', IOSTAT = iosIn) CurrentFile

    IF (iosIn == -1) THEN
      ! Read completed successfully
      EXIT BuildKey
    ELSEIF (iosIn /= 0) THEN
      ! Read failed for some other reason
      EXIT(5)
    END If

    ! Otherwise, we read a line, so increase the counter
    FileCounter = FileCounter + 1
  END DO CountFiles

  ! Now allocate memory for the filenames and periods
  ALLOCATE(Dataset%Filenames(FileCounter), Dataset%StartYear(FileCounter),&
          Dataset%EndYear(FileCounter))

  ! Now rewind the file and read it again
  REWIND(InputUnit)

  BuildDataset: DO Indx = 1, FileCounter
    ! Read in the line
    READ(InputUnit, '(A)', IOSTAT = iosIn) CurrentFile

    ! Place the filename in the dataset structure
    Dataset%Filenames(Indx) = CurrentFile

    ! Reading the start year is easy- just the first 4 characters starting
    ! from the recorded IndxStartDate.
    READ(CurrentFile(IndxStartDate:IndxStartDate+3), *) Dataset%StartYear(Indx)

    ! The end year is not as simple as taking the 4 characters at
    ! IndxEndDate, because the size of the <startdate> string is 3
    ! characters longer than the YYYYMMDD date specifier. Therefore, shift
    ! the index back 3 units and then take the next 4 characters
    READ(CurrentFile(IndxEndDate-3:IndxEndDate), *) Dataset%EndYear(Indx)

  END DO BuildDataset
  !OPEN(OutputUnit, ReservedOutputFile, ios = ios2)

  !! Now iterate through the lines in the input file, inspect the string
  !! to extract the bracketing years. We then append the start and the end
  !! year to the line. We can reuse the CurrentFile string to store the read
  !! line from the temporary file

  !BuildKey: DO
    !! Read in the line and check the status
    !READ(InputUnit, '(A)', IOSTAT = iosIn) CurrentFile

    !IF (iosIn == -1) THEN
      !! Read completed successfully
      !EXIT BuildKey
    !ELSEIF (iosIn /= 0) THEN
      !! Read failed for some other reason
      !EXIT(5)
    !END If

    !! We assume the dates in the filename in the YYYYMMDD format.
    !! For now, we're assuming that the data files contain complete years,
    !! with partial year splits across the files, so all we need is the year.
    !! If in future we have more fine grained data, we can modify this and
    !! the key reader routines.

    !! Reading the start year is easy- just the first 4 characters starting
    !! from the recorded IndxStartDate.
    !StartYear = CurrentFile(IndxStartDate:IndxStartDate+3)

    !! The end year is not as simple as taking the 4 characters at
    !! IndxEndDate, because the size of the <startdate> string is 3
    !! characters longer than the YYYYMMDD date specifier. Therefore, shift
    !! the index back 3 units and then take the next 4 characters
    !EndYear = CurrentFile(IndxEndDate-3:IndxEndDate)

    !! Now we have the information required to build the Spatio_temporal_dataset
    !! structure
    !Dataset%
    !! Now we have everything we need to write back to the key file, in the
    !! format {filename} {StartYear} {EndYear}.

    !!WRITE(OutputUnit, '(A)', IOSTAT = iosOut) TRIM(CurrentFile)//" "//&
      !!StartYear//" "//EndYear

    !!! Any non-zero iostat is bad
    !!IF (iosOut =/ 0) THEN
      !!EXIT(5)
    !!END IF
  !END DO BuildKey

  ! Now we've build the key for the variable, we can remove the temporary
  ! we used to store the filenames
  CALL execute_command_line("rm __tmpFileWithNoClashes__.txt")
END SUBROUTINE build_dataset_key

SUBROUTINE get_cru_CO2(CRU, CurrentYear, CO2air)
  ! Get the atmospheric CO2

  TYPE(CRU_TYPE), INTENT(IN)    :: CRU
  INTEGER, INTENT(IN)           :: CurrentYear
  REAL, INTENT(OUT)             :: CO2air

  ! When we retrieve a specific year's CO2, we want to convert a string to int
  INTEGER                       :: CO2Year

  ! How do we actually get the correct CO2?
  ! First, check the CO2 method
  IF (CRU%CO2Method == "Yearly") THEN
    CO2air = CRU%CO2Vals(CurrentYear)
  ELSEIF (CRU%CO2Method == "Spatial") THEN
    CONTINUE
    ! This clause is not yet implemented.
  ELSE
    ! Get the CO2 from a specific year
    READ(CRU%CO2Method, '(I)') CO2Year
    CO2Air = CRU%CO2Vals(CO2Year)
  END IF
END SUBROUTINE get_cru_co2    
  
SUBROUTINE CRU_GET_SUBDIURNAL_MET(CRU, MET, CurYear, ktau, kend)

  ! Obtain one day of CRU-NCEP meteorology, subdiurnalise it using a weather
  ! and return the result to the CABLE driver.

  USE cable_def_types_mod,   ONLY: MET_TYPE, r_2
  USE cable_IO_vars_module,  ONLY: LANDPT, latitude
  USE cable_common_module,   ONLY: DOYSOD2YMDHMS
  USE cable_weathergenerator,ONLY: WEATHER_GENERATOR_TYPE, WGEN_INIT, &
       WGEN_DAILY_CONSTANTS, WGEN_SUBDIURNAL_MET
  USE mo_utils,              ONLY: eq

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
  REAL      :: CO2air                 ! CO2 concentration in ppm
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
  !CRU%ktau  = ktau     ! ktau is the current timestep in the year.

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
