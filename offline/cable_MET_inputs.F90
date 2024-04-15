! Author: Lachlan Whyborn
! Last Modified: Mon 15 Apr 2024 16:45:50

MODULE MET_INPUTS

IMPLICIT NONE

USE NetCDF
USE CRU,      ONLY: CRU_TYPE

TYPE SPATIO_TEMPORAL_DATASET
  ! When we actually want a piece of data, we pass this along with the year
  ! and day to a read routine. The routine checks the passed year against
  ! StartYear(CurrentFileIndx) and EndYear(CurrentFileIndx), and decides whether
  ! to either a) continue with the current file if it falls within that range,
  !        or b) move to a new file if it falls outside that range.
  ! By storing the ID associated with the currently open file, we only need to
  ! open a new IO stream a couple of times per simulation, rather than every
  ! time we start a new year or a new time step.

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
  INTEGER :: CurrentFileID, CurrentVarId, CurrentFileIndx = 0
END TYPE SPATIO_TEMPORAL_DATASET

CONTAINS

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
  ! TRENDY configuration, so keep it for consistency.

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
  CHARACTER(LEN=128)  :: CO2Method, NDepMethod
  INTEGER             :: MetRecyc
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
                    ReadDiffFrac, CO2Method, NDepMethod, MetRecyc, LeapYears

  ! Get a temporary unique ID and use it to read the namelist
  nmlUnit = 1
  OPEN(nmlUnit, file = "met_config.nml", status = 'old', action = 'read')
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
  CRU%DtSecs = DtHrs * 3600
  CRU%MetRecyc = MetRecyc
  CRU%ReadDiffFrac = ReadDiffFrac
END SUBROUTINE read_MET_namelist_cbl

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

  CALL get_unit(nml_unit)
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

SUBROUTINE prepare_spatiotemporal_dataset(FileTemplate, Dataset)
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

  ! We then create a data structure which effectively holds metadata about
  ! the possibly multi-file dataset.
  ! The data structure contains the list of files, their periods, as well as
  ! a reference to the current open file.
  ! Specify the intent of the arguments
  CHARACTER(len=256), INTENT(IN)  :: FileTemplate
  TYPE(SPATIO_TEMPORAL_DATASET), INTENT(OUT)  :: Dataset

  ! We need to have a definition of the substrings we're going to replace.
  CHARACTER(len=11) :: StartDate = "<startdate>"
  CHARACTER(len=9)  :: EndDate = "<enddate>"

  ! For clarity in the looping, and so we don't mutate the original inputs,
  ! we extract the filename from the array.
  CHARACTER(len=256):: CurrentFile

  ! We write the outputs from the `ls` command to a specific reserved filename
  CHARACTER(len=32) :: ReservedOutputFile = "__FileNameWithNoClashes__.txt"

  ! We read the start and end years to strings before writing them to the key.
  CHARACTER(len=4)  :: StartYear, EndYear

  ! Initialise integers to store the status of the command.
  INTEGER           :: ExStat, CStat

  ! We're going to want to remember where "<startdate>" and "<enddate>" occur
  ! in the filename templates, so we can retrieve the time periods of each
  ! file later.
  INTEGER           :: IndxStartDate = 0, IndxEndDate = 0

  ! INTEGERs for the file ID that reads the output from the 'ls' command.
  INTEGER           :: InputUnit

  ! And status integers for the IO
  INTEGER           :: ios

  ! Need to compute the number of files in the dataset
  INTEGER           :: FileCounter

  ! Iterators
  INTEGER   :: CharIndx, FileIndx

  ! Get a unique ID here.
  CALL get_unit(InputUnit)

  ! First thing we have to do is determine the location of the <startdate>
  ! and <enddate> strings in the supplied filenames. To do that, we compare
  ! substrings of the relevant length, with the substring moving along the
  ! length of the filename string.
  ! We iterate separately for <startdate> and <enddate> so that we make
  ! absolutely sure we never read past the end of the allocated string.
  ! We stop the iteration when we either find <startdate>/<enddate> or we 
  ! reach the character 11/9 elements from the end of the string (11/9 are 
  ! the lengths of the substrings we're comparing against).

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
  IF ((IndxStartDate == 0) .OR. (IndxEndDate == 0)) THEN
    write(*,*) "File for "//CurrentFile//" does not match the template."
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
     __FileNameWithNoClashes__.txt", EXITSTAT = ExStat, CMDSTAT = CStat)

  ! Now we want to read this temporary file back in and inspect the contents
  ! to determine the date ranges for each file. We then write back to a
  ! reserved name file, that now contains the original filenames as well as
  ! the year ranges for each file as described in the opening preamble to
  ! this subroutine.

  OPEN(InputUnit, FILE = "__FileNameWithNoClashes__.txt", IOSTAT = ios)

  ! Now iterate through this temporary file and check how many files are in the
  ! dataset, so we can allocate memory in the SpatioTemporalDataset
  FileCounter = 0
  CountFiles: DO
    ! Read in the line and check the status. We write the line to the variable,
    ! but it's not actually necessary at the moment.
    READ(InputUnit, '(A)', IOSTAT = ios) CurrentFile

    IF (ios == -1) THEN
      ! Read completed successfully
      EXIT CountFiles
    ELSEIF (ios /= 0) THEN
      ! Read failed for some other reason
      CALL EXIT(5)
    END If

    ! Otherwise, we read a line, so increase the counter
    FileCounter = FileCounter + 1
  END DO CountFiles

  ! Now allocate memory for the filenames and periods
  ALLOCATE(Dataset%Filenames(FileCounter), Dataset%StartYear(FileCounter),&
          Dataset%EndYear(FileCounter))

  ! Now rewind the file and read it again
  REWIND(InputUnit)

  BuildDataset: DO FileIndx = 1, FileCounter
    ! Read in the line, this time we do need to store the line
    READ(InputUnit, '(A)', IOSTAT = ios) CurrentFile

    ! Place the filename in the dataset structure
    Dataset%Filenames(FileIndx) = CurrentFile

    ! Reading the start year is easy- just the first 4 characters starting
    ! from the recorded IndxStartDate.
    READ(CurrentFile(IndxStartDate:IndxStartDate+3), *)&
      Dataset%StartYear(FileIndx)

    ! The end year is not as simple as taking the 4 characters at
    ! IndxEndDate, because the size of the <startdate> string is 3
    ! characters longer than the YYYYMMDD date specifier. Therefore, shift
    ! the index back 3 units and then take the next 4 characters
    READ(CurrentFile(IndxEndDate-3:IndxEndDate), *) Dataset%EndYear(FileIndx)

  END DO BuildDataset
  
  ! Finished the work (we assign the VarNames later). Now delete the temporary
  ! file we used to store the `ls` output.
  CALL execute_command_line("rm __tmpFileWithNoClashes__.txt")
END SUBROUTINE prepare_spatiotemporal_dataset

SUBROUTINE cru_read_metvals(STD, MetData, Year, DayOfYear,&
    LeapYears)
  ! Get the data from a day
  TYPE(SPATIO_TEMPORAL_DATASET), INTENT(INOUT)  :: STD
  TYPE(CRU_MET_TYPE), INTENT(INOUT)               :: MetData
  INTEGER, INTENT(INOUT)                        :: DayOfYear
  LOGICAL, INTENT(IN)                           :: LeapYears
  INTEGER, INTENT(INOUT)                        :: Year

  ! Need a temporary array to read in the data from the netcdf file
  REAL, DIMENSION(:, :), ALLOCATABLE            :: TmpArray
  INTEGER   :: LatSize, LonSize
  INTEGER   :: LatID, LonID

  ! Iterators
  INTEGER   :: FileIndx, VarNameIndx, YearIndx, ok

  ! We'll need to compute the record index to grab
  INTEGER     :: TimeIndex

  ! Check 2 conditions- if the file index is 0, it means we 
  ! haven't selected a file for opening yet, or the year is outside the range of
  ! the current file. In either case, we need to choose a new file.
  WRITE(*,*) "Start looking for the file."
  IF ((STD%CurrentFileIndx .EQ. 0) .OR.&
    ((YEAR < STD%StartYear(STD%CurrentFileIndx)) .AND.&
    (YEAR > STD%EndYear(STD%CurrentFileIndx)))) THEN
    ! In either instance, find the relevant file
    ! First, we might be in an instance where we're in no file's time period.
    WRITE(*,*) "Our current index is wrong."
    IF (YEAR < STD%StartYear(1)) THEN
      WRITE(*,*) "The year is before the era of the data."
      ! We're before the era of our data, just get the first day in the dataset
      STD%CurrentFileIndx = 1
      WRITE(*,*) STD%StartYear
      YEAR = STD%StartYear(1)
      WRITE(*,*) "Got year"
      DayOfYear = 1
    ELSEIF (YEAR > STD%EndYear(SIZE(STD%EndYear))) THEN
      WRITE(*,*) "The year is after the era of the data."
      ! We're after the era of our data, get the last day of the last year
      STD%CurrentFileIndx = SIZE(STD%FileNames)
      YEAR = STD%EndYear(SIZE(STD%EndYear))
      ! Need special handling if leapyears
      IF (LEAPYEARS) THEN
        IF (is_leapyear(YEAR)) THEN
          DayOfYear = 366
        ELSE
          DayOfYear = 365
        END IF
      ELSE
        DayOfYear = 365
      END IF
    ELSE
      ! Standard operation, we're in the era of our data
      FindFile: DO FileIndx = 1, SIZE(STD%FileNames)
        IF ((YEAR >= STD%StartYear(FileIndx)) .AND. &
          (YEAR <= STD%EndYear(FileIndx))) THEN
          STD.CurrentFileIndx = FileIndx
          EXIT FindFile
        END IF
      END DO FindFile
    END IF
  END IF
  WRITE(*,*) "Successfully found the file."

    !! We've selected the file, obtain a unique file ID and open the input stream
    !CALL get_unit(STD%CurrentFileID)
    !ok = NF90_OPEN(STD%FileNames(STD%CurrentFileIndx), NF90_NOWRITE,&
      !STD%CurrentFileID)

    !! Find the variable and get it's ID
    !FindVar: DO VarNameIndx = 1, SIZE(STD%VarNames)
      !ok = NF90_INQ_VARID(STD%CurrentFileID,&
        !STD%VarNames(VarNameIndx),&
        !STD%CurrentVarID)
      !IF (ok == NF90_NOERR) THEN
        !EXIT FindVar
      !END IF
    !END DO FindVar
  !END IF

  ! Now read the desired time step
  TimeIndex = DayOfYear

  ! Due to leapyears, we can't just do add 365 * number of years from startyear.
  IF (LeapYears) THEN
    CountDays: DO YearIndx = STD%StartYear(&
      STD%CurrentFileIndx), (Year - 1)
      IF (is_leapyear(YEAR)) THEN
        TimeIndex = TimeIndex + 366
      ELSE
        TimeIndex = TimeIndex + 365
      END IF
    END DO CountDays
  ELSE
    TimeIndex = TimeIndex + 365 * (Year - STD%StartYear(&
    STD%CurrentFileIndx) - 1)
  END IF

  WRITE(*,*) STD%FileNames(STD%CurrentFileIndx), TimeIndex
  !! Need to find out the dimensions of the data
  !ok = NF90_INQ_DIMID(STD%CurrentFileID, "lat", LatID)
  !ok = NF90_INQUIRE_DIMENSION(STD%CurrentFileID, LatID, LEN = LatSize)
  !ok = NF90_INQ_DIMID(STD%CurrentFileID, "lon", LonID)
  !ok = NF90_INQUIRE_DIMENSION(STD%CurrentFileID, LonID, LEN = LonSize)

  !! Now we have the index, we can grab the data
  !ALLOCATE(TmpArray(LatSize, LonSize))

  !ok = NF90_GET_VAR(STD%CurrentFileID,&
    !STD%CurrentVarID, TmpArray, START=(/1, 1, TimeIndex/))

END SUBROUTINE cru_read_metvals
