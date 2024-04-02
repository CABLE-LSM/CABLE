! Author: Lachlan Whyborn
! Last Modified: Wed 03 Apr 2024 08:56:33 AM AEDT

module CRU

  type CRU_FILES
    ! Store the templated strings for each of the meteorology files. The intent
    ! is for the filenames to follow a template that allows simple selection
    ! of the correct file for a given year. The pattern is as follows:
    ! some/long/path/pre_filename_<startdate>_<enddate>_post_filename.nc
    ! This way we can easily extract the time period for a given file.
    ! Each file relates to a local parameter names in the meteorological 
    ! forcing.
    character(len=200)  :: rainFile, lwdnFile, swdnFile, presFile, qairFile, &
                        :: tmaxFile, tminFile, uwindFile, fdiffFile, co2File, &
                        :: ndepFile

  end type CRU_FILES

  type CRU_CONFIG
    ! Contains all the configuration data for the cru met forcing

  subroutine CRU_INIT(CRU)

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
#endif

    implicit none

    ! Type that will contain the filenames for the CRU variables.
    type(CRU_FILES), intent(inout) :: CRUFiles
    
    ! Handle the constructing of the filenames and periods for each file
    call PREPARE_MET_FILES(CRUFILES)






    integer :: ErrStatus  ! Error status returned by nc routines (zero=ok, non-zero=error)
    integer :: nmlunit    ! Unit number for reading namelist file
    integer :: FID        ! NetCDF id for the landmask file
    integer :: latID, lonID  ! NetCDF ids for dimensions in the landmask file
    integer :: landID     ! NetCDF id for the landmask variable in the landmask file
    integer :: landcnt    ! Manually incremented counter for the number of land cells
    integer :: xcol, yrow ! Column and row position in the data file grids
    integer :: imetvar    ! loop counter through met variables

    ! We have two separate namelists within the cru namelist file, one to set
    ! all the filenames, and one to set the configuration options.
    ! Set temporary names to read from the filename namelist. 200 characters
    ! should be sufficient for any file.
    character(len=200)  :: MetPath, rainFile, lwdnFile, swdnFile, presFile, &
                           qairFile, tmaxFile, tminFile, uwindFile, fdiffFile, &
                           co2File, ndepFile

    ! the contents of the namelist will be different
    namelist /crufiles/ MetPath, rainFile, lwdnFile, swdnFile, presFile, &
              qairFile, tmaxFile, tminFile, uwindFile, vwindFile, fdiffFile, &
              co2File, ndepFile

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
  
    namelist /cruconfig/ 

  subroutine PREPARE_MET_FILES(METFILES, METVARIABLES)

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
    ! determine the time period for each file and add to them to the file so
    ! what we end up with is a file that contains:
    ! "<filename_1>.nc" startYear endYear
    ! "<filename_2>.nc" startYear endYear
    ! So for the example above, it would look like
    ! "pre/precip_18500101_18991231.nc" 1850 1899
    ! "pre/precip_19000101_19491231.nc" 1900 1949

    ! We can then read these filelists with a standardised reader.
    ! This approach should be able to be extended to any time series data, but
    ! we're just applying it to the Met data for now.

    ! Specify the intent of the arguments
    character(len=256), dimension(:), intent(in)  :: METFILES
    character(len=8),   dimension(:), intent(in)  :: METVARIABLES

    ! We need to have a definition of the substrings we're going to replace.
    character(len=11) :: StartDate = "<startdate>"
    character(len=9)  :: EndDate = "<enddate>"

    ! For clarity in the looping, and so we don't mutate the original inputs,
    ! we extract the filename from the array.
    character(len=256):: CurrentFile

    ! We write the outputs to a specific reserved filename
    character(len=32) :: ReservedOutputFile

    ! We read the start and end years to strings before writing them to the key.
    character(len=4)  :: StartYear, EndYear

    ! Initialise integers to store the status of the command.
    integer           :: ExStat, CStat

    ! We're going to want to remember where "<startdate>" and "<enddate>" occur
    ! in the filename templates, so we can retrieve the time periods of each
    ! file later.
    integer           :: IndxStartDate = 0, IndxEndDate = 0

    ! integers for the file IDs that read and write the file metadata.
    integer           :: InputUnit, OutputUnit

    ! And status integers for those files
    integer           :: iosIn, iosOut

    ! We want to reuse the units for the input and output files, so get unique
    ! IDs here rather than inside the loop.
    CALL get_unit(InputUnit)
    CALL get_unit(OutputUnit)

    ! Loop over the variable names and process the files for each variable.
    LoopVariables: DO VarIndx = 1, SIZE(METVARIABLES)
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
      ReservedOutputFile = METVARIABLES(VarIndx)//"Key.txt"

      ! Pull out the current file from the array
      CurrentFile = METFILES(VarIndx)

      ! Iterate through the string for <startdate>
      FindStart: DO CharIndx = 1, LEN(CurrentFile) - 11
        ! Check against <startdate>, if found set the index
        IF (CURRENTFILE(CharIndx:CharIndx+10) == StartDate) THEN
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
      OPEN(OutputUnit, ReservedOutputFile, ios = ios2)

      ! Now iterate through the lines in the input file, inspect the string
      ! to extract the bracketing years. We then append the start and the end
      ! year to the line. We can reuse the CurrentFile string to store the read
      ! line from the temporary file

      BuildKey: DO
        ! Read in the line and check the status
        READ(InputUnit, '(A)', IOSTAT = iosIn) CurrentFile

        IF (iosIn == -1) THEN
          ! Read completed successfully
          EXIT BuildKey
        ELSEIF (iosIn /= 0) THEN
          ! Read failed for some other reason
          EXIT(5)
        END If

        ! We assume the dates in the filename in the YYYYMMDD format.
        ! For now, we're assuming that the data files contain complete years,
        ! with partial year splits across the files, so all we need is the year.
        ! If in future we have more fine grained data, we can modify this and
        ! the key reader routines.

        ! Reading the start year is easy- just the first 4 characters starting
        ! from the recorded IndxStartDate.
        StartYear = CurrentFile(IndxStartDate:IndxStartDate+3)

        ! The end year is not as simple as taking the 4 characters at
        ! IndxEndDate, because the size of the <startdate> string is 3
        ! characters longer than the YYYYMMDD date specifier. Therefore, shift
        ! the index back 3 units and then take the next 4 characters
        EndYear = CurrentFile(IndxEndDate-3:IndxEndDate)

        ! Now we have everything we need to write back to the key file, in the
        ! format {filename} {StartYear} {EndYear}.

        WRITE(OutputUnit, '(A)', IOSTAT = iosOut) TRIM(CurrentFile)//" "//&
          StartYear//" "//EndYear

        ! Any non-zero iostat is bad
        IF (iosOut =/ 0) THEN
          EXIT(5)
        END IF

      END DO LoopVariables

      ! Now we've build the key for each variable, we can remove the temporary
      ! we used to store the filenames
      CALL execute_command_line("rm __tmpFileWithNoClashes__.txt")
