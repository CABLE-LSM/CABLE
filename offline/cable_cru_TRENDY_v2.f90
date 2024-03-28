! Author: Lachlan Whyborn
! Last Modified: Thu 28 Mar 2024 17:10:47

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

  subroutine PREPARE_MET_FILES(CRUFILES, CRUVARIABLES)

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

    ! We need to have a definition of the substrings we're going to replace.
    character(len=9)  :: StartDate = "<startdate>"
    character(len=7)  :: EndDate = "<enddate>"

    ! Initialise integers to store the status of the command.
    integer           :: ExitStat, CmdStat

    ! We're going to want to remember where "<startdate>" and "<enddate>" occur
    ! in the filename templates, so we can retrieve the time periods of each
    ! file.
    integer           :: IndxStartDate, IndxEndDate

    ! Loop over the variable names
    DO VarIndx = 1, SIZE(CRUFILES)
      ! First thing we need to do is make the strings compatible with ls, by
      ! replacing <startdate> and <enddate> with wildcard characters.

      ! First replace startdate, so only loop up to startdate's length from the
      ! end of the string. It seems to work fine even if we continue up to the 
      ! end of the string, but accessing past the end of the array is undefined
      ! behaviour and theoertically bad things could happen.
      DO CharIndx = 1, LEN(CRUFiles(VarIndx)) - 9
        ! If the local substring matches StartDate, replace it with the wildcard
        ! and shift the rest of the string 8 characters to the left.
        IF (CRUFiles(VarIndx)(CharIndx:CharIndx + 9) == StartDate) THEN
          CRUFiles(VarIndx)(CharIndx:CharIndx) = "*"
          CRUFiles(VarIndx)(CharIndx+1:) = CRUFiles(VarIndx)(CharIndx+9:)

          ! Remember where the start date starts
          IndxStartDate = CharIndx
        END IF
      END DO

      ! Do the same with enddate.
      DO CharIndx = 1, LEN(CRUFiles(VarIndx)) - 7
        ! Match and replace <enddate>
        IF (CRUFiles(VarIndx)(CharIndx:CharIndx + 7) == EndDate) THEN
          CRUFiles(VarIndx)(CharIndx:CharIndx) = "*"
          CRUFiles(VarIndx)(CharIndx+1:) = CRUFiles(VarIndx)(CharIndx+7:)
          
          ! Remember where the end date starts
          IndxEndDate = CharIndx
        END IF
      END DO

      ! Now we can use execute_command_lane to glob all the matching files and
      ! write them to a specified filename.
      call execute_command_lane("ls -1 "\\CRUFiles(VarIndx)\\" > "\\&
          CRUVariables(VarIndx)\\".txt", exitstat = ExitStat, cmdstat = CmdStat)

      
        





