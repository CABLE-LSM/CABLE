
!###############################################################################
!###############################################################################

MODULE bios_type_def
!------------------------------------------------------------------------------- 
! PRB: 18-06-2000
! This module explicitly defines the sizes of variable types
!-------------------------------------------------------------------------------
  implicit none
  save
! Define integer kind parameters to accommodate the range of numbers usually 
! associated with 4, 2, and 1 byte integers. 
  integer,parameter :: i4b = selected_int_kind(9) 
  integer,parameter :: i2b = selected_int_kind(4)
  integer,parameter :: i1b = selected_int_kind(2)
! Define single and double precision real kind parameters: 
! * Kind(1.0)   defines sp as the machine's default size for single precision
! * Kind(1.0d0) defines dp as the machine's default size for double precision
  integer,parameter :: sp  = kind(1.0)    
  integer,parameter :: dp  = kind(1.0d0)
! lgt is set to the default kind required for representing logical values. 
  integer,parameter :: lgt = kind(.true.)

END MODULE bios_type_def

!###############################################################################
!###############################################################################
    
MODULE bios_date_functions
!------------------------------------------------------------------------------
! Date Arithmetic (PRB, 12-07-2000; MRR, 11-oct-05)
!
! (PRB, 8/9/2015) Use operator overloading to add date arithmetic functionality
! to +, - and all logical operators so that dates of type DMYDATE can be added,
! subtracted, and compared. Note that this will work for logical operators
! regardless of which form is used in the code (e.g. >= or .GE.). Fortran 
! keywords capitalised to meet CABLE coding standards.
!------------------------------------------------------------------------------
USE bios_type_def
IMPLICIT NONE
PUBLIC
SAVE
! Define a TYPE for Australian-style dates built from 3 INTEGER fields.
! Assign as Date = dmydate(iday,imth,iyear)
! Access components as iday=Date%Day, imth=Date%Month, iyear=Date%Year
TYPE dmydate
  INTEGER(i4b):: Day
  INTEGER(i4b):: Month
  INTEGER(i4b):: Year
END TYPE dmydate
! Define interfaces for the +, -, ==, /=, <, >, <=, >= operators to:
! * add, subtract an integer number of days from a dmydate, with AddDay, SubDay; 
! * logically compare two dates, with EQDates, NEDates etc.
! These explicit interfaces are required.
INTERFACE OPERATOR (+)
  MODULE PROCEDURE AddDay
END INTERFACE
INTERFACE OPERATOR (-)
  MODULE PROCEDURE SubDay
END INTERFACE
INTERFACE OPERATOR (==)
  MODULE PROCEDURE EQDates
END INTERFACE
INTERFACE OPERATOR (/=)
  MODULE PROCEDURE NEDates
END INTERFACE
INTERFACE OPERATOR (<)
  MODULE PROCEDURE LTDates
END INTERFACE
INTERFACE OPERATOR (>)
  MODULE PROCEDURE GTDates
END INTERFACE
INTERFACE OPERATOR (<=)
  MODULE PROCEDURE LEDates
END INTERFACE
INTERFACE OPERATOR (>=)
  MODULE PROCEDURE GEDates
END INTERFACE
!-------------------------------------------------------------------------------

CONTAINS 

!*******************************************************************************
  
  FUNCTION AddDay(Today,Days2Add)
!-------------------------------------------------------------------------------
! Extends the '+' operator to enable adding days to a date
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE (dmydate)              :: AddDay    ! Date with days added (FUNCTION name)
  TYPE (dmydate), INTENT(IN)  :: Today     ! Current date
  INTEGER(i4b),   INTENT(IN)  :: Days2Add  ! Days to add to current date
  REAL(dp)                    :: JDay      ! Julian Day
!-------------------------------------------------------------------------------
  JDay = JulianDay(Today)
  JDay = JDay + FLOAT(Days2Add)
  AddDay = GregDate(JDay)
  END FUNCTION AddDay

!*******************************************************************************

  FUNCTION SubDay(Today,Days2Sub)
!-------------------------------------------------------------------------------
! Extends the '-' operator to enable subtracting days from a date
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE (dmydate)              :: SubDay    ! Date with days added (FUNCTION name)
  TYPE (dmydate), INTENT(IN)  :: Today     ! Current date
  INTEGER(i4b),   INTENT(IN)  :: Days2Sub  ! Days to add to current date
  REAL(dp)                    :: JDay      ! Julian Day
!-------------------------------------------------------------------------------
  JDay = JulianDay(Today)
  JDay = JDay - FLOAT(Days2Sub)
  SubDay = GregDate(JDay)
  END FUNCTION SubDay

!*******************************************************************************

  FUNCTION EQDayofYear(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '==' operator to enable comparison of two dates for equality.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: EQDayofYear   ! Result of equality test between dates
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!-------------------------------------------------------------------------------
  EQDayofYear = .true.
  IF (Date1%day   /= Date2%day   .OR. &
      Date1%month /= Date2%month ) EQDayofYear = .false.
  END FUNCTION EQDayofYear
!*******************************************************************************

  FUNCTION EQDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '==' operator to enable comparison of two dates for equality.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: EQDates   ! Result of equality test between dates
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!-------------------------------------------------------------------------------
  EQDates = .true.
  IF (Date1%day   /= Date2%day   .OR. &
      Date1%month /= Date2%month .OR. &
      Date1%year  /= Date2%year) EQDates = .false.
  END FUNCTION EQDates

!*******************************************************************************

  FUNCTION NEDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '/=' operator to enable comparison of two dates for inequality.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: NEDates   ! Result of inequality test between dates
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!-------------------------------------------------------------------------------  
  NEDates = .not. EQDates(Date1,Date2)
  END FUNCTION NEDates

!*******************************************************************************

  FUNCTION LTDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '<' operator to enable comparison of two dates for Date1 < Date2.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: LTDates   ! Result of LT test between dates
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!------------------------------------------------------------------------------- 
  LTDates = (JulianDay(Date1)<JulianDay(Date2))
  END FUNCTION LTDates

!*******************************************************************************

  FUNCTION GTDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '>' operator to enable comparison of two dates for Date1 > Date2.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: GTDates   ! Result of GT test between dates
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!-------------------------------------------------------------------------------  
  GTDates = (JulianDay(Date1)>JulianDay(Date2))
  END FUNCTION GTDates

!*******************************************************************************

  FUNCTION LEDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '<=' operator to enable comparison of two dates Date1 <= Date2.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: LEDates   ! Result of LE test between dates
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!-------------------------------------------------------------------------------  
  LEDates = (LTDates(Date1,Date2) .OR. EQDates(Date1,Date2))
  END FUNCTION LEDates

!*******************************************************************************

  FUNCTION GEDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '>=' operator to enable comparison of two dates Date1 >= Date2.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: GEDates   ! Result of GE test between dates
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!-------------------------------------------------------------------------------  
  GEDates = (GTDates(Date1,Date2) .OR. EQDates(Date1,Date2))
  END FUNCTION GEDates
!*******************************************************************************

  FUNCTION EQMonths(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '==' operator to enable comparison of two time-series months for equality.
! VH, 25-03-11
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: EQMonths   ! Result of equality test between months
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!-------------------------------------------------------------------------------
  EQMonths = .true.
  IF (Date1%month /= Date2%month .OR. &
     Date1%year  /= Date2%year) EQMonths = .false.
  END FUNCTION EQMonths
!*******************************************************************************

  FUNCTION EQMonthofYear(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '==' operator to enable comparison of two months of year for equality.
! VH, 25-03-11
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL                     :: EQMonthofYear   ! Result of equality test between months
  TYPE (dmydate), INTENT(IN)  :: Date1     ! First date for LOGICAL comparison
  TYPE (dmydate), INTENT(IN)  :: Date2     ! Second date for LOGICAL comparison
!-------------------------------------------------------------------------------
  EQMonthofYear = .true.
  IF (Date1%month /= Date2%month) EQMonthofYear = .false.
  END FUNCTION EQMonthofYear
!*******************************************************************************

  FUNCTION JulianDay(GregDate)
!-------------------------------------------------------------------------------
! Returns a REAL(sp) Julian Day when given a Gregorian Date.
! Adapted from Date Algorithms of Peter Baum:
! http://vsg.cape.com/~pbaum/date/date0.htm
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(dp)                   :: JulianDay
  TYPE (dmydate), INTENT(IN) :: GregDate
  REAL(dp)                   :: D,M,Y
!-------------------------------------------------------------------------------
  D = real(GregDate%Day,dp)
  M = real(GregDate%Month,dp)
  Y = real(GregDate%Year,dp)
  IF (M<3.0_dp) THEN 
    M = M + 12.0_dp 
    Y = Y - 1.0_dp 
  END IF 
  JulianDay = D + int((153.0_dp * M - 457.0_dp)/5.0_dp) + (365.0_dp*Y) +     &
              FLOOR(Y/4.0_dp) - FLOOR(Y/100.0_dp) + FLOOR(Y/400.0_dp) + 1721118.5_dp
  END FUNCTION JulianDay

!*******************************************************************************

  FUNCTION GregDate(JulianDay)
!-------------------------------------------------------------------------------
! Returns a Gregorian Date when given a Julian Day 
! Modified from Date Algorithms of Peter Baum 
! http://vsg.cape.com/~pbaum/date/date0.htm
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE (dmydate)             :: GregDate
  REAL(dp), INTENT(IN)       :: JulianDay
  REAL(dp)                   :: D,M,Y
  REAL(dp)                   :: Z,R,A,G,B,C
!-------------------------------------------------------------------------------
  Z = FLOOR(JulianDay - 1721118.5_dp) 
  R = JulianDay - 1721118.5_dp - Z
  G = Z - 0.25_dp 
  A = FLOOR(G/36524.25_dp) 
  B = A - FLOOR(A/4.0_dp) 
  Y = FLOOR((B+G)/365.25_dp) 
  C = B + Z - FLOOR(365.25_dp*Y) 
  M = int(((5.0_dp*C) + 456.0_dp) / 153.0_dp) 
  D = C - int((153.0_dp * M - 457.0_dp) / 5.0_dp) + R  
  IF (M > 12.0_dp) THEN 
    Y = Y + 1.0_dp 
    M = M - 12.0_dp 
  END IF
  GregDate = dmydate(int(D),int(M),int(Y))  ! Keep truncated day number only
  END FUNCTION GregDate

!*******************************************************************************

  FUNCTION LeapDay(Year)
!-------------------------------------------------------------------------------
! Returns 1 if leap year, 0 if not.  Add it to the length of February.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(i4b)               :: LeapDay  
  INTEGER(i4b), INTENT(IN)   :: Year
!-------------------------------------------------------------------------------
  IF (MOD(Year,4)/=0) THEN
    LeapDay = 0
  ELSE IF (MOD(Year,400)==0) THEN
    LeapDay = 1
  ELSE IF (MOD(Year,100)==0) THEN
    LeapDay = 0
  ELSE
    LeapDay = 1
  END IF
  END FUNCTION LeapDay

!*******************************************************************************

FUNCTION DaysInMonth(Date)
!-------------------------------------------------------------------------------
! returns number of days in current month
! MRR, 12-oct-2005: return 0 if month not legal
!-------------------------------------------------------------------------------
TYPE(dmydate),INTENT(IN):: Date
INTEGER(i4b)          :: DaysInMonth  
INTEGER(i4b),PARAMETER:: MonthDays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
!-------------------------------------------------------------------------------
IF (Date%Month >= 1 .AND. Date%Month <= 12) THEN
  IF (Date%Month==2) THEN
      DaysInMonth = MonthDays(2) + LeapDay(Date%Year)
  ELSE 
    DaysInMonth = MonthDays(Date%Month)
  END IF
ELSE
  DaysInMonth = 0
END IF
END FUNCTION DaysInMonth

!*******************************************************************************

  FUNCTION YearDay(Date)
!-------------------------------------------------------------------------------
  TYPE(dmydate), INTENT(IN)   :: Date
  INTEGER(i4b)                :: YearDay 
  INTEGER(i4b), DIMENSION(12) :: MonthDays
!-------------------------------------------------------------------------------
  MonthDays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  MonthDays(2) = 28 + LeapDay(Date%Year)
  IF (Date%Month==1) THEN
    YearDay = Date%Day
  ELSE
    YearDay = SUM(MonthDays(1:Date%Month-1)) + Date%Day
  END IF
  END FUNCTION YearDay

!*******************************************************************************

FUNCTION DayDifference (Date1,Date0)
!-------------------------------------------------------------------------------
! Returns Date1 - Date0 in INTEGER days.
! MRR, 11-oct-05
!-------------------------------------------------------------------------------
USE bios_type_def
IMPLICIT NONE
TYPE(dmydate),INTENT(IN):: Date1, Date0
INTEGER(i4b):: DayDifference
!-------------------------------------------------------------------------------
DayDifference = nint(JulianDay(Date1) - JulianDay(Date0))
END FUNCTION DayDifference

!*******************************************************************************

FUNCTION LegalDate (Date)
!-------------------------------------------------------------------------------
! Returns .true. if date is legal (test D,M only), otherwise .false.
! MRR, 11-oct-05
!-------------------------------------------------------------------------------
USE bios_type_def
IMPLICIT NONE
TYPE(dmydate),INTENT(IN):: Date
LOGICAL(lgt):: LegalDate
!-------------------------------------------------------------------------------
LegalDate = .false.
IF ( (Date%month >= 1 .AND. Date%month <= 12) .AND.         &   ! check month
     (Date%day >= 1 .AND. Date%day <= DaysInMonth(Date)) )  &   ! check day
   LegalDate = .true.
END FUNCTION LegalDate

!*******************************************************************************

FUNCTION EndMonth (Date)
!-------------------------------------------------------------------------------
! Returns .true. if date is last day of month, otherwise .false.
! MRR, 11-oct-05
!-------------------------------------------------------------------------------
USE bios_type_def
IMPLICIT NONE
TYPE(dmydate),INTENT(IN):: Date
LOGICAL(lgt):: EndMonth
!-------------------------------------------------------------------------------
EndMonth = .false.
IF (Date%day == DaysInMonth(Date)) EndMonth = .true.
END FUNCTION EndMonth
!*******************************************************************************

FUNCTION BeginMonth (Date)
!-------------------------------------------------------------------------------
! Returns .true. if date is first day of month, otherwise .false.
! VH, 25-03-11
!-------------------------------------------------------------------------------
USE bios_type_def
IMPLICIT NONE
TYPE(dmydate),INTENT(IN):: Date
LOGICAL(lgt):: BeginMonth
!-------------------------------------------------------------------------------
BeginMonth = .false.
IF (Date%day == 1) BeginMonth = .true.
END FUNCTION BeginMonth

!*******************************************************************************

FUNCTION EndYear (Date)
!-------------------------------------------------------------------------------
! Returns .true. if date is last day of year, otherwise .false.
! MRR, 11-oct-05
!-------------------------------------------------------------------------------
USE bios_type_def
IMPLICIT NONE
TYPE(dmydate),INTENT(IN):: Date
LOGICAL(lgt):: EndYear
!-------------------------------------------------------------------------------
EndYear = .false.
IF (Date%day == 31 .AND. Date%month == 12) EndYear = .true.
END FUNCTION EndYear

!*******************************************************************************

END MODULE bios_date_functions
    
!*******************************************************************************
!*******************************************************************************
    
MODULE bios_misc_io
    
USE bios_type_def

CONTAINS

  SUBROUTINE ReadArcFltHeader(iunit,Filename,Cols,Rows,xLL,yLL,CellSize,NoDataVal)
!-------------------------------------------------------------------------------
! Read the header file (.hdr) associated with an ArcGIS binary grid file (.hdr).
! The elements are req'd in standard order, although ArcGIS is not reliant on this
! order. A seventh line is often found specifying the endianness, but this is not
! read. Return the grid dimensions and no data value.  
!-------------------------------------------------------------------------------
  implicit none
  integer(i4b)    ,intent(in) :: iunit
  character(200)  ,intent(in) :: Filename
  integer(i4b)    ,intent(out):: Cols, Rows
  real(sp)        ,intent(out):: xLL, yLL, CellSize, NoDataVal
  character(12)               :: Head
!-------------------------------------------------------------------------------
  open (unit=iunit,file=Filename,status='old')

  read (iunit,*) Head, Cols      ! Number of columns
  read (iunit,*) Head, Rows      ! Number of rows
  read (iunit,*) Head, xLL       ! Western boundary
  read (iunit,*) Head, yLL       ! Southern boundary
  read (iunit,*) Head, CellSize  ! Resolution (both W-E, N-S)
  read (iunit,*) Head, NoDataVal ! Missing data value
  close (unit=iunit)

  END SUBROUTINE ReadArcFltHeader

END MODULE bios_misc_io

!*******************************************************************************
!*******************************************************************************
    
MODULE cable_bios_met_obs_params

! Input routines for BIOS meteorology, remote sensing observations, and soil  
! parameters for use in the CABLE land surface scheme.

  USE bios_type_def
  USE bios_date_functions
  USE bios_misc_io
  
  USE CABLE_COMMON_MODULE, ONLY: GET_UNIT, cable_user ! Subr assigns an unique unit number for file opens.

  USE cable_IO_vars_module, ONLY: &
      logn, landpt , nmetpatches         ! Unit number for writing logfile entries and land point array
  USE cable_def_types_mod,   ONLY: MET_TYPE, soil_parameter_type, mland
  
  USE cable_weathergenerator,ONLY: WEATHER_GENERATOR_TYPE, WGEN_INIT, &
                                   WGEN_DAILY_CONSTANTS, WGEN_SUBDIURNAL_MET

  IMPLICIT NONE

! Define the filenames from the namelist: pathnames, landmask, daily bios .bin met files, and soil parameter files
! We define them across the whole module because they are read in by cable_bios_init but used later by
! cable_bios_load_params, i.e. outside of the bios initialisation process. We do this separately because in cable_driver  
! parameter initialisation happens after other initialisations. Our choice is to overwrite the defaults, so this has to be
! done after the defaults are read in. 

  CHARACTER(200) :: Run            ! 'spinup', 'premet', 'standard'
                                   ! 'spinup' uses recycled met (from later date), and preind CO2 and Ndep
                                   ! 'premet' uses recycled met (from later date), and actual CO2 and Ndep (e.g. appropriate before met starts in 1900)
                                   ! 'standard' uses actual met, and actual CO2 and Ndep
  CHARACTER(len=15) :: CO2         ! CO2 takes value       : "preind","actual"
  CHARACTER(len=15) :: Ndep        ! Ndep takes value      : "preind","actual"
  CHARACTER(len=15) :: MetForcing  ! MetForcing takes value: "recycled", "actual"


  CHARACTER(200) :: met_path, param_path, landmaskflt_file, landmaskhdr_file, &
       rain_file, swdown_file, tairmax_file, tairmin_file, wind_file, &
       vp0900_file, vp1500_file, co2_file, &
    b1_file, b2_file, bulkdens1_kgm3_file, bulkdens2_kgm3_file, clayfrac1_file, clayfrac2_file, &
    csoil1_file, csoil2_file, depth1_m_file, depth2_m_file, hyk1sat_ms_file, hyk2sat_ms_file, & 
    psie1_m_file, psie2_m_file, siltfrac1_file, siltfrac2_file, wvol1fc_m3m3_file, wvol2fc_m3m3_file, &
    wvol1sat_m3m3_file, wvol2sat_m3m3_file, wvol1w_m3m3_file, wvol2w_m3m3_file, MVG_file , &
    c4frac_file , vegtypeigbp_file, avgannmax_fapar_file 
    !, &
!    slope_deg_file  ! Terrain slope in degrees

! Define the daily bios met variables and file unit numbers which are required by multiple subroutines.
  
  REAL(sp),     ALLOCATABLE :: rain_day(:)          ! Packed vector of daily AWAP/BIOS rain (mm) 
  REAL(sp),     ALLOCATABLE :: swdown_day(:)        ! Packed vector of daily AWAP/BIOS swdown (MJ)
  REAL(sp),     ALLOCATABLE :: tairmax_day(:)       ! Packed vector of daily AWAP/BIOS max air temp (deg C)
  REAL(sp),     ALLOCATABLE :: tairmin_day(:)       ! Packed vector of daily AWAP/BIOS min air temp (deg C)
  REAL(sp),     ALLOCATABLE :: prev_tairmax_day(:)       ! Packed vector of previous day's AWAP/BIOS max air temp (deg C)
  REAL(sp),     ALLOCATABLE :: next_tairmin_day(:)       ! Packed vector of next day's AWAP/BIOS min air temp (deg C)
  REAL(sp),     ALLOCATABLE :: wind_day(:)          ! Packed vector of daily wind (ms-1)
  REAL(sp),     ALLOCATABLE :: vp0900(:)          ! Packed vector of 9am vapour pressure (mb)
  REAL(sp),     ALLOCATABLE :: vp1500(:)          ! Packed vector of 3pm vapour pressure (mb)
  REAL(sp),     ALLOCATABLE :: prev_vp1500(:)          ! Packed vector of 3pm vapour pressure (mb)
  REAL(sp),     ALLOCATABLE :: next_vp0900(:)          ! Packed vector of 9am vapour pressure (mb)

  REAL(sp)                  :: co2air_year          ! Single global value of co2 (ppm)
  
  INTEGER(i4b),  SAVE :: rain_unit, swdown_unit, tairmax_unit, tairmin_unit, co2_unit ! Met file unit numbers
  INTEGER(i4b),  SAVE :: wind_unit, vp1500_unit, vp0900_unit
  TYPE(dmydate), SAVE :: previous_date ! The day before the current date, for noting changes of year
  TYPE(dmydate), SAVE :: bios_rundate  ! The day before the current date, for noting changes of year
  TYPE(dmydate)       :: dummydate     ! Dummy date for when keeping the date is not required
  TYPE(dmydate),SAVE  :: MetDate       ! Date of met to access (equals current date for normals runs, but
                                       ! must be calculated for spinup and initialisation runs (for dates before 1900)
  INTEGER(i4b),PARAMETER :: recycle_met_startdate = 1981 ! range for met to be recycled for spinup and initialisation
  INTEGER(i4b),PARAMETER :: recycle_met_enddate = 2010
  INTEGER(i4b)   :: skipdays                        ! Days of met to skip when user_startdate is after bios_startdate
  TYPE(dmydate), SAVE  :: bios_startdate, bios_enddate    ! First and last dates found in bios met files (read from rain file)
  REAL(sp), PRIVATE, PARAMETER :: SecDay = 86400.
  TYPE(WEATHER_GENERATOR_TYPE), SAVE :: WG

CONTAINS

  SUBROUTINE cable_bios_init(dels,curyear,met,kend,ktauday)

  USE cable_IO_vars_module, ONLY: & 
    latitude, longitude, &  ! Vectors for lat and long of each land cell
    land_y, land_x,       &  ! Vectors for row and col of each land cell
    nmetpatches,         &
    mask,                &  ! Integer land mask
    metGrid, &
    sdoy, smoy, syear, shod, &
    xdimsize, ydimsize, &
    lat_all, lon_all

  IMPLICIT NONE

  REAL, INTENT(INOUT) :: dels  ! time step size in seconds, read from bios.nml namelist
  INTEGER, INTENT(INOUT) :: curyear, kend 
  INTEGER, INTENT(IN) :: ktauday
  TYPE(MET_TYPE), INTENT(INOUT)       :: MET 
  ! Local variables
  
  LOGICAL,  SAVE :: call1 = .TRUE.
  INTEGER(i4b)   :: iunit  
  INTEGER(i4b)   :: MaskCols, MaskRows  ! Landmask col and row dimensions
  REAL(sp)       :: MaskBndW, MaskBndS  ! Landmask outer bound dimensions in decimal degrees (West & South)
  REAL(sp)       :: MaskCtrW, MaskCtrS
  REAL(sp)       :: MaskRes, NoDataVal  ! Landmask resolution (dec deg) and no-data value
  
  INTEGER(i4b)   :: icol, irow, iland   ! Loop counters for cols, rows, land cells
  
  LOGICAL(lgt), ALLOCATABLE :: LandMaskLogical(:,:) ! Logical land mask for vector packing
  REAL(sp),     ALLOCATABLE :: LandMaskReal(:,:)    ! Real land mask as read from landmask file
  INTEGER(i4b), ALLOCATABLE :: ColRowGrid(:,:)      ! Temp grid to hold col or row numbers for packing to land_x or land_y
  
 ! TYPE(dmydate)  :: bios_startdate, bios_enddate    ! First and last dates found in bios met files (read from rain file)
  TYPE(dmydate)  :: user_startdate, user_enddate    ! First and last run dates specified by user (read from cable_user%)
  INTEGER(i4b)   :: bios_rundays                    ! Days in run, from bios_startdate or user_startdate to bios_enddate 
  INTEGER(i4b)   :: co2_startyear, co2_endyear      ! First and last years found in bios global CO2 files 
  INTEGER(i4b)   :: co2_skipyears                   ! Years of annual CO2 to skip to position for reading the current year
  INTEGER(i4b)   :: iday, iyear ! counters
  INTEGER :: error_status
  NAMELIST /biosnml/ Run, met_path, param_path, landmaskflt_file, landmaskhdr_file, &
       rain_file, swdown_file, tairmax_file, tairmin_file, &
       wind_file, &
       vp0900_file, vp1500_file, co2_file, &
    b1_file, b2_file, bulkdens1_kgm3_file, bulkdens2_kgm3_file, clayfrac1_file, clayfrac2_file, &
    csoil1_file, csoil2_file, depth1_m_file, depth2_m_file, hyk1sat_ms_file, hyk2sat_ms_file, & 
    psie1_m_file, psie2_m_file, siltfrac1_file, siltfrac2_file, wvol1fc_m3m3_file, wvol2fc_m3m3_file, &
    wvol1sat_m3m3_file, wvol2sat_m3m3_file, wvol1w_m3m3_file, wvol2w_m3m3_file, MVG_file, &
    c4frac_file, vegtypeigbp_file, avgannmax_fapar_file, &
!    slope_deg_file, &  ! Terrain slope in degrees
    dels               ! time step size in seconds
  
  ! Open and read the BIOS namelist file.
  
  CALL GET_UNIT(iunit)
  OPEN (iunit, FILE="bios.nml")
  READ (iunit, NML=biosnml)
  CLOSE(iunit)

END SUBROUTINE cable_bios_init
!******************************************************************************

SUBROUTINE cable_bios_load_params(soil)

USE cable_def_types_mod,  ONLY: mland, soil_parameter_type

IMPLICIT NONE

INTEGER(i4b) :: is, ie ! Index start/end points within cable spatial vectors
                                ! for the current land-cell's tiles. These are just 
                                ! aliases to improve code readability
INTEGER(i4b) :: iland         ! loop counter through mland land cells
INTEGER(i4b) :: param_unit    ! Unit number for reading (all) parameter files.
INTEGER(i4b) :: error_status  ! Error status returned by OPENs

TYPE(soil_parameter_type), INTENT(INOUT) :: soil


! Temporary soil parameter variables. Varnames match corresponding bios parameter filenames (roughly).
! Dimensions are mland, which will be mapped to the more complicated mland+tiles dimensions of the 
! equivalent cable soil variables. Filenames for these variables are defined at module level
! because they are read in from bios.nml as part of the initialisation but only accessed here,
! after the default cable params are read in by cable_driver. 
REAL(sp), ALLOCATABLE :: b1(:)
REAL(sp), ALLOCATABLE :: b2(:)
REAL(sp), ALLOCATABLE :: bulkdens1_kgm3(:)
REAL(sp), ALLOCATABLE :: bulkdens2_kgm3(:)
REAL(sp), ALLOCATABLE :: clayfrac1(:)
REAL(sp), ALLOCATABLE :: clayfrac2(:)
REAL(sp), ALLOCATABLE :: csoil1(:)
REAL(sp), ALLOCATABLE :: csoil2(:)
REAL(sp), ALLOCATABLE :: depth1_m(:)
REAL(sp), ALLOCATABLE :: depth2_m(:)
REAL(sp), ALLOCATABLE :: hyk1sat_ms(:)
REAL(sp), ALLOCATABLE :: hyk2sat_ms(:)
REAL(sp), ALLOCATABLE :: psie1_m(:)
REAL(sp), ALLOCATABLE :: psie2_m(:)
REAL(sp), ALLOCATABLE :: siltfrac1(:)
REAL(sp), ALLOCATABLE :: siltfrac2(:)
REAL(sp), ALLOCATABLE :: wvol1fc_m3m3(:)
REAL(sp), ALLOCATABLE :: wvol2fc_m3m3(:)
REAL(sp), ALLOCATABLE :: wvol1sat_m3m3(:)
REAL(sp), ALLOCATABLE :: wvol2sat_m3m3(:)
REAL(sp), ALLOCATABLE :: wvol1w_m3m3(:)
REAL(sp), ALLOCATABLE :: wvol2w_m3m3(:)
! REAL(sp), ALLOCATABLE :: slope_deg(:)


! Allocate all bios soil parameters to mland, the number of land cells.
ALLOCATE (b1(mland), b2(mland), bulkdens1_kgm3(mland), bulkdens2_kgm3(mland))
ALLOCATE (clayfrac1(mland), clayfrac2(mland), csoil1(mland), csoil2(mland))
ALLOCATE (depth1_m(mland), depth2_m(mland), hyk1sat_ms(mland), hyk2sat_ms(mland))
ALLOCATE (psie1_m(mland), psie2_m(mland), siltfrac1(mland), siltfrac2(mland))
ALLOCATE (wvol1fc_m3m3(mland), wvol2fc_m3m3(mland), wvol1sat_m3m3(mland), wvol2sat_m3m3(mland))
ALLOCATE (wvol1w_m3m3(mland), wvol2w_m3m3(mland))
!ALLOCATE (slope_deg(mland))

CALL GET_UNIT(param_unit)  ! Obtain an unused unit number for file reading, reused for all soil vars.

! Open, read, and close each soil parameter in turn, stopping for any missing file.
OPEN (param_unit, FILE=TRIM(param_path)//TRIM(b1_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
WRITE(*,*) "Param path: ", TRIM(param_path), " b1 file: ", b1_file
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(b1_file) ; STOP 1
ELSE
  READ (param_unit) b1
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(b2_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(b2_file) ; STOP 2
ELSE
  READ (param_unit) b2
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(bulkdens1_kgm3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(bulkdens1_kgm3_file) ; STOP 3
ELSE
  READ (param_unit) bulkdens1_kgm3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(bulkdens2_kgm3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(bulkdens2_kgm3_file) ; STOP 4
ELSE
  READ (param_unit) bulkdens2_kgm3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(clayfrac1_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(clayfrac1_file) ; STOP 5
ELSE
  READ (param_unit) clayfrac1
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(clayfrac2_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(clayfrac2_file) ; STOP 6
ELSE
  READ (param_unit) clayfrac2
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(csoil1_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(csoil1_file) ; STOP 7
ELSE
  READ (param_unit) csoil1
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(csoil2_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(csoil2_file) ; STOP 8
ELSE
  READ (param_unit) csoil2
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(depth1_m_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(depth1_m_file) ; STOP 9
ELSE
  READ (param_unit) depth1_m
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(depth2_m_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(depth2_m_file) ; STOP 10
ELSE
  READ (param_unit) depth2_m
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(hyk1sat_ms_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(hyk1sat_ms_file) ; STOP 11
ELSE
  READ (param_unit) hyk1sat_ms
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(hyk2sat_ms_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(hyk2sat_ms_file) ; STOP 12
ELSE
  READ (param_unit) hyk2sat_ms
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(psie1_m_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(psie1_m_file); STOP 13
ELSE
  READ (param_unit) psie1_m
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(psie2_m_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(psie2_m_file); STOP 14
ELSE
  READ (param_unit) psie2_m
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(siltfrac1_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(siltfrac1_file); STOP 15
ELSE
  READ (param_unit) siltfrac1
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(siltfrac2_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(siltfrac2_file) ; STOP 16
ELSE
  READ (param_unit) siltfrac2
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol1fc_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol1fc_m3m3_file); STOP 17
ELSE
  READ (param_unit) wvol1fc_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol2fc_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)   
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol2fc_m3m3_file); STOP 18
ELSE
  READ (param_unit) wvol2fc_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol1sat_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol1sat_m3m3_file) ; STOP 19
ELSE
  READ (param_unit) wvol1sat_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol2sat_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol2sat_m3m3_file) ; STOP 20
ELSE
  READ (param_unit) wvol2sat_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol1w_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol1w_m3m3_file) ; STOP 21
ELSE
  READ (param_unit) wvol1w_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol2w_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol2w_m3m3_file) ; STOP 22
ELSE
  READ (param_unit) wvol2w_m3m3
  CLOSE (param_unit)
END IF

!OPEN (param_unit, FILE=TRIM(param_path)//TRIM(slope_deg_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
!IF (error_status > 0) THEN
!  WRITE (*,'("STOP - File not found: ",  LEN_TRIM(TRIM(param_path)//TRIM(slope_deg_file)))') ; STOP -1
!ELSE
!  READ (param_unit) slope_deg
!  CLOSE (param_unit)
!END IF

! Map the values of each soil parameter read from bios parameter files (one value
! for each of mland land cells) onto the equivalent CABLE parameter, overwriting
! the default values and doing unit and other conversions as required. Assign the 
! same bios parameter value for each land cell to all tiles within that land cell. 
! 
! Explanation: CABLE parameters have multiple 'tiles' (eg. veg types) per land cell.
! The number of tiles may vary with each land cell so to save memory they
! are incorporated within the same dimension as mland. Within the larger vector
! an index keeps track of the location of the first and last tile of each land cell.
! This index is landpt, which for each of mland cells has as cstart and cend index.

DO iland = 1,mland ! For each land cell...
  is = landpt(iland)%cstart  ! Index position for the first tile of this land cell.
  ie = landpt(iland)%cend    ! Index position for the last tile of this land cell.
!%%%%%  soil%albsoil(is:ie) = albsoil(iland)
  soil%bch(is:ie)     = min(b1(iland),16.0) 
  soil%silt(is:ie)    = siltfrac1(iland)
  soil%clay(is:ie)    = clayfrac1(iland)
  soil%sand(is:ie)    = 1.0 - soil%silt(is:ie) - soil%clay(is:ie)
  soil%css(is:ie)     = csoil1(iland)
  soil%hyds(is:ie)    = max(hyk1sat_ms(iland),1.0e-8) 
  soil%sfc(is:ie)     = wvol1fc_m3m3(iland)
  soil%ssat(is:ie)    = max(0.4,wvol1sat_m3m3(iland))
  soil%sucs(is:ie)    = max(psie1_m(iland),-2.0)         
  soil%swilt(is:ie)   = min(0.2,wvol1w_m3m3(iland)) 
  soil%rhosoil(is:ie) = bulkdens1_kgm3(iland)

! PRB: Comment this out because CABLE is not using spatial zse. The dimension is ms (soil levels), not mp.
! Remove this section when both horizons are dealt with in the commented out code below 
!  print *,"is, ie ",is,ie
!  soil%zse(is:ie)   = depth1_m(iland) 
!  WHERE (soil%zse(is:ie).lt.0.01)    
!    soil%zse(is:ie) = 0.20 ! for compatibility with AWAP
!  endwhere

!%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSULT WITH VANESSA ABOUT HOW TO DEAL WITH B HORIZON
  ! B horizon parameters (for use in soil-litter)
 ! soil%bchB(is:ie)     = min(b2(iland),16.0)  
 ! soil%siltB(is:ie)    = siltfrac2(iland)
 ! soil%clayB(is:ie)    = clayfrac2(iland)
 ! soil%cssB(is:ie)     = csoil2(iland)
 ! soil%hydsB(is:ie)    = max(hyk2sat_ms(iland),1.0e-8)  ! m s-1  
 ! soil%ssatB(is:ie)    = max(0.4,wvol2sat_m3m3(iland))
 ! soil%sucsB(is:ie)    = max(psie2_m(iland),-2.0)          
 ! soil%swiltB(is:ie)   = min(0.2,wvol2w_m3m3(iland))
 ! soil%rhosoilB(is:ie) = bulkdens2_kgm3(iland)
 ! soil%sfcB(is:ie)     = wvol2fc_m3m3(iland)
 ! soil%depthA(is:ie)   = depth1_m(iland)
 ! WHERE (soil%depthA(is:ie).lt.0.01)    
 !   soil%depthA(is:ie) = 0.20 ! for compatibility with AWAP
 ! endwhere
 ! soil%depthB(is:ie)   = depth2_m(iland)
 ! WHERE (soil%depthB(is:ie).lt.0.01)
 !   soil%depthB(is:ie) = 1.00 ! to avoid very small B-Horizon
 ! endwhere
!%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSULT WITH VANESSA ABOUT HOW TO DEAL WITH B HORIZON

!%%%%%%%%%%%%%%%%%%%%%%%%%%% Assign terrain slope variable to new cable variable (tbd) here
! cablevar%slope_deg(is:ie)  = slope_deg(iland)

END DO

! Deallocate soil parameter variables which have been copied to their CABLE equivalents.
DEALLOCATE (b1, b2, bulkdens1_kgm3, bulkdens2_kgm3)
DEALLOCATE (clayfrac1, clayfrac2, csoil1, csoil2)
DEALLOCATE (depth1_m, depth2_m, hyk1sat_ms, hyk2sat_ms)
DEALLOCATE (psie1_m, psie2_m, siltfrac1, siltfrac2)
DEALLOCATE (wvol1fc_m3m3, wvol2fc_m3m3, wvol1sat_m3m3, wvol2sat_m3m3)
DEALLOCATE (wvol1w_m3m3, wvol2w_m3m3)
!DEALLOCATE (slope_deg)

END SUBROUTINE cable_bios_load_params

!******************************************************************************

SUBROUTINE cable_bios_load_biome(MVG)

  USE cable_def_types_mod,  ONLY: mland


  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: MVG(:) ! climate variables

  INTEGER(i4b) :: param_unit    ! Unit number for reading (all) parameter files.
  INTEGER(i4b) :: error_status  ! Error status returned by OPENs
  REAL(sp), ALLOCATABLE :: tmp(:)

  ! Temporary soil parameter variables. Varnames match corresponding bios parameter filenames (roughly).
  ! Dimensions are mland, which will be mapped to the more complicated mland+tiles dimensions of the 
  ! equivalent cable soil variables. Filenames for these variables are defined at module level
  ! because they are read in from bios.nml as part of the initialisation but only accessed here,
  ! after the default cable params are read in by cable_driver. 

  ALLOCATE (tmp(mland))

  CALL GET_UNIT(param_unit)  ! Obtain an unused unit number for file reading, reused for all soil vars.

  ! Open, read, and close Major Veg Group file
  OPEN (param_unit, FILE=TRIM(param_path)//TRIM(MVG_file), ACCESS='STREAM', &
       FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
  IF (error_status > 0) THEN
     WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(MVG_file) ; STOP -1
  ELSE
     READ (param_unit) tmp
     CLOSE (param_unit)
  END IF

  MVG = int(tmp)

END SUBROUTINE cable_bios_load_biome


!******************************************************************************

SUBROUTINE cable_bios_load_fracC4(fracC4)

  USE cable_def_types_mod,  ONLY: mland


  IMPLICIT NONE

  REAL, INTENT(INOUT) :: fracC4(:) ! climate variables

  INTEGER(i4b) :: param_unit    ! Unit number for reading (all) parameter files.
  INTEGER(i4b) :: error_status  ! Error status returned by OPENs
  REAL(sp), ALLOCATABLE :: tmp(:)

  ! Temporary soil parameter variables. Varnames match corresponding bios parameter filenames (roughly).
  ! Dimensions are mland, which will be mapped to the more complicated mland+tiles dimensions of the 
  ! equivalent cable soil variables. Filenames for these variables are defined at module level
  ! because they are read in from bios.nml as part of the initialisation but only accessed here,
  ! after the default cable params are read in by cable_driver. 

  ALLOCATE (tmp(mland))

  CALL GET_UNIT(param_unit)  ! Obtain an unused unit number for file reading

  ! Open, read, and close Major Veg Group file
  OPEN (param_unit, FILE=TRIM(param_path)//TRIM(c4frac_file), ACCESS='STREAM', &
       FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
  IF (error_status > 0) THEN
     WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(c4frac_file) ; STOP -1
  ELSE
     READ (param_unit) tmp
     CLOSE (param_unit)
  END IF
  
  fracC4 = tmp

END SUBROUTINE cable_bios_load_fracC4

!******************************************************************************
SUBROUTINE cable_bios_load_climate_params(climate)

USE cable_def_types_mod,  ONLY: mland, climate_type

  IMPLICIT NONE


INTEGER(i4b) :: is, ie ! Index start/end points within cable spatial vectors
                                ! for the current land-cell's tiles. These are just 
                                ! aliases to improve code readability
INTEGER(i4b) :: iland         ! loop counter through mland land cells
INTEGER(i4b) :: param_unit    ! Unit number for reading (all) parameter files.
INTEGER(i4b) :: error_status  ! Error status returned by OPENs
TYPE(climate_type), INTENT(INOUT)       :: climate ! climate variables

REAL(sp), ALLOCATABLE :: vegtypeigbp(:)
REAL(sp), ALLOCATABLE :: avgannmax_fapar(:)


ALLOCATE (vegtypeigbp(mland),  avgannmax_fapar(mland))

CALL GET_UNIT(param_unit)  ! Obtain an unused unit number for file reading, reused for all soil vars.

OPEN(param_unit, FILE=TRIM(param_path)//TRIM(vegtypeigbp_file), ACCESS='STREAM', &
     FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(vegtypeigbp_file) ; STOP -1
ELSE
  READ (param_unit) vegtypeigbp
  CLOSE (param_unit)
END IF



OPEN (param_unit, FILE=TRIM(param_path)//TRIM(avgannmax_fapar_file), ACCESS='STREAM', &
     FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(avgannmax_fapar_file) ; STOP -1
ELSE
  READ (param_unit) avgannmax_fapar
  CLOSE (param_unit)
END IF

DO iland = 1,mland ! For each land cell...
   is = landpt(iland)%cstart  ! Index position for the first tile of this land cell.
   ie = landpt(iland)%cend    ! Index position for the last tile of this land cell.
   climate%modis_igbp(is:ie) = INT(vegtypeigbp(iland))
   !climate%AvgAnnRainf(is:ie) = MAP(iland)
   climate%AvgAnnMaxFAPAR(is:ie) = avgannmax_fapar(iland)
ENDDO


END SUBROUTINE cable_bios_load_climate_params
!******************************************************************************
  
END MODULE cable_bios_met_obs_params
    
!******************************************************************************
!******************************************************************************

  
    
    
    
    
