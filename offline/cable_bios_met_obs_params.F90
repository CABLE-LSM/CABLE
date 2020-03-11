
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
  IF (M<3.0_8) THEN 
    M = M + 12.0_8 
    Y = Y - 1.0_8 
  END IF 
  JulianDay = D + int((153.0_8 * M - 457.0_8)/5.0_8) + (365.0_8*Y) +     &
              FLOOR(Y/4.0_8) - FLOOR(Y/100.0_8) + FLOOR(Y/400.0_8) + 1721118.5_8
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
  Z = FLOOR(JulianDay - 1721118.5_8) 
  R = JulianDay - 1721118.5_8 - Z
  G = Z - 0.25_8 
  A = FLOOR(G/36524.25_8) 
  B = A - FLOOR(A/4.0_8) 
  Y = FLOOR((B+G)/365.25_8) 
  C = B + Z - FLOOR(365.25_8*Y) 
  M = int(((5.0_8*C) + 456.0_8) / 153.0_8) 
  D = C - int((153.0_8 * M - 457.0_8) / 5.0_8) + R  
  IF (M > 12.0_8) THEN 
    Y = Y + 1.0_8 
    M = M - 12.0_8 
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
  REAL(sp)       :: MaskCtrW, MaskCtrS , tmp
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

  metgrid = 'land'

!%%%%%%%%%%%%%%%

!print *,Run, met_path, param_path, landmaskflt_file, landmaskhdr_file, &
!    rain_file, swdown_file, tairmax_file, tairmin_file, co2_file, &
!    b1_file, b2_file, bulkdens1_kgm3_file, bulkdens2_kgm3_file, clayfrac1_file, clayfrac2_file, &
!    csoil1_file, csoil2_file, depth1_m_file, depth2_m_file, hyk1sat_ms_file, hyk2sat_ms_file, & 
!    psie1_m_file, psie2_m_file, siltfrac1_file, siltfrac2_file, wvol1fc_m3m3_file, wvol2fc_m3m3_file, &
!    wvol1sat_m3m3_file, wvol2sat_m3m3_file, wvol1w_m3m3_file, wvol2w_m3m3_file

!stop 'end test'

  ! Assign MetForcing, CO2 and Ndep cases based only on the value of Run from bios namelist file
  SELECT CASE (TRIM(Run))
  CASE( "spinup" )       ! corresponds to S0_TRENDY in CRU
    MetForcing = "recycled"
    CO2     = "preind"
    Ndep    = "preind"
    WRITE(*   ,*)"Run = 'spinup': Therefore MetForcing = 'recycled', CO2 = 'preind', Ndep = 'preind'"
    WRITE(logn,*)"Run = 'spinup': Therefore MetForcing = 'recycled', CO2 = 'preind', Ndep = 'preind'"
  CASE( "premet" )    ! corresponds to S1_TRENDY in CRU
    MetForcing = "recycled"
    CO2     = "actual"
    Ndep     = "actual"
    WRITE(*   ,*)"Run = 'premet': Therefore MetForcing = 'recycled', CO2 = 'actual', Ndep = 'actual'"
    WRITE(logn,*)"Run = 'premet': Therefore MetForcing = 'recycled', CO2 = 'actual', Ndep = 'actual'"
  CASE( "standard" )    ! corresponds to S2_TRENDY in CRU
    MetForcing = "actual"
    CO2     = "actual"
    Ndep     = "actual"
    WRITE(*   ,*)"Run = 'standard': Therefore MetForcing = 'actual', CO2 = 'actual', Ndep = 'actual'"
    WRITE(logn,*)"Run = 'standard': Therefore MetForcing = 'actual', CO2 = 'actual', Ndep = 'actual'"
  CASE  default
    WRITE(*   ,*)"Wrong Run in bios nml file: ",Run
    WRITE(*   ,*)"Use: spinup, premet or standard!"
    WRITE(*   ,*)"Program stopped"
    WRITE(logn,*)"Wrong Run in bios nml file: ",Run
    WRITE(logn,*)"Use: spinup, premet or standard!"
    WRITE(logn,*)"Program stopped"
    STOP  
  END SELECT

  ! Print settings
  WRITE(*   ,*)"========================================= BIOS ============"
  WRITE(*   ,*)"BIOS settings chosen:"
  WRITE(*   ,*)" Run      = ",TRIM(Run)
  WRITE(*   ,*)" MetForcing (assigned) : ",TRIM(MetForcing)
  WRITE(*   ,*)" CO2     (assigned) : ",TRIM(CO2)
  WRITE(*   ,*)" Ndep     (assigned) : ",TRIM(Ndep)
  WRITE(*   ,*)" met_path      = ",TRIM(met_path)
  WRITE(*   ,*)" param_path    = ",TRIM(param_path)
  WRITE(*   ,*)" landmaskflt_file  = ",TRIM(landmaskflt_file)
  WRITE(*   ,*)" landmaskhdr_file  = ",TRIM(landmaskhdr_file)
  WRITE(*   ,*)" rain_file     = ",TRIM(rain_file)
  WRITE(*   ,*)" swdown_file   = ",TRIM(swdown_file)
  WRITE(*   ,*)" tairmax_file  = ",TRIM(tairmax_file)
  WRITE(*   ,*)" tairmin_file  = ",TRIM(tairmin_file)
  WRITE(*   ,*)" wind_file  = ",TRIM(wind_file)
  WRITE(*   ,*)" vp0900_file  = ",TRIM(vp0900_file)
  WRITE(*   ,*)" vp1500_file  = ",TRIM(vp1500_file)
  WRITE(*   ,*)" co2_file      = ",TRIM(co2_file)
  WRITE(*   ,*)" b1_file             = ",TRIM(b1_file)
  WRITE(*   ,*)" b2_file             = ",TRIM(b2_file)
  WRITE(*   ,*)" bulkdens1_kgm3_file = ",TRIM(bulkdens1_kgm3_file)
  WRITE(*   ,*)" bulkdens2_kgm3_file = ",TRIM(bulkdens2_kgm3_file)
  WRITE(*   ,*)" clayfrac1_file      = ",TRIM(clayfrac1_file)
  WRITE(*   ,*)" clayfrac2_file      = ",TRIM(clayfrac2_file)
  WRITE(*   ,*)" csoil1_file         = ",TRIM(csoil1_file)
  WRITE(*   ,*)" csoil2_file         = ",TRIM(csoil2_file)
  WRITE(*   ,*)" depth1_m_file       = ",TRIM(depth1_m_file)
  WRITE(*   ,*)" depth2_m_file       = ",TRIM(depth2_m_file)
  WRITE(*   ,*)" hyk1sat_ms_file     = ",TRIM(hyk1sat_ms_file)
  WRITE(*   ,*)" hyk2sat_ms_file     = ",TRIM(hyk2sat_ms_file)
  WRITE(*   ,*)" psie1_m_file        = ",TRIM(psie1_m_file)
  WRITE(*   ,*)" psie2_m_file        = ",TRIM(psie2_m_file)
  WRITE(*   ,*)" siltfrac1_file      = ",TRIM(siltfrac1_file)
  WRITE(*   ,*)" siltfrac2_file      = ",TRIM(siltfrac2_file)
  WRITE(*   ,*)" wvol1fc_m3m3_file   = ",TRIM(wvol1fc_m3m3_file)
  WRITE(*   ,*)" wvol2fc_m3m3_file   = ",TRIM(wvol2fc_m3m3_file)
  WRITE(*   ,*)" wvol1sat_m3m3_file  = ",TRIM(wvol1sat_m3m3_file)
  WRITE(*   ,*)" wvol2sat_m3m3_file  = ",TRIM(wvol2sat_m3m3_file)
  WRITE(*   ,*)" wvol1w_m3m3_file    = ",TRIM(wvol1w_m3m3_file)
  WRITE(*   ,*)" wvol2w_m3m3_file    = ",TRIM(wvol2w_m3m3_file)
  WRITE(*   ,*)" MVG_file    = ",TRIM(MVG_file)
  WRITE(*   ,*)" c4frac_file    = ",TRIM(c4frac_file)
  WRITE(*   ,*)" vegtypeigbp_file    = ",TRIM(vegtypeigbp_file)
  WRITE(*   ,*)" avgannmax_fapar_file   = ",TRIM(avgannmax_fapar_file)
  !WRITE(*   ,*)" slope_deg_file      = ",TRIM(slope_deg_file)
  WRITE(*   ,*)" DT(secs): ",dels

  WRITE(logn,*)"========================================= BIOS ============"
  WRITE(logn,*)"BIOS settings chosen:"
  WRITE(logn,*)" Run      = ",TRIM(Run)
  WRITE(logn,*)" MetForcing (assigned) : ",TRIM(MetForcing)
  WRITE(logn,*)" CO2     (assigned) : ",TRIM(CO2)
  WRITE(logn,*)" Ndep     (assigned) : ",TRIM(Ndep)
  WRITE(logn,*)" met_path      = ",TRIM(met_path)
  WRITE(logn,*)" param_path    = ",TRIM(param_path)
  WRITE(logn,*)" landmaskflt_file  = ",TRIM(landmaskflt_file)
  WRITE(logn,*)" landmaskhdr_file  = ",TRIM(landmaskhdr_file)
  WRITE(logn,*)" rain_file     = ",TRIM(rain_file)
  WRITE(logn,*)" swdown_file   = ",TRIM(swdown_file)
  WRITE(logn,*)" tairmax_file  = ",TRIM(tairmax_file)
  WRITE(logn,*)" tairmin_file  = ",TRIM(tairmin_file)
  WRITE(logn,*)" co2_file      = ",TRIM(co2_file)
  WRITE(logn,*)" b1_file             = ",TRIM(b1_file)
  WRITE(logn,*)" b2_file             = ",TRIM(b2_file)
  WRITE(logn,*)" bulkdens1_kgm3_file = ",TRIM(bulkdens1_kgm3_file)
  WRITE(logn,*)" bulkdens2_kgm3_file = ",TRIM(bulkdens2_kgm3_file)
  WRITE(logn,*)" clayfrac1_file      = ",TRIM(clayfrac1_file)
  WRITE(logn,*)" clayfrac2_file      = ",TRIM(clayfrac2_file)
  WRITE(logn,*)" csoil1_file         = ",TRIM(csoil1_file)
  WRITE(logn,*)" csoil2_file         = ",TRIM(csoil2_file)
  WRITE(logn,*)" depth1_m_file       = ",TRIM(depth1_m_file)
  WRITE(logn,*)" depth2_m_file       = ",TRIM(depth2_m_file)
  WRITE(logn,*)" hyk1sat_ms_file     = ",TRIM(hyk1sat_ms_file)
  WRITE(logn,*)" hyk2sat_ms_file     = ",TRIM(hyk2sat_ms_file)
  WRITE(logn,*)" psie1_m_file        = ",TRIM(psie1_m_file)
  WRITE(logn,*)" psie2_m_file        = ",TRIM(psie2_m_file)
  WRITE(logn,*)" siltfrac1_file      = ",TRIM(siltfrac1_file)
  WRITE(logn,*)" siltfrac2_file      = ",TRIM(siltfrac2_file)
  WRITE(logn,*)" wvol1fc_m3m3_file   = ",TRIM(wvol1fc_m3m3_file)
  WRITE(logn,*)" wvol2fc_m3m3_file   = ",TRIM(wvol2fc_m3m3_file)
  WRITE(logn,*)" wvol1sat_m3m3_file  = ",TRIM(wvol1sat_m3m3_file)
  WRITE(logn,*)" wvol2sat_m3m3_file  = ",TRIM(wvol2sat_m3m3_file)
  WRITE(logn,*)" wvol1w_m3m3_file    = ",TRIM(wvol1w_m3m3_file)
  WRITE(logn,*)" wvol2w_m3m3_file    = ",TRIM(wvol2w_m3m3_file)
  WRITE(logn,*)" MVG_file    = ",TRIM(MVG_file)
  WRITE(logn ,*)" c4frac_file    = ",TRIM(c4frac_file)
  WRITE(logn,*)" vegtypeigbp_file    = ",TRIM(vegtypeigbp_file)
  WRITE(logn,*)" avgannmax_fapar_file   = ",TRIM(avgannmax_fapar_file)
  !WRITE(logn,*)" slope_deg_file      = ",TRIM(slope_deg_file)
  WRITE(logn,*)" timestep in secs  = ",dels
  WRITE(logn,*)" DT(secs): ",dels

  WRITE(*   ,*)"========================================= BIOS ============"
  WRITE(logn,*)"========================================= BIOS ============"
  
  ! Read the header file for the landmask file and parse it for dimensions
  ! and no-data value. Allocate logical, integer (CABLE) and real (temporary
  ! arrays for land masks.
  
  CALL GET_UNIT(iunit)
  CALL ReadArcFltHeader(iunit,landmaskhdr_file,MaskCols,MaskRows,MaskBndW, & 
       MaskBndS,MaskRes,NoDataVal)
  CLOSE (iunit) 
 
  xdimsize = MaskCols
  ydimsize = MaskRows
  write(*,*) ' MaskCols,MaskRows', MaskCols, MaskRows 
  ALLOCATE (LandMaskLogical(MaskCols,MaskRows))
  ALLOCATE (mask           (MaskCols,MaskRows))   ! mask: CABLE's integer land mask
  ALLOCATE (LandMaskReal   (MaskCols,MaskRows))
  
  ! Read the real land mask.
  
  CALL GET_UNIT(iunit)
  OPEN (iunit,file=landmaskflt_file,access='stream',form='unformatted',status='old')
  write(*,*) 'after landmask_flt_file open'
  !READ (iunit) tmp
  !READ (iunit) tmp
  !write(*,*) tmp
  REWIND (iunit) 
  READ (iunit) LandMaskReal
  write(*,*) 'after landmask_flt_file read'
  CLOSE (iunit)
  
  ! Initialise CABLE (mask) and logical landmasks, and assign them from the  
  ! real land mask. 
  
  mask = 0
  LandMaskLogical = .false.
  WHERE (NINT(LandMaskReal) > NINT(NoDataVal)) ! Where land...
      mask = 1
      LandMaskLogical = .true.
  END WHERE
  DEALLOCATE (LandMaskReal)
    
  mland = COUNT(LandMaskLogical)         ! mland: CABLE's count of land points
  nmetpatches = 1
  
  ! Allocate vectors for supplying the lat/long and col/row pairs of each
  ! land cell to cable.
  ALLOCATE( latitude(mland), longitude(mland) )
  ALLOCATE( land_y  (mland), land_x   (mland) )
  
  ! Populate a temporary integer grid with column numbers, then row numbers,
  ! packing them with the landmask into the cable vectors that record the 
  ! column and row numbers of the land cells.
  ALLOCATE (ColRowGrid(MaskCols,MaskRows))
  FORALL (icol=1:MaskCols) ColRowGrid(icol,:) = icol
  land_x = PACK(ColRowGrid,mask=LandMaskLogical)
  FORALL (irow=1:MaskRows) ColRowGrid(:,irow) = irow
  land_y = PACK(ColRowGrid,mask=LandMaskLogical)

! convert latitude and longitude to array
  ALLOCATE (lat_all(MaskCols,MaskRows))
  ALLOCATE (lon_all(MaskCols,MaskRows))
  
  
!
  
! Using the landmask grid boundaries and dimensions, translate the land  
! cell cols and rows of land_x and land_y into corresponding lats and longs.
  MaskCtrW = MaskBndW + (MaskRes / 2.0) ! Convert western and southern 
  MaskCtrS = MaskBndS + (MaskRes / 2.0) ! boundaries to cell centres.
  DO iLand = 1,mland
    longitude(iLand) = MaskRes * real((land_x(iLand) - 1)) + MaskCtrW
    latitude(iLand) = MaskCtrS + (real(MaskRows - land_y(iLand)) * MaskRes)
 END DO

  lat_all = UNPACK(latitude,mask=LandMaskLogical,field=-9999.)
  lon_all = UNPACK(longitude,mask=LandMaskLogical,field=-9999.)

  !Finished reading grids. Only mland vectors from now on.
  DEALLOCATE (ColRowGrid)
  DEALLOCATE (LandMaskLogical)

! If this is the first ever initialisation, open the met files, sort
! out the run's date range with respect to the available met data,
! and allocate memory for met vars. Otherwise, just rewind the  
! already-open met files.
  IF (call1) then
      
    CALL open_bios_met
    ALLOCATE (   rain_day(mland))
    ALLOCATE ( swdown_day(mland))
    ALLOCATE (tairmax_day(mland))
    ALLOCATE (tairmin_day(mland))
    ALLOCATE (prev_tairmax_day(mland))
    ALLOCATE (next_tairmin_day(mland))
    ALLOCATE (wind_day(mland))
    ALLOCATE (vp0900(mland))
    ALLOCATE (vp1500(mland))
    ALLOCATE (prev_vp1500(mland))
    ALLOCATE (next_vp0900(mland))

! Whether start and end dates are user specified or not, we need to know 
! what the date range in the met files is, so dummy-read the rainfall.
    READ (rain_unit,IOSTAT=error_status) bios_startdate, rain_day
    DO WHILE (error_status > 0)
       READ (rain_unit,IOSTAT=error_status) bios_enddate, rain_day
    END DO
    REWIND (rain_unit)
    
! If the user-specified dates are zero, this means we should define the 
! run length by the first and last dates in the met file. 
    IF (cable_user%yearstart == 0 .and. cable_user%yearend == 0) THEN

! Using date arithmetic, calculate the number of days in the bios run. 
      bios_rundays =  DayDifference(bios_enddate, bios_startdate) + 1
      skipdays = 0

! Initialise the previous_date as one day before the bios_startdate,
! using date arithmetic
      previous_date = bios_startdate - 1
      
! Set the CABLE timing vars
      curyear = bios_startdate%year  ! Current CABLE year set to BIOS start year
      kend = ktauday * bios_rundays  ! Run endpoint is the number of timesteps
                                     ! per day (ktauday) * num of days in run.

      shod        = 0.
      sdoy        = 1
      smoy        = 1
      syear       = curyear

!      met%hod (landpt(:)%cstart) = 0 ! Cable run always starts at midnight
!      met%doy (landpt(:)%cstart) = 1 ! Cable run always starts on Jan 1...
!      met%year(landpt(:)%cstart) = curyear ! ...of the current year

      cable_user%yearstart = bios_startdate%year
      cable_user%yearend = bios_enddate%year
      
    ELSE
! Dates are user specified, so they must be tested for compliance with the 
! available met, and the CABLE timing parameters estimated.
        
! User-specified start & end dates are just years: convert them to dmy dates.
      user_startdate = dmydate(1,1,cable_user%yearstart)
      user_enddate   = dmydate(31,12,cable_user%yearend)
      curyear = user_startdate%year  ! The initial current year is set prior to any skipping.
      
! Initialise the previous_date as one day before the user_startdate,
! using date arithmetic.
      previous_date = user_startdate - 1
      shod        = 0.
      sdoy        = 1
      smoy        = 1
      !syear       = 1690
      syear = 1981
      write(*,*) 'prev:',previous_date%year,previous_date%month,previous_date%day
      write(*,*) 'run:',  user_startdate%year, user_startdate%month,  user_startdate%day      
      ! For spinup and initialisation before bios met begins (1900),
      ! calculate the required met year for repeatedly cycling through the
      ! spinup meteorology between recycle_met_startdate - recycle_met_enddate (1900-1929).
      ! For normal runs (e.g. 1900-2016), MetDate%Year = curyear.
      MetDate%day = 1
      MetDate%month = 1
      IF (TRIM(MetForcing) .EQ. 'recycled') THEN
         MetDate%Year = recycle_met_startdate + MOD(curyear-syear,recycle_met_enddate-recycle_met_startdate+1)
         write(*,*) 'metdatestart,: ',  MetDate%Year,  recycle_met_startdate,  curyear, syear, &
                       recycle_met_enddate, recycle_met_startdate,  MOD(curyear-syear,recycle_met_enddate-recycle_met_startdate+1)
      ELSE IF (TRIM(MetForcing) .EQ. 'actual' ) THEN
        MetDate%Year = curyear
      ENDIF

! The user-specified run startdate must be no earlier than the first available met.
! Stop if the user startdate is mis-specified in this way.
write(6,*) 'MetDate, bios_startdate=',MetDate, bios_startdate
      skipdays = DayDifference(MetDate, bios_startdate)
      IF (skipdays < 0) THEN
        STOP "ERROR in cable_bios_init: user_startdate < bios_startdate"
      END IF
      
      ! In this case the bios run is what is left of the bios met
      ! after skipping to the user startdate.
      bios_rundays = DayDifference(bios_enddate , user_startdate) + 1
      kend = bios_rundays * ktauday  ! Final timestep

!      met%hod (landpt(:)%cstart) = 0
!      met%doy (landpt(:)%cstart) = 1
!      met%year(landpt(:)%cstart) = curyear

    END IF  ! End of test for user start/end dates set to zero.
    
! CO2 is handled as a special case because it is annual. 
! Dummy-read the CO2 file to get the date range, of which
! only the year is relevant to preserve.
    READ (co2_unit,IOSTAT=error_status) dummydate, co2air_year
    co2_startyear = dummydate%year

    
    DO WHILE (error_status > 0)
       READ (co2_unit,IOSTAT=error_status) dummydate, co2air_year
       co2_endyear = dummydate%year
    END DO
    REWIND (co2_unit)

! From the curyear (the first year of relevant met), calculate how many
! years of annual co2 records to skip so as to position so that the next
! record to read will be the current year.
    co2_skipyears = curyear - co2_startyear
    
! Initialise Weather Generator
    CALL WGEN_INIT( WG, mland, latitude, dels )
    
! Record that the first call is now over.
    call1 = .FALSE.
     
  ELSE ! Not the first ever initialisation, so just rewind the met files.
    REWIND (rain_unit)
    REWIND (swdown_unit)
    REWIND (tairmax_unit)  
    REWIND (tairmin_unit)
     IF (TRIM(wind_file) .NE. 'none') THEN
        REWIND (wind_unit)
     ENDIF
     IF (TRIM(vp0900_file) .NE. 'none') THEN
        REWIND (vp0900_unit)
     ENDIF
     IF (TRIM(vp1500_file) .NE. 'none') THEN
        REWIND (vp1500_unit)
     ENDIF
    REWIND (co2_unit)
  END IF
  
  ! Regardless of whether this is the first or subsequent initialisation:
  ! If skipdays > 0, the user_startdate is after the start of the met files, so
  ! all the BIOS met files must be repositioned to that date by reading skipdays 
  ! daily met records. 
  DO iday = 1,skipdays
    READ (rain_unit) dummydate, rain_day
    READ (swdown_unit) dummydate, swdown_day
    READ (tairmax_unit) dummydate, tairmax_day
    READ (tairmin_unit) dummydate, tairmin_day
    IF (TRIM(wind_file) .NE. 'none') THEN
       READ (wind_unit) dummydate, wind_day          !
    ENDIF
    IF (TRIM(vp0900_file) .NE. 'none') THEN
       READ (vp0900_unit) dummydate, vp0900          !
    ENDIF
    IF (TRIM(vp1500_file) .NE. 'none') THEN     
       READ (vp1500_unit) dummydate, vp1500
    ENDIF!

  END DO
  dummydate = dummydate + 1
  IF (skipdays .gt. 0) THEN
    write(*,*) 'Met skipped, ready to read from ',dummydate
    write(logn,*) 'Met skipped, ready to read from ',dummydate
  END IF
  
  ! Likewise, skip through the annual CO2 records to position for reading
  ! the first relevant year.
  DO iyear = 1,co2_skipyears
     READ (co2_unit) dummydate, co2air_year
  END DO
  
  CONTAINS

    SUBROUTINE open_bios_met

    IMPLICIT NONE
    INTEGER :: error_status
  
  ! Open each met file for the first time, stopping if the file is not found.
  
    CALL GET_UNIT(rain_unit)  ! Rainfall
    OPEN (rain_unit, FILE=TRIM(met_path)//TRIM(rain_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
    IF (error_status > 0) THEN
      WRITE (*,*) "STOP - File not found: ", TRIM(met_path)//TRIM(rain_file)
      STOP ''
    END IF
  
    CALL GET_UNIT(swdown_unit)  ! Shortwave downward solar radiation
    OPEN (swdown_unit, FILE=TRIM(met_path)//TRIM(swdown_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
         IOSTAT=error_status)
    IF (error_status > 0) THEN
      WRITE (*,*) "STOP - File not found: ", TRIM(met_path)//TRIM(swdown_file)
      STOP ''
    END IF
  
    CALL GET_UNIT(tairmax_unit)  ! Maximum air temperature
    OPEN (tairmax_unit, FILE=TRIM(met_path)//TRIM(tairmax_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
         IOSTAT=error_status)
    IF (error_status > 0) THEN
      WRITE (*,*) "STOP - File not found: ", TRIM(met_path)//TRIM(tairmax_file)
      STOP ''
    END IF   
  
    CALL GET_UNIT(tairmin_unit)  ! Minimum air temperature  
    OPEN (tairmin_unit, FILE=TRIM(met_path)//TRIM(tairmin_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
         IOSTAT=error_status)
    IF (error_status > 0) THEN
      WRITE (*,*) "STOP - File not found: ", TRIM(met_path)//TRIM(tairmin_file)
      STOP ''
   END IF

   IF (TRIM(wind_file) .NE. 'none') THEN
      CALL GET_UNIT(wind_unit)  ! wind speed  
      OPEN (wind_unit, FILE=TRIM(met_path)//TRIM(wind_file), &
           ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
      IF (error_status > 0) THEN
         WRITE (*,*) "STOP - File not found: ", TRIM(met_path)//TRIM(wind_file)
         STOP ''
      END IF
   ENDIF
   
   IF (TRIM(vp0900_file) .NE. 'none') THEN
      CALL GET_UNIT(vp0900_unit)  ! vp 0900  
      OPEN (vp0900_unit, FILE=TRIM(met_path)//TRIM(vp0900_file), &
           ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
      IF (error_status > 0) THEN
         WRITE (*,*) "STOP - File not found: ", TRIM(met_path)//TRIM(vp0900_file)
         STOP ''
      END IF
   ENDIF
   
   IF (TRIM(vp1500_file) .NE. 'none') THEN
      CALL GET_UNIT(vp1500_unit)  ! vp 1500
      OPEN (vp1500_unit, FILE=TRIM(met_path)//TRIM(vp1500_file), &
           ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
      IF (error_status > 0) THEN
         WRITE (*,*) "STOP - File not found: ", TRIM(met_path)//TRIM(vp1500_file)
         STOP ''
      END IF
   ENDIF
  
    CALL GET_UNIT(co2_unit)  ! CO2
    OPEN (co2_unit, FILE=TRIM(met_path)//TRIM(co2_file),&
         ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
    IF (error_status > 0) THEN
      WRITE (*,*) "STOP - File not found: ", TRIM(met_path)//TRIM(co2_file)
      STOP ''
    END IF
  
    END SUBROUTINE open_bios_met 
  
  
  END SUBROUTINE cable_bios_init 

!******************************************************************************
  
  SUBROUTINE cable_bios_read_met(MET, CurYear, ktau, kend, islast, dels )
  
    ! Read a single day of meteorology from all bios met files, updating the bios_rundate
    ! If a change of year has occurred, read an annual CO2 record
    
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: CurYear, ktau, kend
    LOGICAL, INTENT(IN)  :: islast
    REAL, INTENT(IN) :: dels                        ! time step size in seconds
    TYPE(MET_TYPE), INTENT(INOUT)       :: MET

    LOGICAL(lgt)   :: newday
!    real(sp),parameter:: RMW       = 0.018016 ! molecular wt of water     [kg/mol]
!    real(sp),parameter:: RMA       = 0.02897 ! atomic wt of C            [kg/mol]
    real(sp),parameter:: RMWbyRMA  = 0.62188471 ! molecular wt of water [kg/mol] / atomic wt of C [kg/mol]
    integer(i4b)   :: iday
    integer(i4b)   :: iland       ! Loop counter through mland land cells
    integer(i4b)   :: is, ie      ! For each land cell, the start and ending index position within the larger cable spatial
                                  ! vectors of the first and last tile for that land cell. 

    met%hod (landpt(:)%cstart) = REAL(MOD( (ktau-1) * NINT(dels), INT(SecDay)) ) / 3600.
    met%doy (landpt(:)%cstart) = INT(REAL(ktau-1) * dels / SecDay ) + 1
    met%year(landpt(:)%cstart) = Curyear  

    newday = ( met%hod(landpt(1)%cstart).EQ. 0 )
    IF ( newday ) THEN
       ! get current day's met
       READ (rain_unit) bios_rundate, rain_day          ! Packed vector of daily AWAP/BIOS rain (mm) 

       READ (swdown_unit) bios_rundate, swdown_day        ! Packed vector of daily AWAP/BIOS swdown (MJ)
       READ (tairmax_unit) bios_rundate, tairmax_day       ! Packed vector of daily AWAP/BIOS max air temp (deg C)
       READ (tairmin_unit) bios_rundate, tairmin_day       ! Packed vector of daily AWAP/BIOS min air temp (deg C)
       IF (TRIM(wind_file) .NE. 'none') THEN
          READ (wind_unit) bios_rundate, wind_day          !
       ENDIF
       IF (TRIM(vp0900_file) .NE. 'none') THEN
          READ (vp0900_unit) bios_rundate, vp0900          !
       ENDIF
       IF (TRIM(vp1500_file) .NE. 'none') THEN     
          READ (vp1500_unit) bios_rundate, vp1500
       ENDIF

       IF (MetDate /= bios_rundate) THEN
         write(*,*) 'Expecting to read met for ',MetDate, &
              'but actually read met for ',bios_rundate
         write(*,*) 'Program stopped'
         STOP
       END IF

       ! Increment MetDate. 
       MetDate = MetDate + 1
  
       ! If MetForcing is 'recycled', check whether we need to rewind files and skip to correct position
       if (MetForcing .eq. 'recycled') then
         if (MetDate%Year .gt. recycle_met_enddate) then
           MetDate%Day = 1
           MetDate%Month = 1
           MetDate%Year = recycle_met_startdate
           skipdays = DayDifference(MetDate, bios_startdate)
           write(*,*) 'skipdays', skipdays
           REWIND (rain_unit)
           REWIND (swdown_unit)
           REWIND (tairmax_unit)
           REWIND (tairmin_unit)
           IF (TRIM(wind_file) .NE. 'none') THEN
              REWIND (wind_unit)
           ENDIF
           IF (TRIM(vp0900_file) .NE. 'none') THEN
              REWIND (vp0900_unit)
           ENDIF
           IF (TRIM(vp1500_file) .NE. 'none') THEN       
              REWIND (vp1500_unit)
           ENDIF
           !REWIND (co2_unit)
           DO iday = 1,skipdays
             READ (rain_unit) dummydate, rain_day
             READ (swdown_unit) dummydate, swdown_day
             READ (tairmax_unit) dummydate, tairmax_day
             READ (tairmin_unit) dummydate, tairmin_day
             
             IF (TRIM(wind_file) .NE. 'none') THEN
                READ (wind_unit) dummydate, wind_day
             ENDIF
             IF (TRIM(vp0900_file) .NE. 'none') THEN           
                READ (vp0900_unit) dummydate, vp0900
             ENDIF
             IF (TRIM(vp1500_file) .NE. 'none') THEN   
                READ (vp1500_unit) dummydate, vp1500
             ENDIF
          END DO
           dummydate = dummydate + 1
           IF (skipdays .gt. 0) THEN
              write(*,*) 'Met skipped, ready to read from ',dummydate
            
             write(logn,*) 'Met skipped, ready to read from ',dummydate
           END IF
         endif
       endif

       !! NB !! Peter to fix as in Cabledyn !!
       !if (ktau.eq.1) then       
       !    WG%TempMaxDayPrev =  tairmax_day
       !else
          WG%TempMaxDayPrev = WG%TempMaxDay
       !endif

       !if (ktau.ne.kend) then
       !   READ (tairmin_unit) bios_rundate, next_tairmin_day       ! Packed vector of daily AWAP/BIOS min air temp (deg C)
       !   BACKSPACE(tairmin_unit)
       !endif

       next_tairmin_day =   tairmin_day
       prev_vp1500 = vp1500
       next_vp0900 = vp0900

        IF (TRIM(wind_file) .NE. 'none') THEN
           WG%WindDay        = wind_day
        ELSE
           WG%WindDay        = 3.0
        ENDIF
       WG%TempMinDay     = tairmin_day  
       WG%TempMaxDay     = tairmax_day  
       WG%TempMinDayNext =  next_tairmin_day

       WG%VapPmbDay = esatf(tairmin_day)

       IF (TRIM(vp0900_file) .NE. 'none') THEN
          WG%VapPmb0900 = vp0900
          WG%VapPmb1500 = vp1500
          WG%VapPmb1500Prev = prev_vp1500
          WG%VapPmb0900Next = next_vp0900
       
       ELSE
          WG%VapPmb0900 =  WG%VapPmbDay 
          WG%VapPmb1500 =   WG%VapPmbDay 
          WG%VapPmb1500Prev =  WG%VapPmbDay 
          WG%VapPmb0900Next =  WG%VapPmbDay
       ENDIF

       if (swdown_file(1:4) .eq. 'rsds') then
          WG%SolarMJDay     = swdown_day *3600.*24. *1.0e-6 ! Wm-2 -> MJm-2/Day
       else
          WG%SolarMJDay     = swdown_day
       endif
       WG%PrecipDay      = rain_day / 1000. ! ->[m/d]
       where ( WG%TempMinDay.lt.-2.0) 
          WG%SnowDay        =  rain_day / 1000. ! ->[m/d]
          WG%PrecipDay      = 0.0
       elsewhere
          WG%SnowDay        = 0.0
       end where
  
       WG%PmbDay = 1000.0 ! Air pressure in mb fixed in space and time

       IF (previous_date%year .ne. bios_rundate%year) THEN
         IF ( TRIM(CO2) .EQ. "preind") THEN
           met%ca(:) = 286.42 / 1.e+6  ! CO2 for 1860
         ELSE
           READ (co2_unit) dummydate, co2air_year          ! Single global value of co2 (ppm)
           met%ca(:) = CO2air_year / 1.e+6   
        END IF
         ! write(596,*), previous_date%year, bios_rundate%year, dummydate,  CO2air_year,  met%ca(1)*1.0e6
      END IF

      ! Finished operations for this day, so the current date becomes the previous date.  
       previous_date = bios_rundate

       CALL WGEN_DAILY_CONSTANTS( WG, mland, INT(met%doy(1))+1 )
  
    END IF

 CALL WGEN_SUBDIURNAL_MET( WG, mland, NINT(met%hod(1)*3600./dels) )
 
 DO iland = 1,mland ! For each land cell...
   is = landpt(iland)%cstart  ! Index position for the first tile of this land cell.
   ie = landpt(iland)%cend    ! Index position for the last tile of this land cell.

   met%precip(is:ie)    = WG%Precip(iland) 
   met%precip_sn(is:ie) = WG%Snow(iland)
   ! combine Rainf and Snowf inputs
   met%precip(is:ie)    = met%precip(is:ie) + met%precip_sn(is:ie)
   met%fld(is:ie)       = WG%PhiLD(iland)
   met%fsd(is:ie,1)     = WG%PhiSD(iland) * 0.5
   met%fsd(is:ie,2)     = WG%PhiSD(iland) * 0.5
   met%tk(is:ie)        = WG%Temp(iland) + 273.15
   met%ua(is:ie)        = WG%Wind(iland) *2.0  ! factor of 2 to convert from 2m screen height to 40 m zref, assuming logarithmic profile
   met%coszen(is:ie)    = WG%coszen(iland)
   met%qv(is:ie)        = WG%VapPmb(iland)/WG%Pmb(iland)*RMWbyRMA ! specific humidity (kg/kg)
   met%pmb(is:ie)       = WG%Pmb(iland)
   met%rhum(is:ie)  =  WG%VapPmb(iland)/esatf(real(WG%Temp(iland),sp)) *100.0 ! rel humidity (%)
   met%u10(is:ie) = met%ua(is:ie) 
   ! initialise within canopy air temp
   met%tvair(is:ie)     = met%tk(is:ie) 
   met%tvrad(is:ie)     = met%tk(is:ie)
end do

  

!*******************************************************************************
CONTAINS
  ELEMENTAL FUNCTION Esatf(TC)
!-------------------------------------------------------------------------------
! At temperature TC [deg C], return saturation water vapour pressure Esatf [mb] 
! from Teten formula.
! MRR, xx/1987
! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
!                 just like intrinsic functions.
! MRR, 12-mar-02: Convert Qsatf (specific humidity routine) to Esatf
!-------------------------------------------------------------------------------
  implicit none
  real(sp), intent(in):: TC           ! temp [deg C]
  real(sp):: Esatf                    ! saturation vapour pressure [mb]
  real(sp):: TCtmp                    ! local
  real(sp),parameter:: A = 6.106      ! Teten coefficients
  real(sp),parameter:: B = 17.27      ! Teten coefficients
  real(sp),parameter:: C = 237.3      ! Teten coefficients
!-------------------------------------------------------------------------------
  TCtmp = TC                          ! preserve TC
  if (TCtmp.gt.100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
  if (TCtmp.lt.-40.0) TCtmp = -40.0
  Esatf = A*EXP(B*TCtmp/(C+TCtmp))    ! sat vapour pressure (mb)
 
  END FUNCTION Esatf

!*******************************************************************************
 
END SUBROUTINE cable_bios_read_met

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
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(b1_file) ; STOP ''
ELSE
  READ (param_unit) b1
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(b2_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(b2_file) ; STOP ''
ELSE
  READ (param_unit) b2
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(bulkdens1_kgm3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(bulkdens1_kgm3_file) ; STOP ''
ELSE
  READ (param_unit) bulkdens1_kgm3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(bulkdens2_kgm3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(bulkdens2_kgm3_file) ; STOP ''
ELSE
  READ (param_unit) bulkdens2_kgm3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(clayfrac1_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(clayfrac1_file) ; STOP ''
ELSE
  READ (param_unit) clayfrac1
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(clayfrac2_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(clayfrac2_file) ; STOP ''
ELSE
  READ (param_unit) clayfrac2
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(csoil1_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(csoil1_file) ; STOP ''
ELSE
  READ (param_unit) csoil1
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(csoil2_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(csoil2_file) ; STOP ''
ELSE
  READ (param_unit) csoil2
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(depth1_m_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(depth1_m_file) ; STOP ''
ELSE
  READ (param_unit) depth1_m
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(depth2_m_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(depth2_m_file) ; STOP ''
ELSE
  READ (param_unit) depth2_m
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(hyk1sat_ms_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(hyk1sat_ms_file) ; STOP ''
ELSE
  READ (param_unit) hyk1sat_ms
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(hyk2sat_ms_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(hyk2sat_ms_file) ; STOP ''
ELSE
  READ (param_unit) hyk2sat_ms
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(psie1_m_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(psie1_m_file); STOP ''
ELSE
  READ (param_unit) psie1_m
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(psie2_m_file), &
     ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(psie2_m_file); STOP ''
ELSE
  READ (param_unit) psie2_m
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(siltfrac1_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(siltfrac1_file); STOP ''
ELSE
  READ (param_unit) siltfrac1
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(siltfrac2_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(siltfrac2_file) ; STOP ''
ELSE
  READ (param_unit) siltfrac2
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol1fc_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol1fc_m3m3_file); STOP ''
ELSE
  READ (param_unit) wvol1fc_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol2fc_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)   
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol2fc_m3m3_file); STOP ''
ELSE
  READ (param_unit) wvol2fc_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol1sat_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol1sat_m3m3_file) ; STOP ''
ELSE
  READ (param_unit) wvol1sat_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol2sat_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol2sat_m3m3_file) ; STOP ''
ELSE
  READ (param_unit) wvol2sat_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol1w_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol1w_m3m3_file) ; STOP ''
ELSE
  READ (param_unit) wvol1w_m3m3
  CLOSE (param_unit)
END IF

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(wvol2w_m3m3_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD', &
     IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(wvol2w_m3m3_file) ; STOP ''
ELSE
  READ (param_unit) wvol2w_m3m3
  CLOSE (param_unit)
END IF

!OPEN (param_unit, FILE=TRIM(param_path)//TRIM(slope_deg_file), ACCESS='STREAM',FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
!IF (error_status > 0) THEN
!  WRITE (*,'("STOP - File not found: ",  LEN_TRIM(TRIM(param_path)//TRIM(slope_deg_file)))') ; STOP ''
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


INTEGER(i4b) :: is, ie ! Index start/end points within cable spatial vectors
                                ! for the current land-cell's tiles. These are just 
                                ! aliases to improve code readability
INTEGER(i4b) :: iland         ! loop counter through mland land cells
INTEGER(i4b) :: param_unit    ! Unit number for reading (all) parameter files.
INTEGER(i4b) :: error_status  ! Error status returned by OPENs
REAL(sp), ALLOCATABLE :: tmp(:)
INTEGER, INTENT(INOUT)       :: MVG(:) ! climate variables


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
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(MVG_file) ; STOP ''
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


INTEGER(i4b) :: is, ie ! Index start/end points within cable spatial vectors
                                ! for the current land-cell's tiles. These are just 
                                ! aliases to improve code readability
INTEGER(i4b) :: iland         ! loop counter through mland land cells
INTEGER(i4b) :: param_unit    ! Unit number for reading (all) parameter files.
INTEGER(i4b) :: error_status  ! Error status returned by OPENs
REAL(sp), ALLOCATABLE :: tmp(:)
REAL, INTENT(INOUT)       :: fracC4(:) ! climate variables


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
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(c4frac_file) ; STOP ''
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

OPEN (param_unit, FILE=TRIM(param_path)//TRIM(vegtypeigbp_file), ACCESS='STREAM', &
     FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
print*, TRIM(param_path)//TRIM(vegtypeigbp_file)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(vegtypeigbp_file) ; STOP ''
ELSE
  READ (param_unit) vegtypeigbp
  CLOSE (param_unit)
END IF



OPEN (param_unit, FILE=TRIM(param_path)//TRIM(avgannmax_fapar_file), ACCESS='STREAM', &
     FORM='UNFORMATTED', STATUS='OLD',IOSTAT=error_status)
IF (error_status > 0) THEN
  WRITE (*,*) "STOP - File not found: ", TRIM(param_path)//TRIM(avgannmax_fapar_file) ; STOP ''
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

  
    
    
    
    
