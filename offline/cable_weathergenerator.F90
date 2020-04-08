MODULE CABLE_WEATHERGENERATOR
  
  USE CABLE_COMMON_MODULE, ONLY:  cable_user
  
  IMPLICIT NONE

  ! Parameters
  INTEGER, PARAMETER,PRIVATE :: sp = 8 ! this is double precision
  REAL(sp),PARAMETER,PRIVATE :: &
       Pi         = 3.14159265_sp , &            ! Pi
       PiBy2      = 1.57079632_sp , &            ! Pi/2
       SecDay     = 86400.0_sp    , &            ! Seconds/day
       SolarConst = 1370._sp*SecDay/1e6_sp, &    ! Solar constant [MJ/m2/day]
       epsilon    = 0.736_sp        , &
       SBoltz     = 5.67e-8_sp  ! Stefan-Boltzmann constant [W/m2/K4]
  ! Global variables
  TYPE WEATHER_GENERATOR_TYPE
     ! general
     INTEGER :: np
     INTEGER :: ndtime   ! number of subdiurnal steps in 24 hr
     REAL    :: delT     ! subdiurnal timestep (in s)
     REAL(sp),DIMENSION(:),ALLOCATABLE :: LatDeg

     ! Daily input to subdiurnal met calculations
     REAL(sp),DIMENSION(:),ALLOCATABLE :: &
          WindDay        ,&
          TempMinDay     ,&    ! Daily minimum air temp [degC]
          TempMaxDay     ,&    ! Daily maximum air temp [degC]
          TempMinDayNext ,&    ! Daily minimum air temp tomorrow [degC]
          TempMaxDayPrev ,&    ! Daily maximum air temp yesterday [degC]
          SolarMJDay     ,&      ! Total daily down-ward irradiation [MJ/m2/d] (swdown)
          PrecipDay      ,&  ! 24hr-total precipitation [m/day] (current day = 24h from 0900 on previous day)
          SnowDay        ,&  ! 24hr-total precipitation [m/day] (current day = 24h from 0900 on previous day)
          PmbDay         ,&    ! 24hr-av pressure [mb]
          VapPmbDay      ,&           ! 24hr-av water vapour pressure [mb]
!CLNx     PrecipDayNext   ,&  ! 24hr-total precipitation [m/day] (next day = 24 h from 0900 on current day)
          VapPmb0900      ,&  ! 0900 water vapour pressure [mb]
          VapPmb1500      ,&  ! 1500 water vapour pressure [mb]
          VapPmb1500Prev  ,&  ! 1500 (prev day) water vapour pressure [mb]
          VapPmb0900Next    ! 0900(next day) water vapour pressure [mb]

     ! Daily constants
     REAL(sp) :: DecRad                                   ! Declination in radians
     REAL(sp),DIMENSION(:),ALLOCATABLE :: &
          WindDark           ,&  ! Wind [m/s] during night hrs
          WindLite           ,&  ! Wind [m/s] during daylight hrs
          SolarNorm          ,&  ! (daily solar)/(solar const): geometry
          LatRad             ,&  ! Latitude in radians
          DayLength          ,&
          TimeSunsetPrev     ,&
          TimeSunrise        ,&
          TimeMaxTemp        ,&
          TimeSunset         ,&
          TempSunsetPrev     ,&
          TempSunset         ,&
          TempNightRate      ,&
          TempNightRatePrev  ,&
          TempRangeDay       ,&
          TempRangeAft

     !   Solar and Temperature Params
     !   Hourly Met outgoing
     REAL(sp),DIMENSION(:),ALLOCATABLE :: &
          PhiSd           ,&  ! downward solar irradiance [W/m2]
          PhiLd           ,&  ! down longwave irradiance  [W/m2]
          Precip          ,&  ! precip [mm/h]
          Snow            ,&  ! precip [mm/h]
          Wind            ,&  ! wind   [m/s]
          Temp            ,&  ! temp   [degC]
          VapPmb          ,&  ! vapour pressure [mb]
          Pmb             ,&  ! pressure [mb]
          coszen              ! cos(theta)

  END TYPE WEATHER_GENERATOR_TYPE


CONTAINS

!===============================================================================

SUBROUTINE WGEN_INIT( WG, np, latitude, dels )

  IMPLICIT NONE

  TYPE(WEATHER_GENERATOR_TYPE) :: WG
  INTEGER, INTENT(IN) :: np
  REAL,    INTENT(IN) :: latitude(np), dels

  WG%np     = np
  WG%delT   = dels
  WG%ndtime = NINT ( REAL(SecDay)/dels )

  ALLOCATE ( WG%LatDeg             (np) )
  ALLOCATE ( WG%WindDay            (np) )
  ALLOCATE ( WG%TempMinDay         (np) ) ! Daily minimum air temp [degC]
  ALLOCATE ( WG%TempMaxDay         (np) ) ! Daily maximum air temp [degC]
  ALLOCATE ( WG%TempMinDayNext     (np) ) ! Daily minimum air temp tomorrow [degC]
  ALLOCATE ( WG%TempMaxDayPrev     (np) ) ! Daily maximum air temp yesterday [degC]
  ALLOCATE ( WG%SolarMJDay         (np) )
  ALLOCATE ( WG%WindDark           (np) )  ! Wind [m/s] during night hrs
  ALLOCATE ( WG%WindLite           (np) )  ! Wind [m/s] during daylight hrs
  ALLOCATE ( WG%SolarNorm          (np) )  ! (daily solar)/(solar const): geometry
  ALLOCATE ( WG%LatRad             (np) )  ! Latitude in radians
  ALLOCATE ( WG%DayLength          (np) )
  ALLOCATE ( WG%TimeSunsetPrev     (np) )
  ALLOCATE ( WG%TimeSunrise        (np) )
  ALLOCATE ( WG%TimeMaxTemp        (np) )
  ALLOCATE ( WG%TimeSunset         (np) )
  ALLOCATE ( WG%TempSunsetPrev     (np) )
  ALLOCATE ( WG%TempSunset         (np) )
  ALLOCATE ( WG%TempNightRate      (np) )
  ALLOCATE ( WG%TempNightRatePrev  (np) )
  ALLOCATE ( WG%TempRangeDay       (np) )
  ALLOCATE ( WG%TempRangeAft       (np) )
  ALLOCATE ( WG%PrecipDay          (np) )  ! 24hr-total precipitation [m/day] (current day = 24h from 0900 on previous day)
  ALLOCATE ( WG%SnowDay            (np) )  ! [m/d]
  ALLOCATE ( WG%PmbDay             (np) )  ! 24hr-av pressure [mb]
  ALLOCATE ( WG%VapPmbDay          (np) )  ! 24hr-av vapour pressure [mb]
  ALLOCATE ( WG%PhiSd              (np) )  ! downward solar irradiance [W/m2]
  ALLOCATE ( WG%PhiLd              (np) )  ! down longwave irradiance  [W/m2]
  ALLOCATE ( WG%Precip             (np) )  ! precip [mm/h]
  ALLOCATE ( WG%Snow               (np) )  ! precip [mm/h]
  ALLOCATE ( WG%Wind               (np) )  ! wind   [m/s]
  ALLOCATE ( WG%Temp               (np) )  ! temp   [degC]
  ALLOCATE ( WG%VapPmb             (np) )  ! vapour pressure [mb]
  ALLOCATE ( WG%Pmb                (np) )  ! pressure [mb]
  ALLOCATE ( WG%coszen             (np) )  ! cos(theta)
  
  ALLOCATE ( WG%VapPmb0900         (np) )  ! 0900 water vapour pressure [mb]
  ALLOCATE ( WG%VapPmb1500         (np) )  ! 1500 water vapour pressure [mb]
  ALLOCATE ( WG%VapPmb1500Prev     (np) )  ! 1500 (prev day) water vapour pressure [mb]
  ALLOCATE ( WG%VapPmb0900Next     (np) )  ! 0900(next day) water vapour pressure [mb]

  !MCINI
  call wgen_zero(WG)

  WG%LatDeg(:) = latitude(:)

END SUBROUTINE WGEN_INIT


subroutine wgen_zero(wg)

  implicit none

  type(weather_generator_type), intent(inout) :: wg

  wg%LatDeg            = 0.0_sp
  wg%WindDay           = 0.0_sp
  wg%TempMinDay        = 0.0_sp
  wg%TempMaxDay        = 0.0_sp
  wg%TempMinDayNext    = 0.0_sp
  wg%TempMaxDayPrev    = 0.0_sp
  wg%SolarMJDay        = 0.0_sp
  wg%PrecipDay         = 0.0_sp
  wg%SnowDay           = 0.0_sp
  wg%PmbDay            = 0.0_sp
  wg%VapPmbDay         = 0.0_sp
  wg%VapPmb0900        = 0.0_sp
  wg%VapPmb1500        = 0.0_sp
  wg%VapPmb1500Prev    = 0.0_sp
  wg%VapPmb0900Next    = 0.0_sp
  wg%WindDark          = 0.0_sp
  wg%WindLite          = 0.0_sp
  wg%SolarNorm         = 0.0_sp
  wg%LatRad            = 0.0_sp
  wg%DayLength         = 0.0_sp
  wg%TimeSunsetPrev    = 0.0_sp
  wg%TimeSunrise       = 0.0_sp
  wg%TimeMaxTemp       = 0.0_sp
  wg%TimeSunset        = 0.0_sp
  wg%TempSunsetPrev    = 0.0_sp
  wg%TempSunset        = 0.0_sp
  wg%TempNightRate     = 0.0_sp
  wg%TempNightRatePrev = 0.0_sp
  wg%TempRangeDay      = 0.0_sp
  Wg%TempRangeAft      = 0.0_sp
  wg%PhiSd             = 0.0_sp
  wg%PhiLd             = 0.0_sp
  wg%Precip            = 0.0_sp
  wg%Snow              = 0.0_sp
  wg%Wind              = 0.0_sp
  wg%Temp              = 0.0_sp
  wg%VapPmb            = 0.0_sp
  wg%Pmb               = 0.0_sp
  wg%coszen            = 0.0_sp
  
end subroutine wgen_zero


SUBROUTINE WGEN_DAILY_CONSTANTS( WG, np, YearDay )

!-------------------------------------------------------------------------------
! Routine for calculating solar and met parameters that are constant over the
! the course of a 24hr day. These have been split off from version 01 of
! SubDiurnalMet which no longer calculates a day's worth of subdiurnal met in
! one call. Subdiurnal met for ntime times is now calculated using ntime calls,
! using the day-constants calculated here.
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  TYPE(WEATHER_GENERATOR_TYPE) :: WG
  INTEGER, INTENT(IN) :: np, YearDay

  ! Local variables
  REAL(sp) :: YearRad
  REAL(sp),DIMENSION(np) :: LatDeg1
  REAL(sp),DIMENSION(np) :: TanTan
  REAL(sp),DIMENSION(np) :: HDLRad
  REAL(sp),DIMENSION(np) :: TimeSunriseNext
  REAL(sp) :: RatioWindLiteDark

! -------------------------
! Downward solar irradiance
! -------------------------
LatDeg1    = SIGN(MIN(ABS(WG%LatDeg),89.9_sp), WG%LatDeg)   ! avoid singularity at pole
WG%LatRad  = LatDeg1*Pi/180.0_sp                      ! latitude in radians
YearRad    = 2.0_sp*Pi*(YearDay-1)/365.0_sp              ! day of year in radians
!YearRadPrev = 2.0*Pi*(YearDay-2)/365.0          ! previous day of year in radians
                                                ! (Ok for YD - 2 = -1)

! DecRad = Declination in radians (+23.5 deg on 22 June, -23.5 deg on 22 Dec):
WG%DecRad = 0.006918_sp - 0.399912_sp*COS(YearRad) + 0.070257_sp*SIN(YearRad)     &
         - 0.006758_sp*COS(2.0_sp*YearRad) + 0.000907_sp*SIN(2.0_sp*YearRad)      &
         - 0.002697_sp*COS(3.0_sp*YearRad) + 0.001480_sp*SIN(3.0_sp*YearRad)
                                        ! Paltridge and Platt eq [3.7]

! Daylength: HDLRad = Half Day Length in radians (dawn:noon = noon:dusk):
TanTan = -TAN(WG%LatRad)*TAN(WG%DecRad)
WHERE (TanTan .LE. -1.0_sp)
  HDLRad = Pi                           ! polar summer: sun never sets
ELSEWHERE (TanTan .GE. 1.0_sp)
  HDLRad = 0.0_sp                          ! polar winter: sun never rises
ELSEWHERE
  HDLRad = ACOS(TanTan)                 ! Paltridge and Platt eq [3.21]
END WHERE                                  ! (HDLRad = their capital PI)
WG%DayLength = 24.0_sp*2.0_sp*HDLRad / (2.0_sp*Pi)  ! Daylength (dawn:dusk) in hours

! Daily solar irradiance without atmosphere, normalised by solar constant,
! with both energy fluxes in MJ/m2/day, calculated from solar geometry:
WG%SolarNorm =  &                          ! Paltridge and Platt eq [3.22]
  (HDLRad*SIN(WG%LatRad)*SIN(WG%DecRad) + &
  COS(WG%LatRad)*COS(WG%DecRad)*SIN(HDLRad)) / Pi

! ----
! Wind
! ----

RatioWindLiteDark = 3.0_sp                 ! (daytime wind) / (nighttime wind)
WG%WindDark  = WG%WindDay / &
     ( (24.0_sp-WG%DayLength)/24.0_sp + RatioWindLiteDark*WG%DayLength/24.0_sp )
WG%WindLite  = WG%WindDay / &
     ( (1.0_sp/RatioWindLiteDark)*(24.0_sp-WG%DayLength)/24.0_sp + WG%DayLength/24.0_sp )

! -----------
! Temperature
! -----------
! These are parameters required for the calculation of temperature according to
! Cesaraccio et al 2001 for sunrise-to-sunrise-next-day. Because we are calculating temp for
! midnight-to-midnight, we need to calculate midnight-to-sunrise temps using data for the
! previous day (-24h), hence the extra parameters denoted by *'s below which are not
! mentioned per se in Cesaraccio. Cesaraccio symbology for incoming met data is included
! here as comments for completeness:
!                                                                 Sym in Cesaraccio et al 2001
! TempMinDay                                                        Tn
! TempMaxDay                                                        Tx
! TempMinDayNext                                                    Tp
! TempMaxDayPrev  
!!$write(71, "( 1000e16.6)") WG%DecRad
!!$   write(72, "( 1000e16.6)") TAN(WG%DecRad)
!!$   write(73, "( 1000e16.6)") TAN(WG%LatRad)
!!$ write(74, "( 1000e16.6)")  WG%DayLength 
WHERE ( WG%DayLength .LT. 0.01_sp )
   ! Polar night
   WG%TimeSunrise     = 9._sp              ! Hn
   WG%TimeSunset      = 16._sp                            ! Ho
   WG%TimeMaxTemp     = 13._sp                                    ! Hx
ELSEWHERE ( WG%DayLength .GT. 23.0_sp )
   ! Polar day
   WG%TimeSunrise     = 1._sp              ! Hn
   WG%TimeSunset      = 23._sp                            ! Ho
   WG%TimeMaxTemp     = 13._sp                                    ! Hx
ELSEWHERE
   
   WG%TimeSunrise     = (ACOS(min(TAN(WG%LatRad)*TAN(WG%DecRad),0.9999_sp)))*12._sp/Pi              ! Hn
   WG%TimeSunset      = WG%TimeSunrise + WG%DayLength                             ! Ho
   WG%TimeMaxTemp     = WG%TimeSunset - MIN(4._sp,( WG%DayLength * 0.4_sp)) ! Hx
END WHERE
WG%TimeSunsetPrev  = WG%TimeSunset - 24._sp                                 ! * Ho-24h (a negative hour)
TimeSunriseNext    = WG%TimeSunrise + 24._sp                               ! Hp
WG%TempSunset      = WG%TempMaxDay - &
     (0.39_sp * (WG%TempMaxDay - WG%TempMinDayNext))               ! To
WG%TempSunsetPrev  = WG%TempMaxDayPrev - &
     (0.39_sp * (WG%TempMaxDayPrev - WG%TempMinDay))           ! * To-24h
WG%TempRangeDay    = WG%TempMaxDay - WG%TempMinDay                            ! alpha = Tx-Tn
WG%TempRangeAft    = WG%TempMaxDay - WG%TempSunset                            ! R = Tx-To
WG%TempNightRate   = (WG%TempMinDayNext - WG%TempSunset)/ &
     SQRT(TimeSunriseNext-WG%TimeSunset)                  ! b = (Tp-To)/sqrt(Hp-Ho)
WG%TempNightRatePrev = (WG%TempMinDay - WG%TempSunsetPrev)/ &
     SQRT(WG%TimeSunrise-WG%TimeSunsetPrev)                  ! * b-24h = (Tn-(To-24h))/sqrt(Hn-(Ho-24h))

END SUBROUTINE WGEN_DAILY_CONSTANTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE WGEN_SUBDIURNAL_MET(WG, np, itime)

  !-------------------------------------------------------------------------------
  ! Routine for downscaling daily met data to ntime subdiurnal values per day,
  ! at the instants ((it/ntime)*2400, it=1,ntime).
  !
  ! ALGORITHMS:
  ! * Downward solar irradiance:
  !   * Fit a daily irradiance time series based on solar geometry, normalised to
  !     observed daily total SolarMJDay. Reference: Paltridge and Platt (1976).
  ! * Downward longwave irradiance:
  !   * Computed from downscaled air temperature using Swinbank (1963) formula.
  ! * Precipitation:
  !   * Assume steady through 24 hours at average value (needs improvement?)
  ! * Wind:
  !   * Set RatioWindLiteDark = (daytime wind) / (nighttime wind), a predetermined
  !     value, typically 3. Then calculate wind by day and wind by night (both
  !     steady, but different) to make the average work out right.
  ! * Temperature:
  !   * Fit a sine wave from sunrise (TempMinDay) to peak at TempMaxDay, another
  !     from TempMaxDay to sunset, TempMinDay, and a square-root function from
  !     sunset to sunrise the next day (Cesaraccio et al 2001).
  ! * Water Vapour Pressure, Air Pressure:
  !   * Hold constant. Return rel humidity at TempMinDay as a diagnostic on obs.
  !
  ! HISTORY:
  ! * 04-sep-2003 (MRR): 01 Written and tested in Program SubDiurnalWeather
  ! * 30-nov-2007 (PRB): 02 Vectorise with deferred shape arrays, adding np dimension.
  ! * 11-dec-2007 (PRB): 03 Switch to Cesaraccio et al 2001 algorithm for temperature
  ! * 28/02/2012 (VH) :  04 cahnge precip from uniform distribution to evenly distributed over
  ! the periods 0600:0700; 0700:0800; 1500:1600; 1800:1900
  !-------------------------------------------------------------------------------
#ifdef __MPI__
  use mpi, only: MPI_Abort
#endif

  IMPLICIT NONE

  TYPE(WEATHER_GENERATOR_TYPE) :: WG
  INTEGER, INTENT(IN) :: np, itime

  ! Local variables
  REAL(sp) :: TimeNoon, adjust_fac(np)
  REAL(sp) :: TimeRad
  REAL(sp) :: rntime    ! Real version of ntime
  REAL(sp) :: ritime    ! Real version of current time
  REAL(sp),DIMENSION(np):: PhiLd_Swinbank     ! down longwave irradiance  [W/m2]
#ifdef __MPI__
  integer :: ierr
#endif

  !-------------------------------------------------------------------------------

  ritime = REAL(itime)     * WG%delT/3600._sp  ! Convert the current time to real
  rntime = REAL(WG%ndtime) * WG%delT/3600._sp  ! Convert ntime to real

  ! Instantaneous downward hemispheric solar irradiance PhiSd
  TimeNoon = ritime/rntime - 0.5_sp   ! Time in day frac (-0.5 to 0.5, zero at noon)
  TimeRad  = 2.0_sp*Pi*TimeNoon       ! Time in day frac (-Pi to Pi, zero at noon)
  WHERE (ritime >= WG%TimeSunrise .AND. ritime <= WG%TimeSunset&
       .AND.WG%SolarNorm.gt.1.e-3_sp  ) ! Sun is up
     WG%PhiSd = MAX((WG%SolarMJDay/WG%SolarNorm)  &   ! PhiSd [MJ/m2/day]
          * ( SIN(WG%DecRad)*SIN(WG%LatRad) + COS(WG%DecRad)*COS(WG%LatRad) &
          * COS(TimeRad) ), 0.0_sp)  ! fix to avoid negative PhiSD vh 13/05/08
     ! Paltridge and Platt eq [3.4]
     WG%coszen = ( SIN(WG%DecRad)*SIN(WG%LatRad) + COS(WG%DecRad)*COS(WG%LatRad)*COS(TimeRad) )
  ELSEWHERE ! sun is down
     WG%PhiSd  = 0.0_sp
     WG%coszen = 0.0_sp
  END WHERE

  WG%PhiSd    = WG%PhiSd*1.e6_sp/SecDay       ! Convert PhiSd: [MJ/m2/day] to [W/m2]
  ! -------------
  ! Precipitation
  ! -------------

  !Precip = PrecipDay*1000./rntime  ! convert from m/d to mm/h
  !Precip = (PrecipDay*1000.*9./24. + PrecipDayNext*1000.*15./24.)/rntime
  !Precip = (PrecipDay*1000.*9./24. + PrecipDayNext*1000.*15./24.)/24.  ! hourly precip [mm/h]

  IF ( ABS(ritime-REAL(INT(ritime))) .GT. 1.e-7_sp) THEN
     WRITE(*,*) "Only works for integer hourly timestep! tstep = ", ritime
     write(*,*) "cable_weathergenerator.F90!"
#ifdef __MPI__
     call MPI_Abort(0, 83, ierr) ! Do not know comm nor rank here
#else
     stop 83
#endif
  ENDIF

  IF ((WG%ndtime .eq. 24) .OR. (WG%ndtime .eq. 8)) then
     IF ((ritime >= 15._sp .AND. ritime < 16._sp).OR.(ritime >= 18._sp .AND. ritime < 19._sp)) THEN
        WG%Precip = WG%PrecipDay *1000._sp/2._sp
        WG%Snow   = WG%SnowDay*1000._sp/2._sp
     ELSE
        WG%Precip = 0._sp
        WG%Snow   = 0._sp
     ENDIF
  elseif (WG%ndtime .eq. 12) then
     IF ((ritime >= 16._sp .AND. ritime < 17._sp).OR.(ritime >= 18._sp .AND. ritime < 19._sp)) THEN
        WG%Precip = WG%PrecipDay *1000._sp/2._sp
        WG%Snow   = WG%SnowDay*1000._sp/2._sp
     ELSE
        WG%Precip = 0._sp
        WG%Snow   = 0._sp
     ENDIF
     
  endif

  ! ----
  ! Wind
  ! ----

  WHERE (ritime >= WG%TimeSunrise .AND. ritime <= WG%TimeSunset) ! Sun is up
     WG%Wind = WG%WindLite
  ELSEWHERE ! Sun is down
     WG%Wind = WG%WindDark
  END WHERE

  ! -----------
  ! Temperature
  ! -----------
  ! Calculate temperature according to Cesaraccio et al 2001, including midnight to
  ! sunrise period using previous days info, and ignoring the period from after the
  ! following midnight to sunrise the next day, normally calculated by Cesaraccio.
  WHERE ( ritime <= WG%TimeSunrise )
     ! Midnight to sunrise
     WG%Temp = WG%TempSunsetPrev + WG%TempNightRatePrev * SQRT(ritime - WG%TimeSunsetPrev)
  ELSEWHERE (ritime > WG%TimeSunrise .AND. ritime <= WG%TimeMaxTemp)
     ! Sunrise to time of maximum temperature
     WG%Temp = WG%TempMinDay + &
          WG% TempRangeDay *SIN (((ritime-WG%TimeSunrise)/ &
          (WG%TimeMaxTemp-WG%TimeSunrise))*PiBy2)
  ELSEWHERE (ritime > WG%TimeMaxTemp .AND. ritime <= WG%TimeSunset)
     ! Time of maximum temperature to sunset
     WG%Temp = WG%TempSunset + &
          WG%TempRangeAft * SIN(PiBy2 + ((ritime-WG%TimeMaxTemp)/4._sp*PiBy2))
  ELSEWHERE (ritime > WG%TimeSunset )
     ! Sunset to midnight
     WG%Temp = WG%TempSunset + WG%TempNightRate * SQRT(ritime - WG%TimeSunset)
  END WHERE

  ! -----------------------------------
  ! Water Vapour Pressure, Air Pressure
  ! -----------------------------------
  !WG%VapPmb = WG%VapPmbDay
  WG%Pmb    = WG%PmbDay

  IF (TRIM(cable_user%MetType).EQ.'bios') THEN
     
     IF (ritime <= 9._sp) THEN
        ! before 9am
        WG%VapPmb = WG%VapPmb1500Prev + (WG%VapPmb0900 - WG%VapPmb1500Prev) * (9._sp + ritime)/18._sp
     ELSEIF (ritime > 9_sp .AND. ritime <= 15._sp) THEN
        ! between 9am and 15:00
        WG%VapPmb = WG%VapPmb0900 + (WG%VapPmb1500 - WG%VapPmb0900) * (ritime - 9._sp)/(15._sp-9._sp)
     ELSEIF (ritime > 15._sp) THEN
        ! after 15:00
        WG%VapPmb  = WG%VapPmb1500 + (WG%VapPmb0900Next - WG%VapPmb1500) * (ritime - 15._sp)/18._sp
     END IF
   
  ELSE
    
     WG%VapPmb = WG%VapPmbDay
  ENDIF

  

  ! ----------------------------
  ! Downward longwave irradiance
  ! ----------------------------

  PhiLd_Swinbank = 335.97_sp * (((WG%Temp + 273.16_sp) / 293.0_sp)**6)   ! [W/m2] (Swinbank 1963)

  ! -------------------------------
  ! Alternate longwave formulation
  ! ----------------------------

  WG%PhiLd = epsilon * SBoltz * (WG%Temp + 273.16_sp)**4       ! [W/m2] (Brutsaert)

  WHERE (WG%PhiSd.GT.50.0_sp)
     adjust_fac = ((1.17_sp)**(WG%SolarNorm))/1.17_sp
  ELSEWHERE
     adjust_fac = 0.9_sp
  ENDWHERE

  WG%PhiLd = WG%PhiLd /adjust_fac * (1.0_sp + WG%PhiSd/8000._sp)       ! adjustment (formulation from Gab Abramowitz)

  WHERE ((WG%PhiLd.GT.500.0_sp).OR.(WG%PhiLd.LT.100.0_sp))
     WG%PhiLd = PhiLd_Swinbank
  ENDWHERE

  IF (ANY((WG%PhiLd.GT.750.0_sp).OR.(WG%PhiLd.LT.100.0_sp))) THEN
     !write(*,*) 'PhiLD out of range'
  ENDIF

END SUBROUTINE WGEN_SUBDIURNAL_MET

END MODULE CABLE_WEATHERGENERATOR
