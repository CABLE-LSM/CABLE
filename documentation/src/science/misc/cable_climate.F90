!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: tracking climate variables for use in phenology and potential pft modules
!
! Called from: cable_driver
!
! History: Vanessa Haverd Jan 2015

! ==============================================================================
!#define UM_BUILD YES
MODULE cable_climate_mod

  USE cable_def_types_mod, ONLY: met_type, climate_type, canopy_type, mp, &
       r_2, alloc_cbm_var, air_type, radiation_type
  USE TypeDef,              ONLY: i4b, dp
  !CABLE_LSM: see CABLE Ticket#149. yet still inclueded file?? legacy-hack??
# ifndef UM_BUILD
  USE cable_IO_vars_module, ONLY: patch
# endif
  USE CABLE_COMMON_MODULE, ONLY: CurYear, filename, cable_user
  USE casa_ncdf_module, ONLY: HANDLE_ERR

CONTAINS
  ! ==============================================================================


  SUBROUTINE cable_climate(ktau,kstart,kend,ktauday,idoy,LOY,met,climate, canopy, &
       air, rad, dels, np)


    IMPLICIT NONE

    INTEGER,      INTENT(IN) :: ktau ! integration step number
    INTEGER,      INTENT(IN) :: kstart ! starting value of ktau
    INTEGER,      INTENT(IN) :: kend ! total # timesteps in run

    INTEGER,      INTENT(IN)                  :: idoy ,LOY ! day of year (1-365) , Length oy
    INTEGER,      INTENT(IN)                  :: ktauday
    TYPE (met_type), INTENT(IN)       :: met  ! met input variables
    TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
    TYPE (canopy_type), INTENT(IN) :: canopy ! vegetation variables
    TYPE (air_type), INTENT(IN)       :: air
    TYPE (radiation_type), INTENT(IN)  :: rad        ! radiation variables
    REAL, INTENT(IN)               :: dels ! integration time setp (s)
    INTEGER,      INTENT(IN)                  :: np
    INTEGER :: d, y, k
    INTEGER, PARAMETER:: COLDEST_DAY_NHEMISPHERE = 355
    INTEGER, PARAMETER:: COLDEST_DAY_SHEMISPHERE = 172
    REAL, PARAMETER:: CoeffPT = 1.26
    REAL,      DIMENSION(mp)  :: mtemp_last, mmoist_last, phiEq, ppc, EpsA, RhoA
    INTEGER :: startyear
    INTEGER :: MonthDays(12)
    INTEGER::  DaysInMonth, nmonth, tmp
    LOGICAL :: IsLastDay ! last day of month?
    REAL, PARAMETER:: Gaero = 0.015  ! (m s-1) aerodynmaic conductance (for use in PT evap)
    REAL, PARAMETER:: Capp   = 29.09    ! isobaric spec heat air    [J/molA/K]
    REAL, PARAMETER:: SBoltz  = 5.67e-8  ! Stefan-Boltzmann constant [W/m2/K4]
    climate%doy = idoy

!!$! * Find irradiances, available energy, equilibrium latent heat flux
!!$PPc    = Gaero / ( Gaero + 4.0*SBoltz*((TempA+273.16)**3)/(RhoA*Capp) )
!!$                                                    ! PPc = Ga/(Ga+Gr)      [-]
!!$EpsA   = Epsif(TempA, Pmb)                          ! Epsi at TempA         [-]
!!$PhiSd  = SolarMJ * 1.0e6 / (DayltFrac*SecDay)       ! daylt down solar      [W/m2]
!!$PhiLd  = 335.97 * (((TempA + 273.16) / 293.0)**6)   ! daylt down thermal    [W/m2]
!!$                                                    !   (Swinbank formula)
!!$PhiAi  = (1.0-Albedo)*PhiSd + Emis *    &           ! daylt iso-avail engy  [W/m2]
!!$         (PhiLd - SBoltz*((TempA + 273.16)**4))     !   (veg + soil)
!!$PhiEq  = PhiAi * (PPc*EpsA) / (PPc*EpsA + 1.0)      ! equil ltnt heat flux  [W/m2]
!!$PhiEq  = max(PhiEq, 1.0)                            ! PhiEq > +1 W/m2, so non-negative
!!$                                    ! precipitation         [m/day]
!!$FWPT   = CoeffPT * PhiEq * ((DayltFrac*SecDay) / (RhoW*Rlat))

    ! accumulate annual evaporation and potential evaporation
    !ppc = 1.0
    RhoA = met%pmb * 100.0 / (8.314 * (met%tk)) ! air density [molA/m3]
    PPc    = Gaero / ( Gaero + 4.0*SBoltz*((met%tk**3)/(RhoA*Capp) ))

    EpSA = Epsif(met%tk - 273.16, met%pmb)
    phiEq = canopy%rniso * (PPc*EpsA) / (PPc*EpsA + 1.0)      ! equil ltnt heat flux  [W/m2]

    IF (idoy==1 .AND. MOD(ktau,ktauday)==1 ) THEN
       !  climate%evap_PT =  max(phiEq,1.0)*CoeffPT/air%rlam*dels  ! mm
       climate%evap_PT =  phiEq*CoeffPT/air%rlam*dels  ! mm
       !  climate%evap_PT = canopy%epot  ! mm
       !  climate%aevap  =   canopy%fe/air%rlam*dels ! mm
       climate%aevap = met%precip ! mm
    ELSE

       climate%evap_PT = climate%evap_PT + MAX(phiEq,1.0)*CoeffPT/air%rlam*dels  ! mm
       ! climate%evap_PT =climate%evap_PT + canopy%epot  ! mm
       ! climate%aevap =  climate%aevap + canopy%fe/air%rlam*dels ! mm
       climate%aevap = climate%aevap + met%precip ! mm
    ENDIF

    ! accumulate daily temperature, evap and potential evap
    IF(MOD(ktau,ktauday)==1) THEN
       climate%dtemp = met%tk - 273.15
       climate%dmoist = canopy%fwsoil
    ELSE
       climate%dtemp = climate%dtemp + met%tk - 273.15
       climate%dmoist = climate%dmoist + canopy%fwsoil
    ENDIF

    IF(MOD((ktau-kstart+1),ktauday)==0) THEN  ! end of day
       climate%dtemp = climate%dtemp/FLOAT(ktauday)
       climate%dmoist = climate%dmoist/FLOAT(ktauday)
    ENDIF



    IF(MOD((ktau-kstart+1),ktauday)==0) THEN  ! end of day
       ! get month and check if end of month
       IsLastDay = .FALSE.
       MonthDays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
       IF (LOY==366) MonthDays(2) = MonthDays(2) + 1
       nmonth = 1
       tmp = MonthDays(1)
       DO WHILE(tmp.LT.idoy)
          tmp = tmp +  MonthDays(nmonth+1)
          nmonth = nmonth + 1
       ENDDO

       IF (idoy == SUM(MonthDays(1:nmonth))) IsLastDay = .TRUE.

       ! On first day of year ...
       IF (idoy==1) THEN

          ! ... reset annual GDD5 counter
          climate%agdd5=0.0
          climate%agdd0=0.0
          climate%evap_PT = 0     ! annual PT evap [mm]
          climate%aevap  = 0      ! annual evap [mm]


       ENDIF

       !jhan: see CABLE Ticket#149 and above CABLE_LSM: comment
# ifndef UM_BUILD
       WHERE ((patch%latitude>=0.0 .AND. idoy==COLDEST_DAY_NHEMISPHERE).OR. &
            (patch%latitude<0.0 .AND. idoy==COLDEST_DAY_SHEMISPHERE) )

          ! In midwinter, reset GDD counter for summergreen phenology
          climate%gdd5=0.0
          climate%gdd0=0.0
       END WHERE
# endif
       ! Update GDD counters and chill day count
       climate%gdd0 = climate%gdd0 + MAX(0.0,climate%dtemp-0.0)
       climate%agdd0= climate%agdd0 + MAX(0.0,climate%dtemp-0.0)

       climate%gdd5 = climate%gdd5 + MAX(0.0,climate%dtemp-5.0)
       climate%agdd5= climate%agdd5 + MAX(0.0,climate%dtemp-5.0)
       WHERE (climate%dtemp<5.0 .AND. climate%chilldays<=365)
          climate%chilldays = climate%chilldays + 1
       ENDWHERE

       ! Save yesterday's mean temperature for the last month
       mtemp_last=climate%mtemp

       ! Update daily temperatures, and mean overall temperature, for last 31 days

       climate%mtemp=climate%dtemp
       climate%qtemp=climate%dtemp
       climate%mmoist = climate%dmoist
       DO d=1,30
          climate%dtemp_31(:,d)=climate%dtemp_31(:,d+1)
          climate%mtemp= climate%mtemp + climate%dtemp_31(:,d)
          climate%dmoist_31(:,d)=climate%dmoist_31(:,d+1)
          climate%mmoist = climate%mmoist + climate%dmoist_31(:,d)
       ENDDO
       DO d=1,90
          climate%dtemp_91(:,d)=climate%dtemp_91(:,d+1)
          climate%qtemp= climate%qtemp + climate%dtemp_91(:,d)
       ENDDO
       climate%dtemp_31(:,31)=climate%dtemp
       climate%dtemp_91(:,91)=climate%dtemp
       climate%dmoist_31(:,31)=climate%dmoist
       climate%qtemp = climate%qtemp/91.0  ! average temperature over the last quarter
       climate%mtemp = climate%mtemp/31.0 ! average temperature over the last month
       climate%mmoist = climate%mmoist/31.0 ! average moisture index over the last month

       ! Reset GDD and chill day counter if mean monthly temperature falls below base
       ! temperature

       WHERE(mtemp_last>=5.0 .AND. climate%mtemp<5.0)
          climate%gdd5=0.0
          climate%chilldays=0
       ENDWHERE

       ! On last day of month ...

       IF (IsLastDay) THEN

          ! Update mean temperature for the last 12 months
          ! atemp_mean_new = atemp_mean_old * (11/12) + mtemp * (1/12)

          climate%atemp_mean=climate%atemp_mean*(11./12.)+climate%mtemp*(1./12.)

          ! Record minimum and maximum monthly temperatures

          IF (nmonth==1) THEN
             climate%mtemp_min=climate%mtemp;
             climate%mtemp_max=climate%mtemp;
             climate%qtemp_max_last_year =  climate%qtemp_max;
             climate%qtemp_max=climate%qtemp;
          ELSE
             WHERE (climate%mtemp<climate%mtemp_min) &
                  climate%mtemp_min=climate%mtemp
             WHERE (climate%mtemp>climate%mtemp_max) &
                  climate%mtemp_max=climate%mtemp
             WHERE (climate%qtemp>climate%qtemp_max) &
                  climate%qtemp_max=climate%qtemp
          ENDIF  ! first month of year

          ! On 31 December update records of minimum monthly temperatures for the last
          ! 20 years and find minimum monthly temperature for the last 20 years

          IF (nmonth==12) THEN
             climate%nyears = climate%nyears +1


             startyear=20-MIN(19,climate%nyears-1)
             climate%mtemp_min20=0.0
             climate%mtemp_max20=0.0
             climate%alpha_PT20 = 0.0

             IF (startyear<20) THEN
                DO y=startyear,19
                   climate%mtemp_min_20(:,y)=climate%mtemp_min_20(:,y+1)
                   climate%mtemp_min20=climate%mtemp_min20+ &
                        climate%mtemp_min_20(:,y)
                   climate%mtemp_max_20(:,y)=climate%mtemp_max_20(:,y+1)
                   climate%mtemp_max20 =climate%mtemp_max20 + &
                        climate%mtemp_max_20(:,y)
                   climate%alpha_PT_20(:,y)=climate%alpha_PT_20(:,y+1)
                   climate%alpha_PT20 =climate%alpha_PT20 + &
                        climate%alpha_PT_20(:,y)
                ENDDO

                climate%mtemp_min20=climate%mtemp_min20/REAL(20-startyear)
                climate%mtemp_max20=climate%mtemp_max20/REAL(20-startyear)
                climate%alpha_PT20=climate%alpha_PT20/REAL(20-startyear)
             ELSE
                ! only occurs when climate%nyears = 1
                climate%mtemp_min20 = climate%mtemp_min
                climate%mtemp_max20 = climate%mtemp_max
                climate%alpha_PT20 = climate%alpha_PT
             ENDIF

             climate%mtemp_min_20(:,20)=climate%mtemp_min
             climate%mtemp_max_20(:,20)=climate%mtemp_max


             climate%alpha_PT = MAX(climate%aevap/climate%evap_PT, 0.0)     ! ratio of annual evap to annual PT evap
             climate%alpha_PT_20(:,20)=climate%alpha_PT


             CALL biome1_pft(climate,np)

          ENDIF  ! last month of year



       ENDIF     ! last day of month


    ENDIF ! end of day

  END SUBROUTINE cable_climate
  !###############################################################################

  ELEMENTAL FUNCTION Epsif(TC,Pmb)
    !-------------------------------------------------------------------------------
    ! At temperature TC [deg C] and pressure Pmb [mb], return
    ! epsi = (RLAM/CAPP) * d(sat spec humidity)/dT [(kg/kg)/K], from Teten formula.
    ! MRR, xx/1987, 27-jan-94
    ! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
    !                 just like intrinsic functions.
    ! MRR, 28-mar-05: use Rlat [J/molW], Capp [J/molA/K] from MODULE Constants,
    !                 to ensure consistency with other uses of Rlat, Capp
    ! MRR, 28-mar-05: Remove dependence of Rlat (latent heat vaporisation of water)
    !                 on temperature, use value at 20 C
    !-------------------------------------------------------------------------------
    !USE TypeDef
    !USE Constants
    IMPLICIT NONE
    REAL,INTENT(in):: TC, Pmb       ! temp [deg C], pressure [mb]
    REAL:: Epsif                    ! epsi
    REAL:: TCtmp, ES, dESdT         ! local
    REAL,PARAMETER:: A = 6.106      ! Teten coefficients
    REAL,PARAMETER:: B = 17.27      ! Teten coefficients
    REAL,PARAMETER:: C = 237.3      ! Teten coefficients
    REAL,PARAMETER:: Rlat      = 44140.0  ! lat heat evap H2O at 20C  [J/molW]
    REAL,PARAMETER:: Capp      = 29.09    ! isobaric spec heat air    [J/molA/K]
    !-------------------------------------------------------------------------------
    TCtmp = TC                          ! preserve TC
    IF (TCtmp.GT.100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
    IF (TCtmp.LT.-40.0) TCtmp = -40.0
    ES    = A*EXP(B*TCtmp/(C+TCtmp))    ! sat vapour pressure
    dESdT = ES*B*C/(C+TCtmp)**2         ! d(sat VP)/dT: (mb/K)
    Epsif = (Rlat/Capp) * dESdT / Pmb   ! dimensionless (ES/Pmb = molW/molA)

  END FUNCTION Epsif

  !###############################################################################
  ! ==============================================================================
  !
  SUBROUTINE biome1_pft(climate, np)
    IMPLICIT NONE

    TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
    INTEGER, INTENT(IN) :: np
    INTEGER :: k, j, npft
    INTEGER, ALLOCATABLE :: pft_biome1(:,:)
    REAL, ALLOCATABLE:: alpha_PT_scaled(:)

    ! TABLE 1 , Prentice et al. J. Biogeog., 19, 117-134, 1992
    ! pft_biome1: Trees (1)tropical evergreen; (2) tropical raingreen; (3) warm temp evergreen ;
    ! (4) temperate summergreen; (5) cool-temperate conifer; (6) boreal evergreen conifer;
    ! (7) boreal summergreen;
    ! Non-trees: (8) sclerophyll/succulent; (9) warm grass/shrub; (10) cool grass/shrub;
    ! (11) cold grass/shrub; (12) hot desert shrub; (13) cold desert shrub.

    ALLOCATE(pft_biome1(np,4))
    ALLOCATE(alpha_PT_scaled(np))
    alpha_PT_scaled =  MIN(climate%alpha_PT20, 1.0)

    DO k=1,np

       pft_biome1(k,:) = 999

       IF (climate%mtemp_min20(k) .GE. 15.5) THEN
          IF (alpha_PT_scaled(k).GE.0.80) THEN
             pft_biome1(k,1) = 1
             IF (alpha_PT_scaled(k).LE.0.85) THEN
                pft_biome1(k,2) = 2
             ENDIF
          ELSEIF (alpha_PT_scaled(k).GE.0.4 .AND. alpha_PT_scaled(k).LT.0.80) THEN
             pft_biome1(k,1) = 2
          ENDIF
       ENDIF


       IF (climate%mtemp_min20(k).GE.5 .AND.alpha_PT_scaled(k).GE.0.4 &
            .AND. pft_biome1(k,1).EQ.999 ) THEN
          pft_biome1(k,1) = 3
       ENDIF


       IF (climate%mtemp_min20(k).GE.-15 .AND. climate%mtemp_min20(k).LE.15.5 .AND. &
            alpha_PT_scaled(k).GE.0.35 .AND. climate%agdd5(k).GT.1200 & !
            .AND. pft_biome1(k,1).GT.3 ) THEN
          pft_biome1(k,1) = 4
       ENDIF

       IF (climate%mtemp_min20(k).GE.-19 .AND. climate%mtemp_min20(k).LE.5 .AND. &
            alpha_PT_scaled(k).GE.0.35 .AND. climate%agdd5(k).GT.900 )  THEN
          IF (pft_biome1(k,1).GT.4 ) THEN
             pft_biome1(k,1) = 5
          ELSEIF (pft_biome1(k,1).EQ.4 ) THEN
             pft_biome1(k,2) = 5
          ENDIF
       ENDIF

       IF (climate%mtemp_min20(k).GE.-35 .AND. climate%mtemp_min20(k).LE.-2 .AND. &
            alpha_PT_scaled(k).GE.0.35 .AND. climate%agdd5(k).GT.350 )  THEN
          IF (pft_biome1(k,1).EQ.999 ) THEN
             pft_biome1(k,1) = 6
          ELSEIF (pft_biome1(k,2).EQ.999 ) THEN
             pft_biome1(k,2) = 6
          ELSE
             pft_biome1(k,3) = 6
          ENDIF
       ENDIF

       IF ( climate%mtemp_min20(k).LE. 5 .AND. &
            alpha_PT_scaled(k).GE.0.45 .AND. climate%agdd5(k).GT.350 )  THEN
          IF (pft_biome1(k,1).EQ.999 ) THEN
             pft_biome1(k,1) = 7
          ELSEIF (pft_biome1(k,2).EQ.999 ) THEN
             pft_biome1(k,2) = 7
          ELSEIF (pft_biome1(k,3).EQ.999 ) THEN
             pft_biome1(k,3) = 7
          ELSE
             pft_biome1(k,4) = 7
          ENDIF
       ENDIF

       IF (climate%mtemp_min20(k).GE.5 .AND.alpha_PT_scaled(k).GE.0.2 &
            .AND. pft_biome1(k,1).EQ.999 ) THEN
          pft_biome1(k,1) = 8
       ENDIF

       IF (climate%mtemp_max20(k).GE.22 .AND.alpha_PT_scaled(k).GE.0.1 &
            .AND. pft_biome1(k,1).EQ.999 ) THEN
          pft_biome1(k,1) = 9
       ENDIF

       IF (climate%agdd5(k).GE.500 .AND.alpha_PT_scaled(k).GE.0.33 &
            .AND. pft_biome1(k,1).EQ.999 ) THEN
          pft_biome1(k,1) = 10
       ENDIF

       IF (climate%agdd0(k).GE.100 .AND.alpha_PT_scaled(k).GE.0.33) THEN
          IF (pft_biome1(k,1).EQ.999 ) THEN
             pft_biome1(k,1) = 11
          ELSEIF (pft_biome1(k,1).EQ.10) THEN
             pft_biome1(k,2) = 11
          ENDIF
       ENDIF

       IF (climate%mtemp_max20(k).GE.22 .AND. pft_biome1(k,1).EQ.999 ) THEN
          pft_biome1(k,1) = 12
       ENDIF

       IF (climate%agdd0(k).GE.100 .AND. pft_biome1(k,1).EQ.999 ) THEN
          pft_biome1(k,1) = 13
       ENDIF

       ! end of evironmental constraints on pft
       npft = 0
       DO j=1,4
          IF (pft_biome1(k,j).NE.999) npft = npft+1
       ENDDO
       !     MAP to Biome1 biome and CABLE pft
       ! (1) Tropical Rainforest
       IF(pft_biome1(k,1)==1 .AND. npft .EQ.1 ) THEN
          climate%biome(k) = 1
          climate%iveg(k) = 2
       ENDIF

       ! (2) Tropical Seasonal forest
       IF(pft_biome1(k,1)==1 .AND.pft_biome1(k,2)==2.AND. npft .EQ.2 ) THEN
          climate%biome(k) = 2
          climate%iveg(k) = 2
       ENDIF

       ! (3) Tropical dry forest/savanna
       IF(pft_biome1(k,1)==2.AND. npft .EQ.1 ) THEN
          climate%biome(k) = 3
          climate%iveg(k) = 2  ! N.B. need to include c4 grass
       ENDIF


       ! (4) Broad-leaved evergreen/warm mixed-forest
       IF(pft_biome1(k,1)==3.AND. npft .EQ.1 ) THEN
          climate%biome(k) = 4
          climate%iveg(k) = 2
       ENDIF

       ! (5) Temperate deciduous forest
       IF(pft_biome1(k,1)==4.AND.pft_biome1(k,2)==5.AND. &
            pft_biome1(k,3)==7 .AND. npft .EQ.3 ) THEN
          climate%biome(k) = 5
          climate%iveg(k) = 4
       ENDIF

       ! (6) Cool mixed forest
       IF(pft_biome1(k,1)==4.AND.pft_biome1(k,2)==5.AND. &
            pft_biome1(k,3)==6 .AND.  pft_biome1(k,4)==7 &
            .AND. npft .EQ.4 ) THEN
          climate%biome(k) = 6
          climate%iveg(k) = 4
       ENDIF

       ! (7) Cool conifer forest

       IF(pft_biome1(k,1)==5.AND.pft_biome1(k,2)==6.AND. &
            pft_biome1(k,3)==7 .AND. npft .EQ.3 ) THEN
          climate%biome(k) = 7
          climate%iveg(k) = 1
       ENDIF
       ! (8) Taiga
       IF(pft_biome1(k,1)==6.AND.pft_biome1(k,2)==7 .AND. npft .EQ.2 ) THEN
          climate%biome(k) = 8
          climate%iveg(k) = 1
       ENDIF

       ! (9) Cold mixed forest
       IF(pft_biome1(k,1)==5.AND.pft_biome1(k,2)==7 .AND. npft .EQ.2 ) THEN
          climate%biome(k) = 9
          climate%iveg(k) = 1
       ENDIF

       ! (10) Cold deciduous forest
       IF(pft_biome1(k,1)==7 .AND. npft .EQ.1 ) THEN
          climate%biome(k) = 10
          climate%iveg(k) = 3
       ENDIF

       ! (11) Xerophytic woods/scrub
       IF(pft_biome1(k,1)==8 .AND. npft .EQ.1 ) THEN
          climate%biome(k) = 11
          climate%iveg(k) = 5
       ENDIF

       ! (12) Warm grass/shrub
       IF(pft_biome1(k,1)==9 .AND. npft .EQ.1 ) THEN
          climate%biome(k) = 12
          climate%iveg(k) = 5  ! include C4 grass tile ?
       ENDIF

       ! (13) Cool grass/shrub
       IF(pft_biome1(k,1)==10 .AND.pft_biome1(k,2)==11 .AND.  npft .EQ.2 ) THEN
          climate%biome(k) = 13
          climate%iveg(k) = 5  ! include C3 grass tile ?
       ENDIF

       ! (14) Tundra
       IF(pft_biome1(k,1)==11 .AND. npft .EQ.1 ) THEN
          climate%biome(k) = 14
          climate%iveg(k) = 8  !
       ENDIF

       ! (15) Hot desert
       IF(pft_biome1(k,1)==12 .AND. npft .EQ.1 ) THEN
          climate%biome(k) = 15
          climate%iveg(k) = 14  !
       ENDIF

       ! (16) Semidesert
       IF(pft_biome1(k,1)==13 .AND. npft .EQ.1 ) THEN
          climate%biome(k) = 16
          climate%iveg(k) = 5  !
       ENDIF

       ! (17) Ice/polar desert
       IF(climate%biome(k)==999) THEN
          climate%biome(k) = 17
          climate%iveg(k) = 17  !
       ENDIF

       !jhan: see CABLE Ticket#149 and above CABLE_LSM: comment
# ifndef UM_BUILD
       ! check for DBL or NEL in SH: set to EBL instead
       IF ((climate%iveg(k)==1 .OR.climate%iveg(k)==3 .OR. climate%iveg(k)==4) &
            .AND. patch(k)%latitude<0) THEN
          climate%iveg(k) = 2
       ENDIF
# endif

       !"(/grass:1/shrub:2/woody:3"
       !1,3,Evergreen Needleleaf Forest
       !2,3,Evergreen Broadleaf Forest,,,,,,,,,,,,,,,,,,
       !3,3,Deciduous Needleleaf Forest,,,,,,,,,,,,,,,,,,
       !4,3,Deciduous Broadleaf Forest,,,,,,,,,,,,,,,,,,
       !5,2,shrub,,,,,,,,,,,,,,,,,,
       !6,1,C3 grass,,,,,,,,,,,,,,,,,,
       !7,1,C4 grass,,,,,,,,,,,,,,,,,,
       !8,1,tundra,,,,,,,,,
       ! 1.



    END DO


  END SUBROUTINE BIOME1_PFT

  ! ==============================================================================

  SUBROUTINE climate_init ( climate,np, ktauday )
    IMPLICIT NONE

    TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
    INTEGER, INTENT(IN) :: np, ktauday
    INTEGER :: d

    CALL alloc_cbm_var(climate,np)

    IF (cable_user%climate_fromzero .OR. .NOT.cable_user%call_climate) THEN

       ! Maciej
       !   DO d=1,31
       !climate%dtemp_31(:,d)= climate%dtemp
       ! climate%dmoist_31(:,d)= climate%dmoist
       !      climate%dtemp_31(:,d)= 0
       !      climate%dmoist_31(:,d)= 0
       !
       !   ENDDO
       climate%dtemp_31(:,:)= 0
       climate%dtemp_91(:,:)= 0
       climate%dmoist_31(:,:)= 0
       climate%atemp_mean=0



       climate%nyears = 0
       climate%chilldays = 0
       climate%mtemp = 0
       climate%mmoist = 0
       climate%mtemp_min = 0
       climate%mtemp_max=0
       climate%qtemp_max=0
       climate%qtemp_max_last_year = 0
       climate%mtemp_min20 =0
       climate%mtemp_max20=0
       climate%AGDD5=0
       climate%GDD5=0
       climate%AGDD0=0
       climate%GDD0=0
       climate%alpha_PT=0
       climate%evap_PT=0
       climate%aevap=0
       climate%mtemp_min_20=0
       climate%mtemp_max_20=0
       climate%alpha_PT_20=0
       climate%iveg = 999
       climate%biome = 999


    ELSE
       CALL READ_CLIMATE_RESTART_NC (climate, ktauday)

    ENDIF
    !else
    ! CALL READ_CLIMATE_RESTART_NC (climate)

    !endif

  END SUBROUTINE climate_init

  ! ==============================================================================

  SUBROUTINE WRITE_CLIMATE_RESTART_NC ( climate, ktauday )

    USE netcdf


    IMPLICIT NONE

    TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables
    INTEGER, INTENT(IN) :: ktauday
    INTEGER*4 :: mp4
    INTEGER*4, PARAMETER   :: pmp4 =0
    INTEGER, PARAMETER   :: fmp4 = KIND(pmp4)
    INTEGER*4   :: STATUS
    INTEGER*4   :: FILE_ID, land_ID, nyear_ID, nday_ID, ndayq_ID, i
    CHARACTER :: CYEAR*4, FNAME*99,dum*50

    ! 0 dim arrays
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(20) :: A1
    ! 1 dim arrays (integer) (npt )
    CHARACTER(len=20),DIMENSION(3) :: AI1
    ! 2 dim arrays (npt,20)
    CHARACTER(len=20),DIMENSION(3) :: A2
    ! 2 dim arrays (npt,31)
    CHARACTER(len=20),DIMENSION(2) :: A3
    ! 2 dim arrays (npt,91)
    CHARACTER(len=20),DIMENSION(1) :: A4


    INTEGER*4 ::  VID0(SIZE(A0)),VID1(SIZE(A1)),VIDI1(SIZE(AI1)), &
         VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))

    mp4=INT(mp,fmp4)
    A0(1) = 'nyears'
    A0(2) = 'year'

    A1(1) = 'latitude'
    A1(2) = 'longitude'
    A1(3) = 'dtemp'
    A1(4) = 'mtemp'
    A1(5) = 'qtemp'
    A1(6) = 'mtemp_min'
    A1(7) = 'mtemp_max'
    A1(8) = 'qtemp_max'
    A1(9) = 'qtemp_max_last_year'
    A1(10) = 'mtemp_min20'
    A1(11) = 'mtemp_max20'
    A1(12) = 'atemp_mean'
    A1(13) = 'AGDD5'
    A1(14) = 'GDD5'
    A1(15) = 'AGDD0'
    A1(16) = 'GDD0'
    A1(17) = 'mmoist'
    A1(18) = 'evap_PT'
    A1(19) = 'aevap'
    A1(20)  = 'alpha_PT'

    AI1(1) = 'chilldays'
    AI1(2) = 'iveg'
    AI1(3) = 'biome'

    A2(1) = 'mtemp_min_20'
    A2(2) = 'mtemp_max_20'
    A2(3) = 'alpha_PT_20'

    A3(1) = 'dtemp_31'
    A3(2) = 'dmoist_31'

    A4(1) = 'dtemp_91'

# ifndef UM_BUILD

    ! Get File-Name
    WRITE(CYEAR, FMT='(I4)') CurYear + 1
    fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
         '_climate_rst.nc'
    ! Create NetCDF file:
    STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", "01/01/"//CYEAR  )

    ! Define dimensions:
    ! Land (number of points)
    STATUS = NF90_def_dim(FILE_ID, 'land'   , mp4     , land_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    ! number of years (stored for 20 y running means0
    STATUS = NF90_def_dim(FILE_ID, 'nyear' , 20 , nyear_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    ! number of days (stored for 31 day monthly means)
    STATUS = NF90_def_dim(FILE_ID, 'nday' , 31 , nday_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    ! number of days (stored for 91 day quarterly means)
    STATUS = NF90_def_dim(FILE_ID, 'ndayq' , 91 , ndayq_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    DO i = 1, SIZE(A0)
       STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_INT ,VID0(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO


    DO i = 1, SIZE(A1)
       STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(AI1)
       STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID/),VIDI1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A2)
       STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,nyear_ID/),VID2(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A3)
       STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,nday_ID/),VID3(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A4)
       STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,ndayq_ID/),VID4(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! PUT nyears and current year
    STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), climate%nyears )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), CurYear )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT LAT / LON
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(1), patch%latitude )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(2), patch%longitude )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT VARS
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(3), climate%dtemp )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(4), climate%mtemp )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(5), climate%qtemp )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(6), climate%mtemp_min )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(7), climate%mtemp_max )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(8), climate%qtemp_max )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(9), climate%qtemp_max_last_year )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), climate%mtemp_min20 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), climate%mtemp_max20  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), climate%atemp_mean  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), climate%AGDD5  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(14), climate%GDD5  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), climate%AGDD0 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), climate%GDD0 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), climate%mmoist )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), climate%evap_PT )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), climate%aevap )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(20), climate%alpha_PT )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), climate%chilldays )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(2), climate%iveg )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(3), climate%biome )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), climate%mtemp_min_20 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), climate%mtemp_max_20 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), climate%alpha_PT_20 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), climate%dtemp_31 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), climate%dmoist_31 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), climate%dtemp_91 )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

# endif

  END SUBROUTINE WRITE_CLIMATE_RESTART_NC
  ! ==============================================================================

  SUBROUTINE READ_CLIMATE_RESTART_NC ( climate, ktauday )

    USE netcdf


    IMPLICIT NONE

    TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
    INTEGER, INTENT(IN) :: ktauday
    INTEGER*4 :: mp4
    INTEGER*4, PARAMETER   :: pmp4 =0
    INTEGER, PARAMETER   :: fmp4 = KIND(pmp4)
    INTEGER*4   :: STATUS
    INTEGER*4   :: FILE_ID, land_ID, nyear_ID, nday_ID, dID, i, land_dim
    CHARACTER :: CYEAR*4, FNAME*99,dum*50

    ! 0 dim arrays
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(20) :: A1
    ! 1 dim arrays (integer) (npt )
    CHARACTER(len=20),DIMENSION(3) :: AI1
    ! 2 dim arrays (npt,20)
    CHARACTER(len=20),DIMENSION(3) :: A2
    ! 2 dim arrays (npt,31)
    CHARACTER(len=20),DIMENSION(2) :: A3
    ! 2 dim arrays (npt,91)
    CHARACTER(len=20),DIMENSION(1) :: A4

    REAL(r_2), DIMENSION(mp)          :: LAT, LON, TMP
    REAL(r_2)                         :: TMP2(mp,20),TMP3(mp,31),TMP4(mp,91)
    INTEGER*4 :: TMPI(mp), TMPI0
    LOGICAL            ::  EXISTFILE

    mp4=INT(mp,fmp4)
    A0(1) = 'nyears'
    A0(2) = 'year'

    A1(1) = 'latitude'
    A1(2) = 'longitude'
    A1(3) = 'dtemp'
    A1(4) = 'mtemp'
    A1(5) = 'qtemp'
    A1(6) = 'mtemp_min'
    A1(7) = 'mtemp_max'
    A1(8) = 'qtemp_max'
    A1(9) = 'qtemp_max_last_year'
    A1(10) = 'mtemp_min20'
    A1(11) = 'mtemp_max20'
    A1(12) = 'atemp_mean'
    A1(13) = 'AGDD5'
    A1(14) = 'GDD5'
    A1(15) = 'AGDD0'
    A1(16) = 'GDD0'
    A1(17) = 'mmoist'
    A1(18) = 'evap_PT'
    A1(19) = 'aevap'
    A1(20)  = 'alpha_PT'

    AI1(1) = 'chilldays'
    AI1(2) = 'iveg'
    AI1(3) = 'biome'

    A2(1) = 'mtemp_min_20'
    A2(2) = 'mtemp_max_20'
    A2(3) = 'alpha_PT_20'

    A3(1) = 'dtemp_31'
    A3(2) = 'dmoist_31'

    A4(1) = 'dtemp_91'


    ! Get File-Name
    WRITE(CYEAR, FMT='(I4)') CurYear + 1
    fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
         '_climate_rst.nc'

    INQUIRE( FILE=TRIM( fname ), EXIST=EXISTFILE )

    IF ( .NOT.EXISTFILE) WRITE(*,*) fname, ' does not exist!!'

# ifndef UM_BUILD
    ! Open NetCDF file:
    STATUS = NF90_OPEN(fname, NF90_NOWRITE, FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


    ! dimensions:
    ! Land (number of points)
    STATUS = NF90_INQ_DIMID(FILE_ID, 'land'   , dID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    ! number of years (stored for 20 y running means
    STATUS = NF90_INQ_DIMID(FILE_ID, 'nyear'   , dID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    ! number of days (stored for 31 d monthly means
    STATUS = NF90_INQ_DIMID(FILE_ID, 'nday'   , dID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    IF ( land_dim .NE. SIZE(patch%latitude)) THEN
       WRITE(*,*) "Dimension misfit, ", fname
       WRITE(*,*) "land_dim", land_dim, SIZE(patch%latitude)
       !     STOP
    ENDIF

    ! LAT & LON
    STATUS = NF90_INQ_VARID( FILE_ID, A1(1), dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_GET_VAR( FILE_ID, dID, LAT )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


    STATUS = NF90_INQ_VARID( FILE_ID, A1(2), dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_GET_VAR( FILE_ID, dID, LON )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! READ scalar fields
    DO i = 1, SIZE(A0)
       STATUS = NF90_INQ_VARID( FILE_ID, A0(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMPI0 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A0(i)))
       CASE ('nyears'      ) ; climate%nyears      = TMPI0
       END SELECT
    END DO


    ! READ 1-dimensional real fields
    DO i = 3, SIZE(A1)
       STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A1(i)))
       CASE ('mtemp'      ) ; climate%mtemp       = TMP
       CASE ('qtemp'      ) ; climate%mtemp       = TMP
       CASE ('mtemp_min'   ) ; climate%mtemp_min   = TMP
       CASE ('mtemp_max'  ) ; climate%mtemp_max  = TMP
       CASE ('qtemp_max'  ) ; climate%qtemp_max  = TMP
       CASE ('qtemp_max_last_year'  ) ; climate%qtemp_max_last_year  = TMP
       CASE ('mtemp_min20' ) ; climate%mtemp_min20 = TMP
       CASE ('mtemp_max20'  ) ; climate%mtemp_max20  = TMP
       CASE ('atemp_mean'  ) ; climate%atemp_mean  = TMP
       CASE ('AGDD5'  ) ; climate%AGDD5  = TMP
       CASE ('GDD5'  ) ; climate%GDD5  = TMP
       CASE ('AGDD0'  ) ; climate%AGDD0  = TMP
       CASE ('GDD0'  ) ; climate%GDD0  = TMP
       CASE ('evap_PT'  ) ; climate%evap_PT  = TMP
       CASE ('aevap'  ) ; climate%aevap  = TMP
       CASE ('alpha_PT'  ) ; climate%alpha_PT  = TMP
       END SELECT
    END DO

    ! READ 1-dimensional integer fields
    DO i = 1, SIZE(AI1)

       WRITE(*,*)  TRIM(AI1(i))
       STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMPI )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(AI1(i)))
       CASE ('chilldays'      ) ; climate%chilldays      = TMPI
       CASE ('iveg'      ) ; climate%iveg     = TMPI
       CASE ('biome'      ) ; climate%biome     = TMPI
       END SELECT
    END DO


    ! READ 2-dimensional fields (nyear)
    DO i = 1, SIZE(A2)
       WRITE(*,*)  TRIM(A2(i))
       STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A2(i)))
       CASE ('mtemp_min_20' ) ; climate%mtemp_min_20 = TMP2
       CASE ('mtemp_max_20' ) ; climate%mtemp_max_20 = TMP2
       CASE ('alpha_PT_20' ) ; climate%alpha_PT_20 = TMP2
       END SELECT
    END DO

    ! READ 2-dimensional fields (nday)
    DO i = 1, SIZE(A3)
       STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A3(i)))
       CASE ('dtemp_31' ) ; climate%dtemp_31 = TMP3
       CASE ('dmoist_31' ) ; climate%dmoist_31 = TMP3
       END SELECT
    END DO

    ! READ 2-dimensional fields (ndayq)
    DO i = 1, SIZE(A4)
       STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A4(i)))
       CASE ('dtemp_91' ) ; climate%dtemp_91 = TMP4
       END SELECT
    END DO


    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
# endif


  END SUBROUTINE  READ_CLIMATE_RESTART_NC


END MODULE cable_climate_mod
