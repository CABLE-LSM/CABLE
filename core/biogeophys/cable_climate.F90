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
! Aug 2017: additional leaf-level met variables stored for use in optimising ratio
! of Jmax to Vcmax
! ==============================================================================
MODULE cable_climate_mod

 Use cable_def_types_mod, ONLY: met_type, climate_type, canopy_type,soil_snow_type, mp, &
      r_2, alloc_cbm_var, air_type, radiation_type, veg_parameter_type
 USE TypeDef,              ONLY: i4b, dp
 USE cable_IO_vars_module, ONLY: patch
 USE CABLE_COMMON_MODULE, ONLY: CurYear, filename, cable_user, HANDLE_ERR

CONTAINS
! ==============================================================================


SUBROUTINE cable_climate(ktau,kstart,kend,ktauday,idoy,LOY,met,climate, canopy, &
     veg,ssnow,air, rad, dels, np)


  IMPLICIT NONE

  INTEGER,      INTENT(IN) :: ktau ! integration step number
  INTEGER,      INTENT(IN) :: kstart ! starting value of ktau
  INTEGER,      INTENT(IN) :: kend ! total # timesteps in run

  INTEGER,      INTENT(IN)                  :: idoy ,LOY ! day of year (1-365) , Length oy
  INTEGER,      INTENT(IN)                  :: ktauday
  TYPE (met_type), INTENT(IN)       :: met  ! met input variables
  TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
  TYPE (canopy_type), INTENT(IN) :: canopy ! vegetation variables
  TYPE (soil_snow_type), INTENT(IN) :: ssnow
  TYPE (air_type), INTENT(IN)       :: air
  TYPE (radiation_type), INTENT(IN)  :: rad        ! radiation variables
  TYPE (veg_parameter_type), INTENT(IN)  :: veg  ! vegetation parameters
  REAL, INTENT(IN)               :: dels ! integration time setp (s)
  INTEGER,      INTENT(IN)                  :: np
  INTEGER :: d, y, k
  INTEGER, PARAMETER:: COLDEST_DAY_NHEMISPHERE = 355
  INTEGER, PARAMETER:: COLDEST_DAY_SHEMISPHERE = 172
  real, PARAMETER:: CoeffPT = 1.26
  real,      dimension(mp)  :: mtemp_last, mmoist_last, phiEq, ppc, EpsA, RhoA
  integer :: startyear
  integer :: MonthDays(12)
  integer::  DaysInMonth, nmonth, tmp, nsd
  logical :: IsLastDay ! last day of month?
  real, PARAMETER:: Gaero = 0.015  ! (m s-1) aerodynmaic conductance (for use in PT evap)
  real, PARAMETER:: Capp   = 29.09    ! isobaric spec heat air    [J/molA/K]
  real, PARAMETER:: SBoltz  = 5.67e-8  ! Stefan-Boltzmann constant [W/m2/K4]
  real, PARAMETER:: moisture_min = 0.15 ! threshold for setting "growing moisture days", as required for drought-deciduous phenology
  real, PARAMETER:: T1 = 0.0, T2 = -3.0, T3 = -4.0, T6 = -5.0 ! for computing fractional spring recovery
  real, PARAMETER:: ffrost = 0.1, fdorm0 = 0.15  ! for computing fractional spring recovery
  real, PARAMETER:: gdd0_rec0 = 500.0
  real, dimension(mp) :: f1, f2, frec0, dkbdi

  nsd = ktauday*5 ! number of subdirunal time-steps
                  ! to be accumulated for storing variables needed to implement
                  ! coordination of photosynthesis

  climate%doy = idoy

  ! accumulate annual evaporation and potential evaporation
  RhoA = met%pmb * 100.0 / (8.314 * (met%tk)) ! air density [molA/m3]
  PPc    = Gaero / ( Gaero + 4.0*SBoltz*((met%tk**3)/(RhoA*Capp) ))
  
  EpSA = Epsif(met%tk - 273.16, met%pmb)
  phiEq = canopy%rniso * (PPc*EpsA) / (PPc*EpsA + 1.0)      ! equil ltnt heat flux  [W/m2]

  IF (idoy==1 .and. MOD(ktau,ktauday)==1 ) THEN   ! first time step of year
        climate%evap_PT =  phiEq*CoeffPT/2.5014e6*dels  ! mm
        climate%aevap = met%precip ! mm
        climate%fapar_ann_max = 0.0 
  ELSE
     climate%evap_PT = climate%evap_PT + phiEq*CoeffPT/2.5014e6*dels  ! mm
    ! climate%evap_PT = climate%evap_PT + max(phiEq,1.0)*CoeffPT/air%rlam*dels  ! mm
    ! climate%evap_PT =climate%evap_PT + canopy%epot  ! mm
    ! climate%aevap =  climate%aevap + canopy%fe/air%rlam*dels ! mm
     climate%aevap = climate%aevap + met%precip ! mm
  ENDIF

  ! accumulate daily temperature, evap and potential evap
  IF(MOD(ktau,ktauday)==1) THEN
     climate%dtemp = met%tk - 273.15
     climate%dmoist = sum(ssnow%wb(:,:)*veg%froot(:,:),2)
     climate%dtemp_min =  climate%dtemp
     climate%dtemp_max =  climate%dtemp
     climate%drhum = met%rhum
     climate%du10_max = met%u10
     climate%dprecip = met%precip
  ELSE
     climate%dtemp = climate%dtemp + met%tk - 273.15
     climate%dmoist =  climate%dmoist + sum(ssnow%wb(:,:)*veg%froot(:,:),2)
     climate%dtemp_min = min(met%tk - 273.15, climate%dtemp_min)
     climate%dtemp_max = max(met%tk - 273.15, climate%dtemp_max)
     climate%drhum = climate%drhum + met%rhum
     climate%du10_max = max(met%u10, climate%du10_max)
     climate%dprecip = climate%dprecip + met%precip
  ENDIF

  IF(MOD((ktau-kstart+1),ktauday)==0) THEN  ! end of day
     ! compute daily averages
     climate%dtemp = climate%dtemp/FLOAT(ktauday)
     climate%dmoist = climate%dmoist/FLOAT(ktauday)
     climate%drhum = climate%drhum/FLOAT(ktauday)

    IF(cable_user%CALL_BLAZE) then
     ! update days since last rain and precip since last day without rain
     WHERE (climate%dprecip .gt. 0.01)
        WHERE (climate%DSLR .gt. 0)
           climate%last_precip = climate%dprecip
        ELSEWHERE
           climate%last_precip = climate%last_precip + climate%dprecip
        END WHERE
        climate%DSLR = 0  ! reset days since last rain
     ELSEWHERE
        climate%DSLR =climate%DSLR + 1 
     END WHERE

     ! calculate Keetch-Byram Drought Index
     WHERE ( climate%DSLR == 0 ) 
        WHERE ( climate%last_precip > 5. ) 
           dkbdi = 5. - climate%last_precip
        ELSEWHERE
           dkbdi = 0.
        END WHERE
        
     ELSEWHERE
        dkbdi = (( 800. - climate%KBDI ) * (.968 * EXP(.0486 * &
             (climate%dtemp_max * 9./5. + 32.)) &
             - 8.3) / 1000. / (1. + 10.88 * &
             EXP(-.0441 * climate%aprecip_av20 /25.4)) * .254)
     ENDWHERE
     climate%KBDI = max(0.0, dkbdi + climate%KBDI)

     ! calcululate MacArthur Drought-Factor D
     climate%D_MacArthur = .191 * (climate%KBDI + 104.) * &
                            (real(climate%DSLR) + 1.)**1.5 / &
                            ( 3.52 * (real(climate%DSLR) + 1.)**1.5 &
                            + climate%last_precip - 1. )
     climate%D_MacArthur =  MAX(0.,MIN(10.,climate%D_MacArthur))
 

        ! MacArthur FFDI
     climate%FFDI = 2. * EXP( -.45 + .987 * LOG(climate%D_MacArthur + .001) &
           - .03456 *  climate%drhum + .0338 * climate%dtemp_max  + &
           .0234 *  climate%du10_max )

     ! Nesterov Index

     WHERE ( climate%dprecip .GE. 3. .OR. (climate%dtemp_max -climate%dtemp_min) .LT. 4.) 
        climate%Nesterov_Current = 0
     ELSEWHERE
        climate%Nesterov_Current = climate%Nesterov_Current + &
             ( climate%dtemp_max -climate%dtemp_min + 4. ) * climate%dtemp_max 
     END WHERE
     climate%Nesterov_ann_max = MAX(climate%Nesterov_Current,climate%Nesterov_ann_max)
     climate%Nesterov_ann_max = MIN(150000.,climate%Nesterov_ann_max)

     WHERE (climate%NDAY_Nesterov .gt. 365)
        climate%NDAY_Nesterov= 0
        climate%Nesterov_ann_running_max = climate%Nesterov_ann_max       
     ELSEWHERE (climate%NDAY_Nesterov .le. 365)

        WHERE (climate%Nesterov_Current >  climate%Nesterov_ann_running_max)
           climate%Nesterov_ann_running_max = MIN(150000.,climate%Nesterov_Current)
           climate%NDAY_Nesterov = 0
        ELSEWHERE (climate%Nesterov_ann_running_max >= climate%Nesterov_Current)
           climate%NDAY_Nesterov = climate%NDAY_Nesterov + 1
        ENDWHERE
        
     ENDWHERE
  ENDIF
    

!!$     write(3333,"(200f16.6)") real(idoy), real(climate%NDAY_Nesterov(1)), &
!!$          climate%Nesterov_Current(1), &
!!$          climate%Nesterov_ann_max(1), climate%Nesterov_ann_max_last_year(1), &
!!$          climate%Nesterov_ann_running_max(1),  climate%FFDI(1), climate%D_MacArthur(1), &
!!$          climate%KBDI(1), climate%dprecip(1)
  ENDIF

  !  midday fraction of incoming visible radiation absorbed by the canopy
  IF ( mod(ktau,INT(24.0*3600.0/dels)) == INT(24.0*3600.0/dels)/2 &
               ) THEN
    ! climate%fapar_ann_max =   max( (1.-rad%rhocdf(:,1))*(1.-rad%fbeam(:,1)) + &
    !      (1.-rad%rhocbm(:,1))*rad%fbeam(:,1) , climate%fapar_ann_max)
     WHERE (rad%fbeam(:,1).GE.0.01)
        climate%fapar_ann_max = max(1.- rad%extkbm(:,1)*canopy%vlaiw,  &
             climate%fapar_ann_max)
     ENDWHERE
  
  ENDIF

  ! accumulate sub-diurnal sun- and shade-leaf met variables that are relevant for calc of Anet
  climate%APAR_leaf_sun(:,1:nsd-1) =  climate%APAR_leaf_sun(:,2:nsd)
  climate%APAR_leaf_sun(:,nsd) = rad%qcan(:,1,1)*4.6 ! umol m-2 s-1
  climate%APAR_leaf_shade(:,1:nsd-1) = climate%APAR_leaf_shade(:,2:nsd)
  climate%APAR_leaf_shade(:,nsd) = rad%qcan(:,2,1)*4.6 !umol m-2 s-1
  climate%Dleaf_sun(:,1:nsd-1) =  climate%Dleaf_sun(:,2:nsd)
  climate%Dleaf_sun(:,nsd) = canopy%dlf
  climate%fwsoil(:,1:nsd-1) =  climate%fwsoil(:,2:nsd)
  climate%fwsoil(:,nsd) = canopy%fwsoil

  climate%Dleaf_shade(:,1:nsd-1) = climate%Dleaf_shade(:,2:nsd)
  climate%Dleaf_shade(:,nsd) = canopy%dlf

  climate%Tleaf_sun(:,1:nsd-1) = climate%Tleaf_sun(:,2:nsd)
  climate%Tleaf_sun(:,nsd) = canopy%tlf

  climate%Tleaf_shade(:,1:nsd-1) = climate%Tleaf_shade(:,2:nsd)
  climate%Tleaf_shade(:,nsd) = canopy%tlf

  climate%cs_sun(:,1:nsd-1) = climate%cs_sun(:,2:nsd)
  climate%cs_sun(:,nsd) = canopy%cs_sl ! ppm

  climate%cs_shade(:,1:nsd-1) = climate%cs_shade(:,2:nsd)
  climate%cs_shade(:,nsd) = canopy%cs_sh ! ppm

  climate%scalex_sun(:,1:nsd-1) = climate%scalex_sun(:,2:nsd)
  climate%scalex_sun(:,nsd) = rad%scalex(:,1)

  climate%scalex_shade(:,1:nsd-1) = climate%scalex_shade(:,2:nsd)
  climate%scalex_shade(:,nsd) = rad%scalex(:,2)

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

     if (idoy == sum(MonthDays(1:nmonth))) IsLastDay = .TRUE.

     ! On first day of year ...
     IF (idoy==1) THEN

        ! ... reset annual GDD5 counter
        climate%agdd5=0.0
        climate%agdd0=0.0
        climate%evap_PT = 0     ! annual PT evap [mm]
        climate%aevap  = 0      ! annual evap [mm]  

        ! ... reset annual min and max soil moisture
        climate%dmoist_min = climate%dmoist
        climate%dmoist_max = climate%dmoist

        ! ... reset annual precip
        climate%aprecip = 0.0

        ! ...reset annual max nesterov index
        climate%Nesterov_ann_max = 0.0
        climate%fapar_ann_max = 0.0

     ENDIF

     WHERE ((patch%latitude>=0.0 .and. idoy==COLDEST_DAY_NHEMISPHERE).OR. &
          (patch%latitude<0.0 .and. idoy==COLDEST_DAY_SHEMISPHERE) )
        ! In midwinter, reset GDD counter for summergreen phenology
        !climate%gdd5=0.0
        climate%gdd0=0.0
   
        ! reset day degree sum related to spring photosynthetic recovery
        climate%gdd0_rec = 0.0
      END WHERE

     WHERE ((patch%latitude<=0.0 .and. idoy==COLDEST_DAY_NHEMISPHERE).OR. &
          (patch%latitude>0.0 .and. idoy==COLDEST_DAY_SHEMISPHERE) )

        ! In mid-summer, reset dormancy fraction
        climate%fdorm = 1.0

     ENDWHERE

     ! Update GDD counters and chill day count
     climate%gdd0 = climate%gdd0 + max(0.0,climate%dtemp-0.0)
     climate%agdd0= climate%agdd0 + max(0.0,climate%dtemp-0.0)

     climate%gdd5 = climate%gdd5 + max(0.0,climate%dtemp-5.0)
     climate%agdd5= climate%agdd5 + max(0.0,climate%dtemp-5.0)

     ! Update min and max daily soil moisture
     climate%dmoist_min = min(climate%dmoist, climate%dmoist_min)
     climate%dmoist_max = max(climate%dmoist, climate%dmoist_max)

     ! Update annual rainfall total
     climate%aprecip =  climate%aprecip + climate%dprecip

     ! Update dormancy fraction if there has been a frost
     WHERE (climate%dtemp_min .GE. T6 .AND. climate%dtemp_min .LT. 0.0)
        climate%fdorm = max(climate%fdorm - ffrost * climate%dtemp_min/T6 , 0.0)
     ELSEWHERE (climate%dtemp_min .LT. T6)
        climate%fdorm = max(climate%fdorm - ffrost , 0.0)
     ENDWHERE
     
     frec0 = fdorm0 + (1.0 - fdorm0) * climate%fdorm



     WHERE ( climate%dtemp_min .ge. T1 .AND. climate%dtemp_31(:,31) .ge. T1)
        f2 = 1.0
     ELSEWHERE ( climate%dtemp_min .le. T2 .AND. climate%dtemp_31(:,31) .le. T2)
        f2 = 0.0
     ELSEWHERE
        f2 = min( (climate%dtemp_min - T2)/(T1-T2) , &
             ( climate%dtemp_31(:,31) - T2)/(T1-T2))
     ENDWHERE

     WHERE ( climate%dtemp_min .GE. T2)
        f1 = 0.0
     ELSEWHERE ( climate%dtemp_min .LE. T3)
        f1 = 0.3
     ELSEWHERE
        f1 = 0.3*(T2 -  climate%dtemp_min)/ (T2 - T3)
     ENDWHERE

     WHERE (climate%dtemp_min .ge. T2)
       climate%gdd0_rec = max(climate%gdd0_rec + climate%dtemp * f2, 0.0)
     ELSEWHERE (climate%dtemp_min .lt. T2)
        climate%gdd0_rec = max(climate%gdd0_rec*(1. - f1), 0.0)
     ENDWHERE

     WHERE (climate%gdd0_rec .LE. gdd0_rec0)
        climate%frec = frec0 + (1.0 - frec0)* climate%gdd0_rec/gdd0_rec0
     ELSEWHERE
        climate%frec = 1.0
     ENDWHERE

     WHERE (climate%dtemp<5.0 .and. climate%chilldays<=365)
        climate%chilldays = climate%chilldays + 1
     ENDWHERE
     
     ! update GMD (growing moisture day) counter
     where (climate%dmoist .gt. &
          climate%dmoist_min20 + moisture_min*(climate%dmoist_max20 - climate%dmoist_min20))
        climate%gmd = climate%gmd + 1
     elsewhere
        climate%gmd = 0
     endwhere
     
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

     WHERE(mtemp_last>=5.0 .and. climate%mtemp<5.0)
        climate%gdd5=0.0
        climate%chilldays=0
     ENDWHERE

     ! On last day of month ...

     if (IsLastDay) THEN

        ! Update mean temperature for the last 12 months
        ! atemp_mean_new = atemp_mean_old * (11/12) + mtemp * (1/12)

        climate%atemp_mean=climate%atemp_mean*(11./12.)+climate%mtemp*(1./12.)

        ! Record minimum and maximum monthly temperatures

        if (nmonth==1) THEN
           climate%mtemp_min=climate%mtemp;
           climate%mtemp_max=climate%mtemp;
           climate%qtemp_max_last_year =  climate%qtemp_max;
           climate%qtemp_max=climate%qtemp;
        else
           where (climate%mtemp<climate%mtemp_min) &
                climate%mtemp_min=climate%mtemp
           where (climate%mtemp>climate%mtemp_max) &
                climate%mtemp_max=climate%mtemp
           where (climate%qtemp>climate%qtemp_max) &
                climate%qtemp_max=climate%qtemp
        ENDIF  ! first month of year

        ! On 31 December update records of minimum monthly temperatures for the last
        ! 20 years and find minimum monthly temperature for the last 20 years

        if (nmonth==12) THEN
           climate%nyears = climate%nyears +1


           startyear=20-min(19,climate%nyears-1)
           climate%mtemp_min20=0.0
           climate%mtemp_max20=0.0
           climate%alpha_PT20 = 0.0

           climate%dmoist_min20 = 0.0
           climate%dmoist_max20 = 0.0

           climate%aprecip_av20 = 0.0

           climate%fapar_ann_max_last_year =  climate%fapar_ann_max

           climate%Nesterov_ann_max_last_year = climate%Nesterov_ann_max

           if (startyear<20) then
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
                 
                 climate%dmoist_min_20(:,y) = climate%dmoist_min_20(:,y+1)
                 climate%dmoist_max_20(:,y) = climate%dmoist_max_20(:,y+1)
                 climate%dmoist_min20 = climate%dmoist_min20 + climate%dmoist_min_20(:,y)
                 climate%dmoist_max20 = climate%dmoist_max20 + climate%dmoist_max_20(:,y)

                 climate%aprecip_20(:,y) = climate%aprecip_20(:,y+1)
                 climate%aprecip_av20 =  climate%aprecip_av20 + climate%aprecip_20(:,y)
              ENDDO

              climate%mtemp_min20=climate%mtemp_min20/real(20-startyear)
              climate%mtemp_max20=climate%mtemp_max20/real(20-startyear)
              climate%alpha_PT20=climate%alpha_PT20/real(20-startyear)
              climate%dmoist_min20 = climate%dmoist_min20/real(20-startyear)
              climate%dmoist_max20 = climate%dmoist_max20/real(20-startyear)
              climate%aprecip_av20 = climate%aprecip_av20/real(20-startyear)
           else
              ! only occurs when climate%nyears = 1
              climate%mtemp_min20 = climate%mtemp_min
              climate%mtemp_max20 = climate%mtemp_max
              climate%alpha_PT20 = climate%alpha_PT
              climate%dmoist_min20 = climate%dmoist_min
              climate%dmoist_max20 = climate%dmoist_max
              climate%aprecip_av20 = climate%aprecip
           endif

           climate%mtemp_min_20(:,20)=climate%mtemp_min
           climate%mtemp_max_20(:,20)=climate%mtemp_max
           

           climate%alpha_PT = max(climate%aevap/climate%evap_PT, 0.0)     ! ratio of annual evap to annual PT evap
           climate%alpha_PT_20(:,20)=climate%alpha_PT
 
           climate%dmoist_min_20(:,20)= climate%dmoist_min 
           climate%dmoist_max_20(:,20)= climate%dmoist_max

           climate%aprecip_20(:,20) = climate%aprecip

          CALL biome1_pft(climate,np)

        ENDIF  ! last month of year



     ENDIF     ! last day of month


  ENDIF ! end of day
  ! test for seasonal acclimation
  !write(5669,*) climate%qtemp_max_last_year(1), climate%mtemp(1)
  if (cable_user%acclimate_autoresp_seasonal) then
     climate%qtemp_max_last_year =  climate%mtemp
  endif

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
implicit none
real,intent(in):: TC, Pmb       ! temp [deg C], pressure [mb]
real:: Epsif                    ! epsi
real:: TCtmp, ES, dESdT         ! local
real,parameter:: A = 6.106      ! Teten coefficients
real,parameter:: B = 17.27      ! Teten coefficients
real,parameter:: C = 237.3      ! Teten coefficients
real,parameter:: Rlat      = 44140.0  ! lat heat evap H2O at 20C  [J/molW]
real,parameter:: Capp      = 29.09    ! isobaric spec heat air    [J/molA/K]
!-------------------------------------------------------------------------------
TCtmp = TC                          ! preserve TC
if (TCtmp.gt.100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
if (TCtmp.lt.-40.0) TCtmp = -40.0
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
alpha_PT_scaled =  min(climate%alpha_PT20, 1.0)


DO k=1,np

   pft_biome1(k,:) = 999

   IF (climate%mtemp_min20(k) .GE. 15.5) THEN
      IF (alpha_PT_scaled(k).GE.0.85) THEN
         pft_biome1(k,1) = 1
         IF (alpha_PT_scaled(k).LE.0.90) THEN
         !IF (alpha_PT_scaled(k).LE.0.95) THEN
            pft_biome1(k,2) = 2
         ENDIF
      ELSEIF (alpha_PT_scaled(k).GE.0.4 .and. alpha_PT_scaled(k).LT.0.85) THEN
         pft_biome1(k,1) = 2
      ENDIF
   ENDIF

   
   IF (climate%mtemp_min20(k).GE.5 .and.alpha_PT_scaled(k).GE.0.4 &
        .and. pft_biome1(k,1).eq.999 ) THEN
      pft_biome1(k,1) = 3
   ENDIF

   
   IF (climate%mtemp_min20(k).GE.-15 .and. climate%mtemp_min20(k).LE.15.5 .and. &
        alpha_PT_scaled(k).GE.0.35 .and. climate%agdd5(k).gt.1200 & !
        .and. pft_biome1(k,1).gt.3 ) THEN
      pft_biome1(k,1) = 4
   ENDIF
   
   IF (climate%mtemp_min20(k).GE.-19 .and. climate%mtemp_min20(k).LE.5 .and. &
        alpha_PT_scaled(k).GE.0.35 .and. climate%agdd5(k).gt.900 )  THEN
      IF (pft_biome1(k,1).gt.4 ) THEN
         pft_biome1(k,1) = 5
      ELSEIF (pft_biome1(k,1).eq.4 ) THEN
         pft_biome1(k,2) = 5
      ENDIF
   ENDIF
   
   IF (climate%mtemp_min20(k).GE.-35 .and. climate%mtemp_min20(k).LE.-2 .and. &
     alpha_PT_scaled(k).GE.0.35 .and. climate%agdd5(k).gt.550 )  THEN
      IF (pft_biome1(k,1).eq.999 ) THEN
         pft_biome1(k,1) = 6
      ELSEIF (pft_biome1(k,2).eq.999 ) THEN
         pft_biome1(k,2) = 6
      ELSE
         pft_biome1(k,3) = 6
      ENDIF
   ENDIF
   
   IF ( climate%mtemp_min20(k).LE. 5 .and. &
        alpha_PT_scaled(k).GE.0.35 .and. climate%agdd5(k).gt.550 )  THEN
      IF (pft_biome1(k,1).eq.999 ) THEN
         pft_biome1(k,1) = 7
      ELSEIF (pft_biome1(k,2).eq.999 ) THEN
         pft_biome1(k,2) = 7
      ELSEIF (pft_biome1(k,3).eq.999 ) THEN
         pft_biome1(k,3) = 7  
      ELSE
         pft_biome1(k,4) = 7
      ENDIF
   ENDIF
   
   IF (climate%mtemp_min20(k).GE.5 .and.alpha_PT_scaled(k).GE.0.2 &
        .and. pft_biome1(k,1).eq.999 ) THEN
      pft_biome1(k,1) = 8
   ENDIF
   
   IF (climate%mtemp_max20(k).GE.22 .and.alpha_PT_scaled(k).GE.0.1 &
        .and. pft_biome1(k,1).eq.999 ) THEN
      pft_biome1(k,1) = 9
   ENDIF
   
   IF (climate%agdd5(k).GE.500 .and.alpha_PT_scaled(k).GE.0.33 &
   .and. pft_biome1(k,1).eq.999 ) THEN
      pft_biome1(k,1) = 10
   ENDIF
   
   IF (climate%agdd0(k).GE.100 .and.alpha_PT_scaled(k).GE.0.33) THEN 
      IF (pft_biome1(k,1).eq.999 ) THEN
         pft_biome1(k,1) = 11
      ELSEIF (pft_biome1(k,1).eq.10) THEN
         pft_biome1(k,2) = 11
      ENDIF
   ENDIF
   
   IF (climate%mtemp_max20(k).GE.22 .and. pft_biome1(k,1).eq.999 ) THEN
      pft_biome1(k,1) = 12
   ENDIF
   
   IF (climate%agdd0(k).GE.100 .and. pft_biome1(k,1).eq.999 ) THEN
      pft_biome1(k,1) = 13
   ENDIF
   
   ! end of evironmental constraints on pft
  npft = 0
  Do j=1,4
     if (pft_biome1(k,j).NE.999) npft = npft+1
  ENDDO
!     MAP to Biome1 biome and CABLE pft
! (1) Tropical Rainforest
   if(pft_biome1(k,1)==1 .and. npft .eq.1 ) then
       climate%biome(k) = 1
       climate%iveg(k) = 2
   endif

! (2) Tropical Seasonal forest
   if(pft_biome1(k,1)==1 .and.pft_biome1(k,2)==2.and. npft .eq.2 ) then
       climate%biome(k) = 2
       climate%iveg(k) = 2
   endif

! (3) Tropical dry forest/savanna
  if(pft_biome1(k,1)==2.and. npft .eq.1 ) then
       climate%biome(k) = 3
       climate%iveg(k) = 2  ! N.B. need to include c4 grass
  endif


! (4) Broad-leaved evergreen/warm mixed-forest
 if(pft_biome1(k,1)==3.and. npft .eq.1 ) then
       climate%biome(k) = 4
       climate%iveg(k) = 2
 endif

! (5) Temperate deciduous forest
 if(pft_biome1(k,1)==4.and.pft_biome1(k,2)==5.and. &
      pft_biome1(k,3)==7 .and. npft .eq.3 ) then
       climate%biome(k) = 5
       climate%iveg(k) = 4
 endif

! (6) Cool mixed forest
if(pft_biome1(k,1)==4.and.pft_biome1(k,2)==5.and. &
     pft_biome1(k,3)==6 .and.  pft_biome1(k,4)==7 &
     .and. npft .eq.4 ) then
       climate%biome(k) = 6
       climate%iveg(k) = 4
 endif

! (7) Cool conifer forest

 if(pft_biome1(k,1)==5.and.pft_biome1(k,2)==6.and. &
      pft_biome1(k,3)==7 .and. npft .eq.3 ) then
       climate%biome(k) = 7
       climate%iveg(k) = 1
 endif
! (8) Taiga
if(pft_biome1(k,1)==6.and.pft_biome1(k,2)==7 .and. npft .eq.2 ) then
       climate%biome(k) = 8
       climate%iveg(k) = 1
 endif

! (9) Cold mixed forest
if(pft_biome1(k,1)==5.and.pft_biome1(k,2)==7 .and. npft .eq.2 ) then
       climate%biome(k) = 9
       climate%iveg(k) = 1
endif

! (10) Cold deciduous forest
if(pft_biome1(k,1)==7 .and. npft .eq.1 ) then
       climate%biome(k) = 10
       climate%iveg(k) = 3
endif

! (11) Xerophytic woods/scrub
if(pft_biome1(k,1)==8 .and. npft .eq.1 ) then
       climate%biome(k) = 11
       climate%iveg(k) = 5
endif

! (12) Warm grass/shrub
if(pft_biome1(k,1)==9 .and. npft .eq.1 ) then
       climate%biome(k) = 12
       climate%iveg(k) = 5  ! include C4 grass tile ?
endif

! (13) Cool grass/shrub
if(pft_biome1(k,1)==10 .and.pft_biome1(k,2)==11 .and.  npft .eq.2 ) then
       climate%biome(k) = 13
       climate%iveg(k) = 5  ! include C3 grass tile ?
endif

! (14) Tundra
if(pft_biome1(k,1)==11 .and. npft .eq.1 ) then
       climate%biome(k) = 14
       climate%iveg(k) = 8  !
endif

! (15) Hot desert
if(pft_biome1(k,1)==12 .and. npft .eq.1 ) then
       climate%biome(k) = 15
       climate%iveg(k) = 14  !
endif

! (16) Semidesert
if(pft_biome1(k,1)==13 .and. npft .eq.1 ) then
       climate%biome(k) = 16
       climate%iveg(k) = 5  !
endif

! (17) Ice/polar desert
if(climate%biome(k)==999) then
       climate%biome(k) = 17
       climate%iveg(k) = 17  !
endif

! check for DBL or NEL in SH: set to EBL instead
if ((climate%iveg(k)==1 .OR.climate%iveg(k)==3 .OR. climate%iveg(k)==4) &
     .and. patch(k)%latitude<0) THEN
   climate%iveg(k) = 2
   climate%biome(k) = 4
endif


! check for EBL in temperate South America: set to Warm grass/shrub instead.
if (climate%biome(k)==4 .and. &
     (patch(k)%latitude>=-46.25 .and. patch(k)%latitude<= -23.25 &
      .and. patch(k)%longitude>=-65.25 .and. patch(k)%longitude<=-42.75)) then
      climate%biome(k) = 12
      climate%iveg(k) = 5
endif


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

SUBROUTINE climate_init ( climate,np ,ktauday )
IMPLICIT NONE

TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
INTEGER, INTENT(IN) :: np, ktauday
INTEGER :: d

!CALL alloc_cbm_var(climate,np,ktauday)



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
   climate%gmd = 0
   climate%frec = 1.0
   climate%GDD0_rec = 0.0
   climate%fdorm = 1.0
   climate%APAR_leaf_sun =0.0
   climate%APAR_leaf_shade = 0.0
   climate%Dleaf_sun = 0.0
   climate%Dleaf_shade = 0.0
   climate%Tleaf_sun = 0.0
   climate%Tleaf_shade = 0.0
   climate%cs_sun = 0.0
   climate%cs_shade = 0.0
   climate%scalex_sun = 0.0
   climate%scalex_shade = 0.0
   climate%fwsoil = 0.0
   climate%dmoist = 0.0
   climate%dmoist_min = 0.0
   climate%dmoist_max = 0.0
   climate%dmoist_min20 = 0.0
   climate%dmoist_max20 = 0.0
   climate%dmoist_min_20 = 0.0
   climate%dmoist_max_20 = 0.0
   climate%fapar_ann_max_last_year = 0.0

   climate%modis_igbp = 0
   climate%AvgAnnMaxFAPAR = 0.0
   climate%dtemp_max = 0.0
   climate%drhum = 0.0
   climate%du10_max = 0.0
   climate%dprecip = 0.0
   climate%aprecip = 0.0
   climate%last_precip = 0.0
   climate%KBDI = 0.0
   climate%FFDI = 0.0
   climate%D_MacArthur = 0.0
   climate%Nesterov_Current = 0.0

   climate%aprecip_20 = 0.0

   climate%Nesterov_ann_max_last_year = 0.0
   climate%Nesterov_ann_max = 0.0
   climate%Nesterov_ann_running_max = 0.0
   climate%NDAY_Nesterov= 0
   
!if (.not.cable_user%climate_fromzero) then
!   CALL READ_CLIMATE_RESTART_NC (climate, ktauday )
!endif


END SUBROUTINE climate_init

! ==============================================================================

SUBROUTINE WRITE_CLIMATE_RESTART_NC ( climate, ktauday )

  USE netcdf


  IMPLICIT NONE

  TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables
  INTEGER, INTENT(IN) :: ktauday
  INTEGER*4 :: mp4, nsd
  INTEGER*4, parameter   :: pmp4 =0
  INTEGER, parameter   :: fmp4 = kind(pmp4)
  INTEGER*4   :: STATUS
  INTEGER*4   :: FILE_ID, land_ID, nyear_ID, nday_ID, ndayq_ID, nsd_ID, i
    CHARACTER :: CYEAR*4, FNAME*99,dum*50

  ! 0 dim arrays
  CHARACTER(len=20),DIMENSION(2) :: A0
  ! 1 dim arrays (npt )
  CHARACTER(len=30),DIMENSION(33) :: A1
 ! 1 dim arrays (integer) (npt )
  CHARACTER(len=20),DIMENSION(6) :: AI1
  ! 2 dim arrays (npt,20)
  CHARACTER(len=20),DIMENSION(6) :: A2
  ! 2 dim arrays (npt,31)
  CHARACTER(len=20),DIMENSION(2) :: A3
  ! 2 dim arrays (npt,91)
  CHARACTER(len=20),DIMENSION(1) :: A4
  ! 2 dim arrays (npt,ktauday*14)
  CHARACTER(len=20),DIMENSION(11) :: A5


  INTEGER*4 ::  VID0(SIZE(A0)),VID1(SIZE(A1)),VIDI1(SIZE(AI1)), &
       VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4)), VID5(SIZE(A5))

  nsd = int(ktauday*5,fmp4) ! number of sub-diurnal time-steps (for photosynthesis drivers)
  mp4=int(mp,fmp4)
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
  A1(21) = 'GDD0_rec'
  A1(22) = 'frec'
  A1(23) = 'fdorm'
  A1(24) = 'dmoist_min20'
  A1(25) = 'dmoist_max20'
  A1(26) =  'fapar_ann_max_last_year'

  A1(27) = 'last_precip'
  A1(28) = 'Nesterov_Current'
  A1(29) = 'KBDI'
  A1(30) = 'aprecip_av20'
  A1(31) = 'Nesterov_ann_max'
  A1(32) = 'Nesterov_ann_max_last_year'
  A1(33) = 'Nesterov_ann_running_max'
  

  AI1(1) = 'chilldays'
  AI1(2) = 'iveg'
  AI1(3) = 'biome'
  AI1(4) = 'GMD'

  AI1(5) = 'DSLR'
  AI1(6) = 'NDAY_Nesterov'

  A2(1) = 'mtemp_min_20'
  A2(2) = 'mtemp_max_20'
  A2(3) = 'alpha_PT_20'
  A2(4) = 'dmoist_min_20'
  A2(5) = 'dmoist_max_20'

  A2(6) = 'aprecip_20'

  A3(1) = 'dtemp_31'
  A3(2) = 'dmoist_31'

  A4(1) = 'dtemp_91'

  A5(1) = 'APAR_leaf_sun'
  A5(2) = 'APAR_leaf_shade'
  A5(3) = 'Dleaf_sun'
  A5(4) = 'Dleaf_shade'
  A5(5) = 'Tleaf_sun'
  A5(6) = 'Tleaf_shade'
  A5(7) = 'cs_sun'
  A5(8) = 'cs_shade'
  A5(9) = 'scalex_sun'
  A5(10) = 'scalex_shade'
  A5(11) = 'fwsoil'

!# define UM_BUILD YES
# ifndef UM_BUILD

  ! Get File-Name
  IF (LEN_TRIM(cable_user%climate_restart_out).gt.0) THEN
     fname = TRIM(cable_user%climate_restart_out)
  ELSE


     WRITE(CYEAR, FMT='(I4)') CurYear + 1
     fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
          '_climate_rst.nc'
  ENDIF
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
  ! number of sub-diurnal time-steps (stored for leaf photosynthesis drivers)
   STATUS = NF90_def_dim(FILE_ID, 'nsd' , nsd , nsd_ID)
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

  DO i = 1, SIZE(A5)
     STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT,(/land_ID,nsd_ID/),VID5(i))
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

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(21), climate%GDD0_rec )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(22), climate%frec )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(23), climate%fdorm )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(24), climate%dmoist_min20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(25), climate%dmoist_max20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(26), climate%fapar_ann_max_last_year )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(27), climate%last_precip )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(28), climate%Nesterov_Current )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(29), climate%KBDI )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(30), climate%aprecip_av20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(31), climate%Nesterov_ann_max )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(32), climate%Nesterov_ann_max_last_year )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(33), climate%Nesterov_ann_running_max )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), climate%chilldays )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(2), climate%iveg )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(3), climate%biome )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(4), climate%GMD )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(5), climate%DSLR )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(6), climate%NDAY_Nesterov )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


  
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), climate%mtemp_min_20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), climate%mtemp_max_20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), climate%alpha_PT_20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), climate%dmoist_min_20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(5), climate%dmoist_max_20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(6), climate%aprecip_20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), climate%dtemp_31 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), climate%dmoist_31 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), climate%dtemp_91 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), climate%APAR_leaf_sun )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(2), climate%APAR_leaf_shade )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(3), climate%Dleaf_sun )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(4), climate%Dleaf_shade )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(5), climate%Tleaf_sun )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(6), climate%Tleaf_shade )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(7), climate%cs_sun )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(8), climate%cs_shade )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(9), climate%scalex_sun )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(10), climate%scalex_shade )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID5(11), climate%fwsoil )
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
  INTEGER*4 :: mp4, nsd
  INTEGER*4, parameter   :: pmp4 =0
  INTEGER, parameter   :: fmp4 = kind(pmp4)
  INTEGER*4   :: STATUS
  INTEGER*4   :: FILE_ID, land_ID, nyear_ID, nday_ID,nsd_ID, dID, i, land_dim
  CHARACTER :: CYEAR*4, FNAME*99,dum*50

  ! 0 dim arrays
  CHARACTER(len=20),DIMENSION(2) :: A0
  ! 1 dim arrays (npt )
  CHARACTER(len=30),DIMENSION(33) :: A1
 ! 1 dim arrays (integer) (npt )
  CHARACTER(len=20),DIMENSION(6) :: AI1
  ! 2 dim arrays (npt,20)
  CHARACTER(len=20),DIMENSION(6) :: A2
  ! 2 dim arrays (npt,31)
  CHARACTER(len=20),DIMENSION(2) :: A3
 ! 2 dim arrays (npt,91)
  CHARACTER(len=20),DIMENSION(1) :: A4
  ! 2 dim arrays (npt,ktauday*5)
  CHARACTER(len=20),DIMENSION(11) :: A5

  ! MC - climate variables are in single precision
  ! REAL(r_2), DIMENSION(mp)          :: LAT, LON, TMP
  ! REAL(r_2)                         :: TMP2(mp,20),TMP3(mp,31),TMP4(mp,91)
  ! REAL(r_2), ALLOCATABLE:: TMP5(:,:)
  real, dimension(mp) :: LAT, LON, TMP
  real                :: TMP2(mp,20), TMP3(mp,31), TMP4(mp,91)
  real, allocatable   :: TMP5(:,:)
  INTEGER*4 :: TMPI(mp), TMPI0
  LOGICAL            ::  EXISTFILE

  nsd = int(ktauday*5,fmp4) ! number of sub-diurnal time-steps (for photosynthesis drivers)
  ALLOCATE(TMP5(mp, nsd)) 
  mp4=int(mp,fmp4)
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
  A1(21) = 'GDD0_rec'
  A1(22) = 'frec'
  A1(23) = 'fdorm'
  A1(24) = 'dmoist_min20'
  A1(25) = 'dmoist_max20'
  A1(26) =  'fapar_ann_max_last_year'

  A1(27) = 'last_precip'
  A1(28) = 'Nesterov_Current'
  A1(29) = 'KBDI'
  A1(30) = 'aprecip_av20'
  A1(31) = 'Nesterov_ann_max'
  A1(32) = 'Nesterov_ann_max_last_year'
  A1(33) = 'Nesterov_ann_running_max'
  
  AI1(1) = 'chilldays'
  AI1(2) = 'iveg'
  AI1(3) = 'biome'
  AI1(4) = 'GMD'

  AI1(5) = 'DSLR'
  AI1(6) = 'NDAY_Nesterov'

  A2(1) = 'mtemp_min_20'
  A2(2) = 'mtemp_max_20'
  A2(3) = 'alpha_PT_20'
  A2(4) = 'dmoist_min_20'
  A2(5) = 'dmoist_max_20'

  A2(6) = 'aprecip_20'

  A3(1) = 'dtemp_31'
  A3(2) = 'dmoist_31'
  
  A4(1) = 'dtemp_91'

  A5(1) = 'APAR_leaf_sun'
  A5(2) = 'APAR_leaf_shade'
  A5(3) = 'Dleaf_sun'
  A5(4) = 'Dleaf_shade'
  A5(5) = 'Tleaf_sun'
  A5(6) = 'Tleaf_shade'
  A5(7) = 'cs_sun'
  A5(8) = 'cs_shade'
  A5(9) = 'scalex_sun'
  A5(10) = 'scalex_shade'
  A5(11) = 'fwsoil'

  ! Get File-Name
  IF (LEN_TRIM(cable_user%climate_restart_in).gt.0) THEN
     fname = TRIM(cable_user%climate_restart_in)
  ELSE
     WRITE(CYEAR, FMT='(I4)') CurYear + 1
     fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
          '_climate_rst.nc'
  ENDIF
  INQUIRE( FILE=TRIM( fname ), EXIST=EXISTFILE )

        IF ( .NOT.EXISTFILE) write(*,*) fname, ' does not exist!!'

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
 ! number of sub-diurnal time-steps (stored for leaf photosynthesis drivers)
  STATUS = NF90_INQ_DIMID(FILE_ID, 'nsd'   , dID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


write(*,*) 'patch%latitude',  SIZE(patch%latitude)
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
     print*,  TRIM(A1(i))
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
     CASE ('GDD0_rec'  ) ; climate%GDD0_rec  = TMP
     CASE ('frec'  ) ; climate%frec  = TMP
     CASE ('fdorm'  ) ; climate%fdorm  = TMP
     CASE ('dmoist_min20'  ) ; climate%dmoist_min20  = TMP
     CASE ('dmoist_max20'  ) ; climate%dmoist_max20  = TMP
     CASE ('fapar_ann_max_last_year'  ) ; climate%fapar_ann_max_last_year  = TMP
     CASE ('last_precip' ) ; climate%last_precip = TMP
     CASE ('Nesterov_Current' ); climate%Nesterov_Current = TMP
     CASE ('KBDI' ); climate%KBDI = TMP
     CASE ('aprecip_av20' ); climate%aprecip_av20 = TMP
     CASE ('Nesterov_ann_max' ); climate%Nesterov_ann_max = TMP
     CASE ('Nesterov_ann_max_last_year' ); climate%Nesterov_ann_max_last_year = TMP
     CASE ('Nesterov_ann_running_max' ); climate%Nesterov_ann_running_max = TMP
        
     END SELECT
  END DO

! READ 1-dimensional integer fields
  DO i = 1, SIZE(AI1)

     print*,  TRIM(AI1(i))
     STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMPI )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(AI1(i)))
     CASE ('chilldays'      ) ; climate%chilldays      = TMPI
     CASE ('iveg'      ) ; climate%iveg     = TMPI
     CASE ('biome'      ) ; climate%biome     = TMPI
     CASE ('GMD'      ) ; climate%GMD     = TMPI

     CASE ('DSLR' ) ; climate%DSLR = TMPI
     CASE ('NDAY_Nesterov' ) ; climate%NDAY_Nesterov = TMPI
     END SELECT
  END DO


 ! READ 2-dimensional fields (nyear)
  DO i = 1, SIZE(A2)
     write(*,*)  TRIM(A2(i))
     STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A2(i)))
     CASE ('mtemp_min_20' ) ; climate%mtemp_min_20 = TMP2
     CASE ('mtemp_max_20' ) ; climate%mtemp_max_20 = TMP2
     CASE ('alpha_PT_20' ) ; climate%alpha_PT_20 = TMP2
     CASE ('dmoist_min_20' ) ; climate%dmoist_min_20 = TMP2
     CASE ('dmoist_max_20' ) ; climate%dmoist_max_20 = TMP2
     CASE ('aprecip_20' )    ; climate%aprecip_20 = TMP2
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

! READ 2-dimensional fields (nsd)
  DO i = 1, SIZE(A5)
     STATUS = NF90_INQ_VARID( FILE_ID, A5(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP5 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A5(i)))
     CASE ('APAR_leaf_sun' ) ; climate%APAR_leaf_sun  = TMP5
     CASE ('APAR_leaf_shade' ) ; climate%APAR_leaf_shade  = TMP5
     CASE ('Dleaf_sun' ) ; climate%Dleaf_sun  = TMP5
     CASE ('Dleaf_shade' ) ; climate%Dleaf_shade  = TMP5
     CASE ('Tleaf_sun' ) ; climate%Tleaf_sun  = TMP5
     CASE ('Tleaf_shade' ) ; climate%Tleaf_shade  = TMP5
     CASE ('cs_sun' ) ; climate%cs_sun  = TMP5
     CASE ('cs_shade' ) ; climate%cs_shade  = TMP5
     CASE ('scalex_sun' ) ; climate%scalex_sun  = TMP5
     CASE ('scalex_shade' ) ; climate%scalex_shade  = TMP5
     CASE ('fwsoil' ) ; climate%fwsoil  = TMP5


     END SELECT
  END DO


  ! Close NetCDF file:
  STATUS = NF90_close(FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
# endif 


END SUBROUTINE  READ_CLIMATE_RESTART_NC


END MODULE cable_climate_mod
