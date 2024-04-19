! ==============================================================================
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
module cable_climate_mod

  use cable_def_types_mod,  only: met_type, climate_type, canopy_type,soil_snow_type, mp, &
       r_2, radiation_type, veg_parameter_type
  use TypeDef,              only: i4b, dp
  use cable_IO_vars_module, only: patch
  use CABLE_COMMON_MODULE,  only: CurYear, filename, cable_user, HANDLE_ERR

contains

  ! ------------------------------------------------------------------

  subroutine cable_climate(ktau, kstart, ktauday, idoy, LOY, &
       met, climate, canopy, veg, ssnow, rad, &
       dels, np)

    implicit none

    integer,                  intent(IN)    :: ktau      ! integration step number
    integer,                  intent(IN)    :: kstart    ! starting value of ktau
    integer,                  intent(IN)    :: ktauday
    integer,                  intent(IN)    :: idoy, LOY ! day of year (1-365) , Length oy
    type(met_type),           intent(IN)    :: met       ! met input variables
    type(climate_type),       intent(INOUT) :: climate   ! climate variables
    type(canopy_type),        intent(IN)    :: canopy    ! vegetation variables
    type(veg_parameter_type), intent(IN)    :: veg       ! vegetation parameters
    type(soil_snow_type),     intent(IN)    :: ssnow
    type(radiation_type),     intent(IN)    :: rad       ! radiation variables
    real,                     intent(IN)    :: dels      ! integration time setp (s)
    integer,                  intent(IN)    :: np

    integer :: d, y
    integer, parameter :: COLDEST_DAY_NHEMISPHERE = 355
    integer, parameter :: COLDEST_DAY_SHEMISPHERE = 172
    real, parameter :: CoeffPT = 1.26
    real, dimension(mp) :: mtemp_last, phiEq, ppc, EpsA, RhoA
    integer :: startyear
    integer :: MonthDays(12)
    integer :: nmonth, tmp, nsd
    logical :: IsLastDay ! last day of month?
    real, parameter :: Gaero = 0.015  ! (m s-1) aerodynmaic conductance (for use in PT evap)
    real, parameter :: Capp   = 29.09    ! isobaric spec heat air    [J/molA/K]
    real, parameter :: SBoltz  = 5.67e-8  ! Stefan-Boltzmann constant [W/m2/K4]
    ! threshold for setting "growing moisture days", as required for drought-deciduous phenology
    real, parameter :: moisture_min = 0.15
    real, parameter :: T1 = 0.0, T2 = -3.0, T3 = -4.0, T6 = -5.0 ! for computing fractional spring recovery
    real, parameter :: ffrost = 0.1, fdorm0 = 0.15  ! for computing fractional spring recovery
    real, parameter :: gdd0_rec0 = 250.0
    real, dimension(mp) :: f1, f2, frec0, dkbdi

    nsd = ktauday*5 ! number of subdirunal time-steps
    ! to be accumulated for storing variables needed to implement
    ! coordination of photosynthesis

    climate%doy = idoy

    ! accumulate annual evaporation and potential evaporation
    RhoA = met%pmb * 100.0 / (8.314 * (met%tk)) ! air density [molA/m3]
    PPc  = Gaero / ( Gaero + 4.0*SBoltz*(met%tk**3)/(RhoA*Capp) )

    EpSA  = epsif(met%tk - 273.16, met%pmb)
    phiEq = canopy%rniso * (PPc*EpsA) / (PPc*EpsA + 1.0)      ! equil ltnt heat flux  [W/m2]

    if ((idoy == 1) .and. (mod(ktau, ktauday) == 1)) then
       ! first time step of year
       climate%evap_PT       = phiEq * CoeffPT / 2.5014e6 * dels  ! mm
       climate%aevap         = met%precip ! mm
       climate%fapar_ann_max = 0.0
    else
       climate%evap_PT = climate%evap_PT + phiEq * CoeffPT / 2.5014e6 * dels  ! mm
       ! climate%evap_PT = climate%evap_PT + max(phiEq,1.0)*CoeffPT/air%rlam*dels  ! mm
       ! climate%evap_PT =climate%evap_PT + canopy%epot  ! mm
       ! climate%aevap =  climate%aevap + canopy%fe/air%rlam*dels ! mm
       climate%aevap = climate%aevap + met%precip ! mm
    endif

    ! accumulate daily temperature, evap and potential evap
    if (mod(ktau, ktauday) == 1) then
       climate%dtemp     = met%tk - 273.15
       climate%dmoist    = sum(real(ssnow%wb(:,:))*veg%froot(:,:), 2)
       climate%dtemp_min = climate%dtemp
       climate%dtemp_max = climate%dtemp
       climate%drhum     = met%rhum
       climate%du10_max  = met%u10
       climate%dprecip   = met%precip
    else
       climate%dtemp     = climate%dtemp + met%tk - 273.15
       climate%dmoist    = climate%dmoist + sum(real(ssnow%wb(:,:))*veg%froot(:,:), 2)
       climate%dtemp_min = min(met%tk - 273.15, climate%dtemp_min)
       climate%dtemp_max = max(met%tk - 273.15, climate%dtemp_max)
       climate%drhum     = climate%drhum + met%rhum
       climate%du10_max  = max(met%u10, climate%du10_max)
       climate%dprecip   = climate%dprecip + met%precip
    endif

    if (mod((ktau-kstart+1), ktauday) == 0) then  ! end of day
       ! compute daily averages
       climate%dtemp  = climate%dtemp  / real(ktauday)
       climate%dmoist = climate%dmoist / real(ktauday)
       climate%drhum  = climate%drhum  / real(ktauday)

       if (cable_user%CALL_BLAZE) then
          ! update days since last rain and precip since last day without rain
          where (climate%dprecip .gt. 0.01)
             where (climate%DSLR .gt. 0)
                climate%last_precip = climate%dprecip
             elsewhere
                climate%last_precip = climate%last_precip + climate%dprecip
             end where
             climate%DSLR = 0  ! reset days since last rain
          elsewhere
             climate%DSLR =climate%DSLR + 1
          end where
          ! calculate Keetch-Byram Drought Index
          where (climate%DSLR == 0)
             where (climate%last_precip > 5.)
                dkbdi = 5. - climate%last_precip
             elsewhere
                dkbdi = 0.
             end where
          elsewhere
             dkbdi = ((800. - climate%KBDI) * (0.968 * exp(0.0486 * &
                  (climate%dtemp_max * 9./5. + 32.)) &
                  - 8.3) / 1000. / (1. + 10.88 * &
                  exp(-0.0441 * climate%aprecip_av20 /25.4)) * 0.254)
          ENDWHERE
          climate%KBDI = max(0.0, dkbdi + climate%KBDI)

          ! calcululate MacArthur Drought-Factor D
          climate%D_MacArthur = 0.191 * (climate%KBDI + 104.) * &
               (real(climate%DSLR) + 1.)**1.5 / &
               ( 3.52 * (real(climate%DSLR) + 1.)**1.5 &
               + climate%last_precip - 1. )
          climate%D_MacArthur =  max(0., min(10., climate%D_MacArthur))

          ! MacArthur FFDI
          climate%FFDI = 2. * exp( -0.45 + 0.987 * log(climate%D_MacArthur + 0.001) &
               - 0.03456 *  climate%drhum + 0.0338 * climate%dtemp_max  + &
               0.0234 *  climate%du10_max )

          ! Nesterov Index
          where ( climate%dprecip .ge. 3. .or. (climate%dtemp_max -climate%dtemp_min) .lt. 4.)
             climate%Nesterov_Current = 0.
          elsewhere
             climate%Nesterov_Current = climate%Nesterov_Current + &
                  ( climate%dtemp_max -climate%dtemp_min + 4. ) * climate%dtemp_max
          end where
          climate%Nesterov_ann_max = max(climate%Nesterov_Current,climate%Nesterov_ann_max)
          climate%Nesterov_ann_max = min(150000.,climate%Nesterov_ann_max)

          where (climate%NDAY_Nesterov .gt. 365)
             climate%NDAY_Nesterov= 0
             climate%Nesterov_ann_running_max = climate%Nesterov_ann_max
          elsewhere (climate%NDAY_Nesterov .le. 365)
             where (climate%Nesterov_Current >  climate%Nesterov_ann_running_max)
                climate%Nesterov_ann_running_max = min(150000.,climate%Nesterov_Current)
                climate%NDAY_Nesterov = 0
             elsewhere (climate%Nesterov_ann_running_max >= climate%Nesterov_Current)
                climate%NDAY_Nesterov = climate%NDAY_Nesterov + 1
             ENDWHERE
          ENDWHERE
       endif ! call blaze

       ! write(3333,"(200f16.6)") real(idoy), real(climate%NDAY_Nesterov(1)), &
       !      climate%Nesterov_Current(1), &
       !      climate%Nesterov_ann_max(1), climate%Nesterov_ann_max_last_year(1), &
       !      climate%Nesterov_ann_running_max(1),  climate%FFDI(1), climate%D_MacArthur(1), &
       !      climate%KBDI(1), climate%dprecip(1)
    endif

    !  midday fraction of incoming visible radiation absorbed by the canopy
    if (mod(ktau,int(24.0*3600.0/dels)) == int(24.0*3600.0/dels)/2) then
       ! climate%fapar_ann_max =   max( (1.-rad%rhocdf(:,1))*(1.-rad%fbeam(:,1)) + &
       !      (1.-rad%rhocbm(:,1))*rad%fbeam(:,1) , climate%fapar_ann_max)
       where (rad%fbeam(:,1).ge.0.01)
          climate%fapar_ann_max = max(1.- rad%extkbm(:,1)*canopy%vlaiw,  &
               climate%fapar_ann_max)
       ENDWHERE
    endif

    ! accumulate sub-diurnal sun- and shade-leaf met variables that are relevant for calc of Anet
    climate%APAR_leaf_sun(:,1:nsd-1)   = climate%APAR_leaf_sun(:,2:nsd)
    climate%APAR_leaf_sun(:,nsd)       = rad%qcan(:,1,1)*4.6 ! umol m-2 s-1
    climate%APAR_leaf_shade(:,1:nsd-1) = climate%APAR_leaf_shade(:,2:nsd)
    climate%APAR_leaf_shade(:,nsd)     = rad%qcan(:,2,1)*4.6 !umol m-2 s-1
    climate%Dleaf_sun(:,1:nsd-1)       = climate%Dleaf_sun(:,2:nsd)
    climate%Dleaf_sun(:,nsd)           = real(canopy%dlf)
    climate%fwsoil(:,1:nsd-1)          = climate%fwsoil(:,2:nsd)
    climate%fwsoil(:,nsd)              = real(canopy%fwsoil)

    climate%Dleaf_shade(:,1:nsd-1) = climate%Dleaf_shade(:,2:nsd)
    climate%Dleaf_shade(:,nsd)     = real(canopy%dlf)

    climate%Tleaf_sun(:,1:nsd-1) = climate%Tleaf_sun(:,2:nsd)
    climate%Tleaf_sun(:,nsd)     = real(canopy%tlf)

    climate%Tleaf_shade(:,1:nsd-1) = climate%Tleaf_shade(:,2:nsd)
    climate%Tleaf_shade(:,nsd)     = real(canopy%tlf)

    climate%cs_sun(:,1:nsd-1) = climate%cs_sun(:,2:nsd)
    climate%cs_sun(:,nsd)     = real(canopy%cs_sl) ! ppm

    climate%cs_shade(:,1:nsd-1) = climate%cs_shade(:,2:nsd)
    climate%cs_shade(:,nsd)     = real(canopy%cs_sh) ! ppm

    climate%scalex_sun(:,1:nsd-1) = climate%scalex_sun(:,2:nsd)
    climate%scalex_sun(:,nsd)     = rad%scalex(:,1)

    climate%scalex_shade(:,1:nsd-1) = climate%scalex_shade(:,2:nsd)
    climate%scalex_shade(:,nsd)     = rad%scalex(:,2)

    if (mod((ktau-kstart+1), ktauday) == 0) then  ! end of day
       ! get month and check if end of month
       IsLastDay = .false.
       MonthDays = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
       if(LOY==366) MonthDays(2) = MonthDays(2) + 1
       nmonth = 1
       tmp = MonthDays(1)
       do while(tmp.lt.idoy)
          tmp = tmp +  MonthDays(nmonth+1)
          nmonth = nmonth + 1
       enddo

       if (idoy == sum(MonthDays(1:nmonth))) IsLastDay = .true.

       ! On first day of year ...
       if (idoy==1) then
          ! ... reset annual GDD5 counter
          climate%agdd5   = 0.0
          climate%agdd0   = 0.0
          climate%evap_PT = 0.0 ! annual PT evap [mm]
          climate%aevap   = 0.0 ! annual evap [mm]
          ! ... reset annual min and max soil moisture
          climate%dmoist_min = climate%dmoist
          climate%dmoist_max = climate%dmoist
          ! ... reset annual precip
          climate%aprecip = 0.0
          ! ...reset annual max nesterov index
          climate%Nesterov_ann_max = 0.0
          climate%fapar_ann_max = 0.0
       endif

       where ( ((patch%latitude>=0.0) .and. (idoy==COLDEST_DAY_NHEMISPHERE)) .or. &
            ((patch%latitude<0.0) .and. (idoy==COLDEST_DAY_SHEMISPHERE)) )
          ! In midwinter, reset GDD counter for summergreen phenology
          ! climate%gdd5=0.0
          climate%gdd0 = 0.0
          ! reset day degree sum related to spring photosynthetic recovery
          climate%gdd0_rec = 0.0
       end where

       where ( ((patch%latitude<=0.0) .and. (idoy==COLDEST_DAY_NHEMISPHERE)) .or. &
            ((patch%latitude>0.0) .and. (idoy==COLDEST_DAY_SHEMISPHERE)) )
          ! In mid-summer, reset dormancy fraction
          climate%fdorm = 1.0
       ENDWHERE

       ! Update GDD counters and chill day count
       climate%gdd0  = climate%gdd0  + max(0.0, climate%dtemp-0.0)
       climate%agdd0 = climate%agdd0 + max(0.0, climate%dtemp-0.0)

       climate%gdd5  = climate%gdd5  + max(0.0, climate%dtemp-5.0)
       climate%agdd5 = climate%agdd5 + max(0.0, climate%dtemp-5.0)

       ! Update min and max daily soil moisture
       climate%dmoist_min = min(climate%dmoist, climate%dmoist_min)
       climate%dmoist_max = max(climate%dmoist, climate%dmoist_max)

       ! Update annual rainfall total
       climate%aprecip =  climate%aprecip + climate%dprecip

       ! Update dormancy fraction if there has been a frost
       where ((climate%dtemp_min .ge. T6) .and. (climate%dtemp_min .lt. 0.0))
          climate%fdorm = max(climate%fdorm - ffrost * climate%dtemp_min/T6, 0.0)
       elsewhere (climate%dtemp_min .lt. T6)
          climate%fdorm = max(climate%fdorm - ffrost, 0.0)
       ENDWHERE

       frec0 = fdorm0 + (1.0 - fdorm0) * climate%fdorm

       where ((climate%dtemp_min .ge. T1) .and. (climate%dtemp_31(:,31) .ge. T1))
          f2 = 1.0
       elsewhere ((climate%dtemp_min .le. T2) .and. (climate%dtemp_31(:,31) .le. T2))
          f2 = 0.0
       elsewhere
          f2 = min((climate%dtemp_min - T2)/(T1-T2) , &
               (climate%dtemp_31(:,31) - T2)/(T1-T2))
       ENDWHERE

       where (climate%dtemp_min .ge. T2)
          f1 = 0.0
       elsewhere (climate%dtemp_min .le. T3)
          f1 = 0.3
       elsewhere
          f1 = 0.3*(T2 -  climate%dtemp_min)/ (T2 - T3)
       ENDWHERE

       where (climate%dtemp_min .ge. T2)
          climate%gdd0_rec = max(climate%gdd0_rec + climate%dtemp * f2, 0.0)
       elsewhere (climate%dtemp_min .lt. T2)
          climate%gdd0_rec = max(climate%gdd0_rec*(1. - f1), 0.0)
       ENDWHERE

       where (climate%gdd0_rec .le. gdd0_rec0)
          climate%frec = frec0 + (1.0 - frec0) * climate%gdd0_rec / gdd0_rec0
       elsewhere
          climate%frec = 1.0
       ENDWHERE

       where ((climate%dtemp<5.0) .and. (climate%chilldays<=365))
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
       mtemp_last = climate%mtemp

       ! Update daily temperatures, and mean overall temperature, for last 31 days

       climate%mtemp  = climate%dtemp
       climate%qtemp  = climate%dtemp
       climate%mmoist = climate%dmoist
       do d=1, 30
          climate%dtemp_31(:, d)  = climate%dtemp_31(:, d+1)
          climate%mtemp           = climate%mtemp + climate%dtemp_31(:, d)
          climate%dmoist_31(:, d) = climate%dmoist_31(:, d+1)
          climate%mmoist          = climate%mmoist + climate%dmoist_31(:, d)
       enddo
       do d=1, 90
          climate%dtemp_91(:, d) = climate%dtemp_91(:, d+1)
          climate%qtemp          = climate%qtemp + climate%dtemp_91(:, d)
       enddo
       climate%dtemp_31(:, 31)  = climate%dtemp
       climate%dtemp_91(:, 91)  = climate%dtemp
       climate%dmoist_31(:, 31) = climate%dmoist
       climate%qtemp  = climate%qtemp / 91.0  ! average temperature over the last quarter
       climate%mtemp  = climate%mtemp / 31.0  ! average temperature over the last month
       climate%mmoist = climate%mmoist / 31.0 ! average moisture index over the last month

       ! Reset GDD and chill day counter if mean monthly temperature falls below base
       ! temperature

       where ((mtemp_last>=5.0) .and. (climate%mtemp<5.0))
          climate%gdd5      = 0.0
          climate%chilldays = 0
       ENDWHERE

       ! On last day of month ...

       if (IsLastDay) then

          ! Update mean temperature for the last 12 months
          ! atemp_mean_new = atemp_mean_old * (11/12) + mtemp * (1/12)

          climate%atemp_mean = climate%atemp_mean*(11./12.) + climate%mtemp*(1./12.)

          ! Record minimum and maximum monthly temperatures
          if (nmonth==1) then
             climate%mtemp_min = climate%mtemp
             climate%mtemp_max = climate%mtemp
             climate%qtemp_max_last_year = climate%qtemp_max
             climate%qtemp_max = climate%qtemp
          else
             where (climate%mtemp < climate%mtemp_min) &
                  climate%mtemp_min = climate%mtemp
             where (climate%mtemp > climate%mtemp_max) &
                  climate%mtemp_max = climate%mtemp
             where (climate%qtemp > climate%qtemp_max) &
                  climate%qtemp_max = climate%qtemp
          endif  ! first month of year

          ! On 31 December update records of minimum monthly temperatures for the last
          ! 20 years and find minimum monthly temperature for the last 20 years

          if (nmonth==12) then
             climate%nyears = climate%nyears + 1

             startyear = 20 - min(19, climate%nyears-1)
             climate%mtemp_min20 = 0.0
             climate%mtemp_max20 = 0.0
             climate%alpha_PT20  = 0.0

             climate%dmoist_min20 = 0.0
             climate%dmoist_max20 = 0.0

             climate%aprecip_av20 = 0.0

             climate%fapar_ann_max_last_year = climate%fapar_ann_max

             climate%Nesterov_ann_max_last_year = climate%Nesterov_ann_max

             if (startyear<20) then
                do y=startyear, 19
                   climate%mtemp_min_20(:, y) = climate%mtemp_min_20(:, y+1)
                   climate%mtemp_min20 = climate%mtemp_min20+ &
                        climate%mtemp_min_20(:, y)
                   climate%mtemp_max_20(:, y) = climate%mtemp_max_20(:, y+1)
                   climate%mtemp_max20 = climate%mtemp_max20 + &
                        climate%mtemp_max_20(:, y)
                   climate%alpha_PT_20(:, y) = climate%alpha_PT_20(:, y+1)
                   climate%alpha_PT20 = climate%alpha_PT20 + &
                        climate%alpha_PT_20(:, y)

                   climate%dmoist_min_20(:, y) = climate%dmoist_min_20(:, y+1)
                   climate%dmoist_max_20(:, y) = climate%dmoist_max_20(:, y+1)
                   climate%dmoist_min20 = climate%dmoist_min20 + climate%dmoist_min_20(:, y)
                   climate%dmoist_max20 = climate%dmoist_max20 + climate%dmoist_max_20(:, y)

                   climate%aprecip_20(:, y) = climate%aprecip_20(:, y+1)
                   climate%aprecip_av20 = climate%aprecip_av20 + climate%aprecip_20(:, y)
                enddo

                climate%mtemp_min20  = climate%mtemp_min20 / real(20-startyear)
                climate%mtemp_max20  = climate%mtemp_max20 / real(20-startyear)
                climate%alpha_PT20   = climate%alpha_PT20 / real(20-startyear)
                climate%dmoist_min20 = climate%dmoist_min20 / real(20-startyear)
                climate%dmoist_max20 = climate%dmoist_max20 / real(20-startyear)
                climate%aprecip_av20 = climate%aprecip_av20 / real(20-startyear)
             else
                ! only occurs when climate%nyears = 1
                climate%mtemp_min20 = climate%mtemp_min
                climate%mtemp_max20 = climate%mtemp_max
                climate%alpha_PT20 = climate%alpha_PT
                climate%dmoist_min20 = climate%dmoist_min
                climate%dmoist_max20 = climate%dmoist_max
                climate%aprecip_av20 = climate%aprecip
             endif

             climate%mtemp_min_20(:, 20)=climate%mtemp_min
             climate%mtemp_max_20(:, 20)=climate%mtemp_max

             climate%alpha_PT = max(climate%aevap/climate%evap_PT, 0.0) ! ratio of annual evap to annual PT evap
             climate%alpha_PT_20(:, 20) = climate%alpha_PT

             climate%dmoist_min_20(:, 20) = climate%dmoist_min
             climate%dmoist_max_20(:, 20) = climate%dmoist_max

             climate%aprecip_20(:, 20) = climate%aprecip

             call biome1_pft(climate,np)

          endif  ! last month of year

       endif     ! last day of month

    endif ! end of day
    ! test for seasonal acclimation
    !write(5669,*) climate%qtemp_max_last_year(1), climate%mtemp(1)
    if (cable_user%acclimate_autoresp_seasonal) then
       climate%qtemp_max_last_year = climate%mtemp
    endif

  end subroutine cable_climate

  ! ------------------------------------------------------------------

  elemental function epsif(TC, Pmb)
    ! ------------------------------------------------------------------------------
    ! At temperature TC [deg C] and pressure Pmb [mb], return
    ! epsi = (RLAM/CAPP) * d(sat spec humidity)/dT [(kg/kg)/K], from Teten formula.
    ! MRR, xx/1987, 27-jan-94
    ! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
    !                 just like intrinsic functions.
    ! MRR, 28-mar-05: use Rlat [J/molW], Capp [J/molA/K] from MODULE Constants,
    !                 to ensure consistency with other uses of Rlat, Capp
    ! MRR, 28-mar-05: Remove dependence of Rlat (latent heat vaporisation of water)
    !                 on temperature, use value at 20 C
    ! ------------------------------------------------------------------------------
    ! USE TypeDef
    ! USE Constants

    implicit none

    real, intent(in) :: TC, Pmb       ! temp [deg C], pressure [mb]
    real :: epsif                     ! epsi

    real :: TCtmp, ES, dESdT         ! local
    real, parameter:: A = 6.106      ! Teten coefficients
    real, parameter:: B = 17.27      ! Teten coefficients
    real, parameter:: C = 237.3      ! Teten coefficients
    real, parameter:: Rlat      = 44140.0  ! lat heat evap H2O at 20C  [J/molW]
    real, parameter:: Capp      = 29.09    ! isobaric spec heat air    [J/molA/K]

    TCtmp = TC                            ! preserve TC
    if (TCtmp .gt. 100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
    if (TCtmp .lt. -40.0) TCtmp = -40.0
    ES    = A*exp(B*TCtmp/(C+TCtmp))    ! sat vapour pressure
    dESdT = ES*B*C/(C+TCtmp)**2         ! d(sat VP)/dT: (mb/K)
    epsif = (Rlat/Capp) * dESdT / Pmb   ! dimensionless (ES/Pmb = molW/molA)

  end function epsif

  ! ------------------------------------------------------------------

  subroutine biome1_pft(climate, np)
    implicit none

    type (climate_type), intent(INOUT)       :: climate  ! climate variables
    integer, intent(IN) :: np
    integer :: k, j, npft
    integer, allocatable :: pft_biome1(:,:)
    real, allocatable:: alpha_PT_scaled(:)

    ! TABLE 1 , Prentice et al. J. Biogeog., 19, 117-134, 1992
    ! pft_biome1: Trees (1)tropical evergreen; (2) tropical raingreen; (3) warm temp evergreen ;
    ! (4) temperate summergreen; (5) cool-temperate conifer; (6) boreal evergreen conifer;
    ! (7) boreal summergreen;
    ! Non-trees: (8) sclerophyll/succulent; (9) warm grass/shrub; (10) cool grass/shrub;
    ! (11) cold grass/shrub; (12) hot desert shrub; (13) cold desert shrub.

    allocate(pft_biome1(np,4))
    allocate(alpha_PT_scaled(np))
    alpha_PT_scaled =  min(climate%alpha_PT20, 1.0)

    do k=1,np

       pft_biome1(k,:) = 999

       if (climate%mtemp_min20(k) .ge. 15.5) then
          if (alpha_PT_scaled(k).ge.0.85) then
             pft_biome1(k,1) = 1
             if (alpha_PT_scaled(k).le.0.90) then
                !IF (alpha_PT_scaled(k).LE.0.95) THEN
                pft_biome1(k,2) = 2
             endif
          elseif (alpha_PT_scaled(k).ge.0.4 .and. alpha_PT_scaled(k).lt.0.85) then
             pft_biome1(k,1) = 2
          endif
       endif


       if (climate%mtemp_min20(k).ge.5 .and.alpha_PT_scaled(k).ge.0.4 &
            .and. pft_biome1(k,1).eq.999 ) then
          pft_biome1(k,1) = 3
       endif


       if (climate%mtemp_min20(k).ge.-15 .and. climate%mtemp_min20(k).le.15.5 .and. &
            alpha_PT_scaled(k).ge.0.35 .and. climate%agdd5(k).gt.1200 & !
            .and. pft_biome1(k,1).gt.3) then
          pft_biome1(k,1) = 4
       endif

       if (climate%mtemp_min20(k).ge.-19 .and. climate%mtemp_min20(k).le.5 .and. &
            alpha_PT_scaled(k).ge.0.35 .and. climate%agdd5(k).gt.900)  then
          if (pft_biome1(k,1).gt.4) then
             pft_biome1(k,1) = 5
          elseif (pft_biome1(k,1).eq.4) then
             pft_biome1(k,2) = 5
          endif
       endif

       if (climate%mtemp_min20(k).ge.-35 .and. climate%mtemp_min20(k).le.-2 .and. &
            alpha_PT_scaled(k).ge.0.35 .and. climate%agdd5(k).gt.550)  then
          if (pft_biome1(k,1).eq.999) then
             pft_biome1(k,1) = 6
          elseif (pft_biome1(k,2).eq.999) then
             pft_biome1(k,2) = 6
          else
             pft_biome1(k,3) = 6
          endif
       endif

       if ( climate%mtemp_min20(k).le. 5 .and. &
            alpha_PT_scaled(k).ge.0.35 .and. climate%agdd5(k).gt.550 )  then
          if (pft_biome1(k,1).eq.999) then
             pft_biome1(k,1) = 7
          elseif (pft_biome1(k,2).eq.999) then
             pft_biome1(k,2) = 7
          elseif (pft_biome1(k,3).eq.999) then
             pft_biome1(k,3) = 7
          else
             pft_biome1(k,4) = 7
          endif
       endif

       if (climate%mtemp_min20(k).ge.5 .and.alpha_PT_scaled(k).ge.0.2 &
            .and. pft_biome1(k,1).eq.999) then
          pft_biome1(k,1) = 8
       endif

       if (climate%mtemp_max20(k).ge.22 .and.alpha_PT_scaled(k).ge.0.1 &
            .and. pft_biome1(k,1).eq.999) then
          pft_biome1(k,1) = 9
       endif

       if (climate%agdd5(k).ge.500 .and.alpha_PT_scaled(k).ge.0.33 &
            .and. pft_biome1(k,1).eq.999) then
          pft_biome1(k,1) = 10
       endif

       if (climate%agdd0(k).ge.100 .and.alpha_PT_scaled(k).ge.0.33) then
          if (pft_biome1(k,1).eq.999) then
             pft_biome1(k,1) = 11
          elseif (pft_biome1(k,1).eq.10) then
             pft_biome1(k,2) = 11
          endif
       endif

       if (climate%mtemp_max20(k).ge.22 .and. pft_biome1(k,1).eq.999) then
          pft_biome1(k,1) = 12
       endif

       if (climate%agdd0(k).ge.100 .and. pft_biome1(k,1).eq.999) then
          pft_biome1(k,1) = 13
       endif

       ! end of evironmental constraints on pft
       npft = 0
       do j=1,4
          if (pft_biome1(k,j).ne.999) npft = npft+1
       enddo
       !     MAP to Biome1 biome and CABLE pft
       ! (1) Tropical Rainforest
       if (pft_biome1(k,1)==1 .and. npft .eq.1) then
          climate%biome(k) = 1
          climate%iveg(k) = 2
       endif

       ! (2) Tropical Seasonal forest
       if (pft_biome1(k,1)==1 .and.pft_biome1(k,2)==2.and. npft .eq.2) then
          climate%biome(k) = 2
          climate%iveg(k) = 2
       endif

       ! (3) Tropical dry forest/savanna
       if (pft_biome1(k,1)==2.and. npft .eq.1) then
          climate%biome(k) = 3
          climate%iveg(k) = 2  ! N.B. need to include c4 grass
       endif


       ! (4) Broad-leaved evergreen/warm mixed-forest
       if (pft_biome1(k,1)==3.and. npft .eq.1) then
          climate%biome(k) = 4
          climate%iveg(k) = 2
       endif

       ! (5) Temperate deciduous forest
       if (pft_biome1(k,1)==4.and.pft_biome1(k,2)==5.and. &
            pft_biome1(k,3)==7 .and. npft .eq.3) then
          climate%biome(k) = 5
          climate%iveg(k) = 4
       endif

       ! (6) Cool mixed forest
       if (pft_biome1(k,1)==4.and.pft_biome1(k,2)==5.and. &
            pft_biome1(k,3)==6 .and.  pft_biome1(k,4)==7 &
            .and. npft .eq.4) then
          climate%biome(k) = 6
          climate%iveg(k) = 4
       endif

       ! (7) Cool conifer forest
       if (pft_biome1(k,1)==5.and.pft_biome1(k,2)==6.and. &
            pft_biome1(k,3)==7 .and. npft .eq.3) then
          climate%biome(k) = 7
          climate%iveg(k) = 1
       endif

       ! (8) Taiga
       if (pft_biome1(k,1)==6.and.pft_biome1(k,2)==7 .and. npft .eq.2) then
          climate%biome(k) = 8
          climate%iveg(k) = 1
       endif

       ! (9) Cold mixed forest
       if (pft_biome1(k,1)==5.and.pft_biome1(k,2)==7 .and. npft .eq.2) then
          climate%biome(k) = 9
          climate%iveg(k) = 1
       endif

       ! (10) Cold deciduous forest
       if (pft_biome1(k,1)==7 .and. npft .eq.1) then
          climate%biome(k) = 10
          climate%iveg(k) = 3
       endif

       ! (11) Xerophytic woods/scrub
       if (pft_biome1(k,1)==8 .and. npft .eq.1) then
          climate%biome(k) = 11
          climate%iveg(k) = 5
       endif

       ! (12) Warm grass/shrub
       if (pft_biome1(k,1)==9 .and. npft .eq.1) then
          climate%biome(k) = 12
          climate%iveg(k) = 5  ! include C4 grass tile ?
       endif

       ! (13) Cool grass/shrub
       if (pft_biome1(k,1)==10 .and.pft_biome1(k,2)==11 .and.  npft .eq.2) then
          climate%biome(k) = 13
          climate%iveg(k) = 5  ! include C3 grass tile ?
       endif

       ! (14) Tundra
       if (pft_biome1(k,1)==11 .and. npft .eq.1) then
          climate%biome(k) = 14
          climate%iveg(k) = 8  !
       endif

       ! (15) Hot desert
       if (pft_biome1(k,1)==12 .and. npft .eq.1) then
          climate%biome(k) = 15
          climate%iveg(k) = 14  !
       endif

       ! (16) Semidesert
       if (pft_biome1(k,1)==13 .and. npft .eq.1) then
          climate%biome(k) = 16
          climate%iveg(k) = 5  !
       endif

       ! (17) Ice/polar desert
       if (climate%biome(k)==999) then
          climate%biome(k) = 17
          climate%iveg(k) = 17  !
       endif

       ! check for DBL or NEL in SH: set to EBL instead
       if ((climate%iveg(k)==1 .or.climate%iveg(k)==3 .or. climate%iveg(k)==4) &
            .and. patch(k)%latitude<0) then
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

    end do

  end subroutine BIOME1_PFT

  ! ------------------------------------------------------------------

  subroutine climate_init(climate)

    implicit none

    type (climate_type), intent(INOUT) :: climate  ! climate variables

    ! CALL alloc_cbm_var(climate,np,ktauday)

    ! Maciej
    !   DO d=1,31
    !climate%dtemp_31(:,d)= climate%dtemp
    ! climate%dmoist_31(:,d)= climate%dmoist
    !      climate%dtemp_31(:,d)= 0
    !      climate%dmoist_31(:,d)= 0
    !
    !   ENDDO

    climate%nyears = 0
    climate%doy    = 1

    climate%chilldays     = 0
    climate%iveg          = 999
    climate%biome         = 999
    climate%GMD           = 0
    climate%modis_igbp    = 0
    climate%DSLR          = 0
    climate%NDAY_Nesterov = 0

    climate%dtemp                      = 0.0
    climate%dmoist                     = 0.0
    climate%dmoist_min                 = 0.0
    climate%dmoist_min20               = 0.0
    climate%dmoist_max                 = 0.0
    climate%dmoist_max20               = 0.0
    climate%mtemp                      = 0.0
    climate%qtemp                      = 0.0
    climate%mmoist                     = 0.0
    climate%mtemp_min                  = 0.0
    climate%mtemp_max                  = 0.0
    climate%qtemp_max                  = 0.0
    climate%qtemp_max_last_year        = 0.0
    climate%mtemp_min20                = 0.0
    climate%mtemp_max20                = 0.0
    climate%atemp_mean                 = 0.0
    climate%AGDD5                      = 0.0
    climate%GDD5                       = 0.0
    climate%AGDD0                      = 0.0
    climate%GDD0                       = 0.0
    climate%alpha_PT                   = 0.0
    climate%evap_PT                    = 0.0
    climate%aevap                      = 0.0
    climate%alpha_PT20                 = 0.0
    climate%GDD0_rec                   = 0.0
    climate%frec                       = 1.0
    climate%dtemp_min                  = 0.0
    climate%fdorm                      = 1.0
    climate%fapar_ann_max              = 0.0
    climate%fapar_ann_max_last_year    = 0.0
    climate%AvgAnnMaxFAPAR             = 0.0
    climate%dtemp_max                  = 0.0
    climate%drhum                      = 0.0
    climate%du10_max                   = 0.0
    climate%dprecip                    = 0.0
    climate%aprecip                    = 0.0
    climate%aprecip_av20               = 0.0
    climate%last_precip                = 0.0
    climate%KBDI                       = 0.0
    climate%FFDI                       = 0.0
    climate%D_MacArthur                = 0.0
    climate%Nesterov_Current           = 0.0
    climate%Nesterov_ann_max           = 0.0
    climate%Nesterov_ann_max_last_year = 0.0
    climate%Nesterov_ann_running_max   = 0.0

    climate%mtemp_min_20    = 0.0
    climate%mtemp_max_20    = 0.0
    climate%dmoist_min_20   = 0.0
    climate%dmoist_max_20   = 0.0
    climate%dtemp_31        = 0.0
    climate%dmoist_31       = 0.0
    climate%alpha_PT_20     = 0.0
    climate%dtemp_91        = 0.0
    climate%APAR_leaf_sun   = 0.0
    climate%APAR_leaf_shade = 0.0
    climate%Dleaf_sun       = 0.0
    climate%Dleaf_shade     = 0.0
    climate%Tleaf_sun       = 0.0
    climate%Tleaf_shade     = 0.0
    climate%cs_sun          = 0.0
    climate%cs_shade        = 0.0
    climate%scalex_sun      = 0.0
    climate%scalex_shade    = 0.0
    climate%fwsoil          = 0.0
    climate%aprecip_20      = 0.0
    climate%Rd_sun          = 0.0
    climate%Rd_shade        = 0.0

    !if (.not.cable_user%climate_fromzero) then
    !   CALL READ_CLIMATE_RESTART_NC (climate, ktauday)
    !endif

  end subroutine climate_init

  ! ------------------------------------------------------------------

  subroutine WRITE_CLIMATE_RESTART_NC(climate)

    use cable_def_types_mod, only: climate_type, write_netcdf_cbm_var

    implicit none

    type(climate_type), intent(in) :: climate  ! climate variables

    character(len=4)  :: cyear
    character(len=99) :: fname

# ifndef UM_BUILD
    ! get file name
    if (len_trim(cable_user%climate_restart_out) > 0) then
       fname = trim(cable_user%climate_restart_out)
    else
       write(cyear, fmt='(I4)') CurYear + 1
       fname = trim(filename%path) // '/' // trim(cable_user%RunIden) // &
            '_climate_rst.nc'
    endif

    call write_netcdf_cbm_var(trim(fname), climate)
# endif

  end subroutine WRITE_CLIMATE_RESTART_NC

  ! ------------------------------------------------------------------

  subroutine READ_CLIMATE_RESTART_NC(climate)

    use cable_def_types_mod, only: climate_type, read_netcdf_cbm_var

    use netcdf

    implicit none

    type(climate_type), intent(inout) :: climate  ! climate variables

    character(len=4)  :: cyear
    character(len=99) :: fname

    ! get file name
    if (len_trim(cable_user%climate_restart_in) > 0) then
       fname = trim(cable_user%climate_restart_in)
    else
       write(cyear, fmt='(I4)') CurYear + 1
       fname = trim(filename%path) // '/' // trim(cable_user%RunIden) // &
            '_climate_rst.nc'
    endif

    ! read netCDF file
    call read_netcdf_cbm_var(trim(fname), climate)

  end subroutine READ_CLIMATE_RESTART_NC

end module cable_climate_mod
