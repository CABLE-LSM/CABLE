! This file contains Fortran90 code for the POP model,
! a stand-alone tree demography and landscape structure module for Earth system models
! 17-01-2014
! Written by Vanessa Haverd, Ben Smith and Lars Nieradzik
! Report Bugs to Vanessa.Haverd@csiro.au


!CITATION
!--------------------------------------------------------
!When referring to this code in publications, please cite:
! Haverd, V., Smith, B., Cook, G., Briggs, P.R., Nieradzik, L., Roxburgh, S.R., Liedloff, A.,
! Meyer, C.P. and Canadell, J.G., 2013.
! A stand-alone tree demography and landscape structure module for Earth system models.
! Geophysical Research Letters, 40: 1-6.


!DISCLAIMER, COPYRIGHT AND LICENCE

!--------------------------------------------------------

! Use of this code is subject to the Legal Notice and Disclaimer at

! http://www.csiro.au/org/LegalNoticeAndDisclaimer.html

! This code is Copyright, CSIRO, 2014.

! This code is made available under the conditions of the Creative Commons

! Attribution-Share Alike 3.0 License:
! http://creativecommons.org/licenses/by-sa/3.0/
!*******************************************************************************

MODULE TypeDef
  !-------------------------------------------------------------------------------
  ! This module explicitly defines the sizes of variable types
  !-------------------------------------------------------------------------------
  IMPLICIT NONE
  SAVE
  ! Define integer kind parameters to accommodate the range of numbers usually
  ! associated with 4, 2, and 1 byte integers.
  INTEGER,PARAMETER :: i4b = SELECTED_INT_KIND(9)
  INTEGER,PARAMETER :: i2b = SELECTED_INT_KIND(4)
  INTEGER,PARAMETER :: i1b = SELECTED_INT_KIND(2)
  ! Define single and double precision real kind parameters:
  ! * Kind(1.0)   defines sp as the machine's default size for single precision
  ! * Kind(1.0d0) defines dp as the machine's default size for double precision
  INTEGER,PARAMETER :: sp  = KIND(1.0)
  INTEGER,PARAMETER :: dp  = KIND(1.0d0)
  ! lgt is set to the default kind required for representing logical values.
  INTEGER,PARAMETER :: lgt = KIND(.TRUE.)

END MODULE TypeDef


!*******************************************************************************
MODULE POP_Constants
  USE TYPEdef, ONLY: dp, i4b


  REAL(dp),PARAMETER:: FULTON_ALPHA= 3.5 ! recruitment scalar alpha in Fulton (1991)
  REAL(dp),PARAMETER:: DENSINDIV_MAX=0.2   !  Maximum density of individuals within a cohort indiv/m2
  REAL(dp),PARAMETER:: DENSINDIV_MIN=1e-9 !
  REAL(dp),PARAMETER:: Kbiometric=50.0 ! Constant in height-diameter relationship
  REAL(dp),PARAMETER:: WD= 300.0 ! Wood density kgC/m3
  REAL(dp),PARAMETER:: GROWTH_EFFICIENCY_MIN=0.008  ! threshold growth efficiency for enhanced mortality (higher value gives higher biomass turnover)
  REAL(dp),PARAMETER:: Pmort=5.0 ! exponent in mortality formula
  REAL(dp),PARAMETER:: MORT_MAX=0.3 ! upper asymptote for enhanced mortality
  REAL(dp),PARAMETER:: THETA_recruit=0.95 ! shape parameter in recruitment equation
  REAL(dp),PARAMETER:: CMASS_STEM_INIT= 1e-4 ! initial biomass kgC/m2
  REAL(dp),PARAMETER:: POWERbiomass=0.67 ! exponent for biomass in proportion to which cohorts preempt resources
  REAL(dp),PARAMETER:: POWERGrowthEfficiency = 0.67
  REAL(dp),PARAMETER:: CrowdingFactor = 0.029
  REAL(dp),PARAMETER:: ALPHA_CPC = 3.5
  REAL(dp),PARAMETER:: k_allom1 = 200.0 ! crown area =  k_allom1 * diam ** k_rp
  REAL(dp),PARAMETER:: k_rp = 1.67  ! constant in crown area relation to tree diameter
  REAL(dp),PARAMETER:: ksapwood = 0.05 ! rate constant for conversion of sapwood to heartwood (y-1)
  REAL(dp),PARAMETER:: Q=7.0 ! governs rate of increase of mortality with age (2=exponential)
  REAL(dp),PARAMETER:: shootfrac = 0.63
  REAL(dp),PARAMETER:: CtoNw = 400
  REAL(dp),PARAMETER::  CtoNl = 60
  REAL(dp),PARAMETER:: CtoNr = 70
  REAL(dp),PARAMETER:: N_EXTENT = 2 ! multiple of crown diameters within which tree competes with other cohorts
  REAL(dp),PARAMETER:: EPS=1e-12
  INTEGER(i4b),PARAMETER :: NLAYER = 1 ! number of vertical veg layers (1 is currently the only option)
  INTEGER(i4b),PARAMETER :: NCOHORT_MAX = 20 ! maximum number of cohorts
  INTEGER(i4b),PARAMETER :: NDISTURB=1 ! number of disturbance regimes (1 (total only)  or 2 (partial and total))
  INTEGER(i4b),PARAMETER :: PATCH_REPS=10 ! higher number reduces 'noise'
  INTEGER(i4b),PARAMETER :: NAGE_MAX = 1 ! number of maxium ages
  INTEGER(i4b),PARAMETER :: PATCH_REPS1=60 ! number of first dist years
  INTEGER(i4b),PARAMETER :: PATCH_REPS2=1 ! number of second dist years
  INTEGER(i4b),PARAMETER :: NPATCH =PATCH_REPS1*PATCH_REPS2
  INTEGER(i4b),PARAMETER :: NPATCH1D= NPATCH
  INTEGER(i4b),PARAMETER :: NPATCH2D= NPATCH
  INTEGER(i4b),PARAMETER ::  HEIGHT_BINS=12 ! number of height categories to keep track of for diagnostics
  REAL(dp),PARAMETER:: BIN_POWER=1.4 ! bins have muscles
  ! Time base factor (to be multiplied by mean dist interval to give TIMEBASE)
  ! for sampling disturbance probabilities from Poisson distribution
  INTEGER(i4b),PARAMETER :: TIMEBASE_FACTOR=50
  REAL(dp),PARAMETER:: PI=3.14159265358979323846264
  INTEGER(i4b),PARAMETER :: ALLOM_SWITCH = 0 ! 0 == default; 1 = top-end allometry (requires precip as input to POPSTEP)
  ! 0 == binnned max height variable; 1 = continuous (needs lots of memory); 2 = binned by integer heights
  INTEGER(i4b),PARAMETER :: MAX_HEIGHT_SWITCH = 2
  INTEGER(i4b),PARAMETER :: RESOURCE_SWITCH = 1 ! 0 = default; 1  fraction net resource uptake
  INTEGER(i4b),PARAMETER :: RECRUIT_SWITCH = 1 ! 0 = default, 1 = Pgap-dependence
  INTEGER(i4b),PARAMETER :: INTERP_SWITCH = 1 ! 0 = sum over weighted patches, 1 = sum over interpolated patches
  INTEGER(i4b),PARAMETER :: SMOOTH_SWITCH = 0 ! smooth disturbance flux
  INTEGER(i4b),PARAMETER :: NYEAR_SMOOTH = 11 ! smoothing window (y)
  INTEGER(i4b),PARAMETER :: NYEAR_HISTORY =  NYEAR_SMOOTH-NYEAR_SMOOTH/2
  INTEGER(i4b),PARAMETER :: AGEMAX = 1000 
END MODULE POP_Constants
!*******************************************************************************
MODULE POP_Types
  USE TYPEdef, ONLY: dp, i4b
  USE POP_Constants, ONLY: NCOHORT_MAX, NLAYER, HEIGHT_BINS, NDISTURB, NPATCH, NPATCH2D, &
       NYEAR_HISTORY, AGEMAX


  TYPE Cohort
     INTEGER(i4b) :: id
     INTEGER(i4b) :: age ! cohort age
     REAL(dp)     :: biomass ! cohort biomass
     REAL(dp)     :: density ! landscape tree density (weighted mean over patches)
     REAL(dp)     :: frac_resource_uptake
     REAL(dp)     :: frac_light_uptake
     REAL(dp)     :: frac_interception
     REAL(dp)     :: frac_respiration
     REAL(dp)     :: frac_NPP
     REAL(dp)     :: respiration_scalar
     REAL(dp)     :: crown_area
     REAL(dp)     :: Pgap
     REAL(dp)     :: height
     REAL(dp)     :: diameter
     REAL(dp)     :: sapwood
     REAL(dp)     :: heartwood
     REAL(dp)     :: sapwood_area
     REAL(dp)     :: basal_area
     REAL(dp)     :: LAI
     REAL(dp)     :: Cleaf
     REAL(dp)     :: Croot
  END TYPE Cohort

  TYPE Layer
     TYPE (Cohort), DIMENSION(NCOHORT_MAX) :: Cohort
     INTEGER(i4b) :: ncohort ! number of cohorts with density >0
     REAL(dp)    :: biomass ! layer biomass
     REAL(dp)    :: density ! layer tree density
     REAL(dp)     :: hmean ! layer mean tree height (weighted mean over patches)
     REAL(dp)     :: hmax  ! layer max tree height
  END TYPE Layer

  TYPE Patch
     TYPE (Layer), DIMENSION(NLAYER) :: Layer
     REAL(dp)     :: factor_recruit
     REAL(dp)     :: pgap
     REAL(dp)     :: lai
     REAL(dp)     :: biomass ! total biomass in patch
     REAL(dp)     :: biomass_old ! total biomass in patch
     REAL(dp)     :: sapwood ! total sapwood biomass in patch
     REAL(dp)     :: heartwood ! total heartwood biomass in patch
     REAL(dp)     :: sapwood_old ! total sapwood biomass in patch
     REAL(dp)     :: sapwood_area  ! total sapwood area in patch
     REAL(dp)     :: sapwood_area_old ! total sapwood area in patch
     REAL(dp)     :: stress_mortality ! biomass lost in each patch due to stress
     REAL(dp)     :: fire_mortality ! biomass lost in each patch due partial fire disturbance
     REAL(dp)     :: cat_mortality ! biomass lost in each patch due partial fire disturbance
     REAL(dp)     :: crowding_mortality ! biomass lost to crowding mortality
     REAL(dp)     :: cpc
     REAL(dp)     :: mortality !
     REAL(dp)     :: sapwood_loss
     REAL(dp)     :: sapwood_area_loss
     REAL(dp)     :: growth ! biomass growth in each patch due to stem increment
     REAL(dp)     :: area_growth ! basal area growth in each patch due to stem increment
     INTEGER(i4b) :: disturbance_interval(NDISTURB)  ! prescribed disturbance(s) interval for this patch
     INTEGER(i4b) :: first_disturbance_year(NDISTURB)
     INTEGER(i4b) :: age(NDISTURB) ! number of years since last disturbance(s)
     INTEGER(i4b) :: id
     REAL(dp)     :: frac_NPP
     REAL(dp)     :: frac_respiration
     REAL(dp)     :: frac_light_uptake
  END TYPE Patch

  TYPE Landscape
     TYPE (Patch), DIMENSION(NPATCH2D) :: patch
     REAL(dp), DIMENSION(NPATCH2D)     :: freq ! patch weighting
     REAL(dp), DIMENSION(NPATCH2D)     :: freq_old ! patch weighting (previous time-step)
     REAL(dp), DIMENSION(NPATCH2D)     :: fire_freq     !
     REAL(dp), DIMENSION(NPATCH2D)     :: fire_freq_old !
     REAL(dp), DIMENSION(NPATCH2D)     :: cat_freq      !
     REAL(dp), DIMENSION(NPATCH2D)     :: cat_freq_old  !
     REAL(dp), DIMENSION(NPATCH2D,NDISTURB)     :: freq_ranked_age_unique ! unique age weighting
     INTEGER(i4b), DIMENSION(NPATCH2D, NDISTURB)     :: ranked_age_unique ! unique age
     INTEGER(i4b), DIMENSION(NDISTURB)     :: n_age ! number of unique ages
     REAL(dp), DIMENSION(NLAYER)     :: biomass ! landscape stem biomass (weighted mean over patches)
     REAL(dp), DIMENSION(NLAYER)     :: density ! landscape tree density (weighted mean over patches)
     REAL(dp), DIMENSION(NLAYER)     :: hmean ! landscape mean treen height (weighted mean over patches)
     REAL(dp), DIMENSION(NLAYER)     :: hmax  ! landscape max tree height
     REAL(dp), DIMENSION(HEIGHT_BINS)     :: cmass_stem_bin ! biomass by height bin
     REAL(dp), DIMENSION(HEIGHT_BINS)     :: densindiv_bin ! density by height bin
     REAL(dp), DIMENSION(HEIGHT_BINS)     :: height_bin ! mean height in each bin
     REAL(dp), DIMENSION(HEIGHT_BINS)     :: diameter_bin ! mean diameter in each bin
     CHARACTER(100), DIMENSION(HEIGHT_BINS) :: bin_labels ! text strings for bin bounds
     REAL(dp) :: cmass_sum ! landscape biomass
     REAL(dp) :: cmass_sum_old ! landscape biomass
     REAL(dp) :: cheartwood_sum ! landscape biomass (heart wood)
     REAL(dp) :: csapwood_sum ! landscape biomass (sap wood)
     REAL(dp) :: csapwood_sum_old ! landscape biomass
     REAL(dp) :: densindiv ! landscape density of individuals
     REAL(dp) :: height_mean
     REAL(dp) :: height_max
     REAL(dp) :: basal_area
     REAL(dp) :: sapwood_loss ! (kg C m-2 y-1) ! total sapwood loss (turnover + mortality)
     REAL(dp) :: sapwood_area_loss ! ( m2/m-2 y-1) saapwood area loss (mortality only)
     REAL(dp) :: stress_mortality ! (kg C m-2 y-1)
     REAL(dp) :: crowding_mortality ! (kg C m-2 y-1)
     REAL(dp) :: fire_mortality ! (kg C m-2 y-1)
     REAL(dp) :: cat_mortality ! (kg C m-2 y-1)
     REAL(dp) :: res_mortality ! (kg C m-2 y-1)
     REAL(dp) :: growth
     REAL(dp) :: area_growth ! m2/ha
     REAL(dp) :: crown_cover
     REAL(dp) :: crown_area
     REAL(dp) :: crown_volume
     REAL(dp) :: sapwood_area
     REAL(dp) :: sapwood_area_old
     REAL(dp) :: Kclump ! clumping factor
     INTEGER(i4b) :: npatch_active
     INTEGER(i4b) :: LU
     REAL(dp) :: smoothing_buffer
     REAL(dp) :: smoothing_buffer_cat
     REAL(dp) :: fire_mortality_smoothed
     REAL(dp) :: cat_mortality_smoothed
     REAL(dp), DIMENSION(NYEAR_HISTORY) :: fire_mortality_history
     REAL(dp), DIMENSION(NYEAR_HISTORY) :: cat_mortality_history
     REAL(dp), DIMENSION(AGEMAX)     :: freq_age ! age weighting (by age in y: 0:AGE_MAX-1)
     REAL(dp), DIMENSION(AGEMAX)     :: biomass_age
  END TYPE Landscape

  TYPE POP_TYPE
     TYPE(Landscape), DIMENSION(:), ALLOCATABLE :: pop_grid
     INTEGER , DIMENSION(:), Allocatable    :: it_pop
     INTEGER :: np
     INTEGER, DIMENSION(:), Allocatable :: Iwood ! , LU
  END TYPE POP_TYPE

END MODULE POP_Types
!*******************************************************************************

MODULE POPModule
  !-------------------------------------------------------------------------------
  ! * This module contains all subroutines for POP calcs at a single time step.
  !-------------------------------------------------------------------------------
  USE TYPEdef, ONLY: sp, i4b
  USE POP_Types
  USE POP_Constants


CONTAINS


  !*******************************************************************************
  SUBROUTINE ZeroPOP(POP,n)
    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER, OPTIONAL, INTENT(IN) ::n
    INTEGER:: g,k,l,c, np,a,b

    IF ( .NOT. ALLOCATED(pop%pop_grid) ) THEN
       WRITE(*,*)" POP not allocated! Abort in ZeroPOP."
       STOP -1
    ENDIF



    np = SIZE(pop%pop_grid)
! optional interger n intended for zeroing secondary forest tiles
    
    IF (PRESENT(n)) THEN
       a = n
       b = n
       !pop%LU(n) = 2
       POP%pop_grid(n)%LU = 2
    ELSE
       a = 1
       b = np
       !pop%LU = 1
       POP%pop_grid%LU = 1
    endif
 

    DO g=a,b
       POP%pop_grid(g)%freq = 0 ! patch weighting
       POP%pop_grid(g)%freq_old = 0 ! patch weighting
       POP%pop_grid(g)%fire_freq = 0
       POP%pop_grid(g)%fire_freq_old = 0
       POP%pop_grid(g)%cat_freq = 0
       POP%pop_grid(g)%cat_freq_old = 0
       POP%pop_grid(g)%biomass = 0 ! landscape stem biomass (weighted mean over patches)
       POP%pop_grid(g)%density= 0  ! landscape tree density (weighted mean over patches)
       POP%pop_grid(g)%hmean= 0  ! landscape mean treen height (weighted mean over patches)
       POP%pop_grid(g)%hmax= 0   ! landscape max tree height
       POP%pop_grid(g)%cmass_stem_bin= 0  ! biomass by height bin
       POP%pop_grid(g)%densindiv_bin= 0  ! density by height bin
       POP%pop_grid(g)%height_bin= 0  ! mean height in each bin
       POP%pop_grid(g)%diameter_bin= 0  ! mean diameter in each bin
       POP%pop_grid(g)%bin_labels= ' '  ! text strings for bin bounds
       POP%pop_grid(g)%cmass_sum= 0  ! landscape biomass
       POP%pop_grid(g)%csapwood_sum= 0  ! landscape biomass
       POP%pop_grid(g)%cmass_sum_old= 0  ! landscape biomass
       POP%pop_grid(g)%csapwood_sum_old= 0  ! landscape biomass
       POP%pop_grid(g)%cheartwood_sum= 0  ! landscape biomass
       POP%pop_grid(g)%densindiv= 0  ! landscape density of individuals
       POP%pop_grid(g)%height_mean= 0
       POP%pop_grid(g)%height_max= 0
       POP%pop_grid(g)%basal_area= 0
       POP%pop_grid(g)%sapwood_loss = 0
       POP%pop_grid(g)%sapwood_area_loss = 0
       POP%pop_grid(g)%crowding_mortality = 0 ! (kg C m-2 y-1)
       POP%pop_grid(g)%stress_mortality = 0 ! (kg C m-2 y-1)
       POP%pop_grid(g)%fire_mortality = 0 ! (kg C m-2 y-1)
       POP%pop_grid(g)%cat_mortality = 0 ! (kg C m-2 y-1)
       POP%pop_grid(g)%res_mortality = 0 ! (kg C m-2 y-1)
       POP%pop_grid(g)%growth= 0
       POP%pop_grid(g)%area_growth= 0
       POP%pop_grid(g)%crown_cover = 0
       POP%pop_grid(g)%crown_area = 0
       POP%pop_grid(g)%crown_volume = 0
       POP%pop_grid(g)%sapwood_area = 0
       POP%pop_grid(g)%sapwood_area_old = 0
       POP%pop_grid(g)%Kclump = 1
       POP%pop_grid(g)%fire_mortality_smoothed = 0
       POP%pop_grid(g)%cat_mortality_smoothed = 0
       POP%pop_grid(g)%smoothing_buffer = 0
       POP%pop_grid(g)%smoothing_buffer_cat = 0
       POP%pop_grid(g)%fire_mortality_history = 0
       POP%pop_grid(g)%cat_mortality_history = 0
       POP%pop_grid(g)%freq_age = 0
       IF (PRESENT(n)) THEN
          POP%pop_grid(g)%freq_age(1) = 1
       ENDIF
       POP%pop_grid(g)%biomass_age = 0

       DO k=1,NPATCH2D
          POP%pop_grid(g)%patch(k)%factor_recruit= 0
          POP%pop_grid(g)%patch(k)%pgap= 0
          POP%pop_grid(g)%patch(k)%lai= 0
          POP%pop_grid(g)%patch(k)%biomass= 0  ! total biomass in patch
          POP%pop_grid(g)%patch(k)%biomass_old= 0
          POP%pop_grid(g)%patch(k)%sapwood = 0  ! total biomass in patch (sapwood)
          POP%pop_grid(g)%patch(k)%heartwood = 0  ! total biomass in patch (heartwood)
          POP%pop_grid(g)%patch(k)%sapwood_old= 0
          POP%pop_grid(g)%patch(k)%sapwood_area = 0
          POP%pop_grid(g)%patch(k)%sapwood_area_old= 0
          POP%pop_grid(g)%patch(k)%stress_mortality = 0  ! biomass lost in each patch due to stress
          POP%pop_grid(g)%patch(k)%fire_mortality= 0 ! biomass lost in each patch due to fire partial dist
          POP%pop_grid(g)%patch(k)%cat_mortality= 0 ! biomass lost in each patch due to fire partial dist
          POP%pop_grid(g)%patch(k)%crowding_mortality = 0
          POP%pop_grid(g)%patch(k)%cpc = 0
          POP%pop_grid(g)%patch(k)%mortality = 0
          POP%pop_grid(g)%patch(k)%sapwood_loss = 0
          POP%pop_grid(g)%patch(k)%sapwood_area_loss = 0
          POP%pop_grid(g)%patch(k)%growth= 0  ! biomass growth in each patch due stem increment
          POP%pop_grid(g)%patch(k)%area_growth= 0
          POP%pop_grid(g)%patch(k)%disturbance_interval= 0   ! prescribed disturbance(s) interval for this patch
          POP%pop_grid(g)%patch(k)%first_disturbance_year = 0
          POP%pop_grid(g)%patch(k)%age= 0  ! number of years since last disturbance(s)
          POP%pop_grid(g)%patch(k)%id = 0
          POP%pop_grid(g)%patch(k)%frac_NPP = 0
          POP%pop_grid(g)%patch(k)%frac_respiration = 0
          POP%pop_grid(g)%patch(k)%frac_light_uptake = 0
          DO l=1,NLAYER
             POP%pop_grid(g)%patch(k)%Layer(L)%ncohort = 0 ! number of cohorts with density >0
             POP%pop_grid(g)%patch(k)%Layer(L)%biomass = 0 ! layer biomass
             POP%pop_grid(g)%patch(k)%Layer(L)%density= 0  ! layer tree density
             POP%pop_grid(g)%patch(k)%Layer(L)%hmean= 0  ! layer mean tree height (weighted mean over patches)
             POP%pop_grid(g)%patch(k)%Layer(L)%hmax= 0   ! layer max tree height
             DO c = 1,NCOHORT_MAX
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%id = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%age = 0 ! cohort age
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%biomass = 0 ! cohort biomass
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%density = 0 ! landscape tree density (weighted mean over patches)
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_resource_uptake = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_light_uptake = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_interception = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_respiration = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_NPP = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%respiration_scalar = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%crown_area = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%Pgap = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%height = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%diameter = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%sapwood = 0 ! cohort sapwood
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%heartwood = 0 ! cohort heartwood
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%sapwood_area = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%basal_area = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%LAI = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%Cleaf = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%Croot = 0
             ENDDO
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE ZeroPOP
  !*******************************************************************************
 SUBROUTINE InitPOP2D_Poisson(POP, mean_disturbance_interval, m)
  ! Initialises vector of patches with maximum age correpondding to 95% of pdf
  ! Starting year: uniform distribution up to maximum age

  IMPLICIT NONE

  TYPE(POP_TYPE), INTENT(INOUT) :: POP
  INTEGER(i4b), INTENT(IN) ::  mean_disturbance_interval(:,:)
  INTEGER(i4b), INTENT(IN), optional :: m 
  INTEGER(i4b) :: j, k, g, ipatch, idist, p, c, n, i
  INTEGER(i4b) :: disturbance_interval
  INTEGER(i4b):: patch_disturbance_interval_idist(NDISTURB,NPATCH2D)
  INTEGER(i4b):: patch_first_disturbance_year_idist(NDISTURB,NPATCH2D), &
       patch_first_disturbance_year_unique(NDISTURB,NPATCH2D)

  INTEGER(i4b):: Poisson_age(1000),Poisson_freq(1000)
  REAL(dp):: Poisson_weight(1000), CumPoisson_weight(1000)
  INTEGER(i4b):: disturbances_per_timebase, timebase
  INTEGER:: i_min, i_max, age_sample(2,NAGE_MAX**2), tmp(NAGE_MAX**2)
  INTEGER:: age_tmp, tmp_unique(NAGE_MAX**2), n_age, np
  REAL(dp):: disturbance_freq, tmp1
  INTEGER:: tmp2(PATCH_REPS1) ,  n_first_disturbance_year(NDISTURB), tmp3(PATCH_REPS2)
  INTEGER:: a,b
  np = SIZE(POP%pop_grid)
  a = 1
  b = np
  IF (PRESENT(m)) THEN
     a = m
     b = m
  ENDIF



  DO g=a,b

     ! calculate Poisson weights for each of the 2 mean disturbance intervals
     IF (NPATCH.gt.1) THEN
        DO idist=1,NDISTURB
           disturbance_freq=1.0/REAL(mean_disturbance_interval(g,idist))
           DO p = 1,1000
              Poisson_age(p) = p
              Poisson_weight(p) = Exponential(disturbance_freq,p)
              CumPoisson_weight(p) = CumExponential(disturbance_freq,REAL(p,dp))
           ENDDO
           ! set max age to correspond to 95% percentile of cum pdf
           DO k =1,NPATCH2D
              i_max = MAXLOC(Poisson_age,1,CumPoisson_weight.LE.0.95)
              POP%pop_grid(g)%patch(k)%disturbance_interval(idist) = Poisson_age(i_max)
              POP%pop_grid(g)%patch(k)%id = k
              POP%pop_grid(g)%patch(k)%age = 0
           ENDDO
        ENDDO
!!

       DO idist =1,ndisturb
        ! set first disturbance year for first dist interval class
        if (idist .eq. 1) then
        disturbance_interval = POP%pop_grid(g)%patch(1)%disturbance_interval(idist)
        DO c = 1,PATCH_REPS1
           if (c==1) then
              tmp2(1) = 1
           else
              tmp2(1) = max(disturbance_interval*(c-1)/(PATCH_REPS1),1)+1
           endif
           tmp2(2) = max(disturbance_interval*c/(PATCH_REPS1),1)
           tmp2(3) = tmp2(1)
       !    write(*,*) 'tmp2', c, disturbance_interval, tmp2(1),tmp2(2)
           DO j = 1,PATCH_REPS2
              ipatch = (c-1)*PATCH_REPS2 + j
              POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(idist) = tmp2(3)
              tmp2(3)=tmp2(3)+1
              if (tmp2(3)>tmp2(2)) then
                 tmp2(3) = tmp2(1)
              ENDIF
           ENDDO
        ENDDO

!!$        ! set first disturbance year for first dist interval class
!!$        idist = 1
!!$        disturbance_interval = POP%pop_grid(g)%patch(1)%disturbance_interval(idist)
!!$        DO c = 1,PATCH_REPS1
!!$           tmp2(c) = max(disturbance_interval*c/(PATCH_REPS1),1)
!!$        ENDDO
!!$        DO c = 1,PATCH_REPS1
!!$           i = 0
!!$           DO j = 1,PATCH_REPS2
!!$              ipatch = (j-1)*PATCH_REPS1 + c
!!$              i = i+1
!!$              IF (i.gt.PATCH_REPS1) then
!!$                 i = 1
!!$              ENDIF
!!$              do while ((tmp2(i+1).eq. tmp2(i)).and.(i.lt.PATCH_REPS1))
!!$                 i = i+1
!!$                 IF (i.gt.PATCH_REPS1) then
!!$                    i = 1
!!$                 ENDIF
!!$
!!$              ENDDO
!!$
!!$
!!$              write(*,*) i, tmp2(i)
!!$              POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(idist) = tmp2(i)
!!$
!!$           ENDDO
!!$        ENDDO


    ! set first disturbance year for first 2nd interval class
        ELSEIF (idist.eq.2) then
           disturbance_interval = POP%pop_grid(g)%patch(1)%disturbance_interval(idist)
           DO c = 1,PATCH_REPS2
              tmp3(c) = max(disturbance_interval*(c-1)/(PATCH_REPS2),1)
           ENDDO

           DO c = 1,PATCH_REPS2
              i = 0
              DO j = 1,PATCH_REPS1
                  ipatch = (j-1)*PATCH_REPS2 + c
                  POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(idist) = &
                       tmp3(c) +(j-1)*max((tmp3(idist)-tmp3(1))/PATCH_REPS1,1)

               !  i = i+1
               !  if (i.gt.(tmp3(2)-tmp3(1))) i = 0
              ENDDO
           ENDDO
        ENDIF




        ENDDO



     ELSE   ! NPATCH =1 (single patch mode)
        k = 1
        DO idist=1,NDISTURB
           POP%pop_grid(g)%patch(k)%disturbance_interval(idist) = mean_disturbance_interval(g,idist)
           POP%pop_grid(g)%patch(k)%first_disturbance_year(idist) = 113
           POP%pop_grid(g)%patch(k)%age = 0
           POP%pop_grid(g)%patch(k)%id = k
        ENDDO
     ENDIF

     POP%pop_grid(g)%npatch_active = NPATCH

  ENDDO

END SUBROUTINE InitPOP2D_Poisson
!
  !*******************************************************************************

  SUBROUTINE POPStep(POP, StemNPP, disturbance_interval, disturbance_intensity,LAI,Cleaf,Croot, &
       NPPtoGPP, StemNPP_av,frac_intensity1,precip )
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    REAL(dp), INTENT(IN) :: StemNPP(:,:)
    REAL(dp), INTENT(IN) :: disturbance_intensity(:,:)
    INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:)
    REAL(dp), INTENT(IN) ::  LAI(:)
    REAL(dp), INTENT(IN) ::  Cleaf(:)
    REAL(dp), INTENT(IN) ::  Croot(:)
    REAL(dp), INTENT(IN) ::  NPPtoGPP(:)
    REAL(dp), INTENT(IN), OPTIONAL :: frac_intensity1(:), precip(:)
    REAL(dp), INTENT(IN), OPTIONAL :: StemNPP_av(:)
    INTEGER(i4b) :: idisturb,np,g
    INTEGER(i4b), allocatable :: it(:)
    REAL(dp):: dallocW

    !INTEGER, INTENT(IN) :: wlogn
    pop%it_pop = pop%it_pop + 1

    !it = pop%it_pop(1)
    np = SIZE(POP%POP_grid)
    allocate(it(np))

    do g=1,np
       it(g) =  maxval(pop%pop_grid(g)%patch(:)%age(1)) + 1
    enddo
!!$    DO idisturb = 1,NDISTURB
!!$       CALL GetUniqueAgeFrequencies(POP, disturbance_interval, idisturb, it)
!!$    ENDDO
!!$
!!$    CALL GetPatchFrequencies(POP,it)

   !call flush(wlogn)
    IF (PRESENT(precip)) THEN
       IF(PRESENT(StemNPP_av)) THEN
          CALL PatchAnnualDynamics(POP, StemNPP,NPPtoGPP,disturbance_interval, it, precip=precip,StemNPP_av=StemNPP_av)
       ELSE
          CALL PatchAnnualDynamics(POP, StemNPP,NPPtoGPP,disturbance_interval, it, precip=precip)
       ENDIF
    ELSE
       IF(PRESENT(StemNPP_av)) THEN
          CALL PatchAnnualDynamics(POP, StemNPP,NPPtoGPP,disturbance_interval, it,StemNPP_av=StemNPP_av)
       ELSE
          CALL PatchAnnualDynamics(POP, StemNPP,NPPtoGPP,disturbance_interval, it)
       ENDIF
    ENDIF

    IF (NDISTURB.EQ.1) THEN
       IF (PRESENT(precip)) THEN
          !   CALL Patch_disturb(POP,it,1,precip)
          CALL Patch_partial_disturb2(POP,1,precip)
       ELSE
          CALL Patch_disturb(POP,1)
          ! CALL Patch_partial_disturb2(POP,it,1)
       ENDIF
    ELSEIF (NDISTURB.EQ.2) THEN
       IF (PRESENT(frac_intensity1)) THEN
          IF (PRESENT(precip)) THEN
             CALL Patch_partial_disturb(POP,1,disturbance_intensity,precip,frac_intensity1=frac_intensity1)
          ELSE
             CALL Patch_partial_disturb(POP,1,disturbance_intensity,frac_intensity1=frac_intensity1)
          ENDIF
       ELSE
          IF (PRESENT(precip)) THEN
             CALL Patch_partial_disturb(POP,1,disturbance_intensity,precip=precip)
          ELSE
             CALL Patch_partial_disturb(POP,1,disturbance_intensity)
          ENDIF
       ENDIF
       IF (PRESENT(precip)) THEN
          !CALL Patch_partial_disturb2(POP,it,2,precip)
          CALL Patch_disturb(POP,2,precip)
       ELSE
          ! CALL Patch_partial_disturb2(POP,it,2)
          CALL Patch_disturb(POP,2)
       ENDIF
    ENDIF

    DO idisturb = 1,NDISTURB
       CALL GetUniqueAgeFrequencies(POP, disturbance_interval, idisturb, it)
    ENDDO

    CALL GetPatchFrequencies(POP,it)

    !PRINT*,"Get Diags"
    IF (PRESENT(precip)) THEN
       CALL GetDiagnostics(pop, LAI,Cleaf,Croot, disturbance_interval, it,precip)
    ELSE
       CALL GetDiagnostics(pop, LAI,Cleaf,Croot, disturbance_interval, it)
    ENDIF



  END SUBROUTINE POPStep
  !*******************************************************************************
  SUBROUTINE PatchAnnualDynamics(pop, StemNPP,NPPtoGPP, disturbance_interval, it, precip,StemNPP_av)
    IMPLICIT NONE

    TYPE( POP_TYPE ), INTENT(INOUT) :: pop
    REAL(dp), INTENT(IN)            :: StemNPP(:,:)
    REAL(dp), INTENT(IN)            :: NPPtoGPP(:)
    INTEGER(i4b), INTENT(IN)        ::  disturbance_interval(:,:)
    REAL(dp), INTENT(IN), OPTIONAL  :: precip(:)
    REAL(dp), OPTIONAL, INTENT(IN)            :: StemNPP_av(:)
    INTEGER(i4b), INTENT(IN)        :: it(:)

    REAL(dp) :: f, mu, densindiv, cmass
    REAL(dp) :: tmp,tmp_light,tmp_respiration,tmp_fracnpp, cmass_stem_sum,cmass_stem_inc
    INTEGER(i4b) :: j, k,c, ncohort, idist
    INTEGER(i4b) :: ivec(NCOHORT_MAX), nc, tmp1(NPATCH2D), np, idisturb
    REAL(dp) :: growth_efficiency,cmass_stem
    REAL(dp) :: mort, mort_bg, fire_mort
    REAL(dp) :: s2, cpc, crown_area
    REAL(dp) :: mort_cpc
    REAL(dp) :: basal, ht, diam, area_growth_grid , basal_grid, basal_new, basal_old
    REAL(dp) :: tmp2(NCOHORT_MAX), freq,tmp3(NPATCH2D),tmp4(NPATCH2D),tmp5(NPATCH2D)
    idisturb = 1
    np = SIZE(POP%POP_grid)

    ! growth
    ! Distributes layer biomass increment among cohorts and increments age
    ! calculate fractional resource uptake by each cohort
    DO j=1,np
       basal_grid = 0.0
       area_growth_grid = 0.0
       pop%pop_grid(j)%sapwood_area_old      = pop%pop_grid(j)%sapwood_area
       pop%pop_grid(j)%freq_old      = pop%pop_grid(j)%freq
       pop%pop_grid(j)%fire_freq_old = pop%pop_grid(j)%fire_freq
       pop%pop_grid(j)%cat_freq_old  = pop%pop_grid(j)%cat_freq

       ! Get fraction allocation for each patch
       tmp = 0.0
       tmp_light = 0.0
       tmp_respiration = 0.0
       tmp_fracNPP = 0.0
  
       if (NPATCH2D >1.and.it(j) > 1.and. RESOURCE_SWITCH>0) then
          DO k=1,NPATCH2D
             freq =  pop%pop_grid(j)%freq(pop%pop_grid(j)%patch(k)%id)
             nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             DO c=1,nc
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_light_uptake = &
                     pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%frac_interception      ! defined in terms of Pgap
                ! total autotrophic resp, summed over all cohorts and patches
                tmp_respiration = tmp_respiration + &
                     freq*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%respiration_scalar
             ENDDO
             tmp_light = tmp_light + freq*(1. - pop%pop_grid(j)%patch(k)%Pgap)
          ENDDO
          IF (tmp_respiration .gt. 1.e-8 .and. tmp_light .gt. 1.e-8) then
             DO k=1,NPATCH2D
                ! fraction respiration and un-normalised NPP for each patch
                nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
                freq =  pop%pop_grid(j)%freq(pop%pop_grid(j)%patch(k)%id)
                ! frac autotrophic resp
              
                pop%pop_grid(j)%patch(k)%frac_respiration = &
                     sum(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%respiration_scalar)/tmp_respiration

                ! frac gpp
                pop%pop_grid(j)%patch(k)%frac_light_uptake = &
                     (1. - pop%pop_grid(j)%patch(k)%pgap) /tmp_light
                ! frac npp
                pop%pop_grid(j)%patch(k)%frac_NPP = &
                     max(pop%pop_grid(j)%patch(k)%frac_light_uptake*(1/NPPtoGPP(j)) - &
                     pop%pop_grid(j)%patch(k)%frac_respiration*(1/NPPtoGPP(j)-1.0),0.0_dp)

                tmp_fracNPP = tmp_fracNPP + freq*pop%pop_grid(j)%patch(k)%frac_NPP


             ENDDO

             ! normalised fraction NPP
             DO k=1,NPATCH2D
                pop%pop_grid(j)%patch(k)%frac_NPP = &
                     pop%pop_grid(j)%patch(k)%frac_NPP/tmp_fracNPP

             ENDDO
          ELSE
             pop%pop_grid(j)%patch(:)%frac_NPP = 1.0
             pop%pop_grid(j)%patch(:)%frac_respiration = 1.0
             pop%pop_grid(j)%patch(:)%frac_light_uptake = 1.0


          ENDIF
       else
          pop%pop_grid(j)%patch(:)%frac_NPP = 1.0
          pop%pop_grid(j)%patch(:)%frac_respiration = 1.0
          pop%pop_grid(j)%patch(:)%frac_light_uptake = 1.0
       endif
       ! End Get fraction allocation for each patch

       ! Get fraction allocation for each cohort in each patch
       DO k=1,NPATCH2D
          tmp = 0.0
          tmp_light = 0.0
          tmp_respiration = 0.0
          tmp_fracNPP = 0.0
          if (pop%pop_grid(j)%patch(k)%Layer(1)%ncohort>1) then
             nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

                cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
                densindiv = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density

                ! get initial basal area

                IF ( PRESENT(precip) ) THEN
                   CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_old, precip(j))
                ELSE
                   CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_old )
                ENDIF

                if (ALLOM_SWITCH.eq.1) then
                   !! assumes crown radius (m) = 0.1492 * dbh (cm) (from G. Cook, pers. comm.)
                   crown_area = densindiv*PI*(diam*100.*0.1492)**2
                else
                   crown_area = densindiv*PI*(((k_allom1 * diam ** k_rp )/PI)**0.5)**2
                endif

                tmp = tmp + (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass/ &  ! sum over all cohorts
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density)**POWERbiomass * &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density

                tmp_light = tmp_light + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception

                tmp_respiration = tmp_respiration + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%respiration_scalar

                tmp2(c) = sum((pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c:nc)%biomass/ &  ! sum over all cohorts c:nc
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c:nc)%density)**POWERbiomass * &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c:nc)%density)

             ENDDO

             ! un-normalised fractional gross resource uptake: weighted combination of components
             ! where cohort competes with older cohorts and where it does not
             DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

                if (RESOURCE_SWITCH ==1) then
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception = &
                        pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception/tmp_light
                else
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception = 1.0
                endif

             ENDDO

             !normalised fractional gross resource uptake
             nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
                !normalised fractional gross resource uptake
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_light_uptake = &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception/ &
                     sum(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%frac_interception)
             ENDDO


             ! fraction respiration and un-normalised NPP
             DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_respiration = &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%respiration_scalar/tmp_respiration


                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP = &
                     max(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_light_uptake*(1/NPPtoGPP(j)) - &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_respiration*(1/NPPtoGPP(j)-1.0),0.0_dp)


                tmp_fracNPP = tmp_fracNPP +  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP

             ENDDO

             ! normalised fraction NPP
             DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP = &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP/tmp_fracNPP

             ENDDO


             ! fraction net resource uptake

             DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

                if (RESOURCE_SWITCH==0) then

                   ! default net fraction resource uptake
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake = &
                        (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass/ &
                        pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density)**POWERbiomass * &
                        pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density/tmp


                elseif (RESOURCE_SWITCH==1) then

                   ! fraction net resource uptake = fraction NPP

                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake = &
                        pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP * &
                        pop%pop_grid(j)%patch(k)%frac_NPP

                endif

             ENDDO

          else
             c = 1
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP = 1
             tmp_fracNPP = 1
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_respiration = 1
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_light_uptake =1
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake = 1

             if (RESOURCE_SWITCH==1) then
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake = pop%pop_grid(j)%patch(k)%frac_NPP
             endif

          endif


       ENDDO

       tmp = 0
       DO k=1,NPATCH2D
          pop%pop_grid(j)%patch(k)%sapwood_loss = 0.0
          pop%pop_grid(j)%patch(k)%sapwood_area_loss = 0.0
          pop%pop_grid(j)%patch(k)%sapwood_old = pop%pop_grid(j)%patch(k)%sapwood
          pop%pop_grid(j)%patch(k)%sapwood_area_old = pop%pop_grid(j)%patch(k)%sapwood_area
          pop%pop_grid(j)%patch(k)%biomass_old = pop%pop_grid(j)%patch(k)%biomass
          pop%pop_grid(j)%patch(k)%growth = 0.0
          pop%pop_grid(j)%patch(k)%area_growth = 0.0
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          freq =  pop%pop_grid(j)%freq(pop%pop_grid(j)%patch(k)%id)
          DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

             cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
             densindiv = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density

             ! get initial basal area

             IF ( PRESENT(precip) ) THEN
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_old, precip(j))
             ELSE
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_old )
             ENDIF

             ! increment biomass in cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass +  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake*StemNPP(j,1)

             cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
             tmp = tmp + freq*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake

             ! get incremented basal area
             IF ( PRESENT(precip) ) THEN
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_new, precip(j))
             ELSE
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_new )
             ENDIF


             ! increment sapwood in cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood +  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake*StemNPP(j,1)


             ! increment heartwood in cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood + &
                  ksapwood*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood

             ! keep track of patch-level sapwood loss
             pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                  ksapwood*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood

             ! decrease sapwood in cohort (accounting for loss to heartwood)
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood =  &
                  (1. - ksapwood)*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood


             !if ( pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.gt.1e-9) then
             ! patch biomass increment
             pop%pop_grid(j)%patch(k)%growth = pop%pop_grid(j)%patch(k)%growth +  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake*StemNPP(j,1)

             ! patch sapwood area increment
             pop%pop_grid(j)%patch(k)%area_growth = pop%pop_grid(j)%patch(k)%area_growth +  &
                  basal_new-basal_old


             !  endif
             ! increment cohort age
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%age = &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%age + 1

          ENDDO
          ! Layer biomass (summed over cohorts)
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          pop%pop_grid(j)%patch(k)%Layer(1)%biomass = SUM(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%biomass)


       ENDDO


    ENDDO

    ! Mortality
    !Implements resource stress mortality and crowding mortality for all cohorts in layer

    DO j=1,np
       DO k=1,NPATCH2D
          nc = 0
          ivec = 0
          pop%pop_grid(j)%patch(k)%stress_mortality = 0.0
          pop%pop_grid(j)%patch(k)%fire_mortality = 0.0
          pop%pop_grid(j)%patch(k)%crowding_mortality = 0.0
          pop%pop_grid(j)%patch(k)%mortality = 0.0
          crown_area = 0.0
          DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
             cmass_stem_inc=StemNPP(j,1)*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake

             if (present(StemNPP_av)) then
                growth_efficiency=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake* &
                     StemNPP_av(j)  /(cmass_stem**(POWERGrowthEfficiency))
             else
                growth_efficiency=cmass_stem_inc/(cmass_stem**(POWERGrowthEfficiency))
             endif
             ! growth_efficiency=cmass_stem_inc/(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood**(POWERGrowthEfficiency))
             mort=MORT_MAX/(1.0+(growth_efficiency/GROWTH_EFFICIENCY_MIN)**Pmort)

             ! mort = 0 ! test

             pop%pop_grid(j)%patch(k)%stress_mortality = pop%pop_grid(j)%patch(k)%stress_mortality &
                  + mort*cmass_stem
             IF (pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter*100.GT.1.) THEN
                if (ALLOM_SWITCH.eq.1) then
                   ! assumes crown radius (m) = 0.14 * dbh (cm)
                   crown_area = crown_area + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density* &
                        PI*(pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter*100.*0.14)**2
                else
                   crown_area = crown_area + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density* &
                        k_allom1 * pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter ** k_rp
                endif
             ELSE
                crown_area = crown_area + 0.5*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%LAI
             ENDIF

             cpc = 1. - exp(-crown_area)
             pop%pop_grid(j)%patch(k)%cpc = cpc
             if (cpc.gt.1e-3 .and.  alpha_cpc * (1. - 1./cpc).gt.-50.0) then
                mort_cpc = exp(alpha_cpc * (1. - 1./cpc))
             else
                mort_cpc = 0.
             endif

             !mort_cpc = 0 ! test
             pop%pop_grid(j)%patch(k)%crowding_mortality = pop%pop_grid(j)%patch(k)%crowding_mortality + &
                  min((mort_cpc*CrowdingFactor),cmass_stem_inc/cmass_stem)*cmass_stem

             mort = mort + min((mort_cpc*CrowdingFactor),cmass_stem_inc/cmass_stem)

             pop%pop_grid(j)%patch(k)%mortality = pop%pop_grid(j)%patch(k)%mortality + mort*cmass_stem

             pop%pop_grid(j)%patch(k)%sapwood_loss =  pop%pop_grid(j)%patch(k)%sapwood_loss + &
                  mort*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood

             pop%pop_grid(j)%patch(k)%sapwood_area_loss =  pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                  mort*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood_area


             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = cmass_stem*(1.-mort)
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood *(1.-mort)
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood *(1.-mort)

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density*(1.-mort)
             IF (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.LT.DENSINDIV_MIN) THEN
                ! remove cohort
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 0
                pop%pop_grid(j)%patch(k)%stress_mortality = pop%pop_grid(j)%patch(k)%stress_mortality + &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood_area
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood = 0.0

             ELSE
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 1
                !COMMLN Why is id here 1 instead of c or sth useful? Call it differently
                nc = nc+1
                ivec(nc)=c
             ENDIF
          ENDDO
          ! SHUFFLE if necessary to remove zero-density cohorts
          IF (nc.LT.pop%pop_grid(j)%patch(k)%Layer(1)%ncohort) THEN
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ivec(1:nc))
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = nc

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%density = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%id      = 0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%biomass = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood_area = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%basal_area = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%heartwood = 0.0
          ENDIF
          ! Layer biomass (summed over cohorts)
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          pop%pop_grid(j)%patch(k)%Layer(1)%biomass = SUM(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%biomass)
       ENDDO
    ENDDO

    ! recruitment
    IF (PRESENT(precip)) THEN
       CALL layer_recruitment(pop, precip)
    ELSE
       CALL layer_recruitment(pop)
    ENDIF

    ! Update time since last patch disturbance
    DO j=1,np
       DO k=1,NPATCH2D

          DO idist =1, NDISTURB
             pop%pop_grid(j)%patch(k)%age(idist) = pop%pop_grid(j)%patch(k)%age(idist) + 1
          ENDDO

       ENDDO

    ENDDO


  END SUBROUTINE PatchAnnualDynamics
  !*******************************************************************************
  SUBROUTINE GetUniqueAgeFrequencies(pop, disturbance_interval, idisturb, it)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:), idisturb, it(:)
    INTEGER(i4b) :: g, i,j,k,ct,lastct,agecopy,idcopy
    REAL(dp), ALLOCATABLE :: midpoint(:)
    INTEGER(i4b), ALLOCATABLE :: ranked_age(:), ranked_age_init(:)
    INTEGER(i4b) ::  tmp_count, tmp_i, age_tmp
    INTEGER(i4b), ALLOCATABLE :: ranked_age_unique_id(:), ranked_age_id(:), counter(:)
    REAL(dp), ALLOCATABLE :: tmp(:), freq_tmp(:), freq_tmp1(:)
    REAL(dp) :: p,cump,lastcump, freq, tmp1
    INTEGER(i4b) :: n_age ! number of unique ages
    INTEGER(i4b) :: npatch_active ! number of active patches
    REAL(dp):: disturbance_freq
    INTEGER(i4b) :: i_max, age_max, Poisson_age(1000), np
    REAL(dp):: Poisson_weight(1000), CumPoisson_weight(1000)
    INTEGER(i4b), ALLOCATABLE :: bound(:,:), unique_age(:)

    !Fills array freq with weights (frequencies across landscape) for each unique age
    ! given specified mean disturbance interval

    np = SIZE(POP%POP_grid)
    DO g=1,np

       npatch_active = NPATCH2D
       IF (.NOT.ALLOCATED(midpoint)) ALLOCATE(midpoint(npatch_active))
       IF (.NOT.ALLOCATED(counter)) ALLOCATE(counter(npatch_active))
       IF (.NOT.ALLOCATED(ranked_age)) ALLOCATE(ranked_age(npatch_active))
       IF (.NOT.ALLOCATED(ranked_age_init)) ALLOCATE(ranked_age_init(npatch_active))
       IF (.NOT.ALLOCATED(ranked_age_id)) ALLOCATE(ranked_age_id(npatch_active))
       IF (.NOT.ALLOCATED(ranked_age_unique_id)) ALLOCATE(ranked_age_unique_id(npatch_active))
       IF (.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(npatch_active))
       IF (.NOT.ALLOCATED(freq_tmp)) ALLOCATE(freq_tmp(npatch_active))
       IF (.NOT.ALLOCATED(freq_tmp1)) ALLOCATE(freq_tmp1(npatch_active))


       ! rank patches in order of age
       pop%pop_grid(g)%ranked_age_unique(:, idisturb) = 0
       ranked_age_init = pop%pop_grid(g)%patch%age(idisturb)
       ranked_age = pop%pop_grid(g)%patch%age(idisturb)
       ranked_age_id = pop%pop_grid(g)%patch%id
       ranked_age_unique_id = 0
       freq_tmp = 0.0
       freq = 0
       pop%pop_grid(g)%freq_ranked_age_unique(:, idisturb) = 0.0
       midpoint = 0.0


       DO i = 1, npatch_active -1
          DO j = i+1, npatch_active
             IF (ranked_age(i).GT.ranked_age(j)) THEN
                agecopy          = ranked_age(i)
                idcopy           = ranked_age_id(i)
                ranked_age(i)    = ranked_age(j)
                ranked_age_id(i) = ranked_age_id(j)
                ranked_age(j)    = agecopy
                ranked_age_id(j) = idcopy
             ENDIF
          ENDDO
       ENDDO

       ! subset to unique ages
       k=0
       age_tmp = -1
       DO i = 1, npatch_active
          IF (ranked_age(i).NE.age_tmp) k = k+1
          pop%pop_grid(g)%ranked_age_unique(k, idisturb) = ranked_age(i)
          ranked_age_unique_id(k) = ranked_age_id(i)
          age_tmp = ranked_age(i)
          n_age  = k
       ENDDO

       disturbance_freq=1.0/REAL(disturbance_interval(g,idisturb))
       DO i =1,1000
          Poisson_age(i) = i
          CumPoisson_weight(i) = CumExponential(disturbance_freq,REAL(i,dp))

       ENDDO


       ! construct upper and lower bounds for each unique age: these set the range of ages to be
       ! represented by an unique age
       ALLOCATE(bound(n_age,2))
       ALLOCATE (unique_age(n_age))
       bound = 0
       unique_age = pop%pop_grid(g)%ranked_age_unique(1:n_age,idisturb)
       DO i=1,n_age
          IF (unique_age(i).EQ.0) THEN
             bound(i,1) = 0
             bound(i,2) = 0
          ELSEIF ((i.EQ.1).AND.(unique_age(i).GT.0)) THEN
             bound(i,1) = 0
             bound(i,2) = unique_age(i)
          ELSEIF ((unique_age(i).GT.0).AND.(i.GT.1).AND.(unique_age(i-1).EQ.unique_age(i)-1)) THEN
             bound(i,1) = unique_age(i)
             IF (i.LT.n_age) THEN
                bound(i,2) = unique_age(i)
             ELSE
                i_max = MAXLOC(Poisson_age, 1, CumPoisson_weight.LE.0.99)
                bound(i, 2) = Poisson_age(i_max)
             ENDIF
          ELSEIF ((unique_age(i).GT.0).AND.(i.GT.1).AND.(unique_age(i-1).NE.unique_age(i)-1)) THEN
             bound(i,1) = bound(i-1,2)+1
             IF (i.LT.n_age) THEN
                bound(i,2) = (unique_age(i)+ unique_age(i+1))/2
             ELSE
                i_max = MAXLOC(Poisson_age, 1, CumPoisson_weight.LE.0.99)
                bound(i, 2) = Poisson_age(i_max)
             ENDIF
          ENDIF
          

       ENDDO

       ! calculate weighting for each unique age
       DO i=1,n_age
          DO j = bound(i,1),bound(i,2)

             !IF (pop%LU(g)==2) THEN  ! secondary forest
             IF (POP%pop_grid(g)%LU ==2) THEN
                freq_tmp(i) = freq_tmp(i) +  pop%pop_grid(g)%freq_age(j+1)
             ELSE
                freq_tmp(i) = freq_tmp(i) + REALExponential(disturbance_freq,REAL(j,dp))
             ENDIF

          ENDDO
       ENDDO

       pop%pop_grid(g)%freq_ranked_age_unique(1:npatch_active,idisturb) = freq_tmp
       pop%pop_grid(g)%n_age(idisturb) = n_age

       DEALLOCATE (bound)
       DEALLOCATE (unique_age)

    ENDDO

  END SUBROUTINE GetUniqueAgeFrequencies
  !*******************************************************************************
  SUBROUTINE GetPatchFrequencies(pop,it)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b), INTENT(IN) :: it(:)
    INTEGER(i4b) :: n1, n2, g, REPCOUNT, tmp1(NPATCH1D), np, idist
    REAL(dp) ::  tmp2(NPATCH1D), tmp3(NPATCH1D), sum_freq

    np = SIZE(Pop%pop_grid)
    DO g=1,np

       pop%pop_grid(g)%freq = 0.0
       DO idist = 1, NDISTURB
          IF (idist.EQ.1) THEN
             DO n1=1,pop%pop_grid(g)%n_age(1)
                repcount = COUNT(pop%pop_grid(g)%patch(:)%age(1).EQ.pop%pop_grid(g)%ranked_age_unique(n1,1))
                WHERE (pop%pop_grid(g)%patch(:)%age(1).EQ.pop%pop_grid(g)%ranked_age_unique(n1,1))
                   pop%pop_grid(g)%freq = pop%pop_grid(g)%freq_ranked_age_unique(n1,1) /REAL(repcount,dp)
                ENDWHERE
             ENDDO


          ELSEIF (idist.EQ.2) THEN
             ! first calculate weights for patches with age(2)>age(1)
             DO n1=1,pop%pop_grid(g)%n_age(1)


                DO n2=1,pop%pop_grid(g)%n_age(idist)
                   repcount = COUNT((pop%pop_grid(g)%patch(1:NPATCH)%age(1) .EQ. &
                        pop%pop_grid(g)%ranked_age_unique(n1,1)).AND. &
                        (pop%pop_grid(g)%patch(1:NPATCH)%age(idist) .EQ.  &
                        pop%pop_grid(g)%ranked_age_unique(n2,idist)))
                   WHERE ((pop%pop_grid(g)%patch(1:NPATCH)%age(1).EQ.pop%pop_grid(g)%ranked_age_unique(n1,1)).AND. &
                        (pop%pop_grid(g)%patch(1:NPATCH)%age(idist).EQ.pop%pop_grid(g)%ranked_age_unique(n2,idist)))
                      pop%pop_grid(g)%freq(1:NPATCH) = pop%pop_grid(g)%freq_ranked_age_unique(n1,1)* &
                           pop%pop_grid(g)%freq_ranked_age_unique(n2,idist) &
                           /REAL(repcount,dp)

                   ENDWHERE


                ENDDO
             ENDDO


          ENDIF
       ENDDO ! end loop over idist

       sum_freq = SUM(pop%pop_grid(g)%freq)
       if (sum_freq.gt.0.0) pop%pop_grid(g)%freq = pop%pop_grid(g)%freq/sum_freq

    ENDDO


  END SUBROUTINE GetPatchFrequencies

  !*******************************************************************************
  SUBROUTINE GetDiagnostics(pop,LAI,Cleaf,Croot,disturbance_interval, it, precip)
    ! Gets diagnostic data for current landscape structure
    IMPLICIT NONE
    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    REAL(dp), INTENT(IN) ::  LAI(:)
    REAL(dp), INTENT(IN) ::  Cleaf(:)
    REAL(dp), INTENT(IN) ::  Croot(:)
    INTEGER(i4b), INTENT(IN)        ::  disturbance_interval(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    INTEGER(i4b), INTENT(IN) :: it(:)
    INTEGER(i4b) :: P, g,i,j,ct, ct_highres
    REAL(dp) :: limits(HEIGHT_BINS+1)
    REAL(dp) :: ht, htmax, cmass_stem,densindiv, freq, freq_old
    CHARACTER(len=12) :: string1, string2
    CHARACTER(len=9) :: fmt
    INTEGER(i4b) :: npatch_active  ! number of active patches
    INTEGER(i4b) :: np, nc, i_height, id
    REAL(dp) :: diam,basal, cump
    REAL(dp) :: patch_crown_area(NPATCH2D), patch_crown_cover(NPATCH2D)
    REAL(dp), ALLOCATABLE :: height_list(:), height_list_weight(:)
    REAL(dp) :: height_copy, weight_copy, Pwc, FAVD
    INTEGER(i4b), PARAMETER :: HEIGHT_BINS_highres=100 ! bins for assessing height_max
    REAL(dp), ALLOCATABLE :: limits_highres(:), DENSINDIV_HIGHRES(:)
    REAL(dp) :: a, b, C, res_flux, fire_flux, cat_flux, vol, vol1, tmp, tmp2
    integer :: arg1

    fmt = '(f5.1)'
    limits(1) = 0.
    IF(.NOT.ALLOCATED(limits_highres)) ALLOCATE(limits_highres(HEIGHT_BINS_highres+1))
    IF(.NOT.ALLOCATED(DENSINDIV_HIGHRES)) ALLOCATE(DENSINDIV_HIGHRES(HEIGHT_BINS_highres))


    limits_highres(1) = 0.
    np = SIZE(Pop%pop_grid)

    DO g=1,np
       npatch_active = NPATCH2D
       IF (MAX_HEIGHT_SWITCH.EQ.1) THEN
          ALLOCATE(height_list(NPATCH2D*NCOHORT_MAX))
          ALLOCATE(height_list_weight(NPATCH2D*NCOHORT_MAX))
       ENDIF
     !  IF(.NOT.ALLOCATED(MASK)) ALLOCATE(MASK(POP%pop_grid%npatch_active))



       DO i=1,HEIGHT_BINS
          limits(i+1) = BIN_POWER**REAL(i)
          WRITE(string1,fmt) (limits(i))
          WRITE(string2,fmt) (limits(i+1))
          pop%pop_grid(g)%bin_labels(i) = 'Height_'//TRIM(ADJUSTL(string1))//'-'//TRIM(ADJUSTL(string2))//'m'
          pop%pop_grid(g)%cmass_stem_bin(i) = 0.0
          pop%pop_grid(g)%densindiv_bin(i) = 0.0
          pop%pop_grid(g)%cmass_stem_bin(i) = 0.0
          pop%pop_grid(g)%height_bin(i) = REAL(limits(i)+limits(i+1))/2.
          pop%pop_grid(g)%diameter_bin(i) = ((REAL(limits(i))/Kbiometric)**(3/2)+(REAL(limits(i+1))/Kbiometric)**(3/2))/2.
       ENDDO

       DO i=1,HEIGHT_BINS_highres
          limits_highres(i+1) = REAL(i)
       ENDDO

       IF (MAX_HEIGHT_SWITCH.EQ.1) THEN
          height_list = 0.0
          height_list_weight = 0.0
       ENDIF
       i_height = 0
       pop%pop_grid(g)%cmass_sum_old = pop%pop_grid(g)%cmass_sum
       pop%pop_grid(g)%csapwood_sum_old = pop%pop_grid(g)%csapwood_sum
       pop%pop_grid(g)%cmass_sum = 0.0
       pop%pop_grid(g)%csapwood_sum = 0.0
       pop%pop_grid(g)%cheartwood_sum = 0.0
       pop%pop_grid(g)%height_mean = 0.0
       pop%pop_grid(g)%fire_mortality = 0.0
       pop%pop_grid(g)%cat_mortality = 0.0
       pop%pop_grid(g)%res_mortality = 0.0
       pop%pop_grid(g)%stress_mortality = 0.0
       pop%pop_grid(g)%crowding_mortality = 0.0
       pop%pop_grid(g)%sapwood_loss = 0.0
       pop%pop_grid(g)%sapwood_area_loss = 0.0
       pop%pop_grid(g)%growth = 0.0
       pop%pop_grid(g)%area_growth = 0.0
       pop%pop_grid(g)%basal_area = 0.0
       pop%pop_grid(g)%densindiv = 0.0
       pop%pop_grid(g)%height_max = 0.0
       pop%pop_grid(g)%crown_cover = 0.0
       pop%pop_grid(g)%crown_area = 0.0
       pop%pop_grid(g)%sapwood_area = 0.0
       pop%pop_grid(g)%crown_volume = 0.0
       densindiv_highres = 0.0
       ! loop through patches
       DO P = 1, npatch_active
          pop%pop_grid(g)%patch(p)%biomass = 0.0
          pop%pop_grid(g)%patch(p)%sapwood = 0.0
          pop%pop_grid(g)%patch(p)%sapwood_area = 0.0
          pop%pop_grid(g)%patch(p)%heartwood = 0.0
          pop%pop_grid(g)%patch(p)%layer(1)%biomass = 0.0
          pop%pop_grid(g)%patch(p)%layer(1)%density = 0.0
          patch_crown_area(p) = 0.0
          patch_crown_cover(p) = 0.0
          tmp2 = sum(pop%pop_grid(g)%patch(p)%layer(1)%cohort(1:pop%pop_grid(g)%patch(p)%layer(1)%ncohort)%sapwood_area)

          freq =  pop%pop_grid(g)%freq(pop%pop_grid(g)%patch(p)%id)
          freq_old = pop%pop_grid(g)%freq_old(pop%pop_grid(g)%patch(p)%id)

          ! loop through cohorts
          DO i = 1, pop%pop_grid(g)%patch(p)%layer(1)%ncohort
             cmass_stem = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%biomass
             densindiv = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%density

             IF ( PRESENT(precip) ) THEN
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal, precip(g))
             ELSE
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal )

             ENDIF

             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%height   = ht
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter = diam

             ! basal area in each cohort
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%basal_area = basal

             ! sapwood area in each cohort
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area = basal - & ! m2 ha-1
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%heartwood/(ht*WD)*1.e4

             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area = &
                  max(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area, 0.0)

             ! get bin
             ct = 1
             DO j=1,HEIGHT_BINS
                IF (ht.GT.limits(j)) ct = j
             ENDDO ! bins

             ! get high res bin
             ct_highres = 1
             DO j=1,HEIGHT_BINS_highres
                IF (ht.GT.limits_highres(j)) ct_highres = j
             ENDDO ! bins



             pop%pop_grid(g)%patch(p)%layer(1)%biomass = pop%pop_grid(g)%patch(p)%layer(1)%biomass + cmass_stem
             pop%pop_grid(g)%patch(p)%layer(1)%density = pop%pop_grid(g)%patch(p)%layer(1)%density + densindiv


             IF (diam*100.0.GT.1.) THEN
                patch_crown_area(p) = patch_crown_area(p) + densindiv*PI*(diam*100.*0.1492)**2 ! uses GC relationship
                pop%pop_grid(g)%crown_volume = pop%pop_grid(g)%crown_volume + &
                     freq*densindiv*(4./3.)*PI*(diam*100.*0.1492)**2*(1.5*(diam*100.*0.1492))
                ! assumes vertical radius = 1.5 * horizontal radius
             ENDIF

             IF (diam*100..GT.5.) THEN
                if (ALLOM_SWITCH.eq.1) then
                   !! assumes crown radius (m) = 0.1492 * dbh (cm) (from G. Cook, pers. comm.)
                   ! assumes vertical radius = 1.5 * horizontal radius
                   pop%pop_grid(g)%crown_volume = pop%pop_grid(g)%crown_volume + &
                        freq*densindiv*(4./3.)*PI*(diam*100.*0.1492)**2*(1.5*(diam*100.*0.1492))
                else
                   !! global allometry
                   ! assumes vertical radius = 1.5 * horizontal radius
                   pop%pop_grid(g)%crown_volume = pop%pop_grid(g)%crown_volume + &
                        freq*densindiv*(4./3.)*PI*1.5*((k_allom1 * diam ** k_rp )/PI)**1.5
                endif

             ENDIF

             pop%pop_grid(g)%patch(p)%sapwood = pop%pop_grid(g)%patch(p)%sapwood + &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood
             pop%pop_grid(g)%patch(p)%sapwood_area = pop%pop_grid(g)%patch(p)%sapwood_area + &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area
             pop%pop_grid(g)%patch(p)%heartwood = pop%pop_grid(g)%patch(p)%heartwood + &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%heartwood
             pop%pop_grid(g)%patch(p)%biomass = pop%pop_grid(g)%patch(p)%biomass + cmass_stem
             pop%pop_grid(g)%cmass_stem_bin(ct) = pop%pop_grid(g)%cmass_stem_bin(ct) + freq*cmass_stem
             pop%pop_grid(g)%densindiv_bin(ct) = pop%pop_grid(g)%densindiv_bin(ct) + freq*densindiv
             pop%pop_grid(g)%cmass_sum = pop%pop_grid(g)%cmass_sum + freq*cmass_stem
             pop%pop_grid(g)%csapwood_sum = pop%pop_grid(g)%csapwood_sum + &
                  freq*pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood
             pop%pop_grid(g)%cheartwood_sum = pop%pop_grid(g)%cheartwood_sum + &
                  freq*pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%heartwood

             pop%pop_grid(g)%densindiv = pop%pop_grid(g)%densindiv + freq*densindiv
             pop%pop_grid(g)%height_mean = pop%pop_grid(g)%height_mean + ht*freq*densindiv
             pop%pop_grid(g)%basal_area = pop%pop_grid(g)%basal_area +basal*freq
             densindiv_highres(ct_highres) = densindiv_highres(ct_highres) + freq*densindiv
             IF (MAX_HEIGHT_SWITCH.EQ.1) THEN
                i_height = i_height+1
                height_list(i_height) = ht
                height_list_weight(i_height) = densindiv*freq
             ENDIF
          ENDDO ! cohorts

          pop%pop_grid(g)%stress_mortality = pop%pop_grid(g)%stress_mortality + &
               freq*pop%pop_grid(g)%patch(p)%stress_mortality

          pop%pop_grid(g)%crowding_mortality = pop%pop_grid(g)%crowding_mortality + &
               freq*pop%pop_grid(g)%patch(p)%crowding_mortality

          pop%pop_grid(g)%sapwood_area = pop%pop_grid(g)%sapwood_area + &
               freq*sum(pop%pop_grid(g)%patch(p)%layer(1)%cohort(1:pop%pop_grid(g)%patch(p)%layer(1)%ncohort)%sapwood_area)

          pop%pop_grid(g)%sapwood_loss = pop%pop_grid(g)%sapwood_loss + &
               freq*pop%pop_grid(g)%patch(p)%sapwood_loss

          pop%pop_grid(g)%sapwood_loss =  pop%pop_grid(g)%sapwood_loss + &
               pop%pop_grid(g)%patch(p)%sapwood_old*(freq_old-freq)

          pop%pop_grid(g)%sapwood_area_loss = pop%pop_grid(g)%sapwood_area_loss + &
               freq*pop%pop_grid(g)%patch(p)%sapwood_area_loss

          pop%pop_grid(g)%sapwood_area_loss =  pop%pop_grid(g)%sapwood_area_loss + &
               pop%pop_grid(g)%patch(p)%sapwood_area_old*(freq_old-freq)


          pop%pop_grid(g)%cat_mortality  = pop%pop_grid(g)%cat_mortality +  &
               freq*pop%pop_grid(g)%patch(p)%cat_mortality
          pop%pop_grid(g)%res_mortality = pop%pop_grid(g)%res_mortality + &
               pop%pop_grid(g)%patch(p)%biomass_old*(freq_old-freq)
          pop%pop_grid(g)%fire_mortality = pop%pop_grid(g)%fire_mortality + &
               freq*pop%pop_grid(g)%patch(p)%fire_mortality
          pop%pop_grid(g)%growth =  pop%pop_grid(g)%growth + freq_old*pop%pop_grid(g)%patch(p)%growth
          !write(*,*) 'freq_old',freq_old


          pop%pop_grid(g)%area_growth =  pop%pop_grid(g)%area_growth + &
               freq*pop%pop_grid(g)%patch(p)%area_growth

       ENDDO ! patches



IF (INTERP_SWITCH==1.and.NDISTURB.eq.2) then
   !CALL INTERPOLATE_BIOMASS_2D(pop, disturbance_interval,it)
   CALL INTERPOLATE_BIOMASS_2D(pop, disturbance_interval,it(g),g)

ELSEIF (INTERP_SWITCH==1.and.NDISTURB.eq.1) then

   CALL INTERPOLATE_BIOMASS_1D(pop, disturbance_interval,it(g),g)

ENDIF

arg1 = NYEAR_SMOOTH-(NYEAR_SMOOTH/2)
IF (SMOOTH_SWITCH==1) THEN
   IF (it(g).LE.NYEAR_SMOOTH-NYEAR_SMOOTH/2) THEN
      CALL SMOOTH_FLUX(POP,g,it(g))
   ELSE
      CALL SMOOTH_FLUX(POP,g,int(arg1,i4b))
   ENDIF
ENDIF


       ! leaf area index in each cohort
       DO P = 1, npatch_active

          freq =  pop%pop_grid(g)%freq(pop%pop_grid(g)%patch(p)%id)
          ! loop through cohorts
          DO i = 1, pop%pop_grid(g)%patch(p)%layer(1)%ncohort
             cmass_stem = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%biomass
             densindiv = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%density
             basal=PI*(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter/2.0)* &
                  (pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter/2.0)*densindiv*1.e4
             ! leaf area index in each cohort
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%LAI = LAI(g) * &
                  min(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area  &
                  /max(pop%pop_grid(g)%sapwood_area,1e-3), 10.0_dp)
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Cleaf = Cleaf(g) * &
                  min(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area  &
                  /max(pop%pop_grid(g)%sapwood_area,1.e-3), 10.0_dp)
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Croot = Croot(g) * &
                  min(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area  &
                  /max(pop%pop_grid(g)%sapwood_area,1e-3), 10.0_dp)
          ENDDO ! cohorts
          pop%pop_grid(g)%patch(p)%LAI = sum(pop%pop_grid(g)%patch(p)%layer(1)% &
               cohort(1:pop%pop_grid(g)%patch(p)%layer(1)%ncohort)%LAI)
       ENDDO ! patches

       ! PGap = (1-fcover) calculation

       if (pop%pop_grid(g)%crown_volume>0.0) then
          FAVD = LAI(g)/pop%pop_grid(g)%crown_volume ! foliage area volume density
       else
          FAVD = 0.0
       endif

       DO P = 1, npatch_active
          freq =  pop%pop_grid(g)%freq(pop%pop_grid(g)%patch(p)%id)
          nc =  pop%pop_grid(g)%patch(p)%layer(1)%ncohort
          ! loop through cohorts
          DO i = 1, nc
             cmass_stem = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%biomass
             densindiv = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%density

             IF ( PRESENT(precip) ) THEN
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal, precip(g))
             ELSE
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal )
             ENDIF


             IF (diam*100.GT.1.) THEN

                if (ALLOM_SWITCH.eq.1) then
                   !! assumes crown radius (m) = 0.1492 * dbh (cm) (from G. Cook, pers. comm.)
                   pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = densindiv*PI*(diam*100.*0.1492)**2
                   Pwc = EXP(-0.5 * pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%LAI/ &
                        pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area)
                   pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = densindiv*PI*(diam*100.*0.1492)**2*(1.-Pwc)

                else
                   !! global allometry
                   pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = &
                        densindiv*PI*(((k_allom1 * diam ** k_rp )/PI)**0.5)**2
                   Pwc = EXP(max(-0.5 * pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%LAI/ &
                        max(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area,1.e-3),-20.0))

                   pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = &
                        pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area*(1.-Pwc) !*1.4142
                endif

             ELSE
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = &
                     0.5*pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%LAI !*1.4142
             ENDIF
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area= &
                  max(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area,0.01)

             pop%pop_grid(g)%crown_area = pop%pop_grid(g)%crown_area + &
                  freq*pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area

             IF (i.eq.1) THEN ! top cohort
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Pgap = &
                     exp(-pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area)

                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%frac_interception = &
                     1- exp(-pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area)
             ELSE
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Pgap = &
                     pop%pop_grid(g)%patch(p)%layer(1)%cohort(i-1)%Pgap* &
                     exp(-pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area)
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%frac_interception = &
                     pop%pop_grid(g)%patch(p)%layer(1)%cohort(i-1)%Pgap - &
                     pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Pgap
             ENDIF
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%respiration_scalar =  &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood/shootfrac/CtoNw + &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Cleaf/CtoNl + &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Croot/CtoNr

          ENDDO ! cohorts

          IF (nc>0) THEN
             pop%pop_grid(g)%patch(p)%pgap  = &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(nc)%Pgap
          ELSE
             pop%pop_grid(g)%patch(p)%pgap  = 1
          ENDIF


       ENDDO ! patches

       pop%pop_grid(g)%Kclump = max(pop%pop_grid(g)%crown_area/(0.5*LAI(g)),0.1_dp)
       pop%pop_grid(g)%crown_cover = 1.-EXP(-pop%pop_grid(g)%crown_area)

      
       pop%pop_grid(g)%height_mean = pop%pop_grid(g)%height_mean/max(pop%pop_grid(g)%densindiv,1.0e-5)

       ! Height Diagnostics
       IF (MAX_HEIGHT_SWITCH.EQ.0) THEN
          ! Set landscape maximum height to centre of bin with <5% of trees in a bin of higher size classes
          cump = 0.
          j = 1
          DO WHILE (cump.LT.0.95)
             cump = cump + pop%pop_grid(g)%densindiv_bin(j)/max(pop%pop_grid(g)%densindiv,1.0e-5)
             pop%pop_grid(g)%height_max = pop%pop_grid(g)%height_bin(j)
             j = j+1
          ENDDO
       ELSEIF (MAX_HEIGHT_SWITCH.EQ.1) THEN

          ! sort height list
          DO i = 1, i_height -1
             DO j = i+1, i_height
                IF (height_list(i).GT.height_list(j)) THEN
                   height_copy = height_list(i)
                   weight_copy =  height_list_weight(i)
                   height_list(i) = height_list(j)
                   height_list_weight(i) = height_list_weight(j)
                   height_list(j) = height_copy
                   height_list_weight(j) = weight_copy
                ENDIF
             ENDDO
          ENDDO ! end sort height list

          ! normailse height list weights
          height_list_weight=height_list_weight/SUM(height_list_weight(1:i_height))
          cump = 0.
          j = 1
          DO WHILE (cump.LT.0.95)
             cump = cump + height_list_weight(j)
             pop%pop_grid(g)%height_max = height_list(j)
             j = j+1
          ENDDO
          DEALLOCATE(height_list)
          DEALLOCATE(height_list_weight)

       ELSEIF (MAX_HEIGHT_SWITCH.EQ.2) THEN
          cump = 0.
          j = 1
          densindiv_highres= densindiv_highres/max(SUM(densindiv_highres),1.0e-5)
          DO WHILE ((cump.LT.0.95).AND.(j.LE.HEIGHT_BINS_highres))
             cump = cump + densindiv_highres(j)
             pop%pop_grid(g)%height_max = (limits_highres(j+1) + limits_highres(j))/2.
             j = j+1
          ENDDO
       ENDIF
    !deallocate(MASK)

    ENDDO ! end loop over grid cells


  END SUBROUTINE GetDiagnostics
  !*******************************************************************************
  SUBROUTINE Patch_partial_disturb(pop,idisturb,intensity,precip,frac_intensity1)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b), INTENT(IN) ::  idisturb
    REAL(dp), INTENT(IN) :: intensity(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:), frac_intensity1(:)
    INTEGER(i4b) :: j, k, i, g, c, nc, np
    INTEGER(i4b) ::  ivec(NCOHORT_MAX)
    REAL(dp) :: ht, diam
    REAL(dp) :: Psurvival_l, Psurvival_s, Psurvival, char_height

    np = SIZE(Pop%pop_grid)
    !Print*,"CLN Hier PROBLEM PSurv1 unten l1345 wegen PRESENT() "

    ! Kills a fraction of biomass in patch when prescribed disturbance interval is reached
    DO j=1,np
       DO k=1,NPATCH
          pop%pop_grid(j)%patch(k)%fire_mortality = 0.0

          IF (((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).NE.0).AND. &
               (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))).OR. &
               (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))) THEN


             ! loop through cohorts
             ivec = 0
             nc = 0
             DO c = 1, pop%pop_grid(j)%patch(k)%layer(1)%ncohort
                ! kill fraction of each cohort
                char_height = 3.7*(1.-EXP(-0.19*Intensity(j,1)))
                ht = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%height
                diam = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter*100. ! diameter in cm
                IF ((ht.GT.8.5).AND.(ht.GT.char_height)) THEN
                   Psurvival_s  =(-0.0011*Intensity(j,1) -0.00002)*ht &
                        +(0.0075*Intensity(j,1)+1.)

                ELSEIF ((ht.LE.8.5).AND.(ht.GT.char_height)) THEN
                   Psurvival_s =(0.0178*Intensity(j,1) + 0.0144)*ht &
                        + (-0.1174*Intensity(j,1)+0.9158)

                ELSE
                   Psurvival_s = 0.0
                ENDIF
                Psurvival_s = MIN(Psurvival_s,1.0_dp)
                Psurvival_s = MAX(Psurvival_s,1e-3_dp)
                Psurvival = Psurvival_s


                IF (PRESENT(frac_intensity1)) THEN
                   char_height = 3.7*(1.-EXP(-0.19*Intensity(j,2)))
                   ht = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%height
                   IF ((ht.GT.8.5).AND.(ht.GT.char_height)) THEN
                      Psurvival_s  =(-0.0011*Intensity(j,2) -0.00002)*ht &
                           +(0.0075*Intensity(j,2)+1.)
                   ELSEIF ((ht.LE.8.5).AND.(ht.GT.char_height)) THEN
                      Psurvival_s =(0.0178*Intensity(j,2) + 0.0144)*ht &
                           + (-0.1174*Intensity(j,2)+0.9158)
                   ELSE
                      Psurvival_s = 0.0
                   ENDIF
                   Psurvival_s = MIN(Psurvival_s,1.0_dp)
                   Psurvival_s = MAX(Psurvival_s,1e-3_dp)
                   Psurvival = Psurvival_s*(1.-frac_intensity1(j)) + Psurvival*frac_intensity1(j)
                ENDIF
                !Psurvival = 1 ! test
                pop%pop_grid(j)%patch(k)%fire_mortality = pop%pop_grid(j)%patch(k)%fire_mortality + &
                     (1.-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                     (1.-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                     (1.-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood_area
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%heartwood = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%heartwood
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density
                IF (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.LT.DENSINDIV_MIN) THEN
                   ! remove cohort
                   pop%pop_grid(j)%patch(k)%fire_mortality = pop%pop_grid(j)%patch(k)%fire_mortality + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                   pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                   pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood_area
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood_area = 0.0

                ELSE
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 1
                   nc = nc+1
                   ivec(nc)=c
                ENDIF
             ENDDO

             ! SHUFFLE if necessary to remove zero-density cohorts
             IF (nc.LT.pop%pop_grid(j)%patch(k)%Layer(1)%ncohort) THEN
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ivec(1:nc))
                pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = nc

                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%density = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%id = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%biomass = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood_area = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%heartwood = 0.0
             ENDIF

             pop%pop_grid(j)%patch(k)%age(idisturb) = 0
             pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0


          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE Patch_partial_disturb
  !*******************************************************************************
  SUBROUTINE Patch_partial_disturb2(pop,idisturb,precip)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b), INTENT(IN) ::  idisturb
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    INTEGER(i4b) :: j, k, i, g, c, nc, np
    INTEGER(i4b) ::  ivec(NCOHORT_MAX)
    REAL(dp) :: ht, diam
    REAL(dp) :: Psurvival_l, Psurvival_s, Psurvival, char_height, frac_mort, Pmort

    np = SIZE(Pop%pop_grid)

    ! Kills a fraction (80%) biomass in patch when prescribed disturbance interval is reached
    DO j=1,np
       DO k=1,NPATCH2D
          pop%pop_grid(j)%patch(k)%cat_mortality = 0.
          ! Layer biomass (summed over cohorts)
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          pop%pop_grid(j)%patch(k)%Layer(1)%biomass = SUM(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%biomass)

          IF (((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).NE.0).AND. &
               (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))).OR. &
               (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))) THEN


             ! loop through cohorts
             ivec = 0
             nc = 0
             frac_mort = 0.0
             pop%pop_grid(j)%patch(k)%cat_mortality = 0.0
             DO c = 1, pop%pop_grid(j)%patch(k)%layer(1)%ncohort
                ! kill fraction of each cohort, up to 80% of patch biomass

                if (pop%pop_grid(j)%patch(k)%cat_mortality  < 0.8 * pop%pop_grid(j)%patch(k)%Layer(1)%biomass ) then

                   Pmort = min(  (0.8*pop%pop_grid(j)%patch(k)%Layer(1)%biomass-pop%pop_grid(j)%patch(k)%fire_mortality) &
                        /pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass, 1.0_dp)

                else

                   Pmort = 0.0
                endif
                Psurvival = 1- Pmort



                pop%pop_grid(j)%patch(k)%cat_mortality = pop%pop_grid(j)%patch(k)%cat_mortality + &
                     (1.-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                     (1.-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                     (1.-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood_area
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%heartwood = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%heartwood
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density
                IF (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.LT.DENSINDIV_MIN) THEN
                   ! remove cohort
                   pop%pop_grid(j)%patch(k)%cat_mortality = pop%pop_grid(j)%patch(k)%cat_mortality + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                   pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                   pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood_area
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood_area = 0.0

                ELSE
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 1
                   nc = nc+1
                   ivec(nc)=c
                ENDIF
             ENDDO

             ! SHUFFLE if necessary to remove zero-density cohorts
             IF (nc.LT.pop%pop_grid(j)%patch(k)%Layer(1)%ncohort) THEN
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ivec(1:nc))
                pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = nc

                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%density = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%id = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%biomass = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood_area = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%heartwood = 0.0
             ENDIF

             pop%pop_grid(j)%patch(k)%age(idisturb) = 0
             pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0


          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE Patch_partial_disturb2
  !*******************************************************************************

  SUBROUTINE Patch_disturb(pop,idisturb,precip)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT)  :: POP
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    !INTEGER(i4b), INTENT(IN) :: it(:),idisturb
    INTEGER(i4b), INTENT(IN) :: idisturb
    INTEGER(i4b) :: j, k, np, nc

    np = SIZE(Pop%pop_grid)
    ! Kills all biomass in patch when prescribed disturbance interval is reached
    ! Should be called after accounting for this year

    DO j=1,np
       DO k=1,NPATCH2D
          pop%pop_grid(j)%patch(k)%cat_mortality = 0.
          IF (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).NE.0) THEN
             IF ((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb)).or. &
                  (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb)) ) THEN
                ! kill entire layer
                nc = pop%pop_grid(j)%patch(k)%layer(1)%ncohort

                ! pop%pop_grid(j)%patch(k)%fire_mortality = SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%biomass)
                pop%pop_grid(j)%patch(k)%cat_mortality = SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%biomass)
                pop%pop_grid(j)%patch(k)%sapwood_loss =  pop%pop_grid(j)%patch(k)%sapwood_loss + &
                     SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%sapwood)
                pop%pop_grid(j)%patch(k)%sapwood_area_loss =  pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                     SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%sapwood_area)
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%ncohort = 0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%biomass = 0.0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%density = 0.0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmean = 0.0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmax  = 0.0
                pop%pop_grid(j)%patch(k)%age(idisturb) = 0
                pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0

                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%density = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%id = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%biomass = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%sapwood = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%sapwood_area = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%heartwood = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%age = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%lai = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%height = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%pgap = 1.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_resource_uptake = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_light_uptake = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_interception = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_respiration = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_NPP = 0.0
                pop%pop_grid(j)%patch(k)%area_growth = 0.0
                pop%pop_grid(j)%patch(k)%pgap = 1.0
                ! understorey recruitment
                IF (PRESENT(precip)) THEN
                   CALL layer_recruitment_single_patch(pop,k,j,precip)
                ELSE
                   CALL layer_recruitment_single_patch(pop,k,j)

                ENDIF
             ENDIF
          ELSEIF (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb)) THEN
             ! kill entire layer
             nc = pop%pop_grid(j)%patch(k)%layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%sapwood_loss =  pop%pop_grid(j)%patch(k)%sapwood_loss + &
                  SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%sapwood)
             pop%pop_grid(j)%patch(k)%sapwood_area_loss =  pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                  SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%sapwood_area)
             pop%pop_grid(j)%patch(k)%cat_mortality = SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%biomass)
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%ncohort = 0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%biomass = 0.0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%density = 0.0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmean = 0.0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmax  = 0.0
             pop%pop_grid(j)%patch(k)%age(idisturb) = 0
             pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%density = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%id = 0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%biomass = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%sapwood = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%sapwood_area = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%heartwood = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%age = 0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%lai = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%height = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%pgap = 1.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_resource_uptake = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_light_uptake = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_interception = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_respiration = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_NPP = 0.0
             pop%pop_grid(j)%patch(k)%area_growth = 0.0
             pop%pop_grid(j)%patch(k)%pgap = 1.0
             ! understorey recruitment
             IF (PRESENT(precip)) THEN
                CALL layer_recruitment_single_patch(pop,k,j,precip)
             ELSE
                CALL layer_recruitment_single_patch(pop,k,j)

             ENDIF
          ENDIF

       ENDDO

    ENDDO

  END SUBROUTINE Patch_disturb
  !*******************************************************************************
  SUBROUTINE  layer_recruitment(pop,precip)
    IMPLICIT NONE
    TYPE(POP_TYPE), INTENT(INOUT)  :: POP
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    REAL(dp) :: f, mu, densindiv, cmass, ht
    REAL(dp) :: tmp, cmass_stem_sum,cmass_stem_inc
    INTEGER(i4b) :: j, k,c, ncohort, np
    REAL(dp) :: diam,basal

    np = SIZE(Pop%pop_grid)

    DO j=1,np
       DO k=1,NPATCH2D
          IF (RECRUIT_SWITCH==0) THEN
             pop%pop_grid(j)%patch(k)%factor_recruit = EXP(-0.6*((pop%pop_grid(j)%patch(k)%Layer(1)%biomass)**(0.6667)))
          ELSEIF (RECRUIT_SWITCH==1) THEN
             pop%pop_grid(j)%patch(k)%factor_recruit = max(pop%pop_grid(j)%patch(k)%pgap,1e-3)
          ENDIF
          f = pop%pop_grid(j)%patch(k)%factor_recruit
          mu=EXP(max(FULTON_ALPHA*(1.0-2*THETA_recruit/(f+1-SQRT((f+1)*(f+1)-4*THETA_recruit*f))),-50.0));
          densindiv=DENSINDIV_MAX*mu;
          cmass=CMASS_STEM_INIT*densindiv/DENSINDIV_MAX;

          !COMMLN below: should not be cohort +1 or .LE. !
          IF (cmass>EPS*10..AND.densindiv>DENSINDIV_MIN.AND.(pop%pop_grid(j)%patch(k)%Layer(1)%ncohort+1).LT.NCOHORT_MAX) THEN

             ! create a new cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort + 1
             ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%biomass = cmass
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%density = densindiv
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%sapwood = cmass

             IF ( PRESENT(precip) ) THEN
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass, densindiv, ht, diam, basal, precip(j))
             ELSE
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass, densindiv, ht, diam, basal )
             ENDIF

             pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%height   = ht
             pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%diameter = diam

          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE layer_recruitment

  !*******************************************************************************
  SUBROUTINE  layer_recruitment_single_patch(pop, index, grid_index,precip)
    IMPLICIT NONE
    TYPE(POP_TYPE), INTENT(INOUT)  :: POP
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    INTEGER(i4b), INTENT(IN) :: index, grid_index
    REAL(dp) :: f, mu, densindiv, cmass, ht
    REAL(dp) :: tmp, cmass_stem_sum,cmass_stem_inc
    INTEGER(i4b) :: j, k,c, ncohort, np
    REAL(dp) :: diam,basal

    np = SIZE(Pop%pop_grid)
    DO j=grid_index,grid_index
       DO k=index,index
          IF (RECRUIT_SWITCH==0) THEN
             pop%pop_grid(j)%patch(k)%factor_recruit = EXP(-0.6*((pop%pop_grid(j)%patch(k)%Layer(1)%biomass)**(0.6667)))
          ELSEIF (RECRUIT_SWITCH==1) THEN
             !pop%pop_grid(j)%patch(k)%factor_recruit = pop%pop_grid(j)%patch(k)%pgap
             pop%pop_grid(j)%patch(k)%factor_recruit = 1
          ENDIF
          f = pop%pop_grid(j)%patch(k)%factor_recruit
          mu=EXP(FULTON_ALPHA*(1.0-2*THETA_recruit/(f+1-SQRT((f+1)*(f+1)-4*THETA_recruit*f))));
          densindiv=DENSINDIV_MAX*mu;
          cmass=CMASS_STEM_INIT*densindiv/DENSINDIV_MAX;
          ! write(*,*) 'layer_recruitment_single_patch', densindiv, pop%pop_grid(j)%patch(k)%Layer(1)%biomass
          IF (cmass>EPS*10..AND.densindiv>DENSINDIV_MIN.AND.(pop%pop_grid(j)%patch(k)%Layer(1)%ncohort+1).LT.NCOHORT_MAX) THEN
             ! create a new cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort + 1
             ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%biomass = cmass
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%density = densindiv
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%sapwood = cmass

             IF ( PRESENT(precip) ) THEN
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass, densindiv, ht, diam, basal, precip(j))
             ELSE
                CALL GET_ALLOMETRY( ALLOM_SWITCH,  cmass, densindiv, ht, diam, basal )
             ENDIF

             pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%height   = ht
             pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%diameter = diam


          ENDIF


       ENDDO
    ENDDO


  END SUBROUTINE layer_recruitment_single_patch

  !*******************************************************************************
  ! Exponential distribution
  ! Returns probability of a given time-between-events (x)
  ! Given a Poisson process with expected frequency (events per unit time) lambda
  ! Reference: http://en.wikipedia.org/wiki/Exponential_distribution
  ! Use to determine average age (x, years) of patches with a given random disturbance
  ! frequency lambda (disturbances per year)

  REAL(dp) FUNCTION Exponential(lambda, x)
    IMPLICIT NONE
    INTEGER(i4b), INTENT(IN) :: x
    REAL(dp), INTENT(IN) ::  lambda

    IF (x.LT.0) THEN ! Shouldn't happen but ...
       Exponential=0.0
    ELSE
       Exponential=lambda*EXP(-lambda*x)
    ENDIF

  END FUNCTION Exponential
  !*******************************************************************************
  ! Exponential distribution
  ! Returns probability of a given time-between-events (x)
  ! Given a Poisson process with expected frequency (events per unit time) lambda
  ! Reference: http://en.wikipedia.org/wiki/Exponential_distribution
  ! Use to determine average age (x, years) of patches with a given random disturbance
  ! frequency lambda (disturbances per year)

  REAL(dp) FUNCTION REALExponential(lambda, x)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) ::  x
    REAL(dp), INTENT(IN) ::  lambda

    IF (x.LT.0) THEN ! Shouldn't happen but ...
       REALExponential=0.0
    ELSE
       REALExponential=lambda*EXP(-lambda*x)
    ENDIF

  END FUNCTION REALExponential



  !*******************************************************************************
  REAL(dp) FUNCTION CumExponential(lambda, x)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) ::  lambda

    IF (x.LT.0) THEN ! Shouldn't happen but ...
       CumExponential=0.0
    ELSE
       CumExponential=1.-EXP(-lambda*x)
    ENDIF

  END FUNCTION CumExponential


  !*******************************************************************************

  REAL(dp)  FUNCTION Factorial(n)
    IMPLICIT NONE
    INTEGER , INTENT(IN) :: n
    INTEGER :: i
    REAL(dp) :: Ans

    Ans = 1.
    DO i = 1, n
       Ans = Ans *REAL(i,dp)
    END DO

    Factorial = Ans

  END FUNCTION Factorial
  !*******************************************************************************
  ! ALLOMETRY
  !*******************************************************************************
  SUBROUTINE GET_ALLOMETRY( ALLOM_SWITCH,  biomass, density, ht, diam, basal, precip )

    IMPLICIT NONE
    INTEGER(i4b), INTENT(IN) :: ALLOM_SWITCH
    REAL(dp),     INTENT(IN) :: biomass
    REAL(dp),     INTENT(IN) :: density
    REAL(dp),     INTENT(IN), OPTIONAL :: precip
    REAL(dp),     INTENT(OUT):: ht, diam, basal

    ! Standard Allometry
    IF (ALLOM_SWITCH.EQ.0) THEN
       ht   = (Kbiometric**(3.0/4.0))*(4.*biomass/(max(density,1e-5)*WD*PI))**(1.0/4.0)
       diam = (ht/Kbiometric)**(1.5)
       basal= PI * (diam/2.0) * (diam/2.0) * density * 1.e4

       ! Top-End Allometry following G.Cook
    ELSEIF (ALLOM_SWITCH.EQ.1.AND.PRESENT(precip)) THEN
       ht   =GetHeight(precip,biomass,density)
       CALL Allometry(ht,biomass,density,diam,basal)

       ! Allometry following Williams 2005, Model 5b
    ELSEIF ( ALLOM_SWITCH.EQ.2 ) THEN
       CALL Williams_Allometry(biomass,density,ht,diam,basal)

    ELSE
       WRITE(*,*)"Invalid Allometry settings in POP!"
       WRITE(*,*)"ALLOM_SWITCH   = ",ALLOM_SWITCH
       WRITE(*,*)"Precip present = ",PRESENT(precip)
       STOP -1
    ENDIF

  END SUBROUTINE GET_ALLOMETRY
  !*******************************************************************************
  ! TOP-END ALLOMETRY STARTS HERE
  !*******************************************************************************
  ! Tree height based on precipitation and Gary Cook Top-End allometry
  ! Bisection solution for tree height (m) based on modified height-DBH relationship
  ! from Garry Cook (pers. comm. 15/4/2013)
  ! "I have been using H=0.054xExp(0.0014xRF)xD + 4.05*exp(-0.00032*rf) with Rf in mm, D in cm and height in m."
  ! Since the above expression does not go to zero at diameter=0 it is linked to a simple linear equation
  ! for initial height growth (H=50*D with D in m) using a non-rectangular hyperbola to smooth between the two
  ! Mathematical derivation: see POP documentation
  ! Arguments:
  ! precip  = annual precipitation (mm)
  ! biomass = tree stem C biomass across patch (kgC/m2)
  ! density = tree density (indiv/m2)

  REAL(dp) FUNCTION GetHeight(precip,biomass,density)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: precip
    REAL(dp), INTENT(IN) :: biomass
    REAL(dp), INTENT(IN) :: density

    REAL(dp),PARAMETER:: THETA=0.99 ! Shape parameter, should be slightly <1
    REAL(dp),PARAMETER:: HMIN=0.001 ! min bound for tree height
    REAL(dp),PARAMETER:: HMAX=100 ! max bound for tree height
    REAL(dp),PARAMETER:: EPS=0.01 ! precision of the root
    INTEGER(i4b), PARAMETER :: MAXTRIES=25

    REAL(dp) :: alpha,beta,delta,rh,st,x1,x2,rtbis,dx,fmid,xmid,lhs,rhs
    INTEGER(i4b) :: b

    alpha=4.05*EXP(-0.00032*precip)
    beta=5.4*EXP(0.0014*precip)
    delta=2.0*SQRT(biomass/density/WD/PI)

    x1=HMIN
    x2=HMAX
    rtbis=x1
    dx=x2-x1
    b=0
    fmid=EPS+1.0

    DO WHILE (ABS(dx).GT.EPS.AND.b.LE.MAXTRIES)
       b=b+1
       dx=dx*0.5
       xmid=rtbis+dx

       ! Evaluate LHS-RHS at height=xmid
       ! LHS-RHS should increase with increasing height

       lhs=xmid
       rh=1.0/SQRT(xmid)
       st=alpha+beta*delta*rh+100*delta*rh
       rhs=1.0/2.0/THETA* &
            (st-SQRT(st*st-400*THETA*alpha*delta*rh- &
            400*THETA*beta*delta*delta/xmid))
       fmid=lhs-rhs

       IF (fmid.LT.0.0) rtbis=xmid

    ENDDO

    GetHeight=xmid

  END FUNCTION GetHeight
!*******************************************************************************
SUBROUTINE INTERPOLATE_BIOMASS_1D(pop, disturbance_interval,it,g)
IMPLICIT NONE

TYPE(POP_TYPE), INTENT(INOUT) :: POP
INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:)
INTEGER(i4b), INTENT(IN) ::  it,g

INTEGER(i4b) :: nage,iage, i_min, i_max, tmp_array(NPATCH2D)
INTEGER(i4b) :: i_min_growth, i_max_growth
REAL(dp) :: disturbance_freq,tmp_min,tmp_max, tmp1_min, tmp1_max
REAL(dp) :: tmp2_min, tmp2_max
REAL(dp) :: tmp3_min, tmp3_max
REAL(dp) :: tmp4_min, tmp4_max
LOGICAL :: MASK(NPATCH2D)
INTEGER(i4b) :: age_min, age_max
INTEGER(i4b) :: age_min_growth, age_max_growth
INTEGER(i4b), ALLOCATABLE :: age(:)
REAL(dp), ALLOCATABLE ::cmass_age(:), stress_mort_age(:), crowd_mort_age(:)
REAL(dp), ALLOCATABLE ::csapwood_age(:), sapwood_area_age(:), growth_age(:)
REAL(dp), ALLOCATABLE ::freq_age(:)

! get interpolated biomass,sapwood, stress mortality, crowding mortality, disturbance mortality
POP%pop_grid(g)%cmass_sum= 0
POP%pop_grid(g)%stress_mortality = 0
POP%pop_grid(g)%cat_mortality = 0
pop%pop_grid(g)%crowding_mortality = 0
pop%pop_grid(g)%csapwood_sum = 0
pop%pop_grid(g)%sapwood_area = 0
tmp_array = 0
nage =  min(POP%pop_grid(g)%patch(1)%disturbance_interval(1),it)+1 ! maximum age
!nage = maxval(pop%pop_grid(g)%patch(:)%age(1))
disturbance_freq=1.0/REAL(disturbance_interval(g,1))
IF(.NOT.ALLOCATED(age)) ALLOCATE(age(nage))
IF(.NOT.ALLOCATED(freq_age)) ALLOCATE(freq_age(nage))
IF(.NOT.ALLOCATED(cmass_age)) ALLOCATE(cmass_age(nage))
IF(.NOT.ALLOCATED(growth_age)) ALLOCATE(growth_age(nage))
IF(.NOT.ALLOCATED(csapwood_age)) ALLOCATE(csapwood_age(nage))
IF(.NOT.ALLOCATED(sapwood_area_age)) ALLOCATE(sapwood_area_age(nage))
IF(.NOT.ALLOCATED(stress_mort_age)) ALLOCATE(stress_mort_age(nage))
IF(.NOT.ALLOCATED(crowd_mort_age)) ALLOCATE(crowd_mort_age(nage))
!pop%pop_grid(g)%biomass_age(2:agemax) = pop%pop_grid(g)%biomass_age(1:agemax-1)
!pop%pop_grid(g)%biomass_age(1) = 0.0
!cmass_age = pop%pop_grid(g)%biomass_age
tmp_min = 0.0
tmp_max = 0.0
pop%pop_grid(g)%biomass_age = 0.0


!IF (pop%LU(g)==2) THEN  ! secondary forest
IF (POP%pop_grid(g)%LU==2) then ! secondary forest
   DO iage = 1, nage
      age(iage) = iage-1
      freq_age(iage) =  pop%pop_grid(g)%freq_age(iage)
   ENDDO
ELSE
   DO iage = 1, nage
      age(iage) = iage-1
      freq_age(iage) =  REALExponential(disturbance_freq,REAL(age(iage),dp))
      pop%pop_grid(g)%freq_age(iage) = freq_age(iage)
     
   END DO

ENDIF
if (sum(freq_age)>0.0) freq_age = freq_age/sum(freq_age)

DO iage = 1, nage
   ! get nearest ages bracketing age(iage)
   if (any(pop%pop_grid(g)%patch(:)%age(1).LE.age(iage))) then
      age_min = MAXVAL(pop%pop_grid(g)%patch(:)%age(1), &
           pop%pop_grid(g)%patch(:)%age(1).LE.age(iage))
      i_min = MAXLOC(pop%pop_grid(g)%patch(:)%age(1), 1, &
           pop%pop_grid(g)%patch(:)%age(1).LE.age(iage))
   else
      age_min = 0
      i_min = 0
   endif
   if (any(pop%pop_grid(g)%patch(:)%age(1).GE.age(iage))) then
      age_max = MINVAL(pop%pop_grid(g)%patch(:)%age(1), &
           pop%pop_grid(g)%patch(:)%age(1).GE.age(iage))
      i_max = MINLOC(pop%pop_grid(g)%patch(:)%age(1), 1, &
           pop%pop_grid(g)%patch(:)%age(1).GE.age(iage))
   else
      age_max = 0
      i_max = 0
   endif

   age_min_growth = age_min
   age_max_growth = age_max
   i_min_growth = i_min
   i_max_growth = i_max
!write(*,*) 'I0',it, age(iage), age_min, age_max, i_min, i_max 
   if ((i_min.gt.0).and.(i_max.gt.0).and.(age_max.eq.age_min)) then
      ! no need to interpolate

      MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_min
      where (MASK)
         tmp_array= 1
      elsewhere
         tmp_array=0
      endwhere
      cmass_age(iage) =  &
           SUM(pop%pop_grid(g)%patch(:)%layer(1)%biomass,MASK)/SUM(tmp_array)
      growth_age(iage) =  &
           SUM(pop%pop_grid(g)%patch(:)%growth,MASK)/SUM(tmp_array)
      csapwood_age(iage) = SUM(pop%pop_grid(g)%patch(:)%sapwood,MASK)/SUM(tmp_array)
      sapwood_area_age(iage) = & 
           SUM(pop%pop_grid(g)%patch(:)%sapwood_area,MASK)/SUM(tmp_array)
      stress_mort_age(iage)= &
           SUM(pop%pop_grid(g)%patch(:)%stress_mortality,MASK)/SUM(tmp_array)
      crowd_mort_age(iage)= &
           SUM(pop%pop_grid(g)%patch(:)%crowding_mortality,MASK)/SUM(tmp_array)
   else
      ! interpolate or extrapolate
      if ((i_min.eq.0).and.(i_max.gt.0)) then
         ! interpolate to zero
         age_min = 0
         i_min = 0
      elseif  ((i_max.eq.0).and.(i_min.gt.0)) then
         ! extrapolate to higher age
         age_max = age_min
         i_max = i_min
         age_min = MAXVAL(pop%pop_grid(g)%patch(:)%age(1), &
              pop%pop_grid(g)%patch(:)%age(1).LT.age_max)
         i_min = MAXLOC(pop%pop_grid(g)%patch(:)%age(1),1, &
              pop%pop_grid(g)%patch(:)%age(1).LT.age_max)
      endif

      ! interpolate or extrapolate (growth)
      if ((i_min_growth.eq.0).and.(i_max_growth.gt.0).and.age(iage).LE.2) then
         ! interpolate to zero
         age_min_growth = 0
         i_min_growth = 0
      elseif (((age_min_growth.LE.2).OR.(i_min_growth.eq.0)).and. &
           (i_max_growth.gt.0).and.age(iage).GT.2) then
         ! extrapolate to lower age
         age_min_growth = age_max_growth
         i_min_growth = i_max_growth

         age_max_growth = MINVAL(pop%pop_grid(g)%patch(:)%age(1), &
              pop%pop_grid(g)%patch(:)%age(1).GT.age_min_growth)
         i_max_growth = MINLOC(pop%pop_grid(g)%patch(:)%age(1), 1, &
              pop%pop_grid(g)%patch(:)%age(1).GT.age_min_growth)
      elseif  ((i_max_growth.eq.0).and.(i_min_growth.gt.0)) then
         ! extrapolate to higher age
         age_max_growth = age_min_growth
         i_max_growth = i_min_growth

         age_min_growth = MAXVAL(pop%pop_grid(g)%patch(:)%age(1), &
              pop%pop_grid(g)%patch(:)%age(1).LT.age_max_growth)
         i_min_growth = MAXLOC(pop%pop_grid(g)%patch(:)%age(1),1, &
              pop%pop_grid(g)%patch(:)%age(1).LT.age_max_growth)
      endif
     

      if (i_min.ne.0.and.age_min.ne.0) then
         MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_min
         where (MASK)
            tmp_array= 1
         elsewhere
            tmp_array=0
         endwhere
         tmp_min = SUM(pop%pop_grid(g)%patch(:)%layer(1)%biomass,MASK)/SUM(tmp_array)
         tmp1_min = SUM(pop%pop_grid(g)%patch(:)%stress_mortality,MASK)/SUM(tmp_array)
         tmp2_min = SUM(pop%pop_grid(g)%patch(:)%crowding_mortality,MASK)/SUM(tmp_array)
         tmp3_min = SUM(pop%pop_grid(g)%patch(:)%sapwood,MASK)/SUM(tmp_array)
         tmp4_min = SUM(pop%pop_grid(g)%patch(:)%sapwood_area,MASK)/SUM(tmp_array)
      else
         tmp_min = 0.0
         tmp1_min = 0.0
         tmp2_min = 0.0
         tmp3_min = 0.0
         tmp4_min = 0.0
      endif

      MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_max
      where (MASK)
         tmp_array= 1
      elsewhere
         tmp_array=0
      endwhere
      tmp_max = SUM(pop%pop_grid(g)%patch(:)%layer(1)%biomass,MASK)/SUM(tmp_array)
      tmp1_max = SUM(pop%pop_grid(g)%patch(:)%stress_mortality,MASK)/SUM(tmp_array)
      tmp2_max = SUM(pop%pop_grid(g)%patch(:)%crowding_mortality,MASK)/SUM(tmp_array)
      tmp3_max = SUM(pop%pop_grid(g)%patch(:)%sapwood,MASK)/SUM(tmp_array)
      tmp4_max = SUM(pop%pop_grid(g)%patch(:)%sapwood_area,MASK)/SUM(tmp_array)

      cmass_age(iage) = tmp_min + (tmp_max-tmp_min)/real(age_max-age_min)* &
           real(age(iage)-age_min)

      stress_mort_age(iage) = tmp1_min + (tmp1_max-tmp1_min)/real(age_max-age_min)* &
           real(age(iage)-age_min)

      crowd_mort_age(iage) = tmp2_min + (tmp2_max-tmp2_min)/real(age_max-age_min)* &
           real(age(iage)-age_min)

      csapwood_age(iage) = tmp3_min + (tmp3_max-tmp3_min)/real(age_max-age_min)* &
           real(age(iage)-age_min)

      sapwood_area_age(iage) = tmp4_min + (tmp4_max-tmp4_min)/real(age_max-age_min)* &
           real(age(iage)-age_min)


      if (i_min_growth.ne.0.and.age_min_growth.ne.0) then
         MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_min_growth
         where (MASK)
            tmp_array= 1
         elsewhere
            tmp_array=0
         endwhere
         tmp_min = SUM(pop%pop_grid(g)%patch(:)%growth,MASK)/SUM(tmp_array)
      else
         tmp_min = 0.0
      endif

      MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_max_growth
      where (MASK)
         tmp_array= 1
      elsewhere
         tmp_array=0
      endwhere
      tmp_max = SUM(pop%pop_grid(g)%patch(:)%growth,MASK)/SUM(tmp_array)

      growth_age(iage) = tmp_min + (tmp_max-tmp_min)/real(age_max_growth-age_min_growth)* &
           real(age(iage)-age_min_growth)

      !cmass_age(iage) = cmass_age(iage) + growth_age(iage) -  stress_mort_age(iage) - crowd_mort_age(iage)
   endif
   ! write(*,*) 'age_min, age_max 1:', age(iage), age_min, age_max, i_min, i_max

   !if (it.eq.403) then
   !   witer(58,*) age(iage),age_min, age_max, tmp1_min, tmp1_max,  stress_mort_age(iage), i_min
   !endif

   POP%pop_grid(g)%cmass_sum =  POP%pop_grid(g)%cmass_sum + &
        freq_age(iage)*cmass_age(iage)
   POP%pop_grid(g)%csapwood_sum =  POP%pop_grid(g)%csapwood_sum + &
        freq_age(iage)*csapwood_age(iage)
   POP%pop_grid(g)%sapwood_area =  POP%pop_grid(g)%sapwood_area + &
        freq_age(iage)*sapwood_area_age(iage)

!!$if (g==2) then 
!!$write(71, "(2i4, 350e16.6)")  it,  iage, freq_age(iage), cmass_age(iage), growth_age(iage),  stress_mort_age(iage), &
!!$ crowd_mort_age(iage), tmp_min, tmp_max, real(age_max_growth), real(age_min_growth)
!!$endif 

!!$if (g==2) write(72, "(2i4, 350e16.6)")  it,  iage, freq_age(iage), cmass_age(iage)
!!$if (g==1) write(71, "(2i4, 350e16.6)")  it,  iage, freq_age(iage), cmass_age(iage)
   POP%pop_grid(g)%stress_mortality =  POP%pop_grid(g)%stress_mortality + &
        freq_age(iage)*stress_mort_age(iage)
   POP%pop_grid(g)%crowding_mortality =  POP%pop_grid(g)%crowding_mortality + &
        freq_age(iage)*crowd_mort_age(iage)

   POP%pop_grid(g)%cat_mortality = POP%pop_grid(g)%growth - &
        POP%pop_grid(g)%stress_mortality - &
        POP%pop_grid(g)%crowding_mortality - &
        ( POP%pop_grid(g)%cmass_sum- POP%pop_grid(g)%cmass_sum_old)

   POP%pop_grid(g)%fire_mortality = 0.0
   POP%pop_grid(g)%res_mortality = 0.0

   pop%pop_grid(g)%biomass_age(iage) = cmass_age(iage)

enddo
!!$if (it.gt.400) then
!!$   write(*,*) 'it, nage', it, nage
!!$   write(591, "(350e16.6)") freq_age
!!$   write(601,"(350e16.6)") cmass_age
!!$   write(602,"(350e16.6)") stress_mort_age
!!$   write(603,"(350e16.6)") real(age)
!!$endif
DEALLOCATE(age)
DEALLOCATE(freq_age)
DEALLOCATE(cmass_age)
DEALLOCATE(stress_mort_age)




END SUBROUTINE INTERPOLATE_BIOMASS_1D
!*******************************************************************************
SUBROUTINE INTERPOLATE_BIOMASS_2D(pop, disturbance_interval,it,g)
IMPLICIT NONE

TYPE(POP_TYPE), INTENT(INOUT) :: POP
INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:)
INTEGER(i4b), INTENT(IN) ::  it,g
INTEGER(i4b), allocatable :: A1(:), A2(:) ! interpolated ages
INTEGER(i4b), allocatable :: xobs(:), yobs(:)     ! observed ages
REAL(dp), allocatable :: z1obs(:), z2obs(:), z3obs(:) ! observed biomass, stress_mort, crowd_mort
REAL(dp), allocatable :: z1interp(:), z2interp(:), z3interp(:) ! interpolated biomass, stress mortality, crowding mortality
REAL(dp), allocatable :: freq_interp(:) ! weightings for interpolated age pairs
REAL(dp), allocatable :: zp(:)  ! euclidean distance from interpolated age pair to observed age pairs
INTEGER(i4b) :: age_max(2), nrep(NPATCH2D+1)
INTEGER(i4b) :: tmp1, tmp2, I1, I2,I3,I4
INTEGER(i4b) :: x1, x2, x3, x4, y1, y2, y3, y4
INTEGER(i4b) :: p, j, k, n, np, nobs, count_extrap,l, ct
LOGICAL :: flag
REAL(dp) :: biomass(NPATCH2D+1), stress_mort(NPATCH2D+1), crowd_mort(NPATCH2D+1)
INTEGER(i4b) :: age1(NPATCH2D+1),  age2(NPATCH2D+1)
REAL(dp) :: zmin
INTEGER(i4b), allocatable :: interp_case(:), tmp_array(:), tmp(:)
REAL(dp) :: area(4,4), x(4), y(4),  disturbance_freq1,  disturbance_freq2
INTEGER(i4b) :: triangle_points(4,3), I_inside_triangle, Ineighbour(8)
LOGICAL ::  MASK_INSIDE_TRIANGLE(4), IS_NEIGHBOUR(8), tmp_logical
LOGICAL, allocatable :: MASK2(:), MASK3(:), MASK4(:)
INTEGER(i4b), allocatable :: address(:,:)

 POP%pop_grid(g)%cmass_sum = 0.0
 POP%pop_grid(g)%stress_mortality = 0.0
 POP%pop_grid(g)%crowding_mortality = 0.0

! Construct Age  Interpolating Grid
age_max(1) =  min(POP%pop_grid(g)%patch(1)%disturbance_interval(1),it)+1 ! maximum age
age_max(2) =  min(POP%pop_grid(g)%patch(1)%disturbance_interval(NDISTURB),it)+1 ! maximum age

p=0;
DO j=0,age_max(1)
   DO k=0, age_max(2)
      IF (k>j) THEN
         p=p+1
      ENDIF
   ENDDO
ENDDO
np = p

ALLOCATE(A1(np))
ALLOCATE(A2(np))
ALLOCATE(z1interp(np))
ALLOCATE(z2interp(np))
ALLOCATE(z3interp(np))
ALLOCATE(freq_interp(np))
ALLOCATE(interp_case(np))
ALLOCATE(tmp_array(np))
ALLOCATE(tmp(np))
ALLOCATE(address(age_max(1)+1,age_max(2)+1))
p=0;
address = -9999
DO j=0,age_max(1)
   DO k=0, age_max(2)
      IF (k>j) THEN
         p=p+1
         A1(p) = j
         A2(p) = k
         address(j+1,k+1) = p
      ENDIF
   ENDDO
ENDDO

! Construct Age observations

nrep = 1
flag = .false.
j = 1
! create "observed" age pair (0,0) with zero biomass
stress_mort(j) = 0.0
biomass(j) = 0.0
crowd_mort(j) = 0.0
age1(j) = 0
age2(j) = 0
j = j+1


do k=1,NPATCH2D
   tmp1 = pop%pop_grid(g)%patch(k)%age(1)
   tmp2 = pop%pop_grid(g)%patch(k)%age(NDISTURB)
   if (j.gt.1) then
      do n=1,j-1
         if (tmp1 == age1(n) .and. tmp2==age2(n)) then
            flag = .true.
            nrep(n) = nrep(n) + 1
            biomass(n) = biomass(n) + pop%pop_grid(g)%patch(k)%layer(1)%biomass
            stress_mort(n) = stress_mort(n) + pop%pop_grid(g)%patch(k)%stress_mortality
            crowd_mort(n) = crowd_mort(n) +  pop%pop_grid(g)%patch(k)%crowding_mortality
         endif
      enddo
   endif
   if (flag .eqv. .false.) then
       age1(j) = tmp1
       age2(j) = tmp2
       biomass(j) =  pop%pop_grid(g)%patch(k)%layer(1)%biomass
       stress_mort(j) =  pop%pop_grid(g)%patch(k)%stress_mortality
       crowd_mort(j) =  pop%pop_grid(g)%patch(k)%crowding_mortality
       if (k.ne.NPATCH2D) j = j+1
    else
       flag = .false.
    endif
enddo
nobs = j
biomass(1:nobs) = biomass(1:nobs)/nrep(1:nobs)
stress_mort(1:nobs) = stress_mort(1:nobs)/nrep(1:nobs)
crowd_mort(1:nobs) = crowd_mort(1:nobs)/nrep(1:nobs)
ALLOCATE(xobs(nobs))
ALLOCATE(yobs(nobs))
ALLOCATE(z1obs(nobs))
ALLOCATE(z2obs(nobs))
ALLOCATE(z3obs(nobs))
ALLOCATE(zp(nobs))
ALLOCATE(MASK2(nobs))
ALLOCATE(MASK3(nobs))
ALLOCATE(MASK4(nobs))
xobs = age1(1:nobs)
yobs = age2(1:nobs)
z1obs = biomass(1:nobs)
z2obs = stress_mort(1:nobs)
z3obs = crowd_mort(1:nobs)

if (it.eq.500) then
 do k=1,nobs
    write(502, "(2i5,2e16.6,i5)") xobs(k), yobs(k), z1obs(k), z3obs(k), nrep(k)
 enddo
endif

! get weightings for each interpolated age pair
do k = 1, np
   disturbance_freq1=1.0/REAL(disturbance_interval(g,1))
   disturbance_freq2=1.0/REAL(disturbance_interval(g,2))

   freq_interp(k) = REALExponential(disturbance_freq1,REAL(A1(k),dp)) * &
        REALExponential(disturbance_freq2,REAL(A2(k),dp))

ENDDO
freq_interp = freq_interp/sum(freq_interp)

! interpolate
DO p=1,np   ! loop over interpolated age pairs
   ! get distance to all observations
   DO j=1,nobs
      zp(j) = sqrt((real(A1(p))-real(xobs(j)))**2+(real(A2(p))-real(yobs(j)))**2)
   ENDDO

   ! get closest point
   zmin = MINVAL(zp)
   I1 = MINLOC(zp,1)
   x1 = xobs(I1)
   y1 = yobs(I1)


   ! check for obs locations forming a quadrangle around interpolating point
   MASK2 = (sign(1,A1(p)-xobs)== -sign(1,A1(p)-x1)).AND.(sign(1,A2(p)-yobs)== sign(1,A2(p)-y1).and.A1(p).NE.xobs)
   MASK3 = (sign(1,A1(p)-xobs)== sign(1,A1(p)-x1)).AND.(sign(1,A2(p)-yobs)== -sign(1,A2(p)-y1).and.A2(p).NE.yobs)
   MASK4 = (sign(1,A1(p)-xobs)== -sign(1,A1(p)-x1)).AND.(sign(1,A2(p)-yobs)== -sign(1,A2(p)-y1) &
        .and.A1(p).NE.xobs.and.A2(p).NE.yobs)

   IF ((ANY(MASK2)).and.(ANY(MASK3)).and.(ANY(MASK4)))    THEN
      ! get nearest point with opposing sign of x displacement
      I2 =  MINLOC(zp,1, MASK2)
      x2 = xobs(I2)
      y2 = yobs(I2)
      ! get nearest point with opposing sign of y displacement
      I3 =  MINLOC(zp,1, MASK3)
      x3 = xobs(I3)
      y3 = yobs(I3)
      ! get nearest point with opposing sign of x & y displacements
      I4 =  MINLOC(zp,1, MASK4)
      x4 = xobs(I4)
      y4 = yobs(I4)

    tmp_logical = .NOT.( (x2.eq.0.and.y2.eq.0).OR.(x3.eq.0.and.y3.eq.0).OR.(x4.eq.0.and.y4.eq.0))
   ENDIF




   IF ((A1(p)==x1).and.(A2(p)==y1)) THEN
      interp_case(p)=1
   ELSEIF ((ANY(MASK2)).and.(ANY(MASK3)).and.(ANY(MASK4)) .AND.tmp_logical) THEN ! quadrangle (without (0,0)) exists
      interp_case(p) = 2
   ELSE
      interp_case(p) = 3
   ENDIF!j=1



   select case (interp_case(p))

   case(1) ! interpolated point is the same as observation
      z1interp(p) = z1obs(I1)
      z2interp(p) = z2obs(I1)
      z3interp(p) = z3obs(I1)
   case(3) ! extrapolation required: set to value of nearest observation
      z1interp(p) = z1obs(I1)
      z2interp(p) = z2obs(I1)
      z3interp(p) = z3obs(I1)
   case(2)  ! quadrangle
      ! get nearest point with opposing sign of x displacement
      I2 =  MINLOC(zp,1, MASK2)
      x2 = xobs(I2)
      y2 = yobs(I2)
      ! get nearest point with opposing sign of y displacement
      I3 =  MINLOC(zp,1, MASK3)
      x3 = xobs(I3)
      y3 = yobs(I3)
      ! get nearest point with opposing sign of x & y displacements
      I4 =  MINLOC(zp,1, MASK4)
      x4 = xobs(I4)
      y4 = yobs(I4)


      x=real((/x1, x2, x3, x4/),dp)
      y =real((/y1, y2, y3, y4/),dp)

      ! get area of four possible triangles from 4 vertices, and corresponding 3 partial
      ! triangles, with observation as one vertex
      area = 0.0
      triangle_points = 0
      triangle_points(1,:) = (/I1, I2, I3/)
      area(1,1) = area_triangle(x(1), y(1),x(2), y(2),x(3), y(3))
      area(1,2) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(2), y(2), x(3), y(3))
      area(1,3) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(1), y(1), x(3), y(3))
      area(1,4) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(1), y(1), x(2), y(2))

      triangle_points(2,:) =  (/I1, I2, I4/)
      area(2,1) = area_triangle(x(1), y(1), x(2), y(2), x(4), y(4))
      area(2,2) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(2), y(2), x(4), y(4))
      area(2,3) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(1), y(1), x(4), y(4))
      area(2,4) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(1), y(1), x(2), y(2))

      triangle_points(3,:) =  (/I1, I4, I3/)
      area(3,1) = area_triangle(x(1) , y(1), x(3), y(3),x(4), y(4))
      area(3,2) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(4), y(4), x(3), y(3))
      area(3,3) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(1), y(1), x(3), y(3))
      area(3,4) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(1), y(1), x(4), y(4))

      triangle_points(4,:) =  (/I2, I4, I3/)
      area(4,1) = area_triangle(x(2), y(2), x(3), y(3), x(4), y(4));
      area(4,2) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(4), y(4), x(3), y(3))
      area(4,3) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(2), y(2), x(3), y(3))
      area(4,4) = area_triangle(real(A1(p),dp), real(A2(p),dp), x(2), y(2), x(4), y(4))

      MASK_INSIDE_TRIANGLE=(area(:,1)==sum(area(:,2:4),2));
      I_inside_triangle = MINLOC(area(:,1),1,MASK_INSIDE_TRIANGLE)
      z1interp(p) =SUM( z1obs(triangle_points(I_inside_triangle,:))* &
           area(I_inside_triangle,2:4))/sum(area(I_inside_triangle,2:4))
      z2interp(p) =SUM( z2obs(triangle_points(I_inside_triangle,:))* &
           area(I_inside_triangle,2:4))/sum(area(I_inside_triangle,2:4))
      z3interp(p) =SUM( z3obs(triangle_points(I_inside_triangle,:))* &
           area(I_inside_triangle,2:4))/sum(area(I_inside_triangle,2:4))


   case default
      write(*,*) " illegal interpolation case."
      stop
   end select ! interpolation case

ENDDO ! loop over interpolated age pairs

p = 0 ! counter for while loop
count_extrap = 0
! Extrapolation

IF (ANY(interp_case.NE.3)) THEN
   DO WHILE (ANY(interp_case==3))
      DO ct = 3,3
         DO p=1,np
            ! find number of neighbours for each extrapolable point
            is_neighbour = .FALSE.
            Ineighbour = 1
            tmp(p) = 0
            IF (interp_case(p)==3) THEN

               DO k=1,8
                  if(k==1) then
                     tmp1=A1(p)-1
                     tmp2=A2(p)
                  elseif(k==2) then
                     tmp1 = A1(p)+1
                     tmp2 =  A2(p)
                  elseif (k==3) then
                     tmp1 = A1(p)
                     tmp2 =  A2(p)-1
                  elseif (k==4) then
                     tmp1 = A1(p)
                     tmp2 =  A2(p)+1
                  elseif (k==5) then
                     tmp1=A1(p)-1
                     tmp2=A2(p)-1
                  elseif (k==6) then
                     tmp1=A1(p)+1
                     tmp2=A2(p)+1
                  elseif (k==7) then
                     tmp1=A1(p)+1
                     tmp2=A2(p)-1
                  elseif (k==8) then
                     tmp1=A1(p)-1
                     tmp2=A2(p)+1
                  endif
                  IF(tmp1.ge.0.and.tmp1.LE.age_max(1).and.tmp2.ge.0.and.tmp2.LE.age_max(2).and.tmp2.gt.tmp1) THEN
                     is_neighbour(k) = (interp_case(address(tmp1+1,tmp2+1)).ne.3)
                     Ineighbour(k) = address(tmp1+1,tmp2+1)
                  ENDIF

               ENDDO


               tmp(p) = COUNT(is_neighbour)
               if (tmp(p).GE.ct) then
                  z1interp(p) = sum(z1interp(Ineighbour) ,is_neighbour)/real(tmp(p))
                  z2interp(p) = sum(z2interp(Ineighbour) ,is_neighbour)/real(tmp(p))
                  z3interp(p) = sum(z3interp(Ineighbour) ,is_neighbour)/real(tmp(p))
                  interp_case(p) = 2
               endif
            ENDIF
         ENDDO !1, np
      ENDDO
      IF (ALL(tmp.lt.3)) THEN
         ! maximum of one neighbour
         p = MAXLOC(tmp,1)
         is_neighbour = .FALSE.
         Ineighbour = 1
         tmp(p) = 0
         DO k=1,8
            if(k==1) then
               tmp1=A1(p)-1
               tmp2=A2(p)
            elseif(k==2) then
               tmp1 = A1(p)+1
               tmp2 =  A2(p)
            elseif (k==3) then
               tmp1 = A1(p)
               tmp2 =  A2(p)-1
            elseif (k==4) then
               tmp1 = A1(p)
               tmp2 =  A2(p)+1
            elseif (k==5) then
               tmp1=A1(p)-1
               tmp2=A2(p)-1
            elseif (k==6) then
               tmp1=A1(p)+1
               tmp2=A2(p)+1
            elseif (k==7) then
               tmp1=A1(p)+1
               tmp2=A2(p)-1
            elseif (k==8) then
               tmp1=A1(p)-1
               tmp2=A2(p)+1
            endif
            IF(tmp1.ge.0.and.tmp1.LE.age_max(1).and.tmp2.ge.0.and.tmp2.LE.age_max(2).and.tmp2.gt.tmp1) THEN
               is_neighbour(k) = (interp_case(address(tmp1+1,tmp2+1)).ne.3)
               Ineighbour(k) = address(tmp1+1,tmp2+1)
            ENDIF

         ENDDO

         tmp(p) = COUNT(is_neighbour)
         z1interp(p) = sum(z1interp(Ineighbour) ,is_neighbour)/real(tmp(p))
         z2interp(p) = sum(z2interp(Ineighbour) ,is_neighbour)/real(tmp(p))
         z3interp(p) = sum(z3interp(Ineighbour) ,is_neighbour)/real(tmp(p))
         interp_case(p) = 2

      ENDIF  ! ALL(tmp.lt.2)

   ENDDO ! do while
ENDIF


DO p = 1,np
   POP%pop_grid(g)%cmass_sum =  POP%pop_grid(g)%cmass_sum + &
        freq_interp(p)*z1interp(p)
   POP%pop_grid(g)%stress_mortality =  POP%pop_grid(g)%stress_mortality + &
        freq_interp(p)*z2interp(p)
   POP%pop_grid(g)%crowding_mortality =  POP%pop_grid(g)%crowding_mortality + &
        freq_interp(p)*z3interp(p)


ENDDO

POP%pop_grid(g)%cat_mortality = POP%pop_grid(g)%cmass_sum *  disturbance_freq2
POP%pop_grid(g)%fire_mortality = POP%pop_grid(g)%growth - &
     POP%pop_grid(g)%cat_mortality - &
     POP%pop_grid(g)%stress_mortality - &
     POP%pop_grid(g)%crowding_mortality - &
     ( POP%pop_grid(g)%cmass_sum- POP%pop_grid(g)%cmass_sum_old)


DEALLOCATE(xobs)
DEALLOCATE(yobs)
DEALLOCATE(z1obs)
DEALLOCATE(z2obs)
DEALLOCATE(z3obs)
DEALLOCATE(zp)
DEALLOCATE(MASK2)
DEALLOCATE(MASK3)
DEALLOCATE(MASK4)
DEALLOCATE(A1)
DEALLOCATE(A2)
DEALLOCATE(z1interp)
DEALLOCATE(z2interp)
DEALLOCATE(z3interp)
DEALLOCATE(freq_interp)
DEALLOCATE(interp_case)
DEALLOCATE(tmp)
DEALLOCATE(tmp_array)
DEALLOCATE(address)

END SUBROUTINE INTERPOLATE_BIOMASS_2D

!******************************************************************************
SUBROUTINE SMOOTH_FLUX(POP,g,t)
IMPLICIT NONE

TYPE(POP_TYPE), INTENT(INOUT) :: POP
INTEGER(i4b), INTENT(IN) :: g, t
INTEGER(i4b), PARAMETER :: SPAN = NYEAR_SMOOTH/2
REAL(dp) :: x(SPAN+1), y(SPAN+1), a, b, r
REAL(dp) :: sumflux, sumsmooth, flux(NYEAR_HISTORY), smoothed_flux
REAL(dp) :: dbuf
INTEGER(i4b) :: t0, tt, n, k

! update fire_mortality_history

IF (t.gt.NYEAR_HISTORY) THEN
   DO k = 1, NYEAR_HISTORY-1
      POP%pop_grid(g)%fire_mortality_history(k) = POP%pop_grid(g)%fire_mortality_history(k+1)
   ENDDO
ENDIF
POP%pop_grid(g)%fire_mortality_history(t) = POP%pop_grid(g)%fire_mortality

flux = POP%pop_grid(g)%fire_mortality_history
t0 = t-SPAN
n = 0
sumflux = 0.0
sumsmooth = 0.0

DO tt = 1,NYEAR_SMOOTH
   IF ((t0+tt).GE.1 .AND. (t0+tt).LE.t+1) THEN
      sumflux = sumflux + flux(t0+tt)
      y(tt) = flux(t0+tt)
      x(tt) = tt
      n = n+1
      IF ((t0+tt).eq.t+1) THEN
         CALL regress(x,y,n,a,b,r)
      ENDIF
   ELSE
      sumflux = sumflux + (a + b*tt)
      n = n+ 1
   ENDIF
ENDDO

dbuf =POP%pop_grid(g)%smoothing_buffer/(real(NYEAR_SMOOTH)/2.0)
!MC smoothed_flux=max(sumflux/real(n)+dbuf, 0.0)
smoothed_flux=max(sumflux/real(n)+dbuf, 0.0_dp)
POP%pop_grid(g)%smoothing_buffer = POP%pop_grid(g)%smoothing_buffer + flux(t) - smoothed_flux

POP%pop_grid(g)%fire_mortality_smoothed = smoothed_flux



END SUBROUTINE SMOOTH_FLUX
!******************************************************************************
SUBROUTINE SMOOTH_FLUX_cat(POP,g,t)
IMPLICIT NONE

TYPE(POP_TYPE), INTENT(INOUT) :: POP
INTEGER(i4b), INTENT(IN) :: g, t
INTEGER(i4b), PARAMETER :: SPAN = NYEAR_SMOOTH/2
REAL(dp) :: x(SPAN+1), y(SPAN+1), a, b, r
REAL(dp) :: sumflux, sumsmooth, flux(NYEAR_HISTORY), smoothed_flux
REAL(dp) :: dbuf
INTEGER(i4b) :: t0, tt, n, k

! update fire_mortality_history

IF (t.gt.NYEAR_HISTORY) THEN
   DO k = 1, NYEAR_HISTORY-1
      POP%pop_grid(g)%cat_mortality_history(k) = POP%pop_grid(g)%cat_mortality_history(k+1)
   ENDDO
ENDIF
POP%pop_grid(g)%cat_mortality_history(t) = POP%pop_grid(g)%cat_mortality

flux = POP%pop_grid(g)%cat_mortality_history
t0 = t-SPAN
n = 0
sumflux = 0.0
sumsmooth = 0.0

DO tt = 1,NYEAR_SMOOTH
   IF ((t0+tt).GE.1 .AND. (t0+tt).LE.t+1) THEN
      sumflux = sumflux + flux(t0+tt-1)
      y(tt) = flux(t0+tt-1)
      x(tt) = tt
      n = n+1
      IF ((t0+tt).eq.t+1) THEN
         CALL regress(x,y,n,a,b,r)
      ENDIF
   ELSE
      sumflux = sumflux + (a + b*tt)
      n = n+ 1
   ENDIF
ENDDO

dbuf =POP%pop_grid(g)%smoothing_buffer_cat/(real(NYEAR_SMOOTH)/2.0)
!MC smoothed_flux=max(sumflux/real(n)+dbuf, 0.0)
smoothed_flux=max(sumflux/real(n)+dbuf, 0.0_dp)
POP%pop_grid(g)%smoothing_buffer_cat = POP%pop_grid(g)%smoothing_buffer_cat + flux(t) - smoothed_flux

POP%pop_grid(g)%cat_mortality_smoothed = smoothed_flux



END SUBROUTINE SMOOTH_FLUX_cat
!******************************************************************************
SUBROUTINE REGRESS(x, y, n, a, b, r)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: x(:), y(:)
REAL(dp), INTENT(OUT) :: a, b, r
INTEGER(i4b), INTENT(IN) :: n
REAL(dp)::  sx,sy,sxx,sxy,delta,meanx,meany,sdx,sdy
INTEGER(i4b) :: i

!Performs a linear regression of array y on array x (n values)
! returning parameters a and b in the fitted model: y=a+bx
! Source: Press et al 1986, Sect 14.2
!  also returns Pearson r


sx=0.0
sy=0.0
sxx=0.0
sxy=0.0
DO i=1,n
   sx = sx + x(i)
   sy = sy + y(i)
   sxx = sxx + x(i) * y(i)
   sxy = sxy + x(i)*y(i)
ENDDO
delta = real(n)*sxx - sx*sx
a=(sxx*sy-sx*sxy)/delta
b=(real(n)*sxy-sx*sy)/delta
meanx=sx/real(n)
meany=sy/real(n)
sdx = 0.0
sdy = 0.0
DO i=1,n
   sdx = sdx + (x(i)-meanx)*(x(i)-meanx)
   sdy = sdy + (y(i)-meany)*(y(i)-meany)
ENDDO
sdx=sqrt(sdx/real(n-1))
sdy=sqrt(sdy/real(n-1))
!MC if (sdx .eq. 0.0 .or. sdy .eq. 0.0) then
if ((abs(sdx) .lt. tiny(0.0_dp)) .or. (abs(sdy) .lt. tiny(0.0_dp))) then
   r = 0.0_dp
else
   r =b*sdx/sdy
ENDIF


END SUBROUTINE REGRESS
!******************************************************************************
REAL(dp) FUNCTION Area_Triangle(x1,y1,x2,y2,x3,y3)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: x1, y1, x2, y2, x3, y3

area_triangle = abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);

END FUNCTION Area_Triangle


!******************************************************************************
  ! Top-End Allometry
  ! Computes tree stem diameter (m) and basal area (m2/ha)
  ! given height (m), stem biomass (kgC/m2) and tree population density (indiv/m2)

  SUBROUTINE Allometry(height,biomass,density,diam,basal)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: height
    REAL(dp), INTENT(IN) :: biomass
    REAL(dp), INTENT(IN) :: density
    REAL(dp), INTENT(OUT) :: diam
    REAL(dp), INTENT(OUT) :: basal

    REAL(dp) :: delta,rh

    delta=2.0*SQRT(biomass/density/WD/PI)
    rh=1.0/SQRT(height)

    diam=delta*rh
    basal=PI*(diam/2.0)*(diam/2.0)*density*1e4

  END SUBROUTINE Allometry
  !*******************************************************************************
  SUBROUTINE Williams_Allometry(agBiomass, density, height, dbh, basal)

    IMPLICIT NONE
    ! Allometry following Williams 2005, Model 5b (see table 2)
    ! Williams et al., "Allometry for estimating aboveground tree biomass in tropical
    ! and subtropical eucalypt woodlands: towards general predictive equations",
    ! Australian Journal of Botany, 2005, 53, 607-619
    ! INPUT
    ! agbiomass: above ground biomass [kg(C)/m2]
    ! density  : tree population density [m-2]
    ! OUTPUT
    ! height   : tree height [m]
    ! dbh      : Diameter at breast height [m]
    ! basal    : Basal area [m2/ha]
    REAL(dp), INTENT(IN) :: agbiomass, density
    REAL(dp), INTENT(OUT):: height, dbh, basal
    REAL(dp), PARAMETER  :: beta0 = -2.3046
    REAL(dp), PARAMETER  :: beta1 = 2.5243
    REAL(dp), PARAMETER  :: gC2DM = 1./0.49  ! ratio Dry matter mass to g(C)

    ! Compute dbh using model 5b and converting from cm -> m
    dbh    = ( agbiomass * gC2DM / density / exp(beta0) ) ** (1./beta1)
    dbh    = dbh * 0.01
    ! Compute basal area [m^2/m^2]
    basal  = PI * (0.5 * dbh)**2 * density
    ! Compute Height using cylindrical approach
    height = agBiomass / ( WD * basal )
    ! Basal Area [m^2/m^2]->[m^2/ha]
    basal  = basal * 1.e4

  END SUBROUTINE Williams_Allometry

  !*******************************************************************************
  SUBROUTINE POP_init(POP, disturbance_interval,np,Iwood,precip)
    USE POP_types,     ONLY: POP_TYPE
    USE TypeDef,       ONLY: i4b
    IMPLICIT NONE

    INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:), np
    INTEGER(i4b), INTENT(IN) :: Iwood(:)
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    TYPE( POP_TYPE )        , INTENT(INOUT) :: POP
    INTEGER(i4b) :: j, k



 
    CALL alloc_POP(pop,int(np))
    POP%np = np
    POP%Iwood = Iwood
    !POP%LU = 1  ! initialise to primary forest
    POP%pop_grid(:)%LU =1 
    POP%it_pop = 0
    CALL ZeroPOP(pop)
    CALL InitPOP2D_Poisson(pop,INT(disturbance_interval,i4b))

    DO j=1,np
       DO k=1,NPATCH2D
          ! understorey recruitment
          IF (PRESENT(precip)) THEN
             CALL layer_recruitment_single_patch(pop,k,j,precip)
          ELSE
             CALL layer_recruitment_single_patch(pop,k,j)

          ENDIF
       ENDDO
    ENDDO



  END SUBROUTINE POP_init
  !*******************************************************************************
  SUBROUTINE POP_init_single(POP, disturbance_interval,n,precip)
    USE POP_types,     ONLY: POP_TYPE
    USE TypeDef,       ONLY: i4b
    IMPLICIT NONE

    INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:)
    INTEGER(i4b), INTENT(IN) :: n
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    TYPE( POP_TYPE )        , INTENT(INOUT) :: POP
    INTEGER(i4b) :: j, k



 
    
    
    
    POP%it_pop(n) = 0
    CALL ZeroPOP(pop,INT(n))   !mrd561: types do not match unless converted here
    CALL InitPOP2D_Poisson(pop,INT(disturbance_interval,i4b),n)

    DO j=n,n
       DO k=1,NPATCH2D
          ! understorey recruitment
          IF (PRESENT(precip)) THEN
             CALL layer_recruitment_single_patch(pop,k,j,precip)
          ELSE
             CALL layer_recruitment_single_patch(pop,k,j)

          ENDIF
       ENDDO
    ENDDO



  END SUBROUTINE POP_init_single

  !*******************************************************************************
  SUBROUTINE alloc_POP(POP, arraysize)

    USE TypeDef,   Only: i4b, dp
    USE POP_Types, Only: POP_TYPE


    TYPE( POP_TYPE ),INTENT(INOUT) :: POP
    INTEGER,            INTENT(IN) :: arraysize


    IF (.NOT.ALLOCATED(POP%POP_Grid)) ALLOCATE (POP%POP_Grid(arraysize))
    IF (.NOT.ALLOCATED(POP%Iwood)) ALLOCATE (POP%Iwood(arraysize))
    !IF (.NOT.ALLOCATED(POP%LU)) ALLOCATE (POP%LU(arraysize))
    IF (.NOT.ALLOCATED(POP%it_pop)) ALLOCATE (POP%it_pop(arraysize))

  END SUBROUTINE alloc_POP

  !*******************************************************************************
END MODULE POPModule
!*******************************************************************************
