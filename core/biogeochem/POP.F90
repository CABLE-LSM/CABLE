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

module TypeDef
  !-------------------------------------------------------------------------------
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

end module TypeDef


!*******************************************************************************


module POP_Constants

  use TYPEdef, only: dp, i4b

  implicit none

  ! REAL(dp),PARAMETER:: FULTON_ALPHA= 5.6 ! recruitment scalar alpha in Fulton (1991)
  ! REAL(dp),PARAMETER:: DENSINDIV_MAX=2  ! 0.5 !  Maximum density of individuals within a cohort indiv/m2
  ! REAL(dp),PARAMETER:: DENSINDIV_MIN=1e-9 !
  ! REAL(dp),PARAMETER:: Kbiometric=50.0 ! Constant in height-diameter relationship
  ! REAL(dp),PARAMETER:: WD= 300.0 ! Wood density kgC/m3
  ! ! threshold growth efficiency for enhanced mortality (higher value gives higher biomass turnover)
  ! REAL(dp),PARAMETER:: GROWTH_EFFICIENCY_MIN=0.008 ! 0.008
  ! REAL(dp),PARAMETER:: Pmort=5.0 ! exponent in mortality formula
  ! REAL(dp),PARAMETER:: MORT_MAX=0.3 ! upper asymptote for enhanced mortality
  ! REAL(dp),PARAMETER:: THETA_recruit=0.95 ! shape parameter in recruitment equation
  ! REAL(dp),PARAMETER:: CMASS_STEM_INIT= 1e-4 ! initial biomass kgC/m2
  ! REAL(dp),PARAMETER:: POWERbiomass=0.75 ! exponent for biomass in proportion to which cohorts preempt resources
  ! REAL(dp),PARAMETER:: POWERGrowthEfficiency = 0.75
  ! REAL(dp),PARAMETER:: CrowdingFactor = 0.0128
  ! REAL(dp),PARAMETER:: ALPHA_CPC = 3.0

  real(dp), parameter :: FULTON_ALPHA = 3.5_dp ! recruitment scalar alpha in Fulton (1991)
  real(dp), parameter :: DENSINDIV_MAX = 0.2_dp  ! 0.5 !  Maximum density of individuals within a cohort indiv/m2
  real(dp), parameter :: DENSINDIV_MIN = 1.0e-9_dp !
  real(dp), parameter :: Kbiometric = 50.0_dp ! Constant in height-diameter relationship
  real(dp), parameter :: WD = 300.0_dp ! Wood density kgC/m3
  ! threshold growth efficiency for enhanced mortality (higher value gives higher biomass turnover)
  real(dp), parameter :: GROWTH_EFFICIENCY_MIN = 0.009_dp ! 0.0095 ! 0.0089 ! 0.0084
  real(dp), parameter :: Pmort = 5.0_dp ! exponent in mortality formula
  real(dp), parameter :: MORT_MAX = 0.3_dp ! upper asymptote for enhanced mortality
  real(dp), parameter :: THETA_recruit = 0.95_dp ! shape parameter in recruitment equation
  real(dp), parameter :: CMASS_STEM_INIT = 1.0e-4_dp ! initial biomass kgC/m2
  real(dp), parameter :: POWERbiomass = 0.67_dp ! exponent for biomass in proportion to which cohorts preempt resources
  real(dp), parameter :: POWERGrowthEfficiency = 0.67_dp
  real(dp), parameter :: CrowdingFactor = 0.043_dp ! 0.043 ! 0.039  !0.029 ! 0.033
  real(dp), parameter :: ALPHA_CPC = 3.5_dp
  real(dp), parameter :: k_allom1 = 200.0_dp ! crown area =  k_allom1 * diam ** k_rp
  real(dp), parameter :: k_rp = 1.67_dp  ! constant in crown area relation to tree diameter
  real(dp), parameter :: ksapwood = 0.05_dp ! rate constant for conversion of sapwood to heartwood (y-1)
  real(dp), parameter :: Q=7.0_dp ! governs rate of increase of mortality with age (2=exponential)
  real,     parameter :: rshootfrac = 0.63
  real(dp), parameter :: shootfrac = real(rshootfrac,dp)
  real(dp), parameter :: CtoNw = 400.0_dp
  real(dp), parameter ::  CtoNl = 60.0_dp
  real(dp), parameter :: CtoNr = 70.0_dp
  real(dp), parameter :: N_EXTENT = 2.0_dp ! multiple of crown diameters within which tree competes with other cohorts
  real(dp), parameter :: EPS = 1.0e-12_dp
  integer(i4b), parameter :: NLAYER = 1 ! number of vertical veg layers (1 is currently the only option)
  integer(i4b), parameter :: NCOHORT_MAX = 20 ! maximum number of cohorts
  integer(i4b), parameter :: NDISTURB = 1 ! number of disturbance regimes (1 (total only)  or 2 (partial and total))
  integer(i4b), parameter :: PATCH_REPS = 10 ! higher number reduces 'noise'
  integer(i4b), parameter :: NAGE_MAX = 1 ! number of maxium ages
  integer(i4b), parameter :: PATCH_REPS1 = 60 ! number of first dist years
  integer(i4b), parameter :: PATCH_REPS2 = 1 ! number of second dist years
  integer(i4b), parameter :: NPATCH = PATCH_REPS1*PATCH_REPS2
  integer(i4b), parameter :: NPATCH1D = NPATCH
  integer(i4b), parameter :: NPATCH2D = NPATCH
  integer(i4b), parameter ::  HEIGHT_BINS = 12 ! number of height categories to keep track of for diagnostics
  real(dp), parameter :: BIN_POWER = 1.4_dp ! bins have muscles
  ! Time base factor (to be multiplied by mean dist interval to give TIMEBASE)
  ! for sampling disturbance probabilities from Poisson distribution
  integer(i4b), parameter :: TIMEBASE_FACTOR=50
  real(dp), parameter :: PI=3.14159265358979323846264_dp
  ! 0 == default; 1 = top-end allometry (requires precip as input to POPSTEP); 2 = Allometry following Williams 2005, Model 5b
  integer(i4b), parameter :: ALLOM_SWITCH = 2
  ! 0 == binnned max height variable; 1 = continuous (needs lots of memory); 2 = binned by integer heights
  integer(i4b), parameter :: MAX_HEIGHT_SWITCH = 2
  integer(i4b), parameter :: RESOURCE_SWITCH = 1 ! 0 = default; 1  fraction net resource uptake
  integer(i4b), parameter :: RECRUIT_SWITCH = 1 ! 0 = default, 1 = Pgap-dependence
  integer(i4b), parameter :: INTERP_SWITCH = 1 ! 0 = sum over weighted patches, 1 = sum over interpolated patches
  integer(i4b), parameter :: SMOOTH_SWITCH = 0 ! smooth disturbance flux
  integer(i4b), parameter :: NYEAR_WINDOW  = 5                  ! one-side of smoothing window (y)
  integer(i4b), parameter :: NYEAR_SMOOTH  = 2*NYEAR_WINDOW + 1 ! smoothing window (y)
  integer(i4b), parameter :: NYEAR_HISTORY = NYEAR_SMOOTH-NYEAR_WINDOW
  integer(i4b), parameter :: AGEMAX = 1000

end module POP_Constants


!*******************************************************************************


module POP_Types

  use TYPEdef, only: dp, i4b
  use POP_Constants, only: NCOHORT_MAX, NLAYER, HEIGHT_BINS, NDISTURB, NPATCH, NPATCH2D, &
       NYEAR_HISTORY, AGEMAX

  implicit none

  type Cohort
     integer(i4b) :: id
     integer(i4b) :: age ! cohort age
     real(dp)     :: biomass ! cohort biomass
     real(dp)     :: density ! landscape tree density (weighted mean over patches)
     real(dp)     :: frac_resource_uptake
     real(dp)     :: frac_light_uptake
     real(dp)     :: frac_interception
     real(dp)     :: frac_respiration
     real(dp)     :: frac_NPP
     real(dp)     :: respiration_scalar
     real(dp)     :: crown_area
     real(dp)     :: Pgap
     real(dp)     :: height
     real(dp)     :: diameter
     real(dp)     :: sapwood
     real(dp)     :: heartwood
     real(dp)     :: sapwood_area
     real(dp)     :: basal_area
     real(dp)     :: LAI
     real(dp)     :: Cleaf
     real(dp)     :: Croot
  end type Cohort

  type Layer
     type(Cohort), dimension(NCOHORT_MAX) :: Cohort
     integer(i4b) :: ncohort ! number of cohorts with density >0
     real(dp)     :: biomass ! layer biomass
     real(dp)     :: density ! layer tree density
     real(dp)     :: hmean ! layer mean tree height (weighted mean over patches)
     real(dp)     :: hmax  ! layer max tree height
  end type Layer

  type Patch
     type(Layer), dimension(NLAYER) :: Layer
     real(dp)     :: factor_recruit
     real(dp)     :: pgap
     real(dp)     :: lai
     real(dp)     :: biomass ! total biomass in patch
     real(dp)     :: biomass_old ! total biomass in patch
     real(dp)     :: sapwood ! total sapwood biomass in patch
     real(dp)     :: heartwood ! total heartwood biomass in patch
     real(dp)     :: sapwood_old ! total sapwood biomass in patch
     real(dp)     :: sapwood_area  ! total sapwood area in patch
     real(dp)     :: sapwood_area_old ! total sapwood area in patch
     real(dp)     :: stress_mortality ! biomass lost in each patch due to stress
     real(dp)     :: fire_mortality ! biomass lost in each patch due partial fire disturbance
     real(dp)     :: cat_mortality ! biomass lost in each patch due partial fire disturbance
     real(dp)     :: crowding_mortality ! biomass lost to crowding mortality
     real(dp)     :: cpc
     real(dp)     :: mortality !
     real(dp)     :: sapwood_loss
     real(dp)     :: sapwood_area_loss
     real(dp)     :: growth ! biomass growth in each patch due to stem increment
     real(dp)     :: area_growth ! basal area growth in each patch due to stem increment
     integer(i4b) :: disturbance_interval(NDISTURB)  ! prescribed disturbance(s) interval for this patch
     integer(i4b) :: first_disturbance_year(NDISTURB)
     integer(i4b) :: age(NDISTURB) ! number of years since last disturbance(s)
     integer(i4b) :: id
     real(dp)     :: frac_NPP
     real(dp)     :: frac_respiration
     real(dp)     :: frac_light_uptake
     real(dp)     :: fire_top_kill_density
  end type Patch

  type Landscape
     type(Patch), dimension(NPATCH2D) :: patch
     real(dp), dimension(NPATCH2D)    :: freq ! patch weighting
     real(dp), dimension(NPATCH2D)    :: freq_old ! patch weighting (previous time-step)
     real(dp), dimension(NPATCH2D)    :: fire_freq     !
     real(dp), dimension(NPATCH2D)    :: fire_freq_old !
     real(dp), dimension(NPATCH2D)    :: cat_freq      !
     real(dp), dimension(NPATCH2D)    :: cat_freq_old  !
     real(dp), dimension(NPATCH2D,NDISTURB) :: freq_ranked_age_unique ! unique age weighting
     integer(i4b), dimension(NPATCH2D, NDISTURB) :: ranked_age_unique ! unique age
     integer(i4b), dimension(NDISTURB) :: n_age ! number of unique ages
     real(dp), dimension(NLAYER) :: biomass ! landscape stem biomass (weighted mean over patches)
     real(dp), dimension(NLAYER) :: density ! landscape tree density (weighted mean over patches)
     real(dp), dimension(NLAYER) :: hmean ! landscape mean treen height (weighted mean over patches)
     real(dp), dimension(NLAYER) :: hmax  ! landscape max tree height
     real(dp), dimension(HEIGHT_BINS) :: cmass_stem_bin ! biomass by height bin
     real(dp), dimension(HEIGHT_BINS) :: densindiv_bin ! density by height bin
     real(dp), dimension(HEIGHT_BINS) :: height_bin ! mean height in each bin
     real(dp), dimension(HEIGHT_BINS) :: diameter_bin ! mean diameter in each bin
     character(100), dimension(HEIGHT_BINS) :: bin_labels ! text strings for bin bounds
     real(dp) :: cmass_sum ! landscape biomass
     real(dp) :: cmass_sum_old ! landscape biomass
     real(dp) :: cheartwood_sum ! landscape biomass (heart wood)
     real(dp) :: csapwood_sum ! landscape biomass (sap wood)
     real(dp) :: csapwood_sum_old ! landscape biomass
     real(dp) :: densindiv ! landscape density of individuals
     real(dp) :: height_mean
     real(dp) :: height_max
     real(dp) :: basal_area
     real(dp) :: sapwood_loss ! (kg C m-2 y-1) ! total sapwood loss (turnover + mortality)
     real(dp) :: sapwood_area_loss ! ( m2/m-2 y-1) sapwood area loss (mortality only)
     real(dp) :: stress_mortality ! (kg C m-2 y-1)
     real(dp) :: crowding_mortality ! (kg C m-2 y-1)
     real(dp) :: fire_mortality ! (kg C m-2 y-1)
     real(dp) :: cat_mortality ! (kg C m-2 y-1)
     real(dp) :: res_mortality ! (kg C m-2 y-1)
     real(dp) :: growth
     real(dp) :: area_growth ! m2/ha
     real(dp) :: crown_cover
     real(dp) :: crown_area
     real(dp) :: crown_volume
     real(dp) :: sapwood_area
     real(dp) :: sapwood_area_old
     real(dp) :: Kclump ! clumping factor
     integer(i4b) :: npatch_active
     integer(i4b) :: LU
     real(dp) :: smoothing_buffer
     real(dp) :: smoothing_buffer_cat
     real(dp) :: fire_mortality_smoothed
     real(dp) :: cat_mortality_smoothed
     real(dp), dimension(NYEAR_HISTORY) :: fire_mortality_history
     real(dp), dimension(NYEAR_HISTORY) :: cat_mortality_history
     real(dp), dimension(AGEMAX) :: freq_age ! age weighting (by age in y: 0:AGE_MAX-1)
     real(dp), dimension(AGEMAX) :: biomass_age
     real(dp) :: rkill
  end type Landscape

  type POP_TYPE
     type(Landscape), dimension(:), allocatable :: pop_grid
     integer,         dimension(:), allocatable :: it_pop
     integer                                    :: np
     integer,         dimension(:), allocatable :: Iwood ! , LU
  end type POP_TYPE

end module POP_Types


!*******************************************************************************


module POPModule
  !-------------------------------------------------------------------------------
  ! * This module contains all subroutines for POP calcs at a single time step.
  !-------------------------------------------------------------------------------
  use TYPEdef, only: sp, i4b
  use POP_Types
  use POP_Constants

  implicit none

contains

  !*******************************************************************************

  subroutine ZeroPOP(POP,n)

#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif

    implicit none

    type(POP_TYPE),    intent(INOUT) :: POP
    integer, optional, intent(IN)    :: n

    integer:: g, k, l, c, np, a, b
#ifdef __MPI__
    integer :: ierr
#endif

    if (.not. allocated(pop%pop_grid)) then
       write(*,*)" POP not allocated! Abort in ZeroPOP."
#ifdef __MPI__
       call MPI_Abort(0, 84, ierr) ! Do not know comm nor rank here
#else
       stop 84
#endif
    endif

    np = size(pop%pop_grid)

    ! optional integer n intended for zeroing secondary forest tiles
    if (present(n)) then
       a = n
       b = n
       !pop%LU(n) = 2
       POP%pop_grid(n)%LU = 2
    else
       a = 1
       b = np
       !pop%LU = 1
       POP%pop_grid%LU = 1
    endif

    do g=a, b
       POP%pop_grid(g)%freq                    = 0.0_dp ! patch weighting
       POP%pop_grid(g)%freq_old                = 0.0_dp ! patch weighting
       POP%pop_grid(g)%fire_freq               = 0.0_dp
       POP%pop_grid(g)%fire_freq_old           = 0.0_dp
       POP%pop_grid(g)%cat_freq                = 0.0_dp
       POP%pop_grid(g)%cat_freq_old            = 0.0_dp
       POP%pop_grid(g)%freq_ranked_age_unique  = 0.0_dp
       POP%pop_grid(g)%ranked_age_unique       = 0
       POP%pop_grid(g)%n_age                   = 0
       POP%pop_grid(g)%biomass                 = 0.0_dp ! landscape stem biomass (weighted mean over patches)
       POP%pop_grid(g)%density                 = 0.0_dp ! landscape tree density (weighted mean over patches)
       POP%pop_grid(g)%hmean                   = 0.0_dp ! landscape mean treen height (weighted mean over patches)
       POP%pop_grid(g)%hmax                    = 0.0_dp ! landscape max tree height
       POP%pop_grid(g)%cmass_stem_bin          = 0.0_dp ! biomass by height bin
       POP%pop_grid(g)%densindiv_bin           = 0.0_dp ! density by height bin
       POP%pop_grid(g)%height_bin              = 0.0_dp ! mean height in each bin
       POP%pop_grid(g)%diameter_bin            = 0.0_dp ! mean diameter in each bin
       POP%pop_grid(g)%bin_labels              = ' '    ! text strings for bin bounds
       POP%pop_grid(g)%cmass_sum               = 0.0_dp ! landscape biomass
       POP%pop_grid(g)%cmass_sum_old           = 0.0_dp ! landscape biomass
       POP%pop_grid(g)%cheartwood_sum          = 0.0_dp ! landscape biomass
       POP%pop_grid(g)%csapwood_sum            = 0.0_dp ! landscape biomass
       POP%pop_grid(g)%csapwood_sum_old        = 0.0_dp ! landscape biomass
       POP%pop_grid(g)%densindiv               = 0.0_dp ! landscape density of individuals
       POP%pop_grid(g)%height_mean             = 0.0_dp
       POP%pop_grid(g)%height_max              = 0.0_dp
       POP%pop_grid(g)%basal_area              = 0.0_dp
       POP%pop_grid(g)%sapwood_loss            = 0.0_dp
       POP%pop_grid(g)%sapwood_area_loss       = 0.0_dp
       POP%pop_grid(g)%stress_mortality        = 0.0_dp ! (kg C m-2 y-1)
       POP%pop_grid(g)%crowding_mortality      = 0.0_dp ! (kg C m-2 y-1)
       POP%pop_grid(g)%fire_mortality          = 0.0_dp ! (kg C m-2 y-1)
       POP%pop_grid(g)%cat_mortality           = 0.0_dp ! (kg C m-2 y-1)
       POP%pop_grid(g)%res_mortality           = 0.0_dp ! (kg C m-2 y-1)
       POP%pop_grid(g)%growth                  = 0.0_dp
       POP%pop_grid(g)%area_growth             = 0.0_dp
       POP%pop_grid(g)%crown_cover             = 0.0_dp
       POP%pop_grid(g)%crown_area              = 0.0_dp
       POP%pop_grid(g)%crown_volume            = 0.0_dp
       POP%pop_grid(g)%sapwood_area            = 0.0_dp
       POP%pop_grid(g)%sapwood_area_old        = 0.0_dp
       POP%pop_grid(g)%Kclump                  = 1.0_dp
       POP%pop_grid(g)%npatch_active           = 0
       POP%pop_grid(g)%smoothing_buffer        = 0.0_dp
       POP%pop_grid(g)%smoothing_buffer_cat    = 0.0_dp
       POP%pop_grid(g)%fire_mortality_smoothed = 0.0_dp
       POP%pop_grid(g)%cat_mortality_smoothed  = 0.0_dp
       POP%pop_grid(g)%fire_mortality_history  = 0.0_dp
       POP%pop_grid(g)%cat_mortality_history   = 0.0_dp
       POP%pop_grid(g)%freq_age                = 0.0_dp
       POP%pop_grid(g)%rkill                   = 0.0_dp
       if (present(n)) then
          POP%pop_grid(g)%freq_age(1) = 1.0_dp
       endif
       POP%pop_grid(g)%biomass_age = 0.0_dp

       do k=1, NPATCH2D
          POP%pop_grid(g)%patch(k)%factor_recruit         = 0.0_dp
          POP%pop_grid(g)%patch(k)%pgap                   = 0.0_dp
          POP%pop_grid(g)%patch(k)%lai                    = 0.0_dp
          POP%pop_grid(g)%patch(k)%biomass                = 0.0_dp ! total biomass in patch
          POP%pop_grid(g)%patch(k)%biomass_old            = 0.0_dp
          POP%pop_grid(g)%patch(k)%sapwood                = 0.0_dp ! total biomass in patch (sapwood)
          POP%pop_grid(g)%patch(k)%heartwood              = 0.0_dp ! total biomass in patch (heartwood)
          POP%pop_grid(g)%patch(k)%sapwood_old            = 0.0_dp
          POP%pop_grid(g)%patch(k)%sapwood_area           = 0.0_dp
          POP%pop_grid(g)%patch(k)%sapwood_area_old       = 0.0_dp
          POP%pop_grid(g)%patch(k)%stress_mortality       = 0.0_dp ! biomass lost in each patch due to stress
          POP%pop_grid(g)%patch(k)%fire_mortality         = 0.0_dp ! biomass lost in each patch due to fire partial dist
          POP%pop_grid(g)%patch(k)%cat_mortality          = 0.0_dp ! biomass lost in each patch due to fire partial dist
          POP%pop_grid(g)%patch(k)%crowding_mortality     = 0.0_dp
          POP%pop_grid(g)%patch(k)%cpc                    = 0.0_dp
          POP%pop_grid(g)%patch(k)%mortality              = 0.0_dp
          POP%pop_grid(g)%patch(k)%sapwood_loss           = 0.0_dp
          POP%pop_grid(g)%patch(k)%sapwood_area_loss      = 0.0_dp
          POP%pop_grid(g)%patch(k)%growth                 = 0.0_dp ! biomass growth in each patch due stem increment
          POP%pop_grid(g)%patch(k)%area_growth            = 0.0_dp
          POP%pop_grid(g)%patch(k)%disturbance_interval   = 0      ! prescribed disturbance(s) interval for this patch
          POP%pop_grid(g)%patch(k)%first_disturbance_year = 0
          POP%pop_grid(g)%patch(k)%age                    = 0      ! number of years since last disturbance(s)
          POP%pop_grid(g)%patch(k)%id                     = 0
          POP%pop_grid(g)%patch(k)%frac_NPP               = 0.0_dp
          POP%pop_grid(g)%patch(k)%frac_respiration       = 0.0_dp
          POP%pop_grid(g)%patch(k)%frac_light_uptake      = 0.0_dp
          POP%pop_grid(g)%patch(k)%fire_top_kill_density  = 0.0_dp

          do l=1, NLAYER
             POP%pop_grid(g)%patch(k)%Layer(L)%ncohort = 0      ! number of cohorts with density >0.0_dp
             POP%pop_grid(g)%patch(k)%Layer(L)%biomass = 0.0_dp ! layer biomass
             POP%pop_grid(g)%patch(k)%Layer(L)%density = 0.0_dp ! layer tree density
             POP%pop_grid(g)%patch(k)%Layer(L)%hmean   = 0.0_dp ! layer mean tree height (weighted mean over patches)
             POP%pop_grid(g)%patch(k)%Layer(L)%hmax    = 0.0_dp ! layer max tree height

             do c=1, NCOHORT_MAX
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%id                   = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%age                  = 0      ! cohort age
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%biomass              = 0.0_dp ! cohort biomass
                ! landscape tree density (weighted mean over patches)
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%density              = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_resource_uptake = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_light_uptake    = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_interception    = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_respiration     = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_NPP             = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%respiration_scalar   = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%crown_area           = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%Pgap                 = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%height               = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%diameter             = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%sapwood              = 0.0_dp ! cohort sapwood
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%heartwood            = 0.0_dp ! cohort heartwood
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%sapwood_area         = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%basal_area           = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%LAI                  = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%Cleaf                = 0.0_dp
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%Croot                = 0.0_dp
             enddo ! NCOHORT_MAX

          enddo ! NLAYER

       enddo ! NPATCH2D

    enddo ! pop_grid%np

  end subroutine ZeroPOP


  !*******************************************************************************


  subroutine InitPOP2D_Poisson(POP, mean_disturbance_interval, m)
    ! Initialises vector of patches with maximum age correpondding to 95% of pdf
    ! Starting year: uniform distribution up to maximum age

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer(i4b),   intent(IN)    :: mean_disturbance_interval(:,:)
    integer(i4b),   intent(IN), optional :: m

    integer(i4b) :: j, k, g, ipatch, idist, p, c, i
    integer(i4b) :: disturbance_interval
    integer(i4b):: Poisson_age(1000)
    real(dp):: Poisson_weight(1000), CumPoisson_weight(1000)
    integer:: i_max
    integer:: np
    real(dp):: disturbance_freq
    integer:: tmp2(PATCH_REPS1), tmp3(PATCH_REPS2)
    integer:: a,b

    np = size(POP%pop_grid)
    a = 1
    b = np
    if (present(m)) then
       a = m
       b = m
    endif

    do g=a,b
       ! calculate Poisson weights for each of the 2 mean disturbance intervals
       if (NPATCH.gt.1) then
          do idist=1,NDISTURB
             disturbance_freq=1.0_dp/real(mean_disturbance_interval(g,idist),dp)
             do p = 1,1000
                Poisson_age(p) = p
                Poisson_weight(p) = Exponential(disturbance_freq,p)
                CumPoisson_weight(p) = CumExponential(disturbance_freq,real(p,dp))
             enddo
             ! set max age to correspond to 95% percentile of cum pdf
             do k =1,NPATCH2D
                i_max = maxloc(Poisson_age,1,CumPoisson_weight.le.0.95_dp)
                POP%pop_grid(g)%patch(k)%disturbance_interval(idist) = Poisson_age(i_max)
                POP%pop_grid(g)%patch(k)%id = k
                POP%pop_grid(g)%patch(k)%age = 0
             enddo
          enddo

          do idist =1,ndisturb
             ! set first disturbance year for first dist interval class
             if (idist .eq. 1) then
                disturbance_interval = POP%pop_grid(g)%patch(1)%disturbance_interval(idist)
                do c = 1,PATCH_REPS1
                   if (c==1) then
                      tmp2(1) = 1
                   else
                      tmp2(1) = max(disturbance_interval*(c-1)/(PATCH_REPS1),1)+1
                   endif
                   tmp2(2) = max(disturbance_interval*c/(PATCH_REPS1),1)
                   tmp2(3) = tmp2(1)
                   !    write(*,*) 'tmp2', c, disturbance_interval, tmp2(1),tmp2(2)
                   do j = 1,PATCH_REPS2
                      ipatch = (c-1)*PATCH_REPS2 + j
                      POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(idist) = tmp2(3)
                      tmp2(3)=tmp2(3)+1
                      if (tmp2(3)>tmp2(2)) then
                         tmp2(3) = tmp2(1)
                      endif
                   enddo
                enddo

                ! ! set first disturbance year for first dist interval class
                ! idist = 1
                ! disturbance_interval = POP%pop_grid(g)%patch(1)%disturbance_interval(idist)
                ! DO c = 1,PATCH_REPS1
                !    tmp2(c) = max(disturbance_interval*c/(PATCH_REPS1),1)
                ! ENDDO
                ! DO c = 1,PATCH_REPS1
                !    i = 0
                !    DO j = 1,PATCH_REPS2
                !       ipatch = (j-1)*PATCH_REPS1 + c
                !       i = i+1
                !       IF (i.gt.PATCH_REPS1) then
                !          i = 1
                !       ENDIF
                !       do while ((tmp2(i+1).eq. tmp2(i)).and.(i.lt.PATCH_REPS1))
                !          i = i+1
                !          IF (i.gt.PATCH_REPS1) then
                !             i = 1
                !          ENDIF

                !       ENDDO

                !       write(*,*) i, tmp2(i)
                !       POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(idist) = tmp2(i)

                !    ENDDO
                ! ENDDO

                ! set first disturbance year for first 2nd interval class
             elseif (idist.eq.2) then
                disturbance_interval = POP%pop_grid(g)%patch(1)%disturbance_interval(idist)
                do c = 1,PATCH_REPS2
                   tmp3(c) = max(disturbance_interval*(c-1)/(PATCH_REPS2),1)
                enddo

                do c = 1,PATCH_REPS2
                   i = 0
                   do j = 1,PATCH_REPS1
                      ipatch = (j-1)*PATCH_REPS2 + c
                      POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(idist) = &
                           tmp3(c) +(j-1)*max((tmp3(idist)-tmp3(1))/PATCH_REPS1,1)

                      !  i = i+1
                      !  if (i.gt.(tmp3(2)-tmp3(1))) i = 0
                   enddo
                enddo
             endif

          enddo

       else   ! NPATCH =1 (single patch mode)
          k = 1
          do idist=1,NDISTURB
             POP%pop_grid(g)%patch(k)%disturbance_interval(idist) = mean_disturbance_interval(g,idist)
             POP%pop_grid(g)%patch(k)%first_disturbance_year(idist) = 113
             POP%pop_grid(g)%patch(k)%age = 0
             POP%pop_grid(g)%patch(k)%id = k
          enddo
       endif

       POP%pop_grid(g)%npatch_active = NPATCH

    enddo

  end subroutine InitPOP2D_Poisson


  !*******************************************************************************


  subroutine POPStep(POP, StemNPP, disturbance_interval, disturbance_intensity,LAI,Cleaf,Croot, &
       NPPtoGPP, StemNPP_av,frac_intensity1,precip)

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    real(dp), intent(IN) :: StemNPP(:,:)
    real(dp), intent(IN) :: disturbance_intensity(:,:)
    integer(i4b), intent(IN) ::  disturbance_interval(:,:)
    real(dp), intent(IN) ::  LAI(:)
    real(dp), intent(IN) ::  Cleaf(:)
    real(dp), intent(IN) ::  Croot(:)
    real(dp), intent(IN) ::  NPPtoGPP(:)
    real(dp), intent(IN), optional :: frac_intensity1(:), precip(:)
    real(dp), intent(IN), optional :: StemNPP_av(:)

    integer(i4b) :: idisturb,np,g
    integer(i4b), allocatable :: it(:)

    !INTEGER, INTENT(IN) :: wlogn
    pop%it_pop = pop%it_pop + 1
    !it = pop%it_pop(1)
    np = size(POP%POP_grid)
    allocate(it(np))

    do g=1, np
       it(g) = maxval(pop%pop_grid(g)%patch(:)%age(1)) + 1
    enddo
    ! DO idisturb = 1,NDISTURB
    !    CALL GetUniqueAgeFrequencies(POP, disturbance_interval, idisturb)
    ! ENDDO

    ! CALL GetPatchFrequencies(POP)

    !call flush(wlogn)
    if (present(precip)) then
       if(present(StemNPP_av)) then
          call PatchAnnualDynamics(POP, StemNPP, NPPtoGPP, it, precip=precip, StemNPP_av=StemNPP_av)
       else
          call PatchAnnualDynamics(POP, StemNPP, NPPtoGPP, it, precip=precip)
       endif
    else
       if(present(StemNPP_av)) then
          call PatchAnnualDynamics(POP, StemNPP, NPPtoGPP, it, StemNPP_av=StemNPP_av)
       else
          call PatchAnnualDynamics(POP, StemNPP, NPPtoGPP, it)
       endif
    endif

    if (NDISTURB.eq.1) then
       if (present(precip)) then
          !   CALL Patch_disturb(POP,it,1,precip)
          call Patch_partial_disturb2(POP,1)
       else
          call Patch_disturb(POP,1)
          ! CALL Patch_partial_disturb2(POP,it)
       endif
    elseif (NDISTURB.eq.2) then
       if (present(frac_intensity1)) then
          call Patch_partial_disturb(POP,1,disturbance_intensity,frac_intensity1=frac_intensity1)
       else
          call Patch_partial_disturb(POP,1,disturbance_intensity)
       endif
       if (present(precip)) then
          !CALL Patch_partial_disturb2(POP,it,2)
          call Patch_disturb(POP,2,precip)
       else
          ! CALL Patch_partial_disturb2(POP,it,2)
          call Patch_disturb(POP,2)
       endif
    endif

    do idisturb = 1,NDISTURB
       call GetUniqueAgeFrequencies(POP, disturbance_interval, idisturb)
    enddo

    call GetPatchFrequencies(POP)

    if (present(precip)) then
       call GetDiagnostics(pop, LAI,Cleaf,Croot, disturbance_interval, it,precip)
    else
       call GetDiagnostics(pop, LAI,Cleaf,Croot, disturbance_interval, it)
    endif

  end subroutine POPStep


  !*******************************************************************************


  subroutine PatchAnnualDynamics(pop, StemNPP, NPPtoGPP, it, StemNPP_av, precip)

    implicit none

    type( POP_TYPE ), intent(INOUT) :: pop
    real(dp), intent(IN)            :: StemNPP(:,:)
    real(dp), intent(IN)            :: NPPtoGPP(:)
    real(dp), intent(IN), optional  :: precip(:)
    real(dp), optional, intent(IN)            :: StemNPP_av(:)
    integer(i4b), intent(IN)        :: it(:)

    real(dp) :: densindiv
    real(dp) :: tmp,tmp_light,tmp_respiration,tmp_fracnpp, cmass_stem_inc
    integer(i4b) :: j, k,c, idist
    integer(i4b) :: ivec(NCOHORT_MAX), nc, np
    real(dp) :: growth_efficiency,cmass_stem
    real(dp) :: mort
    real(dp) :: cpc, crown_area
    real(dp) :: mort_cpc
    real(dp) :: ht, diam, area_growth_grid , basal_grid, basal_new, basal_old
    real(dp) :: tmp2(NCOHORT_MAX), freq

    np = size(POP%POP_grid)

    ! growth
    ! Distributes layer biomass increment among cohorts and increments age
    ! calculate fractional resource uptake by each cohort
    do j=1,np
       basal_grid = 0.0_dp
       area_growth_grid = 0.0_dp
       pop%pop_grid(j)%sapwood_area_old      = pop%pop_grid(j)%sapwood_area
       pop%pop_grid(j)%freq_old      = pop%pop_grid(j)%freq
       pop%pop_grid(j)%fire_freq_old = pop%pop_grid(j)%fire_freq
       pop%pop_grid(j)%cat_freq_old  = pop%pop_grid(j)%cat_freq


       ! Get fraction allocation for each patch
       tmp = 0.0_dp
       tmp_light = 0.0_dp
       tmp_respiration = 0.0_dp
       tmp_fracNPP = 0.0_dp

       if (NPATCH2D >1 .and. it(j) > 1 .and. RESOURCE_SWITCH>0) then
          do k=1,NPATCH2D
             freq =  pop%pop_grid(j)%freq(pop%pop_grid(j)%patch(k)%id)
             nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             do c=1,nc
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_light_uptake = &
                     pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%frac_interception      ! defined in terms of Pgap
                ! total autotrophic resp, summed over all cohorts and patches
                tmp_respiration = tmp_respiration + &
                     freq*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%respiration_scalar
             enddo
             tmp_light = tmp_light + freq*(1.0_dp - pop%pop_grid(j)%patch(k)%Pgap)
          enddo
          if (tmp_respiration .gt. 1.0e-8_dp .and. tmp_light .gt. 1.0e-8_dp) then
             do k=1,NPATCH2D
                ! fraction respiration and un-normalised NPP for each patch
                nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
                freq =  pop%pop_grid(j)%freq(pop%pop_grid(j)%patch(k)%id)
                ! frac autotrophic resp

                pop%pop_grid(j)%patch(k)%frac_respiration = &
                     sum(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%respiration_scalar)/tmp_respiration

                ! frac gpp
                pop%pop_grid(j)%patch(k)%frac_light_uptake = &
                     (1.0_dp - pop%pop_grid(j)%patch(k)%pgap) /tmp_light
                ! frac npp
                pop%pop_grid(j)%patch(k)%frac_NPP = &
                     max(pop%pop_grid(j)%patch(k)%frac_light_uptake*(1.0_dp/NPPtoGPP(j)) - &
                     pop%pop_grid(j)%patch(k)%frac_respiration*(1.0_dp/NPPtoGPP(j)-1.0_dp),0.0_dp)

                tmp_fracNPP = tmp_fracNPP + freq*pop%pop_grid(j)%patch(k)%frac_NPP

             enddo

             ! normalised fraction NPP
             do k=1,NPATCH2D
                pop%pop_grid(j)%patch(k)%frac_NPP = &
                     pop%pop_grid(j)%patch(k)%frac_NPP/tmp_fracNPP

             enddo
          else
             pop%pop_grid(j)%patch(:)%frac_NPP = 1.0_dp
             pop%pop_grid(j)%patch(:)%frac_respiration = 1.0_dp
             pop%pop_grid(j)%patch(:)%frac_light_uptake = 1.0_dp


          endif
       else
          pop%pop_grid(j)%patch(:)%frac_NPP = 1.0_dp
          pop%pop_grid(j)%patch(:)%frac_respiration = 1.0_dp
          pop%pop_grid(j)%patch(:)%frac_light_uptake = 1.0_dp
       endif
       ! End Get fraction allocation for each patch
       ! Get fraction allocation for each cohort in each patch
       do k=1,NPATCH2D
          tmp = 0.0_dp
          tmp_light = 0.0_dp
          tmp_respiration = 0.0_dp
          tmp_fracNPP = 0.0_dp
          if (pop%pop_grid(j)%patch(k)%Layer(1)%ncohort>1) then
             nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             do c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

                cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
                densindiv = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density

                ! get initial basal area

                if ( present(precip) ) then
                   call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_old, precip(j))
                else
                   call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_old )
                endif

                if (ALLOM_SWITCH.eq.1) then
                   !! assumes crown radius (m) = 0.1492 * dbh (cm) (from G. Cook, pers. comm.)
                   crown_area = densindiv*PI*(diam*100.0_dp*0.1492_dp)**2
                else
                   crown_area = densindiv*PI*(((k_allom1 * diam ** k_rp )/PI)**0.5_dp)**2
                endif

                tmp = tmp + (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass/ &  ! sum over all cohorts
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density)**POWERbiomass * &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density

                tmp_light = tmp_light + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception

                tmp_respiration = tmp_respiration + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%respiration_scalar

                tmp2(c) = sum((pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c:nc)%biomass/ &  ! sum over all cohorts c:nc
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c:nc)%density)**POWERbiomass * &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c:nc)%density)

             enddo

             ! un-normalised fractional gross resource uptake: weighted combination of components
             ! where cohort competes with older cohorts and where it does not
             do c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

                if (RESOURCE_SWITCH ==1) then
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception = &
                        pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception/tmp_light
                else
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception = 1.0_dp
                endif

             enddo

             !normalised fractional gross resource uptake
             nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             do c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
                !normalised fractional gross resource uptake
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_light_uptake = &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_interception/ &
                     sum(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%frac_interception)
             enddo


             ! fraction respiration and un-normalised NPP
             do c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_respiration = &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%respiration_scalar/tmp_respiration


                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP = &
                     max(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_light_uptake*(1.0_dp/NPPtoGPP(j)) - &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_respiration*(1.0_dp/NPPtoGPP(j)-1.0_dp),0.0_dp)


                tmp_fracNPP = tmp_fracNPP +  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP

             enddo

             ! normalised fraction NPP
             do c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP = &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_NPP/tmp_fracNPP

             enddo


             ! fraction net resource uptake

             do c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

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

             enddo

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


       enddo

       tmp = 0
       do k=1,NPATCH2D
          pop%pop_grid(j)%patch(k)%sapwood_loss = 0.0_dp
          pop%pop_grid(j)%patch(k)%sapwood_area_loss = 0.0_dp
          pop%pop_grid(j)%patch(k)%sapwood_old = pop%pop_grid(j)%patch(k)%sapwood
          pop%pop_grid(j)%patch(k)%sapwood_area_old = pop%pop_grid(j)%patch(k)%sapwood_area
          pop%pop_grid(j)%patch(k)%biomass_old = pop%pop_grid(j)%patch(k)%biomass
          pop%pop_grid(j)%patch(k)%growth = 0.0_dp
          pop%pop_grid(j)%patch(k)%area_growth = 0.0_dp
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          freq =  pop%pop_grid(j)%freq(pop%pop_grid(j)%patch(k)%id)
          do c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

             cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
             densindiv = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density

             ! get initial basal area

             if ( present(precip) ) then
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_old, precip(j))
             else
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_old )
             endif

             ! increment biomass in cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass +  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake*StemNPP(j,1)

             cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
             tmp = tmp + freq*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake

             ! get incremented basal area
             if ( present(precip) ) then
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_new, precip(j))
             else

                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal_new )
             endif


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
                  (1.0_dp - ksapwood)*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood


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

          enddo
          ! Layer biomass (summed over cohorts)
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          pop%pop_grid(j)%patch(k)%Layer(1)%biomass = sum(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%biomass)


       enddo


    enddo

    ! Mortality
    !Implements resource stress mortality and crowding mortality for all cohorts in layer

    do j=1,np
       do k=1,NPATCH2D
          nc = 0
          ivec = 0
          pop%pop_grid(j)%patch(k)%stress_mortality = 0.0_dp
          pop%pop_grid(j)%patch(k)%crowding_mortality = 0.0_dp
          pop%pop_grid(j)%patch(k)%mortality = 0.0_dp
          crown_area = 0.0_dp
          do c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
             cmass_stem_inc=StemNPP(j,1)*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake

             if (present(StemNPP_av)) then
                growth_efficiency=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake* &
                     StemNPP_av(j)  /(cmass_stem**(POWERGrowthEfficiency))
             else
                growth_efficiency=cmass_stem_inc/(cmass_stem**(POWERGrowthEfficiency))
             endif
             ! growth_efficiency=cmass_stem_inc/(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood**(POWERGrowthEfficiency))
             mort=MORT_MAX/(1.0_dp+(growth_efficiency/GROWTH_EFFICIENCY_MIN)**Pmort)

             ! mort = 0 ! test

             pop%pop_grid(j)%patch(k)%stress_mortality = pop%pop_grid(j)%patch(k)%stress_mortality &
                  + mort*cmass_stem
             if (pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter*100_dp .gt. 1.0_dp) then
                if (ALLOM_SWITCH.eq.1) then
                   ! assumes crown radius (m) = 0.14 * dbh (cm)
                   crown_area = crown_area + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density* &
                        PI*(pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter*100.0_dp*0.14_dp)**2
                else
                   crown_area = crown_area + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density* &
                        k_allom1 * pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter ** k_rp
                endif
             else
                crown_area = crown_area + 0.5_dp*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%LAI
             endif

             cpc = 1.0_dp - exp(-crown_area)
             pop%pop_grid(j)%patch(k)%cpc = cpc
             if (cpc.gt.1.0e-3_dp .and.  alpha_cpc * (1.0_dp - 1.0_dp/cpc).gt.-50.0_dp) then
                mort_cpc = exp(alpha_cpc * (1.0_dp - 1.0_dp/cpc))
             else
                mort_cpc = 0.0_dp
             endif

             !mort_cpc = 0 ! test
             pop%pop_grid(j)%patch(k)%crowding_mortality = pop%pop_grid(j)%patch(k)%crowding_mortality + &
                  min((mort_cpc*CrowdingFactor),cmass_stem_inc/cmass_stem)*cmass_stem


             pop%pop_grid(j)%patch(k)%mortality = pop%pop_grid(j)%patch(k)%mortality + &
                  mort*cmass_stem

             pop%pop_grid(j)%patch(k)%sapwood_loss =  pop%pop_grid(j)%patch(k)%sapwood_loss + &
                  mort*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood

             pop%pop_grid(j)%patch(k)%sapwood_area_loss =  pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                  mort*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood_area


             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = cmass_stem*(1.0_dp-mort)
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood *(1.0_dp-mort)
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood *(1.0_dp-mort)

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density*(1.0_dp-mort)
             if (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.lt.DENSINDIV_MIN) then
                ! remove cohort
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 0
                pop%pop_grid(j)%patch(k)%stress_mortality = pop%pop_grid(j)%patch(k)%stress_mortality + &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood_area
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood = 0.0_dp

             else
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 1
                !COMMLN Why is id here 1 instead of c or sth useful? Call it differently
                nc = nc+1
                ivec(nc)=c
             endif
          enddo
          ! SHUFFLE if necessary to remove zero-density cohorts
          if (nc.lt.pop%pop_grid(j)%patch(k)%Layer(1)%ncohort) then
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ivec(1:nc))
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = nc

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%density = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%id      = 0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%biomass = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood_area = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%basal_area = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%heartwood = 0.0_dp
          endif
          ! Layer biomass (summed over cohorts)
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          pop%pop_grid(j)%patch(k)%Layer(1)%biomass = sum(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%biomass)
       enddo
    enddo

    ! recruitment
    if (present(precip)) then
       call layer_recruitment(pop, precip)
    else
       call layer_recruitment(pop)
    endif

    ! Update time since last patch disturbance
    do j=1,np
       do k=1,NPATCH2D

          do idist =1, NDISTURB
             pop%pop_grid(j)%patch(k)%age(idist) = pop%pop_grid(j)%patch(k)%age(idist) + 1
          enddo

       enddo

    enddo

  end subroutine PatchAnnualDynamics


  !*******************************************************************************


  subroutine GetUniqueAgeFrequencies(pop, disturbance_interval, idisturb)

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer(i4b), intent(IN) ::  disturbance_interval(:,:), idisturb

    integer(i4b) :: g, i,j,k,agecopy,idcopy
    real(dp), allocatable :: midpoint(:)
    integer(i4b), allocatable :: ranked_age(:), ranked_age_init(:)
    integer(i4b) :: age_tmp
    integer(i4b), allocatable :: ranked_age_unique_id(:), ranked_age_id(:), counter(:)
    real(dp), allocatable :: tmp(:), freq_tmp(:), freq_tmp1(:)
    real(dp) :: freq
    integer(i4b) :: n_age ! number of unique ages
    integer(i4b) :: npatch_active ! number of active patches
    real(dp):: disturbance_freq
    integer(i4b) :: i_max, Poisson_age(1000), np
    real(dp):: CumPoisson_weight(1000)
    integer(i4b), allocatable :: bound(:,:), unique_age(:)

    !Fills array freq with weights (frequencies across landscape) for each unique age
    ! given specified mean disturbance interval

    np = size(POP%POP_grid)
    do g=1,np

       npatch_active = NPATCH2D
       if (.not.allocated(midpoint)) allocate(midpoint(npatch_active))
       if (.not.allocated(counter)) allocate(counter(npatch_active))
       if (.not.allocated(ranked_age)) allocate(ranked_age(npatch_active))
       if (.not.allocated(ranked_age_init)) allocate(ranked_age_init(npatch_active))
       if (.not.allocated(ranked_age_id)) allocate(ranked_age_id(npatch_active))
       if (.not.allocated(ranked_age_unique_id)) allocate(ranked_age_unique_id(npatch_active))
       if (.not.allocated(tmp)) allocate(tmp(npatch_active))
       if (.not.allocated(freq_tmp)) allocate(freq_tmp(npatch_active))
       if (.not.allocated(freq_tmp1)) allocate(freq_tmp1(npatch_active))


       ! rank patches in order of age
       pop%pop_grid(g)%ranked_age_unique(:, idisturb) = 0
       ranked_age_init = pop%pop_grid(g)%patch%age(idisturb)
       ranked_age = pop%pop_grid(g)%patch%age(idisturb)
       ranked_age_id = pop%pop_grid(g)%patch%id
       ranked_age_unique_id = 0
       freq_tmp = 0.0_dp
       freq = 0.0_dp
       pop%pop_grid(g)%freq_ranked_age_unique(:, idisturb) = 0.0_dp
       midpoint = 0.0_dp


       do i = 1, npatch_active -1
          do j = i+1, npatch_active
             if (ranked_age(i).gt.ranked_age(j)) then
                agecopy          = ranked_age(i)
                idcopy           = ranked_age_id(i)
                ranked_age(i)    = ranked_age(j)
                ranked_age_id(i) = ranked_age_id(j)
                ranked_age(j)    = agecopy
                ranked_age_id(j) = idcopy
             endif
          enddo
       enddo

       ! subset to unique ages
       k=0
       age_tmp = -1
       do i = 1, npatch_active
          if (ranked_age(i).ne.age_tmp) k = k+1
          pop%pop_grid(g)%ranked_age_unique(k, idisturb) = ranked_age(i)
          ranked_age_unique_id(k) = ranked_age_id(i)
          age_tmp = ranked_age(i)
          n_age  = k
       enddo

       disturbance_freq=1.0_dp/real(disturbance_interval(g,idisturb),dp)
       do i =1,1000
          Poisson_age(i) = i
          CumPoisson_weight(i) = CumExponential(disturbance_freq,real(i,dp))

       enddo


       ! construct upper and lower bounds for each unique age: these set the range of ages to be
       ! represented by an unique age
       allocate(bound(n_age,2))
       allocate (unique_age(n_age))
       bound = 0
       unique_age = pop%pop_grid(g)%ranked_age_unique(1:n_age,idisturb)
       do i=1,n_age
          if (unique_age(i).eq.0) then
             bound(i,1) = 0
             bound(i,2) = 0
          elseif ((i.eq.1).and.(unique_age(i).gt.0)) then
             bound(i,1) = 0
             bound(i,2) = unique_age(i)
          elseif ((unique_age(i).gt.0).and.(i.gt.1).and.(unique_age(i-1).eq.unique_age(i)-1)) then
             bound(i,1) = unique_age(i)
             if (i.lt.n_age) then
                bound(i,2) = unique_age(i)
             else
                i_max = maxloc(Poisson_age, 1, CumPoisson_weight.le.0.99_dp)
                bound(i, 2) = Poisson_age(i_max)
             endif
          elseif ((unique_age(i).gt.0).and.(i.gt.1).and.(unique_age(i-1).ne.unique_age(i)-1)) then
             bound(i,1) = bound(i-1,2)+1
             if (i.lt.n_age) then
                bound(i,2) = (unique_age(i)+ unique_age(i+1))/2
             else
                i_max = maxloc(Poisson_age, 1, CumPoisson_weight.le.0.99_dp)
                bound(i, 2) = Poisson_age(i_max)
             endif
          endif


       enddo

       ! calculate weighting for each unique age
       do i=1,n_age
          do j = bound(i,1),bound(i,2)

             !IF (pop%LU(g)==2) THEN  ! secondary forest
             if (POP%pop_grid(g)%LU ==2) then
                freq_tmp(i) = freq_tmp(i) +  pop%pop_grid(g)%freq_age(j+1)
             else
                freq_tmp(i) = freq_tmp(i) + REALExponential(disturbance_freq,real(j,dp))
             endif

          enddo
       enddo

       pop%pop_grid(g)%freq_ranked_age_unique(1:npatch_active,idisturb) = freq_tmp
       pop%pop_grid(g)%n_age(idisturb) = n_age

       deallocate (bound)
       deallocate (unique_age)

    enddo

  end subroutine GetUniqueAgeFrequencies


  !*******************************************************************************


  subroutine GetPatchFrequencies(pop)

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP

    integer(i4b) :: n1, n2, g, REPCOUNT, np, idist
    real(dp) ::  sum_freq

    np = size(Pop%pop_grid)
    do g=1,np

       pop%pop_grid(g)%freq = 0.0_dp
       do idist = 1, NDISTURB
          if (idist.eq.1) then
             do n1=1,pop%pop_grid(g)%n_age(1)
                repcount = count(pop%pop_grid(g)%patch(:)%age(1).eq.pop%pop_grid(g)%ranked_age_unique(n1,1))
                where (pop%pop_grid(g)%patch(:)%age(1).eq.pop%pop_grid(g)%ranked_age_unique(n1,1))
                   pop%pop_grid(g)%freq = pop%pop_grid(g)%freq_ranked_age_unique(n1,1) /real(repcount,dp)
                ENDWHERE
             enddo


          elseif (idist.eq.2) then
             ! first calculate weights for patches with age(2)>age(1)
             do n1=1,pop%pop_grid(g)%n_age(1)


                do n2=1,pop%pop_grid(g)%n_age(idist)
                   repcount = count((pop%pop_grid(g)%patch(1:NPATCH)%age(1) .eq. &
                        pop%pop_grid(g)%ranked_age_unique(n1,1)).and. &
                        (pop%pop_grid(g)%patch(1:NPATCH)%age(idist) .eq.  &
                        pop%pop_grid(g)%ranked_age_unique(n2,idist)))
                   where ((pop%pop_grid(g)%patch(1:NPATCH)%age(1).eq.pop%pop_grid(g)%ranked_age_unique(n1,1)).and. &
                        (pop%pop_grid(g)%patch(1:NPATCH)%age(idist).eq.pop%pop_grid(g)%ranked_age_unique(n2,idist)))
                      pop%pop_grid(g)%freq(1:NPATCH) = pop%pop_grid(g)%freq_ranked_age_unique(n1,1)* &
                           pop%pop_grid(g)%freq_ranked_age_unique(n2,idist) &
                           /real(repcount,dp)

                   ENDWHERE


                enddo
             enddo


          endif
       enddo ! end loop over idist

       sum_freq = sum(pop%pop_grid(g)%freq)
       if (sum_freq.gt.0.0_dp) pop%pop_grid(g)%freq = pop%pop_grid(g)%freq/sum_freq

    enddo


  end subroutine GetPatchFrequencies


  !*******************************************************************************


  subroutine GetDiagnostics(pop,LAI,Cleaf,Croot,disturbance_interval, it, precip)
    ! Gets diagnostic data for current landscape structure

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    real(dp), intent(IN) ::  LAI(:)
    real(dp), intent(IN) ::  Cleaf(:)
    real(dp), intent(IN) ::  Croot(:)
    integer(i4b), intent(IN)        ::  disturbance_interval(:,:)
    real(dp), intent(IN), optional :: precip(:)
    integer(i4b), intent(IN) :: it(:)
    integer(i4b) :: P, g,i,j,ct, ct_highres
    real(dp) :: limits(HEIGHT_BINS+1)
    real(dp) :: ht, cmass_stem,densindiv, freq, freq_old
    character(len=12) :: string1, string2
    character(len=9) :: fmt
    integer(i4b) :: npatch_active  ! number of active patches
    integer(i4b) :: np, nc, i_height
    real(dp) :: diam,basal, cump
    real(dp) :: patch_crown_area(NPATCH2D), patch_crown_cover(NPATCH2D)
    real(dp), allocatable :: height_list(:), height_list_weight(:)
    real(dp) :: height_copy, weight_copy, Pwc, FAVD
    integer(i4b), parameter :: HEIGHT_BINS_highres=100 ! bins for assessing height_max
    real(dp), allocatable :: limits_highres(:), DENSINDIV_HIGHRES(:)
    real(dp) :: tmp2
    integer :: arg1

    fmt = '(f5.1)'
    limits(1) = 0.0_dp
    if(.not.allocated(limits_highres)) allocate(limits_highres(HEIGHT_BINS_highres+1))
    if(.not.allocated(DENSINDIV_HIGHRES)) allocate(DENSINDIV_HIGHRES(HEIGHT_BINS_highres))

    limits_highres(1) = 0.0_dp
    np = size(Pop%pop_grid)

    do g=1, np
       npatch_active = NPATCH2D
       if (MAX_HEIGHT_SWITCH.eq.1) then
          allocate(height_list(NPATCH2D*NCOHORT_MAX))
          allocate(height_list_weight(NPATCH2D*NCOHORT_MAX))
       endif
       !  IF(.NOT.ALLOCATED(MASK)) ALLOCATE(MASK(POP%pop_grid%npatch_active))

       do i=1,HEIGHT_BINS
          limits(i+1) = BIN_POWER**real(i,dp)
          write(string1,fmt) (limits(i))
          write(string2,fmt) (limits(i+1))
          pop%pop_grid(g)%bin_labels(i) = 'Height_'//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'m'
          pop%pop_grid(g)%cmass_stem_bin(i) = 0.0_dp
          pop%pop_grid(g)%densindiv_bin(i) = 0.0_dp
          pop%pop_grid(g)%cmass_stem_bin(i) = 0.0_dp
          pop%pop_grid(g)%height_bin(i) = real(limits(i)+limits(i+1),dp)/2.0_dp
          pop%pop_grid(g)%diameter_bin(i) = ( (real(limits(i),dp)/Kbiometric)**(3.0_dp/2.0_dp) + &
               (real(limits(i+1),dp)/Kbiometric)**(3.0_dp/2.0_dp) ) / 2.0_dp
       enddo

       do i=1,HEIGHT_BINS_highres
          limits_highres(i+1) = real(i,dp)
       enddo

       if (MAX_HEIGHT_SWITCH.eq.1) then
          height_list = 0.0_dp
          height_list_weight = 0.0_dp
       endif
       i_height = 0
       pop%pop_grid(g)%cmass_sum_old      = pop%pop_grid(g)%cmass_sum
       pop%pop_grid(g)%csapwood_sum_old   = pop%pop_grid(g)%csapwood_sum
       pop%pop_grid(g)%cmass_sum          = 0.0_dp
       pop%pop_grid(g)%csapwood_sum       = 0.0_dp
       pop%pop_grid(g)%cheartwood_sum     = 0.0_dp
       pop%pop_grid(g)%height_mean        = 0.0_dp
       pop%pop_grid(g)%fire_mortality     = 0.0_dp
       pop%pop_grid(g)%cat_mortality      = 0.0_dp
       pop%pop_grid(g)%res_mortality      = 0.0_dp
       pop%pop_grid(g)%stress_mortality   = 0.0_dp
       pop%pop_grid(g)%crowding_mortality = 0.0_dp
       pop%pop_grid(g)%sapwood_loss       = 0.0_dp
       pop%pop_grid(g)%sapwood_area_loss  = 0.0_dp
       pop%pop_grid(g)%growth             = 0.0_dp
       pop%pop_grid(g)%area_growth        = 0.0_dp
       pop%pop_grid(g)%basal_area         = 0.0_dp
       pop%pop_grid(g)%densindiv          = 0.0_dp
       pop%pop_grid(g)%height_max         = 0.0_dp
       pop%pop_grid(g)%crown_cover        = 0.0_dp
       pop%pop_grid(g)%crown_area         = 0.0_dp
       pop%pop_grid(g)%sapwood_area       = 0.0_dp
       pop%pop_grid(g)%crown_volume       = 0.0_dp
       densindiv_highres = 0.0_dp
       ! loop through patches
       do P = 1, npatch_active
          pop%pop_grid(g)%patch(p)%biomass          = 0.0_dp
          pop%pop_grid(g)%patch(p)%sapwood          = 0.0_dp
          pop%pop_grid(g)%patch(p)%sapwood_area     = 0.0_dp
          pop%pop_grid(g)%patch(p)%heartwood        = 0.0_dp
          pop%pop_grid(g)%patch(p)%layer(1)%biomass = 0.0_dp
          pop%pop_grid(g)%patch(p)%layer(1)%density = 0.0_dp
          patch_crown_area(p)  = 0.0_dp
          patch_crown_cover(p) = 0.0_dp
          tmp2 = sum(pop%pop_grid(g)%patch(p)%layer(1)%cohort(1:pop%pop_grid(g)%patch(p)%layer(1)%ncohort)%sapwood_area)

          freq =  pop%pop_grid(g)%freq(pop%pop_grid(g)%patch(p)%id)
          freq_old = pop%pop_grid(g)%freq_old(pop%pop_grid(g)%patch(p)%id)

          ! loop through cohorts
          do i = 1, pop%pop_grid(g)%patch(p)%layer(1)%ncohort
             cmass_stem = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%biomass
             densindiv = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%density

             if ( present(precip) ) then
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal, precip(g))
             else
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal )

             endif

             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%height   = ht
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter = diam

             ! basal area in each cohort
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%basal_area = basal

             ! sapwood area in each cohort
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area = basal - & ! m2 ha-1
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%heartwood/(ht*WD)*1.0e4_dp

             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area = &
                  max(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area, 0.0_dp)

             ! get bin
             ct = 1
             do j=1,HEIGHT_BINS
                if (ht.gt.limits(j)) ct = j
             enddo ! bins

             ! get high res bin
             ct_highres = 1
             do j=1,HEIGHT_BINS_highres
                if (ht.gt.limits_highres(j)) ct_highres = j
             enddo ! bins

             pop%pop_grid(g)%patch(p)%layer(1)%biomass = pop%pop_grid(g)%patch(p)%layer(1)%biomass + cmass_stem
             pop%pop_grid(g)%patch(p)%layer(1)%density = pop%pop_grid(g)%patch(p)%layer(1)%density + densindiv

             if (diam*100.0_dp .gt. 1.0_dp) then
                patch_crown_area(p) = patch_crown_area(p) + densindiv*PI*(diam*100.0_dp*0.1492_dp)**2 ! uses GC relationship
                pop%pop_grid(g)%crown_volume = pop%pop_grid(g)%crown_volume + &
                     freq*densindiv*(4.0_dp/3.0_dp)*PI*(diam*100.0_dp*0.1492_dp)**2*(1.5_dp*(diam*100.0_dp*0.1492_dp))
                ! assumes vertical radius = 1.5 * horizontal radius
             endif

             if (diam*100.0_dp .gt. 5.0_dp) then
                if (ALLOM_SWITCH.eq.1) then
                   !! assumes crown radius (m) = 0.1492 * dbh (cm) (from G. Cook, pers. comm.)
                   ! assumes vertical radius = 1.5 * horizontal radius
                   pop%pop_grid(g)%crown_volume = pop%pop_grid(g)%crown_volume + &
                        freq*densindiv*(4.0_dp/3.0_dp)*PI*(diam*100.0_dp*0.1492_dp)**2*(1.5_dp*(diam*100.0_dp*0.1492_dp))
                else
                   !! global allometry
                   ! assumes vertical radius = 1.5 * horizontal radius
                   pop%pop_grid(g)%crown_volume = pop%pop_grid(g)%crown_volume + &
                        freq*densindiv*(4.0_dp/3.0_dp)*PI*1.5_dp*((k_allom1 * diam ** k_rp )/PI)**1.5_dp
                endif

             endif

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
             if (MAX_HEIGHT_SWITCH.eq.1) then
                i_height = i_height+1
                height_list(i_height) = ht
                height_list_weight(i_height) = densindiv*freq
             endif
          enddo ! cohorts

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

          pop%pop_grid(g)%area_growth =  pop%pop_grid(g)%area_growth + &
               freq*pop%pop_grid(g)%patch(p)%area_growth

       enddo ! patches

       if (INTERP_SWITCH==1.and.NDISTURB.eq.2) then
          !CALL INTERPOLATE_BIOMASS_2D(pop, disturbance_interval,it)
          call INTERPOLATE_BIOMASS_2D(pop, disturbance_interval,it(g),g)
       elseif (INTERP_SWITCH==1.and.NDISTURB.eq.1) then
          call INTERPOLATE_BIOMASS_1D(pop, disturbance_interval,it(g),g)
       endif

       arg1 = NYEAR_HISTORY
       if (SMOOTH_SWITCH==1) then
          if (it(g).le.NYEAR_HISTORY) then
             call SMOOTH_FLUX(POP,g,it(g))
          else
             call SMOOTH_FLUX(POP,g,int(arg1,i4b))
          endif
       endif

       ! leaf area index in each cohort
       do P = 1, npatch_active

          freq =  pop%pop_grid(g)%freq(pop%pop_grid(g)%patch(p)%id)
          ! loop through cohorts
          do i = 1, pop%pop_grid(g)%patch(p)%layer(1)%ncohort
             cmass_stem = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%biomass
             densindiv = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%density
             basal=PI*(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter/2.0_dp)* &
                  (pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter/2.0_dp)*densindiv*1.0e4_dp
             ! leaf area index in each cohort
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%LAI = LAI(g) * &
                  min(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area  &
                  /max(pop%pop_grid(g)%sapwood_area,1.0e-3_dp), 10.0_dp)
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Cleaf = Cleaf(g) * &
                  min(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area  &
                  /max(pop%pop_grid(g)%sapwood_area,1.0e-3_dp), 10.0_dp)
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Croot = Croot(g) * &
                  min(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood_area  &
                  /max(pop%pop_grid(g)%sapwood_area,1.0e-3_dp), 10.0_dp)
          enddo ! cohorts
          pop%pop_grid(g)%patch(p)%LAI = sum(pop%pop_grid(g)%patch(p)%layer(1)% &
               cohort(1:pop%pop_grid(g)%patch(p)%layer(1)%ncohort)%LAI)
       enddo ! patches

       ! PGap = (1-fcover) calculation

       if (pop%pop_grid(g)%crown_volume>0.0_dp) then
          FAVD = LAI(g)/pop%pop_grid(g)%crown_volume ! foliage area volume density
       else
          FAVD = 0.0_dp
       endif

       do P = 1, npatch_active
          freq =  pop%pop_grid(g)%freq(pop%pop_grid(g)%patch(p)%id)
          nc =  pop%pop_grid(g)%patch(p)%layer(1)%ncohort
          ! loop through cohorts
          do i = 1, nc
             cmass_stem = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%biomass
             densindiv = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%density

             if ( present(precip) ) then
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal, precip(g))
             else
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass_stem, densindiv, ht, diam, basal )
             endif


             if (diam*100.0_dp .gt. 1.0_dp) then

                if (ALLOM_SWITCH.eq.1) then
                   !! assumes crown radius (m) = 0.1492 * dbh (cm) (from G. Cook, pers. comm.)
                   pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = densindiv*PI*(diam*100.0_dp*0.1492_dp)**2
                   Pwc = exp(-0.5_dp * pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%LAI/ &
                        pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area)
                   pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = &
                        densindiv*PI*(diam*100.0_dp*0.1492_dp)**2*(1.0_dp-Pwc)

                else
                   !! global allometry
                   pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = &
                        densindiv*PI*(((k_allom1 * diam ** k_rp )/PI)**0.5_dp)**2
                   Pwc = exp(max(-0.5_dp * pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%LAI/ &
                        max(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area,1.e-3_dp),-20.0_dp))

                   pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = &
                        pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area*(1.0_dp-Pwc) !*1.4142
                endif

             else
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area = &
                     0.5_dp*pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%LAI !*1.4142
             endif
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area= &
                  max(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area,0.01_dp)

             pop%pop_grid(g)%crown_area = pop%pop_grid(g)%crown_area + &
                  freq*pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area

             if (i.eq.1) then ! top cohort
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Pgap = &
                     exp(-pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area)

                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%frac_interception = &
                     1- exp(-pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area)
             else
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Pgap = &
                     pop%pop_grid(g)%patch(p)%layer(1)%cohort(i-1)%Pgap* &
                     exp(-pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%crown_area)
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%frac_interception = &
                     pop%pop_grid(g)%patch(p)%layer(1)%cohort(i-1)%Pgap - &
                     pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Pgap
             endif
             pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%respiration_scalar =  &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%sapwood/shootfrac/CtoNw + &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Cleaf/CtoNl + &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%Croot/CtoNr

          enddo ! cohorts

          if (nc>0) then
             pop%pop_grid(g)%patch(p)%pgap  = &
                  pop%pop_grid(g)%patch(p)%layer(1)%cohort(nc)%Pgap
          else
             pop%pop_grid(g)%patch(p)%pgap  = 1
          endif


       enddo ! patches

       pop%pop_grid(g)%Kclump = max(pop%pop_grid(g)%crown_area/(0.5_dp*LAI(g)),0.1_dp)
       pop%pop_grid(g)%crown_cover = 1.0_dp-exp(-pop%pop_grid(g)%crown_area)


       pop%pop_grid(g)%height_mean = pop%pop_grid(g)%height_mean/max(pop%pop_grid(g)%densindiv,1.0e-5_dp)

       ! Height Diagnostics
       if (MAX_HEIGHT_SWITCH.eq.0) then
          ! Set landscape maximum height to centre of bin with <5% of trees in a bin of higher size classes
          cump = 0.0_dp
          j = 1
          do while (cump.lt.0.95_dp)
             cump = cump + pop%pop_grid(g)%densindiv_bin(j)/max(pop%pop_grid(g)%densindiv,1.0e-5_dp)
             pop%pop_grid(g)%height_max = pop%pop_grid(g)%height_bin(j)
             j = j+1
          enddo
       elseif (MAX_HEIGHT_SWITCH.eq.1) then

          ! sort height list
          do i = 1, i_height -1
             do j = i+1, i_height
                if (height_list(i).gt.height_list(j)) then
                   height_copy = height_list(i)
                   weight_copy =  height_list_weight(i)
                   height_list(i) = height_list(j)
                   height_list_weight(i) = height_list_weight(j)
                   height_list(j) = height_copy
                   height_list_weight(j) = weight_copy
                endif
             enddo
          enddo ! end sort height list

          ! normailse height list weights
          height_list_weight=height_list_weight/sum(height_list_weight(1:i_height))
          cump = 0.0_dp
          j = 1
          do while (cump.lt.0.95_dp)
             cump = cump + height_list_weight(j)
             pop%pop_grid(g)%height_max = height_list(j)
             j = j+1
          enddo
          deallocate(height_list)
          deallocate(height_list_weight)

       elseif (MAX_HEIGHT_SWITCH.eq.2) then
          cump = 0.0_dp
          j = 1
          densindiv_highres= densindiv_highres/max(sum(densindiv_highres),1.0e-5_dp)
          do while ((cump.lt.0.95_dp).and.(j.le.HEIGHT_BINS_highres))
             cump = cump + densindiv_highres(j)
             pop%pop_grid(g)%height_max = (limits_highres(j+1) + limits_highres(j))/2.0_dp
             j = j+1
          enddo
       endif
       !deallocate(MASK)

    enddo ! end loop over grid cells


  end subroutine GetDiagnostics


  !*******************************************************************************


  subroutine Patch_partial_disturb(pop,idisturb,intensity,frac_intensity1)

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer(i4b), intent(IN) ::  idisturb
    real(dp), intent(IN) :: intensity(:,:)
    real(dp), intent(IN), optional :: frac_intensity1(:)

    integer(i4b) :: j, k, c, nc, np
    integer(i4b) ::  ivec(NCOHORT_MAX)
    real(dp) :: ht, diam
    real(dp) :: Psurvival_s, Psurvival, char_height

    np = size(Pop%pop_grid)

    ! Kills a fraction of biomass in patch when prescribed disturbance interval is reached
    do j=1,np
       do k=1,NPATCH
          pop%pop_grid(j)%patch(k)%fire_mortality = 0.0_dp

          if (((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).ne.0).and. &
               (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).eq.pop%pop_grid(j)%patch(k)%age(idisturb))).or. &
               (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).eq.pop%pop_grid(j)%patch(k)%age(idisturb))) then


             ! loop through cohorts
             ivec = 0
             nc = 0
             do c = 1, pop%pop_grid(j)%patch(k)%layer(1)%ncohort
                ! kill fraction of each cohort
                char_height = 3.7_dp*(1.0_dp-exp(-0.19_dp*Intensity(j,1)))
                ht = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%height
                diam = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter*100.0_dp ! diameter in cm
                if ((ht.gt.8.5_dp).and.(ht.gt.char_height)) then
                   Psurvival_s  =(-0.0011_dp*Intensity(j,1) -0.00002_dp)*ht &
                        +(0.0075_dp*Intensity(j,1)+1.0_dp)

                elseif ((ht.le.8.5_dp).and.(ht.gt.char_height)) then
                   Psurvival_s =(0.0178_dp*Intensity(j,1) + 0.0144_dp)*ht &
                        + (-0.1174_dp*Intensity(j,1)+0.9158_dp)

                else
                   Psurvival_s = 0.0_dp
                endif
                Psurvival_s = min(Psurvival_s,1.0_dp)
                Psurvival_s = max(Psurvival_s,1.0e-3_dp)
                Psurvival = Psurvival_s


                if (present(frac_intensity1)) then
                   char_height = 3.7_dp*(1.0_dp-exp(-0.19_dp*Intensity(j,2)))
                   ht = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%height
                   if ((ht.gt.8.5_dp).and.(ht.gt.char_height)) then
                      Psurvival_s  =(-0.0011_dp*Intensity(j,2) -0.00002_dp)*ht &
                           +(0.0075_dp*Intensity(j,2)+1.0_dp)
                   elseif ((ht.le.8.5_dp).and.(ht.gt.char_height)) then
                      Psurvival_s =(0.0178_dp*Intensity(j,2) + 0.0144_dp)*ht &
                           + (-0.1174_dp*Intensity(j,2)+0.9158_dp)
                   else
                      Psurvival_s = 0.0_dp
                   endif
                   Psurvival_s = min(Psurvival_s,1.0_dp)
                   Psurvival_s = max(Psurvival_s,1.0e-3_dp)
                   Psurvival = Psurvival_s*(1.0_dp-frac_intensity1(j)) + Psurvival*frac_intensity1(j)
                endif
                ! Psurvival = 1.0_dp ! test
                pop%pop_grid(j)%patch(k)%fire_mortality = pop%pop_grid(j)%patch(k)%fire_mortality + &
                     (1.0_dp-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                     (1.0_dp-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                     (1.0_dp-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood_area
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%heartwood = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%heartwood
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density
                if (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.lt.DENSINDIV_MIN) then
                   ! remove cohort
                   pop%pop_grid(j)%patch(k)%fire_mortality = pop%pop_grid(j)%patch(k)%fire_mortality + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                   pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                   pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood_area
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = 0.0_dp
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = 0.0_dp
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood = 0.0_dp
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood = 0.0_dp
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood_area = 0.0_dp

                else
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 1
                   nc = nc+1
                   ivec(nc)=c
                endif
             enddo

             ! SHUFFLE if necessary to remove zero-density cohorts
             if (nc.lt.pop%pop_grid(j)%patch(k)%Layer(1)%ncohort) then
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ivec(1:nc))
                pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = nc

                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%density = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%id = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%biomass = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood_area = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%heartwood = 0.0_dp
             endif

             pop%pop_grid(j)%patch(k)%age(idisturb) = 0
             pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0


          endif

       enddo
    enddo

  end subroutine Patch_partial_disturb


  !*******************************************************************************


  subroutine Patch_partial_disturb2(pop,idisturb)

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer(i4b), intent(IN) ::  idisturb

    integer(i4b) :: j, k, c, nc, np
    integer(i4b) ::  ivec(NCOHORT_MAX)
    real(dp) :: Psurvival, frac_mort, Pmort

    np = size(Pop%pop_grid)

    ! Kills a fraction (80%) biomass in patch when prescribed disturbance interval is reached
    do j=1,np
       do k=1,NPATCH2D
          pop%pop_grid(j)%patch(k)%cat_mortality = 0.0_dp
          ! Layer biomass (summed over cohorts)
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          pop%pop_grid(j)%patch(k)%Layer(1)%biomass = sum(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%biomass)

          if (((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).ne.0).and. &
               (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).eq.pop%pop_grid(j)%patch(k)%age(idisturb))).or. &
               (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).eq.pop%pop_grid(j)%patch(k)%age(idisturb))) then


             ! loop through cohorts
             ivec = 0
             nc = 0
             frac_mort = 0.0_dp
             pop%pop_grid(j)%patch(k)%cat_mortality = 0.0_dp
             do c = 1, pop%pop_grid(j)%patch(k)%layer(1)%ncohort
                ! kill fraction of each cohort, up to 80% of patch biomass
                if (pop%pop_grid(j)%patch(k)%cat_mortality  < 0.8_dp * pop%pop_grid(j)%patch(k)%Layer(1)%biomass ) then
                   Pmort = min(  (0.8_dp*pop%pop_grid(j)%patch(k)%Layer(1)%biomass-pop%pop_grid(j)%patch(k)%fire_mortality) &
                        /pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass, 1.0_dp)
                else
                   Pmort = 0.0_dp
                endif
                Psurvival = 1.0_dp - Pmort

                pop%pop_grid(j)%patch(k)%cat_mortality = pop%pop_grid(j)%patch(k)%cat_mortality + &
                     (1.0_dp-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                     (1.0_dp-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                     (1.0_dp-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood_area
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%heartwood = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%heartwood
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density = &
                     Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density
                if (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.lt.DENSINDIV_MIN) then
                   ! remove cohort
                   pop%pop_grid(j)%patch(k)%cat_mortality = pop%pop_grid(j)%patch(k)%cat_mortality + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                   pop%pop_grid(j)%patch(k)%sapwood_loss = pop%pop_grid(j)%patch(k)%sapwood_loss + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood
                   pop%pop_grid(j)%patch(k)%sapwood_area_loss = pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%sapwood_area
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = 0.0_dp
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = 0.0_dp
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%heartwood = 0.0_dp
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood = 0.0_dp
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%sapwood_area = 0.0_dp

                else
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 1
                   nc = nc+1
                   ivec(nc)=c
                endif
             enddo

             ! SHUFFLE if necessary to remove zero-density cohorts
             if (nc.lt.pop%pop_grid(j)%patch(k)%Layer(1)%ncohort) then
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ivec(1:nc))
                pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = nc

                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%density = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%id = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%biomass = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%sapwood_area = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%heartwood = 0.0_dp
             endif

             pop%pop_grid(j)%patch(k)%age(idisturb) = 0
             pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0


          endif

       enddo
    enddo

  end subroutine Patch_partial_disturb2


  !*******************************************************************************


  subroutine Patch_disturb(pop,idisturb,precip)
    implicit none

    type(POP_TYPE), intent(INOUT)  :: POP
    real(dp), intent(IN), optional :: precip(:)
    !INTEGER(i4b), INTENT(IN) :: it(:),idisturb
    integer(i4b), intent(IN) :: idisturb
    integer(i4b) :: j, k, np, nc

    np = size(Pop%pop_grid)
    ! Kills all biomass in patch when prescribed disturbance interval is reached
    ! Should be called after accounting for this year

    do j=1,np
       do k=1,NPATCH2D
          pop%pop_grid(j)%patch(k)%cat_mortality = 0.0_dp
          if (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).ne.0) then
             if ((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).eq.pop%pop_grid(j)%patch(k)%age(idisturb)).or. &
                  (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).eq.pop%pop_grid(j)%patch(k)%age(idisturb)) ) then
                ! kill entire layer
                nc = pop%pop_grid(j)%patch(k)%layer(1)%ncohort

                ! pop%pop_grid(j)%patch(k)%fire_mortality = SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%biomass)
                pop%pop_grid(j)%patch(k)%cat_mortality = sum(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%biomass)
                pop%pop_grid(j)%patch(k)%sapwood_loss =  pop%pop_grid(j)%patch(k)%sapwood_loss + &
                     sum(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%sapwood)
                pop%pop_grid(j)%patch(k)%sapwood_area_loss =  pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                     sum(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%sapwood_area)
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%ncohort = 0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%biomass = 0.0_dp
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%density = 0.0_dp
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmean = 0.0_dp
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmax  = 0.0_dp
                pop%pop_grid(j)%patch(k)%age(idisturb) = 0
                pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0

                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%density = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%id = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%biomass = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%sapwood = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%sapwood_area = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%heartwood = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%age = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%lai = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%height = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%pgap = 1.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_resource_uptake = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_light_uptake = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_interception = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_respiration = 0.0_dp
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_NPP = 0.0_dp
                pop%pop_grid(j)%patch(k)%area_growth = 0.0_dp
                pop%pop_grid(j)%patch(k)%pgap = 1.0_dp
                ! understorey recruitment
                if (present(precip)) then
                   call layer_recruitment_single_patch(pop,k,j,precip)
                else
                   call layer_recruitment_single_patch(pop,k,j)

                endif
             endif
          elseif (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).eq.pop%pop_grid(j)%patch(k)%age(idisturb)) then
             ! kill entire layer
             nc = pop%pop_grid(j)%patch(k)%layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%sapwood_loss =  pop%pop_grid(j)%patch(k)%sapwood_loss + &
                  sum(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%sapwood)
             pop%pop_grid(j)%patch(k)%sapwood_area_loss =  pop%pop_grid(j)%patch(k)%sapwood_area_loss + &
                  sum(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%sapwood_area)
             pop%pop_grid(j)%patch(k)%cat_mortality = sum(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%biomass)
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%ncohort = 0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%biomass = 0.0_dp
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%density = 0.0_dp
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmean = 0.0_dp
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmax  = 0.0_dp
             pop%pop_grid(j)%patch(k)%age(idisturb) = 0
             pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%density = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%id = 0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%biomass = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%sapwood = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%sapwood_area = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%heartwood = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%age = 0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%lai = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%height = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%pgap = 1.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_resource_uptake = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_light_uptake = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_interception = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_respiration = 0.0_dp
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%frac_NPP = 0.0_dp
             pop%pop_grid(j)%patch(k)%area_growth = 0.0_dp
             pop%pop_grid(j)%patch(k)%pgap = 1.0_dp
             ! understorey recruitment
             if (present(precip)) then
                call layer_recruitment_single_patch(pop,k,j,precip)
             else
                call layer_recruitment_single_patch(pop,k,j)

             endif
          endif

       enddo

    enddo

  end subroutine Patch_disturb


  !*******************************************************************************


  subroutine  layer_recruitment(pop,precip)

    implicit none

    type(POP_TYPE), intent(INOUT)  :: POP
    real(dp), intent(IN), optional :: precip(:)

    real(dp) :: f, mu, densindiv, cmass, ht
    integer(i4b) :: j, k, ncohort, np
    real(dp) :: diam, basal

    np = size(Pop%pop_grid)

    do j=1,np
       do k=1,NPATCH2D
          if (RECRUIT_SWITCH==0) then
             pop%pop_grid(j)%patch(k)%factor_recruit = exp(-0.6_dp*((pop%pop_grid(j)%patch(k)%Layer(1)%biomass)**(0.6667_dp)))
          elseif (RECRUIT_SWITCH==1) then
             pop%pop_grid(j)%patch(k)%factor_recruit = max(pop%pop_grid(j)%patch(k)%pgap,1.0e-3_dp)
          endif
          f = pop%pop_grid(j)%patch(k)%factor_recruit
          mu=exp(max(FULTON_ALPHA*(1.0_dp-2.0_dp*THETA_recruit/ &
               (f+1.0_dp-sqrt((f+1.0_dp)*(f+1.0_dp)-4.0_dp*THETA_recruit*f))), &
               -50.0_dp))
          densindiv=DENSINDIV_MAX*mu + pop%pop_grid(j)%patch(k)%fire_top_kill_density
          cmass=CMASS_STEM_INIT*densindiv/DENSINDIV_MAX

          !write(5599,*),  pop%pop_grid(j)%patch(k)%fire_top_kill_density,  densindiv, pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          !COMMLN below: should not be cohort +1 or .LE. !
          if (cmass>EPS*10.0_dp .and. densindiv>DENSINDIV_MIN .and. &
               (pop%pop_grid(j)%patch(k)%Layer(1)%ncohort+1).lt.NCOHORT_MAX) then

             ! create a new cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort + 1
             ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%biomass = cmass
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%density = densindiv
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%sapwood = cmass

             if ( present(precip) ) then
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass, densindiv, ht, diam, basal, precip(j))
             else
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass, densindiv, ht, diam, basal )
             endif

             pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%height   = ht
             pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%diameter = diam

          endif
          pop%pop_grid(j)%patch(k)%fire_top_kill_density = 0.0_dp
       enddo
    enddo

  end subroutine layer_recruitment


  !*******************************************************************************


  subroutine  layer_recruitment_single_patch(pop, index, grid_index,precip)

    implicit none

    type(POP_TYPE), intent(INOUT)  :: POP
    real(dp), intent(IN), optional :: precip(:)
    integer(i4b), intent(IN) :: index, grid_index

    real(dp) :: f, mu, densindiv, cmass, ht
    integer(i4b) :: j, k, ncohort, np
    real(dp) :: diam,basal

    np = size(Pop%pop_grid)
    do j=grid_index,grid_index
       do k=index,index
          if (RECRUIT_SWITCH==0) then
             pop%pop_grid(j)%patch(k)%factor_recruit = exp(-0.6_dp*((pop%pop_grid(j)%patch(k)%Layer(1)%biomass)**(0.6667_dp)))
          elseif (RECRUIT_SWITCH==1) then
             !pop%pop_grid(j)%patch(k)%factor_recruit = pop%pop_grid(j)%patch(k)%pgap
             pop%pop_grid(j)%patch(k)%factor_recruit = 1
          endif
          f = pop%pop_grid(j)%patch(k)%factor_recruit
          mu=exp(FULTON_ALPHA*(1.0_dp-2.0_dp*THETA_recruit/(f+1.0_dp-sqrt((f+1.0_dp)*(f+1.0_dp)-4.0_dp*THETA_recruit*f))))
          densindiv=DENSINDIV_MAX*mu
          cmass=CMASS_STEM_INIT*densindiv/DENSINDIV_MAX

          if (cmass>EPS*10.0_dp .and. densindiv>DENSINDIV_MIN .and. &
               (pop%pop_grid(j)%patch(k)%Layer(1)%ncohort+1).lt.NCOHORT_MAX) then
             ! create a new cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort + 1
             ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%biomass = cmass
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%density = densindiv
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%sapwood = cmass

             if ( present(precip) ) then
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass, densindiv, ht, diam, basal, precip(j))
             else
                call GET_ALLOMETRY( ALLOM_SWITCH,  cmass, densindiv, ht, diam, basal )
             endif

             pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%height   = ht
             pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%diameter = diam


          endif


       enddo
    enddo

  end subroutine layer_recruitment_single_patch

  !*******************************************************************************
  ! Exponential distribution
  ! Returns probability of a given time-between-events (x)
  ! Given a Poisson process with expected frequency (events per unit time) lambda
  ! Reference: http://en.wikipedia.org/wiki/Exponential_distribution
  ! Use to determine average age (x, years) of patches with a given random disturbance
  ! frequency lambda (disturbances per year)

  real(dp) function Exponential(lambda, x)

    implicit none

    integer(i4b), intent(IN) :: x
    real(dp), intent(IN) ::  lambda

    if (x .lt. 0.0_dp) then ! Shouldn't happen but ...
       Exponential = 0.0_dp
    else
       Exponential = lambda*exp(-lambda*x)
    endif

  end function Exponential

  !*******************************************************************************
  ! Exponential distribution
  ! Returns probability of a given time-between-events (x)
  ! Given a Poisson process with expected frequency (events per unit time) lambda
  ! Reference: http://en.wikipedia.org/wiki/Exponential_distribution
  ! Use to determine average age (x, years) of patches with a given random disturbance
  ! frequency lambda (disturbances per year)

  real(dp) function REALExponential(lambda, x)

    implicit none

    real(dp), intent(IN) ::  x
    real(dp), intent(IN) ::  lambda

    if (x .lt. 0.0_dp) then ! Shouldn't happen but ...
       REALExponential = 0.0_dp
    else
       REALExponential = lambda*exp(-lambda*x)
    endif

  end function REALExponential

  !*******************************************************************************

  real(dp) function CumExponential(lambda, x)

    implicit none

    real(dp), intent(IN) :: x
    real(dp), intent(IN) ::  lambda

    if (x .lt. 0.0_dp) then ! Shouldn't happen but ...
       CumExponential = 0.0_dp
    else
       CumExponential = 1.0_dp - exp(-lambda*x)
    endif

  end function CumExponential

  !*******************************************************************************

  real(dp) function Factorial(n)

    implicit none

    integer, intent(IN) :: n

    integer :: i
    real(dp) :: Ans

    Ans = 1.0_dp
    do i = 1, n
       Ans = Ans * real(i,dp)
    end do

    Factorial = Ans

  end function Factorial

  !*******************************************************************************
  ! ALLOMETRY
  !*******************************************************************************

  subroutine GET_ALLOMETRY( ALLOM_SWITCH,  biomass, density, ht, diam, basal, precip )

#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif

    implicit none

    integer(i4b), intent(IN) :: ALLOM_SWITCH
    real(dp),     intent(IN) :: biomass
    real(dp),     intent(IN) :: density
    real(dp),     intent(IN), optional :: precip
    real(dp),     intent(OUT):: ht, diam, basal

#ifdef __MPI__
    integer :: ierr
#endif

    ! Standard Allometry
    if (ALLOM_SWITCH.eq.0) then
       ht   = (Kbiometric**(3.0_dp/4.0_dp))*(4.0_dp*biomass/(max(density,1.0e-5_dp)*WD*PI))**(1.0_dp/4.0_dp)
       diam = (ht/Kbiometric)**(1.5_dp)
       basal= PI * (diam/2.0_dp) * (diam/2.0_dp) * density * 1.0e4_dp

       ! Top-End Allometry following G.Cook
    elseif (ALLOM_SWITCH.eq.1.and.present(precip)) then
       ht   =GetHeight(precip,biomass,density)
       call Allometry(ht,biomass,density,diam,basal)

       ! Allometry following Williams 2005, Model 5b
    elseif ( ALLOM_SWITCH.eq.2 ) then
       call Williams_Allometry(biomass,density,ht,diam,basal)

    else
       write(*,*)"Invalid Allometry settings in POP!"
       write(*,*)"ALLOM_SWITCH   = ",ALLOM_SWITCH
       write(*,*)"Precip present = ",present(precip)
#ifdef __MPI__
       call MPI_Abort(0, 85, ierr) ! Do not know comm nor rank here
#else
       stop 85
#endif
    endif

  end subroutine GET_ALLOMETRY
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

  real(dp) function GetHeight(precip,biomass,density)

    implicit none

    real(dp), intent(IN) :: precip
    real(dp), intent(IN) :: biomass
    real(dp), intent(IN) :: density

    real(dp),parameter:: THETA=0.99_dp ! Shape parameter, should be slightly <1
    real(dp),parameter:: HMIN=0.001_dp ! min bound for tree height
    real(dp),parameter:: HMAX=100.0_dp ! max bound for tree height
    real(dp),parameter:: EPS=0.01_dp ! precision of the root
    integer(i4b), parameter :: MAXTRIES=25

    real(dp) :: alpha,beta,delta,rh,st,x1,x2,rtbis,dx,fmid,xmid,lhs,rhs
    integer(i4b) :: b

    alpha=4.05_dp*exp(-0.00032_dp*precip)
    beta=5.4_dp*exp(0.0014_dp*precip)
    delta=2.0_dp*sqrt(biomass/density/WD/PI)

    x1=HMIN
    x2=HMAX
    rtbis=x1
    dx=x2-x1
    b=0
    fmid=EPS+1.0_dp

    do while (abs(dx).gt.EPS.and.b.le.MAXTRIES)
       b=b+1
       dx=dx*0.5_dp
       xmid=rtbis+dx

       ! Evaluate LHS-RHS at height=xmid
       ! LHS-RHS should increase with increasing height

       lhs=xmid
       rh=1.0_dp/sqrt(xmid)
       st=alpha+beta*delta*rh+100.0_dp*delta*rh
       rhs=1.0_dp/2.0_dp/THETA* &
            (st-sqrt(st*st-400.0_dp*THETA*alpha*delta*rh- &
            400.0_dp*THETA*beta*delta*delta/xmid))
       fmid=lhs-rhs

       if (fmid.lt.0.0_dp) rtbis=xmid

    enddo

    GetHeight=xmid

  end function GetHeight


  !*******************************************************************************


  subroutine INTERPOLATE_BIOMASS_1D(pop, disturbance_interval,it,g)
    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer(i4b), intent(IN) ::  disturbance_interval(:,:)
    integer(i4b), intent(IN) ::  it,g

    integer(i4b) :: nage,iage, i_min, i_max
    integer(i4b) :: i_min_growth, i_max_growth
    real(dp) :: disturbance_freq, tmp_min, tmp_max, tmp1_min, tmp1_max, tmp_array(NPATCH2D)
    real(dp) :: tmp2_min, tmp2_max
    real(dp) :: tmp3_min, tmp3_max
    real(dp) :: tmp4_min, tmp4_max
    logical :: MASK(NPATCH2D)
    integer(i4b) :: age_min, age_max
    integer(i4b) :: age_min_growth, age_max_growth
    integer(i4b), allocatable :: age(:)
    real(dp), allocatable ::cmass_age(:), stress_mort_age(:), crowd_mort_age(:)
    real(dp), allocatable ::csapwood_age(:), sapwood_area_age(:), growth_age(:)
    real(dp), allocatable ::freq_age(:)

    ! get interpolated biomass,sapwood, stress mortality, crowding mortality, disturbance mortality
    POP%pop_grid(g)%cmass_sum= 0.0_dp
    POP%pop_grid(g)%stress_mortality = 0.0_dp
    POP%pop_grid(g)%cat_mortality = 0.0_dp
    pop%pop_grid(g)%crowding_mortality = 0.0_dp
    pop%pop_grid(g)%csapwood_sum = 0.0_dp
    pop%pop_grid(g)%sapwood_area = 0.0_dp
    tmp_array = 0.0_dp
    nage =  min(POP%pop_grid(g)%patch(1)%disturbance_interval(1),it)+1 ! maximum age

    !nage = maxval(pop%pop_grid(g)%patch(:)%age(1))
    if (POP%pop_grid(g)%LU==2) then ! secondary forest
       nage = AGEMAX
       do iage=AGEMAX,1,-1
          if (pop%pop_grid(g)%freq_age(iage)>0) then
             exit
          else
             nage = nage - 1
          endif
       enddo
    endif


    disturbance_freq=1.0_dp/real(disturbance_interval(g,1),dp)
    if(.not.allocated(age)) allocate(age(nage))
    if(.not.allocated(freq_age)) allocate(freq_age(nage))
    if(.not.allocated(cmass_age)) allocate(cmass_age(nage))
    if(.not.allocated(growth_age)) allocate(growth_age(nage))
    if(.not.allocated(csapwood_age)) allocate(csapwood_age(nage))
    if(.not.allocated(sapwood_area_age)) allocate(sapwood_area_age(nage))
    if(.not.allocated(stress_mort_age)) allocate(stress_mort_age(nage))
    if(.not.allocated(crowd_mort_age)) allocate(crowd_mort_age(nage))
    !pop%pop_grid(g)%biomass_age(2:agemax) = pop%pop_grid(g)%biomass_age(1:agemax-1)
    !pop%pop_grid(g)%biomass_age(1) = 0.0
    !cmass_age = pop%pop_grid(g)%biomass_age
    tmp_min = 0.0_dp
    tmp_max = 0.0_dp
    pop%pop_grid(g)%biomass_age = 0.0_dp


    !IF (pop%LU(g)==2) THEN  ! secondary forest
    if (POP%pop_grid(g)%LU==2) then ! secondary forest
       do iage = 1, nage
          age(iage) = iage-1
          freq_age(iage) =  pop%pop_grid(g)%freq_age(iage)
       enddo
    else
       do iage = 1, nage
          age(iage) = iage-1
          freq_age(iage) =  REALExponential(disturbance_freq,real(age(iage),dp))
          pop%pop_grid(g)%freq_age(iage) = freq_age(iage)
       end do
    endif
    if (sum(freq_age)>0.0_dp) freq_age = freq_age/sum(freq_age)

    do iage = 1, nage
       ! get nearest ages bracketing age(iage)
       if (any(pop%pop_grid(g)%patch(:)%age(1).le.age(iage))) then
          age_min = maxval(pop%pop_grid(g)%patch(:)%age(1), &
               pop%pop_grid(g)%patch(:)%age(1).le.age(iage))
          i_min = maxloc(pop%pop_grid(g)%patch(:)%age(1), 1, &
               pop%pop_grid(g)%patch(:)%age(1).le.age(iage))
       else
          age_min = 0
          i_min   = 0
       endif
       if (any(pop%pop_grid(g)%patch(:)%age(1).ge.age(iage))) then
          age_max = minval(pop%pop_grid(g)%patch(:)%age(1), &
               pop%pop_grid(g)%patch(:)%age(1).ge.age(iage))
          i_max = minloc(pop%pop_grid(g)%patch(:)%age(1), 1, &
               pop%pop_grid(g)%patch(:)%age(1).ge.age(iage))
       else
          age_max = 0
          i_max   = 0
       endif

       age_min_growth = age_min
       age_max_growth = age_max
       i_min_growth = i_min
       i_max_growth = i_max

       if ((i_min.gt.0).and.(i_max.gt.0).and.(age_max.eq.age_min)) then
          ! no need to interpolate
          MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_min
          where (MASK)
             tmp_array = 1.0_dp
          elsewhere
             tmp_array = 0.0_dp
          endwhere
          cmass_age(iage) =  &
               sum(pop%pop_grid(g)%patch(:)%layer(1)%biomass,MASK)/sum(tmp_array)
          growth_age(iage) =  &
               sum(pop%pop_grid(g)%patch(:)%growth,MASK)/sum(tmp_array)
          csapwood_age(iage) = sum(pop%pop_grid(g)%patch(:)%sapwood,MASK)/sum(tmp_array)
          sapwood_area_age(iage) = &
               sum(pop%pop_grid(g)%patch(:)%sapwood_area,MASK)/sum(tmp_array)
          stress_mort_age(iage)= &
               sum(pop%pop_grid(g)%patch(:)%stress_mortality,MASK)/sum(tmp_array)
          crowd_mort_age(iage)= &
               sum(pop%pop_grid(g)%patch(:)%crowding_mortality,MASK)/sum(tmp_array)
       else
          ! interpolate or extrapolate
          if ((i_min.eq.0).and.(i_max.gt.0)) then
             ! interpolate to zero
             age_min = 0
             i_min   = 0
          elseif  ((i_max.eq.0).and.(i_min.gt.0)) then
             ! extrapolate to higher age
             age_max = age_min
             i_max   = i_min
             age_min = maxval(pop%pop_grid(g)%patch(:)%age(1), &
                  pop%pop_grid(g)%patch(:)%age(1).lt.age_max)
             i_min = maxloc(pop%pop_grid(g)%patch(:)%age(1),1, &
                  pop%pop_grid(g)%patch(:)%age(1).lt.age_max)
          endif

          ! interpolate or extrapolate (growth)
          if ((i_min_growth.eq.0).and.(i_max_growth.gt.0).and.age(iage).le.2) then
             ! interpolate to zero
             age_min_growth = 0
             i_min_growth   = 0
          elseif (((age_min_growth.le.2).or.(i_min_growth.eq.0)).and. &
               (i_max_growth.gt.0).and.age(iage).gt.2) then
             ! extrapolate to lower age
             age_min_growth = age_max_growth
             i_min_growth   = i_max_growth

             age_max_growth = minval(pop%pop_grid(g)%patch(:)%age(1), &
                  pop%pop_grid(g)%patch(:)%age(1).gt.age_min_growth)
             i_max_growth = minloc(pop%pop_grid(g)%patch(:)%age(1), 1, &
                  pop%pop_grid(g)%patch(:)%age(1).gt.age_min_growth)
          elseif  ((i_max_growth.eq.0).and.(i_min_growth.gt.0)) then
             ! extrapolate to higher age
             age_max_growth = age_min_growth
             i_max_growth   = i_min_growth

             age_min_growth = maxval(pop%pop_grid(g)%patch(:)%age(1), &
                  pop%pop_grid(g)%patch(:)%age(1).lt.age_max_growth)
             i_min_growth = maxloc(pop%pop_grid(g)%patch(:)%age(1),1, &
                  pop%pop_grid(g)%patch(:)%age(1).lt.age_max_growth)
          endif

          if (i_min.ne.0.and.age_min.ne.0) then
             MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_min
             where (MASK)
                tmp_array = 1.0_dp
             elsewhere
                tmp_array = 0.0_dp
             endwhere
             tmp_min = sum(pop%pop_grid(g)%patch(:)%layer(1)%biomass,MASK)/sum(tmp_array)
             tmp1_min = sum(pop%pop_grid(g)%patch(:)%stress_mortality,MASK)/sum(tmp_array)
             tmp2_min = sum(pop%pop_grid(g)%patch(:)%crowding_mortality,MASK)/sum(tmp_array)
             tmp3_min = sum(pop%pop_grid(g)%patch(:)%sapwood,MASK)/sum(tmp_array)
             tmp4_min = sum(pop%pop_grid(g)%patch(:)%sapwood_area,MASK)/sum(tmp_array)
          else
             tmp_min  = 0.0_dp
             tmp1_min = 0.0_dp
             tmp2_min = 0.0_dp
             tmp3_min = 0.0_dp
             tmp4_min = 0.0_dp
          endif

          MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_max
          where (MASK)
             tmp_array = 1.0_dp
          elsewhere
             tmp_array = 0.0_dp
          endwhere
          tmp_max  = sum(pop%pop_grid(g)%patch(:)%layer(1)%biomass,MASK)/sum(tmp_array)
          tmp1_max = sum(pop%pop_grid(g)%patch(:)%stress_mortality,MASK)/sum(tmp_array)
          tmp2_max = sum(pop%pop_grid(g)%patch(:)%crowding_mortality,MASK)/sum(tmp_array)
          tmp3_max = sum(pop%pop_grid(g)%patch(:)%sapwood,MASK)/sum(tmp_array)
          tmp4_max = sum(pop%pop_grid(g)%patch(:)%sapwood_area,MASK)/sum(tmp_array)

          cmass_age(iage) = tmp_min + (tmp_max-tmp_min)/real(age_max-age_min,dp)* &
               real(age(iage)-age_min,dp)

          stress_mort_age(iage) = tmp1_min + (tmp1_max-tmp1_min)/real(age_max-age_min,dp)* &
               real(age(iage)-age_min,dp)

          crowd_mort_age(iage) = tmp2_min + (tmp2_max-tmp2_min)/real(age_max-age_min,dp)* &
               real(age(iage)-age_min,dp)

          csapwood_age(iage) = tmp3_min + (tmp3_max-tmp3_min)/real(age_max-age_min,dp)* &
               real(age(iage)-age_min,dp)

          sapwood_area_age(iage) = tmp4_min + (tmp4_max-tmp4_min)/real(age_max-age_min,dp)* &
               real(age(iage)-age_min,dp)

          if (i_min_growth.ne.0.and.age_min_growth.ne.0) then
             MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_min_growth
             where (MASK)
                tmp_array = 1.0_dp
             elsewhere
                tmp_array = 0.0_dp
             endwhere
             tmp_min = sum(pop%pop_grid(g)%patch(:)%growth,MASK)/sum(tmp_array)
          else
             tmp_min = 0.0_dp
          endif

          MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_max_growth
          where (MASK)
             tmp_array = 1.0_dp
          elsewhere
             tmp_array = 0.0_dp
          endwhere
          tmp_max = sum(pop%pop_grid(g)%patch(:)%growth,MASK)/sum(tmp_array)

          growth_age(iage) = tmp_min + (tmp_max-tmp_min)/real(age_max_growth-age_min_growth,dp)* &
               real(age(iage)-age_min_growth,dp)

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

       pop%pop_grid(g)%biomass_age(iage) = cmass_age(iage)
    enddo


    POP%pop_grid(g)%cat_mortality = POP%pop_grid(g)%growth - &
         POP%pop_grid(g)%stress_mortality - &
         POP%pop_grid(g)%crowding_mortality - &
         ( POP%pop_grid(g)%cmass_sum- POP%pop_grid(g)%cmass_sum_old)


    POP%pop_grid(g)%fire_mortality = 0.0_dp
    POP%pop_grid(g)%res_mortality  = 0.0_dp

!!$if (g==4) then
!!$   write(*,*) 'it, nage, growth', it, nage
!!$write(*,*) 'patch biomass', pop%pop_grid(g)%patch(1:5)%layer(1)%biomass
!!$write(*,*) 'patch growth', pop%pop_grid(g)%patch(1:5)%growth
!!$write(*,*) 'stress mort', pop%pop_grid(g)%patch(1:5)%stress_mortality
!!$   write(591, "(350e16.6)") freq_age
!!$   write(601,"(350e16.6)") cmass_age
!!$   write(602,"(350e16.6)") stress_mort_age
!!$   write(603,"(350e16.6)") real(age)
!!$if (it==5) stop
!!$endif
    deallocate(age)
    deallocate(freq_age)
    deallocate(cmass_age)
    deallocate(stress_mort_age)

  end subroutine INTERPOLATE_BIOMASS_1D


  !*******************************************************************************


  subroutine INTERPOLATE_FIREMORTALITY(pop, disturbance_interval,it,g)
    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer(i4b), intent(IN) ::  disturbance_interval(:,:)
    integer(i4b), intent(IN) ::  it,g

    integer(i4b) :: nage,iage, i_min, i_max
    integer(i4b) :: i_min_growth, i_max_growth
    real(dp) :: disturbance_freq,tmp_min, tmp_max, tmp_array(NPATCH2D)
    real(dp) :: tmp5_min, tmp5_max
    logical :: MASK(NPATCH2D)
    integer(i4b) :: age_min, age_max
    integer(i4b) :: age_min_growth, age_max_growth
    integer(i4b), allocatable :: age(:)
    real(dp), allocatable :: fire_mort_age(:)
    real(dp), allocatable :: freq_age(:)

    ! get interpolated fire mortality
    POP%pop_grid(g)%fire_mortality = 0.0_dp
    tmp_array = 0.0_dp
    nage = min(POP%pop_grid(g)%patch(1)%disturbance_interval(1),it)+1 ! maximum age

    if (POP%pop_grid(g)%LU==2) then ! secondary forest
       nage = AGEMAX
       do iage=AGEMAX,1,-1
          if (pop%pop_grid(g)%freq_age(iage)>0) then
             exit
          else
             nage = nage - 1
          endif
       enddo
    endif


    disturbance_freq=1.0_dp/real(disturbance_interval(g,1),dp)
    if(.not.allocated(age)) allocate(age(nage))
    if(.not.allocated(freq_age)) allocate(freq_age(nage))
    if(.not.allocated(fire_mort_age)) allocate(fire_mort_age(nage))

    tmp_min = 0.0_dp
    tmp_max = 0.0_dp
    pop%pop_grid(g)%biomass_age = 0.0_dp

    if (POP%pop_grid(g)%LU==2) then ! secondary forest
       do iage = 1, nage
          age(iage) = iage-1
          freq_age(iage) =  pop%pop_grid(g)%freq_age(iage)
       enddo
    else
       do iage = 1, nage
          age(iage) = iage-1
          freq_age(iage) =  REALExponential(disturbance_freq,real(age(iage),dp))
          pop%pop_grid(g)%freq_age(iage) = freq_age(iage)
       end do

    endif
    if (sum(freq_age)>0.0_dp) freq_age = freq_age/sum(freq_age)

    do iage = 1, nage
       ! get nearest ages bracketing age(iage)
       if (any(pop%pop_grid(g)%patch(:)%age(1).le.age(iage))) then
          age_min = maxval(pop%pop_grid(g)%patch(:)%age(1), &
               pop%pop_grid(g)%patch(:)%age(1).le.age(iage))
          i_min = maxloc(pop%pop_grid(g)%patch(:)%age(1), 1, &
               pop%pop_grid(g)%patch(:)%age(1).le.age(iage))
       else
          age_min = 0
          i_min = 0
       endif
       if (any(pop%pop_grid(g)%patch(:)%age(1).ge.age(iage))) then
          age_max = minval(pop%pop_grid(g)%patch(:)%age(1), &
               pop%pop_grid(g)%patch(:)%age(1).ge.age(iage))
          i_max = minloc(pop%pop_grid(g)%patch(:)%age(1), 1, &
               pop%pop_grid(g)%patch(:)%age(1).ge.age(iage))
       else
          age_max = 0
          i_max = 0
       endif

       age_min_growth = age_min
       age_max_growth = age_max
       i_min_growth = i_min
       i_max_growth = i_max

       if ((i_min.gt.0).and.(i_max.gt.0).and.(age_max.eq.age_min)) then
          ! no need to interpolate

          MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_min
          where (MASK)
             tmp_array = 1.0_dp
          elsewhere
             tmp_array = 0.0_dp
          endwhere
          fire_mort_age(iage)= &
               sum(pop%pop_grid(g)%patch(:)%fire_mortality,MASK)/sum(tmp_array)
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
             age_min = maxval(pop%pop_grid(g)%patch(:)%age(1), &
                  pop%pop_grid(g)%patch(:)%age(1).lt.age_max)
             i_min = maxloc(pop%pop_grid(g)%patch(:)%age(1),1, &
                  pop%pop_grid(g)%patch(:)%age(1).lt.age_max)
          endif


          tmp5_min = 0.0_dp
          if (i_min.ne.0.and.age_min.ne.0) then
             MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_min
             where (MASK)
                tmp_array = 1.0_dp
             elsewhere
                tmp_array = 0.0_dp
             endwhere
             if ( any(MASK) ) then
                tmp5_min = sum(pop%pop_grid(g)%patch(:)%fire_mortality,MASK)/sum(tmp_array)
             endif
          endif

          tmp5_max = 0.0_dp
          MASK = pop%pop_grid(g)%patch(:)%age(1).eq.age_max
          where (MASK)
             tmp_array = 1.0_dp
          elsewhere
             tmp_array = 0.0_dp
          endwhere
          IF any(MASK) then
             tmp5_max = sum(pop%pop_grid(g)%patch(:)%fire_mortality,MASK)/sum(tmp_array)
          endif

          fire_mort_age(iage) = tmp5_min + (tmp5_max-tmp5_min)/real(age_max-age_min,dp)* &
               real(age(iage)-age_min,dp)

       endif
       POP%pop_grid(g)%fire_mortality =  POP%pop_grid(g)%fire_mortality + &
            freq_age(iage)*fire_mort_age(iage)

    enddo

    deallocate(age)
    deallocate(freq_age)
    deallocate(fire_mort_age)

  end subroutine INTERPOLATE_FIREMORTALITY


  !*******************************************************************************


  subroutine ADJUST_POP_FOR_FIRE(pop,disturbance_interval, burned_area, FLI)
    ! reduces biomass on a cohort basis according to mortality vs dbh function
    ! interpolates patch-based fire mortality to get grid-cell mortality
    implicit none

    type( POP_TYPE ), intent(INOUT)  :: pop
    integer(i4b), intent(IN)        ::  disturbance_interval(:,:)
    real(dp),  intent(IN)            :: burned_area(:), FLI(:)
    integer(i4b) :: g, np, c, k, it, nc
    real(dp) :: mort, cmass_stem, dbh


    np = size(POP%POP_grid)
    mort = 0.0

    do g=1,np
       POP%pop_grid(g)%fire_mortality = 0.0_dp



       if (burned_area(g) > 0.0_dp) then

          it = maxval(pop%pop_grid(g)%patch(:)%age(1)) + 1
          do k=1,NPATCH
             nc = pop%pop_grid(g)%patch(k)%Layer(1)%ncohort

             pop%pop_grid(g)%patch(k)%fire_mortality = 0.0_dp
             do c=1,nc

                dbh = pop%pop_grid(g)%patch(k)%layer(1)%cohort(c)%diameter*100.0_dp
                cmass_stem = pop%pop_grid(g)%patch(k)%layer(1)%cohort(c)%biomass

                mort = TopKill_Collins(dbh, FLI(g)) * burned_area(g)

                pop%pop_grid(g)%patch(k)%fire_mortality = mort* &
                     pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%biomass+ &
                     pop%pop_grid(g)%patch(k)%fire_mortality

                pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%biomass = cmass_stem*(1.0_dp-mort)
                pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%heartwood =  &
                     pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%heartwood *(1.0_dp-mort)
                pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%sapwood =  &
                     pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%sapwood *(1.0_dp-mort)

                pop%pop_grid(g)%patch(k)%fire_top_kill_density = &
                     pop%pop_grid(g)%patch(k)%fire_top_kill_density + &
                     pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%density *mort

                pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%density = &
                     pop%pop_grid(g)%patch(k)%Layer(1)%cohort(c)%density*(1.0_dp-mort)


             enddo


             nc = pop%pop_grid(g)%patch(k)%Layer(1)%ncohort
             pop%pop_grid(g)%patch(k)%biomass_old =  pop%pop_grid(g)%patch(k)%Layer(1)%biomass
             pop%pop_grid(g)%patch(k)%Layer(1)%biomass = &
                  sum(pop%pop_grid(g)%patch(k)%Layer(1)%cohort(1:nc)%biomass)

             ! need to remove cohorts with very low density?
             ! This will get done at the end of the year anyway

          enddo


       endif
       ! INTREPOLATE amongst patches to get total biomass lost to fire
       ! creates new value for  POP%pop_grid(g)%fire_mortality
       call INTERPOLATE_FIREMORTALITY(pop, disturbance_interval,it,g)

       !CLN Kill ratio to be used within BLAZE to compute fluxes
       POP%pop_grid(g)%rkill = 0.
       if ( POP%pop_grid(g)%cmass_sum .gt. 0.) then
          POP%pop_grid(g)%rkill = POP%pop_grid(g)%fire_mortality / POP%pop_grid(g)%cmass_sum
       else
          POP%pop_grid(g)%rkill = 0.
       endif

       POP%pop_grid(g)%cmass_sum = POP%pop_grid(g)%cmass_sum - POP%pop_grid(g)%fire_mortality
       


     enddo

   end subroutine ADJUST_POP_FOR_FIRE


   !*******************************************************************************


subroutine INTERPOLATE_BIOMASS_2D(pop, disturbance_interval,it,g)

     use mo_utils, only: eq
#ifdef __MPI__
     use mpi, only: MPI_Abort
#endif

implicit none

type(POP_TYPE), intent(INOUT) :: POP
integer(i4b), intent(IN) ::  disturbance_interval(:,:)
integer(i4b), intent(IN) ::  it,g
integer(i4b), allocatable :: A1(:), A2(:) ! interpolated ages
integer(i4b), allocatable :: xobs(:), yobs(:)     ! observed ages
real(dp), allocatable :: z1obs(:), z2obs(:), z3obs(:) ! observed biomass, stress_mort, crowd_mort
real(dp), allocatable :: z1interp(:), z2interp(:), z3interp(:) ! interpolated biomass, stress mortality, crowding mortality
real(dp), allocatable :: freq_interp(:) ! weightings for interpolated age pairs
real(dp), allocatable :: zp(:)  ! euclidean distance from interpolated age pair to observed age pairs
integer(i4b) :: age_max(2), nrep(NPATCH2D+1)
integer(i4b) :: tmp1, tmp2, I1, I2,I3,I4
integer(i4b) :: x1, x2, x3, x4, y1, y2, y3, y4
integer(i4b) :: p, j, k, n, np, nobs, count_extrap, ct
logical :: flag
real(dp) :: biomass(NPATCH2D+1), stress_mort(NPATCH2D+1), crowd_mort(NPATCH2D+1)
integer(i4b) :: age1(NPATCH2D+1),  age2(NPATCH2D+1)
real(dp) :: zmin
integer(i4b), allocatable :: interp_case(:), tmp_array(:), tmp(:)
real(dp) :: area(4,4), x(4), y(4),  disturbance_freq1,  disturbance_freq2
integer(i4b) :: triangle_points(4,3), I_inside_triangle, Ineighbour(8)
logical ::  MASK_INSIDE_TRIANGLE(4), IS_NEIGHBOUR(8), tmp_logical
logical, allocatable :: MASK2(:), MASK3(:), MASK4(:)
integer(i4b), allocatable :: address(:,:)
#ifdef __MPI__
 integer :: ierr
#endif

 POP%pop_grid(g)%cmass_sum = 0.0_dp
 POP%pop_grid(g)%stress_mortality = 0.0_dp
 POP%pop_grid(g)%crowding_mortality = 0.0_dp

! Construct Age  Interpolating Grid
age_max(1) =  min(POP%pop_grid(g)%patch(1)%disturbance_interval(1),it)+1 ! maximum age
age_max(2) =  min(POP%pop_grid(g)%patch(1)%disturbance_interval(NDISTURB),it)+1 ! maximum age

p=0
do j=0,age_max(1)
   do k=0, age_max(2)
      if (k>j) then
         p=p+1
      endif
   enddo
enddo
np = p

allocate(A1(np))
allocate(A2(np))
allocate(z1interp(np))
allocate(z2interp(np))
allocate(z3interp(np))
allocate(freq_interp(np))
allocate(interp_case(np))
allocate(tmp_array(np))
allocate(tmp(np))
allocate(address(age_max(1)+1,age_max(2)+1))
p=0
address = -9999
do j=0,age_max(1)
   do k=0, age_max(2)
      if (k>j) then
         p=p+1
         A1(p) = j
         A2(p) = k
         address(j+1,k+1) = p
      endif
   enddo
enddo

! Construct Age observations

nrep = 1
flag = .false.
j = 1
! create "observed" age pair (0,0) with zero biomass
stress_mort(j) = 0.0_dp
biomass(j) = 0.0_dp
crowd_mort(j) = 0.0_dp
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
allocate(xobs(nobs))
allocate(yobs(nobs))
allocate(z1obs(nobs))
allocate(z2obs(nobs))
allocate(z3obs(nobs))
allocate(zp(nobs))
allocate(MASK2(nobs))
allocate(MASK3(nobs))
allocate(MASK4(nobs))
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
   disturbance_freq1=1.0_dp/real(disturbance_interval(g,1),dp)
   disturbance_freq2=1.0_dp/real(disturbance_interval(g,2),dp)

   freq_interp(k) = REALExponential(disturbance_freq1,real(A1(k),dp)) * &
        REALExponential(disturbance_freq2,real(A2(k),dp))

enddo
freq_interp = freq_interp/sum(freq_interp)

! interpolate
do p=1,np   ! loop over interpolated age pairs
   ! get distance to all observations
   do j=1,nobs
      zp(j) = sqrt((real(A1(p),dp)-real(xobs(j),dp))**2+(real(A2(p),dp)-real(yobs(j),dp))**2)
   enddo

   ! get closest point
   zmin = minval(zp)
   I1 = minloc(zp,1)
   x1 = xobs(I1)
   y1 = yobs(I1)


   ! check for obs locations forming a quadrangle around interpolating point
   MASK2 = (sign(1,A1(p)-xobs)== -sign(1,A1(p)-x1)).and.(sign(1,A2(p)-yobs)== sign(1,A2(p)-y1).and.A1(p).ne.xobs)
   MASK3 = (sign(1,A1(p)-xobs)== sign(1,A1(p)-x1)).and.(sign(1,A2(p)-yobs)== -sign(1,A2(p)-y1).and.A2(p).ne.yobs)
   MASK4 = (sign(1,A1(p)-xobs)== -sign(1,A1(p)-x1)).and.(sign(1,A2(p)-yobs)== -sign(1,A2(p)-y1) &
        .and.A1(p).ne.xobs.and.A2(p).ne.yobs)

   if ((any(MASK2)).and.(any(MASK3)).and.(any(MASK4)))    then
      ! get nearest point with opposing sign of x displacement
      I2 =  minloc(zp,1, MASK2)
      x2 = xobs(I2)
      y2 = yobs(I2)
      ! get nearest point with opposing sign of y displacement
      I3 =  minloc(zp,1, MASK3)
      x3 = xobs(I3)
      y3 = yobs(I3)
      ! get nearest point with opposing sign of x & y displacements
      I4 =  minloc(zp,1, MASK4)
      x4 = xobs(I4)
      y4 = yobs(I4)

    tmp_logical = .not.( (x2.eq.0.and.y2.eq.0).or.(x3.eq.0.and.y3.eq.0).or.(x4.eq.0.and.y4.eq.0))
   endif


   if ((A1(p)==x1).and.(A2(p)==y1)) then
      interp_case(p)=1
   elseif ((any(MASK2)).and.(any(MASK3)).and.(any(MASK4)) .and.tmp_logical) then ! quadrangle (without (0,0)) exists
      interp_case(p) = 2
   else
      interp_case(p) = 3
   endif!j=1



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
      I2 =  minloc(zp,1, MASK2)
      x2 = xobs(I2)
      y2 = yobs(I2)
      ! get nearest point with opposing sign of y displacement
      I3 =  minloc(zp,1, MASK3)
      x3 = xobs(I3)
      y3 = yobs(I3)
      ! get nearest point with opposing sign of x & y displacements
      I4 =  minloc(zp,1, MASK4)
      x4 = xobs(I4)
      y4 = yobs(I4)


      x=real((/x1, x2, x3, x4/),dp)
      y =real((/y1, y2, y3, y4/),dp)

      ! get area of four possible triangles from 4 vertices, and corresponding 3 partial
      ! triangles, with observation as one vertex
      area = 0.0_dp
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

      MASK_INSIDE_TRIANGLE = eq(area(:,1), sum(area(:,2:4),2))
      I_inside_triangle = minloc(area(:,1),1,MASK_INSIDE_TRIANGLE)
      z1interp(p) =sum( z1obs(triangle_points(I_inside_triangle,:))* &
           area(I_inside_triangle,2:4))/sum(area(I_inside_triangle,2:4))
      z2interp(p) =sum( z2obs(triangle_points(I_inside_triangle,:))* &
           area(I_inside_triangle,2:4))/sum(area(I_inside_triangle,2:4))
      z3interp(p) =sum( z3obs(triangle_points(I_inside_triangle,:))* &
           area(I_inside_triangle,2:4))/sum(area(I_inside_triangle,2:4))


   case default
      write(*,*) " illegal interpolation case."
#ifdef __MPI__
      call MPI_Abort(0, 86, ierr) ! Do not know comm nor rank here
#else
      stop 86
#endif
   end select ! interpolation case

enddo ! loop over interpolated age pairs

p = 0 ! counter for while loop
count_extrap = 0
! Extrapolation

if (any(interp_case.ne.3)) then
   do while (any(interp_case==3))
      do ct = 3,3
         do p=1,np
            ! find number of neighbours for each extrapolable point
            is_neighbour = .false.
            Ineighbour = 1
            tmp(p) = 0
            if (interp_case(p)==3) then

               do k=1,8
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
                  if(tmp1.ge.0.and.tmp1.le.age_max(1).and.tmp2.ge.0.and.tmp2.le.age_max(2).and.tmp2.gt.tmp1) then
                     is_neighbour(k) = (interp_case(address(tmp1+1,tmp2+1)).ne.3)
                     Ineighbour(k) = address(tmp1+1,tmp2+1)
                  endif

               enddo


               tmp(p) = count(is_neighbour)
               if (tmp(p).ge.ct) then
                  z1interp(p) = sum(z1interp(Ineighbour) ,is_neighbour)/real(tmp(p),dp)
                  z2interp(p) = sum(z2interp(Ineighbour) ,is_neighbour)/real(tmp(p),dp)
                  z3interp(p) = sum(z3interp(Ineighbour) ,is_neighbour)/real(tmp(p),dp)
                  interp_case(p) = 2
               endif
            endif
         enddo !1, np
      enddo
      if (all(tmp.lt.3)) then
         ! maximum of one neighbour
         p = maxloc(tmp,1)
         is_neighbour = .false.
         Ineighbour = 1
         tmp(p) = 0
         do k=1,8
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
            if(tmp1.ge.0.and.tmp1.le.age_max(1).and.tmp2.ge.0.and.tmp2.le.age_max(2).and.tmp2.gt.tmp1) then
               is_neighbour(k) = (interp_case(address(tmp1+1,tmp2+1)).ne.3)
               Ineighbour(k) = address(tmp1+1,tmp2+1)
            endif

         enddo

         tmp(p) = count(is_neighbour)
         z1interp(p) = sum(z1interp(Ineighbour) ,is_neighbour)/real(tmp(p),dp)
         z2interp(p) = sum(z2interp(Ineighbour) ,is_neighbour)/real(tmp(p),dp)
         z3interp(p) = sum(z3interp(Ineighbour) ,is_neighbour)/real(tmp(p),dp)
         interp_case(p) = 2

      endif  ! ALL(tmp.lt.2)

   enddo ! do while
endif


do p = 1,np
   POP%pop_grid(g)%cmass_sum =  POP%pop_grid(g)%cmass_sum + &
        freq_interp(p)*z1interp(p)
   POP%pop_grid(g)%stress_mortality =  POP%pop_grid(g)%stress_mortality + &
        freq_interp(p)*z2interp(p)
   POP%pop_grid(g)%crowding_mortality =  POP%pop_grid(g)%crowding_mortality + &
        freq_interp(p)*z3interp(p)


enddo

POP%pop_grid(g)%cat_mortality = POP%pop_grid(g)%cmass_sum *  disturbance_freq2
POP%pop_grid(g)%fire_mortality = POP%pop_grid(g)%growth - &
     POP%pop_grid(g)%cat_mortality - &
     POP%pop_grid(g)%stress_mortality - &
     POP%pop_grid(g)%crowding_mortality - &
     ( POP%pop_grid(g)%cmass_sum- POP%pop_grid(g)%cmass_sum_old)


deallocate(xobs)
deallocate(yobs)
deallocate(z1obs)
deallocate(z2obs)
deallocate(z3obs)
deallocate(zp)
deallocate(MASK2)
deallocate(MASK3)
deallocate(MASK4)
deallocate(A1)
deallocate(A2)
deallocate(z1interp)
deallocate(z2interp)
deallocate(z3interp)
deallocate(freq_interp)
deallocate(interp_case)
deallocate(tmp)
deallocate(tmp_array)
deallocate(address)

end subroutine INTERPOLATE_BIOMASS_2D


!******************************************************************************


subroutine SMOOTH_FLUX(POP,g,t)

  implicit none

  type(POP_TYPE), intent(INOUT) :: POP
  integer(i4b),   intent(IN)    :: g, t

  integer(i4b), parameter :: SPAN = NYEAR_WINDOW
  real(dp) :: x(NYEAR_SMOOTH), y(NYEAR_SMOOTH), a, b, r
  real(dp) :: sumflux, sumsmooth, flux(NYEAR_HISTORY), smoothed_flux
  real(dp) :: dbuf
  integer(i4b) :: t0, tt, n, k

  ! update fire_mortality_history

  if (t.gt.NYEAR_HISTORY) then
     do k = 1, NYEAR_HISTORY-1
        POP%pop_grid(g)%fire_mortality_history(k) = POP%pop_grid(g)%fire_mortality_history(k+1)
     enddo
  endif
  POP%pop_grid(g)%fire_mortality_history(t) = POP%pop_grid(g)%fire_mortality

  flux = POP%pop_grid(g)%fire_mortality_history
  t0 = t-SPAN
  n = 0
  sumflux   = 0.0_dp
  sumsmooth = 0.0_dp
  do tt=1, NYEAR_SMOOTH
     if ((t0+tt).ge.1 .and. (t0+tt).le.t+1) then
        sumflux = sumflux + flux(t0+tt)
        y(tt) = flux(t0+tt)
        x(tt) = real(tt,dp)
        n = n+1
        if ((t0+tt).eq.t+1) then
           call regress(x,y,n,a,b,r)
        endif
     else
        sumflux = sumflux + (a + b*tt)
        n = n+ 1
     endif
  enddo

  dbuf = POP%pop_grid(g)%smoothing_buffer / (real(NYEAR_SMOOTH,dp)/2.0_dp)
  smoothed_flux = max(sumflux/real(n)+dbuf, 0.0_dp)
  POP%pop_grid(g)%smoothing_buffer = POP%pop_grid(g)%smoothing_buffer + flux(t) - smoothed_flux
  POP%pop_grid(g)%fire_mortality_smoothed = smoothed_flux

end subroutine SMOOTH_FLUX


!******************************************************************************


subroutine SMOOTH_FLUX_cat(POP,g,t)

  implicit none

  type(POP_TYPE), intent(INOUT) :: POP
  integer(i4b), intent(IN) :: g, t
  integer(i4b), parameter :: SPAN = NYEAR_WINDOW
  real(dp) :: x(NYEAR_SMOOTH), y(NYEAR_SMOOTH), a, b, r
  real(dp) :: sumflux, sumsmooth, flux(NYEAR_HISTORY), smoothed_flux
  real(dp) :: dbuf
  integer(i4b) :: t0, tt, n, k

  ! update cat_mortality_history

  if (t.gt.NYEAR_HISTORY) then
     do k = 1, NYEAR_HISTORY-1
        POP%pop_grid(g)%cat_mortality_history(k) = POP%pop_grid(g)%cat_mortality_history(k+1)
     enddo
  endif
  POP%pop_grid(g)%cat_mortality_history(t) = POP%pop_grid(g)%cat_mortality

  flux = POP%pop_grid(g)%cat_mortality_history
  t0 = t-SPAN
  n = 0
  sumflux = 0.0_dp
  sumsmooth = 0.0_dp
  do tt=1, NYEAR_SMOOTH
     if ((t0+tt).ge.1 .and. (t0+tt).le.t+1) then
        sumflux = sumflux + flux(t0+tt-1)
        y(tt) = flux(t0+tt-1)
        x(tt) = real(tt, dp)
        n = n+1
        if ((t0+tt).eq.t+1) then
           call regress(x,y,n,a,b,r)
        endif
     else
        sumflux = sumflux + (a + b*tt)
        n = n+ 1
     endif
  enddo

  dbuf = POP%pop_grid(g)%smoothing_buffer_cat/(real(NYEAR_SMOOTH,dp)/2.0_dp)
  smoothed_flux = max(sumflux/real(n)+dbuf, 0.0_dp)
  POP%pop_grid(g)%smoothing_buffer_cat = POP%pop_grid(g)%smoothing_buffer_cat + flux(t) - smoothed_flux
  POP%pop_grid(g)%cat_mortality_smoothed = smoothed_flux

end subroutine SMOOTH_FLUX_cat


!******************************************************************************


subroutine REGRESS(x, y, n, a, b, r)

  implicit none

  real(dp), intent(IN) :: x(:), y(:)
  real(dp), intent(OUT) :: a, b, r
  integer(i4b), intent(IN) :: n
  real(dp)::  sx,sy,sxx,sxy,delta,meanx,meany,sdx,sdy
  integer(i4b) :: i

  ! Performs a linear regression of array y on array x (n values)
  ! returning parameters a and b in the fitted model: y=a+bx
  ! Source: Press et al 1986, Sect 14.2
  ! also returns Pearson r

  sx=0.0_dp
  sy=0.0_dp
  sxx=0.0_dp
  sxy=0.0_dp
  do i=1,n
     sx = sx + x(i)
     sy = sy + y(i)
     sxx = sxx + x(i) * y(i)
     sxy = sxy + x(i)*y(i)
  enddo
  delta = real(n,dp)*sxx - sx*sx
  a=(sxx*sy-sx*sxy)/delta
  b=(real(n,dp)*sxy-sx*sy)/delta
  meanx=sx/real(n,dp)
  meany=sy/real(n,dp)
  sdx = 0.0_dp
  sdy = 0.0_dp
  do i=1,n
     sdx = sdx + (x(i)-meanx)*(x(i)-meanx)
     sdy = sdy + (y(i)-meany)*(y(i)-meany)
  enddo
  sdx=sqrt(sdx/real(n-1,dp))
  sdy=sqrt(sdy/real(n-1,dp))
  if ((abs(sdx) .lt. tiny(0.0_dp)) .or. (abs(sdy) .lt. tiny(0.0_dp))) then
     r = 0.0_dp
  else
     r = b*sdx/sdy
  endif

end subroutine REGRESS


!******************************************************************************


real(dp) function Area_Triangle(x1,y1,x2,y2,x3,y3)

  implicit none

  real(dp), intent(IN) :: x1, y1, x2, y2, x3, y3

  area_triangle = abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0_dp)

end function Area_Triangle


!******************************************************************************


! Fraction of topkill by DBH , according to Fig. 2 of Collins, J. Ec., 2020
real(dp) function TopKill_Collins(dbh, FLI)

  implicit none

  real(dp), intent(IN) :: dbh, FLI
  real(dp), parameter:: a = 0.79130437_dp, b =  0.07999477_dp, c = 0.06326282_dp

  if ((FLI > 3000.0_dp).and.(FLI < 7000.0_dp)) then
     TopKill_Collins = 0.08_dp
  elseif (FLI > 7000.0_dp) then
     TopKill_Collins = max(a * exp(-b * dbh) + c, 0.08_dp)
  else ! low intensity fires
     TopKill_Collins = 0.0_dp
  endif

end function TopKill_Collins


!******************************************************************************


! Top-End Allometry
! Computes tree stem diameter (m) and basal area (m2/ha)
! given height (m), stem biomass (kgC/m2) and tree population density (indiv/m2)

subroutine Allometry(height,biomass,density,diam,basal)

  implicit none

  real(dp), intent(IN) :: height
  real(dp), intent(IN) :: biomass
  real(dp), intent(IN) :: density
  real(dp), intent(OUT) :: diam
  real(dp), intent(OUT) :: basal

  real(dp) :: delta,rh

  delta=2.0_dp*sqrt(biomass/density/WD/PI)
  rh=1.0_dp/sqrt(height)

  diam=delta*rh
  basal=PI*(diam/2.0_dp)*(diam/2.0_dp)*density*1.0e4_dp

end subroutine Allometry


!*******************************************************************************


  subroutine Williams_Allometry(agBiomass, density, height, dbh, basal)

    implicit none
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
    real(dp), intent(IN) :: agbiomass, density
    real(dp), intent(OUT):: height, dbh, basal

    real(dp), parameter  :: beta0 = -2.3046_dp
    real(dp), parameter  :: beta1 = 2.5243_dp
    real(dp), parameter  :: gC2DM = 1.0_dp/0.49_dp  ! ratio Dry matter mass to g(C)

    ! Compute dbh using model 5b and converting from cm -> m
    dbh    = ( agbiomass * gC2DM / density / exp(beta0) ) ** (1.0_dp/beta1)
    dbh    = dbh * 0.01_dp
    ! Compute basal area [m^2/m^2]
    basal  = PI * (0.5_dp * dbh)**2 * density
    ! Compute Height using cylindrical approach
    height = agBiomass / ( WD * basal )
    ! Basal Area [m^2/m^2]->[m^2/ha]
    basal  = basal * 1.0e4_dp

  end subroutine Williams_Allometry


  !*******************************************************************************


  subroutine POP_init(POP, disturbance_interval, np, Iwood, precip)

    use POP_types, only: POP_TYPE
    use TypeDef,   only: i4b

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer(i4b),   intent(IN)    :: disturbance_interval(:,:)
    integer(i4b),   intent(IN)    :: np
    integer(i4b),   intent(IN)    :: Iwood(:)
    real(dp),       intent(IN), optional :: precip(:)

    integer(i4b) :: j, k

    call alloc_POP(pop,int(np))
    POP%np     = np
    POP%Iwood  = Iwood
    POP%it_pop = 0
    ! POP%LU = 1  ! initialise to primary forest
    POP%pop_grid(:)%LU = 1

    call ZeroPOP(pop)

    call InitPOP2D_Poisson(pop, int(disturbance_interval,i4b))

    do j=1,np
       do k=1,NPATCH2D
          ! understorey recruitment
          if (present(precip)) then
             call layer_recruitment_single_patch(pop, k, j, precip)
          else
             call layer_recruitment_single_patch(pop, k, j)
          endif
       enddo
    enddo

  end subroutine POP_init


  !*******************************************************************************


  subroutine POP_init_single(POP, disturbance_interval, n, precip)

    use POP_types, only: POP_TYPE
    use TypeDef,   only: i4b

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer(i4b),   intent(IN)    :: disturbance_interval(:,:)
    integer(i4b),   intent(IN)    :: n
    real(dp),       intent(IN), optional :: precip(:)

    integer(i4b) :: j, k

    POP%it_pop(n) = 0
    call ZeroPOP(pop, n)

    call InitPOP2D_Poisson(pop, int(disturbance_interval,i4b), n)

    do j=n,n
       do k=1,NPATCH2D
          ! understorey recruitment
          if (present(precip)) then
             call layer_recruitment_single_patch(pop, k, j, precip)
          else
             call layer_recruitment_single_patch(pop, k, j)
          endif
       enddo
    enddo

  end subroutine POP_init_single


  !*******************************************************************************


  subroutine alloc_POP(POP, arraysize)

    use POP_Types, only: POP_TYPE

    implicit none

    type(POP_TYPE), intent(INOUT) :: POP
    integer,        intent(IN)    :: arraysize

    if (.not. allocated(POP%POP_Grid)) allocate(POP%POP_Grid(arraysize))
    if (.not. allocated(POP%Iwood))    allocate(POP%Iwood(arraysize))
    ! IF (.NOT. ALLOCATED(POP%LU)) ALLOCATE (POP%LU(arraysize))
    if (.not. allocated(POP%it_pop))   allocate(POP%it_pop(arraysize))

  end subroutine alloc_POP

  !*******************************************************************************

end module POPModule

!*******************************************************************************
