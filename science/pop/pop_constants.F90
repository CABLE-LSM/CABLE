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
!===============================================================================

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
