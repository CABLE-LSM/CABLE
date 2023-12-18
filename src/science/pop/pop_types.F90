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
     INTEGER , DIMENSION(:), ALLOCATABLE    :: it_pop
     INTEGER :: np
     INTEGER, DIMENSION(:), ALLOCATABLE :: Iwood ! , LU
  END TYPE POP_TYPE

END MODULE POP_Types

