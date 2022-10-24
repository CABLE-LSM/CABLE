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
! Purpose: defines/allocates variables for CASA-CNP
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed for offline CASA-CNP, code revision likely to better
!          suit ACCESS and to merge more consistently with CABLE code
!
!
! ==============================================================================
! casa_variable.f90
!
! the following modules are used when "casacnp" is coupled to "cable"
!   casadimension
!   casaparm
!   casavariable with subroutine alloc_casavariable
!   phenvariable with subroutine alloc_phenvariable

MODULE casaparm
  USE casadimension

  IMPLICIT NONE
  INTEGER, PARAMETER :: initcasa= 1   ! =0 spin; 1 restart file
  INTEGER, PARAMETER :: iceland  = 17 !=13 for casa vegtype =15 for IGBP vegtype
  INTEGER, PARAMETER :: cropland = 9  ! 12 and 14 for IGBP vegtype
  INTEGER, PARAMETER :: croplnd2 =10  ! ditto
  INTEGER, PARAMETER :: forest  = 3
  INTEGER, PARAMETER :: shrub   = 2
  INTEGER, PARAMETER :: grass   = 1
  INTEGER, PARAMETER :: icewater= 0
  INTEGER, PARAMETER :: LEAF    = 1
  INTEGER, PARAMETER :: WOOD    = 2
  INTEGER, PARAMETER :: FROOT   = 3
  !  INTEGER, PARAMETER :: LABILE  = 4
  INTEGER, PARAMETER :: METB    = 1
  INTEGER, PARAMETER :: STR     = 2
  INTEGER, PARAMETER :: CWD     = 3
  INTEGER, PARAMETER :: MIC     = 1
  INTEGER, PARAMETER :: SLOW    = 2
  INTEGER, PARAMETER :: PASS    = 3
  INTEGER, PARAMETER :: PLAB    = 1
  INTEGER, PARAMETER :: PSORB   = 2
  INTEGER, PARAMETER :: POCC    = 3
  !! vh_js !! LALLOC moved to bgcdriver to allow for value to be switchable
  !INTEGER, PARAMETER :: LALLOC  = 0      !=0 constant; 1 variable
  REAL(r_2), PARAMETER :: z30=0.3
  REAL(r_2), PARAMETER :: R0=0.3
  REAL(r_2), PARAMETER :: S0=0.3
  REAL(r_2), PARAMETER :: fixed_stem=1.0/3.0
  REAL(r_2), PARAMETER :: Q10alloc=2.0
  REAL(r_2), PARAMETER :: ratioNCstrfix = 1.0/150.0
  REAL(r_2), PARAMETER :: ratioNPstrfix = 25.0
  REAL(r_2), PARAMETER :: fracCbiomass = 0.50
  REAL(r_2), PARAMETER :: tsoilrefc=25.0
  REAL(r_2), PARAMETER :: tkzeroc=273.15
  REAL(r_2), PARAMETER :: frootparma = 0.3192
  REAL(r_2), PARAMETER :: frootparmb =-0.0485
  REAL(r_2), PARAMETER :: frootparmc = 0.1755
  REAL(r_2), PARAMETER :: xweightalloc = 0.2
  !  REAL(r_2), PARAMETER :: xkplab=0.5*deltcasa
  !  REAL(r_2), PARAMETER :: xkpsorb=0.01*deltcasa
  !  REAL(r_2), PARAMETER :: xkpocc =0.01*deltcasa
END MODULE casaparm
