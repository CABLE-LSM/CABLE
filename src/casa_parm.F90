!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: defines/allocates variables for CASA-CNP
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed for offline CASA-CNP, code revision likely to better
!          suit ACCESS and to merge more consistently with CABLE code
!
!
! ==============================================================================
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
  INTEGER, PARAMETER :: LALLOC  = 0      !=0 constant; 1 variable
  REAL(r_2), PARAMETER :: z30=0.3
  REAL(r_2), PARAMETER :: R0=0.3
  REAL(r_2), PARAMETER :: S0=0.3
  REAL(r_2), PARAMETER :: fixed_stem=1.0/3.0
  REAL(r_2), PARAMETER :: Q10alloc=2.0
  REAL(r_2), PARAMETER :: ratioNCstrfix = 1.0/150.0
  REAL(r_2), PARAMETER :: ratioPCstrfix = ratioNCstrfix/25.0
  REAL(r_2), PARAMETER :: fracCbiomass = 0.50
  REAL(r_2), PARAMETER :: tsoilrefc=25.0
  REAL(r_2), PARAMETER :: tkzeroc=273.15
  REAL(r_2), PARAMETER :: frootparma = 0.3192
  REAL(r_2), PARAMETER :: frootparmb =-0.0485
  REAL(r_2), PARAMETER :: frootparmc = 0.1755
  REAL(r_2), PARAMETER :: xweightalloc = 0.2
  REAL(r_2), PARAMETER :: xkplab=0.5*deltcasa
  REAL(r_2), PARAMETER :: xkpsorb=0.01*deltcasa
  REAL(r_2), PARAMETER :: xkpocc =0.01*deltcasa
END MODULE casaparm

