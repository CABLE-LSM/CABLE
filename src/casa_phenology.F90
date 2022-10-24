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

MODULE phenvariable
  USE casadimension
  IMPLICIT NONE
  TYPE phen_variable
     INTEGER,   DIMENSION(:),  POINTER :: phase
     REAL(r_2), DIMENSION(:),  POINTER :: TKshed
     INTEGER,   DIMENSION(:,:),POINTER :: doyphase
     REAL, DIMENSION(:),  POINTER :: phen   ! fraction of max LAI
     REAL, DIMENSION(:),  POINTER :: aphen  ! annual leaf on sum
     INTEGER,   DIMENSION(:,:),POINTER :: phasespin
     INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_1
     INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_2
     INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_3
     INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_4

  END TYPE phen_variable

CONTAINS

  SUBROUTINE alloc_phenvariable(phen,arraysize)

    IMPLICIT NONE
    TYPE(phen_variable), INTENT(INOUT) :: phen
    INTEGER,             INTENT(IN) :: arraysize

    ALLOCATE(phen%Tkshed(mvtype))
    ALLOCATE(phen%phase(arraysize),         &
         phen%doyphase(arraysize,mphase))
    ALLOCATE(phen%phen(arraysize), &
         phen%aphen(arraysize), &
         phen%phasespin(arraysize,mdyear), &
         phen%doyphasespin_1(arraysize,mdyear), &
         phen%doyphasespin_2(arraysize,mdyear), &
         phen%doyphasespin_3(arraysize,mdyear), &
         phen%doyphasespin_4(arraysize,mdyear))
  END SUBROUTINE alloc_phenvariable

END MODULE phenvariable
