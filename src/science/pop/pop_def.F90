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

