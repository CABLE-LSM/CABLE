!> \file mo_isotope.f90

!> \brief Isotope parameters and helper routines

!> \details Provides isotope parameters such as international standards
!> and convenience routines such as delta value calculations.

!> \author Matthias Cuntz
!> \date Apr 2019

! License
! -------
! This file is part of the JAMS Fortran package, distributed under the MIT License.
!
! Copyright (c) 2019 Matthias Cuntz - mc (at) macu (dot) de
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

MODULE mo_isotope

  use cable_def_types_mod, only: r2
  use mo_constants, only: Mw_r2, Ma_r2, twothird_r2

  implicit none

  private

  ! Public routines
  public :: delta     ! isotopic delta value
  public :: delta1000 ! isotopic delta value in permil
  public :: isoratio  ! ratio of rare to abundant isotope
  public :: isosanity ! set rare and abundant to 0 if any < precision

  ! Public parameters
  ! ratio diffusion of air vs. water vapour ~ 1.6
  real(r2), parameter, public :: ratio_diff_air2vap = Ma_r2 / Mw_r2
  ! ratio of "diffusion" of air vs. water vapour through boundary layer ~ 1.4
  real(r2), parameter, public :: ratio_boundary_air2vap = ratio_diff_air2vap**twothird_r2
  ! ratio diffusion of water vapour vs. air ~ 0.62
  real(r2), parameter, public :: ratio_diff_vap2air = Mw_r2 / Ma_r2
  ! ratio of "diffusion" of water vapour vs. air through boundary layer ~ 0.73
  real(r2), parameter, public :: ratio_boundary_vap2air = ratio_diff_vap2air**twothird_r2
  ! 13C/12C VPDB standard
  real(r2), parameter, public :: vpdbc13 = 0.0112372_r2
  ! molecular weight of 13C (kg mol-1)
  real(r2), parameter, public :: Mc13   = 13.e-3_r2
  ! molecular weight of 13CO2 (kg mol-1)
  real(r2), parameter, public :: Mc13o2 = 45.e-3_r2

contains

  ! ------------------------------------------------------------------

  !     NAME
  !         delta

  !     PURPOSE
  !>        \brief Isotopic delta value.

  !>        \details Elemental routine to calculate the delta value
  !>        from an isotope ratio or the rare and the abundant isotope concentrations
  !>        relative to a given standard.

  !>        The delta value can be set to a default value if the abundant
  !>        isotopic concentration is below a give precision.

  !     CALLING SEQUENCE
  !         del = delta(rare, abundant, standard, default, precision)

  !     PARAMETER
  !>        \param[in] "real(r2) :: rare"                   Rare isotope concentration or isotope ratio
  !>        \param[in] "real(r2), optional :: abundant"     Abundant isotope concentration if rare is rare isotope concentration
  !>                                                        Default: not used.
  !>        \param[in] "real(r2), optional :: standard"     Isotop ratio of the standard material.
  !>                                                        Default: 1.
  !>        \param[in] "real(r2), optional :: default"      Set isoratio default if abundant < precision.
  !>                                                        Default: 0.
  !>        \param[in] "real(r2), optional :: precision"    Set isoratio default if abundant < precision.

  !     RETURN
  !>       \return      real(r2) :: delta &mdash; Isotopic delta value: rare / abundant / standard - 1.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  elemental pure function delta(rare, abundant, standard, default, precision)

    implicit none

    real(r2), intent(in)           :: rare
    real(r2), intent(in), optional :: abundant
    real(r2), intent(in), optional :: standard
    real(r2), intent(in), optional :: default
    real(r2), intent(in), optional :: precision
    real(r2)                       :: delta

    real(r2) :: istandard, idefault, iprecision

    ! default optionals
    if (present(standard)) then
       istandard = standard
    else
       istandard = 0.0_r2
    endif
    if (present(default)) then
       idefault = default
    else
       idefault = 0.0_r2
    endif
    if (present(precision)) then
       iprecision = precision
    else
       iprecision = 0.0_r2
    endif

    if (present(abundant)) then
       if (abs(abundant) > iprecision) then
          delta = rare / abundant / istandard - 1.0_r2
       else
          delta = idefault
       end if
    else
       delta = rare / istandard - 1.0_r2
    end if
    !MCTest
    ! if (abs(delta) < 100.0_r2*epsilon(1.0_r2)) delta = 0.0_r2
    !MCTest

  end function delta


  ! ------------------------------------------------------------------

  !     NAME
  !         delta1000

  !     PURPOSE
  !>        \brief Isotopic delta value in permil.

  !>        \details Elemental routine to calculate the delta value in permil
  !>        from an isotope ratio or the rare and the abundant isotope concentrations
  !>        relative to a given standard.

  !>        The delta value can be set to a default value if the abundant
  !>        isotopic concentration is below a give precision.

  !     CALLING SEQUENCE
  !         del = delta1000(rare, abundant, standard, default, precision)

  !     PARAMETER
  !>        \param[in] "real(r2) :: rare"                   Rare isotope concentration or isotope ratio
  !>        \param[in] "real(r2), optional :: abundant"     Abundant isotope concentration if rare is rare isotope concentration
  !>                                                        Default: not used.
  !>        \param[in] "real(r2), optional :: standard"     Isotop ratio of the standard material.
  !>                                                        Default: 1.
  !>        \param[in] "real(r2), optional :: default"      Set isoratio default if abundant < precision.
  !>                                                        Default: 0.
  !>        \param[in] "real(r2), optional :: precision"    Set isoratio default if abundant < precision.

  !     RETURN
  !>       \return      real(r2) :: delta &mdash; Isotopic delta value in permil: (rare / abundant / standard - 1.) * 1000.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  elemental pure function delta1000(rare, abundant, standard, default, precision)

    implicit none

    real(r2), intent(in)           :: rare
    real(r2), intent(in), optional :: abundant
    real(r2), intent(in), optional :: standard
    real(r2), intent(in), optional :: default
    real(r2), intent(in), optional :: precision
    real(r2)                       :: delta1000

    real(r2) :: istandard, idefault, iprecision

    ! default optionals
    if (present(standard)) then
       istandard = standard
    else
       istandard = 0.0_r2
    endif
    if (present(default)) then
       idefault = default
    else
       idefault = 0.0_r2
    endif
    if (present(precision)) then
       iprecision = precision
    else
       iprecision = 0.0_r2
    endif

    if (present(abundant)) then
       if (abs(abundant) > iprecision) then
          delta1000 = (rare / abundant / istandard - 1.0_r2) * 1000.0_r2
          ! !MCTest
          ! if (abs(delta1000) < 1000.0_r2*epsilon(1.0_r2)*1000.0_r2*(10._r2**abs(log10(abundant)))) delta1000 = 0.0_r2
          ! !MCTest
       else
          delta1000 = idefault
       end if
    else
       delta1000 = (rare / istandard - 1.0_r2) * 1000.0_r2
       ! !MCTest
       ! if (abs(delta1000) < 1000.0_r2*epsilon(1.0_r2)*1000.0_r2*(10._r2**abs(log10(rare)))) delta1000 = 0.0_r2
       ! !MCTest
    end if

  end function delta1000


  ! ------------------------------------------------------------------

  !     NAME
  !         isoratio

  !     PURPOSE
  !>        \brief Ratio of rare to abundant isotopes.

  !>        \details Elemental routine to calculate the ratio of the rare
  !>        and the abundant isotope concentrations.

  !>        The ratio can be set to a default value if the abundant
  !>        isotopic concentration is below a give precision.

  !     CALLING SEQUENCE
  !         ratio = isoratio(rare, abundant, default, precision)

  !     PARAMETER
  !>        \param[in] "real(r2) :: rare"                   Rare isotope concentration
  !>        \param[in] "real(r2) :: abundant"               Abundant isotope concentration
  !>        \param[in] "real(r2), optional :: default"      Set isoratio default if abundant < precision.
  !>                                                        Default: 0.
  !>        \param[in] "real(r2), optional :: precision"    Set isoratio default if abundant < precision.

  !     RETURN
  !>       \return      real(r2) :: isoratio &mdash; Isotope ratio: rare / abundant

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  elemental pure function isoratio(rare, abundant, default, precision)

    implicit none

    real(r2), intent(in)           :: rare
    real(r2), intent(in)           :: abundant
    real(r2), intent(in), optional :: default
    real(r2), intent(in), optional :: precision
    real(r2)                       :: isoratio

    real(r2) :: idefault, iprecision

    ! default optionals
    if (present(default)) then
       idefault = default
    else
       idefault = 0.0_r2
    endif
    if (present(precision)) then
       iprecision = precision
    else
       iprecision = 0.0_r2
    endif

    if (abs(abundant) > iprecision) then
       isoratio = rare / abundant
    else
       isoratio = idefault
    end if

  end function isoratio


  ! ------------------------------------------------------------------

  !     NAME
  !         isosanity

  !     PURPOSE
  !>        \brief Set rare and abundant isotope composition to zero
  !>        if any of the two is less than epsilon.

  !>        \details If either the rare or the abundant concentration is less than
  !>        a given precision, both concentrations are set to zero.

  !>        The default precision is epsilon of the kind.

  !>        The test can be done on the absolute of the concentrations if negative
  !>        compositions are possible.

  !     CALLING SEQUENCE
  !         call isosanity(rare, abundant, precision, absolute)

  !     PARAMETER
  !>        \param[inout] "real(r2) :: rare"                   Rare isotope concentration
  !>        \param[inout] "real(r2) :: abundant"               Abundant isotope concentration
  !>        \param[in] "real(r2), optional :: precision"    Set rare and abundant to 0 if any of the two is < precision.
  !>                                                        Default: epsilon(1.0_r2)
  !>        \param[in] "logical, optional :: absolute"      If .true.: check abs(rare) and abs(abundant) < precision
  !>                                                        Default: .false.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date May 2020
  elemental pure subroutine isosanity(rare, abundant, precision, absolute)

    implicit none

    real(r2), intent(inout)        :: rare
    real(r2), intent(inout)        :: abundant
    real(r2), intent(in), optional :: precision
    logical,  intent(in), optional :: absolute

    real(r2) :: iprecision
    logical  :: iabsolute

    ! default optionals
    if (present(precision)) then
       iprecision = precision
    else
       iprecision = epsilon(1.0_r2)
    endif
    if (present(absolute)) then
       iabsolute = absolute
    else
       iabsolute = .false.
    endif

    if (iabsolute) then
       if ((abs(rare) < iprecision) .or. (abs(abundant) < iprecision)) then
          rare     = 0.0_r2
          abundant = 0.0_r2
       endif
    else
       if ((rare < iprecision) .or. (abundant < iprecision)) then
          rare     = 0.0_r2
          abundant = 0.0_r2
       endif
    endif

  end subroutine isosanity

END MODULE mo_isotope
