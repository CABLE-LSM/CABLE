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
!
! Purpose: Solve quadratic and cubic equations and its derivatives
!
! History:
!    * Copy and combine routines fAn_c3, etc. from cable_canopy
!      and cable_optimiseJVratio, Nov 2025, Matthias Cuntz
!
module cable_polynomial

  implicit none

  public :: solve_quadratic, solve_cubic
  public :: deriv_solve_quadratic, deriv_solve_cubic

  private

  ! ieee_quiet_nan or ieee_signaling_nan
  character(len=*), parameter :: nan = 'ieee_quiet_nan'

contains

  !> \brief Solution of quadratic equation
  !
  !> \details
  !!    The solution of the quadratic equation
  !!    \f[ a x^2 + b x + c = 0 \f]
  !!    is
  !!    \f[ x_{1/2} = (-b \pm \sqrt(b^2 - 4ac)) / (2a)2 \f]
  !!
  !!    The sign can be chosen. The default uses the negative sign.
  !
  !> \param[in] "real(r2) :: a, b, c"         The coefficients of the quadratic equation
  !> \param[in] "character, optional :: sgn"  '+' or '-' for sign in solution (default: '-')
  !> \retval "real(r2) :: x"                   The solution of the quadratic equation
  !
  !> \author Matthias Cuntz
  !> \date Nov 2025
  elemental pure function solve_quadratic(a, b, c, sgn)

    use cable_def_types_mod, only: r2
    use mo_utils, only: special_value

    implicit none

    real(r2),  intent(in) :: a, b, c
    character, intent(in), optional :: sgn
    real(r2) :: solve_quadratic

    real(r2) :: s2
    character :: isgn

    if (present(sgn)) then
       isgn = sgn
    else
       isgn = '-'
    endif
    
    if (abs(c) > tiny(1.0_r2)) then
       if (abs(a) > tiny(1.0_r2)) then
          s2 = b**2 - 4.0_r2 * a * c
          if (s2 > tiny(1.0_r2)) then
             if (isgn == '-') then
                solve_quadratic = (-b - sqrt(s2)) / (2.0_r2 * a)
             else
                solve_quadratic = (-b + sqrt(s2)) / (2.0_r2 * a)
             endif
          else
             solve_quadratic = special_value(c, nan)
          endif
       else
          if (abs(b) > tiny(1.0_r2)) then
             solve_quadratic = -c / b
          else
             solve_quadratic = special_value(c, nan)
          endif
       endif
    else
       if (abs(a) > tiny(1.0_r2)) then
          solve_quadratic = -b / a
       else
          solve_quadratic = 0.0_r2
       endif
    endif

  end function solve_quadratic

  
  !> \brief Derivative of the solution of quadratic equation
  !
  !> \details
  !!    The solution of the quadratic equation
  !!    \f[ a x^2 + b x + c = 0 \f]
  !!    is
  !!    \f[ x_{1/2} = (-b \pm \sqrt(b^2 - 4ac)) / (2a)2 \f]
  !!    and its derivative is
  !!    \f[ dx_{1/2}/dy = (-db/dy \pm (b * db/dy - 2a * dc/dy - 2c * da/dy)
  !!                      / \sqrt(b^2 - 4ac)) / 2a
  !!                      -(-b \pm sqrt(b^2 - 4ac))/(2a**2) * da/dy \f]
  !!
  !!    The sign can be chosen. The default uses the negative sign.
  !
  !> \param[in] "real(r2) :: a, b, c"         The coefficients of the quadratic equation
  !> \param[in] "real(r2) :: da, db, dc"      The derivatives of the coefficients
  !> \param[in] "character, optional :: sgn"  '+' or '-' for sign in solution
  !!     (default: '-')
  !> \retval "real(r2) :: dx"                 The derivative of the solution of
  !!     the quadratic equation
  !
  !> \author Matthias Cuntz
  !> \date Nov 2025
  elemental pure function deriv_solve_quadratic(a, b, c, da, db, dc, sgn)

    use cable_def_types_mod, only: r2
    use mo_utils, only: special_value

    implicit none

    real(r2), intent(in) :: a, b, c, da, db, dc
    character, intent(in), optional :: sgn
    real(r2) :: deriv_solve_quadratic

    real(r2) :: s2, s, p
    character :: isgn

    if (present(sgn)) then
       isgn = sgn
    else
       isgn = '-'
    endif

    if (abs(c) > tiny(1.0_r2)) then
       if (abs(a) > tiny(1.0_r2)) then
          s2 = b**2 - 4.0_r2 * a * c
          if (s2 > tiny(1.0_r2)) then
             s = sqrt(s2)
             p = (2.0_r2 * b * db - 4.0_r2 * c * da - 4.0_r2 * a * dc) &
                  / (2.0_r2 * s)
             if (isgn == '-') then
                deriv_solve_quadratic = (-db - p) / (2.0_r2 * a) &
                     - (-b - s) / (2.0_r2 * a**2) * da
             else
                deriv_solve_quadratic = (-db + p) / (2.0_r2 * a) &
                     - (-b + s) / (2.0_r2 * a**2) * da
             endif
          else
             deriv_solve_quadratic = special_value(c, nan)
          endif
       else
          if (abs(b) > tiny(1.0_r2)) then
             deriv_solve_quadratic = -1.0_r2 / b * dc + c / b**2 * db
          else
             deriv_solve_quadratic = special_value(c, nan)
          endif
       endif
    else
       if (abs(a) > tiny(1.0_r2)) then
          deriv_solve_quadratic = -1.0_r2 / a * db + b / a**2 * da
       else
          deriv_solve_quadratic = 0.0_r2
       endif
    endif

  end function deriv_solve_quadratic

  
  ! Roots of solution of cubic equation: ax**3 + bx**2 + cx + d = 0
  ! p = (3ac - b2) / 3a**2
  ! q = (2b**3 - 9abc + 27a**2d) / 27a**3
  ! x = 2 * sqrt(-p/3) * cos(1/3 acos(3q/2p sqrt(-3/p)) - 2pi/2 * k) - b/3a
  !     with k=0, 1, 2
  ! default: k=1
  elemental pure subroutine solve_cubic_pq(a, b, c, d, p, q)

    use cable_def_types_mod, only: r2
    use mo_utils, only: special_value

    implicit none

    real(r2), intent(in) :: a,b,c,d
    real(r2), intent(out) :: p, q

    if (abs(a) > tiny(1.0_r2)) then
       p = (3.0_r2 * a * c - b**2) / (3.0_r2 * a**2)
       q = (2.0_r2 * b**3 - 9.0_r2 * a * b * c + 27.0_r2 * a**2* d ) &
            / (27.0_r2 * a**3)
    else
       p = special_value(c, nan)
       q = special_value(c, nan)
    endif

  end subroutine solve_cubic_pq

  
  !> \brief Solution of cubic equation
  !
  !> \details
  !!    The solution of the cubic equation
  !!    \f[ a x^3 + b x^2 + c x + d = 0 \f]
  !!    is
  !!    \f[ x = 2 * sqrt(-p/3) * cos(1/3 acos(3q/2p sqrt(-3/p)) - 2\pi/2 * k)
  !!            - b/3a \f]
  !!    with
  !!    \f[ p = (3ac - b2) / 3a^2 \f]
  !!    \f[ q = (2b^3 - 9abc + 27a^2d) / 27a^3 \f]
  !!    and \f[ k=0, 1, 2 \f]
  !!
  !!    The default uses k=1.
  !
  !> \param[in] "real(r2) :: a, b, c, d"  The coefficients of the cubic equation
  !> \param[in] "integer, optional :: k"  k=0, 1, 2 in solution (default: 1)
  !> \retval "real(r2) :: x"              The solution of the cubic equation
  !
  !> \author Matthias Cuntz
  !> \date Nov 2025
  elemental pure function solve_cubic(a, b, c, d, k)

    use cable_def_types_mod, only: r2
    use mo_constants, only: pi => pi_r2
    use mo_utils, only: special_value

    implicit none

    real(r2), intent(in)  :: a, b, c, d
    integer,  intent(in), optional :: k
    real(r2) :: solve_cubic

    real(r2) :: p, q, pq
    real(r2) :: kk

    if (present(k)) then
       kk = real(k, r2)
    else
       kk = 1.0_r2
    endif

    call solve_cubic_pq(a, b, c, d, p, q)

    if (p > tiny(1.0_r2)) then
       ! a must be /= 0 so that p is defined
       pq = 3.0_r2 * q / (2.0_r2 * p) * sqrt(-3.0_r2 / p)
       if (pq < 0.999999999999_r2) then
          solve_cubic = 2.0_r2 * sqrt(-p / 3.0_r2) &
               * cos(acos(pq) / 3.0_r2 - 2.0_r2 * pi * kk /3.0_r2) - b/(3.0_r2 * a)
       else
          solve_cubic = special_value(c, nan)
       endif
    else
       solve_cubic = special_value(c, nan)
    endif

  end function solve_cubic

  
  ! Derivative of roots of solution of cubic equation:
  !    ax**3 + bx**2 + cx + d = 0
  ! default: k=1  
  elemental pure subroutine deriv_solve_cubic_pq(a, b, c, d, da, db, dc, dd, dp, dq)

    use cable_def_types_mod, only: r2
    use mo_utils, only: special_value

    implicit none

    real(r2), intent(in)  :: a, b, c, d, da, db, dc, dd
    real(r2), intent(out) :: dp, dq

    if (abs(a) > tiny(1.0_r2)) then
       dp = (3.0_r2 *da * c + 3.0_r2 * a * dc - 2.0_r2 * b * db)/(3.0_r2 * a**2) &
            - 2.0_r2 * (3.0_r2 * a * c - b**2) / (3.0_r2 * a**3) * da
       dq = (6.0_r2 * b**2 * db - 9.0_r2 * da * b * c - 9.0_r2 * a * db * c &
            - 9.0_r2 * a * b * dc + 54.0_r2 * a * da * d + 27.0_r2 * a**2 * dd) &
            / (27.0_r2 * a**3) &
            - 3.0_r2 * (2.0_r2 * b**3 - 9.0_r2 * a * b * c + 27.0_r2 * a**2 * d) &
            / (27.0_r2 * a**4) * da
    else
       dp = special_value(c, nan)
       dq = special_value(c, nan)
    endif
       
  end subroutine deriv_solve_cubic_pq


  !> \brief Derivative of solution of cubic equation
  !
  !> \details
  !!    The derivative \f[ dx/dy \f] of the solution of the cubic equation
  !!    \f[ a x^3 + b x^2 + c x + d = 0 \f]
  !!
  !!    The default uses k=1 for k in [0, 1, 2].
  !
  !> \param[in] "real(r2) :: a, b, c, d"      The coefficients of the cubic equation
  !> \param[in] "real(r2) :: da, db, dc, dd"  The derivatives of the coefficients
  !> \param[in] "integer, optional :: k"      k=0, 1, 2 in solution (default: 1)
  !> \retval "real(r2) :: dx"                 Derivative of the solution of the
  !!     cubic equation
  !
  !> \author Matthias Cuntz
  !> \date Nov 2025
  elemental pure function deriv_solve_cubic(a, b, c, d, da, db, dc, dd, k)

    use cable_def_types_mod, only: r2
    use mo_constants, only: pi => pi_r2
    use mo_utils, only: special_value

    implicit none

    real(r2), intent(in) :: a, b, c, d, da, db, dc, dd
    integer,  intent(in), optional :: k
    real(r2) :: deriv_solve_cubic

    real(r2) :: p, q, dp, dq
    real(r2) :: pq
    real(r2) :: kk

    if (present(k)) then
       kk = real(k, r2)
    else
       kk = 1.0_r2
    endif

    call solve_cubic_pq(a, b, c, d, p, q)

    if (p > tiny(1.0_r2)) then
       call deriv_solve_cubic_pq(a, b, c, d, da, db, dc, dd, dp, dq)

       pq = 3.0_r2 * q / (2.0_r2 * p) * sqrt(-3.0_r2 / p)
       if (pq < 0.999999999999_r2) then
          deriv_solve_cubic = -1.0_r2 / (3.0_r2 * sqrt(-p / 3.0_r2)) &
               * cos(acos(pq) / 3.0_r2 - 2.0_r2 * pi * kk / 3.0_r2) * dp &
               + 2.0_r2 * sqrt(-p / 3.0_r2) &
               * ( sin(acos(pq) / 3.0_r2 - 2.0_r2 * pi * kk / 3.0_r2) &
               / (3.0_r2 * sqrt(1.0_r2 - pq**2)) &
               * 3.0_r2 / 2.0_r2 &
               * ( dq / p * sqrt(-3.0_r2 / p) &
               - q / p**2 * dp * sqrt(-3.0_r2 / p) &
               + q / p * 3.0_r2 / 2.0_r2 * 1.0_r2 / sqrt(-3.0_r2 / p) &
               * 1.0_r2 / p**2 * dp ) ) &
               - db / (3.0_r2 * a) + b / (3.0_r2 * a**2) * da
       else
          deriv_solve_cubic = special_value(c, nan)
       endif
    else
       deriv_solve_cubic = special_value(c, nan)
    endif

  end function deriv_solve_cubic
  
end module cable_polynomial
