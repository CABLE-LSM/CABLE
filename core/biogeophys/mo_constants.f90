!> \file mo_constants.f90

!> \brief Provides computational, mathematical, physical, and file constants

!> \details Provides computational constants like epsilon, mathematical constants such as Pi,
!> physical constants such as the Stefan-Boltzmann constant, and file units for some standard streams
!> such as standard in.

!> \author Matthias Cuntz
!> \date Nov 2011

module mo_constants

  !  This module contains basic and derived constants
  !
  !  Written  Nov 2011, Matthias Cuntz
  !  Modified Mar 2014, Matthias Cuntz                   - iso_fortran_env
  !           Jun 2018, Matthias Cuntz, Johannes Brenner - constants for standard atmosphere (mtclim)

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2011-2018 Matthias Cuntz - mc (at) macu (dot) de
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

  use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit

  use cable_def_types_mod, only: r1, r2

  implicit none

  ! Computational
  !> epsilon(1.0) in double precision
  real(r2), parameter :: eps_r2 = epsilon(1.0_r2)
  !> epsilon(1.0) in single precision
  real(r1), parameter :: eps_r1 = epsilon(1.0_r1)

  
  ! Mathematical
  !> pi in double precision
  real(r2), parameter :: pi_r2        = 3.141592653589793238462643383279502884197_r2    ! pi
  !> pi in single precision
  real(r1), parameter :: pi_r1        = 3.141592653589793238462643383279502884197_r1
  !> pi/2 in double precision
  real(r2), parameter :: pio2_r2      = 1.57079632679489661923132169163975144209858_r2  ! pi/2
  !> pi/2 in single precision
  real(r1), parameter :: pio2_r1      = 1.57079632679489661923132169163975144209858_r1
  !> 2*pi in double precision
  real(r2), parameter :: twopi_r2     = 6.283185307179586476925286766559005768394_r2    ! 2*pi
  !> 2*pi in single precision
  real(r1), parameter :: twopi_r1     = 6.283185307179586476925286766559005768394_r1
  !> square root of 2 in double precision
  real(r2), parameter :: sqrt2_r2     = 1.41421356237309504880168872420969807856967_r2  ! sqrt(2)
  !> square root of 2 in single precision
  real(r1), parameter :: sqrt2_r1     = 1.41421356237309504880168872420969807856967_r1
  !> 1/3 in double precision
  real(r2), parameter :: onethird_r2  = 0.3333333333333333333333333333333333333_r2      ! 1/3
  !> 1/3 in single precision
  real(r1), parameter :: onethird_r1  = 0.3333333333333333333333333333333333333_r1
  !> 2/3 in double precision
  real(r2), parameter :: twothird_r2  = 1.0_r2 - onethird_r2                            ! 2/3
  !> 2/3 in single precision
  real(r1), parameter :: twothird_r1  = 1.0_r1 - onethird_r1
  !> degree to radian conversion (pi/180) in double precision
  real(r2), parameter :: deg2rad_r2   = pi_r2/180._r2            ! deg2rad
  !> degree to radian conversion (pi/180) in single precision
  real(r1), parameter :: deg2rad_r1   = pi_r1/180._r1
  !> radian to degree conversion (180/pi) in double precision
  real(r2), parameter :: rad2deg_r2   = 180._r2/pi_r2            ! rad2deg
  !> radian to degree conversion (180/pi) in single precision
  real(r1), parameter :: rad2deg_r1   = 180._r1/pi_r1
  !> arc to radian conversion (180/pi) in double precision
  real(r2), parameter :: arc2rad_r2   = pi_r2/12._r2             ! arc2rad
  !> arc to radian conversion (180/pi) in single precision
  real(r1), parameter :: arc2rad_r1   = pi_r1/12._r1
  ! seconds per radian of hour angle in single precision
  real(r2), parameter :: secperrad_r2 = 13750.9871_r2            ! secperrad
  ! seconds per radian of hour angle in single precision
  real(r1), parameter :: secperrad_r1 = 13750.9871_r1
  

  ! Physical
  !> seconds per day [s] in double precision
  real(r2), parameter :: secday_r2      = 86400._r2               ! secday [s]
  !> seconds per day [s] in single precision
  real(r1), parameter :: secday_r1      = 86400._r1
  !> psychrometric constant [kPa K^-1] in double precision
  real(r2), parameter :: psychro_r2     = 0.0646_r2               ! psychrometric constant [kPa c-1]
  !> psychrometric constant [kPa K^-1] in sibgle precision
  real(r1), parameter :: psychro_r1     = 0.0646_r1
  !> gravity accelaration [m^2 s^-1] in double precision
  real(r2), parameter :: gravity_r2     = 9.81_r2                 ! gravity acceleration [m^2/s]
  !> gravity accelaration [m^2 s^-1] in single precision
  real(r1), parameter :: gravity_r1     = 9.81_r1
  !>  solar constant in [J m^-2 s^-1] in double precision
  real(r2), parameter :: solarconst_r2  = 1367._r2                ! solar constant in [W m-2 = kg s-3]
  !>  solar constant in [J m^-2 s^-1] in single precision
  real(r1), parameter :: solarconst_r1  = 1367._r1
  !> specific heat for vaporization of water in [J m-2 mm-1] in double precision
  real(r2), parameter :: specheatet_r2  = 2.45e06_r2              ! specific heat in [W s m-2 mm-1 = kg s-2 mm-1]
  !> specific heat for vaporization of water in [J m-2 mm-1] in single precision
  real(r1), parameter :: specheatet_r1  = 2.45e06_r1
  !> standard temperature [k] in double precision
  real(r2), parameter :: T0_r2          = 273.15_r2               ! celcius <-> kelvin [K]
  !> standard temperature [k] in single precision
  real(r1), parameter :: T0_r1          = 273.15_r1
  !> stefan-boltzmann constant [W m^-2 K^-4] in double precision
  real(r2), parameter :: sigma_r2       = 5.67e-08_r2             ! stefan-boltzmann constant [w/m^2/k^4]
  !> stefan-boltzmann constant [W m^-2 K^-4] in single precision
  real(r1), parameter :: sigma_r1       = 5.67e-08_r1
  ! earth radius [m] in double precision
  real(r2), parameter :: radiusearth_r2 = 6371228._r2
  ! earth radius [m] in single precision
  real(r1), parameter :: radiusearth_r1 = 6371228._r1
  ! radians of earth orbit per julian day
  real(r2), parameter :: rar2erday_r2   = 0.017214_r2
  ! radians of earth orbit per julian day
  real(r1), parameter :: rar2erday_r1   = 0.017214_r1
  ! minimum declination (radians)
  real(r2), parameter :: mindecl_r2     = -0.4092797_r2
  ! minimum declination (radians)
  real(r1), parameter :: mindecl_r1     = -0.4092797_r1
  ! julian day offset of winter solstice
  real(r2), parameter :: daysoff_r2     = 11.25_r2
  ! julian day offset of winter solstice
  real(r1), parameter :: daysoff_r1     = 11.25_r1
  ! molecular weight of air (kg mol-1)
  real(r2), parameter :: Ma_r2          = 28.9644e-3_r2
  ! molecular weight of air (kg mol-1)
  real(r1), parameter :: Ma_r1          = 28.9644e-3_r1
  ! molecular weight of water (kg mol-1)
  real(r2), parameter :: Mw_r2          = 18.0148e-3_r2
  ! molecular weight of water (kg mol-1)
  real(r1), parameter :: Mw_r1          = 18.0148e-3_r1
  ! unitless ratio of molec weights (mw/ma)
  real(r2), parameter :: epsair_r2      = Mw_r2/Ma_r2
  ! unitless ratio of molec weights (mw/ma)
  real(r1), parameter :: epsair_r1      = Mw_r1/Ma_r1
  ! gas law constant (m3 pa mol-1 K-1)
  real(r2), parameter :: R_r2           = 8.3143_r2
  ! gas law constant (m3 pa mol-1 K-1)
  real(r1), parameter :: R_r1           = 8.3143_r1
  ! standard temperature lapse rate (-K m-1)
  real(r2), parameter :: lr_std_r2      = 0.0065_r2
  ! standard temperature lapse rate (-K m-1)
  real(r1), parameter :: lr_std_r1      = 0.0065_r1

  !> Standard atmosphere
  !> standard pressure [Pa] in double precision
  real(r2), parameter :: p0_r2    = 101325._r2      ! standard pressure [Pa]
  !> standard pressure [Pa] in single precision
  real(r1), parameter :: p0_r1    = 101325._r1
  !> standard density  [kg m^-3] in double precision
  real(r2), parameter :: rho0_r2  = 1.225_r2        ! standard air density
  !> standard density  [kg m^-3] in single precision
  real(r1), parameter :: rho0_r1  = 1.225_r1
  !> specific heat capacity of air [J kg^-1 K^-1] in double precision
  real(r2), parameter :: cp0_r2   = 1005.0_r2       ! standard specific heat of air
  !> specific heat capacity of air [J kg^-1 K^-1] in single precision
  real(r1), parameter :: cp0_r1   = 1005.0_r1
  !> standard temp at 0.0 m elevation (K)           ! standard temperature at 0 m asl
  real(r2), parameter :: T_std_r2 = 288.15_r2
  !> standard temp at 0.0 m elevation (K)
  real(r1), parameter :: T_std_r1 = 288.15_r1

  
  ! Numerical Recipes
  !> pi in double precision
  real(r2), parameter :: pi_d    = 3.141592653589793238462643383279502884197_r2      ! pi
  !> pi in single precision
  real(r1), parameter :: pi      = 3.141592653589793238462643383279502884197_r1
  !> pi/2 in double precision
  real(r2), parameter :: pio2_d  = 1.57079632679489661923132169163975144209858_r2    ! pi/2
  !> pi/2 in single precision
  real(r1), parameter :: pio2    = 1.57079632679489661923132169163975144209858_r1
  !> 2*pi in double precision
  real(r2), parameter :: twopi_d = 6.283185307179586476925286766559005768394_r2      ! 2*pi
  !> 2*pi in single precision
  real(r1), parameter :: twopi   = 6.283185307179586476925286766559005768394_r1
  !> square root of 2 in double precision
  real(r2), parameter :: sqrt2_d = 1.41421356237309504880168872420969807856967_r2    ! sqrt(2)
  !> square root of 2 in single precision
  real(r1), parameter :: sqrt2   = 1.41421356237309504880168872420969807856967_r1
  !> Euler's constant in double precision
  real(r2), parameter :: euler_d = 0.5772156649015328606065120900824024310422_r2     ! Euler
  !> Euler's constant in single precision
  real(r1), parameter :: euler   = 0.5772156649015328606065120900824024310422_r1

  ! Standard file units
  !> standard input file unit
  ! integer, parameter :: nin  = 5   ! standard input stream
  integer, parameter :: nin  = input_unit   ! standard input stream
  !> standard output file unit
  ! integer, parameter :: nout = 6   ! standard output stream
  integer, parameter :: nout = output_unit   ! standard output stream
  !> standard error file unit
  ! integer, parameter :: nerr = 0   ! error output stream
  integer, parameter :: nerr = error_unit   ! error output stream
  !> standard file unit for namelist
  integer, parameter :: nnml = 100 ! namelist unit

end module mo_constants
