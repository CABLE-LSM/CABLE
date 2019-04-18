!> \file cable_c13o2_def

!> \brief Constants and types for 13CO2 calculations in CABLE

!> \details All type definitions, allocation and output variables for calculating 13C within Cable.

!> \author Matthias Cuntz, Juergen Knauer
!> \date Apr 2019

MODULE cable_c13o2_def

  use cable_def_types_mod, only: dp => r_2

  implicit none

  private

  public :: c13o2_pool
  ! routines
  public :: alloc_c13o2
  public :: update_sum_c13o2
  public :: zero_sum_c13o2

  ! types
  type c13o2_pool
     integer                           :: nland, nplant, nlitter, nsoil
     real(dp), dimension(:),   pointer :: clabile
     real(dp), dimension(:,:), pointer :: cplant
     real(dp), dimension(:,:), pointer :: clitter
     real(dp), dimension(:,:), pointer :: csoil
  end type c13o2_pool

  ! variables
  type(c13o2_pool), public :: c13o2pools
  type(c13o2_pool), public :: sum_c13o2pools

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------

  ! Allocate all 13CO2 variables
  subroutine alloc_c13o2(c13o2pools,nland)

    use casadimension, only: mplant, mlitter, msoil

    implicit none
    
    type(c13o2_pool), intent(inout) :: c13o2pools
    integer,          intent(in)    :: nland

    c13o2pools%nland   = nland
    c13o2pools%nplant  = mplant
    c13o2pools%nlitter = mlitter
    c13o2pools%nsoil   = msoil
    allocate(c13o2pools%clabile(nland))
    allocate(c13o2pools%cplant(nland,mplant))
    allocate(c13o2pools%clitter(nland,mlitter))
    allocate(c13o2pools%csoil(nland,msoil))
    
  end subroutine alloc_c13o2

  ! ------------------------------------------------------------------

  ! Sum current time step to accumulated output for 13CO2
  ! or divide at end by number of time steps
  subroutine update_sum_c13o2(sum_c13o2pools, c13o2pools, sum_now, average_now, nsteps)

    implicit none
    
    type(c13o2_pool), intent(inout) :: sum_c13o2pools
    type(c13o2_pool), intent(in)    :: c13o2pools
    logical,          intent(in)    :: sum_now, average_now
    integer,          intent(in)    :: nsteps

    real(dp) :: rsteps ! 1/real(nsteps)
    
    if (sum_now) then
       sum_c13o2pools%clabile = sum_c13o2pools%clabile + c13o2pools%clabile
       sum_c13o2pools%cplant  = sum_c13o2pools%cplant  + c13o2pools%cplant
       sum_c13o2pools%clitter = sum_c13o2pools%clitter + c13o2pools%clitter
       sum_c13o2pools%csoil   = sum_c13o2pools%csoil   + c13o2pools%csoil
    else if (average_now) then
       rsteps = 1._dp / real(nsteps,dp)
       sum_c13o2pools%clabile = sum_c13o2pools%clabile * rsteps
       sum_c13o2pools%cplant  = sum_c13o2pools%cplant  * rsteps
       sum_c13o2pools%clitter = sum_c13o2pools%clitter * rsteps
       sum_c13o2pools%csoil   = sum_c13o2pools%csoil   * rsteps
    endif

  end subroutine update_sum_c13o2

  ! ------------------------------------------------------------------
  
  ! Zero the accumulated output for 13CO2
  subroutine zero_sum_c13o2(sum_c13o2pools)

    implicit none
    
    type(c13o2_pool), intent(inout) :: sum_c13o2pools

    sum_c13o2pools%clabile = 0._dp
    sum_c13o2pools%cplant  = 0._dp
    sum_c13o2pools%clitter = 0._dp
    sum_c13o2pools%csoil   = 0._dp

  end subroutine zero_sum_c13o2

END MODULE cable_c13o2_def
