!> \file cable_c13o2_def

!> \brief Constants, types and allocation/zeroing routines for 13CO2 in CABLE

!> \details All type definitions, allocation and output variables for calculating 13C within Cable.

!> \author Matthias Cuntz, Juergen Knauer
!> \date Apr 2019

MODULE cable_c13o2_def

  use cable_def_types_mod, only: dp => r_2

  implicit none

  private

  ! types
  public :: c13o2_flux
  public :: c13o2_pool
  public :: c13o2_luc
  ! routines
  public :: c13o2_alloc_flux
  public :: c13o2_alloc_pools
  public :: c13o2_alloc_luc
  public :: c13o2_zero_flux
  public :: c13o2_update_sum_pools
  public :: c13o2_zero_sum_pools
  public :: c13o2_zero_luc

  ! types
  type c13o2_flux
     integer                         :: ntile
     real(dp), dimension(:), pointer :: ass   ! net assimilation 13CO2 flux [mol/m2s]
  end type c13o2_flux

  type c13o2_pool
     integer                           :: ntile, nplant, nlitter, nsoil, npools
     real(dp), dimension(:,:), pointer :: cplant   ! 13C content in plants: leaves, wood, fine roots
     real(dp), dimension(:,:), pointer :: clitter  ! 13C content in litter: metabolic, fine structural, coarse woody debris
     real(dp), dimension(:,:), pointer :: csoil    ! 13C content in soils: fast, medium, slow
     real(dp), dimension(:),   pointer :: clabile  ! 13C content in excess pool
     real(dp), dimension(:),   pointer :: charvest ! 13C content in agricultural harvest products
  end type c13o2_pool
  
  type c13o2_luc
     integer                           :: nland, nharvest, nclearance, npools
     real(dp), dimension(:,:), pointer :: charvest   ! 13C content in harvest products
     real(dp), dimension(:,:), pointer :: cclearance ! 13C content in clearance products
     real(dp), dimension(:),   pointer :: cagric     ! 13C content in agricultural products
  end type c13o2_luc

  ! variables
  type(c13o2_pool), public :: c13o2pools
  type(c13o2_pool), public :: sum_c13o2pools
  type(c13o2_luc),  public :: c13o2luc

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------

  ! Allocate all 13C fluxes
  subroutine c13o2_alloc_flux(c13o2flux, ntile)

    implicit none
    
    type(c13o2_flux), intent(inout) :: c13o2flux
    integer,          intent(in)    :: ntile

    c13o2flux%ntile = ntile
    allocate(c13o2flux%ass(ntile))
    
  end subroutine c13o2_alloc_flux

  ! ------------------------------------------------------------------

  ! Allocate all 13C Casa pools
  subroutine c13o2_alloc_pools(c13o2pools, ntile)

    use casadimension, only: mplant, mlitter, msoil

    implicit none
    
    type(c13o2_pool), intent(inout) :: c13o2pools
    integer,          intent(in)    :: ntile

    c13o2pools%ntile   = ntile
    c13o2pools%nplant  = mplant
    c13o2pools%nlitter = mlitter
    c13o2pools%nsoil   = msoil
    c13o2pools%npools  = mplant + mlitter + msoil + 1 ! w/o harvest
    allocate(c13o2pools%cplant(ntile,mplant))
    allocate(c13o2pools%clitter(ntile,mlitter))
    allocate(c13o2pools%csoil(ntile,msoil))
    allocate(c13o2pools%clabile(ntile))
    allocate(c13o2pools%charvest(ntile))
    
  end subroutine c13o2_alloc_pools

  ! ------------------------------------------------------------------

  ! Allocate all 13C LUC pools
  subroutine c13o2_alloc_luc(c13o2luc, nland)

    use popluc_module, only: kHarvProd, kClearProd
    
    implicit none
    
    type(c13o2_luc), intent(inout) :: c13o2luc
    integer,         intent(in)    :: nland

    c13o2luc%nland      = nland
    c13o2luc%nharvest   = size(kHarvProd,1)
    c13o2luc%nclearance = size(kClearProd,1)
    c13o2luc%npools     = size(kHarvProd,1) + size(kClearProd,1) + 1
    allocate(c13o2luc%charvest(nland,size(kHarvProd,1)))
    allocate(c13o2luc%cclearance(nland,size(kClearProd,1)))
    allocate(c13o2luc%cagric(nland))

  end subroutine c13o2_alloc_luc

  ! ------------------------------------------------------------------

  subroutine c13o2_zero_flux(c13o2flux)

    use cable_def_types_mod, only: dp => r_2

    implicit none

    type(c13o2_flux), intent(inout) :: c13o2flux

    c13o2flux%ass  = 0._dp

  end subroutine c13o2_zero_flux

  ! ------------------------------------------------------------------

  ! Sum current time step to accumulated output for 13CO2
  ! or divide at end by number of time steps
  subroutine c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, sum_now, average_now, nsteps)

    implicit none
    
    type(c13o2_pool), intent(inout) :: sum_c13o2pools
    type(c13o2_pool), intent(in)    :: c13o2pools
    logical,          intent(in)    :: sum_now, average_now
    integer,          intent(in)    :: nsteps

    real(dp) :: rsteps ! 1/real(nsteps)
    
    if (sum_now) then
       sum_c13o2pools%cplant  =  sum_c13o2pools%cplant   + c13o2pools%cplant
       sum_c13o2pools%clitter  = sum_c13o2pools%clitter  + c13o2pools%clitter
       sum_c13o2pools%csoil    = sum_c13o2pools%csoil    + c13o2pools%csoil
       sum_c13o2pools%clabile  = sum_c13o2pools%clabile  + c13o2pools%clabile
    else if (average_now) then
       rsteps = 1._dp / real(nsteps,dp)
       sum_c13o2pools%cplant   = sum_c13o2pools%cplant  * rsteps
       sum_c13o2pools%clitter  = sum_c13o2pools%clitter * rsteps
       sum_c13o2pools%csoil    = sum_c13o2pools%csoil   * rsteps
       sum_c13o2pools%clabile  = sum_c13o2pools%clabile * rsteps
    endif

  end subroutine c13o2_update_sum_pools

  ! ------------------------------------------------------------------
  
  ! Zero the accumulated output for 13CO2
  subroutine c13o2_zero_sum_pools(sum_c13o2pools)

    implicit none
    
    type(c13o2_pool), intent(inout) :: sum_c13o2pools

    sum_c13o2pools%cplant   = 0._dp
    sum_c13o2pools%clitter  = 0._dp
    sum_c13o2pools%csoil    = 0._dp
    sum_c13o2pools%clabile  = 0._dp

  end subroutine c13o2_zero_sum_pools

  ! ------------------------------------------------------------------

  subroutine c13o2_zero_luc(c13o2luc)

    use cable_def_types_mod, only: dp => r_2

    implicit none

    type(c13o2_luc), intent(inout) :: c13o2luc

    c13o2luc%charvest   = 0._dp
    c13o2luc%cclearance = 0._dp
    c13o2luc%cagric     = 0._dp

  end subroutine c13o2_zero_luc

  ! ------------------------------------------------------------------

END MODULE cable_c13o2_def
