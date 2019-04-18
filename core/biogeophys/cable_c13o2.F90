!> \file cable_c13o2

!> \brief Drivers for 13CO2 calculations in CABLE

!> \details Routines for calculating 13C within Cable.

!> \author Matthias Cuntz, Juergen Knauer
!> \date Apr 2019

MODULE cable_c13o2

  implicit none

  private

  public :: c13o2_init
  public :: c13o2_init_restart
  public :: c13o2_update_pools
  public :: c13o2_save_casapool

  ! ------------------------------------------------------------------

  private :: c13o2_c13o2pools_back
  private :: c13o2_err_handler
  private :: c13o2_fluxmatrix
  private :: c13o2_save_c13o2pools
  private :: c13o2_sinks
  private :: c13o2_sources

contains

  ! ------------------------------------------------------------------

  ! Initialise all 13CO2 pools to 0.
  subroutine c13o2_init(c13o2pools)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool

    implicit none

    type(c13o2_pool), intent(inout) :: c13o2pools

    c13o2pools%clabile = 0._dp
    c13o2pools%cplant  = 0._dp
    c13o2pools%clitter = 0._dp
    c13o2pools%csoil   = 0._dp

  end subroutine c13o2_init

  ! ------------------------------------------------------------------

  ! Read 13CO2 pools from restart file
  subroutine c13o2_init_restart(c13o2_restart_in, c13o2pools)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool
    use netcdf,              only: NF90_OPEN, NF90_NOWRITE, NF90_NOERR, &
         NF90_INQ_VARID, NF90_GET_VAR, NF90_CLOSE

    implicit none

    character(len=*), intent(in)    :: c13o2_restart_in
    type(c13o2_pool), intent(inout) :: c13o2pools

    logical :: iexist
    integer :: status, file_id, var_id
    
    inquire(file=trim(c13o2_restart_in), exist=iexist)
    if (.not. iexist) &
         call c13o2_err_handler('13CO2 restart_in file does not exist: '//trim(c13o2_restart_in))

    write(*,*) 'Read 13CO2 restart_in file: ', trim(c13o2_restart_in)

    status = NF90_OPEN(trim(c13o2_restart_in), NF90_NOWRITE, file_id)
    if (status /= NF90_noerr) call c13o2_err_handler('Error 13CO2 restart_in file: '//trim(c13o2_restart_in))

    !MC13 ToDo: check dimensions
    ! clabile
    status = NF90_INQ_VARID(file_id, 'clabile', var_id)
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Variable clabile not found in restart_in file: '//trim(c13o2_restart_in))
    status = NF90_GET_VAR(file_id, var_id, c13o2pools%clabile, start=(/1/), count=(/c13o2pools%nland/))
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Error reading variable clabile from restart_in file: '//trim(c13o2_restart_in))
    ! cplant
    status = NF90_INQ_VARID(file_id, 'cplant', var_id)
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Variable cplant not found in restart_in file: '//trim(c13o2_restart_in))
    status = NF90_GET_VAR(file_id, var_id, c13o2pools%cplant, &
         start=(/1,1/), count=(/c13o2pools%nland, c13o2pools%nplant/))
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Error reading variable cplant from restart_in file: '//trim(c13o2_restart_in))
    ! clitter
    status = NF90_INQ_VARID(file_id, 'clitter', var_id)
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Variable clitter not found in restart_in file: '//trim(c13o2_restart_in))
    status = NF90_GET_VAR(file_id, var_id, c13o2pools%clitter, &
         start=(/1,1/), count=(/c13o2pools%nland, c13o2pools%nlitter/))
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Error reading variable clitter from restart_in file: '//trim(c13o2_restart_in))
    ! csoil
    status = NF90_INQ_VARID(file_id, 'csoil', var_id)
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Variable csoil not found in restart_in file: '//trim(c13o2_restart_in))
    status = NF90_GET_VAR(file_id, var_id, c13o2pools%csoil, &
         start=(/1,1/), count=(/c13o2pools%nland, c13o2pools%nsoil/))
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Error reading variable csoil from restart_in file: '//trim(c13o2_restart_in))

    status = NF90_CLOSE(file_id)
    if (status /= NF90_noerr) &
         call c13o2_err_handler('Could not close the c13o2 restart_in file: '//trim(c13o2_restart_in))

  end subroutine c13o2_init_restart

  ! ------------------------------------------------------------------

  subroutine c13o2_save_casapool(casapool, casasave)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_pool
    
    implicit none

    type(casa_pool),          intent(in)  :: casapool
    real(dp), dimension(:,:), intent(out) :: casasave

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = size(casapool%cplant,2)
    nlitter = size(casapool%clitter,2)
    nsoil   = size(casapool%csoil,2)

    nstart = 1
    nend   = nplant
    casasave(nstart:nend,:) = casapool%cplant
    nstart = nstart + nplant
    nend   = nend + nlitter
    casasave(nstart:nend,:) = casapool%clitter
    nstart = nstart + nlitter
    nend   = nend + nsoil
    casasave(nstart:nend,:) = casapool%csoil
    nstart = nstart + nsoil
    nend   = nstart
    casasave(nstart,:) = casapool%clabile

  end subroutine c13o2_save_casapool

  ! ------------------------------------------------------------------

  subroutine c13o2_fluxmatrix(casaflux, fluxmatrix)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_flux
    
    implicit none

    type(casa_flux),            intent(in)  :: casaflux
    real(dp), dimension(:,:,:), intent(out) :: fluxmatrix

    integer :: nplant, nlitter, nsoil, nstartr, nendr, nstartc, nendc

    nplant  = size(casaflux%FluxFromPtoL,2)
    nlitter = size(casaflux%FluxFromPtoL,3)
    nsoil   = size(casaflux%FluxFromLtoS,3)

    fluxmatrix = 0._dp
    
    nstartr = 1
    nendr   = nplant
    nstartc = nplant + 1
    nendc   = nplant + nlitter
    fluxmatrix(nstartr:nendr,nstartc:nendc,:) = &
         reshape(casaflux%FluxFromPtoL, shape(fluxmatrix(nstartr:nendr,nstartc:nendc,:)), order=(/3,1,2/))
    nstartr = nstartr + nplant
    nendr   = nendr + nlitter
    nstartc = nstartc + nplant
    nendc   = nendc + nsoil
    fluxmatrix(nstartr:nendr,nstartc:nendc,:) = &
         reshape(casaflux%FluxFromLtoS, shape(fluxmatrix(nstartr:nendr,nstartc:nendc,:)), order=(/3,1,2/))
    nstartr = nstartr + nlitter
    nendr   = nendr + nsoil
    fluxmatrix(nstartr:nendr,nstartc:nendc,:) = &
         reshape(casaflux%FluxFromStoS, shape(fluxmatrix(nstartr:nendr,nstartc:nendc,:)), order=(/3,1,2/))

  end subroutine c13o2_fluxmatrix

  ! ------------------------------------------------------------------

  subroutine c13o2_sources(casaflux, casasources)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_flux
    
    implicit none

    type(casa_flux),          intent(in)  :: casaflux
    real(dp), dimension(:,:), intent(out) :: casasources

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = size(casaflux%fracCalloc,2)
    nlitter = size(casaflux%FluxFromLtoCO2,2)
    nsoil   = size(casaflux%FluxFromStoCO2,2)

    casasources = 0._dp
    casasources(1:nplant,:) = transpose(spread(casaflux%Cnpp, dim=2, ncopies=nplant) * casaflux%fracCalloc)
    nstart = nplant + nlitter + nsoil + 1
    nend   = nstart
    casasources(nstart,:) = casaflux%Cgpp * casaflux%fracClabile

  end subroutine c13o2_sources

  ! ------------------------------------------------------------------

  subroutine c13o2_sinks(casaflux, casasinks)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_flux
    
    implicit none

    type(casa_flux),          intent(in)  :: casaflux
    real(dp), dimension(:,:), intent(out) :: casasinks

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = size(casaflux%FluxFromPtoCO2,2)
    nlitter = size(casaflux%FluxFromLtoCO2,2)
    nsoil   = size(casaflux%FluxFromStoCO2,2)

    nstart = 1
    nend   = nplant
    casasinks(nstart:nend,:) = transpose(casaflux%FluxFromPtoCO2)
    nstart = nstart + nplant
    nend   = nend + nlitter
    casasinks(nstart:nend,:) = transpose(casaflux%FluxFromLtoCO2)
    nstart = nstart + nlitter
    nend   = nend + nsoil
    casasinks(nstart:nend,:) = transpose(casaflux%FluxFromStoCO2)
    nstart = nstart + nsoil
    nend   = nstart
    casasinks(nstart,:) = casaflux%clabloss

  end subroutine c13o2_sinks

  ! ------------------------------------------------------------------

  subroutine c13o2_save_c13o2pools(c13o2pools, c13o2save)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool
    
    implicit none

    type(c13o2_pool),         intent(in)  :: c13o2pools
    real(dp), dimension(:,:), intent(out) :: c13o2save

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil

    nstart = 1
    nend   = nplant
    c13o2save(nstart:nend,:) = transpose(c13o2pools%cplant)
    nstart = nstart + nplant
    nend   = nend + nlitter
    c13o2save(nstart:nend,:) = transpose(c13o2pools%clitter)
    nstart = nstart + nlitter
    nend   = nend + nsoil
    c13o2save(nstart:nend,:) = transpose(c13o2pools%csoil)
    nstart = nstart + nsoil
    nend   = nstart
    c13o2save(nstart,:) = c13o2pools%clabile

  end subroutine c13o2_save_c13o2pools

  ! ------------------------------------------------------------------

  subroutine c13o2_c13o2pools_back(c13o2pools, c13o2save)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool
    
    implicit none

    type(c13o2_pool),       intent(inout) :: c13o2pools
    real(dp), dimension(:,:), intent(in)    :: c13o2save

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil

    nstart = 1
    nend   = nplant
    c13o2pools%cplant = transpose(c13o2save(nstart:nend,:))
    nstart = nstart + nplant
    nend   = nend + nlitter
    c13o2pools%clitter = transpose(c13o2save(nstart:nend,:))
    nstart = nstart + nlitter
    nend   = nend + nsoil
    c13o2pools%csoil = transpose(c13o2save(nstart:nend,:))
    nstart = nstart + nsoil
    nend   = nstart
    c13o2pools%clabile = c13o2save(nstart,:)

  end subroutine c13o2_c13o2pools_back

  ! ------------------------------------------------------------------
  
  subroutine c13o2_update_pools(casasave, casaflux, c13o2pools)

    use cable_def_types_mod,   only: dp => r_2
    use casavariable,          only: casa_flux
    use cable_c13o2_def,       only: c13o2_pool
    use mo_isotope_pool_model, only: isotope_pool_model
    
    implicit none

    real(dp), dimension(:,:), intent(in)    :: casasave
    type(casa_flux),          intent(in)    :: casaflux
    type(c13o2_pool),         intent(inout) :: c13o2pools

    integer :: nland, i
    real(dp), dimension(size(casasave,1),size(casasave,2)) :: c13o2save, casasources, casasinks
    real(dp), dimension(size(casasave,1),size(casasave,1),size(casasave,2)) :: fluxmatrix

    call c13o2_save_c13o2pools(c13o2pools, c13o2save)
    call c13o2_fluxmatrix(casaflux, fluxmatrix)
    call c13o2_sources(casaflux, casasources)
    call c13o2_sinks(casaflux, casasinks)
    call isotope_pool_model(1.0_dp, c13o2save, casasave, fluxmatrix, S=casasources, Si=casasinks)
    call c13o2_c13o2pools_back(c13o2pools, c13o2save)
    
  end subroutine c13o2_update_pools

  ! ------------------------------------------------------------------

  ! Write out error message and stop program
  subroutine c13o2_err_handler(message)
    
    implicit none

    character(len=*), intent(in) :: message

    write(*,*) trim(message)
    stop 9

  end subroutine c13o2_err_handler

     
END MODULE cable_c13o2
