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

  ! ------------------------------------------------------------------

  private :: c13o2_err_handler

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

  ! Write out error message and stop program
  subroutine c13o2_err_handler(message)
    
    implicit none

    character(len=*), intent(in) :: message

    write(*,*) trim(message)
    stop 9

  end subroutine c13o2_err_handler

     
END MODULE cable_c13o2
