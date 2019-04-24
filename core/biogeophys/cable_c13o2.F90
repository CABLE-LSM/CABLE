!> \file cable_c13o2

!> \brief Drivers for 13CO2 calculations in CABLE

!> \details Routines for calculating 13C within Cable.

!> \author Matthias Cuntz, Juergen Knauer
!> \date Apr 2019

MODULE cable_c13o2

  implicit none

  private

  ! ------------------------------------------------------------------
  ! Public
  
  ! CASA and LUC
  public :: c13o2_init           ! Initialise all 13C pools to 0
  public :: c13o2_save_casapool  ! Save Casa pools before Casa update
  public :: c13o2_update_pools   ! Update 13C Casa pools

  ! Output
  public :: c13o2_create_output  ! Create output netcdf file
  public :: c13o2_write_output   ! Write to output netcdf file
  public :: c13o2_close_output   ! Close output netcdf file

  ! Restart
  public :: c13o2_read_restart   ! Read 13CO2 restart_in file
  public :: c13o2_write_restart  ! Write 13CO2 restart_out file

  ! ------------------------------------------------------------------
  ! Private

  ! CASA and LUC
  private :: c13o2_fluxmatrix      ! Flux matrix between Casa pools
  private :: c13o2_sinks           ! Sinks of isotope pool model
  private :: c13o2_sources         ! Sources of isotope pool model
  private :: c13o2_save_c13o2pools ! Save 13C Casa pools before update
  private :: c13o2_c13o2pools_back ! Write 13C Casa pools to c13o2pools after update

  ! General
  private :: c13o2_err_handler     ! Internal error handler

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------
  
  ! Public Routines
  
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  
  ! CASA and LUC
  
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
  
  subroutine c13o2_update_pools(casasave, casaflux, c13o2pools)

    use cable_def_types_mod,   only: dp => r_2
    use casavariable,          only: casa_flux
    use cable_c13o2_def,       only: c13o2_pool
    use mo_isotope_pool_model, only: isotope_pool_model
    
    implicit none

    real(dp), dimension(:,:), intent(in)    :: casasave
    type(casa_flux),          intent(in)    :: casaflux
    type(c13o2_pool),         intent(inout) :: c13o2pools

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
  
  ! Output
  
  ! ------------------------------------------------------------------

  subroutine c13o2_create_output(casamet, c13o2pools, file_id, vars, var_ids)

    use cable_common_module,  only: cable_user, filename
    use cable_io_vars_module, only: timeunits, calendar
    use casavariable,         only: casa_met, casafile
    use cable_c13o2_def,      only: c13o2_pool
    use netcdf,               only: nf90_create, nf90_clobber, nf90_noerr, &
         nf90_redef, nf90_put_att, nf90_global, nf90_def_dim, nf90_unlimited, &
         nf90_def_var, nf90_float, nf90_int, nf90_enddef, nf90_put_var

    implicit none

    type(casa_met),                  intent(in)  :: casamet
    type(c13o2_pool),                intent(in)  :: c13o2pools
    integer,                         intent(out) :: file_id
    integer, parameter :: nvars = 7
    character(len=20), dimension(nvars), intent(out) :: vars
    integer,           dimension(nvars), intent(out) :: var_ids

    ! local variables
    integer :: i, status
    character(len=200) :: fname, dum
    integer :: olen

    ! dimensions
    integer, parameter :: ndims = 5
    character(len=20), dimension(ndims) :: dims    ! names
    integer,           dimension(ndims) :: dim_ids ! ids
    integer,           dimension(ndims) :: ldims   ! lengths

    ! variables
    character(len=40), dimension(nvars) :: lvars ! long names
    character(len=20), dimension(nvars) :: uvars ! units
    integer, dimension(nvars) :: dvars           ! number of dimensions
    integer, dimension(nvars) :: tvars           ! type
    integer, dimension(3)     :: idids           ! tmp for dim ids

    ! dimension names
    dims(1) = 'nland'
    dims(2) = 'nplant'
    dims(3) = 'nlitter'
    dims(4) = 'nsoil'
    dims(5) = 'time'
    ! dimension lengths
    ldims(1) = c13o2pools%nland
    ldims(2) = c13o2pools%nplant
    ldims(3) = c13o2pools%nlitter
    ldims(4) = c13o2pools%nsoil
    ldims(5) = nf90_unlimited
    
    ! variable names
    vars(1) = 'time'
    vars(2) = 'latitude'
    vars(3) = 'longitude'
    vars(4) = 'clabile'
    vars(5) = 'cplant'
    vars(6) = 'clitter'
    vars(7) = 'csoil'
    ! variable long_name
    lvars(1) = 'Time'
    lvars(2) = 'Latitude'
    lvars(3) = 'Longitude'
    lvars(4) = '13C content of excess carbon pool'
    lvars(5) = '13C content of plant pools'
    lvars(6) = '13C content of litter pools'
    lvars(7) = '13C content of soil pools'
    ! variable units
    uvars(1) = trim(timeunits)
    uvars(2) = 'degrees_north'
    uvars(3) = 'degrees_east'
    uvars(4) = 'kg(13C)/m^2'
    uvars(5) = 'kg(13C)/m^2'
    uvars(6) = 'kg(13C)/m^2'
    uvars(7) = 'kg(13C)/m^2'
    ! number of dimensions
    dvars(1) = 1 ! time
    dvars(2) = 1 ! nland
    dvars(3) = 1 ! nland
    dvars(4) = 2 ! nland, time
    dvars(5) = 3 ! nland, nplant, time
    dvars(6) = 3 ! nland, nlitter, time
    dvars(7) = 3 ! nland, nsoil, time
    ! variable type
    tvars(1) = nf90_int
    tvars(2) = nf90_float
    tvars(3) = nf90_float
    tvars(4) = nf90_float
    tvars(5) = nf90_float
    tvars(6) = nf90_float
    tvars(7) = nf90_float

    ! output file name
    if (len_trim(cable_user%c13o2_outfile) > 0) then
       fname = trim(cable_user%c13o2_restart_out)
    else
       if (len_trim(casafile%out) > 0) then
          olen = len_trim(casafile%out)
          fname = casafile%out(1:olen-3)//'_c13o2.nc'
       else
          if (len_trim(cable_user%mettype) > 0) then
             write(dum, fmt="(i4,'_',i4)") cable_user%yearstart, cable_user%yearend
             if ((cable_user%yearstart < 1000) .and. (cable_user%yearend < 1000)) then
                write(dum, fmt="(i3,'_',i3)") cable_user%yearstart, cable_user%yearend
             else if (cable_user%yearstart < 1000) then
                write(dum, fmt="(i3,'_',i4)") cable_user%yearstart, cable_user%yearend
             endif
             fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_'//trim(dum)//'_c13o2_casa_out.nc'
          else
             fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_c13o2_casa_out.nc'
          endif
       endif
    endif

    ! create output file
    status = nf90_create(fname, nf90_clobber, file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not open c13o2 output file: '//trim(fname))
    write(*,*) 'Writing 13CO2 output file'
    status = nf90_redef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not redef c13o2 output file: '//trim(fname))
    
    ! global attributes
    status = nf90_put_att(file_id, nf90_global, "StartYear", cable_user%yearstart)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global attribute StartYear in c13o2 output file: '//trim(fname))
    status = nf90_put_att(file_id, nf90_global, "EndYear", cable_user%yearend)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global attribute EndYear in c13o2 output file: '//trim(fname))
    status = nf90_put_att(file_id, nf90_global, "RunIden", cable_user%RunIden)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global attribute RunIden in c13o2 output file: '//trim(fname))
    
    ! define dimensions
    do i=1, ndims
       status = nf90_def_dim(file_id, trim(dims(i)), ldims(i), dim_ids(i))
       if (status /= nf90_noerr) then
          write(*,*) 'Dim: ', trim(dims(i)), ', Size: ', ldims(i)
          call c13o2_err_handler('Could not define dimension in c13o2 output file: '//trim(fname))
       endif
    end do
    
    ! define variables
    do i=1, nvars
       if (trim(vars(i)) == 'time') then
          idids(1) = dim_ids(5)
       else if (trim(vars(i)) == 'latitude') then
          idids(1) = dim_ids(1)
       else if (trim(vars(i)) == 'longitude') then
          idids(1) = dim_ids(1)
       else if (trim(vars(i)) == 'clabile') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(5)
       else if (trim(vars(i)) == 'cplant') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(2)
          idids(3) = dim_ids(5)
       else if (trim(vars(i)) == 'clitter') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(3)
          idids(3) = dim_ids(5)
       else if (trim(vars(i)) == 'csoil') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(4)
          idids(3) = dim_ids(5)
       else
          call c13o2_err_handler('Dimension not known for variable '//trim(vars(i))//' in c13o2 output file: '//trim(fname))
       endif
       status = nf90_def_var(file_id, trim(vars(i)), tvars(i), idids(1:dvars(i)), var_ids(i))
       if (status /= nf90_noerr) then
          write(*,*) 'Var: ', trim(vars(i)), ', dim_ids: ', idids(1:dvars(i))
          call c13o2_err_handler('Could not define variable in c13o2 output file: '//trim(fname))
       endif
       status = nf90_put_att(file_id, var_ids(i), 'long_name', lvars(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define long_name for '//trim(vars(i))//' in c13o2 output file: '//trim(fname))
       status = nf90_put_att(file_id, var_ids(i), 'units', uvars(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define units for '//trim(vars(i))//' in c13o2 output file: '//trim(fname))
       if (trim(vars(i)) == 'time') then
          status = nf90_put_att(file_id, var_ids(i), 'calendar', calendar)
          if (status /= nf90_noerr) &
               call c13o2_err_handler('Could not define calendar for variable time in c13o2 output file: '//trim(fname))
       endif
    end do ! nvars

    ! end definition phase
    status = nf90_enddef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not end definition phase of c13o2 output file: '//trim(fname))

    ! put static variables lat and lon
    do i=1, nvars
       if (trim(vars(i)) == 'latitude') then
          status = nf90_put_var(file_id, var_ids(i), casamet%lat)
          if (status /= nf90_noerr) &
               call c13o2_err_handler('Could not put latitudes to c13o2 output file: '//trim(fname))
       else if (trim(vars(i)) == 'longitude') then
          status = nf90_put_var(file_id, var_ids(i), casamet%lon)
          if (status /= nf90_noerr) &
               call c13o2_err_handler('Could not put longitudes to c13o2 output file: '//trim(fname))
       endif
    end do

  end subroutine c13o2_create_output

  ! ------------------------------------------------------------------

  subroutine c13o2_close_output(file_id)

    use netcdf, only: nf90_close, nf90_noerr
    
    implicit none

    integer, intent(in) :: file_id

    integer :: status

    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close 13CO2 output file.')

  end subroutine c13o2_close_output
  
  ! ------------------------------------------------------------------

  subroutine c13o2_write_output(file_id, vars, var_ids, timestep, c13o2pools)

    use cable_c13o2_def, only: c13o2_pool
    use netcdf,          only: nf90_put_var, nf90_noerr

    implicit none

    integer,                         intent(in) :: file_id
    character(len=20), dimension(:), intent(in) :: vars
    integer,           dimension(:), intent(in) :: var_ids
    integer,                         intent(in) :: timestep
    type(c13o2_pool),                intent(in) :: c13o2pools

    ! local variables
    integer :: i, status
    integer :: nvars, nland, nplant, nlitter, nsoil

    nvars   = size(vars,1)
    nland   = c13o2pools%nland
    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil
    ! define variables
    do i=1, nvars
       if (trim(vars(i)) == 'time') then
          status = nf90_put_var(file_id, var_ids(i), timestep, &
               start=(/timestep/))
       else if (trim(vars(i)) == 'clabile') then
          status = nf90_put_var(file_id, var_ids(i), c13o2pools%clabile, &
               start=(/1,timestep/), count=(/nland,1/))
       else if (trim(vars(i)) == 'cplant') then
          status = nf90_put_var(file_id, var_ids(i), c13o2pools%cplant, &
               start=(/1,1,timestep/), count=(/nland,nplant,1/))
       else if (trim(vars(i)) == 'clitter') then
          status = nf90_put_var(file_id, var_ids(i), c13o2pools%clitter, &
               start=(/1,1,timestep/), count=(/nland,nlitter,1/))
       else if (trim(vars(i)) == 'csoil') then
          status = nf90_put_var(file_id, var_ids(i), c13o2pools%csoil, &
               start=(/1,1,timestep/), count=(/nland,nsoil,1/))
       endif
       if (status /= nf90_noerr) then
          write(*,*) 'Var: ', trim(vars(i)), ', var_id: ', var_ids(i)
          call c13o2_err_handler('Could not put variable in c13o2 output file')
       endif
    end do ! nvars

  end subroutine c13o2_write_output
  
  ! ------------------------------------------------------------------
  
  ! Restart

  ! ------------------------------------------------------------------

  ! Read 13CO2 pools from restart file
#ifndef UM_BUILD
  subroutine c13o2_read_restart(c13o2_restart_in, c13o2pools)

    use cable_c13o2_def,     only: c13o2_pool
    use netcdf,              only: nf90_open, nf90_nowrite, nf90_noerr, &
         nf90_inq_varid, nf90_get_var, nf90_close

    implicit none

    character(len=*), intent(in)    :: c13o2_restart_in
    type(c13o2_pool), intent(inout) :: c13o2pools

    logical :: iexist
    integer :: status, file_id, var_id
    
    inquire(file=trim(c13o2_restart_in), exist=iexist)
    if (.not. iexist) &
         call c13o2_err_handler('13CO2 restart_in file does not exist: '//trim(c13o2_restart_in))

    write(*,*) 'Read 13CO2 restart_in file: ', trim(c13o2_restart_in)

    status = nf90_open(trim(c13o2_restart_in), nf90_nowrite, file_id)
    if (status /= nf90_noerr) call c13o2_err_handler('Error 13CO2 restart_in file: '//trim(c13o2_restart_in))

    !MC13 ToDo - Code dimension check
    ! clabile
    status = nf90_inq_varid(file_id, 'clabile', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable clabile not found in restart_in file: '//trim(c13o2_restart_in))
    status = nf90_get_var(file_id, var_id, c13o2pools%clabile, start=(/1/), count=(/c13o2pools%nland/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable clabile from restart_in file: '//trim(c13o2_restart_in))
    ! cplant
    status = nf90_inq_varid(file_id, 'cplant', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable cplant not found in restart_in file: '//trim(c13o2_restart_in))
    status = nf90_get_var(file_id, var_id, c13o2pools%cplant, &
         start=(/1,1/), count=(/c13o2pools%nland, c13o2pools%nplant/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable cplant from restart_in file: '//trim(c13o2_restart_in))
    ! clitter
    status = nf90_inq_varid(file_id, 'clitter', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable clitter not found in restart_in file: '//trim(c13o2_restart_in))
    status = nf90_get_var(file_id, var_id, c13o2pools%clitter, &
         start=(/1,1/), count=(/c13o2pools%nland, c13o2pools%nlitter/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable clitter from restart_in file: '//trim(c13o2_restart_in))
    ! csoil
    status = nf90_inq_varid(file_id, 'csoil', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable csoil not found in restart_in file: '//trim(c13o2_restart_in))
    status = nf90_get_var(file_id, var_id, c13o2pools%csoil, &
         start=(/1,1/), count=(/c13o2pools%nland, c13o2pools%nsoil/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable csoil from restart_in file: '//trim(c13o2_restart_in))

    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close the c13o2 restart_in file: '//trim(c13o2_restart_in))

  end subroutine c13o2_read_restart
#endif
  
  ! ------------------------------------------------------------------

#ifndef UM_BUILD
  subroutine c13o2_write_restart(c13o2pools)

    use cable_common_module, only: cable_user, filename, CurYear
    use cable_c13o2_def,     only: c13o2_pool
    use netcdf,              only: nf90_create, nf90_clobber, nf90_noerr, &
         nf90_redef, nf90_put_att, nf90_global, nf90_def_dim, &
         nf90_def_var, nf90_float, nf90_enddef, nf90_put_var, nf90_close

    implicit none

    type(c13o2_pool), intent(in) :: c13o2pools

    ! local variables
    integer :: file_id, land_id, plant_id, litter_id, soil_id
    integer :: i, status
    character(len=200) :: fname, cyear

    ! 1 dim arrays (npt,)
    integer, parameter :: nlandvars = 1
    character(len=20), dimension(nlandvars)  :: landvars
    ! 2 dim arrays (npt,nplant)
    integer, parameter :: nplantvars = 1
    character(len=20), dimension(nplantvars)  :: plantvars
    ! 2 dim arrays (npt,nlitter)
    integer, parameter :: nlittervars = 1
    character(len=20), dimension(nlittervars) :: littervars
    ! 2 dim arrays (npt,nsoil)
    integer, parameter :: nsoilvars = 1
    character(len=20), dimension(nsoilvars)   :: soilvars
    ! variable ids
    integer, dimension(nlandvars)   :: landvars_id
    integer, dimension(nplantvars)  :: plantvars_id
    integer, dimension(nlittervars) :: littervars_id
    integer, dimension(nsoilvars)   :: soilvars_id

    ! Variables
    landvars(1)   = 'clabile'
    plantvars(1)  = 'cplant'
    littervars(1) = 'clitter'
    soilvars(1)   = 'csoil'

    ! restart file name
    if (len_trim(cable_user%c13o2_restart_out) > 0) then
       fname = trim(cable_user%c13o2_restart_out)
    else
       fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_c13o2_rst.nc'
    endif

    ! create restart file
    status = nf90_create(fname, nf90_clobber, file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not open c13o2 restart_out file: '//trim(fname))
    write(*,*) 'Writing 13CO2 restart file'
    status = nf90_redef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not redef c13o2 restart_out file: '//trim(fname))
    
    ! global attributes
    write(cyear, fmt='(I4)') CurYear + 1
    status = nf90_put_att(file_id, nf90_global, "Valid restart date", "01.01."//cyear)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global date attribute in c13o2 restart_out file: '//trim(fname))
    
    ! define dimensions
    status = nf90_def_dim(file_id, 'nland', c13o2pools%nland, land_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define land dimension in c13o2 restart_out file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nplant', c13o2pools%nplant, plant_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define plant dimension in c13o2 restart_out file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nlitter', c13o2pools%nlitter, litter_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define litter dimension in c13o2 restart_out file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nsoil', c13o2pools%nsoil, soil_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define soil dimension in c13o2 restart_out file: '//trim(fname))
    
    ! define variables
    do i=1, nlandvars
       status = nf90_def_var(file_id, trim(landvars(i)), nf90_float, &
            (/land_id/), landvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(landvars(i))//' in c13o2 restart_out file: '//trim(fname))
    end do
    do i=1, nplantvars
       status = nf90_def_var(file_id, trim(plantvars(i)), nf90_float, &
            (/land_id,plant_id/), plantvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(plantvars(i))//' in c13o2 restart_out file: '//trim(fname))
    end do
    do i=1, nlittervars
       status = nf90_def_var(file_id, trim(littervars(i)), nf90_float, &
            (/land_id,litter_id/), littervars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(littervars(i))//' in c13o2 restart_out file: '//trim(fname))
    end do
    do i=1, nsoilvars
       status = nf90_def_var(file_id, trim(soilvars(i)), nf90_float, &
            (/land_id,soil_id/), soilvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(soilvars(i))//' in c13o2 restart_out file: '//trim(fname))
    end do

    status = nf90_enddef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not end definition phase of c13o2 restart_out file: '//trim(fname))
    
    ! put variables
    status = nf90_put_var(file_id, landvars_id(1), c13o2pools%clabile)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(landvars(1))//' to c13o2 restart_out file: '//trim(fname))
    status = nf90_put_var(file_id, plantvars_id(1), c13o2pools%cplant)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(plantvars(1))//' to c13o2 restart_out file: '//trim(fname))
    status = nf90_put_var(file_id, littervars_id(1), c13o2pools%clitter)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(littervars(1))//' to c13o2 restart_out file: '//trim(fname))
    status = nf90_put_var(file_id, soilvars_id(1), c13o2pools%csoil)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(soilvars(1))//' to c13o2 restart_out file: '//trim(fname))

    ! close restart file
    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close c13o2 restart_out file: '//trim(fname))

  end subroutine c13o2_write_restart
#endif

  ! ------------------------------------------------------------------
  
  ! Private Routines
  
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  
  ! CASA and LUC
  
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
  
  ! General
  
  ! ------------------------------------------------------------------

  ! Write out error message and stop program
  subroutine c13o2_err_handler(message)
    
    implicit none

    character(len=*), intent(in) :: message

    write(*,*) trim(message)
    stop 9

  end subroutine c13o2_err_handler
     
END MODULE cable_c13o2
