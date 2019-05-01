!> \file cable_c13o2

!> \brief 13CO2 calculations in CABLE

!> \details Routines for calculating, shuffling, reading and writing 13C within Cable.

!> \author Matthias Cuntz, Juergen Knauer
!> \date Apr 2019

MODULE cable_c13o2

  implicit none

  private

  ! ------------------------------------------------------------------
  ! Public
  
  ! Casa
  public :: c13o2_init_pools    ! Initialise all 13C pools to 0
  public :: c13o2_save_casapool ! Save Casa pools before Casa update
  public :: c13o2_update_pools  ! Update 13C Casa pools

  ! LUC
  public :: c13o2_init_luc      ! Initialise all 13C LUC pools to 0
  public :: c13o2_save_luc      ! Save Casa pools, LUC pools, and land-use fractions
  public :: c13o2_update_luc    ! Update 13C Casa and LUC pools

  ! Output
  public :: c13o2_create_output ! Create output netcdf file
  public :: c13o2_write_output  ! Write to output netcdf file
  public :: c13o2_close_output  ! Close output netcdf file

  ! Restart
  public :: c13o2_read_restart_pools  ! Read 13CO2 restart_in_pools file
  public :: c13o2_write_restart_pools ! Write 13CO2 restart_out_pools file
  public :: c13o2_read_restart_luc    ! Read 13CO2 restart_in_luc file
  public :: c13o2_write_restart_luc   ! Write 13CO2 restart_out_luc file

  ! Generell
  public :: c13o2_print_delta_pools ! Print delta values of all Casa pools on screen
  public :: c13o2_print_delta_luc   ! Print delta values of all LUC pools on screen

  ! ------------------------------------------------------------------
  ! Private

  ! Casa
  private :: c13o2_fluxmatrix_pools ! Flux matrix between Casa pools
  private :: c13o2_sinks_pools      ! Sinks of isotope pool model
  private :: c13o2_sources_pools    ! Sources of isotope pool model
  private :: c13o2_rsources_pools   ! Isotope ratios of sources of isotope pool model
  private :: c13o2_save_c13o2pools  ! Save 13C Casa pools before update
  private :: c13o2_c13o2pools_back  ! Write 13C Casa pools to c13o2pools after update

  ! General
  private :: c13o2_err_handler      ! Internal error handler

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------
  
  ! Public Routines
  
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  
  ! Casa
  
  ! ------------------------------------------------------------------

  ! Initialise all 13CO2 pools to 0.
  subroutine c13o2_init_pools(c13o2pools)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool

    implicit none

    type(c13o2_pool), intent(inout) :: c13o2pools

    c13o2pools%cplant   = 0._dp
    c13o2pools%clitter  = 0._dp
    c13o2pools%csoil    = 0._dp
    c13o2pools%clabile  = 0._dp
    c13o2pools%charvest = 0._dp

  end subroutine c13o2_init_pools

  ! ------------------------------------------------------------------

  ! Save old C concentrations of Casa pools
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

    casasave = 0._dp
    
    ! plant
    nstart = 1
    nend   = nplant
    casasave(nstart:nend,:) = casapool%cplant
    ! litter
    nstart = nstart + nplant
    nend   = nend + nlitter
    casasave(nstart:nend,:) = casapool%clitter
    ! soil
    nstart = nstart + nlitter
    nend   = nend + nsoil
    casasave(nstart:nend,:) = casapool%csoil
    ! labile
    nstart = nstart + nsoil
    nend   = nstart
    casasave(nstart,:) = casapool%clabile

  end subroutine c13o2_save_casapool

  ! ------------------------------------------------------------------
  
  ! Calculate the generic isotope pool model, updating 13C concentrations in Casa pools
  subroutine c13o2_update_pools(casasave, casaflux, c13o2pools)

    use cable_def_types_mod,   only: dp => r_2
    use casavariable,          only: casa_flux
    use cable_c13o2_def,       only: c13o2_pool
    use mo_isotope_pool_model, only: isotope_pool_model
    
    implicit none

    real(dp), dimension(:,:), intent(in)    :: casasave
    type(casa_flux),          intent(in)    :: casaflux
    type(c13o2_pool),         intent(inout) :: c13o2pools

    real(dp), dimension(size(casasave,1),size(casasave,2)) :: c13o2save, casasources, casasinks, rsources
    real(dp), dimension(size(casasave,1),size(casasave,1),size(casasave,2)) :: fluxmatrix

    ! update 13C of cumulative harvest
    where (casasave(1,:) > 0._dp) c13o2pools%charvest = c13o2pools%charvest + &
         casaflux%FluxFromPtoHarvest * c13o2pools%cplant(1,:) / casasave(1,:)

    ! save the old 13C concentration in casa pools
    call c13o2_save_c13o2pools(c13o2pools, c13o2save)

    ! Prepare isotope pool model by getting pool exchange fluxes, sources and sinks
    ! exchange fluxes between pools
    call c13o2_fluxmatrix_pools(casaflux, fluxmatrix)
    ! sources such as photosynthesis
    call c13o2_sources_pools(casaflux, casasources)
    ! isotope ratios of sources
    call c13o2_rsources_pools(c13o2pools, casasave, rsources)
    ! sinks such as respiration
    call c13o2_sinks_pools(casaflux, casasinks)

    ! Calc the isotope pool model
    call isotope_pool_model(1.0_dp, c13o2save, casasave, fluxmatrix, S=casasources, Rs=rsources, Si=casasinks)

    ! put new solution into initial c13o2_pool type
    call c13o2_c13o2pools_back(c13o2pools, c13o2save)    
    
  end subroutine c13o2_update_pools

  ! ------------------------------------------------------------------
  
  ! LUC
  
  ! ------------------------------------------------------------------

  ! Initialise all 13CO2 LUC pools to 0.
  subroutine c13o2_init_luc(c13o2luc, c13o2pools, veg, np)

    use cable_def_types_mod, only: veg_parameter_type, dp => r_2, i4 => i_d
    use cable_common_module, only: cable_user
    use casavariable,        only: casa_pool
    use casaparm,            only: leaf, wood, froot
    use cable_c13o2_def,     only: c13o2_luc, c13o2_pool, c13o2_alloc_luc, c13o2_zero_luc
    use mo_isotope,          only: vpdbc13

    implicit none

    type(c13o2_luc),          intent(inout) :: c13o2luc
    type(c13o2_pool),         intent(inout) :: c13o2pools
    type(veg_parameter_type), intent(in)    :: veg
    integer(i4),              intent(in)    :: np

    call c13o2_alloc_luc(c13o2luc, np)
    call c13o2_zero_luc(c13o2luc)
    
    if (cable_user%POPLUC_RunType == 'init') then
       ! zero biomass in secondary forest tiles (both CASA and POP variables)
       where (veg%iLU(:) == 2)
          c13o2pools%cplant(:,leaf)  = 0.01 * vpdbc13
          c13o2pools%cplant(:,wood)  = 0.01 * vpdbc13
          c13o2pools%cplant(:,froot) = 0.01 * vpdbc13
       end where
    else if (cable_user%POPLUC_RunType == 'restart') then
       call c13o2_read_restart_luc(cable_user%c13o2_restart_in_luc, c13o2luc)
    endif

  end subroutine c13o2_init_luc

  ! ------------------------------------------------------------------

  ! Save old C concentrations of Casa pools
  subroutine c13o2_save_luc(casapool, popluc, casasave, lucsave)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_pool
    use popluc_types,        only: popluc_type

    implicit none

    type(casa_pool),          intent(in)  :: casapool
    type(popluc_type),        intent(in)  :: popluc
    real(dp), dimension(:,:), intent(out) :: casasave
    real(dp), dimension(:,:), intent(out) :: lucsave

    integer :: nharvest, nclearance, nstart, nend

    nharvest   = size(popluc%HarvProd,2)
    nclearance = size(popluc%ClearProd,2)

    call c13o2_save_casapool(casapool, casasave)

    lucsave = 0._dp
    
    ! harvest
    nstart = 1
    nend   = nharvest
    lucsave(nstart:nend,:) = popluc%HarvProd
    ! clearance
    nstart = nstart + nharvest
    nend   = nend + nclearance
    lucsave(nstart:nend,:) = popluc%ClearProd
    ! agric
    nstart = nstart + nclearance
    nend   = nstart
    lucsave(nstart,:) = popluc%AgProd

  end subroutine c13o2_save_luc

  ! ------------------------------------------------------------------
  
  ! Calculate the generic isotope pool model, updating 13C concentrations in Casa pools
  subroutine c13o2_update_luc(casasave, lucsave, popluc, prim_only, c13o2pools, c13o2luc)

    use cable_def_types_mod,   only: dp => r_2
    use cable_io_vars_module,  only: landpt, patch
    use casavariable,          only: casa_flux
    use popluc_types,          only: popluc_type
    use popluc_constants,      only: nLU
    use popluc_module,         only: kHarvProd, kClearProd, kAgProd
    use cable_c13o2_def,       only: c13o2_pool, c13o2_luc
    use mo_isotope_pool_model, only: isotope_luc_model, isotope_pool_model

    implicit none

    real(dp), dimension(:,:), intent(in)    :: casasave
    real(dp), dimension(:,:), intent(in)    :: lucsave
    type(popluc_type),        intent(in)    :: popluc
    logical,  dimension(:),   intent(in)    :: prim_only
    type(c13o2_pool),         intent(inout) :: c13o2pools
    type(c13o2_luc),          intent(inout) :: c13o2luc
    
    real(dp), dimension(size(casasave,1),size(casasave,2))   :: c13o2savepools, rsavepools ! 13C conc, isotope ratio
    real(dp), dimension(nLU)                                 :: A, Anew  ! patch fractions
    real(dp), dimension(nLU,nLU)                             :: dA, dAp  ! area change matrix
    real(dp), dimension(nLU)                                 :: c13sluc, tluc ! 13C source, C sink
    real(dp), dimension(c13o2luc%nharvest)                   :: c13sharv ! 13C source
    real(dp), dimension(c13o2luc%nharvest,c13o2luc%nharvest) :: fharv    ! fake harvest flux matrix
    real(dp), dimension(c13o2luc%nclearance)                     :: c13sclear ! 13C source
    real(dp), dimension(c13o2luc%nclearance,c13o2luc%nclearance) :: fclear    ! fake clearance flux matrix
    real(dp), dimension(1)                                       :: c13sag    ! 13C source
    real(dp), dimension(1,1)                                     :: fag       ! fake clearance flux matrix

    integer :: nplant, nlitter, nsoil, nharvest, nclearance
    integer :: g, j, l, c, cs, ce

    nplant     = c13o2pools%nplant
    nlitter    = c13o2pools%nlitter
    nsoil      = c13o2pools%nsoil
    nharvest   = c13o2luc%nharvest
    nclearance = c13o2luc%nclearance

    fharv  = 0._dp
    fclear = 0._dp
    fag    = 0._dp

    ! save old 13C concentrations in Casa and LUC
    call c13o2_save_c13o2pools(c13o2pools, c13o2savepools)
    ! isotope ratio in old pools
    rsavepools = 0._dp
    where(casasave > 0._dp) rsavepools = c13o2savepools / casasave
    
    do g=1, popluc%np ! loop over popluc gridcells == Cable gridcells
       if (.not. prim_only(g)) then ! only then land-use change is possible
          j = landpt(g)%cstart ! start index of Cable tiles for grid cell g
          l = landpt(g)%cend   ! end index of Cable tiles for grid cell g
          if ((l-j+1) /= nLU) then
             write(*,*) 'I do still not understand the tiling: ', g, j, l, nLU
             stop 9
          endif
          if ((popluc%ptos(g) + popluc%ptos(g) + popluc%ptos(g) + popluc%ptos(g)) > 0._dp) then             
             A       = patch(j:l)%frac
             dA      = 0._dp
             dA(1,2) = popluc%ptos(g)
             dA(1,3) = popluc%ptog(g)
             dA(2,3) = popluc%stog(g)
             dA(3,2) = popluc%gtos(g)
             ! plant
             dAp  = 0.
             Anew = A - sum(dA, dim=2) + sum(dA, dim=1)
             do c=1, nplant
                cs = c
                tluc = casasave(j:l,cs) * sum(dA, dim=2)
                call isotope_luc_model(c13o2pools%cplant(j:l,c), A, dAp, C=casasave(j:l,cs), T=tluc, At=Anew)
             end do
             ! litter
             do c=1, nlitter
                cs = nplant + c
                c13sluc = 0._dp
                if (c==2) then
                   c13sluc(:) = sum(dA*spread(c13o2savepools(j:l,1), 2, nLU), dim=1) + & ! spread(array, dim, ncopies)
                        sum(dA*spread(c13o2savepools(j:l,3), 2, nLU), dim=1)
                else if (c==3) then
                   c13sluc(2) = popluc%FluxPHarvResidtoLitter(g) * rsavepools(j,2) + &
                        popluc%FluxSHarvResidtoLitter(g) * rsavepools(j+1,2)
                   c13sluc(3) = popluc%FluxPClearResidtoLitter(g) * rsavepools(j,2) + &
                        popluc%FluxSClearResidtoLitter(g) * rsavepools(j+1,2)
                endif
                call isotope_luc_model(c13o2pools%clitter(j:l,c), A, dA, C=casasave(j:l,cs), S=c13sluc)
             end do
             ! soil
             do c=1, nsoil
                cs = nplant + nlitter + c
                call isotope_luc_model(c13o2pools%csoil(j:l,c), A, dA, C=casasave(j:l,cs))
             end do
             ! labile
             cs = nplant + nlitter + nsoil + 1
             call isotope_luc_model(c13o2pools%clabile(j:l), A, dA, C=casasave(j:l,cs))
             ! harvest
             cs = 1
             ce = cs + nharvest
             c13sharv = popluc%fracHarvProd(g,:) * sum(popluc%FHarvest(g,:) * rsavepools(j:l,2))
             call isotope_pool_model(1.0_dp, c13o2luc%charvest(g,:), lucsave(g,cs:ce), fharv, S=c13sharv, beta=kHarvProd)
             ! clearance
             cs = nharvest + 1
             ce = cs + nclearance
             c13sclear = popluc%fracClearProd(g,:) * sum(popluc%FClearance(g,:) * rsavepools(j:l,2))
             call isotope_pool_model(1.0_dp, c13o2luc%cclearance(g,:), lucsave(g,cs:ce), fclear, S=c13sclear, beta=kClearProd)
             ! agric
             cs = nharvest + nclearance + 1
             ce = cs
             c13sag = sum(c13o2pools%charvest(j:l) * A)
             call isotope_pool_model(1.0_dp, c13o2luc%cagric(g:g), lucsave(g,cs:ce), fag, S=c13sag, beta=(/kAgProd/))
          end if
       end if
    end do
    
  end subroutine c13o2_update_luc

  ! ------------------------------------------------------------------
  
  ! Output
  
  ! ------------------------------------------------------------------

  ! Create 13C Casa output file
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
    dims(1) = 'ntile'
    dims(2) = 'nplant'
    dims(3) = 'nlitter'
    dims(4) = 'nsoil'
    dims(5) = 'time'
    ! dimension lengths
    ldims(1) = c13o2pools%ntile
    ldims(2) = c13o2pools%nplant
    ldims(3) = c13o2pools%nlitter
    ldims(4) = c13o2pools%nsoil
    ldims(5) = nf90_unlimited
    
    ! variable names
    vars(1) = 'time'
    vars(2) = 'latitude'
    vars(3) = 'longitude'
    vars(4) = 'cplant'
    vars(5) = 'clitter'
    vars(6) = 'csoil'
    vars(7) = 'clabile'
    ! vars(8) = 'charvest'
    ! variable long_name
    lvars(1) = 'Time'
    lvars(2) = 'Latitude'
    lvars(3) = 'Longitude'
    lvars(4) = '13C content of plant pools'
    lvars(5) = '13C content of litter pools'
    lvars(6) = '13C content of soil pools'
    lvars(7) = '13C content of excess carbon pool'
    ! lvars(8) = '13C content of accumulated harvest'
    ! variable units
    uvars(1) = trim(timeunits)
    uvars(2) = 'degrees_north'
    uvars(3) = 'degrees_east'
    uvars(4) = 'kg(13C)/m^2'
    uvars(5) = 'kg(13C)/m^2'
    uvars(6) = 'kg(13C)/m^2'
    uvars(7) = 'kg(13C)/m^2'
    ! uvars(8) = 'kg(13C)/m^2'
    ! number of dimensions
    dvars(1) = 1 ! time
    dvars(2) = 1 ! ntile
    dvars(3) = 1 ! ntile
    dvars(4) = 3 ! ntile, nplant, time
    dvars(5) = 3 ! ntile, nlitter, time
    dvars(6) = 3 ! ntile, nsoil, time
    dvars(7) = 2 ! ntile, time
    ! dvars(8) = 2 ! ntile, time
    ! variable type
    tvars(1) = nf90_int
    tvars(2) = nf90_float
    tvars(3) = nf90_float
    tvars(4) = nf90_float
    tvars(5) = nf90_float
    tvars(6) = nf90_float
    tvars(7) = nf90_float
    ! tvars(8) = nf90_float

    ! output file name
    if (len_trim(cable_user%c13o2_outfile) > 0) then
       fname = trim(cable_user%c13o2_restart_out_pools)
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
             fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_'//trim(dum)//'_casa_out_c13o2.nc'
          else
             fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_casa_out_c13o2.nc'
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
       else if (trim(vars(i)) == 'clabile') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(5)
       ! else if (trim(vars(i)) == 'charvest') then
       !    idids(1) = dim_ids(1)
       !    idids(2) = dim_ids(5)
       else
          call c13o2_err_handler('Dimension not known for variable '//trim(vars(i))// &
               ' in c13o2 output file: '//trim(fname))
       endif
       status = nf90_def_var(file_id, trim(vars(i)), tvars(i), idids(1:dvars(i)), var_ids(i))
       if (status /= nf90_noerr) then
          write(*,*) 'Var: ', trim(vars(i)), ', dim_ids: ', idids(1:dvars(i))
          call c13o2_err_handler('Could not define variable in c13o2 output file: '//trim(fname))
       endif
       status = nf90_put_att(file_id, var_ids(i), 'long_name', lvars(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define long_name for '//trim(vars(i))// &
            ' in c13o2 output file: '//trim(fname))
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

  ! Write into 13C Casa output file
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
    nland   = c13o2pools%ntile
    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil
    ! define variables
    do i=1, nvars
       if (trim(vars(i)) == 'time') then
          status = nf90_put_var(file_id, var_ids(i), timestep, &
               start=(/timestep/))
       else if (trim(vars(i)) == 'cplant') then
          status = nf90_put_var(file_id, var_ids(i), c13o2pools%cplant, &
               start=(/1,1,timestep/), count=(/nland,nplant,1/))
       else if (trim(vars(i)) == 'clitter') then
          status = nf90_put_var(file_id, var_ids(i), c13o2pools%clitter, &
               start=(/1,1,timestep/), count=(/nland,nlitter,1/))
       else if (trim(vars(i)) == 'csoil') then
          status = nf90_put_var(file_id, var_ids(i), c13o2pools%csoil, &
               start=(/1,1,timestep/), count=(/nland,nsoil,1/))
       else if (trim(vars(i)) == 'clabile') then
          status = nf90_put_var(file_id, var_ids(i), c13o2pools%clabile, &
               start=(/1,timestep/), count=(/nland,1/))
       ! else if (trim(vars(i)) == 'charvest') then
       !    status = nf90_put_var(file_id, var_ids(i), c13o2pools%charvest, &
       !         start=(/1,timestep/), count=(/nland,1/))
       endif
       if (status /= nf90_noerr) then
          write(*,*) 'Var: ', trim(vars(i)), ', var_id: ', var_ids(i)
          call c13o2_err_handler('Could not put variable in c13o2 output file')
       endif
    end do ! nvars

  end subroutine c13o2_write_output

  ! ------------------------------------------------------------------

  ! Close 13C Casa output file
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
  
  ! Restart

  ! ------------------------------------------------------------------

#ifndef UM_BUILD
  ! Read 13C Casa pools from restart file
  subroutine c13o2_read_restart_pools(c13o2_restart_in_pools, c13o2pools)

    use cable_c13o2_def,     only: c13o2_pool
    use netcdf,              only: nf90_open, nf90_nowrite, nf90_noerr, &
         nf90_inq_varid, nf90_get_var, nf90_close

    implicit none

    character(len=*), intent(in)    :: c13o2_restart_in_pools
    type(c13o2_pool), intent(inout) :: c13o2pools

    logical :: iexist
    integer :: status, file_id, var_id
    
    inquire(file=trim(c13o2_restart_in_pools), exist=iexist)
    if (.not. iexist) &
         call c13o2_err_handler('13CO2 restart_in_pools file does not exist: '//trim(c13o2_restart_in_pools))

    write(*,*) 'Read 13CO2 restart_in_pools file: ', trim(c13o2_restart_in_pools)

    status = nf90_open(trim(c13o2_restart_in_pools), nf90_nowrite, file_id)
    if (status /= nf90_noerr) call c13o2_err_handler('Error 13CO2 restart_in_pools file: '//trim(c13o2_restart_in_pools))

    !MC13 ToDo - dimension check
    ! cplant
    status = nf90_inq_varid(file_id, 'cplant', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable cplant not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%cplant, &
         start=(/1,1/), count=(/c13o2pools%ntile, c13o2pools%nplant/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable cplant from restart_in_pools file: '//trim(c13o2_restart_in_pools))
    ! clitter
    status = nf90_inq_varid(file_id, 'clitter', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable clitter not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%clitter, &
         start=(/1,1/), count=(/c13o2pools%ntile, c13o2pools%nlitter/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable clitter from restart_in_pools file: '//trim(c13o2_restart_in_pools))
    ! csoil
    status = nf90_inq_varid(file_id, 'csoil', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable csoil not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%csoil, &
         start=(/1,1/), count=(/c13o2pools%ntile, c13o2pools%nsoil/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable csoil from restart_in_pools file: '//trim(c13o2_restart_in_pools))
    ! clabile
    status = nf90_inq_varid(file_id, 'clabile', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable clabile not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%clabile, start=(/1/), count=(/c13o2pools%ntile/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable clabile from restart_in_pools file: '//trim(c13o2_restart_in_pools))
    ! charvest
    status = nf90_inq_varid(file_id, 'charvest', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable charvest not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%charvest, start=(/1/), count=(/c13o2pools%ntile/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable charvest from restart_in_pools file: '//trim(c13o2_restart_in_pools))

    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close the c13o2 restart_in_pools file: '//trim(c13o2_restart_in_pools))

  end subroutine c13o2_read_restart_pools
#endif
  
  ! ------------------------------------------------------------------

#ifndef UM_BUILD
  ! Write 13C Casa pools into restart file
  subroutine c13o2_write_restart_pools(c13o2pools)

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
    integer, parameter :: nlandvars = 2
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
    landvars(2)   = 'charvest'
    plantvars(1)  = 'cplant'
    littervars(1) = 'clitter'
    soilvars(1)   = 'csoil'

    ! restart file name
    if (len_trim(cable_user%c13o2_restart_out_pools) > 0) then
       fname = trim(cable_user%c13o2_restart_out_pools)
    else
       fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_c13o2_pools_rst.nc'
    endif

    ! create restart file
    status = nf90_create(fname, nf90_clobber, file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not open c13o2 restart_out_pools file: '//trim(fname))
    write(*,*) 'Writing 13CO2 restart file'
    status = nf90_redef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not redef c13o2 restart_out_pools file: '//trim(fname))
    
    ! global attributes
    write(cyear, fmt='(I4)') CurYear + 1
    status = nf90_put_att(file_id, nf90_global, "Valid restart date", "01.01."//cyear)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global date attribute in c13o2 restart_out_pools file: '//trim(fname))
    
    ! define dimensions
    status = nf90_def_dim(file_id, 'nland', c13o2pools%ntile, land_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define land dimension in c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nplant', c13o2pools%nplant, plant_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define plant dimension in c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nlitter', c13o2pools%nlitter, litter_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define litter dimension in c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nsoil', c13o2pools%nsoil, soil_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define soil dimension in c13o2 restart_out_pools file: '//trim(fname))
    
    ! define variables
    do i=1, nlandvars
       status = nf90_def_var(file_id, trim(landvars(i)), nf90_float, &
            (/land_id/), landvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(landvars(i))// &
            ' in c13o2 restart_out_pools file: '//trim(fname))
    end do
    do i=1, nplantvars
       status = nf90_def_var(file_id, trim(plantvars(i)), nf90_float, &
            (/land_id,plant_id/), plantvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(plantvars(i))// &
            ' in c13o2 restart_out_pools file: '//trim(fname))
    end do
    do i=1, nlittervars
       status = nf90_def_var(file_id, trim(littervars(i)), nf90_float, &
            (/land_id,litter_id/), littervars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(littervars(i))// &
            ' in c13o2 restart_out_pools file: '//trim(fname))
    end do
    do i=1, nsoilvars
       status = nf90_def_var(file_id, trim(soilvars(i)), nf90_float, &
            (/land_id,soil_id/), soilvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(soilvars(i))// &
            ' in c13o2 restart_out_pools file: '//trim(fname))
    end do

    status = nf90_enddef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not end definition phase of c13o2 restart_out_pools file: '//trim(fname))
    
    ! put variables
    status = nf90_put_var(file_id, landvars_id(1), c13o2pools%clabile)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(landvars(1))//' to c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_put_var(file_id, landvars_id(2), c13o2pools%charvest)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(landvars(2))//' to c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_put_var(file_id, plantvars_id(1), c13o2pools%cplant)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(plantvars(1))//' to c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_put_var(file_id, littervars_id(1), c13o2pools%clitter)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(littervars(1))//' to c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_put_var(file_id, soilvars_id(1), c13o2pools%csoil)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(soilvars(1))//' to c13o2 restart_out_pools file: '//trim(fname))

    ! close restart file
    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close c13o2 restart_out_pools file: '//trim(fname))

  end subroutine c13o2_write_restart_pools
#endif

  ! ------------------------------------------------------------------

  ! Read 13C LUC pools from restart file
  subroutine c13o2_read_restart_luc(c13o2_restart_in_luc, c13o2luc)

    use cable_c13o2_def,     only: c13o2_luc
    use netcdf,              only: nf90_open, nf90_nowrite, nf90_noerr, &
         nf90_inq_varid, nf90_get_var, nf90_close

    implicit none

    character(len=*), intent(in)    :: c13o2_restart_in_luc
    type(c13o2_luc),  intent(inout) :: c13o2luc

    logical :: iexist
    integer :: status, file_id, var_id
    
    inquire(file=trim(c13o2_restart_in_luc), exist=iexist)
    if (.not. iexist) &
         call c13o2_err_handler('13CO2 restart_in_luc file does not exist: '//trim(c13o2_restart_in_luc))

    write(*,*) 'Read 13CO2 restart_in_luc file: ', trim(c13o2_restart_in_luc)

    status = nf90_open(trim(c13o2_restart_in_luc), nf90_nowrite, file_id)
    if (status /= nf90_noerr) call c13o2_err_handler('Error 13CO2 restart_in_luc file: '//trim(c13o2_restart_in_luc))

    !MC13 ToDo - dimension check
    ! charvest
    status = nf90_inq_varid(file_id, 'charvest', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable charvest not found in restart_in_luc file: '//trim(c13o2_restart_in_luc))
    status = nf90_get_var(file_id, var_id, c13o2luc%charvest, &
         start=(/1,1/), count=(/c13o2luc%nland, c13o2luc%nharvest/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable charvest from restart_in_luc file: '//trim(c13o2_restart_in_luc))
    ! cclearance
    status = nf90_inq_varid(file_id, 'cclearance', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable cclearance not found in restart_in_luc file: '//trim(c13o2_restart_in_luc))
    status = nf90_get_var(file_id, var_id, c13o2luc%cclearance, &
         start=(/1,1/), count=(/c13o2luc%nland, c13o2luc%nclearance/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable cclearance from restart_in_luc file: '//trim(c13o2_restart_in_luc))
    ! cagric
    status = nf90_inq_varid(file_id, 'cagric', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable cagric not found in restart_in_luc file: '//trim(c13o2_restart_in_luc))
    status = nf90_get_var(file_id, var_id, c13o2luc%cagric, &
         start=(/1/), count=(/c13o2luc%nland/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable cagric from restart_in_luc file: '//trim(c13o2_restart_in_luc))

    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close the c13o2 restart_in_luc file: '//trim(c13o2_restart_in_luc))

  end subroutine c13o2_read_restart_luc
  
  ! ------------------------------------------------------------------

  ! Write 13C Casa pools into restart file
  subroutine c13o2_write_restart_luc(c13o2luc)

    use cable_common_module, only: cable_user, filename, CurYear
    use cable_c13o2_def,     only: c13o2_luc
    use netcdf,              only: nf90_create, nf90_clobber, nf90_noerr, &
         nf90_redef, nf90_put_att, nf90_global, nf90_def_dim, &
         nf90_def_var, nf90_float, nf90_enddef, nf90_put_var, nf90_close

    implicit none

    type(c13o2_luc), intent(in) :: c13o2luc

    ! local variables
    integer :: file_id, land_id, harvest_id, clearance_id
    integer :: i, status
    character(len=200) :: fname, cyear

    ! 1 dim arrays (npt,)
    integer, parameter :: nlandvars = 1
    character(len=20), dimension(nlandvars)      :: landvars
    ! 2 dim arrays (npt,nharvest)
    integer, parameter :: nharvestvars = 1
    character(len=20), dimension(nharvestvars)   :: harvestvars
    ! 2 dim arrays (npt,nclearance)
    integer, parameter :: nclearancevars = 1
    character(len=20), dimension(nclearancevars) :: clearancevars
    ! variable ids
    integer, dimension(nlandvars)      :: landvars_id
    integer, dimension(nharvestvars)   :: harvestvars_id
    integer, dimension(nclearancevars) :: clearancevars_id

    ! Variables
    landvars(1)      = 'cagric'
    harvestvars(1)   = 'charvest'
    clearancevars(1) = 'cclearance'

    ! restart file name
    if (len_trim(cable_user%c13o2_restart_out_luc) > 0) then
       fname = trim(cable_user%c13o2_restart_out_luc)
    else
       fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_c13o2_luc_rst.nc'
    endif

    ! create restart file
    status = nf90_create(fname, nf90_clobber, file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not open c13o2 restart_out_luc file: '//trim(fname))
    write(*,*) 'Writing 13CO2 restart file'
    status = nf90_redef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not redef c13o2 restart_out_luc file: '//trim(fname))
    
    ! global attributes
    write(cyear, fmt='(I4)') CurYear + 1
    status = nf90_put_att(file_id, nf90_global, "Valid restart date", "01.01."//cyear)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global date attribute in c13o2 restart_out_luc file: '//trim(fname))
    
    ! define dimensions
    status = nf90_def_dim(file_id, 'nland', c13o2luc%nland, land_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define land dimension in c13o2 restart_out_luc file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nharvest', c13o2luc%nharvest, harvest_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define harvest dimension in c13o2 restart_out_luc file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nclearance', c13o2luc%nclearance, clearance_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define clearance dimension in c13o2 restart_out_luc file: '//trim(fname))
    
    ! define variables
    do i=1, nlandvars
       status = nf90_def_var(file_id, trim(landvars(i)), nf90_float, &
            (/land_id/), landvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(landvars(i))// &
            ' in c13o2 restart_out_luc file: '//trim(fname))
    end do
    do i=1, nharvestvars
       status = nf90_def_var(file_id, trim(harvestvars(i)), nf90_float, &
            (/land_id,harvest_id/), harvestvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(harvestvars(i))// &
            ' in c13o2 restart_out_luc file: '//trim(fname))
    end do
    do i=1, nclearancevars
       status = nf90_def_var(file_id, trim(clearancevars(i)), nf90_float, &
            (/land_id,clearance_id/), clearancevars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(clearancevars(i))// &
            ' in c13o2 restart_out_luc file: '//trim(fname))
    end do

    status = nf90_enddef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not end definition phase of c13o2 restart_out_luc file: '//trim(fname))
    
    ! put variables
    status = nf90_put_var(file_id, landvars_id(1), c13o2luc%cagric)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(landvars(1))//' to c13o2 restart_out_luc file: '//trim(fname))
    status = nf90_put_var(file_id, harvestvars_id(1), c13o2luc%charvest)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(harvestvars(1))//' to c13o2 restart_out_luc file: '//trim(fname))
    status = nf90_put_var(file_id, clearancevars_id(1), c13o2luc%cclearance)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(clearancevars(1))//' to c13o2 restart_out_luc file: '//trim(fname))

    ! close restart file
    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close c13o2 restart_out_luc file: '//trim(fname))

  end subroutine c13o2_write_restart_luc

  ! ------------------------------------------------------------------
  
  ! Generell
  
  ! ------------------------------------------------------------------

  ! Print 13C delta values of Casa pools on screen
  subroutine c13o2_print_delta_pools(casapool, casaflux, c13o2pools)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_pool
    use casavariable,        only: casa_flux
    use cable_c13o2_def,     only: c13o2_pool
    use mo_isotope,          only: delta1000, vpdbc13

    implicit none

    type(casa_pool),  intent(in) :: casapool
    type(casa_flux),  intent(in) :: casaflux
    type(c13o2_pool), intent(in) :: c13o2pools

    ! integer :: nplant, nlitter, nsoil

    ! nplant  = size(casapool%cplant,2)
    ! nlitter = size(casapool%clitter,2)
    ! nsoil   = size(casapool%csoil,2)

    write(*,*) 'delta-13C of Casa pools'
    write(*,*) '    plant:    ', delta1000(c13o2pools%cplant,   casapool%cplant,   vpdbc13, -1000._dp, tiny(1.0_dp))
    write(*,*) '    litter:   ', delta1000(c13o2pools%clitter,  casapool%clitter,  vpdbc13, -1000._dp, tiny(1.0_dp))
    write(*,*) '    soil:     ', delta1000(c13o2pools%csoil,    casapool%csoil,    vpdbc13, -1000._dp, tiny(1.0_dp))
    write(*,*) '    labile:   ', delta1000(c13o2pools%clabile,  casapool%clabile,  vpdbc13, -1000._dp, tiny(1.0_dp))
    write(*,*) '    Charvest: ', delta1000(c13o2pools%charvest, casaflux%Charvest, vpdbc13, -1000._dp, tiny(1.0_dp))

  end subroutine c13o2_print_delta_pools

  ! ------------------------------------------------------------------

  ! Print 13C delta values of LUC pools on screen
  subroutine c13o2_print_delta_luc(popluc, c13o2luc)

    use cable_def_types_mod, only: dp => r_2
    use popluc_types,        only: popluc_type
    use cable_c13o2_def,     only: c13o2_luc
    use mo_isotope,          only: delta1000, vpdbc13

    implicit none

    type(popluc_type), intent(in) :: popluc
    type(c13o2_luc),   intent(in) :: c13o2luc

    ! integer :: nharvest, nclearance

    ! nharvest   = size(popluc%HarvProd,2)
    ! nclearance = size(popluc%ClearProd,2)

    write(*,*) 'delta-13C of LUC pools'
    write(*,*) '    harvest:   ', delta1000(c13o2luc%charvest,   popluc%HarvProd,  vpdbc13, -1000._dp, tiny(1.0_dp))
    write(*,*) '    clearance: ', delta1000(c13o2luc%cclearance, popluc%ClearProd, vpdbc13, -1000._dp, tiny(1.0_dp))
    write(*,*) '    agric:     ', delta1000(c13o2luc%cagric,     popluc%AgProd,    vpdbc13, -1000._dp, tiny(1.0_dp))

  end subroutine c13o2_print_delta_luc

  ! ------------------------------------------------------------------
  
  ! Private Routines
  
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  
  ! Casa
  
  ! ------------------------------------------------------------------

  ! Exchange matrix of all fluxes between Casa pools
  subroutine c13o2_fluxmatrix_pools(casaflux, fluxmatrix)

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

    ! plant to litter
    nstartr = 1
    nendr   = nplant
    nstartc = nplant + 1
    nendc   = nplant + nlitter
    fluxmatrix(nstartr:nendr,nstartc:nendc,:) = &
         reshape(casaflux%FluxFromPtoL, shape(fluxmatrix(nstartr:nendr,nstartc:nendc,:)), order=(/3,1,2/))
    ! litter to soil
    nstartr = nstartr + nplant
    nendr   = nendr + nlitter
    nstartc = nstartc + nlitter
    nendc   = nendc + nsoil
    fluxmatrix(nstartr:nendr,nstartc:nendc,:) = &
         reshape(casaflux%FluxFromLtoS, shape(fluxmatrix(nstartr:nendr,nstartc:nendc,:)), order=(/3,1,2/))
    ! soil to soil
    nstartr = nstartr + nlitter
    nendr   = nendr + nsoil
    nstartc = nstartc
    nendc   = nendc
    fluxmatrix(nstartr:nendr,nstartc:nendc,:) = &
         reshape(casaflux%FluxFromStoS, shape(fluxmatrix(nstartr:nendr,nstartc:nendc,:)), order=(/3,1,2/))
    ! labile = 0.

  end subroutine c13o2_fluxmatrix_pools

  ! ------------------------------------------------------------------

  ! Sinks of Casa pools others than exchanges to other Casa pools, such as respiration
  subroutine c13o2_sinks_pools(casaflux, casasinks)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_flux
    
    implicit none

    type(casa_flux),          intent(in)  :: casaflux
    real(dp), dimension(:,:), intent(out) :: casasinks

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = size(casaflux%FluxFromPtoCO2,2)
    nlitter = size(casaflux%FluxFromLtoCO2,2)
    nsoil   = size(casaflux%FluxFromStoCO2,2)

    casasinks = 0._dp
    
    ! plant
    nstart = 1
    nend   = nplant
    casasinks(nstart:nend,:) = transpose(casaflux%FluxFromPtoCO2)
    ! litter
    nstart = nstart + nplant
    nend   = nend + nlitter
    casasinks(nstart:nend,:) = transpose(casaflux%FluxFromLtoCO2)
    ! soil
    nstart = nstart + nlitter
    nend   = nend + nsoil
    casasinks(nstart:nend,:) = transpose(casaflux%FluxFromStoCO2)
    ! labile
    nstart = nstart + nsoil
    nend   = nstart
    casasinks(nstart,:) = casaflux%clabloss

  end subroutine c13o2_sinks_pools

  ! ------------------------------------------------------------------

  ! Sources of Casa pools others than exchanges to other Casa pools, such as from photosynthesis
  subroutine c13o2_sources_pools(casaflux, casasources)

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

    ! plant
    casasources(1:nplant,:) = transpose(spread(casaflux%Cnpp, dim=2, ncopies=nplant) * casaflux%fracCalloc)
    ! litter = 0.
    ! soil = 0.
    ! labile
    nstart = nplant + nlitter + nsoil + 1
    nend   = nstart
    casasources(nstart,:) = casaflux%Cgpp * casaflux%fracClabile

  end subroutine c13o2_sources_pools

  ! ------------------------------------------------------------------

  ! Isotope ratios of sources of Casa pools others than exchanges to other Casa pools, such as from photosynthesis
  subroutine c13o2_rsources_pools(c13o2pools, casasave, rsources)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool
    
    implicit none

    type(c13o2_pool),         intent(in)  :: c13o2pools
    real(dp), dimension(:,:), intent(in)  :: casasave
    real(dp), dimension(:,:), intent(out) :: rsources

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil

    rsources = 1._dp

    ! plant
    nstart = 1
    nend   = nplant
    rsources(nstart:nend,:) = 1._dp !MC13 ToDo - discrimination of photosynthesis <- or put into sources directly
    ! litter = 1.
    ! soil = 1.
    ! labile
    nstart = nstart + nplant + nlitter + nsoil
    nend   = nstart
    rsources(nstart,:) = 1._dp      !MC13 ToDo - discrimination of photosynthesis <- or put into sources directly

  end subroutine c13o2_rsources_pools

  ! ------------------------------------------------------------------

  ! Save old 13C concentration of Casa pools
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

    ! plant
    nstart = 1
    nend   = nplant
    c13o2save(nstart:nend,:) = transpose(c13o2pools%cplant)
    ! litter
    nstart = nstart + nplant
    nend   = nend + nlitter
    c13o2save(nstart:nend,:) = transpose(c13o2pools%clitter)
    ! soil
    nstart = nstart + nlitter
    nend   = nend + nsoil
    c13o2save(nstart:nend,:) = transpose(c13o2pools%csoil)
    ! labile
    nstart = nstart + nsoil
    nend   = nstart
    c13o2save(nstart,:) = c13o2pools%clabile

  end subroutine c13o2_save_c13o2pools

  ! ------------------------------------------------------------------

  ! Write back new 13C concentration of Casa pools c13o2_pool type
  subroutine c13o2_c13o2pools_back(c13o2pools, c13o2save)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool
    
    implicit none

    type(c13o2_pool),         intent(inout) :: c13o2pools
    real(dp), dimension(:,:), intent(in)    :: c13o2save

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil

    ! plant
    nstart = 1
    nend   = nplant
    c13o2pools%cplant = transpose(c13o2save(nstart:nend,:))
    ! litter
    nstart = nstart + nplant
    nend   = nend + nlitter
    c13o2pools%clitter = transpose(c13o2save(nstart:nend,:))
    ! soil
    nstart = nstart + nlitter
    nend   = nend + nsoil
    c13o2pools%csoil = transpose(c13o2save(nstart:nend,:))
    ! labile
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
