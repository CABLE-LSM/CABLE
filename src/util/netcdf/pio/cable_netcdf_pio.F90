module cable_netcdf_pio_mod
  use cable_netcdf_mod

  use cable_mpi_mod, only: mpi_grp_t
  use cable_abort_module, only: cable_abort

  use pio, only: pio_file_desc_t => file_desc_t
  use pio, only: pio_iosystem_desc_t => iosystem_desc_t
  use pio, only: pio_io_desc_t => io_desc_t
  use pio, only: pio_var_desc_t => var_desc_t
  use pio, only: pio_init
  use pio, only: pio_initdecomp
  use pio, only: pio_createfile
  use pio, only: pio_openfile
  use pio, only: pio_closefile
  use pio, only: pio_syncfile
  use pio, only: pio_def_dim
  use pio, only: pio_def_var
  use pio, only: pio_put_att
  use pio, only: pio_get_att
  use pio, only: pio_put_var
  use pio, only: pio_get_var
  use pio, only: pio_setframe
  use pio, only: pio_write_darray
  use pio, only: pio_read_darray
  use pio, only: pio_strerror
  use pio, only: pio_enddef
  use pio, only: pio_redef
  use pio, only: pio_inq_dimid
  use pio, only: pio_inquire_dimension
  use pio, only: pio_inquire_variable
  use pio, only: pio_inq_varid
  use pio, only: pio_finalize
  use pio, only: PIO_MAX_NAME
  use pio, only: PIO_OFFSET_KIND
  use pio, only: PIO_INT
  use pio, only: PIO_REAL
  use pio, only: PIO_DOUBLE
  use pio, only: PIO_REARR_BOX
  use pio, only: PIO_IOTYPE_NETCDF
  use pio, only: PIO_IOTYPE_NETCDF4C
  use pio, only: PIO_IOTYPE_NETCDF4P
  use pio, only: PIO_CLOBBER
  use pio, only: PIO_NOCLOBBER
  use pio, only: PIO_WRITE
  use pio, only: PIO_NOWRITE
  use pio, only: PIO_UNLIMITED
  use pio, only: PIO_NOERR
  use pio, only: PIO_GLOBAL
  implicit none

  private
  public :: cable_netcdf_pio_io_t

  type, extends(cable_netcdf_decomp_t) :: cable_netcdf_pio_decomp_t
    type(pio_io_desc_t), private :: pio_io_desc
  end type

  type, extends(cable_netcdf_io_t) :: cable_netcdf_pio_io_t
    private
    type(mpi_grp_t) :: mpi_grp
    type(pio_iosystem_desc_t) :: pio_iosystem_desc
  contains
    procedure :: init => cable_netcdf_pio_io_init
    procedure :: finalise => cable_netcdf_pio_io_finalise
    procedure :: create_file => cable_netcdf_pio_io_create_file
    procedure :: open_file => cable_netcdf_pio_io_open_file
    procedure :: create_decomp => cable_netcdf_pio_io_create_decomp
  end type

  interface cable_netcdf_pio_io_t
    procedure cable_netcdf_pio_io_constructor
  end interface

  type, extends(cable_netcdf_file_t) :: cable_netcdf_pio_file_t
    type(pio_file_desc_t), private :: pio_file_desc
  contains
    procedure :: close => cable_netcdf_pio_file_close
    procedure :: end_def => cable_netcdf_pio_file_end_def
    procedure :: redef => cable_netcdf_pio_file_redef
    procedure :: sync => cable_netcdf_pio_file_sync
    procedure :: def_dims => cable_netcdf_pio_file_def_dims
    procedure :: def_var => cable_netcdf_pio_file_def_var
    procedure :: put_att_global_string => cable_netcdf_pio_file_put_att_global_string
    procedure :: put_att_global_int32 => cable_netcdf_pio_file_put_att_global_int32
    procedure :: put_att_global_real32 => cable_netcdf_pio_file_put_att_global_real32
    procedure :: put_att_global_real64 => cable_netcdf_pio_file_put_att_global_real64
    procedure :: put_att_var_string => cable_netcdf_pio_file_put_att_var_string
    procedure :: put_att_var_int32 => cable_netcdf_pio_file_put_att_var_int32
    procedure :: put_att_var_real32 => cable_netcdf_pio_file_put_att_var_real32
    procedure :: put_att_var_real64 => cable_netcdf_pio_file_put_att_var_real64
    procedure :: get_att_global_string => cable_netcdf_pio_file_get_att_global_string
    procedure :: get_att_global_int32 => cable_netcdf_pio_file_get_att_global_int32
    procedure :: get_att_global_real32 => cable_netcdf_pio_file_get_att_global_real32
    procedure :: get_att_global_real64 => cable_netcdf_pio_file_get_att_global_real64
    procedure :: get_att_var_string => cable_netcdf_pio_file_get_att_var_string
    procedure :: get_att_var_int32 => cable_netcdf_pio_file_get_att_var_int32
    procedure :: get_att_var_real32 => cable_netcdf_pio_file_get_att_var_real32
    procedure :: get_att_var_real64 => cable_netcdf_pio_file_get_att_var_real64
    procedure :: inq_dim_len => cable_netcdf_pio_file_inq_dim_len
    procedure :: inq_var_ndims => cable_netcdf_pio_file_inq_var_ndims
    procedure :: put_var_int32_0d => cable_netcdf_pio_file_put_var_int32_0d
    procedure :: put_var_int32_1d => cable_netcdf_pio_file_put_var_int32_1d
    procedure :: put_var_int32_2d => cable_netcdf_pio_file_put_var_int32_2d
    procedure :: put_var_int32_3d => cable_netcdf_pio_file_put_var_int32_3d
    procedure :: put_var_real32_0d => cable_netcdf_pio_file_put_var_real32_0d
    procedure :: put_var_real32_1d => cable_netcdf_pio_file_put_var_real32_1d
    procedure :: put_var_real32_2d => cable_netcdf_pio_file_put_var_real32_2d
    procedure :: put_var_real32_3d => cable_netcdf_pio_file_put_var_real32_3d
    procedure :: put_var_real64_0d => cable_netcdf_pio_file_put_var_real64_0d
    procedure :: put_var_real64_1d => cable_netcdf_pio_file_put_var_real64_1d
    procedure :: put_var_real64_2d => cable_netcdf_pio_file_put_var_real64_2d
    procedure :: put_var_real64_3d => cable_netcdf_pio_file_put_var_real64_3d
    procedure :: write_darray_int32_1d => cable_netcdf_pio_file_write_darray_int32_1d
    procedure :: write_darray_int32_2d => cable_netcdf_pio_file_write_darray_int32_2d
    procedure :: write_darray_int32_3d => cable_netcdf_pio_file_write_darray_int32_3d
    procedure :: write_darray_real32_1d => cable_netcdf_pio_file_write_darray_real32_1d
    procedure :: write_darray_real32_2d => cable_netcdf_pio_file_write_darray_real32_2d
    procedure :: write_darray_real32_3d => cable_netcdf_pio_file_write_darray_real32_3d
    procedure :: write_darray_real64_1d => cable_netcdf_pio_file_write_darray_real64_1d
    procedure :: write_darray_real64_2d => cable_netcdf_pio_file_write_darray_real64_2d
    procedure :: write_darray_real64_3d => cable_netcdf_pio_file_write_darray_real64_3d
    procedure :: get_var_int32_0d => cable_netcdf_pio_file_get_var_int32_0d
    procedure :: get_var_int32_1d => cable_netcdf_pio_file_get_var_int32_1d
    procedure :: get_var_int32_2d => cable_netcdf_pio_file_get_var_int32_2d
    procedure :: get_var_int32_3d => cable_netcdf_pio_file_get_var_int32_3d
    procedure :: get_var_real32_0d => cable_netcdf_pio_file_get_var_real32_0d
    procedure :: get_var_real32_1d => cable_netcdf_pio_file_get_var_real32_1d
    procedure :: get_var_real32_2d => cable_netcdf_pio_file_get_var_real32_2d
    procedure :: get_var_real32_3d => cable_netcdf_pio_file_get_var_real32_3d
    procedure :: get_var_real64_0d => cable_netcdf_pio_file_get_var_real64_0d
    procedure :: get_var_real64_1d => cable_netcdf_pio_file_get_var_real64_1d
    procedure :: get_var_real64_2d => cable_netcdf_pio_file_get_var_real64_2d
    procedure :: get_var_real64_3d => cable_netcdf_pio_file_get_var_real64_3d
    procedure :: read_darray_int32_1d => cable_netcdf_pio_file_read_darray_int32_1d
    procedure :: read_darray_int32_2d => cable_netcdf_pio_file_read_darray_int32_2d
    procedure :: read_darray_int32_3d => cable_netcdf_pio_file_read_darray_int32_3d
    procedure :: read_darray_real32_1d => cable_netcdf_pio_file_read_darray_real32_1d
    procedure :: read_darray_real32_2d => cable_netcdf_pio_file_read_darray_real32_2d
    procedure :: read_darray_real32_3d => cable_netcdf_pio_file_read_darray_real32_3d
    procedure :: read_darray_real64_1d => cable_netcdf_pio_file_read_darray_real64_1d
    procedure :: read_darray_real64_2d => cable_netcdf_pio_file_read_darray_real64_2d
    procedure :: read_darray_real64_3d => cable_netcdf_pio_file_read_darray_real64_3d
  end type

contains

  function type_pio(basetype)
    integer, intent(in) :: basetype
    integer :: type_pio
    select case(basetype)
    case(CABLE_NETCDF_INT)
      type_pio = PIO_INT
    case(CABLE_NETCDF_FLOAT)
      type_pio = PIO_REAL
    case(CABLE_NETCDF_DOUBLE)
      type_pio = PIO_DOUBLE
    case default
      call cable_abort("cable_netcdf_pio_mod: Error: type not supported")
    end select
  end function type_pio

  function iotype_pio(iotype)
    integer, intent(in) :: iotype
    integer :: iotype_pio
    select case(iotype)
    case(CABLE_NETCDF_IOTYPE_CLASSIC)
      iotype_pio = PIO_IOTYPE_NETCDF
    case(CABLE_NETCDF_IOTYPE_NETCDF4C)
      iotype_pio = PIO_IOTYPE_NETCDF4C
    case(CABLE_NETCDF_IOTYPE_NETCDF4P)
      iotype_pio = PIO_IOTYPE_NETCDF4P
    case default
      call cable_abort("cable_netcdf_pio_mod: Error: iotype not supported")
    end select
  end function iotype_pio

  function mode_pio(mode)
    integer, intent(in), optional :: mode
    integer :: mode_pio

    if (.not. present(mode)) then
      mode_pio = PIO_WRITE
      return
    end if

    select case(mode)
    case(CABLE_NETCDF_MODE_CLOBBER)
      mode_pio = PIO_CLOBBER
    case(CABLE_NETCDF_MODE_NOCLOBBER)
      mode_pio = PIO_NOCLOBBER
    case(CABLE_NETCDF_MODE_WRITE)
      mode_pio = PIO_WRITE
    case(CABLE_NETCDF_MODE_NOWRITE)
      mode_pio = PIO_NOWRITE
    case default
      call cable_abort("Error: mode not supported", __FILE__, __LINE__)
    end select

  end function mode_pio

  subroutine check_pio(status)
    integer, intent(in) :: status
    integer :: strerror_status
    character(len=PIO_MAX_NAME) :: err_msg
    if (status /= PIO_NOERR) then
        strerror_status = pio_strerror(status, err_msg)
        call cable_abort(trim(err_msg), file=__FILE__, line=__LINE__)
    end if
  end subroutine check_pio

  subroutine get_start_count_nonoptionals(start_nonopt, count_nonopt, shape, start, count)
    integer, allocatable, intent(out) :: start_nonopt(:), count_nonopt(:)
    integer, intent(in) :: shape(:)
    integer, optional, intent(in) :: start(:), count(:)
    if (present(start)) then
      start_nonopt = start
    else
      allocate(start_nonopt, mold=shape)
      start_nonopt = 1
    end if
    if (present(count)) then
      count_nonopt = count
    else
      allocate(count_nonopt, source=shape)
    end if
  end subroutine

  function cable_netcdf_pio_io_constructor(mpi_grp) result(this)
    type(cable_netcdf_pio_io_t) :: this
    type(mpi_grp_t), intent(in) :: mpi_grp
    this%mpi_grp = mpi_grp
  end function

  subroutine cable_netcdf_pio_io_init(this)
    class(cable_netcdf_pio_io_t), intent(inout) :: this
    ! TODO: get PIO configuration settings
    call pio_init( &
      comp_rank=this%mpi_grp%rank, &
      comp_comm=this%mpi_grp%comm%mpi_val, &
      num_iotasks=1, &
      num_aggregator=0, & ! This argument is obsolete (see https://github.com/NCAR/ParallelIO/issues/1888)
      stride=1, &
      rearr=PIO_REARR_BOX, &
      iosystem=this%pio_iosystem_desc, &
      base=1 &
    )
  end subroutine

  subroutine cable_netcdf_pio_io_finalise(this)
    class(cable_netcdf_pio_io_t), intent(inout) :: this
    integer :: status
    call pio_finalize(this%pio_iosystem_desc, status)
    call check_pio(status)
  end subroutine


  function cable_netcdf_pio_io_create_file(this, path, iotype, mode) result(file)
    class(cable_netcdf_pio_io_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    integer, intent(in) :: iotype
    integer, intent(in), optional :: mode
    class(cable_netcdf_file_t), allocatable :: file
    type(pio_file_desc_t) :: pio_file_desc
    call check_pio(pio_createfile(this%pio_iosystem_desc, pio_file_desc, iotype_pio(iotype), path, mode_pio(mode)))
    file = cable_netcdf_pio_file_t(pio_file_desc)
  end function

  function cable_netcdf_pio_io_open_file(this, path, iotype, mode) result(file)
    class(cable_netcdf_pio_io_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    integer, intent(in) :: iotype
    integer, intent(in), optional :: mode
    class(cable_netcdf_file_t), allocatable :: file
    type(pio_file_desc_t) :: pio_file_desc
    call check_pio(pio_openfile(this%pio_iosystem_desc, pio_file_desc, iotype_pio(iotype), path, mode_pio(mode)))
    file = cable_netcdf_pio_file_t(pio_file_desc)
  end function

  function cable_netcdf_pio_io_create_decomp(this, compmap, dims, type) result(decomp)
    class(cable_netcdf_pio_io_t), intent(inout) :: this
    integer, intent(in) :: compmap(:), dims(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp
    type(pio_io_desc_t) :: pio_io_desc
    call pio_initdecomp( &
      this%pio_iosystem_desc, &
      type_pio(type), &
      dims, &
      compmap, &
      pio_io_desc &
    )
    allocate(decomp, source=cable_netcdf_pio_decomp_t(compmap, dims, type, pio_io_desc=pio_io_desc))
  end function

  subroutine cable_netcdf_pio_file_close(this)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    call pio_closefile(this%pio_file_desc)
  end subroutine

  subroutine cable_netcdf_pio_file_end_def(this)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    call check_pio(pio_enddef(this%pio_file_desc))
  end subroutine

  subroutine cable_netcdf_pio_file_redef(this)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    call check_pio(pio_redef(this%pio_file_desc))
  end subroutine

  subroutine cable_netcdf_pio_file_sync(this)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    call pio_syncfile(this%pio_file_desc)
  end subroutine

  subroutine cable_netcdf_pio_file_def_dims(this, dim_names, dim_lens)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: dim_names(:)
    integer, intent(in) :: dim_lens(:)
    integer :: i, tmp
    do i = 1, size(dim_names)
      if (dim_lens(i) == CABLE_NETCDF_UNLIMITED) then
        call check_pio(pio_def_dim(this%pio_file_desc, dim_names(i), PIO_UNLIMITED, tmp))
      else
        call check_pio(pio_def_dim(this%pio_file_desc, dim_names(i), dim_lens(i), tmp))
      end if
    end do
  end subroutine

  subroutine cable_netcdf_pio_file_def_var(this, var_name, dim_names, type)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: dim_names(:)
    integer, intent(in) :: type
    integer, allocatable :: dimids(:)
    integer :: i
    type(pio_var_desc_t) :: tmp

    if (.not. present(dim_names)) then
      call check_pio(pio_def_var(this%pio_file_desc, var_name, type_pio(type), tmp))
      return
    end if

    allocate(dimids(size(dim_names)))
    do i = 1, size(dimids)
      call check_pio(pio_inq_dimid(this%pio_file_desc, dim_names(i), dimids(i)))
    end do
    call check_pio(pio_def_var(this%pio_file_desc, var_name, type_pio(type), dimids, tmp))

  end subroutine

  subroutine cable_netcdf_pio_file_put_att_global_string(this, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name, att_value
    call check_pio(pio_put_att(this%pio_file_desc, PIO_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_put_att_global_int32(this, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
    call check_pio(pio_put_att(this%pio_file_desc, PIO_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_put_att_global_real32(this, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
    call check_pio(pio_put_att(this%pio_file_desc, PIO_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_put_att_global_real64(this, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
    call check_pio(pio_put_att(this%pio_file_desc, PIO_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_put_att_var_string(this, var_name, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name, att_value
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_att(this%pio_file_desc, var_desc, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_put_att_var_int32(this, var_name, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_att(this%pio_file_desc, var_desc, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_put_att_var_real32(this, var_name, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_att(this%pio_file_desc, var_desc, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_put_att_var_real64(this, var_name, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_att(this%pio_file_desc, var_desc, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_get_att_global_string(this, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    character(len=*), intent(out) :: att_value
    call check_pio(pio_get_att(this%pio_file_desc, PIO_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_get_att_global_int32(this, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
    call check_pio(pio_get_att(this%pio_file_desc, PIO_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_get_att_global_real32(this, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
    call check_pio(pio_get_att(this%pio_file_desc, PIO_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_get_att_global_real64(this, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
    call check_pio(pio_get_att(this%pio_file_desc, PIO_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_get_att_var_string(this, var_name, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    character(len=*), intent(out) :: att_value
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_att(this%pio_file_desc, var_desc, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_get_att_var_int32(this, var_name, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_att(this%pio_file_desc, var_desc, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_get_att_var_real32(this, var_name, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_att(this%pio_file_desc, var_desc, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_get_att_var_real64(this, var_name, att_name, att_value)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_att(this%pio_file_desc, var_desc, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_pio_file_inq_dim_len(this, dim_name, dim_len)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: dim_name
    integer, intent(out) :: dim_len
    integer :: dimid
    call check_pio(pio_inq_dimid(this%pio_file_desc, dim_name, dimid))
    call check_pio(pio_inquire_dimension(this%pio_file_desc, dimid, len=dim_len))
  end subroutine

  subroutine cable_netcdf_pio_file_inq_var_ndims(this, var_name, ndims)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer, intent(out) :: ndims
    integer :: varid
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, varid))
    call check_pio(pio_inquire_variable(this%pio_file_desc, varid, ndims=ndims))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_int32_0d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, [1], start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_int32_1d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_int32_2d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_int32_3d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_real32_0d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, [1], start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_real32_1d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_real32_2d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_real32_3d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_real64_0d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, [1], start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_real64_1d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_real64_2d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_put_var_real64_3d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_put_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_int32(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: status
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    if (present(frame)) then
      call pio_setframe(this%pio_file_desc, var_desc, int(frame, PIO_OFFSET_KIND))
    end if
    select type(decomp)
    type is (cable_netcdf_pio_decomp_t)
      select rank(values)
      rank(1)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      rank(2)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      rank(3)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      end select
      call check_pio(status)
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_pio_decomp_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_int32_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_int32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_int32_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_int32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_int32_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_int32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_real32(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: status
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    if (present(frame)) then
      call pio_setframe(this%pio_file_desc, var_desc, int(frame, PIO_OFFSET_KIND))
    end if
    select type(decomp)
    type is (cable_netcdf_pio_decomp_t)
      select rank(values)
      rank(1)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      rank(2)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      rank(3)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      end select
      call check_pio(status)
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_pio_decomp_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_real32_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_real32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_real32_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_real32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_real32_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_real32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_real64(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: status
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    if (present(frame)) then
      call pio_setframe(this%pio_file_desc, var_desc, int(frame, PIO_OFFSET_KIND))
    end if
    select type(decomp)
    type is (cable_netcdf_pio_decomp_t)
      select rank(values)
      rank(1)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      rank(2)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      rank(3)
        call pio_write_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status, fill_value)
      end select
      call check_pio(status)
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_pio_decomp_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_real64_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_real64(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_real64_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_real64(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_write_darray_real64_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_write_darray_real64(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_int32_0d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, [1], start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_int32_1d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_int32_2d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_int32_3d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_real32_0d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, [1], start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_real32_1d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_real32_2d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_real32_3d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_real64_0d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, [1], start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_real64_1d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_real64_2d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_get_var_real64_3d(this, var_name, values, start, count)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer, allocatable :: start_nonopt(:), count_nonopt(:)
    type(pio_var_desc_t) :: var_desc
    call get_start_count_nonoptionals(start_nonopt, count_nonopt, shape(values), start, count)
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    call check_pio(pio_get_var(this%pio_file_desc, var_desc, start_nonopt, count_nonopt, values))
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_int32(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: status
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    if (present(frame)) then
      call pio_setframe(this%pio_file_desc, var_desc, int(frame, PIO_OFFSET_KIND))
    end if
    select type(decomp)
    type is (cable_netcdf_pio_decomp_t)
      select rank(values)
      rank(1)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      rank(2)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      rank(3)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      end select
      call check_pio(status)
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_pio_decomp_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_int32_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_int32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_int32_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_int32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_int32_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_int32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_real32(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: status
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    if (present(frame)) then
      call pio_setframe(this%pio_file_desc, var_desc, int(frame, PIO_OFFSET_KIND))
    end if
    select type(decomp)
    type is (cable_netcdf_pio_decomp_t)
      select rank(values)
      rank(1)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      rank(2)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      rank(3)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      end select
      call check_pio(status)
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_pio_decomp_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_real32_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_real32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_real32_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_real32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_real32_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_real32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_real64(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: status
    type(pio_var_desc_t) :: var_desc
    call check_pio(pio_inq_varid(this%pio_file_desc, var_name, var_desc))
    if (present(frame)) then
      call pio_setframe(this%pio_file_desc, var_desc, int(frame, PIO_OFFSET_KIND))
    end if
    select type(decomp)
    type is (cable_netcdf_pio_decomp_t)
      select rank(values)
      rank(1)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      rank(2)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      rank(3)
        call pio_read_darray(this%pio_file_desc, var_desc, decomp%pio_io_desc, values, status)
      end select
      call check_pio(status)
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_pio_decomp_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_real64_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_real64(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_real64_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_real64(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_pio_file_read_darray_real64_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_pio_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_pio_file_read_darray_real64(this, var_name, values, decomp, frame)
  end subroutine

end module cable_netcdf_pio_mod
