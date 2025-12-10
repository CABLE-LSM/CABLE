module cable_restart_mod
  use iso_fortran_env, only: int32, real32, real64

  use cable_abort_module, only: cable_abort

  use cable_def_types_mod, only: mp, mp_global
  use cable_def_types_mod, only: mland_global
  use cable_def_types_mod, only: ms, msn, nrb, ncp, ncs

  use cable_io_vars_module, only: patch_decomp_start
  use cable_io_vars_module, only: patch
  use cable_io_vars_module, only: timeunits, calendar, time_coord

  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: cable_netcdf_file_t
  use cable_netcdf_mod, only: cable_netcdf_create_file
  use cable_netcdf_mod, only: CABLE_NETCDF_IOTYPE_CLASSIC
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE

  use cable_netcdf_decomp_util_mod, only: dim_spec_t
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_patch
  implicit none
  private

  public :: cable_restart_mod_init
  public :: cable_restart_mod_end
  public :: cable_restart_variable_write
  public :: cable_restart_variable_write_darray
  public :: cable_restart_write_time

  ! TODO(Sean): is an interface overkill here? It does make things more intuitive for distributed I/O

  interface cable_restart_variable_write_darray
    module procedure cable_restart_variable_write_darray_int32_1d
    module procedure cable_restart_variable_write_darray_int32_2d
    module procedure cable_restart_variable_write_darray_int32_3d
    module procedure cable_restart_variable_write_darray_real32_1d
    module procedure cable_restart_variable_write_darray_real32_2d
    module procedure cable_restart_variable_write_darray_real32_3d
    module procedure cable_restart_variable_write_darray_real64_1d
    module procedure cable_restart_variable_write_darray_real64_2d
    module procedure cable_restart_variable_write_darray_real64_3d
  end interface cable_restart_variable_write_darray

  interface cable_restart_variable_write
    module procedure cable_restart_variable_write_int32_1d
    module procedure cable_restart_variable_write_int32_2d
    module procedure cable_restart_variable_write_int32_3d
    module procedure cable_restart_variable_write_real32_1d
    module procedure cable_restart_variable_write_real32_2d
    module procedure cable_restart_variable_write_real32_3d
    module procedure cable_restart_variable_write_real64_1d
    module procedure cable_restart_variable_write_real64_2d
    module procedure cable_restart_variable_write_real64_3d
  end interface cable_restart_variable_write

  class(cable_netcdf_file_t), allocatable :: restart_output_file

  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_int32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_real32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_real64
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_soil_int32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_soil_real32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_soil_real64
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_snow_int32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_snow_real32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_snow_real64
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_rad_int32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_rad_real32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_rad_real64
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_plantcarbon_int32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_plantcarbon_real32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_plantcarbon_real64
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_soilcarbon_int32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_soilcarbon_real32
  class(cable_netcdf_decomp_t), allocatable, target :: decomp_patch_soilcarbon_real64

contains

  subroutine cable_restart_mod_init()

    type(dim_spec_t), allocatable :: mem_shape_patch(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_soil(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_snow(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_rad(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_plantcarbon(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_soilcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_patch(:)
    type(dim_spec_t), allocatable :: var_shape_patch_soil(:)
    type(dim_spec_t), allocatable :: var_shape_patch_snow(:)
    type(dim_spec_t), allocatable :: var_shape_patch_rad(:)
    type(dim_spec_t), allocatable :: var_shape_patch_plantcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_patch_soilcarbon(:)

    mem_shape_patch             = [dim_spec_t('patch', mp)]
    mem_shape_patch_soil        = [dim_spec_t('patch', mp), dim_spec_t('soil', ms)]
    mem_shape_patch_snow        = [dim_spec_t('patch', mp), dim_spec_t('snow', msn)]
    mem_shape_patch_rad         = [dim_spec_t('patch', mp), dim_spec_t('rad', nrb)]
    mem_shape_patch_plantcarbon = [dim_spec_t('patch', mp), dim_spec_t('plantcarbon', ncp)]
    mem_shape_patch_soilcarbon  = [dim_spec_t('patch', mp), dim_spec_t('soilcarbon', ncs)]

    var_shape_patch             = [dim_spec_t('patch', mp_global)]
    var_shape_patch_soil        = [dim_spec_t('patch', mp_global), dim_spec_t('soil', ms)]
    var_shape_patch_snow        = [dim_spec_t('patch', mp_global), dim_spec_t('snow', msn)]
    var_shape_patch_rad         = [dim_spec_t('patch', mp_global), dim_spec_t('rad', nrb)]
    var_shape_patch_plantcarbon = [dim_spec_t('patch', mp_global), dim_spec_t('plantcarbon', ncp)]
    var_shape_patch_soilcarbon  = [dim_spec_t('patch', mp_global), dim_spec_t('soilcarbon', ncs)]

    decomp_patch_int32              = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch, var_shape_patch, CABLE_NETCDF_INT)
    decomp_patch_real32             = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch, var_shape_patch, CABLE_NETCDF_FLOAT)
    decomp_patch_real64             = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch, var_shape_patch, CABLE_NETCDF_DOUBLE)
    decomp_patch_soil_int32         = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soil, var_shape_patch_soil, CABLE_NETCDF_INT)
    decomp_patch_soil_real32        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soil, var_shape_patch_soil, CABLE_NETCDF_FLOAT)
    decomp_patch_soil_real64        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soil, var_shape_patch_soil, CABLE_NETCDF_DOUBLE)
    decomp_patch_snow_int32         = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_snow, var_shape_patch_snow, CABLE_NETCDF_INT)
    decomp_patch_snow_real32        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_snow, var_shape_patch_snow, CABLE_NETCDF_FLOAT)
    decomp_patch_snow_real64        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_snow, var_shape_patch_snow, CABLE_NETCDF_DOUBLE)
    decomp_patch_rad_int32          = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_rad, var_shape_patch_rad, CABLE_NETCDF_INT)
    decomp_patch_rad_real32         = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_rad, var_shape_patch_rad, CABLE_NETCDF_FLOAT)
    decomp_patch_rad_real64         = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_rad, var_shape_patch_rad, CABLE_NETCDF_DOUBLE)
    decomp_patch_plantcarbon_int32  = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_plantcarbon, var_shape_patch_plantcarbon, CABLE_NETCDF_INT)
    decomp_patch_plantcarbon_real32 = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_plantcarbon, var_shape_patch_plantcarbon, CABLE_NETCDF_FLOAT)
    decomp_patch_plantcarbon_real64 = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_plantcarbon, var_shape_patch_plantcarbon, CABLE_NETCDF_DOUBLE)
    decomp_patch_soilcarbon_int32   = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soilcarbon, var_shape_patch_soilcarbon, CABLE_NETCDF_INT)
    decomp_patch_soilcarbon_real32  = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soilcarbon, var_shape_patch_soilcarbon, CABLE_NETCDF_FLOAT)
    decomp_patch_soilcarbon_real64  = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soilcarbon, var_shape_patch_soilcarbon, CABLE_NETCDF_DOUBLE)

    restart_output_file = cable_netcdf_create_file("test_restart.nc", iotype=CABLE_NETCDF_IOTYPE_CLASSIC) ! TODO(Sean): use filename from namelist

    call restart_output_file%def_dims(["mland"], [mland_global])
    call restart_output_file%def_dims(["mp"], [mp_global])
    call restart_output_file%def_dims(["soil"], [ms])
    call restart_output_file%def_dims(["snow"], [msn])
    call restart_output_file%def_dims(["rad"], [nrb])
    call restart_output_file%def_dims(["soil_carbon_pools"], [ncs])
    call restart_output_file%def_dims(["plant_carbon_pools"], [ncp])
    call restart_output_file%def_dims(["time"], [1])

    call restart_output_file%end_def()

  end subroutine cable_restart_mod_init

  subroutine cable_restart_mod_end()

    if (allocated(restart_output_file)) call restart_output_file%close()

  end subroutine cable_restart_mod_end

  subroutine define_variable(output_file, var_name, var_dims, var_type, long_name, units)
    class(cable_netcdf_file_t), intent(inout) :: output_file
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call output_file%redef()
    call output_file%def_var(var_name, var_dims, var_type)
    call output_file%put_att(var_name, "long_name", long_name)
    call output_file%put_att(var_name, "units", units)
    call output_file%end_def()

  end subroutine define_variable

  subroutine associate_decomp_int32(var_name, decomp, data_shape)
    character(len=*), intent(in) :: var_name
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    integer, dimension(:), intent(in) :: data_shape

    if (all(data_shape == [mp])) then
      decomp => decomp_patch_int32
    else if (all(data_shape == [mp, ms])) then
      decomp => decomp_patch_soil_int32
    else if (all(data_shape == [mp, msn])) then
      decomp => decomp_patch_snow_int32
    else if (all(data_shape == [mp, nrb])) then
      decomp => decomp_patch_rad_int32
    else if (all(data_shape == [mp, ncp])) then
      decomp => decomp_patch_plantcarbon_int32
    else if (all(data_shape == [mp, ncs])) then
      decomp => decomp_patch_soilcarbon_int32
    else
      call cable_abort("Unexpected data shape for variable " // var_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_int32

  subroutine associate_decomp_real32(var_name, decomp, data_shape)
    character(len=*), intent(in) :: var_name
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    integer, dimension(:), intent(in) :: data_shape

    if (all(data_shape == [mp])) then
      decomp => decomp_patch_real32
    else if (all(data_shape == [mp, ms])) then
      decomp => decomp_patch_soil_real32
    else if (all(data_shape == [mp, msn])) then
      decomp => decomp_patch_snow_real32
    else if (all(data_shape == [mp, nrb])) then
      decomp => decomp_patch_rad_real32
    else if (all(data_shape == [mp, ncp])) then
      decomp => decomp_patch_plantcarbon_real32
    else if (all(data_shape == [mp, ncs])) then
      decomp => decomp_patch_soilcarbon_real32
    else
      call cable_abort("Unexpected data shape for variable " // var_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_real32

  subroutine associate_decomp_real64(var_name, decomp, data_shape)
    character(len=*), intent(in) :: var_name
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    integer, dimension(:), intent(in) :: data_shape

    if (all(data_shape == [mp])) then
      decomp => decomp_patch_real64
    else if (all(data_shape == [mp, ms])) then
      decomp => decomp_patch_soil_real64
    else if (all(data_shape == [mp, msn])) then
      decomp => decomp_patch_snow_real64
    else if (all(data_shape == [mp, nrb])) then
      decomp => decomp_patch_rad_real64
    else if (all(data_shape == [mp, ncp])) then
      decomp => decomp_patch_plantcarbon_real64
    else if (all(data_shape == [mp, ncs])) then
      decomp => decomp_patch_soilcarbon_real64
    else
      call cable_abort("Unexpected data shape for variable " // var_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_real64

  subroutine cable_restart_variable_write_darray_int32_1d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    integer(kind=int32), intent(in) :: data(:)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_int32(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_int32_1d

  subroutine cable_restart_variable_write_darray_int32_2d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    integer(kind=int32), intent(in) :: data(:, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_int32(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_int32_2d

  subroutine cable_restart_variable_write_darray_int32_3d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    integer(kind=int32), intent(in) :: data(:, :, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_int32(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_int32_3d

  subroutine cable_restart_variable_write_darray_real32_1d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real32), intent(in) :: data(:)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_real32(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_real32_1d

  subroutine cable_restart_variable_write_darray_real32_2d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real32), intent(in) :: data(:, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_real32(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_real32_2d

  subroutine cable_restart_variable_write_darray_real32_3d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real32), intent(in) :: data(:, :, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_real32(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_real32_3d

  subroutine cable_restart_variable_write_darray_real64_1d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real64), intent(in) :: data(:)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_real32(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_real64_1d

  subroutine cable_restart_variable_write_darray_real64_2d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real64), intent(in) :: data(:, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_real64(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_real64_2d

  subroutine cable_restart_variable_write_darray_real64_3d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real64), intent(in) :: data(:, :, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units
    class(cable_netcdf_decomp_t), pointer :: decomp

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call associate_decomp_real64(var_name, decomp, shape(data))
    call restart_output_file%write_darray(var_name, data, decomp)

  end subroutine cable_restart_variable_write_darray_real64_3d

  subroutine cable_restart_variable_write_int32_1d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    integer(kind=int32), intent(in) :: data(:)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_int32_1d

  subroutine cable_restart_variable_write_int32_2d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    integer(kind=int32), intent(in) :: data(:, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_int32_2d

  subroutine cable_restart_variable_write_int32_3d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    integer(kind=int32), intent(in) :: data(:, :, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_int32_3d

  subroutine cable_restart_variable_write_real32_1d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real32), intent(in) :: data(:)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_real32_1d

  subroutine cable_restart_variable_write_real32_2d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real32), intent(in) :: data(:, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_real32_2d

  subroutine cable_restart_variable_write_real32_3d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real32), intent(in) :: data(:, :, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_real32_3d

  subroutine cable_restart_variable_write_real64_1d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real64), intent(in) :: data(:)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_real64_1d

  subroutine cable_restart_variable_write_real64_2d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real64), intent(in) :: data(:, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_real64_2d

  subroutine cable_restart_variable_write_real64_3d(var_name, var_dims, data, var_type, long_name, units)
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: var_dims(:)
    real(kind=real64), intent(in) :: data(:, :, :)
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: units

    call define_variable(restart_output_file, var_name, var_dims, var_type, long_name, units)
    call restart_output_file%put_var(var_name, data)

  end subroutine cable_restart_variable_write_real64_3d

  subroutine cable_restart_write_time(time_value)
    real, intent(in) :: time_value

    call restart_output_file%redef()
    call restart_output_file%def_var("time", ["time"], CABLE_NETCDF_DOUBLE)
    call restart_output_file%put_att("time", "units", timeunits)
    call restart_output_file%put_att("time", "coordinate", time_coord)
    call restart_output_file%put_att("time", "calendar", calendar)
    call restart_output_file%end_def()
    call restart_output_file%put_var("time", [time_value])

  end subroutine cable_restart_write_time

end module
