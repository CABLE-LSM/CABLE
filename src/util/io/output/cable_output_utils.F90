module cable_output_utils_mod

  use iso_fortran_env, only: int32, real32, real64

  use cable_common_module, only: filename

  use cable_def_types_mod, only: mp
  use cable_def_types_mod, only: mp_global
  use cable_def_types_mod, only: mland
  use cable_def_types_mod, only: mland_global
  use cable_def_types_mod, only: ms
  use cable_def_types_mod, only: msn
  use cable_def_types_mod, only: nrb
  use cable_def_types_mod, only: ncs
  use cable_def_types_mod, only: ncp

  use cable_io_vars_module, only: output
  use cable_io_vars_module, only: metgrid
  use cable_io_vars_module, only: xdimsize
  use cable_io_vars_module, only: ydimsize
  use cable_io_vars_module, only: max_vegpatches
  use cable_io_vars_module, only: timeunits
  use cable_io_vars_module, only: time_coord
  use cable_io_vars_module, only: calendar

  use cable_abort_module, only: cable_abort

  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: cable_netcdf_file_t
  use cable_netcdf_mod, only: MAX_LEN_DIM => CABLE_NETCDF_MAX_STR_LEN_DIM
  use cable_netcdf_mod, only: CABLE_NETCDF_UNLIMITED
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE

  use cable_io_decomp_mod, only: io_decomp_t

  use cable_output_types_mod, only: cable_output_dim_t
  use cable_output_types_mod, only: cable_output_variable_t
  use cable_output_types_mod, only: cable_output_profile_t
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_PATCH
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SOIL
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SNOW
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_RAD
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_PLANTCARBON
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SOILCARBON
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_LAND
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_LAND_GLOBAL
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_X
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_Y
  use cable_output_types_mod, only: FILL_VALUE_INT32
  use cable_output_types_mod, only: FILL_VALUE_REAL32
  use cable_output_types_mod, only: FILL_VALUE_REAL64

  implicit none
  private

  public :: init_decomp_pointers
  public :: allocate_grid_reduction_buffers
  public :: deallocate_grid_reduction_buffers
  public :: requires_x_y_output_grid
  public :: requires_land_output_grid
  public :: check_invalid_frequency
  public :: dim_size
  public :: infer_dim_names
  public :: define_variables
  public :: set_global_attributes
  public :: associate_decomp_int32
  public :: associate_decomp_real32
  public :: associate_decomp_real64
  public :: associate_temp_buffer_int32
  public :: associate_temp_buffer_real32
  public :: associate_temp_buffer_real64

  ! Decomposition pointers for each variable class and data type
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soil_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soil_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soil_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_snow_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_snow_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_snow_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_rad_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_rad_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_rad_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_plantcarbon_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_plantcarbon_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_plantcarbon_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soilcarbon_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soilcarbon_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soilcarbon_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soil_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soil_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soil_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_snow_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_snow_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_snow_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_rad_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_rad_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_rad_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_plantcarbon_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_plantcarbon_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_plantcarbon_real64
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soilcarbon_int32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soilcarbon_real32
  class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soilcarbon_real64

  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_int32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_real32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_real64
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_soil_int32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_soil_real32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_soil_real64
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_snow_int32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_snow_real32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_snow_real64
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_rad_int32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_rad_real32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_rad_real64
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_plantcarbon_int32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_plantcarbon_real32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_plantcarbon_real64
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_soilcarbon_int32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_soilcarbon_real32
  class(cable_netcdf_decomp_t), pointer :: restart_decomp_patch_soilcarbon_real64

  ! Temporary buffers for computing grid-cell averages for each variable class
  integer(kind=int32), allocatable, target :: temp_buffer_land_int32(:)
  real(kind=real32),   allocatable, target :: temp_buffer_land_real32(:)
  real(kind=real64),   allocatable, target :: temp_buffer_land_real64(:)
  integer(kind=int32), allocatable, target :: temp_buffer_land_soil_int32(:, :)
  real(kind=real32),   allocatable, target :: temp_buffer_land_soil_real32(:, :)
  real(kind=real64),   allocatable, target :: temp_buffer_land_soil_real64(:, :)
  integer(kind=int32), allocatable, target :: temp_buffer_land_snow_int32(:, :)
  real(kind=real32),   allocatable, target :: temp_buffer_land_snow_real32(:, :)
  real(kind=real64),   allocatable, target :: temp_buffer_land_snow_real64(:, :)
  integer(kind=int32), allocatable, target :: temp_buffer_land_rad_int32(:, :)
  real(kind=real32),   allocatable, target :: temp_buffer_land_rad_real32(:, :)
  real(kind=real64),   allocatable, target :: temp_buffer_land_rad_real64(:, :)
  integer(kind=int32), allocatable, target :: temp_buffer_land_plantcarbon_int32(:, :)
  real(kind=real32),   allocatable, target :: temp_buffer_land_plantcarbon_real32(:, :)
  real(kind=real64),   allocatable, target :: temp_buffer_land_plantcarbon_real64(:, :)
  integer(kind=int32), allocatable, target :: temp_buffer_land_soilcarbon_int32(:, :)
  real(kind=real32),   allocatable, target :: temp_buffer_land_soilcarbon_real32(:, :)
  real(kind=real64),   allocatable, target :: temp_buffer_land_soilcarbon_real64(:, :)

  type(io_decomp_t) :: io_decomp_global

contains

  subroutine init_decomp_pointers(io_decomp)
    type(io_decomp_t), intent(in), target :: io_decomp

    if (requires_x_y_output_grid(output%grid, metGrid)) then
      output_decomp_base_int32                    => io_decomp%land_to_x_y_int32
      output_decomp_base_real32                   => io_decomp%land_to_x_y_real32
      output_decomp_base_real64                   => io_decomp%land_to_x_y_real64
      output_decomp_base_soil_int32               => io_decomp%land_soil_to_x_y_soil_int32
      output_decomp_base_soil_real32              => io_decomp%land_soil_to_x_y_soil_real32
      output_decomp_base_soil_real64              => io_decomp%land_soil_to_x_y_soil_real64
      output_decomp_base_snow_int32               => io_decomp%land_snow_to_x_y_snow_int32
      output_decomp_base_snow_real32              => io_decomp%land_snow_to_x_y_snow_real32
      output_decomp_base_snow_real64              => io_decomp%land_snow_to_x_y_snow_real64
      output_decomp_base_rad_int32                => io_decomp%land_rad_to_x_y_rad_int32
      output_decomp_base_rad_real32               => io_decomp%land_rad_to_x_y_rad_real32
      output_decomp_base_rad_real64               => io_decomp%land_rad_to_x_y_rad_real64
      output_decomp_base_plantcarbon_int32        => io_decomp%land_plantcarbon_to_x_y_plantcarbon_int32
      output_decomp_base_plantcarbon_real32       => io_decomp%land_plantcarbon_to_x_y_plantcarbon_real32
      output_decomp_base_plantcarbon_real64       => io_decomp%land_plantcarbon_to_x_y_plantcarbon_real64
      output_decomp_base_soilcarbon_int32         => io_decomp%land_soilcarbon_to_x_y_soilcarbon_int32
      output_decomp_base_soilcarbon_real32        => io_decomp%land_soilcarbon_to_x_y_soilcarbon_real32
      output_decomp_base_soilcarbon_real64        => io_decomp%land_soilcarbon_to_x_y_soilcarbon_real64
      output_decomp_base_patch_int32              => io_decomp%patch_to_x_y_patch_int32
      output_decomp_base_patch_real32             => io_decomp%patch_to_x_y_patch_real32
      output_decomp_base_patch_real64             => io_decomp%patch_to_x_y_patch_real64
      output_decomp_base_patch_soil_int32         => io_decomp%patch_soil_to_x_y_patch_soil_int32
      output_decomp_base_patch_soil_real32        => io_decomp%patch_soil_to_x_y_patch_soil_real32
      output_decomp_base_patch_soil_real64        => io_decomp%patch_soil_to_x_y_patch_soil_real64
      output_decomp_base_patch_snow_int32         => io_decomp%patch_snow_to_x_y_patch_snow_int32
      output_decomp_base_patch_snow_real32        => io_decomp%patch_snow_to_x_y_patch_snow_real32
      output_decomp_base_patch_snow_real64        => io_decomp%patch_snow_to_x_y_patch_snow_real64
      output_decomp_base_patch_rad_int32          => io_decomp%patch_rad_to_x_y_patch_rad_int32
      output_decomp_base_patch_rad_real32         => io_decomp%patch_rad_to_x_y_patch_rad_real32
      output_decomp_base_patch_rad_real64         => io_decomp%patch_rad_to_x_y_patch_rad_real64
      output_decomp_base_patch_plantcarbon_int32  => io_decomp%patch_plantcarbon_to_x_y_patch_plantcarbon_int32
      output_decomp_base_patch_plantcarbon_real32 => io_decomp%patch_plantcarbon_to_x_y_patch_plantcarbon_real32
      output_decomp_base_patch_plantcarbon_real64 => io_decomp%patch_plantcarbon_to_x_y_patch_plantcarbon_real64
      output_decomp_base_patch_soilcarbon_int32   => io_decomp%patch_soilcarbon_to_x_y_patch_soilcarbon_int32
      output_decomp_base_patch_soilcarbon_real32  => io_decomp%patch_soilcarbon_to_x_y_patch_soilcarbon_real32
      output_decomp_base_patch_soilcarbon_real64  => io_decomp%patch_soilcarbon_to_x_y_patch_soilcarbon_real64
    else if (requires_land_output_grid(output%grid, metGrid)) then
      output_decomp_base_int32                    => io_decomp%land_to_land_int32
      output_decomp_base_real32                   => io_decomp%land_to_land_real32
      output_decomp_base_real64                   => io_decomp%land_to_land_real64
      output_decomp_base_soil_int32               => io_decomp%land_soil_to_land_soil_int32
      output_decomp_base_soil_real32              => io_decomp%land_soil_to_land_soil_real32
      output_decomp_base_soil_real64              => io_decomp%land_soil_to_land_soil_real64
      output_decomp_base_snow_int32               => io_decomp%land_snow_to_land_snow_int32
      output_decomp_base_snow_real32              => io_decomp%land_snow_to_land_snow_real32
      output_decomp_base_snow_real64              => io_decomp%land_snow_to_land_snow_real64
      output_decomp_base_rad_int32                => io_decomp%land_rad_to_land_rad_int32
      output_decomp_base_rad_real32               => io_decomp%land_rad_to_land_rad_real32
      output_decomp_base_rad_real64               => io_decomp%land_rad_to_land_rad_real64
      output_decomp_base_plantcarbon_int32        => io_decomp%land_plantcarbon_to_land_plantcarbon_int32
      output_decomp_base_plantcarbon_real32       => io_decomp%land_plantcarbon_to_land_plantcarbon_real32
      output_decomp_base_plantcarbon_real64       => io_decomp%land_plantcarbon_to_land_plantcarbon_real64
      output_decomp_base_soilcarbon_int32         => io_decomp%land_soilcarbon_to_land_soilcarbon_int32
      output_decomp_base_soilcarbon_real32        => io_decomp%land_soilcarbon_to_land_soilcarbon_real32
      output_decomp_base_soilcarbon_real64        => io_decomp%land_soilcarbon_to_land_soilcarbon_real64
      output_decomp_base_patch_int32              => io_decomp%patch_to_land_patch_int32
      output_decomp_base_patch_real32             => io_decomp%patch_to_land_patch_real32
      output_decomp_base_patch_real64             => io_decomp%patch_to_land_patch_real64
      output_decomp_base_patch_soil_int32         => io_decomp%patch_soil_to_land_patch_soil_int32
      output_decomp_base_patch_soil_real32        => io_decomp%patch_soil_to_land_patch_soil_real32
      output_decomp_base_patch_soil_real64        => io_decomp%patch_soil_to_land_patch_soil_real64
      output_decomp_base_patch_snow_int32         => io_decomp%patch_snow_to_land_patch_snow_int32
      output_decomp_base_patch_snow_real32        => io_decomp%patch_snow_to_land_patch_snow_real32
      output_decomp_base_patch_snow_real64        => io_decomp%patch_snow_to_land_patch_snow_real64
      output_decomp_base_patch_rad_int32          => io_decomp%patch_rad_to_land_patch_rad_int32
      output_decomp_base_patch_rad_real32         => io_decomp%patch_rad_to_land_patch_rad_real32
      output_decomp_base_patch_rad_real64         => io_decomp%patch_rad_to_land_patch_rad_real64
      output_decomp_base_patch_plantcarbon_int32  => io_decomp%patch_plantcarbon_to_land_patch_plantcarbon_int32
      output_decomp_base_patch_plantcarbon_real32 => io_decomp%patch_plantcarbon_to_land_patch_plantcarbon_real32
      output_decomp_base_patch_plantcarbon_real64 => io_decomp%patch_plantcarbon_to_land_patch_plantcarbon_real64
      output_decomp_base_patch_soilcarbon_int32   => io_decomp%patch_soilcarbon_to_land_patch_soilcarbon_int32
      output_decomp_base_patch_soilcarbon_real32  => io_decomp%patch_soilcarbon_to_land_patch_soilcarbon_real32
      output_decomp_base_patch_soilcarbon_real64  => io_decomp%patch_soilcarbon_to_land_patch_soilcarbon_real64
    else
      call cable_abort("Error: Unable to determine output grid type", __FILE__, __LINE__)
    end if

    restart_decomp_patch_int32                => io_decomp%patch_to_patch_int32
    restart_decomp_patch_real32               => io_decomp%patch_to_patch_real32
    restart_decomp_patch_real64               => io_decomp%patch_to_patch_real64
    restart_decomp_patch_soil_int32           => io_decomp%patch_soil_to_patch_soil_int32
    restart_decomp_patch_soil_real32          => io_decomp%patch_soil_to_patch_soil_real32
    restart_decomp_patch_soil_real64          => io_decomp%patch_soil_to_patch_soil_real64
    restart_decomp_patch_snow_int32           => io_decomp%patch_snow_to_patch_snow_int32
    restart_decomp_patch_snow_real32          => io_decomp%patch_snow_to_patch_snow_real32
    restart_decomp_patch_snow_real64          => io_decomp%patch_snow_to_patch_snow_real64
    restart_decomp_patch_rad_int32            => io_decomp%patch_rad_to_patch_rad_int32
    restart_decomp_patch_rad_real32           => io_decomp%patch_rad_to_patch_rad_real32
    restart_decomp_patch_rad_real64           => io_decomp%patch_rad_to_patch_rad_real64
    restart_decomp_patch_plantcarbon_int32    => io_decomp%patch_plantcarbon_to_patch_plantcarbon_int32
    restart_decomp_patch_plantcarbon_real32   => io_decomp%patch_plantcarbon_to_patch_plantcarbon_real32
    restart_decomp_patch_plantcarbon_real64   => io_decomp%patch_plantcarbon_to_patch_plantcarbon_real64
    restart_decomp_patch_soilcarbon_int32     => io_decomp%patch_soilcarbon_to_patch_soilcarbon_int32
    restart_decomp_patch_soilcarbon_real32    => io_decomp%patch_soilcarbon_to_patch_soilcarbon_real32
    restart_decomp_patch_soilcarbon_real64    => io_decomp%patch_soilcarbon_to_patch_soilcarbon_real64

  end subroutine

  subroutine allocate_grid_reduction_buffers()

    allocate(temp_buffer_land_int32(mland))
    allocate(temp_buffer_land_real32(mland))
    allocate(temp_buffer_land_real64(mland))
    allocate(temp_buffer_land_soil_int32(mland, ms))
    allocate(temp_buffer_land_soil_real32(mland, ms))
    allocate(temp_buffer_land_soil_real64(mland, ms))
    allocate(temp_buffer_land_snow_int32(mland, msn))
    allocate(temp_buffer_land_snow_real32(mland, msn))
    allocate(temp_buffer_land_snow_real64(mland, msn))
    allocate(temp_buffer_land_rad_int32(mland, nrb))
    allocate(temp_buffer_land_rad_real32(mland, nrb))
    allocate(temp_buffer_land_rad_real64(mland, nrb))
    allocate(temp_buffer_land_plantcarbon_int32(mland, ncp))
    allocate(temp_buffer_land_plantcarbon_real32(mland, ncp))
    allocate(temp_buffer_land_plantcarbon_real64(mland, ncp))
    allocate(temp_buffer_land_soilcarbon_int32(mland, ncs))
    allocate(temp_buffer_land_soilcarbon_real32(mland, ncs))
    allocate(temp_buffer_land_soilcarbon_real64(mland, ncs))

  end subroutine

  subroutine deallocate_grid_reduction_buffers()

    deallocate(temp_buffer_land_int32)
    deallocate(temp_buffer_land_real32)
    deallocate(temp_buffer_land_real64)
    deallocate(temp_buffer_land_soil_int32)
    deallocate(temp_buffer_land_soil_real32)
    deallocate(temp_buffer_land_soil_real64)
    deallocate(temp_buffer_land_snow_int32)
    deallocate(temp_buffer_land_snow_real32)
    deallocate(temp_buffer_land_snow_real64)
    deallocate(temp_buffer_land_rad_int32)
    deallocate(temp_buffer_land_rad_real32)
    deallocate(temp_buffer_land_rad_real64)
    deallocate(temp_buffer_land_plantcarbon_int32)
    deallocate(temp_buffer_land_plantcarbon_real32)
    deallocate(temp_buffer_land_plantcarbon_real64)
    deallocate(temp_buffer_land_soilcarbon_int32)
    deallocate(temp_buffer_land_soilcarbon_real32)
    deallocate(temp_buffer_land_soilcarbon_real64)

  end subroutine

  logical function requires_x_y_output_grid(output_grid, met_grid)
    character(len=*), intent(in) :: output_grid
    character(len=*), intent(in) :: met_grid
    requires_x_y_output_grid = (( &
      output_grid == "default" .AND. met_grid == "mask" &
    ) .OR. ( &
      output_grid == "mask" .OR. output_grid == "ALMA" &
    ))
  end function

  logical function requires_land_output_grid(output_grid, met_grid)
    character(len=*), intent(in) :: output_grid
    character(len=*), intent(in) :: met_grid
    requires_land_output_grid = ( &
      output_grid == "land" .OR. (output_grid == "default" .AND. met_grid == "land") &
    )
  end function

  elemental integer function dim_size(dim)
    type(cable_output_dim_t), intent(in) :: dim

    select case (dim%value)
    case (CABLE_OUTPUT_DIM_PATCH%value)
      dim_size = mp
    case (CABLE_OUTPUT_DIM_SOIL%value)
      dim_size = ms
    case (CABLE_OUTPUT_DIM_SNOW%value)
      dim_size = msn
    case (CABLE_OUTPUT_DIM_RAD%value)
      dim_size = nrb
    case (CABLE_OUTPUT_DIM_PLANTCARBON%value)
      dim_size = ncp
    case (CABLE_OUTPUT_DIM_SOILCARBON%value)
      dim_size = ncs
    case (CABLE_OUTPUT_DIM_LAND%value)
      dim_size = mland
    case (CABLE_OUTPUT_DIM_LAND_GLOBAL%value)
      dim_size = mland_global
    case (CABLE_OUTPUT_DIM_X%value)
      dim_size = xdimsize
    case (CABLE_OUTPUT_DIM_Y%value)
      dim_size = ydimsize
    case default
      dim_size = -1 ! Unknown dimension
    end select

  end function

  subroutine check_invalid_frequency(sampling_frequency, accumulation_frequency, var_name, file_name)
    character(len=*), intent(in) :: sampling_frequency
    character(len=*), intent(in) :: accumulation_frequency
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in) :: file_name

    integer :: sampling_period_in_hours, accumulation_period_in_hours

    character(len=256) :: err_message

    err_message = ( &
      "Invalid combination of sampling frequency '" // sampling_frequency // &
      "' with accumulation frequency '" // accumulation_frequency // "' for variable '" // &
      var_name // "' in file '" // file_name // "'" &
    )

    select case (sampling_frequency)
    case ("all")
      if (accumulation_frequency /= "all") call cable_abort(err_message, __FILE__, __LINE__)
    case ("user")
      read(sampling_frequency(5:7), *) sampling_period_in_hours
      if (accumulation_frequency == "user") then
        read(accumulation_frequency(5:7), *) accumulation_period_in_hours
        if (sampling_period_in_hours < accumulation_period_in_hours) then
          call cable_abort(err_message, __FILE__, __LINE__)
        end if
      else if (accumulation_frequency /= "all") then
        call cable_abort(err_message, __FILE__, __LINE__)
      end if
    case ("daily")
      if (.not. any(accumulation_frequency == ["all", "daily", "user"])) then
        call cable_abort(err_message, __FILE__, __LINE__)
      end if
    case ("monthly")
      if (.not. any(accumulation_frequency == ["all", "daily", "user", "monthly"])) then
        call cable_abort(err_message, __FILE__, __LINE__)
      end if
    case default
      call cable_abort("Invalid sampling frequency '" // sampling_frequency // &
        "' for variable '" // var_name // "' in file '" // file_name // "'", __FILE__, __LINE__)
    end select

  end subroutine check_invalid_frequency

  function infer_dim_names(output_variable, restart) result(dim_names)
    type(cable_output_variable_t), intent(in) :: output_variable
    logical, intent(in), optional :: restart

    character(MAX_LEN_DIM), allocatable :: dim_names(:)
    logical :: restart_local
    integer :: j

    restart_local = .false.
    if (present(restart)) restart_local = restart

    allocate(dim_names(0))
    do j = 1, size(output_variable%data_shape)
      if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_PATCH) then
        if (restart_local) then
          dim_names = [dim_names, "mp"]
        else if (requires_land_output_grid(output%grid, metgrid)) then
          if (output_variable%reduction_method == "none") then
            dim_names = [dim_names, "land", "patch"]
          else
            dim_names = [dim_names, "land"]
          end if
        else if (requires_x_y_output_grid(output%grid, metgrid)) then
          if (output_variable%reduction_method == "none") then
            dim_names = [dim_names, "x", "y", "patch"]
          else
            dim_names = [dim_names, "x", "y"]
          end if
        end if
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_LAND) then
        if (restart_local) then
          dim_names = [dim_names, "mland"]
        else if (requires_land_output_grid(output%grid, metgrid)) then
          dim_names = [dim_names, "land"]
        else if (requires_x_y_output_grid(output%grid, metgrid)) then
          dim_names = [dim_names, "x", "y"]
        end if
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_LAND_GLOBAL) then
        if (restart_local) then
          dim_names = [dim_names, "mland"]
        else
          dim_names = [dim_names, "land"]
        end if
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_SOIL) then
        dim_names = [dim_names, "soil"]
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_SNOW) then
        dim_names = [dim_names, "snow"]
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_RAD) then
        dim_names = [dim_names, "rad"]
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_PLANTCARBON) then
        if (restart_local) then
          dim_names = [dim_names, "plant_carbon_pools"]
        else
          dim_names = [dim_names, "plantcarbon"]
        end if
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_SOILCARBON) then
        if (restart_local) then
          dim_names = [dim_names, "soil_carbon_pools"]
        else
          dim_names = [dim_names, "soilcarbon"]
        end if
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_X) then
        dim_names = [dim_names, "x"]
      else if (output_variable%data_shape(j) == CABLE_OUTPUT_DIM_Y) then
        dim_names = [dim_names, "y"]
      else
        call cable_abort("Unexpected data shape for variable " // output_variable%name, __FILE__, __LINE__)
      end if
    end do

    if (.not. restart_local .and. .not. output_variable%parameter) dim_names = [dim_names, "time"]

  end function

  function infer_cell_methods(output_variable) result(cell_methods)
    type(cable_output_variable_t), intent(in) :: output_variable
    character(len=256) :: cell_methods

    if (.not. output_variable%parameter) then
      select case (output_variable%aggregation_method)
      case ("point")
        cell_methods = "time: point"
      case ("mean")
        cell_methods = "time: mean"
      case ("sum")
        cell_methods = "time: sum"
      case ("min")
        cell_methods = "time: minimum"
      case ("max")
        cell_methods = "time: maximum"
      case default
        call cable_abort("Unexpected aggregation method '" // output_variable%aggregation_method // &
          "' for variable '" // output_variable%name // "'", __FILE__, __LINE__)
      end select
    end if

    select case (output_variable%reduction_method)
    case ("none", "first_patch_in_grid_cell")
      ! no additional cell methods
    case ("grid_cell_average")
      if (len_trim(cell_methods) > 0) then
        cell_methods = cell_methods // " area: mean"
      else
        cell_methods = "area: mean"
      end if
    case default
      call cable_abort("Unexpected reduction method '" // output_variable%reduction_method // &
        "' for variable '" // output_variable%name // "'", __FILE__, __LINE__)
    end select

  end function

  subroutine define_variables(output_file, output_variables, restart)
    class(cable_netcdf_file_t), intent(inout) :: output_file
    type(cable_output_variable_t), intent(in) :: output_variables(:)
    logical, intent(in), optional :: restart

    integer :: i, j
    logical :: restart_local

    character(MAX_LEN_DIM), allocatable :: required_dimensions(:), dim_names(:)

    restart_local = .false.
    if (present(restart)) restart_local = restart

    do i = 1, size(output_variables)
      associate(output_var => output_variables(i))
        if (restart_local .and. .not. output_var%restart) cycle
        if (.not. allocated(output_var%data_shape)) cycle
        dim_names = infer_dim_names(output_var, restart_local)
        if (.not. allocated(required_dimensions)) then
          required_dimensions = dim_names
        else
          required_dimensions = [ &
            required_dimensions, &
            pack(dim_names, [( &
              .not. any(dim_names(j) == required_dimensions), &
              j = 1, &
              size(dim_names) &
            )]) &
          ]
        end if
      end associate
    end do

    do i = 1, size(required_dimensions)
      select case (required_dimensions(i))
      case ("mp")
        call output_file%def_dims(["mp"], [mp_global])
      case ("mland")
        call output_file%def_dims(["mland"], [mland_global])
      case ("land")
        call output_file%def_dims(["land"], [mland_global])
      case ("x")
        call output_file%def_dims(["x"], [xdimsize])
      case ("y")
        call output_file%def_dims(["y"], [ydimsize])
      case ("patch")
        call output_file%def_dims(["patch"], [max_vegpatches])
      case ("soil")
        call output_file%def_dims(["soil"], [ms])
      case ("rad")
        call output_file%def_dims(["rad"], [nrb])
      case ("soil_carbon_pools")
        call output_file%def_dims(["soil_carbon_pools"], [ncs])
      case ("soilcarbon")
        call output_file%def_dims(["soilcarbon"], [ncs])
      case ("plant_carbon_pools")
        call output_file%def_dims(["plant_carbon_pools"], [ncp])
      case ("plantcarbon")
        call output_file%def_dims(["plantcarbon"], [ncp])
      case ("time")
        ! time dimension defined separately below
      case default
        call cable_abort("Unexpected dimension name '" // required_dimensions(i) // "'", __FILE__, __LINE__)
      end select
    end do

    if (restart_local) then
      call output_file%def_dims(["time"], [1])
    else
      call output_file%def_dims(["time"], [CABLE_NETCDF_UNLIMITED])
    end if

    call output_file%def_var("time", ["time"], CABLE_NETCDF_DOUBLE)
    call output_file%put_att("time", "units", timeunits)
    call output_file%put_att("time", "coordinate", time_coord)
    call output_file%put_att("time", "calendar", calendar)

    if (.not. restart_local) then
      call output_file%def_dims(["nv"], [2])
      call output_file%def_var("time_bnds", ["nv", "time"], CABLE_NETCDF_DOUBLE)
      call output_file%put_att("time", "bounds", "time_bnds")
    end if

    do i = 1, size(output_variables)
      associate(output_var => output_variables(i))
        if (restart_local .and. .not. output_var%restart) cycle
        call output_file%def_var( &
          var_name=output_var%name, &
          dim_names=infer_dim_names(output_var, restart_local), &
          type=output_var%var_type &
        )
        call output_file%put_att(output_var%name, "units", output_var%units)
        call output_file%put_att(output_var%name, "long_name", output_var%long_name)
        select case (output_var%var_type)
        case (CABLE_NETCDF_INT)
          call output_file%put_att(output_var%name, "_FillValue", FILL_VALUE_INT32)
          call output_file%put_att(output_var%name, "missing_value", FILL_VALUE_INT32)
        case (CABLE_NETCDF_FLOAT)
          call output_file%put_att(output_var%name, "_FillValue", FILL_VALUE_REAL32)
          call output_file%put_att(output_var%name, "missing_value", FILL_VALUE_REAL32)
        case (CABLE_NETCDF_DOUBLE)
          call output_file%put_att(output_var%name, "_FillValue", FILL_VALUE_REAL64)
          call output_file%put_att(output_var%name, "missing_value", FILL_VALUE_REAL64)
        end select
        call output_file%put_att(output_var%name, "cell_methods", infer_cell_methods(output_var))
      end associate
    end do

  end subroutine define_variables

  subroutine set_global_attributes(output_profile)
    type(cable_output_profile_t), intent(inout) :: output_profile

    character(32) :: todaydate, nowtime

    call date_and_time(todaydate, nowtime)
    todaydate = todaydate(1:4) // "/" // todaydate(5:6) // "/" // todaydate(7:8)
    nowtime = nowtime(1:2) // ":" // nowtime(3:4) // ":" // nowtime(5:6)
    call output_profile%output_file%put_att("Production", trim(todaydate) // " at " // trim(nowtime))
    call output_profile%output_file%put_att("Source", "CABLE LSM output file")
    call output_profile%output_file%put_att("CABLE_input_file", trim(filename%met))

    select case (output_profile%sampling_frequency)
    case ("user")
       call output_profile%output_file%put_att("Output_averaging", TRIM(output_profile%sampling_frequency(5:7)) // "-hourly output")
    case ("all")
       call output_profile%output_file%put_att("Output_averaging", "all timesteps recorded")
    case ("daily")
       call output_profile%output_file%put_att("Output_averaging", "daily")
    case ("monthly")
       call output_profile%output_file%put_att("Output_averaging", "monthly")
    case default
       call cable_abort("Invalid sampling frequency '" // output_profile%sampling_frequency // "'", __FILE__, __LINE__)
    end select

  end subroutine set_global_attributes

  subroutine associate_decomp_int32(output_var, decomp, restart)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    logical, intent(in), optional :: restart

    if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_int32
      else
        decomp => output_decomp_base_int32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_int32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_soil_int32
      else
        decomp => output_decomp_base_soil_int32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_soil_int32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_snow_int32
      else
        decomp => output_decomp_base_snow_int32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_snow_int32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_rad_int32
      else
        decomp => output_decomp_base_rad_int32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_rad_int32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_plantcarbon_int32
      else
        decomp => output_decomp_base_plantcarbon_int32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_plantcarbon_int32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_soilcarbon_int32
      else
        decomp => output_decomp_base_soilcarbon_int32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_soilcarbon_int32
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_int32

  subroutine associate_decomp_real32(output_var, decomp, restart)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    logical, intent(in), optional :: restart

    if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_real32
      else
        decomp => output_decomp_base_real32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_real32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_soil_real32
      else
        decomp => output_decomp_base_soil_real32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_soil_real32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_snow_real32
      else
        decomp => output_decomp_base_snow_real32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_snow_real32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_rad_real32
      else
        decomp => output_decomp_base_rad_real32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_rad_real32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_plantcarbon_real32
      else
        decomp => output_decomp_base_plantcarbon_real32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_plantcarbon_real32
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_soilcarbon_real32
      else
        decomp => output_decomp_base_soilcarbon_real32
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_soilcarbon_real32
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_real32

  subroutine associate_decomp_real64(output_var, decomp, restart)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    logical, intent(in), optional :: restart

    if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_real64
      else
        decomp => output_decomp_base_real64
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_real64
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_soil_real64
      else
        decomp => output_decomp_base_soil_real64
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_soil_real64
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_snow_real64
      else
        decomp => output_decomp_base_snow_real64
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_snow_real64
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_rad_real64
      else
        decomp => output_decomp_base_rad_real64
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_rad_real64
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_plantcarbon_real64
      else
        decomp => output_decomp_base_plantcarbon_real64
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_plantcarbon_real64
      end if
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp_base_patch_soilcarbon_real64
      else
        decomp => output_decomp_base_soilcarbon_real64
      end if
      if (present(restart)) then
        if (restart) decomp => restart_decomp_patch_soilcarbon_real64
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_real64

  subroutine associate_temp_buffer_int32(output_var, temp_buffer_int32_1d, temp_buffer_int32_2d, temp_buffer_int32_3d)
    type(cable_output_variable_t), intent(inout) :: output_var
    integer(kind=int32), pointer, intent(inout), optional :: temp_buffer_int32_1d(:)
    integer(kind=int32), pointer, intent(inout), optional :: temp_buffer_int32_2d(:,:)
    integer(kind=int32), pointer, intent(inout), optional :: temp_buffer_int32_3d(:,:,:)

    if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH])) then
      if (.not. present(temp_buffer_int32_1d)) call cable_abort( &
        "temp_buffer_int32_1d must be provided for 1D data shape", __FILE__, __LINE__)
      temp_buffer_int32_1d => temp_buffer_land_int32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (.not. present(temp_buffer_int32_2d)) call cable_abort( &
        "temp_buffer_int32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_int32_2d => temp_buffer_land_soil_int32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (.not. present(temp_buffer_int32_2d)) call cable_abort( &
        "temp_buffer_int32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_int32_2d => temp_buffer_land_rad_int32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (.not. present(temp_buffer_int32_2d)) call cable_abort( &
        "temp_buffer_int32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_int32_2d => temp_buffer_land_snow_int32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (.not. present(temp_buffer_int32_2d)) call cable_abort( &
        "temp_buffer_int32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_int32_2d => temp_buffer_land_rad_int32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (.not. present(temp_buffer_int32_2d)) call cable_abort( &
        "temp_buffer_int32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_int32_2d => temp_buffer_land_plantcarbon_int32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (.not. present(temp_buffer_int32_2d)) call cable_abort( &
        "temp_buffer_int32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_int32_2d => temp_buffer_land_soilcarbon_int32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_int32

  subroutine associate_temp_buffer_real32(output_var, temp_buffer_real32_1d, temp_buffer_real32_2d, temp_buffer_real32_3d)
    type(cable_output_variable_t), intent(inout) :: output_var
    real(kind=real32), pointer, intent(inout), optional :: temp_buffer_real32_1d(:)
    real(kind=real32), pointer, intent(inout), optional :: temp_buffer_real32_2d(:,:)
    real(kind=real32), pointer, intent(inout), optional :: temp_buffer_real32_3d(:,:,:)

    if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH])) then
      if (.not. present(temp_buffer_real32_1d)) call cable_abort( &
        "temp_buffer_real32_1d must be provided for 1D data shape", __FILE__, __LINE__)
      temp_buffer_real32_1d => temp_buffer_land_real32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (.not. present(temp_buffer_real32_2d)) call cable_abort( &
        "temp_buffer_real32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real32_2d => temp_buffer_land_soil_real32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (.not. present(temp_buffer_real32_2d)) call cable_abort( &
        "temp_buffer_real32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real32_2d => temp_buffer_land_rad_real32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (.not. present(temp_buffer_real32_2d)) call cable_abort( &
        "temp_buffer_real32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real32_2d => temp_buffer_land_snow_real32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (.not. present(temp_buffer_real32_2d)) call cable_abort( &
        "temp_buffer_real32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real32_2d => temp_buffer_land_rad_real32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (.not. present(temp_buffer_real32_2d)) call cable_abort( &
        "temp_buffer_real32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real32_2d => temp_buffer_land_plantcarbon_real32
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (.not. present(temp_buffer_real32_2d)) call cable_abort( &
        "temp_buffer_real32_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real32_2d => temp_buffer_land_soilcarbon_real32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_real32

  subroutine associate_temp_buffer_real64(output_var, temp_buffer_real64_1d, temp_buffer_real64_2d, temp_buffer_real64_3d)
    type(cable_output_variable_t), intent(inout) :: output_var
    real(kind=real64), pointer, intent(inout), optional :: temp_buffer_real64_1d(:)
    real(kind=real64), pointer, intent(inout), optional :: temp_buffer_real64_2d(:,:)
    real(kind=real64), pointer, intent(inout), optional :: temp_buffer_real64_3d(:,:,:)

    if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH])) then
      if (.not. present(temp_buffer_real64_1d)) call cable_abort( &
        "temp_buffer_real64_1d must be provided for 1D data shape", __FILE__, __LINE__)
      temp_buffer_real64_1d => temp_buffer_land_real64
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (.not. present(temp_buffer_real64_2d)) call cable_abort( &
        "temp_buffer_real64_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real64_2d => temp_buffer_land_soil_real64
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (.not. present(temp_buffer_real64_2d)) call cable_abort( &
        "temp_buffer_real64_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real64_2d => temp_buffer_land_rad_real64
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (.not. present(temp_buffer_real64_2d)) call cable_abort( &
        "temp_buffer_real64_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real64_2d => temp_buffer_land_snow_real64
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (.not. present(temp_buffer_real64_2d)) call cable_abort( &
        "temp_buffer_real64_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real64_2d => temp_buffer_land_rad_real64
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (.not. present(temp_buffer_real64_2d)) call cable_abort( &
        "temp_buffer_real64_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real64_2d => temp_buffer_land_plantcarbon_real64
    else if (all(output_var%data_shape == [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (.not. present(temp_buffer_real64_2d)) call cable_abort( &
        "temp_buffer_real64_2d must be provided for 2D data shape", __FILE__, __LINE__)
      temp_buffer_real64_2d => temp_buffer_land_soilcarbon_real64
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_real64

end module
