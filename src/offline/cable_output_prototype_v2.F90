module cable_output_prototype_v2_mod
  use iso_fortran_env, only: int32, real32, real64

  use cable_def_types_mod, only: mp, mp_global
  use cable_def_types_mod, only: mland
  use cable_def_types_mod, only: ms
  use cable_def_types_mod, only: nrb
  use cable_def_types_mod, only: ncs
  use cable_def_types_mod, only: ncp

  use cable_abort_module, only: cable_abort

  use cable_io_vars_module, only: metGrid, patch_type, land_type, xdimsize, ydimsize, max_vegpatches
  use cable_io_vars_module, only: timeunits, calendar, time_coord

  use cable_io_vars_module, only: output_options => output, patchout_options => patchout

  use cable_io_decomp_mod, only: io_decomp_t

  use cable_timing_utils_mod, only: time_step_matches

  use aggregator_mod, only: aggregator_mod_init
  use aggregator_mod, only: aggregator_mod_end
  use aggregator_mod, only: aggregator_t
  use aggregator_mod, only: aggregator_handle_t
  use aggregator_mod, only: aggregator_int32_1d_t
  use aggregator_mod, only: aggregator_int32_2d_t
  use aggregator_mod, only: aggregator_int32_3d_t
  use aggregator_mod, only: aggregator_real32_1d_t
  use aggregator_mod, only: aggregator_real32_2d_t
  use aggregator_mod, only: aggregator_real32_3d_t
  use aggregator_mod, only: aggregator_real64_1d_t
  use aggregator_mod, only: aggregator_real64_2d_t
  use aggregator_mod, only: aggregator_real64_3d_t
  use aggregator_mod, only: store_aggregator

  use cable_netcdf_mod, only: cable_netcdf_file_t
  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: cable_netcdf_create_file
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE
  use cable_netcdf_mod, only: CABLE_NETCDF_UNLIMITED
  use cable_netcdf_mod, only: MAX_LEN_VAR => CABLE_NETCDF_MAX_STR_LEN_VAR
  use cable_netcdf_mod, only: MAX_LEN_DIM => CABLE_NETCDF_MAX_STR_LEN_DIM

  use cable_output_utils_mod

  implicit none
  private

  public :: cable_output_mod_init
  public :: cable_output_mod_end
  public :: cable_output_add_variable
  public :: cable_output_aggregator_t
  public :: cable_output_add_aggregator
  public :: cable_output_commit
  public :: cable_output_update
  public :: output_options
  public :: patchout_options
  public :: requires_x_y_output_grid
  public :: requires_land_output_grid

  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_UNDEFINED              = 0
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE                   = 1
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_SOIL              = 2
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_SNOW              = 3
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_RAD               = 4
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_PLANTCARBON       = 5
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_SOILCARBON        = 6
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH             = 7
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SOIL        = 8
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SNOW        = 9
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_RAD         = 10
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_PLANTCARBON = 11
  integer, parameter, public :: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SOILCARBON  = 12

  integer(kind=int32), parameter :: FILL_VALUE_INT32  = -9999_int32
  real(kind=real32),   parameter :: FILL_VALUE_REAL32 = -1.0e+33_real32
  real(kind=real64),   parameter :: FILL_VALUE_REAL64 = -1.0e+33_real64

  type :: cable_output_aggregator_t
    type(aggregator_handle_t) :: aggregator_handle
    character(len=20) :: accumulation_frequency
    character(len=20) :: aggregation_frequency
  end type

  type cable_output_variable_t
    character(len=MAX_LEN_VAR) :: name
    character(len=MAX_LEN_DIM), allocatable :: dims(:)
    integer :: var_type
    character(len=50) :: units
    character(len=100) :: long_name
    character(len=100) :: cell_methods
    logical :: active
    logical :: grid_cell_averaging
    integer :: shape_type = CABLE_OUTPUT_SHAPE_TYPE_UNDEFINED
    real, dimension(2) :: range
    type(cable_output_aggregator_t) :: output_aggregator
    class(cable_netcdf_decomp_t), pointer :: decomp => null()
    real(kind=real32), pointer :: temp_buffer_real32_1d(:)       => null()
    real(kind=real32), pointer :: temp_buffer_real32_2d(:, :)    => null()
    real(kind=real32), pointer :: temp_buffer_real32_3d(:, :, :) => null()
    real(kind=real64), pointer :: temp_buffer_real64_1d(:)       => null()
    real(kind=real64), pointer :: temp_buffer_real64_2d(:, :)    => null()
    real(kind=real64), pointer :: temp_buffer_real64_3d(:, :, :) => null()
  end type

  type cable_output_profile_t
    real :: previous_write_time = 0.0
    integer :: frame = 0
    class(cable_netcdf_file_t), allocatable :: output_file
    !> List of output aggregators sorted in decreasing accumulation_frequency,
    ! then aggregation_frequency. Sorting the aggregators this way ensures that
    ! intermediate aggregators are updated before any aggregators which may be
    ! dependent on them.
    type(cable_output_aggregator_t), allocatable :: output_aggregators(:)
    type(cable_output_variable_t), allocatable :: output_variables(:)
  end type

  ! Decomposition mappings for each variable class and type
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

  ! Temporary buffers for computing grid-cell averages for each variable class
  real(kind=real32), allocatable, target :: temp_buffer_land_real32(:)
  real(kind=real64), allocatable, target :: temp_buffer_land_real64(:)
  real(kind=real32), allocatable, target :: temp_buffer_land_soil_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_soil_real64(:, :)
  real(kind=real32), allocatable, target :: temp_buffer_land_snow_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_snow_real64(:, :)
  real(kind=real32), allocatable, target :: temp_buffer_land_rad_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_rad_real64(:, :)
  real(kind=real32), allocatable, target :: temp_buffer_land_plantcarbon_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_plantcarbon_real64(:, :)
  real(kind=real32), allocatable, target :: temp_buffer_land_soilcarbon_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_soilcarbon_real64(:, :)

  ! TODO(Sean): once cable_write.F90 is removed, move the output_inclusion_type
  ! from cable_io_vars_module to here (as this would no longer introduce a cyclic
  ! module dependency). Then uncomment declarations below:
  ! type(output_inclusion_t) :: output_options
  ! type(output_inclusion_t) :: patchout_options ! do we want patch-specific info

  type(cable_output_profile_t), allocatable :: global_profile

contains

  logical function requires_x_y_output_grid(output_grid, met_grid)
    character(len=*), intent(in) :: output_grid
    character(len=*), intent(in) :: met_grid
    requires_x_y_output_grid = (( &
      output_grid == 'default' .AND. met_grid == 'mask' &
    ) .OR. ( &
      output_grid == 'mask' .OR. output_grid == 'ALMA' &
    ))
  end function

  logical function requires_land_output_grid(output_grid, met_grid)
    character(len=*), intent(in) :: output_grid
    character(len=*), intent(in) :: met_grid
    requires_land_output_grid = ( &
      output_grid == 'land' .OR. (output_grid == 'default' .AND. met_grid == 'land') &
    )
  end function

  function compare_aggregators_by_frequency(a, b) result(is_less)
    type(cable_output_aggregator_t), intent(in) :: a, b
    logical :: is_less

    ! TODO(Sean): sort frequency by decreasing accumulation_frequency first, then decreasing aggregation_frequency

    is_less = .false.

  end function

  subroutine sort_aggregators_by_frequency(output_aggregators)
    type(cable_output_aggregator_t), intent(inout) :: output_aggregators(:)
    integer :: i, j
    type(cable_output_aggregator_t) :: temp

    do i = 1, size(output_aggregators) - 1
      do j = i + 1, size(output_aggregators)
        if (compare_aggregators_by_frequency(output_aggregators(i), output_aggregators(j))) then
          temp = output_aggregators(i)
          output_aggregators(i) = output_aggregators(j)
          output_aggregators(j) = temp
        end if
      end do
    end do

  end subroutine

  subroutine cable_output_mod_init(io_decomp)
    type(io_decomp_t), intent(in), target :: io_decomp
    class(cable_netcdf_file_t), allocatable :: output_file

    if (requires_x_y_output_grid(output_options%grid, metGrid)) then
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
    else if (requires_land_output_grid(output_options%grid, metGrid)) then
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
      call cable_abort("Unable to determine output I/O decomposition", __FILE__, __LINE__)
    end if

    ! Initialize temporary buffers for grid-cell averaging
    allocate(temp_buffer_land_real32(mland))
    allocate(temp_buffer_land_real64(mland))
    allocate(temp_buffer_land_soil_real32(mland, mp))
    allocate(temp_buffer_land_soil_real64(mland, mp))
    allocate(temp_buffer_land_snow_real32(mland, mp))
    allocate(temp_buffer_land_snow_real64(mland, mp))
    allocate(temp_buffer_land_rad_real32(mland, mp))
    allocate(temp_buffer_land_rad_real64(mland, mp))
    allocate(temp_buffer_land_plantcarbon_real32(mland, mp))
    allocate(temp_buffer_land_plantcarbon_real64(mland, mp))
    allocate(temp_buffer_land_soilcarbon_real32(mland, mp))
    allocate(temp_buffer_land_soilcarbon_real64(mland, mp))

    call aggregator_mod_init()

    allocate(cable_output_profile_t::global_profile)

  end subroutine

  subroutine cable_output_mod_end()

    if (allocated(global_profile%output_file)) call global_profile%output_file%close()

    deallocate(global_profile)

    call aggregator_mod_end()

    if (associated(output_decomp_base_int32))                    nullify(output_decomp_base_int32)
    if (associated(output_decomp_base_real32))                   nullify(output_decomp_base_real32)
    if (associated(output_decomp_base_real64))                   nullify(output_decomp_base_real64)
    if (associated(output_decomp_base_soil_int32))               nullify(output_decomp_base_soil_int32)
    if (associated(output_decomp_base_soil_real32))              nullify(output_decomp_base_soil_real32)
    if (associated(output_decomp_base_soil_real64))              nullify(output_decomp_base_soil_real64)
    if (associated(output_decomp_base_snow_int32))               nullify(output_decomp_base_snow_int32)
    if (associated(output_decomp_base_snow_real32))              nullify(output_decomp_base_snow_real32)
    if (associated(output_decomp_base_snow_real64))              nullify(output_decomp_base_snow_real64)
    if (associated(output_decomp_base_rad_int32))                nullify(output_decomp_base_rad_int32)
    if (associated(output_decomp_base_rad_real32))               nullify(output_decomp_base_rad_real32)
    if (associated(output_decomp_base_rad_real64))               nullify(output_decomp_base_rad_real64)
    if (associated(output_decomp_base_plantcarbon_int32))        nullify(output_decomp_base_plantcarbon_int32)
    if (associated(output_decomp_base_plantcarbon_real32))       nullify(output_decomp_base_plantcarbon_real32)
    if (associated(output_decomp_base_plantcarbon_real64))       nullify(output_decomp_base_plantcarbon_real64)
    if (associated(output_decomp_base_soilcarbon_int32))         nullify(output_decomp_base_soilcarbon_int32)
    if (associated(output_decomp_base_soilcarbon_real32))        nullify(output_decomp_base_soilcarbon_real32)
    if (associated(output_decomp_base_soilcarbon_real64))        nullify(output_decomp_base_soilcarbon_real64)
    if (associated(output_decomp_base_patch_int32))              nullify(output_decomp_base_patch_int32)
    if (associated(output_decomp_base_patch_real32))             nullify(output_decomp_base_patch_real32)
    if (associated(output_decomp_base_patch_real64))             nullify(output_decomp_base_patch_real64)
    if (associated(output_decomp_base_patch_soil_int32))         nullify(output_decomp_base_patch_soil_int32)
    if (associated(output_decomp_base_patch_soil_real32))        nullify(output_decomp_base_patch_soil_real32)
    if (associated(output_decomp_base_patch_soil_real64))        nullify(output_decomp_base_patch_soil_real64)
    if (associated(output_decomp_base_patch_snow_int32))         nullify(output_decomp_base_patch_snow_int32)
    if (associated(output_decomp_base_patch_snow_real32))        nullify(output_decomp_base_patch_snow_real32)
    if (associated(output_decomp_base_patch_snow_real64))        nullify(output_decomp_base_patch_snow_real64)
    if (associated(output_decomp_base_patch_rad_int32))          nullify(output_decomp_base_patch_rad_int32)
    if (associated(output_decomp_base_patch_rad_real32))         nullify(output_decomp_base_patch_rad_real32)
    if (associated(output_decomp_base_patch_rad_real64))         nullify(output_decomp_base_patch_rad_real64)
    if (associated(output_decomp_base_patch_plantcarbon_int32))  nullify(output_decomp_base_patch_plantcarbon_int32)
    if (associated(output_decomp_base_patch_plantcarbon_real32)) nullify(output_decomp_base_patch_plantcarbon_real32)
    if (associated(output_decomp_base_patch_plantcarbon_real64)) nullify(output_decomp_base_patch_plantcarbon_real64)
    if (associated(output_decomp_base_patch_soilcarbon_int32))   nullify(output_decomp_base_patch_soilcarbon_int32)
    if (associated(output_decomp_base_patch_soilcarbon_real32))  nullify(output_decomp_base_patch_soilcarbon_real32)
    if (associated(output_decomp_base_patch_soilcarbon_real64))  nullify(output_decomp_base_patch_soilcarbon_real64)

    deallocate(temp_buffer_land_real32)
    deallocate(temp_buffer_land_real64)
    deallocate(temp_buffer_land_soil_real32)
    deallocate(temp_buffer_land_soil_real64)
    deallocate(temp_buffer_land_snow_real32)
    deallocate(temp_buffer_land_snow_real64)
    deallocate(temp_buffer_land_rad_real32)
    deallocate(temp_buffer_land_rad_real64)
    deallocate(temp_buffer_land_plantcarbon_real32)
    deallocate(temp_buffer_land_plantcarbon_real64)
    deallocate(temp_buffer_land_soilcarbon_real32)
    deallocate(temp_buffer_land_soilcarbon_real64)

  end subroutine

  subroutine cable_output_add_variable( &
    name, dims, var_type, units, long_name, active, grid_cell_averaging, &
    shape_type, range, accumulation_frequency, aggregation_frequency, aggregator &
  )
    character(len=*), intent(in) :: name
    character(len=*), dimension(:), intent(in) :: dims
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: long_name
    logical, intent(in) :: active
    logical, intent(in) :: grid_cell_averaging
    integer, intent(in) :: shape_type
    real, dimension(2), intent(in) :: range
    character(len=*), intent(in) :: accumulation_frequency
    character(len=*), intent(in) :: aggregation_frequency
    class(aggregator_t), intent(in) :: aggregator

    type(cable_output_variable_t) :: output_var

    if (grid_cell_averaging) then
      select type (aggregator)
      type is (aggregator_real32_1d_t)
        if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid cell averaging", __FILE__, __LINE__)
      type is (aggregator_real32_2d_t)
        if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid cell averaging", __FILE__, __LINE__)
      type is (aggregator_real32_3d_t)
        if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid cell averaging", __FILE__, __LINE__)
      type is (aggregator_real64_1d_t)
        if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid cell averaging", __FILE__, __LINE__)
      type is (aggregator_real64_2d_t)
        if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid cell averaging", __FILE__, __LINE__)
      type is (aggregator_real64_3d_t)
        if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid cell averaging", __FILE__, __LINE__)
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    end if

    ! TODO(Sean): determine cell_methods based on grid_cell_averaging and aggregator method

    output_var%name = trim(adjustl(name))
    output_var%dims = dims
    output_var%units = trim(adjustl(units))
    output_var%long_name = trim(adjustl(long_name))
    output_var%active = active
    output_var%grid_cell_averaging = grid_cell_averaging
    output_var%range = range
    output_var%shape_type = shape_type
    output_var%var_type = var_type

    if (active) then
      call cable_output_add_aggregator( &
        aggregator=aggregator, &
        accumulation_frequency=accumulation_frequency, &
        aggregation_frequency=aggregation_frequency, &
        output_aggregator=output_var%output_aggregator &
      )
    end if

    if (grid_cell_averaging) then
      select case(shape_type)
      case (CABLE_OUTPUT_SHAPE_TYPE_BASE)
        select type(aggregator)
        type is (aggregator_real32_1d_t)
          output_var%temp_buffer_real32_1d => temp_buffer_land_real32
        type is (aggregator_real64_1d_t)
          output_var%temp_buffer_real64_1d => temp_buffer_land_real64
        class default
          call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
        end select
      case (CABLE_OUTPUT_SHAPE_TYPE_BASE_SOIL)
        select type(aggregator)
        type is (aggregator_real32_2d_t)
          output_var%temp_buffer_real32_2d => temp_buffer_land_soil_real32
        type is (aggregator_real64_2d_t)
          output_var%temp_buffer_real64_2d => temp_buffer_land_soil_real64
        class default
          call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
        end select
      case (CABLE_OUTPUT_SHAPE_TYPE_BASE_SNOW)
        select type(aggregator)
        type is (aggregator_real32_2d_t)
          output_var%temp_buffer_real32_2d => temp_buffer_land_snow_real32
        type is (aggregator_real64_2d_t)
          output_var%temp_buffer_real64_2d => temp_buffer_land_snow_real64
        class default
          call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
        end select
      case (CABLE_OUTPUT_SHAPE_TYPE_BASE_RAD)
        select type(aggregator)
        type is (aggregator_real32_2d_t)
          output_var%temp_buffer_real32_2d => temp_buffer_land_rad_real32
        type is (aggregator_real64_2d_t)
          output_var%temp_buffer_real64_2d => temp_buffer_land_rad_real64
        class default
          call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
        end select
      case (CABLE_OUTPUT_SHAPE_TYPE_BASE_PLANTCARBON)
        select type(aggregator)
        type is (aggregator_real32_2d_t)
          output_var%temp_buffer_real32_2d => temp_buffer_land_plantcarbon_real32
        type is (aggregator_real64_2d_t)
          output_var%temp_buffer_real64_2d => temp_buffer_land_plantcarbon_real64
        class default
          call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
        end select
      case (CABLE_OUTPUT_SHAPE_TYPE_BASE_SOILCARBON)
        select type(aggregator)
        type is (aggregator_real32_2d_t)
          output_var%temp_buffer_real32_2d => temp_buffer_land_soilcarbon_real32
        type is (aggregator_real64_2d_t)
          output_var%temp_buffer_real64_2d => temp_buffer_land_soilcarbon_real64
        class default
          call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
        end select
      case default
        call cable_abort("Unexpected shape_type", __FILE__, __LINE__)
      end select
    end if

    select case(shape_type)
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE)
      select type(aggregator)
      type is (aggregator_int32_1d_t)
        output_var%decomp => output_decomp_base_int32
      type is (aggregator_real32_1d_t)
        output_var%decomp => output_decomp_base_real32
      type is (aggregator_real64_1d_t)
        output_var%decomp => output_decomp_base_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_SOIL)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_soil_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_soil_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_soil_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_SNOW)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_snow_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_snow_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_snow_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_RAD)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_rad_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_rad_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_rad_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_PLANTCARBON)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_plantcarbon_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_plantcarbon_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_plantcarbon_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_SOILCARBON)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_soilcarbon_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_soilcarbon_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_soilcarbon_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH)
      select type(aggregator)
      type is (aggregator_int32_1d_t)
        output_var%decomp => output_decomp_base_patch_int32
      type is (aggregator_real32_1d_t)
        output_var%decomp => output_decomp_base_patch_real32
      type is (aggregator_real64_1d_t)
        output_var%decomp => output_decomp_base_patch_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SOIL)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_patch_soil_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_patch_soil_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_patch_soil_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SNOW)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_patch_snow_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_patch_snow_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_patch_snow_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_RAD)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_patch_rad_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_patch_rad_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_patch_rad_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_PLANTCARBON)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_patch_plantcarbon_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_patch_plantcarbon_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_patch_plantcarbon_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case (CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SOILCARBON)
      select type(aggregator)
      type is (aggregator_int32_2d_t)
        output_var%decomp => output_decomp_base_patch_soilcarbon_int32
      type is (aggregator_real32_2d_t)
        output_var%decomp => output_decomp_base_patch_soilcarbon_real32
      type is (aggregator_real64_2d_t)
        output_var%decomp => output_decomp_base_patch_soilcarbon_real64
      class default
        call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
      end select
    case default
      call cable_abort("Unexpected shape_type", __FILE__, __LINE__)
    end select

    if (.not. allocated(global_profile%output_variables)) then
      global_profile%output_variables = [output_var]
    else
      global_profile%output_variables = [global_profile%output_variables, output_var]
    end if

  end subroutine cable_output_add_variable

  subroutine cable_output_add_aggregator(aggregator, accumulation_frequency, aggregation_frequency, output_aggregator)
    class(aggregator_t), intent(in) :: aggregator
    character(len=*), intent(in) :: accumulation_frequency
    character(len=*), intent(in) :: aggregation_frequency
    type(cable_output_aggregator_t), intent(out) :: output_aggregator

    output_aggregator = cable_output_aggregator_t( &
      accumulation_frequency=accumulation_frequency, &
      aggregation_frequency=aggregation_frequency, &
      aggregator_handle=store_aggregator(aggregator) &
    )

    if (.not. allocated(global_profile%output_aggregators)) then
      global_profile%output_aggregators = [output_aggregator]
    else
      global_profile%output_aggregators = [global_profile%output_aggregators, output_aggregator]
    end if

  end subroutine cable_output_add_aggregator

  subroutine cable_output_commit()
    class(cable_netcdf_file_t), allocatable :: output_file
    integer :: i

    output_file = cable_netcdf_create_file("test_output.nc") ! TODO(Sean): use filename from namelist

    call output_file%def_dims(["x", "y"], [xdimsize, ydimsize])
    call output_file%def_dims(["patch"], [max_vegpatches])
    call output_file%def_dims(["soil"], [ms])
    call output_file%def_dims(["rad"], [nrb])
    call output_file%def_dims(["soil_carbon_pools"], [ncs])
    call output_file%def_dims(["plant_carbon_pools"], [ncp])
    call output_file%def_dims(["time"], [CABLE_NETCDF_UNLIMITED])
    call output_file%def_dims(["nv"], [2])

    if (requires_x_y_output_grid(output_options%grid, metgrid)) then
      call output_file%def_dims(["z"], [1]) ! Atmospheric 'z' dim of size 1 to comply with ALMA grid type
    else if (requires_land_output_grid(output_options%grid, metgrid)) then
      call output_file%def_dims(["land"], [mland])
      call output_file%def_var("local_lat", ["land"], CABLE_NETCDF_FLOAT)
      call output_file%put_att("local_lat", "units", "degrees_north")
      call output_file%def_var("local_lon", ["land"], CABLE_NETCDF_FLOAT)
      call output_file%put_att("local_lon", "units", "degrees_east")
    else
      call cable_abort("Error: Unable to determine output grid type", __FILE__, __LINE__)
    end if

    call output_file%def_var("time", ["time"], CABLE_NETCDF_DOUBLE)
    call output_file%put_att("time", "units", timeunits)
    call output_file%put_att("time", "coordinate", time_coord)
    call output_file%put_att("time", "calendar", calendar)
    call output_file%put_att("time", "bounds", "time_bnds")
    call output_file%def_var("time_bnds", ["nv", "time"], CABLE_NETCDF_DOUBLE)

    ! Define latitude and longitude variable (ALMA):
    call output_file%def_var("latitude", ["x", "y"], CABLE_NETCDF_FLOAT)
    call output_file%put_att("latitude", "units", "degrees_north")
    call output_file%def_var("longitude", ["x", "y"], CABLE_NETCDF_FLOAT)
    call output_file%put_att("longitude", "units", "degrees_east")

    ! Write "cordinate variables" to enable reading by GrADS:
    call output_file%def_var("x", ["x"], CABLE_NETCDF_FLOAT)
    call output_file%put_att("x", "units", "degrees_east")
    call output_file%put_att("x", "comment", "x coordinate variable for GrADS compatibility")
    call output_file%def_var("y", ["y"], CABLE_NETCDF_FLOAT)
    call output_file%put_att("y", "units", "degrees_north")
    call output_file%put_att("y", "comment", "y coordinate variable for GrADS compatibility")

    ! TODO(Sean): define remaining coordinate variables

    ! TODO(Sean): add global attributes

    global_profile%output_variables = pack(global_profile%output_variables, global_profile%output_variables(:)%active)

    do i = 1, size(global_profile%output_variables)
      associate(output_var => global_profile%output_variables(i))
        call output_file%def_var( &
          var_name=output_var%name, &
          dim_names=output_var%dims, &
          type=output_var%var_type &
        )
        call output_file%put_att(output_var%name, 'units', output_var%units)
        call output_file%put_att(output_var%name, 'long_name', output_var%long_name)
        select case (output_var%var_type)
        case (CABLE_NETCDF_INT)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_INT32)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_INT32)
        case (CABLE_NETCDF_FLOAT)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_REAL32)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_REAL32)
        case (CABLE_NETCDF_DOUBLE)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_REAL64)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_REAL64)
        end select
        ! TODO(Sean): set cell_methods attribute
      end associate
    end do

    global_profile%output_file = output_file

    call sort_aggregators_by_frequency(global_profile%output_aggregators)

    ! Initialize all aggregators
    do i = 1, size(global_profile%output_aggregators)
      associate(aggregator => global_profile%output_aggregators(i)%aggregator_handle%aggregator)
        call aggregator%init()
      end associate
    end do

  end subroutine

  subroutine cable_output_update(time_index, dels, leaps, start_year, patch, landpt)
    integer, intent(in) :: time_index
    real, intent(in) :: dels
    logical, intent(in) :: leaps
    integer, intent(in) :: start_year
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)

    real :: current_time
    integer :: i

    do i = 1, size(global_profile%output_aggregators)
      associate(output_aggregator => global_profile%output_aggregators(i))
        if (time_step_matches(dels, time_index, output_aggregator%accumulation_frequency, leaps, start_year)) then
          call output_aggregator%aggregator_handle%accumulate()
        end if
        if (time_step_matches(dels, time_index, output_aggregator%aggregation_frequency, leaps, start_year)) then
          call output_aggregator%aggregator_handle%normalise()
        end if
      end associate
    end do

    if (time_step_matches(dels, time_index, output_options%averaging, leaps, start_year)) then

      do i = 1, size(global_profile%output_variables)
        associate(output_variable => global_profile%output_variables(i))
          if (output_variable%grid_cell_averaging) then
            call write_variable_grid_cell_average(output_variable, global_profile%output_file, global_profile%frame + 1, patch, landpt)
          else
            call write_variable(output_variable, global_profile%output_file, global_profile%frame + 1)
          end if
        end associate
      end do

      current_time = time_index * dels
      call global_profile%output_file%put_var("time", (current_time + global_profile%previous_write_time) / 2.0, start=[global_profile%frame + 1])
      call global_profile%output_file%put_var("time_bnds", [global_profile%previous_write_time, current_time], start=[1, global_profile%frame + 1])
      global_profile%previous_write_time = current_time
      global_profile%frame = global_profile%frame + 1

    end if

    do i = 1, size(global_profile%output_aggregators)
      associate(output_aggregator => global_profile%output_aggregators(i))
        if (time_step_matches(dels, time_index, output_aggregator%aggregation_frequency, leaps, start_year)) then
          call output_aggregator%aggregator_handle%reset()
        end if
      end associate
    end do

  end subroutine cable_output_update

  subroutine write_variable(output_variable, output_file, time_index)
    type(cable_output_variable_t), intent(inout) :: output_variable
    class(cable_netcdf_file_t), intent(inout) :: output_file
    integer, intent(in) :: time_index

    select type (aggregator => output_variable%output_aggregator%aggregator_handle%aggregator)
    type is (aggregator_int32_1d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            frame=time_index)
    type is (aggregator_int32_2d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            frame=time_index)
    type is (aggregator_int32_3d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            frame=time_index)
    type is (aggregator_real32_1d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real32_2d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real32_3d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real64_1d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    type is (aggregator_real64_2d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    type is (aggregator_real64_3d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%storage, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine write_variable

  subroutine write_variable_grid_cell_average(output_variable, output_file, time_index, patch, landpt)
    type(cable_output_variable_t), intent(inout) :: output_variable
    class(cable_netcdf_file_t), intent(inout) :: output_file
    integer, intent(in) :: time_index
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)

    select type (aggregator => output_variable%output_aggregator%aggregator_handle%aggregator)
    type is (aggregator_real32_1d_t)
      call grid_cell_average( &
            input_array=aggregator%storage, &
            output_array=output_variable%temp_buffer_real32_1d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real32_1d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real32_2d_t)
      call grid_cell_average( &
            input_array=aggregator%storage, &
            output_array=output_variable%temp_buffer_real32_2d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real32_2d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real32_3d_t)
      call grid_cell_average( &
            input_array=aggregator%storage, &
            output_array=output_variable%temp_buffer_real32_3d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real32_3d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real64_1d_t)
      call grid_cell_average( &
            input_array=aggregator%storage, &
            output_array=output_variable%temp_buffer_real64_1d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real64_1d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    type is (aggregator_real64_2d_t)
      call grid_cell_average( &
            input_array=aggregator%storage, &
            output_array=output_variable%temp_buffer_real64_2d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real64_2d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    type is (aggregator_real64_3d_t)
      call grid_cell_average( &
            input_array=aggregator%storage, &
            output_array=output_variable%temp_buffer_real64_3d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real64_3d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine write_variable_grid_cell_average

end module
