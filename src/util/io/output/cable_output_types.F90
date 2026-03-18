module cable_output_types_mod

  use iso_fortran_env, only: int32, real32, real64

  use aggregator_mod, only: aggregator_t

  use cable_netcdf_mod, only: cable_netcdf_file_t

  implicit none
  private

  integer, parameter, public :: CABLE_OUTPUT_VAR_TYPE_UNDEFINED = -1

  integer(kind=int32), parameter, public :: FILL_VALUE_INT32  = -9999999_int32
  real(kind=real32),   parameter, public :: FILL_VALUE_REAL32 = -1.0e+33_real32
  real(kind=real64),   parameter, public :: FILL_VALUE_REAL64 = -1.0e+33_real64

  character(32), parameter, public :: allowed_reduction_methods(3) = [ &
    "none                    ", &
    "grid_cell_average       ", &
    "first_patch_in_grid_cell" &
  ]

  character(32), parameter, public :: allowed_aggregation_methods(5) = [ &
      "point", &
      "mean ", &
      "max  ", &
      "min  ", &
      "sum  " &
  ]

  type, public :: cable_output_dim_t
    character(len=64) :: name
    integer :: size
  end type

  type, public :: cable_output_attribute_t
    character(len=64) :: name
    character(len=256) :: value
  end type

  type, public :: cable_output_variable_t
    character(len=64) :: field_name
    character(len=64) :: netcdf_name = ""
    character(len=64) :: accumulation_frequency = "all"
    character(len=64) :: reduction_method = "none"
    character(len=64) :: aggregation_method = "point"
    logical :: active = .true.
    logical :: parameter = .false.
    logical :: distributed = .true.
    logical :: restart = .false.
    logical :: patchout = .false.
    integer :: var_type = CABLE_OUTPUT_VAR_TYPE_UNDEFINED
    real, dimension(2) :: range = [-huge(0.0), huge(0.0)]
    real :: scale_by = 1.0
    real :: divide_by = 1.0
    real :: offset_by = 0.0
    type(cable_output_dim_t), allocatable :: data_shape(:)
    class(aggregator_t), allocatable :: aggregator
    type(cable_output_attribute_t), allocatable :: metadata(:)
  end type

  type, public :: cable_output_profile_t
    real :: previous_write_time = 0.0
    integer :: frame = 0
    character(len=64) :: sampling_frequency
    character(len=64) :: grid_type
    character(len=256) :: file_name
    class(cable_netcdf_file_t), allocatable :: output_file
    type(cable_output_variable_t), allocatable :: coordinate_variables(:)
    type(cable_output_variable_t), allocatable :: output_variables(:)
    type(cable_output_attribute_t), allocatable :: metadata(:)
  end type

end module
