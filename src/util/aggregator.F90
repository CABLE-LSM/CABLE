module aggregator_mod

  use iso_fortran_env, only: int32, real32, real64
  use cable_abort_module, only: cable_abort

  use aggregator_types_mod

  implicit none
  private

  public :: aggregator_t
  public :: aggregator_handle_t
  public :: aggregator_int32_1d_t
  public :: aggregator_int32_2d_t
  public :: aggregator_int32_3d_t
  public :: aggregator_real32_1d_t
  public :: aggregator_real32_2d_t
  public :: aggregator_real32_3d_t
  public :: aggregator_real64_1d_t
  public :: aggregator_real64_2d_t
  public :: aggregator_real64_3d_t
  public :: aggregator_mod_init
  public :: aggregator_mod_end
  public :: new_aggregator
  public :: store_aggregator

  type aggregator_store_t
    class(aggregator_t), allocatable :: aggregator
  end type aggregator_store_t

  interface new_aggregator
    module procedure new_aggregator_int32_1d_t
    module procedure new_aggregator_int32_2d_t
    module procedure new_aggregator_int32_3d_t
    module procedure new_aggregator_real32_1d
    module procedure new_aggregator_real32_2d
    module procedure new_aggregator_real32_3d
    module procedure new_aggregator_real64_1d
    module procedure new_aggregator_real64_2d
    module procedure new_aggregator_real64_3d
  end interface

  integer, parameter :: DEFAULT_MAX_AGGREGATORS = 1000
  integer :: num_aggregators = 0
  type(aggregator_store_t), allocatable, target :: aggregators(:)


contains

  subroutine aggregator_mod_init()
    integer :: max_aggregators
    max_aggregators = DEFAULT_MAX_AGGREGATORS ! TODO(Sean): Make this configurable
    allocate(aggregators(max_aggregators))
    num_aggregators = 0
  end subroutine

  subroutine aggregator_mod_end()
    if (allocated(aggregators)) deallocate(aggregators)
  end subroutine

  function store_aggregator(aggregator) result(aggregator_handle)
    class(aggregator_t), intent(in) :: aggregator
    type(aggregator_handle_t) :: aggregator_handle
    integer :: index

    num_aggregators = num_aggregators + 1
    if (num_aggregators > size(aggregators)) then
      ! Note: we cannot resize the aggregators array as this would
      ! invalidate any pointers to its elements elsewhere.
      ! TODO(Sean): provide a recommendation to increase the max number of aggregators
      call cable_abort("Exceeded maximum number of aggregators.", __FILE__, __LINE__)
    end if

    index = num_aggregators

    ! Copy the aggregator into the store
    aggregators(index)%aggregator = aggregator

    ! Create a handle pointing to the stored aggregator
    aggregator_handle%aggregator => aggregators(index)%aggregator

  end function store_aggregator

  function new_aggregator_int32_1d_t(source_data, method) result(agg)
    integer(kind=int32), dimension(:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_int32_1d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_int32_1d_t

  function new_aggregator_int32_2d_t(source_data, method) result(agg)
    integer(kind=int32), dimension(:,:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_int32_2d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_int32_2d_t

  function new_aggregator_int32_3d_t(source_data, method) result(agg)
    integer(kind=int32), dimension(:,:,:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_int32_3d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_int32_3d_t

  function new_aggregator_real32_1d(source_data, method) result(agg)
    real(kind=real32), dimension(:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_real32_1d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_real32_1d

  function new_aggregator_real32_2d(source_data, method) result(agg)
    real(kind=real32), dimension(:,:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_real32_2d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_real32_2d

  function new_aggregator_real32_3d(source_data, method) result(agg)
    real(kind=real32), dimension(:,:,:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_real32_3d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_real32_3d

  function new_aggregator_real64_1d(source_data, method) result(agg)
    real(kind=real64), dimension(:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_real64_1d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_real64_1d

  function new_aggregator_real64_2d(source_data, method) result(agg)
    real(kind=real64), dimension(:,:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_real64_2d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_real64_2d

  function new_aggregator_real64_3d(source_data, method) result(agg)
    real(kind=real64), dimension(:,:,:), intent(inout), target :: source_data
    character(len=*), intent(in) :: method
    type(aggregator_real64_3d_t) :: agg

    agg%source_data => source_data
    call agg%set_method(method)

  end function new_aggregator_real64_3d

end module
