module aggregator_types_mod
  use iso_fortran_env, only: int32, real32, real64
  use cable_abort_module, only: cable_abort
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

  type, abstract :: aggregator_t
    integer :: counter = 0
    procedure(accumulate_data), pointer :: accumulate
    procedure(normalise_data), pointer :: normalise
    procedure(reset_data), pointer :: reset
  contains
    procedure :: init => aggregator_init
    procedure :: set_method => aggregator_set_method
  end type aggregator_t

  abstract interface
    subroutine accumulate_data(this)
      import aggregator_t
      class(aggregator_t), intent(inout) :: this
    end subroutine accumulate_data
    subroutine normalise_data(this)
      import aggregator_t
      class(aggregator_t), intent(inout) :: this
    end subroutine normalise_data
    subroutine reset_data(this)
      import aggregator_t
      class(aggregator_t), intent(inout) :: this
    end subroutine reset_data
  end interface

  type aggregator_handle_t
    class(aggregator_t), pointer :: aggregator => null()
  contains
    procedure :: init => aggregator_handle_init
    procedure :: accumulate => aggregator_handle_accumulate
    procedure :: normalise => aggregator_handle_normalise
    procedure :: reset => aggregator_handle_reset
  end type aggregator_handle_t

  type, extends(aggregator_t) :: aggregator_int32_1d_t
    integer(kind=int32), dimension(:), allocatable :: storage
    integer(kind=int32), dimension(:), pointer :: source_data => null()
  end type aggregator_int32_1d_t

  type, extends(aggregator_t) :: aggregator_int32_2d_t
    integer(kind=int32), dimension(:,:), allocatable :: storage
    integer(kind=int32), dimension(:,:), pointer :: source_data => null()
  end type aggregator_int32_2d_t

  type, extends(aggregator_t) :: aggregator_int32_3d_t
    integer(kind=int32), dimension(:,:,:), allocatable :: storage
    integer(kind=int32), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_int32_3d_t

  type, extends(aggregator_t) :: aggregator_real32_1d_t
    real(kind=real32), dimension(:), allocatable :: storage
    real(kind=real32), dimension(:), pointer :: source_data => null()
  end type aggregator_real32_1d_t

  type, extends(aggregator_t) :: aggregator_real32_2d_t
    real(kind=real32), dimension(:,:), allocatable :: storage
    real(kind=real32), dimension(:,:), pointer :: source_data => null()
  end type aggregator_real32_2d_t

  type, extends(aggregator_t) :: aggregator_real32_3d_t
    real(kind=real32), dimension(:,:,:), allocatable :: storage
    real(kind=real32), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_real32_3d_t

  type, extends(aggregator_t) :: aggregator_real64_1d_t
    real(kind=real64), dimension(:), allocatable :: storage
    real(kind=real64), dimension(:), pointer :: source_data => null()
  end type aggregator_real64_1d_t

  type, extends(aggregator_t) :: aggregator_real64_2d_t
    real(kind=real64), dimension(:,:), allocatable :: storage
    real(kind=real64), dimension(:,:), pointer :: source_data => null()
  end type aggregator_real64_2d_t

  type, extends(aggregator_t) :: aggregator_real64_3d_t
    real(kind=real64), dimension(:,:,:), allocatable :: storage
    real(kind=real64), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_real64_3d_t

contains

  subroutine aggregator_handle_init(this)
    class(aggregator_handle_t), intent(inout) :: this

    call this%aggregator%init()

  end subroutine aggregator_handle_init

  subroutine aggregator_handle_accumulate(this)
    class(aggregator_handle_t), intent(inout) :: this

    call this%aggregator%accumulate()

  end subroutine aggregator_handle_accumulate

  subroutine aggregator_handle_normalise(this)
    class(aggregator_handle_t), intent(inout) :: this

    call this%aggregator%normalise()

  end subroutine aggregator_handle_normalise

  subroutine aggregator_handle_reset(this)
    class(aggregator_handle_t), intent(inout) :: this

    call this%aggregator%reset()

  end subroutine aggregator_handle_reset

  subroutine aggregator_init(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    type is (aggregator_int32_2d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    type is (aggregator_int32_3d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    type is (aggregator_real32_1d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    type is (aggregator_real32_2d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    type is (aggregator_real32_3d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    type is (aggregator_real64_1d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    type is (aggregator_real64_2d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    type is (aggregator_real64_3d_t)
      if (.not. allocated(this%storage)) allocate(this%storage, mold=this%source_data)
    end select

    call this%reset()

  end subroutine aggregator_init

  subroutine aggregator_set_method(this, method)
    class(aggregator_t), intent(inout) :: this
    character(len=*), intent(in) :: method

    if (method == "mean") then
      this%accumulate => sum_accumulate
      this%normalise => mean_normalise
      this%reset => other_reset
    elseif (method == "sum") then
      this%accumulate => sum_accumulate
      this%normalise => other_normalise
      this%reset => other_reset
    elseif (method == "point") then
      this%accumulate => point_accumulate
      this%normalise => other_normalise
      this%reset => point_reset
    elseif (method == "min") then
      this%accumulate => min_accumulate
      this%normalise => other_normalise
      this%reset => min_reset
    elseif (method == "max") then
      this%accumulate => max_accumulate
      this%normalise => other_normalise
      this%reset => max_reset
    else
      call cable_abort("Aggregation method "//method//" is invalid.")
    endif

  end subroutine aggregator_set_method

  subroutine sum_accumulate(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      this%storage = this%storage + this%source_data
    type is (aggregator_int32_2d_t)
      this%storage = this%storage + this%source_data
    type is (aggregator_int32_3d_t)
      this%storage = this%storage + this%source_data
    type is (aggregator_real32_1d_t)
      this%storage = this%storage + this%source_data
    type is (aggregator_real32_2d_t)
      this%storage = this%storage + this%source_data
    type is (aggregator_real32_3d_t)
      this%storage = this%storage + this%source_data
    type is (aggregator_real64_1d_t)
      this%storage = this%storage + this%source_data
    type is (aggregator_real64_2d_t)
      this%storage = this%storage + this%source_data
    type is (aggregator_real64_3d_t)
      this%storage = this%storage + this%source_data
    end select

    this%counter = this%counter + 1

  end subroutine sum_accumulate

  subroutine point_accumulate(this)
    class(aggregator_t), intent(inout) :: this
  end subroutine point_accumulate

  subroutine min_accumulate(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      this%storage = min(this%storage, this%source_data)
    type is (aggregator_int32_2d_t)
      this%storage = min(this%storage, this%source_data)
    type is (aggregator_int32_3d_t)
      this%storage = min(this%storage, this%source_data)
    type is (aggregator_real32_1d_t)
      this%storage = min(this%storage, this%source_data)
    type is (aggregator_real32_2d_t)
      this%storage = min(this%storage, this%source_data)
    type is (aggregator_real32_3d_t)
      this%storage = min(this%storage, this%source_data)
    type is (aggregator_real64_1d_t)
      this%storage = min(this%storage, this%source_data)
    type is (aggregator_real64_2d_t)
      this%storage = min(this%storage, this%source_data)
    type is (aggregator_real64_3d_t)
      this%storage = min(this%storage, this%source_data)
    end select

    this%counter = this%counter + 1

  end subroutine min_accumulate

  subroutine max_accumulate(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      this%storage = max(this%storage, this%source_data)
    type is (aggregator_int32_2d_t)
      this%storage = max(this%storage, this%source_data)
    type is (aggregator_int32_3d_t)
      this%storage = max(this%storage, this%source_data)
    type is (aggregator_real32_1d_t)
      this%storage = max(this%storage, this%source_data)
    type is (aggregator_real32_2d_t)
      this%storage = max(this%storage, this%source_data)
    type is (aggregator_real32_3d_t)
      this%storage = max(this%storage, this%source_data)
    type is (aggregator_real64_1d_t)
      this%storage = max(this%storage, this%source_data)
    type is (aggregator_real64_2d_t)
      this%storage = max(this%storage, this%source_data)
    type is (aggregator_real64_3d_t)
      this%storage = max(this%storage, this%source_data)
    end select

    this%counter = this%counter + 1

  end subroutine max_accumulate

  subroutine mean_normalise(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      this%storage = this%storage / this%counter
    type is (aggregator_int32_2d_t)
      this%storage = this%storage / this%counter
    type is (aggregator_int32_3d_t)
      this%storage = this%storage / this%counter
    type is (aggregator_real32_1d_t)
      this%storage = this%storage / this%counter
    type is (aggregator_real32_2d_t)
      this%storage = this%storage / this%counter
    type is (aggregator_real32_3d_t)
      this%storage = this%storage / this%counter
    type is (aggregator_real64_1d_t)
      this%storage = this%storage / this%counter
    type is (aggregator_real64_2d_t)
      this%storage = this%storage / this%counter
    type is (aggregator_real64_3d_t)
      this%storage = this%storage / this%counter
    end select

  end subroutine mean_normalise

  subroutine point_normalise(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      this%storage = this%source_data
    type is (aggregator_int32_2d_t)
      this%storage = this%source_data
    type is (aggregator_int32_3d_t)
      this%storage = this%source_data
    type is (aggregator_real32_1d_t)
      this%storage = this%source_data
    type is (aggregator_real32_2d_t)
      this%storage = this%source_data
    type is (aggregator_real32_3d_t)
      this%storage = this%source_data
    type is (aggregator_real64_1d_t)
      this%storage = this%source_data
    type is (aggregator_real64_2d_t)
      this%storage = this%source_data
    type is (aggregator_real64_3d_t)
      this%storage = this%source_data
    end select

  end subroutine point_normalise

  subroutine other_normalise(this)
    class(aggregator_t), intent(inout) :: this
  end subroutine other_normalise

  subroutine point_reset(this)
    class(aggregator_t), intent(inout) :: this
  end subroutine point_reset

  subroutine min_reset(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      this%storage = huge(int(0_int32))
    type is (aggregator_int32_2d_t)
      this%storage = huge(int(0_int32))
    type is (aggregator_int32_3d_t)
      this%storage = huge(int(0_int32))
    type is (aggregator_real32_1d_t)
      this%storage = huge(real(0.0_real32))
    type is (aggregator_real32_2d_t)
      this%storage = huge(real(0.0_real32))
    type is (aggregator_real32_3d_t)
      this%storage = huge(real(0.0_real32))
    type is (aggregator_real64_1d_t)
      this%storage = huge(real(0.0_real64))
    type is (aggregator_real64_2d_t)
      this%storage = huge(real(0.0_real64))
    type is (aggregator_real64_3d_t)
      this%storage = huge(real(0.0_real64))
    end select

    this%counter = 0

  end subroutine min_reset

  subroutine max_reset(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      this%storage = -huge(int(0_int32))
    type is (aggregator_int32_2d_t)
      this%storage = -huge(int(0_int32))
    type is (aggregator_int32_3d_t)
      this%storage = -huge(int(0_int32))
    type is (aggregator_real32_1d_t)
      this%storage = -huge(real(0.0_real32))
    type is (aggregator_real32_2d_t)
      this%storage = -huge(real(0.0_real32))
    type is (aggregator_real32_3d_t)
      this%storage = -huge(real(0.0_real32))
    type is (aggregator_real64_1d_t)
      this%storage = -huge(real(0.0_real64))
    type is (aggregator_real64_2d_t)
      this%storage = -huge(real(0.0_real64))
    type is (aggregator_real64_3d_t)
      this%storage = -huge(real(0.0_real64))
    end select

    this%counter = 0

  end subroutine max_reset

  subroutine other_reset(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_1d_t)
      this%storage = 0_int32
    type is (aggregator_int32_2d_t)
      this%storage = 0_int32
    type is (aggregator_int32_3d_t)
      this%storage = 0_int32
    type is (aggregator_real32_1d_t)
      this%storage = 0.0_real32
    type is (aggregator_real32_2d_t)
      this%storage = 0.0_real32
    type is (aggregator_real32_3d_t)
      this%storage = 0.0_real32
    type is (aggregator_real64_1d_t)
      this%storage = 0.0_real64
    type is (aggregator_real64_2d_t)
      this%storage = 0.0_real64
    type is (aggregator_real64_3d_t)
      this%storage = 0.0_real64
    end select

    this%counter = 0

  end subroutine other_reset

end module
