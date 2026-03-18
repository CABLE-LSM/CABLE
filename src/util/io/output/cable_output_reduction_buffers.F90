module cable_output_reduction_buffers_mod

  use cable_error_handler_mod, only: cable_abort

  use cable_array_utils_mod, only: array_eq

  use iso_fortran_env, only: int32, real32, real64

  use cable_def_types_mod, only: mland
  use cable_def_types_mod, only: mp
  use cable_def_types_mod, only: ms
  use cable_def_types_mod, only: msn
  use cable_def_types_mod, only: nrb
  use cable_def_types_mod, only: ncs
  use cable_def_types_mod, only: ncp

  use cable_output_types_mod, only: cable_output_variable_t

  implicit none
  private

  public :: allocate_grid_reduction_buffers
  public :: deallocate_grid_reduction_buffers
  public :: associate_temp_buffer

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

  interface associate_temp_buffer
    module procedure associate_temp_buffer_1d_int32
    module procedure associate_temp_buffer_1d_real32
    module procedure associate_temp_buffer_1d_real64
    module procedure associate_temp_buffer_2d_int32
    module procedure associate_temp_buffer_2d_real32
    module procedure associate_temp_buffer_2d_real64
    module procedure associate_temp_buffer_3d_int32
    module procedure associate_temp_buffer_3d_real32
    module procedure associate_temp_buffer_3d_real64
  end interface

contains

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

  subroutine associate_temp_buffer_1d_int32(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    integer(kind=int32), pointer, intent(inout) :: temp_buffer(:)

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      temp_buffer => temp_buffer_land_int32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_1d_int32

  subroutine associate_temp_buffer_1d_real32(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    real(kind=real32), pointer, intent(inout) :: temp_buffer(:)

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      temp_buffer => temp_buffer_land_real32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_1d_real32

  subroutine associate_temp_buffer_1d_real64(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    real(kind=real64), pointer, intent(inout) :: temp_buffer(:)

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      temp_buffer => temp_buffer_land_real64
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_1d_real64

  subroutine associate_temp_buffer_2d_int32(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    integer(kind=int32), pointer, intent(inout) :: temp_buffer(:,:)

    if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      temp_buffer => temp_buffer_land_soil_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      temp_buffer => temp_buffer_land_snow_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      temp_buffer => temp_buffer_land_rad_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      temp_buffer => temp_buffer_land_plantcarbon_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      temp_buffer => temp_buffer_land_soilcarbon_int32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_2d_int32

  subroutine associate_temp_buffer_2d_real32(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    real(kind=real32), pointer, intent(inout) :: temp_buffer(:,:)

    if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      temp_buffer => temp_buffer_land_soil_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      temp_buffer => temp_buffer_land_snow_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      temp_buffer => temp_buffer_land_rad_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      temp_buffer => temp_buffer_land_plantcarbon_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      temp_buffer => temp_buffer_land_soilcarbon_real32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_2d_real32

  subroutine associate_temp_buffer_2d_real64(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    real(kind=real64), pointer, intent(inout) :: temp_buffer(:,:)

    if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      temp_buffer => temp_buffer_land_soil_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      temp_buffer => temp_buffer_land_snow_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      temp_buffer => temp_buffer_land_rad_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      temp_buffer => temp_buffer_land_plantcarbon_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      temp_buffer => temp_buffer_land_soilcarbon_real64
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine associate_temp_buffer_2d_real64

  subroutine associate_temp_buffer_3d_int32(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    integer(kind=int32), pointer, intent(inout) :: temp_buffer(:,:,:)

    call cable_abort("Grid reduction buffers not implemented for this data type.", __FILE__, __LINE__)

  end subroutine associate_temp_buffer_3d_int32

  subroutine associate_temp_buffer_3d_real32(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    real(kind=real32), pointer, intent(inout) :: temp_buffer(:,:,:)

    call cable_abort("Grid reduction buffers not implemented for this data type.", __FILE__, __LINE__)

  end subroutine associate_temp_buffer_3d_real32

  subroutine associate_temp_buffer_3d_real64(output_var, temp_buffer)
    type(cable_output_variable_t), intent(inout) :: output_var
    real(kind=real64), pointer, intent(inout) :: temp_buffer(:,:,:)

    call cable_abort("Grid reduction buffers not implemented for this data type.", __FILE__, __LINE__)

  end subroutine associate_temp_buffer_3d_real64

end module
