! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

submodule (cable_output_mod:cable_output_common_smod) cable_output_reduction_buffers_smod
  !* Implementation of procedures used for managing the grid reduction buffers
  ! used in the output system.

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

  implicit none

  integer(kind=int32), allocatable, target :: temp_buffer_land_int32(:)
    !! Grid reduction buffer for data with shape `[mp]` and type `integer(int32)`.
  real(kind=real32),   allocatable, target :: temp_buffer_land_real32(:)
    !! Grid reduction buffer for data with shape `[mp]` and type `real(real32)`.
  real(kind=real64),   allocatable, target :: temp_buffer_land_real64(:)
    !! Grid reduction buffer for data with shape `[mp]` and type `real(real64)`.
  integer(kind=int32), allocatable, target :: temp_buffer_land_soil_int32(:, :)
    !! Grid reduction buffer for data with shape `[mp, ms]` and type `integer(int32)`.
  real(kind=real32),   allocatable, target :: temp_buffer_land_soil_real32(:, :)
    !! Grid reduction buffer for data with shape `[mp, ms]` and type `real(real32)`.
  real(kind=real64),   allocatable, target :: temp_buffer_land_soil_real64(:, :)
    !! Grid reduction buffer for data with shape `[mp, ms]` and type `real(real64)`.
  integer(kind=int32), allocatable, target :: temp_buffer_land_snow_int32(:, :)
    !! Grid reduction buffer for data with shape `[mp, msn]` and type `integer(int32)`.
  real(kind=real32),   allocatable, target :: temp_buffer_land_snow_real32(:, :)
    !! Grid reduction buffer for data with shape `[mp, msn]` and type `real(real32)`.
  real(kind=real64),   allocatable, target :: temp_buffer_land_snow_real64(:, :)
    !! Grid reduction buffer for data with shape `[mp, msn]` and type `real(real64)`.
  integer(kind=int32), allocatable, target :: temp_buffer_land_rad_int32(:, :)
    !! Grid reduction buffer for data with shape `[mp, nrb]` and type `integer(int32)`.
  real(kind=real32),   allocatable, target :: temp_buffer_land_rad_real32(:, :)
    !! Grid reduction buffer for data with shape `[mp, nrb]` and type `real(real32)`.
  real(kind=real64),   allocatable, target :: temp_buffer_land_rad_real64(:, :)
    !! Grid reduction buffer for data with shape `[mp, nrb]` and type `real(real64)`.
  integer(kind=int32), allocatable, target :: temp_buffer_land_plantcarbon_int32(:, :)
    !! Grid reduction buffer for data with shape `[mp, ncp]` and type `integer(int32)`.
  real(kind=real32),   allocatable, target :: temp_buffer_land_plantcarbon_real32(:, :)
    !! Grid reduction buffer for data with shape `[mp, ncp]` and type `real(real32)`.
  real(kind=real64),   allocatable, target :: temp_buffer_land_plantcarbon_real64(:, :)
    !! Grid reduction buffer for data with shape `[mp, ncp]` and type `real(real64)`.
  integer(kind=int32), allocatable, target :: temp_buffer_land_soilcarbon_int32(:, :)
    !! Grid reduction buffer for data with shape `[mp, ncs]` and type `integer(int32)`.
  real(kind=real32),   allocatable, target :: temp_buffer_land_soilcarbon_real32(:, :)
    !! Grid reduction buffer for data with shape `[mp, ncs]` and type `real(real32)`.
  real(kind=real64),   allocatable, target :: temp_buffer_land_soilcarbon_real64(:, :)
    !! Grid reduction buffer for data with shape `[mp, ncs]` and type `real(real64)`.

contains

  module subroutine cable_output_reduction_buffers_init()
    !! Initialises the buffers used for performing grid reductions in the output system.

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

  module subroutine cable_output_reduction_buffers_free()
    !! Deallocates the buffers used for performing grid reductions in the output system.

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

  module subroutine cable_output_reduction_buffers_associate_1d_int32(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 1D 32-bit integer variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    integer(kind=int32), pointer, intent(inout) :: temp_buffer(:)
      !! The pointer array to associate with the appropriate reduction buffer.

    if (array_eq(output_var%data_shape(:)%size(), [mp])) then
      temp_buffer => temp_buffer_land_int32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine cable_output_reduction_buffers_associate_1d_int32

  module subroutine cable_output_reduction_buffers_associate_1d_real32(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 1D 32-bit real variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    real(kind=real32), pointer, intent(inout) :: temp_buffer(:)
      !! The pointer array to associate with the appropriate reduction buffer.

    if (array_eq(output_var%data_shape(:)%size(), [mp])) then
      temp_buffer => temp_buffer_land_real32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine cable_output_reduction_buffers_associate_1d_real32

  module subroutine cable_output_reduction_buffers_associate_1d_real64(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 1D 64-bit real variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    real(kind=real64), pointer, intent(inout) :: temp_buffer(:)
      !! The pointer array to associate with the appropriate reduction buffer.

    if (array_eq(output_var%data_shape(:)%size(), [mp])) then
      temp_buffer => temp_buffer_land_real64
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine cable_output_reduction_buffers_associate_1d_real64

  module subroutine cable_output_reduction_buffers_associate_2d_int32(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 2D 32-bit integer variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    integer(kind=int32), pointer, intent(inout) :: temp_buffer(:,:)
      !! The pointer array to associate with the appropriate reduction buffer.

    if (array_eq(output_var%data_shape(:)%size(), [mp, ms])) then
      temp_buffer => temp_buffer_land_soil_int32
    else if (array_eq(output_var%data_shape(:)%size(), [mp, msn])) then
      temp_buffer => temp_buffer_land_snow_int32
    else if (array_eq(output_var%data_shape(:)%size(), [mp, nrb])) then
      temp_buffer => temp_buffer_land_rad_int32
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncp])) then
      temp_buffer => temp_buffer_land_plantcarbon_int32
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncs])) then
      temp_buffer => temp_buffer_land_soilcarbon_int32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine cable_output_reduction_buffers_associate_2d_int32

  module subroutine cable_output_reduction_buffers_associate_2d_real32(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 2D 32-bit real variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    real(kind=real32), pointer, intent(inout) :: temp_buffer(:,:)
      !! The pointer array to associate with the appropriate reduction buffer.

    if (array_eq(output_var%data_shape(:)%size(), [mp, ms])) then
      temp_buffer => temp_buffer_land_soil_real32
    else if (array_eq(output_var%data_shape(:)%size(), [mp, msn])) then
      temp_buffer => temp_buffer_land_snow_real32
    else if (array_eq(output_var%data_shape(:)%size(), [mp, nrb])) then
      temp_buffer => temp_buffer_land_rad_real32
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncp])) then
      temp_buffer => temp_buffer_land_plantcarbon_real32
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncs])) then
      temp_buffer => temp_buffer_land_soilcarbon_real32
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine cable_output_reduction_buffers_associate_2d_real32

  module subroutine cable_output_reduction_buffers_associate_2d_real64(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 2D 64-bit real variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    real(kind=real64), pointer, intent(inout) :: temp_buffer(:,:)
      !! The pointer array to associate with the appropriate reduction buffer.

    if (array_eq(output_var%data_shape(:)%size(), [mp, ms])) then
      temp_buffer => temp_buffer_land_soil_real64
    else if (array_eq(output_var%data_shape(:)%size(), [mp, msn])) then
      temp_buffer => temp_buffer_land_snow_real64
    else if (array_eq(output_var%data_shape(:)%size(), [mp, nrb])) then
      temp_buffer => temp_buffer_land_rad_real64
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncp])) then
      temp_buffer => temp_buffer_land_plantcarbon_real64
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncs])) then
      temp_buffer => temp_buffer_land_soilcarbon_real64
    else
      call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
    end if

  end subroutine cable_output_reduction_buffers_associate_2d_real64

  module subroutine cable_output_reduction_buffers_associate_3d_int32(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 3D 32-bit integer variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    integer(kind=int32), pointer, intent(inout) :: temp_buffer(:,:,:)
      !! The pointer array to associate with the appropriate reduction buffer.

    call cable_abort("Grid reduction buffers not implemented for this data type.", __FILE__, __LINE__)

  end subroutine cable_output_reduction_buffers_associate_3d_int32

  module subroutine cable_output_reduction_buffers_associate_3d_real32(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 3D 32-bit real variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    real(kind=real32), pointer, intent(inout) :: temp_buffer(:,:,:)
      !! The pointer array to associate with the appropriate reduction buffer.

    call cable_abort("Grid reduction buffers not implemented for this data type.", __FILE__, __LINE__)

  end subroutine cable_output_reduction_buffers_associate_3d_real32

  module subroutine cable_output_reduction_buffers_associate_3d_real64(output_var, temp_buffer)
    !! The reduction buffer association subroutine for 3D 64-bit real variables.
    type(cable_output_variable_t), intent(inout) :: output_var
      !! The output variable for which to associate the reduction buffer.
    real(kind=real64), pointer, intent(inout) :: temp_buffer(:,:,:)
      !! The pointer array to associate with the appropriate reduction buffer.

    call cable_abort("Grid reduction buffers not implemented for this data type.", __FILE__, __LINE__)

  end subroutine cable_output_reduction_buffers_associate_3d_real64

end submodule cable_output_reduction_buffers_smod
