module cable_grid_reductions_mod

  use iso_fortran_env, only: int32, real32, real64

  use cable_io_vars_module, only: patch_type, land_type

  implicit none
  private

  public :: grid_cell_average
  public :: first_patch_in_grid_cell

  interface first_patch_in_grid_cell
    module procedure first_patch_in_grid_cell_int32_1d
    module procedure first_patch_in_grid_cell_int32_2d
    module procedure first_patch_in_grid_cell_int32_3d
    module procedure first_patch_in_grid_cell_real32_1d
    module procedure first_patch_in_grid_cell_real32_2d
    module procedure first_patch_in_grid_cell_real32_3d
    module procedure first_patch_in_grid_cell_real64_1d
    module procedure first_patch_in_grid_cell_real64_2d
    module procedure first_patch_in_grid_cell_real64_3d
  end interface

  interface grid_cell_average
    module procedure grid_cell_average_real32_1d
    module procedure grid_cell_average_real32_2d
    module procedure grid_cell_average_real32_3d
    module procedure grid_cell_average_real64_1d
    module procedure grid_cell_average_real64_2d
    module procedure grid_cell_average_real64_3d
  end interface

contains

  subroutine grid_cell_average_real32_1d(input_array, output_array, patch, landpt)
    real(kind=real32), intent(in) :: input_array(:)
    real(kind=real32), intent(out) :: output_array(:)
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, patch_index

    do land_index = 1, size(output_array)
      output_array(land_index) = 0.0_real32
      do patch_index = landpt(land_index)%cstart, landpt(land_index)%cend
        output_array(land_index) = output_array(land_index) + &
              input_array(patch_index) * patch(patch_index)%frac
      end do
    end do

  end subroutine

  subroutine grid_cell_average_real32_2d(input_array, output_array, patch, landpt)
    real(kind=real32), intent(in) :: input_array(:, :)
    real(kind=real32), intent(out) :: output_array(:, :)
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, patch_index, j

    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        output_array(land_index, j) = 0.0_real32
        do patch_index = landpt(land_index)%cstart, landpt(land_index)%cend
          output_array(land_index, j) = ( &
            output_array(land_index, j) + input_array(patch_index, j) * patch(patch_index)%frac &
          )
        end do
      end do
    end do

  end subroutine

  subroutine grid_cell_average_real32_3d(input_array, output_array, patch, landpt)
    real(kind=real32), intent(in) :: input_array(:, :, :)
    real(kind=real32), intent(out) :: output_array(:, :, :)
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, patch_index, j, k

    do k = 1, size(output_array, 3)
      do j = 1, size(output_array, 2)
        do land_index = 1, size(output_array, 1)
          output_array(land_index, j, k) = 0.0_real32
          do patch_index = landpt(land_index)%cstart, landpt(land_index)%cend
            output_array(land_index, j, k) = ( &
              output_array(land_index, j, k) + &
              input_array(patch_index, j, k) * patch(patch_index)%frac &
            )
          end do
        end do
      end do
    end do

  end subroutine

  subroutine grid_cell_average_real64_1d(input_array, output_array, patch, landpt)
    real(kind=real64), intent(in) :: input_array(:)
    real(kind=real64), intent(out) :: output_array(:)
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, patch_index

    do land_index = 1, size(output_array)
      output_array(land_index) = 0.0_real64
      do patch_index = landpt(land_index)%cstart, landpt(land_index)%cend
        output_array(land_index) = output_array(land_index) + &
              input_array(patch_index) * patch(patch_index)%frac
      end do
    end do

  end subroutine

  subroutine grid_cell_average_real64_2d(input_array, output_array, patch, landpt)
    real(kind=real64), intent(in) :: input_array(:, :)
    real(kind=real64), intent(out) :: output_array(:, :)
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, patch_index, j

    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        output_array(land_index, j) = 0.0_real64
        do patch_index = landpt(land_index)%cstart, landpt(land_index)%cend
          output_array(land_index, j) = ( &
            output_array(land_index, j) + input_array(patch_index, j) * patch(patch_index)%frac &
          )
        end do
      end do
    end do

  end subroutine

  subroutine grid_cell_average_real64_3d(input_array, output_array, patch, landpt)
    real(kind=real64), intent(in) :: input_array(:, :, :)
    real(kind=real64), intent(out) :: output_array(:, :, :)
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, patch_index, j, k

    do k = 1, size(output_array, 3)
      do j = 1, size(output_array, 2)
        do land_index = 1, size(output_array, 1)
          output_array(land_index, j, k) = 0.0_real64
          do patch_index = landpt(land_index)%cstart, landpt(land_index)%cend
            output_array(land_index, j, k) = ( &
              output_array(land_index, j, k) + &
              input_array(patch_index, j, k) * patch(patch_index)%frac &
            )
          end do
        end do
      end do
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_int32_1d(input_array, output_array, landpt)
    integer(kind=int32), intent(in) :: input_array(:)
    integer(kind=int32), intent(out) :: output_array(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index

    do land_index = 1, size(output_array)
      output_array(land_index) = input_array(landpt(land_index)%cstart)
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_int32_2d(input_array, output_array, landpt)
    integer(kind=int32), intent(in) :: input_array(:, :)
    integer(kind=int32), intent(out) :: output_array(:, :)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, j

    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        output_array(land_index, j) = input_array(landpt(land_index)%cstart, j)
      end do
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_int32_3d(input_array, output_array, landpt)
    integer(kind=int32), intent(in) :: input_array(:, :, :)
    integer(kind=int32), intent(out) :: output_array(:, :, :)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, j, k

    do k = 1, size(output_array, 3)
      do j = 1, size(output_array, 2)
        do land_index = 1, size(output_array, 1)
          output_array(land_index, j, k) = input_array(landpt(land_index)%cstart, j, k)
        end do
      end do
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real32_1d(input_array, output_array, landpt)
    real(kind=real32), intent(in) :: input_array(:)
    real(kind=real32), intent(out) :: output_array(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index

    do land_index = 1, size(output_array)
      output_array(land_index) = input_array(landpt(land_index)%cstart)
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real32_2d(input_array, output_array, landpt)
    real(kind=real32), intent(in) :: input_array(:, :)
    real(kind=real32), intent(out) :: output_array(:, :)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, j

    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        output_array(land_index, j) = input_array(landpt(land_index)%cstart, j)
      end do
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real32_3d(input_array, output_array, landpt)
    real(kind=real32), intent(in) :: input_array(:, :, :)
    real(kind=real32), intent(out) :: output_array(:, :, :)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, j, k

    do k = 1, size(output_array, 3)
      do j = 1, size(output_array, 2)
        do land_index = 1, size(output_array, 1)
          output_array(land_index, j, k) = input_array(landpt(land_index)%cstart, j, k)
        end do
      end do
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real64_1d(input_array, output_array, landpt)
    real(kind=real64), intent(in) :: input_array(:)
    real(kind=real64), intent(out) :: output_array(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index

    do land_index = 1, size(output_array)
      output_array(land_index) = input_array(landpt(land_index)%cstart)
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real64_2d(input_array, output_array, landpt)
    real(kind=real64), intent(in) :: input_array(:, :)
    real(kind=real64), intent(out) :: output_array(:, :)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, j

    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        output_array(land_index, j) = input_array(landpt(land_index)%cstart, j)
      end do
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real64_3d(input_array, output_array, landpt)
    real(kind=real64), intent(in) :: input_array(:, :, :)
    real(kind=real64), intent(out) :: output_array(:, :, :)
    type(land_type), intent(in) :: landpt(:)
    integer :: land_index, j, k

    do k = 1, size(output_array, 3)
      do j = 1, size(output_array, 2)
        do land_index = 1, size(output_array, 1)
          output_array(land_index, j, k) = input_array(landpt(land_index)%cstart, j, k)
        end do
      end do
    end do

  end subroutine

end module
