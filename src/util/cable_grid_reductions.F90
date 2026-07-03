! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_grid_reductions_mod
  !* This module provides procedures for performing various grid cell reductions
  ! for data along some dimension. This is commonly used for reducing data along
  ! the tile/patch dimension to a per grid cell value.

  use iso_fortran_env, only: int32, real32, real64

  use cable_io_vars_module, only: patch_type, land_type

  implicit none
  private

  public :: grid_cell_average
  public :: first_patch_in_grid_cell

  interface grid_cell_average
    !* Interface for computing the area weighted average over the patch/tile
    ! dimension for various data types and array ranks.
    module procedure grid_cell_average_real32_1d
    module procedure grid_cell_average_real32_2d
    module procedure grid_cell_average_real32_3d
    module procedure grid_cell_average_real64_1d
    module procedure grid_cell_average_real64_2d
    module procedure grid_cell_average_real64_3d
  end interface

  interface first_patch_in_grid_cell
    !* Interface for extracting the value from the first patch/tile in each grid
    ! cell for various data types and array ranks. This is useful for arrays where
    ! averaging along the patch/tile dimension does not make sense, or where the
    ! array contains the same value everywhere along the patch/tile dimension.
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

contains

  subroutine grid_cell_average_real32_1d(input_array, output_array, patch, landpt)
    !* Computes the area weighted average over the patch/tile dimension for a 1D
    ! 32-bit real array.
    real(kind=real32), intent(in) :: input_array(:)
      !* The input array to be reduced. The first (i.e. fastest varying)
      ! dimension of this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:)
      !* The output array containing the grid cell averaged values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !* The `patch_type` instance describing the area fraction of each active
      ! patch/tile dimension.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
    !* Computes the area weighted average over the patch/tile dimension for a 2D
    ! 32-bit real array.
    real(kind=real32), intent(in) :: input_array(:, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:, :)
      !* The output array containing the grid cell averaged values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the
      ! number of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !* The `patch_type` instance describing the area fraction of each active
      ! patch/tile dimension.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
    !* Computes the area weighted average over the patch/tile dimension for a 3D
    ! 32-bit real array.
    real(kind=real32), intent(in) :: input_array(:, :, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:, :, :)
      !* The output array containing the grid cell averaged values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the
      ! number of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !* The `patch_type` instance describing the area fraction of each active
      ! patch/tile dimension.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
    !* Computes the area weighted average over the patch/tile dimension for a 1D
    ! 64-bit real array.
    real(kind=real64), intent(in) :: input_array(:)
      !* The input array to be reduced. The first (i.e. fastest varying)
      ! dimension of this array must be the patch/tile dimension being reduced.
    real(kind=real64), intent(out) :: output_array(:)
      !* The output array containing the grid cell averaged values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !* The `patch_type` instance describing the area fraction of each active
      ! patch/tile dimension.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
    !* Computes the area weighted average over the patch/tile dimension for a 2D
    ! 64-bit real array.
    real(kind=real64), intent(in) :: input_array(:, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real64), intent(out) :: output_array(:, :)
      !* The output array containing the grid cell averaged values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the
      ! number of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !* The `patch_type` instance describing the area fraction of each active
      ! patch/tile dimension.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
    !* Computes the area weighted average over the patch/tile dimension for a 3D
    ! 64-bit real array.
    real(kind=real64), intent(in) :: input_array(:, :, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real64), intent(out) :: output_array(:, :, :)
      !* The output array containing the grid cell averaged values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the
      ! number of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !* The `patch_type` instance describing the area fraction of each active
      ! patch/tile dimension.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
    !! Extracts the first patch value for each grid cell from a 1D integer array.
    integer(kind=int32), intent(in) :: input_array(:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    integer(kind=int32), intent(out) :: output_array(:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index

    do land_index = 1, size(output_array)
      output_array(land_index) = input_array(landpt(land_index)%cstart)
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_int32_2d(input_array, output_array, landpt)
    !! Extracts the first patch value for each grid cell from a 2D integer array.
    integer(kind=int32), intent(in) :: input_array(:, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    integer(kind=int32), intent(out) :: output_array(:, :)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, j

    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        output_array(land_index, j) = input_array(landpt(land_index)%cstart, j)
      end do
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_int32_3d(input_array, output_array, landpt)
    !! Extracts the first patch value for each grid cell from a 3D integer array.
    integer(kind=int32), intent(in) :: input_array(:, :, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    integer(kind=int32), intent(out) :: output_array(:, :, :)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
    !! Extracts the first patch value for each grid cell from a 1D 32-bit real array.
    real(kind=real32), intent(in) :: input_array(:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index

    do land_index = 1, size(output_array)
      output_array(land_index) = input_array(landpt(land_index)%cstart)
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real32_2d(input_array, output_array, landpt)
    !! Extracts the first patch value for each grid cell from a 2D 32-bit real array.
    real(kind=real32), intent(in) :: input_array(:, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:, :)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, j

    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        output_array(land_index, j) = input_array(landpt(land_index)%cstart, j)
      end do
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real32_3d(input_array, output_array, landpt)
    !! Extracts the first patch value for each grid cell from a 3D 32-bit real array.
    real(kind=real32), intent(in) :: input_array(:, :, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:, :, :)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
    !! Extracts the first patch value for each grid cell from a 1D 64-bit real array.
    real(kind=real64), intent(in) :: input_array(:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real64), intent(out) :: output_array(:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index

    do land_index = 1, size(output_array)
      output_array(land_index) = input_array(landpt(land_index)%cstart)
    end do

  end subroutine

  subroutine first_patch_in_grid_cell_real64_2d(input_array, output_array, landpt)
    !! Extracts the first patch value for each grid cell from a 2D 64-bit real array.
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
    !! Extracts the first patch value for each grid cell from a 3D 64-bit real array.
    real(kind=real64), intent(in) :: input_array(:, :, :)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real64), intent(out) :: output_array(:, :, :)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
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
