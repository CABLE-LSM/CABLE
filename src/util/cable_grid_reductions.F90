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
  public :: first_tile_on_cell
  public :: dominant_tile

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

  interface first_tile_on_cell
    !* Interface for extracting the value from the first patch/tile in each grid
    ! cell for various data types and array ranks. This is useful for arrays where
    ! averaging along the patch/tile dimension does not make sense, or where the
    ! array contains the same value everywhere along the patch/tile dimension.
    module procedure first_tile_on_cell_int32_1d
    module procedure first_tile_on_cell_int32_2d
    module procedure first_tile_on_cell_int32_3d
    module procedure first_tile_on_cell_real32_1d
    module procedure first_tile_on_cell_real32_2d
    module procedure first_tile_on_cell_real32_3d
    module procedure first_tile_on_cell_real64_1d
    module procedure first_tile_on_cell_real64_2d
    module procedure first_tile_on_cell_real64_3d
  end interface

  interface dominant_tile
    !* Interface for extracting the value from the tile with the largest area
    ! fraction in each grid cell for various data types and array ranks. If
    ! several tiles share the largest fraction, the first of them is used.
    module procedure dominant_tile_int32_1d
    module procedure dominant_tile_int32_2d
    module procedure dominant_tile_int32_3d
    module procedure dominant_tile_real32_1d
    module procedure dominant_tile_real32_2d
    module procedure dominant_tile_real32_3d
    module procedure dominant_tile_real64_1d
    module procedure dominant_tile_real64_2d
    module procedure dominant_tile_real64_3d
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

  subroutine first_tile_on_cell_int32_1d(input_array, output_array, landpt)
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

  subroutine first_tile_on_cell_int32_2d(input_array, output_array, landpt)
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

  subroutine first_tile_on_cell_int32_3d(input_array, output_array, landpt)
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

  subroutine first_tile_on_cell_real32_1d(input_array, output_array, landpt)
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

  subroutine first_tile_on_cell_real32_2d(input_array, output_array, landpt)
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

  subroutine first_tile_on_cell_real32_3d(input_array, output_array, landpt)
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

  subroutine first_tile_on_cell_real64_1d(input_array, output_array, landpt)
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

  subroutine first_tile_on_cell_real64_2d(input_array, output_array, landpt)
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

  subroutine first_tile_on_cell_real64_3d(input_array, output_array, landpt)
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

  pure function dominant_tile_index(land_entry, patch) result(dominant_index)
    !! Index of the tile with the largest area fraction in the grid cell
    !! described by `land_entry`. Ties go to the first such tile.
    type(land_type), intent(in) :: land_entry
    type(patch_type), intent(in) :: patch(:)
    integer :: dominant_index
    integer :: patch_index

    ! Start with the cell's first tile as the leader and let each later tile
    ! take over only if it is strictly larger. Because a later tile that merely
    ! equals the leader does not replace it, ties go to the first tile, which
    ! keeps the choice stable.
    dominant_index = land_entry%cstart
    do patch_index = land_entry%cstart + 1, land_entry%cend
      if (patch(patch_index)%frac > patch(dominant_index)%frac) dominant_index = patch_index
    end do
  end function dominant_tile_index

  subroutine dominant_tile_int32_1d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 1D int32 array.
    integer(kind=int32), intent(in) :: input_array(:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    integer(kind=int32), intent(out) :: output_array(:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do land_index = 1, size(output_array, 1)
      dominant_index = dominant_tile_index(landpt(land_index), patch)
      output_array(land_index) = input_array(dominant_index)
    end do
  end subroutine dominant_tile_int32_1d

  subroutine dominant_tile_int32_2d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 2D int32 array.
    integer(kind=int32), intent(in) :: input_array(:,:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    integer(kind=int32), intent(out) :: output_array(:,:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index, j

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        dominant_index = dominant_tile_index(landpt(land_index), patch)
        output_array(land_index, j) = input_array(dominant_index, j)
      end do
    end do
  end subroutine dominant_tile_int32_2d

  subroutine dominant_tile_int32_3d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 3D int32 array.
    integer(kind=int32), intent(in) :: input_array(:,:,:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    integer(kind=int32), intent(out) :: output_array(:,:,:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index, j, k

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do k = 1, size(output_array, 3)
      do j = 1, size(output_array, 2)
        do land_index = 1, size(output_array, 1)
          dominant_index = dominant_tile_index(landpt(land_index), patch)
          output_array(land_index, j, k) = input_array(dominant_index, j, k)
        end do
      end do
    end do
  end subroutine dominant_tile_int32_3d

  subroutine dominant_tile_real32_1d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 1D real32 array.
    real(kind=real32), intent(in) :: input_array(:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do land_index = 1, size(output_array, 1)
      dominant_index = dominant_tile_index(landpt(land_index), patch)
      output_array(land_index) = input_array(dominant_index)
    end do
  end subroutine dominant_tile_real32_1d

  subroutine dominant_tile_real32_2d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 2D real32 array.
    real(kind=real32), intent(in) :: input_array(:,:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:,:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index, j

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        dominant_index = dominant_tile_index(landpt(land_index), patch)
        output_array(land_index, j) = input_array(dominant_index, j)
      end do
    end do
  end subroutine dominant_tile_real32_2d

  subroutine dominant_tile_real32_3d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 3D real32 array.
    real(kind=real32), intent(in) :: input_array(:,:,:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real32), intent(out) :: output_array(:,:,:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index, j, k

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do k = 1, size(output_array, 3)
      do j = 1, size(output_array, 2)
        do land_index = 1, size(output_array, 1)
          dominant_index = dominant_tile_index(landpt(land_index), patch)
          output_array(land_index, j, k) = input_array(dominant_index, j, k)
        end do
      end do
    end do
  end subroutine dominant_tile_real32_3d

  subroutine dominant_tile_real64_1d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 1D real64 array.
    real(kind=real64), intent(in) :: input_array(:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real64), intent(out) :: output_array(:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do land_index = 1, size(output_array, 1)
      dominant_index = dominant_tile_index(landpt(land_index), patch)
      output_array(land_index) = input_array(dominant_index)
    end do
  end subroutine dominant_tile_real64_1d

  subroutine dominant_tile_real64_2d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 2D real64 array.
    real(kind=real64), intent(in) :: input_array(:,:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real64), intent(out) :: output_array(:,:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index, j

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do j = 1, size(output_array, 2)
      do land_index = 1, size(output_array, 1)
        dominant_index = dominant_tile_index(landpt(land_index), patch)
        output_array(land_index, j) = input_array(dominant_index, j)
      end do
    end do
  end subroutine dominant_tile_real64_2d

  subroutine dominant_tile_real64_3d(input_array, output_array, patch, landpt)
    !! Extracts the value of the largest-area tile in each grid cell from a 3D real64 array.
    real(kind=real64), intent(in) :: input_array(:,:,:)
      !* The input array to be reduced. The first (i.e. fastest varying) dimension of
      ! this array must be the patch/tile dimension being reduced.
    real(kind=real64), intent(out) :: output_array(:,:,:)
      !* The output array containing the reduced per grid cell values. The first
      ! (i.e. fastest varying) dimension of this array must be equal to the number
      ! of grid cells.
    type(patch_type), intent(in) :: patch(:)
      !! The `patch_type` instance describing the area fraction of each patch/tile.
    type(land_type), intent(in) :: landpt(:)
      !* The `land_type` instance describing the starting and ending patch/tile
      ! indexes in the input array for each grid cell.
    integer :: land_index, dominant_index, j, k

    ! Each output point is one grid cell. Find that cell's largest tile, then copy
    ! that tile's value across. Any further dimensions (soil layers and so on)
    ! are simply looped over, so the same tile is chosen for every layer.
    do k = 1, size(output_array, 3)
      do j = 1, size(output_array, 2)
        do land_index = 1, size(output_array, 1)
          dominant_index = dominant_tile_index(landpt(land_index), patch)
          output_array(land_index, j, k) = input_array(dominant_index, j, k)
        end do
      end do
    end do
  end subroutine dominant_tile_real64_3d

end module
