! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_cable_grid_reductions
  !! Tests for the grid cell reductions.
  use iso_fortran_env, only: int32, real32, real64
  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_case
  use fortuno_interface_mod, only: check
  use cable_io_vars_module, only: patch_type, land_type
  use cable_grid_reductions_mod, only: grid_cell_average
  use cable_grid_reductions_mod, only: first_tile_on_cell
  use cable_grid_reductions_mod, only: dominant_tile
  implicit none
  private

  public :: cable_grid_reductions_test_list

contains

  function cable_grid_reductions_test_list() result(test_list)
    type(test_list_t) :: test_list

    test_list = test_list_t([ &
      test_case("test_reductions_first_tile_on_cell", test_first_tile_on_cell), &
      test_case("test_reductions_dominant_tile", test_dominant_tile), &
      test_case("test_reductions_dominant_tile_tie_takes_first", test_dominant_tile_tie), &
      test_case("test_reductions_dominant_tile_2d_and_3d", test_dominant_tile_2d_and_3d), &
      test_case("test_reductions_dominant_tile_integer", test_dominant_tile_integer), &
      test_case("test_reductions_grid_cell_average", test_grid_cell_average) &
    ])
  end function cable_grid_reductions_test_list

  subroutine make_grid(patch, landpt)
    !! Three cells: tiles 1-3 (fractions .2 .5 .3), tile 4 (1.0), tiles 5-6 (.6 .4).
    type(patch_type), intent(out) :: patch(6)
    type(land_type), intent(out) :: landpt(3)

    patch(:)%frac = [0.2, 0.5, 0.3, 1.0, 0.6, 0.4]
    landpt(1)%cstart = 1; landpt(1)%cend = 3
    landpt(2)%cstart = 4; landpt(2)%cend = 4
    landpt(3)%cstart = 5; landpt(3)%cend = 6
  end subroutine make_grid

  subroutine test_first_tile_on_cell()
    type(patch_type) :: patch(6)
    type(land_type) :: landpt(3)
    real(kind=real32) :: values(6), reduced(3)

    call make_grid(patch, landpt)
    values = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    call first_tile_on_cell(values, reduced, landpt=landpt)
    call check(all(reduced == [10.0, 40.0, 50.0]), msg="first tile of each cell")
  end subroutine test_first_tile_on_cell

  subroutine test_dominant_tile()
    type(patch_type) :: patch(6)
    type(land_type) :: landpt(3)
    real(kind=real32) :: values(6), reduced(3)
    real(kind=real64) :: values64(6), reduced64(3)

    call make_grid(patch, landpt)
    values = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    call dominant_tile(values, reduced, patch, landpt)
    call check(all(reduced == [20.0, 40.0, 50.0]), msg="largest tile of each cell")
    values64 = real(values, real64)
    call dominant_tile(values64, reduced64, patch, landpt)
    call check(all(reduced64 == [20.0_real64, 40.0_real64, 50.0_real64]), msg="real64")
  end subroutine test_dominant_tile

  subroutine test_dominant_tile_tie()
    type(patch_type) :: patch(6)
    type(land_type) :: landpt(3)
    real(kind=real32) :: values(6), reduced(3)

    call make_grid(patch, landpt)
    patch(1:3)%frac = [0.4, 0.4, 0.2]
    values = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    call dominant_tile(values, reduced, patch, landpt)
    call check(reduced(1) == 10.0, msg="the first of equal tiles is used")
  end subroutine test_dominant_tile_tie

  subroutine test_dominant_tile_2d_and_3d()
    type(patch_type) :: patch(6)
    type(land_type) :: landpt(3)
    real(kind=real32) :: values2(6, 2), reduced2(3, 2), values3(6, 2, 2), reduced3(3, 2, 2)
    integer :: j, k

    call make_grid(patch, landpt)
    do j = 1, 2
      values2(:, j) = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0] * j
      do k = 1, 2
        values3(:, j, k) = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0] * j * k
      end do
    end do
    call dominant_tile(values2, reduced2, patch, landpt)
    call check(all(reduced2(:, 1) == [20.0, 40.0, 50.0]) .and. all(reduced2(:, 2) == [40.0, 80.0, 100.0]), msg="2d")
    call dominant_tile(values3, reduced3, patch, landpt)
    call check(all(reduced3(:, 2, 2) == [80.0, 160.0, 200.0]), msg="3d")
  end subroutine test_dominant_tile_2d_and_3d

  subroutine test_dominant_tile_integer()
    type(patch_type) :: patch(6)
    type(land_type) :: landpt(3)
    integer(kind=int32) :: values(6), reduced(3)

    call make_grid(patch, landpt)
    values = [1, 2, 3, 4, 5, 6]
    call dominant_tile(values, reduced, patch, landpt)
    call check(all(reduced == [2, 4, 5]), msg="integers are reduced without averaging")
  end subroutine test_dominant_tile_integer

  subroutine test_grid_cell_average()
    type(patch_type) :: patch(6)
    type(land_type) :: landpt(3)
    real(kind=real32) :: values(6), reduced(3)

    call make_grid(patch, landpt)
    values = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    call grid_cell_average(values, reduced, patch, landpt)
    call check(abs(reduced(1) - (2.0 + 10.0 + 9.0)) < 1.0e-5, msg="area-weighted average of cell 1")
    call check(abs(reduced(2) - 40.0) < 1.0e-5, msg="single tile")
    call check(abs(reduced(3) - (30.0 + 24.0)) < 1.0e-5, msg="area-weighted average of cell 3")
  end subroutine test_grid_cell_average

end module test_cable_grid_reductions
