! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_array_utils_mod
  !! Utility procedures for working with arrays.
  implicit none
  private

  public array_offset
  public array_index
  public array_partition

contains

  !> Calculate the memory offset corresponding to a given index in a multi-dimensional array.
  function array_offset(index, shape) result(offset)
    integer, intent(in) :: index(:) !! The multi-dimensional index for which to calculate the offset.
    integer, intent(in) :: shape(:) !! The shape of the multi-dimensional array.
    integer :: offset !! The calculated memory offset corresponding to the given index.
    integer :: i, scale
    offset = 1; scale = 1
    do i = 1, size(index)
      offset = offset + (index(i) - 1) * scale
      scale = scale * shape(i)
    end do
  end function

  !> Calculate the index corresponding to a given memory offset in a multi-dimensional array.
  subroutine array_index(offset_in, shape, index)
    integer, intent(in) :: offset_in !! The memory offset for which to calculate the multi-dimensional index.
    integer, intent(in) :: shape(:) !! The shape of the multi-dimensional array.
    integer, intent(inout) :: index(:) !! The calculated multi-dimensional index corresponding to the given offset.
    integer :: i, offset
    offset = offset_in
    do i = 1, size(shape)
      index(i) = mod(offset - 1, shape(i)) + 1
      offset = (offset - 1) / shape(i) + 1
    end do
  end subroutine

  subroutine array_partition(n, k, p, start, count)
    !* Compute start and count for the `p`'th partition of an array of size `n`
    ! where `p = 0, 1, ... , k - 1`.
    !
    ! For `k` partitions, an array of `n` elements can be partitioned into `r`
    ! partitions of length `q + 1`, and `k - r` partitions of length `q` where `q`
    ! and `r` are the quotient and remainder of `n` divided by `k` (i.e. `n = q *
    ! k + r`).  Note, we assume that the `r` partitions of length `q + 1` precede
    ! the `k - r` partitions of length `q` in the array.
    integer, intent(in) :: n !! The total number of elements in the array to be partitioned.
    integer, intent(in) :: k !! The total number of partitions.
    integer, intent(in) :: p !! The index of the partition for which to compute the start and count (0-based).
    integer, intent(out) :: start !! The starting index (1-based) of the `p`'th partition in the array.
    integer, intent(out) :: count !! The number of elements in the `p`'th partition.
    integer :: q, r

    q = n / k
    r = mod(n, k)

    if (p < r) then
      count = q + 1
      start = 1 + (q + 1) * p
    else
      count = q
      start = 1 + (q + 1) * r + q * (p - r)
    end if

  end subroutine array_partition

end module cable_array_utils_mod
