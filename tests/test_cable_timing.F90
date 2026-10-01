! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_cable_timing
  !! Tests for the output frequencies in cable_timing_mod.
  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_case
  use fortuno_interface_mod, only: check
  use cable_common_module, only: CurYear
  use cable_io_vars_module, only: leaps
  use cable_timing_mod, only: cable_timing_frequency_matches
  use cable_timing_mod, only: cable_timing_frequency_rank
  use cable_timing_mod, only: cable_timing_date_string
  use cable_timing_mod, only: cable_timing_set_start_year
  implicit none
  private

  public :: cable_timing_test_list

  real, parameter :: dels = 1800.0
    !! Half-hour time steps: 48 per day, 6 per 3 hours.

contains

  function cable_timing_test_list() result(test_list)
    type(test_list_t) :: test_list

    test_list = test_list_t([ &
      test_case("test_timing_timestep", test_timestep), &
      test_case("test_timing_3hrly_and_daily", test_three_hourly_and_daily), &
      test_case("test_timing_monthly", test_monthly), &
      test_case("test_timing_yearly", test_yearly), &
      test_case("test_timing_frequency_rank", test_frequency_rank), &
      test_case("test_timing_date_string", test_date_string) &
    ])
  end function cable_timing_test_list

  subroutine test_timestep()
    call check(cable_timing_frequency_matches(dels, 1, "timestep"), msg="every step matches")
    call check(cable_timing_frequency_matches(dels, 12345, "timestep"), msg="any step matches")
  end subroutine test_timestep

  subroutine test_three_hourly_and_daily()
    call check(cable_timing_frequency_matches(dels, 6, "3hrly"), msg="first 3 hours")
    call check(cable_timing_frequency_matches(dels, 12, "3hrly"), msg="second 3 hours")
    call check(.not. cable_timing_frequency_matches(dels, 7, "3hrly"), msg="between")
    call check(cable_timing_frequency_matches(dels, 48, "daily"), msg="end of day 1")
    call check(.not. cable_timing_frequency_matches(dels, 47, "daily"), msg="before end of day 1")
    call check(cable_timing_frequency_matches(dels, 96, "daily"), msg="end of day 2")
  end subroutine test_three_hourly_and_daily

  subroutine test_monthly()
    call cable_timing_set_start_year(2001)
    CurYear = 2001
    leaps = .false.
    call check(cable_timing_frequency_matches(dels, 31 * 48, "monthly"), msg="end of January")
    call check(.not. cable_timing_frequency_matches(dels, 31 * 48 - 1, "monthly"), msg="before the end of January")
    call check(cable_timing_frequency_matches(dels, 59 * 48, "monthly"), msg="end of February, no leap year")
    call check(cable_timing_frequency_matches(dels, 365 * 48, "monthly"), msg="end of December")
    CurYear = 2002
    call check(cable_timing_frequency_matches(dels, (365 + 31) * 48, "monthly"), msg="end of January in the second year")
  end subroutine test_monthly

  subroutine test_yearly()
    call cable_timing_set_start_year(2001)
    CurYear = 2001
    leaps = .false.
    call check(cable_timing_frequency_matches(dels, 365 * 48, "yearly"), msg="end of the first year")
    call check(.not. cable_timing_frequency_matches(dels, 365 * 48 - 1, "yearly"), msg="just before the end of the year")
    call check(.not. cable_timing_frequency_matches(dels, 31 * 48, "yearly"), msg="not at the end of January")
    CurYear = 2002
    call check(cable_timing_frequency_matches(dels, 730 * 48, "yearly"), msg="end of the second year")
    ! With leap years the year 2004 has 366 days.
    leaps = .true.
    CurYear = 2004
    call cable_timing_set_start_year(2004)
    call check(cable_timing_frequency_matches(dels, 366 * 48, "yearly"), msg="leap year has 366 days")
    leaps = .false.
  end subroutine test_yearly

  subroutine test_frequency_rank()
    call check(cable_timing_frequency_rank("timestep") == 1, msg="timestep is finest")
    call check(cable_timing_frequency_rank("3hrly") == 2, msg="3hrly")
    call check(cable_timing_frequency_rank("daily") == 3, msg="daily")
    call check(cable_timing_frequency_rank("monthly") == 4, msg="monthly")
    call check(cable_timing_frequency_rank("yearly") == 5, msg="yearly is coarsest")
    call check(cable_timing_frequency_rank("all") == 0, msg="the old name is not accepted")
    call check(cable_timing_frequency_rank("user006") == 0, msg="custom intervals are not accepted")
  end subroutine test_frequency_rank

  subroutine test_date_string()
    call check(cable_timing_date_string(2000, 1, 0, .true.) == "2000-01-01", msg="start date")
    call check(cable_timing_date_string(2000, 1, 30, .true.) == "2000-01-31", msg="end of January")
    call check(cable_timing_date_string(2000, 1, 31, .true.) == "2000-02-01", msg="start of February")
    call check(cable_timing_date_string(2000, 60, 0, .true.) == "2000-02-29", msg="leap day")
    call check(cable_timing_date_string(2000, 60, 0, .false.) == "2000-03-01", msg="no leap day without leap years")
    call check(cable_timing_date_string(2001, 365, 0, .true.) == "2001-12-31", msg="last day of a common year")
    call check(cable_timing_date_string(2001, 365, 1, .true.) == "2002-01-01", msg="into the next year")
    call check(cable_timing_date_string(2000, 1, 366, .true.) == "2001-01-01", msg="a leap year has 366 days")
    call check(cable_timing_date_string(2000, 1, 365, .false.) == "2001-01-01", msg="a year has 365 days without leap years")
    call check(cable_timing_date_string(2002, 1, 16 * 365 + 4 - 1, .true.) == "2017-12-31", msg="16 years including leap days")
  end subroutine test_date_string

end module test_cable_timing
