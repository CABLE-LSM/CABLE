! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_timing_mod
  !! Module for handling timing in CABLE.
  use cable_error_handler_mod, only: cable_abort
  use cable_common_module, only: is_leapyear, current_year => CurYear
  use cable_io_vars_module, only: leaps
  implicit none
  private

  public :: cable_timing_frequency_matches
  public :: cable_timing_frequency_is_greater_than
  public :: cable_timing_set_start_year

  integer, parameter, public :: seconds_per_hour = 3600
  integer, parameter, public :: hours_per_day = 24
  integer, parameter, public :: seconds_per_day = 86400
  integer, parameter, public :: months_in_year = 12

  !> Cumulative day of year at the end of each month for a non-leap year.
  integer, parameter, dimension(months_in_year) :: last_day = [ &
    31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 &
  ]

  !> Cumulative day of year at the end of each month for a leap year.
  integer, parameter, dimension(months_in_year) :: last_day_leap = [ &
    31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 &
  ]

  integer, parameter :: START_YEAR_UNDEFINED = -1

  !> Start year of the simulation.
  integer :: start_year = START_YEAR_UNDEFINED

contains

  subroutine cable_timing_set_start_year(year)
    !* Set the start year of the simulation. This is used for calculating monthly
    ! timing.
    integer, intent(in) :: year
    start_year = year
  end subroutine

  function cable_timing_frequency_matches(dels, ktau, frequency) result(match)
    !! Determines whether the current time step matches the specified frequency.
    real, intent(in) :: dels !! Model time step in seconds
    integer, intent(in) :: ktau !! Current time step index
    character(len=*), intent(in) :: frequency !! Frequency string: 'all', 'user', 'daily', 'monthly'
    logical :: match
    integer :: i, time_steps_per_interval, interval_in_hours
    integer :: last_day_of_month_in_total_elapsed_days(months_in_year)

    select case (frequency)
    case ('user')
      read(frequency(5:7), *) interval_in_hours
      time_steps_per_interval = seconds_per_hour * interval_in_hours / int(dels)
      match = mod(ktau, time_steps_per_interval) == 0
    case ('all')
      match = .true.
    case ('daily')
      time_steps_per_interval = seconds_per_hour * hours_per_day / int(dels)
      match = mod(ktau, time_steps_per_interval) == 0
    case ('monthly')
      if (start_year == START_YEAR_UNDEFINED) then
        call cable_abort('start_year undefined for monthly frequency', __FILE__, __LINE__)
      end if
      last_day_of_month_in_total_elapsed_days = 0
      do i = start_year, current_year - 1
        if (leaps .and. is_leapyear(i)) then
          last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + 366
        else
          last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + 365
        end if
      end do
      if (leaps .and. is_leapyear(current_year)) then
        last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + last_day_leap
      else
        last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + last_day
      end if
      match = any(int(real(last_day_of_month_in_total_elapsed_days) * seconds_per_day / dels) == ktau)
    case default
      call cable_abort('Error: unknown frequency "' // trim(adjustl(frequency)) // '"', __FILE__, __LINE__)
    end select

  end function

  logical function cable_timing_frequency_is_greater_than(freq_a, freq_b, dels) result(freq_a_greater_than_b)
    !* Utility function to determine whether one frequency is greater than
    ! another following the ordering "all" > "user" > "daily" > "monthly".
    character(len=*), intent(in) :: freq_a
      !! The first frequency to compare, one of "all", "user", "daily", or "monthly".
    character(len=*), intent(in) :: freq_b
      !! The second frequency to compare, one of "all", "user", "daily", or "monthly".
    real, intent(in) :: dels
      !! Model time step in seconds, used for comparing against "user" frequencies.

    integer :: period_in_hours_a, period_in_hours_b

    select case (freq_a)
    case ("all")
      if (freq_b == "all") then
        freq_a_greater_than_b = .false.
      else if (freq_b == "user") then
        read(freq_b(5:7), *) period_in_hours_b
        freq_a_greater_than_b = dels / seconds_per_hour < period_in_hours_b
      else
        freq_a_greater_than_b = .true.
      end if
    case ("user")
      read(freq_a(5:7), *) period_in_hours_a
      if (freq_b == "user") then
        read(freq_b(5:7), *) period_in_hours_b
        freq_a_greater_than_b = period_in_hours_a < period_in_hours_b
      else if (freq_b == "all") then
        freq_a_greater_than_b = period_in_hours_a < dels / seconds_per_hour
      else
        freq_a_greater_than_b = .true.
      end if
    case ("daily")
      freq_a_greater_than_b = freq_b == "monthly"
    case ("monthly")
      freq_a_greater_than_b = .false.
    case default
      call cable_abort("Unexpected sampling frequency '" // freq_a, __FILE__, __LINE__)
    end select

  end function

end module
