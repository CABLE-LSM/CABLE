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
  public :: cable_timing_frequency_rank
  public :: cable_timing_date_string
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

  !> Frequencies at which output can be written, from finest to coarsest.
  character(len=16), parameter, public :: cable_timing_frequencies(5) = [character(len=16) :: &
    "timestep", "3hrly", "daily", "monthly", "yearly" &
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

  pure function cable_timing_frequency_rank(frequency) result(rank)
    !! Position of `frequency` in `cable_timing_frequencies`, from finest (1) to
    !! coarsest, or 0 if it is not a known frequency.
    character(len=*), intent(in) :: frequency
    integer :: rank
    integer :: candidate

    ! The list is ordered finest to coarsest, so a frequency's position in it
    ! is its rank. Comparing two ranks tells which of two frequencies is finer.
    ! Zero means "not in the list", which callers use to spot a bad name.
    rank = 0
    do candidate = 1, size(cable_timing_frequencies)
      if (trim(cable_timing_frequencies(candidate)) == trim(frequency)) then
        rank = candidate
        return
      end if
    end do
  end function cable_timing_frequency_rank

  function cable_timing_date_string(start_year, start_day_of_year, elapsed_days, leap_years) result(date)
    !! The date `elapsed_days` after the given day of the given year, as `YYYY-MM-DD`.
    integer, intent(in) :: start_year !! Year of the starting date.
    integer, intent(in) :: start_day_of_year !! Day of the year of the starting date, from 1.
    integer, intent(in) :: elapsed_days !! Days after the starting date; not negative.
    logical, intent(in) :: leap_years !! Whether leap years have 366 days.
    character(len=10) :: date
    integer :: year, day_of_year, month, cumulative_days(months_in_year)

    if (elapsed_days < 0) call cable_abort("Elapsed days must not be negative", __FILE__, __LINE__)
    ! Work with a day-of-year count. Add the elapsed days, which may run past
    ! the end of the year, and then roll whole years off the front until what
    ! is left fits inside a single year.
    year = start_year
    day_of_year = start_day_of_year + elapsed_days
    do while (day_of_year > days_in_year(year, leap_years))
      day_of_year = day_of_year - days_in_year(year, leap_years)
      year = year + 1
    end do

    ! Turn the day of the year into month and day. `cumulative_days(m)` is how
    ! many days have passed by the end of month m, so the month is the first
    ! one whose end is on or after the day. The day of the month is what is
    ! left after taking off the days of all the earlier months.
    cumulative_days = last_day
    if (leap_years .and. is_leapyear(year)) cumulative_days = last_day_leap
    month = 1
    do while (day_of_year > cumulative_days(month))
      month = month + 1
    end do
    if (month > 1) day_of_year = day_of_year - cumulative_days(month - 1)
    ! i4.4 / i2.2 pad with zeros: year 2000, month 1, day 5 -> 2000-01-05.
    write(date, '(i4.4,"-",i2.2,"-",i2.2)') year, month, day_of_year
  end function cable_timing_date_string

  pure function days_in_year(year, leap_years) result(days)
    integer, intent(in) :: year
    logical, intent(in) :: leap_years
    integer :: days

    ! A year is 365 days unless leap years are switched on and this is one.
    days = 365
    if (leap_years .and. is_leapyear(year)) days = 366
  end function days_in_year

  function cable_timing_frequency_matches(dels, ktau, frequency) result(match)
    !! Determines whether the current time step is the end of a period of the
    !! given frequency.
    real, intent(in) :: dels !! Model time step in seconds
    integer, intent(in) :: ktau !! Current time step index
    character(len=*), intent(in) :: frequency
      !! One of `cable_timing_frequencies`: `timestep`, `3hrly`, `daily`, `monthly` or `yearly`.
    logical :: match
    integer :: i, time_steps_per_interval
    integer :: last_day_of_month_in_total_elapsed_days(months_in_year)
    integer :: elapsed_days_at_end_of_year

    select case (trim(frequency))
    ! A period ends when the step count is a whole multiple of the number of
    ! steps the period contains (or, for months and years, on a known step).
    case ('timestep')
      ! Every step is the end of a period.
      match = .true.
    case ('3hrly')
      ! Steps per 3 hours = seconds in 3 hours / seconds per step.
      time_steps_per_interval = seconds_per_hour * 3 / int(dels)
      match = mod(ktau, time_steps_per_interval) == 0
    case ('daily')
      time_steps_per_interval = seconds_per_hour * hours_per_day / int(dels)
      match = mod(ktau, time_steps_per_interval) == 0
    case ('monthly')
      if (start_year == START_YEAR_UNDEFINED) then
        call cable_abort('start_year undefined for monthly frequency', __FILE__, __LINE__)
      end if
      ! Build the 12 month-end positions of the current year as a count of days
      ! since the start of the run: first the days in all the earlier years...
      last_day_of_month_in_total_elapsed_days = 0
      do i = start_year, current_year - 1
        if (leaps .and. is_leapyear(i)) then
          last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + 366
        else
          last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + 365
        end if
      end do
      ! ...then add the running total of days at the end of each month of this
      ! year (a whole array is added at once, one value per month).
      if (leaps .and. is_leapyear(current_year)) then
        last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + last_day_leap
      else
        last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + last_day
      end if
      ! This step ends a month if it is the step that lands on any of those days.
      match = any(int(real(last_day_of_month_in_total_elapsed_days) * seconds_per_day / dels) == ktau)
    case ('yearly')
      if (start_year == START_YEAR_UNDEFINED) then
        call cable_abort('start_year undefined for yearly frequency', __FILE__, __LINE__)
      end if
      ! Count the days from the start of the run to the end of the current year
      ! by adding up the length of each year in turn (366 for leap years when
      ! they are in use). The year ends on the step that reaches that many days.
      elapsed_days_at_end_of_year = 0
      do i = start_year, current_year
        if (leaps .and. is_leapyear(i)) then
          elapsed_days_at_end_of_year = elapsed_days_at_end_of_year + 366
        else
          elapsed_days_at_end_of_year = elapsed_days_at_end_of_year + 365
        end if
      end do
      ! Convert days to a step number: days * seconds per day / seconds per step.
      match = int(real(elapsed_days_at_end_of_year) * seconds_per_day / dels) == ktau
    case default
      call cable_abort('Error: unknown frequency "' // trim(adjustl(frequency)) // '"', __FILE__, __LINE__)
    end select

  end function cable_timing_frequency_matches

end module
