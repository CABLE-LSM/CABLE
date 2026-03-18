module cable_timing_utils_mod
  use cable_error_handler_mod, only: cable_abort
  use cable_common_module, only: is_leapyear, current_year => CurYear
  implicit none
  private

  public :: time_step_matches
  public :: seconds_per_hour
  public :: hours_per_day
  public :: seconds_per_day

  integer, parameter :: seconds_per_hour = 3600
  integer, parameter :: hours_per_day = 24
  integer, parameter :: seconds_per_day = 86400
  integer, parameter :: months_in_year = 12
  integer, parameter, dimension(months_in_year) :: &
       daysm = [31,28,31,30,31,30,31,31,30,31,30,31], &
       daysml = [31,29,31,30,31,30,31,31,30,31,30,31], &
       lastday = [31,59,90,120,151,181,212,243,273,304,334,365], &
       lastdayl = [31,60,91,121,152,182,213,244,274,305,335,366]

contains

  function time_step_matches(dels, ktau, frequency, leaps, start_year) result(match)
    real, intent(in) :: dels !! Model time step in seconds
    integer, intent(in) :: ktau !! Current time step index
    character(len=*), intent(in) :: frequency !! Frequency string: 'all', 'daily', 'monthly'
    logical, intent(in) :: leaps !! Are we using leap years?
    integer, intent(in) :: start_year !! Start year of the simulation
    logical :: match
    integer :: i
    integer :: time_steps_per_interval
    integer :: interval_in_hours
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
      last_day_of_month_in_total_elapsed_days = 0
      do i = start_year, current_year - 1
        if (leaps .and. is_leapyear(i)) then
          last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + 366
        else
          last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + 365
        end if
      end do
      if (leaps .and. is_leapyear(current_year)) then
        last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + lastdayl
      else
        last_day_of_month_in_total_elapsed_days = last_day_of_month_in_total_elapsed_days + lastday
      end if
      match = any(int(real(last_day_of_month_in_total_elapsed_days) * seconds_per_day / dels) == ktau)
    case default
      call cable_abort('Error: unknown frequency "' // trim(adjustl(frequency)) // '"', __FILE__, __LINE__)
    end select

  end function

end module
