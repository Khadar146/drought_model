! =====================================================================
! File: src/drought_analysis_mod.f90
! Purpose: Drought analysis tools for Somaliland climate research
! Author: Khadar - Advanced drought characterization and trends
! =====================================================================
module drought_analysis_mod
  use iso_fortran_env, only: wp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use config_mod, only: drought_config
  implicit none
  private
  public :: analyze_drought_characteristics, compute_drought_statistics
  public :: drought_event, drought_statistics

  type :: drought_event
    integer :: start_time    ! Time index when drought started
    integer :: end_time      ! Time index when drought ended  
    integer :: duration      ! Duration in months
    real(wp) :: severity     ! Average SPI during event
    real(wp) :: intensity    ! Minimum SPI reached
    real(wp) :: magnitude    ! Sum of negative SPI values
  end type drought_event

  type :: drought_statistics
    integer :: total_events
    integer :: extreme_events    ! Count of extreme droughts (SPI <= -2.0)
    real(wp) :: avg_duration     ! Average drought duration (months)
    real(wp) :: max_duration     ! Longest drought (months)
    real(wp) :: avg_severity     ! Average drought severity
    real(wp) :: trend_slope      ! Linear trend in SPI (per decade)
    real(wp) :: trend_p_value    ! P-value for trend significance
    real(wp) :: frequency        ! Droughts per decade
  end type drought_statistics

contains

  subroutine analyze_drought_characteristics(spi_ts, year_idx, month_idx, cfg, &
                                           events, stats, region_name)
    real(wp), intent(in) :: spi_ts(:)              ! SPI time series
    integer, intent(in) :: year_idx(:), month_idx(:)
    type(drought_config), intent(in) :: cfg
    type(drought_event), allocatable, intent(out) :: events(:)
    type(drought_statistics), intent(out) :: stats
    character(len=*), intent(in), optional :: region_name
    
    logical, allocatable :: is_drought(:)
    integer :: t, n_events, event_start, event_end
    integer :: ntime
    character(len=50) :: region
    
    ntime = size(spi_ts)
    region = "Grid Point"
    if (present(region_name)) region = region_name
    
    ! Identify drought periods (SPI <= -1.0 for moderate drought)
    allocate(is_drought(ntime))
    is_drought = (spi_ts <= cfg%moderate_drought) .and. (.not. ieee_is_nan(spi_ts))
    
    ! Count and extract drought events
    call extract_drought_events(spi_ts, is_drought, events)
    
    ! Compute statistics
    call compute_drought_statistics(spi_ts, events, year_idx, stats)
    
    ! Print summary for this location/region
    if (allocated(events) .and. size(events) > 0) then
      print *, "=== DROUGHT ANALYSIS FOR ", trim(region), " ==="
      print *, "Total drought events: ", stats%total_events
      print *, "Extreme drought events (SPI <= -2.0): ", stats%extreme_events
      print *, "Average drought duration: ", stats%avg_duration, " months"
      print *, "Longest drought: ", stats%max_duration, " months"
      print *, "Drought frequency: ", stats%frequency, " events per decade"
      print *, "Average severity: ", stats%avg_severity
      if (abs(stats%trend_p_value) < 0.05) then
        print *, "Significant trend detected (p <", stats%trend_p_value, ")"
        print *, "SPI trend: ", stats%trend_slope, " per decade"
      end if
      print *, "================================================="
    end if
  end subroutine analyze_drought_characteristics

  subroutine extract_drought_events(spi_ts, is_drought, events)
    real(wp), intent(in) :: spi_ts(:)
    logical, intent(in) :: is_drought(:)
    type(drought_event), allocatable, intent(out) :: events(:)
    
    integer :: t, ntime, n_events, event_idx
    integer :: event_start, event_end
    logical :: in_drought
    type(drought_event), allocatable :: temp_events(:)
    
    ntime = size(spi_ts)
    
    ! First pass: count events
    n_events = 0
    in_drought = .false.
    
    do t = 1, ntime
      if (is_drought(t) .and. .not. in_drought) then
        ! Start of new drought event
        n_events = n_events + 1
        in_drought = .true.
      else if (.not. is_drought(t) .and. in_drought) then
        ! End of drought event
        in_drought = .false.
      end if
    end do
    
    if (n_events == 0) then
      allocate(events(0))
      return
    end if
    
    ! Second pass: extract events
    allocate(temp_events(n_events))
    event_idx = 0
    in_drought = .false.
    
    do t = 1, ntime
      if (is_drought(t) .and. .not. in_drought) then
        ! Start of new drought event
        event_idx = event_idx + 1
        event_start = t
        in_drought = .true.
      else if ((.not. is_drought(t) .or. t == ntime) .and. in_drought) then
        ! End of drought event
        event_end = t - 1
        if (t == ntime .and. is_drought(t)) event_end = t
        
        ! Characterize the event
        call characterize_drought_event(spi_ts, event_start, event_end, &
                                       temp_events(event_idx))
        in_drought = .false.
      end if
    end do
    
    ! Copy to output
    events = temp_events
  end subroutine extract_drought_events

  subroutine characterize_drought_event(spi_ts, start_t, end_t, event)
    real(wp), intent(in) :: spi_ts(:)
    integer, intent(in) :: start_t, end_t
    type(drought_event), intent(out) :: event
    
    integer :: t, duration
    real(wp) :: total_deficit
    
    duration = end_t - start_t + 1
    
    event%start_time = start_t
    event%end_time = end_t
    event%duration = duration
    event%intensity = minval(spi_ts(start_t:end_t))
    event%severity = sum(spi_ts(start_t:end_t)) / real(duration, wp)
    
    ! Magnitude: sum of negative departures from zero
    total_deficit = 0.0_wp
    do t = start_t, end_t
      if (spi_ts(t) < 0.0_wp) then
        total_deficit = total_deficit + abs(spi_ts(t))
      end if
    end do
    event%magnitude = total_deficit
  end subroutine characterize_drought_event

  subroutine compute_drought_statistics(spi_ts, events, year_idx, stats)
    real(wp), intent(in) :: spi_ts(:)
    type(drought_event), intent(in) :: events(:)
    integer, intent(in) :: year_idx(:)
    type(drought_statistics), intent(out) :: stats
    
    integer :: i, n_years
    real(wp) :: total_duration, total_severity
    
    stats%total_events = size(events)
    
    if (stats%total_events == 0) then
      stats%extreme_events = 0
      stats%avg_duration = 0.0_wp
      stats%max_duration = 0.0_wp
      stats%avg_severity = 0.0_wp
      stats%frequency = 0.0_wp
      stats%trend_slope = 0.0_wp
      stats%trend_p_value = 1.0_wp
      return
    end if
    
    ! Count extreme events
    stats%extreme_events = count(events%intensity <= -2.0_wp)
    
    ! Duration statistics
    total_duration = sum(real(events%duration, wp))
    stats%avg_duration = total_duration / real(stats%total_events, wp)
    stats%max_duration = real(maxval(events%duration), wp)
    
    ! Severity statistics
    total_severity = sum(events%severity)
    stats%avg_severity = total_severity / real(stats%total_events, wp)
    
    ! Frequency (events per decade)
    n_years = year_idx(size(year_idx)) - year_idx(1) + 1
    stats%frequency = real(stats%total_events, wp) * 10.0_wp / real(n_years, wp)
    
    ! Trend analysis (simplified linear trend)
    call compute_linear_trend(spi_ts, stats%trend_slope, stats%trend_p_value)
  end subroutine compute_drought_statistics

  subroutine compute_linear_trend(data, slope, p_value)
    real(wp), intent(in) :: data(:)
    real(wp), intent(out) :: slope, p_value
    
    integer :: n, i
    real(wp) :: x_mean, y_mean, num, den
    real(wp), allocatable :: x(:), y(:)
    
    n = size(data)
    allocate(x(n), y(n))
    
    ! Remove NaN values
    y = data
    x = [(real(i, wp), i=1,n)]
    
    ! Simple linear regression
    x_mean = sum(x) / real(n, wp)
    y_mean = sum(y, mask=(.not. ieee_is_nan(y))) / real(count(.not. ieee_is_nan(y)), wp)
    
    num = sum((x - x_mean) * (y - y_mean), mask=(.not. ieee_is_nan(y)))
    den = sum((x - x_mean)**2)
    
    if (den > 0.0_wp) then
      slope = num / den
      ! Convert to per-decade trend
      slope = slope * 120.0_wp  ! 120 months = 10 years
    else
      slope = 0.0_wp
    end if
    
    ! Simplified p-value (would need proper statistical test for accuracy)
    p_value = 0.5_wp  ! Placeholder - implement proper significance test if needed
  end subroutine compute_linear_trend

end module drought_analysis_mod
