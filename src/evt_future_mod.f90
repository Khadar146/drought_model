!> @file evt_future_mod.f90
!! @brief Extreme Value Theory module for future climate projections
!! 
!! This module applies EVT analysis to future climate scenarios to assess
!! changes in extreme drought characteristics under climate change.
!! Compares future extreme events with historical baseline.
!!
!! @author Khadar Daahir
!! @date August 2025
!!
module evt_future_mod
  use iso_fortran_env, only: wp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use fsml_dst, only: f_dst_gpd_cdf
  implicit none
  private
  
  public :: compute_evt_future, compare_evt_changes, project_future_extremes
  public :: assess_risk_changes, validate_future_evt
  
  !> Parameters for future EVT analysis
  real(wp), parameter :: CLIMATE_CHANGE_THRESHOLD = 1.5_wp  ! Significant change threshold
  integer,  parameter :: MIN_FUTURE_YEARS = 30              ! Minimum years for robust analysis
  
contains

  !> @brief Compute EVT analysis for future climate scenarios
  !! @param spi_future Future SPI data [nlon, nlat, ntime, nperiods]
  !! @param historical_evt Historical EVT parameters for comparison
  !! @param future_start_year Starting year of future period
  !! @param[out] future_evt Future EVT parameters
  !! @param[out] change_signals Detected changes in extreme characteristics
  subroutine compute_evt_future(spi_future, nlon, nlat, ntime, nperiods, &
                               historical_evt, future_start_year, &
                               future_evt, change_signals)
    integer,  intent(in)  :: nlon, nlat, ntime, nperiods, future_start_year
    real(wp), intent(in)  :: spi_future(nlon, nlat, ntime, nperiods)
    real(wp), intent(in)  :: historical_evt(nlon, nlat, nperiods, 3) ! threshold, shape, scale
    real(wp), intent(out) :: future_evt(nlon, nlat, nperiods, 3)
    logical,  intent(out) :: change_signals(nlon, nlat, nperiods, 4) ! threshold, shape, scale, frequency
    
    integer :: i, j, ap
    real(wp), allocatable :: extremes(:)
    real(wp) :: threshold, shape, scale
    logical :: fit_ok
    integer :: n_years
    
    ! Initialize outputs
    future_evt = ieee_value(0.0_wp, ieee_quiet_nan)
    change_signals = .false.
    
    n_years = ntime / 12  ! Assuming monthly data
    
    if (n_years < MIN_FUTURE_YEARS) then
      print *, "Warning: Future period too short for robust EVT analysis"
      return
    end if
    
    do ap = 1, nperiods
      do i = 1, nlon
        do j = 1, nlat
          ! Extract valid SPI values for this grid point and period
          call extract_valid_spi_future(spi_future(i,j,:,ap), ntime, extremes)
          
          if (size(extremes) > 50) then  ! Need sufficient data for future analysis
            ! Use historical threshold for consistency
            threshold = historical_evt(i,j,ap,1)
            
            if (.not. ieee_is_nan(threshold)) then
              future_evt(i,j,ap,1) = threshold
              
              ! Fit GPD to future exceedances
              call fit_gpd_to_exceedances_future(extremes, threshold, shape, scale, fit_ok)
              
              if (fit_ok) then
                future_evt(i,j,ap,2) = shape
                future_evt(i,j,ap,3) = scale
                
                ! Detect significant changes
                call detect_evt_changes(historical_evt(i,j,ap,:), future_evt(i,j,ap,:), &
                                       extremes, threshold, change_signals(i,j,ap,:))
              end if
            end if
          end if
          
          if (allocated(extremes)) deallocate(extremes)
        end do
      end do
    end do
    
    print *, "> Future EVT analysis completed for", nperiods, "SPI period(s)"
  end subroutine compute_evt_future

  !> @brief Compare EVT changes between historical and future periods
  subroutine compare_evt_changes(historical_evt, future_evt, &
                                nlon, nlat, nperiods, &
                                return_periods, nreturn, &
                                change_metrics)
    integer,  intent(in)  :: nlon, nlat, nperiods, nreturn
    real(wp), intent(in)  :: historical_evt(nlon, nlat, nperiods, 3)
    real(wp), intent(in)  :: future_evt(nlon, nlat, nperiods, 3)
    real(wp), intent(in)  :: return_periods(nreturn)
    real(wp), intent(out) :: change_metrics(nlon, nlat, nperiods, nreturn, 3) ! absolute, relative, significance
    
    integer :: i, j, ap, ir
    real(wp) :: hist_level, future_level, abs_change, rel_change
    
    change_metrics = ieee_value(0.0_wp, ieee_quiet_nan)
    
    do ap = 1, nperiods
      do i = 1, nlon
        do j = 1, nlat
          do ir = 1, nreturn
            ! Compute historical return level
            call compute_return_level_future(historical_evt(i,j,ap,1), &
                                           historical_evt(i,j,ap,2), &
                                           historical_evt(i,j,ap,3), &
                                           return_periods(ir), hist_level)
            
            ! Compute future return level
            call compute_return_level_future(future_evt(i,j,ap,1), &
                                           future_evt(i,j,ap,2), &
                                           future_evt(i,j,ap,3), &
                                           return_periods(ir), future_level)
            
            if (.not. ieee_is_nan(hist_level) .and. .not. ieee_is_nan(future_level)) then
              ! Absolute change
              abs_change = future_level - hist_level
              change_metrics(i,j,ap,ir,1) = abs_change
              
              ! Relative change (%)
              if (abs(hist_level) > 1.0e-6_wp) then
                rel_change = 100.0_wp * abs_change / abs(hist_level)
                change_metrics(i,j,ap,ir,2) = rel_change
              end if
              
              ! Significance assessment (simplified)
              if (abs(rel_change) > 20.0_wp) then  ! >20% change considered significant
                change_metrics(i,j,ap,ir,3) = 1.0_wp
              else
                change_metrics(i,j,ap,ir,3) = 0.0_wp
              end if
            end if
          end do
        end do
      end do
    end do
    
    print *, "> EVT change analysis completed"
  end subroutine compare_evt_changes

  !> @brief Project future extreme events and their probabilities
  subroutine project_future_extremes(future_evt, nlon, nlat, nperiods, &
                                    projection_years, extreme_thresholds, &
                                    projected_frequencies, confidence_intervals)
    integer,  intent(in)  :: nlon, nlat, nperiods
    real(wp), intent(in)  :: future_evt(nlon, nlat, nperiods, 3)
    real(wp), intent(in)  :: projection_years(:)
    real(wp), intent(in)  :: extreme_thresholds(:)
    real(wp), intent(out) :: projected_frequencies(nlon, nlat, nperiods, size(extreme_thresholds))
    real(wp), intent(out) :: confidence_intervals(nlon, nlat, nperiods, size(extreme_thresholds), 2)
    
    integer :: i, j, ap, ith
    real(wp) :: threshold, shape, scale, prob, frequency
    real(wp) :: lower_ci, upper_ci
    integer :: n_thresholds
    
    n_thresholds = size(extreme_thresholds)
    projected_frequencies = ieee_value(0.0_wp, ieee_quiet_nan)
    confidence_intervals = ieee_value(0.0_wp, ieee_quiet_nan)
    
    do ap = 1, nperiods
      do i = 1, nlon
        do j = 1, nlat
          threshold = future_evt(i,j,ap,1)
          shape = future_evt(i,j,ap,2)
          scale = future_evt(i,j,ap,3)
          
          if (.not. ieee_is_nan(threshold) .and. .not. ieee_is_nan(shape) .and. &
              .not. ieee_is_nan(scale) .and. scale > 0.0_wp) then
            
            do ith = 1, n_thresholds
              if (extreme_thresholds(ith) > threshold) then
                ! Compute exceedance probability
                call compute_exceedance_probability(extreme_thresholds(ith), &
                                                  threshold, shape, scale, prob)
                
                ! Convert to frequency (events per year)
                frequency = prob  ! Simplified: assuming annual exceedance rate = 1
                projected_frequencies(i,j,ap,ith) = frequency
                
                ! Rough confidence intervals (simplified)
                call compute_frequency_confidence(frequency, 30.0_wp, lower_ci, upper_ci)
                confidence_intervals(i,j,ap,ith,1) = lower_ci
                confidence_intervals(i,j,ap,ith,2) = upper_ci
              end if
            end do
          end if
        end do
      end do
    end do
    
    print *, "> Future extreme projections completed"
  end subroutine project_future_extremes

  !> @brief Assess changes in drought risk under future climate
  subroutine assess_risk_changes(historical_frequencies, future_frequencies, &
                                nlon, nlat, nperiods, n_thresholds, &
                                risk_ratios, risk_categories)
    integer,  intent(in)  :: nlon, nlat, nperiods, n_thresholds
    real(wp), intent(in)  :: historical_frequencies(nlon, nlat, nperiods, n_thresholds)
    real(wp), intent(in)  :: future_frequencies(nlon, nlat, nperiods, n_thresholds)
    real(wp), intent(out) :: risk_ratios(nlon, nlat, nperiods, n_thresholds)
    integer,  intent(out) :: risk_categories(nlon, nlat, nperiods, n_thresholds)
    
    integer :: i, j, ap, ith
    real(wp) :: ratio
    
    risk_ratios = ieee_value(0.0_wp, ieee_quiet_nan)
    risk_categories = 0  ! 0=no change, 1=low increase, 2=moderate increase, 3=high increase, -1=decrease
    
    do ap = 1, nperiods
      do i = 1, nlon
        do j = 1, nlat
          do ith = 1, n_thresholds
            if (.not. ieee_is_nan(historical_frequencies(i,j,ap,ith)) .and. &
                .not. ieee_is_nan(future_frequencies(i,j,ap,ith)) .and. &
                historical_frequencies(i,j,ap,ith) > 1.0e-6_wp) then
              
              ratio = future_frequencies(i,j,ap,ith) / historical_frequencies(i,j,ap,ith)
              risk_ratios(i,j,ap,ith) = ratio
              
              ! Categorize risk change
              if (ratio > 2.0_wp) then
                risk_categories(i,j,ap,ith) = 3  ! High increase (>100%)
              else if (ratio > 1.5_wp) then
                risk_categories(i,j,ap,ith) = 2  ! Moderate increase (50-100%)
              else if (ratio > 1.2_wp) then
                risk_categories(i,j,ap,ith) = 1  ! Low increase (20-50%)
              else if (ratio < 0.8_wp) then
                risk_categories(i,j,ap,ith) = -1  ! Decrease (>20% reduction)
              else
                risk_categories(i,j,ap,ith) = 0   ! No significant change
              end if
            end if
          end do
        end do
      end do
    end do
    
    print *, "> Risk assessment completed"
  end subroutine assess_risk_changes

  !> @brief Validate future EVT results
  logical function validate_future_evt(future_evt, spi_future, nlon, nlat, ntime, nperiods) result(is_valid)
    integer,  intent(in) :: nlon, nlat, ntime, nperiods
    real(wp), intent(in) :: future_evt(nlon, nlat, nperiods, 3)
    real(wp), intent(in) :: spi_future(nlon, nlat, ntime, nperiods)
    
    integer :: i, j, ap, valid_fits, total_points
    real(wp) :: shape, scale
    
    is_valid = .true.
    valid_fits = 0
    total_points = 0
    
    do ap = 1, nperiods
      do i = 1, nlon
        do j = 1, nlat
          total_points = total_points + 1
          shape = future_evt(i,j,ap,2)
          scale = future_evt(i,j,ap,3)
          
          if (.not. ieee_is_nan(shape) .and. .not. ieee_is_nan(scale)) then
            ! Check parameter validity
            if (scale > 0.0_wp .and. abs(shape) < 0.5_wp) then
              valid_fits = valid_fits + 1
            else
              is_valid = .false.
            end if
          end if
        end do
      end do
    end do
    
    ! Require at least 50% successful fits
    if (real(valid_fits, wp) / real(total_points, wp) < 0.5_wp) then
      is_valid = .false.
    end if
    
    print *, "> Future EVT validation:", merge("PASSED", "FAILED", is_valid)
    print *, "> Valid fits:", valid_fits, "/", total_points
  end function validate_future_evt

  ! Helper subroutines
  
  subroutine extract_valid_spi_future(spi_series, ntime, valid_data)
    integer,  intent(in)  :: ntime
    real(wp), intent(in)  :: spi_series(ntime)
    real(wp), allocatable, intent(out) :: valid_data(:)
    
    integer :: t, count_valid
    
    count_valid = 0
    do t = 1, ntime
      if (.not. ieee_is_nan(spi_series(t))) count_valid = count_valid + 1
    end do
    
    if (count_valid > 0) then
      allocate(valid_data(count_valid))
      count_valid = 0
      do t = 1, ntime
        if (.not. ieee_is_nan(spi_series(t))) then
          count_valid = count_valid + 1
          valid_data(count_valid) = spi_series(t)
        end if
      end do
    else
      allocate(valid_data(0))
    end if
  end subroutine extract_valid_spi_future

  subroutine fit_gpd_to_exceedances_future(data, threshold, shape, scale, success)
    real(wp), intent(in)  :: data(:)
    real(wp), intent(in)  :: threshold
    real(wp), intent(out) :: shape, scale
    logical,  intent(out) :: success
    
    real(wp), allocatable :: exceedances(:)
    integer :: i, n_exceed
    
    success = .false.
    shape = ieee_value(0.0_wp, ieee_quiet_nan)
    scale = ieee_value(0.0_wp, ieee_quiet_nan)
    
    n_exceed = count(data > threshold)
    if (n_exceed < 20) return  ! Need sufficient exceedances
    
    allocate(exceedances(n_exceed))
    n_exceed = 0
    do i = 1, size(data)
      if (data(i) > threshold) then
        n_exceed = n_exceed + 1
        exceedances(n_exceed) = data(i) - threshold
      end if
    end do
    
    ! Use method of moments for fitting
    call fit_gpd_moments_future(exceedances, shape, scale, success)
    
    deallocate(exceedances)
  end subroutine fit_gpd_to_exceedances_future

  subroutine fit_gpd_moments_future(exceedances, shape, scale, success)
    real(wp), intent(in)  :: exceedances(:)
    real(wp), intent(out) :: shape, scale
    logical,  intent(out) :: success
    
    real(wp) :: mean_exc, var_exc
    integer :: n
    
    n = size(exceedances)
    success = .false.
    
    if (n < 2) return
    
    mean_exc = sum(exceedances) / real(n, wp)
    var_exc = sum((exceedances - mean_exc)**2) / real(n - 1, wp)
    
    if (mean_exc > 0.0_wp .and. var_exc > 0.0_wp) then
      shape = 0.5_wp * ((mean_exc**2) / var_exc - 1.0_wp)
      scale = 0.5_wp * mean_exc * ((mean_exc**2) / var_exc + 1.0_wp)
      
      if (scale > 0.0_wp .and. abs(shape) < 0.5_wp) then
        success = .true.
      end if
    end if
  end subroutine fit_gpd_moments_future

  subroutine detect_evt_changes(historical_params, future_params, extremes, threshold, changes)
    real(wp), intent(in)  :: historical_params(3), future_params(3)
    real(wp), intent(in)  :: extremes(:)
    real(wp), intent(in)  :: threshold
    logical,  intent(out) :: changes(4)
    
    real(wp) :: shape_change, scale_change, freq_change
    integer :: n_exceed_hist, n_exceed_fut
    
    changes = .false.
    
    ! Shape parameter change
    if (.not. ieee_is_nan(historical_params(2)) .and. .not. ieee_is_nan(future_params(2))) then
      shape_change = abs(future_params(2) - historical_params(2))
      if (shape_change > 0.1_wp) changes(2) = .true.
    end if
    
    ! Scale parameter change
    if (.not. ieee_is_nan(historical_params(3)) .and. .not. ieee_is_nan(future_params(3))) then
      if (historical_params(3) > 0.0_wp) then
        scale_change = abs((future_params(3) - historical_params(3)) / historical_params(3))
        if (scale_change > 0.2_wp) changes(3) = .true.
      end if
    end if
    
    ! Frequency change (simplified)
    n_exceed_fut = count(extremes > threshold)
    if (n_exceed_fut > 0) then
      changes(4) = .true.  ! Simplified: any exceedances indicate potential change
    end if
  end subroutine detect_evt_changes

  subroutine compute_return_level_future(threshold, shape, scale, return_period, return_level)
    real(wp), intent(in)  :: threshold, shape, scale, return_period
    real(wp), intent(out) :: return_level
    
    real(wp) :: lambda, prob
    
    return_level = ieee_value(0.0_wp, ieee_quiet_nan)
    
    if (scale <= 0.0_wp .or. return_period <= 0.0_wp) return
    
    lambda = 1.0_wp  ! Annual exceedance rate
    prob = 1.0_wp - 1.0_wp / (lambda * return_period)
    
    if (abs(shape) < 1.0e-6_wp) then
      return_level = threshold - scale * log(1.0_wp - prob)
    else
      return_level = threshold + (scale / shape) * ((1.0_wp - prob)**(-shape) - 1.0_wp)
    end if
  end subroutine compute_return_level_future

  subroutine compute_exceedance_probability(value, threshold, shape, scale, prob)
    real(wp), intent(in)  :: value, threshold, shape, scale
    real(wp), intent(out) :: prob
    
    if (value <= threshold .or. scale <= 0.0_wp) then
      prob = 0.0_wp
      return
    end if
    
    if (abs(shape) < 1.0e-6_wp) then
      prob = exp(-(value - threshold) / scale)
    else
      prob = (1.0_wp + shape * (value - threshold) / scale)**(-1.0_wp/shape)
    end if
    
    prob = max(0.0_wp, min(1.0_wp, prob))
  end subroutine compute_exceedance_probability

  subroutine compute_frequency_confidence(frequency, n_years, lower_ci, upper_ci)
    real(wp), intent(in)  :: frequency, n_years
    real(wp), intent(out) :: lower_ci, upper_ci
    
    real(wp) :: margin
    
    ! Simplified confidence interval (assuming normal approximation)
    margin = 1.96_wp * sqrt(frequency / n_years)  ! 95% CI
    lower_ci = max(0.0_wp, frequency - margin)
    upper_ci = frequency + margin
  end subroutine compute_frequency_confidence

end module evt_future_mod
