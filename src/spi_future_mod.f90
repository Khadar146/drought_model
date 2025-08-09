!> @file spi_future_mod.f90
!! @brief SPI module for future climate projections
!! 
!! This module computes SPI for future climate scenarios using historical reference periods.
!! Designed to detect changes in drought characteristics under climate change.
!! Uses the same methodology as historical SPI but with different data sources.
!!
!! @author Khadar Daahir
!! @date August 2025
!!
module spi_future_mod
  use iso_fortran_env, only: wp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use fsml_dst, only: f_dst_gammai_core, f_dst_norm_ppf_core
  implicit none
  private
  
  public :: compute_spi_future, apply_bias_correction, validate_future_spi
  public :: compare_with_historical
  
contains

  !> @brief Compute SPI for future climate projections
  !! @param future_precip Future precipitation data [nlon, nlat, ntime]
  !! @param historical_params Historical gamma parameters for reference [nlon, nlat, 12, nperiods, 3]
  !! @param start_year Starting year of future data
  !! @param start_month Starting month of future data
  !! @param accum_periods Array of accumulation periods
  !! @param[out] spi_future Future SPI values [nlon, nlat, ntime, nperiods]
  !! @param[out] anomaly_flags Flags for extreme anomalies [nlon, nlat, ntime, nperiods]
  subroutine compute_spi_future(future_precip, nlon, nlat, ntime, &
                               historical_params, &
                               start_year, start_month, &
                               accum_periods, num_periods, &
                               spi_future, anomaly_flags, &
                               year_idx, month_idx)
    integer,  intent(in) :: nlon, nlat, ntime, num_periods
    real(wp), intent(in) :: future_precip(nlon, nlat, ntime)
    real(wp), intent(in) :: historical_params(nlon, nlat, 12, num_periods, 3) ! shape, scale, p_zero
    integer,  intent(in) :: start_year, start_month
    integer,  intent(in) :: accum_periods(num_periods)
    real(wp), intent(out):: spi_future(nlon, nlat, ntime, num_periods)
    logical,  intent(out):: anomaly_flags(nlon, nlat, ntime, num_periods)
    integer,  intent(in), optional :: year_idx(ntime), month_idx(ntime)

    integer :: ap, i, j, k
    real(wp), allocatable :: acc_series(:,:,:)
    logical :: have_calendar
    real(wp), parameter :: EXTREME_THRESHOLD = 3.0_wp

    spi_future = ieee_value(0.0_wp, ieee_quiet_nan)
    anomaly_flags = .false.
    have_calendar = present(year_idx) .and. present(month_idx)

    do ap = 1, num_periods
      k = accum_periods(ap)
      call accumulate_field_future(future_precip, nlon, nlat, ntime, k, acc_series)

      do i = 1, nlon
        do j = 1, nlat
          if (have_calendar) then
            call spi_for_cell_future_arrays(acc_series(i,j,:), ntime, &
                                           historical_params(i,j,:,ap,:), &
                                           month_idx, k, spi_future(i,j,:,ap))
          else
            call spi_for_cell_future_linear(acc_series(i,j,:), ntime, &
                                           historical_params(i,j,:,ap,:), &
                                           start_month, k, spi_future(i,j,:,ap))
          end if
          
          ! Flag extreme anomalies
          call flag_extreme_anomalies(spi_future(i,j,:,ap), ntime, EXTREME_THRESHOLD, &
                                    anomaly_flags(i,j,:,ap))
        end do
      end do

      deallocate(acc_series)
    end do

    print *, "> Future SPI computed for", num_periods, "accumulation period(s)."
  end subroutine compute_spi_future

  !> @brief Apply bias correction to future precipitation data
  subroutine apply_bias_correction(raw_future, historical_obs, historical_model, &
                                 nlon, nlat, ntime_fut, ntime_hist, &
                                 corrected_future)
    integer,  intent(in)  :: nlon, nlat, ntime_fut, ntime_hist
    real(wp), intent(in)  :: raw_future(nlon, nlat, ntime_fut)
    real(wp), intent(in)  :: historical_obs(nlon, nlat, ntime_hist)
    real(wp), intent(in)  :: historical_model(nlon, nlat, ntime_hist)
    real(wp), intent(out) :: corrected_future(nlon, nlat, ntime_fut)
    
    integer :: i, j, t, m
    real(wp) :: obs_mean(12), model_mean(12), bias_factor(12)
    real(wp) :: obs_sum, model_sum
    integer :: obs_count, model_count
    
    corrected_future = raw_future
    
    do i = 1, nlon
      do j = 1, nlat
        ! Compute monthly bias correction factors
        do m = 1, 12
          obs_sum = 0.0_wp; model_sum = 0.0_wp
          obs_count = 0; model_count = 0
          
          ! Calculate monthly means from historical data
          do t = 1, ntime_hist
            if (mod(t-1, 12) + 1 == m) then
              if (.not. ieee_is_nan(historical_obs(i,j,t))) then
                obs_sum = obs_sum + historical_obs(i,j,t)
                obs_count = obs_count + 1
              end if
              if (.not. ieee_is_nan(historical_model(i,j,t))) then
                model_sum = model_sum + historical_model(i,j,t)
                model_count = model_count + 1
              end if
            end if
          end do
          
          if (obs_count > 0 .and. model_count > 0) then
            obs_mean(m) = obs_sum / real(obs_count, wp)
            model_mean(m) = model_sum / real(model_count, wp)
            
            if (model_mean(m) > 0.0_wp) then
              bias_factor(m) = obs_mean(m) / model_mean(m)
            else
              bias_factor(m) = 1.0_wp
            end if
          else
            bias_factor(m) = 1.0_wp
          end if
        end do
        
        ! Apply correction to future data
        do t = 1, ntime_fut
          m = mod(t-1, 12) + 1
          if (.not. ieee_is_nan(raw_future(i,j,t))) then
            corrected_future(i,j,t) = raw_future(i,j,t) * bias_factor(m)
          end if
        end do
      end do
    end do
    
    print *, "> Bias correction applied to future precipitation data"
  end subroutine apply_bias_correction

  !> @brief Compare future SPI with historical baseline
  subroutine compare_with_historical(spi_future, spi_historical, &
                                   nlon, nlat, ntime_fut, ntime_hist, nperiods, &
                                   change_metrics)
    integer,  intent(in)  :: nlon, nlat, ntime_fut, ntime_hist, nperiods
    real(wp), intent(in)  :: spi_future(nlon, nlat, ntime_fut, nperiods)
    real(wp), intent(in)  :: spi_historical(nlon, nlat, ntime_hist, nperiods)
    real(wp), intent(out) :: change_metrics(nlon, nlat, nperiods, 5) ! mean, std, extreme_freq, drought_freq, trend
    
    integer :: i, j, ap, t
    real(wp) :: future_mean, hist_mean, future_std, hist_std
    real(wp) :: future_extreme_freq, hist_extreme_freq
    real(wp) :: future_drought_freq, hist_drought_freq
    integer :: valid_count_fut, valid_count_hist
    integer :: extreme_count_fut, extreme_count_hist
    integer :: drought_count_fut, drought_count_hist
    real(wp), parameter :: EXTREME_THRESHOLD = 2.0_wp
    real(wp), parameter :: DROUGHT_THRESHOLD = -1.0_wp
    
    change_metrics = ieee_value(0.0_wp, ieee_quiet_nan)
    
    do ap = 1, nperiods
      do i = 1, nlon
        do j = 1, nlat
          ! Calculate statistics for future period
          call calculate_spi_statistics(spi_future(i,j,:,ap), ntime_fut, &
                                      future_mean, future_std, valid_count_fut)
          
          ! Calculate statistics for historical period
          call calculate_spi_statistics(spi_historical(i,j,:,ap), ntime_hist, &
                                      hist_mean, hist_std, valid_count_hist)
          
          if (valid_count_fut > 0 .and. valid_count_hist > 0) then
            ! Mean change
            change_metrics(i,j,ap,1) = future_mean - hist_mean
            
            ! Standard deviation change
            change_metrics(i,j,ap,2) = future_std - hist_std
            
            ! Extreme frequency change
            extreme_count_fut = count_extremes(spi_future(i,j,:,ap), ntime_fut, EXTREME_THRESHOLD)
            extreme_count_hist = count_extremes(spi_historical(i,j,:,ap), ntime_hist, EXTREME_THRESHOLD)
            
            future_extreme_freq = real(extreme_count_fut, wp) / real(valid_count_fut, wp)
            hist_extreme_freq = real(extreme_count_hist, wp) / real(valid_count_hist, wp)
            change_metrics(i,j,ap,3) = future_extreme_freq - hist_extreme_freq
            
            ! Drought frequency change
            drought_count_fut = count_droughts(spi_future(i,j,:,ap), ntime_fut, DROUGHT_THRESHOLD)
            drought_count_hist = count_droughts(spi_historical(i,j,:,ap), ntime_hist, DROUGHT_THRESHOLD)
            
            future_drought_freq = real(drought_count_fut, wp) / real(valid_count_fut, wp)
            hist_drought_freq = real(drought_count_hist, wp) / real(valid_count_hist, wp)
            change_metrics(i,j,ap,4) = future_drought_freq - hist_drought_freq
            
            ! Trend analysis (simplified)
            call calculate_trend(spi_future(i,j,:,ap), ntime_fut, change_metrics(i,j,ap,5))
          end if
        end do
      end do
    end do
    
    print *, "> Comparison with historical baseline completed"
  end subroutine compare_with_historical

  !> @brief Validate future SPI results
  logical function validate_future_spi(spi_future, nlon, nlat, ntime, nperiods) result(is_valid)
    integer,  intent(in) :: nlon, nlat, ntime, nperiods
    real(wp), intent(in) :: spi_future(nlon, nlat, ntime, nperiods)
    
    integer :: i, j, ap, valid_count, total_count
    real(wp) :: mean_val, std_val
    
    is_valid = .true.
    total_count = 0
    valid_count = 0
    
    do ap = 1, nperiods
      do i = 1, nlon
        do j = 1, nlat
          call calculate_spi_statistics(spi_future(i,j,:,ap), ntime, mean_val, std_val, valid_count)
          total_count = total_count + 1
          
          ! Check for reasonable statistics
          if (valid_count > 0) then
            if (abs(mean_val) > 2.0_wp .or. std_val < 0.5_wp .or. std_val > 3.0_wp) then
              is_valid = .false.
              return
            end if
          end if
        end do
      end do
    end do
    
    print *, "> Future SPI validation:", merge("PASSED", "FAILED", is_valid)
  end function validate_future_spi

  ! Helper subroutines (similar to historical SPI module)
  
  subroutine accumulate_field_future(x, nlon, nlat, ntime, k, acc)
    real(wp), intent(in)  :: x(nlon, nlat, ntime)
    integer,  intent(in)  :: nlon, nlat, ntime, k
    real(wp), allocatable, intent(out) :: acc(:,:,:)
    integer :: t, i, j

    allocate(acc(nlon, nlat, ntime))
    acc = ieee_value(0.0_wp, ieee_quiet_nan)

    do t = k, ntime
      do i = 1, nlon
        do j = 1, nlat
          if (all(.not. ieee_is_nan(x(i,j,t-k+1:t)))) then
            acc(i,j,t) = sum(x(i,j,t-k+1:t))
          end if
        end do
      end do
    end do
  end subroutine accumulate_field_future

  subroutine spi_for_cell_future_arrays(acc, ntime, historical_params, month_idx, k, spi_cell)
    real(wp), intent(in)  :: acc(ntime)
    integer,  intent(in)  :: ntime, k
    real(wp), intent(in)  :: historical_params(12, 3) ! [month, param] where param: 1=shape, 2=scale, 3=p_zero
    integer,  intent(in)  :: month_idx(ntime)
    real(wp), intent(out) :: spi_cell(ntime)

    integer :: t, m

    spi_cell = ieee_value(0.0_wp, ieee_quiet_nan)

    do t = k, ntime
      m = month_idx(t)
      if (.not. ieee_is_nan(acc(t)) .and. historical_params(m,2) > 0.0_wp) then
        spi_cell(t) = gamma_to_spi_future(acc(t), historical_params(m,1), &
                                        historical_params(m,2), historical_params(m,3))
      end if
    end do
  end subroutine spi_for_cell_future_arrays

  subroutine spi_for_cell_future_linear(acc, ntime, historical_params, start_month, k, spi_cell)
    real(wp), intent(in)  :: acc(ntime)
    integer,  intent(in)  :: ntime, start_month, k
    real(wp), intent(in)  :: historical_params(12, 3)
    real(wp), intent(out) :: spi_cell(ntime)

    integer :: t, m

    spi_cell = ieee_value(0.0_wp, ieee_quiet_nan)

    do t = k, ntime
      m = month_of_index_future(t, start_month)
      if (.not. ieee_is_nan(acc(t)) .and. historical_params(m,2) > 0.0_wp) then
        spi_cell(t) = gamma_to_spi_future(acc(t), historical_params(m,1), &
                                        historical_params(m,2), historical_params(m,3))
      end if
    end do
  end subroutine spi_for_cell_future_linear

  pure function gamma_to_spi_future(value, shape, scale, p_zero) result(spi_val)
    real(wp), intent(in) :: value, shape, scale, p_zero
    real(wp) :: spi_val, Pg, P
    real(wp), parameter :: eps = 1.0e-12_wp

    if (shape <= 0.0_wp .or. scale <= 0.0_wp) then
      spi_val = ieee_value(0.0_wp, ieee_quiet_nan)
      return
    end if

    if (value <= 0.0_wp) then
      Pg = 0.0_wp
    else
      Pg = f_dst_gammai_core(value/scale, shape) / gamma(shape)
    end if

    P = p_zero + (1.0_wp - p_zero) * Pg
    P = max(eps, min(1.0_wp - eps, P))

    spi_val = f_dst_norm_ppf_core(P, 0.0_wp, 1.0_wp)
  end function gamma_to_spi_future

  pure function month_of_index_future(t, start_month) result(m)
    integer, intent(in) :: t, start_month
    integer :: m, z
    z = (start_month - 1) + (t - 1)
    m = 1 + mod(z, 12)
  end function month_of_index_future

  subroutine flag_extreme_anomalies(spi_series, ntime, threshold, flags)
    real(wp), intent(in)  :: spi_series(ntime)
    integer,  intent(in)  :: ntime
    real(wp), intent(in)  :: threshold
    logical,  intent(out) :: flags(ntime)
    
    integer :: t
    
    do t = 1, ntime
      flags(t) = .false.
      if (.not. ieee_is_nan(spi_series(t))) then
        flags(t) = abs(spi_series(t)) > threshold
      end if
    end do
  end subroutine flag_extreme_anomalies

  subroutine calculate_spi_statistics(spi_series, ntime, mean_val, std_val, valid_count)
    real(wp), intent(in)  :: spi_series(ntime)
    integer,  intent(in)  :: ntime
    real(wp), intent(out) :: mean_val, std_val
    integer,  intent(out) :: valid_count
    
    integer :: t
    real(wp) :: sum_val, sum_sq
    
    valid_count = 0
    sum_val = 0.0_wp
    sum_sq = 0.0_wp
    
    do t = 1, ntime
      if (.not. ieee_is_nan(spi_series(t))) then
        valid_count = valid_count + 1
        sum_val = sum_val + spi_series(t)
        sum_sq = sum_sq + spi_series(t)**2
      end if
    end do
    
    if (valid_count > 0) then
      mean_val = sum_val / real(valid_count, wp)
      if (valid_count > 1) then
        std_val = sqrt((sum_sq - sum_val**2/real(valid_count, wp)) / real(valid_count - 1, wp))
      else
        std_val = 0.0_wp
      end if
    else
      mean_val = ieee_value(0.0_wp, ieee_quiet_nan)
      std_val = ieee_value(0.0_wp, ieee_quiet_nan)
    end if
  end subroutine calculate_spi_statistics

  integer function count_extremes(spi_series, ntime, threshold) result(count)
    real(wp), intent(in) :: spi_series(ntime)
    integer,  intent(in) :: ntime
    real(wp), intent(in) :: threshold
    
    integer :: t
    
    count = 0
    do t = 1, ntime
      if (.not. ieee_is_nan(spi_series(t))) then
        if (abs(spi_series(t)) > threshold) count = count + 1
      end if
    end do
  end function count_extremes

  integer function count_droughts(spi_series, ntime, threshold) result(count)
    real(wp), intent(in) :: spi_series(ntime)
    integer,  intent(in) :: ntime
    real(wp), intent(in) :: threshold
    
    integer :: t
    
    count = 0
    do t = 1, ntime
      if (.not. ieee_is_nan(spi_series(t))) then
        if (spi_series(t) < threshold) count = count + 1
      end if
    end do
  end function count_droughts

  subroutine calculate_trend(spi_series, ntime, trend_slope)
    real(wp), intent(in)  :: spi_series(ntime)
    integer,  intent(in)  :: ntime
    real(wp), intent(out) :: trend_slope
    
    ! Simplified trend calculation using least squares
    real(wp) :: sum_x, sum_y, sum_xy, sum_x2, denom
    integer :: t, valid_count
    
    sum_x = 0.0_wp; sum_y = 0.0_wp; sum_xy = 0.0_wp; sum_x2 = 0.0_wp
    valid_count = 0
    
    do t = 1, ntime
      if (.not. ieee_is_nan(spi_series(t))) then
        valid_count = valid_count + 1
        sum_x = sum_x + real(t, wp)
        sum_y = sum_y + spi_series(t)
        sum_xy = sum_xy + real(t, wp) * spi_series(t)
        sum_x2 = sum_x2 + real(t, wp)**2
      end if
    end do
    
    if (valid_count > 1) then
      denom = real(valid_count, wp) * sum_x2 - sum_x**2
      if (abs(denom) > 1.0e-12_wp) then
        trend_slope = (real(valid_count, wp) * sum_xy - sum_x * sum_y) / denom
      else
        trend_slope = 0.0_wp
      end if
    else
      trend_slope = ieee_value(0.0_wp, ieee_quiet_nan)
    end if
  end subroutine calculate_trend

end module spi_future_mod
