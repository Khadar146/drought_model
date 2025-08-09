!> @file evt_mod.f90
!! @brief Extreme Value Theory module for historical drought analysis
!! 
!! This module implements extreme value theory analysis for historical drought events
!! using Generalized Pareto Distribution (GPD) and Peak-Over-Threshold methods.
!! Designed to establish baseline extreme event statistics for future comparison.
!!
!! @author Khadar Daahir
!! @date August 2025
!!
module evt_mod
  use iso_fortran_env, only: wp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use fsml_dst, only: f_dst_gpd_cdf
  implicit none
  private
  
  public :: compute_evt_analysis, fit_gpd_threshold, extract_extremes
  public :: compute_return_periods, validate_gpd_fit
  
  !> Parameters for EVT analysis
  real(wp), parameter :: DEFAULT_THRESHOLD_QUANTILE = 0.95_wp
  integer,  parameter :: MIN_EXCEEDANCES = 20
  
contains

  !> @brief Main EVT analysis routine for SPI extreme values
  !! @param spi_data SPI time series data [nlon, nlat, ntime, nperiods]
  !! @param threshold_quantile Quantile for threshold selection (default 0.95)
  !! @param return_periods Array of return periods to compute [years]
  !! @param[out] thresholds Computed thresholds [nlon, nlat, nperiods]
  !! @param[out] gpd_shape GPD shape parameters [nlon, nlat, nperiods]
  !! @param[out] gpd_scale GPD scale parameters [nlon, nlat, nperiods]  
  !! @param[out] return_levels Return levels [nlon, nlat, nperiods, nreturn]
  subroutine compute_evt_analysis(spi_data, nlon, nlat, ntime, nperiods, &
                                  threshold_quantile, return_periods, nreturn, &
                                  thresholds, gpd_shape, gpd_scale, return_levels)
    integer,  intent(in)  :: nlon, nlat, ntime, nperiods, nreturn
    real(wp), intent(in)  :: spi_data(nlon, nlat, ntime, nperiods)
    real(wp), intent(in)  :: threshold_quantile
    real(wp), intent(in)  :: return_periods(nreturn)
    real(wp), intent(out) :: thresholds(nlon, nlat, nperiods)
    real(wp), intent(out) :: gpd_shape(nlon, nlat, nperiods)
    real(wp), intent(out) :: gpd_scale(nlon, nlat, nperiods)
    real(wp), intent(out) :: return_levels(nlon, nlat, nperiods, nreturn)
    
    integer :: i, j, ap, ir
    real(wp), allocatable :: extremes(:)
    real(wp) :: threshold, shape, scale
    logical :: fit_ok
    
    ! Initialize outputs
    thresholds = ieee_value(0.0_wp, ieee_quiet_nan)
    gpd_shape = ieee_value(0.0_wp, ieee_quiet_nan)
    gpd_scale = ieee_value(0.0_wp, ieee_quiet_nan)
    return_levels = ieee_value(0.0_wp, ieee_quiet_nan)
    
    do ap = 1, nperiods
      do i = 1, nlon
        do j = 1, nlat
          ! Extract valid SPI values for this grid point and period
          call extract_valid_spi(spi_data(i,j,:,ap), ntime, extremes)
          
          if (size(extremes) > MIN_EXCEEDANCES) then
            ! Determine threshold
            call fit_gpd_threshold(extremes, threshold_quantile, threshold, fit_ok)
            
            if (fit_ok) then
              thresholds(i,j,ap) = threshold
              
              ! Fit GPD to exceedances
              call fit_gpd_to_exceedances(extremes, threshold, shape, scale, fit_ok)
              
              if (fit_ok) then
                gpd_shape(i,j,ap) = shape
                gpd_scale(i,j,ap) = scale
                
                ! Compute return levels
                do ir = 1, nreturn
                  call compute_return_level(threshold, shape, scale, &
                                          return_periods(ir), size(extremes), &
                                          return_levels(i,j,ap,ir))
                end do
              end if
            end if
          end if
          
          if (allocated(extremes)) deallocate(extremes)
        end do
      end do
    end do
    
    print *, "> EVT analysis completed for", nperiods, "SPI period(s)"
  end subroutine compute_evt_analysis

  !> @brief Extract valid (non-NaN) SPI values
  subroutine extract_valid_spi(spi_series, ntime, valid_data)
    integer,  intent(in)  :: ntime
    real(wp), intent(in)  :: spi_series(ntime)
    real(wp), allocatable, intent(out) :: valid_data(:)
    
    integer :: t, count_valid
    
    ! Count valid values
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
  end subroutine extract_valid_spi

  !> @brief Determine threshold for GPD fitting using quantile method
  subroutine fit_gpd_threshold(data, quantile, threshold, success)
    real(wp), intent(in)  :: data(:)
    real(wp), intent(in)  :: quantile
    real(wp), intent(out) :: threshold
    logical,  intent(out) :: success
    
    real(wp), allocatable :: sorted_data(:)
    integer :: n, threshold_idx
    
    n = size(data)
    success = .false.
    threshold = ieee_value(0.0_wp, ieee_quiet_nan)
    
    if (n < MIN_EXCEEDANCES) return
    
    allocate(sorted_data(n))
    sorted_data = data
    call quicksort(sorted_data, 1, n)
    
    threshold_idx = int(quantile * real(n, wp))
    threshold_idx = max(1, min(n, threshold_idx))
    
    threshold = sorted_data(threshold_idx)
    
    ! Check if enough exceedances
    if ((n - threshold_idx + 1) >= MIN_EXCEEDANCES) then
      success = .true.
    end if
    
    deallocate(sorted_data)
  end subroutine fit_gpd_threshold

  !> @brief Fit GPD to exceedances above threshold
  subroutine fit_gpd_to_exceedances(data, threshold, shape, scale, success)
    real(wp), intent(in)  :: data(:)
    real(wp), intent(in)  :: threshold
    real(wp), intent(out) :: shape, scale
    logical,  intent(out) :: success
    
    real(wp), allocatable :: exceedances(:)
    integer :: i, n_exceed
    
    success = .false.
    shape = ieee_value(0.0_wp, ieee_quiet_nan)
    scale = ieee_value(0.0_wp, ieee_quiet_nan)
    
    ! Extract exceedances
    n_exceed = count(data > threshold)
    if (n_exceed < MIN_EXCEEDANCES) return
    
    allocate(exceedances(n_exceed))
    n_exceed = 0
    do i = 1, size(data)
      if (data(i) > threshold) then
        n_exceed = n_exceed + 1
        exceedances(n_exceed) = data(i) - threshold
      end if
    end do
    
    ! Fit GPD using method of moments as fallback
    call fit_gpd_moments(exceedances, shape, scale, success)
    
    deallocate(exceedances)
  end subroutine fit_gpd_to_exceedances

  !> @brief Fit GPD using method of moments
  subroutine fit_gpd_moments(exceedances, shape, scale, success)
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
      
      ! Basic validation
      if (scale > 0.0_wp .and. abs(shape) < 0.5_wp) then
        success = .true.
      end if
    end if
  end subroutine fit_gpd_moments

  !> @brief Compute return level for given return period
  subroutine compute_return_level(threshold, shape, scale, return_period, n_years, return_level)
    real(wp), intent(in)  :: threshold, shape, scale, return_period
    integer,  intent(in)  :: n_years
    real(wp), intent(out) :: return_level
    
    real(wp) :: lambda, prob
    
    return_level = ieee_value(0.0_wp, ieee_quiet_nan)
    
    if (scale <= 0.0_wp .or. n_years <= 0) return
    
    ! Estimate exceedance rate
    lambda = 1.0_wp  ! Assuming annual data
    
    ! Compute return level
    prob = 1.0_wp - 1.0_wp / (lambda * return_period)
    
    if (abs(shape) < 1.0e-6_wp) then
      ! Exponential case
      return_level = threshold - scale * log(1.0_wp - prob)
    else
      ! General GPD case
      return_level = threshold + (scale / shape) * ((1.0_wp - prob)**(-shape) - 1.0_wp)
    end if
  end subroutine compute_return_level

  !> @brief Extract extreme values (both minima and maxima)
  subroutine extract_extremes(spi_data, ntime, block_size, extremes_min, extremes_max)
    real(wp), intent(in)  :: spi_data(ntime)
    integer,  intent(in)  :: ntime, block_size
    real(wp), allocatable, intent(out) :: extremes_min(:), extremes_max(:)
    
    integer :: n_blocks, block, start_idx, end_idx, i
    real(wp) :: block_min, block_max
    logical :: has_valid
    
    n_blocks = ntime / block_size
    allocate(extremes_min(n_blocks), extremes_max(n_blocks))
    
    extremes_min = ieee_value(0.0_wp, ieee_quiet_nan)
    extremes_max = ieee_value(0.0_wp, ieee_quiet_nan)
    
    do block = 1, n_blocks
      start_idx = (block - 1) * block_size + 1
      end_idx = min(block * block_size, ntime)
      
      has_valid = .false.
      block_min = huge(1.0_wp)
      block_max = -huge(1.0_wp)
      
      do i = start_idx, end_idx
        if (.not. ieee_is_nan(spi_data(i))) then
          has_valid = .true.
          block_min = min(block_min, spi_data(i))
          block_max = max(block_max, spi_data(i))
        end if
      end do
      
      if (has_valid) then
        extremes_min(block) = block_min
        extremes_max(block) = block_max
      end if
    end do
  end subroutine extract_extremes

  !> @brief Validate GPD fit quality
  logical function validate_gpd_fit(exceedances, shape, scale) result(is_valid)
    real(wp), intent(in) :: exceedances(:), shape, scale
    
    real(wp) :: ks_stat
    integer :: n
    
    is_valid = .false.
    n = size(exceedances)
    
    if (n < MIN_EXCEEDANCES .or. scale <= 0.0_wp) return
    
    ! Basic parameter checks
    if (abs(shape) > 0.5_wp) return  ! Conservative shape bound
    
    ! Could add Kolmogorov-Smirnov test here
    ! For now, use basic validation
    is_valid = .true.
  end function validate_gpd_fit

  !> @brief Compute return periods for given values
  subroutine compute_return_periods(values, threshold, shape, scale, n_years, return_periods)
    real(wp), intent(in)  :: values(:)
    real(wp), intent(in)  :: threshold, shape, scale
    integer,  intent(in)  :: n_years
    real(wp), intent(out) :: return_periods(size(values))
    
    integer :: i
    real(wp) :: lambda, prob
    
    lambda = 1.0_wp  ! Annual exceedance rate
    
    do i = 1, size(values)
      if (values(i) > threshold .and. scale > 0.0_wp) then
        if (abs(shape) < 1.0e-6_wp) then
          prob = exp(-(values(i) - threshold) / scale)
        else
          prob = (1.0_wp + shape * (values(i) - threshold) / scale)**(-1.0_wp/shape)
        end if
        return_periods(i) = 1.0_wp / (lambda * (1.0_wp - prob))
      else
        return_periods(i) = ieee_value(0.0_wp, ieee_quiet_nan)
      end if
    end do
  end subroutine compute_return_periods

  !> @brief Simple quicksort implementation
  recursive subroutine quicksort(array, low, high)
    real(wp), intent(inout) :: array(:)
    integer,  intent(in)    :: low, high
    integer :: pivot_idx
    
    if (low < high) then
      call partition(array, low, high, pivot_idx)
      call quicksort(array, low, pivot_idx - 1)
      call quicksort(array, pivot_idx + 1, high)
    end if
  end subroutine quicksort

  !> @brief Partition for quicksort
  subroutine partition(array, low, high, pivot_idx)
    real(wp), intent(inout) :: array(:)
    integer,  intent(in)    :: low, high
    integer,  intent(out)   :: pivot_idx
    
    real(wp) :: pivot, temp
    integer :: i, j
    
    pivot = array(high)
    i = low - 1
    
    do j = low, high - 1
      if (array(j) <= pivot) then
        i = i + 1
        temp = array(i)
        array(i) = array(j)
        array(j) = temp
      end if
    end do
    
    temp = array(i + 1)
    array(i + 1) = array(high)
    array(high) = temp
    
    pivot_idx = i + 1
  end subroutine partition

end module evt_mod
