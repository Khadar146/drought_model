!> @file spi_mod.f90
!! @brief Computes Standardized Precipitation Index (SPI) on gridded data.
!! @details Implements SPI for multiple accumulation periods using gamma fitting
!!          and standard normal transformation, across all grid cells.
!!          Best-practice approach following WMO guidelines.
!!
!! @date 2025
!! @author Khadar

module spi_mod
  use iso_fortran_env, only: wp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private
  public :: compute_spi

contains

  subroutine compute_spi(precip_mm, nlon, nlat, ntime, accum_periods, num_periods, spi)
    !! Computes SPI for all grid cells and accumulation periods
    real(wp), intent(in) :: precip_mm(nlon, nlat, ntime)
    integer,  intent(in) :: nlon, nlat, ntime
    integer,  intent(in) :: accum_periods(:)      ! E.g., [1, 3, 6, 12]
    integer,  intent(in) :: num_periods
    real(wp), intent(out):: spi(nlon, nlat, ntime, num_periods)

    integer :: i, j, t, k, ap, m
    real(wp), allocatable :: acc_p(:)
    real(wp), allocatable :: shape(:,:,:), scale(:,:,:)

    allocate(shape(nlon, nlat, 12))
    allocate(scale(nlon, nlat, 12))
    shape = ieee_value(0.0_wp, ieee_quiet_nan)
    scale = ieee_value(0.0_wp, ieee_quiet_nan)
    spi   = ieee_value(0.0_wp, ieee_quiet_nan)

    do ap = 1, num_periods
      do i = 1, nlon
        do j = 1, nlat
          ! Accumulated precipitation
          call accumulate_series(precip_mm(i,j,:), accum_periods(ap), acc_p)

          ! Fit gamma parameters by calendar month (m=1..12)
          do m = 1, 12
            call fit_gamma_by_month(acc_p, m, ntime, shape(i,j,m), scale(i,j,m))
          end do

          ! Transform to SPI
          do t = accum_periods(ap), ntime
            m = mod(t-1, 12) + 1
            if (.not. ieee_is_nan(acc_p(t))) then
              spi(i,j,t,ap) = gamma_to_spi(acc_p(t), shape(i,j,m), scale(i,j,m))
            end if
          end do
        end do
      end do
    end do

    deallocate(acc_p, shape, scale)
    print *, "> SPI computation complete for all grid cells and time steps."

  end subroutine compute_spi

  subroutine accumulate_series(series, period, result)
    !! Computes moving sums (accumulated precip)
    real(wp), intent(in)  :: series(:)
    integer,  intent(in)  :: period
    real(wp), allocatable, intent(out) :: result(:)

    integer :: t
    integer :: ntime
    ntime = size(series)
    allocate(result(ntime))
    result = ieee_value(0.0_wp, ieee_quiet_nan)

    do t = period, ntime
      if (all(.not. ieee_is_nan(series(t-period+1:t)))) then
        result(t) = sum(series(t-period+1:t))
      end if
    end do
  end subroutine accumulate_series

  subroutine fit_gamma_by_month(data, month_idx, ntime, shape, scale)
    !! Fits gamma distribution using method of moments for a given month
    real(wp), intent(in) :: data(:)
    integer,  intent(in) :: month_idx, ntime
    real(wp), intent(out):: shape, scale

    real(wp) :: mean, var, logmean
    integer :: t, count
    real(wp), allocatable :: month_vals(:)

    allocate(month_vals(ntime / 12 + 1))
    count = 0

    do t = month_idx, ntime, 12
      if (.not. ieee_is_nan(data(t))) then
        count = count + 1
        month_vals(count) = data(t)
      end if
    end do

    if (count >= 5) then
      mean = sum(month_vals(1:count)) / count
      var  = sum((month_vals(1:count) - mean)**2) / (count - 1)
      if (var > 0.0_wp .and. mean > 0.0_wp) then
        shape = (mean**2) / var
        scale = var / mean
      else
        shape = ieee_value(0.0_wp, ieee_quiet_nan)
        scale = shape
      end if
    else
      shape = ieee_value(0.0_wp, ieee_quiet_nan)
      scale = shape
    end if

    deallocate(month_vals)
  end subroutine fit_gamma_by_month

  function gamma_to_spi(value, shape, scale) result(spi_val)
    !! Computes SPI from gamma parameters and value
    real(wp), intent(in) :: value, shape, scale
    real(wp) :: spi_val, prob

    if (ieee_is_nan(shape) .or. ieee_is_nan(scale)) then
      spi_val = ieee_value(0.0_wp, ieee_quiet_nan)
    else if (value == 0.0_wp) then
      spi_val = -2.0_wp ! SPI near lowest bound
    else
      prob = gammainc(value / scale, shape)
      spi_val = inverse_normal(prob)
    end if
  end function gamma_to_spi

 !> @file spi_mod.f90
!! @brief Computes Standardized Precipitation Index (SPI) on gridded data.
!! @details Implements SPI for multiple accumulation periods using gamma fitting
!!          and standard normal transformation, across all grid cells.
!!          Best-practice approach following WMO guidelines.
!!
!! @date 2025
!! @author Khadar

module spi_mod
  use iso_fortran_env, only: wp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use lib_math_mod, only: gammainc_regularized, norminv
  implicit none
  private
  public :: compute_spi

contains

  subroutine compute_spi(precip_mm, nlon, nlat, ntime, accum_periods, num_periods, spi)
    !! Computes SPI for all grid cells and accumulation periods
    real(wp), intent(in) :: precip_mm(nlon, nlat, ntime)
    integer,  intent(in) :: nlon, nlat, ntime
    integer,  intent(in) :: accum_periods(:)
    integer,  intent(in) :: num_periods
    real(wp), intent(out):: spi(nlon, nlat, ntime, num_periods)

    integer :: i, j, t, k, ap, m
    real(wp), allocatable :: acc_p(:)
    real(wp), allocatable :: shape(:,:,:), scale(:,:,:)

    allocate(shape(nlon, nlat, 12))
    allocate(scale(nlon, nlat, 12))
    shape = ieee_value(0.0_wp, ieee_quiet_nan)
    scale = ieee_value(0.0_wp, ieee_quiet_nan)
    spi   = ieee_value(0.0_wp, ieee_quiet_nan)

    do ap = 1, num_periods
      do i = 1, nlon
        do j = 1, nlat
          call accumulate_series(precip_mm(i,j,:), accum_periods(ap), acc_p)

          do m = 1, 12
            call fit_gamma_by_month(acc_p, m, ntime, shape(i,j,m), scale(i,j,m))
          end do

          do t = accum_periods(ap), ntime
            m = mod(t-1, 12) + 1
            if (.not. ieee_is_nan(acc_p(t))) then
              spi(i,j,t,ap) = gamma_to_spi(acc_p(t), shape(i,j,m), scale(i,j,m))
            end if
          end do
        end do
      end do
    end do

    deallocate(acc_p, shape, scale)
    print *, "> SPI computation complete for all grid cells and time steps."

  end subroutine compute_spi

  subroutine accumulate_series(series, period, result)
    real(wp), intent(in)  :: series(:)
    integer,  intent(in)  :: period
    real(wp), allocatable, intent(out) :: result(:)

    integer :: t, ntime
    ntime = size(series)
    allocate(result(ntime))
    result = ieee_value(0.0_wp, ieee_quiet_nan)

    do t = period, ntime
      if (all(.not. ieee_is_nan(series(t-period+1:t)))) then
        result(t) = sum(series(t-period+1:t))
      end if
    end do
  end subroutine accumulate_series

  subroutine fit_gamma_by_month(data, month_idx, ntime, shape, scale)
    real(wp), intent(in) :: data(:)
    integer,  intent(in) :: month_idx, ntime
    real(wp), intent(out):: shape, scale

    real(wp) :: mean, var
    integer :: t, count
    real(wp), allocatable :: month_vals(:)

    allocate(month_vals(ntime / 12 + 1))
    count = 0

    do t = month_idx, ntime, 12
      if (.not. ieee_is_nan(data(t))) then
        count = count + 1
        month_vals(count) = data(t)
      end if
    end do

    if (count >= 5) then
      mean = sum(month_vals(1:count)) / count
      var  = sum((month_vals(1:count) - mean)**2) / (count - 1)
      if (var > 0.0_wp .and. mean > 0.0_wp) then
        shape = (mean**2) / var
        scale = var / mean
      else
        shape = ieee_value(0.0_wp, ieee_quiet_nan)
        scale = shape
      end if
    else
      shape = ieee_value(0.0_wp, ieee_quiet_nan)
      scale = shape
    end if

    deallocate(month_vals)
  end subroutine fit_gamma_by_month

  function gamma_to_spi(value, shape, scale) result(spi_val)
    real(wp), intent(in) :: value, shape, scale
    real(wp) :: spi_val, prob

    if (ieee_is_nan(shape) .or. ieee_is_nan(scale)) then
      spi_val = ieee_value(0.0_wp, ieee_quiet_nan)
    else if (value == 0.0_wp) then
      spi_val = -2.0_wp
    else
      prob = gammainc_regularized(value / scale, shape)
      spi_val = norminv(prob)
    end if
  end function gamma_to_spi

end module spi_mod