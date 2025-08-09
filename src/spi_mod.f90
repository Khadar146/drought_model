!> @file spi_mod.f90
!! @brief Gridded SPI using Gamma fit (per month) and FSML transforms.
!! @details
!!  - Accumulates monthly precip (mm) over k = 1/3/6/12...
!!  - Fits Gamma(shape, scale) per calendar month using ONLY the reference period.
!!  - Mixture with zero-precip probability p_zero.
!!  - Uses FSML (fsml_dst): f_dst_gammai_core (lower incomplete gamma),
!!                          f_dst_norm_ppf_core (inverse normal CDF).
!!
!! Inputs:
!!   precip_mm(nlon,nlat,ntime)    : monthly precip in mm (NaNs allowed)
!!   start_year, start_month       : index 1 corresponds to this (e.g., 1980,1)
!!   ref_start_year, ref_end_year  : fit window (e.g., 1991–2020)
!!   accum_periods(:)              : e.g., [1,3,6,12]
!!
!! Output:
!!   spi(nlon,nlat,ntime,num_periods) : SPI values (NaN where insufficient data)
!!
!! @author Khadar
!! @date 2025


module spi_mod
  use iso_fortran_env, only: wp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use fsml_dst,  only: f_dst_gammai_core, f_dst_norm_ppf_core
  implicit none
  private
  public :: compute_spi
contains

  !------------------------------!
  ! Public driver: all periods   !
  !------------------------------!
  subroutine compute_spi(precip_mm, nlon, nlat, ntime,                 &
                         start_year, start_month,                       &
                         ref_start_year, ref_end_year,                  &
                         accum_periods, num_periods,                    &
                         spi,                                           &
                         year_idx, month_idx)
    integer,  intent(in) :: nlon, nlat, ntime
    real(wp), intent(in) :: precip_mm(nlon, nlat, ntime)
    integer,  intent(in) :: start_year, start_month
    integer,  intent(in) :: ref_start_year, ref_end_year
    integer,  intent(in) :: accum_periods(:), num_periods
    real(wp), intent(out):: spi(nlon, nlat, ntime, num_periods)
    ! Optional calendar arrays from preprocess (preferred if present)
    integer,  intent(in), optional :: year_idx(ntime), month_idx(ntime)

    integer :: ap, i, j, k
    real(wp), allocatable :: acc_series(:,:,:)
    logical,  allocatable :: ref_mask(:)
    logical :: have_calendar

    spi = ieee_value(0.0_wp, ieee_quiet_nan)

    have_calendar = present(year_idx) .and. present(month_idx)

    allocate(ref_mask(ntime))
    if (have_calendar) then
      call build_reference_mask_from_arrays(ntime, year_idx, ref_start_year, ref_end_year, ref_mask)
    else
      call build_reference_mask_linear(ntime, start_year, start_month, ref_start_year, ref_end_year, ref_mask)
    end if

    do ap = 1, num_periods
      k = accum_periods(ap)
      call accumulate_field(precip_mm, nlon, nlat, ntime, k, acc_series)

      do i = 1, nlon
        do j = 1, nlat
          if (have_calendar) then
            call spi_for_cell_with_arrays(acc_series(i,j,:), ntime, ref_mask, month_idx, k, spi(i,j,:,ap))
          else
            call spi_for_cell_linear     (acc_series(i,j,:), ntime, ref_mask, start_month, k, spi(i,j,:,ap))
          end if
        end do
      end do

      deallocate(acc_series)
    end do

    deallocate(ref_mask)
    print *, "> SPI computed for", num_periods, "accumulation period(s)."
  end subroutine compute_spi

  !-----------------------------------------!
  ! Accumulate (moving sum) the whole field !
  !-----------------------------------------!
  subroutine accumulate_field(x, nlon, nlat, ntime, k, acc)
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
  end subroutine accumulate_field

  !----------------------------------------------!
  ! SPI cell using explicit month_idx (preferred) !
  !----------------------------------------------!
  subroutine spi_for_cell_with_arrays(acc, ntime, ref_mask, month_idx, k, spi_cell)
    real(wp), intent(in)  :: acc(ntime)
    integer,  intent(in)  :: ntime, k
    logical,  intent(in)  :: ref_mask(ntime)
    integer,  intent(in)  :: month_idx(ntime)
    real(wp), intent(out) :: spi_cell(ntime)

    real(wp) :: shape(12), scale(12), p_zero(12)
    logical  :: okfit(12)
    integer  :: m, t

    spi_cell = ieee_value(0.0_wp, ieee_quiet_nan)
    shape    = ieee_value(0.0_wp, ieee_quiet_nan)
    scale    = ieee_value(0.0_wp, ieee_quiet_nan)
    p_zero   = 0.0_wp
    okfit    = .false.

    do m = 1, 12
      call fit_gamma_ref_month_arrays(acc, ntime, ref_mask, month_idx, m, k, shape(m), scale(m), p_zero(m), okfit(m))
    end do

    do t = k, ntime
      m = month_idx(t)
      if (.not. ieee_is_nan(acc(t)) .and. okfit(m)) then
        spi_cell(t) = gamma_to_spi(acc(t), shape(m), scale(m), p_zero(m))
      end if
    end do
  end subroutine spi_for_cell_with_arrays

  !---------------------------------------------------!
  ! SPI cell using linear month from start_month      !
  ! (kept for backward compatibility)                 !
  !---------------------------------------------------!
  subroutine spi_for_cell_linear(acc, ntime, ref_mask, start_month, k, spi_cell)
    real(wp), intent(in)  :: acc(ntime)
    integer,  intent(in)  :: ntime, start_month, k
    logical,  intent(in)  :: ref_mask(ntime)
    real(wp), intent(out) :: spi_cell(ntime)

    real(wp) :: shape(12), scale(12), p_zero(12)
    logical  :: okfit(12)
    integer  :: m, t

    spi_cell = ieee_value(0.0_wp, ieee_quiet_nan)
    shape    = ieee_value(0.0_wp, ieee_quiet_nan)
    scale    = ieee_value(0.0_wp, ieee_quiet_nan)
    p_zero   = 0.0_wp
    okfit    = .false.

    do m = 1, 12
      call fit_gamma_ref_month_linear(acc, ntime, ref_mask, m, start_month, k, shape(m), scale(m), p_zero(m), okfit(m))
    end do

    do t = k, ntime
      m = month_of_index(t, start_month)
      if (.not. ieee_is_nan(acc(t)) .and. okfit(m)) then
        spi_cell(t) = gamma_to_spi(acc(t), shape(m), scale(m), p_zero(m))
      end if
    end do
  end subroutine spi_for_cell_linear

  !-------------------- fits (calendar arrays) --------------------!
  subroutine fit_gamma_ref_month_arrays(acc, ntime, ref_mask, month_idx, m_target, k, &
                                        shape, scale, p_zero, okfit)
    real(wp), intent(in)  :: acc(ntime)
    integer,  intent(in)  :: ntime, month_idx(ntime), m_target, k
    logical,  intent(in)  :: ref_mask(ntime)
    real(wp), intent(out) :: shape, scale, p_zero
    logical,  intent(out) :: okfit

    real(wp) :: mean_pos, var_pos, v
    integer  :: count_all, count_pos, t

    shape = ieee_value(0.0_wp, ieee_quiet_nan)
    scale = ieee_value(0.0_wp, ieee_quiet_nan)
    p_zero = 0.0_wp
    okfit  = .false.

    count_all = 0; count_pos = 0
    mean_pos  = 0.0_wp; var_pos = 0.0_wp

    do t = k, ntime
      if (.not. ref_mask(t)) cycle
      if (month_idx(t) /= m_target) cycle
      if (ieee_is_nan(acc(t))) cycle

      count_all = count_all + 1
      v = acc(t)
      if (v > 0.0_wp) then
        count_pos = count_pos + 1
        mean_pos  = mean_pos + v
      end if
    end do

    if (count_all > 0) p_zero = real(count_all - count_pos, wp) / real(count_all, wp)

    if (count_pos >= 5) then
      mean_pos = mean_pos / real(count_pos, wp)
      do t = k, ntime
        if (.not. ref_mask(t)) cycle
        if (month_idx(t) /= m_target) cycle
        if (ieee_is_nan(acc(t))) cycle
        v = acc(t)
        if (v > 0.0_wp) var_pos = var_pos + (v - mean_pos)**2
      end do
      if (count_pos > 1) var_pos = var_pos / real(count_pos - 1, wp)

      if (mean_pos > 0.0_wp .and. var_pos > 0.0_wp) then
        shape = (mean_pos**2) / var_pos
        scale = var_pos / mean_pos
        okfit = .true.
      end if
    end if
  end subroutine fit_gamma_ref_month_arrays

  !-------------------- fits (linear calendar) --------------------!
  subroutine fit_gamma_ref_month_linear(acc, ntime, ref_mask, m_target, start_month, k, &
                                        shape, scale, p_zero, okfit)
    real(wp), intent(in)  :: acc(ntime)
    integer,  intent(in)  :: ntime, m_target, start_month, k
    logical,  intent(in)  :: ref_mask(ntime)
    real(wp), intent(out) :: shape, scale, p_zero
    logical,  intent(out) :: okfit

    real(wp) :: mean_pos, var_pos, v
    integer  :: count_all, count_pos, t, m

    shape = ieee_value(0.0_wp, ieee_quiet_nan)
    scale = ieee_value(0.0_wp, ieee_quiet_nan)
    p_zero = 0.0_wp
    okfit  = .false.

    count_all = 0; count_pos = 0
    mean_pos  = 0.0_wp; var_pos = 0.0_wp

    do t = k, ntime
      if (.not. ref_mask(t)) cycle
      m = month_of_index(t, start_month)
      if (m /= m_target) cycle
      if (ieee_is_nan(acc(t))) cycle

      count_all = count_all + 1
      v = acc(t)
      if (v > 0.0_wp) then
        count_pos = count_pos + 1
        mean_pos  = mean_pos + v
      end if
    end do

    if (count_all > 0) p_zero = real(count_all - count_pos, wp) / real(count_all, wp)

    if (count_pos >= 5) then
      mean_pos = mean_pos / real(count_pos, wp)
      do t = k, ntime
        if (.not. ref_mask(t)) cycle
        m = month_of_index(t, start_month)
        if (m /= m_target) cycle
        if (ieee_is_nan(acc(t))) cycle
        v = acc(t)
        if (v > 0.0_wp) var_pos = var_pos + (v - mean_pos)**2
      end do
      if (count_pos > 1) var_pos = var_pos / real(count_pos - 1, wp)

      if (mean_pos > 0.0_wp .and. var_pos > 0.0_wp) then
        shape = (mean_pos**2) / var_pos
        scale = var_pos / mean_pos
        okfit = .true.
      end if
    end if
  end subroutine fit_gamma_ref_month_linear

  !----------------------------------------------!
  ! Gamma → probability → Normal quantile (SPI)  !
  !----------------------------------------------!
  pure function gamma_to_spi(value, shape, scale, p_zero) result(spi_val)
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
  end function gamma_to_spi

  !---------------------------!
  ! Build reference time mask !
  !---------------------------!
  subroutine build_reference_mask_from_arrays(ntime, year_idx, y1, y2, ref_mask)
    integer, intent(in)  :: ntime, year_idx(ntime), y1, y2
    logical, intent(out) :: ref_mask(ntime)
    integer :: t
    do t = 1, ntime
      ref_mask(t) = (year_idx(t) >= y1 .and. year_idx(t) <= y2)
    end do
  end subroutine build_reference_mask_from_arrays

  subroutine build_reference_mask_linear(ntime, start_year, start_month, y1, y2, ref_mask)
    integer, intent(in)  :: ntime, start_year, start_month, y1, y2
    logical, intent(out) :: ref_mask(ntime)
    integer :: t, y, m
    do t = 1, ntime
      call year_month_of_index(t, start_year, start_month, y, m)
      ref_mask(t) = (y >= y1 .and. y <= y2)
    end do
  end subroutine build_reference_mask_linear

  !------------------------------!
  ! Helpers: calendar indexing   !
  !------------------------------!
  pure subroutine year_month_of_index(t, start_year, start_month, y, m)
    integer, intent(in)  :: t, start_year, start_month
    integer, intent(out) :: y, m
    integer :: z
    z = (start_month - 1) + (t - 1)
    y = start_year + z / 12
    m = 1 + mod(z, 12)
  end subroutine year_month_of_index

  pure function month_of_index(t, start_month) result(m)
    integer, intent(in) :: t, start_month
    integer :: m, z
    z = (start_month - 1) + (t - 1)
    m = 1 + mod(z, 12)
  end function month_of_index

end module spi_mod