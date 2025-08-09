!> @file prep_mod.f90
!! Preprocess precipitation: units, calendar, gap check, gap filling, rolling sums.
!! Works with io_mod outputs:
!!   - Preferred: pass valid_time(:), time_units, calendar
!!   - Fallback : pass start_year, start_month
module prep_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use iso_fortran_env, only: wp => real64, int64
  implicit none
  private
  public :: preprocess_precip

contains

  subroutine preprocess_precip(precip_in, nlon, nlat, ntime,                   &
                               precip_mm,                                      &
                               year_idx, month_idx,                            &
                               valid_time, time_units, calendar,               &
                               start_year, start_month,                        &
                               agg3, agg12)
    ! ---- Required ----
    integer,  intent(in)  :: nlon, nlat, ntime
    real(wp), intent(in)  :: precip_in(nlon, nlat, ntime)     ! m/month
    real(wp), intent(out) :: precip_mm(nlon, nlat, ntime)     ! mm/month

    ! ---- Optional time info from io_mod ----
    integer(int64), intent(in),  optional :: valid_time(ntime)  ! seconds since base date
    character(len=*), intent(in), optional :: time_units        ! "seconds since YYYY-MM-DD"
    character(len=*), intent(in), optional :: calendar          ! e.g., "proleptic_gregorian"

    ! ---- Optional fallback origin ----
    integer, intent(in),  optional :: start_year, start_month

    ! ---- Optional outputs ----
    integer, intent(out), optional :: year_idx(ntime), month_idx(ntime)
    real(wp), intent(out), optional :: agg3(nlon, nlat, ntime)
    real(wp), intent(out), optional :: agg12(nlon, nlat, ntime)

    ! ---- Locals ----
    integer :: i, j, t, m, gaps
    real(wp), allocatable :: clim_mean(:,:,:)
    integer :: y0, m0
    logical :: have_time_meta, have_y_m

    ! 1) Units: m -> mm
    precip_mm = precip_in * 1000.0_wp

    ! 2) Calendar build
    have_time_meta = present(valid_time) .and. present(time_units)
    if (have_time_meta) then
      call build_calendar_from_units(valid_time, time_units, calendar, ntime,    &
                                     y0, m0, year_idx, month_idx)
      write(*,*) "  Calendar: built from valid_time + units (", trim(time_units), ")."
    else
      if (present(start_year)) then
        y0 = start_year
      else
        y0 = 1980
      end if
      if (present(start_month)) then
        m0 = start_month
      else
        m0 = 1
      end if
      if (present(year_idx) .and. present(month_idx)) then
        do t = 1, ntime
          year_idx(t)  = y0 + ((m0 - 1) + (t - 1)) / 12
          month_idx(t) = 1  + mod((m0 - 1) + (t - 1), 12)
        end do
      end if
      write(*,*) "  Calendar: fallback from start_year/start_month."
    end if
    have_y_m = present(year_idx) .and. present(month_idx)

    ! 3) Gap check (only if Y/M available)
    gaps = 0
    if (have_y_m) then
      do t = 2, ntime
        if (.not. is_next_month(year_idx(t-1), month_idx(t-1), year_idx(t), month_idx(t))) then
          gaps = gaps + 1
          write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
             "  WARNING: gap near t=", t, " (", year_idx(t-1), "-", month_idx(t-1), &
             ") -> (", year_idx(t), "-", month_idx(t), ")"
        end if
      end do
      if (gaps == 0) then
        write(*,*) "  Gap check: no obvious calendar gaps."
      else
        write(*,*) "  Gap check: gaps found =", gaps, "(SPI windows skip NaN windows)."
      end if
    else
      write(*,*) "  Gap check: skipped (no explicit Y/M)."
    end if

    ! 4) Monthly climatology per (lon,lat,month)
    allocate(clim_mean(nlon, nlat, 12)); clim_mean = 0.0_wp
    if (have_y_m) then
      do m = 1, 12
        do i = 1, nlon
          do j = 1, nlat
            call monthly_mean_by_index(precip_mm(i,j,:), ntime, month_idx, m, clim_mean(i,j,m))
          end do
        end do
      end do
    else
      write(*,*) "  NOTE: climatology uses modulo(12) since month_idx not provided."
      do m = 1, 12
        do i = 1, nlon
          do j = 1, nlat
            call monthly_mean_mod12(precip_mm(i,j,:), ntime, m, clim_mean(i,j,m))
          end do
        end do
      end do
    end if

    ! 5) Fill NaNs with monthly climatology
    do t = 1, ntime
      if (have_y_m) then
        m = month_idx(t)
      else
        m = mod(t-1, 12) + 1
      end if
      do i = 1, nlon
        do j = 1, nlat
          if (ieee_is_nan(precip_mm(i,j,t))) precip_mm(i,j,t) = clim_mean(i,j,m)
        end do
      end do
    end do
    deallocate(clim_mean)

    ! 6) Optional rolling sums
    if (present(agg3))  call rolling_sum_mm(precip_mm, nlon, nlat, ntime, 3,  agg3)
    if (present(agg12)) call rolling_sum_mm(precip_mm, nlon, nlat, ntime, 12, agg12)

    write(*,*) "> Preprocessing complete: units, calendar, gap check, fill, (optional) aggregates."
  end subroutine preprocess_precip

  ! ---------- helpers ----------

  subroutine monthly_mean_by_index(x, ntime, month_idx, m_target, mean_out)
    real(wp), intent(in)  :: x(ntime)
    integer,  intent(in)  :: ntime, month_idx(ntime), m_target
    real(wp), intent(out) :: mean_out
    integer :: t, c
    real(wp) :: s
    s = 0.0_wp; c = 0
    do t = 1, ntime
      if (month_idx(t) == m_target .and. .not. ieee_is_nan(x(t))) then
        s = s + x(t); c = c + 1
      end if
    end do
    if (c > 0) then
      mean_out = s / real(c, wp)
    else
      mean_out = ieee_value(0.0_wp, ieee_quiet_nan)
    end if
  end subroutine monthly_mean_by_index

  subroutine monthly_mean_mod12(x, ntime, month1_12, mean_out)
    real(wp), intent(in)  :: x(ntime)
    integer,  intent(in)  :: ntime, month1_12
    real(wp), intent(out) :: mean_out
    integer :: t, c
    real(wp) :: s
    s = 0.0_wp; c = 0
    do t = month1_12, ntime, 12
      if (.not. ieee_is_nan(x(t))) then
        s = s + x(t); c = c + 1
      end if
    end do
    if (c > 0) then
      mean_out = s / real(c, wp)
    else
      mean_out = ieee_value(0.0_wp, ieee_quiet_nan)
    end if
  end subroutine monthly_mean_mod12

  subroutine rolling_sum_mm(x, nlon, nlat, ntime, k, outsum)
    real(wp), intent(in)  :: x(nlon, nlat, ntime)
    integer,  intent(in)  :: nlon, nlat, ntime, k
    real(wp), intent(out) :: outsum(nlon, nlat, ntime)
    integer :: i, j, t
    outsum = ieee_value(0.0_wp, ieee_quiet_nan)
    if (ntime < k) return
    do t = k, ntime
      do i = 1, nlon
        do j = 1, nlat
          if (all(.not. ieee_is_nan(x(i,j,t-k+1:t)))) then
            outsum(i,j,t) = sum(x(i,j,t-k+1:t))
          end if
        end do
      end do
    end do
  end subroutine rolling_sum_mm

  logical pure function is_next_month(y1, m1, y2, m2)
    integer, intent(in) :: y1, m1, y2, m2
    is_next_month = .false.
    if (y2 == y1 .and. m2 == m1 + 1) is_next_month = .true.
    if (m1 == 12 .and. m2 == 1 .and. y2 == y1 + 1) is_next_month = .true.
  end function is_next_month

  subroutine build_calendar_from_units(valid_time, units, cal, ntime, y0, m0, year_idx, month_idx)
    integer(int64), intent(in) :: valid_time(ntime)
    character(len=*), intent(in) :: units
    character(len=*), intent(in), optional :: cal
    integer, intent(in)  :: ntime
    integer, intent(out) :: y0, m0
    integer, intent(out), optional :: year_idx(ntime), month_idx(ntime)
    integer :: ybase, mb, db, y, m, d, t
    call parse_units_epoch(units, ybase, mb, db)
    y0 = ybase; m0 = mb
    if (present(year_idx) .and. present(month_idx)) then
      do t = 1, ntime
        call epoch_seconds_to_ymd(ybase, mb, db, valid_time(t), y, m, d)
        year_idx(t)  = y
        month_idx(t) = m
      end do
    end if
  end subroutine build_calendar_from_units

  subroutine parse_units_epoch(units, y0, m0, d0)
    character(len=*), intent(in)  :: units
    integer, intent(out) :: y0, m0, d0
    integer :: p, ios
    y0 = 1970; m0 = 1; d0 = 1
    p = index(units, "since")
    if (p > 0) then
      read(units(p+6:), '(I4,1x,I2,1x,I2)', iostat=ios) y0, m0, d0
      if (ios /= 0) then
        y0 = 1970; m0 = 1; d0 = 1
      end if
    end if
  end subroutine parse_units_epoch

  subroutine epoch_seconds_to_ymd(ybase, mbase, dbase, seconds, y, m, d)
    integer, intent(in)        :: ybase, mbase, dbase
    integer(int64), intent(in) :: seconds
    integer, intent(out)       :: y, m, d
    integer(int64) :: days_left
    integer :: dim(12)
    dim = [31,28,31,30,31,30,31,31,30,31,30,31]
    y = ybase; m = mbase; d = dbase
    days_left = seconds / 86400_int64
    do while (days_left > 0_int64)
      if (is_leap(y)) then
        dim(2) = 29
      else
        dim(2) = 28
      end if
      if (d < dim(m)) then
        d = d + 1
      else
        d = 1
        m = m + 1
        if (m > 12) then
          m = 1; y = y + 1
        end if
      end if
      days_left = days_left - 1_int64
    end do
  end subroutine epoch_seconds_to_ymd

  logical pure function is_leap(y)
    integer, intent(in) :: y
    is_leap = (mod(y,4) == 0 .and. (mod(y,100) /= 0 .or. mod(y,400) == 0))
  end function is_leap

end module prep_mod