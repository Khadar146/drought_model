    ! =====================================================================
! File: src/prep_mod.f90
! Purpose: PREP stage — units + calendar + sanity checks. No filling.
! API: preprocess_precip(...) matches main_historical.
! =====================================================================
module prep_mod
  use iso_fortran_env, only: wp => real64, int64
  implicit none
  private
  public :: preprocess_precip
  public :: build_calendar, check_calendar  ! exposed for optional QA

contains

  subroutine preprocess_precip(precip_in, nlon, nlat, ntime,             &
                               precip_mm,                                 &
                               year_idx, month_idx,                       &
                               valid_time, time_units, calendar)
    ! Required inputs
    integer,          intent(in)  :: nlon, nlat, ntime
    real(wp),         intent(in)  :: precip_in(nlon, nlat, ntime)   ! meters/month
    integer(int64),   intent(in)  :: valid_time(ntime)
    character(len=*), intent(in)  :: time_units
    character(len=*), intent(in)  :: calendar                        ! informational

    ! Outputs
    real(wp),         intent(out) :: precip_mm(nlon, nlat, ntime)    ! mm/month
    integer,          intent(out) :: year_idx(ntime), month_idx(ntime)

    integer :: t

    ! 1) Units
    precip_mm = precip_in * 1000.0_wp

    ! 2) Calendar from valid_time + units
    call build_calendar(valid_time, time_units, year_idx, month_idx)

    ! 3) Sanity: first/last, consecutive months
    call check_calendar(year_idx, month_idx)

    ! 4) Minimal echo so user sees alignment quickly
    write(*,'(A,1x,I0,1x,A,1x,I0)') 'Data start:', year_idx(1), '-', month_idx(1)
    write(*,'(A,1x,I0,1x,A,1x,I0)') 'Data end  :', year_idx(ntime), '-', month_idx(ntime)
  end subroutine preprocess_precip

  ! --------------------------------------------------------------
  ! Calendar decode: CF-like "X since YYYY-MM-DD" + int64 values
  ! --------------------------------------------------------------
  subroutine build_calendar(time_val, time_units, year, month)
    integer(int64), intent(in) :: time_val(:)
    character(len=*), intent(in) :: time_units
    integer, intent(out) :: year(size(time_val)), month(size(time_val))

    integer :: base_y, base_m, base_d, sec_per_unit, t
    integer(int64) :: secs
    call parse_units(time_units, base_y, base_m, base_d, sec_per_unit)
    do t = 1, size(time_val)
      secs = time_val(t) * int(sec_per_unit, int64)
      call epoch_seconds_to_ym(secs, base_y, base_m, base_d, year(t), month(t))
    end do
  end subroutine build_calendar

  subroutine check_calendar(year, month)
    integer, intent(in) :: year(:), month(:)
    integer :: t, gaps
    gaps = 0

    write(*,'(A,1x,I0)') 'Calendar length (months) =', size(year)
    write(*,'(A)') 'First 6 (y-m):'
    do t=1, min(6, size(year)); write(*,'(I4,"-",I2.2)') year(t), month(t); end do
    write(*,'(A)') 'Last 6 (y-m):'
    do t=max(1, size(year)-5), size(year); write(*,'(I4,"-",I2.2)') year(t), month(t); end do

    do t=2, size(year)
      if (.not. is_next_month(year(t-1), month(t-1), year(t), month(t))) gaps = gaps + 1
    end do
    if (gaps > 0) then
      write(*,'(A,I0)') 'FAIL(prep): calendar gaps found =', gaps   ! why: SPI windows rely on month alignment
      stop 1
    else
      write(*,'(A)') 'Calendar OK: consecutive months with no gaps.'
    end if
  end subroutine check_calendar

  ! ----------------- helpers -----------------

  subroutine parse_units(units, y0, m0, d0, sec_per_unit)
    character(len=*), intent(in)  :: units
    integer, intent(out)          :: y0, m0, d0, sec_per_unit
    integer :: p, ios
    y0=1970; m0=1; d0=1; sec_per_unit=1
    if (index(units,'hours since')>0)   sec_per_unit = 3600
    if (index(units,'minutes since')>0) sec_per_unit = 60
    if (index(units,'seconds since')>0) sec_per_unit = 1
    p = index(units, 'since')
    if (p>0) then
      read(units(p+6:), '(I4,1x,I2,1x,I2)', iostat=ios) y0, m0, d0
      if (ios/=0) then; y0=1970; m0=1; d0=1; end if   ! why: conservative fallback
    end if
  end subroutine parse_units

  subroutine epoch_seconds_to_ym(sec_from_base, yb, mb, db, y, m)
    integer(int64), intent(in) :: sec_from_base
    integer, intent(in)        :: yb, mb, db
    integer, intent(out)       :: y, m
    integer(int64) :: days_left
    integer :: dim(12), d, mm, yy
    dim=[31,28,31,30,31,30,31,31,30,31,30,31]   ! why: simple civil calendar; calendar attr is informational here
    yy=yb; mm=mb; d=db
    days_left = sec_from_base / 86400_int64
    do while (days_left>0_int64)
      dim(2) = merge(29,28, is_leap(yy))
      if (d < dim(mm)) then
        d=d+1
      else
        d=1; mm=mm+1
        if (mm>12) then; mm=1; yy=yy+1; end if
      end if
      days_left = days_left - 1_int64
    end do
    y=yy; m=mm
  contains
    pure logical function is_leap(yy)
      integer, intent(in) :: yy
      is_leap = (mod(yy,4)==0 .and. (mod(yy,100)/=0 .or. mod(yy,400)==0))
    end function is_leap
  end subroutine epoch_seconds_to_ym

  pure logical function is_next_month(y1,m1,y2,m2)
    integer, intent(in) :: y1,m1,y2,m2
    is_next_month = ( (y2==y1 .and. m2==m1+1) .or. (m1==12 .and. m2==1 .and. y2==y1+1) )
  end function is_next_month

end module prep_mod