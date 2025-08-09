program main_historical
  use iso_fortran_env, only: wp => real64, int64
  use io_mod,     only: read_precip_from_netcdf
  use prep_mod,   only: preprocess_precip
  use spi_mod,    only: compute_spi
  implicit none

  ! ---- I/O paths ----
  character(len=*), parameter :: filename = "data/historical/precip_monthly.nc"
  character(len=*), parameter :: varname  = "tp"
  character(len=*), parameter :: outdir   = "output/spi/"

  ! ---- Discovered sizes ----
  integer :: nlon, nlat, ntime

  ! ---- Time / reference period ----
  integer :: start_year, start_month        ! will derive from valid_time
  integer, parameter :: ref_start = 1991
  integer, parameter :: ref_end   = 2020

  ! ---- Accumulation periods ----
  integer, parameter :: num_periods = 3
  integer, parameter :: accum_periods(num_periods) = [1, 3, 12]

  ! ---- Data arrays ----
  real(wp), allocatable :: precip(:,:,:)
  real(wp), allocatable :: precip_mm(:,:,:)
  real(wp), allocatable :: spi_all(:,:,:,:)
  real(wp), allocatable :: tmp3d(:,:,:)

  ! ---- Time metadata from IO ----
  integer(int64),   allocatable :: valid_time(:)
  character(len=:), allocatable :: time_units, calendar

  ! ---- Calendar indices from PREP ----
  integer, allocatable :: year_idx(:), month_idx(:)

  integer :: p, t

  ! Ensure output directory exists
  call execute_command_line('mkdir -p ' // outdir)

  ! ---- Step 1: Read precip + time metadata from NetCDF ----
  call read_precip_from_netcdf(filename, varname,              &
                               precip, nlon, nlat, ntime,      &
                               valid_time, time_units, calendar)

  ! Allocate downstream arrays
  allocate(precip_mm(nlon, nlat, ntime))
  allocate(spi_all(nlon, nlat, ntime, num_periods))
  allocate(tmp3d(nlon, nlat, ntime))
  allocate(year_idx(ntime), month_idx(ntime))

! ---- Step 2: Preprocess (units, calendar from valid_time, fill, optional sums) ----
call preprocess_precip(precip, nlon, nlat, ntime,                     &
                       precip_mm,                                     &
                       year_idx, month_idx,                           &
                       valid_time, time_units, calendar)

! Derive start_year/start_month from the first timestamp for SPI alignment
start_year  = year_idx(1)
start_month = month_idx(1)

! --- Quick check ---
print *, "Data starts at year=", start_year, " month=", start_month
print *, "First 3 timestamps:"
do t = 1, min(3, ntime)
    print *, "  t=", t, " -> ", year_idx(t), "-", month_idx(t)
end do
print *, "Last 3 timestamps:"
do t = max(ntime-2, 1), ntime
    print *, "  t=", t, " -> ", year_idx(t), "-", month_idx(t)
end do

  ! ---- Step 3: SPI for [1,3,12] months (Gamma fit per calendar month over 1991–2020) ----
  call compute_spi(precip_mm, nlon, nlat, ntime,                        &
                   start_year, start_month,                             &
                   ref_start, ref_end,                                  &
                   accum_periods, num_periods,                          &
                   spi_all)

  ! ---- Step 4: Save each period to NetCDF and CSV ----
  do p = 1, num_periods
     tmp3d = spi_all(:,:,:,p)
     select case (accum_periods(p))
     case (1)
        call save_to_netcdf(tmp3d, nlon, nlat, ntime, outdir//'spi_01.nc', 'SPI_1')
        call save_summary_csv(tmp3d, nlon, nlat, ntime, outdir//'spi_01_summary.csv')
     case (3)
        call save_to_netcdf(tmp3d, nlon, nlat, ntime, outdir//'spi_03.nc', 'SPI_3')
        call save_summary_csv(tmp3d, nlon, nlat, ntime, outdir//'spi_03_summary.csv')
     case (12)
        call save_to_netcdf(tmp3d, nlon, nlat, ntime, outdir//'spi_12.nc', 'SPI_12')
        call save_summary_csv(tmp3d, nlon, nlat, ntime, outdir//'spi_12_summary.csv')
     end select
  end do

  print *, "> Main historical SPI pipeline complete. Outputs saved to:", trim(outdir)

contains

  ! ---------------------------
  ! Simple NetCDF writer (3D)
  ! ---------------------------
  subroutine save_to_netcdf(var3d, nlon, nlat, ntime, path, vname)
    use netcdf
    use iso_fortran_env, only: wp => real64
    implicit none
    integer,          intent(in) :: nlon, nlat, ntime
    real(wp),         intent(in) :: var3d(nlon, nlat, ntime)
    character(len=*), intent(in) :: path, vname
    integer :: ncid, lon_dim, lat_dim, time_dim, varid, rv
    integer :: dimids(3)

    print *, '>> Writing NetCDF: ', trim(path)

    rv = nf90_create(trim(path), NF90_CLOBBER, ncid);                 call chk(rv,'create')
    rv = nf90_def_dim(ncid, "lon",  nlon,  lon_dim);                  call chk(rv,'def_dim lon')
    rv = nf90_def_dim(ncid, "lat",  nlat,  lat_dim);                  call chk(rv,'def_dim lat')
    rv = nf90_def_dim(ncid, "time", ntime, time_dim);                 call chk(rv,'def_dim time')

    dimids = [lon_dim, lat_dim, time_dim]
    rv = nf90_def_var(ncid, trim(vname), NF90_DOUBLE, dimids, varid); call chk(rv,'def_var data')

    rv = nf90_enddef(ncid);                                           call chk(rv,'enddef')
    rv = nf90_put_var(ncid, varid, var3d);                            call chk(rv,'put data')
    rv = nf90_close(ncid);                                            call chk(rv,'close')
    print *, '>> Wrote: ', trim(path)
  end subroutine save_to_netcdf

  ! ---------------------------
  ! Tiny CSV summary (per grid)
  ! ---------------------------
  subroutine save_summary_csv(var3d, nlon, nlat, ntime, path)
    use iso_fortran_env, only: wp => real64
    implicit none
    integer,          intent(in) :: nlon, nlat, ntime
    real(wp),         intent(in) :: var3d(nlon, nlat, ntime)
    character(len=*), intent(in) :: path
    integer :: u, i, j, t, c
    real(wp) :: s, mn, mx

    open(newunit=u, file=path, status="replace", action="write")
    write(u,'(A)') "i,j,count,min,max,mean"
    do i = 1, nlon
      do j = 1, nlat
        s = 0.0_wp; c = 0; mn = huge(0.0_wp); mx = -huge(0.0_wp)
        do t = 1, ntime
          if (.not. (var3d(i,j,t) /= var3d(i,j,t))) then  ! NaN check via self-inequality
            s  = s + var3d(i,j,t)
            c  = c + 1
            mn = min(mn, var3d(i,j,t))
            mx = max(mx, var3d(i,j,t))
          end if
        end do
        if (c > 0) then
          write(u,'(I0,1x,I0,1x,I0,1x,ES12.5,1x,ES12.5,1x,ES12.5)') i, j, c, mn, mx, s/real(c,wp)
        else
          write(u,'(I0,1x,I0,1x,I0,1x,A,1x,A,1x,A)') i, j, 0, "NaN", "NaN", "NaN"
        end if
      end do
    end do
    close(u)
  end subroutine save_summary_csv

  subroutine chk(status, msg)
    use netcdf
    integer, intent(in) :: status
    character(len=*), intent(in) :: msg
    if (status /= nf90_noerr) then
      print *, "NetCDF error in ", trim(msg), ": ", trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine chk

end program main_historical