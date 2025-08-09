!> NetCDF Input/Output Module for Drought Analysis
!! 
!! This module provides robust NetCDF file reading capabilities specifically
!! designed for precipitation data analysis in climate research applications.
!! 
!! Key Features:
!! - Automatic dimension detection and reordering
!! - Comprehensive missing value handling
!! - Data quality validation and reporting
!! - Memory-efficient operations for large datasets
!!
!! @author Khadar Daahir
!! @date 2025-08-09
!! @version 1.0.0
module io_mod
  use iso_fortran_env, only: wp => real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use netcdf
  implicit none
  private
  public :: read_precip_from_netcdf

contains

  !> Read precipitation data from NetCDF files with comprehensive validation
  !!
  !! This subroutine reads precipitation data from CF-compliant NetCDF files,
  !! automatically detecting dimension order and coordinate variables. It handles
  !! multiple missing value conventions and performs extensive data quality checks.
  !!
  !! The routine automatically reorders data to the standard (longitude, latitude, time)
  !! format regardless of the input file's dimension order, making downstream
  !! processing more consistent and reliable.
  !!
  !! @param[in]     filename    Path to the NetCDF input file
  !! @param[in]     var_name    Name of the precipitation variable to read
  !! @param[out]    precip      Precipitation data array (lon,lat,time) in meters
  !! @param[out]    lon         Longitude coordinate array in degrees East
  !! @param[out]    lat         Latitude coordinate array in degrees North  
  !! @param[out]    nlon        Number of longitude points
  !! @param[out]    nlat        Number of latitude points
  !! @param[out]    ntime       Number of time points
  !! @param[out]    valid_time  Time coordinate values (seconds since epoch)
  !! @param[out]    time_units  Time units string from NetCDF file
  !! @param[out]    calendar    Calendar type string from NetCDF file
  !!
  !! @note All missing values (GRIB_missingValue, _FillValue) are converted to IEEE NaN
  !! @note Input data units are preserved (typically meters for ERA5 data)
  !! @note Comprehensive validation report is written to 'io_read.log'
  !!
  !! @warning Stops execution on NetCDF errors or dimension inconsistencies
  !!
  !! Example usage:
  !! ```fortran
  !! use io_mod, only: read_precip_from_netcdf
  !! real(wp), allocatable :: precip(:,:,:), lon(:), lat(:)
  !! integer :: nlon, nlat, ntime
  !! call read_precip_from_netcdf('data.nc', 'tp', precip, lon, lat, &
  !!                              nlon, nlat, ntime, valid_time, time_units, calendar)
  !! ```
  subroutine read_precip_from_netcdf(filename, var_name,                 &
                                     precip, lon, lat,                   &
                                     nlon, nlat, ntime,                  &
                                     valid_time, time_units, calendar)
    ! --- inputs ---
    character(len=*), intent(in)  :: filename, var_name
    ! --- outputs ---
    real(wp), allocatable, intent(out) :: precip(:,:,:)   ! (lon,lat,time) in file units (m)
    real(wp), allocatable, intent(out) :: lon(:), lat(:)
    integer,  intent(out) :: nlon, nlat, ntime
    integer(int64), allocatable, intent(out) :: valid_time(:)
    character(len=:),   allocatable, intent(out) :: time_units, calendar

    ! --- locals ---
    integer :: ncid, varid, ndims, dimids(3), i
    character(len=NF90_MAX_NAME) :: dname
    integer :: dlen(3)
    integer :: pos_time, pos_lat, pos_lon
    integer :: id_time, id_lon, id_lat
    real(wp), allocatable :: raw(:,:,:)
    real(wp) :: grib_missing, fill_r
    real     :: fill_f
    integer  :: rv, has_grib, has_fill_r, has_fill_f
    integer  :: logu
    integer  :: total_count, missing_count, neg_count
    real(wp) :: vmin, vmax, vmean, miss_pct

    print *, '>>> IO: Reading file=', trim(filename), ' var=', trim(var_name)

    ! open + locate var
    call chk(nf90_open(trim(filename), NF90_NOWRITE, ncid), 'open')
    call chk(nf90_inq_varid(ncid, trim(var_name), varid), 'inq_varid data')
    call chk(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), 'inquire_variable data')
    if (ndims /= 3) stop 'FAIL(io): expected 3D precip variable'

    ! map dimension positions by name
    pos_time = -1; pos_lat = -1; pos_lon = -1
    do i=1,3
      call chk(nf90_inquire_dimension(ncid, dimids(i), name=dname, len=dlen(i)), 'inquire_dimension')
      select case (trim(dname))
      case ('time','valid_time','time_counter','t'); pos_time = i
      case ('lat','latitude','y');                   pos_lat  = i
      case ('lon','longitude','x');                  pos_lon  = i
      end select
    end do
    if (pos_time<0 .or. pos_lat<0 .or. pos_lon<0) stop 'FAIL(io): could not resolve time/lat/lon dims'

    ntime = dlen(pos_time); nlat = dlen(pos_lat); nlon = dlen(pos_lon)

    ! coord var ids (try common names sequentially)
    id_time = -1
    rv = nf90_inq_varid(ncid,'valid_time',id_time)
    if (rv /= nf90_noerr) rv = nf90_inq_varid(ncid,'time',id_time)
    if (rv /= nf90_noerr) rv = nf90_inq_varid(ncid,'time_counter',id_time)
    if (rv /= nf90_noerr) rv = nf90_inq_varid(ncid,'t',id_time)

    id_lat  = -1
    rv = nf90_inq_varid(ncid,'lat',id_lat)
    if (rv /= nf90_noerr) rv = nf90_inq_varid(ncid,'latitude',id_lat)
    if (rv /= nf90_noerr) rv = nf90_inq_varid(ncid,'y',id_lat)

    id_lon  = -1
    rv = nf90_inq_varid(ncid,'lon',id_lon)
    if (rv /= nf90_noerr) rv = nf90_inq_varid(ncid,'longitude',id_lon)
    if (rv /= nf90_noerr) rv = nf90_inq_varid(ncid,'x',id_lon)

    if (id_time<0 .or. id_lat<0 .or. id_lon<0) stop 'FAIL(io): missing coordinate variable(s)'

    ! read coords and time attrs
    allocate(valid_time(ntime), lat(nlat), lon(nlon))
    call chk(nf90_get_var(ncid, id_time, valid_time), 'get time')
    call chk(nf90_get_var(ncid, id_lat,  lat),        'get lat')
    call chk(nf90_get_var(ncid, id_lon,  lon),        'get lon')
    call get_att_str(ncid, id_time, 'units',    time_units)
    call get_att_str(ncid, id_time, 'calendar', calendar)

    ! read data in file order
    allocate(raw(dlen(1), dlen(2), dlen(3)))
    call chk(nf90_get_var(ncid, varid, raw), 'get data')

    ! missing value attributes
    has_grib  = nf90_get_att(ncid, varid, 'GRIB_missingValue', grib_missing)
    has_fill_r= nf90_get_att(ncid, varid, '_FillValue',        fill_r)
    has_fill_f= nf90_get_att(ncid, varid, '_FillValue',        fill_f)

    call chk(nf90_close(ncid), 'close')

    ! reorder to (lon,lat,time)
    allocate(precip(nlon, nlat, ntime))
    call reorder_to_lonlatime(raw, pos_lon, pos_lat, pos_time, precip)
    deallocate(raw)

    ! map all known missing sentinels to NaN
    if (has_grib == nf90_noerr) then
      where (precip == grib_missing) precip = ieee_value(0.0_wp, ieee_quiet_nan)
    end if
    if (has_fill_r == nf90_noerr) then
      where (precip == fill_r) precip = ieee_value(0.0_wp, ieee_quiet_nan)
    else if (has_fill_f == nf90_noerr) then
      where (precip == real(fill_f, wp)) precip = ieee_value(0.0_wp, ieee_quiet_nan)
    end if

    ! -------- sanity + log --------
    open(newunit=logu, file='io_read.log', status='replace', action='write')
    write(logu,'(A)') 'IO sanity report'
    write(logu,'(A)') '----------------'
    call print_dim('time', ntime, logu)
    call print_dim('lat ', nlat,  logu)
    call print_dim('lon ', nlon,  logu)
    call print_minmax('lon', lon, logu)
    call print_minmax('lat', lat, logu)
    write(logu,'(A,1x,A)') 'time units:', trim(time_units)
    write(logu,'(A,1x,A)') 'calendar   :', trim(calendar)

    total_count   = size(precip)
    missing_count = total_count - count(isfinite_r(precip))
    neg_count     = count( (precip < 0.0_wp) .and. isfinite_r(precip) )
    miss_pct      = 100.0_wp * real(missing_count,wp) / real(total_count,wp)

    if (missing_count < total_count) then
      vmin  = minval(precip, mask=isfinite_r(precip))
      vmax  = maxval(precip, mask=isfinite_r(precip))
      vmean = sum(   precip, mask=isfinite_r(precip)) / &
              real(total_count - missing_count, wp)
    else
      vmin=0.0_wp; vmax=0.0_wp; vmean=0.0_wp
    end if

    ! short writes to avoid line truncation
    write(*,'(A,I0)')        'NaN count =', missing_count
    write(*,'(A,F6.2,A)')    'NaN percent =', miss_pct, '%'
    write(*,'(A,I0)')        'Negative values count =', neg_count
    write(*,'(A,1x,ES12.4,1x,ES12.4,1x,ES12.4)') 'Min/Max/Mean (m):', vmin, vmax, vmean
    write(*,'(A,1x,ES12.5,1x,ES12.5)') 'lon min/max:', minval(lon), maxval(lon)
    write(*,'(A,1x,ES12.5,1x,ES12.5)') 'lat min/max:', minval(lat), maxval(lat)
    write(*,'(A,1x,A)') 'time units:', trim(time_units)

    write(logu,'(A,I0)')     'NaN count =', missing_count
    write(logu,'(A,F6.2,A)') 'NaN percent =', miss_pct, '%'
    write(logu,'(A,I0)')     'Negative values count =', neg_count
    write(logu,'(A,1x,ES12.4,1x,ES12.4,1x,ES12.4)') 'Min/Max/Mean (m):', vmin, vmax, vmean
    close(logu)
    print *, '>>> IO OK. Sanity written to io_read.log'

  contains
    subroutine chk(status, where)
      integer, intent(in) :: status
      character(len=*), intent(in) :: where
      if (status /= nf90_noerr) then
        write(*,*) 'NetCDF error @', trim(where), ': ', trim(nf90_strerror(status))
        stop 1
      end if
    end subroutine chk

    subroutine get_att_str(ncid, varid, aname, out)
      integer, intent(in) :: ncid, varid
      character(len=*), intent(in) :: aname
      character(len=:), allocatable, intent(out) :: out
      character(len=256) :: buf
      if (nf90_get_att(ncid, varid, trim(aname), buf) == nf90_noerr) then
        out = trim(buf)
      else
        out = ''
      end if
    end subroutine get_att_str

    subroutine reorder_to_lonlatime(raw, plon, plat, ptime, outv)
      ! why: normalize order once for downstream math
      real(wp), intent(in)  :: raw(:,:,:)
      integer,  intent(in)  :: plon, plat, ptime
      real(wp), intent(out) :: outv(:,:,:)
      integer :: i,j,t, il,jl,tl
      do t=1,size(outv,3)
        do j=1,size(outv,2)
          do i=1,size(outv,1)
            il = pick_index(plon, i,j,t)
            jl = pick_index(plat, i,j,t)
            tl = pick_index(ptime, i,j,t)
            outv(i,j,t) = raw(il,jl,tl)
          end do
        end do
      end do
    end subroutine reorder_to_lonlatime

    pure integer function pick_index(which, ilon, jlat, ttime) result(ix)
      integer, intent(in) :: which, ilon, jlat, ttime
      select case (which)
      case (1); ix = ilon
      case (2); ix = jlat
      case (3); ix = ttime
      case default; ix = 1
      end select
    end function pick_index

    pure logical function is_monotonic_real(a, decreasing)
      real(wp), intent(in) :: a(:)
      logical,  intent(in) :: decreasing
      integer :: k
      is_monotonic_real = .true.
      do k=2,size(a)
        if (decreasing) then
          if (a(k) > a(k-1)) then; is_monotonic_real=.false.; return; end if
        else
          if (a(k) < a(k-1)) then; is_monotonic_real=.false.; return; end if
        end if
      end do
    end function is_monotonic_real

    pure logical function is_monotonic_int64(a)
      integer(int64), intent(in) :: a(:)
      integer :: k
      is_monotonic_int64 = .true.
      do k=2,size(a)
        if (a(k) <= a(k-1)) then; is_monotonic_int64=.false.; return; end if
      end do
    end function is_monotonic_int64

    elemental logical function isfinite_r(x)
      real(wp), intent(in) :: x
      isfinite_r = (x == x) .and. (abs(x) < huge(1.0_wp))
    end function isfinite_r

    subroutine print_dim(name, n, u)
      character(len=*), intent(in) :: name
      integer, intent(in) :: n, u
      write(*,'(A,1x,I0)')  trim(name)//' length =', n
      write(u,'(A,1x,I0)')  trim(name)//' length =', n
    end subroutine print_dim

    subroutine print_minmax(name, a, u)
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: a(:)
      integer, intent(in) :: u
      write(*,'(A,1x,ES12.5,1x,ES12.5)') trim(name)//' min/max =', minval(a), maxval(a)
      write(u,'(A,1x,ES12.5,1x,ES12.5)')  trim(name)//' min/max =', minval(a), maxval(a)
    end subroutine print_minmax

  end subroutine read_precip_from_netcdf
end module io_mod