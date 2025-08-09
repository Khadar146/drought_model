module io_mod
  use netcdf
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use iso_fortran_env, only: wp => real64, int64
  implicit none
  private
  public :: read_precip_from_netcdf
contains

  subroutine read_precip_from_netcdf(filename, var_name,                 &
                                     precip, nlon, nlat, ntime,          &
                                     valid_time, time_units, calendar)
    ! inputs
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: var_name
    ! outputs
    integer,          intent(out) :: nlon, nlat, ntime
    real(wp), allocatable, intent(out) :: precip(:,:,:)
    integer(int64),   allocatable, intent(out) :: valid_time(:)
    character(len=:), allocatable, intent(out) :: time_units, calendar

    ! handles / metadata
    integer :: ncid, varid, tvarid
    integer :: ndims, dimids(3), dimlen
    character(len=NF90_MAX_NAME) :: dimname, time_var_name
    integer :: dimlen_lon, dimlen_lat, dimlen_time
    integer :: i, rv
    integer :: start3(3), count3(3)     ! <-- moved up (must be declared before exec stmts)

    ! missing values
    real(wp) :: fill_value, grib_missing
    logical  :: has_fill, has_grib

    ! stats
    integer :: missing_count, neg_count, extreme_count, total_count
    logical, allocatable :: valid_mask(:,:,:)
    real(wp) :: min_val, max_val, mean_val
    integer :: log_unit
    character(len=*), parameter :: log_file = "io_stats.log"

    ! temp read buffer in FILE order (t, lat, lon)
    real(wp), allocatable :: raw(:,:,:)

    print *, "> Reading precipitation from NetCDF file:", trim(filename)

    ! open file
    call chk(nf90_open(trim(filename), NF90_NOWRITE, ncid), "open file")

    ! precip var id
    call chk(nf90_inq_varid(ncid, trim(var_name), varid), "inq_varid "//trim(var_name))

    ! inquire variable dims (order matters!)
    call chk(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), "inquire_variable "//trim(var_name))
    if (ndims /= 3) stop "Expected 3D precip variable"

    ! resolve sizes (recognize multiple names)
    dimlen_lon  = -1; dimlen_lat = -1; dimlen_time = -1
    do i = 1, ndims
      call chk(nf90_inquire_dimension(ncid, dimids(i), name=dimname, len=dimlen), "inquire_dimension")
      select case (trim(dimname))
      case ("lon","longitude","x")
        dimlen_lon = dimlen
      case ("lat","latitude","y")
        dimlen_lat = dimlen
      case ("time","valid_time","time_counter","t")
        dimlen_time = dimlen
      case default
        print *, "  WARNING: Unrecognized dimension name in var: ", trim(dimname)
      end select
    end do
    if (dimlen_lon<0 .or. dimlen_lat<0 .or. dimlen_time<0) stop "Could not resolve lon/lat/time sizes"

    nlon  = dimlen_lon
    nlat  = dimlen_lat
    ntime = dimlen_time

    ! find time variable
    if (nf90_inq_varid(ncid, "valid_time", tvarid) == nf90_noerr) then
      time_var_name = "valid_time"
    else if (nf90_inq_varid(ncid, "time", tvarid) == nf90_noerr) then
      time_var_name = "time"
    else if (nf90_inq_varid(ncid, "time_counter", tvarid) == nf90_noerr) then
      time_var_name = "time_counter"
    else
      stop "No time variable named 'valid_time'/'time'/'time_counter' found"
    end if

    ! read time coord and attributes
    allocate(valid_time(ntime))
    call chk(nf90_get_var(ncid, tvarid, valid_time), "get_var "//trim(time_var_name))
    block
      character(len=256) :: tu, cal
      tu = ""; cal = ""
      if (nf90_get_att(ncid, tvarid, "units", tu)     == nf90_noerr) then
        time_units = trim(tu)
      else
        time_units = ""
      end if
      if (nf90_get_att(ncid, tvarid, "calendar", cal) == nf90_noerr) then
        calendar = trim(cal)
      else
        calendar = ""  ! CF default is "gregorian"
      end if
    end block

    ! allocate raw (time,lat,lon)
    allocate(raw(ntime, nlat, nlon))

    ! missing-value attrs
    rv = nf90_get_att(ncid, varid, "_FillValue",        fill_value);   has_fill = (rv == nf90_noerr)
    rv = nf90_get_att(ncid, varid, "GRIB_missingValue", grib_missing); has_grib = (rv == nf90_noerr)

    ! build start/count arrays that match dimids(1:3)
    start3 = [1,1,1]
    do i = 1, ndims
      call chk(nf90_inquire_dimension(ncid, dimids(i), name=dimname, len=dimlen), "inquire_dimension")
      select case (trim(dimname))
      case ("time","valid_time","time_counter","t")
        count3(i) = ntime
      case ("lat","latitude","y")
        count3(i) = nlat
      case ("lon","longitude","x")
        count3(i) = nlon
      case default
        count3(i) = dimlen
      end select
    end do

    call chk(nf90_get_var(ncid, varid, raw, start=start3, count=count3), "get_var "//trim(var_name))

    ! close file
    call chk(nf90_close(ncid), "close file")

    ! reorder to (lon,lat,time)
    allocate(precip(nlon, nlat, ntime))
    do i = 1, ntime
      precip(:,:,i) = transpose(raw(i,:,:))
    end do
    deallocate(raw)

    ! apply missing masks
    if (has_grib) where (precip == grib_missing) precip = ieee_value(0.0_wp, ieee_quiet_nan)
    if (has_fill  .and. .not. ieee_is_nan(fill_value)) &
         where (precip == fill_value) precip = ieee_value(0.0_wp, ieee_quiet_nan)

    ! stats
    allocate(valid_mask(nlon, nlat, ntime))
    valid_mask = .not. ieee_is_nan(precip)

    missing_count = count(.not. valid_mask)
    total_count   = size(precip)
    print *, "  Missing values found:", missing_count, "(",  &
             100.0_wp * real(missing_count,wp) / real(total_count,wp), "%)"

    if (any(valid_mask)) then
      min_val  = minval(precip, mask=valid_mask)
      max_val  = maxval(precip, mask=valid_mask)
      mean_val = sum(precip, mask=valid_mask) / count(valid_mask)
      print *, "  Min (m/month):", min_val
      print *, "  Max (m/month):", max_val
      print *, "  Mean (m/month):", mean_val
    else
      print *, "  No valid precipitation values found."
    end if

    neg_count     = count((precip < 0.0_wp) .and. valid_mask)
    extreme_count = count((precip > 2.0_wp) .and. valid_mask)
    if (neg_count > 0)     print *, "  WARNING: Negative precipitation values detected!"
    if (extreme_count > 0) print *, "  WARNING: Extreme precipitation values detected!"

    ! log
    open(newunit=log_unit, file=log_file, status="replace", action="write")
    write(log_unit,*) "NetCDF file: ", trim(filename)
    write(log_unit,*) "Variable: ", trim(var_name)
    write(log_unit,*) "Dimensions (lon x lat x time): ", nlon, " x ", nlat, " x ", ntime
    write(log_unit,*) "Missing values: ", missing_count, " (",  &
                      100.0_wp * real(missing_count,wp) / real(total_count,wp), "%)"
    write(log_unit,*) "Negative values: ", neg_count
    write(log_unit,*) "Extreme values (>2.0 m/month): ", extreme_count
    write(log_unit,*) "time var: ", trim(time_var_name)
    write(log_unit,*) "time units: ", trim(time_units)
    write(log_unit,*) "calendar: ", trim(calendar)
    close(log_unit)

    deallocate(valid_mask)
    print *, "> Successfully read precipitation + time metadata."
    print *, "> Statistics written to:", log_file

  contains
    subroutine chk(status, msg)
      integer, intent(in) :: status
      character(len=*), intent(in) :: msg
      if (status /= nf90_noerr) then
        print *, "NetCDF error in ", trim(msg), ": ", trim(nf90_strerror(status))
        stop 1
      end if
    end subroutine chk
  end subroutine read_precip_from_netcdf

end module io_mod