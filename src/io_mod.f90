!> @file io_mod.f90
!! @brief Reads precipitation data from NetCDF into 3D arrays with cleaning and logs statistics.
!! @details Handles _FillValue and GRIB_missingValue replacement with NaN,
!!          computes basic statistics, flags anomalies, and saves them to a log file.
!!
!! @author Khadar
!! @date 2025

module io_mod
  use netcdf
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use iso_fortran_env, only: wp => real64
  implicit none
  private
  public :: read_precip_from_netcdf

contains

  subroutine read_precip_from_netcdf(filename, var_name, precip, nlon, nlat, ntime)
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: var_name
    integer,          intent(in)  :: nlon, nlat, ntime
    real(wp),         intent(out) :: precip(nlon, nlat, ntime)

    integer :: ncid, varid, rv
    real(wp) :: fill_value, grib_missing
    logical :: has_fill, has_grib
    integer :: missing_count, neg_count, extreme_count
    integer :: total_count
    logical, allocatable :: valid_mask(:,:,:)
    real(wp) :: min_val, max_val, mean_val
    character(len=*), parameter :: log_file = "io_stats.log"
    integer :: log_unit

    print *, "> Reading precipitation from NetCDF file:", trim(filename)

    ! --- Open NetCDF file ---
    rv = nf90_open(filename, NF90_NOWRITE, ncid)
    if (rv /= nf90_noerr) stop "Error opening file"

    ! --- Get variable ID ---
    rv = nf90_inq_varid(ncid, var_name, varid)
    if (rv /= nf90_noerr) stop "Error finding variable"

    ! --- Check for missing value attributes ---
    rv = nf90_get_att(ncid, varid, "_FillValue", fill_value)
    has_fill = (rv == nf90_noerr)

    rv = nf90_get_att(ncid, varid, "GRIB_missingValue", grib_missing)
    has_grib = (rv == nf90_noerr)

    ! --- Read precipitation data ---
    rv = nf90_get_var(ncid, varid, precip)
    if (rv /= nf90_noerr) stop "Error reading variable"

    ! --- Close file ---
    rv = nf90_close(ncid)

    ! --- Replace GRIB missing values with NaN ---
    if (has_grib) where (precip == grib_missing) precip = ieee_value(precip, ieee_quiet_nan)

    ! --- Replace _FillValue with NaN ---
    if (has_fill .and. .not. ieee_is_nan(fill_value)) then
      where (precip == fill_value) precip = ieee_value(precip, ieee_quiet_nan)
    end if

    ! --- Create validity mask ---
    allocate(valid_mask(nlon, nlat, ntime))
    valid_mask = .not. ieee_is_nan(precip)

    ! --- Count missing values ---
    missing_count = count(.not. valid_mask)
    total_count = size(precip)
    print *, "  Missing values found:", missing_count, &
             " (", 100.0_wp * missing_count / total_count, "%)"

    ! --- Compute basic stats ---
    if (any(valid_mask)) then
      min_val = minval(precip, mask=valid_mask)
      max_val = maxval(precip, mask=valid_mask)
      mean_val = sum(precip, mask=valid_mask) / count(valid_mask)
      print *, "  Min precipitation (m/month):", min_val
      print *, "  Max precipitation (m/month):", max_val
      print *, "  Mean precipitation (m/month):", mean_val
    else
      min_val = ieee_value(0.0_wp, ieee_quiet_nan)
      max_val = min_val
      mean_val = min_val
      print *, "  No valid precipitation values found."
    end if

    ! --- Check anomalies ---
    neg_count = count((precip < 0.0_wp) .and. valid_mask)
    extreme_count = count((precip > 2.0_wp) .and. valid_mask)
    print *, "  Negative values count:", neg_count
    print *, "  Extreme values (>2.0 m/month) count:", extreme_count

    if (neg_count > 0) then
      print *, "  WARNING: Negative precipitation values detected!"
    else if (extreme_count > 0) then
      print *, "  WARNING: Extreme precipitation values detected!"
    else
      print *, "  Data ranges appear within typical bounds."
    end if

    ! --- Write log file ---
    open(newunit=log_unit, file=log_file, status="replace", action="write")
    write(log_unit,*) "NetCDF file: ", trim(filename)
    write(log_unit,*) "Variable: ", trim(var_name)
    write(log_unit,*) "Dimensions: ", nlon, " x ", nlat, " x ", ntime
    write(log_unit,*) "Missing values: ", missing_count, " (", 100.0_wp * missing_count / total_count, "%)"
    write(log_unit,*) "Min (m/month): ", min_val
    write(log_unit,*) "Max (m/month): ", max_val
    write(log_unit,*) "Mean (m/month): ", mean_val
    write(log_unit,*) "Negative values: ", neg_count
    write(log_unit,*) "Extreme values (>2.0 m/month): ", extreme_count
    close(log_unit)

    ! --- Clean up ---
    deallocate(valid_mask)

    print *, "> Successfully read precipitation data."
    print *, "> Statistics written to:", log_file

  end subroutine read_precip_from_netcdf

end module io_mod