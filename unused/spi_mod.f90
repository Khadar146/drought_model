module spi_mod
  ! -------------------------------------------------------------------
  ! SPI Module – Drought Index Construction
  ! Description:
  !   Reads monthly precipitation data from NetCDF for Somaliland,
  !   preprocesses it, and computes SPI at multiple timescales.
  !
  ! Dependencies: netcdf
  ! -------------------------------------------------------------------
  use netcdf
  implicit none

  private
  public :: run_spi_module

contains

  ! -------------------------------------------------------------------
  subroutine check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, "NetCDF Error: ", nf90_strerror(status)
      stop 1
    end if
  end subroutine check

  ! -------------------------------------------------------------------
  subroutine run_spi_module()
    character(len=*), parameter :: file_path = "data/historical/precip_monthly.nc"
    character(len=*), parameter :: varname   = "tp"

    ! Declare arrays for input data
    real, allocatable :: precip(:,:,:)
    real, allocatable :: time(:)
    real, allocatable :: lat(:), lon(:)

    ! Step 1: Load monthly precipitation from NetCDF
    call read_precip_from_netcdf(file_path, varname, precip, time, lat, lon)

    ! Step 2: Preprocess precipitation (e.g. unit conversion, missing data)
    call preprocess_precip(precip)

    ! Step 3: Compute SPI at multiple timescales (e.g. 1, 3, 6 months)
    call calculate_spi(precip, time, lat, lon)

    print *, "SPI module completed successfully."
  end subroutine run_spi_module

  ! -------------------------------------------------------------------
  subroutine read_precip_from_netcdf(path, varname, precip, time, lat, lon)
    use netcdf  

    character(len=*), intent(in)  :: path, varname
    real, allocatable, intent(out):: precip(:,:,:), time(:), lat(:), lon(:)
    integer :: ncid, varid, dimid_time, dimid_lat, dimid_lon
    integer :: ntime, nlat, nlon
    integer :: retval

    ! Open NetCDF file
    call check(nf90_open(trim(path), NF90_NOWRITE, ncid))

    ! Get dimensions
    call check(nf90_inq_dimid(ncid, "valid_time", dimid_time))
    call check(nf90_inq_dimlen(ncid, dimid_time, ntime))
    call check(nf90_inq_dimid(ncid, "latitude", dimid_lat))
    call check(nf90_inq_dimlen(ncid, dimid_lat, nlat))
    call check(nf90_inq_dimid(ncid, "longitude", dimid_lon))
    call check(nf90_inq_dimlen(ncid, dimid_lon, nlon))

    ! Allocate arrays
    allocate(precip(ntime, nlat, nlon))
    allocate(time(ntime), lat(nlat), lon(nlon))

    ! Read precipitation data
    call check(nf90_inq_varid(ncid, varname, varid))
    call check(nf90_get_var(ncid, varid, precip))

    ! Read coordinates
    call check(nf90_inq_varid(ncid, "valid_time", varid))
    call check(nf90_get_var(ncid, varid, time))
    call check(nf90_inq_varid(ncid, "latitude", varid))
    call check(nf90_get_var(ncid, varid, lat))
    call check(nf90_inq_varid(ncid, "longitude", varid))
    call check(nf90_get_var(ncid, varid, lon))

    ! Close file
    call check(nf90_close(ncid))

    print *, "Loaded precipitation data from ", trim(path)
    print *, "Shape: ", ntime, "x", nlat, "x", nlon
  end subroutine read_precip_from_netcdf

  ! -------------------------------------------------------------------
  subroutine preprocess_precip(precip)
    ! Converts units and handles missing values
    real, intent(inout) :: precip(:,:,:)
    integer :: i, j, k
    real, parameter :: missing_value = -9999.0  ! Assumed missing value
    integer :: ntime, nlat, nlon

    ntime = size(precip, 1)
    nlat  = size(precip, 2)
    nlon  = size(precip, 3)

    do i = 1, ntime
      do j = 1, nlat
        do k = 1, nlon
          if (precip(i,j,k) < 0.0 .or. precip(i,j,k) == missing_value) then
            precip(i,j,k) = 0.0
          else
            precip(i,j,k) = precip(i,j,k) * 1000.0  ! Convert m -> mm if needed
          end if
        end do
      end do
    end do

    print *, "Preprocessing complete: unit conversion and missing values handled."
  end subroutine preprocess_precip

  ! -------------------------------------------------------------------
  subroutine calculate_spi(precip, time, lat, lon)
    ! Simple SPI-1 calculation (demo at 1 grid point: center)
    real, intent(in) :: precip(:,:,:), time(:), lat(:), lon(:)
    integer :: ntime, nlat, nlon, i, j, k
    real, allocatable :: series(:)
    real :: mean_val, std_val
    real, allocatable :: spi(:)

    ntime = size(precip, 1)
    nlat  = size(precip, 2)
    nlon  = size(precip, 3)

    ! Select grid point (center point for demo)
    j = nlat / 2
    k = nlon / 2

    ! Extract 1D time series at that grid point
    allocate(series(ntime))
    do i = 1, ntime
      series(i) = precip(i, j, k)
    end do

    ! Compute mean and standard deviation
    mean_val = sum(series) / ntime
    std_val  = sqrt(sum((series - mean_val)**2) / (ntime - 1))

    ! Standardize to get SPI-1
    allocate(spi(ntime))
    do i = 1, ntime
      if (std_val > 0.0) then
        spi(i) = (series(i) - mean_val) / std_val
      else
        spi(i) = 0.0
      end if
    end do

    ! Save SPI-1 and validate against known events
    call save_spi_output(spi, time, "results/spi1_center.csv")
    call validate_historical_droughts(spi, time)

    print *, "SPI-1 computed at center point (j=", j, ", k=", k, ")"
  end subroutine calculate_spi

  ! -------------------------------------------------------------------
  subroutine save_spi_output(spi, time, out_path)
    real, intent(in) :: spi(:), time(:)
    character(len=*), intent(in) :: out_path
    integer :: i, n, unit

    n = size(spi)
    unit = 20  ! Arbitrary unit number

    open(unit=unit, file=trim(out_path), status='replace')
    write(unit,*) "time,sci_spi1"
    do i = 1, n
      write(unit, '(F12.2,1x,F8.4)') time(i), spi(i)
    end do
    close(unit)

    print *, "Saved SPI-1 time series to ", trim(out_path)
  end subroutine save_spi_output

  ! -------------------------------------------------------------------
  subroutine validate_historical_droughts(spi, time)
    real, intent(in) :: spi(:), time(:)
    integer :: i, n
    real :: year

    n = size(spi)
    print *, "Historical Event Validation (SPI vs Known Drought Years)"
    print *, "-----------------------------------------------"
    print *, " Index   Year     SPI"

    do i = 1, n
      year = 1970.0 + time(i) / (365.25 * 24.0 * 3600.0)  ! convert from seconds since epoch
      if (abs(year - 2011.0) < 1.0 .or. &
          abs(year - 2016.5) < 1.0 .or. &
          abs(year - 2021.0) < 1.0 .or. &
          abs(year - 2022.5) < 1.0) then
        write(*,'(I5,2X,F7.2,2X,F8.3)') i, year, spi(i)
      end if
    end do

    print *, "If SPI < -1.0, conditions likely aligned with historical droughts."
  end subroutine validate_historical_droughts

end module spi_mod



