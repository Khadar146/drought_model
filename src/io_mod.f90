!> @file io_mod.f90
!! @brief Module for reading precipitation data from NetCDF files.
!! @details This I/O module loads 3D climate data fields (e.g. total precipitation)
!!          from NetCDF format into memory using high-precision real types.
!!          Includes metadata and data consistency checks (e.g. units, lat-lon values).
!!
!! @author Khadar Daahir
!! @affiliation University of Glasgow, Climate Dynamics Lab
!! @date 2025

module io_mod
  use netcdf
  use iso_fortran_env, only: wp => real64
  implicit none

  private
  public :: read_precip_from_netcdf

contains

  subroutine read_precip_from_netcdf(filename, precip, lon, lat, time, status)
    character(len=*), intent(in)         :: filename
    real(wp), allocatable, intent(out)   :: precip(:,:,:)
    real(wp), allocatable, intent(out)   :: lon(:), lat(:), time(:)
    integer, intent(out)                 :: status

    integer :: ncid
    integer :: varid_tp, varid_lon, varid_lat, varid_time
    integer :: dimid_lon, dimid_lat, dimid_time
    integer :: nlon, nlat, ntime
    integer :: retval, att_len
    character(len=100) :: units_str
    real(wp), allocatable :: lat_raw(:), precip_raw(:,:,:)

    ! Open NetCDF file
    retval = nf90_open(trim(filename), nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      status = retval
      return
    end if

    ! Dimension lengths
    call nf90_inq_dimid(ncid, "longitude", dimid_lon)
    call nf90_inq_dimlen(ncid, dimid_lon, nlon)

    call nf90_inq_dimid(ncid, "latitude", dimid_lat)
    call nf90_inq_dimlen(ncid, dimid_lat, nlat)

    call nf90_inq_dimid(ncid, "valid_time", dimid_time)
    call nf90_inq_dimlen(ncid, dimid_time, ntime)

    ! Allocate
    allocate(lon(nlon), lat_raw(nlat), time(ntime))
    allocate(precip_raw(nlon, nlat, ntime))

    ! Variable IDs
    call nf90_inq_varid(ncid, "longitude", varid_lon)
    call nf90_inq_varid(ncid, "latitude", varid_lat)
    call nf90_inq_varid(ncid, "valid_time", varid_time)
    call nf90_inq_varid(ncid, "tp", varid_tp)

    ! Read variables
    call nf90_get_var(ncid, varid_lon, lon)
    call nf90_get_var(ncid, varid_lat, lat_raw)
    call nf90_get_var(ncid, varid_time, time)
    call nf90_get_var(ncid, varid_tp, precip_raw)

    ! Attribute checks (latitude)
    call nf90_inquire_attribute(ncid, varid_lat, "units", len=att_len)
    call nf90_get_att(ncid, varid_lat, "units", units_str)
    if (trim(units_str(1:att_len)) /= "degrees_north") then
      print *, "Latitude units not as expected."
      status = -1
      call nf90_close(ncid)
      return
    end if

    ! Attribute checks (longitude)
    call nf90_inquire_attribute(ncid, varid_lon, "units", len=att_len)
    call nf90_get_att(ncid, varid_lon, "units", units_str)
    if (trim(units_str(1:att_len)) /= "degrees_east") then
      print *, "Longitude units not as expected."
      status = -1
      call nf90_close(ncid)
      return
    end if

    ! Close file
    call nf90_close(ncid)

    ! Flip and convert
    allocate(lat(nlat))
    allocate(precip(nlon, nlat, ntime))
    lat = lat_raw(nlat:1:-1)
    precip = precip_raw(:, nlat:1:-1, :)

    precip = precip * 1000.0_wp  ! meters to mm
    status = 0
  end subroutine read_precip_from_netcdf

end module io_mod