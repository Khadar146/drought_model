program main
  use io_mod, only: read_precip_from_netcdf
  use iso_fortran_env, only: wp => real64
  implicit none

  ! --- Declare dimensions and data array ---
  integer, parameter :: nlat = 41, nlon = 66, ntime = 540
  real(wp)           :: precip(nlon, nlat, ntime)
  character(len=128) :: filename

  ! File location (adjust as needed)
  filename = "data/historical/precip_monthly.nc"

  ! Read precipitation variable ("tp")
  call read_precip_from_netcdf(filename, "tp", precip, nlon, nlat, ntime)

end program main