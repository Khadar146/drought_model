program test_netcdf
  use netcdf
  implicit none
  integer :: ncid, status

  status = nf90_open("dummy.nc", nf90_nowrite, ncid)
  print *, "NetCDF open status:", status
end program test_netcdf