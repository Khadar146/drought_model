program test1
  use netcdf
  implicit none

  integer, parameter :: t_dim = 540, lat_dim = 41, lon_dim = 66
  integer :: ncid, varid, retval, t
  real, allocatable :: tp(:,:,:), ts(:)
  character(len=*), parameter :: filename = "data/Climate/ERA5_Land/Precipitation/tp_ERA5Land_1981-2024.nc"

  ! Open NetCDF file
  retval = nf90_open(filename, NF90_NOWRITE, ncid)
  if (retval /= nf90_noerr) then
     print *, "Error opening file:", nf90_strerror(retval)
     stop
  end if

  ! Get variable ID
  retval = nf90_inq_varid(ncid, "tp", varid)
  if (retval /= nf90_noerr) then
     print *, "Error getting varid:", nf90_strerror(retval)
     stop
  end if

  ! Allocate and read full precipitation data
  allocate(tp(lon_dim, lat_dim, t_dim))
  retval = nf90_get_var(ncid, varid, tp)
  if (retval /= nf90_noerr) then
     print *, "Error reading variable:", nf90_strerror(retval)
     stop
  end if

  ! Close NetCDF file
  retval = nf90_close(ncid)
  if (retval /= nf90_noerr) then
     print *, "Error closing file:", nf90_strerror(retval)
     stop
  end if

  ! Allocate and extract full time series for grid point (10,35)
  allocate(ts(t_dim))
  do t = 1, t_dim
     ts(t) = tp(10,35,t) * 1000.0  ! Convert m to mm
  end do

  ! Write time series to CSV
  open(unit=10, file="output/test/precip_timeseries.csv", status="replace")
  write(10,*) "Month,Precip_mm"
  do t = 1, t_dim
     write(10,'(I4,1x,F10.2)') t, ts(t)
  end do
  close(10)

  print *, "✅ Saved precipitation time series for (10,35) to CSV."

end program test1