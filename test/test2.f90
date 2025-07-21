! app/test2.f90
program test2
  use, intrinsic :: iso_fortran_env, only: compiler_version
  use netcdf  ! This is the critical test
  implicit none
  
  print *, "----------------------------------------"
  print *, "NetCDF Module Verification Program"
  print *, "----------------------------------------"
  print *, "Compiler: ", compiler_version()
  print *, "NetCDF Module Path: /Users/kdc/miniconda3/envs/drought_env/include"
  print *, "NF90_NOERR constant value: ", NF90_NOERR
  print *, "----------------------------------------"
  print *, "SUCCESS: NetCDF module is properly accessible"
  print *, "----------------------------------------"
end program test2