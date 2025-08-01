!> @file main.f90
!! @brief Main application for drought modelling in Somaliland
!! @details Loads climate data and runs the drought analysis pipeline
!!          including SPI, EVT, and future forecasting modules.
!!
!!          Currently active: NetCDF data loading using io_mod.
!!          SPI and EVT modules are scaffolded for later integration.
!!
!! @author Khadar Daahir
!! @affiliation University of Glasgow, Climate Dynamics Lab
!! @date July 2025

program main
  use io_mod
  implicit none

  character(len=100) :: filename
  real(wp), allocatable :: precip(:,:,:), lon(:), lat(:), time(:)
  integer :: status

  print *, "Reading NetCDF precipitation data..."
  filename = "data/historical/precip_monthly.nc" 

  call read_precip_from_netcdf(trim(filename), precip, lon, lat, time, status)

  if (status == 0) then
    print *, "File read successfully!"
    print *, "Dimensions:"
    print *, "  Longitude points:", size(lon)
    print *, "  Latitude points :", size(lat)
    print *, "  Time steps      :", size(time)
  else
    print *, "NetCDF read failed. Status code =", status
  end if
end program main