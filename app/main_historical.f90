program main_historical
  use io_mod,    only: read_precip_from_netcdf
  use prep_mod,  only: preprocess_precip
  use spi_mod,   only: compute_spi
  use output_mod,only: save_spi_to_netcdf, save_spi_summary_csv
  implicit none

  ! Parameters
  integer, parameter :: nlon = 120, nlat = 80, ntime = 540 ! adjust as needed
  character(len=*), parameter :: filename = "data/historical/precip_monthly.nc"
  character(len=*), parameter :: varname = "tp"
  character(len=*), parameter :: outdir = "output/spi/"

  ! Data Arrays
  real(8) :: precip_raw(nlon, nlat, ntime)
  real(8) :: precip_mm(nlon, nlat, ntime)
  real(8) :: spi03(nlon, nlat, ntime)
  real(8) :: spi06(nlon, nlat, ntime)
  real(8) :: spi12(nlon, nlat, ntime)

  ! Step 1: Read and clean input
  call read_precip_from_netcdf(filename, varname, precip_raw, nlon, nlat, ntime)

  ! Step 2: Preprocessing (units, fill, aggregation)
  call preprocess_precip(precip_raw, precip_mm, nlon, nlat, ntime)

  ! Step 3: Compute SPI for 3, 6, 12 months
  call compute_spi(precip_mm, 3,  spi03, nlon, nlat, ntime)
  call compute_spi(precip_mm, 6,  spi06, nlon, nlat, ntime)
  call compute_spi(precip_mm, 12, spi12, nlon, nlat, ntime)

  ! Step 4: Save to NetCDF and CSV
  call save_spi_to_netcdf(spi03, nlon, nlat, ntime, outdir//'spi_03.nc', 'SPI_3')
  call save_spi_to_netcdf(spi06, nlon, nlat, ntime, outdir//'spi_06.nc', 'SPI_6')
  call save_spi_to_netcdf(spi12, nlon, nlat, ntime, outdir//'spi_12.nc', 'SPI_12')

  call save_spi_summary_csv(spi03, nlon, nlat, ntime, outdir//'spi_03_summary.csv')
  call save_spi_summary_csv(spi06, nlon, nlat, ntime, outdir//'spi_06_summary.csv')
  call save_spi_summary_csv(spi12, nlon, nlat, ntime, outdir//'spi_12_summary.csv')

  print *, "> Main historical SPI pipeline complete. Outputs saved to:", trim(outdir)

end program main_historical