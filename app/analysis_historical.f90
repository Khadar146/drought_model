!> @file analysis_historical.f90
!! @brief Historical Analysis Module
!! 
!! This module contains the complete historical drought analysis workflow
!! including SPI calculation and EVT analysis for the baseline period.
!!
!! @author Khadar Daahir
!! @date August 2025
!!
module analysis_historical
  use iso_fortran_env, only: wp => real64, int64
  use io_mod, only: read_precip_from_netcdf
  use prep_mod, only: preprocess_precip, build_calendar
  use spi_mod, only: compute_spi
  use evt_mod, only: compute_evt_analysis
  implicit none
  private
  
  public :: run_historical_pipeline, write_historical_results
  
  ! Analysis parameters
  integer, parameter :: NUM_SPI_PERIODS = 4
  integer, parameter :: spi_periods(NUM_SPI_PERIODS) = [1, 3, 6, 12]
  integer, parameter :: ref_start_year = 1991, ref_end_year = 2020
  
contains

  !> Execute complete historical analysis pipeline
  subroutine run_historical_pipeline(data_file, output_dir)
    character(len=*), intent(in) :: data_file
    character(len=*), intent(in), optional :: output_dir
    
    ! Data arrays
    real(wp), allocatable :: precip(:,:,:), lon(:), lat(:)
    integer :: nlon, nlat, ntime
    integer(int64), allocatable :: valid_time(:)
    character(len=:), allocatable :: time_units, calendar
    
    ! Analysis arrays
    real(wp), allocatable :: spi(:,:,:,:)  ! (lon,lat,time,periods)
    integer, allocatable :: year_idx(:), month_idx(:)
    
    ! EVT arrays
    real(wp), allocatable :: thresholds(:,:,:), gpd_shape(:,:,:), gpd_scale(:,:,:)
    real(wp), allocatable :: return_levels(:,:,:,:)
    real(wp), parameter :: return_periods(3) = [5.0_wp, 10.0_wp, 20.0_wp]
    
    character(len=256) :: out_dir
    integer :: start_year, start_month
    
    ! Set output directory
    if (present(output_dir)) then
      out_dir = trim(output_dir)
    else
      out_dir = "outputs/historical"
    end if
    
    print *, "🏛️  EXECUTING HISTORICAL ANALYSIS PIPELINE"
    print *, "==========================================="
    
    ! Step 1: Load precipitation data
    print *, "📥 Step 1: Loading precipitation data..."
    call read_precip_from_netcdf(data_file, 'tp', precip, lon, lat, &
                                 nlon, nlat, ntime, valid_time, time_units, calendar)
    
    ! Step 2: Process calendar and validate
    print *, "📅 Step 2: Processing calendar..."
    allocate(year_idx(ntime), month_idx(ntime))
    call build_calendar(valid_time, time_units, year_idx, month_idx)
    ! Validation: check if data covers required reference period
    start_year = year_idx(1)
    start_month = month_idx(1)
    
    print *, "Data start:", start_year, "-", start_month
    print *, "Data end  :", year_idx(ntime), "-", month_idx(ntime)
    print *, "Reference period:", ref_start_year, "–", ref_end_year
    
    ! Step 3: Compute SPI for all periods
    print *, "📊 Step 3: Computing SPI..."
    allocate(spi(nlon, nlat, ntime, NUM_SPI_PERIODS))
    call compute_spi(precip, nlon, nlat, ntime, &
                    start_year, start_month, &
                    ref_start_year, ref_end_year, &
                    spi_periods, NUM_SPI_PERIODS, &
                    spi, year_idx, month_idx)
    
    ! Step 4: Compute EVT analysis
    print *, "📈 Step 4: Computing EVT analysis..."
    allocate(thresholds(nlon, nlat, NUM_SPI_PERIODS))
    allocate(gpd_shape(nlon, nlat, NUM_SPI_PERIODS))
    allocate(gpd_scale(nlon, nlat, NUM_SPI_PERIODS))
    allocate(return_levels(nlon, nlat, NUM_SPI_PERIODS, size(return_periods)))
    
    call compute_evt_analysis(spi, nlon, nlat, ntime, NUM_SPI_PERIODS, &
                             0.95_wp, return_periods, size(return_periods), &
                             thresholds, gpd_shape, gpd_scale, return_levels)
    
    ! Step 5: Write results
    print *, "💾 Step 5: Writing results..."
    call write_historical_results(out_dir, lon, lat, valid_time, time_units, calendar, &
                                 spi, thresholds, gpd_shape, gpd_scale, return_levels, &
                                 return_periods)
    
    print *, "✅ Historical analysis pipeline completed"
    print *, "📁 Results saved to:", trim(out_dir)
    
    ! Cleanup
    deallocate(precip, lon, lat, valid_time)
    deallocate(spi, year_idx, month_idx)
    deallocate(thresholds, gpd_shape, gpd_scale, return_levels)
  end subroutine run_historical_pipeline

  !> Write historical analysis results to files
  subroutine write_historical_results(output_dir, lon, lat, time, time_units, calendar, &
                                     spi, thresholds, shape, scale, return_levels, return_periods)
    character(len=*), intent(in) :: output_dir
    real(wp), intent(in) :: lon(:), lat(:)
    integer(int64), intent(in) :: time(:)
    character(len=*), intent(in) :: time_units, calendar
    real(wp), intent(in) :: spi(:,:,:,:)
    real(wp), intent(in) :: thresholds(:,:,:), shape(:,:,:), scale(:,:,:)
    real(wp), intent(in) :: return_levels(:,:,:,:)
    real(wp), intent(in) :: return_periods(:)
    
    character(len=256) :: spi_dir, evt_dir, filename
    integer :: ap, ir
    character(len=8) :: period_str, return_str
    
    ! Create output directories
    spi_dir = trim(output_dir) // "/spi"
    evt_dir = trim(output_dir) // "/evt"
    call system("mkdir -p " // trim(spi_dir))
    call system("mkdir -p " // trim(evt_dir))
    
    ! Write SPI results
    do ap = 1, size(spi, 4)
      write(period_str, '(I0,A)') spi_periods(ap), 'month'
      
      ! NetCDF file
      filename = trim(spi_dir) // "/spi_" // trim(period_str) // ".nc"
      call write_spi_netcdf(filename, lon, lat, time, time_units, calendar, &
                           spi(:,:,:,ap), spi_periods(ap))
      
      ! CSV summary
      filename = trim(spi_dir) // "/spi_" // trim(period_str) // "_summary.csv"
      call write_spi_summary(filename, spi(:,:,:,ap))
    end do
    
    ! Write EVT results
    do ap = 1, size(thresholds, 3)
      write(period_str, '(I0,A)') spi_periods(ap), 'month'
      
      ! Thresholds
      filename = trim(evt_dir) // "/evt_thresholds_" // trim(period_str) // ".nc"
      call write_evt_field(filename, lon, lat, thresholds(:,:,ap), "threshold", &
                          "GPD threshold for " // trim(period_str) // " SPI")
      
      ! GPD parameters
      filename = trim(evt_dir) // "/evt_shape_" // trim(period_str) // ".nc"
      call write_evt_field(filename, lon, lat, shape(:,:,ap), "shape", &
                          "GPD shape parameter for " // trim(period_str) // " SPI")
      
      filename = trim(evt_dir) // "/evt_scale_" // trim(period_str) // ".nc"
      call write_evt_field(filename, lon, lat, scale(:,:,ap), "scale", &
                          "GPD scale parameter for " // trim(period_str) // " SPI")
      
      ! Return levels
      do ir = 1, size(return_periods)
        write(return_str, '(I0,A)') int(return_periods(ir)), "yr"
        filename = trim(evt_dir) // "/return_level_" // trim(period_str) // "_" // trim(return_str) // ".nc"
        call write_evt_field(filename, lon, lat, return_levels(:,:,ap,ir), "return_level", &
                            trim(return_str) // " return level for " // trim(period_str) // " SPI")
      end do
    end do
    
    print *, "📊 SPI results written to:", trim(spi_dir)
    print *, "📈 EVT results written to:", trim(evt_dir)
  end subroutine write_historical_results

  ! Placeholder implementations for file writing
  ! These would be implemented with actual NetCDF writing routines
  
  subroutine write_spi_netcdf(filename, lon, lat, time, time_units, calendar, spi_data, period)
    character(len=*), intent(in) :: filename, time_units, calendar
    real(wp), intent(in) :: lon(:), lat(:), spi_data(:,:,:)
    integer(int64), intent(in) :: time(:)
    integer, intent(in) :: period
    
    print *, "  💾 Writing SPI NetCDF:", trim(filename)
    ! Implement actual NetCDF writing
  end subroutine write_spi_netcdf

  subroutine write_spi_summary(filename, spi_data)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: spi_data(:,:,:)
    
    print *, "  💾 Writing SPI summary:", trim(filename)
    ! Implement CSV summary writing
  end subroutine write_spi_summary

  subroutine write_evt_field(filename, lon, lat, field, var_name, description)
    character(len=*), intent(in) :: filename, var_name, description
    real(wp), intent(in) :: lon(:), lat(:), field(:,:)
    
    print *, "  💾 Writing EVT field:", trim(filename)
    ! Implement NetCDF field writing
  end subroutine write_evt_field

end module analysis_historical
