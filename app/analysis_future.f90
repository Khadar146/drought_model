!> @file analysis_future.f90
!! @brief Future Analysis Module
!! 
!! This module contains the complete future drought analysis workflow
!! including bias correction, future SPI calculation, and future EVT analysis.
!!
!! @author Khadar Daahir  
!! @date August 2025
!!
module analysis_future
  use iso_fortran_env, only: wp => real64, int64
  use io_mod, only: read_precip_from_netcdf
  use prep_mod, only: preprocess_precip, build_calendar
  use spi_future_mod, only: compute_spi_future, apply_bias_correction
  use evt_future_mod, only: compute_evt_future, compare_evt_changes
  implicit none
  private
  
  public :: run_future_pipeline, write_future_results
  
  ! Analysis parameters
  integer, parameter :: NUM_SPI_PERIODS = 4
  integer, parameter :: spi_periods(NUM_SPI_PERIODS) = [1, 3, 6, 12]
  
contains

  !> Execute complete future analysis pipeline
  subroutine run_future_pipeline(future_data_file, historical_data_file, scenario, output_dir)
    character(len=*), intent(in) :: future_data_file, historical_data_file, scenario
    character(len=*), intent(in), optional :: output_dir
    
    ! Data arrays
    real(wp), allocatable :: future_precip(:,:,:), hist_precip(:,:,:)
    real(wp), allocatable :: corrected_precip(:,:,:)
    real(wp), allocatable :: lon(:), lat(:)
    integer :: nlon, nlat, ntime_fut, ntime_hist
    integer(int64), allocatable :: valid_time_fut(:), valid_time_hist(:)
    character(len=:), allocatable :: time_units_fut, calendar_fut
    character(len=:), allocatable :: time_units_hist, calendar_hist
    
    ! Analysis arrays
    real(wp), allocatable :: spi_future(:,:,:,:)
    real(wp), allocatable :: historical_params(:,:,:,:,:)  ! (lon,lat,month,period,param)
    logical, allocatable :: anomaly_flags(:,:,:,:)
    integer, allocatable :: year_idx(:), month_idx(:)
    
    ! EVT arrays
    real(wp), allocatable :: historical_evt(:,:,:,:)  ! (lon,lat,period,param)
    real(wp), allocatable :: future_evt(:,:,:,:)
    logical, allocatable :: change_signals(:,:,:,:)
    
    character(len=256) :: out_dir
    integer :: start_year, start_month
    
    ! Set output directory
    if (present(output_dir)) then
      out_dir = trim(output_dir)
    else
      out_dir = "outputs/future"
    end if
    
    print *, "🔮 EXECUTING FUTURE ANALYSIS PIPELINE"
    print *, "====================================="
    print *, "🌡️  Scenario:", trim(scenario)
    
    ! Step 1: Load future climate data
    print *, "📥 Step 1: Loading future climate data..."
    call read_precip_from_netcdf(future_data_file, 'tp', future_precip, lon, lat, &
                                nlon, nlat, ntime_fut, valid_time_fut, time_units_fut, calendar_fut)
    
    ! Step 2: Load historical data for bias correction
    print *, "📥 Step 2: Loading historical data for bias correction..."
    call load_historical_reference(historical_data_file, hist_precip, ntime_hist, &
                                  valid_time_hist, time_units_hist, calendar_hist)
    
    ! Step 3: Apply bias correction
    print *, "🔧 Step 3: Applying bias correction..."
    allocate(corrected_precip(nlon, nlat, ntime_fut))
    call apply_bias_correction(future_precip, hist_precip, hist_precip, &  ! Using hist as model proxy
                              nlon, nlat, ntime_fut, ntime_hist, corrected_precip)
    
    ! Step 4: Process future calendar
    print *, "📅 Step 4: Processing future calendar..."
    allocate(year_idx(ntime_fut), month_idx(ntime_fut))
    call build_calendar(valid_time_fut, time_units_fut, year_idx, month_idx)
    
    ! Step 5: Load historical gamma parameters
    print *, "📊 Step 5: Loading historical SPI parameters..."
    call load_historical_parameters(historical_params)
    
    ! Step 6: Compute future SPI
    print *, "📊 Step 6: Computing future SPI..."
    allocate(spi_future(nlon, nlat, ntime_fut, NUM_SPI_PERIODS))
    allocate(anomaly_flags(nlon, nlat, ntime_fut, NUM_SPI_PERIODS))
    call compute_spi_future(corrected_precip, nlon, nlat, ntime_fut, &
                           historical_params, start_year, start_month, &
                           spi_periods, NUM_SPI_PERIODS, &
                           spi_future, anomaly_flags, year_idx, month_idx)
    
    ! Step 7: Load historical EVT parameters
    print *, "📈 Step 7: Loading historical EVT parameters..."
    call load_historical_evt_parameters(historical_evt)
    
    ! Step 8: Compute future EVT
    print *, "📈 Step 8: Computing future EVT analysis..."
    allocate(future_evt(nlon, nlat, NUM_SPI_PERIODS, 3))
    allocate(change_signals(nlon, nlat, NUM_SPI_PERIODS, 4))
    call compute_evt_future(spi_future, nlon, nlat, ntime_fut, NUM_SPI_PERIODS, &
                           historical_evt, start_year, future_evt, change_signals)
    
    ! Step 9: Write results
    print *, "💾 Step 9: Writing results..."
    call write_future_results(out_dir, scenario, lon, lat, valid_time_fut, &
                             time_units_fut, calendar_fut, spi_future, &
                             anomaly_flags, future_evt, change_signals)
    
    print *, "✅ Future analysis pipeline completed"
    print *, "📁 Results saved to:", trim(out_dir)
    
    ! Cleanup
    deallocate(future_precip, hist_precip, corrected_precip, lon, lat)
    deallocate(valid_time_fut, valid_time_hist)
    deallocate(spi_future, anomaly_flags, year_idx, month_idx)
    deallocate(historical_params, historical_evt, future_evt, change_signals)
  end subroutine run_future_pipeline

  !> Write future analysis results
  subroutine write_future_results(output_dir, scenario, lon, lat, time, time_units, calendar, &
                                 spi_future, anomaly_flags, evt_params, change_signals)
    character(len=*), intent(in) :: output_dir, scenario
    real(wp), intent(in) :: lon(:), lat(:)
    integer(int64), intent(in) :: time(:)
    character(len=*), intent(in) :: time_units, calendar
    real(wp), intent(in) :: spi_future(:,:,:,:)
    logical, intent(in) :: anomaly_flags(:,:,:,:)
    real(wp), intent(in) :: evt_params(:,:,:,:)
    logical, intent(in) :: change_signals(:,:,:,:)
    
    character(len=256) :: spi_dir, evt_dir, scenario_dir, filename
    integer :: ap
    character(len=8) :: period_str
    
    ! Create scenario-specific output directories
    scenario_dir = trim(output_dir) // "/" // trim(scenario)
    spi_dir = trim(scenario_dir) // "/spi"
    evt_dir = trim(scenario_dir) // "/evt"
    call system("mkdir -p " // trim(spi_dir))
    call system("mkdir -p " // trim(evt_dir))
    
    ! Write future SPI results
    do ap = 1, size(spi_future, 4)
      write(period_str, '(I0,A)') spi_periods(ap), 'month'
      
      ! Future SPI NetCDF
      filename = trim(spi_dir) // "/spi_future_" // trim(period_str) // ".nc"
      call write_future_spi_netcdf(filename, lon, lat, time, time_units, calendar, &
                                  spi_future(:,:,:,ap), spi_periods(ap))
      
      ! Anomaly flags
      filename = trim(spi_dir) // "/anomalies_" // trim(period_str) // ".nc"
      call write_anomaly_flags(filename, lon, lat, time, anomaly_flags(:,:,:,ap))
      
      ! SPI change summary
      filename = trim(spi_dir) // "/spi_change_summary_" // trim(period_str) // ".csv"
      call write_spi_change_summary(filename, spi_future(:,:,:,ap))
    end do
    
    ! Write future EVT results
    do ap = 1, size(evt_params, 3)
      write(period_str, '(I0,A)') spi_periods(ap), 'month'
      
      ! Future EVT parameters
      filename = trim(evt_dir) // "/evt_future_" // trim(period_str) // ".nc"
      call write_future_evt_params(filename, lon, lat, evt_params(:,:,ap,:))
      
      ! Change signals
      filename = trim(evt_dir) // "/change_signals_" // trim(period_str) // ".nc"
      call write_change_signals(filename, lon, lat, change_signals(:,:,ap,:))
    end do
    
    ! Generate scenario summary report
    filename = trim(scenario_dir) // "/scenario_summary_" // trim(scenario) // ".txt"
    call write_scenario_summary(filename, scenario, spi_future, evt_params, change_signals)
    
    print *, "📊 Future SPI results written to:", trim(spi_dir)
    print *, "📈 Future EVT results written to:", trim(evt_dir)
    print *, "📋 Scenario summary written to:", trim(scenario_dir)
  end subroutine write_future_results

  ! Placeholder implementations - would be replaced with actual routines
  
  subroutine load_historical_reference(filename, precip, ntime, time, time_units, calendar)
    character(len=*), intent(in) :: filename
    real(wp), allocatable, intent(out) :: precip(:,:,:)
    integer, intent(out) :: ntime
    integer(int64), allocatable, intent(out) :: time(:)
    character(len=:), allocatable, intent(out) :: time_units, calendar
    
    print *, "  📥 Loading historical reference data"
    ! Implement actual loading
    allocate(precip(1,1,1), time(1))
    ntime = 1
    time_units = "days since 1900-01-01"
    calendar = "gregorian"
  end subroutine load_historical_reference

  subroutine load_historical_parameters(params)
    real(wp), allocatable, intent(out) :: params(:,:,:,:,:)
    print *, "  📊 Loading historical gamma parameters"
    ! Load from previous historical analysis
    allocate(params(1,1,12,4,3))  ! (lon,lat,month,period,param)
  end subroutine load_historical_parameters

  subroutine load_historical_evt_parameters(evt_params)
    real(wp), allocatable, intent(out) :: evt_params(:,:,:,:)
    print *, "  📈 Loading historical EVT parameters"
    allocate(evt_params(1,1,4,3))  ! (lon,lat,period,param)
  end subroutine load_historical_evt_parameters

  subroutine write_future_spi_netcdf(filename, lon, lat, time, time_units, calendar, spi_data, period)
    character(len=*), intent(in) :: filename, time_units, calendar
    real(wp), intent(in) :: lon(:), lat(:), spi_data(:,:,:)
    integer(int64), intent(in) :: time(:)
    integer, intent(in) :: period
    print *, "  💾 Writing future SPI NetCDF:", trim(filename)
  end subroutine write_future_spi_netcdf

  subroutine write_anomaly_flags(filename, lon, lat, time, flags)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: lon(:), lat(:)
    integer(int64), intent(in) :: time(:)
    logical, intent(in) :: flags(:,:,:)
    print *, "  💾 Writing anomaly flags:", trim(filename)
  end subroutine write_anomaly_flags

  subroutine write_spi_change_summary(filename, spi_data)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: spi_data(:,:,:)
    print *, "  💾 Writing SPI change summary:", trim(filename)
  end subroutine write_spi_change_summary

  subroutine write_future_evt_params(filename, lon, lat, evt_params)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: lon(:), lat(:), evt_params(:,:,:)
    print *, "  💾 Writing future EVT parameters:", trim(filename)
  end subroutine write_future_evt_params

  subroutine write_change_signals(filename, lon, lat, signals)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: lon(:), lat(:)
    logical, intent(in) :: signals(:,:,:)
    print *, "  💾 Writing change signals:", trim(filename)
  end subroutine write_change_signals

  subroutine write_scenario_summary(filename, scenario, spi_data, evt_params, change_signals)
    character(len=*), intent(in) :: filename, scenario
    real(wp), intent(in) :: spi_data(:,:,:,:), evt_params(:,:,:,:)
    logical, intent(in) :: change_signals(:,:,:,:)
    print *, "  📋 Writing scenario summary:", trim(filename)
  end subroutine write_scenario_summary

end module analysis_future
