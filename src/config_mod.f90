! =====================================================================
! File: src/config_mod.f90
! Purpose: Centralized configuration for drought analysis pipeline
! Author: Khadar - Somaliland Drought Analysis
! =====================================================================
module config_mod
  use iso_fortran_env, only: wp => real64
  implicit none
  private
  public :: drought_config, load_config, print_config

  type :: drought_config
    ! Input data
    character(len=256) :: input_file = "data/historical/precip_monthly.nc"
    character(len=20)  :: var_name = "tp"
    character(len=256) :: output_dir = "output/spi/"
    
    ! Reference period for SPI calculation (critical for drought analysis)
    integer :: ref_start_year = 1991
    integer :: ref_end_year   = 2020
    
    ! SPI accumulation periods (months)
    integer :: num_periods = 4
    integer :: accum_periods(10) = [1, 3, 6, 12, 0, 0, 0, 0, 0, 0]
    
    ! Quality control
    real(wp) :: min_data_coverage = 0.70_wp    ! Require 70% valid data
    real(wp) :: spi_extreme_cap = 5.0_wp       ! Cap extreme SPI values
    
    ! Drought thresholds for Somaliland context
    real(wp) :: moderate_drought = -1.0_wp     ! SPI <= -1.0
    real(wp) :: severe_drought   = -1.5_wp     ! SPI <= -1.5
    real(wp) :: extreme_drought  = -2.0_wp     ! SPI <= -2.0
    
    ! Regional settings
    character(len=50) :: region = "Somaliland"
    character(len=50) :: analysis_type = "Historical Drought Analysis"
  end type drought_config

contains

  subroutine load_config(cfg, config_file)
    type(drought_config), intent(out) :: cfg
    character(len=*), intent(in), optional :: config_file
    
    ! Default configuration - can be overridden by file if needed
    cfg = drought_config()
    
    ! Future enhancement: read from config file if provided
    if (present(config_file)) then
      ! Could implement TOML/JSON reading here
      print *, "Loading configuration from: ", trim(config_file)
    end if
  end subroutine load_config

  subroutine print_config(cfg)
    type(drought_config), intent(in) :: cfg
    integer :: i
    
    print *, ""
    print *, "=== DROUGHT ANALYSIS CONFIGURATION ==="
    print *, "Region: ", trim(cfg%region)
    print *, "Analysis: ", trim(cfg%analysis_type)
    print *, "Input file: ", trim(cfg%input_file)
    print *, "Variable: ", trim(cfg%var_name)
    print *, "Output directory: ", trim(cfg%output_dir)
    print *, "Reference period: ", cfg%ref_start_year, "-", cfg%ref_end_year
    print *, "SPI periods (months):", (cfg%accum_periods(i), i=1,cfg%num_periods)
    print *, "Drought thresholds:"
    print *, "  Moderate: SPI <=", cfg%moderate_drought
    print *, "  Severe:   SPI <=", cfg%severe_drought  
    print *, "  Extreme:  SPI <=", cfg%extreme_drought
    print *, "======================================="
    print *, ""
  end subroutine print_config

end module config_mod
