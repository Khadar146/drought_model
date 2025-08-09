!> @file main_drought_analysis.f90
!! @brief Main Drought Analysis Application
!! 
!! Unified application for Somaliland drought analysis supporting both
!! historical baseline establishment and future climate projections.
!! This is the single entry point for all drought analysis workflows.
!!
!! @author Khadar Daahir
!! @date August 2025
!!
program main_drought_analysis
  use analysis_historical, only: run_historical_pipeline
  use analysis_future, only: run_future_pipeline  
  use analysis_comparison, only: run_comparison_pipeline
  implicit none
  
  ! Command line arguments
  character(len=256) :: analysis_type, data_file, scenario
  integer :: num_args, status
  
  ! Print header
  call print_header()
  
  ! Parse command line arguments
  call parse_arguments(analysis_type, data_file, scenario, status)
  if (status /= 0) then
    call print_usage()
    stop 1
  end if
  
  ! Route to appropriate analysis
  select case (trim(adjustl(analysis_type)))
  case ('historical', 'hist', 'h')
    call run_historical_pipeline(data_file)
  case ('future', 'fut', 'f') 
    call run_future_analysis(data_file, scenario)
  case ('comparison', 'compare', 'c')
    call run_comparison_analysis(data_file, scenario)
  case ('all', 'complete')
    call run_complete_analysis(data_file, scenario)
  case default
    print *, "ERROR: Unknown analysis type '", trim(analysis_type), "'"
    call print_usage()
    stop 1
  end select
  
  print *, ""
  print *, "🎉 Drought analysis completed successfully!"
  print *, "📊 Check outputs/ directory for results"

contains

  !> Print application header
  subroutine print_header()
    print *, "================================================================"
    print *, "    🌍 SOMALILAND DROUGHT ANALYSIS SYSTEM"
    print *, "================================================================"
    print *, "Research Hypothesis: 'Future droughts will become more extreme"
    print *, "                     and frequent due to climate change'"
    print *, ""
    print *, "Author: Khadar Daahir (University of Glasgow)"
    print *, "Date:   August 2025"
    print *, "================================================================"
    print *, ""
  end subroutine print_header

  !> Parse command line arguments
  subroutine parse_arguments(analysis_type, data_file, scenario, status)
    character(len=*), intent(out) :: analysis_type, data_file, scenario
    integer, intent(out) :: status
    
    status = 0
    analysis_type = ""
    data_file = ""
    scenario = ""
    
    num_args = command_argument_count()
    
    if (num_args < 1) then
      status = 1
      return
    end if
    
    ! First argument: analysis type
    call get_command_argument(1, analysis_type)
    
    ! Second argument: data file (optional, has defaults)
    if (num_args >= 2) then
      call get_command_argument(2, data_file)
    else
      ! Set default data file based on analysis type
      select case (trim(adjustl(analysis_type)))
      case ('historical', 'hist', 'h')
        data_file = "data/climate/historical/precip_monthly.nc"
      case ('future', 'fut', 'f')
        data_file = "data/climate/future/precip_future.nc"
      case default
        data_file = "data/climate/historical/precip_monthly.nc"
      end select
    end if
    
    ! Third argument: scenario (for future analysis)
    if (num_args >= 3) then
      call get_command_argument(3, scenario)
    else
      scenario = "ssp585"  ! Default scenario
    end if
  end subroutine parse_arguments

  !> Print usage information
  subroutine print_usage()
    print *, "USAGE:"
    print *, "  fpm run DroughtAnalysis -- <analysis_type> [data_file] [scenario]"
    print *, ""
    print *, "ANALYSIS TYPES:"
    print *, "  historical (h)  - Historical baseline analysis (1980-2024)"
    print *, "  future (f)      - Future projections analysis (2025-2100)"
    print *, "  comparison (c)  - Compare historical vs future"
    print *, "  all            - Complete analysis workflow"
    print *, ""
    print *, "EXAMPLES:"
    print *, "  fpm run DroughtAnalysis -- historical"
    print *, "  fpm run DroughtAnalysis -- future data/climate/future/ssp585.nc"
    print *, "  fpm run DroughtAnalysis -- comparison"
    print *, "  fpm run DroughtAnalysis -- all"
    print *, ""
    print *, "OUTPUTS:"
    print *, "  📊 outputs/historical/ - Historical SPI and EVT results"
    print *, "  📈 outputs/future/     - Future projection results"
    print *, "  📉 outputs/comparison/ - Change analysis results"
  end subroutine print_usage

  !> Run future analysis workflow  
  subroutine run_future_analysis(data_file, scenario)
    character(len=*), intent(in) :: data_file, scenario
    logical :: file_exists
    
    print *, "🔮 FUTURE ANALYSIS WORKFLOW"
    print *, "==========================="
    print *, "📂 Data file: ", trim(data_file)
    print *, "🌡️  Scenario: ", trim(scenario)
    print *, "📅 Projection period: 2025-2100"
    print *, ""
    
    ! Check if future data file exists
    inquire(file=trim(data_file), exist=file_exists)
    if (.not. file_exists) then
      print *, "❌ ERROR: Future data file not found: ", trim(data_file)
      print *, ""
      print *, "📥 To get SSP climate projections:"
      print *, "   1. Visit: https://esgf-node.llnl.gov/search/cmip6/"
      print *, "   2. Search for: precipitation, ", trim(scenario), ", monthly"
      print *, "   3. Download data for Somaliland region (5-12°N, 42-51°E)"
      print *, "   4. Place in: data/climate/future/"
      print *, ""
      print *, "📋 Required SSP scenarios:"
      print *, "   - SSP1-2.6: Low emissions (Paris Agreement)"
      print *, "   - SSP2-4.5: Medium emissions (likely outcome)"  
      print *, "   - SSP3-7.0: High emissions (regional rivalry)"
      print *, "   - SSP5-8.5: Very high emissions (fossil development)"
      stop 1
    end if
    
    ! Call the future analysis module
    call run_future_pipeline(data_file, "data/climate/historical/precip_monthly.nc", &
                            scenario, "outputs/future")
    
    print *, "✅ Future analysis completed"
    print *, "📁 Results saved to: outputs/future/"
  end subroutine run_future_analysis

  !> Run comparison analysis
  subroutine run_comparison_analysis(data_file, scenario)
    character(len=*), intent(in) :: data_file, scenario
    logical :: hist_exists, fut_exists
    
    print *, "📊 COMPARISON ANALYSIS WORKFLOW"
    print *, "==============================="
    print *, "🔬 Testing hypothesis: Future droughts more extreme/frequent"
    print *, ""
    
    ! Check if both historical and future results exist
    inquire(file="outputs/historical/spi_1month.nc", exist=hist_exists)
    inquire(file="outputs/future/spi_1month.nc", exist=fut_exists)
    
    if (.not. hist_exists) then
      print *, "❌ ERROR: Historical results not found"
      print *, "   Please run: fpm run DroughtAnalysis -- historical"
      stop 1
    end if
    
    if (.not. fut_exists) then
      print *, "❌ ERROR: Future results not found"
      print *, "   Please run: fpm run DroughtAnalysis -- future"
      stop 1
    end if
    
    ! Call the comparison analysis module
    call run_comparison_pipeline("outputs/historical", "outputs/future", &
                                "outputs/comparison")
    
    print *, "✅ Comparison analysis completed"
    print *, "📁 Results saved to: outputs/comparison/"
  end subroutine run_comparison_analysis

  !> Run complete analysis workflow
  subroutine run_complete_analysis(data_file, scenario)
    character(len=*), intent(in) :: data_file, scenario
    
    print *, "🚀 COMPLETE ANALYSIS WORKFLOW"
    print *, "============================="
    print *, "🔬 Full hypothesis testing pipeline"
    print *, ""
    
    ! Run historical analysis first
    print *, "Step 1: Historical Analysis"
    print *, "---------------------------"
    call run_historical_pipeline(data_file, "outputs/historical")
    
    print *, ""
    print *, "Step 2: Future Analysis"
    print *, "-----------------------"
    call run_future_analysis(data_file, scenario)
    
    print *, ""
    print *, "Step 3: Comparison Analysis"
    print *, "---------------------------"
    call run_comparison_analysis(data_file, scenario)
    
    print *, "✅ Complete analysis workflow finished"
    print *, "🎯 Hypothesis testing results available in outputs/"
  end subroutine run_complete_analysis

end program main_drought_analysis
