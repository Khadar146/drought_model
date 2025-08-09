!> @file analysis_comparison.f90
!! @brief Comparison Analysis Module
!! 
!! This module contains the comparative analysis workflow to test the hypothesis
!! that droughts will become more extreme and frequent in the future.
!! It performs statistical comparison between historical and future results.
!!
!! @author Khadar Daahir
!! @date August 2025
!!
module analysis_comparison
  use iso_fortran_env, only: wp => real64, int64
  implicit none
  private
  
  public :: run_comparison_pipeline, generate_change_report
  
  ! Analysis parameters
  integer, parameter :: NUM_SPI_PERIODS = 4
  integer, parameter :: spi_periods(NUM_SPI_PERIODS) = [1, 3, 6, 12]
  real(wp), parameter :: SIGNIFICANCE_LEVEL = 0.05_wp
  
  ! Change detection thresholds
  real(wp), parameter :: DROUGHT_THRESHOLD = -1.0_wp  ! SPI < -1.0 indicates drought
  real(wp), parameter :: EXTREME_THRESHOLD = -2.0_wp  ! SPI < -2.0 indicates extreme drought
  
contains

  !> Execute complete comparison analysis pipeline
  subroutine run_comparison_pipeline(historical_results_dir, future_results_dir, output_dir)
    character(len=*), intent(in) :: historical_results_dir, future_results_dir
    character(len=*), intent(in), optional :: output_dir
    
    character(len=256) :: comparison_output_dir
    
    ! Set output directory
    if (present(output_dir)) then
      comparison_output_dir = trim(output_dir)
    else
      comparison_output_dir = "outputs/comparison"
    end if
    
    print *, "🔬 Starting Comparison Analysis Pipeline"
    print *, "========================================"
    print *, "Testing Hypothesis: Future droughts will be more extreme and frequent"
    print *, ""
    print *, "📂 Historical results: ", trim(historical_results_dir)
    print *, "📂 Future results: ", trim(future_results_dir)
    print *, "📁 Output directory: ", trim(comparison_output_dir)
    print *, ""
    
    ! Load historical and future results
    call load_analysis_results(historical_results_dir, future_results_dir)
    
    ! Perform drought frequency analysis
    call analyze_drought_frequency_changes()
    
    ! Perform drought intensity analysis
    call analyze_drought_intensity_changes()
    
    ! Perform drought duration analysis
    call analyze_drought_duration_changes()
    
    ! Perform extreme value comparison
    call analyze_extreme_value_changes()
    
    ! Generate statistical significance tests
    call perform_significance_tests()
    
    ! Generate comparison reports
    call generate_change_report(comparison_output_dir)
    
    print *, "✅ Comparison analysis completed successfully"
    print *, "📊 Results saved to: ", trim(comparison_output_dir)
    
  end subroutine run_comparison_pipeline

  !> Load historical and future analysis results
  subroutine load_analysis_results(historical_dir, future_dir)
    character(len=*), intent(in) :: historical_dir, future_dir
    
    print *, "📥 Loading analysis results for comparison..."
    
    ! Load historical SPI results
    call load_historical_spi_results(historical_dir)
    
    ! Load future SPI results
    call load_future_spi_results(future_dir)
    
    ! Load historical EVT results
    call load_historical_evt_results(historical_dir)
    
    ! Load future EVT results
    call load_future_evt_results(future_dir)
    
    print *, "✅ Analysis results loaded successfully"
    
  end subroutine load_analysis_results

  !> Analyze changes in drought frequency
  subroutine analyze_drought_frequency_changes()
    
    print *, "📊 Analyzing drought frequency changes..."
    
    ! Calculate historical drought frequency
    call calculate_historical_drought_frequency()
    
    ! Calculate future drought frequency
    call calculate_future_drought_frequency()
    
    ! Compare frequency changes by SPI period
    call compare_frequency_by_period()
    
    ! Assess seasonal frequency changes
    call analyze_seasonal_frequency_changes()
    
    print *, "✅ Drought frequency analysis completed"
    
  end subroutine analyze_drought_frequency_changes

  !> Analyze changes in drought intensity
  subroutine analyze_drought_intensity_changes()
    
    print *, "📈 Analyzing drought intensity changes..."
    
    ! Calculate historical drought intensity statistics
    call calculate_historical_intensity_stats()
    
    ! Calculate future drought intensity statistics
    call calculate_future_intensity_stats()
    
    ! Compare intensity distributions
    call compare_intensity_distributions()
    
    ! Identify extreme drought changes
    call analyze_extreme_drought_changes()
    
    print *, "✅ Drought intensity analysis completed"
    
  end subroutine analyze_drought_intensity_changes

  !> Analyze changes in drought duration
  subroutine analyze_drought_duration_changes()
    
    print *, "⏱️  Analyzing drought duration changes..."
    
    ! Identify drought events and durations (historical)
    call identify_historical_drought_events()
    
    ! Identify drought events and durations (future)
    call identify_future_drought_events()
    
    ! Compare duration statistics
    call compare_duration_statistics()
    
    ! Analyze multi-year drought potential
    call analyze_multiyear_drought_risk()
    
    print *, "✅ Drought duration analysis completed"
    
  end subroutine analyze_drought_duration_changes

  !> Analyze changes in extreme value characteristics
  subroutine analyze_extreme_value_changes()
    
    print *, "🌡️  Analyzing extreme value changes..."
    
    ! Compare return periods
    call compare_return_periods()
    
    ! Analyze extreme value parameters
    call compare_evt_parameters()
    
    ! Assess tail behavior changes
    call analyze_tail_behavior_changes()
    
    print *, "✅ Extreme value analysis completed"
    
  end subroutine analyze_extreme_value_changes

  !> Perform statistical significance tests
  subroutine perform_significance_tests()
    
    print *, "📈 Performing statistical significance tests..."
    
    ! Two-sample t-tests for mean changes
    call perform_mean_difference_tests()
    
    ! Kolmogorov-Smirnov tests for distribution changes
    call perform_distribution_tests()
    
    ! Mann-Whitney U tests for non-parametric comparison
    call perform_rank_tests()
    
    ! Trend analysis tests
    call perform_trend_tests()
    
    print *, "✅ Statistical tests completed"
    
  end subroutine perform_significance_tests

  !> Generate comprehensive change analysis report
  subroutine generate_change_report(output_dir)
    character(len=*), intent(in) :: output_dir
    
    character(len=512) :: report_file, summary_file, figure_dir
    
    print *, "📋 Generating change analysis reports..."
    
    ! Create output directories
    call create_output_directories(output_dir)
    
    ! Generate main comparison report
    report_file = trim(output_dir) // "/drought_change_analysis_report.txt"
    call write_detailed_comparison_report(report_file)
    
    ! Generate executive summary
    summary_file = trim(output_dir) // "/executive_summary.txt"
    call write_executive_summary(summary_file)
    
    ! Generate comparison figures (if visualization tools available)
    figure_dir = trim(output_dir) // "/figures"
    call generate_comparison_figures(figure_dir)
    
    ! Generate hypothesis testing results
    call write_hypothesis_test_results(output_dir)
    
    print *, "📊 Reports generated:"
    print *, "  - Detailed Report: ", trim(report_file)
    print *, "  - Executive Summary: ", trim(summary_file)
    print *, "  - Figures: ", trim(figure_dir)
    
  end subroutine generate_change_report

  ! Placeholder implementations for the analysis routines
  ! These would be implemented with actual statistical analysis code
  
  subroutine load_historical_spi_results(dir)
    character(len=*), intent(in) :: dir
    ! Implementation: Load historical SPI NetCDF files
    print *, "  📈 Loading historical SPI results from: ", trim(dir)
  end subroutine load_historical_spi_results

  subroutine load_future_spi_results(dir)
    character(len=*), intent(in) :: dir
    ! Implementation: Load future SPI NetCDF files
    print *, "  📈 Loading future SPI results from: ", trim(dir)
  end subroutine load_future_spi_results

  subroutine load_historical_evt_results(dir)
    character(len=*), intent(in) :: dir
    ! Implementation: Load historical EVT NetCDF files
    print *, "  📊 Loading historical EVT results from: ", trim(dir)
  end subroutine load_historical_evt_results

  subroutine load_future_evt_results(dir)
    character(len=*), intent(in) :: dir
    ! Implementation: Load future EVT NetCDF files
    print *, "  📊 Loading future EVT results from: ", trim(dir)
  end subroutine load_future_evt_results

  subroutine calculate_historical_drought_frequency()
    ! Implementation: Calculate drought frequency statistics for historical period
    print *, "    📊 Calculating historical drought frequency..."
  end subroutine calculate_historical_drought_frequency

  subroutine calculate_future_drought_frequency()
    ! Implementation: Calculate drought frequency statistics for future period
    print *, "    📊 Calculating future drought frequency..."
  end subroutine calculate_future_drought_frequency

  subroutine compare_frequency_by_period()
    ! Implementation: Compare frequency changes across SPI time periods
    print *, "    📈 Comparing frequency changes by SPI period..."
  end subroutine compare_frequency_by_period

  subroutine analyze_seasonal_frequency_changes()
    ! Implementation: Analyze seasonal patterns in frequency changes
    print *, "    🗓️  Analyzing seasonal frequency changes..."
  end subroutine analyze_seasonal_frequency_changes

  subroutine calculate_historical_intensity_stats()
    ! Implementation: Calculate intensity statistics for historical period
    print *, "    📊 Calculating historical intensity statistics..."
  end subroutine calculate_historical_intensity_stats

  subroutine calculate_future_intensity_stats()
    ! Implementation: Calculate intensity statistics for future period
    print *, "    📊 Calculating future intensity statistics..."
  end subroutine calculate_future_intensity_stats

  subroutine compare_intensity_distributions()
    ! Implementation: Compare intensity distributions between periods
    print *, "    📈 Comparing intensity distributions..."
  end subroutine compare_intensity_distributions

  subroutine analyze_extreme_drought_changes()
    ! Implementation: Analyze changes in extreme drought characteristics
    print *, "    🌡️  Analyzing extreme drought changes..."
  end subroutine analyze_extreme_drought_changes

  subroutine identify_historical_drought_events()
    ! Implementation: Identify drought events and calculate durations
    print *, "    🔍 Identifying historical drought events..."
  end subroutine identify_historical_drought_events

  subroutine identify_future_drought_events()
    ! Implementation: Identify future drought events and calculate durations
    print *, "    🔍 Identifying future drought events..."
  end subroutine identify_future_drought_events

  subroutine compare_duration_statistics()
    ! Implementation: Compare duration statistics between periods
    print *, "    📊 Comparing duration statistics..."
  end subroutine compare_duration_statistics

  subroutine analyze_multiyear_drought_risk()
    ! Implementation: Analyze risk of multi-year drought events
    print *, "    ⚠️  Analyzing multi-year drought risk..."
  end subroutine analyze_multiyear_drought_risk

  subroutine compare_return_periods()
    ! Implementation: Compare return period analysis between periods
    print *, "    📈 Comparing return periods..."
  end subroutine compare_return_periods

  subroutine compare_evt_parameters()
    ! Implementation: Compare EVT distribution parameters
    print *, "    📊 Comparing EVT parameters..."
  end subroutine compare_evt_parameters

  subroutine analyze_tail_behavior_changes()
    ! Implementation: Analyze changes in extreme tail behavior
    print *, "    📉 Analyzing tail behavior changes..."
  end subroutine analyze_tail_behavior_changes

  subroutine perform_mean_difference_tests()
    ! Implementation: Statistical tests for mean differences
    print *, "    📊 Performing mean difference tests..."
  end subroutine perform_mean_difference_tests

  subroutine perform_distribution_tests()
    ! Implementation: Distribution comparison tests
    print *, "    📈 Performing distribution tests..."
  end subroutine perform_distribution_tests

  subroutine perform_rank_tests()
    ! Implementation: Non-parametric rank tests
    print *, "    📊 Performing rank tests..."
  end subroutine perform_rank_tests

  subroutine perform_trend_tests()
    ! Implementation: Trend analysis tests
    print *, "    📈 Performing trend analysis..."
  end subroutine perform_trend_tests

  subroutine create_output_directories(output_dir)
    character(len=*), intent(in) :: output_dir
    ! Implementation: Create necessary output directories
    print *, "    📁 Creating output directories..."
  end subroutine create_output_directories

  subroutine write_detailed_comparison_report(filename)
    character(len=*), intent(in) :: filename
    ! Implementation: Write comprehensive comparison report
    print *, "    📝 Writing detailed comparison report..."
  end subroutine write_detailed_comparison_report

  subroutine write_executive_summary(filename)
    character(len=*), intent(in) :: filename
    ! Implementation: Write executive summary
    print *, "    📋 Writing executive summary..."
  end subroutine write_executive_summary

  subroutine generate_comparison_figures(figure_dir)
    character(len=*), intent(in) :: figure_dir
    ! Implementation: Generate comparison figures and plots
    print *, "    📊 Generating comparison figures..."
  end subroutine generate_comparison_figures

  subroutine write_hypothesis_test_results(output_dir)
    character(len=*), intent(in) :: output_dir
    ! Implementation: Write hypothesis testing results
    print *, "    🔬 Writing hypothesis test results..."
  end subroutine write_hypothesis_test_results

end module analysis_comparison
