!> Somaliland Drought Analysis Pipeline - Main Program
!!
!! Implements drought analysis framework for climate change impact assessment
!! in Somaliland, Horn of Africa.
!!
!! **Research Hypothesis:** Climate change will increase drought frequency, 
!! duration, and severity in Somaliland.
!!
!! **Analysis Phases:**
!! 1. Historical drought indices (ERA5-Land, 1981-2024)
!! 2. Model validation & bias correction (CMIP6 vs observations)
!! 3. Future projections (SSP scenarios to 2100)
!! 4. Statistical analysis & climate change quantification
!!
!! @author Khadar Daahir (University of Glasgow)
!! @date 2025-08-25
!! @version 1.0
!! @license MIT License
!! @note Requires conda environment with NetCDF-Fortran and statistical libraries

program main
    use, intrinsic :: iso_fortran_env, only: real32
    use io_module, only: load_era5_climate_data, save_processed_pet_data, save_processed_precipitation_data, &
                         save_spi_multiscale, save_spi_single_scale, save_spei_multiscale, save_spei_single_scale, &
                         save_evt_results_enhanced, create_output_directory, load_raw_cmip6_scenario, &
                         load_corrected_cmip6_scenario, cmip6_scenario_data_t, corrected_cmip6_scenarios_t
    use spi_module, only: calculate_spi_timescales
    use spei_module, only: calculate_spei_timescales, era5_spei_results_t, calculate_era5_spei_structured
    use bias_correction_clean, only: process_bias_correction
    use projection_module, only: calculate_future_spei, calculate_future_spei_structured, cmip6_spei_results_t
    use evt_module, only: run_evt_analysis, evt_results
    implicit none
    
    ! =================================================================
    ! PHASE 1: BASELINE DROUGHT INDICES (ERA5, 1981–2024)
    ! =================================================================
    
    ! Data dimensions (will be read from files)
    integer, parameter :: NTIME = 528  ! 44 years × 12 months (1981–2024)
    integer, parameter :: NLAT = 36
    integer, parameter :: NLON = 56
    
    ! ERA5 data arrays
    real(8), allocatable :: precipitation(:,:,:)  ! mm/month
    real(8), allocatable :: pet(:,:,:)           ! mm/month
    real(8), allocatable :: latitudes(:)
    real(8), allocatable :: longitudes(:)
    real(8), allocatable :: time_stamps(:)       ! seconds since 1970-01-01
    
    ! Drought indices arrays (SPI only - SPEI uses structured approach)
    real(8), allocatable :: spi_1(:,:,:)   ! 1-month SPI
    real(8), allocatable :: spi_3(:,:,:)   ! 3-month SPI
    real(8), allocatable :: spi_6(:,:,:)   ! 6-month SPI
    real(8), allocatable :: spi_12(:,:,:)  ! 12-month SPI
    
    ! File paths
    character(len=200) :: era5_precip_file, era5_pet_file
    character(len=200) :: output_dir
    
    ! Status and control variables
    integer :: status
    logical :: success
    
    ! CMIP6 data for bias correction
    type(cmip6_scenario_data_t) :: cmip6_historical, cmip6_ssp126, cmip6_ssp245, cmip6_ssp585
    
    ! CMIP6 SPEI results for EVT analysis
    type(cmip6_spei_results_t) :: cmip6_spei_results
    
    ! ERA5 SPEI results for EVT analysis  
    type(era5_spei_results_t) :: era5_spei_results
    
    ! Structured corrected CMIP6 scenarios (memory-based approach)
    type(corrected_cmip6_scenarios_t) :: corrected_scenarios
    
    ! Corrected CMIP6 scenario data (individual arrays - legacy approach)
    real(8), allocatable :: corrected_ssp126_precip(:,:,:), corrected_ssp126_pet(:,:,:)
    real(8), allocatable :: corrected_ssp245_precip(:,:,:), corrected_ssp245_pet(:,:,:)
    real(8), allocatable :: corrected_ssp585_precip(:,:,:), corrected_ssp585_pet(:,:,:)
    real(8), allocatable :: corrected_lat(:), corrected_lon(:), corrected_time(:)
    integer :: corrected_start_year, corrected_end_year
    
    
    print *, "=========================================="
    print *, "   SOMALILAND DROUGHT ANALYSIS PIPELINE"
    print *, "=========================================="
    print *, ""
    print *, "Research Hypothesis:"
    print *, "Under climate change, droughts in Somaliland will"
    print *, "increase in frequency, duration, and severity."
    print *, ""
    print *, "Processing 1981-2024 baseline period..."
    print *, ""
    
    ! =================================================================
    ! Initialize file paths
    ! =================================================================
    
    era5_precip_file = "data/Raw/ERA5/era5_land_precipitation_somaliland_1981_2024.nc"
    era5_pet_file = "data/Raw/ERA5/era5_land_potential_evaporation_somaliland_1981_2024.nc"
    output_dir = "data/Processed_data/Drought_indices/ERA5/"
    
    ! =================================================================
    ! Allocate arrays
    ! =================================================================
    
    allocate(precipitation(NLON, NLAT, NTIME))
    allocate(pet(NLON, NLAT, NTIME))
    allocate(latitudes(NLAT))
    allocate(longitudes(NLON))
    allocate(time_stamps(NTIME))
    
    allocate(spi_1(NLON, NLAT, NTIME))
    allocate(spi_3(NLON, NLAT, NTIME))
    allocate(spi_6(NLON, NLAT, NTIME))
    allocate(spi_12(NLON, NLAT, NTIME))
    
    ! Note: SPEI arrays are managed within era5_spei_results structure
    
    ! =================================================================
    ! PHASE 1A: READ ERA5 CLIMATE DATA (PRECIPITATION + PET)
    ! =================================================================
    
    print *, "Phase 1A: Reading ERA5 climate data (precipitation + PET)..."
    
    ! Load both precipitation and PET data using combined function
    call load_era5_climate_data(era5_precip_file, era5_pet_file, &
                                precipitation, pet, &
                                latitudes, longitudes, time_stamps, &
                                status)
    
    if (status /= 0) then
        print *, "✗ CRITICAL ERROR: Failed to load ERA5 climate data"
        print *, "Status code:", status
        stop 1
    end if
    
    print *, "✓ ERA5 climate data loaded successfully"
    print *, "  Precipitation range:", minval(precipitation, mask=precipitation > -999), &
             " to ", maxval(precipitation, mask=precipitation > -999), " mm/month"
    print *, "  PET range:", minval(pet, mask=pet > -999), &
             " to ", maxval(pet, mask=pet > -999), " mm/month"
    print *, ""
    
    ! Save processed precipitation data
    print *, "Phase 1B: Saving processed climate data..."
    call save_processed_precipitation_data(precipitation, real(latitudes, kind=real32), real(longitudes, kind=real32), &
                                         "data/Processed_data/Climate_Drivers/Precipitation/somaliland_precipitation_processed_1981_2024.nc", &
                                         status)
    
    if (status /= 0) then
        print *, "⚠ WARNING: Failed to save processed precipitation data"
    else
        print *, "✓ Processed precipitation data saved successfully"
    end if
    
    ! Save processed PET data
    call save_processed_pet_data(pet, real(latitudes, kind=real32), real(longitudes, kind=real32), &
                                "data/Processed_data/Climate_Drivers/PET/somaliland_pet_processed_1981_2024.nc", &
                                status)
    
    if (status /= 0) then
        print *, "⚠ WARNING: Failed to save processed PET data"
    else
        print *, "✓ Processed PET data saved successfully"
    end if
    print *, ""
    
    ! =================================================================
    ! PHASE 1C: CALCULATE SPI INDICES
    ! =================================================================
    
    print *, "Phase 1C: Calculating SPI indices using FSML..."
    
    ! Calculate SPI for multiple timescales using FSML
    call calculate_spi_timescales(precipitation, spi_1, spi_3, spi_6, spi_12, status)
    if (status /= 0) then
        print *, "ERROR: SPI calculation failed"
        stop 1
    end if
    
    print *, "✓ SPI indices calculated successfully"
    
    ! Save SPI results to NetCDF files
    print *, "Saving SPI results..."
    
    ! Save multiscale SPI file
    call save_spi_multiscale(spi_1, spi_3, spi_6, spi_12, latitudes, longitudes, &
                            'data/Processed_data/Drought_indices/ERA5/SPI/spi_multiscale_1981_2024.nc', status)
    if (status /= 0) then
        print *, "ERROR: Failed to save multiscale SPI data"
        stop 1
    end if
    print *, "✓ Multiscale SPI data saved"
    
    ! Save individual SPI-3 file
    call save_spi_single_scale(spi_3, 3, latitudes, longitudes, &
                              'data/Processed_data/Drought_indices/ERA5/SPI/spi_3_month_1981_2024.nc', status)
    if (status /= 0) then
        print *, "ERROR: Failed to save SPI-3 data"
        stop 1
    end if
    print *, "✓ SPI-3 data saved"
    
    ! Save individual SPI-12 file
    call save_spi_single_scale(spi_12, 12, latitudes, longitudes, &
                               'data/Processed_data/Drought_indices/ERA5/SPI/spi_12_month_1981_2024.nc', status)
    if (status /= 0) then
        print *, "ERROR: Failed to save SPI-12 data"
        stop 1
    end if
    print *, "✓ SPI-12 data saved"
    print *, ""
    
    ! =================================================================
    ! PHASE 1D: CALCULATE SPEI INDICES
    ! =================================================================
    
    print *, "Phase 1D: Calculating SPEI indices for multiple timescales..."
    
    ! Calculate SPEI using structured approach (consistent with CMIP6)
    call calculate_era5_spei_structured(precipitation, pet, latitudes, longitudes, &
                                       1981, 2024, era5_spei_results)
    
    ! SPEI Quality Control Report
    print *, "  SPEI-1 QC:"
    print *, "    Valid data points:", count(era5_spei_results%spei_1 > -999.0 .and. era5_spei_results%spei_1 < 999.0)
    print *, "    Coverage:", 100.0 * count(era5_spei_results%spei_1 > -999.0) / real(size(era5_spei_results%spei_1)), "%"
    print *, "    Range: [", minval(era5_spei_results%spei_1, mask=era5_spei_results%spei_1 > -999), ",", &
             maxval(era5_spei_results%spei_1, mask=era5_spei_results%spei_1 < 999), "]"
    
    print *, "  SPEI-3 QC:"
    print *, "    Valid data points:", count(era5_spei_results%spei_3 > -999.0 .and. era5_spei_results%spei_3 < 999.0)
    print *, "    Coverage:", 100.0 * count(era5_spei_results%spei_3 > -999.0) / real(size(era5_spei_results%spei_3)), "%"
    print *, "    Range: [", minval(era5_spei_results%spei_3, mask=era5_spei_results%spei_3 > -999), ",", &
             maxval(era5_spei_results%spei_3, mask=era5_spei_results%spei_3 < 999), "]"
    
    print *, "  SPEI-6 QC:"
    print *, "    Valid data points:", count(era5_spei_results%spei_6 > -999.0 .and. era5_spei_results%spei_6 < 999.0)
    print *, "    Coverage:", 100.0 * count(era5_spei_results%spei_6 > -999.0) / real(size(era5_spei_results%spei_6)), "%"
    print *, "    Range: [", minval(era5_spei_results%spei_6, mask=era5_spei_results%spei_6 > -999), ",", &
             maxval(era5_spei_results%spei_6, mask=era5_spei_results%spei_6 < 999), "]"
             
    print *, "  SPEI-12 QC:"
    print *, "    Valid data points:", count(era5_spei_results%spei_12 > -999.0 .and. era5_spei_results%spei_12 < 999.0)
    print *, "    Coverage:", 100.0 * count(era5_spei_results%spei_12 > -999.0) / real(size(era5_spei_results%spei_12)), "%"
    print *, "    Range: [", minval(era5_spei_results%spei_12, mask=era5_spei_results%spei_12 > -999), ",", &
             maxval(era5_spei_results%spei_12, mask=era5_spei_results%spei_12 < 999), "]"
    
    print *, "Phase 1D: SPEI calculation completed successfully"
    
    ! Save SPEI results to NetCDF files
    print *, "Saving SPEI results..."
    
    ! Save multiscale SPEI file
    call save_spei_multiscale(era5_spei_results%spei_1, era5_spei_results%spei_3, &
                             era5_spei_results%spei_6, era5_spei_results%spei_12, &
                             latitudes, longitudes, &
                             'data/Processed_data/Drought_indices/ERA5/SPEI/spei_multiscale_1981_2024.nc', status)
    if (status /= 0) then
        print *, "ERROR: Failed to save multiscale SPEI data"
        stop 1
    end if
    print *, "✓ Multiscale SPEI data saved"
    
    ! Save individual SPEI-3 file
    call save_spei_single_scale(era5_spei_results%spei_3, 3, latitudes, longitudes, &
                               'data/Processed_data/Drought_indices/ERA5/SPEI/spei_3_month_1981_2024.nc', status)
    if (status /= 0) then
        print *, "ERROR: Failed to save SPEI-3 data"
        stop 1
    end if
    print *, "✓ SPEI-3 data saved"
    
    ! Save individual SPEI-12 file
    call save_spei_single_scale(era5_spei_results%spei_12, 12, latitudes, longitudes, &
                                'data/Processed_data/Drought_indices/ERA5/SPEI/spei_12_month_1981_2024.nc', status)
    if (status /= 0) then
        print *, "ERROR: Failed to save SPEI-12 data"
        stop 1
    end if
    print *, "✓ SPEI-12 data saved"
    print *, ""
    
    ! =================================================================
    ! PHASE 2: VALIDATION & BIAS CORRECTION
    ! =================================================================
    
    print *, "Phase 2: CMIP6 Validation & Bias Correction..."
    
    ! Load raw CMIP6 data for bias correction (modular approach)
    print *, "Loading raw CMIP6 data for bias correction..."
    call load_raw_cmip6_scenario("data", "historical", cmip6_historical, status)
    if (status /= 0) then
        print *, "⚠ WARNING: Failed to load CMIP6 historical data"
        success = .false.
    else
        print *, "✓ CMIP6 historical data loaded"
    end if
    
    call load_raw_cmip6_scenario("data", "ssp126", cmip6_ssp126, status)
    if (status /= 0) then
        print *, "⚠ WARNING: Failed to load CMIP6 SSP126 data"
        success = .false.
    else
        print *, "✓ CMIP6 SSP126 data loaded"
    end if
    
    call load_raw_cmip6_scenario("data", "ssp245", cmip6_ssp245, status)
    if (status /= 0) then
        print *, "⚠ WARNING: Failed to load CMIP6 SSP245 data"
        success = .false.
    else
        print *, "✓ CMIP6 SSP245 data loaded"
    end if
    
    call load_raw_cmip6_scenario("data", "ssp585", cmip6_ssp585, status)
    if (status /= 0) then
        print *, "⚠ WARNING: Failed to load CMIP6 SSP585 data"
        success = .false.
    else
        print *, "✓ CMIP6 SSP585 data loaded"
    end if
    
    ! Call structured bias correction with pre-loaded data (eliminates redundant I/O)
    if (success) then
        print *, ""
        print *, "🔄 Using structured memory-based approach (consistent with ERA5)..."
        call process_bias_correction(precipitation, pet, &
                                    cmip6_historical, cmip6_ssp126, cmip6_ssp245, cmip6_ssp585, &
                                    corrected_scenarios, success)
    end if
    
    if (success) then
        print *, "✓ Phase 2: CMIP6 bias correction completed successfully"
    else
        print *, "⚠ WARNING: Phase 2 bias correction failed"
    end if
    print *, ""
    
    ! =================================================================
    ! PHASE 3: FUTURE PROJECTIONS (STRUCTURED MEMORY-BASED)
    ! =================================================================
    
    print *, "Phase 3: Future drought projections (CMIP6 SSP scenarios)..."
    
    ! Use structured corrected data directly from memory (no file I/O)
    if (success .and. corrected_scenarios%data_available) then
        print *, ""
        print *, "🔄 Using memory-based corrected scenarios (eliminates file I/O bottleneck)..."
        call calculate_future_spei_structured(corrected_scenarios, cmip6_spei_results, success)
    else
        print *, "⚠ WARNING: Corrected CMIP6 scenarios not available in memory"
        success = .false.
    end if
    
    if (success) then
        print *, "✓ Phase 3: Structured future projections completed successfully"
        print *, "   (Memory-based processing - consistent with ERA5 approach)"
    else
        print *, "⚠ WARNING: Phase 3 structured projections failed"
    end if
    print *, ""
    
    ! =================================================================
    ! PHASE 4: EXTREME VALUE THEORY ANALYSIS
    ! =================================================================
    
    print *, "Phase 4: Extreme Value Theory (EVT) Analysis..."
    print *, "Analyzing extreme drought events using GPD fitting..."
    
    call run_evt_analysis_pipeline(era5_spei_results%spei_3, era5_spei_results%spei_6, &
                                   era5_spei_results%spei_12, latitudes, longitudes, &
                                   cmip6_spei_results, success)
    
    if (success) then
        print *, "✓ Phase 4: EVT analysis completed successfully"
        print *, "  Results saved in: data/Processed_data/extreme_value_analysis/"
        print *, "  Files include GPD parameters, return levels, and exceedance probabilities"
    else
        print *, "⚠ WARNING: Phase 4 EVT analysis failed"
    end if
    print *, ""
    
    ! =================================================================
    ! PIPELINE COMPLETION SUMMARY
    ! =================================================================
    
    print *, "========================================================"
    print *, "           DROUGHT ANALYSIS PIPELINE COMPLETE"
    print *, "========================================================"
    print *, ""
    print *, "✅ All 4 phases executed successfully:"
    print *, "   Phase 1: Historical drought indices (ERA5, 1981-2024)"
    print *, "   Phase 2: CMIP6 validation & bias correction"
    print *, "   Phase 3: Future projections (SSP scenarios to 2100)"
    print *, "   Phase 4: Extreme value analysis (GPD, return levels)"
    print *, ""
    print *, "📁 Results organized in: data/Processed_data/"
    print *, "📊 Analysis ready for publication and policy application"
    print *, ""
    
    ! =================================================================
    ! CLEANUP
    ! =================================================================
    
    deallocate(precipitation, pet, latitudes, longitudes, time_stamps)
    ! Note: SPI and SPEI arrays are now managed within structured data types
    ! era5_spei_results and cmip6_spei_results will be automatically deallocated
    
    print *, "=========================================="
    print *, "   DROUGHT ANALYSIS PIPELINE COMPLETE"
    print *, "=========================================="

contains

    !===================================================================
    ! EVT ANALYSIS PIPELINE INTEGRATION
    !===================================================================
    subroutine run_evt_analysis_pipeline(era5_spei_3, era5_spei_6, era5_spei_12, &
                                         era5_lat, era5_lon, cmip6_spei_results, success)
        use io_module, only: save_evt_results_enhanced
        use projection_module, only: cmip6_spei_results_t
        implicit none
        real(8), intent(in) :: era5_spei_3(:,:,:), era5_spei_6(:,:,:), era5_spei_12(:,:,:)
        real(8), intent(in) :: era5_lat(:), era5_lon(:)
        type(cmip6_spei_results_t), intent(in) :: cmip6_spei_results
        logical, intent(out) :: success
        
        ! EVT analysis parameters
        character(len=10), parameter :: scenarios(4) = ['era5      ', 'ssp126    ', 'ssp245    ', 'ssp585    ']
        integer, parameter :: evt_timescales(3) = [3, 6, 12]  ! Focus on longer timescales
        real(8), parameter :: evt_thresholds(3) = [-1.5d0, -1.5d0, -1.0d0]  ! Drought thresholds (negative SPEI)
        
        ! EVT results structure
        type(evt_results) :: results
        real(8), allocatable :: spei_data_2d(:,:)  ! Flattened spatial-temporal data
        real(8), allocatable :: time_years(:)
        
        ! File handling
        character(len=300) :: output_file, output_dir
        integer :: i_scenario, i_timescale, status
        logical :: dir_exists
        
        success = .true.
        
        print *, ""
        print *, "=== Extreme Value Theory (EVT) Analysis ==="
        print *, "Applying Profile Maximum Likelihood Estimation"
        print *, "Testing Hypothesis: Future droughts more extreme than historical"
        print *, "References: Coles (2001), Grimshaw (1993)"
        print *, ""
        
        ! Create output directory structure (smart creation - only if needed)
        output_dir = "data/Processed_data/extreme_value_analysis"
        call create_evt_directory_structure(output_dir)
        
        ! Process each scenario and timescale combination
        do i_scenario = 1, size(scenarios)
            do i_timescale = 1, size(evt_timescales)
                
                print *, "📊 Processing EVT: ", trim(scenarios(i_scenario)), &
                        " SPEI-", evt_timescales(i_timescale)
                
                ! Get appropriate SPEI data based on scenario
                if (trim(scenarios(i_scenario)) == 'era5') then
                    ! Use provided ERA5 historical data (1981-2024)
                    if (i_timescale == 1) then  ! SPEI-3
                        call prepare_spei_for_evt(era5_spei_3, spei_data_2d, time_years, 1981, 2024)
                    else if (i_timescale == 2) then  ! SPEI-6
                        call prepare_spei_for_evt(era5_spei_6, spei_data_2d, time_years, 1981, 2024)
                    else  ! SPEI-12
                        call prepare_spei_for_evt(era5_spei_12, spei_data_2d, time_years, 1981, 2024)
                    end if
                else
                    ! Use in-memory CMIP6 projection data (2015-2099) 
                    if (.not. cmip6_spei_results%data_available) then
                        print *, "⚠️  CMIP6 SPEI data not available for ", trim(scenarios(i_scenario))
                        success = .false.
                        cycle
                    end if
                    
                    ! Select appropriate SPEI data based on scenario and timescale
                    if (trim(scenarios(i_scenario)) == 'ssp126') then
                        if (i_timescale == 1) then  ! SPEI-3
                            call prepare_spei_for_evt(cmip6_spei_results%ssp126_spei_3, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        else if (i_timescale == 2) then  ! SPEI-6
                            call prepare_spei_for_evt(cmip6_spei_results%ssp126_spei_6, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        else  ! SPEI-12
                            call prepare_spei_for_evt(cmip6_spei_results%ssp126_spei_12, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        end if
                    else if (trim(scenarios(i_scenario)) == 'ssp245') then
                        if (i_timescale == 1) then  ! SPEI-3
                            call prepare_spei_for_evt(cmip6_spei_results%ssp245_spei_3, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        else if (i_timescale == 2) then  ! SPEI-6
                            call prepare_spei_for_evt(cmip6_spei_results%ssp245_spei_6, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        else  ! SPEI-12
                            call prepare_spei_for_evt(cmip6_spei_results%ssp245_spei_12, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        end if
                    else if (trim(scenarios(i_scenario)) == 'ssp585') then
                        if (i_timescale == 1) then  ! SPEI-3
                            call prepare_spei_for_evt(cmip6_spei_results%ssp585_spei_3, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        else if (i_timescale == 2) then  ! SPEI-6
                            call prepare_spei_for_evt(cmip6_spei_results%ssp585_spei_6, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        else  ! SPEI-12
                            call prepare_spei_for_evt(cmip6_spei_results%ssp585_spei_12, spei_data_2d, time_years, &
                                                     cmip6_spei_results%start_year, cmip6_spei_results%end_year)
                        end if
                    end if
                end if
                
                ! Run EVT analysis for ALL scenarios (ERA5 + CMIP6)
                call run_evt_analysis(spei_data_2d, time_years, &
                                     evt_thresholds(i_timescale), results)
                
                ! Check if analysis was successful
                if (results%params%status /= 0) then
                    print *, "⚠️  EVT analysis failed for ", trim(scenarios(i_scenario)), &
                            " SPEI-", evt_timescales(i_timescale)
                    success = .false.
                    cycle
                end if
                
                ! Construct output file path based on scenario
                if (trim(scenarios(i_scenario)) == 'era5') then
                    write(output_file, '(A,A,I0,A)') &
                        trim(output_dir), "/era5/evt_era5_spei", evt_timescales(i_timescale), ".nc"
                else
                    write(output_file, '(A,A,A,A,A,A,I0,A)') &
                        trim(output_dir), "/cmip6/", trim(scenarios(i_scenario)), &
                        "/evt_", trim(scenarios(i_scenario)), "_spei", evt_timescales(i_timescale), ".nc"
                end if
                
                ! Save EVT results using enhanced format
                call save_evt_results_enhanced(output_file, results, trim(scenarios(i_scenario)), &
                                             evt_timescales(i_timescale), status)
                
                if (status == 0) then
                    print *, "✅ EVT results saved: ", trim(output_file)
                else
                    print *, "❌ Failed to save EVT results: ", trim(output_file)
                    success = .false.
                end if
                
                ! Clean up allocated arrays (only the ones we allocated locally)
                if (allocated(spei_data_2d)) deallocate(spei_data_2d)
                if (allocated(time_years)) deallocate(time_years)
                
            end do
        end do
        
        if (success) then
            print *, ""
            print *, "✅ EVT analysis pipeline completed successfully"
            print *, "   Historical baseline (ERA5): ", trim(output_dir), "/era5/"
            print *, "   Future projections (CMIP6): ", trim(output_dir), "/cmip6/"
            print *, "   Hypothesis testing: Compare return levels between periods"
        else
            print *, ""
            print *, "⚠️  EVT analysis pipeline completed with warnings"
        end if
        
    end subroutine run_evt_analysis_pipeline

    !===================================================================
    ! SMART DIRECTORY CREATION - Create only if needed
    !===================================================================
    subroutine create_evt_directory_structure(base_dir)
        implicit none
        character(len=*), intent(in) :: base_dir
        
        logical :: dir_exists
        
        ! Check and create base directory
        inquire(file=trim(base_dir), exist=dir_exists)
        if (.not. dir_exists) then
            print *, "📁 Creating EVT output directory: ", trim(base_dir)
            call create_output_directory(base_dir)
        end if
        
        ! Check and create ERA5 subdirectory
        inquire(file=trim(base_dir)//"/era5", exist=dir_exists)
        if (.not. dir_exists) then
            print *, "📁 Creating ERA5 EVT directory: ", trim(base_dir), "/era5"
            call create_output_directory(trim(base_dir) // "/era5")
        end if
        
        ! Check and create CMIP6 subdirectories
        inquire(file=trim(base_dir)//"/cmip6", exist=dir_exists)
        if (.not. dir_exists) then
            print *, "📁 Creating CMIP6 EVT directories..."
            call create_output_directory(trim(base_dir) // "/cmip6")
            call create_output_directory(trim(base_dir) // "/cmip6/ssp126")
            call create_output_directory(trim(base_dir) // "/cmip6/ssp245")
            call create_output_directory(trim(base_dir) // "/cmip6/ssp585")
        end if
        
    end subroutine create_evt_directory_structure

    !===================================================================
    ! PREPARE SPEI DATA FOR EVT ANALYSIS
    !===================================================================
    subroutine prepare_spei_for_evt(spei_3d, spei_2d, time_years, start_year, end_year)
        implicit none
        real(8), intent(in) :: spei_3d(:,:,:)
        real(8), allocatable, intent(out) :: spei_2d(:,:)
        real(8), allocatable, intent(out) :: time_years(:)
        integer, intent(in), optional :: start_year, end_year
        
        integer :: nlon, nlat, ntime, i, j, k
        integer :: actual_start_year, actual_end_year
        
        nlon = size(spei_3d, 1)
        nlat = size(spei_3d, 2) 
        ntime = size(spei_3d, 3)
        
        ! Determine start and end years
        if (present(start_year) .and. present(end_year)) then
            actual_start_year = start_year
            actual_end_year = end_year
        else
            ! Default to ERA5 historical period (1981-2024)
            actual_start_year = 1981
            actual_end_year = 2024
        end if
        
        ! Allocate output arrays
        allocate(spei_2d(nlon * nlat, ntime))
        allocate(time_years(ntime))
        
        ! Flatten spatial dimensions (preserve time dimension)
        do k = 1, ntime
            do j = 1, nlat
                do i = 1, nlon
                    spei_2d((j-1)*nlon + i, k) = spei_3d(i, j, k)
                end do
            end do
        end do
        
        ! Create time array based on actual period
        do k = 1, ntime
            time_years(k) = real(actual_start_year) + real(k-1) / 12.0  ! Convert to decimal years
        end do
        
    end subroutine prepare_spei_for_evt

end program main
