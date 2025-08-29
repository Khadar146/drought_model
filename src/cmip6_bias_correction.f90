module cmip6_bias_correction
    !===============================================================================
    ! CMIP6 BIAS CORRECTION MODULE -  DRIVER-BASED APPROACH
    !===============================================================================
    !
    ! Purpose: Proper bias correction following climate science best practices
    !
    !  WORKFLOW (PET-from-drivers + direct PET correction):
    !   Phase 2A - Historical Calibration (1981-2014):
    !     1. Load CMIP6-historical P, Tmax, Tmin for overlap period
    !     2. Calculate CMIP6-historical PET from Tmax/Tmin (Hargreaves method)
    !     3. Load ERA5 P and ERA5 PET for same period
    !     4. Compute bias correction factors:
    !        • For P: CMIP6-historical P vs ERA5 P (multiplicative)
    !        • For PET: CMIP6-historical PET vs ERA5 PET (multiplicative)
    !     5. Apply corrections to CMIP6-historical P and PET
    !     6. Save corrected historical data (P_corrected, PET_corrected)
    !
    !   Phase 2B - Future Application (2015-2099):
    !     1. Load CMIP6-future P, Tmax, Tmin (scenarios: ssp126, ssp245, ssp585)
    
    !     2. Calculate CMIP6-future PET from Tmax/Tmin (same method as Phase 2A)
    !     3. Apply SAME bias correction factors from Phase 2A:
    !        • Correct P_future using P bias factors
    !        • Correct PET_future using PET bias factors
    !     4. Save corrected future data (P_corrected, PET_corrected)
    !
    ! Key Principle: PET-from-drivers + direct PET correction approach
    !               (since we have ERA5 PET but not ERA5 Tmax/Tmin)
    !
    ! Author: Khadar Daahir
    ! Date: August 2025
    !===============================================================================
    
    use iso_fortran_env, only: real64, output_unit
    use ieee_arithmetic  ! For NaN handling
    use netcdf
    use io_module, only: save_corrected_cmip6_precipitation, save_corrected_cmip6_pet_hargreaves, &
                         cmip6_scenario_data_t, load_raw_cmip6_scenario, save_corrected_cmip6_scenario, &
                         corrected_cmip6_scenarios_t
    implicit none
    
    private
    
    ! Public interfaces
    public :: bias_correction_factors_t
    public :: process_cmip6_bias_correction_proper, process_cmip6_bias_correction_structured
    public :: load_cmip6_scenario_data
    public :: calculate_hargreaves_pet_from_cmip6_drivers
    
    ! Data structures
    
    type :: bias_correction_factors_t
        ! Bias correction factors (monthly)
        real(real64), allocatable :: precip_bias_factor(:,:,:)    ! multiplicative [12 months]
        real(real64), allocatable :: pet_bias_factor(:,:,:)       ! multiplicative [12 months]
        
        ! Variance scaling factors (monthly)
        real(real64), allocatable :: precip_var_factor(:,:,:)     ! variance scaling [12 months]
        real(real64), allocatable :: pet_var_factor(:,:,:)        ! variance scaling [12 months]
        
        ! Quality metrics
        real(real64), allocatable :: precip_correlation(:,:)      ! spatial correlation
        real(real64), allocatable :: pet_correlation(:,:)         ! spatial correlation
        real(real64), allocatable :: precip_rmse(:,:)             ! RMSE [mm/month]
        real(real64), allocatable :: pet_rmse(:,:)                ! RMSE [mm/month]
        
        ! Metadata
        integer :: nlons, nlats
        logical :: is_computed
        character(len=50) :: reference_period                     ! "1981-2014"
    end type bias_correction_factors_t
    
    ! Module parameters
    integer, parameter :: HISTORICAL_START_YEAR = 1981
    integer, parameter :: HISTORICAL_END_YEAR = 2014
    integer, parameter :: OVERLAP_MONTHS = (HISTORICAL_END_YEAR - HISTORICAL_START_YEAR + 1) * 12
    real(real64), parameter :: MISSING_VALUE = -999.0_real64
    real(real64), parameter :: MIN_PRECIP_THRESHOLD = 0.1_real64  ! mm/month
    real(real64), parameter :: KELVIN_TO_CELSIUS = 273.15_real64
    real(real64), parameter :: LAND_THRESHOLD = 0.75_real64     ! 75% valid data threshold for land pixels
    real(real64), parameter :: SECONDS_PER_MONTH = 2.628e6_real64  ! Conversion factor kg/m²/s to mm/month
    
    ! Solar radiation constants for proper Hargreaves calculation
    real(real64), parameter :: PI = 3.14159265359_real64
    real(real64), parameter :: SOLAR_CONSTANT = 1367.0_real64  ! W/m² (solar constant)
    real(real64), parameter :: HARGREAVES_COEFF = 0.0023_real64
    
contains

    !===========================================================================
    ! STRUCTURED BIAS CORRECTION: Memory-based approach for consistency
    !===========================================================================
    subroutine process_cmip6_bias_correction_structured(era5_precip, era5_pet, era5_lats, era5_lons, &
                                                       cmip6_historical, cmip6_ssp126, cmip6_ssp245, cmip6_ssp585, &
                                                       corrected_scenarios, success)
        !===========================================================================
        ! Memory-based CMIP6 bias correction - consistent with ERA5 approach
        ! Accepts pre-loaded CMIP6 data, returns structured corrected data
        !===========================================================================
        real(real64), intent(in) :: era5_precip(:,:,:)    ! ERA5 reference [mm/month]
        real(real64), intent(in) :: era5_pet(:,:,:)       ! ERA5 reference [mm/month] 
        real(real64), intent(in) :: era5_lats(:)          ! ERA5 latitudes
        real(real64), intent(in) :: era5_lons(:)          ! ERA5 longitudes
        type(cmip6_scenario_data_t), intent(in) :: cmip6_historical, cmip6_ssp126, cmip6_ssp245, cmip6_ssp585
        type(corrected_cmip6_scenarios_t), intent(out) :: corrected_scenarios
        logical, intent(out) :: success
        
        ! Working copies for bias correction processing
        type(cmip6_scenario_data_t) :: work_historical, work_ssp126, work_ssp245, work_ssp585
        type(bias_correction_factors_t) :: bias_factors
        logical :: phase_success
        
        success = .false.
        
        print *, ""
        print *, "=========================================="
        print *, "   STRUCTURED CMIP6 BIAS CORRECTION"
        print *, "=========================================="
        print *, ""
        print *, "Approach: Memory-based structured data processing"
        print *, "  1. Use pre-loaded CMIP6 data (eliminate redundant I/O)"
        print *, "  2. Calculate CMIP6 PET from Tmax/Tmin (Hargreaves)"
        print *, "  3. Compute bias correction factors"
        print *, "  4. Apply corrections to all scenarios"
        print *, "  5. Return structured corrected data"
        print *, ""
        
        ! Make working copies of input data (to avoid modifying originals)
        work_historical = cmip6_historical
        work_ssp126 = cmip6_ssp126
        work_ssp245 = cmip6_ssp245
        work_ssp585 = cmip6_ssp585
        
        ! Calculate PET for all scenarios using Hargreaves method
        call calculate_hargreaves_pet_from_cmip6_drivers(work_historical, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to calculate historical PET from drivers"
            return
        end if
        
        call calculate_hargreaves_pet_from_cmip6_drivers(work_ssp126, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to calculate SSP126 PET from drivers"
            return
        end if
        
        call calculate_hargreaves_pet_from_cmip6_drivers(work_ssp245, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to calculate SSP245 PET from drivers"
            return
        end if
        
        call calculate_hargreaves_pet_from_cmip6_drivers(work_ssp585, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to calculate SSP585 PET from drivers"
            return
        end if
        
        ! Compute bias correction factors using historical period
        call compute_bias_correction_factors_correct(work_historical, era5_precip, era5_pet, &
                                                    era5_lats, era5_lons, bias_factors, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to compute bias correction factors"
            return
        end if
        
        ! Apply bias corrections to all scenarios
        call apply_precipitation_pet_bias_correction(work_historical, bias_factors, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to apply bias correction to historical"
            return
        end if
        
        call apply_precipitation_pet_bias_correction(work_ssp126, bias_factors, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to apply bias correction to SSP126"
            return
        end if
        
        call apply_precipitation_pet_bias_correction(work_ssp245, bias_factors, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to apply bias correction to SSP245"
            return
        end if
        
        call apply_precipitation_pet_bias_correction(work_ssp585, bias_factors, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to apply bias correction to SSP585"
            return
        end if
        
        ! Package results into structured output (memory-based)
        call package_corrected_scenarios(work_ssp126, work_ssp245, work_ssp585, corrected_scenarios, phase_success)
        if (.not. phase_success) then
            print *, "❌ Failed to package corrected scenarios"
            return
        end if
        
        success = .true.
        print *, "✅ Structured CMIP6 bias correction completed successfully"
        print *, "   Data available in memory for direct processing"
        print *, ""
        
    end subroutine process_cmip6_bias_correction_structured

    subroutine process_cmip6_bias_correction_proper(data_directory, era5_precip, era5_pet, &
                                                   era5_lats, era5_lons, success)
        !===========================================================================
        ! Main driver for PROPER CMIP6 bias correction process
        !===========================================================================
        character(len=*), intent(in) :: data_directory
        real(real64), intent(in) :: era5_precip(:,:,:)    ! ERA5 reference [mm/month]
        real(real64), intent(in) :: era5_pet(:,:,:)       ! ERA5 reference [mm/month] 
        real(real64), intent(in) :: era5_lats(:)          ! ERA5 latitudes
        real(real64), intent(in) :: era5_lons(:)          ! ERA5 longitudes
        logical, intent(out) :: success
        
        ! CMIP6 scenario data
        type(cmip6_scenario_data_t) :: cmip6_historical, cmip6_ssp126, cmip6_ssp245, cmip6_ssp585
        type(bias_correction_factors_t) :: bias_factors
        
        logical :: phase_success
        
        success = .false.
        
        write(*,*) ""
        write(*,*) "=========================================="
        write(*,*) "   PHASE 2: CORRECTED CMIP6 BIAS CORRECTION"
        write(*,*) "=========================================="
        write(*,*) ""
        write(*,*) "Approach: PET-from-drivers + direct PET correction"
        write(*,*) "  1. Calculate CMIP6 PET from Tmax/Tmin (Hargreaves)"
        write(*,*) "  2. Compare CMIP6 PET vs ERA5 PET"
        write(*,*) "  3. Bias correct P and PET directly"
        write(*,*) "  4. Historical calibration (1981-2014)"
        write(*,*) "  5. Apply to future scenarios (2015-2099)"
        write(*,*) ""
        
        ! ===================================================================
        ! PHASE 2A: HISTORICAL CALIBRATION (1981-2014)
        ! ===================================================================
        
        write(*,*) "PHASE 2A: Historical Calibration (1981-2014)"
        write(*,*) "============================================="
        
        ! Step 1: Load CMIP6 historical data
        call load_cmip6_scenario_data(data_directory, "historical", cmip6_historical, phase_success)
        if (.not. phase_success) then
            write(*,*) "❌ Failed to load CMIP6 historical data"
            return
        end if
        
        ! Step 2: Calculate CMIP6 PET from drivers (Hargreaves method)
        call calculate_hargreaves_pet_from_cmip6_drivers(cmip6_historical, phase_success)
        if (.not. phase_success) then
            write(*,*) "❌ Failed to calculate CMIP6 historical PET from drivers"
            return
        end if
        
        ! Step 3: Compute bias correction factors (P vs ERA5_P, PET vs ERA5_PET)
        call compute_bias_correction_factors_correct(cmip6_historical, era5_precip, era5_pet, era5_lats, era5_lons, &
                                                    bias_factors, phase_success)
        if (.not. phase_success) then
            write(*,*) "❌ Failed to compute bias correction factors"
            return
        end if
        
        ! Step 4: Apply bias corrections to historical P and PET
        call apply_precipitation_pet_bias_correction(cmip6_historical, bias_factors, phase_success)
        if (.not. phase_success) then
            write(*,*) "❌ Failed to apply bias correction to historical P and PET"
            return
        end if
        
        ! Step 5: Save corrected historical data using io_module
        call save_corrected_cmip6_scenario(cmip6_historical, "historical", phase_success)
        if (.not. phase_success) then
            write(*,*) "❌ Failed to save corrected historical data"
            return
        end if
        
        write(*,*) "✅ PHASE 2A: Historical calibration completed successfully"
        write(*,*) ""
        
        ! ===================================================================
        ! PHASE 2B: FUTURE PROJECTIONS APPLICATION (2015-2099)
        ! ===================================================================
        
        write(*,*) "PHASE 2B: Future Projections Application (2015-2099)"
        write(*,*) "===================================================="
        
        ! Process each future scenario using the same bias correction factors
        call process_future_scenario("ssp126", cmip6_ssp126, data_directory, bias_factors, phase_success)
        if (.not. phase_success) then
            write(*,*) "❌ Failed to process SSP1-2.6 scenario"
            return
        end if
        
        call process_future_scenario("ssp245", cmip6_ssp245, data_directory, bias_factors, phase_success)
        if (.not. phase_success) then
            write(*,*) "❌ Failed to process SSP2-4.5 scenario"
            return
        end if
        
        call process_future_scenario("ssp585", cmip6_ssp585, data_directory, bias_factors, phase_success)
        if (.not. phase_success) then
            write(*,*) "❌ Failed to process SSP5-8.5 scenario"
            return
        end if
        
        write(*,*) "✅ PHASE 2B: Future projections completed successfully"
        write(*,*) ""
        
        write(*,*) "=========================================="
        write(*,*) "   PHASE 2 COMPLETE: CORRECTED BIAS CORRECTION"
        write(*,*) "=========================================="
        write(*,*) ""
        write(*,*) "Approach Used: PET-from-drivers + direct PET correction"
        write(*,*) "  ✓ CMIP6 drivers → CMIP6 PET (Hargreaves method)"
        write(*,*) "  ✓ Bias correction: P vs ERA5_P, PET vs ERA5_PET"
        write(*,*) "  ✓ Historical calibration (1981-2014) applied to futures"
        write(*,*) ""
        write(*,*) "Output saved to:"
        write(*,*) "  data/Processed_data/Climate_Drivers/Corrected_CMIP6/"
        write(*,*) "    ├── cmip6_historical/precipitation/ & PET/"
        write(*,*) "    ├── cmip6_ssp126/precipitation/ & PET/"
        write(*,*) "    ├── cmip6_ssp245/precipitation/ & PET/"
        write(*,*) "    └── cmip6_ssp585/precipitation/ & PET/"
        
        success = .true.
        
    end subroutine process_cmip6_bias_correction_proper

    subroutine process_future_scenario(scenario_name, cmip6_data, data_directory, bias_factors, success)
        !===========================================================================
        ! Process a future scenario using historical bias correction factors
        !===========================================================================
        character(len=*), intent(in) :: scenario_name
        type(cmip6_scenario_data_t), intent(out) :: cmip6_data
        character(len=*), intent(in) :: data_directory
        type(bias_correction_factors_t), intent(in) :: bias_factors
        logical, intent(out) :: success
        
        success = .false.
        
        write(*,*) "📂 Processing ", trim(scenario_name), " scenario..."
        
        ! Step 1: Load future scenario data
        call load_cmip6_scenario_data(data_directory, scenario_name, cmip6_data, success)
        if (.not. success) then
            write(*,*) "❌ Failed to load ", trim(scenario_name), " data"
            return
        end if
        
        ! Step 2: Calculate CMIP6 PET from drivers (same method as historical)
        call calculate_hargreaves_pet_from_cmip6_drivers(cmip6_data, success)
        if (.not. success) then
            write(*,*) "❌ Failed to calculate CMIP6 ", trim(scenario_name), " PET from drivers"
            return
        end if
        
        ! Step 3: Apply SAME bias correction factors from historical period
        call apply_precipitation_pet_bias_correction(cmip6_data, bias_factors, success)
        if (.not. success) then
            write(*,*) "❌ Failed to apply bias correction to ", trim(scenario_name), " P and PET"
            return
        end if
        
        ! Step 4: Save corrected data using io_module
        call save_corrected_cmip6_scenario(cmip6_data, scenario_name, success)
        if (.not. success) then
            write(*,*) "❌ Failed to save corrected ", trim(scenario_name), " data"
            return
        end if
        
        write(*,*) "✅ ", trim(scenario_name), " scenario processed successfully"
        
    end subroutine process_future_scenario

    subroutine load_cmip6_scenario_data(data_directory, scenario, cmip6_data, success)
        !===========================================================================
        ! Load CMIP6 scenario data - now delegated to io_module
        !===========================================================================
        character(len=*), intent(in) :: data_directory
        character(len=*), intent(in) :: scenario
        type(cmip6_scenario_data_t), intent(out) :: cmip6_data
        logical, intent(out) :: success
        
        integer :: status
        
        ! Delegate to io_module for centralized I/O
        call load_raw_cmip6_scenario(data_directory, scenario, cmip6_data, status)
        success = (status == 0)
        
    end subroutine load_cmip6_scenario_data

    subroutine build_cmip6_file_paths(base_dir, scenario, pr_file, tasmax_file, tasmin_file)
        !===========================================================================
        ! Build file paths for CMIP6 data based on actual directory structure
        !===========================================================================
        character(len=*), intent(in) :: base_dir, scenario
        character(len=*), intent(out) :: pr_file, tasmax_file, tasmin_file
        
        character(len=256) :: scenario_dir
        
        ! Build scenario directory path
        scenario_dir = trim(base_dir) // "/Raw/CMIP6/cmip6_" // trim(scenario)
        
        ! Build file paths based on actual structure
        if (trim(scenario) == "historical") then
            pr_file = trim(scenario_dir) // "/precipitation/pr_Amon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_19810116-20141216.nc"
            tasmax_file = trim(scenario_dir) // "/PET/Daily maximum near-surface air temperature/tasmax_Amon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_19810116-20141216.nc"
            tasmin_file = trim(scenario_dir) // "/PET/Daily minimum near-surface air temperature/tasmin_Amon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_19810116-20141216.nc"
        else
            ! For SSP scenarios (2015-2099)
            pr_file = trim(scenario_dir) // "/precipitation/pr_Amon_MPI-ESM1-2-LR_" // trim(scenario) // "_r1i1p1f1_gn_20150116-20991216.nc"
            tasmax_file = trim(scenario_dir) // "/PET/Daily maximum near-surface air temperature/tasmax_Amon_MPI-ESM1-2-LR_" // trim(scenario) // "_r1i1p1f1_gn_20150116-20991216.nc"
            tasmin_file = trim(scenario_dir) // "/PET/Daily minimum near-surface air temperature/tasmin_Amon_MPI-ESM1-2-LR_" // trim(scenario) // "_r1i1p1f1_gn_20150116-20991216.nc"
        end if
                      
    end subroutine build_cmip6_file_paths

    !===========================================================================
    ! Apply quality control and land masking to CMIP6 data (similar to ERA5)
    !===========================================================================
    subroutine apply_quality_control_cmip6(cmip6_data)
        type(cmip6_scenario_data_t), intent(inout) :: cmip6_data
        
        integer :: i, j, t, valid_count
        real(real64) :: valid_frac
        real(real64), parameter :: LAND_THRESHOLD = 0.75_real64  ! From ERA5 processing
        real(real64), parameter :: MIN_VALID_PRECIP = 0.0_real64
        real(real64), parameter :: MAX_VALID_PRECIP = 1000.0_real64  ! mm/month
        real(real64), parameter :: MIN_VALID_TEMP = 200.0_real64     ! K
        real(real64), parameter :: MAX_VALID_TEMP = 330.0_real64     ! K
        
        write(*,*) "🔍 Applying quality control to CMIP6 data..."
        
        ! Initialize quality arrays
        cmip6_data%land_mask = .true.
        cmip6_data%data_quality = 1.0_real64
        
        ! Create land mask based on temporal data coverage
        do i = 1, cmip6_data%nlons
            do j = 1, cmip6_data%nlats
                valid_count = 0
                
                ! Count valid data points over time
                do t = 1, cmip6_data%ntimes
                    ! Check if all variables are valid (within reasonable bounds)
                    if (cmip6_data%precipitation(i, j, t) >= MIN_VALID_PRECIP .and. &
                        cmip6_data%precipitation(i, j, t) <= MAX_VALID_PRECIP .and. &
                        cmip6_data%tasmax(i, j, t) >= MIN_VALID_TEMP .and. &
                        cmip6_data%tasmax(i, j, t) <= MAX_VALID_TEMP .and. &
                        cmip6_data%tasmin(i, j, t) >= MIN_VALID_TEMP .and. &
                        cmip6_data%tasmin(i, j, t) <= MAX_VALID_TEMP .and. &
                        cmip6_data%tasmax(i, j, t) > cmip6_data%tasmin(i, j, t)) then
                        valid_count = valid_count + 1
                    end if
                end do
                
                ! Calculate valid data fraction
                valid_frac = real(valid_count, real64) / real(cmip6_data%ntimes, real64)
                cmip6_data%data_quality(i, j) = valid_frac
                
                ! Apply land threshold (following ERA5 approach)
                if (valid_frac < LAND_THRESHOLD) then
                    cmip6_data%land_mask(i, j) = .false.
                end if
            end do
        end do
        
        ! Apply land mask to data - set ocean points to missing
        do i = 1, cmip6_data%nlons
            do j = 1, cmip6_data%nlats
                if (.not. cmip6_data%land_mask(i, j)) then
                    cmip6_data%precipitation(i, j, :) = ieee_value(0.0_real64, ieee_quiet_nan)
                    cmip6_data%tasmax(i, j, :) = ieee_value(0.0_real64, ieee_quiet_nan)
                    cmip6_data%tasmin(i, j, :) = ieee_value(0.0_real64, ieee_quiet_nan)
                end if
            end do
        end do
        
        cmip6_data%has_quality_control = .true.
        
        write(*,*) "✅ Quality control applied:"
        write(*,*) "   Land fraction: ", real(count(cmip6_data%land_mask)) / real(size(cmip6_data%land_mask)) * 100.0, "%"
        write(*,*) "   Average data quality: ", sum(cmip6_data%data_quality) / size(cmip6_data%data_quality) * 100.0, "%"
        
    end subroutine apply_quality_control_cmip6
    
    !===========================================================================
    ! Helper function to check for NaN values
    !===========================================================================
    logical function is_nan(value)
        real(real64), intent(in) :: value
        is_nan = (value /= value)  ! NaN property: NaN != NaN
    end function is_nan

    subroutine compute_bias_correction_factors_correct(cmip6_hist, era5_precip, era5_pet, era5_lats, era5_lons, &
                                                      bias_factors, success)
        !===========================================================================
        ! Compute real bias correction factors using statistical comparison
        ! Compare CMIP6-calculated vs ERA5 variables for overlap period (1981-2014)
        !===========================================================================
        type(cmip6_scenario_data_t), intent(in) :: cmip6_hist
        real(real64), intent(in) :: era5_precip(:,:,:)    ! ERA5 reference [mm/month]
        real(real64), intent(in) :: era5_pet(:,:,:)       ! ERA5 reference [mm/month]
        real(real64), intent(in) :: era5_lats(:)          ! ERA5 latitudes
        real(real64), intent(in) :: era5_lons(:)          ! ERA5 longitudes
        type(bias_correction_factors_t), intent(out) :: bias_factors
        logical, intent(out) :: success
        
        integer :: i, j, m, t, month_idx, year_start_idx, valid_count
        integer :: era5_overlap_start, era5_overlap_end
        real(real64) :: era5_monthly_mean, cmip6_monthly_mean, era5_monthly_var, cmip6_monthly_var
        real(real64) :: era5_annual_mean, cmip6_annual_mean, correlation_sum, correlation_count
        real(real64) :: era5_regional_mean  ! For spatial averaging
        real(real64), allocatable :: era5_month_values(:), cmip6_month_values(:)
        real(real64), allocatable :: era5_annual_values(:), cmip6_annual_values(:)
        logical :: has_valid_data
        
        success = .false.
        
        write(*,*) "📊 Computing real bias correction factors using statistical comparison..."
        write(*,*) "   Reference period: ", HISTORICAL_START_YEAR, "-", HISTORICAL_END_YEAR
        write(*,*) "   ERA5 grid: ", size(era5_lons), " x ", size(era5_lats)
        write(*,*) "   CMIP6 grid: ", cmip6_hist%nlons, " x ", cmip6_hist%nlats
        write(*,*) "   Computing actual ERA5/CMIP6 statistics..."
        
        ! Check if CMIP6 PET has been calculated
        if (.not. allocated(cmip6_hist%pet_raw)) then
            write(*,*) "❌ CMIP6 PET not calculated yet - call calculate_hargreaves_pet_from_cmip6_drivers first"
            return
        end if
        
        ! Initialize bias factors structure
        bias_factors%nlons = cmip6_hist%nlons
        bias_factors%nlats = cmip6_hist%nlats
        bias_factors%reference_period = "1981-2014"
        
        allocate(bias_factors%precip_bias_factor(cmip6_hist%nlons, cmip6_hist%nlats, 12))
        allocate(bias_factors%pet_bias_factor(cmip6_hist%nlons, cmip6_hist%nlats, 12))
        allocate(bias_factors%precip_var_factor(cmip6_hist%nlons, cmip6_hist%nlats, 12))
        allocate(bias_factors%pet_var_factor(cmip6_hist%nlons, cmip6_hist%nlats, 12))
        allocate(bias_factors%precip_correlation(cmip6_hist%nlons, cmip6_hist%nlats))
        allocate(bias_factors%pet_correlation(cmip6_hist%nlons, cmip6_hist%nlats))
        allocate(bias_factors%precip_rmse(cmip6_hist%nlons, cmip6_hist%nlats))
        allocate(bias_factors%pet_rmse(cmip6_hist%nlons, cmip6_hist%nlats))
        
        ! Calculate indices for ERA5 overlap period (1981-2014)
        era5_overlap_start = 1
        era5_overlap_end = min(OVERLAP_MONTHS, size(era5_precip, 3))  ! 34 years * 12 months
        
        write(*,*) "   ERA5 overlap indices: ", era5_overlap_start, " to ", era5_overlap_end
        write(*,*) "   CMIP6 historical length: ", cmip6_hist%ntimes
        
        ! Allocate temporary arrays for monthly statistics
        allocate(era5_month_values(HISTORICAL_END_YEAR - HISTORICAL_START_YEAR + 1))
        allocate(cmip6_month_values(HISTORICAL_END_YEAR - HISTORICAL_START_YEAR + 1))
        allocate(era5_annual_values(HISTORICAL_END_YEAR - HISTORICAL_START_YEAR + 1))
        allocate(cmip6_annual_values(HISTORICAL_END_YEAR - HISTORICAL_START_YEAR + 1))
        
        ! Initialize with neutral values
        bias_factors%precip_bias_factor = 1.0_real64
        bias_factors%pet_bias_factor = 1.0_real64
        bias_factors%precip_var_factor = 1.0_real64
        bias_factors%pet_var_factor = 1.0_real64
        bias_factors%precip_correlation = 0.0_real64
        bias_factors%pet_correlation = 0.0_real64
        bias_factors%precip_rmse = 0.0_real64
        bias_factors%pet_rmse = 0.0_real64
        
        ! For proper spatial comparison, use regional averages instead of single points
        ! (in practice, you'd do proper spatial interpolation between grids)
        write(*,*) "   Computing grid-based bias factors with regional averaging..."
        
        ! Compute monthly bias factors for each grid point
        do m = 1, 12  ! For each month
            do j = 1, cmip6_hist%nlats
                do i = 1, cmip6_hist%nlons
                    
                    ! Only process land points
                    if (.not. cmip6_hist%land_mask(i, j)) cycle
                    
                    ! Extract monthly values for this month across all years
                    valid_count = 0
                    era5_monthly_mean = 0.0_real64
                    cmip6_monthly_mean = 0.0_real64
                    
                    ! Use regional average from ERA5 instead of single point
                    do year_start_idx = 1, cmip6_hist%ntimes, 12
                        month_idx = year_start_idx + m - 1
                        if (month_idx <= cmip6_hist%ntimes .and. month_idx <= era5_overlap_end) then
                            ! Calculate regional ERA5 average (more robust than single point)
                            era5_regional_mean = sum(era5_precip(:, :, month_idx)) / size(era5_precip(:, :, month_idx))
                            
                            if (era5_regional_mean > 0.01_real64 .and. &
                                cmip6_hist%precipitation(i, j, month_idx) > 0.01_real64) then
                                
                                valid_count = valid_count + 1
                                era5_monthly_mean = era5_monthly_mean + era5_regional_mean
                                cmip6_monthly_mean = cmip6_monthly_mean + cmip6_hist%precipitation(i, j, month_idx)
                            end if
                        end if
                    end do
                    
                    ! Compute precipitation bias factor
                    if (valid_count >= 10) then  ! Require at least 10 valid years
                        era5_monthly_mean = era5_monthly_mean / real(valid_count, real64)
                        cmip6_monthly_mean = cmip6_monthly_mean / real(valid_count, real64)
                        
                        if (cmip6_monthly_mean > 0.01_real64) then
                            bias_factors%precip_bias_factor(i, j, m) = era5_monthly_mean / cmip6_monthly_mean
                            ! Constrain bias factors to reasonable range
                            bias_factors%precip_bias_factor(i, j, m) = max(0.5_real64, &
                                min(2.0_real64, bias_factors%precip_bias_factor(i, j, m)))
                        end if
                    end if
                    
                    ! Compute PET bias factor (similar approach)
                    valid_count = 0
                    era5_monthly_mean = 0.0_real64
                    cmip6_monthly_mean = 0.0_real64
                    
                    do year_start_idx = 1, cmip6_hist%ntimes, 12
                        month_idx = year_start_idx + m - 1
                        if (month_idx <= cmip6_hist%ntimes .and. month_idx <= era5_overlap_end) then
                            ! Use regional ERA5 average for PET comparison too
                            era5_regional_mean = sum(era5_pet(:, :, month_idx)) / size(era5_pet(:, :, month_idx))
                            
                            if (era5_regional_mean > 0.01_real64 .and. &
                                cmip6_hist%pet_raw(i, j, month_idx) > 0.01_real64) then
                                
                                valid_count = valid_count + 1
                                era5_monthly_mean = era5_monthly_mean + era5_regional_mean
                                cmip6_monthly_mean = cmip6_monthly_mean + cmip6_hist%pet_raw(i, j, month_idx)
                            end if
                        end if
                    end do
                    
                    if (valid_count >= 10) then  ! Require at least 10 valid years
                        era5_monthly_mean = era5_monthly_mean / real(valid_count, real64)
                        cmip6_monthly_mean = cmip6_monthly_mean / real(valid_count, real64)
                        
                        if (cmip6_monthly_mean > 0.01_real64) then
                            bias_factors%pet_bias_factor(i, j, m) = era5_monthly_mean / cmip6_monthly_mean
                            ! Constrain bias factors to reasonable range
                            bias_factors%pet_bias_factor(i, j, m) = max(0.5_real64, &
                                min(2.0_real64, bias_factors%pet_bias_factor(i, j, m)))
                        end if
                    end if
                    
                end do
            end do
        end do
        
        ! Compute real correlation and RMSE statistics for all land grid cells
        do j = 1, cmip6_hist%nlats
            do i = 1, cmip6_hist%nlons
                
                ! Only process land points
                if (.not. cmip6_hist%land_mask(i, j)) cycle
                
                ! Calculate precipitation correlation and RMSE
                valid_count = 0
                era5_annual_mean = 0.0_real64
                cmip6_annual_mean = 0.0_real64
                
                ! First pass: calculate means
                do t = 1, min(cmip6_hist%ntimes, era5_overlap_end)
                    era5_regional_mean = sum(era5_precip(:, :, t)) / size(era5_precip(:, :, t))
                    if (era5_regional_mean > 0.0_real64 .and. cmip6_hist%precipitation(i, j, t) > 0.0_real64) then
                        era5_annual_mean = era5_annual_mean + era5_regional_mean
                        cmip6_annual_mean = cmip6_annual_mean + cmip6_hist%precipitation(i, j, t)
                        valid_count = valid_count + 1
                    end if
                end do
                
                if (valid_count > 10) then  ! Require sufficient data
                    era5_annual_mean = era5_annual_mean / real(valid_count, real64)
                    cmip6_annual_mean = cmip6_annual_mean / real(valid_count, real64)
                    
                    ! Second pass: calculate correlation and RMSE
                    correlation_sum = 0.0_real64
                    correlation_count = 0.0_real64
                    era5_monthly_var = 0.0_real64
                    cmip6_monthly_var = 0.0_real64
                    bias_factors%precip_rmse(i, j) = 0.0_real64
                    
                    do t = 1, min(cmip6_hist%ntimes, era5_overlap_end)
                        era5_regional_mean = sum(era5_precip(:, :, t)) / size(era5_precip(:, :, t))
                        if (era5_regional_mean > 0.0_real64 .and. cmip6_hist%precipitation(i, j, t) > 0.0_real64) then
                            ! Correlation calculation
                            correlation_sum = correlation_sum + &
                                (era5_regional_mean - era5_annual_mean) * (cmip6_hist%precipitation(i, j, t) - cmip6_annual_mean)
                            era5_monthly_var = era5_monthly_var + (era5_regional_mean - era5_annual_mean)**2
                            cmip6_monthly_var = cmip6_monthly_var + (cmip6_hist%precipitation(i, j, t) - cmip6_annual_mean)**2
                            
                            ! RMSE calculation
                            bias_factors%precip_rmse(i, j) = bias_factors%precip_rmse(i, j) + &
                                (era5_regional_mean - cmip6_hist%precipitation(i, j, t))**2
                            correlation_count = correlation_count + 1.0_real64
                        end if
                    end do
                    
                    ! Finalize correlation
                    if (era5_monthly_var > 0.0_real64 .and. cmip6_monthly_var > 0.0_real64) then
                        bias_factors%precip_correlation(i, j) = correlation_sum / sqrt(era5_monthly_var * cmip6_monthly_var)
                        bias_factors%precip_correlation(i, j) = max(-1.0_real64, min(1.0_real64, bias_factors%precip_correlation(i, j)))
                    else
                        bias_factors%precip_correlation(i, j) = 0.0_real64
                    end if
                    
                    ! Finalize RMSE
                    if (correlation_count > 0.0_real64) then
                        bias_factors%precip_rmse(i, j) = sqrt(bias_factors%precip_rmse(i, j) / correlation_count)
                    end if
                end if
                
                ! Calculate PET correlation and RMSE
                valid_count = 0
                era5_annual_mean = 0.0_real64
                cmip6_annual_mean = 0.0_real64
                
                ! First pass: calculate PET means
                do t = 1, min(cmip6_hist%ntimes, era5_overlap_end)
                    era5_regional_mean = sum(era5_pet(:, :, t)) / size(era5_pet(:, :, t))
                    if (era5_regional_mean > 0.0_real64 .and. cmip6_hist%pet_raw(i, j, t) > 0.0_real64) then
                        era5_annual_mean = era5_annual_mean + era5_regional_mean
                        cmip6_annual_mean = cmip6_annual_mean + cmip6_hist%pet_raw(i, j, t)
                        valid_count = valid_count + 1
                    end if
                end do
                
                if (valid_count > 10) then  ! Require sufficient data
                    era5_annual_mean = era5_annual_mean / real(valid_count, real64)
                    cmip6_annual_mean = cmip6_annual_mean / real(valid_count, real64)
                    
                    ! Second pass: calculate PET correlation and RMSE
                    correlation_sum = 0.0_real64
                    correlation_count = 0.0_real64
                    era5_monthly_var = 0.0_real64
                    cmip6_monthly_var = 0.0_real64
                    bias_factors%pet_rmse(i, j) = 0.0_real64
                    
                    do t = 1, min(cmip6_hist%ntimes, era5_overlap_end)
                        era5_regional_mean = sum(era5_pet(:, :, t)) / size(era5_pet(:, :, t))
                        if (era5_regional_mean > 0.0_real64 .and. cmip6_hist%pet_raw(i, j, t) > 0.0_real64) then
                            ! PET correlation calculation
                            correlation_sum = correlation_sum + &
                                (era5_regional_mean - era5_annual_mean) * (cmip6_hist%pet_raw(i, j, t) - cmip6_annual_mean)
                            era5_monthly_var = era5_monthly_var + (era5_regional_mean - era5_annual_mean)**2
                            cmip6_monthly_var = cmip6_monthly_var + (cmip6_hist%pet_raw(i, j, t) - cmip6_annual_mean)**2
                            
                            ! PET RMSE calculation
                            bias_factors%pet_rmse(i, j) = bias_factors%pet_rmse(i, j) + &
                                (era5_regional_mean - cmip6_hist%pet_raw(i, j, t))**2
                            correlation_count = correlation_count + 1.0_real64
                        end if
                    end do
                    
                    ! Finalize PET correlation
                    if (era5_monthly_var > 0.0_real64 .and. cmip6_monthly_var > 0.0_real64) then
                        bias_factors%pet_correlation(i, j) = correlation_sum / sqrt(era5_monthly_var * cmip6_monthly_var)
                        bias_factors%pet_correlation(i, j) = max(-1.0_real64, min(1.0_real64, bias_factors%pet_correlation(i, j)))
                    else
                        bias_factors%pet_correlation(i, j) = 0.0_real64
                    end if
                    
                    ! Finalize PET RMSE
                    if (correlation_count > 0.0_real64) then
                        bias_factors%pet_rmse(i, j) = sqrt(bias_factors%pet_rmse(i, j) / correlation_count)
                    end if
                end if
                
            end do
        end do
        
        ! Clean up
        deallocate(era5_month_values, cmip6_month_values, era5_annual_values, cmip6_annual_values)
        
        bias_factors%is_computed = .true.
        
        write(*,*) "✅ Real bias correction factors computed from data statistics"
        write(*,*) "   Sample precipitation bias factors (month 1): ", &
                   minval(bias_factors%precip_bias_factor(:,:,1)), " to ", &
                   maxval(bias_factors%precip_bias_factor(:,:,1))
        write(*,*) "   Sample PET bias factors (month 1): ", &
                   minval(bias_factors%pet_bias_factor(:,:,1)), " to ", &
                   maxval(bias_factors%pet_bias_factor(:,:,1))
        
        ! ===================================================================
        ! SUMMARY STATISTICS FOR THESIS REPORT
        ! ===================================================================
        block
            real(real64) :: mean_precip_bias, mean_pet_bias, mean_r_precip, mean_r_pet, mean_rmse_precip, mean_rmse_pet
            integer :: valid_precip_count, valid_pet_count, valid_r_precip_count, valid_r_pet_count
            integer :: valid_rmse_precip_count, valid_rmse_pet_count
            
            ! Calculate mean bias factors (only valid land points)
            valid_precip_count = count(bias_factors%precip_bias_factor > 0.0_real64 .and. &
                                     bias_factors%precip_bias_factor < 10.0_real64)
            valid_pet_count = count(bias_factors%pet_bias_factor > 0.0_real64 .and. &
                                  bias_factors%pet_bias_factor < 10.0_real64)
            
            if (valid_precip_count > 0) then
                mean_precip_bias = sum(bias_factors%precip_bias_factor, &
                                     mask=(bias_factors%precip_bias_factor > 0.0_real64 .and. &
                                           bias_factors%precip_bias_factor < 10.0_real64)) / real(valid_precip_count, real64)
            else
                mean_precip_bias = 1.0_real64
            end if
            
            if (valid_pet_count > 0) then
                mean_pet_bias = sum(bias_factors%pet_bias_factor, &
                                  mask=(bias_factors%pet_bias_factor > 0.0_real64 .and. &
                                        bias_factors%pet_bias_factor < 10.0_real64)) / real(valid_pet_count, real64)
            else
                mean_pet_bias = 1.0_real64
            end if
            
            ! Calculate mean correlations (only valid correlations)
            valid_r_precip_count = count(abs(bias_factors%precip_correlation) > 0.01_real64 .and. &
                                       abs(bias_factors%precip_correlation) <= 1.0_real64)
            valid_r_pet_count = count(abs(bias_factors%pet_correlation) > 0.01_real64 .and. &
                                    abs(bias_factors%pet_correlation) <= 1.0_real64)
            
            if (valid_r_precip_count > 0) then
                mean_r_precip = sum(bias_factors%precip_correlation, &
                                  mask=(abs(bias_factors%precip_correlation) > 0.01_real64 .and. &
                                        abs(bias_factors%precip_correlation) <= 1.0_real64)) / real(valid_r_precip_count, real64)
            else
                mean_r_precip = 0.0_real64
            end if
            
            if (valid_r_pet_count > 0) then
                mean_r_pet = sum(bias_factors%pet_correlation, &
                               mask=(abs(bias_factors%pet_correlation) > 0.01_real64 .and. &
                                     abs(bias_factors%pet_correlation) <= 1.0_real64)) / real(valid_r_pet_count, real64)
            else
                mean_r_pet = 0.0_real64
            end if
            
            ! Calculate mean RMSE (only valid RMSE values)
            valid_rmse_precip_count = count(bias_factors%precip_rmse > 0.0_real64 .and. &
                                          bias_factors%precip_rmse < 1000.0_real64)
            valid_rmse_pet_count = count(bias_factors%pet_rmse > 0.0_real64 .and. &
                                       bias_factors%pet_rmse < 1000.0_real64)
            
            if (valid_rmse_precip_count > 0) then
                mean_rmse_precip = sum(bias_factors%precip_rmse, &
                                     mask=(bias_factors%precip_rmse > 0.0_real64 .and. &
                                           bias_factors%precip_rmse < 1000.0_real64)) / real(valid_rmse_precip_count, real64)
            else
                mean_rmse_precip = 0.0_real64
            end if
            
            if (valid_rmse_pet_count > 0) then
                mean_rmse_pet = sum(bias_factors%pet_rmse, &
                                  mask=(bias_factors%pet_rmse > 0.0_real64 .and. &
                                        bias_factors%pet_rmse < 1000.0_real64)) / real(valid_rmse_pet_count, real64)
            else
                mean_rmse_pet = 0.0_real64
            end if
            
            write(*,*) ""
            write(*,*) "=========================================="
            write(*,*) "  BIAS CORRECTION SUMMARY (1981–2014)"
            write(*,*) "=========================================="
            write(*,*) "  Mean Precip Bias Factor: ", mean_precip_bias
            write(*,*) "  Mean PET Bias Factor:    ", mean_pet_bias
            write(*,*) "  Mean Correlation P/PET:  ", mean_r_precip, "/", mean_r_pet
            write(*,*) "  Mean RMSE P/PET (mm):    ", mean_rmse_precip, "/", mean_rmse_pet
            write(*,*) "  Valid grid cells (P/PET):", valid_precip_count, "/", valid_pet_count
            write(*,*) "=========================================="
            write(*,*) ""
            
            ! Save summary table for thesis reporting
            call save_bias_correction_summary_table(mean_precip_bias, mean_pet_bias, &
                                                   mean_r_precip, mean_r_pet, &
                                                   mean_rmse_precip, mean_rmse_pet, &
                                                   valid_precip_count, valid_pet_count)
            
            ! Save summary statistics to table file for thesis
            call save_bias_correction_summary_table(mean_precip_bias, mean_pet_bias, &
                                                   mean_r_precip, mean_r_pet, &
                                                   mean_rmse_precip, mean_rmse_pet, &
                                                   valid_precip_count, valid_pet_count)
        end block
        
        success = .true.
        
    end subroutine compute_bias_correction_factors_correct

    subroutine apply_precipitation_pet_bias_correction(cmip6_data, bias_factors, success)
        !===========================================================================
        ! Apply bias correction to precipitation and PET directly (not drivers)
        !===========================================================================
        type(cmip6_scenario_data_t), intent(inout) :: cmip6_data
        type(bias_correction_factors_t), intent(in) :: bias_factors
        logical, intent(out) :: success
        
        integer :: i, j, k, month
        
        success = .false.
        
        if (.not. bias_factors%is_computed) then
            write(*,*) "❌ Bias correction factors not computed"
            return
        end if
        
        if (.not. allocated(cmip6_data%pet_corrected)) then
            write(*,*) "❌ CMIP6 PET not calculated yet - call calculate_hargreaves_pet_from_cmip6_drivers first"
            return
        end if
        
        write(*,*) "🔧 Applying P and PET bias correction to ", trim(cmip6_data%scenario), "..."
        
        ! Apply bias correction to each time step
        do k = 1, cmip6_data%ntimes
            month = mod(k-1, 12) + 1  ! Get month (1-12)
            
            do j = 1, cmip6_data%nlats
                do i = 1, cmip6_data%nlons
                    
                    ! Correct precipitation (multiplicative)
                    cmip6_data%precipitation_corrected(i,j,k) = cmip6_data%precipitation(i,j,k) * &
                                                               bias_factors%precip_bias_factor(i,j,month)
                    cmip6_data%precipitation_corrected(i,j,k) = max(0.0_real64, cmip6_data%precipitation_corrected(i,j,k))
                    
                    ! Clean up tiny values (< 1e-6 mm/month)
                    if (cmip6_data%precipitation_corrected(i,j,k) < 1e-6_real64) then
                        cmip6_data%precipitation_corrected(i,j,k) = 0.0_real64
                    end if
                    
                    ! Correct PET (multiplicative) - using CMIP6-calculated PET raw data
                    cmip6_data%pet_corrected(i,j,k) = cmip6_data%pet_raw(i,j,k) * &
                                                      bias_factors%pet_bias_factor(i,j,month)
                    cmip6_data%pet_corrected(i,j,k) = max(0.0_real64, cmip6_data%pet_corrected(i,j,k))
                    
                end do
            end do
        end do
        
        cmip6_data%is_corrected = .true.
        
        write(*,*) "✅ P and PET bias correction applied to ", trim(cmip6_data%scenario)
        write(*,*) "   Corrected precipitation range: ", minval(cmip6_data%precipitation_corrected), " to ", &
                   maxval(cmip6_data%precipitation_corrected), " mm/month"
        write(*,*) "   Corrected PET range: ", minval(cmip6_data%pet_corrected), " to ", &
                   maxval(cmip6_data%pet_corrected), " mm/month"
        
        success = .true.
        
    end subroutine apply_precipitation_pet_bias_correction

    subroutine calculate_hargreaves_pet_from_cmip6_drivers(cmip6_data, success)
        !===========================================================================
        ! Calculate Hargreaves PET from CMIP6 drivers with variable solar radiation
        ! Similar to ERA5 approach but using CMIP6 temperature drivers
        !===========================================================================
        type(cmip6_scenario_data_t), intent(inout) :: cmip6_data
        logical, intent(out), optional :: success
        
        real(real64), allocatable :: tmean(:,:,:), trange(:,:,:)
        real(real64) :: Ra_lat_month, lat_rad, day_of_year, solar_declination, &
                        sunset_hour_angle, dr, lat_deg, tmean_val, trange_val
        real(real64), parameter :: hargreaves_coeff = 0.0023_real64
        real(real64), parameter :: PI = 3.14159265359_real64
        real(real64), parameter :: DAYS_PER_MONTH = 30.44_real64  ! From ERA5 processing
        real(real64), parameter :: SOLAR_CONSTANT = 0.0820_real64  ! MJ m-2 day-1
        integer :: i, j, k, month, year_offset, current_year, current_month
        logical :: local_success
        
        local_success = .false.
        
        write(*,*) "🌡️  Calculating Hargreaves PET with variable Ra from CMIP6 drivers for ", trim(cmip6_data%scenario), "..."
        
        ! Allocate temporary arrays
        allocate(tmean(cmip6_data%nlons, cmip6_data%nlats, cmip6_data%ntimes))
        allocate(trange(cmip6_data%nlons, cmip6_data%nlats, cmip6_data%ntimes))
        
        ! Make sure pet_raw array is allocated (will store CMIP6-calculated PET before bias correction)
        if (.not. allocated(cmip6_data%pet_raw)) then
            allocate(cmip6_data%pet_raw(cmip6_data%nlons, cmip6_data%nlats, cmip6_data%ntimes))
        end if
        
        ! Convert from Kelvin to Celsius and calculate components
        tmean = (cmip6_data%tasmax + cmip6_data%tasmin) / 2.0_real64 - KELVIN_TO_CELSIUS
        trange = sqrt(abs(cmip6_data%tasmax - cmip6_data%tasmin))  ! Ensure positive
        
        ! Calculate PET using Hargreaves formula with variable Ra
        do k = 1, cmip6_data%ntimes
            ! Determine the month and year for this time step
            current_month = mod(k - 1, 12) + 1
            current_year = cmip6_data%start_year + (k - 1) / 12
            
            ! Calculate representative day of year for this month (middle of month)
            day_of_year = real(current_month - 1, real64) * DAYS_PER_MONTH + DAYS_PER_MONTH/2.0_real64
            
            do j = 1, cmip6_data%nlats
                ! Get latitude in degrees
                lat_deg = cmip6_data%latitude(j)
                lat_rad = lat_deg * PI / 180.0_real64
                
                ! Calculate solar radiation components
                solar_declination = 0.4093_real64 * sin((2.0_real64 * PI * day_of_year / 365.0_real64) - 1.405_real64)
                sunset_hour_angle = acos(-tan(lat_rad) * tan(solar_declination))
                dr = 1.0_real64 + 0.033_real64 * cos(2.0_real64 * PI * day_of_year / 365.0_real64)
                
                ! Calculate extraterrestrial radiation (Ra) in mm/day
                Ra_lat_month = (24.0_real64 * 60.0_real64 / PI) * SOLAR_CONSTANT * dr * &
                              (sunset_hour_angle * sin(lat_rad) * sin(solar_declination) + &
                               cos(lat_rad) * cos(solar_declination) * sin(sunset_hour_angle))
                
                ! Convert to mm/month using ERA5 approach
                Ra_lat_month = Ra_lat_month * DAYS_PER_MONTH
                
                do i = 1, cmip6_data%nlons
                    ! Only calculate for land points
                    if (cmip6_data%land_mask(i, j)) then
                        tmean_val = tmean(i, j, k)
                        trange_val = trange(i, j, k)
                        
                        ! Check for valid temperature data
                        if (tmean_val > -50.0_real64 .and. tmean_val < 60.0_real64 .and. trange_val >= 0.0_real64) then
                            
                            ! Calculate PET using Hargreaves formula
                            cmip6_data%pet_raw(i, j, k) = hargreaves_coeff * Ra_lat_month * &
                                                         (tmean_val + 17.8_real64) * trange_val
                            
                            ! Ensure non-negative PET
                            if (cmip6_data%pet_raw(i, j, k) < 0.0_real64) then
                                cmip6_data%pet_raw(i, j, k) = 0.0_real64
                            end if
                        else
                            ! Set invalid data to zero for now
                            cmip6_data%pet_raw(i, j, k) = 0.0_real64
                        end if
                    else
                        ! Ocean points
                        cmip6_data%pet_raw(i, j, k) = 0.0_real64
                    end if
                end do
            end do
        end do
        
        ! Initialize pet_corrected to pet_raw (will be modified during bias correction)
        cmip6_data%pet_corrected = cmip6_data%pet_raw
        
        ! Cleanup
        deallocate(tmean, trange)
        
        write(*,*) "✅ Enhanced Hargreaves PET calculated with variable Ra for ", trim(cmip6_data%scenario)
        write(*,*) "   CMIP6-calculated PET range: ", minval(cmip6_data%pet_raw), &
                   " to ", maxval(cmip6_data%pet_raw), " mm/month"
        write(*,*) "   Valid PET array allocated with dimensions: ", shape(cmip6_data%pet_raw)
        
        local_success = .true.
        if (present(success)) success = local_success
        
    end subroutine calculate_hargreaves_pet_from_cmip6_drivers

    subroutine save_bias_correction_summary_table(mean_precip_bias, mean_pet_bias, &
                                                 mean_r_precip, mean_r_pet, &
                                                 mean_rmse_precip, mean_rmse_pet, &
                                                 valid_precip_count, valid_pet_count)
        !===========================================================================
        ! Save bias correction summary statistics to table file for thesis
        !===========================================================================
        real(real64), intent(in) :: mean_precip_bias, mean_pet_bias
        real(real64), intent(in) :: mean_r_precip, mean_r_pet
        real(real64), intent(in) :: mean_rmse_precip, mean_rmse_pet
        integer, intent(in) :: valid_precip_count, valid_pet_count
        
        character(len=512) :: table_file
        integer :: unit_num, status
        character(len=20) :: date_str, time_str
        
        ! Create output file path
        table_file = "data/Final_results/Tables/bias_correction_summary_1981_2014.csv"
        
        ! Create directory if it doesn't exist
        call system("mkdir -p data/Final_results/Tables")
        
        ! Get current date and time
        call date_and_time(date_str, time_str)
        
        ! Open file for writing
        open(newunit=unit_num, file=trim(table_file), status='replace', iostat=status)
        
        if (status /= 0) then
            write(*,*) "⚠ Warning: Could not save bias correction table to ", trim(table_file)
            return
        end if
        
        ! Write CSV header
        write(unit_num, '(A)') "# CMIP6 Bias Correction Summary Statistics (1981-2014)"
        write(unit_num, '(A)') "# Generated on " // date_str(1:4) // "-" // date_str(5:6) // "-" // date_str(7:8) // &
                              " at " // time_str(1:2) // ":" // time_str(3:4) // ":" // time_str(5:6)
        write(unit_num, '(A)') "# Reference: ERA5-Land vs CMIP6 MPI-ESM1-2-LR historical"
        write(unit_num, '(A)') "# Methodology: Regional averaging with constrained bias factors (0.5-2.0x)"
        write(unit_num, '(A)') ""
        write(unit_num, '(A)') "Variable,Mean_Bias_Factor,Mean_Correlation,Mean_RMSE_mm,Valid_Grid_Cells,Performance_Rating"
        
        ! Write precipitation statistics
        write(unit_num, '(A,F8.4,A,F6.3,A,F8.2,A,I0,A,A)') &
            "Precipitation,", mean_precip_bias, ",", mean_r_precip, ",", mean_rmse_precip, ",", &
            valid_precip_count, ",", trim(assess_performance(mean_r_precip, mean_precip_bias))
        
        ! Write PET statistics
        write(unit_num, '(A,F8.4,A,F6.3,A,F8.2,A,I0,A,A)') &
            "PET,", mean_pet_bias, ",", mean_r_pet, ",", mean_rmse_pet, ",", &
            valid_pet_count, ",", trim(assess_performance(mean_r_pet, mean_pet_bias))
        
        ! Write summary information
        write(unit_num, '(A)') ""
        write(unit_num, '(A)') "# Summary Statistics"
        write(unit_num, '(A,F6.3)') "# Overall Mean Correlation (P & PET): ", (mean_r_precip + mean_r_pet) / 2.0_real64
        write(unit_num, '(A,F6.3)') "# Overall Mean Bias Factor: ", (mean_precip_bias + mean_pet_bias) / 2.0_real64
        write(unit_num, '(A,I0)') "# Total Valid Grid Cells: ", valid_precip_count + valid_pet_count
        write(unit_num, '(A)') "# "
        write(unit_num, '(A)') "# Performance Rating:"
        write(unit_num, '(A)') "#   Excellent: r > 0.8, bias 0.9-1.1"
        write(unit_num, '(A)') "#   Good:      r > 0.6, bias 0.8-1.2"
        write(unit_num, '(A)') "#   Fair:      r > 0.4, bias 0.7-1.3"
        write(unit_num, '(A)') "#   Poor:      r < 0.4, bias < 0.7 or > 1.3"
        
        close(unit_num)
        
        write(*,*) "📊 Bias correction summary table saved to:"
        write(*,*) "   ", trim(table_file)
        
    end subroutine save_bias_correction_summary_table
    
    function assess_performance(correlation, bias_factor) result(rating)
        !===========================================================================
        ! Assess bias correction performance based on correlation and bias factor
        !===========================================================================
        real(real64), intent(in) :: correlation, bias_factor
        character(len=10) :: rating
        
        if (abs(correlation) > 0.8_real64 .and. bias_factor >= 0.9_real64 .and. bias_factor <= 1.1_real64) then
            rating = "Excellent"
        else if (abs(correlation) > 0.6_real64 .and. bias_factor >= 0.8_real64 .and. bias_factor <= 1.2_real64) then
            rating = "Good"
        else if (abs(correlation) > 0.4_real64 .and. bias_factor >= 0.7_real64 .and. bias_factor <= 1.3_real64) then
            rating = "Fair"
        else
            rating = "Poor"
        end if
        
    end function assess_performance

    !===========================================================================
    ! PACKAGE_CORRECTED_SCENARIOS: Memory-based data packaging
    !===========================================================================
    subroutine package_corrected_scenarios(ssp126_data, ssp245_data, ssp585_data, &
                                          corrected_scenarios, success)
        !===========================================================================
        ! Package bias-corrected data into structured output for memory-based processing
        !===========================================================================
        type(cmip6_scenario_data_t), intent(in) :: ssp126_data, ssp245_data, ssp585_data
        type(corrected_cmip6_scenarios_t), intent(out) :: corrected_scenarios
        logical, intent(out) :: success
        
        integer :: nlon, nlat, ntime
        
        success = .false.
        
        ! Verify all scenarios have consistent dimensions
        if (ssp126_data%nlons /= ssp245_data%nlons .or. ssp245_data%nlons /= ssp585_data%nlons .or. &
            ssp126_data%nlats /= ssp245_data%nlats .or. ssp245_data%nlats /= ssp585_data%nlats .or. &
            ssp126_data%ntimes /= ssp245_data%ntimes .or. ssp245_data%ntimes /= ssp585_data%ntimes) then
            print *, "❌ Inconsistent dimensions across CMIP6 scenarios"
            return
        end if
        
        ! Get dimensions from SSP126 (all should be the same)
        nlon = ssp126_data%nlons
        nlat = ssp126_data%nlats
        ntime = ssp126_data%ntimes
        
        ! Allocate arrays in structured output
        allocate(corrected_scenarios%ssp126_precip(nlon, nlat, ntime))
        allocate(corrected_scenarios%ssp126_pet(nlon, nlat, ntime))
        allocate(corrected_scenarios%ssp245_precip(nlon, nlat, ntime))
        allocate(corrected_scenarios%ssp245_pet(nlon, nlat, ntime))
        allocate(corrected_scenarios%ssp585_precip(nlon, nlat, ntime))
        allocate(corrected_scenarios%ssp585_pet(nlon, nlat, ntime))
        
        allocate(corrected_scenarios%latitudes(nlat))
        allocate(corrected_scenarios%longitudes(nlon))
        allocate(corrected_scenarios%time_values(ntime))
        
        ! Copy bias-corrected data
        corrected_scenarios%ssp126_precip = ssp126_data%precipitation_corrected
        corrected_scenarios%ssp126_pet = ssp126_data%pet_corrected
        corrected_scenarios%ssp245_precip = ssp245_data%precipitation_corrected
        corrected_scenarios%ssp245_pet = ssp245_data%pet_corrected
        corrected_scenarios%ssp585_precip = ssp585_data%precipitation_corrected
        corrected_scenarios%ssp585_pet = ssp585_data%pet_corrected
        
        ! Copy coordinate arrays (from SSP126, all should be identical)
        corrected_scenarios%latitudes = ssp126_data%latitude
        corrected_scenarios%longitudes = ssp126_data%longitude
        corrected_scenarios%time_values = ssp126_data%time_values
        
        ! Set metadata
        corrected_scenarios%start_year = ssp126_data%start_year
        corrected_scenarios%end_year = ssp126_data%end_year
        corrected_scenarios%data_available = .true.
        
        success = .true.
        print *, "✅ Bias-corrected scenarios packaged into structured data"
        print *, "   Dimensions: ", nlon, "x", nlat, "x", ntime
        print *, "   Period: ", corrected_scenarios%start_year, "-", corrected_scenarios%end_year
        
    end subroutine package_corrected_scenarios

end module cmip6_bias_correction
