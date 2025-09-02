module bias_correction_clean
    !===============================================================================
    ! STREAMLINED QUANTILE MAPPING BIAS CORRECTION MODULE
    !===============================================================================
    ! Purpose: Statistical bias correction for CMIP6 climate data using quantile mapping
    ! Method: Distribution-based quantile mapping (Themeßl et al., 2011; Cannon et al., 2015)
    ! 
    ! Scientific Approach:
    !   1. Fit probability distributions to historical observations (ERA5) and model data (CMIP6)
    !   2. For each future CMIP6 value: 
    !      a) Calculate quantile in CMIP6 distribution: q = F_model(x)
    !      b) Map to observation distribution: x_corrected = F_obs^(-1)(q)
    !   3. Precipitation: Gamma distribution (shape, scale parameters)
    !   4. PET: Normal distribution (mean, standard deviation)
    !
    ! References:
    !   - Themeßl et al. (2011): Empirical-statistical downscaling and error correction
    !   - Cannon et al. (2015): Bias correction of GCM precipitation by quantile mapping
    !
    ! Author: Khadar Daahir (University of Glasgow)
    ! Date: August 2025
    !===============================================================================
    
    use iso_fortran_env, only: real64, output_unit, error_unit
    use ieee_arithmetic, only: ieee_is_finite
    use io_module, only: cmip6_scenario_data_t, corrected_cmip6_scenarios_t, FILL_VALUE
    use fsml, only: fsml_mean, fsml_std, fsml_var, fsml_pcc
    use fsml_dst, only: f_dst_gamma_cdf_core, f_dst_gamma_ppf_core, &
                        f_dst_norm_cdf_core, f_dst_norm_ppf_core
    implicit none
    
    private
    public :: bias_correction_parameters_t, apply_quantile_mapping_correction
    public :: process_bias_correction, calculate_hargreaves_pet
    public :: generate_validation_tables, bias_correction_validation_t
    public :: generate_bias_correction_validation_tables
    public :: generate_spei12_comparison_files
    
    ! Clean data structure for bias correction parameters
    type :: bias_correction_parameters_t
        ! Monthly distribution parameters [lon, lat, month]
        real(real64), allocatable :: era5_precip_shape(:,:,:)     ! Gamma shape
        real(real64), allocatable :: era5_precip_scale(:,:,:)     ! Gamma scale
        real(real64), allocatable :: cmip6_precip_shape(:,:,:)    ! Gamma shape
        real(real64), allocatable :: cmip6_precip_scale(:,:,:)    ! Gamma scale
        
        real(real64), allocatable :: era5_pet_mean(:,:,:)         ! Normal mean
        real(real64), allocatable :: era5_pet_std(:,:,:)          ! Normal std
        real(real64), allocatable :: cmip6_pet_mean(:,:,:)        ! Normal mean
        real(real64), allocatable :: cmip6_pet_std(:,:,:)         ! Normal std
        
        integer :: nlons, nlats
        logical :: fitted = .false.
    end type bias_correction_parameters_t
    
    ! Validation metrics type for table generation
    type :: bias_correction_validation_t
        real(real64) :: bias_before, bias_after               ! Difference bias (mm/month)
        real(real64) :: bias_ratio_before, bias_ratio_after   ! Ratio bias (model/obs)
        real(real64) :: rmse_before, rmse_after
        real(real64) :: correlation_before, correlation_after
        real(real64) :: mean_obs, mean_corrected, mean_original
        real(real64) :: std_obs, std_corrected, std_original
    end type bias_correction_validation_t
    
    ! Constants
    real(real64), parameter :: MIN_PRECIP = 0.01_real64
    real(real64), parameter :: KELVIN_TO_CELSIUS = 273.15_real64
    real(real64), parameter :: PI = 3.14159265359_real64
    
    ! Historical period for parameter fitting
    integer, parameter :: HIST_START_YEAR = 1981
    integer, parameter :: HIST_END_YEAR = 2014
    
contains
    
    !===============================================================================
    ! Calculate bias correction validation metrics using FSML statistical functions
    !===============================================================================
    function calculate_bias_metrics(obs_data, original_data, corrected_data) result(metrics)
        implicit none
        real(real64), intent(in) :: obs_data(:)
        real(real64), intent(in) :: original_data(:)
        real(real64), intent(in) :: corrected_data(:)
        type(bias_correction_validation_t) :: metrics
        
        ! Calculate means using FSML
        metrics%mean_obs = fsml_mean(obs_data)
        metrics%mean_original = fsml_mean(original_data)
        metrics%mean_corrected = fsml_mean(corrected_data)
        
        ! Calculate standard deviations using FSML
        metrics%std_obs = fsml_std(obs_data)
        metrics%std_original = fsml_std(original_data)
        metrics%std_corrected = fsml_std(corrected_data)
        
        ! Calculate biases as differences (model - observations) in mm/month
        ! This is the standard for precipitation bias correction studies
        metrics%bias_before = metrics%mean_original - metrics%mean_obs
        metrics%bias_after = metrics%mean_corrected - metrics%mean_obs
        
        ! Calculate bias ratios (model / observations) for supplementary context
        ! Avoid division by zero by checking for very small values
        if (abs(metrics%mean_obs) > 1.0e-10_real64) then
            metrics%bias_ratio_before = metrics%mean_original / metrics%mean_obs
            metrics%bias_ratio_after = metrics%mean_corrected / metrics%mean_obs
        else
            metrics%bias_ratio_before = 1.0_real64  ! No bias if observations are effectively zero
            metrics%bias_ratio_after = 1.0_real64
        end if
        
        ! Calculate RMSE
        metrics%rmse_before = sqrt(fsml_mean((original_data - obs_data)**2))
        metrics%rmse_after = sqrt(fsml_mean((corrected_data - obs_data)**2))
        
        ! Calculate correlations using FSML
        metrics%correlation_before = fsml_pcc(obs_data, original_data)
        metrics%correlation_after = fsml_pcc(obs_data, corrected_data)
        
    end function calculate_bias_metrics

    subroutine process_bias_correction(era5_precip, era5_pet, &
                                     historical_data, ssp126_data, ssp245_data, ssp585_data, &
                                     corrected_scenarios, success)
        !===========================================================================
        ! Main bias correction workflow - clean and efficient
        !===========================================================================
        real(real64), intent(in) :: era5_precip(:,:,:), era5_pet(:,:,:)
        type(cmip6_scenario_data_t), intent(inout) :: historical_data, ssp126_data, ssp245_data, ssp585_data
        type(corrected_cmip6_scenarios_t), intent(out) :: corrected_scenarios
        logical, intent(out) :: success
        
        type(bias_correction_parameters_t) :: bc_params
        
        success = .false.
        
        write(output_unit, '(A)') '🔧 Starting quantile mapping bias correction...'
        
        ! Step 1: Calculate PET for all CMIP6 scenarios
        call calculate_hargreaves_pet(historical_data, success)
        if (.not. success) return
        
        call calculate_hargreaves_pet(ssp126_data, success)
        if (.not. success) return
        
        call calculate_hargreaves_pet(ssp245_data, success)
        if (.not. success) return
        
        call calculate_hargreaves_pet(ssp585_data, success)
        if (.not. success) return
        
        ! Step 2: Fit distribution parameters using historical period
        call fit_distribution_parameters(era5_precip, era5_pet, historical_data, bc_params, success)
        if (.not. success) return
        
        ! Step 3: Apply quantile mapping to all scenarios
        call apply_quantile_mapping_correction(historical_data, bc_params, success)
        if (.not. success) return
        
        call apply_quantile_mapping_correction(ssp126_data, bc_params, success)
        if (.not. success) return
        
        call apply_quantile_mapping_correction(ssp245_data, bc_params, success)
        if (.not. success) return
        
        call apply_quantile_mapping_correction(ssp585_data, bc_params, success)
        if (.not. success) return
        
        ! Step 4: Package results
        call package_corrected_scenarios(ssp126_data, ssp245_data, ssp585_data, corrected_scenarios, success)
        
        write(output_unit, '(A)') '✅ Bias correction completed successfully'
        
    end subroutine process_bias_correction

    subroutine calculate_hargreaves_pet(cmip6_data, success)
        !===========================================================================
        ! Calculate PET using Hargreaves-Samani method (FAO-56 based)
        !===========================================================================
        type(cmip6_scenario_data_t), intent(inout) :: cmip6_data
        logical, intent(out) :: success
        
        integer :: i, j, k
        real(real64) :: tmax_c, tmin_c, tmean_c, temp_range
        real(real64) :: lat_rad, Ra, pet_value
        
        success = .false.
        
        if (.not. allocated(cmip6_data%pet_raw)) then
            allocate(cmip6_data%pet_raw(cmip6_data%nlons, cmip6_data%nlats, cmip6_data%ntimes))
        end if
        
        ! Calculate PET for each grid point and time step
        do k = 1, cmip6_data%ntimes
            do j = 1, cmip6_data%nlats
                lat_rad = cmip6_data%latitude(j) * PI / 180.0_real64
                Ra = calculate_extraterrestrial_radiation(lat_rad, k)
                
                do i = 1, cmip6_data%nlons
                    if (cmip6_data%land_mask(i,j)) then
                        if (cmip6_data%tasmax(i,j,k) /= FILL_VALUE .and. &
                            cmip6_data%tasmin(i,j,k) /= FILL_VALUE .and. &
                            cmip6_data%tasmax(i,j,k) > 200.0_real64 .and. &
                            cmip6_data%tasmin(i,j,k) > 200.0_real64) then
                            
                            ! Convert to Celsius
                            tmax_c = cmip6_data%tasmax(i,j,k) - KELVIN_TO_CELSIUS
                            tmin_c = cmip6_data%tasmin(i,j,k) - KELVIN_TO_CELSIUS
                            tmean_c = (tmax_c + tmin_c) / 2.0_real64
                            temp_range = max(1.0_real64, tmax_c - tmin_c)
                            
                            ! Hargreaves-Samani formula: PET = 0.0023 * Ra * (Tmean + 17.8) * sqrt(DTR)
                            pet_value = 0.0023_real64 * Ra * (tmean_c + 17.8_real64) * sqrt(temp_range)
                            pet_value = max(0.1_real64, min(1000.0_real64, pet_value))  ! Reasonable bounds
                            
                            cmip6_data%pet_raw(i,j,k) = pet_value
                        else
                            cmip6_data%pet_raw(i,j,k) = FILL_VALUE
                        end if
                    else
                        cmip6_data%pet_raw(i,j,k) = FILL_VALUE
                    end if
                end do
            end do
        end do
        
        success = .true.
        
    end subroutine calculate_hargreaves_pet

    function calculate_extraterrestrial_radiation(lat_rad, time_index) result(Ra)
        !===========================================================================
        ! Calculate extraterrestrial radiation using FAO-56 method
        !===========================================================================
        real(real64), intent(in) :: lat_rad
        integer, intent(in) :: time_index
        real(real64) :: Ra
        
        real(real64) :: day_of_year, dr, delta, omega_s
        real(real64), parameter :: GSC = 0.0820_real64  ! Solar constant [MJ m-2 min-1]
        
        ! Approximate day of year from monthly time index
        day_of_year = real(mod(time_index-1, 12) * 30 + 15, real64)
        
        ! Inverse relative distance Earth-Sun
        dr = 1.0_real64 + 0.033_real64 * cos(2.0_real64 * PI * day_of_year / 365.0_real64)
        
        ! Solar declination
        delta = 0.409_real64 * sin(2.0_real64 * PI * day_of_year / 365.0_real64 - 1.39_real64)
        
        ! Sunset hour angle
        omega_s = acos(-tan(lat_rad) * tan(delta))
        
        ! Extraterrestrial radiation [MJ m-2 day-1]
        Ra = (24.0_real64 * 60.0_real64 / PI) * GSC * dr * &
             (omega_s * sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * sin(omega_s))
        
        ! Ensure reasonable bounds for Horn of Africa latitude
        Ra = max(15.0_real64, min(40.0_real64, Ra))
        
    end function calculate_extraterrestrial_radiation

    subroutine fit_distribution_parameters(era5_precip, era5_pet, cmip6_data, bc_params, success)
        !===========================================================================
        ! Fit monthly distribution parameters for quantile mapping
        !===========================================================================
        real(real64), intent(in) :: era5_precip(:,:,:), era5_pet(:,:,:)
        type(cmip6_scenario_data_t), intent(in) :: cmip6_data
        type(bias_correction_parameters_t), intent(out) :: bc_params
        logical, intent(out) :: success
        
        integer :: i, j, month, overlap_months
        integer :: fitted_cells, total_cells
        
        success = .false.
        fitted_cells = 0
        total_cells = 0
        
        ! Initialize parameters structure
        bc_params%nlons = cmip6_data%nlons
        bc_params%nlats = cmip6_data%nlats
        
        allocate(bc_params%era5_precip_shape(bc_params%nlons, bc_params%nlats, 12))
        allocate(bc_params%era5_precip_scale(bc_params%nlons, bc_params%nlats, 12))
        allocate(bc_params%cmip6_precip_shape(bc_params%nlons, bc_params%nlats, 12))
        allocate(bc_params%cmip6_precip_scale(bc_params%nlons, bc_params%nlats, 12))
        
        allocate(bc_params%era5_pet_mean(bc_params%nlons, bc_params%nlats, 12))
        allocate(bc_params%era5_pet_std(bc_params%nlons, bc_params%nlats, 12))
        allocate(bc_params%cmip6_pet_mean(bc_params%nlons, bc_params%nlats, 12))
        allocate(bc_params%cmip6_pet_std(bc_params%nlons, bc_params%nlats, 12))
        
        ! Calculate overlap period (1981-2014)
        overlap_months = min((HIST_END_YEAR - HIST_START_YEAR + 1) * 12, &
                            size(era5_precip, 3), cmip6_data%ntimes)
        
        write(output_unit, '(A,I0,A)') '📊 Fitting distributions for ', overlap_months, ' months overlap'
        
        ! Fit parameters for each grid cell and month
        do j = 1, bc_params%nlats
            do i = 1, bc_params%nlons
                if (cmip6_data%land_mask(i,j)) then
                    total_cells = total_cells + 1
                    
                    call fit_cell_parameters(era5_precip(i,j,:), era5_pet(i,j,:), &
                                           cmip6_data%precipitation(i,j,:), cmip6_data%pet_raw(i,j,:), &
                                           overlap_months, bc_params, i, j, fitted_cells)
                end if
            end do
        end do
        
        ! Check if enough cells were fitted
        if (fitted_cells >= total_cells / 10) then  ! At least 10% success rate
            bc_params%fitted = .true.
            success = .true.
            write(output_unit, '(A,I0,A,I0,A)') '✅ Successfully fitted ', fitted_cells, '/', total_cells, ' land cells'
        else
            write(error_unit, '(A,I0,A,I0,A)') '❌ Insufficient fits: ', fitted_cells, '/', total_cells, ' cells'
        end if
        
    end subroutine fit_distribution_parameters

    subroutine fit_cell_parameters(era5_p, era5_pet, cmip6_p, cmip6_pet, overlap_months, &
                                  bc_params, i, j, fitted_count)
        !===========================================================================
        ! Fit distribution parameters for a single grid cell
        !===========================================================================
        real(real64), intent(in) :: era5_p(:), era5_pet(:), cmip6_p(:), cmip6_pet(:)
        integer, intent(in) :: overlap_months, i, j
        type(bias_correction_parameters_t), intent(inout) :: bc_params
        integer, intent(inout) :: fitted_count
        
        integer :: month, year_idx, month_idx, valid_years
        real(real64), allocatable :: era5_month_p(:), cmip6_month_p(:)
        real(real64), allocatable :: era5_month_pet(:), cmip6_month_pet(:)
        logical :: month_fitted
        
        month_fitted = .false.
        
        ! Fit parameters for each month
        do month = 1, 12
            ! Extract monthly values for the overlap period
            valid_years = overlap_months / 12
            allocate(era5_month_p(valid_years), cmip6_month_p(valid_years))
            allocate(era5_month_pet(valid_years), cmip6_month_pet(valid_years))
            
            call extract_monthly_values(era5_p, cmip6_p, era5_pet, cmip6_pet, &
                                      month, overlap_months, &
                                      era5_month_p, cmip6_month_p, &
                                      era5_month_pet, cmip6_month_pet, valid_years)
            
            if (valid_years >= 5) then  ! Minimum years for reliable fitting
                ! Fit precipitation (gamma distribution)
                call fit_gamma_distribution(era5_month_p(1:valid_years), &
                                          bc_params%era5_precip_shape(i,j,month), &
                                          bc_params%era5_precip_scale(i,j,month))
                
                call fit_gamma_distribution(cmip6_month_p(1:valid_years), &
                                          bc_params%cmip6_precip_shape(i,j,month), &
                                          bc_params%cmip6_precip_scale(i,j,month))
                
                ! Fit PET (normal distribution)
                call fit_normal_distribution(era5_month_pet(1:valid_years), &
                                           bc_params%era5_pet_mean(i,j,month), &
                                           bc_params%era5_pet_std(i,j,month))
                
                call fit_normal_distribution(cmip6_month_pet(1:valid_years), &
                                           bc_params%cmip6_pet_mean(i,j,month), &
                                           bc_params%cmip6_pet_std(i,j,month))
                
                month_fitted = .true.
            else
                ! Set default parameters if insufficient data
                bc_params%era5_precip_shape(i,j,month) = 2.0_real64
                bc_params%era5_precip_scale(i,j,month) = 1.0_real64
                bc_params%cmip6_precip_shape(i,j,month) = 2.0_real64
                bc_params%cmip6_precip_scale(i,j,month) = 1.0_real64
                
                bc_params%era5_pet_mean(i,j,month) = 100.0_real64
                bc_params%era5_pet_std(i,j,month) = 50.0_real64
                bc_params%cmip6_pet_mean(i,j,month) = 100.0_real64
                bc_params%cmip6_pet_std(i,j,month) = 50.0_real64
            end if
            
            deallocate(era5_month_p, cmip6_month_p, era5_month_pet, cmip6_month_pet)
        end do
        
        if (month_fitted) fitted_count = fitted_count + 1
        
    end subroutine fit_cell_parameters

    subroutine extract_monthly_values(era5_p, cmip6_p, era5_pet, cmip6_pet, target_month, &
                                     overlap_months, era5_month_p, cmip6_month_p, &
                                     era5_month_pet, cmip6_month_pet, valid_count)
        !===========================================================================
        ! Extract values for a specific month from time series
        !===========================================================================
        real(real64), intent(in) :: era5_p(:), cmip6_p(:), era5_pet(:), cmip6_pet(:)
        integer, intent(in) :: target_month, overlap_months
        real(real64), intent(out) :: era5_month_p(:), cmip6_month_p(:)
        real(real64), intent(out) :: era5_month_pet(:), cmip6_month_pet(:)
        integer, intent(inout) :: valid_count
        
        integer :: year_start, month_idx, year_count
        
        valid_count = 0
        
        do year_start = 1, overlap_months, 12
            month_idx = year_start + target_month - 1
            if (month_idx <= overlap_months) then
                if (era5_p(month_idx) > MIN_PRECIP .and. cmip6_p(month_idx) > MIN_PRECIP .and. &
                    era5_pet(month_idx) > 0.0_real64 .and. cmip6_pet(month_idx) > 0.0_real64) then
                    
                    valid_count = valid_count + 1
                    if (valid_count <= size(era5_month_p)) then
                        era5_month_p(valid_count) = era5_p(month_idx)
                        cmip6_month_p(valid_count) = cmip6_p(month_idx)
                        era5_month_pet(valid_count) = era5_pet(month_idx)
                        cmip6_month_pet(valid_count) = cmip6_pet(month_idx)
                    end if
                end if
            end if
        end do
        
    end subroutine extract_monthly_values

    subroutine fit_gamma_distribution(data, shape, scale)
        !===========================================================================
        ! Fit gamma distribution using method of moments
        !===========================================================================
        real(real64), intent(in) :: data(:)
        real(real64), intent(out) :: shape, scale
        
        real(real64) :: data_mean, data_var
        integer :: n
        
        n = size(data)
        
        ! Calculate sample statistics
        data_mean = sum(data) / real(n, real64)
        data_var = sum((data - data_mean)**2) / real(max(1, n-1), real64)
        
        ! Method of moments estimators for gamma distribution
        if (data_var > 1e-10_real64 .and. data_mean > 1e-10_real64) then
            scale = data_var / data_mean           ! Scale parameter
            shape = data_mean / scale              ! Shape parameter
            
            ! Ensure reasonable parameter bounds
            shape = max(0.1_real64, min(50.0_real64, shape))
            scale = max(0.01_real64, min(1000.0_real64, scale))
        else
            ! Default parameters for degenerate cases
            shape = 2.0_real64
            scale = 1.0_real64
        end if
        
    end subroutine fit_gamma_distribution

    subroutine fit_normal_distribution(data, mean, std)
        !===========================================================================
        ! Fit normal distribution parameters
        !===========================================================================
        real(real64), intent(in) :: data(:)
        real(real64), intent(out) :: mean, std
        
        integer :: n
        
        n = size(data)
        
        ! Calculate sample statistics
        mean = sum(data) / real(n, real64)
        std = sqrt(sum((data - mean)**2) / real(max(1, n-1), real64))
        
        ! Ensure positive standard deviation
        std = max(0.1_real64, std)
        
    end subroutine fit_normal_distribution

    subroutine apply_quantile_mapping_correction(cmip6_data, bc_params, success)
        !===========================================================================
        ! Apply quantile mapping using fitted distribution parameters
        !===========================================================================
        type(cmip6_scenario_data_t), intent(inout) :: cmip6_data
        type(bias_correction_parameters_t), intent(in) :: bc_params
        logical, intent(out) :: success
        
        integer :: i, j, k, month
        real(real64) :: quantile, corrected_value
        
        success = .false.
        
        if (.not. bc_params%fitted) then
            write(error_unit, '(A)') '❌ Cannot apply correction: parameters not fitted'
            return
        end if
        
        ! Allocate corrected arrays if not already allocated
        if (.not. allocated(cmip6_data%precipitation_corrected)) then
            allocate(cmip6_data%precipitation_corrected(cmip6_data%nlons, cmip6_data%nlats, cmip6_data%ntimes))
        end if
        
        if (.not. allocated(cmip6_data%pet_corrected)) then
            allocate(cmip6_data%pet_corrected(cmip6_data%nlons, cmip6_data%nlats, cmip6_data%ntimes))
        end if
        
        ! Apply quantile mapping to each time step
        do k = 1, cmip6_data%ntimes
            month = mod(k-1, 12) + 1
            
            do j = 1, cmip6_data%nlats
                do i = 1, cmip6_data%nlons
                    if (cmip6_data%land_mask(i,j)) then
                        
                        ! Correct precipitation (gamma distribution)
                        if (cmip6_data%precipitation(i,j,k) > MIN_PRECIP) then
                            quantile = gamma_cdf(cmip6_data%precipitation(i,j,k), &
                                               bc_params%cmip6_precip_shape(i,j,month), &
                                               bc_params%cmip6_precip_scale(i,j,month))
                            
                            corrected_value = gamma_ppf(quantile, &
                                                       bc_params%era5_precip_shape(i,j,month), &
                                                       bc_params%era5_precip_scale(i,j,month))
                            
                            cmip6_data%precipitation_corrected(i,j,k) = max(0.0_real64, corrected_value)
                        else
                            cmip6_data%precipitation_corrected(i,j,k) = 0.0_real64
                        end if
                        
                        ! Correct PET (normal distribution)
                        if (cmip6_data%pet_raw(i,j,k) > 0.0_real64) then
                            quantile = normal_cdf(cmip6_data%pet_raw(i,j,k), &
                                                bc_params%cmip6_pet_mean(i,j,month), &
                                                bc_params%cmip6_pet_std(i,j,month))
                            
                            corrected_value = normal_ppf(quantile, &
                                                        bc_params%era5_pet_mean(i,j,month), &
                                                        bc_params%era5_pet_std(i,j,month))
                            
                            cmip6_data%pet_corrected(i,j,k) = max(0.1_real64, corrected_value)
                        else
                            cmip6_data%pet_corrected(i,j,k) = FILL_VALUE
                        end if
                        
                    else
                        cmip6_data%precipitation_corrected(i,j,k) = FILL_VALUE
                        cmip6_data%pet_corrected(i,j,k) = FILL_VALUE
                    end if
                end do
            end do
        end do
        
        success = .true.
        
    end subroutine apply_quantile_mapping_correction

    function gamma_cdf(x, shape, scale) result(cdf)
        !===========================================================================
        ! Gamma cumulative distribution function using FSML (proper implementation)
        !===========================================================================
        real(real64), intent(in) :: x, shape, scale
        real(real64) :: cdf
        
        if (x <= 0.0_real64) then
            cdf = 0.0_real64
        else
            ! Use FSML gamma CDF with proper parameterization
            cdf = f_dst_gamma_cdf_core(x, alpha=shape, beta=scale, loc=0.0_real64, tail="left")
        end if
        
    end function gamma_cdf

    function gamma_ppf(p, shape, scale) result(x)
        !===========================================================================
        ! Gamma percent point function using FSML (proper implementation)
        !===========================================================================
        real(real64), intent(in) :: p, shape, scale
        real(real64) :: x
        
        if (p <= 0.0_real64) then
            x = 0.0_real64
        else if (p >= 1.0_real64) then
            x = shape * scale * 10.0_real64  ! Large value for extreme quantiles
        else
            ! Use FSML gamma PPF with proper parameterization
            x = f_dst_gamma_ppf_core(p, alpha=shape, beta=scale, loc=0.0_real64)
            x = max(0.0_real64, x)  ! Ensure non-negative
        end if
        
    end function gamma_ppf

    function normal_cdf(x, mean, std) result(cdf)
        !===========================================================================
        ! Normal cumulative distribution function using FSML
        !===========================================================================
        real(real64), intent(in) :: x, mean, std
        real(real64) :: cdf
        
        cdf = f_dst_norm_cdf_core(x, mu=mean, sigma=std, tail="left")
        
    end function normal_cdf

    function normal_ppf(p, mean, std) result(x)
        !===========================================================================
        ! Normal percent point function using FSML
        !===========================================================================
        real(real64), intent(in) :: p, mean, std
        real(real64) :: x
        
        if (p <= 0.0_real64) then
            x = mean - 6.0_real64 * std
        else if (p >= 1.0_real64) then
            x = mean + 6.0_real64 * std
        else
            x = f_dst_norm_ppf_core(p, mu=mean, sigma=std)
        end if
        
    end function normal_ppf

    subroutine package_corrected_scenarios(ssp126_data, ssp245_data, ssp585_data, corrected_scenarios, success)
        !===========================================================================
        ! Package corrected data into structured format for projection module
        !===========================================================================
        type(cmip6_scenario_data_t), intent(in) :: ssp126_data, ssp245_data, ssp585_data
        type(corrected_cmip6_scenarios_t), intent(out) :: corrected_scenarios
        logical, intent(out) :: success
        
        success = .false.
        
        ! Allocate coordinate arrays
        allocate(corrected_scenarios%longitudes(ssp126_data%nlons))
        allocate(corrected_scenarios%latitudes(ssp126_data%nlats))
        corrected_scenarios%longitudes = ssp126_data%longitude
        corrected_scenarios%latitudes = ssp126_data%latitude
        
        ! Allocate scenario arrays
        allocate(corrected_scenarios%ssp126_precip(ssp126_data%nlons, ssp126_data%nlats, ssp126_data%ntimes))
        allocate(corrected_scenarios%ssp126_pet(ssp126_data%nlons, ssp126_data%nlats, ssp126_data%ntimes))
        
        allocate(corrected_scenarios%ssp245_precip(ssp245_data%nlons, ssp245_data%nlats, ssp245_data%ntimes))
        allocate(corrected_scenarios%ssp245_pet(ssp245_data%nlons, ssp245_data%nlats, ssp245_data%ntimes))
        
        allocate(corrected_scenarios%ssp585_precip(ssp585_data%nlons, ssp585_data%nlats, ssp585_data%ntimes))
        allocate(corrected_scenarios%ssp585_pet(ssp585_data%nlons, ssp585_data%nlats, ssp585_data%ntimes))
        
        ! Copy corrected data
        corrected_scenarios%ssp126_precip = ssp126_data%precipitation_corrected
        corrected_scenarios%ssp126_pet = ssp126_data%pet_corrected
        
        corrected_scenarios%ssp245_precip = ssp245_data%precipitation_corrected
        corrected_scenarios%ssp245_pet = ssp245_data%pet_corrected
        
        corrected_scenarios%ssp585_precip = ssp585_data%precipitation_corrected
        corrected_scenarios%ssp585_pet = ssp585_data%pet_corrected
        
        ! Set metadata
        corrected_scenarios%start_year = 2015  ! Start of future projections
        corrected_scenarios%end_year = 2100    ! End of future projections
        corrected_scenarios%data_available = .true.
        
        success = .true.
        
    end subroutine package_corrected_scenarios
    
    !===============================================================================
    ! Generate professional validation tables for bias correction analysis
    !===============================================================================
    subroutine generate_validation_tables(obs_precip, obs_pet, &
                                         original_precip, original_pet, &
                                         corrected_precip, corrected_pet, &
                                         scenario_precip, scenario_pet)
        implicit none
        real(real64), intent(in) :: obs_precip(:), obs_pet(:)
        real(real64), intent(in) :: original_precip(:), original_pet(:)
        real(real64), intent(in) :: corrected_precip(:), corrected_pet(:)
        real(real64), intent(in) :: scenario_precip(:,:), scenario_pet(:,:)  ! [time, scenario]
        
        type(bias_correction_validation_t) :: precip_metrics, pet_metrics
        integer :: unit_num, ios, i, valid_count
        real :: scenario_means(4), scenario_stds(4), coverage_percent
        character(len=256) :: filename
        
        ! Calculate validation metrics
        precip_metrics = calculate_bias_metrics(obs_precip, original_precip, corrected_precip)
        pet_metrics = calculate_bias_metrics(obs_pet, original_pet, corrected_pet)
        
        ! Table 1: Historical Validation Results (matching existing format)
        filename = 'data/Final_results/Tables/bias_correction/bias_correction_validation_historical.csv'
        open(newunit=unit_num, file=trim(filename), status='replace', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Error opening file: ', trim(filename)
            return
        end if
        
        write(unit_num, '(A)') '# Bias Correction Validation - Historical Period (1981-2014)'
        write(unit_num, '(A)') '# Comparison: ERA5 vs CMIP6 Raw vs CMIP6 Corrected'
        write(unit_num, '(A)') '# Generated by Clean Bias Correction Module'
        write(unit_num, '(A)') 'Variable,Mean_Obs,Mean_Model_Before,Mean_Model_After,Bias_Diff_Before,Bias_Diff_After,Bias_Ratio_Before,Bias_Ratio_After,RMSE_Before,RMSE_After,Correlation_Before,Correlation_After'
        
        ! Write precipitation data with appropriate precision
        write(unit_num, '(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F6.3,A,F6.3,A,F8.3,A,F8.3,A,F6.3,A,F6.3)') &
            'Precipitation,', &
            precip_metrics%mean_obs, ',', &
            precip_metrics%mean_original, ',', &
            precip_metrics%mean_corrected, ',', &
            precip_metrics%bias_before, ',', precip_metrics%bias_after, ',', &
            precip_metrics%bias_ratio_before, ',', precip_metrics%bias_ratio_after, ',', &
            precip_metrics%rmse_before, ',', precip_metrics%rmse_after, ',', &
            precip_metrics%correlation_before, ',', precip_metrics%correlation_after
        
        write(unit_num, '(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F6.3,A,F6.3,A,F8.3,A,F8.3,A,F6.3,A,F6.3)') &
            'PET,', &
            pet_metrics%mean_obs, ',', &
            pet_metrics%mean_original, ',', &
            pet_metrics%mean_corrected, ',', &
            pet_metrics%bias_before, ',', pet_metrics%bias_after, ',', &
            pet_metrics%bias_ratio_before, ',', pet_metrics%bias_ratio_after, ',', &
            pet_metrics%rmse_before, ',', pet_metrics%rmse_after, ',', &
            pet_metrics%correlation_before, ',', pet_metrics%correlation_after
        
        close(unit_num)
        write(output_unit, '(A,A)') 'Generated historical validation table: ', trim(filename)
        
        ! Table 2: Future Scenario Diagnostics (matching existing format)
        filename = 'data/Final_results/Tables/bias_correction/bias_correction_scenarios_diagnostics.csv'
        open(newunit=unit_num, file=trim(filename), status='replace', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') 'Error opening file: ', trim(filename)
            return
        end if
        
        write(unit_num, '(A)') '# Future Scenarios Coverage and Diagnostics (2015-2099)'
        write(unit_num, '(A)') '# Quality control metrics for bias-corrected CMIP6 data'
        write(unit_num, '(A)') '# Generated by Clean Bias Correction Module'
        write(unit_num, '(A)') 'Scenario,Variable,Valid_Coverage_Percent,Mean_mm_month,Std_mm_month,Notes'
        
        ! Calculate scenario statistics using FSML
        do i = 1, min(3, size(scenario_precip, 2))  ! SSP126, SSP245, SSP585
            if (i <= size(scenario_precip, 2)) then
                scenario_means(i) = real(fsml_mean(scenario_precip(:, i)))
                scenario_stds(i) = real(fsml_std(scenario_precip(:, i)))
                
                ! Calculate coverage (assuming non-missing values)
                valid_count = count(ieee_is_finite(scenario_precip(:, i)))
                coverage_percent = 100.0 * real(valid_count) / real(size(scenario_precip, 1))
                
                ! Write precipitation data
                if (i == 1) then
                    write(unit_num, '(A,F0.2,A,F0.6,A,F0.6,A)') 'SSP126,Precipitation,', &
                        coverage_percent, ',', scenario_means(i), ',', scenario_stds(i), ',Bias-corrected'
                else if (i == 2) then
                    write(unit_num, '(A,F0.2,A,F0.6,A,F0.6,A)') 'SSP245,Precipitation,', &
                        coverage_percent, ',', scenario_means(i), ',', scenario_stds(i), ',Bias-corrected'
                else if (i == 3) then
                    write(unit_num, '(A,F0.2,A,F0.6,A,F0.6,A)') 'SSP585,Precipitation,', &
                        coverage_percent, ',', scenario_means(i), ',', scenario_stds(i), ',Bias-corrected'
                end if
            end if
        end do
        
        ! Calculate PET scenario statistics
        do i = 1, min(3, size(scenario_pet, 2))
            if (i <= size(scenario_pet, 2)) then
                scenario_means(i) = real(fsml_mean(scenario_pet(:, i)))
                scenario_stds(i) = real(fsml_std(scenario_pet(:, i)))
                
                valid_count = count(ieee_is_finite(scenario_pet(:, i)))
                coverage_percent = 100.0 * real(valid_count) / real(size(scenario_pet, 1))
                
                ! Write PET data
                if (i == 1) then
                    write(unit_num, '(A,F0.2,A,F0.6,A,F0.6,A)') 'SSP126,PET,', &
                        coverage_percent, ',', scenario_means(i), ',', scenario_stds(i), ',Bias-corrected'
                else if (i == 2) then
                    write(unit_num, '(A,F0.2,A,F0.6,A,F0.6,A)') 'SSP245,PET,', &
                        coverage_percent, ',', scenario_means(i), ',', scenario_stds(i), ',Bias-corrected'
                else if (i == 3) then
                    write(unit_num, '(A,F0.2,A,F0.6,A,F0.6,A)') 'SSP585,PET,', &
                        coverage_percent, ',', scenario_means(i), ',', scenario_stds(i), ',Bias-corrected'
                end if
            end if
        end do
        
        close(unit_num)
        write(output_unit, '(A,A)') 'Generated scenario diagnostics table: ', trim(filename)
        
    end subroutine generate_validation_tables
    
    !===============================================================================
    ! Wrapper for generating validation tables with real 3D climate data
    !===============================================================================
    subroutine generate_bias_correction_validation_tables(era5_precip, era5_pet, &
                                                         cmip6_historical, &
                                                         corrected_scenarios, success)
        implicit none
        real(real64), intent(in) :: era5_precip(:,:,:), era5_pet(:,:,:)
        type(cmip6_scenario_data_t), intent(in) :: cmip6_historical
        type(corrected_cmip6_scenarios_t), intent(in) :: corrected_scenarios
        logical, intent(out) :: success
        
        ! Local arrays for spatial-temporal averages
        integer :: total_points, valid_points, i, j, k
        integer :: n_times_era5, n_times_cmip6, min_times
        integer :: era5_nlons, era5_nlats, cmip6_nlons, cmip6_nlats
        real(real64), allocatable :: obs_precip_1d(:), obs_pet_1d(:)
        real(real64), allocatable :: original_precip_1d(:), original_pet_1d(:)
        real(real64), allocatable :: corrected_precip_1d(:), corrected_pet_1d(:)
        real(real64), allocatable :: scenario_precip_2d(:,:), scenario_pet_2d(:,:)
        ! Spatial averaging variables
        real(real64) :: sum_era5_precip, sum_era5_pet, sum_cmip6_precip, sum_cmip6_pet
        real(real64) :: sum_corr_precip, sum_corr_pet, weight_sum, weight_sum_scenarios_precip, weight_sum_scenarios_pet, weight, lat_rad
        ! Validation filtering variables
        real(real64), allocatable :: valid_obs_precip(:), valid_obs_pet(:)
        real(real64), allocatable :: valid_orig_precip(:), valid_orig_pet(:)
        real(real64), allocatable :: valid_corr_precip(:), valid_corr_pet(:)
        integer :: valid_precip_count, valid_pet_count, idx
        
        success = .false.
        
        ! Get and verify dimensions
        era5_nlons = size(era5_precip, 1)
        era5_nlats = size(era5_precip, 2)
        n_times_era5 = size(era5_precip, 3)
        
        cmip6_nlons = cmip6_historical%nlons
        cmip6_nlats = cmip6_historical%nlats
        n_times_cmip6 = cmip6_historical%ntimes
        
        write(output_unit, '(A,I0,A,I0)') '   ERA5 grid dimensions: ', era5_nlons, ' x ', era5_nlats
        write(output_unit, '(A,I0,A,I0)') '   CMIP6 grid dimensions: ', cmip6_nlons, ' x ', cmip6_nlats
        write(output_unit, '(A,I0,A,I0)') '   Time dimensions - ERA5: ', n_times_era5, ', CMIP6: ', n_times_cmip6
        
        ! Check latitude ordering (should be north to south for proper indexing)
        if (cmip6_historical%latitude(1) < cmip6_historical%latitude(cmip6_nlats)) then
            write(output_unit, '(A)') '   ✓ Latitude ordering: South to North (standard)'
        else
            write(output_unit, '(A)') '   ✓ Latitude ordering: North to South (standard)'
        end if
        write(output_unit, '(A,F8.2,A,F8.2)') '   Latitude range: ', &
            cmip6_historical%latitude(1), ' to ', cmip6_historical%latitude(cmip6_nlats)
        
        ! Verify grids match (both should be regridded to same grid)
        if (era5_nlons /= cmip6_nlons .or. era5_nlats /= cmip6_nlats) then
            write(error_unit, '(A)') 'ERROR: ERA5 and CMIP6 grid dimensions do not match!'
            write(error_unit, '(A,I0,A,I0)') '  ERA5: ', era5_nlons, ' x ', era5_nlats
            write(error_unit, '(A,I0,A,I0)') '  CMIP6: ', cmip6_nlons, ' x ', cmip6_nlats
            return
        end if
        
        min_times = min(n_times_era5, n_times_cmip6)
        total_points = cmip6_nlons * cmip6_nlats
        valid_points = count(cmip6_historical%land_mask)
        
        write(output_unit, '(A,I0,A,I0,A)') '   Processing ', valid_points, ' land points out of ', total_points, ' total points'
        write(output_unit, '(A,I0,A)') '   Using ', min_times, ' overlapping time steps'
        
        ! Allocate 1D arrays for spatial-temporal averages
        allocate(obs_precip_1d(min_times))
        allocate(obs_pet_1d(min_times))
        allocate(original_precip_1d(min_times))
        allocate(original_pet_1d(min_times))
        allocate(corrected_precip_1d(min_times))
        allocate(corrected_pet_1d(min_times))
        allocate(scenario_precip_2d(min_times, 3))  ! SSP126, SSP245, SSP585
        allocate(scenario_pet_2d(min_times, 3))
        
        ! Calculate spatial averages for each time step (land points only, latitude-weighted)
        do k = 1, min_times
            ! Initialize sums and weights
            sum_era5_precip = 0.0_real64
            sum_era5_pet = 0.0_real64
            sum_cmip6_precip = 0.0_real64
            sum_cmip6_pet = 0.0_real64
            sum_corr_precip = 0.0_real64
            sum_corr_pet = 0.0_real64
            weight_sum = 0.0_real64
            
            do j = 1, cmip6_nlats
                ! Calculate latitude weight (cos weighting for proper spatial averaging)
                lat_rad = cmip6_historical%latitude(j) * PI / 180.0_real64
                weight = cos(lat_rad)
                
                do i = 1, cmip6_nlons
                    if (cmip6_historical%land_mask(i,j)) then
                        ! Check ERA5 data validity (finite and non-fill)
                        if (ieee_is_finite(era5_precip(i,j,k)) .and. era5_precip(i,j,k) >= 0.0_real64 .and. &
                            era5_precip(i,j,k) /= FILL_VALUE .and. &
                            ieee_is_finite(era5_pet(i,j,k)) .and. era5_pet(i,j,k) > 0.0_real64 .and. &
                            era5_pet(i,j,k) /= FILL_VALUE) then
                            
                            ! Check CMIP6 original data validity
                            if (ieee_is_finite(cmip6_historical%precipitation(i,j,k)) .and. &
                                cmip6_historical%precipitation(i,j,k) >= 0.0_real64 .and. &
                                cmip6_historical%precipitation(i,j,k) /= FILL_VALUE .and. &
                                ieee_is_finite(cmip6_historical%pet_raw(i,j,k)) .and. &
                                cmip6_historical%pet_raw(i,j,k) > 0.0_real64 .and. &
                                cmip6_historical%pet_raw(i,j,k) /= FILL_VALUE) then
                                
                                ! Check CMIP6 corrected data validity
                                if (ieee_is_finite(cmip6_historical%precipitation_corrected(i,j,k)) .and. &
                                    cmip6_historical%precipitation_corrected(i,j,k) >= 0.0_real64 .and. &
                                    cmip6_historical%precipitation_corrected(i,j,k) /= FILL_VALUE .and. &
                                    ieee_is_finite(cmip6_historical%pet_corrected(i,j,k)) .and. &
                                    cmip6_historical%pet_corrected(i,j,k) > 0.0_real64 .and. &
                                    cmip6_historical%pet_corrected(i,j,k) /= FILL_VALUE) then
                                    
                                    ! Accumulate weighted sums
                                    sum_era5_precip = sum_era5_precip + weight * era5_precip(i,j,k)
                                    sum_era5_pet = sum_era5_pet + weight * era5_pet(i,j,k)
                                    
                                    sum_cmip6_precip = sum_cmip6_precip + weight * cmip6_historical%precipitation(i,j,k)
                                    sum_cmip6_pet = sum_cmip6_pet + weight * cmip6_historical%pet_raw(i,j,k)
                                    
                                    sum_corr_precip = sum_corr_precip + weight * cmip6_historical%precipitation_corrected(i,j,k)
                                    sum_corr_pet = sum_corr_pet + weight * cmip6_historical%pet_corrected(i,j,k)
                                    
                                    weight_sum = weight_sum + weight
                                end if
                            end if
                        end if
                    end if
                end do
            end do
            
            ! Calculate weighted spatial averages
            if (weight_sum > 0.0_real64) then
                obs_precip_1d(k) = sum_era5_precip / weight_sum
                obs_pet_1d(k) = sum_era5_pet / weight_sum
                original_precip_1d(k) = sum_cmip6_precip / weight_sum
                original_pet_1d(k) = sum_cmip6_pet / weight_sum
                corrected_precip_1d(k) = sum_corr_precip / weight_sum
                corrected_pet_1d(k) = sum_corr_pet / weight_sum
            else
                ! No valid data for this time step
                obs_precip_1d(k) = FILL_VALUE
                obs_pet_1d(k) = FILL_VALUE
                original_precip_1d(k) = FILL_VALUE
                original_pet_1d(k) = FILL_VALUE
                corrected_precip_1d(k) = FILL_VALUE
                corrected_pet_1d(k) = FILL_VALUE
            end if
        end do
        
        ! Extract scenario data (using same time period for consistency)
        do k = 1, min_times
            ! Calculate spatial averages manually for scenarios (with latitude weighting)
            valid_points = 0
            weight_sum_scenarios_precip = 0.0_real64
            weight_sum_scenarios_pet = 0.0_real64
            scenario_precip_2d(k, 1) = 0.0_real64
            scenario_precip_2d(k, 2) = 0.0_real64
            scenario_precip_2d(k, 3) = 0.0_real64
            scenario_pet_2d(k, 1) = 0.0_real64
            scenario_pet_2d(k, 2) = 0.0_real64
            scenario_pet_2d(k, 3) = 0.0_real64
            
            ! Check if we have scenario data available
            if (k <= size(corrected_scenarios%ssp126_precip, 3)) then
                do j = 1, size(corrected_scenarios%ssp126_precip, 2)
                    ! Calculate latitude weight (same as historical)
                    lat_rad = cmip6_historical%latitude(j) * PI / 180.0_real64
                    weight = cos(lat_rad)
                    
                    do i = 1, size(corrected_scenarios%ssp126_precip, 1)
                        ! Check precipitation validity and add to precip sums
                        if (ieee_is_finite(corrected_scenarios%ssp126_precip(i,j,k)) .and. &
                            corrected_scenarios%ssp126_precip(i,j,k) >= 0.0_real64) then
                            
                            scenario_precip_2d(k, 1) = scenario_precip_2d(k, 1) + weight * corrected_scenarios%ssp126_precip(i,j,k)
                            scenario_precip_2d(k, 2) = scenario_precip_2d(k, 2) + weight * corrected_scenarios%ssp245_precip(i,j,k)
                            scenario_precip_2d(k, 3) = scenario_precip_2d(k, 3) + weight * corrected_scenarios%ssp585_precip(i,j,k)
                            
                            weight_sum_scenarios_precip = weight_sum_scenarios_precip + weight
                            valid_points = valid_points + 1
                        end if
                        
                        ! Check PET validity separately and add to PET sums
                        if (ieee_is_finite(corrected_scenarios%ssp126_pet(i,j,k)) .and. &
                            corrected_scenarios%ssp126_pet(i,j,k) > 0.0_real64) then
                            
                            scenario_pet_2d(k, 1) = scenario_pet_2d(k, 1) + weight * corrected_scenarios%ssp126_pet(i,j,k)
                            scenario_pet_2d(k, 2) = scenario_pet_2d(k, 2) + weight * corrected_scenarios%ssp245_pet(i,j,k)
                            scenario_pet_2d(k, 3) = scenario_pet_2d(k, 3) + weight * corrected_scenarios%ssp585_pet(i,j,k)
                            
                            weight_sum_scenarios_pet = weight_sum_scenarios_pet + weight
                        end if
                    end do
                end do
                
                ! Calculate weighted averages
                if (weight_sum_scenarios_precip > 0.0_real64) then
                    scenario_precip_2d(k, 1) = scenario_precip_2d(k, 1) / weight_sum_scenarios_precip
                    scenario_precip_2d(k, 2) = scenario_precip_2d(k, 2) / weight_sum_scenarios_precip
                    scenario_precip_2d(k, 3) = scenario_precip_2d(k, 3) / weight_sum_scenarios_precip
                end if
                
                if (weight_sum_scenarios_pet > 0.0_real64) then
                    scenario_pet_2d(k, 1) = scenario_pet_2d(k, 1) / weight_sum_scenarios_pet
                    scenario_pet_2d(k, 2) = scenario_pet_2d(k, 2) / weight_sum_scenarios_pet
                    scenario_pet_2d(k, 3) = scenario_pet_2d(k, 3) / weight_sum_scenarios_pet
                end if
            end if
        end do
        
        ! Filter out FILL_VALUE entries and calculate statistics on valid data only
        ! Count valid time steps SEPARATELY for precipitation and PET
        valid_precip_count = count(obs_precip_1d /= FILL_VALUE .and. &
                                  original_precip_1d /= FILL_VALUE .and. &
                                  corrected_precip_1d /= FILL_VALUE)
        
        valid_pet_count = count(obs_pet_1d /= FILL_VALUE .and. &
                               original_pet_1d /= FILL_VALUE .and. &
                               corrected_pet_1d /= FILL_VALUE)
        
        write(output_unit, '(A,I0,A,I0,A)') '   Valid precip time steps: ', valid_precip_count, ' out of ', min_times
        write(output_unit, '(A,I0,A,I0,A)') '   Valid PET time steps: ', valid_pet_count, ' out of ', min_times
        
        if (valid_precip_count < 10) then
            write(error_unit, '(A)') 'ERROR: Too few valid precipitation time steps for reliable validation'
            return
        end if
        
        if (valid_pet_count < 10) then
            write(error_unit, '(A)') 'ERROR: Too few valid PET time steps for reliable validation'
            return
        end if
        
        ! Extract valid data points SEPARATELY for each variable
        allocate(valid_obs_precip(valid_precip_count), valid_orig_precip(valid_precip_count), valid_corr_precip(valid_precip_count))
        allocate(valid_obs_pet(valid_pet_count), valid_orig_pet(valid_pet_count), valid_corr_pet(valid_pet_count))
        
        ! Fill precipitation arrays
        idx = 0
        do k = 1, min_times
            if (obs_precip_1d(k) /= FILL_VALUE .and. &
                original_precip_1d(k) /= FILL_VALUE .and. &
                corrected_precip_1d(k) /= FILL_VALUE) then
                idx = idx + 1
                valid_obs_precip(idx) = obs_precip_1d(k)
                valid_orig_precip(idx) = original_precip_1d(k)
                valid_corr_precip(idx) = corrected_precip_1d(k)
            end if
        end do
        
        ! Fill PET arrays
        idx = 0
        do k = 1, min_times
            if (obs_pet_1d(k) /= FILL_VALUE .and. &
                original_pet_1d(k) /= FILL_VALUE .and. &
                corrected_pet_1d(k) /= FILL_VALUE) then
                idx = idx + 1
                valid_obs_pet(idx) = obs_pet_1d(k)
                valid_orig_pet(idx) = original_pet_1d(k)
                valid_corr_pet(idx) = corrected_pet_1d(k)
            end if
        end do
        
        write(output_unit, '(A,F8.2,A)') '   ERA5 precipitation mean: ', real(fsml_mean(valid_obs_precip)), ' mm/month'
        write(output_unit, '(A,F8.2,A)') '   CMIP6 original precipitation mean: ', real(fsml_mean(valid_orig_precip)), ' mm/month'
        write(output_unit, '(A,F8.2,A)') '   CMIP6 corrected precipitation mean: ', real(fsml_mean(valid_corr_precip)), ' mm/month'
        write(output_unit, '(A,F8.2,A)') '   ERA5 PET mean: ', real(fsml_mean(valid_obs_pet)), ' mm/month'
        
        ! Debug scenario values 
        
        ! Check for any remaining negative values (should not happen now)
        if (minval(valid_obs_precip) < 0.0_real64) then
            write(error_unit, '(A,F8.2)') 'WARNING: Negative ERA5 precipitation values found, min: ', minval(valid_obs_precip)
        end if
        if (minval(valid_obs_pet) < 0.0_real64) then
            write(error_unit, '(A,F8.2)') 'WARNING: Negative ERA5 PET values found, min: ', minval(valid_obs_pet)
        end if
        
        ! Generate validation tables using the processed valid arrays
        call generate_validation_tables(valid_obs_precip, valid_obs_pet, &
                                       valid_orig_precip, valid_orig_pet, &
                                       valid_corr_precip, valid_corr_pet, &
                                       scenario_precip_2d, scenario_pet_2d)
        
        ! Clean up
        deallocate(obs_precip_1d, obs_pet_1d)
        deallocate(original_precip_1d, original_pet_1d)
        deallocate(corrected_precip_1d, corrected_pet_1d)
        deallocate(scenario_precip_2d, scenario_pet_2d)
        
        success = .true.
        
    end subroutine generate_bias_correction_validation_tables

    !===========================================================================
    ! GENERATE_SPEI12_COMPARISON_FILES: Create SPEI-12 NetCDF files for bias correction validation
    !===========================================================================
    subroutine generate_spei12_comparison_files(era5_precip, era5_pet, &
                                               corrected_scenarios, latitudes, longitudes, &
                                               era5_time_values, success, era5_spei12_output)
        !===========================================================================
        ! Purpose: Generate SPEI-12 comparison files for bias correction validation
        ! Uses the same SPEI calculation as the main pipeline for consistency
        ! 
        ! Outputs:
        !   1. era5_historical_spei12_1981_2014.nc (reference)
        !   2. cmip6_historical_spei12_corrected_1981_2014.nc (bias-corrected)
        !   3. era5_spei12_output: ERA5 SPEI-12 data for FSML comparison (optional)
        !===========================================================================
        use spei_module, only: calculate_single_spei
        use io_module, only: save_spei_single_scale, create_output_directory
        
        ! Input arguments
        real(real64), intent(in) :: era5_precip(:,:,:)     ! ERA5 precipitation [lon,lat,time]
        real(real64), intent(in) :: era5_pet(:,:,:)        ! ERA5 PET [lon,lat,time]
        type(corrected_cmip6_scenarios_t), intent(in) :: corrected_scenarios
        real(real64), intent(in) :: latitudes(:), longitudes(:)
        real(real64), intent(in) :: era5_time_values(:)
        logical, intent(out) :: success
        
        ! Optional output: ERA5 SPEI-12 for alignment with bias correction
        real(real64), allocatable, intent(out), optional :: era5_spei12_output(:,:,:)
        
        ! Local variables
        integer :: nlon, nlat, ntime_era5, ntime_overlap
        integer :: overlap_start_era5, overlap_start_cmip6
        integer :: i, j, status
        
        ! SPEI-12 calculation arrays (only SPEI-12, not all timescales)
        real(real64), allocatable :: spei12_era5(:,:,:), spei12_cmip6_corrected(:,:,:)
        real(real64), allocatable :: water_balance_era5(:,:,:), water_balance_cmip6(:,:,:)
        
        ! Time slicing arrays
        real(real64), allocatable :: era5_time_overlap(:)
        
        success = .false.
        
        write(output_unit, '(A)') ''
        write(output_unit, '(A)') '🔍 GENERATING SPEI-12 COMPARISON FILES FOR BIAS CORRECTION VALIDATION'
        write(output_unit, '(A)') '================================================================='
        
        ! Get dimensions
        nlon = size(era5_precip, 1)
        nlat = size(era5_precip, 2)
        ntime_era5 = size(era5_precip, 3)
        
        write(output_unit, '(A,I0,A,I0,A,I0,A)') '   Grid dimensions: ', nlon, ' × ', nlat, ' × ', ntime_era5, ' (ERA5)'
        
        ! Determine overlap period (1981-2014: 34 years = 408 months)
        ntime_overlap = 34 * 12  ! 1981-2014
        overlap_start_era5 = 1   ! ERA5 starts in 1981
        overlap_start_cmip6 = 1  ! Bias-corrected data should start from overlap period
        
        write(output_unit, '(A,I0,A)') '   Using overlap period: 1981-2014 (', ntime_overlap, ' months)'
        
        ! Validate time dimensions
        if (ntime_era5 < ntime_overlap) then
            write(error_unit, '(A,I0,A,I0)') 'ERROR: ERA5 time dimension too small: ', ntime_era5, ' < ', ntime_overlap
            return
        end if
        
        if (size(corrected_scenarios%ssp126_precip, 3) < ntime_overlap) then
            write(error_unit, '(A,I0,A,I0)') 'ERROR: Corrected CMIP6 time dimension too small: ', &
                size(corrected_scenarios%ssp126_precip, 3), ' < ', ntime_overlap
            return
        end if
        
        ! Create output directory
        call create_output_directory('data/Final_results/validation/')
        
        ! ================================================================
        ! STEP 1: CALCULATE SPEI-12 FOR ERA5 DATA (OVERLAP PERIOD)
        ! ================================================================
        write(output_unit, '(A)') '   📊 Calculating SPEI-12 for ERA5 data (1981-2014)...'
        
        ! Allocate arrays for SPEI-12 calculation only
        allocate(water_balance_era5(nlon, nlat, ntime_overlap))
        allocate(spei12_era5(nlon, nlat, ntime_overlap))
        
        ! Calculate water balance (P - PET) for overlap period
        do i = 1, nlon
            do j = 1, nlat
                water_balance_era5(i, j, :) = &
                    era5_precip(i, j, overlap_start_era5:overlap_start_era5+ntime_overlap-1) - &
                    era5_pet(i, j, overlap_start_era5:overlap_start_era5+ntime_overlap-1)
            end do
        end do
        
        ! Calculate SPEI-12 only using FSML
        call calculate_single_spei(water_balance_era5, 12, spei12_era5, status)
        if (status /= 0) then
            write(error_unit, '(A)') 'ERROR: SPEI-12 calculation failed for ERA5 data'
            return
        end if
        
        write(output_unit, '(A)') '      ✓ ERA5 SPEI-12 calculated successfully using FSML'
        
        ! ================================================================
        ! STEP 2: CALCULATE SPEI-12 FOR BIAS-CORRECTED CMIP6 DATA
        ! ================================================================
        write(output_unit, '(A)') '   📊 Calculating SPEI-12 for bias-corrected CMIP6 data (1981-2014)...'
        
        ! Allocate arrays for SPEI-12 calculation only
        allocate(water_balance_cmip6(nlon, nlat, ntime_overlap))
        allocate(spei12_cmip6_corrected(nlon, nlat, ntime_overlap))
        
        ! Calculate water balance (P - PET) for overlap period
        do i = 1, nlon
            do j = 1, nlat
                water_balance_cmip6(i, j, :) = &
                    corrected_scenarios%ssp126_precip(i, j, overlap_start_cmip6:overlap_start_cmip6+ntime_overlap-1) - &
                    corrected_scenarios%ssp126_pet(i, j, overlap_start_cmip6:overlap_start_cmip6+ntime_overlap-1)
            end do
        end do
        
        ! Calculate SPEI-12 only using FSML
        call calculate_single_spei(water_balance_cmip6, 12, spei12_cmip6_corrected, status)
        if (status /= 0) then
            write(error_unit, '(A)') 'ERROR: SPEI-12 calculation failed for bias-corrected CMIP6 data'
            return
        end if
        
        write(output_unit, '(A)') '      ✓ Bias-corrected CMIP6 SPEI-12 calculated successfully using FSML'
        
        ! ================================================================
        ! STEP 3: SAVE SPEI-12 NETCDF FILES
        ! ================================================================
        write(output_unit, '(A)') '   💾 Saving SPEI-12 comparison files...'
        
        ! Extract time values for overlap period
        allocate(era5_time_overlap(ntime_overlap))
        era5_time_overlap = era5_time_values(overlap_start_era5:overlap_start_era5+ntime_overlap-1)
        
        ! Save ERA5 SPEI-12
        call save_spei_single_scale(spei12_era5, 12, latitudes, longitudes, &
                                   'data/Final_results/validation/era5_historical_spei12_1981_2014.nc', status)
        if (status /= 0) then
            write(error_unit, '(A)') 'ERROR: Failed to save ERA5 SPEI-12 file'
            return
        end if
        
        write(output_unit, '(A)') '      ✓ ERA5 SPEI-12 saved: era5_historical_spei12_1981_2014.nc'
        
        ! Save bias-corrected CMIP6 SPEI-12
        call save_spei_single_scale(spei12_cmip6_corrected, 12, corrected_scenarios%latitudes, &
                                   corrected_scenarios%longitudes, &
                                   'data/Final_results/validation/cmip6_historical_spei12_corrected_1981_2014.nc', status)
        if (status /= 0) then
            write(error_unit, '(A)') 'ERROR: Failed to save bias-corrected CMIP6 SPEI-12 file'
            return
        end if
        
        write(output_unit, '(A)') '      ✓ Corrected CMIP6 SPEI-12 saved: cmip6_historical_spei12_corrected_1981_2014.nc'
        
        ! ================================================================
        ! STEP 4: SUMMARY
        ! ================================================================
        write(output_unit, '(A)') ''
        write(output_unit, '(A)') '   ✅ SPEI-12 COMPARISON FILES GENERATED SUCCESSFULLY!'
        write(output_unit, '(A)') '   📁 Location: data/Final_results/validation/'
        write(output_unit, '(A)') '   📊 Files created:'
        write(output_unit, '(A)') '      1. era5_historical_spei12_1981_2014.nc (reference)'
        write(output_unit, '(A)') '      2. cmip6_historical_spei12_corrected_1981_2014.nc (bias-corrected)'
        write(output_unit, '(A)') '   🎯 Use these files to validate bias correction effectiveness in drought indices'
        write(output_unit, '(A)') ''
        
        ! Optionally return ERA5 SPEI-12 data for FSML comparison alignment
        if (present(era5_spei12_output)) then
            allocate(era5_spei12_output(nlon, nlat, ntime_overlap))
            era5_spei12_output = spei12_era5
            write(output_unit, '(A)') '   📤 ERA5 SPEI-12 data (1981-2014) exported for FSML comparison'
        end if
        
        ! Cleanup
        if (allocated(water_balance_era5)) deallocate(water_balance_era5)
        if (allocated(water_balance_cmip6)) deallocate(water_balance_cmip6)
        if (allocated(spei12_era5)) deallocate(spei12_era5)
        if (allocated(spei12_cmip6_corrected)) deallocate(spei12_cmip6_corrected)
        if (allocated(era5_time_overlap)) deallocate(era5_time_overlap)
        
        success = .true.
        
    end subroutine generate_spei12_comparison_files

end module bias_correction_clean
