!> Climate Data I/O Module for Somaliland Drought Analysis
!! @author Khadar Daahir (University of Glasgow)
!! @date 2025-08-25

module io_module
    use, intrinsic :: iso_fortran_env, only: real32, real64, int64
    use netcdf
    implicit none
    
    private
    public :: load_era5_climate_data, save_processed_precipitation_data, save_processed_pet_data, &
              create_land_mask, apply_land_mask, save_spi_multiscale, save_spi_single_scale, &
              save_spei_multiscale, save_spei_single_scale, save_corrected_cmip6_precipitation, &
              save_corrected_cmip6_pet_hargreaves, &
              save_corrected_cmip6_scenario, load_corrected_cmip6_scenario, load_raw_cmip6_scenario, &
              save_future_drought_indices, save_evt_results_nc, save_evt_results_enhanced, &
              create_output_directory, get_evt_filename, get_processed_filename, save_evt_results, &
              load_netcdf_spei_data, corrected_cmip6_scenarios_t, add_bias_correction_metadata, &
              bias_correction_validation_t, save_corrected_cmip6_precipitation_enhanced, &
              save_corrected_cmip6_pet_enhanced, build_regridded_cmip6_file_paths, &
              read_regridded_cmip6_netcdf_data, load_regridded_cmip6_scenario, FILL_VALUE
    
    integer, parameter :: dp = real64
    real(dp), parameter :: FILL_VALUE = -999.0_dp
    character(len=*), parameter :: SOFTWARE_VERSION = "v1.0.0"
    character(len=*), parameter :: PROCESSING_HISTORY = "Created by Somaliland Drought Analysis Pipeline"
    real(dp), parameter :: LAT_MIN = 8.0_dp, LAT_MAX = 11.5_dp
    real(dp), parameter :: LON_MIN = 43.0_dp, LON_MAX = 48.5_dp
    real(dp), parameter :: M_TO_MM = 1000.0_dp
    real(dp), parameter :: DAYS_PER_MONTH = 30.44_dp
    real(dp), parameter :: LAND_THRESHOLD = 0.75_dp
    
    !> CMIP6 scenario data structure
    type, public :: cmip6_scenario_data_t
        ! Raw CMIP6 drivers
        real(dp), allocatable :: precipitation(:,:,:)  ! (lon, lat, time) [mm/month]
        real(dp), allocatable :: tasmax(:,:,:)         ! (lon, lat, time) [K]
        real(dp), allocatable :: tasmin(:,:,:)         ! (lon, lat, time) [K]
        
        ! PET arrays (raw and corrected)
        real(dp), allocatable :: pet_raw(:,:,:)               ! [mm/month] - CMIP6-calculated PET
        real(dp), allocatable :: pet_corrected(:,:,:)         ! [mm/month] - bias-corrected PET
        real(dp), allocatable :: precipitation_corrected(:,:,:) ! [mm/month] - bias-corrected P
        
        ! Coordinate arrays
        real(dp), allocatable :: longitude(:)          ! longitude values [degrees]
        real(dp), allocatable :: latitude(:)           ! latitude values [degrees]
        real(dp), allocatable :: time_values(:)        ! time values [days since ref]
        
        ! Data quality arrays
        logical, allocatable :: land_mask(:,:)             ! land/ocean mask
        real(dp), allocatable :: data_quality(:,:)     ! fraction of valid data per pixel
        
        ! Metadata
        integer :: nlons, nlats, ntimes                     ! dimensions
        character(len=50) :: scenario                       ! scenario name
        integer :: start_year, end_year                     ! temporal coverage
        logical :: is_corrected                             ! bias correction status
    end type cmip6_scenario_data_t
    
    !> Bias correction validation metrics structure
    type, public :: bias_correction_validation_t
        real(dp) :: bias            ! Mean bias
        real(dp) :: rmse            ! Root mean square error
        real(dp) :: correlation     ! Pearson correlation coefficient
        real(dp) :: ks_statistic    ! Kolmogorov-Smirnov statistic
        real(dp) :: mae             ! Mean absolute error
        logical :: validation_passed ! Overall validation status
    end type bias_correction_validation_t
    
    !> Corrected CMIP6 scenarios data structure (memory-based approach)
    type, public :: corrected_cmip6_scenarios_t
        ! Bias-corrected precipitation and PET for historical and future scenarios [lon, lat, time]
        real(dp), allocatable :: historical_precip(:,:,:), historical_pet(:,:,:)
        real(dp), allocatable :: ssp126_precip(:,:,:), ssp126_pet(:,:,:)
        real(dp), allocatable :: ssp245_precip(:,:,:), ssp245_pet(:,:,:)
        real(dp), allocatable :: ssp585_precip(:,:,:), ssp585_pet(:,:,:)
        
        ! Coordinate arrays
        real(dp), allocatable :: latitudes(:), longitudes(:), time_values(:)
        
        ! Metadata
        integer :: start_year, end_year
        logical :: data_available
        logical :: historical_available
    end type corrected_cmip6_scenarios_t
    
    interface save_evt_results
        module procedure save_evt_results_nc
        module procedure save_evt_results_enhanced
    end interface save_evt_results

contains

!> Create output directory if it doesn't exist
subroutine create_output_directory(dir_path)
    character(len=*), intent(in) :: dir_path
    character(len=500) :: mkdir_cmd
    
    mkdir_cmd = "mkdir -p " // trim(dir_path)
    call system(trim(mkdir_cmd))
end subroutine create_output_directory

function get_evt_filename(base_dir, scenario, timescale) result(filepath)
    character(len=*), intent(in) :: base_dir, scenario
    integer, intent(in) :: timescale
    character(len=512) :: filepath
    character(len=10) :: tscale_str
    
    write(tscale_str, '(I0)') timescale
    filepath = trim(base_dir) // "/EVT/" // trim(scenario) // "_evt_spei" // &
               trim(tscale_str) // ".nc"
end function get_evt_filename

function get_processed_filename(base_dir, data_type, scenario, timescale) result(filepath)
    character(len=*), intent(in) :: base_dir, data_type, scenario
    integer, intent(in), optional :: timescale
    character(len=512) :: filepath
    character(len=10) :: tscale_str
    
    if (present(timescale)) then
        write(tscale_str, '(I0)') timescale
        filepath = trim(base_dir) // "/" // trim(data_type) // "/" // &
                   trim(scenario) // "_" // trim(data_type) // "_" // &
                   trim(tscale_str) // "month.nc"
    else
        filepath = trim(base_dir) // "/" // trim(data_type) // "/" // &
                   trim(scenario) // "_" // trim(data_type) // "_processed.nc"
    end if
end function get_processed_filename

!> Add enhanced metadata to NetCDF files
subroutine add_enhanced_metadata(ncid, data_type, scenario, period_start, period_end)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: data_type, scenario
    integer, intent(in), optional :: period_start, period_end
    character(len=25) :: timestamp
    character(len=100) :: spatial_coverage, temporal_coverage
    
    ! Get timestamp
    call get_timestamp(timestamp)
    
    ! Spatial coverage metadata
    write(spatial_coverage, '(A,F5.1,A,F5.1,A,F5.1,A,F5.1,A)') &
        "Longitude: ", LON_MIN, " to ", LON_MAX, &
        ", Latitude: ", LAT_MIN, " to ", LAT_MAX, " degrees"
    
    ! Temporal coverage metadata
    if (present(period_start) .and. present(period_end)) then
        write(temporal_coverage, '(I4,A,I4)') period_start, "-", period_end
    else
        temporal_coverage = "1981-2024 baseline + 2025-2100 projections"
    end if
    
    ! Add enhanced global attributes
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "title", &
        "Somaliland Drought Analysis - " // trim(data_type)))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "institution", &
        "University of Glasgow"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", &
        "Khadar Daahir"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "source", &
        "Drought Analysis Pipeline - " // trim(scenario)))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "history", &
        PROCESSING_HISTORY))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "software_version", &
        SOFTWARE_VERSION))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "creation_date", &
        timestamp))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "spatial_coverage", &
        spatial_coverage))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "temporal_coverage", &
        temporal_coverage))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", &
        "CF-1.8"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "references", &
        "Coles (2001), Vicente-Serrano et al. (2010), McKee et al. (1993)"))
end subroutine add_enhanced_metadata

!> Add comprehensive bias correction metadata to NetCDF files
subroutine add_bias_correction_metadata(ncid, scenario, variable_type, validation_metrics)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: scenario, variable_type
    type(bias_correction_validation_t), intent(in), optional :: validation_metrics
    
    character(len=25) :: timestamp
    character(len=200) :: history_string
    character(len=100) :: method_description
    
    ! Get timestamp
    call get_timestamp(timestamp)
    
    ! ✅ General provenance
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "title", &
        "CMIP6 bias-corrected dataset - " // trim(variable_type) // " - " // trim(scenario)))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "institution", &
        "University of Glasgow"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "source", &
        "MPI-ESM1-2-LR CMIP6, corrected against ERA5-Land"))
    
    ! Build comprehensive history string
    write(history_string, '(A,A,A)') &
        "Bias correction applied ", trim(timestamp), &
        " using FSML-based Quantile Mapping (Gamma/Normal distributions)"
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "history", history_string))
    
    ! ✅ Methodology details
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "bias_correction_method", &
        "Quantile Mapping with distribution-specific parameters"))
    
    if (trim(variable_type) == "precipitation") then
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "distribution_type", "Gamma"))
        method_description = "Gamma distribution quantile mapping with shape/scale parameters"
    else if (trim(variable_type) == "pet" .or. trim(variable_type) == "PET") then
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "distribution_type", "Normal"))
        method_description = "Normal distribution quantile mapping with mean/std parameters"
    end if
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "method_description", method_description))
    
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "regridding_method", &
        "Bilinear interpolation (Fortran implementation)"))
    
    ! ✅ Dataset references
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "reference_data", &
        "ERA5-Land (ECMWF), 0.1° resolution"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "model_data", &
        "CMIP6 MPI-ESM1-2-LR, 1.25° resolution"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "scenario", scenario))
    
    ! Set time coverage based on scenario
    if (trim(scenario) == "historical") then
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_start", "1981-01"))
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_end", "2014-12"))
    else
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_start", "2015-01"))
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_end", "2099-12"))
    end if
    
    ! ✅ Quality control
    if (present(validation_metrics)) then
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "validation_performed", "true"))
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "mean_bias", validation_metrics%bias))
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "rmse", validation_metrics%rmse))
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "correlation", validation_metrics%correlation))
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "ks_statistic", validation_metrics%ks_statistic))
    else
        call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "validation_performed", "false"))
    end if
    
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "missing_value", "-999.0"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "conventions", "CF-1.8"))
    
    ! Processing metadata
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "processing_software", &
        "Fortran drought pipeline with FSML statistical library"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "creation_date", timestamp))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", "Khadar Daahir"))
    call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "project", &
        "Somaliland Drought Analysis - University of Glasgow"))
    
end subroutine add_bias_correction_metadata

!> Load ERA5-Land climate data with unit conversion and land masking
!! Converts ERA5 tp (m/day) and pev (m/day) to mm/month for drought analysis
subroutine load_era5_climate_data(precip_file, pet_file, precipitation, pet, &
                                 latitudes, longitudes, time_stamps, status)
    character(len=*), intent(in) :: precip_file, pet_file
    real(dp), intent(out) :: precipitation(:,:,:)  ! (lon,lat,time) in mm/month
    real(dp), intent(out) :: pet(:,:,:)           ! (lon,lat,time) in mm/month
    real(dp), intent(out) :: latitudes(:)         ! degrees_north
    real(dp), intent(out) :: longitudes(:)        ! degrees_east  
    real(dp), intent(out) :: time_stamps(:)       ! seconds since 1970-01-01
    integer, intent(out) :: status
    
    ! Local variables
    integer :: ncid_p, ncid_pet
    integer :: varid_tp, varid_pev, varid_lat, varid_lon, varid_time
    integer :: nlon, nlat, ntime
    real(real32), allocatable :: temp_data(:,:,:)
    integer :: i, j, t
    
    status = 0
    
    ! Open precipitation file
    status = nf90_open(precip_file, nf90_nowrite, ncid_p)
    if (status /= nf90_noerr) then
        print *, "ERROR: Failed to open precipitation file: ", trim(precip_file)
        print *, "NetCDF error: ", trim(nf90_strerror(status))
        return
    end if
    
    ! Open PET file
    status = nf90_open(pet_file, nf90_nowrite, ncid_pet)
    if (status /= nf90_noerr) then
        print *, "ERROR: Failed to open PET file: ", trim(pet_file)
        print *, "NetCDF error: ", trim(nf90_strerror(status))
        status = nf90_close(ncid_p)
        return
    end if
    
    ! Validate dimensions
    
    call validate_dimensions(ncid_p, nlon, nlat, ntime, status)
    if (status /= 0) then
        status = nf90_close(ncid_p)
        status = nf90_close(ncid_pet)
        return
    end if
    
    ! Check array sizes match file dimensions
    if (size(precipitation, 1) /= nlon .or. size(precipitation, 2) /= nlat .or. &
        size(precipitation, 3) /= ntime) then
        print *, "ERROR: Array dimensions mismatch"
        print *, "Expected: ", nlon, "x", nlat, "x", ntime
        print *, "Got: ", size(precipitation, 1), "x", size(precipitation, 2), "x", size(precipitation, 3)
        status = -1
        status = nf90_close(ncid_p)
        status = nf90_close(ncid_pet)
        return
    end if
    
    ! Read coordinates
    
    call read_coordinates(ncid_p, latitudes, longitudes, time_stamps, status)
    if (status /= 0) then
        status = nf90_close(ncid_p)
        status = nf90_close(ncid_pet)
        return
    end if
    
    ! Read precipitation data (tp)
    
    status = nf90_inq_varid(ncid_p, 'tp', varid_tp)
    if (status /= nf90_noerr) then
        print *, "ERROR: Variable 'tp' not found in precipitation file"
        status = nf90_close(ncid_p)
        status = nf90_close(ncid_pet)
        return
    end if
    
    allocate(temp_data(nlon, nlat, ntime))
    
    ! Read precipitation data
    status = nf90_get_var(ncid_p, varid_tp, temp_data, &
                         start=[1, 1, 1], count=[nlon, nlat, ntime])
    if (status /= nf90_noerr) then
        print *, "ERROR: Failed to read precipitation data"
        print *, "NetCDF error: ", trim(nf90_strerror(status))
        deallocate(temp_data)
        status = nf90_close(ncid_p)
        status = nf90_close(ncid_pet)
        return
    end if
    
    ! Convert precipitation: m/day → mm/month
    do t = 1, ntime
        do j = 1, nlat
            do i = 1, nlon
                if (temp_data(i,j,t) /= temp_data(i,j,t) .or. &  ! NaN check
                    abs(temp_data(i,j,t)) > 1e10) then            ! Large value check
                    precipitation(i,j,t) = FILL_VALUE
                else
                    precipitation(i,j,t) = real(temp_data(i,j,t), dp) * M_TO_MM * DAYS_PER_MONTH
                end if
            end do
        end do
    end do
    
    status = nf90_close(ncid_p)
    
    ! Read PET data (pev)
    
    status = nf90_inq_varid(ncid_pet, 'pev', varid_pev)
    if (status /= nf90_noerr) then
        print *, "ERROR: Variable 'pev' not found in PET file"
        deallocate(temp_data)
        status = nf90_close(ncid_pet)
        return
    end if
    
    ! Read PET data
    status = nf90_get_var(ncid_pet, varid_pev, temp_data, &
                         start=[1, 1, 1], count=[nlon, nlat, ntime])
    if (status /= nf90_noerr) then
        print *, "ERROR: Failed to read PET data"
        print *, "NetCDF error: ", trim(nf90_strerror(status))
        deallocate(temp_data)
        status = nf90_close(ncid_pet)
        return
    end if
    
    ! Convert PET: m/day → mm/month with sign flip
    do t = 1, ntime
        do j = 1, nlat
            do i = 1, nlon
                if (temp_data(i,j,t) /= temp_data(i,j,t) .or. &  ! NaN check
                    abs(temp_data(i,j,t)) > 1e10) then            ! Large value check
                    pet(i,j,t) = FILL_VALUE
                else
                    pet(i,j,t) = -real(temp_data(i,j,t), dp) * M_TO_MM * DAYS_PER_MONTH
                end if
            end do
        end do
    end do
    
    status = nf90_close(ncid_pet)
    deallocate(temp_data)
    
    ! Create and apply land mask
    
    call create_and_apply_land_mask(precipitation, pet, status)
    if (status /= 0) then
        print *, "ERROR: Failed to create land mask"
        return
    end if
    
    ! Validate loaded data
    call validate_climate_data(precipitation, pet, status)
    if (status /= 0) then
        print *, "ERROR: Climate data validation failed"
        return
    end if
    
    ! Apply land mask to remove ocean pixels
    print *, "Applying land mask to remove ocean pixels..."
    call create_and_apply_land_mask(precipitation, pet, status)
    if (status /= 0) then
        print *, "ERROR: Land masking failed"
        return
    end if
    
    print *, "ERA5 climate data loaded successfully with land masking applied"
    
end subroutine load_era5_climate_data

subroutine validate_dimensions(ncid, nlon, nlat, ntime, status)
    integer, intent(in) :: ncid
    integer, intent(out) :: nlon, nlat, ntime, status
    
    integer :: dimid
    
    status = 0
    
    ! Get longitude dimension
    status = nf90_inq_dimid(ncid, 'longitude', dimid)
    if (status /= nf90_noerr) then
        print *, "ERROR: 'longitude' dimension not found"
        return
    end if
    status = nf90_inquire_dimension(ncid, dimid, len=nlon)
    if (status /= nf90_noerr) return
    
    ! Get latitude dimension  
    status = nf90_inq_dimid(ncid, 'latitude', dimid)
    if (status /= nf90_noerr) then
        print *, "ERROR: 'latitude' dimension not found"
        return
    end if
    status = nf90_inquire_dimension(ncid, dimid, len=nlat)
    if (status /= nf90_noerr) return
    
    ! Get time dimension
    status = nf90_inq_dimid(ncid, 'valid_time', dimid)
    if (status /= nf90_noerr) then
        print *, "ERROR: 'valid_time' dimension not found"
        return
    end if
    status = nf90_inquire_dimension(ncid, dimid, len=ntime)
    if (status /= nf90_noerr) return
    
    print *, "  Dimensions: ", nlon, " lon ×", nlat, " lat ×", ntime, " time"
    
end subroutine validate_dimensions

subroutine read_coordinates(ncid, latitudes, longitudes, time_stamps, status)
    integer, intent(in) :: ncid
    real(dp), intent(out) :: latitudes(:), longitudes(:), time_stamps(:)
    integer, intent(out) :: status
    
    integer :: varid
    integer(int64), allocatable :: time_int(:)
    
    status = 0
    
    ! Read latitudes
    status = nf90_inq_varid(ncid, 'latitude', varid)
    if (status /= nf90_noerr) then
        print *, "ERROR: 'latitude' variable not found"
        return
    end if
    status = nf90_get_var(ncid, varid, latitudes)
    if (status /= nf90_noerr) then
        print *, "ERROR: Failed to read latitude data"
        return
    end if
    
    ! Read longitudes
    status = nf90_inq_varid(ncid, 'longitude', varid)
    if (status /= nf90_noerr) then
        print *, "ERROR: 'longitude' variable not found"
        return
    end if
    status = nf90_get_var(ncid, varid, longitudes)
    if (status /= nf90_noerr) then
        print *, "ERROR: Failed to read longitude data"
        return
    end if
    
    ! Read time stamps (int64 → real64)
    allocate(time_int(size(time_stamps)))
    status = nf90_inq_varid(ncid, 'valid_time', varid)
    if (status /= nf90_noerr) then
        print *, "ERROR: 'valid_time' variable not found"
        return
    end if
    status = nf90_get_var(ncid, varid, time_int)
    if (status /= nf90_noerr) then
        print *, "ERROR: Failed to read time data"
        deallocate(time_int)
        return
    end if
    
    time_stamps = real(time_int, dp)
    deallocate(time_int)
    
    print *, "  Latitude range: ", minval(latitudes), " to ", maxval(latitudes), "°N"
    print *, "  Longitude range: ", minval(longitudes), " to ", maxval(longitudes), "°E"
    print *, "  Time span: ", size(time_stamps), " months"
    
end subroutine read_coordinates

subroutine validate_climate_data(precipitation, pet, status)
    real(dp), intent(in) :: precipitation(:,:,:), pet(:,:,:)
    integer, intent(out) :: status
    
    real(dp) :: p_min, p_max, p_mean, pet_min, pet_max, pet_mean
    integer :: valid_p, valid_pet, total_points
    
    status = 0
    total_points = size(precipitation)
    
    ! Calculate statistics for valid precipitation data
    valid_p = count(precipitation > FILL_VALUE)
    if (valid_p > 0) then
        p_min = minval(precipitation, mask=precipitation > FILL_VALUE)
        p_max = maxval(precipitation, mask=precipitation > FILL_VALUE)
        p_mean = sum(precipitation, mask=precipitation > FILL_VALUE) / real(valid_p, dp)
    else
        p_min = FILL_VALUE
        p_max = FILL_VALUE  
        p_mean = FILL_VALUE
    end if
    
    ! Calculate statistics for valid PET data
    valid_pet = count(pet > FILL_VALUE)
    if (valid_pet > 0) then
        pet_min = minval(pet, mask=pet > FILL_VALUE)
        pet_max = maxval(pet, mask=pet > FILL_VALUE)
        pet_mean = sum(pet, mask=pet > FILL_VALUE) / real(valid_pet, dp)
    else
        pet_min = FILL_VALUE
        pet_max = FILL_VALUE
        pet_mean = FILL_VALUE
    end if
    
    ! Print validation summary
    print *, "  Precipitation validation:"
    print *, "    Valid points: ", valid_p, "/", total_points, &
             " (", 100.0 * valid_p / total_points, "%)"
    if (valid_p > 0) then
        print *, "    Range: ", p_min, " to ", p_max, " mm/month"
        print *, "    Mean: ", p_mean, " mm/month"
    end if
    
    print *, "  PET validation:"
    print *, "    Valid points: ", valid_pet, "/", total_points, &
             " (", 100.0 * valid_pet / total_points, "%)"
    if (valid_pet > 0) then
        print *, "    Range: ", pet_min, " to ", pet_max, " mm/month" 
        print *, "    Mean: ", pet_mean, " mm/month"
    end if
    
    ! Sanity checks for arid/semi-arid climate (Somaliland)
    if (valid_p > 0 .and. (p_mean < 0.0 .or. p_mean > 200.0)) then
        print *, "  ⚠ WARNING: Precipitation mean outside expected range (0-200 mm/month) for arid climate"
    end if
    
    if (valid_pet > 0 .and. (pet_mean < 60.0 .or. pet_mean > 500.0)) then
        print *, "  ⚠ WARNING: PET mean outside expected range (60-500 mm/month) for arid climate"
    end if
    
end subroutine validate_climate_data

subroutine save_processed_precipitation_data(precipitation, latitudes, longitudes, &
                                          output_file, status)
    real(dp), intent(in) :: precipitation(:,:,:)
    real(real32), intent(in) :: latitudes(:), longitudes(:)
    character(len=*), intent(in) :: output_file
    integer, intent(out) :: status
    
    integer :: ncid, dimid_lon, dimid_lat, dimid_time
    integer :: varid_lon, varid_lat, varid_time, varid_precip
    integer :: nlon, nlat, ntime, i
    
    nlon = size(precipitation, 1)
    nlat = size(precipitation, 2) 
    ntime = size(precipitation, 3)
    
    ! Create NetCDF file
    status = nf90_create(output_file, nf90_clobber, ncid)
    if (status /= nf90_noerr) then
        print *, "ERROR: Cannot create file ", trim(output_file)
        return
    end if
    
    ! Define dimensions
    status = nf90_def_dim(ncid, 'longitude', nlon, dimid_lon)
    status = nf90_def_dim(ncid, 'latitude', nlat, dimid_lat)
    status = nf90_def_dim(ncid, 'time', ntime, dimid_time)
    
    ! Define coordinate variables
    status = nf90_def_var(ncid, 'longitude', nf90_float, [dimid_lon], varid_lon)
    status = nf90_put_att(ncid, varid_lon, 'units', 'degrees_east')
    status = nf90_put_att(ncid, varid_lon, 'long_name', 'longitude')
    
    status = nf90_def_var(ncid, 'latitude', nf90_float, [dimid_lat], varid_lat)
    status = nf90_put_att(ncid, varid_lat, 'units', 'degrees_north')
    status = nf90_put_att(ncid, varid_lat, 'long_name', 'latitude')
    
    status = nf90_def_var(ncid, 'time', nf90_int, [dimid_time], varid_time)
    status = nf90_put_att(ncid, varid_time, 'units', 'months since 1981-01-01')
    
    ! Define precipitation variable
    status = nf90_def_var(ncid, 'precipitation', nf90_double, &
                         [dimid_lon, dimid_lat, dimid_time], varid_precip)
    status = nf90_put_att(ncid, varid_precip, 'units', 'mm/month')
    status = nf90_put_att(ncid, varid_precip, 'long_name', 'monthly precipitation')
    status = nf90_put_att(ncid, varid_precip, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_precip, 'standard_name', 'precipitation_amount')
    status = nf90_put_att(ncid, varid_precip, 'comment', 'ERA5-Land total precipitation converted to mm/month')
    
    ! Global attributes
    status = nf90_put_att(ncid, nf90_global, 'title', 'Processed ERA5-Land Precipitation for Somaliland')
    status = nf90_put_att(ncid, nf90_global, 'institution', 'University of Glasgow')
    status = nf90_put_att(ncid, nf90_global, 'source', 'ERA5-Land total precipitation with unit conversion')
    status = nf90_put_att(ncid, nf90_global, 'creator_name', 'Khadar Daahir')
    status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.7')
    status = nf90_put_att(ncid, nf90_global, 'date_created', '2025-08-26')
    status = nf90_put_att(ncid, nf90_global, 'summary', 'Monthly precipitation data for Somaliland drought analysis')
    status = nf90_put_att(ncid, nf90_global, 'keywords', 'precipitation, climate, ERA5-Land, Somaliland, Horn of Africa')
    status = nf90_put_att(ncid, nf90_global, 'references', 'Muñoz-Sabater et al. (2021), doi:10.5194/essd-13-4349-2021')
    
    status = nf90_enddef(ncid)
    
    ! Write data
    status = nf90_put_var(ncid, varid_lon, longitudes)
    status = nf90_put_var(ncid, varid_lat, latitudes)
    status = nf90_put_var(ncid, varid_time, [(i, i=1,ntime)])
    status = nf90_put_var(ncid, varid_precip, precipitation)
    
    status = nf90_close(ncid)
    
end subroutine save_processed_precipitation_data

subroutine save_processed_pet_data(pet, latitudes, longitudes, output_file, status)
    real(dp), intent(in) :: pet(:,:,:)
    real(real32), intent(in) :: latitudes(:), longitudes(:)
    character(len=*), intent(in) :: output_file
    integer, intent(out) :: status
    
    integer :: ncid, dimid_lon, dimid_lat, dimid_time
    integer :: varid_lon, varid_lat, varid_time, varid_pet
    integer :: nlon, nlat, ntime, i
    
    nlon = size(pet, 1)
    nlat = size(pet, 2)
    ntime = size(pet, 3)
    
    ! Create NetCDF file
    status = nf90_create(output_file, nf90_clobber, ncid)
    if (status /= nf90_noerr) then
        print *, "ERROR: Cannot create file ", trim(output_file)
        return
    end if
    
    ! Define dimensions
    status = nf90_def_dim(ncid, 'longitude', nlon, dimid_lon)
    status = nf90_def_dim(ncid, 'latitude', nlat, dimid_lat)
    status = nf90_def_dim(ncid, 'time', ntime, dimid_time)
    
    ! Define coordinate variables
    status = nf90_def_var(ncid, 'longitude', nf90_float, [dimid_lon], varid_lon)
    status = nf90_put_att(ncid, varid_lon, 'units', 'degrees_east')
    status = nf90_put_att(ncid, varid_lon, 'long_name', 'longitude')
    
    status = nf90_def_var(ncid, 'latitude', nf90_float, [dimid_lat], varid_lat)
    status = nf90_put_att(ncid, varid_lat, 'units', 'degrees_north')
    status = nf90_put_att(ncid, varid_lat, 'long_name', 'latitude')
    
    status = nf90_def_var(ncid, 'time', nf90_int, [dimid_time], varid_time)
    status = nf90_put_att(ncid, varid_time, 'units', 'months since 1981-01-01')
    
    ! Define PET variable
    status = nf90_def_var(ncid, 'potential_evapotranspiration', nf90_double, &
                         [dimid_lon, dimid_lat, dimid_time], varid_pet)
    status = nf90_put_att(ncid, varid_pet, 'units', 'mm/month')
    status = nf90_put_att(ncid, varid_pet, 'long_name', 'monthly potential evapotranspiration')
    status = nf90_put_att(ncid, varid_pet, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_pet, 'standard_name', 'water_potential_evaporation_amount')
    status = nf90_put_att(ncid, varid_pet, 'comment', 'ERA5-Land potential evaporation converted to mm/month')
    
    ! Global attributes
    status = nf90_put_att(ncid, nf90_global, 'title', 'Processed ERA5-Land Potential Evapotranspiration for Somaliland')
    status = nf90_put_att(ncid, nf90_global, 'institution', 'University of Glasgow')
    status = nf90_put_att(ncid, nf90_global, 'source', 'ERA5-Land potential evaporation with unit conversion')
    status = nf90_put_att(ncid, nf90_global, 'creator_name', 'Khadar Daahir')
    status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.7')
    status = nf90_put_att(ncid, nf90_global, 'date_created', '2025-08-26')
    status = nf90_put_att(ncid, nf90_global, 'summary', 'Monthly potential evapotranspiration data for Somaliland drought analysis')
    status = nf90_put_att(ncid, nf90_global, 'keywords', 'evapotranspiration, climate, ERA5-Land, Somaliland, Horn of Africa')
    status = nf90_put_att(ncid, nf90_global, 'references', 'Muñoz-Sabater et al. (2021), doi:10.5194/essd-13-4349-2021')
    
    status = nf90_enddef(ncid)
    
    ! Write data
    status = nf90_put_var(ncid, varid_lon, longitudes)
    status = nf90_put_var(ncid, varid_lat, latitudes)
    status = nf90_put_var(ncid, varid_time, [(i, i=1,ntime)])
    status = nf90_put_var(ncid, varid_pet, pet)
    
    status = nf90_close(ncid)
    
end subroutine save_processed_pet_data

subroutine create_and_apply_land_mask(precipitation, pet, status)
    real(dp), intent(inout) :: precipitation(:,:,:), pet(:,:,:)
    integer, intent(out) :: status
    
    logical, allocatable :: land_mask(:,:)
    integer :: nlon, nlat, ntime, i, j, t
    integer :: valid_p_count, valid_pet_count
    real(dp) :: valid_fraction_p, valid_fraction_pet
    integer :: total_masked, original_valid_p, original_valid_pet
    
    nlon = size(precipitation, 1)
    nlat = size(precipitation, 2)
    ntime = size(precipitation, 3)
    
    allocate(land_mask(nlon, nlat))
    land_mask = .false.
    
    ! Count original valid points
    original_valid_p = count(precipitation > FILL_VALUE)
    original_valid_pet = count(pet > FILL_VALUE)
    
    print *, "  Creating land mask..."
    
    ! Create land mask based on data availability
    ! A pixel is considered "land" if it has valid data for a sufficient fraction of time steps
    do j = 1, nlat
        do i = 1, nlon
            valid_p_count = count(precipitation(i,j,:) > FILL_VALUE)
            valid_pet_count = count(pet(i,j,:) > FILL_VALUE)
            
            valid_fraction_p = real(valid_p_count, dp) / real(ntime, dp)
            valid_fraction_pet = real(valid_pet_count, dp) / real(ntime, dp)
            
            ! Pixel is land if both P and PET have sufficient valid data
            if (valid_fraction_p >= LAND_THRESHOLD .and. valid_fraction_pet >= LAND_THRESHOLD) then
                land_mask(i,j) = .true.
            end if
        end do
    end do
    
    ! Apply land mask - set ocean/invalid pixels to fill value
    total_masked = 0
    do j = 1, nlat
        do i = 1, nlon
            if (.not. land_mask(i,j)) then
                ! Mask all time steps for this pixel
                do t = 1, ntime
                    if (precipitation(i,j,t) > FILL_VALUE) then
                        precipitation(i,j,t) = FILL_VALUE
                        total_masked = total_masked + 1
                    end if
                    if (pet(i,j,t) > FILL_VALUE) then
                        pet(i,j,t) = FILL_VALUE
                    end if
                end do
            end if
        end do
    end do
    
    print *, "    Land pixels identified: ", count(land_mask), " out of ", nlon*nlat, &
             " (", 100.0*count(land_mask)/(nlon*nlat), "%)"
    print *, "    Ocean pixels masked: ", total_masked, " data points"
    print *, "    Final valid P points: ", count(precipitation > FILL_VALUE), &
             " (was ", original_valid_p, ")"
    print *, "    Final valid PET points: ", count(pet > FILL_VALUE), &
             " (was ", original_valid_pet, ")"
    
    deallocate(land_mask)
    status = 0
    
end subroutine create_and_apply_land_mask

!> Apply pre-computed land mask to data
subroutine apply_land_mask(precipitation, pet, land_mask, status)
    real(dp), intent(inout) :: precipitation(:,:,:), pet(:,:,:)
    logical, intent(in) :: land_mask(:,:)
    integer, intent(out) :: status
    
    integer :: i, j, t, nlon, nlat, ntime
    integer :: original_valid_p, original_valid_pet, masked_valid_p, masked_valid_pet
    
    nlon = size(precipitation, 1)
    nlat = size(precipitation, 2)
    ntime = size(precipitation, 3)
    
    ! Count original valid points
    original_valid_p = count(precipitation > FILL_VALUE)
    original_valid_pet = count(pet > FILL_VALUE)
    
    ! Apply land mask - set ocean pixels to fill value
    do j = 1, nlat
        do i = 1, nlon
            if (.not. land_mask(i,j)) then
                precipitation(i,j,:) = FILL_VALUE
                pet(i,j,:) = FILL_VALUE
            end if
        end do
    end do
    
    ! Count masked valid points
    masked_valid_p = count(precipitation > FILL_VALUE)
    masked_valid_pet = count(pet > FILL_VALUE)
    
    print *, "  Land masking applied:"
    print *, "    Land pixels: ", count(land_mask), "/", size(land_mask), &
             " (", 100.0 * count(land_mask) / size(land_mask), "%)"
    print *, "    Valid P points: ", masked_valid_p, " (was ", original_valid_p, ")"
    print *, "    Valid PET points: ", masked_valid_pet, " (was ", original_valid_pet, ")"
    
    status = 0
    
end subroutine apply_land_mask

!> Create land mask for external use
subroutine create_land_mask(precipitation, pet, land_mask, status)
    real(dp), intent(in) :: precipitation(:,:,:), pet(:,:,:)
    logical, intent(out) :: land_mask(:,:)
    integer, intent(out) :: status
    
    integer :: nlon, nlat, ntime, i, j
    integer :: valid_p_count, valid_pet_count
    real(dp) :: valid_fraction_p, valid_fraction_pet
    
    nlon = size(precipitation, 1)
    nlat = size(precipitation, 2)
    ntime = size(precipitation, 3)
    
    status = 0
    land_mask = .false.
    
    ! Check array dimensions
    if (size(land_mask, 1) /= nlon .or. size(land_mask, 2) /= nlat) then
        print *, "ERROR: Land mask dimensions mismatch"
        status = -1
        return
    end if
    
    ! Create land mask based on data availability
    do j = 1, nlat
        do i = 1, nlon
            valid_p_count = count(precipitation(i,j,:) > FILL_VALUE)
            valid_pet_count = count(pet(i,j,:) > FILL_VALUE)
            
            valid_fraction_p = real(valid_p_count, dp) / real(ntime, dp)
            valid_fraction_pet = real(valid_pet_count, dp) / real(ntime, dp)
            
            ! Pixel is land if both P and PET have sufficient valid data
            if (valid_fraction_p >= LAND_THRESHOLD .and. valid_fraction_pet >= LAND_THRESHOLD) then
                land_mask(i,j) = .true.
            end if
        end do
    end do
    
end subroutine create_land_mask

!> Save all SPI timescales to one NetCDF file
subroutine save_spi_multiscale(spi_1, spi_3, spi_6, spi_12, latitudes, longitudes, output_file, status)
    real(dp), intent(in) :: spi_1(:,:,:), spi_3(:,:,:), spi_6(:,:,:), spi_12(:,:,:)
    real(dp), intent(in) :: latitudes(:), longitudes(:)
    character(len=*), intent(in) :: output_file
    integer, intent(out) :: status
    
    integer :: ncid, dimid_lon, dimid_lat, dimid_time
    integer :: varid_lon, varid_lat, varid_time
    integer :: varid_spi1, varid_spi3, varid_spi6, varid_spi12
    integer :: nlon, nlat, ntime
    
    nlon = size(spi_1, 1)
    nlat = size(spi_1, 2)
    ntime = size(spi_1, 3)
    
    status = 0
    
    ! Create NetCDF file
    status = nf90_create(output_file, nf90_clobber, ncid)
    if (status /= nf90_noerr) return
    
    ! Define dimensions
    status = nf90_def_dim(ncid, 'longitude', nlon, dimid_lon)
    status = nf90_def_dim(ncid, 'latitude', nlat, dimid_lat)
    status = nf90_def_dim(ncid, 'time', ntime, dimid_time)
    
    ! Define coordinate variables
    status = nf90_def_var(ncid, 'longitude', nf90_float, dimid_lon, varid_lon)
    status = nf90_put_att(ncid, varid_lon, 'units', 'degrees_east')
    status = nf90_put_att(ncid, varid_lon, 'long_name', 'longitude')
    
    status = nf90_def_var(ncid, 'latitude', nf90_float, dimid_lat, varid_lat)
    status = nf90_put_att(ncid, varid_lat, 'units', 'degrees_north')
    status = nf90_put_att(ncid, varid_lat, 'long_name', 'latitude')
    
    status = nf90_def_var(ncid, 'time', nf90_int, dimid_time, varid_time)
    status = nf90_put_att(ncid, varid_time, 'units', 'months since 1981-01-01')
    
    ! Define SPI variables (Fortran order: lon, lat, time)
    status = nf90_def_var(ncid, 'spi_1', nf90_double, [dimid_lon, dimid_lat, dimid_time], varid_spi1)
    status = nf90_put_att(ncid, varid_spi1, 'units', '1')
    status = nf90_put_att(ncid, varid_spi1, 'long_name', '1-month Standardized Precipitation Index')
    status = nf90_put_att(ncid, varid_spi1, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_spi1, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, varid_spi1, 'standard_name', 'precipitation_index')
    status = nf90_put_att(ncid, varid_spi1, 'comment', 'Standardized using gamma distribution fitting (WMO guidelines)')
    
    status = nf90_def_var(ncid, 'spi_3', nf90_double, [dimid_lon, dimid_lat, dimid_time], varid_spi3)
    status = nf90_put_att(ncid, varid_spi3, 'units', '1')
    status = nf90_put_att(ncid, varid_spi3, 'long_name', '3-month Standardized Precipitation Index')
    status = nf90_put_att(ncid, varid_spi3, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_spi3, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, varid_spi3, 'standard_name', 'precipitation_index')
    status = nf90_put_att(ncid, varid_spi3, 'comment', 'First 2 months are undefined (set to FillValue). Standardized using gamma distribution fitting (WMO guidelines)')
    
    status = nf90_def_var(ncid, 'spi_6', nf90_double, [dimid_lon, dimid_lat, dimid_time], varid_spi6)
    status = nf90_put_att(ncid, varid_spi6, 'units', '1')
    status = nf90_put_att(ncid, varid_spi6, 'long_name', '6-month Standardized Precipitation Index')
    status = nf90_put_att(ncid, varid_spi6, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_spi6, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, varid_spi6, 'standard_name', 'precipitation_index')
    status = nf90_put_att(ncid, varid_spi6, 'comment', 'First 5 months are undefined (set to FillValue). Standardized using gamma distribution fitting (WMO guidelines)')
    
    status = nf90_def_var(ncid, 'spi_12', nf90_double, [dimid_lon, dimid_lat, dimid_time], varid_spi12)
    status = nf90_put_att(ncid, varid_spi12, 'units', '1')
    status = nf90_put_att(ncid, varid_spi12, 'long_name', '12-month Standardized Precipitation Index')
    status = nf90_put_att(ncid, varid_spi12, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_spi12, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, varid_spi12, 'standard_name', 'precipitation_index')
    status = nf90_put_att(ncid, varid_spi12, 'comment', 'First 11 months are undefined (set to FillValue). Standardized using gamma distribution fitting (WMO guidelines)')
    
    ! Global attributes
    status = nf90_put_att(ncid, nf90_global, 'title', 'Standardized Precipitation Index (SPI) for Somaliland')
    status = nf90_put_att(ncid, nf90_global, 'institution', 'University of Glasgow')
    status = nf90_put_att(ncid, nf90_global, 'source', 'ERA5-Land precipitation with FSML statistical processing')
    status = nf90_put_att(ncid, nf90_global, 'creator_name', 'Khadar Daahir')
    status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.7')
    status = nf90_put_att(ncid, nf90_global, 'date_created', '2025-08-26')
    status = nf90_put_att(ncid, nf90_global, 'summary', 'Multiscale SPI drought indices for Somaliland calculated using WMO-compliant gamma distribution fitting and FSML statistical library')
    status = nf90_put_att(ncid, nf90_global, 'keywords', 'drought, SPI, standardized precipitation index, climate, Somaliland, Horn of Africa')
    status = nf90_put_att(ncid, nf90_global, 'geospatial_lat_min', 8.0_dp)
    status = nf90_put_att(ncid, nf90_global, 'geospatial_lat_max', 11.5_dp)
    status = nf90_put_att(ncid, nf90_global, 'geospatial_lon_min', 43.0_dp)
    status = nf90_put_att(ncid, nf90_global, 'geospatial_lon_max', 48.5_dp)
    status = nf90_put_att(ncid, nf90_global, 'time_coverage_start', '1981-01-01')
    status = nf90_put_att(ncid, nf90_global, 'time_coverage_end', '2024-12-31')
    status = nf90_put_att(ncid, nf90_global, 'methodology', 'Gamma distribution fitting followed by normal standardization (WMO SPI guidelines)')
    
    status = nf90_enddef(ncid)
    
    ! Write data
    status = nf90_put_var(ncid, varid_lon, longitudes)
    status = nf90_put_var(ncid, varid_lat, latitudes)
    status = nf90_put_var(ncid, varid_spi1, spi_1)
    if (status /= nf90_noerr) print *, "ERROR writing SPI-1: ", trim(nf90_strerror(status))
    status = nf90_put_var(ncid, varid_spi3, spi_3)
    if (status /= nf90_noerr) print *, "ERROR writing SPI-3: ", trim(nf90_strerror(status))
    status = nf90_put_var(ncid, varid_spi6, spi_6)
    if (status /= nf90_noerr) print *, "ERROR writing SPI-6: ", trim(nf90_strerror(status))
    status = nf90_put_var(ncid, varid_spi12, spi_12)
    if (status /= nf90_noerr) print *, "ERROR writing SPI-12: ", trim(nf90_strerror(status))
    
    status = nf90_close(ncid)
    
end subroutine save_spi_multiscale

!> Save individual SPI timescale to NetCDF file
subroutine save_spi_single_scale(spi_data, timescale, latitudes, longitudes, output_file, status)
    real(dp), intent(in) :: spi_data(:,:,:)
    integer, intent(in) :: timescale
    real(dp), intent(in) :: latitudes(:), longitudes(:)
    character(len=*), intent(in) :: output_file
    integer, intent(out) :: status
    
    integer :: ncid, dimid_lon, dimid_lat, dimid_time
    integer :: varid_lon, varid_lat, varid_time, varid_spi
    integer :: nlon, nlat, ntime
    character(len=100) :: var_name, long_name
    character(len=200) :: comment_text  ! Separate variable for comments
    
    nlon = size(spi_data, 1)
    nlat = size(spi_data, 2)
    ntime = size(spi_data, 3)
    
    status = 0
    
    ! Create variable names
    write(var_name, '("spi_", I0)') timescale
    write(long_name, '(I0, "-month Standardized Precipitation Index")') timescale
    
    ! Create NetCDF file
    status = nf90_create(output_file, nf90_clobber, ncid)
    if (status /= nf90_noerr) return
    
    ! Define dimensions
    status = nf90_def_dim(ncid, 'longitude', nlon, dimid_lon)
    status = nf90_def_dim(ncid, 'latitude', nlat, dimid_lat)
    status = nf90_def_dim(ncid, 'time', ntime, dimid_time)
    
    ! Define coordinate variables
    status = nf90_def_var(ncid, 'longitude', nf90_float, dimid_lon, varid_lon)
    status = nf90_put_att(ncid, varid_lon, 'units', 'degrees_east')
    status = nf90_put_att(ncid, varid_lon, 'long_name', 'longitude')
    
    status = nf90_def_var(ncid, 'latitude', nf90_float, dimid_lat, varid_lat)
    status = nf90_put_att(ncid, varid_lat, 'units', 'degrees_north')
    status = nf90_put_att(ncid, varid_lat, 'long_name', 'latitude')
    
    status = nf90_def_var(ncid, 'time', nf90_int, dimid_time, varid_time)
    status = nf90_put_att(ncid, varid_time, 'units', 'months since 1981-01-01')
    
    ! Define SPI variable (Fortran order: lon, lat, time)
    status = nf90_def_var(ncid, trim(var_name), nf90_double, [dimid_lon, dimid_lat, dimid_time], varid_spi)
    status = nf90_put_att(ncid, varid_spi, 'units', '1')
    status = nf90_put_att(ncid, varid_spi, 'long_name', trim(long_name))
    status = nf90_put_att(ncid, varid_spi, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_spi, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, varid_spi, 'standard_name', 'precipitation_index')
    
    ! Add timescale-specific comment
    if (timescale > 1) then
        write(comment_text, '("First ", I0, " months are undefined (set to FillValue). Standardized using gamma distribution fitting (WMO guidelines)")') timescale-1
        status = nf90_put_att(ncid, varid_spi, 'comment', trim(comment_text))
    else
        status = nf90_put_att(ncid, varid_spi, 'comment', 'Standardized using gamma distribution fitting (WMO guidelines)')
    end if
    
    ! Global attributes
    write(long_name, '(I0, "-month Standardized Precipitation Index")') timescale  ! Reset long_name
    status = nf90_put_att(ncid, nf90_global, 'title', trim(long_name)//' for Somaliland')
    status = nf90_put_att(ncid, nf90_global, 'institution', 'University of Glasgow')
    status = nf90_put_att(ncid, nf90_global, 'source', 'ERA5-Land precipitation with FSML statistical processing')
    status = nf90_put_att(ncid, nf90_global, 'creator_name', 'Khadar Daahir')
    status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.7')
    status = nf90_put_att(ncid, nf90_global, 'date_created', '2025-08-26')
    status = nf90_put_att(ncid, nf90_global, 'summary', trim(long_name)//' drought index for Somaliland calculated using WMO-compliant gamma distribution fitting')
    status = nf90_put_att(ncid, nf90_global, 'keywords', 'drought, SPI, standardized precipitation index, climate, Somaliland, Horn of Africa')
    status = nf90_put_att(ncid, nf90_global, 'methodology', 'Gamma distribution fitting followed by normal standardization (WMO SPI guidelines)')
    
    status = nf90_enddef(ncid)
    
    ! Write data
    status = nf90_put_var(ncid, varid_lon, longitudes)
    status = nf90_put_var(ncid, varid_lat, latitudes)
    status = nf90_put_var(ncid, varid_spi, spi_data)
    
    status = nf90_close(ncid)
    
end subroutine save_spi_single_scale

!> Save all SPEI timescales to NetCDF file
subroutine save_spei_multiscale(spei_1, spei_3, spei_6, spei_12, latitudes, longitudes, output_file, status)
    real(dp), intent(in) :: spei_1(:,:,:), spei_3(:,:,:), spei_6(:,:,:), spei_12(:,:,:)
    real(dp), intent(in) :: latitudes(:), longitudes(:)
    character(len=*), intent(in) :: output_file
    integer, intent(out) :: status
    
    integer :: ncid, lat_dimid, lon_dimid, time_dimid
    integer :: lat_varid, lon_varid, time_varid
    integer :: spei_1_varid, spei_3_varid, spei_6_varid, spei_12_varid
    integer :: dims(3)
    integer :: nlon, nlat, ntime
    integer :: i
    
    nlon = size(longitudes)
    nlat = size(latitudes) 
    ntime = size(spei_1, 3)
    
    ! Create NetCDF file
    status = nf90_create(output_file, NF90_CLOBBER, ncid)
    if (status /= nf90_noerr) return
    
    ! Define dimensions
    status = nf90_def_dim(ncid, "longitude", nlon, lon_dimid)
    if (status /= nf90_noerr) return
    status = nf90_def_dim(ncid, "latitude", nlat, lat_dimid)
    if (status /= nf90_noerr) return
    status = nf90_def_dim(ncid, "time", ntime, time_dimid)
    if (status /= nf90_noerr) return
    
    ! Define coordinate variables
    status = nf90_def_var(ncid, 'longitude', nf90_float, lon_dimid, lon_varid)
    status = nf90_put_att(ncid, lon_varid, 'units', 'degrees_east')
    status = nf90_put_att(ncid, lon_varid, 'long_name', 'longitude')
    
    status = nf90_def_var(ncid, 'latitude', nf90_float, lat_dimid, lat_varid)
    status = nf90_put_att(ncid, lat_varid, 'units', 'degrees_north')
    status = nf90_put_att(ncid, lat_varid, 'long_name', 'latitude')
    
    status = nf90_def_var(ncid, 'time', nf90_int, time_dimid, time_varid)
    status = nf90_put_att(ncid, time_varid, 'units', 'months since 1981-01-01')
    
    ! Define SPEI variables
    dims = [lon_dimid, lat_dimid, time_dimid]
    status = nf90_def_var(ncid, "spei_1", nf90_double, dims, spei_1_varid)
    status = nf90_def_var(ncid, "spei_3", nf90_double, dims, spei_3_varid)
    status = nf90_def_var(ncid, "spei_6", nf90_double, dims, spei_6_varid)
    status = nf90_def_var(ncid, "spei_12", nf90_double, dims, spei_12_varid)
    
    ! Add attributes for SPEI-1
    status = nf90_put_att(ncid, spei_1_varid, 'units', '1')
    status = nf90_put_att(ncid, spei_1_varid, "long_name", "1-month Standardized Precipitation-Evapotranspiration Index")
    status = nf90_put_att(ncid, spei_1_varid, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, spei_1_varid, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, spei_1_varid, 'standard_name', 'precipitation_evapotranspiration_index')
    status = nf90_put_att(ncid, spei_1_varid, 'comment', 'Standardized using gamma distribution fitting with data shifting')
    
    ! Add attributes for SPEI-3
    status = nf90_put_att(ncid, spei_3_varid, 'units', '1')
    status = nf90_put_att(ncid, spei_3_varid, "long_name", "3-month Standardized Precipitation-Evapotranspiration Index")
    status = nf90_put_att(ncid, spei_3_varid, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, spei_3_varid, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, spei_3_varid, 'standard_name', 'precipitation_evapotranspiration_index')
    status = nf90_put_att(ncid, spei_3_varid, 'comment', 'First 2 months are undefined (set to FillValue). Standardized using gamma distribution fitting with data shifting')
    
    ! Add attributes for SPEI-6
    status = nf90_put_att(ncid, spei_6_varid, 'units', '1')
    status = nf90_put_att(ncid, spei_6_varid, "long_name", "6-month Standardized Precipitation-Evapotranspiration Index")
    status = nf90_put_att(ncid, spei_6_varid, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, spei_6_varid, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, spei_6_varid, 'standard_name', 'precipitation_evapotranspiration_index')
    status = nf90_put_att(ncid, spei_6_varid, 'comment', 'First 5 months are undefined (set to FillValue). Standardized using gamma distribution fitting with data shifting')
    
    ! Add attributes for SPEI-12
    status = nf90_put_att(ncid, spei_12_varid, 'units', '1')
    status = nf90_put_att(ncid, spei_12_varid, "long_name", "12-month Standardized Precipitation-Evapotranspiration Index")
    status = nf90_put_att(ncid, spei_12_varid, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, spei_12_varid, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, spei_12_varid, 'standard_name', 'precipitation_evapotranspiration_index')
    status = nf90_put_att(ncid, spei_12_varid, 'comment', 'First 11 months are undefined (set to FillValue). Standardized using gamma distribution fitting with data shifting')
    
    ! Global attributes
    status = nf90_put_att(ncid, nf90_global, 'title', 'Multi-scale Standardized Precipitation-Evapotranspiration Index for Somaliland')
    status = nf90_put_att(ncid, nf90_global, 'institution', 'University of Glasgow')
    status = nf90_put_att(ncid, nf90_global, 'source', 'ERA5-Land precipitation and PET with FSML statistical processing')
    status = nf90_put_att(ncid, nf90_global, 'creator_name', 'Khadar Daahir')
    status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.7')
    status = nf90_put_att(ncid, nf90_global, 'date_created', '2025-08-26')
    status = nf90_put_att(ncid, nf90_global, 'summary', 'Multi-scale SPEI drought indices for Somaliland calculated using gamma distribution fitting with data shifting')
    status = nf90_put_att(ncid, nf90_global, 'keywords', 'drought, SPEI, standardized precipitation evapotranspiration index, climate, Somaliland, Horn of Africa')
    status = nf90_put_att(ncid, nf90_global, 'methodology', 'Gamma distribution fitting with data shifting followed by normal standardization')
    
    ! End definition mode
    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) return
    
    ! Write coordinate data
    status = nf90_put_var(ncid, lon_varid, longitudes)
    status = nf90_put_var(ncid, lat_varid, latitudes)
    status = nf90_put_var(ncid, time_varid, [(i, i=1,ntime)])
    
    ! Write SPEI data
    status = nf90_put_var(ncid, spei_1_varid, spei_1)
    if (status /= nf90_noerr) print *, "ERROR writing SPEI-1: ", trim(nf90_strerror(status))
    
    status = nf90_put_var(ncid, spei_3_varid, spei_3)
    if (status /= nf90_noerr) print *, "ERROR writing SPEI-3: ", trim(nf90_strerror(status))
    
    status = nf90_put_var(ncid, spei_6_varid, spei_6)
    if (status /= nf90_noerr) print *, "ERROR writing SPEI-6: ", trim(nf90_strerror(status))
    
    status = nf90_put_var(ncid, spei_12_varid, spei_12)
    if (status /= nf90_noerr) print *, "ERROR writing SPEI-12: ", trim(nf90_strerror(status))
    
    status = nf90_close(ncid)
    
end subroutine save_spei_multiscale

!> Save individual SPEI timescale to NetCDF file
subroutine save_spei_single_scale(spei_data, timescale, latitudes, longitudes, output_file, status)
    real(dp), intent(in) :: spei_data(:,:,:)
    integer, intent(in) :: timescale
    real(dp), intent(in) :: latitudes(:), longitudes(:)
    character(len=*), intent(in) :: output_file
    integer, intent(out) :: status
    
    integer :: ncid, dimid_lat, dimid_lon, dimid_time
    integer :: varid_lat, varid_lon, varid_time, varid_spei
    integer :: nlon, nlat, ntime
    integer :: i
    character(len=80) :: var_name, long_name
    character(len=200) :: comment_text
    
    nlon = size(longitudes)
    nlat = size(latitudes)
    ntime = size(spei_data, 3)
    
    ! Create variable names
    write(var_name, '("spei_", I0)') timescale
    write(long_name, '(I0, "-month Standardized Precipitation-Evapotranspiration Index")') timescale
    
    ! Create NetCDF file
    status = nf90_create(output_file, nf90_clobber, ncid)
    if (status /= nf90_noerr) return
    
    ! Define dimensions
    status = nf90_def_dim(ncid, 'longitude', nlon, dimid_lon)
    status = nf90_def_dim(ncid, 'latitude', nlat, dimid_lat)
    status = nf90_def_dim(ncid, 'time', ntime, dimid_time)
    
    ! Define coordinate variables
    status = nf90_def_var(ncid, 'longitude', nf90_float, dimid_lon, varid_lon)
    status = nf90_put_att(ncid, varid_lon, 'units', 'degrees_east')
    status = nf90_put_att(ncid, varid_lon, 'long_name', 'longitude')
    
    status = nf90_def_var(ncid, 'latitude', nf90_float, dimid_lat, varid_lat)
    status = nf90_put_att(ncid, varid_lat, 'units', 'degrees_north')
    status = nf90_put_att(ncid, varid_lat, 'long_name', 'latitude')
    
    status = nf90_def_var(ncid, 'time', nf90_int, dimid_time, varid_time)
    status = nf90_put_att(ncid, varid_time, 'units', 'months since 1981-01-01')
    
    ! Define SPEI variable (Fortran order: lon, lat, time)
    status = nf90_def_var(ncid, trim(var_name), nf90_double, [dimid_lon, dimid_lat, dimid_time], varid_spei)
    status = nf90_put_att(ncid, varid_spei, 'units', '1')
    status = nf90_put_att(ncid, varid_spei, 'long_name', trim(long_name))
    status = nf90_put_att(ncid, varid_spei, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_spei, 'valid_range', (/-3.5_dp, 3.5_dp/))
    status = nf90_put_att(ncid, varid_spei, 'standard_name', 'precipitation_evapotranspiration_index')
    
    ! Add timescale-specific comment
    if (timescale > 1) then
        write(comment_text, '("First ", I0, " months are undefined (set to FillValue). Standardized using gamma distribution fitting with data shifting")') timescale-1
        status = nf90_put_att(ncid, varid_spei, 'comment', trim(comment_text))
    else
        status = nf90_put_att(ncid, varid_spei, 'comment', 'Standardized using gamma distribution fitting with data shifting')
    end if
    
    ! Global attributes
    write(long_name, '(I0, "-month Standardized Precipitation-Evapotranspiration Index")') timescale  ! Reset long_name
    status = nf90_put_att(ncid, nf90_global, 'title', trim(long_name)//' for Somaliland')
    status = nf90_put_att(ncid, nf90_global, 'institution', 'University of Glasgow')
    status = nf90_put_att(ncid, nf90_global, 'source', 'ERA5-Land precipitation and PET with FSML statistical processing')
    status = nf90_put_att(ncid, nf90_global, 'creator_name', 'Khadar Daahir')
    status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.7')
    status = nf90_put_att(ncid, nf90_global, 'date_created', '2025-08-26')
    status = nf90_put_att(ncid, nf90_global, 'summary', trim(long_name)//' drought index for Somaliland calculated using gamma distribution fitting with data shifting')
    status = nf90_put_att(ncid, nf90_global, 'keywords', 'drought, SPEI, standardized precipitation evapotranspiration index, climate, Somaliland, Horn of Africa')
    status = nf90_put_att(ncid, nf90_global, 'methodology', 'Gamma distribution fitting with data shifting followed by normal standardization')
    
    status = nf90_enddef(ncid)
    
    ! Write data
    status = nf90_put_var(ncid, varid_lon, longitudes)
    status = nf90_put_var(ncid, varid_lat, latitudes)
    status = nf90_put_var(ncid, varid_time, [(i, i=1,ntime)])
    status = nf90_put_var(ncid, varid_spei, spei_data)
    
    status = nf90_close(ncid)
    
end subroutine save_spei_single_scale

!> Export bias-corrected CMIP6 precipitation
subroutine save_corrected_cmip6_precipitation(precipitation, latitudes, longitudes, &
                                            scenario, model, start_year, end_year, &
                                            output_file, status, bias_factor, correlation, rmse)
    real(dp), intent(in) :: precipitation(:,:,:)
    real(dp), intent(in) :: latitudes(:), longitudes(:)
    character(len=*), intent(in) :: scenario, model
    integer, intent(in) :: start_year, end_year
    character(len=*), intent(in) :: output_file
    integer, intent(out) :: status
    real(dp), intent(in), optional :: bias_factor, correlation, rmse
    
    integer :: ncid, dimid_lon, dimid_lat, dimid_time
    integer :: varid_lon, varid_lat, varid_time, varid_precip
    integer :: nlon, nlat, ntime, i
    character(len=50) :: time_units
    
    nlon = size(precipitation, 1)
    nlat = size(precipitation, 2) 
    ntime = size(precipitation, 3)
    
    ! Create NetCDF file
    status = nf90_create(output_file, nf90_clobber, ncid)
    if (status /= nf90_noerr) then
        print *, "ERROR: Cannot create file ", trim(output_file)
        return
    end if
    
    ! Define dimensions
    status = nf90_def_dim(ncid, 'longitude', nlon, dimid_lon)
    status = nf90_def_dim(ncid, 'latitude', nlat, dimid_lat)
    status = nf90_def_dim(ncid, 'time', ntime, dimid_time)
    
    ! Define coordinate variables
    status = nf90_def_var(ncid, 'longitude', nf90_double, [dimid_lon], varid_lon)
    status = nf90_put_att(ncid, varid_lon, 'units', 'degrees_east')
    status = nf90_put_att(ncid, varid_lon, 'long_name', 'longitude')
    status = nf90_put_att(ncid, varid_lon, 'standard_name', 'longitude')
    
    status = nf90_def_var(ncid, 'latitude', nf90_double, [dimid_lat], varid_lat)
    status = nf90_put_att(ncid, varid_lat, 'units', 'degrees_north')
    status = nf90_put_att(ncid, varid_lat, 'long_name', 'latitude')
    status = nf90_put_att(ncid, varid_lat, 'standard_name', 'latitude')
    
    ! Create time units string
    write(time_units, '("months since ", I0, "-01-01")') start_year
    status = nf90_def_var(ncid, 'time', nf90_int, [dimid_time], varid_time)
    status = nf90_put_att(ncid, varid_time, 'units', trim(time_units))
    status = nf90_put_att(ncid, varid_time, 'long_name', 'time')
    status = nf90_put_att(ncid, varid_time, 'standard_name', 'time')
    
    ! Define precipitation variable
    status = nf90_def_var(ncid, 'precipitation', nf90_double, &
                         [dimid_lon, dimid_lat, dimid_time], varid_precip)
    status = nf90_put_att(ncid, varid_precip, 'units', 'mm/month')
    status = nf90_put_att(ncid, varid_precip, 'long_name', 'bias-corrected monthly precipitation')
    status = nf90_put_att(ncid, varid_precip, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_precip, 'standard_name', 'precipitation_amount')
    status = nf90_put_att(ncid, varid_precip, 'comment', 'CMIP6 precipitation bias-corrected using ERA5-Land reference')
    
    ! Global attributes
    status = nf90_put_att(ncid, nf90_global, 'title', 'Bias-Corrected CMIP6 Precipitation for Somaliland')
    status = nf90_put_att(ncid, nf90_global, 'institution', 'University of Glasgow')
    status = nf90_put_att(ncid, nf90_global, 'source', 'CMIP6 ' // trim(model) // ' ' // trim(scenario))
    status = nf90_put_att(ncid, nf90_global, 'creator_name', 'Khadar Daahir')
    status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.7')
    status = nf90_put_att(ncid, nf90_global, 'date_created', '2025-08-26')
    status = nf90_put_att(ncid, nf90_global, 'summary', 'Bias-corrected monthly precipitation for Somaliland drought analysis')
    status = nf90_put_att(ncid, nf90_global, 'keywords', 'precipitation, CMIP6, bias correction, Somaliland, drought')
    status = nf90_put_att(ncid, nf90_global, 'scenario', trim(scenario))
    status = nf90_put_att(ncid, nf90_global, 'model', trim(model))
    status = nf90_put_att(ncid, nf90_global, 'time_coverage_start', trim(time_units))
    status = nf90_put_att(ncid, nf90_global, 'references', 'Eyring et al. (2016), doi:10.5194/gmd-9-1937-2016')
    status = nf90_put_att(ncid, nf90_global, 'processing_note', 'Bias correction applied using ERA5-Land 1981-2014 reference period')
    
    ! Add bias correction metadata attributes
    if (present(bias_factor)) then
        status = nf90_put_att(ncid, nf90_global, 'bias_correction_factor', bias_factor)
    end if
    if (present(correlation)) then
        status = nf90_put_att(ncid, nf90_global, 'bias_correction_correlation', correlation)
    end if
    if (present(rmse)) then
        status = nf90_put_att(ncid, nf90_global, 'bias_correction_rmse_mm_month', rmse)
    end if
    status = nf90_put_att(ncid, nf90_global, 'bias_correction_reference', 'ERA5-Land 1981-2014')
    status = nf90_put_att(ncid, nf90_global, 'bias_correction_method', 'Regional averaging with constrained multiplicative factors (0.5-2.0x)')
    
    status = nf90_enddef(ncid)
    
    ! Write data
    status = nf90_put_var(ncid, varid_lon, longitudes)
    status = nf90_put_var(ncid, varid_lat, latitudes)
    status = nf90_put_var(ncid, varid_time, [(i, i=1,ntime)])
    status = nf90_put_var(ncid, varid_precip, precipitation)
    
    status = nf90_close(ncid)
    
end subroutine save_corrected_cmip6_precipitation

!> Export bias-corrected CMIP6 Hargreaves PET
subroutine save_corrected_cmip6_pet_hargreaves(pet, latitudes, longitudes, &
                                             scenario, model, start_year, end_year, &
                                             output_file, status, bias_factor, correlation, rmse)
    real(dp), intent(in) :: pet(:,:,:)
    real(dp), intent(in) :: latitudes(:), longitudes(:)
    character(len=*), intent(in) :: scenario, model
    integer, intent(in) :: start_year, end_year
    character(len=*), intent(in) :: output_file
    integer, intent(out) :: status
    real(dp), intent(in), optional :: bias_factor, correlation, rmse
    
    integer :: ncid, dimid_lon, dimid_lat, dimid_time
    integer :: varid_lon, varid_lat, varid_time, varid_pet
    integer :: nlon, nlat, ntime, i
    character(len=50) :: time_units
    
    nlon = size(pet, 1)
    nlat = size(pet, 2)
    ntime = size(pet, 3)
    
    ! Create NetCDF file
    status = nf90_create(output_file, nf90_clobber, ncid)
    if (status /= nf90_noerr) then
        print *, "ERROR: Cannot create file ", trim(output_file)
        return
    end if
    
    ! Define dimensions
    status = nf90_def_dim(ncid, 'longitude', nlon, dimid_lon)
    status = nf90_def_dim(ncid, 'latitude', nlat, dimid_lat)
    status = nf90_def_dim(ncid, 'time', ntime, dimid_time)
    
    ! Define coordinate variables
    status = nf90_def_var(ncid, 'longitude', nf90_double, [dimid_lon], varid_lon)
    status = nf90_put_att(ncid, varid_lon, 'units', 'degrees_east')
    status = nf90_put_att(ncid, varid_lon, 'long_name', 'longitude')
    status = nf90_put_att(ncid, varid_lon, 'standard_name', 'longitude')
    
    status = nf90_def_var(ncid, 'latitude', nf90_double, [dimid_lat], varid_lat)
    status = nf90_put_att(ncid, varid_lat, 'units', 'degrees_north')
    status = nf90_put_att(ncid, varid_lat, 'long_name', 'latitude')
    status = nf90_put_att(ncid, varid_lat, 'standard_name', 'latitude')
    
    ! Create time units string
    write(time_units, '("months since ", I0, "-01-01")') start_year
    status = nf90_def_var(ncid, 'time', nf90_int, [dimid_time], varid_time)
    status = nf90_put_att(ncid, varid_time, 'units', trim(time_units))
    status = nf90_put_att(ncid, varid_time, 'long_name', 'time')
    status = nf90_put_att(ncid, varid_time, 'standard_name', 'time')
    
    ! Define PET variable
    status = nf90_def_var(ncid, 'potential_evapotranspiration_hargreaves', nf90_double, &
                         [dimid_lon, dimid_lat, dimid_time], varid_pet)
    status = nf90_put_att(ncid, varid_pet, 'units', 'mm/month')
    status = nf90_put_att(ncid, varid_pet, 'long_name', 'bias-corrected monthly potential evapotranspiration (Hargreaves method)')
    status = nf90_put_att(ncid, varid_pet, '_FillValue', FILL_VALUE)
    status = nf90_put_att(ncid, varid_pet, 'standard_name', 'water_potential_evaporation_amount')
    status = nf90_put_att(ncid, varid_pet, 'method', 'Hargreaves equation: PET = 0.0023 * Ra * (Tmean + 17.8) * sqrt(Tmax - Tmin)')
    status = nf90_put_att(ncid, varid_pet, 'comment', 'CMIP6 Hargreaves PET bias-corrected using ERA5-Land reference')
    
    ! Global attributes
    status = nf90_put_att(ncid, nf90_global, 'title', 'Bias-Corrected CMIP6 Hargreaves PET for Somaliland')
    status = nf90_put_att(ncid, nf90_global, 'institution', 'University of Glasgow')
    status = nf90_put_att(ncid, nf90_global, 'source', 'CMIP6 ' // trim(model) // ' ' // trim(scenario))
    status = nf90_put_att(ncid, nf90_global, 'creator_name', 'Khadar Daahir')
    status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.7')
    status = nf90_put_att(ncid, nf90_global, 'date_created', '2025-08-26')
    status = nf90_put_att(ncid, nf90_global, 'summary', 'Bias-corrected monthly Hargreaves PET for Somaliland drought analysis')
    status = nf90_put_att(ncid, nf90_global, 'keywords', 'evapotranspiration, Hargreaves, CMIP6, bias correction, Somaliland, drought')
    status = nf90_put_att(ncid, nf90_global, 'scenario', trim(scenario))
    status = nf90_put_att(ncid, nf90_global, 'model', trim(model))
    status = nf90_put_att(ncid, nf90_global, 'time_coverage_start', trim(time_units))
    status = nf90_put_att(ncid, nf90_global, 'references', 'Hargreaves & Samani (1985); Eyring et al. (2016)')
    status = nf90_put_att(ncid, nf90_global, 'processing_note', 'Bias correction applied using ERA5-Land 1981-2014 reference period')
    
    ! Add bias correction metadata attributes
    if (present(bias_factor)) then
        status = nf90_put_att(ncid, nf90_global, 'bias_correction_factor', bias_factor)
    end if
    if (present(correlation)) then
        status = nf90_put_att(ncid, nf90_global, 'bias_correction_correlation', correlation)
    end if
    if (present(rmse)) then
        status = nf90_put_att(ncid, nf90_global, 'bias_correction_rmse_mm_month', rmse)
    end if
    status = nf90_put_att(ncid, nf90_global, 'bias_correction_reference', 'ERA5-Land 1981-2014')
    status = nf90_put_att(ncid, nf90_global, 'bias_correction_method', 'Regional averaging with constrained multiplicative factors (0.5-2.0x)')
    
    status = nf90_enddef(ncid)
    
    ! Write data
    status = nf90_put_var(ncid, varid_lon, longitudes)
    status = nf90_put_var(ncid, varid_lat, latitudes)
    status = nf90_put_var(ncid, varid_time, [(i, i=1,ntime)])
    status = nf90_put_var(ncid, varid_pet, pet)
    
    status = nf90_close(ncid)
    
end subroutine save_corrected_cmip6_pet_hargreaves

    !-------------------------------------------------------------------
    ! BUILD CMIP6 FILE PATHS
    !-------------------------------------------------------------------
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

    !-------------------------------------------------------------------
    ! READ CMIP6 NETCDF DATA - Full NetCDF implementation
    !-------------------------------------------------------------------
    subroutine read_cmip6_netcdf_data(pr_file, tasmax_file, tasmin_file, cmip6_data, success)
        character(len=*), intent(in) :: pr_file, tasmax_file, tasmin_file
        type(cmip6_scenario_data_t), intent(inout) :: cmip6_data
        logical, intent(out) :: success
        
        integer :: ncid, varid, dimid
        integer :: status
        integer :: nlons, nlats, ntimes
        real(dp), allocatable :: temp_pr(:,:,:)
        
        success = .false.
        
        ! Read precipitation file first to get dimensions
        status = nf90_open(pr_file, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error opening precipitation file: ", trim(pr_file)
            return
        end if
        
        ! Get dimensions
        status = nf90_inq_dimid(ncid, "longitude", dimid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error getting longitude dimension"
            status = nf90_close(ncid)
            return
        end if
        status = nf90_inquire_dimension(ncid, dimid, len=nlons)
        
        status = nf90_inq_dimid(ncid, "latitude", dimid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error getting latitude dimension"
            status = nf90_close(ncid)
            return
        end if
        status = nf90_inquire_dimension(ncid, dimid, len=nlats)
        
        status = nf90_inq_dimid(ncid, "time", dimid)
        status = nf90_inquire_dimension(ncid, dimid, len=ntimes)
        
        ! Store dimensions
        cmip6_data%nlons = nlons
        cmip6_data%nlats = nlats
        cmip6_data%ntimes = ntimes
        
        ! Allocate arrays for raw data
        allocate(cmip6_data%precipitation(nlons, nlats, ntimes))
        allocate(cmip6_data%tasmax(nlons, nlats, ntimes))
        allocate(cmip6_data%tasmin(nlons, nlats, ntimes))
        allocate(cmip6_data%longitude(nlons))
        allocate(cmip6_data%latitude(nlats))
        allocate(cmip6_data%time_values(ntimes))
        allocate(temp_pr(nlons, nlats, ntimes))
        
        ! Allocate arrays for processed data
        allocate(cmip6_data%pet_raw(nlons, nlats, ntimes))
        allocate(cmip6_data%pet_corrected(nlons, nlats, ntimes))
        allocate(cmip6_data%precipitation_corrected(nlons, nlats, ntimes))
        allocate(cmip6_data%land_mask(nlons, nlats))
        allocate(cmip6_data%data_quality(nlons, nlats))
        
        ! Initialize status flags
        cmip6_data%land_mask = .true.  ! Initialize as all land (will be refined later)
        cmip6_data%data_quality = 1.0_dp  ! Initialize as perfect quality
        
        ! Read coordinate variables
        status = nf90_inq_varid(ncid, "longitude", varid)
        status = nf90_get_var(ncid, varid, cmip6_data%longitude)
        
        status = nf90_inq_varid(ncid, "latitude", varid)
        status = nf90_get_var(ncid, varid, cmip6_data%latitude)
        
        status = nf90_inq_varid(ncid, "time", varid)
        status = nf90_get_var(ncid, varid, cmip6_data%time_values)
        
        ! Read precipitation data
        status = nf90_inq_varid(ncid, "pr", varid)
        status = nf90_get_var(ncid, varid, temp_pr)
        
        ! Convert precipitation from kg/m²/s to mm/month
        ! 1 kg/m²/s = 1 mm/s, so multiply by seconds per month (approx 2.628e6)
        cmip6_data%precipitation = temp_pr * 2.628e6_dp
        
        status = nf90_close(ncid)
        
        ! Read tasmax file
        status = nf90_open(tasmax_file, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error opening tasmax file: ", trim(tasmax_file)
            return
        end if
        
        status = nf90_inq_varid(ncid, "tasmax", varid)
        status = nf90_get_var(ncid, varid, cmip6_data%tasmax)
        status = nf90_close(ncid)
        
        ! Read tasmin file
        status = nf90_open(tasmin_file, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error opening tasmin file: ", trim(tasmin_file)
            return
        end if
        
        status = nf90_inq_varid(ncid, "tasmin", varid)
        status = nf90_get_var(ncid, varid, cmip6_data%tasmin)
        status = nf90_close(ncid)
        
        deallocate(temp_pr)
        
        ! Apply quality control and land masking (call bias correction module for now)
        ! TODO: Move quality control to io_module for complete separation
        
        success = .true.
        
        write(*,*) "✅ NetCDF data files read successfully:"
        write(*,*) "   Dimensions: ", nlons, " x ", nlats, " x ", ntimes
        write(*,*) "   Precipitation range: ", minval(cmip6_data%precipitation), " to ", maxval(cmip6_data%precipitation)
        write(*,*) "   Temperature range (K): ", minval(cmip6_data%tasmin), " to ", maxval(cmip6_data%tasmax)
        write(*,*) "   Note: Quality control and PET calculation will be done by bias correction module"
        
    end subroutine read_cmip6_netcdf_data

    !-------------------------------------------------------------------
    ! LOAD RAW CMIP6 SCENARIO DATA
    !-------------------------------------------------------------------
    subroutine load_raw_cmip6_scenario(data_directory, scenario, cmip6_data, status)
        !===========================================================================
        ! Load raw CMIP6 scenario data (pr, tasmax, tasmin) from NetCDF files
        !===========================================================================
        character(len=*), intent(in) :: data_directory
        character(len=*), intent(in) :: scenario
        type(cmip6_scenario_data_t), intent(out) :: cmip6_data
        integer, intent(out) :: status
        
        character(len=512) :: pr_file, tasmax_file, tasmin_file
        logical :: read_success
        
        status = -1  ! Initialize as error
        cmip6_data%scenario = trim(scenario)
        cmip6_data%is_corrected = .false.
        
        ! Set temporal metadata based on scenario
        if (trim(scenario) == "historical") then
            cmip6_data%start_year = 1981
            cmip6_data%end_year = 2014
        else
            cmip6_data%start_year = 2015
            cmip6_data%end_year = 2099
        end if
        
        ! Build file paths based on actual directory structure
        call build_cmip6_file_paths(data_directory, scenario, pr_file, tasmax_file, tasmin_file)
        
        write(*,*) "📂 Loading raw CMIP6 ", trim(scenario), " data..."
        
        ! Read the actual NetCDF files
        call read_cmip6_netcdf_data(pr_file, tasmax_file, tasmin_file, cmip6_data, read_success)
        
        if (read_success) then
            status = 0  ! Success
            write(*,*) "✅ Raw CMIP6 ", trim(scenario), " data loaded successfully"
            write(*,*) "   Dimensions: ", cmip6_data%nlons, " x ", cmip6_data%nlats, " x ", cmip6_data%ntimes
        else
            status = -1  ! Error
            write(*,*) "❌ Failed to load raw CMIP6 ", trim(scenario), " data"
        end if
        
    end subroutine load_raw_cmip6_scenario

!> Save bias-corrected CMIP6 scenario data (precipitation and PET)
!!
!! Complete I/O operation for saving bias-corrected CMIP6 climate data
!! with proper directory structure, metadata, and error handling.
!! 
!! @param[in] cmip6_data Complete CMIP6 scenario data structure
!! @param[in] scenario Scenario name (historical, ssp126, ssp245, ssp585)
!! @param[out] success Operation success status
    subroutine save_corrected_cmip6_scenario(cmip6_data, scenario, success)
        type(cmip6_scenario_data_t), intent(in) :: cmip6_data
        character(len=*), intent(in) :: scenario
        logical, intent(out) :: success
        
        character(len=512) :: output_base_dir, precip_output_dir, pet_output_dir
        character(len=512) :: precip_file, pet_file
        integer :: status
        logical :: precip_success, pet_success
        
        ! TODO: These should be passed from the bias correction module as actual computed metrics
        type(bias_correction_validation_t) :: precip_validation, pet_validation
        
        success = .false.
        
        if (.not. cmip6_data%is_corrected) then
            write(*,*) "❌ Cannot save: data not bias-corrected yet"
            return
        end if
        
        ! Initialize validation metrics (placeholder - should come from bias correction module)
        precip_validation%bias = 0.0_dp
        precip_validation%rmse = 0.0_dp
        precip_validation%correlation = 0.0_dp
        precip_validation%ks_statistic = 0.0_dp
        precip_validation%mae = 0.0_dp
        precip_validation%validation_passed = .true.
        
        pet_validation%bias = 0.0_dp
        pet_validation%rmse = 0.0_dp
        pet_validation%correlation = 0.0_dp
        pet_validation%ks_statistic = 0.0_dp
        pet_validation%mae = 0.0_dp
        pet_validation%validation_passed = .true.
        
        ! Create output directory paths following established structure
        output_base_dir = "data/Processed_data/Climate_Drivers/Corrected_CMIP6/cmip6_" // trim(scenario)
        precip_output_dir = trim(output_base_dir) // "/precipitation"
        pet_output_dir = trim(output_base_dir) // "/PET"
        
        ! Create output file paths
        precip_file = trim(precip_output_dir) // "/corrected_precipitation_" // trim(scenario) // "_1981_2014.nc"
        pet_file = trim(pet_output_dir) // "/corrected_pet_hargreaves_" // trim(scenario) // "_1981_2014.nc"
        
        ! Adjust file names for future scenarios
        if (trim(scenario) /= "historical") then
            precip_file = trim(precip_output_dir) // "/corrected_precipitation_" // trim(scenario) // "_2015_2099.nc"
            pet_file = trim(pet_output_dir) // "/corrected_pet_hargreaves_" // trim(scenario) // "_2015_2099.nc"
        end if
        
        write(*,*) "💾 Saving corrected ", trim(scenario), " data with comprehensive metadata..."
        write(*,*) "   Precipitation: ", trim(precip_file)
        write(*,*) "   PET: ", trim(pet_file)
        
        ! Create directories
        call create_output_directory(precip_output_dir)
        call create_output_directory(pet_output_dir)
        
        ! Save using enhanced functions with bias correction metadata
        call save_corrected_cmip6_precipitation_enhanced(cmip6_data%precipitation_corrected, &
                                                        cmip6_data%latitude, cmip6_data%longitude, &
                                                        trim(scenario), "MPI-ESM1-2-LR", &
                                                        cmip6_data%start_year, cmip6_data%end_year, &
                                                        precip_file, status, precip_validation)
        precip_success = (status == 0)
        
        call save_corrected_cmip6_pet_enhanced(cmip6_data%pet_corrected, &
                                              cmip6_data%latitude, cmip6_data%longitude, &
                                              trim(scenario), "MPI-ESM1-2-LR", &
                                              cmip6_data%start_year, cmip6_data%end_year, &
                                              pet_file, status, pet_validation)
        pet_success = (status == 0)
        
        if (precip_success .and. pet_success) then
            write(*,*) "✅ Corrected ", trim(scenario), " data saved with comprehensive metadata"
            success = .true.
        else
            write(*,*) "❌ Failed to save corrected ", trim(scenario), " data"
            if (.not. precip_success) write(*,*) "   Precipitation save failed"
            if (.not. pet_success) write(*,*) "   PET save failed"
        end if
        
    end subroutine save_corrected_cmip6_scenario

!> Load bias-corrected CMIP6 scenario data for projections
!!
!! Loads the bias-corrected precipitation and PET data that was saved by
!! the bias correction module for use in future SPEI/SPI calculations.
!!
!! @param[in] scenario Scenario name (ssp126, ssp245, ssp585)
!! @param[out] precipitation Corrected precipitation data [mm/month]
!! @param[out] pet Corrected PET data [mm/month]  
!! @param[out] latitudes Latitude coordinates [degrees]
!! @param[out] longitudes Longitude coordinates [degrees]
!! @param[out] time_values Time coordinate values
!! @param[out] start_year Starting year of data
!! @param[out] end_year Ending year of data
!! @param[out] status Operation status (0=success, /=0=error)
    subroutine load_corrected_cmip6_scenario(scenario, precipitation, pet, latitudes, longitudes, &
                                            time_values, start_year, end_year, status)
        character(len=*), intent(in) :: scenario
        real(dp), allocatable, intent(out) :: precipitation(:,:,:)
        real(dp), allocatable, intent(out) :: pet(:,:,:)
        real(dp), allocatable, intent(out) :: latitudes(:)
        real(dp), allocatable, intent(out) :: longitudes(:)
        real(dp), allocatable, intent(out) :: time_values(:)
        integer, intent(out) :: start_year, end_year
        integer, intent(out) :: status
        
        character(len=512) :: precip_file, pet_file, base_dir
        integer :: ncid_p, ncid_pet, varid, dimid
        integer :: nlons, nlats, ntimes
        
        status = -1
        
        ! Build file paths for corrected data
        base_dir = "data/Processed_data/Climate_Drivers/Corrected_CMIP6/cmip6_" // trim(scenario)
        
        ! Set temporal metadata and file paths based on scenario
        if (trim(scenario) == "historical") then
            start_year = 1981
            end_year = 2014
            precip_file = trim(base_dir) // "/precipitation/corrected_precipitation_" // trim(scenario) // "_1981_2014.nc"
            pet_file = trim(base_dir) // "/PET/corrected_pet_hargreaves_" // trim(scenario) // "_1981_2014.nc"
        else
            start_year = 2015
            end_year = 2099
            precip_file = trim(base_dir) // "/precipitation/corrected_precipitation_" // trim(scenario) // "_2015_2099.nc"
            pet_file = trim(base_dir) // "/PET/corrected_pet_hargreaves_" // trim(scenario) // "_2015_2099.nc"
        end if
        
        write(*,*) "📂 Loading corrected CMIP6 ", trim(scenario), " data..."
        write(*,*) "   Precipitation: ", trim(precip_file)
        write(*,*) "   PET: ", trim(pet_file)
        
        ! Open precipitation file and get dimensions
        if (nf90_open(precip_file, NF90_NOWRITE, ncid_p) /= NF90_NOERR) then
            write(*,*) "❌ Error opening precipitation file: ", trim(precip_file)
            return
        end if
        
        ! Get dimensions from precipitation file
        if (nf90_inq_dimid(ncid_p, "longitude", dimid) /= NF90_NOERR) then
            write(*,*) "❌ Error getting longitude dimension"
            status = nf90_close(ncid_p)
            return
        end if
        status = nf90_inquire_dimension(ncid_p, dimid, len=nlons)
        
        if (nf90_inq_dimid(ncid_p, "latitude", dimid) /= NF90_NOERR) then
            write(*,*) "❌ Error getting latitude dimension"
            status = nf90_close(ncid_p)
            return
        end if
        status = nf90_inquire_dimension(ncid_p, dimid, len=nlats)
        
        if (nf90_inq_dimid(ncid_p, "time", dimid) /= NF90_NOERR) then
            write(*,*) "❌ Error getting time dimension"
            status = nf90_close(ncid_p)
            return
        end if
        status = nf90_inquire_dimension(ncid_p, dimid, len=ntimes)
        
        ! Allocate arrays
        allocate(precipitation(nlons, nlats, ntimes))
        allocate(pet(nlons, nlats, ntimes))
        allocate(latitudes(nlats))
        allocate(longitudes(nlons))
        allocate(time_values(ntimes))
        
        ! Read coordinate variables
        if (nf90_inq_varid(ncid_p, "longitude", varid) /= NF90_NOERR) then
            write(*,*) "❌ Error finding longitude variable"
            status = nf90_close(ncid_p)
            return
        end if
        status = nf90_get_var(ncid_p, varid, longitudes)
        
        if (nf90_inq_varid(ncid_p, "latitude", varid) /= NF90_NOERR) then
            write(*,*) "❌ Error finding latitude variable"
            status = nf90_close(ncid_p)
            return
        end if
        status = nf90_get_var(ncid_p, varid, latitudes)
        
        if (nf90_inq_varid(ncid_p, "time", varid) /= NF90_NOERR) then
            write(*,*) "❌ Error finding time variable"
            status = nf90_close(ncid_p)
            return
        end if
        status = nf90_get_var(ncid_p, varid, time_values)
        
        ! Read precipitation data
        if (nf90_inq_varid(ncid_p, "precipitation", varid) /= NF90_NOERR) then
            write(*,*) "❌ Error finding precipitation variable"
            status = nf90_close(ncid_p)
            return
        end if
        status = nf90_get_var(ncid_p, varid, precipitation)
        status = nf90_close(ncid_p)
        
        ! Open and read PET file
        if (nf90_open(pet_file, NF90_NOWRITE, ncid_pet) /= NF90_NOERR) then
            write(*,*) "❌ Error opening PET file: ", trim(pet_file)
            return
        end if
        
        if (nf90_inq_varid(ncid_pet, "potential_evapotranspiration_hargreaves", varid) /= NF90_NOERR) then
            write(*,*) "❌ Error finding PET variable 'potential_evapotranspiration_hargreaves' in file: ", trim(pet_file)
            status = nf90_close(ncid_pet)
            return
        else
            write(*,*) "✅ Found PET variable 'potential_evapotranspiration_hargreaves'"
        end if
        status = nf90_get_var(ncid_pet, varid, pet)
        status = nf90_close(ncid_pet)
        
        status = 0
        write(*,*) "✅ Corrected CMIP6 data loaded successfully:"
        write(*,*) "   Dimensions: ", nlons, " x ", nlats, " x ", ntimes
        write(*,*) "   Period: ", start_year, "-", end_year
        write(*,*) "   Precipitation range: ", minval(precipitation), " to ", maxval(precipitation), " mm/month"
        write(*,*) "   PET range: ", minval(pet), " to ", maxval(pet), " mm/month"
        
    end subroutine load_corrected_cmip6_scenario

!> Save future drought indices for CMIP6 scenarios
!!
!! Creates the specific directory structure for CMIP6 drought indices:
!! data/Final_results/Drought_Indices/CMIP6/{scenario}/
!! Saves both multiscale and individual timescale files
!!
!! @param[in] scenario_name SSP scenario (ssp126, ssp245, ssp585)
!! @param[in] spi_* SPI arrays for different timescales
!! @param[in] spei_* SPEI arrays for different timescales
!! @param[in] latitudes Latitude coordinate array
!! @param[in] longitudes Longitude coordinate array
!! @param[in] time_values Time coordinate array
!! @param[in] start_year Start year of the data
!! @param[in] end_year End year of the data
!! @param[out] status Operation status (0=success, 1=error)
    subroutine save_future_drought_indices(scenario_name, spi_1, spi_3, spi_6, spi_12, &
                                          spei_1, spei_3, spei_6, spei_12, &
                                          latitudes, longitudes, time_values, &
                                          start_year, end_year, status)
        character(len=*), intent(in) :: scenario_name
        real(dp), intent(in) :: spi_1(:,:,:), spi_3(:,:,:), spi_6(:,:,:), spi_12(:,:,:)
        real(dp), intent(in) :: spei_1(:,:,:), spei_3(:,:,:), spei_6(:,:,:), spei_12(:,:,:)
        real(dp), intent(in) :: latitudes(:), longitudes(:), time_values(:)
        integer, intent(in) :: start_year, end_year
        integer, intent(out) :: status
        
        character(len=512) :: output_dir, file_path
        character(len=50) :: period_str
        
        status = 1  ! Default to error
        
        ! Create period string
        write(period_str, '(I0,A,I0)') start_year, "_", end_year
        
        ! Create base output directory for this scenario (following documented standards)
        output_dir = "data/Processed_data/Drought_indices/CMIP6/" // trim(scenario_name)
        
        write(*,*) "💾 Saving future drought indices for ", trim(scenario_name), "..."
        write(*,*) "   Output directory: ", trim(output_dir)
        write(*,*) "   Period: ", trim(period_str)
        
        ! Create directory structure
        call create_output_directory(output_dir)
        
        ! Save SPEI multiscale (all timescales in one file)
        file_path = trim(output_dir) // "/spei_multiscale_" // trim(scenario_name) // "_" // trim(period_str) // ".nc"
        call save_spei_multiscale(spei_1, spei_3, spei_6, spei_12, latitudes, longitudes, &
                                 file_path, status)
        if (status /= 0) then
            write(*,*) "❌ Failed to save SPEI multiscale for ", trim(scenario_name)
            return
        end if
        write(*,*) "   ✅ SPEI multiscale saved: ", trim(file_path)
        
        ! Save individual SPEI timescales
        ! SPEI-1
        file_path = trim(output_dir) // "/spei_1_" // trim(scenario_name) // "_" // trim(period_str) // ".nc"
        call save_spei_single_scale(spei_1, 1, latitudes, longitudes, file_path, status)
        if (status /= 0) then
            write(*,*) "❌ Failed to save SPEI-1 for ", trim(scenario_name)
            return
        end if
        
        ! SPEI-3
        file_path = trim(output_dir) // "/spei_3_" // trim(scenario_name) // "_" // trim(period_str) // ".nc"
        call save_spei_single_scale(spei_3, 3, latitudes, longitudes, file_path, status)
        if (status /= 0) then
            write(*,*) "❌ Failed to save SPEI-3 for ", trim(scenario_name)
            return
        end if
        
        ! SPEI-6
        file_path = trim(output_dir) // "/spei_6_" // trim(scenario_name) // "_" // trim(period_str) // ".nc"
        call save_spei_single_scale(spei_6, 6, latitudes, longitudes, file_path, status)
        if (status /= 0) then
            write(*,*) "❌ Failed to save SPEI-6 for ", trim(scenario_name)
            return
        end if
        
        ! SPEI-12
        file_path = trim(output_dir) // "/spei_12_" // trim(scenario_name) // "_" // trim(period_str) // ".nc"
        call save_spei_single_scale(spei_12, 12, latitudes, longitudes, file_path, status)
        if (status /= 0) then
            write(*,*) "❌ Failed to save SPEI-12 for ", trim(scenario_name)
            return
        end if
        
      
        status = 0  ! Success
        
        write(*,*) "✅ All future drought indices saved successfully for ", trim(scenario_name)
        write(*,*) "   Directory structure created:"
        write(*,*) "   ", trim(output_dir), "/"
        write(*,*) "   ├── spei_multiscale_", trim(scenario_name), "_", trim(period_str), ".nc"
        write(*,*) "   ├── spei_1_", trim(scenario_name), "_", trim(period_str), ".nc"
        write(*,*) "   ├── spei_3_", trim(scenario_name), "_", trim(period_str), ".nc"
        write(*,*) "   ├── spei_6_", trim(scenario_name), "_", trim(period_str), ".nc"
        write(*,*) "   └── spei_12_", trim(scenario_name), "_", trim(period_str), ".nc"
        
    end subroutine save_future_drought_indices

    !===============================================================================
    ! SAVE EVT RESULTS TO NETCDF
    !===============================================================================
    !> Save Extreme Value Theory (EVT) analysis results to NetCDF
    !!
    !! Saves GPD parameters, return levels, and exceedance probabilities from
    !! extreme drought analysis using FSML-based Generalized Pareto Distribution fitting.
    !!
    !! @param[in] filename Output NetCDF file path
    !! @param[in] gpd_params GPD parameters [location(μ), scale(σ), shape(ξ)]
    !! @param[in] return_levels Return period levels [10yr, 20yr, 50yr]
    !! @param[in] exceed_prob Probability of exceeding historical 20-year event
    !! @param[in] scenario Climate scenario name (era5, ssp126, ssp245, ssp585)
    !! @param[in] timescale SPEI timescale used (1, 3, 6, 12 months)
    !! @param[in] threshold GPD threshold value used for exceedances
    !! @param[out] status Success/failure status (0=success)
    subroutine save_evt_results_nc(filename, gpd_params, return_levels, exceed_prob, &
                                   scenario, timescale, threshold, status)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none
        
        character(len=*), intent(in) :: filename, scenario
        integer, intent(in) :: timescale
        real(real64), intent(in) :: gpd_params(3), return_levels(3), exceed_prob, threshold
        integer, intent(out) :: status
        
        integer :: ncid, dimid_scalar, dimid_params, dimid_returns
        integer :: varid_gpd, varid_returns, varid_exceed, varid_threshold
        character(len=10) :: timescale_str
        character(len=25) :: timestamp
        
        ! Convert timescale to string
        write(timescale_str, '(I0)') timescale
        
        ! Get current timestamp
        call get_timestamp(timestamp)
        
        ! Create NetCDF file
        status = nf90_create(filename, NF90_CLOBBER, ncid)
        if (status /= NF90_NOERR) then
            write(*,*) "ERROR: Failed to create EVT results file: ", trim(filename)
            return
        end if
        
        ! Define dimensions
        status = nf90_def_dim(ncid, "scalar", 1, dimid_scalar)
        status = nf90_def_dim(ncid, "gpd_parameters", 3, dimid_params)      ! [μ, σ, ξ]
        status = nf90_def_dim(ncid, "return_periods", 3, dimid_returns)     ! [10, 20, 50 years]
        
        ! Define variables
        status = nf90_def_var(ncid, "gpd_parameters", NF90_DOUBLE, [dimid_params], varid_gpd)
        status = nf90_def_var(ncid, "return_levels", NF90_DOUBLE, [dimid_returns], varid_returns)
        status = nf90_def_var(ncid, "exceed_probability", NF90_DOUBLE, [dimid_scalar], varid_exceed)
        status = nf90_def_var(ncid, "threshold", NF90_DOUBLE, [dimid_scalar], varid_threshold)
        
        ! Add variable attributes
        status = nf90_put_att(ncid, varid_gpd, "long_name", "Generalized Pareto Distribution Parameters")
        status = nf90_put_att(ncid, varid_gpd, "units", "1")
        status = nf90_put_att(ncid, varid_gpd, "description", "Location (mu), Scale (sigma), Shape (xi) parameters")
        status = nf90_put_att(ncid, varid_gpd, "_FillValue", FILL_VALUE)
        
        status = nf90_put_att(ncid, varid_returns, "long_name", "Return Period Levels")
        status = nf90_put_att(ncid, varid_returns, "units", "SPEI units")
        status = nf90_put_att(ncid, varid_returns, "description", "10-year, 20-year, 50-year return levels")
        status = nf90_put_att(ncid, varid_returns, "_FillValue", FILL_VALUE)
        
        status = nf90_put_att(ncid, varid_exceed, "long_name", "Exceedance Probability")
        status = nf90_put_att(ncid, varid_exceed, "units", "1")
        status = nf90_put_att(ncid, varid_exceed, "description", "Probability of exceeding historical 20-year event")
        
        status = nf90_put_att(ncid, varid_threshold, "long_name", "GPD Threshold")
        status = nf90_put_att(ncid, varid_threshold, "units", "SPEI units")
        status = nf90_put_att(ncid, varid_threshold, "description", "Threshold for peak-over-threshold analysis")
        
        ! Add global attributes
        status = nf90_put_att(ncid, NF90_GLOBAL, "title", &
                             "Extreme Value Theory Analysis for Drought Events")
        status = nf90_put_att(ncid, NF90_GLOBAL, "institution", "University of Glasgow")
        status = nf90_put_att(ncid, NF90_GLOBAL, "source", "Somaliland Drought Analysis Pipeline")
        status = nf90_put_att(ncid, NF90_GLOBAL, "method", "FSML-GPD Peak-over-Threshold Analysis")
        status = nf90_put_att(ncid, NF90_GLOBAL, "scenario", trim(scenario))
        status = nf90_put_att(ncid, NF90_GLOBAL, "spei_timescale", trim(timescale_str) // " months")
        status = nf90_put_att(ncid, NF90_GLOBAL, "created_on", trim(timestamp))
        status = nf90_put_att(ncid, NF90_GLOBAL, "contact", "khadardahir146@gmail.com")
        status = nf90_put_att(ncid, NF90_GLOBAL, "references", &
                             "Coles (2001) Introduction to Statistical Modeling of Extreme Values")
        status = nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8")
        
        ! End definitions
        status = nf90_enddef(ncid)
        if (status /= NF90_NOERR) then
            write(*,*) "ERROR: Failed to end definitions in EVT file"
            return
        end if
        
        ! Write data
        status = nf90_put_var(ncid, varid_gpd, gpd_params)
        status = nf90_put_var(ncid, varid_returns, return_levels)
        status = nf90_put_var(ncid, varid_exceed, exceed_prob)
        status = nf90_put_var(ncid, varid_threshold, threshold)
        
        ! Close file
        status = nf90_close(ncid)
        
        if (status == NF90_NOERR) then
            write(*,*) "✅ EVT results saved to ", trim(filename)
            status = 0  ! Success
        else
            write(*,*) "❌ Failed to write EVT results to ", trim(filename)
            status = 1  ! Failure
        end if
        
    end subroutine save_evt_results_nc
    
    !-------------------------------------------------------------------
    ! SAVE ENHANCED EVT RESULTS with Profile MLE output
    ! Following academic standards for EVT analysis reporting
    !-------------------------------------------------------------------
    subroutine save_evt_results_enhanced(filename, results, scenario, timescale, status)
        use evt_module, only: evt_results
        use netcdf
        character(len=*), intent(in) :: filename, scenario
        type(evt_results), intent(in) :: results
        integer, intent(in) :: timescale
        integer, intent(out) :: status
        
        integer :: ncid, varid_xi, varid_sigma, varid_mu, varid_loglik
        integer :: varid_xi_se, varid_sigma_se, varid_rl10, varid_rl20, varid_rl50
        integer :: varid_ci_lower, varid_ci_upper, varid_exceed_prob, varid_threshold
        integer :: varid_n_exceed, varid_mean_excess, dim_return_period
        character(len=25) :: timestamp
        
        status = 0
        
        ! Extract directory from filename and ensure it exists
        call create_output_directory(filename(:index(filename, '/', back=.true.)-1))
        
        ! Create NetCDF file
        if (nf90_create(filename, NF90_CLOBBER, ncid) == NF90_NOERR) then
            
            ! Define dimensions
            call check_nc(nf90_def_dim(ncid, "return_period", 3, dim_return_period))
            
            ! Define variables with enhanced metadata
            call check_nc(nf90_def_var(ncid, "gpd_shape", NF90_DOUBLE, varid_xi))
            call check_nc(nf90_put_att(ncid, varid_xi, "long_name", &
                "GPD shape parameter (xi) from Profile Maximum Likelihood"))
            call check_nc(nf90_put_att(ncid, varid_xi, "units", "dimensionless"))
            call check_nc(nf90_put_att(ncid, varid_xi, "references", &
                "Coles (2001), Grimshaw (1993)"))
            
            call check_nc(nf90_def_var(ncid, "gpd_scale", NF90_DOUBLE, varid_sigma))
            call check_nc(nf90_put_att(ncid, varid_sigma, "long_name", &
                "GPD scale parameter (sigma) from Profile MLE"))
            call check_nc(nf90_put_att(ncid, varid_sigma, "units", "SPEI units"))
            
            call check_nc(nf90_def_var(ncid, "gpd_location", NF90_DOUBLE, varid_mu))
            call check_nc(nf90_put_att(ncid, varid_mu, "long_name", &
                "GPD location parameter (threshold)"))
            call check_nc(nf90_put_att(ncid, varid_mu, "units", "SPEI units"))
            
            call check_nc(nf90_def_var(ncid, "log_likelihood", NF90_DOUBLE, varid_loglik))
            call check_nc(nf90_put_att(ncid, varid_loglik, "long_name", &
                "Profile log-likelihood value"))
            
            call check_nc(nf90_def_var(ncid, "xi_standard_error", NF90_DOUBLE, varid_xi_se))
            call check_nc(nf90_put_att(ncid, varid_xi_se, "long_name", &
                "Standard error of shape parameter"))
            call check_nc(nf90_put_att(ncid, varid_xi_se, "method", &
                "Fisher Information Matrix"))
            
            call check_nc(nf90_def_var(ncid, "sigma_standard_error", NF90_DOUBLE, varid_sigma_se))
            call check_nc(nf90_put_att(ncid, varid_sigma_se, "long_name", &
                "Standard error of scale parameter"))
            
            ! Return levels
            call check_nc(nf90_def_var(ncid, "return_levels", NF90_DOUBLE, &
                [dim_return_period], varid_rl10))
            call check_nc(nf90_put_att(ncid, varid_rl10, "long_name", &
                "Return levels for 10, 20, 50-year return periods"))
            call check_nc(nf90_put_att(ncid, varid_rl10, "units", "SPEI units"))
            
            ! Confidence intervals
            call check_nc(nf90_def_var(ncid, "return_levels_ci_lower", NF90_DOUBLE, &
                [dim_return_period], varid_ci_lower))
            call check_nc(nf90_put_att(ncid, varid_ci_lower, "long_name", &
                "95% confidence interval lower bounds"))
            
            call check_nc(nf90_def_var(ncid, "return_levels_ci_upper", NF90_DOUBLE, &
                [dim_return_period], varid_ci_upper))
            call check_nc(nf90_put_att(ncid, varid_ci_upper, "long_name", &
                "95% confidence interval upper bounds"))
            
            ! Additional diagnostics
            call check_nc(nf90_def_var(ncid, "exceedance_probability", NF90_DOUBLE, varid_exceed_prob))
            call check_nc(nf90_put_att(ncid, varid_exceed_prob, "long_name", &
                "Probability of exceeding 20-year return level"))
            
            call check_nc(nf90_def_var(ncid, "threshold", NF90_DOUBLE, varid_threshold))
            call check_nc(nf90_put_att(ncid, varid_threshold, "long_name", &
                "Peak-over-threshold value"))
            
            call check_nc(nf90_def_var(ncid, "n_exceedances", NF90_INT, varid_n_exceed))
            call check_nc(nf90_put_att(ncid, varid_n_exceed, "long_name", &
                "Number of threshold exceedances"))
            
            call check_nc(nf90_def_var(ncid, "mean_excess", NF90_DOUBLE, varid_mean_excess))
            call check_nc(nf90_put_att(ncid, varid_mean_excess, "long_name", &
                "Mean excess above threshold"))
            
            ! Enhanced global attributes using centralized metadata function
            call add_enhanced_metadata(ncid, "Extreme Value Theory Analysis", scenario, 1981, 2100)
            
            ! EVT-specific metadata
            call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "spei_timescale", timescale))
            call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "methodology", &
                "Profile Maximum Likelihood Estimation with Golden Section Optimization"))
            call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "evt_references", &
                "Coles (2001), Grimshaw (1993), Hosking & Wallis (1987)"))
            
            ! End definition mode
            call check_nc(nf90_enddef(ncid))
            
            ! Write data
            call check_nc(nf90_put_var(ncid, varid_xi, results%params%xi))
            call check_nc(nf90_put_var(ncid, varid_sigma, results%params%sigma))
            call check_nc(nf90_put_var(ncid, varid_mu, results%params%mu))
            call check_nc(nf90_put_var(ncid, varid_loglik, results%params%loglik))
            call check_nc(nf90_put_var(ncid, varid_xi_se, results%params%xi_se))
            call check_nc(nf90_put_var(ncid, varid_sigma_se, results%params%sigma_se))
            call check_nc(nf90_put_var(ncid, varid_rl10, results%return_levels))
            call check_nc(nf90_put_var(ncid, varid_ci_lower, results%return_levels_ci(:,1)))
            call check_nc(nf90_put_var(ncid, varid_ci_upper, results%return_levels_ci(:,2)))
            call check_nc(nf90_put_var(ncid, varid_exceed_prob, results%exceed_prob))
            call check_nc(nf90_put_var(ncid, varid_threshold, results%threshold))
            call check_nc(nf90_put_var(ncid, varid_n_exceed, results%n_exceedances))
            call check_nc(nf90_put_var(ncid, varid_mean_excess, results%mean_excess_function))
            
            ! Close file
            call check_nc(nf90_close(ncid))
            
            write(*,*) "✅ Enhanced EVT results saved to ", trim(filename)
            
        else
            write(*,*) "❌ Failed to create enhanced EVT file ", trim(filename)
            status = 1
        end if
        
    end subroutine save_evt_results_enhanced

    !> Check NetCDF operation status and handle errors
    !! 
    !! Utility subroutine for robust NetCDF error handling throughout
    !! the I/O operations, ensuring data integrity and proper debugging.
    !!
    !! @param[in] status NetCDF operation return status
    subroutine check_nc(status)
        integer, intent(in) :: status
        
        if (status /= nf90_noerr) then
            print *, "NetCDF Error: ", trim(nf90_strerror(status))
            stop "NetCDF operation failed"
        end if
    end subroutine check_nc

    !> Get current timestamp for metadata
    subroutine get_timestamp(timestamp)
        character(len=25), intent(out) :: timestamp
        character(len=8) :: date
        character(len=10) :: time
        call date_and_time(date, time)
        timestamp = date(1:4) // "-" // date(5:6) // "-" // date(7:8) // "T" // &
                   time(1:2) // ":" // time(3:4) // ":" // time(5:6)
    end subroutine get_timestamp

    !> Load NetCDF SPEI data for EVT analysis
    !! 
    !! Load SPEI indices from NetCDF files for extreme value analysis.
    !! This subroutine handles both ERA5 and CMIP6 SPEI data files.
    !!
    !! @param[in] filename Path to NetCDF file containing SPEI data
    !! @param[out] spei_data SPEI values (lon, lat, time)
    !! @param[out] lat Latitude values
    !! @param[out] lon Longitude values  
    !! @param[out] time Time values
    !! @param[out] status Operation status (0=success)
    subroutine load_netcdf_spei_data(filename, spei_data, lat, lon, time, status)
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: spei_data(:,:,:)
        real(dp), allocatable, intent(out) :: lat(:), lon(:), time(:)
        integer, intent(out) :: status
        
        integer :: ncid, varid_spei, varid_lat, varid_lon, varid_time
        integer :: nlon, nlat, ntime
        integer :: dimid_lon, dimid_lat, dimid_time
        integer :: nc_status, i, j, k
        character(len=32) :: spei_var_name
        real(dp), allocatable :: temp_spei(:,:,:)  ! Temporary array with file dimensions
        
        status = 0
        
        ! Open NetCDF file
        nc_status = nf90_open(filename, nf90_nowrite, ncid)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot open file: ", trim(filename)
            status = 1
            return
        end if
        
        ! Get dimension IDs and sizes
        nc_status = nf90_inq_dimid(ncid, "longitude", dimid_lon)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot find longitude dimension in file"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_inq_dimid(ncid, "latitude", dimid_lat)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot find latitude dimension in file"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_inq_dimid(ncid, "time", dimid_time)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot find time dimension in file"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_inquire_dimension(ncid, dimid_lon, len=nlon)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot get longitude dimension size"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_inquire_dimension(ncid, dimid_lat, len=nlat)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot get latitude dimension size"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_inquire_dimension(ncid, dimid_time, len=ntime)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot get time dimension size"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        ! Allocate arrays
        allocate(spei_data(nlon, nlat, ntime))
        allocate(lat(nlat))
        allocate(lon(nlon))
        allocate(time(ntime))
        
        ! For SPEI data, we need a temporary array with NetCDF Fortran dimensions
        ! File shows: spei_3(time, latitude, longitude) 
        ! Fortran reads: temp_spei(longitude, latitude, time)
        allocate(temp_spei(nlon, nlat, ntime))
        
        ! Determine SPEI variable name from filename
        if (index(filename, "spei_3") > 0) then
            spei_var_name = "spei_3"
        else if (index(filename, "spei_6") > 0) then
            spei_var_name = "spei_6"
        else if (index(filename, "spei_12") > 0) then
            spei_var_name = "spei_12"
        else
            spei_var_name = "spei"  ! Default fallback
        end if
        
        ! Get variable IDs
        nc_status = nf90_inq_varid(ncid, trim(spei_var_name), varid_spei)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot find ", trim(spei_var_name), " variable in file"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_inq_varid(ncid, "latitude", varid_lat)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot find latitude variable in file"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_inq_varid(ncid, "longitude", varid_lon)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot find longitude variable in file"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_inq_varid(ncid, "time", varid_time)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot find time variable in file"
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        ! Read data with proper dimensions - NetCDF Fortran reads in Fortran order 
        ! File: spei_3(time, latitude, longitude) -> Fortran: temp_spei(longitude, latitude, time)
        nc_status = nf90_get_var(ncid, varid_spei, temp_spei)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot read SPEI data from file: ", nf90_strerror(nc_status)
            print *, "   Expected dimensions: ", ntime, "×", nlat, "×", nlon
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_get_var(ncid, varid_lat, lat)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot read latitude data from file: ", nf90_strerror(nc_status)
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_get_var(ncid, varid_lon, lon)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot read longitude data from file: ", nf90_strerror(nc_status)
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        nc_status = nf90_get_var(ncid, varid_time, time)
        if (nc_status /= nf90_noerr) then
            print *, "❌ Error: Cannot read time data from file: ", nf90_strerror(nc_status)
            nc_status = nf90_close(ncid)
            status = 1
            return
        end if
        
        ! Copy data directly since dimensions are already correct (lon, lat, time)
        spei_data = temp_spei
        
        ! Clean up
        deallocate(temp_spei)
        
        ! Close file
        nc_status = nf90_close(ncid)
        if (nc_status /= nf90_noerr) then
            print *, "⚠️  Warning: Problem closing NetCDF file"
        end if
        
        print *, "✅ SPEI data loaded successfully: ", trim(spei_var_name), " ", nlon, "×", nlat, "×", ntime
        
    end subroutine load_netcdf_spei_data

    !> Enhanced save function for bias-corrected CMIP6 precipitation with comprehensive metadata
    subroutine save_corrected_cmip6_precipitation_enhanced(precipitation, lats, lons, scenario, &
                                                          model, start_year, end_year, filename, &
                                                          status, validation_metrics)
        real(dp), intent(in) :: precipitation(:,:,:)
        real(dp), intent(in) :: lats(:), lons(:)
        character(len=*), intent(in) :: scenario, model, filename
        integer, intent(in) :: start_year, end_year
        integer, intent(out) :: status
        type(bias_correction_validation_t), intent(in) :: validation_metrics
        
        integer :: ncid, lon_dimid, lat_dimid, time_dimid
        integer :: lon_varid, lat_varid, time_varid, precip_varid
        integer :: nlons, nlats, ntimes
        
        nlons = size(lons)
        nlats = size(lats)
        ntimes = size(precipitation, 3)
        
        ! Create NetCDF file
        status = nf90_create(filename, NF90_CLOBBER, ncid)
        if (status /= nf90_noerr) return
        
        ! Define dimensions
        status = nf90_def_dim(ncid, "longitude", nlons, lon_dimid)
        if (status /= nf90_noerr) return
        status = nf90_def_dim(ncid, "latitude", nlats, lat_dimid)
        if (status /= nf90_noerr) return
        status = nf90_def_dim(ncid, "time", ntimes, time_dimid)
        if (status /= nf90_noerr) return
        
        ! Define coordinate variables
        status = nf90_def_var(ncid, "longitude", NF90_DOUBLE, lon_dimid, lon_varid)
        if (status /= nf90_noerr) return
        status = nf90_def_var(ncid, "latitude", NF90_DOUBLE, lat_dimid, lat_varid)
        if (status /= nf90_noerr) return
        status = nf90_def_var(ncid, "time", NF90_DOUBLE, time_dimid, time_varid)
        if (status /= nf90_noerr) return
        
        ! Define data variable
        status = nf90_def_var(ncid, "precipitation", NF90_DOUBLE, &
                             [lon_dimid, lat_dimid, time_dimid], precip_varid)
        if (status /= nf90_noerr) return
        
        ! Add variable attributes
        call check_nc(nf90_put_att(ncid, precip_varid, "units", "mm/month"))
        call check_nc(nf90_put_att(ncid, precip_varid, "long_name", &
            "Bias-corrected monthly precipitation"))
        call check_nc(nf90_put_att(ncid, precip_varid, "standard_name", "precipitation_amount"))
        call check_nc(nf90_put_att(ncid, precip_varid, "_FillValue", -999.0_dp))
        call check_nc(nf90_put_att(ncid, precip_varid, "missing_value", -999.0_dp))
        
        ! Add coordinate attributes
        call check_nc(nf90_put_att(ncid, lon_varid, "units", "degrees_east"))
        call check_nc(nf90_put_att(ncid, lat_varid, "units", "degrees_north"))
        call check_nc(nf90_put_att(ncid, time_varid, "units", "months since 1981-01-01"))
        
        ! Add comprehensive bias correction metadata
        call add_bias_correction_metadata(ncid, scenario, "precipitation", validation_metrics)
        
        ! End define mode
        status = nf90_enddef(ncid)
        if (status /= nf90_noerr) return
        
        ! Write data
        status = nf90_put_var(ncid, lon_varid, lons)
        if (status /= nf90_noerr) return
        status = nf90_put_var(ncid, lat_varid, lats)
        if (status /= nf90_noerr) return
        status = nf90_put_var(ncid, precip_varid, precipitation)
        if (status /= nf90_noerr) return
        
        ! Close file
        status = nf90_close(ncid)
        
    end subroutine save_corrected_cmip6_precipitation_enhanced

    !> Enhanced save function for bias-corrected CMIP6 PET with comprehensive metadata
    subroutine save_corrected_cmip6_pet_enhanced(pet, lats, lons, scenario, model, &
                                                start_year, end_year, filename, status, &
                                                validation_metrics)
        real(dp), intent(in) :: pet(:,:,:)
        real(dp), intent(in) :: lats(:), lons(:)
        character(len=*), intent(in) :: scenario, model, filename
        integer, intent(in) :: start_year, end_year
        integer, intent(out) :: status
        type(bias_correction_validation_t), intent(in) :: validation_metrics
        
        integer :: ncid, lon_dimid, lat_dimid, time_dimid
        integer :: lon_varid, lat_varid, time_varid, pet_varid
        integer :: nlons, nlats, ntimes
        
        nlons = size(lons)
        nlats = size(lats)
        ntimes = size(pet, 3)
        
        ! Create NetCDF file
        status = nf90_create(filename, NF90_CLOBBER, ncid)
        if (status /= nf90_noerr) return
        
        ! Define dimensions
        status = nf90_def_dim(ncid, "longitude", nlons, lon_dimid)
        if (status /= nf90_noerr) return
        status = nf90_def_dim(ncid, "latitude", nlats, lat_dimid)
        if (status /= nf90_noerr) return
        status = nf90_def_dim(ncid, "time", ntimes, time_dimid)
        if (status /= nf90_noerr) return
        
        ! Define coordinate variables
        status = nf90_def_var(ncid, "longitude", NF90_DOUBLE, lon_dimid, lon_varid)
        if (status /= nf90_noerr) return
        status = nf90_def_var(ncid, "latitude", NF90_DOUBLE, lat_dimid, lat_varid)
        if (status /= nf90_noerr) return
        status = nf90_def_var(ncid, "time", NF90_DOUBLE, time_dimid, time_varid)
        if (status /= nf90_noerr) return
        
        ! Define data variable
        status = nf90_def_var(ncid, "potential_evapotranspiration", NF90_DOUBLE, &
                             [lon_dimid, lat_dimid, time_dimid], pet_varid)
        if (status /= nf90_noerr) return
        
        ! Add variable attributes
        call check_nc(nf90_put_att(ncid, pet_varid, "units", "mm/month"))
        call check_nc(nf90_put_att(ncid, pet_varid, "long_name", &
            "Bias-corrected monthly potential evapotranspiration (Hargreaves method)"))
        call check_nc(nf90_put_att(ncid, pet_varid, "standard_name", &
            "water_potential_evapotranspiration_amount"))
        call check_nc(nf90_put_att(ncid, pet_varid, "_FillValue", -999.0_dp))
        call check_nc(nf90_put_att(ncid, pet_varid, "missing_value", -999.0_dp))
        call check_nc(nf90_put_att(ncid, pet_varid, "method", "Hargreaves-Samani"))
        
        ! Add coordinate attributes
        call check_nc(nf90_put_att(ncid, lon_varid, "units", "degrees_east"))
        call check_nc(nf90_put_att(ncid, lat_varid, "units", "degrees_north"))
        call check_nc(nf90_put_att(ncid, time_varid, "units", "months since 1981-01-01"))
        
        ! Add comprehensive bias correction metadata
        call add_bias_correction_metadata(ncid, scenario, "pet", validation_metrics)
        
        ! End define mode
        status = nf90_enddef(ncid)
        if (status /= nf90_noerr) return
        
        ! Write data
        status = nf90_put_var(ncid, lon_varid, lons)
        if (status /= nf90_noerr) return
        status = nf90_put_var(ncid, lat_varid, lats)
        if (status /= nf90_noerr) return
        status = nf90_put_var(ncid, pet_varid, pet)
        if (status /= nf90_noerr) return
        
        ! Close file
        status = nf90_close(ncid)
        
    end subroutine save_corrected_cmip6_pet_enhanced

    !-------------------------------------------------------------------
    ! REGRIDDED CMIP6 DATA FUNCTIONS
    !-------------------------------------------------------------------
    
    !> Build file paths for regridded CMIP6 data
    subroutine build_regridded_cmip6_file_paths(base_dir, scenario, pr_file, tasmax_file, tasmin_file)
        !===========================================================================
        ! Build file paths for regridded CMIP6 data in Regridded_data structure
        ! Updated to match actual CDO regridding script output
        !===========================================================================
        character(len=*), intent(in) :: base_dir, scenario
        character(len=*), intent(out) :: pr_file, tasmax_file, tasmin_file
        
        character(len=256) :: scenario_dir
        
        ! Build scenario directory path for regridded data (matches CDO script structure)
        scenario_dir = trim(base_dir) // "/Processed_data/Regridded_data/CMIP6/" // trim(scenario)
        
        ! Build file paths for regridded files (simple naming from CDO script)
        pr_file = trim(scenario_dir) // "/precipitation_regridded.nc"
        tasmax_file = trim(scenario_dir) // "/tasmax_regridded.nc"
        tasmin_file = trim(scenario_dir) // "/tasmin_regridded.nc"
                      
    end subroutine build_regridded_cmip6_file_paths

    !> Read regridded CMIP6 NetCDF data
    subroutine read_regridded_cmip6_netcdf_data(pr_file, tasmax_file, tasmin_file, cmip6_data, success)
        !===========================================================================
        ! Read regridded CMIP6 data from Regridded_data structure
        ! This data is already on ERA5 0.1° grid and ready for bias correction
        !===========================================================================
        character(len=*), intent(in) :: pr_file, tasmax_file, tasmin_file
        type(cmip6_scenario_data_t), intent(inout) :: cmip6_data
        logical, intent(out) :: success
        
        integer :: ncid, varid, dimid
        integer :: status
        integer :: nlons, nlats, ntimes
        real(dp), allocatable :: temp_pr(:,:,:)
        
        success = .false.
        
        ! Read precipitation file first to get dimensions
        status = nf90_open(pr_file, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error opening regridded precipitation file: ", trim(pr_file)
            return
        end if
        
        ! Get dimensions (should match ERA5 grid: 56x36) - handle both naming conventions
        ! Try "longitude" first (regridded files), then "lon" (original files)
        status = nf90_inq_dimid(ncid, "longitude", dimid)
        if (status /= NF90_NOERR) then
            status = nf90_inq_dimid(ncid, "lon", dimid)
            if (status /= NF90_NOERR) then
                write(*,*) "❌ Error getting longitude dimension from regridded data"
                write(*,*) "   Tried both 'longitude' and 'lon' dimension names"
                status = nf90_close(ncid)
                return
            end if
        end if
        status = nf90_inquire_dimension(ncid, dimid, len=nlons)
        
        ! Try "latitude" first (regridded files), then "lat" (original files)
        status = nf90_inq_dimid(ncid, "latitude", dimid)
        if (status /= NF90_NOERR) then
            status = nf90_inq_dimid(ncid, "lat", dimid)
            if (status /= NF90_NOERR) then
                write(*,*) "❌ Error getting latitude dimension from regridded data"
                write(*,*) "   Tried both 'latitude' and 'lat' dimension names"
                status = nf90_close(ncid)
                return
            end if
        end if
        status = nf90_inquire_dimension(ncid, dimid, len=nlats)
        
        status = nf90_inq_dimid(ncid, "time", dimid)
        status = nf90_inquire_dimension(ncid, dimid, len=ntimes)
        
        ! Verify we have ERA5 grid dimensions (56 lon x 36 lat)
        if (nlons /= 56 .or. nlats /= 36) then
            write(*,*) "⚠️  WARNING: Regridded data dimensions don't match ERA5 grid"
            write(*,*) "   Expected: 56 x 36, Got: ", nlons, " x ", nlats
        end if
        
        ! Store dimensions
        cmip6_data%nlons = nlons
        cmip6_data%nlats = nlats
        cmip6_data%ntimes = ntimes
        
        ! Allocate arrays for raw data
        allocate(cmip6_data%precipitation(nlons, nlats, ntimes))
        allocate(cmip6_data%tasmax(nlons, nlats, ntimes))
        allocate(cmip6_data%tasmin(nlons, nlats, ntimes))
        allocate(cmip6_data%longitude(nlons))
        allocate(cmip6_data%latitude(nlats))
        allocate(cmip6_data%time_values(ntimes))
        allocate(temp_pr(nlons, nlats, ntimes))
        
        ! Allocate arrays for processed data
        allocate(cmip6_data%pet_raw(nlons, nlats, ntimes))
        allocate(cmip6_data%pet_corrected(nlons, nlats, ntimes))
        allocate(cmip6_data%precipitation_corrected(nlons, nlats, ntimes))
        allocate(cmip6_data%land_mask(nlons, nlats))
        allocate(cmip6_data%data_quality(nlons, nlats))
        
        ! Initialize status flags
        cmip6_data%land_mask = .true.  ! Initialize as all land (will be refined later)
        cmip6_data%data_quality = 1.0_dp  ! Initialize as perfect quality
        
        ! Read coordinate variables - handle both naming conventions
        ! Try "longitude" first (regridded files), then "lon" (original files)
        status = nf90_inq_varid(ncid, "longitude", varid)
        if (status /= NF90_NOERR) then
            status = nf90_inq_varid(ncid, "lon", varid)
        end if
        status = nf90_get_var(ncid, varid, cmip6_data%longitude)
        
        ! Try "latitude" first (regridded files), then "lat" (original files)
        status = nf90_inq_varid(ncid, "latitude", varid)
        if (status /= NF90_NOERR) then
            status = nf90_inq_varid(ncid, "lat", varid)
        end if
        status = nf90_get_var(ncid, varid, cmip6_data%latitude)
        
        status = nf90_inq_varid(ncid, "time", varid)
        status = nf90_get_var(ncid, varid, cmip6_data%time_values)
        
        ! Read precipitation data (in kg m-2 s-1, need to convert to mm/month)
        status = nf90_inq_varid(ncid, "pr", varid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error finding 'pr' variable in regridded data"
            status = nf90_close(ncid)
            return
        end if
        status = nf90_get_var(ncid, varid, temp_pr)
        
        ! First mask fill values before conversion
        where (temp_pr >= 1.e+19_dp) temp_pr = FILL_VALUE
        
        ! Convert from kg m-2 s-1 to mm/month 
        ! kg m-2 s-1 * 86400 s/day * ~30.44 days/month = mm/month
        ! Using 30.44 as average days per month (365.25/12)
        where (temp_pr /= FILL_VALUE)
            cmip6_data%precipitation = temp_pr * 86400.0_real64 * 30.44_real64
        elsewhere
            cmip6_data%precipitation = FILL_VALUE
        end where
        
        write(*,*) "   💧 Converted precipitation: kg/m²/s → mm/month (with fill value masking)"
        write(*,*) "      Valid data count:", count(temp_pr /= FILL_VALUE), "/", size(temp_pr)
        
        status = nf90_close(ncid)
        
        ! Read tasmax file
        status = nf90_open(tasmax_file, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error opening regridded tasmax file: ", trim(tasmax_file)
            return
        end if
        
        status = nf90_inq_varid(ncid, "tasmax", varid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error finding 'tasmax' variable in regridded data"
            status = nf90_close(ncid)
            return
        end if
        status = nf90_get_var(ncid, varid, cmip6_data%tasmax)
        
        ! Mask fill values in tasmax data
        where (cmip6_data%tasmax >= 1.e+19_dp) cmip6_data%tasmax = FILL_VALUE
        
        status = nf90_close(ncid)
        
        ! Read tasmin file
        status = nf90_open(tasmin_file, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error opening regridded tasmin file: ", trim(tasmin_file)
            return
        end if
        
        status = nf90_inq_varid(ncid, "tasmin", varid)
        if (status /= NF90_NOERR) then
            write(*,*) "❌ Error finding 'tasmin' variable in regridded data"
            status = nf90_close(ncid)
            return
        end if
        status = nf90_get_var(ncid, varid, cmip6_data%tasmin)
        
        ! Mask fill values in tasmin data
        where (cmip6_data%tasmin >= 1.e+19_dp) cmip6_data%tasmin = FILL_VALUE
        
        status = nf90_close(ncid)
        
        write(*,*) "✅ Successfully loaded regridded CMIP6 data"
        write(*,*) "   Dimensions: ", nlons, " x ", nlats, " x ", ntimes
        write(*,*) "   Longitude range: ", minval(cmip6_data%longitude), " to ", maxval(cmip6_data%longitude)
        write(*,*) "   Latitude range: ", minval(cmip6_data%latitude), " to ", maxval(cmip6_data%latitude)
        
        ! Create proper land mask based on valid precipitation data
        call create_cmip6_land_mask(cmip6_data)
        
        success = .true.
        
        deallocate(temp_pr)
        
    end subroutine read_regridded_cmip6_netcdf_data

    !> Create land mask for CMIP6 data based on valid precipitation values
    subroutine create_cmip6_land_mask(cmip6_data)
        type(cmip6_scenario_data_t), intent(inout) :: cmip6_data
        
        integer :: i, j, t, valid_count, total_count
        real(dp) :: mean_precip
        
        ! Initialize as no land
        cmip6_data%land_mask = .false.
        
        ! For each grid cell, check if there's reasonable precipitation data
        do j = 1, cmip6_data%nlats
            do i = 1, cmip6_data%nlons
                valid_count = 0
                total_count = 0
                mean_precip = 0.0_dp
                
                ! Check first 12 months of data for valid precipitation
                do t = 1, min(12, cmip6_data%ntimes)
                    if (cmip6_data%precipitation(i,j,t) /= FILL_VALUE .and. &     ! Not our fill value
                        cmip6_data%precipitation(i,j,t) >= 0.0_dp .and. &         ! Non-negative
                        cmip6_data%precipitation(i,j,t) < 1000.0_dp) then         ! Reasonable upper bound
                        valid_count = valid_count + 1
                        mean_precip = mean_precip + cmip6_data%precipitation(i,j,t)
                    end if
                    total_count = total_count + 1
                end do
                
                ! Mark as land if:
                ! 1. At least 50% of data is valid
                ! 2. Mean precipitation is reasonable (between 0.1 and 500 mm/month)
                if (valid_count >= total_count/2 .and. valid_count > 0) then
                    mean_precip = mean_precip / real(valid_count, dp)
                    if (mean_precip >= 0.1_dp .and. mean_precip <= 500.0_dp) then
                        cmip6_data%land_mask(i,j) = .true.
                    end if
                end if
            end do
        end do
        
        write(*,*) "   🗺️  Created CMIP6 land mask:"
        write(*,*) "      Land pixels: ", count(cmip6_data%land_mask), " / ", &
                   size(cmip6_data%land_mask), " (", &
                   100.0 * count(cmip6_data%land_mask) / size(cmip6_data%land_mask), "%)"
        
    end subroutine create_cmip6_land_mask

    !> Load regridded CMIP6 scenario data
    subroutine load_regridded_cmip6_scenario(data_directory, scenario, cmip6_data, status)
        !===========================================================================
        ! Load regridded CMIP6 scenario data (pr, tasmax, tasmin) from regridded NetCDF files
        ! This function mirrors load_raw_cmip6_scenario but uses regridded data
        !===========================================================================
        character(len=*), intent(in) :: data_directory
        character(len=*), intent(in) :: scenario
        type(cmip6_scenario_data_t), intent(out) :: cmip6_data
        integer, intent(out) :: status
        
        character(len=512) :: pr_file, tasmax_file, tasmin_file
        logical :: read_success
        
        status = -1  ! Initialize as error
        cmip6_data%scenario = trim(scenario)
        cmip6_data%is_corrected = .false.
        
        ! Set temporal metadata based on scenario
        if (trim(scenario) == "historical") then
            cmip6_data%start_year = 1981
            cmip6_data%end_year = 2014
        else
            cmip6_data%start_year = 2015
            cmip6_data%end_year = 2099
        end if
        
        ! Build file paths for regridded data
        call build_regridded_cmip6_file_paths(data_directory, scenario, pr_file, tasmax_file, tasmin_file)
        
        write(*,*) "📂 Loading regridded CMIP6 ", trim(scenario), " data..."
        write(*,*) "   Files:"
        write(*,*) "   - ", trim(pr_file)
        write(*,*) "   - ", trim(tasmax_file) 
        write(*,*) "   - ", trim(tasmin_file)
        
        ! Read the regridded NetCDF files (ERA5 grid compatible)
        call read_regridded_cmip6_netcdf_data(pr_file, tasmax_file, tasmin_file, cmip6_data, read_success)
        
        if (read_success) then
            status = 0  ! Success
            write(*,*) "✅ Regridded CMIP6 ", trim(scenario), " data loaded successfully"
            write(*,*) "   Dimensions: ", cmip6_data%nlons, " x ", cmip6_data%nlats, " x ", cmip6_data%ntimes
            write(*,*) "   📍 Data is on ERA5 0.1° grid for consistent bias correction"
        else
            status = -1  ! Error
            write(*,*) "❌ Failed to load regridded CMIP6 ", trim(scenario), " data"
            write(*,*) "   💡 Make sure regridding has been completed using:"
            write(*,*) "      ./regrid_cmip6_complete.sh"
        end if
        
    end subroutine load_regridded_cmip6_scenario

end module io_module
