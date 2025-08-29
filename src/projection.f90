!=======================================================================
! PROJECTION MODULE FOR FUTURE DROUGHT ANALYSIS
!=======================================================================
! Purpose:
!   Computes future drought projections using bias-corrected CMIP6 data
!   under multiple emission scenarios (SSP1-2.6, SSP2-4.5, SSP5-8.5).
!   Generates SPEI indices for climate change impact assessment.
!
! Methodology:
!   1. Load bias-corrected CMIP6 precipitation and temperature data
!   2. Calculate potential evapotranspiration using temperature-based methods
!   3. Compute water balance (P - PET) for each scenario
!   4. Calculate SPEI at multiple timescales using reference period parameters
!   5. Analyse trends and changes in drought characteristics
!
! Emission Scenarios:
!   SSP1-2.6: Low emissions, strong mitigation (~1.5°C warming)
!   SSP2-4.5: Moderate emissions, middle-of-the-road (~2.5°C warming)
!   SSP5-8.5: High emissions, fossil fuel development (~4.5°C warming)
!
! Time Periods:
!   Historical: 1981-2014 (reference for bias correction)
!   Near-term: 2021-2040
!   Mid-century: 2041-2070
!   End-century: 2071-2100
!
! Author: Khadar
! Date: August 2025
!=======================================================================

module projection_module
    use, intrinsic :: iso_fortran_env, only: real64
    use io_module, only: load_corrected_cmip6_scenario, save_spei_multiscale, save_spi_multiscale, &
                         save_future_drought_indices, corrected_cmip6_scenarios_t
    use spei_module, only: calculate_spei_timescales
    use spi_module, only: calculate_spi_timescales
    implicit none
    
    private
    public :: calculate_future_spei, calculate_future_spei_structured, cmip6_spei_results_t
    
    ! Precision parameter
    integer, parameter :: dp = real64
    
    ! Emission scenarios
    integer, parameter :: n_scenarios = 3
    character(len=10), parameter :: scenarios(n_scenarios) = &
        ['ssp126    ', 'ssp245    ', 'ssp585    ']
    
    ! Time periods for analysis
    integer, parameter :: historical_start = 1981
    integer, parameter :: historical_end = 2014
    integer, parameter :: near_term_start = 2021
    integer, parameter :: near_term_end = 2040
    integer, parameter :: mid_century_start = 2041
    integer, parameter :: mid_century_end = 2070
    integer, parameter :: end_century_start = 2071
    integer, parameter :: end_century_end = 2100
    
    ! SPEI timescales
    integer, parameter :: n_timescales = 4
    integer, parameter :: spei_timescales(n_timescales) = [1, 3, 6, 12]

    ! Data structure to hold CMIP6 SPEI results for EVT analysis
    type, public :: cmip6_spei_results_t
        ! SPEI arrays for each scenario [lon, lat, time]
        real(dp), allocatable :: ssp126_spei_3(:,:,:), ssp126_spei_6(:,:,:), ssp126_spei_12(:,:,:)
        real(dp), allocatable :: ssp245_spei_3(:,:,:), ssp245_spei_6(:,:,:), ssp245_spei_12(:,:,:)
        real(dp), allocatable :: ssp585_spei_3(:,:,:), ssp585_spei_6(:,:,:), ssp585_spei_12(:,:,:)
        
        ! Coordinate arrays
        real(dp), allocatable :: latitudes(:), longitudes(:), time_values(:)
        
        ! Metadata
        integer :: start_year, end_year
        logical :: data_available
    end type cmip6_spei_results_t

contains

    !-------------------------------------------------------------------
    ! CALCULATE FUTURE SPEI - STRUCTURED VERSION (consistent with ERA5)
    !-------------------------------------------------------------------
    subroutine calculate_future_spei_structured(corrected_scenarios, cmip6_results, success)
        !===========================================================================
        ! Memory-based SPEI calculation for all CMIP6 scenarios
        ! Accepts structured corrected data, returns structured SPEI results
        ! Consistent with ERA5 memory-based approach
        !===========================================================================
        type(corrected_cmip6_scenarios_t), intent(in) :: corrected_scenarios
        type(cmip6_spei_results_t), intent(out) :: cmip6_results
        logical, intent(out) :: success
        
        logical :: ssp126_success, ssp245_success, ssp585_success
        integer :: nlon, nlat, ntime
        
        success = .true.
        cmip6_results%data_available = .false.
        
        print *, ""
        print *, "=============================================="
        print *, "   STRUCTURED FUTURE DROUGHT PROJECTIONS"
        print *, "=============================================="
        print *, "Computing future SPEI for emission scenarios:"
        print *, "(Using memory-based structured CMIP6 data)"
        print *, ""
        
        ! Verify input data is available
        if (.not. corrected_scenarios%data_available) then
            print *, "❌ ERROR: Corrected CMIP6 scenarios data not available"
            success = .false.
            return
        end if
        
        ! Get dimensions
        nlon = size(corrected_scenarios%ssp126_precip, 1)
        nlat = size(corrected_scenarios%ssp126_precip, 2) 
        ntime = size(corrected_scenarios%ssp126_precip, 3)
        
        ! Allocate SPEI result arrays
        allocate(cmip6_results%ssp126_spei_3(nlon, nlat, ntime))
        allocate(cmip6_results%ssp126_spei_6(nlon, nlat, ntime))
        allocate(cmip6_results%ssp126_spei_12(nlon, nlat, ntime))
        allocate(cmip6_results%ssp245_spei_3(nlon, nlat, ntime))
        allocate(cmip6_results%ssp245_spei_6(nlon, nlat, ntime))
        allocate(cmip6_results%ssp245_spei_12(nlon, nlat, ntime))
        allocate(cmip6_results%ssp585_spei_3(nlon, nlat, ntime))
        allocate(cmip6_results%ssp585_spei_6(nlon, nlat, ntime))
        allocate(cmip6_results%ssp585_spei_12(nlon, nlat, ntime))
        
        ! Copy coordinate arrays
        allocate(cmip6_results%latitudes(size(corrected_scenarios%latitudes)))
        allocate(cmip6_results%longitudes(size(corrected_scenarios%longitudes)))
        allocate(cmip6_results%time_values(size(corrected_scenarios%time_values)))
        cmip6_results%latitudes = corrected_scenarios%latitudes
        cmip6_results%longitudes = corrected_scenarios%longitudes
        cmip6_results%time_values = corrected_scenarios%time_values
        
        ! Set metadata
        cmip6_results%start_year = corrected_scenarios%start_year
        cmip6_results%end_year = corrected_scenarios%end_year
        
        ! Process each scenario using the structured data
        call process_single_scenario_structured("ssp126", &
                                               corrected_scenarios%ssp126_precip, corrected_scenarios%ssp126_pet, &
                                               cmip6_results%ssp126_spei_3, cmip6_results%ssp126_spei_6, cmip6_results%ssp126_spei_12, &
                                               corrected_scenarios%start_year, corrected_scenarios%end_year, ssp126_success)
        
        call process_single_scenario_structured("ssp245", &
                                               corrected_scenarios%ssp245_precip, corrected_scenarios%ssp245_pet, &
                                               cmip6_results%ssp245_spei_3, cmip6_results%ssp245_spei_6, cmip6_results%ssp245_spei_12, &
                                               corrected_scenarios%start_year, corrected_scenarios%end_year, ssp245_success)
        
        call process_single_scenario_structured("ssp585", &
                                               corrected_scenarios%ssp585_precip, corrected_scenarios%ssp585_pet, &
                                               cmip6_results%ssp585_spei_3, cmip6_results%ssp585_spei_6, cmip6_results%ssp585_spei_12, &
                                               corrected_scenarios%start_year, corrected_scenarios%end_year, ssp585_success)
        
        ! Check overall success
        success = ssp126_success .and. ssp245_success .and. ssp585_success
        
        if (success) then
            cmip6_results%data_available = .true.
            print *, ""
            print *, "=============================================="
            print *, "   PHASE 3 COMPLETE: FUTURE PROJECTIONS"
            print *, "=============================================="
            print *, ""
        else
            print *, "❌ ERROR: Some CMIP6 scenarios failed to process"
            cmip6_results%data_available = .false.
        end if
        
    end subroutine calculate_future_spei_structured

    !-------------------------------------------------------------------
    ! CALCULATE FUTURE SPEI FOR ALL SCENARIOS - MODULAR VERSION
    !-------------------------------------------------------------------
    subroutine calculate_future_spei(ssp126_precip, ssp126_pet, ssp245_precip, ssp245_pet, &
                                    ssp585_precip, ssp585_pet, corrected_lat, corrected_lon, &
                                    corrected_time, corrected_start_year, corrected_end_year, &
                                    cmip6_results, success)
        implicit none
        real(dp), intent(in) :: ssp126_precip(:,:,:), ssp126_pet(:,:,:)
        real(dp), intent(in) :: ssp245_precip(:,:,:), ssp245_pet(:,:,:)
        real(dp), intent(in) :: ssp585_precip(:,:,:), ssp585_pet(:,:,:)
        real(dp), intent(in) :: corrected_lat(:), corrected_lon(:), corrected_time(:)
        integer, intent(in) :: corrected_start_year, corrected_end_year
        type(cmip6_spei_results_t), intent(out) :: cmip6_results
        logical, intent(out) :: success
        
        logical :: ssp126_success, ssp245_success, ssp585_success
        
        success = .true.  ! Start optimistic
        cmip6_results%data_available = .false.
        
        ! Initialize metadata
        cmip6_results%start_year = corrected_start_year
        cmip6_results%end_year = corrected_end_year
        
        ! Copy coordinate arrays
        allocate(cmip6_results%latitudes(size(corrected_lat)))
        allocate(cmip6_results%longitudes(size(corrected_lon)))
        allocate(cmip6_results%time_values(size(corrected_time)))
        cmip6_results%latitudes = corrected_lat
        cmip6_results%longitudes = corrected_lon
        cmip6_results%time_values = corrected_time
        
        write(*,*) ""
        write(*,*) "=============================================="
        write(*,*) "   PHASE 3: FUTURE DROUGHT PROJECTIONS"
        write(*,*) "=============================================="
        write(*,*) "Computing future SPEI for emission scenarios:"
        write(*,*) "(Using pre-loaded bias-corrected CMIP6 data)"
        write(*,*) ""
        
        ! Process SSP126 scenario with provided data
        write(*,*) "📊 Processing scenario: ssp126"
        call process_single_scenario_modular("ssp126", ssp126_precip, ssp126_pet, &
                                            corrected_lat, corrected_lon, corrected_time, &
                                            corrected_start_year, corrected_end_year, &
                                            cmip6_results%ssp126_spei_3, cmip6_results%ssp126_spei_6, &
                                            cmip6_results%ssp126_spei_12, ssp126_success)
        if (ssp126_success) then
            write(*,*) "✅ ssp126 completed successfully"
        else
            write(*,*) "❌ ssp126 processing failed"
            success = .false.
        end if
        write(*,*) ""
        
        ! Process SSP245 scenario with provided data
        write(*,*) "📊 Processing scenario: ssp245"
        call process_single_scenario_modular("ssp245", ssp245_precip, ssp245_pet, &
                                            corrected_lat, corrected_lon, corrected_time, &
                                            corrected_start_year, corrected_end_year, &
                                            cmip6_results%ssp245_spei_3, cmip6_results%ssp245_spei_6, &
                                            cmip6_results%ssp245_spei_12, ssp245_success)
        if (ssp245_success) then
            write(*,*) "✅ ssp245 completed successfully"
        else
            write(*,*) "❌ ssp245 processing failed"
            success = .false.
        end if
        write(*,*) ""
        
        ! Process SSP585 scenario with provided data
        write(*,*) "📊 Processing scenario: ssp585"
        call process_single_scenario_modular("ssp585", ssp585_precip, ssp585_pet, &
                                            corrected_lat, corrected_lon, corrected_time, &
                                            corrected_start_year, corrected_end_year, &
                                            cmip6_results%ssp585_spei_3, cmip6_results%ssp585_spei_6, &
                                            cmip6_results%ssp585_spei_12, ssp585_success)
        if (ssp585_success) then
            write(*,*) "✅ ssp585 completed successfully"
        else
            write(*,*) "❌ ssp585 processing failed"
            success = .false.
        end if
        write(*,*) ""
        
        ! Mark data as available if all scenarios succeeded
        cmip6_results%data_available = success
        
        write(*,*) "=============================================="
        write(*,*) "   PHASE 3 COMPLETE: FUTURE PROJECTIONS"
        write(*,*) "=============================================="
        write(*,*) ""
        
    end subroutine calculate_future_spei

    !-------------------------------------------------------------------
    ! PROCESS SINGLE SCENARIO - STRUCTURED VERSION
    !-------------------------------------------------------------------
    subroutine process_single_scenario_structured(scenario_name, precip_data, pet_data, &
                                                 spei_3, spei_6, spei_12, start_year, end_year, success)
        !===========================================================================
        ! Process a single scenario using pre-loaded structured data
        ! Eliminates file I/O during projection calculations
        !===========================================================================
        character(len=*), intent(in) :: scenario_name
        real(dp), dimension(:,:,:), intent(in) :: precip_data, pet_data
        real(dp), dimension(:,:,:), intent(out) :: spei_3, spei_6, spei_12
        integer, intent(in) :: start_year, end_year
        logical, intent(out) :: success
        
        print *, "  ➤ Processing ", trim(scenario_name), " (", start_year, "-", end_year, ")"
        
        ! Calculate SPEI using the structured data directly
        ! Need to declare spei_1 (not used but required by interface)
        block
            real(dp), allocatable :: spei_1(:,:,:)
            integer :: nlon, nlat, ntime
            
            nlon = size(precip_data, 1)
            nlat = size(precip_data, 2)
            ntime = size(precip_data, 3)
            
            allocate(spei_1(nlon, nlat, ntime))
            
            call calculate_spei_timescales(precip_data, pet_data, nlon, nlat, ntime, &
                                         spei_1, spei_3, spei_6, spei_12)
            
            deallocate(spei_1)
            success = .true.
        end block
        
        if (success) then
            print *, "    ✓ ", trim(scenario_name), " SPEI calculation completed"
        else
            print *, "    ❌ ", trim(scenario_name), " SPEI calculation failed"
        end if
        
    end subroutine process_single_scenario_structured

    !-------------------------------------------------------------------
    ! PROCESS SINGLE EMISSION SCENARIO - MODULAR VERSION (USES PROVIDED DATA)
    !-------------------------------------------------------------------
    subroutine process_single_scenario_modular(scenario_name, precipitation, pet, &
                                              latitudes, longitudes, time_values, &
                                              start_year, end_year, &
                                              spei_3_out, spei_6_out, spei_12_out, success)
        implicit none
        character(len=*), intent(in) :: scenario_name
        real(dp), intent(in) :: precipitation(:,:,:), pet(:,:,:)
        real(dp), intent(in) :: latitudes(:), longitudes(:), time_values(:)
        integer, intent(in) :: start_year, end_year
        real(dp), allocatable, intent(out) :: spei_3_out(:,:,:), spei_6_out(:,:,:), spei_12_out(:,:,:)
        logical, intent(out) :: success
        
        integer :: nlons, nlats, ntimes, status
        
        ! Drought indices arrays
        real(dp), allocatable :: spi_1(:,:,:), spi_3(:,:,:), spi_6(:,:,:), spi_12(:,:,:)
        real(dp), allocatable :: spei_1(:,:,:), spei_3(:,:,:), spei_6(:,:,:), spei_12(:,:,:)
        
        success = .false.
        
        ! Get dimensions from provided data
        nlons = size(precipitation, 1)
        nlats = size(precipitation, 2)
        ntimes = size(precipitation, 3)
        
        write(*,*) "   ✅ Using pre-loaded CMIP6 data: ", nlons, "×", nlats, "×", ntimes
        write(*,*) "   Time period: ", start_year, "-", end_year
        
        ! Allocate drought indices arrays
        allocate(spi_1(nlons, nlats, ntimes))
        allocate(spi_3(nlons, nlats, ntimes))
        allocate(spi_6(nlons, nlats, ntimes))
        allocate(spi_12(nlons, nlats, ntimes))
        
        allocate(spei_1(nlons, nlats, ntimes))
        allocate(spei_3(nlons, nlats, ntimes))
        allocate(spei_6(nlons, nlats, ntimes))
        allocate(spei_12(nlons, nlats, ntimes))
        
        ! Calculate SPI indices using provided precipitation data
        write(*,*) "   📊 Calculating SPI indices..."
        call calculate_spi_timescales(precipitation, spi_1, spi_3, spi_6, spi_12, status)
        
        if (status /= 0) then
            write(*,*) "   ❌ Failed to calculate SPI for ", trim(scenario_name)
            return
        end if
        write(*,*) "   ✅ SPI calculated successfully"
        
        ! Calculate SPEI indices using provided precipitation and PET data
        write(*,*) "   📊 Calculating SPEI indices..."
        call calculate_spei_timescales(precipitation, pet, nlons, nlats, ntimes, &
                                      spei_1, spei_3, spei_6, spei_12)
        write(*,*) "   ✅ SPEI calculated successfully"
        
        ! Save drought indices
        write(*,*) "   💾 Saving future drought indices for ", trim(scenario_name), "..."
        call save_future_drought_indices(scenario_name, spi_1, spi_3, spi_6, spi_12, &
                                        spei_1, spei_3, spei_6, spei_12, &
                                        latitudes, longitudes, time_values, &
                                        start_year, end_year, status)
        
        if (status == 0) then
            write(*,*) "   ✅ Future drought indices saved successfully"
            success = .true.
        else
            write(*,*) "   ❌ Failed to save drought indices"
            return
        end if
        
        ! Return SPEI data for EVT analysis (allocate and copy)
        allocate(spei_3_out(nlons, nlats, ntimes))
        allocate(spei_6_out(nlons, nlats, ntimes))
        allocate(spei_12_out(nlons, nlats, ntimes))
        
        spei_3_out = spei_3
        spei_6_out = spei_6 
        spei_12_out = spei_12
        
        ! Cleanup SPI and SPEI-1 data (we don't need them for EVT)
        deallocate(spi_1, spi_3, spi_6, spi_12)
        deallocate(spei_1)
        ! Keep spei_3, spei_6, spei_12 allocated as they're returned
        
    end subroutine process_single_scenario_modular

    !-------------------------------------------------------------------
    ! PROCESS SINGLE EMISSION SCENARIO - LEGACY VERSION (LOADS DATA INTERNALLY)
    !-------------------------------------------------------------------
    subroutine process_single_scenario(scenario_name, success)
        implicit none
        character(len=*), intent(in) :: scenario_name
        logical, intent(out) :: success
        
        ! CMIP6 data arrays
        real(dp), allocatable :: precipitation(:,:,:), pet(:,:,:)
        real(dp), allocatable :: latitudes(:), longitudes(:), time_values(:)
        integer :: start_year, end_year, status
        integer :: nlons, nlats, ntimes
        
        ! Drought indices arrays
        real(dp), allocatable :: spi_1(:,:,:), spi_3(:,:,:), spi_6(:,:,:), spi_12(:,:,:)
        real(dp), allocatable :: spei_1(:,:,:), spei_3(:,:,:), spei_6(:,:,:), spei_12(:,:,:)
        
        success = .false.
        
        ! Step 1: Load bias-corrected CMIP6 data
        write(*,*) "   📂 Loading bias-corrected CMIP6 data..."
        call load_corrected_cmip6_scenario(scenario_name, precipitation, pet, latitudes, longitudes, &
                                          time_values, start_year, end_year, status)
        
        if (status /= 0) then
            write(*,*) "   ❌ Failed to load CMIP6 data for ", trim(scenario_name)
            return
        end if
        
        ! Get dimensions
        nlons = size(precipitation, 1)
        nlats = size(precipitation, 2)
        ntimes = size(precipitation, 3)
        
        write(*,*) "   ✅ CMIP6 data loaded: ", nlons, "×", nlats, "×", ntimes
        write(*,*) "   Time period: ", start_year, "-", end_year
        
        ! Step 2: Allocate drought indices arrays
        allocate(spi_1(nlons, nlats, ntimes))
        allocate(spi_3(nlons, nlats, ntimes))
        allocate(spi_6(nlons, nlats, ntimes))
        allocate(spi_12(nlons, nlats, ntimes))
        
        allocate(spei_1(nlons, nlats, ntimes))
        allocate(spei_3(nlons, nlats, ntimes))
        allocate(spei_6(nlons, nlats, ntimes))
        allocate(spei_12(nlons, nlats, ntimes))
        
        ! Step 3: Calculate SPI indices
        write(*,*) "   📊 Calculating SPI indices..."
        call calculate_spi_timescales(precipitation, spi_1, spi_3, spi_6, spi_12, status)
        
        if (status /= 0) then
            write(*,*) "   ❌ Failed to calculate SPI for ", trim(scenario_name)
            return
        end if
        write(*,*) "   ✅ SPI calculated successfully"
        
        ! Step 4: Calculate SPEI indices
        write(*,*) "   📊 Calculating SPEI indices..."
        call calculate_spei_timescales(precipitation, pet, nlons, nlats, ntimes, &
                                      spei_1, spei_3, spei_6, spei_12)
        write(*,*) "   ✅ SPEI calculated successfully"
        
        ! Step 5: Save drought indices (placeholder for now)
        write(*,*) "   💾 Saving future drought indices for ", trim(scenario_name), "..."
        call save_future_drought_indices(scenario_name, spi_1, spi_3, spi_6, spi_12, &
                                        spei_1, spei_3, spei_6, spei_12, &
                                        latitudes, longitudes, time_values, &
                                        start_year, end_year, status)
        
        if (status == 0) then
            write(*,*) "   ✅ Future drought indices saved successfully"
            success = .true.
        else
            write(*,*) "   ❌ Failed to save drought indices"
        end if
        
        ! Cleanup
        deallocate(precipitation, pet, latitudes, longitudes, time_values)
        deallocate(spi_1, spi_3, spi_6, spi_12)
        deallocate(spei_1, spei_3, spei_6, spei_12)
        
    end subroutine process_single_scenario

    !-------------------------------------------------------------------
    ! COMPUTE FUTURE WATER BALANCE (to be implemented)
    !-------------------------------------------------------------------
    ! subroutine compute_future_water_balance(precip, pet, water_balance)
    !     implicit none
    !     real(dp), intent(in) :: precip(:,:,:)
    !     real(dp), intent(in) :: pet(:,:,:)
    !     real(dp), intent(out) :: water_balance(:,:,:)
    !     
    !     ! Calculate water balance for future projections
    !     
    ! end subroutine compute_future_water_balance
    
    !-------------------------------------------------------------------
    ! CALCULATE FUTURE SPEI USING REFERENCE PARAMETERS (to be implemented)
    !-------------------------------------------------------------------
    ! subroutine calculate_future_spei_values(water_balance, spei_values)
    !     implicit none
    !     real(dp), intent(in) :: water_balance(:,:,:)
    !     real(dp), intent(out) :: spei_values(:,:,:,:)  ! (time, lat, lon, timescale)
    !     
    !     ! Use distribution parameters from historical period
    !     ! Apply to future water balance data
    !     
    ! end subroutine calculate_future_spei_values
    
    !-------------------------------------------------------------------
    ! ANALYSE DROUGHT TRENDS (to be implemented)
    !-------------------------------------------------------------------
    ! subroutine analyse_drought_trends(spei_data, scenario_name)
    !     implicit none
    !     real(dp), intent(in) :: spei_data(:,:,:,:)
    !     character(len=*), intent(in) :: scenario_name
    !     
    !     ! Analyse trends in drought frequency, severity, duration
    !     ! Compare different time periods
    !     ! Generate trend statistics
    !     
    ! end subroutine analyse_drought_trends

end module projection_module
