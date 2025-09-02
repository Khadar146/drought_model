# Somaliland Drought Analysis Pipeline

## Author Information
**Author**: Khadar Dahir Abdisalan  
**Supervisor**: Sebastian Mutz ## Installation & Usage

### Prerequisites
```bash
conda create -n Drought fortran-compiler netcdf-fortran
conda activate Drought
conda install fpm
```

### Build & Execute
```bash
git clone https://github.com/Khadar146/drought_model.git
cd drought_model
fpm build --flag "-I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lnetcdff -lnetcdf"
fpm run --flag "-I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lnetcdff -lnetcdf"
```

## Naming Conventions

### File Naming Standards
**Pattern**: `{variable}_{dataset}_{scenario}_{timeperiod}.{ext}`

**Examples**:
- `spei_multiscale_era5_1981-2024.nc` - SPEI indices from ERA5 reference data
- `precipitation_corrected_ssp245_2015-2099.nc` - Bias-corrected precipitation projections  
- `evt_parameters_ensemble_1981-2014.nc` - Extreme value model parameters

### Variable Naming Conventions

#### Climate Variables
```fortran
! Precipitation and evapotranspiration
real(dp) :: precipitation(:,:,:)        ! Monthly precipitation [mm]
real(dp) :: pet(:,:,:)                   ! Potential evapotranspiration [mm]
real(dp) :: temperature_max(:,:,:)       ! Maximum temperature [K]
real(dp) :: temperature_min(:,:,:)       ! Minimum temperature [K]
```

#### Drought Indices
```fortran
! Multi-timescale drought indices
real(dp) :: spi_1(:,:,:)                 ! 1-month SPI
real(dp) :: spi_3(:,:,:)                 ! 3-month SPI
real(dp) :: spi_6(:,:,:)                 ! 6-month SPI
real(dp) :: spi_12(:,:,:)                ! 12-month SPI

real(dp) :: spei_1(:,:,:)                ! 1-month SPEI
real(dp) :: spei_3(:,:,:)                ! 3-month SPEI
real(dp) :: spei_6(:,:,:)                ! 6-month SPEI
real(dp) :: spei_12(:,:,:)               ! 12-month SPEI
```

#### Coordinate Arrays
```fortran
! Spatial and temporal coordinates
real(dp) :: latitudes(:)                 ! Latitude values [degrees_north]
real(dp) :: longitudes(:)                ! Longitude values [degrees_east]
real(dp) :: time_stamps(:)               ! Time values [seconds since 1970-01-01]
```

#### EVT Parameters
```fortran
! Extreme value theory parameters
real(dp) :: xi                           ! GPD shape parameter
real(dp) :: sigma                        ! GPD scale parameter  
real(dp) :: mu                           ! GPD location parameter (threshold)
real(dp) :: return_levels(3)             ! 10, 20, 50-year return levels
```

### Module Naming Standards
**Pattern**: `{purpose}_module`

**Examples**:
- `spi_module` - Standardized Precipitation Index calculations
- `spei_module` - Standardized Precipitation Evapotranspiration Index
- `io_module` - Input/output and data management
- `evt_fsml_comparison` - Extreme value theory analysis
- `bias_correction_clean` - Quantile mapping bias correction
- `projection_module` - Future scenario processing

### Function Naming Standards
**Pattern**: `{action}_{object}`

**Examples**:
```fortran
! Data loading functions
subroutine load_era5_climate_data(...)
subroutine load_corrected_cmip6_scenario(...)

! Calculation functions  
subroutine calculate_spi_timescales(...)
subroutine calculate_spei_timescales(...)

! Analysis functions
subroutine run_fsml_evt_comparison(...)
subroutine apply_quantile_mapping_correction(...)

! Save functions
subroutine save_processed_precipitation_data(...)
subroutine save_drought_indices(...)
```

### Data Type Naming Standards
**Pattern**: `{descriptor}_t`

**Examples**:
```fortran
! Structured data types
type :: era5_spei_results_t
type :: cmip6_spei_results_t  
type :: bias_correction_parameters_t
type :: fsml_comparison_results_t
type :: corrected_cmip6_scenarios_t
```

### Constants and Parameters
```fortran
! Precision and fill values
integer, parameter :: dp = real64
real(dp), parameter :: FILL_VALUE = -999.0_dp

! Physical constants
real(dp), parameter :: DAYS_PER_MONTH = 30.44_dp
real(dp), parameter :: M_TO_MM = 1000.0_dp

! Domain boundaries
real(dp), parameter :: LAT_MIN = 8.0_dp, LAT_MAX = 11.5_dp
real(dp), parameter :: LON_MIN = 43.0_dp, LON_MAX = 48.5_dp
```*: University of Glasgow, School of Geographical Sciences  
**Email**: khadardahir146@gmail.com

## Abstract

This pipeline implements a comprehensive climate analysis framework for assessing drought patterns and climate change impacts in Somaliland, Horn of Africa. The system processes historical climate data (1981-2024) and generates future drought projections (2015-2100) under multiple emission scenarios using standardized drought indices (SPI/SPEI) and extreme value theory.

### Research Objective
Quantify climate change impacts on drought frequency, duration, and severity in Somaliland through statistical analysis of historical observations and bias-corrected climate model projections.

### Study Domain
- **Region**: Somaliland (43.0°-48.5°E, 8.0°-11.5°N)
- **Area**: ~176,000 km² arid/semi-arid climate
- **Resolution**: 0.25° × 0.25° grid (~25 km)
- **Temporal Coverage**: 1981-2099 (44 years historical + 85 years projections)

## Methodology

### Core Statistical Framework
The pipeline implements a rigorous statistical approach using the **FSML (Fortran Statistical and Machine Learning) library** for all probability distributions and hypothesis testing.

### Data Sources & Processing
- **ERA5-Land**: ECMWF reanalysis (1981-2024) - 0.25° resolution reference dataset
- **CMIP6**: Multi-model ensemble (historical + SSP1-2.6, SSP2-4.5, SSP5-8.5)
- **Variables**: Precipitation, temperature (max/min), potential evapotranspiration (Hargreaves method)

### Analytical Workflow

#### 1. Drought Index Calculation
**SPI (Standardized Precipitation Index)**
- **Distribution**: Gamma distribution fitting using FSML maximum likelihood estimation
- **Timescales**: 1, 3, 6, 12-month aggregations
- **Standardization**: CDF → Normal inverse transformation
- **Implementation**: `spi_module.f90` with FSML statistical functions

**SPEI (Standardized Precipitation Evapotranspiration Index)**
- **Water Balance**: P - PET (precipitation minus potential evapotranspiration)
- **Distribution**: Gamma distribution for positive water balance, empirical for negative
- **PET Method**: Hargreaves temperature-based calculation
- **Implementation**: `spei_module.f90` with structured data types

#### 2. Bias Correction
**Quantile Mapping Method** (Themeßl et al., 2011; Cannon et al., 2015)
- **Precipitation**: Gamma distribution quantile mapping (shape/scale parameters)
- **PET**: Normal distribution quantile mapping (mean/std parameters)
- **Reference Period**: 1981-2014 (ERA5-CMIP6 overlap)
- **Validation**: FSML correlation, RMSE, bias metrics
- **Implementation**: `bias_correction_clean.f90` with FSML distribution functions

#### 3. Future Projections
**Climate Scenarios**
- **SSP1-2.6**: Low emissions (~1.5°C warming by 2100)
- **SSP2-4.5**: Moderate emissions (~2.5°C warming)
- **SSP5-8.5**: High emissions (~4.5°C warming)
- **Time Periods**: 2015-2099 bias-corrected projections
- **Implementation**: `projection_module.f90` with structured scenario handling

#### 4. Extreme Value Theory (EVT)
**FSML-based Climate Change Testing**
- **Distribution**: Generalized Pareto Distribution (GPD) for drought exceedances
- **Parameter Estimation**: Method of Moments with FSML validation
- **Statistical Tests**: 
  - Paired t-tests for severity changes
  - Two-sample t-tests for frequency changes
  - Confidence intervals for parameter reliability
  - Extreme drought intensification tests
- **Implementation**: `evt_fsml_comparison.f90` with comprehensive hypothesis testing

### Technical Standards
- **Programming**: Fortran 2008+ with modular architecture
- **Statistics**: FSML library for all probability distributions and tests
- **Data Format**: NetCDF-4 with CF-1.8 conventions
- **Quality Control**: Comprehensive validation and error handling

## Analysis Pipeline

### Phase 1: Baseline Drought Indices (ERA5, 1981–2024)
- **1A**: ERA5 climate data loading (precipitation + PET)
- **1B**: Data preprocessing and quality control
- **1C**: SPI calculation using FSML statistical library
- **1D**: SPEI calculation with Hargreaves PET method

### Phase 2: Validation & Bias Correction
- **2A**: Historical calibration using CMIP6 vs ERA5 (1981-2014)
- **2B**: Future scenario processing (SSP1-2.6, SSP2-4.5, SSP5-8.5)

### Phase 3: Future Projections (2015-2100)
- Bias-corrected CMIP6 climate scenarios
- Multi-scale drought index calculations
- Climate change signal analysis

### Phase 4: Statistical Analysis (Planned)
- Trend analysis (Mann-Kendall, Sen's slope)
- Change point detection (PETTITT tests)
- Regional vulnerability assessment

## Core Modules

### Statistical Computing
- **`spi_module.f90`**: SPI calculation using FSML gamma distribution fitting
- **`spei_module.f90`**: SPEI calculation with structured data types and water balance
- **`bias_correction_clean.f90`**: Quantile mapping with FSML distribution functions
- **`evt_fsml_comparison.f90`**: Extreme value analysis with comprehensive hypothesis testing

### Data Management
- **`io_module.f90`**: NetCDF I/O with CF-1.8 compliance and metadata handling
- **`projection_module.f90`**: Future scenario processing and CMIP6 data management

### Pipeline Execution
- **`app/main.f90`**: Complete analysis workflow (ERA5 → bias correction → projections → EVT)

## Key Features

### Statistical Rigor
- **FSML Integration**: Complete statistical library for probability distributions
- **Hypothesis Testing**: Four-test framework for climate change detection
- **Validation Metrics**: Correlation, RMSE, bias assessment with confidence intervals
- **Quality Assurance**: Comprehensive error handling and data validation

### Climate Science Standards
- **CF Conventions**: Fully compliant NetCDF metadata
- **WMO Guidelines**: Standard drought index calculation methods
- **CMIP6 Compatibility**: Direct processing of climate model outputs
- **Bias Correction**: State-of-the-art quantile mapping techniques

## Installation & Usage

### Prerequisites
```bash
conda create -n Drought fortran-compiler netcdf-fortran
conda activate Drought
conda install fpm
```

### Build & Execute
```bash
git clone https://github.com/Khadar146/drought_model.git
cd drought_model
fpm build --flag "-I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lnetcdff -lnetcdf"
fpm run --flag "-I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lnetcdff -lnetcdf"
```

## Results Summary

### Analytical Capabilities

#### Drought Index Products
- **SPI**: Multi-timescale (1,3,6,12-month) with FSML gamma distribution fitting
- **SPEI**: Water balance approach with temperature-based PET calculation
- **Coverage**: Complete Somaliland domain (43.0°-48.5°E, 8.0°-11.5°N)
- **Resolution**: 0.25° grid (~25 km spacing)

#### Climate Change Detection
- **Historical Trends**: 44-year ERA5 baseline (1981-2024)
- **Future Projections**: Bias-corrected CMIP6 scenarios (2015-2099)
- **Statistical Testing**: FSML-based hypothesis testing framework
- **Extreme Events**: GPD-based return period analysis

#### Quality Metrics
| Component | Validation | Standard |
|-----------|------------|----------|
| SPI/SPEI Calculation | FSML statistical functions | WMO guidelines |
| Bias Correction | r > 0.7, RMSE < 20% | ERA5 reference |
| EVT Analysis | 4-test hypothesis framework | Climate change detection |
| Data Coverage | >95% spatial completeness | CF-1.8 compliance |

### Research Applications

#### Climate Impact Assessment
- **Drought Intensification**: Quantified changes in severity and frequency
- **Regional Vulnerability**: Spatial patterns of drought response
- **Scenario Analysis**: Multi-emission pathway comparison
- **Extreme Events**: Return period estimation for drought planning

#### Technical Validation
- **Reference Dataset**: ERA5-Land as observational baseline
- **Statistical Robustness**: FSML library ensures computational accuracy
- **Uncertainty Quantification**: Confidence intervals and ensemble statistics
- **Reproducibility**: Fully documented methodology and open-source code

## License & Citation

**License**: MIT License
 https://github.com/Khadar146/drought_model

## Acknowledgments

- ECMWF ERA5-Land reanalysis data
- CMIP6 modeling groups and data providers
- FSML statistical library (Sebastian Mutz)

---
**Version**: 1.0.0 | **Status**: CONTINOUS DEVELOPMENT | **Updated**: September 2025
