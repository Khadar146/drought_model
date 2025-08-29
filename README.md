# Somaliland Drought Analysis Pipeline

[![License: MIT]

## Project Overview

The Somaliland Drought Analysis Pipeline is a comprehensive climate analysis framework designed to assess drought patterns and climate change impacts in Somaliland, Horn of Africa. The pipeline processes historical climate data (1981-2024) and generates future drought projections (2015-2100) under multiple emission scenarios.

### Research Hypothesis
**"Climate change will increase drought frequency, duration, and severity in Somaliland."**

### Geographic Domain
- **Region**: Somaliland (43.0°-48.5°E, 8.0°-11.5°N)
- **Area**: ~176,000 km² of arid/semi-arid climate
- **Resolution**: 0.1° × 0.1° grid (~11 km spacing)
- **Time Period**: 44 years (1981-2024) + projections to 2100

## Quick Start

### Prerequisites
```bash
# Create conda environment
conda create -n Drought fortran-compiler netcdf-fortran
conda activate Drought

# Install FPM (Fortran Package Manager)
conda install fpm
```

### Build and Run
```bash
# Clone repository
git clone https://github.com/Khadar146/drought_model.git
cd drought_model

# Build pipeline

fpm build 
or
fpm build --flag "-I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lnetcdff -lnetcdf"

# Run full analysis
fpm run
or
fpm run --flag "-I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lnetcdff -lnetcdf"
```

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

## Documentation

### Comprehensive Guides
- **[Bias Correction and Projection Methods](docs/BIAS_CORRECTION_AND_PROJECTION_METHODS.md)** - Detailed methodology, statistical methods, and limitations
- **[Technical Specifications](docs/TECHNICAL_SPECIFICATIONS.md)** - Data formats, coordinate systems, and processing parameters
- **[Quick Reference Tables](docs/QUICK_REFERENCE.md)** - Summary tables, error codes, and troubleshooting

###  Specialized Documentation
- **[SPEI Documentation](docs/SPEI_DOCUMENTATION.md)** - Standardized Precipitation Evapotranspiration Index
- **[SPI Documentation](docs/SPI_DOCUMENTATION.md)** - Standardized Precipitation Index
- **[IO Module Documentation](docs/IO_Module_Documentation.md)** - NetCDF data handling
- **[NetCDF Data Layout](docs/NetCDF_Data_Layout_Explanation.md)** - File formats and structure

###  Validation Reports
- **[SPI Validation Report](docs/SPI_VALIDATION_REPORT.md)** - Quality control and verification
- **[Drought Indices Validation](docs/DROUGHT_INDICES_VALIDATION_REPORT.md)** - Cross-validation results
- **[Missing Value Analysis](docs/Missing_Value_Analysis.md)** - Data quality assessment

## Technical Stack

### Core Libraries
| Library | Purpose | Version |
|---------|---------|---------|
| **[FSML](https://github.com/sebastian-mutz/fsml)** | Statistical methods (Gamma distribution fitting) | Latest |
| **NetCDF-Fortran** | Climate data I/O | ≥ 4.5.0 |
| **Fortran Standard Library** | Utility functions | ≥ 0.2.0 |

### Key Features
- **Modular Design**: Separation of I/O, computation, and analysis
- **Statistical Rigor**: FSML library for robust distribution fitting
- **Climate Standards**: CF-1.8 compliant NetCDF output
- **Bias Correction**: Multiplicative scaling with ERA5 reference
- **Quality Control**: Comprehensive validation and error handling

## 📈 Data Products

### Drought Indices
- **SPI** (Standardized Precipitation Index): 1, 3, 6, 12-month timescales
- **SPEI** (Standardized Precipitation Evapotranspiration Index): 1, 3, 6, 12-month timescales

### Climate Scenarios
- **ERA5-Land** (1981-2024): Historical baseline
- **SSP1-2.6** (2015-2099): Low emissions (~1.5°C warming)
- **SSP2-4.5** (2015-2099): Moderate emissions (~2.5°C warming)
- **SSP5-8.5** (2015-2099): High emissions (~4.5°C warming)

### Output Structure
```
data/Final_results/Drought_Indices/
├── ERA5/                        # Historical baseline (1981-2024)
│   ├── SPI/
│   └── SPEI/
└── CMIP6/                       # Future projections (2015-2099)
    ├── ssp126/
    ├── ssp245/
    └── ssp585/
```

## 🔬 Methodology Highlights

### Bias Correction
- **Driver-based approach**: Calculate PET from temperature using Hargreaves method
- **Direct PET correction**: Bias-correct both precipitation and PET separately
- **Multiplicative scaling**: Preserve relative climate change signals
- **Reference period**: 1981-2014 overlap between CMIP6 and ERA5

### Statistical Methods
- **Gamma distribution fitting**: Maximum likelihood estimation using FSML
- **Standardization**: Normal inverse transformation for drought indices
- **Quality control**: Range validation, temporal continuity checks
- **Missing data handling**: Spatial and temporal interpolation

### Limitations and Justifications
- **Hargreaves PET**: Temperature-only method justified by data availability
- **Stationarity assumption**: Bias factors constant over time
- **CMIP6 resolution**: ~50 km grid, bias-corrected to 0.1° resolution
- **Model uncertainty**: Addressed through ensemble approach

## 🎯 Results Summary

### Key Findings (Preliminary)
- **Historical trends**: Increasing drought frequency observed in ERA5 data
- **Future projections**: Clear progression of drought severity across emission scenarios
- **Climate signal**: SSP5-8.5 shows extreme drought conditions (+3.09 SPEI consistently)
- **Regional patterns**: Spatial variability in drought response across Somaliland

### Data Quality
- **Coverage**: >95% spatial completeness
- **Validation**: Strong correlation with ERA5 reference (r > 0.7)
- **Physical realism**: SPEI values within expected range [-4, +4]
- **Temporal continuity**: No gaps or discontinuities in time series

## 🤝 Contributing

### Development Environment
```bash
# Set up development environment
git clone https://github.com/Khadar146/drought_model.git
cd drought_model
conda env create -f environment.yml
conda activate Drought
```

### Code Standards
- **Fortran 2008+**: Modern Fortran practices
- **Modular design**: Separate modules for I/O, computation, analysis
- **Documentation**: Comprehensive inline comments and external docs
- **Testing**: Validation against reference datasets

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact

**Khadar Daahir**  
PhD Student, School of Geographical and Earth Sciences  
University of Glasgow  
📧 Email: khadardahir146@gmail.com  
🔗 GitHub: [@Khadar146](https://github.com/Khadar146)

## 🙏 Acknowledgments

- **CMIP6 Modeling Groups**: Climate model outputs
- **ECMWF**: ERA5-Land reanalysis data
- **FSML Library**: Statistical methods implementation
- **University of Glasgow**: Research support and supervision
- **Climate Research Community**: Methodological foundations

## 📊 Project Metrics

| Metric | Value | Notes |
|--------|-------|-------|
| **Lines of Code** | ~3,500 | Fortran source code |
| **Data Volume** | ~4.9 GB | Complete pipeline outputs |
| **Processing Time** | ~45 minutes | Full pipeline on 8-core system |
| **Spatial Points** | 2,016 | 0.1° grid over Somaliland |
| **Temporal Coverage** | 119 years | 1981-2099 analysis period |
| **Scenarios Analyzed** | 4 | ERA5 + 3 SSP scenarios |

---

**Pipeline Version**: 1.0.0  
**Documentation Updated**: August 26, 2025  
**Status**: ✅ Production Ready
