# Somaliland Drought Analysis System

## 🌍 Project Overview

A comprehensive Fortran-based system for analyzing drought patterns in Somaliland, designed to test the hypothesis that **future droughts will become more extreme and frequent due to climate change**.

## 🎯 Research Hypothesis

**"Droughts in Somaliland will become more extreme and frequent in the future based on climate projections"**

This system provides:
- ✅ **Historical drought baseline** (1980-2024) - COMPLETED
- 🔄 **Future drought projections** (2025-2100) - IN DEVELOPMENT
- 📊 **Comparative change analysis** - PLANNED

## 🏗️ Project Structure

```
📁 droughtmodel/
├── 📁 src/                    # Core analysis modules
├── 📁 app/
│   ├── 📁 historical/         # Historical analysis (SPI + EVT)
│   └── 📁 future/            # Future projections (SPI + EVT)
├── 📁 outputs/
│   ├── 📁 historical/         # Historical results
│   ├── 📁 future/            # Future projection results  
│   └── 📁 comparison/        # Change analysis
├── 📁 data/
│   ├── 📁 climate/historical/ # Historical precipitation data
│   └── 📁 climate/future/    # Climate projection data
├── 📁 test/                   # Validation and testing
└── 📁 docs/                   # Documentation
```

## 🔬 Scientific Methods

### Historical Analysis ✅
- **SPI Calculation**: Standardized Precipitation Index using FSML statistical routines
- **EVT Analysis**: Extreme Value Theory with Generalized Pareto Distribution
- **Reference Period**: 1991-2020 (WMO standard)
- **Temporal Resolution**: Monthly data, multiple accumulation periods (1, 3, 6, 12 months)

### Future Projections 🔄
- **Climate Models**: CMIP6 projections (SSP scenarios)
- **Bias Correction**: Statistical adjustment using historical relationships
- **Future SPI**: Using historical gamma parameters for consistency
- **Change Detection**: Quantitative assessment of drought frequency/intensity changes

## 🚀 Quick Start

### Prerequisites
- Fortran compiler (gfortran recommended)
- NetCDF libraries (`brew install netcdf-fortran` on macOS)
- FPM (Fortran Package Manager)

### Build & Run
```bash
# Clone and navigate to project
cd droughtmodel

# Build the system
fpm build --flag "$(nf-config --fflags)" --link-flag "$(nf-config --flibs)"

# Run historical analysis (all SPI timescales)
fpm run HistoricalAnalysis

# Check outputs
ls outputs/historical/spi/
```

## 📊 Current Results

### Historical Analysis (VALIDATED ✅)
- **Dataset**: 540 months (1980-2024)
- **Spatial Coverage**: 66×41 grid points (Somaliland region)
- **SPI Range**: -5.2 to +7.5 (no artificial capping)
- **Statistical Properties**: Mean ≈ 0, appropriate distribution

## 🤝 Author & Contact

**Khadar Daahir Abdisalan**  
Email: khadar@glasgow.ac.uk  
Institution: University of Glasgow

---

**Status**: Historical analysis complete ✅ | Future projections in development 🔄  
**Last Updated**: August 2025
