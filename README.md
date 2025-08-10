# Somaliland Drought Analysis

Fortran-based drought monitoring system for Somaliland using Standardized Precipitation Index (SPI) analysis.

## Requirements
- Fortran compiler (gfortran)
- NetCDF-Fortran library
- Fortran Package Manager (fpm)

## Build
```bash
chmod +x build.sh
./build.sh
```

## Usage
```bash
fpm run HistoricalSPI
```

## Structure
```
app/                    # Main application
src/                    # Core modules (io_mod, prep_mod, spi_mod)
data/                   # Input data
outputs/                # Results
```

## Outputs
- NetCDF files with SPI values (1, 3, 6, 12 months)
- Statistical summaries (CSV)

## Study Area
Somaliland region using ERA5-Land precipitation data (0.1° resolution, 1980-2024)

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
