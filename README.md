<<<<<<< HEAD
# drought_model
Somaliland. Computes SPI from NetCDF precipitation, includes preprocessing, historical validation, and EVT-based severity analysis. Designed for academic use in climate risk and drought forecasting.
=======
# droughtmodel

A modular Fortran project to compute SPI/SPEI drought indices and apply robust feature selection for drought forecasting in data-scarce regions.

## Author
Khadar Daahir Abdisalan 
Contact: khadardahir146@gmail.com

## Structure
- `src/`: Main modules (SPI, SPEI, filtering, FS)
- `app/`: Execution logic
- `data/`, `output/`: I/O folders
- `docs/`: Documentation (for use with Ford)


## Main Modules
- `spi_mod.f90`: SPI computation
- `evt.f90`: EVT and GEV fitting
- `main.f90`: Main program


## Build & Run

```bash
fpm build
fpm run
>>>>>>> 99852a9 (Initial commit: Fortran-based SPI and EVT drought model)
