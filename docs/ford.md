---
project: DroughtModel
version: 1.0.0
license: MIT
summary: >
  A comprehensive Fortran-based drought analysis pipeline for Somaliland using 
  historical precipitation data (1980-2024). Implements Standardized Precipitation 
  Index (SPI) calculations across multiple time scales with robust NetCDF I/O, 
  data preprocessing, and extreme value analysis for climate research applications.

author: Khadar Daahir
email: khadar@glasgow.ac.uk
github: https://github.com/username/droughtmodel
institution: University of Glasgow, Climate Dynamics Lab
date: 2025-08-09

# Documentation Configuration
src_dir: src
output_dir: docs/html
exclude_dir: build
exclude_dir: unused
exclude: fpm.toml

# Display Options
display: public
sort: alpha
graph: true
search: true
mathjax: true
warn: true
coloured_edges: true

# File Extensions
extensions: f90 F90 f95 F95 f03 F03 f08 F08

# Project Information
project_description: >
  ## Project Overview
  
  The Somaliland Drought Analysis Model is a scientific computing pipeline designed 
  for comprehensive drought monitoring and analysis in the Horn of Africa region. 
  The model processes historical precipitation data to generate multi-scale drought 
  indices and provides statistical analysis of drought events.
  
  ## Key Features
  
  - **Multi-scale SPI Analysis**: 1, 3, 6, and 12-month SPI calculations
  - **Robust NetCDF I/O**: Handles ERA5 reanalysis data with proper coordinate transformations
  - **Data Quality Assurance**: Comprehensive validation and error checking
  - **Extreme Value Analysis**: Statistical characterization of drought events
  - **Climate Research Focus**: Designed for scientific applications and publications
  
  ## Technical Architecture
  
  ```
  NetCDF Input → IO Module → PREP Module → SPI Module → Analysis & Output
  ```
  
  ### Module Description
  
  - **IO Module**: NetCDF reading, coordinate handling, missing value processing
  - **PREP Module**: Unit conversion, temporal aggregation, quality control
  - **SPI Module**: Standardized Precipitation Index calculations
  - **Analysis Module**: Drought event detection and statistical analysis
  
  ## Data Requirements
  
  - **Input**: Monthly precipitation data in NetCDF format
  - **Coverage**: Somaliland region (42-48.5°E, 8.5-12.5°N)
  - **Resolution**: 0.1° spatial grid (~11km)
  - **Temporal**: 1980-2024 (540 months)
  
  ## Usage
  
  ```bash
  # Compile with FPM
  fpm build
  
  # Run historical analysis
  fpm run main_historical
  
  # Validate IO module
  ./validate_io.sh
  ```

# External Links
github: https://github.com/username/droughtmodel
website: https://climate-dynamics.glasgow.ac.uk

# Documentation Pages
page_dir: docs/pages
media_dir: docs/media

# Code Documentation Standards
docmark: !
predocmark: >
creation_date: 2025-08-09
---
