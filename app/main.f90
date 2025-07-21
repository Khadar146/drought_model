program main
  ! -------------------------------------------------------------------
  ! MSc Drought Modelling – Main Program
  ! Purpose: Executes the drought analysis workflow for Somaliland
  ! using SPI, EVT, and forecasting modules with NetCDF input.
  !
  ! Author: Khadar
  ! Date: June - August 2025
  ! -------------------------------------------------------------------
  use spi_mod         ! Computes SPI from NetCDF data (includes loading + preprocessing)
  ! use evt_mod       ! Applies Extreme Value Theory to SPI results (future implementation)
  ! use forecast_mod  ! Performs drought forecasting (e.g. regression) (future implementation)
  implicit none

  ! Welcome message
  print *, "---------------------------------------------"
  print *, " Running Drought Modelling Pipeline..."
  print *, "---------------------------------------------"

  ! Step 1: Compute SPI from historical and/or future data
  call run_spi_module()

  ! Step 2: Apply EVT to analyse severity and frequency
  ! call run_evt_module()

  ! Step 3: Run drought forecasting using climate predictors
  ! call run_forecast_module()

  ! Completion message
  print *, "---------------------------------------------"
  print *, " Pipeline completed successfully."
  print *, " Outputs saved to /output."
  print *, "---------------------------------------------"

end program main