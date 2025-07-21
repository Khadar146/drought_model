program drought_model
  use extract               ! Stage 1 – Extraction
  use preprocess_mod            ! Stage 2 – Cleaning
  use spi_mod                   ! Stage 3a – SPI
  use spei_mod                  ! Stage 3b – SPEI
  use feature_selection_mod     ! Stage 4 – Feature Selection
  use regression_mod            ! Stage 5 – Modelling
  use evt_mod                   ! Stage 6 – EVT
  use validation_mod            ! Stage 7 – Validation
  implicit none

  !───────────────
  ! Stage 1: NetCDF to CSV
  !───────────────
  call run_extraction()

  !───────────────
  ! Stage 2: Data Cleaning + Preprocessing
  !───────────────
  call clean_data()

  !───────────────
  ! Stage 3: SPI and SPEI Index Calculation
  !───────────────
  call compute_spi()
  call compute_spei()

  !───────────────
  ! Stage 4: Feature Selection
  !───────────────
  call run_selection()

  !───────────────
  ! Stage 5: Regression Model Training (OLS/Ridge)
  !───────────────
  call train_model()

  !───────────────
  ! Stage 6: Extreme Value Theory Modelling
  !───────────────
  call run_evt()

  !───────────────
  ! Stage 7: Validation + Cross Checks
  !───────────────
  call validate_against_historical()

  print *, "✅ All stages completed successfully."
end program drought_model