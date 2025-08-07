!> @file prep_mod.f90
!! @brief Preprocesses precipitation data for drought analysis.
!! @details Converts from m/month to mm/month,
!!          fills missing values using monthly climatology,
!!          and optionally aggregates to seasonal (3-month)
!!          and annual (12-month) scales.
!!
!! @author Khadar
!! @affiliation University of Glasgow, Climate Dynamics Lab
!! @date 2025

module prep_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use iso_fortran_env, only: wp => real64
  implicit none
  private
  public :: preprocess_precip

contains

  subroutine preprocess_precip(precip_in, precip_mm, nlon, nlat, ntime, agg3, agg12)
    real(wp), intent(in)  :: precip_in(nlon, nlat, ntime)
    real(wp), intent(out) :: precip_mm(nlon, nlat, ntime)
    integer,  intent(in)  :: nlon, nlat, ntime
    real(wp), intent(out), optional :: agg3(nlon, nlat, ntime)
    real(wp), intent(out), optional :: agg12(nlon, nlat, ntime)

    integer :: i, j, t, m
    real(wp), allocatable :: clim_mean(:,:,:)
    real(wp) :: sum_val
    integer  :: count_val

    ! --- 1. Unit conversion (m/month to mm/month) ---
    precip_mm = precip_in * 1000.0_wp

    ! --- 2. Allocate monthly climatology array ---
    allocate(clim_mean(nlon, nlat, 12))
    clim_mean = 0.0_wp

    ! --- 3. Compute monthly climatology (mean over years, ignoring NaNs) ---
    do m = 1, 12
      do i = 1, nlon
        do j = 1, nlat
          sum_val = 0.0_wp
          count_val = 0
          do t = m, ntime, 12
            if (.not. ieee_is_nan(precip_mm(i,j,t))) then
              sum_val = sum_val + precip_mm(i,j,t)
              count_val = count_val + 1
            end if
          end do
          if (count_val > 0) then
            clim_mean(i,j,m) = sum_val / count_val
          else
            clim_mean(i,j,m) = ieee_value(0.0_wp, ieee_quiet_nan)
          end if
        end do
      end do
    end do

    ! --- 4. Fill missing values using climatology ---
    do t = 1, ntime
      m = mod(t - 1, 12) + 1
      do i = 1, nlon
        do j = 1, nlat
          if (ieee_is_nan(precip_mm(i,j,t))) then
            precip_mm(i,j,t) = clim_mean(i,j,m)
          end if
        end do
      end do
    end do

    ! --- 5. Compute 3-month seasonal aggregates ---
    if (present(agg3)) then
      agg3 = ieee_value(0.0_wp, ieee_quiet_nan)
      do t = 3, ntime
        do i = 1, nlon
          do j = 1, nlat
            if (all(.not. ieee_is_nan(precip_mm(i,j,t-2:t)))) then
              agg3(i,j,t) = sum(precip_mm(i,j,t-2:t))
            end if
          end do
        end do
      end do
    end if

    ! --- 6. Compute 12-month annual aggregates ---
    if (present(agg12)) then
      agg12 = ieee_value(0.0_wp, ieee_quiet_nan)
      do t = 12, ntime
        do i = 1, nlon
          do j = 1, nlat
            if (all(.not. ieee_is_nan(precip_mm(i,j,t-11:t)))) then
              agg12(i,j,t) = sum(precip_mm(i,j,t-11:t))
            end if
          end do
        end do
      end do
    end if

    ! --- 7. Clean up ---
    deallocate(clim_mean)

    print *, "> Preprocessing complete: units converted, gaps filled, aggregation done."

  end subroutine preprocess_precip

end module prep_mod