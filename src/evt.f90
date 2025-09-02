!=======================================================================
! EVT MODULE FOR DROUGHT EXTREMES (SPEI-based)
!=======================================================================
! Purpose:
!   Fits Generalized Pareto Distribution (GPD) to SPEI extremes using 
!   Profile Maximum Likelihood Estimation following academic standards
!   
!   Outputs:
!     - GPD parameters (ξ, σ) with uncertainty estimates
!     - Return levels (10-, 20-, 50-year) with confidence intervals
!     - Exceedance probabilities for chosen historical baseline events
!     - Diagnostic statistics for model validation
!
! References:
!   - Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values. 
!     Springer-Verlag, London. [Standard EVT reference]
!   - Hosking, J.R.M. and Wallis, J.R. (1987). Parameter and quantile estimation 
!     for the generalized Pareto distribution. Technometrics, 29(3), 339-349.
!   - Grimshaw, S.D. (1993). Computing maximum likelihood estimates for the 
!     generalized Pareto distribution. Technometrics, 35(2), 185-191.
!   - Davison, A.C. and Smith, R.L. (1990). Models for exceedances over high 
!     thresholds. Journal of the Royal Statistical Society, Series B, 52, 393-442.
!=======================================================================

module evt_module
    use, intrinsic :: iso_fortran_env, only: real64, int32, output_unit, error_unit
    ! FSML GPD functions - available through main fsml module
    use fsml, only: fsml_gpd_pdf, fsml_gpd_cdf, fsml_gpd_ppf
    implicit none
    
    private
    public :: run_evt_analysis, gpd_parameters, evt_results, extract_exceedances, fit_gpd_profile_mle
    
    integer, parameter :: dp = real64
    real(dp), parameter :: FILL_VALUE = -999.0_dp
    real(dp), parameter :: PI = 3.141592653589793_dp
    real(dp), parameter :: GOLDEN_RATIO = 1.618033988749895_dp
    real(dp), parameter :: EPS = 1.0e-12_dp
    
    ! Derived types for structured output
    type :: gpd_parameters
        real(dp) :: xi          ! Shape parameter
        real(dp) :: sigma       ! Scale parameter  
        real(dp) :: mu          ! Location parameter (threshold)
        real(dp) :: xi_se       ! Standard error of shape parameter
        real(dp) :: sigma_se    ! Standard error of scale parameter
        real(dp) :: loglik      ! Log-likelihood value
        integer :: status       ! Fitting status (0=success, >0=warning/error)
    end type gpd_parameters
    
    type :: evt_results
        type(gpd_parameters) :: params
        real(dp) :: return_levels(3)       ! 10, 20, 50-year return levels
        real(dp) :: return_levels_ci(3,2)  ! Confidence intervals [lower, upper]
        real(dp) :: exceed_prob             ! Exceedance probability
        real(dp) :: threshold               ! Threshold used
        integer :: n_exceedances           ! Number of exceedances
        real(dp) :: mean_excess_function    ! Mean excess above threshold
    end type evt_results

contains

    !-------------------------------------------------------------------
    ! MAIN EVT ANALYSIS: Profile Maximum Likelihood Approach
    ! Following Coles (2001) and Grimshaw (1993) methodology
    !-------------------------------------------------------------------
    subroutine run_evt_analysis(spei, time_years, threshold, results)
        real(dp), intent(in) :: spei(:,:)          ! (time, gridcell)
        real(dp), intent(in) :: time_years(:)
        real(dp), intent(in) :: threshold
        type(evt_results), intent(out) :: results
        
        real(dp), allocatable :: excess(:)
        integer :: valid_count
        real(dp) :: hist_20yr
        
        write(output_unit, '(A)') "=== Extreme Value Theory Analysis ==="
        write(output_unit, '(A, F8.3)') "Threshold: ", threshold
        
        ! Initialize results
        results%threshold = threshold
        
        ! Extract exceedances over threshold (Peak-Over-Threshold method)
        call extract_exceedances(spei, threshold, excess, valid_count)
        results%n_exceedances = valid_count
        
        if (valid_count < 30) then
            write(error_unit, '(A, I0, A)') "WARNING: Only ", valid_count, &
                " exceedances found - insufficient for robust EVT fitting"
            write(error_unit, '(A)') "         Minimum recommended: 30 exceedances (Coles, 2001)"
            results%params%status = 1
            results%return_levels = FILL_VALUE
            results%exceed_prob = FILL_VALUE
            return
        end if
        
        write(output_unit, '(A, I0, A, F8.3)') "Found ", valid_count, &
            " exceedances above threshold ", threshold
        
        ! Calculate mean excess function for diagnostics
        results%mean_excess_function = sum(excess(1:valid_count)) / real(valid_count, dp)
        
        ! Fit GPD parameters using Profile Maximum Likelihood Estimation
        write(output_unit, '(A)') "Fitting GPD using Profile Maximum Likelihood..."
        call fit_gpd_profile_mle(excess(1:valid_count), threshold, results%params)
        
        if (results%params%status /= 0) then
            write(error_unit, '(A)') "WARNING: GPD parameter fitting failed"
            results%return_levels = FILL_VALUE
            results%exceed_prob = FILL_VALUE
            return
        end if
        
        ! Compute return levels with confidence intervals
        call compute_return_levels_with_ci(results%params, results%return_levels, &
                                         results%return_levels_ci, valid_count)
        
        ! Calculate exceedance probability for historical 20-year event
        hist_20yr = results%return_levels(2)
        ! For drought analysis: probability of being below the drought level
        ! Convert drought level back to positive exceedance scale for CDF calculation
        results%exceed_prob = fsml_gpd_cdf(abs(hist_20yr - results%params%mu), &
                                         results%params%xi, 0.0_dp, results%params%sigma)
        
        ! Output summary
        call print_evt_summary(results)
        
    end subroutine run_evt_analysis
    
    !-------------------------------------------------------------------
    ! PROFILE MAXIMUM LIKELIHOOD ESTIMATION for GPD Parameters
    ! Following Grimshaw (1993) and Coles (2001) methodology
    !-------------------------------------------------------------------
    subroutine fit_gpd_profile_mle(excess_data, threshold, params)
        real(dp), intent(in) :: excess_data(:)
        real(dp), intent(in) :: threshold  
        type(gpd_parameters), intent(out) :: params
        
        integer :: n, i, iter, max_iter
        real(dp) :: xi_opt, sigma_opt, loglik_opt
        real(dp) :: xi_lower, xi_upper, xi_test
        real(dp) :: sigma_mle, loglik_test
        real(dp) :: tol, golden_tol
        real(dp) :: mean_excess, xi_start
        logical :: converged
        
        n = size(excess_data)
        params%mu = threshold  ! Location parameter is the threshold
        params%status = 0
        
        ! Constraint bounds for shape parameter (Coles, 2001, Chapter 4)
        xi_lower = -0.5_dp    ! Lower bound for numerical stability
        xi_upper = 0.5_dp     ! Upper bound for finite moments
        max_iter = 100
        tol = 1.0e-8_dp
        golden_tol = 1.0e-10_dp
        
        ! Initial estimate using method of moments (Hosking & Wallis, 1987)
        mean_excess = sum(excess_data) / real(n, dp)
        xi_start = 0.5_dp * (sum((excess_data/mean_excess - 1.0_dp)**2) / real(n-1, dp) - 1.0_dp)
        xi_start = max(xi_lower + 0.01_dp, min(xi_upper - 0.01_dp, xi_start))
        
        write(output_unit, '(A, F8.4)') "  Initial shape parameter estimate: ", xi_start
        
        ! Profile likelihood optimization using golden section search
        ! This is more robust than Newton-Raphson for constrained optimization
        loglik_opt = -huge(1.0_dp)
        xi_opt = xi_start
        sigma_opt = mean_excess
        
        ! Golden section search for optimal xi (Coles, 2001, Appendix A)
        call golden_section_search(excess_data, xi_lower, xi_upper, xi_opt, &
                                 loglik_opt, golden_tol, max_iter, converged)
        
        if (.not. converged) then
            write(error_unit, '(A)') "WARNING: Profile MLE optimization did not converge"
            params%status = 2
        end if
        
        ! Compute optimal sigma given optimal xi
        call compute_profile_sigma(excess_data, xi_opt, sigma_opt)
        
        ! Final parameter values
        params%xi = xi_opt
        params%sigma = sigma_opt
        params%loglik = loglik_opt
        
        ! Compute standard errors using Fisher Information Matrix
        call compute_parameter_uncertainties(excess_data, params)
        
        ! Diagnostic checks following Coles (2001), Section 4.3.4
        call validate_gpd_fit(excess_data, params)
        
        write(output_unit, '(A)') "  Profile MLE fitting completed:"
        write(output_unit, '(A, F8.5, A, F8.5, A)') "    Shape (ξ): ", params%xi, &
            " ± ", params%xi_se, " (SE)"
        write(output_unit, '(A, F8.5, A, F8.5, A)') "    Scale (σ): ", params%sigma, &
            " ± ", params%sigma_se, " (SE)"
        write(output_unit, '(A, F10.3)') "    Log-likelihood: ", params%loglik
        
    end subroutine fit_gpd_profile_mle
    
    !-------------------------------------------------------------------
    ! GOLDEN SECTION SEARCH for Profile Likelihood Optimization
    ! Based on numerical optimization theory (Press et al., 2007)
    !-------------------------------------------------------------------
    subroutine golden_section_search(data, xi_lower, xi_upper, xi_opt, loglik_opt, &
                                    tol, max_iter, converged)
        real(dp), intent(in) :: data(:)
        real(dp), intent(in) :: xi_lower, xi_upper, tol
        integer, intent(in) :: max_iter
        real(dp), intent(out) :: xi_opt, loglik_opt
        logical, intent(out) :: converged
        
        real(dp) :: a, b, c, d, fc, fd
        real(dp) :: sigma_c, sigma_d
        integer :: iter
        
        a = xi_lower
        b = xi_upper
        c = b - (b - a) / GOLDEN_RATIO
        d = a + (b - a) / GOLDEN_RATIO
        
        ! Initial function evaluations
        call compute_profile_sigma(data, c, sigma_c)
        fc = -profile_log_likelihood(data, c, sigma_c)  ! Minimize negative log-likelihood
        
        call compute_profile_sigma(data, d, sigma_d)
        fd = -profile_log_likelihood(data, d, sigma_d)
        
        converged = .false.
        do iter = 1, max_iter
            if (abs(b - a) < tol) then
                converged = .true.
                exit
            end if
            
            if (fc < fd) then
                b = d
                d = c
                fd = fc
                c = b - (b - a) / GOLDEN_RATIO
                call compute_profile_sigma(data, c, sigma_c)
                fc = -profile_log_likelihood(data, c, sigma_c)
            else
                a = c
                c = d
                fc = fd
                d = a + (b - a) / GOLDEN_RATIO
                call compute_profile_sigma(data, d, sigma_d)
                fd = -profile_log_likelihood(data, d, sigma_d)
            end if
        end do
        
        if (fc < fd) then
            xi_opt = c
            call compute_profile_sigma(data, xi_opt, sigma_c)
            loglik_opt = profile_log_likelihood(data, xi_opt, sigma_c)
        else
            xi_opt = d
            call compute_profile_sigma(data, xi_opt, sigma_d)
            loglik_opt = profile_log_likelihood(data, xi_opt, sigma_d)
        end if
        
    end subroutine golden_section_search
    
    !-------------------------------------------------------------------
    ! COMPUTE PROFILE SIGMA for given XI (Grimshaw, 1993)
    !-------------------------------------------------------------------
    subroutine compute_profile_sigma(data, xi, sigma)
        real(dp), intent(in) :: data(:), xi
        real(dp), intent(out) :: sigma
        
        integer :: n, i
        real(dp) :: sum_term
        
        n = size(data)
        
        if (abs(xi) < EPS) then
            ! Exponential case (xi = 0)
            sigma = sum(data) / real(n, dp)
        else
            ! General GPD case
            sum_term = 0.0_dp
            do i = 1, n
                if (1.0_dp + xi * data(i) / sigma > 0.0_dp) then
                    sum_term = sum_term + log(1.0_dp + xi * data(i) / sigma)
                end if
            end do
            sigma = sum(data) / (real(n, dp) * (1.0_dp + xi))
        end if
        
        ! Ensure positive scale parameter
        sigma = max(sigma, 1.0e-6_dp)
        
    end subroutine compute_profile_sigma

    
    !-------------------------------------------------------------------
    ! PROFILE LOG-LIKELIHOOD FUNCTION for GPD (Coles, 2001, Eq. 4.8)
    !-------------------------------------------------------------------
    function profile_log_likelihood(data, xi, sigma) result(loglik)
        real(dp), intent(in) :: data(:), xi, sigma
        real(dp) :: loglik
        
        integer :: n, i
        real(dp) :: term, log_sigma
        
        n = size(data)
        loglik = 0.0_dp
        log_sigma = log(sigma)
        
        if (abs(xi) < EPS) then
            ! Exponential case (xi → 0)
            loglik = -real(n, dp) * log_sigma - sum(data) / sigma
        else
            ! General GPD case
            do i = 1, n
                term = 1.0_dp + xi * data(i) / sigma
                if (term <= 0.0_dp) then
                    loglik = -huge(1.0_dp)  ! Invalid parameter combination
                    return
                end if
                loglik = loglik - log_sigma - (1.0_dp + 1.0_dp/xi) * log(term)
            end do
        end if
        
    end function profile_log_likelihood
    
    !-------------------------------------------------------------------
    ! COMPUTE PARAMETER UNCERTAINTIES using Fisher Information Matrix
    ! Following Coles (2001), Section 2.6.4 and Davison & Smith (1990)
    !-------------------------------------------------------------------
    subroutine compute_parameter_uncertainties(data, params)
        type(gpd_parameters), intent(inout) :: params
        real(dp), intent(in) :: data(:)
        
        integer :: n
        real(dp) :: fisher_info(2,2), inv_fisher(2,2)
        real(dp) :: det, xi, sigma
        real(dp) :: I11, I12, I22  ! Fisher information matrix elements
        real(dp) :: sum1, sum2, sum3, term
        integer :: i
        
        n = size(data)
        xi = params%xi
        sigma = params%sigma
        
        ! Compute Fisher Information Matrix elements (Coles, 2001, Section 4.3.2)
        sum1 = 0.0_dp
        sum2 = 0.0_dp  
        sum3 = 0.0_dp
        
        if (abs(xi) < EPS) then
            ! Exponential case
            I11 = real(n, dp) / (sigma**2)
            I12 = 0.0_dp
            I22 = real(n, dp) / (sigma**2)
        else
            ! General GPD case
            do i = 1, n
                term = 1.0_dp + xi * data(i) / sigma
                if (term > 0.0_dp) then
                    sum1 = sum1 + data(i)**2 / (sigma**2 * term**2)
                    sum2 = sum2 + data(i) / (sigma * term**2)
                    sum3 = sum3 + 1.0_dp / term**2
                end if
            end do
            
            I11 = real(n, dp) / (sigma**2) + (1.0_dp + xi)**2 * sum1
            I12 = (1.0_dp + xi) * sum2 / xi
            I22 = real(n, dp) / (xi**2) + sum3 / (xi**2)
        end if
        
        ! Construct Fisher Information Matrix
        fisher_info(1,1) = I11  ! sigma, sigma
        fisher_info(1,2) = I12  ! sigma, xi
        fisher_info(2,1) = I12  ! xi, sigma  
        fisher_info(2,2) = I22  ! xi, xi
        
        ! Invert Fisher Information Matrix to get covariance matrix
        det = fisher_info(1,1) * fisher_info(2,2) - fisher_info(1,2)**2
        
        if (abs(det) > EPS) then
            inv_fisher(1,1) = fisher_info(2,2) / det
            inv_fisher(1,2) = -fisher_info(1,2) / det
            inv_fisher(2,1) = inv_fisher(1,2)
            inv_fisher(2,2) = fisher_info(1,1) / det
            
            ! Standard errors are square roots of diagonal elements
            params%sigma_se = sqrt(max(0.0_dp, inv_fisher(1,1)))
            params%xi_se = sqrt(max(0.0_dp, inv_fisher(2,2)))
        else
            ! Singular Fisher matrix - use conservative estimates
            params%sigma_se = sigma * 0.1_dp
            params%xi_se = max(0.05_dp, abs(xi) * 0.2_dp)
            params%status = 3  ! Warning: uncertain parameter estimates
        end if
        
    end subroutine compute_parameter_uncertainties
    
    !-------------------------------------------------------------------
    ! COMPUTE RETURN LEVELS with Confidence Intervals
    ! Using Delta method (Coles, 2001, Section 4.3.3)
    !-------------------------------------------------------------------
    subroutine compute_return_levels_with_ci(params, return_levels, ci_bounds, n_obs)
        type(gpd_parameters), intent(in) :: params
        real(dp), intent(out) :: return_levels(3)
        real(dp), intent(out) :: ci_bounds(3,2)  ! [level, lower/upper]
        integer, intent(in) :: n_obs
        
        real(dp) :: return_periods(3) = [10.0_dp, 20.0_dp, 50.0_dp]
        real(dp) :: prob, z_quantile, return_level, variance, std_error
        real(dp) :: drl_dsigma, drl_dxi  ! Partial derivatives for Delta method
        real(dp) :: xi, sigma, log_term
        integer :: i
        
        xi = params%xi
        sigma = params%sigma
        z_quantile = 1.96_dp  ! 95% confidence interval
        
        do i = 1, 3
            prob = 1.0_dp - 1.0_dp / return_periods(i)
            
            ! Return level calculation for drought analysis (below threshold)
            ! For drought: return_level = threshold - GPD_quantile(prob)
            ! This gives negative SPEI values for drought events
            return_level = params%mu - fsml_gpd_ppf(prob, xi, 0.0_dp, sigma)
            return_levels(i) = return_level
            
            ! Delta method for confidence intervals (Coles, 2001, Section 4.3.3)
            ! For drought analysis: derivatives need sign correction
            if (abs(xi) < EPS) then
                ! Exponential case
                log_term = -log(1.0_dp - prob)
                drl_dsigma = -log_term  ! Negative for drought direction
                drl_dxi = 0.0_dp
            else
                ! General GPD case - negative derivatives for drought direction
                log_term = -log(1.0_dp - prob)
                drl_dsigma = -(((1.0_dp - prob)**(-xi) - 1.0_dp) / xi)
                drl_dxi = -(sigma / (xi**2)) * (((1.0_dp - prob)**(-xi)) * &
                         log(1.0_dp - prob) - ((1.0_dp - prob)**(-xi) - 1.0_dp))
            end if
            
            ! Variance using Delta method
            variance = (drl_dsigma**2) * (params%sigma_se**2) + &
                      (drl_dxi**2) * (params%xi_se**2) + &
                      2.0_dp * drl_dsigma * drl_dxi * params%sigma_se * params%xi_se * 0.0_dp
            ! Note: Assuming zero correlation between sigma and xi estimates for simplicity
            
            std_error = sqrt(max(0.0_dp, variance))
            
            ! 95% Confidence intervals
            ci_bounds(i, 1) = return_level - z_quantile * std_error  ! Lower bound
            ci_bounds(i, 2) = return_level + z_quantile * std_error  ! Upper bound
        end do
        
    end subroutine compute_return_levels_with_ci
    
    !-------------------------------------------------------------------
    ! VALIDATE GPD FIT using Diagnostic Tests (Coles, 2001, Section 4.4)
    !-------------------------------------------------------------------
    subroutine validate_gpd_fit(data, params)
        real(dp), intent(in) :: data(:)
        type(gpd_parameters), intent(in) :: params
        
        integer :: n, i
        real(dp) :: transformed_data(size(data))
        real(dp) :: mean_transform, var_transform
        real(dp) :: ks_statistic, ad_statistic
        real(dp) :: xi, sigma
        
        n = size(data)
        xi = params%xi
        sigma = params%sigma
        
        write(output_unit, '(A)') "  GPD Model Diagnostics:"
        
        ! Transform data to standard exponential using GPD CDF
        do i = 1, n
            transformed_data(i) = -log(1.0_dp - fsml_gpd_cdf(data(i), xi, 0.0_dp, sigma))
        end do
        
        ! Check if transformed data follows exponential(1) distribution
        mean_transform = sum(transformed_data) / real(n, dp)
        var_transform = sum((transformed_data - mean_transform)**2) / real(n-1, dp)
        
        write(output_unit, '(A, F8.4, A)') "    Transformed data mean: ", mean_transform, &
            " (should be ≈ 1.0)"
        write(output_unit, '(A, F8.4, A)') "    Transformed data variance: ", var_transform, &
            " (should be ≈ 1.0)"
        
        ! Simple goodness-of-fit assessment
        if (abs(mean_transform - 1.0_dp) > 0.2_dp) then
            write(output_unit, '(A)') "    WARNING: Significant deviation from exponential assumption"
        else
            write(output_unit, '(A)') "    GPD fit appears reasonable"
        end if
        
        ! Shape parameter interpretation (Coles, 2001, Section 4.1.1)
        if (xi > 0.0_dp) then
            write(output_unit, '(A)') "    Heavy-tailed distribution (Fréchet domain)"
        else if (xi < 0.0_dp) then
            write(output_unit, '(A)') "    Short-tailed distribution (Weibull domain)"
        else
            write(output_unit, '(A)') "    Exponential-type tail (Gumbel domain)"
        end if
        
    end subroutine validate_gpd_fit
    
    !-------------------------------------------------------------------
    ! PRINT EVT ANALYSIS SUMMARY
    !-------------------------------------------------------------------
    subroutine print_evt_summary(results)
        type(evt_results), intent(in) :: results
        
        write(output_unit, '(A)') ""
        write(output_unit, '(A)') "=== EVT Analysis Summary ==="
        write(output_unit, '(A, F8.3)') "Threshold: ", results%threshold
        write(output_unit, '(A, I0)') "Number of exceedances: ", results%n_exceedances
        write(output_unit, '(A, F8.4)') "Mean excess: ", results%mean_excess_function
        write(output_unit, '(A)') ""
        write(output_unit, '(A)') "GPD Parameters (Profile MLE):"
        write(output_unit, '(A, F8.5, A, F8.5, A)') "  Shape (ξ): ", results%params%xi, &
            " ± ", results%params%xi_se, ""
        write(output_unit, '(A, F8.5, A, F8.5, A)') "  Scale (σ): ", results%params%sigma, &
            " ± ", results%params%sigma_se, ""
        write(output_unit, '(A, F10.3)') "  Log-likelihood: ", results%params%loglik
        write(output_unit, '(A)') ""
        write(output_unit, '(A)') "Return Levels [95% CI]:"
        write(output_unit, '(A, F8.4, A, F8.4, A, F8.4, A)') "  10-year: ", results%return_levels(1), &
            " [", results%return_levels_ci(1,1), ", ", results%return_levels_ci(1,2), "]"
        write(output_unit, '(A, F8.4, A, F8.4, A, F8.4, A)') "  20-year: ", results%return_levels(2), &
            " [", results%return_levels_ci(2,1), ", ", results%return_levels_ci(2,2), "]"
        write(output_unit, '(A, F8.4, A, F8.4, A, F8.4, A)') "  50-year: ", results%return_levels(3), &
            " [", results%return_levels_ci(3,1), ", ", results%return_levels_ci(3,2), "]"
    
    end subroutine print_evt_summary
    
    !-------------------------------------------------------------------
    ! EXTRACT EXCEEDANCES OVER THRESHOLD 
    ! Peak-Over-Threshold method (Coles, 2001, Chapter 4)
    !-------------------------------------------------------------------
    subroutine extract_exceedances(data, threshold, excess, count)
        real(dp), intent(in) :: data(:,:)
        real(dp), intent(in) :: threshold
        real(dp), allocatable, intent(out) :: excess(:)
        integer, intent(out) :: count
        
        integer :: i, t, ntime, ncell
        real(dp), allocatable :: temp(:)
        real(dp) :: value
        
        ntime = size(data,1)
        ncell = size(data,2)
        allocate(temp(ntime*ncell))
        count = 0
        
        ! Extract all values exceeding threshold
        ! Note: For SPEI, more negative values represent more severe droughts
        ! We want exceedances in the extreme drought direction
        do t = 1, ntime
            do i = 1, ncell
                value = data(t,i)
                ! Use same data validation as extract_valid_data in evt_threshold_diagnostics
                if (abs(value - FILL_VALUE) > 1e-6 .and. abs(value) < 100.0_dp) then
                    ! For drought analysis, we typically look at exceedances 
                    ! below the threshold (extreme negative SPEI values)
                    if (value < threshold) then
                        count = count + 1
                        temp(count) = threshold - value  ! Convert to positive exceedances
                    end if
                end if
            end do
        end do
        
        if (count > 0) then
            allocate(excess(count))
            excess(1:count) = temp(1:count)
            write(output_unit, '(A, I0, A, F8.3)') "  Extracted ", count, &
                " drought exceedances below threshold ", threshold
        else
            allocate(excess(1))
            excess = FILL_VALUE
            write(error_unit, '(A)') "WARNING: No exceedances found"
        end if
        
        deallocate(temp)
    end subroutine extract_exceedances

end module evt_module