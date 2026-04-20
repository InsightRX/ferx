#' Fit a nonlinear mixed effects model
#'
#' Fits a NLME model using FOCE or FOCEI estimation with a Rust backend
#' powered by automatic differentiation (Enzyme).
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV file (ID, TIME, DV, EVID, AMT, CMT, ...)
#' @param method Estimation method(s). Either a single string or a character
#'   vector of methods to run in sequence. Each stage is seeded with the
#'   previous stage's converged parameters, and only the final stage produces
#'   the reported covariance/diagnostics. Supported methods: \code{"foce"},
#'   \code{"focei"}, \code{"saem"}, \code{"gn"} (Gauss-Newton / BHHH), or
#'   \code{"gn_hybrid"} (Gauss-Newton followed by a FOCEI polish step).
#'   Example chain: \code{c("saem", "focei")}.
#' @param maxiter Maximum number of outer optimization iterations
#' @param covariance Logical; compute the covariance step for standard errors
#' @param verbose Logical; print progress during estimation
#' @param bloq_method Handling of observations below the lower limit of
#'   quantification. \code{NULL} (default) keeps whatever the model file
#'   specified; \code{"m3"} enables Beal's M3 likelihood (requires a
#'   \code{CENS} column in the data, with \code{DV} carrying the LLOQ value
#'   on \code{CENS=1} rows); \code{"drop"} disables M3.
#' @param optimizer Outer optimizer used during estimation. One of
#'   \code{"slsqp"} (default), \code{"lbfgs"}, \code{"mma"},
#'   \code{"bobyqa"}, or \code{"trust_region"}. \code{"slsqp"} and
#'   \code{"lbfgs"} use gradient information from Enzyme autodiff and are
#'   fastest on smooth surfaces. \code{"bobyqa"} is derivative-free and more
#'   robust on discontinuous or noisy objectives. \code{"trust_region"} uses a
#'   second-order trust-region method. \code{"mma"} is a gradient-based
#'   method-of-moving-asymptotes approach.
#' @param inner_maxiter Maximum number of inner (individual) optimization
#'   iterations per outer step. Default is \code{200}. Reduce for speed at the
#'   cost of inner convergence accuracy; increase if individual fits are not
#'   converging.
#' @param inner_tol Convergence tolerance for the inner (individual)
#'   optimizer. Default is \code{1e-6}. Looser values (e.g. \code{1e-4})
#'   speed up each outer iteration; tighter values improve accuracy.
#' @param mu_referencing Logical. If \code{TRUE} (default), automatically
#'   detects mu-referencing from the model structure for faster and more
#'   accurate convergence. Applies to all estimation methods. Set to
#'   \code{FALSE} to disable for comparison purposes. Detection works
#'   automatically for standard parameterizations such as
#'   \code{PARAM = THETA * exp(ETA)}; unusual parameterizations fall back
#'   silently to zero-centred ETA initialisation with no error. No changes
#'   to the \code{.ferx} model file are needed. Check \code{fit$warnings}
#'   to see which ETAs were detected.
#' @param steihaug_max_iters Maximum conjugate-gradient iterations for the
#'   Steihaug subproblem solver, used only when \code{optimizer = "trust_region"}.
#'   Default is \code{50}, which covers most population PK models (1- and
#'   2-compartment with covariates). Increase for complex models with more than
#'   50 parameters; the theoretical maximum needed is \code{n_params}.
#'
#' @return A list with components:
#'   \item{converged}{Logical; did the optimizer converge}
#'   \item{method}{Estimation method used}
#'   \item{n_iterations}{Number of outer iterations run}
#'   \item{ofv}{Objective function value (-2 log-likelihood)}
#'   \item{aic}{Akaike Information Criterion}
#'   \item{bic}{Bayesian Information Criterion}
#'   \item{theta}{Named numeric vector of fixed effect estimates}
#'   \item{omega}{Between-subject variability covariance matrix}
#'   \item{sigma}{Residual error parameter estimates}
#'   \item{se_theta}{Standard errors for theta (NULL if covariance step failed)}
#'   \item{se_omega}{Standard errors for omega diagonal}
#'   \item{se_sigma}{Standard errors for sigma}
#'   \item{sdtab}{Data frame with ID, TIME, DV, PRED, IPRED, CWRES, IWRES, ETA1..n}
#'   \item{warnings}{Character vector of warnings}
#'   \item{sir_ess}{SIR effective sample size (NULL if SIR not run)}
#'   \item{sir_ci_theta, sir_ci_omega, sir_ci_sigma}{SIR 95\% CI matrices
#'     with columns \code{lower} and \code{upper} (NULL if SIR not run)}
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # Default — mu-referencing on
#' fit <- ferx_fit(ex$model, ex$data)
#'
#' # Compare with mu-referencing off
#' fit_no_mu <- ferx_fit(ex$model, ex$data, mu_referencing = FALSE)
#'
#' # Check which ETAs were detected
#' fit$warnings
#'
#' # Default SLSQP optimizer
#' fit <- ferx_fit(ex$model, ex$data)
#'
#' # Derivative-free BOBYQA — more robust on difficult surfaces
#' fit_bobyqa <- ferx_fit(ex$model, ex$data, optimizer = "bobyqa")
#'
#' # Second-order trust region
#' fit_tr <- ferx_fit(ex$model, ex$data, optimizer = "trust_region")
#'
#' # Trust region with more CG iterations for a complex model (25+ parameters)
#' fit_tr2 <- ferx_fit(ex$model, ex$data, optimizer = "trust_region",
#'                     steihaug_max_iters = 100L)
#'
#' # Fine-tune inner loop speed
#' fit_fast <- ferx_fit(ex$model, ex$data,
#'                      optimizer = "bobyqa",
#'                      inner_maxiter = 100,
#'                      inner_tol = 1e-6)
#'
#' # Chain SAEM to FOCEI (SAEM explores, FOCEI polishes):
#' result <- ferx_fit("warfarin.ferx", "warfarin.csv",
#'                    method = c("saem", "focei"))
#'
#' # Likelihood-based BLOQ handling (M3):
#' bloq <- ferx_example("warfarin_bloq")
#' result <- ferx_fit(bloq$model, bloq$data, method = "focei", bloq_method = "m3")
#' }
#'
#' @export
ferx_fit <- function(model, data,
                     method = "foce",
                     maxiter = 500L,
                     covariance = TRUE,
                     verbose = TRUE,
                     bloq_method = NULL,
                     optimizer = "slsqp",
                     inner_maxiter = 200L,
                     inner_tol = 1e-6,
                     mu_referencing = TRUE,
                     steihaug_max_iters = 50L) {
  stopifnot(file.exists(model), file.exists(data))
  if (!is.logical(mu_referencing) || length(mu_referencing) != 1L || is.na(mu_referencing)) {
    stop("`mu_referencing` must be TRUE or FALSE")
  }
  if (!is.character(method) || length(method) == 0L) {
    stop("`method` must be a non-empty character vector")
  }
  method <- vapply(
    method,
    function(m) {
      match.arg(
        tolower(gsub("[^a-z0-9]", "_", m)),
        c("foce", "focei", "saem", "gn", "gn_hybrid")
      )
    },
    character(1L),
    USE.NAMES = FALSE
  )
  if (is.null(bloq_method)) {
    bloq_arg <- ""
  } else {
    bloq_arg <- match.arg(tolower(bloq_method), c("drop", "m3"))
  }
  optimizer <- match.arg(
    tolower(optimizer),
    c("slsqp", "lbfgs", "mma", "bobyqa", "trust_region")
  )

  raw <- ferx_rust_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    method = method,
    maxiter = as.integer(maxiter),
    covariance = covariance,
    verbose = verbose,
    bloq_method = bloq_arg,
    optimizer = optimizer,
    inner_maxiter = as.integer(inner_maxiter),
    inner_tol = as.double(inner_tol),
    mu_referencing = mu_referencing,
    steihaug_max_iters = as.integer(steihaug_max_iters)
  )

  if (length(raw) == 0) {
    stop("Model fitting failed. Check model and data files.")
  }

  # Structure the result
  result <- raw

  # Name the theta vector
  if (length(result$theta) > 0 && length(result$theta_names) > 0) {
    names(result$theta) <- result$theta_names
  }

  # Reshape omega into a matrix
  if (length(result$omega) > 0 && !is.null(result$omega_dim)) {
    d <- result$omega_dim
    result$omega <- matrix(result$omega, nrow = d, ncol = d)
  }

  # Name SE vectors
  if (length(result$se_theta) > 0 && length(result$theta_names) > 0) {
    names(result$se_theta) <- result$theta_names
  }
  if (length(result$se_theta) == 0) result$se_theta <- NULL
  if (length(result$se_omega) == 0) result$se_omega <- NULL
  if (length(result$se_sigma) == 0) result$se_sigma <- NULL

  # SIR: NaN ess => not computed; flat [lo, hi, lo, hi, ...] => (n, 2) matrix
  if (is.null(result$sir_ess) || !is.finite(result$sir_ess)) {
    result$sir_ess <- NULL
  }
  reshape_ci <- function(v, row_names = NULL) {
    if (length(v) == 0) return(NULL)
    m <- matrix(v, ncol = 2, byrow = TRUE,
                dimnames = list(row_names, c("lower", "upper")))
    m
  }
  result$sir_ci_theta <- reshape_ci(result$sir_ci_theta, result$theta_names)
  n_eta <- if (is.null(dim(result$omega))) NULL else nrow(result$omega)
  eta_names <- if (!is.null(n_eta)) paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")") else NULL
  result$sir_ci_omega <- reshape_ci(result$sir_ci_omega, eta_names)
  sig_names <- paste0("SIGMA(", seq_along(result$sigma), ")")
  result$sir_ci_sigma <- reshape_ci(result$sir_ci_sigma, sig_names)

  # Clean up internal fields
  result$theta_names <- NULL
  result$omega_dim <- NULL

  # Print mu-referencing detections as informational messages
  mu_ref_warnings <- grep("mu-ref", result$warnings, value = TRUE)
  for (w in mu_ref_warnings) {
    eta_names <- sub("^mu-ref:\\s*", "", w)
    message("Mu-referencing detected for: ", eta_names)
  }

  class(result) <- "ferx_fit"
  result
}

#' @export
print.ferx_fit <- function(x, ...) {
  bar <- strrep("=", 60)
  cat(bar, "\n", sep = "")
  cat("NONLINEAR MIXED EFFECTS MODEL ESTIMATION\n")
  cat(bar, "\n\n", sep = "")

  cat("Converged: ", if (isTRUE(x$converged)) "YES" else "NO", "\n", sep = "")
  if (!is.null(x$method_chain) && length(x$method_chain) > 1) {
    cat("Estimation chain:  ", paste(x$method_chain, collapse = " -> "), "\n", sep = "")
  } else {
    cat("Estimation method: ", x$method, "\n", sep = "")
  }
  if (!is.null(x$n_iterations)) {
    cat("Iterations: ", x$n_iterations, "\n", sep = "")
  }

  cat("\n--- Objective Function ---\n")
  cat(sprintf("OFV:  %.4f\n", x$ofv))
  cat(sprintf("AIC:  %.4f\n", x$aic))
  cat(sprintf("BIC:  %.4f\n", x$bic))

  n_par <- if (is.null(x$n_parameters)) NA_integer_ else x$n_parameters
  cat(sprintf(
    "\nSubjects: %d  Observations: %d  Parameters: %s\n",
    x$n_subjects, x$n_obs,
    if (is.na(n_par)) "?" else format(n_par)
  ))

  # THETA
  cat("\n--- THETA Estimates ---\n")
  cat(sprintf("%-16s %12s %12s %10s\n", "Parameter", "Estimate", "SE", "%RSE"))
  cat(strrep("-", 52), "\n", sep = "")
  theta_names <- names(x$theta)
  if (is.null(theta_names)) theta_names <- paste0("THETA", seq_along(x$theta))
  for (i in seq_along(x$theta)) {
    est <- x$theta[i]
    if (!is.null(x$se_theta) && length(x$se_theta) >= i) {
      se_val <- x$se_theta[i]
      rse <- if (abs(est) > 1e-12) abs(se_val / est) * 100 else NaN
      se_str <- sprintf("%.6f", se_val)
      rse_str <- sprintf("%.1f", rse)
    } else {
      se_str <- "N/A"; rse_str <- "N/A"
    }
    cat(sprintf("%-16s %12.6f %12s %10s\n", theta_names[i], est, se_str, rse_str))
  }

  # OMEGA
  cat("\n--- OMEGA Estimates ---\n")
  om <- x$omega
  if (is.null(dim(om))) om <- matrix(om, 1, 1)
  n_eta <- nrow(om)
  has_offdiag <- FALSE
  for (i in seq_len(n_eta)) {
    var_ii <- om[i, i]
    cv_pct <- if (var_ii > 0) sqrt(var_ii) * 100 else 0
    se_str <- if (!is.null(x$se_omega) && length(x$se_omega) >= i) {
      sprintf("%.6f", x$se_omega[i])
    } else {
      "N/A"
    }
    cat(sprintf(
      "  OMEGA(%d,%d) = %.6f  (CV%% = %.1f)  SE = %s\n",
      i, i, var_ii, cv_pct, se_str
    ))
    for (j in seq_len(i - 1L)) {
      if (abs(om[i, j]) > 1e-15) has_offdiag <- TRUE
    }
  }
  if (has_offdiag) {
    cat("  --- Correlations ---\n")
    for (i in seq_len(n_eta)) {
      for (j in seq_len(i - 1L)) {
        cov_ij <- om[i, j]
        var_i <- om[i, i]; var_j <- om[j, j]
        corr <- if (var_i > 0 && var_j > 0) cov_ij / (sqrt(var_i) * sqrt(var_j)) else 0
        cat(sprintf(
          "  OMEGA(%d,%d) = %.6f  (corr = %.4f)\n",
          i, j, cov_ij, corr
        ))
      }
    }
  }

  # SIGMA
  cat("\n--- SIGMA Estimates ---\n")
  for (i in seq_along(x$sigma)) {
    se_str <- if (!is.null(x$se_sigma) && length(x$se_sigma) >= i) {
      sprintf("%.6f", x$se_sigma[i])
    } else {
      "N/A"
    }
    cat(sprintf("  SIGMA(%d) = %.6f  SE = %s\n", i, x$sigma[i], se_str))
  }

  # SIR uncertainty
  if (!is.null(x$sir_ess)) {
    cat("\n--- SIR Uncertainty (95% CI) ---\n")
    cat(sprintf("Effective sample size: %.1f\n", x$sir_ess))
    print_ci <- function(m) {
      if (is.null(m)) return(invisible())
      for (i in seq_len(nrow(m))) {
        cat(sprintf("  %s : [%.6f, %.6f]\n", rownames(m)[i], m[i, 1], m[i, 2]))
      }
    }
    print_ci(x$sir_ci_theta)
    print_ci(x$sir_ci_omega)
    print_ci(x$sir_ci_sigma)
  }

  if (length(x$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in x$warnings) cat("  *", w, "\n")
  }

  cat(bar, "\n", sep = "")
  invisible(x)
}
