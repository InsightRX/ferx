#' Fit a nonlinear mixed effects model
#'
#' Fits a NLME model using FOCE or FOCEI estimation with a Rust backend
#' powered by automatic differentiation (Enzyme).
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV file (ID, TIME, DV, EVID, AMT, CMT, ...)
#' @param method Estimation method. One of \code{"foce"}, \code{"focei"},
#'   \code{"saem"}, \code{"gn"} (Gauss-Newton / BHHH), or \code{"gn_hybrid"}
#'   (Gauss-Newton followed by a FOCEI polish step).
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
  method <- match.arg(
    tolower(gsub("[^a-z0-9]", "_", method)),
    c("foce", "focei", "saem", "gn", "gn_hybrid")
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
  cat("NLME Fit Result (", x$method, ")\n", sep = "")
  cat("  Converged:", x$converged, "\n")
  cat("  OFV:", round(x$ofv, 4), " AIC:", round(x$aic, 4), " BIC:", round(x$bic, 4), "\n")
  cat("  Subjects:", x$n_subjects, " Observations:", x$n_obs, "\n\n")

  cat("Theta estimates:\n")
  if (!is.null(x$se_theta)) {
    df <- data.frame(
      Estimate = x$theta,
      SE = x$se_theta,
      RSE = paste0(round(abs(x$se_theta / x$theta) * 100, 1), "%")
    )
  } else {
    df <- data.frame(Estimate = x$theta)
  }
  print(df)

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  *", w, "\n")
  }

  invisible(x)
}
