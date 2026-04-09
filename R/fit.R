#' Fit a nonlinear mixed effects model
#'
#' Fits a NLME model using FOCE or FOCEI estimation with a Rust backend
#' powered by automatic differentiation (Enzyme).
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV file (ID, TIME, DV, EVID, AMT, CMT, ...)
#' @param method Estimation method: "foce" or "focei"
#' @param maxiter Maximum number of outer optimization iterations
#' @param covariance Logical; compute the covariance step for standard errors
#' @param verbose Logical; print progress during estimation
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
#' result <- ferx_fit("warfarin.ferx", "warfarin.csv")
#' result$theta
#' head(result$sdtab)
#' }
#'
#' @export
ferx_fit <- function(model, data,
                     method = "foce",
                     maxiter = 500L,
                     covariance = TRUE,
                     verbose = TRUE) {
  stopifnot(file.exists(model), file.exists(data))
  method <- match.arg(tolower(method), c("foce", "focei"))

  raw <- ferx_rust_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    method = method,
    maxiter = as.integer(maxiter),
    covariance = covariance,
    verbose = verbose
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
