#' Simulate from a NLME model
#'
#' Simulates observations from a parsed model with between-subject variability
#' and residual error.
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV (provides population structure: doses, obs times)
#' @param n_sim Number of simulation replicates
#' @param seed Random seed for reproducibility
#'
#' @return A data.frame with columns: SIM, ID, TIME, IPRED, DV_SIM
#'
#' @examples
#' \dontrun{
#' sim <- ferx_simulate("warfarin.ferx", "warfarin.csv", n_sim = 100, seed = 42)
#' }
#'
#' @export
ferx_simulate <- function(model, data, n_sim = 1L, seed = 42L) {
  stopifnot(file.exists(model), file.exists(data))

  ferx_rust_simulate(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    n_sim = as.integer(n_sim),
    seed = as.integer(seed)
  )
}

#' Population predictions from a NLME model
#'
#' Computes population-level predictions (eta = 0) for all subjects.
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV
#'
#' @return A data.frame with columns: ID, TIME, PRED
#'
#' @examples
#' \dontrun{
#' preds <- ferx_predict("warfarin.ferx", "warfarin.csv")
#' }
#'
#' @export
ferx_predict <- function(model, data) {
  stopifnot(file.exists(model), file.exists(data))

  ferx_rust_predict(
    model_path = normalizePath(model),
    data_path = normalizePath(data)
  )
}
