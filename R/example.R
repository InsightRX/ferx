#' Get paths to example model and data files
#'
#' Returns file paths to bundled example models and datasets.
#' Called without arguments, lists available example names.
#'
#' @param name Name of the example (e.g. "warfarin"). If NULL, returns
#'   a character vector of available example names.
#'
#' @return If \code{name} is NULL, a character vector of available examples.
#'   Otherwise, a list with components:
#'   \item{model}{Path to the .ferx model file}
#'   \item{data}{Path to the NONMEM-format CSV data file}
#'
#' @examples
#' ferx_example()
#' ex <- ferx_example("warfarin")
#' ex$model
#' ex$data
#'
#' @export
ferx_example <- function(name = NULL) {
  models_dir <- system.file("examples", "models", package = "ferx")
  if (models_dir == "") {
    stop("Example files not found. Is ferx installed?")
  }

  available <- tools::file_path_sans_ext(list.files(models_dir, pattern = "\\.ferx$"))

  if (is.null(name)) {
    return(available)
  }

  if (!name %in% available) {
    stop(
      "Example '", name, "' not found. Available examples: ",
      paste(available, collapse = ", ")
    )
  }

  list(
    model = system.file("examples", "models", paste0(name, ".ferx"), package = "ferx"),
    data = system.file("examples", "data", paste0(name, ".csv"), package = "ferx")
  )
}
