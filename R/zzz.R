.onAttach <- function(libname, pkgname) {
  enabled <- tryCatch(ferx_rust_autodiff_enabled(), error = function(e) NA)
  if (isFALSE(enabled)) {
    packageStartupMessage(
      "ferx: built WITHOUT autodiff. Gradient-based fits will be unavailable ",
      "or fall back to finite differences. Rebuild without FERX_NO_AUTODIFF=1 ",
      "(requires the Enzyme Rust toolchain) to enable autodiff."
    )
  }
}
