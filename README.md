# ferx

Fast nonlinear mixed effects (NLME) modeling in R, powered by a Rust backend with [Enzyme](https://enzyme.mit.edu/) automatic differentiation for exact gradients.

## Features

- **FOCE/FOCEI estimation** with automatic differentiation
- **Analytical PK models**: 1- and 2-compartment (oral/IV)
- **ODE-based models**: Dormand-Prince RK45 solver for general ODEs
- **NONMEM-compatible**: reads standard NONMEM CSV datasets
- **Model DSL**: define models in `.ferx` text files

## Installation

ferx depends on the experimental Rust autodiff feature, which requires a nightly rustc plus the Enzyme LLVM plugin. As of 2026, Enzyme is not yet shipped by rustup, so a one-time plugin build is needed.

### Prerequisites

1. **rustup + upstream nightly**
   ```bash
   # If you had snap's rustup, remove it first: sudo snap remove rustup
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
   source "$HOME/.cargo/env"
   rustup toolchain install nightly
   ```

2. **The `enzyme` toolchain must be registered.** The ferx build system pins to a toolchain named `enzyme`. Either:

   **Option A — simple (single user, dev machine):** use your nightly as the enzyme toolchain
   ```bash
   rustup toolchain link enzyme "$(rustup which --toolchain nightly rustc | xargs dirname | xargs dirname)"
   ```

   **Option B — shared (multi-user server):** a sysadmin stages nightly + the Enzyme plugin in `/opt/rust-nightly`, each user links that:
   ```bash
   rustup toolchain link enzyme /opt/rust-nightly
   ```
   See [INSTALL-SYSADMIN.md](INSTALL-SYSADMIN.md) for the sysadmin side: installing matching LLVM, building the Enzyme plugin, and dropping `libEnzyme-<N>.so` in the right sysroot location.

3. **R environment** — tell R to find rustup and pick the enzyme toolchain. Add to `~/.Renviron`:
   ```
   PATH=/opt/rust-nightly/bin:${HOME}/.cargo/bin:${PATH}
   RUSTUP_TOOLCHAIN=enzyme
   ```
   (Drop the `/opt/rust-nightly/bin` entry if you're using Option A above.) Restart R.

### Install the package

From R:
```r
devtools::install_github("InsightRX/ferx")
```

Or from a local clone:
```bash
R CMD INSTALL .
```

### Verifying the Enzyme plugin is available

Before installing, this should print `LLVM version: <N>.x.y` and **not** error:
```bash
rustc +enzyme -Z autodiff=Enable - </dev/null 2>&1 | head
# Expected output: "error[E0601]: `main` function not found"
# That error is fine — it means rustc + Enzyme loaded successfully.
# The real problem is if it says "autodiff backend not found in the sysroot".
```

If you see "autodiff backend not found", the `libEnzyme-<N>.so` file is missing or in the wrong place. See the sysadmin install guide for how to build and place it.

## Quick Start

```r
library(ferx)

# Get bundled example paths
ex <- ferx_example("warfarin")

# Fit a one-compartment oral PK model
result <- ferx_fit(ex$model, ex$data, method = "focei")
result

# Simulate for VPC
sim <- ferx_simulate(ex$model, ex$data, n_sim = 100, seed = 42)

# Population predictions
preds <- ferx_predict(ex$model, ex$data)
```

## Model Specification

Models are defined in `.ferx` files:

```
[parameters]
  theta TVCL(0.2, 0.001, 10.0)   # name(initial, lower, upper)
  theta TVV(10.0, 0.1, 500.0)
  theta TVKA(1.5, 0.01, 50.0)

  omega ETA_CL ~ 0.09            # between-subject variability (variance)
  omega ETA_V  ~ 0.04
  omega ETA_KA ~ 0.30

  sigma PROP_ERR ~ 0.02

[individual_parameters]
  CL = TVCL * exp(ETA_CL)
  V  = TVV  * exp(ETA_V)
  KA = TVKA * exp(ETA_KA)

[structural_model]
  pk one_cpt_oral(cl=CL, v=V, ka=KA)

[error_model]
  DV ~ proportional(PROP_ERR)
```

See `ferx_example()` for available bundled examples.

## API Reference

| Function | Description |
|---|---|
| `ferx_fit()` | Fit a NLME model (FOCE/FOCEI) |
| `ferx_simulate()` | Simulate replicates with BSV and residual error |
| `ferx_predict()` | Population predictions (ETA = 0) |
| `ferx_example()` | Get paths to bundled example models and data |

## License

MIT
