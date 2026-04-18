# ferx

Fast nonlinear mixed effects (NLME) modeling in R, powered by a Rust backend with [Enzyme](https://enzyme.mit.edu/) automatic differentiation for exact gradients.

## Features

- **FOCE/FOCEI estimation** with automatic differentiation
- **Analytical PK models**: 1- and 2-compartment (oral/IV)
- **ODE-based models**: Dormand-Prince RK45 solver for general ODEs
- **NONMEM-compatible**: reads standard NONMEM CSV datasets
- **Model DSL**: define models in `.ferx` text files

## Installation

A Rust installation with the Enzyme AutoDifferentiation engine is required for FeRx to compute gradients. Most likely
you will need to build Rust from source, which may take an hour or so.

See [documentation](https://insightrx.github.io/ferx-nlme/installation.html) for installation instructions. 

### Install the package

After installing Rust, in R run:

```r
devtools::install_github("InsightRX/ferx")
```

Or from a local clone:
```bash
R CMD INSTALL .
```

### Docker

A Docker image is available that bundles the Enzyme toolchain (built from source), ferx CLI, the ferx R package, and RStudio Server — no local Rust/Enzyme setup required. If you're on Windows, this is currently the only supported way of running FeRx.

```bash
# Build (first build takes ~45-60 min; cached after that)
docker build -t ferx:latest .

# Run RStudio Server
docker run --rm -p 8787:8787 -e PASSWORD=ferx ferx:latest
# -> http://localhost:8787   user: rstudio   password: ferx
```

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

MIT — see [LICENSE.md](LICENSE.md).

Copyright © 2026 InsightRX, Inc.
