# ferx

Fast nonlinear mixed effects (NLME) modeling in R, powered by a Rust backend with [Enzyme](https://enzyme.mit.edu/) automatic differentiation for exact gradients.

## Features

- **FOCE/FOCEI estimation** with automatic differentiation
- **Analytical PK models**: 1- and 2-compartment (oral/IV)
- **ODE-based models**: Dormand-Prince RK45 solver for general ODEs
- **NONMEM-compatible**: reads standard NONMEM CSV datasets
- **Model DSL**: define models in `.ferx` text files

## Installation

Requires the [Enzyme Rust toolchain](https://enzyme.mit.edu/rust/) (a custom Rust fork with autodiff support).

```r
# Install from source
install.packages("ferx", repos = NULL, type = "source")
```

Or from the command line:

```bash
R CMD INSTALL .
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
  theta TVCL(0.134, 0.001, 10.0)
  theta TVV(8.1, 0.1, 500.0)
  theta TVKA(1.0, 0.01, 50.0)
  omega ETA_CL ~ 0.07
  sigma PROP_ERR ~ 0.01

[individual_parameters]
  CL = TVCL * exp(ETA_CL)
  V  = TVV  * exp(ETA_V)
  KA = TVKA * exp(ETA_KA)

[structural_model]
  pk one_cpt_oral(cl=CL, v=V, ka=KA)

[error_model]
  DV ~ proportional(PROP_ERR)

[initial_values]
  theta = [0.2, 10.0, 1.5]
  omega = [0.09, 0.04, 0.30]
  sigma = [0.02]
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
