# ferx development reference

This file summarises all changes made to ferx-nlme and ferx, the build
setup on this machine, and guidance for future development sessions.

---

## Machine setup (MacBook, username: teun)

### Folder structure
    ~/ferx-nlme     Rust engine (main computation)
    ~/ferx          R package (user-facing wrapper)

Both are cloned from InsightRX and have a myfork remote pointing to
YOUR_USERNAME's GitHub forks.

### Build constraint — limited RAM
The full autodiff/release build hangs due to low available RAM.
Always use this flag combination for all cargo commands:

    --no-default-features --features ci

This skips Enzyme autodiff and uses finite differences instead.
It does not affect the optimizer or mu-referencing code.

### Standard build commands

    # Check compilation (fastest)
    cargo check --no-default-features --features ci

    # Build
    cargo build --release --no-default-features --features ci

    # Run warfarin example
    cargo run --release --no-default-features --features ci --bin ferx \
      -- examples/warfarin.ferx --data data/warfarin.csv

    # Run two-compartment example
    cargo run --release --no-default-features --features ci --bin ferx \
      -- examples/two_cpt_oral_cov.ferx --data data/warfarin.csv

### Install R package after changes
    cd ~/ferx
    R CMD INSTALL .

### Patched files — do not revert these
These were changed during setup to work around the RAM constraint:

    ferx/src/Makevars
      — removed autodiff flags, uses ci build

    ferx/src/rust/Cargo.toml
      — ferx-nlme dependency: default-features = false, features = ci
      — lto = thin (was fat, caused memory hang)

### Enzyme toolchain (installed but not used in dev builds)
    Location: /Users/teun/.rustup/toolchains/nightly-x86_64-apple-darwin
    Plugin:   lib/rustlib/x86_64-apple-darwin/lib/libEnzyme-22.dylib
    Verify:   rustc +enzyme -Z autodiff=Enable - </dev/null 2>&1 | head
    Expected: error[E0601] — means Enzyme loaded correctly

---

## Changes made to ferx-nlme

### 1. New outer optimizers (src/estimation/)

Added two new optimizers selectable via optimizer = X in [fit_options]:

**BOBYQA** (optimizer = bobyqa)
  - Derivative-free quadratic interpolation via NLopt
  - No new dependencies — NLopt already included
  - Robust on poorly-scaled objective surfaces
  - Does not use gradient callback (NLopt handles this automatically
    via the existing if let Some(g) = grad pattern)
  - Files changed: src/types.rs, src/parser/model_parser.rs,
    src/estimation/outer_optimizer.rs

**Newton Trust Region** (optimizer = trust_region)
  - Second-order method via argmin crate + Steihaug CG inner solver
  - Uses exact AD gradients + finite-difference Hessian (first pass)
  - Exact Enzyme Hessian can replace FD Hessian in future
  - Files changed: src/estimation/trust_region.rs (new),
    src/estimation/outer_optimizer.rs, src/estimation/mod.rs,
    Cargo.toml (added argmin = 0.10, argmin-math = 0.4)

### 2. Speed improvements

**Parallel subject evaluation** (src/stats/likelihood.rs)
  - foce_population_nll changed from .iter() to .par_iter()
  - Benefits all optimizers — each subject is independent
  - rayon was already in Cargo.toml

**Warm-start ETAs for trust region** (src/estimation/trust_region.rs)
  - Added cached_etas: Mutex<Vec<DVector<f64>>> to FoceiProblem
  - Mirrors the NloptState.cached_etas pattern already in SLSQP/BOBYQA
  - Mutex used instead of RefCell because argmin traits take &self

### 3. New user-facing settings in [fit_options]

    optimizer     = slsqp        # default (also: lbfgs, mma, bobyqa, trust_region)
    maxiter       = 500          # outer iteration limit
    inner_maxiter = 200          # ETA estimation iterations per subject
    inner_tol     = 1e-8         # ETA convergence tolerance

    Files changed: src/parser/model_parser.rs, src/types.rs,
    src/estimation/inner_optimizer.rs

### 4. Automatic mu-referencing (all estimation methods)

Mu-referencing is a model reparameterization that benefits FOCE, FOCEI,
SAEM, Gauss-Newton, trust region, and BOBYQA simultaneously. It centres
the ETA posterior at zero regardless of theta values, improving
convergence speed and Laplace approximation accuracy.

**Stage 1 — generic ETA initialization** (src/estimation/parameterization.rs)
  - get_eta_init(n_eta, warm_start, mu_refs) shared by all methods
  - Priority: warm-start > mu_k > zero
  - Called by: inner_optimizer.rs, saem.rs, gauss_newton.rs,
    trust_region.rs

**Stage 2 — automatic pattern detection** (src/parser/model_parser.rs)
  - Detects ETA anchor patterns in [individual_parameters]:
      PARAM = THETA * exp(ETA)              → mu = log(THETA)
      PARAM = exp(log(THETA) + ETA)         → mu = log(THETA)
      PARAM = THETA + ETA                   → mu = THETA
      PARAM = THETA * exp(ETA) * covariates → mu = log(THETA)
  - Stores MuRef mappings in CompiledModel.mu_refs (HashMap)
  - Falls back to mu = 0 silently for unrecognised patterns
  - MuRef struct added to src/types.rs

**Stage 3 — apply mu_k at each outer iteration**
  - Computes mu_k from current theta at start of each outer step
  - All inner loops receive mu_k via get_eta_init()
  - SAEM MH proposal centred on mu_k during exploration phase
  - Output ETAs unshifted before writing to sdtab (mean-zero,
    NONMEM-compatible)
  - Files changed: src/api.rs or src/estimation/outer_optimizer.rs,
    src/io/output.rs

### 5. New example file
    examples/warfarin_bobyqa.ferx
      — documented example showing optimizer, inner_maxiter, inner_tol

---

## Changes made to ferx (R package)

### New arguments in ferx_fit()

    ferx_fit(
      model,
      data,
      method         = "focei",
      optimizer      = "slsqp",      # new
      inner_maxiter  = 200,          # new
      inner_tol      = 1e-8,         # new
      mu_referencing = TRUE,         # new
      ...
    )

    Valid optimizers: slsqp, lbfgs, mma, bobyqa, trust_region

### Mu-referencing message
When mu_referencing = TRUE, ferx_fit() prints a message() showing
which ETAs were automatically detected, e.g.:
  "Mu-referencing detected for: ETA_CL, ETA_V, ETA_KA"

### Files changed in ferx
    R/ferx_fit.R              new arguments, validation, mu message, roxygen docs
    src/rust/src/lib.rs       new arguments passed to FitOptions
    R/extendr-wrappers.R      auto-regenerated by rextendr::document()

---

## Testing in RStudio

```r
library(ferx)
ex <- ferx_example("warfarin")

# All optimizers — OFVs should be similar
fit_slsqp <- ferx_fit(ex$model, ex$data, optimizer = "slsqp")
fit_bobyqa <- ferx_fit(ex$model, ex$data, optimizer = "bobyqa")
fit_tr     <- ferx_fit(ex$model, ex$data, optimizer = "trust_region")

# Timing comparison
system.time(ferx_fit(ex$model, ex$data, optimizer = "slsqp"))
system.time(ferx_fit(ex$model, ex$data, optimizer = "bobyqa"))
system.time(ferx_fit(ex$model, ex$data, optimizer = "trust_region"))

# Mu-referencing comparison
fit_mu    <- ferx_fit(ex$model, ex$data, mu_referencing = TRUE)
fit_no_mu <- ferx_fit(ex$model, ex$data, mu_referencing = FALSE)

# Fine-tuned inner loop
fit_fast <- ferx_fit(ex$model, ex$data,
                     optimizer     = "bobyqa",
                     inner_maxiter = 100,
                     inner_tol     = 1e-6)
```

---

## GitHub branches

Both repos have a feature branch pushed to your forks:

    feature/bobyqa-trust-region

Pull requests open against InsightRX/ferx-nlme:main and InsightRX/ferx:main.

To push new changes:

    cd ~/ferx-nlme
    git add -A
    git commit -m "your message"
    git push myfork feature/bobyqa-trust-region

    cd ~/ferx
    git add -A
    git commit -m "your message"
    git push myfork feature/bobyqa-trust-region

---

## Future work not yet implemented

### Exact Hessian via Enzyme (trust region)
Currently trust region uses finite-difference Hessian (n gradient calls
per iteration). Replacing with Enzyme forward-over-reverse AD would give
exact Hessian in ~2x gradient cost regardless of parameter count.
Requires full autodiff build — needs machine with more RAM or CI.
File to change: src/estimation/trust_region.rs Hessian trait impl.

### BOBYQA interpolation point tuning
NLopt BOBYQA default uses 2n+1 interpolation points. Reducing to n+2
makes each iteration cheaper. Check nlopt crate docs for API.

### Expose trust region radius in [fit_options]
Currently hardcoded: initial = 1.0, max = 10.0.
Adding trust_region_radius and trust_region_max_radius as fit_options
keys would give users control over step aggressiveness.

### NUTS sampler (full Bayesian)
nuts-rs crate provides a pure-Rust NUTS sampler. Would add full
posterior estimation comparable to NONMEM BAYES method.

---

## Starting a new Claude Code session

Always start with:

    cd ~/ferx-nlme   (or ~/ferx)
    claude

Opening instruction for ferx-nlme work:

    Read CLAUDE.md for context.
    Use --no-default-features --features ci in all cargo commands.
    Do not touch any autodiff or Enzyme code.

Opening instruction for ferx work:

    Read CLAUDE.md for context.
    Do not change src/rust/Cargo.toml — keep ferx-nlme with
    default-features = false, features = ci, and lto = thin.
    After any Rust changes run rextendr::document() then R CMD INSTALL .
