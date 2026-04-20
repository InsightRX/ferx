# Trust region subproblem solver — findings and next steps

---

## Which subproblem solver is used

The ferx-nlme trust region implementation uses **Steihaug-CG** via the
argmin crate. The relevant code in src/estimation/trust_region.rs is:

```rust
let subproblem = Steihaug::new().with_max_iters(20);
let solver = TrustRegion::new(subproblem)
    .with_radius(1.0)
    .with_max_radius(10.0);
```

---

## Why this matters for mu-referencing

The three subproblem solver options and their sensitivity to conditioning:

| Solver | Conditioning sensitivity | Mu-referencing benefit |
|---|---|---|
| Steihaug-CG | Highest | Largest |
| Dogleg | Medium | Moderate (has Cauchy fallback) |
| Exact | Lowest | Smallest |

Steihaug-CG is the most sensitive to Hessian conditioning — it gains
the most from mu-referencing because a well-conditioned Hessian allows
the CG iterations to converge faster and more reliably.

Dogleg has an explicit fallback to the steepest descent (Cauchy)
direction when conditioning is poor, so it degrades more gracefully
but gains less from mu-referencing.

## Conclusion

Choosing Steihaug-CG and implementing mu-referencing is the **best
possible pairing** from a numerical standpoint. The two features
reinforce each other:

- Mu-referencing centres the ETA posterior at zero, giving Steihaug-CG
  a well-conditioned Hessian to work with
- Steihaug-CG exploits good conditioning to converge the subproblem
  in fewer CG iterations per outer step
- Result: faster and more reliable trust region convergence compared
  to either feature alone

---

## Current limitation — Steihaug max iterations too low

The current `with_max_iters(20)` may truncate CG too early for models
with more than 20 parameters. Steihaug-CG needs at most n iterations
where n is the number of parameters (thetas + omegas + sigmas).

For a typical 2-compartment model with covariates, n can reach 15-25.
For complex models it can be higher.

### Fix — expose steihaug_max_iters as a tunable setting

**In ferx-nlme:**

Ask Claude Code in ~/ferx-nlme:

```
Read src/estimation/trust_region.rs, src/types.rs,
and src/parser/model_parser.rs.

Add steihaug_max_iters as an optional field in FitOptions with default 50.
Pass it to Steihaug::new().with_max_iters() in the trust region solver.
Expose it as an optional key in [fit_options] in the .ferx parser.
Verify with: cargo check --no-default-features --features ci
```

**In ferx (R package):**

After ferx-nlme change is done, ask Claude Code in ~/ferx:

```
Add steihaug_max_iters as a formal argument to ferx_fit() with default 50.
Pass it through src/rust/src/lib.rs to FitOptions.
Add it to the roxygen documentation alongside inner_maxiter and inner_tol.
Run rextendr::document() then R CMD INSTALL .
```

**Usage in .ferx file after change:**

```
[fit_options]
  optimizer          = trust_region
  maxiter            = 300
  inner_maxiter      = 100
  inner_tol          = 1e-6
  steihaug_max_iters = 50
```

**Usage in R after change:**

```r
fit <- ferx_fit(ex$model, ex$data,
                optimizer          = "trust_region",
                steihaug_max_iters = 50)
```

---

## Recommended default values by model complexity

| Model type | n params | steihaug_max_iters |
|---|---|---|
| 1-compartment, 2-3 ETAs | ~8 | 20 (current default fine) |
| 2-compartment, 3-4 ETAs | ~15 | 30 |
| 2-compartment + covariates | ~20 | 50 (new default) |
| Complex ODE, 5+ ETAs | 25+ | n_params × 2 |

Setting default to 50 covers most real population PK models safely.

---

## Push changes to GitHub

After implementing the fix in both repos:

```bash
# ferx-nlme
cd ~/ferx-nlme
git add -A
git commit -m "feat: expose steihaug_max_iters as tunable fit option

Default 50 (was hardcoded 20). Prevents early CG truncation
on models with more than 20 parameters."
git push myfork feature/bobyqa-trust-region

# ferx
cd ~/ferx
git add -A
git commit -m "feat: add steihaug_max_iters argument to ferx_fit()

Exposes trust region CG iteration limit with default 50.
Documents alongside inner_maxiter and inner_tol."
git push myfork feature/bobyqa-trust-region
```

---

## Summary of all tunable trust region settings after this change

| Setting | Default | Controls |
|---|---|---|
| maxiter | 500 | Outer Newton iterations |
| inner_maxiter | 200 | ETA BFGS iterations per subject |
| inner_tol | 1e-8 | ETA convergence tolerance |
| steihaug_max_iters | 50 | CG iterations for subproblem |
| trust_region_radius | 1.0 | Initial step radius (hardcoded) |
| trust_region_max_radius | 10.0 | Maximum step radius (hardcoded) |

The last two (radius settings) are still hardcoded and could be exposed
as future fit_options keys if needed.
