use extendr_api::prelude::*;
use ferx_nlme::types::*;
use nalgebra::DMatrix;
use std::path::Path;

/// Fit a NLME model to data.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param method Estimation method: "foce" or "focei"
/// @param maxiter Maximum outer iterations
/// @param covariance Run covariance step (TRUE/FALSE)
/// @param verbose Print progress (TRUE/FALSE)
/// @param bloq_method BLOQ handling: "drop", "m3", or "" to use the model default
/// @param optimizer Outer optimizer: "slsqp", "lbfgs", "mma", "bobyqa", or "trust_region"
/// @param inner_maxiter Maximum inner (individual) optimizer iterations
/// @param inner_tol Convergence tolerance for the inner optimizer
/// @param mu_referencing Use mu-referencing for ETA initialisation (TRUE/FALSE)
/// @param steihaug_max_iters Maximum CG iterations for Steihaug subproblem (trust_region only)
/// @return Named list with fit results
/// @export
#[extendr]
fn ferx_rust_fit(
    model_path: &str,
    data_path: &str,
    method: &str,
    maxiter: i32,
    covariance: bool,
    verbose: bool,
    bloq_method: &str,
    optimizer: &str,
    inner_maxiter: i32,
    inner_tol: f64,
    mu_referencing: bool,
    steihaug_max_iters: i32,
) -> List {
    let mut parsed =
        match ferx_nlme::parser::model_parser::parse_full_model_file(Path::new(model_path)) {
            Ok(p) => p,
            Err(e) => {
                rprintln!("Error parsing model: {}", e);
                return List::new(0);
            }
        };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return List::new(0);
        }
    };

    // Build fit options
    let mut opts = parsed.fit_options.clone();
    let m = method.to_lowercase();
    if m.contains("saem") {
        opts.method = EstimationMethod::Saem;
        opts.interaction = false;
    } else if m.contains("hybrid") {
        opts.method = EstimationMethod::FoceGnHybrid;
        opts.interaction = false;
    } else if m.contains("gn") || m.contains("gauss") {
        opts.method = EstimationMethod::FoceGn;
        opts.interaction = false;
    } else if m.contains("focei") || m.contains("foce-i") || m.contains("interaction") {
        opts.method = EstimationMethod::FoceI;
        opts.interaction = true;
    } else {
        opts.method = EstimationMethod::Foce;
        opts.interaction = false;
    }
    opts.outer_maxiter = maxiter as usize;
    opts.run_covariance_step = covariance;
    opts.verbose = verbose;

    // Optional R-side override for BLOQ handling. Empty string → keep whatever
    // the model file specified.
    match bloq_method.trim().to_lowercase().as_str() {
        "" => {}
        "m3" => opts.bloq_method = BloqMethod::M3,
        "drop" | "none" | "ignore" => opts.bloq_method = BloqMethod::Drop,
        other => {
            rprintln!(
                "Unknown bloq_method '{}' — expected 'm3' or 'drop' (falling back to model default)",
                other
            );
        }
    }
    // Mirror onto the compiled model so likelihood functions pick it up.
    parsed.model.bloq_method = opts.bloq_method;

    // Inner optimizer settings
    opts.inner_maxiter = inner_maxiter as usize;
    opts.inner_tol = inner_tol;
    opts.mu_referencing = mu_referencing;
    opts.steihaug_max_iters = steihaug_max_iters as usize;

    // Outer optimizer override
    match optimizer.trim().to_lowercase().as_str() {
        "" | "slsqp" => opts.optimizer = Optimizer::Slsqp,
        "lbfgs" => opts.optimizer = Optimizer::Lbfgs,
        "mma" => opts.optimizer = Optimizer::Mma,
        "bobyqa" => opts.optimizer = Optimizer::Bobyqa,
        "trust_region" => opts.optimizer = Optimizer::TrustRegion,
        other => {
            rprintln!(
                "Unknown optimizer '{}' — expected slsqp, lbfgs, mma, bobyqa, or trust_region (falling back to slsqp)",
                other
            );
            opts.optimizer = Optimizer::Slsqp;
        }
    }

    // Build initial parameters (initial values now live in [parameters] block)
    let init_params = parsed.model.default_params.clone();

    // Fit
    let result = match ferx_nlme::fit(&parsed.model, &population, &init_params, &opts) {
        Ok(r) => r,
        Err(e) => {
            rprintln!("Fit error: {}", e);
            return List::new(0);
        }
    };

    // Convert to R list
    fit_result_to_list(&result, &population)
}

/// Simulate from a NLME model.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV (for population structure)
/// @param n_sim Number of simulations
/// @param seed Random seed
/// @return Data frame with ID, TIME, IPRED, DV_SIM columns
/// @export
#[extendr]
fn ferx_rust_simulate(
    model_path: &str,
    data_path: &str,
    n_sim: i32,
    seed: i32,
) -> Robj {
    let model = match ferx_nlme::parse_model_file(Path::new(model_path)) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let results = ferx_nlme::simulate_with_seed(
        &model,
        &population,
        &model.default_params,
        n_sim as usize,
        seed as u64,
    );

    sim_results_to_df(&results)
}

/// Simulate using fitted parameters.
///
/// @param model_path Path to .ferx model file (used for structure + default metadata)
/// @param data_path Path to NONMEM-format CSV (for population structure)
/// @param theta Fitted theta vector
/// @param omega_flat Row-major flattened omega matrix
/// @param omega_dim Side length of the omega matrix
/// @param sigma Fitted sigma vector
/// @param n_sim Number of simulations
/// @param seed Random seed
/// @return Data frame with SIM, ID, TIME, IPRED, DV_SIM columns
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn ferx_rust_simulate_from_fit(
    model_path: &str,
    data_path: &str,
    theta: Vec<f64>,
    omega_flat: Vec<f64>,
    omega_dim: i32,
    sigma: Vec<f64>,
    n_sim: i32,
    seed: i32,
) -> Robj {
    let model = match ferx_nlme::parse_model_file(Path::new(model_path)) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let params = match params_from_fit(&model, &theta, &omega_flat, omega_dim, &sigma) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };

    let results =
        ferx_nlme::simulate_with_seed(&model, &population, &params, n_sim as usize, seed as u64);
    sim_results_to_df(&results)
}

/// Population predictions from a NLME model.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @return Data frame with ID, TIME, PRED columns
/// @export
#[extendr]
fn ferx_rust_predict(
    model_path: &str,
    data_path: &str,
) -> Robj {
    let model = match ferx_nlme::parse_model_file(Path::new(model_path)) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let results = ferx_nlme::predict(&model, &population, &model.default_params);

    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let pred: Vec<f64> = results.iter().map(|r| r.pred).collect();

    data_frame!(ID = id, TIME = time, PRED = pred).into()
}

/// Population predictions using fitted parameters.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param theta Fitted theta vector
/// @param omega_flat Row-major flattened omega matrix
/// @param omega_dim Side length of the omega matrix
/// @param sigma Fitted sigma vector
/// @return Data frame with ID, TIME, PRED columns
/// @export
#[extendr]
fn ferx_rust_predict_from_fit(
    model_path: &str,
    data_path: &str,
    theta: Vec<f64>,
    omega_flat: Vec<f64>,
    omega_dim: i32,
    sigma: Vec<f64>,
) -> Robj {
    let model = match ferx_nlme::parse_model_file(Path::new(model_path)) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let params = match params_from_fit(&model, &theta, &omega_flat, omega_dim, &sigma) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };

    let results = ferx_nlme::predict(&model, &population, &params);

    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let pred: Vec<f64> = results.iter().map(|r| r.pred).collect();

    data_frame!(ID = id, TIME = time, PRED = pred).into()
}

// -- Helper: materialize ModelParameters from R-side theta/omega/sigma --

fn params_from_fit(
    model: &CompiledModel,
    theta: &[f64],
    omega_flat: &[f64],
    omega_dim: i32,
    sigma: &[f64],
) -> std::result::Result<ModelParameters, String> {
    let template = &model.default_params;

    if theta.len() != template.theta.len() {
        return Err(format!(
            "Fit error: theta length {} does not match model ({} expected)",
            theta.len(),
            template.theta.len()
        ));
    }
    if sigma.len() != template.sigma.values.len() {
        return Err(format!(
            "Fit error: sigma length {} does not match model ({} expected)",
            sigma.len(),
            template.sigma.values.len()
        ));
    }
    let d = omega_dim as usize;
    if d != template.omega.dim() {
        return Err(format!(
            "Fit error: omega dim {} does not match model ({} expected)",
            d,
            template.omega.dim()
        ));
    }
    if omega_flat.len() != d * d {
        return Err(format!(
            "Fit error: omega_flat length {} does not match dim²={}",
            omega_flat.len(),
            d * d
        ));
    }

    // R/Rust both use row-major-from-column-major interchangeably here since
    // the matrix is symmetric in the intended use (covariance); reconstruct
    // via DMatrix::from_row_slice to be explicit.
    let omega_mat = DMatrix::from_row_slice(d, d, omega_flat);
    let omega = OmegaMatrix::from_matrix(
        omega_mat,
        template.omega.eta_names.clone(),
        template.omega.diagonal,
    );

    Ok(ModelParameters {
        theta: theta.to_vec(),
        theta_names: template.theta_names.clone(),
        theta_lower: template.theta_lower.clone(),
        theta_upper: template.theta_upper.clone(),
        omega,
        sigma: SigmaVector {
            values: sigma.to_vec(),
            names: template.sigma.names.clone(),
        },
    })
}

// -- Helper: SimulationResult slice → R data frame --

fn sim_results_to_df(results: &[ferx_nlme::api::SimulationResult]) -> Robj {
    let sim: Vec<i32> = results.iter().map(|r| r.sim as i32).collect();
    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let ipred: Vec<f64> = results.iter().map(|r| r.ipred).collect();
    let dv_sim: Vec<f64> = results.iter().map(|r| r.dv_sim).collect();

    data_frame!(
        SIM = sim,
        ID = id,
        TIME = time,
        IPRED = ipred,
        DV_SIM = dv_sim
    )
    .into()
}

// -- Helper: convert FitResult + Population to R named list --

fn fit_result_to_list(result: &FitResult, population: &Population) -> List {
    // Theta
    let theta_names: Vec<String> = result.theta_names.clone();
    let theta_values: Vec<f64> = result.theta.clone();

    // Omega as flat matrix (row-major)
    let n_eta = result.omega.nrows();
    let mut omega_flat: Vec<f64> = Vec::with_capacity(n_eta * n_eta);
    for i in 0..n_eta {
        for j in 0..n_eta {
            omega_flat.push(result.omega[(i, j)]);
        }
    }

    // SE vectors
    let se_theta: Vec<f64> = result.se_theta.clone().unwrap_or_default();
    let se_omega: Vec<f64> = result.se_omega.clone().unwrap_or_default();
    let se_sigma: Vec<f64> = result.se_sigma.clone().unwrap_or_default();

    // SDTAB as data frame
    let sdtab_cols = ferx_nlme::io::output::sdtab(result, population);
    let sdtab = sdtab_to_dataframe(&sdtab_cols);

    // Warnings
    let warnings: Vec<String> = result.warnings.clone();

    let method_label = match result.method {
        EstimationMethod::Saem => "SAEM",
        EstimationMethod::FoceI => "FOCEI",
        EstimationMethod::Foce => "FOCE",
        EstimationMethod::FoceGn => "FOCE-GN",
        EstimationMethod::FoceGnHybrid => "FOCE-GN-Hybrid",
    };

    list!(
        converged = result.converged,
        method = method_label,
        ofv = result.ofv,
        aic = result.aic,
        bic = result.bic,
        n_subjects = result.n_subjects as i32,
        n_obs = result.n_obs as i32,
        n_parameters = result.n_parameters as i32,
        theta = theta_values,
        theta_names = theta_names,
        omega = omega_flat,
        omega_dim = n_eta as i32,
        sigma = result.sigma.clone(),
        se_theta = se_theta,
        se_omega = se_omega,
        se_sigma = se_sigma,
        sdtab = sdtab,
        warnings = warnings
    )
}

fn sdtab_to_dataframe(cols: &[(String, Vec<f64>)]) -> Robj {
    if cols.is_empty() {
        return ().into();
    }

    let mut pairs: Vec<(&str, Robj)> = Vec::new();
    for (name, values) in cols {
        pairs.push((name.as_str(), values.clone().into()));
    }

    // Build data frame via list + class attribute
    let mut df = List::from_pairs(pairs);
    df.set_class(&["data.frame"]).unwrap();

    let n_rows = cols[0].1.len();
    let row_names: Vec<i32> = (1..=n_rows as i32).collect();
    df.set_attrib("row.names", row_names).unwrap();

    df.into()
}

extendr_module! {
    mod ferx;
    fn ferx_rust_fit;
    fn ferx_rust_simulate;
    fn ferx_rust_simulate_from_fit;
    fn ferx_rust_predict;
    fn ferx_rust_predict_from_fit;
}
