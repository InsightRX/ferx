use extendr_api::prelude::*;
use ferx_nlme::types::*;
use std::path::Path;

/// Fit a NLME model to data.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param method Estimation method: "foce" or "focei"
/// @param maxiter Maximum outer iterations
/// @param covariance Run covariance step (TRUE/FALSE)
/// @param verbose Print progress (TRUE/FALSE)
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
) -> List {
    let parsed = match ferx_nlme::parser::model_parser::parse_full_model_file(Path::new(model_path)) {
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
    opts.interaction = method.to_lowercase().contains("focei")
        || method.to_lowercase().contains("interaction");
    opts.outer_maxiter = maxiter as usize;
    opts.run_covariance_step = covariance;
    opts.verbose = verbose;

    // Build initial parameters
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

    // Convert to R data frame
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

    list!(
        converged = result.converged,
        method = if result.interaction { "FOCEI" } else { "FOCE" },
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
    fn ferx_rust_predict;
}
