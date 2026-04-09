## Initialize
library(ferx)
library(vpc)
library(dplyr)

ex <- ferx_example("warfarin")

## Fit
result <- ferx_fit(
  model = ex$model, 
  data = ex$data, 
  method = "focei"
)

# Simulate and create VPC
sim <- ferx_simulate(ex$model, ex$data, n_sim = 100, seed = 42)
obs <- read.csv(file = ex$data)
vpc(
  obs = obs |> 
    mutate(DV = as.numeric(DV)), 
  sim = sim, 
  obs_cols = list(), 
  sim_cols = list(dv = "DV_SIM")
)
