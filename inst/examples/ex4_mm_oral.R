## Initialize
library(ferx)
library(vpc)
library(dplyr)

ex <- ferx_example("mm_oral")

## Fit
result <- ferx_fit(
  model = ex$model, 
  data = ex$data, 
  method = "saem"
)

# Simulate and create VPC
sim <- ferx_simulate(
  model = ex$model, 
  data = ex$data, 
  n_sim = 200, 
  seed = 42,
  fit = result
)
obs <- read.csv(file = ex$data)
vpc(
  obs = obs |> 
    mutate(DV = as.numeric(DV)), 
  sim = sim, 
  obs_cols = list(), 
  sim_cols = list(dv = "DV_SIM")
)
