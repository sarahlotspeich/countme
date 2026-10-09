# Reproducibility --------------------------------------------------------------
## Random seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array unique seed using the array ID
sim_seed = as.integer(args)
## Be reproducible queens
set.seed(sim_seed)

# Setup ------------------------------------------------------------------------
## Load libraries
library(bizicount) ## for zic.reg
library(countme) ## for smle_nb and smle_zi_nb
## Source data generating/simulation running functions
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/countme/refs/heads/main/sims/run_sett_vary_zero_infl.R")
source("~/Documents/countme/sims/run_sett_vary_zero_infl.R")

# Simulate multiple replications -----------------------------------------------
print(Sys.time()) ## print start time (for reference)
run_sett_vary_zero_infl(
  eta0 = -4, ## intercept in zero inflation model
  nrep = 2 ## number of simulated replicates
)
