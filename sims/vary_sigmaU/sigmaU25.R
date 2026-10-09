# Reproducibility --------------------------------------------------------------
## Random seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array unique seed using the array ID
sim_seed = as.integer(args)
## Be reproducible queens
set.seed(sim_seed)

# Setup ------------------------------------------------------------------------
## Load libraries
library(MASS) ## for glm.nb
library(countme) ## for smle_nb
## Source data generating/simulation running functions
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/countme/refs/heads/main/sims/run_sett_vary_sigma.R")

# Simulate multiple replications -----------------------------------------------
print(Sys.time()) ## print start time (for reference)
run_sett_vary_dispersion(
  sigmaU = 0.25, ## error standard deviation
  nrep = 50 ## number of simulated replicates
)
