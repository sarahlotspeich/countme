# Write function to simulate data ----------------------------------------------
sim_data = function(eta0, n = 2000, sn = 35, sigmaU = 0.25, pv = 0.15, k = 2, beta = matrix(data = c(5, 1), ncol = 1)) {
  ## Generate error-free covariate
  x1 = x1f = rnorm(n = n)
  ## Generate additional error-free covariate
  z = rbinom(n = n,
             size = 1,
             prob = 0.2)
  ### Design matrix (add intercept column)
  x = data.matrix(data.frame(int = 1, x1 = x1)) ## n x (p + 1) matrix
  ### Mean parameters for Y|X
  mu = exp(x %*% beta)
  ## Generate outcome
  y = rnbinom(n = n,
              size = k,
              prob = (k / (mu + k)))
  ### Add potential zero inflation
  g = rbinom(n = n,
             size = 1,
             prob = 1 / (1 + exp(- (eta0 + z))))
  y[g == 1] = 0 #### if Z = 1, force Y = 0
  ## Generate error-prone covariate
  x1star = x1 + rnorm(n = 100)
  ## Generate validation indicators
  v = sample(x = c(0, 1),
             size = n,
             replace = TRUE,
             #prob = c(0.65, 0.35))
             prob = c((1 - pv), pv))
  x1[v == 0] = NA
  ## Setup new B-splines
  B = splines::bs(x = x1star, ## Error-prone ALI (from EHR)
                  df = sn,
                  Boundary.knots = range(x1star),
                  intercept = TRUE)
  colnames(B) = paste0("bs", seq(1, sn))

  ## Build dataset
  data = data.frame(y, x1f, x1, x1star, z, B)
  return(data)
}

# Write function to run multiple reps of this setting --------------------------
run_sett_vary_zero_infl = function(eta0, sn = 35, nrep = 50) {
  ## Initialize empty dataframe for results
  res = data.frame(rep = 1:nrep,
                   eta0,
                   beta0_gs = NA,
                   beta1_gs = NA,
                   eta0_gs = NA,
                   eta1_gs = NA,
                   theta_gs = NA,
                   se_beta0_gs = NA,
                   se_beta1_gs = NA,
                   se_eta0_gs = NA,
                   se_eta1_gs = NA,
                   se_theta_gs = NA,
                   beta0_n = NA,
                   beta1_n = NA,
                   eta0_n = NA,
                   eta1_n = NA,
                   theta_n = NA,
                   se_beta0_n = NA,
                   se_beta1_n = NA,
                   se_eta0_n = NA,
                   se_eta1_n = NA,
                   se_theta_n = NA,
                   beta0_cc = NA,
                   beta1_cc = NA,
                   eta0_cc = NA,
                   eta1_cc = NA,
                   theta_cc = NA,
                   se_beta0_cc = NA,
                   se_beta1_cc = NA,
                   se_eta0_cc = NA,
                   se_eta1_cc = NA,
                   se_theta_cc = NA,
                   beta0_smle = NA,
                   beta1_smle = NA,
                   theta_smle = NA,
                   se_beta0_smle = NA,
                   se_beta1_smle = NA,
                   se_theta_smle = NA,
                   beta0_smle_zi = NA,
                   beta1_smle_zi = NA,
                   eta0_smle_zi = NA,
                   eta1_smle_zi = NA,
                   theta_smle_zi = NA,
                   se_beta0_smle_zi = NA,
                   se_beta1_smle_zi = NA,
                   se_eta0_smle_zi = NA,
                   se_eta1_smle_zi = NA,
                   se_theta_smle_zi = NA)

  ## Loop over replications ----------------------------------------------------
  print(Sys.time()) ## print start time (for reference)
  for (r in 1:nrep) {
    ### Simulate data
    rdat = sim_data(eta0 = eta0)

    ## Fit gold standard model
    gs_fit = zic.reg(
      fmla = y ~ x1f | z,
      data = rdat,
      dist = "nbinom",
      optimizer = "nlm"
    )

    ### Save estimates to res
    res[r, c("beta0_gs", "beta1_gs", "eta0_gs", "eta1_gs", "theta_gs")] = gs_fit$coef
    res[r, c("se_beta0_gs","se_beta1_gs", "se_eta0_gs", "se_eta1_gs", "se_theta_gs")] = sqrt(diag(vcov(gs_fit)))

    ## Fit naive model
    n_fit = zic.reg(
      fmla = y ~ x1star | z,
      data = rdat,
      dist = "nbinom",
      optimizer = "nlm"
    )

    ### Save estimates to res
    res[r, c("beta0_n", "beta1_n", "eta0_n", "eta1_n", "theta_n")] = n_fit$coef
    res[r, c("se_beta0_n","se_beta1_n", "se_eta0_n", "se_eta1_n", "se_theta_n")] = sqrt(diag(vcov(n_fit)))

    ## Fit complete-case model
    rdat_cc = subset(rdat, !is.na(x1))
    cc_fit = zic.reg(
      fmla = y ~ x1 | z,
      data = rdat_cc,
      dist = "nbinom",
      optimizer = "nlm"
    )

    ### Save estimates to res
    res[r, c("beta0_cc", "beta1_cc", "eta0_cc", "eta1_cc", "theta_cc")] = cc_fit$coef
    res[r, c("se_beta0_cc","se_beta1_cc", "se_eta0_cc", "se_eta1_cc", "se_theta_cc")] = sqrt(diag(vcov(cc_fit)))

    ## Fit SMLE model
    smle_fit = smle_nb(analysis_formula = y ~ x1,
                       error_formula = paste("x1 ~", paste(paste0("bs", 1:sn), collapse = "+")),
                       data = rdat,
                       no_se = FALSE)
    ### Save estimates to res
    res[r, c("beta0_smle", "beta1_smle", "theta_smle")] = smle_fit$coefficients$Estimate
    res[r, c("se_beta0_smle", "se_beta1_smle", "se_theta_smle")] = smle_fit$coefficients$`Std. Error`

    ## Fit SMLE model
    smle_fit = smle_zi_nb(analysis_formula = y ~ x1 | z,
                          error_formula = paste("x1 ~", paste(paste0("bs", 1:sn), collapse = "+")),
                          data = rdat,
                          no_se = FALSE)

    ### Save estimates to res
    res[r, c("beta0_smle_zi", "beta1_smle_zi", "eta0_smle_zi", "eta1_smle_zi", "theta_smle_zi")] = with(smle_fit, coefficients$Estimate)
    res[r, c("se_beta0_smle_zi", "se_beta1_smle_zi", "se_eta0_smle_zi", "se_eta1_smle_zi", "se_theta_smle_zi")] = with(smle_fit, coefficients$`Std. Error`)

    ### Track progress
    if (r %% 25 == 0) {
      print(r)
      print(Sys.time())
    }

    ### Save results
    res |>
      write.csv(paste0("vary_zero_infl_neg", abs(eta0), "_seed", sim_seed, ".csv"),
                row.names = FALSE)
  }
}
