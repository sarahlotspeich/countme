# Write function to simulate data ----------------------------------------------
sim_data = function(n, sn, sigmaU = 0.25, pv = 0.15, k = 0.5, beta = matrix(data = c(5, 1), ncol = 1)) {
  ## Generate error-free covariate
  x1 = x1f = rnorm(n = n)
  ### Design matrix (add intercept column)
  x = data.matrix(data.frame(int = 1, x1 = x1)) ## n x (p + 1) matrix
  ### Mean parameters for Y|X
  mu = exp(x %*% beta)
  ## Generate outcome
  y = rnbinom(n = n,
              size = k,
              prob = (k / (mu + k)))
  ## Generate error-prone covariate
  x1star = x1 + rnorm(n = n, sd = sigmaU)
  ## Generate validation indicators
  v = sample(x = c(0, 1),
             size = n,
             replace = TRUE,
             prob = c(1 - pv, pv))
  x1[v == 0] = NA
  ## Setup new B-splines
  B = splines::bs(x = x1star, ## Error-prone ALI (from EHR)
                  df = sn,
                  Boundary.knots = range(x1star),
                  intercept = TRUE)
  colnames(B) = paste0("bs", seq(1, sn))

  ## Build dataset
  data = data.frame(y, x1f, x1, x1star, B)
  return(data)
}

# Write function to run multiple reps of this setting --------------------------
run_sett_vary_n_pv = function(n, pv, sn, nrep = 50) {
  ## Initialize empty dataframe for results
  res = data.frame(rep = 1:nrep,
                   n, pv, sn,
                   beta0_gs = NA,
                   beta1_gs = NA,
                   theta_gs = NA,
                   se_beta0_gs = NA,
                   se_beta1_gs = NA,
                   se_theta_gs = NA,
                   beta0_n = NA,
                   beta1_n = NA,
                   theta_n = NA,
                   se_beta0_n = NA,
                   se_beta1_n = NA,
                   se_theta_n = NA,
                   beta0_cc = NA,
                   beta1_cc = NA,
                   theta_cc = NA,
                   se_beta0_cc = NA,
                   se_beta1_cc = NA,
                   se_theta_cc = NA,
                   beta0_smle = NA,
                   beta1_smle = NA,
                   theta_smle = NA,
                   se_beta0_smle = NA,
                   se_beta1_smle = NA,
                   se_theta_smle = NA)

  ## Loop over replications ----------------------------------------------------
  print(Sys.time()) ## print start time (for reference)
  for (r in 1:nrep) {
    ### Simulate data
    rdat = sim_data(n = n, sn = sn)

    ### Fit gold standard model
    gs_fit = glm.nb(formula = y ~ x1f,
                    data = rdat)
    #### Save estimates to res
    res[r, c("beta0_gs", "beta1_gs", "theta_gs")] = with(gs_fit, c(coefficients, theta))
    res[r, c("se_beta0_gs", "se_beta1_gs", "se_theta_gs")] = c(sqrt(diag(vcov(gs_fit))), gs_fit$SE.theta)

    ### Fit naive model
    n_fit = glm.nb(formula = y ~ x1star,
                   data = rdat)
    #### Save estimates to res
    res[r, c("beta0_n", "beta1_n", "theta_n")] = with(n_fit, c(coefficients, theta))
    res[r, c("se_beta0_n", "se_beta1_n", "se_theta_n")] = c(sqrt(diag(vcov(n_fit))), n_fit$SE.theta)

    ### Fit complete-case model
    cc_rdat = rdat[!is.na(rdat$x1), ]
    cc_fit = glm.nb(formula = y ~ x1,
                    data = cc_rdat)
    #### Save estimates to res
    res[r, c("beta0_cc", "beta1_cc", "theta_cc")] = with(cc_fit, c(coefficients, theta))
    res[r, c("se_beta0_cc", "se_beta1_cc", "se_theta_cc")] = c(sqrt(diag(vcov(cc_fit))), cc_fit$SE.theta)

    ### Fit SMLE model
    smle_fit = smle_nb(analysis_formula = y ~ x1,
                       error_formula = paste("x1 ~", paste(paste0("bs", 1:sn), collapse = "+")),
                       data = rdat,
                       no_se = FALSE)
    #### Save estimates to res
    res[r, c("beta0_smle", "beta1_smle", "theta_smle")] = smle_fit$coefficients$Estimate
    res[r, c("se_beta0_smle", "se_beta1_smle", "se_theta_smle")] = smle_fit$coefficients$`Std. Error`

    ### Track progress
    if (r %% 25 == 0) {
      print(r)
      print(Sys.time())
    }

    ### Save results
    res |>
      write.csv(paste0("vary_n", n, "_pv", (pv * 100), "_seed", sim_seed, ".csv"),
                row.names = FALSE)
  }
}
