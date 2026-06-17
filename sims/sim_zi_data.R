# Write function to simulate zero-inflated data --------------------------------
sim_zi_data = function(eta0 = -7.4, eta1 = 6.6, n = 2000, sn = 35, gammaU = 0.25, pv = 0.15, k = 0.3, 
                       beta = matrix(data = c(-1.7, 0.2, 0.8), ncol = 1), 
                       pS = c(0.2500000, 0.9870130, 0.4549098, 0.1450000, 0.0580000,
                              0.2490119, 0.3138501, 0.3316391, 0.3111111, 0.0000000),
                       tprS = 0.95, fprS = 0.05,
                       pM = c(0.996, 0.153, 0.002, 0.000, 0.000,
                              0.494, 0.213, 0.213, 0.955, 0.983)) {
  ## Generate error-free covariate (validated ALI)
  ### Begin with stress indicators (10 per person) from Bernoulli (pS)
  S = rbinom(n = n * 10,
             size = 1,
             prob = rep(x = pS, times = n))
  S_mat = matrix(data = S,
                 nrow = n,
                 ncol = 10,
                 byrow = TRUE)
  x1 = x1f = rowMeans(S_mat)
  
  ## Generate additional error-free covariate (sex)
  z = rbinom(
    n = n,
    size = 1,
    prob = 0.4
    )
  ### Design matrix (add intercept column)
  x = data.matrix(data.frame(int = 1, x1 = x1, z = z)) ## n x (p + 1) matrix
  ### Mean parameters for Y|X
  mu = exp(x %*% beta)
  
  ## Generate outcome
  y = rnbinom(n = n,
              size = k,
              prob = (k / (mu + k)))
  ### Add potential zero inflation
  g = rbinom(n = n,
             size = 1,
             prob = 1 / (1 + exp(- (eta0 + eta1 * z))))
  y[g == 1] = 0 #### if Z = 1, force Y = 0
  
  ## Generate error-prone covariate (EHR ALI)
  ### Begin with stress indicators (10 per person) from Bernoulli (1 / (1 + exp(-(gamma0 + gamma1 S))))
  gamma0 = - log((1 - fprS) / fprS) ### define intercept such that P(S* = 0|S = 0) = FPR
  gamma1 = - log((1 - tprS) / tprS) - gamma0 ### define slope such that P(S* = 1|S = 1) = TPR
  Sstar = rbinom(n = n * 10,
                 size = 1,
                 prob = 1 / (1 + exp(- (gamma0 + gamma1 * S))))
  ### Simulate missingness in EHR version of the covariate
  Sstar_miss = rbinom(n = n * 10,
                      size = 1,
                      prob = rep(x = pM, times = n))
  Sstar[which(Sstar_miss == 1)] = NA
  Sstar_mat = matrix(data = Sstar,
                     nrow = n,
                     ncol = 10,
                     byrow = TRUE)
  x1star = rowMeans(Sstar_mat,
                    na.rm = TRUE)
  
  ## Generate validation indicators
  v = sample(x = c(0, 1),
             size = n,
             replace = TRUE,
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
