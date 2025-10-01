# this script contains the negative log-likelihood function useful for the estimation of the FESD-GPD model
# GPD.GAS.nllk.i (internal): returns the negative log-likelihood for a specific unit across time
# GPD.GAS.nllk.fe (external): returns the pooled negative log-likelihood for all units across time

GPD.GAS.nllk.i <- function(par, i, ngroups, x, tails) {
  # par: vector of static parameters of the FESD-GPD model
  # i: string indicating label of cross-sectional unit
  # ngroups: total number of units
  # x: vector of (total observations - threshold): only the one where x>0 will be considered for the estimation
  # tails: tail type assumption, "heavy"
  
  # assign parameters
  omega1 <- par[i] # intercept for the scale dynamic
  omega2 <- par[i + ngroups] # intercept for the shape dynamic
  A11 <- exp(par[ngroups * 2 + 1]) # parameter associated to the score of the scale
  A22 <- exp(par[ngroups * 2 + 2]) # parameter associated to the score of the shape
  B11 <- 1 / (1 + exp(-par[ngroups * 2 + 3])) # autoregressive parameter for the scale dynamic
  B22 <- 1 / (1 + exp(-par[ngroups * 2 + 4])) # autoregressive parameter for the shape dynamic

  # settings
  x <- as.numeric(x)
  N <- length(x)
  n_pos <- sum(ifelse((x > 0), 1, 0))
  s <- s_tilde <- matrix(NA, N, 2)
  f <- f_tilde <- matrix(NA, N, 2)

  # initialization
  if (B11 != 1 & B22 != 1) {
    f_tilde[1, 1] = omega1 / (1 - B11)
    f_tilde[1, 2] = omega2 / (1 - B22)
  } else {
    f_tilde[1, ] = 0
  }

  f[1, 1] = exp(f_tilde[1, 1]) # initialize scale
  f[1, 2] = param(f_tilde[1, 2], tails) # initialize shape

  s_tilde[1, ] = s[1, ] = 0

  # score-based recursions
  for (t in 2:N) {
    # unconstrained scale and shape dynamics
    f_tilde[t, 1] <- omega1 + A11 * s_tilde[t - 1, 1] + B11 * f_tilde[t - 1, 1] # scale
    f_tilde[t, 2] <- omega2 + A22 * s_tilde[t - 1, 2] + B22 * f_tilde[t - 1, 2] # shape

    # reparameterized dynamics
    f[t, 1] <- exp(f_tilde[t, 1]) # scale
    f[t, 2] <- param(f_tilde[t, 2], tails) # shape

    # score scale
    grad <- 1 / f[t, 1]
    s[t, 1] <- GPD.scaled.score(f[t, 2], x[t], f[t, 1], dyn_par = "scale")
    s_tilde[t, 1] <- grad * s[t, 1]
    # score shape
    grad <- grad.rep(f[t, 2], tails)
    s[t, 2] <- GPD.scaled.score(f[t, 1], x[t], f[t, 2], dyn_par = "shape")
    s_tilde[t, 2] <- grad * s[t, 2]
  }

  # negative log-likelihood
  sum = 0
  for (t in 1:N) {
    if (x[t] > 0) {
      sum = sum +
        (-log(f[t, 1]) - (1 / f[t, 2] + 1) * log(1 + f[t, 2] / f[t, 1] * x[t]))
    }
  }
  return(-sum)
}

GPD.GAS.nllk.fe <- function(par, x, tails, groups) {
  # par: vector of static parameters of the FESD-GPD model
  # x: vector of (total observations - threshold): only the one where x>0 will be considered for the estimation
  # tails: tail type assumption, "heavy"
  # groups: character vector containing unique labels of cross-sectional units
  
  ngroups <- length(groups)
  nllk <- function(i) {
    GPD.GAS.nllk.i(
      par,
      i,
      length(groups),
      x[x$ID == groups[i], which(colnames(x) == "exc")],
      tails
    )
  }

  # pooled negative log-likelihood to be minimized
  res <- sum(unlist(mclapply(
    1:length(unique(groups)),
    nllk,
    mc.cores = length(unique(groups))
  )))
  return(res)
}
