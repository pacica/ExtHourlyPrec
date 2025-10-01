# par.dyn.GPD returns the constrained dynamics for the
# score-driven GPD scale and shape parameters for a specific cross-sectional unit
# the scale lies on the positive real line, while the shape ranges in (0,1)

par.dyn.GPD <- function(est_par, x, tails) {
  # est_par: vector of constrained static parameters associated to a specific cross-sectional unit
  #           estimated with the FESD-GPD model
  # x: vector of (total observations - threshold) associated to a specific cross-sectional unit
  # tails: character string indicating the type of tail behavior
  
  # settings
  N <- length(x)
  dyn.par <- 2

  omega1 <- est_par[1]
  omega2 <- est_par[2]
  A11 <- est_par[3]
  A22 <- est_par[4]
  B11 <- est_par[5]
  B22 <- est_par[6]

  f_est <- f_tilde_est <- matrix(NA, N + 1, dyn.par)
  s_est <- s_tilde_est <- matrix(NA, N + 1, dyn.par)

  # initialization
  if (B11 != 1 & B22 != 1) {
    f_tilde_est[1, 1] <- omega1 / (1 - B11)
    f_tilde_est[1, 2] <- omega2 / (1 - B22)
  } else {
    f_tilde_est[1, ] <- 0
  }
  f_est[1, ] <- exp(f_tilde_est[1, ])

  for (t in 1:N) {
    # score for scale parameter
    grad <- 1 / f_est[t, 1]
    s_est[t, 1] <- GPD.scaled.score(
      f_est[t, 2],
      x[t],
      f_est[t, 1],
      dyn_par = "scale"
    )
    s_tilde_est[t, 1] <- grad * s_est[t, 1]
    # score for shape parameter
    grad <- grad.rep(f_est[t, 2], tails)
    s_est[t, 2] <- GPD.scaled.score(
      f_est[t, 1],
      x[t],
      f_est[t, 2],
      dyn_par = "shape"
    )
    s_tilde_est[t, 2] <- grad * s_est[t, 2]

    # scale dynamic
    f_tilde_est[t + 1, 1] <- omega1 +
      A11 * s_tilde_est[t, 1] +
      B11 * f_tilde_est[t, 1]
    f_est[t + 1, 1] <- exp(f_tilde_est[t + 1, 1])
    # shape dynamic
    f_tilde_est[t + 1, 2] <- omega2 +
      A22 * s_tilde_est[t, 2] +
      B22 * f_tilde_est[t, 2]
    f_est[t + 1, 2] <- param(f_tilde_est[t + 1, 2], tails)
  }
  return(f_est[-c(1, nrow(f_est)), ])
}
