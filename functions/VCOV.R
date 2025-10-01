# vcovdiff.fe estimates the sandwich variance-covariance matrix
# of the static parameters of the FESD-GPD model
# based on the numerical differentiation of the negative log-likelihood
# it returns a list with both the variance-covariance matrix for the
# unconstrained and constrained static parameters
# and the tolerance value used in the numerical approximation of the derivatives

vcovdiff.fe <- function(func, theta, obs, tails, groups) {
  # func: the function to differentiate numerically, the negative log-likelihood of the FESD-GPD model
  # theta: the unconstrained static parameter
  # obs: a data matrix containing the underlying observations (obs$HourlyPrecipitation) and the threshold (obs$tvth_mv)
  # tails: character string indicating the type of tail behavior, "heavy"
  # groups: character vector containing unique labels of cross-sectional units
  
  # function that computes the approximate first derivative
  fstdr <- function(func, theta, dim, obs, tails, groups, eps) {
    grad <- numeric(dim)
    sel <- numeric(dim)

    for (i in 1:dim) {
      sel[i] <- 1

      a <- theta + eps * sel
      fa <- -func(a, obs, tails, groups)
      b <- theta - eps * sel
      fb <- -func(b, obs, tails, groups)
      grad[i] <- (fa - fb) / (2 * eps[i])

      sel[i] <- 0
    }
    return(grad)
  }

  # function that computes the approximate second derivative
  snddr <- function(func, theta, dim, obs, tails, groups, eps) {
    hess <- matrix(NA, dim, dim)
    sel1 <- sel2 <- numeric(dim)

    for (i in 1:dim) {
      for (j in 1:dim) {
        sel1[i] <- 1
        sel2[j] <- 1

        a <- theta + eps * sel1 + eps * sel2
        fa <- -func(a, obs, tails, groups)
        b <- theta
        fb <- -func(b, obs, tails, groups)
        c <- theta - eps * sel1 - eps * sel2
        fc <- -func(c, obs, tails, groups)

        hess[i, j] <- (fa - 2 * fb + fc) / eps^2

        sel1[i] <- 0
        sel2[j] <- 0
      }
    }
    return(hess)
  }

  # function that computes the derivative of reparameterized static parameters
  jacobian <- function(theta, groups) {
    jac <- matrix(0, length(theta), length(theta))

    for (i in 1:length(theta)) {
      if (i %in% 1:(length(groups) * 2)) {
        jac[i, i] <- 1
      } else {
        if (i == (length(groups) * 2 + 1) | i == (length(groups) * 2 + 2)) {
          jac[i, i] <- exp(theta[i])
        } else {
          jac[i, i] <- exp(-theta[i]) / (1 + exp(-theta[i]))^2
        }
      }
    }
    return(jac)
  }

  ############

  n <- nrow(obs[obs$HourlyPrecipitation > obs$tvth_mv, ])
  dim <- length(theta)
  eps_grad <- c(abs(theta[1:(length(groups) * 2)]) * 0.0001, rep(10^(-4), 4))
  eps_hess <- min(eps_grad)

  grad <- fstdr(func, theta, dim, obs, tails, groups, eps_grad)
  hess <- snddr(func, theta, dim, obs, tails, groups, eps_hess)

  # estimated Fisher Info matrix
  Iinv <- solve(-hess / n)
  # estimated outer product
  OP <- grad %*% t(grad) / n
  # sandwich variance covariance matrix
  SNDW <- Iinv %*% OP %*% Iinv
  # delta method
  jacob <- jacobian(theta, groups)
  SNDW_rep <- t(jacob) %*% SNDW %*% jacob

  return(list(
    "SNDW" = SNDW,
    "SNDW_rep" = SNDW_rep,
    "eps_g" = eps_grad,
    "eps_h" = eps_hess
  ))
}
