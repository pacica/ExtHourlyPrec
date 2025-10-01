# GPD.scaled.score computes the scaled score of the GPD scale and shape parameters

GPD.scaled.score <- function(par, y, f, dyn_par) {
  # par: scalar at time t of the GPD parameter of which the score is not computed
  # y: scalar of (total observations - threshold) at time t
  # f: scalar indicating the value of the dynamic parameter at time t for which the score has to be computed
  # dyn_par: string indicating the dynamic parameter for which the score has to be computed, ("scale", "shape")

  if (dyn_par == "scale") {
    xi <- par[1]
    sigma <- f

    if (y > 0) {
      # scaling factor
      I <- -(-1 / sigma^2 * (1 - (2 * xi * (xi + 1)) / (2 * xi^2 + 3 * xi + 1)))
      scale <- (I)^(-1)

      # score
      z <- 1 + xi / sigma * y
      p <- -(1 / xi + 1) * 1 / z
      nabla_z_sigma <- -xi * y / sigma^2
      q <- -1 / sigma
      s <- p * nabla_z_sigma + q
    } else {
      s <- scale <- 0
    }
  } else if (dyn_par == "shape") {
    sigma <- par[1]
    xi <- f

    if (y > 0) {
      # scaling factor
      I <- -(1 /
        (xi^2 * (xi + 1)) +
        (xi + 1) / xi * 2 / (2 * xi^2 + 3 * xi + 1) -
        1 / xi^2 * (2 - 1 / (1 + xi)))
      scale <- (I)^(-1)

      # score
      z <- 1 + xi / sigma * y
      p <- -(1 / xi + 1) * 1 / z
      nabla_z_xi <- y / sigma
      q <- 1 / xi^2 * log(z)
      s <- p * nabla_z_xi + q
    } else {
      s <- scale <- 0
    }
  }

  # scaled score
  scaled_score <- s * scale

  return(scaled_score)
}
