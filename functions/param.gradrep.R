# this script contains the functions used in the score-driven recursion of the shape parameter
# param (internal): if tails == "heavy", returns the dynamic of the shape parameter constrained between (0,1)
# grad.rep (internal): returns the derivative of the unconstrained dynamic with respect to the constrained dynamic

param <- function(f_tilde, tails) {
  # f_tilde: vector of unconstrained shape dynamic defined on the real line
  # tails: character string indicating the type of tail behavior
  
  if (tails == "heavy") {
    f <- 1 / (1 + exp(-f_tilde))
  } else {
    f <- f_tilde
  }
  return(f)
}


grad.rep <- function(f, tails) {
  # f: vector of constrained shape dynamic
  # tails: character string indicating the type of tail behavior
  
  if (tails == "light") {
    der <- 1
  } else if (tails == "norm") {
    der <- 1
  } else {
    if (tails == "heavy") {
      der <- 1 / (f * (1 - f))
    }
  }
  return(der)
}
