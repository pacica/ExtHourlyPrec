# GPD.SD.VaR returns a vector containing the dynamic (1-prob.level)% extreme quantile

GPD.SD.VaR <- function(est.sigma, est.xi, tv.th, prob.level, data) {
  # est.sigma: vector of estimated scale dynamic
  # est.xi: vector of estimated shape dynamic
  # tv.th: vector of time varying threshold
  # prob.level: scalar indicating level at which the extreme quantile is estimated
  #             it is the value corresponding to the probability on the left of the quantile, so prob.level > threshold level
  #             for instance, prob.level can be equal to 0.999
  # data: vector of underlying observations

  tv.th <- tv.th[-1] # remove first obs since filtered estimates not available at time t=1
  data <- data[-1] # remove first obs since filtered estimates not available at time t=1
  p <- numeric(length(data))

  # empirical estimate for the probability of exceeding the threshold
  for (t in 1:length(data)) {
    p[t] <- cumsum(data > tv.th)[t] / t
  }
  quant <- tv.th + est.sigma / est.xi * (((1 - prob.level) / (p))^(-est.xi) - 1)

  return(quant)
}
