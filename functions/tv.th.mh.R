tv.th.mh <- function(xx, quant.level, lag, y_burn, freq, symmetric = T) {
  if (freq == "month") {
    # threshold by month
    xx_burn <- subset(xx, xx$year %in% y_burn)
    quant_burn <- as.vector(by(
      xx_burn$HourlyPrecipitation,
      xx_burn$month,
      function(x) try(as.numeric(quantile(x[(!is.na(x))], quant.level)))
    ))

    years <- sort(unique(xx$year))
    tvth_mv <- tvth_exp <- quant_burn

    tvth_y <- rep(as.numeric(y_burn[length(y_burn)]) + 1, length(quant_burn))
    tvth_m <- months <- sort(unique(xx$month))

    if (symmetric == T) {
      # symmetric time window

      for (y in (length(y_burn) + 2):length(years)) {
        if (y <= (length(years) - lag)) {
          # years before the lag
          for (m in months) {
            tvth_mv <- c(
              tvth_mv,
              as.numeric(quantile(
                subset(
                  xx$HourlyPrecipitation,
                  xx$year %in% c(years[y - lag]:years[y + lag]) & xx$month == m
                ),
                quant.level
              ))
            )
            tvth_exp <- c(
              tvth_exp,
              as.numeric(quantile(
                subset(
                  xx$HourlyPrecipitation,
                  xx$year %in% c(min(years):years[y - 1]) & xx$month == m
                ),
                quant.level
              ))
            )
            tvth_m <- c(tvth_m, m)
          }
        } else {
          # years after the lag
          for (m in months) {
            tvth_mv <- c(
              tvth_mv,
              as.numeric(quantile(
                subset(
                  xx$HourlyPrecipitation,
                  xx$year %in%
                    c(years[y - lag]:years[length(years)]) &
                    xx$month == m
                ),
                quant.level
              ))
            )
            tvth_exp <- c(
              tvth_exp,
              as.numeric(quantile(
                subset(
                  xx$HourlyPrecipitation,
                  xx$year %in% c(min(years):years[y - 1]) & xx$month == m
                ),
                quant.level
              ))
            )
            tvth_m <- c(tvth_m, m)
          }
        }
        tvth_y <- c(tvth_y, rep(years[y], length(quant_burn)))
      }
    } else {
      # time window only past
      for (y in (length(y_burn) + 2):length(years)) {
        for (m in months) {
          tvth_mv <- c(
            tvth_mv,
            as.numeric(quantile(
              subset(
                xx$HourlyPrecipitation,
                xx$year %in% c(years[y - lag]:years[y - 1]) & xx$month == m
              ),
              quant.level
            ))
          )
          tvth_exp <- c(
            tvth_exp,
            as.numeric(quantile(
              subset(
                xx$HourlyPrecipitation,
                xx$year %in% c(min(years):years[y - 1]) & xx$month == m
              ),
              quant.level
            ))
          )
          tvth_m <- c(tvth_m, m)
        }
        tvth_y <- c(tvth_y, rep(years[y], length(quant_burn)))
      }
    }
    tvth <- as.data.frame(cbind(tvth_y, tvth_m, tvth_exp, tvth_mv))
    xx_1 <- merge(
      xx,
      tvth,
      by.x = c("year", "month"),
      by.y = c("tvth_y", "tvth_m")
    )
    xx_1 <- xx_1[with(xx_1, order(as.Date(DATE))), ]
  } else if (freq == "hour") {
    # threshold by month and hour
    xx$mh <- paste0("M", substr(xx$DATE, 6, 7), "H", substr(xx$DATE, 12, 13))
    xx_burn <- subset(xx, xx$year %in% y_burn)
    quant_burn <- as.vector(by(
      xx_burn$HourlyPrecipitation,
      xx_burn$mh,
      function(x) try(as.numeric(quantile(x[(!is.na(x))], quant.level)))
    ))
    years <- sort(unique(xx$year))
    tvth_mh <- mh <- sort(unique(xx$mh))
    tvth_mv <- tvth_exp <- quant_burn
    tvth_y <- rep(as.numeric(y_burn[length(y_burn)]) + 1, length(quant_burn))

    if (symmetric == T) {
      # symmetric time window
      for (y in (length(y_burn) + 2):length(years)) {
        if (y <= (length(years) - lag)) {
          for (m in mh) {
            # tvth_mv <- c(tvth_mv,as.numeric(quantile(subset(xx$HourlyPrecipitation,xx$year %in% c(years[y-lag]:years[y-1]) & xx$mh==m),quant.level)))
            tvth_mv <- c(
              tvth_mv,
              as.numeric(quantile(
                subset(
                  xx$HourlyPrecipitation,
                  xx$year %in% c(years[y - lag]:years[y + lag]) & xx$mh == m
                ),
                quant.level
              ))
            )
            tvth_exp <- c(
              tvth_exp,
              as.numeric(quantile(
                subset(
                  xx$HourlyPrecipitation,
                  xx$year %in% c(min(years):years[y - 1]) & xx$mh == m
                ),
                quant.level
              ))
            )
            tvth_mh <- c(tvth_mh, m)
          }
        } else {
          for (m in mh) {
            # tvth_mv <- c(tvth_mv,as.numeric(quantile(subset(xx$HourlyPrecipitation,xx$year %in% c(years[y-lag]:years[y-1]) & xx$mh==m),quant.level)))
            tvth_mv <- c(
              tvth_mv,
              as.numeric(quantile(
                subset(
                  xx$HourlyPrecipitation,
                  xx$year %in%
                    c(years[y - lag]:years[length(years)]) &
                    xx$mh == m
                ),
                quant.level
              ))
            )
            tvth_exp <- c(
              tvth_exp,
              as.numeric(quantile(
                subset(
                  xx$HourlyPrecipitation,
                  xx$year %in% c(min(years):years[y - 1]) & xx$mh == m
                ),
                quant.level
              ))
            )
            tvth_mh <- c(tvth_mh, m)
          }
        }
        tvth_y <- c(tvth_y, rep(years[y], length(quant_burn)))
      }
    } else {
      # time window only past
      for (y in (length(y_burn) + 2):length(years)) {
        for (m in mh) {
          tvth_mv <- c(
            tvth_mv,
            as.numeric(quantile(
              subset(
                xx$HourlyPrecipitation,
                xx$year %in% c(years[y - lag]:years[y - 1]) & xx$mh == m
              ),
              quant.level
            ))
          )
          tvth_exp <- c(
            tvth_exp,
            as.numeric(quantile(
              subset(
                xx$HourlyPrecipitation,
                xx$year %in% c(min(years):years[y - 1]) & xx$mh == m
              ),
              quant.level
            ))
          )
          tvth_mh <- c(tvth_mh, m)
        }
        tvth_y <- c(tvth_y, rep(years[y], length(quant_burn)))
      }
    }
    tvth <- as.data.frame(cbind(tvth_y, tvth_mh, tvth_exp, tvth_mv))
    tvth$tvth_y <- as.numeric(tvth$tvth_y)
    tvth$tvth_exp <- as.numeric(tvth$tvth_exp)
    tvth$tvth_mv <- as.numeric(tvth$tvth_mv)
    xx_1 <- merge(
      xx,
      tvth,
      by.x = c("year", "mh"),
      by.y = c("tvth_y", "tvth_mh")
    )
    xx_1 <- xx_1[with(xx_1, order(as.Date(DATE))), ]
  }
  return(list("quant.level" = quant.level, "dt" = xx_1))
}
