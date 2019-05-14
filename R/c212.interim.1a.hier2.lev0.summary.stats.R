# c212.BB.summary
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.interim.1a.hier2.lev0.summary.stats.R,v 1.5 2019/05/05 13:18:12 clb13102 Exp clb13102 $"

c212.interim.1a.hier2.lev0.summary.stats <- function(raw)
{
	if (is.null(raw)) {
		print("NULL raw data");
		return(NULL)
	}

    if (M_global$INTERIM_check_summ_name_1a_2(raw)) {
        print("Missing names");
        return(NULL)
    }

	# Check which variables we are monitoring
	monitor = raw$monitor
	theta_mon = monitor[monitor$variable == "theta",]$monitor
	gamma_mon = monitor[monitor$variable == "gamma",]$monitor
	mu.theta_mon = monitor[monitor$variable == "mu.theta",]$monitor
	mu.gamma_mon = monitor[monitor$variable == "mu.gamma",]$monitor
	sigma2.theta_mon = monitor[monitor$variable == "sigma2.theta",]$monitor
	sigma2.gamma_mon = monitor[monitor$variable == "sigma2.gamma",]$monitor

	model = attr(raw, "model")
	if (is.null(model)) {
		print("Model attribute missing");
		return(NULL)
	}

	nchains = raw$chains

	gamma_summ = data.frame(interval = character(0), B = character(0), AE = character(0),
													mean = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	theta_summ = data.frame(interval = character(0), B = character(0), AE = character(0),
													mean = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	mu.gamma_summ = data.frame(interval = character(0), B = character(0), mean = numeric(0),
													hpi_lower = numeric(0), hpi_upper = numeric(0),
													SD = numeric(0), SE = numeric(0))
	mu.theta_summ = data.frame(interval = character(0), B = character(0), mean = numeric(0),
													hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	sigma2.gamma_summ = data.frame(interval = character(0), B = character(0), mean = numeric(0),
													hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	sigma2.theta_summ = data.frame(interval = character(0), B = character(0), mean = numeric(0),
													hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))

	samples_combined = rep(NA, (raw$iter - raw$burnin)*nchains)

	for (i in 1:raw$nIntervals) {
		for (b in 1:raw$nBodySys[i]) {
			bs = raw$B[i,b]
			for (j in 1:raw$nAE[i, b]) {
				AE = raw$AE[i,b,j]

				# gamma
				if (gamma_mon == 1) {
					s = M_global$summaryStats(raw$gamma[, i, b, j, ], nchains)
					row <- data.frame(interval = raw$Intervals[i], B = raw$B[i, b],
								AE = raw$AE[i, b,j], mean = s[1],
								hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
								SE = s[6])
					gamma_summ = rbind(gamma_summ, row)
				}

				# theta
				if (theta_mon == 1) {
					s = M_global$summaryStats(raw$theta[, i, b, j, ], nchains)
					row <- data.frame(interval = raw$Intervals[i], B = raw$B[i, b],
								AE = raw$AE[i, b,j], mean = s[1],
								hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
								SE = s[6])
					theta_summ = rbind(theta_summ, row)
				}
			}

			# mu.gamma
			if (mu.gamma_mon == 1) {
				s = M_global$summaryStats(raw$mu.gamma[, i, b, ], nchains)
				row <- data.frame(interval = raw$Intervals[i], B = raw$B[i, b], mean = s[1],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5], SE = s[6])
				mu.gamma_summ = rbind(mu.gamma_summ, row)
			}

			# mu.theta
			if (mu.theta_mon == 1) {
				s = M_global$summaryStats(raw$mu.theta[, i, b, ], nchains)
				row <- data.frame(interval = raw$Intervals[i], B = raw$B[i, b], mean = s[1],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5], SE = s[6])
				mu.theta_summ = rbind(mu.theta_summ, row)
			}

			# sigma2.theta
			if (sigma2.theta_mon == 1) {
				s = M_global$summaryStats(raw$sigma2.theta[, i, b, ], nchains)
				row <- data.frame(interval = raw$Intervals[i], B = raw$B[i, b], mean = s[1],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5], SE = s[6])
				sigma2.theta_summ = rbind(sigma2.theta_summ, row)
			}

			# sigma2.gamma
			if (sigma2.gamma_mon == 1) {
				s = M_global$summaryStats(raw$sigma2.gamma[, i, b, ], nchains)
				row <- data.frame(interval = raw$Intervals[i], B = raw$B[i, b], mean = s[1],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5], SE = s[6])
				sigma2.gamma_summ = rbind(sigma2.gamma_summ, row)
			}
		}
	}

	rownames(gamma_summ) <- NULL
	rownames(theta_summ) <- NULL
	rownames(mu.gamma_summ) <- NULL
	rownames(mu.theta_summ) <- NULL
	rownames(sigma2.gamma_summ) <- NULL
	rownames(sigma2.theta_summ) <- NULL

	summary.stats = list(theta.summary = theta_summ, gamma.summary = gamma_summ,
								mu.gamma.summary = mu.gamma_summ,
                                mu.theta.summary = mu.theta_summ,
                                sigma2.gamma.summary = sigma2.gamma_summ,
                                sigma2.theta.summary = sigma2.theta_summ,
								monitor = monitor)

	attr(summary.stats, "model") = model

	return(summary.stats)
}

c212.interim.1a.hier2.lev0.print.summary.stats <- function(summ)
{
	if (is.null(summ)) {
		print("NULL summary data");
		return(NULL)
	}

	# Check which variables we are monitoring
	monitor = summ$monitor
	theta_mon = monitor[monitor$variable == "theta",]$monitor
	gamma_mon = monitor[monitor$variable == "gamma",]$monitor
	mu.theta_mon = monitor[monitor$variable == "mu.theta",]$monitor
	mu.gamma_mon = monitor[monitor$variable == "mu.gamma",]$monitor
	sigma2.theta_mon = monitor[monitor$variable == "sigma2.theta",]$monitor
	sigma2.gamma_mon = monitor[monitor$variable == "sigma2.gamma",]$monitor

	model = attr(summ, "model")
	if (is.null(model)) {
		print("Missing model attribute");
		return(NULL)
	}

	if (theta_mon == 1 && !("theta.summary" %in% names(summ))) {
		print("Missing theta.summary data");
		return(NULL)
	}
	if (gamma_mon == 1 && !("gamma.summary" %in% names(summ))) {
		print("Missing gamma.summary data");
		return(NULL)
	}
	if (mu.gamma_mon == 1 && !("mu.gamma.summary" %in% names(summ))) {
		print("Missing mu.gamma.summary data");
		return(NULL)
	}
	if (mu.theta_mon == 1 && !("mu.theta.summary" %in% names(summ))) {
		print("Missing mu.theta.summary data");
		return(NULL)
	}
	if (sigma2.gamma_mon == 1 && !("sigma2.gamma.summary" %in% names(summ))) {
		print("Missing sigma2.gamma.summary data");
		return(NULL)
	}
	if (sigma2.theta_mon == 1 && !("sigma2.theta.summary" %in% names(summ))) {
		print("Missing sigma2.theta.summary data");
		return(NULL)
	}

	cat(sprintf("Variable          Mean        (95%% HPI)          SD       SE\n"))
	cat(sprintf("=============================================================\n"))
	if (gamma_mon == 1) {
		for (i in 1:nrow(summ$gamma.summary)) {
			row = summ$gamma.summary[i, ]
			cat(sprintf("gamma[%s, %s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n", row$interval,
					row$B, row$AE, row$mean, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
		}
	}
	if (theta_mon == 1) {
		for (i in 1:nrow(summ$theta.summary)) {
			row = summ$theta.summary[i, ]
			cat(sprintf("theta[%s, %s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n", row$interval,
					row$B, row$AE, row$mean, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
		}
	}
	if (mu.gamma_mon == 1) {
		for (i in 1:nrow(summ$mu.gamma.summary)) {
			row = summ$mu.gamma.summary[i, ]
			cat(sprintf("mu.gamma[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n", row$interval,
					row$B, row$mean, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
		}
	}
	if (mu.theta_mon == 1) {
		for (i in 1:nrow(summ$mu.theta.summary)) {
			row = summ$mu.theta.summary[i, ]
			cat(sprintf("mu.theta[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n", row$interval,
					row$B, row$mean, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
		}
	}
	if (sigma2.gamma_mon == 1) {
		for (i in 1:nrow(summ$sigma2.gamma.summary)) {
			row = summ$sigma2.gamma.summary[i, ]
			cat(sprintf("sigma2.gamma[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n", row$interval,
					row$B, row$mean, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
		}
	}
	if (sigma2.theta_mon == 1) {
		for (i in 1:nrow(summ$sigma2.theta.summary)) {
			row = summ$sigma2.theta.summary[i, ]
			cat(sprintf("sigma2.theta[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n", row$interval,
					row$B, row$mean, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
		}
	}
}
