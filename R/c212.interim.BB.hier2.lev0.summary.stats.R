# c212.BB.summary
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.interim.BB.hier2.lev0.summary.stats.R,v 1.6 2019/05/05 13:18:12 clb13102 Exp clb13102 $"

c212.interim.BB.hier2.lev0.summary.stats <- function(raw)
{
	s_base = c212.interim.1a.hier2.lev0.summary.stats(raw)

	# Check which variables we are monitoring
	monitor = raw$monitor
	pi_mon = monitor[monitor$variable == "pi",]$monitor

	if (is.null(s_base)) {
		return(NULL)
	}

	if (pi_mon == 1 && !("pi" %in% names(raw))) {
		print("Missing pi data");
		return(NULL)
	}

	nchains = raw$chains

	pi_summ = data.frame(interval = character(0), B = character(0), mean = numeric(0),
													hpi_lower = numeric(0), hpi_upper = numeric(0),
													SD = numeric(0), SE = numeric(0))

	samples_combined <- rep(NA, (raw$iter - raw$burnin)*nchains)

	if (pi_mon == 1) {
		for (i in 1:raw$nIntervals) {
			for (b in 1:raw$nBodySys[i]) {
				bs = raw$B[i,b]

				s = M_global$summaryStats(raw$pi[, i, b, ], nchains)
				row <- data.frame(interval = raw$Intervals[i], B = raw$B[i, b], mean = s[1],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5], SE = s[6])
				pi_summ = rbind(pi_summ, row)
			}
		}	
	}

	rownames(pi_summ) <- NULL

	s_BB = list(pi.summary = pi_summ)

	summary.stats= c(s_base, s_BB)

	attr(summary.stats, "model") = attr(raw, "model")

	return(summary.stats)
}

c212.interim.BB.hier2.lev0.print.summary.stats <- function(summ)
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
	pi_mon = monitor[monitor$variable == "pi",]$monitor

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

	if (pi_mon == 1 && !("pi.summary" %in% names(summ))) {
		print("Missing pi.summary data");
		return(NULL)
	}

	cat(sprintf("Variable          Mean        (95%% HPI)          SD       SE\n"))
	cat(sprintf("=============================================================\n"))
	if (gamma_mon == 1) {
		for (i in 1:nrow(summ$gamma.summary)) {
			row = summ$gamma.summary[i, ]
			cat(sprintf("gamma[%s, %s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$interval, row$B, row$AE, row$mean, row$hpi_lower,
					row$hpi_upper, row$SD, row$SE))
		}
	}

	if (theta_mon == 1) {
		for (i in 1:nrow(summ$theta.summary)) {
			row = summ$theta.summary[i, ]
			cat(sprintf("theta[%s, %s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$interval, row$B, row$AE, row$mean, row$hpi_lower,
					row$hpi_upper, row$SD, row$SE))
		}
	}
	if (mu.gamma_mon == 1) {
		for (i in 1:nrow(summ$mu.gamma.summary)) {
			row = summ$mu.gamma.summary[i, ]
			cat(sprintf("mu.gamma[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$interval,
					row$B, row$mean, row$hpi_lower,
					row$hpi_upper, row$SD, row$SE))
		}
	}
	if (mu.theta_mon == 1) {
		for (i in 1:nrow(summ$mu.theta.summary)) {
			row = summ$mu.theta.summary[i, ]
			cat(sprintf("mu.theta[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$interval,
					row$B, row$mean, row$hpi_lower,
					row$hpi_upper, row$SD, row$SE))
		}
	}
	if (sigma2.gamma_mon == 1) {
		for (i in 1:nrow(summ$sigma2.gamma.summary)) {
			row = summ$sigma2.gamma.summary[i, ]
			cat(sprintf("sigma2.gamma[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$interval,
					row$B, row$mean, row$hpi_lower,
					row$hpi_upper, row$SD, row$SE))
		}
	}
	if (sigma2.theta_mon == 1) {
		for (i in 1:nrow(summ$sigma2.theta.summary)) {
			row = summ$sigma2.theta.summary[i, ]
			cat(sprintf("sigma2.theta[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$interval,
					row$B, row$mean, row$hpi_lower,
					row$hpi_upper, row$SD, row$SE))
		}
	}
	if (pi_mon == 1) {
		for (i in 1:nrow(summ$pi.summary)) {
			row = summ$pi.summary[i, ]
			cat(sprintf("pi[%s, %s]: %0.6f (%0.6f %0.6f) %0.6f %0.6f\n", row$interval,
						row$B, row$mean, row$hpi_lower,
					row$hpi_upper, row$SD, row$SE))
		}
	}
}
