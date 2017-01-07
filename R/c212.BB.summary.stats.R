# c212.BB.summary
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.BB.summary.stats.R,v 1.4 2016/10/11 12:34:17 clb13102 Exp clb13102 $"

c212.BB.summary.stats <- function(raw)
{
	Id <- "$Id: c212.BB.summary.stats.R,v 1.4 2016/10/11 12:34:17 clb13102 Exp clb13102 $"

	if (is.null(raw)) {
		print("NULL raw data");
		return(NULL)
	}

	model = attr(raw, "model")
	if (is.null(model)) {
		print("Model attribute missing");
		return(NULL)
	}

	n = c("chains", "nBodySys", "maxAEs", "nAE", "theta", "B", "AE", "gamma", "mu.gamma", "mu.theta", "mu.gamma.0", "mu.theta.0",
			"tau2.gamma.0", "tau2.theta.0", "sigma2.theta", "sigma2.gamma", "iter", "burnin")
	if (model == "BB") {
		n = c(n, c("pi", "alpha.pi", "beta.pi"))
	}

	if (M_global$checkNames(n, raw)) {
		print("Missing names");
		return(NULL)
	}

	nchains = raw$chains

	gamma_summ = data.frame(B = character(0), AE = character(0), mean = numeric(0), 
									median = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	theta_summ = data.frame(B = character(0), AE = character(0), mean = numeric(0), 
									median = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	mu.gamma_summ = data.frame(B = character(0), mean = numeric(0), 
									median = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	mu.theta_summ = data.frame(B = character(0), mean = numeric(0), 
									median = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	sigma2.gamma_summ = data.frame(B = character(0), mean = numeric(0), 
									median = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	sigma2.theta_summ = data.frame(B = character(0), mean = numeric(0), 
									median = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
	mu.gamma.0_summ = data.frame(mean = numeric(0), median = numeric(0), hpi_lower = numeric(0), hpi_upper = numeric(0),
													SD = numeric(0), SE = numeric(0))
	mu.theta.0_summ = data.frame(mean = numeric(0), median = numeric(0), hpi_lower = numeric(0), hpi_upper = numeric(0),
													SD = numeric(0), SE = numeric(0))
	tau2.theta.0_summ = data.frame(mean = numeric(0), median = numeric(0), hpi_lower = numeric(0), hpi_upper = numeric(0),
													SD = numeric(0), SE = numeric(0))
	tau2.gamma.0_summ = data.frame(mean = numeric(0), median = numeric(0), hpi_lower = numeric(0), hpi_upper = numeric(0),
													SD = numeric(0), SE = numeric(0))

	if (model == "BB") {
		pi_summ = data.frame(B = character(0), mean = numeric(0), median = numeric(0), hpi_lower = numeric(0),
													hpi_upper = numeric(0), SD = numeric(0), SE = numeric(0))
		pi_combined <- array(NA, dim=c(raw$nBodySys, (raw$iter - raw$burnin)*nchains))
		alpha.pi_summ = data.frame(mean = numeric(0), median = numeric(0), hpi_lower = numeric(0), hpi_upper = numeric(0),
													SD = numeric(0), SE = numeric(0))
		beta.pi_summ = data.frame(mean = numeric(0), median = numeric(0), hpi_lower = numeric(0), hpi_upper = numeric(0),
													SD = numeric(0), SE = numeric(0))
	}

	gamma_combined <- array(NA, dim=c(raw$nBodySys, raw$maxAEs, (raw$iter - raw$burnin)*nchains))
	theta_combined <- array(NA, dim=c(raw$nBodySys, raw$maxAEs, (raw$iter - raw$burnin)*nchains))
	mu.gamma_combined <- array(NA, dim=c(raw$nBodySys, (raw$iter - raw$burnin)*nchains))
	mu.theta_combined <- array(NA, dim=c(raw$nBodySys, (raw$iter - raw$burnin)*nchains))
	sigma2.theta_combined <- array(NA, dim=c(raw$nBodySys, (raw$iter - raw$burnin)*nchains))
	sigma2.gamma_combined <- array(NA, dim=c(raw$nBodySys, (raw$iter - raw$burnin)*nchains))

	for (b in 1:raw$nBodySys) {
		bs = raw$B[b]
		for (j in 1:raw$nAE[b]) {
			AE = raw$AE[b,j]

			# gamma
			s = M_global$summaryStats(raw$gamma[, b, j, ], nchains)
			row <- data.frame(B = raw$B[b], AE = raw$AE[b,j], mean = s[1], median = s[2],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
							SE = s[6])
			gamma_summ = rbind(gamma_summ, row)


			# theta
			s = M_global$summaryStats(raw$theta[, b, j, ], nchains)
			row <- data.frame(B = raw$B[b], AE = raw$AE[b,j], mean = s[1], median = s[2],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
							SE = s[6])
			theta_summ = rbind(theta_summ, row)
		}

		# mu.gamma
		s = M_global$summaryStats(raw$mu.gamma[, b, ], nchains)
		row <- data.frame(B = raw$B[b], mean = s[1], median = s[2],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
							SE = s[6])
		mu.gamma_summ = rbind(mu.gamma_summ, row)

		# mu.theta
		s = M_global$summaryStats(raw$mu.theta[, b, ], nchains)
		row <- data.frame(B = raw$B[b], mean = s[1], median = s[2],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
							SE = s[6])
		mu.theta_summ = rbind(mu.theta_summ, row)

		# sigma2.theta
		s = M_global$summaryStats(raw$sigma2.theta[, b, ], nchains)
		row <- data.frame(B = raw$B[b], mean = s[1], median = s[2],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
							SE = s[6])
		sigma2.theta_summ = rbind(sigma2.theta_summ, row)

		# sigma2.gamma
		s = M_global$summaryStats(raw$sigma2.gamma[, b, ], nchains)
		row <- data.frame(B = raw$B[b], mean = s[1], median = s[2],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
							SE = s[6])
		sigma2.gamma_summ = rbind(sigma2.gamma_summ, row)

		if (model == "BB") {
			s = M_global$summaryStats(raw$pi[, b, ], nchains)
			row <- data.frame(B = raw$B[b], mean = s[1], median = s[2],
							hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
							SE = s[6])
			pi_summ = rbind(pi_summ, row)
		}
	}

	# mu.gamma.0
	s = M_global$summaryStats(raw$mu.gamma.0[,], nchains)
	row <- data.frame(mean = s[1], median = s[2],
					hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
					SE = s[6])
	mu.gamma.0_summ = rbind(mu.gamma.0_summ, row)
	
	# mu.theta.0
	s = M_global$summaryStats(raw$mu.theta.0[,], nchains)
	row <- data.frame(mean = s[1], median = s[2],
					hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
					SE = s[6])
	mu.theta.0_summ = rbind(mu.theta.0_summ, row)

	# tau2.gamma.0
	s = M_global$summaryStats(raw$tau2.gamma.0[,], nchains)
	row <- data.frame(mean = s[1], median = s[2],
					hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
					SE = s[6])
	tau2.gamma.0_summ = rbind(tau2.gamma.0_summ, row)

	# tau2.theta.0
	s = M_global$summaryStats(raw$tau2.theta.0[,], nchains)
	row <- data.frame(mean = s[1], median = s[2],
					hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
					SE = s[6])
	tau2.theta.0_summ = rbind(tau2.theta.0_summ, row)

	if (model == "BB") {
		s = M_global$summaryStats(raw$alpha.pi[,], nchains)
		row <- data.frame(mean = s[1], median = s[2],
					hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
					SE = s[6])
		alpha.pi_summ = rbind(alpha.pi_summ, row)

		s = M_global$summaryStats(raw$beta.pi[,], nchains)
		row <- data.frame(mean = s[1], median = s[2],
					hpi_lower = s[3], hpi_upper = s[4], SD = s[5],
					SE = s[6])
		beta.pi_summ = rbind(beta.pi_summ, row)
	}

	rownames(gamma_summ) <- NULL
	rownames(theta_summ) <- NULL
	rownames(mu.gamma_summ) <- NULL
	rownames(mu.theta_summ) <- NULL
	rownames(sigma2.gamma_summ) <- NULL
	rownames(sigma2.theta_summ) <- NULL
	rownames(mu.gamma.0_summ) <- NULL
	rownames(mu.theta.0_summ) <- NULL
	rownames(tau2.gamma.0_summ) <- NULL
	rownames(tau2.theta.0_summ) <- NULL

	summary.stats = list(theta.summary = theta_summ, gamma.summary = gamma_summ,
								mu.gamma.summary = mu.gamma_summ,
                                mu.theta.summary = mu.theta_summ,
                                sigma2.gamma.summary = sigma2.gamma_summ,
                                sigma2.theta.summary = sigma2.theta_summ,
                                mu.gamma.0.summary = mu.gamma.0_summ,
                                mu.theta.0.summary = mu.theta.0_summ,
                                tau2.gamma.0.summary = tau2.gamma.0_summ,
                                tau2.theta.0.summary = tau2.theta.0_summ)

	if (model == "BB") {
		rownames(pi_summ) <- NULL
		rownames(alpha.pi_summ) <- NULL
		rownames(alpha.pi_summ) <- NULL

		summary.stats$pi.summary <- pi_summ
		summary.stats$alpha.pi.summary <- alpha.pi_summ
		summary.stats$beta.pi.summary <- beta.pi_summ
	}

	attr(summary.stats, "model") = model

	return(summary.stats)
}

c212.BB.print.summary.stats <- function(summ)
{
	if (is.null(summ)) {
		print("NULL summary data");
		return(NULL)
	}

	model = attr(summ, "model")
	if (is.null(model)) {
		print("Missing model attribute");
		return(NULL)
	}

	n = c("theta.summary", "gamma.summary", "mu.gamma.summary", "mu.theta.summary", "sigma2.gamma.summary", "sigma2.theta.summary",
			"mu.gamma.0.summary", "mu.theta.0.summary", "tau2.theta.0.summary", "tau2.gamma.0.summary")
	if (model == "BB") {
		n = c(n, c("pi.summary", "alpha.pi.summary", "beta.pi.summary"))
	}

	if (M_global$checkNames(n, summ)) {
		print("Missing names");
		return(NULL)
	}

	cat(sprintf("Variable          Mean        Median     (95%% HPI)          SD       SE\n"))
	cat(sprintf("========================================================================\n"))
	for (i in 1:nrow(summ$gamma.summary)) {
		row = summ$gamma.summary[i, ]
		cat(sprintf("gamma[%s, %s]: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$B, row$AE, row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	for (i in 1:nrow(summ$theta.summary)) {
		row = summ$theta.summary[i, ]
		cat(sprintf("theta[%s, %s]: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$B, row$AE, row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	for (i in 1:nrow(summ$mu.gamma.summary)) {
		row = summ$mu.gamma.summary[i, ]
		cat(sprintf("mu.gamma[%s]: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$B, row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	for (i in 1:nrow(summ$mu.theta.summary)) {
		row = summ$mu.theta.summary[i, ]
		cat(sprintf("mu.theta[%s]: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$B, row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	for (i in 1:nrow(summ$sigma2.gamma.summary)) {
		row = summ$sigma2.gamma.summary[i, ]
		cat(sprintf("sigma2.gamma[%s]: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$B, row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	for (i in 1:nrow(summ$sigma2.theta.summary)) {
		row = summ$sigma2.theta.summary[i, ]
		cat(sprintf("sigma2.theta[%s]: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$B, row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	if (model == "BB") {
		for (i in 1:nrow(summ$pi.summary)) {
			row = summ$pi.summary[i, ]
			cat(sprintf("pi[%s]: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
						row$B, row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
		}
	}
	for (i in 1:nrow(summ$mu.gamma.0.summary)) {
		row = summ$mu.gamma.0.summary[i, ]
		cat(sprintf("mu.gamma.0: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	for (i in 1:nrow(summ$mu.theta.0.summary)) {
		row = summ$mu.theta.0.summary[i, ]
		cat(sprintf("mu.theta.0: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	for (i in 1:nrow(summ$tau2.gamma.0.summary)) {
		row = summ$tau2.gamma.0.summary[i, ]
		cat(sprintf("tau2.gamma.0: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	for (i in 1:nrow(summ$tau2.theta.0.summary)) {
		row = summ$tau2.theta.0.summary[i, ]
		cat(sprintf("tau2.theta.0: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
	if (model == "BB") {
		row = summ$alpha.pi.summary[i, ]
		cat(sprintf("alpha.pi: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
		row = summ$beta.pi.summary[i, ]
		cat(sprintf("beta.pi: %0.6f %0.6f (%0.6f %0.6f) %0.6f %0.6f\n",
					row$mean, row$median, row$hpi_lower, row$hpi_upper, row$SD, row$SE))
	}
}
