#c212.BB.convergence.diag # Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014
#
# If the MCMC simulation has been run for more than one chain report the Gelman-Rubin statistic.
# If the MCMC simulation has been run for only one chain report the Geweke diagnostic (Z-score)
#
# The long-handed implementation here is actually quicker than using apply on each set of data:
#
# dd = function(x, nchains) {
#    mcmc_obj <- list(NA)
#
#    for (c in 1:nchains) {
#        mcmc_obj[[c]] = mcmc(x[c,])
#    }
#
#    mlist <- mcmc.list(mcmc_obj)
#    g <- gelman.diag(mlist)
#    return(c(g$psrf[1], g$psrf[2]))
# }
#
# th = apply(raw$theta, c(2,3), dd, nchains)
# ga = apply(raw$gamma, c(2,3), dd, nchains)
# mg = apply(raw$mu.gamma, 2, dd, nchains)
# mt = apply(raw$mu.theta, 2, dd, nchains)
# sg = apply(raw$sigma2.gamma, 2, dd, nchains)
# st = apply(raw$sigma2.theta, 2, dd, nchains)
#
# This really needs to be moved to C++ code.
#

Id <- "$Id: c212.BB.convergence.R,v 1.6 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.BB.convergence.diag <- function(raw)
{
	if (is.null(raw)) {
		print("NULL raw data")
		return(NULL)
	}

	model = attr(raw, "model")
	if (is.null(model)) {
		print("Simulation model attribute missing")
		return(NULL)
	}

	n = c("chains", "nBodySys", "maxAEs", "nAE", "theta", "B", "AE", "gamma", "mu.gamma", "theta",
			"mu.theta", "mu.gamma.0", "mu.theta.0", "tau2.gamma.0", "tau2.theta.0", "sigma2.theta", "sigma2.gamma")

	if (model == "BB") {
		n = c(n, c("alpha.pi", "beta.pi", "pi"))
	}

	if (raw$sim_type == "MH") {
		n = c(n, c("theta_acc", "gamma_acc"))
		if (model == "BB") {
			n = c(n, c("alpha.pi_acc", "beta.pi_acc", "theta_zero_prop", "theta_zero_acc"))
		}
	}
	else {
		if (model == "BB") {
			n = c(n, c("theta_acc", "theta_zero_prop", "theta_zero_acc"))
		}
	}

	if (M_global$checkNames(n, raw)) {
		print("Missing names");
		return(NULL)
	}

	nchains = raw$chains

	len = sum(raw$nAE)
	gamma_conv = data.frame(B = character(len), AE = character(len), stat = numeric(len), upper_ci = numeric(len),
								stringsAsFactors=FALSE)
	theta_conv = data.frame(B = character(len), AE = character(len), stat = numeric(len), upper_ci = numeric(len),
								stringsAsFactors=FALSE)
	mu.gamma_conv = data.frame(B = character(raw$nBodySys), stat = numeric(raw$nBodySys), upper_ci = numeric(raw$nBodySys),
								stringsAsFactors=FALSE)
	mu.theta_conv = data.frame(B = character(raw$nBodySys), stat = numeric(raw$nBodySys), upper_ci = numeric(raw$nBodySys),
								stringsAsFactors=FALSE)
	sigma2.gamma_conv = data.frame(B = character(raw$nBodySys), stat = numeric(raw$nBodySys), upper_ci = numeric(raw$nBodySys),
								stringsAsFactors=FALSE)
	sigma2.theta_conv = data.frame(B = character(raw$nBodySys), stat = numeric(raw$nBodySys), upper_ci = numeric(raw$nBodySys),
								stringsAsFactors=FALSE)
	mu.gamma.0_conv = data.frame(stat = numeric(1), upper_ci = numeric(1), stringsAsFactors=FALSE)
	mu.theta.0_conv = data.frame(stat = numeric(1), upper_ci = numeric(1), stringsAsFactors=FALSE)
	tau2.theta.0_conv = data.frame(stat = numeric(1), upper_ci = numeric(1), stringsAsFactors=FALSE)
	tau2.gamma.0_conv = data.frame(stat = numeric(1), upper_ci = numeric(1), stringsAsFactors=FALSE)

	if (model == "BB") {
		pi_conv = data.frame(B = character(raw$nBodySys), stat = numeric(raw$nBodySys), upper_ci = numeric(raw$nBodySys),
								stringsAsFactors=FALSE)
		alpha.pi_conv = data.frame(stat = numeric(1), upper_ci = numeric(1), stringsAsFactors=FALSE)
		beta.pi_conv = data.frame(stat = numeric(1), upper_ci = numeric(1), stringsAsFactors=FALSE)
	}

	type <- NA

	if (nchains > 1) {
		# Gelman-Rubin Statistics

		type = "Gelman-Rubin"

		indx = 1
		for (b in 1:raw$nBodySys) {
			for (j in 1:raw$nAE[b]) {
				# theta
				g = M_global$GelmanRubin(raw$theta[ , b, j, ], nchains)
				theta_conv[indx,] = data.frame(B = raw$B[b], AE = raw$AE[b,j],
												stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

				# gamma
				g = M_global$GelmanRubin(raw$gamma[ , b, j, ], nchains)
				gamma_conv[indx,] = data.frame(B = raw$B[b], AE = raw$AE[b,j],
												stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)
				indx = indx + 1
			}

			# mu.gamma
			g = M_global$GelmanRubin(raw$mu.gamma[ , b, ], nchains)
			mu.gamma_conv[b,] = data.frame(B = raw$B[b], stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

	 		# mu.theta
			g = M_global$GelmanRubin(raw$mu.theta[ , b, ], nchains)
			mu.theta_conv[b,] = data.frame(B = raw$B[b], stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

			# sigma2.theta
			g = M_global$GelmanRubin(raw$sigma2.theta[ , b, ], nchains)
			sigma2.theta_conv[b,] = data.frame(B = raw$B[b], stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

			# sigma2.gamma
			g = M_global$GelmanRubin(raw$sigma2.gamma[ , b, ], nchains)
			sigma2.gamma_conv[b,] = data.frame(B = raw$B[b], stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

			if (model == "BB") {
				# pi
				g = M_global$GelmanRubin(raw$pi[ , b, ], nchains)
				pi_conv[b,] = data.frame(B = raw$B[b], stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)
			}
		}

		# mu.gamma.0
		g = M_global$GelmanRubin(raw$mu.gamma.0, nchains)
		mu.gamma.0_conv[1,] = data.frame(stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

		# mu.theta.0
		g = M_global$GelmanRubin(raw$mu.theta.0, nchains)
		mu.theta.0_conv[1,] = data.frame(stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

		# tau2.gamma.0
		g = M_global$GelmanRubin(raw$tau2.gamma.0, nchains)
		tau2.gamma.0_conv[1,] = data.frame(stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

		# tau2.theta.0
		g = M_global$GelmanRubin(raw$tau2.theta.0, nchains)
		tau2.theta.0_conv[1,] = data.frame(stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

		if (model == "BB") {
			# alpha.pi
			g = M_global$GelmanRubin(raw$alpha.pi, nchains)
			alpha.pi_conv[1,] = data.frame(stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)

			# beta.pi
			g = M_global$GelmanRubin(raw$beta.pi, nchains)
			beta.pi_conv[1,] = data.frame(stat = g$psrf[1], upper_ci =  g$psrf[2], stringsAsFactors=FALSE)
		}
	}
	else {
		# Geweke Diagnostic

		type = "Geweke"

		indx = 1
		for (b in 1:raw$nBodySys) {
			for (j in 1:raw$nAE[b]) {
				# theta
				g = M_global$Geweke(raw$theta[1, b, j, ])
				theta_conv[indx,] = data.frame(B = raw$B[b], AE = raw$AE[b,j], stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)

				# gamma
				g = M_global$Geweke(raw$gamma[1, b, j, ])
				gamma_conv[indx,] = data.frame(B = raw$B[b], AE = raw$AE[b,j], stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)

				indx = indx + 1
			}

			# mu.gamma
			g = M_global$Geweke(raw$mu.gamma[1, b, ])
			mu.gamma_conv[b,] = data.frame(B = raw$B[b], stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)

			# mu.theta
			g = M_global$Geweke(raw$mu.theta[1, b, ])
			mu.theta_conv[b,] = data.frame(B = raw$B[b], stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)

			# sigma2.theta
			g = M_global$Geweke(raw$sigma2.theta[1, b, ])
			sigma2.theta_conv[b,] = data.frame(B = raw$B[b], stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)
			
			# sigma2.gamma
			g = M_global$Geweke(raw$sigma2.gamma[1, b, ])
			sigma2.gamma_conv[b,] = data.frame(B = raw$B[b], stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)

			if (model == "BB") {
				# pi
				g = M_global$Geweke(raw$pi[1, b, ])
				pi_conv[b,] = data.frame(B = raw$B[b], stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)
			}
		}

		g = M_global$Geweke(raw$mu.gamma.0[1,])
		mu.gamma.0_conv[1,] = data.frame(stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)
	
		g = M_global$Geweke(raw$mu.theta.0[1,])
		mu.theta.0_conv[1,] = data.frame(stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)
	
		g = M_global$Geweke(raw$tau2.gamma.0[1,])
		tau2.gamma.0_conv[1,] = data.frame(stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)
		
		g = M_global$Geweke(raw$tau2.theta.0[1,])
		tau2.theta.0_conv[1,] = data.frame(stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)

		if (model == "BB") {
			g = M_global$Geweke(raw$alpha.pi[1,])
			alpha.pi_conv[1,] = data.frame(stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)

			g = M_global$Geweke(raw$beta.pi[1,])
			beta.pi_conv[1,] = data.frame(stat = g$z, upper_ci = NA, stringsAsFactors=FALSE)
		}
	}

	theta_acc = data.frame(chain = numeric(0), B = character(0), AE = character(0), rate = numeric(0), stringsAsFactors=FALSE)
	gamma_acc = data.frame(chain = numeric(0), B = character(0), AE = character(0), rate = numeric(0), stringsAsFactors=FALSE)

	if (model == "BB") {
		theta_zero_prop <- data.frame(chain = numeric(0), B = character(0),
										AE = character(0), rate = numeric(0))
		theta_zero_acc <- data.frame(chain = numeric(0), B = character(0),
										AE = character(0), rate = numeric(0))

		alpha.pi_acc <- rep(NA, nchains)
		beta.pi_acc <- rep(NA, nchains)
	}

	if (raw$sim_type == "MH") {
		for (b in 1:raw$nBodySys) {
			for (j in 1:raw$nAE[b]) {
				for (c in 1:nchains) {
					rate <- raw$theta_acc[c, b, j]/raw$iter
					row = data.frame(chain = c, B = raw$B[b], AE = raw$AE[b,j], rate = rate, stringsAsFactors=FALSE)
					theta_acc = rbind(theta_acc, row)
				}
			}
		}

		for (b in 1:raw$nBodySys) {
			for (j in 1:raw$nAE[b]) {
				for (c in 1:nchains) {
					rate <- raw$gamma_acc[c, b, j]/raw$iter
					row = data.frame(chain = c, B = raw$B[b], AE = raw$AE[b,j], rate = rate, stringsAsFactors=FALSE)
					gamma_acc = rbind(gamma_acc, row)
				}
			}
		}

		if (model == "BB") {
			for (c in 1:nchains) {
				alpha.pi_acc[c] <- raw$alpha.pi_acc[c]/raw$iter
				beta.pi_acc[c] <- raw$beta.pi_acc[c]/raw$iter
			}
			for (b in 1:raw$nBodySys) {
				for (j in 1:raw$nAE[b]) {
					for (c in 1:nchains) {
						rate <- raw$theta_zero_prop[c, b, j]/raw$iter
						row = data.frame(chain = c, B = raw$B[b], AE = raw$AE[b,j], rate = rate, stringsAsFactors=FALSE)
						theta_zero_prop = rbind(theta_zero_prop, row)

						rate <- raw$theta_zero_acc[c, b, j]/raw$theta_zero_prop[c, b, j]
						row = data.frame(chain = c, B = raw$B[b], AE = raw$AE[b,j], rate = rate)
						theta_zero_acc = rbind(theta_zero_acc, row)
					}
				}
			}
		}
	}
	else {
		if (model == "BB") {
			for (b in 1:raw$nBodySys) {
				for (j in 1:raw$nAE[b]) {
					for (c in 1:nchains) {
						rate = raw$theta_acc[c, b, j]/raw$iter
						row = data.frame(chain = c, B = raw$B[b], AE = raw$AE[b,j], rate = rate, stringsAsFactors=FALSE)
						theta_acc = rbind(theta_acc, row)

						rate = raw$theta_zero_prop[c, b, j]/raw$iter
						row = data.frame(chain = c, B = raw$B[b], AE = raw$AE[b,j], rate = rate, stringsAsFactors=FALSE)
						theta_zero_prop = rbind(theta_zero_prop, row)

						rate = raw$theta_zero_acc[c, b, j]/raw$theta_zero_prop[c, b, j]
						row = data.frame(chain = c, B = raw$B[b], AE = raw$AE[b,j], rate = rate, stringsAsFactors=FALSE)
						theta_zero_acc = rbind(theta_zero_acc, row)
					}
				}
			}
		}
	}
	
	rownames(gamma_conv) <- NULL
	rownames(theta_conv) <- NULL
	rownames(mu.gamma_conv) <- NULL
	rownames(mu.theta_conv) <- NULL
	rownames(sigma2.gamma_conv) <- NULL
	rownames(sigma2.theta_conv) <- NULL

	rownames(mu.gamma.0_conv) <- NULL
	rownames(mu.theta.0_conv) <- NULL
	rownames(tau2.theta.0_conv) <- NULL
	rownames(tau2.gamma.0_conv) <- NULL

	rownames(gamma_acc) <- NULL
	rownames(theta_acc) <- NULL

	if (model == "BB") {
		rownames(pi_conv) <- NULL
		rownames(alpha.pi_conv) <- NULL
		rownames(beta.pi_conv) <- NULL
		rownames(alpha.pi_acc) <- NULL
		rownames(beta.pi_acc) <- NULL
		rownames(theta_zero_prop) <- NULL
		rownames(theta_zero_acc) <- NULL
	}

	conv.diag = list(sim_type = raw$sim_type, type = type, gamma.conv.diag = gamma_conv,
								theta.conv.diag = theta_conv,
                                mu.gamma.conv.diag = mu.gamma_conv,
                                mu.theta.conv.diag = mu.theta_conv,
                                sigma2.gamma.conv.diag = sigma2.gamma_conv,
                                sigma2.theta.conv.diag = sigma2.theta_conv,
                                mu.gamma.0.conv.diag = mu.gamma.0_conv,
                                mu.theta.0.conv.diag = mu.theta.0_conv,
                                tau2.gamma.0.conv.diag = tau2.gamma.0_conv,
                                tau2.theta.0.conv.diag = tau2.theta.0_conv,
                                gamma_acc = gamma_acc,
                                theta_acc = theta_acc)

	if (model == "BB") {
		conv.diag$pi.conv.diag = pi_conv
		conv.diag$alpha.pi.conv.diag = alpha.pi_conv
		conv.diag$beta.pi.conv.diag = beta.pi_conv

		conv.diag$alpha.pi_acc = alpha.pi_acc
		conv.diag$beta.pi_acc = beta.pi_acc

		conv.diag$theta_zero_prop = theta_zero_prop
		conv.diag$theta_zero_acc = theta_zero_acc
	}

	attr(conv.diag, "model") = attr(raw, "model")
	return(conv.diag)
}

c212.BB.print.convergence.summary <- function(conv) {

	if (is.null(conv)) {
		print("NULL conv data")
		return(NULL)
	}

	model = attr(conv, "model")
	if (is.null(model)) {
		print("Convergence model attribute missing")
		return(NULL)
	}

	n = c("gamma.conv.diag", "theta.conv.diag", "mu.gamma.conv.diag", "mu.theta.conv.diag", "sigma2.gamma.conv.diag",
			"sigma2.theta.conv.diag", "mu.gamma.0.conv.diag", "mu.theta.0.conv.diag", "tau2.gamma.0.conv.diag",
			"tau2.theta.0.conv.diag", "gamma_acc", "theta_acc")

	if (model == "BB") {
		n = c(n, c("pi.conv.diag", "alpha.pi.conv.diag", "beta.pi.conv.diag", "alpha.pi_acc", "beta.pi_acc", "theta_zero_prop",
			"theta_zero_acc"))
	}

	if (M_global$checkNames(n, conv)) {
		print("Missing names");
		return(NULL)
	}

	cat(sprintf("Summary Convergence Diagnostics:\n"))
	cat(sprintf("================================\n"))

	chk = 0
	text_out = ""
	if (conv$type == "Gelman-Rubin") {
		chk = 0
		text_out = "Gelman-Rubin diagnostic"
	}
	else {
		chk = 1
		text_out = "Geweke statistic"
	}

	cat(sprintf("theta:\n"))
	cat(sprintf("------\n"))
	M_global$EOTprintConvSummLev1(conv$theta.conv.diag, text_out, chk)

	cat(sprintf("gamma:\n"))
	cat(sprintf("------\n"))
	M_global$EOTprintConvSummLev1(conv$gamma.conv.diag, text_out, chk)

	cat(sprintf("mu.gamma:\n"))
	cat(sprintf("---------\n"))
	M_global$EOTprintConvSummLev2(conv$mu.gamma.conv.diag, text_out, chk)

	cat(sprintf("mu.theta:\n"))
	cat(sprintf("---------\n"))
	M_global$EOTprintConvSummLev2(conv$mu.theta.conv.diag, text_out, chk)

	cat(sprintf("sigma2.gamma:\n"))
	cat(sprintf("-------------\n"))
	M_global$EOTprintConvSummLev2(conv$sigma2.gamma.conv.diag, text_out, chk)

	cat(sprintf("sigma2.theta:\n"))
	cat(sprintf("-------------\n"))
	M_global$EOTprintConvSummLev2(conv$sigma2.theta.conv.diag, text_out, chk)

	if (model == "BB") {
		cat(sprintf("pi:\n"))
		cat(sprintf("---\n"))
		M_global$EOTprintConvSummLev2(conv$pi.conv.diag, text_out, chk)
	}

	cat(sprintf("mu.gamma.0:\n"))
	cat(sprintf("-----------\n"))
	M_global$EOTprintConvSummLev3(conv$mu.gamma.0.conv.diag, text_out, chk)

	cat(sprintf("mu.theta.0:\n"))
	cat(sprintf("-----------\n"))
	M_global$EOTprintConvSummLev3(conv$mu.theta.0.conv.diag, text_out, chk)

	cat(sprintf("tau2.gamma.0:\n"))
	cat(sprintf("-------------\n"))
	M_global$EOTprintConvSummLev3(conv$tau2.gamma.0.conv.diag, text_out, chk)

	cat(sprintf("tau2.theta.0:\n"))
	cat(sprintf("-------------\n"))
	M_global$EOTprintConvSummLev3(conv$tau2.theta.0.conv.diag, text_out, chk)

	if (model == "BB") {
		cat(sprintf("alpha.pi:\n"))
		cat(sprintf("---------\n"))
		M_global$EOTprintConvSummLev3(conv$alpha.pi.conv.diag, text_out, chk)

		cat(sprintf("beta.pi:\n"))
		cat(sprintf("--------\n"))
		M_global$EOTprintConvSummLev3(conv$beta.pi.conv.diag, text_out, chk)
	}

	if (conv$sim_type == "MH") {
		cat("\nSampling Acceptance Rates:\n")
		cat("==========================\n")
		cat("theta:\n")
		cat("------\n")
		#print(sprintf("Min: %0.6f, Max: %0.6f", min(conv$theta_acc[!is.na(conv$theta_acc)]),
		#										max(conv$theta_acc[!is.na(conv$theta_acc)])))
		print(sprintf("Min: %0.6f, Max: %0.6f", min(conv$theta_acc$rate),
												max(conv$theta_acc$rate)))
		if (model == "BB") {
			#print(sprintf("Zero Acceptance Min: %0.6f, Max: %0.6f",
			#			min(conv$theta_zero_acc[!is.na(conv$theta_zero_acc)]),
			#			max(conv$theta_zero_acc[!is.na(conv$theta_zero_acc)])))
			#print(sprintf("Zero Acceptance Min: %0.6f, Max: %0.6f",
			#			min(conv$theta_zero_acc$rate),
			#			max(conv$theta_zero_acc$rate)))
		}

		cat("gamma:\n")
		cat("------\n")
		print(sprintf("Min: %0.6f, Max: %0.6f", min(conv$gamma_acc$rate),
												max(conv$gamma_acc$rate)))

		if (model == "BB") {
			cat("alpha.pi:\n")
			cat("---------\n")
			print(sprintf("Min: %0.6f, Max: %0.6f", min(conv$alpha.pi_acc[!is.na(conv$alpha.pi_acc)]),
                                                max(conv$alpha.pi_acc[!is.na(conv$alpha.pi_acc)])))

			cat("beta.pi:\n")
			cat("--------\n")
			print(sprintf("Min: %0.6f, Max: %0.6f", min(conv$beta.pi_acc[!is.na(conv$beta.pi_acc)]),
                                                max(conv$beta.pi_acc[!is.na(conv$beta.pi_acc)])))
		}
	}
	else {
		if (model == "BB") {
			cat("\nSampling Acceptance Rates:\n")
			cat("==========================\n")
			cat("theta:\n")
			cat("------\n")
			print(sprintf("Min: %0.6f, Max: %0.6f", min(conv$theta_acc$rate),
												max(conv$theta_acc$rate)))
		}
	}
}
