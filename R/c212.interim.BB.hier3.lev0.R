# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 05/06/2015


Mi1 <- new.env()

Mi1$Id <- "$Id: c212.interim.BB.hier3.lev0.R,v 1.13 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.interim.BB.indep <- function(trial.data, sim_type = "SLICE", burnin = 10000, iter = 60000, nchains = 5,
	theta_algorithm = "MH",
	global.sim.params = data.frame(type = c("MH", "MH", "MH", "MH", "SLICE", "SLICE", "SLICE"),
                            param = c("sigma_MH_alpha", "sigma_MH_beta", "sigma_MH_gamma", "sigma_MH_theta",
                            "w_alpha", "w_beta", "w_gamma"),
                            value = c(3, 3, 0.2, 0.25, 1, 1, 1), control = c(0, 0, 0, 0, 6, 6, 6),
							stringsAsFactors = FALSE),
	sim.params = NULL,
	monitor = data.frame(variable = c("theta", "gamma", "mu.gamma", "mu.theta",
					"sigma2.theta", "sigma2.gamma",
		            "mu.theta.0", "mu.gamma.0", "tau2.theta.0", "tau2.gamma.0",
					"pi", "alpha.pi", "beta.pi"),
					monitor = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
					stringsAsFactors = FALSE),
	initial_values = NULL,
	hyper_params = list(mu.gamma.0.0 = 0, tau2.gamma.0.0 = 10,
	mu.theta.0.0 = 0, tau2.theta.0.0 = 10, alpha.gamma.0.0 = 3, beta.gamma.0.0 = 1, alpha.theta.0.0 = 3,
	beta.theta.0.0 = 1, alpha.gamma = 3, beta.gamma = 1, alpha.theta = 3, beta.theta = 1, lambda.alpha = 1.0,
	lambda.beta = 1.0),
	global.pm.weight = 0.5,
	pm.weights = NULL,
	adapt_phase=1, memory_model = "HIGH")
{
	interim = M_global$INTERIMdata(Mi1, trial.data, iter, nchains, burnin, initial_values)

	if (is.null(interim)) {
		return(NULL)
	}

	trial.data = interim$trial.data
	cntrl.data = interim$cntrl.data

	Mi1$Algo <- theta_algorithm

	Mi1$sim_type <- sim_type

	if (nrow(global.sim.params[global.sim.params$type == sim_type,]) == 0) {
		print("Missing simulation parametetrs");
		return(NULL)
	}

	if (!all(global.sim.params$value > 0)) {
		print("Invalid simulation parameter value");
		return(NULL)
	}

	Mi1$global.sim.params <- global.sim.params

	Mi1$Level = 0

	sp = M_global$INTERIM_sim_paramsBB_3(Mi1, sim.params, pm.weights, sim_type, trial.data, cntrl.data)

	sim.params = sp$sim.params
	pm.weights = sp$pm.weights

	monitor = M_global$INTERIM_monitor_BB_3(monitor)

	# Initialise the hyper-parameters
	Mi1$mu.gamma.0.0 <- hyper_params$mu.gamma.0.0
	Mi1$tau2.gamma.0.0 <- hyper_params$tau2.gamma.0.0
	Mi1$alpha.gamma <- hyper_params$alpha.gamma
	Mi1$beta.gamma <- hyper_params$beta.gamma
	Mi1$alpha.gamma.0.0 <- hyper_params$alpha.gamma.0.0
	Mi1$beta.gamma.0.0 <- hyper_params$beta.gamma.0.0

	Mi1$mu.theta.0.0 <- hyper_params$mu.theta.0.0
	Mi1$tau2.theta.0.0 <- hyper_params$tau2.theta.0.0
	Mi1$alpha.theta <- hyper_params$alpha.theta
	Mi1$beta.theta <- hyper_params$beta.theta
	Mi1$alpha.theta.0.0 <- hyper_params$alpha.theta.0.0
	Mi1$beta.theta.0.0 <- hyper_params$beta.theta.0.0

	# BB2004 parameters
	Mi1$lambda.alpha <- hyper_params$lambda.alpha
	Mi1$lambda.beta <- hyper_params$lambda.beta

	algo = 1
	if (Mi1$Algo == "BB2004") {
		algo <- 1;
	} else {
		if (Mi1$Algo == "MH") {
			algo <- 2;
		} else {
			if (Mi1$Algo == "Adapt") {
				algo <- 3;
			} else {
				if (Mi1$Algo == "Indep") {
					algo <- 4;
				} else {
					algo <- 1;
				}
			}
		}
	}


	Ret2 = .Call("c212BB_poisson_mc_exec", as.integer(nchains), as.integer(burnin),
					as.integer(iter), Mi1$sim_type,
					memory_model, Mi1$global.sim.params,
					sim.params,
					as.numeric(global.pm.weight),
					pm.weights,
					monitor,
					as.integer(Mi1$numIntervals), as.integer(Mi1$Level),
					Mi1$maxBs, as.integer(Mi1$numB), as.integer(Mi1$maxAEs),
					as.integer(t(Mi1$nAE)), as.integer(aperm(Mi1$x)), as.integer(aperm(Mi1$y)),
					as.integer(aperm(Mi1$C)),
					as.integer(aperm(Mi1$T)),
					as.numeric(aperm(Mi1$theta)),
					as.numeric(aperm(Mi1$gamma)),
					as.numeric(Mi1$mu.gamma.0.0),
					as.numeric(Mi1$tau2.gamma.0.0),
					as.numeric(Mi1$mu.theta.0.0),
					as.numeric(Mi1$tau2.theta.0.0),
					as.numeric(Mi1$alpha.gamma.0.0),
					as.numeric(Mi1$beta.gamma.0.0),
					as.numeric(Mi1$alpha.theta.0.0),
					as.numeric(Mi1$beta.theta.0.0),
					as.numeric(Mi1$alpha.gamma),
					as.numeric(Mi1$beta.gamma),
					as.numeric(Mi1$alpha.theta),
					as.numeric(Mi1$beta.theta),
					as.numeric(aperm(Mi1$mu.gamma.0)),
					as.numeric(aperm(Mi1$tau2.gamma.0)),
					as.numeric(aperm(Mi1$mu.theta.0)),
					as.numeric(aperm(Mi1$tau2.theta.0)),
					as.numeric(aperm(Mi1$mu.gamma)),
					as.numeric(aperm(Mi1$mu.theta)),
					as.numeric(aperm(Mi1$sigma2.gamma)),
					as.numeric(aperm(Mi1$sigma2.theta)),
					as.numeric(aperm(Mi1$pi)),
					as.numeric(aperm(Mi1$alpha.pi)),
					as.numeric(aperm(Mi1$beta.pi)),
					as.numeric(Mi1$lambda.alpha),
					as.numeric(Mi1$lambda.beta),
					as.integer(algo),
					as.integer(adapt_phase))

	mu.gamma.0_samples = NULL
	if (monitor[monitor$variable == "mu.gamma.0", ]$monitor == 1) {
		mu.gamma.0_samples <- .Call("getMuGamma0SamplesInterimAll")
		mu.gamma.0_samples = aperm(mu.gamma.0_samples)
	}

	mu.theta.0_samples = NULL
	if (monitor[monitor$variable == "mu.theta.0", ]$monitor == 1) {
		mu.theta.0_samples <- .Call("getMuTheta0SamplesInterimAll")
		mu.theta.0_samples = aperm(mu.theta.0_samples)
	}

	tau2.gamma.0_samples = NULL
	if (monitor[monitor$variable == "tau2.gamma.0", ]$monitor == 1) {
		tau2.gamma.0_samples <- .Call("getTau2Gamma0SamplesInterimAll")
		tau2.gamma.0_samples = aperm(tau2.gamma.0_samples)
	}

	tau2.theta.0_samples = NULL
	if (monitor[monitor$variable == "tau2.theta.0", ]$monitor == 1) {
		tau2.theta.0_samples <- .Call("getTau2Theta0SamplesInterimAll")
		tau2.theta.0_samples = aperm(tau2.theta.0_samples)
	}

	mu.theta_samples = NULL
	if (monitor[monitor$variable == "mu.theta", ]$monitor == 1) {
		mu.theta_samples <- .Call("getMuThetaSamplesInterimAll")
		mu.theta_samples <- aperm(mu.theta_samples)
	}

	mu.gamma_samples = NULL
	if (monitor[monitor$variable == "mu.gamma", ]$monitor == 1) {
		mu.gamma_samples <- .Call("getMuGammaSamplesInterimAll")
		mu.gamma_samples <- aperm(mu.gamma_samples)
	}

	sigma2.theta_samples = NULL
	if (monitor[monitor$variable == "sigma2.theta", ]$monitor == 1) {
		sigma2.theta_samples <- .Call("getSigma2ThetaSamplesInterimAll")
		sigma2.theta_samples <- aperm(sigma2.theta_samples)
	}

	sigma2.gamma_samples = NULL
	if (monitor[monitor$variable == "sigma2.gamma", ]$monitor == 1) {
		sigma2.gamma_samples <- .Call("getSigma2GammaSamplesInterimAll")
		sigma2.gamma_samples <- aperm(sigma2.gamma_samples)
	}

	pi_samples = NULL
	if (monitor[monitor$variable == "pi", ]$monitor == 1) {
		pi_samples = .Call("getPiSamplesInterimAll")
		pi_samples <- aperm(pi_samples)
	}

	alpha.pi_samples = NULL
	alpha.pi_acc = NULL
	if (monitor[monitor$variable == "alpha.pi", ]$monitor == 1) {
		alpha.pi_samples = .Call("getAlphaPiSamplesInterimAll")
		alpha.pi_samples = aperm(alpha.pi_samples)

		alpha.pi_acc = .Call("getAlphaPiAcceptInterimAll")
		alpha.pi_acc = aperm(alpha.pi_acc)
	}

	beta.pi_samples = NULL
	beta.pi_acc = NULL
	if (monitor[monitor$variable == "beta.pi", ]$monitor == 1) {
		beta.pi_samples = .Call("getBetaPiSamplesInterimAll")
		beta.pi_samples = aperm(beta.pi_samples)

		beta.pi_acc = .Call("getBetaPiAcceptInterimAll")
		beta.pi_acc = aperm(beta.pi_acc)
	}

	gamma_samples = NULL
	gamma_acc = NULL
	if (monitor[monitor$variable == "gamma", ]$monitor == 1) {
		gamma_samples = .Call("getGammaSamplesInterimAll")
		gamma_samples = aperm(gamma_samples)

		gamma_acc = .Call("getGammaAcceptInterimAll")
		gamma_acc <- aperm(gamma_acc)
	}

	theta_samples = NULL
	theta_acc = NULL
	if (monitor[monitor$variable == "theta", ]$monitor == 1) {
		theta_samples = .Call("getThetaSamplesInterimAll")
		theta_samples = aperm(theta_samples)

		theta_acc = .Call("getThetaAcceptInterimAll")
		theta_acc <- aperm(theta_acc)
	}

	.C("Release_Interim")

	model_fit = list(id = Mi1$Id, sim_type = Mi1$sim_type, chains = nchains, nIntervals = Mi1$numIntervals,
			Intervals = Mi1$Intervals, nBodySys = Mi1$numB, maxBs = Mi1$maxBs,
			maxAEs = Mi1$maxAEs, nAE = Mi1$nAE, AE=Mi1$AE, B = Mi1$B,
			burnin = burnin, iter = iter,
			monitor = monitor,
			gamma = gamma_samples,
			theta = theta_samples,
			mu.gamma = mu.gamma_samples,
			mu.theta = mu.theta_samples,
			sigma2.gamma = sigma2.gamma_samples,
			sigma2.theta = sigma2.theta_samples,
			pi = pi_samples,
			alpha.pi = alpha.pi_samples,
			beta.pi = beta.pi_samples,
			alpha.pi_acc = alpha.pi_acc,
			beta.pi_acc = beta.pi_acc,
			mu.gamma.0 = mu.gamma.0_samples,
			mu.theta.0 = mu.theta.0_samples,
			tau2.gamma.0 = tau2.gamma.0_samples,
			tau2.theta.0 = tau2.theta.0_samples,
			gamma_acc = gamma_acc,
			theta_acc = theta_acc)
			
	# Model is poisson with BB hierarchy and independent intervals
	attr(model_fit, "model") = "BB_pois_indep"

	return(model_fit)
}

Mi1$initVars = function() {

    # Data Structure
    Mi1$B <- c()
    Mi1$numB <- NA
    Mi1$numIntervals <- NA
    Mi1$nAE <- c()
    Mi1$maxAEs <- NA

    # Trial Event Data
    Mi1$x <- array()
    Mi1$C <- array()
    Mi1$y <- array()
    Mi1$T <- array()

    # Hyperparameters
    Mi1$mu.gamma.0.0 <- NA
    Mi1$tau2.gamma.0.0 <- NA
    Mi1$mu.theta.0.0 <- NA
    Mi1$tau2.theta.0.0 <- NA
    Mi1$alpha.gamma.0.0 <- NA
    Mi1$beta.gamma.0.0 <- NA
    Mi1$alpha.theta.0.0 <- NA
    Mi1$beta.theta.0.0 <- NA
    Mi1$alpha.gamma <- NA
    Mi1$beta.gamma <- NA
    Mi1$alpha.theta <- NA
    Mi1$beta.theta <- NA

	# Parameters/Simulated values
    # Stage 3
    Mi1$mu.gamma.0 <- c()
    Mi1$tau2.gamma.0 <- c()
    Mi1$mu.theta.0 <- c()
    Mi1$tau2.theta.0 <- c()

    # Stage 2
    Mi1$mu.gamma <- array()
    Mi1$mu.theta <- array()
    Mi1$sigma2.gamma <- array()
    Mi1$sigma2.theta <- array()

    # Stage 1
    Mi1$theta <- array()
    Mi1$gamma <- array()

    # Acceptance rates

	# BB2004 parameters
	Mi1$lambda.alpha <- NA
	Mi1$lambda.beta <- NA
	Mi1$alpha.pi <- NA
	Mi1$beta.pi <- NA
	Mi1$pi <- NA
}

Mi1$initChains = function(c) {
	# Choose random values for gamma and theta
	for (i in 1:Mi1$numIntervals) {
		numB = Mi1$numB[i]
		for (b in 1:numB) {
			Mi1$gamma[c, i, b, 1:Mi1$nAE[i, b]] <- runif(Mi1$nAE[i, b], -10, 10)
			Mi1$theta[c, i, b, 1:Mi1$nAE[i, b]] <- runif(Mi1$nAE[i, b], -10, 10)

			Mi1$theta[c, i, b, ][is.infinite(Mi1$theta[c, i, b, ])] = -10
			Mi1$gamma[c, i, b, ][is.infinite(Mi1$gamma[c, i, b, ])] = -10

			Mi1$theta[1, i, b, ][is.nan(Mi1$theta[1, i, b, ])] = -10 # -1000
			Mi1$gamma[1, i, b, ][is.nan(Mi1$gamma[1, i, b, ])] = -10 # -1000
		}

		Mi1$mu.gamma[c, i, 1:numB] = runif(numB, -10, 10)
		Mi1$mu.theta[c, i, 1:numB] = runif(numB, -10, 10)
		Mi1$sigma2.gamma[c, i, 1:numB] = runif(numB, 5, 20)
		Mi1$sigma2.theta[c, i, 1:numB] = runif(numB, 5, 20)

		Mi1$pi[c, i, 1:numB] = runif(numB, 0, 1)

		Mi1$mu.gamma.0[c, i] = runif(1, -10, 10)
		Mi1$tau2.gamma.0[c, i] = runif(1, 5, 20)
		Mi1$mu.theta.0[c, i] = runif(1, -10, 10)
		Mi1$tau2.theta.0[c, i] = runif(1, 5, 20)

		Mi1$alpha.pi[c, i] = runif(1, 1.25, 100)
		Mi1$beta.pi[c, i] = runif(1, 1.25, 100)
	}
}

Mi1$initialiseChains = function(initial_values, nchains) {

	Mi1$theta = array(0, dim=c(nchains, Mi1$numIntervals, Mi1$maxBs, Mi1$maxAEs))
	Mi1$gamma = array(0, dim=c(nchains, Mi1$numIntervals, Mi1$maxBs, Mi1$maxAEs))

	if (is.null(initial_values)) {

		# Initialise the first chain with the data
		for (i in 1:Mi1$numIntervals) {
			numB = Mi1$numB[i]
			for (b in 1:numB) {
				Mi1$gamma[1, i, b, ] <- log(Mi1$x[i, b,]/Mi1$C[i, b, ])
				Mi1$theta[1, i, b, ] <- log(Mi1$y[i, b,]/Mi1$T[i]) - Mi1$gamma[1, i, b, ]

				Mi1$theta[1, i, b, ][is.infinite(Mi1$theta[1, i, b, ])] = -10 # -1000
				Mi1$gamma[1, i, b, ][is.infinite(Mi1$gamma[1, i, b, ])] = -10 # -1000

				Mi1$theta[1, i, b, ][is.nan(Mi1$theta[1, i, b, ])] = -10 # -1000
				Mi1$gamma[1, i, b, ][is.nan(Mi1$gamma[1, i, b, ])] = -10 # -1000
			}
		}

		Mi1$mu.gamma <- array(0, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))
		Mi1$mu.theta <- array(0, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))
		Mi1$sigma2.gamma <- array(10, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))
		Mi1$sigma2.theta <- array(10, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))

		Mi1$pi <- array(0.5, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))

		Mi1$mu.gamma.0 <- array(0, dim = c(nchains, Mi1$numIntervals))
		Mi1$tau2.gamma.0 <- array(10, dim = c(nchains, Mi1$numIntervals))
		Mi1$mu.theta.0 <- array(0, dim = c(nchains, Mi1$numIntervals))
		Mi1$tau2.theta.0 <- array(10, dim = c(nchains, Mi1$numIntervals))

		Mi1$alpha.pi <- array(1.5, dim = c(nchains, Mi1$numIntervals))
		Mi1$beta.pi <- array(1.5, dim = c(nchains, Mi1$numIntervals))

		if (nchains > 1) {
			for (c in 2:nchains) {
				Mi1$initChains(c)
			}
		}
	}
	else {

		Mi1$mu.gamma.0 <- array(0, dim = c(nchains, Mi1$numIntervals))
		Mi1$tau2.gamma.0 <- array(10, dim = c(nchains, Mi1$numIntervals))
		Mi1$mu.theta.0 <- array(0, dim = c(nchains, Mi1$numIntervals))
		Mi1$tau2.theta.0 <- array(10, dim = c(nchains, Mi1$numIntervals))

		Mi1$alpha.pi <- array(10, dim = c(nchains, Mi1$numIntervals))
		Mi1$beta.pi <- array(10, dim = c(nchains, Mi1$numIntervals))

		for (c in 1:nchains) {
			for (i in 1:Mi1$numIntervals) {
				interval = Mi1$Intervals[i]
				data = initial_values$mu.gamma.0[initial_values$mu.gamma.0$chain == c &
												initial_values$mu.gamma.0$Interval == interval, ]
				Mi1$mu.gamma.0[c, i] = data$value
				data = initial_values$mu.theta.0[initial_values$mu.theta.0$chain == c &
												initial_values$mu.theta.0$Interval == interval, ]
				Mi1$mu.theta.0[c, i] = data$value

				data = initial_values$tau2.gamma.0[initial_values$tau2.gamma.0$chain == c &
												initial_values$tau2.gamma.0$Interval == interval, ]
				Mi1$tau2.gamma.0[c, i] = data$value
				data = initial_values$tau2.theta.0[initial_values$tau2.theta.0$chain == c &
												initial_values$tau2.theta.0$Interval == interval, ]
				Mi1$tau2.theta.0[c, i] = data$value

				data = initial_values$alpha.pi[initial_values$alpha.pi$chain == c &
                                                initial_values$alpha.pi$Interval == interval, ]
				Mi1$alpha.pi[c, i] = data$value
				data = initial_values$beta.pi[initial_values$beta.pi$chain == c &
                                                initial_values$beta.pi$Interval == interval, ]
				Mi1$beta.pi[c, i] = data$value
			}
		}

		Mi1$mu.gamma <- array(0, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))
		Mi1$mu.theta <- array(0, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))
		Mi1$sigma2.gamma <- array(0, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))
		Mi1$sigma2.theta <- array(0, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))

		Mi1$pi <- array(0.5, dim = c(nchains, Mi1$numIntervals, Mi1$maxBs))

		for (c in 1:nchains) {
			for (i in 1:Mi1$numIntervals) {
				interval = Mi1$Intervals[i]
				for (b in 1:Mi1$numB[i]) {
					data = initial_values$mu.gamma[initial_values$mu.gamma$chain == c &
									initial_values$mu.gamma$Interval == interval
												& initial_values$mu.gamma$B == Mi1$B[i, b],]
					Mi1$mu.gamma[c, i, b] = data$value

					data = initial_values$mu.theta[initial_values$mu.theta$chain == c &
									initial_values$mu.theta$Interval == interval
												& initial_values$mu.theta$B == Mi1$B[i, b],]
					Mi1$mu.theta[c, i, b] = data$value

					data = initial_values$sigma2.gamma[initial_values$sigma2.gamma$chain == c &
									initial_values$sigma2.gamma$Interval == interval
												& initial_values$sigma2.gamma$B == Mi1$B[i, b],]
					Mi1$sigma2.gamma[c, i, b] = data$value

					data = initial_values$sigma2.theta[initial_values$sigma2.theta$chain == c &
									initial_values$sigma2.theta$Interval == interval
												& initial_values$sigma2.theta$B == Mi1$B[i, b],]
					Mi1$sigma2.theta[c, i, b] = data$value

					data = initial_values$pi[initial_values$pi$chain == c &
									initial_values$pi$Interval == interval
												& initial_values$pi$B == Mi1$B[i, b],]
					Mi1$pi[c, i, b] = data$value
				}
			}
		}

		for (c in 1:nchains) {
			for (i in 1:Mi1$numIntervals) {
				interval = Mi1$Intervals[i]
				for (b in 1:Mi1$numB[i]) {
					for (j in 1:Mi1$nAE[i, b]) {
						ae = Mi1$AE[i, b, j]
						data = initial_values$gamma[initial_values$gamma$chain == c
										& initial_values$gamma$Interval == interval
										& initial_values$gamma$B == Mi1$B[i, b]
										& initial_values$gamma$AE == ae,]
						Mi1$gamma[c, i, b, j] = data$value

						data = initial_values$theta[initial_values$theta$chain == c
										& initial_values$theta$Interval == interval
										& initial_values$theta$B == Mi1$B[i, b]
										& initial_values$theta$AE == ae,]
						Mi1$theta[c, i, b, j] = data$value
					}
				}
			}
		}
	}
}
