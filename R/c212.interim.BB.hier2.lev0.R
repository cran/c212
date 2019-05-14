# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 05/06/2015


Mi_BB_h2_l0 <- new.env()

Mi_BB_h2_l0$Id <- "$Id: c212.interim.BB.hier2.lev0.R,v 1.9 2019/05/05 13:18:12 clb13102 Exp clb13102 $"

c212.interim.BB.hier2.lev0 <- function(trial.data, sim_type = "SLICE", burnin = 10000, iter = 60000, nchains = 5,
	theta_algorithm = "MH",
	global.sim.params = data.frame(type = c("MH", "MH", "MH", "MH", "SLICE", "SLICE", "SLICE"),
                            param = c("sigma_MH_alpha", "sigma_MH_beta", "sigma_MH_gamma", "sigma_MH_theta",
                            "w_alpha", "w_beta", "w_gamma"),
                            value = c(3, 3, 0.2, 0.5, 1, 1, 1), control = c(0, 0, 0, 0, 6, 6, 6),
							stringsAsFactors = FALSE),
	sim.params = NULL,
	monitor = data.frame(variable = c("theta", "gamma", "mu.gamma", "mu.theta",
		"sigma2.theta", "sigma2.gamma", "pi"),
		monitor = c(1, 1, 1, 1, 1, 1, 1), stringsAsFactors = FALSE),
	initial_values = NULL,
	hyper_params = list(mu.gamma.0 = 0, tau2.gamma.0 = 10,
	mu.theta.0 = 0, tau2.theta.0 = 10, alpha.gamma = 3, beta.gamma = 1, alpha.theta = 3,
	beta.theta = 1, alpha.pi = 1.1, beta.pi = 1.1),
	global.pm.weight = 0.5,
	pm.weights = NULL,
	adapt_phase=1, memory_model = "HIGH")
{
	interim = M_global$INTERIMdata(Mi_BB_h2_l0, trial.data, iter, nchains, burnin, initial_values)

	if (is.null(interim)) {
		return(NULL)
	}

	trial.data = interim$trial.data
	cntrl.data = interim$cntrl.data

	Mi_BB_h2_l0$Algo <- theta_algorithm

	Mi_BB_h2_l0$sim_type <- sim_type

	if (nrow(global.sim.params[global.sim.params$type == sim_type,]) == 0) {
		print("Missing simulation parametetrs");
		return(NULL)
	}

	if (!all(global.sim.params$value > 0)) {
		print("Invalid simulation parameter value");
		return(NULL)
	}

	Mi_BB_h2_l0$global.sim.params <- global.sim.params

	Mi_BB_h2_l0$Level = 0

	sp = M_global$INTERIM_sim_paramsBB_2(Mi_BB_h2_l0, sim.params, pm.weights, sim_type, trial.data, cntrl.data)

	sim.params = sp$sim.params
	pm.weights = sp$pm.weights

	monitor = M_global$INTERIM_monitor_BB_2(monitor)

	# Initialise the hyper-parameters
	Mi_BB_h2_l0$mu.gamma.0 <- hyper_params$mu.gamma.0
	Mi_BB_h2_l0$tau2.gamma.0 <- hyper_params$tau2.gamma.0
	Mi_BB_h2_l0$alpha.gamma <- hyper_params$alpha.gamma
	Mi_BB_h2_l0$beta.gamma <- hyper_params$beta.gamma

	Mi_BB_h2_l0$mu.theta.0 <- hyper_params$mu.theta.0
	Mi_BB_h2_l0$tau2.theta.0 <- hyper_params$tau2.theta.0
	Mi_BB_h2_l0$alpha.theta <- hyper_params$alpha.theta
	Mi_BB_h2_l0$beta.theta <- hyper_params$beta.theta

	Mi_BB_h2_l0$alpha.pi <- hyper_params$alpha.pi
	Mi_BB_h2_l0$beta.pi <- hyper_params$beta.pi

	algo = 1
	if (Mi_BB_h2_l0$Algo == "BB2004") {
		algo <- 1;
	} else {
		if (Mi_BB_h2_l0$Algo == "MH") {
			algo <- 2;
		} else {
			if (Mi_BB_h2_l0$Algo == "Adapt") {
				algo <- 3;
			} else {
				if (Mi_BB_h2_l0$Algo == "Indep") {
					algo <- 4;
				} else {
					algo <- 1;
				}
			}
		}
	}

	Ret2 = .Call("c212BB_interim_hier2_exec", as.integer(nchains), as.integer(burnin),
					as.integer(iter), Mi_BB_h2_l0$sim_type, memory_model,
					Mi_BB_h2_l0$global.sim.params,
					sim.params,
					as.numeric(global.pm.weight),
					pm.weights,
					monitor,
					as.integer(Mi_BB_h2_l0$numIntervals),
					as.integer(Mi_BB_h2_l0$Level),
					Mi_BB_h2_l0$maxBs,
					as.integer(Mi_BB_h2_l0$numB),
					as.integer(Mi_BB_h2_l0$maxAEs),
					as.integer(t(Mi_BB_h2_l0$nAE)),
					as.integer(aperm(Mi_BB_h2_l0$x)),
					as.integer(aperm(Mi_BB_h2_l0$y)),
					as.numeric(aperm(Mi_BB_h2_l0$C)),
					as.numeric(aperm(Mi_BB_h2_l0$T)),
					as.numeric(aperm(Mi_BB_h2_l0$theta)),
					as.numeric(aperm(Mi_BB_h2_l0$gamma)),
					as.numeric(Mi_BB_h2_l0$alpha.gamma),
					as.numeric(Mi_BB_h2_l0$beta.gamma),
					as.numeric(Mi_BB_h2_l0$alpha.theta),
					as.numeric(Mi_BB_h2_l0$beta.theta),
					as.numeric(Mi_BB_h2_l0$mu.gamma.0),
					as.numeric(Mi_BB_h2_l0$tau2.gamma.0),
					as.numeric(Mi_BB_h2_l0$mu.theta.0),
					as.numeric(Mi_BB_h2_l0$tau2.theta.0),
					as.numeric(aperm(Mi_BB_h2_l0$mu.gamma)),
					as.numeric(aperm(Mi_BB_h2_l0$mu.theta)),
					as.numeric(aperm(Mi_BB_h2_l0$sigma2.gamma)),
					as.numeric(aperm(Mi_BB_h2_l0$sigma2.theta)),
					as.numeric(aperm(Mi_BB_h2_l0$pi)),
					as.numeric(Mi_BB_h2_l0$alpha.pi),
					as.numeric(Mi_BB_h2_l0$beta.pi),
					as.integer(algo),
					as.integer(adapt_phase))

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

	model_fit = list(id = Mi_BB_h2_l0$Id, sim_type = Mi_BB_h2_l0$sim_type, chains = nchains, nIntervals = Mi_BB_h2_l0$numIntervals,
			Intervals = Mi_BB_h2_l0$Intervals, nBodySys = Mi_BB_h2_l0$numB, maxBs = Mi_BB_h2_l0$maxBs,
			maxAEs = Mi_BB_h2_l0$maxAEs, nAE = Mi_BB_h2_l0$nAE, AE=Mi_BB_h2_l0$AE, B = Mi_BB_h2_l0$B,
			burnin = burnin, iter = iter,
			monitor = monitor,
			gamma = gamma_samples,
			theta = theta_samples,
			mu.gamma = mu.gamma_samples,
			mu.theta = mu.theta_samples,
			sigma2.gamma = sigma2.gamma_samples,
			sigma2.theta = sigma2.theta_samples,
			pi = pi_samples,
			gamma_acc = gamma_acc,
			theta_acc = theta_acc)
			
	# Model is poisson with BB hierarchy and independent intervals
	attr(model_fit, "model") = "BB_pois_h2_l0"

	return(model_fit)
}

Mi_BB_h2_l0$initVars = function() {

    # Data Structure
    Mi_BB_h2_l0$B <- c()
    Mi_BB_h2_l0$numB <- NA
    Mi_BB_h2_l0$numIntervals <- NA
    Mi_BB_h2_l0$nAE <- c()
    Mi_BB_h2_l0$maxAEs <- NA

    # Trial Event Data
    Mi_BB_h2_l0$x <- array()
    Mi_BB_h2_l0$C <- array()
    Mi_BB_h2_l0$y <- array()
    Mi_BB_h2_l0$T <- array()

    # Hyperparameters
    Mi_BB_h2_l0$mu.gamma.0 <- NA
    Mi_BB_h2_l0$tau2.gamma.0 <- NA
    Mi_BB_h2_l0$mu.theta.0 <- NA
    Mi_BB_h2_l0$tau2.theta.0 <- NA
    Mi_BB_h2_l0$alpha.gamma <- NA
    Mi_BB_h2_l0$beta.gamma <- NA
    Mi_BB_h2_l0$alpha.theta <- NA
    Mi_BB_h2_l0$beta.theta <- NA

    # Stage 2
    Mi_BB_h2_l0$mu.gamma <- array()
    Mi_BB_h2_l0$mu.theta <- array()
    Mi_BB_h2_l0$sigma2.gamma <- array()
    Mi_BB_h2_l0$sigma2.theta <- array()

    # Stage 1
    Mi_BB_h2_l0$theta <- array()
    Mi_BB_h2_l0$gamma <- array()

	# BB2004 parameters
	Mi_BB_h2_l0$alpha.pi <- NA
	Mi_BB_h2_l0$beta.pi <- NA
	Mi_BB_h2_l0$pi <- NA
}

Mi_BB_h2_l0$initChains = function(c) {
	# Choose random values for gamma and theta
	for (i in 1:Mi_BB_h2_l0$numIntervals) {
		numB = Mi_BB_h2_l0$numB[i]
		for (b in 1:numB) {
			Mi_BB_h2_l0$gamma[c, i, b, 1:Mi_BB_h2_l0$nAE[i, b]] <- runif(Mi_BB_h2_l0$nAE[i, b], -10, 10)
			Mi_BB_h2_l0$theta[c, i, b, 1:Mi_BB_h2_l0$nAE[i, b]] <- runif(Mi_BB_h2_l0$nAE[i, b], -10, 10)

			Mi_BB_h2_l0$theta[c, i, b, ][is.infinite(Mi_BB_h2_l0$theta[c, i, b, ])] = -10
			Mi_BB_h2_l0$gamma[c, i, b, ][is.infinite(Mi_BB_h2_l0$gamma[c, i, b, ])] = -10

			Mi_BB_h2_l0$theta[c, i, b, ][is.nan(Mi_BB_h2_l0$theta[c, i, b, ])] = -10 # -1000
			Mi_BB_h2_l0$gamma[c, i, b, ][is.nan(Mi_BB_h2_l0$gamma[c, i, b, ])] = -10 # -1000
		}

		Mi_BB_h2_l0$mu.gamma[c, i, 1:numB] = runif(numB, -10, 10)
		Mi_BB_h2_l0$mu.theta[c, i, 1:numB] = runif(numB, -10, 10)
		Mi_BB_h2_l0$sigma2.gamma[c, i, 1:numB] = runif(numB, 5, 20)
		Mi_BB_h2_l0$sigma2.theta[c, i, 1:numB] = runif(numB, 5, 20)

		Mi_BB_h2_l0$pi[c, i, 1:numB] = runif(numB, 0, 1)
	}
}

Mi_BB_h2_l0$initialiseChains = function(initial_values, nchains) {

	Mi_BB_h2_l0$theta = array(0, dim=c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs, Mi_BB_h2_l0$maxAEs))
	Mi_BB_h2_l0$gamma = array(0, dim=c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs, Mi_BB_h2_l0$maxAEs))

	if (is.null(initial_values)) {

		# Initialise the first chain with the data
		for (i in 1:Mi_BB_h2_l0$numIntervals) {
			numB = Mi_BB_h2_l0$numB[i]
			for (b in 1:numB) {
				Mi_BB_h2_l0$gamma[1, i, b, ] <- log(Mi_BB_h2_l0$x[i, b,]/Mi_BB_h2_l0$C[i, b, ])
				Mi_BB_h2_l0$theta[1, i, b, ] <- log(Mi_BB_h2_l0$y[i, b,]/Mi_BB_h2_l0$T[i, b, ]) - Mi_BB_h2_l0$gamma[1, i, b, ]

				Mi_BB_h2_l0$theta[1, i, b, ][is.infinite(Mi_BB_h2_l0$theta[1, i, b, ])] = -10 # -1000
				Mi_BB_h2_l0$gamma[1, i, b, ][is.infinite(Mi_BB_h2_l0$gamma[1, i, b, ])] = -10 # -1000

				Mi_BB_h2_l0$theta[1, i, b, ][is.nan(Mi_BB_h2_l0$theta[1, i, b, ])] = -10 # -1000
				Mi_BB_h2_l0$gamma[1, i, b, ][is.nan(Mi_BB_h2_l0$gamma[1, i, b, ])] = -10 # -1000
			}
		}

		Mi_BB_h2_l0$mu.gamma <- array(0, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))
		Mi_BB_h2_l0$mu.theta <- array(0, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))
		Mi_BB_h2_l0$sigma2.gamma <- array(10, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))
		Mi_BB_h2_l0$sigma2.theta <- array(10, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))

		Mi_BB_h2_l0$pi <- array(0.5, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))

		if (nchains > 1) {
			for (c in 2:nchains) {
				Mi_BB_h2_l0$initChains(c)
			}
		}
	}
	else {

		Mi_BB_h2_l0$mu.gamma <- array(0, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))
		Mi_BB_h2_l0$mu.theta <- array(0, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))
		Mi_BB_h2_l0$sigma2.gamma <- array(0, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))
		Mi_BB_h2_l0$sigma2.theta <- array(0, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))

		Mi_BB_h2_l0$pi <- array(0.5, dim = c(nchains, Mi_BB_h2_l0$numIntervals, Mi_BB_h2_l0$maxBs))

		for (c in 1:nchains) {
			for (i in 1:Mi_BB_h2_l0$numIntervals) {
				interval = Mi_BB_h2_l0$Intervals[i]
				for (b in 1:Mi_BB_h2_l0$numB[i]) {
					data = initial_values$mu.gamma[initial_values$mu.gamma$chain == c &
												initial_values$mu.gamma$Interval == interval
												& initial_values$mu.gamma$B == Mi_BB_h2_l0$B[1, b],]
					Mi_BB_h2_l0$mu.gamma[c, i, b] = data$value

					data = initial_values$mu.theta[initial_values$mu.theta$chain == c &
												initial_values$mu.theta$Interval == interval
												& initial_values$mu.theta$B == Mi_BB_h2_l0$B[1, b],]
					Mi_BB_h2_l0$mu.theta[c, i, b] = data$value

					data = initial_values$sigma2.gamma[initial_values$sigma2.gamma$chain == c &
												initial_values$sigma2.gamma$Interval == interval
												& initial_values$sigma2.gamma$B == Mi_BB_h2_l0$B[1, b],]
					Mi_BB_h2_l0$sigma2.gamma[c, i, b] = data$value

					data = initial_values$sigma2.theta[initial_values$sigma2.theta$chain == c &
												initial_values$sigma2.theta$Interval == interval
												& initial_values$sigma2.theta$B == Mi_BB_h2_l0$B[1, b],]
					Mi_BB_h2_l0$sigma2.theta[c, i, b] = data$value

					data = initial_values$pi[initial_values$pi$chain == c &
												initial_values$pi$Interval == interval
												& initial_values$pi$B == Mi_BB_h2_l0$B[1, b],]
					Mi_BB_h2_l0$pi[c, i, b] = data$value
				}
			}
		}

		for (c in 1:nchains) {
			for (i in 1:Mi_BB_h2_l0$numIntervals) {
				interval = Mi_BB_h2_l0$Intervals[i]
				for (b in 1:Mi_BB_h2_l0$numB[i]) {
					for (j in 1:Mi_BB_h2_l0$nAE[i, b]) {
						ae = Mi_BB_h2_l0$AE[i, b, j]
						data = initial_values$gamma[initial_values$gamma$chain == c
										& initial_values$gamma$Interval == interval
										& initial_values$gamma$B == Mi_BB_h2_l0$B[i, b]
										& initial_values$gamma$AE == ae,]
						Mi_BB_h2_l0$gamma[c, i, b, j] = data$value

						data = initial_values$theta[initial_values$theta$chain == c
										& initial_values$theta$Interval == interval
										& initial_values$theta$B == Mi_BB_h2_l0$B[i, b]
										& initial_values$theta$AE == ae,]
						Mi_BB_h2_l0$theta[c, i, b, j] = data$value
					}
				}
			}
		}
	}
}
