# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 28/04/2015


Mi <- new.env()

Mi$Id <- "$Id: c212.interim.1a.hier3.lev0.R,v 1.12 2018/10/03 15:40:55 clb13102 Exp clb13102 $"

c212.interim.1a.indep <- function(trial.data, sim_type = "SLICE", burnin = 10000, iter = 40000, nchains = 3,
	global.sim.params = data.frame(type = c("MH", "SLICE"), param = c("sigma_MH", "w"), value = c(0.2,1),
	control = c(0,6)),
	sim.params = NULL,
	monitor = data.frame(variable = c("theta", "gamma", "mu.gamma", "mu.theta",
					"sigma2.theta", "sigma2.gamma",
		            "mu.theta.0", "mu.gamma.0", "tau2.theta.0", "tau2.gamma.0"),
					monitor = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
					stringsAsFactors = FALSE),
	initial_values = NULL,
	hyper_params = list(mu.gamma.0.0 = 0, tau2.gamma.0.0 = 10,
	mu.theta.0.0 = 0, tau2.theta.0.0 = 10, alpha.gamma.0.0 = 3, beta.gamma.0.0 = 1, alpha.theta.0.0 = 3,
	beta.theta.0.0 = 1, alpha.gamma = 3, beta.gamma = 1, alpha.theta = 3, beta.theta = 1), memory_model = "HIGH")
{

	interim = M_global$INTERIMdata(Mi, trial.data, iter, nchains, burnin, initial_values)

	if (is.null(interim)) {
		return(NULL)
	}

	trial.data = interim$trial.data
	cntrl.data = interim$cntrl.data

	Mi$sim_type <- sim_type

	if (nrow(global.sim.params[global.sim.params$type == sim_type,]) != 1) {
		print("Missing simulation parametetrs");
		return(NULL)
	}

	Mi$global.sim.param <- global.sim.params[global.sim.params$type == sim_type,]$value
	Mi$global.sim.param_ctrl <- global.sim.params[global.sim.params$type == sim_type,]$control

	if (Mi$global.sim.param <= 0) {
		print("Invalid simulation parametetr value");
		return(NULL)
	}

	Mi$level = 0

	sim.params = M_global$INTERIM_sim_params1a(Mi, sim.params, sim_type, trial.data, cntrl.data)

	monitor = M_global$INTERIM_monitor_1a_3(monitor)

	# Initialise the hyper-parameters
	Mi$mu.gamma.0.0 <- hyper_params$mu.gamma.0.0
	Mi$tau2.gamma.0.0 <- hyper_params$tau2.gamma.0.0
	Mi$alpha.gamma <- hyper_params$alpha.gamma
	Mi$beta.gamma <- hyper_params$beta.gamma
	Mi$alpha.gamma.0.0 <- hyper_params$alpha.gamma.0.0
	Mi$beta.gamma.0.0 <- hyper_params$beta.gamma.0.0

	Mi$mu.theta.0.0 <- hyper_params$mu.theta.0.0
	Mi$tau2.theta.0.0 <- hyper_params$tau2.theta.0.0
	Mi$alpha.theta <- hyper_params$alpha.theta
	Mi$beta.theta <- hyper_params$beta.theta
	Mi$alpha.theta.0.0 <- hyper_params$alpha.theta.0.0
	Mi$beta.theta.0.0 <- hyper_params$beta.theta.0.0

	Ret2 = .Call("c2121a_poisson_mc_exec", as.integer(nchains), as.integer(burnin),
					as.integer(iter), Mi$sim_type,
					memory_model,
					as.numeric(Mi$global.sim.param),
					as.numeric(Mi$global.sim.param_ctrl),
					sim.params,
					monitor,
					as.integer(Mi$numIntervals), as.integer(Mi$level),
					Mi$maxBs, as.integer(Mi$numB), as.integer(Mi$maxAEs),
					as.integer(t(Mi$nAE)), as.integer(aperm(Mi$x)),
					as.integer(aperm(Mi$y)),
					as.numeric(aperm(Mi$C)),
					as.numeric(aperm(Mi$T)),
					as.numeric(aperm(Mi$theta)),
					as.numeric(aperm(Mi$gamma)),
					as.numeric(Mi$mu.gamma.0.0),
					as.numeric(Mi$tau2.gamma.0.0),
					as.numeric(Mi$mu.theta.0.0),
					as.numeric(Mi$tau2.theta.0.0),
					as.numeric(Mi$alpha.gamma.0.0),
					as.numeric(Mi$beta.gamma.0.0),
					as.numeric(Mi$alpha.theta.0.0),
					as.numeric(Mi$beta.theta.0.0),
					as.numeric(Mi$alpha.gamma),
					as.numeric(Mi$beta.gamma),
					as.numeric(Mi$alpha.theta),
					as.numeric(Mi$beta.theta),
					as.numeric(aperm(Mi$mu.gamma.0)),
					as.numeric(aperm(Mi$tau2.gamma.0)),
					as.numeric(aperm(Mi$mu.theta.0)),
					as.numeric(aperm(Mi$tau2.theta.0)),
					as.numeric(aperm(Mi$mu.gamma)),
					as.numeric(aperm(Mi$mu.theta)),
					as.numeric(aperm(Mi$sigma2.gamma)),
					as.numeric(aperm(Mi$sigma2.theta)))

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

	model_fit = list(id = Mi$Id, sim_type = Mi$sim_type, chains = nchains, nIntervals = Mi$numIntervals,
			Intervals = Mi$Intervals, nBodySys = Mi$numB, maxBs = Mi$maxBs,
			maxAEs = Mi$maxAEs, nAE = Mi$nAE, AE=Mi$AE, B = Mi$B,
			burnin = burnin, iter = iter,
			monitor = monitor,
			gamma = gamma_samples,
			theta = theta_samples,
			mu.gamma = mu.gamma_samples,
			mu.theta = mu.theta_samples,
			sigma2.gamma = sigma2.gamma_samples,
			sigma2.theta = sigma2.theta_samples,
			mu.gamma.0 = mu.gamma.0_samples,
			mu.theta.0 = mu.theta.0_samples,
			tau2.gamma.0 = tau2.gamma.0_samples,
			tau2.theta.0 = tau2.theta.0_samples,
			gamma_acc = gamma_acc,
			theta_acc = theta_acc)
			
	# Model is poisson with BB1a hierarchy and independent intervals
	attr(model_fit, "model") = "1a_pois_indep"

	return(model_fit)
}

Mi$initVars = function() {

    # Data Structure
    Mi$B <- c()
    Mi$numB <- NA
    Mi$numIntervals <- NA
    Mi$nAE <- c()
    Mi$maxAEs <- NA

    # Trial Event Data
    Mi$x <- array()
    Mi$C <- array()
    Mi$y <- array()
    Mi$T <- array()

    # Hyperparameters
    Mi$mu.gamma.0.0 <- NA
    Mi$tau2.gamma.0.0 <- NA
    Mi$mu.theta.0.0 <- NA
    Mi$tau2.theta.0.0 <- NA
    Mi$alpha.gamma.0.0 <- NA
    Mi$beta.gamma.0.0 <- NA
    Mi$alpha.theta.0.0 <- NA
    Mi$beta.theta.0.0 <- NA
    Mi$alpha.gamma <- NA
    Mi$beta.gamma <- NA
    Mi$alpha.theta <- NA
    Mi$beta.theta <- NA

    # Parameters/Simulated values
    # Stage 3
    Mi$mu.gamma.0 <- c()
    Mi$tau2.gamma.0 <- c()
    Mi$mu.theta.0 <- c()
    Mi$tau2.theta.0 <- c()

    # Stage 2
    Mi$mu.gamma <- array()
    Mi$mu.theta <- array()
    Mi$sigma2.gamma <- array()
    Mi$sigma2.theta <- array()

    # Stage 1
    Mi$theta <- array()
    Mi$gamma <- array()
}

Mi$initChains = function(c) {
	# Choose random values for gamma and theta
	for (i in 1:Mi$numIntervals) {
		numB = Mi$numB[i]
		for (b in 1:numB) {
			Mi$gamma[c, i, b, 1:Mi$nAE[i, b]] <- runif(Mi$nAE[i, b], -10, 10)
			Mi$theta[c, i, b, 1:Mi$nAE[i, b]] <- runif(Mi$nAE[i, b], -10, 10)

			Mi$theta[c, i, b, ][is.infinite(Mi$theta[c, i, b, ])] = -10
			Mi$gamma[c, i, b, ][is.infinite(Mi$gamma[c, i, b, ])] = -10

			Mi$theta[c, i, b, ][is.nan(Mi$theta[c, i, b, ])] = -10 # -1000
			Mi$gamma[c, i, b, ][is.nan(Mi$gamma[c, i, b, ])] = -10 # -1000
		}

		Mi$mu.gamma[c, i, 1:numB] = runif(numB, -10, 10)
		Mi$mu.theta[c, i, 1:numB] = runif(numB, -10, 10)
		Mi$sigma2.gamma[c, i, 1:numB] = runif(numB, 5, 20)
		Mi$sigma2.theta[c, i, 1:numB] = runif(numB, 5, 20)

		Mi$mu.gamma.0[c, i] = runif(1, -10, 10)
		Mi$tau2.gamma.0[c, i] = runif(1, 5, 20)
		Mi$mu.theta.0[c, i] = runif(1, -10, 10)
		Mi$tau2.theta.0[c, i] = runif(1, 5, 20)
	}
}

Mi$initialiseChains = function(initial_values, nchains) {

	Mi$theta = array(0, dim=c(nchains, Mi$numIntervals, Mi$maxBs, Mi$maxAEs))
	Mi$gamma = array(0, dim=c(nchains, Mi$numIntervals, Mi$maxBs, Mi$maxAEs))

	if (is.null(initial_values)) {

		# Initialise the first chain with the data
		for (i in 1:Mi$numIntervals) {
			numB = Mi$numB[i]
			for (b in 1:numB) {
				Mi$gamma[1, i, b, ] <- log(Mi$x[i, b,]/Mi$C[i, b, ])
				Mi$theta[1, i, b, ] <- log(Mi$y[i, b,]/Mi$T[i, b, ]) - Mi$gamma[1, i, b, ]

				Mi$theta[1, i, b, ][is.infinite(Mi$theta[1, i, b, ])] = -10 # -1000
				Mi$gamma[1, i, b, ][is.infinite(Mi$gamma[1, i, b, ])] = -10 # -1000

				Mi$theta[1, i, b, ][is.nan(Mi$theta[1, i, b, ])] = -10 # -1000
				Mi$gamma[1, i, b, ][is.nan(Mi$gamma[1, i, b, ])] = -10 # -1000
			}
		}

		Mi$mu.gamma <- array(0, dim = c(nchains, Mi$numIntervals, Mi$maxBs))
		Mi$mu.theta <- array(0, dim = c(nchains, Mi$numIntervals, Mi$maxBs))
		Mi$sigma2.gamma <- array(10, dim = c(nchains, Mi$numIntervals, Mi$maxBs))
		Mi$sigma2.theta <- array(10, dim = c(nchains, Mi$numIntervals, Mi$maxBs))

		Mi$mu.gamma.0 <- array(0, dim = c(nchains, Mi$numIntervals))
		Mi$tau2.gamma.0 <- array(10, dim = c(nchains, Mi$numIntervals))
		Mi$mu.theta.0 <- array(0, dim = c(nchains, Mi$numIntervals))
		Mi$tau2.theta.0 <- array(10, dim = c(nchains, Mi$numIntervals))

		if (nchains > 1) {
			for (c in 2:nchains) {
				Mi$initChains(c)
			}

		}
	}
	else {

		Mi$mu.gamma.0 <- array(0, dim = c(nchains, Mi$numIntervals))
		Mi$tau2.gamma.0 <- array(10, dim = c(nchains, Mi$numIntervals))
		Mi$mu.theta.0 <- array(0, dim = c(nchains, Mi$numIntervals))
		Mi$tau2.theta.0 <- array(10, dim = c(nchains, Mi$numIntervals))

		for (c in 1:nchains) {
			for (i in 1:Mi$numIntervals) {
				interval = Mi$Intervals[i]
				data = initial_values$mu.gamma.0[initial_values$mu.gamma.0$chain == c &
												initial_values$mu.gamma.0$Interval == interval, ]
				Mi$mu.gamma.0[c, i] = data$value
				data = initial_values$mu.theta.0[initial_values$mu.theta.0$chain == c &
												initial_values$mu.theta.0$Interval == interval, ]
				Mi$mu.theta.0[c, i] = data$value

				data = initial_values$tau2.gamma.0[initial_values$tau2.gamma.0$chain == c &
												initial_values$tau2.gamma.0$Interval == interval, ]
				Mi$tau2.gamma.0[c, i] = data$value
				data = initial_values$tau2.theta.0[initial_values$tau2.theta.0$chain == c &
												initial_values$tau2.theta.0$Interval == interval, ]
				Mi$tau2.theta.0[c, i] = data$value

			}
		}

		Mi$mu.gamma <- array(0, dim = c(nchains, Mi$numIntervals, Mi$maxBs))
		Mi$mu.theta <- array(0, dim = c(nchains, Mi$numIntervals, Mi$maxBs))
		Mi$sigma2.gamma <- array(0, dim = c(nchains, Mi$numIntervals, Mi$maxBs))
		Mi$sigma2.theta <- array(0, dim = c(nchains, Mi$numIntervals, Mi$maxBs))

		for (c in 1:nchains) {
			for (i in 1:Mi$numIntervals) {
				interval = Mi$Intervals[i]
				for (b in 1:Mi$numB[i]) {
					data = initial_values$mu.gamma[initial_values$mu.gamma$chain == c &
									initial_values$mu.gamma$Interval == interval
												& initial_values$mu.gamma$B == Mi$B[i, b],]
					Mi$mu.gamma[c, i, b] = data$value

					data = initial_values$mu.theta[initial_values$mu.theta$chain == c &
									initial_values$mu.theta$Interval == interval
												& initial_values$mu.theta$B == Mi$B[i, b],]
					Mi$mu.theta[c, i, b] = data$value

					data = initial_values$sigma2.gamma[initial_values$sigma2.gamma$chain == c &
									initial_values$sigma2.gamma$Interval == interval
												& initial_values$sigma2.gamma$B == Mi$B[i, b],]
					Mi$sigma2.gamma[c, i, b] = data$value

					data = initial_values$sigma2.theta[initial_values$sigma2.theta$chain == c &
									initial_values$sigma2.theta$Interval == interval
												& initial_values$sigma2.theta$B == Mi$B[i, b],]
					Mi$sigma2.theta[c, i, b] = data$value
				}
			}
		}

		for (c in 1:nchains) {
			for (i in 1:Mi$numIntervals) {
				interval = Mi$Intervals[i]
				for (b in 1:Mi$numB[i]) {
					for (j in 1:Mi$nAE[i, b]) {
						ae = Mi$AE[i, b, j]
						data = initial_values$gamma[initial_values$gamma$chain == c
										& initial_values$gamma$Interval == interval
										& initial_values$gamma$B == Mi$B[i, b]
										& initial_values$gamma$AE == ae,]
						Mi$gamma[c, i, b, j] = data$value

						data = initial_values$theta[initial_values$theta$chain == c
										& initial_values$theta$Interval == interval
										& initial_values$theta$B == Mi$B[i, b]
										& initial_values$theta$AE == ae,]
						Mi$theta[c, i, b, j] = data$value
					}
				}
			}
		}
	}
}
