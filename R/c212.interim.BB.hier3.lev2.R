# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 10/06/2015


Md1 <- new.env()

Md1$Id <- "$Id: c212.interim.BB.hier3.lev2.R,v 1.17 2018/10/03 15:40:56 clb13102 Exp clb13102 $"

c212.interim.BB.dep.lev2 <- function(trial.data, sim_type = "SLICE", burnin = 10000, iter = 60000, nchains = 5,
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
	interim = M_global$INTERIMdata(Md1, trial.data, iter, nchains, burnin, initial_values)

	if (is.null(interim)) {
		return(NULL)
	}

	trial.data = interim$trial.data
	cntrl.data = interim$cntrl.data

	Md1$Algo <- theta_algorithm

	Md1$sim_type <- sim_type

	if (nrow(global.sim.params[global.sim.params$type == sim_type,]) == 0) {
		print("Missing simulation parametetrs");
		return(NULL)
	}

	if (!all(global.sim.params$value > 0)) {
		print("Invalid simulation parameter value");
		return(NULL)
	}

	Md1$global.sim.params <- global.sim.params

	Md1$Level = 2

	sp = M_global$INTERIM_sim_paramsBB_3(Md1, sim.params, pm.weights, sim_type, trial.data, cntrl.data)

	sim.params = sp$sim.params
	pm.weights = sp$pm.weights

	monitor = M_global$INTERIM_monitor_BB_3(monitor)

	# Initialise the hyper-parameters
	Md1$mu.gamma.0.0 <- hyper_params$mu.gamma.0.0
	Md1$tau2.gamma.0.0 <- hyper_params$tau2.gamma.0.0
	Md1$alpha.gamma <- hyper_params$alpha.gamma
	Md1$beta.gamma <- hyper_params$beta.gamma
	Md1$alpha.gamma.0.0 <- hyper_params$alpha.gamma.0.0
	Md1$beta.gamma.0.0 <- hyper_params$beta.gamma.0.0

	Md1$mu.theta.0.0 <- hyper_params$mu.theta.0.0
	Md1$tau2.theta.0.0 <- hyper_params$tau2.theta.0.0
	Md1$alpha.theta <- hyper_params$alpha.theta
	Md1$beta.theta <- hyper_params$beta.theta
	Md1$alpha.theta.0.0 <- hyper_params$alpha.theta.0.0
	Md1$beta.theta.0.0 <- hyper_params$beta.theta.0.0

	# BB2004 parameters
	Md1$lambda.alpha <- hyper_params$lambda.alpha
	Md1$lambda.beta <- hyper_params$lambda.beta

	algo = 1
	if (Md1$Algo == "BB2004") {
		algo <- 1;
	} else {
		if (Md1$Algo == "MH") {
			algo <- 2;
		} else {
			if (Md1$Algo == "Adapt") {
				algo <- 3;
			} else {
				if (Md1$Algo == "Indep") {
					algo <- 4;
				} else {
					algo <- 1;
				}
			}
		}
	}


	Ret2 = .Call("c212BB_poisson_mc_exec", as.integer(nchains), as.integer(burnin),
					as.integer(iter), Md1$sim_type,
					memory_model, Md1$global.sim.params,
					sim.params,
					as.numeric(global.pm.weight),
					pm.weights,
					monitor,
					as.integer(Md1$numIntervals), as.integer(Md1$Level),
					Md1$maxBs, as.integer(Md1$numB), as.integer(Md1$maxAEs),
					as.integer(t(Md1$nAE)), as.integer(aperm(Md1$x)), as.integer(aperm(Md1$y)),
					as.numeric(aperm(Md1$C)),
					as.numeric(aperm(Md1$T)),
					as.numeric(aperm(Md1$theta)),
					as.numeric(aperm(Md1$gamma)),
					as.numeric(Md1$mu.gamma.0.0),
					as.numeric(Md1$tau2.gamma.0.0),
					as.numeric(Md1$mu.theta.0.0),
					as.numeric(Md1$tau2.theta.0.0),
					as.numeric(Md1$alpha.gamma.0.0),
					as.numeric(Md1$beta.gamma.0.0),
					as.numeric(Md1$alpha.theta.0.0),
					as.numeric(Md1$beta.theta.0.0),
					as.numeric(Md1$alpha.gamma),
					as.numeric(Md1$beta.gamma),
					as.numeric(Md1$alpha.theta),
					as.numeric(Md1$beta.theta),
					as.numeric(Md1$mu.gamma.0),
					as.numeric(Md1$tau2.gamma.0),
					as.numeric(Md1$mu.theta.0),
					as.numeric(Md1$tau2.theta.0),
					as.numeric(aperm(Md1$mu.gamma)),
					as.numeric(aperm(Md1$mu.theta)),
					as.numeric(aperm(Md1$sigma2.gamma)),
					as.numeric(aperm(Md1$sigma2.theta)),
					as.numeric(aperm(Md1$pi)),
					as.numeric(Md1$alpha.pi),
					as.numeric(Md1$beta.pi),
					as.numeric(Md1$lambda.alpha),
					as.numeric(Md1$lambda.beta),
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

	model_fit = list(id = Md1$Id, sim_type = Md1$sim_type, chains = nchains, nIntervals = Md1$numIntervals,
			Intervals = Md1$Intervals, nBodySys = Md1$numB, maxBs = Md1$maxBs,
			maxAEs = Md1$maxAEs, nAE = Md1$nAE, AE=Md1$AE, B = Md1$B,
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
			
	# Model is poisson with BB hierarchy and dependent intervals
	attr(model_fit, "model") = "BB_pois_dep_lev2"

	return(model_fit)
}

Md1$initVars = function() {

    # Data Structure
    Md1$B <- c()
    Md1$numB <- NA
    Md1$numIntervals <- NA
    Md1$nAE <- c()
    Md1$maxAEs <- NA

    # Trial Event Data
    Md1$x <- array()
    Md1$C <- array()
    Md1$y <- array()
    Md1$T <- array()

    # Hyperparameters
    Md1$mu.gamma.0.0 <- NA
    Md1$tau2.gamma.0.0 <- NA
    Md1$mu.theta.0.0 <- NA
    Md1$tau2.theta.0.0 <- NA
    Md1$alpha.gamma.0.0 <- NA
    Md1$beta.gamma.0.0 <- NA
    Md1$alpha.theta.0.0 <- NA
    Md1$beta.theta.0.0 <- NA
    Md1$alpha.gamma <- NA
    Md1$beta.gamma <- NA
    Md1$alpha.theta <- NA
    Md1$beta.theta <- NA

	# Parameters/Simulated values
    # Stage 3
    Md1$mu.gamma.0 <- c()
    Md1$tau2.gamma.0 <- c()
    Md1$mu.theta.0 <- c()
    Md1$tau2.theta.0 <- c()

    # Stage 2
    Md1$mu.gamma <- array()
    Md1$mu.theta <- array()
    Md1$sigma2.gamma <- array()
    Md1$sigma2.theta <- array()

    # Stage 1
    Md1$theta <- array()
    Md1$gamma <- array()

	# BB2004 parameters
	Md1$lambda.alpha <- NA
	Md1$lambda.beta <- NA
	Md1$alpha.pi <- NA
	Md1$beta.pi <- NA
	Md1$pi <- NA
}

Md1$initChains = function(c) {
	# Choose random values for gamma and theta
	for (i in 1:Md1$numIntervals) {
		numB = Md1$numB[i]
		for (b in 1:numB) {
			Md1$gamma[c, i, b, 1:Md1$nAE[i, b]] <- runif(Md1$nAE[i, b], -10, 10)
			Md1$theta[c, i, b, 1:Md1$nAE[i, b]] <- runif(Md1$nAE[i, b], -10, 10)

			Md1$theta[c, i, b, ][is.infinite(Md1$theta[c, i, b, ])] = -10
			Md1$gamma[c, i, b, ][is.infinite(Md1$gamma[c, i, b, ])] = -10

			Md1$theta[c, i, b, ][is.nan(Md1$theta[c, i, b, ])] = -10 # -1000
			Md1$gamma[c, i, b, ][is.nan(Md1$gamma[c, i, b, ])] = -10 # -1000
		}

		Md1$mu.gamma[c, i, 1:numB] = runif(numB, -10, 10)
		Md1$mu.theta[c, i, 1:numB] = runif(numB, -10, 10)
		Md1$sigma2.gamma[c, i, 1:numB] = runif(numB, 5, 20)
		Md1$sigma2.theta[c, i, 1:numB] = runif(numB, 5, 20)

		Md1$pi[c, i, 1:numB] = runif(numB, 0, 1)

		Md1$mu.gamma.0[c] = runif(1, -10, 10)
		Md1$tau2.gamma.0[c] = runif(1, 5, 20)
		Md1$mu.theta.0[c] = runif(1, -10, 10)
		Md1$tau2.theta.0[c] = runif(1, 5, 20)

		Md1$alpha.pi[c] = runif(1, 1.25, 100)
		Md1$beta.pi[c] = runif(1, 1.25, 100)
	}
}

Md1$initialiseChains = function(initial_values, nchains) {

	Md1$theta = array(0, dim=c(nchains, Md1$numIntervals, Md1$maxBs, Md1$maxAEs))
	Md1$gamma = array(0, dim=c(nchains, Md1$numIntervals, Md1$maxBs, Md1$maxAEs))

	if (is.null(initial_values)) {

		# Initialise the first chain with the data
		for (i in 1:Md1$numIntervals) {
			numB = Md1$numB[i]
			for (b in 1:numB) {
				Md1$gamma[1, i, b, ] <- log(Md1$x[i, b,]/Md1$C[i, b, ])
				Md1$theta[1, i, b, ] <- log(Md1$y[i, b,]/Md1$T[i, b, ]) - Md1$gamma[1, i, b, ]

				Md1$theta[1, i, b, ][is.infinite(Md1$theta[1, i, b, ])] = -10 # -1000
				Md1$gamma[1, i, b, ][is.infinite(Md1$gamma[1, i, b, ])] = -10 # -1000

				Md1$theta[1, i, b, ][is.nan(Md1$theta[1, i, b, ])] = -10 # -1000
				Md1$gamma[1, i, b, ][is.nan(Md1$gamma[1, i, b, ])] = -10 # -1000
			}
		}

		Md1$mu.gamma <- array(0, dim = c(nchains, Md1$numIntervals, Md1$maxBs))
		Md1$mu.theta <- array(0, dim = c(nchains, Md1$numIntervals, Md1$maxBs))
		Md1$sigma2.gamma <- array(10, dim = c(nchains, Md1$numIntervals, Md1$maxBs))
		Md1$sigma2.theta <- array(10, dim = c(nchains, Md1$numIntervals, Md1$maxBs))

		Md1$pi <- array(0.5, dim = c(nchains, Md1$numIntervals, Md1$maxBs))

		Md1$mu.gamma.0 <- rep(0, nchains)
		Md1$tau2.gamma.0 <- rep(10, nchains)
		Md1$mu.theta.0 <- rep(0, nchains)
		Md1$tau2.theta.0 <- rep(10, nchains)

		Md1$alpha.pi <- rep(1.5, nchains)
		Md1$beta.pi <- rep(1.5, nchains)

		if (nchains > 1) {
			for (c in 2:nchains) {
				Md1$initChains(c)
			}
		}
	}
	else {

		Md1$mu.gamma.0 <- rep(0, nchains)
		Md1$tau2.gamma.0 <- rep(10, nchains)
		Md1$mu.theta.0 <- rep(0, nchains)
		Md1$tau2.theta.0 <- rep(10, nchains)

		Md1$alpha.pi <- rep(10, nchains)
		Md1$beta.pi <- rep(10, nchains)

		for (c in 1:nchains) {
			Md1$mu.gamma.0[c] = initial_values$mu.gamma.0[c]
			Md1$tau2.gamma.0[c] = initial_values$tau2.gamma.0[c]
			Md1$mu.theta.0[c] = initial_values$mu.theta.0[c]
			Md1$tau2.theta.0[c] = initial_values$tau2.theta.0[c]

			Md1$alpha.pi[c] = initial_values$alpha.pi[c]
			Md1$beta.pi[c] = initial_values$beta.pi[c]
		}

		Md1$mu.gamma <- array(0, dim = c(nchains, Md1$numIntervals, Md1$maxBs))
		Md1$mu.theta <- array(0, dim = c(nchains, Md1$numIntervals, Md1$maxBs))
		Md1$sigma2.gamma <- array(0, dim = c(nchains, Md1$numIntervals, Md1$maxBs))
		Md1$sigma2.theta <- array(0, dim = c(nchains, Md1$numIntervals, Md1$maxBs))

		Md1$pi <- array(0.5, dim = c(nchains, Md1$numIntervals, Md1$maxBs))

		for (c in 1:nchains) {
			for (i in 1:Md1$numIntervals) {
				interval = Md1$Intervals[i]
				for (b in 1:Md1$numB[i]) {
					data = initial_values$mu.gamma[initial_values$mu.gamma$chain == c &
									initial_values$mu.gamma$Interval == interval
												& initial_values$mu.gamma$B == Md1$B[i, b],]
					Md1$mu.gamma[c, i, b] = data$value

					data = initial_values$mu.theta[initial_values$mu.theta$chain == c &
									initial_values$mu.theta$Interval == interval
												& initial_values$mu.theta$B == Md1$B[i, b],]
					Md1$mu.theta[c, i, b] = data$value

					data = initial_values$sigma2.gamma[initial_values$sigma2.gamma$chain == c &
									initial_values$sigma2.gamma$Interval == interval
												& initial_values$sigma2.gamma$B == Md1$B[i, b],]
					Md1$sigma2.gamma[c, i, b] = data$value

					data = initial_values$sigma2.theta[initial_values$sigma2.theta$chain == c &
									initial_values$sigma2.theta$Interval == interval
												& initial_values$sigma2.theta$B == Md1$B[i, b],]
					Md1$sigma2.theta[c, i, b] = data$value

					data = initial_values$pi[initial_values$pi$chain == c &
									initial_values$pi$Interval == interval
												& initial_values$pi$B == Md1$B[i, b],]
					Md1$pi[c, i, b] = data$value
				}
			}
		}

		for (c in 1:nchains) {
			for (i in 1:Md1$numIntervals) {
				interval = Md1$Intervals[i]
				for (b in 1:Md1$numB[i]) {
					for (j in 1:Md1$nAE[i, b]) {
						ae = Md1$AE[i, b, j]
						data = initial_values$gamma[initial_values$gamma$chain == c
										& initial_values$gamma$Interval == interval
										& initial_values$gamma$B == Md1$B[i, b]
										& initial_values$gamma$AE == ae,]
						Md1$gamma[c, i, b, j] = data$value

						data = initial_values$theta[initial_values$theta$chain == c
										& initial_values$theta$Interval == interval
										& initial_values$theta$B == Md1$B[i, b]
										& initial_values$theta$AE == ae,]
						Md1$theta[c, i, b, j] = data$value
					}
				}
			}
		}
	}
}
