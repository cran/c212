# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 28/04/2015


Md <- new.env()

Md$Id <- "$Id: c212.interim.1a.hier3.lev2.R,v 1.9 2018/10/03 15:40:55 clb13102 Exp clb13102 $"

c212.interim.1a.dep.lev2 <- function(trial.data, sim_type = "SLICE", burnin = 10000, iter = 40000, nchains = 3,
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
	interim = M_global$INTERIMdata(Md, trial.data, iter, nchains, burnin, initial_values)

	if (is.null(interim)) {
		return(NULL)
	}

	trial.data = interim$trial.data
	cntrl.data = interim$cntrl.data

	Md$sim_type <- sim_type

	if (nrow(global.sim.params[global.sim.params$type == sim_type,]) != 1) {
		print("Missing simulation parametetrs");
		return(NULL)
	}

	Md$global.sim.param <- global.sim.params[global.sim.params$type == sim_type,]$value
	Md$global.sim.param_ctrl <- global.sim.params[global.sim.params$type == sim_type,]$control

	if (Md$global.sim.param <= 0) {
		print("Invalid simulation parametetr value");
		return(NULL)
	}

	Md$level = 2

	sim.params = M_global$INTERIM_sim_params1a(Md, sim.params, sim_type, trial.data, cntrl.data)

	monitor = M_global$INTERIM_monitor_1a_3(monitor)

	# Initialise the hyper-parameters
	Md$mu.gamma.0.0 <- hyper_params$mu.gamma.0.0
	Md$tau2.gamma.0.0 <- hyper_params$tau2.gamma.0.0
	Md$alpha.gamma <- hyper_params$alpha.gamma
	Md$beta.gamma <- hyper_params$beta.gamma
	Md$alpha.gamma.0.0 <- hyper_params$alpha.gamma.0.0
	Md$beta.gamma.0.0 <- hyper_params$beta.gamma.0.0

	Md$mu.theta.0.0 <- hyper_params$mu.theta.0.0
	Md$tau2.theta.0.0 <- hyper_params$tau2.theta.0.0
	Md$alpha.theta <- hyper_params$alpha.theta
	Md$beta.theta <- hyper_params$beta.theta
	Md$alpha.theta.0.0 <- hyper_params$alpha.theta.0.0
	Md$beta.theta.0.0 <- hyper_params$beta.theta.0.0

	Ret2 = .Call("c2121a_poisson_mc_exec", as.integer(nchains), as.integer(burnin),
					as.integer(iter), Md$sim_type,
					memory_model,
					as.numeric(Md$global.sim.param),
					as.numeric(Md$global.sim.param_ctrl),
					sim.params,
					monitor,
					as.integer(Md$numIntervals), as.integer(Md$level),
					Md$maxBs, as.integer(Md$numB), as.integer(Md$maxAEs),
					as.integer(t(Md$nAE)), as.integer(aperm(Md$x)), as.integer(aperm(Md$y)),
					as.numeric(aperm(Md$C)),
					as.numeric(aperm(Md$T)),
					as.numeric(aperm(Md$theta)),
					as.numeric(aperm(Md$gamma)),
					as.numeric(Md$mu.gamma.0.0),
					as.numeric(Md$tau2.gamma.0.0),
					as.numeric(Md$mu.theta.0.0),
					as.numeric(Md$tau2.theta.0.0),
					as.numeric(Md$alpha.gamma.0.0),
					as.numeric(Md$beta.gamma.0.0),
					as.numeric(Md$alpha.theta.0.0),
					as.numeric(Md$beta.theta.0.0),
					as.numeric(Md$alpha.gamma),
					as.numeric(Md$beta.gamma),
					as.numeric(Md$alpha.theta),
					as.numeric(Md$beta.theta),
					as.numeric(Md$mu.gamma.0),
					as.numeric(Md$tau2.gamma.0),
					as.numeric(Md$mu.theta.0),
					as.numeric(Md$tau2.theta.0),
					as.numeric(aperm(Md$mu.gamma)),
					as.numeric(aperm(Md$mu.theta)),
					as.numeric(aperm(Md$sigma2.gamma)),
					as.numeric(aperm(Md$sigma2.theta)))

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

	model_fit = list(id = Md$Id, sim_type = Md$sim_type, chains = nchains, nIntervals = Md$numIntervals,
			Intervals = Md$Intervals, nBodySys = Md$numB, maxBs = Md$maxBs,
			maxAEs = Md$maxAEs, nAE = Md$nAE, AE=Md$AE, B = Md$B,
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
	attr(model_fit, "model") = "1a_pois_dep_lev2"

	return(model_fit)
}

Md$initVars = function() {

    # Data Structure
    Md$B <- c()
    Md$numB <- NA
    Md$numIntervals <- NA
    Md$nAE <- c()
    Md$maxAEs <- NA

    # Trial Event Data
    Md$x <- array()
    Md$C <- array()
    Md$y <- array()
    Md$T <- array()

    # Hyperparameters
    Md$mu.gamma.0.0 <- NA
    Md$tau2.gamma.0.0 <- NA
    Md$mu.theta.0.0 <- NA
    Md$tau2.theta.0.0 <- NA
    Md$alpha.gamma.0.0 <- NA
    Md$beta.gamma.0.0 <- NA
    Md$alpha.theta.0.0 <- NA
    Md$beta.theta.0.0 <- NA
    Md$alpha.gamma <- NA
    Md$beta.gamma <- NA
    Md$alpha.theta <- NA
    Md$beta.theta <- NA

	# Parameters/Simulated values
    # Stage 3
    Md$mu.gamma.0 <- c()
    Md$tau2.gamma.0 <- c()
    Md$mu.theta.0 <- c()
    Md$tau2.theta.0 <- c()

    # Stage 2
    Md$mu.gamma <- array()
    Md$mu.theta <- array()
    Md$sigma2.gamma <- array()
    Md$sigma2.theta <- array()

    # Stage 1
    Md$theta <- array()
    Md$gamma <- array()
}

Md$initChains = function(c) {
	# Choose random values for gamma and theta
	for (i in 1:Md$numIntervals) {
		numB = Md$numB[i]
		for (b in 1:numB) {
			Md$gamma[c, i, b, 1:Md$nAE[i, b]] <- runif(Md$nAE[i, b], -10, 10)
			Md$theta[c, i, b, 1:Md$nAE[i, b]] <- runif(Md$nAE[i, b], -10, 10)

			Md$theta[c, i, b, ][is.infinite(Md$theta[c, i, b, ])] = -10
			Md$gamma[c, i, b, ][is.infinite(Md$gamma[c, i, b, ])] = -10

			Md$theta[c, i, b, ][is.nan(Md$theta[c, i, b, ])] = -10 # -1000
			Md$gamma[c, i, b, ][is.nan(Md$gamma[c, i, b, ])] = -10 # -1000
		}

		Md$mu.gamma[c, i, 1:numB] = runif(numB, -10, 10)
		Md$mu.theta[c, i, 1:numB] = runif(numB, -10, 10)
		Md$sigma2.gamma[c, i, 1:numB] = runif(numB, 5, 20)
		Md$sigma2.theta[c, i, 1:numB] = runif(numB, 5, 20)

		Md$mu.gamma.0[c] = runif(1, -10, 10)
		Md$tau2.gamma.0[c] = runif(1, 5, 20)
		Md$mu.theta.0[c] = runif(1, -10, 10)
		Md$tau2.theta.0[c] = runif(1, 5, 20)
	}
}

Md$initialiseChains = function(initial_values, nchains) {

	Md$theta = array(0, dim=c(nchains, Md$numIntervals, Md$maxBs, Md$maxAEs))
	Md$gamma = array(0, dim=c(nchains, Md$numIntervals, Md$maxBs, Md$maxAEs))

	if (is.null(initial_values)) {

		# Initialise the first chain with the data
		for (i in 1:Md$numIntervals) {
			numB = Md$numB[i]
			for (b in 1:numB) {
				Md$gamma[1, i, b, ] <- log(Md$x[i, b,]/Md$C[i, b, ])
				Md$theta[1, i, b, ] <- log(Md$y[i, b,]/Md$T[i, b, ]) - Md$gamma[1, i, b, ]

				Md$theta[1, i, b, ][is.infinite(Md$theta[1, i, b, ])] = -10 # -1000
				Md$gamma[1, i, b, ][is.infinite(Md$gamma[1, i, b, ])] = -10 # -1000

				Md$theta[1, i, b, ][is.nan(Md$theta[1, i, b, ])] = -10 # -1000
				Md$gamma[1, i, b, ][is.nan(Md$gamma[1, i, b, ])] = -10 # -1000
			}
		}

		Md$mu.gamma <- array(0, dim = c(nchains, Md$numIntervals, Md$maxBs))
		Md$mu.theta <- array(0, dim = c(nchains, Md$numIntervals, Md$maxBs))
		Md$sigma2.gamma <- array(10, dim = c(nchains, Md$numIntervals, Md$maxBs))
		Md$sigma2.theta <- array(10, dim = c(nchains, Md$numIntervals, Md$maxBs))

		Md$mu.gamma.0 <- rep(0, nchains)
		Md$tau2.gamma.0 <- rep(10, nchains)
		Md$mu.theta.0 <- rep(0, nchains)
		Md$tau2.theta.0 <- rep(10, nchains)

		if (nchains > 1) {
			for (c in 2:nchains) {
				Md$initChains(c)
			}

		}
	}
	else {

		Md$mu.gamma.0 <- rep(0, nchains)
		Md$tau2.gamma.0 <- rep(10, nchains)
		Md$mu.theta.0 <- rep(0, nchains)
		Md$tau2.theta.0 <- rep(10, nchains)

		for (c in 1:nchains) {
			Md$mu.gamma.0[c] = initial_values$mu.gamma.0[c]
			Md$mu.theta.0[c] = initial_values$mu.theta.0[c]
			Md$tau2.gamma.0[c] = initial_values$tau2.gamma.0[c]
			Md$tau2.theta.0[c] = initial_values$tau2.theta.0[c]
		}

		Md$mu.gamma <- array(0, dim = c(nchains, Md$numIntervals, Md$maxBs))
		Md$mu.theta <- array(0, dim = c(nchains, Md$numIntervals, Md$maxBs))
		Md$sigma2.gamma <- array(0, dim = c(nchains, Md$numIntervals, Md$maxBs))
		Md$sigma2.theta <- array(0, dim = c(nchains, Md$numIntervals, Md$maxBs))

		for (c in 1:nchains) {
			for (i in 1:Md$numIntervals) {
				interval = Md$Intervals[i]
				for (b in 1:Md$numB[i]) {
					data = initial_values$mu.gamma[initial_values$mu.gamma$chain == c &
									initial_values$mu.gamma$Interval == interval
												& initial_values$mu.gamma$B == Md$B[i, b],]
					Md$mu.gamma[c, i, b] = data$value

					data = initial_values$mu.theta[initial_values$mu.theta$chain == c &
									initial_values$mu.theta$Interval == interval
												& initial_values$mu.theta$B == Md$B[i, b],]
					Md$mu.theta[c, i, b] = data$value

					data = initial_values$sigma2.gamma[initial_values$sigma2.gamma$chain == c &
									initial_values$sigma2.gamma$Interval == interval
												& initial_values$sigma2.gamma$B == Md$B[i, b],]
					Md$sigma2.gamma[c, i, b] = data$value

					data = initial_values$sigma2.theta[initial_values$sigma2.theta$chain == c &
									initial_values$sigma2.theta$Interval == interval
												& initial_values$sigma2.theta$B == Md$B[i, b],]
					Md$sigma2.theta[c, i, b] = data$value
				}
			}
		}

		for (c in 1:nchains) {
			for (i in 1:Md$numIntervals) {
				interval = Md$Intervals[i]
				for (b in 1:Md$numB[i]) {
					for (j in 1:Md$nAE[i, b]) {
						ae = Md$AE[i, b, j]
						data = initial_values$gamma[initial_values$gamma$chain == c
										& initial_values$gamma$Interval == interval
										& initial_values$gamma$B == Md$B[i, b]
										& initial_values$gamma$AE == ae,]
						Md$gamma[c, i, b, j] = data$value

						data = initial_values$theta[initial_values$theta$chain == c
										& initial_values$theta$Interval == interval
										& initial_values$theta$B == Md$B[i, b]
										& initial_values$theta$AE == ae,]
						Md$theta[c, i, b, j] = data$value
					}
				}
			}
		}
	}
}
