# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 28/04/2015


Mi_1a_h2_l1 <- new.env()

Mi_1a_h2_l1$Id <- "$Id: c212.interim.1a.hier2.lev1.R,v 1.9 2018/10/03 15:40:56 clb13102 Exp clb13102 $"

c212.interim.1a.hier2.lev1 <- function(trial.data, sim_type = "SLICE", burnin = 10000,
	iter = 40000, nchains = 3,
	global.sim.params = data.frame(type = c("MH", "SLICE"), param = c("sigma_MH", "w"),
	value = c(0.2,1),
	control = c(0,6)),
	sim.params = NULL,
	monitor = data.frame(variable = c("theta", "gamma", "mu.gamma", "mu.theta",
		"sigma2.theta", "sigma2.gamma"),
		monitor = c(1, 1, 1, 1, 1, 1), stringsAsFactors = FALSE),
	initial_values = NULL,
	hyper_params = list(mu.gamma.0 = 0, tau2.gamma.0 = 10,
	mu.theta.0 = 0, tau2.theta.0 = 10, alpha.gamma = 3, beta.gamma = 1, alpha.theta = 3,
	beta.theta = 1), memory_model = "HIGH")
{
	interim = M_global$INTERIMdata(Mi_1a_h2_l1, trial.data, iter, nchains, burnin, initial_values)

	if (is.null(interim)) {
		return(NULL)
	}

	trial.data = interim$trial.data
	cntrl.data = interim$cntrl.data

	if (M_global$checkBS(Mi_1a_h2_l1, cntrl.data)) {
		return(NULL)
	}

	Mi_1a_h2_l1$sim_type <- sim_type

	if (nrow(global.sim.params[global.sim.params$type == sim_type,]) != 1) {
		print("Missing simulation parametetrs");
		return(NULL)
	}

	Mi_1a_h2_l1$global.sim.param <- global.sim.params[global.sim.params$type == sim_type,]$value
	Mi_1a_h2_l1$global.sim.param_ctrl <- global.sim.params[global.sim.params$type == sim_type,]$control

	if (Mi_1a_h2_l1$global.sim.param <= 0) {
		print("Invalid simulation parametetr value");
		return(NULL)
	}

	Mi_1a_h2_l1$level = 1

	sim.params = M_global$INTERIM_sim_params1a(Mi_1a_h2_l1, sim.params, sim_type, trial.data, cntrl.data)

	monitor = M_global$INTERIM_monitor_1a_2(monitor)

	# Initialise the hyper-parameters
	Mi_1a_h2_l1$mu.gamma.0 <- hyper_params$mu.gamma.0
	Mi_1a_h2_l1$tau2.gamma.0 <- hyper_params$tau2.gamma.0
	Mi_1a_h2_l1$alpha.gamma <- hyper_params$alpha.gamma
	Mi_1a_h2_l1$beta.gamma <- hyper_params$beta.gamma

	Mi_1a_h2_l1$mu.theta.0 <- hyper_params$mu.theta.0
	Mi_1a_h2_l1$tau2.theta.0 <- hyper_params$tau2.theta.0
	Mi_1a_h2_l1$alpha.theta <- hyper_params$alpha.theta
	Mi_1a_h2_l1$beta.theta <- hyper_params$beta.theta

	Ret2 = .Call("c2121a_interim_hier2_exec", as.integer(nchains), as.integer(burnin),
					as.integer(iter), Mi_1a_h2_l1$sim_type, 
					memory_model,
					as.numeric(Mi_1a_h2_l1$global.sim.param),
					as.numeric(Mi_1a_h2_l1$global.sim.param_ctrl),
					sim.params,
					monitor,
					as.integer(Mi_1a_h2_l1$numIntervals),
					as.integer(Mi_1a_h2_l1$level),
					Mi_1a_h2_l1$maxBs, as.integer(Mi_1a_h2_l1$numB),
					as.integer(Mi_1a_h2_l1$maxAEs),
					as.integer(t(Mi_1a_h2_l1$nAE)),
					as.integer(aperm(Mi_1a_h2_l1$x)),
					as.integer(aperm(Mi_1a_h2_l1$y)),
					as.numeric(aperm(Mi_1a_h2_l1$C)),
					as.numeric(aperm(Mi_1a_h2_l1$T)),
					as.numeric(aperm(Mi_1a_h2_l1$theta)),
					as.numeric(aperm(Mi_1a_h2_l1$gamma)),
					as.numeric(Mi_1a_h2_l1$mu.gamma.0),
					as.numeric(Mi_1a_h2_l1$tau2.gamma.0),
					as.numeric(Mi_1a_h2_l1$mu.theta.0),
					as.numeric(Mi_1a_h2_l1$tau2.theta.0),
					as.numeric(Mi_1a_h2_l1$alpha.gamma),
					as.numeric(Mi_1a_h2_l1$beta.gamma),
					as.numeric(Mi_1a_h2_l1$alpha.theta),
					as.numeric(Mi_1a_h2_l1$beta.theta),
					as.numeric(aperm(Mi_1a_h2_l1$mu.gamma)),
					as.numeric(aperm(Mi_1a_h2_l1$mu.theta)),
					as.numeric(aperm(Mi_1a_h2_l1$sigma2.gamma)),
					as.numeric(aperm(Mi_1a_h2_l1$sigma2.theta)))

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

	model_fit = list(id = Mi_1a_h2_l1$Id,
			sim_type = Mi_1a_h2_l1$sim_type,
			chains = nchains,
			nIntervals = Mi_1a_h2_l1$numIntervals,
			Intervals = Mi_1a_h2_l1$Intervals,
			nBodySys = Mi_1a_h2_l1$numB,
			maxBs = Mi_1a_h2_l1$maxBs,
			maxAEs = Mi_1a_h2_l1$maxAEs,
			nAE = Mi_1a_h2_l1$nAE,
			AE=Mi_1a_h2_l1$AE,
			B = Mi_1a_h2_l1$B,
			burnin = burnin,
			iter = iter,
			monitor = monitor,
			gamma = gamma_samples,
			theta = theta_samples,
			mu.gamma = mu.gamma_samples,
			mu.theta = mu.theta_samples,
			sigma2.gamma = sigma2.gamma_samples,
			sigma2.theta = sigma2.theta_samples,
			gamma_acc = gamma_acc,
			theta_acc = theta_acc)
			
	# Model is poisson with BB1a hierarchy and independent intervals
	attr(model_fit, "model") = "1a_pois_h2_l1"

	return(model_fit)
}

Mi_1a_h2_l1$initVars = function() {

    # Data Structure
    Mi_1a_h2_l1$B <- c()
    Mi_1a_h2_l1$numB <- NA
    Mi_1a_h2_l1$numIntervals <- NA
    Mi_1a_h2_l1$nAE <- c()
    Mi_1a_h2_l1$maxAEs <- NA

    # Trial Event Data
    Mi_1a_h2_l1$x <- array()
    Mi_1a_h2_l1$C <- array()
    Mi_1a_h2_l1$y <- array()
    Mi_1a_h2_l1$T <- array()

    # Hyperparameters
    Mi_1a_h2_l1$mu.gamma.0 <- NA
    Mi_1a_h2_l1$tau2.gamma.0 <- NA
    Mi_1a_h2_l1$mu.theta.0 <- NA
    Mi_1a_h2_l1$tau2.theta.0 <- NA
    Mi_1a_h2_l1$alpha.gamma <- NA
    Mi_1a_h2_l1$beta.gamma <- NA
    Mi_1a_h2_l1$alpha.theta <- NA
    Mi_1a_h2_l1$beta.theta <- NA

    # Parameters/Simulated values
    # Stage 2
    Mi_1a_h2_l1$mu.gamma <- array()
    Mi_1a_h2_l1$mu.theta <- array()
    Mi_1a_h2_l1$sigma2.gamma <- array()
    Mi_1a_h2_l1$sigma2.theta <- array()

    # Stage 1
    Mi_1a_h2_l1$theta <- array()
    Mi_1a_h2_l1$gamma <- array()
}

Mi_1a_h2_l1$initChains = function(c) {
	# Choose random values for gamma and theta
	for (i in 1:Mi_1a_h2_l1$numIntervals) {
		numB = Mi_1a_h2_l1$numB[i]
		for (b in 1:numB) {
			Mi_1a_h2_l1$gamma[c, i, b, 1:Mi_1a_h2_l1$nAE[i, b]] <- runif(Mi_1a_h2_l1$nAE[i, b], -10, 10)
			Mi_1a_h2_l1$theta[c, i, b, 1:Mi_1a_h2_l1$nAE[i, b]] <- runif(Mi_1a_h2_l1$nAE[i, b], -10, 10)

			Mi_1a_h2_l1$theta[c, i, b, ][is.infinite(Mi_1a_h2_l1$theta[c, i, b, ])] = -10
			Mi_1a_h2_l1$gamma[c, i, b, ][is.infinite(Mi_1a_h2_l1$gamma[c, i, b, ])] = -10

			Mi_1a_h2_l1$theta[c, i, b, ][is.nan(Mi_1a_h2_l1$theta[c, i, b, ])] = -10 # -1000
			Mi_1a_h2_l1$gamma[c, i, b, ][is.nan(Mi_1a_h2_l1$gamma[c, i, b, ])] = -10 # -1000
		}
	}

	Mi_1a_h2_l1$mu.theta[c, 1:numB] = runif(numB, -10, 10)
	Mi_1a_h2_l1$mu.gamma[c, 1:numB] = runif(numB, -10, 10)
	Mi_1a_h2_l1$sigma2.gamma[c, 1:numB] = runif(numB, 5, 20)
	Mi_1a_h2_l1$sigma2.theta[c, 1:numB] = runif(numB, 5, 20)
}

Mi_1a_h2_l1$initialiseChains = function(initial_values, nchains) {

	Mi_1a_h2_l1$theta = array(0, dim=c(nchains, Mi_1a_h2_l1$numIntervals, Mi_1a_h2_l1$maxBs, Mi_1a_h2_l1$maxAEs))
	Mi_1a_h2_l1$gamma = array(0, dim=c(nchains, Mi_1a_h2_l1$numIntervals, Mi_1a_h2_l1$maxBs, Mi_1a_h2_l1$maxAEs))

	if (is.null(initial_values)) {

		# Initialise the first chain with the data
		for (i in 1:Mi_1a_h2_l1$numIntervals) {
			numB = Mi_1a_h2_l1$numB[i]
			for (b in 1:numB) {
				Mi_1a_h2_l1$gamma[1, i, b, ] <- log(Mi_1a_h2_l1$x[i, b,]/Mi_1a_h2_l1$C[i, b, ])
				Mi_1a_h2_l1$theta[1, i, b, ] <- log(Mi_1a_h2_l1$y[i, b,]/Mi_1a_h2_l1$T[i, b, ]) - Mi_1a_h2_l1$gamma[1, i, b, ]

				Mi_1a_h2_l1$theta[1, i, b, ][is.infinite(Mi_1a_h2_l1$theta[1, i, b, ])] = -10 # -1000
				Mi_1a_h2_l1$gamma[1, i, b, ][is.infinite(Mi_1a_h2_l1$gamma[1, i, b, ])] = -10 # -1000

				Mi_1a_h2_l1$theta[1, i, b, ][is.nan(Mi_1a_h2_l1$theta[1, i, b, ])] = -10 # -1000
				Mi_1a_h2_l1$gamma[1, i, b, ][is.nan(Mi_1a_h2_l1$gamma[1, i, b, ])] = -10 # -1000
			}
		}

		Mi_1a_h2_l1$mu.theta <- array(0, dim = c(nchains, Mi_1a_h2_l1$maxBs))


		Mi_1a_h2_l1$mu.gamma <- array(0, dim = c(nchains, Mi_1a_h2_l1$maxBs))
		Mi_1a_h2_l1$sigma2.gamma <- array(10, dim = c(nchains, Mi_1a_h2_l1$maxBs))
		Mi_1a_h2_l1$sigma2.theta <- array(10, dim = c(nchains, Mi_1a_h2_l1$maxBs))

		if (nchains > 1) {
			for (c in 2:nchains) {
				Mi_1a_h2_l1$initChains(c)
			}

		}
	}
	else {

		Mi_1a_h2_l1$mu.theta <- array(0, dim = c(nchains, Mi_1a_h2_l1$maxBs))

		Mi_1a_h2_l1$mu.gamma <- array(0, dim = c(nchains, Mi_1a_h2_l1$maxBs))
		Mi_1a_h2_l1$sigma2.gamma <- array(0, dim = c(nchains, Mi_1a_h2_l1$maxBs))
		Mi_1a_h2_l1$sigma2.theta <- array(0, dim = c(nchains, Mi_1a_h2_l1$maxBs))

		for (c in 1:nchains) {
			# Have the same 
			for (b in 1:length(Mi_1a_h2_l1$B[1,])) {
				data = initial_values$mu.theta[initial_values$mu.theta$chain == c
												& initial_values$mu.theta$B == Mi_1a_h2_l1$B[1, b],]
				Mi_1a_h2_l1$mu.theta[c, b] = data$value

				data = initial_values$mu.gamma[initial_values$mu.gamma$chain == c
												& initial_values$mu.gamma$B == Mi_1a_h2_l1$B[1, b],]
				Mi_1a_h2_l1$mu.gamma[c, b] = data$value

				data = initial_values$sigma2.theta[initial_values$sigma2.theta$chain == c
												& initial_values$sigma2.theta$B == Mi_1a_h2_l1$B[1, b],]
				Mi_1a_h2_l1$sigma2.theta[c, b] = data$value

				data = initial_values$sigma2.gamma[initial_values$sigma2.gamma$chain == c
												& initial_values$sigma2.gamma$B == Mi_1a_h2_l1$B[1, b],]
				Mi_1a_h2_l1$sigma2.gamma[c, b] = data$value
			}
		}

		for (c in 1:nchains) {
			for (i in 1:Mi_1a_h2_l1$numIntervals) {
				interval = Mi_1a_h2_l1$Intervals[i]
				for (b in 1:Mi_1a_h2_l1$numB[i]) {
					for (j in 1:Mi_1a_h2_l1$nAE[i, b]) {
						ae = Mi_1a_h2_l1$AE[i, b, j]
						data = initial_values$gamma[initial_values$gamma$chain == c
										& initial_values$gamma$Interval == interval
										& initial_values$gamma$B == Mi_1a_h2_l1$B[i, b]
										& initial_values$gamma$AE == ae,]
						Mi_1a_h2_l1$gamma[c, i, b, j] = data$value

						data = initial_values$theta[initial_values$theta$chain == c
										& initial_values$theta$Interval == interval
										& initial_values$theta$B == Mi_1a_h2_l1$B[i, b]
										& initial_values$theta$AE == ae,]
						Mi_1a_h2_l1$theta[c, i, b, j] = data$value
					}
				}
			}
		}
	}
}
