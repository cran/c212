# c212.1a
# Case 2/12 Model c212.1a
# R. Carragher
# Date: 01/10/2013


M <- new.env()

M$Id <- "$Id: c212.1a.R,v 1.16 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.1a <- function(trial.data, sim_type = "SLICE", burnin = 10000, iter = 40000,
	nchains = 3,
	global.sim.params = data.frame(type = c("MH", "SLICE"), param = c("sigma_MH", "w"),
				value = c(0.35,1), control = c(0,6), stringsAsFactors = FALSE),
	sim.params = NULL,
	initial_values = NULL,
	hyper_params = list(mu.gamma.0.0 = 0, tau2.gamma.0.0 = 10,
		mu.theta.0.0 = 0, tau2.theta.0.0 = 10, alpha.gamma.0.0 = 3, beta.gamma.0.0 = 1,
		alpha.theta.0.0 = 3, beta.theta.0.0 = 1, alpha.gamma = 3, beta.gamma = 1,
		alpha.theta = 3, beta.theta = 1))
{
	eot = M_global$EOTdata(M, trial.data, iter, nchains, burnin, initial_values)

	if (is.null(eot)) {
		return(NULL)
	}

	trial.data = eot$trial.data
	cntrl.data = eot$cntrl.data

	M$sim_type <- sim_type

	if (nrow(global.sim.params[global.sim.params$type == sim_type,]) != 1) {
		print("Missing default simulation parametetrs");
		return(NULL)
	}

	M$global.sim.param <- global.sim.params[global.sim.params$type == sim_type,]$value
	M$global.sim.param_ctrl <- global.sim.params[global.sim.params$type ==
																	sim_type,]$control

	if (M$global.sim.param <= 0) {
		print("Invalid simulation parametetr value");
		return(NULL)
	}

	sim.params = M_global$sim_params1a(M, sim.params, sim_type, trial.data, cntrl.data)

	# Hyperparameters
	if (!is.null(hyper_params$mu.gamma.0.0))
		M$mu.gamma.0.0 <- hyper_params$mu.gamma.0.0
	if (!is.null(hyper_params$tau2.gamma.0.0))
		M$tau2.gamma.0.0 <- hyper_params$tau2.gamma.0.0
	if (!is.null(hyper_params$mu.theta.0.0))
		M$mu.theta.0.0 <- hyper_params$mu.theta.0.0
	if (!is.null(hyper_params$tau2.theta.0.0))
		M$tau2.theta.0.0 <- hyper_params$tau2.theta.0.0
	if (!is.null(hyper_params$alpha.gamma.0.0))
		M$alpha.gamma.0.0 <- hyper_params$alpha.gamma.0.0
	if (!is.null(hyper_params$beta.gamma.0.0))
		M$beta.gamma.0.0 <- hyper_params$beta.gamma.0.0
	if (!is.null(hyper_params$alpha.theta.0.0))
		M$alpha.theta.0.0 <- hyper_params$alpha.theta.0.0
	if (!is.null(hyper_params$beta.theta.0.0))
		M$beta.theta.0.0 <- hyper_params$beta.theta.0.0
	if (!is.null(hyper_params$alpha.gamma))
		M$alpha.gamma <- hyper_params$alpha.gamma
	if (!is.null(hyper_params$beta.gamma))
		M$beta.gamma <- hyper_params$beta.gamma
	if (!is.null(hyper_params$alpha.theta))
		M$alpha.theta <- hyper_params$alpha.theta
	if (!is.null(hyper_params$beta.theta))
		M$beta.theta <- hyper_params$beta.theta

	x <- M$x
	y <- M$y
	x[is.na(x)] <- 0
	y[is.na(y)] <- 0
	nc <- M$NC
	nt <- M$NT
	nc[is.na(nc)] <- 0
	nt[is.na(nt)] <- 0

	Ret2 = .Call("c2121a_exec", as.integer(nchains),
					as.integer(burnin), as.integer(iter),
					as.integer(M$NumBodySys), as.integer(M$maxAEs),
					as.integer(M$nAE), M$sim_type, as.numeric(M$global.sim.param),
					as.numeric(M$global.sim.param_ctrl),
					sim.params,
					as.vector(as.integer(t(x))),
					as.vector(as.integer(t(y))),
					as.vector(as.integer(t(nc))),
					as.vector(as.integer(t(nt))),
					as.vector(aperm(M$theta)),
					as.vector(aperm(M$gamma)),
					as.numeric(M$mu.gamma.0.0),
					as.numeric(M$tau2.gamma.0.0),
					as.numeric(M$mu.theta.0.0),
					as.numeric(M$tau2.theta.0.0),
					as.numeric(M$alpha.gamma.0.0),
					as.numeric(M$beta.gamma.0.0),
					as.numeric(M$alpha.theta.0.0),
					as.numeric(M$beta.theta.0.0),
					as.numeric(M$alpha.gamma),
					as.numeric(M$beta.gamma),
					as.numeric(M$alpha.theta),
					as.numeric(M$beta.theta),
					as.vector(as.numeric(M$mu.gamma.0)),
					as.vector(as.numeric(M$tau2.gamma.0)),
					as.vector(as.numeric(M$mu.theta.0)),
					as.vector(as.numeric(M$tau2.theta.0)),
					as.numeric(aperm(M$mu.gamma)),
					as.numeric(aperm(M$mu.theta)),
					as.numeric(aperm(M$sigma2.gamma)),
					as.numeric(aperm(M$sigma2.theta))
				)

	# Getting the samples in this order reduces the maximum memory used.
	mu.gamma.0_samples = .Call("getMuGamma0SamplesAll")
	mu.gamma.0_samples = aperm(mu.gamma.0_samples)

	mu.theta.0_samples = .Call("getMuTheta0SamplesAll")
	mu.theta.0_samples = aperm(mu.theta.0_samples)

	tau2.theta.0_samples = .Call("getTau2Theta0SamplesAll")
	tau2.theta.0_samples = aperm(tau2.theta.0_samples)

	tau2.gamma.0_samples = .Call("getTau2Gamma0SamplesAll")
	tau2.gamma.0_samples = aperm(tau2.gamma.0_samples)

	mu.gamma_samples = .Call("getMuGammaSamplesAll")
	mu.gamma_samples = aperm(mu.gamma_samples)

	mu.theta_samples = .Call("getMuThetaSamplesAll")
	mu.theta_samples = aperm(mu.theta_samples)

	sigma2.gamma_samples = .Call("getSigma2GammaSamplesAll")
	sigma2.gamma_samples = aperm(sigma2.gamma_samples)

	sigma2.theta_samples = .Call("getSigma2ThetaSamplesAll")
	sigma2.theta_samples = aperm(sigma2.theta_samples)

	theta_samples = .Call("getThetaSamplesAll")
	theta_samples = aperm(theta_samples)

	gamma_samples = .Call("getGammaSamplesAll")
	gamma_samples = aperm(gamma_samples)

	theta_acc = .Call("getThetaAcceptAll")
	theta_acc = aperm(theta_acc)

	gamma_acc = .Call("getGammaAcceptAll")
	gamma_acc = aperm(gamma_acc)

	.C("Release")

	print("MCMC fitting complete.")

	model_fit = list(id = M$Id, sim_type = M$sim_type, chains = nchains,
							nBodySys = M$NumBodySys,
							maxAEs = M$maxAEs,
							nAE = M$nAE, AE = M$AE,
							B = M$BodySys,
							burnin = burnin,
							iter = iter,
							mu.gamma.0 = mu.gamma.0_samples,
							mu.theta.0 = mu.theta.0_samples,
							tau2.gamma.0 = tau2.gamma.0_samples,
							tau2.theta.0 = tau2.theta.0_samples,
							mu.gamma = mu.gamma_samples,
							mu.theta = mu.theta_samples,
							sigma2.gamma = sigma2.gamma_samples,
							sigma2.theta = sigma2.theta_samples,
							gamma = gamma_samples,
							theta = theta_samples,
							gamma_acc = gamma_acc,
							theta_acc = theta_acc)

	attr(model_fit, "model") = "1a"

	return(model_fit)
}

M$initVars <- function() {

	# Data Structure
	M$BodySys <- c()
	M$NumBodySys <- NA
	M$maxAEs <- NA
	M$nAE <- c()

	# Trial Event Data
	M$x <- matrix()
	M$y <- matrix()
	M$NT <- matrix()
	M$NC <- matrix()

	# Initial Values
	M$x_init <- array()
	M$y_init <- array()

	# Hyperparameters
	M$mu.gamma.0.0 <- 0
	M$tau2.gamma.0.0 <- 10
	M$mu.theta.0.0 <- 0
	M$tau2.theta.0.0 <- 10
	M$alpha.gamma.0.0 <- 3
	M$beta.gamma.0.0 <- 1
	M$alpha.theta.0.0 <- 3
	M$beta.theta.0.0 <- 1
	M$alpha.gamma <- 3
	M$beta.gamma <- 1
	M$alpha.theta <- 3
	M$beta.theta <- 1

	# Parameters/Simulated values 
	# Stage 3
	M$mu.gamma.0 <- c()
	M$tau2.gamma.0 <- c()
	M$mu.theta.0 <- c()
	M$tau2.theta.0 <- c()


	# Stage 2
	M$mu.gamma <- array()
	M$mu.theta <- array()
	M$sigma2.gamma <- array()
	M$sigma2.theta <- array()

	# Stage 1
	M$theta <- array()
	M$gamma <- array()
}

M$initChain <- function(chain) {
	x_chain <- matrix(NA, nrow = M$NumBodySys, ncol = M$maxAEs)
	y_chain <- matrix(NA, nrow = M$NumBodySys, ncol = M$maxAEs)

	for (b in 1:M$NumBodySys) {
		for (j in 1:M$nAE[b]) {
			ctrl <- seq(0, M$NC[b,j], 1)
			trt <- seq(0, M$NT[b,j], 1)
			x_chain[b, j] <- sample(ctrl, 1)
			y_chain[b, j] <- sample(trt, 1)
		}
	}

	x_chain_init <- x_chain/M$NC
	y_chain_init <- y_chain/M$NT
	x_chain_init[x_chain_init == 0] = 1 / max(M$NC[!is.na(M$NC)])
	y_chain_init[y_chain_init == 0] = 1 / max(M$NT[!is.na(M$NT)])
	x_chain_init[x_chain_init == 1] =
						(max(M$NC[!is.na(M$NC)]) - 1) / max(M$NC[!is.na(M$NC)])
	y_chain_init[y_chain_init == 1] =
						(max(M$NT[!is.na(M$NT)]) - 1) / max(M$NT[!is.na(M$NT)])
	M$x_init[chain,,] <- x_chain_init
	M$y_init[chain,,] <- y_chain_init

	# Stage 1
	M$gamma[chain,,] <- M_global$logit(M$x_init[chain,,])
	M$theta[chain,,] <- M_global$logit(M$y_init[chain,,]) - M$gamma[chain,,]

	# Stage 2
	u <- runif(M$NumBodySys, -50, 50)
	M$mu.gamma[chain, ] <- u
	u <- runif(M$NumBodySys, -50, 50)
	M$mu.theta[chain,] <- u

	u <- runif(M$NumBodySys, 20, 50)
	M$sigma2.gamma[chain,] <- u
	u <- runif(M$NumBodySys, 20, 50)
	M$sigma2.theta[chain,] <- u

	# Stage 3
	u <- runif(1, -50, 50)
	M$mu.gamma.0[chain] <- u
	u <- runif(1, -10, 10)
	M$mu.theta.0[chain] <- u
	u <- runif(1, 20, 50)
	M$tau2.gamma.0[chain] <- u
	u <- runif(1, 20, 50)
	M$tau2.theta.0[chain] <- u
}

M$initialiseChains <- function(initial_values, nchains) {

	M$theta <- array(NA, dim=c(nchains, M$NumBodySys, M$maxAEs))
	M$gamma <- array(NA, dim=c(nchains, M$NumBodySys, M$maxAEs))

	if (is.null(initial_values)) {
		# Default initialisation
		# Initial Values for the simulation - first chain
		M$x_init <- array(NA, dim=c(nchains, M$NumBodySys, M$maxAEs))
		M$y_init <- array(NA, dim=c(nchains, M$NumBodySys, M$maxAEs))
		x_init <- M$x/M$NC
		y_init <- M$y/M$NT
		x_init[x_init == 0] = 1 / max(M$NC[!is.na(M$NC)])
		y_init[y_init == 0] = 1 / max(M$NT[!is.na(M$NT)])
		x_init[x_init == 1] = (max(M$NC[!is.na(M$NC)]) - 1) / max(M$NC[!is.na(M$NC)])
		y_init[y_init == 1] = (max(M$NT[!is.na(M$NT)]) - 1) / max(M$NT[!is.na(M$NT)])
		M$x_init[1,,] <- x_init
		M$y_init[1,,] <- y_init

		# Stage 1
		M$gamma[1,,] <- M_global$logit(M$x_init[1,,])
		M$theta[1,,] <- M_global$logit(M$y_init[1,,]) - M$gamma[1,,]

		# Stage 2
		M$mu.gamma <- array(0, dim=c(nchains, M$NumBodySys))
		M$sigma2.gamma <- array(10, dim=c(nchains, M$NumBodySys))
		M$mu.theta <- array(0, dim=c(nchains, M$NumBodySys))
		M$sigma2.theta <- array(10, dim=c(nchains, M$NumBodySys))

		# Stage 3
		M$mu.gamma.0 <- rep(0, nchains)
		M$tau2.gamma.0 <- rep(10, nchains)
		M$mu.theta.0 <- rep(0, nchains)
		M$tau2.theta.0 <- rep(10, nchains)

		# Initialise any further chains
		if (nchains > 1) {
			for (i in 2:nchains) {
				M$initChain(i)
			}
		}
	}
	else {
		# Use values passed in by caller

		# Stage 1
		for (c in 1:nchains) {
			init_gamma_vals = initial_values$gamma[initial_values$gamma$chain == c,,
																		drop = FALSE]
			init_theta_vals = initial_values$theta[initial_values$theta$chain == c,,
																		drop = FALSE]
			for (b in 1:M$NumBodySys) {
				for (j in 1:M$nAE[b]) {
					bs <- M$BodySys[b]
					val = init_gamma_vals[init_gamma_vals$B == bs &
										init_gamma_vals$AE == M$AE[b,j],,
										drop = FALSE]$value
					M$gamma[c, b, j] = val
					val = init_theta_vals[init_theta_vals$B == bs &
										init_theta_vals$AE == M$AE[b,j],,
										drop = FALSE]$value
					M$theta[c, b, j] = val
				}
			}
		}


		# Stage 2
		M$mu.gamma <- array(NA, dim=c(nchains, M$NumBodySys))
		M$sigma2.gamma <- array(NA, dim=c(nchains, M$NumBodySys))
		M$mu.theta <- array(NA, dim=c(nchains, M$NumBodySys))
		M$sigma2.theta <- array(NA, dim=c(nchains, M$NumBodySys))

		for (c in 1:nchains) {
			for (b in 1:M$NumBodySys) {
				bs <- M$BodySys[b]
				v1 = initial_values$mu.gamma[initial_values$mu.gamma$chain == c &
											initial_values$mu.gamma$B == bs,,
											drop = FALSE]$value
				v2 = initial_values$mu.theta[initial_values$mu.theta$chain == c &
											initial_values$mu.theta$B == bs,,
											drop = FALSE]$value
				v3 = initial_values$sigma2.gamma[initial_values$sigma2.gamma$chain ==
											c & initial_values$sigma2.gamma$B == bs,,
											drop = FALSE]$value
				v4 = initial_values$sigma2.theta[initial_values$sigma2.theta$chain ==
											c & initial_values$sigma2.theta$B == bs,,
											drop = FALSE]$value
				M$mu.gamma[c,b] = v1
				M$mu.theta[c,b]  = v2
				M$sigma2.gamma[c,b]  = v3
				M$sigma2.theta[c,b] = v4
			}
		}

		# Stage 3
		M$mu.gamma.0 <- initial_values$mu.gamma.0
		M$mu.theta.0 <- initial_values$mu.theta.0
		M$tau2.gamma.0 <- initial_values$tau2.gamma.0
		M$tau2.theta.0 <- initial_values$tau2.theta.0
	}
}
