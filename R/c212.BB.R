# c212.BB
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 05/08/2013
#
# In coding terms the main difference between this and c212.1a is a change to
# sample_theta and the addition of new parameter pi, alpha_pi, beta_pi and their priors
# The implementation of sample_theta follows [BB2004] - which has a typo in its
# algorithm for sampling theta for the case theta(C) != 0, theta(0) = 0

M1 <- new.env()

M1$Id <- "$Id: c212.BB.R,v 1.14 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.BB <- function(trial.data, burnin = 20000, iter = 60000, nchains = 3,
	theta_algorithm = "MH", sim_type = "SLICE",
	global.sim.params = data.frame(type = c("MH", "MH", "MH", "MH",
							"SLICE", "SLICE", "SLICE"),
		param = c("sigma_MH_alpha", "sigma_MH_beta", "sigma_MH_gamma", "sigma_MH_theta",
				"w_alpha", "w_beta", "w_gamma"),
		value = c(3, 3, 0.2, 0.2, 1, 1, 1), control = c(0, 0, 0, 0, 6, 6, 6),
		stringsAsFactors = FALSE),
	sim.params = NULL,
	initial_values = NULL, hyper_params = list(mu.gamma.0.0 = 0,
		tau2.gamma.0.0 = 10, mu.theta.0.0 = 0, tau2.theta.0.0 = 10, alpha.gamma.0.0 = 3,
		beta.gamma.0.0 = 1, alpha.theta.0.0 = 3,beta.theta.0.0 = 1, alpha.gamma = 3,
		beta.gamma = 1, alpha.theta = 3, beta.theta = 1,
		lambda.alpha = 1.0, lambda.beta = 1.0),
	global.pm.weight = 0.5,
	pm.weights = NULL,
	adapt_params = data.frame(min_w = 0.25, chains = 3, burnin = 20000 , iter = 40000),
	adapt_phase=0)
{
	eot = M_global$EOTdata(M1, trial.data, iter, nchains, burnin, initial_values)

	if (is.null(eot)) {
		return(NULL)
	}

	trial.data = eot$trial.data
	cntrl.data = eot$cntrl.data

	M1$Algo <- theta_algorithm

	M1$sim_type <- sim_type

	if (nrow(global.sim.params[global.sim.params$type == sim_type,]) == 0) {
		print("Missing simulation parameters");
		return(NULL)
	}

	if (!all(global.sim.params$value > 0)) {
			print("Invalid simulation parameter value");
			return(NULL)
	}

	M1$global.sim.params <- global.sim.params

	sp = M_global$sim_paramsBB(M1, sim.params, pm.weights, sim_type, trial.data, cntrl.data)

	sim.params = sp$sim.params
	pm.weights = sp$pm.weights

	# Hyperparameters
	M1$mu.gamma.0.0 <- hyper_params$mu.gamma.0.0
	M1$tau2.gamma.0.0 <- hyper_params$tau2.gamma.0.0
	M1$mu.theta.0.0 <- hyper_params$mu.theta.0.0
	M1$tau2.theta.0.0 <- hyper_params$tau2.theta.0.0
	M1$alpha.gamma.0.0 <- hyper_params$alpha.gamma.0.0
	M1$beta.gamma.0.0 <- hyper_params$beta.gamma.0.0
	M1$alpha.theta.0.0 <- hyper_params$alpha.theta.0.0
	M1$beta.theta.0.0 <- hyper_params$beta.theta.0.0
	M1$alpha.gamma <- hyper_params$alpha.gamma
	M1$beta.gamma <- hyper_params$beta.gamma
	M1$alpha.theta <- hyper_params$alpha.theta
	M1$beta.theta <- hyper_params$beta.theta

	# BB2004 parameters
	M1$lambda.alpha <- hyper_params$lambda.alpha
	M1$lambda.beta <- hyper_params$lambda.beta

	x <- M1$x
	y <- M1$y
	x[is.na(x)] <- 0
	y[is.na(y)] <- 0
	nc <- M1$NC
	nt <- M1$NT
	nc[is.na(nc)] <- 0
	nt[is.na(nt)] <- 0

	algo <- 1
	if (M1$Algo == "MH") {
		algo <- 1;
	} else if (M1$Algo == "MH_Adapt") {
		algo <- 2;
	} else if (M1$Algo == "MIS_Adapt") {
		algo <- 3;
	} else if (M1$Algo == "Indep") {
		algo <- 4;
	} else {
		algo <- 1;
	}

	Ret2 = .Call("c212BB_exec", as.integer(nchains),
					as.integer(burnin), as.integer(iter),
					as.integer(M1$NumBodySys), as.integer(M1$maxAEs),
					as.integer(M1$nAE),
					as.vector(as.integer(t(x))), as.vector(as.integer(t(y))), 
					as.vector(as.integer(t(nc))),
					as.vector(as.integer(t(nt))),
					as.vector(aperm(M1$theta)),
					as.vector(aperm(M1$gamma)),
					as.numeric(M1$mu.gamma.0.0),
					as.numeric(M1$tau2.gamma.0.0),
					as.numeric(M1$mu.theta.0.0),
					as.numeric(M1$tau2.theta.0.0),
					as.numeric(M1$alpha.gamma.0.0),
					as.numeric(M1$beta.gamma.0.0),
					as.numeric(M1$alpha.theta.0.0),
					as.numeric(M1$beta.theta.0.0),
					as.numeric(M1$alpha.gamma),
					as.numeric(M1$beta.gamma),
					as.numeric(M1$alpha.theta),
					as.numeric(M1$beta.theta),
					as.vector(as.numeric(M1$mu.gamma.0)),
					as.vector(as.numeric(M1$tau2.gamma.0)),
					as.vector(as.numeric(M1$mu.theta.0)),
					as.vector(as.numeric(M1$tau2.theta.0)),
					as.numeric(aperm(M1$mu.gamma)),
					as.numeric(aperm(M1$mu.theta)),
					as.numeric(aperm(M1$sigma2.gamma)),
					as.numeric(aperm(M1$sigma2.theta)),
					as.vector(as.numeric(M1$alpha.pi)),
					as.vector(as.numeric(M1$beta.pi)),
					as.numeric(M1$lambda.alpha),
					as.numeric(M1$lambda.beta),
					as.numeric(aperm(M1$pi)),
					as.integer(algo),
					as.integer(adapt_phase),
					M1$sim_type,
					M1$global.sim.params,
					sim.params,
					as.numeric(global.pm.weight),
					pm.weights,
					as.numeric(adapt_params$min_w),
					as.integer(adapt_params$chains),
					as.integer(adapt_params$burnin),
					as.integer(adapt_params$iter))

	alpha.pi_samples = .Call("getAlphaPiSamplesAll")
	alpha.pi_samples = aperm(alpha.pi_samples)

	beta.pi_samples = .Call("getBetaPiSamplesAll")
	beta.pi_samples = aperm(beta.pi_samples)

	alpha.pi_acc <- .Call("getAlphaPiAcceptAll")
	beta.pi_acc <- .Call("getBetaPiAcceptAll")

	mu.gamma.0_samples = .Call("getMuGamma0SamplesAll")
	mu.gamma.0_samples = aperm(mu.gamma.0_samples)

	mu.theta.0_samples = .Call("getMuTheta0SamplesAll")
	mu.theta.0_samples = aperm(mu.theta.0_samples)

	tau2.gamma.0_samples = .Call("getTau2Gamma0SamplesAll")
	tau2.gamma.0_samples = aperm(tau2.gamma.0_samples)

	tau2.theta.0_samples = .Call("getTau2Theta0SamplesAll")
	tau2.theta.0_samples = aperm(tau2.theta.0_samples)

	mu.gamma_samples = .Call("getMuGammaSamplesAll")
	mu.gamma_samples = aperm(mu.gamma_samples)

	mu.theta_samples = .Call("getMuThetaSamplesAll")
	mu.theta_samples = aperm(mu.theta_samples)

	sigma2.gamma_samples = .Call("getSigma2GammaSamplesAll")
	sigma2.gamma_samples = aperm(sigma2.gamma_samples)

	sigma2.theta_samples = .Call("getSigma2ThetaSamplesAll")
	sigma2.theta_samples = aperm(sigma2.theta_samples)

	pi_samples = .Call("getPiSamplesAll")
	pi_samples = aperm(pi_samples)

	theta_samples = .Call("getThetaSamplesAll")
	theta_samples = aperm(theta_samples)

	gamma_samples = .Call("getGammaSamplesAll")
	gamma_samples = aperm(gamma_samples)

	theta_acc = .Call("getThetaAcceptAll")
	theta_acc = aperm(theta_acc)
	
	gamma_acc = .Call("getGammaAcceptAll")
	gamma_acc = aperm(gamma_acc)

	theta_zero_prop = .Call("getThetaZeroPropAll")
	theta_zero_prop = aperm(theta_zero_prop)

	theta_zero_acc = .Call("getThetaZeroAcceptAll")
	theta_zero_acc = aperm(theta_zero_acc)

	.C("Release")

	print("MCMC fitting complete.")

	model_fit =  list(id = M1$Id, theta_alg = M1$Algo, sim_type = M1$sim_type,
							chains = nchains, nBodySys = M1$NumBodySys, maxAEs = M1$maxAEs,
							nAE = M1$nAE, AE = M1$AE, B = M1$BodySys,
							burnin = burnin, iter = iter,
							mu.gamma.0 = mu.gamma.0_samples,
							mu.theta.0 = mu.theta.0_samples,
							tau2.gamma.0 = tau2.gamma.0_samples,
							tau2.theta.0 = tau2.theta.0_samples,
							mu.gamma = mu.gamma_samples,
							mu.theta = mu.theta_samples,
							sigma2.gamma = sigma2.gamma_samples,
							sigma2.theta = sigma2.theta_samples,
							pi = pi_samples,
							alpha.pi = alpha.pi_samples,
							beta.pi = beta.pi_samples,
							alpha.pi_acc = alpha.pi_acc,
							beta.pi_acc = beta.pi_acc,
							gamma = gamma_samples,
							theta = theta_samples,
							gamma_acc = gamma_acc,
							theta_acc = theta_acc,
							theta_zero_prop = theta_zero_prop,
							theta_zero_acc = theta_zero_acc)

		attr(model_fit, "model") = "BB"

		return(model_fit)
}

M1$initVars <- function() {

	# Data Structure
	M1$BodySys <- c()
	M1$NumBodySys <- NA
	M1$maxAEs <- NA
	M1$nAE <- c()

	# Trial Event Data
	M1$x <- matrix()
	M1$y <- matrix()
	M1$NT <- matrix()
	M1$NC <- matrix()

	# Initial Values
	M1$x_init <- array()
	M1$y_init <- array()

	# Hyperparameters
	M1$mu.gamma.0.0 <- NA
	M1$tau2.gamma.0.0 <- NA
	M1$mu.theta.0.0 <- NA
	M1$tau2.theta.0.0 <- NA
	M1$alpha.gamma.0.0 <- NA
	M1$beta.gamma.0.0 <- NA
	M1$alpha.theta.0.0 <- NA
	M1$beta.theta.0.0 <- NA
	M1$alpha.gamma <- NA
	M1$beta.gamma <- NA
	M1$alpha.theta <- NA
	M1$beta.theta <- NA

	# Parameters/Simulated values 
	# Stage 3
	M1$mu.gamma.0 <- c()
	M1$tau2.gamma.0 <- c()
	M1$mu.theta.0 <- c()
	M1$tau2.theta.0 <- c()

	# Stage 2
	M1$mu.gamma <- array()
	M1$mu.theta <- array()
	M1$sigma2.gamma <- array()
	M1$sigma2.theta <- array()

	# Stage 1
	M1$theta <- array()
	M1$gamma <- array()

	# BB2004 parameters
	M1$lambda.alpha <- NA
	M1$lambda.beta <- NA
	M1$alpha.pi <- NA
	M1$beta.pi <- NA
	M1$pi <- NA
}

M1$initChain <- function(chain) {
	x_chain <- matrix(NA, nrow = M1$NumBodySys, ncol = M1$maxAEs)
	y_chain <- matrix(NA, nrow = M1$NumBodySys, ncol = M1$maxAEs)

	for (b in 1:M1$NumBodySys) {
		for (j in 1:M1$nAE[b]) {
			ctrl <- seq(0, M1$NC[b,j], 1)
			trt <- seq(0, M1$NT[b,j], 1)
			x_chain[b, j] <- sample(ctrl, 1)
			y_chain[b, j] <- sample(trt, 1)
		}
	}

	x_chain_init <- x_chain/M1$NC
	y_chain_init <- y_chain/M1$NT
	x_chain_init[x_chain_init == 0] = 1 / max(M1$NC[!is.na(M1$NC)])
	y_chain_init[y_chain_init == 0] = 1 / max(M1$NT[!is.na(M1$NT)])
	x_chain_init[x_chain_init == 1] =
							(max(M1$NC[!is.na(M1$NC)]) - 1) / max(M1$NC[!is.na(M1$NC)])
	y_chain_init[y_chain_init == 1] =
							(max(M1$NT[!is.na(M1$NT)]) - 1) / max(M1$NT[!is.na(M1$NT)])
	M1$x_init[chain,,] <- x_chain_init
	M1$y_init[chain,,] <- y_chain_init

	d <- M1$x_init[chain,,][!is.na(M1$x_init[chain,,])]
	d <- d[d >= 1]
	e <- M1$y_init[chain,,][!is.na(M1$y_init[chain,,])]
	e <- e[e <= 0]
	if (length(d) > 0) {
		browser()
	}
	if (length(e) > 0) {
		browser()
	}

	# Stage 1
	M1$gamma[chain,,] <- M_global$logit(M1$x_init[chain,,])
	M1$theta[chain,,] <- M_global$logit(M1$y_init[chain,,]) - M1$gamma[chain,,]

	# Stage 2
	u <- runif(M1$NumBodySys, -50, 50)
	M1$mu.gamma[chain, ] <- u
	u <- runif(M1$NumBodySys, -50, 50)
	M1$mu.theta[chain,] <- u

	u <- runif(M1$NumBodySys, 20, 50)
	M1$sigma2.gamma[chain,] <- u
	u <- runif(M1$NumBodySys, 20, 50)
	M1$sigma2.theta[chain,] <- u

	# Stage 3
	u <- runif(1, -50, 50)
	M1$mu.gamma.0[chain] <- u
	u <- runif(1, -50, 50)
	M1$mu.theta.0[chain] <- u
	u <- runif(1, 20, 50)
	M1$tau2.gamma.0[chain] <- u
	u <- runif(1, 20, 50)
	M1$tau2.theta.0[chain] <- u

	# BB2004 parameters
	u <- runif(M1$NumBodySys, 0, 1)
	M1$pi[chain, ] <- u
	u <- runif(1, 1.25,100)
	M1$alpha.pi[chain] <- u
	u <- runif(1, 1.25,100)
	M1$beta.pi[chain] <- u
}

M1$initialiseChains <- function(initial_values, nchains) {

	M1$theta <- array(NA, dim=c(nchains, M1$NumBodySys, M1$maxAEs))
	M1$gamma <- array(NA, dim=c(nchains, M1$NumBodySys, M1$maxAEs))

	if (is.null(initial_values)) {

		M1$x_init <- array(NA, dim=c(nchains, M1$NumBodySys, M1$maxAEs))
		M1$y_init <- array(NA, dim=c(nchains, M1$NumBodySys, M1$maxAEs))
		x_init <- M1$x/M1$NC
		y_init <- M1$y/M1$NT
		x_init[x_init == 0] = 1 / max(M1$NC[!is.na(M1$NC)])
		y_init[y_init == 0] = 1 / max(M1$NT[!is.na(M1$NT)])
		x_init[x_init == 1] = (max(M1$NC[!is.na(M1$NC)]) - 1) / max(M1$NC[!is.na(M1$NC)])
		y_init[y_init == 1] = (max(M1$NT[!is.na(M1$NT)]) - 1) / max(M1$NT[!is.na(M1$NT)])
		M1$x_init[1,,] <- x_init
		M1$y_init[1,,] <- y_init


		# Stage 1
		M1$gamma[1,,] <- M_global$logit(M1$x_init[1,,])
		M1$theta[1,,] <- M_global$logit(M1$y_init[1,,]) - M1$gamma[1,,]

		# Stage 2
		M1$mu.gamma <- array(0, dim=c(nchains, M1$NumBodySys))
		M1$sigma2.gamma <- array(10, dim=c(nchains, M1$NumBodySys))
		M1$mu.theta <- array(0, dim=c(nchains, M1$NumBodySys))
		M1$sigma2.theta <- array(10, dim=c(nchains, M1$NumBodySys))

		# BB2004 parameters
		M1$pi <- array(0.5, dim=c(nchains, M1$NumBodySys))

		# Stage 3
		M1$mu.gamma.0 <- rep(0, nchains)
		M1$tau2.gamma.0 <- rep(10, nchains)
		M1$mu.theta.0 <- rep(0, nchains)
		M1$tau2.theta.0 <- rep(10, nchains)

		# BB2004 parameter
		M1$alpha.pi <- rep(1.5, nchains)
		M1$beta.pi <- rep(1.5, nchains)

		# Initialise any further chains
		if (nchains > 1) {
			for (i in 2:nchains) {
				M1$initChain(i)
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
			for (b in 1:M1$NumBodySys) {
				for (j in 1:M1$nAE[b]) {
					bs <- M1$BodySys[b]
					val = init_gamma_vals[init_gamma_vals$B == bs &
										init_gamma_vals$AE == M1$AE[b,j],,
																	drop = FALSE]$value
					M1$gamma[c, b, j] = val
					val = init_theta_vals[init_theta_vals$B == bs &
										init_theta_vals$AE == M1$AE[b,j],,
																	drop = FALSE]$value
					M1$theta[c, b, j] = val
				}
			}
		}

        # Stage 2
		M1$mu.gamma <- array(NA, dim=c(nchains, M1$NumBodySys))
		M1$sigma2.gamma <- array(NA, dim=c(nchains, M1$NumBodySys))
		M1$mu.theta <- array(NA, dim=c(nchains, M1$NumBodySys))
		M1$sigma2.theta <- array(NA, dim=c(nchains, M1$NumBodySys))
		M1$pi <- array(NA, dim=c(nchains, M1$NumBodySys))


		# BB2004 parameters
		M1$pi <- array(0.5, dim=c(nchains, M1$NumBodySys))
		for (c in 1:nchains) {
			for (b in 1:M1$NumBodySys) {
				bs <- M1$BodySys[b]
				v1 = initial_values$mu.gamma[initial_values$mu.gamma$chain == c &
                                            initial_values$mu.gamma$B == bs,,
											drop = FALSE]$value
				v2 = initial_values$mu.theta[initial_values$mu.theta$chain ==
											c & initial_values$mu.theta$B == bs,,
											drop = FALSE]$value
				v3 = initial_values$sigma2.gamma[initial_values$sigma2.gamma$chain ==
											c & initial_values$sigma2.gamma$B == bs,,
											drop = FALSE]$value
				v4 = initial_values$sigma2.theta[initial_values$sigma2.theta$chain ==
											c & initial_values$sigma2.theta$B == bs,,
											drop = FALSE]$value
				v5 = initial_values$pi[initial_values$pi$chain ==
											c & initial_values$pi$B == bs,,
											drop = FALSE]$value
				M1$mu.gamma[c,b] = v1
				M1$mu.theta[c,b]  = v2
				M1$sigma2.gamma[c,b]  = v3
				M1$sigma2.theta[c,b] = v4
				M1$pi[c, b] = v5
            }
		}

		M1$pi <- array(0.5, dim=c(nchains, M1$NumBodySys))

		# Stage 3
		M1$mu.gamma.0 <- initial_values$mu.gamma.0
		M1$mu.theta.0 <- initial_values$mu.theta.0
		M1$tau2.gamma.0 <- initial_values$tau2.gamma.0
		M1$tau2.theta.0 <- initial_values$tau2.theta.0

		M1$alpha.pi <- initial_values$alpha.pi
		M1$beta.pi <- initial_values$beta.pi
	}
}
