c212.gen.initial.values = function(trial.data, nchains = 3, model = "1a", hier = 3, level = 0) {

	if (nchains < 1)
		return(NULL)

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors = FALSE)
	}

	facs <- sapply(trial.data, is.factor)
	trial.data[facs] <- lapply(trial.data[facs], as.character)

	initial_values = NULL

	if ("Total" %in% names(trial.data)) {
		initial_values = c212.eot.gen.initial.values(trial.data, nchains, model)
	}
	else {
		initial_values = c212.interim.gen.initial.values(trial.data, nchains, model, hier, level)
	}

	initial_values
}

c212.eot.gen.initial.values = function(trial.data, nchains, model) {

	initial_values = NULL

	# Check the correct columns are defined
	cols = c("B", "AE", "Count", "Group", "Total")
	if (M_global$checkCols(cols, trial.data)) {
		print("Missing columns");
		return(NULL)
	}

	ordered.data <- trial.data[order(trial.data$B, trial.data$AE, trial.data$Group),,
                                                                            drop=FALSE]

	# Check both contain the same data
	cntrl.data = ordered.data[ordered.data$Group == 1,]
	treat.data = ordered.data[ordered.data$Group == 2,]

	if(!identical(cntrl.data$B, treat.data$B)) {
		print("Mismatced body-system data");
		return(NULL)
	}

	if(!identical(cntrl.data$AE, treat.data$AE)) {
		print("Mismatced adverse event data");
		return(NULL)
	}


	n = nrow(cntrl.data)
	B = as.character(unique(ordered.data$B))
	numBS = length(B)

	# The data structures
	# First chain used fixed values based on the original hyper-params
	mu.gamma.0 <- rep(0, nchains)
	mu.theta.0 <- rep(0, nchains)
	tau2.gamma.0 <- rep(10, nchains)
	tau2.theta.0 <- rep(10, nchains)
	alpha.pi <- rep(1.5, nchains)
	beta.pi <- rep(1.5, nchains)

	mu.gamma <- data.frame(chain = numeric(nchains*numBS), B = character(nchains*numBS),
											value = numeric(nchains*numBS), stringsAsFactors = FALSE)
	mu.theta <- data.frame(chain = numeric(nchains*numBS), B = character(nchains*numBS),
											value = numeric(nchains*numBS), stringsAsFactors = FALSE)
	sigma2.gamma <- data.frame(chain = numeric(nchains*numBS), B = character(nchains*numBS),
											value = numeric(nchains*numBS), stringsAsFactors = FALSE)
	sigma2.theta <- data.frame(chain = numeric(nchains*numBS), B = character(nchains*numBS),
											value = numeric(nchains*numBS), stringsAsFactors = FALSE)
	pi <- data.frame(chain = numeric(nchains*numBS), B = character(nchains*numBS),
											value = numeric(nchains*numBS), stringsAsFactors = FALSE)

	mu.gamma[1:numBS,]$chain = 1
	mu.gamma[1:numBS,]$B = B
	mu.gamma[1:numBS,]$value = 0
	mu.theta[1:numBS,]$chain = 1
	mu.theta[1:numBS,]$B = B
	mu.theta[1:numBS,]$value = 0
	sigma2.gamma[1:numBS,]$chain = 1
	sigma2.gamma[1:numBS,]$B = B
	sigma2.gamma[1:numBS,]$value = 10
	sigma2.theta[1:numBS,]$chain = 1
	sigma2.theta[1:numBS,]$B = B
	sigma2.theta[1:numBS,]$value = 10

	pi[1:numBS,]$chain = 1
	pi[1:numBS,]$B = B
	pi[1:numBS,]$value = 0.5

	# initial values for gamma/theta
	data = cntrl.data[,colnames(cntrl.data) %in% c("B", "AE")]
	chain = rep(1, n)
	
	t = data
	for (c in 2:nchains) {
		t = rbind(t, data)
		chain = c(chain, rep(c, n))
	}
	t = cbind(t, chain)

	gamma = t
	theta = t

	val = rep(0, nrow(gamma))
	gamma = cbind(gamma, val)
	theta = cbind(theta, val)

	x_chain <- rep(0,n)
	y_chain <- rep(0,n)

	# First chain - derive the theta/gamma from the data
	x_chain = cntrl.data$Count/cntrl.data$Total
	y_chain = treat.data$Count/treat.data$Total

	Nc = max(cntrl.data$Total)
	Nt = max(treat.data$Total)

	x_chain[x_chain == 0] <- 1/Nc
  	y_chain[y_chain == 0] <- 1/Nt

	x_chain[x_chain == 1] <- (Nc - 1)/Nc
	y_chain[y_chain == 1] <- (Nt - 1)/Nt

	ga = log(x_chain/(1 - x_chain))
	th = log(y_chain/(1 - y_chain)) - ga

	gamma$value[1:n] = ga
	theta$value[1:n] = th

	if (nchains > 1) {
		for (c in 2:nchains) {
			# gamma/theta
			for (a in 1:n) {
				x_chain[a] <- (sample(0:cntrl.data[a,]$Total, 1))/cntrl.data[a,]$Total
				y_chain[a] <- (sample(0:treat.data[a,]$Total, 1))/treat.data[a,]$Total
			}

			x_chain[x_chain == 0] <- 1/Nc
   			y_chain[y_chain == 0] <- 1/Nt

			x_chain[x_chain == 1] <- (Nc - 1)/Nc
			y_chain[y_chain == 1] <- (Nt - 1)/Nt

			ga = log(x_chain/(1-x_chain))
			th = log(y_chain/(1-y_chain)) - ga

			gamma$value[1:n + a * (c - 1)] = ga
			theta$value[1:n + a * (c - 1)] = th
	
			# mu.gamma etc.
			u <- runif(numBS, -50, 50)
			mu.gamma[1:numBS + numBS * (c -1),]$chain = c
			mu.gamma[1:numBS + numBS * (c -1),]$value = u
			mu.gamma[1:numBS + numBS * (c -1),]$B = B
			u <- runif(numBS, -50, 50)
			mu.theta[1:numBS + numBS * (c -1),]$chain = c
			mu.theta[1:numBS + numBS * (c -1),]$value = u
			mu.theta[1:numBS + numBS * (c -1),]$B = B
			u <- runif(numBS, 20, 50)
			sigma2.gamma[1:numBS + numBS * (c -1),]$chain = c
			sigma2.gamma[1:numBS + numBS * (c -1),]$value = u
			sigma2.gamma[1:numBS + numBS * (c -1),]$B = B
			u <- runif(numBS, 20, 50)
			sigma2.theta[1:numBS + numBS * (c -1),]$chain = c
			sigma2.theta[1:numBS + numBS * (c -1),]$value = u
			sigma2.theta[1:numBS + numBS * (c -1),]$B = B
	
			u <- runif(numBS, 0 , 1)
			pi[1:numBS + numBS * (c -1),]$chain = c
			pi[1:numBS + numBS * (c -1),]$value = u
			pi[1:numBS + numBS * (c -1),]$B = B
		
	
			# mu.gamma.0 etc.
			u <- runif(1, -50, 50)
			mu.gamma.0[c] <- u
			u <- runif(1, -50, 50)
			mu.theta.0[c] <- u
			u <- runif(1, 5, 20)
			tau2.gamma.0[c] <- u
			u <- runif(1, 5, 20)
			tau2.theta.0[c] <- u
	
			u <- runif(1, 1.25,100)
			alpha.pi[c] <- u
			u <- runif(1, 1.25,100)
			beta.pi[c] <- u
		}
	}

	initial_values = list(gamma = gamma, theta = theta, mu.gamma = mu.gamma,
						mu.theta = mu.theta, sigma2.gamma = sigma2.gamma,
						sigma2.theta = sigma2.theta, mu.gamma.0 = mu.gamma.0,
						mu.theta.0 = mu.theta.0, tau2.gamma.0 = tau2.gamma.0,
						tau2.theta.0 = tau2.theta.0)

	if (model == "BB") {
		bb_init = list(pi = pi, alpha.pi = alpha.pi, beta.pi = beta.pi)
		initial_values = c(initial_values, bb_init)
	}

	initial_values
}

c212.interim.gen.initial.values = function(trial.data, nchains, model, hier, level) {

	initial_values = NULL

	# Check the correct columns are defined
	cols = c("B", "AE", "Count", "Group", "Interval", "Exposure", "I_index")
	if (M_global$checkCols(cols, trial.data)) {
		print("Missing columns");
		return(NULL)
	}

   # Order by body-system, adverse event, interval and group
    ordered.data <- trial.data[order(trial.data$I_index, trial.data$B, trial.data$AE, trial.data$Group),, drop=FALSE]

    cntrl.data <- ordered.data[ordered.data$Group == 1, ]
    treat.data <- ordered.data[ordered.data$Group == 2, ]

    # The data is ordered so a straight comparison is possible
	if(!identical(cntrl.data$B, treat.data$B)) {
		print("Mismatced body-system data");
		return(NULL)
	}
	if(!identical(cntrl.data$AE, treat.data$AE)) {
		print("Mismatced adverse event data");
		return(NULL)
	}
	if(!identical(cntrl.data$Interval, treat.data$Interval)) {
		print("Mismatced interval data");
		return(NULL)
	}

	Intervals = unique(cntrl.data$Interval)
	numIntervals = length(Intervals)
	numB = c()

	for (i in 1:numIntervals) {
		c = cntrl.data[cntrl.data$Interval == Intervals[i],]
		b = unique(c$B)
		numB[i] = length(b)
	}

	maxBs = max(numB)
	B = array(NA, dim = c(numIntervals, maxBs))
	for (i in 1:numIntervals) {
		c = cntrl.data[cntrl.data$Interval == Intervals[i],]
		b = unique(c$B)
		B[i, 1:length(b)] = b
	}

	nAE = array(0, dim = c(numIntervals, maxBs))
	for (i in 1:numIntervals) {
		n = numB[i]
		for (b in 1:n) {
			nAE[i, b] = length(unique(cntrl.data[cntrl.data$B == B[i, b], ]$AE))
		}
	}

    maxAEs <- max(nAE)

	n = nrow(cntrl.data)

	data = cntrl.data[,colnames(cntrl.data) %in% c("Interval", "I_index", "B", "AE")]
	chain = rep(1, n)
	
	t = data
	if (nchains > 1) {
		for (c in 2:nchains) {
			t = rbind(t, data)
			chain = c(chain, rep(c, n))
		}
	}
	t = cbind(t, chain)

	gamma = t
	theta = t

	value = rep(0, nrow(gamma))
	theta = cbind(theta, value)
	gamma = cbind(gamma, value)

	# First chain - derive the theta/gamma from the data
	x_chain = cntrl.data$Count/cntrl.data$Exposure
	y_chain = treat.data$Count/treat.data$Exposure

	ga = log(x_chain)
	th = log(y_chain) - ga

	ga[is.infinite(ga)] = -10 
	th[is.infinite(th)] = -10 

	gamma$value[1:n] = ga
	theta$value[1:n] = th

	if (nchains > 1) {
		for (c in 2:nchains) {
			u = runif(n, -10, 10)
			gamma$value[1:n + n * (c - 1)] = u
			u = runif(n, -10, 10)
			theta$value[1:n + n * (c - 1)] = u
		}
	}

	if (level == 1) {
		sz = length(B[1,])
		mu.gamma <- data.frame(chain = numeric(nchains*sz), B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)
		mu.theta <- data.frame(chain = numeric(nchains*sz), B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)
		sigma2.gamma <- data.frame(chain = numeric(nchains*sz), B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)
		sigma2.theta <- data.frame(chain = numeric(nchains*sz), B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)

		pi <- data.frame(chain = numeric(nchains*sz), B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)

		mu.gamma[1:sz, ]$chain = 1
		mu.gamma[1:sz, ]$B = B[1,]
		mu.gamma[1:sz, ]$value = 0
		mu.theta[1:sz, ]$chain = 1
		mu.theta[1:sz, ]$B = B[1,]
		mu.theta[1:sz, ]$value = 0
		sigma2.gamma[1:sz, ]$chain = 1
		sigma2.gamma[1:sz, ]$B = B[1,]
		sigma2.gamma[1:sz, ]$value = 10
		sigma2.theta[1:sz, ]$chain = 1
		sigma2.theta[1:sz, ]$B = B[1,]
		sigma2.theta[1:sz, ]$value = 10
		pi[1:sz, ]$chain = 1
		pi[1:sz, ]$B = B[1,]
		pi[1:sz, ]$value = 0.5

		if (nchains > 1) {
			offset = sz
			for (c in 2:nchains) {
				mu.gamma[offset + 1:sz, ]$chain = c
				mu.gamma[offset + 1:sz, ]$B = B[1,]
				mu.gamma[offset + 1:sz, ]$value = runif(sz, -10, 10)
				mu.theta[offset + 1:sz, ]$chain = c
				mu.theta[offset + 1:sz, ]$B = B[1,]
				mu.theta[offset + 1:sz, ]$value = runif(sz, -10, 10)
				sigma2.gamma[offset + 1:sz, ]$chain = c
				sigma2.gamma[offset + 1:sz, ]$B = B[1,]
				sigma2.gamma[offset + 1:sz, ]$value = runif(sz, 5, 20)
				sigma2.theta[offset + 1:sz, ]$chain = c
				sigma2.theta[offset + 1:sz, ]$B = B[1,]
				sigma2.theta[offset + 1:sz, ]$value = runif(sz, 5, 20)
	
				pi[offset + 1:sz, ]$chain = c
				pi[offset + 1:sz, ]$B = B[1,]
				pi[offset + 1:sz, ]$value = runif(sz, 0, 1)
				offset = offset + sz
			}
		}
	}
	else {

		# 1a hier3 lev 0, 2
		sz = sum(numB)
		mu.gamma <- data.frame(chain = numeric(nchains*sz), Interval = character(nchains*sz),
											B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)
		mu.theta <- data.frame(chain = numeric(nchains*sz), Interval = character(nchains*sz),
											B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)
		sigma2.gamma <- data.frame(chain = numeric(nchains*sz), Interval = character(nchains*sz),
											B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)
		sigma2.theta <- data.frame(chain = numeric(nchains*sz), Interval = character(nchains*sz),
											B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)
		pi <- data.frame(chain = numeric(nchains*sz), Interval = character(nchains*sz),
											B = character(nchains*sz),
											value = numeric(nchains*sz), stringsAsFactors = FALSE)

		# First chain used fixed values based on the original hyper-params

		offset = 0
		for (i in 1:numIntervals) {
			mu.gamma[offset + 1:numB[i],]$chain = 1
			mu.gamma[offset + 1:numB[i],]$Interval = Intervals[i]
			mu.gamma[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
			mu.gamma[offset + 1:numB[i],]$value = 0
			mu.theta[offset + 1:numB[i],]$chain = 1
			mu.theta[offset + 1:numB[i],]$Interval = Intervals[i]
			mu.theta[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
			mu.theta[offset + 1:numB[i],]$value = 0
			sigma2.gamma[offset + 1:numB[i],]$chain = 1
			sigma2.gamma[offset + 1:numB[i],]$Interval = Intervals[i]
			sigma2.gamma[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
			sigma2.gamma[offset + 1:numB[i],]$value = 0
			sigma2.theta[offset + 1:numB[i],]$chain = 1
			sigma2.theta[offset + 1:numB[i],]$Interval = Intervals[i]
			sigma2.theta[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
			sigma2.theta[offset + 1:numB[i],]$value = 0

			pi[offset + 1:numB[i],]$chain = 1
			pi[offset + 1:numB[i],]$Interval = Intervals[i]
			pi[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
			pi[offset + 1:numB[i],]$value = 0.5

			offset = offset + numB[i]
		}

		if (nchains > 1) {
			for (c in 2:nchains) {
				for (i in 1:numIntervals) {
					mu.gamma[offset + 1:numB[i],]$chain = c
					mu.gamma[offset + 1:numB[i],]$Interval = Intervals[i]
					mu.gamma[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
					mu.gamma[offset + 1:numB[i],]$value = runif(numB[i], -10, 10)
					mu.theta[offset + 1:numB[i],]$chain = c
					mu.theta[offset + 1:numB[i],]$Interval = Intervals[i]
					mu.theta[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
					mu.theta[offset + 1:numB[i],]$value = runif(numB[i], -10, 10)
					sigma2.gamma[offset + 1:numB[i],]$chain = c
					sigma2.gamma[offset + 1:numB[i],]$Interval = Intervals[i]
					sigma2.gamma[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
					sigma2.gamma[offset + 1:numB[i],]$value = runif(numB[i], 5, 20)
					sigma2.theta[offset + 1:numB[i],]$chain = c
					sigma2.theta[offset + 1:numB[i],]$Interval = Intervals[i]
					sigma2.theta[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
					sigma2.theta[offset + 1:numB[i],]$value = runif(numB[i], 5, 20)
	
					pi[offset + 1:numB[i],]$chain = c
					pi[offset + 1:numB[i],]$Interval = Intervals[i]
					pi[offset + 1:numB[i],]$B = B[i, 1:numB[i]]
					pi[offset + 1:numB[i],]$value = runif(numB[i], 0, 1)
	
					offset = offset + numB[i]
				}
			}
		}
	}

	if (level == 0) {
		mu.gamma.0 <- data.frame(chain = numeric(nchains*numIntervals), Interval = character(nchains*numIntervals),
                                            value = numeric(nchains*numIntervals), stringsAsFactors = FALSE)
		mu.theta.0 <- data.frame(chain = numeric(nchains*numIntervals), Interval = character(nchains*numIntervals),
                                            value = numeric(nchains*numIntervals), stringsAsFactors = FALSE)
		tau2.gamma.0 <- data.frame(chain = numeric(nchains*numIntervals), Interval = character(nchains*numIntervals),
                                            value = numeric(nchains*numIntervals), stringsAsFactors = FALSE)
		tau2.theta.0 <- data.frame(chain = numeric(nchains*numIntervals), Interval = character(nchains*numIntervals),
                                            value = numeric(nchains*numIntervals), stringsAsFactors = FALSE)
		alpha.pi <- data.frame(chain = numeric(nchains*numIntervals), Interval = character(nchains*numIntervals),
                                            value = numeric(nchains*numIntervals), stringsAsFactors = FALSE)
		beta.pi <- data.frame(chain = numeric(nchains*numIntervals), Interval = character(nchains*numIntervals),
                                            value = numeric(nchains*numIntervals), stringsAsFactors = FALSE)

		for (i in 1:numIntervals) {
			mu.gamma.0[i,]$chain = 1
			mu.gamma.0[i,]$Interval = Intervals[i]
			mu.gamma.0[i,]$value = 0
			mu.theta.0[i,]$chain = 1
			mu.theta.0[i,]$Interval = Intervals[i]
			mu.theta.0[i,]$value = 0

			tau2.gamma.0[i,]$chain = 1
			tau2.gamma.0[i,]$Interval = Intervals[i]
			tau2.gamma.0[i,]$value = 0
			tau2.theta.0[i,]$chain = 1
			tau2.theta.0[i,]$Interval = Intervals[i]
			tau2.theta.0[i,]$value = 0

			alpha.pi[i,]$chain = 1
			alpha.pi[i,]$Interval = Intervals[i]
			alpha.pi[i,]$value = 1.5
			beta.pi[i,]$chain = 1
			beta.pi[i,]$Interval = Intervals[i]
			beta.pi[i,]$value = 1.5
		}
	
		if (nchains > 1) {
			for (c in 2:nchains) {
				for (i in 1:numIntervals) {
					mu.gamma.0[numIntervals * (c - 1) + i,]$chain = c
					mu.gamma.0[numIntervals * (c - 1) + i,]$Interval = Intervals[i]
					mu.gamma.0[numIntervals * (c - 1) + i,]$value = runif(1, -10, 10)
					mu.theta.0[numIntervals * (c - 1) + i,]$chain = c
					mu.theta.0[numIntervals * (c - 1) + i,]$Interval = Intervals[i]
					mu.theta.0[numIntervals * (c - 1) + i,]$value = runif(1, -10, 10)
	
					tau2.gamma.0[numIntervals * (c - 1) + i,]$chain = c
					tau2.gamma.0[numIntervals * (c - 1) + i,]$Interval = Intervals[i]
					tau2.gamma.0[numIntervals * (c - 1) + i,]$value = runif(1, 5, 20)
					tau2.theta.0[numIntervals * (c - 1) + i,]$chain = c
					tau2.theta.0[numIntervals * (c - 1) + i,]$Interval = Intervals[i]
					tau2.theta.0[numIntervals * (c - 1) + i,]$value = runif(1, 5, 20)
	
					alpha.pi[numIntervals * (c - 1) + i,]$chain = c
					alpha.pi[numIntervals * (c - 1) + i,]$Interval = Intervals[i]
					alpha.pi[numIntervals * (c - 1) + i,]$value = runif(1, 1.25, 100)
					beta.pi[numIntervals * (c - 1) + i,]$chain = c
					beta.pi[numIntervals * (c - 1) + i,]$Interval = Intervals[i]
					beta.pi[numIntervals * (c - 1) + i,]$value = runif(1, 1.25, 100)
				}
			}
		}
	}
	else {
		# level 1, 2
		mu.gamma.0 = rep(0, nchains)
		mu.theta.0 = rep(0, nchains)
		tau2.gamma.0 = rep(10, nchains)
		tau2.theta.0 = rep(10, nchains)

		alpha.pi = rep(1.5, nchains)
		beta.pi = rep(1.5, nchains)

		if (nchains > 1) {
			for (c in 2:nchains) {
				mu.gamma.0[c] = runif(1, -10, 10)
				mu.theta.0[c] = runif(1, -10, 10)
				tau2.gamma.0[c] = runif(1, 5, 20)
				tau2.theta.0[c] = runif(1, 5, 20)
				alpha.pi[c] = runif(1, 1.25, 100)
				beta.pi[c] = runif(1, 1.25, 100)
			}
		}
	}

	initial_values = list(gamma = gamma, theta = theta, mu.gamma = mu.gamma,
						mu.theta = mu.theta, sigma2.gamma = sigma2.gamma,
						sigma2.theta = sigma2.theta)

	if (model == "BB") {
		bb2_init = list(pi = pi)
		initial_values = c(initial_values, bb2_init)
	}

	if (hier == 3) {
		h3_init = list(mu.gamma.0 = mu.gamma.0, mu.theta.0 = mu.theta.0,
						tau2.gamma.0 = tau2.gamma.0,
                        tau2.theta.0 = tau2.theta.0)
		initial_values = c(initial_values, h3_init)
		if (model == "BB") {
			bb3_init = list(alpha.pi = alpha.pi, beta.pi = beta.pi)
			initial_values = c(initial_values, bb3_init)
		}
	}

	initial_values
}
