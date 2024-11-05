# Global functions etc

M_global <- new.env()

M_global$logit <- function (p) {
	return (log(p/(1-p)))
}

M_global$checkCols <- function(cols, table) {

	for (i in 1:length(cols)) {

		if (!(cols[i] %in% colnames(table))) {
			print(sprintf("Missing %s column", cols[i]));
			return(1)
		}
	}
	0
}

M_global$checkNames <- function(n, table) {

	for (i in 1:length(n)) {

		if (!(n[i] %in% names(table))) {
			print(sprintf("Missing %s data", n[i]));
			return(1)
		}
	}
	0
}

# Code shared between c212.1a and c212.BB for end-of-trial data
M_global$EOTdata <- function(M_env, trial.data, iter, nchains, burnin, initial_values) {

	if (iter <= burnin || nchains < 1 || iter < 0 || burnin < 0) {
		print("Invalid simulation setup parametetrs");
		return(NULL)
	}

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors = FALSE)
	}

	if (is.null(trial.data)) {
		print("NULL trial data")
		return(NULL)
	}

	# Perform some validation checks
	if ((is.null(nrow(trial.data))) || (nrow(trial.data) == 0)) {
		print("Missing trial data set");
		return(NULL)
	}

	# Check the correct columns are defined
	cols = c("B", "AE", "Count", "Group", "Total")
	if (M_global$checkCols(cols, trial.data)) {
		print("Missing columns");
		return(NULL)
	}

	# Convert any factor columns to character (not needed if the data is in a file but if a data.frame
	# is passed to the function then we need to convert to strings)
	facs <- sapply(trial.data, is.factor)
	trial.data[facs] <- lapply(trial.data[facs], as.character)

    # Order by body-system, adverse event and group
	ordered.data <- trial.data[order(trial.data$B, trial.data$AE, trial.data$Group),,
																			drop=FALSE]

	# Assume 2 groups initially
	#num.groups <- max(trial.data$Group)

	cntrl.data <- ordered.data[ordered.data$Group == 1, ]
	treat.data <- ordered.data[ordered.data$Group == 2, ]

	# Check that we have matching body-systems and AEs in the control and treatment
	# groups
	# The data is ordered so a straight comparison is possible
	if(!identical(cntrl.data$B, treat.data$B)) {
		print("Mismatced body-system data");
		return(NULL)
	}
	if(!identical(cntrl.data$AE, treat.data$AE)) {
		print("Mismatced adverse event data");
		return(NULL)
	}

	M_env$initVars()

	# Size the data structures
	# BodySystems - unique preserves the order
	M_env$BodySys <- as.character(unique(cntrl.data$B))
	M_env$NumBodySys <- length(M_env$BodySys)

	# AEs
	M_env$nAE <- as.integer(by(cntrl.data, cntrl.data$B, nrow))
	maxAEs <- max(M_env$nAE)
	M_env$maxAEs <- maxAEs

	# Event Count Data
	M_env$x <- matrix(NA, nrow = M_env$NumBodySys, ncol = maxAEs)
	M_env$NC <- matrix(NA, nrow = M_env$NumBodySys, ncol = maxAEs)
	M_env$y <- matrix(NA, nrow = M_env$NumBodySys, ncol = maxAEs)
	M_env$NT <- matrix(NA, nrow = M_env$NumBodySys, ncol = maxAEs)
	M_env$AE = matrix(NA, nrow = M_env$NumBodySys, ncol = maxAEs)

	for (b in 1:M_env$NumBodySys) {
		bs <- M_env$BodySys[b]
		c <- cntrl.data[cntrl.data$B == bs, ]
		t <- treat.data[treat.data$B == bs, ]
		M_env$x[b,1:length(c$Count)] <- c$Count
		M_env$NC[b,1:length(c$Total)] <- c$Total
		M_env$y[b,1:length(t$Count)] <- t$Count
		M_env$NT[b,1:length(t$Total)] <- t$Total
		M_env$AE[b,1:length(c$AE)] <- as.character(c$AE)
	}

	M_env$initialiseChains(initial_values, nchains)

	out = list(trial.data = trial.data, cntrl.data = cntrl.data, treat.data = treat.data)

	return(out)
}

M_global$sim_params1a <- function(M_env, sim.params, sim_type, trial.data, cntrl.data) {

	# Have any of the global simulation parameters been overridden
	if (!is.null(sim.params)) {

		# 1. Check if we have a full set for the simulation type
		sim.params = sim.params[sim.params$type == sim_type, ]
	
		if (nrow(sim.params) > 0) {
			sim.params <- sim.params[order(sim.params$B, sim.params$AE),,drop=FALSE]
			gamma.params = sim.params[sim.params$variable == "gamma",]
			theta.params = sim.params[sim.params$variable == "theta",]

			if (!identical(cntrl.data$B, gamma.params$B) ||
					!identical(cntrl.data$AE, gamma.params$AE) ||
					!identical(cntrl.data$B, theta.params$B) ||
					!identical(cntrl.data$AE, theta.params$AE)) {
				# If we don't have a full set then we need to merge the values

				params = c212.sim.control.params(trial.data, "1a")
				params = params[params$type == sim_type, ]

				params = merge(params, sim.params, by = c("type", "variable", "B", "AE",
									"param"), all.x = T)

				params = params[, !(names(params) %in% c("value.x", "control.x"))]
				names(params)[names(params) == "value.y"] = "value"
				names(params)[names(params) == "control.y"] = "control"
				params <- params[order(params$B, params$AE),, drop=FALSE]

				sim.params = params
			}

			B = rep(1:M_env$NumBodySys)
			B = rep(B, M_env$nAE)
			j = sequence(M_env$nAE)

			# Only the indices are used in the C++ sampler
			B = as.integer(B)
			j = as.integer(j)

			gamma.params = sim.params[sim.params$variable == "gamma",]
			theta.params = sim.params[sim.params$variable == "theta",]

			gamma.params = gamma.params[, !(names(gamma.params) %in% c("B", "AE"))]
			gamma.params = cbind(gamma.params, B)
			gamma.params = cbind(gamma.params, j)

			theta.params = theta.params[, !(names(theta.params) %in% c("B", "AE"))]
			theta.params = cbind(theta.params, B)
			theta.params = cbind(theta.params, j)

			sim.params = rbind(gamma.params, theta.params)
			sim.params = sim.params[!is.na(sim.params$value), ]

			sim.params$value = as.numeric(sim.params$value)
			sim.params$control = as.numeric(sim.params$control)
		}
	}

	return(sim.params)
}


M_global$sim_paramsBB <- function(M_env, sim.params, pm.weights, sim_type, trial.data, cntrl.data) {

	# Have any global parameters been overridden?
	# 1. Simulation parameters
	# 2. Point mass weights

	# Only the indices are used in the C++ sampler
	B = rep(1:M_env$NumBodySys)
	B = rep(B, M_env$nAE)
	j = sequence(M_env$nAE)

	B = as.integer(B)
	j = as.integer(j)

	# Have any of the global simulation parameters been overridden
	if (!is.null(sim.params)) {

		# 1. Check if we have a full set for the simulation type
		sim.params <- sim.params[order(sim.params$B, sim.params$AE),,drop=FALSE]
		sim.params = sim.params[sim.params$type ==
						sim_type | sim.params$variable == "theta", ]

		if (nrow(sim.params) > 0) {

			gamma.params = sim.params[sim.params$variable == "gamma",]
			theta.params = sim.params[sim.params$variable == "theta",]

			if (!identical(cntrl.data$B, gamma.params$B) ||
				!identical(cntrl.data$AE, gamma.params$AE) ||
				!identical(cntrl.data$B, theta.params$B) ||
				!identical(cntrl.data$AE, theta.params$AE)) {
				# If we don't have a full set then we need to merge the values

				params = c212.sim.control.params(trial.data, "BB")
				params = params[params$type == sim_type | params$variable == "theta", ]

				params = merge(params, sim.params, by = c("type", "variable", "B", "AE",
																"param"), all.x = T)

				params = params[, !(names(params) %in% c("value.x", "control.x"))]
				names(params)[names(params) == "value.y"] = "value"
				names(params)[names(params) == "control.y"] = "control"
				params <- params[order(params$B, params$AE),, drop=FALSE]

				sim.params = params
			}

			gamma.params = sim.params[sim.params$variable == "gamma",]
			theta.params = sim.params[sim.params$variable == "theta",]

			gamma.params = gamma.params[, !(names(gamma.params) %in% c("B", "AE"))]
			gamma.params = cbind(gamma.params, B)
			gamma.params = cbind(gamma.params, j)

			theta.params = theta.params[, !(names(theta.params) %in% c("B", "AE"))]
			theta.params = cbind(theta.params, B)
			theta.params = cbind(theta.params, j)

			sim.params = rbind(gamma.params, theta.params)
			sim.params = sim.params[!is.na(sim.params$value), ]

			sim.params$value = as.numeric(sim.params$value)
			sim.params$control = as.numeric(sim.params$control)
		}
	}

	# Have any of the default point-mass weights been overridden
	if (!is.null(pm.weights)) {
		pm.weights <- pm.weights[order(pm.weights$B, pm.weights$AE),,drop=FALSE]

		if (!identical(cntrl.data$B, pm.weights$B) ||
				!identical(cntrl.data$AE, pm.weights$AE)) {
			# If we don't have a full set then we need to merge the values

			w = c212.pointmass.weights(trial.data)

			w = merge(w, pm.weights, by = c("B", "AE"), all.x = T)

			w = w[, !(names(w) %in% c("weight_pm.x"))]
			names(w)[names(w) == "weight_pm.y"] = "weight_pm"
			w <- w[order(w$B, w$AE),, drop=FALSE]

			pm.weights = w
		}

		pm.weights = pm.weights[, !(names(pm.weights) %in% c("B", "AE")), drop = FALSE]

		pm.weights = cbind(pm.weights, B)
		pm.weights = cbind(pm.weights, j)

		pm.weights = pm.weights[!is.na(pm.weights$weight_pm), ]
	}

	out = list(sim.params = sim.params, pm.weights = pm.weights)

	return(out)
}

M_global$INTERIMdata <- function(M_env, trial.data, iter, nchains, burnin, initial_values) {

	if (iter <= burnin || nchains < 1 || iter < 0 || burnin < 0) {
		print("Invalid simulation setup parametetrs");
		return(NULL)
	}

	M_env$initVars()

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors = FALSE)
	}

	# Perform some validation checks
	if (is.null(trial.data)) {
		print("NULL trial data")
		return(NULL)
	}

	if ((is.null(nrow(trial.data))) || (nrow(trial.data) == 0)) {
		print("Missing trial data set");
		return(NULL)
	}

	# Check the correct columns are defined
	cols = c("B", "AE", "Count", "Group", "Interval", "Exposure", "I_index")
	if (M_global$checkCols(cols, trial.data)) {
		print("Missing columns");
		return(NULL)
	}

	# Convert any factor columns to character (not needed if the data is in a file but if a data.frame
	# is passed to the function then we need to convert to strings)
	facs <- sapply(trial.data, is.factor)
	trial.data[facs] <- lapply(trial.data[facs], as.character)

	# Order by body-system, adverse event, interval and group
	ordered.data <- trial.data[order(trial.data$I_index, trial.data$B, trial.data$AE, trial.data$Group),, drop=FALSE]

	cntrl.data <- ordered.data[ordered.data$Group == 1, ]
	treat.data <- ordered.data[ordered.data$Group == 2, ]

	# Check that cntrl, treat data are matched 1 to 1
	# Check that we have matching body-systems, AEs and intervals in the control and
	# treatment groups
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

	# Size the data structures
	
	# At this stage the groups are matched exactly by Interval, BS, and AE
	# However it is possible that not all Body-systems or adverse events will appear
	# in all groups

	M_env$Intervals = unique(cntrl.data$Interval)
	M_env$numIntervals = length(M_env$Intervals)
	M_env$numB = c()
	for (i in 1:M_env$numIntervals) {
		c = cntrl.data[cntrl.data$Interval == M_env$Intervals[i],]
		b = unique(c$B)
		M_env$numB[i] = length(b)
	}
	M_env$maxBs = max(M_env$numB)
	M_env$B = array(NA, dim = c(M_env$numIntervals, M_env$maxBs))
	for (i in 1:M_env$numIntervals) {
		c = cntrl.data[cntrl.data$Interval == M_env$Intervals[i],]
		b = unique(c$B)
		M_env$B[i, 1:length(b)] = b
	}

	# AEs
	M_env$nAE = array(0, dim = c(M_env$numIntervals, M_env$maxBs))
	for (i in 1:M_env$numIntervals) {
		numB = M_env$numB[i]
		for (b in 1:numB) {
			M_env$nAE[i, b] = length(unique(cntrl.data[cntrl.data$B == M_env$B[i, b], ]$AE))
		}
	}

	maxAEs <- max(M_env$nAE)
	M_env$maxAEs <- maxAEs
	M_env$AE = array(NA, dim = c(M_env$numIntervals, M_env$maxBs, maxAEs))
	for (i in 1:M_env$numIntervals) {
		numB = M_env$numB[i]
		for (b in 1:numB) {
			c = cntrl.data[cntrl.data$Interval == M_env$Intervals[i] & cntrl.data$B == M_env$B[i, b],]
			# Do we need a unique here - each AE should only apper once in a body-system
			# in an interval - Let's add a check. Probably need to add it to c212.BB
			# etc. as well.
			if (length(unique(c$AE)) != length(c$AE)) {
				print("Multipe adverse events in body-system");
				return(NULL)
			}

			M_env$AE[i, b,1:length(c$AE)] = c$AE
		}
	}

	# Trial Data
	# Control: counts and exposure
	M_env$x = array(NA, dim = c(M_env$numIntervals, M_env$maxBs, maxAEs))
	M_env$C <- array(NA, dim = c(M_env$numIntervals, M_env$maxBs, maxAEs))

	# Treatment: counts and exposure
	M_env$y = array(NA, dim = c(M_env$numIntervals, M_env$maxBs, maxAEs))
	M_env$T <- array(NA, dim = c(M_env$numIntervals, M_env$maxBs, maxAEs))

	for (i in 1:M_env$numIntervals) {
		numB = M_env$numB[i]
		for (b in 1:numB) {
			control = cntrl.data[cntrl.data$Interval == M_env$Intervals[i] & cntrl.data$B == M_env$B[i, b], ]
			M_env$x[i, b, 1:length(control$Count)] <- control$Count
			M_env$C[i, b, 1:length(control$Count)] <- control$Exposure

			treatment = treat.data[treat.data$Interval == M_env$Intervals[i] & treat.data$B == M_env$B[i, b], ]
			M_env$y[i, b, 1:length(control$Count)] <- treatment$Count
			M_env$T[i, b, 1:length(control$Count)] <- treatment$Exposure
		}
	}

	M_env$initialiseChains(initial_values, nchains)

	out = list(trial.data = trial.data, cntrl.data = cntrl.data, treat.data = treat.data)

	return(out)
}

# Check that all intervals contain the same body-systems
# We only need to check the control data as, from above, the treatment data has the exact same structure
# This check is only used in level 1 models
M_global$checkBS <- function(M_env, cntrl.data) {

	c1 = cntrl.data[cntrl.data$Interval == M_env$Intervals[1],]
	b1 = c1$B
	if (M_env$numIntervals > 1) {
		for (i in 2:M_env$numIntervals) {
			c = cntrl.data[cntrl.data$Interval == M_env$Intervals[i],]
			b = c$B

			if (!identical(b1, b)) {
				print("Intervals contain different body-systems");
				return(1)
			}
		}
	}
	0
}

M_global$INTERIM_sim_params1a <- function(M_env, sim.params, sim_type, trial.data, cntrl.data) {

	# Have any of the global simulation parameters been overridden
	if (!is.null(sim.params)) {
		# 1. Check if we have a full set for the simulation type
		sim.params = sim.params[sim.params$type == sim_type, ]
		sim.params <- sim.params[order(sim.params$I_index, sim.params$B, sim.params$AE),,drop=FALSE]
		gamma.params = sim.params[sim.params$variable == "gamma",]
		theta.params = sim.params[sim.params$variable == "theta",]

		if (!identical(cntrl.data$B, gamma.params$B) ||
			!identical(cntrl.data$AE, gamma.params$AE) ||
			!identical(cntrl.data$Interval, gamma.params$Interval) ||
			!identical(cntrl.data$B, theta.params$B) ||
			!identical(cntrl.data$AE, theta.params$AE) || 
			!identical(cntrl.data$Interval, theta.params$Interval)) {

			# If we don't have a full set then we need to merge the values
			params = c212.sim.control.params(trial.data, "1a")
			params = params[params$type == sim_type, ]

			params = merge(params, sim.params, by = c("type", "variable", "Interval", "I_index", "B", "AE", "param"), all.x = T)

			params = params[, !(names(params) %in% c("value.x", "control.x"))]
			names(params)[names(params) == "value.y"] = "value"
			names(params)[names(params) == "control.y"] = "control"
			params <- params[order(params$I_index, params$B, params$AE),, drop=FALSE]

			sim.params = params
		}

		# Add in the indices for intervals, body-systems and aes
		I_index = c()
		B = c()
		j = c()

		x = apply(M_env$nAE, 1, sum)
		I_index = rep(1:M_env$numIntervals, x)
		d = cbind(M_env$nAE, M_env$numB)
		B = as.vector(apply(d, 1, function(x) { n = length(x); rep(1:x[n], x[1:x[n]])})) # Should handle missing BS's - maybe
		j = as.vector(apply(M_env$nAE,1, sequence))

		I_index = as.integer(I_index)
		B = as.integer(B)
		j = as.integer(j)

		gamma.params = sim.params[sim.params$variable == "gamma",]
		theta.params = sim.params[sim.params$variable == "theta",]

		gamma.params = gamma.params[, !(names(gamma.params) %in% c("Interval", "I_index","B", "AE"))]
		gamma.params = cbind(gamma.params, I_index)
		gamma.params = cbind(gamma.params, B)
		gamma.params = cbind(gamma.params, j)

		theta.params = theta.params[, !(names(theta.params) %in% c("Interval", "I_index","B", "AE"))]
		theta.params = cbind(theta.params, I_index)
		theta.params = cbind(theta.params, B)
		theta.params = cbind(theta.params, j)

		sim.params = rbind(gamma.params, theta.params)
		sim.params = sim.params[!is.na(sim.params$value), ]

		sim.params$value = as.numeric(sim.params$value)
		sim.params$control = as.numeric(sim.params$control)
	}

	return(sim.params)
}

M_global$INTERIM_sim_paramsBB_2 <- function(M_env, sim.params, pm.weights, sim_type, trial.data, cntrl.data) {

	# Have any of the global parameters been overridden?
	# 1. Simulation parameters
	# 2. Point mass weights

	# Indices for changed parameters
	I_index = c()
	B = c()
	j = c()

	x = apply(M_env$nAE, 1, sum)
	I_index = rep(1:M_env$numIntervals, x)
	d = cbind(M_env$nAE, M_env$numB)
	B = as.vector(apply(d, 1, function(x) { n = length(x); rep(1:x[n], x[1:x[n]])})) # Should handle missing BS's - maybe
	j = as.vector(apply(M_env$nAE,1, sequence))

	I_index = as.integer(I_index)
	B = as.integer(B)
	j = as.integer(j)

	if ((!is.null(sim.params)) && (nrow(sim.params) > 0)) {

		# 1. Check if we have a full set for the simulation type
		sim.params <- sim.params[order(sim.params$I_index, sim.params$B, sim.params$AE),,drop=FALSE]
		sim.params = sim.params[sim.params$type == sim_type | sim.params$variable == "theta", ]

		if (nrow(sim.params) > 0) {

			gamma.params = sim.params[sim.params$variable == "gamma",]
			theta.params = sim.params[sim.params$variable == "theta",]

			if (!identical(cntrl.data$B, gamma.params$B) ||
				!identical(cntrl.data$AE, gamma.params$AE) ||
				!identical(cntrl.data$Interval, gamma.params$Interval) ||
				!identical(cntrl.data$B, theta.params$B) ||
				!identical(cntrl.data$AE, theta.params$AE) || 
				!identical(cntrl.data$Interval, theta.params$Interval)) {

				# If we don't have a full set then we need to merge the values
				params = c212.sim.control.params(trial.data, "BB")
				params = params[params$variable == "theta" | params$variable == "gamma",]

				params = params[params$type == sim_type | params$variable == "theta", ]

				params = merge(params, sim.params, by = c("type", "variable", "Interval", "I_index", "B", "AE", "param"), all.x = T)

				params = params[, !(names(params) %in% c("value.x", "control.x"))]
				names(params)[names(params) == "value.y"] = "value"
				names(params)[names(params) == "control.y"] = "control"
				params <- params[order(params$I_index, params$B, params$AE),, drop=FALSE]
	
				sim.params = params
			}

			# Add in the indices for intervals, body-systems and aes

			gamma.params = sim.params[sim.params$variable == "gamma",]
			theta.params = sim.params[sim.params$variable == "theta",]

			gamma.params = gamma.params[, !(names(gamma.params) %in% c("Interval", "I_index","B", "AE"))]
			gamma.params = cbind(gamma.params, I_index)
			gamma.params = cbind(gamma.params, B)
			gamma.params = cbind(gamma.params, j)

			theta.params = theta.params[, !(names(theta.params) %in% c("Interval", "I_index","B", "AE"))]
			theta.params = cbind(theta.params, I_index)
			theta.params = cbind(theta.params, B)
			theta.params = cbind(theta.params, j)

			sim.params = rbind(gamma.params, theta.params)
			sim.params = sim.params[!is.na(sim.params$value), ]

			sim.params$value = as.numeric(sim.params$value)
			sim.params$control = as.numeric(sim.params$control)
		}
	}

	pm.weights <- M_global$INTERIM_pm_weights(M_env, pm.weights, trial.data, cntrl.data, I_index, B, j)

	out = list(sim.params = sim.params, pm.weights = pm.weights)

	return(out)
}

M_global$INTERIM_sim_paramsBB_3 <- function(M_env, sim.params, pm.weights, sim_type, trial.data, cntrl.data) {

	# Have any of the global parameters been overridden?
	# 1. Simulation parameters
	# 2. Point mass weights

	# Indices for changed parameters
	I_index = c()
	B = c()
	j = c()

	x = apply(M_env$nAE, 1, sum)
	I_index = rep(1:M_env$numIntervals, x)
	d = cbind(M_env$nAE, M_env$numB)
	B = as.vector(apply(d, 1, function(x) { n = length(x); rep(1:x[n], x[1:x[n]])})) # Should handle missing BS's - maybe
	j = as.vector(apply(M_env$nAE,1, sequence))

	I_index = as.integer(I_index)
	B = as.integer(B)
	j = as.integer(j)

	if ((!is.null(sim.params)) && (nrow(sim.params) > 0)) {

		# 1. Check if we have a full set for the simulation type
		sim.params <- sim.params[order(sim.params$I_index, sim.params$B, sim.params$AE), ,drop=FALSE]
		sim.params = sim.params[sim.params$type == sim_type | sim.params$variable == "theta", ]

		if (nrow(sim.params) > 0) {

			gamma.params = sim.params[sim.params$variable == "gamma",]
			theta.params = sim.params[sim.params$variable == "theta",]
			alpha.params = sim.params[sim.params$variable == "alpha",]
			beta.params = sim.params[sim.params$variable == "beta",]

			# If used at all, alpha/beta apply to Intervals only
			if (nrow(alpha.params) > 0) {
				alpha.params$B <- NA
				alpha.params$AE <- NA
			}
			if (nrow(beta.params) > 0) {
				beta.params$B <- NA
				beta.params$AE <- NA
			}

			if (!identical(cntrl.data$B, gamma.params$B) ||
				!identical(cntrl.data$AE, gamma.params$AE) ||
				!identical(cntrl.data$Interval, gamma.params$Interval) ||
				!identical(cntrl.data$B, theta.params$B) ||
				!identical(cntrl.data$AE, theta.params$AE) || 
				!identical(cntrl.data$Interval, theta.params$Interval) ||
				!identical(unique(cntrl.data$Interval), alpha.params$Interval) ||
				!identical(unique(cntrl.data$Interval), beta.params$Interval)) {

				# If we don't have a full set then we need to merge the values
				params = c212.sim.control.params(trial.data, "BB")

				params = params[params$type == sim_type | params$variable == "theta", ]

				params = merge(params, sim.params, by = c("type", "variable", "Interval", "I_index", "B", "AE", "param"), all.x = T)

				params = params[, !(names(params) %in% c("value.x", "control.x"))]
				names(params)[names(params) == "value.y"] = "value"
				names(params)[names(params) == "control.y"] = "control"
				params <- params[order(params$I_index, params$B, params$AE),, drop=FALSE]
	
				sim.params = params
			}

			# Add in the indices for intervals, body-systems and aes

			gamma.params = sim.params[sim.params$variable == "gamma",]
			theta.params = sim.params[sim.params$variable == "theta",]
			alpha.params = sim.params[sim.params$variable == "alpha",]
			beta.params = sim.params[sim.params$variable == "beta",]

			gamma.params = gamma.params[, !(names(gamma.params) %in% c("Interval", "I_index","B", "AE"))]
			gamma.params = cbind(gamma.params, I_index)
			gamma.params = cbind(gamma.params, B)
			gamma.params = cbind(gamma.params, j)

			theta.params = theta.params[, !(names(theta.params) %in% c("Interval", "I_index","B", "AE"))]
			theta.params = cbind(theta.params, I_index)
			theta.params = cbind(theta.params, B)
			theta.params = cbind(theta.params, j)

			alpha.params = alpha.params[, !(names(alpha.params) %in% c("Interval", "I_index","B", "AE"))]
			alpha.params = cbind(alpha.params, I_index = unique(I_index))
			alpha.params$B <- NA
			alpha.params$j <- NA

			beta.params = beta.params[, !(names(beta.params) %in% c("Interval", "I_index","B", "AE"))]
			beta.params = cbind(beta.params, I_index = unique(I_index))
			beta.params$B <- NA
			beta.params$j <- NA

			sim.params = rbind(gamma.params, theta.params, alpha.params, beta.params)
			sim.params = sim.params[!is.na(sim.params$value), ]

			sim.params$value = as.numeric(sim.params$value)
			sim.params$control = as.numeric(sim.params$control)
		}
	}

	pm.weights <- M_global$INTERIM_pm_weights(M_env, pm.weights, trial.data, cntrl.data, I_index, B, j)

	out = list(sim.params = sim.params, pm.weights = pm.weights)

	return(out)
}

M_global$INTERIM_pm_weights <- function(M_env, pm.weights, trial.data, cntrl.data, I_index, B, j) {

	# Have any of the default point-mass weights been overridden
	if ((!is.null(pm.weights)) && (nrow(pm.weights) > 0)) {

		pm.weights <- pm.weights[order(pm.weights$I_index, pm.weights$B, pm.weights$AE),,drop=FALSE]

		if (!identical(cntrl.data$B, pm.weights$B) ||
			!identical(cntrl.data$AE, pm.weights$AE)) {
			# If we don't have a full set then we need to merge the values

			w = c212.pointmass.weights(trial.data)

			w = merge(w, pm.weights, by = c("Interval", "I_index", "B", "AE"), all.x = T)

			w = w[, !(names(w) %in% c("weight_pm.x"))]
			names(w)[names(w) == "weight_pm.y"] = "weight_pm"
			w <- w[order(w$I_index, w$B, w$AE),, drop=FALSE]

			pm.weights = w
		}

		pm.weights = pm.weights[, !(names(pm.weights) %in% c("Interval", "I_index", "B", "AE")), drop = FALSE]

		pm.weights = cbind(pm.weights, I_index)
		pm.weights = cbind(pm.weights, B)
		pm.weights = cbind(pm.weights, j)

		pm.weights = pm.weights[!is.na(pm.weights$weight_pm), ]
	}

	return(pm.weights)
}

M_global$INTERIM_monitor_1a_2 <- function(monitor) {

	# monitor null means monitor all variables
	if (is.null(monitor)) {
		v = c("theta", "gamma", "mu.gamma", "mu.theta", "sigma2.theta", "sigma2.gamma")
		s = rep(1, length(v))
		monitor = data.frame(variable = v, monitor = s, stringsAsFactors = FALSE)
	}
	else {
		# Does monitor include values for all possible variables?
		if (nrow(monitor[monitor$variable == "theta", ]) == 0) {
			monitor = rbind(monitor, c("theta", 0))
		}
		if (nrow(monitor[monitor$variable == "gamma", ]) == 0) {
			monitor = rbind(monitor, c("gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.theta", ]) == 0) {
			monitor = rbind(monitor, c("mu.theta", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.gamma", ]) == 0) {
			monitor = rbind(monitor, c("mu.gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "sigma2.theta", ]) == 0) {
			monitor = rbind(monitor, c("sigma2.theta", 0))
		}
		if (nrow(monitor[monitor$variable == "sigma2.gamma", ]) == 0) {
			monitor = rbind(monitor, c("sigma2.gamma", 0))
		}
	}
	# Coerce any non-correct type
	monitor$monitor = as.integer(monitor$monitor)

	return(monitor)
}

M_global$INTERIM_monitor_1a_3 <- function(monitor) {

	# monitor null means monitor all variables
	if (is.null(monitor)) {
		v = c("theta", "gamma", "mu.gamma", "mu.theta", "sigma2.theta", "sigma2.gamma",
				"mu.theta.0", "mu.gamma.0", "tau2.theta.0", "tau2.gamma.0")
		s = rep(1, length(v))
		monitor = data.frame(variable = v, monitor = s, stringsAsFactors = FALSE)
	}
	else {
		# Does monitor include values for all possible variables?
		if (nrow(monitor[monitor$variable == "theta", ]) == 0) {
			monitor = rbind(monitor, c("theta", 0))
		}
		if (nrow(monitor[monitor$variable == "gamma", ]) == 0) {
			monitor = rbind(monitor, c("gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.theta", ]) == 0) {
			monitor = rbind(monitor, c("mu.theta", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.gamma", ]) == 0) {
			monitor = rbind(monitor, c("mu.gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "sigma2.theta", ]) == 0) {
			monitor = rbind(monitor, c("sigma2.theta", 0))
		}
		if (nrow(monitor[monitor$variable == "sigma2.gamma", ]) == 0) {
			monitor = rbind(monitor, c("sigma2.gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.theta.0", ]) == 0) {
			monitor = rbind(monitor, c("mu.theta.0", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.gamma.0", ]) == 0) {
			monitor = rbind(monitor, c("mu.gamma.0", 0))
		}
		if (nrow(monitor[monitor$variable == "tau2.gamma.0", ]) == 0) {
			monitor = rbind(monitor, c("tau2.gamma.0", 0))
		}
		if (nrow(monitor[monitor$variable == "tau2.theta.0", ]) == 0) {
			monitor = rbind(monitor, c("tau2.theta.0", 0))
		}
	}
	# Coerce any non-correct type
	monitor$monitor = as.integer(monitor$monitor)

	return(monitor)
}

M_global$INTERIM_monitor_BB_2 <- function(monitor) {

	# monitor null means monitor all variables
	if (is.null(monitor)) {
		v = c("theta", "gamma", "mu.gamma", "mu.theta", "sigma2.theta", "sigma2.gamma",
				"pi")
		s = rep(1, length(v))
		monitor = data.frame(variable = v, monitor = s, stringsAsFactors = FALSE)
	}
	else {
		# Does monitor include values for all possible variables?
		if (nrow(monitor[monitor$variable == "theta", ]) == 0) {
			monitor = rbind(monitor, c("theta", 0))
		}
		if (nrow(monitor[monitor$variable == "gamma", ]) == 0) {
			monitor = rbind(monitor, c("gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.theta", ]) == 0) {
			monitor = rbind(monitor, c("mu.theta", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.gamma", ]) == 0) {
			monitor = rbind(monitor, c("mu.gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "sigma2.theta", ]) == 0) {
			monitor = rbind(monitor, c("sigma2.theta", 0))
		}
		if (nrow(monitor[monitor$variable == "sigma2.gamma", ]) == 0) {
			monitor = rbind(monitor, c("sigma2.gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "pi", ]) == 0) {
			monitor = rbind(monitor, c("pi", 0))
		}
	}
	# Coerce any non-correct type
	monitor$monitor = as.integer(monitor$monitor)

	return(monitor)
}

M_global$INTERIM_monitor_BB_3 <- function(monitor) {

	# monitor null means monitor all variables
	if (is.null(monitor)) {
		v = c("theta", "gamma", "mu.gamma", "mu.theta", "sigma2.theta", "sigma2.gamma",
				"mu.theta.0", "mu.gamma.0", "tau2.theta.0", "tau2.gamma.0",
				"pi", "alpha.pi", "beta.pi")
		s = rep(1, length(v))
		monitor = data.frame(variable = v, monitor = s, stringsAsFactors = FALSE)
	}
	else {
		# Does monitor include values for all possible variables?
		if (nrow(monitor[monitor$variable == "theta", ]) == 0) {
			monitor = rbind(monitor, c("theta", 0))
		}
		if (nrow(monitor[monitor$variable == "gamma", ]) == 0) {
			monitor = rbind(monitor, c("gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.theta", ]) == 0) {
			monitor = rbind(monitor, c("mu.theta", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.gamma", ]) == 0) {
			monitor = rbind(monitor, c("mu.gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "sigma2.theta", ]) == 0) {
			monitor = rbind(monitor, c("sigma2.theta", 0))
		}
		if (nrow(monitor[monitor$variable == "sigma2.gamma", ]) == 0) {
			monitor = rbind(monitor, c("sigma2.gamma", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.theta.0", ]) == 0) {
			monitor = rbind(monitor, c("mu.theta.0", 0))
		}
		if (nrow(monitor[monitor$variable == "mu.gamma.0", ]) == 0) {
			monitor = rbind(monitor, c("mu.gamma.0", 0))
		}
		if (nrow(monitor[monitor$variable == "tau2.gamma.0", ]) == 0) {
			monitor = rbind(monitor, c("tau2.gamma.0", 0))
		}
		if (nrow(monitor[monitor$variable == "tau2.theta.0", ]) == 0) {
			monitor = rbind(monitor, c("tau2.theta.0", 0))
		}
		if (nrow(monitor[monitor$variable == "pi", ]) == 0) {
			monitor = rbind(monitor, c("pi", 0))
		}
		if (nrow(monitor[monitor$variable == "alpha.pi", ]) == 0) {
			monitor = rbind(monitor, c("alpha.pi", 0))
		}
		if (nrow(monitor[monitor$variable == "beta.pi", ]) == 0) {
			monitor = rbind(monitor, c("beta.pi", 0))
		}
	}
	# Coerce any non-correct type
	monitor$monitor = as.integer(monitor$monitor)

	return(monitor)
}

M_global$INTERIM_check_conv_name_1a_2 <- function(raw) {

	# Check which variables we are monitoring
	monitor = raw$monitor
	theta_mon = monitor[monitor$variable == "theta",]$monitor
	gamma_mon = monitor[monitor$variable == "gamma",]$monitor
	mu.theta_mon = monitor[monitor$variable == "mu.theta",]$monitor
	mu.gamma_mon = monitor[monitor$variable == "mu.gamma",]$monitor
	sigma2.theta_mon = monitor[monitor$variable == "sigma2.theta",]$monitor
	sigma2.gamma_mon = monitor[monitor$variable == "sigma2.gamma",]$monitor

	n = c("chains", "nIntervals", "Intervals", "nBodySys", "maxBs", "maxAEs", "nAE", "B", "AE", "iter")

	if (theta_mon == 1) {
		n = c(n, "theta")
		if (raw$sim_type == "MH") {
			n = c(n, "theta_acc")
		}
	}
	if (gamma_mon == 1 ) {
		n = c(n, "gamma")
		if (raw$sim_type == "MH") {
			n = c(n, "gamma_acc")
		}
	}
	if (mu.gamma_mon == 1) {
		n = c(n, "mu.gamma")
	}
	if (mu.theta_mon == 1) {
		n = c(n, "mu.theta")
	}
	if (sigma2.theta_mon) {
		n = c(n, "sigma2.theta")
	}
	if (sigma2.gamma_mon) {
		n = c(n, "sigma2.gamma")
	}

	if (M_global$checkNames(n, raw)) {
		print("Missing names");
		return(1)
	}
	0
}

M_global$INTERIM_check_conv_name_1a_3 <- function(raw) {

	if (M_global$INTERIM_check_conv_name_1a_2(raw)) {
		return(1)
	}

	monitor = raw$monitor
	mu.theta.0_mon = monitor[monitor$variable == "mu.theta.0",]$monitor
	mu.gamma.0_mon = monitor[monitor$variable == "mu.gamma.0",]$monitor
	tau2.theta.0_mon = monitor[monitor$variable == "tau2.theta.0",]$monitor
	tau2.gamma.0_mon = monitor[monitor$variable == "tau2.gamma.0",]$monitor

	n = c()

    if (mu.gamma.0_mon == 1) {
		n = c(n, "mu.gamma.0")
    }
    if (mu.theta.0_mon == 1) {
		n = c(n, "mu.theta.0")
    }
    if (tau2.gamma.0_mon == 1) {
		n = c(n, "tau2.gamma.0")
    }
    if (tau2.theta.0_mon == 1) {
		n = c(n, "tau2.theta.0")
    }

	if (length(n) > 0) {
		if (M_global$checkNames(n, raw)) {
			return(1)
		}
    }

	0
}

M_global$INTERIM_check_summ_name_1a_2 <- function(raw) {

	# Check which variables we are monitoring
	monitor = raw$monitor
	theta_mon = monitor[monitor$variable == "theta",]$monitor
	gamma_mon = monitor[monitor$variable == "gamma",]$monitor
	mu.theta_mon = monitor[monitor$variable == "mu.theta",]$monitor
	mu.gamma_mon = monitor[monitor$variable == "mu.gamma",]$monitor
	sigma2.theta_mon = monitor[monitor$variable == "sigma2.theta",]$monitor
	sigma2.gamma_mon = monitor[monitor$variable == "sigma2.gamma",]$monitor

	n = c("chains", "nIntervals", "nBodySys", "maxBs", "maxAEs", "nAE", "B", "AE", "iter", "burnin")

	if (theta_mon == 1) {
		 n = c(n, "theta")
	}
	if (gamma_mon == 1) {
		 n = c(n, "gamma")
	}
	if (mu.gamma_mon == 1) {
		 n = c(n, "mu.gamma")
	}
	if (theta_mon == 1) {
		 n = c(n, "theta")
	}
	if (mu.theta_mon == 1) {
		 n = c(n, "mu.theta")
	}
	if (sigma2.theta_mon == 1) {
		 n = c(n, "sigma2.theta")
	}
	if (sigma2.gamma_mon == 1) {
		 n = c(n, "sigma2.gamma")
	}

	if (length(n) > 0) {
		if (M_global$checkNames(n, raw)) {
			return(1)
		}
	}

    0
}

M_global$INTERIM_check_summ_name_1a_3 <- function(raw) {

	if (M_global$INTERIM_check_summ_name_1a_2(raw)) {
		return(1)
	}

	monitor = raw$monitor
	mu.theta.0_mon = monitor[monitor$variable == "mu.theta.0",]$monitor
	mu.gamma.0_mon = monitor[monitor$variable == "mu.gamma.0",]$monitor
	tau2.theta.0_mon = monitor[monitor$variable == "tau2.theta.0",]$monitor
	tau2.gamma.0_mon = monitor[monitor$variable == "tau2.gamma.0",]$monitor

	n = c()

    if (mu.gamma.0_mon == 1) {
		n = c(n, "mu.gamma.0")
    }
    if (mu.theta.0_mon == 1) {
		n = c(n, "mu.theta.0")
    }
    if (tau2.gamma.0_mon == 1) {
		n = c(n, "tau2.gamma.0")
    }
    if (tau2.theta.0_mon == 1) {
		n = c(n, "tau2.theta.0")
    }

	if (length(n) > 0) {
		if (M_global$checkNames(n, raw)) {
			return(1)
		}
    }

	0
}

M_global$Geweke <- function(x) {
	mcmc_obj <- mcmc(x)
	g <- geweke.diag(mcmc_obj)

	return(g)
}

M_global$GelmanRubin <- function(x, nchains) {

	mcmc_obj <- list(NA)

	x1 = split(x, row(x))
	x2 = lapply(x1, mcmc)
	mlist <- mcmc.list(x2)

	g <- gelman.diag(mlist)

	return(g)
}	

M_global$GelmanRubin_new <- function(x, nchains) {

	mcmc_obj <- list(NA)

	x1 = split(x, row(x))
	x2 = lapply(x1, mcmc)
	mlist <- mcmc.list(x2)

	g <- gelman.diag(mlist)

	return(c(g$psrf[1], g$psrf[2]))
}	

M_global$summaryStats <- function(x, nchains) {

	if (nchains == 1) {
		x = as.matrix(t(x))
	}

	x1 = c(x[1:nchains,])
	m <- mcmc(x1)
	h <- HPDinterval(m)
	m = c(mean(x1), median(x1))

	mcmc_obj <- list(NA)
	x1 = split(x, row(x))
	x2 = lapply(x1, mcmc)
	mlist <- mcmc.list(x2)
	stats = summary(mlist)

	return(c(m[1], m[2], h[1], h[2], stats$statistics["SD"], stats$statistics["Time-series SE"]))
}

chk_val <- function(val, q = 0.975) {
	if (abs(val) > qnorm(q)) {
		return("*")
	}
	else {
		return("-")
	}
}


M_global$EOTprintConvSummLev1 <- function(x, text, chk = 0) {
	max_t = head(x[x$stat == max(x$stat),, drop = FALSE], 1)
	min_t = head(x[x$stat == min(x$stat),, drop = FALSE], 1)

	if (chk == 0) {
		cat(sprintf("Max %s (%s %s): %0.6f\n", text, max_t$B, max_t$AE, max_t$stat))
		cat(sprintf("Min %s (%s %s): %0.6f\n", text, min_t$B, min_t$AE, min_t$stat))
	}
	else {
		cat(sprintf("Max %s (%s %s): %0.6f (%s)\n", text, max_t$B, max_t$AE, max_t$stat, chk_val(max_t$stat)))
		cat(sprintf("Min %s (%s %s): %0.6f (%s)\n", text, min_t$B, min_t$AE, min_t$stat, chk_val(min_t$stat)))
	}
}

M_global$EOTprintConvSummLev2 <- function(x, text, chk = 0) {
	max_t = head(x[x$stat == max(x$stat),, drop = FALSE], 1)
	min_t = head(x[x$stat == min(x$stat),, drop = FALSE], 1)

	if (chk == 0) {
		cat(sprintf("Max %s (%s): %0.6f\n", text, max_t$B, max_t$stat))
		cat(sprintf("Min %s (%s): %0.6f\n", text, min_t$B, min_t$stat))
	}
	else {
		cat(sprintf("Max %s (%s): %0.6f (%s)\n", text, max_t$B, max_t$stat, chk_val(max_t$stat)))
		cat(sprintf("Min %s (%s): %0.6f (%s)\n", text, min_t$B, min_t$stat, chk_val(min_t$stat)))
	}
}

M_global$EOTprintConvSummLev3 <- function(x, text, chk = 0) {

	if (chk == 0) {
		cat(sprintf("%s: %0.6f\n", text, x$stat))
	}
	else {
		cat(sprintf("%s: %0.6f (%s)\n", text, x$stat, chk_val(x$stat)))
	}
}
