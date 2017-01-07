# c212.GBH
# Case 2/12 - Group BH procedure
# R. Carragher
# Date: 25/11/2013

Id <- "$Id: c212GBH.R,v 1.3 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.GBH <- function(trial.data, pi0 = "TST", alpha = 0.05) {

	if (is.null(trial.data)) {
		print("NULL trial data")
		return(NULL)
	}

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
	}

	# Check the correct columns are defined
	if (!("p" %in% colnames(trial.data))) {
		print("p column missing")
		return(NULL);
	}
	if (!("B" %in% colnames(trial.data))) {
		print("B column missing")
		return(NULL);
	}

    # No data in data frame - just return it
    if (is.null(nrow(trial.data)) || (nrow(trial.data) == 0)) {
        return(trial.data)
    }

	if (pi0 == "TST") {
		pi0_est <- c212.TST(trial.data, alpha)
	}
	else {
		pi0_est <- c212.LSL(trial.data)
	}

	if (is.null(pi0_est)) {
		return(NULL)
	}

	# If all estimated proportions of true hypotheses are 1, return no significant hypotheses
	t <- rle(pi0_est$pi)$values
	if (length(t) == 1 && t == 1) {
		return(trial.data[0,, drop = FALSE])
	}


	sig <- GBH(trial.data, pi0_est, alpha)

	row.names(sig) = NULL

	return(sig)
}

GBH <- function(trial.data, pi0_est, alpha) {

	B <- unique(trial.data$B)
	numB <- length(B)

	m <- nrow(trial.data)

	p_weighted <- rep(0, m)

	for (i in 1:m) {
		b <- trial.data$B[i]
		pi0 <- pi0_est[B == b, ]$pi0
		p_weighted[i] <- trial.data$p[i] * ((pi0)/(1 - pi0))
	}
	trial.data <- cbind(trial.data, p_weighted)

	ordered.data <- trial.data[order(trial.data$p_weighted),, drop=FALSE]

	pi0 <- 0
	for (i in 1:numB) {
		ng <- nrow(trial.data[trial.data$B == B[i],, drop = FALSE])
		pi0 <- pi0 + ng * pi0_est[B == B[i], ]$pi0
	}
	pi0 <- pi0/m

	alpha_w <- alpha/(1 - pi0)

	k <- 0

	for (i in 1:m) {
		if (ordered.data$p_weighted[i] <= i * alpha_w/m) {
			k <- i
		}
	}

	if (k > 0) {
		sig <- ordered.data[1:k,, drop = FALSE]
	}
	else {
		sig <- ordered.data[0,, drop = FALSE]
	}

	return(sig)
}
