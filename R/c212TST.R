# c212.TST
# Case 2/12 - Estimators of the proportions of true hypotheses
# R. Carragher
# Date: 29/11/2013

Id <- "$Id: c212TST.R,v 1.2 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

# Two stage method
c212.TST <- function(trial.data, alpha = 0.05) {

	# Check the correct columns are defined
	if (!("p" %in% colnames(trial.data))) {
		print("p column missing")
		return(NULL);
	}
	if (!("B" %in% colnames(trial.data))) {
		print("B column missing")
		return(NULL);
	}

	# No data in data frame - return NULL
	if (is.null(nrow(trial.data)) || (nrow(trial.data) == 0)) {
		print("No data found")
		return(NULL)
	}

	alpha_prime <- alpha/(1 + alpha)

	B <- unique(trial.data$B)
	numB <- length(B)

	# Apply the BH at level alpha_prime to the groups
	gamma <- rep(0, numB)

	for (i in 1:numB) {

		g <- trial.data[trial.data$B == B[i],, drop = FALSE]
		ng <- nrow(g)
		g_bh <- c212.BH(g, alpha_prime)

		rg = nrow(g_bh)

		gamma[i] <- (ng - rg)/ng
	}

	gamma_TST <- data.frame(B, pi0 = gamma)

	return(gamma_TST)
}
