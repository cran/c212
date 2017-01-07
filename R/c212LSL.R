# c212.LSL
# Case 2/12 - Estimators of the proportions of true hypotheses
# R. Carragher
# Date: 29/11/2013


Id <- "$Id: c212LSL.R,v 1.2 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

# Least slope method
c212.LSL <- function(trial.data) {

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

	B <- unique(trial.data$B)
	numB <- length(B)

	gamma <- rep(0, numB)

	for (i in 1:numB) {
		g <- trial.data[trial.data$B == B[i],, drop = FALSE]
		g_ordered <- g[order(g$p),, drop=FALSE]
		ng <- nrow(g)

		j = 1
		l = (ng:1)/(1 - g_ordered$p)

		if (ng > 1) {
			for (j in 2:ng) {
				if (l[j] > l[j - 1]) {
					break;
				}
			}
		}

		gamma[i] <- min((floor(l[j]) + 1)/ng, 1)
	}

	gamma_LSL <- data.frame(B, pi0 = gamma)

	return(gamma_LSL)
}
