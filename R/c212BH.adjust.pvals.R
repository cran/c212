# c212.BH.adjust.pvals
# Case 2/12 - Adjust p values for use with the BH procedure
# R. Carragher
# Date: 25/11/2013

Id <- "$Id: c212BH.adjust.pvals.R,v 1.3 2016/10/14 10:39:03 clb13102 Exp clb13102 $"

c212.BH.adjust.pvals <- function(trial.data)
{
	# Check the correct columns are defined
	if (!("p" %in% colnames(trial.data))) {
		print("p column missing")
		return(NULL);
	}

	#trial.data <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
	ordered.data <- trial.data[order(trial.data$p),, drop = FALSE]

	m <- nrow(ordered.data)

	p_adj <- rep(NA, m)

	p_adj[m] = ordered.data$p[m]

	j <- m - 1
	while (j >= 1) {
		p_val <- ordered.data$p[j]

		p_cand <- (m/j)*p_val

		p_adj[j] = min(p_adj[j + 1], p_cand)

		j <- j - 1
	}

	ordered.data <- cbind(ordered.data, p_adj)

	row.names(ordered.data) = NULL

	return(ordered.data)
}
