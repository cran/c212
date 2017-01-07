# c212.DFDR
# Case 2/12 - Double False Discovery Rate
# R. Carragher
# Date: 25/11/2013

Id <- "$Id: c212DFDR.R,v 1.3 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.DFDR <- function(trial.data, alpha = 0.05)
{
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

	F <- trial.data[0,]

	# Body systems
	B <- unique(trial.data$B)
	numB <- length(B)

	# Find the representative p-value for the group
	# This is the minimum of the within group adjusted p-values for each group

	p_val <- rep(0, numB)

	for (i in 1:numB) {
		g <- trial.data[trial.data$B == B[i],, drop = FALSE]

		# Find the BH-adjusted p-values in the group
		p_grp <- c212.BH.adjust.pvals(g)

		p_val[i] <- min(p_grp$p_adj)
	}

	p_rep <- data.frame(B, p = p_val)

	p_rep_adj <- c212.BH.adjust.pvals(p_rep)

	# Create the family F

	for (i in 1:numB) {
		if (p_rep_adj$p_adj[i] <= alpha) {
			#print(sprintf("BS: %d added to F", B[i]))
			F <- rbind(F, trial.data[trial.data$B == p_rep_adj$B[i],, drop = FALSE])
		}
	}

	if (nrow(F) == 0) {
		return(F)
	}

	# Apply the BH procedure to F
	F_BH <- c212.BH(F, alpha)

	return(F_BH)
}
