# c212.BONF
# Case 2/12 - Bonferroni corretion
# R. Carragher
# Date: 25/11/2013

Id <- "$Id: c212BONF.R,v 1.3 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.BONF <- function(trial.data, alpha = 0.05)
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

	n <- nrow(trial.data)

	if (is.null(n)) {
		print("No data found");
		return(NULL);
	}

	# Empty data set - just return it
	if (n == 0) {
		return (trial.data)
	}

	sig <- trial.data[trial.data$p <= alpha/n,, drop = FALSE ]

	sig <- sig[order(sig$p),, drop=FALSE]

	row.names(sig) = NULL

	return(sig);
}
