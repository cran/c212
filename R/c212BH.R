# c212.BH
# Case 2/12 - Implementation of BH procedure for controlling the FDR
# R. Carragher
# Date: 25/11/2013

Id <- "$Id: c212BH.R,v 1.3 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.BH <- function(trial.data, alpha = 0.05) {

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

	m <- nrow(trial.data)

	if (is.null(m)) {
		print("No data found");
		return(NULL);
	}

	# No data in data frame - just return it
	if (m == 0) {
		return(trial.data)
	}

	ordered.data <- trial.data[order(trial.data$p),, drop=FALSE]

	k <- 0

	for (j in 1:m) {
		if (ordered.data$p[j] <= j*alpha/m) {
			k <- j
		}
	}

	if (k > 0) {
		sig <- ordered.data[1:k,, drop=FALSE]
	}
	else {
		sig <- ordered.data[0,, drop=FALSE]
	}

	row.names(sig) = NULL
	return(sig)
}
