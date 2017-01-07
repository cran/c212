# c212.ssBH
# Case 2/12 subset BH ([Y2007])
# R. Carragher
# Date: 25/11/2013


Id <- "$Id: c212ssBH.R,v 1.3 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.ssBH <- function(trial.data, alpha = 0.05) {

	if (is.null(trial.data)) {
		print("NULL trial data")
		return(NULL)
	}

	if (is.character(trial.data)) {
		file = trial.data
		trial.data = read.table(file, header = TRUE, stringsAsFactors=FALSE)
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

	m <- nrow(trial.data)

	if (is.null(m)) {
		print("No data found")
		return(NULL)
	}

	# No data in data frame - just return it
	if (m == 0) {
		return(trial.data)
	}

	rej <- trial.data[0,, drop = FALSE]

	B <- unique(trial.data$B)
	numB <- length(B)


	mg <- rep(0, numB)

	for (i in 1:numB) {
		g <- trial.data[trial.data$B == B[i],, drop = FALSE]
		mg[i] <- nrow(g)

		level <- (mg[i]/m)*alpha

		sig <- c212.BH(g, alpha = level)

		rej <- rbind(rej, sig)
	}

	return(rej)
}
