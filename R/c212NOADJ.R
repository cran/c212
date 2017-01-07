# c212.NOADJ
# Case 2/12 - Unadjusted p-values
# R. Carragher
# Date: 25/11/2013

Id <- "$Id: c212NOADJ.R,v 1.2 2016/08/25 15:13:40 clb13102 Exp clb13102 $"

c212.NOADJ <- function(trial.data, alpha = 0.05)
{
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

	sig <- trial.data[trial.data$p <= alpha,, drop = FALSE]

	sig <- sig[order(sig$p),, drop=FALSE]

	row.names(sig) = NULL
	return(sig);
}
