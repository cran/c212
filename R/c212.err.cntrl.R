# c212.err.cntrl
# Case 2/12
# R. Carragher
# Date: 19/07/2019

Id <- "$Id: c212.err.cntrl.R,v 1.2 2020/08/31 15:12:07 clb13102 Exp clb13102 $"


c212.err.cntrl <- function(trial.data, alpha = 0.05, method = "NOADJ",...)
{
	if (is.null(trial.data)) {
		print("NULL trial data")
		return(NULL)
	}

	if (is.character(trial.data)) {
		file = trial.data
		trial.data = read.table(file, header = TRUE, stringsAsFactors=FALSE)
	}

	if (!is.numeric(alpha)) {
		print("alpha value must be numberic")
		return(NULL)
	}

	err.cntrl.grp <- c("GBH", "ssBH", "DFDR")
	err.cntrl <- c("NOADJ", "BONF", "BH", err.cntrl.grp)
	if (!(method %in% err.cntrl)) {
		print(sprintf("Invalid method: %s", method))
		return(NULL)
	}

	# Check the correct columns are defined
	cols <- colnames(trial.data)
	if (!("p" %in% cols)) {
		print("p column missing")
		return(NULL);
	}

	if ((method %in% err.cntrl.grp) & (!("B" %in% cols))) {
		print("B column missing")
		return(NULL);
	}

	sig <- NA
	if (method == "NOADJ")
		sig <- c212.NOADJ(trial.data = trial.data, alpha = alpha)
	else if (method == "BONF")
		sig <- c212.BONF(trial.data = trial.data, alpha = alpha)
	else if (method == "BH")
		sig <- c212.BH(trial.data = trial.data, alpha = alpha)
	else if (method == "GBH") {
		arg <- list(...)
		if (is.null(arg$pi0))
			sig <- c212.GBH(trial.data = trial.data, alpha = alpha)
		else
			sig <- c212.GBH(trial.data = trial.data, pi0 = arg$pi0, alpha = alpha)
	}
	else if (method == "ssBH")
		sig <- c212.ssBH(trial.data = trial.data, alpha = alpha)
	else if (method == "DFDR")
		sig <- c212.DFDR(trial.data = trial.data, alpha = alpha)
	else if (method == "BONF")
		sig <- c212.BONF(trial.data = trial.data, alpha = alpha)

	row.names(sig) = NULL
	return(sig);
}
