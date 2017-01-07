# c212.BB.ptheta
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.ptheta.R,v 1.1 2015/07/07 11:55:20 clb13102 Exp clb13102 $"

c212.ptheta <- function(raw)
{
	if (is.null(raw)) {
		print("NULL raw data");
		return(NULL)
	}

	model = attr(raw, "model")

	if (is.null(model)) {
		print("Missing model attribute");
		return(NULL)
	}

	 if (model == "1a" || model == "BB") {
		summ = c212.BB.ptheta(raw)
	}
	else {
		summ = c212.interim.ptheta(raw)
	}

	return(summ)
}
