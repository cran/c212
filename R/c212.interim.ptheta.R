# c212.BB.ptheta
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.interim.ptheta.R,v 1.5 2019/05/05 13:18:12 clb13102 Exp clb13102 $"

c212.interim.ptheta <- function(raw)
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

	if (!("chains" %in% names(raw))) {
		print("Missing chains data");
		return(NULL)
	}
	if (!("maxBs" %in% names(raw))) {
		print("Missing chains data");
		return(NULL)
	}
	if (!("nBodySys" %in% names(raw))) {
		print("Missing chains data");
		return(NULL)
	}
	if (!("maxAEs" %in% names(raw))) {
		print("Missing chains data");
		return(NULL)
	}
	if (!("nAE" %in% names(raw))) {
		print("Missing nAE data");
		return(NULL)
	}
	if (!("theta" %in% names(raw))) {
		print("Missing theta data");
		return(NULL)
	}
	if (!("B" %in% names(raw))) {
		print("Missing B data");
		return(NULL)
	}
	if (!("AE" %in% names(raw))) {
		print("Missing AE data");
		return(NULL)
	}

	nchains = raw$chains

	summ <- data.frame(interval = character(0), B = character(0), AE = character(0), ptheta = numeric(0))

	# Inference - based on combined chains:\n")
	samples_combined = rep(NA, (raw$iter - raw$burnin)*nchains)

	for (i in 1:raw$nIntervals) {
		for (b in 1:raw$nBodySys[i]) {
			for (j in 1:raw$nAE[i, b]) {
				mcmc_obj <- list(NA)
				for (c in 1:nchains) {
					mcmc_obj[[c]] <- mcmc(raw$theta[c, i, b, j, ])
				}
				mlist <- mcmc.list(mcmc_obj)

				samples_combined <- c(raw$theta[1:nchains, i, b, j, ])
				s <- ecdf(samples_combined)
				th <-  1 - s(0)
				row <- data.frame(interval = raw$Intervals[i], B = raw$B[i, b], AE = raw$AE[i, b, j], ptheta = th)
				summ = rbind(summ, row)
			}
		}
	}

	rownames(summ) <- NULL
	return(summ)
}
