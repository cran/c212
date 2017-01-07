# c212.BB.ptheta
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.BB.ptheta.R,v 1.5 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.BB.ptheta <- function(raw)
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
	if (!("iter" %in% names(raw))) {
		print("Missing iter data");
		return(NULL)
	}
	if (!("burnin" %in% names(raw))) {
		print("Missing burnin data");
		return(NULL)
	}
	if (model == "BB") {
		if (!("pi" %in% names(raw))) {
			print("Missing pi data");
			return(NULL)
		}
		if (!("alpha.pi" %in% names(raw))) {
			print("Missing alpha.pi data");
			return(NULL)
		}
		if (!("beta.pi" %in% names(raw))) {
			print("Missing beta.pi data");
			return(NULL)
		}
	}

	nchains = raw$chains

	summ <- data.frame(B = character(0), AE = character(0), ptheta = numeric(0))

	# Inference - based on combined chains:\n")
	th <- array(NA, dim = c(raw$nBodySys, raw$maxAEs))
	theta_combined <- array(NA, dim=c(raw$nBodySys, raw$maxAEs, (raw$iter - raw$burnin)*nchains))
	for (b in 1:raw$nBodySys) {
		for (j in 1:raw$nAE[b]) {
			mcmc_obj <- list(NA)
			for (c in 1:nchains) {
				mcmc_obj[[c]] <- mcmc(raw$theta[c, b, j, ])
			}
			mlist <- mcmc.list(mcmc_obj)

			theta_combined[b,j,] <- c(raw$theta[1:nchains, b, j, ])
			s <- ecdf(theta_combined[b, j, ])
			th[b, j] <-  1 - s(0)
			row <- data.frame(B = raw$B[b], AE = raw$AE[b, j], ptheta = th[b, j])
			summ = rbind(summ, row)
		}
	}

	rownames(summ) <- NULL
	return(summ)
}
