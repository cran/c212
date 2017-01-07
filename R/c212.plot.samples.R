# c212.plot.posterior
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 08/12/2014

Id <- "$Id: c212.plot.samples.R,v 1.2 2015/06/05 13:44:18 clb13102 Exp clb13102 $"


c212.plot.samples <- function(samples, title) {
	if (is.null(samples)) {
		print("NULL sample");
		return(NULL)
	}

	chains = 0
	d = dim(samples)

	mcmc_obj <- list(NA)

	if (is.null(d)) {
		mcmc_obj[[1]] <- mcmc(samples)
	}
	else {

		l = length(d)

		if (l > 2) {
			print("Dimension length error");
			return(NULL)
		}

		chains = d[1]

		for (i in 1:chains) {
			mcmc_obj[[i]] <- mcmc(samples[i,])
		}
	}

	m = mcmc.list(mcmc_obj)
	plot(m, main = title)
}
