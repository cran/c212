#c212.BB.convergence.diag # Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014
#
# If the MCMC simulation has been run for more than one chain report the Gelman-Rubin statistic.
# If the MCMC simulation has been run for only one chain report the Geweke diagnostic (Z-score)
#

Id <- "$Id: c212.convergence.R,v 1.11 2016/10/14 10:39:04 clb13102 Exp clb13102 $"

c212.convergence.diag <- function(raw, debug_diagnostic = FALSE)
{
	if (is.null(raw)) {
		print("NULL raw data")
		return(NULL)
	}

	model = attr(raw, "model")
	if (is.null(model)) {
		print("Simulation model attribute missing")
		return(NULL)
	}

	conv.diag = list()

	if (model == "1a" || model == "BB") {
		conv.diag = c212.BB.convergence.diag(raw)
	} else if (model == "1a_pois_indep") {
		conv.diag = c212.interim.1a.indep.convergence.diag(raw, debug_diagnostic);
	} else if (model == "1a_pois_dep_lev2") {
		conv.diag = c212.interim.1a.dep.lev2.convergence.diag(raw, debug_diagnostic);
	} else if (model == "1a_pois_dep_lev1") {
		conv.diag = c212.interim.1a.dep.lev1.convergence.diag(raw, debug_diagnostic);
	} else if (model == "BB_pois_indep") {
		conv.diag = c212.interim.BB.indep.convergence.diag(raw, debug_diagnostic);
	} else if (model == "BB_pois_dep_lev2") {
		conv.diag = c212.interim.BB.dep.lev2.convergence.diag(raw, debug_diagnostic);
	} else if (model == "BB_pois_dep_lev1") {
		conv.diag = c212.interim.BB.dep.lev1.convergence.diag(raw, debug_diagnostic);
	} else if (model == "1a_pois_h2_l0") {
		conv.diag = c212.interim.1a.hier2.lev0.convergence.diag(raw, debug_diagnostic);
	} else if (model == "1a_pois_h2_l1") {
		conv.diag = c212.interim.1a.hier2.lev1.convergence.diag(raw, debug_diagnostic);
	} else if (model == "BB_pois_h2_l0") {
		conv.diag = c212.interim.BB.hier2.lev0.convergence.diag(raw, debug_diagnostic);
	} else if (model == "BB_pois_h2_l1") {
		conv.diag = c212.interim.BB.hier2.lev1.convergence.diag(raw, debug_diagnostic);
	}
	else {
		conv.diag = NULL
	}

	return(conv.diag)
}

c212.print.convergence.summary <- function(conv) {

	if (is.null(conv)) {
		print("NULL conv data")
		return(NULL)
	}

	model = attr(conv, "model")
	if (is.null(model)) {
		print("Convergence model attribute missing")
		return(NULL)
	}

	if (model == "1a" || model == "BB") {
		c212.BB.print.convergence.summary(conv)
	} else if (model == "1a_pois_indep") {
		c212.interim.1a.indep.print.convergence.summary(conv)
	} else if (model == "1a_pois_dep_lev2") {
		c212.interim.1a.dep.lev2.print.convergence.summary(conv)
	} else if (model == "1a_pois_dep_lev1") {
		c212.interim.1a.dep.lev1.print.convergence.summary(conv)
	} else if (model == "BB_pois_indep") {
		c212.interim.BB.indep.print.convergence.summary(conv)
	} else if (model == "BB_pois_dep_lev2") {
		c212.interim.BB.dep.lev2.print.convergence.summary(conv)
	} else if (model == "BB_pois_dep_lev1") {
		c212.interim.BB.dep.lev1.print.convergence.summary(conv)
	} else if (model == "1a_pois_h2_l0") {
		c212.interim.1a.hier2.lev0.print.convergence.summary(conv)
	} else if (model == "1a_pois_h2_l1") {
		c212.interim.1a.hier2.lev1.print.convergence.summary(conv)
	} else if (model == "BB_pois_h2_l0") {
		c212.interim.BB.hier2.lev0.print.convergence.summary(conv)
	} else if (model == "BB_pois_h2_l1") {
		c212.interim.BB.hier2.lev1.print.convergence.summary(conv)
	}
}
