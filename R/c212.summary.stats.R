# c212.BB.summary
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.summary.stats.R,v 1.6 2016/03/18 15:39:32 clb13102 Exp clb13102 $"

c212.summary.stats <- function(raw)
{
	if (is.null(raw)) {
		print("NULL raw data");
		return(NULL)
	}

	model = attr(raw, "model")
	if (is.null(model)) {
		print("Model attribute missing");
		return(NULL)
	}

	summary.stats = list()

	if (model == "1a" || model == "BB") {
		summary.stats = c212.BB.summary.stats(raw)
	} else if (model == "1a_pois_indep") {
		summary.stats = c212.interim.1a.indep.summary.stats(raw)
	} else if (model == "1a_pois_dep_lev2") {
		summary.stats = c212.interim.1a.dep.lev2.summary.stats(raw)
	} else if (model == "1a_pois_dep_lev1") {
		summary.stats = c212.interim.1a.dep.lev1.summary.stats(raw)
	} else if (model == "BB_pois_indep") {
		summary.stats = c212.interim.BB.indep.summary.stats(raw)
	} else if (model == "BB_pois_dep_lev2") {
		summary.stats = c212.interim.BB.dep.lev2.summary.stats(raw)
	} else if (model == "BB_pois_dep_lev1") {
		summary.stats = c212.interim.BB.dep.lev1.summary.stats(raw)
	} else if (model == "1a_pois_h2_l0") {
		summary.stats = c212.interim.1a.hier2.lev0.summary.stats(raw)
	} else if (model == "1a_pois_h2_l1") {
		summary.stats = c212.interim.1a.hier2.lev1.summary.stats(raw)
	} else if (model == "BB_pois_h2_l0") {
		summary.stats = c212.interim.BB.hier2.lev0.summary.stats(raw)
	} else if (model == "BB_pois_h2_l1") {
		summary.stats = c212.interim.BB.hier2.lev1.summary.stats(raw)
	}
	return(summary.stats)
}

c212.print.summary.stats <- function(summ)
{
	if (is.null(summ)) {
		print("NULL summary data");
		return(NULL)
	}

	model = attr(summ, "model")
	if (is.null(model)) {
		print("Missing model attribute");
		return(NULL)
	}

	if (model == "1a" || model == "BB") {
		c212.BB.print.summary.stats(summ)
	} else if (model == "1a_pois_indep") {
		c212.interim.1a.indep.print.summary.stats(summ)
	} else if (model == "1a_pois_dep_lev2") {
		c212.interim.1a.dep.lev2.print.summary.stats(summ)
	} else if (model == "1a_pois_dep_lev1") {
		c212.interim.1a.dep.lev1.print.summary.stats(summ)
	} else if (model == "BB_pois_indep") {
		c212.interim.BB.indep.print.summary.stats(summ)
	} else if (model == "BB_pois_dep_lev2") {
		c212.interim.BB.dep.lev2.print.summary.stats(summ)
	} else if (model == "BB_pois_dep_lev1") {
		c212.interim.BB.dep.lev1.print.summary.stats(summ)
	} else if (model == "1a_pois_h2_l0") {
		c212.interim.1a.hier2.lev0.print.summary.stats(summ)
	} else if (model == "1a_pois_h2_l1") {
		c212.interim.1a.hier2.lev1.print.summary.stats(summ)
	} else if (model == "BB_pois_h2_l0") {
		c212.interim.BB.hier2.lev0.print.summary.stats(summ)
	} else if (model == "BB_pois_h2_l1") {
		c212.interim.BB.hier2.lev1.print.summary.stats(summ)
	}
}
