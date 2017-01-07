# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 28/04/2015


Id <- "$Id: c212.interim.1a.hier3.R,v 1.16 2016/12/08 13:55:26 clb13102 Exp clb13102 $"

c212.interim.1a.hier3 <- function(trial.data, sim_type = "SLICE", burnin = 10000,
	iter = 40000, nchains = 3,
	global.sim.params = data.frame(type = c("MH", "SLICE"), param = c("sigma_MH", "w"),
							value = c(0.2,1), control = c(0,6), stringsAsFactors = FALSE),
	sim.params = NULL,
	monitor = data.frame(variable = c("theta", "gamma", "mu.gamma", "mu.theta",
					"sigma2.theta", "sigma2.gamma",
		            "mu.theta.0", "mu.gamma.0", "tau2.theta.0", "tau2.gamma.0"),
					monitor = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
					stringsAsFactors = FALSE),
	initial_values = NULL, level = 1,
	hyper_params = list(mu.gamma.0.0 = 0, tau2.gamma.0.0 = 10, mu.theta.0.0 = 0,
		tau2.theta.0.0 = 10, alpha.gamma.0.0 = 3, beta.gamma.0.0 = 1,
		alpha.theta.0.0 = 3, beta.theta.0.0 = 1, alpha.gamma = 3, beta.gamma = 1,
		alpha.theta = 3, beta.theta = 1), memory_model = "HIGH"
	)
{
	if (level == 0) {
		model_fit = c212.interim.1a.indep(trial.data, sim_type, burnin,
						iter, nchains, global.sim.params, sim.params, monitor,
						initial_values, hyper_params, memory_model)
	}
	else if (level == 1) {
		model_fit = c212.interim.1a.dep.lev1(trial.data, sim_type, burnin,
							iter, nchains, global.sim.params, sim.params, monitor,
							initial_values, hyper_params, memory_model)
	} else if (level == 2) {
			model_fit = c212.interim.1a.dep.lev2(trial.data, sim_type, burnin,
							iter, nchains, global.sim.params, sim.params, monitor,
							initial_values, hyper_params, memory_model)
	} else {
		return(NULL)
	}

	return(model_fit)
}
