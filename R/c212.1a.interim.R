# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 03/08/2020


Id <- "$Id: c212.1a.interim.R,v 1.1 2020/08/31 10:04:49 clb13102 Exp clb13102 $"

c212.1a.interim <- function(trial.data, sim_type = "SLICE", burnin = 10000,
	iter = 40000, nchains = 3,
	global.sim.params = NULL,
	sim.params = NULL,
	monitor = NULL,
	initial_values = NULL, hier = 3, level = 1,
	hyper_params = NULL,
	 memory_model = "HIGH"
	)
{
	g.s.p <-global.sim.params
	mon <- monitor
	h.p <- hyper_params
	model.fit <- NULL

	if (is.null(g.s.p)) {
		g.s.p  = c212.global.sim.params(trial.data, model = "1a", hier = hier)
	}
	if (is.null(mon)) {
		mon <- c212.monitor.samples(model = "1a", hier = hier)
		mon$monitor <- 1
	}
	if (is.null(h.p)) {
		h.p <- c212.hyper.params(trial.data, model = "1a", hier = hier)
	}

	if (hier == 3) {
		model.fit <- c212.interim.1a.hier3(trial.data, sim_type, burnin,
						iter, nchains,
						global.sim.params = g.s.p,
						sim.params,
						monitor = mon,
						initial_values,
						level,
						hyper_params = h.p,
						memory_model)
	}
	else if (hier == 2) {
		model.fit <- c212.interim.1a.hier2(trial.data, sim_type, burnin,
						iter, nchains,
						global.sim.params = g.s.p,
						sim.params,
						monitor = mon,
						initial_values,
						level,
						hyper_params = h.p,
						memory_model)
	}

	return(model.fit)
}
