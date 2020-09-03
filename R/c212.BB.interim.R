# c212.interim
# Case 2/12: Interim Analysis wrapper
# R. Carragher
# Date: 05/06/2015


Id <- "$Id: c212.BB.interim.R,v 1.1 2020/08/31 10:04:55 clb13102 Exp clb13102 $"

c212.BB.interim <- function(trial.data,
	sim_type = "SLICE", burnin = 20000, iter = 60000, nchains = 5,
	theta_algorithm = "MH",
	global.sim.params = NULL,
	sim.params = NULL,
	monitor = NULL,
	initial_values = NULL,
	hier = 3,
	level = 1,
	hyper_params = NULL,
	global.pm.weight = 0.5,
	pm.weights = NULL,
	adapt_phase=1, memory_model = "HIGH")
{
	g.s.p <-global.sim.params
 	mon <- monitor
 	h.p <- hyper_params
	model.fit <- NULL

	if (is.null(g.s.p)) {
		g.s.p  = c212.global.sim.params(trial.data, model = "BB", hier = hier)
	}
	if (is.null(mon)) {
		mon <- c212.monitor.samples(model = "BB", hier = hier)
		mon$monitor <- 1
    }
    if (is.null(h.p)) {
		h.p <- c212.hyper.params(trial.data, model = "BB", hier = hier)
    }

    if (hier == 3) {
		model.fit <- c212.interim.BB.hier3(trial.data, sim_type,
						burnin, iter, nchains, theta_algorithm,
						global.sim.params = g.s.p,
						sim.params,
						monitor = mon,
						initial_values,
						level,
						hyper_params = h.p,
						global.pm.weight,
						pm.weights,
						adapt_phase, memory_model)
	}
	else if (hier == 2) {
		model.fit <- c212.interim.BB.hier2(trial.data, sim_type,
						burnin, iter, nchains, theta_algorithm,
						global.sim.params = g.s.p,
						sim.params,
						monitor = mon,
						initial_values,
						level,
						hyper_params = h.p,
						global.pm.weight,
						pm.weights,
						adapt_phase, memory_model)
	}
 
	return(model.fit)
}
