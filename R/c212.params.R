c212.sim.control.params = function(trial.data, model = "1a") {

	if (model == "1a") {
		sim_param = c212.1a.sim.control.params(trial.data)
	}
	else {
		sim_param = c212.BB.sim.control.params(trial.data)
	}

	sim_param
}

c212.1a.sim.control.params = function(trial.data) {

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors = FALSE)
	}

	trial.data = trial.data[trial.data$Group == 2,]
	n = nrow(trial.data)

	if ("Interval" %in% names(trial.data)) {

		p_SLICE = "w"
		v_SLICE = 1.0
		c_SLICE = 6

		p_MH = "sigma_MH"
		v_MH = 0.2
		c_MH = 0.0

		sp1 = data.frame(type = "SLICE", variable = "gamma",
				Interval = trial.data$Interval,
				I_index = trial.data$I_index, B = trial.data$B,
				AE = trial.data$AE,
				param = p_SLICE, value = v_SLICE, control = c_SLICE,
				stringsAsFactors = FALSE)

		sp2 = data.frame(type = "SLICE", variable = "theta",
				Interval = trial.data$Interval,
				I_index = trial.data$I_index, B = trial.data$B,
				AE = trial.data$AE,
				param = p_SLICE, value = v_SLICE, control = c_SLICE,
				stringsAsFactors = FALSE)

		sp3 = data.frame(type = "MH", variable = "gamma",
				Interval = trial.data$Interval,
				I_index = trial.data$I_index, B = trial.data$B,
				AE = trial.data$AE, 
				param = p_MH, value = v_MH, control = c_MH,
				stringsAsFactors = FALSE)

		sp4 = data.frame(type = "MH", variable = "theta",
				Interval = trial.data$Interval,
				I_index = trial.data$I_index, B = trial.data$B,
				AE = trial.data$AE, 
				param = p_MH, value = v_MH, control = c_MH,
				stringsAsFactors = FALSE)
	}
	else {

		p_SLICE = "w"
		v_SLICE = 1.0
		c_SLICE = 6

		p_MH = "sigma_MH"
		v_MH = 0.35
		c_MH = 0.0

		sp1 = data.frame(type = "SLICE", variable = "gamma",
								B = trial.data$B,
								AE = trial.data$AE, 
								param = p_SLICE, value = v_SLICE, control = c_SLICE,
								stringsAsFactors = FALSE)

		sp2 = data.frame(type = "SLICE", variable = "theta",
								B = trial.data$B,
								AE = trial.data$AE, 
								param = p_SLICE, value = v_SLICE, control = c_SLICE,
								stringsAsFactors = FALSE)

		sp3 = data.frame(type = "MH", variable = "gamma",
								B = trial.data$B,
								AE = trial.data$AE, 
								param = p_MH, value = v_MH, control = c_MH,
								stringsAsFactors = FALSE)

		sp4 = data.frame(type = "MH", variable = "theta",
								B = trial.data$B,
								AE = trial.data$AE, 
								param = p_MH, value = v_MH, control = c_MH,
								stringsAsFactors = FALSE)
	}
	sim_params = rbind(sp1, sp2, sp3, sp4)

	sim_params
}

c212.BB.sim.control.params = function(trial.data) {

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors = FALSE)
	}

	trial.data = trial.data[trial.data$Group == 2,]

	if ("Interval" %in% names(trial.data)) {

		v_sigma_MH_gamma = 0.2
		v_sigma_MH_theta = 0.15

		v_w = 1
		v_w_control = 6
		v_sigma = 3

		sp1 = data.frame(type = "MH", variable = "gamma",
								Interval = trial.data$Interval,
								I_index = trial.data$I_index,
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "sigma_MH_gamma", value = v_sigma_MH_gamma,
								control = 0,
								stringsAsFactors = FALSE)

		sp2 = data.frame(type = "MH", variable = "theta",
								Interval = trial.data$Interval,
								I_index = trial.data$I_index,
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "sigma_MH_theta", value = v_sigma_MH_theta,
								control = 0,
								stringsAsFactors = FALSE)

		sp3 = data.frame(type = "SLICE", variable = "gamma",
								Interval = trial.data$Interval,
								I_index = trial.data$I_index,
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "w_gamma", value = v_w,
								control = v_w_control,
								stringsAsFactors = FALSE)

		sp4 = data.frame(type = "MH", variable = "alpha",
								Interval = trial.data$Interval,
								I_index = trial.data$I_index,
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "sigma_MH_alpha", value = v_sigma,
								control = 0,
								stringsAsFactors = FALSE)

		sp5 = data.frame(type = "SLICE", variable = "alpha",
								Interval = trial.data$Interval,
								I_index = trial.data$I_index,
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "w_alpha", value = v_w,
								control = v_w_control,
								stringsAsFactors = FALSE)

		sp6 = data.frame(type = "MH", variable = "beta",
								Interval = trial.data$Interval,
								I_index = trial.data$I_index,
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "sigma_MH_beta", value = v_sigma,
								control = 0,
								stringsAsFactors = FALSE)

		sp7 = data.frame(type = "SLICE", variable = "beta",
								Interval = trial.data$Interval,
								I_index = trial.data$I_index,
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "w_beta", value = v_w,
								control = v_w_control,
								stringsAsFactors = FALSE)

		sim_params = rbind(sp1, sp2, sp3, sp4, sp5, sp6, sp7)

	}
	else {

		v_sigma_MH_gamma = 0.2
		v_sigma_MH_theta = 0.2

		v_w_gamma = 1
		v_w_gamma_control = 6

		sp1 = data.frame(type = "MH", variable = "gamma",
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "sigma_MH_gamma", value = v_sigma_MH_gamma,
								control = 0,
								stringsAsFactors = FALSE)

		sp2 = data.frame(type = "MH", variable = "theta",
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "sigma_MH_theta", value = v_sigma_MH_theta,
								control = 0,
								stringsAsFactors = FALSE)

		sp3 = data.frame(type = "SLICE", variable = "gamma",
								B = trial.data$B,
								AE = trial.data$AE, 
								param = "w_gamma", value = v_w_gamma,
								control = v_w_gamma_control,
								stringsAsFactors = FALSE)

		sim_params = rbind(sp1, sp2, sp3)
	}

	sim_params
}

