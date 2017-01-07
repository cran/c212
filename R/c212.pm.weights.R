c212.pointmass.weights = function(trial.data) {

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors = FALSE)
	}

	trial.data = trial.data[trial.data$Group == 2,]
	n = nrow(trial.data)

	if ("Interval" %in% names(trial.data)) {
		weights = data.frame(Interval = trial.data$Interval,
				I_index = trial.data$I_index, B = trial.data$B, AE = trial.data$AE,
				weight_pm = rep(0.5,n), stringsAsFactors = FALSE)
	}
	else {
		weights = data.frame(B = trial.data$B, AE = trial.data$AE,
								weight_pm = rep(0.5,n),
								stringsAsFactors = FALSE)
	}

	weights
}
