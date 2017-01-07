# c212.interim
# Case 2/12: Interim Analysis MLE
# R. Carragher
# Date: 28/04/2015



c212.interim.MLE <- function(trial.data)
{
	Id <- "$Id: c212.interim.MLE.R,v 1.2 2015/07/28 13:46:34 clb13102 Exp clb13102 $"

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors = FALSE)
	}

	if (is.null(trial.data)) {
		print("NULL trial data")
		return(NULL)
	}

	# Perform some validation checks
	if ((is.null(nrow(trial.data))) || (nrow(trial.data) == 0)) {
		print("Missing trial data set");
		return(NULL)
	}
	# Check the correct columns are defined
	if (!("B" %in% colnames(trial.data))) {
		print("Missing B column");
		return(NULL)
	}
	if (!("AE" %in% colnames(trial.data))) {
		print("Missing AE column");
		return(NULL)
	}
	if (!("Count" %in% colnames(trial.data))) {
		print("Missing Count column");
		return(NULL)
	}
	if (!("Group" %in% colnames(trial.data))) {
		print("Missing Group column");
		return(NULL)
	}
	if (!("Interval" %in% colnames(trial.data))) {
		print("Missing Interval column");
		return(NULL)
	}
	if (!("I_index" %in% colnames(trial.data))) {
		print("Missing Interval column");
		return(NULL)
	}
	if (!("Exposure" %in% colnames(trial.data))) {
		print("Missing Exposure column");
		return(NULL)
	}

	# Order by body-system, adverse event, interval and group
	ordered.data <- trial.data[order(trial.data$I_index, trial.data$B, trial.data$AE,
							trial.data$Group),, drop=FALSE]

	cntrl.data <- ordered.data[ordered.data$Group == 1, ]
	treat.data <- ordered.data[ordered.data$Group == 2, ]

	# Check that cntrl, treat data are matched 1 to 1
	# Check that we have matching body-systems, AEs and intervals in the control and treatment groups
	# The data is ordered so a straight comparison is possible
	if(!identical(cntrl.data$B, treat.data$B)) {
		print("Mismatced body-system data");
		return(NULL)
	}
	if(!identical(cntrl.data$AE, treat.data$AE)) {
		print("Mismatced adverse event data");
		return(NULL)
	}
	if(!identical(cntrl.data$Interval, treat.data$Interval)) {
		print("Mismatced interval data");
		return(NULL)
	}

	# Size the data structures
	sz = nrow(cntrl.data)
	gamma_summ = data.frame(interval = rep(NA, sz), B = rep(NA, sz), AE = rep(NA, sz),
								estimate = rep(0, sz), lower = rep(0, sz),
								upper = rep(0, sz));

	theta_summ = data.frame(interval = rep(NA, sz), B = rep(NA, sz), AE = rep(NA, sz),
								estimate = rep(0, sz), lower = rep(0, sz),
								upper = rep(0, sz));

	control_rate_summ = data.frame(interval = rep(NA, sz), B = rep(NA, sz), AE = rep(NA, sz),
								estimate = rep(0, sz), lower = rep(0, sz),
								upper = rep(0, sz));

	treatment_rate_summ = data.frame(interval = rep(NA, sz), B = rep(NA, sz), AE = rep(NA, sz),
								estimate = rep(0, sz), lower = rep(0, sz),
								upper = rep(0, sz));


	for (i in 1:nrow(cntrl.data)) {
		estimate = cntrl.data[i,]$Count/cntrl.data[i,]$Exposure
		control_rate_summ[i,]$interval = cntrl.data[i,]$Interval
		control_rate_summ[i,]$B = cntrl.data[i,]$B
		control_rate_summ[i,]$AE = cntrl.data[i,]$AE
		control_rate_summ[i,]$estimate = estimate

		estimate = log(cntrl.data[i,]$Count/cntrl.data[i,]$Exposure)
		gamma_summ[i,]$interval = cntrl.data[i,]$Interval
		gamma_summ[i,]$B = cntrl.data[i,]$B
		gamma_summ[i,]$AE = cntrl.data[i,]$AE
		gamma_summ[i,]$estimate = estimate

		estimate = treat.data[i,]$Count/treat.data[i,]$Exposure
		treatment_rate_summ[i,]$interval = treat.data[i,]$Interval
		treatment_rate_summ[i,]$B = treat.data[i,]$B
		treatment_rate_summ[i,]$AE = treat.data[i,]$AE
		treatment_rate_summ[i,]$estimate = estimate

		estimate = log(treat.data[i,]$Count/treat.data[i,]$Exposure) - log(cntrl.data[i,]$Count/cntrl.data[i,]$Exposure)
		theta_summ[i,]$interval = treat.data[i,]$Interval
		theta_summ[i,]$B = treat.data[i,]$B
		theta_summ[i,]$AE = treat.data[i,]$AE
		theta_summ[i,]$estimate = estimate
	}

	model_fit = list(id = Id,
			gamma = gamma_summ,
			theta = theta_summ,
			control_rate = control_rate_summ,
			treatment_rate = treatment_rate_summ)
			
	attr(model_fit, "model") = "pois_MLE"

	return(model_fit)
}
