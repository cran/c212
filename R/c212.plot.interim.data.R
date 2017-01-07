# c212.plot.interim.data
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.plot.interim.data.R,v 1.1 2016/08/25 15:13:41 clb13102 Exp clb13102 $"

c212.plot.interim.data <- function(trial.data, body_sys, cex = 0.8, title = NULL) {

	if (is.null(trial.data)) {
		print("NULL trial data")
		return(NULL)
	}

	if (is.character(trial.data)) {
		file = trial.data
		trial.data = read.table(file, header=TRUE, stringsAsFactors = FALSE)
	}

	if (is.null(nrow(trial.data)) || (nrow(trial.data) == 0)) {
		print("Missing trial data");
		return(NULL)
	}

	# Some validation
	n = c("Interval", "I_index", "B", "AE", "Group", "Count", "Exposure")

	if (M_global$checkCols(n, trial.data)) {
		print("Missing required columns")
		return(NULL)
	}

	trial.data = trial.data[trial.data$B == body_sys, ]

	ordered.data <- trial.data[order(trial.data$I_index, trial.data$B, trial.data$AE, trial.data$Group),, drop=FALSE]

	cntrl.data <- ordered.data[ordered.data$Group == 1, ]
	treat.data <- ordered.data[ordered.data$Group == 2, ]

	aes = as.character(unique(cntrl.data$AE))
	intervals = as.character(unique(cntrl.data$Interval))

	linetype <- c(1:length(aes))
	linetype_cntrl <- rep(2, length(aes))
	linetype_trt <- rep(1, length(aes))
	colors_cntrl <- c(1:length(aes))
	colors_trt <- c(1:length(aes))

	control = matrix(cntrl.data$Count, nrow = length(aes), ncol = length(intervals), byrow = FALSE)
	treatment = matrix(treat.data$Count, nrow = length(aes), ncol = length(intervals), byrow = FALSE)

	max_count = max(control, treatment)

	if (is.null(title)) {
		title = sprintf("Adverse Events Counts by Interval\n%s", body_sys)
	}

	plot(1:length(intervals), 1:length(intervals), type="n", xaxt = 'n', ylab="Adverse Event Counts", xlab = "Intervals",
			main = title, , ylim=c(0, max(control, treatment)))

	for (i in 1:length(aes)) {
		cntrl = control[i,]
		points(1:length(intervals), cntrl, lwd = 0.4, col=colors_cntrl[i], pch = 1)
		lines(1:length(intervals), cntrl, lwd = 0.4, col=colors_cntrl[i], lty = linetype_cntrl[i])
	}
	for (i in 1:length(aes)) {
		trt = treatment[i,]
		points(1:length(intervals), trt, lwd = 0.4, col=colors_trt[i], pch = 5)
		lines(1:length(intervals), trt, lwd = 0.4, col=colors_trt[i], lty = linetype_trt[i])
	}
	axis(1, at=1:length(intervals), labels=intervals)

	L = legend("topright",
		c(aes),
		lty=c(linetype_trt),
		pch=rep(5, length(aes)),
		lwd=rep(0.4,length(aes)),col=colors_trt, bg="white", cex = 0.8, title="Treatment")
	legend(x = L$rect$left - L$rect$w, y = L$rect$top,
		c(aes),
		lty=c(linetype_cntrl),
		pch=rep(1, length(aes)),
		lwd=rep(0.4,length(aes)),col=colors_cntrl, bg="white", cex = 0.8, title="Control")
}
