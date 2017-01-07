# c212.BB.plot.eot.data
# Case 2/12 Model c212.BB
# R. Carragher
# Date: 28/11/2014

Id <- "$Id: c212.plot.eot.data.R,v 1.5 2016/10/14 10:39:05 clb13102 Exp clb13102 $"


c212.plot.eot.data <- function(trial.data, legend = TRUE, interactive = FALSE, cex = 0.5)
{

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
	# Check the correct columns are defined
	if (!("B" %in% colnames(trial.data))) {
		print("Missing B data");
		return(NULL)
	}
	if (!("AE" %in% colnames(trial.data))) {
		print("Missing AE data");
		return(NULL)
	}
	if (!("Count" %in% colnames(trial.data))) {
		print("Missing Count data");
		return(NULL)
	}
	if (!("Group" %in% colnames(trial.data))) {
		print("Missing Group data");
		return(NULL)
	}

	# Order by body-system, AE and group
	ordered.data <- trial.data[order(trial.data$B, trial.data$AE, trial.data$Group),, drop=FALSE]

	cntrl.data <- ordered.data[ordered.data$Group == 1, ]
	treat.data <- ordered.data[ordered.data$Group == 2, ]

	# Some further validation
	# Check that we have matching body-systems and AEs in the control and treatment groups
	# The data is ordered so a straight comparison is possible
	if(!identical(cntrl.data$B, treat.data$B)) {
		print("Mismatch in body-systems data");
		return(NULL)
	}
	if(!identical(cntrl.data$AE, treat.data$AE)) {
		print("Mismatch in adverse events data");
		return(NULL)
	}


	B = unique(ordered.data$B)

	cntrl_b_tot = rep(0, length(B))
	treat_b_tot = rep(0, length(B))

	for (i in 1:length(B)) {
		c = cntrl.data[cntrl.data$B == B[i],, drop = FALSE]
		t = treat.data[treat.data$B == B[i],, drop = FALSE]
		cntrl_b_tot[i] = sum(c$Count)
		treat_b_tot[i] = sum(t$Count)
	}

	old.par = par(mfrow=c(2,1))
	#par(mfrow=c(3,1))

	plot(1:length(B), cntrl_b_tot, xaxt = 'n', xlab = "Body Systems", ylab="AE Incidence Counts", col="green", bty='L',
	main = "AE Incidence Counts by Body-System", ylim=c(0, max(treat_b_tot, cntrl_b_tot)))
	points(1:length(B), cntrl_b_tot, lwd = 0.4, col="green")
	lines(1:length(B), cntrl_b_tot, lwd = 0.4, col="green")
	points(1:length(B), treat_b_tot, lwd = 0.4, col="red")
	lines(1:length(B), treat_b_tot, lwd = 0.4, col="red")
	axis(1, at=1:length(B), labels=B)


	#legend(4.5,120,
	if (legend) {
		legend("topright",
			c("Control","Treatment"),
			lty=c(1,1),
			lwd=c(0.4,0.4),col=c("green","red"), bg="white")
	}

	# Plot with lines per body-system
	plot(1:length(cntrl.data$Count), cntrl.data$Count, xaxt = 'n', xlab = "Adverse Events", ylab="AE Incidence Counts", col="green", bty='L', main = "Individual AE Incidence Counts", type='n', ylim=c(0, max(cntrl.data$Count, treat.data$Count)))

	cnt = 0
	abline(v = cnt + 0.5, lwd = 0.5, lty=2, col="black")
	for (i in 1:length(B)) {
		c = cntrl.data[cntrl.data$B == B[i],,drop=FALSE]
		t = treat.data[treat.data$B == B[i],,drop=FALSE]
		points((cnt +1):(cnt + length(c$Count)), c$Count, lwd = 0.4, col="green")
		lines((cnt+1):(cnt + length(c$Count)), c$Count, lwd = 0.4, col="green")
		points((cnt +1):(cnt + length(t$Count)), t$Count, lwd = 0.4, col="red")
		lines((cnt+1):(cnt + length(t$Count)), t$Count, lwd = 0.4, col="red")
		cnt = cnt + length(c$Count)
		abline(v = cnt + 0.5, lwd = 0.5, lty=2, col="black")
	}
	axis(1, at=1:length(cntrl.data$Count), labels=FALSE)

	max_text = max(nchar(as.character(cntrl.data$AE)))
	labels = as.character(cntrl.data$AE)

	x0 = 0
	par_usr = 0.2
	if (max_text >= 10) {
		labels = substr(labels, 1, 10)
		substr(labels, 8, 10) <- "..."
		par_usr = 4.0
		x0 = 0
	}

	text(1:length(cntrl.data$Count) - x0, labels= labels, srt=80, pos = 1, par("usr")[3] - par_usr, xpd=TRUE, cex = cex)

	if (legend) {
			legend("topright",
			c("Control","Treatment"),
			lty=c(1,1),
			lwd=c(0.4,0.4),col=c("green","red"), bg="white")
	}
	
	if (interactive) {
		x = c(1:length(cntrl.data$Count), 1:length(treat.data$Count))
		y = c(cntrl.data$Count, treat.data$Count)
		labels = c(as.character(cntrl.data$AE), as.character(treat.data$AE))
		n = identify(x, y, labels = labels, cex = 0.75)
	}

	par(old.par)
}
