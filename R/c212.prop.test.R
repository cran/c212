# R. Carragher
# Date: 09/02/2015
#

Id <- "$Id: c212.prop.test.R,v 1.3 2017/01/04 12:20:56 clb13102 Exp clb13102 $"

c212.bin.test <- function(trial.data, alternative = "two.sided", correct = TRUE)
{

	if (is.character(trial.data)) {
		file = trial.data
		trial.data <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
	}

	if (is.null(trial.data)) {
		print("NULL trial data")
		return(NULL)
	}

	# Perform some validation checks
	if (is.null(nrow(trial.data)) || (nrow(trial.data) == 0)) {
		print("Missing trial data");
		return (NULL)
	}

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
	if (!("Total" %in% colnames(trial.data))) {
		print("Missing Total data");
		return(NULL)
	}

	ordered.data <- trial.data[order(trial.data$B, trial.data$AE, trial.data$Group),, drop=FALSE]
	cntrl.data <- ordered.data[ordered.data$Group == 1, ]
	treat.data <- ordered.data[ordered.data$Group == 2, ]

	if(!identical(cntrl.data$B, treat.data$B)) {
		print("Mismatch in body-systems data");
		return(NULL)
	}
	if(!identical(cntrl.data$AE, treat.data$AE)) {
		print("Mismatch in adverse events data");
		return(NULL)
	}

	n = nrow(cntrl.data)
	f <- matrix(nrow = 2, ncol = 2)

	prop = data.frame(B = character(0), j = integer(0), AE = character(0), p = numeric(0), stringsAsFactors=FALSE)

	B = treat.data[1, ]$B
	j = 1

	for (i in 1:n) {

		if (B != treat.data[i,]$B) {
			B = treat.data[i,]$B
			j = 1
		}
		

		f[1,1] <- treat.data[i,]$Count 
		f[1,2] <- treat.data[i,]$Total - treat.data[i,]$Count
		f[2,1] <- cntrl.data[i,]$Count
		f[2,2] <- cntrl.data[i,]$Total - cntrl.data[i,]$Count
		res = prop.test(f, alternative = alternative, correct = correct)

		row = data.frame(B = treat.data[i,]$B, j = j, AE = treat.data[i,]$AE, p = res$p.value, stringsAsFactors=FALSE)
		j = j + 1

		prop = rbind(prop, row)

	}

	return(prop)
}
