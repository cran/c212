\name{c212.ptheta}
\alias{c212.ptheta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reports the posterior probability that theta (the increase in the log-odds) is greater than zero for each Adverse Event}
\description{
This function reports the posterior probability that theta is positive, i.e. that there is an increase in log
odds of an adverse event being associated with treatment.
}
\usage{
	c212.ptheta(raw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
The output from a call to c212.BB.
}
}
\value{
A data frame containing the columns: \emph{interval} if analysing interim data, \emph{B}: body system, \emph{AE}: adverse event and \emph{ptheta}, the posterior probability that theta is positive.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
%%trial.data <- data.frame(B = c("BS6", "BS6", "BS7", "BS7", "BS6", "BS6",
%%"BS5", "BS7", "BS6", "BS2", "BS2", "BS5", "BS4", "BS1", "BS6", "BS5",
%%"BS4", "BS4", "BS2", "BS5", "BS1", "BS6", "BS8", "BS6", "BS3", "BS8",
%%"BS8", "BS5", "BS6", "BS8", "BS5", "BS2", "BS6", "BS6", "BS3", "BS3",
%%"BS8", "BS8", "BS6", "BS4", "BS5", "BS2", "BS3", "BS4", "BS5", "BS2",
%%"BS3", "BS5", "BS6", "BS4", "BS8", "BS2", "BS6", "BS5", "BS5", "BS2",
%%"BS3", "BS8", "BS6", "BS5", "BS8", "BS7", "BS3", "BS6", "BS4", "BS6",
%%"BS6", "BS3", "BS4", "BS6", "BS3", "BS5", "BS5", "BS5", "BS6", "BS3",
%%"BS6", "BS3", "BS7", "BS3", "BS3", "BS5", "BS8", "BS8", "BS4", "BS7",
%%"BS6", "BS4", "BS3", "BS8"),
%%AE = c("AE30", "AE33", "AE38", "AE39", "AE31", "AE31", "AE25", "AE37",
%%"AE28", "AE3", "AE3", "AE24", "AE15", "AE1", "AE35", "AE18", "AE16",
%%"AE17", "AE2", "AE23", "AE1", "AE35", "AE42", "AE27", "AE7", "AE43",
%%"AE42", "AE24", "AE36", "AE40", "AE20", "AE4", "AE29", "AE27", "AE912",
%%"AE910", "AE41", "AE45", "AE26", "AE16", "AE21", "AE5", "AE912", "AE14",
%%"AE21", "AE5", "AE911", "AE23", "AE33", "AE14", "AE40", "AE4", "AE36",
%%"AE22", "AE22", "AE2", "AE9", "AE45", "AE32", "AE25", "AE44", "AE39",
%%"AE7", "AE34", "AE17", "AE34", "AE28", "AE6", "AE15", "AE29", "AE8",
%%"AE19", "AE19", "AE20", "AE26", "AE8", "AE30", "AE910", "AE38", "AE6",
%%"AE9", "AE18", "AE43", "AE44", "AE13", "AE37", "AE32", "AE13", "AE911",
%%"AE41"),
%%Group = c(2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2, 1, 2,
%%1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 2, 1,
%%2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2,
%%1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 2,
%%1, 2, 2, 1, 2, 2, 1, 1, 1, 1),
%%Count = c(20, 15, 25, 22, 24, 19, 49, 16, 12, 23, 18, 52, 19, 21, 18,
%%51, 24, 22, 48, 57, 17, 20, 24, 20, 57, 20, 31, 56, 25, 19, 46, 25, 24,
%%19, 22, 29, 20, 31, 28, 19, 66, 14, 60, 18, 51, 22, 52, 58, 22, 17, 12,
%%14, 23, 57, 52, 28, 55, 18, 22, 52, 22, 24, 23, 27, 23, 26, 28, 22, 25,
%%19, 23, 54, 54, 49, 22, 58, 21, 63, 22, 56, 23, 59, 25, 23, 27, 32, 13,
%%20, 24, 23),
%%Total = c(450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450))

\dontrun{
data(c212.trial.data)
raw = c212.BB(c212.trial.data)
}

\dontrun{
p = c212.BB.ptheta(rm)

head(p)
          B       AE    ptheta
1 Bdy-sys_1 Adv-Ev_1 0.2560500
2 Bdy-sys_2 Adv-Ev_2 0.9426417
3 Bdy-sys_2 Adv-Ev_3 0.8751500
4 Bdy-sys_2 Adv-Ev_4 0.1154917
5 Bdy-sys_2 Adv-Ev_5 0.2317417
6 Bdy-sys_3 Adv-Ev_6 1.0000000
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.ptheta}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
