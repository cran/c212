\name{c212.err.cntrl}
\alias{c212.err.cntrl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implementaion of Group Bonferroni-Hochberg procedure for control of the False Discovery Rate}
\description{
Common interface to the error controlling methods: Unadjutsed hypothesis testing (NOADJ), Bonferroni correction (BONF), Benjamini-Hochberg procedure (BH),
Group Benjamini-Hochberg (GBH),
Double False Discover Rate (DFDR), subset Benjamini-Hochberg (ssBH).
}
\usage{
c212.err.cntrl(trial.data, alpha = 0.05, method = "NOADJ",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
File or data frame containing the p-values for the hypotheses being tested. 
The data must contain the following columns: \emph{B}: the index or name of the groupings; \emph{p}: the p-values of the hypotheses.
}
\item{alpha}{
The level for error control. E.g. 0.05.
}
\item{method}{
The error control procedure to be applied:
"NOAD" - unadjusted testing,
"BONF" - Bonferroni correction
"BH" -  Benjamini-Hochberg procedure
"GBH" - Group Benjamini-Hochberg
"DFDR" - Double False Discover Rate
"ssBH" - subset Benjamini-Hochberg.
}
\item{...}{
Additional optional parameter for the GBH method: \emph{pi0}. 
}
}
\value{
The subset of hypotheses in \emph{file} or \emph{trial.data} deemed significant by the Group Benjamini-Hochberg process.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
trial.data <- data.frame(B = c(1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
j = c(1, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5),
AE = c("AE1", "AE2", "AE3", "AE4", "AE5", "AE6", "AE7", "AE8", "AE9", "AE10", "AE11",
"AE12", "AE13", "AE14", "AE15", "AE16", "AE17"),
p = c(0.135005, 0.010000, 0.001000, 0.005000, 0.153501, 0.020000, 0.0013, 0.0023,
0.011, 0.023000, 0.016, 0.0109, 0.559111, 0.751986, 0.308339, 0.837154, 0.325882))


c212.err.cntrl(trial.data = trial.data, alpha = 0.05, method = "GBH")

\dontrun{
   B j   AE        p   p_weighted
1  3 1  AE6 0.020000 0.0000000000
2  3 2  AE7 0.001300 0.0000000000
3  3 3  AE8 0.002300 0.0000000000
4  3 4  AE9 0.011000 0.0000000000
5  3 5 AE10 0.023000 0.0000000000
6  3 6 AE11 0.016000 0.0000000000
7  3 7 AE12 0.010900 0.0000000000
8  2 2  AE3 0.001000 0.0003333333
9  2 3  AE4 0.005000 0.0016666667
10 2 1  AE2 0.010000 0.0033333333
11 2 4  AE5 0.153501 0.0511670000
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.err.cntrl}
\keyword{BONF}
\keyword{BH}
\keyword{GBH}
\keyword{DFDR}
\keyword{ssBH}
