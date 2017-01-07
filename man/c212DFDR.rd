\name{c212.DFDR}
\alias{c212.DFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implementation of the Double False Discovery Rate for controlling the False Discovery Rate.}
\description{
The Double False Discovery Rate is designed to take advantage of possible groupings which may exist within 
sets of hypotheses. It applies the BH-procedure twice. Once at the group level, to identify sets of hypotheses
which may contain significant hypotheses. It then groups these hypotheses together to form a single family and
applies the BH-procedure again to declare hypotheses significant.}
\usage{
c212.DFDR(trial.data, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
File or data frame containing the p-values for the hypotheses being tested. 
The data must contain the following columns: \emph{B}: the index or name of the groupings; \emph{p}: the p-values of the hypotheses.
}
\item{alpha}{
The level for FDR control. E.g. 0.05.
}
}
\value{
The subset of hypotheses in \emph{file} or \emph{trial.data} deemed significant by the Double False Discovery Rate process.
}
\references{
Mehrotra, D. V. and Adewale, A. J. (2012). Flagging clinical adverse experiences: reducing false discoveries without materially compromising power for detecting true signals. Stat Med, 31(18):1918-30.
}
\author{
R. Carragher
}

\examples{
trial.data <- data.frame(B = c(1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
AE = c("AE1", "AE2", "AE3", "AE4", "AE5", "AE6", "AE7", "AE8", "AE9", "AE10", "AE11",
"AE12", "AE13", "AE14", "AE15", "AE16", "AE17"),
p = c(0.135005, 0.010000, 0.001000, 0.005000, 0.153501, 0.020000, 0.0013, 0.0023,
0.011, 0.023000, 0.016, 0.0109, 0.559111, 0.751986, 0.308339, 0.837154, 0.325882))


c212.DFDR(trial.data, 0.05)


\dontrun{
   B j   AE      p
1  2 2  AE3 0.0010
2  3 2  AE7 0.0013
3  3 3  AE8 0.0023
4  2 3  AE4 0.0050
5  2 1  AE2 0.0100
6  3 7 AE12 0.0109
7  3 4  AE9 0.0110
8  3 6 AE11 0.0160
9  3 1  AE6 0.0200
10 3 5 AE10 0.0230
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.DFDR}
\keyword{DFDR}
\keyword{Double False Discovery Rate}
