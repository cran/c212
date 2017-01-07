\name{c212.ssBH}
\alias{c212.ssBH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implementation of Subset Benjamini-Hochberg for False Discover Rate control}
\description{
The Subset Benjamini-Hochberg allows for the use of subsets to allow the extension of the Benjamini-Hochberg
procedure to types of non-positively dependent regression statistics.
}
\usage{
c212.ssBH(trial.data, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{trial.data}{File or data frame containing the p-values for the hypotheses being tested. 
The data must contain the following columns: \emph{B}: the index of the groupings; \emph{p}: the p-values of the hypotheses.}
\item{alpha}{
The level for FDR control. E.g. 0.05.
}
}
\value{
The subset of hypotheses in \emph{file} or \emph{trial.data} deemed significant by the Subset Benjamini-Hochberg process.
}
\references{
Yekutieli, Daniel (2008). False discovery rate control for non-positively regression dependent test statistics. Journal of Statistical Planning and Inference, 138(2):405-415.
}
\author{
R. Carragher
}
\note{
This process is at most as powerful as the Benjamini-Hochberg procedure.
The subsets do not have to be disjoint.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
trial.data <- data.frame(B = c(1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
j = c(1, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5),
AE = c("AE1", "AE2", "AE3", "AE4", "AE5", "AE6", "AE7", "AE8", "AE9", "AE10", "AE11",
"AE12", "AE13", "AE14", "AE15", "AE16", "AE17"),
p = c(0.135005, 0.010000, 0.001000, 0.005000, 0.153501, 0.020000, 0.0013, 0.0023,
0.011, 0.023000, 0.016, 0.0109, 0.559111, 0.751986, 0.308339, 0.837154, 0.325882))

c212.ssBH(trial.data, alpha=0.05)

\dontrun{
  B j   AE      p
1 2 2  AE3 0.0010
2 2 3  AE4 0.0050
3 3 2  AE7 0.0013
4 3 3  AE8 0.0023
5 3 7 AE12 0.0109
6 3 4  AE9 0.0110
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.ssBH}
\keyword{ssBH}
\keyword{Subset Benjamini-Hochberg}
