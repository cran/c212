\name{c212.BONF}
\alias{c212.BONF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implementation of Bonferroni correction for error control}
\description{
The Bonferroni correction controls the Familywise Error Rate.}
\usage{
c212.BONF(trial.data, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{trial.data}{
File or data frame containing the p-values for the hypotheses being tested. The data must include a column called \emph{p} which contains
the p-values of the hypotheses.
}
\item{alpha}{
The value for error control, e.g. 0.05.
}
}
\value{
The subset of hypotheses in \emph{file} or \emph{trial.data} deemed significant.
}
\references{
Matthews, John N. S. (2006) Introduction to Randomized Controlled Clinical Trials, Second Edition. Chapman & Hall/CRC Texts in Statistical Science.
}
\author{
R. Carragher
}

\note{
No check is made for duplicate rows in the input file or data frame.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
trial.data <- data.frame(B = c(1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
j = c(1, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5),
AE = c("AE1", "AE2", "AE3", "AE4", "AE5", "AE6", "AE7", "AE8", "AE9", "AE10", "AE11",
"AE12", "AE13", "AE14", "AE15", "AE16", "AE17"),
p = c(0.135005, 0.010000, 0.001000, 0.005000, 0.153501, 0.020000, 0.0013, 0.0023, 0.011,
0.023000, 0.016, 0.0109, 0.559111, 0.751986, 0.308339, 0.837154, 0.325882))


c212.BONF(trial.data, 0.05)


\dontrun{
  B j  AE      p
1 2 2 AE3 0.0010
2 3 2 AE7 0.0013
3 3 3 AE8 0.0023
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.BONF}
\keyword{Bonferroni}
\keyword{Familywise error rate}
\keyword{FWER}
