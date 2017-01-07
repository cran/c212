\name{c212.print.summary.stats}
\alias{c212.print.summary.stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print the Summary Statistics of Posterior Distributions}
\description{
The function prints the variable names, the mean, the 95% HPI interval, the standard distribtion and the
MCMC standard error for the simulated sample.}
\usage{
	c212.print.summary.stats(summ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summ}{
The output from a call to c212.summary.stats
}
}
\value{
Nothing
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
summ = c212.summary.stats(rm)
c212.print.summary.stats(summ)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.print.summary.stats}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry}
