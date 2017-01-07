\name{c212.pointmass.weights}
\alias{c212.pointmass.weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a template for the point-mass weightings.}
\description{
This function generate a template for weights for the proposal distribution used to sample \emph{theta} variables in models which 
use a point-mass.
}
%\details{
%}

\usage{
	c212.pointmass.weights(trial.data)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{trial.data}{
A file or data frame containing the trial data, either for end of trial or interim analysis.
}
}


\value{
A dataframe containing the weightings template for each Body-system, adverse event and, if required, interval.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
data(c212.trial.data)
c212.pointmass.weights(c212.trial.data)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.pointmass.weights}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
