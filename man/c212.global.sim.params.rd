\name{c212.global.sim.params}
\alias{c212.global.sim.params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a template for the individual model parameter simulation control parameters.}
\description{
This function generates the default global simulation parameters used by the model simulation functions (e.g. c212.BB).
}
%\details{
%}

\usage{
	c212.global.sim.params(trial.data, model = "BB", hier = 3)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{trial.data}{
A file or data frame containing the trial data, either for end of trial or interim analysis.
}
  \item{model}{
Allowed values are "BB" and "1a" for point-mass and non-point-mass models respectively.
}
  \item{hier}{
Generate parameters for a two level or three level hierarchy. Allowed values are 2 and 3 respectively.
}
}


\value{
A dataframe containing the global simulation parameters.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(c212.trial.data)
global.sim.prams <- c212.global.sim.params(c212.trial.data)
\dontrun{
data(c212.trial.data)
global.sim.prams <- c212.global.sim.params(c212.trial.data)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.global.sim.params}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
