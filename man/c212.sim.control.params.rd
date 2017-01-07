\name{c212.sim.control.params}
\alias{c212.sim.control.params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a template for the individual model parameter simulation control parameters.}
\description{
This function generates a template for overriding the global simulation parameters used by the model simulation functions (e.g. c212.BB).
}
%\details{
%}

\usage{
	c212.sim.control.params(trial.data, model = "1a")
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{trial.data}{
A file or data frame containing the trial data, either for end of trial or interim analysis.
}
  \item{model}{
Allowed values are "BB" and "1a" for point-mass and non-point-mass models respectively.
}
}


\value{
A dataframe containing the simulation parameters template.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
data(c212.trial.data)
c212.sim.control.params(c212.trial.data)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.sim.control.params}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
