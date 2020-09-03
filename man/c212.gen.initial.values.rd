\name{c212.gen.initial.values}
\alias{c212.gen.initial.values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a template simulation initial values.}
\description{
This function generates a template for the initial values to be used to start the simulation. They can be updated by the caller and passed to the simulation function.
}
%\details{
%}

\usage{
	c212.gen.initial.values(trial.data, nchains = 3,
						model = "1a", hier = 3, level = 0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{trial.data}{
A file or data frame containing the trial data, either for end of trial or interim analysis.
}
\item{nchains}{
The number of chains in the simulation.
}
\item{model}{
The model type: "BB" for point-mass models, "1a" for non-point-mass models.
}
\item{hier}{
Allowed values are 2, 3 indicating two or three level hierarchies.
}
\item{level}{
The dependency level in the model: 0, 1, 2.
}
}

\value{
A dataframe containing the template of initial values.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(c212.trial.data)
init.vals <- c212.gen.initial.values(c212.trial.data)
print(init.vals$mu.gamma.0)
\dontrun{
data(c212.trial.data)
init.vals <- c212.gen.initial.values(c212.trial.data)
print(init.vals$mu.gamma.0)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.gen.initial.values}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
