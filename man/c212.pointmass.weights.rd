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
data(c212.trial.data)
pmw <- c212.pointmass.weights(c212.trial.data)
head(pmw)
\dontrun{
data(c212.trial.data)
pmw <- c212.pointmass.weights(c212.trial.data)
head(pmw)
          B        AE weight_pm
1 Bdy-sys_2  Adv-Ev_5       0.5
2 Bdy-sys_5 Adv-Ev_24       0.5
3 Bdy-sys_6 Adv-Ev_31       0.5
4 Bdy-sys_8 Adv-Ev_42       0.5
5 Bdy-sys_7 Adv-Ev_39       0.5
6 Bdy-sys_6 Adv-Ev_34       0.5
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.pointmass.weights}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
