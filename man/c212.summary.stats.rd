\name{c212.summary.stats}
\alias{c212.summary.stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary Statistics for the Posterior Distributions in the model.}
\description{
Returns the Summary Statistics for the Posterior Distributions in the model.
}
\details{
The function reports the mean, upper and lower bounds of the 95% HPI (higest probabily interval), the 
standard deviation and MCMC standard error.
}
\usage{
	c212.summary.stats(raw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
The output from a model simulation (e.g. c212.BB).
}
}
\value{
Returns a list of the summary statistics for each sampled variable.
Each element of the list is a data.frame containing at least the columns \emph{mean}, \emph{hpi_lower},
\emph{hpi_upper}, \emph{SD} and \emph{SE}.
For the simulation return by \emph{c212.1a} the output is as follows:
\preformatted{
list(theta.summary, gamma.summary, mu.gamma.summary,
								mu.theta.summary = mu.theta_summ,
								sigma2.gamma.summary, sigma2.theta.summary,
								mu.gamma.0.summary, mu.theta.0.summary,
								tau2.gamma.0.summary, tau2.theta.0.summary)
}
Additional columns which may be used to indentify the individual variables are \emph{B}, the body-system, and
\emph{AE}, the Adverse Event and \emph{interval}.
}
\note{
The MCMC error is found using the 'coda' summary function.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
data(c212.trial.data)
raw = c212.BB(c212.trial.data)
summ = c212.summary.stats(raw)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.summary.stats}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
