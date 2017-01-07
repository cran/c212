\name{c212.convergence.diag}
\alias{c212.convergence.diag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convergence Diagnostics of the Simulation}
\description{
The function applies either Gelman-Rubin or the Geweke diagnostic to the raw output of model simulation (e.g. c212.BB).
It returns the convergence diagnostics and, if applicable, the acceptance rates for the sampling distributions.
}
\details{
parameter is time consuming. This function applies one of two convergence diagnostics to the raw output of a model simulation
in order to allow convergence to be assessed. The two diagnotics are:

i) Gelman-Rubin diagnostic - used when there is more than one chain. A value close to 1 is consistent with
an MCMC simulation which has converged. The `coda' diagnostic returns a point estimate and upper confidence
limits.

ii) Geweke diagnostic - used when there is a single chain. A Z-score which is consistent with a standard normal
distribution is expected from an MCMC simulation which has converged.

The raw sample data is converted to `coda' format (mcmc objects) and the `coda' methods gelman.diag and
geweke.diag are used to perform the checks.
}
\usage{
	c212.convergence.diag(raw, debug_diagnostic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
The output from a model simulation.
}
  \item{debug_diagnostic}{
Unused parameter.
}
}
\value{
Returns a list of the diagnostics for each sampled variable. Each individual element of the list is a
data.frame containing at least the columns \emph{type}, which is the type of diagnostic
(`Gelman-Rubin' or `Geweke'), \emph{stat}, which is the value of the dignostic, and \emph{upper_ci} which is
the upper confidence interval for the Gelman-Rubin diagnostic. For the Geweke diagnostic \emph{upper_ci}
contains the value NA. Depending on the simulation performed the return from \emph{c212.convergence.diag} will contain different
variables. The return for a simulation from \emph{c212.1a} is as follows:
\preformatted{
list(sim_type, gamma.conv.diag, theta.conv.diag, mu.gamma.conv.diag,
                       mu.theta.conv.diag, sigma2.gamma.conv.diag,
                       sigma2.theta.conv.diag, mu.gamma.0.conv.diag,
                       mu.theta.0.conv.diag, tau2.gamma.0.conv.diag,
                       tau2.theta.0.conv.diag)
}
Additional columns which may be used to indentify the individual samples are \emph{B}, the body-system, and
\emph{AE}, the Adverse Event and \emph{interval}.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

\dontrun{
data(c212.trial.data)
raw = c212.BB(c212.trial.data)
conv = c212.convergence.diag(raw)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.convergence.diag}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Gelman-Rubin}
\keyword{Berry and Berry}
