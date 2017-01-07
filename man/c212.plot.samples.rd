\name{c212.plot.samples}
\alias{c212.plot.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Posterior Distribution}
\description{
This function plots a graph of the sampled posterior distribution.
}
\details{
Two graphs are displayed on the same panel. The left graph is the traceplot of the chains. The right graph is
a plot of the distribution.
}
\usage{
	c212.plot.samples(samples, title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{
An array of samples indexed by \emph{chain}.
}
  \item{title}{
The graph title.
}
}
\value{
Nothing is returned.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
data(c212.trial.data)
raw = c212.1a(c212.trial.data)
sample = raw$theta[,2,2,]
c212.plot.samples(sample, sprintf("\%s: \%s \%s", "theta", raw$B[2], raw$AE[2]))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.plot.samples}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
