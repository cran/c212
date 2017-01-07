\name{c212.plot.eot.data}
\alias{c212.plot.eot.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Adverse Event Incidence Data}
\description{
This function plots a graph of the total adverse event incidence counts by body-system and by individual
adverse event for end of trial data. }
\details{
Two graphs are displayed on the same panel. The top graph is of AE Incidence Counts by Body-System. The
lower graphs is of the Individual AE Incidence Counts.
}
\usage{
	c212.plot.eot.data(trial.data, legend = TRUE, interactive = FALSE,
		cex = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
A file or data frame containing the trial data. The data frame must contain the columns \emph{B} (body-system), \emph{AE} (adverse event), \emph{Group} (1 - control, 2 treatment), \emph{Count} (total number of events), \emph{Total} (total number of participants).
}
  \item{legend}{
Boolean. If TRUE print a legend.
}
  \item{interactive}{
Boolean. If TRUE allow the user to identify individual adverse events on the individual adverse events graph.
}
  \item{cex}{
Font size of the labels on the Adverse Event counts graph.
}
}
\value{
Nothing is returned.
}
\author{
R. Carragher
}
\note{
The legend may obscure a portion of the graph. In this case the legend may be supressed by choosing \emph{legend = FALSE} when calling the function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
data(c212.trial.data)
c212.plot.eot.data(c212.trial.data)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.plot.eot.data}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
