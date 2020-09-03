\name{c212.plot.interim.data.rd}
\alias{c212.plot.interim.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Adverse Event Count Data for a Body-system by Interval}
\description{
Plot adverse event interval data for a body-system.
}
\details{
This function plots a graph of the count of adverse events which have occurred in an interval for a particular body-system by interval.
}
\usage{
	c212.plot.interim.data(trial.data, body_sys, cex = 0.8, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
A file or data frame containing the trial data. The data frame must contain the columns \emph{I_index} (interval), \emph{B} (body-system), \emph{AE} (adverse event), \emph{Group} (1 - control, 2 treatment), \emph{Count} (total number of events), \emph{Exposure} (total time of participants spent in the interval).
}
  \item{body_sys}{
The body-system for which to plot the events.
}
  \item{cex}{
Font size of the labels on the Adverse Event counts graph.
}
  \item{title}{
Main title of the graph.
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
data(c212.trial.interval.data1)
c212.plot.interim.data(c212.trial.interval.data1, "Bdy-sys_3")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.plot.interim.data}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
