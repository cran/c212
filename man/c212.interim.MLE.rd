\name{c212.interim.MLE}
\alias{c212.interim.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Poisson Maximim Likelihood Estimator}
\description{
Calculate the Poisson Maximim Likelihood Estimator for interim analysis data.
}

\usage{
	c212.interim.MLE(trial.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
A file or data frame containing the trial data. It must contain must contain the columns \emph{I_index} (interval index), \emph{B} (body-system), \emph{AE} (adverse event), \emph{Group} (1 - control, 2 treatment), \emph{Count} (total number of events), \emph{Total} (total number of participants in the trial
arm).
}
}

%\details{
%}
\value{
The maximum likelihood estimators and summary statistics.
}
%\references{
%}
\author{
R. Carragher
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
data(c212.trial.interval.data1)
raw = c212.interim.MLE(c212.trial.interval.data1)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.interim.MLE}
\keyword{Poisson MLE} % __ONLY ONE__ keyword per line
