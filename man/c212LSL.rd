\name{c212.LSL}
\alias{c212.LSL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implementaion of the least-slope estimator estimator (LSL) for the proprtion of true null hypotheses.
}
\description{
The least-slope estimator estimator (LSL) is one of a number of estimators of the proportion of true null hypotheses. This implementation assumes a grouped structure for the data.
}
\usage{
c212.LSL(trial.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{trial.data}{
Data frame containing the p-values for the hypotheses being tested. The data must contain the following columns: \emph{B}: the index or name of the groupings; \emph{p}: the p-values of the hypotheses.
}
}
\value{
An estimate of the proportion of true null hypotheses.
}
\references{
Hu, J. X. and Zhao, H. and Zhou, H. H. (2010). False Discovery Rate Control With Groups. J Am Stat Assoc, 105(491):1215-1227.

Benjamini Y, Hochberg Y. (2000). On the Adaptive Control of the False Discovery Rate in Multiple Testing
With Independent Statistics. Journal of Educational and Behavioral Statistics, 25(1):60–83.
}
\author{
R. Carragher<raymond.carragher@strath.ac.uk>
}

\note{
The implementation is that described in Hu, J. X. and Zhao, H. and Zhou, H. H. (2010).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%+\examples{
%%+##---- Should be DIRECTLY executable !! ----
%%+##-- ==>  Define data, use random,
%%+##--	or do  help(data=index)  for the standard data sets.
%%+
%%+## The function is currently defined as
%%+c212.LSL(trial.data, alpha) {
%%+}
%%+}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.LSL}
