\name{c212.TST}
\alias{c212.TST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implementaion of the two-stage estimator (TST) for the proportion of true null hypotheses.
}
\description{
The two-stage estimator (TST) is one of a number of estimators of the proportion of true null hypotheses. It
uses the Benjamini-Hochberg procedure at a reduced level to make the estimate. This implementation assumes
a grouped structure for the data.
}
\usage{
c212.TST(trial.data, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{trial.data}{
Data frame containing the p-values for the hypotheses being tested. The data must contain the following columns: \emph{B}: the index or name of the groupings; \emph{p}: the p-values of the hypotheses.
}
	\item{alpha}{
The level for FDR control. E.g. 0.05.
}
}
\value{
An estimate of the proportion of true null hypotheses.
}
\references{
Hu, J. X. and Zhao, H. and Zhou, H. H. (2010). False Discovery Rate Control With Groups. J Am Stat Assoc, 105(491):1215-1227.

Y. Benjamini, A. M. Krieger, and D. Yekutieli (2006). Adaptive linear step-up procedures that
control the false discovery rate. Biometrika, 93(3):491–507.
}
\author{
R. Carragher<raymond.carragher@strath.ac.uk>
}

\note{
The implementation is that described in Hu, J. X. and Zhao, H. and Zhou, H. H. (2010).
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(c212.FDR.data)
c212.TST(c212.FDR.data)
\dontrun{
          B pi0
1 Bdy-sys_5 1.0
2 Bdy-sys_6 1.0
3 Bdy-sys_7 1.0
4 Bdy-sys_8 1.0
5 Bdy-sys_2 0.5
6 Bdy-sys_3 0.0
7 Bdy-sys_4 1.0
8 Bdy-sys_1 1.0
}
}

%%+\examples{
%%+##---- Should be DIRECTLY executable !! ----
%%+##-- ==>  Define data, use random,
%%+##--	or do  help(data=index)  for the standard data sets.
%%+
%%+## The function is currently defined as
%%+c212.TST(file, trial.data, use_file, alpha) {
%%+}
%%+}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.TST}
