\name{c212.fisher.test}
\alias{c212.fisher.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fisher Exact Test}
\description{
Perform a Fisher exact test on clinical trial data.
}
\details{
Perform a Fisher exact test on clinical trial data.
}
\usage{
	c212.fisher.test(trial.data, alternative = "two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
A file or data frame containing the trial data. The data frame must contain the columns \emph{B} (body-system), \emph{AE} (adverse event), \emph{Group} (1 - control, 2 treatment), \emph{Count} (total number of events), \emph{Total} (total number of participants).
}
  \item{alternative}{
Alternative hypothesis may be "two.sided", "greater" or "less". The default is "two.sided".
}
}

\value{
Dataframe containing the results of the test. A copy of the input dataframe with an additional column \emph{p} containing the
p-value frmo the test.
}
\author{
R. Carragher
}

\note{
Wrapper for the R function 'fisher.test' in package 'stats'.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

\dontrun{
data(c212.trial.data)
f = c212.fisher.test(c212.trial.data)
head(f)
          B j       AE            p
1 Bdy-sys_1 1 Adv-Ev_1 2.892876e-01
2 Bdy-sys_2 1 Adv-Ev_2 5.333164e-03
3 Bdy-sys_2 2 Adv-Ev_3 1.601311e-02
4 Bdy-sys_2 3 Adv-Ev_4 6.502108e-01
5 Bdy-sys_2 4 Adv-Ev_5 7.437946e-01
6 Bdy-sys_3 1 Adv-Ev_6 3.746249e-08

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.fisher.test}
\keyword{Fisher} % __ONLY ONE__ keyword per line
\keyword{Exact} % __ONLY ONE__ keyword per line
\keyword{Test} % __ONLY ONE__ keyword per line
