\name{c212.bin.test}
\alias{c212.bin.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Raw Adverse Event Incidence Data}
\description{
Test the hypothesis that the proportions in two groups are the same.
adverse event. }
\details{
Test the hypothesis that the proportions in two groups are the same.
}
\usage{
	c212.bin.test(trial.data, alternative = "two.sided", correct = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
A file or data frame containing the trial data. The data frame must contain the columns \emph{B} (body-system), \emph{AE} (adverse event), \emph{Group} (1 - control, 2 treatment), \emph{Count} (total number of events), \emph{Total} (total number of participants).
}
  \item{alternative}{
Alternative hypothesis may be "two.sided", "greater" or "less". The default is "two.sided".
}
  \item{correct}{
Apply a continuity correction.
}
}

\value{
Dataframe containing the results of the test. A copy of the input dataframe with an additional column \emph{p} containing the
p-value from the test.
}
\author{
R. Carragher
}

\note{
Wrapper for the R function 'prop.test' in package 'stats'.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(c212.trial.data)
pr = c212.bin.test(c212.trial.data)
head(pr)

\dontrun{
          B j       AE            p
1 Bdy-sys_1 1 Adv-Ev_1 2.893605e-01
2 Bdy-sys_2 1 Adv-Ev_2 5.711463e-03
3 Bdy-sys_2 2 Adv-Ev_3 1.655715e-02
4 Bdy-sys_2 3 Adv-Ev_4 6.497695e-01
5 Bdy-sys_2 4 Adv-Ev_5 7.433433e-01
6 Bdy-sys_3 1 Adv-Ev_6 8.419469e-08

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.bin.test}
\keyword{Proportion} % __ONLY ONE__ keyword per line
\keyword{Hypothesis} % __ONLY ONE__ keyword per line
\keyword{Test} % __ONLY ONE__ keyword per line
