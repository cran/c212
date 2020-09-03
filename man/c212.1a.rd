\name{c212.1a}
\alias{c212.1a}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implementation of the Berry and Berry Three-Level Hierarchical Model without Point-Mass.}
\description{
Implementaion of Berry and Berry model without the point-mass (Model 1a Xia et al (2011))}
\usage{
	c212.1a(trial.data, sim_type = "SLICE", burnin = 10000, iter = 40000,
	nchains = 3,
	global.sim.params = data.frame(type = c("MH", "SLICE"),
	param = c("sigma_MH", "w"),
	value = c(0.35,1), control = c(0,6), stringsAsFactors = FALSE),
	sim.params = NULL,
	initial_values = NULL,
	hyper_params = list(mu.gamma.0.0 = 0, tau2.gamma.0.0 = 10,
	mu.theta.0.0 = 0, tau2.theta.0.0 = 10, alpha.gamma.0.0 = 3,
	beta.gamma.0.0 = 1, alpha.theta.0.0 = 3, beta.theta.0.0 = 1,
	alpha.gamma = 3, beta.gamma = 1,
	alpha.theta = 3, beta.theta = 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
A file or data frame containing the trial data. It must contain must contain the columns \emph{B} (body-system), \emph{AE} (adverse event), \emph{Group} (1 - control, 2 treatment), \emph{Count} (total number of events), \emph{Total} (total number of participants in the trial
arm).
}
  \item{sim_type}{
The type of MCMC method to use for simulating from non-standard distributions. Allowed values are \emph{"MH"}
and \emph{"SLICE"} for Metropolis_Hastings and Slice sampling respectively.
}
  \item{burnin}{
The burnin period for the monte-carlo simulation. These are discarded from the returned samples.
}
  \item{iter}{
The total number of iterations for which the monte-carlo simulation is run. This includes the burnin period.
The total number of samples returned is \emph{iter - burnin}
}
  \item{nchains}{
The number of independent chains to run.
}
\item{global.sim.params}{
A data frame containing the parameters for the simulation type \emph{sim_type}. For \emph{"MH"} the parameter
is the variance of the normal distribution used to simulate the next candidate value centred on the current
value. For \emph{"SLICE"} the parameters are the estimated width of the slice and a value limiting the search for the next sample.
}
\item{sim.params}{
A dataframe containing simulation parameters which override the global simulation parameters (\emph{global.sim.params}) for particular model
parameters. \emph{sim.params} must contain the following columns: type: the simulation type ("MH" or "SLICE"); variable: the model parameter 
for which the simulation parameters are being overridden; B: the body-system (if applicable); AE: the adverse event (if applicable);
param: the simulation parameter; value: the overridden value; control: the overridden control value.

The function \emph{c212.sim.control.params} generates a template for \emph{sim.params} which can be edited by the user.
}
  \item{initial_values}{
The initial values for starting the chains. If NULL (the default) is passed the function generates the initial
values for the chains. initial_values is a list with the following format:
\preformatted{
list(gamma, theta, mu.gamma, mu.theta, sigma2.gamma,
sigma2.theta, mu.gamma.0, mu.theta.0, tau2.gamma.0,
tau2.theta.0)
}
where each element of the list is either a dataframe or array.
The function \emph{c212.gen.initial.values} can be used to generate a template for the list which can be updated by the user if required.
The formats of the list elements are as follows:

\emph{gamma, theta}: dataframe with columns \emph{B}, \emph{AE}, \emph{chain}, \emph{value}

\emph{mu.gamma, mu.theta, sigma2.gamma, sigma2.theta}: dataframe with columns \emph{B}, \emph{chain}, \emph{value}

\emph{mu.gamma.0, mu.theta.0, tau2.gamma.0, tau2.theta.0}: array of size \emph{chain}.
}
  \item{hyper_params}{
The hyperparameters for the model. The default hyperparameters are those given in Berry and Berry 2004.
}
}
\details{
The model is fitted by a Gibbs sampler. The details of the complete conditional densities are given in Berry
and Berry (2004). The posterior distributions for \emph{gamma} and \emph{theta} are sampled with either a Metropolis-Hastings step or a
slice sampler.
}
\value{
The output from the simulation including all the sampled values is as follows:
\preformatted{
list(id, sim_type, chains, nBodySys, maxAEs, nAE, AE, B, burnin,
	iter, mu.gamma.0, mu.theta.0, tau2.gamma.0, tau2.theta.0,
	mu.gamma, mu.theta, sigma2.gamma, sigma2.theta, gamma,
	theta, gamma_acc, theta_acc)
}
where

\emph{id} - a string identifying the version of the function

\emph{sim_type} - an string identifying the sampling method used for non-standard distributions, either \emph{"MH"} or \emph{"SLICE"}

\emph{chains} - the number of chains for which the simulation was run

\emph{nBodySys} - the number of body-systems

\emph{maxAEs} - the maximum number of AEs in a body-system

\emph{nAE} - an array. The number of AEs in each body-system.

\emph{AE} - an array of dimension \emph{nBodySys}, \emph{maxAEs}. The Adverse Events.

\emph{B} - an array. The body-systems.

\emph{burnin} - the burnin period for the simulation.

\emph{iter} - the total number of iterations in the simulation.

\emph{mu.gamma.0} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{mu.theta.0} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{tau2.gamma.0} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{tau2.theta.0} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{mu.gamma} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}

\emph{mu.theta} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}

\emph{sigma2.gamma} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}

\emph{sigma2.theta} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}

\emph{gamma} - array of samples of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}, \emph{iter - burnin}

\emph{theta} - array of samples of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}, \emph{iter - burnin}

\emph{gamma_acc} - the acceptance rate for the gamma samples if a Metropolis-Hastings method is used. An array of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}

\emph{theta_acc} - the acceptance rate for the theta samples if a Metropolis-Hastings method is used. An array of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}
}
\references{
S. M. Berry and D. A. Berry (2004). Accounting for multiplicities in assessing drug safety: a three-
level hierarchical mixture model.
Biometrics, 60(2):418-26.

H. Amy Xia, Haijun Ma, and Bradley P. Carlin (2011). Bayesian hierarchical modelling for
detecting safety signals in clinical trials. Journal of Biopharmaceutical Statistics, 21(5):1006–
1029.

Scott M. Berry, Bradley P. Carlin, J. Jack Lee, and Peter M¨ller (2010). Bayesian adaptive
methods for clinical trials. CRC Press.
}
\author{
R. Carragher
}
\note{
The function performs the simulation and returns the raw output. No checks for convergence are performed. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(c212.trial.data)
raw = c212.1a(c212.trial.data, burnin = 100, iter = 200)
\dontrun{
data(c212.trial.data)
raw = c212.1a(c212.trial.data)

raw$B
[1] "Bdy-sys_1" "Bdy-sys_2" "Bdy-sys_3" "Bdy-sys_4" "Bdy-sys_5" "Bdy-sys_6"
[7] "Bdy-sys_7" "Bdy-sys_8"

mean(rm$theta[2, 3,1,])
[1] 1.306362

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.1a}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
