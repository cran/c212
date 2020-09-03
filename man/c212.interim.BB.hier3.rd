\name{c212.interim.BB.hier3}
\alias{c212.interim.BB.hier3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Three-Level Hierarchical Body-system based Model for interim analysis with Point-Mass.}
\description{
Implementation of a Three-Level Hierarchical Body-system based Model for interim analysis with Point-Mass.}
\usage{
	c212.interim.BB.hier3(trial.data, sim_type = "SLICE", burnin = 20000,
	iter = 60000, nchains = 5, theta_algorithm = "MH",
	global.sim.params = data.frame(type = c("MH", "MH", "MH", "MH",
	"SLICE", "SLICE", "SLICE"),
	param = c("sigma_MH_alpha", "sigma_MH_beta", "sigma_MH_gamma",
	"sigma_MH_theta", "w_alpha", "w_beta", "w_gamma"),
	value = c(3, 3, 0.2, 0.25, 1, 1, 1), control = c(0, 0, 0, 0, 6, 6, 6),
	stringsAsFactors = FALSE),
	sim.params = NULL,
	monitor = data.frame(variable = c("theta", "gamma", "mu.gamma", "mu.theta",
	"sigma2.theta", "sigma2.gamma",
	"mu.theta.0", "mu.gamma.0", "tau2.theta.0", "tau2.gamma.0",
	"pi", "alpha.pi", "beta.pi"),
	monitor = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
	stringsAsFactors = FALSE),
	initial_values = NULL, level = 1, hyper_params = list(mu.gamma.0.0 = 0,
	tau2.gamma.0.0 = 10, mu.theta.0.0 = 0, tau2.theta.0.0 = 10,
	alpha.gamma.0.0 = 3, beta.gamma.0.0 = 1, alpha.theta.0.0 = 3,
	beta.theta.0.0 = 1, alpha.gamma = 3,
	beta.gamma = 1, alpha.theta = 3, beta.theta = 1,
	lambda.alpha = 1.0, lambda.beta = 1.0),
	global.pm.weight = 0.5,
	pm.weights = NULL,
	adapt_phase=1, memory_model = "HIGH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trial.data}{
A file or data frame containing the trial data. It must contain must contain the columns \emph{B} (body-system), \emph{AE} (adverse event), \emph{Group} (1 - control, 2 treatment), \emph{Count} (total number of events), \emph{Total} (total number of participants).
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
  \item{theta_algorithm}{
MCMC algorithm used to sample the theta variables. "MH" is the only currently supported stable algorithm.
}
  \item{sim_type}{
The type of MCMC method to use for simulating from non-standard distributions apart from theta. Allowed values are \emph{"MH"} and \emph{"SLICE"} for Metropolis_Hastings and Slice sampling respectively.
}

\item{monitor}{
A dataframe indicating which sets of variables to monitor.
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
	sigma2.theta, pi, mu.gamma.0, mu.theta.0,
	tau2.gamma.0, tau2.theta.0, alpha.pi, beta.pi)
}
The function \emph{c212.gen.initial.values} can be used to generate a template for the list which can be updated by the user if required.
The formats of the list elements are as follows:

\emph{gamma, theta}: dataframe with columns \emph{B}, \emph{AE}, \emph{chain}, \emph{value}

\emph{mu.gamma, mu.theta, sigma2.gamma, sigma2.theta, pi}: dataframe with columns \emph{B}, \emph{chain}, \emph{value}

\emph{mu.gamma.0, mu.theta.0, tau2.gamma.0, tau2.theta.0, alpha.pi, beta.pi}: array of size \emph{chain}.

}
  \item{level}{
Allowed values are 0, 1, 2. Respectively these indicate independent intervals, common body-system means across the intervals and weak relationships between the intervals.
}
  \item{hyper_params}{
The hyperparameters for the model. The default hyperparameters are those given in Berry and Berry 2004.
}

\item{global.pm.weight}{A global weighting for the proposal distribution used to sample theta.}
\item{pm.weights}{Override global.pm.weight for specific adverse events.}

  \item{adapt_phase}{
Unused parameter.
}

\item{memory_model}{
Allowed values are "HIGH" and "LOW". "HIGH" means use as much memory as possible. "LOW" means use the minimum amount of memory.
}
}
\details{
The model is fitted by a Gibbs sampler. The details of the complete conditional densities are given in Berry
and Berry (2004).
}
\value{
The output from the simulation including all the sampled values is as follows:
\preformatted{
list(id, theta_alg, sim_type, chains, nIntervals, Intervals, nBodySys,
	maxBs, maxAEs, nAE, AE, B, burnin,
	iter, monitor, mu.gamma.0, mu.theta.0, tau2.gamma.0, tau2.theta.0,
	mu.gamma, mu.theta, sigma2.gamma, sigma2.theta, pi, alpha.pi, beta.pi,
	alpha.pi_acc, beta.pi_acc, gamma, theta, gamma_acc, theta_acc)
}
where

\emph{id} - a string identifying the version of the function

\emph{theta_alg} - an string identifying the algorithm used to sample theta

\emph{sim_type} - an string identifying the sampling method used for non-standard distributions, either \emph{"MH"} or \emph{"SLICE"}

\emph{chains} - the number of chains for which the simulation was run

\emph{nIntervals} - the number of intervals in the simulation

\emph{Intervals} - an array. The intervals.

\emph{nBodySys} - the number of body-systems

\emph{maxBs} - the maximum number of body-systems in an interval

\emph{maxAEs} - the maximum number of AEs in a body-system

\emph{nAE} - an array. The number of AEs in each body-system.

\emph{AE} - an array of dimension \emph{nBodySys}, \emph{maxAEs}. The Adverse Events.

\emph{B} - an array. The body-systems.

\emph{burnin} - burnin used for the simulation.

\emph{iter} - the total number of iterations in the simulation.

\emph{monitor} - the variables being monitored. A dataframe.

\emph{mu.gamma.0} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{mu.theta.0} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{tau2.gamma.0} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{tau2.theta.0} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{mu.gamma} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}

\emph{mu.theta} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}

\emph{sigma2.gamma} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}

\emph{sigma2.theta} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}

\emph{pi} - array of samples of dimension \emph{chains}, \emph{nBodySys} \emph{iter - burnin}
\emph{alpha.pi} - array of samples of dimension \emph{chains}, \emph{iter - burnin}
\emph{beta.pi} - array of samples of dimension \emph{chains}, \emph{iter - burnin}

\emph{alpha.pi_acc} - the acceptance rate for the alpha.pi samples if a Metropolis-Hastings method is used. An array of dimension \emph{chains}, \emph{maxAEs}

\emph{beta.pi_acc} - the acceptance rate for the beta.pi samples if a Metropolis-Hastings method is used. An array of dimension \emph{chains}, \emph{maxAEs}


\emph{gamma} - array of samples of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}, \emph{iter - burnin}

\emph{theta} - array of samples of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}, \emph{iter - burnin}

\emph{gamma_acc} - the acceptance rate for the gamma samples if a Metropolis-Hastings method is used. An array of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}

\emph{theta_acc} - the acceptance rate for the theta samples. An array of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}

}
%\references{
%S. M. Berry and D. A. Berry (2004). Accounting for multiplicities in assessing drug safety: a three-
%level hierarchical mixture model.
%Biometrics, 60(2):418-26.
%}
\author{
R. Carragher
}
\note{
The function performs the simulation and returns the raw output. No checks for convergence are performed. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(c212.trial.interval.data1)
raw = c212.interim.BB.hier3(c212.trial.interval.data1, level = 1, burnin = 100, iter = 200)

\dontrun{
data(c212.trial.interval.data1)
raw = c212.interim.BB.hier3(c212.trial.interval.data1, level = 1)

raw$B
     [,1]        [,2]         [,3]         [,4]         [,5]        
[1,] "Bdy-sys_1" "Bdy-sys_10" "Bdy-sys_11" "Bdy-sys_12" "Bdy-sys_13"
[2,] "Bdy-sys_1" "Bdy-sys_10" "Bdy-sys_11" "Bdy-sys_12" "Bdy-sys_13"
[3,] "Bdy-sys_1" "Bdy-sys_10" "Bdy-sys_11" "Bdy-sys_12" "Bdy-sys_13"
[4,] "Bdy-sys_1" "Bdy-sys_10" "Bdy-sys_11" "Bdy-sys_12" "Bdy-sys_13"
[5,] "Bdy-sys_1" "Bdy-sys_10" "Bdy-sys_11" "Bdy-sys_12" "Bdy-sys_13"
[6,] "Bdy-sys_1" "Bdy-sys_10" "Bdy-sys_11" "Bdy-sys_12" "Bdy-sys_13"
     [,6]         [,7]         [,8]        [,9]        [,10]       [,11]      
[1,] "Bdy-sys_14" "Bdy-sys_15" "Bdy-sys_2" "Bdy-sys_3" "Bdy-sys_4" "Bdy-sys_5"
[2,] "Bdy-sys_14" "Bdy-sys_15" "Bdy-sys_2" "Bdy-sys_3" "Bdy-sys_4" "Bdy-sys_5"
[3,] "Bdy-sys_14" "Bdy-sys_15" "Bdy-sys_2" "Bdy-sys_3" "Bdy-sys_4" "Bdy-sys_5"
[4,] "Bdy-sys_14" "Bdy-sys_15" "Bdy-sys_2" "Bdy-sys_3" "Bdy-sys_4" "Bdy-sys_5"
[5,] "Bdy-sys_14" "Bdy-sys_15" "Bdy-sys_2" "Bdy-sys_3" "Bdy-sys_4" "Bdy-sys_5"
[6,] "Bdy-sys_14" "Bdy-sys_15" "Bdy-sys_2" "Bdy-sys_3" "Bdy-sys_4" "Bdy-sys_5"
     [,12]       [,13]       [,14]       [,15]      
[1,] "Bdy-sys_6" "Bdy-sys_7" "Bdy-sys_8" "Bdy-sys_9"
[2,] "Bdy-sys_6" "Bdy-sys_7" "Bdy-sys_8" "Bdy-sys_9"
[3,] "Bdy-sys_6" "Bdy-sys_7" "Bdy-sys_8" "Bdy-sys_9"
[4,] "Bdy-sys_6" "Bdy-sys_7" "Bdy-sys_8" "Bdy-sys_9"
[5,] "Bdy-sys_6" "Bdy-sys_7" "Bdy-sys_8" "Bdy-sys_9"
[6,] "Bdy-sys_6" "Bdy-sys_7" "Bdy-sys_8" "Bdy-sys_9"

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.interim.BB.hier3}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Point-mass} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
