\name{c212.BB}
\alias{c212.BB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implementation of the Berry and Berry Three-Level Hierarchical Model.}
\description{
Implementaion of Berry and Berry model (2004), also model 1b from Xia et al (2011).}
\usage{
	c212.BB(trial.data, burnin = 20000, iter = 60000, nchains = 3,
	theta_algorithm = "MH", sim_type = "SLICE",
	global.sim.params = data.frame(type = c("MH", "MH", "MH", "MH",
	"SLICE", "SLICE", "SLICE"),
	param = c("sigma_MH_alpha", "sigma_MH_beta", "sigma_MH_gamma",
	"sigma_MH_theta", "w_alpha", "w_beta", "w_gamma"),
	value = c(3, 3, 0.2, 0.2, 1, 1, 1), control = c(0, 0, 0, 0, 6, 6, 6),
	stringsAsFactors = FALSE),
	sim.params = NULL,
	initial_values = NULL,
	hyper_params = list(mu.gamma.0.0 = 0,
	tau2.gamma.0.0 = 10, mu.theta.0.0 = 0, tau2.theta.0.0 = 10,
	alpha.gamma.0.0 = 3, beta.gamma.0.0 = 1, alpha.theta.0.0 = 3,
	beta.theta.0.0 = 1, alpha.gamma = 3,
	beta.gamma = 1, alpha.theta = 3, beta.theta = 1,
	lambda.alpha = 1.0, lambda.beta = 1.0),
	global.pm.weight = 0.5,
	pm.weights = NULL,
	adapt_params = data.frame(min_w = 0.25, chains = 3, burnin = 20000,
	iter = 40000),
	adapt_phase=0)
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
The type of MCMC method to use for simulating from non-standard distributions apart from theta. Allowed values are \emph{"MH"} and \emph{"SLICE"} for Metropis_Hastings and Slice sampling respectively.
}
\item{global.sim.params}{
A data frame containing the parameters for the simuation type \emph{sim_type}. For \emph{"MH"} the parameter
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
  \item{hyper_params}{
The hyperparameters for the model. The default hyperparamters are those given in Berry and Berry 2004.
}

\item{global.pm.weight}{A global weighting for the proposal distribution used to sample theta.}
\item{pm.weights}{Override global.pm.weight for specific adverse events.}

  \item{adapt_params}{
Unused parameter.
}
  \item{adapt_phase}{
Unused parameter.
}
}
\details{
The model is fitted by a Gibbs sampler. The details of the complete conditional densities are given in Berry
and Berry (2004).
}
\value{
The output from the simulation including all the sampled values is as follows:
\preformatted{
list(id, theta_alg, sim_type, chains, nBodySys, maxAEs, nAE, AE, B,
	burnin, iter, mu.gamma.0, mu.theta.0, tau2.gamma.0,
	tau2.theta.0, mu.gamma, mu.theta, sigma2.gamma, sigma2.theta,
	pi, alpha.pi, beta.pi, alpha.pi_acc, beta.pi_acc, gamma, theta,
	gamma_acc, theta_acc, theta_zero_prop, theta_zero_acc)
}
where

\emph{id} - a string identifying the verions of the function

\emph{theta_alg} - an string identifying the algorithm used to smaple theta

\emph{sim_type} - an string identifying the samlping method used for non-standard distributions, either \emph{"MH"} or \emph{"SLICE"}

\emph{chains} - the number of chains for which the simulation was run

\emph{nBodySys} - the number of body-systems

\emph{maxAEs} - the maximum number of AEs in a body-system

\emph{nAE} - an array. The number of AEs in each body-system.

\emph{AE} - an array of dimension \emph{nBodySys}, \emph{maxAEs}. The Adverse Events.

\emph{B} - an array. The body-systems.

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

\emph{theta_zero_prop} - the number of zeros proposed in theta sampling. An array of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}
\emph{theta_zero_acc} - the acceptance rate for zeros for the theta samples. An array of dimension \emph{chains}, \emph{nBodySys}, \emph{maxAEs}
}
\references{
S. M. Berry and D. A. Berry (2004). Accounting for multiplicities in assessing drug safety: a three-
level hierarchical mixture model.
Biometrics, 60(2):418-26.

H. Amy Xia, Haijun Ma, and Bradley P. Carlin (2011). Bayesian hierarchical modeling for
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
%%trial.data <- data.frame(B = c("BS6", "BS6", "BS7", "BS7", "BS6", "BS6",
%%"BS5", "BS7", "BS6", "BS2", "BS2", "BS5", "BS4", "BS1", "BS6", "BS5",
%%"BS4", "BS4", "BS2", "BS5", "BS1", "BS6", "BS8", "BS6", "BS3", "BS8",
%%"BS8", "BS5", "BS6", "BS8", "BS5", "BS2", "BS6", "BS6", "BS3", "BS3",
%%"BS8", "BS8", "BS6", "BS4", "BS5", "BS2", "BS3", "BS4", "BS5", "BS2",
%%"BS3", "BS5", "BS6", "BS4", "BS8", "BS2", "BS6", "BS5", "BS5", "BS2",
%%"BS3", "BS8", "BS6", "BS5", "BS8", "BS7", "BS3", "BS6", "BS4", "BS6",
%%"BS6", "BS3", "BS4", "BS6", "BS3", "BS5", "BS5", "BS5", "BS6", "BS3",
%%"BS6", "BS3", "BS7", "BS3", "BS3", "BS5", "BS8", "BS8", "BS4", "BS7",
%%"BS6", "BS4", "BS3", "BS8"),
%%AE = c("AE30", "AE33", "AE38", "AE39", "AE31", "AE31", "AE25", "AE37",
%%"AE28", "AE3", "AE3", "AE24", "AE15", "AE1", "AE35", "AE18", "AE16",
%%"AE17", "AE2", "AE23", "AE1", "AE35", "AE42", "AE27", "AE7", "AE43",
%%"AE42", "AE24", "AE36", "AE40", "AE20", "AE4", "AE29", "AE27", "AE912",
%%"AE910", "AE41", "AE45", "AE26", "AE16", "AE21", "AE5", "AE912", "AE14",
%%"AE21", "AE5", "AE911", "AE23", "AE33", "AE14", "AE40", "AE4", "AE36",
%%"AE22", "AE22", "AE2", "AE9", "AE45", "AE32", "AE25", "AE44", "AE39",
%%"AE7", "AE34", "AE17", "AE34", "AE28", "AE6", "AE15", "AE29", "AE8",
%%"AE19", "AE19", "AE20", "AE26", "AE8", "AE30", "AE910", "AE38", "AE6",
%%"AE9", "AE18", "AE43", "AE44", "AE13", "AE37", "AE32", "AE13", "AE911",
%%"AE41"),
%%Group = c(2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2, 1, 2,
%%1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 2, 1,
%%2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2,
%%1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 2,
%%1, 2, 2, 1, 2, 2, 1, 1, 1, 1),
%%Count = c(20, 15, 25, 22, 24, 19, 49, 16, 12, 23, 18, 52, 19, 21, 18,
%%51, 24, 22, 48, 57, 17, 20, 24, 20, 57, 20, 31, 56, 25, 19, 46, 25, 24,
%%19, 22, 29, 20, 31, 28, 19, 66, 14, 60, 18, 51, 22, 52, 58, 22, 17, 12,
%%14, 23, 57, 52, 28, 55, 18, 22, 52, 22, 24, 23, 27, 23, 26, 28, 22, 25,
%%19, 23, 54, 54, 49, 22, 58, 21, 63, 22, 56, 23, 59, 25, 23, 27, 32, 13,
%%20, 24, 23),
%%Total = c(450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450,
%%450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450, 450))

\dontrun{
%%rm = c212.BB(trial.data, nchains = 3)
data(c212.trial.data)
raw = c212.BB(c212.trial.data)

raw$B
[1] "Bdy-sys_1" "Bdy-sys_2" "Bdy-sys_3" "Bdy-sys_4" "Bdy-sys_5" "Bdy-sys_6"
[7] "Bdy-sys_7" "Bdy-sys_8"

mean(raw$theta[2, 1,1,])
[1] 0.1088401

median(raw$theta[2, 1,1,])
[1] 0

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{c212.BB}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Point-mass} % __ONLY ONE__ keyword per line
\keyword{Berry and Berry} % __ONLY ONE__ keyword per line
