#include<cstdio>
#include<cstdlib>
#include<cmath>

#include "c212_Rdefines.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>

#include <Rdefines.h>
#include <Rinternals.h>

#include "c2121a.h"

using namespace std;

//static const char *rcsId = "$Id: c2121a.cpp,v 1.19 2020/09/04 12:49:13 clb13102 Exp clb13102 $";

const char* c2121a::sColType = "type";
const char* c2121a::sColVariable = "variable";
const char* c2121a::sColParam = "param";
const char* c2121a::sColValue = "value";
const char* c2121a::sColControl = "control";
const char* c2121a::sColB = "B";
const char* c2121a::sColj = "j";

const char* c2121a::sParam_w = "w";
const char* c2121a::sParam_sigma_MH = "sigma_MH";
const char* c2121a::sVariable_gamma = "gamma";
const char* c2121a::sVariable_theta = "theta";

c2121a::c2121a()
{
	//Rprintf("c2121a::c2121a: Default constructor\n");
	gChains = 0;
	gBurnin = 0;
	gIter = 0;
	sim_type = NULL;
	gNAE = NULL;
	gNumBodySys = 0;
	gMaxAEs = 0;
	gSim_Param = 0.0;
	gSim_Param_cntrl = 0.0;

	mu_theta_0_0 = 0.0;
	mu_gamma_0_0 = 0.0;
	tau2_theta_0_0 = 0.0;
	tau2_gamma_0_0 = 0.0;
	alpha_gamma_0_0 = 0.0;
	beta_gamma_0_0 = 0.0;
	alpha_theta_0_0 = 0.0;
	beta_theta_0_0 = 0.0;
	alpha_gamma = 0.0;
	beta_gamma = 0.0;
	alpha_theta = 0.0;
	beta_theta = 0.0;

	mu_theta_0 = NULL;
	mu_gamma_0 = NULL;
	tau2_theta_0 = NULL;
	tau2_gamma_0 = NULL;

	mu_theta = NULL;
	mu_gamma = NULL;
	sigma2_theta = NULL;
	sigma2_gamma = NULL;

	gTheta = NULL;
	gGamma = NULL;
	gTheta_acc = NULL;
	gGamma_acc = NULL;

	x = NULL;
	y = NULL;
	NC = NULL;
	NT = NULL;

	gTheta_samples = NULL;
	gGamma_samples = NULL;
	mu_theta_samples = NULL;
	mu_gamma_samples = NULL;
	sigma2_theta_samples = NULL;
	sigma2_gamma_samples = NULL;
	mu_theta_0_samples = NULL;
	mu_gamma_0_samples = NULL;
	tau2_theta_0_samples = NULL;
	tau2_gamma_0_samples = NULL;

	gW_gamma = NULL;
	gW_theta = NULL;
	gW_gamma_control = NULL;
	gW_theta_control = NULL;
	gSigma_MH_gamma = NULL;
	gSigma_MH_theta = NULL;
}

c2121a::c2121a(SEXP pChains, SEXP pBurnin, SEXP pIter, SEXP pNumBodySys, SEXP pMaxAEs,
                    SEXP pNAE, SEXP pSim_Type, SEXP pGlobal_Sim_Param,
					SEXP pGlobal_Sim_Param_Cntrl,
					SEXP pSim_Param,
                    SEXP pX, SEXP pY,
                    SEXP pNC, SEXP pNT,
                    SEXP ptheta, SEXP pgamma,
                    SEXP pmu_gamma_0_0,
                    SEXP ptau2_gamma_0_0,
                    SEXP pmu_theta_0_0,
                    SEXP ptau2_theta_0_0,
                    SEXP palpha_gamma_0_0,
                    SEXP pbeta_gamma_0_0,
                    SEXP palpha_theta_0_0,
                    SEXP pbeta_theta_0_0,
                    SEXP palpha_gamma,
                    SEXP pbeta_gamma,
                    SEXP palpha_theta,
                    SEXP pbeta_theta,
                    SEXP pmu_gamma_0,
                    SEXP ptau2_gamma_0,
                    SEXP pmu_theta_0,
                    SEXP ptau2_theta_0,
                    SEXP pmu_gamma,
                    SEXP pmu_theta,
                    SEXP psigma2_gamma,
                    SEXP psigma2_theta)
{
	//Rprintf("c2121a::c2121a: constructor\n");
	gChains = 0;
	gBurnin = 0;
	gIter = 0;
	sim_type = NULL;
	gNAE = NULL;
	gNumBodySys = 0;
	gMaxAEs = 0;
	gSim_Param = 0.0;
	gSim_Param_cntrl = 0.0;

	mu_theta_0_0 = 0.0;
	mu_gamma_0_0 = 0.0;
	tau2_theta_0_0 = 0.0;
	tau2_gamma_0_0 = 0.0;
	alpha_gamma_0_0 = 0.0;
	beta_gamma_0_0 = 0.0;
	alpha_theta_0_0 = 0.0;
	beta_theta_0_0 = 0.0;
	alpha_gamma = 0.0;
	beta_gamma = 0.0;
	alpha_theta = 0.0;
	beta_theta = 0.0;

	mu_theta_0 = NULL;
	mu_gamma_0 = NULL;
	tau2_theta_0 = NULL;
	tau2_gamma_0 = NULL;

	mu_theta = NULL;
	mu_gamma = NULL;
	sigma2_theta = NULL;
	sigma2_gamma = NULL;

	gTheta = NULL;
	gGamma = NULL;
	gTheta_acc = NULL;
	gGamma_acc = NULL;

	x = NULL;
	y = NULL;
	NC = NULL;
	NT = NULL;

	gTheta_samples = NULL;
	gGamma_samples = NULL;
	mu_theta_samples = NULL;
	mu_gamma_samples = NULL;
	sigma2_theta_samples = NULL;
	sigma2_gamma_samples = NULL;
	mu_theta_0_samples = NULL;
	mu_gamma_0_samples = NULL;
	tau2_theta_0_samples = NULL;
	tau2_gamma_0_samples = NULL;

	gW_gamma = NULL;
	gW_theta = NULL;
	gW_gamma_control = NULL;
	gW_theta_control = NULL;
	gSigma_MH_gamma = NULL;
	gSigma_MH_theta = NULL;

	init(pChains, pBurnin, pIter, pNumBodySys, pMaxAEs, pNAE, pSim_Type,
			pGlobal_Sim_Param, pGlobal_Sim_Param_Cntrl, pSim_Param, pX, pY, pNC, pNT,
			ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0, pmu_theta_0_0,
			ptau2_theta_0_0,
			palpha_gamma_0_0, pbeta_gamma_0_0, palpha_theta_0_0, pbeta_theta_0_0,
			palpha_gamma, pbeta_gamma, palpha_theta, pbeta_theta, pmu_gamma_0,
			ptau2_gamma_0,
            pmu_theta_0, ptau2_theta_0, pmu_gamma, pmu_theta, psigma2_gamma,
			psigma2_theta);
}

c2121a::~c2121a()
{
	//Rprintf("c2121a::c2121a - destructor\n");
	release();
}

void c2121a::gibbs_sampler()
{
	if (strcmp(sim_type, "MH") == 0) {
		simulate_MH();
	}
	else {
		simulate_SLICE();
	}

	return;
}

void c2121a::simulate_MH()
{
	int i = 0, c = 0;

	for (c = 0; c < gChains; c++) {
		for (i = 0; i < gIter; i++) {
#ifndef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			sample_mu_gamma_0(c, gBurnin, i);
			sample_mu_theta_0(c, gBurnin, i);
			sample_tau2_gamma_0(c, gBurnin, i);
			sample_tau2_theta_0(c, gBurnin, i);
			sample_mu_gamma(c, gBurnin, i);
			sample_mu_theta(c, gBurnin, i);
			sample_sigma2_gamma(c, gBurnin, i);
			sample_sigma2_theta(c, gBurnin, i);
			sample_gamma_MH(c, gBurnin, i);
			sample_theta_MH(c, gBurnin, i);
#ifndef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			if ((i + 1)%1000 == 0) {
				Rprintf("%d iterations...\n", i + 1);
			}
		}
		Rprintf("MCMC chain fitting complete.\n");
	}
}

void c2121a::simulate_SLICE()
{
	int i = 0, c = 0;

	for (c = 0; c < gChains; c++) {
		for (i = 0; i < gIter; i++) {
#ifndef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			sample_mu_gamma_0(c, gBurnin, i);
			sample_mu_theta_0(c, gBurnin, i);
			sample_tau2_gamma_0(c, gBurnin, i);
			sample_tau2_theta_0(c, gBurnin, i);
			sample_mu_gamma(c, gBurnin, i);
			sample_mu_theta(c, gBurnin, i);
			sample_sigma2_gamma(c, gBurnin, i);
			sample_sigma2_theta(c, gBurnin, i);
			sample_gamma_SLICE(c, gBurnin, i);
			sample_theta_SLICE(c, gBurnin, i);
#ifndef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			if ((i + 1)%1000 == 0) {
				Rprintf("%d iterations...\n", i + 1);
			}
		}
		Rprintf("MCMC chain fitting complete.\n");
	}
}

void c2121a::sample_mu_gamma_0(int c, int burnin, int iter)
{
	double denom = tau2_gamma_0[c] + tau2_gamma_0_0 * ((double)gNumBodySys);

	double mu_gamma_tot = 0.0;
	int i = 0;

	for (i = 0; i < gNumBodySys; i++) {
		mu_gamma_tot += mu_gamma[c][i];
	}

	double mean = (tau2_gamma_0[c] * mu_gamma_0_0 + tau2_gamma_0_0 * mu_gamma_tot)/denom;
	double var = (tau2_gamma_0[c] * tau2_gamma_0_0) / denom;

	double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
	GetRNGstate();
#endif
	mu_gamma_0[c] = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
	PutRNGstate();
#endif

	if (iter >= burnin) {
		mu_gamma_0_samples[c][iter - burnin] = mu_gamma_0[c];
	}
}

void c2121a::sample_mu_theta_0(int c, int burnin, int iter)
{

	double denom = tau2_theta_0[c] + tau2_theta_0_0 * ((double)gNumBodySys);

	double mu_theta_tot = 0.0;
	int i = 0;

	for (i = 0; i < gNumBodySys; i++) {
		mu_theta_tot += mu_theta[c][i];
	}

	double mean = (tau2_theta_0[c] * mu_theta_0_0 + tau2_theta_0_0 * mu_theta_tot)/denom;
	double var = (tau2_theta_0[c] * tau2_theta_0_0) / denom;

	double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
	GetRNGstate();
#endif
	mu_theta_0[c] = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
	PutRNGstate();
#endif

	if (iter >= burnin) {
		mu_theta_0_samples[c][iter - burnin] = mu_theta_0[c];
	}
}

void c2121a::sample_tau2_gamma_0(int c, int burnin, int iter)
{
	double s = alpha_gamma_0_0 + ((double)gNumBodySys)/2.0;
	double r = 0.0;
	double isum = 0.0;

	int i = 0;
	for (i = 0; i < gNumBodySys; i++) {
		isum += (pow((mu_gamma[c][i] - mu_gamma_0[c]), 2.0));
	}
	r = beta_gamma_0_0 + 0.5*isum;

	// In the C API the gamma distribution is defined to be shape/scale rather than shape/rate
#ifdef INDIVIDUAL_RNG
	GetRNGstate();
#endif
	double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
	PutRNGstate();
#endif
	
	tau2_gamma_0[c] = 1/cand;

	if (iter >= burnin) {
		tau2_gamma_0_samples[c][iter - burnin] = tau2_gamma_0[c];
	}
}

void c2121a::sample_tau2_theta_0(int c, int burnin, int iter)
{

	double s = alpha_theta_0_0 + ((double)gNumBodySys)/2.0;
	double isum = 0.0;
	double r = 0.0;

	int i = 0;
	for (i = 0; i < gNumBodySys; i++) {
		isum += (pow((mu_theta[c][i] - mu_theta_0[c]), 2.0));
	}

	r = beta_theta_0_0 + 0.5 * isum;

#ifdef INDIVIDUAL_RNG
	GetRNGstate();
#endif
	double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
	PutRNGstate();
#endif
	
	tau2_theta_0[c] = 1/cand;

	if (iter >= burnin) {
		tau2_theta_0_samples[c][iter - burnin] = tau2_theta_0[c];
	}
}

void c2121a::sample_mu_gamma(int c, int burnin, int iter)
{
	int b = 0;

	for (b = 0; b < gNumBodySys; b++) {
		double denom = sigma2_gamma[c][b] + ((double)gNAE[b])*tau2_gamma_0[c];

		double t = 0.0;
		int j = 0;
		for (j = 0; j < gNAE[b]; j++) {
			t += gGamma[c][b][j];
		}

		double mean = (sigma2_gamma[c][b] * mu_gamma_0[c] + tau2_gamma_0[c] * t)/denom;

		double var = (sigma2_gamma[c][b]*tau2_gamma_0[c])/denom;

		double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		double cand = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		mu_gamma[c][b] = cand;

		if (iter >= burnin) {
			mu_gamma_samples[c][b][iter - burnin] = mu_gamma[c][b];
		}
	}
}

void c2121a::sample_mu_theta(int c, int burnin, int iter)
{
	int b = 0;

	for (b = 0; b < gNumBodySys; b++) {
		double denom = sigma2_theta[c][b] + ((double)gNAE[b])*tau2_theta_0[c];

		double t = 0.0;
		int j = 0;
		for (j = 0; j < gNAE[b]; j++) {
			t += gTheta[c][b][j];
		}

		double mean = (sigma2_theta[c][b] * mu_theta_0[c] + tau2_theta_0[c] * t)/denom;

		double var = (sigma2_theta[c][b]*tau2_theta_0[c])/denom;

		double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		double cand = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		mu_theta[c][b] = cand;

		if (iter >= burnin) {
			mu_theta_samples[c][b][iter - burnin] = mu_theta[c][b];
		}

	}
}

void c2121a::sample_sigma2_gamma(int c, int burnin, int iter)
{
	int b = 0;
	for (b = 0; b < gNumBodySys; b++) {

		double s = alpha_gamma + ((double)gNAE[b])/2.0;

		double t = 0.0;
		int j = 0;
		for (j = 0; j < gNAE[b]; j++) {
			t += (pow(gGamma[c][b][j] - mu_gamma[c][b],2.0));
		}

		double r = beta_gamma + t/2.0;

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		sigma2_gamma[c][b] = 1/cand;

		if (iter >= burnin) {
			sigma2_gamma_samples[c][b][iter - burnin] = sigma2_gamma[c][b];
		}
	}
}

void c2121a::sample_sigma2_theta(int c, int burnin, int iter)
{

	int b = 0;

	for (b = 0; b < gNumBodySys; b++) {

		double s = alpha_theta + ((double)gNAE[b])/2.0;


		double t = 0;
		int j = 0;
		for (j = 0; j < gNAE[b]; j++) {
			t += (pow((gTheta[c][b][j] - mu_theta[c][b]),2.0));
		}

		double r = beta_theta + t/2.0;

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		sigma2_theta[c][b] = 1/cand;

		if (iter >= burnin) {
			sigma2_theta_samples[c][b][iter - burnin] = sigma2_theta[c][b];
		}
	}

}

double c2121a::log_f_gamma(int c, int b, int j, double gamm)
{
	// Drop the constant term from N(mu.gamma[b], sigma^2[b])
	double q = ((-1.0/2.0)*(pow((gamm - mu_gamma[c][b]),2.0))/sigma2_gamma[c][b]);

	double f = q + ((double)x[b][j]) * gamm  - ((double)NC[b][j]) *(log(1 + exp(gamm))) + ((double)y[b][j]) * (gamm + gTheta[c][b][j]) - ((double)NT[b][j]) *(log(1.0 + exp(gamm + gTheta[c][b][j])));

	return(f);
}

void c2121a::sample_gamma_MH(int c, int burnin, int iter)
{
	int b = 0, j = 0;

	for (b = 0; b < gNumBodySys; b++) {
		for (j = 0; j < gNAE[b]; j++) {


#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double cand = rnorm(gGamma[c][b][j], gSigma_MH_gamma[b][j]);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			double f1 = log_f_gamma(c, b, j , cand);
			double f2 = log_f_gamma(c, b, j , gGamma[c][b][j]);

			double ratio = exp(f1 - f2);

			ratio = cMIN(ratio, 1.0);

			if (u <= ratio) {
				gGamma[c][b][j] = cand;
				gGamma_acc[c][b][j] = gGamma_acc[c][b][j] + 1;
			}

			if (iter >= burnin) {
				gGamma_samples[c][b][j][iter - burnin] = gGamma[c][b][j];
			}

		}
	}
}

void c2121a::sample_gamma_SLICE(int c, int burnin, int iter)
{
	int K = 0, J = 0;
	int b = 0, j = 0;

	for (b = 0; b < gNumBodySys; b++) {
		for (j = 0; j < gNAE[b]; j++) {

			int m = (int)gW_gamma_control[b][j];

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			J = floor(runif(0,m));
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif
			K = (m-1) - J;

			double l = 0.0, r = 0.0;
			double g = log_f_gamma(c, b, j, gGamma[c][b][j]);
			double logy = 0.0;
#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double e = rexp(1);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif
			logy = g - e;

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double u = runif(0, gW_gamma[b][j]);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			l = gGamma[c][b][j] - u;
			r = gGamma[c][b][j] + (gW_gamma[b][j] - u);

			while (J > 0) {
				if (logy >= log_f_gamma(c, b, j, l)) {
					break;
				}
				l = l - gW_gamma[b][j];
				J--;
			}

			while (K > 0) {
				if (logy >= log_f_gamma(c, b, j, r)) {
					break;
				}
				r = r + gW_gamma[b][j];
				K--;
			}

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double cand = runif(l, r);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			while (logy >= log_f_gamma(c, b, j, cand)) {
				if (cand < gGamma[c][b][j]) {
					l = cand;
				}
				else {
					r = cand;
				}
#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				cand = runif(l, r);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif
			}

			gGamma[c][b][j] = cand;

			if (iter >= burnin) {
				gGamma_samples[c][b][j][iter - burnin] = gGamma[c][b][j];
			}
		}
	}
}

double c2121a::log_f_theta(int c, int b, int j, double theta)
{
	// Drop the constant term from N(mu_theta[c][b], sigma^2.theta[b])
	double q = ((-1.0/2.0)*(pow((theta - mu_theta[c][b]),2.0))/sigma2_theta[c][b]);

	double f = q + ((double)y[b][j]) * (gGamma[c][b][j] + theta) - ((double)NT[b][j]) *(log(1.0 + exp(gGamma[c][b][j] + theta)));

	return(f);
}

void c2121a::sample_theta_MH(int c, int burnin, int iter)
{
	int b = 0, j = 0;
	for (b = 0; b < gNumBodySys; b++) {
		for ( j = 0; j < gNAE[b]; j++) {

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double cand = rnorm(gTheta[c][b][j], gSigma_MH_theta[b][j]);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			double f1 = log_f_theta(c, b, j, cand);
			double f2 = log_f_theta(c, b, j, gTheta[c][b][j]);

			double ratio = exp(f1 - f2);

			ratio = cMIN(ratio, 1);

			if (u <= ratio) {
				gTheta[c][b][j] = cand;
				gTheta_acc[c][b][j] = gTheta_acc[c][b][j] + 1;
			}

			if (iter >= burnin) {
				gTheta_samples[c][b][j][iter - burnin] = gTheta[c][b][j];
			}
		}
	}
}

void c2121a::sample_theta_SLICE(int c, int burnin, int iter)
{
	int K = 0, J = 0;
	int b = 0, j = 0;

	for (b = 0; b < gNumBodySys; b++) {
		for (j = 0; j < gNAE[b]; j++) {

			int m = (int)gW_theta_control[b][j];

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			J = floor(runif(0,m));
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif
			K = (m-1) - J;

			double l = 0.0, r = 0.0;
			double g = log_f_theta(c, b, j, gTheta[c][b][j]);
			double logy = 0.0;
#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double e = rexp(1);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif
			logy = g - e;

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double u = runif(0, gW_theta[b][j]);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			l = gTheta[c][b][j] - u;
			r = gTheta[c][b][j] + (gW_theta[b][j] - u);

			while (J > 0) {
				if (logy >= log_f_theta(c, b, j, l)) {
					break;
				}
				l = l - gW_theta[b][j];
				J--;
			}

			while (K > 0) {
				if (logy >= log_f_theta(c, b, j, r)) {
					break;
				}
				r = r + gW_theta[b][j];
				K--;
			}

			
#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double cand = runif(l, r);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			// Stopping condition:
			// Paper:
			//		logy < log f(x1)
			// Neal's code:
			//		log f(x1) >= logy or logy <= log f(x1)
			// So a difference here, don't think it matters much but I've used the
			// original paper's condition:
			// We only stop if logy < log f(x1) - otherwise we continue.
			// I'm not sure if this causes an infinite loop at some stage(?)
			// if the the original sampled value is at the boundary, i.e. if
			// y ~ runif(0, f(x0)) give f(x0) as its value?

			while (logy >= log_f_theta(c, b, j, cand)) {
				if (cand < gTheta[c][b][j]) {
					l = cand;
				}
				else {
					r = cand;
				}
#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				cand = runif(l, r);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif
			}

			gTheta[c][b][j] = cand;

			if (iter >= burnin) {
				gTheta_samples[c][b][j][iter - burnin] = gTheta[c][b][j];
			}
		}
	}
}

double c2121a::cMIN(double a, double b)
{
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

void c2121a::init(SEXP pChains, SEXP pBurnin, SEXP pIter, SEXP pNumBodySys,
                            SEXP pMaxAEs, SEXP pNAE,
							SEXP pSim_Type,
							SEXP pGlobal_Sim_Param, SEXP pGlobal_Sim_Param_Cntrl,
							SEXP pSim_Param,
							SEXP pX,
                            SEXP pY, SEXP pNC, SEXP pNT,
                            SEXP ptheta, SEXP pgamma, SEXP pmu_gamma_0_0,
                            SEXP ptau2_gamma_0_0,
                            SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0,
                            SEXP palpha_gamma_0_0,
                            SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0,
                            SEXP pbeta_theta_0_0,
                            SEXP palpha_gamma, SEXP pbeta_gamma,
                            SEXP palpha_theta,
                            SEXP pbeta_theta, SEXP pmu_gamma_0,
                            SEXP ptau2_gamma_0,
                            SEXP pmu_theta_0, SEXP ptau2_theta_0,
                            SEXP pmu_gamma, SEXP pmu_theta,
                            SEXP psigma2_gamma, SEXP psigma2_theta)
{
	release();

	gChains = *(INTEGER(pChains));
	gBurnin = *(INTEGER(pBurnin));
	gIter = *(INTEGER(pIter));
	gNumBodySys = *(INTEGER(pNumBodySys));
	gMaxAEs = *(INTEGER(pMaxAEs));

	int l = 0;
	gNAE = (int *)malloc(gNumBodySys * sizeof(int));
	for (l = 0; l < gNumBodySys; l++) {
		gNAE[l] = (INTEGER(pNAE))[l];
	}

	alpha_gamma_0_0 = *(REAL(palpha_gamma_0_0));
	beta_gamma_0_0 = *(REAL(pbeta_gamma_0_0));
	alpha_theta_0_0 = *(REAL(palpha_theta_0_0));
	beta_theta_0_0 = *(REAL(pbeta_theta_0_0));

	alpha_gamma = *(REAL(palpha_gamma));
	beta_gamma = *(REAL(pbeta_gamma));
	alpha_theta = *(REAL(palpha_theta));
	beta_theta = *(REAL(pbeta_theta));

	mu_theta_0_0 = *(REAL(pmu_theta_0_0));
	mu_gamma_0_0 = *(REAL(pmu_gamma_0_0));
	tau2_theta_0_0 = *(REAL(ptau2_theta_0_0));
	tau2_gamma_0_0 = *(REAL(ptau2_gamma_0_0));

	int c = 0;
	mu_gamma_0 = (double *)malloc(gChains * sizeof(double));
	mu_theta_0 = (double *)malloc(gChains * sizeof(double));
	tau2_gamma_0 = (double *)malloc(gChains * sizeof(double));
	tau2_theta_0 = (double *)malloc(gChains * sizeof(double));

	double *vmu_gamma_0 = REAL(pmu_gamma_0);
	double *vmu_theta_0 = REAL(pmu_theta_0);
	double *vtau2_gamma_0 = REAL(ptau2_gamma_0);
	double *vtau2_theta_0 = REAL(ptau2_theta_0);
	for (c = 0; c < gChains; c++) {
		mu_gamma_0[c] = *vmu_gamma_0++;
		mu_theta_0[c] = *vmu_theta_0++;
		tau2_gamma_0[c] = *vtau2_gamma_0++;
		tau2_theta_0[c] = *vtau2_theta_0++;
	}


	// Initialise the data
	x = (int **)malloc(gNumBodySys * sizeof(int*));
	y = (int **)malloc(gNumBodySys * sizeof(int*));
	NC = (int **)malloc(gNumBodySys * sizeof(int*));
	NT = (int **)malloc(gNumBodySys * sizeof(int*));

	gTheta = (double ***)malloc(gChains * sizeof(double**));
	gGamma = (double ***)malloc(gChains * sizeof(double**));
	gTheta_acc = (int ***)malloc(gChains * sizeof(int**));
	gGamma_acc = (int ***)malloc(gChains * sizeof(int**));

	for (c = 0; c < gChains; c++) {
		gTheta[c] = (double **)malloc(gNumBodySys * sizeof(double*));
		gGamma[c] = (double **)malloc(gNumBodySys * sizeof(double*));
		gTheta_acc[c] = (int **)malloc(gNumBodySys * sizeof(int*));
		gGamma_acc[c] = (int **)malloc(gNumBodySys * sizeof(int*));
		int i = 0;
		for (i = 0; i < gNumBodySys; i++) {
			gTheta[c][i] = (double *)malloc(gNAE[i] * sizeof(double));
			gGamma[c][i] = (double *)malloc(gNAE[i] * sizeof(double));
			gTheta_acc[c][i] = (int *)malloc(gNAE[i] * sizeof(int));
			gGamma_acc[c][i] = (int *)malloc(gNAE[i] * sizeof(int));
		}
	}

	int *vNC = INTEGER(pNC);
	int *vNT = INTEGER(pNT);
	int *vX = INTEGER(pX);
	int *vY = INTEGER(pY);
	double *vtheta = REAL(ptheta);
	double *vgamma = REAL(pgamma);

	for (c = 0; c < gChains; c++) {
		int i = 0, j = 0;
		for (i = 0; i < gNumBodySys; i++) {
			for (j = 0; j < gMaxAEs; j++) {
				if (j < gNAE[i]) {
					gTheta[c][i][j] = *vtheta;
					gGamma[c][i][j] = *vgamma;
					gGamma_acc[c][i][j] = 0;
					gTheta_acc[c][i][j] = 0;
				}
				vtheta++;
				vgamma++;
			}
		}
	}

	int i = 0, j = 0;
	for (i = 0; i < gNumBodySys; i++) {
		x[i] = (int *)malloc(gNAE[i] * sizeof(int));
		y[i] = (int *)malloc(gNAE[i] * sizeof(int));
		NC[i] = (int *)malloc(gNAE[i] * sizeof(int));
		NT[i] = (int *)malloc(gNAE[i] * sizeof(int));
		for (j = 0; j < gMaxAEs; j++) {
			if (j < gNAE[i]) {
				NC[i][j] = *vNC;
				NT[i][j] = *vNT;
				x[i][j] = *vX;
				y[i][j] = *vY;
			}
			vNC++;
			vNT++;
			vX++;
			vY++;
		}
	}

	mu_theta = (double **)malloc(gChains * sizeof(double*));
	mu_gamma = (double **)malloc(gChains * sizeof(double*));
	sigma2_theta = (double **)malloc(gChains * sizeof(double*));
	sigma2_gamma = (double **)malloc(gChains * sizeof(double*));
	for (c = 0; c < gChains; c++) {
		mu_theta[c] = (double *)malloc(gNumBodySys * sizeof(double));
		mu_gamma[c] = (double *)malloc(gNumBodySys * sizeof(double));
		sigma2_theta[c] = (double *)malloc(gNumBodySys * sizeof(double));
		sigma2_gamma[c] = (double *)malloc(gNumBodySys * sizeof(double));
	}

	double *vmu_theta = REAL(pmu_theta);
	double *vmu_gamma = REAL(pmu_gamma);
	double *vsigma2_gamma = REAL(psigma2_gamma);
	double *vsigma2_theta = REAL(psigma2_theta);

	for (c = 0; c < gChains; c++) {
		for (i = 0; i < gNumBodySys; i++) {
			mu_theta[c][i] = *vmu_theta++;
			mu_gamma[c][i] = *vmu_gamma++;
			sigma2_theta[c][i] = *vsigma2_theta++;
			sigma2_gamma[c][i] = *vsigma2_gamma++;
		}
	}

	// The samples
	mu_gamma_0_samples = (double **)malloc(gChains *sizeof(double*));
	mu_theta_0_samples = (double **)malloc(gChains *sizeof(double*));
	tau2_theta_0_samples = (double **)malloc(gChains *sizeof(double*));
	tau2_gamma_0_samples = (double **)malloc(gChains *sizeof(double*));
	for (c = 0; c < gChains; c++) {
		mu_gamma_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		mu_theta_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		tau2_theta_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		tau2_gamma_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
	}

	mu_theta_samples = (double ***)malloc(gChains *sizeof(double**));
	mu_gamma_samples = (double ***)malloc(gChains *sizeof(double**));
	sigma2_theta_samples = (double ***)malloc(gChains *sizeof(double**));
	sigma2_gamma_samples = (double ***)malloc(gChains *sizeof(double**));
	for (c = 0; c < gChains; c++) {
		mu_theta_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		mu_gamma_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		sigma2_theta_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		sigma2_gamma_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));

		for (i = 0; i < gNumBodySys; i++) {
			mu_theta_samples[c][i]
					= (double *)malloc((gIter - gBurnin) *sizeof(double));
			mu_gamma_samples[c][i]
					= (double *)malloc((gIter - gBurnin) *sizeof(double));
			sigma2_theta_samples[c][i]
					= (double *)malloc((gIter - gBurnin) *sizeof(double));
			sigma2_gamma_samples[c][i]
					= (double *)malloc((gIter - gBurnin) *sizeof(double));
		}
	}

	gTheta_samples = (double ****)malloc(gChains *sizeof(double***));
	gGamma_samples = (double ****)malloc(gChains *sizeof(double***));

	for (c = 0; c < gChains; c++) {
		gTheta_samples[c] = (double ***)malloc(gNumBodySys *sizeof(double**));
		gGamma_samples[c] = (double ***)malloc(gNumBodySys *sizeof(double**));
		for (i = 0; i < gNumBodySys; i++) {
			gTheta_samples[c][i] = (double **)malloc(gNAE[i] *sizeof(double*));
			gGamma_samples[c][i] = (double **)malloc(gNAE[i] *sizeof(double*));
			for (j = 0; j < gNAE[i]; j++) {
				gTheta_samples[c][i][j]
						= (double *)malloc((gIter - gBurnin) *sizeof(double));
				gGamma_samples[c][i][j]
						= (double *)malloc((gIter - gBurnin) *sizeof(double));
			}
		}
	}

	// Simulation Parameters

	// Global Parameters
	gSim_Param = *(REAL(pGlobal_Sim_Param));
	gSim_Param_cntrl = *(REAL(pGlobal_Sim_Param_Cntrl));
	l = strlen(CHAR(STRING_ELT(pSim_Type, 0)));
	sim_type = (char *)malloc((l + 1)*sizeof(char));
	if (sim_type) {
		strcpy(sim_type, CHAR(STRING_ELT(pSim_Type, 0)));
		sim_type[l] = 0;
	}

	// Individual Variable Parameters
	initSimParams(pSim_Param);

	// Report the global parameters
	Rprintf("Global Simulation Parameters:\n");
	Rprintf("\tSimulation Type: %s\n", sim_type);
	if (0 == strcmp("SLICE", sim_type)) {
		Rprintf("\tw (width): %0.6f\n", gSim_Param);
		Rprintf("\tm (control): %0.6f\n", gSim_Param_cntrl);
	}
	else {
		Rprintf("\tsigma_MH: %0.6f\n", gSim_Param);
	}
}

void c2121a::initSimParams(SEXP sim_params)
{
	gW_gamma = (double **)malloc(gNumBodySys * sizeof(double*));
	gW_theta = (double **)malloc(gNumBodySys * sizeof(double*));
	gW_gamma_control = (int **)malloc(gNumBodySys * sizeof(int*));
	gW_theta_control = (int **)malloc(gNumBodySys * sizeof(int*));
	gSigma_MH_gamma = (double **)malloc(gNumBodySys * sizeof(double*));
	gSigma_MH_theta = (double **)malloc(gNumBodySys * sizeof(double*));

	int i = 0, j = 0;
	for (i = 0; i < gNumBodySys; i++) {
		gW_gamma[i] = (double *)malloc(gNAE[i] * sizeof(double));
		gW_theta[i] = (double *)malloc(gNAE[i] * sizeof(double));
		gW_gamma_control[i] = (int *)malloc(gNAE[i] * sizeof(int));
		gW_theta_control[i] = (int *)malloc(gNAE[i] * sizeof(int));
		gSigma_MH_gamma[i] = (double *)malloc(gNAE[i] * sizeof(double));
		gSigma_MH_theta[i] = (double *)malloc(gNAE[i] * sizeof(double));
		for (j = 0; j < gNAE[i]; j++) {
			gW_gamma[i][j] = gSim_Param;
			gW_theta[i][j] = gSim_Param;
			gW_gamma_control[i][j] = (int)gSim_Param_cntrl;
			gW_theta_control[i][j] = (int)gSim_Param_cntrl;
			gSigma_MH_gamma[i][j] = gSim_Param;
			gSigma_MH_theta[i][j] = gSim_Param;
		}
	}

	int len = Rf_length(sim_params);

	//SEXP sType = R_NilValue;
	SEXP sVariables = R_NilValue;
	SEXP sParams = R_NilValue;
	SEXP sValues = R_NilValue;
	SEXP sControl = R_NilValue;
	SEXP sB = R_NilValue;
	SEXP sj = R_NilValue;

	if (len && Rf_isNewList(sim_params)) {

		SEXP names = Rf_getAttrib(sim_params, R_NamesSymbol);

		for (i = 0; i < len; i++) {
			if (strcmp(sColValue, CHAR(STRING_ELT(names, i))) == 0) {
				sValues = VECTOR_ELT(sim_params, i);
			}
			if (strcmp(sColParam, CHAR(STRING_ELT(names, i))) == 0) {
				sParams = VECTOR_ELT(sim_params, i);
			}
			if (strcmp(sColControl, CHAR(STRING_ELT(names, i))) == 0) {
				sControl = VECTOR_ELT(sim_params, i);
			}
			if (strcmp(sColVariable, CHAR(STRING_ELT(names, i))) == 0) {
				sVariables = VECTOR_ELT(sim_params, i);
			}
			if (strcmp(sColB, CHAR(STRING_ELT(names, i))) == 0) {
				sB = VECTOR_ELT(sim_params, i);
			}
			if (strcmp(sColj, CHAR(STRING_ELT(names, i))) == 0) {
				sj = VECTOR_ELT(sim_params, i);
			}
		}

		len = Rf_length(sParams);

		if (len > 0) {

			double* vals = REAL(sValues);
			double* cntrl = REAL(sControl);
			int* B = INTEGER(sB);
			int* j = INTEGER(sj);
		
			for (i = 0; i < len; i++) {
				const char *var = CHAR(STRING_ELT(sVariables, i));
				const char *param = CHAR(STRING_ELT(sParams, i));
	
				int b = B[i] - 1;
				int a = j[i] - 1;
				if (0 == strcmp(sVariable_gamma, var)) {
					if (0 == strcmp(param, sParam_w)) {
						gW_gamma[b][a] = vals[i];
						gW_gamma_control[b][a] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH)) {
						gSigma_MH_gamma[b][a] = vals[i];
					}
				}
				else if (0 == strcmp(sVariable_theta, var)) {
					if (0 == strcmp(param, sParam_w)) {
						gW_theta[b][a] = vals[i];
						gW_theta_control[b][a] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH)) {
						gSigma_MH_theta[b][a] = vals[i];
					}
				}
			}
		}
	}
}

void c2121a::release()
{
	int i = 0, j = 0, c = 0;

	if (x != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(x[i]);
		}
		free(x);
		x = NULL;
	}

	if (y != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(y[i]);
		}
		free(y);
		y = NULL;
	}

	if (NC != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(NC[i]);
		}
		free(NC);
		NC = NULL;
	}

	if (NT != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(NT[i]);
		}
		free(NT);
		NT = NULL;
	}

	if (gTheta != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(gTheta[c][i]);
			}
			free(gTheta[c]);
		}
		free(gTheta);
		gTheta = NULL;
	}

	if (gGamma != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(gGamma[c][i]);
			}
			free(gGamma[c]);
		}
		free(gGamma);
		gGamma = NULL;
	}

	if (gTheta_acc != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(gTheta_acc[c][i]);
			}
			free(gTheta_acc[c]);
		}
		free(gTheta_acc);
		gTheta_acc = NULL;
	}

	if (gGamma_acc != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(gGamma_acc[c][i]);
			}
			free(gGamma_acc[c]);
		}
		free(gGamma_acc);
		gGamma_acc = NULL;
	}

	if (mu_theta != NULL) {
		for (c = 0; c < gChains; c++) {
			free(mu_theta[c]);
		}
		free(mu_theta);
		mu_theta = NULL;
	}

	if (mu_gamma != NULL) {
		for (c = 0; c < gChains; c++) {
			free(mu_gamma[c]);
		}
		free(mu_gamma);
		mu_gamma = NULL;
	}

	if (sigma2_theta != NULL) {
		for (c = 0; c < gChains; c++) {
			free(sigma2_theta[c]);
		}
		free(sigma2_theta);
		sigma2_theta = NULL;
	}

	if (sigma2_gamma != NULL) {
		for (c = 0; c < gChains; c++) {
			free(sigma2_gamma[c]);
		}
		free(sigma2_gamma);
		sigma2_gamma = NULL;
	}

	if (mu_gamma_0 != NULL) {
		free(mu_gamma_0);
		mu_gamma_0 = NULL;
	}

	if (mu_theta_0 != NULL) {
		free(mu_theta_0);
		mu_theta = NULL;
	}

	if (tau2_theta_0 != NULL) {
		free(tau2_theta_0);
		tau2_theta_0 = NULL;
	}

	if (tau2_gamma_0 != NULL) {
		free(tau2_gamma_0);
		tau2_gamma_0 = NULL;
	}

	// The samples
	if (gTheta_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				for (j = 0; j < gNAE[i]; j++) {
					free(gTheta_samples[c][i][j]);
				}
				free(gTheta_samples[c][i]);
			}
			free(gTheta_samples[c]);
		}
		free(gTheta_samples);
		gTheta_samples = NULL;
	}

	if (gGamma_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				for (j = 0; j < gNAE[i]; j++) {
					free(gGamma_samples[c][i][j]);
				}
				free(gGamma_samples[c][i]);
			}
			free(gGamma_samples[c]);
		}
		free(gGamma_samples);
		gGamma_samples = NULL;
	}

	if (mu_theta_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(mu_theta_samples[c][i]);
			}
			free(mu_theta_samples[c]);
		}
		free(mu_theta_samples);
		mu_theta_samples = NULL;
	}

	if (mu_gamma_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(mu_gamma_samples[c][i]);
			}
			free(mu_gamma_samples[c]);
		}
		free(mu_gamma_samples);
		mu_gamma_samples = NULL;
	}

	if (sigma2_theta_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(sigma2_theta_samples[c][i]);
			}
			free(sigma2_theta_samples[c]);
		}
		free(sigma2_theta_samples);
		sigma2_theta_samples = NULL;
	}

	if (sigma2_gamma_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(sigma2_gamma_samples[c][i]);
			}
			free(sigma2_gamma_samples[c]);
		}
		free(sigma2_gamma_samples);
		sigma2_gamma_samples = NULL;
	}

	if (mu_gamma_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			free(mu_gamma_0_samples[c]);
		}
		free(mu_gamma_0_samples);
		mu_gamma_0_samples = NULL;
	}

	if (mu_theta_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			free(mu_theta_0_samples[c]);
		}
		free(mu_theta_0_samples);
		mu_theta_0_samples = NULL;
	}

	if (tau2_theta_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			free(tau2_theta_0_samples[c]);
		}
		free(tau2_theta_0_samples);
		tau2_theta_0_samples = NULL;
	}

	if (tau2_gamma_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			free(tau2_gamma_0_samples[c]);
		}
		free(tau2_gamma_0_samples);
		tau2_gamma_0_samples = NULL;
	}

	if (gNAE != NULL) {
		free(gNAE);
		gNAE = NULL;
	}

	// Free the simulation parameters
	if (sim_type != NULL) {
		free(sim_type);
		sim_type = NULL;
	}

	if (gW_gamma != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gW_gamma[i]);
		}
		free(gW_gamma);
		gW_gamma = NULL;
	}

	if (gW_theta != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gW_theta[i]);
		}
		free(gW_theta);
		gW_theta = NULL;
	}

	if (gW_gamma_control != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gW_gamma_control[i]);
		}
		free(gW_gamma_control);
		gW_gamma_control = NULL;
	}

	if (gW_theta_control != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gW_theta_control[i]);
		}
		free(gW_theta_control);
		gW_theta_control = NULL;
	}

	if (gSigma_MH_gamma != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gSigma_MH_gamma[i]);
		}
		free(gSigma_MH_gamma);
		gSigma_MH_gamma = NULL;
	}

	if (gSigma_MH_theta != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gSigma_MH_theta[i]);
		}
		free(gSigma_MH_theta);
		gSigma_MH_theta = NULL;
	}
}

SEXP c2121a::getThetaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL1Samples(gTheta_samples);

	return samples;
}

SEXP c2121a::getGammaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL1Samples(gGamma_samples);

	return samples;
}

SEXP c2121a::getMuThetaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(mu_theta_samples);

	return samples;
}

SEXP c2121a::getMuGammaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(mu_gamma_samples);

	return samples;
}

SEXP c2121a::getSigma2ThetaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(sigma2_theta_samples);

	return samples;
}

SEXP c2121a::getSigma2GammaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(sigma2_gamma_samples);

	return samples;
}

SEXP c2121a::getL1Samples(double**** &data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = Rf_allocVector(REALSXP, gChains * gNumBodySys * gMaxAEs * (gIter - gBurnin)));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		int b = 0;
		for (b = 0; b < gNumBodySys; b++) {
			int j = 0;
			for (j = 0; j < gMaxAEs; j++) {
				if (j < gNAE[b]) {
					memcpy(REAL(samples) + i, data[c][b][j], (gIter - gBurnin)*sizeof(double));
					free(data[c][b][j]);
					data[c][b][j] = NULL;
				}
				i += (gIter - gBurnin);
			}
			free(data[c][b]);
			data[c][b] = NULL;
		}
		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = Rf_allocVector(INTSXP, 4));

	INTEGER(dim)[0] = (gIter - gBurnin);
	INTEGER(dim)[1] = gMaxAEs;
	INTEGER(dim)[2] = gNumBodySys;
	INTEGER(dim)[3] = gChains;

	Rf_setAttrib(samples, R_DimSymbol, dim);

	UNPROTECT(2);

	return samples;
}

SEXP c2121a::getL2Samples(double*** & data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = Rf_allocVector(REALSXP, gChains * gNumBodySys * (gIter - gBurnin)));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		int b = 0;
		for (b = 0; b < gNumBodySys; b++) {
			memcpy(REAL(samples) + i, data[c][b], (gIter - gBurnin)*sizeof(double));
			i += (gIter - gBurnin);
			free(data[c][b]);
			data[c][b] = NULL;
		}
		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = Rf_allocVector(INTSXP, 3));

	INTEGER(dim)[0] = (gIter - gBurnin);
	INTEGER(dim)[1] = gNumBodySys;
	INTEGER(dim)[2] = gChains;

	Rf_setAttrib(samples, R_DimSymbol, dim);

	UNPROTECT(2);

	return samples;
}

SEXP c2121a::getL3Samples(double** &data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = Rf_allocVector(REALSXP, gChains * (gIter - gBurnin)));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		memcpy(REAL(samples) + i, data[c], (gIter - gBurnin)*sizeof(double));
		i += (gIter - gBurnin);
		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = Rf_allocVector(INTSXP, 2));

	INTEGER(dim)[0] = (gIter - gBurnin);
	INTEGER(dim)[1] = gChains;

	Rf_setAttrib(samples, R_DimSymbol, dim);

	UNPROTECT(2);

	return samples;
}

SEXP c2121a::getMuTheta0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(mu_theta_0_samples);

	return samples;
}

SEXP c2121a::getTau2Gamma0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(tau2_gamma_0_samples);

	return samples;
}

SEXP c2121a::getMuGamma0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(mu_gamma_0_samples);

	return samples;
}

SEXP c2121a::getTau2Theta0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(tau2_theta_0_samples);

	return samples;
}

SEXP c2121a::getL1Accept(int*** &data)
{
	SEXP acc = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(acc = Rf_allocVector(INTSXP, gChains * gNumBodySys * gMaxAEs));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		int b = 0;
		for (b = 0; b < gNumBodySys; b++) {
			memcpy(INTEGER(acc) + i, data[c][b], gNAE[b]*sizeof(int));
			i += gMaxAEs;
			free(data[c][b]);
			data[c][b] = NULL;
		}
		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = Rf_allocVector(INTSXP, 3));

	INTEGER(dim)[0] = gMaxAEs;
	INTEGER(dim)[1] = gNumBodySys;
	INTEGER(dim)[2] = gChains;

	Rf_setAttrib(acc, R_DimSymbol, dim);

	UNPROTECT(2);

	return acc;
}

SEXP c2121a::getThetaAccept()
{
	SEXP acc = R_NilValue;

	acc = getL1Accept(gTheta_acc);

	return acc;
}

SEXP c2121a::getGammaAccept()
{
	SEXP acc = R_NilValue;

	acc = getL1Accept(gGamma_acc);

	return acc;
}

void c2121a::getThetaSamples(int* c, int* b, int* j, double* theta_samples)
{
	int C = (*c) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

	memcpy(theta_samples, gTheta_samples[C][B][J], (gIter - gBurnin)*sizeof(double));
}

double* c2121a::getThetaSamples(int c, int b, int j)
{
	return(gTheta_samples[c][b][j]);
}

void c2121a::getGammaSamples(int* c, int* b, int* j, double* gamma_samples)
{
	int C = (*c) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

	memcpy(gamma_samples, gGamma_samples[C][B][J], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getMuThetaSamples(int* c, int* b, double* mu_theta)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	memcpy(mu_theta, mu_theta_samples[C][B], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getMuGammaSamples(int* c, int* b, double* mu_gamma)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	memcpy(mu_gamma, mu_gamma_samples[C][B], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getSigma2ThetaSamples(int* c, int* b, double* sigma2)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	memcpy(sigma2, sigma2_theta_samples[C][B], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getSigma2GammaSamples(int* c, int* b, double* sigma2)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	memcpy(sigma2, sigma2_gamma_samples[C][B], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getMuGamma0Samples(int *c, double* mu)
{
	int C = (*c) - 1;

	memcpy(mu, mu_gamma_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getMuTheta0Samples(int* c, double* mu)
{
	int C = (*c) - 1;

	memcpy(mu, mu_theta_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getTau2Gamma0Samples(int* c, double* tau2)
{
	int C = (*c) - 1;

	memcpy(tau2, tau2_gamma_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getTau2Theta0Samples(int* c, double* tau2)
{
	int C = (*c) - 1;

	memcpy(tau2, tau2_theta_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c2121a::getThetaAccept(int* c, int* b, int* j, double* theta_acc)
{
	int C = (*c) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

	*theta_acc = gTheta_acc[C][B][J];
}

void c2121a::getGammaAccept(int* c, int* b, int* j, double* gamma_acc)
{
	int C = (*c) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

	*gamma_acc = gGamma_acc[C][B][J];
}
