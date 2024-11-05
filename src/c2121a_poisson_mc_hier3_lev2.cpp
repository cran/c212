#include<cstdio>
#include<cstdlib>
#include <cstring>
#include<cmath>

#include "c212_Rdefines.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <Rinternals.h>


#include "c2121a_poisson_mc_hier2_lev0.h"
#include "c2121a_poisson_mc_hier3_lev0.h"
#include "c2121a_poisson_mc_hier3_lev2.h"

using namespace std;

//static const char *rcsId = "$Id: c2121a_poisson_mc_hier3_lev2.cpp,v 1.17 2018/10/03 15:40:27 clb13102 Exp clb13102 $";

c2121a_poisson_mc_hier3_lev2::c2121a_poisson_mc_hier3_lev2()
{
	//Rprintf("c2121a_poisson_mc_hier3_lev2::c2121a_poisson_mc_hier3_lev2: Default constructor\n");

	mu_theta_0 = NULL;
	mu_gamma_0 = NULL;
	tau2_theta_0 = NULL;
	tau2_gamma_0 = NULL;

	mu_theta_0_samples = NULL;
	mu_gamma_0_samples = NULL;
	tau2_theta_0_samples = NULL;
	tau2_gamma_0_samples = NULL;
}

c2121a_poisson_mc_hier3_lev2::c2121a_poisson_mc_hier3_lev2(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_cntrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE, SEXP pX,
					SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma, SEXP pmu_gamma_0_0,
					SEXP ptau2_gamma_0_0, SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0, SEXP palpha_gamma_0_0,
					SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0, SEXP pbeta_theta_0_0, SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta, SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0, SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta)

{
	mu_theta_0 = NULL;
	mu_gamma_0 = NULL;
	tau2_theta_0 = NULL;
	tau2_gamma_0 = NULL;

	mu_theta_0_samples = NULL;
	mu_gamma_0_samples = NULL;
	tau2_theta_0_samples = NULL;
	tau2_gamma_0_samples = NULL;

	init(sChains, sBurnin, sIter, sSim_Type, sMem_Model, sGlobal_Sim_Param,
				sGlobal_Sim_Param_cntrl,
				sSim_Param,
				sMonitor,
				sNumIntervals, sMaxBs, sNumBodySys, sMaxAEs, sNAE,
				pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0, pmu_theta_0_0, ptau2_theta_0_0,
				palpha_gamma_0_0, pbeta_gamma_0_0, palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma,
				pbeta_gamma, palpha_theta, pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0,
				ptau2_theta_0, pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta);

}

void c2121a_poisson_mc_hier3_lev2::clear()
{
	release();
	c2121a_poisson_mc_hier3_lev0::release();
	c2121a_poisson_mc_hier2_lev0::release();
}

void c2121a_poisson_mc_hier3_lev2::initL3Samples()
{
	int c = 0;

	// The samples
	if (retainSamples(iMonitor_mu_gamma_0))
		mu_gamma_0_samples = (double **)malloc(gChains * sizeof(double*));
	if (retainSamples(iMonitor_mu_theta_0))
		mu_theta_0_samples = (double **)malloc(gChains * sizeof(double*));
	if (retainSamples(iMonitor_tau2_gamma_0))
		tau2_gamma_0_samples = (double **)malloc(gChains * sizeof(double*));
	if (retainSamples(iMonitor_tau2_theta_0))
		tau2_theta_0_samples = (double **)malloc(gChains * sizeof(double*));

	for (c = 0; c < gChains; c++) {
		if (retainSamples(iMonitor_mu_gamma_0))
			mu_gamma_0_samples[c] = (double *)malloc((gIter - gBurnin)* sizeof(double));
		if (retainSamples(iMonitor_mu_theta_0))
			mu_theta_0_samples[c] = (double *)malloc((gIter - gBurnin)* sizeof(double));
		if (retainSamples(iMonitor_tau2_gamma_0))
			tau2_gamma_0_samples[c] =
								(double *)malloc((gIter - gBurnin)* sizeof(double));
		if (retainSamples(iMonitor_tau2_theta_0))
			tau2_theta_0_samples[c] =
								(double *)malloc((gIter - gBurnin)* sizeof(double));
	}
}

void c2121a_poisson_mc_hier3_lev2::initL3Variables(SEXP pmu_gamma_0,
            SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0)
{
	int c = 0;

	mu_gamma_0 = (double*)malloc(gChains * sizeof(double));
	double *vmu_gamma_0 = REAL(pmu_gamma_0);
	for (c = 0; c < gChains; c++) {
		mu_gamma_0[c] = *vmu_gamma_0;
		vmu_gamma_0++;
	}

	mu_theta_0 = (double*)malloc(gChains * sizeof(double));
	double *vmu_theta_0 = REAL(pmu_theta_0);
	for (c = 0; c < gChains; c++) {
		mu_theta_0[c] = *vmu_theta_0;
		vmu_theta_0++;
	}

	tau2_gamma_0 = (double*)malloc(gChains * sizeof(double));
	double *vtau2_gamma_0 = REAL(ptau2_gamma_0);
	for (c = 0; c < gChains; c++) {
		tau2_gamma_0[c] = *vtau2_gamma_0;
		vtau2_gamma_0++;
	}

	tau2_theta_0 = (double*)malloc(gChains * sizeof(double));
	double *vtau2_theta_0 = REAL(ptau2_theta_0);
	for (c = 0; c < gChains; c++) {
		tau2_theta_0[c] = *vtau2_theta_0;
		vtau2_theta_0++;
	}
}

c2121a_poisson_mc_hier3_lev2::~c2121a_poisson_mc_hier3_lev2()
{
	//Rprintf("c2121a_poisson_mc_hier3_lev2::c2121a_poisson_mc_hier3_lev2 - destructor\n");
	release();
}

void c2121a_poisson_mc_hier3_lev2::gibbs_sampler()
{
	if (strcmp(sim_type, "MH") == 0) {
		simulate_MH();
	}
	else {
		simulate_SLICE();
	}

	return;
}

void c2121a_poisson_mc_hier3_lev2::simulate_MH()
{
	int i = 0;

	for (i = 0; i < gIter; i++) {
#ifndef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		sample_mu_gamma_0(gBurnin, i);
		sample_mu_theta_0(gBurnin, i);
		sample_tau2_gamma_0(gBurnin, i);
		sample_tau2_theta_0(gBurnin, i);
		sample_mu_gamma(gBurnin, i);
		sample_mu_theta(gBurnin, i);
		sample_sigma2_gamma(gBurnin, i);
		sample_sigma2_theta(gBurnin, i);
		sample_gamma_MH(gBurnin, i);
		sample_theta_MH(gBurnin, i);
#ifndef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		if ((i + 1)%1000 == 0) {
			Rprintf("%d iterations...\n", i + 1);
		}
	}
	Rprintf("MCMC fitting complete.\n");
}

void c2121a_poisson_mc_hier3_lev2::simulate_SLICE()
{
	int i = 0;

	for (i = 0; i < gIter; i++) {
#ifndef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		sample_mu_gamma_0(gBurnin, i);
		sample_mu_theta_0(gBurnin, i);
		sample_tau2_gamma_0(gBurnin, i);
		sample_tau2_theta_0(gBurnin, i);
		sample_mu_gamma(gBurnin, i);
		sample_mu_theta(gBurnin, i);
		sample_sigma2_gamma(gBurnin, i);
		sample_sigma2_theta(gBurnin, i);
		sample_gamma_SLICE(gBurnin, i);
		sample_theta_SLICE(gBurnin, i);
#ifndef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		if ((i + 1)%1000 == 0) {
			Rprintf("%d iterations...\n", i + 1);
		}
	}
	Rprintf("MCMC fitting complete.\n");
}

void c2121a_poisson_mc_hier3_lev2::sample_mu_gamma_0(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {

		int sum_B = 0;
		double mu_gamma_tot = 0.0;

		for (l = 0; l < gNumIntervals; l++) {
			sum_B = sum_B + gNumBodySys[l];

			int b = 0;
			for (b = 0; b < gNumBodySys[l]; b++) {
				mu_gamma_tot += mu_gamma[c][l][b];
			}
		}

		double denom = tau2_gamma_0[c] + tau2_gamma_0_0 * (double)sum_B;

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

		if (iter >= burnin  && retainSamples(iMonitor_mu_gamma_0)) {
			mu_gamma_0_samples[c][iter - burnin] = mu_gamma_0[c];
		}
	}
}

void c2121a_poisson_mc_hier3_lev2::sample_mu_theta_0 (int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {

		int sum_B = 0;
		double mu_theta_tot = 0.0;

		for (l = 0; l < gNumIntervals; l++) {
			sum_B = sum_B + gNumBodySys[l];

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {
				mu_theta_tot += mu_theta[c][l][b];
			}
		}

		double denom = tau2_theta_0[c] + tau2_theta_0_0 * ((double)sum_B);


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

		if (iter >= burnin  && retainSamples(iMonitor_mu_theta_0)) {
			mu_theta_0_samples[c][iter - burnin] = mu_theta_0[c];
		}
	}
}

void c2121a_poisson_mc_hier3_lev2::sample_tau2_gamma_0(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {
		int sum_B = 0;
		double isum = 0.0;
		for (l = 0; l < gNumIntervals; l++) {
			sum_B = sum_B + gNumBodySys[l];

			int b = 0;
			for (b = 0; b < gNumBodySys[l]; b++) {
				isum += (pow((mu_gamma[c][l][b] - mu_gamma_0[c]), 2.0));
			}
		}

		double s = alpha_gamma_0_0 + ((double)sum_B)/2.0;
		double r = 0.0;

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

		if (iter >= burnin && retainSamples(iMonitor_tau2_gamma_0)) {
			tau2_gamma_0_samples[c][iter - burnin] = tau2_gamma_0[c];
		}
	}
}

void c2121a_poisson_mc_hier3_lev2::sample_tau2_theta_0(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {
		int sum_B = 0;
		double isum = 0.0;
		for (l = 0; l < gNumIntervals; l++) {
			sum_B = sum_B + gNumBodySys[l];

			int b = 0;
			for (b = 0; b < gNumBodySys[l]; b++) {
				isum += (pow((mu_theta[c][l][b] - mu_theta_0[c]), 2.0));
			}
		}

		double s = alpha_theta_0_0 + ((double)sum_B)/2.0;
		double r = 0.0;


		r = beta_theta_0_0 + 0.5 * isum;

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		tau2_theta_0[c] = 1/cand;


		if (iter >= burnin && retainSamples(iMonitor_tau2_theta_0)) {
			tau2_theta_0_samples[c][iter - burnin] = tau2_theta_0[c];
		}
	}
}

void c2121a_poisson_mc_hier3_lev2::sample_mu_gamma(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {
				double denom = sigma2_gamma[c][l][b] + ((double)gNAE[l][b])*tau2_gamma_0[c];


				double t = 0.0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					t += gGamma[c][l][b][j];
				}


				double mean = (sigma2_gamma[c][l][b] * mu_gamma_0[c] + tau2_gamma_0[c] * t)/denom;

				double var = (sigma2_gamma[c][l][b]*tau2_gamma_0[c])/denom;

				double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				double cand = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif

				mu_gamma[c][l][b] = cand;

				if (iter >= burnin && retainSamples(iMonitor_mu_gamma)) {
					mu_gamma_samples[c][l][b][iter - burnin] = mu_gamma[c][l][b];
				}
			}
		}
	}
}

void c2121a_poisson_mc_hier3_lev2::sample_mu_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {
				double denom = sigma2_theta[c][l][b] + ((double)gNAE[l][b])*tau2_theta_0[c];

				double t = 0.0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					t += gTheta[c][l][b][j];
				}

				double mean = (sigma2_theta[c][l][b] * mu_theta_0[c] + tau2_theta_0[c] * t)/denom;

				double var = (sigma2_theta[c][l][b]*tau2_theta_0[c])/denom;

				double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				double cand = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif

				mu_theta[c][l][b] = cand;

				if (iter >= burnin && retainSamples(iMonitor_mu_theta)) {
					mu_theta_samples[c][l][b][iter - burnin] = mu_theta[c][l][b];
				}

			}
		}
	}
}

void c2121a_poisson_mc_hier3_lev2::sample_sigma2_gamma(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;
			for (b = 0; b < gNumBodySys[l]; b++) {

				double s = alpha_gamma + ((double)gNAE[l][b])/2.0;

				double t = 0.0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					t += (pow(gGamma[c][l][b][j] - mu_gamma[c][l][b],2.0));
				}

				double r = beta_gamma + t/2.0;


#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif

				sigma2_gamma[c][l][b] = 1/cand;

				if (iter >= burnin && retainSamples(iMonitor_sigma2_gamma)) {
					sigma2_gamma_samples[c][l][b][iter - burnin] = sigma2_gamma[c][l][b];
				}
			}
		}
	}
}

void c2121a_poisson_mc_hier3_lev2::sample_sigma2_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {

				double s = alpha_theta + ((double)gNAE[l][b])/2.0;


				double t = 0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					t += (pow((gTheta[c][l][b][j] - mu_theta[c][l][b]),2.0));
				}

				double r = beta_theta + t/2.0;

#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif

				sigma2_theta[c][l][b] = 1/cand;


				if (iter >= burnin && retainSamples(iMonitor_sigma2_theta)) {
					sigma2_theta_samples[c][l][b][iter - burnin] = sigma2_theta[c][l][b];
				}
			}
		}
	}
}

double c2121a_poisson_mc_hier3_lev2::cMIN(double a, double b)
{
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

void c2121a_poisson_mc_hier3_lev2::releaseL3Variables()
{
	if (mu_theta_0 != NULL) {
		free(mu_theta_0);
		mu_theta_0 = NULL;
	}

	if (mu_gamma_0 != NULL) {
		free(mu_gamma_0);
		mu_gamma_0 = 0;
	}

	if (tau2_theta_0 != NULL) {
		free(tau2_theta_0);
		tau2_theta_0 = NULL;
	}

	if (tau2_gamma_0 != NULL) {
		free(tau2_gamma_0);
		tau2_gamma_0 = NULL;
	}
}

void c2121a_poisson_mc_hier3_lev2::releaseL3Samples()
{
	int c = 0;

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
	if (tau2_gamma_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			free(tau2_gamma_0_samples[c]);
		}
		free(tau2_gamma_0_samples);
		tau2_gamma_0_samples = NULL;
	}
	if (tau2_theta_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			free(tau2_theta_0_samples[c]);
		}
		free(tau2_theta_0_samples);
		tau2_theta_0_samples = NULL;
	}
}

void c2121a_poisson_mc_hier3_lev2::release()
{
	releaseL3Variables();

	releaseL3Samples();
}

SEXP c2121a_poisson_mc_hier3_lev2::getL3Samples(double** &data)
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

SEXP c2121a_poisson_mc_hier3_lev2::getMuGamma0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(mu_gamma_0_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier3_lev2::getMuTheta0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(mu_theta_0_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier3_lev2::getTau2Gamma0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(tau2_gamma_0_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier3_lev2::getTau2Theta0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(tau2_theta_0_samples);

	return samples;
}


void c2121a_poisson_mc_hier3_lev2::getMuGamma0Samples(int *c, int *l, double* mu)
{
	int C = (*c) - 1;

	if (mu_gamma_0_samples)
		memcpy(mu, mu_gamma_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier3_lev2::getMuTheta0Samples(int *c, int *l, double* mu)
{
	int C = (*c) - 1;

	if (mu_theta_0_samples)
		memcpy(mu, mu_theta_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier3_lev2::getTau2Gamma0Samples(int *c, int *l, double* tau2)
{
	int C = (*c) - 1;

	if (tau2_gamma_0_samples)
		memcpy(tau2, tau2_gamma_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier3_lev2::getTau2Theta0Samples(int *c, int *l, double* tau2)
{
	int C = (*c) - 1;

	if (tau2_theta_0_samples)
		memcpy(tau2, tau2_theta_0_samples[C], (gIter - gBurnin)*sizeof(double));
}
