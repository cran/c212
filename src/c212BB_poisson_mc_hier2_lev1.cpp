#include<cstdio>
#include<cstdlib>

#include<cstring>
#include<cmath>

#include "c212_Rdefines.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <Rinternals.h>


#include "c2121a_poisson_mc_hier2_lev0.h"
#include "c212BB_poisson_mc_hier2_lev0.h"
#include "c212BB_poisson_mc_hier2_lev1.h"

using namespace std;

//static const char *rcsId = "$Id: c212BB_poisson_mc_hier2_lev1.cpp,v 1.8 2018/10/03 15:40:28 clb13102 Exp clb13102 $";

c212BB_poisson_mc_hier2_lev1::c212BB_poisson_mc_hier2_lev1()
{
	//Rprintf("c212BB_poisson_mc_hier2_lev1::c212BB_poisson_mc_hier2_lev1: Default constructor\n");
	mu_theta = NULL;
	mu_gamma = NULL;
	sigma2_theta = NULL;
	sigma2_gamma = NULL;
	gPi = NULL;

	mu_theta_samples = NULL;
	mu_gamma_samples = NULL;
	sigma2_theta_samples = NULL;
	sigma2_gamma_samples = NULL;
	gPi_samples = NULL;
}

c212BB_poisson_mc_hier2_lev1::c212BB_poisson_mc_hier2_lev1(SEXP sChains, SEXP sBurnin,
		SEXP sIter, SEXP sSim_Type,
		SEXP sMem_Model,
		SEXP sGlobal_Sim_Params,
		SEXP sSim_Params,
		SEXP MH_weight,
		SEXP pm_weights,
		SEXP sMonitor,
		SEXP sNumIntervals, SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
		SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
		SEXP palpha_gamma,
		SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta, SEXP pmu_gamma_0,
		SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0, SEXP pmu_gamma,
		SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
		SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi,
		SEXP palgo, SEXP padapt_phase)

{
	mu_theta = NULL;
	mu_gamma = NULL;
	sigma2_theta = NULL;
	sigma2_gamma = NULL;
	gPi = NULL;

	mu_theta_samples = NULL;
	mu_gamma_samples = NULL;
	gPi_samples = NULL;
	sigma2_theta_samples = NULL;
	sigma2_gamma_samples = NULL;

	init(sChains, sBurnin, sIter, sSim_Type, sMem_Model, sGlobal_Sim_Params,
				sSim_Params,
				MH_weight,
				pm_weights,
				sMonitor,
				sNumIntervals, sMaxBs, sNumBodySys, sMaxAEs, sNAE,
				pX, pY, pC, pT, ptheta, pgamma,
				palpha_gamma,
				pbeta_gamma, palpha_theta, pbeta_theta, pmu_gamma_0, ptau2_gamma_0,
				pmu_theta_0,
				ptau2_theta_0, pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
				pPi, palpha_pi, pbeta_pi,
				palgo, padapt_phase);
}

void c212BB_poisson_mc_hier2_lev1::clear()
{
	release();
	c212BB_poisson_mc_hier2_lev0::release();
	c2121a_poisson_mc_hier2_lev0::release();
}

void c212BB_poisson_mc_hier2_lev1::initL2Variables(SEXP pmu_gamma, SEXP pmu_theta,
									SEXP psigma2_gamma, SEXP psigma2_theta, SEXP pPi)
{
	int c = 0, b = 0;

	double* vmu_gamma = REAL(pmu_gamma);
	mu_gamma = (double**)malloc(gChains * sizeof(double*));
	for (c = 0; c < gChains; c++) {
		mu_gamma[c] = (double*)malloc(gMaxBs * sizeof(double));
		for (b = 0; b < gMaxBs; b++) {
			mu_gamma[c][b] = *vmu_gamma;
			vmu_gamma++;
		}
	}

	double* vmu_theta = REAL(pmu_theta);
	mu_theta = (double**)malloc(gChains * sizeof(double*));
	for (c = 0; c < gChains; c++) {
		mu_theta[c] = (double*)malloc(gMaxBs * sizeof(double));
		for (b = 0; b < gMaxBs; b++) {
			mu_theta[c][b] = *vmu_theta;
			vmu_theta++;
		}
	}

	double* vsigma2_gamma = REAL(psigma2_gamma);
	sigma2_gamma = (double**)malloc(gChains * sizeof(double*));
	for (c = 0; c < gChains; c++) {
		sigma2_gamma[c] = (double*)malloc(gMaxBs * sizeof(double));
		for (b = 0; b < gMaxBs; b++) {
			sigma2_gamma[c][b] = *vsigma2_gamma;
			vsigma2_gamma++;
		}
	}

	double* vsigma2_theta = REAL(psigma2_theta);
	sigma2_theta = (double**)malloc(gChains * sizeof(double*));
	for (c = 0; c < gChains; c++) {
		sigma2_theta[c] = (double*)malloc(gMaxBs * sizeof(double));
		for (b = 0; b < gMaxBs; b++) {
			sigma2_theta[c][b] = *vsigma2_theta;
			vsigma2_theta++;
		}
	}

	double* vpi = REAL(pPi);
	gPi = (double**)malloc(gChains * sizeof(double*));
	for (c = 0; c < gChains; c++) {
		gPi[c] = (double*)malloc(gMaxBs * sizeof(double));
		for (b = 0; b < gMaxBs; b++) {
			gPi[c][b] = *vpi;
			vpi++;
		}
	}
}

void c212BB_poisson_mc_hier2_lev1::releaseL2Variables()
{
	int c = 0;

	if (gPi != NULL) {
		for (c = 0; c < gChains; c++) {
			free(gPi[c]);
		}
		free(gPi);
		gPi = 0;
	}

	if (mu_gamma != NULL) {
		for (c = 0; c < gChains; c++) {
			free(mu_gamma[c]);
		}
		free(mu_gamma);
		mu_gamma = 0;
	}

	if (mu_theta != NULL) {
		for (c = 0; c < gChains; c++) {
			free(mu_theta[c]);
		}
		free(mu_theta);
		mu_theta = 0;
	}

	if (sigma2_gamma != NULL) {
		for (c = 0; c < gChains; c++) {
			free(sigma2_gamma[c]);
		}
		free(sigma2_gamma);
		sigma2_gamma = 0;
	}

	if (sigma2_theta != NULL) {
		for (c = 0; c < gChains; c++) {
			free(sigma2_theta[c]);
		}
		free(sigma2_theta);
		sigma2_theta = 0;
	}
}

void c212BB_poisson_mc_hier2_lev1::initL2Samples()
{
	// The samples
	int c = 0, b = 0, l = 0;

	if (retainSamples(iMonitor_mu_theta))
		mu_theta_samples = (double ***)malloc(gChains *sizeof(double**));
	if (retainSamples(iMonitor_mu_gamma))
		mu_gamma_samples = (double ***)malloc(gChains *sizeof(double**));
	if (retainSamples(iMonitor_sigma2_theta))
		sigma2_theta_samples = (double ***)malloc(gChains *sizeof(double**));
	if (retainSamples(iMonitor_sigma2_gamma))
		sigma2_gamma_samples = (double ***)malloc(gChains *sizeof(double**));
	if (retainSamples(iMonitor_pi))
		gPi_samples = (double ***)malloc(gChains *sizeof(double**));

	for (c = 0; c < gChains; c++) {
		if (retainSamples(iMonitor_mu_theta))
			mu_theta_samples[c] = (double **)malloc(gMaxBs *sizeof(double*));
		if (retainSamples(iMonitor_mu_gamma))
			mu_gamma_samples[c] = (double **)malloc(gMaxBs *sizeof(double*));
		if (retainSamples(iMonitor_sigma2_theta))
			sigma2_theta_samples[c] = (double **)malloc(gMaxBs *sizeof(double*));
		if (retainSamples(iMonitor_sigma2_gamma))
			sigma2_gamma_samples[c] = (double **)malloc(gMaxBs *sizeof(double*));
		if (retainSamples(iMonitor_pi))
			gPi_samples[c] = (double **)malloc(gMaxBs *sizeof(double*));

		for (b = 0; b < gNumBodySys[l]; b++) {
			if (retainSamples(iMonitor_mu_theta))
				mu_theta_samples[c][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));
			if (retainSamples(iMonitor_mu_gamma))
				mu_gamma_samples[c][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));

			if (retainSamples(iMonitor_sigma2_theta))
				sigma2_theta_samples[c][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));

			if (retainSamples(iMonitor_sigma2_gamma))
				sigma2_gamma_samples[c][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));

			if (retainSamples(iMonitor_pi))
				gPi_samples[c][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));
		}
	}
}

void c212BB_poisson_mc_hier2_lev1::releaseL2Samples()
{
	int c = 0, b = 0, l = 0;

	l = 0;
	if (gPi_samples) {
		for (c = 0; c < gChains; c++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(gPi_samples[c][b]);
			}
			free(gPi_samples[c]);
		}
		free(gPi_samples);
		gPi_samples = NULL;
	}

	if (mu_theta_samples) {
		for (c = 0; c < gChains; c++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(mu_theta_samples[c][b]);
			}
			free(mu_theta_samples[c]);
		}
		free(mu_theta_samples);
		mu_theta_samples = NULL;
	}
	if (mu_gamma_samples != NULL) {
		for (c = 0; c < gChains; c++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					free(mu_gamma_samples[c][b]);
				}
			free(mu_gamma_samples[c]);
		}
		free(mu_gamma_samples);
		mu_gamma_samples = NULL;
	}
	if (sigma2_theta_samples != NULL) {
		for (c = 0; c < gChains; c++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					free(sigma2_theta_samples[c][b]);
				}
			free(sigma2_theta_samples[c]);
		}
		free(sigma2_theta_samples);
		sigma2_theta_samples = NULL;
	}
	if (sigma2_gamma_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(sigma2_gamma_samples[c][b]);
			}
			free(sigma2_gamma_samples[c]);
		}
		free(sigma2_gamma_samples);
		sigma2_gamma_samples = NULL;
	}
}

c212BB_poisson_mc_hier2_lev1::~c212BB_poisson_mc_hier2_lev1()
{
	//Rprintf("c212BB_poisson_mc_hier2_lev1::c212BB_poisson_mc_hier2_lev1 - destructor\n");
	release();
}

void c212BB_poisson_mc_hier2_lev1::gibbs_sampler()
{
	switch(gSimType) {
		case eSim_Type_MH:
			simulate_MH();
		break;

		case eSim_Type_SLICE:
			simulate_SLICE();
		break;

		default:
			simulate_MH();
		break;
	}
		
	return;
}

void c212BB_poisson_mc_hier2_lev1::simulate_MH()
{
	int i = 0;

	for (i = 0; i < gIter; i++) {
#ifndef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		sample_pi(gBurnin, i);
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

void c212BB_poisson_mc_hier2_lev1::simulate_SLICE()
{
	int i = 0;

	for (i = 0; i < gIter; i++) {
#ifndef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		sample_pi(gBurnin, i);
		sample_mu_gamma(gBurnin, i);
		sample_mu_theta(gBurnin, i);
		sample_sigma2_gamma(gBurnin, i);
		sample_sigma2_theta(gBurnin, i);

		sample_gamma_SLICE(gBurnin, i);
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

void c212BB_poisson_mc_hier2_lev1::sample_pi(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		int b = 0;
		for (b = 0; b < gNumBodySys[0]; b++) {

			int theta_zero_count = 0;
			int nae = 0;
			for (l = 0; l < gNumIntervals; l++) {
				int j = 0;
				for (j = 0; j< gNAE[l][b]; j++) {
					if (gTheta[c][l][b][j] == 0.0) {
						theta_zero_count++;
					}
					nae++;
				}
			}

			double shape1 = alpha_pi + (double)theta_zero_count;
			double shape2 = beta_pi + (double)nae - (double)theta_zero_count;

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			gPi[c][b] = rbeta(shape1, shape2);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			if (iter >= burnin && retainSamples(iMonitor_pi)) {
				gPi_samples[c][b][iter - burnin] = gPi[c][b];
			}
		}
	}
}

void c212BB_poisson_mc_hier2_lev1::sample_mu_gamma(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {

		int b = 0;

		for (b = 0; b < gNumBodySys[0]; b++) {

			int nae = 0;
			for (l = 0; l < gNumIntervals; l++) {
				nae = nae + gNAE[l][b];
			}


			double t = 0.0;
			int j = 0;
			for (l = 0; l < gNumIntervals; l++) {
				for (j = 0; j < gNAE[l][b]; j++) {
					t += gGamma[c][l][b][j];
				}
			}

			double denom = sigma2_gamma[c][b] + ((double)nae)*tau2_gamma_0;

			double mean = (sigma2_gamma[c][b] * mu_gamma_0 + tau2_gamma_0 * t)/denom;

			double var = (sigma2_gamma[c][b]*tau2_gamma_0)/denom;

			double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double cand = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			mu_gamma[c][b] = cand;

			if (iter >= burnin && retainSamples(iMonitor_mu_gamma)) {
				mu_gamma_samples[c][b][iter - burnin] = mu_gamma[c][b];
			}
		}
	}
}

void c212BB_poisson_mc_hier2_lev1::sample_mu_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {

		int b = 0;

		for (b = 0; b < gNumBodySys[0]; b++) {

			double t = 0.0;
			int j = 0;
			int Kb = 0;
			for (l = 0; l < gNumIntervals; l++) {
				for (j = 0; j < gNAE[l][b]; j++) {
					if (gTheta[c][l][b][j] != 0.0) {
						Kb++;
					}
					t += gTheta[c][l][b][j];
				}
			}

			double denom = sigma2_theta[c][b] + ((double)Kb)*tau2_theta_0;

			double mean = (sigma2_theta[c][b] * mu_theta_0 + tau2_theta_0 * t)/denom;

			double var = (sigma2_theta[c][b]*tau2_theta_0)/denom;

			double sd = sqrt(var);


#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double cand = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			mu_theta[c][b] = cand;

			if (iter >= burnin && retainSamples(iMonitor_mu_theta)) {
				mu_theta_samples[c][b][iter - burnin] = mu_theta[c][b];
			}

		}
	}
}

void c212BB_poisson_mc_hier2_lev1::sample_sigma2_gamma(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {

		int b = 0;
		for (b = 0; b < gNumBodySys[0]; b++) {

			int nae = 0;
			for (l = 0; l < gNumIntervals; l++) {
				nae = nae + gNAE[l][b];
			}

			double s = alpha_gamma + ((double)nae)/2.0;

			double t = 0.0;
			int j = 0;
			for (l = 0; l < gNumIntervals; l++) {
				for (j = 0; j < gNAE[l][b]; j++) {
					t += (pow(gGamma[c][l][b][j] - mu_gamma[c][b],2.0));
				}
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

			if (iter >= burnin && retainSamples(iMonitor_sigma2_gamma)) {
				sigma2_gamma_samples[c][b][iter - burnin] = sigma2_gamma[c][b];
			}
		}
	}
}

void c212BB_poisson_mc_hier2_lev1::sample_sigma2_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {

		int b = 0;

		for (b = 0; b < gNumBodySys[0]; b++) {

			double t = 0;
			int j = 0;
			int Kb = 0;
			for (l = 0; l < gNumIntervals; l++) {
				for (j = 0; j < gNAE[l][b]; j++) {
					if (gTheta[c][l][b][j] != 0.0) {
						Kb++;
						t += (pow((gTheta[c][l][b][j] - mu_theta[c][b]),2.0));
					}
				}
			}

			double s = alpha_theta + ((double)Kb)/2.0;
			double r = beta_theta + t/2.0;

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			sigma2_theta[c][b] = 1/cand;

			if (iter >= burnin && retainSamples(iMonitor_sigma2_theta)) {
				sigma2_theta_samples[c][b][iter - burnin] = sigma2_theta[c][b];
			}
		}
	}
}

double c212BB_poisson_mc_hier2_lev1::log_f_gamma(int c, int i, int b, int j, double gamm)
{
	double f1 = 0.0, f2 = 0.0, f3 = 0.0, f4 = 0.0, f5 = 0.0;

	f1 = ((double)x[i][b][j]) * gamm;
	f2 = -(exp(gamm)) * ((double)C[i][b][j]);
	f3 = ((double)y[i][b][j]) * (gamm + gTheta[c][i][b][j]);
	f4 = -(exp(gamm + gTheta[c][i][b][j]))*((double)T[i][b][j]);
	f5 = -(pow((gamm - mu_gamma[c][b]), 2.0))/(2.0 * sigma2_gamma[c][b]);

	double f = f1 + f2 + f3 + f4 + f5;

	return(f);
}

void c212BB_poisson_mc_hier2_lev1::sample_gamma_MH(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0, j = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {
				for (j = 0; j < gNAE[l][b]; j++) {

#ifdef INDIVIDUAL_RNG
					GetRNGstate();
#endif
					double cand = rnorm(gGamma[c][l][b][j], gSigma_MH_gamma[l][b][j]);
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

					double f1 = log_f_gamma(c, l, b, j , cand);
					double f2 = log_f_gamma(c, l, b, j , gGamma[c][l][b][j]);

					double ratio = exp(f1 - f2);

					ratio = cMIN(ratio, 1.0);

					if (u <= ratio) {
						gGamma[c][l][b][j] = cand;
						gGamma_acc[c][l][b][j] = gGamma_acc[c][l][b][j] + 1;
					}


					if (iter >= burnin && retainSamples(iMonitor_gamma)) {
						gGamma_samples[c][l][b][j][iter - burnin] = gGamma[c][l][b][j];
					}
				}
			}
		}
	}
}

void c212BB_poisson_mc_hier2_lev1::sample_gamma_SLICE(int burnin, int iter)
{
	int c = 0, i = 0;
	int K = 0, J = 0;

	for (c = 0; c < gChains; c++) {
		for (i = 0; i < gNumIntervals; i++) {

			int b = 0, j = 0;
			double cand = 0.0;

			for (b = 0; b < gNumBodySys[i]; b++) {
				for (j = 0; j < gNAE[i][b]; j++) {

					int m = gW_gamma_control[i][b][j];

#ifdef INDIVIDUAL_RNG
					GetRNGstate();
#endif
					J = floor(runif(0,m));
#ifdef INDIVIDUAL_RNG
					PutRNGstate();
#endif
				    K = (m-1) - J;

					double l = 0.0, r = 0.0;
					double g = log_f_gamma(c, i, b, j, gGamma[c][i][b][j]);
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
					double u = runif(0, gW_gamma[i][b][j]);
#ifdef INDIVIDUAL_RNG
					PutRNGstate();
#endif

					l = gGamma[c][i][b][j] - u;
					r = gGamma[c][i][b][j] + (gW_gamma[i][b][j] - u);

					while (J > 0) {
						if (logy >= log_f_gamma(c, i, b, j, l)) {
							break;
						}
						l = l - gW_gamma[i][b][j];
						J--;
					}

					while (K > 0) {
						if (logy >= log_f_gamma(c, i, b, j, r)) {
							break;
						}
						r = r + gW_gamma[i][b][j];
						K--;
					}
		
#ifdef INDIVIDUAL_RNG
					GetRNGstate();
#endif
					cand = runif(l, r);
#ifdef INDIVIDUAL_RNG
					PutRNGstate();
#endif

					while (logy >= log_f_gamma(c, i, b, j, cand)) {
						if (cand < gGamma[c][i][b][j]) {
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

					gGamma[c][i][b][j] = cand;

					if (iter >= burnin && retainSamples(iMonitor_gamma)) {
						gGamma_samples[c][i][b][j][iter - burnin] = gGamma[c][i][b][j];
					}
				}
			}
		}
	}
}

double c212BB_poisson_mc_hier2_lev1::log_q_theta(int l, int b, int j, double p, double theta, double mean)
{
	double f = 0.0;

	if (theta == 0.0) {
		f = log(p);
	}
	else {
		f = log(1 - p) + log((1.0/(gSigma_MH_theta[l][b][j] * sqrt(2.0 * M_PI))))  + (-1.0/(2.0*gSigma_MH_theta[l][b][j]*gSigma_MH_theta[l][b][j])) * pow((theta - mean), 2.0);
	}

	return f;
}


double c212BB_poisson_mc_hier2_lev1::log_f_theta(int c, int i, int b, int j, double theta)
{
	double f1 = 0.0, f2 = 0.0;

	f1 = (((double)y[i][b][j]) * theta) - (exp(gGamma[c][i][b][j] + theta)) * ((double)T[i][b][j]);

	if (theta == 0.0) {
		f2 = log(gPi[c][b]);
	}
	else {
		f2 = log(1 - gPi[c][b]) + log(1.0/sqrt(2.0 * M_PI*sigma2_theta[c][b]))
				+ ((-1.0/2.0)*(pow(theta -mu_theta[c][b], 2.0))/sigma2_theta[c][b]);
	}

	double f = f1 + f2;

	return(f);
}

/*
* Sample theta using a MH step as detailed in: Gottardo, Raftery - Markov Chain Monte Carlo
* With Mixtures of Mutually Singular Distributions
* Perform an adaption step to get a good candidate density
*/
void c212BB_poisson_mc_hier2_lev1::sample_theta_MH(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0, j = 0;
			for (b = 0; b < gNumBodySys[l]; b++) {
				for ( j = 0; j < gNAE[l][b]; j++) {

#ifdef INDIVIDUAL_RNG
					GetRNGstate();
#endif
					double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
					PutRNGstate();
#endif

					double cand = 0.0;

					if (u < gWp[l][b][j]) {
						cand = 0.0;
					}
					else {
#ifdef INDIVIDUAL_RNG
						GetRNGstate();
#endif
						cand = rnorm(gTheta[c][l][b][j], gSigma_MH_theta[l][b][j]);
#ifdef INDIVIDUAL_RNG
						PutRNGstate();
#endif
					}

					double f_cand = log_f_theta(c, l, b, j, cand);
					double f_prev = log_f_theta(c, l, b, j, gTheta[c][l][b][j]);

					double q_cand = log_q_theta(l, b, j, gWp[l][b][j],
														cand, gTheta[c][l][b][j]);
					double q_prev = log_q_theta(l, b, j, gWp[l][b][j],
														gTheta[c][l][b][j], cand);

					double lratio = f_cand - f_prev + q_prev - q_cand;

					double ratio = exp(lratio);

#ifdef INDIVIDUAL_RNG
					GetRNGstate();
#endif
					u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
					PutRNGstate();
#endif
					if (u <= ratio) {
						gTheta[c][l][b][j] = cand;
						gTheta_acc[c][l][b][j] = gTheta_acc[c][l][b][j] + 1;
						//if (iter >= burnin)
						//	gTheta_acc[c][l][b][j] = gTheta_acc[c][l][b][j] + 1;
					}

					if (iter >= burnin && retainSamples(iMonitor_theta)) {
						gTheta_samples[c][l][b][j][iter - burnin] = gTheta[c][l][b][j];
					}

//					if (cand == 0.0 && gTheta[c][l][b][j] == 0.0) {
//						gTheta[c][l][b][j] = cand;
//						if (iter >= burnin)
//							gTheta_acc[c][l][b][j] = gTheta_acc[c][l][b][j] + 1;
//					}
//					else {
//						if (u <= ratio) {
//							gTheta[c][l][b][j] = cand;
//							if (iter >= burnin)
//								gTheta_acc[c][l][b][j] = gTheta_acc[c][l][b][j] + 1;
//						}
//					}
//
//					if (iter >= burnin && retainSamples(iMonitor_theta)) {
//						gTheta_samples[c][l][b][j][iter - burnin] = gTheta[c][l][b][j];
//					}
				}
			}
		}
	}
}

double c212BB_poisson_mc_hier2_lev1::cMIN(double a, double b)
{
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

void c212BB_poisson_mc_hier2_lev1::release()
{
	releaseL2Variables();

	releaseL2Samples();
}

SEXP c212BB_poisson_mc_hier2_lev1::getL2Samples(double*** &data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = Rf_allocVector(REALSXP, gChains * gMaxBs * (gIter - gBurnin)));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		int b = 0;
		for (b = 0; b < gMaxBs; b++) {
//			if (b < gNumBodySys[l]) {
				memcpy(REAL(samples) + i, data[c][b],
											(gIter - gBurnin)*sizeof(double));
//			}
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
	INTEGER(dim)[1] = gMaxBs;
	INTEGER(dim)[2] = gChains;

	Rf_setAttrib(samples, R_DimSymbol, dim);

	UNPROTECT(2);

	return samples;
}

SEXP c212BB_poisson_mc_hier2_lev1::getMuThetaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(mu_theta_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier2_lev1::getMuGammaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(mu_gamma_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier2_lev1::getSigma2ThetaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(sigma2_theta_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier2_lev1::getSigma2GammaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(sigma2_gamma_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier2_lev1::getPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(gPi_samples);

	return samples;
}


void c212BB_poisson_mc_hier2_lev1::getMuThetaSamples(int *c, int *l, int* b, double* mu_theta)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	if (mu_theta_samples)
		memcpy(mu_theta, mu_theta_samples[C][B], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier2_lev1::getMuGammaSamples(int *c, int *l, int* b, double* mu_gamma)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	if (mu_gamma_samples)
		memcpy(mu_gamma, mu_gamma_samples[C][B], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier2_lev1::getSigma2ThetaSamples(int *c, int *l, int* b, double* sigma2)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	if (sigma2_theta_samples)
		memcpy(sigma2, sigma2_theta_samples[C][B], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier2_lev1::getSigma2GammaSamples(int *c, int *l, int* b, double* sigma2)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	if (sigma2_gamma_samples)
		memcpy(sigma2, sigma2_gamma_samples[C][B], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier2_lev1::getPiSamples(int *c, int *l, int* b, double* pi)
{
	int C = (*c) - 1;
	int B = (*b) - 1;

	if (gPi_samples)
		memcpy(pi, gPi_samples[C][B], (gIter - gBurnin)*sizeof(double));
}
