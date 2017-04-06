#include<cstdio>
#include<cstdlib>

#include<cstring>
#include<cmath>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <Rinternals.h>


#include "c2121a_poisson_mc_hier2_lev0.h"
#include "c2121a_poisson_mc_hier3_lev0.h"
#include "c212BB_poisson_mc_hier3_lev0.h"
#include "c212BB_poisson_mc_hier3_lev2.h"

using namespace std;

static const char *rcsId = "$Id: c212BB_poisson_mc_hier3_lev2.cpp,v 1.16 2017/03/22 16:12:09 clb13102 Exp clb13102 $";

c212BB_poisson_mc_hier3_lev2::c212BB_poisson_mc_hier3_lev2()
{
	//Rprintf("c212BB_poisson_mc_hier3_lev2::c212BB_poisson_mc_hier3_lev2: Default constructor\n");

	mu_theta_0 = NULL;
	mu_gamma_0 = NULL;
	tau2_theta_0 = NULL;
	tau2_gamma_0 = NULL;
	alpha_pi = NULL;
	beta_pi = NULL;

	alpha_pi_acc = NULL;
	beta_pi_acc = NULL;

	mu_theta_0_samples = NULL;
	mu_gamma_0_samples = NULL;
	tau2_theta_0_samples = NULL;
	tau2_gamma_0_samples = NULL;
	alpha_pi_samples = NULL;
	beta_pi_samples = NULL;
}

c212BB_poisson_mc_hier3_lev2::c212BB_poisson_mc_hier3_lev2(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model, SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE, SEXP pX,
					SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma, SEXP pmu_gamma_0_0,
					SEXP ptau2_gamma_0_0, SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0, SEXP palpha_gamma_0_0,
					SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0, SEXP pbeta_theta_0_0, SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta, SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0, SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi, SEXP plambda_alpha, SEXP plambda_beta,
					SEXP palgo, SEXP padapt_phase)

{
	mu_theta_0 = NULL;
	mu_gamma_0 = NULL;
	tau2_theta_0 = NULL;
	tau2_gamma_0 = NULL;
	alpha_pi = NULL;
	beta_pi = NULL;

	alpha_pi_acc = NULL;
	beta_pi_acc = NULL;

	mu_theta_0_samples = NULL;
	mu_gamma_0_samples = NULL;
	tau2_theta_0_samples = NULL;
	tau2_gamma_0_samples = NULL;
	alpha_pi_samples = NULL;
	beta_pi_samples = NULL;

	init(sChains, sBurnin, sIter, sSim_Type, sMem_Model, sGlobal_Sim_Params,
				sSim_Params, MH_weight, pm_weights,
				sMonitor,
				sNumIntervals, sMaxBs, sNumBodySys, sMaxAEs, sNAE,
				pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0, pmu_theta_0_0, ptau2_theta_0_0,
				palpha_gamma_0_0, pbeta_gamma_0_0, palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma,
				pbeta_gamma, palpha_theta, pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0,
				ptau2_theta_0, pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
				pPi, palpha_pi, pbeta_pi, plambda_alpha, plambda_beta,
				palgo, padapt_phase);
}

void c212BB_poisson_mc_hier3_lev2::init(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type, SEXP sMem_Model, SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP sMonitor,
					SEXP sNumIntervals,
					SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP pmu_gamma_0_0,
					SEXP ptau2_gamma_0_0, SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0, SEXP palpha_gamma_0_0,
					SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0, SEXP pbeta_theta_0_0, SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta, SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0, SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi, SEXP plambda_alpha, SEXP plambda_beta,
					SEXP palgo, SEXP padapt_phase)

{
	release();
	c212BB_poisson_mc_hier3_lev0::release();
	c2121a_poisson_mc_hier3_lev0::release();
	c2121a_poisson_mc_hier2_lev0::release();

	initMonitor(sMonitor);

	int l = 0, b = 0, c = 0, a = 0;

	gChains = *(INTEGER(sChains));
	gBurnin = *(INTEGER(sBurnin));
	gIter = *(INTEGER(sIter));
	gNumIntervals = *(INTEGER(sNumIntervals));
	gMaxBs = *(INTEGER(sMaxBs));

	// Body-system Data
	gNumBodySys = (int *)malloc(gNumIntervals * sizeof(int));
	for (l = 0; l < gNumIntervals; l++) {
		gNumBodySys[l] = (INTEGER(sNumBodySys))[l];
	}

	// AE data
	gMaxAEs = *(INTEGER(sMaxAEs));
	gNAE = (int **)malloc(gNumIntervals * sizeof(int *));
	for (l = 0; l < gNumIntervals; l++) {
		gNAE[l] = (int *)malloc(gMaxBs * sizeof(int));
	}

	int indx = 0;
	for (l = 0; l < gNumIntervals; l++) {
		for (b = 0; b < gMaxBs; b++) {
			gNAE[l][b] = (INTEGER(sNAE))[indx++];
		}
	}

	l = strlen(CHAR(STRING_ELT(sMem_Model, 0)));
	char *mem_model = (char *)malloc((l + 1)*sizeof(char));
	if (mem_model) {
		strcpy(mem_model, CHAR(STRING_ELT(sMem_Model, 0)));
		mem_model[l] = 0;

		Rprintf("Memory Model: %s\n", mem_model);

		if (0 == strcmp("LOW", mem_model)) {
			eMemory_Model = LOW;
		}
		else {
			eMemory_Model = HIGH;
		}

		free(mem_model);
		mem_model = NULL;
	}

	x = (int***)malloc(gNumIntervals * sizeof(int**));
	y = (int***)malloc(gNumIntervals * sizeof(int**));
	C = (int***)malloc(gNumIntervals * sizeof(int**));
	T = (int***)malloc(gNumIntervals * sizeof(int**));
	for (l = 0; l < gNumIntervals; l++) {
		x[l] = (int**)malloc(gMaxBs * sizeof(int*));
		y[l] = (int**)malloc(gMaxBs * sizeof(int*));
		C[l] = (int**)malloc(gMaxBs * sizeof(int*));
		T[l] = (int**)malloc(gMaxBs * sizeof(int*));
		for (b = 0; b < gMaxBs; b++) {
			x[l][b] = (int*)malloc(gMaxAEs * sizeof(int));
			y[l][b] = (int*)malloc(gMaxAEs * sizeof(int));
			C[l][b] = (int*)malloc(gMaxAEs * sizeof(int));
			T[l][b] = (int*)malloc(gMaxAEs * sizeof(int));
		}
	}

	int *vX = INTEGER(pX);
	int *vY = INTEGER(pY);
	int *vC = INTEGER(pC);
	int *vT = INTEGER(pT);
	for (l = 0; l < gNumIntervals; l++) {
		int b = 0;
		for (b = 0; b < gMaxBs; b++) {
			for (a = 0; a < gMaxAEs; a++) {
				x[l][b][a] = *vX;
				y[l][b][a] = *vY;
				C[l][b][a] = *vC;
				T[l][b][a] = *vT;
				vX++;
				vY++;
				vC++;
				vT++;
			}
		}
	}

	gTheta = (double****)malloc(gChains * sizeof(double***));
	gGamma = (double****)malloc(gChains * sizeof(double***));
	for (c = 0; c < gChains; c++) {
		gTheta[c] = (double***)malloc(gNumIntervals * sizeof(double**));
		gGamma[c] = (double***)malloc(gNumIntervals * sizeof(double**));
		for (l = 0; l < gNumIntervals; l++) {
			gTheta[c][l] = (double**)malloc(gMaxBs * sizeof(double*));
			gGamma[c][l] = (double**)malloc(gMaxBs * sizeof(double*));
			for (b = 0; b < gMaxBs; b++) {
				gTheta[c][l][b] = (double*)malloc(gMaxAEs * sizeof(double));
				gGamma[c][l][b] = (double*)malloc(gMaxAEs * sizeof(double));
			}
		}
	}

	double* vtheta = REAL(ptheta);
	double* vgamma = REAL(pgamma);
	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gMaxBs; b++) {
				for (a = 0; a < gMaxAEs; a++) {
					gTheta[c][l][b][a] = *vtheta;
					gGamma[c][l][b][a] = *vgamma;
					vtheta++;
					vgamma++;
				}
			}
		}
	}

	mu_gamma_0_0 = *(REAL(pmu_gamma_0_0));
	tau2_gamma_0_0 = *(REAL(ptau2_gamma_0_0));
	mu_theta_0_0 = *(REAL(pmu_theta_0_0));
	tau2_theta_0_0 = *(REAL(ptau2_theta_0_0));
	alpha_gamma_0_0 = *(REAL(palpha_gamma_0_0));
	beta_gamma_0_0 = *(REAL(pbeta_gamma_0_0));
	alpha_theta_0_0 = *(REAL(palpha_theta_0_0));
	beta_theta_0_0 = *(REAL(pbeta_theta_0_0));
	alpha_gamma = *(REAL(palpha_gamma));
	beta_gamma = *(REAL(pbeta_gamma));
	alpha_theta = *(REAL(palpha_theta));
	beta_theta = *(REAL(pbeta_theta));
	lambda_alpha = *(REAL(plambda_alpha));
	lambda_beta = *(REAL(plambda_beta));

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

	alpha_pi = (double*)malloc(gChains * sizeof(double));
	double *valpha_pi = REAL(palpha_pi);
	for (c = 0; c < gChains; c++) {
		alpha_pi[c] = *valpha_pi;
		valpha_pi++;
	}

	beta_pi = (double*)malloc(gChains * sizeof(double));
	double *vbeta_pi = REAL(pbeta_pi);
	for (c = 0; c < gChains; c++) {
		beta_pi[c] = *vbeta_pi;
		vbeta_pi++;
	}

	double* vmu_gamma = REAL(pmu_gamma);
	mu_gamma = (double***)malloc(gChains * sizeof(double**));
	for (c = 0; c < gChains; c++) {
		mu_gamma[c] = (double**)malloc(gNumIntervals * sizeof(double*));
		for (l = 0; l < gNumIntervals; l++) {
			mu_gamma[c][l] = (double*)malloc(gMaxBs * sizeof(double));
			for (b = 0; b < gMaxBs; b++) {
				mu_gamma[c][l][b] = *vmu_gamma;
				vmu_gamma++;
			}
		}
	}

	double* vmu_theta = REAL(pmu_theta);
	mu_theta = (double***)malloc(gChains * sizeof(double**));
	for (c = 0; c < gChains; c++) {
		mu_theta[c] = (double**)malloc(gNumIntervals * sizeof(double*));
		for (l = 0; l < gNumIntervals; l++) {
			mu_theta[c][l] = (double*)malloc(gMaxBs * sizeof(double));
			for (b = 0; b < gMaxBs; b++) {
				mu_theta[c][l][b] = *vmu_theta;
				vmu_theta++;
			}
		}
	}

	double* vsigma2_gamma = REAL(psigma2_gamma);
	sigma2_gamma = (double***)malloc(gChains * sizeof(double**));
	for (c = 0; c < gChains; c++) {
		sigma2_gamma[c] = (double**)malloc(gNumIntervals * sizeof(double*));
		for (l = 0; l < gNumIntervals; l++) {
			sigma2_gamma[c][l] = (double*)malloc(gMaxBs * sizeof(double));
			for (b = 0; b < gMaxBs; b++) {
				sigma2_gamma[c][l][b] = *vsigma2_gamma;
				vsigma2_gamma++;
			}
		}
	}

	double* vsigma2_theta = REAL(psigma2_theta);
	sigma2_theta = (double***)malloc(gChains * sizeof(double**));
	for (c = 0; c < gChains; c++) {
		sigma2_theta[c] = (double**)malloc(gNumIntervals * sizeof(double*));
		for (l = 0; l < gNumIntervals; l++) {
			sigma2_theta[c][l] = (double*)malloc(gMaxBs * sizeof(double));
			for (b = 0; b < gMaxBs; b++) {
				sigma2_theta[c][l][b] = *vsigma2_theta;
				vsigma2_theta++;
			}
		}
	}

	double* vpi = REAL(pPi);
	gPi = (double***)malloc(gChains * sizeof(double**));
	for (c = 0; c < gChains; c++) {
		gPi[c] = (double**)malloc(gNumIntervals * sizeof(double*));
		for (l = 0; l < gNumIntervals; l++) {
			gPi[c][l] = (double*)malloc(gMaxBs * sizeof(double));
			for (b = 0; b < gMaxBs; b++) {
				gPi[c][l][b] = *vpi;
				vpi++;
			}
		}
	}

	// The samples
	if (retainSamples(iMonitor_mu_gamma_0))
		mu_gamma_0_samples = (double **)malloc(gChains * sizeof(double*));
	if (retainSamples(iMonitor_mu_theta_0))
		mu_theta_0_samples = (double **)malloc(gChains * sizeof(double*));
	if (retainSamples(iMonitor_tau2_gamma_0))
		tau2_gamma_0_samples = (double **)malloc(gChains * sizeof(double*));
	if (retainSamples(iMonitor_tau2_theta_0))
		tau2_theta_0_samples = (double **)malloc(gChains * sizeof(double*));
	if (retainSamples(iMonitor_alpha_pi))
		alpha_pi_samples = (double **)malloc(gChains * sizeof(double*));
	if (retainSamples(iMonitor_beta_pi))
		beta_pi_samples = (double **)malloc(gChains * sizeof(double*));

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
		if (retainSamples(iMonitor_alpha_pi))
			alpha_pi_samples[c] =
									(double *)malloc((gIter - gBurnin)* sizeof(double));
		if (retainSamples(iMonitor_beta_pi))
			beta_pi_samples[c] =
									(double *)malloc((gIter - gBurnin)* sizeof(double));
	}

	if (retainSamples(iMonitor_mu_theta))
		mu_theta_samples = (double ****)malloc(gChains *sizeof(double***));
	if (retainSamples(iMonitor_mu_gamma))
		mu_gamma_samples = (double ****)malloc(gChains *sizeof(double***));
	if (retainSamples(iMonitor_sigma2_theta))
		sigma2_theta_samples = (double ****)malloc(gChains *sizeof(double***));
	if (retainSamples(iMonitor_sigma2_gamma))
		sigma2_gamma_samples = (double ****)malloc(gChains *sizeof(double***));
	if (retainSamples(iMonitor_pi))
		gPi_samples = (double ****)malloc(gChains *sizeof(double***));

	for (c = 0; c < gChains; c++) {
		if (retainSamples(iMonitor_mu_theta))
			mu_theta_samples[c] = (double ***)malloc(gNumIntervals *sizeof(double**));
		if (retainSamples(iMonitor_mu_gamma))
			mu_gamma_samples[c] = (double ***)malloc(gNumIntervals *sizeof(double**));
		if (retainSamples(iMonitor_sigma2_theta))
			sigma2_theta_samples[c] =
								(double ***)malloc(gNumIntervals *sizeof(double**));
		if (retainSamples(iMonitor_sigma2_gamma))
			sigma2_gamma_samples[c] =
								(double ***)malloc(gNumIntervals *sizeof(double**));
		if (retainSamples(iMonitor_pi))
			gPi_samples[c] = (double ***)malloc(gNumIntervals *sizeof(double**));

		for (l = 0; l < gNumIntervals; l++) {
			if (retainSamples(iMonitor_mu_theta))
				mu_theta_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));
			if (retainSamples(iMonitor_mu_gamma))
				mu_gamma_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));
			if (retainSamples(iMonitor_sigma2_theta))
				sigma2_theta_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));
			if (retainSamples(iMonitor_sigma2_gamma))
				sigma2_gamma_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));
			if (retainSamples(iMonitor_pi))
				gPi_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));

			for (b = 0; b < gNumBodySys[l]; b++) {
				if (retainSamples(iMonitor_mu_theta))
					mu_theta_samples[c][l][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));
				if (retainSamples(iMonitor_mu_gamma))
					mu_gamma_samples[c][l][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));
				if (retainSamples(iMonitor_sigma2_theta))
				sigma2_theta_samples[c][l][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));
				if (retainSamples(iMonitor_sigma2_gamma))
					sigma2_gamma_samples[c][l][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));
				if (retainSamples(iMonitor_pi))
					gPi_samples[c][l][b] =
									(double *)malloc((gIter - gBurnin)*sizeof(double));
			}
		}
	}

	if (retainSamples(iMonitor_theta))
		gTheta_samples = (double *****)malloc(gChains*sizeof(double****));
	if (retainSamples(iMonitor_gamma))
		gGamma_samples = (double *****)malloc(gChains*sizeof(double****));

	for (c = 0; c < gChains; c++) {
		if (retainSamples(iMonitor_theta))
			gTheta_samples[c] = (double ****)malloc(gNumIntervals *sizeof(double***));
		if (retainSamples(iMonitor_gamma))
			gGamma_samples[c] = (double ****)malloc(gNumIntervals *sizeof(double***));

		for (l = 0; l < gNumIntervals; l++) {
			if (retainSamples(iMonitor_theta))
				gTheta_samples[c][l] =
								(double ***)malloc(gNumBodySys[l] *sizeof(double**));
			if (retainSamples(iMonitor_gamma))
				gGamma_samples[c][l] =
								(double ***)malloc(gNumBodySys[l] *sizeof(double**));

			for (b = 0; b < gNumBodySys[l]; b++) {
				if (retainSamples(iMonitor_theta))
					gTheta_samples[c][l][b] =
									(double **)malloc(gNAE[l][b] *sizeof(double*));
				if (retainSamples(iMonitor_gamma))
					gGamma_samples[c][l][b] =
									(double **)malloc(gNAE[l][b] *sizeof(double*));

				for (a = 0; a < gNAE[l][b]; a++) {
					if (retainSamples(iMonitor_theta))
						gTheta_samples[c][l][b][a] =
									(double *)malloc((gIter - gBurnin) *sizeof(double));
					if (retainSamples(iMonitor_gamma))
						gGamma_samples[c][l][b][a] =
									(double *)malloc((gIter - gBurnin) *sizeof(double));
				}
			}
		}
	}

	gTheta_acc = (int****)malloc(gChains * sizeof(int***));
	gGamma_acc = (int****)malloc(gChains * sizeof(int***));
	for (c = 0; c < gChains; c++) {
		gTheta_acc[c] = (int***)malloc(gNumIntervals * sizeof(int**));
		gGamma_acc[c] = (int***)malloc(gNumIntervals * sizeof(int**));
		for (l = 0; l < gNumIntervals; l++) {
			gTheta_acc[c][l] = (int**)malloc(gMaxBs * sizeof(int*));
			gGamma_acc[c][l] = (int**)malloc(gMaxBs * sizeof(int*));
			for (b = 0; b < gMaxBs; b++) {
				gTheta_acc[c][l][b] = (int*)malloc(gMaxAEs * sizeof(int));
				gGamma_acc[c][l][b] = (int*)malloc(gMaxAEs * sizeof(int));
				for (a = 0; a < gMaxAEs; a++) {
					gTheta_acc[c][l][b][a] = 0;
					gGamma_acc[c][l][b][a] = 0;
				}
			}
		}
	}

	alpha_pi_acc = (int*)malloc(gChains * sizeof(int));
	beta_pi_acc = (int*)malloc(gChains * sizeof(int));
	for (c = 0; c < gChains; c++) {
		alpha_pi_acc[c] = 0;
		beta_pi_acc[c] = 0;
	}

	// Global simulation parameters
	initGlobalSimParams(sSim_Type, sGlobal_Sim_Params);

	// Individual simulation parameters
	initSimParams(sSim_Params);

	// MH point-mass weights
	gMH_weight = *(REAL(MH_weight));
	initPMWeights(pm_weights);

	//Rprintf("Simultion Data: %s (%0.6f) %d, %d, %d\n",
	//							sim_type, gSim_Param, gChains, gBurnin, gIter);
	//
	//Rprintf("Intervals: %d\n", gNumIntervals);
	//for (l = 0; l < gNumIntervals; l++)
	//	Rprintf("\tInterval %d: Contains %d body-systems\n", l, gNumBodySys[l]);
	//
 	//Rprintf("\tMaxBs: %d\n", gMaxBs);
 	//Rprintf("\tMaxAEs: %d\n", gMaxAEs);
	//
	//for (l = 0; l < gNumIntervals; l++) {
	//	for (b = 0; b < gMaxBs; b++) {
	//		Rprintf("\tAEs in Interval: %d, BS: %d = %d\n", l, b, gNAE[l][b]);
	//	}
	//}
	//
	//Rprintf("Control Count Data:\n");
	//for (l = 0; l < gNumIntervals; l++) {
	//	Rprintf("\tInterval: %d\n", l);
	//	int b = 0;
	//	for (b = 0; b < gNumBodySys[l]; b++) {
	//		Rprintf("\t\tBody-system: %d\n", b);
	//		int a = 0;
	//		for (a = 0; a < gNAE[l][b]; a++) {
	//			Rprintf("\t\t\tAE: %d - Count: %d\n", a, x[l][b][a]);
	//		}
	//	}
	//}
	//
	//Rprintf("Treatment Count Data:\n");
	//for (l = 0; l < gNumIntervals; l++) {
	//	Rprintf("\tInterval: %d\n", l);
	//	for (b = 0; b < gNumBodySys[l]; b++) {
	//		Rprintf("\t\tBody-system: %d\n", b);
	//		int a = 0;
	//		for (a = 0; a < gNAE[l][b]; a++) {
	//			Rprintf("\tAE: %d - Count: %d\n", a, y[l][b][a]);
	//		}
	//	}
	//}
	//
	//Rprintf("Control Exposure Data:\n");
	//for (l = 0; l < gNumIntervals; l++) {
	//	Rprintf("\tInterval: %d\n", l);
	//	for (b = 0; b < gNumBodySys[l]; b++) {
	//		Rprintf("\t\tBody-system: %d\n", b);
	//		int a = 0;
	//		for (a = 0; a < gNAE[l][b]; a++) {
	//			Rprintf("\t\t\tAE: %d - Count: %d\n", a, C[l][b][a]);
	//		}
	//	}
	//}
	//
	//Rprintf("Treatment Exposure Data:\n");
	//for (l = 0; l < gNumIntervals; l++) {
	//	Rprintf("\tInterval: %d\n", l);
	//	for (b = 0; b < gNumBodySys[l]; b++) {
	//		Rprintf("\t\tBody-system: %d\n", b);
	//		int a = 0;
	//		for (a = 0; a < gNAE[l][b]; a++) {
	//			Rprintf("\t\t\tAE: %d - Count: %d\n", a, T[l][b][a]);
	//		}
	//	}
	//}
	//
	//Rprintf("Theta initialised values:\n");
	//for (c = 0; c < gChains; c++) {
	//	Rprintf("\tChain: %d\n", c);
	//	for (l = 0; l < gNumIntervals; l++) {
	//		Rprintf("\tInterval: %d\n", l);
	//		for (b = 0; b < gNumBodySys[l]; b++) {
	//			Rprintf("\t\tBody-system: %d\n", b);
	//			int a = 0;
	//			for (a = 0; a < gNAE[l][b]; a++) {
	//				Rprintf("\t\t\t\t: %0.6f\n", gTheta[c][l][b][a]);
	//			}
	//		}
	//	}
	//}
	//
	//Rprintf("Gamma initialised values:\n");
	//for (c = 0; c < gChains; c++) {
	//	Rprintf("\tChain: %d\n", c);
	//	for (l = 0; l < gNumIntervals; l++) {
	//		Rprintf("\tInterval: %d\n", l);
	//		for (b = 0; b < gNumBodySys[l]; b++) {
	//			Rprintf("\t\tBody-system: %d\n", b);
	//			int a = 0;
	//			for (a = 0; a < gNAE[l][b]; a++) {
	//				Rprintf("\t\t\t\t: %0.6f\n", gGamma[c][l][b][a]);
	//			}
	//		}
	//	}
	//}
	//
	//Rprintf("mu.gamma.0 initialised values:\n");
	//for (c = 0; c < gChains; c++) {
	//	for (l = 0; l < gNumIntervals; l++) {
	//		Rprintf("\t%0.6f\n", mu_gamma_0[c][l]);
	//	}
	//}
	//
	//Rprintf("mu.theta.0 initialised values:\n");
	//for (c = 0; c < gChains; c++) {
	//	for (l = 0; l < gNumIntervals; l++) {
	//		Rprintf("\t%0.6f\n", mu_theta_0[c][l]);
	//	}
	//}
	//
	//Rprintf("tau2.gamma.0 initialised values:\n");
	//for (c = 0; c < gChains; c++) {
	//	for (l = 0; l < gNumIntervals; l++) {
	//		Rprintf("\t%0.6f\n", tau2_gamma_0[c][l]);
	//	}
	//}
	//
	//Rprintf("tau2.theta.0 initialised values:\n");
	//for (c = 0; c < gChains; c++) {
	//	for (l = 0; l < gNumIntervals; l++) {
	//		Rprintf("\t%0.6f\n", tau2_theta_0[c][l]);
	//	}
	//}
	//
	//Rprintf("sigma2_gamma initialised values:\n");
	//for (c = 0; c < gChains; c++) {
	//	for (l = 0; l < gNumIntervals; l++) {
	//		for (b = 0; b < gNumBodySys[l]; b++) {
	//			Rprintf("\t%d %d %d: %0.6f\n", c, l , b, sigma2_gamma[c][l][b]);
	//		}
	//	}
	//}
	//
	//Rprintf("mu.gamma.0.0: %0.6f\n", mu_gamma_0_0);
	//Rprintf("tau2.gamma.0.0: %0.6f\n", tau2_gamma_0_0);
	//Rprintf("mu.theta.0.0: %0.6f\n", mu_theta_0_0);
	//Rprintf("tau2.theta.0.0: %0.6f\n", tau2_theta_0_0);
	//Rprintf("alpha.gamma.0.0: %0.6f\n", alpha_gamma_0_0);
	//Rprintf("beta.gamma.0.0: %0.6f\n", beta_gamma_0_0);
	//Rprintf("alpha.theta.0.0: %0.6f\n", alpha_theta_0_0);
	//Rprintf("beta.theta.0.0: %0.6f\n", beta_theta_0_0);
	//Rprintf("alpha.gamma: %0.6f\n", alpha_gamma);
	//Rprintf("beta.gamma: %0.6f\n", beta_gamma);
	//Rprintf("alpha.theta: %0.6f\n", alpha_theta);
	//Rprintf("beta.theta: %0.6f\n", beta_theta);
}

c212BB_poisson_mc_hier3_lev2::~c212BB_poisson_mc_hier3_lev2()
{
	//Rprintf("c212BB_poisson_mc_hier3_lev2::c212BB_poisson_mc_hier3_lev2 - destructor\n");
	release();
}

void c212BB_poisson_mc_hier3_lev2::gibbs_sampler()
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

void c212BB_poisson_mc_hier3_lev2::simulate_MH()
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

		sample_alpha_pi_MH(gBurnin, i);
		sample_beta_pi_MH(gBurnin, i);
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

void c212BB_poisson_mc_hier3_lev2::simulate_SLICE()
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

		sample_alpha_pi_SLICE(gBurnin, i);
		sample_beta_pi_SLICE(gBurnin, i);
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

void c212BB_poisson_mc_hier3_lev2::sample_mu_gamma_0(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		int sum_B = 0;
		double mu_gamma_tot = 0.0;

		for (l = 0; l < gNumIntervals; l++) {
			sum_B = sum_B + gNumBodySys[l];

			int b = 0;
			for (b = 0; b < gNumBodySys[l]; b++) {
				mu_gamma_tot += mu_gamma[c][l][b];
			}
		}

		double denom = tau2_gamma_0[c] + tau2_gamma_0_0 * ((double)sum_B);

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

		if (iter >= burnin && retainSamples(iMonitor_mu_gamma_0)) {
			mu_gamma_0_samples[c][iter - burnin] = mu_gamma_0[c];
		}
	}
}

void c212BB_poisson_mc_hier3_lev2::sample_mu_theta_0 (int burnin, int iter)
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

		if (iter >= burnin && retainSamples(iMonitor_mu_theta_0)) {
			mu_theta_0_samples[c][iter - burnin] = mu_theta_0[c];
		}
	}
}

void c212BB_poisson_mc_hier3_lev2::sample_tau2_gamma_0(int burnin, int iter)
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

void c212BB_poisson_mc_hier3_lev2::sample_tau2_theta_0(int burnin, int iter)
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

double c212BB_poisson_mc_hier3_lev2::log_f_alpha_pi(int c, double alpha)
{
	double f = 0.0;
	int l = 0;

	for (l = 0; l < gNumIntervals; l++) {
		int b = 0;
		double log_pi_sum = 0.0;
		for (b = 0; b < gNumBodySys[l]; b++) {
			log_pi_sum += log(gPi[c][l][b]);
		}

		f = f + ((double)gNumBodySys[l]) * (lgammafn(alpha + beta_pi[c]) - lgammafn(alpha));

		f = f + (alpha - 1.0)*log_pi_sum;
	}

	f = f - alpha * lambda_alpha;

	return(f);
}

void c212BB_poisson_mc_hier3_lev2::sample_alpha_pi_MH(int burnin, int iter)
{
	int c = 0;

	for (c = 0; c< gChains; c++) {

		double cand = 0;

		// alpha_pi is restricted to being greater than zero
		// This is rejection sampling of a normal distribution truncated at 1.
		// See c212BB.cpp
	    while (cand <= 1.0) {
#ifdef INDIVIDUAL_RNG
	        GetRNGstate();
#endif
	        cand = rnorm(alpha_pi[c], gDefault_Sigma_MH_alpha);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif
		}

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		double f1 = log_f_alpha_pi(c, cand);
		double f2 = log_f_alpha_pi(c, alpha_pi[c]);

		double q1 = pnorm((alpha_pi[c] - 1)/gDefault_Sigma_MH_alpha, 0, 1, 1, 0);
		double q2 = pnorm((cand - 1)/gDefault_Sigma_MH_alpha, 0, 1, 1, 0);

		double ratio = (exp(f1 - f2)) * q1/q2;
		ratio = cMIN(ratio, 1);

	    if (u <= ratio) {
	        alpha_pi[c] = cand;
			//if (iter >= burnin)
				alpha_pi_acc[c] = alpha_pi_acc[c] + 1;
		}

		if (iter >= burnin && retainSamples(iMonitor_alpha_pi)) {
			alpha_pi_samples[c][iter - burnin] = alpha_pi[c];
		}
	}
}

void c212BB_poisson_mc_hier3_lev2::sample_alpha_pi_SLICE(int burnin, int iter)
{
	int c = 0;
	int m = gDefault_W_alpha_control, K = 0, J = 0;

	for (c = 0; c < gChains; c++) {

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		J = floor(runif(0,m));
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif
	    K = (m-1) - J;

		double cand = 0.0;
		double l = 0.0, r = 0.0;

		double g = log_f_alpha_pi(c, alpha_pi[c]);
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
		double u = runif(0, gDefault_W_alpha);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		l = alpha_pi[c] - u;
		r = alpha_pi[c] + (gDefault_W_alpha - u);

		while (J > 0) {
			if (l <= 1.0)
				break;

			if (logy >= log_f_alpha_pi(c, l)) {
				break;
			}
			l = l - gDefault_W_alpha;

			J--;
		}

		while (K > 0) {
			if (logy >= log_f_alpha_pi(c, r)) {
				break;
			}
			r = r + gDefault_W_alpha;
			K--;
		}

		if (l <= 1.0) {
			l = 1.0;
		}

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		cand = runif(l, r);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		while (logy >= log_f_alpha_pi(c, cand)) {
			if (cand < alpha_pi[c]) {
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

		alpha_pi[c] = cand;

		if (iter >= burnin && retainSamples(iMonitor_alpha_pi)) {
			alpha_pi_samples[c][iter - burnin] = alpha_pi[c];
		}
	}
}

double c212BB_poisson_mc_hier3_lev2::log_f_beta_pi(int c, double beta)
{
	double f = 0.0;
	int l = 0;

	for (l = 0; l < gNumIntervals; l++) {
		int b = 0;
		double log_sum = 0.0;
		for (b = 0; b < gNumBodySys[l]; b++) {
			log_sum += log(1 - gPi[c][l][b]);
		}

		f = f + ((double)gNumBodySys[l]) * (lgammafn(alpha_pi[c] + beta) - lgammafn(beta));

		f = f + (beta - 1.0)*log_sum;
	}

	f = f - beta * lambda_alpha;

	return(f);
}

void c212BB_poisson_mc_hier3_lev2::sample_beta_pi_MH(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {

		double cand = 0.0;

		while (cand <= 1.0) {
#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			cand = rnorm(beta_pi[l], gDefault_Sigma_MH_beta);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif
		}

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		double f1 = log_f_beta_pi(c, cand);
		double f2 = log_f_beta_pi(c, beta_pi[c]);

		double q1 = pnorm((beta_pi[c] - 1)/gDefault_Sigma_MH_beta, 0, 1, 1, 0);
		double q2 = pnorm((cand - 1)/gDefault_Sigma_MH_beta, 0, 1, 1, 0);

		double ratio = (exp(f1 - f2)) * (q1/q2);

		ratio = cMIN(ratio, 1);

		if (u <= ratio) {
			beta_pi[c] = cand;
			//if (iter >= burnin)
				beta_pi_acc[c] = beta_pi_acc[c] + 1;
		}

		if (iter >= burnin && retainSamples(iMonitor_beta_pi)) {
			beta_pi_samples[c][iter - burnin] = beta_pi[c];
		}
	}
}

void c212BB_poisson_mc_hier3_lev2::sample_beta_pi_SLICE(int burnin, int iter)
{
	int c = 0;
	int m = gDefault_W_beta_control, K = 0, J = 0;

	for (c = 0; c < gChains; c++) {

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		J = floor(runif(0,m));
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif
		K = (m-1) - J;

		double l = 0.0, r = 0.0;
		double cand = 0.0;

		double g = log_f_beta_pi(c, beta_pi[c]);
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
		double u = runif(0, gDefault_W_beta);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		l = beta_pi[c] - u;
		r = beta_pi[c] + (gDefault_W_beta - u);


		// beta is retricted to being greater than 1
		// need a do - while loop
		while (J > 0) {
			if (l <= 1.0)
				break;

			if (logy >= log_f_beta_pi(c, l)) {
				break;
			}
			l = l - gDefault_W_beta;
			J--;
		}

		while (K > 0) {
			if (logy >= log_f_beta_pi(c, r)) {
				break;
			}
			r = r + gDefault_W_beta;
			K--;
		}

		if (l <= 1.0) {
			l = 1.0;
		}


#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		cand = runif(l, r);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

		while (logy >= log_f_beta_pi(c, cand)) {
			if (cand < beta_pi[c]) {
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

		beta_pi[c] = cand;

		if (iter >= burnin && retainSamples(iMonitor_beta_pi)) {
			beta_pi_samples[c][iter - burnin] = beta_pi[c];
		}
	}
}

void c212BB_poisson_mc_hier3_lev2::sample_pi(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;
			for (b = 0; b < gNumBodySys[l]; b++) {
				int theta_zero_count = 0;

				int j = 0;
				for (j = 0; j< gNAE[l][b]; j++) {
					if (gTheta[c][l][b][j] == 0.0) {
						theta_zero_count++;
					}
				}

				double shape1 = alpha_pi[c] + (double)theta_zero_count;
				double shape2 = beta_pi[c] + (double)gNAE[l][b] - (double)theta_zero_count;

#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				gPi[c][l][b] = rbeta(shape1, shape2);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif

				if (iter >= burnin && retainSamples(iMonitor_pi)) {
					gPi_samples[c][l][b][iter - burnin] = gPi[c][l][b];
				}
			}
		}

	}
}

void c212BB_poisson_mc_hier3_lev2::sample_mu_gamma(int burnin, int iter)
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

void c212BB_poisson_mc_hier3_lev2::sample_mu_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {

				double t = 0.0;
				int Kb = 0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					if (gTheta[c][l][b][j] != 0.0) {
						Kb++;
					}
					t += gTheta[c][l][b][j];
				}

				double denom = sigma2_theta[c][l][b] + ((double)Kb)*tau2_theta_0[c];

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

void c212BB_poisson_mc_hier3_lev2::sample_sigma2_gamma(int burnin, int iter)
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

void c212BB_poisson_mc_hier3_lev2::sample_sigma2_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {


				double t = 0;
				int j = 0;
				int Kb = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					if (gTheta[c][l][b][j] != 0.0) {
						Kb++;
						t += (pow((gTheta[c][l][b][j] - mu_theta[c][l][b]),2.0));
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

				sigma2_theta[c][l][b] = 1/cand;

				if (iter >= burnin && retainSamples(iMonitor_sigma2_theta)) {
					sigma2_theta_samples[c][l][b][iter - burnin] = sigma2_theta[c][l][b];
				}
			}
		}
	}
}

double c212BB_poisson_mc_hier3_lev2::cMIN(double a, double b)
{
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

void c212BB_poisson_mc_hier3_lev2::release()
{
	int c = 0;

	if (alpha_pi != NULL) {
		free(alpha_pi);
		alpha_pi = NULL;
	}

	if (alpha_pi_acc != NULL) {
		free(alpha_pi_acc);
		alpha_pi_acc = NULL;
	}

	if (beta_pi != NULL) {
		free(beta_pi);
		beta_pi = NULL;
	}

	if (beta_pi_acc != NULL) {
		free(beta_pi_acc);
		beta_pi_acc = NULL;
	}

	if (alpha_pi_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			free(alpha_pi_samples[c]);
		}
		free(alpha_pi_samples);
		alpha_pi_samples = NULL;
	}

	if (beta_pi_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			free(beta_pi_samples[c]);
		}
		free(beta_pi_samples);
		beta_pi_samples = NULL;
	}

	if (mu_theta_0 != NULL) {
		free(mu_theta_0);
		mu_theta_0 = NULL;
	}

	if (mu_gamma_0 != NULL) {
		free(mu_gamma_0);
		mu_gamma_0 = NULL;
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

SEXP c212BB_poisson_mc_hier3_lev2::getL3Samples(double** &data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = allocVector(REALSXP, gChains * (gIter - gBurnin)));

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

	PROTECT(dim = allocVector(INTSXP, 2));

	INTEGER(dim)[0] = (gIter - gBurnin);
	INTEGER(dim)[1] = gChains;

	setAttrib(samples, R_DimSymbol, dim);

	UNPROTECT(2);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev2::getMuGamma0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(mu_gamma_0_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev2::getMuTheta0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(mu_theta_0_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev2::getTau2Gamma0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(tau2_gamma_0_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev2::getTau2Theta0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(tau2_theta_0_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev2::getAlphaPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(alpha_pi_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev2::getBetaPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(beta_pi_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev2::getL3Accept(int* &data)
{
	SEXP acc = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(acc = allocVector(INTSXP, gChains));
	memcpy(INTEGER(acc), data, gChains*sizeof(int));

	free(data);
	data = NULL;

	PROTECT(dim = allocVector(INTSXP, 1));

	INTEGER(dim)[0] = gChains;
	setAttrib(acc, R_DimSymbol, dim);

	UNPROTECT(2);

	return acc;
}

SEXP c212BB_poisson_mc_hier3_lev2::getAlphaPiAccept()
{
	SEXP acc = R_NilValue;

	acc = getL3Accept(alpha_pi_acc);

	return acc;
}

SEXP c212BB_poisson_mc_hier3_lev2::getBetaPiAccept()
{
	SEXP acc = R_NilValue;

	acc = getL3Accept(beta_pi_acc);

	return acc;
}

void c212BB_poisson_mc_hier3_lev2::getMuGamma0Samples(int *c, int *l, double* mu)
{
	int C = (*c) - 1;

	if (mu_gamma_0_samples)
		memcpy(mu, mu_gamma_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev2::getMuTheta0Samples(int *c, int *l, double* mu)
{
	int C = (*c) - 1;

	if (mu_theta_0_samples)
		memcpy(mu, mu_theta_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev2::getTau2Gamma0Samples(int *c, int *l, double* tau2)
{
	int C = (*c) - 1;

	if (tau2_gamma_0_samples)
		memcpy(tau2, tau2_gamma_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev2::getTau2Theta0Samples(int *c, int *l, double* tau2)
{
	int C = (*c) - 1;

	if (tau2_theta_0_samples)
		memcpy(tau2, tau2_theta_0_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev2::getAlphaPiSamples(int *c, int *l, double* alpha_pi)
{
	int C = (*c) - 1;

	if (alpha_pi_samples)
		memcpy(alpha_pi, alpha_pi_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev2::getBetaPiSamples(int *c, int *l, double* beta_pi)
{
	int C = (*c) - 1;

	if (beta_pi_samples)
		memcpy(beta_pi, beta_pi_samples[C], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev2::getAlphaPiAccept(int *c, int *l, double* acc)
{
	int C = (*c) - 1;

	*acc = alpha_pi_acc[C];
}

void c212BB_poisson_mc_hier3_lev2::getBetaPiAccept(int *c, int* l,  double* acc)
{
	int C = (*c) - 1;

	*acc = beta_pi_acc[C];
}
