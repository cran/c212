#include<cstdio>
#include<cstdlib>
#include <cstring>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "c2121a_poisson_mc_hier2_lev0.h"
#include "c2121a_poisson_mc_hier3_lev0.h"

static const char *rcsId = "$Id: c2121a_poisson_mc_hier3_lev0.cpp,v 1.15 2016/08/25 15:13:14 clb13102 Exp clb13102 $";

const char* c2121a_poisson_mc_hier3_lev0::sMonitor_mu_theta_0 = "mu.theta.0";
const char* c2121a_poisson_mc_hier3_lev0::sMonitor_mu_gamma_0 = "mu.gamma.0";
const char* c2121a_poisson_mc_hier3_lev0::sMonitor_tau2_theta_0 = "tau2.theta.0";
const char* c2121a_poisson_mc_hier3_lev0::sMonitor_tau2_gamma_0 = "tau2.gamma.0";

c2121a_poisson_mc_hier3_lev0::c2121a_poisson_mc_hier3_lev0()
{
	//Rprintf("c2121a_poisson_mc_hier3_lev0::c2121a_poisson_mc_hier3_lev0: Default constructor\n");

	iMonitor_mu_theta_0 = 0;
	iMonitor_mu_gamma_0 = 0;
	iMonitor_tau2_theta_0 = 0;
	iMonitor_tau2_gamma_0 = 0;

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

	mu_theta_0_samples = NULL;
	mu_gamma_0_samples = NULL;
	tau2_theta_0_samples = NULL;
	tau2_gamma_0_samples = NULL;
}

c2121a_poisson_mc_hier3_lev0::c2121a_poisson_mc_hier3_lev0(SEXP sChains, SEXP sBurnin, SEXP sIter,
					SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_cntrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs,
					SEXP sNAE, SEXP pX,
					SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP pmu_gamma_0_0,
					SEXP ptau2_gamma_0_0, SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0,
					SEXP palpha_gamma_0_0,
					SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0, SEXP pbeta_theta_0_0,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta)
{
	iMonitor_mu_theta_0 = 0;
	iMonitor_mu_gamma_0 = 0;
  	iMonitor_tau2_theta_0 = 0;
  	iMonitor_tau2_gamma_0 = 0;

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

	mu_theta_0_samples = NULL;
	mu_gamma_0_samples = NULL;
	tau2_theta_0_samples = NULL;
	tau2_gamma_0_samples = NULL;

	init(sChains, sBurnin, sIter, sSim_Type, sMem_Model, sGlobal_Sim_Param,
				sGlobal_Sim_Param_cntrl, sSim_Param,
				sMonitor,
				sNumIntervals, sMaxBs, sNumBodySys, sMaxAEs, sNAE,
				pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0,
				pmu_theta_0_0, ptau2_theta_0_0,
				palpha_gamma_0_0, pbeta_gamma_0_0, palpha_theta_0_0, pbeta_theta_0_0,
				palpha_gamma,
				pbeta_gamma, palpha_theta, pbeta_theta, pmu_gamma_0, ptau2_gamma_0,
				pmu_theta_0,
				ptau2_theta_0, pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta);

}

void c2121a_poisson_mc_hier3_lev0::init(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_cntrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals,
					SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP pmu_gamma_0_0,
					SEXP ptau2_gamma_0_0, SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0,
					SEXP palpha_gamma_0_0,
					SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0, SEXP pbeta_theta_0_0,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta)

{
	release();
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

	mu_gamma_0 = (double**)malloc(gChains * sizeof(double*));
	double *vmu_gamma_0 = REAL(pmu_gamma_0);
	for (c = 0; c < gChains; c++) {
		mu_gamma_0[c] = (double *)malloc(gNumIntervals * sizeof(double));
		for (l = 0; l < gNumIntervals; l++) {
			mu_gamma_0[c][l] = *vmu_gamma_0;
			vmu_gamma_0++;
		}
	}

	mu_theta_0 = (double**)malloc(gChains * sizeof(double*));
	double *vmu_theta_0 = REAL(pmu_theta_0);
	for (c = 0; c < gChains; c++) {
		mu_theta_0[c] = (double *)malloc(gNumIntervals * sizeof(double));
		for (l = 0; l < gNumIntervals; l++) {
			mu_theta_0[c][l] = *vmu_theta_0;
			vmu_theta_0++;
		}
	}

	tau2_gamma_0 = (double**)malloc(gChains * sizeof(double*));
	double *vtau2_gamma_0 = REAL(ptau2_gamma_0);
	for (c = 0; c < gChains; c++) {
		tau2_gamma_0[c] = (double *)malloc(gNumIntervals * sizeof(double));
		for (l = 0; l < gNumIntervals; l++) {
			tau2_gamma_0[c][l] = *vtau2_gamma_0;
			vtau2_gamma_0++;
		}
	}

	tau2_theta_0 = (double**)malloc(gChains * sizeof(double*));
	double *vtau2_theta_0 = REAL(ptau2_theta_0);
	for (c = 0; c < gChains; c++) {
		tau2_theta_0[c] = (double *)malloc(gNumIntervals * sizeof(double));
		for (l = 0; l < gNumIntervals; l++) {
			tau2_theta_0[c][l] = *vtau2_theta_0;
			vtau2_theta_0++;
		}
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

	// The samples
	if (retainSamples(iMonitor_mu_gamma_0))
		mu_gamma_0_samples = (double ***)malloc(gChains * sizeof(double**));
	if (retainSamples(iMonitor_mu_theta_0))
		mu_theta_0_samples = (double ***)malloc(gChains * sizeof(double**));
	if (retainSamples(iMonitor_tau2_gamma_0))
		tau2_gamma_0_samples = (double ***)malloc(gChains * sizeof(double**));
	if (retainSamples(iMonitor_tau2_theta_0))
		tau2_theta_0_samples = (double ***)malloc(gChains * sizeof(double**));

	for (c = 0; c < gChains; c++) {
		if (retainSamples(iMonitor_mu_gamma_0))
			mu_gamma_0_samples[c] = (double **)malloc(gNumIntervals* sizeof(double *));
		if (retainSamples(iMonitor_mu_theta_0))
			mu_theta_0_samples[c] = (double **)malloc(gNumIntervals* sizeof(double *));
		if (retainSamples(iMonitor_tau2_gamma_0))
			tau2_gamma_0_samples[c] =
									(double **)malloc(gNumIntervals* sizeof(double *));
		if (retainSamples(iMonitor_tau2_theta_0))
			tau2_theta_0_samples[c] =
									(double **)malloc(gNumIntervals* sizeof(double *));
		for (l = 0; l < gNumIntervals; l++) {
			if (retainSamples(iMonitor_mu_gamma_0))
				mu_gamma_0_samples[c][l] =
									(double *)malloc((gIter - gBurnin)* sizeof(double));
			if (retainSamples(iMonitor_mu_theta_0))
				mu_theta_0_samples[c][l] =
									(double *)malloc((gIter - gBurnin)* sizeof(double));
			if (retainSamples(iMonitor_tau2_gamma_0))
				tau2_gamma_0_samples[c][l] =
									(double *)malloc((gIter - gBurnin)* sizeof(double));
			if (retainSamples(iMonitor_tau2_theta_0))
				tau2_theta_0_samples[c][l] =
									(double *)malloc((gIter - gBurnin)* sizeof(double));
		}
	}

	if (retainSamples(iMonitor_mu_theta))
		mu_theta_samples = (double ****)malloc(gChains *sizeof(double***));
	if (retainSamples(iMonitor_mu_gamma))
		mu_gamma_samples = (double ****)malloc(gChains *sizeof(double***));
	if (retainSamples(iMonitor_sigma2_theta))
		sigma2_theta_samples = (double ****)malloc(gChains *sizeof(double***));
	if (retainSamples(iMonitor_sigma2_gamma))
		sigma2_gamma_samples = (double ****)malloc(gChains *sizeof(double***));

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
		for (l = 0; l < gNumIntervals; l++) {
			if (retainSamples(iMonitor_mu_theta))
				mu_theta_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));
			if (retainSamples(iMonitor_mu_gamma))
				mu_gamma_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));
			if (retainSamples(iMonitor_sigma2_theta))
				sigma2_theta_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));
			if (retainSamples(iMonitor_sigma2_gamma))
				sigma2_gamma_samples[c][l] = (double **)malloc(gMaxBs *sizeof(double*));

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

	// Simulation Parameters
	l = strlen(CHAR(STRING_ELT(sSim_Type, 0)));
	sim_type = (char *)malloc((l + 1)*sizeof(char));
	if (sim_type) {
		strcpy(sim_type, CHAR(STRING_ELT(sSim_Type, 0)));
		sim_type[l] = 0;
	}

	// Global Simulation Parameters
	gSim_Param = *REAL(sGlobal_Sim_Param);
	gSim_Param_cntrl = *REAL(sGlobal_Sim_Param_cntrl);

	// Individual Simulation Parameters
	initSimParams(sSim_Param);

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

void c2121a_poisson_mc_hier3_lev0::initMonitor(SEXP sMonitor)
{
    int len = Rf_length(sMonitor);

    SEXP sVariables = R_NilValue;
    SEXP sValues = R_NilValue;

    if (len > 0 && isNewList(sMonitor)) {

        SEXP names = getAttrib(sMonitor, R_NamesSymbol);

        int i = 0;

        for (i = 0; i < len; i++) {
            if (strcmp(sColMonitorVariables, CHAR(STRING_ELT(names, i))) == 0) {
                sVariables = VECTOR_ELT(sMonitor, i);
            }
            if (strcmp(sColMonitorValues, CHAR(STRING_ELT(names, i))) == 0) {
                sValues = VECTOR_ELT(sMonitor, i);
            }
        }
            
        len = Rf_length(sVariables);

		if (len > 0) {

	        int* vals = INTEGER(sValues);
	
	        for (i = 0; i < len; i++) {
				const char *t = CHAR(STRING_ELT(sVariables, i));
	
	            if (0 == strcmp(t, sMonitor_theta)) {
	                iMonitor_theta = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_gamma)) {
	                iMonitor_gamma = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_mu_theta)) {
	                iMonitor_mu_theta = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_mu_gamma)) {
	                iMonitor_mu_gamma = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_sigma2_theta)) {
	                iMonitor_sigma2_theta = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_sigma2_gamma)) {
	                iMonitor_sigma2_gamma = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_mu_theta_0)) {
	                iMonitor_mu_theta_0 = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_mu_gamma_0)) {
	                iMonitor_mu_gamma_0 = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_tau2_gamma_0)) {
	                iMonitor_tau2_gamma_0 = vals[i];
	            }
	            if (0 == strcmp(t, sMonitor_tau2_theta_0)) {
	                iMonitor_tau2_theta_0 = vals[i];
	            }
			}
		}
	}
}

c2121a_poisson_mc_hier3_lev0::~c2121a_poisson_mc_hier3_lev0()
{
	release();
}

void c2121a_poisson_mc_hier3_lev0::gibbs_sampler()
{
	if (strcmp(sim_type, "MH") == 0) {
		simulate_MH();
	}
	else {
		simulate_SLICE();
	}

	return;
}

void c2121a_poisson_mc_hier3_lev0::simulate_MH()
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

void c2121a_poisson_mc_hier3_lev0::simulate_SLICE()
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

void c2121a_poisson_mc_hier3_lev0::sample_mu_gamma_0(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			double denom = tau2_gamma_0[c][l] + tau2_gamma_0_0 * ((double)gNumBodySys[l]);

			double mu_gamma_tot = 0.0;
			int i = 0;

			for (i = 0; i < gNumBodySys[l]; i++) {
				mu_gamma_tot += mu_gamma[c][l][i];
			}

			double mean = (tau2_gamma_0[c][l] * mu_gamma_0_0 + tau2_gamma_0_0 * mu_gamma_tot)/denom;
			double var = (tau2_gamma_0[c][l] * tau2_gamma_0_0) / denom;

			double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			mu_gamma_0[c][l] = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			if (iter >= burnin && retainSamples(iMonitor_mu_gamma_0)) {
				mu_gamma_0_samples[c][l][iter - burnin] = mu_gamma_0[c][l];
			}
		}
	}
}

void c2121a_poisson_mc_hier3_lev0::sample_mu_theta_0 (int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			double denom = tau2_theta_0[c][l] + tau2_theta_0_0 * ((double)gNumBodySys[l]);

			double mu_theta_tot = 0.0;
			int i = 0;

			for (i = 0; i < gNumBodySys[l]; i++) {
				mu_theta_tot += mu_theta[c][l][i];
			}

			double mean = (tau2_theta_0[c][l] * mu_theta_0_0 + tau2_theta_0_0 * mu_theta_tot)/denom;
			double var = (tau2_theta_0[c][l] * tau2_theta_0_0) / denom;

			double sd = sqrt(var);

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			mu_theta_0[c][l] = rnorm(mean, sd);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			if (iter >= burnin && retainSamples(iMonitor_mu_theta_0)) {
				mu_theta_0_samples[c][l][iter - burnin] = mu_theta_0[c][l];
			}
		}
	}
}

void c2121a_poisson_mc_hier3_lev0::sample_tau2_gamma_0(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			double s = alpha_gamma_0_0 + ((double)gNumBodySys[l])/2.0;
			double r = 0.0;
			double isum = 0.0;

			int i = 0;
			for (i = 0; i < gNumBodySys[l]; i++) {
				isum += (pow((mu_gamma[c][l][i] - mu_gamma_0[c][l]), 2.0));
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
	
			tau2_gamma_0[c][l] = 1/cand;

			if (iter >= burnin && retainSamples(iMonitor_tau2_gamma_0)) {
				tau2_gamma_0_samples[c][l][iter - burnin] = tau2_gamma_0[c][l];
			}
		}
	}
}

void c2121a_poisson_mc_hier3_lev0::sample_tau2_theta_0(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			double s = alpha_theta_0_0 + ((double)gNumBodySys[l])/2.0;
			double isum = 0.0;
			double r = 0.0;

			int i = 0;
			for (i = 0; i < gNumBodySys[l]; i++) {
				isum += (pow((mu_theta[c][l][i] - mu_theta_0[c][l]), 2.0));
			}

			r = beta_theta_0_0 + 0.5 * isum;

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double cand = rgamma(s, 1/r);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			tau2_theta_0[c][l] = 1/cand;

			if (iter >= burnin && retainSamples(iMonitor_tau2_theta_0)) {
				tau2_theta_0_samples[c][l][iter - burnin] = tau2_theta_0[c][l];
			}
		}
	}
}

void c2121a_poisson_mc_hier3_lev0::sample_mu_gamma(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {
				double denom = sigma2_gamma[c][l][b] + ((double)gNAE[l][b])*tau2_gamma_0[c][l];


				double t = 0.0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					t += gGamma[c][l][b][j];
				}

				double mean = (sigma2_gamma[c][l][b] * mu_gamma_0[c][l] + tau2_gamma_0[c][l] * t)/denom;

				double var = (sigma2_gamma[c][l][b]*tau2_gamma_0[c][l])/denom;

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

void c2121a_poisson_mc_hier3_lev0::sample_mu_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {
				double denom = sigma2_theta[c][l][b] + ((double)gNAE[l][b])*tau2_theta_0[c][l];

				double t = 0.0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					t += gTheta[c][l][b][j];
				}

				double mean = (sigma2_theta[c][l][b] * mu_theta_0[c][l] + tau2_theta_0[c][l] * t)/denom;

				double var = (sigma2_theta[c][l][b]*tau2_theta_0[c][l])/denom;

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

void c2121a_poisson_mc_hier3_lev0::sample_sigma2_gamma(int burnin, int iter)
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

void c2121a_poisson_mc_hier3_lev0::sample_sigma2_theta(int burnin, int iter)
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

double c2121a_poisson_mc_hier3_lev0::cMIN(double a, double b)
{
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

void c2121a_poisson_mc_hier3_lev0::release()
{
	int c = 0, l = 0;

	if (mu_theta_0 != NULL) {
		for (c = 0; c < gChains; c++) {
			free(mu_theta_0[c]);
		}
		free(mu_theta_0);
		mu_theta_0 = NULL;
	}

	if (mu_gamma_0 != NULL) {
		for (c = 0; c < gChains; c++) {
			free(mu_gamma_0[c]);
		}
		free(mu_gamma_0);
		mu_gamma_0 = 0;
	}

	if (tau2_theta_0 != NULL) {
		for (c = 0; c < gChains; c++) {
			free(tau2_theta_0[c]);
		}
		free(tau2_theta_0);
		tau2_theta_0 = NULL;
	}

	if (tau2_gamma_0 != NULL) {
		for (c = 0; c < gChains; c++) {
			free(tau2_gamma_0[c]);
		}
		free(tau2_gamma_0);
		tau2_gamma_0 = NULL;
	}

	if (mu_gamma_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(mu_gamma_0_samples[c][l]);
			}
			free(mu_gamma_0_samples[c]);
		}
		free(mu_gamma_0_samples);
		mu_gamma_0_samples = NULL;
	}
	if (mu_theta_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(mu_theta_0_samples[c][l]);
			}
			free(mu_theta_0_samples[c]);
		}
		free(mu_theta_0_samples);
		mu_theta_0_samples = NULL;
	}
	if (tau2_gamma_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(tau2_gamma_0_samples[c][l]);
			}
			free(tau2_gamma_0_samples[c]);
		}
		free(tau2_gamma_0_samples);
		tau2_gamma_0_samples = NULL;
	}
	if (tau2_theta_0_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(tau2_theta_0_samples[c][l]);
			}
			free(tau2_theta_0_samples[c]);
		}
		free(tau2_theta_0_samples);
		tau2_theta_0_samples = NULL;
	}
}

SEXP c2121a_poisson_mc_hier3_lev0::getL3Samples(double*** &data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = allocVector(REALSXP, gChains * gNumIntervals * (gIter - gBurnin)));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		int l = 0;
		for (l = 0; l < gNumIntervals; l++) {
			memcpy(REAL(samples) + i, data[c][l], (gIter - gBurnin)*sizeof(double));

			i += (gIter - gBurnin);
			free(data[c][l]);
			data[c][l] = NULL;
		}
		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = allocVector(INTSXP, 3));

	INTEGER(dim)[0] = (gIter - gBurnin);
	INTEGER(dim)[1] = gNumIntervals;
	INTEGER(dim)[2] = gChains;

	setAttrib(samples, R_DimSymbol, dim);

	UNPROTECT(2);

	return samples;
}

SEXP c2121a_poisson_mc_hier3_lev0::getMuGamma0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(mu_gamma_0_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier3_lev0::getMuTheta0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(mu_theta_0_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier3_lev0::getTau2Gamma0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(tau2_gamma_0_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier3_lev0::getTau2Theta0Samples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(tau2_theta_0_samples);

	return samples;
}

void c2121a_poisson_mc_hier3_lev0::getMuGamma0Samples(int *c, int *l, double* mu)
{
	int C = (*c) - 1;
	int L = (*l) - 1;

	if (mu_gamma_0_samples)
		memcpy(mu, mu_gamma_0_samples[C][L], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier3_lev0::getMuTheta0Samples(int *c, int *l, double* mu)
{
	int C = (*c) - 1;
	int L = (*l) - 1;

	if (mu_theta_0_samples)
		memcpy(mu, mu_theta_0_samples[C][L], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier3_lev0::getTau2Gamma0Samples(int *c, int *l, double* tau2)
{
	int C = (*c) - 1;
	int L = (*l) - 1;

	if (tau2_gamma_0_samples)
		memcpy(tau2, tau2_gamma_0_samples[C][L], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier3_lev0::getTau2Theta0Samples(int *c, int *l, double* tau2)
{
	int C = (*c) - 1;
	int L = (*l) - 1;

	if (tau2_theta_0_samples)
		memcpy(tau2, tau2_theta_0_samples[C][L], (gIter - gBurnin)*sizeof(double));
}
