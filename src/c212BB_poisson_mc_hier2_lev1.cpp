#include<cstdio>
#include<cstdlib>

#include<cstring>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <Rinternals.h>


#include "c2121a_poisson_mc_hier2_lev0.h"
#include "c212BB_poisson_mc_hier2_lev0.h"
#include "c212BB_poisson_mc_hier2_lev1.h"

static const char *rcsId = "$Id: c212BB_poisson_mc_hier2_lev1.cpp,v 1.6 2016/12/02 13:48:53 clb13102 Exp clb13102 $";

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

void c212BB_poisson_mc_hier2_lev1::init(SEXP sChains, SEXP sBurnin, SEXP sIter,
		SEXP sSim_Type, SEXP sMem_Model, SEXP sGlobal_Sim_Params, 
		SEXP sSim_Params,
		SEXP MH_weight,
		SEXP pm_weights,
		SEXP sMonitor,
		SEXP sNumIntervals,
		SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
		SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
		SEXP palpha_gamma,
		SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta, SEXP pmu_gamma_0,
		SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0, SEXP pmu_gamma,
		SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
		SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi,
		SEXP palgo, SEXP padapt_phase)

{
	release();
	c212BB_poisson_mc_hier2_lev0::release();
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

	mu_gamma_0 = *(REAL(pmu_gamma_0));
	tau2_gamma_0 = *(REAL(ptau2_gamma_0));
	mu_theta_0 = *(REAL(pmu_theta_0));
	tau2_theta_0 = *(REAL(ptau2_theta_0));
	alpha_gamma = *(REAL(palpha_gamma));
	beta_gamma = *(REAL(pbeta_gamma));
	alpha_theta = *(REAL(palpha_theta));
	beta_theta = *(REAL(pbeta_theta));

	alpha_pi = *(REAL(palpha_pi));
	beta_pi = *(REAL(pbeta_pi));

	//Rprintf("Setting alpha_pi = %f, beta_pi = %f\n", alpha_pi, beta_pi);

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

	// The samples
	l = 0;
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

	// Global Simulation Parameters
	initGlobalSimParams(sSim_Type, sGlobal_Sim_Params);

	// Individual Simulation Parameters
	initSimParams(sSim_Params);

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
	//			Rprintf("\t%d %d %d: %0.6f\n", c, l , b, sigma2_gamma[c][b]);
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
	int c = 0, l = 0, b = 0;

	if (gPi != NULL) {
		for (c = 0; c < gChains; c++) {
			free(gPi[c]);
		}
		free(gPi);
		gPi = 0;
	}

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

	l = 0;
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

SEXP c212BB_poisson_mc_hier2_lev1::getL2Samples(double*** &data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = allocVector(REALSXP, gChains * gMaxBs * (gIter - gBurnin)));

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

	PROTECT(dim = allocVector(INTSXP, 3));

	INTEGER(dim)[0] = (gIter - gBurnin);
	INTEGER(dim)[1] = gMaxBs;
	INTEGER(dim)[2] = gChains;

	setAttrib(samples, R_DimSymbol, dim);

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
