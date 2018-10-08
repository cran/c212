#include<cstdio>
#include<cstdlib>
#include <cstring>
#include<cmath>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "c2121a_poisson_mc_hier2_lev0.h"

using namespace std;

static const char *rcsId = "$Id: c2121a_poisson_mc_hier2_lev0.cpp,v 1.7 2018/10/03 15:40:28 clb13102 Exp clb13102 $";

const char* c2121a_poisson_mc_hier2_lev0::sColType = "type";
const char* c2121a_poisson_mc_hier2_lev0::sColVariable = "variable";
const char* c2121a_poisson_mc_hier2_lev0::sColParam = "param";
const char* c2121a_poisson_mc_hier2_lev0::sColValue = "value";
const char* c2121a_poisson_mc_hier2_lev0::sColControl = "control";
const char* c2121a_poisson_mc_hier2_lev0::sColB = "B";
const char* c2121a_poisson_mc_hier2_lev0::sColj = "j";
const char* c2121a_poisson_mc_hier2_lev0::sColI_index = "I_index";

const char* c2121a_poisson_mc_hier2_lev0::sParam_w = "w";
const char* c2121a_poisson_mc_hier2_lev0::sParam_sigma_MH = "sigma_MH";
const char* c2121a_poisson_mc_hier2_lev0::sVariable_gamma = "gamma";
const char* c2121a_poisson_mc_hier2_lev0::sVariable_theta = "theta";

const char* c2121a_poisson_mc_hier2_lev0::sColMonitorVariables = "variable";
const char* c2121a_poisson_mc_hier2_lev0::sColMonitorValues = "monitor";

const char* c2121a_poisson_mc_hier2_lev0::sMonitor_theta = "theta";
const char* c2121a_poisson_mc_hier2_lev0::sMonitor_gamma = "gamma";
const char* c2121a_poisson_mc_hier2_lev0::sMonitor_mu_theta = "mu.theta";
const char* c2121a_poisson_mc_hier2_lev0::sMonitor_mu_gamma = "mu.gamma";
const char* c2121a_poisson_mc_hier2_lev0::sMonitor_sigma2_theta = "sigma2.theta";
const char* c2121a_poisson_mc_hier2_lev0::sMonitor_sigma2_gamma = "sigma2.gamma";

c2121a_poisson_mc_hier2_lev0::c2121a_poisson_mc_hier2_lev0()
{
	//Rprintf("c2121a_poisson_mc_hier2_lev0::c2121a_poisson_mc_hier2_lev0: Default constructor\n");
	gChains = 0;
	gBurnin = 0;
	gIter = 0;
	gMaxBs = 0;
	sim_type = NULL;
	gNAE = NULL;
	gNumIntervals = 0;
	gNumBodySys = NULL;
	gMaxAEs = 0;
	gSim_Param = 0.0;
	gSim_Param_cntrl = 0.0;
	sim_type = NULL;

 	gW_gamma = NULL;
	gW_theta = NULL;
	gW_gamma_control = NULL;
	gW_theta_control = NULL;
	gSigma_MH_gamma = NULL;
	gSigma_MH_theta = NULL;

	eMemory_Model = HIGH;
	iMonitor_theta = 0;
	iMonitor_gamma = 0;
	iMonitor_mu_theta = 0;
	iMonitor_mu_gamma = 0;
	iMonitor_sigma2_theta = 0;
	iMonitor_sigma2_gamma = 0;

	mu_theta_0 = 0.0;
	mu_gamma_0 = 0.0;
	tau2_theta_0 = 0.0;
	tau2_gamma_0 = 0.0;
	alpha_gamma = 0.0;
	beta_gamma = 0.0;
	alpha_theta = 0.0;
	beta_theta = 0.0;

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
	C = NULL;
	T = NULL;

	gTheta_samples = NULL;
	gGamma_samples = NULL;
	mu_theta_samples = NULL;
	mu_gamma_samples = NULL;
	sigma2_theta_samples = NULL;
	sigma2_gamma_samples = NULL;
}

c2121a_poisson_mc_hier2_lev0::c2121a_poisson_mc_hier2_lev0(SEXP sChains, SEXP sBurnin,
					SEXP sIter,
					SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_cntrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs,
					SEXP sNAE, SEXP pX,
					SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta)
{
	gChains = 0;
	gBurnin = 0;
	gIter = 0;
	sim_type = NULL;
	gNAE = NULL;
	gNumIntervals = 0;
	gMaxBs = 0;
	gNumBodySys = NULL;
	gMaxAEs = 0;
	gSim_Param = 0.0;
	gSim_Param_cntrl = 0.0;
	sim_type = NULL;

 	gW_gamma = NULL;
	gW_theta = NULL;
	gW_gamma_control = NULL;
	gW_theta_control = NULL;
	gSigma_MH_gamma = NULL;
	gSigma_MH_theta = NULL;

	eMemory_Model = HIGH;
	iMonitor_theta = 0;
	iMonitor_gamma = 0;
	iMonitor_mu_theta = 0;
	iMonitor_mu_gamma = 0;
	iMonitor_sigma2_theta = 0;
	iMonitor_sigma2_gamma = 0;

	mu_theta_0 = 0.0;
	mu_gamma_0 = 0.0;
	tau2_theta_0 = 0.0;
	tau2_gamma_0 = 0.0;
	alpha_gamma = 0.0;
	beta_gamma = 0.0;
	alpha_theta = 0.0;
	beta_theta = 0.0;

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
	C = NULL;
	T = NULL;

	gTheta_samples = NULL;
	gGamma_samples = NULL;
	mu_theta_samples = NULL;
	mu_gamma_samples = NULL;
	sigma2_theta_samples = NULL;
	sigma2_gamma_samples = NULL;

	init(sChains, sBurnin, sIter, sSim_Type, sMem_Model, sGlobal_Sim_Param,
				sGlobal_Sim_Param_cntrl,
				sSim_Param,
				sMonitor,
				sNumIntervals, sMaxBs, sNumBodySys, sMaxAEs, sNAE,
				pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0, ptau2_gamma_0,
				pmu_theta_0, ptau2_theta_0,
				palpha_gamma, pbeta_gamma, palpha_theta, pbeta_theta,
				pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta);

}

void c2121a_poisson_mc_hier2_lev0::init(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_cntrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals,
					SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta)
{
	clear();

	initMonitor(sMonitor);

	initBaselineVariables(sChains, sBurnin, sIter,
				sMem_Model, sNumIntervals, sMaxBs, sNumBodySys, sMaxAEs, sNAE);

	initDataVariables(pX, pY, pC, pT);

	initL1Variables(ptheta, pgamma);

	initL2Params(pmu_gamma_0, ptau2_gamma_0, pmu_theta_0, ptau2_theta_0,
						palpha_gamma, pbeta_gamma, palpha_theta, pbeta_theta);

	initL2Variables(pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta);

	initL2Samples();

	initL1Samples();

	initGlobalSimParams(sSim_Type, sGlobal_Sim_Param, sGlobal_Sim_Param_cntrl);

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
	//Rprintf("mu.gamma.0: %0.6f\n", mu_gamma_0);
	//Rprintf("tau2.gamma.0: %0.6f\n", tau2_gamma_0);
	//Rprintf("mu.theta.0: %0.6f\n", mu_theta_0);
	//Rprintf("tau2.theta.0: %0.6f\n", tau2_theta_0);
	//Rprintf("alpha.gamma: %0.6f\n", alpha_gamma);
	//Rprintf("beta.gamma: %0.6f\n", beta_gamma);
	//Rprintf("alpha.theta: %0.6f\n", alpha_theta);
	//Rprintf("beta.theta: %0.6f\n", beta_theta);
	//Rprintf("alpha.gamma: %0.6f\n", alpha_gamma);
	//Rprintf("beta.gamma: %0.6f\n", beta_gamma);
	//Rprintf("alpha.theta: %0.6f\n", alpha_theta);
	//Rprintf("beta.theta: %0.6f\n", beta_theta);
}


void c2121a_poisson_mc_hier2_lev0::initBaselineVariables(SEXP sChains, SEXP sBurnin, SEXP sIter,
                    SEXP sMem_Model, SEXP sNumIntervals, SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE)
{
	int l = 0, b = 0;

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
}

void c2121a_poisson_mc_hier2_lev0::clear()
{
	release();
}

void c2121a_poisson_mc_hier2_lev0::releaseBaselineVariables()
{
	int l = 0;

	if (gNumBodySys != NULL) {
		free(gNumBodySys);
		gNumBodySys = NULL;
	}

	if (gNAE != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			free(gNAE[l]);
		}
		free(gNAE);
		gNAE = NULL;
	}
}

void c2121a_poisson_mc_hier2_lev0::initGlobalSimParams(SEXP sSim_Type, SEXP sGlobal_Sim_Param, SEXP sGlobal_Sim_Param_cntrl)
{
	int l = 0;

	// Simulations Parameters
	l = strlen(CHAR(STRING_ELT(sSim_Type, 0)));
	sim_type = (char *)malloc((l + 1)*sizeof(char));
	if (sim_type) {
		strcpy(sim_type, CHAR(STRING_ELT(sSim_Type, 0)));
		sim_type[l] = 0;
	}

	// Global Simulation Parameters
	gSim_Param = *REAL(sGlobal_Sim_Param);
	gSim_Param_cntrl = *REAL(sGlobal_Sim_Param_cntrl);
}

void c2121a_poisson_mc_hier2_lev0::releaseGlobalSimParams()
{
	if (sim_type != NULL) {
		free(sim_type);
		sim_type = NULL;
	}
}

void c2121a_poisson_mc_hier2_lev0::initSimParams(SEXP sSim_Params)
{
	gW_gamma = (double ***)malloc(gNumIntervals * sizeof(double**));
	gW_theta = (double***)malloc(gNumIntervals * sizeof(double**));
	gW_gamma_control = (int ***)malloc(gNumIntervals * sizeof(int**));
	gW_theta_control = (int ***)malloc(gNumIntervals * sizeof(int**));
	gSigma_MH_gamma = (double ***)malloc(gNumIntervals * sizeof(double**));
	gSigma_MH_theta = (double ***)malloc(gNumIntervals * sizeof(double**));

	int i = 0, b = 0, j = 0;

	for (i = 0; i < gNumIntervals; i++) {

		gW_gamma[i] = (double**)malloc(gNumBodySys[i] * sizeof(double*));
		gW_theta[i] = (double**)malloc(gNumBodySys[i] * sizeof(double*));
		gW_gamma_control[i] = (int**)malloc(gNumBodySys[i] * sizeof(int*));
		gW_theta_control[i] = (int**)malloc(gNumBodySys[i] * sizeof(int*));
		gSigma_MH_gamma[i] = (double**)malloc(gNumBodySys[i] * sizeof(double*));
		gSigma_MH_theta[i] = (double**)malloc(gNumBodySys[i] * sizeof(double*));

		for (b = 0; b < gNumBodySys[i]; b++) {

			gW_gamma[i][b] = (double*)malloc(gNAE[i][b] * sizeof(double));
			gW_theta[i][b] = (double*)malloc(gNAE[i][b] * sizeof(double));
			gW_gamma_control[i][b] = (int*)malloc(gNAE[i][b] * sizeof(int));
			gW_theta_control[i][b] = (int*)malloc(gNAE[i][b] * sizeof(int));
			gSigma_MH_gamma[i][b] = (double*)malloc(gNAE[i][b] * sizeof(double));
			gSigma_MH_theta[i][b] = (double*)malloc(gNAE[i][b] * sizeof(double));

			for (j = 0; j < gNAE[i][b]; j++) {
				gW_gamma[i][b][j] = gSim_Param;
				gW_theta[i][b][j] = gSim_Param;
				gW_gamma_control[i][b][j] = (int)gSim_Param_cntrl;
				gW_theta_control[i][b][j] = (int)gSim_Param_cntrl;
				gSigma_MH_gamma[i][b][j] = gSim_Param;
				gSigma_MH_theta[i][b][j] = gSim_Param;
			}
		}
	}

	int len = Rf_length(sSim_Params);

	if (len && isNewList(sSim_Params)) {

		SEXP sVariables = R_NilValue;
		SEXP sParams = R_NilValue;
		SEXP sValues = R_NilValue;
		SEXP sControl = R_NilValue;
		SEXP sB = R_NilValue;
		SEXP sj = R_NilValue;
		//SEXP sIntervals = R_NilValue;
		SEXP sI_index = R_NilValue;

		SEXP names = getAttrib(sSim_Params, R_NamesSymbol);

		for (i = 0; i < len; i++) {
			if (strcmp(sColValue, CHAR(STRING_ELT(names, i))) == 0) {
				sValues = VECTOR_ELT(sSim_Params, i);
			}
			if (strcmp(sColParam, CHAR(STRING_ELT(names, i))) == 0) {
				sParams = VECTOR_ELT(sSim_Params, i);
			}
			if (strcmp(sColControl, CHAR(STRING_ELT(names, i))) == 0) {
				sControl = VECTOR_ELT(sSim_Params, i);
			}
			if (strcmp(sColVariable, CHAR(STRING_ELT(names, i))) == 0) {
				sVariables = VECTOR_ELT(sSim_Params, i);
			}
			if (strcmp(sColB, CHAR(STRING_ELT(names, i))) == 0) {
				sB = VECTOR_ELT(sSim_Params, i);
			}
			if (strcmp(sColj, CHAR(STRING_ELT(names, i))) == 0) {
				sj = VECTOR_ELT(sSim_Params, i);
			}
			if (strcmp(sColI_index, CHAR(STRING_ELT(names, i))) == 0) {
				sI_index = VECTOR_ELT(sSim_Params, i);
			}
		}

		len = Rf_length(sParams);

		if (len > 0) {
			double* vals = REAL(sValues);
			double* cntrl = REAL(sControl);
			int* B = INTEGER(sB);
			int* j = INTEGER(sj);
			int* i_index = INTEGER(sI_index);
	
			for (i = 0; i < len; i++) {
				const char *var = CHAR(STRING_ELT(sVariables, i));
				const char *param = CHAR(STRING_ELT(sParams, i));
	
				int l = i_index[i] - 1;
				int b = B[i] - 1;
				int a = j[i] - 1;
				if (0 == strcmp(sVariable_gamma, var)) {
					if (0 == strcmp(param, sParam_w)) {
						gW_gamma[l][b][a] = vals[i];
						gW_gamma_control[l][b][a] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH)) {
						gSigma_MH_gamma[l][b][a] = vals[i];
					}
				}
				else if (0 == strcmp(sVariable_theta, var)) {
					if (0 == strcmp(param, sParam_w)) {
						gW_theta[l][b][a] = vals[i];
						gW_theta_control[l][b][a] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH)) {
						gSigma_MH_theta[l][b][a] = vals[i];
					}
				}
			}
		}
	}
}

void c2121a_poisson_mc_hier2_lev0::releaseSimParams()
{
	int l = 0, b = 0;

	if (gW_gamma != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(gW_gamma[l][b]);
			}
			free(gW_gamma[l]);
		}
		free(gW_gamma);
		gW_gamma = NULL;
	}

	if (gW_theta != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(gW_theta[l][b]);
			}
			free(gW_theta[l]);
		}
		free(gW_theta);
		gW_theta = NULL;
	}

	if (gW_gamma_control != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(gW_gamma_control[l][b]);
			}
			free(gW_gamma_control[l]);
		}
		free(gW_gamma_control);
		gW_gamma_control = NULL;
	}

	if (gW_theta_control != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(gW_theta_control[l][b]);
			}
			free(gW_theta_control[l]);
		}
		free(gW_theta_control);
		gW_theta_control = NULL;
	}

	if (gSigma_MH_gamma != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(gSigma_MH_gamma[l][b]);
			}
			free(gSigma_MH_gamma[l]);
		}
		free(gSigma_MH_gamma);
		gSigma_MH_gamma = NULL;
	}

	if (gSigma_MH_theta != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(gSigma_MH_theta[l][b]);
			}
			free(gSigma_MH_theta[l]);
		}
		free(gSigma_MH_theta);
		gSigma_MH_theta = NULL;
	}


}

void c2121a_poisson_mc_hier2_lev0::initMonitor(SEXP sMonitor)
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
			}
		}
	}
}

void c2121a_poisson_mc_hier2_lev0::initL2Params(SEXP pmu_gamma_0,
						SEXP ptau2_gamma_0, SEXP pmu_theta_0,
						SEXP ptau2_theta_0, SEXP palpha_gamma,
						SEXP pbeta_gamma, SEXP palpha_theta,
						SEXP pbeta_theta)
{
	mu_gamma_0 = *(REAL(pmu_gamma_0));
	tau2_gamma_0 = *(REAL(ptau2_gamma_0));
	mu_theta_0 = *(REAL(pmu_theta_0));
	tau2_theta_0 = *(REAL(ptau2_theta_0));
	alpha_gamma = *(REAL(palpha_gamma));
	beta_gamma = *(REAL(pbeta_gamma));
	alpha_theta = *(REAL(palpha_theta));
	beta_theta = *(REAL(pbeta_theta));
}

void c2121a_poisson_mc_hier2_lev0::initL2Variables(SEXP pmu_gamma, SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta)
{
	int c = 0, l = 0, b = 0;

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
}

void c2121a_poisson_mc_hier2_lev0::releaseL2Variables()
{
	int c = 0, l = 0;

	if (mu_gamma != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(mu_gamma[c][l]);
			}
			free(mu_gamma[c]);
		}
		free(mu_gamma);
		mu_gamma = 0;
	}

	if (mu_theta != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(mu_theta[c][l]);
			}
			free(mu_theta[c]);
		}
		free(mu_theta);
		mu_theta = 0;
	}

	if (sigma2_gamma != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(sigma2_gamma[c][l]);
			}
			free(sigma2_gamma[c]);
		}
		free(sigma2_gamma);
		sigma2_gamma = 0;
	}

	if (sigma2_theta != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(sigma2_theta[c][l]);
			}
			free(sigma2_theta[c]);
		}
		free(sigma2_theta);
		sigma2_theta = 0;
	}
}

void c2121a_poisson_mc_hier2_lev0::initL1Variables(SEXP ptheta, SEXP pgamma)
{
	int c = 0, l = 0, b = 0, a = 0;

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
}

void c2121a_poisson_mc_hier2_lev0::releaseL1Variables()
{
	int c = 0, l = 0, b = 0;

	if (gTheta != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gMaxBs; b++) {
					free(gTheta[c][l][b]);
				}
				free(gTheta[c][l]);
			}
			free(gTheta[c]);
		}
		free(gTheta);
		gTheta = NULL;
	}
	if (gGamma != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gMaxBs; b++) {
					free(gGamma[c][l][b]);
				}
				free(gGamma[c][l]);
			}
			free(gGamma[c]);
		}
		free(gGamma);
		gGamma = NULL;
	}
}

void c2121a_poisson_mc_hier2_lev0::initDataVariables(SEXP pX, SEXP pY, SEXP pC, SEXP pT)
{
	int l = 0, b = 0, a = 0;

	x = (int***)malloc(gNumIntervals * sizeof(int**));
	y = (int***)malloc(gNumIntervals * sizeof(int**));
	C = (double***)malloc(gNumIntervals * sizeof(double**));
	T = (double***)malloc(gNumIntervals * sizeof(double**));
	for (l = 0; l < gNumIntervals; l++) {
		x[l] = (int**)malloc(gMaxBs * sizeof(int*));
		y[l] = (int**)malloc(gMaxBs * sizeof(int*));
		C[l] = (double**)malloc(gMaxBs * sizeof(double*));
		T[l] = (double**)malloc(gMaxBs * sizeof(double*));
		for (b = 0; b < gMaxBs; b++) {
			x[l][b] = (int*)malloc(gMaxAEs * sizeof(int));
			y[l][b] = (int*)malloc(gMaxAEs * sizeof(int));
			C[l][b] = (double*)malloc(gMaxAEs * sizeof(double));
			T[l][b] = (double*)malloc(gMaxAEs * sizeof(double));
		}
	}

	int *vX = INTEGER(pX);
	int *vY = INTEGER(pY);
	double *vC = REAL(pC);
	double *vT = REAL(pT);
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
}

void c2121a_poisson_mc_hier2_lev0::releaseDataVariables()
{
	int l = 0, b = 0;

	if (x != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gMaxBs; b++) {
				free(x[l][b]);
			}
			free(x[l]);
		}
		free(x);
		x = NULL;
	}
	if (y != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gMaxBs; b++) {
				free(y[l][b]);
			}
			free(y[l]);
		}
		free(y);
		y = NULL;
	}
	if (C != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gMaxBs; b++) {
				free(C[l][b]);
			}
			free(C[l]);
		}
		free(C);
		C = NULL;
	}
	if (T != NULL) {
		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gMaxBs; b++) {
				free(T[l][b]);
			}
			free(T[l]);
		}
		free(T);
		T = NULL;
	}
}

void c2121a_poisson_mc_hier2_lev0::initL2Samples()
{
	int c = 0, l = 0, b = 0;

	// The samples
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
}

void c2121a_poisson_mc_hier2_lev0::releaseL2Samples()
{
	int c = 0, l = 0, b = 0;

	if (mu_theta_samples) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					free(mu_theta_samples[c][l][b]);
				}
				free(mu_theta_samples[c][l]);
			}
			free(mu_theta_samples[c]);
		}
		free(mu_theta_samples);
		mu_theta_samples = NULL;
	}
	if (mu_gamma_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					free(mu_gamma_samples[c][l][b]);
				}
				free(mu_gamma_samples[c][l]);
			}
			free(mu_gamma_samples[c]);
		}
		free(mu_gamma_samples);
		mu_gamma_samples = NULL;
	}
	if (sigma2_theta_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					free(sigma2_theta_samples[c][l][b]);
				}
				free(sigma2_theta_samples[c][l]);
			}
			free(sigma2_theta_samples[c]);
		}
		free(sigma2_theta_samples);
		sigma2_theta_samples = NULL;
	}
	if (sigma2_gamma_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					free(sigma2_gamma_samples[c][l][b]);
				}
				free(sigma2_gamma_samples[c][l]);
			}
			free(sigma2_gamma_samples[c]);
		}
		free(sigma2_gamma_samples);
		sigma2_gamma_samples = NULL;
	}
}

void c2121a_poisson_mc_hier2_lev0::initL1Samples()
{
	int c = 0, l = 0, b = 0, a = 0;

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
}

void c2121a_poisson_mc_hier2_lev0::releaseL1Samples()
{
	int c = 0, l = 0, b = 0, a = 0;

	if (gTheta_samples) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					for (a = 0; a < gNAE[l][b]; a++) {
						free(gTheta_samples[c][l][b][a]);
					}
					free(gTheta_samples[c][l][b]);
				}
				free(gTheta_samples[c][l]);
			}
			free(gTheta_samples[c]);
		}
		free(gTheta_samples);
		gTheta_samples = NULL;
	}
	if (gGamma_samples) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					for (a = 0; a < gNAE[l][b]; a++) {
						free(gGamma_samples[c][l][b][a]);
					}
					free(gGamma_samples[c][l][b]);
				}
				free(gGamma_samples[c][l]);
			}
			free(gGamma_samples[c]);
		}
		free(gGamma_samples);
		gGamma_samples = NULL;
	}

	if (gTheta_acc != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gMaxBs; b++) {
					free(gTheta_acc[c][l][b]);
				}
				free(gTheta_acc[c][l]);
			}
			free(gTheta_acc[c]);
		}
		free(gTheta_acc);
		gTheta_acc = NULL;
	}
	if (gGamma_acc != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gMaxBs; b++) {
					free(gGamma_acc[c][l][b]);
				}
				free(gGamma_acc[c][l]);
			}
			free(gGamma_acc[c]);
		}
		free(gGamma_acc);
		gGamma_acc = NULL;
	}
}

c2121a_poisson_mc_hier2_lev0::~c2121a_poisson_mc_hier2_lev0()
{
	//Rprintf("c2121a_poisson_mc_hier2_lev0::c2121a_poisson_mc_hier2_lev0 - destructor\n");
	release();
}

void c2121a_poisson_mc_hier2_lev0::gibbs_sampler()
{
	if (strcmp(sim_type, "MH") == 0) {
		simulate_MH();
	}
	else {
		simulate_SLICE();
	}

	return;
}

void c2121a_poisson_mc_hier2_lev0::simulate_MH()
{
	int i = 0;

	for (i = 0; i < gIter; i++) {
#ifndef INDIVIDUAL_RNG
		GetRNGstate();
#endif
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

void c2121a_poisson_mc_hier2_lev0::simulate_SLICE()
{
	int i = 0;

	for (i = 0; i < gIter; i++) {
#ifndef INDIVIDUAL_RNG
		GetRNGstate();
#endif
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

void c2121a_poisson_mc_hier2_lev0::sample_mu_gamma(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {
				double denom = sigma2_gamma[c][l][b] + ((double)gNAE[l][b])*tau2_gamma_0;

				double t = 0.0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					t += gGamma[c][l][b][j];
				}

				double mean = (sigma2_gamma[c][l][b] * mu_gamma_0 + tau2_gamma_0 * t)/denom;

				double var = (sigma2_gamma[c][l][b]*tau2_gamma_0)/denom;

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

void c2121a_poisson_mc_hier2_lev0::sample_mu_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {
				double denom = sigma2_theta[c][l][b] + ((double)gNAE[l][b])*tau2_theta_0;

				double t = 0.0;
				int j = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					t += gTheta[c][l][b][j];
				}

				double mean = (sigma2_theta[c][l][b] * mu_theta_0 + tau2_theta_0 * t)/denom;

				double var = (sigma2_theta[c][l][b]*tau2_theta_0)/denom;

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

void c2121a_poisson_mc_hier2_lev0::sample_sigma2_gamma(int burnin, int iter)
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
					t += (pow(gGamma[c][l][b][j] - mu_gamma[c][l][b], 2.0));
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

void c2121a_poisson_mc_hier2_lev0::sample_sigma2_theta(int burnin, int iter)
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
					t += (pow((gTheta[c][l][b][j] - mu_theta[c][l][b]), 2.0));
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

double c2121a_poisson_mc_hier2_lev0::log_f_gamma(int c, int i, int b, int j, double gamm)
{
	double f1 = 0.0, f2 = 0.0, f3 = 0.0, f4 = 0.0, f5 = 0.0;

	f1 = ((double)x[i][b][j]) * gamm;
	f2 = -(exp(gamm)) * ((double)C[i][b][j]);
	f3 = ((double)y[i][b][j]) * (gamm + gTheta[c][i][b][j]);
	f4 = -(exp(gamm + gTheta[c][i][b][j]))*((double)T[i][b][j]);
	f5 = -(pow((gamm - mu_gamma[c][i][b]), 2.0))/(2.0 * sigma2_gamma[c][i][b]);

	double f = f1 + f2 + f3 + f4 + f5;

	return(f);
}

void c2121a_poisson_mc_hier2_lev0::sample_gamma_MH(int burnin, int iter)
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

void c2121a_poisson_mc_hier2_lev0::sample_gamma_SLICE(int burnin, int iter)
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

double c2121a_poisson_mc_hier2_lev0::log_f_theta(int c, int i, int b, int j, double theta)
{
	double f1 = 0.0, f2 = 0.0, f3 = 0.0;

	f1 = ((double)y[i][b][j]) * (gGamma[c][i][b][j] + theta);
	f2 = -(exp(gGamma[c][i][b][j] + theta)) * ((double)T[i][b][j]);
	f3 = - ((pow(theta - mu_theta[c][i][b], 2.0)))/(2.0 * sigma2_theta[c][i][b]);

	double f = f1 + f2 + f3;

	return(f);
}

void c2121a_poisson_mc_hier2_lev0::sample_theta_MH(int burnin, int iter)
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
					double cand = rnorm(gTheta[c][l][b][j], gSigma_MH_theta[l][b][j]);
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

					double f1 = log_f_theta(c, l, b, j, cand);
					double f2 = log_f_theta(c, l, b, j, gTheta[c][l][b][j]);

					double ratio = exp(f1 - f2);

					ratio = cMIN(ratio, 1);

					if (u <= ratio) {
						gTheta[c][l][b][j] = cand;
						gTheta_acc[c][l][b][j] = gTheta_acc[c][l][b][j] + 1;
					}

					if (iter >= burnin && retainSamples(iMonitor_theta)) {
						gTheta_samples[c][l][b][j][iter - burnin] = gTheta[c][l][b][j];
					}
				}
			}
		}
	}
}

void c2121a_poisson_mc_hier2_lev0::sample_theta_SLICE(int burnin, int iter)
{
	int c = 0, i = 0;
	//int m = gSim_Param_cntrl,
	int K = 0, J = 0;

	for (c = 0; c < gChains; c++) {
		for (i = 0; i < gNumIntervals; i++) {


			int b = 0, j = 0;
			double cand = 0.0;

			for (b = 0; b < gNumBodySys[i]; b++) {
				for (j = 0; j < gNAE[i][b]; j++) {

					int m = gW_theta_control[i][b][j];
#ifdef INDIVIDUAL_RNG
					GetRNGstate();
#endif
					J = floor(runif(0,m));
#ifdef INDIVIDUAL_RNG
					PutRNGstate();
#endif
					K = (m-1) - J;

					double l = 0.0, r = 0.0;
					double g = log_f_theta(c, i, b, j, gTheta[c][i][b][j]);
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
					double u = runif(0, gW_theta[i][b][j]);
#ifdef INDIVIDUAL_RNG
					PutRNGstate();
#endif

					l = gTheta[c][i][b][j] - u;
					r = gTheta[c][i][b][j] + (gW_theta[i][b][j] - u);

					while (J > 0) {
						if (logy >= log_f_theta(c, i, b, j, l)) {
							break;
						}
						l = l - gW_theta[i][b][j];
						J--;
					}

					while (K > 0) {
						if (logy >= log_f_theta(c, i, b, j, r)) {
							break;
						}
						r = r + gW_theta[i][b][j];
						K--;
					}
			
#ifdef INDIVIDUAL_RNG
					GetRNGstate();
#endif
					cand = runif(l, r);
#ifdef INDIVIDUAL_RNG
					PutRNGstate();
#endif

					while (logy >= log_f_theta(c, i, b, j, cand)) {
						if (cand < gTheta[c][i][b][j]) {
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

					gTheta[c][i][b][j] = cand;

					if (iter >= burnin && retainSamples(iMonitor_theta)) {
						gTheta_samples[c][i][b][j][iter - burnin] = gTheta[c][i][b][j];
					}
				}
			}
		}
	}
}

double c2121a_poisson_mc_hier2_lev0::cMIN(double a, double b)
{
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

void c2121a_poisson_mc_hier2_lev0::release()
{
	releaseGlobalSimParams();

	releaseDataVariables();

	releaseL1Variables();

	releaseL2Variables();

	releaseL2Samples();

	releaseL1Samples();

	releaseSimParams();

	releaseBaselineVariables();
}

SEXP c2121a_poisson_mc_hier2_lev0::getL1Samples(double***** &data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = allocVector(REALSXP, gChains * gNumIntervals * gMaxBs * gMaxAEs * (gIter - gBurnin)));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		int l = 0;
		for (l = 0; l < gNumIntervals; l++) {
			int b = 0;
			for (b = 0; b < gMaxBs; b++) {
				int j = 0;
				if (b < gNumBodySys[l]) {
					for (j = 0; j < gMaxAEs; j++) {
						if (j < gNAE[l][b]) {
							memcpy(REAL(samples) + i, data[c][l][b][j],
														(gIter - gBurnin)*sizeof(double));
							free(data[c][l][b][j]);
							data[c][l][b][j] = NULL;
						}
						i += (gIter - gBurnin);
					}
					free(data[c][l][b]);
					data[c][l][b] = NULL;
				}
				else {
					i += gMaxAEs*(gIter - gBurnin);
				}
			}
			free(data[c][l]);
			data[c][l] = NULL;
		}
		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = allocVector(INTSXP, 5));

	INTEGER(dim)[0] = (gIter - gBurnin);
	INTEGER(dim)[1] = gMaxAEs;
	INTEGER(dim)[2] = gMaxBs;
	INTEGER(dim)[3] = gNumIntervals;
	INTEGER(dim)[4] = gChains;

	setAttrib(samples, R_DimSymbol, dim);

	UNPROTECT(2);

	return samples;
}

SEXP c2121a_poisson_mc_hier2_lev0::getL2Samples(double**** &data)
{
	SEXP samples = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(samples = allocVector(REALSXP, gChains * gNumIntervals * gMaxBs * (gIter - gBurnin)));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		int l = 0;
		for (l = 0; l < gNumIntervals; l++) {
			int b = 0;
			for (b = 0; b < gMaxBs; b++) {
				if (b < gNumBodySys[l]) {
					memcpy(REAL(samples) + i, data[c][l][b],
														(gIter - gBurnin)*sizeof(double));
				}
				i += (gIter - gBurnin);
				free(data[c][l][b]);
				data[c][l][b] = NULL;
			}
			free(data[c][l]);
			data[c][l] = NULL;
		}
		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = allocVector(INTSXP, 4));

	INTEGER(dim)[0] = (gIter - gBurnin);
	INTEGER(dim)[1] = gMaxBs;
	INTEGER(dim)[2] = gNumIntervals;
	INTEGER(dim)[3] = gChains;

	setAttrib(samples, R_DimSymbol, dim);

	UNPROTECT(2);

	return samples;
}

SEXP c2121a_poisson_mc_hier2_lev0::getThetaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL1Samples(gTheta_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier2_lev0::getGammaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL1Samples(gGamma_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier2_lev0::getMuThetaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(mu_theta_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier2_lev0::getMuGammaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(mu_gamma_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier2_lev0::getSigma2ThetaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(sigma2_theta_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier2_lev0::getSigma2GammaSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(sigma2_gamma_samples);

	return samples;
}

SEXP c2121a_poisson_mc_hier2_lev0::getL1Accept(int**** &data)
{
	SEXP acc = R_NilValue;
	SEXP dim = R_NilValue;

   PROTECT(acc = allocVector(INTSXP, gChains * gNumIntervals * gMaxBs * gMaxAEs));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		int l = 0;
		for (l = 0; l < gNumIntervals; l++) {
			int b = 0;
			for (b = 0; b < gMaxBs; b++) {
				if (b < gNumBodySys[l]) {
					memcpy(INTEGER(acc) + i, data[c][l][b], gMaxAEs*sizeof(int));
				}
				i += gMaxAEs;
				free(data[c][l][b]);
				data[c][l][b] = NULL;
			}
			free(data[c][l]);
			data[c][l] = NULL;
		}
		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = allocVector(INTSXP, 4));

	INTEGER(dim)[0] = gMaxAEs;
	INTEGER(dim)[1] = gMaxBs;
	INTEGER(dim)[2] = gNumIntervals;
	INTEGER(dim)[3] = gChains;

	setAttrib(acc, R_DimSymbol, dim);

	UNPROTECT(2);

	return acc;
}

SEXP c2121a_poisson_mc_hier2_lev0::getThetaAccept()
{
	SEXP acc = R_NilValue;

	acc = getL1Accept(gTheta_acc);

	return acc;
}

SEXP c2121a_poisson_mc_hier2_lev0::getGammaAccept()
{
	SEXP acc = R_NilValue;

	acc = getL1Accept(gGamma_acc);

	return acc;
}

void c2121a_poisson_mc_hier2_lev0::getThetaSamples(int *c, int*l, int* b, int* j, double* theta_samples)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

	if (gTheta_samples)
		memcpy(theta_samples,
					gTheta_samples[C][L][B][J], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier2_lev0::getGammaSamples(int *c, int *l, int* b, int* j, double* gamma_samples)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

	if (gGamma_samples)
		memcpy(gamma_samples,
					gGamma_samples[C][L][B][J], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier2_lev0::getMuThetaSamples(int *c, int *l, int* b, double* mu_theta)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;

	if (mu_theta_samples)
		memcpy(mu_theta, mu_theta_samples[C][L][B], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier2_lev0::getMuGammaSamples(int *c, int *l, int* b, double* mu_gamma)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;

	if (mu_gamma_samples)
		memcpy(mu_gamma, mu_gamma_samples[C][L][B], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier2_lev0::getSigma2ThetaSamples(int *c, int *l, int* b, double* sigma2)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;

	if (sigma2_theta_samples)
		memcpy(sigma2, sigma2_theta_samples[C][L][B], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier2_lev0::getSigma2GammaSamples(int *c, int *l, int* b, double* sigma2)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;

	if (sigma2_gamma_samples)
		memcpy(sigma2, sigma2_gamma_samples[C][L][B], (gIter - gBurnin)*sizeof(double));
}

void c2121a_poisson_mc_hier2_lev0::getThetaAccept(int *c, int *l, int* b, int* j, double* theta_acc)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

	*theta_acc = gTheta_acc[C][L][B][J];
}

void c2121a_poisson_mc_hier2_lev0::getGammaAccept(int *c, int *l, int* b, int* j, double* gamma_acc)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

	*gamma_acc = gGamma_acc[C][L][B][J];
}
