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
#include "c2121a_poisson_mc_hier3_lev0.h"
#include "c212BB_poisson_mc_hier3_lev0.h"

using namespace std;

static const char *rcsId = "$Id: c212BB_poisson_mc_hier3_lev0.cpp,v 1.16 2017/03/22 16:12:09 clb13102 Exp clb13102 $";

const char* c212BB_poisson_mc_hier3_lev0::sMonitor_pi = "pi";
const char* c212BB_poisson_mc_hier3_lev0::sMonitor_alpha_pi = "alpha.pi";
const char* c212BB_poisson_mc_hier3_lev0::sMonitor_beta_pi = "beta.pi";

const char* c212BB_poisson_mc_hier3_lev0::sColType = "type";
const char* c212BB_poisson_mc_hier3_lev0::sColParam = "param";
const char* c212BB_poisson_mc_hier3_lev0::sColValue = "value";
const char* c212BB_poisson_mc_hier3_lev0::sColControl = "control";

const char* c212BB_poisson_mc_hier3_lev0::sColPMweight = "weight_pm";

const char* c212BB_poisson_mc_hier3_lev0::sParam_sigma_MH_gamma = "sigma_MH_gamma";
const char* c212BB_poisson_mc_hier3_lev0::sParam_sigma_MH_theta = "sigma_MH_theta";
const char* c212BB_poisson_mc_hier3_lev0::sParam_sigma_MH_alpha = "sigma_MH_alpha";
const char* c212BB_poisson_mc_hier3_lev0::sParam_sigma_MH_beta = "sigma_MH_beta";
const char* c212BB_poisson_mc_hier3_lev0::sParam_w_gamma = "w_gamma";
const char* c212BB_poisson_mc_hier3_lev0::sParam_w_theta = "w_theta";
const char* c212BB_poisson_mc_hier3_lev0::sParam_w_alpha = "w_alpha";
const char* c212BB_poisson_mc_hier3_lev0::sParam_w_beta = "w_beta";

const char* c212BB_poisson_mc_hier3_lev0::sVariable_alpha = "alpha";
const char* c212BB_poisson_mc_hier3_lev0::sVariable_beta = "beta";


c212BB_poisson_mc_hier3_lev0::c212BB_poisson_mc_hier3_lev0()
{
	//Rprintf("c212BB_poisson_mc_hier3_lev0::c212BB_poisson_mc_hier3_lev0: Default constructor\n");

	gW_alpha = NULL;
	gW_beta = NULL;
	gW_alpha_control = NULL;
	gW_beta_control = NULL;
	gSigma_MH_alpha = NULL;
	gSigma_MH_beta = NULL;

	gWp = NULL;
	gMH_weight = 0.5;

	iMonitor_pi = 0;
	iMonitor_alpha_pi = 0;
	iMonitor_beta_pi = 0;

	gSimType = eSim_Type_SLICE;

	gDefault_Sigma_MH_alpha = 1.0;
	gDefault_Sigma_MH_beta = 1.0;
	gDefault_Sigma_MH_gamma = 1.0;
	gDefault_Sigma_MH_theta = 1.0;
	gDefault_W_alpha = 1.0;
	gDefault_W_beta = 1.0;
	gDefault_W_gamma = 1.0;
	gDefault_W_alpha_control = 6.0;
	gDefault_W_beta_control = 6.0;
	gDefault_W_gamma_control = 6.0;

	lambda_alpha = 0.0;
	lambda_beta = 0.0;

	alpha_pi = NULL;
	beta_pi = NULL;

	gPi = NULL;

	alpha_pi_acc = NULL;
	beta_pi_acc = NULL;

	gPi_samples = NULL;
	alpha_pi_samples = NULL;
	beta_pi_samples = NULL;
}

c212BB_poisson_mc_hier3_lev0::c212BB_poisson_mc_hier3_lev0(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model, SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
                    SEXP MH_weight,
                    SEXP pm_weights,
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
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi, SEXP plambda_alpha,
					SEXP plambda_beta,
					SEXP palgo, SEXP padapt_phase)

{
	gW_alpha = NULL;
	gW_beta = NULL;
	gW_alpha_control = NULL;
	gW_beta_control = NULL;
	gSigma_MH_alpha = NULL;
	gSigma_MH_beta = NULL;

	gWp = NULL;
	gMH_weight = 0.5;

	iMonitor_pi = 0;
	iMonitor_alpha_pi = 0;
    iMonitor_beta_pi = 0;

	gSimType = eSim_Type_SLICE;

	gDefault_Sigma_MH_alpha = 1.0;
	gDefault_Sigma_MH_beta = 1.0;
	gDefault_Sigma_MH_gamma = 1.0;
	gDefault_Sigma_MH_theta = 1.0;
	gDefault_W_alpha = 1.0;
	gDefault_W_beta = 1.0;
	gDefault_W_gamma = 1.0;
	gDefault_W_alpha_control = 6.0;
	gDefault_W_beta_control = 6.0;
	gDefault_W_gamma_control = 6.0;

	lambda_alpha = 0.0;
	lambda_beta = 0.0;

	alpha_pi = NULL;
	beta_pi = NULL;

	gPi = NULL;

	alpha_pi_acc = NULL;
	beta_pi_acc = NULL;

	gPi_samples = NULL;
	alpha_pi_samples = NULL;
	beta_pi_samples = NULL;

	init(sChains, sBurnin, sIter, sSim_Type, sMem_Model, sGlobal_Sim_Params,
				sSim_Params, MH_weight, pm_weights,
				sMonitor,
				sNumIntervals, sMaxBs, sNumBodySys, sMaxAEs, sNAE,
				pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0,
				pmu_theta_0_0, ptau2_theta_0_0,
				palpha_gamma_0_0, pbeta_gamma_0_0, palpha_theta_0_0, pbeta_theta_0_0,
				palpha_gamma,
				pbeta_gamma, palpha_theta, pbeta_theta, pmu_gamma_0, ptau2_gamma_0,
				pmu_theta_0,
				ptau2_theta_0, pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
				pPi, palpha_pi, pbeta_pi, plambda_alpha, plambda_beta,
				palgo, padapt_phase);
}

void c212BB_poisson_mc_hier3_lev0::initPMWeights(SEXP pm_weights)
{
	gWp = (double ***)malloc(gNumIntervals * sizeof(double**));

	int i = 0, b = 0, j = 0;
	for (i = 0; i < gNumIntervals; i++) {
		gWp[i] = (double**)malloc(gNumBodySys[i] * sizeof(double*));
		for (b = 0; b < gNumBodySys[i]; b++) {
			gWp[i][b] = (double*)malloc(gNAE[i][b] * sizeof(double));
			for (j = 0; j < gNAE[i][b]; j++) {
				gWp[i][b][j] = gMH_weight;
			}
		}
	}

	int len = Rf_length(pm_weights);

	SEXP sPM_Weights = R_NilValue;
	SEXP sI_index = R_NilValue;
	SEXP sB = R_NilValue;
	SEXP sj = R_NilValue;

	if (len && isNewList(pm_weights)) {

		SEXP names = getAttrib(pm_weights, R_NamesSymbol);

		for (i = 0; i < len; i++) {
			if (strcmp(sColPMweight, CHAR(STRING_ELT(names, i))) == 0) {
				sPM_Weights = VECTOR_ELT(pm_weights, i);
			}
			if (strcmp(sColI_index, CHAR(STRING_ELT(names, i))) == 0) {
				sI_index = VECTOR_ELT(pm_weights, i);
			}
			if (strcmp(sColB, CHAR(STRING_ELT(names, i))) == 0) {
				sB = VECTOR_ELT(pm_weights, i);
			}
			if (strcmp(sColj, CHAR(STRING_ELT(names, i))) == 0) {
				sj = VECTOR_ELT(pm_weights, i);
			}
		}

		len = Rf_length(sPM_Weights);

		if (len > 0) {
			double* weights = REAL(sPM_Weights);
			int* L = INTEGER(sI_index);
			int* B = INTEGER(sB);
			int* j = INTEGER(sj);

			for (i = 0; i < len; i++) {
				int l = L[i] - 1;
				int b = B[i] - 1;
				int a = j[i] - 1;
				gWp[l][b][a] = weights[i];
			}
		}
	}
}

void c212BB_poisson_mc_hier3_lev0::initSimParams(SEXP sSim_Params)
{
	gW_gamma = (double ***)malloc(gNumIntervals * sizeof(double**));
	gW_gamma_control = (int ***)malloc(gNumIntervals * sizeof(int**));
	gSigma_MH_gamma = (double ***)malloc(gNumIntervals * sizeof(double**));
	gSigma_MH_theta = (double ***)malloc(gNumIntervals * sizeof(double**));
	gW_alpha = (double *)malloc(gNumIntervals * sizeof(double));
	gW_beta = (double *)malloc(gNumIntervals * sizeof(double));
	gW_alpha_control = (double *)malloc(gNumIntervals * sizeof(double));
	gW_beta_control = (double *)malloc(gNumIntervals * sizeof(double));
	gSigma_MH_alpha = (double *)malloc(gNumIntervals * sizeof(double));
	gSigma_MH_beta = (double *)malloc(gNumIntervals * sizeof(double));

	int i = 0, b = 0, j = 0;

	for (i = 0; i < gNumIntervals; i++) {

		gW_gamma[i] = (double**)malloc(gNumBodySys[i] * sizeof(double*));
		gW_gamma_control[i] = (int**)malloc(gNumBodySys[i] * sizeof(int*));
		gSigma_MH_gamma[i] = (double**)malloc(gNumBodySys[i] * sizeof(double*));
		gSigma_MH_theta[i] = (double**)malloc(gNumBodySys[i] * sizeof(double*));

		gW_alpha[i] = gDefault_W_alpha;
		gW_beta[i] = gDefault_W_alpha;
		gW_alpha_control[i] = gDefault_W_alpha_control;
		gW_beta_control[i] = gDefault_W_beta_control;
		gSigma_MH_alpha[i] = gDefault_Sigma_MH_alpha;
		gSigma_MH_beta[i] = gDefault_Sigma_MH_beta;

		for (b = 0; b < gNumBodySys[i]; b++) {

			gW_gamma[i][b] = (double*)malloc(gNAE[i][b] * sizeof(double));
			gW_gamma_control[i][b] = (int*)malloc(gNAE[i][b] * sizeof(int));
			gSigma_MH_gamma[i][b] = (double*)malloc(gNAE[i][b] * sizeof(double));
			gSigma_MH_theta[i][b] = (double*)malloc(gNAE[i][b] * sizeof(double));

			for (j = 0; j < gNAE[i][b]; j++) {
				gW_gamma[i][b][j] = gDefault_W_gamma;
				gW_gamma_control[i][b][j] = (int)gDefault_W_gamma_control;
				gSigma_MH_gamma[i][b][j] = gDefault_Sigma_MH_gamma;
				gSigma_MH_theta[i][b][j] = gDefault_Sigma_MH_theta;
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
					if (0 == strcmp(param, sParam_w_gamma)) {
						gW_gamma[l][b][a] = vals[i];
						gW_gamma_control[l][b][a] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH_gamma)) {
						gSigma_MH_gamma[l][b][a] = vals[i];
					}
				}
				else if (0 == strcmp(sVariable_theta, var)) {
					if (0 == strcmp(param, sParam_w_theta)) {
						gW_theta[l][b][a] = vals[i];
						gW_theta_control[l][b][a] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH_theta)) {
						gSigma_MH_theta[l][b][a] = vals[i];
					}
				}
				else if (0 == strcmp(sVariable_alpha, var)) {
					if (0 == strcmp(param, sParam_w_alpha)) {
						gW_alpha[l] = vals[i];
						gW_alpha_control[l] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH_alpha)) {
						gSigma_MH_alpha[l] = vals[i];
					}
				}
				else if (0 == strcmp(sVariable_beta, var)) {
					if (0 == strcmp(param, sParam_w_beta)) {
						gW_beta[l] = vals[i];
						gW_beta_control[l] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH_beta)) {
						gSigma_MH_beta[l] = vals[i];
					}
				}
			}
		}
	}
}

void c212BB_poisson_mc_hier3_lev0::initGlobalSimParams(SEXP sSim_Type, SEXP sGlobal_Sim_Params)
{
	int len = Rf_length(sGlobal_Sim_Params);

	SEXP sParams = R_NilValue;
	SEXP sValues = R_NilValue;
	SEXP sControl = R_NilValue;

	if (strcmp("MH", CHAR(STRING_ELT(sSim_Type, 0))) == 0) {
		gSimType = eSim_Type_MH;
	}
	else {
		gSimType = eSim_Type_SLICE;
	}

    if (len > 0 && isNewList(sGlobal_Sim_Params)) {

        SEXP names = getAttrib(sGlobal_Sim_Params, R_NamesSymbol);

        int i = 0;

        for (i = 0; i < len; i++) {
            if (strcmp(sColValue, CHAR(STRING_ELT(names, i))) == 0) {
                sValues = VECTOR_ELT(sGlobal_Sim_Params, i);
            }
            //if (strcmp(sColType, CHAR(STRING_ELT(names, i))) == 0) {
            //  sType = VECTOR_ELT(sGlobal_Sim_Params, i);
            //}
            if (strcmp(sColParam, CHAR(STRING_ELT(names, i))) == 0) {
                sParams = VECTOR_ELT(sGlobal_Sim_Params, i);
            }
            if (strcmp(sColControl, CHAR(STRING_ELT(names, i))) == 0) {
                sControl = VECTOR_ELT(sGlobal_Sim_Params, i);
            }
        }

        len = Rf_length(sParams);

		if (len > 0) {

			double* vals = REAL(sValues);
			double* cntrl = REAL(sControl);

			for (i = 0; i < len; i++) {
				const char *t = CHAR(STRING_ELT(sParams, i));

				if (0 == strcmp(t, sParam_sigma_MH_gamma)) {
					gDefault_Sigma_MH_gamma = vals[i];
				}
				if (0 == strcmp(t, sParam_sigma_MH_theta)) {
					gDefault_Sigma_MH_theta = vals[i];
				}
				if (0 == strcmp(t, sParam_sigma_MH_alpha)) {
					gDefault_Sigma_MH_alpha = vals[i];
				}
				if (0 == strcmp(t, sParam_sigma_MH_beta)) {
					gDefault_Sigma_MH_beta = vals[i];
				}
				if (0 == strcmp(t, sParam_w_gamma)) {
					gDefault_W_gamma = vals[i];
					gDefault_W_gamma_control = cntrl[i];
				}
				if (0 == strcmp(t, sParam_w_alpha)) {
					gDefault_W_alpha = vals[i];
					gDefault_W_alpha_control = cntrl[i];
				}
				if (0 == strcmp(t, sParam_w_beta)) {
					gDefault_W_beta = vals[i];
					gDefault_W_beta_control = cntrl[i];
				}
			}
    	}
    }
}

void c212BB_poisson_mc_hier3_lev0::init(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type, SEXP sMem_Model, SEXP sGlobal_Sim_Params,
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

	alpha_pi = (double**)malloc(gChains * sizeof(double*));
	double *valpha_pi = REAL(palpha_pi);
	for (c = 0; c < gChains; c++) {
		alpha_pi[c] = (double *)malloc(gNumIntervals * sizeof(double));
		for (l = 0; l < gNumIntervals; l++) {
			alpha_pi[c][l] = *valpha_pi;
			valpha_pi++;
		}
	}

	beta_pi = (double**)malloc(gChains * sizeof(double*));
	double *vbeta_pi = REAL(pbeta_pi);
	for (c = 0; c < gChains; c++) {
		beta_pi[c] = (double *)malloc(gNumIntervals * sizeof(double));
		for (l = 0; l < gNumIntervals; l++) {
			beta_pi[c][l] = *vbeta_pi;
			vbeta_pi++;
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
		mu_gamma_0_samples = (double ***)malloc(gChains * sizeof(double**));
	if (retainSamples(iMonitor_mu_theta_0))
		mu_theta_0_samples = (double ***)malloc(gChains * sizeof(double**));
	if (retainSamples(iMonitor_tau2_gamma_0))
		tau2_gamma_0_samples = (double ***)malloc(gChains * sizeof(double**));
	if (retainSamples(iMonitor_tau2_theta_0))
		tau2_theta_0_samples = (double ***)malloc(gChains * sizeof(double**));
	if (retainSamples(iMonitor_alpha_pi))
		alpha_pi_samples = (double ***)malloc(gChains * sizeof(double**));
	if (retainSamples(iMonitor_beta_pi))
		beta_pi_samples = (double ***)malloc(gChains * sizeof(double**));

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
		if (retainSamples(iMonitor_alpha_pi))
			alpha_pi_samples[c] = (double **)malloc(gNumIntervals* sizeof(double *));
		if (retainSamples(iMonitor_beta_pi))
			beta_pi_samples[c] = (double **)malloc(gNumIntervals* sizeof(double *));

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
			if (retainSamples(iMonitor_alpha_pi))
				alpha_pi_samples[c][l] =
									(double *)malloc((gIter - gBurnin)* sizeof(double));
			if (retainSamples(iMonitor_beta_pi))
				beta_pi_samples[c][l] =
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
					mu_gamma_samples[c][l][b]
								= (double *)malloc((gIter - gBurnin)*sizeof(double));
				if (retainSamples(iMonitor_sigma2_theta))
					sigma2_theta_samples[c][l][b]
								= (double *)malloc((gIter - gBurnin)*sizeof(double));
				if (retainSamples(iMonitor_sigma2_gamma))
					sigma2_gamma_samples[c][l][b]
								= (double *)malloc((gIter - gBurnin)*sizeof(double));
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

	alpha_pi_acc = (int**)malloc(gChains * sizeof(int*));
	beta_pi_acc = (int**)malloc(gChains * sizeof(int*));
	for (c = 0; c < gChains; c++) {
		alpha_pi_acc[c] = (int*)malloc(gNumIntervals * sizeof(int));
		beta_pi_acc[c] = (int*)malloc(gNumIntervals * sizeof(int));
		for (l = 0; l < gNumIntervals; l++) {
			alpha_pi_acc[c][l] = 0;
			beta_pi_acc[c][l] = 0;
		}
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

void c212BB_poisson_mc_hier3_lev0::initMonitor(SEXP sMonitor)
{
    int len = Rf_length(sMonitor);

    SEXP sVariables = R_NilValue;
    SEXP sValues = R_NilValue;

    if (isNewList(sMonitor)) {

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
			if (0 == strcmp(t, sMonitor_pi)) {
				iMonitor_pi = vals[i];
			}
			if (0 == strcmp(t, sMonitor_alpha_pi)) {
				iMonitor_alpha_pi = vals[i];
			}
			if (0 == strcmp(t, sMonitor_beta_pi)) {
				iMonitor_beta_pi = vals[i];
			}
		}
	}
}

c212BB_poisson_mc_hier3_lev0::~c212BB_poisson_mc_hier3_lev0()
{
	//Rprintf("c212BB_poisson_mc_hier3_lev0::c212BB_poisson_mc_hier3_lev0 - destructor\n");
	release();
}

void c212BB_poisson_mc_hier3_lev0::gibbs_sampler()
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

void c212BB_poisson_mc_hier3_lev0::simulate_MH()
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

void c212BB_poisson_mc_hier3_lev0::simulate_SLICE()
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

double c212BB_poisson_mc_hier3_lev0::log_f_alpha_pi(int c, int l, double alpha)
{
	double f = 0.0;
	double log_pi_sum = 0.0;

	int b = 0;
	for (b = 0; b < gNumBodySys[l]; b++) {
		log_pi_sum += log(gPi[c][l][b]);
	}

	f = ((double)gNumBodySys[l]) * (lgammafn(alpha + beta_pi[c][l]) - lgammafn(alpha))
							+ (alpha - 1.0) * log_pi_sum - alpha * lambda_alpha; 
	return(f);
}

void c212BB_poisson_mc_hier3_lev0::sample_alpha_pi_MH(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			double cand = 0;

		    // alpha_pi is restricted to being greater than zero
		    while (cand <= 1.0) {
#ifdef INDIVIDUAL_RNG
		        GetRNGstate();
#endif
		        //cand = rnorm(alpha_pi[c][l], gDefault_Sigma_MH_alpha);
		        cand = rnorm(alpha_pi[c][l], gSigma_MH_alpha[l]);
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

			double f1 = log_f_alpha_pi(c, l, cand);
			double f2 = log_f_alpha_pi(c, l, alpha_pi[c][l]);

			double q1 = pnorm((alpha_pi[c][l] - 1)/gSigma_MH_alpha[l], 0, 1, 1, 0);
			double q2 = pnorm((cand - 1)/gSigma_MH_alpha[l], 0, 1, 1, 0);

			double ratio = (exp(f1 - f2))* q1/q2;
			ratio = cMIN(ratio, 1);

		    if (u <= ratio) {
		        alpha_pi[c][l] = cand;
				//if (iter >= burnin)
					alpha_pi_acc[c][l] = alpha_pi_acc[c][l] + 1;
			}

			if (iter >= burnin && retainSamples(iMonitor_alpha_pi)) {
				alpha_pi_samples[c][l][iter - burnin] = alpha_pi[c][l];
			}
		}
	}
}

void c212BB_poisson_mc_hier3_lev0::sample_alpha_pi_SLICE(int burnin, int iter)
{
	int c = 0, i = 0;
	int K = 0, J = 0;

	for (c = 0; c< gChains; c++) {
		for (i = 0; i < gNumIntervals; i++) {

			int m = gW_alpha_control[i];

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

			double g = log_f_alpha_pi(c, i , alpha_pi[c][i]);
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
			double u = runif(0, gW_alpha[i]);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			l = alpha_pi[c][i] - u;
			r = alpha_pi[c][i] + (gW_alpha[i] - u);

			while (J > 0) {
				if (l <= 1.0)
					break;

				if (logy >= log_f_alpha_pi(c, i, l)) {
					break;
				}
				l = l - gW_alpha[i];

				J--;
			}

			while (K > 0) {
				if (logy >= log_f_alpha_pi(c, i, r)) {
                            break;
				}
				r = r + gW_alpha[i];
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

			while (logy >= log_f_alpha_pi(c, i, cand)) {
				if (cand < alpha_pi[c][i]) {
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

			alpha_pi[c][i] = cand;

			if (iter >= burnin && retainSamples(iMonitor_alpha_pi)) {
				alpha_pi_samples[c][i][iter - burnin] = alpha_pi[c][i];
			}
		}
	}
}

double c212BB_poisson_mc_hier3_lev0::log_f_beta_pi(int c, int l, double beta)
{
    double f = 0.0;
    double log_sum = 0.0;

    int b = 0;
    for (b = 0; b < gNumBodySys[l]; b++) {
        log_sum += log(1.0 - gPi[c][l][b]);
    }

    f = ((double)gNumBodySys[l]) * (lgammafn(alpha_pi[c][l] + beta) - lgammafn(beta)) + (beta - 1.0) * log_sum - beta * lambda_beta; 
    return(f);
}

void c212BB_poisson_mc_hier3_lev0::sample_beta_pi_MH(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c< gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			double cand = 0.0;

			while (cand <= 1.0) {
#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				cand = rnorm(beta_pi[c][l], gSigma_MH_beta[l]);
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

			double f1 = log_f_beta_pi(c, l, cand);
			double f2 = log_f_beta_pi(c, l, beta_pi[c][l]);

			double q1 = pnorm((beta_pi[c][l] - 1)/gSigma_MH_beta[l], 0, 1, 1, 0);
			double q2 = pnorm((cand - 1)/gSigma_MH_beta[l], 0, 1, 1, 0);

			double ratio = (exp(f1 - f2)) * (q1/q2);

			ratio = cMIN(ratio, 1);

			if (u <= ratio) {
				beta_pi[c][l] = cand;
				//if (iter >= burnin)
					beta_pi_acc[c][l] = beta_pi_acc[c][l] + 1;
			}

			if (iter >= burnin && retainSamples(iMonitor_beta_pi)) {
				beta_pi_samples[c][l][iter - burnin] = beta_pi[c][l];
			}
		}
	}
}

void c212BB_poisson_mc_hier3_lev0::sample_beta_pi_SLICE(int burnin, int iter)
{
	int c = 0, i = 0;
	int K = 0, J = 0;

	for (c = 0; c< gChains; c++) {
		for (i = 0; i < gNumIntervals; i++) {

			int m = gW_beta_control[i];

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

			double g = log_f_beta_pi(c, i, beta_pi[c][i]);
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
			double u = runif(0, gW_beta[i]);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			l = beta_pi[c][i] - u;
			r = beta_pi[c][i] + (gW_beta[i] - u);

			// beta is retricted to being greater than 1
			// need a do - while loop
			while (J > 0) {
				if (l <= 1.0) {
					l = 1.0;
					break;
				}

				if (logy >= log_f_beta_pi(c, i, l)) {
					break;
				}
				l = l - gW_beta[i];
				J--;
			}

			while (K > 0) {
				if (logy >= log_f_beta_pi(c, i, r)) {
					break;
				}
				r = r + gW_beta[i];
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

			while (logy >= log_f_beta_pi(c, i, cand)) {
				if (cand < beta_pi[c][i]) {
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

			beta_pi[c][i] = cand;

			if (iter >= burnin && retainSamples(iMonitor_beta_pi)) {
				beta_pi_samples[c][i][iter - burnin] = beta_pi[c][i];
			}
		}
	}
}

void c212BB_poisson_mc_hier3_lev0::sample_pi(int burnin, int iter)
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

				double shape1 = alpha_pi[c][l] + (double)theta_zero_count;
				double shape2 = beta_pi[c][l] + (double)gNAE[l][b] - (double)theta_zero_count;

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

void c212BB_poisson_mc_hier3_lev0::sample_mu_theta(int burnin, int iter)
{
	int c = 0, l = 0;

	for (c = 0; c < gChains; c++) {
		for (l = 0; l < gNumIntervals; l++) {

			int b = 0;

			for (b = 0; b < gNumBodySys[l]; b++) {

				double t = 0.0;
				int j = 0;
				int Kb = 0;
				for (j = 0; j < gNAE[l][b]; j++) {
					if (gTheta[c][l][b][j] != 0.0) {
						Kb++;
					}
					t += gTheta[c][l][b][j];
				}

				double denom = sigma2_theta[c][l][b] + ((double)Kb)*tau2_theta_0[c][l];

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

void c212BB_poisson_mc_hier3_lev0::sample_sigma2_theta(int burnin, int iter)
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

double c212BB_poisson_mc_hier3_lev0::log_q_theta(int l, int b, int j, double p, double theta, double mean)
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


double c212BB_poisson_mc_hier3_lev0::log_f_theta(int c, int i, int b, int j, double theta)
{
	double f1 = 0.0, f2 = 0.0;

	f1 = (((double)y[i][b][j]) * theta) - (exp(gGamma[c][i][b][j] + theta)) * ((double)T[i][b][j]);

	if (theta == 0) {
		f2 = log(gPi[c][i][b]);
	}
	else {
		f2 = log(1 - gPi[c][i][b]) + log(1.0/sqrt(2.0 * M_PI*sigma2_theta[c][i][b]))
				+ ((-1.0/2.0)*(pow(theta -mu_theta[c][i][b], 2.0))/sigma2_theta[c][i][b]);
	}

	double f = f1 + f2;

	return(f);
}

/*
* Sample theta using a MH step as detailed in: Gottardo, Raftery - Markov Chain Monte Carlo
* With Mixtures of Mutually Singular Distributions
* Perform an adaption step to get a good candidate density
*/
void c212BB_poisson_mc_hier3_lev0::sample_theta_MH(int burnin, int iter)
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
					}

					if (iter >= burnin && retainSamples(iMonitor_theta)) {
						gTheta_samples[c][l][b][j][iter - burnin] = gTheta[c][l][b][j];
					}
				}
			}
		}
	}
}

double c212BB_poisson_mc_hier3_lev0::cMIN(double a, double b)
{
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

void c212BB_poisson_mc_hier3_lev0::release()
{
	int c = 0, l = 0, b = 0;

	if (alpha_pi != NULL) {
		for (c = 0; c < gChains; c++) {
			free(alpha_pi[c]);
		}
		free(alpha_pi);
		alpha_pi = NULL;
	}

	if (alpha_pi_acc != NULL) {
		for (c = 0; c < gChains; c++) {
			free(alpha_pi_acc[c]);
		}

		free(alpha_pi_acc);
		alpha_pi_acc = NULL;
	}

	if (beta_pi != NULL) {
		for (c = 0; c < gChains; c++) {
			free(beta_pi[c]);
		}
		free(beta_pi);
		beta_pi = NULL;
	}

	if (beta_pi_acc != NULL) {
		for (c = 0; c < gChains; c++) {
			free(beta_pi_acc[c]);
		}

		free(beta_pi_acc);
		beta_pi_acc = NULL;
	}

	if (gPi != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(gPi[c][l]);
			}
			free(gPi[c]);
		}
		free(gPi);
		gPi = 0;
	}

	if (alpha_pi_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(alpha_pi_samples[c][l]);
			}
			free(alpha_pi_samples[c]);
		}
		free(alpha_pi_samples);
		alpha_pi_samples = NULL;
	}

	if (beta_pi_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				free(beta_pi_samples[c][l]);
			}
			free(beta_pi_samples[c]);
		}
		free(beta_pi_samples);
		beta_pi_samples = NULL;
	}

	if (gPi_samples) {
		for (c = 0; c < gChains; c++) {
			for (l = 0; l < gNumIntervals; l++) {
				for (b = 0; b < gNumBodySys[l]; b++) {
					free(gPi_samples[c][l][b]);
				}
				free(gPi_samples[c][l]);
			}
			free(gPi_samples[c]);
		}
		free(gPi_samples);
		gPi_samples = NULL;
	}

	if (gWp != NULL) {

		for (l = 0; l < gNumIntervals; l++) {
			for (b = 0; b < gNumBodySys[l]; b++) {
				free(gWp[l][b]);
			}
			free(gWp[l]);
		}

		free(gWp);
		gWp = NULL;
	}

	if (gW_alpha != NULL) {
		free(gW_alpha);
		gW_alpha = NULL;
	}
	if (gW_beta != NULL) {
		free(gW_beta);
		gW_beta = NULL;
	}
	if (gW_alpha_control != NULL) {
		free(gW_alpha_control);
		gW_alpha_control = NULL;
	}
	if (gW_beta_control != NULL) {
		free(gW_beta_control);
		gW_beta_control = NULL;
	}
	if (gSigma_MH_alpha != NULL) {
		free(gSigma_MH_alpha);
		gSigma_MH_alpha = NULL;
	}
	if (gSigma_MH_beta != NULL) {
		free(gSigma_MH_beta);
		gSigma_MH_beta = NULL;
	}
}

SEXP c212BB_poisson_mc_hier3_lev0::getPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(gPi_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev0::getAlphaPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(alpha_pi_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev0::getBetaPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(beta_pi_samples);

	return samples;
}

SEXP c212BB_poisson_mc_hier3_lev0::getAlphaPiAccept()
{
	SEXP acc = R_NilValue;

	acc = getL3Accept(alpha_pi_acc);

	return acc;
}

SEXP c212BB_poisson_mc_hier3_lev0::getL3Accept(int** &data)
{
	SEXP acc = R_NilValue;
	SEXP dim = R_NilValue;

	PROTECT(acc = allocVector(INTSXP, gChains * gNumIntervals));

	int i = 0;
	int c = 0;
	for (c = 0; c < gChains; c++) {
		memcpy(INTEGER(acc) + i, data[c], gNumIntervals * sizeof(int));

		//int l = 0;
		//for (l = 0; l < gNumIntervals; l++) {
		//	//memcpy(INTEGER(acc) + i, data[c][l], sizeof(int));
		//	(INTEGER(acc)[i] = data[c][l];
		//	i++;
		//}

		free(data[c]);
		data[c] = NULL;
	}
	free(data);
	data = NULL;

	PROTECT(dim = allocVector(INTSXP, 2));

	INTEGER(dim)[0] = gNumIntervals;
	INTEGER(dim)[1] = gChains;

	setAttrib(acc, R_DimSymbol, dim);

	UNPROTECT(2);

	return acc;
}

SEXP c212BB_poisson_mc_hier3_lev0::getBetaPiAccept()
{
	SEXP acc = R_NilValue;

	acc = getL3Accept(beta_pi_acc);

	return acc;
}

void c212BB_poisson_mc_hier3_lev0::getPiSamples(int *c, int *l, int* b, double* pi)
{
	int C = (*c) - 1;
	int L = (*l) - 1;
	int B = (*b) - 1;

	memcpy(pi, gPi_samples[C][L][B], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev0::getAlphaPiSamples(int *c, int *l, double* alpha_pi)
{
	int C = (*c) - 1;
	int L = (*l) - 1;

	memcpy(alpha_pi, alpha_pi_samples[C][L], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev0::getBetaPiSamples(int *c, int *l, double* beta_pi)
{
	int C = (*c) - 1;
	int L = (*l) - 1;

	memcpy(beta_pi, beta_pi_samples[C][L], (gIter - gBurnin)*sizeof(double));
}

void c212BB_poisson_mc_hier3_lev0::getAlphaPiAccept(int *c, int *l, double* acc)
{
	int C = (*c) - 1;
	int L = (*l) - 1;

	*acc = alpha_pi_acc[C][L];
}

void c212BB_poisson_mc_hier3_lev0::getBetaPiAccept(int *c, int* l,  double* acc)
{
	int C = (*c) - 1;
	int L = (*l) - 1;

	*acc = beta_pi_acc[C][L];
}
