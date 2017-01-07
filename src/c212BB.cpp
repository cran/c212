#include<cstdio>
#include<cstdlib>

#include<map>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "c2121a.h"
#include "c212BB.h"

static const char *rcsId = "$Id: c212BB.cpp,v 1.18 2017/01/06 10:26:20 clb13102 Exp clb13102 $";

const char* c212BB::sColPMweight = "weight_pm";

const char* c212BB::sParam_sigma_MH_gamma = "sigma_MH_gamma";
const char* c212BB::sParam_sigma_MH_theta = "sigma_MH_theta";
const char* c212BB::sParam_sigma_MH_alpha = "sigma_MH_alpha";
const char* c212BB::sParam_sigma_MH_beta = "sigma_MH_beta";
const char* c212BB::sParam_w_gamma = "w_gamma";
const char* c212BB::sParam_w_theta = "w_theta";
const char* c212BB::sParam_w_alpha = "w_alpha";
const char* c212BB::sParam_w_beta = "w_beta";


c212BB::c212BB(SEXP sChains, SEXP sBurnin, SEXP sIter,  SEXP sNumBodySys, SEXP sMaxAEs,
					SEXP pNAE,
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
                    SEXP psigma2_theta,
					SEXP palpha_pi,
					SEXP pbeta_pi,
					SEXP plambda_alpha,
					SEXP plambda_beta,
					SEXP pPi,
					SEXP palgo,
					SEXP padapt_phase,
					SEXP sim_type,
					SEXP global_sim_params,
					SEXP sim_params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP adapt_min_w,
					SEXP adapt_chains,
					SEXP adapt_burnin,
					SEXP adapt_iter,
					bool on_screen)
{
	gScreen = on_screen;

	// Give these default values to allow the Gibbs simulator to run even if no or
	// partial values are supplied.
	gSigma_MH_alpha = 3.0;
	gSigma_MH_beta = 3.0;
	gDefault_Sigma_MH_gamma = 0.2;
	gDefault_Sigma_MH_theta = 0.2;
	gW_alpha = 1.0;
	gW_beta = 1.0;
	gW_alpha_control = 6.0;
	gW_beta_control = 6.0;
	gDefault_W_gamma = 1.0;
	gDefault_W_gamma_control = 1.0;

	gWp = NULL;
	gMH_weight = 0.5;
	gMHAdaptParams.w_min = 0.25;
	gMHAdaptParams.chains = 1;
	gMHAdaptParams.burnin = 1000;
	gMHAdaptParams.iter = 10000;

	gTheta_zero_prop = NULL;
	gTheta_zero_acc = NULL;

	alpha_pi  = NULL;  // Current value of the sampled distribution
	beta_pi  = NULL;  // Current value of the sampled distribution
	lambda_alpha = 0.0;
	lambda_beta = 0.0;
	gPi = NULL;
	alpha_pi_acc  = NULL;  // Acceptance rate
	beta_pi_acc  = NULL;  // Acceptance rate

	alpha_pi_samples = NULL;
	beta_pi_samples = NULL;
	gPi_samples = NULL;

	gAlgo = (eAlgoType)(*(INTEGER(palgo)));
	gSimType = eSim_Type_SLICE;

	// Adapt MH phase
	gAdapt_Phase_alpha = 0;
	gAdapt_Phase_beta = 0;
	alpha_pi_acc_adapt  = 0;
	beta_pi_acc_adapt  = 0;

	gAdapt_Phase_theta = NULL;
	theta_acc_adapt = NULL;
	theta_mix_p = NULL;
	theta_max_p = NULL;

	// Adaptive MCMC
	gM = 1;
	gW0 = NULL;
	gW = NULL;
	gMU = NULL;
	gSIGMA2 = NULL;
	//gMU_tilde = NULL;
	//gSIGMA2_tilde = NULL;
	//gLAMBDA = NULL;

	// Adapt for some of the MH steps
	in_apapt_phase = 0;

	init(sChains, sBurnin, sIter, sNumBodySys, sMaxAEs, pNAE, pX, pY, pNC, pNT,
		ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0, pmu_theta_0_0, ptau2_theta_0_0,
        palpha_gamma_0_0, pbeta_gamma_0_0, palpha_theta_0_0, pbeta_theta_0_0,
		palpha_gamma, pbeta_gamma, palpha_theta, pbeta_theta, pmu_gamma_0,
		ptau2_gamma_0, pmu_theta_0, ptau2_theta_0, pmu_gamma, pmu_theta, psigma2_gamma,
		psigma2_theta, palpha_pi, pbeta_pi, plambda_alpha, plambda_beta, pPi,
		padapt_phase, sim_type, global_sim_params, sim_params, MH_weight, pm_weights,
		adapt_min_w, adapt_chains, adapt_burnin,
		adapt_iter);

	if (gScreen) {
		// Report the global parameters
		Rprintf("Global Simulation Parameters:\n");
		Rprintf("\tSimulation Type: %d\n", gAlgo);
		Rprintf("\tw_alpha (width): %0.6f\n", gW_alpha);
		Rprintf("\tm alpha (control): %0.6f\n", gW_alpha_control);
		Rprintf("\tw_beta (width): %0.6f\n", gW_beta);
		Rprintf("\tm beta (control): %0.6f\n", gW_beta_control);
		Rprintf("\tw_gamma (width): %0.6f\n", gDefault_W_gamma);
		Rprintf("\tm gamma (control): %0.6f\n", gDefault_W_gamma_control);
		Rprintf("\tsigma_MH_alpha: %0.6f\n", gSigma_MH_alpha);
		Rprintf("\tsigma_MH_beta: %0.6f\n", gSigma_MH_beta);
		Rprintf("\tsigma_MH_gamma: %0.6f\n", gDefault_Sigma_MH_gamma);
		Rprintf("\tsigma_MH_theta: %0.6f\n", gDefault_Sigma_MH_theta);
		Rprintf("\tdefault weight: %0.6f\n", gMH_weight);
	}
}

c212BB::c212BB(int sChains, int sBurnin, int sIter,  int sNumBodySys, int sMaxAEs,
					int* pNAE,
					int** pX, int** pY,
					int** pNC, int** pNT,
					double*** ptheta, double*** pgamma,
                    double pmu_gamma_0_0,
                    double ptau2_gamma_0_0,
                    double pmu_theta_0_0,
                    double ptau2_theta_0_0,
                    double palpha_gamma_0_0,
                    double pbeta_gamma_0_0,
                    double palpha_theta_0_0,
                    double pbeta_theta_0_0,
                    double palpha_gamma,
                    double pbeta_gamma,
                    double palpha_theta,
                    double pbeta_theta,
                    double* pmu_gamma_0,
					double* ptau2_gamma_0,
                    double* pmu_theta_0,
                    double* ptau2_theta_0,
                    double** pmu_gamma,
                    double** pmu_theta,
                    double** psigma2_gamma,
                    double** psigma2_theta,
					double* palpha_pi,
					double* pbeta_pi,
					double plambda_alpha,
					double plambda_beta,
					double** pPi,
					int palgo,
					int padapt_phase,
					const char* sim_type,
					std::map<const char*, gSimParams>& mGlobalSimParams,
					double** gamma_params,
					int** gamma_ctrl,
					double** theta_params,
					double MH_weight,
					double** pm_weights,
					MHAdaptParams& adapt_params, bool on_screen)
{
	gScreen = on_screen;

	gBurnin = 0;
	gChains = 0;
	gIter = 0;
	gNAE = NULL;
	gNumBodySys = 0;
	gMaxAEs = 0;

	gSim_Param = 0.0;
	gSim_Param_cntrl = 0.0;
	gW_gamma = NULL;
	gW_theta = NULL;
	gW_gamma_control = NULL;
	gW_theta_control = NULL;
	gSigma_MH_gamma = NULL;
	gSigma_MH_theta = NULL;

	// Give these default values to allow the Gibbs simulator to run even if no or
	// partial values are supplied.
	gSigma_MH_alpha = 3.0;
	gSigma_MH_beta = 3.0;
	gDefault_Sigma_MH_gamma = 0.2;
	gDefault_Sigma_MH_theta = 0.2;
	gW_alpha = 1.0;
	gW_beta = 1.0;
	gW_alpha_control = 6.0;
	gW_beta = 1.0;
	gW_beta_control = 6.0;
	gDefault_W_gamma = 1.0;
	gDefault_W_gamma_control = 6.0;

	gWp = NULL;
	gMH_weight = 0.5;
	gMHAdaptParams.w_min = 0.25;
	gMHAdaptParams.chains = 1;
	gMHAdaptParams.burnin = 1000;
	gMHAdaptParams.iter = 10000;

	mu_theta_0_0 = 0.0;   // Fixed hyper-parameter value
	mu_gamma_0_0  = 0.0; // Fixed hyper-parameter value
	tau2_theta_0_0 = 0.0; // Fixed hyper-parameter value
	tau2_gamma_0_0 = 0.0; // Fixed hyper-parameter value
	alpha_gamma_0_0 = 0.0; // Fixed hype-parameter value
	beta_gamma_0_0 = 0.0; // Fixed hype-parameter value
	alpha_theta_0_0 = 0.0; // Fixed hype-parameter value
	beta_theta_0_0 = 0.0; // Fixed hype-parameter value
	alpha_gamma = 0.0; // Fixed hype-parameter value
	beta_gamma = 0.0; // Fixed hype-parameter value
	alpha_theta = 0.0; // Fixed hype-parameter value
	beta_theta = 0.0; // Fixed hype-parameter value

	mu_theta_0 = NULL; // Current value of the sampled distribution
	mu_gamma_0 = NULL; // Current value of the sampled distribution
	tau2_theta_0 = NULL; // Current value of the sampled distribution
	tau2_gamma_0 = NULL; // Current value of the sampled distribution

	gTheta_zero_prop = NULL;
	gTheta_zero_acc = NULL;

	alpha_pi  = NULL;  // Current value of the sampled distribution
	beta_pi  = NULL;  // Current value of the sampled distribution
	lambda_alpha = 0.0;
	lambda_beta = 0.0;
	gPi = NULL;
	alpha_pi_acc  = NULL;  // Acceptance rate
	beta_pi_acc  = NULL;  // Acceptance rate

	mu_theta = NULL;
	mu_gamma = NULL;
	sigma2_theta = NULL;
	sigma2_gamma = NULL;

	gTheta = NULL;
	gGamma = NULL;
	gTheta_acc = NULL;
	gGamma_acc = NULL;

	// Data values
	x = NULL;
	y = NULL;
	NC = NULL;
	NT = NULL;

	// Samples
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

	alpha_pi_samples = NULL;
	beta_pi_samples = NULL;
	gPi_samples = NULL;

	gSimType = eSim_Type_SLICE;
	gAlgo = (eAlgoType)(palgo);

	// Adapt MH phase
	gAdapt_Phase_alpha = 0;
	gAdapt_Phase_beta = 0;
	alpha_pi_acc_adapt  = 0;
	beta_pi_acc_adapt  = 0;

	gAdapt_Phase_theta = NULL;
	theta_acc_adapt = NULL;
	theta_mix_p = NULL;
	theta_max_p = NULL;

	// Adaptive MCMC
	gM = 1;
	gW0 = NULL;
	gW = NULL;
	gMU = NULL;
	gSIGMA2 = NULL;
	//gMU_tilde = NULL;
	//gSIGMA2_tilde = NULL;
	//gLAMBDA = NULL;

	// Adapt for some of the MH steps
	in_apapt_phase = 0;

	init(sChains, sBurnin, sIter, sNumBodySys, sMaxAEs, pNAE, pX, pY, pNC, pNT,
		ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0, pmu_theta_0_0, ptau2_theta_0_0,
        palpha_gamma_0_0, pbeta_gamma_0_0, palpha_theta_0_0, pbeta_theta_0_0,
		palpha_gamma, pbeta_gamma, palpha_theta, pbeta_theta, pmu_gamma_0,
		ptau2_gamma_0, pmu_theta_0, ptau2_theta_0, pmu_gamma, pmu_theta, psigma2_gamma,
		psigma2_theta, palpha_pi, pbeta_pi, plambda_alpha, plambda_beta, pPi,
		padapt_phase, sim_type, mGlobalSimParams,
		gamma_params, gamma_ctrl, theta_params,
		MH_weight,
		pm_weights,
		adapt_params);
}

c212BB::~c212BB()
{
	//Rprintf("c212BB::c212BB - destructor\n");
	release();
}

void c212BB::gibbs_sampler()
{
	if (gAlgo == MH_ADAPT) {
		adaptPhaseMH();
	}

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
}

void c212BB::initGlobalSimParams(SEXP sim_type, SEXP sim_params)
{
	int len = Rf_length(sim_params);

	//SEXP sType = R_NilValue;
	SEXP sParams = R_NilValue;
	SEXP sValues = R_NilValue;
	SEXP sControl = R_NilValue;

	if (strcmp("MH", CHAR(STRING_ELT(sim_type, 0))) == 0) {
		gSimType = eSim_Type_MH;
	}
	else {
		gSimType = eSim_Type_SLICE;
	}

	if (isNewList(sim_params)) {

		SEXP names = getAttrib(sim_params, R_NamesSymbol);

		int i = 0;

		for (i = 0; i < len; i++) {
			if (strcmp(sColValue, CHAR(STRING_ELT(names, i))) == 0) {
				sValues = VECTOR_ELT(sim_params, i);
			}
			//if (strcmp(sColType, CHAR(STRING_ELT(names, i))) == 0) {
			//	sType = VECTOR_ELT(sim_params, i);
			//}
			if (strcmp(sColParam, CHAR(STRING_ELT(names, i))) == 0) {
				sParams = VECTOR_ELT(sim_params, i);
			}
			if (strcmp(sColControl, CHAR(STRING_ELT(names, i))) == 0) {
				sControl = VECTOR_ELT(sim_params, i);
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
					mGlobalSimParams[sParam_sigma_MH_gamma]  = gSimParams(vals[i], 0.0);
				}
				if (0 == strcmp(t, sParam_sigma_MH_theta)) {
					gDefault_Sigma_MH_theta = vals[i];
					mGlobalSimParams[sParam_sigma_MH_theta]  = gSimParams(vals[i], 0.0);
				}
				if (0 == strcmp(t, sParam_sigma_MH_alpha)) {
					gSigma_MH_alpha = vals[i];
					mGlobalSimParams[sParam_sigma_MH_alpha]  = gSimParams(vals[i], 0.0);
				}
				if (0 == strcmp(t, sParam_sigma_MH_beta)) {
					gSigma_MH_beta = vals[i];
					mGlobalSimParams[sParam_sigma_MH_beta]  = gSimParams(vals[i], 0.0);
				}
				if (0 == strcmp(t, sParam_w_gamma)) {
					gDefault_W_gamma = vals[i];
					gDefault_W_gamma_control = cntrl[i];
					mGlobalSimParams[sParam_w_gamma]  = gSimParams(vals[i], cntrl[i]);
				}
				if (0 == strcmp(t, sParam_w_alpha)) {
					gW_alpha = vals[i];
					gW_alpha_control = cntrl[i];
					mGlobalSimParams[sParam_w_alpha]  = gSimParams(vals[i], cntrl[i]);
				}
				if (0 == strcmp(t, sParam_w_beta)) {
					gW_beta = vals[i];
					gW_beta_control = cntrl[i];
					mGlobalSimParams[sParam_w_beta]  = gSimParams(vals[i], cntrl[i]);
				}
			}
		}
	}
}

void c212BB::initGlobalSimParams(const char* sim_type, std::map<const char*, gSimParams>& mGlobalSimParams)
{
	gDefault_Sigma_MH_gamma = mGlobalSimParams[sParam_sigma_MH_gamma].value;
	gDefault_Sigma_MH_theta = mGlobalSimParams[sParam_sigma_MH_theta].value;
	gSigma_MH_alpha = mGlobalSimParams[sParam_sigma_MH_alpha].value;
	gSigma_MH_beta = mGlobalSimParams[sParam_sigma_MH_beta].value;

	gDefault_W_gamma = mGlobalSimParams[sParam_w_gamma].value;
	gDefault_W_gamma_control = mGlobalSimParams[sParam_w_gamma].control;

	gW_alpha = mGlobalSimParams[sParam_w_alpha].value;
	gW_alpha_control = mGlobalSimParams[sParam_w_alpha].control;

	gW_beta = mGlobalSimParams[sParam_w_beta].value;
	gW_beta_control = mGlobalSimParams[sParam_w_beta].control;

	if (strcmp("MH", sim_type) == 0) {
		gSimType = eSim_Type_MH;
	}
	else {
		gSimType = eSim_Type_SLICE;
	}
}

void c212BB::initSimParams(SEXP sim_params)
{
	gW_gamma = (double **)malloc(gNumBodySys * sizeof(double*));
	gW_gamma_control = (int **)malloc(gNumBodySys * sizeof(int*));
	gSigma_MH_gamma = (double **)malloc(gNumBodySys * sizeof(double*));
	gSigma_MH_theta = (double **)malloc(gNumBodySys * sizeof(double*));

	int i = 0, j = 0;
	for (i = 0; i < gNumBodySys; i++) {
		gW_gamma[i] = (double *)malloc(gNAE[i] * sizeof(double));
		gW_gamma_control[i] = (int *)malloc(gNAE[i] * sizeof(int));
		gSigma_MH_gamma[i] = (double *)malloc(gNAE[i] * sizeof(double));
		gSigma_MH_theta[i] = (double *)malloc(gNAE[i] * sizeof(double));
		for (j = 0; j < gNAE[i]; j++) {
			gW_gamma[i][j] = gDefault_W_gamma;
			gW_gamma_control[i][j] = (int)gDefault_W_gamma_control;
			gSigma_MH_gamma[i][j] = gDefault_Sigma_MH_gamma;
			gSigma_MH_theta[i][j] = gDefault_Sigma_MH_theta;
		}
	}

	int len = Rf_length(sim_params);

	SEXP sVariables = R_NilValue;
	SEXP sParams = R_NilValue;
	SEXP sValues = R_NilValue;
	SEXP sControl = R_NilValue;
	SEXP sB = R_NilValue;
	SEXP sj = R_NilValue;

	if (len && isNewList(sim_params)) {

		SEXP names = getAttrib(sim_params, R_NamesSymbol);

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
					if (0 == strcmp(param, sParam_w_gamma)) {
						gW_gamma[b][a] = vals[i];
						gW_gamma_control[b][a] = (int)cntrl[i];
					}
					else if (0 == strcmp(param, sParam_sigma_MH_gamma)) {
						gSigma_MH_gamma[b][a] = vals[i];
					}
				}
				else if (0 == strcmp(sVariable_theta, var)) {
					if (0 == strcmp(param, sParam_sigma_MH_theta)) {
						gSigma_MH_theta[b][a] = vals[i];
					}
				}
			}
		}
	}
}

void c212BB::initSimParams(double** gamma_params, int** gamma_cntrl, double** theta_params)
{
	gW_gamma = (double **)malloc(gNumBodySys * sizeof(double*));
	gW_gamma_control = (int **)malloc(gNumBodySys * sizeof(int*));
	gSigma_MH_gamma = (double **)malloc(gNumBodySys * sizeof(double*));
	gSigma_MH_theta = (double **)malloc(gNumBodySys * sizeof(double*));

	int i = 0, j = 0;
	for (i = 0; i < gNumBodySys; i++) {
		gW_gamma[i] = (double *)malloc(gNAE[i] * sizeof(double));
		gW_gamma_control[i] = (int *)malloc(gNAE[i] * sizeof(int));
		gSigma_MH_gamma[i] = (double *)malloc(gNAE[i] * sizeof(double));
		gSigma_MH_theta[i] = (double *)malloc(gNAE[i] * sizeof(double));
		for (j = 0; j < gNAE[i]; j++) {
			if (gSimType == eSim_Type_SLICE) {
				gW_gamma[i][j] = gamma_params[i][j];
				gW_gamma_control[i][j] = (int)gamma_cntrl[i][j];
			}
			else {
				gSigma_MH_gamma[i][j] = gamma_params[i][j];
			}
			gSigma_MH_theta[i][j] = theta_params[i][j];
		}
	}
}

void c212BB::initPMWeights(SEXP pm_weights)
{
	gWp = (double**)malloc(gNumBodySys *sizeof(double*));

	int i = 0, j = 0;
	for (i = 0; i < gNumBodySys; i++) {
		gWp[i] = (double *)malloc(gNAE[i] * sizeof(double));
		for (j = 0; j < gNAE[i]; j++) {
			gWp[i][j] = gMH_weight;
		}
	}

	int len = Rf_length(pm_weights);

	SEXP sPM_Weights = R_NilValue;
	SEXP sB = R_NilValue;
	SEXP sj = R_NilValue;

	if (len && isNewList(pm_weights)) {

		SEXP names = getAttrib(pm_weights, R_NamesSymbol);

		for (i = 0; i < len; i++) {
			if (strcmp(sColPMweight, CHAR(STRING_ELT(names, i))) == 0) {
				sPM_Weights = VECTOR_ELT(pm_weights, i);
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
			int* B = INTEGER(sB);
			int* j = INTEGER(sj);

			for (i = 0; i < len; i++) {
				int b = B[i] - 1;
				int a = j[i] - 1;
				gWp[b][a] = weights[i];
			}
		}
	}
}

void c212BB::initPMWeights(double** pm_weights)
{
	int i = 0, j = 0;

	gWp = (double**)malloc(gNumBodySys *sizeof(double*));

	for (i = 0; i < gNumBodySys; i++) {
		gWp[i] = (double *)malloc(gNAE[i] * sizeof(double));
		for (j = 0; j < gNAE[i]; j++) {
			gWp[i][j] = pm_weights[i][j];
		}
	}
}

void c212BB::adaptPhaseMH()
{
	if (gScreen)
		Rprintf("Adaptive phase...\n");

	const char *sAdapt_Sim_Type = "SLICE";

	c212BB adapt(gMHAdaptParams.chains, gMHAdaptParams.burnin, gMHAdaptParams.iter,
					gNumBodySys, gMaxAEs,
					gNAE,
					x, y,
					NC, NT,
					gTheta, gGamma,
                    mu_gamma_0_0,
                    tau2_gamma_0_0,
                    mu_theta_0_0,
                    tau2_theta_0_0,
                    alpha_gamma_0_0,
                    beta_gamma_0_0,
                    alpha_theta_0_0,
                    beta_theta_0_0,
                    alpha_gamma,
                    beta_gamma,
                    alpha_theta,
                    beta_theta,
                    mu_gamma_0,
					tau2_gamma_0,
                    mu_theta_0,
                    tau2_theta_0,
                    mu_gamma,
                    mu_theta,
                    sigma2_gamma,
                    sigma2_theta,
					alpha_pi,
					beta_pi,
					lambda_alpha,
					lambda_beta,
					gPi,
					MH,
					gAdapt_Phase_alpha,
					sAdapt_Sim_Type,
					mGlobalSimParams,
					gW_gamma,
					gW_gamma_control,
					gSigma_MH_theta,
					gMH_weight,
					gWp,
					gMHAdaptParams,
					FALSE);

	adapt.gibbs_sampler();

	int c = 0, i = 0, b = 0, j = 0;
	int sz = gMHAdaptParams.iter - gMHAdaptParams.burnin;
	for (b = 0; b < gNumBodySys; b++) {
		for (j = 0; j < gNAE[b]; j++) {
			int z_zero = 0;
			for (c = 0; c < gMHAdaptParams.chains; c++) {
				double *theta_samples = adapt.getThetaSamples(c, b, j);
				for (i = 0; i < sz; i++) {
					if (theta_samples[i] == 0.0)
						z_zero++;
				}
			}

			double w = ((double)z_zero)/((double)(c * sz));
			gWp[b][j] = w;

			if (gWp[b][j] < gMHAdaptParams.w_min)
				gWp[b][j] = gMHAdaptParams.w_min;
			if (gWp[b][j] > 1 - gMHAdaptParams.w_min)
				gWp[b][j] = 1 - gMHAdaptParams.w_min;
		}
	}

	if (gScreen)
		Rprintf("Complete.\n");
}

void c212BB::simulate_MH()
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

			sample_alpha_pi_MH(c, gBurnin, i);
			sample_beta_pi_MH(c, gBurnin, i);
			sample_pi(c, gBurnin, i);

			sample_mu_gamma(c, gBurnin, i);
			sample_mu_theta(c, gBurnin, i);
			sample_sigma2_gamma(c, gBurnin, i);
			sample_sigma2_theta(c, gBurnin, i);
			sample_gamma_MH(c, gBurnin, i);

			switch(gAlgo) {
				//case BB2004:
				//	sample_theta_BB2004(c, gBurnin, i);
				//break;

				case MH:
				case MH_ADAPT:
					sample_theta_MH(c, gBurnin, i);
				break;

				case MIS_ADAPT:
					sample_theta_MIS_Adapt(c, gBurnin, i);
				break;

				case INDEP:
					sample_theta_Independent_MH(c, gBurnin, i);
				break;

				default:
					sample_theta_MH(c, gBurnin, i);
				break;
			}
#ifndef INDIVIDUAL_RNG
			PutRNGstate();
#endif
			if (((i + 1)%1000 == 0) && gScreen) {
				Rprintf("%d iterations...\n", i + 1);
			}
		}
	}
	if (gScreen)
		Rprintf("MCMC chain fitting complete.\n");
}

void c212BB::simulate_SLICE()
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

			sample_alpha_pi_SLICE(c, gBurnin, i);
			sample_beta_pi_SLICE(c, gBurnin, i);
			sample_pi(c, gBurnin, i);

			sample_mu_gamma(c, gBurnin, i);
			sample_mu_theta(c, gBurnin, i);
			sample_sigma2_gamma(c, gBurnin, i);
			sample_sigma2_theta(c, gBurnin, i);
			sample_gamma_SLICE(c, gBurnin, i);

			switch(gAlgo) {
				//case BB2004:
				//	sample_theta_BB2004(c, gBurnin, i);
				//break;

				case MH:
				case MH_ADAPT:
					sample_theta_MH(c, gBurnin, i);
				break;

				case MIS_ADAPT:
					sample_theta_MIS_Adapt(c, gBurnin, i);
				break;

				case INDEP:
					sample_theta_Independent_MH(c, gBurnin, i);
				break;

				default:
					sample_theta_MH(c, gBurnin, i);
				break;
			}
#ifndef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			if (((i + 1)%1000 == 0) && gScreen) {
				Rprintf("%d iterations...\n", i + 1);
			}
		}
	}
	if (gScreen)
		Rprintf("MCMC chain fitting complete.\n");
}

double c212BB::log_f_alpha_pi(int c, double alpha)
{
	double f = 0.0;
	double log_pi_sum = 0.0;

    int b = 0;
	for (b = 0; b < gNumBodySys; b++) {
		 log_pi_sum += log(gPi[c][b]);
	}

	f = ((double)gNumBodySys) * (lgamma(alpha + beta_pi[c]) - lgamma(alpha)) + (alpha - 1.0) * log_pi_sum - alpha * lambda_alpha;

    return(f);
}

void c212BB::sample_alpha_pi_MH(int c, int burnin, int iter)
{
	double cand = 0;

	// Original implementation was incorrect. We need to include an extra term in
	// the ratio corresponding to the fact that this is a rejection sampling of 
	// rnorm trucated at 1.
	// alpha_pi is restricted to being greater than one

	while (cand <= 1.0) {
#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		cand = rnorm(alpha_pi[c], gSigma_MH_alpha);
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

	// PHI((mu - 1)/sigma)
	double q1 = pnorm((alpha_pi[c] - 1)/gSigma_MH_alpha, 0, 1, 1, 0);
	double q2 = pnorm((cand - 1)/gSigma_MH_alpha, 0, 1, 1, 0);

	double ratio = (exp(f1 - f2)) * q1/q2;

	ratio = cMIN(ratio, 1);

	if (u <= ratio) {
		alpha_pi[c] = cand;
		//if (iter >= burnin)
			alpha_pi_acc[c] = alpha_pi_acc[c] + 1;

        alpha_pi_acc_adapt = alpha_pi_acc_adapt + 1;
	}

	// Attempt to adjust the "step" size to get an acceptable acceptance rate
//	if (gAdapt_Phase_alpha && (iter < burnin) && iter && (iter % 2000 == 0)) {
//		// If the acceptance rate is too small we need to decrease the "step"
//		if (alpha_pi_acc_adapt/2000.0 < 0.25) {
//			if (gSigma_MH_alpha > 0.25) {
//				gSigma_MH_alpha -= 0.25;
//			}
//			else {
//				gSigma_MH_alpha = gSigma_MH_alpha/2.0;
//			}
//		}
//		else if (alpha_pi_acc_adapt/2000.0 > 0.55) {
//			// If the acceptance rate is too large we need to increase the "step"
//			gSigma_MH_alpha += 0.25;
//		}
//		else {
//			// Acceptable rate
//			gAdapt_Phase_alpha = 0;
//		}
 //       alpha_pi_acc_adapt = 0;
//	}

	if (iter >= burnin) {
		alpha_pi_samples[c][iter - burnin] = alpha_pi[c];
	}
}

void c212BB::sample_alpha_pi_SLICE(int c, int burnin, int iter)
{
	int m = (int)gW_alpha_control, K = 0, J = 0;

#ifdef INDIVIDUAL_RNG
	GetRNGstate();
#endif
	J = floor(runif(0, m));
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
	double u = runif(0, gW_alpha);
#ifdef INDIVIDUAL_RNG
	PutRNGstate();
#endif

	l = alpha_pi[c] - u;
	r = alpha_pi[c] + (gW_alpha - u);

	while (J > 0) {
		if (l <= 1.0)
			break;

		if (logy >= log_f_alpha_pi(c, l)) {
			break;
		}
		l = l - gW_alpha;

		J--;
	}

	while (K > 0) {
		if (logy >= log_f_alpha_pi(c, r)) {
			break;
		}
		r = r + gW_alpha;
		K--;
	}

	// 1.0 is the lowerbound for alpha_pi, i.e. alpha_pi > 1. The probability of
	// sampling 1.0 is zero.
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

	if (iter >= burnin) {
		alpha_pi_samples[c][iter - burnin] = alpha_pi[c];
	}
}

double c212BB::log_f_beta_pi(int c, double beta)
{
	double f = 0.0;
    double log_sum = 0.0;

    int b = 0;
	for (b = 0; b < gNumBodySys; b++) {
    	log_sum += log(1.0 - gPi[c][b]);
	}

    f = ((double)gNumBodySys) * (lgamma(alpha_pi[c] + beta) - lgamma(beta)) + (beta - 1.0) * log_sum - beta * lambda_beta;

    return(f);
}

void c212BB::sample_beta_pi_MH(int c, int burnin, int iter)
{
    double cand = 0.0;

    // beta_pi is restricted to being greater than zero
    while (cand <= 1.0) {
#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
        cand = rnorm(beta_pi[c], gSigma_MH_beta);
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

	double q1 = pnorm((beta_pi[c] - 1)/gSigma_MH_beta, 0, 1, 1, 0);
	double q2 = pnorm((cand - 1)/gSigma_MH_beta, 0, 1, 1, 0);

    double ratio = (exp(f1 - f2)) * (q1/q2);

    ratio = cMIN(ratio, 1);

    if (u <= ratio) {
        beta_pi[c] = cand;
		//if (iter >= burnin)
        	beta_pi_acc[c] = beta_pi_acc[c] + 1;

        beta_pi_acc_adapt = beta_pi_acc_adapt + 1;
    }

	// Attempt to adjust the "step" size to get an acceptable acceptance rate
//	if (gAdapt_Phase_beta && (iter < burnin) && iter && (iter % 2000 == 0)) {
//		// If the acceptance rate is too small we need to decrease the "step"
//		if (beta_pi_acc_adapt/2000.0 < 0.25) {
//			if (gSigma_MH_beta > 0.5) {
//				gSigma_MH_beta -= 0.5;
//			}
//			else {
//				gSigma_MH_beta = gSigma_MH_beta/2.0;
//			}
//		}
//		else if (beta_pi_acc_adapt/2000.0 > 0.55) {
//			// If the acceptance rate is too large we need to increase the "step"
//			gSigma_MH_beta += 0.5;
//		}
//		else {
//			gAdapt_Phase_beta = 0;
//		}
 //       beta_pi_acc_adapt = 0;
//	}

    if (iter >= burnin) {
        beta_pi_samples[c][iter - burnin] = beta_pi[c];
    }
}

void c212BB::sample_beta_pi_SLICE(int c, int burnin, int iter)
{
	int m = gW_beta_control, K = 0, J = 0;

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
	double u = runif(0, gW_beta);
#ifdef INDIVIDUAL_RNG
	PutRNGstate();
#endif

	l = beta_pi[c] - u;
	r = beta_pi[c] + (gW_beta - u);

	// beta is retricted to being greater than 1
	while (J > 0) {
		if (l <= 1.0)
			break;

		if (logy >= log_f_beta_pi(c, l)) {
			break;
		}
		l = l - gW_beta;
		J--;
	}

	while (K > 0) {
		if (logy >= log_f_beta_pi(c, r)) {
			break;
		}
		r = r + gW_beta;
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

    if (iter >= burnin) {
        beta_pi_samples[c][iter - burnin] = beta_pi[c];
    }
}

void c212BB::sample_pi(int c, int burnin, int iter)
{
	int b = 0;
    for (b = 0; b < gNumBodySys; b++) {
		int theta_zero_count = 0;

		int j = 0;
		for (j = 0; j< gNAE[b]; j++) {
			if (gTheta[c][b][j] == 0.0) {
        		theta_zero_count++;
			}
		}

        double shape1 = alpha_pi[c] + (double)theta_zero_count;
        double shape2 = beta_pi[c] + (double)gNAE[b] - (double)theta_zero_count;

#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
        gPi[c][b] = rbeta(shape1, shape2);
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif

        if (iter >= burnin) {
            gPi_samples[c][b][iter - burnin] = gPi[c][b];
        }
    }
}

void c212BB::sample_mu_theta(int c, int burnin, int iter)
{
	int b = 0;

	for (b = 0; b < gNumBodySys; b++) {

		double t = 0.0;
		int K_b = 0;
		int j = 0;
		for (j = 0; j < gNAE[b]; j++) {
			if (gTheta[c][b][j] != 0.0) {
				K_b++;
			}
			t += gTheta[c][b][j];
		}

		double denom = sigma2_theta[c][b] + ((double)K_b)*tau2_theta_0[c];

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

void c212BB::sample_sigma2_theta(int c, int burnin, int iter)
{

	int b = 0;

	for (b = 0; b < gNumBodySys; b++) {

		double t = 0;
		int K_b = 0;
		int j = 0;
		for (j = 0; j < gNAE[b]; j++) {
			if (gTheta[c][b][j] != 0.0) {
				K_b++;
				t += (pow((gTheta[c][b][j] - mu_theta[c][b]),2.0));
			}
		}

		double s = alpha_theta + ((double)K_b)/2.0;

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

double c212BB::log_f_theta(int c, int b, int j, double theta)
{
	double f1 = ((double)y[b][j]) * (gGamma[c][b][j] + theta) - ((double)NT[b][j]) *(log(1 + exp(gGamma[c][b][j] + theta)));

	double f2 = 0;

	if (theta == 0) {
		f2 = log(gPi[c][b]);
	}

	if (theta != 0) {
		f2 = log(1 - gPi[c][b]) + log(1.0/sqrt(2.0*M_PI*sigma2_theta[c][b])) + ((-1.0/2.0)*(pow(theta -mu_theta[c][b], 2.0))/sigma2_theta[c][b]);
	}

	double f = f1 + f2;

	return(f);
}

// Implementation as described in BB2004.
// There were a number of issues:
// 1. There is a typo in the appendix of the original paper:
//		For theta^C <> 0, theta^0 = 0 the denominator should use theta^C in the
//		exponential and not theta^0
// 2. Using r(theta^C, theta^0) etc can cause the ratio to become infinite. So
//		we replace them with log until the comparison.
// Despite what the appendix says this is just a standard MH step with mixture proposal 
// distribution, actually the same as we implemented in sample_theta_MH. So this
// function is no longer needed.
//void c212BB::sample_theta_BB2004(int c, int burnin, int iter)
//{
//#ifdef __C212DEBUG
//	Rprintf("11a - Samples:\n");
//#endif
//
//	int b = 0, j = 0;
//	for (b = 0; b < gNumBodySys; b++) {
//		for ( j = 0; j < gNAE[b]; j++) {
//
//			// Following BB2004 Appendix
//			// lambda = 0.5 and w_{m,n} = 0 giving:
//			// q(y|x) ~ 0.5 I(theta = 0) + (0.5) N(x, M$sigma_MH_theta^2)
//			// This is NOT an adaptive MCMC algorithm.
//
//#ifdef INDIVIDUAL_RNG
//			GetRNGstate();
//#endif
//			double u = runif(0, 1);
//#ifdef INDIVIDUAL_RNG
//			PutRNGstate();
//#endif
//#ifdef __C212DEBUG
//			Rprintf("\t\tu: %0.6f\n", u);
//#endif
//
//			double cand = 0.0;
//
//			if (u < 0.5) {
//				cand = 0.0;
//			}
//			else {
//#ifdef INDIVIDUAL_RNG
//				GetRNGstate();
//#endif
//				cand = rnorm(gTheta[c][b][j], gSigma_MH_theta[b][j]);
//#ifdef INDIVIDUAL_RNG
//				PutRNGstate();
//#endif
//			}
//#ifdef __C212DEBUG
//			Rprintf("\t\tcand: %0.6f\n", cand);
//#endif
//
//
//			double f1 = log_f_theta(c, b, j, cand);
//			double f2 = log_f_theta(c, b, j, gTheta[c][b][j]);
//
//			double ratio = exp(f1 - f2);
//			double lratio = f1 - f2;
//#ifdef __C212DEBUG
//			Rprintf("\t\tf1 f2: %0.6f %0.6f\n", f1, f2);
//#endif
//
//#ifdef __C212DEBUG
//			Rprintf("\t\tratio: %0.6f\n", ratio);
//#endif
//
//			//ratio = cMIN(ratio, 1);
//
//#ifdef INDIVIDUAL_RNG
//			GetRNGstate();
//#endif
//			u = runif(0, 1);
//#ifdef INDIVIDUAL_RNG
//			PutRNGstate();
//#endif
//
//#ifdef __C212DEBUG
//			Rprintf("\t\tu: %0.6f\n", u);
//			Rprintf("\t\tcurr: %0.6f\n", gTheta[c][b][j]);
//#endif
//
//			double m_old = (1.0/(gSigma_MH_theta[b][j] * sqrt(2.0*M_PI))) * exp((-1.0/(2.0*gSigma_MH_theta[b][j]*gSigma_MH_theta[b][j])) * pow(gTheta[c][b][j], 2.0));
//			double m_cand = (1.0/(gSigma_MH_theta[b][j] * sqrt(2.0*M_PI))) * exp((-1.0/(2.0*gSigma_MH_theta[b][j]*gSigma_MH_theta[b][j])) * pow(cand, 2.0));
//
//			double l_old = log(1) - log((gSigma_MH_theta[b][j] * sqrt(2.0*M_PI))) + ((-1.0/(2.0*gSigma_MH_theta[b][j]*gSigma_MH_theta[b][j])) * pow(gTheta[c][b][j], 2.0));
//			double l_cand = log(1) - log((gSigma_MH_theta[b][j] * sqrt(2.0*M_PI))) + ((-1.0/(2.0*gSigma_MH_theta[b][j]*gSigma_MH_theta[b][j])) * pow(cand, 2.0));
//
//#ifdef __C212DEBUG
//			Rprintf("\t\tm_old m_cand: %0.6f %0.6f\n", m_old, m_cand);
//#endif
//
//			if (cand == 0.0 && gTheta[c][b][j] == 0.0) {
//				gTheta[c][b][j] = cand;
//				if (iter >= burnin)
//					gTheta_acc[c][b][j] = gTheta_acc[c][b][j] + 1;
//#ifdef __C212DEBUG
//				Rprintf("\t\ttheta - ACC1\n");
//#endif
//			}
//			else {
//
//				if (cand == 0 && gTheta[c][b][j] != 0.0) {
//#ifdef __C212DEBUG
//					double xxx = m_old*ratio;
//					Rprintf("\t\tACC2: theta - %0.6f %0.6f\n", u, xxx);
//#endif
//					//if (u <= m_old * ratio) {
//					if (u <= exp(l_old + lratio)) {
//						gTheta[c][b][j] = cand;
//						if (iter >= burnin)
//							gTheta_acc[c][b][j] = gTheta_acc[c][b][j] + 1;
//#ifdef __C212DEBUG
//						Rprintf("\t\ttheta - ACC2\n");
//#endif
//					}
//				}
//				else {
//
//					if (cand != 0.0 && gTheta[c][b][j] == 0.0) {
//						// This looks to be a mistake in the Appendix!!!
//						//if (u <= ratio / m_old) {
//						//if (u <= ratio / m_cand) {
//						if (u <= exp(lratio - l_cand)) {
//							gTheta[c][b][j] = cand;
//							if (iter >= burnin)
//								gTheta_acc[c][b][j] = gTheta_acc[c][b][j] + 1;
//#ifdef __C212DEBUG
//							Rprintf("\t\ttheta - ACC3\n");
//#endif
//						}
//					}
//					else {
//
//						if (cand != 0.0 && gTheta[c][b][j] != 0.0) {
//							if (u <= ratio) {
//								gTheta[c][b][j] = cand;
//								if (iter >= burnin)
//									gTheta_acc[c][b][j] = gTheta_acc[c][b][j] + 1;
//#ifdef __C212DEBUG
//								Rprintf("\t\ttheta - ACC4\n");
//#endif
//							}
//						}
//
//
//					}
//				}
//			}
//
//			if (iter >= burnin) {
//				gTheta_samples[c][b][j][iter - burnin] = gTheta[c][b][j];
//			}
//
//#ifdef __C212DEBUG
//			Rprintf("\tVal: %0.6f\n", gTheta[b][j]);
//#endif
//		}
//	}
//}

double c212BB::log_q_theta(int c, int b, int j, double p, double theta, double mean)
{
	double f = 0.0;

	if (theta == 0.0) {
		f = log(p);
	}
	else {
		f = log(1 - p) + log((1.0/(gSigma_MH_theta[b][j] * sqrt(2.0 * M_PI))))  + (-1.0/(2.0*gSigma_MH_theta[b][j]*gSigma_MH_theta[b][j])) * pow((theta - mean), 2.0);
	}

	return f;
}

/*
* Sample theta using a MH step as detailed in:
* Gottardo, Raftery - Markov Chain Monte Carlo
* With Mixtures of Mutually Singular Distributions
*/
void c212BB::sample_theta_MH(int c, int burnin, int iter)
{
	int b = 0, j = 0;
	for (b = 0; b < gNumBodySys; b++) {
		for ( j = 0; j < gNAE[b]; j++) {

			// q(y|x) ~ 0.5 I(theta = 0) + (0.5) N(x, M$sigma_MH_theta^2)

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			double cand = 0.0;

			if (u < gWp[b][j]) {
				cand = 0.0;
				gTheta_zero_prop[c][b][j]++;
			}
			else {
#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				cand = rnorm(gTheta[c][b][j], gSigma_MH_theta[b][j]);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif
			}

			double f_cand = log_f_theta(c, b, j, cand);
			double f_prev = log_f_theta(c, b, j, gTheta[c][b][j]);

			double q_cand = log_q_theta(c, b, j, gWp[b][j], cand, gTheta[c][b][j]);
			double q_prev = log_q_theta(c, b, j, gWp[b][j], gTheta[c][b][j], cand);

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
				gTheta[c][b][j] = cand;
				gTheta_acc[c][b][j] = gTheta_acc[c][b][j] + 1;
				if (cand == 0)
					gTheta_zero_acc[c][b][j]++;
			}

			if (iter >= burnin) {
				gTheta_samples[c][b][j][iter - burnin] = gTheta[c][b][j];
			}
		}
	}
}

/*
* This is an adaptive MMC algorithm - Algorithm 3 from Ji and Schmidler:
* "Adaptive Markov Chain Monte Carlo for Bayesian Variable Selection"
*/
void c212BB::sample_theta_MIS_Adapt(int c, int burnin, int iter)
{
	int b = 0, j = 0;
	for (b = 0; b < gNumBodySys; b++) {
		for (j = 0; j < gNAE[b]; j++) {

			double w0 = gW0[b][j];
			//double* w = gW[b][j];
			//double* mu = gMU[b][j];
			//double* sigma2 = gSIGMA2[b][j];
			//double lambda = gLAMBDA[b][j];

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			double cand = sample_qn(b, j);

			double ratio = 0;

			double f_cand = log_f_theta(c, b, j, cand);
			double f_prev = log_f_theta(c, b, j, gTheta[c][b][j]);

			//double p = exp(f_cand - f_prev);

			double sn_cand = sn(cand, b, j);
			double sn_prev = sn(gTheta[c][b][j], b, j);

			if (gTheta[c][b][j] == 0.0 && cand == 0.0) {
				ratio = 1;
			} else {
				if (gTheta[c][b][j] != 0.0 && cand == 0.0) {
					ratio = exp(f_cand - f_prev  + log(sn_prev) - log(w0));

				} else {
					if (gTheta[c][b][j] == 0.0 && cand != 0.0) {
						ratio = exp(f_cand - f_prev  + w0 - sn_cand);
					} else {
						ratio = exp(f_cand - f_prev  + sn_prev - sn_cand);
					}
				}
			}


			if (u <= ratio) {
				gTheta[c][b][j] = cand;
				if (iter >= burnin)
					gTheta_acc[c][b][j] = gTheta_acc[c][b][j] + 1;
			}
			else {
			}

			if (iter >= burnin) {
				gTheta_samples[c][b][j][iter - burnin] = gTheta[c][b][j];
			}

			// Update the parameters
			update_params(gTheta[c][b][j], b, j, iter);
		}
	}
}

double c212BB::sample_qn(int b, int j)
{
//	int m = gM;
	double w0 = gW0[b][j];
	double* w = gW[b][j];
	double* mu = gMU[b][j];
	double* sigma2 = gSIGMA2[b][j];

	double cumsum = w0;

	//double lambda = gLAMBDA[b][j];
	//double mu_tilde = gMU_tilde[b][j];
	//double sigma2_tilde = gSIGMA2_tilde[b][j];

	double q = 0.0;

#ifdef INDIVIDUAL_RNG
	GetRNGstate();
#endif
    double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
	PutRNGstate();
#endif

	// Use "<=" here rather than the "<" used in the rest of the code, it makes the code
	// a bit shorter and is correct.
	if (u <= w0) {
		q = 0.0;
	}
	else {
#ifdef INDIVIDUAL_RNG
		GetRNGstate();
#endif
		int i = 0;
		for (i = 0; i < gM; i++) {
			cumsum = cumsum + w[i];
			if (u <= cumsum) {
				q = rnorm(mu[i], sqrt(sigma2[i]));
				break;
			}
		}
#ifdef INDIVIDUAL_RNG
		PutRNGstate();
#endif
	}

	return(q);
}

double c212BB::sn(double x, int b, int j)
{
	int m = gM;
//	double w0 = gW0[b][j];
	double* w = gW[b][j];
	double* mu = gMU[b][j];
	double* sigma2 = gSIGMA2[b][j];

	//double lambda = gLAMBDA[b][j];
	//double mu_tilde = gMU_tilde[b][j];
	//double sigma2_tilde = gSIGMA2_tilde[b][j];

    double q = 0.0;

    //double s[m];
	double* s = (double *)malloc(m * sizeof(double));

	int i = 0;
	for (i = 0; i < m; i++) {
    	s[i] = (1.0/sqrt(2.0*M_PI*sigma2[i]))
					*exp((-1.0/2.0)*(pow((x - mu[i]),2.0))/(sigma2[i]));
	}

	for (i = 0; i < m ; i++) {
    	q +=  w[i]*s[i];
	}

	free(s);

    return(q);
}


double c212BB::phi (double val, double mu, double sigma2)
{

	double s = (1.0/sqrt(2.0*M_PI*sigma2))*exp((-1.0/2.0)*(pow((val - mu),2.0))/(sigma2));
	return(s);
}

void c212BB::update_params(double x, int b, int j, int n)
{
	int m = gM;
	double w0 = gW0[b][j];
	double* w = gW[b][j];
	double* mu = gMU[b][j];
	double* sigma2 = gSIGMA2[b][j];

	double O0 = 0.0;
	double* O = (double *)malloc(m * sizeof(double));
	double O_bar = 0.0;
	double* k = (double *)malloc(m * sizeof(double));

	double r = 0.1/((double)n + 1.0);

	double s = 0.0;

	int i = 0;
	for (i = 0; i < m; i++) {
		s += w[i]*phi(x, mu[i], sigma2[i]);
	}

	if (x == 0.0) {
		O0 = 1/w0;
	}
	else {
		O0 = 0.0;
	}

	for (i = 0; i < m; i++) {
		if (x == 0.0) {
			O[i] = 0.0;
		} else {
			O[i] = phi(x, mu[i], sigma2[i])/s;
		}
	}

	double t = 0;
	for (i = 0; i < m;  i++) {
		t += O[i];
	}

	O_bar = (1.0/((double)m + 1.0)) * (O0 + t);

	for (i = 0; i < m;  i++) {
		k[i] = r * w[i] * O[i];
	}

	w0 = w0 + r*(O0 - O_bar);

	for (i = 0; i < m;  i++) {
		w[i] = w[i] + r*(O[i] - O_bar);
		if (x == 0.0) {
			mu[i] = mu[i];
			sigma2[i] = sigma2[i];
		} else {
			sigma2[i] = sigma2[i] + k[i]*(pow((x - mu[i]),2.0) - sigma2[i]);
			mu[i] = mu[i] + k[i]*(x - mu[i]);
		}
	}

	gW0[b][j] = w0;
	for (i = 0; i < m; i++) {
		gW[b][j][i] = w[i];
		gMU[b][j][i] = mu[i];
		gSIGMA2[b][j][i] = sigma2[i];
	}
	free(O);
	free(k);

}

double c212BB::log_g(double val, int b, int j) {

    double f = 0;

    if (val == 0) {
        f = log(0.5);
    } else {
        f = log(0.5) + log(1.0/(gSigma_MH_theta[b][j]*sqrt(2.0*M_PI))) + (-1.0/2.0)*(pow(val,2.0))/(pow(gSigma_MH_theta[b][j],2.0));
    }

    return(f);
}


void c212BB::sample_theta_Independent_MH(int c, int burnin, int iter)
{
	int b = 0, j = 0;
	for (b = 0; b < gNumBodySys; b++) {
		for ( j = 0; j < gNAE[b]; j++) {

			// Independent MH Algorithm
			// g(y) ~ 0.5 I(theta = 0) + (0.5) N(0, 10^2)

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			double u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif

			double cand = 0.0;

			if (u < 0.5) {
				cand = 0.0;
			}
			else {
#ifdef INDIVIDUAL_RNG
				GetRNGstate();
#endif
				cand = rnorm(0, gSigma_MH_theta[b][j]);
#ifdef INDIVIDUAL_RNG
				PutRNGstate();
#endif
			}

			double f1 = log_f_theta(c, b, j, cand);
			double f2 = log_f_theta(c, b, j, gTheta[c][b][j]);

			double q1 = log_g(cand, b, j);
			double q2 = log_g(gTheta[c][b][j], b, j);

			double lratio = f1 - f2 + q2 - q1;
			double ratio = exp(lratio);

#ifdef INDIVIDUAL_RNG
			GetRNGstate();
#endif
			u = runif(0, 1);
#ifdef INDIVIDUAL_RNG
			PutRNGstate();
#endif


			if (u <= ratio) {
				gTheta[c][b][j] = cand;
				if (iter >= burnin)
					gTheta_acc[c][b][j] = gTheta_acc[c][b][j] + 1;
			}

			if (iter >= burnin) {
				gTheta_samples[c][b][j][iter - burnin] = gTheta[c][b][j];
			}
		}
	}
}

void c212BB::init(SEXP pChains, SEXP pBurnin, SEXP pIter, SEXP pNumBodySys,
					SEXP pMaxAEs, SEXP pNAE,
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
                    SEXP psigma2_theta,
					SEXP palpha_pi,
					SEXP pbeta_pi,
					SEXP plambda_alpha,
					SEXP plambda_beta,
					SEXP pPi,
					SEXP padapt_phase,
					SEXP sim_type,
					SEXP global_sim_params,
					SEXP sim_params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP adapt_min_w,
					SEXP adapt_chains,
					SEXP adapt_burnin,
					SEXP adapt_iter)
{
	release();
	c2121a::release();

	gChains = *(INTEGER(pChains));
	gBurnin = *(INTEGER(pBurnin));
	gIter = *(INTEGER(pIter));
	gNumBodySys = *(INTEGER(pNumBodySys));
	gMaxAEs = *(INTEGER(pMaxAEs));

	alpha_pi_acc = NULL;
	beta_pi_acc = NULL;
	alpha_pi_acc_adapt = 0;
	beta_pi_acc_adapt = 0;

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

	alpha_pi = (double *)malloc(gChains * sizeof(double));
	beta_pi = (double *)malloc(gChains * sizeof(double));
	alpha_pi_acc = (int *)malloc(gChains * sizeof(int));
	beta_pi_acc = (int *)malloc(gChains * sizeof(int));
	lambda_alpha = *(REAL(plambda_alpha));
	lambda_beta = *(REAL(plambda_beta));

	double *valpha_pi = REAL(palpha_pi);
	double *vbeta_pi = REAL(pbeta_pi);

	for (c = 0; c < gChains; c++) {
		alpha_pi[c] = *valpha_pi++;
		beta_pi[c] = *vbeta_pi++;
		alpha_pi_acc[c] = 0;
		beta_pi_acc[c] = 0;
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

	gTheta_zero_prop = (int ***)malloc(gChains * sizeof(int**));
	gTheta_zero_acc = (int ***)malloc(gChains * sizeof(int**));

	for (c = 0; c < gChains; c++) {
		gTheta[c] = (double **)malloc(gNumBodySys * sizeof(double*));
		gGamma[c] = (double **)malloc(gNumBodySys * sizeof(double*));
		gTheta_acc[c] = (int **)malloc(gNumBodySys * sizeof(int*));
		gGamma_acc[c] = (int **)malloc(gNumBodySys * sizeof(int*));

		gTheta_zero_prop[c] = (int **)malloc(gNumBodySys * sizeof(int*));
		gTheta_zero_acc[c] = (int **)malloc(gNumBodySys * sizeof(int*));

		int i = 0;
		for (i = 0; i < gNumBodySys; i++) {
			gTheta[c][i] = (double *)malloc(gNAE[i] * sizeof(double));
			gGamma[c][i] = (double *)malloc(gNAE[i] * sizeof(double));
			gTheta_acc[c][i] = (int *)malloc(gNAE[i] * sizeof(int));
			gGamma_acc[c][i] = (int *)malloc(gNAE[i] * sizeof(int));

			gTheta_zero_prop[c][i] = (int *)malloc(gNAE[i] * sizeof(int));
			gTheta_zero_acc[c][i] = (int *)malloc(gNAE[i] * sizeof(int));
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

					gTheta_zero_prop[c][i][j] = 0;
					gTheta_zero_acc[c][i][j] = 0;
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

	gPi = (double **)malloc(gChains * sizeof(double*));

	for (c = 0; c < gChains; c++) {
		mu_theta[c] = (double *)malloc(gNumBodySys * sizeof(double));
		mu_gamma[c] = (double *)malloc(gNumBodySys * sizeof(double));
		sigma2_theta[c] = (double *)malloc(gNumBodySys * sizeof(double));
		sigma2_gamma[c] = (double *)malloc(gNumBodySys * sizeof(double));

		gPi[c] = (double *)malloc(gNumBodySys * sizeof(double));
	}

	double *vmu_theta = REAL(pmu_theta);
	double *vmu_gamma = REAL(pmu_gamma);
	double *vsigma2_gamma = REAL(psigma2_gamma);
	double *vsigma2_theta = REAL(psigma2_theta);
	double *vPi = REAL(pPi);

	for (c = 0; c < gChains; c++) {
		for (i = 0; i < gNumBodySys; i++) {
			mu_theta[c][i] = *vmu_theta++;
			mu_gamma[c][i] = *vmu_gamma++;
			gPi[c][i] = *vPi++;
			sigma2_theta[c][i] = *vsigma2_theta++;
			sigma2_gamma[c][i] = *vsigma2_gamma++;
		}
	}

	// The samples
	mu_gamma_0_samples = (double **)malloc(gChains * sizeof(double*));
	mu_theta_0_samples = (double **)malloc(gChains *sizeof(double*));
	tau2_theta_0_samples = (double **)malloc(gChains *sizeof(double*));
	tau2_gamma_0_samples = (double **)malloc(gChains *sizeof(double*));
	alpha_pi_samples = (double **)malloc(gChains *sizeof(double*));
	beta_pi_samples = (double **)malloc(gChains *sizeof(double*));

	for (c = 0; c < gChains; c++) {
		mu_gamma_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		mu_theta_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		tau2_theta_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		tau2_gamma_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		alpha_pi_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		beta_pi_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
	}

	mu_theta_samples = (double ***)malloc(gChains *sizeof(double**));
	mu_gamma_samples = (double ***)malloc(gChains *sizeof(double**));
	sigma2_theta_samples = (double ***)malloc(gChains *sizeof(double**));
	sigma2_gamma_samples = (double ***)malloc(gChains *sizeof(double**));

	gPi_samples = (double ***)malloc(gChains *sizeof(double**));

	for (c = 0; c < gChains; c++) {
		mu_theta_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		mu_gamma_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		sigma2_theta_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		sigma2_gamma_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		gPi_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		for (i = 0; i < gNumBodySys; i++) {
			mu_theta_samples[c][i] =
						(double *)malloc((gIter - gBurnin) *sizeof(double));
			mu_gamma_samples[c][i] =
						(double *)malloc((gIter - gBurnin) *sizeof(double));
			sigma2_theta_samples[c][i] =
						(double *)malloc((gIter - gBurnin) *sizeof(double));
			sigma2_gamma_samples[c][i] =
						(double *)malloc((gIter - gBurnin) *sizeof(double));

			gPi_samples[c][i] = (double *)malloc((gIter - gBurnin) *sizeof(double));
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
				gTheta_samples[c][i][j] = (double *)malloc((gIter - gBurnin) *sizeof(double));
				gGamma_samples[c][i][j] = (double *)malloc((gIter - gBurnin) *sizeof(double));
			}
		}
	}

	// Simulation Parameters

	// Global Parameters
	initGlobalSimParams(sim_type, global_sim_params);

	initSimParams(sim_params);

	// MH point-mass weights
	gMH_weight = *(REAL(MH_weight));
	gMHAdaptParams.w_min = *(REAL(adapt_min_w));;
	gMHAdaptParams.chains = *(INTEGER(adapt_chains));
	gMHAdaptParams.burnin = *(INTEGER(adapt_burnin));
	gMHAdaptParams.iter = *(INTEGER(adapt_iter));

	// The adaptive chains must have at most the same number of chains as the
	// run being adapted
	if (gMHAdaptParams.chains > gChains)
		gMHAdaptParams.chains = gChains;

	initPMWeights(pm_weights);

	// Adaptive MCMC
	gW0 = (double **)malloc(gNumBodySys *sizeof(double*));
	for (i = 0; i < gNumBodySys; i++) {
		gW0[i] = (double *)malloc(gNAE[i] * sizeof(double));
		for (j = 0; j < gNAE[i]; j++) {
			gW0[i][j] = 0.5;
		}
	}

	gW = (double ***)malloc(gNumBodySys *sizeof(double**));
	for (i = 0; i < gNumBodySys; i++) {
		gW[i] = (double **)malloc(gNAE[i] * sizeof(double*));
		for (j = 0; j < gNAE[i]; j++) {
			gW[i][j] = (double *)malloc(gM*sizeof(double));
			int k = 0;
			for (k = 0; k < gM; k++) {
				// Equally wieghted over the M normal distributions
				// currently M (gM) is 1.
				gW[i][j][k] = 0.5/gM;
			}
		}
	}

	gMU = (double ***)malloc(gNumBodySys *sizeof(double**));
	for (i = 0; i < gNumBodySys; i++) {
		gMU[i] = (double **)malloc(gNAE[i] * sizeof(double*));
		for (j = 0; j < gNAE[i]; j++) {
			gMU[i][j] = (double *)malloc(gM*sizeof(double));
			int k = 0;
			// Currently M (gM) = 1 so this gMU is 0 for all AEs
			for (k = 0; k < gM; k++) {
				//gMU[i][j][k] = 0;
				gMU[i][j][k] = k + 1;
			}
		}
	}

	gSIGMA2 = (double ***)malloc(gNumBodySys *sizeof(double**));
	for (i = 0; i < gNumBodySys; i++) {
		gSIGMA2[i] = (double **)malloc(gNAE[i] * sizeof(double*));
		for (j = 0; j < gNAE[i]; j++) {
			gSIGMA2[i][j] = (double *)malloc(gM*sizeof(double));
			int k = 0;
			for (k = 0; k < gM; k++) {
				gSIGMA2[i][j][k] = 100;
				gSIGMA2[i][j][k] = 10;
				gSIGMA2[i][j][k] = gDefault_Sigma_MH_theta;
				//gSIGMA2[i][j][k] = 0.5;
				//gSIGMA2[i][j][k] = 0.2;
				//gSIGMA2[i][j][k] = 1.0;
			}
		}
	}

	//gMU_tilde = (double **)malloc(gNumBodySys *sizeof(double*));
	//for (i = 0; i < gNumBodySys; i++) {
	//	gMU_tilde[i] = (double *)malloc(pNAE[i] * sizeof(double));
	//	for (j = 0; j < pNAE[i]; j++) {
	//			gMU_tilde[i][j] = -1.0;
	//	}
	//}

	//gSIGMA2_tilde = (double **)malloc(gNumBodySys *sizeof(double*));
	//for (i = 0; i < gNumBodySys; i++) {
	//	gSIGMA2_tilde[i] = (double *)malloc(pNAE[i] * sizeof(double));
	//	for (j = 0; j < pNAE[i]; j++) {
	//			gSIGMA2_tilde[i][j] = 0.5;
	//			gSIGMA2_tilde[i][j] = 10;
	//	}
	//}

	//gLAMBDA = (double **)malloc(gNumBodySys *sizeof(double*));
	//for (i = 0; i < gNumBodySys; i++) {
	//	gLAMBDA[i] = (double *)malloc(pNAE[i] * sizeof(double));
	//	for (j = 0; j < pNAE[i]; j++) {
	//			gLAMBDA[i][j] = 0.1;
	//			//gLAMBDA[i][j] = 0.0;
	//	}
	//}

	gAdapt_Phase_alpha = *(INTEGER(padapt_phase));
	gAdapt_Phase_beta = *(INTEGER(padapt_phase));
}

void c212BB::init(int pChains, int pBurnin, int pIter, int pNumBodySys,
					int pMaxAEs, int* pNAE,
					int** pX, int** pY,
					int** pNC, int** pNT,
					double*** ptheta, double*** pgamma,
                    double pmu_gamma_0_0,
                    double ptau2_gamma_0_0,
                    double pmu_theta_0_0,
                    double ptau2_theta_0_0,
                    double palpha_gamma_0_0,
                    double pbeta_gamma_0_0,
                    double palpha_theta_0_0,
                    double pbeta_theta_0_0,
                    double palpha_gamma,
                    double pbeta_gamma,
                    double palpha_theta,
                    double pbeta_theta,
                    double* pmu_gamma_0,
					double* ptau2_gamma_0,
                    double* pmu_theta_0,
                    double* ptau2_theta_0,
                    double** pmu_gamma,
                    double** pmu_theta,
                    double** psigma2_gamma,
                    double** psigma2_theta,
					double* palpha_pi,
					double* pbeta_pi,
					double plambda_alpha,
					double plambda_beta,
					double** pPi,
					int padapt_phase,
					const char* sim_type,
					std::map<const char*, gSimParams>& mGlobalSimParams,
					double** gamma_params,
					int** gamma_ctrl,
					double** theta_params,
					double MH_weight,
					double** pm_weights,
					MHAdaptParams& adapt_params)
{
	release();
	c2121a::release();

	gChains = pChains;
	gBurnin = pBurnin;
	gIter = pIter;

	gNumBodySys = pNumBodySys;
	gMaxAEs = pMaxAEs;

	alpha_pi_acc = NULL;
	beta_pi_acc = NULL;
	alpha_pi_acc_adapt = 0;
	beta_pi_acc_adapt = 0;

	int l = 0;
	gNAE = (int *)malloc(gNumBodySys * sizeof(int));
	for (l = 0; l < gNumBodySys; l++) {
		gNAE[l] = pNAE[l];
	}

	alpha_gamma_0_0 = palpha_gamma_0_0;
	beta_gamma_0_0 = pbeta_gamma_0_0;
	alpha_theta_0_0 = palpha_theta_0_0;
	beta_theta_0_0 = pbeta_theta_0_0;

	alpha_gamma = palpha_gamma;
	beta_gamma = pbeta_gamma;
	alpha_theta = palpha_theta;
	beta_theta = pbeta_theta;

	mu_theta_0_0 = pmu_theta_0_0;
	mu_gamma_0_0 = pmu_gamma_0_0;
	tau2_theta_0_0 = ptau2_theta_0_0;
	tau2_gamma_0_0 = ptau2_gamma_0_0;

	int c = 0;
	mu_gamma_0 = (double *)malloc(gChains * sizeof(double));
	mu_theta_0 = (double *)malloc(gChains * sizeof(double));
	tau2_gamma_0 = (double *)malloc(gChains * sizeof(double));
	tau2_theta_0 = (double *)malloc(gChains * sizeof(double));


	for (c = 0; c < gChains; c++) {
		mu_gamma_0[c] = pmu_gamma_0[c];
		mu_theta_0[c] = pmu_theta_0[c];
		tau2_gamma_0[c] = ptau2_gamma_0[c];
		tau2_theta_0[c] = ptau2_theta_0[c];
	}

	alpha_pi = (double *)malloc(gChains * sizeof(double));
	beta_pi = (double *)malloc(gChains * sizeof(double));
	alpha_pi_acc = (int *)malloc(gChains * sizeof(int));
	beta_pi_acc = (int *)malloc(gChains * sizeof(int));
	lambda_alpha = plambda_alpha;
	lambda_beta = plambda_beta;

	for (c = 0; c < gChains; c++) {
		alpha_pi[c] = palpha_pi[c];
		beta_pi[c] = pbeta_pi[c];
		alpha_pi_acc[c] = 0;
		beta_pi_acc[c] = 0;
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

	gTheta_zero_prop = (int ***)malloc(gChains * sizeof(int**));
	gTheta_zero_acc = (int ***)malloc(gChains * sizeof(int**));

	for (c = 0; c < gChains; c++) {
		gTheta[c] = (double **)malloc(gNumBodySys * sizeof(double*));
		gGamma[c] = (double **)malloc(gNumBodySys * sizeof(double*));
		gTheta_acc[c] = (int **)malloc(gNumBodySys * sizeof(int*));
		gGamma_acc[c] = (int **)malloc(gNumBodySys * sizeof(int*));

		gTheta_zero_prop[c] = (int **)malloc(gNumBodySys * sizeof(int*));
		gTheta_zero_acc[c] = (int **)malloc(gNumBodySys * sizeof(int*));

		int i = 0;
		for (i = 0; i < gNumBodySys; i++) {
			gTheta[c][i] = (double *)malloc(gNAE[i] * sizeof(double));
			gGamma[c][i] = (double *)malloc(gNAE[i] * sizeof(double));
			gTheta_acc[c][i] = (int *)malloc(gNAE[i] * sizeof(int));
			gGamma_acc[c][i] = (int *)malloc(gNAE[i] * sizeof(int));

			gTheta_zero_prop[c][i] = (int *)malloc(gNAE[i] * sizeof(int));
			gTheta_zero_acc[c][i] = (int *)malloc(gNAE[i] * sizeof(int));
		}
	}

	for (c = 0; c < gChains; c++) {
		int i = 0, j = 0;
		for (i = 0; i < gNumBodySys; i++) {
			for (j = 0; j < gMaxAEs; j++) {
				if (j < gNAE[i]) {
					gTheta[c][i][j] = ptheta[c][i][j];
					gGamma[c][i][j] = pgamma[c][i][j];
					gGamma_acc[c][i][j] = 0;
					gTheta_acc[c][i][j] = 0;

					gTheta_zero_prop[c][i][j] = 0;
					gTheta_zero_acc[c][i][j] = 0;
				}
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
				NC[i][j] = pNC[i][j];
				NT[i][j] = pNT[i][j];
				x[i][j] = pX[i][j];
				y[i][j] = pY[i][j];
			}
		}
	}

	mu_theta = (double **)malloc(gChains * sizeof(double*));
	mu_gamma = (double **)malloc(gChains * sizeof(double*));
	sigma2_theta = (double **)malloc(gChains * sizeof(double*));
	sigma2_gamma = (double **)malloc(gChains * sizeof(double*));

	gPi = (double **)malloc(gChains * sizeof(double*));

	for (c = 0; c < gChains; c++) {
		mu_theta[c] = (double *)malloc(gNumBodySys * sizeof(double));
		mu_gamma[c] = (double *)malloc(gNumBodySys * sizeof(double));
		sigma2_theta[c] = (double *)malloc(gNumBodySys * sizeof(double));
		sigma2_gamma[c] = (double *)malloc(gNumBodySys * sizeof(double));

		gPi[c] = (double *)malloc(gNumBodySys * sizeof(double));
	}

	for (c = 0; c < gChains; c++) {
		for (i = 0; i < gNumBodySys; i++) {
			mu_theta[c][i] = pmu_theta[c][i];
			mu_gamma[c][i] = pmu_gamma[c][i];
			gPi[c][i] = pPi[c][i];
			sigma2_theta[c][i] = psigma2_theta[c][i];
			sigma2_gamma[c][i] = psigma2_gamma[c][i];
		}
	}

	// The samples
	mu_gamma_0_samples = (double **)malloc(gChains * sizeof(double*));
	mu_theta_0_samples = (double **)malloc(gChains *sizeof(double*));
	tau2_theta_0_samples = (double **)malloc(gChains *sizeof(double*));
	tau2_gamma_0_samples = (double **)malloc(gChains *sizeof(double*));
	alpha_pi_samples = (double **)malloc(gChains *sizeof(double*));
	beta_pi_samples = (double **)malloc(gChains *sizeof(double*));

	for (c = 0; c < gChains; c++) {
		mu_gamma_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		mu_theta_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		tau2_theta_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		tau2_gamma_0_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		alpha_pi_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
		beta_pi_samples[c] = (double *)malloc((gIter - gBurnin) *sizeof(double));
	}

	mu_theta_samples = (double ***)malloc(gChains *sizeof(double**));
	mu_gamma_samples = (double ***)malloc(gChains *sizeof(double**));
	sigma2_theta_samples = (double ***)malloc(gChains *sizeof(double**));
	sigma2_gamma_samples = (double ***)malloc(gChains *sizeof(double**));

	gPi_samples = (double ***)malloc(gChains *sizeof(double**));

	for (c = 0; c < gChains; c++) {
		mu_theta_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		mu_gamma_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		sigma2_theta_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		sigma2_gamma_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		gPi_samples[c] = (double **)malloc(gNumBodySys *sizeof(double*));
		for (i = 0; i < gNumBodySys; i++) {
			mu_theta_samples[c][i] =
						(double *)malloc((gIter - gBurnin) *sizeof(double));
			mu_gamma_samples[c][i] =
						(double *)malloc((gIter - gBurnin) *sizeof(double));
			sigma2_theta_samples[c][i] =
						(double *)malloc((gIter - gBurnin) *sizeof(double));
			sigma2_gamma_samples[c][i] =
						(double *)malloc((gIter - gBurnin) *sizeof(double));

			gPi_samples[c][i] = (double *)malloc((gIter - gBurnin) *sizeof(double));
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
				gTheta_samples[c][i][j] = (double *)malloc((gIter - gBurnin) *sizeof(double));
				gGamma_samples[c][i][j] = (double *)malloc((gIter - gBurnin) *sizeof(double));
			}
		}
	}

	// Simulation Parameters

	// Global Parameters
	initGlobalSimParams(sim_type, mGlobalSimParams);

	initSimParams(gamma_params, gamma_ctrl, theta_params);

	// MH
	gMH_weight = MH_weight;
	gMHAdaptParams = adapt_params;

	initPMWeights(pm_weights);

	// Adaptive MCMC
	gW0 = (double **)malloc(gNumBodySys *sizeof(double*));
	for (i = 0; i < gNumBodySys; i++) {
		gW0[i] = (double *)malloc(gNAE[i] * sizeof(double));
		for (j = 0; j < gNAE[i]; j++) {
			gW0[i][j] = 0.5;
		}
	}

	gW = (double ***)malloc(gNumBodySys *sizeof(double**));
	for (i = 0; i < gNumBodySys; i++) {
		gW[i] = (double **)malloc(gNAE[i] * sizeof(double*));
		for (j = 0; j < gNAE[i]; j++) {
			gW[i][j] = (double *)malloc(gM*sizeof(double));
			int k = 0;
			for (k = 0; k < gM; k++) {
				// Equally wieghted over the M normal distributions
				// currently M (gM) is 1.
				gW[i][j][k] = 0.5/gM;
			}
		}
	}

	gMU = (double ***)malloc(gNumBodySys *sizeof(double**));
	for (i = 0; i < gNumBodySys; i++) {
		gMU[i] = (double **)malloc(gNAE[i] * sizeof(double*));
		for (j = 0; j < gNAE[i]; j++) {
			gMU[i][j] = (double *)malloc(gM*sizeof(double));
			int k = 0;
			// Currently M (gM) = 1 so this gMU is 0 for all AEs
			for (k = 0; k < gM; k++) {
				//gMU[i][j][k] = 0;
				gMU[i][j][k] = k + 1;
			}
		}
	}

	gSIGMA2 = (double ***)malloc(gNumBodySys *sizeof(double**));
	for (i = 0; i < gNumBodySys; i++) {
		gSIGMA2[i] = (double **)malloc(gNAE[i] * sizeof(double*));
		for (j = 0; j < gNAE[i]; j++) {
			gSIGMA2[i][j] = (double *)malloc(gM*sizeof(double));
			int k = 0;
			for (k = 0; k < gM; k++) {
				gSIGMA2[i][j][k] = 100;
				gSIGMA2[i][j][k] = 10;
				gSIGMA2[i][j][k] = gDefault_Sigma_MH_theta;
				//gSIGMA2[i][j][k] = 0.5;
				//gSIGMA2[i][j][k] = 0.2;
				//gSIGMA2[i][j][k] = 1.0;
			}
		}
	}

	//gMU_tilde = (double **)malloc(gNumBodySys *sizeof(double*));
	//for (i = 0; i < gNumBodySys; i++) {
	//	gMU_tilde[i] = (double *)malloc(pNAE[i] * sizeof(double));
	//	for (j = 0; j < pNAE[i]; j++) {
	//			gMU_tilde[i][j] = -1.0;
	//	}
	//}

	//gSIGMA2_tilde = (double **)malloc(gNumBodySys *sizeof(double*));
	//for (i = 0; i < gNumBodySys; i++) {
	//	gSIGMA2_tilde[i] = (double *)malloc(pNAE[i] * sizeof(double));
	//	for (j = 0; j < pNAE[i]; j++) {
	//			gSIGMA2_tilde[i][j] = 0.5;
	//			gSIGMA2_tilde[i][j] = 10;
	//	}
	//}

	//gLAMBDA = (double **)malloc(gNumBodySys *sizeof(double*));
	//for (i = 0; i < gNumBodySys; i++) {
	//	gLAMBDA[i] = (double *)malloc(pNAE[i] * sizeof(double));
	//	for (j = 0; j < pNAE[i]; j++) {
	//			gLAMBDA[i][j] = 0.1;
	//			//gLAMBDA[i][j] = 0.0;
	//	}
	//}

	gAdapt_Phase_alpha = padapt_phase;
	gAdapt_Phase_beta = padapt_phase;
}

void c212BB::release()
{
	int i = 0, j = 0, c = 0;

	if (gPi != NULL) {
		for (c = 0; c < gChains; c++) {
			free(gPi[c]);
		}
		free(gPi);
		gPi = NULL;
	}

	if (alpha_pi != NULL) {
		free(alpha_pi);
		alpha_pi = NULL;
	}

	if (beta_pi != NULL) {
		free(beta_pi);
		beta_pi = NULL;
	}

	if (alpha_pi_acc != NULL) {
		free(alpha_pi_acc);
		alpha_pi_acc = NULL;
	}

	if (beta_pi_acc != NULL) {
		free(beta_pi_acc);
		beta_pi_acc = NULL;
	}

	if (gTheta_zero_prop != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(gTheta_zero_prop[c][i]);
			}
			free(gTheta_zero_prop[c]);
		}
		free(gTheta_zero_prop);
		gTheta_zero_prop = NULL;
	}

	if (gTheta_zero_acc != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(gTheta_zero_acc[c][i]);
			}
			free(gTheta_zero_acc[c]);
		}
		free(gTheta_zero_acc);
		gTheta_zero_acc = NULL;
	}


	// The samples
	if (gPi_samples != NULL) {
		for (c = 0; c < gChains; c++) {
			for (i = 0; i < gNumBodySys; i++) {
				free(gPi_samples[c][i]);
			}
			free(gPi_samples[c]);
		}
		free(gPi_samples);
		gPi_samples = NULL;
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

	// MH
	if (gWp != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gWp[i]);
		}
		free(gWp);
	}

	// Adaptive MIS
	if (gW0 != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gW0[i]);
		}
		free(gW0);
		gW0 = NULL;
	}

	if (gW != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			for (j = 0; j < gNAE[i]; j++) {
				free(gW[i][j]);
			}
			free(gW[i]);
		}
		free(gW);
		gW = NULL;
	}

	if (gMU != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			for (j = 0; j < gNAE[i]; j++) {
				free(gMU[i][j]);
			}
			free(gMU[i]);
		}
		free(gMU);
		gMU = NULL;
	}

	if (gSIGMA2 != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			for (j = 0; j < gNAE[i]; j++) {
				free(gSIGMA2[i][j]);
			}
			free(gSIGMA2[i]);
		}
		free(gSIGMA2);
		gSIGMA2 = NULL;
	}

	//if (gMU_tilde != NULL) {
	//	for (i = 0; i < gNumBodySys; i++) {
	//		free(gMU_tilde[i]);
	//	}
	//	free(gMU_tilde);
	//	gMU_tilde = NULL;
	//}

	//if (gSIGMA2_tilde != NULL) {
	//	for (i = 0; i < gNumBodySys; i++) {
	//		free(gSIGMA2_tilde[i]);
	//	}
	//	free(gSIGMA2_tilde);
	//	gSIGMA2_tilde = NULL;
	//}

	//if (gLAMBDA != NULL) {
	//	for (i = 0; i < gNumBodySys; i++) {
	//		free(gLAMBDA[i]);
	//	}
	//	free(gLAMBDA);
	//	gLAMBDA = NULL;
	//}

	if (gAdapt_Phase_theta != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(gAdapt_Phase_theta[i]);
		}
		free(gAdapt_Phase_theta);
		gAdapt_Phase_theta = NULL;
	}

	if (theta_acc_adapt != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(theta_acc_adapt[i]);
		}
		free(theta_acc_adapt);
		theta_acc_adapt = NULL;
	}

	if (theta_mix_p != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(theta_mix_p[i]);
		}
		free(theta_mix_p);
		theta_mix_p = NULL;
	}

	if (theta_max_p != NULL) {
		for (i = 0; i < gNumBodySys; i++) {
			free(theta_max_p[i]);
		}
		free(theta_max_p);
		theta_max_p = NULL;
	}
}

SEXP c212BB::getPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL2Samples(gPi_samples);

	return samples;
}

SEXP c212BB::getAlphaPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(alpha_pi_samples);

	return samples;
}

SEXP c212BB::getBetaPiSamples()
{
	SEXP samples = R_NilValue;

	samples = getL3Samples(beta_pi_samples);

	return samples;
}

SEXP c212BB::getL3Accept(int* &data)
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

SEXP c212BB::getAlphaPiAccept()
{
	SEXP acc = R_NilValue;

	acc = getL3Accept(alpha_pi_acc);

	return acc;
}

SEXP c212BB::getBetaPiAccept()
{
	SEXP acc = R_NilValue;

	acc = getL3Accept(beta_pi_acc);

	return acc;
}

void c212BB::getPiSamples(int* c, int* b, double* p)
{
	int i = 0;
	int C = (*c) - 1;
	int B = (*b) - 1;
	for (i = 0; i < (gIter - gBurnin); i++) {
		p[i] = gPi_samples[C][B][i];
	}
}

void c212BB::getAlphaPiSamples(int* c, double* alpha_pi)
{
	int i = 0;
	int C = (*c) - 1;
	for (i = 0; i < (gIter - gBurnin); i++) {
		alpha_pi[i] = alpha_pi_samples[C][i];
	}
}

void c212BB::getBetaPiSamples(int* c, double* beta_pi)
{
	int i = 0;
	int C = (*c) - 1;
	for (i = 0; i < (gIter - gBurnin); i++) {
		beta_pi[i] = beta_pi_samples[C][i];
	}
}

void c212BB::getAlphaPiAccept(int* c, double* acc)
{
	int C = (*c) - 1;
	*acc = alpha_pi_acc[C];
}

void c212BB::getBetaPiAccept(int* c, double* acc)
{
	int C = (*c) - 1;
	*acc = beta_pi_acc[C];
}

SEXP c212BB::getThetaZeroAccept()
{
	SEXP acc = R_NilValue;

	acc = getL1Accept(gTheta_zero_acc);

	return acc;
}

SEXP c212BB::getThetaZeroProp()
{
	SEXP acc = R_NilValue;

	acc = getL1Accept(gTheta_zero_prop);

	return acc;
}

void c212BB::getThetaZeroAccept(int* c, int* b, int* j, double* zero_prop, double* zero_acc)
{
	int C = (*c) - 1;
	int B = (*b) - 1;
	int J = (*j) - 1;

    *zero_prop = gTheta_zero_prop[C][B][J];
    *zero_acc = gTheta_zero_acc[C][B][J];
}
