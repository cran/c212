#ifndef C212_EXEC_H
#define C212_EXEC_H

extern "C" {

SEXP c2121a_exec(SEXP pChains, SEXP pBurnin, SEXP pIter, SEXP pNumBodySys,
					SEXP pMaxAEs,
					SEXP pNAE, SEXP sim_type,
					SEXP pGlobal_Sim_Param, SEXP pGlobal_Sim_Param_Cntrl,
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
                    SEXP psigma2_theta);

void Release();

SEXP getThetaSamplesAll();
SEXP getGammaSamplesAll();
SEXP getMuThetaSamplesAll();
SEXP getMuGammaSamplesAll();
SEXP getSigma2ThetaSamplesAll();
SEXP getSigma2GammaSamplesAll();
SEXP getMuGamma0SamplesAll();
SEXP getMuTheta0SamplesAll();
SEXP getTau2Gamma0SamplesAll();
SEXP getTau2Theta0SamplesAll();
SEXP getThetaAcceptAll();
SEXP getGammaAcceptAll();
SEXP getThetaZeroPropAll();
SEXP getThetaZeroAcceptAll();

void getThetaSamples(int* c, int* b, int* j, double* theta_samples);
void getGammaSamples(int* c, int* b, int* j, double* gamma_samples);
void getMuThetaSamples(int* c, int* b, double* mu_theta);
void getMuGammaSamples(int* c, int* b, double* mu_gamma);
void getSigma2ThetaSamples(int* c, int* b, double* mu_theta);
void getSigma2GammaSamples(int* c, int* b, double* mu_gamma);
void getMuGamma0Samples(int* c, double* mu);
void getMuTheta0Samples(int* c, double* mu);
void getTau2Gamma0Samples(int* c, double* tau2);
void getTau2Theta0Samples(int* c, double* tau2);
void getThetaAccept(int* c, int* b, int* j, double* theta_acc);
void getGammaAccept(int* c, int* b, int* j, double* gamma_acc);

void getThetaZeroAccept(int* c, int* b, int* j, double* zero_prop, double* zero_acc);

SEXP c212BB_exec(SEXP sChains, SEXP sBurnin, SEXP sIter,  SEXP sNumBodySys,
					SEXP sMaxAEs, SEXP pNAE,
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
					SEXP adapt_iter);

SEXP getPiSamplesAll();
SEXP getAlphaPiSamplesAll();
SEXP getBetaPiSamplesAll();
SEXP getAlphaPiAcceptAll();
SEXP getBetaPiAcceptAll();

void getPiSamples(int* c, int* b, double* pi);
void getAlphaPiSamples(int* c, double* alpha_pi);
void getBetaPiSamples(int* c, double* beta_pi);
void getAlphaPiAccept(int* c, double* acc);
void getBetaPiAccept(int* c, double* acc);

SEXP c2121a_poisson_mc_exec(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_ctrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sIndep_Intervals,
					SEXP sMaxBs, SEXP sNumBodySys,
					SEXP sMaxAEs, SEXP sNAE,
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
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta);

SEXP c2121a_interim_hier2_exec(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_ctrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sIndep_Intervals,
					SEXP sMaxBs, SEXP sNumBodySys,
					SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta);


SEXP getThetaSamplesInterimAll();
SEXP getGammaSamplesInterimAll();
SEXP getMuThetaSamplesInterimAll();
SEXP getMuGammaSamplesInterimAll();
SEXP getSigma2ThetaSamplesInterimAll();
SEXP getSigma2GammaSamplesInterimAll();
SEXP getThetaAcceptInterimAll();
SEXP getGammaAcceptInterimAll();
SEXP getMuGamma0SamplesInterimAll();
SEXP getMuTheta0SamplesInterimAll();
SEXP getTau2Gamma0SamplesInterimAll();
SEXP getTau2Theta0SamplesInterimAll();

void getThetaSamplesInterim(int *c, int* l, int* b, int* j, double* theta_samples);
void getGammaSamplesInterim(int *c, int* l, int* b, int* j, double* theta_samples);
void getMuThetaSamplesInterim(int *c, int *l, int* b, double* mu_theta);
void getMuGammaSamplesInterim(int *c, int *l, int* b, double* mu_gamma);
void getSigma2ThetaSamplesInterim(int *c, int *l, int* b, double* mu_theta);
void getSigma2GammaSamplesInterim(int *c, int *l, int* b, double* mu_gamma);
void getMuGamma0SamplesInterim(int *c, int *l, double* mu);
void getMuTheta0SamplesInterim(int *c, int *l, double* mu);
void getTau2Gamma0SamplesInterim(int *c, int *l, double* tau2);
void getTau2Theta0SamplesInterim(int *c, int *l, double* tau2);
void getThetaAcceptInterim(int *c, int *l, int* b, int* j, double* theta_acc);
void getGammaAcceptInterim(int *c, int *l, int* b, int* j, double* gamma_acc);
void Release_Interim();

SEXP c212BB_poisson_mc_exec(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model, SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sIndep_Intervals,
					SEXP sMaxBs, SEXP sNumBodySys,
					SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma, SEXP pmu_gamma_0_0,
					SEXP ptau2_gamma_0_0, SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0, SEXP palpha_gamma_0_0,
					SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0, SEXP pbeta_theta_0_0, SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta, SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0, SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi,
					SEXP palpha_pi,
					SEXP pbeta_pi,
					SEXP plambda_alpha,
					SEXP plambda_beta,
					SEXP palgo,
					SEXP padapt_phase);

SEXP c212BB_interim_hier2_exec(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobalSim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sIndep_Intervals,
					SEXP sMaxBs, SEXP sNumBodySys,
					SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta, SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0, SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi,
					SEXP palpha_pi,
					SEXP pbeta_pi,
					SEXP palgo,
					SEXP padapt_phase);

SEXP getPiSamplesInterimAll();
SEXP getAlphaPiSamplesInterimAll();
SEXP getBetaPiSamplesInterimAll();
SEXP getAlphaPiAcceptInterimAll();
SEXP getBetaPiAcceptInterimAll();

void getPiSamplesInterim(int* c, int* l, int* b, double* pi);
void getAlphaPiSamplesInterim(int *c, int* l, double* alpha_pi);
void getBetaPiSamplesInterim(int *c, int* l, double* beta_pi);
void getAlphaPiAcceptInterim(int *c, int* l, double* acc);
void getBetaPiAcceptInterim(int *c, int* l, double* acc);

}

#endif
