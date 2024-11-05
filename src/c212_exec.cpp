#include<cstdio>
#include<cstdlib>

#include<map>

#include "c212_Rdefines.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "c212_exec.h"
#include "c2121a.h"
#include "c212BB.h"
#include "c2121a_poisson_mc_hier2_lev0.h"
#include "c2121a_poisson_mc_hier2_lev1.h"
#include "c212BB_poisson_mc_hier2_lev0.h"
#include "c212BB_poisson_mc_hier2_lev1.h"
#include "c2121a_poisson_mc_hier3_lev0.h"
#include "c2121a_poisson_mc_hier3_lev2.h"
#include "c2121a_poisson_mc_hier3_lev1.h"
#include "c212BB_poisson_mc_hier3_lev0.h"
#include "c212BB_poisson_mc_hier3_lev2.h"
#include "c212BB_poisson_mc_hier3_lev1.h"

//#include <pthread.h>
//#include "c2121a_poisson_mt.h"

using namespace std;

//static const char *rcsId = "$Id: c212_exec.cpp,v 1.26 2018/09/21 10:30:19 clb13102 Exp clb13102 $";

// These should really come from a common base class or be static within the class
static c2121a* model = NULL;
static c2121a_poisson_mc_hier2_lev0* model_interim = NULL;

SEXP c2121a_exec(SEXP pChains, SEXP pBurnin, SEXP pIter, SEXP pNumBodySys,
				SEXP pMaxAEs, SEXP pNAE, SEXP sim_type,
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
                    SEXP psigma2_theta)
{
	try {

		if (model) {
			delete model;
			model = NULL;
		}

		model = new c2121a(pChains, pBurnin, pIter, pNumBodySys, pMaxAEs, pNAE,
					sim_type, pGlobal_Sim_Param, pGlobal_Sim_Param_Cntrl,
					pSim_Param,
					pX, pY,
					pNC, pNT,
					ptheta, pgamma,
					pmu_gamma_0_0,
                    ptau2_gamma_0_0,
					pmu_theta_0_0,
					ptau2_theta_0_0,
					palpha_gamma_0_0,
					pbeta_gamma_0_0,
					palpha_theta_0_0,
					pbeta_theta_0_0,
					palpha_gamma,
					pbeta_gamma,
					palpha_theta,
					pbeta_theta,
					pmu_gamma_0,
					ptau2_gamma_0,
					pmu_theta_0,
					ptau2_theta_0,
                    pmu_gamma,
					pmu_theta,
					psigma2_gamma,
					psigma2_theta); 

		model->gibbs_sampler();
	}
	catch(...) {
		Rprintf("Unexpected exception c2121a_exec\n");
	}
	return R_NilValue;
}

SEXP getThetaSamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getThetaSamples();

	return samples;
}

SEXP getGammaSamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getGammaSamples();

	return samples;
}
SEXP getMuGammaSamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getMuGammaSamples();

	return samples;
}

SEXP getMuThetaSamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getMuThetaSamples();

	return samples;
}
SEXP getSigma2GammaSamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getSigma2GammaSamples();

	return samples;
}

SEXP getSigma2ThetaSamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getSigma2ThetaSamples();

	return samples;
}

SEXP getMuGamma0SamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getMuGamma0Samples();

	return samples;
}

SEXP getMuTheta0SamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getMuTheta0Samples();

	return samples;
}

SEXP getTau2Gamma0SamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getTau2Gamma0Samples();

	return samples;
}

SEXP getTau2Theta0SamplesAll()
{
	SEXP samples = R_NilValue;

	if (model)
		samples = model->getTau2Theta0Samples();

	return samples;
}

SEXP getThetaAcceptAll()
{
	SEXP acc = R_NilValue;

	if (model)
		acc = model->getThetaAccept();

	return acc;
}

SEXP getGammaAcceptAll()
{
	SEXP acc = R_NilValue;

	if (model)
		acc = model->getGammaAccept();

	return acc;
}

SEXP getThetaZeroAcceptAll()
{
	SEXP acc = R_NilValue;

	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		acc = modelBB->getThetaZeroAccept();

	return acc;
}

SEXP getThetaZeroPropAll()
{
	SEXP acc = R_NilValue;

	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		acc = modelBB->getThetaZeroProp();

	return acc;
}

void Release()
{
	if (model)
		delete model;
	model = NULL;
}

void getThetaSamples(int* c, int* b, int* j, double* theta_samples)
{
	if (model)
		model->getThetaSamples(c, b ,j, theta_samples);
}

void getGammaSamples(int* c, int* b, int* j, double* gamma_samples)
{
	if (model)
		model->getGammaSamples(c, b, j, gamma_samples);
}

void getMuThetaSamples(int* c, int* b, double* mu_theta)
{
	if (model)
		model->getMuThetaSamples(c, b, mu_theta);
}

void getMuGammaSamples(int* c, int* b, double* mu_gamma)
{
	if (model)
		model->getMuGammaSamples(c, b, mu_gamma);
}

void getSigma2ThetaSamples(int* c, int* b, double* sigma2)
{
	if (model)
		model->getSigma2ThetaSamples(c, b, sigma2);
}

void getSigma2GammaSamples(int* c, int* b, double* sigma2)
{
	if (model)
		model->getSigma2GammaSamples(c, b, sigma2);
}

void getMuGamma0Samples(int* c, double* mu)
{
	if (model)
		model->getMuGamma0Samples(c, mu);
}

void getMuTheta0Samples(int* c, double* mu)
{
	if (model)
		model->getMuTheta0Samples(c, mu);
}

void getTau2Gamma0Samples(int*c, double* tau2)
{
	if (model)
		model->getTau2Gamma0Samples(c, tau2);
}

void getTau2Theta0Samples(int* c, double* tau2)
{
	if (model)
		model->getTau2Theta0Samples(c, tau2);
}

void getThetaAccept(int* c, int* b, int* j, double* theta_acc)
{
	if (model)
		model->getThetaAccept(c, b ,j, theta_acc);
}

void getGammaAccept(int* c, int* b, int* j, double* gamma_acc)
{
	if (model)
		model->getGammaAccept(c, b ,j, gamma_acc);
}

void getThetaZeroAccept(int* c, int* b, int* j, double* zero_prop, double* zero_acc)
{
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		modelBB->getThetaZeroAccept(c, b ,j, zero_prop, zero_acc);
}

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
					SEXP mh_weight,
					SEXP pm_weights,
					SEXP adapt_min_w,
					SEXP adapt_chains,
					SEXP adapt_burnin,
					SEXP adapt_iter)
{

	try {

		if (model) {
			delete model;
			model = NULL;
		}

		model = new c212BB(sChains, sBurnin, sIter,  sNumBodySys, sMaxAEs, pNAE,
					pX, pY,
					pNC, pNT,
					ptheta, pgamma,
                    pmu_gamma_0_0,
                    ptau2_gamma_0_0,
                    pmu_theta_0_0,
                    ptau2_theta_0_0,
                    palpha_gamma_0_0,
                    pbeta_gamma_0_0,
                    palpha_theta_0_0,
                    pbeta_theta_0_0,
                    palpha_gamma,
                    pbeta_gamma,
                    palpha_theta,
                    pbeta_theta,
                    pmu_gamma_0,
					ptau2_gamma_0,
                    pmu_theta_0,
                    ptau2_theta_0,
                    pmu_gamma,
                    pmu_theta,
                    psigma2_gamma,
                    psigma2_theta,
					palpha_pi,
					pbeta_pi,
					plambda_alpha,
					plambda_beta,
					pPi,
					palgo,
					padapt_phase,
					sim_type,
					global_sim_params,
					sim_params,
					mh_weight,
					pm_weights,
					adapt_min_w,
					adapt_chains,
					adapt_burnin,
					adapt_iter);

		model->gibbs_sampler();
	}
	catch(...) {
		Rprintf("Unexpected exception c212BB_exec\n");
	}

	return R_NilValue;
}

SEXP getPiSamplesAll()
{
	SEXP samples = R_NilValue;
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		samples = modelBB->getPiSamples();

	return samples;
}

SEXP getAlphaPiSamplesAll()
{
	SEXP samples = R_NilValue;
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		samples = modelBB->getAlphaPiSamples();

	return samples;
}

SEXP getBetaPiSamplesAll()
{
	SEXP samples = R_NilValue;
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		samples = modelBB->getBetaPiSamples();

	return samples;
}

SEXP getAlphaPiAcceptAll()
{
	SEXP acc = R_NilValue;
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		acc = modelBB->getAlphaPiAccept();

	return acc;
}

SEXP getBetaPiAcceptAll()
{
	SEXP acc = R_NilValue;
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		acc = modelBB->getBetaPiAccept();

	return acc;
}

void getPiSamples(int*c, int* b, double* pi)
{
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		modelBB->getPiSamples(c, b, pi);
}

void getAlphaPiSamples(int* c, double* alpha_pi)
{
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		modelBB->getAlphaPiSamples(c, alpha_pi);
}

void getBetaPiSamples(int* c, double* beta_pi)
{
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		modelBB->getBetaPiSamples(c, beta_pi);
}

void getAlphaPiAccept(int* c, double* acc)
{
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		modelBB->getAlphaPiAccept(c, acc);
}

void getBetaPiAccept(int* c, double* acc)
{
	c212BB* modelBB = dynamic_cast<c212BB *>(model);

	if (modelBB)
		modelBB->getBetaPiAccept(c, acc);
}

SEXP c2121a_interim_hier2_exec(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_cntrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sLevel,
					SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta)
{
	try {

		if (model_interim) {
			delete model_interim;
			model_interim = NULL;
		}

		int level = *(INTEGER(sLevel));

		switch (level) {

			case 0:
				model_interim = new c2121a_poisson_mc_hier2_lev0(sChains, sBurnin,
										sIter,
										sSim_Type, sMem_Model, sGlobal_Sim_Param,
										sGlobal_Sim_Param_cntrl,
										sSim_Param,
										sMonitor, sNumIntervals, sMaxBs,
										sNumBodySys,
										sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,	
										pmu_gamma_0, ptau2_gamma_0,
										pmu_theta_0, ptau2_theta_0,
										palpha_gamma,
										pbeta_gamma,
										palpha_theta, pbeta_theta,
										pmu_gamma, pmu_theta, psigma2_gamma,
										psigma2_theta);
			break;

			case 1:
				model_interim = new c2121a_poisson_mc_hier2_lev1(sChains, sBurnin,
										sIter,
										sSim_Type, sMem_Model, sGlobal_Sim_Param,
										sGlobal_Sim_Param_cntrl,
										sSim_Param,
										sMonitor, sNumIntervals, sMaxBs,
										sNumBodySys,
										sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,
										pmu_gamma_0, ptau2_gamma_0,
										pmu_theta_0, ptau2_theta_0,
										palpha_gamma, pbeta_gamma,
										palpha_theta, pbeta_theta,
										pmu_gamma, pmu_theta, psigma2_gamma,
										psigma2_theta);
			break;

			default: // level 0
				model_interim = new c2121a_poisson_mc_hier2_lev0(sChains, sBurnin,
										sIter,
										sSim_Type, sMem_Model, sGlobal_Sim_Param,
										sGlobal_Sim_Param_cntrl,
										sSim_Param,
										sMonitor, sNumIntervals, sMaxBs,
										sNumBodySys,
										sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,
										pmu_gamma_0, ptau2_gamma_0,
										pmu_theta_0, ptau2_theta_0,
										palpha_gamma, pbeta_gamma,
										palpha_theta, pbeta_theta,
										pmu_gamma, pmu_theta, psigma2_gamma,
										psigma2_theta);
			break;
		}

		model_interim->gibbs_sampler();
		//model_hier2->gibbs_sampler();
	}
	catch(...) {
		Rprintf("Unexpected exception c2121a_poisson_mc\n");
	}
	return R_NilValue;
}

// Interim Analysis Models
SEXP c2121a_poisson_mc_exec(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,	
					SEXP sGlobal_Sim_Param,
					SEXP sGlobal_Sim_Param_cntrl,
					SEXP sSim_Param,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sLevel,
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
	try {

		if (model_interim) {
			delete model_interim;
			model_interim = NULL;
		}

		int level = *(INTEGER(sLevel));

		switch (level) {

			case 0:
				model_interim = new c2121a_poisson_mc_hier3_lev0(sChains, sBurnin,
										sIter,
										sSim_Type,
										sMem_Model,
										sGlobal_Sim_Param,
										sGlobal_Sim_Param_cntrl,
										sSim_Param,
										sMonitor,
										sNumIntervals, sMaxBs,
										sNumBodySys,
										sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,	
										pmu_gamma_0_0, ptau2_gamma_0_0,
										pmu_theta_0_0, ptau2_theta_0_0,
										palpha_gamma_0_0,
										pbeta_gamma_0_0,
										palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma,
										pbeta_gamma, palpha_theta,
										pbeta_theta, pmu_gamma_0, ptau2_gamma_0,
										pmu_theta_0,
										ptau2_theta_0,
										pmu_gamma, pmu_theta, psigma2_gamma,
										psigma2_theta);
			break;

			case 1:
				model_interim = new c2121a_poisson_mc_hier3_lev1(sChains, sBurnin, sIter,
										sSim_Type,
										sMem_Model, 
										sGlobal_Sim_Param,
										sGlobal_Sim_Param_cntrl,
										sSim_Param,
										sMonitor,
										sNumIntervals, sMaxBs,
										sNumBodySys,
										sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,
										pmu_gamma_0_0, ptau2_gamma_0_0,
										pmu_theta_0_0, ptau2_theta_0_0,
										palpha_gamma_0_0, pbeta_gamma_0_0,
										palpha_theta_0_0, pbeta_theta_0_0,
										palpha_gamma, pbeta_gamma, palpha_theta,
										pbeta_theta, pmu_gamma_0, ptau2_gamma_0,
										pmu_theta_0, ptau2_theta_0,
										pmu_gamma, pmu_theta, psigma2_gamma,
										psigma2_theta);
			break;

			case 2:
				model_interim = new c2121a_poisson_mc_hier3_lev2(sChains, sBurnin, sIter,
										sSim_Type,
										sMem_Model,
										sGlobal_Sim_Param,
										sGlobal_Sim_Param_cntrl,
										sSim_Param,
										sMonitor,
										sNumIntervals, sMaxBs,
										sNumBodySys,
										sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,
										pmu_gamma_0_0, ptau2_gamma_0_0,
										pmu_theta_0_0, ptau2_theta_0_0,
										palpha_gamma_0_0, pbeta_gamma_0_0,
										palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma,
										pbeta_gamma, palpha_theta,
										pbeta_theta, pmu_gamma_0, ptau2_gamma_0,
										pmu_theta_0, ptau2_theta_0,
										pmu_gamma, pmu_theta, psigma2_gamma,
										psigma2_theta);
			break;

			default: // level 0
				model_interim = new c2121a_poisson_mc_hier3_lev0(sChains, sBurnin,
										sIter,
										sSim_Type,
										sMem_Model,
										sGlobal_Sim_Param,
										sGlobal_Sim_Param_cntrl,
										sSim_Param,
										sMonitor,
										sNumIntervals, sMaxBs,
										sNumBodySys,
										sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,
										pmu_gamma_0_0, ptau2_gamma_0_0,
										pmu_theta_0_0, ptau2_theta_0_0,
										palpha_gamma_0_0, pbeta_gamma_0_0,
										palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma,
										pbeta_gamma, palpha_theta,
										pbeta_theta, pmu_gamma_0, ptau2_gamma_0,
										pmu_theta_0, ptau2_theta_0,
										pmu_gamma, pmu_theta, psigma2_gamma,
										psigma2_theta);
			break;
		}

		model_interim->gibbs_sampler();
	}
	catch(...) {
		Rprintf("Unexpected exception c2121a_poisson_mc\n");
	}
	return R_NilValue;
}

SEXP getThetaSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	if (model_interim)
		samples = model_interim->getThetaSamples();

	return samples;
}

SEXP getGammaSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	if (model_interim)
		samples = model_interim->getGammaSamples();

	return samples;
}

SEXP getMuThetaSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	if (model_interim)
		samples = model_interim->getMuThetaSamples();

	return samples;
}

SEXP getMuGammaSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	if (model_interim)
		samples = model_interim->getMuGammaSamples();

	return samples;
}

SEXP getSigma2GammaSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	if (model_interim)
		samples = model_interim->getSigma2GammaSamples();

	return samples;
}

SEXP getSigma2ThetaSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	if (model_interim)
		samples = model_interim->getSigma2ThetaSamples();

	return samples;
}

SEXP getGammaAcceptInterimAll()
{
	SEXP acc = R_NilValue;

	if (model_interim)
		acc = model_interim->getGammaAccept();

	return acc;
}

SEXP getThetaAcceptInterimAll()
{
	SEXP acc = R_NilValue;

	if (model_interim)
		acc = model_interim->getThetaAccept();

	return acc;
}

SEXP getMuGamma0SamplesInterimAll()
{
	SEXP samples = R_NilValue;

	c2121a_poisson_mc_hier3_lev0* model1a = dynamic_cast<c2121a_poisson_mc_hier3_lev0 *>(model_interim);

	if (model1a)
		samples = model1a->getMuGamma0Samples();

	return samples;
}

SEXP getMuTheta0SamplesInterimAll()
{
	SEXP samples = R_NilValue;

	c2121a_poisson_mc_hier3_lev0* model1a = dynamic_cast<c2121a_poisson_mc_hier3_lev0 *>(model_interim);

	if (model1a)
		samples = model1a->getMuTheta0Samples();

	return samples;
}

SEXP getTau2Gamma0SamplesInterimAll()
{
	SEXP samples = R_NilValue;

	c2121a_poisson_mc_hier3_lev0* model1a = dynamic_cast<c2121a_poisson_mc_hier3_lev0 *>(model_interim);

	if (model1a)
		samples = model1a->getTau2Gamma0Samples();

	return samples;
}

SEXP getTau2Theta0SamplesInterimAll()
{
	SEXP samples = R_NilValue;

	c2121a_poisson_mc_hier3_lev0* model1a = dynamic_cast<c2121a_poisson_mc_hier3_lev0 *>(model_interim);

	if (model1a)
		samples = model1a->getTau2Theta0Samples();

	return samples;
}

void getThetaSamplesInterim(int *c, int* l, int* b, int* j, double* theta_samples)
{
	if (model_interim)
		model_interim->getThetaSamples(c, l, b ,j, theta_samples);
}

void getGammaSamplesInterim(int *c, int* l, int* b, int* j, double* gamma_samples)
{
	if (model_interim)
		model_interim->getGammaSamples(c, l, b, j, gamma_samples);
}

void getMuThetaSamplesInterim(int *c, int* l, int* b, double* mu_theta)
{
	if (model_interim)
		model_interim->getMuThetaSamples(c, l, b, mu_theta);
}

void getMuGammaSamplesInterim(int *c, int* l, int* b, double* mu_gamma)
{
	if (model_interim)
		model_interim->getMuGammaSamples(c, l, b, mu_gamma);
}

void getSigma2ThetaSamplesInterim(int *c, int* l, int* b, double* sigma2)
{
	if (model_interim)
		model_interim->getSigma2ThetaSamples(c, l, b, sigma2);
}

void getSigma2GammaSamplesInterim(int *c, int* l, int* b, double* sigma2)
{
	if (model_interim)
		model_interim->getSigma2GammaSamples(c, l, b, sigma2);
}

void getMuGamma0SamplesInterim(int *c, int* l, double* mu)
{
	c2121a_poisson_mc_hier3_lev0* model1a = dynamic_cast<c2121a_poisson_mc_hier3_lev0 *>(model_interim);

	if (model1a)
		model1a->getMuGamma0Samples(c, l, mu);
}

void getMuTheta0SamplesInterim(int *c, int* l, double* mu)
{
	c2121a_poisson_mc_hier3_lev0* model1a = dynamic_cast<c2121a_poisson_mc_hier3_lev0 *>(model_interim);

	if (model1a)
		model1a->getMuTheta0Samples(c, l, mu);
}

void getTau2Gamma0SamplesInterim(int *c, int* l, double* tau2)
{
	c2121a_poisson_mc_hier3_lev0* model1a = dynamic_cast<c2121a_poisson_mc_hier3_lev0 *>(model_interim);

	if (model1a)
		model1a->getTau2Gamma0Samples(c, l, tau2);
}

void getTau2Theta0SamplesInterim(int *c, int* l, double* tau2)
{
	c2121a_poisson_mc_hier3_lev0* model1a = dynamic_cast<c2121a_poisson_mc_hier3_lev0 *>(model_interim);

	if (model1a)
		model1a->getTau2Theta0Samples(c, l, tau2);
}

void getThetaAcceptInterim(int *c, int* l, int* b, int* j, double* theta_acc)
{
	if (model_interim)
		model_interim->getThetaAccept(c, l, b ,j, theta_acc);
}

void getGammaAcceptInterim(int *c, int* l, int* b, int* j, double* gamma_acc)
{
	if (model_interim)
		model_interim->getGammaAccept(c, l, b ,j, gamma_acc);
}

void Release_Interim()
{
	if (model_interim)
		delete model_interim;
	model_interim = NULL;
}

SEXP c212BB_poisson_mc_exec(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model, SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sLevel,
					SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
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
					SEXP padapt_phase)
{
	try {

		if (model_interim) {
			delete model_interim;
			model_interim = NULL;
		}

		int level = *(INTEGER(sLevel));

		switch(level) {
			case 0:
				model_interim = new c212BB_poisson_mc_hier3_lev0(sChains, sBurnin, sIter, sSim_Type,
								sMem_Model, sGlobal_Sim_Params,
								sSim_Params,
								MH_weight,
								pm_weights,
								sMonitor,
								sNumIntervals, sMaxBs, sNumBodySys,
								sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0,
								pmu_theta_0_0, ptau2_theta_0_0, palpha_gamma_0_0, pbeta_gamma_0_0,
								palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma, pbeta_gamma, palpha_theta,
								pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0, ptau2_theta_0,
								pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
								pPi,
								palpha_pi,
								pbeta_pi,
								plambda_alpha,
								plambda_beta,
								palgo,
								padapt_phase);
			break;

			case 1:
				model_interim = new c212BB_poisson_mc_hier3_lev1(sChains, sBurnin, sIter, sSim_Type,
								sMem_Model, sGlobal_Sim_Params,
								sSim_Params,
								MH_weight,
								pm_weights,
								sMonitor,
								sNumIntervals, sMaxBs, sNumBodySys,
								sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0,
								pmu_theta_0_0, ptau2_theta_0_0, palpha_gamma_0_0, pbeta_gamma_0_0,
								palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma, pbeta_gamma, palpha_theta,
								pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0, ptau2_theta_0,
								pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
								pPi,
								palpha_pi,
								pbeta_pi,
								plambda_alpha,
								plambda_beta,
								palgo,
								padapt_phase);
			break;

			case 2:
				model_interim = new c212BB_poisson_mc_hier3_lev2(sChains, sBurnin, sIter, sSim_Type,
								sMem_Model, sGlobal_Sim_Params,
								sSim_Params,
								MH_weight,
								pm_weights,
								sMonitor,
								sNumIntervals, sMaxBs, sNumBodySys,
								sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0,
								pmu_theta_0_0, ptau2_theta_0_0, palpha_gamma_0_0, pbeta_gamma_0_0,
								palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma, pbeta_gamma, palpha_theta,
								pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0, ptau2_theta_0,
								pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
								pPi,
								palpha_pi,
								pbeta_pi,
								plambda_alpha,
								plambda_beta,
								palgo,
								padapt_phase);
			break;

			default:
				model_interim = new c212BB_poisson_mc_hier3_lev0(sChains, sBurnin, sIter, sSim_Type,
								sMem_Model,
								sGlobal_Sim_Params,
								sSim_Params,
								MH_weight,
								pm_weights,
								sMonitor,
								sNumIntervals, sMaxBs, sNumBodySys,
								sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma, pmu_gamma_0_0, ptau2_gamma_0_0,
								pmu_theta_0_0, ptau2_theta_0_0, palpha_gamma_0_0, pbeta_gamma_0_0,
								palpha_theta_0_0, pbeta_theta_0_0, palpha_gamma, pbeta_gamma, palpha_theta,
								pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0, ptau2_theta_0,
								pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
								pPi,
								palpha_pi,
								pbeta_pi,
								plambda_alpha,
								plambda_beta,
								palgo,
								padapt_phase);
			break;
		}

		model_interim->gibbs_sampler();
	}
	catch(...) {
		Rprintf("Unexpected exception c2121a_poisson_mc\n");
	}
	return R_NilValue;
}

SEXP c212BB_interim_hier2_exec(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobalSim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP sMonitor,
					SEXP sNumIntervals, SEXP sLevel,
					SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi,
					SEXP palpha_pi,
					SEXP pbeta_pi,
					SEXP palgo,
					SEXP padapt_phase)
{
	try {

		if (model_interim) {
			delete model_interim;
			model_interim = NULL;
		}

		int level = *(INTEGER(sLevel));

		switch(level) {
			case 0:
				model_interim = new c212BB_poisson_mc_hier2_lev0(sChains, sBurnin, sIter, sSim_Type,
								sMem_Model,
								sGlobalSim_Params,
								sSim_Params,
								MH_weight,
								pm_weights,
								sMonitor,
								sNumIntervals, sMaxBs, sNumBodySys,
								sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,
								palpha_gamma, pbeta_gamma, palpha_theta,
								pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0, ptau2_theta_0,
								pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
								pPi,
								palpha_pi,
								pbeta_pi,
								palgo,
								padapt_phase);
			break;

			case 1:
				model_interim = new c212BB_poisson_mc_hier2_lev1(sChains, sBurnin, sIter, sSim_Type,
								sMem_Model,
								sGlobalSim_Params,
								sSim_Params,
								MH_weight,
								pm_weights,
								sMonitor,
								sNumIntervals, sMaxBs, sNumBodySys,
								sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,
								palpha_gamma, pbeta_gamma, palpha_theta,
								pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0, ptau2_theta_0,
								pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
								pPi,
								palpha_pi,
								pbeta_pi,
								palgo,
								padapt_phase);
			break;

			default:
				model_interim = new c212BB_poisson_mc_hier2_lev0(sChains, sBurnin, sIter, sSim_Type,
								sMem_Model,
								sGlobalSim_Params,
								sSim_Params,
								MH_weight,
								pm_weights,
								sMonitor,
								sNumIntervals, sMaxBs, sNumBodySys,
								sMaxAEs, sNAE, pX, pY, pC, pT, ptheta, pgamma,
								palpha_gamma, pbeta_gamma, palpha_theta,
								pbeta_theta, pmu_gamma_0, ptau2_gamma_0, pmu_theta_0, ptau2_theta_0,
								pmu_gamma, pmu_theta, psigma2_gamma, psigma2_theta,
								pPi,
								palpha_pi,
								pbeta_pi,
								palgo,
								padapt_phase);
			break;
		}

		model_interim->gibbs_sampler();
	}
	catch(...) {
		Rprintf("Unexpected exception c2121a_poisson_mc\n");
	}
	return R_NilValue;
}

SEXP getPiSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	c212BB_poisson_mc_hier2_lev0* modelBB_h2 =
				dynamic_cast<c212BB_poisson_mc_hier2_lev0 *>(model_interim);
	c212BB_poisson_mc_hier3_lev0* modelBB_h3 = NULL;

	if (modelBB_h2)
		samples = modelBB_h2->getPiSamples();
	else {
		modelBB_h3 = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);
		if (modelBB_h3)
			samples = modelBB_h3->getPiSamples();
	}

	return samples;
}

SEXP getAlphaPiSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	c212BB_poisson_mc_hier3_lev0* modelBB = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);

	if (modelBB)
		samples = modelBB->getAlphaPiSamples();

	return samples;
}

SEXP getBetaPiSamplesInterimAll()
{
	SEXP samples = R_NilValue;

	c212BB_poisson_mc_hier3_lev0* modelBB = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);

	if (modelBB)
		samples = modelBB->getBetaPiSamples();

	return samples;
}

SEXP getAlphaPiAcceptInterimAll()
{
	SEXP acc = R_NilValue;

	c212BB_poisson_mc_hier3_lev0* modelBB = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);

	if (modelBB)
		acc = modelBB->getAlphaPiAccept();

	return acc;
}

SEXP getBetaPiAcceptInterimAll()
{
	SEXP acc = R_NilValue;

	c212BB_poisson_mc_hier3_lev0* modelBB = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);

	if (modelBB)
		acc = modelBB->getBetaPiAccept();

	return acc;
}

void getPiSamplesInterim(int *c, int* l, int* b, double* pi)
{
	c212BB_poisson_mc_hier2_lev0* modelBB_h2 =
				dynamic_cast<c212BB_poisson_mc_hier2_lev0 *>(model_interim);
	c212BB_poisson_mc_hier3_lev0* modelBB_h3 = NULL;

	if (modelBB_h2)
		modelBB_h2->getPiSamples(c, l, b, pi);
	else {
		modelBB_h3 = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);
		if (modelBB_h3)
			modelBB_h3->getPiSamples(c, l, b, pi);
	}
}

void getAlphaPiSamplesInterim(int *c, int* l, double* alpha_pi)
{
	c212BB_poisson_mc_hier3_lev0* modelBB = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);

	if (modelBB)
		modelBB->getAlphaPiSamples(c, l, alpha_pi);
}

void getBetaPiSamplesInterim(int *c, int* l, double* beta_pi)
{
	c212BB_poisson_mc_hier3_lev0* modelBB = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);

	if (modelBB)
		modelBB->getBetaPiSamples(c, l, beta_pi);
}

void getAlphaPiAcceptInterim(int *c, int* l, double* acc)
{
	c212BB_poisson_mc_hier3_lev0* modelBB = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);

	if (modelBB)
		modelBB->getAlphaPiAccept(c, l, acc);
}

void getBetaPiAcceptInterim(int *c, int* l, double* acc)
{
	c212BB_poisson_mc_hier3_lev0* modelBB = dynamic_cast<c212BB_poisson_mc_hier3_lev0 *>(model_interim);

	if (modelBB)
		modelBB->getBetaPiAccept(c, l, acc);
}
