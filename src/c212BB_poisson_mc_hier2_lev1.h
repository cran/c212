#ifndef C212BB_POISSON_MC_HIER2_LEV1_H
#define C212BB_POISSON_MC_HIER2_LEV1_H

class c212BB_poisson_mc_hier2_lev1 : public c212BB_poisson_mc_hier2_lev0 {
	public:
		c212BB_poisson_mc_hier2_lev1();

		c212BB_poisson_mc_hier2_lev1(SEXP sChains, SEXP sBurnin, SEXP sIter,
					SEXP sim_type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP sMonitor,
					SEXP sNumIntervals,
					SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi,
					SEXP palgo, SEXP padapt_phase);

		virtual ~c212BB_poisson_mc_hier2_lev1();

		virtual void gibbs_sampler();

	protected:
		virtual void release();
		virtual void simulate_MH();
		virtual void simulate_SLICE();

		virtual void sample_pi(int gBurnin, int i);

		virtual double log_f_theta(int c, int i, int b, int j, double theta);
		virtual double log_q_theta(int l, int b, int j, double p, double theta, double mean);

		virtual void sample_mu_gamma(int burnin, int iter);
		virtual void sample_mu_theta(int burnin, int iter);
		virtual void sample_sigma2_gamma(int burnin, int iter);
		virtual void sample_sigma2_theta(int burnin, int iter);
		virtual double log_f_gamma(int c, int i, int b, int j, double gamm);
		virtual void sample_gamma_MH(int burnin, int iter);
		virtual void sample_gamma_SLICE(int burnin, int iter);
		virtual void sample_theta_MH(int burnin, int iter);
		double cMIN(double a, double b);

		virtual void init(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP sMonitor,
					SEXP sNumIntervals,
					SEXP sMaxBs, SEXP sNumBodySys, SEXP sMaxAEs, SEXP sNAE,
					SEXP pX, SEXP pY, SEXP pC, SEXP pT, SEXP ptheta, SEXP pgamma,
					SEXP palpha_gamma,
					SEXP pbeta_gamma, SEXP palpha_theta, SEXP pbeta_theta,
					SEXP pmu_gamma_0,
					SEXP ptau2_gamma_0, SEXP pmu_theta_0, SEXP ptau2_theta_0,
					SEXP pmu_gamma,
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi,
					SEXP palgo, SEXP padapt_phase);

		virtual SEXP getL2Samples(double*** &data);

	public:
		virtual SEXP getMuThetaSamples();
		virtual SEXP getMuGammaSamples();
		virtual SEXP getSigma2ThetaSamples();
		virtual SEXP getSigma2GammaSamples();
		virtual SEXP getPiSamples();

		virtual void getMuThetaSamples(int *c, int *l, int* b, double* mu_theta);
		virtual void getMuGammaSamples(int *c, int *l, int* b, double* mu_gamma);
		virtual void getSigma2ThetaSamples(int *c, int *l, int* b, double* sigma2);
		virtual void getSigma2GammaSamples(int *c, int *l, int* b, double* sigma2);
		virtual void getPiSamples(int *c, int *l, int* b, double* pi);

	protected:

		double** gPi;
		double*** gPi_samples;

		double** mu_theta;
		double** mu_gamma;
		double** sigma2_theta;
		double** sigma2_gamma;

		double*** mu_theta_samples;
		double*** mu_gamma_samples;
		double*** sigma2_theta_samples;
		double*** sigma2_gamma_samples;

		// Adapt for some of the MH steps
		int in_apapt_phase;
};

#endif


