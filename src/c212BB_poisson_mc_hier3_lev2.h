#ifndef C212BB_POISSON_MC_HIER3_LEV2_H
#define C212BB_POISSON_MC_HIER3_LEV2_H

class c212BB_poisson_mc_hier3_lev2 : public c212BB_poisson_mc_hier3_lev0 {
	public:
		c212BB_poisson_mc_hier3_lev2();

		c212BB_poisson_mc_hier3_lev2(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sim_type,
					SEXP sMem_Model,
					SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
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
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi, SEXP plambda_alpha,
					SEXP plambda_beta,
					SEXP palgo, SEXP padapt_phase);

		virtual ~c212BB_poisson_mc_hier3_lev2();

		virtual void gibbs_sampler();

	protected:
		virtual void release();
		virtual void simulate_MH();
		virtual void simulate_SLICE();

		virtual void sample_mu_gamma_0(int burnin, int iter);
		virtual void sample_mu_theta_0(int burnin, int iter);
		virtual void sample_tau2_gamma_0(int burnin, int iter);
		virtual void sample_tau2_theta_0(int burnin, int iter);

		virtual void sample_alpha_pi_MH(int gBurnin, int i);
		virtual void sample_alpha_pi_SLICE(int gBurnin, int i);
		double log_f_alpha_pi(int c, double alpha);

		virtual void sample_beta_pi_SLICE(int gBurnin, int i);
		virtual void sample_beta_pi_MH(int gBurnin, int i);
		virtual double log_f_beta_pi(int c, double alpha);

		virtual void sample_pi(int gBurnin, int i);


		virtual void sample_mu_gamma(int burnin, int iter);
		virtual void sample_mu_theta(int burnin, int iter);
		virtual void sample_sigma2_gamma(int burnin, int iter);
		virtual void sample_sigma2_theta(int burnin, int iter);
		double cMIN(double a, double b);

		virtual void init(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sSim_Type,
					SEXP sMem_Model, SEXP sGlobal_Sim_Params,
					SEXP sSim_Params,
					SEXP MH_weight,
					SEXP pm_weights,
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
					SEXP pmu_theta, SEXP psigma2_gamma, SEXP psigma2_theta,
					SEXP pPi, SEXP palpha_pi, SEXP pbeta_pi, SEXP plambda_alpha,
					SEXP plambda_beta,
					SEXP palgo, SEXP padapt_phase);

		virtual SEXP getL3Samples(double** &data);
		virtual SEXP getL3Accept(int* &data);

	public:
		virtual SEXP getMuGamma0Samples();
		virtual SEXP getMuTheta0Samples();
		virtual SEXP getTau2Gamma0Samples();
		virtual SEXP getTau2Theta0Samples();
		virtual SEXP getAlphaPiSamples();
		virtual SEXP getBetaPiSamples();
		virtual SEXP getAlphaPiAccept();
		virtual SEXP getBetaPiAccept();

		virtual void getMuGamma0Samples(int *c, int *l, double* mu);
		virtual void getMuTheta0Samples(int *c, int *l, double* mu);
		virtual void getTau2Gamma0Samples(int *c, int *l, double* tau2);
		virtual void getTau2Theta0Samples(int *c, int *l, double* tau2);

		virtual void getAlphaPiSamples(int *c, int*l, double* alpha_pi);
		virtual void getBetaPiSamples(int *c, int*l, double* beta_pi);
		virtual void getAlphaPiAccept(int *c, int*l, double* acc);
		virtual void getBetaPiAccept(int *c, int*l, double* acc);

		double *mu_theta_0;   // Current value of the sampled distribution - updated constantly
		double *mu_gamma_0;  // Current value of the sampled distribution - updated constantly
		double *tau2_theta_0;   // Current value of the sampled distribution - updated constantly
		double *tau2_gamma_0;   // Current value of the sampled distribution - updated constantly

		double** mu_theta_0_samples;
		double** mu_gamma_0_samples;
		double** tau2_theta_0_samples;
		double** tau2_gamma_0_samples;

		double* alpha_pi;  // Current value of the sampled distribution - updated constantly
		double** alpha_pi_samples;
		int* alpha_pi_acc;

		double* beta_pi;  // Current value of the sampled distribution - updated constantly
		double** beta_pi_samples;
		int* beta_pi_acc;

	protected:

		// Adapt for some of the MH steps
		int in_apapt_phase;
};

#endif


