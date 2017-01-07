#ifndef C2121A_H
#define C2121A_H

class c2121a {
	public:
		c2121a();
		c2121a(SEXP pChains, SEXP pBurnin, SEXP pIter, SEXP pNumBodySys, SEXP pMaxAEs,
					SEXP pNAE, SEXP sim_type, SEXP pGlobal_Sim_Param,
					SEXP pGlobal_Sim_Param_Cntrl,
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

		virtual ~c2121a();

		virtual void gibbs_sampler();

	protected:
		virtual void release();
		virtual void simulate_MH();
		virtual void simulate_SLICE();

		virtual void sample_mu_gamma_0(int c,int burnin, int iter);
		virtual void sample_mu_theta_0(int c,int burnin, int iter);
		virtual void sample_tau2_gamma_0(int c,int burnin, int iter);
		virtual void sample_tau2_theta_0(int c,int burnin, int iter);
		virtual void sample_mu_gamma(int c,int burnin, int iter);
		virtual void sample_mu_theta(int c,int burnin, int iter);
		virtual void sample_sigma2_gamma(int c,int burnin, int iter);
		virtual void sample_sigma2_theta(int c,int burnin, int iter);
		virtual double log_f_gamma(int c, int b, int j, double gamm);
		virtual void sample_gamma_MH(int c,int burnin, int iter);
		virtual void sample_gamma_SLICE(int c,int burnin, int iter);
		virtual double log_f_theta(int c, int b, int j, double theta);
		virtual void sample_theta_MH(int c,int burnin, int iter);
		virtual void sample_theta_SLICE(int c,int burnin, int iter);
		virtual double cMIN(double a, double b);

		virtual void init(SEXP Chains, SEXP pBurnin, SEXP pIter, SEXP pNumBodySys,
							SEXP pMaxAEs, SEXP pNAE, 
							SEXP pSim_Type,
							SEXP pGlobal_Sim_Param, SEXP pGlobal_Sim_Param_cntrl,
							SEXP pSim_Param,
							SEXP pX,
							SEXP pY, SEXP pNC, SEXP pNT,
							SEXP ptheta, SEXP pgamma, SEXP pmu_gamma_0_0,
							SEXP ptau2_gamma_0_0,
							SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0,
							SEXP palpha_gamma_0_0,
							SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0,
							SEXP pbeta_theta_0_0,
							SEXP palpha_gamma, SEXP pbeta_gamma,
							SEXP palpha_theta,
							SEXP pbeta_theta, SEXP pmu_gamma_0,
							SEXP ptau2_gamma_0,
							SEXP pmu_theta_0, SEXP ptau2_theta_0,
							SEXP pmu_gamma, SEXP pmu_theta,
							SEXP psigma2_gamma, SEXP psigma2_theta);

		virtual void initSimParams(SEXP sim_params);

		virtual SEXP getL1Samples(double**** &samples);
		virtual SEXP getL2Samples(double*** &samples);
		virtual SEXP getL3Samples(double** &data);

		virtual SEXP getL1Accept(int*** &samples);

	public:
		virtual SEXP getThetaSamples();
		virtual SEXP getGammaSamples();
		virtual SEXP getMuThetaSamples();
		virtual SEXP getMuGammaSamples();
		virtual SEXP getSigma2ThetaSamples();
		virtual SEXP getSigma2GammaSamples();
		virtual SEXP getMuGamma0Samples();
		virtual SEXP getMuTheta0Samples();
		virtual SEXP getTau2Gamma0Samples();
		virtual SEXP getTau2Theta0Samples();
		virtual SEXP getThetaAccept();
		virtual SEXP getGammaAccept();

		virtual void getThetaSamples(int* c, int* b, int* j, double* theta_samples);
		virtual double* getThetaSamples(int c, int b, int j);
		virtual void getGammaSamples(int* c, int* b, int* j, double* gamma_samples);
		virtual void getMuThetaSamples(int* c, int* b, double* mu_theta);
		virtual void getMuGammaSamples(int* c, int* b, double* mu_gamma);
		virtual void getSigma2ThetaSamples(int* c, int* b, double* sigma2);
		virtual void getSigma2GammaSamples(int* c, int* b, double* sigma2);
		virtual void getMuGamma0Samples(int* c, double* mu);
		virtual void getMuTheta0Samples(int* c, double* mu);
		virtual void getTau2Gamma0Samples(int* c, double* tau2);
		virtual void getTau2Theta0Samples(int* c, double* tau2);
		virtual void getThetaAccept(int* c, int* b, int* j, double* theta_acc);
		virtual void getGammaAccept(int* c, int* b, int* j, double* gamma_acc);

	protected:
		int gChains;
		int gBurnin;
		int gIter;
		char* sim_type;
		int* gNAE;
		int gNumBodySys;
		int gMaxAEs;
		// Global simulation parameter values
		double gSim_Param;
		double gSim_Param_cntrl;

		// Individual simulation parameter values
		double** gW_gamma;
		double** gW_theta;
		int** gW_gamma_control;
		int** gW_theta_control;
		double** gSigma_MH_gamma;
		double** gSigma_MH_theta;

		static const char* sColType;
		static const char* sColVariable;
		static const char* sColParam;
		static const char* sColValue;
		static const char* sColControl;
		static const char* sColB;
		static const char* sColj;

		static const char* sParam_w;
		static const char* sParam_sigma_MH;
		static const char* sVariable_gamma;
		static const char* sVariable_theta;

		// Hyper-parameters
		double mu_theta_0_0;   // Fixed hyper-parameter value
		double mu_gamma_0_0; // Fixed hyper-parameter value
		double tau2_theta_0_0; // Fixed hyper-parameter value
		double tau2_gamma_0_0; // Fixed hyper-parameter value
		double alpha_gamma_0_0; // Fixed hype-parameter value
		double beta_gamma_0_0; // Fixed hype-parameter value
		double alpha_theta_0_0; // Fixed hype-parameter value
		double beta_theta_0_0; // Fixed hype-parameter value
		double alpha_gamma; // Fixed hype-parameter value
		double beta_gamma; // Fixed hype-parameter value
		double alpha_theta; // Fixed hype-parameter value
		double beta_theta; // Fixed hype-parameter value

		double *mu_theta_0;   // Current value of the sampled distribution - updated constantly
		double *mu_gamma_0;  // Current value of the sampled distribution - updated constantly
		double *tau2_theta_0;   // Current value of the sampled distribution - updated constantly
		double *tau2_gamma_0;   // Current value of the sampled distribution - updated constantly

		// All chains
		double** mu_theta;
		double** mu_gamma;
		double** sigma2_theta;
		double** sigma2_gamma;

		double*** gTheta;
		double*** gGamma;
		int*** gTheta_acc;
		int*** gGamma_acc;

		// Data values
		int **x;
		int **y;
		int **NC;
		int **NT;

		// Samples
		double**** gTheta_samples;
		double**** gGamma_samples;
		double** mu_theta_0_samples;
		double** mu_gamma_0_samples;
		double** tau2_theta_0_samples;
		double** tau2_gamma_0_samples;
		double*** mu_theta_samples;
		double*** mu_gamma_samples;
		double*** sigma2_theta_samples;
		double*** sigma2_gamma_samples;
};

#endif
