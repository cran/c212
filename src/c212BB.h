#ifndef C212BB_H
#define C212BB_H

struct gSimParams {
	double value;
	double control;
	gSimParams(double v, double c) { value = v; control = c; }
	gSimParams() { value = 0.0; control = 0.0; }
};

struct MHAdaptParams {
	double w_min;
	int chains;
	int burnin;
	int iter;
};

class c212BB : public c2121a {

	public:
		c212BB(SEXP sChains, SEXP sBurnin, SEXP sIter, SEXP sNumBodySys, SEXP sMaxAEs,
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
					bool on_screen = TRUE);

		c212BB(int sChains, int sBurnin, int sIter, int sNumBodySys, int sMaxAEs,
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
					std::map<const char*, gSimParams>& mSimParams,
					double** gamma_params,
					int** gamma_ctrl,
					double** theta_params,
					double MH_weight,
					double** pm_weights,
					MHAdaptParams& adapt_params,
					bool on_screen = TRUE);
		~c212BB();

	protected:
		void simulate_MH();
		void simulate_SLICE();
		
		double log_f_theta(int c, int b, int j, double theta);
		//void sample_theta_BB2004(int c, int burnin, int iter);
		double log_q_theta(int c, int b, int j, double p, double theta, double mean);
		void sample_theta_MH(int c, int burnin, int iter);
		void sample_theta_MIS_Adapt(int c, int burnin, int iter);
		void sample_theta_Independent_MH(int c, int burnin, int iter);
		
		double log_f_beta_pi(int c, double beta);
		void sample_alpha_pi_MH(int c, int burnin, int iter);
		void sample_alpha_pi_SLICE(int c, int burnin, int iter);
		double log_f_alpha_pi(int c, double alpha);
		void sample_beta_pi_MH(int c, int burnin, int iter);
		void sample_beta_pi_SLICE(int c, int burnin, int iter);
		void sample_pi(int c, int burnin, int iter);

		virtual void sample_mu_theta(int c,int burnin, int iter);
		virtual void sample_sigma2_theta(int c,int burnin, int iter);

		void init(SEXP pChains, SEXP pBurnin, SEXP pIter, SEXP pNumBodySys,
					SEXP pMaxAEs, SEXP pNAE, 
					SEXP pX, SEXP pY, SEXP pNC, SEXP pNT,
					SEXP ptheta, SEXP pgamma, SEXP pmu_gamma_0_0, SEXP ptau2_gamma_0_0,
					SEXP pmu_theta_0_0, SEXP ptau2_theta_0_0, SEXP palpha_gamma_0_0,
                    SEXP pbeta_gamma_0_0, SEXP palpha_theta_0_0, SEXP pbeta_theta_0_0,
					SEXP palpha_gamma, SEXP pbeta_gamma, SEXP palpha_theta,
					SEXP pbeta_theta, SEXP pmu_gamma_0, SEXP ptau2_gamma_0,
					SEXP pmu_theta_0, SEXP ptau2_theta_0, SEXP pmu_gamma,
					SEXP pmu_theta,
                    SEXP psigma2_gamma, SEXP psigma2_theta,
                    SEXP palpha_pi, SEXP pbeta_pi, SEXP plambda_alpha,
					SEXP plambda_beta, SEXP pPi, SEXP padapt_phase, SEXP sim_type,
					SEXP global_sim_params,
					SEXP sim_params,
					SEXP MH_weight,
					SEXP pm_weights,
					SEXP adapt_min_w,
					SEXP adapt_chains, SEXP adapt_burnin, SEXP adapt_iter);

		void init(int pChains, int pBurnin, int pIter, int pNumBodySys,
					int pMaxAEs, int* pNAE, 
					int** pX, int** pY, int** pNC, int** pNT,
					double*** ptheta, double*** pgamma, double pmu_gamma_0_0,
					double ptau2_gamma_0_0,
					double pmu_theta_0_0, double ptau2_theta_0_0,
					double palpha_gamma_0_0,
                    double pbeta_gamma_0_0, double palpha_theta_0_0,
					double pbeta_theta_0_0,
					double palpha_gamma, double pbeta_gamma, double palpha_theta,
					double pbeta_theta, double* pmu_gamma_0, double* ptau2_gamma_0,
					double* pmu_theta_0, double* ptau2_theta_0, double** pmu_gamma,
					double** pmu_theta,
                    double** psigma2_gamma, double** psigma2_theta,
                    double* palpha_pi, double* pbeta_pi, double plambda_alpha,
					double plambda_beta, double** pPi, int padapt_phase,
					const char* sim_type,
					std::map<const char*, gSimParams>& mSimParams,
					double** gamma_params,
					int** gamma_ctrl,
					double** theta_params,
					double MH_weight,
					double** pm_weights,
					MHAdaptParams& adapt_params);

		void initGlobalSimParams(SEXP sim_type, SEXP sim_params);
		void initGlobalSimParams(const char* sim_type,
									std::map<const char*, gSimParams>& mSimParams);

		void initSimParams(SEXP sim_params);
		void initSimParams(double** gamma_params, int** gamma_cntrl,
													double** theta_params);

		void initPMWeights(SEXP pm_weights);
		void initPMWeights(double** pm_weights);

		void adaptPhaseMH();

		void release();

		double log_g(double val, int b, int j);

		virtual SEXP getL3Accept(int* &data);

	public:
		void gibbs_sampler();

		virtual SEXP getPiSamples();
		virtual SEXP getAlphaPiSamples();
		virtual SEXP getBetaPiSamples();
		virtual SEXP getAlphaPiAccept();
		virtual SEXP getBetaPiAccept();

		SEXP getThetaZeroAccept();
		SEXP getThetaZeroProp();

		virtual void getPiSamples(int* c, int* b, double* pi);
		double sample_qn(int b, int j);
		double phi (double val, double mu, double sigma2);
		double sn(double x, int b, int j);
		void update_params(double x, int b, int j, int n);
		virtual void getAlphaPiSamples(int* c, double* alpha_pi);
		virtual void getBetaPiSamples(int* c, double* beta_pi);
		virtual void getAlphaPiAccept(int* c, double* acc);
		virtual void getBetaPiAccept(int* c, double* acc);

		virtual void getThetaZeroAccept(int* c, int* b, int* j, double* zero_prop, double* zero_acc);
		
	protected:
		static const char* sColPMweight;

		static const char* sParam_sigma_MH_gamma;
		static const char* sParam_sigma_MH_theta;
		static const char* sParam_sigma_MH_alpha;
		static const char* sParam_sigma_MH_beta;
		static const char* sParam_w_gamma;
		static const char* sParam_w_theta;
		static const char* sParam_w_alpha;
		static const char* sParam_w_beta;

		bool gScreen;
		
		//typedef enum {BB2004 = 1, MH, ADAPT, INDEP} eAlgoType;
		typedef enum {MH = 1, MH_ADAPT, MIS_ADAPT, INDEP} eAlgoType;
		typedef enum {eSim_Type_MH = 1, eSim_Type_SLICE} eSimType;
		
		eAlgoType gAlgo;
		eSimType gSimType;

		// Global Simulation Parameters
		// Give these default values to allow the Gibbs simulator to run even if no or
		// partial values are supplied.
		double gSigma_MH_alpha;
		double gSigma_MH_beta;
		double gDefault_Sigma_MH_gamma;
		double gDefault_Sigma_MH_theta;
		double gW_alpha;
		double gW_alpha_control;
		double gW_beta;
		double gW_beta_control;
		double gDefault_W_gamma;
		double gDefault_W_gamma_control;

		// Weights for MH sampler - these can be adapted
		double **gWp;
		std::map<const char*, gSimParams> mGlobalSimParams;
		double gMH_weight;
		MHAdaptParams gMHAdaptParams;

		double lambda_alpha;
		double lambda_beta;

		int*** gTheta_zero_prop; // Zero proposal and acceptance
		int*** gTheta_zero_acc;

		double* alpha_pi;  // Current value of the sampled distribution
		double* beta_pi;  // Current value of the sampled distribution
		double** gPi;
		
		int* alpha_pi_acc;  // Current value of the sampled distribution
		int* beta_pi_acc;  // Current value of the sampled distribution

		double** alpha_pi_samples;
		double** beta_pi_samples;
		double*** gPi_samples;

		// Adapt MH phase
		int gAdapt_Phase_alpha;
		int gAdapt_Phase_beta;
		int alpha_pi_acc_adapt;
		int beta_pi_acc_adapt;

		int** gAdapt_Phase_theta;
		int** theta_acc_adapt;
		double** theta_mix_p;
		double** theta_max_p;

		// Adaptive MCMC - Algorithm 3 - Ji and Schmidler:
		// Adaptive Markov Chain Monte Carlo for Bayesian Variable Selection
		int gM;
		double** gW0;
		double*** gW;
		double*** gMU;
		double*** gSIGMA2;
		// Leave these out for the time being
		//double** gMU_tilde;
		//double** gSIGMA2_tilde;
		//double** gLAMBDA;

		// Adapt for some of the MH steps
		int in_apapt_phase;
};

#endif
