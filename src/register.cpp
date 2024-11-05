#include<cstdlib>
#include<cstdio>

#include "c212_Rdefines.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
#include <R_ext/Visibility.h>

#include "register.h"
#include "c212_exec.h"

//
// Register native methods with R (entry-points into compiled code)
//

//static const char *rcsId = "$Id: register.cpp,v 1.2 2017/03/28 15:38:23 clb13102 Exp clb13102 $";

// .C methods
static R_CMethodDef cMethods[] = {
	{"Release", (DL_FUNC)&Release, 0, NULL},
	{"Release_Interim", (DL_FUNC)&Release_Interim, 0, NULL},
	{"getAlphaPiAcceptInterim", (DL_FUNC)&getAlphaPiAcceptInterim, 3, NULL},
	{"getAlphaPiSamplesInterim", (DL_FUNC)&getAlphaPiSamplesInterim, 3, NULL},
	{"getBetaPiAcceptInterim", (DL_FUNC)&getBetaPiAcceptInterim, 3, NULL},
	{"getBetaPiSamplesInterim", (DL_FUNC)&getBetaPiSamplesInterim, 3, NULL},
	{"getGammaAcceptInterim", (DL_FUNC)&getGammaAcceptInterim, 5, NULL},
	{"getGammaSamplesInterim", (DL_FUNC)&getGammaSamplesInterim, 5, NULL},
	{"getMuGamma0SamplesInterim", (DL_FUNC)&getMuGamma0SamplesInterim, 3, NULL},
	{"getMuGammaSamplesInterim", (DL_FUNC)&getMuGammaSamplesInterim, 4, NULL},
	{"getMuTheta0SamplesInterim", (DL_FUNC)&getMuTheta0SamplesInterim, 3, NULL},
	{"getMuThetaSamplesInterim", (DL_FUNC)&getMuThetaSamplesInterim, 4, NULL},
	{"getPiSamplesInterim", (DL_FUNC)&getPiSamplesInterim, 4, NULL},
	{"getSigma2GammaSamplesInterim", (DL_FUNC)&getSigma2GammaSamplesInterim, 4, NULL},
	{"getSigma2ThetaSamplesInterim", (DL_FUNC)&getSigma2ThetaSamplesInterim, 4, NULL},
	{"getTau2Gamma0SamplesInterim", (DL_FUNC)&getTau2Gamma0SamplesInterim, 3, NULL},
	{"getTau2Theta0SamplesInterim", (DL_FUNC)&getTau2Theta0SamplesInterim, 3, NULL},
	{"getThetaAcceptInterim", (DL_FUNC)&getThetaAcceptInterim, 5, NULL},
	{"getThetaSamplesInterim", (DL_FUNC)&getThetaSamplesInterim, 5, NULL},
	{NULL, NULL, 0, NULL},
};

// .Call methods
static R_CallMethodDef callMethods[] = {
	{"c2121a_exec", (DL_FUNC)&c2121a_exec, 36},											// 36
	{"c2121a_interim_hier2_exec", (DL_FUNC)&c2121a_interim_hier2_exec, 33},				// 33
	{"c2121a_poisson_mc_exec", (DL_FUNC)&c2121a_poisson_mc_exec, 41},					// 41
	{"c212BB_exec", (DL_FUNC)&c212BB_exec, -1},											// 48
	{"c212BB_interim_hier2_exec", (DL_FUNC)&c212BB_interim_hier2_exec, 39},				// 39
	{"c212BB_poisson_mc_exec", (DL_FUNC)&c212BB_poisson_mc_exec, 49},					// 49
	{"getAlphaPiAcceptAll", (DL_FUNC)&getAlphaPiAcceptAll, 0},							// 0
	{"getAlphaPiAcceptInterimAll", (DL_FUNC)&getAlphaPiAcceptInterimAll, 0},			// 0
	{"getAlphaPiSamplesAll", (DL_FUNC)&getAlphaPiSamplesAll, 0},						// 0
	{"getAlphaPiSamplesInterimAll", (DL_FUNC)&getAlphaPiSamplesInterimAll, 0},			// 0
	{"getBetaPiAcceptAll", (DL_FUNC)&getBetaPiAcceptAll, 0},							// 0
	{"getBetaPiAcceptInterimAll", (DL_FUNC)&getBetaPiAcceptInterimAll, 0},				// 0
	{"getBetaPiSamplesAll", (DL_FUNC)&getBetaPiSamplesAll, 0},							// 0
	{"getBetaPiSamplesInterimAll", (DL_FUNC)&getBetaPiSamplesInterimAll, 0},			// 0
	{"getGammaAcceptAll", (DL_FUNC)&getGammaAcceptAll, 0},								// 0
	{"getGammaAcceptInterimAll", (DL_FUNC)&getGammaAcceptInterimAll, 0},				// 0
	{"getGammaSamplesAll", (DL_FUNC)&getGammaSamplesAll, 0},							// 0
	{"getGammaSamplesInterimAll", (DL_FUNC)&getGammaSamplesInterimAll, 0},				// 0
	{"getMuGamma0SamplesAll", (DL_FUNC)&getMuGamma0SamplesAll, 0},						// 0
	{"getMuGamma0SamplesInterimAll", (DL_FUNC)&getMuGamma0SamplesInterimAll, 0},		// 0
	{"getMuGammaSamplesAll", (DL_FUNC)&getMuGammaSamplesAll, 0},						// 0
	{"getMuGammaSamplesInterimAll", (DL_FUNC)&getMuGammaSamplesInterimAll, 0},			// 0
	{"getMuTheta0SamplesAll", (DL_FUNC)&getMuTheta0SamplesAll, 0},						// 0
	{"getMuTheta0SamplesInterimAll", (DL_FUNC)&getMuTheta0SamplesInterimAll, 0},		// 0
	{"getMuThetaSamplesAll", (DL_FUNC)&getMuThetaSamplesAll, 0},						// 0
	{"getMuThetaSamplesInterimAll", (DL_FUNC)&getMuThetaSamplesInterimAll, 0},			// 0
	{"getPiSamplesAll", (DL_FUNC)&getPiSamplesAll, 0},									// 0
	{"getPiSamplesInterimAll", (DL_FUNC)&getPiSamplesInterimAll, 0},					// 0
	{"getSigma2GammaSamplesAll", (DL_FUNC)&getSigma2GammaSamplesAll, 0},				// 0
	{"getSigma2GammaSamplesInterimAll", (DL_FUNC)&getSigma2GammaSamplesInterimAll, 0},	// 0
	{"getSigma2ThetaSamplesAll", (DL_FUNC)&getSigma2ThetaSamplesAll, 0},				// 0
	{"getSigma2ThetaSamplesInterimAll", (DL_FUNC)&getSigma2ThetaSamplesInterimAll, 0},	// 0
	{"getTau2Gamma0SamplesAll", (DL_FUNC)&getTau2Gamma0SamplesAll, 0},					// 0
	{"getTau2Gamma0SamplesInterimAll", (DL_FUNC)&getTau2Gamma0SamplesInterimAll, 0},	// 0
	{"getTau2Theta0SamplesAll", (DL_FUNC)&getTau2Theta0SamplesAll, 0},					// 0
	{"getTau2Theta0SamplesInterimAll", (DL_FUNC)&getTau2Theta0SamplesInterimAll, 0},	// 0
	{"getThetaAcceptAll", (DL_FUNC)&getThetaAcceptAll, 0},								// 0
	{"getThetaAcceptInterimAll", (DL_FUNC)&getThetaAcceptInterimAll, 0},				// 0
	{"getThetaSamplesAll", (DL_FUNC)&getThetaSamplesAll, 0},							// 0
	{"getThetaSamplesInterimAll", (DL_FUNC)&getThetaSamplesInterimAll, 0},				// 0
	{"getThetaZeroAcceptAll", (DL_FUNC)&getThetaZeroAcceptAll, 0},						// 0
	{"getThetaZeroPropAll", (DL_FUNC)&getThetaZeroPropAll, 0},							// 0
    {NULL, NULL, 0}
};

// Register the methods with R
void attribute_visible R_init_c212(DllInfo *info)
{
	// Register the entry-points with R
	R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
	// Ensure only registered entry-points can be called by name e.g. .C("Release")
	R_useDynamicSymbols(info, FALSE);
}
