#include "mex.h"
#include <math.h>
#include "../../../../../NdgMath/NdgMemory.h"
#include "../../../../../NdgMath/NdgMath.h"
#include "eqstate.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double *BaroclinicT, *BaroclinicS, *BaroclinicDTS;

extern char *BaroDensityInitialized;

void MyExit()
{
    if (!strcmp("True", BaroDensityInitialized)) {
        BaroDensityMemoryDeAllocation();
        BaroDensityInitialized = "False";
    }
    return;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    mexAtExit(&MyExit);
    double *height = mxGetPr(prhs[0]);
    double *hT = mxGetPr(prhs[1]);
    double *hS = mxGetPr(prhs[2]);
    double *z = mxGetPr(prhs[3]);
    double hcrit = mxGetScalar(prhs[4]);
    int Np = (int)mxGetScalar(prhs[5]);
    int K = (int)mxGetScalar(prhs[6]);
	//This is for linear EOS
	double rho0 = mxGetScalar(prhs[7]);
	//This is for the thermal expansion coefficient
	double alphaT = mxGetScalar(prhs[8]);
	//This is for the salinity expansion coefficient
	double betaS = mxGetScalar(prhs[9]);
	//This is for linear EOS
	double T0 = mxGetScalar(prhs[10]);
	//This is for linear EOS
	double S0 = mxGetScalar(prhs[11]);

	char *EosType;
	EosType = mxArrayToString(prhs[12]);
    
    plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
    double *rho = mxGetPr(plhs[0]);
    
    if (!strcmp("False", BaroDensityInitialized)) {
        BaroDensityMemoryAllocation(Np, K);
    }
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int i = 0; i < K; i++) {

		DotCriticalDivide(BaroclinicT + i*Np, hT + i*Np, &hcrit, height + i*Np, Np);

		DotCriticalDivide(BaroclinicS + i*Np, hS + i*Np, &hcrit, height + i*Np, Np);

		if (!strcmp(EosType, "Jackett05")) {
			for (int p = 0; p < Np; p++) {
				EosByFeistel(rho + i*Np + p, max(*(BaroclinicT + i*Np + p), 0.0), max(*(BaroclinicS + i*Np + p), 0.0) );
			}
		}
		else if (!strcmp(EosType, "UNESCO83")) {
			for (int p = 0; p < Np; p++) {
				EosByUNESCO(rho + i*Np + p, max(*(BaroclinicT + i*Np + p), 0.0), max(*(BaroclinicS + i*Np + p), 0.0));
			}
		}
		else if (!strcmp(EosType, "Linear")) {
			for (int p = 0; p < Np; p++) {
				EosByLinear(rho + i*Np + p, max(*(BaroclinicT + i*Np + p), 0.0), max(*(BaroclinicS + i*Np + p), 0.0), rho0, T0, S0, alphaT, betaS);
			}
		}
		else {
			printf("Equation of state(EOS) needs to be pointed for this part!\n");
			break;
		}
    }
}