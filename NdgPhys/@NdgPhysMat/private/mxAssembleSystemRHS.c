#include <mex.h>
#include "../../../NdgMath/NdgMath.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int Np = (int)mxGetScalar(prhs[0]);
	int K = (int)mxGetScalar(prhs[1]);
	int Nvar = (int)mxGetScalar(prhs[2]);
	double *EXa = mxGetPr(prhs[3]);
	int EStage = (int)mxGetNumberOfElements(prhs[3]);
	double *IMa = mxGetPr(prhs[4]);
	int IStage = (int)mxGetNumberOfElements(prhs[4]);
	double *Tempfphys = mxGetPr(prhs[5]);
	double *ExplicitRHS = mxGetPr(prhs[6]);
	double *ImplicitRHS = mxGetPr(prhs[7]);
	double dt = mxGetScalar(prhs[8]);
	const size_t NdimOut = 3;
	const mwSize dimOut[3] = { Np, K, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *SystemRHS = mxGetPr(plhs[0]);
	double *TempRHS = malloc(Np*K*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int i = 0; i < Nvar; i++){
			Add(SystemRHS + i*Np*K + k*Np, SystemRHS + i*Np*K + k*Np, Tempfphys + i*Np*K + k*Np, Np);
			for (int j = 0; j < IStage; j++){
				MultiplyByConstant(TempRHS + k*Np, ExplicitRHS + (i*EStage + j)*Np*K + k*Np, dt*EXa[j], Np);
				Add(SystemRHS + i*Np*K + k*Np, SystemRHS + i*Np*K + k*Np, TempRHS + k*Np, Np);
				MultiplyByConstant(TempRHS + k*Np, ImplicitRHS + (i*IStage + j)*Np*K + k*Np, dt*IMa[j], Np);
				Add(SystemRHS + i*Np*K + k*Np, SystemRHS + i*Np*K + k*Np, TempRHS + k*Np, Np);
			}
			MultiplyByConstant(TempRHS + k*Np, ExplicitRHS + ((Nvar-1)*EStage + EStage-1)*Np*K + k*Np, dt*EXa[EStage-1], Np);
		}
	}
	free(TempRHS);
}