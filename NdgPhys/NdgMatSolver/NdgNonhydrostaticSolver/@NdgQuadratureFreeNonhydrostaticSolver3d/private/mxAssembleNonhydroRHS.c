#include "SWENonhydrostatic3d.h"
//#include "../../../../../NdgMath/NdgMemory.h"

//extern double *InvSHeight;

void GetInverseHeight(double *, double *, double, int);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *PUPX = mxGetPr(prhs[0]);
	int Np = (int)mxGetM(prhs[0]);
	int  K = (int)mxGetN(prhs[0]);
	double *PUPS = mxGetPr(prhs[1]);
	double *PVPY = mxGetPr(prhs[2]);
	double *PVPS = mxGetPr(prhs[3]);
	double *PWPS = mxGetPr(prhs[4]);
	double *PSPX = mxGetPr(prhs[5]);
	double *PSPY = mxGetPr(prhs[6]);
	double *Height = mxGetPr(prhs[7]);
	double dt = mxGetScalar(prhs[8]);
	double rho = mxGetScalar(prhs[9]);
	double Hcrit = mxGetScalar(prhs[10]);

	plhs[0] = mxCreateDoubleMatrix(Np*K, 1, mxREAL);
	double *RHS = mxGetPr(plhs[0]);

	double *InvSHeight = malloc(Np*K*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		GetInverseHeight(InvSHeight + k*Np, Height+k*Np, Hcrit, Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int p = 0; p < Np; p++){
			RHS[k*Np + p] = rho / dt*(PUPX[k*Np + p] + PUPS[k*Np + p] * PSPX[k*Np + p] + \
				PVPY[k*Np + p] + PVPS[k*Np + p] * PSPY[k*Np + p] + \
				InvSHeight[k*Np + p] * PWPS[k*Np + p]);
		}
	}
	free(InvSHeight);

}

void GetInverseHeight(double *dest, double *source, double hcrit, int Np){
	for (int i = 0; i < Np; i++){
		if (source[i] >= hcrit) {
			dest[i] = 1.0 / source[i];
		}
		else{
			dest[i] = 0;
		}	
	}
}