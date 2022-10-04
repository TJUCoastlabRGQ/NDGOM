#include "SWENonhydrostatic3d.h"
//#include "../../../../../NdgMath/NdgMemory.h"

//extern double *InvSHeight;

void GetInverseHeight(double *, double *, double, int);

/*The input variables are organized as follows:
PUPX: $\frac{\partial u}{\partial x}$, indexed as 0.
PUPS: $\frac{\partial u}{\partial \sigma}$, indexed as 1.
PVPY: $\frac{\partial v}{\partial y}$, indexed as 2.
PVPS: $\frac{\partial v}{\partial \sigma}$, indexed as 3.
PWPS: $\frac{\partial w}{\partial \sigma}$, indexed as 4.
PSPX: $\frac{\partial \sigma}{\partial x}$, indexed as 5.
PSPY: $\frac{\partial \sigma}{\partial y}$, indexed as 6.
   h: The three dimensional water depth, indexed as 7.
  dt: The time step, indexed as 8.
 rho: The water density, indexed as 9.
Hcrit: The critical water depth, indexed as 10.
*/

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
	double *J = mxGetPr(prhs[11]);
	double *M3d = mxGetPr(prhs[12]);

	plhs[0] = mxCreateDoubleMatrix(Np*K, 1, mxREAL);
	double *RHS = mxGetPr(plhs[0]);

	double *InvSHeight = malloc(Np*K*sizeof(double));
	double *TempRHS = malloc(Np*K*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		GetInverseHeight(InvSHeight + k*Np, Height+k*Np, Hcrit, Np);
	}

	ptrdiff_t One = 1;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		double *EleMass3d = malloc(Np*Np*sizeof(double));
		DiagMultiply(EleMass3d, M3d, J+k*Np, Np);
		for (int p = 0; p < Np; p++){
			TempRHS[k*Np + p] = rho / dt*(PUPX[k*Np + p] + PUPS[k*Np + p] * PSPX[k*Np + p] + \
				PVPY[k*Np + p] + PVPS[k*Np + p] * PSPY[k*Np + p] + \
				InvSHeight[k*Np + p] * PWPS[k*Np + p]);
		}
		MatrixMultiply("N", "N", (ptrdiff_t)Np, One, (ptrdiff_t)Np, 1.0, EleMass3d,
			(ptrdiff_t)Np, TempRHS + k*Np, (ptrdiff_t)Np, 0.0, RHS + k*Np, (ptrdiff_t)Np);
		free(EleMass3d);
	}
	free(InvSHeight);
	free(TempRHS);
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