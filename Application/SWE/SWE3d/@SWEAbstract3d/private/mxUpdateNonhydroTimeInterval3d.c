#include "mex.h"
#include "mxSWE3d.h"
#include <math.h>
#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgMemory.h"

extern double *TimeIntervalu, *TimeIntervalv, *TimeIntervalw, *TimeIntervalt;

extern char *TimeIntervalInitialized;

void MyExit()
{
	if (!strcmp("True", TimeIntervalInitialized)){
		SWENH3dTimeIntervalMemoryDeAllocation();
		TimeIntervalInitialized = "False";
	}
	return;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  mexAtExit(&MyExit);
  double gra = mxGetScalar(prhs[0]);
  double N = mxGetScalar(prhs[1]);
  double Nz = mxGetScalar(prhs[2]);
  int Np = (int)mxGetScalar(prhs[3]);
  double* dx = mxGetPr(prhs[4]);
  signed char* regionType = (signed char*)mxGetPr(prhs[5]);
  // double* fphys = mxGetPr(prhs[5]);
  double *Hu = mxGetPr(prhs[6]);
  double *Hv = mxGetPr(prhs[7]);
  double *H = mxGetPr(prhs[8]);
  double *Hw = mxGetPr(prhs[9]);
  int NLayer = (int)mxGetScalar(prhs[10]);
  int K2d = (int)mxGetScalar(prhs[11]);

  int K3d = NLayer*K2d;

  if (!strcmp("False", TimeIntervalInitialized)){
	  SWENH3dTimeIntervalMemoryAllocation(Np, K3d, K2d);
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < K2d; k++){
	  double spe;
	  TimeIntervalt[k] = pow(10.0, 6.0);
	  NdgRegionType type = (NdgRegionType)regionType[k];
	  if (type == NdgRegionDry) {
		  continue;
	  }
	  for (int L = 0; L < NLayer; L++){
		  /*$u$*/
		  DotDivide(TimeIntervalu + k*NLayer*Np + L*Np, Hu + k*NLayer*Np + L*Np, H + k*NLayer*Np + L*Np, Np);
		  /*$v$*/
		  DotDivide(TimeIntervalv + k*NLayer*Np + L*Np, Hv + k*NLayer*Np + L*Np, H + k*NLayer*Np + L*Np, Np);
		  /*$w$*/
		  DotDivide(TimeIntervalw + k*NLayer*Np + L*Np, Hw + k*NLayer*Np + L*Np, H + k*NLayer*Np + L*Np, Np);
		  for (int p = 0; p < Np; p++){
			  /*$\sqrt{(u^2+v^2)}$*/
			  spe = sqrt(pow(TimeIntervalu[k*NLayer*Np + L*Np + p], 2.0) + pow(TimeIntervalv[k*NLayer*Np + L*Np + p], 2.0));
			  /*$\frac{dx}{(\sqrt{(u^2+v^2)} + \sqrt{gh})(2N+1)}$*/
			  TimeIntervalt[k] = min(TimeIntervalt[k], dx[k] / (spe + sqrt(gra * H[k*NLayer*Np + L*Np + p])) / (2 * N + 1)* pow(2.0, -1.0 / N));
			  /*$\frac{\Delta \sigma h}{w(2N_z+1)}$*/
			  TimeIntervalt[k] = min(TimeIntervalt[k], 1.0 / NLayer *H[k*NLayer*Np + L*Np + p] / fabs(TimeIntervalw[k*NLayer*Np + L*Np + p])/(2.0*Nz+1) * pow(2.0,-1.0/Nz));
		  }
	  }
  }

  plhs[0] = mxCreateDoubleScalar(1e6);
  for (int k = 0; k < K2d; k++){
	  *mxGetPr(plhs[0]) = min(*mxGetPr(plhs[0]), TimeIntervalt[k]);
  }

  return;
}