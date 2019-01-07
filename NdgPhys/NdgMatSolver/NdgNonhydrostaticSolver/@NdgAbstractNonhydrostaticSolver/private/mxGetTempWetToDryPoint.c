# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 4
#define NLHS 1


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  /* check input & output */
  if (nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NLHS);
  }
  
  int Np = (int)mxGetScalar(prhs[0]);
  int Nfp = (int)mxGetScalar(prhs[1]);
  double *NewWetDyrFace = mxGetPr(prhs[2]);
  mwSize M = mxGetM(prhs[2]);
  double *Fmask = mxGetPr(prhs[3]);

  plhs[0] = mxCreateDoubleMatrix(M*Nfp, 1, mxREAL);
  double *TempWetToDryPoint = mxGetPr(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  // get the global point index of the the point that located on the face of wet cell and adjacent to a dry cell 
for(mwIndex i = 0; i<M; i++){
	int cell = (int)NewWetDyrFace[i];
	int face = (int)NewWetDyrFace[M+i];
	for(mwIndex j= 0; j<Nfp; j++){
		TempWetToDryPoint[i*Nfp+j] = (cell-1)*Np + Fmask[(face-1)*Nfp+j];
	}
  }
}