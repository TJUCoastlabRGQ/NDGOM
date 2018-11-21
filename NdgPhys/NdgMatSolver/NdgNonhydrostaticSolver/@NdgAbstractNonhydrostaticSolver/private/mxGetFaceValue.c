# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 5
#define NLHS 2

typedef enum {
  NdgWetDryFaceZero = 0,
  NdgWetDryFaceZeroGrad = 1,
} NdgWetDryFaceBoundaryType;

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
  
  double *fm = mxGetPr(prhs[0]);
  double *fp = mxGetPr(prhs[1]);
  double *NonhydroFmPoint = mxGetPr(prhs[2]);
  double *NonhydroFpPoint = mxGetPr(prhs[3]);
  signed char *ftype = (signed char *)mxGetData(prhs[4]);
  NdgWetDryFaceBoundaryType type = (NdgWetDryFaceBoundaryType)ftype[0];
  mwSize M = mxGetM(prhs[1]);
  mwSize N = mxGetN(prhs[1]);
  plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
  double *LocalValue = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(M, N, mxREAL);
  double *AdjacentValue = mxGetPr(plhs[1]);
  
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif	
for(mwIndex i=0; i<N; i++){
	for(mwIndex j=0; j<M; j++){
		/*Program for obtaining the local value with zero and zero grad condition considered*/
		/*if the first point of a given edge is the given face type, then the whole face is of that type*/
	if(type == NdgWetDryFaceZero && NonhydroFmPoint[i*M] == -1)
		LocalValue[i*M + j] = -1* fp[i*M + j];	
    else if(type == NdgWetDryFaceZeroGrad && NonhydroFmPoint[i*M] == -1)
		LocalValue[i*M + j] = 1* fp[i*M + j];
	else
		LocalValue[i*M + j] = fm[i*M + j];
       /*Program for obtaining the adjacent value with zero and zero grad condition considered*/
	if(type == NdgWetDryFaceZero && NonhydroFpPoint[i*M] == -1)
		AdjacentValue[i*M + j] = -1* fm[i*M + j];	
    else if(type == NdgWetDryFaceZeroGrad && NonhydroFpPoint[i*M] == -1)
		AdjacentValue[i*M + j] = 1* fm[i*M + j];
	else
		AdjacentValue[i*M + j] = fp[i*M + j];	
	}
}
}