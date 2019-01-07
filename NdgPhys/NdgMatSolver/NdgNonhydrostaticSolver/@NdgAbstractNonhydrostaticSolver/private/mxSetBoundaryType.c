# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 2
#define NLHS 1

typedef enum {
  NdgEdgeInner = 0,
  NdgEdgeGaussEdge = 1,
  NdgEdgeSlipWall = 2,
  NdgEdgeNonSlipWall = 3,
  NdgEdgeZeroGrad = 4,
  NdgEdgeClamped = 5,
  NdgEdgeClampedDepth = 6,
  NdgEdgeClampedVel = 7,
  NdgEdgeFlather = 8,
  NdgEdgeNonLinearFlather = 9,
  NdgEdgeNonLinearFlatherFlow = 10,
  NdgEdgeNonReflectingFlux = 11
} NdgEdgeType;

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
  
  signed char *ftype = (signed char *)mxGetData(prhs[0]);
  double *FToN1 = mxGetPr(prhs[1]);
  mwSize M = mxGetM(prhs[1]);
  mwSize N = mxGetN(prhs[1]);
  plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
  double *EidBoundaryType = mxGetPr(plhs[0]);
  
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif	  
 
 for(mwIndex e = 0; e< N; e++){
	 NdgEdgeType type = (NdgEdgeType)ftype[e];
	 if (type == NdgEdgeSlipWall || type == NdgEdgeNonSlipWall){
		 for (mwIndex n = 0; n < M; n++){
			 EidBoundaryType[e*M + n] = 1;
		 }
	 }
	 else if (type == NdgEdgeZeroGrad || type == NdgEdgeClamped || type == NdgEdgeClampedDepth || type == NdgEdgeClampedVel){
		 for (mwIndex n = 0; n < M; n++){
			 EidBoundaryType[e*M + n] = -1;
		 }
	 }
		 else{
             mexPrintf("Unknown or Unrealized boundary type\n");			 
		 }
	 }
 }