# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 4
#define NLHS 1

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
  signed char *ftype = (signed char *)mxGetData(prhs[2]);
  NdgWetDryFaceBoundaryType type = (NdgWetDryFaceBoundaryType)ftype[0];
  double *EidBoundaryType = mxGetPr(prhs[3]);
  mwSize M = mxGetM(prhs[0]);
  mwSize N = mxGetN(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
  double *AdjacentValue = mxGetPr(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  if (type == NdgWetDryFaceZero){	
	  for (mwIndex i = 0; i<N; i++){
		  if ( EidBoundaryType[i*M] == -1){
			  for (mwIndex j = 0; j < M; j++)
				  /*Program for obtaining the local value with zero and zero grad condition considered*/
				  /*if the first point of a given edge is the given face type, then the whole face is of that type*/
				  AdjacentValue[i*M + j] = -1 * fm[i*M + j];
		  }
		  else if ( EidBoundaryType[i*M] != -1){
			  for (mwIndex j = 0; j<M; j++)
				  AdjacentValue[i*M + j] = fm[i*M + j];
		  }

	  }
  }
  else if (type == NdgWetDryFaceZeroGrad){
	  for (mwIndex i = 0; i < N; i++){
		  if ( EidBoundaryType[i*M] == -1){
			  for (mwIndex j = 0; j < M; j++)
				  /*Program for obtaining the local value with zero and zero grad condition considered*/
				  /*if the first point of a given edge is the given face type, then the whole face is of that type*/
				  AdjacentValue[i*M + j] = fm[i*M + j];
		  }
		  else if (EidBoundaryType[i*M] != -1){
			  for (mwIndex j = 0; j<M; j++)
				  AdjacentValue[i*M + j] = -1 * fm[i*M + j];
		  }
	  }
  }
 }