# include "mex.h"
# include <math.h>
#include <malloc.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 6
#define NLHS 1

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	  /* check input & output */
  if (nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NLHS);
  }
  double *ZeroFluxBoundary = mxGetPr(prhs[0]); 
  int ZeroFluxBoundaryIndex = (int)mxGetScalar(prhs[1]);
  double *AdjacentDryCellAndFace = mxGetPr(prhs[2]);
  int Nfp = (int)mxGetScalar(prhs[3]);
  double *FluxTerm = mxGetPr(prhs[4]);
  mwSize row = mxGetM(prhs[4]); 
  mwSize col = mxGetN(prhs[4]);
  int TNfp = (int)mxGetScalar(prhs[5]);
  plhs[0] = mxCreateDoubleMatrix(row,col,mxREAL);
  double *NewFluxTerm = mxGetPr(plhs[0]);
  
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif
 
 for(mwIndex i =0; i < ZeroFluxBoundaryIndex; i++){
	 int cell = (int)ZeroFluxBoundary[i];
	 int DryCell = (int)AdjacentDryCellAndFace[i];
	 int edge = (int)ZeroFluxBoundary[i + ZeroFluxBoundaryIndex];
	 int DryEdge = (int)AdjacentDryCellAndFace[i + ZeroFluxBoundaryIndex];
	 for (mwIndex j = 0; j < Nfp; j++){
		 FluxTerm[(cell - 1)*TNfp + (edge - 1)*Nfp + j] = 0;
		 FluxTerm[(DryCell - 1)*TNfp + (DryEdge - 1)*Nfp + j] = 0;
	 }
		 
 }
 
 /*data copy*/
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif
  for(mwIndex i =0; i < row*col; i++)
    NewFluxTerm[i] = FluxTerm[i];
}