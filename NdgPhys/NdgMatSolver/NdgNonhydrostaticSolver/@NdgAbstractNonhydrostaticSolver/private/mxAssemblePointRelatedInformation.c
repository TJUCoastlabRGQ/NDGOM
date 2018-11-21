# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 4
#define NLHS 2

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
	double *ZeroFluxBoundary = mxGetPr(prhs[0]);
	double *AdjacentDryCellAndFace = mxGetPr(prhs[1]);
	double *FToE = mxGetPr(prhs[2]);
	double *FToN1 = mxGetPr(prhs[3]);
	mwSize M = mxGetM(prhs[3]);
	mwSize N = mxGetN(prhs[3]);
	mwSize NumCell = mxGetM(prhs[1]);
	plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(M, N, mxREAL);
    double *NonhydroFmPoint = mxGetPr(plhs[0]);
	double *NonhydroFpPoint = mxGetPr(plhs[1]);
	mwSize Ne = mxGetN(prhs[2]);
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif	
  for(mwIndex i = 0; i<M*N; i++){
     NonhydroFmPoint[i] = 1;
	 NonhydroFpPoint[i] = 1;
  }
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif		
	for(mwIndex i = 0; i<NumCell; i++){
		int ele1 = (int)ZeroFluxBoundary[i];
		int ele2 = (int)AdjacentDryCellAndFace[i];
		for(mwIndex j = 0; j<Ne; j++){
			/*For this situation, the local element is wet, while the adjacent cell is dry*/
			if((ele1 == FToE[j*2])&&(ele2 == FToE[j*2+1])){
				for(mwIndex k=0; k < M; k++)
					NonhydroFpPoint[j*M+k] = -1;
			}
			/*For this situation, the adjacent element is wet, while the local cell is dry*/
			else if((ele2 == FToE[j*2])&&(ele1 == FToE[j*2+1])){
				for(mwIndex k=0; k < M; k++)
					NonhydroFmPoint[j*M+k] = -1;
			}
		}
	}
}