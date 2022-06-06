# include "mex.h"
# include <math.h>

#define NRHS 2
#define NLHS 2

#ifdef _OPENMP
#include <omp.h>
#endif


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
	double Num = mxGetScalar(prhs[0]);
	mwIndex *jc = mxGetJc(prhs[1]);
	mwIndex *jr = mxGetIr(prhs[1]);

	/*Program for obtaining the JcsGlobalStiffMatrix and JrsGlobalStiffMatrix*/
	plhs[0] = mxCreateDoubleMatrix((int)Num + 1, 1, mxREAL);
	double *JcStiffMatrix = mxGetPr(plhs[0]);

  int Uplimit = (int)Num+1;	
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif
	for (mwIndex i = 0; i < Uplimit; i++){
		JcStiffMatrix[i] = (double)jc[i];
	}	
	plhs[1] = mxCreateDoubleMatrix(jc[(int)Num], 1, mxREAL);
	double *JrStiffMatrix = mxGetPr(plhs[1]);
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif
	for (mwIndex i = 0; i < jc[(int)Num]; i++){
		JrStiffMatrix[i] = (double)jr[i];
	}
}