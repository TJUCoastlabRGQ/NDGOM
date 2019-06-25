# include "mex.h"
# include <math.h>

#include "../../../SWE2d/@SWEAbstract2d/private/mxSWE2d.h"

#define NRHS 9
#define NLHS 1

#ifdef _OPENMP
#include <omp.h>
#endif

/* This function is used to calculate the local or the adjacent surface flux
 * for the vertical momentum equation included in the nonhydrostatic model
 */

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
    signed char* regType = (signed char*)mxGetData(prhs[0]);
    //double *FToE = mxGetPr(prhs[1]);
    double hcrit = mxGetScalar(prhs[2]);
    double *h = mxGetPr(prhs[3]);
    double *hu = mxGetPr(prhs[4]);
    double *hv = mxGetPr(prhs[5]);
    double *hw = mxGetPr(prhs[6]);
    double *nx = mxGetPr(prhs[7]);
    double *ny = mxGetPr(prhs[8]);
    
    mwSize M = mxGetM(prhs[3]);
    mwSize N = mxGetN(prhs[3]);
    /*Program for obtaining the JcsGlobalStiffMatrix and JrsGlobalStiffMatrix*/
    plhs[0] = mxCreateDoubleMatrix( M, N, mxREAL);
    double *surfFlux = mxGetPr(plhs[0]);
    
    mxArray *TempDryFaceFlag = mxCreateDoubleMatrix(1, N, mxREAL);
    double *DryFaceFlag = mxGetPr(TempDryFaceFlag);
    
    evaluateWetDryInterface( regType, prhs[1], DryFaceFlag);
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (mwIndex k = 0; k<N; k++){
        if (DryFaceFlag[k] == 1) {
            // the considered face is one between two adjacent dry cells, or one dry cell and another wet cell, for this situation, do nothing
            continue;
        }
        else{
            for (mwIndex j = 0; j < M; j++){
            	mwIndex i = k*M + j;
                    surfFlux[i] = hu[i] * hw[i]/h[i] * nx[i] + hv[i] * hw[i]/h[i] * ny[i];
            }
        }
    }
}
