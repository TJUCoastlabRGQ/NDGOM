#include "mex.h"
                
#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 5
#define NLHS 2

/*Function used to get the volume flux for the vertical momentum included in 
 * the non-hydrostatic model. It's noted that only the wet cell is conidered, 
 * and the volume flux for other cells are set to be zero*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    /* check input & output */
    if (nrhs != NRHS) {
        mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NRHS);
    }
    
    if (nlhs != NLHS) {
        mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NLHS);
    }
    
    double *status = mxGetPr(prhs[0]);
    double *hu = mxGetPr(prhs[1]);
    double *hv = mxGetPr(prhs[2]);
    double *hw = mxGetPr(prhs[3]);
    double *h = mxGetPr(prhs[4]);
    
    mwSize Np = mxGetM(prhs[4]);
    mwSize Ne = mxGetN(prhs[4]);
    
    plhs[0] = mxCreateDoubleMatrix(Np, Ne, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Np, Ne, mxREAL);
    double *huw = mxGetPr(plhs[0]);
    double *hvw = mxGetPr(plhs[1]);
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for(mwIndex i = 0; i< Ne; i++){
        if (status[i] == 4){ // The considered cell is a wet cell
            for(mwIndex j = i*Np; j< (i+1)*Np; j++ ){
                huw[j] = hu[j] * hw[j] / h[j];
                hvw[j] = hv[j] * hw[j] / h[j];
            }
        }
    }
    
    return;
}