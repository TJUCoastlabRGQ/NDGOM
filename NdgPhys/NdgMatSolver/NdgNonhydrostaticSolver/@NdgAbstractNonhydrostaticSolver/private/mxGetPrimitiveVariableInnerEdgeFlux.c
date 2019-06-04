# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 3
#define NLHS 1

/* This function is used to set the flux for the primitive at the wet-dry interface to be zero following the numerical flux 
 for the LDG method given by Peraire and Perssony (2008). This flux is set to be zero because the homogeneous Dirichlet boundary 
 conditon for non-hydrostatic pressure is imposed*/

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
    double *WetDryFaceOrder = mxGetPr(prhs[0]);
    mwSize Num = mxGetM(prhs[0]);
    
    double *fluxS = mxGetPr(prhs[1]);
    int Nfp = (int)mxGetScalar(prhs[2]);
    
    mwSize M = mxGetM(prhs[1]);
    mwSize N = mxGetN(prhs[1]);
    
    
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    double *NewFluxS = mxGetPr(plhs[0]);
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (mwIndex i = 0; i<M*N; i++){
        NewFluxS[i] = fluxS[i];
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif    
    for (mwIndex i = 0; i<Num; i++){
        for ( mwIndex j = 0; j<Nfp; j++ ){
            NewFluxS[ ( (int)WetDryFaceOrder[i] - 1 ) * Nfp + j ] = 0;
        }
    }
    
}