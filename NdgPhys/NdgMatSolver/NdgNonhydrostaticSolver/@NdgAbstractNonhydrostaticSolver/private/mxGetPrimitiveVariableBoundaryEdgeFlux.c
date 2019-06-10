# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 3
#define NLHS 1

/* This function is used to set the flux for the primitive variable at the boundary edge following the numerical flux
 * for the LDG method given by Peraire and Perssony (2008). This flux is set to be zero at the out flow boundary
 * and set to be homogeneous Newmann boundary at other boundaries
 */
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
    double *EidBoundaryType = mxGetPr(prhs[0]);
    double *vector = mxGetPr(prhs[1]);
    double *Um = mxGetPr(prhs[2]);
    
    mwSize M = mxGetM(prhs[2]);
    mwSize N = mxGetN(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    double *FluxS = mxGetPr(plhs[0]);
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (mwIndex i = 0; i < M*N; i++){
        if( EidBoundaryType[i] == -1 ){
            FluxS[i] = 0;
        }
        else{
            FluxS[i] = vector[i] * Um[i];
        }
    }
}