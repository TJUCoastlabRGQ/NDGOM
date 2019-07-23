# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 5
#define NLHS 1

/* This function is used to set the flux for the auxialary variable at the wet-dry interface following the numerical flux
 * for the LDG method given by Peraire and Perssony (2008). For this part only the homogeneous dirichlet boundary condition
 * considered
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
    mwSize Nfp = mxGetM(prhs[0]);
    mwSize Ne = mxGetN(prhs[0]);
    double Tau = mxGetScalar(prhs[1]);
    double *Um = mxGetPr(prhs[2]);
    double *Sigmam = mxGetPr(prhs[3]);
    double *vector = mxGetPr(prhs[4]);
    
    plhs[0] = mxCreateDoubleMatrix(Nfp, Ne, mxREAL);
    double *fluxS = mxGetPr(plhs[0]);
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    
    for(mwIndex i = 0; i < Ne; i++){
        // The first point of the studied boundary edge is flaged as a Dirichlet boundary
        if ( (int)EidBoundaryType[i*Nfp] == -1 ){
            for( mwIndex j = i*Nfp; j< (i+1)*Nfp; j++ ){
                // Homogeneous Dirichlet boundary condition for pressure
                fluxS[j] = ( Sigmam[j] - Tau * Um[j] * vector[j] ) * vector[j];
            }
        }
    }
    
}