# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 6
#define NLHS 3

/* This function is used to set the flux term at the wet-dry interface and the local value of the wet cell adjacent to
 * the dry cell to be zero such that the derivative of the given flow variable at the wet-dry interface is calculated
 * in a local manner.
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
    double *WetDryFaceOrder = mxGetPr(prhs[0]);
    mwSize Num = mxGetM(prhs[0]);
    
    double *NonhydroFmPoint = mxGetPr(prhs[1]);
    mwSize Nfp = mxGetM(prhs[1]);
    double *NonhydroFpPoint = mxGetPr(prhs[2]);
    
    double *Um = mxGetPr(prhs[3]);
    double *Up = mxGetPr(prhs[4]);
    double *fluxS = mxGetPr(prhs[5]);
    mwSize M = mxGetM(prhs[5]);
    mwSize N = mxGetN(prhs[5]);
    
    
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    double *NewUm = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(M, N, mxREAL);
    double *NewUp = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(M, N, mxREAL);
    double *NewFluxS = mxGetPr(plhs[2]);
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (mwIndex i = 0; i<M*N; i++){
        NewUm[i] = Um[i];
        NewUp[i] = Up[i];
        NewFluxS[i] = fluxS[i];
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (mwIndex i = 0; i<Num; i++){
        for ( mwIndex j = 0; j<Nfp; j++ ){
            NewFluxS[ ( (int)WetDryFaceOrder[i] - 1 ) * Nfp + j ] = 0;
        }
        if (NonhydroFmPoint[ ( (int)WetDryFaceOrder[i] - 1 ) * Nfp ] == -1 ){
            for ( mwIndex j = 0; j<Nfp; j++ ){
                NewUm[ ( (int)WetDryFaceOrder[i] - 1 ) * Nfp + j ] = 0;
            }
        }
        else if (NonhydroFpPoint[ ( (int)WetDryFaceOrder[i] - 1 ) * Nfp ] == -1 ) {
            for ( mwIndex j = 0; j<Nfp; j++ ){
                NewUp[ ( (int)WetDryFaceOrder[i] - 1 ) * Nfp + j ] = 0;
            }            
        }
        else{
           mexPrintf("Matlab:%s:InvalidWetDryStatus,\n", __FILE__); 
        }
            
    }
    
}