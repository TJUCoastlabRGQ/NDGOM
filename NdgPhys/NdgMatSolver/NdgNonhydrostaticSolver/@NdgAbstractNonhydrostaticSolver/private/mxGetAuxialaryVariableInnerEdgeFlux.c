# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 10
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
    double *WetDryFaceOrder = mxGetPr(prhs[0]);
    mwSize NumFace = mxGetM(prhs[0]);
    double *fluxS = mxGetPr(prhs[1]);
    double Tau = mxGetScalar(prhs[2]);
    double *NonhydroFmPoint = mxGetPr(prhs[3]);
    double *NonhydroFpPoint = mxGetPr(prhs[4]);
    mwSize Nfp = mxGetM(prhs[2]);
    mwSize Ne = mxGetN(prhs[2]);
    double *Um = mxGetPr(prhs[5]);
    double *Up = mxGetPr(prhs[6]);
    double *Sigmam = mxGetPr(prhs[7]);
    double *Sigmap = mxGetPr(prhs[8]);  
    double *vector = mxGetPr(prhs[9]);
    
    plhs[0] = mxCreateDoubleMatrix(Nfp, Ne, mxREAL);
    double *NewfluxS = mxGetPr(plhs[0]);
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif    
    for(mwIndex i = 0; i<Nfp*Ne; i++){
        NewfluxS[i] = fluxS[i];
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for(mwIndex i = 0; i < NumFace; i++){
        for( mwIndex j = (int)( WetDryFaceOrder[i] - 1 )*Nfp; j<(int)WetDryFaceOrder[i]*Nfp; j++ ){
            if( NonhydroFmPoint[j] == NonhydroFpPoint[j] )
                mexPrintf("Matlab:%s:InvalidWetDryStatus,\n", __FILE__);
            if(NonhydroFmPoint[j] == -1)
                // Homogeneous Dirichlet boundary condition for pressure 
                NewfluxS[j] = -1 * ( Sigmap[j] + Tau * Up[j] * vector[j] ) * vector[j];
            if(NonhydroFpPoint[j] == -1)
                NewfluxS[j] = ( Sigmam[j] - Tau * Um[j] * vector[j] ) * vector[j];
            
        }
    }
    
}