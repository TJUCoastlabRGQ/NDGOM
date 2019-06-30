#include "mex.h"
#include "../../../SWE1d/@SWEAbstract1d/private/mxSWE1d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 9
#define NLHS 1

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
    
    signed char* regType = (signed char*)mxGetData(prhs[0]);
    //double *FToE = mxGetPr(prhs[1]);
    double *hm = mxGetPr(prhs[2]);
    double *hum = mxGetPr(prhs[3]);
    //double *hvm = mxGetPr(prhs[4]);
    double *hwm = mxGetPr(prhs[4]);
    double *hp = mxGetPr(prhs[5]);
    double *hup = mxGetPr(prhs[6]);
    //double *hvp = mxGetPr(prhs[8]);
    double *hwp = mxGetPr(prhs[7]);
    double *nx = mxGetPr(prhs[8]);
    //double *ny = mxGetPr(prhs[11]);
    
    mwSize Nfp = mxGetM(prhs[2]);
    mwSize Ne = mxGetN(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(Nfp, Ne, mxREAL);
    double *fluxS = mxGetPr(plhs[0]);
    mxArray *TempDryFaceFlag = mxCreateDoubleMatrix(1, Ne, mxREAL);
    double *DryFaceFlag = mxGetPr(TempDryFaceFlag);
    
    evaluateWetDryInterface( regType, prhs[1], DryFaceFlag);
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (mwIndex k = 0; k < Ne; k++) {
        if (DryFaceFlag[k] == 1) {
            // the considered face is one between two adjacent dry cells, or one dry cell and another wet cell, for this situation, do nothing
            continue;
        }
        else{
            for(mwIndex i = 0; i < Nfp; i++){
                mwIndex j  = k*Nfp + i;
                if(hum[j] * nx[j] > 0 && -hup[j] * nx[j] <= 0)
                    fluxS[j] = hwm[j] * ( hum[j] * nx[j] )/hm[j];
                else if(hum[j] * nx[j] <= 0 && -hup[j] * nx[j] > 0)
                    fluxS[j] = hwp[j] * ( hup[j] * nx[j] )/hp[j];
                else
                    fluxS[j] = ( hwm[j] * hum[j]/hm[j] + hwp[j] * hup[j]/hp[j] ) * nx[j] / 2 ;
            }
        }
    }
    
    return;
}