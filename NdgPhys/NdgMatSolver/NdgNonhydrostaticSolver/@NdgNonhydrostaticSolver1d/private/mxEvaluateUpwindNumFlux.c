#include "mex.h"
#include "../../../../../Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 7
#define NLHS 1

/*
* Purpose: This function is used to calculate the partial derivative about the physical variable
* appeared in the quadrature free two dimensional non-hydrostatic model in the upwind manner.
*
* Input:
*       signed char[Ne] status the cell status of the studied mesh, used to identify a cell is dry or wet, so as to determine the numerical flux at the wet-dry and dry-dry interface
*       double[nfp x Ne]  hum  local flux term in x direction, used together with variables hvm, hup and hvp to determine the upwind direction
* 		double[nfp x Ne]  hup  adjacent flux term in x direction, used together with variables hum, hvm and hvp to determine the upwind direction
* 		double[nfp x Ne]  variablexm local variable in x direction, used together with nx to determine the flux in x direction when the flux direction is outward
* 		double[nfp x Ne]  variablexp adjacent variable in x direction, used together with nx to determine the flux in x direction when the flux direction is inward 
*       double[nfp x Ne]  nx the direction vector in x direction of the local cell 
* Output:
* 		double[nfp x Ne] 	UpwindFluxX the upwind flux term in x direction
*/

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
    double *hum = mxGetPr(prhs[2]);
    double *hup = mxGetPr(prhs[3]);
	double *variablexm = mxGetPr(prhs[4]);
	double *variablexp = mxGetPr(prhs[5]);
    double *nx = mxGetPr(prhs[6]);
    
    mwSize Nfp = mxGetM(prhs[2]);
    mwSize Ne = mxGetN(prhs[2]);
	
    plhs[0] = mxCreateDoubleMatrix(Nfp, Ne, mxREAL);
    double *UpwindFluxX = mxGetPr(plhs[0]);
	
    mxArray *TempDryFaceFlag = mxCreateDoubleMatrix(1, Ne, mxREAL);
    double *DryFaceFlag = mxGetPr(TempDryFaceFlag);
    
    evaluateWetDryInterface( regType, prhs[1], DryFaceFlag);
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (mwIndex k = 0; k < Ne; k++) {
        if (DryFaceFlag[k] == 1) {
            /* the considered face is one between two adjacent dry cells, or one dry cell and another wet cell, for this situation, do nothing.
             * numerical flux at the wet-dry interface is set to zero for the current version */
            continue;
        }
        else{
            for(mwIndex i = 0; i < Nfp; i++){
                mwIndex j  = k*Nfp + i;
                if(hum[j] * nx[j]  > 0 && -hup[j] * nx[j]  <= 0){
                	UpwindFluxX[j] = variablexm[j]*nx[j];
				}
                else if(hum[j] * nx[j] <= 0 && -hup[j] * nx[j] > 0){
                	UpwindFluxX[j] = variablexp[j]*nx[j];
				}
				else{
                	UpwindFluxX[j] = ( variablexm[j] + variablexp[j]) *nx[j]/2;
				}
            }
        }
    }
    
    return;			    
    
}
