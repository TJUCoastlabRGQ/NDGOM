#include "mex.h"
#include "../../../../../Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 6
#define NLHS 2

/*
* Purpose: This function is used to get the local or adjacent flux term about the physical variable
* appeared in the quadrature free two dimensional non-hydrostatic model in the downwind manner.
*
* Input:
*       signed char[Ne] status the cell status of the studied mesh, used to identify a cell is dry or wet, so as to determine the numerical flux at the wet-dry and dry-dry interface
* 		double[nfp x Ne]  variablexmp local or adjacent variable in x direction
* 		double[nfp x Ne]  variableymp local or adjacent variable in y direction
*       double[nfp x Ne]  nx the direction vector in x direction of the local cell 
* 		double[nfp x Ne]  ny the direction vector in y direction of the local cell
* Output:
* 		double[nfp x Ne] 	LocalOrAdjacentFluxX the local or adjacent flux term in x direction
* 		double[nfp x Ne] 	LocalOrAdjacentFluxY the local or adjacent flux term in Y direction
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
	double *variablexmp = mxGetPr(prhs[2]);  
	double *variableymp = mxGetPr(prhs[3]);
	double *nx = mxGetPr(prhs[4]);
	double *ny = mxGetPr(prhs[5]);
	
    mwSize Nfp = mxGetM(prhs[2]);
    mwSize Ne = mxGetN(prhs[2]);	
	
    plhs[0] = mxCreateDoubleMatrix(Nfp, Ne, mxREAL);
    double *LocalOrAdjacentFluxX = mxGetPr(plhs[0]);	
	
    plhs[1] = mxCreateDoubleMatrix(Nfp, Ne, mxREAL);
    double *LocalOrAdjacentFluxY = mxGetPr(plhs[1]);
	
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
                LocalOrAdjacentFluxX[j] = variablexmp[j] * nx[j];
                LocalOrAdjacentFluxY[j] = variableymp[j] * ny[j];
            }
        }
    }
    
    return;			    
}
