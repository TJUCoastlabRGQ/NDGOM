#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 12
#define NLHS 3

/*
 * Purpose: This function is used to get the normal face direction vector at the studied
 * element.
 *
 * Input:
 *       double[1]  ele  index of the studied local element
 * 		double[1]  adjacentEle index of the element adjacent to the studied element
 *       double[1] face face index of the studied local element
 *       double[2 x NOBE ]   BoundaryEdge_FToF the face to face topological relations for boundary edge
 * 		double[nfp x NOBE]  BoundaryEdge_nx  the normal face direction vector of boundary edge in x direction
 * 		double[nfp x NOBE]  BoundaryEdge_ny  the normal face direction vector of boundary edge in y direction
 *       double[2 x NOBE ]   BoundaryEdge_FToE the face to element topological relations for boundary edge
 * 		double[nfp x NOIE]  InnerEdge_nx the normal face direction vector of inner edge in x direction
 * 		double[nfp x NOIE]  InnerEdge_ny the normal face direction vector of inner edge in y direction
 *       double[2 x NOIE ]   InnerEdge_Js the face jacobian for inner edge
 *      double[2 x NOBE ]   BoundaryEdge_Js the face jacobian for boundary edge
 *
 * Output:
 * 		double[1] 	nx the normal face direction vector of the studied local element in x direction
 * 		double[1] 	ny the normal face direction vector of the studied local element in y direction
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
    
    //double *FToE = mxGetPr(prhs[1]);
    double ele = mxGetScalar(prhs[0]);
    double adjacentEle = mxGetScalar(prhs[1]);
    double face = mxGetScalar(prhs[2]);
    double *BoundaryEdge_FToF = mxGetPr(prhs[3]);
    double *BoundaryEdge_nx = mxGetPr(prhs[4]);
    double *BoundaryEdge_ny = mxGetPr(prhs[5]);
    double *BoundaryEdge_FToE = mxGetPr(prhs[6]);
    double *InnerEdge_nx = mxGetPr(prhs[7]);
    double *InnerEdge_ny = mxGetPr(prhs[8]);
    double *InnerEdge_FToE = mxGetPr(prhs[9]);
    double *BoundaryEdge_Js = mxGetPr(prhs[10]);
    double *InnerEdge_Js = mxGetPr(prhs[11]);
    
    mwSize BoundaryEdge_Nfp = mxGetM(prhs[4]);
    mwSize BoundaryEdge_Ne = mxGetN(prhs[4]);
    
    mwSize InnerEdge_Nfp = mxGetM(prhs[7]);
    mwSize InnerEdge_Ne = mxGetN(prhs[7]);
    
    if ( ele == adjacentEle ){ /*the studied element is adjacent to the boundary edge*/
        for(mwIndex i = 0; i < BoundaryEdge_Ne; i++){
            if ( ele == BoundaryEdge_FToE[2*i] && face == BoundaryEdge_FToF[2*i] ){
				plhs[0] = mxCreateDoubleScalar(BoundaryEdge_nx[BoundaryEdge_Nfp*i]);
				plhs[1] = mxCreateDoubleScalar(BoundaryEdge_ny[BoundaryEdge_Nfp*i]);
                plhs[2] = mxCreateDoubleScalar(BoundaryEdge_Js[BoundaryEdge_Nfp*i]);
                break;
            }
        }
    } 
    else{
        for(mwIndex i = 0; i < InnerEdge_Ne; i++ ){
            mwIndex ele1 = (mwIndex)InnerEdge_FToE[2*i];
            mwIndex ele2 = (mwIndex)InnerEdge_FToE[2*i+1];
            if( ele == ele1 && adjacentEle == ele2 ){
			  plhs[0] = mxCreateDoubleScalar(InnerEdge_nx[InnerEdge_Nfp*i]);
			  plhs[1] = mxCreateDoubleScalar(InnerEdge_ny[InnerEdge_Nfp*i]);
              plhs[2] = mxCreateDoubleScalar(InnerEdge_Js[InnerEdge_Nfp*i]);
              break;
            }
            else if ( ele == ele2 && adjacentEle == ele1 ){
			  plhs[0] = mxCreateDoubleScalar(-1 * InnerEdge_nx[InnerEdge_Nfp*i]);
			  plhs[1] = mxCreateDoubleScalar(-1 * InnerEdge_ny[InnerEdge_Nfp*i]);
              plhs[2] = mxCreateDoubleScalar(InnerEdge_Js[InnerEdge_Nfp*i]);
              break;                
            }
        }
    }
        return;
}
