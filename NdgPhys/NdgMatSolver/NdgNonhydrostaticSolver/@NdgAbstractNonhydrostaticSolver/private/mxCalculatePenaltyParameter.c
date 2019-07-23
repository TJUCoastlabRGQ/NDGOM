# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 7
#define NLHS 1

/*
 * Purpose: This function is used to get the global penalty paremeter for the IPDG method.
 *This function is programmed according to ( Shahbazi, 2005 ).
 *
 * Input:
 *       double[K] LAV the legth, area or volume of the studied one, two or three dimensional mesh.
 * 		double[NIe x 2]  InnerEdgeFToE face to element relation for the inner edge
 * 		double[NIe]  InnerEdgeLAV length or area of the inner edge for the studied one, two or three dimensional mesh, this value is set to 1 for the one dimensional mesh
 * 		double[NBe x 2]  BoundaryEdgeFToE face to element relation for the boundary edge
 * 		double[NIe]  BoundaryEdgeLAV length or area of the boundary edge for the studied one, two or three dimensional mesh, this value is set to 1 for the one dimensional mesh
 *      signed char[1] MeshType mesh type of the studied mesh, one, two or three according to the studied mesh
 *      double[1]   N the approximation order
 * Output:
 * 		double[1] 	Tau the global penalty paremeter
 */

typedef enum {
    One = 1,
    Two = 2,
    Three = 3
} Dimension;

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
    double * LAV = mxGetPr(prhs[0]);
    mwSize K = mxGetNumberOfElements(prhs[0]);
    double * InnerEdgeFToE = mxGetPr(prhs[1]);
    double * InnerEdgeLAV = mxGetPr(prhs[2]);
    mwSize NOIE = mxGetNumberOfElements(prhs[2]);
    double * BoundaryEdgeFToE = mxGetPr(prhs[3]);
    double * BoundaryEdgeLAV = mxGetPr(prhs[4]);
    mwSize NOBE = mxGetNumberOfElements(prhs[4]);
    signed char * MeshType = (signed char*)mxGetData(prhs[5]);
    double  N = mxGetScalar(prhs[6]);
    mxArray *PenaltyParameter = mxCreateDoubleMatrix(1, K, mxREAL);
    double *TempPenaltyParameter = mxGetPr(PenaltyParameter);
    /*Loop for all the elements adjacent to the inner edge*/
	for (mwIndex i = 0; i <  NOIE; i++){
        /*Length or Area of the studied inner edge*/
        double IELAV = InnerEdgeLAV[i];
        /*Local cell index adjacent to the inner edge*/
        mwIndex LocalEle = (mwIndex)InnerEdgeFToE[2*i];
		/*Adjacent cell index adjacent to the inner edge*/
		mwIndex AdjacentEle = (mwIndex)InnerEdgeFToE[2*i + 1];
        /*Length, area or volume of the local cell adjacent to the studied inner edge */
		double LocalCellLAV = LAV[LocalEle - 1];
		/*Length, area or volume of the adjacent cell adjacent to the studied inner edge */
		double AdjacentCellLAV = LAV[AdjacentEle - 1];
		TempPenaltyParameter[LocalEle - 1] += (IELAV / 2) / LocalCellLAV;
		TempPenaltyParameter[AdjacentEle - 1] += (IELAV / 2) / AdjacentCellLAV;
    }
    /*Loop for all the elements adjacent to the boundary edge*/
	for (mwIndex i = 0; i < NOBE; i++){
        /*Length or Area of the studied boundary edge*/
        double BELAV = BoundaryEdgeLAV[i];
        /*Cell adjacent to the boundary edge, for the boundary edge, both the local and adjacent cell are the same one*/
        mwIndex ele = (mwIndex)BoundaryEdgeFToE[2*i];
        /*Length, area or volume of the cell adjacent to the studied inner edge */
        double CellLAV = LAV[ele-1];
        TempPenaltyParameter[ele-1] += BELAV/CellLAV;
    }
    /*multiply the penalty parameter with coefficient related to the dimension and approximation order*/
    Dimension meshType = (Dimension)MeshType[0];
        switch(meshType){
            case One: 
				for (mwIndex i = 0; i < K; i++)
                    TempPenaltyParameter[i] *= (N+1) * (N+1)/1;
                break;
            case Two:
				for (mwIndex i = 0; i < K; i++)
                    TempPenaltyParameter[i] *= (N+1) * (N+2)/2;
                break;                
            case Three:
				for (mwIndex i = 0; i < K; i++)
                    TempPenaltyParameter[i] *= (N+1) * (N+3)/3;
                break;
        }
    /*find the global penalty parameter*/    
 double GlobalMaximum = 0;
 for (mwIndex i = 0; i < K; i++){
  if (TempPenaltyParameter[i] >= GlobalMaximum)
      GlobalMaximum = TempPenaltyParameter[i];
 }
  plhs[0] = mxCreateDoubleScalar( GlobalMaximum );  
 
}