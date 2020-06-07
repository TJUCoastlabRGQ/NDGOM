# include "mex.h"
# include <math.h>
# include "mxNonhydroSWE.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 4
#define NLHS 1

/*
 * Purpose: This function is used to define the edge type of the computation mesh, such that different nonhydro-pressure related boundary condition can be considered.
 *Edge of the computation mesh is categarized into inner edge, Neumann edge and Dirichlet edge.
 *
 * Input:
 * 		double[Nf x K]  meshEToE all the elements that adjacent to a studied mesh
 *      double[NBe x 2]  BoundaryEdgeFToF order of face that located one the boundary
 * 		double[NBe x 2]  BoundaryEdgeFToE elements adjacent to a given boundary edge
 * 		signed char[NBe]  BoundaryEdgeFtype face type of the boundary edge
 * Output:
 * 		double[Nf x K] 	EdgeType the edge type of the studied mesh
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
    double * meshEToE = mxGetPr(prhs[0]);
    mwSize Nf = mxGetM(prhs[0]);
    mwSize K = mxGetN(prhs[0]);
    double * BoundaryEdgeFToF = mxGetPr(prhs[1]);
    double * BoundaryEdgeFToE = mxGetPr(prhs[2]);
    signed char * BoundaryEdgeFtype = (signed char*)mxGetPr(prhs[3]);
    mwSize NBe = mxGetNumberOfElements(prhs[3]);
    /*Face type of inner edge is set to 0 by default*/
    plhs[0] = mxCreateDoubleMatrix(Nf, K, mxREAL);
    double * EdgeType = mxGetPr(plhs[0]);
    /*Program to define the face type of the computation mesh at the boundary*/
    for(mwIndex i = 0; i< NBe; i++){
        mwIndex ele = (mwIndex)BoundaryEdgeFToE[2*i];
        mwIndex face = (mwIndex)BoundaryEdgeFToF[2*i];
        NdgNonhydroEdgeType Boundarytype = (NdgNonhydroEdgeType)BoundaryEdgeFtype[i];
        mxSetEdgeType( Boundarytype, EdgeType + (ele - 1)*Nf + face - 1);
    }
    
}