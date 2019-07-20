# include "mex.h"
# include <math.h>

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

typedef enum {
        Inner = 0,
        GaussEdge = 1,
        SlipWall = 2,
        NonSlipWall = 3,
        ZeroGrad = 4,
        Clamped = 5,
        ClampedDepth = 6,
        ClampedVel = 7,
        Flather = 8,
        NonLinearFlatherDepth = 9,
        NonLinearFlatherFlow = 10,
        NonReflectFlux = 11,
        BottomBoundary  = 12,
        UpperSurfaceBoundary =13
} NdgEdgeType;

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
       NdgEdgeType Boundarytype = (NdgEdgeType)BoundaryEdgeFtype[i];
       switch(Boundarytype){
           /*Face type at the wall boundary is set to be the Neumann boundary, and is flagged by 1*/
	      case SlipWall:
               EdgeType[(ele - 1)*Nf + face - 1] = 1;
               break;
	      case NonSlipWall:
               EdgeType[(ele - 1)*Nf + face - 1] = 1;
               break;               
           /*Face type at other boundary is set to be the homogeneous Dirichlet boundary, and is flagged by 2*/
           default:
               EdgeType[(ele - 1)*Nf + face - 1] = 2;
       }
    }
    
}