# include "mex.h"
# include <math.h>
# include <malloc.h>
# include "../../@NdgAbstractNonhydrostaticSolver/private/mxNonhydroSWE.h"
#include "../../../../../Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.h"

#define NRHS 6
#define NLHS 1

/*
 * Purpose: This function is used to identify the wet-dry interface during the calculation
 *
 * Input:
 *      signed char[K] regType status of the mesh cell, used to identify a cell is wet or not
 * 		double[Nf x K]  EToE the mesh topological relation, used to study the status of adjacent cells
 *      double[2 x NOIE]  InnerEdgeFToE elements adjacent to a given inner edge, used to find the cells adjacent to a given inner edge
 * 		double[2 x NOBE]  BoundaryEdgeFToE elements adjacent to a given boundary edge, used to find the cells adjacent to a given boundary edge
 * 		signed char[NOBE]  BoundaryEdgeFtype face type of the boundary edge, used to impose boundary condition when an edge of the studied cell belong to the boundary
 *      double[2 x NOBE] BoundaryEdgeFToF the local face order of the studied boundary edge, used to identify whether the studied face of a cell adjacent to the boudary belongs to the boundary
 * Output:
 * 		double[(Nf + 1) x K_Interface] 	EdgeType the edge type of the studied mesh
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check input & output */
    if (nrhs != NRHS) {
        mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NRHS);
    }
    
    if (nlhs != NLHS) {
        mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NLHS);
    }
    /* Status of the studied cell of the given mesh, wet or dry*/
    signed char* regType = (signed char*)mxGetData(prhs[0]);
    /* Element to element conection relationship*/
    double *EToE = mxGetPr(prhs[1]);
    /*Number of faces for the studied mesh cell*/
    int Nf = (int)mxGetM(prhs[1]);
    double *InnerEdgeFToE = mxGetPr(prhs[2]);
    /*Number of inner edges*/
    mwSize NOIE = (mwSize)mxGetN(prhs[2]);
    double *BoundaryEdgeFToE = mxGetPr(prhs[3]);
    signed char* BoundaryEdgeFtype = (signed char*)mxGetData(prhs[4]);
    double *BoundaryEdgeFToF = mxGetPr(prhs[5]);
    /*Number of boundary edges*/
    mwSize NOBE = (mwSize)mxGetN(prhs[5]);
    
    mxArray *TempDryFaceFlag = mxCreateDoubleMatrix(1, NOIE, mxREAL);
    double *DryFaceFlag = mxGetPr(TempDryFaceFlag);
    
    evaluateWetDryInterface( regType, prhs[2], DryFaceFlag);
    /*Total number of elements at the wet-dry interface*/
    int total_element = 1;
    
    double *TempWetDryCell = (double *)mxCalloc((Nf + 1) * total_element, sizeof(double));
    double *ptrWetDryCell = NULL;
    
    for (mwIndex i = 0; i < NOIE; i++){
        /*This part is used to find the cells located at wet-dry interface*/
        if (DryFaceFlag[i] == 1){
            /*The local cell is wet, while the adjacent one is not*/
            if ( (NdgRegionType)regType[ (int)InnerEdgeFToE[2*i] - 1 ] == NdgRegionWet && (NdgRegionType)regType[ (int)InnerEdgeFToE[2*i + 1] - 1 ] != NdgRegionWet ){
                TempWetDryCell[(total_element-1) * (Nf + 1)] = InnerEdgeFToE[2*i];
            }
            /*The adjacent cell is wet, while the local one is not*/
            else if( (NdgRegionType)regType[ (int)InnerEdgeFToE[2*i + 1] - 1 ] == NdgRegionWet && (NdgRegionType)regType[ (int)InnerEdgeFToE[2*i] - 1 ] != NdgRegionWet ){
                TempWetDryCell[(total_element-1) * (Nf + 1)] = InnerEdgeFToE[2*i + 1];
            }
            else
            {
                /*The studied inner face is flaged as a face located besides a wet-dry interface, but this is not reflected in the status checking in this part*/
                continue;
            }
            /*Index of the cell located besides the wet-dry interface*/
            mwIndex ele = (mwIndex)TempWetDryCell[(total_element-1) * (Nf + 1)];
            /*This part is used to check the status of cell adjacent to the studied wet cell*/
            for(mwIndex j=0; j<Nf; j++){
                mwIndex AdjacentEle = (mwIndex)EToE[(ele - 1)*Nf + j];
                if ( (NdgRegionType)regType[ AdjacentEle - 1 ] == NdgRegionWet && AdjacentEle != ele )
                    /*if the adjacent cell is a wet cell and is not the cell itself, then the edge is an inner edge */
                    TempWetDryCell[(total_element-1) * (Nf + 1) + j + 1] = 0;
                else if ( AdjacentEle == ele ){
                    /*if the adjcent cell is the cell itself, the the cell is located besides
                     * the boundary, and the boudary type should be defined based on the boudary type
                     */
                    for (mwIndex k = 0 ; k < NOBE; k++ ){
                        /*if the cell adjacent to the studied boudary edge is the studied cell,
                         * and the face index belonging to the BoundaryEdgeFToF property is equal
                         * to the local face index, then the boundary type is set based on the
                         * BoundaryEdgeFtype property.
                         */
                        if ( BoundaryEdgeFToE[2*k] == AdjacentEle && BoundaryEdgeFToF[2*k] == j ){
                            NdgNonhydroEdgeType Boundarytype = (NdgNonhydroEdgeType)BoundaryEdgeFtype[k];
                            mxSetEdgeType( Boundarytype, TempWetDryCell + (total_element-1) * (Nf + 1) + j + 1);
                            break;
                        }
                    }
                }
                else if ( (NdgRegionType)regType[ AdjacentEle - 1  ] != NdgRegionWet && AdjacentEle != ele )
                    /*if the adjacent cell is not a wet cell, then the edge is an Dirichlet edge */
                    TempWetDryCell[(total_element-1) * (Nf + 1) + j + 1] = 2;
                
            }
            total_element++;
            /*Reallocate memory space for the TempWetDryCell data*/
            ptrWetDryCell = (double *)mxRealloc(TempWetDryCell, total_element * (Nf + 1) * sizeof(double));
            if (ptrWetDryCell == NULL )
            {
                mexPrintf("can't expand the table!\n");
                break;
            }
            else
            {
                // if succeded, reset the address of the enlarged space to the pointer establised at the beginning
                TempWetDryCell = ptrWetDryCell;
            }
        }
    }
    
    plhs[0] = mxCreateDoubleMatrix( Nf + 1 , total_element - 1 , mxREAL);
    double * WetToDryCell = mxGetPr(plhs[0]);
    for(mwIndex i = 0; i < total_element - 1; i++){
        WetToDryCell[i*(Nf+1)] = TempWetDryCell[i * (Nf + 1)];
        for (mwIndex j =0; j<Nf; j++){
            WetToDryCell[i*(Nf+1) + j + 1] = TempWetDryCell[i * (Nf + 1) + j + 1];
        }
    }
}