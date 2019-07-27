# include "mex.h"
# include <math.h>
# include <malloc.h>
# include "mxNonhydroSWE.h"

#define NRHS 6
#define NLHS 5

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
  // Status of the studied cell of the given mesh, wet or dry
  double *Status = mxGetPr(prhs[0]);  
  // Element to element conection relationship
  double *EToE = mxGetPr(prhs[1]);
  // Number of face for the studied cell
  int Nface = (int)mxGetScalar(prhs[2]);
  // Boundary node order of the studied cell
  double *Fmask = mxGetPr(prhs[3]);
 // Element number of the given mesh
  mwSize Ncell = mxGetN(prhs[1]);
  // Number of point per face of the studied cell
  int PointPerface = (int)mxGetScalar(prhs[4]);
  // Number of interpolation point per element
  int Np = (int)mxGetScalar(prhs[5]);
  size_t current_element = 0;
  size_t total_element = 1;
  int WetNum = 0;
  // Allocate space for storing the outcome and return the correspoing pointer to the allocated address
  int *TempWetCellIndex = (int *)mxCalloc(total_element, sizeof(int));
  int *TempWetDryFace = (int *)mxCalloc(total_element, sizeof(int));
  int *TempAdjacentDryCell = (int *)mxCalloc(total_element, sizeof(int));
  int *TempAdjacentFace = (int *)mxCalloc(total_element, sizeof(int));
  int *TempWetDryPoint = (int *)mxCalloc(total_element*PointPerface, sizeof(int));
  int *ptrWetCellIndex = NULL, *ptrWetDryFace = NULL, *ptrAdjacentDryCell = NULL, \
      *ptrAdjacentFace = NULL, *ptrWetDryPoint = NULL;
  for(mwIndex i=0;i<Ncell;i++)
  {
	  if (Status[i] == 4)   //if the cell is a wet cell
	  {
		  WetNum++;        // Number of the wet cell increases
		  for (mwIndex j = 0; j<Nface; j++)
		  {
              // Index of the adjacent cell of the studied wet cell
			  int AdjacentCell = (int)EToE[i*Nface+j];  
			  if (Status[AdjacentCell-1] == 4)
			  {   //if the adjacent cell of the wet cell is a wet cell, doing nothing
			       continue;
			  }
              else // the adjacent cell of the wet cell is a dry cell
			  {
				  int dryface;
				  int flag = 0;  
				  for (mwIndex k = 0; k<Nface; k++)
				   { //to find the adjacent dry face of the corresponding adjacent dry cell				  
					  if ((int)EToE[(AdjacentCell - 1)*Nface + k] == i+1)  // the adjcact cell of the studied dre cell is the wet cell we studied
					  {
                          // the face index that adjacent to the studied wet cell of the studied dry cell
						  dryface = k + 1;
						  flag = 1;
					  }
				   }
				  if (flag == 0)
				   {
                      // None of the face of the studied dry cell adjacent to the studied wet cell, wrong!
					  mexPrintf("error when look for the adjacent dry face!\n");
				   }
				if(current_element==total_element-1)
				{
                    // enlarge space for the afore allocated space when the space is not enough
                    total_element+=1;
                    // enlarge the space for the afore allocated space and return this pointer to the allocated space to a new pointer
                    ptrWetCellIndex=(int *)mxRealloc(TempWetCellIndex,total_element*sizeof(int));
					ptrWetDryFace = (int *)mxRealloc(TempWetDryFace, total_element*sizeof(int));
					ptrAdjacentDryCell = (int *)mxRealloc(TempAdjacentDryCell, total_element*sizeof(int));
					ptrAdjacentFace = (int *)mxRealloc(TempAdjacentFace, total_element*sizeof(int));
					ptrWetDryPoint = (int *)mxRealloc(TempWetDryPoint, total_element*PointPerface*sizeof(int));
					if (ptrWetCellIndex == NULL || ptrWetDryFace == NULL || ptrAdjacentDryCell == NULL\
						|| ptrAdjacentFace == NULL || ptrWetDryPoint == NULL)
					{
							mexPrintf("can't expand the table!\n");
					}
					else
					{
                        // if succeded, reset the address of the enlarged space to the pointer establised at the beginning  
						TempWetCellIndex = ptrWetCellIndex;
						TempWetDryFace = ptrWetDryFace;
						TempAdjacentDryCell = ptrAdjacentDryCell;
						TempAdjacentFace = ptrAdjacentFace;
						TempWetDryPoint = ptrWetDryPoint;
					}
                }
        // the studied wet cell index stored
        TempWetCellIndex[current_element] = i+1;
        // the face of the studied wet cell that adjacent to the studied dry cell 
		TempWetDryFace[current_element] = j+1;
        // the dry cell that adjacent to the studied wet cell
		TempAdjacentDryCell[current_element] = AdjacentCell;
        // the face of the studied dry cell that adjacent to the studied wet cell 
		TempAdjacentFace[current_element] = dryface;
        // the global point index of the point belonging to the wet cell located at the wet and dry interface
		for (mwIndex k = 0; k<PointPerface; k++)
			TempWetDryPoint[current_element*PointPerface + k] = i*Np + (int)Fmask[j*PointPerface+k];
		current_element++;
			  }				  
		  }
	  }  
  }
  
  double *ZeroFluxBoundary;
  double *AdjacentDryCellAndFace;
  double *WetDryPoint;

  plhs[0] = mxCreateDoubleMatrix(current_element, 2, mxREAL);
  ZeroFluxBoundary = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(current_element, 2, mxREAL);
  AdjacentDryCellAndFace = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleMatrix(current_element * PointPerface, 1, mxREAL);
  WetDryPoint = mxGetPr(plhs[2]);
  plhs[3] = mxCreateDoubleScalar(current_element);
  plhs[4] = mxCreateDoubleScalar(WetNum);
  for (mwIndex i = 0; i < current_element; i++)
  {
	  ZeroFluxBoundary[i] = TempWetCellIndex[i];
	  ZeroFluxBoundary[i + current_element] = TempWetDryFace[i];
	  AdjacentDryCellAndFace[i] = TempAdjacentDryCell[i];
	  AdjacentDryCellAndFace[i + current_element] = TempAdjacentFace[i];
	  for (mwIndex j = 0; j < PointPerface; j++)
		  WetDryPoint[i*PointPerface + j] = TempWetDryPoint[i*PointPerface + j];
  } 
  mxFree(TempWetCellIndex);
  TempWetCellIndex = NULL;
  mxFree(TempWetDryFace);
  TempWetDryFace = NULL;
  mxFree(TempAdjacentDryCell);
  TempAdjacentDryCell = NULL;
  mxFree(TempAdjacentFace);
  TempAdjacentFace = NULL;
  mxFree(TempWetDryPoint);
  TempWetDryPoint = NULL;
}