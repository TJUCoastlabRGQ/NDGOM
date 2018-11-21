# include "mex.h"
# include <math.h>
#include <malloc.h>

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
  double *Status = mxGetPr(prhs[0]);  
  double *EToE = mxGetPr(prhs[1]);
  // Nface = mxGetScalar(prhs[2]);
  int Nface = (int)mxGetScalar(prhs[2]);
  double *Fmask = mxGetPr(prhs[3]);
 // int Ncell = mxGetN(prhs[1]);
  mwSize Ncell = mxGetN(prhs[1]);
  int PointPerface = (int)mxGetScalar(prhs[4]);
  int Np = (int)mxGetScalar(prhs[5]);
  size_t current_element = 0;
  size_t total_element = 1;
  int WetNum = 0;
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
		  WetNum++;
		  for (mwIndex j = 0; j<Nface; j++)
		  {
			  int AdjacentCell = (int)EToE[i*Nface+j];
			  if (Status[AdjacentCell-1] == 4)
			  {   //if the adjacent cell of the wet cell is a wet cell, doing nothing
			       continue;
			  }
			  else
			  {
				  int dryface;
				  int flag = 0;  
				  for (mwIndex k = 0; k<Nface; k++)
				   { //to find the adjacent dry face of the corresponding adjacent dry cell				  
					  if ((int)EToE[(AdjacentCell - 1)*Nface + k] == i+1)
					  {
						  dryface = k + 1;
						  flag = 1;
					  }
				   }
				  if (flag == 0)
				   {
					  mexPrintf("error when look for the adjacent dry face!\n");
				   }
				if(current_element==total_element-1)
				{
                    total_element+=1;
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
						TempWetCellIndex = ptrWetCellIndex;
						TempWetDryFace = ptrWetDryFace;
						TempAdjacentDryCell = ptrAdjacentDryCell;
						TempAdjacentFace = ptrAdjacentFace;
						TempWetDryPoint = ptrWetDryPoint;
					}
                }
        TempWetCellIndex[current_element] = i+1;    //The wet cell index
		TempWetDryFace[current_element] = j+1;
		TempAdjacentDryCell[current_element] = AdjacentCell;
		TempAdjacentFace[current_element] = dryface;
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