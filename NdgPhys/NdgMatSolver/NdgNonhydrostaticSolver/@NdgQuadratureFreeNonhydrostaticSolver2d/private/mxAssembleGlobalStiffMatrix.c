# include "mex.h"
# include <math.h>
# include <string.h>
# include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 14
#define NLHS 1

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#define IsZero(d) ((d)==0)
#endif
/*The following function is used to find the order of the input parameter aim in matrix a*/
int bin(int aim, int low, int high, double *a)
{
	int mid;
	while (low <= high)
	{
		mid = (low + high) / 2;
		if ((int)a[mid] == aim)return mid;
		else if ((int)a[mid]<aim)low = mid + 1;
		else high = mid - 1;
	}
	return -1;//Not Find
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check input & output */
  if (nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NLHS);
  }
  
  double dt = mxGetScalar(prhs[0]);
  double *height = mxGetPr(prhs[1]);
  
  double *TempSecondOrderTerm = mxGetPr(prhs[2]);
  mwIndex *jcTempSecondOrderTerm = mxGetJc(prhs[2]);
  mwIndex *irTempSecondOrderTerm = mxGetIr(prhs[2]);
  
  
  double *TempPNPX = mxGetPr(prhs[3]);
  mwIndex *jcTempPNPX = mxGetJc(prhs[3]);
  mwIndex *irTempPNPX = mxGetIr(prhs[3]);
  
  double *fhx = mxGetPr(prhs[4]);
    
  double *TempPNPY = mxGetPr(prhs[5]);
  mwIndex *jcTempPNPY = mxGetJc(prhs[5]);
  mwIndex *irTempPNPY = mxGetIr(prhs[5]);
  
  double *fhy = mxGetPr(prhs[6]);
  
  double *NP = mxGetPr(prhs[7]);
  mwIndex *jcNp = mxGetJc(prhs[7]);
  mwIndex *irNp = mxGetIr(prhs[7]);
  
  double *H2Bx = mxGetPr(prhs[8]);
  
  double *H2By = mxGetPr(prhs[9]);
  
  double *HBxSquare = mxGetPr(prhs[10]);
  
  double *HBySquare = mxGetPr(prhs[11]);
  
  int row,col;
  
  row = (int)mxGetM(prhs[7]);
  col = (int)mxGetN(prhs[7]);
  
  double *sr;
  mwIndex *irs,*jcs;

  double *JcStiffMatrix = mxGetPr(prhs[12]);
  double *JrStiffMatrix = mxGetPr(prhs[13]);

  plhs[0] = mxCreateSparse(row, col, JcStiffMatrix[col], mxREAL);
  sr = mxGetPr(plhs[0]);
  irs = mxGetIr(plhs[0]);
  jcs = mxGetJc(plhs[0]);

//printf("The total threads out of the parallel domain is:%d\n",omp_get_num_threads());

 #ifdef _OPENMP
 #pragma omp parallel for num_threads(omp_get_max_threads())
 #endif
	  
  for (int i = 0; i < col; i++)
  {
   //   printf("The total threads inside the parallel domain is:%d\n",omp_get_num_threads());
   //   printf("The current thread number is:%d\n",omp_get_thread_num());
	  int Num = (int)JcStiffMatrix[i + 1] - (int)JcStiffMatrix[i];
	  double *temprhsu = sr + (int)JcStiffMatrix[i];
	  memset(temprhsu, 0, Num*sizeof(double));
	  /*How many non-zero Np included in column i*/
	  int NumNp = (int)jcNp[i + 1] - (int)jcNp[i];
	  for (int j = 0; j < NumNp; j++){
		  /*The position of the studied NP element in the input sparse matrix corresponding to non-hydrostatic pressure*/
		  int Index = (int)jcNp[i] + j;
		  /*Index of the row that contains the studied non-zero Np in the studied column*/
		  int rowIndex = (int)irNp[Index];
		  /*Return the position of the studied NP element in the column of the output stiff matrix.
		  Here rowIndex may not equal to Order, since the final output stiff matrix is decided by 
		  several part.*/
		  int Order = bin(rowIndex, 0, NumNp - 1, JrStiffMatrix + (int)JcStiffMatrix[i]);
		  temprhsu[Order] = dt * NP[Index] * \
			  (height[rowIndex] * (H2Bx[rowIndex] + H2By[rowIndex]) - \
			  (HBxSquare[rowIndex] + HBySquare[rowIndex] + 4));
	  }

	  int NumSecondOrderTerm = (int)jcTempSecondOrderTerm[i + 1] - (int)jcTempSecondOrderTerm[i];
	  for (int j = 0; j < NumSecondOrderTerm; j++){
		  int rowIndex = (int)irTempSecondOrderTerm[j + (int)jcTempSecondOrderTerm[i]];
		  int Order = bin(rowIndex, 0, NumSecondOrderTerm - 1, JrStiffMatrix + (int)JcStiffMatrix[i]);
		  int Index = (int)jcTempSecondOrderTerm[i] + j;
		  temprhsu[Order] = temprhsu[Order] + dt * height[rowIndex] * (height[rowIndex] * TempSecondOrderTerm[Index]);
	  }

	  int NumPNPX = jcTempPNPX[i + 1] - jcTempPNPX[i];
	  for (int j = 0; j < NumPNPX; j++){
		  int rowIndex = (int)irTempPNPX[j + (int)jcTempPNPX[i]];
		  int Order = bin(rowIndex, 0, NumPNPX - 1, JrStiffMatrix + (int)JcStiffMatrix[i]);
		  int Index = (int)jcTempPNPX[i] + j;
		  temprhsu[Order] = temprhsu[Order] + dt * height[rowIndex] * (TempPNPX[Index] * fhx[rowIndex]);
	  }

	  int NumPNPY = jcTempPNPY[i + 1] - jcTempPNPY[i];
	  for (int j = 0; j < NumPNPY; j++){
		  int rowIndex = (int)irTempPNPY[j + (int)jcTempPNPY[i]];
		  int Order = bin(rowIndex, 0, NumPNPY - 1, JrStiffMatrix + (int)JcStiffMatrix[i]);
		  int Index = (int)jcTempPNPY[i] + j;
		  temprhsu[Order] = temprhsu[Order] + dt * height[rowIndex] * (TempPNPY[Index] * fhy[rowIndex]);
	  }
  }
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(omp_get_max_threads())
 #endif  
  for (int i = 0; i < col + 1; i++)
	  jcs[i] = JcStiffMatrix[i];

 #ifdef _OPENMP
 #pragma omp parallel for num_threads(omp_get_max_threads())
 #endif  
  for (int i = 0; i < (int)JcStiffMatrix[col]; i++)
	  irs[i] = JrStiffMatrix[i];
  return;
}
