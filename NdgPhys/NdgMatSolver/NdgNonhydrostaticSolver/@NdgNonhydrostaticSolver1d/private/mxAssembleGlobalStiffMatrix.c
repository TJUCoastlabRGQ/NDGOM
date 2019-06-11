# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 10
#define NLHS 1

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#define IsZero(d) ((d)==0)
#endif


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
  
  double *TempSPNPX = mxGetPr(prhs[2]);
  mwIndex *jcTempSPNPX = mxGetJc(prhs[2]);
  mwIndex *irTempSPNPX = mxGetIr(prhs[2]);
  
//   double *TempSPNPY = mxGetPr(prhs[3]);
//   mwIndex *jcTempSPNPY = mxGetJc(prhs[3]);
//   mwIndex *irTempSPNPY = mxGetIr(prhs[3]);  
  
  double *TempPNPX = mxGetPr(prhs[3]);
  mwIndex *jcTempPNPX = mxGetJc(prhs[3]);
  mwIndex *irTempPNPX = mxGetIr(prhs[3]);
  
  double *fhx = mxGetPr(prhs[4]);
    
//   double *TempPNPY = mxGetPr(prhs[6]);
//   mwIndex *jcTempPNPY = mxGetJc(prhs[6]);
//   mwIndex *irTempPNPY = mxGetIr(prhs[6]);
  
//   double *fhy = mxGetPr(prhs[7]);
  
  double *NP = mxGetPr(prhs[5]);
  mwIndex *jcNp = mxGetJc(prhs[5]);
  mwIndex *irNp = mxGetIr(prhs[5]);
  
  double *H2Bx = mxGetPr(prhs[6]);
  
//   double *H2By = mxGetPr(prhs[10]);
  
  double *HBxSquare = mxGetPr(prhs[7]);
  
//   double *HBySquare = mxGetPr(prhs[12]);
  
  mwSize row,col;
  
  row = mxGetM(prhs[5]);
  col = mxGetN(prhs[5]);
  
  double *pi,*sr;
  mwIndex *irs,*jcs;

  double *JcStiffMatrix = mxGetPr(prhs[8]);
  double *JrStiffMatrix = mxGetPr(prhs[9]);

  plhs[0] = mxCreateSparse(row, col, JcStiffMatrix[col], 0);
  sr = mxGetPr(plhs[0]);
  irs = mxGetIr(plhs[0]);
  jcs = mxGetJc(plhs[0]);

  int cmplx;


 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif
	  
  for (mwIndex i = 0; i < col; i++)
  {
	  mxArray *tempdata = mxCreateDoubleMatrix(row, 1, mxREAL);
	  double *temprhsu = mxGetPr(tempdata);

	  for (mwIndex j = jcNp[i]; j<jcNp[i+1] && jcNp[i+1]-jcNp[i]>0; j++)
	  {
		  size_t rowIndex = irNp[j];
		  temprhsu[rowIndex] =  dt *NP[j]*\
              ( H2Bx[rowIndex]  - 1 / height[rowIndex]\
                  * ( HBxSquare[rowIndex] + 2 ) );
	  }

	  for (mwIndex j = jcTempSPNPX[i]; j<jcTempSPNPX[i + 1] && jcTempSPNPX[i + 1] - jcTempSPNPX[i]>0; j++)
	  {
		  size_t rowIndex = irTempSPNPX[j];
		  temprhsu[rowIndex] = temprhsu[rowIndex] + dt * ( height[rowIndex] * TempSPNPX[j] );
	  }

// 	  for (mwIndex j = jcTempSPNPY[i]; j<jcTempSPNPY[i + 1] && jcTempSPNPY[i + 1]-jcTempSPNPY[i]>0; j++)
// 	  {
// 		  size_t rowIndex = irTempSPNPY[j];
// 		  temprhsu[rowIndex] = temprhsu[rowIndex] + dt  *( height[rowIndex] * TempSPNPY[j] );
// 	  }

	  for (mwIndex j = jcTempPNPX[i]; j<jcTempPNPX[i + 1] && jcTempPNPX[i + 1]-jcTempPNPX[i]>0; j++)
	  {
		  size_t rowIndex = irTempPNPX[j];
		  temprhsu[rowIndex] = temprhsu[rowIndex] +  dt *( TempPNPX[j] * fhx[rowIndex] );
	  }
 
// 	  for (mwIndex j = jcTempPNPY[i]; j<jcTempPNPY[i + 1] && jcTempPNPY[i + 1] - jcTempPNPY[i]>0; j++)
// 	  {
// 		  size_t rowIndex = irTempPNPY[j];
// 		  temprhsu[rowIndex] = temprhsu[rowIndex] +  dt  *( TempPNPY[j] * fhy[rowIndex] );
// 	  }
      
	    pi = mxGetPi(tempdata);
        cmplx = (pi==NULL ? 0 : 1);
		if (cmplx)
		{
			mexPrintf("%Complex number detected, problematic!.\n" );
		}
		for (mwIndex rowflag = 0; rowflag<JcStiffMatrix[i + 1] - JcStiffMatrix[i]; rowflag++)
        {
			mwIndex index = rowflag + JcStiffMatrix[i];
			sr[index] = temprhsu[(mwIndex)JrStiffMatrix[index]];
			
        }
		mxDestroyArray(tempdata);
		temprhsu = NULL;
  }
 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif  
  for (mwIndex i = 0; i < col + 1; i++)
	  jcs[i] = JcStiffMatrix[i];

 #ifdef _OPENMP
 #pragma omp parallel for num_threads(DG_THREADS)
 #endif  
  for (mwIndex i = 0; i < (int)JcStiffMatrix[col]; i++)
	  irs[i] = JrStiffMatrix[i];
  return;
}



