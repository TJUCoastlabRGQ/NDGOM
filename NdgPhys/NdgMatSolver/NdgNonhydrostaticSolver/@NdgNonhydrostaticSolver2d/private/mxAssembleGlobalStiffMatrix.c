# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 16
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
  double rho = mxGetScalar(prhs[1]);
  double *NP = mxGetPr(prhs[2]);
  mwIndex *jcNp = mxGetJc(prhs[2]);
  mwIndex *irNp = mxGetIr(prhs[2]);

  double *height = mxGetPr(prhs[3]);

  double *TempPNPX = mxGetPr(prhs[4]);
  mwIndex *jcTempPNPX = mxGetJc(prhs[4]);
  mwIndex *irTempPNPX = mxGetIr(prhs[4]);

  double *bx = mxGetPr(prhs[5]);

  double *TempPNPY = mxGetPr(prhs[6]);
  mwIndex *jcTempPNPY = mxGetJc(prhs[6]);
  mwIndex *irTempPNPY = mxGetIr(prhs[6]);

  double *by = mxGetPr(prhs[7]);

  double *TempSPNPX = mxGetPr(prhs[8]);
  mwIndex *jcTempSPNPX = mxGetJc(prhs[8]);
  mwIndex *irTempSPNPX = mxGetIr(prhs[8]);

  double *TempSPNPY = mxGetPr(prhs[9]);
  mwIndex *jcTempSPNPY = mxGetJc(prhs[9]);
  mwIndex *irTempSPNPY = mxGetIr(prhs[9]);

  double *TempFNPBX = mxGetPr(prhs[10]);
  mwIndex *jcTempFNPBX = mxGetJc(prhs[10]);
  mwIndex *irTempFNPBX = mxGetIr(prhs[10]);

  double *TempFNPBY = mxGetPr(prhs[11]);
  mwIndex *jcTempFNPBY = mxGetJc(prhs[11]);
  mwIndex *irTempFNPBY = mxGetIr(prhs[11]);

  double *fhx = mxGetPr(prhs[12]);
  double *fhy = mxGetPr(prhs[13]);
  
  mwSize row,col;
  
  row = mxGetM(prhs[2]);
  col = mxGetN(prhs[2]);
  
  double *pi,*sr;
  mwIndex *irs,*jcs;

  double *JcStiffMatrix = mxGetPr(prhs[14]);
  double *JrStiffMatrix = mxGetPr(prhs[15]);

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
		  temprhsu[rowIndex] = 2 * dt / rho*NP[j]\
			  + 2 * dt / rho *NP[j ] * bx[rowIndex] * bx[rowIndex]\
			  + 2 * dt / rho *NP[j] * by[rowIndex] * by[rowIndex]\
			  + dt / rho *(fhx[rowIndex] * NP[j] * bx[rowIndex])\
			  + dt / rho *(fhy[rowIndex] * NP[j] * by[rowIndex]);
	  }

	  for (mwIndex j = jcTempPNPX[i]; j<jcTempPNPX[i + 1] && jcTempPNPX[i + 1]-jcTempPNPX[i]>0; j++)
	  {
		  size_t rowIndex = irTempPNPX[j];
		  temprhsu[rowIndex] = temprhsu[rowIndex] + 2 * dt / rho *(height[i] * TempPNPX[j]) * bx[rowIndex]\
			  + dt / rho*(fhx[rowIndex] * (height[i] * TempPNPX[j]));
	  }
 
	  for (mwIndex j = jcTempPNPY[i]; j<jcTempPNPY[i + 1] && jcTempPNPY[i + 1] - jcTempPNPY[i]>0; j++)
	  {
		  size_t rowIndex = irTempPNPY[j];
		  temprhsu[rowIndex] = temprhsu[rowIndex] + 2 * dt / rho *(height[i] * TempPNPY[j]) * by[rowIndex]\
			  + dt / rho*(fhy[rowIndex] * (height[i] * TempPNPY[j]));
	  }

	  for (mwIndex j = jcTempSPNPX[i]; j<jcTempSPNPX[i + 1] && jcTempSPNPX[i + 1] - jcTempSPNPX[i]>0; j++)
	  {
		  size_t rowIndex = irTempSPNPX[j];
		  temprhsu[rowIndex] = temprhsu[rowIndex] - dt / rho *(height[i] * TempSPNPX[j]) * height[rowIndex];
	  }

	  for (mwIndex j = jcTempSPNPY[i]; j<jcTempSPNPY[i + 1] && jcTempSPNPY[i + 1]-jcTempSPNPY[i]>0; j++)
	  {
		  size_t rowIndex = irTempSPNPY[j];
		  temprhsu[rowIndex] = temprhsu[rowIndex] - dt / rho *(height[i] * TempSPNPY[j]) * height[rowIndex];
	  }

	  for (mwIndex j = jcTempFNPBX[i]; j<jcTempFNPBX[i + 1] && jcTempFNPBX[i + 1]-jcTempFNPBX[i]>0; j++)
	  {
		  size_t rowIndex = irTempFNPBX[j];
		  temprhsu[rowIndex] = temprhsu[rowIndex] - dt / rho *(height[i] * TempFNPBX[j]);
	  }

	  for (mwIndex j = jcTempFNPBY[i]; j<jcTempFNPBY[i + 1] && jcTempFNPBY[i + 1] - jcTempFNPBY[i]>0; j++)
	  {
		  size_t rowIndex = irTempFNPBY[j];
		  temprhsu[rowIndex] = temprhsu[rowIndex] - dt / rho *(height[i] * TempFNPBY[j]);
	  }
      
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



