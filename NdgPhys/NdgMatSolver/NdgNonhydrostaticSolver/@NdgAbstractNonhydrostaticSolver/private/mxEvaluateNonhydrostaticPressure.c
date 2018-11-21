#if !defined(_WIN32)
#define dgetrf dgetrf_
#define dgetrs dgetrs_
#endif

#include <math.h>
#include "mex.h"
#include "lapack.h"

#define NRHS 2
#define NLHS 1

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

  char  *chn = "N";
  ptrdiff_t row, col, lda, ldb, Nrhs;
  ptrdiff_t *ipiv;
  ptrdiff_t info;
  row = mxGetM(prhs[0]);
  col = mxGetN(prhs[0]); 
  lda = row;
  ldb = row;
  Nrhs = mxGetN(prhs[1]);
  ipiv = ( ptrdiff_t* )malloc( sizeof(ptrdiff_t) * row );
  double *pA = mxGetPr(prhs[0]);
  double *pb = mxGetPr(prhs[1]);
  double *pa;
  plhs[0] = mxCreateDoubleMatrix((int)row, 1, mxREAL);
  pa = mxGetPr(plhs[0]);
  				
  dgetrf(&row, &col, pA, &lda, ipiv, &info);

  dgetrs(chn, &col, &Nrhs, pA, &lda, ipiv, pb, &ldb, &info);

  for(int i=0; i<col; i++)
	  pa[i] = pb[i];
  free(ipiv);
  
}