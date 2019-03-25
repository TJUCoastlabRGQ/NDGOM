
#include "mex.h"

// #if !defined(_WIN32)
// #define dgemm dgemm_
// #endif

#define DEBUG 0

#define NRHS 1
#define NLHS 1

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

  double *a = mxGetPr(prhs[0]);

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {1, 1, 1};

  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  double *frhs = mxGetPr(plhs[0]);
}