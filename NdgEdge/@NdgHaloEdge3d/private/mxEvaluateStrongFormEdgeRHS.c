#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"
#include "blas.h"
#include "..\..\..\NdgMath\NdgMath.h"

// #if !defined(_WIN32)
// #define dgemm dgemm_
// #endif

#define DEBUG 0

#define NRHS 8
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

  double *invM = mxGetPr(prhs[0]);
  double *Mb = mxGetPr(prhs[1]);
  double *FToE = mxGetPr(prhs[2]);
  double *FToN1 = mxGetPr(prhs[3]);
  double *Js = mxGetPr(prhs[4]);
  double *J = mxGetPr(prhs[5]);
  double *fluxM = mxGetPr(prhs[6]);
  double *fluxS = mxGetPr(prhs[7]);

  // dims = mxGetDimensions(prhs[6]);
  const int Np = mxGetM(prhs[5]);  // num of interp nodes
  const int K = mxGetN(prhs[5]);   // num of elements
  const mwSize *dims = mxGetDimensions(prhs[6]);
  const int Nfp = dims[0];
  const int Ne = dims[1];  // num of edges
  int Nfield;

  if (mxGetNumberOfDimensions(prhs[6]) > 2) {
    Nfield = dims[2];
  } else {
    Nfield = 1;  // fluxM is a 2D matrix
  }

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {Np, K, Nfield};

  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  double *frhs = mxGetPr(plhs[0]);

  char *chn = "N";
  double one = 1.0, zero = 0.0;
  ptrdiff_t oneI = 1;
  ptrdiff_t KI = K;
  ptrdiff_t np = Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(12)
#endif
  for (int fld = 0; fld < Nfield; fld++) {
    double *rhs = frhs + Np * K * fld;
    double *fluxM_ = fluxM + Nfp * Ne * fld;
    // double *fluxP_ = fluxP + Nfp * Ne * fld;
    double *fluxS_ = fluxS + Nfp * Ne * fld;

	for (int k = 0; k < Ne; k++){
		StrongFormBoundaryEdgeRHS(k, FToE, Np, Nfp, FToN1, fluxM_, fluxS_, Js, Mb, rhs);
	}

	double *temp = malloc(Np*K*sizeof(double));
	dgemm(chn, chn, &np, &KI, &np, &one, invM, &np, rhs, &np, &zero, temp,
		&np);

	for (int k = 0; k < K; k++){
		DotDivide(rhs + k*Np, temp + k*Np, J + k*Np, Np);
	}
	free(temp);
  }
  return;
}