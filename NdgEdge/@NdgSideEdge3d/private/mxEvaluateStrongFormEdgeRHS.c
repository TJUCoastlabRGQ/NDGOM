#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"
#include "blas.h"
#include <stdlib.h>
#include "../../../NdgMath/NdgMath.h"

#if !defined(_WIN32)
#define dgemm dgemm_
#endif
// #if !defined(_WIN32)
// #define dgemm dgemm_
// #endif

#define DEBUG 0

#define NRHS 10
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
  double *FToF = mxGetPr(prhs[3]);
  double *FToN1 = mxGetPr(prhs[4]);
  double *FToN2 = mxGetPr(prhs[5]);
  double *Js = mxGetPr(prhs[6]);
  double *J = mxGetPr(prhs[7]);
  double *fluxM = mxGetPr(prhs[8]);
  double *fluxP = mxGetPr(prhs[9]);
  double *fluxS = mxGetPr(prhs[10]);
  int Nface = (int)mxGetScalar(prhs[11]);

  // dims = mxGetDimensions(prhs[6]);
  const int Np = mxGetM(prhs[7]);  // num of interp nodes
  const int K = mxGetN(prhs[7]);   // num of elements
  const mwSize *dims = mxGetDimensions(prhs[8]);
  const int Nfp = dims[0];
  const int Ne = dims[1];  // num of edges
  int Nfield;

  if (mxGetNumberOfDimensions(prhs[7]) > 2) {
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
  
double *ERHS = malloc(Np*K*Nfield*Nface*sizeof(double));
memset(ERHS, 0, Np*K*Nfield*Nface*sizeof(double));
double *TempFacialIntegral = malloc(Np*K*Nfield*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int face = 0; face < Ne; face++){
        for (int field = 0; field < Nfield; field++){
			StrongFormInnerEdgeRHS(face, FToE, FToF, Np, K, Nfp, FToN1, FToN2, fluxM + field*Ne*Nfp,\
				fluxP + field*Ne*Nfp, fluxS + field*Ne*Nfp, Js, Mb, ERHS + field*Np*K*Nface);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int k=0; k < K; k++){
        for(int field=0;field<Nfield;field++){
            for(int face=1;face<Nface;face++){
                Add( ERHS + field*Np*K*Nface+k*Np, ERHS + field*Np*K*Nface + k*Np, ERHS + field*Np*K*Nface + face*Np*K + k*Np, Np);
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nfield; field++){
			MultiEdgeContributionByLiftOperator(ERHS + field*Np*K*Nface + k*Np, TempFacialIntegral + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
            for(int n=0; n<Np;n++){
               frhs[field*Np*K+k*Np+n] = ERHS[field*Np*K*Nface + k*Np + n]; 
            }
            
		}
	}

free(ERHS);
free(TempFacialIntegral);
}