#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"
#include "matrix.h"
#include "blas.h"
#include <stdlib.h>
#include <string.h>


typedef enum {
	One = 1,
	Two = 2,
	Three = 3
} NdgMeshType;

void Minus(double *dest, double *sourcea, double *sourceb, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] - sourceb[i];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	mxArray *TempFToN1 = mxGetField(prhs[0],0,"FToN1");
	double *FToN1 = mxGetPr(TempFToN1);
	mxArray *TempFToN2 = mxGetField(prhs[0], 0, "FToN2");
	double *FToN2 = mxGetPr(TempFToN2);
	int Ne = mxGetScalar(mxGetField(prhs[0], 0, "Ne"));
//	double *dt = mxGetPr(prhs[1]);
	//double dt = mxGetScalar(prhs[1]);
	//plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	//double *data = mxGetPr(plhs[0]);
	//*dt = *dt + 1;
	int Nz;
    if ((Nz = mxGetFieldNumber(prhs[0],"Nz"))!=-1){
        printf("There is property Nz contained in the edge object!\n");
    }
    else{
		printf("No property named Nz is contained in the edge object!\n");
    }
    
    signed char *ftype = (signed char *) mxGetData(prhs[2]);
	NdgMeshType type = (NdgMeshType)ftype[0];
	if (type == Two){
		printf("This is a two-dimensional mesh!\n");
	}
	else if (type == Three){
		printf("This is a three-dimensional mesh!\n");
	}

	double *a = mxGetPr(prhs[3]);
	double *b = mxGetPr(prhs[4]);
	int N = mxGetNumberOfElements(prhs[3]);
	double *copya = malloc(N*sizeof(double));
	memcpy(copya, a, N*sizeof(double));
	Minus(copya, copya, b, N);
	for (int i = 0; i < N; i++)
		printf("Result is: %f\n", copya[i]);
}