#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"
#include "blas.h"
#include "matrix.h"
#include "/../../../../../NdgMath/NdgMath.c"
#include <stdlib.h>


#if !defined(_WIN32)
#define dgemm dgemm_
#endif

typedef enum {
	One = 1,
	Two = 2,
	Three = 3
} NdgMeshType;

void volumeIntegral(){

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//Order of the input variable is hcrit, meshUnion.type, prantl, InnerEdge, BoundaryEdge, nv, frhs, fphys, varIndex, cell, mesh, BoundaryEdgefp
	//BoundaryEdgefm, InnerEdgefm and InnerEdgefp have been deleted, because we can fetch them in part.
	int Np, K, Nvar;
	const size_t *PRHS;
	PRHS = mxGetDimensions(prhs[6]);
	Np = (int)PRHS[0];
	K = (int)PRHS[1];
	Nvar = (int)PRHS[2];
	double *Dr = mxGetField(prhs[9], 0, "Dr");
	double *Ds = mxGetField(prhs[9], 0, "Ds");
	double *rx = mxGetField(prhs[10], 0, "rx");
	double *sx = mxGetField(prhs[10], 0, "sx");
	double *ry = mxGetField(prhs[10], 0, "ry");
	double *sy = mxGetField(prhs[10], 0, "sy");
	double *J = mxGetFiled(prhs[10], 0, "J");
	double *varIndex = mxGetPr(prhs[8]);
	double *fphys = mxGetPr(prhs[7]);
	double *Var = malloc(Np*K*(Nvar - 1)*sizeof(double));
	double *Qx  = malloc(Np*K*(Nvar - 1)*sizeof(double));
	double *TempQx = malloc(Np*K*(Nvar - 1)*sizeof(double));
	double *Qy  = malloc(Np*K*(Nvar - 1)*sizeof(double));
	double *TempQy = malloc(Np*K*(Nvar - 1)*sizeof(double));
	char *chn = "N";
	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

	double hcrit = mxGetScalar(prhs[0]);
	signed char *ftype = (signed char *)mxGetData(prhs[1]);
	NdgMeshType type = (NdgMeshType)ftype[0];
	/*If the studied problem is two-dimensional*/
	if (type == Two){
		for (int field = 0; field < Nvar - 1; field++){  //the water depth field has to be excluded from this part
			//parellel computation to be added here
			for (int k = 0; k < K; k++){
				/*calculate variable first. In matlab, the index begins with one, and not zero. The index for H in matlab is one and diffusion is not included for H, so we need to escape this field */
				DotCriticalDivide(Var + field*K*Np + k*Np, fphys + (int)varIndex[field]*K*Np + k*Np, hcrit, fphys + k*Np, Np);
			}

			for (int k = 0; k < K; k++){
				/*$Qr = Dr*fphys$*/
				dgemm(chn, chn, &np, &oneI, &np, &one, Dr, &np, Var + field*K*Np + k*Np, &np, &zero, Qx + field*K*Np + k*Np, &np);
				/*$Qx1 = rx .* Qr$*/
				DotProduct( Qx + field*K*Np + k*Np, rx + k*Np, Np);
				/*$Qs = Ds*fphys$*/
				dgemm(chn, chn, &np, &oneI, &np, &one, Ds, &np, Var + field*K*Np + k*Np, &np, &zero, TempQx + field*K*Np + k*Np, &np);
				/*$Qx2 = sx .* Qs$*/
				DotProduct(TempQx + field*K*Np + k*Np, sx + k*Np, Np);
				/*$Qx = rx .* (Dr*fphys) + sx .* (Ds*fphys)$*/
				Add(Qx + field*K*Np + k*Np, TempQx + field*K*Np + k*Np, Np);
			}
		}

		for (int var = 0;var)
	}
	else if (type == Three){

	}
}