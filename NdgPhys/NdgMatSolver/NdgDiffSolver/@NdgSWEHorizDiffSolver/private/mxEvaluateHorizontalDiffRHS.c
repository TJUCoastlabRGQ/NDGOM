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

void GetFacialFluxTerm(double *dest, double *source, double *n, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = source[i] * n[i];
}
void GetCentralNumFlux(double *dest, double *fm, double *fp, double *n, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = n[i] * (fm[i] + fp[i]) / 2;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//Order of the input variable is 0 hcrit, 1 meshUnion.type, 2 prantl, 3 InnerEdge, 4 BoundaryEdge, 5 nv, 6 frhs, 7 fphys, 8 varIndex, 9 cell, 10 mesh, 11 BoundaryEdgefp
	//BoundaryEdgefm, InnerEdgefm and InnerEdgefp have been deleted, because we can fetch them in part.
	int Np, K, Nvar;
	const size_t *PRHS;
	double *Tempnv = mxGetPr(prhs[5]);
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
	double *J = mxGetField(prhs[10], 0, "J");
	double *varIndex = mxGetPr(prhs[8]);
	double *fphys = mxGetPr(prhs[7]);
	double *TempIENe = mxGetField(prhs[3],0,"Ne");
	int IENe = (int)TempIENe[0];
	double *TempIENfp = mxGetField(prhs[3], 0, "Nfp");
	int IENfp = (int)TempIENfp[0];
	double *IEMb = mxGetField(prhs[3], 0, "M");
	double *IEJs = mxGetField(prhs[3], 0, "Js");
	double *IEnx = mxGetField(prhs[3], 0, "nx");
	double *IEny = mxGetField(prhs[3], 0, "ny");
	double *TempBENe = mxGetField(prhs[4], 0, "Ne");
	int BENe = (int)TempBENe[0];
	double *TempBENfp = mxGetField(prhs[4], 0, "Nfp");
	int BENfp = (int)TempBENfp[0];
	double *BEMb = mxGetField(prhs[4], 0, "M");
	double *BEJs = mxGetField(prhs[4], 0, "Js");
	double *BEnx = mxGetField(prhs[4], 0, "nx");
	double *BEny = mxGetField(prhs[4], 0, "ny");

	int Nfield;
	double *varIndex = mxGetPr(prhs[8]);
	double *fphys = mxGetPr(prhs[7]);
	double *Hcrit = mxGetPr(prhs[0]);
	double *BEfp = mxGetPr(prhs[11]);
	signed char *Dimension = (signed char *)mxGetData(prhs[1]);
	NdgMeshType type = (NdgMeshType)Dimension[0];
	double *variable = NULL, *TempBEfp = NULL;
	double *nv = malloc(Np*K*sizeof(double));
	/*Allocate memory and calcualte variable u, v, $\theta$*/
	if (type == Two){
		/*For 2d shallow water problem, no horizontal diffusion terms are included in the governing equation for water depth $H$*/
		Nfield = Nvar - 1;
		variable = malloc(Np*K*(Nvar - 1)*sizeof(double));
		TempBEfp = malloc(BENfp*BENe*(Nvar - 1)*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K; k++){
			/*set the diffusion coefficient $H\nv$*/
			DotProduct(nv + k*Np, Tempnv + k*Np, fphys + (int)(varIndex[0] - 1)*Np*K + k*Np, Np);
			for (int field = 0; field < Nvar - 1; field++){
				DotCriticalDivide(variable + field*Np*K + k*Np, \
					fphys + (int)(varIndex[field + 1] - 1)*Np*K + k*Np, Hcrit, \
					fphys + (int)(varIndex[0] - 1)*Np*K + k*Np, Np);//For 2d shallow water problem, variable about height is organized as the first variable
			}
		}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int e = 0; e < BENe; e++){
			for (int field = 0; field < Nvar - 1; field++){
				DotCriticalDivide(TempBEfp + field*BENfp*BENe + e*BENfp, \
					BEfp + (int)(varIndex[field + 1] - 1)*BENfp*BENe + e*BENfp, Hcrit, \
					BEfp + (int)(varIndex[0] - 1)*BENfp*BENe + e*BENfp, BENfp);
			}
		}
	else if (type == Three){
		Nfield = Nvar;
		variable = malloc(Np*K*Nvar*sizeof(double));
		TempBEfp = malloc(BENfp*BENe*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K; k++){
			DotProduct(nv + k*Np, Tempnv + k*Np, fphys + 3 * Np*K + k*Np, Np);
			for (int field = 0; field < Nvar - 1; field++){
				DotCriticalDivide(variable + field*Np*K + k*Np, \
					fphys + (int)(varIndex[field] - 1)*Np*K + k*Np, Hcrit, \
					fphys + 3*Np*K + k*Np, Np);//For 3d shallow water problem, variable about height is organized as the forth variable
			}
		}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int e = 0; e < BENe; e++){
			for (int field = 0; field < Nvar; field++){
				DotCriticalDivide(TempBEfp + field*BENfp*BENe + e*BENfp, \
					BEfp + (int)(varIndex[field] - 1)*BENfp*BENe + e*BENfp, Hcrit, \
					BEfp + 3*BENfp*BENe + e*BENfp, BENfp);
			}
		}
	}
	}
	/*Allocate memory for the following computation*/
	/*Allocate memory for auxiallary variable $q_x=\frac{\partial u(v,\theta)}{\partial x}$*/
	double *AVx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for auxiallary variable $q_y=\frac{\partial u(v,\theta)}{\partial y}$*/
	double *AVy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the volumn integral of auxiallary variable $q_{x1}=\bold{r_x}\cdot (D_r*u(v,\theta))$*/
	double *VAVx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the volumn integral of auxiallary variable $q_{x2}=\bold{s_x}\cdot (D_s*u(v,\theta))$*/
	double *TempVAVx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the volumn integral of auxiallary variable $q_{y1}=\bold{r_y}\cdot (D_r*u(v,\theta))$*/
	double *VAVy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the volumn integral of auxiallary variable $q_{y2}=\bold{s_y}\cdot (D_s*u(v,\theta))$*/
	double *TempVAVy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the local face value at inner edge$*/
	double *IEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent face value at inner edge$*/
	double *IEfp = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local flux term at inner edge$*/
	double *IEFluxM = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent flux term at inner edge$*/
	double *IEFluxP = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at inner edge$*/
	double *IEFluxS = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value at boundary edge$*/
	double *BEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local flux term at boundary edge$*/
	double *BEFluxM = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at boundary edge$*/
	double *BEFluxS = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for interior edge contribution to AVx*/
	double *IERHSX = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for interior edge contribution to AVy*/
	double *IERHSY = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for boundary edge contribution to AVx*/
	double *BERHSX = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for boundary edge contribution to AVy*/
	double *BERHSY = malloc(Np*K*Nfield*sizeof(double));


	char *chn = "N";
	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;
	/*Volume integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			/*$Dr*u(v,\theta)$*/
			dgemm("N", "N", &np, &oneI, &np, &one, &Dr, &np, fphys + field*Np*K + k*Np, &np, &zero, VAVx + field*Np*K + k*Np, &np);
			/*$Ds*u(v,\theta)$*/
			dgemm("N", "N", &np, &oneI, &np, &one, &Dr, &np, fphys + field*Np*K + k*Np, &np, &zero, TempVAVx + field*Np*K + k*Np, &np);
			/*$rx\cdot Dr*u(v,\theta)$*/
			DotProduct(VAVx + field*Np*K + k*Np, VAVx + field*Np*K + k*Np, rx + k*Np, Np);
			/*$sx\cdot Ds*u(v,\theta)$*/
			DotProduct(TempVAVx + field*Np*K + k*Np, TempVAVx + field*Np*K + k*Np, sx + k*Np, Np);
			/*rx\cdot Dr*u(v,\theta) + sx\cdot Ds*u(v,\theta)*/
			Add(VAVx + field*Np*K + k*Np, VAVx + field*Np*K + k*Np, TempVAVx + field*Np*K + k*Np, Np);
			/*$Dr*u(v,\theta)$*/
			dgemm("N", "N", &np, &oneI, &np, &one, &Dr, &np, fphys + field*Np*K + k*Np, &np, &zero, VAVy + field*Np*K + k*Np, &np);
			/*$Ds*u(v,\theta)$*/
			dgemm("N", "N", &np, &oneI, &np, &one, &Dr, &np, fphys + field*Np*K + k*Np, &np, &zero, TempVAVy + field*Np*K + k*Np, &np);
			/*$ry\cdot Dr*u(v,\theta)$*/
			DotProduct(VAVy + field*Np*K + k*Np, VAVy + field*Np*K + k*Np, ry + k*Np, Np);
			/*$sy\cdot Ds*u(v,\theta)$*/
			DotProduct(TempVAVy + field*Np*K + k*Np, TempVAVy + field*Np*K + k*Np, sy + k*Np, Np);
			/*ry\cdot Dr*u(v,\theta) + sy\cdot Ds*u(v,\theta)*/
			Add(VAVy + field*Np*K + k*Np, VAVy + field*Np*K + k*Np, TempVAVy + field*Np*K + k*Np, Np);
		}
	}
/*Allocate memory for local primitive diffusion term in both x and y direction, $\nv\Nabla u(v,\theta)$*/
/*This part is used because we need to extract the trace of the local derivative operator to compute the numerical flux*/
	double *LocalPrimitiveDiffTermX = malloc(Np*K*Nfield);
	double *LocalPrimitiveDiffTermY = malloc(Np*K*Nfield);
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			DotProduct(LocalPrimitiveDiffTermX + field*Np*K + k*Np, VAVx + field*Np*K + k*Np, nv + k*Np, Np);
			DotProduct(LocalPrimitiveDiffTermY + field*Np*K + k*Np, VAVy + field*Np*K + k*Np, nv + k*Np, Np);
		}
		/*for substance transport, prantl number is considered*/
		for (int field = 2; field < Nfield; field++){
			DotDivideByConstant(LocalPrimitiveDiffTermX + field*Np*K + k*Np, VAVx + field*Np*K + k*Np, Prantl, Np);
			DotDivideByConstant(LocalPrimitiveDiffTermY + field*Np*K + k*Np, VAVy + field*Np*K + k*Np, Prantl, Np);
		}
	}

/*Facial integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		for (int field = 0; field < Nfield; field++){
			FetchInnerEdgeFacialValue(IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, Np, IENfp);
			GetFacialFluxTerm(IEFluxM + field*IENe*IENfp + face*IENfp, IEfm + field*IENe*IENfp + face*IENfp, IEnx + face*IENfp,IENfp);
			GetFacialFluxTerm(IEFluxP + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, IEnx + face*IENfp,IENfp);
			GetCentralNumFlux(IEFluxS + field*IENe*IENfp, IEfm + field*IENe*IENfp, IEfp + field*IENe*IENfp + face*IENfp, IEnx + face*IENfp, IENfp);
		}
	}
	
}