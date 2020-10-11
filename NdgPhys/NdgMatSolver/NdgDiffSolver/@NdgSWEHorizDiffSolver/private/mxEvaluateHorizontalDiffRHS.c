#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"
#include "blas.h"
#include "matrix.h"
#include "..\..\..\..\..\NdgMath\NdgMath.h"
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
void GetIPNumFlux(double *dest, double *fm, double *fp, double *n, double *vfm,\
	double *vfp, double *Jumpn, int Nfp, double *Tau, double Coefficient ){
	for (int i = 0; i < Nfp; i++)
		dest[i] = n[i] * (fm[i] + fp[i]) / 2 - n[i] * Tau[i] / Coefficient*(Jumpn[i] * (vfm[i] - vfp[i]));
}

void GetIPBoundaryNumFlux(double *dest, double *fm, double *n, double *vfm, \
	double *vfp, double *Jumpn, int Nfp, double *Tau, double Coefficient){
	for (int i = 0; i < Nfp; i++)
		dest[i] = n[i] * fm[i] - n[i] * Tau[i] / Coefficient*(Jumpn[i] * (vfm[i] - vfp[i]));
}

void GetVolumnIntegral(double *dest, double *tempdest, ptrdiff_t *RowOPA, ptrdiff_t *ColOPB, ptrdiff_t *ColOPA, double *alpha,\
	double *Dr, double *Ds, ptrdiff_t *LDA, double *B, ptrdiff_t *LDB, double *Beta, ptrdiff_t *LDC, double *rx, double *sx){
	/*$Dr*u(v,\theta)$*/
	dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dr, LDA, B, LDB, Beta, dest, LDC);
	/*$Ds*u(v,\theta)$*/
	dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Ds, LDA, B, LDB, Beta, tempdest, LDC);
	/*$rx\cdot Dr*u(v,\theta)$*/
	DotProduct(dest, dest, rx, (int)(*LDC));
	/*$sx\cdot Ds*u(v,\theta)$*/
	DotProduct(tempdest, tempdest, sx, (int)(*LDC));
	/*rx\cdot Dr*u(v,\theta) + sx\cdot Ds*u(v,\theta)*/
	Add(dest, dest, tempdest, (int)(*LDC));
}
/*
* Purpose: This function is used to calculate the inner edge contribution to the auxialary variable
* Input:
* 		double[Np x 1]  dest the position to put the contribution due to the face integral on the studied inner face
*       int[ 1 ]  face index of the studied inner face
* 		int[ 1 ]  Nfp number of the interpolation points on the studied inner face
*       int[ 1 ]  field index of studied physical field, i.e. index corresponds to $u,v$ and $\theta$
*       int[Nfp x 1] fm the local physical value on the studied inner edge
*       int[Nfp x 1] fp the adjacent physical value on the studied inner edge
*       double[2 x 1] IEFToE the topology of the inner edge, it holds information about how the element adjacent to the inner edge
*       double[Nfp x 1] IEFToN1 index of the local facial interpolation points
*       double[Nfp x 1] IEFToN2 index of the adjacent facial interpolation points
* 		int[ 1 ]  Np number of the interpolation points of the computation cell
*       double[Nfp x 1] FluxM the local flux term, in this function, only the adress is required, this term is calculated in this function
*       double[Nfp x 1] FluxP the adjacent flux term, in this function, only the adress is required, this term is calculated in this function
*       double[Nfp x 1] FluxS the numerical flux term, in this funciton, only the address is required, this term is calculated in this function
*       double[Nfp x 1] n the direction vector of the studied face
*       double[Nfp x Nfp] Mb the mass matrix of the standard cell
*       double[Nfp x 1]  Js the faical Jacobian of the studied face
* Output:
* 		double[Np x 1] 	dest the contribution due to the face integral on the studied inner face
*/
void GetIEContributionToAuxialaryVariable(double *dest, int face, int Nfp, int field, double *fm, double *fp, double *IEFToE,\
	double *IEFToN1, double *IEFToN2, int Np, double *FluxM, double *FluxP, double *FluxS, double *n, double *Mb, double *Js){
	/*This */
	GetFacialFluxTerm(FluxM, fm, n, Nfp);
	GetFacialFluxTerm(FluxP, fp, n, Nfp);
	GetCentralNumFlux(FluxS, fm, fp, n, Nfp);
	StrongFormInnerEdgeRHS(face, IEFToE, Np, Nfp, IEFToN1, IEFToN2, FluxM, FluxP, FluxS, Js, Mb, dest);
}
/*
* Purpose: This function is used to calculate the boundary edge contribution to the auxialary variable
* Input:
* 		double[Np x 1]  dest the position to put the contribution due to the face integral on the studied inner face
*       int[ 1 ]  face index of the studied boundary  face
* 		int[ 1 ]  Nfp number of the interpolation points on the studied inner face
*       int[ 1 ]  field index of studied physical field, i.e. index corresponds to $u,v$ and $\theta$
*       int[Nfp x 1] fm the local physical value on the studied inner edge
*       int[Nfp x 1] fp the boundary condition on the studied boundary edge
*       double[2 x 1] BEFToE the topology of the inner edge, it holds information about how the element adjacent to the inner edge
*       double[Nfp x 1] BEFToN1 index of the local facial interpolation points
* 		int[ 1 ]  Np number of the interpolation points of the computation cell
*       double[Nfp x 1] FluxM the local flux term, in this function, only the adress is required, this term is calculated in this function
*       double[Nfp x 1] FluxS the numerical flux term, in this funciton, only the address is required, this term is calculated in this function
*       double[Nfp x 1] n the direction vector of the studied face
*       double[Nfp x Nfp] Mb the mass matrix of the standard cell
*       double[Nfp x 1]  Js the faical Jacobian of the studied face
* Output:
* 		double[Np x 1] 	dest the contribution due to the face integral on the studied inner face
*/
void GetBEContributionToAuxialaryVariable(double *dest, int face, int Nfp, int field, double *fm, double *fp, double *BEFToE, \
	double *BEFToN1, int Np, double *FluxM, double *FluxS, double *n, double *Mb, double *Js){
	GetFacialFluxTerm(FluxM, fm, n, Nfp);
	GetFacialFluxTerm(FluxS, fp, n, Nfp);
	StrongFormBoundaryEdgeRHS(face, BEFToE, Np, Nfp, BEFToN1, FluxM, FluxS, Js, Mb, dest);
}

/*
* Purpose: This function is used to calculate the inner edge contribution to the right hand side
* Input:
* 		double[Np x 1]  dest the position to put the contribution due to the face integral on the studied inner face
*       double[Nfp x 1]  LPDTfm local value of the local partial diffrential term on the studied inner face
*       double[Nfp x 1]  LPDTfp adjacent value of the local partial diffrential term on the studied inner face
* 		int[ 1 ]  face index of the studied inner face
*       double[Np x K]  LPDiffTerm the local partial differential term
*       double[2 x IENe] FToE the topology of the inner edge, it holds information about how the element adjacent to the inner edge
*       double[Nfp x 1] FToN1 index of the local facial interpolation points
*       double[Nfp x 1] FToN2 index of the adjacent facial interpolation points
* 		int[ 1 ]  Np number of the interpolation points of the computation cell
* 		int[ 1 ]  Nfp number of the interpolation points on the studied inner face
*       double[Nfp x 1] FluxS the numerical flux term, in this funciton, only the address is required, this term is calculated in this function
*       double[Nfp x 1] n direction vector used to multiply the numerical flux
*       double[Nfp x 1] fm the local physical value on the studied inner edge, used to calculate the jump term
*       double[Nfp x 1] fp the boundary condition on the studied boundary edge, used to calculate the jump term
*       double[Nfp x 1] Jumpn direction vector used to calculate the jump term
*       double[Nfp x 1] Tau the penalty parameter used in interior penalty flux
*       double[ 1 ] Coefficient the parameter stands for prantl number or 1.0
*       double[Nfp x 1] AVfm the local value of the auxialary variable
*       double[Nfp x 1] AVfp the adjacent value of the auxialary variable
*       double[Np x K] AV the auxialary variable
*       double[Nfp x 1] FluxM the local flux term relates to AVfm
*       double[Nfp x 1] FluxP the adjacent flux term relates to AVfp
*       double[Nfp x 1] Js the facial Jacobian parameter
*       double[Nfp x Nfp] M the mass matrix of the master cell corresponds to the studied face 
* Output:
* 		double[Np x 1] 	dest the contribution due to the face integral on the studied inner face
*/
void GetIEContributionToRHS(double *dest, double *LPDTfm, double *LPDTfp, int face, double *LPDiffTerm, double *FToE, double *FToN1, \
	double *FToN2, int Np, int Nfp, double *FluxS, double *n, double *fm, double *fp, double *Jumpn, double *Tau, double Coefficient, \
	double *AVfm, double *AVfp, double *AV, double *FluxM, double *FluxP, double *Js, double *M){
	/*FToE is the start stress of the property FToE contained in inner edge*/
	FetchInnerEdgeFacialValue(LPDTfm, LPDTfp, LPDiffTerm, FToE+2*face, FToN1, FToN2, Np, Nfp);
	GetIPNumFlux(FluxS, LPDTfm, LPDTfp, n, fm, fp, Jumpn, Nfp, Tau, Coefficient);
	FetchInnerEdgeFacialValue(AVfm, AVfp, AV, FToE+2*face, FToN1, FToN2, Np, Nfp);
	GetFacialFluxTerm(FluxM, AVfm, n, Nfp);
	GetFacialFluxTerm(FluxP, AVfp, n, Nfp);
	StrongFormInnerEdgeRHS(face, FToE, Np, Nfp, FToN1, FToN2, FluxM, FluxP, FluxS, Js, M, dest);
}

/*
* Purpose: This function is used to calculate the boundary edge contribution to the right hand side
* Input:
* 		double[Np x 1]  dest the position to put the contribution due to the face integral on the studied inner face
*       double[Nfp x 1]  LPDTfm local value of the local partial diffrential term on the studied inner face
* 		int[ 1 ]  face index of the studied inner face
*       double[Np x K]  LPDiffTerm the local partial differential term
*       double[2 x IENe] FToE the topology of the inner edge, it holds information about how the element adjacent to the inner edge
*       double[Nfp x 1] FToN1 index of the local facial interpolation points
* 		int[ 1 ]  Np number of the interpolation points of the computation cell
* 		int[ 1 ]  Nfp number of the interpolation points on the studied inner face
*       double[Nfp x 1] FluxS the numerical flux term, in this funciton, only the address is required, this term is calculated in this function
*       double[Nfp x 1] n direction vector used to multiply the numerical flux
*       double[Nfp x 1] fm the local physical value on the studied inner edge, used to calculate the jump term
*       double[Nfp x 1] fp the boundary condition on the studied boundary edge, used to calculate the jump term
*       double[Nfp x 1] Jumpn direction vector used to calculate the jump term
*       double[Nfp x 1] Tau the penalty parameter used in interior penalty flux
*       double[ 1 ] Coefficient the parameter stands for prantl number or 1.0
*       double[Nfp x 1] AVfm the local value of the auxialary variable
*       double[Np x K] AV the auxialary variable
*       double[Nfp x 1] FluxM the local flux term relates to AVfm
*       double[Nfp x 1] Js the facial Jacobian parameter
*       double[Nfp x Nfp] M the mass matrix of the master cell corresponds to the studied face
* Output:
* 		double[Np x 1] 	dest the contribution due to the face integral on the studied inner face
*/

void GetBEContributionToRHS(double *dest, double *LPDTfm, int face, double *LPDiffTerm, double *FToE, double *FToN1, \
    int Np, int Nfp, double *FluxS, double *n, double *fm, double *fp, double *Jumpn, double *Tau, double Coefficient, \
	double *AVfm, double *AV, double *FluxM, double *Js, double *M){
	/*FToE is the start adress of the property FToE contained in inner edge*/
	FetchBoundaryEdgeFacialValue(LPDTfm, LPDiffTerm, FToE + 2 * face, FToN1, Np, Nfp);
	GetIPBoundaryNumFlux(FluxS, LPDTfm, n, fm, fp, Jumpn, Nfp, Tau, Coefficient);
	FetchBoundaryEdgeFacialValue(AVfm, AV, FToE + 2 * face, FToN1, Np, Nfp);
	GetFacialFluxTerm(FluxM, AVfm, n, Nfp);
	StrongFormBoundaryEdgeRHS(face, FToE, Np, Nfp, FToN1, FluxM, FluxS, Js, M, dest);	
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//Order of the input variable is 0 hcrit, 1 meshUnion.type, 2 prantl, 3 InnerEdge, 4 BoundaryEdge, 5 nv, 6 frhs, 7 fphys, 8 varIndex, 9 cell, 10 mesh, 11 BoundaryEdgefp
	double *Hcrit = mxGetPr(prhs[0]);
	signed char *Dimension = (signed char *)mxGetData(prhs[1]);
	NdgMeshType type = (NdgMeshType)Dimension[0];
	double Prantl = mxGetScalar(prhs[2]);
	mxArray *TempIENe = mxGetField(prhs[3], 0, "Ne");
	int IENe = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp = mxGetField(prhs[3], 0, "Nfp");
	int IENfp = (int)mxGetScalar(TempIENfp);
	mxArray *TempIEMb = mxGetField(prhs[3], 0, "M");
	double *IEMb = mxGetPr(TempIEMb);
	mxArray *TempIEJs = mxGetField(prhs[3], 0, "Js");
	double *IEJs = mxGetPr(TempIEJs);
	mxArray *TempIEnx = mxGetField(prhs[3], 0, "nx");
	double *IEnx = mxGetPr(TempIEnx);
	mxArray *TempIEny = mxGetField(prhs[3], 0, "ny");
	double *IEny = mxGetPr(TempIEny);
	mxArray *TempIELAV = mxGetField(prhs[3], 0, "LAV");
	double *IELAV = mxGetPr(TempIELAV);
	mxArray *TempIEFToE = mxGetField(prhs[3], 0, "FToE");
	double *IEFToE = mxGetPr(TempIEFToE);
	mxArray *TempIEFToN1 = mxGetField(prhs[3], 0, "FToN1");
	double *IEFToN1 = mxGetPr(TempIEFToN1);
	mxArray *TempIEFToN2 = mxGetField(prhs[3], 0, "FToN2");
	double *IEFToN2 = mxGetPr(TempIEFToN2);

	mxArray *TempBENe = mxGetField(prhs[4], 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBENfp = mxGetField(prhs[4], 0, "Nfp");
	int BENfp = mxGetScalar(TempBENfp);
	mxArray *TempBEMb = mxGetField(prhs[4], 0, "M");
	double *BEMb = mxGetPr(TempBEMb);
	mxArray *TempBEJs = mxGetField(prhs[4], 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	mxArray *TempBEnx = mxGetField(prhs[4], 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(prhs[4], 0, "ny");
	double *BEny = mxGetPr(TempBEny);
	mxArray *TempBELAV = mxGetField(prhs[4], 0, "LAV");
	double *BELAV = mxGetPr(TempBELAV);
	mxArray *TempBEFToE = mxGetField(prhs[4], 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToN1 = mxGetField(prhs[4], 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);



	double *Tempnv = mxGetPr(prhs[5]);
	int Np, K, Nvar;
	const size_t *PRHS;
	PRHS = mxGetDimensions(prhs[6]);
	Np = (int)PRHS[0];
	K = (int)PRHS[1];
	Nvar = (int)PRHS[2];
	double *InputRHS = mxGetPr(prhs[6]);
	double *fphys = mxGetPr(prhs[7]);
	double *varIndex = mxGetPr(prhs[8]);
	mxArray *TempDr = mxGetField(prhs[9], 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(prhs[9], 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempNface = mxGetField(prhs[9], 0, "Nface");
	int Nface = mxGetScalar(TempNface);
	mxArray *Temprx = mxGetField(prhs[10], 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(prhs[10], 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(prhs[10], 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *Tempsy = mxGetField(prhs[10], 0, "sy");
	double *sy = mxGetPr(Tempsy);
	mxArray *TempJ = mxGetField(prhs[10], 0, "J");
	double *J = mxGetPr(TempJ);
	mxArray *TempLAV = mxGetField(prhs[10], 0, "LAV");
	double *LAV = mxGetPr(TempLAV);
	double *TempBEfp = mxGetPr(prhs[11]);
	/*Set the output right hand side*/
	const size_t NdimOut = 3;
	const mwSize dimOut[3] = { Np, K, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputRHS = mxGetPr(plhs[0]);
	memcpy(OutputRHS, InputRHS, Np*K*Nvar*sizeof(double));

	int Nfield;
	double *variable = NULL, *BEfp = NULL;
	mxArray *TempOrder;
	double *nv = malloc(Np*K*sizeof(double));
	int MeshType, Order;
	/*Allocate memory and calcualte variable u, v, $\theta$*/
	if (type == Two){
		/*For 2d shallow water problem, no horizontal diffusion terms are included in the governing equation for water depth $H$*/
		Nfield = Nvar - 1;
		variable = malloc(Np*K*(Nvar - 1)*sizeof(double));
		BEfp = malloc(BENfp*BENe*(Nvar - 1)*sizeof(double));
		MeshType = 2;
		TempOrder = mxGetField(prhs[9], 0, "N");
		Order = mxGetScalar(TempOrder);
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
				DotCriticalDivide(BEfp + field*BENfp*BENe + e*BENfp, \
					TempBEfp + (int)(varIndex[field + 1] - 1)*BENfp*BENe + e*BENfp, Hcrit, \
					TempBEfp + (int)(varIndex[0] - 1)*BENfp*BENe + e*BENfp, BENfp);
			}
		}
	}
	else if (type == Three){
		Nfield = Nvar;
		variable = malloc(Np*K*Nvar*sizeof(double));
		BEfp = malloc(BENfp*BENe*Nvar*sizeof(double));
		MeshType = 3;
		TempOrder = mxGetField(prhs[9], 0, "N");
		Order = (int)mxGetScalar(TempOrder);
		TempOrder = mxGetField(prhs[9], 0, "Nz");
		Order = max(Order, (int)mxGetScalar(TempOrder));
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
	/*Allocate memory for the following computation*/
	/*Allocate memory for auxiallary variable $q_x=\frac{\partial u(v,\theta)}{\partial x}$*/
	double *AVx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for auxiallary variable $q_y=\frac{\partial u(v,\theta)}{\partial y}$*/
	double *AVy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for part of the volumn integral of auxiallary variable $q_{x1}=\bold{r_x}\cdot (D_r*u(v,\theta))$,
	and for part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	double *Vx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the rest part of the volumn integral of auxiallary variable $q_{x2}=\bold{r_x}\cdot (D_r*u(v,\theta))$,
	and for the rest part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	double *TempVx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for part of the volumn integral of auxiallary variable $q_{y1}=\bold{r_y}\cdot (D_r*u(v,\theta))$,
	and for part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	double *Vy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the rest part of the volumn integral of auxiallary variable $q_{y2}=\bold{s_y}\cdot (D_s*u(v,\theta))$,
	and for the rest part of the volumn integral part of the second order operator $\frac{\partial q}{\partial y}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$*/
	double *TempVy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the local face value at inner edge*/
	double *IEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at inner edge*/
	double *AVIEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent face value at inner edge*/
	double *IEfp = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent face value of the auxialary variable at inner edge*/
	double *AVIEfp = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local flux term at inner edge*/
	double *IEFluxM = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent flux term at inner edge*/
	double *IEFluxP = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at inner edge*/
	double *IEFluxS = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value at boundary edge*/
	double *BEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at boundary edge*/
	double *AVBEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local flux term at boundary edge*/
	double *BEFluxM = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at boundary edge*/
	double *BEFluxS = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for interior edge contribution to right hand side in x direction*/
	double *IERHSX = malloc(Np*K*Nfield*sizeof(double));
	memset(IERHSX, 0, Np*K*Nfield*sizeof(double));
	/*Allocate memory for interior edge contribution to right hand side in y direction*/
	double *IERHSY = malloc(Np*K*Nfield*sizeof(double));
	memset(IERHSY, 0, Np*K*Nfield*sizeof(double));
	/*Allocate memory for boundary edge contribution to right hand side in x direction*/
	double *BERHSX = malloc(Np*K*Nfield*sizeof(double));
	memset(BERHSX, 0, Np*K*Nfield*sizeof(double));
	/*Allocate memory for boundary edge contribution to right hand side in y direction*/
	double *BERHSY = malloc(Np*K*Nfield*sizeof(double));
	memset(BERHSY, 0, Np*K*Nfield*sizeof(double));

	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;
	/*Volume integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			/*$\bold{r_x}\cdot (Dr*u(v,\theta))+\bold{s_x}\cdot (Ds*u(v,\theta))$*/			
			GetVolumnIntegral(Vx + field*Np*K + k*Np, TempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, variable + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
			/*$\bold{r_y}\cdot (Dr*u(v,\theta))+\bold{s_y}\cdot (Ds*u(v,\theta))$*/
			GetVolumnIntegral(Vy + field*Np*K + k*Np, TempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, variable + field*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
 		}
	}
/*Allocate memory for local primitive diffusion term in both x and y direction, $\nv\Nabla u(v,\theta)$*/
/*This part is used because we need to extract the trace of the local derivative operator to compute the numerical flux*/
	double *LocalPrimitiveDiffTermX = malloc(Np*K*Nfield*sizeof(double));
	double *LocalPrimitiveDiffTermY = malloc(Np*K*Nfield*sizeof(double));
	double *LPDTIEfm = malloc(IENfp*IENe*sizeof(double));
	double *LPDTIEfp = malloc(IENfp*IENe*sizeof(double));
	double *LPDTBEfm = malloc(BENfp*BENe*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			DotProduct(LocalPrimitiveDiffTermX + field*Np*K + k*Np, Vx + field*Np*K + k*Np, nv + k*Np, Np);
			DotProduct(LocalPrimitiveDiffTermY + field*Np*K + k*Np, Vy + field*Np*K + k*Np, nv + k*Np, Np);
		}
		/*for substance transport, prantl number is considered*/
		for (int field = 2; field < Nfield; field++){
			DotDivideByConstant(LocalPrimitiveDiffTermX + field*Np*K + k*Np, Vx + field*Np*K + k*Np, Prantl, Np);
			DotDivideByConstant(LocalPrimitiveDiffTermY + field*Np*K + k*Np, Vy + field*Np*K + k*Np, Prantl, Np);
		}
	}

/*Inner edge facial integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int field = 0; field < Nfield; field++){
		for (int face = 0; face < IENe; face++){
			FetchInnerEdgeFacialValue(IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, variable + field*Np*K, IEFToE + 2 * face, \
				IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
			/*Inner edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$*/
			GetIEContributionToAuxialaryVariable(IERHSX + field*Np*K, face, IENfp, field, IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, \
				IEFToE, IEFToN1, IEFToN2, Np, IEFluxM + field*IENe*IENfp + face*IENfp, IEFluxP + field*IENe*IENfp + face*IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
				IEnx + face*IENfp, IEMb, IEJs + face*IENfp);
			/*Inner edge contribution to RHSY of $\frac{\partial u(v,\theta)}{\partial y}$*/
			GetIEContributionToAuxialaryVariable(IERHSY + field*Np*K, face, IENfp, field, IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, \
				IEFToE, IEFToN1, IEFToN2, Np, IEFluxM + field*IENe*IENfp + face*IENfp, IEFluxP + field*IENe*IENfp + face*IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
				IEny + face*IENfp, IEMb, IEJs + face*IENfp);
		}
	}
	/*Boundary edge facial integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int field = 0; field < Nfield; field++){
		for (int face = 0; face < BENe; face++){
			FetchBoundaryEdgeFacialValue(BEfm + field*BENe*BENfp + face*BENfp, variable + field*Np*K, BEFToE + 2 * face, \
				BEFToN1 + BENfp*face, Np, BENfp);
			/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$, we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToAuxialaryVariable(BERHSX + field*Np*K, face, BENfp, field, BEfm + field*BENe*BENfp + face*BENfp, \
				BEfp + field*BENe*BENfp + face*BENfp, BEFToE, BEFToN1, Np, BEFluxM + field*BENe*BENfp + face*BENfp,\
				BEFluxS + field*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BEMb, BEJs + face*BENfp);

			/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial y}$, we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToAuxialaryVariable(BERHSY + field*Np*K, face, BENfp, field, BEfm + field*BENe*BENfp + face*BENfp, \
				BEfp + field*BENe*BENfp + face*BENfp, BEFToE, BEFToN1, Np, BEFluxM + field*BENe*BENfp + face*BENfp,\
				BEFluxS + field*BENe*BENfp + face*BENfp, BEny + face*BENfp, BEMb, BEJs + face*BENfp);
		}
	}
	/*Next, sum contribution from volume integral, inner edge contribution and boundary edge contribution into AVx and Avy*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			Minus(AVx + field*Np*K + k*Np, \
				Vx + field*Np*K + k*Np, IERHSX + field*Np*K + k*Np, Np);
			Minus(AVx + field*Np*K + k*Np, \
				AVx + field*Np*K + k*Np, BERHSX + field*Np*K + k*Np, Np);
			Minus(AVy + field*Np*K + k*Np, \
				Vy + field*Np*K + k*Np, IERHSY + field*Np*K + k*Np, Np);
			Minus(AVy + field*Np*K + k*Np, \
				AVy + field*Np*K + k*Np, BERHSY + field*Np*K + k*Np, Np);
		}
	}
	/*Finally, multiply each component of AVx and AVy by its diffusion coefficient to get the final auxialary variable,
	this is conducted according to $(Q,v)=(\nv\hat Q,v)$,where $\hat Q = \frac{\partial u(v,\theta)}{\partial (x,y)}$,
	we note that in this projection precedure, aliasing error is introduced.*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			DotProduct(AVx + field*Np*K + k*Np, AVx + field*Np*K + k*Np, nv + k*Np, Np);
			DotProduct(AVy + field*Np*K + k*Np, AVy + field*Np*K + k*Np, nv + k*Np, Np);
		}
		for (int field = 2; field < Nfield; field++){
			DotDivideByConstant(AVx + field*Np*K + k*Np, AVx + field*Np*K + k*Np, Prantl, Np);
			DotDivideByConstant(AVy + field*Np*K + k*Np, AVy + field*Np*K + k*Np, Prantl, Np);
		}
	}
	double *InnerEdgeTau = malloc(IENe*IENfp*sizeof(double));
	double *BoundaryEdgeTau = malloc(BENe*BENfp*sizeof(double));
	double *IEnvfm = malloc(IENe*IENfp*sizeof(double));
	double *IEnvfp = malloc(IENe*IENfp*sizeof(double));
	double *BEnvfm = malloc(BENe*BENfp*sizeof(double));
	/*Calculate the contribution to the right hand side due to the auxialary variable AVx and AVy with IPDG.*/
	/*Calculate the penalty parameter $\tau$ first, this parameter is calculated as $\tau=\frac{(N+1)(N+d)}{d}\frac{n_0}{2}\frac{A}{V}\nv$*/
	/*Inner edge first*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		FetchInnerEdgeFacialValue(IEnvfm + face*IENfp, IEnvfp + face*IENfp, nv, \
			IEFToE + face * 2, IEFToN1 + face*IENfp, IEFToN2 + face*IENfp, Np, IENfp);
		double localRatio = IELAV[face] / LAV[(int)IEFToE[face * 2]];
		double adjacentRatio = IELAV[face] / LAV[(int)IEFToE[face * 2 + 1]];
		for (int p = 0; p < IENfp; p++){
			InnerEdgeTau[face*IENfp + p] = max((Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * IEnvfm[face*IENfp + p], \
				(Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * IEnvfm[face*IENfp + p]);
		}
	}
	/*Boundary edge next*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		FetchBoundaryEdgeFacialValue(BEnvfm + face*IENfp, nv, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
		double localRatio = BELAV[face] / LAV[(int)BEFToE[face * 2]];
		for (int p = 0; p < BENfp; p++){
			BoundaryEdgeTau[face*BENfp + p] = (Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * BEnvfm[face*IENfp + p];
		}
	}
	/*Volume integral of the second order operator*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		/*Prantl number is not considered here, since we have considered this parameter when calculating the auxialary variable,
		it is not need anymore in the volumn integral part when considering the contribution of second order operator to the RHS*/
		for (int field = 0; field < Nfield; field++){
			/*$\bold{r_x}\cdot (Dr*Q_x)+\bold{s_x}\cdot (Ds*Q_x)$*/	
			GetVolumnIntegral(Vx + field*Np*K + k*Np, TempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, AVx + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
			/*$\bold{r_y}\cdot (Dr*Q_y)+\bold{s_y}\cdot (Ds*Q_y)$*/
			GetVolumnIntegral(Vy + field*Np*K + k*Np, TempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, AVy + field*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
			if (type == Two){
				/*The water depth field is excluded from this part*/
				Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					Vx + field*Np*K + k*Np, Np);
				Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					Vy + field*Np*K + k*Np, Np);
			}
			else if (type == Three){
				Add(OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					Vx + field*Np*K + k*Np, Np);
				Add(OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					Vy + field*Np*K + k*Np, Np);
			}
		}
		/*If Mellor and Blumberg's method in "Modeling vertical and horizontal diffusivities with the sigma coordinate system" is
		adopted, the next line to the further 44 lines should be ignored( current from 393 to 437)*/
		int field = 0;
		/*$\bold{r_x}\cdot (Dr*Q_x)+\bold{s_x}\cdot (Ds*Q_x)$, with $*Q_x=\nv H\frac{\partial u}{\partial x}*$*/
		GetVolumnIntegral(Vx + field*Np*K + k*Np, TempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, AVx + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
		/*$\bold{r_y}\cdot (Dr*Q_x)+\bold{s_y}\cdot (Ds*Q_x)$, with $*Q_x=\nv H\frac{\partial v}{\partial x}*$*/
		GetVolumnIntegral(Vy + field*Np*K + k*Np, TempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, AVx + (field+1)*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		field = 1;
		/*$\bold{r_x}\cdot (Dr*Q_y)+\bold{s_x}\cdot (Ds*Q_y)$, with $*Q_y=\nv H\frac{\partial u}{\partial y}*$*/
		GetVolumnIntegral(Vx + field*Np*K + k*Np, TempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, AVy + (field-1)*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
		/*$\bold{r_y}\cdot (Dr*Q_y)+\bold{s_y}\cdot (Ds*Q_y)$, with $*Q_y=\nv H\frac{\partial v}{\partial y}*$*/
		GetVolumnIntegral(Vy + field*Np*K + k*Np, TempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, AVy + field*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		for (int field = 0; field < 2; field++){
			if (type == Two){
				/*The water depth field is excluded from this part*/
				Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					Vx + field*Np*K + k*Np, Np);
				Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					Vy + field*Np*K + k*Np, Np);
			}
			else if (type == Three){
				Add(OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					Vx + field*Np*K + k*Np, Np);
				Add(OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					Vy + field*Np*K + k*Np, Np);
			}
		}
	}
	/*Reset all the data contained in IERHSX and IERHSY to zero, becaused these space contains the data left when calculating the auxialary variable*/
	memset(IERHSX, 0, Np*K*Nfield*sizeof(double));
	memset(IERHSY, 0, Np*K*Nfield*sizeof(double));
	/*Surface integral of the second order operator*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int field = 0; field < 2; field++){
		for (int face = 0; face < IENe; face++){
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
			GetIEContributionToRHS(IERHSX + field*Np*K, LPDTIEfm + field*IENe*IENfp + face*IENfp, LPDTIEfp + field*IENe*IENfp + face*IENfp, \
				face, LocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
				IEnx + face*IENfp, IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, IEnx + face*IENfp, InnerEdgeTau + face*IENfp, 1.0, \
				AVIEfm + field*IENe*IENfp + face*IENfp, AVIEfp + field*IENe*IENfp + face*IENfp, AVx + field*Np*K, IEFluxM + field*IENe*IENfp + face*IENfp, \
				IEFluxP + field*IENe*IENfp + face*IENfp, IEJs + face*IENfp, IEMb);
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial y}$*/
			GetIEContributionToRHS(IERHSY + field*Np*K, LPDTIEfm + field*IENe*IENfp + face*IENfp, LPDTIEfp + field*IENe*IENfp + face*IENfp, \
				face, LocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
				IEny + face*IENfp, IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, IEny + face*IENfp, InnerEdgeTau + face*IENfp, 1.0, \
				AVIEfm + field*IENe*IENfp + face*IENfp, AVIEfp + field*IENe*IENfp + face*IENfp, AVy + field*Np*K, IEFluxM + field*IENe*IENfp + face*IENfp, \
				IEFluxP + field*IENe*IENfp + face*IENfp, IEJs + face*IENfp, IEMb);

		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int field = 2; field < Nfield; field++){
		for (int face = 0; face < IENe; face++){
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$*/
			GetIEContributionToRHS(IERHSX + field*Np*K, LPDTIEfm + field*IENe*IENfp + face*IENfp, LPDTIEfp + field*IENe*IENfp + face*IENfp, \
				face, LocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
				IEnx + face*IENfp, IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, IEnx + face*IENfp, InnerEdgeTau + face*IENfp, Prantl, \
				AVIEfm + field*IENe*IENfp + face*IENfp, AVIEfp + field*IENe*IENfp + face*IENfp, AVx + field*Np*K, IEFluxM + field*IENe*IENfp + face*IENfp, \
				IEFluxP + field*IENe*IENfp + face*IENfp, IEJs + face*IENfp, IEMb);
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$*/
			GetIEContributionToRHS(IERHSY + field*Np*K, LPDTIEfm + field*IENe*IENfp + face*IENfp, LPDTIEfp + field*IENe*IENfp + face*IENfp, \
				face, LocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
				IEny + face*IENfp, IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, IEny + face*IENfp, InnerEdgeTau + face*IENfp, Prantl, \
				AVIEfm + field*IENe*IENfp + face*IENfp, AVIEfp + field*IENe*IENfp + face*IENfp, AVy + field*Np*K, IEFluxM + field*IENe*IENfp + face*IENfp, \
				IEFluxP + field*IENe*IENfp + face*IENfp, IEJs + face*IENfp, IEMb);

		}
	}
	memset(BERHSX, 0, Np*K*Nfield*sizeof(double));
	memset(BERHSY, 0, Np*K*Nfield*sizeof(double));
	/*Boundary edge facial integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int field = 0; field < 2; field++){
		for (int face = 0; face < BENe; face++){
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
			GetBEContributionToRHS(BERHSX + field*Np*K, LPDTBEfm + field*BENe*BENfp + face*BENfp, face, LocalPrimitiveDiffTermX + field*Np*K, \
				BEFToE, BEFToN1 + BENfp*face, Np, BENfp, BEFluxS + field*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BEfm + field*BENe*BENfp + face*BENfp, \
				BEfp + field*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BoundaryEdgeTau + face*IENfp, 1.0, AVBEfm + field*BENe*BENfp + face*IENfp, AVx + field*Np*K, \
				BEFluxM + field*BENe*BENfp + face*BENfp, BEJs + face*BENfp, BEMb);
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial (u,v)}{\partial y}$*/
			GetBEContributionToRHS(BERHSY + field*Np*K, LPDTBEfm + field*BENe*BENfp + face*BENfp, face, LocalPrimitiveDiffTermY + field*Np*K, \
				BEFToE, BEFToN1 + BENfp*face, Np, BENfp, BEFluxS + field*BENe*BENfp + face*BENfp, BEny + face*BENfp, BEfm + field*BENe*BENfp + face*BENfp, \
				BEfp + field*BENe*BENfp + face*BENfp, BEny + face*BENfp, BoundaryEdgeTau + face*IENfp, 1.0, AVBEfm + field*BENe*BENfp + face*IENfp, AVy + field*Np*K, \
				BEFluxM + field*BENe*BENfp + face*BENfp, BEJs + face*BENfp, BEMb);

		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int field = 2; field < Nfield; field++){
		for (int face = 0; face < BENe; face++){
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$,
			we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToRHS(BERHSX + field*Np*K, LPDTBEfm + field*BENe*BENfp + face*BENfp, face, LocalPrimitiveDiffTermX + field*Np*K, \
				BEFToE, BEFToN1 + BENfp*face, Np, BENfp, BEFluxS + field*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BEfm + field*BENe*BENfp + face*BENfp, \
				BEfp + field*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BoundaryEdgeTau + face*IENfp, Prantl, AVBEfm + field*BENe*BENfp + face*IENfp, AVx + field*Np*K, \
				BEFluxM + field*BENe*BENfp + face*BENfp, BEJs + face*BENfp, BEMb);
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$,
			we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToRHS(BERHSY + field*Np*K, LPDTBEfm + field*BENe*BENfp + face*BENfp, face, LocalPrimitiveDiffTermY + field*Np*K, \
				BEFToE, BEFToN1 + BENfp*face, Np, BENfp, BEFluxS + field*BENe*BENfp + face*BENfp, BEny + face*BENfp, BEfm + field*BENe*BENfp + face*BENfp, \
				BEfp + field*BENe*BENfp + face*BENfp, BEny + face*BENfp, BoundaryEdgeTau + face*IENfp, Prantl, AVBEfm + field*BENe*BENfp + face*IENfp, AVy + field*Np*K, \
				BEFluxM + field*BENe*BENfp + face*BENfp, BEJs + face*BENfp, BEMb);
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, IERHSX + field*Np*K + k*Np, Np);
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, BERHSX + field*Np*K + k*Np, Np);
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, IERHSY + field*Np*K + k*Np, Np);
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, BERHSY + field*Np*K + k*Np, Np);
		}
	}
	/*If Mellor and Blumberg's method in "Modeling vertical and horizontal diffusivities with the sigma coordinate system" is
	adopted, the following part about the facial integral should be ignored*/
	memset(IERHSX, 0, Np*K * 2 * sizeof(double));
	memset(IERHSY, 0, Np*K * 2 * sizeof(double));
	for (int face = 0; face < IENe; face++){
		int field = 0;
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial u}{\partial x}$*/
		GetIEContributionToRHS(IERHSX + field*Np*K, LPDTIEfm + field*IENe*IENfp + face*IENfp, LPDTIEfp + field*IENe*IENfp + face*IENfp, \
			face, LocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
			IEnx + face*IENfp, IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, IEnx + face*IENfp, InnerEdgeTau + face*IENfp, 1.0, \
			AVIEfm + field*IENe*IENfp + face*IENfp, AVIEfp + field*IENe*IENfp + face*IENfp, AVx + field*Np*K, IEFluxM + field*IENe*IENfp + face*IENfp, \
			IEFluxP + field*IENe*IENfp + face*IENfp, IEJs + face*IENfp, IEMb);
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial y}$, here $Q_x=\nv H\frac{\partial v}{\partial x}$*/
		GetIEContributionToRHS(IERHSY + field*Np*K, LPDTIEfm + field*IENe*IENfp + face*IENfp, LPDTIEfp + field*IENe*IENfp + face*IENfp, \
			face, LocalPrimitiveDiffTermX + (field+1)*Np*K, IEFToE, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
			IEny + face*IENfp, IEfm + (field+1)*IENe*IENfp + face*IENfp, IEfp + (field+1)*IENe*IENfp + face*IENfp, IEnx + face*IENfp, InnerEdgeTau + face*IENfp, 1.0, \
			AVIEfm + field*IENe*IENfp + face*IENfp, AVIEfp + field*IENe*IENfp + face*IENfp, AVx + (field+1)*Np*K, IEFluxM + field*IENe*IENfp + face*IENfp, \
			IEFluxP + field*IENe*IENfp + face*IENfp, IEJs + face*IENfp, IEMb);

		field = 1;
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial x}$, here $Q_y=\nv H\frac{\partial u}{\partial y}$*/
		GetIEContributionToRHS(IERHSX + field*Np*K, LPDTIEfm + field*IENe*IENfp + face*IENfp, LPDTIEfp + field*IENe*IENfp + face*IENfp, \
			face, LocalPrimitiveDiffTermY + (field-1)*Np*K, IEFToE, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
			IEnx + face*IENfp, IEfm + (field-1)*IENe*IENfp + face*IENfp, IEfp + (field-1)*IENe*IENfp + face*IENfp, IEny + face*IENfp, InnerEdgeTau + face*IENfp, 1.0, \
			AVIEfm + field*IENe*IENfp + face*IENfp, AVIEfp + field*IENe*IENfp + face*IENfp, AVy + (field-1)*Np*K, IEFluxM + field*IENe*IENfp + face*IENfp, \
			IEFluxP + field*IENe*IENfp + face*IENfp, IEJs + face*IENfp, IEMb);
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial v}{\partial y}$*/
		GetIEContributionToRHS(IERHSY + field*Np*K, LPDTIEfm + field*IENe*IENfp + face*IENfp, LPDTIEfp + field*IENe*IENfp + face*IENfp, \
			face, LocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp, IEFluxS + field*IENe*IENfp + face*IENfp, \
			IEny + face*IENfp, IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, IEny + face*IENfp, InnerEdgeTau + face*IENfp, 1.0, \
			AVIEfm + field*IENe*IENfp + face*IENfp, AVIEfp + field*IENe*IENfp + face*IENfp, AVy + field*Np*K, IEFluxM + field*IENe*IENfp + face*IENfp, \
			IEFluxP + field*IENe*IENfp + face*IENfp, IEJs + face*IENfp, IEMb);

	}
	memset(BERHSX, 0, Np*K * 2 * sizeof(double));
	memset(BERHSY, 0, Np*K * 2 * sizeof(double));
	for (int face = 0; face < BENe; face++){
			int field = 0;
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial u}{\partial x}$*/
			GetBEContributionToRHS(BERHSX + field*Np*K, LPDTBEfm + field*BENe*BENfp + face*BENfp, face, LocalPrimitiveDiffTermX + field*Np*K, \
				BEFToE, BEFToN1 + BENfp*face, Np, BENfp, BEFluxS + field*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BEfm + field*BENe*BENfp + face*BENfp, \
				BEfp + field*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BoundaryEdgeTau + face*BENfp, 1.0, AVBEfm + field*BENe*BENfp + face*IENfp, AVx + field*Np*K, \
				BEFluxM + field*BENe*BENfp + face*BENfp, BEJs + face*BENfp, BEMb);
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial y}$, here $Q_x=\nv H\frac{\partial v}{\partial x}$*/
			GetBEContributionToRHS(BERHSY + field*Np*K, LPDTBEfm + field*BENe*BENfp + face*BENfp, face, LocalPrimitiveDiffTermX + (field+1)*Np*K, \
				BEFToE, BEFToN1 + BENfp*face, Np, BENfp, BEFluxS + field*BENe*BENfp + face*BENfp, BEny + face*BENfp, BEfm + (field+1)*BENe*BENfp + face*BENfp, \
				BEfp + (field+1)*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BoundaryEdgeTau + face*BENfp, 1.0, AVBEfm + field*BENe*BENfp + face*BENfp, AVx + (field+1)*Np*K, \
				BEFluxM + field*BENe*BENfp + face*BENfp, BEJs + face*BENfp, BEMb);

			field = 1;
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial x}$, here $Q_y=\nv H\frac{\partial u}{\partial y}$*/
			GetBEContributionToRHS(BERHSX + field*Np*K, LPDTBEfm + field*BENe*BENfp + face*BENfp, face, LocalPrimitiveDiffTermY + (field - 1)*Np*K, \
				BEFToE, BEFToN1 + BENfp*face, Np, BENfp, BEFluxS + field*BENe*BENfp + face*BENfp, BEnx + face*BENfp, BEfm + (field - 1)*BENe*BENfp + face*BENfp, \
				BEfp + (field - 1)*BENe*BENfp + face*BENfp, BEny + face*BENfp, BoundaryEdgeTau + face*BENfp, 1.0, AVBEfm + field*BENe*BENfp + face*BENfp, AVy + (field - 1)*Np*K, \
				BEFluxM + field*BENe*BENfp + face*BENfp, BEJs + face*BENfp, BEMb);
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial v}{\partial y}$*/
			GetBEContributionToRHS(BERHSY + field*Np*K, LPDTBEfm + field*BENe*BENfp + face*BENfp, face, LocalPrimitiveDiffTermY + field*Np*K, \
				BEFToE, BEFToN1 + BENfp*face, Np, BENfp, BEFluxS + field*BENe*BENfp + face*BENfp, BEny + face*BENfp, BEfm + field*BENe*BENfp + face*BENfp, \
				BEfp + field*BENe*BENfp + face*BENfp, BEny + face*BENfp, BoundaryEdgeTau + face*BENfp, 1.0, AVBEfm + field*BENe*BENfp + face*BENfp, AVy + field*Np*K, \
				BEFluxM + field*BENe*BENfp + face*BENfp, BEJs + face*BENfp, BEMb);
		}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, IERHSX + field*Np*K + k*Np, Np);
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, BERHSX + field*Np*K + k*Np, Np);
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, IERHSY + field*Np*K + k*Np, Np);
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, BERHSY + field*Np*K + k*Np, Np);
		}
	}

	free(nv);
	free(variable);
	free(BEfp);
	free(AVx);
	free(AVy);
	free(Vx);
	free(TempVx);
	free(Vy);
	free(TempVy);
	free(IEfm);
	free(AVIEfm);
	free(IEfp);
	free(AVIEfp);
	free(IEFluxM);
	free(IEFluxP);
	free(IEFluxS);
	free(BEfm);
	free(AVBEfm);
	free(BEFluxM);
	free(BEFluxS);
	free(IERHSX);
	free(IERHSY);
	free(BERHSX);
	free(BERHSY);
	free(LocalPrimitiveDiffTermX);
	free(LocalPrimitiveDiffTermY);
	free(LPDTIEfm);
	free(LPDTIEfp);
	free(LPDTBEfm);
	free(InnerEdgeTau);
	free(BoundaryEdgeTau);
	free(IEnvfm);
	free(IEnvfp);
	free(BEnvfm);
}