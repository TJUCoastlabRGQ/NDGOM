#include "..\..\..\..\..\NdgMath\NdgMath.h"
#include "HorizontalDiffusion.h"
#include "..\..\..\..\..\NdgMath\NdgSWE.h"

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

	mxArray *TempInvM = mxGetField(prhs[9], 0, "invM");
	double *invM = mxGetPr(TempInvM);

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
	signed char *ftype = (signed char *)mxGetData(prhs[11]);
	double gra = mxGetScalar(prhs[12]);
	double *fext = mxGetPr(prhs[13]);


	/*Set the output right hand side*/
	const size_t NdimOut = 3;
	const mwSize dimOut[3] = { Np, K, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputRHS = mxGetPr(plhs[0]);
	memcpy(OutputRHS, InputRHS, Np*K*Nvar*sizeof(double));
	int Nfield;
	mxArray *TempOrder;
	double *variable = NULL, *TempBEfp = NULL, *BEfp = NULL, *TempBEfm = NULL;
	double *nv = malloc(Np*K*sizeof(double));
	int MeshType, Order;
	double *hu = NULL, *hv = NULL, *h = NULL, *z = NULL;
	double *zM = NULL, *zP = NULL;
	double *huM = NULL, *hvM = NULL, *hM = NULL;
	/*Allocate memory and calcualte variable u, v, $\theta$. Then impose the boudary condition and calculate u, v, $\theta$ defined over the ghost edge*/
	if (type == Two){
		/*For 2d shallow water problem, no horizontal diffusion terms are included in the governing equation for water depth $H$*/
		Nfield = Nvar - 1;
		/*Allocate memory for the original variable $u,v$ and $\theta$*/
		variable = malloc(Np*K*Nfield*sizeof(double));
		/*Allocate memory for the original variable over boundary edge*/
		BEfp = malloc(BENfp*BENe*Nfield*sizeof(double));
		MeshType = 2;
		TempOrder = mxGetField(prhs[9], 0, "N");
		Order = mxGetScalar(TempOrder);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
		for (int k = 0; k < K; k++){
			/*set the diffusion coefficient $H\nv$*/
			DotProduct(nv + k*Np, Tempnv + k*Np, fphys + (int)(varIndex[0] - 1)*Np*K + k*Np, Np);
			for (int field = 0; field < Nvar - 1; field++){
				/*For 2d shallow water problem, variable about height is organized as the first variable*/
				DotCriticalDivide(variable + field*Np*K + k*Np, \
					fphys + (int)(varIndex[field + 1] - 1)*Np*K + k*Np, Hcrit, \
					fphys + (int)(varIndex[0] - 1)*Np*K + k*Np, Np);
			}
		}

		zM = malloc(BENfp*BENe*sizeof(double));
		zP = malloc(BENfp*BENe*sizeof(double));
		/*Allocate memory for the local face value at boundary edge*/
		TempBEfp = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
		TempBEfm = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
		huM = TempBEfm, hvM = TempBEfm + BENe*BENfp, hM = TempBEfm + 2 * BENe*BENfp;

		h = fphys;
		hu = fphys;
		hv = fphys + 2 * Np*K;
		z = fphys + 3 * Np*K;

		/*Fetch variable fm and fp first, then impose boundary condition and conduct hydrostatic reconstruction.*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
		for (int face = 0; face < BENe; face++){
			NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
			FetchBoundaryEdgeFacialValue(huM + face*BENfp, hu, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hvM + face*BENfp, hv, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(zM + face*BENfp, z, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
			here 1 stands for the memory occupied by water depth h*/
			for (int field = 2; field < Nfield; field++){
				/*Index for water depth field needs to be considered*/
				FetchBoundaryEdgeFacialValue(TempBEfm + (field + 1)*BENe*BENfp + face*BENfp, \
					fphys + ((int)varIndex[field + 1] - 1)*Np*K, \
					BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
			}
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			ImposeBoundaryCondition(&gra, type, BEnx + face*BENfp, BEny + face*BENfp, TempBEfm + face*BENfp, TempBEfp + face*BENfp, \
				zM + face*BENfp, zP + face*BENfp, fext + face*BENfp, BENfp, Nfield + 1, BENe);
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			EvaluateHydroStaticReconstructValue(*Hcrit, TempBEfm + face*BENfp, TempBEfp + face*BENfp, zM + face*BENfp, zP + face*BENfp, BENfp, Nfield + 1, BENe);
			/*We divide the variable by water depth to get the original variable*/
			for (int field = 0; field < 2; field++){
				DotCriticalDivide(BEfp + field*BENfp*BENe + face*BENfp, \
					TempBEfp + field*BENfp*BENe + face*BENfp, Hcrit, \
					TempBEfp + 2 * BENfp*BENe + face*BENfp, BENfp);
			}

			for (int field = 2; field < Nfield; field++){
				DotCriticalDivide(BEfp + field*BENfp*BENe + face*BENfp, \
					TempBEfp + (field + 1)*BENfp*BENe + face*BENfp, Hcrit, \
					TempBEfp + 2 * BENfp*BENe + face*BENfp, BENfp);
			}
		}
	}
	else if (type == Three){
		Nfield = Nvar;
		/*Allocate memory for the original variable $u,v$ and $\theta$*/
		variable = malloc(Np*K*Nfield*sizeof(double));
		MeshType = 3;
		TempOrder = mxGetField(prhs[9], 0, "N");
		Order = (int)mxGetScalar(TempOrder);
		TempOrder = mxGetField(prhs[9], 0, "Nz");
		Order = max(Order, (int)mxGetScalar(TempOrder));
		
		
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
		for (int k = 0; k < K; k++){
			DotProduct(nv + k*Np, Tempnv + k*Np, fphys + 3 * Np*K + k*Np, Np);
			for (int field = 0; field < Nvar; field++){
				DotCriticalDivide(variable + field*Np*K + k*Np, \
					fphys + (int)(varIndex[field] - 1)*Np*K + k*Np, Hcrit, \
					fphys + 3*Np*K + k*Np, Np);//For 3d shallow water problem, variable about height is organized as the forth variable
			}
		}

		zM = malloc(BENfp*BENe*sizeof(double));
		zP = malloc(BENfp*BENe*sizeof(double));
		TempBEfp = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
		TempBEfm = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
		/*Allocate memory for the original variable over boundary edge*/
		BEfp = malloc(BENfp*BENe*Nfield*sizeof(double));
		huM = TempBEfm, hvM = TempBEfm + BENe*BENfp, hM = TempBEfm + 2 * BENe*BENfp;
		hu = fphys;
		hv = fphys + Np*K;
		h = fphys + 3 * Np*K;
		z = fphys + 5 * Np*K;
		/*Fetch variable fm and fp first, then impose boundary condition and conduct hydrostatic reconstruction.
		Finally, calculate local flux term, adjacent flux term and numerical flux term*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
		for (int face = 0; face < BENe; face++){
			NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
			FetchBoundaryEdgeFacialValue(huM + face*BENfp, hu, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hvM + face*BENfp, hv, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(zM + face*BENfp, z, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
			here 1 stands for the memory occupied by water depth h*/
			for (int field = 2; field < Nvar; field++){
				FetchBoundaryEdgeFacialValue(TempBEfm + (field + 1)*BENe*BENfp + face*BENfp, \
					fphys + ((int)varIndex[field] - 1)*Np*K, \
					BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
			}
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			ImposeBoundaryCondition(&gra, type, BEnx + face*BENfp, BEny + face*BENfp, TempBEfm + face*BENfp, TempBEfp + face*BENfp, \
				zM + face*BENfp, zP + face*BENfp, fext + face*BENfp, BENfp, Nfield + 1, BENe);
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			EvaluateHydroStaticReconstructValue(*Hcrit, TempBEfm + face*BENfp, TempBEfp + face*BENfp, zM + face*BENfp, zP + face*BENfp, BENfp, Nfield + 1, BENe);
			/*We divide the variable by water depth to get the original variable*/
			for (int field = 0; field < 2; field++){
				DotCriticalDivide(BEfp + field*BENfp*BENe + face*BENfp, \
					TempBEfp + field*BENfp*BENe + face*BENfp, Hcrit, \
					TempBEfp + 2 * BENfp*BENe + face*BENfp, BENfp);
			}

			for (int field = 2; field < Nfield; field++){
				DotCriticalDivide(BEfp + field*BENfp*BENe + face*BENfp, \
					TempBEfp + (field + 1)*BENfp*BENe + face*BENfp, Hcrit, \
					TempBEfp + 2 * BENfp*BENe + face*BENfp, BENfp);
			}

		}
	}

	free(TempBEfp);
	free(TempBEfm);
	free(zM);
	free(zP);

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
	/*Allocate memory for the local face value at boundary edge*/
	double *BEfm = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at inner edge*/
	double *IEFluxS = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at boundary edge*/
	double *AVBEfm = malloc(BENe*BENfp*Nfield*sizeof(double));
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
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			/*$\bold{r_x}\cdot (Dr*u(v,\theta))+\bold{s_x}\cdot (Ds*u(v,\theta))$*/			
			GetVolumnIntegral2d(Vx + field*Np*K + k*Np, TempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, variable + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
			/*$\bold{r_y}\cdot (Dr*u(v,\theta))+\bold{s_y}\cdot (Ds*u(v,\theta))$*/
			GetVolumnIntegral2d(Vy + field*Np*K + k*Np, TempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, variable + field*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
 		}
	}
/*Allocate memory for local primitive diffusion term in both x and y direction, $\nv\Nabla u(v,\theta)$*/
/*This part is used because we need to extract the trace of the local derivative operator to compute the numerical flux*/
	double *LocalPrimitiveDiffTermX = malloc(Np*K*Nfield*sizeof(double));
	double *LocalPrimitiveDiffTermY = malloc(Np*K*Nfield*sizeof(double));
	double *LPDTIEfm = malloc(IENfp*IENe*Nfield*sizeof(double));
	double *LPDTIEfp = malloc(IENfp*IENe*Nfield*sizeof(double));
	double *LPDTBEfm = malloc(BENfp*BENe*Nfield*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			DotProduct(LocalPrimitiveDiffTermX + field*Np*K + k*Np, Vx + field*Np*K + k*Np, nv + k*Np, Np);
			DotProduct(LocalPrimitiveDiffTermY + field*Np*K + k*Np, Vy + field*Np*K + k*Np, nv + k*Np, Np);
		}
		/*for substance transport, prantl number is considered*/
		for (int field = 2; field < Nfield; field++){
			DotDivideByConstant(LocalPrimitiveDiffTermX + field*Np*K + k*Np, LocalPrimitiveDiffTermX + field*Np*K + k*Np, Prantl, Np);
			DotDivideByConstant(LocalPrimitiveDiffTermY + field*Np*K + k*Np, LocalPrimitiveDiffTermY + field*Np*K + k*Np, Prantl, Np);
		}
	}

/*Inner edge facial integral part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nfield; field++){
		for (int face = 0; face < IENe; face++){
			FetchInnerEdgeFacialValue(IEfm + field*IENe*IENfp + face*IENfp, IEfp + field*IENe*IENfp + face*IENfp, variable + field*Np*K, IEFToE + 2 * face, \
				IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
			/*Inner edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$*/
			GetIEContributionToAuxialaryVariable(IERHSX + field*Np*K, face, IENe, IENfp, field, IEfm, IEfp, IEFToE, IEFToN1, IEFToN2, Np, IEFluxM, IEFluxP, IEFluxS, IEnx, IEMb, IEJs);
			/*Inner edge contribution to RHSY of $\frac{\partial u(v,\theta)}{\partial y}$*/
			GetIEContributionToAuxialaryVariable(IERHSY + field*Np*K, face, IENe, IENfp, field, IEfm, IEfp, IEFToE, IEFToN1, IEFToN2, Np, IEFluxM, IEFluxP, IEFluxS, IEny, IEMb, IEJs);
		}
	}
	/*Boundary edge facial integral part*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nfield; field++){
		for (int face = 0; face < BENe; face++){
			FetchBoundaryEdgeFacialValue(BEfm + field*BENe*BENfp + face*BENfp, variable + field*Np*K, BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
			/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$, we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToAuxialaryVariable(BERHSX + field*Np*K, face, BENe, BENfp, field, BEfm, BEfp, BEFToE, BEFToN1, Np, BEFluxM, BEFluxS, BEnx, BEMb, BEJs);

			/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial y}$, we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToAuxialaryVariable(BERHSY + field*Np*K, face, BENe, BENfp, field, BEfm, BEfp, BEFToE, BEFToN1, Np, BEFluxM, BEFluxS, BEny, BEMb, BEJs);
		}
	}


	double *TempFacialIntegralX = malloc(Np*K*Nfield*sizeof(double));
	double *TempFacialIntegralY = malloc(Np*K*Nfield*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nfield; field++){

			MultiEdgeContributionByLiftOperator(IERHSX + field*Np*K + k*Np, TempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J+k*Np, Np);

			MultiEdgeContributionByLiftOperator(IERHSY + field*Np*K + k*Np, TempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nfield; field++){

			MultiEdgeContributionByLiftOperator(BERHSX + field*Np*K + k*Np, TempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);

			MultiEdgeContributionByLiftOperator(BERHSY + field*Np*K + k*Np, TempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}



	/*Next, sum contribution from volume integral, inner edge contribution and boundary edge contribution into AVx and Avy*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
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
#pragma omp parallel for num_threads(omp_get_max_threads())
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
	/*
	plhs[1] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *px = mxGetPr(plhs[1]);
	memcpy(px, AVx, Np*K*Nfield*sizeof(double));
	plhs[2] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *py = mxGetPr(plhs[2]);
	memcpy(py, AVy, Np*K*Nfield*sizeof(double));
	*/
	double *InnerEdgeTau = malloc(IENe*IENfp*sizeof(double));
	double *BoundaryEdgeTau = malloc(BENe*BENfp*sizeof(double));
	double *IEnvfm = malloc(IENe*IENfp*sizeof(double));
	double *IEnvfp = malloc(IENe*IENfp*sizeof(double));
	double *BEnvfm = malloc(BENe*BENfp*sizeof(double));
	/*Calculate the contribution to the right hand side due to the auxialary variable AVx and AVy with IPDG.*/
	/*Calculate the penalty parameter $\tau$ first, this parameter is calculated as $\tau=\frac{(N+1)(N+d)}{d}\frac{n_0}{2}\frac{A}{V}\nv$*/
	/*Inner edge first*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < IENe; face++){
		FetchInnerEdgeFacialValue(IEnvfm + face*IENfp, IEnvfp + face*IENfp, nv, \
			IEFToE + face * 2, IEFToN1 + face*IENfp, IEFToN2 + face*IENfp, Np, IENfp);
		double localRatio = IELAV[face] / LAV[(int)IEFToE[face * 2]-1];
		double adjacentRatio = IELAV[face] / LAV[(int)IEFToE[face * 2 + 1]-1];
		for (int p = 0; p < IENfp; p++){
			InnerEdgeTau[face*IENfp + p] = max(localRatio*(Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * IEnvfm[face*IENfp + p], \
				adjacentRatio*(Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * IEnvfp[face*IENfp + p]);
		}
	}
	/*Boundary edge next*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < BENe; face++){
		FetchBoundaryEdgeFacialValue(BEnvfm + face*BENfp, nv, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
		double localRatio = BELAV[face] / LAV[(int)BEFToE[face * 2]-1];
		for (int p = 0; p < BENfp; p++){
			BoundaryEdgeTau[face*BENfp + p] = localRatio * (Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * BEnvfm[face*IENfp + p];
		}
	}

	/*Volume integral of the second order operator*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		/*Prantl number is not considered here, since we have considered this parameter when calculating the auxialary variable,
		it is not need anymore in the volumn integral part when considering the contribution of second order operator to the RHS*/
		for (int field = 0; field < Nfield; field++){
			/*$\bold{r_x}\cdot (Dr*Q_x)+\bold{s_x}\cdot (Ds*Q_x)$*/
			GetVolumnIntegral2d(Vx + field*Np*K + k*Np, TempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, AVx + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
			/*$\bold{r_y}\cdot (Dr*Q_y)+\bold{s_y}\cdot (Ds*Q_y)$*/
			GetVolumnIntegral2d(Vy + field*Np*K + k*Np, TempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
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
/************************************************************************************************************************************/
		/*If Mellor and Blumberg's method in "Modeling vertical and horizontal diffusivities with the sigma coordinate system" is
		adopted, the next part should be ignored*/
/************************************************************************************************************************************/
		int field = 0;
		/*$\bold{r_x}\cdot (Dr*Q_x)+\bold{s_x}\cdot (Ds*Q_x)$, with $*Q_x=\nv H\frac{\partial u}{\partial x}*$*/
		GetVolumnIntegral2d(Vx + field*Np*K + k*Np, TempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, AVx + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
		/*$\bold{r_y}\cdot (Dr*Q_x)+\bold{s_y}\cdot (Ds*Q_x)$, with $*Q_x=\nv H\frac{\partial v}{\partial x}*$*/
		GetVolumnIntegral2d(Vy + field*Np*K + k*Np, TempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, AVx + (field + 1)*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		field = 1;
		/*$\bold{r_x}\cdot (Dr*Q_y)+\bold{s_x}\cdot (Ds*Q_y)$, with $*Q_y=\nv H\frac{\partial u}{\partial y}*$*/
		GetVolumnIntegral2d(Vx + field*Np*K + k*Np, TempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, AVy + (field - 1)*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
		/*$\bold{r_y}\cdot (Dr*Q_y)+\bold{s_y}\cdot (Ds*Q_y)$, with $*Q_y=\nv H\frac{\partial v}{\partial y}*$*/
		GetVolumnIntegral2d(Vy + field*Np*K + k*Np, TempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
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
/***********************************************************************************************************************************************/
	/*Reset all the data contained in IERHSX and IERHSY to zero, becaused these space contains the data left when calculating the auxialary variable*/
	memset(IERHSX, 0, Np*K*Nfield*sizeof(double));
	memset(IERHSY, 0, Np*K*Nfield*sizeof(double));
	/*Surface integral of the second order operator*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < 2; field++){
		for (int face = 0; face < IENe; face++){
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
			GetIEContributionToRHS(IERHSX + field*Np*K, LPDTIEfm, LPDTIEfp, \
				face, LocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToN1, IEFToN2, Np, IENe, IENfp, field, IEFluxS, \
				IEnx, IEfm + field*IENe*IENfp, IEfp + field*IENe*IENfp, IEnx, InnerEdgeTau, 1.0, \
				AVIEfm, AVIEfp, AVx + field*Np*K, IEFluxM, \
				IEFluxP, IEJs, IEMb);
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial y}$*/
			GetIEContributionToRHS(IERHSY + field*Np*K, LPDTIEfm, LPDTIEfp, \
				face, LocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToN1, IEFToN2, Np, IENe, IENfp, field, IEFluxS, \
				IEny, IEfm + field*IENe*IENfp, IEfp + field*IENe*IENfp, IEny, InnerEdgeTau, 1.0, \
				AVIEfm, AVIEfp, AVy + field*Np*K, IEFluxM, \
				IEFluxP, IEJs, IEMb);

		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 2; field < Nfield; field++){
		for (int face = 0; face < IENe; face++){
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$*/
			GetIEContributionToRHS(IERHSX + field*Np*K, LPDTIEfm, LPDTIEfp, \
				face, LocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToN1, IEFToN2, Np, IENe, IENfp, field, IEFluxS, \
				IEnx, IEfm + field*IENe*IENfp, IEfp + field*IENe*IENfp, IEnx, InnerEdgeTau, Prantl, \
				AVIEfm, AVIEfp, AVx + field*Np*K, IEFluxM, \
				IEFluxP, IEJs, IEMb);
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$*/
			GetIEContributionToRHS(IERHSY + field*Np*K, LPDTIEfm, LPDTIEfp, \
				face, LocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToN1, IEFToN2, Np, IENe, IENfp, field, IEFluxS, \
				IEny, IEfm + field*IENe*IENfp, IEfp + field*IENe*IENfp, IEny, InnerEdgeTau, Prantl, \
				AVIEfm, AVIEfp, AVy + field*Np*K, IEFluxM, \
				IEFluxP, IEJs, IEMb);

		}
	}
	memset(BERHSX, 0, Np*K*Nfield*sizeof(double));
	memset(BERHSY, 0, Np*K*Nfield*sizeof(double));
	/*Boundary edge facial integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < 2; field++){
		for (int face = 0; face < BENe; face++){
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
			GetBEContributionToRHS(BERHSX + field*Np*K, LPDTBEfm, face, LocalPrimitiveDiffTermX + field*Np*K, \
				BEFToE, BEFToN1, Np, BENe, BENfp, field, BEFluxS, BEnx, BEfm + field*BENe*BENfp, \
				BEfp + field*BENe*BENfp, BEnx, BoundaryEdgeTau, 1.0, AVBEfm, AVx + field*Np*K, \
				BEFluxM, BEJs, BEMb);
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial (u,v)}{\partial y}$*/
			GetBEContributionToRHS(BERHSY + field*Np*K, LPDTBEfm, face, LocalPrimitiveDiffTermY + field*Np*K, \
				BEFToE, BEFToN1, Np, BENe, BENfp, field, BEFluxS, BEny, BEfm + field*BENe*BENfp, \
				BEfp + field*BENe*BENfp, BEny, BoundaryEdgeTau, 1.0, AVBEfm, AVy + field*Np*K, \
				BEFluxM, BEJs, BEMb);

		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 2; field < Nfield; field++){
		for (int face = 0; face < BENe; face++){
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$,
			we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToRHS(BERHSX + field*Np*K, LPDTBEfm, face, LocalPrimitiveDiffTermX + field*Np*K, \
				BEFToE, BEFToN1, Np, BENe, BENfp, field, BEFluxS, BEnx, BEfm + field*BENe*BENfp, \
				BEfp + field*BENe*BENfp, BEnx, BoundaryEdgeTau, Prantl, AVBEfm, AVx + field*Np*K, \
				BEFluxM, BEJs, BEMb);
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$,
			we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToRHS(BERHSY + field*Np*K, LPDTBEfm, face, LocalPrimitiveDiffTermY + field*Np*K, \
				BEFToE, BEFToN1, Np, BENe, BENfp, field, BEFluxS, BEny, BEfm + field*BENe*BENfp, \
				BEfp + field*BENe*BENfp, BEny, BoundaryEdgeTau, Prantl, AVBEfm, AVy + field*Np*K, \
				BEFluxM, BEJs, BEMb);
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nfield; field++){

			MultiEdgeContributionByLiftOperator(IERHSX + field*Np*K + k*Np, TempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);

			MultiEdgeContributionByLiftOperator(IERHSY + field*Np*K + k*Np, TempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nfield; field++){

			MultiEdgeContributionByLiftOperator(BERHSX + field*Np*K + k*Np, TempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);

			MultiEdgeContributionByLiftOperator(BERHSY + field*Np*K + k*Np, TempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
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
	/****************************************************************************************************************************/
	/*If Mellor and Blumberg's method in "Modeling vertical and horizontal diffusivities with the sigma coordinate system" is
	adopted, the following part about the facial integral should be ignored*/
	memset(IERHSX, 0, Np*K * 2 * sizeof(double));
	memset(IERHSY, 0, Np*K * 2 * sizeof(double));
	for (int face = 0; face < IENe; face++){
		int field = 0;
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial u}{\partial x}$*/
		GetIEContributionToRHS(IERHSX + field*Np*K, LPDTIEfm, LPDTIEfp, \
			face, LocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToN1, IEFToN2, Np, IENe, IENfp, field, IEFluxS, \
			IEnx, IEfm + field*IENe*IENfp, IEfp + field*IENe*IENfp, IEnx, InnerEdgeTau, 1.0, \
			AVIEfm, AVIEfp, AVx + field*Np*K, IEFluxM, \
			IEFluxP, IEJs, IEMb);
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial y}$, here $Q_x=\nv H\frac{\partial v}{\partial x}$*/
		GetIEContributionToRHS(IERHSY + field*Np*K, LPDTIEfm, LPDTIEfp, \
			face, LocalPrimitiveDiffTermX + (field + 1)*Np*K, IEFToE, IEFToN1, IEFToN2, Np, IENe, IENfp, field, IEFluxS, \
			IEny, IEfm + (field + 1)*IENe*IENfp, IEfp + (field + 1)*IENe*IENfp, IEnx, InnerEdgeTau, 1.0, \
			AVIEfm, AVIEfp, AVx + (field + 1)*Np*K, IEFluxM, \
			IEFluxP, IEJs, IEMb);

		field = 1;
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial x}$, here $Q_y=\nv H\frac{\partial u}{\partial y}$*/
		GetIEContributionToRHS(IERHSX + field*Np*K, LPDTIEfm, LPDTIEfp, \
			face, LocalPrimitiveDiffTermY + (field - 1)*Np*K, IEFToE, IEFToN1, IEFToN2, Np, IENe, IENfp, field, IEFluxS, \
			IEnx, IEfm + (field - 1)*IENe*IENfp, IEfp + (field - 1)*IENe*IENfp, IEny, InnerEdgeTau, 1.0, \
			AVIEfm, AVIEfp, AVy + (field - 1)*Np*K, IEFluxM, \
			IEFluxP, IEJs, IEMb);
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial v}{\partial y}$*/
		GetIEContributionToRHS(IERHSY + field*Np*K, LPDTIEfm, LPDTIEfp, \
			face, LocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToN1, IEFToN2, Np, IENe, IENfp, field, IEFluxS, \
			IEny, IEfm + field*IENe*IENfp, IEfp + field*IENe*IENfp, IEny, InnerEdgeTau, 1.0, \
			AVIEfm, AVIEfp, AVy + field*Np*K, IEFluxM, \
			IEFluxP, IEJs, IEMb);

	}
	memset(BERHSX, 0, Np*K * 2 * sizeof(double));
	memset(BERHSY, 0, Np*K * 2 * sizeof(double));
	for (int face = 0; face < BENe; face++){
		int field = 0;
		/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial u}{\partial x}$*/
		GetBEContributionToRHS(BERHSX + field*Np*K, LPDTBEfm, face, LocalPrimitiveDiffTermX + field*Np*K, \
			BEFToE, BEFToN1, Np, BENe, BENfp, field, BEFluxS, BEnx, BEfm + field*BENe*BENfp, \
			BEfp + field*BENe*BENfp, BEnx, BoundaryEdgeTau, 1.0, AVBEfm, AVx + field*Np*K, \
			BEFluxM, BEJs, BEMb);
		/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial y}$, here $Q_x=\nv H\frac{\partial v}{\partial x}$*/
		GetBEContributionToRHS(BERHSY + field*Np*K, LPDTBEfm, face, LocalPrimitiveDiffTermX + (field + 1)*Np*K, \
			BEFToE, BEFToN1, Np, BENe, BENfp, field, BEFluxS, BEny, BEfm + (field + 1)*BENe*BENfp, \
			BEfp + (field + 1)*BENe*BENfp, BEnx, BoundaryEdgeTau, 1.0, AVBEfm, AVx + (field + 1)*Np*K, \
			BEFluxM, BEJs, BEMb);

		field = 1;
		/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial x}$, here $Q_y=\nv H\frac{\partial u}{\partial y}$*/
		GetBEContributionToRHS(BERHSX + field*Np*K, LPDTBEfm, face, LocalPrimitiveDiffTermY + (field - 1)*Np*K, \
			BEFToE, BEFToN1, Np, BENe, BENfp, field, BEFluxS, BEnx, BEfm + (field - 1)*BENe*BENfp, \
			BEfp + (field - 1)*BENe*BENfp, BEny, BoundaryEdgeTau, 1.0, AVBEfm, AVy + (field - 1)*Np*K, \
			BEFluxM, BEJs, BEMb);
		/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial v}{\partial y}$*/
		GetBEContributionToRHS(BERHSY + field*Np*K, LPDTBEfm, face, LocalPrimitiveDiffTermY + field*Np*K, \
			BEFToE, BEFToN1, Np, BENe, BENfp, field, BEFluxS, BEny, BEfm + field*BENe*BENfp, \
			BEfp + field*BENe*BENfp, BEny, BoundaryEdgeTau, 1.0, AVBEfm, AVy + field*Np*K, \
			BEFluxM, BEJs, BEMb);
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < 2; field++){

			MultiEdgeContributionByLiftOperator(IERHSX + field*Np*K + k*Np, TempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);

			MultiEdgeContributionByLiftOperator(IERHSY + field*Np*K + k*Np, TempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < 2; field++){

			MultiEdgeContributionByLiftOperator(BERHSX + field*Np*K + k*Np, TempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);

			MultiEdgeContributionByLiftOperator(BERHSY + field*Np*K + k*Np, TempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
/******************************************************************************************************************************/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
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
	free(TempFacialIntegralX);
	free(TempFacialIntegralY);
}