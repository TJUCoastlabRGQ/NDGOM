#include "..\..\..\..\..\NdgMath\NdgMath.h"
#include "HorizontalDiffusion.h"
#include "..\..\..\..\..\NdgMath\NdgSWE.h"
#include "..\..\..\..\..\NdgMath\NdgMemory.h"

extern double *HorDiffnv, *HorDiffvariable, *HorDiffBEfp, *HorDiffzM, \
*HorDiffzP, *HorDiffTempBEfp, *HorDiffTempBEfm, *HorDiffAVx, \
*HorDiffAVy, *HorDiffVx, *HorDiffTempVx, *HorDiffVy, *HorDiffTempVy, \
*HorDiffIEfm, *HorDiffAVIEfm, *HorDiffIEfp, *HorDiffAVIEfp, *HorDiffIEFluxM, \
*HorDiffIEFluxP, *HorDiffBEfm, *HorDiffIEFluxS, *HorDiffAVBEfm, *HorDiffBEFluxM, \
*HorDiffBEFluxS, *HorDiffERHSX, *HorDiffERHSY, *HorDiffLocalPrimitiveDiffTermX, \
*HorDiffLocalPrimitiveDiffTermY, *HorDiffLPDTIEfm, *HorDiffLPDTIEfp, *HorDiffLPDTBEfm, \
*HorDiffTempFacialIntegralX, *HorDiffTempFacialIntegralY, *HorDiffInnerEdgeTau, \
*HorDiffBoundaryEdgeTau, *HorDiffIEnvfm, *HorDiffIEnvfp, *HorDiffBEnvfm;

extern char *HorDiffInitialized;

void MyExit()
{
	if (HorDiffInitialized == "true"){
		HorizDiffMemoryDeAllocation();
		HorDiffInitialized = "false";
	}
	return;
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
	mxArray *TempIEFToF = mxGetField(prhs[3], 0, "FToF");
	double *IEFToF = mxGetPr(TempIEFToF);    
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
	mxArray *TempBEFToF = mxGetField(prhs[4], 0, "FToF");
	double *BEFToF = mxGetPr(TempBEFToF);    
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
	int tempNface = (int)mxGetScalar(TempNface);
    int Nface;

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
	int MeshType, Order;
	double *hu = NULL, *hv = NULL, *h = NULL, *z = NULL;
	double *huM = NULL, *hvM = NULL, *hM = NULL;

	if (HorDiffInitialized = "false"){

		HorizDiffMemoryAllocation( type, Np, K, Nvar, tempNface, BENfp, BENe, IENfp, IENe);

	}

	/*Allocate memory and calcualte variable u, v, $\theta$. Then impose the boudary condition and calculate u, v, $\theta$ defined over the ghost edge*/
	if (type == Two){
		/*For 2d shallow water problem, no horizontal diffusion terms are included in the governing equation for water depth $H$*/
		Nfield = Nvar - 1;
        /*For 2d shallow water problem, the face number is equal to TempNface, since there is no surface edge and bottom edge*/
        Nface = tempNface;

		MeshType = 2;
		TempOrder = mxGetField(prhs[9], 0, "N");
		Order = mxGetScalar(TempOrder);  

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
		for (int k = 0; k < K; k++){
			/*set the diffusion coefficient $H\nv$*/
			DotProduct(HorDiffnv + k*Np, Tempnv + k*Np, fphys + (int)(varIndex[0] - 1)*Np*K + k*Np, Np);
			for (int field = 0; field < Nvar - 1; field++){
				/*For 2d shallow water problem, variable about height is organized as the first variable*/
				DotCriticalDivide(HorDiffvariable + field*Np*K + k*Np, \
					fphys + (int)(varIndex[field + 1] - 1)*Np*K + k*Np, Hcrit, \
					fphys + (int)(varIndex[0] - 1)*Np*K + k*Np, Np);
			}
		}
		huM = HorDiffTempBEfm, hvM = HorDiffTempBEfm + BENe*BENfp, hM = HorDiffTempBEfm + 2 * BENe*BENfp;

		h = fphys;
		hu = fphys;
		hv = fphys + 2 * Np*K;
		z = fphys + 3 * Np*K;

		/*Fetch variable fm and fp first, then impose boundary condition and conduct hydrostatic reconstruction.*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
		for (int face = 0; face < BENe; face++){
			NdgEdgeType facetype = (NdgEdgeType)ftype[face];  // boundary condition
			FetchBoundaryEdgeFacialValue(huM + face*BENfp, hu, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hvM + face*BENfp, hv, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(HorDiffzM + face*BENfp, z, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
			here 1 stands for the memory occupied by water depth h*/
			for (int field = 2; field < Nfield; field++){
				/*Index for water depth field needs to be considered*/
				FetchBoundaryEdgeFacialValue(HorDiffTempBEfm + (field + 1)*BENe*BENfp + face*BENfp, \
					fphys + ((int)varIndex[field + 1] - 1)*Np*K, \
					BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
			}
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			ImposeBoundaryCondition(&gra, facetype, BEnx + face*BENfp, BEny + face*BENfp, HorDiffTempBEfm + face*BENfp, HorDiffTempBEfp + face*BENfp, \
				HorDiffzM + face*BENfp, HorDiffzP + face*BENfp, fext + face*BENfp, BENfp, Nfield + 1, BENe);
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			EvaluateHydroStaticReconstructValue(*Hcrit, HorDiffTempBEfm + face*BENfp, HorDiffTempBEfp + face*BENfp, HorDiffzM + face*BENfp, HorDiffzP + face*BENfp, BENfp, Nfield + 1, BENe);
			/*We divide the variable by water depth to get the original variable*/
			for (int field = 0; field < 2; field++){
				DotCriticalDivide(HorDiffBEfp + field*BENfp*BENe + face*BENfp, \
					HorDiffTempBEfp + field*BENfp*BENe + face*BENfp, Hcrit, \
					HorDiffTempBEfp + 2 * BENfp*BENe + face*BENfp, BENfp);
			}

			for (int field = 2; field < Nfield; field++){
				DotCriticalDivide(HorDiffBEfp + field*BENfp*BENe + face*BENfp, \
					HorDiffTempBEfp + (field + 1)*BENfp*BENe + face*BENfp, Hcrit, \
					HorDiffTempBEfp + 2 * BENfp*BENe + face*BENfp, BENfp);
			}
		}
	}
	else if (type == Three){
		Nfield = Nvar;
        /*For 3d shallow water problem, the face number is equal to TempNface - 2, since there 
         * is surface edge and bottom edge is not considered for horizontal diffusion term*/
        Nface = tempNface - 2;        
		MeshType = 3;
		TempOrder = mxGetField(prhs[9], 0, "N");
		Order = (int)mxGetScalar(TempOrder);
		TempOrder = mxGetField(prhs[9], 0, "Nz");
		Order = max(Order, (int)mxGetScalar(TempOrder));	
		
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
		for (int k = 0; k < K; k++){
			DotProduct(HorDiffnv + k*Np, Tempnv + k*Np, fphys + 3 * Np*K + k*Np, Np);
			for (int field = 0; field < Nvar; field++){
				DotCriticalDivide(HorDiffvariable + field*Np*K + k*Np, \
					fphys + (int)(varIndex[field] - 1)*Np*K + k*Np, Hcrit, \
					fphys + 3*Np*K + k*Np, Np);//For 3d shallow water problem, variable about height is organized as the forth variable
			}
		}
		huM = HorDiffTempBEfm, hvM = HorDiffTempBEfm + BENe*BENfp, hM = HorDiffTempBEfm + 2 * BENe*BENfp;
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
			NdgEdgeType facetype = (NdgEdgeType)ftype[face];  // boundary condition
			FetchBoundaryEdgeFacialValue(huM + face*BENfp, hu, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hvM + face*BENfp, hv, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(HorDiffzM + face*BENfp, z, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
			/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
			here 1 stands for the memory occupied by water depth h*/
			for (int field = 2; field < Nvar; field++){
				FetchBoundaryEdgeFacialValue(HorDiffTempBEfm + (field + 1)*BENe*BENfp + face*BENfp, \
					fphys + ((int)varIndex[field] - 1)*Np*K, \
					BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
			}
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			ImposeBoundaryCondition(&gra, facetype, BEnx + face*BENfp, BEny + face*BENfp, HorDiffTempBEfm + face*BENfp, HorDiffTempBEfp + face*BENfp, \
				HorDiffzM + face*BENfp, HorDiffzP + face*BENfp, fext + face*BENfp, BENfp, Nfield + 1, BENe);
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			EvaluateHydroStaticReconstructValue(*Hcrit, HorDiffTempBEfm + face*BENfp, HorDiffTempBEfp + face*BENfp, HorDiffzM + face*BENfp, HorDiffzP + face*BENfp, BENfp, Nfield + 1, BENe);
			/*We divide the variable by water depth to get the original variable*/
			for (int field = 0; field < 2; field++){
				DotCriticalDivide(HorDiffBEfp + field*BENfp*BENe + face*BENfp, \
					HorDiffTempBEfp + field*BENfp*BENe + face*BENfp, Hcrit, \
					HorDiffTempBEfp + 2 * BENfp*BENe + face*BENfp, BENfp);
			}
           /*Water depth is stored as the third variable*/
			for (int field = 2; field < Nfield; field++){
				DotCriticalDivide(HorDiffBEfp + field*BENfp*BENe + face*BENfp, \
					HorDiffTempBEfp + (field + 1)*BENfp*BENe + face*BENfp, Hcrit, \
					HorDiffTempBEfp + 2 * BENfp*BENe + face*BENfp, BENfp);
			}

		}
	}

	memset(HorDiffERHSX, 0, Np*K*Nfield*Nface*sizeof(double));
	memset(HorDiffERHSY, 0, Np*K*Nfield*Nface*sizeof(double));
   
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
			GetVolumnIntegral2d(HorDiffVx + field*Np*K + k*Np, HorDiffTempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, HorDiffvariable + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
			/*$\bold{r_y}\cdot (Dr*u(v,\theta))+\bold{s_y}\cdot (Ds*u(v,\theta))$*/
			GetVolumnIntegral2d(HorDiffVy + field*Np*K + k*Np, HorDiffTempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, HorDiffvariable + field*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
 		}
	}
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			DotProduct(HorDiffLocalPrimitiveDiffTermX + field*Np*K + k*Np, HorDiffVx + field*Np*K + k*Np, HorDiffnv + k*Np, Np);
			DotProduct(HorDiffLocalPrimitiveDiffTermY + field*Np*K + k*Np, HorDiffVy + field*Np*K + k*Np, HorDiffnv + k*Np, Np);
		}
		/*for substance transport, prantl number is considered*/
		for (int field = 2; field < Nfield; field++){
			DotDivideByConstant(HorDiffLocalPrimitiveDiffTermX + field*Np*K + k*Np, HorDiffLocalPrimitiveDiffTermX + field*Np*K + k*Np, Prantl, Np);
			DotDivideByConstant(HorDiffLocalPrimitiveDiffTermY + field*Np*K + k*Np, HorDiffLocalPrimitiveDiffTermY + field*Np*K + k*Np, Prantl, Np);
		}
	}

/*Inner edge facial integral part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int face = 0; face < IENe; face++){
        for (int field = 0; field < Nfield; field++){
			FetchInnerEdgeFacialValue(HorDiffIEfm + field*IENe*IENfp + face*IENfp, HorDiffIEfp + field*IENe*IENfp + face*IENfp, HorDiffvariable + field*Np*K, IEFToE + 2 * face, \
				IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
			/*Inner edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$*/
			GetIEContributionToAuxialaryVariable(HorDiffERHSX + field*Np*K*Nface, face, IENe, IENfp, field, HorDiffIEfm, HorDiffIEfp, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, HorDiffIEFluxM, HorDiffIEFluxP, HorDiffIEFluxS, IEnx, IEMb, IEJs);
			/*Inner edge contribution to RHSY of $\frac{\partial u(v,\theta)}{\partial y}$*/
			GetIEContributionToAuxialaryVariable(HorDiffERHSY + field*Np*K*Nface, face, IENe, IENfp, field, HorDiffIEfm, HorDiffIEfp, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, HorDiffIEFluxM, HorDiffIEFluxP, HorDiffIEFluxS, IEny, IEMb, IEJs);
		}
	}
	/*Boundary edge facial integral part*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int face = 0; face < BENe; face++){
        for (int field = 0; field < Nfield; field++){
			FetchBoundaryEdgeFacialValue(HorDiffBEfm + field*BENe*BENfp + face*BENfp, HorDiffvariable + field*Np*K, BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
			/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$, we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToAuxialaryVariable(HorDiffERHSX + field*Np*K*Nface, face, BENe, BENfp, field, HorDiffBEfm, HorDiffBEfp, BEFToE, BEFToF, BEFToN1, Np, K, HorDiffBEFluxM, HorDiffBEFluxS, BEnx, BEMb, BEJs);
			/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial y}$, we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToAuxialaryVariable(HorDiffERHSY + field*Np*K*Nface, face, BENe, BENfp, field, HorDiffBEfm, HorDiffBEfp, BEFToE, BEFToF, BEFToN1, Np, K, HorDiffBEFluxM, HorDiffBEFluxS, BEny, BEMb, BEJs);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int k=0; k<K; k++){
        for(int field=0;field<Nfield;field++){
            for(int face=1;face<Nface;face++){
                Add( HorDiffERHSX + field*Np*K*Nface+k*Np, HorDiffERHSX + field*Np*K*Nface + k*Np, HorDiffERHSX + field*Np*K*Nface + face*Np*K + k*Np, Np);
                Add( HorDiffERHSY + field*Np*K*Nface+k*Np, HorDiffERHSY + field*Np*K*Nface + k*Np, HorDiffERHSY + field*Np*K*Nface + face*Np*K + k*Np, Np);
            }
        }
    }

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nfield; field++){

			MultiEdgeContributionByLiftOperator(HorDiffERHSX + field*Np*K*Nface + k*Np, HorDiffTempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J+k*Np, Np);

			MultiEdgeContributionByLiftOperator(HorDiffERHSY + field*Np*K*Nface + k*Np, HorDiffTempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
	/*Next, sum contribution from volume integral, inner edge contribution and boundary edge contribution into HorDiffAVx and HorDiffAVy*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			Minus(HorDiffAVx + field*Np*K + k*Np, \
				HorDiffVx + field*Np*K + k*Np, HorDiffERHSX + field*Np*K*Nface + k*Np, Np);
			Minus(HorDiffAVy + field*Np*K + k*Np, \
				HorDiffVy + field*Np*K + k*Np, HorDiffERHSY + field*Np*K*Nface + k*Np, Np);
		}
	}
	/*Finally, multiply each component of HorDiffAVx and HorDiffAVy by its diffusion coefficient to get the final auxialary variable,
	this is conducted according to $(Q,v)=(\nv\hat Q,v)$,where $\hat Q = \frac{\partial u(v,\theta)}{\partial (x,y)}$,
	we note that in this projection precedure, aliasing error is introduced.*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			DotProduct(HorDiffAVx + field*Np*K + k*Np, HorDiffAVx + field*Np*K + k*Np, HorDiffnv + k*Np, Np);
			DotProduct(HorDiffAVy + field*Np*K + k*Np, HorDiffAVy + field*Np*K + k*Np, HorDiffnv + k*Np, Np);
		}
		for (int field = 2; field < Nfield; field++){
			DotDivideByConstant(HorDiffAVx + field*Np*K + k*Np, HorDiffAVx + field*Np*K + k*Np, Prantl, Np);
			DotDivideByConstant(HorDiffAVy + field*Np*K + k*Np, HorDiffAVy + field*Np*K + k*Np, Prantl, Np);
		}
	}

	/*Calculate the contribution to the right hand side due to the auxialary variable HorDiffAVx and HorDiffAVy with IPDG.*/
	/*Calculate the penalty parameter $\tau$ first, this parameter is calculated as $\tau=\frac{(N+1)(N+d)}{d}\frac{n_0}{2}\frac{A}{V}\nv$*/
	/*Inner edge first*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < IENe; face++){
		FetchInnerEdgeFacialValue(HorDiffIEnvfm + face*IENfp, HorDiffIEnvfp + face*IENfp, HorDiffnv, \
			IEFToE + face * 2, IEFToN1 + face*IENfp, IEFToN2 + face*IENfp, Np, IENfp);
		double localRatio = IELAV[face] / LAV[(int)IEFToE[face * 2]-1];
		double adjacentRatio = IELAV[face] / LAV[(int)IEFToE[face * 2 + 1]-1];
		for (int p = 0; p < IENfp; p++){
			HorDiffInnerEdgeTau[face*IENfp + p] = max(localRatio*(Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * HorDiffIEnvfm[face*IENfp + p], \
				adjacentRatio*(Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * HorDiffIEnvfp[face*IENfp + p]);
		}
	}
	/*Boundary edge next*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < BENe; face++){
		FetchBoundaryEdgeFacialValue(HorDiffBEnvfm + face*BENfp, HorDiffnv, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
		double localRatio = BELAV[face] / LAV[(int)BEFToE[face * 2]-1];
		for (int p = 0; p < BENfp; p++){
			HorDiffBoundaryEdgeTau[face*BENfp + p] = localRatio * (Order + 1)*(Order + MeshType) / MeshType*Nface / 2 * HorDiffBEnvfm[face*IENfp + p];
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
			GetVolumnIntegral2d(HorDiffVx + field*Np*K + k*Np, HorDiffTempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, HorDiffAVx + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
			/*$\bold{r_y}\cdot (Dr*Q_y)+\bold{s_y}\cdot (Ds*Q_y)$*/
			GetVolumnIntegral2d(HorDiffVy + field*Np*K + k*Np, HorDiffTempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, HorDiffAVy + field*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
			if (type == Two){
				/*The water depth field is excluded from this part*/
				Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					HorDiffVx + field*Np*K + k*Np, Np);
				Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					HorDiffVy + field*Np*K + k*Np, Np);
			}
			else if (type == Three){
				Add(OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					HorDiffVx + field*Np*K + k*Np, Np);
				Add(OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					HorDiffVy + field*Np*K + k*Np, Np);
			}
		}
    }
/************************************************************************************************************************************/
		/*If Mellor and Blumberg's method in "Modeling vertical and horizontal diffusivities with the sigma coordinate system" is
		adopted, the next part should be ignored*/
/************************************************************************************************************************************/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		int field = 0;
		/*$\bold{r_x}\cdot (Dr*Q_x)+\bold{s_x}\cdot (Ds*Q_x)$, with $*Q_x=\nv H\frac{\partial u}{\partial x}*$*/
		GetVolumnIntegral2d(HorDiffVx + field*Np*K + k*Np, HorDiffTempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, HorDiffAVx + field*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
		/*$\bold{r_y}\cdot (Dr*Q_x)+\bold{s_y}\cdot (Ds*Q_x)$, with $*Q_x=\nv H\frac{\partial v}{\partial x}*$*/
		GetVolumnIntegral2d(HorDiffVy + field*Np*K + k*Np, HorDiffTempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, HorDiffAVx + (field + 1)*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		field = 1;
		/*$\bold{r_x}\cdot (Dr*Q_y)+\bold{s_x}\cdot (Ds*Q_y)$, with $*Q_y=\nv H\frac{\partial u}{\partial y}*$*/
		GetVolumnIntegral2d(HorDiffVx + field*Np*K + k*Np, HorDiffTempVx + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, HorDiffAVy + (field - 1)*Np*K + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
		/*$\bold{r_y}\cdot (Dr*Q_y)+\bold{s_y}\cdot (Ds*Q_y)$, with $*Q_y=\nv H\frac{\partial v}{\partial y}*$*/
		GetVolumnIntegral2d(HorDiffVy + field*Np*K + k*Np, HorDiffTempVy + field*Np*K + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, HorDiffAVy + field*Np*K + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		for (int field = 0; field < 2; field++){
			if (type == Two){
				/*The water depth field is excluded from this part*/
				Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					HorDiffVx + field*Np*K + k*Np, Np);
				Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
					HorDiffVy + field*Np*K + k*Np, Np);
			}
			else if (type == Three){
				Add(OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					HorDiffVx + field*Np*K + k*Np, Np);
				Add(OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					OutputRHS + ((int)varIndex[field] - 1)*Np*K + k*Np, \
					HorDiffVy + field*Np*K + k*Np, Np);
			}
		}
	}
/***********************************************************************************************************************************************/
	/*Reset all the data contained in HorDiffERHSX and HorDiffERHSY to zero, becaused these space contains the data left when calculating the auxialary variable*/
	memset(HorDiffERHSX, 0, Np*K*Nfield*Nface*sizeof(double));
	memset(HorDiffERHSY, 0, Np*K*Nfield*Nface*sizeof(double));
	/*Surface integral of the second order operator, this part is calculated seperately for field with index less than 2 and field with index larger than 2 since prantl number should be considered*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int face = 0; face < IENe; face++){
        for (int field = 0; field < 2; field++){
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
			GetIEContributionToRHS(HorDiffERHSX + field*Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
				face, HorDiffLocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
				IEnx, HorDiffIEfm + field*IENe*IENfp, HorDiffIEfp + field*IENe*IENfp, IEnx, HorDiffInnerEdgeTau, 1.0, \
				HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVx + field*Np*K, HorDiffIEFluxM, \
				HorDiffIEFluxP, IEJs, IEMb);
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial y}$*/
			GetIEContributionToRHS(HorDiffERHSY + field*Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
				face, HorDiffLocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
				IEny, HorDiffIEfm + field*IENe*IENfp, HorDiffIEfp + field*IENe*IENfp, IEny, HorDiffInnerEdgeTau, 1.0, \
				HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVy + field*Np*K, HorDiffIEFluxM, \
				HorDiffIEFluxP, IEJs, IEMb);

		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int face = 0; face < IENe; face++){
        for (int field = 2; field < Nfield; field++){
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$*/
			GetIEContributionToRHS(HorDiffERHSX + field*Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
				face, HorDiffLocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
				IEnx, HorDiffIEfm + field*IENe*IENfp, HorDiffIEfp + field*IENe*IENfp, IEnx, HorDiffInnerEdgeTau, Prantl, \
				HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVx + field*Np*K, HorDiffIEFluxM, \
				HorDiffIEFluxP, IEJs, IEMb);
			/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$*/
			GetIEContributionToRHS(HorDiffERHSY + field*Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
				face, HorDiffLocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
				IEny, HorDiffIEfm + field*IENe*IENfp, HorDiffIEfp + field*IENe*IENfp, IEny, HorDiffInnerEdgeTau, Prantl, \
				HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVy + field*Np*K, HorDiffIEFluxM, \
				HorDiffIEFluxP, IEJs, IEMb);

		}
	}
	/*Boundary edge facial integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int face = 0; face < BENe; face++){
        for (int field = 0; field < 2; field++){
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
			GetBEContributionToRHS(HorDiffERHSX + field*Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermX + field*Np*K, \
				BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEnx, HorDiffBEfm + field*BENe*BENfp, \
				HorDiffBEfp + field*BENe*BENfp, BEnx, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVx + field*Np*K, \
				HorDiffBEFluxM, BEJs, BEMb);
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial (u,v)}{\partial y}$*/
			GetBEContributionToRHS(HorDiffERHSY + field*Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermY + field*Np*K, \
				BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEny, HorDiffBEfm + field*BENe*BENfp, \
				HorDiffBEfp + field*BENe*BENfp, BEny, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVy + field*Np*K, \
				HorDiffBEFluxM, BEJs, BEMb);
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int face = 0; face < BENe; face++){
        for (int field = 2; field < Nfield; field++){
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$,
			we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToRHS(HorDiffERHSX + field*Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermX + field*Np*K, \
				BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEnx, HorDiffBEfm + field*BENe*BENfp, \
				HorDiffBEfp + field*BENe*BENfp, BEnx, HorDiffBoundaryEdgeTau, Prantl, HorDiffAVBEfm, HorDiffAVx + field*Np*K, \
				HorDiffBEFluxM, BEJs, BEMb);
			/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$,
			we note that the numerical flux at the boundary is set directly to the outer value*/
			GetBEContributionToRHS(HorDiffERHSY + field*Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermY + field*Np*K, \
				BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEny, HorDiffBEfm + field*BENe*BENfp, \
				HorDiffBEfp + field*BENe*BENfp, BEny, HorDiffBoundaryEdgeTau, Prantl, HorDiffAVBEfm, HorDiffAVy + field*Np*K, \
				HorDiffBEFluxM, BEJs, BEMb);
		}
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int k=0; k<K; k++){
        for(int field=0;field<Nfield;field++){
            for(int face=1;face<Nface;face++){
                Add( HorDiffERHSX + field*Np*K*Nface+k*Np, HorDiffERHSX + field*Np*K*Nface + k*Np, HorDiffERHSX + field*Np*K*Nface + face*Np*K + k*Np, Np);
                Add( HorDiffERHSY + field*Np*K*Nface+k*Np, HorDiffERHSY + field*Np*K*Nface + k*Np, HorDiffERHSY + field*Np*K*Nface + face*Np*K + k*Np, Np);
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nfield; field++){

			MultiEdgeContributionByLiftOperator(HorDiffERHSX + field*Np*K*Nface + k*Np, HorDiffTempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);

			MultiEdgeContributionByLiftOperator(HorDiffERHSY + field*Np*K*Nface + k*Np, HorDiffTempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < Nfield; field++){
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, HorDiffERHSX + field*Np*K*Nface + k*Np, Np);
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, HorDiffERHSY + field*Np*K*Nface + k*Np, Np);
		}
	}
	/****************************************************************************************************************************/
	/*If Mellor and Blumberg's method in "Modeling vertical and horizontal diffusivities with the sigma coordinate system" is
	adopted, the following part about the facial integral should be ignored*/
	memset(HorDiffERHSX, 0, Np * K * 2 * Nface * sizeof(double));
	memset(HorDiffERHSY, 0, Np * K * 2 * Nface * sizeof(double));
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif    
	for (int face = 0; face < IENe; face++){
		int field = 0;
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial u}{\partial x}$*/
		GetIEContributionToRHS(HorDiffERHSX + field*Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
			face, HorDiffLocalPrimitiveDiffTermX + field*Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
			IEnx, HorDiffIEfm + field*IENe*IENfp, HorDiffIEfp + field*IENe*IENfp, IEnx, HorDiffInnerEdgeTau, 1.0, \
			HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVx + field*Np*K, HorDiffIEFluxM, \
			HorDiffIEFluxP, IEJs, IEMb);
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial y}$, here $Q_x=\nv H\frac{\partial v}{\partial x}$*/
		GetIEContributionToRHS(HorDiffERHSY + field*Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
			face, HorDiffLocalPrimitiveDiffTermX + (field + 1)*Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
			IEny, HorDiffIEfm + (field + 1)*IENe*IENfp, HorDiffIEfp + (field + 1)*IENe*IENfp, IEnx, HorDiffInnerEdgeTau, 1.0, \
			HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVx + (field + 1)*Np*K, HorDiffIEFluxM, \
			HorDiffIEFluxP, IEJs, IEMb);

		field = 1;
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial x}$, here $Q_y=\nv H\frac{\partial u}{\partial y}$*/
		GetIEContributionToRHS(HorDiffERHSX + field*Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
			face, HorDiffLocalPrimitiveDiffTermY + (field - 1)*Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
			IEnx, HorDiffIEfm + (field - 1)*IENe*IENfp, HorDiffIEfp + (field - 1)*IENe*IENfp, IEny, HorDiffInnerEdgeTau, 1.0, \
			HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVy + (field - 1)*Np*K, HorDiffIEFluxM, \
			HorDiffIEFluxP, IEJs, IEMb);
		/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial v}{\partial y}$*/
		GetIEContributionToRHS(HorDiffERHSY + field*Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
			face, HorDiffLocalPrimitiveDiffTermY + field*Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
			IEny, HorDiffIEfm + field*IENe*IENfp, HorDiffIEfp + field*IENe*IENfp, IEny, HorDiffInnerEdgeTau, 1.0, \
			HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVy + field*Np*K, HorDiffIEFluxM, \
			HorDiffIEFluxP, IEJs, IEMb);

	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif    
	for (int face = 0; face < BENe; face++){
		int field = 0;
		/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial u}{\partial x}$*/
		GetBEContributionToRHS(HorDiffERHSX + field*Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermX + field*Np*K, \
			BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEnx, HorDiffBEfm + field*BENe*BENfp, \
			HorDiffBEfp + field*BENe*BENfp, BEnx, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVx + field*Np*K, \
			HorDiffBEFluxM, BEJs, BEMb);
		/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial y}$, here $Q_x=\nv H\frac{\partial v}{\partial x}$*/
		GetBEContributionToRHS(HorDiffERHSY + field*Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermX + (field + 1)*Np*K, \
			BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEny, HorDiffBEfm + (field + 1)*BENe*BENfp, \
			HorDiffBEfp + (field + 1)*BENe*BENfp, BEnx, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVx + (field + 1)*Np*K, \
			HorDiffBEFluxM, BEJs, BEMb);

		field = 1;
		/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial x}$, here $Q_y=\nv H\frac{\partial u}{\partial y}$*/
		GetBEContributionToRHS(HorDiffERHSX + field*Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermY + (field - 1)*Np*K, \
			BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEnx, HorDiffBEfm + (field - 1)*BENe*BENfp, \
			HorDiffBEfp + (field - 1)*BENe*BENfp, BEny, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVy + (field - 1)*Np*K, \
			HorDiffBEFluxM, BEJs, BEMb);
		/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial v}{\partial y}$*/
		GetBEContributionToRHS(HorDiffERHSY + field*Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermY + field*Np*K, \
			BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEny, HorDiffBEfm + field*BENe*BENfp, \
			HorDiffBEfp + field*BENe*BENfp, BEny, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVy + field*Np*K, \
			HorDiffBEFluxM, BEJs, BEMb);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k<K; k++){
		for (int field = 0; field<2; field++){
			for (int face = 1; face<Nface; face++){
				Add(HorDiffERHSX + field*Np*K*Nface + k*Np, HorDiffERHSX + field*Np*K*Nface + k*Np, HorDiffERHSX + field*Np*K*Nface + face*Np*K + k*Np, Np);
				Add(HorDiffERHSY + field*Np*K*Nface + k*Np, HorDiffERHSY + field*Np*K*Nface + k*Np, HorDiffERHSY + field*Np*K*Nface + face*Np*K + k*Np, Np);
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < 2; field++){

			MultiEdgeContributionByLiftOperator(HorDiffERHSX + field*Np*K*Nface + k*Np, HorDiffTempFacialIntegralX + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);

			MultiEdgeContributionByLiftOperator(HorDiffERHSY + field*Np*K*Nface + k*Np, HorDiffTempFacialIntegralY + field*Np*K + k*Np, &np, &oneI, &np, \
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
				OutputRHS + field*Np*K + k*Np, HorDiffERHSX + field*Np*K*Nface + k*Np, Np);
			Minus(OutputRHS + field*Np*K + k*Np, \
				OutputRHS + field*Np*K + k*Np, HorDiffERHSY + field*Np*K*Nface + k*Np, Np);
		}
	}
}