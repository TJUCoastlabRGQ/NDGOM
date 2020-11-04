#include "..\..\..\..\..\NdgMath\NdgMath.h"
#include "..\..\..\..\..\NdgMath\NdgSWE.h"
#include "..\..\..\..\..\NdgMath\NdgSWE3D.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	/*Parameter reading part*/
	/*For mesh object*/
	mxArray *Temprx = mxGetField(prhs[0], 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(prhs[0], 0, "sx");
	double *sx = mxGetPr(Tempsx);
//	mxArray *Temptx = mxGetField(prhs[0], 0, "tx");
//	double *tx = mxGetPr(Temptx);
	mxArray *Tempry = mxGetField(prhs[0], 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *Tempsy = mxGetField(prhs[0], 0, "sy");
	double *sy = mxGetPr(Tempsy);
//	mxArray *Tempty = mxGetField(prhs[0], 0, "ty");
//	double *ty = mxGetPr(Tempty);
//	mxArray *Temprz = mxGetField(prhs[0], 0, "rz");
//	double *rz = mxGetPr(Temprz);
//	mxArray *Tempsz = mxGetField(prhs[0], 0, "sz");
//	double *sz = mxGetPr(Tempsz);
	mxArray *Temptz = mxGetField(prhs[0], 0, "tz");
	double *tz = mxGetPr(Temptz);
	mxArray *TempJ = mxGetField(prhs[0], 0, "J");
	double *J = mxGetPr(TempJ);
	int Np = mxGetM(TempJ);
	int K = mxGetN(TempJ);
	/*For cell object*/
	mxArray *TempDr = mxGetField(prhs[1], 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(prhs[1], 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(prhs[1], 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempinvM = mxGetField(prhs[1], 0, "invM");
	double *invM = mxGetPr(TempinvM);
	/*For inner edge object*/
	mxArray *TempIENe = mxGetField(prhs[2], 0, "Ne");
	int IENe = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp = mxGetField(prhs[2], 0, "Nfp");
	int IENfp = (int)mxGetScalar(TempIENfp);
	mxArray *TempIEMb = mxGetField(prhs[2], 0, "M");
	double *IEMb = mxGetPr(TempIEMb);
	mxArray *TempIEJs = mxGetField(prhs[2], 0, "Js");
	double *IEJs = mxGetPr(TempIEJs);
	mxArray *TempIEnx = mxGetField(prhs[2], 0, "nx");
	double *IEnx = mxGetPr(TempIEnx);
	mxArray *TempIEny = mxGetField(prhs[2], 0, "ny");
	double *IEny = mxGetPr(TempIEny);
//	mxArray *TempIEnz = mxGetField(prhs[2], 0, "nz");
//	double *IEnz = mxGetPr(TempIEnz);
	mxArray *TempIEFToE = mxGetField(prhs[2], 0, "FToE");
	double *IEFToE = mxGetPr(TempIEFToE);
	mxArray *TempIEFToN1 = mxGetField(prhs[2], 0, "FToN1");
	double *IEFToN1 = mxGetPr(TempIEFToN1);
	mxArray *TempIEFToN2 = mxGetField(prhs[2], 0, "FToN2");
	double *IEFToN2 = mxGetPr(TempIEFToN2);
	/*For boundary edge object*/
	mxArray *TempBENe = mxGetField(prhs[3], 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBENfp = mxGetField(prhs[3], 0, "Nfp");
	int BENfp = mxGetScalar(TempBENfp);
	mxArray *TempBEMb = mxGetField(prhs[3], 0, "M");
	double *BEMb = mxGetPr(TempBEMb);
	mxArray *TempBEJs = mxGetField(prhs[3], 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	mxArray *TempBEnx = mxGetField(prhs[3], 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(prhs[3], 0, "ny");
	double *BEny = mxGetPr(TempBEny);
//	mxArray *TempBEnz = mxGetField(prhs[3], 0, "nz");
//	double *BEnz = mxGetPr(TempBEnz);
//	mxArray *TempBELAV = mxGetField(prhs[3], 0, "LAV");
//	double *BELAV = mxGetPr(TempBELAV);
	mxArray *TempBEFToE = mxGetField(prhs[3], 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToN1 = mxGetField(prhs[3], 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);
	/*For bottom edge object*/
	mxArray *TempBotNe = mxGetField(prhs[4], 0, "Ne");
	int BotENe = (int)mxGetScalar(TempBotNe);
	mxArray *TempBotNfp = mxGetField(prhs[4], 0, "Nfp");
	int BotENfp = (int)mxGetScalar(TempBotNfp);
	mxArray *TempBotMb = mxGetField(prhs[4], 0, "M");
	double *BotEMb = mxGetPr(TempBotMb);
	mxArray *TempBotJs = mxGetField(prhs[4], 0, "Js");
	double *BotEJs = mxGetPr(TempBotJs);
//	mxArray *TempBotEnx = mxGetField(prhs[4], 0, "nx");
//	double *BotEnx = mxGetPr(TempBotEnx);
//	mxArray *TempBotEny = mxGetField(prhs[4], 0, "ny");
//	double *BotEny = mxGetPr(TempBotEny);
	mxArray *TempBotEnz = mxGetField(prhs[4], 0, "nz");
	double *BotEnz = mxGetPr(TempBotEnz);
	mxArray *TempBotEFToE = mxGetField(prhs[4], 0, "FToE");
	double *BotEFToE = mxGetPr(TempBotEFToE);
	mxArray *TempBotEFToN1 = mxGetField(prhs[4], 0, "FToN1");
	double *BotEFToN1 = mxGetPr(TempBotEFToN1);
	mxArray *TempBotEFToN2 = mxGetField(prhs[4], 0, "FToN2");
	double *BotEFToN2 = mxGetPr(TempBotEFToN2);
	/*For bottom boundary edge object*/
	mxArray *TempBotBENe = mxGetField(prhs[5], 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBENfp = mxGetField(prhs[5], 0, "Nfp");
	int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	mxArray *TempBotBEMb = mxGetField(prhs[5], 0, "M");
	double *BotBEMb = mxGetPr(TempBotBEMb);
	mxArray *TempBotBEJs = mxGetField(prhs[5], 0, "Js");
	double *BotBEJs = mxGetPr(TempBotBEJs);
//	mxArray *TempBotBEnx = mxGetField(prhs[5], 0, "nx");
//	double *BotBEnx = mxGetPr(TempBotBEnx);
//	mxArray *TempBotBEny = mxGetField(prhs[5], 0, "ny");
//	double *BotBEny = mxGetPr(TempBotBEny);
	mxArray *TempBotBEnz = mxGetField(prhs[5], 0, "nz");
	double *BotBEnz = mxGetPr(TempBotBEnz);
	mxArray *TempBotBEFToE = mxGetField(prhs[5], 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToN1 = mxGetField(prhs[5], 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);
    /*For surface boundary edge object*/
	mxArray *TempSurfBENe = mxGetField(prhs[6], 0, "Ne");
	int SurfBENe = (int)mxGetScalar(TempSurfBENe);
	mxArray *TempSurfBENfp = mxGetField(prhs[6], 0, "Nfp");
	int SurfBENfp = (int)mxGetScalar(TempSurfBENfp);
	mxArray *TempSurfBEMb = mxGetField(prhs[6], 0, "M");
	double *SurfBEMb = mxGetPr(TempSurfBEMb);
	mxArray *TempSurfBEJs = mxGetField(prhs[6], 0, "Js");
	double *SurfBEJs = mxGetPr(TempSurfBEJs);
//	mxArray *TempSurfBEnx = mxGetField(prhs[6], 0, "nx");
//	double *SurfBEnx = mxGetPr(TempSurfBEnx);
//	mxArray *TempSurfBEny = mxGetField(prhs[6], 0, "ny");
//	double *SurfBEny = mxGetPr(TempSurfBEny);
	mxArray *TempSurfBEnz = mxGetField(prhs[6], 0, "nz");
	double *SurfBEnz = mxGetPr(TempSurfBEnz);
	mxArray *TempSurfBEFToE = mxGetField(prhs[6], 0, "FToE");
	double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
	mxArray *TempSurfBEFToN1 = mxGetField(prhs[6], 0, "FToN1");
	double *SurfBEFToN1 = mxGetPr(TempSurfBEFToN1);

	double *varFieldIndex = mxGetPr(prhs[7]);
	int Nvar = (int)mxGetNumberOfElements(prhs[7]);

	double *fphys = mxGetPr(prhs[8]);
	double *h = fphys + 3 * Np*K, *hu = fphys, *hv = fphys + Np*K, \
		*omega = fphys + 2 * Np*K, *z = fphys + 5 * Np*K;
	double *fext = mxGetPr(prhs[9]);
	double gra = mxGetScalar(prhs[10]);
	double Hcrit = mxGetScalar(prhs[11]);

	signed char *ftype = (signed char *)mxGetData(prhs[12]);


    size_t NdimOut = 3;
	mwSize dimOut[3] = { Np, K, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputRHS = mxGetPr(plhs[0]);
	                                                                               /************************  Face Integral Part  ****************************/
/*************************************************************************************************************************************/
	/**************************************Inner Edge Part*******************************************************/
	            /********************************************************************/
	/*Allocate memory for fm and fp defined over inner edges. Here, variables correspond to hu, hv, hT, hS, sediment and other passive transport material, vertical velocity omega is not included*/
	double *IEfm = malloc(IENfp*IENe*(Nvar + 1)*sizeof(double));
	double *huM = IEfm, *hvM = IEfm + IENfp*IENe, *hM = IEfm + 2 * IENfp*IENe;
	double *IEfp = malloc(IENfp*IENe*(Nvar + 1)*sizeof(double));
	double *huP = IEfp, *hvP = IEfp + IENfp*IENe, *hP = IEfp + 2 * IENfp*IENe;
	/*Allocate memory for fluxM, fluxP and fluxS, and calculate these flux term*/
	double *IEFluxM = malloc(IENfp*IENe*Nvar*sizeof(double));
	double *IEFluxP = malloc(IENfp*IENe*Nvar*sizeof(double));
	double *IEFluxS = malloc(IENfp*IENe*Nvar*sizeof(double));
	memset(IEFluxS, 0, IENfp*IENe*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < IENe; face++){
		/*Fetch variable fm and fp first*/
		FetchInnerEdgeFacialValue(hM + face*IENfp, hP + face*IENfp, h, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		FetchInnerEdgeFacialValue(huM + face*IENfp, huP + face*IENfp, hu, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		FetchInnerEdgeFacialValue(hvM + face*IENfp, hvP + face*IENfp, hv, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
		here 1 stands for the space occupied by water depth h.*/
		for (int field = 2; field < Nvar; field++){
			FetchInnerEdgeFacialValue(IEfm + (field + 1)*IENe*IENfp + face*IENfp, \
				IEfp + (field + 1)*IENe*IENfp + face*IENfp, fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				IEFToE + 2 * face, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		}
		EvaluateVerticalFaceSurfFlux(IEFluxM + face*IENfp, IEfm + face*IENfp, IEnx + face*IENfp, IEny + face*IENfp, &gra, Hcrit, IENfp, Nvar, IENe);
		EvaluateVerticalFaceSurfFlux(IEFluxP + face*IENfp, IEfp + face*IENfp, IEnx + face*IENfp, IEny + face*IENfp, &gra, Hcrit, IENfp, Nvar, IENe);
		EvaluateVerticalFaceNumFlux(IEFluxS + face*IENfp, IEfm + face*IENfp, IEfp + face*IENfp, \
			IEnx + face*IENfp, IEny + face*IENfp, &gra, Hcrit, IENfp, Nvar, IENe);
	}

/*Allocate memory for contribution to RHS due to inner edge facial integral, and
 calculate contribution to RHS due to inner edge facial integral in strong form manner*/
	double *IERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(IERHS, 0, Np*K*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < IENe; face++){
			StrongFormInnerEdgeRHS(face, IEFToE, Np, IENfp, IEFToN1, IEFToN2, IEFluxM + field*IENe*IENfp,\
				IEFluxP + field*IENe*IENfp, IEFluxS + field*IENe*IENfp, IEJs, IEMb, IERHS + field*Np*K);
		}
	}
/*Multiply the contribution to RHS due to inner edge facial integral by the inverse mass matrix*/
	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;
	double *TempFacialIntegral = malloc(Np*K*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nvar; field++){
			MultiEdgeContributionByLiftOperator(IERHS + field*Np*K + k*Np, TempFacialIntegral + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
	free(IEfm);
	free(IEfp);
	free(IEFluxM);
	free(IEFluxP);
	free(IEFluxS);
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	            /**************************************Boundary Edge Part*******************************************************/
	                           /********************************************************************/
	/*Allocate memory for fm and fp defined over boundary edges. Here, variables correspond to hu, hv, h, hT, hS, 
	sediment and other passive transport material, vertical velocity omega is not included for boundary edges*/
	double *BEfm = malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
	huM = BEfm, hvM = BEfm + BENfp*BENe, hM = BEfm + 2 * BENfp*BENe;
	double *BEfp = malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
	huP = BEfp, hvP = BEfp + BENfp*BENe, hP = BEfp + 2 * BENfp*BENe;
	double *zM = malloc(BENfp*BENe*sizeof(double));
	double *zP = malloc(BENfp*BENe*sizeof(double));
	double *BEFluxM = malloc(BENfp*BENe*Nvar*sizeof(double));
	double *BEFluxS = malloc(BENfp*BENe*Nvar*sizeof(double));
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
			FetchBoundaryEdgeFacialValue(BEfm + (field + 1)*BENe*BENfp + face*BENfp, \
				 fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
		}
		ImposeBoundaryCondition(&gra, type, BEnx + face*BENfp, BEny + face*BENfp, BEfm + face*BENfp, BEfp + face*BENfp, \
			zM + face*BENfp, zP + face*BENfp, fext + face*BENfp, BENfp, Nvar + 1, BENe);
		EvaluateHydroStaticReconstructValue(Hcrit, BEfm + face*BENfp, BEfp + face*BENfp, zM + face*BENfp, zP + face*BENfp, BENfp, Nvar + 1, BENe);
		EvaluateVerticalFaceSurfFlux(BEFluxM + face*BENfp, BEfm + face*BENfp, BEnx + face*BENfp, BEny + face*BENfp, &gra, Hcrit, BENfp, Nvar, BENe);
		EvaluateVerticalFaceNumFlux(BEFluxS + face*BENfp, BEfm + face*BENfp, BEfp + face*BENfp, \
			BEnx + face*BENfp, BEny + face*BENfp, &gra, Hcrit, BENfp, Nvar, BENe);
	}
	/*Allocate memory for contribution to RHS due to boundary edge facial integral, and calculate 
	contribution to RHS due to boundary edge facial integral in strong form manner*/
	double *BERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(BERHS, 0, Np*K*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < BENe; face++){
			StrongFormBoundaryEdgeRHS(face, BEFToE, Np, BENfp, BEFToN1, BEFluxM + field*BENe*BENfp, BEFluxS + field*BENe*BENfp, BEJs, BEMb, BERHS + field*Np*K);
		}
	}
	/*Multiply the contribution to RHS due to boundary edge facial integral by the inverse mass matrix*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nvar; field++){
			MultiEdgeContributionByLiftOperator(BERHS + field*Np*K + k*Np, TempFacialIntegral + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
	free(BEfm);
	free(BEfp);
	free(BEFluxM);
	free(BEFluxS);
	free(zM);
	free(zP);
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	           /**************************************Bottom Edge Part*******************************************************/
	                              /********************************************************************/
	double *BotEfm = malloc(BotENfp*BotENe*(Nvar + 2)*sizeof(double));
	double *omegaM = BotEfm + 2 * BotENfp*BotENe;
	huM = BotEfm, hvM = BotEfm + BotENfp*BotENe, hM = BotEfm + 3 * BotENfp*BotENe;
	double *BotEfp = malloc(BotENfp*BotENe*(Nvar + 2)*sizeof(double));
	double *omegaP = BotEfp + 2 * BotENfp*BotENe;
	huP = BotEfp, hvP = BotEfp + BotENfp*BotENe, hP = BotEfp + 3 * BotENfp*BotENe;
	/*Allocate memory for fluxM, fluxP and fluxS, and calculate these flux term*/
	double *BotEFluxM = malloc(BotENfp*BotENe*Nvar*sizeof(double));
	double *BotEFluxP = malloc(BotENfp*BotENe*Nvar*sizeof(double));
	double *BotEFluxS = malloc(BotENfp*BotENe*Nvar*sizeof(double));
	memset(BotEFluxS, 0, BotENfp*BotENe*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < BotENe; face++){
		/*Fetch variable fm and fp first*/
		FetchInnerEdgeFacialValue(hM + face*BotENfp, hP + face*BotENfp, h, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		FetchInnerEdgeFacialValue(huM + face*BotENfp, huP + face*BotENfp, hu, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		FetchInnerEdgeFacialValue(hvM + face*BotENfp, hvP + face*BotENfp, hv, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		FetchInnerEdgeFacialValue(omegaM + face*BotENfp, omegaP + face*BotENfp, omega, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included*/
		for (int field = 2; field < Nvar; field++){
			FetchInnerEdgeFacialValue(BotEfm + (field+2)*IENe*IENfp + face*IENfp, \
				IEfp + (field + 2)*IENe*IENfp + face*IENfp, fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				IEFToE + 2 * face, IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		}
		EvaluateHorizontalFaceSurfFlux(BotEFluxM + face*BotENfp, BotEfm + face*BotENfp, BotEnz + face*BotENfp, Hcrit, BotENfp, Nvar, BotENe);
		EvaluateHorizontalFaceSurfFlux(BotEFluxP + face*BotENfp, BotEfp + face*BotENfp, BotEnz + face*BotENfp, Hcrit, BotENfp, Nvar, BotENe);
		EvaluateHorizontalFaceNumFlux(BotEFluxS + face*BotENfp, BotEfm + face*BotENfp, BotEfp + face*BotENfp, \
			BotEnz + face*BotENfp, Hcrit, BotENfp, Nvar, BotENe);
	}


	/*Allocate memory for contribution to RHS due to bottom edge facial integral, and
	calculate contribution to RHS due to bottom edge facial integral in strong form manner*/
	double *BotERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(BotERHS, 0, Np*K*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < BotENe; face++){
			StrongFormInnerEdgeRHS(face, BotEFToE, Np, BotENfp, BotEFToN1, BotEFToN2, BotEFluxM + field*BotENe*BotENfp, \
				BotEFluxP + field*BotENe*BotENfp, BotEFluxS + field*BotENe*BotENfp, BotEJs, BotEMb, BotERHS + field*Np*K);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nvar; field++){
			MultiEdgeContributionByLiftOperator(BotERHS + field*Np*K + k*Np, TempFacialIntegral + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}

	free(BotEfm);
	free(BotEfp);
	free(BotEFluxM);
	free(BotEFluxP);
	free(BotEFluxS);
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	             /**************************************Bottom Boundary Edge Part*******************************************************/
	                            /********************************************************************/
	double *BotBEfm = malloc(BotBENfp*BotBENe*(Nvar + 2)*sizeof(double));
	huM = BotBEfm;
	hvM = BotBEfm + BotBENfp*BotBENe;
	omegaM = BotBEfm + 2 * BotBENfp*BotBENe;
	hM = BotBEfm + 3 * BotBENfp*BotBENe;

	/*Allocate memory for fluxM, fluxP and fluxS, and calculate these flux term*/
	double *BotBEFluxM = malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
	double *BotBEFluxS = malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
	memset(BotBEFluxS, 0, BotBENfp*BotBENe*Nvar*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < BotBENe; face++){
		FetchBoundaryEdgeFacialValue(huM + face*BotBENfp, hu, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*BotBENfp, hv, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(omegaM + face*BotBENfp, omega, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(hM + face*BotBENfp, h, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included
		1 stands for the memory occupied by water depth h*/
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(BotBEfm + (field + 1)*BotBENe*BotBENfp + face*BotBENfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				BotBEFToE + 2 * face, BotBEFToN1 + BotBENfp*face, Np, BotBENfp);
		}

		EvaluateHorizontalFaceSurfFlux(BotBEFluxM + face*BotBENfp, BotBEfm + face*BotBENfp, BotBEnz + face*BotBENfp, Hcrit, BotBENfp, Nvar, BotBENe);
	}

	/*Allocate memory for contribution to RHS due to boundary edge facial integral, and calculate
	contribution to RHS due to boundary edge facial integral in strong form manner*/
	double *BotBERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(BotBERHS, 0, Np*K*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < BotBENe; face++){
			StrongFormBoundaryEdgeRHS(face, BotBEFToE, Np, BotBENfp, BotBEFToN1, BotBEFluxM + field*BotBENe*BotBENfp, BotBEFluxS + field*BotBENe*BotBENfp, BotBEJs, BotBEMb, BotBERHS + field*Np*K);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nvar; field++){
			MultiEdgeContributionByLiftOperator(BotBERHS + field*Np*K + k*Np, TempFacialIntegral + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
	free(BotBEfm);
	free(BotBEFluxM);
	free(BotBEFluxS);
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	         /**************************************Surface Boundary Edge Part*******************************************************/
	                               /********************************************************************/
	double *SurfBEfm = malloc(SurfBENfp*SurfBENe*(Nvar + 2)*sizeof(double));
	huM = SurfBEfm;
	hvM = SurfBEfm + SurfBENfp*SurfBENe;
	omegaM = SurfBEfm + 2 * SurfBENfp*SurfBENe;
	hM = SurfBEfm + 3 * SurfBENfp*SurfBENe;

	/*Allocate memory for fluxM, fluxP and fluxS, and calculate these flux term*/
	double *SurfBEFluxM = malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
	/*WE NOTE THAT, FOR THE CONVERGENCE TEST, THIS FLUX IS NOT TAKEN AS ZERO AND SHOULD BE TAKEN FORM THE INPUT*/
	double *SurfBEFluxS = malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
	memset(SurfBEFluxS, 0, SurfBENfp*SurfBENe*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < SurfBENe; face++){
		FetchBoundaryEdgeFacialValue(huM + face*SurfBENfp, hu, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*SurfBENfp, hv, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(omegaM + face*SurfBENfp, omega, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(hM + face*SurfBENfp, h, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);

		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included
		1 stands for the memory occupied by water depth h*/
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(SurfBEfm + (field + 1)*SurfBENe*SurfBENfp + face*SurfBENfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				SurfBEFToE + 2 * face, SurfBEFToN1 + SurfBENfp*face, Np, SurfBENfp);
		}

		EvaluateHorizontalFaceSurfFlux(SurfBEFluxM + face*SurfBENfp, SurfBEfm + face*SurfBENfp, SurfBEnz + face*SurfBENfp, Hcrit, SurfBENfp, Nvar, SurfBENe);
	}

	/*Allocate memory for contribution to RHS due to boundary edge facial integral, and calculate
	contribution to RHS due to boundary edge facial integral in strong form manner*/
	double *SurfBERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(SurfBERHS, 0, Np*K*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < SurfBENe; face++){
			StrongFormBoundaryEdgeRHS(face, SurfBEFToE, Np, SurfBENfp, SurfBEFToN1, SurfBEFluxM + field*SurfBENe*SurfBENfp, SurfBEFluxS + field*SurfBENe*SurfBENfp, SurfBEJs, SurfBEMb, SurfBERHS + field*Np*K);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nvar; field++){
			MultiEdgeContributionByLiftOperator(SurfBERHS + field*Np*K + k*Np, TempFacialIntegral + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}
	free(SurfBEfm);
	free(SurfBEFluxM);
	free(SurfBEFluxS);
	free(TempFacialIntegral);
	/**********************
	/***************************************************************************************************************/


	                                                                                    /************************  Volume Integral Part  ****************************/
	/*Allocate memory for E, G and H, and calculate these volume flux term*/
	double *E = malloc(Np*K*Nvar*sizeof(double));
	double *G = malloc(Np*K*Nvar*sizeof(double));
	double *H = malloc(Np*K*Nvar*sizeof(double));
	double *TempVolumeIntegral = malloc(Np*K*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		EvaluatePrebalanceVolumeTerm(E + k*Np, G + k*Np, H + k*Np, fphys + k*Np, \
			varFieldIndex, Nvar, &gra, Np, K, Hcrit);

		GetVolumnIntegral3d(OutputRHS + k*Np, TempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, Dt, E + k*Np, G + k*Np, H + k*Np, &np, &np, &zero, \
			&np, rx + k*Np, sx + k*Np, ry + k*Np, sy + k*Np, tz + k*Np, Nvar, Np, K);
	}
	free(E);
	free(G);
	free(H);
	free(TempVolumeIntegral);

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		for (int n = 0; n < Nvar; n++){
			Minus(OutputRHS + n*Np*K + k*Np, \
				IERHS + n*Np*K + k*Np, OutputRHS + n*Np*K + k*Np, Np);
			Add(OutputRHS + n*Np*K + k*Np, \
				BERHS + n*Np*K + k*Np, OutputRHS + n*Np*K + k*Np, Np);
			Add(OutputRHS + n*Np*K + k*Np, \
				BotERHS + n*Np*K + k*Np, OutputRHS + n*Np*K + k*Np, Np);
			Add(OutputRHS + n*Np*K + k*Np, \
				BotBERHS + n*Np*K + k*Np, OutputRHS + n*Np*K + k*Np, Np);
			Add(OutputRHS + n*Np*K + k*Np, \
				SurfBERHS + n*Np*K + k*Np, OutputRHS + n*Np*K + k*Np, Np);
		}
	}

	free(IERHS);
	free(BERHS);
	free(BotERHS);
	free(BotBERHS);
	free(SurfBERHS);

}