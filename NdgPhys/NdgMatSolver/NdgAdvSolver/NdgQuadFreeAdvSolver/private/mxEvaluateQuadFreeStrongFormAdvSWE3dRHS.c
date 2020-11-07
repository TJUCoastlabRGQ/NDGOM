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
	mxArray *TempBEFToE = mxGetField(prhs[3], 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToN1 = mxGetField(prhs[3], 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);
	/*For bottom edge object*/
	mxArray *TempBotENe = mxGetField(prhs[4], 0, "Ne");
	int BotENe = (int)mxGetScalar(TempBotENe);
	mxArray *TempBotENfp = mxGetField(prhs[4], 0, "Nfp");
	int BotENfp = (int)mxGetScalar(TempBotENfp);
	mxArray *TempBotEMb = mxGetField(prhs[4], 0, "M");
	double *BotEMb = mxGetPr(TempBotEMb);
	mxArray *TempBotEJs = mxGetField(prhs[4], 0, "Js");
	double *BotEJs = mxGetPr(TempBotEJs);
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
	double *hu = fphys, *hv = fphys + Np*K, *omega = fphys + 2 * Np*K,  \
		*h = fphys + 3 * Np*K, *z = fphys + 5 * Np*K;
	double *fext = mxGetPr(prhs[9]);
	double gra = mxGetScalar(prhs[10]);
	double Hcrit = mxGetScalar(prhs[11]);

	signed char *ftype = (signed char *)mxGetData(prhs[12]);


    size_t NdimOut = 3;
	mwSize dimOut[3] = { Np, K, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputRHS = mxGetPr(plhs[0]);
	plhs[1] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputE = mxGetPr(plhs[1]);
	plhs[2] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputG = mxGetPr(plhs[2]);
	plhs[3] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputH = mxGetPr(plhs[3]);
	plhs[4] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputVolumn = mxGetPr(plhs[4]);

	mwSize IEdimOut[3] = { IENfp, IENe, Nvar };
	plhs[5] = mxCreateNumericArray(NdimOut, IEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputIEFluxM = mxGetPr(plhs[5]);
	plhs[6] = mxCreateNumericArray(NdimOut, IEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputIEFluxP = mxGetPr(plhs[6]);
	plhs[7] = mxCreateNumericArray(NdimOut, IEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputIEFluxS = mxGetPr(plhs[7]);

	mwSize BEdimOut[3] = { BENfp, BENe, Nvar };
	plhs[8] = mxCreateNumericArray(NdimOut, BEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputBEFluxM = mxGetPr(plhs[8]);
	plhs[9] = mxCreateNumericArray(NdimOut, BEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputBEFluxP = mxGetPr(plhs[9]);
	plhs[10] = mxCreateNumericArray(NdimOut, BEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputBEFluxS = mxGetPr(plhs[10]);

	mwSize BotEdimOut[3] = { BotENfp, BotENe, Nvar };
	plhs[11] = mxCreateNumericArray(NdimOut, BotEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputBotEFluxM = mxGetPr(plhs[11]);
	plhs[12] = mxCreateNumericArray(NdimOut, BotEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputBotEFluxP = mxGetPr(plhs[12]);
	plhs[13] = mxCreateNumericArray(NdimOut, BotEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputBotEFluxS = mxGetPr(plhs[13]);

	mwSize BotBEdimOut[3] = { BotBENfp, BotBENe, Nvar };
	plhs[14] = mxCreateNumericArray(NdimOut, BotBEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputBotBEFluxM = mxGetPr(plhs[14]);
	plhs[15] = mxCreateNumericArray(NdimOut, BotBEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputBotBEFluxS = mxGetPr(plhs[15]);

	mwSize SurfBEdimOut[3] = { SurfBENfp, SurfBENe, Nvar };
	plhs[16] = mxCreateNumericArray(NdimOut, SurfBEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputSurfBEFluxM = mxGetPr(plhs[16]);
	plhs[17] = mxCreateNumericArray(NdimOut, SurfBEdimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputSurfBEFluxS = mxGetPr(plhs[17]);

	
	int Ne, Nfp;

	double *FToN1, *FToN2, *FToE, *nx, *ny, *nz;

	double *fm, *fp, *FluxM, *FluxP, *FluxS;

	                                                                               /************************  Face Integral Part  ****************************/
/*************************************************************************************************************************************/
	/**************************************Inner Edge Part*******************************************************/
	            /********************************************************************/
	/*Allocate memory for fm and fp defined over inner edges. Here, variables correspond to hu, hv, hT, hS, sediment and other passive transport material, vertical velocity omega is not included*/
	Ne = IENe, Nfp = IENfp;
	FToN1 = IEFToN1, FToN2 = IEFToN2, FToE = IEFToE, nx = IEnx, ny = IEny;
	fm = malloc(Nfp*Ne*(Nvar + 1)*sizeof(double));
	double *huM = fm, *hvM = fm + Nfp*Ne, *hM = fm + 2 * Nfp*Ne;
	fp = malloc(Nfp*Ne*(Nvar + 1)*sizeof(double));
	double *huP = fp, *hvP = fp + Nfp*Ne, *hP = fp + 2 * Nfp*Ne;
	/*Allocate memory for fluxM, fluxP and fluxS, and calculate these flux term*/
	FluxM = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxM, 0, Nfp*Ne*Nvar*sizeof(double));
	FluxP = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxP, 0, Nfp*Ne*Nvar*sizeof(double));
	FluxS = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxS, 0, Nfp*Ne*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < Ne; face++){
		/*Fetch variable fm and fp first*/
		FetchInnerEdgeFacialValue(hM + face*Nfp, hP + face*Nfp, h, FToE + 2 * face, \
			FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		FetchInnerEdgeFacialValue(huM + face*Nfp, huP + face*Nfp, hu, FToE + 2 * face, \
			FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		FetchInnerEdgeFacialValue(hvM + face*Nfp, hvP + face*Nfp, hv, FToE + 2 * face, \
			FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
		here 1 stands for the space occupied by water depth h.*/
		for (int field = 2; field < Nvar; field++){
			FetchInnerEdgeFacialValue(fm + (field + 1)*Ne*Nfp + face*Nfp, \
				fp + (field + 1)*Ne*Nfp + face*Nfp, fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				FToE + 2 * face, FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		}
		EvaluateVerticalFaceSurfFlux(FluxM + face*Nfp, fm + face*Nfp, nx + face*Nfp, ny + face*Nfp, &gra, Hcrit, Nfp, Nvar, Ne);
		EvaluateVerticalFaceSurfFlux(FluxP + face*Nfp, fp + face*Nfp, nx + face*Nfp, ny + face*Nfp, &gra, Hcrit, Nfp, Nvar, Ne);
		EvaluateVerticalFaceNumFlux(FluxS + face*Nfp, fm + face*Nfp, fp + face*Nfp, \
			nx + face*Nfp, ny + face*Nfp, &gra, Hcrit, Nfp, Nvar, Ne);
	}
	memcpy(OutputIEFluxM, FluxM, Nfp*Ne*Nvar*sizeof(double));
	memcpy(OutputIEFluxP, FluxP, Nfp*Ne*Nvar*sizeof(double));
	memcpy(OutputIEFluxS, FluxS, Nfp*Ne*Nvar*sizeof(double));


/*Allocate memory for contribution to RHS due to inner edge facial integral, and
 calculate contribution to RHS due to inner edge facial integral in strong form manner*/
	double *IERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(IERHS, 0, Np*K*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < Ne; face++){
			StrongFormInnerEdgeRHS(face, FToE, Np, Nfp, FToN1, FToN2, FluxM + field*Ne*Nfp,\
				FluxP + field*Ne*Nfp, FluxS + field*Ne*Nfp, IEJs, IEMb, IERHS + field*Np*K);
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
	free(fm);
	free(fp);
	free(FluxM);
	free(FluxP);
	free(FluxS);
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	            /**************************************Boundary Edge Part*******************************************************/
	                           /********************************************************************/
	/*Allocate memory for fm and fp defined over boundary edges. Here, variables correspond to hu, hv, h, hT, hS, 
	sediment and other passive transport material, vertical velocity omega is not included for boundary edges*/
	Ne = BENe, Nfp = BENfp;
	FToN1 = BEFToN1, FToE = BEFToE, nx = BEnx, ny = BEny;
	fm = malloc(Nfp*Ne*(Nvar + 1)*sizeof(double));
	huM = fm, hvM = fm + Nfp*Ne, hM = fm + 2 * Nfp*Ne;
	fp = malloc(Nfp*Ne*(Nvar + 1)*sizeof(double));
	huP = fp, hvP = fp + Nfp*Ne, hP = fp + 2 * Nfp*Ne;
	double *zM = malloc(Nfp*Ne*sizeof(double));
	double *zP = malloc(Nfp*Ne*sizeof(double));
	FluxM = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxM, 0, Nfp*Ne*Nvar*sizeof(double));
	FluxS = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxS, 0, Nfp*Ne*Nvar*sizeof(double));

	/*Fetch variable fm and fp first, then impose boundary condition and conduct hydrostatic reconstruction.
	Finally, calculate local flux term, adjacent flux term and numerical flux term*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < Ne; face++){
		NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
		FetchBoundaryEdgeFacialValue(huM + face*Nfp, hu, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(hvM + face*Nfp, hv, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(hM + face*Nfp, h, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(zM + face*Nfp, z, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
		here 1 stands for the memory occupied by water depth h*/
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(fm + (field + 1)*Ne*Nfp + face*Nfp, \
				 fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				FToE + 2 * face, FToN1 + Nfp*face, Np, Nfp);
		}
		ImposeBoundaryCondition(&gra, type, nx + face*Nfp, ny + face*Nfp, fm + face*Nfp, fp + face*Nfp, \
			zM + face*Nfp, zP + face*Nfp, fext + face*Nfp, Nfp, Nvar + 1, Ne);
//		EvaluateHydroStaticReconstructValue(Hcrit, fm + face*Nfp, fp + face*Nfp, zM + face*Nfp, zP + face*Nfp, Nfp, Nvar + 1, Ne);
		EvaluateVerticalFaceSurfFlux(FluxM + face*Nfp, fm + face*Nfp, nx + face*Nfp, ny + face*Nfp, &gra, Hcrit, Nfp, Nvar, Ne);
		EvaluateVerticalFaceNumFlux(FluxS + face*Nfp, fm + face*Nfp, fp + face*Nfp, \
			nx + face*Nfp, ny + face*Nfp, &gra, Hcrit, Nfp, Nvar, Ne);
	}
	memcpy(OutputBEFluxM, FluxM, Nfp*Ne*Nvar*sizeof(double));
	memcpy(OutputBEFluxS, FluxS, Nfp*Ne*Nvar*sizeof(double));

	/*Allocate memory for contribution to RHS due to boundary edge facial integral, and calculate 
	contribution to RHS due to boundary edge facial integral in strong form manner*/
	double *BERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(BERHS, 0, Np*K*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < Ne; face++){
			StrongFormBoundaryEdgeRHS(face, FToE, Np, Nfp, FToN1, FluxM + field*Ne*Nfp, FluxS + field*Ne*Nfp, BEJs, BEMb, BERHS + field*Np*K);
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
	free(fm);
	free(fp);
	free(FluxM);
	free(FluxS);
	free(zM);
	free(zP);
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	           /**************************************Bottom Edge Part*******************************************************/
	                              /********************************************************************/
	Ne = BotENe, Nfp = BotENfp;
	FToN1 = BotEFToN1, FToN2 = BotEFToN2, FToE = BotEFToE, nz = BotEnz;
	fm = malloc(Nfp*Ne*(Nvar + 2)*sizeof(double));
	double *omegaM = fm + 2 * Nfp*Ne;
	huM = fm, hvM = fm + Nfp*Ne, hM = fm + 3 * Nfp*Ne;
	fp = malloc(Nfp*Ne*(Nvar + 2)*sizeof(double));
	double *omegaP = fp + 2 * Nfp*Ne;
	huP = fp, hvP = fp + Nfp*Ne, hP = fp + 3 * Nfp*Ne;
	/*Allocate memory for fluxM, fluxP and fluxS, and calculate these flux term*/
	FluxM = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxM, 0, Nfp*Ne*Nvar*sizeof(double));
	FluxP = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxP, 0, Nfp*Ne*Nvar*sizeof(double));
	FluxS = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxS, 0, Nfp*Ne*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < Ne; face++){
		/*Fetch variable fm and fp first*/
		FetchInnerEdgeFacialValue(hM + face*Nfp, hP + face*Nfp, h, FToE + 2 * face, \
			FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		FetchInnerEdgeFacialValue(huM + face*Nfp, huP + face*Nfp, hu, FToE + 2 * face, \
			FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		FetchInnerEdgeFacialValue(hvM + face*Nfp, hvP + face*Nfp, hv, FToE + 2 * face, \
			FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		FetchInnerEdgeFacialValue(omegaM + face*Nfp, omegaP + face*Nfp, omega, FToE + 2 * face, \
			FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
		here 2 stands for the memory occupied by water depth h and omega*/
		for (int field = 2; field < Nvar; field++){
			FetchInnerEdgeFacialValue(fm + (field + 2)*Ne*Nfp + face*Nfp, \
				fp + (field + 2)*Ne*Nfp + face*Nfp, fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				FToE + 2 * face, FToN1 + Nfp*face, FToN2 + Nfp*face, Np, Nfp);
		}
		EvaluateHorizontalFaceSurfFlux(FluxM + face*Nfp, fm + face*Nfp, nz + face*Nfp, Hcrit, Nfp, Nvar, Ne);
		EvaluateHorizontalFaceSurfFlux(FluxP + face*Nfp, fp + face*Nfp, nz + face*Nfp, Hcrit, Nfp, Nvar, Ne);
		EvaluateHorizontalFaceNumFlux(FluxS + face*Nfp, fm + face*Nfp, fp + face*Nfp, \
			nz + face*Nfp, Hcrit, Nfp, Nvar, Ne);
	}
	memcpy(OutputBotEFluxM, FluxM, Nfp*Ne*Nvar*sizeof(double));
	memcpy(OutputBotEFluxP, FluxP, Nfp*Ne*Nvar*sizeof(double));
	memcpy(OutputBotEFluxS, FluxS, Nfp*Ne*Nvar*sizeof(double));
	/*Allocate memory for contribution to RHS due to bottom edge facial integral, and
	calculate contribution to RHS due to bottom edge facial integral in strong form manner*/
	double *BotERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(BotERHS, 0, Np*K*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < Ne; face++){
			StrongFormInnerEdgeRHS(face, FToE, Np, Nfp, FToN1, FToN2, FluxM + field*Ne*Nfp, \
				FluxP + field*Ne*Nfp, FluxS + field*Ne*Nfp, BotEJs, BotEMb, BotERHS + field*Np*K);
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

	free(fm);
	free(fp);
	free(FluxM);
	free(FluxP);
	free(FluxS);
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	             /**************************************Bottom Boundary Edge Part*******************************************************/
	                            /********************************************************************/
	Ne = BotBENe, Nfp = BotBENfp;
	FToN1 = BotBEFToN1, FToE = BotBEFToE, nz = BotBEnz;
	fm = malloc(Nfp*Ne*(Nvar + 2)*sizeof(double));
	huM = fm;
	hvM = fm + Nfp*Ne;
	omegaM = fm + 2 * Nfp*Ne;
	hM = fm + 3 * Nfp*Ne;

	/*Allocate memory for fluxM, fluxP and fluxS, and calculate these flux term*/
	FluxM = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxM, 0, Nfp*Ne*Nvar*sizeof(double));
	FluxS = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxS, 0, Nfp*Ne*Nvar*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < Ne; face++){
		FetchBoundaryEdgeFacialValue(huM + face*Nfp, hu, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(hvM + face*Nfp, hv, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(omegaM + face*Nfp, omega, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(hM + face*Nfp, h, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);

		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included
		2 stands for the memory occupied by water depth h and omega*/
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(fm + (field + 2)*Ne*Nfp + face*Nfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				FToE + 2 * face, FToN1 + Nfp*face, Np, Nfp);
		}

		EvaluateHorizontalFaceSurfFlux(FluxM + face*Nfp, fm + face*Nfp, nz + face*Nfp, Hcrit, Nfp, Nvar, Ne);
	}

	memcpy(OutputBotBEFluxM, FluxM, Nfp*Ne*Nvar*sizeof(double));
	memcpy(OutputBotBEFluxS, FluxS, Nfp*Ne*Nvar*sizeof(double));

	/*Allocate memory for contribution to RHS due to boundary edge facial integral, and calculate
	contribution to RHS due to boundary edge facial integral in strong form manner*/
	double *BotBERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(BotBERHS, 0, Np*K*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < Ne; face++){
			StrongFormBoundaryEdgeRHS(face, FToE, Np, Nfp, FToN1, FluxM + field*Ne*Nfp, FluxS + field*Ne*Nfp, BotBEJs, BotBEMb, BotBERHS + field*Np*K);
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
	free(fm);
	free(FluxM);
	free(FluxS);
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	         /**************************************Surface Boundary Edge Part*******************************************************/
	                               /********************************************************************/

	Ne = SurfBENe, Nfp = SurfBENfp;
	FToN1 = SurfBEFToN1, FToE = SurfBEFToE, nz = SurfBEnz;

	fm = malloc(Nfp*Ne*(Nvar + 2)*sizeof(double));
	huM = fm;
	hvM = fm + Nfp*Ne;
	omegaM = fm + 2 * Nfp*Ne;
	hM = fm + 3 * Nfp*Ne;

	/*Allocate memory for fluxM, fluxP and fluxS, and calculate these flux term*/
	FluxM = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxM, 0, Nfp*Ne*Nvar*sizeof(double));
	/*WE NOTE THAT, FOR THE CONVERGENCE TEST, THIS FLUX IS NOT TAKEN AS ZERO AND SHOULD BE TAKEN FORM THE INPUT*/
	FluxS = malloc(Nfp*Ne*Nvar*sizeof(double));
	memset(FluxS, 0, Nfp*Ne*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int face = 0; face < Ne; face++){
		FetchBoundaryEdgeFacialValue(huM + face*Nfp, hu, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(hvM + face*Nfp, hv, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(omegaM + face*Nfp, omega, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);
		FetchBoundaryEdgeFacialValue(hM + face*Nfp, h, FToE + 2 * face, FToN1 + face*Nfp, Np, Nfp);

		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included
		2 stands for the memory occupied by water depth h and omega*/
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(fm + (field + 2)*Ne*Nfp + face*Nfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				FToE + 2 * face, FToN1 + Nfp*face, Np, Nfp);
		}

		EvaluateHorizontalFaceSurfFlux(FluxM + face*Nfp, fm + face*Nfp, nz + face*Nfp, Hcrit, Nfp, Nvar, Ne);
	}

	memcpy(OutputSurfBEFluxM, FluxM, Nfp*Ne*Nvar*sizeof(double));
	memcpy(OutputSurfBEFluxS, FluxS, Nfp*Ne*Nvar*sizeof(double));

	/*Allocate memory for contribution to RHS due to boundary edge facial integral, and calculate
	contribution to RHS due to boundary edge facial integral in strong form manner*/
	double *SurfBERHS = malloc(Np*K*Nvar*sizeof(double));
	memset(SurfBERHS, 0, Np*K*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int field = 0; field < Nvar; field++){
		for (int face = 0; face < Ne; face++){
			StrongFormBoundaryEdgeRHS(face, FToE, Np, Nfp, FToN1, FluxM + field*Ne*Nfp, FluxS + field*Ne*Nfp, SurfBEJs, SurfBEMb, SurfBERHS + field*Np*K);
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
	free(fm);
	free(FluxM);
	free(FluxS);
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

	memcpy(OutputVolumn, OutputRHS, Np*K*Nvar*sizeof(double));

	memcpy(OutputE, E, Np*K*Nvar*sizeof(double));
	memcpy(OutputG, G, Np*K*Nvar*sizeof(double));
	memcpy(OutputH, H, Np*K*Nvar*sizeof(double));

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