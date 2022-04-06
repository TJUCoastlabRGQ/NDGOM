#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE.h"
#include "../../../../../NdgMath/NdgSWE3D.h"
#include "../../../../../NdgMath/NdgMemory.h"
#include "stdio.h"
#include <omp.h>

extern double *TempFacialIntegral, *IEfm, *IEfp, *IEFluxM, *IEFluxP, \
*IEFluxS, *ERHS, *BEfm, *BEfp, *AdvzM, *AdvzP, *BEFluxM, \
*BEFluxS, *BotEfm, *BotEfp, *BotEFluxM, *BotEFluxP, *BotEFluxS, \
*BotBEfm, *BotBEFluxM, *BotBEFluxS, *SurfBEfm, *SurfBEFluxM, \
*SurfBEFluxS, *E, *G, *H, *TempVolumeIntegral;

extern char *AdvInitialized;

int timepoint = 0;

void MyExit()
{
	if (!strcmp("True", AdvInitialized)){
		AdvMemoryDeAllocation();
		AdvInitialized = "False";
	}
	return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	mexAtExit(&MyExit);
	/*Parameter reading part*/
	/*For mesh object*/
	mxArray *Temprx = mxGetField(prhs[0], 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(prhs[0], 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(prhs[0], 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *Tempsy = mxGetField(prhs[0], 0, "sy");
	double *sy = mxGetPr(Tempsy);
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
	mxArray *TempNface = mxGetField(prhs[1], 0, "Nface");
	int Nface = mxGetScalar(TempNface);    
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
	mxArray *TempIEFToE = mxGetField(prhs[2], 0, "FToE");
	double *IEFToE = mxGetPr(TempIEFToE);
	mxArray *TempIEFToF = mxGetField(prhs[2], 0, "FToF");
	double *IEFToF = mxGetPr(TempIEFToF);    
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
	mxArray *TempBEFToE = mxGetField(prhs[3], 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToF = mxGetField(prhs[3], 0, "FToF");
	double *BEFToF = mxGetPr(TempBEFToF);    
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
	mxArray *TempBotEnz = mxGetField(prhs[4], 0, "nz");
	double *BotEnz = mxGetPr(TempBotEnz);
	mxArray *TempBotEFToE = mxGetField(prhs[4], 0, "FToE");
	double *BotEFToE = mxGetPr(TempBotEFToE);
    mxArray *TempBotEFToF = mxGetField(prhs[4], 0, "FToF");
	double *BotEFToF = mxGetPr(TempBotEFToF);
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
	mxArray *TempBotBEnz = mxGetField(prhs[5], 0, "nz");
	double *BotBEnz = mxGetPr(TempBotBEnz);
	mxArray *TempBotBEFToE = mxGetField(prhs[5], 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToF = mxGetField(prhs[5], 0, "FToF");
	double *BotBEFToF = mxGetPr(TempBotBEFToF);    
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
	mxArray *TempSurfBEnz = mxGetField(prhs[6], 0, "nz");
	double *SurfBEnz = mxGetPr(TempSurfBEnz);
	mxArray *TempSurfBEFToE = mxGetField(prhs[6], 0, "FToE");
	double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
    mxArray *TempSurfBEFToF = mxGetField(prhs[6], 0, "FToF");
	double *SurfBEFToF = mxGetPr(TempSurfBEFToF);
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

	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;
    
    size_t NdimOut = 3;
	mwSize dimOut[3] = { Np, K, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputRHS = mxGetPr(plhs[0]);

	if (!strcmp("False", AdvInitialized)){
		AdvMemoryAllocation(Np, K, Nvar, IENfp, IENe, Nface, BENfp, BENe, BotENfp, \
		  BotENe, BotBENfp, BotBENe, SurfBENfp, SurfBENe);
	}
	                                                                               /************************  Face Integral Part  ****************************/
/*************************************************************************************************************************************/
	/**************************************Inner Edge Part*******************************************************/
	            /********************************************************************/
	/*Allocate memory for fm and fp defined over inner edges. Here, variables correspond to hu, hv, hT, hS, sediment and other passive transport material, vertical velocity omega is not included*/
	
	double *huM = IEfm, *hvM = IEfm + IENfp*IENe, *hM = IEfm + 2 * IENfp*IENe;
	double *huP = IEfp, *hvP = IEfp + IENfp*IENe, *hP = IEfp + 2 * IENfp*IENe;
	/*Allocate memory for IEFluxM, IEFluxP and IEFluxS, and calculate these flux term*/
	memset(IEFluxM, 0, IENfp*IENe*Nvar*sizeof(double));
	memset(IEFluxP, 0, IENfp*IENe*Nvar*sizeof(double));
	memset(IEFluxS, 0, IENfp*IENe*Nvar*sizeof(double));
    
    //printf("Number of threads is:%d\n",omp_get_max_threads());
//	FILE *fp2;
//	fp2 = fopen("D:\\Sharewithpc\\研究工作\\20220404\\IEPureAdv3d.txt", "a");
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
	for (int face = 0; face < IENe; face++){
 //       printf("The order of thread is:%d\n",omp_get_thread_num());  // This is not thread-safe, can not be used!!
		/*Fetch variable IEfm and IEfp first*/
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
		EvaluateVerticalFaceNumFlux_HLLC_LAI(IEFluxS + face*IENfp, IEfm + face*IENfp, IEfp + face*IENfp, \
			IEnx + face*IENfp, IEny + face*IENfp, &gra, Hcrit, IENfp, Nvar, IENe);
//		for (int p = 0; p < IENfp; p++) {
//			fprintf(fp2, "%16.12f \n", (*(IEFluxS + 2 * IENe*IENfp + face*IENfp + p)) / 10.0);
//		}
	}
//	fclose(fp2);

/*Allocate memory for contribution to RHS due to inner edge facial integral, and
 calculate contribution to RHS due to inner edge facial integral in strong form manner*/
	memset(ERHS, 0, Np*K*Nvar*Nface*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < IENe; face++){
        for (int field = 0; field < Nvar; field++){
			StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, IEFluxM + field*IENe*IENfp,\
				IEFluxP + field*IENe*IENfp, IEFluxS + field*IENe*IENfp, IEJs, IEMb, ERHS + field*Np*K*Nface);
		}
	}
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	            /**************************************Boundary Edge Part*******************************************************/
	                           /********************************************************************/
	/*Allocate memory for fm and fp defined over boundary edges. Here, variables correspond to hu, hv, h, hT, hS, 
	sediment and other passive transport material, vertical velocity omega is not included for boundary edges*/
	huM = BEfm, hvM = BEfm + BENfp*BENe, hM = BEfm + 2 * BENfp*BENe;
	huP = BEfp, hvP = BEfp + BENfp*BENe, hP = BEfp + 2 * BENfp*BENe;
	memset(BEFluxM, 0, BENfp*BENe*Nvar*sizeof(double));
	memset(BEFluxS, 0, BENfp*BENe*Nvar*sizeof(double));

//	FILE *fp, *fp1;
//	fp = fopen("D:\\Sharewithpc\\研究工作\\20220404\\Adv3d.txt","a");
//	fp1 = fopen("D:\\Sharewithpc\\研究工作\\20220404\\PureAdv3d.txt", "a");
//	fprintf(fp, "For time points %d:\n", timepoint);
//	timepoint = timepoint + 1;
	/*Fetch variable BEfm and BEfp first, then impose boundary condition and conduct hydrostatic reconstruction.
	Finally, calculate local flux term, adjacent flux term and numerical flux term*/
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
	for (int face = 0; face < BENe; face++){
		NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
		FetchBoundaryEdgeFacialValue(huM + face*BENfp, hu, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*BENfp, hv, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(AdvzM + face*BENfp, z, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
		here 1 stands for the memory occupied by water depth h*/
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(BEfm + (field + 1)*BENe*BENfp + face*BENfp, \
				 fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
		}
		ImposeBoundaryCondition(&gra, type, BEnx + face*BENfp, BEny + face*BENfp, BEfm + face*BENfp, BEfp + face*BENfp, \
			AdvzM + face*BENfp, AdvzP + face*BENfp, fext + face*BENfp, BENfp, Nvar + 1, BENe);
		EvaluateHydroStaticReconstructValue(Hcrit, BEfm + face*BENfp, BEfp + face*BENfp, AdvzM + face*BENfp, AdvzP + face*BENfp, BENfp, Nvar + 1, BENe);
		EvaluateVerticalFaceSurfFlux(BEFluxM + face*BENfp, BEfm + face*BENfp, BEnx + face*BENfp, BEny + face*BENfp, &gra, Hcrit, BENfp, Nvar, BENe);
		EvaluateVerticalFaceNumFlux_HLLC_LAI(BEFluxS + face*BENfp, BEfm + face*BENfp, BEfp + face*BENfp, \
			BEnx + face*BENfp, BEny + face*BENfp, &gra, Hcrit, BENfp, Nvar, BENe);

//		if (type == NdgEdgeClampedVel) {
//			fprintf(fp,"For face %d:\n", face);
//			fprintf(fp,"For hT, the numerical flux is: \n");
//			for (int p = 0; p < BENfp; p++) {
//				fprintf(fp,"%16.12f \n", (*(BEFluxS + 2 * BENe*BENfp + face*BENfp + p))/10.0);
//				fprintf(fp1, "%16.12f \n", (*(BEFluxS + 2 * BENe*BENfp + face*BENfp + p)) / 10.0);
//			}
//		}
	}
//	fclose(fp);
//	fclose(fp1);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		for (int field = 0; field < Nvar; field++){
			StrongFormBoundaryEdgeRHS(face, BEFToE, BEFToF, Np, K, BENfp, BEFToN1, BEFluxM + field*BENe*BENfp, BEFluxS + field*BENe*BENfp, BEJs, BEMb, ERHS + field*Np*K*Nface);
		}
	}
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	           /**************************************Bottom Edge Part*******************************************************/
	                              /********************************************************************/
	double *omegaM = BotEfm + 2 * BotENfp*BotENe;
	huM = BotEfm, hvM = BotEfm + BotENfp*BotENe, hM = BotEfm + 3 * BotENfp*BotENe;
	double *omegaP = BotEfp + 2 * BotENfp*BotENe;
	huP = BotEfp, hvP = BotEfp + BotENfp*BotENe, hP = BotEfp + 3 * BotENfp*BotENe;
	/*Allocate memory for BotEFluxM, BotEFluxP and BotEFluxS, and calculate these flux term*/
	memset(BotEFluxM, 0, BotENfp*BotENe*Nvar*sizeof(double));
	memset(BotEFluxP, 0, BotENfp*BotENe*Nvar*sizeof(double));
	memset(BotEFluxS, 0, BotENfp*BotENe*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++){
		/*Fetch variable BotEfm and BotEfp first*/
		FetchInnerEdgeFacialValue(hM + face*BotENfp, hP + face*BotENfp, h, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		FetchInnerEdgeFacialValue(huM + face*BotENfp, huP + face*BotENfp, hu, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		FetchInnerEdgeFacialValue(hvM + face*BotENfp, hvP + face*BotENfp, hv, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		FetchInnerEdgeFacialValue(omegaM + face*BotENfp, omegaP + face*BotENfp, omega, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
		here 2 stands for the memory occupied by water depth h and omega*/
		for (int field = 2; field < Nvar; field++){
			FetchInnerEdgeFacialValue(BotEfm + (field + 2)*BotENe*BotENfp + face*BotENfp, \
				BotEfp + (field + 2)*BotENe*BotENfp + face*BotENfp, fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				BotEFToE + 2 * face, BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		}
		EvaluateHorizontalFaceSurfFlux(BotEFluxM + face*BotENfp, BotEfm + face*BotENfp, BotEnz + face*BotENfp, Hcrit, BotENfp, Nvar, BotENe);
		EvaluateHorizontalFaceSurfFlux(BotEFluxP + face*BotENfp, BotEfp + face*BotENfp, BotEnz + face*BotENfp, Hcrit, BotENfp, Nvar, BotENe);
		EvaluateHorizontalFaceNumFlux(BotEFluxS + face*BotENfp, BotEfm + face*BotENfp, BotEfp + face*BotENfp, \
			BotEnz + face*BotENfp, Hcrit, BotENfp, Nvar, BotENe);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < BotENe; face++){
        for (int field = 0; field < Nvar; field++){
			StrongFormInnerEdgeRHS(face, BotEFToE, BotEFToF, Np, K, BotENfp, BotEFToN1, BotEFToN2, BotEFluxM + field*BotENe*BotENfp, \
				BotEFluxP + field*BotENe*BotENfp, BotEFluxS + field*BotENe*BotENfp, BotEJs, BotEMb, ERHS + field*Np*K*Nface );
		}
	}
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	             /**************************************Bottom Boundary Edge Part*******************************************************/
	                            /********************************************************************/
	huM = BotBEfm;
	hvM = BotBEfm + BotBENfp*BotBENe;
	omegaM = BotBEfm + 2 * BotBENfp*BotBENe;
	hM = BotBEfm + 3 * BotBENfp*BotBENe;

	/*Allocate memory for BotBEFluxM, fluxP and BotBEFluxS, and calculate these flux term*/
	memset(BotBEFluxM, 0, BotBENfp*BotBENe*Nvar*sizeof(double));
	memset(BotBEFluxS, 0, BotBENfp*BotBENe*Nvar*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		FetchBoundaryEdgeFacialValue(huM + face*BotBENfp, hu, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*BotBENfp, hv, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(omegaM + face*BotBENfp, omega, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(hM + face*BotBENfp, h, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included
		2 stands for the memory occupied by water depth h and omega*/
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(BotBEfm + (field + 2)*BotBENe*BotBENfp + face*BotBENfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				BotBEFToE + 2 * face, BotBEFToN1 + BotBENfp*face, Np, BotBENfp);
		}

		EvaluateHorizontalFaceSurfFlux(BotBEFluxM + face*BotBENfp, BotBEfm + face*BotBENfp, BotBEnz + face*BotBENfp, Hcrit, BotBENfp, Nvar, BotBENe);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		for (int field = 0; field < Nvar; field++){
			StrongFormBoundaryEdgeRHS(face, BotBEFToE, BotBEFToF, Np, K, BotBENfp, BotBEFToN1, BotBEFluxM + field*BotBENe*BotBENfp, BotBEFluxS + field*BotBENe*BotBENfp, BotBEJs, BotBEMb, ERHS + field*Np*K*Nface);
		}
	}
	/*************************************************************************************************************************************/

	/*************************************************************************************************************************************/
	         /**************************************Surface Boundary Edge Part*******************************************************/
	                               /********************************************************************/
	huM = SurfBEfm;
	hvM = SurfBEfm + SurfBENfp*SurfBENe;
	omegaM = SurfBEfm + 2 * SurfBENfp*SurfBENe;
	hM = SurfBEfm + 3 * SurfBENfp*SurfBENe;

	/*Allocate memory for SurfBEFluxM, fluxP and SurfBEFluxS, and calculate these flux term*/
	memset(SurfBEFluxM, 0, SurfBENfp*SurfBENe*Nvar*sizeof(double));
	/*WE NOTE THAT, FOR THE CONVERGENCE TEST, THIS FLUX IS NOT TAKEN AS ZERO AND SHOULD BE TAKEN FORM THE INPUT*/
	//double *SurfFluxS = mxGetPr(prhs[13]);
	memset(SurfBEFluxS, 0, SurfBENfp*SurfBENe*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		FetchBoundaryEdgeFacialValue(huM + face*SurfBENfp, hu, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*SurfBENfp, hv, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(omegaM + face*SurfBENfp, omega, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(hM + face*SurfBENfp, h, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);

		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included
		2 stands for the memory occupied by water depth h and omega*/
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(SurfBEfm + (field + 2)*SurfBENe*SurfBENfp + face*SurfBENfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				SurfBEFToE + 2 * face, SurfBEFToN1 + SurfBENfp*face, Np, SurfBENfp);
		}

		EvaluateHorizontalFaceSurfFlux(SurfBEFluxM + face*SurfBENfp, SurfBEfm + face*SurfBENfp, SurfBEnz + face*SurfBENfp, Hcrit, SurfBENfp, Nvar, SurfBENe);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < SurfBENe; face++){
        for (int field = 0; field < Nvar; field++){
			StrongFormBoundaryEdgeRHS(face, SurfBEFToE, SurfBEFToF, Np, K, SurfBENfp, SurfBEFToN1, SurfBEFluxM + field*SurfBENe*SurfBENfp, SurfBEFluxS + field*SurfBENe*SurfBENfp, SurfBEJs, SurfBEMb, ERHS + field*Np*K*Nface);
		}
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k < K; k++){
        for(int field=0;field<Nvar;field++){
            for(int face=1;face<Nface;face++){
                Add( ERHS + field*Np*K*Nface+k*Np, ERHS + field*Np*K*Nface + k*Np, ERHS + field*Np*K*Nface + face*Np*K + k*Np, Np);
            }
        }
    }    
    

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < Nvar; field++){
			MultiEdgeContributionByLiftOperator(ERHS + field*Np*K*Nface + k*Np, TempFacialIntegral + field*Np*K + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}

	/***************************************************************************************************************/


	                                                                                    /************************  Volume Integral Part  ****************************/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		EvaluatePrebalanceVolumeTerm(E + k*Np, G + k*Np, H + k*Np, fphys + k*Np, \
			varFieldIndex, Nvar, &gra, Np, K, Hcrit);

		GetVolumnIntegral3d(OutputRHS + k*Np, TempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, Dt, E + k*Np, G + k*Np, H + k*Np, &np, &np, &zero, \
			&np, rx + k*Np, sx + k*Np, ry + k*Np, sy + k*Np, tz + k*Np, Nvar, Np, K);
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int n = 0; n < Nvar; n++){
			Minus(OutputRHS + n*Np*K + k*Np, \
				ERHS + n*Np*K*Nface + k*Np, OutputRHS + n*Np*K + k*Np, Np);
		}
	}

}
