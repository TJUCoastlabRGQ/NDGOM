#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE.h"
#include "../../../../../NdgMath/NdgSWE3D.h"
#include "../../../../../NdgMath/NdgMemory.h"
#include <stdio.h>

extern double *PCEUpdatedIEfm2d, *PCEUpdatedIEfp2d, *PCEUpdatedIEFluxM2d, *PCEUpdatedIEFluxP2d, *PCEUpdatedIEFluxS2d, \
*PCEUpdatedERHS2d, *PCEUpdatedVolumeIntegralX, *PCEUpdatedTempVolumeIntegralX, *PCEUpdatedVolumeIntegralY, \
*PCEUpdatedTempVolumeIntegralY, *PCEUpdatedBEfm2d, *PCEUpdatedBEzM2d, *PCEUpdatedBEfp2d, *PCEUpdatedBEzP2d, \
*PCEUpdatedBEFluxS2d, *PCEUpdatedBEFluxM2d, *PCEUpdatedPCETempFacialIntegral, *PCEUpdatedIEfmod, *PCEUpdatedBEfmod,\
*PCEUpdatedIEfm3d, *PCEUpdatedIEfp3d, *PCEUpdatedIEFluxS3d, *PCEUpdatedBEFluxS3d, *PCEUpdatedBEfm3d, *PCEUpdatedBEfp3d, \
*PCEUpdatedBEzM3d, *PCEUpdatedBEzP3d;

extern char *PCEUpdatedInitialized;

//int timepoint = 0;

void MyExit()
{
	if (!strcmp("True", PCEUpdatedInitialized)){
		PCEUpdatedMemoryDeAllocation();
		PCEUpdatedInitialized = "False";
	}
	return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexAtExit(&MyExit);
	/*Properties contained in mesh2d*/
	const mxArray *mesh2d = prhs[0];
	mxArray *Temprx2d = mxGetField(mesh2d, 0, "rx");
	double *rx2d = mxGetPr(Temprx2d);
	int Np2d = (int)mxGetM(Temprx2d);
	int K2d = (int)mxGetN(Temprx2d);
	mxArray *Tempsx2d = mxGetField(mesh2d, 0, "sx");
	double *sx2d = mxGetPr(Tempsx2d);
	mxArray *Tempry2d = mxGetField(mesh2d, 0, "ry");
	double *ry2d = mxGetPr(Tempry2d);
	mxArray *Tempsy2d = mxGetField(mesh2d, 0, "sy");
	double *sy2d = mxGetPr(Tempsy2d);
	mxArray *TempJ2d = mxGetField(mesh2d, 0, "J");
	double *J2d = mxGetPr(TempJ2d);

	const mxArray *mesh3d = prhs[1];
	mxArray *TempK3d = mxGetField(mesh3d, 0, "K");
	int K3d = (int)mxGetScalar(TempK3d);
	mxArray *TempNLayer = mxGetField(mesh3d, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);

	const mxArray *cell3d = prhs[2];
	mxArray *TempNp3d = mxGetField(cell3d, 0, "Np");
	int Np3d = (int)mxGetScalar(TempNp3d);

	/*Properties contained in two dimensional inner edge*/
	const mxArray *InnerEdge2d = prhs[3];
	mxArray *TempIENe2d = mxGetField(InnerEdge2d, 0, "Ne");
	int IENe2d = (int)mxGetScalar(TempIENe2d);
	mxArray *TempIENfp2d = mxGetField(InnerEdge2d, 0, "Nfp");
	int IENfp2d = (int)mxGetScalar(TempIENfp2d);
	mxArray *TempIEFToE2d = mxGetField(InnerEdge2d, 0, "FToE");
	double *IEFToE2d = mxGetPr(TempIEFToE2d);
	mxArray *TempIEFToF2d = mxGetField(InnerEdge2d, 0, "FToF");
	double *IEFToF2d = mxGetPr(TempIEFToF2d);    
	mxArray *TempIEFToN12d = mxGetField(InnerEdge2d, 0, "FToN1");
	double *IEFToN12d = mxGetPr(TempIEFToN12d);
	mxArray *TempIEFToN22d = mxGetField(InnerEdge2d, 0, "FToN2");
	double *IEFToN22d = mxGetPr(TempIEFToN22d);
	mxArray *TempIEnx2d = mxGetField(InnerEdge2d, 0, "nx");
	double *IEnx2d = mxGetPr(TempIEnx2d);
	mxArray *TempIEny2d = mxGetField(InnerEdge2d, 0, "ny");
	double *IEny2d = mxGetPr(TempIEny2d);
	mxArray *TempIEJs2d = mxGetField(InnerEdge2d, 0, "Js");
	double *IEJs2d = mxGetPr(TempIEJs2d);
	mxArray *TempIEMb2d = mxGetField(InnerEdge2d, 0, "M");
	double *IEMb2d = mxGetPr(TempIEMb2d);

	/*Properties contained in two dimensional boundary edge*/
	const mxArray *BoundaryEdge2d = prhs[4];
	mxArray *TempBENe2d = mxGetField(BoundaryEdge2d, 0, "Ne");
	int BENe2d = (int)mxGetScalar(TempBENe2d);
	mxArray *TempBENfp2d = mxGetField(BoundaryEdge2d, 0, "Nfp");
	int BENfp2d = (int)mxGetScalar(TempBENfp2d);
	mxArray *TempBEFToE2d = mxGetField(BoundaryEdge2d, 0, "FToE");
	double *BEFToE2d = mxGetPr(TempBEFToE2d);
    mxArray *TempBEFToF2d = mxGetField(BoundaryEdge2d, 0, "FToF");
	double *BEFToF2d = mxGetPr(TempBEFToF2d);
	mxArray *TempBEFToN12d = mxGetField(BoundaryEdge2d, 0, "FToN1");
	double *BEFToN12d = mxGetPr(TempBEFToN12d);
	mxArray *TempBEnx2d = mxGetField(BoundaryEdge2d, 0, "nx");
	double *BEnx2d = mxGetPr(TempBEnx2d);
	mxArray *TempBEny2d = mxGetField(BoundaryEdge2d, 0, "ny");
	double *BEny2d = mxGetPr(TempBEny2d);
	mxArray *TempBEJs2d = mxGetField(BoundaryEdge2d, 0, "Js");
	double *BEJs2d = mxGetPr(TempBEJs2d);
	mxArray *TempBEMb2d = mxGetField(BoundaryEdge2d, 0, "M");
	double *BEMb2d = mxGetPr(TempBEMb2d);

	/*Data contained in two-dimensional standard cell*/
	const mxArray *cell2d = prhs[5];
	mxArray *TempDr2d = mxGetField(cell2d, 0, "Dr");
	double *Dr2d = mxGetPr(TempDr2d);
	mxArray *TempDs2d = mxGetField(cell2d, 0, "Ds");
	double *Ds2d = mxGetPr(TempDs2d);
    mxArray *TempNface = mxGetField(cell2d, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempinvM2d = mxGetField(cell2d, 0, "invM");
	double *invM2d = mxGetPr(TempinvM2d);

	signed char *ftype2d = (signed char *)mxGetData(prhs[6]);

	signed char *ftype3d = (signed char *)mxGetData(prhs[7]);

	/*Data contained in two dimensional physical field*/
	double *fphys2d = mxGetPr(prhs[8]);
	double *h2d = fphys2d;
	double *hu2d = fphys2d + Np2d*K2d;
	double *hv2d = fphys2d + 2 * Np2d*K2d;
	double *z2d = fphys2d + 3 * Np2d*K2d;

	double *fphys3d = mxGetPr(prhs[9]);
	double *hu3d = fphys3d;
	double *hv3d = fphys3d + Np3d * K3d;
	double *h3d = fphys3d + 3 * Np3d*K3d;
	double *z3d = fphys3d + 5 * Np3d*K3d;

	/*The two-dimensional external value*/
	double *fext2d = mxGetPr(prhs[10]);
	double *fext3d = mxGetPr(prhs[11]);
//	double *BEhE2d = fext2d, *BEhuE2d = fext2d + BENe2d * BENfp2d, *BEhvE2d = fext2d + 2 * BENe2d * BENfp2d;

	double gra = mxGetScalar(prhs[12]);
	double Hcrit = mxGetScalar(prhs[13]);

	const mxArray *InnerEdge3d = prhs[14];
	/*Properties contained in three dimensional inner edge*/
	mxArray *TempIENe = mxGetField(InnerEdge3d, 0, "Ne");
	int IENe3d = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp3d = mxGetField(InnerEdge3d, 0, "Nfp");
	int IENfp3d = (int)mxGetScalar(TempIENfp3d);
	mxArray *TempIEFToE3d = mxGetField(InnerEdge3d, 0, "FToE");
	double *IEFToE3d = mxGetPr(TempIEFToE3d);
	mxArray *TempIEFToF3d = mxGetField(InnerEdge3d, 0, "FToF");
	double *IEFToF3d = mxGetPr(TempIEFToF3d);
	mxArray *TempIEFToN13d = mxGetField(InnerEdge3d, 0, "FToN1");
	double *IEFToN13d = mxGetPr(TempIEFToN13d);
	mxArray *TempIEFToN23d = mxGetField(InnerEdge3d, 0, "FToN2");
	double *IEFToN23d = mxGetPr(TempIEFToN23d);
	mxArray *TempIEnx3d = mxGetField(InnerEdge3d, 0, "nx");
	double *IEnx3d = mxGetPr(TempIEnx3d);
	mxArray *TempIEny3d = mxGetField(InnerEdge3d, 0, "ny");
	double *IEny3d = mxGetPr(TempIEny3d);
	mxArray *TempIEMb3d = mxGetField(InnerEdge3d, 0, "M");
	double *IEMb3d = mxGetPr(TempIEMb3d);
	mxArray *TempIEJz3d = mxGetField(InnerEdge3d, 0, "Jz");
	double *IEJz3d = mxGetPr(TempIEJz3d);
	mxArray *TempV1d = mxGetField(InnerEdge3d, 0, "V1d");
	double *V1d = mxGetPr(TempV1d);
	mxArray *TempV2d = mxGetField(InnerEdge3d, 0, "V2d");
	double *V2d = mxGetPr(TempV2d);

	double *InvV2d = malloc(IENfp3d*IENfp3d*sizeof(double));
	memcpy(InvV2d, V2d, IENfp3d*IENfp3d*sizeof(double));
	MatrixInverse(InvV2d, (ptrdiff_t)IENfp3d);

	const mxArray *BoundaryEdge3d = prhs[15];
	/*Properties contained in three dimensional boundary edge*/
	mxArray *TempBENe3d = mxGetField(BoundaryEdge3d, 0, "Ne");
	int BENe3d = (int)mxGetScalar(TempBENe3d);
	mxArray *TempBENfp3d = mxGetField(BoundaryEdge3d, 0, "Nfp");
	int BENfp3d = (int)mxGetScalar(TempBENfp3d);
	mxArray *TempBEFToE3d = mxGetField(BoundaryEdge3d, 0, "FToE");
	double *BEFToE3d = mxGetPr(TempBEFToE3d);
	mxArray *TempBEFToF3d = mxGetField(BoundaryEdge3d, 0, "FToF");
	double *BEFToF3d = mxGetPr(TempBEFToF3d);
	mxArray *TempBEFToN13d = mxGetField(BoundaryEdge3d, 0, "FToN1");
	double *BEFToN13d = mxGetPr(TempBEFToN13d);
	mxArray *TempBEnx3d = mxGetField(BoundaryEdge3d, 0, "nx");
	double *BEnx3d = mxGetPr(TempBEnx3d);
	mxArray *TempBEny3d = mxGetField(BoundaryEdge3d, 0, "ny");
	double *BEny3d = mxGetPr(TempBEny3d);
	mxArray *TempBEMb3d = mxGetField(BoundaryEdge3d, 0, "M");
	double *BEMb3d = mxGetPr(TempBEMb3d);
	mxArray *TempBEJz3d = mxGetField(BoundaryEdge3d, 0, "Jz");
	double *BEJz3d = mxGetPr(TempBEJz3d);


	plhs[0] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
	double *RHS = mxGetPr(plhs[0]);

	if (!strcmp("False", PCEUpdatedInitialized)){
		PCEUpdatedMemoryAllocation(IENfp2d, IENe2d, Np2d, K2d, Nface, BENe2d, BENfp2d,\
			IENfp3d, IENe3d, BENe3d, BENfp3d);
	}


	double *IEhuM2d = PCEUpdatedIEfm2d, *IEhvM2d = PCEUpdatedIEfm2d + IENfp2d * IENe2d, \
		*IEhM2d = PCEUpdatedIEfm2d + 2 * IENfp2d*IENe2d;
	
	double *IEhuP2d = PCEUpdatedIEfp2d, *IEhvP2d = PCEUpdatedIEfp2d + IENfp2d * IENe2d, \
		*IEhP2d = PCEUpdatedIEfp2d + 2 * IENfp2d*IENe2d;
	
	memset(PCEUpdatedIEFluxM2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(PCEUpdatedIEFluxP2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(PCEUpdatedIEFluxS2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(PCEUpdatedERHS2d, 0, Np2d*K2d*Nface*sizeof(double));

	ptrdiff_t np = Np2d, oneI = 1;
	double one = 1.0, zero = 0.0;
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		//$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$
		GetVolumnIntegral2d(PCEUpdatedVolumeIntegralX + k*Np2d, PCEUpdatedTempVolumeIntegralX + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hu2d + k*Np2d, &np, &zero, &np, rx2d + k*Np2d, sx2d + k*Np2d);
		//$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$
		GetVolumnIntegral2d(PCEUpdatedVolumeIntegralY + k*Np2d, PCEUpdatedTempVolumeIntegralY + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hv2d + k*Np2d, &np, &zero, &np, ry2d + k*Np2d, sy2d + k*Np2d);

		Add(RHS + k*Np2d, PCEUpdatedVolumeIntegralX + k*Np2d, PCEUpdatedVolumeIntegralY + k*Np2d, Np2d);
	}
	
	/*Two dimensional inner edge flux part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		FetchInnerEdgeFacialValue(IEhM2d + e*IENfp2d, IEhP2d + e*IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhuM2d + e*IENfp2d, IEhuP2d + e*IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhvM2d + e*IENfp2d, IEhvP2d + e*IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		GetFacialFluxTerm2d(PCEUpdatedIEFluxM2d + e*IENfp2d, IEhuM2d + e*IENfp2d, IEhvM2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetFacialFluxTerm2d(PCEUpdatedIEFluxP2d + e*IENfp2d, IEhuP2d + e*IENfp2d, IEhvP2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
	}

	
	double *BEhuM2d = PCEUpdatedBEfm2d, *BEhvM2d = PCEUpdatedBEfm2d + BENe2d * BENfp2d, \
		*BEhM2d = PCEUpdatedBEfm2d + 2 * BENe2d * BENfp2d;

	memset(PCEUpdatedBEFluxS2d, 0, BENe2d*BENfp2d*sizeof(double));
	
	memset(PCEUpdatedBEFluxM2d, 0, BENe2d*BENfp2d*sizeof(double));

	int Nfield = 3;
	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
	/*The following void pointer is added on 08/25/2021 to accomadate the usage of function ImposeBoundaryCondition*/
	double *varFieldIndex = NULL;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhM2d + e*BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhuM2d + e*BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhvM2d + e*BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(PCEUpdatedBEzM2d + e*BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		ImposeBoundaryCondition(&gra, type, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, PCEUpdatedBEfm2d + e*BENfp2d, PCEUpdatedBEfp2d + e*BENfp2d, \
			PCEUpdatedBEzM2d + e*BENfp2d, PCEUpdatedBEzP2d + e*BENfp2d, fext2d + e*BENfp2d, BENfp2d, Nfield, BENe2d, varFieldIndex);
		EvaluateHydroStaticReconstructValue(Hcrit, PCEUpdatedBEfm2d + e*BENfp2d, PCEUpdatedBEfp2d + e*BENfp2d, PCEUpdatedBEzM2d + e*BENfp2d, PCEUpdatedBEzP2d + e*BENfp2d, BENfp2d, Nfield, BENe2d);
		GetFacialFluxTerm2d(PCEUpdatedBEFluxM2d + e*BENfp2d, BEhuM2d + e*BENfp2d, BEhvM2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BENfp2d);
	}

	double *IEhuM3d = PCEUpdatedIEfm3d, *IEhvM3d = PCEUpdatedIEfm3d + IENfp3d * IENe3d, \
		*IEhM3d = PCEUpdatedIEfm3d + 2 * IENfp3d*IENe3d;
	double *IEhuP3d = PCEUpdatedIEfp3d, *IEhvP3d = PCEUpdatedIEfp3d + IENfp3d * IENe3d, \
		*IEhP3d = PCEUpdatedIEfp3d + 2 * IENfp3d*IENe3d;
	memset(PCEUpdatedIEFluxS3d, 0, IENfp3d*IENe3d*sizeof(double));

//	FILE *fp2;
//	fp2 = fopen("D:\\Sharewithpc\\研究工作\\20220404\\IEPurePCE.txt", "a");

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		FetchInnerEdgeFacialValue(IEhM3d + e*IENfp3d, IEhP3d + e*IENfp3d, h3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		FetchInnerEdgeFacialValue(IEhuM3d + e*IENfp3d, IEhuP3d + e*IENfp3d, hu3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		FetchInnerEdgeFacialValue(IEhvM3d + e*IENfp3d, IEhvP3d + e*IENfp3d, hv3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		GetPCENumericalFluxTerm_HLLC_LAI(PCEUpdatedIEFluxS3d + e*IENfp3d, PCEUpdatedIEfm3d + e*IENfp3d, PCEUpdatedIEfp3d + e*IENfp3d, IEnx3d + e*IENfp3d, IEny3d + e*IENfp3d, &gra, Hcrit, IENfp3d, IENe3d);
//		for (int p = 0; p < IENfp3d; p++) {
//			fprintf(fp2, "%16.12f \n", *(PCEUpdatedIEFluxS3d + e*IENfp3d + p));
//		}
	}
//	fclose(fp2);

	memset(PCEUpdatedIEfmod, 0, IENe2d*IENfp3d*sizeof(double));

	//void VerticalFaceColumnIntegral(double *dest, double *source, double *fmod, double *InvV2d, int Nfp2d, double *Jz, int Nlayer, double *V1d, int LNfp2d, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		VerticalFaceColumnIntegral(PCEUpdatedIEFluxS2d + e*IENfp2d, PCEUpdatedIEFluxS3d + e*NLayer*IENfp3d, PCEUpdatedIEfmod + e*IENfp3d, InvV2d, (ptrdiff_t)IENfp3d, IEJz3d + e*NLayer*IENfp3d, NLayer, V1d, (ptrdiff_t)IENfp2d, (int)(*(IEFToF3d + e*NLayer * 2)));
	}


	double *BEhuM3d = PCEUpdatedBEfm3d, *BEhvM3d = PCEUpdatedBEfm3d + BENe3d * BENfp3d, \
		*BEhM3d = PCEUpdatedBEfm3d + 2 * BENe3d * BENfp3d;

//	FILE *fp, *fp1;
//	fp = fopen("D:\\Sharewithpc\\研究工作\\20220404\\PCE.txt", "a");
//	fp1 = fopen("D:\\Sharewithpc\\研究工作\\20220404\\PurePCE.txt", "a");
//	fprintf(fp, "For time points %d:\n", timepoint);
//	timepoint = timepoint + 1;

	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype3d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhuM3d + e*BENfp3d, hu3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(BEhvM3d + e*BENfp3d, hv3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(BEhM3d + e*BENfp3d, h3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(PCEUpdatedBEzM3d + e*BENfp3d, z3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);

		ImposeBoundaryCondition(&gra, type, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, PCEUpdatedBEfm3d + e*BENfp3d, PCEUpdatedBEfp3d + e*BENfp3d, \
			PCEUpdatedBEzM3d + e*BENfp3d, PCEUpdatedBEzP3d + e*BENfp3d, fext3d + e*BENfp3d, BENfp3d, Nfield, BENe3d, varFieldIndex);
		EvaluateHydroStaticReconstructValue(Hcrit, PCEUpdatedBEfm3d + e*BENfp3d, PCEUpdatedBEfp3d + e*BENfp3d, PCEUpdatedBEzM3d + e*BENfp3d, PCEUpdatedBEzP3d + e*BENfp3d, BENfp3d, Nfield, BENe3d);
		GetPCENumericalFluxTerm_HLLC_LAI(PCEUpdatedBEFluxS3d + e*BENfp3d, PCEUpdatedBEfm3d + e*BENfp3d, PCEUpdatedBEfp3d + e*BENfp3d, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, &gra, Hcrit, BENfp3d, BENe3d);
		
//		if (type == NdgEdgeClampedVel) {
//			fprintf(fp, "For face %d:\n", e);
//			fprintf(fp, "For h, the numerical flux is: \n");
//			for (int p = 0; p < BENfp3d; p++) {
//				fprintf(fp, "%16.12f \n", *(PCEUpdatedBEFluxS3d + e*BENfp3d + p));
//				fprintf(fp1, "%16.12f \n", *(PCEUpdatedBEFluxS3d + e*BENfp3d + p));
//			}
//		}
	}
//	fclose(fp);
//	fclose(fp1);

	memset(PCEUpdatedBEfmod, 0, BENe2d*BENfp3d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		VerticalFaceColumnIntegral(PCEUpdatedBEFluxS2d + e*BENfp2d, PCEUpdatedBEFluxS3d + e*NLayer*BENfp3d, PCEUpdatedBEfmod + e*BENfp3d, InvV2d, (ptrdiff_t)BENfp3d, BEJz3d + e*NLayer*BENfp3d, NLayer, V1d, (ptrdiff_t)BENfp2d, (int)(*(BEFToF3d + e*NLayer * 2)));
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, PCEUpdatedIEFluxM2d, PCEUpdatedIEFluxP2d, PCEUpdatedIEFluxS2d, IEJs2d, IEMb2d, PCEUpdatedERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif    
	for (int e = 0; e < BENe2d; e++){
		StrongFormBoundaryEdgeRHS(e, BEFToE2d, BEFToF2d, Np2d, K2d, BENfp2d, BEFToN12d, PCEUpdatedBEFluxM2d, PCEUpdatedBEFluxS2d, BEJs2d, BEMb2d, PCEUpdatedERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k < K2d; k++){
		for (int face = 1; face<Nface; face++){
			Add(PCEUpdatedERHS2d + k*Np2d, PCEUpdatedERHS2d + k*Np2d, PCEUpdatedERHS2d + face*Np2d*K2d + k*Np2d, Np2d);
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		MultiEdgeContributionByLiftOperator(PCEUpdatedERHS2d + k*Np2d, PCEUpdatedPCETempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
	}

	/*Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
//		memset(PCEUpdatedERHS2d, 0.0, K2d*Np2d * sizeof(double));
		Minus(RHS + k*Np2d, PCEUpdatedERHS2d + k*Np2d, RHS + k*Np2d, Np2d);
	}

	free(InvV2d);

}