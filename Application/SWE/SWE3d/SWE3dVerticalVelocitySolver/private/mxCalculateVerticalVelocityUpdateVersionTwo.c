
#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE.h"
#include "../../../../../NdgMath/NdgSWE3D.h"
#include "../../../../../NdgMath/NdgMemory.h"

extern double *UpdatedVSrhs2d, *UpdatedVSIEfm2d, *UpdatedVSIEfp2d, *UpdatedVSIEFluxM2d, \
*UpdatedVSIEFluxP2d, *UpdatedVSIEFluxS2d, *UpdatedVSERHS2d, *UpdatedVSVolumeIntegralX, \
*UpdatedVSTempVolumeIntegralX, *UpdatedVSVolumeIntegralY, *UpdatedVSTempVolumeIntegralY, \
*UpdatedVSBEfm2d, *UpdatedVSBEzM2d, *UpdatedVSBEfp2d, *UpdatedVSBEzP2d, *UpdatedVSBEFluxS2d, \
*UpdatedVSBEFluxM2d, *UpdatedVSTempFacialIntegral, *UpdatedVSfield2d, *UpdatedVSrhs3d, \
*UpdatedVSIEfm3d, *UpdatedVSIEfp3d, *UpdatedVSIEFluxM3d, *UpdatedVSIEFluxP3d, *UpdatedVSIEFluxS3d, \
*UpdatedVSERHS3d, *UpdatedVSVolumeIntegralX3d, *UpdatedVSTempVolumeIntegralX3d, \
*UpdatedVSVolumeIntegralY3d, *UpdatedVSTempVolumeIntegralY3d, *UpdatedVSBEfm3d, \
*UpdatedVSBEzM3d, *UpdatedVSBEfp3d, *UpdatedVSBEzP3d, *UpdatedVSBEFluxS3d, *UpdatedVSBEFluxM3d, \
*UpdatedVSTempFacialIntegral3d, *UpdatedVSIEfmod, *UpdatedVSBEfmod, *Updatedfmod;

extern char *UpdatedVertVelocityInitialized;

void MyExit()
{
	if (!strcmp("True", UpdatedVertVelocityInitialized)){
		UpdatedVertVelocitySolverMemoryDeAllocation();
		UpdatedVertVelocityInitialized = "False";
	}
	return;
}

void FetchBoundaryData(double *dest, double *source, double *destIndex, double *sourceIndex, int size)
{
	for (int i = 0; i < size; i++)
		dest[(int)destIndex[i] - 1] = source[(int)sourceIndex[i] - 1];
}

/*
* The numerical flux for the 2d part is calculated by integration of the 3d part,
* The final velocity is calculated weakly (i.e. calculated from bottom
* using the boundary condition and the information of the bottom element)
* from bottom
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	mexAtExit(&MyExit);
	/*Properties contained in mesh2d*/
	mxArray *Temprx2d = mxGetField(prhs[0], 0, "rx");
	double *rx2d = mxGetPr(Temprx2d);
	int Np2d = (int)mxGetM(Temprx2d);
	int K2d = (int)mxGetN(Temprx2d);
	mxArray *Tempsx2d = mxGetField(prhs[0], 0, "sx");
	double *sx2d = mxGetPr(Tempsx2d);
	mxArray *Tempry2d = mxGetField(prhs[0], 0, "ry");
	double *ry2d = mxGetPr(Tempry2d);
	mxArray *Tempsy2d = mxGetField(prhs[0], 0, "sy");
	double *sy2d = mxGetPr(Tempsy2d);
	mxArray *TempJ2d = mxGetField(prhs[0], 0, "J");
	double *J2d = mxGetPr(TempJ2d);

	/*Properties contained in mesh3d*/
	mxArray *Temprx3d = mxGetField(prhs[1], 0, "rx");
	double *rx3d = mxGetPr(Temprx3d);
	int Np3d = (int)mxGetM(Temprx3d);
	int K3d = (int)mxGetN(Temprx3d);
	mxArray *Tempsx3d = mxGetField(prhs[1], 0, "sx");
	double *sx3d = mxGetPr(Tempsx3d);
	mxArray *Tempry3d = mxGetField(prhs[1], 0, "ry");
	double *ry3d = mxGetPr(Tempry3d);
	mxArray *Tempsy3d = mxGetField(prhs[1], 0, "sy");
	double *sy3d = mxGetPr(Tempsy3d);
	mxArray *TempJ3d = mxGetField(prhs[1], 0, "J");
	double *J3d = mxGetPr(TempJ3d);
	mxArray *TempNLayer = mxGetField(prhs[1], 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);
	mxArray *TempJz = mxGetField(prhs[1], 0, "Jz");
	double *Jz = mxGetPr(TempJz);

	/*Properties contained in two dimensional inner edge*/
	mxArray *TempIENe2d = mxGetField(prhs[2], 0, "Ne");
	int IENe2d = (int)mxGetScalar(TempIENe2d);
	mxArray *TempIENfp2d = mxGetField(prhs[2], 0, "Nfp");
	int IENfp2d = (int)mxGetScalar(TempIENfp2d);
	mxArray *TempIEFToE2d = mxGetField(prhs[2], 0, "FToE");
	double *IEFToE2d = mxGetPr(TempIEFToE2d);
	mxArray *TempIEFToF2d = mxGetField(prhs[2], 0, "FToF");
	double *IEFToF2d = mxGetPr(TempIEFToF2d);    
	mxArray *TempIEFToN12d = mxGetField(prhs[2], 0, "FToN1");
	double *IEFToN12d = mxGetPr(TempIEFToN12d);
	mxArray *TempIEFToN22d = mxGetField(prhs[2], 0, "FToN2");
	double *IEFToN22d = mxGetPr(TempIEFToN22d);
	mxArray *TempIEnx2d = mxGetField(prhs[2], 0, "nx");
	double *IEnx2d = mxGetPr(TempIEnx2d);
	mxArray *TempIEny2d = mxGetField(prhs[2], 0, "ny");
	double *IEny2d = mxGetPr(TempIEny2d);
	mxArray *TempIEJs2d = mxGetField(prhs[2], 0, "Js");
	double *IEJs2d = mxGetPr(TempIEJs2d);
	mxArray *TempIEMb2d = mxGetField(prhs[2], 0, "M");
	double *IEMb2d = mxGetPr(TempIEMb2d);

	/*Properties contained in two dimensional boundary edge*/
	mxArray *TempBENe2d = mxGetField(prhs[3], 0, "Ne");
	int BENe2d = (int)mxGetScalar(TempBENe2d);
	mxArray *TempBENfp2d = mxGetField(prhs[3], 0, "Nfp");
	int BENfp2d = (int)mxGetScalar(TempBENfp2d);
	mxArray *TempBEFToE2d = mxGetField(prhs[3], 0, "FToE");
	double *BEFToE2d = mxGetPr(TempBEFToE2d);
	mxArray *TempBEFToF2d = mxGetField(prhs[3], 0, "FToF");
	double *BEFToF2d = mxGetPr(TempBEFToF2d);    
	mxArray *TempBEFToN12d = mxGetField(prhs[3], 0, "FToN1");
	double *BEFToN12d = mxGetPr(TempBEFToN12d);
	mxArray *TempBEnx2d = mxGetField(prhs[3], 0, "nx");
	double *BEnx2d = mxGetPr(TempBEnx2d);
	mxArray *TempBEny2d = mxGetField(prhs[3], 0, "ny");
	double *BEny2d = mxGetPr(TempBEny2d);
	mxArray *TempBEJs2d = mxGetField(prhs[3], 0, "Js");
	double *BEJs2d = mxGetPr(TempBEJs2d);
	mxArray *TempBEMb2d = mxGetField(prhs[3], 0, "M");
	double *BEMb2d = mxGetPr(TempBEMb2d);

	/*Properties contained in three dimensional inner edge*/
	mxArray *TempIENe = mxGetField(prhs[4], 0, "Ne");
	int IENe3d = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp3d = mxGetField(prhs[4], 0, "Nfp");
	int IENfp3d = (int)mxGetScalar(TempIENfp3d);
	mxArray *TempIEFToE3d = mxGetField(prhs[4], 0, "FToE");
	double *IEFToE3d = mxGetPr(TempIEFToE3d);
	mxArray *TempIEFToF3d = mxGetField(prhs[4], 0, "FToF");
	double *IEFToF3d = mxGetPr(TempIEFToF3d);    
	mxArray *TempIEFToN13d = mxGetField(prhs[4], 0, "FToN1");
	double *IEFToN13d = mxGetPr(TempIEFToN13d);
	mxArray *TempIEFToN23d = mxGetField(prhs[4], 0, "FToN2");
	double *IEFToN23d = mxGetPr(TempIEFToN23d);
	mxArray *TempIEnx3d = mxGetField(prhs[4], 0, "nx");
	double *IEnx3d = mxGetPr(TempIEnx3d);
	mxArray *TempIEny3d = mxGetField(prhs[4], 0, "ny");
	double *IEny3d = mxGetPr(TempIEny3d);
	mxArray *TempIEMb3d = mxGetField(prhs[4], 0, "M");
	double *IEMb3d = mxGetPr(TempIEMb3d);
	mxArray *TempIEJs3d = mxGetField(prhs[4], 0, "Js");
	double *IEJs3d = mxGetPr(TempIEJs3d);
	mxArray *TempIEJz3d = mxGetField(prhs[4], 0, "Jz");
	double *IEJz3d = mxGetPr(TempIEJz3d);
	mxArray *TempV1d = mxGetField(prhs[4], 0, "V1d");
	double *V1d = mxGetPr(TempV1d);
	mxArray *TempV2d = mxGetField(prhs[4], 0, "V2d");
	double *V2d = mxGetPr(TempV2d);

	double *InvV2d = malloc(IENfp3d*IENfp3d*sizeof(double));
	memcpy(InvV2d, V2d, IENfp3d*IENfp3d*sizeof(double));
	MatrixInverse(InvV2d, (ptrdiff_t)IENfp3d);

	/*Properties contained in three dimensional boundary edge*/
	mxArray *TempBENe3d = mxGetField(prhs[5], 0, "Ne");
	int BENe3d = (int)mxGetScalar(TempBENe3d);
	mxArray *TempBENfp3d = mxGetField(prhs[5], 0, "Nfp");
	int BENfp3d = (int)mxGetScalar(TempBENfp3d);
	mxArray *TempBEFToE3d = mxGetField(prhs[5], 0, "FToE");
	double *BEFToE3d = mxGetPr(TempBEFToE3d);
	mxArray *TempBEFToF3d = mxGetField(prhs[5], 0, "FToF");
	double *BEFToF3d = mxGetPr(TempBEFToF3d);    
	mxArray *TempBEFToN13d = mxGetField(prhs[5], 0, "FToN1");
	double *BEFToN13d = mxGetPr(TempBEFToN13d);
	mxArray *TempBEnx3d = mxGetField(prhs[5], 0, "nx");
	double *BEnx3d = mxGetPr(TempBEnx3d);
	mxArray *TempBEny3d = mxGetField(prhs[5], 0, "ny");
	double *BEny3d = mxGetPr(TempBEny3d);
	mxArray *TempBEMb3d = mxGetField(prhs[5], 0, "M");
	double *BEMb3d = mxGetPr(TempBEMb3d);
	mxArray *TempBEJs3d = mxGetField(prhs[5], 0, "Js");
	double *BEJs3d = mxGetPr(TempBEJs3d);
	mxArray *TempBEJz3d = mxGetField(prhs[5], 0, "Jz");
	double *BEJz3d = mxGetPr(TempBEJz3d);

	/*Data contained in two dimensional physical field and three dimensional physical field*/
	double *fphys2d = mxGetPr(prhs[6]);
	double *h2d = fphys2d;
	double *hu2d = fphys2d + Np2d*K2d;
	double *hv2d = fphys2d + 2 * Np2d*K2d;
	double *z2d = fphys2d + 3 * Np2d*K2d;
	double *fphys3d = mxGetPr(prhs[7]);
	double *hu3d = fphys3d;
	double *hv3d = fphys3d + Np3d * K3d;
	double *h3d = fphys3d + 3 * Np3d*K3d;
	double *z3d = fphys3d + 5 * Np3d*K3d;

	/*The critical water depth*/
	double Hcrit = (double)mxGetScalar(prhs[8]);
	
	/*Data contained in two-dimensional standard cell*/
	mxArray *TempDr2d = mxGetField(prhs[9], 0, "Dr");
	double *Dr2d = mxGetPr(TempDr2d);
	mxArray *TempDs2d = mxGetField(prhs[9], 0, "Ds");
	double *Ds2d = mxGetPr(TempDs2d);
    mxArray *TempNface = mxGetField(prhs[9], 0, "Nface");
	int Nface = mxGetScalar(TempNface);
	mxArray *TempinvM2d = mxGetField(prhs[9], 0, "invM");
	double *invM2d = mxGetPr(TempinvM2d);

	/*Data contained in three-dimensional standard cell*/
	mxArray *TempDr3d = mxGetField(prhs[10], 0, "Dr");
	double *Dr3d = mxGetPr(TempDr3d);
	mxArray *TempDs3d = mxGetField(prhs[10], 0, "Ds");
	double *Ds3d = mxGetPr(TempDs3d);
	mxArray *TempinvM3d = mxGetField(prhs[10], 0, "invM");
	double *invM3d = mxGetPr(TempinvM3d);
	/*Approximation order in vertical direction*/
	mxArray *TempNz = mxGetField(prhs[10], 0, "Nz");
	int Nz = (int)mxGetScalar(TempNz);
	mxArray *TempNpz = mxGetField(prhs[10], 0, "Npz");
	int Npz = (int)mxGetScalar(TempNpz);
	mxArray *TempVint = mxGetField(prhs[10], 0, "Vint");
	double *Vint = mxGetPr(TempVint);
	mxArray *TempV3d = mxGetField(prhs[10], 0, "V");
	double *V3d = mxGetPr(TempV3d);

	double *InvV3d = malloc(Np3d*Np3d*sizeof(double));
	memcpy(InvV3d, V3d, Np3d*Np3d*sizeof(double));
	MatrixInverse(InvV3d, (ptrdiff_t)Np3d);

	double gra = mxGetScalar(prhs[11]);

	/*The two-dimensional external value*/
	double *fext2d = mxGetPr(prhs[12]);
//	double *BEhE2d = fext2d, *BEhuE2d = fext2d + BENe2d * BENfp2d, *BEhvE2d = fext2d + 2 * BENe2d * BENfp2d;
	
	/*The three-dimensional external value*/	
	double *fext3d = mxGetPr(prhs[13]);

	signed char *ftype2d = (signed char *)mxGetData(prhs[14]);
	signed char *ftype3d = (signed char *)mxGetData(prhs[15]);
    
	if (!strcmp("False", UpdatedVertVelocityInitialized)){

		UpdatedVertVelocitySolverMemoryAllocation(Np2d, K2d, IENfp2d, IENe2d, Nface, BENe2d, BENfp2d, Np3d, \
			  K3d, IENfp3d, IENe3d, BENe3d, BENfp3d);

	 }
     
	plhs[0] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
	double *VerticalVelocity = mxGetPr(plhs[0]);

	/**********************************************************  Three Dimensional Part  *******************************************************************************/

	memset(UpdatedVSrhs3d, 0, Np3d*K3d*sizeof(double));
	double *IEhuM3d = UpdatedVSIEfm3d, *IEhvM3d = UpdatedVSIEfm3d + IENfp3d * IENe3d, \
		*IEhM3d = UpdatedVSIEfm3d + 2 * IENfp3d*IENe3d;
	double *IEhuP3d = UpdatedVSIEfp3d, *IEhvP3d = UpdatedVSIEfp3d + IENfp3d * IENe3d, \
		*IEhP3d = UpdatedVSIEfp3d + 2 * IENfp3d*IENe3d;
	memset(UpdatedVSIEFluxM3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSIEFluxP3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSIEFluxS3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSERHS3d, 0, Np3d*K3d*Nface*sizeof(double));

	ptrdiff_t np = Np3d;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++){
		/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
		GetVolumnIntegral2d(UpdatedVSVolumeIntegralX3d + k*Np3d, UpdatedVSTempVolumeIntegralX3d + k*Np3d, &np, &oneI, &np, &one, \
			Dr3d, Ds3d, &np, hu3d + k*Np3d, &np, &zero, &np, rx3d + k*Np3d, sx3d + k*Np3d);
		/*$\bold{r_y}\cdot (Dr*hv)+\bold{s_y}\cdot (Ds*hv)$*/
		GetVolumnIntegral2d(UpdatedVSVolumeIntegralY3d + k*Np3d, UpdatedVSTempVolumeIntegralY3d + k*Np3d, &np, &oneI, &np, &one, \
			Dr3d, Ds3d, &np, hv3d + k*Np3d, &np, &zero, &np, ry3d + k*Np3d, sy3d + k*Np3d);

		Add(UpdatedVSrhs3d + k*Np3d, UpdatedVSVolumeIntegralX3d + k*Np3d, UpdatedVSVolumeIntegralY3d + k*Np3d, Np3d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		FetchInnerEdgeFacialValue(IEhM3d + e*IENfp3d, IEhP3d + e*IENfp3d, h3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		FetchInnerEdgeFacialValue(IEhuM3d + e*IENfp3d, IEhuP3d + e*IENfp3d, hu3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		FetchInnerEdgeFacialValue(IEhvM3d + e*IENfp3d, IEhvP3d + e*IENfp3d, hv3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		GetFacialFluxTerm2d(UpdatedVSIEFluxM3d + e*IENfp3d, IEhuM3d + e*IENfp3d, IEhvM3d + e*IENfp3d, IEnx3d + e*IENfp3d, IEny3d + e*IENfp3d, IENfp3d);
		GetFacialFluxTerm2d(UpdatedVSIEFluxP3d + e*IENfp3d, IEhuP3d + e*IENfp3d, IEhvP3d + e*IENfp3d, IEnx3d + e*IENfp3d, IEny3d + e*IENfp3d, IENfp3d);
		GetPCENumericalFluxTerm_HLLC_LAI(UpdatedVSIEFluxS3d + e*IENfp3d, UpdatedVSIEfm3d + e*IENfp3d, UpdatedVSIEfp3d + e*IENfp3d, IEnx3d + e*IENfp3d, IEny3d + e*IENfp3d, &gra, Hcrit, IENfp3d, IENe3d);
	}

	memset(UpdatedVSIEfmod, 0, IENe2d*IENfp3d*sizeof(double));

	//void VerticalFaceColumnIntegral(double *dest, double *source, double *fmod, double *InvV2d, int Nfp2d, double *Jz, int Nlayer, double *V1d, int LNfp2d, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		VerticalFaceColumnIntegral(UpdatedVSIEFluxS2d + e*IENfp2d, UpdatedVSIEFluxS3d + e*NLayer*IENfp3d, UpdatedVSIEfmod + e*IENfp3d, InvV2d, (ptrdiff_t)IENfp3d, IEJz3d + e*NLayer*IENfp3d, NLayer, V1d, (ptrdiff_t)IENfp2d, (int)(*(IEFToF3d + e*NLayer * 2)));
	}


	double *BEhuM3d = UpdatedVSBEfm3d, *BEhvM3d = UpdatedVSBEfm3d + BENe3d * BENfp3d, \
		*BEhM3d = UpdatedVSBEfm3d + 2 * BENe3d * BENfp3d;

	int Nfield = 3;

	/*The following void pointer is added on 08/25/2021 to accomadate the usage of function ImposeBoundaryCondition*/
	double *varFieldIndex = NULL;

	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype3d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhuM3d + e*BENfp3d, hu3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(BEhvM3d + e*BENfp3d, hv3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(BEhM3d + e*BENfp3d, h3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(UpdatedVSBEzM3d + e*BENfp3d, z3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);

		ImposeBoundaryCondition(&gra, type, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, UpdatedVSBEfm3d + e*BENfp3d, UpdatedVSBEfp3d + e*BENfp3d, \
			UpdatedVSBEzM3d + e*BENfp3d, UpdatedVSBEzP3d + e*BENfp3d, fext3d + e*BENfp3d, BENfp3d, Nfield, BENe3d, varFieldIndex);
		EvaluateHydroStaticReconstructValue(Hcrit, UpdatedVSBEfm3d + e*BENfp3d, UpdatedVSBEfp3d + e*BENfp3d, UpdatedVSBEzM3d + e*BENfp3d, UpdatedVSBEzP3d + e*BENfp3d, BENfp3d, Nfield, BENe3d);
		GetFacialFluxTerm2d(UpdatedVSBEFluxM3d + e*BENfp3d, BEhuM3d + e*BENfp3d, BEhvM3d + e*BENfp3d, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, BENfp3d);
		GetPCENumericalFluxTerm_HLLC_LAI(UpdatedVSBEFluxS3d + e*BENfp3d, UpdatedVSBEfm3d + e*BENfp3d, UpdatedVSBEfp3d + e*BENfp3d, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, &gra, Hcrit, BENfp3d, BENe3d);
	}

	memset(UpdatedVSBEfmod, 0, BENe2d*BENfp3d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		VerticalFaceColumnIntegral(UpdatedVSBEFluxS2d + e*BENfp2d, UpdatedVSBEFluxS3d + e*NLayer*BENfp3d, UpdatedVSBEfmod + e*BENfp3d, InvV2d, (ptrdiff_t)BENfp3d, BEJz3d + e*NLayer*BENfp3d, NLayer, V1d, (ptrdiff_t)BENfp2d, (int)(*(BEFToF3d + e*NLayer * 2)));
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE3d, IEFToF3d, Np3d, K3d, IENfp3d, IEFToN13d, IEFToN23d, UpdatedVSIEFluxM3d, UpdatedVSIEFluxP3d, UpdatedVSIEFluxS3d, IEJs3d, IEMb3d, UpdatedVSERHS3d);
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		StrongFormBoundaryEdgeRHS(e, BEFToE3d, BEFToF3d, Np3d, K3d, BENfp3d, BEFToN13d, UpdatedVSBEFluxM3d, UpdatedVSBEFluxS3d, BEJs3d, BEMb3d, UpdatedVSERHS3d);
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k < K3d; k++){
        for(int face=1;face<Nface;face++){
			Add(UpdatedVSERHS3d + k*Np3d, UpdatedVSERHS3d + k*Np3d, UpdatedVSERHS3d + face*Np3d*K3d + k*Np3d, Np3d);
        }        
    }
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		MultiEdgeContributionByLiftOperator(UpdatedVSERHS3d + k*Np3d, UpdatedVSTempFacialIntegral3d + k*Np3d, &np, &oneI, &np, \
			&one, invM3d, &np, &np, &zero, &np, J3d + k*Np3d, Np3d);
	}

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++){
		Minus(UpdatedVSrhs3d + k*Np3d, UpdatedVSERHS3d + k*Np3d, UpdatedVSrhs3d + k*Np3d, Np3d);
	}

	/**********************************************************  Three Dimensional Part Finished  *******************************************************************************/

	/**********************************************************  Two Dimensional Part  *******************************************************************************/
	memset(UpdatedVSrhs2d, 0, Np2d*K2d*sizeof(double));
	double *IEhuM2d = UpdatedVSIEfm2d, *IEhvM2d = UpdatedVSIEfm2d + IENfp2d * IENe2d, \
		*IEhM2d = UpdatedVSIEfm2d + 2 * IENfp2d*IENe2d;
	double *IEhuP2d = UpdatedVSIEfp2d, *IEhvP2d = UpdatedVSIEfp2d + IENfp2d * IENe2d, \
		*IEhP2d = UpdatedVSIEfp2d + 2 * IENfp2d*IENe2d;
	memset(UpdatedVSIEFluxM2d, 0, IENfp2d*IENe2d*sizeof(double));
	memset(UpdatedVSIEFluxP2d, 0, IENfp2d*IENe2d*sizeof(double));
	memset(UpdatedVSERHS2d, 0, Np2d*K2d*Nface*sizeof(double));

	np = Np2d;
	oneI = 1;
	one = 1.0, zero = 0.0;
/*
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		//$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$
		GetVolumnIntegral2d(UpdatedVSVolumeIntegralX + k*Np2d, UpdatedVSTempVolumeIntegralX + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hu2d + k*Np2d, &np, &zero, &np, rx2d + k*Np2d, sx2d + k*Np2d);
		//$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$
		GetVolumnIntegral2d(UpdatedVSVolumeIntegralY + k*Np2d, UpdatedVSTempVolumeIntegralY + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hv2d + k*Np2d, &np, &zero, &np, ry2d + k*Np2d, sy2d + k*Np2d);

		Add(UpdatedVSrhs2d + k*Np2d, UpdatedVSVolumeIntegralX + k*Np2d, UpdatedVSVolumeIntegralY + k*Np2d, Np2d);
	}
	*/

	double *Tempfield2d = malloc(Np2d*K2d * sizeof(double));
	memset(Tempfield2d, 0, Np2d*K2d * sizeof(double));
	double *Tempfield3d = malloc(Np3d*K3d * sizeof(double));
	double *fmod = malloc(Np3d*K3d * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++) {
		VerticalColumnIntegralField3d(UpdatedVSrhs2d + i*Np2d, Np2d, V2d, Tempfield2d + i*Np2d, \
			Tempfield3d + i*Np3d*NLayer, UpdatedVSrhs3d + i*Np3d*NLayer, Jz + i*Np3d*NLayer, \
			fmod + i*Np3d*NLayer, InvV3d, Np3d, NLayer);
	}
	free(Tempfield3d);
	free(Tempfield2d);
	free(fmod);

	/*Two dimensional inner edge flux part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		FetchInnerEdgeFacialValue(IEhM2d + e*IENfp2d, IEhP2d + e*IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhuM2d + e*IENfp2d, IEhuP2d + e*IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhvM2d + e*IENfp2d, IEhvP2d + e*IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		GetFacialFluxTerm2d(UpdatedVSIEFluxM2d + e*IENfp2d, IEhuM2d + e*IENfp2d, IEhvM2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetFacialFluxTerm2d(UpdatedVSIEFluxP2d + e*IENfp2d, IEhuP2d + e*IENfp2d, IEhvP2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
	}

	double *BEhuM2d = UpdatedVSBEfm2d, *BEhvM2d = UpdatedVSBEfm2d + BENe2d * BENfp2d, \
		*BEhM2d = UpdatedVSBEfm2d + 2 * BENe2d * BENfp2d;
	memset(UpdatedVSBEFluxM2d, 0, BENe2d*BENfp2d*sizeof(double));

	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhM2d + e*BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhuM2d + e*BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhvM2d + e*BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(UpdatedVSBEzM2d + e*BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		ImposeBoundaryCondition(&gra, type, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, UpdatedVSBEfm2d + e*BENfp2d, UpdatedVSBEfp2d + e*BENfp2d, \
			UpdatedVSBEzM2d + e*BENfp2d, UpdatedVSBEzP2d + e*BENfp2d, fext2d + e*BENfp2d, BENfp2d, Nfield, BENe2d, varFieldIndex);
		EvaluateHydroStaticReconstructValue(Hcrit, UpdatedVSBEfm2d + e*BENfp2d, UpdatedVSBEfp2d + e*BENfp2d, UpdatedVSBEzM2d + e*BENfp2d, UpdatedVSBEzP2d + e*BENfp2d, BENfp2d, Nfield, BENe2d);
		GetFacialFluxTerm2d(UpdatedVSBEFluxM2d + e*BENfp2d, BEhuM2d + e*BENfp2d, BEhvM2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BENfp2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, UpdatedVSIEFluxM2d, UpdatedVSIEFluxP2d, UpdatedVSIEFluxS2d, IEJs2d, IEMb2d, UpdatedVSERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		StrongFormBoundaryEdgeRHS(e, BEFToE2d, BEFToF2d, Np2d, K2d, BENfp2d, BEFToN12d, UpdatedVSBEFluxM2d, UpdatedVSBEFluxS2d, BEJs2d, BEMb2d, UpdatedVSERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		for (int face = 1; face<Nface; face++){
			Add(UpdatedVSERHS2d + k*Np2d, UpdatedVSERHS2d + k*Np2d, UpdatedVSERHS2d + face*Np2d*K2d + k*Np2d, Np2d);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		MultiEdgeContributionByLiftOperator(UpdatedVSERHS2d + k*Np2d, UpdatedVSTempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
	}

	/*Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		Minus(UpdatedVSrhs2d + k*Np2d, UpdatedVSERHS2d + k*Np2d, UpdatedVSrhs2d + k*Np2d, Np2d);
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		NdgExtend2dField(UpdatedVSfield2d, UpdatedVSrhs2d, Np2d, k, Np3d, NLayer, Nz);
	}

	/**********************************************************  Two Dimensional Part Finished  *******************************************************************************/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++){
		/*Substract the VSfield2d from VSrhs3d to assemble the final right hand side*/
		Minus(UpdatedVSrhs3d + k*Np3d, UpdatedVSrhs3d + k*Np3d, UpdatedVSfield2d + k*Np3d, Np3d);
	}


	/*The coefficient matrix corresponds to the right hand side of horizontal partial derivative term*/
	double *RHSCoeMatrix = mxGetPr(prhs[16]);
	/*The coefficient matrix corresponds to the bottom vertical velocity, this part also contributes to the right hand side*/
	double *VertCoeMatrix = mxGetPr(prhs[17]);

	/*The interpolation point index at the bottom most face*/
	double *BotEidM = mxGetPr(prhs[18]);
	/*The interpolation point index at the up most face*/
	double *UpEidM = mxGetPr(prhs[19]);

	double *VSBotVertVelocity = malloc(Np3d*K2d * sizeof(double));
	memset(VSBotVertVelocity, 0, Np3d*K2d * sizeof(double));

	double *VSTempVerticalVelocity = malloc(Np3d*K3d * sizeof(double));
	memset(VSTempVerticalVelocity, 0, Np3d*K3d * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		dgemm("N", "N", &np, &oneI, &np, &one, RHSCoeMatrix + k*NLayer*Np3d*Np3d + (NLayer - 1)*Np3d*Np3d, \
			&np, UpdatedVSrhs3d + k*NLayer*Np3d + (NLayer - 1)*Np3d, &np, &zero, VerticalVelocity + k*NLayer*Np3d + (NLayer - 1)*Np3d, \
			&np);
		for (int L = 1; L < NLayer; L++) {
			FetchBoundaryData(VSBotVertVelocity + k*Np3d, VerticalVelocity + k*NLayer*Np3d + (NLayer - L)*Np3d, BotEidM, UpEidM, Np2d);
			dgemm("N", "N", &np, &oneI, &np, &one, RHSCoeMatrix + k*NLayer*Np3d*Np3d + (NLayer - L - 1)*Np3d*Np3d, \
				&np, UpdatedVSrhs3d + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, &np, &zero, VerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, \
				&np);
			dgemm("N", "N", &np, &oneI, &np, &one, VertCoeMatrix + k*NLayer*Np3d*Np3d + (NLayer - L - 1)*Np3d*Np3d, \
				&np, VSBotVertVelocity + k*Np3d, &np, &zero, VSTempVerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, \
				&np);
			Add(VerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, VerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, \
				VSTempVerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, Np3d);

		}
	}

	free(VSBotVertVelocity);
	free(VSTempVerticalVelocity);
	free(InvV2d);
	free(InvV3d);
}