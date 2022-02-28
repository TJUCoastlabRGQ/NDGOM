
#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE.h"
#include "../../../../../NdgMath/NdgMemory.h"

extern double *VSrhs2d, *VSIEfm2d, *VSIEfp2d, *VSIEFluxM2d, \
*VSIEFluxP2d, *VSIEFluxS2d, *VSERHS2d, *VSVolumeIntegralX, \
*VSTempVolumeIntegralX, *VSVolumeIntegralY, *VSTempVolumeIntegralY, \
*VSBEfm2d, *VSBEzM2d, *VSBEfp2d, *VSBEzP2d, *VSBEFluxS2d, \
*VSBEFluxM2d, *VSTempFacialIntegral, *VSfield2d, *VSrhs3d, \
*VSIEfm3d, *VSIEfp3d, *VSIEFluxM3d, *VSIEFluxP3d, *VSIEFluxS3d, \
*VSERHS3d, *VSVolumeIntegralX3d, *VSTempVolumeIntegralX3d, \
*VSVolumeIntegralY3d, *VSTempVolumeIntegralY3d, *VSBEfm3d, \
*VSBEzM3d, *VSBEfp3d, *VSBEzP3d, *VSBEFluxS3d, *VSBEFluxM3d, \
*VSTempFacialIntegral3d, *VSTempVerticalVelocity, *VSBotVertVelocity ;

extern char *VertVelocityInitialized;

void MyExit()
{
	if (!strcmp("True", VertVelocityInitialized)){
		VertVelocitySolverMemoryDeAllocation();
		VertVelocityInitialized = "False";
	}
	return;
}

void FetchBoundaryData(double *, double *, double *, double *, int);

/* 
* The numerical flux for the 2d part is based on the two-dimensional physical variable, 
* i.e. the 2d model and the 3d model are independent. The final velocity is calculated 
* weakly (i.e. calculated from bottom using the boundary condition and the information 
* of the bottom element) from bottom
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
//	mxArray *Tempftype2d = mxGetField(prhs[3], 0, "ftype");
//	signed char *ftype2d = (signed char *)mxGetData(Tempftype2d);

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
//	mxArray *Tempftype3d = mxGetField(prhs[5], 0, "ftype");
//	signed char *ftype3d = (signed char *)mxGetData(Tempftype3d);
	mxArray *TempBEJs3d = mxGetField(prhs[5], 0, "Js");
	double *BEJs3d = mxGetPr(TempBEJs3d);

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

	double gra = mxGetScalar(prhs[11]);

	/*The two-dimensional external value*/
	double *fext2d = mxGetPr(prhs[12]);
//	double *BEhE2d = fext2d, *BEhuE2d = fext2d + BENe2d * BENfp2d, *BEhvE2d = fext2d + 2 * BENe2d * BENfp2d;
	
	/*The three-dimensional external value*/	
	double *fext3d = mxGetPr(prhs[13]);
//	double *BEhuE3d = fext3d, *BEhvE3d = fext3d + BENe3d*BENfp3d, *BEhE3d = fext3d + 3 * BENe3d*BENfp3d;

	/*The coefficient matrix corresponds to the right hand side of horizontal partial derivative term*/
	double *RHSCoeMatrix = mxGetPr(prhs[14]);
	/*The coefficient matrix corresponds to the bottom vertical velocity, this part also contributes to the right hand side*/
	double *VertCoeMatrix = mxGetPr(prhs[15]);

	/*The interpolation point index at the bottom most face*/
	double *BotEidM = mxGetPr(prhs[16]);
	/*The interpolation point index at the up most face*/
	double *UpEidM = mxGetPr(prhs[17]);
	signed char *ftype2d = (signed char *)mxGetData(prhs[18]);
	signed char *ftype3d = (signed char *)mxGetData(prhs[19]);
    
	 if (!strcmp("False", VertVelocityInitialized)){

		 VertVelocitySolverMemoryAllocation( Np2d, K2d, IENfp2d, IENe2d, Nface, BENe2d, BENfp2d, Np3d, \
			  K3d, IENfp3d, IENe3d, BENe3d, BENfp3d);

	 }
     
	plhs[0] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
	double *VerticalVelocity = mxGetPr(plhs[0]);

/**********************************************************  Two Dimensional Part  *******************************************************************************/
	memset(VSrhs2d, 0, Np2d*K2d*sizeof(double));
	double *IEhuM2d = VSIEfm2d, *IEhvM2d = VSIEfm2d + IENfp2d * IENe2d, \
		*IEhM2d = VSIEfm2d + 2 * IENfp2d*IENe2d;
	double *IEhuP2d = VSIEfp2d, *IEhvP2d = VSIEfp2d + IENfp2d * IENe2d, \
		*IEhP2d = VSIEfp2d + 2 * IENfp2d*IENe2d;
	memset(VSIEFluxM2d, 0, IENfp2d*IENe2d*sizeof(double)); 
	memset(VSIEFluxP2d, 0, IENfp2d*IENe2d*sizeof(double));
	memset(VSIEFluxS2d, 0, IENfp2d*IENe2d*sizeof(double)); 
	memset(VSERHS2d, 0, Np2d*K2d*Nface*sizeof(double));

	ptrdiff_t np = Np2d;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
			/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
		GetVolumnIntegral2d(VSVolumeIntegralX + k*Np2d, VSTempVolumeIntegralX + k*Np2d, &np, &oneI, &np, &one, \
				Dr2d, Ds2d, &np, hu2d + k*Np2d, &np, &zero, &np, rx2d + k*Np2d, sx2d + k*Np2d);
			/*$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$*/
		GetVolumnIntegral2d(VSVolumeIntegralY + k*Np2d, VSTempVolumeIntegralY + k*Np2d, &np, &oneI, &np, &one, \
				Dr2d, Ds2d, &np, hv2d + k*Np2d, &np, &zero, &np, ry2d + k*Np2d, sy2d + k*Np2d);

		Add(VSrhs2d + k*Np2d, VSVolumeIntegralX + k*Np2d, VSVolumeIntegralY + k*Np2d, Np2d);
	}
/*Two dimensional inner edge flux part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		FetchInnerEdgeFacialValue(IEhM2d + e*IENfp2d, IEhP2d + e*IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhuM2d + e*IENfp2d, IEhuP2d + e*IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhvM2d + e*IENfp2d, IEhvP2d + e*IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		GetFacialFluxTerm2d(VSIEFluxM2d + e*IENfp2d, IEhuM2d + e*IENfp2d, IEhvM2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetFacialFluxTerm2d(VSIEFluxP2d + e*IENfp2d, IEhuP2d + e*IENfp2d, IEhvP2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetPCENumericalFluxTerm_HLLC_LAI(VSIEFluxS2d + e*IENfp2d, VSIEfm2d + e*IENfp2d, VSIEfp2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, &gra, Hcrit, IENfp2d, IENe2d);
	}

	double *BEhuM2d = VSBEfm2d, *BEhvM2d = VSBEfm2d + BENe2d * BENfp2d, \
		*BEhM2d = VSBEfm2d + 2 * BENe2d * BENfp2d;
//	double *BEhuP2d = VSBEfp2d, *BEhvP2d = VSBEfp2d + BENe2d * BENfp2d, \
		*BEhP2d = VSBEfp2d + 2 * BENe2d * BENfp2d;
	memset(VSBEFluxS2d, 0, BENe2d*BENfp2d*sizeof(double));
	memset(VSBEFluxM2d, 0, BENe2d*BENfp2d*sizeof(double));
	/*The number of variables is three, since we only consider hu, hv and h when adding the boundary condition*/
	int Nfield = 3;
	/*The following void pointer is added on 08/25/2021 to accomadate the usage of function ImposeBoundaryCondition*/
	double *varFieldIndex = NULL;
	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhM2d + e*BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhuM2d + e*BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhvM2d + e*BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(VSBEzM2d + e*BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		ImposeBoundaryCondition(&gra, type, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, VSBEfm2d + e*BENfp2d, VSBEfp2d + e*BENfp2d, \
			VSBEzM2d + e*BENfp2d, VSBEzP2d + e*BENfp2d, fext2d + e*BENfp2d, BENfp2d, Nfield, BENe2d, varFieldIndex);
		EvaluateHydroStaticReconstructValue(Hcrit, VSBEfm2d + e*BENfp2d, VSBEfp2d + e*BENfp2d, VSBEzM2d + e*BENfp2d, VSBEzP2d + e*BENfp2d, BENfp2d, Nfield, BENe2d);
		GetFacialFluxTerm2d(VSBEFluxM2d + e*BENfp2d, BEhuM2d + e*BENfp2d, BEhvM2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BENfp2d);
		GetPCENumericalFluxTerm_HLLC_LAI(VSBEFluxS2d + e*BENfp2d, VSBEfm2d + e*BENfp2d, VSBEfp2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, &gra, Hcrit, BENfp2d, BENe2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, VSIEFluxM2d, VSIEFluxP2d, VSIEFluxS2d, IEJs2d, IEMb2d, VSERHS2d);
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		StrongFormBoundaryEdgeRHS(e, BEFToE2d, BEFToF2d, Np2d, K2d, BENfp2d, BEFToN12d, VSBEFluxM2d, VSBEFluxS2d, BEJs2d, BEMb2d, VSERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k < K2d; k++){
            for(int face=1; face<Nface; face++){
                Add( VSERHS2d + k*Np2d, VSERHS2d + k*Np2d, VSERHS2d + face*Np2d*K2d + k*Np2d, Np2d);
            }
    }   
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		MultiEdgeContributionByLiftOperator(VSERHS2d + k*Np2d, VSTempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
	}

/*Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		Minus(VSrhs2d + k*Np2d, VSERHS2d + k*Np2d, VSrhs2d + k*Np2d, Np2d);
	}

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		NdgExtend2dField(VSfield2d, VSrhs2d, Np2d, k, Np3d, NLayer, Nz);
	}     

	/**********************************************************  Two Dimensional Part Finished  *******************************************************************************/

	/**********************************************************  Three Dimensional Part  *******************************************************************************/

	memset(VSrhs3d, 0, Np3d*K3d*sizeof(double));
	double *IEhuM3d = VSIEfm3d, *IEhvM3d = VSIEfm3d + IENfp3d * IENe3d, \
		*IEhM3d = VSIEfm3d + 2 * IENfp3d*IENe3d;
	double *IEhuP3d = VSIEfp3d, *IEhvP3d = VSIEfp3d + IENfp3d * IENe3d, \
		*IEhP3d = VSIEfp3d + 2 * IENfp3d*IENe3d;
	memset(VSIEFluxM3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(VSIEFluxP3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(VSIEFluxS3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(VSERHS3d, 0, Np3d*K3d*Nface*sizeof(double));

	np = Np3d;
	oneI = 1;
	one = 1.0, zero = 0.0;

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++){
//		printf("The number of threads in the for loop is: %d\n", omp_get_num_threads());
//		printf("The current thread order is:%d\n", omp_get_thread_num());
		/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
		GetVolumnIntegral2d(VSVolumeIntegralX3d + k*Np3d, VSTempVolumeIntegralX3d + k*Np3d, &np, &oneI, &np, &one, \
			Dr3d, Ds3d, &np, hu3d + k*Np3d, &np, &zero, &np, rx3d + k*Np3d, sx3d + k*Np3d);
		/*$\bold{r_y}\cdot (Dr*hv)+\bold{s_y}\cdot (Ds*hv)$*/
		GetVolumnIntegral2d(VSVolumeIntegralY3d + k*Np3d, VSTempVolumeIntegralY3d + k*Np3d, &np, &oneI, &np, &one, \
			Dr3d, Ds3d, &np, hv3d + k*Np3d, &np, &zero, &np, ry3d + k*Np3d, sy3d + k*Np3d);

		Add(VSrhs3d + k*Np3d, VSVolumeIntegralX3d + k*Np3d, VSVolumeIntegralY3d + k*Np3d, Np3d);
	}
	/*Two dimensional inner edge flux part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		FetchInnerEdgeFacialValue(IEhM3d + e*IENfp3d, IEhP3d + e*IENfp3d, h3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		FetchInnerEdgeFacialValue(IEhuM3d + e*IENfp3d, IEhuP3d + e*IENfp3d, hu3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		FetchInnerEdgeFacialValue(IEhvM3d + e*IENfp3d, IEhvP3d + e*IENfp3d, hv3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		GetFacialFluxTerm2d(VSIEFluxM3d + e*IENfp3d, IEhuM3d + e*IENfp3d, IEhvM3d + e*IENfp3d, IEnx3d + e*IENfp3d, IEny3d + e*IENfp3d, IENfp3d);
		GetFacialFluxTerm2d(VSIEFluxP3d + e*IENfp3d, IEhuP3d + e*IENfp3d, IEhvP3d + e*IENfp3d, IEnx3d + e*IENfp3d, IEny3d + e*IENfp3d, IENfp3d);
		GetPCENumericalFluxTerm_HLLC_LAI(VSIEFluxS3d + e*IENfp3d, VSIEfm3d + e*IENfp3d, VSIEfp3d + e*IENfp3d, IEnx3d + e*IENfp3d, IEny3d + e*IENfp3d, &gra, Hcrit, IENfp3d, IENe3d);
	}


	double *BEhuM3d = VSBEfm3d, *BEhvM3d = VSBEfm3d + BENe3d * BENfp3d, \
		*BEhM3d = VSBEfm3d + 2 * BENe3d * BENfp3d;

	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype3d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhuM3d + e*BENfp3d, hu3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(BEhvM3d + e*BENfp3d, hv3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(BEhM3d + e*BENfp3d, h3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(VSBEzM3d + e*BENfp3d, z3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);

		ImposeBoundaryCondition(&gra, type, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, VSBEfm3d + e*BENfp3d, VSBEfp3d + e*BENfp3d, \
			VSBEzM3d + e*BENfp3d, VSBEzP3d + e*BENfp3d, fext3d + e*BENfp3d, BENfp3d, Nfield, BENe3d, varFieldIndex);
		EvaluateHydroStaticReconstructValue(Hcrit, VSBEfm3d + e*BENfp3d, VSBEfp3d + e*BENfp3d, VSBEzM3d + e*BENfp3d, VSBEzP3d + e*BENfp3d, BENfp3d, Nfield, BENe3d);
		GetFacialFluxTerm2d(VSBEFluxM3d + e*BENfp3d, BEhuM3d + e*BENfp3d, BEhvM3d + e*BENfp3d, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, BENfp3d);
		GetPCENumericalFluxTerm_HLLC_LAI(VSBEFluxS3d + e*BENfp3d, VSBEfm3d + e*BENfp3d, VSBEfp3d + e*BENfp3d, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, &gra, Hcrit, BENfp3d, BENe3d);
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE3d, IEFToF3d, Np3d, K3d, IENfp3d, IEFToN13d, IEFToN23d, VSIEFluxM3d, VSIEFluxP3d, VSIEFluxS3d, IEJs3d, IEMb3d, VSERHS3d);
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		StrongFormBoundaryEdgeRHS(e, BEFToE3d, BEFToF3d, Np3d, K3d, BENfp3d, BEFToN13d, VSBEFluxM3d, VSBEFluxS3d, BEJs3d, BEMb3d, VSERHS3d);
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k < K3d; k++){
        for(int face=1;face<Nface;face++){
            Add( VSERHS3d + k*Np3d, VSERHS3d + k*Np3d, VSERHS3d + face*Np3d*K3d + k*Np3d, Np3d);
        }        
    }

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		MultiEdgeContributionByLiftOperator(VSERHS3d + k*Np3d, VSTempFacialIntegral3d + k*Np3d, &np, &oneI, &np, \
			&one, invM3d, &np, &np, &zero, &np, J3d + k*Np3d, Np3d);
	}

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++){
		Minus(VSrhs3d + k*Np3d, VSERHS3d + k*Np3d, VSrhs3d + k*Np3d, Np3d);
		/*Substract the VSfield2d from VSrhs3d to assemble the final right hand side*/
		Minus(VSrhs3d + k*Np3d, VSrhs3d + k*Np3d, VSfield2d + k*Np3d, Np3d);
	}

	/**********************************************************  Three Dimensional Part Finished  *******************************************************************************/

	memset(VSBotVertVelocity, 0, Np3d*K2d*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		dgemm("N", "N", &np, &oneI, &np, &one, RHSCoeMatrix + k*NLayer*Np3d*Np3d + (NLayer - 1)*Np3d*Np3d, \
			&np, VSrhs3d + k*NLayer*Np3d + (NLayer - 1)*Np3d, &np, &zero, VerticalVelocity + k*NLayer*Np3d + (NLayer - 1)*Np3d, \
			&np);
		for (int L = 1; L < NLayer; L++){
			FetchBoundaryData(VSBotVertVelocity + k*Np3d, VerticalVelocity + k*NLayer*Np3d + (NLayer - L)*Np3d, BotEidM, UpEidM, Np2d);
			dgemm("N", "N", &np, &oneI, &np, &one, RHSCoeMatrix + k*NLayer*Np3d*Np3d + (NLayer - L-1)*Np3d*Np3d, \
				&np, VSrhs3d + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, &np, &zero, VerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, \
				&np);
			dgemm("N", "N", &np, &oneI, &np, &one, VertCoeMatrix + k*NLayer*Np3d*Np3d + (NLayer - L - 1)*Np3d*Np3d, \
				&np, VSBotVertVelocity + k*Np3d, &np, &zero, VSTempVerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, \
				&np);
			Add(VerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, VerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, \
				VSTempVerticalVelocity + k*NLayer*Np3d + (NLayer - L - 1)*Np3d, Np3d);

		}
	}
}

void FetchBoundaryData(double *dest, double *source, double *destIndex, double *sourceIndex, int size)
{
	for (int i = 0; i < size; i++)
		dest[(int)destIndex[i] - 1] = source[(int)sourceIndex[i] - 1];
}