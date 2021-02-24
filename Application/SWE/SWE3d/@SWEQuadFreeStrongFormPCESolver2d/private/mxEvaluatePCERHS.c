#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE.h"
#include "../../../../../NdgMath/NdgMemory.h"

extern double *IEfm2d, *IEfp2d, *IEFluxM2d, *IEFluxP2d, *IEFluxS2d, \
*ERHS2d, *PCEVolumeIntegralX, *PCETempVolumeIntegralX, *PCEVolumeIntegralY, \
*PCETempVolumeIntegralY, *BEfm2d, *BEzM2d, *BEfp2d, *BEzP2d, \
*BEFluxS2d, *BEFluxM2d, *PCETempFacialIntegral;

extern char *PCEInitialized;

void MyExit()
{
	if (!strcmp("True", PCEInitialized)){
		PCEMemoryDeAllocation();
		PCEInitialized = "False";
	}
	return;
}

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

	/*Properties contained in two dimensional inner edge*/
	mxArray *TempIENe2d = mxGetField(prhs[1], 0, "Ne");
	int IENe2d = (int)mxGetScalar(TempIENe2d);
	mxArray *TempIENfp2d = mxGetField(prhs[1], 0, "Nfp");
	int IENfp2d = (int)mxGetScalar(TempIENfp2d);
	mxArray *TempIEFToE2d = mxGetField(prhs[1], 0, "FToE");
	double *IEFToE2d = mxGetPr(TempIEFToE2d);
	mxArray *TempIEFToF2d = mxGetField(prhs[1], 0, "FToF");
	double *IEFToF2d = mxGetPr(TempIEFToF2d);    
	mxArray *TempIEFToN12d = mxGetField(prhs[1], 0, "FToN1");
	double *IEFToN12d = mxGetPr(TempIEFToN12d);
	mxArray *TempIEFToN22d = mxGetField(prhs[1], 0, "FToN2");
	double *IEFToN22d = mxGetPr(TempIEFToN22d);
	mxArray *TempIEnx2d = mxGetField(prhs[1], 0, "nx");
	double *IEnx2d = mxGetPr(TempIEnx2d);
	mxArray *TempIEny2d = mxGetField(prhs[1], 0, "ny");
	double *IEny2d = mxGetPr(TempIEny2d);
	mxArray *TempIEJs2d = mxGetField(prhs[1], 0, "Js");
	double *IEJs2d = mxGetPr(TempIEJs2d);
	mxArray *TempIEMb2d = mxGetField(prhs[1], 0, "M");
	double *IEMb2d = mxGetPr(TempIEMb2d);

	/*Properties contained in two dimensional boundary edge*/
	mxArray *TempBENe2d = mxGetField(prhs[2], 0, "Ne");
	int BENe2d = (int)mxGetScalar(TempBENe2d);
	mxArray *TempBENfp2d = mxGetField(prhs[2], 0, "Nfp");
	int BENfp2d = (int)mxGetScalar(TempBENfp2d);
	mxArray *TempBEFToE2d = mxGetField(prhs[2], 0, "FToE");
	double *BEFToE2d = mxGetPr(TempBEFToE2d);
    mxArray *TempBEFToF2d = mxGetField(prhs[2], 0, "FToF");
	double *BEFToF2d = mxGetPr(TempBEFToF2d);
	mxArray *TempBEFToN12d = mxGetField(prhs[2], 0, "FToN1");
	double *BEFToN12d = mxGetPr(TempBEFToN12d);
	mxArray *TempBEnx2d = mxGetField(prhs[2], 0, "nx");
	double *BEnx2d = mxGetPr(TempBEnx2d);
	mxArray *TempBEny2d = mxGetField(prhs[2], 0, "ny");
	double *BEny2d = mxGetPr(TempBEny2d);
	mxArray *TempBEJs2d = mxGetField(prhs[2], 0, "Js");
	double *BEJs2d = mxGetPr(TempBEJs2d);
	mxArray *TempBEMb2d = mxGetField(prhs[2], 0, "M");
	double *BEMb2d = mxGetPr(TempBEMb2d);

	/*Data contained in two-dimensional standard cell*/
	mxArray *TempDr2d = mxGetField(prhs[3], 0, "Dr");
	double *Dr2d = mxGetPr(TempDr2d);
	mxArray *TempDs2d = mxGetField(prhs[3], 0, "Ds");
	double *Ds2d = mxGetPr(TempDs2d);
    mxArray *TempNface = mxGetField(prhs[3], 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempinvM2d = mxGetField(prhs[3], 0, "invM");
	double *invM2d = mxGetPr(TempinvM2d);

	signed char *ftype2d = (signed char *)mxGetData(prhs[4]);

	/*Data contained in two dimensional physical field*/
	double *fphys2d = mxGetPr(prhs[5]);
	double *h2d = fphys2d;
	double *hu2d = fphys2d + Np2d*K2d;
	double *hv2d = fphys2d + 2 * Np2d*K2d;
	double *z2d = fphys2d + 3 * Np2d*K2d;

	/*The two-dimensional external value*/
	double *fext2d = mxGetPr(prhs[6]);
//	double *BEhE2d = fext2d, *BEhuE2d = fext2d + BENe2d * BENfp2d, *BEhvE2d = fext2d + 2 * BENe2d * BENfp2d;

	double gra = mxGetScalar(prhs[7]);
	double Hcrit = (double)mxGetScalar(prhs[8]);


	plhs[0] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
	double *RHS = mxGetPr(plhs[0]);

	if (!strcmp("False", PCEInitialized)){

		PCEMemoryAllocation( IENfp2d, IENe2d, Np2d, K2d, Nface, BENe2d, BENfp2d);

	}


	double *IEhuM2d = IEfm2d, *IEhvM2d = IEfm2d + IENfp2d * IENe2d, \
		*IEhM2d = IEfm2d + 2 * IENfp2d*IENe2d;
	
	double *IEhuP2d = IEfp2d, *IEhvP2d = IEfp2d + IENfp2d * IENe2d, \
		*IEhP2d = IEfp2d + 2 * IENfp2d*IENe2d;
	
	memset(IEFluxM2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(IEFluxP2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(IEFluxS2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(ERHS2d, 0, Np2d*K2d*Nface*sizeof(double));


	ptrdiff_t np = Np2d;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
		GetVolumnIntegral2d(PCEVolumeIntegralX + k*Np2d, PCETempVolumeIntegralX + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hu2d + k*Np2d, &np, &zero, &np, rx2d + k*Np2d, sx2d + k*Np2d);
		/*$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$*/
		GetVolumnIntegral2d(PCEVolumeIntegralY + k*Np2d, PCETempVolumeIntegralY + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hv2d + k*Np2d, &np, &zero, &np, ry2d + k*Np2d, sy2d + k*Np2d);

		Add(RHS + k*Np2d, PCEVolumeIntegralX + k*Np2d, PCEVolumeIntegralY + k*Np2d, Np2d);
	}
	/*Two dimensional inner edge flux part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		FetchInnerEdgeFacialValue(IEhM2d + e*IENfp2d, IEhP2d + e*IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhuM2d + e*IENfp2d, IEhuP2d + e*IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhvM2d + e*IENfp2d, IEhvP2d + e*IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		GetFacialFluxTerm2d(IEFluxM2d + e*IENfp2d, IEhuM2d + e*IENfp2d, IEhvM2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetFacialFluxTerm2d(IEFluxP2d + e*IENfp2d, IEhuP2d + e*IENfp2d, IEhvP2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetPCENumericalFluxTerm_HLLC_LAI(IEFluxS2d + e*IENfp2d, IEfm2d + e*IENfp2d, IEfp2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, &gra, Hcrit, IENfp2d, IENe2d);
	}

	
	double *BEhuM2d = BEfm2d, *BEhvM2d = BEfm2d + BENe2d * BENfp2d, \
		*BEhM2d = BEfm2d + 2 * BENe2d * BENfp2d;

	memset(BEFluxS2d, 0, BENe2d*BENfp2d*sizeof(double));
	
	memset(BEFluxM2d, 0, BENe2d*BENfp2d*sizeof(double));
	int Nfield = 2;
	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhM2d + e*BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhuM2d + e*BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhvM2d + e*BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEzM2d + e*BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		ImposeBoundaryCondition(&gra, type, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BEfm2d + e*BENfp2d, BEfp2d + e*BENfp2d, \
			BEzM2d + e*BENfp2d, BEzP2d + e*BENfp2d, fext2d + e*BENfp2d, BENfp2d, Nfield, BENe2d);
		EvaluateHydroStaticReconstructValue(Hcrit, BEfm2d + e*BENfp2d, BEfp2d + e*BENfp2d, BEzM2d + e*BENfp2d, BEzP2d + e*BENfp2d, BENfp2d, Nfield, BENe2d);
		GetFacialFluxTerm2d(BEFluxM2d + e*BENfp2d, BEhuM2d + e*BENfp2d, BEhvM2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BENfp2d);
		GetPCENumericalFluxTerm_HLLC_LAI(BEFluxS2d + e*BENfp2d, BEfm2d + e*BENfp2d, BEfp2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, &gra, Hcrit, BENfp2d, BENe2d);
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, IEFluxM2d, IEFluxP2d, IEFluxS2d, IEJs2d, IEMb2d, ERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif    
	for (int e = 0; e < BENe2d; e++){
		StrongFormBoundaryEdgeRHS(e, BEFToE2d, BEFToF2d, Np2d, K2d, BENfp2d, BEFToN12d, BEFluxM2d, BEFluxS2d, BEJs2d, BEMb2d, ERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k < K2d; k++){
            for(int face=1;face<Nface;face++){
                Add( ERHS2d + k*Np2d, ERHS2d + k*Np2d, ERHS2d + face*Np2d*K2d + k*Np2d, Np2d);
            }
    }   
    

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		MultiEdgeContributionByLiftOperator(ERHS2d + k*Np2d, PCETempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
	}

	/*Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		Minus(RHS + k*Np2d, ERHS2d + k*Np2d, RHS + k*Np2d, Np2d);
	}
}