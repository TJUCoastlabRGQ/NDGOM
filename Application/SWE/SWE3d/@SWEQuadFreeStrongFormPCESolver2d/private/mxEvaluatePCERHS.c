#include "..\..\..\..\..\NdgMath\NdgMath.h"
#include "..\..\..\..\..\NdgMath\NdgSWE.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
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
	double *BEhE2d = fext2d, *BEhuE2d = fext2d + BENe2d * BENfp2d, *BEhvE2d = fext2d + 2 * BENe2d * BENfp2d;

	double gra = mxGetScalar(prhs[7]);
	double Hcrit = (double)mxGetScalar(prhs[8]);


	plhs[0] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
	double *RHS = mxGetPr(plhs[0]);

	double *IEfm2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	double *IEhM2d = IEfm2d, *IEhuM2d = IEfm2d + IENfp2d * IENe2d, \
		*IEhvM2d = IEfm2d + 2 * IENfp2d*IENe2d;
	double *IEfp2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	double *IEhP2d = IEfp2d, *IEhuP2d = IEfp2d + IENfp2d * IENe2d, \
		*IEhvP2d = IEfp2d + 2 * IENfp2d*IENe2d;
	double *IEFluxM2d = malloc(IENfp2d*IENe2d*sizeof(double));
	double *IEFluxP2d = malloc(IENfp2d*IENe2d*sizeof(double));
	double *IEFluxS2d = malloc(IENfp2d*IENe2d*sizeof(double));
	double *IERHS2d = malloc(Np2d*K2d*sizeof(double));
	memset(IERHS2d, 0, Np2d*K2d*sizeof(double));
	double *BERHS2d = malloc(Np2d*K2d*sizeof(double));
	memset(BERHS2d, 0, Np2d*K2d*sizeof(double));
	double *VolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	double *TempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	double *VolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	double *TempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));

	ptrdiff_t np = Np2d;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K2d; k++){
		/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
		GetVolumnIntegral2d(VolumeIntegralX + k*Np2d, TempVolumeIntegralX + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hu2d + k*Np2d, &np, &zero, &np, rx2d + k*Np2d, sx2d + k*Np2d);
		/*$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$*/
		GetVolumnIntegral2d(VolumeIntegralY + k*Np2d, TempVolumeIntegralY + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hv2d + k*Np2d, &np, &zero, &np, ry2d + k*Np2d, sy2d + k*Np2d);

		Add(RHS + k*Np2d, VolumeIntegralX + k*Np2d, VolumeIntegralY + k*Np2d, Np2d);
	}
	/*Two dimensional inner edge flux part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int e = 0; e < IENe2d; e++){
		FetchInnerEdgeFacialValue(IEhM2d + e*IENfp2d, IEhP2d + e*IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhuM2d + e*IENfp2d, IEhuP2d + e*IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhvM2d + e*IENfp2d, IEhvP2d + e*IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		GetFacialFluxTerm2d(IEFluxM2d + e*IENfp2d, IEhuM2d + e*IENfp2d, IEhvM2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetFacialFluxTerm2d(IEFluxP2d + e*IENfp2d, IEhuP2d + e*IENfp2d, IEhvP2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetPCENumericalFluxTerm(IEFluxS2d + e*IENfp2d, IEhM2d + e*IENfp2d, IEhuM2d + e*IENfp2d, IEhvM2d + e*IENfp2d, IEhP2d + e*IENfp2d, IEhuP2d + e*IENfp2d, \
			IEhvP2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d, Hcrit, gra);
	}

	double *BEfm2d = malloc(BENe2d * BENfp2d * 4 * sizeof(double));
	double *BEhM2d = BEfm2d, *BEhuM2d = BEfm2d + BENe2d * BENfp2d, \
		*BEhvM2d = BEfm2d + 2 * BENe2d * BENfp2d, *BEzM2d = BEfm2d + 3 * BENe2d * BENfp2d;
	double *BEfp2d = malloc(BENe2d * BENfp2d * 4 * sizeof(double));
	double *BEhP2d = BEfp2d, *BEhuP2d = BEfp2d + BENe2d * BENfp2d, \
		*BEhvP2d = BEfp2d + 2 * BENe2d * BENfp2d, *BEzP2d = BEfp2d + 3 * BENe2d * BENfp2d;
	double *BEFluxS2d = malloc(BENe2d*BENfp2d*sizeof(double));
	double *BEFluxM2d = malloc(BENe2d*BENfp2d*sizeof(double));
	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int e = 0; e < BENe2d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhM2d + e*BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhuM2d + e*BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhvM2d + e*BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEzM2d + e*BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		for (int p = 0; p < IENfp2d; p++){
			ImposeBoundaryCondition(&gra, type, BEnx2d + e*BENfp2d + p, BEny2d + e*BENfp2d + p, BEhM2d + e*BENfp2d + p, \
				BEhuM2d + e*BENfp2d + p, BEhvM2d + e*BENfp2d + p, BEzM2d + e*BENfp2d + p, BEhE2d + e*BENfp2d + p, \
				BEhuE2d + e*BENfp2d + p, BEhvE2d + e*BENfp2d + p, BEhP2d + e*BENfp2d + p, BEhuP2d + e*BENfp2d + p, \
				BEhvP2d + e*BENfp2d + p, BEzP2d + e*BENfp2d + p);
			EvaluateHydroStaticReconstructValue(Hcrit, BEhM2d + e*BENfp2d + p, BEhuM2d + e*BENfp2d + p, BEhvM2d + e*BENfp2d + p, \
				BEzM2d + e*BENfp2d + p, BEhP2d + e*BENfp2d + p, BEhuP2d + e*BENfp2d + p, BEhvP2d + e*BENfp2d + p, BEzP2d + e*BENfp2d + p);
		}
		GetFacialFluxTerm2d(BEFluxM2d + e*IENfp2d, BEhuM2d + e*BENfp2d, BEhvM2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BENfp2d);
		GetPCENumericalFluxTerm(BEFluxS2d + e*BENfp2d, BEhM2d + e*BENfp2d, BEhuM2d + e*BENfp2d, BEhvM2d + e*BENfp2d, BEhP2d + e*BENfp2d, BEhuP2d + e*BENfp2d, \
			BEhvP2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BENfp2d, Hcrit, gra);
	}

	for (int e = 0; e < IENe2d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE2d, Np2d, IENfp2d, IEFToN12d, IEFToN22d, IEFluxM2d, IEFluxP2d, IEFluxS2d, IEJs2d, IEMb2d, IERHS2d);
	}

	for (int e = 0; e < BENe2d; e++){
		StrongFormBoundaryEdgeRHS(e, BEFToE2d, Np2d, BENfp2d, BEFToN12d, BEFluxM2d, BEFluxS2d, BEJs2d, BEMb2d, BERHS2d);
	}


	double *TempFacialIntegral = malloc(Np2d*K2d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K2d; k++) {
		MultiEdgeContributionByLiftOperator(IERHS2d + k*Np2d, TempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
		MultiEdgeContributionByLiftOperator(BERHS2d + k*Np2d, TempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
	}

	/*Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K2d; k++){
		Minus(RHS + k*Np2d, IERHS2d + k*Np2d, RHS + k*Np2d, Np2d);
		Add(RHS + k*Np2d, BERHS2d + k*Np2d, RHS + k*Np2d, Np2d);
	}

	free(IEfm2d);
	free(IEfp2d);
	free(BEfm2d);
	free(BEfp2d);
	free(IEFluxM2d);
	free(IEFluxP2d);
	free(IEFluxS2d);
	free(BEFluxS2d);
	free(BEFluxM2d);
	free(IERHS2d);
	free(BERHS2d);
	free(VolumeIntegralX);
	free(VolumeIntegralY);
	free(TempVolumeIntegralX);
	free(TempVolumeIntegralY);
	free(TempFacialIntegral);
}