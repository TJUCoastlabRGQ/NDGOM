#include "mex.h"
#include "mxSWE3d.h"
#include <math.h>
#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE3D.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	double *hurhs = mxGetPr(prhs[0]);
	double *hvrhs = mxGetPr(prhs[1]);
	const mxArray *InnerEdge = prhs[2];
	const mxArray *BoundaryEdge = prhs[3];
	const mxArray *BottomEdge = prhs[4];
	const mxArray *cell = prhs[5];
	const mxArray *mesh = prhs[6];
	double *h = mxGetPr(prhs[7]);
	double *rho = mxGetPr(prhs[8]);
	double gra = mxGetScalar(prhs[9]);
	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *OutRHSInX = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *OutRHSInY = mxGetPr(plhs[1]);
	double *PDPX = malloc(Np*K * sizeof(double));
	double *PDPY = malloc(Np*K * sizeof(double));
	double *PRHOPX = malloc(Np*K * sizeof(double));
	double *PRHOPY = malloc(Np*K * sizeof(double));
	double *PRHOPS = malloc(Np*K * sizeof(double));
	GetFirstOrderPartialDerivativeInHorizontalDirection(PDPX, PDPY, PRHOPX, PRHOPY, h, rho, InnerEdge, cell, mesh);
	GetFirstOrderPartialDerivativeInVerticalDirection(PRHOPS, BottomEdge, h, rho, cell, mesh);

	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);
	mxArray *TempNLayer = mxGetField(mesh, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempZ = mxGetField(mesh, 0, "z");
	double *z = mxGetPr(TempZ);
	mxArray *TempJz = mxGetField(mesh, 0, "Jz");
	double *Jz = mxGetPr(TempJz);

	mxArray *TempNpz = mxGetField(cell, 0, "Npz");
	int Npz = (int)mxGetScalar(TempNpz);
	mxArray *TempNph = mxGetField(cell, 0, "Nph");
	int Nph = (int)mxGetScalar(TempNph);
	mxArray *TempVintU = mxGetField(cell, 0, "VintU");
	double *VintU = mxGetPr(TempVintU);
	mxArray *TempV = mxGetField(cell, 0, "V");
	double *V = mxGetPr(TempV);

	mxArray *TempK2d = mxGetField(BottomEdge, 0, "Ne");
	int K2d = (int)mxGetScalar(TempK2d);

	double *InvV = malloc(Np*Np * sizeof(double));
	memcpy(InvV, V, Np*Np * sizeof(double));
	MatrixInverse(InvV, (ptrdiff_t)Np);

	double *BaroclinicInXPartOne = malloc(Np*K * sizeof(double));
	double *BaroclinicInXPartTwo = malloc(Np*K * sizeof(double));
	double *BaroclinicInYPartOne = malloc(Np*K * sizeof(double));
	double *BaroclinicInYPartTwo = malloc(Np*K * sizeof(double));
	double *BaroclinicInXTempRHS = malloc(Np*K * sizeof(double));
	double *BaroclinicInYTempRHS = malloc(Np*K * sizeof(double));
	double *Baroclinicfmod = malloc(Np*K2d * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int p = 0; p < Np; p++) {
			DotProduct(BaroclinicInXPartOne + k*Np, h + k*Np, PRHOPX + k*Np, Np);
			DotProduct(BaroclinicInXPartTwo + k*Np, z + k*Np, PDPX + k*Np, Np);
			DotProduct(BaroclinicInXPartTwo + k*Np, BaroclinicInXPartTwo + k*Np, PRHOPS + k*Np, Np);
			Minus(BaroclinicInXPartOne + k*Np, BaroclinicInXPartOne + k*Np, BaroclinicInXPartTwo + k*Np, Np);
			MultiplyByConstant(BaroclinicInXPartOne + k*Np, BaroclinicInXPartOne + k*Np, -1*gra/1000.0, Np);

			DotProduct(BaroclinicInYPartOne + k*Np, h + k*Np, PRHOPY + k*Np, Np);
			DotProduct(BaroclinicInYPartTwo + k*Np, z + k*Np, PDPY + k*Np, Np);
			DotProduct(BaroclinicInYPartTwo + k*Np, BaroclinicInYPartTwo + k*Np, PRHOPS + k*Np, Np);
			Minus(BaroclinicInYPartOne + k*Np, BaroclinicInYPartOne + k*Np, BaroclinicInYPartTwo + k*Np, Np);
			MultiplyByConstant(BaroclinicInYPartOne + k*Np, BaroclinicInYPartOne + k*Np, -1 * gra / 1000.0, Np);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		VerticalIntegralFromBottom(BaroclinicInXTempRHS + k*NLayer*Np, BaroclinicInXPartOne + k*NLayer*Np, Jz + k*NLayer*Np, Baroclinicfmod + k*Np, NLayer, (ptrdiff_t)Np, InvV, Nph, Npz, VintU);
		VerticalIntegralFromBottom(BaroclinicInYTempRHS + k*NLayer*Np, BaroclinicInYPartOne + k*NLayer*Np, Jz + k*NLayer*Np, Baroclinicfmod + k*Np, NLayer, (ptrdiff_t)Np, InvV, Nph, Npz, VintU);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int p = 0; p < Np; p++)
		{
			BaroclinicInXTempRHS[k*Np + p] = h[k*Np + p] * BaroclinicInXTempRHS[k*Np + p];
			OutRHSInX[k*Np + p] = hurhs[k*Np + p] + BaroclinicInXTempRHS[k*Np + p];
			BaroclinicInYTempRHS[k*Np + p] = h[k*Np + p] * BaroclinicInYTempRHS[k*Np + p];
			OutRHSInY[k*Np + p] = hvrhs[k*Np + p] + BaroclinicInYTempRHS[k*Np + p];
		}
	}
	free(BaroclinicInXPartOne);
	free(BaroclinicInXPartTwo);
	free(BaroclinicInYPartOne);
	free(BaroclinicInYPartTwo);
	free(BaroclinicInXTempRHS);
	free(BaroclinicInYTempRHS);
	free(Baroclinicfmod);
	free(InvV);
}

void GetFirstOrderPartialDerivativeInVerticalDirection(double *PRHOPS, const mxArray *BottomEdge, double *h, double *rho,\
	const mxArray *cell, const mxArray *mesh) {
	mxArray *Temptz = mxGetField(mesh, 0, "tz");
	double *tz = mxGetPr(Temptz);
	mxArray *TempJ = mxGetField(mesh, 0, "J");
	double *J = mxGetPr(TempJ);
	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempInvM = mxGetField(cell, 0, "invM");
	double *invM = mxGetPr(TempInvM);

	/*For bottom edge object*/
	mxArray *TempBotENe = mxGetField(BottomEdge, 0, "Ne");
	int BotENe = (int)mxGetScalar(TempBotENe);
	mxArray *TempBotENfp = mxGetField(BottomEdge, 0, "Nfp");
	int BotENfp = (int)mxGetScalar(TempBotENfp);
	mxArray *TempBotEMb = mxGetField(BottomEdge, 0, "M");
	double *BotEMb = mxGetPr(TempBotEMb);
	mxArray *TempBotEJs = mxGetField(BottomEdge, 0, "Js");
	double *BotEJs = mxGetPr(TempBotEJs);
	mxArray *TempBotEnz = mxGetField(BottomEdge, 0, "nz");
	double *BotEnz = mxGetPr(TempBotEnz);
	mxArray *TempBotEFToE = mxGetField(BottomEdge, 0, "FToE");
	double *BotEFToE = mxGetPr(TempBotEFToE);
	mxArray *TempBotEFToF = mxGetField(BottomEdge, 0, "FToF");
	double *BotEFToF = mxGetPr(TempBotEFToF);
	mxArray *TempBotEFToN1 = mxGetField(BottomEdge, 0, "FToN1");
	double *BotEFToN1 = mxGetPr(TempBotEFToN1);
	mxArray *TempBotEFToN2 = mxGetField(BottomEdge, 0, "FToN2");
	double *BotEFToN2 = mxGetPr(TempBotEFToN2);

	double *BaroclinicBotEfm = malloc(BotENe*BotENfp * sizeof(double));
	double *BaroclinicBotEfp = malloc(BotENe*BotENfp * sizeof(double));
	double *BaroclinicBotEFluxM = malloc(BotENe*BotENfp * sizeof(double));
	double *BaroclinicBotEFluxP = malloc(BotENe*BotENfp * sizeof(double));
	double *BaroclinicBotEFluxS = malloc(BotENe*BotENfp * sizeof(double));
	double *BaroclinicERHS = malloc(Np*K * Nface * sizeof(double));
	double *BaroclinicTempFacialIntegral = malloc(Np*K * sizeof(double));

	memset(BaroclinicERHS, 0, Np*K * Nface * sizeof(double));

	double *rhoM = BaroclinicBotEfm, *rhoP = BaroclinicBotEfp, *rhoFluxM = BaroclinicBotEFluxM, \
		*rhoFluxP = BaroclinicBotEFluxP, *rhoFluxS = BaroclinicBotEFluxS;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++) {
		/*Fetch variable BotEfm and BotEfp first*/
		FetchInnerEdgeFacialValue(rhoM + face*BotENfp, rhoP + face*BotENfp, rho, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);

		EvaluateVerticalFaceSurfFlux(rhoFluxM + face*BotENfp, rhoM + face*BotENfp, BotEnz + face*BotENfp, BotENfp);

		EvaluateVerticalFaceSurfFlux(rhoFluxP + face*BotENfp, rhoP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);

		EvaluateVerticalFaceNumFlux_Central(rhoFluxS + face*BotENfp, rhoM + face*BotENfp, rhoP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++) {
		for (int field = 0; field < 1; field++) {
			StrongFormInnerEdgeRHS(face, BotEFToE, BotEFToF, Np, K, BotENfp, BotEFToN1, BotEFToN2, BaroclinicBotEFluxM + field*BotENe*BotENfp, \
				BaroclinicBotEFluxP + field*BotENe*BotENfp, BaroclinicBotEFluxS + field*BotENe*BotENfp, BotEJs, BotEMb, BaroclinicERHS + field*Np*K*Nface);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field<1; field++) {
			for (int face = 1; face<Nface; face++) {
				Add(BaroclinicERHS + field*Np*K*Nface + k*Np, BaroclinicERHS + field*Np*K*Nface + k*Np, BaroclinicERHS + field*Np*K*Nface + face*Np*K + k*Np, Np);
			}
		}
	}

	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < 1; field++) {
			MultiEdgeContributionByLiftOperator(BaroclinicERHS + field*Np*K*Nface + k*Np, BaroclinicTempFacialIntegral + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		/*$\bold{t_z}\cdot (Dt*rho)$*/
		GetVolumnIntegral1d(PRHOPS + k*Np, &np, &oneI, &np, &one, \
			Dt, &np, u + k*Np, &np, &zero, &np, tz + k*Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		Minus(PRHOPS + k*Np, PRHOPS + k*Np, BaroclinicERHS + k*Np, Np);
	}

	free(BaroclinicBotEfm);
	free(BaroclinicBotEfp);
	free(BaroclinicBotEFluxM);
	free(BaroclinicBotEFluxP);
	free(BaroclinicBotEFluxS);
	free(BaroclinicERHS);
	free(BaroclinicTempFacialIntegral);
}

void GetFirstOrderPartialDerivativeInHorizontalDirection(double *PDPX, double *PDPY, double *PRHOPX, double *PRHOPY,\
	double *h, double *rho, const mxArray *InnerEdge, const mxArray *cell, const mxArray *mesh) {
	mxArray *Temprx = mxGetField(mesh, 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(mesh, 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(mesh, 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *Tempsy = mxGetField(mesh, 0, "sy");
	double *sy = mxGetPr(Tempsy);
	mxArray *TempJ = mxGetField(mesh, 0, "J");
	double *J = mxGetPr(TempJ);
	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	mxArray *TempDr = mxGetField(cell, 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(cell, 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempInvM = mxGetField(cell, 0, "invM");
	double *invM = mxGetPr(TempInvM);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);

	mxArray *TempIENe = mxGetField(InnerEdge, 0, "Ne");
	int IENe = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp = mxGetField(InnerEdge, 0, "Nfp");
	int IENfp = (int)mxGetScalar(TempIENfp);
	mxArray *TempIEMb = mxGetField(InnerEdge, 0, "M");
	double *IEMb = mxGetPr(TempIEMb);
	mxArray *TempIEJs = mxGetField(InnerEdge, 0, "Js");
	double *IEJs = mxGetPr(TempIEJs);
	mxArray *TempIEnx = mxGetField(InnerEdge, 0, "nx");
	double *IEnx = mxGetPr(TempIEnx);
	mxArray *TempIEny = mxGetField(InnerEdge, 0, "ny");
	double *IEny = mxGetPr(TempIEny);
	mxArray *TempIELAV = mxGetField(InnerEdge, 0, "LAV");
	double *IELAV = mxGetPr(TempIELAV);
	mxArray *TempIEFToE = mxGetField(InnerEdge, 0, "FToE");
	double *IEFToE = mxGetPr(TempIEFToE);
	mxArray *TempIEFToF = mxGetField(InnerEdge, 0, "FToF");
	double *IEFToF = mxGetPr(TempIEFToF);
	mxArray *TempIEFToN1 = mxGetField(InnerEdge, 0, "FToN1");
	double *IEFToN1 = mxGetPr(TempIEFToN1);
	mxArray *TempIEFToN2 = mxGetField(InnerEdge, 0, "FToN2");
	double *IEFToN2 = mxGetPr(TempIEFToN2);

	double *BaroclinicIEfm = malloc(IENe*IENfp * 2 * sizeof(double));
	double *BaroclinicIEfp = malloc(IENe*IENfp * 2 * sizeof(double));
	double *BaroclinicIEfluxM = malloc(IENe*IENfp * 4 * sizeof(double));
	double *BaroclinicIEfluxP = malloc(IENe*IENfp * 4 * sizeof(double));
	double *BaroclinicIEfluxS = malloc(IENe*IENfp * 4 * sizeof(double));
	double *BaroclinicERHS = malloc(4*Np*K*(Nface-2)*sizeof(double));
	double *BaroclinicTempFacialIntegral = malloc(Np*K * sizeof(double));
	double *BaroclinicTempVolumeIntegral = malloc(Np*K * sizeof(double));

	double *hM = BaroclinicIEfm, *rhoM = BaroclinicIEfm + IENe*IENfp;
	double *hP = BaroclinicIEfp, *rhoP = BaroclinicIEfp + IENe*IENfp;
	double *HIEfluxMx = BaroclinicIEfluxM, *rhoIEfluxMx = BaroclinicIEfluxM + IENe*IENfp, \
		*HIEfluxMy = BaroclinicIEfluxM + 2*IENe*IENfp, *rhoIEfluxMy = BaroclinicIEfluxM + 3*IENe*IENfp;
	double *HIEfluxPx = BaroclinicIEfluxP, *rhoIEfluxPx = BaroclinicIEfluxP + IENe*IENfp, \
		*HIEfluxPy = BaroclinicIEfluxP + 2 * IENe*IENfp, *rhoIEfluxPy = BaroclinicIEfluxP + 3 * IENe*IENfp;
	double *HIEfluxSx = BaroclinicIEfluxSx, *rhoIEfluxSx = BaroclinicIEfluxSx + IENe*IENfp, \
	    *HIEfluxSy = BaroclinicIEfluxSy + 2 * IENe*IENfp, *rhoIEfluxSy = BaroclinicIEfluxSy + 3 * IENe*IENfp;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		/*Fetch variable IEfm and IEfp first*/
		FetchInnerEdgeFacialValue(hM + face*IENfp, hP + face*IENfp, h, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		FetchInnerEdgeFacialValue(rhoM + face*IENfp, rhoP + face*IENfp, rho, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);

		EvaluateVerticalFaceSurfFlux(HIEfluxMx + face*IENfp, hM + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateVerticalFaceSurfFlux(HIEfluxPx + face*IENfp, hP + face*IENfp, IEnx + face*IENfp, IENfp);

		EvaluateVerticalFaceSurfFlux(HIEfluxMy + face*IENfp, hM + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateVerticalFaceSurfFlux(HIEfluxPy + face*IENfp, hP + face*IENfp, IEny + face*IENfp, IENfp);

		EvaluateVerticalFaceSurfFlux(rhoIEfluxMx + face*IENfp, rhoM + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateVerticalFaceSurfFlux(rhoIEfluxPx + face*IENfp, rhoP + face*IENfp, IEnx + face*IENfp, IENfp);

		EvaluateVerticalFaceSurfFlux(rhoIEfluxMy + face*IENfp, rhoM + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateVerticalFaceSurfFlux(rhoIEfluxPy + face*IENfp, rhoP + face*IENfp, IEny + face*IENfp, IENfp);

		EvaluateVerticalFaceNumFlux_Central(HIEfluxSx + face*IENfp, hM + face*IENfp, hP + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateVerticalFaceNumFlux_Central(HIEfluxSy + face*IENfp, hM + face*IENfp, hP + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateVerticalFaceNumFlux_Central(rhoIEfluxSx + face*IENfp, rhoM + face*IENfp, rhoP + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateVerticalFaceNumFlux_Central(rhoIEfluxSy + face*IENfp, rhoM + face*IENfp, rhoP + face*IENfp, IEny + face*IENfp, IENfp);
	}

	memset(BaroclinicERHS, 0, Np * K * 4 * (Nface - 2) * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		for (int field = 0; field < 4; field++) {
			StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, BaroclinicIEfluxM + field*IENe*IENfp, \
				BaroclinicIEfluxP + field*IENe*IENfp, BaroclinicIEfluxS + field*IENe*IENfp, IEJs, IEMb, BaroclinicERHS + field*Np*K*(Nface - 2));
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field<4; field++) {
			for (int face = 1; face<Nface - 2; face++) {
				Add(BaroclinicERHS + field*Np*K*(Nface - 2) + k*Np, BaroclinicERHS + field*Np*K*(Nface - 2) + k*Np, BaroclinicERHS + field*Np*K*(Nface - 2) + face*Np*K + k*Np, Np);
			}
		}
	}

	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < 4; field++) {
			MultiEdgeContributionByLiftOperator(BaroclinicERHS + field*Np*K*(Nface - 2) + k*Np, BaroclinicTempFacialIntegral + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		/*$\bold{r_x}\cdot (Dr*h)+\bold{s_x}\cdot (Ds*h)$*/
		GetVolumnIntegral2d(PDPX + k*Np, BaroclinicTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, h + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);

		GetVolumnIntegral2d(PDPY + k*Np, BaroclinicTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, h + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		GetVolumnIntegral2d(PRHOPX + k*Np, BaroclinicTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, rho + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);

		GetVolumnIntegral2d(PRHOPY + k*Np, BaroclinicTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, rho + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {

		Minus(PDPX + k*Np, PDPX + k*Np, BaroclinicERHS + k*Np, Np);

		Minus(PDPY + k*Np, PDPY + k*Np, BaroclinicERHS + 2*Np*K*(Nface - 2) + k*Np, Np);

		Minus(PRHOPX + k*Np, PRHOPX + k*Np, BaroclinicERHS + Np*K*(Nface - 2) + k*Np, Np);

		Minus(PRHOPY + k*Np, PRHOPY + k*Np, BaroclinicERHS + 4 * Np*K*(Nface - 2) + k*Np, Np);

	}

	free(BaroclinicIEfm);
	free(BaroclinicIEfp);
	free(BaroclinicIEfluxM);
	free(BaroclinicIEfluxP);
	free(BaroclinicIEfluxS);
	free(BaroclinicERHS);
	free(BaroclinicTempFacialIntegral);
	free(BaroclinicTempVolumeIntegral);
}

void EvaluateVerticalFaceSurfFlux(double *dest, double *source, double *vector, int size) {
	for (int i = 0; i < size; i++)
	{
		dest[i] = source[i] * vector[i];
	}
}

void EvaluateVerticalFaceNumFlux_Central(double *dest, double *sourcefm, double *sourcefm, double *vector, int size) {
	for (int i = 0; i < size; i++)
	{
		dest[i] = (sourcefm[i] + sourcefp[i]) / 2.0 * vector[i];
	}
}