#include "mex.h"
#include <math.h>
#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE3D.h"
#include "../../../../../NdgMath/NdgSWE.h"
#include "../../../../../NdgMath/NdgMemory.h"
#include "eqstate.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double *BaroclinicPDPX, *BaroclinicPDPY, *BaroclinicPRHOPX, *BaroclinicPRHOPY, \
*BaroclinicPRHOPS,*BaroclinicInXPartOne, *BaroclinicInXPartTwo, *BaroclinicInYPartOne, \
*BaroclinicInYPartTwo, *BaroclinicInXTempRHS,*BaroclinicInYTempRHS, *Baroclinicfmod, \
*BaroclinicBotEfm, *BaroclinicBotEfp, *BaroclinicBotEFluxM,*BaroclinicBotEFluxP, \
*BaroclinicBotEFluxS, *BaroclinicIEfm, *BaroclinicIEfp,*BaroclinicIEfluxM, \
*BaroclinicIEfluxP, *BaroclinicIEfluxS, *BaroclinicERHS, *BaroclinicTempFacialIntegral,\
*BaroclinicTempVolumeIntegral, *BaroclinicBEfm, *BaroclinicBEfp, *BaroclinicBEfluxM,\
*BaroclinicBEfluxS;

extern char *BaroclinicPartInitialized;

void GetFirstOrderPartialDerivativeInVerticalDirection(double *, const mxArray *,\
	double *, const mxArray *, const mxArray *);

void GetFirstOrderPartialDerivativeInHorizontalDirection(double *, double *, \
	double *, double *, double *, double *, const mxArray *, const mxArray *,\
	const mxArray *, const mxArray *, double *, char *, double , signed char *,\
	double , double , double , double , double );

void ImposeBcsForRhoAndH(double *, double *, double *, double *, double *, double *, \
	int, int, NdgEdgeType, signed char *, double, double, double, double, double,\
	double);

void EvaluateFaceSurfFlux(double *, double *, double *, int );

void EvaluateFaceNumFlux_Central(double *, double *, double *, double *, int );

void MyExit()
{
	if (!strcmp("True", BaroclinicPartInitialized)) {
		BaroclinicPartMemoryDeAllocation();
		BaroclinicPartInitialized = "False";
	}
	return;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	mexAtExit(&MyExit);
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
	double rho0 = mxGetScalar(prhs[10]);
	char *EosType;
	EosType = mxArrayToString(prhs[11]);
	double *fext3d = mxGetPr(prhs[12]);
	double Hcrit = (double)mxGetScalar(prhs[13]);
	signed char *ftype = (signed char *)mxGetData(prhs[14]);
	double S0 = (double)mxGetScalar(prhs[15]);
	double T0 = (double)mxGetScalar(prhs[16]);
	double alphaT = (double)mxGetScalar(prhs[17]);
	double betaS = (double)mxGetScalar(prhs[18]);

	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);
	mxArray *TempNLayer = mxGetField(mesh, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);
	mxArray *TempZ = mxGetField(mesh, 0, "z");
	double *z = mxGetPr(TempZ);
	mxArray *TempJz = mxGetField(mesh, 0, "Jz");
	double *Jz = mxGetPr(TempJz);

	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);

	mxArray *TempIENe = mxGetField(InnerEdge, 0, "Ne");
	int IENe = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp = mxGetField(InnerEdge, 0, "Nfp");
	int IENfp = (int)mxGetScalar(TempIENfp);

	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBENfp = mxGetField(BoundaryEdge, 0, "Nfp");
	int BENfp = (int)mxGetScalar(TempBENfp);

	mxArray *TempNpz = mxGetField(cell, 0, "Npz");
	int Npz = (int)mxGetScalar(TempNpz);
	mxArray *TempNph = mxGetField(cell, 0, "Nph");
	int Nph = (int)mxGetScalar(TempNph);
	mxArray *TempVintU = mxGetField(cell, 0, "VintU");
	double *VintU = mxGetPr(TempVintU);
	mxArray *TempV = mxGetField(cell, 0, "V");
	double *V = mxGetPr(TempV);

	int K2d = K/ NLayer;

	mxArray *TempBotENe = mxGetField(BottomEdge, 0, "Ne");
	int BotENe = (int)mxGetScalar(TempBotENe);

	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *OutRHSInX = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *OutRHSInY = mxGetPr(plhs[1]);

	if (!strcmp("False", BaroclinicPartInitialized)) {
		BaroclinicPartMemoryAllocation( Np, K, K2d, BotENe, Nph, IENe, IENfp, Nface, BENe, BENfp);
	}

	GetFirstOrderPartialDerivativeInHorizontalDirection(BaroclinicPDPX, BaroclinicPDPY, BaroclinicPRHOPX, BaroclinicPRHOPY, h, rho, InnerEdge, \
		cell, mesh, BoundaryEdge, fext3d, EosType, Hcrit, ftype, rho0, S0, T0, alphaT, betaS);
	GetFirstOrderPartialDerivativeInVerticalDirection(BaroclinicPRHOPS, BottomEdge, rho, cell, mesh);

	double *InvV = malloc(Np*Np * sizeof(double));
	memcpy(InvV, V, Np*Np * sizeof(double));
	MatrixInverse(InvV, (ptrdiff_t)Np);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int p = 0; p < Np; p++) {
			DotProduct(BaroclinicInXPartOne + k*Np, h + k*Np, BaroclinicPRHOPX + k*Np, Np);
			DotProduct(BaroclinicInXPartTwo + k*Np, z + k*Np, BaroclinicPDPX + k*Np, Np);
			DotProduct(BaroclinicInXPartTwo + k*Np, BaroclinicInXPartTwo + k*Np, BaroclinicPRHOPS + k*Np, Np);
			Minus(BaroclinicInXPartOne + k*Np, BaroclinicInXPartOne + k*Np, BaroclinicInXPartTwo + k*Np, Np);
			MultiplyByConstant(BaroclinicInXPartOne + k*Np, BaroclinicInXPartOne + k*Np, -1*gra/ rho0, Np);

			DotProduct(BaroclinicInYPartOne + k*Np, h + k*Np, BaroclinicPRHOPY + k*Np, Np);
			DotProduct(BaroclinicInYPartTwo + k*Np, z + k*Np, BaroclinicPDPY + k*Np, Np);
			DotProduct(BaroclinicInYPartTwo + k*Np, BaroclinicInYPartTwo + k*Np, BaroclinicPRHOPS + k*Np, Np);
			Minus(BaroclinicInYPartOne + k*Np, BaroclinicInYPartOne + k*Np, BaroclinicInYPartTwo + k*Np, Np);
			MultiplyByConstant(BaroclinicInYPartOne + k*Np, BaroclinicInYPartOne + k*Np, -1 * gra / rho0, Np);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		VerticalIntegralFromSurface(BaroclinicInXTempRHS + k*NLayer*Np, BaroclinicInXPartOne + k*NLayer*Np, Jz + k*NLayer*Np, Baroclinicfmod + k*Np, NLayer, (ptrdiff_t)Np, InvV, Nph, Npz, VintU);
		VerticalIntegralFromSurface(BaroclinicInYTempRHS + k*NLayer*Np, BaroclinicInYPartOne + k*NLayer*Np, Jz + k*NLayer*Np, Baroclinicfmod + k*Np, NLayer, (ptrdiff_t)Np, InvV, Nph, Npz, VintU);
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
	free(InvV);
}

void GetFirstOrderPartialDerivativeInVerticalDirection(double *BaroclinicPRHOPS, const mxArray *BottomEdge, double *rho,\
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

		EvaluateFaceSurfFlux(rhoFluxM + face*BotENfp, rhoM + face*BotENfp, BotEnz + face*BotENfp, BotENfp);

		EvaluateFaceSurfFlux(rhoFluxP + face*BotENfp, rhoP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);

		EvaluateFaceNumFlux_Central(rhoFluxS + face*BotENfp, rhoM + face*BotENfp, rhoP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
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
		GetVolumnIntegral1d(BaroclinicPRHOPS + k*Np, &np, &oneI, &np, &one, \
			Dt, &np, rho + k*Np, &np, &zero, &np, tz + k*Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		Minus(BaroclinicPRHOPS + k*Np, BaroclinicPRHOPS + k*Np, BaroclinicERHS + k*Np, Np);
	}

}

void GetFirstOrderPartialDerivativeInHorizontalDirection(double *BaroclinicPDPX, double *BaroclinicPDPY, \
	double *BaroclinicPRHOPX, double *BaroclinicPRHOPY, double *h, double *rho, const mxArray *InnerEdge,\
	const mxArray *cell, const mxArray *mesh, const mxArray *BoundaryEdge, double *fext, char *EosType,\
	double Hcrit, signed char *ftype, double rho0, double S0, double T0, double alphaT, double betaS)
{
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

	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBENfp = mxGetField(BoundaryEdge, 0, "Nfp");
	int BENfp = (int)mxGetScalar(TempBENfp);
	mxArray *TempBEMb = mxGetField(BoundaryEdge, 0, "M");
	double *BEMb = mxGetPr(TempBEMb);
	mxArray *TempBEJs = mxGetField(BoundaryEdge, 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	mxArray *TempBEnx = mxGetField(BoundaryEdge, 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(BoundaryEdge, 0, "ny");
	double *BEny = mxGetPr(TempBEny);
	mxArray *TempBELAV = mxGetField(BoundaryEdge, 0, "LAV");
	double *BELAV = mxGetPr(TempBELAV);
	mxArray *TempBEFToE = mxGetField(BoundaryEdge, 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToF = mxGetField(BoundaryEdge, 0, "FToF");
	double *BEFToF = mxGetPr(TempBEFToF);
	mxArray *TempBEFToN1 = mxGetField(BoundaryEdge, 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);

	double *hM = BaroclinicIEfm, *rhoM = BaroclinicIEfm + IENe*IENfp;
	double *hP = BaroclinicIEfp, *rhoP = BaroclinicIEfp + IENe*IENfp;
	double *HIEfluxMx = BaroclinicIEfluxM, *rhoIEfluxMx = BaroclinicIEfluxM + IENe*IENfp, \
		*HIEfluxMy = BaroclinicIEfluxM + 2*IENe*IENfp, *rhoIEfluxMy = BaroclinicIEfluxM + 3*IENe*IENfp;
	double *HIEfluxPx = BaroclinicIEfluxP, *rhoIEfluxPx = BaroclinicIEfluxP + IENe*IENfp, \
		*HIEfluxPy = BaroclinicIEfluxP + 2 * IENe*IENfp, *rhoIEfluxPy = BaroclinicIEfluxP + 3 * IENe*IENfp;
	double *HIEfluxSx = BaroclinicIEfluxS, *rhoIEfluxSx = BaroclinicIEfluxS + IENe*IENfp, \
	    *HIEfluxSy = BaroclinicIEfluxS + 2 * IENe*IENfp, *rhoIEfluxSy = BaroclinicIEfluxS + 3 * IENe*IENfp;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		/*Fetch variable IEfm and IEfp first*/
		FetchInnerEdgeFacialValue(hM + face*IENfp, hP + face*IENfp, h, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		FetchInnerEdgeFacialValue(rhoM + face*IENfp, rhoP + face*IENfp, rho, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);

		EvaluateFaceSurfFlux(HIEfluxMx + face*IENfp, hM + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateFaceSurfFlux(HIEfluxPx + face*IENfp, hP + face*IENfp, IEnx + face*IENfp, IENfp);

		EvaluateFaceSurfFlux(HIEfluxMy + face*IENfp, hM + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateFaceSurfFlux(HIEfluxPy + face*IENfp, hP + face*IENfp, IEny + face*IENfp, IENfp);

		EvaluateFaceSurfFlux(rhoIEfluxMx + face*IENfp, rhoM + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateFaceSurfFlux(rhoIEfluxPx + face*IENfp, rhoP + face*IENfp, IEnx + face*IENfp, IENfp);

		EvaluateFaceSurfFlux(rhoIEfluxMy + face*IENfp, rhoM + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateFaceSurfFlux(rhoIEfluxPy + face*IENfp, rhoP + face*IENfp, IEny + face*IENfp, IENfp);

		EvaluateFaceNumFlux_Central(HIEfluxSx + face*IENfp, hM + face*IENfp, hP + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateFaceNumFlux_Central(HIEfluxSy + face*IENfp, hM + face*IENfp, hP + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateFaceNumFlux_Central(rhoIEfluxSx + face*IENfp, rhoM + face*IENfp, rhoP + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateFaceNumFlux_Central(rhoIEfluxSy + face*IENfp, rhoM + face*IENfp, rhoP + face*IENfp, IEny + face*IENfp, IENfp);
	}

	hM = BaroclinicBEfm, rhoM = BaroclinicBEfm + BENe*BENfp, hP = BaroclinicBEfp;
	// The second part space store hT first, then the density rho next.
	double *hTrhoP = BaroclinicBEfp + BENe*BENfp, *hSp = BaroclinicBEfp + 2*BENe*BENfp;
	double *HBEfluxMx = BaroclinicBEfluxM, *rhoBEfluxMx = BaroclinicBEfluxM + BENe*BENfp, \
		*HBEfluxMy = BaroclinicBEfluxM + 2 * BENe*BENfp, *rhoBEfluxMy = BaroclinicBEfluxM + 3 * BENe*BENfp,\
		*HBEfluxSx = BaroclinicBEfluxS, *rhoBEfluxSx = BaroclinicBEfluxS + BENe*BENfp, \
		*HBEfluxSy = BaroclinicBEfluxS + 2 * BENe*BENfp, *rhoBEfluxSy = BaroclinicBEfluxS + 3 * BENe*BENfp;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++) {
		NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
		FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(rhoM + face*BENfp, rho, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		ImposeBcsForRhoAndH(hTrhoP + face*BENfp, hSp + face*BENfp, hM + face*BENfp, hP + face*BENfp, rhoM + face*BENfp, \
			fext + face*BENfp, BENfp, BENe, type, EosType, Hcrit, rho0, S0, T0, alphaT, betaS);
		/*Fetch variable IEfm and IEfp first*/
		EvaluateFaceSurfFlux(HBEfluxMx + face*BENfp, hM + face*BENfp, BEnx + face*BENfp, BENfp);
		EvaluateFaceSurfFlux(HBEfluxMy + face*BENfp, hM + face*BENfp, BEny + face*BENfp, BENfp);
		EvaluateFaceSurfFlux(rhoBEfluxMx + face*BENfp, rhoM + face*BENfp, BEnx + face*BENfp, BENfp);
		EvaluateFaceSurfFlux(rhoBEfluxMy + face*BENfp, rhoM + face*BENfp, BEny + face*BENfp, BENfp);
		EvaluateFaceNumFlux_Central(HBEfluxSx + face*BENfp, hM + face*BENfp, hP + face*BENfp, BEnx + face*BENfp, BENfp);
		EvaluateFaceNumFlux_Central(HBEfluxSy + face*BENfp, hM + face*BENfp, hP + face*BENfp, BEny + face*BENfp, BENfp);
		EvaluateFaceNumFlux_Central(rhoBEfluxSx + face*BENfp, rhoM + face*BENfp, hTrhoP + face*BENfp, BEnx + face*BENfp, BENfp);
		EvaluateFaceNumFlux_Central(rhoBEfluxSy + face*BENfp, rhoM + face*BENfp, hTrhoP + face*BENfp, BEny + face*BENfp, BENfp);
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
	for (int face = 0; face < BENe; face++) {
		for (int field = 0; field < 4; field++) {
			StrongFormBoundaryEdgeRHS(face, BEFToE, BEFToF, Np, K, BENfp, BEFToN1, BaroclinicBEfluxM + field*BENe*BENfp,\
				BaroclinicBEfluxS + field*BENe*BENfp, BEJs, BEMb, BaroclinicERHS + field*Np*K*(Nface - 2));
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
		GetVolumnIntegral2d(BaroclinicPDPX + k*Np, BaroclinicTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, h + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);

		GetVolumnIntegral2d(BaroclinicPDPY + k*Np, BaroclinicTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, h + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		GetVolumnIntegral2d(BaroclinicPRHOPX + k*Np, BaroclinicTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, rho + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);

		GetVolumnIntegral2d(BaroclinicPRHOPY + k*Np, BaroclinicTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, rho + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {

		Minus(BaroclinicPDPX + k*Np, BaroclinicPDPX + k*Np, BaroclinicERHS + k*Np, Np);

		Minus(BaroclinicPRHOPX + k*Np, BaroclinicPRHOPX + k*Np, BaroclinicERHS + Np*K*(Nface - 2) + k*Np, Np);

		Minus(BaroclinicPDPY + k*Np, BaroclinicPDPY + k*Np, BaroclinicERHS + 2*Np*K*(Nface - 2) + k*Np, Np);

		Minus(BaroclinicPRHOPY + k*Np, BaroclinicPRHOPY + k*Np, BaroclinicERHS + 3 * Np*K*(Nface - 2) + k*Np, Np);

	}

}

/*Note: Only part of the boundary type is implemented at present*/
void ImposeBcsForRhoAndH(double *srcAndDest, double *hS, double *hM, double *hP, double *rhoM, double *fext, \
	int BENfp, int BENe, NdgEdgeType type, signed char *EosType, double Hcrit, double rho0, double S0,\
	double T0, double alphaT, double betaS) 
{
	if (type == NdgEdgeZeroGrad) {
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
	}
	else if (type == NdgEdgeClampedVel) {
		//T
		DotCriticalDivide(srcAndDest, fext + 3*BENe*BENfp, &Hcrit, hM, BENfp);
		//S
		DotCriticalDivide(hS, fext + 4 * BENe*BENfp, &Hcrit, hM, BENfp);
		if (!strcmp(EosType, "Jackett05")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = hM[i];
				EosByFeistel(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0));
			}
		}
		else if (!strcmp(EosType, "UNESCO83")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = hM[i];
				EosByUNESCO(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0));
			}
		}
		else if (!strcmp(EosType, "Linear")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = hM[i];
				EosByLinear(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0), rho0, T0, S0, alphaT, betaS);
			}
		}
		else {
			printf("Equation of state(EOS) needs to be pointed for this part!\n");
		}
	}
	else if (type == NdgEdgeClampedDepth) {
		//T. The exterior value is stored as hu, hv, h, ht, and hs
		DotCriticalDivide(srcAndDest, fext + 3 * BENe*BENfp, &Hcrit, fext + 2 * BENe*BENfp, BENfp);
		//S
		DotCriticalDivide(hS, fext + 4 * BENe*BENfp, &Hcrit, fext + 2 * BENe*BENfp, BENfp);
		if (!strcmp(EosType, "Jackett05")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = fext[2 * BENe*BENfp +i];
				EosByFeistel(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0));
			}
		}
		else if (!strcmp(EosType, "UNESCO83")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = fext[2 * BENe*BENfp + i];
				EosByUNESCO(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0));
			}
		}
		else if (!strcmp(EosType, "Linear")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = fext[2 * BENe*BENfp + i];
				EosByLinear(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0), rho0, T0, S0, alphaT, betaS);
			}
		}
		else {
			printf("Equation of state(EOS) needs to be pointed for this part!\n");
		}
	}
	else if (type == NdgEdgeClamped) {
		//T
		DotCriticalDivide(srcAndDest, fext + 3 * BENe*BENfp, &Hcrit, fext + 2*BENfp*BENe, BENfp);
		//S
		DotCriticalDivide(hS, fext + 4 * BENe*BENfp, &Hcrit, fext + 2 * BENfp*BENe, BENfp);
		if (!strcmp(EosType, "Jackett05")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = fext[2 * BENfp*BENe + i];
				EosByFeistel(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0));
			}
		}
		else if (!strcmp(EosType, "UNESCO83")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = fext[2 * BENfp*BENe + i];
				EosByUNESCO(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0));
			}
		}
		else if (!strcmp(EosType, "Linear")) {
			for (int i = 0; i < BENfp; i++) {
				hP[i] = fext[2 * BENfp*BENe + i];
				EosByLinear(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0), rho0, T0, S0, alphaT, betaS);
			}
		}
		else {
			printf("Equation of state(EOS) needs to be pointed for this part!\n");
		}
	}
	else if (type == NdgEdgeSlipWall) {
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
	}
	else if (type == NdgEdgeNonSlipWall) {
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
	}
}

void EvaluateFaceSurfFlux(double *dest, double *source, double *vector, int size) {
	for (int i = 0; i < size; i++)
	{
		dest[i] = source[i] * vector[i];
	}
}

void EvaluateFaceNumFlux_Central(double *dest, double *sourcefm, double *sourcefp, double *vector, int size) {
	for (int i = 0; i < size; i++)
	{
		dest[i] = (sourcefm[i] + sourcefp[i]) / 2.0 * vector[i];
	}
}