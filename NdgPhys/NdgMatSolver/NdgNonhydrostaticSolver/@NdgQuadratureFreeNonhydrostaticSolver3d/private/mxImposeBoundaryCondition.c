
#include "SWENonhydrostatic3d.h"

#include "stdio.h"

extern double *ImposeBCsK33, *ImposeBCsInvSquaHeight, *BETau, *ImposeBCsNewmannData, *ImposeBCsWx, *ImposeBCsWy, \
*ImposeBCsWxRHS2d, *ImposeBCsWyRHS2d, *ImposeBCsWIEFluxMx2d, *ImposeBCsWIEFluxMy2d, *ImposeBCsWIEFluxPx2d, *ImposeBCsWIEFluxPy2d, \
*ImposeBCsWIEFluxSx2d, *ImposeBCsWIEFluxSy2d, *ImposeBCsVolumeIntegralX, *ImposeBCsTempVolumeIntegralX, *ImposeBCsVolumeIntegralY, \
*ImposeBCsTempVolumeIntegralY, *ImposeBCsIEfm, *ImposeBCsIEfp, *ImposeBCsERHSx, *ImposeBCsERHSy, *ImposeBCsTempFacialIntegral, \
*ImposeBCsBotBEU, *ImposeBCsBotBEV, *ImposeBCsBotBEH, *ImposeBCsBotBEPSPX, *ImposeBCsBotBEPSPY, *ImposeBCsBotBEPWPS, *ImposeBCsCombinedTerms,\
*ImposeBCsWs;

extern char *ImposeBoundaryInitialized;

void MyExit()
{
	if (!strcmp("True", ImposeBoundaryInitialized)){
		SWENH3dImposeBoundaryMemoryDeAllocation();
		ImposeBoundaryInitialized = "False";
	}
	return;
}

void GetInverseSquareHeight(double *, double *, double , int );

void GetPenaltyParameter(double *, double , double , int , int , int );


void ImposeNewmannBoundaryCondition(double *, double *, mwIndex *, mwIndex *, int, \
	double *, double *, double *, double *, double *, double *, int, int, \
	double *, double *, double *, int, double *, double *);

void ImposeDirichletBoundaryCondition(double *, mwIndex *, mwIndex *, int, \
	int, int, double *, double *, double *, double *, double *, double *, double *, double *, \
    double *, double *, double *, double *, double *, double *, \
	int, double *, double *, double *, double *, double);

void GetScalarCentralFluxTerm2d(double *, double *, double *, double *, int );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexAtExit(&MyExit);
	double *StiffMatrix = mxGetPr(prhs[0]);
	mwIndex *Tempjcs = mxGetJc(prhs[0]);
	mwIndex *Tempirs = mxGetIr(prhs[0]);
	int row, col;
	row = (int)mxGetM(prhs[0]);
	col = (int)mxGetN(prhs[0]);
	double *sr;
	mwIndex *irs, *jcs;
	plhs[0] = mxCreateSparse(row, col, Tempjcs[col], mxREAL);
	sr = mxGetPr(plhs[0]);
	irs = mxGetIr(plhs[0]);
	jcs = mxGetJc(plhs[0]);
	memcpy(sr, StiffMatrix, Tempjcs[col] * sizeof(double));
	memcpy(irs, Tempirs, Tempjcs[col] * sizeof(mwIndex));
	memcpy(jcs, Tempjcs, (col + 1)*sizeof(mwIndex));

	double *K13 = mxGetPr(prhs[1]);
	double *K23 = mxGetPr(prhs[2]);
	double Hcrit = mxGetScalar(prhs[3]);
	double *Height = mxGetPr(prhs[4]);

	const mxArray *mesh = prhs[5];
	const mxArray *cell = prhs[6];
	const mxArray *BottomBoundaryEdge = prhs[7];
	const mxArray *BoundaryEdge = prhs[8];
	signed char *ftype = (signed char *)mxGetData(prhs[9]);
	const mxArray *mesh2d = prhs[10];
	const mxArray *cell2d = prhs[11];
	const mxArray *InnerEdge2d = prhs[12];
	const mxArray *BoundaryEdge2d = prhs[13];
	double *Wold = mxGetPr(prhs[14]);
	double *Wnew = mxGetPr(prhs[15]);
	double deltatime = mxGetScalar(prhs[16]);
	double rho = mxGetScalar(prhs[17]);
	double *Hu = mxGetPr(prhs[18]);
	double *Hv = mxGetPr(prhs[19]);
	double *RHS = mxGetPr(prhs[20]);
	double *PWPS = mxGetPr(prhs[21]);

	mxArray *TempFmask = mxGetField(cell, 0, "Fmask");
	double *Fmask = mxGetPr(TempFmask);
	int maxNfp = (int)mxGetM(TempFmask);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempMass3d = mxGetField(cell, 0, "M");
	double  *Mass3d = mxGetPr(TempMass3d);
	mxArray *TempDr = mxGetField(cell, 0, "Dr");
	double  *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(cell, 0, "Ds");
	double  *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double  *Dt = mxGetPr(TempDt);
	mxArray *TempN = mxGetField(cell, 0, "N");
	int N = (int)mxGetScalar(TempN);
	mxArray *TempNz = mxGetField(cell, 0, "Nz");
	int Nz = (int)mxGetScalar(TempNz);

	mxArray *TempNlayer = mxGetField(mesh, 0, "Nz");
	int Nlayer = (int)mxGetScalar(TempNlayer);
	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);
	mxArray *TempEToE = mxGetField(mesh, 0, "EToE");
	double *EToE = mxGetPr(TempEToE);
	mxArray *TempJ = mxGetField(mesh, 0, "J");
	double *J = mxGetPr(TempJ);
	mxArray *TempLAV = mxGetField(mesh, 0, "LAV");
	double *LAV = mxGetPr(TempLAV);
	mxArray *Temprx = mxGetField(mesh, 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(mesh, 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(mesh, 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *Tempsy = mxGetField(mesh, 0, "sy");
	double *sy = mxGetPr(Tempsy);
//	mxArray *Temprz = mxGetField(mesh, 0, "rz");
//	double *rz = mxGetPr(Temprz);
	mxArray *Temptz = mxGetField(mesh, 0, "tz");
	double *tz = mxGetPr(Temptz);

	plhs[1] = mxCreateDoubleMatrix(Np*K, 1, mxREAL);
	double *OutRHS = mxGetPr(plhs[1]);
	memcpy(OutRHS, RHS, Np*K*sizeof(double));

	mxArray *TempBEJs = mxGetField(BoundaryEdge, 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	int BENfp = (int)mxGetM(TempBEJs);
	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBEFToF = mxGetField(BoundaryEdge, 0, "FToF");
	double *BEFToF = mxGetPr(TempBEFToF);
	mxArray *TempBEFToE = mxGetField(BoundaryEdge, 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToN1 = mxGetField(BoundaryEdge, 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);
	mxArray *TempBEFToN2 = mxGetField(BoundaryEdge, 0, "FToN2");
	double *BEFToN2 = mxGetPr(TempBEFToN2);
	mxArray *TempBEnx = mxGetField(BoundaryEdge, 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(BoundaryEdge, 0, "ny");
	double *BEny = mxGetPr(TempBEny);
	mxArray *TempBEnz = mxGetField(BoundaryEdge, 0, "nz");
	double *BEnz = mxGetPr(TempBEnz);
	mxArray *TempBELMass2d = mxGetField(BoundaryEdge, 0, "M");
	double  *BELMass2d = mxGetPr(TempBELMass2d);
	mxArray *TempBELAV = mxGetField(BoundaryEdge, 0, "LAV");
	double *BELAV = mxGetPr(TempBELAV);

	mxArray *TempBotBEJs = mxGetField(BottomBoundaryEdge, 0, "Js");
	double *BotBEJs = mxGetPr(TempBotBEJs);
	int BotBENfp = (int)mxGetM(TempBotBEJs);
	mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBEFToF = mxGetField(BottomBoundaryEdge, 0, "FToF");
	double *BotBEFToF = mxGetPr(TempBotBEFToF);
	mxArray *TempBotBEFToE = mxGetField(BottomBoundaryEdge, 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToN1 = mxGetField(BottomBoundaryEdge, 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);
//	mxArray *TempBotBEnx = mxGetField(BottomBoundaryEdge, 0, "nx");
//	double *BotBEnx = mxGetPr(TempBotBEnx);
//	mxArray *TempBotBEny = mxGetField(BottomBoundaryEdge, 0, "ny");
//	double *BotBEny = mxGetPr(TempBotBEny);
//	mxArray *TempBotBEnz = mxGetField(BottomBoundaryEdge, 0, "nz");
//	double *BotBEnz = mxGetPr(TempBotBEnz);
	mxArray *TempBotBELMass2d = mxGetField(BottomBoundaryEdge, 0, "M");
	double  *BotBELMass2d = mxGetPr(TempBotBELMass2d);
	mxArray *TempBotBEFLAV = mxGetField(BottomBoundaryEdge, 0, "LAV");
	double *BotBEFLAV = mxGetPr(TempBotBEFLAV);

	mxArray *TempNp2d = mxGetField(cell2d, 0, "Np");
	int Np2d = (int)mxGetScalar(TempNp2d);
	mxArray *TempDr2d = mxGetField(cell2d, 0, "Dr");
	double *Dr2d = mxGetPr(TempDr2d);
	mxArray *TempDs2d = mxGetField(cell2d, 0, "Ds");
	double *Ds2d = mxGetPr(TempDs2d);
	mxArray *TempNface2d = mxGetField(cell2d, 0, "Nface");
	int Nface2d = (int)mxGetScalar(TempNface2d);
	mxArray *TempinvM2d = mxGetField(cell2d, 0, "invM");
	double *invM2d = mxGetPr(TempinvM2d);

	mxArray *TempK2d = mxGetField(mesh2d, 0, "K");
	int K2d = (int)mxGetScalar(TempK2d);
	mxArray *Temprx2d = mxGetField(mesh2d, 0, "rx");
	double *rx2d = mxGetPr(Temprx2d);
	mxArray *Tempsx2d = mxGetField(mesh2d, 0, "sx");
	double *sx2d = mxGetPr(Tempsx2d);
	mxArray *Tempry2d = mxGetField(mesh2d, 0, "ry");
	double *ry2d = mxGetPr(Tempry2d);
	mxArray *Tempsy2d = mxGetField(mesh2d, 0, "sy");
	double *sy2d = mxGetPr(Tempsy2d);
	mxArray *TempJ2d = mxGetField(mesh2d, 0, "J");
	double *J2d = mxGetPr(TempJ2d);

	mxArray *TempIENe2d = mxGetField(InnerEdge2d, 0, "Ne");
	int IENe2d = (int)mxGetScalar(TempIENe2d);
	mxArray *TempIENfp2d = mxGetField(InnerEdge2d, 0, "Nfp");
	int IENfp2d = (int)mxGetScalar(TempIENfp2d);
	mxArray *TempIEnx2d = mxGetField(InnerEdge2d, 0, "nx");
	double *IEnx2d = mxGetPr(TempIEnx2d);
	mxArray *TempIEny2d = mxGetField(InnerEdge2d, 0, "ny");
	double *IEny2d = mxGetPr(TempIEny2d);
	mxArray *TempIEFToN12d = mxGetField(InnerEdge2d, 0, "FToN1");
	double* IEFToN12d = mxGetPr(TempIEFToN12d);
	mxArray *TempIEFToN22d = mxGetField(InnerEdge2d, 0, "FToN2");
	double* IEFToN22d = mxGetPr(TempIEFToN22d);
	mxArray *TempIEMb2d = mxGetField(InnerEdge2d, 0, "M");
	double* IEMb2d = mxGetPr(TempIEMb2d);
	mxArray *TempIEJs2d = mxGetField(InnerEdge2d, 0, "Js");
	double *IEJs2d = mxGetPr(TempIEJs2d);
	mxArray *TempIEFToE2d = mxGetField(InnerEdge2d, 0, "FToE");
	double *IEFToE2d = mxGetPr(TempIEFToE2d);
	mxArray *TempIEFToF2d = mxGetField(InnerEdge2d, 0, "FToF");
	double *IEFToF2d = mxGetPr(TempIEFToF2d);

	if (!strcmp("False", ImposeBoundaryInitialized)){
		SWENH3dImposeBoundaryMemoryAllocation(Np, K, BENe, Np2d, K2d, Nface2d, IENfp2d, IENe2d, BotBENe, BotBENfp);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		GetInverseSquareHeight(ImposeBCsInvSquaHeight + Np*k, Height + Np*k, Hcrit, Np);
		for (int p = 0; p < Np; p++)
			ImposeBCsK33[k*Np + p] = pow(K13[k*Np + p], 2) + pow(K23[k*Np + p], 2) + \
			ImposeBCsInvSquaHeight[k*Np + p];
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int edge = 0; edge < BENe; edge++){
		if ((NdgEdgeType)ftype[edge] == NdgEdgeSlipWall){
			//Newmann boundary Doing Nothing
		}
		else{
			/*
			double *FpIndex = malloc(BENfp*sizeof(double));

			for (int p = 0; p < BENfp; p++){
				FpIndex[p] = BEFToN1[BENfp*edge + p];
			}

			CalculatePenaltyParameter(BETau, BEFToE, BEFToN1, BEFToN2, Np, BENfp, \
				edge, K13, K23, ImposeBCsK33, BELAV, LAV, max(N, Nz), Nface);
			

			int LocalEle;
			LocalEle = (int)BEFToE[2 * edge];

			double *TempEToE , *TempJ = NULL, *TempJs = NULL;
			TempEToE = EToE + (LocalEle - 1)*Nface;
			TempJ = J + (LocalEle - 1)*Np;
			TempJs = BEJs + edge * BENfp;

			ImposeDirichletBoundaryCondition(sr, irs, jcs, LocalEle, \
				Np, BENfp, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, \
				sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
				Dr, Ds, Dt, BEnx + edge * BENfp, BEny + edge * BENfp, BEnz + edge*BENfp, \
				TempJs, BELMass2d, TempEToE, Nface, FpIndex, \
				K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, ImposeBCsK33 + (LocalEle - 1)*Np, *(BETau + edge));

			free(FpIndex);
			*/
		}
	}
	/*The following part is used to calculate the Neumann boundary condition according to 
	$\frac{\partial w}{\partial t}+u\frac{\partial w}{\partial x}+v\frac{\partial w}{\partial y} + w\frac{\partial w}{\partial z}=-\frac{1}{\rho D} \frac{\partial p}{\partial \sigma}$*/
	ptrdiff_t np = Np2d;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
		GetVolumnIntegral2d(ImposeBCsVolumeIntegralX + k*Np2d, ImposeBCsTempVolumeIntegralX + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, Wnew + k*Np2d, &np, &zero, &np, rx2d + k*Np2d, sx2d + k*Np2d);
		/*$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$*/
		GetVolumnIntegral2d(ImposeBCsVolumeIntegralY + k*Np2d, ImposeBCsTempVolumeIntegralY + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, Wnew + k*Np2d, &np, &zero, &np, ry2d + k*Np2d, sy2d + k*Np2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		FetchInnerEdgeFacialValue(ImposeBCsIEfm + e*IENfp2d, ImposeBCsIEfp + e*IENfp2d, Wnew, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		GetScalarFluxTerm2d(ImposeBCsWIEFluxMx2d + e*IENfp2d, ImposeBCsIEfm + e*IENfp2d, IEnx2d + e*IENfp2d, IENfp2d);
		GetScalarFluxTerm2d(ImposeBCsWIEFluxPx2d + e*IENfp2d, ImposeBCsIEfp + e*IENfp2d, IEnx2d + e*IENfp2d, IENfp2d);
		GetScalarFluxTerm2d(ImposeBCsWIEFluxMy2d + e*IENfp2d, ImposeBCsIEfm + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetScalarFluxTerm2d(ImposeBCsWIEFluxPy2d + e*IENfp2d, ImposeBCsIEfp + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetScalarCentralFluxTerm2d(ImposeBCsWIEFluxSx2d + e*IENfp2d, ImposeBCsIEfm + e*IENfp2d, ImposeBCsIEfp + e*IENfp2d, IEnx2d + e*IENfp2d, IENfp2d);
		GetScalarCentralFluxTerm2d(ImposeBCsWIEFluxSy2d + e*IENfp2d, ImposeBCsIEfm + e*IENfp2d, ImposeBCsIEfp + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
	}

	memset(ImposeBCsERHSx, 0, Np2d*K2d*Nface2d*sizeof(double));
	memset(ImposeBCsERHSy, 0, Np2d*K2d*Nface2d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, ImposeBCsWIEFluxMx2d, ImposeBCsWIEFluxPx2d, ImposeBCsWIEFluxSx2d, IEJs2d, IEMb2d, ImposeBCsERHSx);
		StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, ImposeBCsWIEFluxMy2d, ImposeBCsWIEFluxPy2d, ImposeBCsWIEFluxSy2d, IEJs2d, IEMb2d, ImposeBCsERHSy);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		for (int face = 1; face<Nface2d; face++){
			Add(ImposeBCsERHSx + k*Np2d, ImposeBCsERHSx + k*Np2d, ImposeBCsERHSx + face*Np2d*K2d + k*Np2d, Np2d);
			Add(ImposeBCsERHSy + k*Np2d, ImposeBCsERHSy + k*Np2d, ImposeBCsERHSy + face*Np2d*K2d + k*Np2d, Np2d);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		MultiEdgeContributionByLiftOperator(ImposeBCsERHSx + k*Np2d, ImposeBCsTempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
		MultiEdgeContributionByLiftOperator(ImposeBCsERHSy + k*Np2d, ImposeBCsTempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		Minus(ImposeBCsWx + k*Np2d, ImposeBCsVolumeIntegralX + k*Np2d, ImposeBCsERHSx + k*Np2d, Np2d);
		Minus(ImposeBCsWy + k*Np2d, ImposeBCsVolumeIntegralY + k*Np2d, ImposeBCsERHSy + k*Np2d, Np2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		/*$\frac{\partial w}{\partial \sigma}$*/
		FetchBoundaryEdgeFacialValue(ImposeBCsBotBEPWPS + face*BotBENfp, PWPS, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*$\frac{\partial \sigma}{\partial x}$*/
		FetchBoundaryEdgeFacialValue(ImposeBCsBotBEPSPX + face*BotBENfp, K13, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*$\frac{\partial \sigma}{\partial y}$*/
		FetchBoundaryEdgeFacialValue(ImposeBCsBotBEPSPY + face*BotBENfp, K23, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*Hu*/
		FetchBoundaryEdgeFacialValue(ImposeBCsBotBEU + face*BotBENfp, Hu, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*Hv*/
		FetchBoundaryEdgeFacialValue(ImposeBCsBotBEV + face*BotBENfp, Hv, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*H*/
		FetchBoundaryEdgeFacialValue(ImposeBCsBotBEH + face*BotBENfp, Height, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*$u=\frac{Hu}{H}$*/
		DotCriticalDivide(ImposeBCsBotBEU + face*BotBENfp, ImposeBCsBotBEU + face*BotBENfp, &Hcrit, ImposeBCsBotBEH + face*BotBENfp, BotBENfp);
		/*$v=\frac{Hv}{H}$*/
		DotCriticalDivide(ImposeBCsBotBEV + face*BotBENfp, ImposeBCsBotBEV + face*BotBENfp, &Hcrit, ImposeBCsBotBEH + face*BotBENfp, BotBENfp);
		/*\frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial x}*/
		DotProduct(ImposeBCsCombinedTerms + face*BotBENfp, ImposeBCsBotBEPWPS + face*BotBENfp, ImposeBCsBotBEPSPX + face*BotBENfp, BotBENfp);
		/*\frac{\partial w}{\partial x} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial x}*/
		Add(ImposeBCsWx + face*BotBENfp, ImposeBCsWx + face*BotBENfp, ImposeBCsCombinedTerms + face*BotBENfp, BotBENfp);
		/*$u\left (\frac{\partial w}{\partial x} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial x}\right )$*/
		DotProduct(ImposeBCsWx + face*BotBENfp, ImposeBCsWx + face*BotBENfp, ImposeBCsBotBEU + face*BotBENfp, BotBENfp);

		/*\frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial y}*/
		DotProduct(ImposeBCsCombinedTerms + face*BotBENfp, ImposeBCsBotBEPWPS + face*BotBENfp, ImposeBCsBotBEPSPY + face*BotBENfp, BotBENfp);
		/*\frac{\partial w}{\partial y} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial y}*/
		Add(ImposeBCsWy + face*BotBENfp, ImposeBCsWy + face*BotBENfp, ImposeBCsCombinedTerms + face*BotBENfp, BotBENfp);
		/*$v\left (\frac{\partial w}{\partial y} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial y}\right )$*/
		DotProduct(ImposeBCsWy + face*BotBENfp, ImposeBCsWy + face*BotBENfp, ImposeBCsBotBEV + face*BotBENfp, BotBENfp);

        /*$\frac{w}{D}$*/
		DotCriticalDivide(ImposeBCsWs + face*BotBENfp, Wnew + face*BotBENfp, &Hcrit, ImposeBCsBotBEH + face*BotBENfp, BotBENfp);
		/*$\frac{w}{D}\frac{\partial w}{\partial \sigma}$*/
		DotProduct(ImposeBCsWs + face*BotBENfp, ImposeBCsWs + face*BotBENfp, ImposeBCsBotBEPWPS + face*BotBENfp, BotBENfp);

		/*$w_{new} - w_{old}$*/
		Minus(ImposeBCsNewmannData + face*BotBENfp, Wnew + face*BotBENfp, Wold + face*BotBENfp, BotBENfp);
		/*$\frac{w_{new} - w_{old}}{\Delta t}$*/
		DotDivideByConstant(ImposeBCsNewmannData + face*BotBENfp, ImposeBCsNewmannData + face*BotBENfp, deltatime, BotBENfp);
		/*$\frac{w_{new} - w_{old}}{\Delta t} + u\left (\frac{\partial w}{\partial x} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial x}\right ) $*/
		Add(ImposeBCsNewmannData + face*BotBENfp, ImposeBCsNewmannData + face*BotBENfp, ImposeBCsWx + face*BotBENfp, BotBENfp);
		/*$\frac{w_{new} - w_{old}}{\Delta t} + u\left (\frac{\partial w}{\partial x} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial x}\right ) + v\left (\frac{\partial w}{\partial y} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial y}\right )$*/
		Add(ImposeBCsNewmannData + face*BotBENfp, ImposeBCsNewmannData + face*BotBENfp, ImposeBCsWy + face*BotBENfp, BotBENfp);

		/*$\frac{w_{new} - w_{old}}{\Delta t} + u\left (\frac{\partial w}{\partial x} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial x}\right ) + v\left (\frac{\partial w}{\partial y} + \frac{\partial w}{\partial \sigma}\frac{\partial \sigma}{\partial y}\right ) + \frac{w}{D}\frac{\partial w}{\partial \sigma}$*/
		Add(ImposeBCsNewmannData + face*BotBENfp, ImposeBCsNewmannData + face*BotBENfp, ImposeBCsWs + face*BotBENfp, BotBENfp);
		/*$\rho D\left (\frac{w_{new} - w_{old}}{\Delta t} + u\frac{\partial w}{\partial x} + v\frac{\partial w}{\partial y}\right )$*/
		MultiplyByConstant(ImposeBCsNewmannData + face*BotBENfp, ImposeBCsNewmannData + face*BotBENfp, rho, BotBENfp);
		/*$\frac{1}{\rho D}\left (\frac{w_{new} - w_{old}}{\Delta t} + u\frac{\partial w}{\partial x} + v\frac{\partial w}{\partial y}\right )$*/
		DotProduct(ImposeBCsNewmannData + face*BotBENfp, ImposeBCsNewmannData + face*BotBENfp, ImposeBCsBotBEH + face*BotBENfp, BotBENfp);
		/*$n_{\sigma}\rho D\left (\frac{w_{new} - w_{old}}{\Delta t} + u\frac{\partial w}{\partial x} + v\frac{\partial w}{\partial y}\right )$*/
		MultiplyByConstant(ImposeBCsNewmannData + face*BotBENfp, ImposeBCsNewmannData + face*BotBENfp, -1.0, BotBENfp);
		/*$\left (\frac{\partial \sigma}{\partial x}\right )^2$*/
		DotProduct(ImposeBCsBotBEPSPX + face*BotBENfp, ImposeBCsBotBEPSPX + face*BotBENfp, ImposeBCsBotBEPSPX + face*BotBENfp, BotBENfp);
		/*$\left (\frac{\partial \sigma}{\partial y}\right )^2$*/
		DotProduct(ImposeBCsBotBEPSPY + face*BotBENfp, ImposeBCsBotBEPSPY + face*BotBENfp, ImposeBCsBotBEPSPY + face*BotBENfp, BotBENfp);
		/*$\frac{1}{D}$*/
		ReverseValue(ImposeBCsBotBEH + face*BotBENfp, ImposeBCsBotBEH + face*BotBENfp, BotBENfp);
		/*$\frac{1}{D^2}$*/
		DotProduct(ImposeBCsBotBEH + face*BotBENfp, ImposeBCsBotBEH + face*BotBENfp, ImposeBCsBotBEH + face*BotBENfp, BotBENfp);
		/*$\left (\frac{\partial \sigma}{\partial x}\right )^2 + \left (\frac{\partial \sigma}{\partial y}\right )^2$*/
		Add(ImposeBCsBotBEPSPX + face*BotBENfp, ImposeBCsBotBEPSPX + face*BotBENfp, ImposeBCsBotBEPSPY + face*BotBENfp, BotBENfp);
		/*$\left (\frac{\partial \sigma}{\partial x}\right )^2 + \left (\frac{\partial \sigma}{\partial y}\right )^2 + \frac{1}{D^2}$*/
		Add(ImposeBCsBotBEPSPX + face*BotBENfp, ImposeBCsBotBEPSPX + face*BotBENfp, ImposeBCsBotBEH + face*BotBENfp, BotBENfp);

		DotProduct(ImposeBCsNewmannData + face*BotBENfp, ImposeBCsNewmannData + face*BotBENfp, ImposeBCsBotBEPSPX + face*BotBENfp, BotBENfp);
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int edge = 0; edge < BotBENe; edge++){

		int LocalEle;
		LocalEle = (int)BotBEFToE[2 * edge];
		double *DxBuff = malloc(Np*Np*sizeof(double));
		DiagMultiply(DxBuff, Dr, rx + (LocalEle - 1)*Np, Np);
		double *Dx = malloc(Np*Np*sizeof(double));
		DiagMultiply(Dx, Ds, sx + (LocalEle - 1)*Np, Np);
		Add(Dx, Dx, DxBuff, Np*Np);

		double *Dz = malloc(Np*Np*sizeof(double));
		DiagMultiply(Dz, Dt, tz + (LocalEle - 1)*Np, Np);

		double *DyBuff = malloc(Np*Np*sizeof(double));
		DiagMultiply(DyBuff, Dr, ry + (LocalEle - 1)*Np, Np);
		double *Dy = malloc(Np*Np*sizeof(double));
		DiagMultiply(Dy, Ds, sy + (LocalEle - 1)*Np, Np);
		Add(Dy, Dy, DyBuff, Np*Np);

		double *TempEToE = NULL, *TempJs = NULL;
		TempEToE = EToE + (LocalEle - 1)*Nface;
		TempJs = BotBEJs + edge * BotBENfp;

		double *FpIndex = malloc(BotBENfp*sizeof(double));

		for (int p = 0; p < BotBENfp; p++){
			FpIndex[p] = BotBEFToN1[BotBENfp*edge + p];
		}

		ImposeNewmannBoundaryCondition( OutRHS, sr, irs, jcs, LocalEle, \
			K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, ImposeBCsK33 + (LocalEle - 1)*Np, \
			Dx, Dy, Dz, Np, BotBENfp, \
			TempJs, BotBELMass2d, TempEToE, Nface, FpIndex, ImposeBCsNewmannData + edge * BotBENfp);

		free(DxBuff);
		free(Dx);
		free(Dz);
		free(DyBuff);
		free(Dy);
		free(FpIndex);
	}
}

void GetInverseSquareHeight(double *dest, double *source, double Hcrit, int Np){
	for (int i = 0; i < Np; i++){
		if (source[i] >= Hcrit){
			dest[i] = 1.0 / source[i] / source[i];
		}
		else{
			dest[i] = 0;
		}
	}
}

void ImposeNewmannBoundaryCondition(double *RHSdest, double *dest, mwIndex *Irs, mwIndex *Jcs, int LocalEle, \
	double *K31, double *K32, double *K33, double *Dx, double *Dy, double *Dz, int Np, int Nfp,\
	double *Js, double *M2d, double *EToE, int Nface, double *FpIndex, double *NewmannData){

	double *TempK31 = malloc(Np*Np*sizeof(double));
	double *TempK32 = malloc(Np*Np*sizeof(double));
	double *TempK33 = malloc(Np*Np*sizeof(double));
	double *TempCoe = malloc(Np*Np*sizeof(double));

	DiagMultiply(TempK31, Dx, K31, Np);

	DiagMultiply(TempK32, Dy, K32, Np);

	DiagMultiply(TempK33, Dz, K33, Np);

	Add(TempCoe, TempK31, TempK32, Np*Np);

	double *TempRHSBuff = malloc(Np * 1 * sizeof(double));

	memset(TempRHSBuff, 0, Np * 1 * sizeof(double));

	double *TempRHSFacialData = malloc(Nfp * 1 * sizeof(double));
//	Add(TempCoe, TempCoe, TempK33, Np*Np);

	/*Withdraw the data in $\frac{\partial p}{\partial \sigma}$, and store them in TempPNPS*/
	double *TempFacialData = malloc(Nfp*sizeof(double));
	double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(EleMass2d, M2d, Js, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	int UniNum = 0, StartPoint;

	/*Find the exact place where to fill in the data, and the place is stored in StartPoint*/
	double *TempEToE = malloc((Nface + 1)*sizeof(double));

	FindUniqueElementAndSortOrder(TempEToE, EToE, &UniNum, Nface, LocalEle);

	int NonzeroPerColumn = Jcs[(LocalEle - 1)*Np + 1] - Jcs[(LocalEle - 1)*Np];

	for (int j = 0; j < UniNum; j++){
		if ((int)TempEToE[j] == LocalEle)
			StartPoint = (int)Jcs[(LocalEle - 1)*Np] + j*Np;
	}
	double *ContributionPerPoint = malloc(Nfp*sizeof(double));

	ptrdiff_t One = 1;
	/*
	for (int j = 0; j < Nfp; j++){

		FetchFacialData(TempFacialData, TempCoe + ((int)FpIndex[j] - 1)*Np, FpIndex, Nfp);

		MatrixMultiply("N", "N", (ptrdiff_t)Nfp, One, (ptrdiff_t)Nfp, 1.0, EleMass2d,
			(ptrdiff_t)Nfp, TempFacialData, (ptrdiff_t)Nfp, 0.0, ContributionPerPoint, (ptrdiff_t)Nfp);

		MultiplyByConstant(ContributionPerPoint, ContributionPerPoint, -1.0, Nfp);

		AssembleDataIntoPoint(TempContribution + ((int)FpIndex[j] - 1)*Np, ContributionPerPoint, FpIndex, Nfp);
	}
	*/
	
	
	for (int j = 0; j < Np; j++){

		FetchFacialData(TempFacialData, TempCoe + j*Np, FpIndex, Nfp);

		MatrixMultiply("N", "N", (ptrdiff_t)Nfp, One, (ptrdiff_t)Nfp, 1.0, EleMass2d,
			(ptrdiff_t)Nfp, TempFacialData, (ptrdiff_t)Nfp, 0.0, ContributionPerPoint, (ptrdiff_t)Nfp);

		MultiplyByConstant(ContributionPerPoint, ContributionPerPoint, -1.0, Nfp);

		AssembleDataIntoPoint(TempContribution + j*Np, ContributionPerPoint, FpIndex, Nfp);
	}
	
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, One, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, NewmannData, (ptrdiff_t)Nfp, 0.0, TempRHSFacialData, (ptrdiff_t)Nfp);

	AssembleDataIntoPoint(TempRHSBuff, TempRHSFacialData, FpIndex, Nfp);

	MultiplyByConstant(TempRHSBuff, TempRHSBuff, -1.0, Np);

	Add(RHSdest + (LocalEle - 1)*Np, RHSdest + (LocalEle - 1)*Np, TempRHSBuff, Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempContribution, NonzeroPerColumn, Np);

	free(TempK31);
	free(TempK32);
	free(TempK33);
	free(TempCoe);
	free(TempFacialData);
	free(EleMass2d);
	free(TempContribution);
	free(TempEToE);
	free(ContributionPerPoint);
	free(TempRHSBuff);
	free(TempRHSFacialData);
}

void ImposeDirichletBoundaryCondition(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, \
	int Np, int Nfp, double *rx, double *sx, double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, \
	double *nx, double *ny, double *nz, double *Js, double *Mass2d, double *EToE, \
	int Nface, double *FpIndex, double *K13, double *K23, double *K33, double Tau){
	double *DxBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DxBuff, Dr, rx, Np);
	double *Dx = malloc(Np*Np*sizeof(double));
	DiagMultiply(Dx, Ds, sx, Np);
	Add(Dx, Dx, DxBuff, Np*Np);

	double *Dz = malloc(Np*Np*sizeof(double));
	DiagMultiply(Dz, Dt, tz, Np);

	double *DyBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DyBuff, Dr, ry, Np);
	double *Dy = malloc(Np*Np*sizeof(double));
	DiagMultiply(Dy, Ds, sy, Np);
	Add(Dy, Dy, DyBuff, Np*Np);

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(EleMass2d, Mass2d, Js, Nfp);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));
	double *EdgeContribution = malloc(Np*Nfp*sizeof(double));

	/*For fifth term part*/
	/* For term $\int_{\partial \Omega^D}u_h\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	/*For $k_{11}\frac{\partial v}{\partial x}n_xp$*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, Dx, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, nx[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/* For term $k_{11}\frac{\partial p}{\partial x}n_xv$, x direction first*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/*For $k_{13}\frac{\partial v}{\partial \sigma}n_xp$*/

	DiagMultiply(TempDiffMatrix, Dz, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, TempDiffMatrix, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, nx[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/* For term $k_{13}\frac{\partial p}{\partial \sigma}n_xv$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/*For $k_{22}\frac{\partial v}{\partial y}n_yp$*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, Dy, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, ny[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/* For term $k_{22}\frac{\partial p}{\partial y}n_yv$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/*For $k_{23}\frac{\partial v}{\partial \sigma}n_yp$*/
	DiagMultiply(TempDiffMatrix, Dz, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, TempDiffMatrix, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, ny[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/* For term $k_{23}\frac{\partial p}{\partial \sigma}n_yv$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/*For $k_{31}\frac{\partial v}{\partial x}n_{\sigma}p$*/
	DiagMultiply(TempDiffMatrix, Dx, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, TempDiffMatrix, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, nz[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/* For term $k_{31}\frac{\partial p}{\partial x}n_{\sigma}v$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/*For $k_{32}\frac{\partial v}{\partial y}n_{\sigma}p$*/
	DiagMultiply(TempDiffMatrix, Dy, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, TempDiffMatrix, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, nz[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/* For term $k_{32}\frac{\partial p}{\partial y}n_{\sigma}v$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/*For $k_{33}\frac{\partial v}{\partial \sigma}n_{\sigma}p$*/
	DiagMultiply(TempDiffMatrix, Dz, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, TempDiffMatrix, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, nz[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/* For term $k_{33}\frac{\partial p}{\partial \sigma}n_{\sigma}v$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));
	MultiplyByConstant(TempMass2d, EleMass2d, Tau, Nfp*Nfp);
	/*For term $-\int_{\partial \Omega^d}\tau^k s u_hd\boldsymbol{x}$*/
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, FpIndex, FpIndex, Np, Nfp, -1.0);

	int UniNum = 0, StartPoint;
	/*Find the exact place where to fill in the data, and the place is stored in StartPoint*/
	double *TempEToE = malloc((Nface + 1)*sizeof(double));
	FindUniqueElementAndSortOrder(TempEToE, EToE, &UniNum, Nface, LocalEle);
	int NonzeroPerColumn = jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np];
	for (int j = 0; j < UniNum; j++){
		if ((int)TempEToE[j] == LocalEle){
			StartPoint = jcs[(LocalEle - 1)*Np] + j*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, TempContribution, NonzeroPerColumn, Np);
			break;
		}
	}

	free(DxBuff);
	free(Dx);
	free(DyBuff);
	free(Dy);
	free(Dz);
	free(EleMass2d);
	free(TempContribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempEToE);
	free(TempMass2d);
	free(TempDiffMatrix);
}

void GetScalarCentralFluxTerm2d(double *dest, double *fm, double *fp, double *Vector, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = 0.5*Vector[i] * (fm[i] + fp[i]);
}