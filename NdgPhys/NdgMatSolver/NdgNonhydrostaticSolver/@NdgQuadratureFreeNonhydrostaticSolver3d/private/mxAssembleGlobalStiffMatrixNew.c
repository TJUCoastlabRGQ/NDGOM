#include "SWENonhydrostatic3d.h"

extern double *K33, *InvSquaHeight, *InnerEdgeTau, *BottomEdgeTau, *SurfaceEdgeTau;
extern double *IEK33M, *IEK33P;
extern double *BotEK33M, *BotEK33P;
extern double *SurfEK33M;
extern double *IEWeight1, *IEWeight2;
extern double *BotEWeight1, *BotEWeight2;
extern double *BEK33M;
extern double *BEWeight1, *BEWeight2;
extern double *BoundaryEdgeTau;
/*
 * Actually, SurfWeight1 and SurfWeight2 are not needed, 
 * we used it here just to accommodate
 * the usage of function CalculatePenaltyParameterAndWeight
 */ 
extern double *SurfEWeight1, *SurfEWeight2;
extern int *NonIr, *NonJc;
extern char *GlobalStiffMatrixInitialized;

void MyExit()
{
	if (!strcmp("True", GlobalStiffMatrixInitialized)){
		GlobalStiffMatrixMemoryDeAllocation();
		GlobalStiffMatrixInitialized = "False";
	}
	return;
}

void GetInverseSquareHeight(double *, double *, double, int); 

void GetLocalVolumnIntegralTerm(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, double *M3d, \
	double *Dr, double *Ds, double *Dt, double *rx, double *sx, double *ry, double *sy, double *tz, \
	double *J, double *K13, double *K23, double *K33, int Localele);

void GetLocalToDownFacialContribution(double *dest, mwIndex *Irs, mwIndex *jcs, int LocalEle, int DownEle, int Np, \
	int Np2d, double *M2d, double *J2d, double *LocalEid, double *DownEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, \
	double *Adjry, double *Adjsy, double *Adjtz, double *Dr, double *Ds, double *Dt, \
	double *LocalK13, double *LocalK23, double *LocalK33, \
	double *AdjK13, double *AdjK23, double *AdjK33, \
	double Tau, int face);

void GetLocalToUpFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int UpEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *UpEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33, double Tau, int face);

void GetLocalDownFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, \
	double *K33, double Tau, int face);

void GetLocalUpFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33, \
	double Tau, int face);

void ImposeDirichletBoundaryCondition(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, \
	double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, \
	double *K23, double *K33, double Tau, int LocalEle);

void GetLocalToAdjacentFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *M2d, double *J2d, double *LocalEid, double *AdjEid, double *Dr, double *Ds, double *Dt, \
	double *LocalRx, double *LocalSx, double *LocalRy, double *LocalSy, double *LocalTz, \
	double *AdjRx, double *AdjSx, double *AdjRy, double *AdjSy, double *AdjTz, double *nx, double *ny, \
	double *LocalK13, double *LocalK23, double *AdjK13, double *AdjK23, double Tau, int LocalEle, \
	int AdjEle, int Face, int Flag);

void GetLocalFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *M2d, double *J2d, double *LocalEid, double *Dr, double *Ds, double *Dt, double *rx, double *sx, \
	double *ry, double *sy, double *tz, double *nx, double *ny, double *K13, double *K23, double Tau, \
	int LocalEle, int Face, int Flag);

void ImposeHorizontalDirichletBoundaryCondition(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, \
	double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, \
	double *K13, double *K23, double Tau, int LocalEle, double *nx, double *ny);

/*This function is used to assemble the global stiff matrix used in the 
three-dimensional non-hydrostatic solver, the input parameters are organized as follows:
PNPS: $\frac{\partial q}{\partial \sigma}$, indexed as 0;
SPNPX: $\frac{\partial^2 q}{\partial x^2}$, indexed as 1;
SPNPY: $\frac{\partial^2 q}{\partial y^2}$, indexed as 2;
SPNPS: $\frac{\partial^2 q}{\partial \sigma^2}$, indexed as 3;
MSPNPX: $\frac{\partial}{\partial x}\left (\frac{\partial q}{\partial \sigma}\right )$, indexed as 4;
MSPNPY: $\frac{\partial}{\partial y}\left (\frac{\partial q}{\partial \sigma}\right )$, indexed as 5;
PSPX: $\frac{\partial \sigma}{\partial x}$, indexed as 6;
PSPY: $\frac{\partial \sigma}{\partial y}$, indexed as 7;
SPSPX: $\frac{\partial^2 \sigma}{\partial x^2}$, indexed as 8;
SPSPY: $\frac{\partial^2 \sigma}{\partial y^2}$, indexed as 9;
SQPSPX: $\left (\frac{\partial \sigma}{\partial x}\right )^2$, indexed as 10;
SQPSPY: $\left (\frac{\partial \sigma}{\partial y}\right )^2$, indexed as 11;
Height: The water depth, indexed as 12;
Hcrit: The critical water depth, indexed as 13;
BoundaryEdge2d: The two-dimensional boundary edge, indexed as 14;
mesh: The three-dimensional mesh object, indexed as 15;
cell: The three-dimensional master cell, indexed as 16;
BoundaryEdge: The three-dimensional mesh object, indexed as 17;
ftype2d: The two-dimensional boundary edge type, indexed as 18;
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexAtExit(&MyExit);

	/*$\frac{\partial \sigma}{\partial x^*}$*/
	double *K13 = mxGetPr(prhs[0]);
	/*$\frac{\partial \sigma}{\partial y^*}$*/
	double *K23 = mxGetPr(prhs[1]);

	double *SQPSPX = mxGetPr(prhs[2]);
	double *SQPSPY = mxGetPr(prhs[3]);

	double Hcrit = mxGetScalar(prhs[4]);
	double *Height = mxGetPr(prhs[5]);

	const mxArray *mesh = prhs[6];

	mxArray *Temprx = mxGetField(mesh, 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(mesh, 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(mesh, 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *TempLAV = mxGetField(mesh, 0, "LAV");
	double *LAV = mxGetPr(TempLAV);
	mxArray *Tempsy = mxGetField(mesh, 0, "sy");
	double *sy = mxGetPr(Tempsy);
	mxArray *Temptz = mxGetField(mesh, 0, "tz");
	double *tz = mxGetPr(Temptz);
	mxArray *TempJ = mxGetField(mesh, 0, "J");
	double *J = mxGetPr(TempJ);
	mxArray *TempEToE = mxGetField(mesh, 0, "EToE");
	double *EToE = mxGetPr(TempEToE);
	mxArray *TempNlayer = mxGetField(mesh, 0, "Nz");
	int Nlayer = (int)mxGetScalar(TempNlayer);
	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);


	const mxArray *cell = prhs[7];
	mxArray *TempFmask = mxGetField(cell, 0, "Fmask");
	double *Fmask = mxGetPr(TempFmask);
	int maxNfp = (int)mxGetM(TempFmask);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	int Nface2d = Nface - 2;
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempN = mxGetField(cell, 0, "N");
	int N = (int)mxGetScalar(TempN);
	mxArray *TempNz = mxGetField(cell, 0, "Nz");
	int Nz = (int)mxGetScalar(TempNz);
	mxArray *TempMass3d = mxGetField(cell, 0, "M");
	double  *M3d = mxGetPr(TempMass3d);
	mxArray *TempDr = mxGetField(cell, 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(cell, 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempNfp = mxGetField(cell, 0, "Nfp");
	double *Nfp = mxGetPr(TempNfp);
	int Np2d = (int)Nfp[Nface - 1];
	double *UpEidM = malloc(Np2d*sizeof(double));
	double *BotEidM = malloc(Np2d*sizeof(double));
	for (int i = 0; i < Np2d; i++){
		UpEidM[i] = Fmask[(Nface - 1)*maxNfp + i];
		BotEidM[i] = Fmask[(Nface - 2)*maxNfp + i];
	}

	const mxArray *InnerEdge = prhs[8];
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
	mxArray *TempIEnz = mxGetField(InnerEdge, 0, "nz");
	double *IEnz = mxGetPr(TempIEnz);
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

	const mxArray *BottomEdge = prhs[9];
	mxArray *TempBotENe = mxGetField(BottomEdge, 0, "Ne");
	int BotENe = (int)mxGetScalar(TempBotENe);
	mxArray *TempBotENfp = mxGetField(BottomEdge, 0, "Nfp");
	int BotENfp = (int)mxGetScalar(TempBotENfp);
	mxArray *TempBotEnx = mxGetField(BottomEdge, 0, "nx");
	double *BotEnx = mxGetPr(TempBotEnx);
	mxArray *TempBotEny = mxGetField(BottomEdge, 0, "ny");
	double *BotEny = mxGetPr(TempBotEny);
	mxArray *TempBotEnz = mxGetField(BottomEdge, 0, "nz");
	double *BotEnz = mxGetPr(TempBotEnz);
	mxArray *TempBotEFToE = mxGetField(BottomEdge, 0, "FToE");
	double *BotEFToE = mxGetPr(TempBotEFToE);
	mxArray *TempBotEFToN1 = mxGetField(BottomEdge, 0, "FToN1");
	double *BotEFToN1 = mxGetPr(TempBotEFToN1);
	mxArray *TempBotEFToN2 = mxGetField(BottomEdge, 0, "FToN2");
	double *BotEFToN2 = mxGetPr(TempBotEFToN2);
	mxArray *TempBotELAV = mxGetField(BottomEdge, 0, "LAV");
	double *BotELAV = mxGetPr(TempBotELAV);

	const mxArray *SurfaceBoundaryEdge = prhs[10];
	mxArray *TempSurfBENe = mxGetField(SurfaceBoundaryEdge, 0, "Ne");
	int SurfBENe = (int)mxGetScalar(TempSurfBENe);
	mxArray *TempSurfBENfp = mxGetField(SurfaceBoundaryEdge, 0, "Nfp");
	int SurfBENfp = (int)mxGetScalar(TempSurfBENfp);
	mxArray *TempSurfBEnx = mxGetField(SurfaceBoundaryEdge, 0, "nx");
	double *SurfBEnx = mxGetPr(TempSurfBEnx);
	mxArray *TempSurfBEny = mxGetField(SurfaceBoundaryEdge, 0, "ny");
	double *SurfBEny = mxGetPr(TempSurfBEny);
	mxArray *TempSurfBEnz = mxGetField(SurfaceBoundaryEdge, 0, "nz");
	double *SurfBEnz = mxGetPr(TempSurfBEnz);
	mxArray *TempSurfBEFToE = mxGetField(SurfaceBoundaryEdge, 0, "FToE");
	double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
	mxArray *TempSurfBEFToN1 = mxGetField(SurfaceBoundaryEdge, 0, "FToN1");
	double *SurfBEFToN1 = mxGetPr(TempSurfBEFToN1);
	mxArray *TempSurfBEFToN2 = mxGetField(SurfaceBoundaryEdge, 0, "FToN2");
	double *SurfBEFToN2 = mxGetPr(TempSurfBEFToN2);
	mxArray *TempSurfBELAV = mxGetField(SurfaceBoundaryEdge, 0, "LAV");
	double *SurfBELAV = mxGetPr(TempSurfBELAV);


	double *M2d = mxGetPr(prhs[11]);
	double *J2d = mxGetPr(prhs[12]);
	int K2d = (int)mxGetScalar(prhs[13]);

	char* BoundaryType;
	BoundaryType = mxArrayToString(prhs[14]);

	double *SortedIEnx = mxGetPr(prhs[15]);
	double *SortedIEny = mxGetPr(prhs[16]);
	double *SortedIEGlobalFace = mxGetPr(prhs[17]);
	double *SortedIEAdjEle = mxGetPr(prhs[18]);
	double *SortedIEReverseFlag = mxGetPr(prhs[19]);
	double *SortedIEInternalFace = mxGetPr(prhs[20]);

	const mxArray *BoundaryEdge = prhs[21];
	signed char *ftype = (signed char *)mxGetData(prhs[22]);
	mxArray *TempBEMb = mxGetField(BoundaryEdge, 0, "M");
	double *BEMb = mxGetPr(TempBEMb);
	mxArray *TempBEJs = mxGetField(BoundaryEdge, 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	mxArray *TempBEFToE = mxGetField(BoundaryEdge, 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToN1 = mxGetField(BoundaryEdge, 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);
	mxArray *TempBENfp = mxGetField(BoundaryEdge, 0, "Nfp");
	int BENfp = (int)mxGetScalar(TempBENfp);
	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBEnx = mxGetField(BoundaryEdge, 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(BoundaryEdge, 0, "ny");
	double *BEny = mxGetPr(TempBEny);
	mxArray *TempBEnz = mxGetField(BoundaryEdge, 0, "nz");
	double *BEnz = mxGetPr(TempBEnz);
	mxArray *TempBELAV = mxGetField(BoundaryEdge, 0, "LAV");
	double *BELAV = mxGetPr(TempBELAV);

	int Nonzero = K * Np*Np + 2 * IENe * (2 * Np*IENfp - IENfp*IENfp) + \
		2 * BotENe * (2 * Np*BotENfp - BotENfp*BotENfp);

	if (!strcmp("False", GlobalStiffMatrixInitialized)){

		GlobalStiffMatrixMemoryAllocation(Np, K, IENe, BotENe, BENe, SurfBENe, BENfp, IENfp, BotENfp, Nonzero);

		GetSparsePattern(NonIr, NonJc, EToE, IEFToE, IEFToN1, IEFToN2, BotEFToE, \
			BotEFToN1, BotEFToN2, Nface, IENfp, BotENfp, Np, K, IENe, BotENe);

	}

	int i, k, p, ele, L, face;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (k = 0; k < K; k++){
		GetInverseSquareHeight(InvSquaHeight + Np*k, Height + Np*k, Hcrit, Np);
		for (p = 0; p < Np; p++)
			K33[k*Np + p] = SQPSPX[k*Np + p] + SQPSPY[k*Np + p] + \
			InvSquaHeight[k*Np + p];
	}

//	int row, col;
	double *sr;
	mwIndex *irs, *jcs;

	plhs[0] = mxCreateSparse(Np*K, Np*K, Nonzero, mxREAL);
	sr = mxGetPr(plhs[0]);
	irs = mxGetIr(plhs[0]);
	jcs = mxGetJc(plhs[0]);
 
	for (i = 0; i < K*Np + 1; i++){
		jcs[i] = (mwIndex)(NonJc[i]);
	}
  
	for (i = 0; i < Nonzero; i++){
		irs[i] = (mwIndex)(NonIr[i]);
        sr[i] += 1.0e-16;
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (i = 0; i < IENe; i++){
		FetchInnerEdgeFacialValue(IEK33M + i*IENfp, IEK33P + i*IENfp, K33, \
			IEFToE + i * 2, IEFToN1 + i*IENfp, IEFToN2 + i*IENfp, Np, IENfp);
		CalculatePenaltyParameterAndWeight(InnerEdgeTau + i, IEWeight1 + i*IENfp, IEWeight2 + i*IENfp,\
			IEK33M + i*IENfp, IEK33P + i*IENfp, IEnx + i*IENfp, IEny + i*IENfp, IEnz + i*IENfp, IENfp, \
			IELAV[i], LAV[(int)IEFToE[2*i]-1], LAV[(int)IEFToE[2 * i + 1]-1], max(N, Nz), Nface);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (i = 0; i < BENe; i++) {
		FetchBoundaryEdgeFacialValue(BEK33M + i*BENfp, K33, \
			BEFToE + i * 2, BEFToN1 + i*BENfp, Np, BENfp);
		CalculatePenaltyParameterAndWeight(BoundaryEdgeTau + i, BEWeight1 + i*BENfp, BEWeight2 + i*BENfp, \
			BEK33M + i*BENfp, BEK33M + i*BENfp, BEnx + i*BENfp, BEny + i*BENfp, BEnz + i*BENfp, BENfp, \
			BELAV[i], LAV[(int)BEFToE[2 * i] - 1], LAV[(int)BEFToE[2 * i + 1] - 1], max(N, Nz), Nface);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (i = 0; i < BotENe; i++){
		FetchInnerEdgeFacialValue(BotEK33M + i*BotENfp, BotEK33P + i*BotENfp, K33, \
			BotEFToE + i * 2, BotEFToN1 + i*BotENfp, BotEFToN2 + i*BotENfp, Np, BotENfp);
		CalculatePenaltyParameterAndWeight(BottomEdgeTau + i, BotEWeight1 + i*BotENfp, BotEWeight2 + i*BotENfp, \
			BotEK33M + i*BotENfp, BotEK33P + i*BotENfp, BotEnx + i*BotENfp, BotEny + i*BotENfp, BotEnz + i*BotENfp,\
			BotENfp, BotELAV[i], LAV[(int)BotEFToE[2 * i]-1], LAV[(int)BotEFToE[2 * i + 1]-1], max(N, Nz), Nface);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (i = 0; i < SurfBENe; i++){
		FetchBoundaryEdgeFacialValue(SurfEK33M + i*SurfBENfp, K33, \
			SurfBEFToE + i * 2, SurfBEFToN1 + i*SurfBENfp, Np, SurfBENfp);
		CalculatePenaltyParameterAndWeight(SurfaceEdgeTau + i, SurfEWeight1 + i*SurfBENfp, \
			SurfEWeight2 + i*SurfBENfp, SurfEK33M + i*SurfBENfp, SurfEK33M + i*SurfBENfp, \
			SurfBEnx + i*SurfBENfp, SurfBEny + i*SurfBENfp, SurfBEnz + i*SurfBENfp, SurfBENfp,\
			SurfBELAV[i], LAV[(int)SurfBEFToE[2 * i]-1], LAV[(int)SurfBEFToE[2 * i + 1]-1], max(N, Nz), Nface);
	}
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (ele = 0; ele < K; ele++){

		GetLocalVolumnIntegralTerm(sr, irs, jcs, \
			Np, M3d, Dr, Ds, Dt, rx + ele*Np, sx + ele*Np, ry + ele*Np, \
			sy + ele*Np, tz + ele*Np, J + ele*Np, K13 + ele*Np, \
			K23 + ele*Np, K33 + ele*Np, ele + 1);

		double *Facialnx = SortedIEnx + ele * Nface2d * IENfp;
		double *Facialny = SortedIEny + ele * Nface2d * IENfp; 
		double *GlobalFace = SortedIEGlobalFace + ele * Nface2d;
		double *AdjEle = SortedIEAdjEle + ele * Nface2d;
		double *ReverseFlag = SortedIEReverseFlag + ele * Nface2d;
		int InternalFace = (int)SortedIEInternalFace[ele];
		double *LocalEidM = malloc(IENfp*sizeof(double));
		double *AdjEidM = malloc(IENfp*sizeof(double));

		for (i = 0; i < InternalFace; i++){
			if ((int)ReverseFlag[i] == 0){
				for (int p = 0; p < IENfp; p++){
					LocalEidM[p] = IEFToN1[(int)GlobalFace[i] * IENfp + p];
					AdjEidM[p] = IEFToN2[(int)GlobalFace[i] * IENfp + p];
				}
			}
			else{
				for (p = 0; p < IENfp; p++){
					LocalEidM[p] = IEFToN2[(int)GlobalFace[i] * IENfp + p];
					AdjEidM[p] = IEFToN1[(int)GlobalFace[i] * IENfp + p];
				}
			}

			GetLocalFacialContributionInHorizontalDirection(sr, irs, jcs, Np, IENfp, \
				IEMb, IEJs + (int)GlobalFace[i] * IENfp, LocalEidM, Dr, Ds, Dt, rx + ele*Np, sx + ele*Np, \
				ry + ele*Np, sy + ele*Np, tz + ele*Np, Facialnx + i * IENfp, Facialny + i*IENfp, K13 + ele*Np, \
				K23 + ele*Np, *(InnerEdgeTau + (int)GlobalFace[i]), ele + 1, (int)GlobalFace[i], (int)ReverseFlag[i]);

			if (AdjEle[i] != ele + 1) {
				GetLocalToAdjacentFacialContributionInHorizontalDirection(sr, irs, jcs, Np, IENfp, \
					IEMb, IEJs + (int)GlobalFace[i] * IENfp, LocalEidM, AdjEidM, Dr, Ds, Dt, rx + ele*Np, sx + ele*Np, \
					ry + ele*Np, sy + ele*Np, tz + ele*Np, rx + (int)(AdjEle[i] - 1)*Np, sx + (int)(AdjEle[i] - 1)*Np, \
					ry + (int)(AdjEle[i] - 1)*Np, sy + (int)(AdjEle[i] - 1)*Np, tz + (int)(AdjEle[i] - 1)*Np, Facialnx + i*IENfp,\
					Facialny + i*IENfp, K13 + ele*Np, K23 + ele*Np, K13 + (int)(AdjEle[i] - 1)*Np, K23 + (int)(AdjEle[i] - 1)*Np, \
					*(InnerEdgeTau + (int)GlobalFace[i]), ele + 1, AdjEle[i], (int)GlobalFace[i], (int)ReverseFlag[i]);
			}

		}
		free(LocalEidM);
		free(AdjEidM);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (face = 0; face < BENe; face++) {
		NdgEdgeType type = (NdgEdgeType)ftype[face];
		//if ((type == NdgEdgeClampedDepth) || (type == NdgEdgeClampedVel) || (type == NdgEdgeZeroGrad)) {
		//if ((type == NdgEdgeClampedVel) || (type == NdgEdgeClampedDepth)) {
		if ((type == NdgEdgeClampedVel)) {
			ImposeHorizontalDirichletBoundaryCondition(sr, irs, jcs, Np, BENfp, \
				BEMb, BEJs + face*BENfp, BEFToN1 + face*BENfp, rx + (int)(BEFToE[2 * face] - 1)*Np, \
				sx + (int)(BEFToE[2 * face] - 1)*Np, ry + (int)(BEFToE[2 * face] - 1)*Np, \
				sy + (int)(BEFToE[2 * face] - 1)*Np, tz + (int)(BEFToE[2 * face] - 1)*Np, \
				Dr, Ds, Dt, K13 + (int)(BEFToE[2 * face] - 1)*Np, K23 + (int)(BEFToE[2 * face] - 1)*Np, \
				*(BoundaryEdgeTau + face), (int)BEFToE[2 * face], BEnx + face*BENfp, BEny + face*BENfp);
		}
	}

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (ele = 0; ele < K2d; ele++){
		for (L = 0; L < Nlayer; L++){
			//Index of the studied element
			int LocalEle = ele*Nlayer + L + 1;
			//Index of the element that located upside of the studied element 
			int UpEle = (int)EToE[ele*Nlayer*Nface + L*Nface + Nface - 1];
			//Index of the element that located downside of the studied element
			int DownEle = (int)EToE[ele*Nlayer*Nface + L*Nface + Nface - 2];

			if (UpEle == DownEle){//Only one layer in the vertical direction, impose the Dirichlet boundary condition
				if (!strcmp(BoundaryType, "Dirichlet")){
					ImposeDirichletBoundaryCondition(sr, irs, jcs, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, \
						K33 + (LocalEle - 1)*Np, *(SurfaceEdgeTau + ele), LocalEle);
				}
			}
			else {//two or more layers included in vertical direction
				if (LocalEle == UpEle) {//This is the top most cell, and L=0
					if (!strcmp(BoundaryType, "Dirichlet")) {
						ImposeDirichletBoundaryCondition(sr, irs, jcs, Np, Np2d, M2d, \
							J2d + ele*Np2d, UpEidM, rx + (LocalEle - 1)*Np, \
							sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
							tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, \
							K33 + (LocalEle - 1)*Np, *(SurfaceEdgeTau + ele), LocalEle);
					}
					GetLocalDownFacialContribution(sr, irs, jcs, LocalEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, BotEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, \
						K33 + (LocalEle - 1)*Np, *(BottomEdgeTau + L*K2d + ele), L*K2d + ele);
					GetLocalToDownFacialContribution(sr, irs, jcs, LocalEle, DownEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, BotEidM, UpEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (DownEle - 1)*Np, sx + (DownEle - 1)*Np, \
						ry + (DownEle - 1)*Np, sy + (DownEle - 1)*Np, tz + (DownEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (DownEle - 1)*Np, K23 + (DownEle - 1)*Np, K33 + (DownEle - 1)*Np, \
						*(BottomEdgeTau + L*K2d + ele), L*K2d + ele);
				}
				else if (LocalEle == DownEle) {// This is the bottom most cell.
					GetLocalUpFacialContribution(sr, irs, jcs, LocalEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, \
						K33 + (LocalEle - 1)*Np, *(BottomEdgeTau + (L - 1)*K2d + ele), (L - 1)*K2d + ele);

					GetLocalToUpFacialContribution(sr, irs, jcs, LocalEle, UpEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, BotEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (UpEle - 1)*Np, sx + (UpEle - 1)*Np, \
						ry + (UpEle - 1)*Np, sy + (UpEle - 1)*Np, tz + (UpEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (UpEle - 1)*Np, K23 + (UpEle - 1)*Np, K33 + (UpEle - 1)*Np, \
						*(BottomEdgeTau + (L - 1)*K2d + ele), (L - 1)*K2d + ele);
				}
				else {
					GetLocalUpFacialContribution(sr, irs, jcs, LocalEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, \
						K33 + (LocalEle - 1)*Np, *(BottomEdgeTau + (L - 1)*K2d + ele), (L - 1)*K2d + ele);

					GetLocalToUpFacialContribution(sr, irs, jcs, LocalEle, UpEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, BotEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (UpEle - 1)*Np, sx + (UpEle - 1)*Np, \
						ry + (UpEle - 1)*Np, sy + (UpEle - 1)*Np, tz + (UpEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (UpEle - 1)*Np, K23 + (UpEle - 1)*Np, K33 + (UpEle - 1)*Np,\
						*(BottomEdgeTau + (L - 1)*K2d + ele), (L - 1)*K2d + ele);

					GetLocalDownFacialContribution(sr, irs, jcs, LocalEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, BotEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, \
						K33 + (LocalEle - 1)*Np, *(BottomEdgeTau + L*K2d + ele), L*K2d + ele);

					GetLocalToDownFacialContribution(sr, irs, jcs, LocalEle, DownEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, BotEidM, UpEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (DownEle - 1)*Np, sx + (DownEle - 1)*Np, \
						ry + (DownEle - 1)*Np, sy + (DownEle - 1)*Np, tz + (DownEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (DownEle - 1)*Np, K23 + (DownEle - 1)*Np, K33 + (DownEle - 1)*Np,\
						*(BottomEdgeTau + L*K2d + ele), L*K2d + ele);
				}
			}					
		}
	}
	free(UpEidM);
	free(BotEidM);
	
}


void GetLocalVolumnIntegralTerm(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, double *M3d, \
	double *Dr, double *Ds, double *Dt, double *rx, double *sx, double *ry, double *sy, double *tz,\
	double *J, double *K13, double *K23, double *K33, int LocalEle){
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));
	double *VertDiffMatrix = malloc(Np*Np*sizeof(double));
	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, M3d, J, Np);
	double *TempContributionBuff = malloc(Np*Np*sizeof(double));
	double *ContributionBuff = malloc(Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	memset(Contribution, 0, Np*Np*sizeof(double));

	/*For term $$k_{11}\frac{\partial v}{\partial x}\frac{\partial p_h}{\partial x}$$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);
	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{13}\frac{\partial q_h}{\partial \sigma}\frac{\partial v}{\partial x}$$.*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);
	DiagMultiply(VertDiffMatrix, VertDiffMatrix, K13, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VertDiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{22}\frac{\partial v}{\partial y}\frac{\partial p_h}{\partial y}$$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);
	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{23}\frac{\partial p_h}{\partial \sigma}\frac{\partial v}{\partial y}$$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);
	DiagMultiply(VertDiffMatrix, VertDiffMatrix, K23, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VertDiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{31}\frac{\partial p_h}{\partial x}\frac{\partial v}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, VertDiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{32}\frac{\partial q_h}{\partial y}\frac{\partial v}{\partial \sigma}$$. Here, $k_{32}=k_{23}$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, VertDiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{33}\frac{\partial p_h}{\partial \sigma}\frac{\partial v}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, VertDiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	AssembleVolumnContributionIntoSparseMatrix(dest, Ir, Jc, Np, Contribution, LocalEle);

	free(DiffMatrix);
	free(TempDiffMatrix);
	free(VertDiffMatrix);
	free(Mass3d);
	free(TempContributionBuff);
	free(ContributionBuff);
	free(Contribution);
}


void GetLocalToDownFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int DownEle, int Np,\
	int Np2d, double *M2d, double *J2d, double *LocalEid, double *DownEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, \
	double *Adjry, double *Adjsy, double *Adjtz, double *Dr, double *Ds, double *Dt,\
	double *LocalK13, double *LocalK23, double *LocalK33, \
	double *AdjK13, double *AdjK23, double *AdjK33,\
	double Tau, int face){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Np2d, Np);
	// Weight defined on the adjacent face used, since the test function v is defined over the adjacent face
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p}{\partial x}\right\}_{\omega}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial p}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Weight defined on the local face used, since the non-hydrostatic pressure p is defined over the local face
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Np2d, Np);
	// Weight defined on the adjacent face used, since the test function v is defined over the adjacent face
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p}{\partial y}\right\}_{\omega}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial p}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Weight defined on the local face used, since the non-hydrostatic pressure p is defined over the local face
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Adjtz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Np2d, Np);
	// Weight defined on the adjacent face is used, since v is defined over the adjacent face
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p}{\partial \sigma}\right\}_{\omega}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Localtz, Np);
	/*For term $k_{33}\frac{\partial p}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Weight defined on the local face is used, since p is defined over the local face
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, DownEid, LocalEid, Np, Np2d, 1.0);

	double *SortedDownEid = malloc(Np2d*sizeof(double));
	memcpy(SortedDownEid, DownEid, Np2d*sizeof(double));
	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Np2d*sizeof(double));
	Sort(SortedDownEid, Np2d);
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedLocalEid, SortedDownEid, Np, Np2d, TempContribution, LocalEle, DownEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempMass2d);
	free(SortedDownEid);
	free(SortedLocalEid);
}


void GetLocalToUpFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int UpEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *UpEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33, double Tau, int face){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Np2d, Np);
	// Weight defined on the local face is used, since v is defined over the local face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p}{\partial x}\right\}_{\omega}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial p}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Weight defined on the adjacent face is used, since p is defined over the adjacent face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Np2d, Np);
	// Weight defined on the local face is used, since v is defined over the local face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p}{\partial y}\right\}_{\omega}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial p}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Weight defined on the adjacent face is used, since p is defined over the adjacent face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Adjtz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Np2d, Np);
	// Weight defined on the local face is used, since v is defined over the local face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p}{\partial \sigma}\right\}_{\omega}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Localtz, Np);
	/*For term $k_{33}\frac{\partial p}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Weight defined on the adjacent face is used, since p is defined over the adjacent face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, UpEid, LocalEid, Np, Np2d, 1.0);

	double *SortedUpEid = malloc(Np2d*sizeof(double));
	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedUpEid, UpEid, Np2d*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Np2d*sizeof(double));
	Sort(SortedUpEid, Np2d);
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedLocalEid, SortedUpEid, Np, Np2d, TempContribution, LocalEle, UpEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempMass2d);
	free(SortedUpEid);
	free(SortedLocalEid);
}


void GetLocalDownFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23,\
	double *K33, double Tau, int face){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Both p and v are defined over the local face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);


	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}_{\omega}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Both p and v are defined over the local face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}_{\omega}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Both p and v are defined over the local face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight1 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}_{\omega}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1.0, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Np2d*sizeof(double));
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedLocalEid, SortedLocalEid, Np, Np2d, TempContribution, LocalEle, LocalEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempMass2d);
	free(SortedLocalEid);
}

void GetLocalUpFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33, \
	double Tau, int face){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Both p and v are defined over the adjacent face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}_{\omega}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Both p and v are defined over the adjacent face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}_{\omega}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}_{\omega}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	// Both p and v are defined over the adjacent face defined according to FToE
	DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, BotEWeight2 + face * Np2d, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);


	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}_{\omega}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Np2d*sizeof(double));
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedLocalEid, SortedLocalEid, Np, Np2d, TempContribution, LocalEle, LocalEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempMass2d);
	free(SortedLocalEid);
}

void ImposeDirichletBoundaryCondition(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx,\
	double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, \
	double *K23, double *K33, double Tau, int LocalEle){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *BoundaryEdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);

	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Np2d*sizeof(double));
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, Ir, Jc, SortedLocalEid, SortedLocalEid, Np, Np2d, TempContribution, LocalEle, LocalEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(BoundaryEdgeContribution);
	free(TempMass2d);
	free(SortedLocalEid);
}

void GetLocalToAdjacentFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *M2d, double *J2d, double *LocalEid, double *AdjEid, double *Dr, double *Ds, double *Dt, \
	double *LocalRx, double *LocalSx, double *LocalRy, double *LocalSy, double *LocalTz, \
	double *AdjRx, double *AdjSx, double *AdjRy, double *AdjSy, double *AdjTz, double *nx, double *ny,\
	double *LocalK13, double *LocalK23, double *AdjK13, double *AdjK23, double Tau, int LocalEle,\
	int AdjEle, int Face, int Flag){
	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *TempAdjDiff = malloc(Np*Np * sizeof(double));
	double *AdjDiff = malloc(Np*Np*sizeof(double));
	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));
	double *InnerEdgeContribution = malloc(Np*Nfp*sizeof(double));

	/*For term $$ \left \{k_{11}\frac{\partial v}{\partial x}\right\}_{\omega}[p_h]_x $$*/
	DiagMultiply(TempAdjDiff, Dr, AdjRx, Np);
	DiagMultiply(AdjDiff, Ds, AdjSx, Np);
	Add(AdjDiff, AdjDiff, TempAdjDiff, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, AdjDiff, AdjEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on adjacent element defined by FToE
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, nx[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{13}\frac{\partial v}{\partial \sigma}\right\}_{\omega}[p_h]_x $$*/
	DiagMultiply(AdjDiff, Dt, AdjTz, Np);
	/*For term $k_{13}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(AdjDiff, AdjDiff, AdjK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, AdjDiff, AdjEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on adjacent element defined by FToE
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, nx[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{22}\frac{\partial v}{\partial y}\right\}_{\omega}[p_h]_y $$*/
	DiagMultiply(TempAdjDiff, Dr, AdjRy, Np);
	DiagMultiply(AdjDiff, Ds, AdjSy, Np);
	Add(AdjDiff, AdjDiff, TempAdjDiff, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, AdjDiff, AdjEid, Nfp, Np);

	if (Flag == 0) {
		// test function v lies on adjacent element defined by FToE
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, ny[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);


	/*For term $$ \left \{k_{23}\frac{\partial v}{\partial \sigma}\right\}_{\omega}[p_h]_y $$*/
	DiagMultiply(AdjDiff, Dt, AdjTz, Np);
	/*For term $k_{23}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(AdjDiff, AdjDiff, AdjK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, AdjDiff, AdjEid, Nfp, Np);

	if (Flag == 0) {
		// test function v lies on adjacent element defined by FToE
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, ny[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	double *TempLocalDiff = malloc(Np*Np * sizeof(double));
	double *LocalDiff = malloc(Np*Np * sizeof(double));

	/*For term $$ \left \{k_{11}\frac{\partial p_h}{\partial x}\right\}_{\omega}[v]_x $$*/
	DiagMultiply(TempLocalDiff, Dr, LocalRx, Np);
	DiagMultiply(LocalDiff, Ds, LocalSx, Np);
	Add(LocalDiff, LocalDiff, TempLocalDiff, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);

	if (Flag == 0) {
		// test function v lies on adjacent element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}

	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -1.0*nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEid, Np, Nfp);


	/*For term $$k_{13}\frac{\partial p_h}{\partial \sigma}_{\omega} [v]_x$$*/
	DiagMultiply(LocalDiff, Dt, LocalTz, Np);
	/*For term $k_{13}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on adjacent element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}

	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -1.0*nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEid, Np, Nfp);

	/*For term $$ \left \{k_{22}\frac{\partial p_h}{\partial y}\right\}_{\omega}[v]_y $$*/
	DiagMultiply(TempLocalDiff, Dr, LocalRy, Np);
	DiagMultiply(LocalDiff, Ds, LocalSy, Np);
	Add(LocalDiff, LocalDiff, TempLocalDiff, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on adjacent element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}

	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -1.0*ny[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEid, Np, Nfp);

	/*For term $$k_{23}\frac{\partial p_h}{\partial \sigma}_{\omega} [v]_y$$*/
	DiagMultiply(LocalDiff, Dt, LocalTz, Np);
	/*For term $k_{23}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on adjacent element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}

	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);

	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -1.0*ny[0], Np*Nfp);

	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEid, Np, Nfp);

	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));

	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Nfp*Nfp);

	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, AdjEid, LocalEid, Np, Nfp, 1.0);

	double *SortedAdjEid = malloc(Nfp*sizeof(double));
	memcpy(SortedAdjEid, AdjEid, Nfp*sizeof(double));
	double *SortedLocalEid = malloc(Nfp*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Nfp*sizeof(double));
	Sort(SortedAdjEid, Nfp);
	Sort(SortedLocalEid, Nfp);

	AssembleFacialContributionIntoSparseMatrix(dest, Ir, Jc, SortedLocalEid, SortedAdjEid, Np, Nfp, TempContribution, LocalEle, AdjEle);

	free(FacialMass2d);
	free(TempContribution);
	free(TempAdjDiff);
	free(AdjDiff);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(TempLocalDiff);
	free(LocalDiff);
	free(TempMass2d);
	free(SortedLocalEid);
	free(SortedAdjEid);
}

void GetLocalFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *M2d, double *J2d, double *LocalEid, double *Dr, double *Ds, double *Dt, double *rx, double *sx,\
	double *ry, double *sy, double *tz, double *nx, double *ny, double *K13, double *K23, double Tau,\
	int LocalEle, int Face, int Flag){

	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *TempLocalDiff = malloc(Np*Np * sizeof(double));
	double *LocalDiff = malloc(Np*Np*sizeof(double));
	double *FacialDiffMatrix = malloc(Np*Nfp * sizeof(double));
	double *InnerEdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{11}\frac{\partial v}{\partial x}\right\}_{\omega}[p_h]_x $$*/
	/*For term $k_{11}\frac{\partial v}{\partial x}$, here $k_{11}=1$*/
	DiagMultiply(TempLocalDiff, Dr, rx, Np);
	DiagMultiply(LocalDiff, Ds, sx, Np);
	Add(LocalDiff, LocalDiff, TempLocalDiff, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on local element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, nx[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{13}\frac{\partial v}{\partial \sigma}\right\}_{\omega}[p_h]_x $$*/
	/*For term $k_{13}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	DiagMultiply(LocalDiff, LocalDiff, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on local element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, nx[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{22}\frac{\partial v}{\partial y}\right\}_{\omega}[p_h]_y $$*/
	/*For term $k_{22}\frac{\partial v}{\partial y}$, here $k_{22}=1$*/
	DiagMultiply(TempLocalDiff, Dr, ry, Np);
	DiagMultiply(LocalDiff, Ds, sy, Np);
	Add(LocalDiff, LocalDiff, TempLocalDiff, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on local element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, ny[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial v}{\partial \sigma}\right\}_{\omega}[p_h]_y $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{23}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on local element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, ny[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{11}\frac{\partial q_h}{\partial x}\right\}_{\omega}[v]_x $$*/
	/*For term $k_{11}\frac{\partial q_h}{\partial x}$, here $k_{11}=1$*/
	DiagMultiply(TempLocalDiff, Dr, rx, Np);
	DiagMultiply(LocalDiff, Ds, sx, Np);
	Add(LocalDiff, LocalDiff, TempLocalDiff, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on local element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{13}\frac{\partial p_h}{\partial \sigma}\right\}_{\omega}[v]_x $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{13}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on local element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{22}\frac{\partial q_h}{\partial y}\right\}_{\omega}[v]_y $$*/
	/*For term $k_{22}\frac{\partial q_h}{\partial y}$, here $k_{11}=1$*/
	DiagMultiply(TempLocalDiff, Dr, ry, Np);
	DiagMultiply(LocalDiff, Ds, sy, Np);
	Add(LocalDiff, LocalDiff, TempLocalDiff, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on local element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, ny[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial p_h}{\partial \sigma}\right\}_{\omega}[v]_y $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{23}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	if (Flag == 0) {
		// test function v lies on local element defined by FToE, p_h on local element, IEWeight1 used
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight1 + Face * Nfp, Nfp, Np);
	}
	else {
		DiagLeftMultiplyUnsymmetric(FacialDiffMatrix, FacialDiffMatrix, IEWeight2 + Face * Nfp, Nfp, Np);
	}
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, ny[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*The follow part is for term $ - \int_{ \epsilon_i }\tau^k[p_h][s]d\boldsymbol{ x }$*/
	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Nfp*Nfp);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Nfp, -1.0);

	double *SortedLocalEidM = malloc(Nfp*sizeof(double));
	memcpy(SortedLocalEidM, LocalEid, Nfp*sizeof(double));
	Sort(SortedLocalEidM, Nfp);
	// Why
	AssembleFacialContributionIntoSparseMatrix(dest, Ir, Jc, SortedLocalEidM, SortedLocalEidM, Np, Nfp, TempContribution, LocalEle, LocalEle);
	free(FacialMass2d);
	free(TempContribution);
	free(TempLocalDiff);
	free(LocalDiff);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(TempMass2d);
	free(SortedLocalEidM);
}

void ImposeHorizontalDirichletBoundaryCondition(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, \
	double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, \
	double *K13, double *K23, double Tau, int LocalEle, double *nx, double *ny) {

	double FacialMass2d[Np2d*Np2d];
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double TempContribution[Np*Np];
	memset(TempContribution, 0, Np*Np * sizeof(double));

	double DiffMatrix[Np*Np];
	double TempDiffMatrix[Np*Np];
	double FacialDiffMatrix[Np*Np2d];
	double BoundaryEdgeContribution[Np*Np2d];

	/*For term $$ \left \{\frac{\partial v}{\partial x}\right\}p_h n_x $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(BoundaryEdgeContribution, BoundaryEdgeContribution, nx[0], Np*Np2d);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);
	/*For term $$ \left \{\frac{\partial p_h}{\partial x}\right\}v n_x $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);
	MultiplyByConstant(BoundaryEdgeContribution, BoundaryEdgeContribution, nx[0], Np*Np2d);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$\frac{\partial v}{\partial \sigma}\frac{\partial\sigma}{\partial x}q_h n_x $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(BoundaryEdgeContribution, BoundaryEdgeContribution, nx[0], Np*Np2d);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);
	/*For term $$\frac{\partial p_h}{\partial \sigma}\frac{\partial\sigma}{\partial x}v n_x $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);
	MultiplyByConstant(BoundaryEdgeContribution, BoundaryEdgeContribution, nx[0], Np*Np2d);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{\frac{\partial v}{\partial x}\right\}p_h n_y $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(BoundaryEdgeContribution, BoundaryEdgeContribution, ny[0], Np*Np2d);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);
	/*For term $$ \left \{\frac{\partial p_h}{\partial y}\right\}v n_y $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);
	MultiplyByConstant(BoundaryEdgeContribution, BoundaryEdgeContribution, ny[0], Np*Np2d);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$\frac{\partial v}{\partial \sigma}\frac{\partial\sigma}{\partial y}q_h n_y $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(BoundaryEdgeContribution, BoundaryEdgeContribution, ny[0], Np*Np2d);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);
	/*For term $$\frac{\partial p_h}{\partial \sigma}\frac{\partial\sigma}{\partial x}v n_y $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);
	MultiplyByConstant(BoundaryEdgeContribution, BoundaryEdgeContribution, ny[0], Np*Np2d);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*The penalty term*/
	double TempMass2d[Np2d*Np2d];
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	double SortedLocalEid[Np2d];
	memcpy(SortedLocalEid, LocalEid, Np2d * sizeof(double));
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, Ir, Jc, SortedLocalEid, SortedLocalEid, Np, Np2d, TempContribution, LocalEle, LocalEle);
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
