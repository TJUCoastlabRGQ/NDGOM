#include "SWENonhydrostatic3d.h"

void GetInverseSquareHeight(double *, double *, double , int );

void ImposeNewmannBoundaryCondition(double *RHSdest, int LocalEle, \
	int Np, int Nfp, double *Js, double *M2d, double *FpIndex, double *NewmannData);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *K13 = mxGetPr(prhs[0]);
	double *K23 = mxGetPr(prhs[1]);
	double Hcrit = mxGetScalar(prhs[2]);

	const mxArray *mesh = prhs[3];
	const mxArray *cell = prhs[4];
	const mxArray *BottomBoundaryEdge = prhs[5];
	double *Uold = mxGetPr(prhs[6]);
	double *Unew = mxGetPr(prhs[7]);

	double *Vold = mxGetPr(prhs[8]);
	double *Vnew = mxGetPr(prhs[9]);

	double *Wold = mxGetPr(prhs[10]);
	double *Wnew = mxGetPr(prhs[11]);

	double deltatime = mxGetScalar(prhs[12]);
	double rho = mxGetScalar(prhs[13]);
	double *varIndex = mxGetPr(prhs[14]);
	double *RHS = mxGetPr(prhs[15]);
	double *frhs = mxGetPr(prhs[16]);
	double *fphys = mxGetPr(prhs[17]);

	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);

	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	double *H3d = fphys + (int)(varIndex[3]-1) * Np*K;

	plhs[0] = mxCreateDoubleMatrix(Np*K, 1, mxREAL);
	double *OutRHS = mxGetPr(plhs[0]);
	memcpy(OutRHS, RHS, Np*K*sizeof(double));

	mxArray *TempBotBEJs = mxGetField(BottomBoundaryEdge, 0, "Js");
	double *BotBEJs = mxGetPr(TempBotBEJs);
	mxArray *TempBotBENfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
	int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBEFToE = mxGetField(BottomBoundaryEdge, 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToN1 = mxGetField(BottomBoundaryEdge, 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);
	mxArray *TempBotBEMass2d = mxGetField(BottomBoundaryEdge, 0, "M");
	double  *BotBEMass2d = mxGetPr(TempBotBEMass2d);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		double H2d[BotBENfp], Urhs[BotBENfp], Vrhs[BotBENfp], Wrhs[BotBENfp], \
			PsPx[BotBENfp], PsPy[BotBENfp], PqPs[BotBENfp], PqPx[BotBENfp], PqPy[BotBENfp];
		FetchBoundaryEdgeFacialValue(H2d, H3d, BotBEFToE + face * 2, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(Urhs,frhs, BotBEFToE + face * 2, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(Vrhs,frhs + Np*K,BotBEFToE + face * 2, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(Wrhs,frhs + 2*Np*K,BotBEFToE + face * 2, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(PsPx,K13,BotBEFToE + face * 2, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(PsPy,K23,BotBEFToE + face * 2, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		for (int p = 0; p < BotBENfp; p++) {
			/*
			PqPs[p] = -1.0*(-1.0*Wrhs[p])*rho*H2d[p];
			PqPx[p] = -1.0*(-1.0*Urhs[p])*rho - \
				H2d[p] * PqPs[p] * PsPx[p];
			PqPy[p] = -1.0*(-1.0*Vrhs[p])*rho - \
				H2d[p] * PqPs[p] * PsPy[p];
			*/
			PqPs[p] = -1.0*(-1.0*Wrhs[p] + (Wnew[face*BotBENfp + p] - Wold[face*BotBENfp + p]) / deltatime)*rho*H2d[p];
			PqPx[p] = -1.0*(-1.0*Urhs[p] + (Unew[face*BotBENfp + p] - Uold[face*BotBENfp + p]) / deltatime)*rho - \
				H2d[p] * PqPs[p] * PsPx[p];
			PqPy[p] = -1.0*(-1.0*Vrhs[p] + (Vnew[face*BotBENfp + p] - Vold[face*BotBENfp + p]) / deltatime)*rho - \
				H2d[p] * PqPs[p] * PsPy[p];
		}
		DotCriticalDivide(PqPx, PqPx, &Hcrit, H2d, BotBENfp);
		DotCriticalDivide(PqPy, PqPy, &Hcrit, H2d, BotBENfp);

		double HInvSquare[BotBENfp], NewmannData[BotBENfp];
		GetInverseSquareHeight(HInvSquare, H2d, Hcrit, BotBENfp);
		for (int p = 0; p < BotBENfp; p++) {
			/*The outward vector and the sign in the primal form considered here*/
			NewmannData[p] = (-1.0) * (-1.0) * (PsPx[p] * PqPx[p] + PsPy[p] * PqPy[p] + \
				(pow(PsPx[p], 2.0) + pow(PsPy[p], 2.0) + HInvSquare[p])*PqPs[p]);
		}
		ImposeNewmannBoundaryCondition(OutRHS + (int)(BotBEFToE[face * 2] - 1)*Np, (int)BotBEFToE[face * 2], \
			Np, BotBENfp, BotBEJs + face*BotBENfp, BotBEMass2d, BotBEFToN1 + face*BotBENfp, NewmannData);
	}
}

void GetInverseSquareHeight(double *Squadest, double *source, double Hcrit, int Np){
	for (int i = 0; i < Np; i++){
		if (source[i] >= Hcrit){
			Squadest[i] = 1.0 / source[i] / source[i];
		}
		else{
			Squadest[i] = 0;
		}
	}
}

void ImposeNewmannBoundaryCondition(double *RHSdest, int LocalEle, \
	int Np, int Nfp, double *Js, double *M2d, double *FpIndex, double *NewmannData){

	double TempRHSBuff[Np];

	memset(TempRHSBuff, 0, Np * sizeof(double));

	double TempRHSFacialData[Nfp];

	double EleMass2d[Nfp*Nfp];

	DiagMultiply(EleMass2d, M2d, Js, Nfp);

	ptrdiff_t One = 1;
	
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, One, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, NewmannData, (ptrdiff_t)Nfp, 0.0, TempRHSFacialData, (ptrdiff_t)Nfp);

	AssembleDataIntoPoint(TempRHSBuff, TempRHSFacialData, FpIndex, Nfp);

	Add(RHSdest, RHSdest, TempRHSBuff, Np);
}