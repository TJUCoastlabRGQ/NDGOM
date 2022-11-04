#include "SWENonhydrostatic3d.h"

void ImposeNewmannBoundaryCondition(double *RHSdest, int LocalEle, \
	int Np, int Nfp, double *Js, double *M2d, double *FpIndex, double *NewmannData);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *K13 = mxGetPr(prhs[0]);
	double *K23 = mxGetPr(prhs[1]);
	double Hcrit = mxGetScalar(prhs[2]);

	const mxArray *mesh = prhs[3];
	const mxArray *cell = prhs[4];
	const mxArray *BoundaryEdge = prhs[5];
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

	signed char *ftype = (signed char *)mxGetData(prhs[18]);


	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);

	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	double *H3d = fphys + (int)(varIndex[3]-1) * Np*K;

	plhs[0] = mxCreateDoubleMatrix(Np*K, 1, mxREAL);
	double *OutRHS = mxGetPr(plhs[0]);
	memcpy(OutRHS, RHS, Np*K*sizeof(double));

	mxArray *TempBEJs = mxGetField(BoundaryEdge, 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	mxArray *TempBENfp = mxGetField(BoundaryEdge, 0, "Nfp");
	int BENfp = (int)mxGetScalar(TempBENfp);
	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBEFToE = mxGetField(BoundaryEdge, 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToN1 = mxGetField(BoundaryEdge, 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);
	mxArray *TempBEMass2d = mxGetField(BoundaryEdge, 0, "M");
	double  *BEMass2d = mxGetPr(TempBEMass2d);
	mxArray *TempBEnx = mxGetField(BoundaryEdge, 0, "nx");
	double  *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(BoundaryEdge, 0, "ny");
	double  *BEny = mxGetPr(TempBEny);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		NdgEdgeType type = (NdgEdgeType)ftype[face];
		if ((type == NdgEdgeClampedDepth) || (type == NdgEdgeClampedVel)) {
		//if ((type == NdgEdgeClampedDepth)) {
			double H2d[BENfp], Urhs[BENfp], Vrhs[BENfp], Wrhs[BENfp], \
				PsPx[BENfp], PsPy[BENfp], PqPs[BENfp], PqPx[BENfp], PqPy[BENfp];
			FetchBoundaryEdgeFacialValue(H2d, H3d, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(Urhs, frhs, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(Vrhs, frhs + Np*K, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(Wrhs, frhs + 2 * Np*K, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(PsPx, K13, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(PsPy, K23, BEFToE + face * 2, BEFToN1 + face*BENfp, Np, BENfp);
			for (int p = 0; p < BENfp; p++) {
				/*
				PqPs[p] = -1.0*(-1.0*Wrhs[p])*rho*H2d[p];
				PqPx[p] = -1.0*(-1.0*Urhs[p])*rho - \
					H2d[p] * PqPs[p] * PsPx[p];
				PqPy[p] = -1.0*(-1.0*Vrhs[p])*rho - \
					H2d[p] * PqPs[p] * PsPy[p];
				*/
				PqPs[p] = -1.0*(-1.0*Wrhs[p] + (Wnew[face*BENfp + p] - Wold[face*BENfp + p]) / deltatime)*rho*H2d[p];
				PqPx[p] = -1.0*(-1.0*Urhs[p] + (Unew[face*BENfp + p] - Uold[face*BENfp + p]) / deltatime)*rho - \
					H2d[p] * PqPs[p] * PsPx[p];
				PqPy[p] = -1.0*(-1.0*Vrhs[p] + (Vnew[face*BENfp + p] - Vold[face*BENfp + p]) / deltatime)*rho - \
					H2d[p] * PqPs[p] * PsPy[p];
			}
			DotCriticalDivide(PqPx, PqPx, &Hcrit, H2d, BENfp);
			DotCriticalDivide(PqPy, PqPy, &Hcrit, H2d, BENfp);

			double NewmannData[BENfp];
			for (int p = 0; p < BENfp; p++) {
				/*The outward vector and the sign in the primal form considered here*/
				NewmannData[p] = (-1.0) * BEnx[face*BENfp + p] * (PqPx[p] + \
					PqPs[p] * PsPx[p]) + (-1.0) * BEny[face*BENfp + p] * (PqPy[p] + \
						PqPs[p] * PsPy[p]);
			}
			ImposeNewmannBoundaryCondition(OutRHS + (int)(BEFToE[face * 2] - 1)*Np, (int)BEFToE[face * 2], \
				Np, BENfp, BEJs + face*BENfp, BEMass2d, BEFToN1 + face*BENfp, NewmannData);
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