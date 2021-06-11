
#include "SWENonhydrostatic3d.h"

void GetInverseSquareHeight(double *, double *, double , int );

void GetPenaltyParameter(double *, double , double , int , int , int );

void ImposeNewmannBoundaryCondition(double *, mwIndex *, mwIndex *, int , \
	double *, double *, double *, double *, double *, double *, int , int , \
	double *, double *, double *, double *, double *, int , double *);

void ImposeDirichletBoundaryCondition(double *, mwIndex *, mwIndex *, int, \
	int, int, double *, double *, double *, double *, double *, double *, double *, double *, \
	double *, double *, double *, double *, double *, double *, double *, double *, double *, \
	int, double *, double *, double *, double *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
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
	mxArray *TempP = mxGetField(cell, 0, "N");
	int P = (int)mxGetScalar(TempP);

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
	mxArray *Temprz = mxGetField(mesh, 0, "rz");
	double *rz = mxGetPr(Temprz);
	mxArray *Temptz = mxGetField(mesh, 0, "tz");
	double *tz = mxGetPr(Temptz);

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
	mxArray *TempBEnx = mxGetField(BoundaryEdge, 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(BoundaryEdge, 0, "ny");
	double *BEny = mxGetPr(TempBEny);
	mxArray *TempBEnz = mxGetField(BoundaryEdge, 0, "nz");
	double *BEnz = mxGetPr(TempBEnz);
	mxArray *TempBELMass2d = mxGetField(BoundaryEdge, 0, "M");
	double  *BELMass2d = mxGetPr(TempBELMass2d);
	mxArray *TempBEFLAV = mxGetField(BoundaryEdge, 0, "LAV");
	double *BEFLAV = mxGetPr(TempBEFLAV);

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

	double *InvSquaHeight = malloc(Np*K*sizeof(double));
	double *K33 = malloc(Np*K*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		GetInverseSquareHeight(InvSquaHeight + Np*k, Height + Np*k, Hcrit, Np);
		for (int p = 0; p < Np; p++)
			K33[k*Np + p] = pow(K13[k*Np + p], 2) + pow(K23[k*Np + p], 2) + \
			InvSquaHeight[k*Np + p];
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int edge = 0; edge < BENe; edge++){
		if ((NdgEdgeType)ftype[edge] == NdgEdgeSlipWall){
			//Newmann boundary Doing Nothing
		}
		else{

			double *FpIndex = malloc(BENfp*sizeof(double));

			for (int p = 0; p < BENfp; p++){
				FpIndex[p] = BEFToN1[BENfp*edge + p];
			}

			int LocalEle;
			LocalEle = (int)BEFToE[2 * edge];
			double *Tau = malloc(BENfp*sizeof(double));
			GetPenaltyParameter(Tau, LAV[LocalEle - 1], BEFLAV[edge], P, Nface, BENfp);

			double *TempEToE = NULL, *TempJ = NULL, *TempJs = NULL;
			TempEToE = EToE + (LocalEle - 1)*Nface;
			TempJ = J + (LocalEle - 1)*Np;
			TempJs = BEJs + edge * BENfp;

			ImposeDirichletBoundaryCondition(sr, irs, jcs, LocalEle, \
				Np, BENfp, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, \
				sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
				Dr, Ds, Dt, Tau, BEnx + edge * BENfp, BEny + edge * BENfp, BEnz + edge*BENfp, \
				Mass3d, TempJ, TempJs, BELMass2d, TempEToE, Nface, FpIndex, \
				K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np);

			free(FpIndex);
			free(Tau);
		}
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

		double *TempEToE = NULL, *TempJ = NULL, *TempJs = NULL;
		TempEToE = EToE + (LocalEle - 1)*Nface;
		TempJ = J + (LocalEle - 1)*Np;
		TempJs = BotBEJs + edge * BotBENfp;

		double *FpIndex = malloc(BotBENfp*sizeof(double));

		for (int p = 0; p < BotBENfp; p++){
			FpIndex[p] = BotBEFToN1[BotBENfp*edge + p];
		}

		ImposeNewmannBoundaryCondition(sr, irs, jcs, LocalEle, \
			K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np,\
			Dx, Dy, Dz, Np, BotBENfp, \
			Mass3d, TempJ, TempJs, BotBELMass2d, TempEToE, Nface, FpIndex);

		free(DxBuff);
		free(Dx);
		free(Dz);
		free(DyBuff);
		free(Dy);
		free(FpIndex);
	}

	free(InvSquaHeight);
	free(K33);
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

void GetPenaltyParameter(double *dest, double LAV, double FLAV, int P, int Nface, int Nfp){
	for (int i = 0; i < Nfp; i++){
		dest[i] = (P + 1)*(P + 3) / 3.0 * Nface / 2.0 * FLAV / LAV;
		//	dest[i] = 2 * 1.0 / sqrt(FLAV); //This parameter is doubled at the boundary
	}
}

void ImposeNewmannBoundaryCondition(double *dest, mwIndex *Irs, mwIndex *Jcs, int LocalEle, \
	double *K31, double *K32, double *K33, double *Dx, double *Dy, double *Dz, int Np, int Nfp,\
	double *M3d, double *J, double *Js, double *M2d, double *EToE, int Nface, double *FpIndex){

	double *TempK31 = malloc(Np*Np*sizeof(double));
	double *TempK32 = malloc(Np*Np*sizeof(double));
	double *TempK33 = malloc(Np*Np*sizeof(double));
	double *TempCoe = malloc(Np*Np*sizeof(double));

	DiagMultiply(TempK31, Dx, K31, Np);

	DiagMultiply(TempK32, Dy, K32, Np);

	DiagMultiply(TempK33, Dz, K33, Np);

	Add(TempCoe, TempK31, TempK32, Np*Np);

//	Add(TempCoe, TempCoe, TempK33, Np*Np);

	/*Withdraw the data in $\frac{\partial p}{\partial \sigma}$, and store them in TempPNPS*/
	double *TempFacialData = malloc(Nfp*sizeof(double));
	double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(EleMass2d, M2d, Js, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	double *InvEleMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(InvEleMass3d, M3d, J, Np);
	MatrixInverse(InvEleMass3d, (ptrdiff_t)Np);

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

	for (int j = 0; j < Nfp; j++){

		FetchFacialData(TempFacialData, TempCoe + ((int)FpIndex[j] - 1)*Np, FpIndex, Nfp);

		MatrixMultiply("N", "N", (ptrdiff_t)Nfp, One, (ptrdiff_t)Nfp, 1.0, EleMass2d,
			(ptrdiff_t)Nfp, TempFacialData, (ptrdiff_t)Nfp, 0.0, ContributionPerPoint, (ptrdiff_t)Nfp);

		MultiplyByConstant(ContributionPerPoint, ContributionPerPoint, -1.0, Nfp);

		AssembleDataIntoPoint(TempContribution + ((int)FpIndex[j] - 1)*Np, ContributionPerPoint, FpIndex, Nfp);
	}

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvEleMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(TempK31);
	free(TempK32);
	free(TempK33);
	free(TempCoe);
	free(TempFacialData);
	free(EleMass2d);
	free(TempContribution);
	free(Contribution);
	free(InvEleMass3d);
	free(TempEToE);
	free(ContributionPerPoint);
}

void ImposeDirichletBoundaryCondition(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, \
	int Np, int Nfp, double *rx, double *sx, double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, \
	double *Tau, double *nx, double *ny, double *nz, double *Mass3d, double *J, double *Js, double *Mass2d, double *EToE, \
	int Nface, double *FpIndex, double *K13, double *K23, double *K33){
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

	double *EleMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(EleMass3d, Mass3d, J, Np);
	double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(EleMass2d, Mass2d, Js, Nfp);
	double *InvEleMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvEleMass3d, EleMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvEleMass3d, (ptrdiff_t)Np);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

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
	DiagMultiply(TempMass2d, EleMass2d, Tau, Nfp);
	/*For term $-\int_{\partial \Omega^d}\tau^k s u_hd\boldsymbol{x}$*/
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, FpIndex, FpIndex, Np, Nfp, -1.0);

	/*Multiply the contribution by inverse matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvEleMass3d, \
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	int UniNum = 0, StartPoint;
	/*Find the exact place where to fill in the data, and the place is stored in StartPoint*/
	double *TempEToE = malloc((Nface + 1)*sizeof(double));
	FindUniqueElementAndSortOrder(TempEToE, EToE, &UniNum, Nface, LocalEle);
	int NonzeroPerColumn = jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np];
	for (int j = 0; j < UniNum; j++){
		if ((int)TempEToE[j] == LocalEle){
			StartPoint = jcs[(LocalEle - 1)*Np] + j*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}

	free(DxBuff);
	free(Dx);
	free(DyBuff);
	free(Dy);
	free(Dz);
	free(EleMass3d);
	free(EleMass2d);
	free(InvEleMass3d);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempEToE);
	free(TempMass2d);
	free(TempDiffMatrix);
}