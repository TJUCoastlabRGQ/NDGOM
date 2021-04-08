#include "../../../../../../NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver3d/private/SWENonhydrostatic3d.h"

void GetFaceTypeAndFaceOrder(int *, int *, int *, double *, double *, signed char *, int);

void GetPenaltyParameter(double *, double , double , double, int, int);

void ImposeDirichletBoundaryCondition(double *, double *, mwIndex *, mwIndex *, int , \
	int , int , double *, double *, double *, double *, double *, double *, double *, \
	double *, double *, double *, double *, double *, double *, double *, \
	int , double *, double *);

void AssembleDataIntoPoint(double *, double *, double *, int );

void SumInColumn(double *, double *, int );

void SumInRow(double *, double *, int , int );

int GetGlobalFace(int , int , double *, double *, int );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *TempSPNPX = mxGetPr(prhs[0]);
	mwIndex *Tempjcs = mxGetJc(prhs[0]);
	mwIndex *Tempirs = mxGetIr(prhs[0]);
	int row, col;
	row = (int)mxGetM(prhs[0]);
	col = (int)mxGetN(prhs[0]);
	double *TempRHS = mxGetPr(prhs[1]);
	double *DirichDataValue = mxGetPr(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(row, 1, mxREAL);
	double *OutRHS = mxGetPr(plhs[0]);
	memcpy(OutRHS, TempRHS, col*sizeof(double));
	double *sr;
	mwIndex *irs, *jcs;
	plhs[1] = mxCreateSparse(row, col, Tempjcs[col], mxREAL);
	sr = mxGetPr(plhs[1]);
	irs = mxGetIr(plhs[1]);
	jcs = mxGetJc(plhs[1]);
	memcpy(sr, TempSPNPX, Tempjcs[col] * sizeof(double));
	memcpy(irs, Tempirs, Tempjcs[col] * sizeof(mwIndex));
	memcpy(jcs, Tempjcs, (col + 1)*sizeof(mwIndex));

	const mxArray *BoundaryEdge2d = prhs[3];
	const mxArray *BoundaryEdge = prhs[4];
	const mxArray *cell = prhs[5];
	const mxArray *mesh = prhs[6];

	signed char *ftype2d = (signed char *)mxGetData(prhs[7]);

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
	mxArray *TempP = mxGetField(cell, 0, "N");
	double P = mxGetScalar(TempP);


	mxArray *TempBENe2d = mxGetField(BoundaryEdge2d, 0, "Ne");
	int BENe2d = (int)mxGetScalar(TempBENe2d);

	mxArray *TempNlayer = mxGetField(mesh, 0, "Nz");
	int Nlayer = (int)mxGetScalar(TempNlayer);
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

	mxArray *TempJs = mxGetField(BoundaryEdge, 0, "Js");
	double *Js = mxGetPr(TempJs);
	int Nfp = (int)mxGetM(TempJs);
	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempFToF = mxGetField(BoundaryEdge, 0, "FToF");
	double *FToF = mxGetPr(TempFToF);
	mxArray *TempFToE = mxGetField(BoundaryEdge, 0, "FToE");
	double *FToE = mxGetPr(TempFToE);
	mxArray *TempFToN1 = mxGetField(BoundaryEdge, 0, "FToN1");
	double *FToN1 = mxGetPr(TempFToN1);
	mxArray *Tempnx = mxGetField(BoundaryEdge, 0, "nx");
	double *nx = mxGetPr(Tempnx);
	mxArray *Tempny = mxGetField(BoundaryEdge, 0, "ny");
	double *ny = mxGetPr(Tempny);
	mxArray *TempLMass2d = mxGetField(BoundaryEdge, 0, "M");
	double  *LMass2d = mxGetPr(TempLMass2d);
	mxArray *TempFLAV = mxGetField(BoundaryEdge, 0, "LAV");
	double *FLAV = mxGetPr(TempFLAV);

	mxArray *TempFToF2d = mxGetField(BoundaryEdge2d, 0, "FToF");
	double *FToF2d = mxGetPr(TempFToF2d);
	mxArray *TempFToE2d = mxGetField(BoundaryEdge2d, 0, "FToE");
	double *FToE2d = mxGetPr(TempFToE2d);

	/*
	  This part can not parallized with OpemMP, since we may alter the the result at the same time
	  through different threads.
	*/
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
	for (int edge = 0; edge < BENe; edge++){

		double *FpIndex = malloc(Nfp*sizeof(double));

		for (int p = 0; p < Nfp; p++){
			FpIndex[p] = FToN1[Nfp*edge + p];
		}
		
		int LocalEle;
		LocalEle = (int)FToE[2 * edge];
		double *Tau = malloc(Nfp*sizeof(double));
		GetPenaltyParameter(Tau, LAV[LocalEle - 1], FLAV[edge], P, Nface, Nfp);
		
		double *TempEToE = NULL, *TempJ = NULL, *TempJs = NULL;
		TempEToE = EToE + (LocalEle - 1)*Nface;
		TempJ = J + (LocalEle - 1)*Np;
		TempJs = Js + edge * Nfp;

		ImposeDirichletBoundaryCondition(sr, OutRHS + (LocalEle - 1)*Np, irs, jcs, LocalEle, \
			Np, Nfp, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, \
			sy + (LocalEle - 1)*Np, Dr, Ds, Tau, nx + edge * Nfp, ny + edge * Nfp, \
			Mass3d, TempJ, TempJs, LMass2d, TempEToE, Nface, FpIndex, DirichDataValue + edge * Nfp);

		free(FpIndex);
		free(Tau);
	}


	/*
	for (int edge = 0; edge < BENe2d; edge++){
		int Flag = 0;
		int face2d = 0;
		int ele2d = 0;
		int Nface2d = Nface - 2;

		double *FpIndex = malloc(Nfp*sizeof(double));

		GetFaceTypeAndFaceOrder(&Flag, &face2d, &ele2d, FToF2d, FToE2d, ftype2d, edge);

		double *Tau = malloc(Nfp*sizeof(double));

		if (Flag == 1){//here 1 stand for Dirichlet boundary condition
			double *TempEToE = NULL, *TempJ = NULL, *TempJs = NULL;

			int GlobalFace, LocalEle;

			for (int L = 0; L < Nlayer; L++){

				LocalEle = (ele2d - 1)*Nlayer + L + 1;

				TempEToE = EToE + (ele2d - 1)*Nlayer*Nface + L*Nface;

				TempJ = J + (ele2d - 1)*Nlayer*Np + L*Np;

				GlobalFace = GetGlobalFace(face2d, BENe, FToE, FToF, LocalEle);

				for (int p = 0; p < Nfp; p++){
					FpIndex[p] = FToN1[Nfp*GlobalFace + p];
				}

				TempJs = Js + GlobalFace * Nfp;

				GetPenaltyParameter(Tau, LAV[LocalEle - 1], FLAV[GlobalFace], P, Nface, Nfp);

				ImposeDirichletBoundaryCondition(sr, OutRHS + (LocalEle - 1)*Np, irs, jcs, LocalEle, \
					Np, Nfp, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np,\
					sy + (LocalEle - 1)*Np, Dr, Ds, Tau, nx + GlobalFace * Nfp, ny + GlobalFace * Nfp,\
					Mass3d, TempJ, TempJs, LMass2d, TempEToE, Nface, FpIndex, DirichDataValue + GlobalFace * Nfp);

			}
		}
		free(FpIndex);
		free(Tau);
	}
	*/
}

void ImposeDirichletBoundaryCondition(double *dest, double *InputRHS, mwIndex *irs, mwIndex *jcs, int LocalEle, \
	int Np, int Nfp, double *rx, double *sx, double *ry, double *sy, double *Dr, double *Ds, double *Tau, \
	double *nx, double *ny, double *Mass3d, double *J, double *Js, double *Mass2d, double *EToE, \
	int Nface, double *FpIndex, double *DirichData){
	double *DxBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DxBuff, Dr, rx, Np);
	double *Dx = malloc(Np*Np*sizeof(double));
	DiagMultiply(Dx, Ds, sx, Np);
	Add(Dx, Dx, DxBuff, Np*Np);

	double *DyBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DyBuff, Dr, ry, Np);
	double *Dy = malloc(Np*Np*sizeof(double));
	DiagMultiply(Dy, Ds, sy, Np);
	Add(Dy, Dy, DyBuff, Np*Np);

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
	double *WeightedEleMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(WeightedEleMass2d, EleMass2d, DirichData, Nfp);
	double *TempRHS = malloc(Np*sizeof(double));
	memset(TempRHS, 0, Np*sizeof(double));
	double *TempRHSBuff = malloc(Np*sizeof(double));
	memset(TempRHSBuff, 0, Np*sizeof(double));
	double *DirichEdge2d = malloc(Nfp*sizeof(double));
	memset(DirichEdge2d, 0, Nfp*sizeof(double));
	double *DirichEdgeBuff = malloc(Nfp*Np*sizeof(double));

	/* For term $\int_{\partial \Omega^D}u_h\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, Dx, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, nx[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);
	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,\
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);
	/* For term $\int_{\partial \Omega^D}s\boldsymbol{n}\cdot \nabla_hu_h d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d,\
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/* For term $\int_{\partial \Omega^D}s\boldsymbol{n}\cdot\nabla_h u_hd\boldsymbol{x}$, y direction follows*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, Dy, FpIndex, Nfp, Np);
	/*The vector is a constant for a given face*/
	MultiplyByConstant(FacialDiffMatrix, FacialDiffMatrix, ny[0], Np*Nfp);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,\
		(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, y direction follows*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,\
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);

	/* For term $\int_{\partial \Omega^D}s\boldsymbol{n}\cdot\nabla_h u_hd\boldsymbol{x}$, y direction follows*/
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d,\
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, FpIndex, Np, Nfp);

	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(TempMass2d, EleMass2d, Tau, Nfp);
	/*For term $-\int_{\partial \Omega^d}\tau^k s u_hd\boldsymbol{x}$*/
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, FpIndex, FpIndex, Np, Nfp, -1.0);

	/*For term $-\int_{\partial \Omega^D}\tau^ksu_Dd\boldsymbol{x}$*/
	DiagMultiply(WeightedEleMass2d, WeightedEleMass2d, Tau, Nfp);

	SumInColumn(DirichEdge2d, WeightedEleMass2d, Nfp);

	MultiplyByConstant(DirichEdge2d, DirichEdge2d, -1, Nfp);

	AssembleDataIntoPoint(TempRHSBuff, DirichEdge2d, FpIndex, Nfp);

	/*Multiply the contribution by inverse matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvEleMass3d,\
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);
	/*Multiply the contribution due to Dirichlet boundary condition by inverse matrix*/
	ptrdiff_t Col = 1;
	MatrixMultiply("N", "N", (ptrdiff_t)Np, Col, (ptrdiff_t)Np, 1.0, InvEleMass3d,\
		(ptrdiff_t)Np, TempRHSBuff, (ptrdiff_t)Np, 0.0, TempRHS, (ptrdiff_t)Np);

	int UniNum = 0, StartPoint;
	/*Find the exact place where to fill in the data, and the place is stored in StartPoint*/
	double *TempEToE = malloc((Nface + 1)*sizeof(double));
	FindUniqueElementAndSortOrder(TempEToE, EToE, &UniNum, Nface - 2, LocalEle);
	int NonzeroPerColumn = jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np];
	for (int j = 0; j < UniNum; j++){
		if ((int)TempEToE[j] == LocalEle){
			StartPoint = (int)jcs[(LocalEle - 1)*Np] + j*Np;
			break;
		}	
	}

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	Add(InputRHS, InputRHS, TempRHS, Np);

	free(DxBuff);
	free(Dx);
	free(DyBuff);
	free(Dy);
	free(EleMass3d);
	free(EleMass2d);
	free(InvEleMass3d);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(WeightedEleMass2d);
	free(TempRHS);
	free(TempRHSBuff);
	free(DirichEdge2d);
	free(DirichEdgeBuff);
	free(TempEToE);
	free(TempMass2d);
}

int GetGlobalFace(int face, int Ne, double *FToE, double *FToF, int LocalElement){
	for (int i = 0; i < Ne; i++){
		if ((int)FToF[2 * i] == face && (int)FToE[2 * i] == LocalElement){
			return i;
			break;
		}
	}
	return -1; //failed
}

void AssembleDataIntoPoint(double *dest, double *source, double *Index, int Nfp){
	for (int i = 0; i < Nfp; i++){
		dest[(int)Index[i] - 1] += source[i];
	}
}

void SumInColumn(double *dest, double *Source, int Np){
	for (int col = 0; col < Np; col++){
		for (int Row = 0; Row < Np; Row++){
			dest[col] += Source[col*Np + Row];
		}
	}
}

void SumInRow(double *dest, double *Source, int Np, int ColNum){
	for (int row = 0; row < Np; row++){
		for (int col = 0; col < ColNum; col++){
			dest[row] += Source[col*Np + row];
		}
	}
}

void GetPenaltyParameter(double *dest, double LAV, double FLAV, double P, int Nface, int Nfp){
	for (int i = 0; i < Nfp; i++){
		dest[i] = (P + 1)*(P + 3) / 3.0 * Nface / 2.0 * FLAV / LAV;
	}
}

void GetFaceTypeAndFaceOrder(int *flag, int *face, int *ele, double *FToF, double *FToE, signed char *ftype, int edge){
	(*ele) = (int)FToE[2 * edge];
	(*face) = (int)FToF[2 * edge];
	/*At present, only Newmann bounary is considered*/
	if ((NdgEdgeType)ftype[edge] == NdgEdgeSlipWall){
		(*flag) = 1;
	}
	else{
		(*flag) = 0;
	}
}