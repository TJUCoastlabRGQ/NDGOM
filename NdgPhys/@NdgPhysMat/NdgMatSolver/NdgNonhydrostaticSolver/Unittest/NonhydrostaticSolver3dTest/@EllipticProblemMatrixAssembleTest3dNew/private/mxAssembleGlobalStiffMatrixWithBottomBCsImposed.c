#include "../../../../../../../NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver3d/private/SWENonhydrostatic3d.h"
#include <stdio.h>

void GetFaceTypeAndFaceOrder(int *, int *, int *, double *, double *, signed char *, int);

void ImposeDirichletBoundaryCondition(double *, double *, mwIndex *, mwIndex *, int, \
	int, int, double *, double *, double *, double *, double *, double *, double *, double *, \
	double, double *, double *, double *, double *, double *, double *, double *, \
	int, double *, double *, double *, double *, double *);

void ImposeNewmannBoundaryCondition(double *, int , int , int , double *, \
	double *, double *, double *, double *, double *);

void SumInColumn(double *, double *, int );

void SumInRow(double *, double *, int , int );

int GetGlobalFace(int , int , double *, double *, int );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *StiffMatrix = mxGetPr(prhs[0]);
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
	memcpy(sr, StiffMatrix, Tempjcs[col] * sizeof(double));
	memcpy(irs, Tempirs, Tempjcs[col] * sizeof(mwIndex));
	memcpy(jcs, Tempjcs, (col + 1)*sizeof(mwIndex));

	const mxArray *BoundaryEdge = prhs[3];
	const mxArray *cell = prhs[4];
	const mxArray *mesh = prhs[5];

	double *K13 = mxGetPr(prhs[6]);
	double *K23 = mxGetPr(prhs[7]);
	double *K33 = mxGetPr(prhs[8]);

    double *NewmannData = mxGetPr(prhs[9]);

    char* BoundaryType;
    BoundaryType = mxArrayToString(prhs[10]); 

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
	mxArray *Tempnz = mxGetField(BoundaryEdge, 0, "nz");
	double *nz = mxGetPr(Tempnz);
	mxArray *TempLMass2d = mxGetField(BoundaryEdge, 0, "M");
	double  *LMass2d = mxGetPr(TempLMass2d);
	mxArray *TempFLAV = mxGetField(BoundaryEdge, 0, "LAV");
	double *FLAV = mxGetPr(TempFLAV);

	double *Tau = malloc(BENe * sizeof(double));
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

		CalculatePenaltyParameter(Tau, FToE, FToN1, FToN1, Np, Nfp, \
			edge, K13, K23, K33, FLAV, LAV, max(N, Nz), Nface);
		
		double *TempJ = NULL, *TempJs = NULL;
		TempJ = J + (LocalEle - 1)*Np;
		TempJs = Js + edge * Nfp;

		if (!strcmp(BoundaryType, "Dirichlet")){
			ImposeDirichletBoundaryCondition(sr, OutRHS + (LocalEle - 1)*Np, irs, jcs, LocalEle, \
				Np, Nfp, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, \
				sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
				Dr, Ds, Dt, *(Tau+edge), nx + edge * Nfp, ny + edge * Nfp, nz + edge*Nfp, \
				Mass3d, TempJ, TempJs, LMass2d, Nface, FpIndex, DirichDataValue + edge * Nfp, \
				K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np);
		}
		else if (!strcmp(BoundaryType, "Newmann")){
			ImposeNewmannBoundaryCondition(OutRHS + (LocalEle - 1)*Np, LocalEle, Np, Nfp, Mass3d, \
				TempJ, TempJs, LMass2d, FpIndex, NewmannData + edge * Nfp);
		}

		free(FpIndex);
	}
	free(Tau);
}

void ImposeNewmannBoundaryCondition(double *InputRHS, int LocalEle, int Np, int Nfp, double *Mass3d, \
	double *J, double *Js, double *Mass2d, double *FpIndex, double *NewmannDataValue){

	double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(EleMass2d, Mass2d, Js, Nfp);

	double *TempRHSBuff = malloc(Np * 1 * sizeof(double));

	memset(TempRHSBuff, 0, Np * 1 * sizeof(double));

	double *TempFacialData = malloc(Nfp * 1 * sizeof(double));

	ptrdiff_t Col = 1;

	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, Col, (ptrdiff_t)Nfp, 1.0, EleMass2d, \
		(ptrdiff_t)Nfp, NewmannDataValue, (ptrdiff_t)Nfp, 0.0, TempFacialData, (ptrdiff_t)Nfp);

	AssembleDataIntoPoint(TempRHSBuff, TempFacialData, FpIndex, Nfp);

	MultiplyByConstant(TempRHSBuff, TempRHSBuff, -1.0, Np);

	Add(InputRHS, InputRHS, TempRHSBuff, Np);


	free(EleMass2d);

	free(TempRHSBuff);

	free(TempFacialData);

}

void ImposeDirichletBoundaryCondition(double *dest, double *InputRHS, mwIndex *irs, mwIndex *jcs, int LocalEle, \
	int Np, int Nfp, double *rx, double *sx, double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt,\
	double Tau, double *nx, double *ny, double *nz, double *Mass3d, double *J, double *Js, double *Mass2d, \
	int Nface, double *FpIndex, double *DirichData, double *K13, double *K23, double *K33){
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
	double *Contribution = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));
	double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
	double *WeightedEleMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(WeightedEleMass2d, EleMass2d, DirichData, Nfp);
	double *TempRHSBuff = malloc(Np*sizeof(double));
	memset(TempRHSBuff, 0, Np*sizeof(double));
	double *DirichEdge2d = malloc(Nfp*sizeof(double));
	memset(DirichEdge2d, 0, Nfp*sizeof(double));
	double *DirichEdgeBuff = malloc(Nfp*Np*sizeof(double));

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

	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix, \
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);

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

	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix, \
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);

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

	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix, \
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);

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

	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix, \
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);

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


	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix, \
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);

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

	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix, \
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);

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

	/*For term $\int_{\partial \Omega^D}u_D\nabla_h s\cdot\boldsymbol{n}d\boldsymbol{x}$, x direction first*/
	MatrixMultiply("T", "T", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix, \
		(ptrdiff_t)Nfp, WeightedEleMass2d, (ptrdiff_t)Nfp, 0.0, DirichEdgeBuff, (ptrdiff_t)Np);

	SumInRow(TempRHSBuff, DirichEdgeBuff, Np, Nfp);

	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));
	MultiplyByConstant(TempMass2d, EleMass2d, Tau, Nfp * Nfp);
	/*For term $-\int_{\partial \Omega^d}\tau^k s u_hd\boldsymbol{x}$*/
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, FpIndex, FpIndex, Np, Nfp, -1.0);

	/*For term $-\int_{\partial \Omega^D}\tau^ksu_Dd\boldsymbol{x}$*/
	MultiplyByConstant(WeightedEleMass2d, WeightedEleMass2d, Tau, Nfp * Nfp);

	SumInColumn(DirichEdge2d, WeightedEleMass2d, Nfp);

	MultiplyByConstant(DirichEdge2d, DirichEdge2d, -1.0, Nfp);

	AssembleDataIntoPoint(TempRHSBuff, DirichEdge2d, FpIndex, Nfp);

	double *SortedEid = malloc(Nfp*sizeof(double));
	memcpy(SortedEid, FpIndex, Nfp*sizeof(double));
	Sort(SortedEid, Nfp);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedEid, SortedEid, Np, Nfp, TempContribution, LocalEle, LocalEle);

	Add(InputRHS, InputRHS, TempRHSBuff, Np);

	free(SortedEid);
	free(DxBuff);
	free(Dx);
	free(DyBuff);
	free(Dy);
	free(Dz);
	free(EleMass2d);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(WeightedEleMass2d);
	free(TempRHSBuff);
	free(DirichEdge2d);
	free(DirichEdgeBuff);
	free(TempMass2d);
	free(TempDiffMatrix);
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
