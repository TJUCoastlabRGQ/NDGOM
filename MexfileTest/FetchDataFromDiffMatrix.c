#include "../NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver3d/private/SWENonhydrostatic3d.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//The first dimension of the sparse matrix
	double *Row = mxGetPr(prhs[0]);

	mxArray *TempJ3d = mxGetField(prhs[1], 0, "J");
	double *J3d = mxGetPr(TempJ3d);
	mxArray *Temprx3d = mxGetField(prhs[1], 0, "rx");
	double *rx3d = mxGetPr(Temprx3d);
	mxArray *Tempsx3d = mxGetField(prhs[1], 0, "sx");
	double *sx3d = mxGetPr(Tempsx3d);
	mxArray *Tempry3d = mxGetField(prhs[1], 0, "ry");
	double *ry3d = mxGetPr(Tempry3d);
	mxArray *Tempsy3d = mxGetField(prhs[1], 0, "sy");
	double *sy3d = mxGetPr(Tempsy3d);
	mxArray *Temptz3d = mxGetField(prhs[1], 0, "tz");
	double *tz3d = mxGetPr(Temptz3d);

	mxArray *TempDr3d = mxGetField(prhs[2], 0, "Dr");
	double *Dr3d = mxGetPr(TempDr3d);
	mxArray *TempDs3d = mxGetField(prhs[2], 0, "Ds");
	double *Ds3d = mxGetPr(TempDs3d);
	mxArray *TempDt3d = mxGetField(prhs[2], 0, "Dt");
	double *Dt3d = mxGetPr(TempDt3d);
	mxArray *TempNp3d = mxGetField(prhs[2], 0, "Np");
	int Np3d = (int)mxGetScalar(TempNp3d);
	mxArray *TempMass3d = mxGetField(prhs[2], 0, "M");
	double *mass3d = mxGetPr(TempMass3d);

	mxArray *TempJ2d = mxGetField(prhs[3], 0, "J");
	double *J2d = mxGetPr(TempJ2d);

    //The interpolation point for 2d master cell
	mxArray *TempNp2d = mxGetField(prhs[4], 0, "Np");
	int Np2d = (int)mxGetScalar(TempNp2d);
	mxArray *TempMass2d = mxGetField(prhs[4], 0, "M");
	double *mass2d = mxGetPr(TempMass2d);

	int Ele = (int)mxGetScalar(prhs[5]);

	double *DiffBuff = malloc(Np3d*Np3d*sizeof(double));
	/*Test case for $D_x$*/
	plhs[0] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *DiffMatrixX = mxGetPr(plhs[0]);

	DiagMultiply(DiffBuff, Dr3d, rx3d + (Ele-1)*Np3d, Np3d);
	DiagMultiply(DiffMatrixX, Ds3d, sx3d + (Ele - 1)*Np3d, Np3d);
	Add(DiffMatrixX, DiffMatrixX, DiffBuff, Np3d*Np3d);
	/*Test case for $D_y$*/
	plhs[1] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *DiffMatrixY = mxGetPr(plhs[1]);
	DiagMultiply(DiffBuff, Dr3d, ry3d + (Ele - 1)*Np3d, Np3d);
	DiagMultiply(DiffMatrixY, Ds3d, sy3d + (Ele - 1)*Np3d, Np3d);
	Add(DiffMatrixY, DiffMatrixY, DiffBuff, Np3d*Np3d);
	/*Test case for $D_z$*/
	plhs[2] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *DiffMatrixZ = mxGetPr(plhs[2]);
	DiagMultiply(DiffMatrixZ, Dt3d, tz3d + (Ele - 1)*Np3d, Np3d);
	/*Test case for fetching facial diff matrix D_{x2d}*/
	plhs[3] = mxCreateDoubleMatrix(Np2d, Np3d, mxREAL);
	double *Dx2d = mxGetPr(plhs[3]);
	AssembleFacialDiffMatrix(Dx2d, DiffMatrixX, Row, Np2d, Np3d);
	/*Test case for fetching facial diff matrix D_{y2d}*/
	plhs[4] = mxCreateDoubleMatrix(Np2d, Np3d, mxREAL);
	double *Dy2d = mxGetPr(plhs[4]);
	AssembleFacialDiffMatrix(Dy2d, DiffMatrixY, Row, Np2d, Np3d);
	/*Test case for three-dimensional mass matrix $M_{3d}$*/
	plhs[5] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *M3d = mxGetPr(plhs[5]);
	DiagMultiply(M3d, mass3d, J3d + (Ele - 1)*Np3d, Np3d);
	/*Test case for three-dimensional inverse mass matrix $M^{-1}_{3d}$*/
	plhs[6] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *InvM3d = mxGetPr(plhs[6]);
	memcpy(InvM3d, M3d, Np3d*Np3d*sizeof(double));
	MatrixInverse(InvM3d, (ptrdiff_t)Np3d);

	/*Test case for matrix multiplication $Dz^TM^k_{3d}$*/
	plhs[7] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *TempContributionBuff = mxGetPr(plhs[7]);
	MatrixMultiply("T", "N", (ptrdiff_t)Np3d, (ptrdiff_t)Np3d, (ptrdiff_t)Np3d, 1.0, DiffMatrixZ,
		(ptrdiff_t)Np3d, M3d, (ptrdiff_t)Np3d, 0.0, TempContributionBuff, (ptrdiff_t)Np3d);

	double *TempContribution = malloc(Np3d*Np3d*sizeof(double));
	/*Test case for matrix multiplication $-(M^k_{3d})^{-1}Dz^TM^k_{3d}Dz$*/
	plhs[8] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *Contribution = mxGetPr(plhs[8]);
	MatrixMultiply("N", "N", (ptrdiff_t)Np3d, (ptrdiff_t)Np3d, (ptrdiff_t)Np3d, 1.0, TempContributionBuff,
		(ptrdiff_t)Np3d, DiffMatrixZ, (ptrdiff_t)Np3d, 0.0, TempContribution, (ptrdiff_t)Np3d);
	MultiplyByConstant(TempContribution, TempContribution, -1, Np3d*Np3d);
	MatrixMultiply("N", "N", (ptrdiff_t)Np3d, (ptrdiff_t)Np3d, (ptrdiff_t)Np3d, 1.0, InvM3d,
		(ptrdiff_t)Np3d, TempContribution, (ptrdiff_t)Np3d, 0.0, Contribution, (ptrdiff_t)Np3d);

	plhs[9] = mxCreateDoubleMatrix(Np2d, Np2d, mxREAL);
	double *M2d = mxGetPr(plhs[9]);
	DiagMultiply(M2d, mass2d, J2d + (Ele - 1)*Np2d, Np2d);

	plhs[10] = mxCreateDoubleMatrix(Np3d, Np2d, mxREAL);
	double *FacialDiffxMultM2d = mxGetPr(plhs[10]);
	MatrixMultiply("T", "N", (ptrdiff_t)Np3d, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, Dx2d,
		(ptrdiff_t)Np2d, M2d, (ptrdiff_t)Np2d, 0.0, FacialDiffxMultM2d, (ptrdiff_t)Np3d);

	plhs[11] = mxCreateDoubleMatrix(Np3d, Np2d, mxREAL);
	double *FacialDiffyMultM2d = mxGetPr(plhs[11]);
	MatrixMultiply("T", "N", (ptrdiff_t)Np3d, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, Dy2d,
		(ptrdiff_t)Np2d, M2d, (ptrdiff_t)Np2d, 0.0, FacialDiffyMultM2d, (ptrdiff_t)Np3d);

	plhs[12] = mxCreateDoubleMatrix(Np2d, Np3d, mxREAL);
	double *M2dMultFacialDiffx = mxGetPr(plhs[12]);
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np3d, (ptrdiff_t)Np2d, 1.0, M2d,
		(ptrdiff_t)Np2d, Dx2d, (ptrdiff_t)Np2d, 0.0, M2dMultFacialDiffx, (ptrdiff_t)Np2d);

	plhs[13] = mxCreateDoubleMatrix(Np2d, Np3d, mxREAL);
	double *M2dMultFacialDiffy = mxGetPr(plhs[13]);
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np3d, (ptrdiff_t)Np2d, 1.0, M2d,
		(ptrdiff_t)Np2d, Dy2d, (ptrdiff_t)Np2d, 0.0, M2dMultFacialDiffy, (ptrdiff_t)Np2d);

	plhs[14] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *TestTempContributionx = mxGetPr(plhs[14]);
	AssembleContributionIntoRow(TestTempContributionx, M2dMultFacialDiffx, Row, Np3d, Np2d);

	plhs[15] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *TestTempContributiony = mxGetPr(plhs[15]);
	AssembleContributionIntoRow(TestTempContributiony, M2dMultFacialDiffy, Row, Np3d, Np2d);

	plhs[16] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *DirichletBuff2 = mxGetPr(plhs[16]);
	AssembleContributionIntoRowAndColumn(DirichletBuff2, M2d, Row, Row, Np3d, Np2d, -1);

	plhs[17] = mxCreateDoubleMatrix(Np3d, Np3d, mxREAL);
	double *TempContributionBuff2 = mxGetPr(plhs[17]);
	double *InnerEdgeContribution = malloc(Np3d*Np2d*sizeof(double));
	MatrixMultiply("T", "N", (ptrdiff_t)Np3d, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, Dx2d,
		(ptrdiff_t)Np2d, M2d, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np3d);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5, Np3d*Np2d);
	AssembleContributionIntoColumn(TempContributionBuff2, InnerEdgeContribution, Row, Np3d, Np2d);

	free(DiffBuff);
	free(TempContribution);
	free(InnerEdgeContribution);
}