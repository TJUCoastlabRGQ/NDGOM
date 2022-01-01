#if !defined(_WIN32)
#define dgesv dgesv_
#endif

#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgMemory.h"

#include <stdio.h>

#include <stdlib.h>

#include <string.h>

extern double *Tau, *u2d, *v2d;

extern char *VertDiffInitialized;

void FetchBoundaryData(double *dest, double *source, const int Np2d, double *Eid)
{
    for (int i = 0; i < Np2d; i++){
        dest[i] = source[(int)(Eid[i]) - 1];
    }
}

void CalculatePenaltyParameter(double *dest, const int Np2d, const int Np3d, double *UpEidM, double *BotEidM, double *nv, \
        int Nz, int P, int Nface)
{
    double nvM[Np2d], nvP[Np2d];
    //Note: the penalty parameter for the topmost face of each column is not needed, so we leave them undefined
    for (int Layer = 1; Layer < Nz; Layer++)
    {
        FetchBoundaryData(nvM, nv + (Layer-1)*Np3d, Np2d, BotEidM);
        FetchBoundaryData(nvP, nv + Layer*Np3d, Np2d, UpEidM);
        for (int p = 0; p < Np2d; p++)
        {
            dest[Layer*Np2d + p] = ((P + 1)*(P + 3) / 3.0)*(Nface / 2.0)*Nz*max(nvM[p], nvP[p]);
        }
    }
    FetchBoundaryData(nvM, nv + (Nz - 1)*Np3d, Np2d, BotEidM);
    for (int p = 0; p < Np2d; p++)
    {
        dest[Nz*Np2d + p] = ((P + 1)*(P + 3) / 3.0)*(Nface / 2.0)*Nz*nvM[p];
    }
}

/*Note: This function is used to assemble the local columns index and local rows index of the studied cell, and has been verified*/
void GetLocalRowsAndColumns(double *LocalRows, double *LocalColumns, int Index, int Np)
{
    for (int i = 0; i < Np; i++){
        LocalRows[i] = Index * Np + i + 1;
        LocalColumns[i] = Index * Np + i + 1;
    }
}
/*Note: This function is used to calculate the volumn integral contained in primal form*/
void VolumnIntegral(double *dest, double *Dz, double *Mass3d, double *diff, ptrdiff_t Np)
{
    double zero = 0.0;
    double Alpha1 = -1.0;
    double Alpha2 = 1.0;
    double *Tempdest = malloc(Np*Np*sizeof(double));
    dgemm("t", "n", &Np, &Np, &Np, &Alpha1, Dz, &Np, Mass3d, &Np, &zero, Tempdest, &Np);
    dgemm("n", "n", &Np, &Np, &Np, &Alpha2, Tempdest, &Np, diff, &Np, &zero, dest, &Np);
    free(Tempdest);
}

/*This function is used to impose the Newmann boundary at the surface and bottom boundary*/
void ImposeNewmannBoundary(double *Eid, double *mass2d, double* InvMassMatrix3d, double dt, double Imparam, \
        const double* BoundNewmannData, ptrdiff_t Np2d, int K2d, double *SystemRHS, ptrdiff_t Np3d, int K3d, double *StiffMatrix, int Nvar)
{
    double *tempRHS = malloc(Np3d*sizeof(double));
    double *RHS2d = malloc(Np2d*sizeof(double));
    for (int i = 0; i < Nvar; i++){
        memset(tempRHS, 0, Np3d*sizeof(double));
        memset(RHS2d, 0, Np2d*sizeof(double));
        ptrdiff_t colB = 1;
        double Alpha = 1.0;
        double Beta = 0.0;
        dgemm("n", "n", &Np2d, &colB, &Np2d, &Alpha, mass2d, &Np2d, BoundNewmannData + i*Np2d*K2d, &Np2d, &Beta, RHS2d, &Np2d);
        for (int j = 0; j < Np2d; j++)
            tempRHS[(int)Eid[j] - 1] = RHS2d[j];
        dgemm("n", "n", &Np3d, &colB, &Np3d, &Alpha, InvMassMatrix3d, &Np3d, tempRHS, &Np3d, &Beta, StiffMatrix + i*Np3d, &Np3d);
        for (int j = 0; j < Np3d; j++)
            SystemRHS[i*K3d*Np3d + j] += dt*Imparam*(*(StiffMatrix + i*Np3d + j));
    }
    free(tempRHS);
    free(RHS2d);
}
/*Note: This function is used to get the row index of cell bellow the studied one*/
void GetBottomRows(double *BottomRows, double *LocalRows, int Np3d)
{
    for (int i = 0; i < Np3d; i++)
        BottomRows[i] = LocalRows[i] + Np3d;
}
/*Note: This function is used to get the row index of cell on top of the studied one*/
void GetUpRows(double *UpRows, double *LocalRows, int Np3d)
{
    for (int i = 0; i < Np3d; i++)
        UpRows[i] = LocalRows[i] - Np3d;
}
/*Note: This function is used to get the facial diff matrix, and has been checked*/
void AssembleFacialDiffMatrix(double *dest, double *source, double *Eid, int Np2d, int Np3d)
{
    int Index = 0;
    for (int colI = 0; colI < Np3d; colI++){
        for (int RowI = 0; RowI < Np2d; RowI++)
        {
            dest[Index] = source[colI*Np3d + (int)Eid[RowI] - 1];
            Index++;
        }
    }
}

/*Note: This function is used to assemble the global stiff matrix*/
void AssembleGlobalStiffMatrix(double *dest, double *invMass, double *OP, double *rows, double *columns, double alpha, int StiffRow, ptrdiff_t Np)
{
    double *TempDest = malloc(Np * Np * sizeof(double));
    double Beta = 0;
    dgemm("n", "n", &Np, &Np, &Np, &alpha, invMass, &Np, OP, &Np, &Beta, TempDest, &Np);
    //here Np stands Np3d, StiffRow is the row of the stiffmatrix
    AssembleContributionIntoRowAndColumn(dest, TempDest, rows, columns, StiffRow, (int)Np, 1);
    free(TempDest);
}
/*Note: This function is used to solve the equation AX=b to get term X, and has been checked*/
void EquationSolve(double *dest, ptrdiff_t m, ptrdiff_t n, double *CoeMatrix, double dt, double ImplicitParam, int Nvar, int K3d, int Np)
{
    double *TempCoeMatrix = malloc(m*m*sizeof(double));
    ptrdiff_t *iPivot = malloc(m*sizeof(ptrdiff_t));
    ptrdiff_t info;
    ptrdiff_t p = m;
    for (int var = 0; var < Nvar; var++){
        memcpy(TempCoeMatrix, CoeMatrix + var*m*m, m*m*sizeof(double));
        for (int col = 0; col < m; col++){
            for (int row = 0; row < m; row++){
                TempCoeMatrix[col*m + row] = -1 * dt*ImplicitParam*TempCoeMatrix[col*m + row];
            }
            TempCoeMatrix[col*m + col] += 1;
        }
        dgesv(&m, &n, TempCoeMatrix, &m, iPivot, dest + var*K3d*Np, &p, &info);
    }
    free(iPivot);
    free(TempCoeMatrix);
}
/*Note: This function is used to calculate the local boundary contribution to the local stiff operator OP11*/
void LocalBoundaryIntegral(double *eid, double *DiffMatrix, double *mass2d, double *TempTau, double *OP11, ptrdiff_t Np3d, ptrdiff_t Np2d, int Flag, double epsilon)
{
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, DiffMatrix, eid, (int)Np2d, (int)Np3d);
    double Alpha = -1 * epsilon*Flag*0.5, Beta = 0.0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));
    dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
    AssembleContributionIntoColumn(OP11, EdgeContribution, eid, (int)Np3d, (int)Np2d);
    Alpha = 0.5*Flag;
    dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
    AssembleContributionIntoRow(OP11, EdgeContribution, eid, (int)Np3d, (int)Np2d);
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
    DiagMultiply(DoubleJump, mass2d, TempTau, (int)Np2d);
    AssembleContributionIntoRowAndColumn(OP11, DoubleJump, eid, eid, (int)Np3d, (int)Np2d, -1);
    free(FDiffMatrix);
    free(EdgeContribution);
    free(DoubleJump);
}
/*Note: This function is used to impose the homogeneous Dirichlet boundary condition at the bottom boundary*/
void ImposeDirichletBoundary(double *eid, double *DiffMatrix, double *mass2d, double *TempTau, double *OP11, ptrdiff_t Np3d, ptrdiff_t Np2d, int Flag, const double epsilon)
{
    /*	double *Tau = malloc(Np2d*sizeof(double));
    for (int i = 0; i < Np2d; i++)
        Tau[i] = 2 * TempTau[i];
     */
    //LocalBoundaryIntegral(eid, DiffMatrix, mass2d, Tau, OP11, Np3d, Np2d, Flag, epsilon);
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, DiffMatrix, eid, (int)Np2d, (int)Np3d);
    double Alpha = -1 * epsilon*Flag*0.5*2, Beta = 0.0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));
    dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
    AssembleContributionIntoColumn(OP11, EdgeContribution, eid, (int)Np3d, (int)Np2d);
    Alpha = 0.5*Flag*2;
    dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
    AssembleContributionIntoRow(OP11, EdgeContribution, eid, (int)Np3d, (int)Np2d);
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
    //DiagMultiply(DoubleJump, mass2d, Tau, (int)Np2d);
    DiagMultiply(DoubleJump, mass2d, TempTau, (int)Np2d);
    AssembleContributionIntoRowAndColumn(OP11, DoubleJump, eid, eid, (int)Np3d, (int)Np2d, -1);
    free(FDiffMatrix);
    free(EdgeContribution);
    free(DoubleJump);
//	free(Tau);
}
/*Note: This function is used to calculate the adjacent boundary contribution to the adjacent stiff operator OP12, here adjacent means the test function is defined over the adjacent cell and not the local one*/
void AdjacentBoundaryIntegral(double *eidM, double *eidP, double *LocalDiff, double *AdjacentDiff, double *mass2d, double *TempTau, double *OP12, ptrdiff_t Np3d, \
        ptrdiff_t Np2d, int Flag, const double epsilon)
{
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, AdjacentDiff, eidP, (int)Np2d, (int)Np3d);
    double Alpha = -1 * epsilon*Flag*0.5, Beta = 0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));
    dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
    AssembleContributionIntoColumn(OP12, EdgeContribution, eidM, (int)Np3d, (int)Np2d);
    //Alpha = 0.5*Flag;
    Alpha = -1*0.5*Flag;
    AssembleFacialDiffMatrix(FDiffMatrix, LocalDiff, eidM, (int)Np2d, (int)Np3d);
    dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
    AssembleContributionIntoRow(OP12, EdgeContribution, eidP, (int)Np3d, (int)Np2d);
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
    DiagMultiply(DoubleJump, mass2d, TempTau, (int)Np2d);
    AssembleContributionIntoRowAndColumn(OP12, DoubleJump, eidP, eidM, (int)Np3d, (int)Np2d, 1);
    free(FDiffMatrix);
    free(EdgeContribution);
    free(DoubleJump);
}
/*Note: This function is used to add the Newmman boundary part back to the implicit right hand side*/
void AssembleBoundaryContribution(double *dest, double *source, int Np, int K3d, int Nvar)
{
    for (int var = 0; var < Nvar; var++){
        for (int i = 0; i < Np; i++)
            dest[var*K3d*Np + i] += source[var*Np+i];
    }
}

void CalculateVelocityAtBottomCellCenter(double *ucenter, double *vcenter, double *hubot, double *hvbot, \
	ptrdiff_t *RowVCV, ptrdiff_t *ColVCV, double *hcrit, double *h2d, double *VCV) {
	char *chn = "N";
	double alpha = 1.0;
	double beta = 0.0;
	ptrdiff_t Col = 1;
	dgemm(chn, chn, RowVCV, &Col, ColVCV, &alpha, VCV, RowVCV, hubot, ColVCV, &beta, ucenter, RowVCV);
	dgemm(chn, chn, RowVCV, &Col, ColVCV, &alpha, VCV, RowVCV, hvbot, ColVCV, &beta, vcenter, RowVCV);
	for (int p = 0; p < (int)(*RowVCV); p++) {
		if (h2d[p] >= (*hcrit)) {
			ucenter[p] = ucenter[p] / h2d[p];
			vcenter[p] = vcenter[p] / h2d[p];
		}
		else {
			ucenter[p] = 0.0;
			vcenter[p] = 0.0;
		}
	}
}

void ImposeImplicitNeumannBoundary(double *dest, double *EidM, double *Cf, double dt, double ImplicitParam, double *Depth, \
	double *EleMass2d, double *u2d, double *v2d, int Np2d, int Np3d, double hcrit) {
	double *Coe = malloc(Np2d * sizeof(double));
	double *TempOP11 = malloc(Np2d*Np2d * sizeof(double));
	memset(TempOP11, 0, Np2d*Np2d * sizeof(double));
	for (int p = 0; p < Np2d; p++) {
		if (Depth[p] >= hcrit) {
			Coe[p] = -1.0*dt*ImplicitParam*Cf[p] * sqrt(u2d[p] * u2d[p] + v2d[p] * v2d[p]) / Depth[p];
		}
		else {
			Coe[p] = 0.0;
		}
	}
	DiagRightMultiply(TempOP11, EleMass2d, Coe, Np2d);

	AssembleContributionIntoRowAndColumn(dest, TempOP11, EidM, EidM, Np3d, Np2d, -1);

	free(Coe);

	free(TempOP11);
}

void GetImplicitBoundaryContribution(double *dest, double *Cf, double *u2d, double *v2d, double *Height, \
	double *hu3d, double *hv3d, double *VCV, ptrdiff_t RowVCV, \
	ptrdiff_t ColVCV, int Np, int Np2d, double *Mass2d, double *InvMass3d, double *EidM, double hcrit) {

	double *TempRHS = malloc(Np * 2 * sizeof(double));
	double *Drag2d = malloc(Np2d * 2 * sizeof(double));
	double *TempRHS2d = malloc(Np2d * 2 * sizeof(double));
	double *ucenter = malloc(Np2d * sizeof(double));
	double *vcenter = malloc(Np2d * sizeof(double));

	CalculateVelocityAtBottomCellCenter(ucenter, vcenter, hu3d, hv3d, \
		&RowVCV, &ColVCV, &hcrit, Height, VCV);

	for (int p = 0; p < Np2d; p++) {
		Drag2d[p] = Cf[p] * sqrt(pow(u2d[p], 2.0) + pow(v2d[p], 2.0))*ucenter[p];
		Drag2d[Np2d + p] = Cf[p] * sqrt(pow(u2d[p], 2.0) + pow(v2d[p], 2.0))*vcenter[p];
	}

	char *chn = "N";
	double alpha = 1.0;
	double beta = 0.0;
	ptrdiff_t Col = 1;
	ptrdiff_t RowMass2d = (ptrdiff_t)Np2d;
	ptrdiff_t ColMass2d = (ptrdiff_t)Np2d;

	ptrdiff_t RowInvMass3d = (ptrdiff_t)Np;
	ptrdiff_t ColInvMass3d = (ptrdiff_t)Np;

	dgemm(chn, chn, &RowMass2d, &Col, &RowMass2d, &alpha, Mass2d, \
		&RowMass2d, Drag2d, &RowMass2d, &beta, TempRHS2d, &RowMass2d);

	dgemm(chn, chn, &RowMass2d, &Col, &RowMass2d, &alpha, Mass2d, \
		&RowMass2d, Drag2d + Np2d, &RowMass2d, &beta, TempRHS2d + Np2d, &RowMass2d);

	AssembleDataIntoPoint(TempRHS, TempRHS2d, EidM, Np2d);

	AssembleDataIntoPoint(TempRHS + Np, TempRHS2d + Np2d, EidM, Np2d);

	dgemm(chn, chn, &RowInvMass3d, &Col, &RowInvMass3d, &alpha, InvMass3d, \
		&RowInvMass3d, TempRHS, &RowInvMass3d, &beta, dest, &RowInvMass3d);

	dgemm(chn, chn, &RowInvMass3d, &Col, &RowInvMass3d, &alpha, InvMass3d, \
		&RowInvMass3d, TempRHS + Np, &RowInvMass3d, &beta, dest + Np, &RowInvMass3d);

	free(TempRHS);
	free(Drag2d);
	free(TempRHS2d);
	free(ucenter);
	free(vcenter);
}

void MyExit()
{
    if (!strcmp("True", VertDiffInitialized)){
        VertDiffMemoryDeAllocation();
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mexAtExit(&MyExit);
    const double *J2d = mxGetPr(prhs[0]);
    const double *J3d = mxGetPr(prhs[1]);
    const double *M2d = mxGetPr(prhs[2]);
    const double *M3d = mxGetPr(prhs[3]);
    const double *tz = mxGetPr(prhs[4]);
    const double *Dt = mxGetPr(prhs[5]);
    double *Diff = mxGetPr(prhs[6]);
    const double *surf = mxGetPr(prhs[7]);
    const double *bot = mxGetPr(prhs[8]);
    const double dt = mxGetScalar(prhs[9]);
    const double ImplicitParam = mxGetScalar(prhs[10]);
    const double *RHS = mxGetPr(prhs[11]);
    const double Prantl = mxGetScalar(prhs[12]);
    const int K2d = (int)mxGetN(prhs[0]);
    const int K3d = (int)mxGetN(prhs[1]);
    const int Nz = K3d / K2d;
    const int Np = (int)mxGetM(prhs[3]);
    const int Np2d = (int)mxGetM(prhs[2]);
    double *UpEidM = mxGetPr(prhs[13]);
    double *BotEidM = mxGetPr(prhs[14]);
    const int P = (int)mxGetScalar(prhs[15]);
    const int Nface = (int)mxGetScalar(prhs[16]);

    char* BoundaryType;
    BoundaryType = mxArrayToString(prhs[17]);
	// Add variables here. Added on 20211231 by RGQ.
	double *huv3d = mxGetPr(prhs[18]);
	double *hu3d = huv3d;
	double *hv3d = huv3d + Np*K3d;
	double *h2d = mxGetPr(prhs[19]);
	double hcrit = mxGetScalar(prhs[20]);
	double *VCV = mxGetPr(prhs[21]);
	double *Cf = mxGetPr(prhs[22]);
//    printf("%s\n", BoundaryType);
    
    mwSize DimOfRHS = mxGetNumberOfDimensions(prhs[11]);
    int Nvar;
    const size_t *PRHS;
    PRHS = mxGetDimensions(prhs[11]);
    if (DimOfRHS == 2)
        Nvar = 1;
    else{
        Nvar = (int)PRHS[2];
    }
//	double *TempRHS = malloc(Np*K3d*Nvar*sizeof(double));
//	memcpy(TempRHS, RHS, Np*K3d*Nvar*sizeof(double));
    
    const size_t NdimOut = 3;
    const mwSize dimOut[3] = {Np,K3d,Nvar};
    plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

    double *fphys = mxGetPr(plhs[0]);
    memcpy(fphys, RHS, Np*K3d*Nvar*sizeof(double));
    double *ImplicitRHS = mxGetPr(plhs[1]);
    /*Here, epsilon is set to be -1, and this correspoinding to the SIPG. This value can be set to 1(IIPG or NIPG) or 0(IIPG or NIPG)*/
    const double epsilon = -1;
    
    /*If not initialized, initialize first*/
    if (!strcmp("False", VertDiffInitialized))
    {
        VertDiffMemoryAllocation(Np2d, K2d, Nz);
    }
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int i = 0; i < K2d; i++){
        CalculatePenaltyParameter(Tau + i*Np2d*(Nz + 1), Np2d, Np, UpEidM, BotEidM, Diff + i*Np*Nz, Nz, P, Nface);
    }

	ptrdiff_t RowVCV = (ptrdiff_t)Np2d;
	ptrdiff_t ColVCV = (ptrdiff_t)Np;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++) {
		CalculateVelocityAtBottomCellCenter(u2d + i*Np2d, v2d + i*Np2d, hu3d + (i + 1)*Nz*Np - Np, hv3d + (i+1)*Nz*Np - Np, \
			&RowVCV, &ColVCV, &hcrit, h2d + i*Np2d, VCV);
	}

	for (int i = 0; i < K2d; i++) {
		CalculatePenaltyParameter(Tau + i*Np2d*(Nz + 1), Np2d, Np, UpEidM, BotEidM, Diff + i*Np*Nz, Nz, P, Nface);
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    
    for (int i = 0; i < K2d; i++){
        double *StiffMatrix = malloc(Np*Nz*Np*Nz*Nvar*sizeof(double));
        memset(StiffMatrix, 0, Np*Nz*Np*Nz*Nvar*sizeof(double));
        double *EleMass3d = malloc(Np*Np*sizeof(double));
        double *InvEleMass3d = malloc(Np*Np*sizeof(double));
        double *EleMass2d = malloc(Np2d*Np2d*sizeof(double));
        double *LocalPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
        double *Dz = malloc(Np*Np*sizeof(double));
        double *OP11 = malloc(Np*Np*sizeof(double));
        double *LocalRows = malloc(Np*sizeof(double));
        double *LocalColumns = malloc(Np*sizeof(double));
        DiagMultiply(EleMass3d, M3d, J3d + i*Nz*Np, Np);
        memcpy(InvEleMass3d, EleMass3d, Np*Np*sizeof(double));
        MatrixInverse(InvEleMass3d, (ptrdiff_t)Np);
        DiagMultiply(EleMass2d, M2d, J2d + i*Np2d, Np2d);
        DiagMultiply(Dz, Dt, tz + i*Nz*Np, Np);
        DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np, Np);
        GetLocalRowsAndColumns(LocalRows, LocalColumns, 0, Np);
        /*Calculate the volumn integral for the first cell, and impose the surface Newmann boundary condition*/
        VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
        double *SurfBoundStiffTerm = malloc(Np*Nvar*sizeof(double));
		memset(SurfBoundStiffTerm, 0, Np*Nvar * sizeof(double));
        double *BotBoundStiffTerm = malloc(Np*Nvar*sizeof(double));
		memset(BotBoundStiffTerm, 0, Np*Nvar * sizeof(double));
        ImposeNewmannBoundary(UpEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, surf + i*Np2d, (ptrdiff_t)Np2d, K2d, fphys + i*Np*Nz, (ptrdiff_t)Np, K3d, SurfBoundStiffTerm, Nvar);
        if (Nz != 1){
            /*When vertical layers greater than one, we calculate the bottom boundary integral part */
            double *BottomAdjacentRows = malloc(Np*sizeof(double));
            double *UpAdjacentRows = malloc(Np*sizeof(double));
            GetBottomRows(BottomAdjacentRows, LocalRows, Np);
            double *BottomPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
            double *UpPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
            DiagMultiply(BottomPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + Np, Np);
            LocalBoundaryIntegral(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + 2 - 1), OP11, (ptrdiff_t)Np, (ptrdiff_t)Np2d, -1, epsilon);
            double *OP12 = malloc(Np*Np*sizeof(double));
            memset(OP12, 0, Np*Np*sizeof(double));
            AdjacentBoundaryIntegral(BotEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + 2 - 1), OP12, Np, Np2d, -1, epsilon);
            
            for (int var = 0; var < 2; var++){
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, 1.0, Np*Nz, Np);
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP12, BottomAdjacentRows, LocalColumns, 1.0, Np*Nz, Np);
            }
            for (int var = 2; var < Nvar; var++){
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, Prantl, Np*Nz, Np);
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP12, BottomAdjacentRows, LocalColumns, Prantl, Np*Nz, Np);
            }
            for (int j = 1; j < Nz-1 ; j++){
                /*Calculate both the volume and surface integral for the second to the last but one cell*/
                GetLocalRowsAndColumns(LocalRows, LocalColumns, j, Np);
                GetBottomRows(BottomAdjacentRows, LocalRows, Np);
                GetUpRows(UpAdjacentRows, LocalRows, Np);
                DiagMultiply(UpPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (j-1)*Np, Np);
                DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + j*Np, Np);
                DiagMultiply(BottomPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (j + 1)*Np, Np);
                VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
                LocalBoundaryIntegral(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + j+1), OP11, (ptrdiff_t)Np, (ptrdiff_t)Np2d, -1, epsilon);
                LocalBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + j), OP11, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1, epsilon);
                memset(OP12, 0, Np*Np*sizeof(double));
                AdjacentBoundaryIntegral(UpEidM, BotEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + j), OP12, Np, Np2d, 1, epsilon);
                
                for (int var = 0; var < 2; var++){
                    AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, 1, Np*Nz, Np);
                    AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP12, UpAdjacentRows, LocalColumns, 1, Np*Nz, Np);
                }
                for (int var = 2; var < Nvar; var++){
                    AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, Prantl, Np*Nz, Np);
                    AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP12, UpAdjacentRows, LocalColumns, Prantl, Np*Nz, Np);
                }
                memset(OP12, 0, Np*Np*sizeof(double));
                AdjacentBoundaryIntegral(BotEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + j + 1), OP12, Np, Np2d, -1, epsilon);
                for (int var = 0; var < 2 && var < Nvar; var++){
                    AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP12, BottomAdjacentRows, LocalColumns, 1, Np*Nz, Np);
                }
                for (int var = 2; var < Nvar; var++){
                    AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP12, BottomAdjacentRows, LocalColumns, Prantl, Np*Nz, Np);
                }
            }
            /*For the bottom most cell, calculate the volumn integral and impose Neumann boundary condition or the Dirichlet boundary condition*/
            GetLocalRowsAndColumns(LocalRows, LocalColumns, Nz - 1, Np);
            GetUpRows(UpAdjacentRows, LocalRows, Np);
            DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (Nz-1)*Np, Np);
            DiagMultiply(UpPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (Nz - 1 - 1)*Np, Np);
            VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
            LocalBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + Nz-1), OP11, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1, epsilon);
            /*For passive transport substances, we only impose Neumann boundary condition, so this has no effect on the stiff matrix*/
            for (int var = 2; var < Nvar; var++){
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, Prantl, Np*Nz, Np);
            }
            /*We note that, the Neumann data is zero by default, so the impositon of Neumann BC for hu and hv will not affect the Dirichlet boundary condition for hu and hv*/
            ImposeNewmannBoundary(BotEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, bot + i*Np2d, (ptrdiff_t)Np2d, K2d, fphys + i*Np*Nz + (Nz - 1)*Np, (ptrdiff_t)Np, K3d, BotBoundStiffTerm, Nvar);
			// Impose implicit Neumann boundary condition here, this part is added on 20211231
			ImposeImplicitNeumannBoundary(OP11, BotEidM, Cf + i*Np2d, dt, ImplicitParam, h2d + i*Np2d, EleMass2d, u2d + i*Np2d, v2d + i*Np2d, Np2d, Np, hcrit);
            /*The following is used to add homogeneous dirichlet boundary for hu and hv*/
            if (!strcmp(BoundaryType, "Dirichlet")){
                ImposeDirichletBoundary(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + Nz), OP11, (ptrdiff_t)Np, (ptrdiff_t)Np2d, -1, epsilon);
            }
            memset(OP12, 0, Np*Np*sizeof(double));
            AdjacentBoundaryIntegral(UpEidM, BotEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + Nz - 1), OP12, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, epsilon);
            for (int var = 0; var < 2; var++){
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP12, UpAdjacentRows, LocalColumns, 1.0, Np*Nz, Np);
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, 1.0, Np*Nz, Np);
            }
            for (int var = 2; var < Nvar; var++){
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP12, UpAdjacentRows, LocalColumns, Prantl, Np*Nz, Np);
            }
            free(BottomAdjacentRows);
            free(UpAdjacentRows);
            free(BottomPhysicalDiffMatrix);
            free(UpPhysicalDiffMatrix);
            free(OP12);
        }
        else{/*If only one layer included in vertical direction, we need to consider the bottom face, and the Upper face has been considered in the first part, from 291-312*/
            /*For passive transport substances, we only impose Neumann boundary condition*/
            for (int var = 2; var < Nvar; var++){
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, Prantl, Np*Nz, Np);
            }
            /*We note that, the Neumann data is zero by default, so the impositon of Neumann BC for hu and hv will not affect the Dirichlet boundary condition for hu and hv*/
            ImposeNewmannBoundary(BotEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, bot + i*Np2d, (ptrdiff_t)Np2d, K2d, fphys + i*Np*Nz + (Nz - 1)*Np, (ptrdiff_t)Np, K3d, BotBoundStiffTerm, Nvar);
			// Impose implicit Neumann boundary condition here, this part is added on 20211231
			ImposeImplicitNeumannBoundary(OP11, BotEidM, Cf + i*Np2d, dt, ImplicitParam, h2d + i*Np2d, EleMass2d, u2d + i*Np2d, v2d + i*Np2d, Np2d, Np, hcrit);
			/* The following is used to add homogeneous dirichlet boundary for hu and hv*/
            if (!strcmp(BoundaryType, "Dirichlet")){
                ImposeDirichletBoundary(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d*(i*(Nz + 1) + Nz), OP11, (ptrdiff_t)Np, (ptrdiff_t)Np2d, -1, epsilon);
            }
            for (int var = 0; var < 2; var++){
                AssembleGlobalStiffMatrix(StiffMatrix + var*Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, 1, Np*Nz, Np);
            }
        }
        
        EquationSolve(fphys + i*Nz*Np, Nz*Np, 1, StiffMatrix, dt, ImplicitParam, Nvar, K3d, Np);
        
        const ptrdiff_t dimension = Np*Nz;
        const ptrdiff_t colB = 1;
        double Alpha = 1.0, Beta = 0.0;

        for (int var = 0; var < Nvar; var++){
            dgemm("n", "n", &dimension, &colB, &dimension, &Alpha, StiffMatrix + var*Np*Nz*Np*Nz, &dimension, fphys + var*K3d*Np + i*Nz*Np, &dimension, &Beta, ImplicitRHS + var*K3d*Np + i*Nz*Np, &dimension);
        }
        AssembleBoundaryContribution(ImplicitRHS  + i*Nz*Np, SurfBoundStiffTerm, Np, K3d, Nvar);
		// BotBoundStiffTerm with index 1 and 2 stands for hu and hv respectively, they are zero since they are treated implicitly.
		// We add them to the ImplicitRHS since they have no effect on the final result. 
        AssembleBoundaryContribution(ImplicitRHS  + i*Nz*Np + (Nz - 1)*Np, BotBoundStiffTerm, Np, K3d, Nvar);
		// We first need to calculate the contribution to the right hand side due to the implicit bottom condition, the add they to the RHS
		GetImplicitBoundaryContribution(BotBoundStiffTerm, Cf + i*Np2d, u2d + i*Np2d, v2d + i*Np2d, h2d + i*Np2d, \
			fphys + (i + 1)*Np*Nz - Np, fphys + Np*K3d + (i + 1)*Np*Nz - Np, VCV, RowVCV, ColVCV, Np, Np2d, EleMass2d, InvEleMass3d, \
			BotEidM, hcrit);
		//We add the right hand side due to the implicit bottom friction term back, only hu and hv considered. Added on 20211231 by RGQ 
		AssembleBoundaryContribution(ImplicitRHS + i*Nz*Np + (Nz - 1)*Np, BotBoundStiffTerm, Np, K3d, 2);

        free(EleMass3d);
        free(InvEleMass3d);
        free(EleMass2d);
        free(LocalPhysicalDiffMatrix);
        free(Dz);
        free(OP11);
        free(LocalRows);
        free(LocalColumns);
        free(SurfBoundStiffTerm);
        free(BotBoundStiffTerm);
        free(StiffMatrix);
    } 
}