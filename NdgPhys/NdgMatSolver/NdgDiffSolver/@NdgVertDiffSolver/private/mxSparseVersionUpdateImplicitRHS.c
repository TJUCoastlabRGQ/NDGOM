
#include "../../../../../NdgMath/NdgMemory.h"

#include "mxImplicitVerticalEddyViscosity.h"

extern int *Ir, *Jc;

extern int NNZ;

extern double *ImTau, *Imu2d, *Imv2d , *GlobalSystemRHS;

extern char *ImVertDiffInitialized, *ImVertEddyInitialized;

void FetchBoundaryData(double *dest, double *source, const int Np2d, double *Eid)
{
    for (int i = 0; i < Np2d; i++){
        dest[i] = source[(int)(Eid[i]) - 1];
    }
}

void SumInColumn( double *dest, double *source, int row, int column){
    for(int col = 0;col < column; col++){
        for(int r = 0; r < row; r++){
            dest[col] += source[col*row + r]; 
        }
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
void VolumnIntegral(double *dest, double *Dz, double *Mass3d, double *diff, int Np)
{
    double zero = 0.0;
    double Alpha1 = -1.0;
    double Alpha2 = 1.0;
    double *Tempdest = malloc(Np*Np*sizeof(double));
	int Np3d = Np;
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
		Np3d, Np3d, Np3d, Alpha1, Dz, Np3d, Mass3d, Np3d, zero, Tempdest, Np3d);
//    dgemm("t", "n", &Np, &Np, &Np, &Alpha1, Dz, &Np, Mass3d, &Np, &zero, Tempdest, &Np);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		Np3d, Np3d, Np3d, Alpha2, Tempdest, Np3d, diff, Np3d, zero, dest, Np3d);
//    dgemm("n", "n", &Np, &Np, &Np, &Alpha2, Tempdest, &Np, diff, &Np, &zero, dest, &Np);
    free(Tempdest);
}

/*This function is used to impose the Newmann boundary at the surface and bottom boundary*/
void ImposeNewmannBoundary(double *Eid, double *mass2d, double* InvMassMatrix3d, double dt, double Imparam, \
        const double* BoundNewmannData, int Np2d, int K2d, double *SystemRHS, int Np3d, int K3d, double *StiffMatrix, int Nvar)
{
    double *tempRHS = malloc(Np3d*sizeof(double));
    double *RHS2d = malloc(Np2d*sizeof(double));
    for (int i = 0; i < Nvar; i++){
        memset(tempRHS, 0, Np3d*sizeof(double));
        memset(RHS2d, 0, Np2d*sizeof(double));
		int colB = 1;
        double Alpha = 1.0;
        double Beta = 0.0;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,\
			Np2d, colB, Np2d, Alpha, mass2d, Np2d, BoundNewmannData + i*Np2d*K2d, Np3d, Beta, RHS2d, Np2d);
        //dgemm("n", "n", &Np2d, &colB, &Np2d, &Alpha, mass2d, &Np2d, BoundNewmannData + i*Np2d*K2d, &Np2d, &Beta, RHS2d, &Np2d);
        for (int j = 0; j < Np2d; j++)
            tempRHS[(int)Eid[j] - 1] = RHS2d[j];
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
			Np3d, colB, Np3d, Alpha, InvMassMatrix3d, Np3d, tempRHS, Np3d, Beta, StiffMatrix + i*Np3d, Np3d);
        //dgemm("n", "n", &Np3d, &colB, &Np3d, &Alpha, InvMassMatrix3d, &Np3d, tempRHS, &Np3d, &Beta, StiffMatrix + i*Np3d, &Np3d);
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

void AssembleLocalToGlobalContribution(double *dest, double *Finaldest, double *invMass, double dt, \
	double ImplicitParam, double *OP11, double Coe, int StartPoint, int NonzeroNum, int Np) {
	double *TempDest = malloc(Np * Np * sizeof(double));
	double Beta = 0.0;
	double alpha = 1.0;
	// Multiply by the inverse mass matrix first
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np, Np, Np, alpha, invMass, Np, OP11, Np, Beta, TempDest, Np);
	//dgemm("n", "n", &Np, &Np, &Np, &alpha, invMass, &Np, OP11, &Np, &Beta, TempDest, &Np);
	// Form the final stiff matrix first
	AssembleContributionIntoSparseMatrix(Finaldest + StartPoint, TempDest, NonzeroNum, Np);
	// Multiply the stiff matrix by parameter -1*dt*Coe*ImplicitParam
	MultiplyByConstant(TempDest, TempDest, -1*dt*Coe*ImplicitParam, Np*Np);
	// Add the diagonal data by 1
	for (int row = 0; row < Np; row++) {
		TempDest[row*Np + row] += 1;
	}
	// Form the stiff matrix used to calculate the final result
	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempDest, NonzeroNum, Np);
	free(TempDest);
}

void AssembleLocalAdjacentToGlobalContribution(double *dest, double *Finaldest, double *AdjInvMass, double dt, \
	double ImplicitParam, double *OP12, double Coe, int StartPoint, int NonzeroNum, int Np) {
	int Np3d = Np;
	double *TempDest = malloc(Np3d * Np3d * sizeof(double));
	double Beta = 0;
	double alpha = 1.0;
	// Multiply by the inverse mass matrix first
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np, Np, Np, alpha, AdjInvMass, Np, OP12, Np, Beta, TempDest, Np);
	//dgemm("n", "n", &Np, &Np, &Np, &alpha, AdjInvMass, &Np, OP12, &Np, &Beta, TempDest, &Np);
	// Form the final stiff matrix first
	AssembleContributionIntoSparseMatrix(Finaldest + StartPoint, TempDest, NonzeroNum, Np3d);
	// Multiply the stiff matrix by parameter -1*dt*ImplicitParam
	MultiplyByConstant(TempDest, TempDest, -1 * dt* Coe *ImplicitParam, Np3d*Np3d);
	// Form the stiff matrix used to calculate the final result
	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempDest, NonzeroNum, Np3d);
	free(TempDest);
}

/*Note: This function is used to calculate the local boundary contribution to the local stiff operator OP11*/
void LocalBoundaryIntegral(double *eid, double *DiffMatrix, double *mass2d, double *TempTau, double *OP11, int Np3d, int Np2d, int Flag, double epsilon)
{
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, DiffMatrix, eid, Np2d, Np3d);
    double Alpha = -1 * epsilon*Flag*0.5, Beta = 0.0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
		Np3d, Np2d, Np2d, Alpha, FDiffMatrix, Np2d, mass2d, Np2d, Beta, EdgeContribution, Np3d);
    //dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
    AssembleContributionIntoColumn(OP11, EdgeContribution, eid, Np3d, Np2d);
    Alpha = 0.5*Flag;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np2d, Np3d, Np2d, Alpha, mass2d, Np2d, FDiffMatrix, Np2d, Beta, EdgeContribution, Np2d);
//    dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
    AssembleContributionIntoRow(OP11, EdgeContribution, eid, Np3d, Np2d);
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
    DiagMultiply(DoubleJump, mass2d, TempTau, Np2d);
    AssembleContributionIntoRowAndColumn(OP11, DoubleJump, eid, eid, Np3d, Np2d, -1);
    free(FDiffMatrix);
    free(EdgeContribution);
    free(DoubleJump);
}
/*Note: This function is used to impose the homogeneous Dirichlet boundary condition at the bottom boundary*/
void ImposeDirichletBoundary(double *eid, double *DiffMatrix, double *mass2d, double *TempTau, double *OP11, int Np3d, int Np2d, int Flag, const double epsilon)
{
    /*	double *Tau = malloc(Np2d*sizeof(double));
     * for (int i = 0; i < Np2d; i++)
     * Tau[i] = 2 * TempTau[i];
     */
    //LocalBoundaryIntegral(eid, DiffMatrix, mass2d, Tau, OP11, Np3d, Np2d, Flag, epsilon);
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, DiffMatrix, eid, Np2d, Np3d);
    double Alpha = -1 * epsilon*Flag*0.5*2, Beta = 0.0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
		Np3d, Np2d, Np2d, Alpha, FDiffMatrix, Np2d, mass2d, Np2d, Beta, EdgeContribution, Np3d);
    //dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
    AssembleContributionIntoColumn(OP11, EdgeContribution, eid, Np3d, Np2d);
    Alpha = 0.5*Flag*2;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np2d, Np3d, Np2d, Alpha, mass2d, Np2d, FDiffMatrix, Np2d, Beta, EdgeContribution, Np2d);
//    dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
    AssembleContributionIntoRow(OP11, EdgeContribution, eid, Np3d, Np2d);
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
    //DiagMultiply(DoubleJump, mass2d, Tau, (int)Np2d);
    DiagMultiply(DoubleJump, mass2d, TempTau, Np2d);
    AssembleContributionIntoRowAndColumn(OP11, DoubleJump, eid, eid, Np3d, Np2d, -1);
    free(FDiffMatrix);
    free(EdgeContribution);
    free(DoubleJump);
//	free(Tau);
}
/*Note: This function is used to calculate the adjacent boundary contribution to the adjacent stiff operator OP12, here adjacent means the test function is defined over the adjacent cell and not the local one*/
void AdjacentBoundaryIntegral(double *eidM, double *eidP, double *LocalDiff, double *AdjacentDiff, double *mass2d, double *TempTau, double *OP12, int Np3d, \
        int Np2d, int Flag, const double epsilon)
{
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, AdjacentDiff, eidP, Np2d, Np3d);
    double Alpha = -1 * epsilon*Flag*0.5, Beta = 0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));

	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
		Np3d, Np2d, Np2d, Alpha, FDiffMatrix, Np2d, mass2d, Np2d, Beta, EdgeContribution, Np3d);
   // dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
    AssembleContributionIntoColumn(OP12, EdgeContribution, eidM, Np3d, Np2d);
    //Alpha = 0.5*Flag;
    Alpha = -1*0.5*Flag;
    AssembleFacialDiffMatrix(FDiffMatrix, LocalDiff, eidM, Np2d, Np3d);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np2d, Np3d, Np2d, Alpha, mass2d, Np2d, FDiffMatrix, Np2d, Beta, EdgeContribution, Np2d);
//    dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
    AssembleContributionIntoRow(OP12, EdgeContribution, eidP, Np3d, Np2d);
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
    DiagMultiply(DoubleJump, mass2d, TempTau, Np2d);
    AssembleContributionIntoRowAndColumn(OP12, DoubleJump, eidP, eidM, Np3d, Np2d, 1);
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
	int RowVCV, int ColVCV, double *hcrit, double *h2d, double *VCV) {
    char *chn = "N";
    double alpha = 1.0;
    double beta = 0.0;
    int Col = 1;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		RowVCV, Col, ColVCV, alpha, VCV, RowVCV, hubot, ColVCV, beta, ucenter, RowVCV);
   // dgemm(chn, chn, RowVCV, &Col, ColVCV, &alpha, VCV, RowVCV, hubot, ColVCV, &beta, ucenter, RowVCV);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		RowVCV, Col, ColVCV, alpha, VCV, RowVCV, hvbot, ColVCV, beta, vcenter, RowVCV);
   // dgemm(chn, chn, RowVCV, &Col, ColVCV, &alpha, VCV, RowVCV, hvbot, ColVCV, &beta, vcenter, RowVCV);
    for (int p = 0; p < RowVCV; p++) {
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

void ImposeImplicitNeumannBoundary(double *dest, double *EidM, double *Cf, double *Depth, \
	double *EleMass2d, double *u2d, double *v2d, int Np2d, int Np3d, double hcrit, double *VCV, double *FacialElemass3d) {
	double *Coe = malloc(Np2d * sizeof(double));
	double *CoeVCV = malloc(Np2d * Np3d * sizeof(double));
	double *TempOP11 = malloc(Np2d * Np3d * sizeof(double));
	double *FinalCoe = malloc(Np3d * sizeof(double));
	memset(FinalCoe, 0, Np3d * sizeof(double));
	for (int p = 0; p < Np2d; p++) {
		if (Depth[p] >= hcrit) {
			//      Coe[p] = -1.0*dt*ImplicitParam*Cf[p] * sqrt(u2d[p] * u2d[p] + v2d[p] * v2d[p]) / Depth[p];
				Coe[p] = -1.0*Cf[p] * sqrt(u2d[p] * u2d[p] + v2d[p] * v2d[p]) / Depth[p];
			//Coe[p] = -1.0*0.005 / Depth[p];
		}
		else {
			Coe[p] = 0.0;
		}
	}

	DiagLeftMultiplyUnsymmetric(CoeVCV, VCV, Coe, Np2d, Np3d);

	SumInColumn(FinalCoe, CoeVCV, Np2d, Np3d);

	DiagRightMultiplyUnsymmetric(TempOP11, FacialElemass3d, FinalCoe, Np2d, Np3d);

	MultiplyByConstant(TempOP11, TempOP11, -1.0, Np3d*Np2d);

	//AssembleContributionIntoRowAndColumn(dest, TempOP11, EidM, EidM, Np3d, Np2d, -1);
	AssembleContributionIntoRow(dest, TempOP11, EidM, Np3d, Np2d);

	free(Coe);

	free(CoeVCV);
	free(TempOP11);
	free(FinalCoe);

}

//This part is not needed, since the implicit part should be returned as OP11
/*
void ImposeImplicitNeumannBoundary(double *dest, double *Finaldest, double *EidM, double *Cf, double dt, \
	double ImplicitParam, double *Depth, double *EleMass2d, double *u2d, double *v2d, int Np2d, int Np3d, \
	double hcrit, double *VCV, double *FacialElemass3d, double *InvLocalElemass, \
	int StartPoint, int NonzeroNum) {

	ptrdiff_t Np = (ptrdiff_t)Np3d;

	double alpha = 1.0, Beta = 0.0;

    double *Coe = malloc(Np2d * sizeof(double)); 

    double *CoeVCV = malloc(Np2d * Np3d *sizeof(double));

    double *TempOP11 = malloc(Np2d * Np3d *sizeof(double));

	double *OP11 = malloc(Np3d * Np3d * sizeof(double));

	memset(OP11, 0, Np3d * Np3d * sizeof(double));

	double *InvOP11 = malloc(Np3d * Np3d * sizeof(double));

	memset(InvOP11, 0, Np3d * Np3d * sizeof(double));

    double *FinalCoe = malloc(Np3d*sizeof(double));

    memset(FinalCoe, 0, Np3d*sizeof(double));

    for (int p = 0; p < Np2d; p++) {
        if (Depth[p] >= hcrit) {
		//	Coe[p] = -1.0*0.005 / Depth[p];
			Coe[p] = -1.0*Cf[p] * sqrt(u2d[p] * u2d[p] + v2d[p] * v2d[p]) / Depth[p];
        }
        else {
            Coe[p] = 0.0;
        }
    }
    
    DiagLeftMultiplyUnsymmetric(CoeVCV, VCV, Coe, Np2d, Np3d);
    
    SumInColumn(FinalCoe, CoeVCV, Np2d, Np3d);
    
    DiagRightMultiplyUnsymmetric(TempOP11, FacialElemass3d, FinalCoe, Np2d, Np3d);
    
    MultiplyByConstant(TempOP11, TempOP11, -1.0, Np3d*Np2d);
      
    //AssembleContributionIntoRowAndColumn(dest, TempOP11, EidM, EidM, Np3d, Np2d, -1);
    AssembleContributionIntoRow(OP11, TempOP11, EidM, Np3d, Np2d);
	// Multiply by the inverse mass matrix first
	dgemm("n", "n", &Np, &Np, &Np, &alpha, InvLocalElemass, &Np, OP11, &Np, &Beta, InvOP11, &Np);
	// Form the final stiff matrix first
	AssembleContributionIntoSparseMatrix(Finaldest + StartPoint, InvOP11, NonzeroNum, Np3d);
	MultiplyByConstant(InvOP11, InvOP11, -1 * dt *ImplicitParam, Np3d*Np3d);
	// Form the stiff matrix used to calculate the final result
	AssembleContributionIntoSparseMatrix(dest + StartPoint, InvOP11, NonzeroNum, Np3d);
    
    free(Coe);
    
    free(CoeVCV);
    free(TempOP11);
	free(OP11);
	free(InvOP11);
    free(FinalCoe);   
}
*/

void GetImplicitBoundaryContribution(double *dest, double *Cf, double *u2d, double *v2d, double *Height, \
        double *hu3d, double *hv3d, double *VCV, int RowVCV, \
        int ColVCV, int Np, int Np2d, double *Mass2d, double *InvMass3d, double *EidM, double hcrit) {
    
    double *TempRHS = malloc(Np * 2 * sizeof(double));
    memset(TempRHS, 0, Np * 2 * sizeof(double));
    double *Drag2d = malloc(Np2d * 2 * sizeof(double));
    double *TempRHS2d = malloc(Np2d * 2 * sizeof(double));
    double *ucenter = malloc(Np2d * sizeof(double));
    double *vcenter = malloc(Np2d * sizeof(double));
    
    CalculateVelocityAtBottomCellCenter(ucenter, vcenter, hu3d, hv3d, \
            RowVCV, ColVCV, &hcrit, Height, VCV);
    
    for (int p = 0; p < Np2d; p++) {
        Drag2d[p] = Cf[p] * sqrt(pow(u2d[p], 2.0) + pow(v2d[p], 2.0))*ucenter[p];
        Drag2d[Np2d + p] = Cf[p] * sqrt(pow(u2d[p], 2.0) + pow(v2d[p], 2.0))*vcenter[p];
    }
    
    char *chn = "N";
    double alpha = 1.0;
    double beta = 0.0;
    int Col = 1;
	int RowMass2d = Np2d;
	int ColMass2d = Np2d;
    
	int RowInvMass3d = Np;
	int ColInvMass3d = Np;

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		RowMass2d, Col, ColMass2d, alpha, Mass2d, RowMass2d, Drag2d, ColMass2d, beta, TempRHS2d, RowMass2d);
    //dgemm(chn, chn, &RowMass2d, &Col, &RowMass2d, &alpha, Mass2d, \
            &RowMass2d, Drag2d, &RowMass2d, &beta, TempRHS2d, &RowMass2d);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		RowMass2d, Col, ColMass2d, alpha, Mass2d, RowMass2d, Drag2d, ColMass2d, beta, TempRHS2d + Np2d, RowMass2d);
    //dgemm(chn, chn, &RowMass2d, &Col, &RowMass2d, &alpha, Mass2d, \
            &RowMass2d, Drag2d + Np2d, &RowMass2d, &beta, TempRHS2d + Np2d, &RowMass2d);
    
    AssembleDataIntoPoint(TempRHS, TempRHS2d, EidM, Np2d);
    
    AssembleDataIntoPoint(TempRHS + Np, TempRHS2d + Np2d, EidM, Np2d);
    
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		RowInvMass3d, Col, ColInvMass3d, alpha, InvMass3d, RowInvMass3d, TempRHS, ColInvMass3d, beta, dest, RowInvMass3d);
    //dgemm(chn, chn, &RowInvMass3d, &Col, &RowInvMass3d, &alpha, InvMass3d, \
            &RowInvMass3d, TempRHS, &RowInvMass3d, &beta, dest, &RowInvMass3d);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		RowInvMass3d, Col, ColInvMass3d, alpha, InvMass3d, RowInvMass3d, TempRHS + Np, ColInvMass3d, beta, dest, RowInvMass3d);
    //dgemm(chn, chn, &RowInvMass3d, &Col, &RowInvMass3d, &alpha, InvMass3d, \
            &RowInvMass3d, TempRHS + Np, &RowInvMass3d, &beta, dest + Np, &RowInvMass3d);
    
    free(TempRHS);
    free(Drag2d);
    free(TempRHS2d);
    free(ucenter);
    free(vcenter);
}

void MyExit()
{
    if ( (!strcmp("True", ImVertDiffInitialized)) && (!strcmp("True", ImVertEddyInitialized)) ){
        ImVertDiffMemoryDeAllocation();
		ImEddyVisInVertDeAllocation();
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
    
    const size_t NdimOut = 3;
    const mwSize dimOut[3] = {Np,K3d,Nvar};
    plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
    
    double *fphys = mxGetPr(plhs[0]);
//    memcpy(fphys, RHS, Np*K3d*Nvar*sizeof(double));
    double *ImplicitRHS = mxGetPr(plhs[1]);
    /*Here, epsilon is set to be -1, and this correspoinding to the SIPG. This value can be set to 1(IIPG or NIPG) or 0(IIPG or NIPG)*/
    const double epsilon = -1;
    
    /*If not initialized, initialize first*/
    if ( (!strcmp("False", ImVertDiffInitialized)) && (!strcmp("False", ImVertEddyInitialized)) )
    {
        ImVertDiffMemoryAllocation(Np2d, K2d, Nz, Np, Nvar);
		ImEddyVisInVertAllocation(Np, Nz);
    }

	memcpy(GlobalSystemRHS, RHS, Np*K3d*Nvar * sizeof(double));
    
    int RowVCV = Np2d;
    int ColVCV = Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int i = 0; i < K2d; i++) {
        CalculateVelocityAtBottomCellCenter(Imu2d + i*Np2d, Imv2d + i*Np2d, hu3d + (i + 1)*Nz*Np - Np, hv3d + (i+1)*Nz*Np - Np, \
                RowVCV, ColVCV, &hcrit, h2d + i*Np2d, VCV);
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif    
    for (int i = 0; i < K2d; i++) {
        CalculatePenaltyParameter(ImTau + i*Np2d*(Nz + 1), Np2d, Np, UpEidM, BotEidM, Diff + i*Np*Nz, Nz, P, Nface);
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int i = 0; i < K2d; i++){
        double *StiffMatrix = malloc(NNZ*Nvar*sizeof(double));
        memset(StiffMatrix, 0, NNZ*Nvar*sizeof(double));
		double *FinalStiffMatrix = malloc(NNZ*Nvar * sizeof(double));
		memset(FinalStiffMatrix, 0, NNZ*Nvar * sizeof(double));
        double *EleMass3d = malloc(Np*Np*sizeof(double));
        double *InvEleMass3d = malloc(Np*Np*sizeof(double));
        double *EleMass2d = malloc(Np2d*Np2d*sizeof(double));
        double *FacialElemass3d = malloc(Np2d*Np*sizeof(double));
        double *LocalPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
        double *Dz = malloc(Np*Np*sizeof(double));
        double *OP11 = malloc(Np*Np*sizeof(double));
        double *LocalRows = malloc(Np*sizeof(double));
        double *LocalColumns = malloc(Np*sizeof(double));
        DiagMultiply(EleMass3d, M3d, J3d + i*Nz*Np, Np);
        /*Fetch the mass matrix corresponding to the facial point from the three dimensional mass matrix*/
        AssembleFacialDiffMatrix(FacialElemass3d, EleMass3d, BotEidM, Np2d, Np);
        memcpy(InvEleMass3d, EleMass3d, Np*Np*sizeof(double));

		ImMatrixInverse(InvEleMass3d, (lapack_int)Np);

        DiagMultiply(EleMass2d, M2d, J2d + i*Np2d, Np2d);
        DiagMultiply(Dz, Dt, tz + i*Nz*Np, Np);
        DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np, Np);
        GetLocalRowsAndColumns(LocalRows, LocalColumns, 0, Np);
        /*Calculate the volumn integral for the first cell, and impose the surface Newmann boundary condition*/
        VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);  //Volume integral part, we need to consider the diagonal contribution
        double *SurfBoundStiffTerm = malloc(Np*Nvar*sizeof(double));
        memset(SurfBoundStiffTerm, 0, Np*Nvar * sizeof(double));
        double *BotBoundStiffTerm = malloc(Np*Nvar*sizeof(double));
        memset(BotBoundStiffTerm, 0, Np*Nvar * sizeof(double));
		/*Treat the surface Neumann boundary condition explicitly*/
        ImposeNewmannBoundary(UpEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, surf + i*Np2d, Np2d, K2d, GlobalSystemRHS + i*Np*Nz, Np, K3d, SurfBoundStiffTerm, Nvar);
		int LocalStartPoint;
		/*The local start point, for pardiso, jc starts from one, so we have to delete this*/
		LocalStartPoint = Jc[0*Np] - 1;
		if (Nz != 1){
            /*When vertical layers greater than one, we calculate the bottom boundary integral part */
            double *BottomAdjacentRows = malloc(Np*sizeof(double));
            double *UpAdjacentRows = malloc(Np*sizeof(double));
            GetBottomRows(BottomAdjacentRows, LocalRows, Np);
            double *BottomPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
            double *UpPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
            DiagMultiply(BottomPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + Np, Np);
            LocalBoundaryIntegral(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + 2 - 1), OP11, Np, Np2d, -1, epsilon);
            double *OP12 = malloc(Np*Np*sizeof(double));
            memset(OP12, 0, Np*Np*sizeof(double));
            AdjacentBoundaryIntegral(BotEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + 2 - 1), OP12, Np, Np2d, -1, epsilon);
            // Local and Local to bottom
            for (int var = 0; var < 2; var++){
				AssembleLocalToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP11, 1.0, LocalStartPoint, Jc[0*Np + 1] - Jc[0 * Np + 0], Np);
				AssembleLocalAdjacentToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP12, 1.0, LocalStartPoint + Np, Jc[0 * Np + 1] - Jc[0 * Np + 0], Np);
            }
			// Local and Local to bottom
            for (int var = 2; var < Nvar; var++){
				AssembleLocalToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP11, Prantl, LocalStartPoint, Jc[0 * Np + 1] - Jc[0 * Np + 0], Np);
				AssembleLocalAdjacentToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP12, Prantl, LocalStartPoint + Np, Jc[0 * Np + 1] - Jc[0 * Np + 0], Np);
			}
            for (int j = 1; j < Nz-1 ; j++){
				/*The local start point, for pardiso, jc starts from one, so we have to delete this*/
				LocalStartPoint = Jc[j*Np] - 1 + Np;
                /*Calculate both the volume and surface integral for the second to the last but one cell*/
                GetLocalRowsAndColumns(LocalRows, LocalColumns, j, Np);
                GetBottomRows(BottomAdjacentRows, LocalRows, Np);
                GetUpRows(UpAdjacentRows, LocalRows, Np);
                DiagMultiply(UpPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (j-1)*Np, Np);
                DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + j*Np, Np);
                DiagMultiply(BottomPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (j + 1)*Np, Np);
                VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
                LocalBoundaryIntegral(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + j+1), OP11, Np, Np2d, -1, epsilon);
                LocalBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + j), OP11, Np, Np2d, 1, epsilon);
                memset(OP12, 0, Np*Np*sizeof(double));
                AdjacentBoundaryIntegral(UpEidM, BotEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + j), OP12, Np, Np2d, 1, epsilon);
				// Local and Local to up
                for (int var = 0; var < 2; var++){
					AssembleLocalToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
						dt, ImplicitParam, OP11, 1.0, LocalStartPoint, Jc[j * Np + 1] - Jc[j * Np + 0], Np);
					AssembleLocalAdjacentToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
						dt, ImplicitParam, OP12, 1.0, LocalStartPoint - Np, Jc[j * Np + 1] - Jc[j * Np + 0], Np);
                }
				// Local and Local to up
                for (int var = 2; var < Nvar; var++){
					AssembleLocalToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
						dt, ImplicitParam, OP11, Prantl, LocalStartPoint, Jc[j * Np + 1] - Jc[j * Np + 0], Np);
					AssembleLocalAdjacentToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
						dt, ImplicitParam, OP12, Prantl, LocalStartPoint - Np, Jc[j * Np + 1] - Jc[j * Np + 0], Np);
                }
                memset(OP12, 0, Np*Np*sizeof(double));
                AdjacentBoundaryIntegral(BotEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + j + 1), OP12, Np, Np2d, -1, epsilon);
				// Local to bottom
				for (int var = 0; var < 2 && var < Nvar; var++){
					AssembleLocalAdjacentToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
						dt, ImplicitParam, OP12, 1.0, LocalStartPoint + Np, Jc[j * Np + 1] - Jc[j * Np + 0], Np);
                }
                for (int var = 2; var < Nvar; var++){
					AssembleLocalAdjacentToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
						dt, ImplicitParam, OP12, Prantl, LocalStartPoint + Np, Jc[j * Np + 1] - Jc[j * Np + 0], Np);
                }
            }
			/*The local start point, for pardiso, jc starts from one, so we have to delete this*/
			LocalStartPoint = Jc[(Nz - 1)*Np] - 1 + Np;
            /*For the bottom most cell, calculate the volumn integral and impose Neumann boundary condition or the Dirichlet boundary condition*/
            GetLocalRowsAndColumns(LocalRows, LocalColumns, Nz - 1, Np);
            GetUpRows(UpAdjacentRows, LocalRows, Np);
            DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (Nz-1)*Np, Np);
            DiagMultiply(UpPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (Nz - 1 - 1)*Np, Np);
            VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
            LocalBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + Nz-1), OP11, Np, Np2d, 1, epsilon);
            /*For passive transport substances, we only impose Neumann boundary condition, so this has no effect on the stiff matrix*/
            // Local for passive transport substances
			for (int var = 2; var < Nvar; var++){
				AssembleLocalToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP11, Prantl, LocalStartPoint, Jc[(Nz-1) * Np + 1] - Jc[(Nz - 1) * Np + 0], Np);
            }
            /*We note that, the Neumann data is zero by default, so the impositon of Neumann BC for hu and hv will not affect the Dirichlet boundary condition for hu and hv*/
            ImposeNewmannBoundary(BotEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, bot + i*Np2d, Np2d, K2d, GlobalSystemRHS + i*Np*Nz + (Nz - 1)*Np, Np, K3d, BotBoundStiffTerm, Nvar);
            /*The following is used to add homogeneous dirichlet boundary for hu and hv*/
            if (!strcmp(BoundaryType, "Dirichlet")) {
                ImposeDirichletBoundary(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + Nz), OP11, Np, Np2d, -1, epsilon);
            }
            else {
                // Impose implicit Neumann boundary condition here, this part is added on 20211231
                ImposeImplicitNeumannBoundary(OP11, BotEidM, Cf + i*Np2d, h2d + i*Np2d, EleMass2d, Imu2d + i*Np2d, Imv2d + i*Np2d, Np2d, Np, hcrit, VCV, FacialElemass3d);
            }
            
            memset(OP12, 0, Np*Np*sizeof(double));
            AdjacentBoundaryIntegral(UpEidM, BotEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + Nz - 1), OP12, Np, Np2d, 1.0, epsilon);
			//Local and Local to up for hu and hv
			for (int var = 0; var < 2; var++) {
				AssembleLocalToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP11, 1.0, LocalStartPoint, Jc[(Nz - 1) * Np + 1] - Jc[(Nz - 1) * Np + 0], Np);
				AssembleLocalAdjacentToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP12, 1.0, LocalStartPoint - Np, Jc[(Nz - 1) * Np + 1] - Jc[(Nz - 1) * Np + 0], Np);
			}
			// Local to up for passive transport substances
			for (int var = 2; var < Nvar; var++) {
				AssembleLocalAdjacentToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP12, Prantl, LocalStartPoint - Np, Jc[(Nz - 1) * Np + 1] - Jc[(Nz - 1) * Np + 0], Np);
			}
            free(BottomAdjacentRows);
            free(UpAdjacentRows);
            free(BottomPhysicalDiffMatrix);
            free(UpPhysicalDiffMatrix);
            free(OP12);
        }
        else{/*If only one layer included in vertical direction, we need to consider the bottom face, and the Upper face has been considered in the first part, from 291-312*/
            
			 /*For passive transport substances, we only impose Neumann boundary condition*/
			 // Local for passive transport substances
			for (int var = 2; var < Nvar; var++) {
				AssembleLocalToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP11, Prantl, LocalStartPoint, Jc[0 * Np + 1] - Jc[0 * Np + 0], Np);
			}
            /*We note that, the Neumann data is zero by default, so the impositon of Neumann BC for hu and hv will not affect the Dirichlet boundary condition for hu and hv*/
            ImposeNewmannBoundary(BotEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, bot + i*Np2d, Np2d, K2d, GlobalSystemRHS + i*Np*Nz + (Nz - 1)*Np, Np, K3d, BotBoundStiffTerm, Nvar);
            /*The following is used to add homogeneous dirichlet boundary for hu and hv*/
            if (!strcmp(BoundaryType, "Dirichlet")) {
                ImposeDirichletBoundary(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + Nz), OP11, Np, Np2d, -1, epsilon);
            }
            else {
                // Impose implicit Neumann boundary condition here, this part is added on 20211231.
                ImposeImplicitNeumannBoundary(OP11, BotEidM, Cf + i*Np2d, h2d + i*Np2d, EleMass2d, Imu2d + i*Np2d, Imv2d + i*Np2d, Np2d, Np, hcrit, VCV, FacialElemass3d);
            }
			//Local for hu and hv
			for (int var = 0; var < 2; var++) {
				AssembleLocalToGlobalContribution(StiffMatrix + var*NNZ, FinalStiffMatrix + var*NNZ, InvEleMass3d, \
					dt, ImplicitParam, OP11, 1.0, LocalStartPoint, Jc[0 * Np + 1] - Jc[0 * Np + 0], Np);
			}
        }

		/*Invoke pardiso from mkl to solve equation Ax = b*/
		SparseEquationSolve(fphys, Nz*Np, StiffMatrix, GlobalSystemRHS, Nz, Np, i, Nvar, K3d);

		for (int var = 0; var < Nvar; var++) {
			SparseMatrixMultiply(ImplicitRHS + var*K3d*Np + i*Nz*Np, FinalStiffMatrix + var*NNZ, fphys + var * Np*K3d + i*Nz*Np, Nz * Np, Jc, Ir);
		}

        AssembleBoundaryContribution(ImplicitRHS  + i*Nz*Np, SurfBoundStiffTerm, Np, K3d, Nvar);
        // BotBoundStiffTerm with index 1 and 2 stands for hu and hv respectively, they are zero since they are treated implicitly.
        // We add them to the ImplicitRHS since they have no effect on the final result.
        AssembleBoundaryContribution(ImplicitRHS  + i*Nz*Np + (Nz - 1)*Np, BotBoundStiffTerm, Np, K3d, Nvar);
        // We first need to calculate the contribution to the right hand side due to the implicit bottom condition, then add they to the RHS
        GetImplicitBoundaryContribution(BotBoundStiffTerm, Cf + i*Np2d, Imu2d + i*Np2d, Imv2d + i*Np2d, h2d + i*Np2d, \
                fphys + (i + 1)*Np*Nz - Np, fphys + Np*K3d + (i + 1)*Np*Nz - Np, VCV, RowVCV, ColVCV, Np, Np2d, EleMass2d, InvEleMass3d, \
                BotEidM, hcrit);
        
        //We add the right hand side due to the implicit bottom friction term back, only hu and hv considered. Added on 20211231 by RGQ
        AssembleBoundaryContribution(ImplicitRHS + i*Nz*Np + (Nz - 1)*Np, BotBoundStiffTerm, Np, K3d, 2);

        free(EleMass3d);
        free(InvEleMass3d);
        free(FacialElemass3d);
        free(EleMass2d);
        free(LocalPhysicalDiffMatrix);
        free(Dz);
        free(OP11);
        free(LocalRows);
        free(LocalColumns);
        free(SurfBoundStiffTerm);
        free(BotBoundStiffTerm);
        free(StiffMatrix);
		free(FinalStiffMatrix);
    }
}