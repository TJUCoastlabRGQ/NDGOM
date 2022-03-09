
#include "../../../../../NdgMath/NdgMemory.h"

#include "mxImplicitVerticalEddyViscosity.h"

/*We note that, the discretization presented in this function is for the following second order operator
* $-\nabla\cdot\left(\bold K\nabla q\right) = f$, and the primal form for this operator is as follows,
$\int_{\Omega}\nabla_h v\cdot\left(\bold K\nabla_h q_h\right)d\boldsymbol x+\\
\int_{\epsilon_I}\tau[q_h]\cdot [v]d\boldsymbol x+ \int_{\partial \Omega^D}\tau u_h vd\boldsymbol x-\\
\int_{\epsilon_I}\{\bold K\nabla_h q_h\}\cdot [v]d\boldsymbol x-\\
\int_{\partial \Omega^D}v\bold K\nabla_hq_h\cdot\boldsymbol nd\boldsymbol x\\
+\epsilon\int_{\epsilon_I}[q_h]\cdot\{\bold K\nabla_h v\}d\boldsymbol x +\\
\epsilon \int_{\partial \Omega^D}q_h\bold K\nabla_hv\cdot\boldsymbol nd\boldsymbol x = \\
\int_{\Omega}fvd\boldsymbol x+\int_{\partial \Omega^D}\tau q_Dvd\boldsymbol x-\\
\int_{\partial \Omega^D}q_D\bold K\nabla_h v\cdot\boldsymbol nd\boldsymbol x + \\
\int_{\partial \Omega^N}vq_Nd\boldsymbol x$
*/

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


/*Note: This function is used to calculate the volumn integral contained in primal form*/
/*For term $\int_{\Omega}\nabla_h v\cdot (\bold K\nabla_h q_h)d\boldsymbol x$, added on 20220217 by RGQ*/
void VolumnIntegral(double *dest, double *Dz, double *Mass3d, double *diff, int Np)
{
    double zero = 0.0;
    double Alpha1 = 1.0;
    double Alpha2 = 1.0;
    double *Tempdest = malloc(Np*Np*sizeof(double));
	int Np3d = Np;
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
		Np3d, Np3d, Np3d, Alpha1, Dz, Np3d, Mass3d, Np3d, zero, Tempdest, Np3d);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		Np3d, Np3d, Np3d, Alpha2, Tempdest, Np3d, diff, Np3d, zero, dest, Np3d);
    free(Tempdest);
}

/*This function is used to impose the Newmann boundary explicitly at the surface and bottom boundary*/
/*For term $dt\times Imparam\times\int_{\partial \Omega^N}vq_Nd\boldsymbol x$, added on 20220217 by RGQ*/
void ImposeNewmannBoundary(double *Eid, double *mass2d, double* InvMassMatrix3d, double dt, double Imparam, \
        const double* BoundNewmannData, int Np2d, int K2d, double *SystemRHS, int Np3d, int K3d, double *StiffMatrix,\
	int Nvar)
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

        for (int j = 0; j < Np2d; j++)
            tempRHS[(int)Eid[j] - 1] = RHS2d[j];
		/*RHS is needed to assemble the final rhs corresponding to implicit part*/
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
			Np3d, colB, Np3d, Alpha, InvMassMatrix3d, Np3d, tempRHS, Np3d, Beta, StiffMatrix + i*Np3d, Np3d);
        for (int j = 0; j < Np3d; j++)
            SystemRHS[i*K3d*Np3d + j] += dt*Imparam*(*(StiffMatrix + i*Np3d + j));
    }
    free(tempRHS);
    free(RHS2d);
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
	/* Form the final stiff matrix first*/
	AssembleContributionIntoSparseMatrix(Finaldest + StartPoint, TempDest, NonzeroNum, Np);
	/*Multiply the stiff matrix by parameter dt*Coe*ImplicitParam */
	MultiplyByConstant(TempDest, TempDest, dt*Coe*ImplicitParam, Np*Np);
	/* Add the diagonal data by 1 */
	for (int row = 0; row < Np; row++) {
		TempDest[row*Np + row] += 1;
	}
	/* Form the stiff matrix used to calculate the final result */
	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempDest, NonzeroNum, Np);
	free(TempDest);
}

void AssembleLocalAdjacentToGlobalContribution(double *dest, double *Finaldest, double *AdjInvMass, double dt, \
	double ImplicitParam, double *OP12, double Coe, int StartPoint, int NonzeroNum, int Np) {
	double *TempDest = malloc(Np * Np * sizeof(double));
	double Beta = 0;
	double alpha = 1.0;
	/* Multiply by the inverse mass matrix first*/
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np, Np, Np, alpha, AdjInvMass, Np, OP12, Np, Beta, TempDest, Np);
	/* Form the final stiff matrix first*/
	AssembleContributionIntoSparseMatrix(Finaldest + StartPoint, TempDest, NonzeroNum, Np);
	/* Multiply the stiff matrix by parameter dt*ImplicitParam */
	MultiplyByConstant(TempDest, TempDest, dt* Coe *ImplicitParam, Np*Np);
	/* Form the stiff matrix used to calculate the final result */
	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempDest, NonzeroNum, Np);
	free(TempDest);
}


void LocalBoundaryIntegral(double *eid, double *DiffMatrix, double *mass2d, double *TempTau, double *OP11, int Np3d, int Np2d, double nz, double epsilon)
{
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, DiffMatrix, eid, Np2d, Np3d);
    double Alpha = epsilon*nz*0.5, Beta = 0.0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));
	/*For term $+\int_{\epsilon_I}[q_h]\left \{\bold K\nabla_h v\right \}d\boldsymbol x$*/
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
		Np3d, Np2d, Np2d, Alpha, FDiffMatrix, Np2d, mass2d, Np2d, Beta, EdgeContribution, Np3d);
    AssembleContributionIntoColumn(OP11, EdgeContribution, eid, Np3d, Np2d);
    Alpha = -0.5*nz;
	/*For term $-\int_{\epsilon_I}\left{\bold K\nabla_hq_h\right}\cdot [v]$*/
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np2d, Np3d, Np2d, Alpha, mass2d, Np2d, FDiffMatrix, Np2d, Beta, EdgeContribution, Np2d);
    AssembleContributionIntoRow(OP11, EdgeContribution, eid, Np3d, Np2d);
	/*For term $\int_{\epsilon_I}\tau[q_h]\cdot [v]d\boldsymbol x$, right multiply the mass2d with the penalty parameter*/
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
	DiagRightMultiply(DoubleJump, mass2d, TempTau, Np2d);
    AssembleContributionIntoRowAndColumn(OP11, DoubleJump, eid, eid, Np3d, Np2d, 1);
    free(FDiffMatrix);
    free(EdgeContribution);
    free(DoubleJump);
}
/*Note: This function is used to impose the homogeneous Dirichlet boundary condition at the bottom boundary*/
void ImposeDirichletBoundary(double *eid, double *DiffMatrix, double *mass2d, double *TempTau, double *OP11, int Np3d, int Np2d, double nz, const double epsilon)
{
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, DiffMatrix, eid, Np2d, Np3d);
    double Alpha = epsilon*nz*0.5*2, Beta = 0.0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));
	/*For term $\epsilon \int_{\partial \Omega^D}q_h\bold K\nabla_h v\cdot\boldsymbol n d\boldsymbol x$*/
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
		Np3d, Np2d, Np2d, Alpha, FDiffMatrix, Np2d, mass2d, Np2d, Beta, EdgeContribution, Np3d);
    AssembleContributionIntoColumn(OP11, EdgeContribution, eid, Np3d, Np2d);
    Alpha = -1.0*nz;
	/*For term $-\int_{\partial \Omega^D}v\bold K\nabla_hq_h\cdot\boldsymbol nd\boldsymbol x$*/
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np2d, Np3d, Np2d, Alpha, mass2d, Np2d, FDiffMatrix, Np2d, Beta, EdgeContribution, Np2d);
    AssembleContributionIntoRow(OP11, EdgeContribution, eid, Np3d, Np2d);
	/*For term $\int_{\partial \Omega^D}\tau q_h vd\boldsymbol x$*/
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
    DiagRightMultiply(DoubleJump, mass2d, TempTau, Np2d);
    AssembleContributionIntoRowAndColumn(OP11, DoubleJump, eid, eid, Np3d, Np2d, 1);
    free(FDiffMatrix);
    free(EdgeContribution);
    free(DoubleJump);
}
/*Note: This function is used to calculate the adjacent boundary contribution to the adjacent stiff operator OP12, here adjacent means the test function is defined over the adjacent cell and not the local one
* Also note here nz is the outward normal from the local element
*/
void AdjacentBoundaryIntegral(double *eidM, double *eidP, double *LocalDiff, double *AdjacentDiff, double *mass2d, double *TempTau, double *OP12, int Np3d, \
        int Np2d, double nz, const double epsilon)
{
    double *FDiffMatrix = malloc(Np2d*Np3d*sizeof(double));
    AssembleFacialDiffMatrix(FDiffMatrix, AdjacentDiff, eidP, Np2d, Np3d);
    double Alpha = epsilon*nz*0.5, Beta = 0;
    double *EdgeContribution = malloc(Np2d*Np3d*sizeof(double));
	/*For term $\epsilon\int_{\epsilon_I}[q_h]\cdot\{\bold K\nabla_h v\}$*/
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
		Np3d, Np2d, Np2d, Alpha, FDiffMatrix, Np2d, mass2d, Np2d, Beta, EdgeContribution, Np3d);
    AssembleContributionIntoColumn(OP12, EdgeContribution, eidM, Np3d, Np2d);
    Alpha = (-1.0)*0.5*(-1.0*nz);
    AssembleFacialDiffMatrix(FDiffMatrix, LocalDiff, eidM, Np2d, Np3d);
	/*For term $-\int_{\epsilon_I}\left{\bold K\nabla_h q_h\right}\cdot [v]d\boldsymbol x$*/
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np2d, Np3d, Np2d, Alpha, mass2d, Np2d, FDiffMatrix, Np2d, Beta, EdgeContribution, Np2d);
    AssembleContributionIntoRow(OP12, EdgeContribution, eidP, Np3d, Np2d);
    double *DoubleJump = malloc(Np2d*Np2d*sizeof(double));
	/*For term $\int_{\epsilon_I}\tau [q_h]\cdot [v]d\boldsymbol x$*/
    DiagRightMultiply(DoubleJump, mass2d, TempTau, Np2d);
    AssembleContributionIntoRowAndColumn(OP12, DoubleJump, eidP, eidM, Np3d, Np2d, -1);
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
    double alpha = 1.0;
    double beta = 0.0;
    int Col = 1;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		RowVCV, Col, ColVCV, alpha, VCV, RowVCV, hubot, ColVCV, beta, ucenter, RowVCV);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		RowVCV, Col, ColVCV, alpha, VCV, RowVCV, hvbot, ColVCV, beta, vcenter, RowVCV);
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
	double *EleMass2d, double *u2d, double *v2d, int Np2d, int Np3d, double hcrit, double *VCV) {
	double *Coe = malloc(Np2d * sizeof(double));
	double *CoeVCV = malloc(Np2d * Np3d * sizeof(double));
	double alpha = 1.0;
	double beta = 0.0;
	double *TempOP11 = malloc(Np2d * Np3d * sizeof(double));

	for (int p = 0; p < Np2d; p++) {
		if (Depth[p] >= hcrit) {
			// the unit outward normal -1.0 first, then move to the left hand side another -1.0
			Coe[p] = (-1.0)*(-1.0)*Cf[p] * sqrt(u2d[p] * u2d[p] + v2d[p] * v2d[p]) / Depth[p];
			//Coe[p] = (-1.0)*(-1.0)*0.005 / Depth[p];
		}
		else {
			Coe[p] = 0.0;
		}
	}

	DiagLeftMultiplyUnsymmetric(CoeVCV, VCV, Coe, Np2d, Np3d);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, \
		Np2d, Np3d, Np2d, alpha, EleMass2d, Np2d, CoeVCV, Np2d, beta, TempOP11, Np2d);

	AssembleContributionIntoRow(dest, TempOP11, EidM, Np3d, Np2d);
	free(Coe);
	free(CoeVCV);
	free(TempOP11);
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
    const double epsilon = 0.0;
    
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
		//FinalStiffMatrix is needed, when we calculate the influence of the diffusion term when IMEXRK time stepping method is used
		double *FinalStiffMatrix = malloc(NNZ*Nvar * sizeof(double));
		memset(FinalStiffMatrix, 0, NNZ*Nvar * sizeof(double));
        double *EleMass3d = malloc(Np*Np*sizeof(double));
        double *InvEleMass3d = malloc(Np*Np*sizeof(double));
        double *EleMass2d = malloc(Np2d*Np2d*sizeof(double));

        double *LocalPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
        double *Dz = malloc(Np*Np*sizeof(double));
        double *OP11 = malloc(Np*Np*sizeof(double));
        DiagMultiply(EleMass3d, M3d, J3d + i*Nz*Np, Np);
        /*Fetch the mass matrix corresponding to the facial point from the three dimensional mass matrix*/
        memcpy(InvEleMass3d, EleMass3d, Np*Np*sizeof(double));

		ImMatrixInverse(InvEleMass3d, (lapack_int)Np);

        DiagMultiply(EleMass2d, M2d, J2d + i*Np2d, Np2d);
        DiagMultiply(Dz, Dt, tz + i*Nz*Np, Np);
        DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np, Np);
        /*Calculate the volumn integral for the first cell, and impose the surface Newmann boundary condition*/
        VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);  
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
            double *BottomPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
            double *UpPhysicalDiffMatrix = malloc(Np*Np*sizeof(double));
            DiagMultiply(BottomPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + Np, Np);
            LocalBoundaryIntegral(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + 2 - 1), OP11, Np, Np2d, -1.0, epsilon);
            double *OP12 = malloc(Np*Np*sizeof(double));
            memset(OP12, 0, Np*Np*sizeof(double));
            AdjacentBoundaryIntegral(BotEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + 2 - 1), OP12, Np, Np2d, -1.0, epsilon);
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
                DiagMultiply(UpPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (j-1)*Np, Np);
                DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + j*Np, Np);
                DiagMultiply(BottomPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (j + 1)*Np, Np);
                VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
				/*The lower surface*/
                LocalBoundaryIntegral(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + j+1), OP11, Np, Np2d, -1.0, epsilon);
                /*The upper surface*/
				LocalBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + j), OP11, Np, Np2d, 1.0, epsilon);
                memset(OP12, 0, Np*Np*sizeof(double));
				/*Local element to up element*/
                AdjacentBoundaryIntegral(UpEidM, BotEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + j), OP12, Np, Np2d, 1.0, epsilon);
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
				/*Local element to bottom element*/
                AdjacentBoundaryIntegral(BotEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + j + 1), OP12, Np, Np2d, -1.0, epsilon);
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
            DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (Nz-1)*Np, Np);
            DiagMultiply(UpPhysicalDiffMatrix, Dz, Diff + i*Nz*Np + (Nz - 1 - 1)*Np, Np);
            VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
			/*The upper surface*/
            LocalBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + Nz-1), OP11, Np, Np2d, 1.0, epsilon);
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
                ImposeDirichletBoundary(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + Nz), OP11, Np, Np2d, -1.0, epsilon);
            }
            else {
                // Impose implicit Neumann boundary condition here, this part is added on 20211231
                ImposeImplicitNeumannBoundary(OP11, BotEidM, Cf + i*Np2d, h2d + i*Np2d, EleMass2d, Imu2d + i*Np2d, Imv2d + i*Np2d, Np2d, Np, hcrit, VCV);
            }
            
            memset(OP12, 0, Np*Np*sizeof(double));
			/*Local element to up element*/
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
                ImposeDirichletBoundary(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, ImTau + Np2d*(i*(Nz + 1) + Nz), OP11, Np, Np2d, -1.0, epsilon);
            }
            else {
                // Impose implicit Neumann boundary condition here, this part is added on 20211231.
                ImposeImplicitNeumannBoundary(OP11, BotEidM, Cf + i*Np2d, h2d + i*Np2d, EleMass2d, Imu2d + i*Np2d, Imv2d + i*Np2d, Np2d, Np, hcrit, VCV);
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

		free(StiffMatrix);
		free(FinalStiffMatrix);
        free(EleMass3d);
        free(InvEleMass3d);
        free(EleMass2d);
        free(LocalPhysicalDiffMatrix);
        free(Dz);
        free(OP11);
        free(SurfBoundStiffTerm);
        free(BotBoundStiffTerm);
    }
}