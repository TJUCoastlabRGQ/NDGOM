#include "SWENonhydrostatic3d.h"

/*This function is called at the initialization stage to calcualte the second order derivative
about nonhydrostatic pressure in vertical direction, i.e. $\frac{\partial^2 p}{\partial \sigma^2}$.
For this term, the primal form is given as:
$$\begin {split}-\int_{\Omega}\nabla_h s\cdot \nabla_h p_hd\boldsymbol{x}+\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}-
\int_{\partial \Omega^D}(p_D-p_h)\nabla_hs\cdot \boldsymbol{n}d\boldsymbol{x} =\\-\int_{\epsilon_i}(\left\{\nabla_hp_h\right\}-\tau^k[p_h])
[s]d\boldsymbol{x}-\int_{\partial\Omega^D}s\boldsymbol{n}\cdot\nabla_hp_hd\boldsymbol{x}+\int_{\partial \Omega^D}
\tau^ksp_hd\boldsymbol{x} -\int_{\partial \Omega^D}\tau^ksp_Dd\boldsymbol{x} -
\int_{\partial \Omega^N}sp_Nd\boldsymbol{x}\end {split}$$
Detail about this form can be found through any markdown editor.
*/


/*This following funciton is used to assemble term corresponds to
$$-\int_{\Omega}\nabla_h s\cdot \nabla_h p_hd\boldsymbol{x}$$
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
Startpoint: Number of nonzeros before the current studied element
in the sparse matrix
Np: Number of interpolation point for the 3d element
NonzeroPerColumn: Number of nonzeros in each column of the studied scope,
i.e. nonzeros in the column corresponds to global index
of the studied point
mass3d: A pointer to the start of the mass matrix for the master cell
Dt: A pointer to the start of the differential matrix for the master cell
tz: A pointer to the start of the transformation Jacobian of the studied element
J: A pointer to the start of the Jacobian of the studied element
Layer: The index of the layer of the studied cell. This parameter is used here,
for we need it to decide the exact start address in sparse matrix. If
layer = 0, we start from dest + Startpoint, otherwise we start from
dest + Startpoint + Np.
*/
void GetLocalVolumuIntegralTermForSecondOrderTerm(double *dest, int Startpoint, \
	int Np, int NonzeroPerColumn, double *mass3d, \
	double *Dt, double *tz, double *J, int Layer){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvMass3d, (ptrdiff_t)Np);
	double *DiffSigma = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffSigma, Dt, tz, Np);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	double *TempContributionBuff = malloc(Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffSigma,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffSigma, (ptrdiff_t)Np, 0.0, TempContribution, (ptrdiff_t)Np);
	MultiplyByConstant(TempContribution, TempContribution, -1, Np*Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	if (Layer == 0){
		AssembleContributionIntoSparseMatrix(dest + Startpoint, Contribution, NonzeroPerColumn, Np);
	}
	else{
		AssembleContributionIntoSparseMatrix(dest + Startpoint + Np, Contribution, NonzeroPerColumn, Np);
	}
	free(Mass3d);
	free(InvMass3d);
	free(DiffSigma);
	free(TempContribution);
	free(TempContributionBuff);
	free(Contribution);
}

/*The following funciton is used to assemble terms corresponding to
$$-\int_{\partial \Omega^D}(p_D-p_h)\nabla_hs\cdot \boldsymbol{n}d\boldsymbol{x} +
\int_{\partial\Omega^D}s\boldsymbol{n}\cdot\nabla_hp_hd\boldsymbol{x}-\int_{\partial \Omega^D}
\tau^ksp_hd\boldsymbol{x} +\int_{\partial \Omega^D}\tau^ksp_Dd\boldsymbol{x}$$
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForSecondOrderTerm.
The rest are as follows:
Np2d: Number of interpolation point for the 2d element.
For the current case, the cell on top of the 3d element.
mass2d: A pointer to the start of the mass matrix for the 2 dimensional master cell
J2d: A pointer to the start of the Jacobian of the 2 dimensional studied element
Eid: A pointer to the start of the interpolation point index on the top of the master cell
Tau: A pointer to the penalty parameter of the upmost face
*/
void ImposeSecondOrderNonhydroDirichletBoundaryCondition(double *dest, int StartPoint, int Np, \
	int Np2d, int NonzeroPerColumn, double *mass3d, double *mass2d, double *J, \
	double *J2d, double *Eid, double *Tau, double *Dt, double *tz){
	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *DiffSigma = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffSigma, Dt, tz, Np);
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffSigma, Eid, Np2d, Np);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	/*Allocate memory for term $-\int_{\partial \Omega^D}(p_D-p_h)\nabla_hs\cdot \boldsymbol{n}d\boldsymbol{x}$, here P_D = 0*/
	double *DirichletBuff1 = malloc(Np2d * Np * sizeof(double));
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, DirichletBuff1, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, DirichletBuff1, Eid, Np, Np2d);

	/*The follow part is for term $\int_{\partial\Omega^D}s\boldsymbol{n}\cdot\nabla_hp_hd\boldsymbol{x}$*/
	/*("N","N", rows  of OP(A), columns of OP(B), columns of op(A), ALPHA, A, the first dimension of A, B,
	the first dimension of B, BETA, C, the first dimension of C )
	*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, Mass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, DirichletBuff1, (ptrdiff_t)Np2d);
	AssembleContributionIntoRow(TempContribution, DirichletBuff1, Eid, Np, Np2d);

	/*The follow part is for term $-\int_{\partial \Omega^D}\tau^ksp_hd\boldsymbol{x}$*/
	double *DirichletBuff2 = malloc(Np2d * Np2d * sizeof(double));
	DiagMultiply(DirichletBuff2, Mass2d, Tau, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, DirichletBuff2, Eid, Eid, Np, Np2d, -1);

	/*Multiply the contribution by inverse matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(Mass3d);
	free(InvMass3d);
	free(Mass2d);
	free(DiffSigma);
	free(FacialDiffMatrix);
	free(TempContribution);
	free(Contribution);
	free(DirichletBuff1);
	free(DirichletBuff2);
}


/*The following funciton is used to assemble terms corresponding to
$$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}
+\int_{\epsilon_i}(\left\{\nabla_hp_h\right\}-\tau^k[p_h])[s]d\boldsymbol{x}$$
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForSecondOrderTerm
and ImposeSecondOrderNonhydroDirichletBoundaryCondition.
The rest are as follows:
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
UpEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
Uptz: A pointer to the start of the transformation Jacobian of the element located upside of the studied cell
Localtz: A pointer to the start of the transformation Jacobian of the element of the studied cell
Tau: A pointer to the penalty parameter of the up face
*/

void GetLocalToUpElementContributionForSecondOrderTerm(double *dest, int StartPoint, int Np, int Np2d, int NonzeroPerColumn, \
	double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *UpEid, double *Dt, double *Uptz,\
	double *Localtz, double *Tau){
	double *UpMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(UpMass3d, mass3d, J, Np);
	double *InvUpMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvUpMass3d, UpMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvUpMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	double *UpDiffSigma = malloc(Np*Np*sizeof(double));
	DiagMultiply(UpDiffSigma, Dt, Uptz, Np);
	double *LocalDiffSigma = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiffSigma, Dt, Localtz, Np);
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));


	/*The follow part is for term $\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$*/
	double *InnerEdgeContribution = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, UpDiffSigma, UpEid, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5, Np*Np2d);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Np2d);


	/*The follow part is for term $\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiffSigma, LocalEid, Np2d, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, Mass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np2d);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, UpEid, Np, Np2d);


	/*The follow part is for term $ - \int_{ \epsilon_i }\tau^k[p_h][s]d\boldsymbol{ x }$*/
	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(TempMass2d, mass2d, Tau, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, UpEid, LocalEid, Np, Np2d, 1);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvUpMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);
	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(UpMass3d);
	free(InvUpMass3d);
	free(Mass2d);
	free(TempContribution);
	free(Contribution);
	free(UpDiffSigma);
	free(LocalDiffSigma);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(TempMass2d);
}


/*The following funciton is used to assemble terms corresponding to
$$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}
+\int_{\epsilon_i}(\left\{\nabla_hp_h\right\}-\tau^k[p_h])[s]d\boldsymbol{x}$$
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForSecondOrderTerm
and ImposeSecondOrderNonhydroDirichletBoundaryCondition.
The rest are as follows:
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
BotEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
Downtz: A pointer to the start of the transformation Jacobian of the element located downside of the studied cell
Localtz: A pointer to the start of the transformation Jacobian of the element of the studied cell
Tau: A pointer to the penalty parameter of the down face
*/

void GetLocalToDownElementContributionForSecondOrderTerm(double *dest, int StartPoint, int Np, int Np2d, int NonzeroPerColumn, \
	double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *BotEid, double *Dt, double *Downtz, \
	double *Localtz, double *Tau){
	double *DownMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(DownMass3d, mass3d, J, Np);
	double *InvDownMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvDownMass3d, DownMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvDownMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	double *DownDiffSigma = malloc(Np*Np*sizeof(double));
	DiagMultiply(DownDiffSigma, Dt, Downtz, Np);
	double *LocalDiffSigma = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiffSigma, Dt, Localtz, Np);
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	/*The follow part is for term $\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$*/
	double *InnerEdgeContribution = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DownDiffSigma, BotEid, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5, Np*Np2d);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Np2d);

	/*The follow part is for term $\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiffSigma, LocalEid, Np2d, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, Mass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np2d);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, BotEid, Np, Np2d);

	/*The follow part is for term $ - \int_{ \epsilon_i }\tau^k[p_h][s]d\boldsymbol{ x }$*/
	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(TempMass2d, mass2d, Tau, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, BotEid, LocalEid, Np, Np2d, 1);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvDownMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(DownMass3d);
	free(InvDownMass3d);
	free(Mass2d);
	free(TempContribution);
	free(Contribution);
	free(DownDiffSigma);
	free(LocalDiffSigma);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(TempMass2d);
}

/*The following funciton is used to assemble terms corresponding to
$$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}
+\int_{\epsilon_i}(\left\{\nabla_hp_h\right\}-\tau^k[p_h])[s]d\boldsymbol{x}$$
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForSecondOrderTerm
and ImposeSecondOrderNonhydroDirichletBoundaryCondition.
The rest are as follows:
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
Localtz: A pointer to the start of the transformation Jacobian of the element of the studied cell
Tau: A pointer to the penalty parameter of the common face shared by local element and downside element
*/

void GetLocalDownContributionForSecondOrderTerm(double *dest, int StartPoint, int Np, int Np2d, int NonzeroPerColumn, \
	double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *Dt, double *Localtz,\
	double *Tau){
	double *LocalMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalMass3d, mass3d, J, Np);
	double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	double *LocalDiffSigma = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiffSigma, Dt, Localtz, Np);
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	/*The follow part is for term $\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$*/
	double *InnerEdgeContribution = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiffSigma, LocalEid, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5, Np*Np2d);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Np2d);

	/*The follow part is for term $\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiffSigma, LocalEid, Np2d, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, Mass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np2d);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEid, Np, Np2d);

	/*The follow part is for term $ - \int_{ \epsilon_i }\tau^k[p_h][s]d\boldsymbol{ x }$*/
	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(TempMass2d, mass2d, Tau, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(LocalMass3d);
	free(InvLocalMass3d);
	free(Mass2d);
	free(TempContribution);
	free(Contribution);
	free(LocalDiffSigma);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(TempMass2d);
}

/*The following funciton is used to assemble terms corresponding to
$$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}
+\int_{\epsilon_i}(\left\{\nabla_hp_h\right\}-\tau^k[p_h])[s]d\boldsymbol{x}$$
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForSecondOrderTerm
and ImposeSecondOrderNonhydroDirichletBoundaryCondition.
The rest are as follows:
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
Localtz: A pointer to the start of the transformation Jacobian of the element of the studied cell
Tau: A pointer to the penalty parameter of the common face shared by local element and upside element
*/

void GetLocalUpContributionForSecondOrderTerm(double *dest, int StartPoint, int Np, int Np2d, int NonzeroPerColumn, \
	double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *Dt, double *Localtz, \
	double *Tau){
	double *LocalMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalMass3d, mass3d, J, Np);
	double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	double *LocalDiffSigma = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiffSigma, Dt, Localtz, Np);
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	/*The follow part is for term $\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$*/
	double *InnerEdgeContribution = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiffSigma, LocalEid, Np2d, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5, Np*Np2d);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Np2d);

	/*The follow part is for term $\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiffSigma, LocalEid, Np2d, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, Mass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, InnerEdgeContribution, (ptrdiff_t)Np2d);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEid, Np, Np2d);

	/*The follow part is for term $ - \int_{ \epsilon_i }\tau^k[p_h][s]d\boldsymbol{ x }$*/
	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(TempMass2d, mass2d, Tau, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(LocalMass3d);
	free(InvLocalMass3d);
	free(Mass2d);
	free(TempContribution);
	free(Contribution);
	free(LocalDiffSigma);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(TempMass2d);
}


void CalculatePenaltyParameter(double *dest, const int Np2d, int Nz, int P, int Nface)
{
	//Note: the penalty parameter for the bottom most face of each column is not needed, so we leave them undefined 
	for (int Layer = 0; Layer < Nz; Layer++)
	{
		for (int p = 0; p < Np2d; p++)
		{
			dest[Layer*Np2d + p] = ((P + 1)*(P + 3) / 3.0)*(Nface / 2.0)*Nz;
		}
	}
}

/*The input parameters are organized as follows:
Np: Number of interpolation points, indexed as 0.
Ele3d: Number of elements in three-dimensional mesh object, indexed as 1.
Nlayer: Number of layers in vertical direction, indexed as 2.
Ele2d: Number of elements in two-dimensional mesh object, indexed as 3.
Nface: Number of face numbers for the three-dimensional master cell, indexed as 4.
EToE: Element to element topological relation of the 3d mesh, indexed as 5.
tz: Jacobian coefficient in vertical direction, indexed as 6.
Dt: Differential matrix of the 3d master cell in direction t, indexed as 7.
J: Jacobian coefficient of the three-dimensional mesh object, indexed as 8.
J2d: Jacobian coefficient of the 2d mesh, indexed as 9.
M2d: Mass matrix of the 2d master cell of the 2d mesh, indexed as 10.
M3d: Mass matrix of the 3d master cell, indexed as 11.
Fmask: Index of the interpolation point on each face of the master cell, indexed as 12.
Nfp: Number of interpolation point on each face of the master cell, indexed as 13.
P: Approximation order of the master cell, indexed as 14.
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int Np = (int)mxGetScalar(prhs[0]);
	int Ele3d = (int)mxGetScalar(prhs[1]);
	int Nlayer = (int)mxGetScalar(prhs[2]);
	int Ele2d = (int)mxGetScalar(prhs[3]);
	int Nface = (int)mxGetScalar(prhs[4]);
	double *EToE = mxGetPr(prhs[5]);
	int TotalNonzero = Ele3d * 3 * Np*Np - Ele2d * 2 * Np*Np;
	mwIndex *TempIr = malloc(TotalNonzero*sizeof(mwIndex));
	mwIndex *TempJc = malloc((Np*Ele3d + 1)*sizeof(mwIndex));
	TempJc[0] = 0;
	double *tz = mxGetPr(prhs[6]);
	double *Dt = mxGetPr(prhs[7]);
	//double *Mass3d = malloc(Np*Np*sizeof(double));
	double *J = mxGetPr(prhs[8]);
	double *J2d = mxGetPr(prhs[9]);
	double *M2d = mxGetPr(prhs[10]);
	double *M3d = mxGetPr(prhs[11]);
	double *Fmask = mxGetPr(prhs[12]);
	int maxNfp = mxGetM(prhs[12]);
	double *Nfp = mxGetPr(prhs[13]);
	int P = (int)mxGetScalar(prhs[14]);
	int Np2d = (int)Nfp[Nface - 2];
	double *UpEidM = malloc(Np2d*sizeof(double));
	double *BotEidM = malloc(Np2d*sizeof(double));
	for (int i = 0; i < Np2d; i++){
		UpEidM[i] = Fmask[(Nface - 1)*maxNfp + i];
		BotEidM[i] = Fmask[(Nface - 2)*maxNfp + i];
	}

	GetSparsePatternInVerticalDirection(TempIr, TempJc, Np, Nlayer, Ele2d);
	plhs[0] = mxCreateSparse(Np*Ele3d, Np*Ele3d, TotalNonzero, mxREAL);
	double *sr = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);
	memcpy(irs, TempIr, TotalNonzero*sizeof(mwIndex));
	mwIndex *jcs = mxGetJc(plhs[0]);
	memcpy(jcs, TempJc, (Np*Ele3d + 1)*sizeof(mwIndex));

	double *Tau = malloc(sizeof(double)*(Np2d*Ele2d*(Nlayer + 1)));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < Ele2d; i++){
		CalculatePenaltyParameter(Tau + i*Np2d*(Nlayer + 1), Np2d, Nlayer + 1, P, Nface);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < Ele2d; e++){
		int StartPoint;

		for (int L = 0; L < Nlayer; L++){
			/* Index of the studied element */
			int LocalEle = e*Nlayer + L + 1;
			/* Index of the element that located upside of the studied element */ 
			int UpEle = (int)EToE[e*Nlayer*Nface + L*Nface + Nface - 1];
			/* Index of the element that located downside of the studied element */
			int DownEle = (int)EToE[e*Nlayer*Nface + L*Nface + Nface - 2];

			StartPoint = (int)jcs[e*Nlayer*Np + L*Np];

			GetLocalVolumuIntegralTermForSecondOrderTerm(sr, StartPoint, \
				Np, jcs[e*Nlayer*Np + L*Np + 1] - jcs[e*Nlayer*Np + L*Np], M3d, \
				Dt, tz + e*Nlayer*Np + L*Np, J + e*Nlayer*Np + L*Np, L);
			
			if (UpEle == DownEle){ 
				/* we have only one layer in vertical direction, and L=0 */
				StartPoint = jcs[e*Nlayer*Np + L*Np];

				ImposeSecondOrderNonhydroDirichletBoundaryCondition(sr, StartPoint, Np, \
					Np2d, jcs[e * Nlayer * Np + 1] - jcs[e * Nlayer * Np], M3d, M2d, J + e * Nlayer * Np + L*Np, \
					J2d + e*Np2d, UpEidM, Tau + e*(Nlayer + 1)*Np2d, Dt, tz + e * Nlayer * Np + L * Np);

			}
			else{//two or more layers included in vertical direction
				if (LocalEle == UpEle){//This is the top most cell, and L=0
					StartPoint = (int)jcs[e*Nlayer*Np + L*Np];

					ImposeSecondOrderNonhydroDirichletBoundaryCondition(sr, StartPoint, Np, \
						Np2d, jcs[e * Nlayer * Np + 1] - jcs[e * Nlayer * Np], M3d, M2d, J + e * Nlayer * Np + L*Np, \
						J2d + e*Np2d, UpEidM, Tau + e*(Nlayer + 1)*Np2d, Dt, tz + e * Nlayer * Np + L * Np);


					StartPoint = (int)jcs[e*Nlayer*Np + L*Np] + Np;

					GetLocalToDownElementContributionForSecondOrderTerm(sr, StartPoint, Np, Np2d, jcs[e * Nlayer * Np + L*Np + 1] - jcs[e * Nlayer * Np + L*Np], \
						M3d, M2d, J + e * Nlayer * Np + L*Np + Np, J2d + e*Np2d, BotEidM, UpEidM, Dt, tz + e * Nlayer * Np + L*Np + Np, \
						tz + e * Nlayer * Np + L*Np, Tau + e*(Nlayer + 1)*Np2d + L*Np2d + Np2d);

					GetLocalDownContributionForSecondOrderTerm(sr, StartPoint - Np, Np, Np2d, jcs[e * Nlayer * Np + L*Np + 1] - jcs[e * Nlayer * Np + L*Np], \
						M3d, M2d, J + e * Nlayer * Np + L*Np, J2d + e*Np2d, BotEidM, Dt, tz + e * Nlayer * Np + L*Np, \
						Tau + e*(Nlayer + 1)*Np2d + L*Np2d + Np2d);
				}
				else if (LocalEle == DownEle){
					/* 
					This is the bottom most cell, jcs[e*Nlayer*Np + L*Np + 1] means the total nozero contained in the sparse matrix before the current point, including the 
					 first point of the studied element, (int)jcs[e*Nlayer*Np + L*Np + 1] - Np means the start of the studied local element, while (int)jcs[e*Nlayer*Np + L*Np + 1] - 2*Np
					the start of the element located upside of the studied local element
					*/
					StartPoint = (int)jcs[e*Nlayer*Np + L*Np + 1] - 2*Np;

					GetLocalToUpElementContributionForSecondOrderTerm(sr, StartPoint, Np, Np2d, jcs[e * Nlayer * Np + L*Np + 1] - jcs[e * Nlayer * Np + L*Np], \
						M3d, M2d, J + e * Nlayer * Np + L*Np - Np, J2d + e*Np2d, UpEidM, BotEidM, Dt, tz + e * Nlayer * Np + L*Np - Np, \
						tz + e * Nlayer * Np + L*Np, Tau + e*(Nlayer + 1)*Np2d + L*Np2d);

					GetLocalUpContributionForSecondOrderTerm(sr, StartPoint + Np, Np, Np2d, jcs[e * Nlayer * Np + L*Np + 1] - jcs[e * Nlayer * Np + L*Np], \
						M3d, M2d, J + e * Nlayer * Np + L*Np, J2d + e*Np2d, UpEidM, Dt, tz + e * Nlayer * Np + L*Np, \
						Tau + e*(Nlayer + 1)*Np2d + L*Np2d);
					/**********************************************Newmann Boundary is not considered here***********************************************************/

				}
				else{
					/*
						Element locates in the middle of the studied column, and there are elements located both upside and downside
					*/
					StartPoint = (int)jcs[e*Nlayer*Np + L*Np];

					GetLocalToUpElementContributionForSecondOrderTerm(sr, StartPoint, Np, Np2d, jcs[e * Nlayer * Np + L*Np + 1] - jcs[e * Nlayer * Np + L*Np], \
						M3d, M2d, J + e * Nlayer * Np + L*Np - Np, J2d + e*Np2d, UpEidM, BotEidM, Dt, tz + e * Nlayer * Np + L*Np - Np, \
						tz + e * Nlayer * Np + L*Np, Tau + e*(Nlayer + 1)*Np2d + L*Np2d);

					GetLocalUpContributionForSecondOrderTerm(sr, StartPoint + Np, Np, Np2d, jcs[e * Nlayer * Np + L*Np + 1] - jcs[e * Nlayer * Np + L*Np], \
						M3d, M2d, J + e * Nlayer * Np + L*Np, J2d + e*Np2d, UpEidM, Dt, tz + e * Nlayer * Np + L*Np, \
						Tau + e*(Nlayer + 1)*Np2d + L*Np2d);

					StartPoint = (int)jcs[e*Nlayer*Np + L*Np] + 2*Np;

					GetLocalToDownElementContributionForSecondOrderTerm(sr, StartPoint, Np, Np2d, jcs[e * Nlayer * Np + L*Np + 1] - jcs[e * Nlayer * Np + L*Np], \
						M3d, M2d, J + e * Nlayer * Np + L*Np + Np, J2d + e*Np2d, BotEidM, UpEidM, Dt, tz + e * Nlayer * Np + L*Np + Np, \
						tz + e * Nlayer * Np + L*Np, Tau + e*(Nlayer + 1)*Np2d + L*Np2d + Np2d);

					GetLocalDownContributionForSecondOrderTerm(sr, StartPoint - Np, Np, Np2d, jcs[e * Nlayer * Np + L*Np + 1] - jcs[e * Nlayer * Np + L*Np], \
						M3d, M2d, J + e * Nlayer * Np + L*Np, J2d + e*Np2d, BotEidM, Dt, tz + e * Nlayer * Np + L*Np, \
						Tau + e*(Nlayer + 1)*Np2d + L*Np2d + Np2d);

				}
			}
		}
	}

	free(TempIr);
	free(TempJc);
	free(UpEidM);
	free(BotEidM);
	free(Tau);
}