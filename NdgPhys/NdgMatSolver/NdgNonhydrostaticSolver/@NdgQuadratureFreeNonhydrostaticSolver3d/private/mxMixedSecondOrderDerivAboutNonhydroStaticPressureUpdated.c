#include "SWENonhydrostatic3d.h"

/*This function is called at the initialization stage to calculate the mixed second order derivative
about nonhydrostatic pressure in horizontal direction, i.e. $\frac{\partial}{\partial x}(\frac{\partial q}{\partial \sigma})$
and $\frac{\partial}{\partial y}(\frac{\partial q}{\partial \sigma})$. The following primal form is given for the first term.
For this term, the primal form is given as:
$$-\int_{\Omega}\frac{\partial v}{\partial x}\frac{\partial p}{\partial \sigma}d\Omega + \int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}\\
+ \int_{\epsilon_D}\frac{\partial v}{\partial x}n_{\sigma}pd\Omega + \int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega + \\
\int_{\epsilon_D}\frac{\partial p}{\partial \sigma}n_xvd\Omega - \int_{\epsilon_i} \tau [p]_x[v]_xd\Omega - \int_{\epsilon_D} \tau n_x^2pvd\Omega \\
= \int_{\Omega}fvd\Omega + \int_{\epsilon_D}\frac{\partial v}{\partial x}n_{\sigma}q_Dd\Omega - \int_{\epsilon_D}\tau n_x^2vq_D d\Omega - \int_{\epsilon_N} n_x^2vq_N d\Omega$$*/

/*This following funciton is used to assemble term corresponds to
$$-\int_{\Omega}\frac{\partial v}{\partial x}\frac{\partial p}{\partial \sigma}d\Omega$$
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
Startpoint: Number of nonzeros before the current studied element
in the sparse matrix
Np: Number of interpolation point for the 3d element
NonzeroPerColumn: Number of nonzeros in each column of the studied scope,
i.e. nonzeros in the column corresponds to global index
of the studied point
mass3d: A pointer to the start of the mass matrix for the master cell
Dr: A pointer to the start of the differential matrix for the master cell in direction r
Ds: A pointer to the start of the differential matrix for the master cell in direction s
rd: A pointer to the start of the transformation Jacobian of the studied element in direction r, here d can be x or y
sd: A pointer to the start of the transformation Jacobian of the studied element in direction s, here d can be x or y
tz: A pointer to the start of the transformation Jacobian of the studied element in direction t
J: A pointer to the start of the Jacobian of the studied element
*/

void GetLocalVolumuIntegralTermForMixedSecondOrderTerm(double *dest, int StartPoint, int Np, int NonzeroPerColumn, \
	double *mass3d, double *Dr, double *Ds, double *Dt, double *rd, double *sd, double *tz, double *J){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvMass3d, (ptrdiff_t)Np);
	double *DiffMatrixBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrixBuff, Dr, rd, Np);
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrix, Ds, sd, Np);
	Add(DiffMatrix, DiffMatrix, DiffMatrixBuff, Np*Np);
	double *VerticalDiffMatrix = malloc(Np*Np*sizeof(double));
	DiagMultiply(VerticalDiffMatrix, Dt, tz, Np);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	double *TempContributionBuff = malloc(Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VerticalDiffMatrix, (ptrdiff_t)Np, 0.0, TempContribution, (ptrdiff_t)Np);
	MultiplyByConstant(TempContribution, TempContribution, -1.0, Np*Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(Mass3d);
	free(InvMass3d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(VerticalDiffMatrix);
	free(TempContribution);
	free(TempContributionBuff);
	free(Contribution);
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma} + \\
\int_{\epsilon_D}\frac{\partial v}{\partial x}n_{\sigma}pd\Omega + \\
\int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega + \\
\int_{\epsilon_D}\frac{\partial p}{\partial \sigma}n_xvd\Omega - \\
\int_{\epsilon_i} \tau [p][v]d\Omega - \int_{\epsilon_D} \tau pvd\Omega$$.
For this primal form, both $$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}$$ and 
$$\int_{\epsilon_D}\frac{\partial v}{\partial x}n_{\sigma}pd\Omega$$ equal to zero
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
LocalEle: The studied local element
Nface2d: Number of faces for 2d element
Nface: Number of faces for 3d element
EToE: A pointer to an array that contains the topological relation of local element
mass2d: A pointer to the mass matrix of the 2d latetal master face
mass3d: A pointer to the mass matrix of the master element
Dt: A pointer to the differential matrix in vertical direction of the master element
tz: A pointer to the start of the transformation Jacobian in vertical direction of the local element
Js: A pointer to the facial Jacobian of the lateral face
J: A pointer to the Jacobian of the studied mesh
TempEToE: A pointer to the sorted arry for the topological relation of the studied element
FToE: A pointer to the start of the topological relation for face to element relation of the inner edge
FToF: A pointer to the start of the topological relation for face to face relation of the inner edge
Fmask: A pointer to the start of the facial interpolation point index
maxNfp: Maximum number of the  facial interpolation point
Nfp: Number of facial interpolation point for lateral face
Ne: Number of edges for inner edge
Np: Number of interpolation point for the master cell
Vector: A pointer to the start of the facial direction vector of the inner edge
Jcs: A pointer to the array that contains the number of nonzero element in each column
*/

void GetLocalToAdjacentElementContributionInHorizontalDirection(double *dest, int LocalEle, int *AdjEle, double *FacialVector, int *GlobalFace,\
	int InternalFace, int Nface, double *mass2d, double *mass3d, double *Dt, double *tz, double *Js, double *J, double *TempEToE, \
	int *Flag, double *FToN1, double *FToN2, int Np, int Nfp, double *Tau, mwIndex *Jcs){

	double *LocalEidM = malloc(Nfp*sizeof(double));
	double *AdjEidM = malloc(Nfp*sizeof(double));
	for (int i = 0; i < InternalFace; i++){
		if (Flag[i] == 0){//The local and adjacent facial point are the same with the data stored in FToN1 and FToN2
			for (int p = 0; p < Nfp; p++){
				LocalEidM[p] = FToN1[GlobalFace[i] * Nfp + p];
				AdjEidM[p] = FToN2[GlobalFace[i] * Nfp + p];
			}
		}
		else{
			for (int p = 0; p < Nfp; p++){//The local and adjacent facial point are reversed in FToN1 and FToN2
				AdjEidM[p] = FToN1[GlobalFace[i] * Nfp + p];
				LocalEidM[p] = FToN2[GlobalFace[i] * Nfp + p];
			}
		}

		double *LocalDiff = malloc(Np*Np*sizeof(double));
		DiagMultiply(LocalDiff, Dt, tz + (LocalEle - 1)*Np, Np);

		double *AdjMass3d = malloc(Np*Np*sizeof(double));
		DiagMultiply(AdjMass3d, mass3d, J + Np*(AdjEle[i] - 1), Np);
		double *InvAdjMass3d = malloc(Np*Np*sizeof(double));
		memcpy(InvAdjMass3d, AdjMass3d, Np*Np*sizeof(double));
		MatrixInverse(InvAdjMass3d, (ptrdiff_t)Np);
		double *Mass2d = malloc(Nfp*Nfp*sizeof(double));
		DiagMultiply(Mass2d, mass2d, Js + Nfp*GlobalFace[i], Nfp);
		double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));

		double *TempContribution = malloc(Np*Np*sizeof(double));
		memset(TempContribution, 0, Np*Np*sizeof(double));
		double *Contribution = malloc(Np*Np*sizeof(double));

		double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));
		// $$\int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega$$
		double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
		AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEidM, Nfp, Np);
		DiagMultiply(EleMass2d, Mass2d, FacialVector + i*Nfp, Nfp);

		MultiplyByConstant(EleMass2d, EleMass2d, -0.5, Nfp*Nfp);

		MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, EleMass2d,
			(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

		AssembleContributionIntoRow(TempContribution, EdgeContribution, AdjEidM, Np, Nfp);

		/*We delete this part, because $n_{\sigma}$ is zero on the  lateral face*/
		/*$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}$$*/
/*		
        DiagMultiply(Mass2d, mass2d, Js + Nfp*GlobalFace[i], Nfp);
		DiagMultiply(EleMass2d, Mass2d, FacialVector + i*Nfp, Nfp);
		MultiplyByConstant(EleMass2d, EleMass2d, 0.5, Nfp*Nfp);
		MatrixMultiply("T", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,\
			(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);
		AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEidM, Np, Nfp);
*/

		/*$$-\tau\left (  q^--q^+, v^--v^+\right)_{\epsilon_i}$$*/
		DiagMultiply(Mass2d, mass2d, Js + Nfp*GlobalFace[i], Nfp);
		DiagMultiply(Mass2d, Mass2d, Tau + Nfp*GlobalFace[i], Nfp);
		AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, AdjEidM, LocalEidM, Np, Nfp, 1.0);

		MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvAdjMass3d,
			(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

		for (int j = 0; j < Nface + 1; j++){
			if (AdjEle[i] == TempEToE[j]){
				int StartPoint = Jcs[(LocalEle - 1)*Np] + j*Np;
				int NonzeroPerColumn = Jcs[(LocalEle - 1)*Np + 1] - Jcs[(LocalEle - 1)*Np];
				AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
				break;
			}
		}
		free(LocalDiff);
		free(AdjMass3d);
		free(InvAdjMass3d);
		free(Mass2d);
		free(EleMass2d);
		free(TempContribution);
		free(Contribution);
		free(FacialDiffMatrix);
		free(EdgeContribution);
	}
	free(LocalEidM);
	free(AdjEidM);
}


/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma} + \\
\int_{\epsilon_D}\frac{\partial v}{\partial x}n_{\sigma}pd\Omega + \\
\int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega + \\
\int_{\epsilon_D}\frac{\partial p}{\partial \sigma}n_xvd\Omega - \\
\int_{\epsilon_i} \tau [p][v]d\Omega - \int_{\epsilon_D} \tau pvd\Omega$$.
For this primal form, both $$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}$$ and
$$\int_{\epsilon_D}\frac{\partial v}{\partial x}n_{\sigma}pd\Omega$$ equal to zero
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
LocalEle: The studied local element
Nface2d: Number of faces for 2d element
Nface: Number of faces for 3d element
EToE: A pointer to an array that contains the topological relation of local element
mass2d: A pointer to the mass matrix of the 2d latetal master face
mass3d: A pointer to the mass matrix of the master element
Dt: A pointer to the differential matrix in vertical direction of the master element
tz: A pointer to the start of the transformation Jacobian in vertical direction of the local element
Js: A pointer to the facial Jacobian of the lateral face
J: A pointer to the Jacobian of the studied mesh
TempEToE: A pointer to the sorted arry for the topological relation of the studied element
FToE: A pointer to the start of the topological relation for face to element relation of the inner edge
FToF: A pointer to the start of the topological relation for face to face relation of the inner edge
Fmask: A pointer to the start of the facial interpolation point index
maxNfp: Maximum number of the  facial interpolation point
Nfp: Number of facial interpolation point for lateral face
Ne: Number of edges for inner edge
Np: Number of interpolation point for the master cell
Vector: A pointer to the start of the facial direction vector of the inner edge
Jcs: A pointer to the array that contains the number of nonzero element in each column
*/

void GetLocalFacialContributionInHorizontalDirection(double *dest, int LocalEle, int StartPoint, double *Vector, int *GlobalFace, \
	int InternalFace, double *mass2d, double *mass3d, double *Dt, double *tz, double *Js, double *J, int *Flag, double *FToN1, \
	double *FToN2, int Np, int Nfp, double *Tau, mwIndex *Jcs){

	double *LocalEidM = malloc(Nfp*sizeof(double));
	for (int i = 0; i < InternalFace; i++){
		if (Flag[i] == 0){
			for (int p = 0; p < Nfp; p++){
				LocalEidM[p] = FToN1[GlobalFace[i] * Nfp + p];
			}
		}
		else{
			for (int p = 0; p < Nfp; p++){
				LocalEidM[p] = FToN2[GlobalFace[i] * Nfp + p];
			}
		}

		double *LocalDiff = malloc(Np*Np*sizeof(double));
		DiagMultiply(LocalDiff, Dt, tz + (int)(LocalEle - 1)*Np, Np);

		double *LocalMass3d = malloc(Np*Np*sizeof(double));
		DiagMultiply(LocalMass3d, mass3d, J + Np*(LocalEle - 1), Np);
		double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
		memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
		MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
		double *Mass2d = malloc(Nfp*Nfp*sizeof(double));
		DiagMultiply(Mass2d, mass2d, Js + Nfp*GlobalFace[i], Nfp);
		double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));

		double *TempContribution = malloc(Np*Np*sizeof(double));
		memset(TempContribution, 0, Np*Np*sizeof(double));
		double *Contribution = malloc(Np*Np*sizeof(double));

		double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));
		// $$\int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega$$
		double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
		AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEidM, Nfp, Np);
		DiagMultiply(Mass2d, Mass2d, Vector + i*Nfp, Nfp);

		MultiplyByConstant(Mass2d, Mass2d, 0.5, Nfp*Nfp);
		
		MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, Mass2d,
			(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

		AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEidM, Np, Nfp);

		/***********We delete this part because nz = 0**********************/
		/*$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}$$*/
//		DiagMultiply(Mass2d, mass2d, Js + Nfp*GlobalFace[i], Nfp);
//		DiagMultiply(EleMass2d, Mass2d, Vector + i*Nfp, Nfp);
//		MultiplyByConstant(EleMass2d, EleMass2d, 0.5, Nfp*Nfp);
//		MatrixMultiply("T", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,\
			(ptrdiff_t)Nfp, EleMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);
//		AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEidM, Np, Nfp);

		//$$-\tau\left (  q^--q^+, v^--v^+\right)_{\epsilon_i}$$
		DiagMultiply(Mass2d, mass2d, Js + Nfp*GlobalFace[i], Nfp);
		DiagMultiply(Mass2d, Mass2d, Tau + Nfp*GlobalFace[i], Nfp);
		AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, LocalEidM, LocalEidM, Np, Nfp, -1.0);

		MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
			(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

		AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, Jcs[(LocalEle-1)*Np + 1] - Jcs[(LocalEle-1)*Np], Np);

		free(LocalDiff);
		free(LocalMass3d);
		free(InvLocalMass3d);
		free(Mass2d);
		free(EleMass2d);
		free(TempContribution);
		free(Contribution);
		free(FacialDiffMatrix);
		free(EdgeContribution);
	}
	free(LocalEidM);
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_D}q_h\nabla_x vn_{\sigma}d\boldsymbol{x} + \\
\int_{\epsilon_D}\nabla_{\sigma} q_h vn_xd\boldsymbol{x} - \\
\int_{\epsilon_D}\tau vq_h$$, here, the second term is zero since $n_h = 0$ on boundary.
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
Startpoint: Number of nonzeros before the current studied element in the sparse matrix, here it's the point to write the data
Np: Number of interpolation point for the 3d element
Np2d: Number of interpolation point for the 2d element. For the current case, the cell on top of the 3d element.
NonzeroPerColumn: Number of nonzeros in each column of the studied scope, i.e. nonzeros in the column corresponds to global index of the studied point
mass3d: A pointer to the start of the mass matrix for the master cell
mass2d: A pointer to the start of the mass matrix for the 2 dimensional master cell
J: A pointer to the start of the Jacobian of the studied element
J2d: A pointer to the start of the Jacobian of the 2 dimensional studied element
Eid: A pointer to the start of the interpolation point index on the top of the master cell
*/

void ImposeSecondOrderNonhydroDirichletBoundaryCondition(double *dest, int StartPoint, int Np, int Np2d, \
	int NonzeroPerColumn, double *mass3d, double *mass2d, double *J, double *J2d, double *Eid, double *Dr, \
	double *Ds, double *rd, double *sd, double *Tau){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *DiffMatrixBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrixBuff, Dr, rd, Np);
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrix, Ds, sd, Np);
	Add(DiffMatrix, DiffMatrix, DiffMatrixBuff, Np*Np);


	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, Eid, Np2d, Np);
	//$$\int_{\epsilon_D}q_h\nabla_x vn_{\sigma}d\boldsymbol{x}$$
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, Eid, Np, Np2d);
	//$$-\tau\left (  q^--q^+, v^--v^+\right)_{\epsilon_i}$$
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	DiagMultiply(Mass2d, Mass2d, Tau, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, Eid, Eid, Np, Np2d, -1.0);

	/*Multiply the contribution by inverse matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(Mass3d);
	free(InvMass3d);
	free(Mass2d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_{bh}}\nabla_vq_hvn_hd\boldsymbol{x}$$
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
LocalEle: The studied local element
Nface2d: Number of faces for 2d element
Nface: Number of faces for 3d element
EToE: A pointer to an array that contains the topological relation of local element
mass2d: A pointer to the mass matrix of the 2d latetal master face
mass3d: A pointer to the mass matrix of the master element
Dt: A pointer to the differential matrix in vertical direction of the master element
tz: A pointer to the start of the transformation Jacobian in vertical direction of the local element
Js: A pointer to the facial Jacobian of the lateral face
J: A pointer to the Jacobian of the studied mesh
TempEToE: A pointer to the sorted arry for the topological relation of the studied element
FToE: A pointer to the start of the topological relation for face to element relation of the inner edge
FToF: A pointer to the start of the topological relation for face to face relation of the inner edge
Fmask: A pointer to the start of the facial interpolation point index
maxNfp: Maximum number of the  facial interpolation point
Nfp: Number of facial interpolation point for lateral face
Ne: Number of edges for inner edge
Np: Number of interpolation point for the master cell
Vector: A pointer to the start of the facial direction vector of the inner edge
Jcs: A pointer to the array that contains the number of nonzero element in each column
*/
//At present, we impose Neumann boundary here, and no Dirichlet boundary considered.
void ImposeBoundaryConditionInHorizontalDirection(double *dest, int LocalEle, int Nface2d, int Nface, \
	double *EToE, double *mass2d, double *mass3d, double *Dt, double *tz, double *Js, double *J, double *TempEToE, \
	double *FToE, double *FToF, double *Fmask, int maxNfp, int Nfp, int Ne, int Np, double *Vector, mwIndex *Jcs){

	for (int i = 0; i < Nface2d; i++){
		if (LocalEle == (int)EToE[i]){
			int GlobalFace = 0, LocalFace  = 0;
			double *FacialVector = malloc(Nfp * sizeof(double));

			FindFaceAndDirectionVectorAtBoundary(FacialVector, &GlobalFace, &LocalFace, Nfp, LocalEle, FToE, FToF, i + 1, \
				Vector, Ne);

			double *LocalEidM = malloc(Nfp*sizeof(double));

			for (int p = 0; p < Nfp; p++){
				LocalEidM[p] = Fmask[(LocalFace - 1)*maxNfp + p];
			}

			double *LocalDiff = malloc(Np*Np*sizeof(double));

			DiagMultiply(LocalDiff, Dt, tz + (int)(LocalEle - 1)*Np, Np);

			double *LocalMass3d = malloc(Np*Np*sizeof(double));
			DiagMultiply(LocalMass3d, mass3d, J + Np*(int)(EToE[i] - 1), Np);
			double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
			memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
			MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
			double *Mass2d = malloc(Nfp*Nfp*sizeof(double));
			DiagMultiply(Mass2d, mass2d, Js + Nfp*GlobalFace, Nfp);

			double *TempContribution = malloc(Np*Np*sizeof(double));
			memset(TempContribution, 0, Np*Np*sizeof(double));
			double *Contribution = malloc(Np*Np*sizeof(double));

			double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

			double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
			AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEidM, Nfp, Np);
			DiagMultiply(Mass2d, Mass2d, FacialVector, Nfp);

			MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, Mass2d,
				(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

			AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEidM, Np, Nfp);

			MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
				(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

			for (int j = 0; j < Nface + 1; j++){
				if (LocalEle == (int)TempEToE[j] && (j+1) == LocalFace){
					int StartPoint = Jcs[(LocalEle - 1)*Np] + j*Np;
					int NonzeroPerColumn = Jcs[(LocalEle - 1)*Np + 1] - Jcs[(LocalEle - 1)*Np];
					AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
					break;
				}
			}

			free(FacialVector);
			free(LocalEidM);
			free(LocalDiff);
			free(LocalMass3d);
			free(InvLocalMass3d);
			free(Mass2d);
			free(TempContribution);
			free(Contribution);
			free(FacialDiffMatrix);
			free(EdgeContribution);

		}
		else{
			continue;
		}
	}
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma} + \\
\int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega - \\
\int_{\epsilon_i} \tau [p][v]d\Omega$$.
For this part, the second term equals to zero, since $n_x=0$ on the horizontal face
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
Np: Number of interpolation point for the 3d element
Np2d: Number of interpolation point for the 2d element. For the current case, the cell on top of the 3d element.
NonzeroPerColumn: Number of nonzeros in each column of the studied scope, i.e. nonzeros in the column corresponds to global index of the studied point
mass3d: A pointer to the start of the mass matrix for the master cell
mass2d: A pointer to the start of the mass matrix for the 2 dimensional master cell
J: A pointer to the start of the Jacobian of the studied element
J2d: A pointer to the start of the Jacobian of the 2 dimensional studied element
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
BotEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
Dr: A pointer to the start of the differential matrix for the master cell in direction r
Ds: A pointer to the start of the differential matrix for the master cell in direction s
rd: A pointer to the start of the transformation Jacobian of the studied element in direction r, here d can be x or y
sd: A pointer to the start of the transformation Jacobian of the studied element in direction s, here d can be x or y
*/

void GetLocalToDownElementContributionForMixedSecondOrderTerm(double *dest, double *EToE, double *TempEToE, double *FToE, double *FToF, int Nface, int LocalEle, \
	int Np, int Np2d, int Ne, int NonzeroPerColumn, double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *BotEid, double *Dr,\
	double *Ds, double *rd, double *sd, double *Tau, mwIndex *Jcs){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvDownMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvDownMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvDownMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *DiffMatrixBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrixBuff, Dr, rd, Np);
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrix, Ds, sd, Np);
	Add(DiffMatrix, DiffMatrix, DiffMatrixBuff, Np*Np);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	/*$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}$*/
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, BotEid, Np2d, Np);
	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(EdgeContribution, EdgeContribution, -0.5, Np2d*Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*$\int_{\epsilon_i} \tau [p][v]d\Omega$*/
	int GlobalFace = 0;
	FindGlobalBottomEdgeFace(&GlobalFace, FToE, LocalEle, (int)EToE[Nface-2], Ne);
	DiagMultiply(Mass2d, Mass2d, Tau + Np2d*GlobalFace, Np2d);
	//Local to adjacent, positive, while local to local, minus.
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, BotEid, LocalEid, Np, Np2d, 1.0);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvDownMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	for (int i = 0; i < (Nface+1); i++){
		if (EToE[Nface - 2] == TempEToE[i]){
			int StartPoint = Jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}
	
	free(Mass3d);
	free(InvDownMass3d);
	free(Mass2d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma} + \\
\int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega - \\
\int_{\epsilon_i} \tau [p][v]d\Omega$$.
For this part, the second term equals to zero, since $n_x=0$ on the horizontal face
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
Np: Number of interpolation point for the 3d element
Np2d: Number of interpolation point for the 2d element. For the current case, the cell on top of the 3d element.
NonzeroPerColumn: Number of nonzeros in each column of the studied scope, i.e. nonzeros in the column corresponds to global index of the studied point
mass3d: A pointer to the start of the mass matrix for the master cell
mass2d: A pointer to the start of the mass matrix for the 2 dimensional master cell
J: A pointer to the start of the Jacobian of the studied element
J2d: A pointer to the start of the Jacobian of the 2 dimensional studied element
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
BotEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
Dr: A pointer to the start of the differential matrix for the master cell in direction r
Ds: A pointer to the start of the differential matrix for the master cell in direction s
rd: A pointer to the start of the transformation Jacobian of the studied element in direction r, here d can be x or y
sd: A pointer to the start of the transformation Jacobian of the studied element in direction s, here d can be x or y
*/

void GetLocalDownFacialContributionForMixedSecondOrderTerm(double *dest, double *EToE, double *TempEToE, double *FToE, double *FToF, int Nface, int LocalEle, \
	int StartPoint, int NonzeroPerColumn, int Np, int Np2d, int Ne, double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, \
	double *Dr, double *Ds, double *rd, double *sd, double *Tau, mwIndex *Jcs){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvLocalMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *DiffMatrixBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrixBuff, Dr, rd, Np);
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrix, Ds, sd, Np);
	Add(DiffMatrix, DiffMatrix, DiffMatrixBuff, Np*Np);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	/*$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}$*/
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(EdgeContribution, EdgeContribution, -0.5, Np2d*Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*$\int_{\epsilon_i} \tau [p][v]d\Omega$*/
	int GlobalFace = 0;
	FindGlobalBottomEdgeFace(&GlobalFace, FToE, LocalEle, (int)EToE[Nface - 2], Ne);
	DiagMultiply(Mass2d, Mass2d, Tau + Np2d*GlobalFace, Np2d);
	//Local to adjacent, positive, while local to local, minus.
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(Mass3d);
	free(InvLocalMass3d);
	free(Mass2d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma} + \\
\int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega - \\
\int_{\epsilon_i} \tau [p][v]d\Omega$$.
For this part, the second term equals to zero, since $n_x=0$ on the horizontal face.
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
Np: Number of interpolation point for the 3d element
Np2d: Number of interpolation point for the 2d element. For the current case, the cell on top of the 3d element.
NonzeroPerColumn: Number of nonzeros in each column of the studied scope, i.e. nonzeros in the column corresponds to global index of the studied point
mass3d: A pointer to the start of the mass matrix for the master cell
mass2d: A pointer to the start of the mass matrix for the 2 dimensional master cell
J: A pointer to the start of the Jacobian of the studied element
J2d: A pointer to the start of the Jacobian of the 2 dimensional studied element
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
UpEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
Dr: A pointer to the start of the differential matrix for the master cell in direction r
Ds: A pointer to the start of the differential matrix for the master cell in direction s
rd: A pointer to the start of the transformation Jacobian of the studied element in direction r, here d can be x or y
sd: A pointer to the start of the transformation Jacobian of the studied element in direction s, here d can be x or y
*/

void GetLocalToUpElementContributionForMixedSecondOrderTerm(double *dest, double *EToE, double *TempEToE, double *FToE, double *FToF, int Nface, int LocalEle, \
int Np, int Np2d, int Ne, int NonzeroPerColumn, double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *UpEid, double *Dr, \
double *Ds, double *rd, double *sd, double *Tau, mwIndex *Jcs){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvUpMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvUpMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvUpMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *DiffMatrixBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrixBuff, Dr, rd, Np);
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrix, Ds, sd, Np);
	Add(DiffMatrix, DiffMatrix, DiffMatrixBuff, Np*Np);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	/*$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}$*/
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Np2d, Np);
	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np2d*Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*$\int_{\epsilon_i} \tau [p][v]d\Omega$*/
	int GlobalFace = 0;
	FindGlobalBottomEdgeFace(&GlobalFace, FToE, LocalEle, (int)EToE[Nface - 1], Ne);
	DiagMultiply(Mass2d, Mass2d, Tau + Np2d*GlobalFace, Np2d);
	//Local to adjacent, positive, while local to local, minus.
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, UpEid, LocalEid, Np, Np2d, 1.0);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvUpMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	for (int i = 0; i < (Nface + 1); i++){
		if (EToE[Nface - 1] == TempEToE[i]){
			int StartPoint = Jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}
	
	free(Mass3d);
	free(InvUpMass3d);
	free(Mass2d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma} + \\
\int_{\epsilon_i}\left\{\frac{\partial p}{\partial \sigma}\right\}[v]_xd\Omega - \\
\int_{\epsilon_i} \tau [p][v]d\Omega$$.
For this part, the second term equals to zero, since $n_x=0$ on the horizontal face.
The input parameter are as follows:
dest: A pointer to the start of the sparse matrix
Np: Number of interpolation point for the 3d element
Np2d: Number of interpolation point for the 2d element. For the current case, the cell on top of the 3d element.
NonzeroPerColumn: Number of nonzeros in each column of the studied scope, i.e. nonzeros in the column corresponds to global index of the studied point
mass3d: A pointer to the start of the mass matrix for the master cell
mass2d: A pointer to the start of the mass matrix for the 2 dimensional master cell
J: A pointer to the start of the Jacobian of the studied element
J2d: A pointer to the start of the Jacobian of the 2 dimensional studied element
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
Dr: A pointer to the start of the differential matrix for the master cell in direction r
Ds: A pointer to the start of the differential matrix for the master cell in direction s
rd: A pointer to the start of the transformation Jacobian of the studied element in direction r, here d can be x or y
sd: A pointer to the start of the transformation Jacobian of the studied element in direction s, here d can be x or y
*/

void GetLocalUpFacialContributionForMixedSecondOrderTerm(double *dest, double *EToE, double *TempEToE, double *FToE, double *FToF, int Nface, \
	int LocalEle, int StartPoint, int NonzeroPerColumn, int Np, int Np2d, int Ne, double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, \
	double *Dr, double *Ds, double *rd, double *sd, double *Tau, mwIndex *Jcs){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvLocalMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *DiffMatrixBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrixBuff, Dr, rd, Np);
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrix, Ds, sd, Np);
	Add(DiffMatrix, DiffMatrix, DiffMatrixBuff, Np*Np);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	/*$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial x}\right\}[p]_{\sigma}$*/
	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);
	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np2d*Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*$\int_{\epsilon_i} \tau [p][v]d\Omega$*/
	int GlobalFace = 0;
	FindGlobalBottomEdgeFace(&GlobalFace, FToE, LocalEle, (int)EToE[Nface - 1], Ne);
	DiagMultiply(Mass2d, Mass2d, Tau + Np2d*GlobalFace, Np2d);
	//Local to adjacent, positive, while local to local, minus.
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);


	free(Mass3d);
	free(InvLocalMass3d);
	free(Mass2d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}

/*The input parameters are organized as follows:
Np: Number of interpolation points, indexed as 0.
Ele3d: Number of elements in three-dimensional mesh object, indexed as 1.
Nlayer: Number of layers in vertical direction, indexed as 2.
Ele2d: Number of elements in two-dimensional mesh object, indexed as 3.
Nface: Number of face numbers for the three-dimensional master cell, indexed as 4.
EToE: Element to element topological relation of the 3d mesh, indexed as 5.
FToE: Face to element topological relation of the inner edge, indexed as 6.
FToF: Face to face topological relation of the inner edge, indexed as 7.
J: Jacobian coefficient of the three-dimensional mesh object, indexed as 8.
J2d: Jacobian coefficient of the two-dimensional mesh object, indexed as 9.
Js: Jacobian coefficient of the inner edge, indexed as 10.
M3d: Mass matrix of the 3d master cell, indexed as 11.
LM2d: Mass matrix of the 2d master cell of the inner edge, indexed as 12.
M2d: Mass matrix of the 2d master cell, indexed as 13.
BENe: Number of boundary edge, indexed as 14.
IENe: Number of inner edge, indexed as 15.
Fmask: Index of the interpolation point on each face of the master cell, indexed as 16.
Nfp: Number of interpolation point on each face of the master cell, indexed as 17.
P: Approximation order of the master cell, indexed as 18.
vector: The direction vector of the inner edge object, indexed as 19.
rd: Differential coefficient of the 3d mesh object in direction r, indexed as 20.
sd: Differential coefficient of the 3d mesh object in direction s, indexed as 21.
tz: Jacobian coefficient in vertical direction, indexed as 22.
Dr: Differential matrix of the 3d master cell in direction r, indexed as 23.
Ds: Differential matrix of the 3d master cell in direction s, indexed as 24.
Dt: Differential matrix of the 3d master cell in direction t, indexed as 25.
BEFToF: Face to face topological relation of the boundary edge, indexed as 26.
BEFToE: Face to element topological relation of the inner edge, indexed as 27.
BEVector: The direction vector of the boundary edge object, indexed as 28.
*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int Np = (int)mxGetScalar(prhs[0]);
	int Ele3d = (int)mxGetScalar(prhs[1]);
	int Nlayer = (int)mxGetScalar(prhs[2]);
	int Ele2d = (int)mxGetScalar(prhs[3]);
	int Nface = (int)mxGetScalar(prhs[4]);
	int Nface2d = Nface - 2;
	double *EToE = mxGetPr(prhs[5]);
	double *FToE = mxGetPr(prhs[6]);
	double *FToF = mxGetPr(prhs[7]);
	double *FToN1 = mxGetPr(prhs[8]);
	double *FToN2 = mxGetPr(prhs[9]);

	double *LAV = mxGetPr(prhs[10]);
	double *J = mxGetPr(prhs[11]);
	double *J2d = mxGetPr(prhs[12]);
	double *Js = mxGetPr(prhs[13]);
	double *M3d = mxGetPr(prhs[14]);
	double *LM2d = mxGetPr(prhs[15]);
	double *M2d = mxGetPr(prhs[16]);
	int BENe = mxGetScalar(prhs[17]);
	int IENe = mxGetScalar(prhs[18]);
	int TotalNonzero = Ele3d * (Nface + 1) * Np*Np - BENe * Np*Np - 2 * Ele2d*Np*Np;
	mwIndex *TempIr = malloc(TotalNonzero*sizeof(mwIndex));
	mwIndex *TempJc = malloc((Np*Ele3d + 1)*sizeof(mwIndex));
	TempJc[0] = 0;
	GetSparsePatternInHorizontalDirection(TempIr, TempJc, EToE, Nface, Nface, Ele3d, Np);

	double *IELAV = mxGetPr(prhs[19]);
	double *Fmask = mxGetPr(prhs[20]);
	int maxNfp = mxGetM(prhs[20]);
	double *Nfp = mxGetPr(prhs[21]);
	int HorNfp = (int)Nfp[0];
	int VertNfp = (int)Nfp[Nface - 1];
	double *UpEidM = malloc(VertNfp*sizeof(double));
	double *BotEidM = malloc(VertNfp*sizeof(double));
	for (int i = 0; i < VertNfp; i++){
		UpEidM[i] = Fmask[(Nface - 1)*maxNfp + i];
		BotEidM[i] = Fmask[(Nface - 2)*maxNfp + i];
	}

	int P = (int)mxGetScalar(prhs[22]);
	double *Vector = mxGetPr(prhs[23]);
	double *rd = mxGetPr(prhs[24]);
	double *sd = mxGetPr(prhs[25]);
	double *tz = mxGetPr(prhs[26]);
	double *Dr = mxGetPr(prhs[27]);
	double *Ds = mxGetPr(prhs[28]);
	double *Dt = mxGetPr(prhs[29]);

	double *BEFToF = mxGetPr(prhs[30]);
	double *BEFToE = mxGetPr(prhs[31]);
	double *BEVector = mxGetPr(prhs[32]);
	double *BEJs = mxGetPr(prhs[33]);

	double *BotELAV = mxGetPr(prhs[34]);
	int BotENe = (int)mxGetScalar(prhs[35]);
	double *BotEFToE = mxGetPr(prhs[36]);
	double *BotEFToF = mxGetPr(prhs[37]);

	double *SurfELAV = mxGetPr(prhs[38]);
	int SurfENe = (int)mxGetScalar(prhs[39]);
	double *SurfEFToE = mxGetPr(prhs[40]);
	
	plhs[0] = mxCreateSparse(Np*Ele3d, Np*Ele3d, TotalNonzero, mxREAL);
	double *sr = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);
	memcpy(irs, TempIr, TotalNonzero*sizeof(mwIndex));
	mwIndex *jcs = mxGetJc(plhs[0]);
	memcpy(jcs, TempJc, (Np*Ele3d + 1)*sizeof(mwIndex));

	double *Tau = malloc(IENe*HorNfp*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		double localRatio = IELAV[face] / LAV[(int)FToE[face * 2] - 1];
		double adjacentRatio = IELAV[face] / LAV[(int)FToE[face * 2 + 1] - 1];
		for (int p = 0; p < HorNfp; p++){
			Tau[face*HorNfp + p] =  max(localRatio*(P + 1)*(P + 3) / 3 * Nface / 2, \
				adjacentRatio*(P + 1)*(P + 3) / 3 * Nface / 2);
		}
	}

	double *BotETau = malloc(BotENe*VertNfp*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++){
		double localRatio = BotELAV[face] / LAV[(int)BotEFToE[face * 2] - 1];
		double adjacentRatio = BotELAV[face] / LAV[(int)BotEFToE[face * 2 + 1] - 1];
		for (int p = 0; p < VertNfp; p++){
			BotETau[face*VertNfp + p] = max(localRatio*(P + 1)*(P + 3) / 3 * Nface / 2, \
				adjacentRatio*(P + 1)*(P + 3) / 3 * Nface / 2);
		}
	}

	double *SurfETau = malloc(SurfENe*VertNfp*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfENe; face++){
		double localRatio = SurfELAV[face] / LAV[(int)SurfEFToE[face * 2] - 1];
		double adjacentRatio = SurfELAV[face] / LAV[(int)SurfEFToE[face * 2 + 1] - 1];
		for (int p = 0; p < VertNfp; p++){
			SurfETau[face*VertNfp + p] = max(localRatio*(P + 1)*(P + 3) / 3 * Nface / 2, \
				adjacentRatio*(P + 1)*(P + 3) / 3 * Nface / 2);
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int ele = 0; ele < Ele2d; ele++){
		int StartPoint;
		int LocalStartPoint;
		for (int L = 0; L < Nlayer; L++){
			//Index of the studied element
			int LocalEle = ele*Nlayer + L + 1;
			//Index of the element that located upside of the studied element 
			int UpEle = (int)EToE[ele*Nlayer*Nface + L*Nface + Nface - 1];
			//Index of the element that located downside of the studied element
			int DownEle = (int)EToE[ele*Nlayer*Nface + L*Nface + Nface - 2];

			int EleNumber = 0;

			double *TempEToE = malloc((Nface + 1)*sizeof(double));

			FindUniqueElementAndSortOrder(TempEToE, EToE + (LocalEle - 1)*Nface, &EleNumber, Nface, LocalEle);

			for (int i = 0; i < EleNumber; i++){
				if (LocalEle == TempEToE[i]){
					LocalStartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
					GetLocalVolumuIntegralTermForMixedSecondOrderTerm(sr, LocalStartPoint, Np, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], \
						M3d, Dr, Ds, Dt, rd + (LocalEle - 1)*Np, sd + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, J + (LocalEle - 1)*Np);
					break;
				}
			}

			double *FacialVector = malloc(Nface2d*HorNfp*sizeof(double));
			int *GlobalFace = malloc(Nface2d*sizeof(int));
			int *AdjEle = malloc(Nface2d*sizeof(int));
			int *ReverseFlag = malloc(Nface2d*sizeof(int));
			int InternalFace = 0;
			
			FindFaceAndDirectionVector(FacialVector, GlobalFace, AdjEle, \
				&InternalFace, ReverseFlag, HorNfp, LocalEle, FToE, FToF, Vector, IENe, Nface2d);

			GetLocalToAdjacentElementContributionInHorizontalDirection(sr, LocalEle, AdjEle, FacialVector, GlobalFace, InternalFace, \
				Nface, LM2d, M3d, Dt, tz, Js, J, TempEToE, ReverseFlag,  FToN1, FToN2, Np, HorNfp, Tau, jcs);

			GetLocalFacialContributionInHorizontalDirection(sr, LocalEle, LocalStartPoint, FacialVector, GlobalFace, InternalFace, \
				LM2d, M3d, Dt, tz, Js, J, ReverseFlag, FToN1, FToN2, Np, HorNfp, Tau, jcs);
			
//			ImposeBoundaryConditionInHorizontalDirection(sr, LocalEle, Nface2d, Nface, \
				EToE + (LocalEle - 1)*Nface, LM2d, M3d, Dt, tz, BEJs, J, TempEToE, \
				BEFToE, BEFToF, Fmask, maxNfp, HorNfp, BENe, Np, BEVector, jcs);

			if (UpEle == DownEle){//Only one layer in the vertical direction, impose the Dirichlet boundary condition
				for (int i = 0; i < EleNumber; i++){
					if (LocalEle == TempEToE[i]){
						StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
						ImposeSecondOrderNonhydroDirichletBoundaryCondition(sr, StartPoint, Np, VertNfp, \
							jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, \
							Dr, Ds, rd + ele * Nlayer * Np + L*Np, sd + ele * Nlayer * Np + L*Np, SurfETau + ele*VertNfp);
						break;
					}
				}
			}
			else{//two or more layers included in vertical direction
				if (LocalEle == UpEle){//This is the top most cell, and L=0

					for (int i = 0; i < EleNumber; i++){
						if (LocalEle == TempEToE[i]){
							StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
							ImposeSecondOrderNonhydroDirichletBoundaryCondition(sr, StartPoint, Np, VertNfp, \
								jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, \
								Dr, Ds, rd + ele * Nlayer * Np + L*Np, sd + ele * Nlayer * Np + L*Np, SurfETau + ele*VertNfp);
						}
					}

					GetLocalToDownElementContributionForMixedSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, BotEFToE, BotEFToF, Nface,\
						LocalEle, Np, VertNfp, BotENe, jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + ele * Nlayer * Np + L*Np + Np,\
						J2d + ele*VertNfp, BotEidM, UpEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np + Np, sd + ele * Nlayer * Np + L*Np + Np, BotETau, jcs);

					GetLocalDownFacialContributionForMixedSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, BotEFToE, BotEFToF, Nface, LocalEle, \
						LocalStartPoint, jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], Np, VertNfp, BotENe, M3d, M2d, J + ele * Nlayer * Np + L*Np,\
						J2d + ele*VertNfp, BotEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np, sd + ele * Nlayer * Np + L*Np, BotETau, jcs);

				}
				else if (LocalEle == DownEle){// This is the bottom most cell.

					GetLocalToUpElementContributionForMixedSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, BotEFToE, BotEFToF, Nface,\
						LocalEle, Np, VertNfp, BotENe, jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + ele * Nlayer * Np + L*Np - Np, \
						J2d + ele*VertNfp, UpEidM, BotEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np - Np, sd + ele * Nlayer * Np + L*Np - Np, BotETau, jcs);

					GetLocalUpFacialContributionForMixedSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, BotEFToE, BotEFToF, Nface, LocalEle, \
						LocalStartPoint, jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], Np, VertNfp, BotENe, M3d, M2d, J + ele * Nlayer * Np + L*Np, \
						J2d + ele*VertNfp, UpEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np, sd + ele * Nlayer * Np + L*Np, BotETau, jcs);

				}
				else{
					GetLocalToUpElementContributionForMixedSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, BotEFToE, BotEFToF, Nface, \
						LocalEle, Np, VertNfp, BotENe, jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + ele * Nlayer * Np + L*Np - Np, \
						J2d + ele*VertNfp, UpEidM, BotEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np - Np, sd + ele * Nlayer * Np + L*Np - Np, BotETau, jcs);

					GetLocalUpFacialContributionForMixedSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, BotEFToE, BotEFToF, Nface, LocalEle, \
						LocalStartPoint, jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], Np, VertNfp, BotENe, M3d, M2d, J + ele * Nlayer * Np + L*Np, \
						J2d + ele*VertNfp, UpEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np, sd + ele * Nlayer * Np + L*Np, BotETau, jcs);

					GetLocalToDownElementContributionForMixedSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, BotEFToE, BotEFToF, Nface, \
						LocalEle, Np, VertNfp, BotENe, jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + ele * Nlayer * Np + L*Np + Np, \
						J2d + ele*VertNfp, BotEidM, UpEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np + Np, sd + ele * Nlayer * Np + L*Np + Np, BotETau, jcs);

					GetLocalDownFacialContributionForMixedSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, BotEFToE, BotEFToF, Nface, LocalEle, \
						LocalStartPoint, jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], Np, VertNfp, BotENe, M3d, M2d, J + ele * Nlayer * Np + L*Np, \
						J2d + ele*VertNfp, BotEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np, sd + ele * Nlayer * Np + L*Np, BotETau, jcs);
				}
			}
			free(TempEToE);
			free(FacialVector);
			free(GlobalFace);
			free(AdjEle);
			free(ReverseFlag);
		}
	}
	free(TempIr);
	free(TempJc);
	free(UpEidM);
	free(BotEidM);
}