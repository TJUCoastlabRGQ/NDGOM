#include "SWENonhydrostatic3d.h"


void GetLocalVolumuIntegralTermForFirstOrderHorizontalTerm(double *dest, mwIndex *Ir, mwIndex *Jc, \
	int Np, double *Dr, double *Ds, double *rd, double *sd, int LocalEle) {

	double *DiffMatrixBuff = malloc(Np*Np * sizeof(double));
	DiagMultiply(DiffMatrixBuff, Dr, rd, Np);
	double *DiffMatrix = malloc(Np*Np * sizeof(double));
	DiagMultiply(DiffMatrix, Ds, sd, Np);
	Add(DiffMatrix, DiffMatrix, DiffMatrixBuff, Np*Np);

	AssembleVolumnContributionIntoSparseMatrix(dest, Ir, Jc, Np, DiffMatrix, LocalEle);

	free(DiffMatrixBuff);
	free(DiffMatrix);
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial \sigma}\right \}[p]_xd\Omega+ \\
\int_{\epsilon_i}\left \{\frac{\partial v}{\partial x}\right \}[p]_{\sigma}d\Omega + \\
\int_{\epsilon_i}\left \{\frac{\partial p}{\partial \sigma} \right \}[v]_xd\Omega+ \\
\int_{\epsilon_i}\left \{\frac{\partial p}{\partial x} \right \}[v]_{\sigma}d\Omega  - \\
\int_{\epsilon_i}\tau[v]_x[p]_xd\Omega - \int_{\epsilon_i}\tau[v]_{\sigma}[p]_{\sigma}d\Omega$$.
For this primal form, the second, forth and the last terms, i.e. $\int_{\epsilon_i}\left \{\frac{\partial v}{\partial x}\right \}[p]_{\sigma}d\Omega$,
$\int_{\epsilon_i}\left \{\frac{\partial p}{\partial x} \right \}[v]_{\sigma}d\Omega$ and $\int_{\epsilon_i}\tau[v]_{\sigma}[p]_{\sigma}d\Omega$ 
equal to zero since $n_{\sigma}=0$ on the lateral face.
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

void GetLocalToAdjacentElementContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, double *mass3d, \
	double *mass2d, double *J, double *Js, int LocalEle, double *AdjEle, double *FacialVector, double *GlobalFace,\
	double InternalFace, int Nface, double *TempEToE, double *Flag, double *FToN1, double *FToN2){

	double *LocalEidM = malloc(Nfp*sizeof(double));
	double *AdjEidM = malloc(Nfp*sizeof(double));
	for (int i = 0; i < (int)InternalFace; i++){
		if ((int)Flag[i] == 0){//The local and adjacent facial point are the same with the data stored in FToN1 and FToN2
			for (int p = 0; p < Nfp; p++){
				LocalEidM[p] = FToN1[(int)GlobalFace[i] * Nfp + p];
				AdjEidM[p] = FToN2[(int)GlobalFace[i] * Nfp + p];
			}
		}
		else{
			for (int p = 0; p < Nfp; p++){//The local and adjacent facial point are reversed in FToN1 and FToN2
				AdjEidM[p] = FToN1[(int)GlobalFace[i] * Nfp + p];
				LocalEidM[p] = FToN2[(int)GlobalFace[i] * Nfp + p];
			}
		}

		double *AdjMass3d = malloc(Np*Np*sizeof(double));
		DiagMultiply(AdjMass3d, mass3d, J + Np*((int)AdjEle[i] - 1), Np);
		double *InvAdjMass3d = malloc(Np*Np*sizeof(double));
		memcpy(InvAdjMass3d, AdjMass3d, Np*Np*sizeof(double));
		MatrixInverse(InvAdjMass3d, (ptrdiff_t)Np);
		double *Mass2d = malloc(Nfp*Nfp*sizeof(double));
		DiagMultiply(Mass2d, mass2d, Js + Nfp*(int)GlobalFace[i], Nfp);
		double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));

		double *TempContribution = malloc(Np*Np*sizeof(double));
		memset(TempContribution, 0, Np*Np*sizeof(double));
		double *Contribution = malloc(Np*Np*sizeof(double));

		DiagMultiply(EleMass2d, Mass2d, FacialVector + i*Nfp, Nfp);

		MultiplyByConstant(EleMass2d, EleMass2d, 0.5, Nfp*Nfp);
		
		AssembleContributionIntoRowAndColumn(TempContribution, EleMass2d, AdjEidM, LocalEidM, Np, Nfp, -1.0);
		
		MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvAdjMass3d,
			(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

		AssembleFacialContributionForFirstOrderTermIntoSparseMatrix(dest, Ir, Jc, LocalEidM, Np, Nfp, Contribution, LocalEle, (int)AdjEle[i]);

		free(AdjMass3d);
		free(InvAdjMass3d);
		free(Mass2d);
		free(EleMass2d);
		free(TempContribution);
		free(Contribution);
	}
	free(LocalEidM);
	free(AdjEidM);
}


/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_i}\left\{\frac{\partial v}{\partial \sigma}\right \}[p]_xd\Omega+ \\
\int_{\epsilon_i}\left \{\frac{\partial v}{\partial x}\right \}[p]_{\sigma}d\Omega + \\
\int_{\epsilon_i}\left \{\frac{\partial p}{\partial \sigma} \right \}[v]_xd\Omega+ \\
\int_{\epsilon_i}\left \{\frac{\partial p}{\partial x} \right \}[v]_{\sigma}d\Omega  - \\
\int_{\epsilon_i}\tau[v]_x[p]_xd\Omega - \int_{\epsilon_i}\tau[v]_{\sigma}[p]_{\sigma}d\Omega$$.
For this primal form, the second, forth and the last terms, i.e. $\int_{\epsilon_i}\left \{\frac{\partial v}{\partial x}\right \}[p]_{\sigma}d\Omega$,
$\int_{\epsilon_i}\left \{\frac{\partial p}{\partial x} \right \}[v]_{\sigma}d\Omega$ and $\int_{\epsilon_i}\tau[v]_{\sigma}[p]_{\sigma}d\Omega$
equal to zero since $n_{\sigma}=0$ on the lateral face.
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

void GetLocalFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *mass3d, double *mass2d, double *J, double *Js, int LocalEle, double *Vector, double *GlobalFace, \
	double InternalFace, double *Flag, double *FToN1, double *FToN2){

	double *LocalEidM = malloc(Nfp*sizeof(double));
	for (int i = 0; i < (int)InternalFace; i++){
		if (Flag[i] == 0){
			for (int p = 0; p < Nfp; p++){
				LocalEidM[p] = FToN1[(int)GlobalFace[i] * Nfp + p];
			}
		}
		else{
			for (int p = 0; p < Nfp; p++){
				LocalEidM[p] = FToN2[(int)GlobalFace[i] * Nfp + p];
			}
		}

		double *LocalMass3d = malloc(Np*Np*sizeof(double));
		DiagMultiply(LocalMass3d, mass3d, J + Np*(LocalEle - 1), Np);
		double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
		memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
		MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);

		double *Mass2d = malloc(Nfp*Nfp*sizeof(double));
		DiagMultiply(Mass2d, mass2d, Js + Nfp*(int)GlobalFace[i], Nfp);

		double *TempContribution = malloc(Np*Np*sizeof(double));
		memset(TempContribution, 0, Np*Np*sizeof(double));
		double *Contribution = malloc(Np*Np*sizeof(double));

		DiagMultiply(Mass2d, Mass2d, Vector + i*Nfp, Nfp);

		MultiplyByConstant(Mass2d, Mass2d, 0.5, Nfp*Nfp);

		AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, LocalEidM, LocalEidM, Np, Nfp, -1.0);

		MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
			(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

		AssembleFacialContributionForFirstOrderTermIntoSparseMatrix(dest, Ir, Jc, LocalEidM, Np, Nfp, Contribution, LocalEle, LocalEle);

		free(LocalMass3d);
		free(InvLocalMass3d);
		free(Mass2d);
		free(TempContribution);
		free(Contribution);
	}
	free(LocalEidM);
}

void ImposeHorizonFirstOrderNonhydroDirichletBoundaryCondition(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, \
	int Np2d, int Ele, double *mass3d, double *mass2d, double *J, \
	double *J2d, double *Eid, double *Vector) {

	double *Mass3d = malloc(Np*Np * sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvMass3d = malloc(Np*Np * sizeof(double));
	memcpy(InvMass3d, Mass3d, Np*Np * sizeof(double));
	MatrixInverse(InvMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d * sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np * sizeof(double));
	memset(TempContribution, 0, Np*Np * sizeof(double));
	double *Contribution = malloc(Np*Np * sizeof(double));

	MultiplyByConstant(Mass2d, Mass2d, Vector[0], Np2d*Np2d);

	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, Eid, Eid, Np, Np2d, -1.0);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleFacialContributionForFirstOrderTermIntoSparseMatrix(dest, Ir, Jc, Eid, Np, Np2d, Contribution, Ele, Ele);

	free(Mass3d);

	free(InvMass3d);

	free(Mass2d);

	free(TempContribution);

	free(Contribution);
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
	double *J = mxGetPr(prhs[10]);
	double *Js = mxGetPr(prhs[11]);
	double *M3d = mxGetPr(prhs[12]);
	double *LM2d = mxGetPr(prhs[13]);
	int IENe = mxGetScalar(prhs[14]);

	double *Nfp = mxGetPr(prhs[15]);
	int HorNfp = (int)Nfp[0];

	double *Vector = mxGetPr(prhs[16]);
	double *rd = mxGetPr(prhs[17]);
	double *sd = mxGetPr(prhs[18]);
	double *Dr = mxGetPr(prhs[19]);
	double *Ds = mxGetPr(prhs[20]);
	double *Fmask = mxGetPr(prhs[21]);
	int MaxNfp = (int)mxGetM(prhs[21]);

	const mxArray *BoundaryEdge = prhs[22];
	mxArray *TempBEMb = mxGetField(BoundaryEdge, 0, "M");
	double *BEMb = mxGetPr(TempBEMb);
	mxArray *TempBEJs = mxGetField(BoundaryEdge, 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	mxArray *TempBEFToE = mxGetField(BoundaryEdge, 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToN1 = mxGetField(BoundaryEdge, 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);
	mxArray *TempBENfp = mxGetField(BoundaryEdge, 0, "Nfp");
	int BENfp = (int)mxGetScalar(TempBENfp);
	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);

	signed char *ftype = (signed char *)mxGetData(prhs[23]);
	double *BEVector = mxGetPr(prhs[24]);

	int TotalNonzero = Ele3d * Np * Np + 2 * IENe * HorNfp * Np;
	mwIndex *TempIr = malloc(TotalNonzero*sizeof(mwIndex));
	mwIndex *TempJc = malloc((Np*Ele3d + 1)*sizeof(mwIndex));
	TempJc[0] = 0;

	GetSparsePatternForHorizontalFirstOrderTerm(TempIr, TempJc, EToE, FToE, FToN1, FToN2, \
		Nface, HorNfp, MaxNfp, Np, Ele3d, IENe, Fmask);
	
	plhs[0] = mxCreateSparse(Np*Ele3d, Np*Ele3d, TotalNonzero, mxREAL);
	double *sr = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);
	memcpy(irs, TempIr, TotalNonzero*sizeof(mwIndex));
	mwIndex *jcs = mxGetJc(plhs[0]);
	memcpy(jcs, TempJc, (Np*Ele3d + 1)*sizeof(mwIndex));

	int ele, L, face, i;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (ele = 0; ele < Ele2d; ele++){
		for ( L = 0; L < Nlayer; L++){
			//Index of the studied element
			int LocalEle = ele*Nlayer + L + 1;

			double EleNumber = 0;

			double *TempEToE = malloc((Nface + 1)*sizeof(double));

			FindUniqueElementAndSortOrder(TempEToE, EToE + (LocalEle - 1)*Nface, &EleNumber, Nface2d, LocalEle);

			for (i = 0; i < EleNumber; i++){
				if (LocalEle == TempEToE[i]){

					GetLocalVolumuIntegralTermForFirstOrderHorizontalTerm(sr, irs, jcs, \
						Np, Dr, Ds, rd + (LocalEle - 1)*Np, sd + (LocalEle - 1)*Np, LocalEle);

					break;
				}
			}

			double *FacialVector = malloc(Nface2d*HorNfp*sizeof(double));
			double *GlobalFace = malloc(Nface2d*sizeof(double));
			double *AdjEle = malloc(Nface2d*sizeof(double));
			double *ReverseFlag = malloc(Nface2d*sizeof(double));
			double InternalFace = 0;
			
			FindFaceAndDirectionVector(FacialVector, GlobalFace, AdjEle, \
				&InternalFace, ReverseFlag, HorNfp, LocalEle, FToE, FToF, Vector, IENe, Nface2d);

			GetLocalToAdjacentElementContributionInHorizontalDirection(sr, irs, jcs, Np, HorNfp, M3d, LM2d, J, Js, LocalEle, AdjEle, \
				FacialVector, GlobalFace, InternalFace, Nface, TempEToE, ReverseFlag, FToN1, FToN2);

			GetLocalFacialContributionInHorizontalDirection(sr, irs, jcs, Np, HorNfp, M3d, LM2d, J, Js, LocalEle, \
				FacialVector, GlobalFace, InternalFace, ReverseFlag, FToN1, FToN2);

			free(TempEToE);
			free(FacialVector);
			free(GlobalFace);
			free(AdjEle);
			free(ReverseFlag);
		}
	}
	free(TempIr);
	free(TempJc);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (face = 0; face < BENe; face++) {
		NdgEdgeType type = (NdgEdgeType)ftype[face];
		//if ((type == NdgEdgeClampedDepth) || (type == NdgEdgeClampedVel) || (type == NdgEdgeZeroGrad)) {
		//if ((type == NdgEdgeClampedVel) || (type == NdgEdgeClampedDepth)) {
		if ((type == NdgEdgeClampedVel)) {
			ImposeHorizonFirstOrderNonhydroDirichletBoundaryCondition(sr, irs, jcs, Np, \
				BENfp, (int)BEFToE[2 * face], M3d, BEMb, J + ((int)BEFToE[2 * face] - 1)*Np, \
				BEJs + face*BENfp, BEFToN1 + face*BENfp, BEVector + face*BENfp);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (i = 0; i < TotalNonzero; i++) {
		sr[i] += 1.0e-16;
	}
}