#include "SWENonhydrostatic3d.h"

/*This function is called at the initialization stage to calcualte the mixed second order derivative
about nonhydrostatic pressure in horizontal direction, i.e. $\frac{\partial}{\partial x}(\frac{\partial q}{\partial \sigma})$
and $\frac{\partial}{\partial y}(\frac{\partial q}{\partial \sigma})$.
For this term, the primal form is given as:
$$-\int_{\Omega_e}\nabla_hv\nabla_{\sigma}q_hd\boldsymbol{x}+\int_{\epsilon_{iv}}[q_h]\left\{\nabla_hv\right\}
d\boldsymbol{x}=-\int_{\epsilon_{ih}}\hat Q_h[v]d\boldsymbol{x}-\int_{\epsilon_{bh}}\hat Q_hn_hvd\boldsymbol{x}
\\=-\int_{\epsilon_{ih}}\left\{\nabla_hq_h\right\}[v]d\boldsymbol{x} - \int_{\epsilon_{bh}}\hat Q_hn_hvd\boldsymbol{x}$$
At present, the last term in the primal form is neglected, i.e. the numerical flux $\hat Q_hn_h$is set to zero at the lateral boundary.
*/


/*This following funciton is used to assemble term corresponds to
$$-\int_{\Omega}\nabla_h v \nabla_{\sigma} q_hd\boldsymbol{x}$$
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

void GetLocalVolumuIntegralTermForSecondOrderTerm(double *dest, int StartPoint, int Np, int NonzeroPerColumn, \
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
	MultiplyByConstant(TempContribution, TempContribution, -1, Np*Np);
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
$$\int_{\epsilon_{ih}}\left\{\nabla_vq_h\right\}[v]d\boldsymbol{x}$$
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
maxNfp: Number of the maximum facial interpolation point
Nfp: Number of facial interpolation point for lateral face
Ne: Number of edges for inner edge
Np: Number of interpolation point for the master cell
Vector: A pointer to the start of the facial direction vector of the inner edge
Jcs: A pointer to the array that contains the number of nonzero element in each column
*/

void GetLocalToAdjacentElementContributionInHorizontalDirection(double *dest, int LocalEle, int Nface2d, int Nface, \
	double *EToE, double *mass2d, double *mass3d, double *Dt, double *tz, double *Js, double *J, double *TempEToE, \
	double *FToE, double *FToF, double *Fmask, int maxNfp, int Nfp, int Ne, int Np, double *Vector, mwIndex *Jcs){
	
	for (int i = 0; i < Nface2d; i++){
		if (LocalEle == (int)EToE[i]) continue;
		else{
			int GlobalFace = 0, LocalFace = 0, AdjFace = 0;
			double *FacialVector = malloc(Nfp * sizeof(double));
			FindFaceAndDirectionVector(FacialVector, &GlobalFace, &LocalFace, &AdjFace, \
				Nfp, LocalEle, EToE[i], FToE, FToF, Vector, Ne);
			double *LocalEidM = malloc(Nfp*sizeof(double));
			double *AdjEidM = malloc(Nfp*sizeof(double));

			for (int p = 0; p < Nfp; p++){
				LocalEidM[p] = Fmask[(LocalFace - 1)*maxNfp + p];
				AdjEidM[p] = Fmask[(AdjFace - 1)*maxNfp + p];
			}

			double *AdjDiff = malloc(Np*Np*sizeof(double));
			DiagMultiply(AdjDiff, Dt, tz + (int)(EToE[i] - 1)*Np, Np);

			double *AdjMass3d = malloc(Np*Np*sizeof(double));
			DiagMultiply(AdjMass3d, mass3d, J + Np*(int)(EToE[i] - 1), Np);
			double *InvAdjMass3d = malloc(Np*Np*sizeof(double));
			memcpy(InvAdjMass3d, AdjMass3d, Np*Np*sizeof(double));
			MatrixInverse(InvAdjMass3d, (ptrdiff_t)Np);
			double *Mass2d = malloc(Nfp*Nfp*sizeof(double));
			DiagMultiply(Mass2d, mass2d, Js + Nfp*GlobalFace, Nfp);

			double *TempContribution = malloc(Np*Np*sizeof(double));
			memset(TempContribution, 0, Np*Np*sizeof(double));
			double *Contribution = malloc(Np*Np*sizeof(double));

			double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

			double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
			AssembleFacialDiffMatrix(FacialDiffMatrix, AdjDiff, AdjEidM, Nfp, Np);
			DiagMultiply(Mass2d, Mass2d, FacialVector, Nfp);
			MultiplyByConstant(Mass2d, Mass2d, -0.5, Nfp*Nfp);
			MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, Mass2d,
				(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Nfp);

			AssembleContributionIntoRow(TempContribution, EdgeContribution, AdjEidM, Np, Nfp);

			MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvAdjMass3d,
				(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

			for (int j = 0; j < Nface + 1; j++){
				if (EToE[i] == TempEToE[j]){
					int StartPoint = Jcs[(LocalEle-1)*Np] + j*Np;
					int NonzeroPerColumn = Jcs[(LocalEle - 1)*Np + 1] - Jcs[(LocalEle - 1)*Np];
					AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
					break;
				}
			}

			free(FacialVector);
			free(LocalEidM);
			free(AdjEidM);
			free(AdjDiff);
			free(AdjMass3d);
			free(InvAdjMass3d);
			free(Mass2d);
			free(TempContribution);
			free(Contribution);
			free(FacialDiffMatrix);
			free(EdgeContribution);
		}
	}
}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_{bv}}^{top}q_hvn_{\sigma}d\boldsymbol{x}$$
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
	int NonzeroPerColumn, double *mass3d, double *mass2d, double *J, double *J2d, double *Eid){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvMass3d, (ptrdiff_t)Np);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, Eid, Eid, Np, Np2d, 1);
	/*Multiply the contribution by inverse matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(Mass3d);
	free(InvMass3d);
	free(Mass2d);
	free(TempContribution);
	free(Contribution);

}

/*The following part is used to assemble term corresponding to:
$$\int_{\epsilon_{iv}}[q_h]\left\{\nabla_hv\right\}d\boldsymbol{x}$$
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

void GetLocalToDownElementContributionForSecondOrderTerm(double *dest, double *EToE, double *TempEToE, int Nface, int LocalEle, int Np, int Np2d,\
	int NonzeroPerColumn, double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *BotEid, double *Dr, double *Ds, double *rd, \
	double *sd, mwIndex *Jcs){

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

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, BotEid, Np2d, Np);
	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(EdgeContribution, EdgeContribution, -0.5, Np2d*Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, BotEid, Np, Np2d);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvDownMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	for (int i = 0; i < (Nface+1); i++){
		if (EToE[Nface - 2] == TempEToE[i]){
			int StartPoint = Jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}
	free(Mass2d);
	free(Mass3d);
	free(InvDownMass3d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}


void GetLocalToUpElementContributionForSecondOrderTerm(double *dest, double *EToE, double *TempEToE, int Nface, int LocalEle, int Np, int Np2d, \
	int NonzeroPerColumn, double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *UpEid, double *Dr, double *Ds, double *rd, \
	double *sd, mwIndex *Jcs){

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

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Np2d, Np);
	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, Mass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np2d*Np);
	AssembleContributionIntoColumn(TempContribution, EdgeContribution, UpEid, Np, Np2d);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvUpMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	for (int i = 0; i < (Nface + 1); i++){
		if (EToE[Nface - 1] == TempEToE[i]){
			int StartPoint = Jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}
	free(Mass2d);
	free(Mass3d);
	free(InvUpMass3d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(TempContribution);
	free(Contribution);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}


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
//	double *LAV = mxGetPr(prhs[8]);
	double *J = mxGetPr(prhs[8]);
	double *J2d = mxGetPr(prhs[9]);
	double *Js = mxGetPr(prhs[10]);
	double *M3d = mxGetPr(prhs[11]);
	double *LM2d = mxGetPr(prhs[12]);
	double *M2d = mxGetPr(prhs[13]);
	int BENe = mxGetScalar(prhs[14]);
	int IENe = mxGetScalar(prhs[15]);
	int TotalNonzero = Ele3d * (Nface + 1) * Np*Np - BENe * Np*Np - 2 * Ele2d*Np*Np;
	mwIndex *TempIr = malloc(TotalNonzero*sizeof(mwIndex));
	mwIndex *TempJc = malloc((Np*Ele3d + 1)*sizeof(mwIndex));
	TempJc[0] = 0;
	GetSparsePatternInHorizontalDirection(TempIr, TempJc, EToE, Nface, Nface, Ele3d, Np);


	double *Fmask = mxGetPr(prhs[16]);
	int maxNfp = mxGetM(prhs[16]);
	double *Nfp = mxGetPr(prhs[17]);
	int HorNfp = Nfp[0];
	int VertNfp = Nfp[Nface - 1];
	double *UpEidM = malloc(VertNfp*sizeof(double));
	double *BotEidM = malloc(VertNfp*sizeof(double));
	for (int i = 0; i < VertNfp; i++){
		UpEidM[i] = Fmask[(Nface - 1)*maxNfp + i];
		BotEidM[i] = Fmask[(Nface - 2)*maxNfp + i];
	}

	int P = (int)mxGetScalar(prhs[18]);
	double *Vector = mxGetPr(prhs[19]);
	double *rd = mxGetPr(prhs[20]);
	double *sd = mxGetPr(prhs[21]);
	double *tz = mxGetPr(prhs[22]);
	double *Dr = mxGetPr(prhs[23]);
	double *Ds = mxGetPr(prhs[24]);
	double *Dt = mxGetPr(prhs[25]);

	plhs[0] = mxCreateSparse(Np*Ele3d, Np*Ele3d, TotalNonzero, mxREAL);
	double *sr = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);
	memcpy(irs, TempIr, TotalNonzero*sizeof(mwIndex));
	mwIndex *jcs = mxGetJc(plhs[0]);
	memcpy(jcs, TempJc, (Np*Ele3d + 1)*sizeof(mwIndex));


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int ele = 0; ele < Ele2d; ele++){
		int StartPoint;
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
					StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
					GetLocalVolumuIntegralTermForSecondOrderTerm(sr, StartPoint, Np, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], \
						M3d, Dr, Ds, Dt, rd + (LocalEle - 1)*Np, sd + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, J + (LocalEle - 1)*Np);
					break;
				}
			}

			GetLocalToAdjacentElementContributionInHorizontalDirection(sr, LocalEle, Nface2d, Nface, \
				EToE + (LocalEle-1)*Nface, LM2d, M3d, Dt, tz, Js, J, TempEToE, FToE, FToF, Fmask, maxNfp, HorNfp, IENe, Np, Vector, jcs);

			if (UpEle == DownEle){//Only one layer in the vertical direction, impose the Dirichlet boundary condition
				StartPoint = jcs[(LocalEle - 1)*Np];
				ImposeSecondOrderNonhydroDirichletBoundaryCondition(sr, StartPoint, Np, VertNfp, \
					jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM);
			}
			else{//two or more layers included in vertical direction
				if (LocalEle == UpEle){//This is the top most cell, and L=0
					StartPoint = jcs[(LocalEle - 1)*Np + L*Np];

					ImposeSecondOrderNonhydroDirichletBoundaryCondition(sr, StartPoint, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM);

					GetLocalToDownElementContributionForSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, Nface, LocalEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + ele * Nlayer * Np + L*Np + Np, J2d + ele*VertNfp, \
						BotEidM, UpEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np + Np, sd + ele * Nlayer * Np + L*Np + Np, jcs);
				}
				else if (LocalEle == DownEle){// This is the bottom most cell.

					GetLocalToUpElementContributionForSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, Nface, LocalEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + ele * Nlayer * Np + L*Np - Np, J2d + ele*VertNfp, \
						UpEidM, BotEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np - Np, sd + ele * Nlayer * Np + L*Np - Np, jcs);
				}
				else{
					GetLocalToUpElementContributionForSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, Nface, LocalEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + ele * Nlayer * Np + L*Np - Np, J2d + ele*VertNfp, \
						UpEidM, BotEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np - Np, sd + ele * Nlayer * Np + L*Np - Np, jcs);

					GetLocalToDownElementContributionForSecondOrderTerm(sr, EToE + Nface*(LocalEle - 1), TempEToE, Nface, LocalEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, J + ele * Nlayer * Np + L*Np + Np, J2d + ele*VertNfp, \
						BotEidM, UpEidM, Dr, Ds, rd + ele * Nlayer * Np + L*Np + Np, sd + ele * Nlayer * Np + L*Np + Np, jcs);
				}
			}
			free(TempEToE);
		}
	}
	free(UpEidM);
	free(BotEidM);
	free(TempIr);
	free(TempJc);
}