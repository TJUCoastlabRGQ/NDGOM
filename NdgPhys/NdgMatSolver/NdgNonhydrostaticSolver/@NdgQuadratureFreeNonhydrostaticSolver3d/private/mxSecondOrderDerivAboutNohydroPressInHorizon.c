#include "SWENonhydrostatic3d.h"

/*This function is called at the initialization stage to calcualte the second order derivative
about nonhydrostatic pressure in horizontal direction, i.e. $\frac{\partial^2 p}{\partial x^2}$
and $\frac{\partial^2 p}{\partial y^2}$.
For this term, the primal form is given as:
$$-\int_{\Omega}\nabla_h s\cdot \nabla_h p_hd\boldsymbol{x}+\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}-
\int_{\partial \Omega^D}(p_D-p_h)\nabla_hs\cdot \boldsymbol{n}d\boldsymbol{x} = -\int_{\epsilon_i}(\left\{\nabla_hp_h\right\}-\tau^k[p_h])
[s]d\boldsymbol{x}-\int_{\partial\Omega^D}s\boldsymbol{n}\cdot\nabla_hp_hd\boldsymbol{x}+\int_{\partial \Omega^D}
\tau^ksp_hd\boldsymbol{x} -\int_{\partial \Omega^D}\tau^ksp_Dd\boldsymbol{x} -
\int_{\partial \Omega^N}sp_Nd\boldsymbol{x}$$
We note that at present we don't consider the boundary condition. Boundary condition about this part is imposed when we 
assemble the global stiff matrix.
Detail about this form can be found through any markdown editor.
*/



/*This following function is used to assemble term corresponds to
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
Dr: A pointer to the start of the differential matrix for the master cell in direction r
Ds: A pointer to the start of the differential matrix for the master cell in direction s
Dt: A pointer to the start of the differential matrix for the master cell in direction t
rd: A pointer to the start of the transformation Jacobian of the studied element in direction r, here d can be x or y
sd: A pointer to the start of the transformation Jacobian of the studied element in direction s, here d can be x or y
J: A pointer to the start of the Jacobian of the studied element
*/

void GetLocalVolumuIntegralTermForSecondOrderTerm(double *dest, int StartPoint, int Np, int NonzeroPerColumn,\
	double *mass3d, double *Dr, double *Ds, double *rd, double *sd, double *J)
{
	double *LocalMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalMass3d, mass3d, J, Np);
	double *LumpedMass3d = malloc(Np*Np*sizeof(double));
	MassLumping(LumpedMass3d, LocalMass3d, Np);

	double *DiffMatrixBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrixBuff, Dr, rd, Np);
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	DiagMultiply(DiffMatrix, Ds, sd, Np);
	Add(DiffMatrix, DiffMatrix, DiffMatrixBuff, Np*Np);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	double *TempContributionBuff = malloc(Np*Np*sizeof(double));

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, LumpedMass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, TempContribution, (ptrdiff_t)Np);
	MultiplyByConstant(TempContribution, TempContribution, -1.0, Np*Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempContribution, NonzeroPerColumn, Np);

	free(LocalMass3d);
	free(LumpedMass3d);
	free(DiffMatrixBuff);
	free(DiffMatrix);
	free(TempContribution);
	free(TempContributionBuff);
}


/*The following function is used to assemble terms corresponding to
$$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x} +
\int_{\epsilon_i}(\left\{\nabla_hp_h\right\}-\tau^k[p_h])[s]d\boldsymbol{x}$$
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForSecondOrderTerm.
The rest are as follows:
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
UpEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
Uptz: A pointer to the start of the transformation Jacobian of the element located upside of the studied cell
Localtz: A pointer to the start of the transformation Jacobian of the element of the studied cell
Tau: A pointer to the penalty parameter of the up face
*/

//void GetLocalToAdjacentElementContributionForSecondOrderTerm(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEidM, double *AdjEidM, double *Dr, double *Ds, \
	double *Localrd, double *Adjrd, double *Localsd, double *Adjsd, double *Tau, double *Vector){
void GetLocalToAdjacentElementContributionForSecondOrderTerm(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *mass2d, double *J2d, double *LocalEidM, double *AdjEidM, double *Dr, double *Ds, \
	double *Localrd, double *Adjrd, double *Localsd, double *Adjsd, double *Vector){

	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, mass2d, J2d, Nfp);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *AdjDiffBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(AdjDiffBuff, Dr, Adjrd, Np);
	double *AdjDiff = malloc(Np*Np*sizeof(double));
	DiagMultiply(AdjDiff, Ds, Adjsd, Np);
	Add(AdjDiff, AdjDiff, AdjDiffBuff, Np*Np);

	double *LocalDiffBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiffBuff, Dr, Localrd, Np);
	double *LocalDiff = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiff, Ds, Localsd, Np);
	Add(LocalDiff, LocalDiff, LocalDiffBuff, Np*Np);

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	/*The follow part is for term $\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$*/
	double *InnerEdgeContribution = malloc(Np*Nfp*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, AdjDiff, AdjEidM, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*Vector[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*The follow part is for term $\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEidM, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5*Vector[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEidM, Np, Nfp);

	/*The follow part is for term $ - \int_{ \epsilon_i }\tau^k[p_h][s]d\boldsymbol{ x }$*/
/*	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(TempMass2d, FacialMass2d, Tau, Nfp);
	DiagMultiply(TempMass2d, TempMass2d, Vector, Nfp);
	DiagMultiply(TempMass2d, TempMass2d, Vector, Nfp);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, AdjEidM, LocalEidM, Np, Nfp, 1.0);
	*/
	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempContribution, NonzeroPerColumn, Np);

	free(FacialMass2d);
	free(TempContribution);
	free(AdjDiffBuff);
	free(AdjDiff);
	free(LocalDiffBuff);
	free(LocalDiff);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
//	free(TempMass2d);
}

/*The following function is used to assemble terms corresponding to
$$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x} +
\int_{\epsilon_i}(\left\{\nabla_hp_h\right\}-\tau^k[p_h])[s]d\boldsymbol{x}$$
This part corresponds to the effect of local boundary integral have on the local element.
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForSecondOrderTerm.
The rest are as follows:
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
UpEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
Uptz: A pointer to the start of the transformation Jacobian of the element located upside of the studied cell
Localtz: A pointer to the start of the transformation Jacobian of the element of the studied cell
Tau: A pointer to the penalty parameter of the up face
*/

//void GetLocalFacialContributionForSecondOrderTerm(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEidM, double *Dr, double *Ds, \
	double *Localrd, double *Localsd, double *Tau, double *Vector){
void GetLocalFacialContributionForSecondOrderTerm(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *mass2d, double *J2d, double *LocalEidM, double *Dr, double *Ds, \
	double *Localrd, double *Localsd, double *Vector){

	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, mass2d, J2d, Nfp);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *LocalDiffBuff = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiffBuff, Dr, Localrd, Np);
	double *LocalDiff = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiff, Ds, Localsd, Np);
	Add(LocalDiff, LocalDiff, LocalDiffBuff, Np*Np);

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	/*The follow part is for term $\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$*/
	double *InnerEdgeContribution = malloc(Np*Nfp*sizeof(double));
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEidM, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*Vector[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*The follow part is for term $\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$*/
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEidM, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d,\
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	//DiagMultiply(InnerEdgeContribution, InnerEdgeContribution, Vector, Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*Vector[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*The follow part is for term $ - \int_{ \epsilon_i }\tau^k[p_h][s]d\boldsymbol{ x }$*/
/*	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(TempMass2d, FacialMass2d, Tau, Nfp);
	DiagMultiply(TempMass2d, TempMass2d, Vector, Nfp);
	DiagMultiply(TempMass2d, TempMass2d, Vector, Nfp);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEidM, LocalEidM, Np, Nfp, -1.0);
	*/
	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempContribution, NonzeroPerColumn, Np);

	free(FacialMass2d);
	free(TempContribution);
	free(LocalDiffBuff);
	free(LocalDiff);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
//	free(TempMass2d);
}

void GetLocalToAdjacentContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int AdjEle, int Np, int Np2d, \
	int NonzeroPerColumn, double *mass2d, \
	double *J2d, double *AdjEid, double *LocalEid, double *FToE, int Ne, double *Tau){

	int GlobalFace = 0;
	FindGlobalBottomEdgeFace(&GlobalFace, FToE, LocalEle, AdjEle, Ne);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	DiagMultiply(Mass2d, Mass2d, Tau + Np2d * GlobalFace, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, AdjEid, LocalEid, Np, Np2d, 1.0);

	for (int i = 0; i < EleNumber; i++){
		if (AdjEle == (int)TempEToE[i]){
			int StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, TempContribution, NonzeroPerColumn, Np);
			break;
		}
	}
	free(Mass2d);
	free(TempContribution);
}

void GetLocalPenaltyContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int AdjEle, int Np, int Np2d, \
	int NonzeroPerColumn, double *mass2d, \
	double *J2d, double *Eid, double* FToE, int Ne, double *Tau){

	int GlobalFace = 0;
	FindGlobalBottomEdgeFace(&GlobalFace, FToE, LocalEle, AdjEle, Ne);
	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);

	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	DiagMultiply(Mass2d, Mass2d, Tau + Np2d * GlobalFace, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, Eid, Eid, Np, Np2d, -1.0);

	for (int i = 0; i < EleNumber; i++){
		if (LocalEle == (int)TempEToE[i]){
			int StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, TempContribution, NonzeroPerColumn, Np);
			break;
		}
	}
	free(Mass2d);
	free(TempContribution);
}

void ImposeDirichletBoundaryCondition(double *dest, int StartPoint, int Np, int Np2d, \
	int NonzeroPerColumn, double *mass2d, double *J2d, double *Eid, double *Tau){

	double *Mass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(Mass2d, mass2d, J2d, Np2d);
	/*Since nx = ny = 0 on surface, we only consider the penalty part here*/
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	DiagMultiply(Mass2d, Mass2d, Tau, Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, Eid, Eid, Np, Np2d, -1.0);


	AssembleContributionIntoSparseMatrix(dest + StartPoint, TempContribution, NonzeroPerColumn, Np);

	free(Mass2d);
	free(TempContribution);
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
LAV: Volume of the three-dimensional mesh object, indexed as 8.
J: Jacobian coefficient of the three-dimensional mesh object, indexed as 9.
Js: Jacobian coefficient of the inner edge, indexed as 10.
M3d: Mass matrix of the 3d master cell, indexed as 11.
M: Mass matrix of the 2d master cell of the inner edge, indexed as 12.
BENe: Number of boundary edge, indexed as 13.
IENe: Number of inner edge, indexed as 14.
IELAV: Area of the inner edge object, indexed as 15.
Fmask: Index of the interpolation point on each face of the master cell, indexed as 16.
Nfp: Number of interpolation point on each face of the master cell, indexed as 17.
P: Approximation order of the master cell, indexed as 18.
vector: The direction vector of the inner edge object, indexed as 19.
rd: Differential coefficient of the 3d mesh object in direction r, indexed as 20.
sd: Differential coefficient of the 3d mesh object in direction s, indexed as 21.
Dr: Differential matrix of the 3d master cell in direction r, indexed as 22.
Ds: Differential matrix of the 3d master cell in direction s, indexed as 23.
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
	double *Js = mxGetPr(prhs[11]);
	double *M = mxGetPr(prhs[12]);
	int BENe = mxGetScalar(prhs[13]);
	int IENe = mxGetScalar(prhs[14]);
	double *IELAV = mxGetPr(prhs[15]);

	int TotalNonzero = Ele3d * (Nface + 1) * Np*Np - BENe * Np*Np - 2 * Ele2d*Np*Np;
	mwIndex *TempIr = malloc(TotalNonzero*sizeof(mwIndex));
	mwIndex *TempJc = malloc((Np*Ele3d + 1)*sizeof(mwIndex));
	TempJc[0] = 0;

	double *Fmask = mxGetPr(prhs[16]);
	int maxNfp = mxGetM(prhs[16]);
	double *Nfp = mxGetPr(prhs[17]);
	int IENfp = (int)Nfp[0];
	int VertNfp = (int)Nfp[Nface - 1];
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
	double *Dr = mxGetPr(prhs[22]);
	double *Ds = mxGetPr(prhs[23]);

	double *BotELAV = mxGetPr(prhs[24]);
	int BotENe = (int)mxGetScalar(prhs[25]);
	double *BotEFToE = mxGetPr(prhs[26]);

	double *SurfELAV = mxGetPr(prhs[27]);
	int SurfENe = (int)mxGetScalar(prhs[28]);
	double *SurfEFToE = mxGetPr(prhs[29]);
	double *M2d = mxGetPr(prhs[30]);
	double *J2d = mxGetPr(prhs[31]);

	char* BoundaryType;
	BoundaryType = mxArrayToString(prhs[32]);

	double *J = mxGetPr(prhs[33]);
	double *M3d = mxGetPr(prhs[34]);

	GetSparsePatternInHorizontalDirection(TempIr, TempJc, EToE, Nface, Nface, Ele3d, Np);
	plhs[0] = mxCreateSparse(Np*Ele3d, Np*Ele3d, TotalNonzero, mxREAL);
	double *sr = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);
	memcpy(irs, TempIr, TotalNonzero*sizeof(mwIndex));
	mwIndex *jcs = mxGetJc(plhs[0]);
	memcpy(jcs, TempJc, (Np*Ele3d + 1)*sizeof(mwIndex));

	free(TempIr);
	free(TempJc);
	/*
	
	double *Tau = malloc(IENe*IENfp*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		double localRatio = IELAV[face] / LAV[(int)FToE[face * 2] - 1];
		double adjacentRatio = IELAV[face] / LAV[(int)FToE[face * 2 + 1] - 1];
		for (int p = 0; p < IENfp; p++){
			//Tau[face*IENfp + p] = 2000000;
			Tau[face*IENfp + p] = 100 * max(localRatio*(P + 1)*(P + 3) / 3 * Nface / 2, \
				adjacentRatio*(P + 1)*(P + 3) / 3 * Nface / 2);
//			Tau[face*IENfp + p] = 1.0 / sqrt(IELAV[face]);
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
			BotETau[face*VertNfp + p] = 100 * 0.5 * max(localRatio*(P + 1)*(P + 3) / 3 * Nface / 2, \
				adjacentRatio*(P + 1)*(P + 3) / 3 * Nface / 2);
//			BotETau[face*VertNfp + p] = 0.5 * 1.0 / sqrt(BotELAV[face]);
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
			SurfETau[face*VertNfp + p] = 100 * 0.5 * max(localRatio*(P + 1)*(P + 3) / 3 * Nface / 2, \
				adjacentRatio*(P + 1)*(P + 3) / 3 * Nface / 2);
//			SurfETau[face*VertNfp + p] = 1.0 / sqrt(SurfELAV[face]); //This parameter is doubled
		}
	}
	*/


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int ele = 0; ele < Ele3d; ele++){
		int EleNumber = 0;
		int StartPoint;
		int LocalStartPoint;
		double *TempEToE = malloc((Nface+1)*sizeof(double));
		FindUniqueElementAndSortOrder(TempEToE, EToE + ele*Nface, &EleNumber, Nface, ele + 1);

		for (int i = 0; i < EleNumber; i++){
			if (ele + 1 == (int)TempEToE[i]){
				LocalStartPoint = jcs[ele*Np] + i*Np;
				GetLocalVolumuIntegralTermForSecondOrderTerm(sr, LocalStartPoint, \
					Np, jcs[ele*Np + 1] - jcs[ele*Np], \
					M3d, Dr, Ds, rd + ele*Np, sd + ele*Np, J + ele*Np);
				break;
			}
		}

		double *FacialVector = malloc(Nface2d*IENfp*sizeof(double));
		int *GlobalFace = malloc(Nface2d*sizeof(int));
		int *AdjEle = malloc(Nface2d*sizeof(int));
		int *ReverseFlag = malloc(Nface2d*sizeof(int));
		int InternalFace = 0;

		FindFaceAndDirectionVector(FacialVector, GlobalFace, AdjEle, \
			&InternalFace, ReverseFlag, IENfp, ele + 1, FToE, FToF, Vector, IENe, Nface2d);

		double *LocalEidM = malloc(IENfp*sizeof(double));

		double *AdjEidM = malloc(IENfp*sizeof(double));

		for (int i = 0; i < InternalFace; i++){
			if (ReverseFlag[i] == 0){
				for (int p = 0; p < IENfp; p++){
					LocalEidM[p] = FToN1[GlobalFace[i] * IENfp + p];
					AdjEidM[p] = FToN2[GlobalFace[i] * IENfp + p];
				}
			}
			else{
				for (int p = 0; p < IENfp; p++){
					LocalEidM[p] = FToN2[GlobalFace[i] * IENfp + p];
					AdjEidM[p] = FToN1[GlobalFace[i] * IENfp + p];
				}
			}

			GetLocalFacialContributionForSecondOrderTerm(sr, LocalStartPoint, Np, IENfp, jcs[ele*Np + 1] - jcs[ele*Np], \
				 M, Js + GlobalFace[i] * IENfp, LocalEidM, Dr, Ds, \
				rd + ele*Np, sd + ele*Np, FacialVector + i * IENfp);

			for (int j = 0; j < EleNumber; j++){
				if ((int)TempEToE[j] == AdjEle[i]){
					StartPoint = jcs[ele*Np] + j*Np;
					GetLocalToAdjacentElementContributionForSecondOrderTerm(sr, StartPoint, Np, IENfp, jcs[ele*Np + 1] - jcs[ele*Np], \
					    M, Js + GlobalFace[i] * IENfp, LocalEidM, AdjEidM, Dr, Ds, \
						rd + ele*Np, rd + (int)(TempEToE[j] - 1)*Np, sd + ele*Np, sd + (int)(TempEToE[j] - 1)*Np, \
						FacialVector + i*IENfp);
				}
			}
		}
		free(TempEToE);
		free(FacialVector);
		free(GlobalFace);
		free(AdjEle);
		free(ReverseFlag);
		free(LocalEidM);
		free(AdjEidM);
	}

//	free(Tau);
	/*
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

			if (UpEle == DownEle){//Only one layer in the vertical direction, impose the Dirichlet boundary condition
				for (int i = 0; i < EleNumber; i++){
					if (LocalEle == (int)TempEToE[i]){
						StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
						if (!strcmp(BoundaryType, "Dirichlet")){
							ImposeDirichletBoundaryCondition(sr, StartPoint, Np, VertNfp, \
								jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
								J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, SurfETau + ele*VertNfp);
						}
						break;
					}
				}
			}
			else{//two or more layers included in vertical direction
				if (LocalEle == UpEle){//This is the top most cell, and L=0
					for (int i = 0; i < EleNumber; i++){
						if (LocalEle == (int)TempEToE[i]){
							StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
							if (!strcmp(BoundaryType, "Dirichlet")){
								ImposeDirichletBoundaryCondition(sr, StartPoint, Np, VertNfp, \
									jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
									J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, SurfETau + ele*VertNfp);
							}
							break;
						}
					}

					GetLocalPenaltyContribution(sr, jcs, TempEToE, EleNumber, LocalEle, DownEle, Np, VertNfp,\
						 jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						 J + (LocalEle - 1)*Np, J2d + ele*VertNfp, BotEidM, BotEFToE, BotENe, BotETau);

					GetLocalToAdjacentContribution(sr, jcs, TempEToE, EleNumber, LocalEle, DownEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, BotEidM, BotEFToE, BotENe, BotETau);
				}
				else if (LocalEle == DownEle){// This is the bottom most cell.

					GetLocalPenaltyContribution(sr, jcs, TempEToE, EleNumber, LocalEle, UpEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, BotEFToE, BotENe, BotETau);

					GetLocalToAdjacentContribution(sr, jcs, TempEToE, EleNumber, LocalEle, UpEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, BotEidM, UpEidM, BotEFToE, BotENe, BotETau);
				}
				else{
					// Bottom surface
					GetLocalPenaltyContribution(sr, jcs, TempEToE, EleNumber, LocalEle, DownEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, BotEidM, BotEFToE, BotENe, BotETau);

					GetLocalToAdjacentContribution(sr, jcs, TempEToE, EleNumber, LocalEle, DownEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, BotEidM, BotEFToE, BotENe, BotETau);

					//Upper surface
					GetLocalPenaltyContribution(sr, jcs, TempEToE, EleNumber, LocalEle, UpEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, BotEFToE, BotENe, BotETau);

					GetLocalToAdjacentContribution(sr, jcs, TempEToE, EleNumber, LocalEle, UpEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, BotEidM, UpEidM, BotEFToE, BotENe, BotETau);

				}

			}
			free(TempEToE);
		}
	}
	*/

	free(UpEidM);
	free(BotEidM);
//	free(BotETau);
//	free(SurfETau);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < TotalNonzero; i++){
		sr[i] = sr[i] + pow(10.0,-16);
	}
}
