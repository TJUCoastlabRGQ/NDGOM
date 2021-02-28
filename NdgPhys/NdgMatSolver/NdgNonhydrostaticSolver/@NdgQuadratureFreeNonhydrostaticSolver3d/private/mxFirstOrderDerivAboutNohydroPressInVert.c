#include "SWENonhydrostatic3d.h"
/*This function is called at the initialization stage to calcualte the first order derivative 
 about nonhydrostatic pressure in vertical direction, i.e. $\frac{\partial p}{\partial \sigma}$.
 For this term, the primal form is given as:  
 $\int_{\Omega}\boldsymbol{\sigma}_h\cdot\boldsymbol{v}d\boldsymbol{x} = 
\int_{\Omega}\boldsymbol{v}\cdot \nabla_h p_h d\boldsymbol{x}  - 
\int_{{\epsilon}_i} [p_h] \cdot \boldsymbol{\left\{v\right\}}d\boldsymbol{x} +
\int_{\partial \Omega^D} p_D \boldsymbol{v}\cdot\boldsymbol{n}d\boldsymbol{x}-
\int_{\partial \Omega^D}p_h\boldsymbol{v}\cdot\boldsymbol{n}d\boldsymbol{x}$
Detail about this form can be found through any markdown editor.
*/


/*This following function is used to assemble term corresponds to 
$\int_{\Omega}\boldsymbol{v}\cdot \nabla_h p_h d\boldsymbol{x}$.
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
void GetLocalVolumuIntegralTermForFirstOrderTerm(double *dest, int Startpoint, \
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
	double *Contribution = malloc(Np*Np*sizeof(double));
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, Mass3d,
		(ptrdiff_t)Np, DiffSigma, (ptrdiff_t)Np, 0.0, TempContribution, (ptrdiff_t)Np);
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
	free(Contribution);
}


/*The following funciton is used to assemble terms corresponding to 
$-\int_{\partial \Omega^D}p_h\boldsymbol{v}\cdot\boldsymbol{n}d\boldsymbol{x}$.
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForFirstOrderTerm.
The rest are as follows:
Np2d: Number of interpolation point for the 2d element.
      For the current case, the cell on top of the 3d element.
mass2d: A pointer to the start of the mass matrix for the 2 dimensional master cell
J2d: A pointer to the start of the Jacobian of the 2 dimensional studied element
Eid: A pointer to the start of the interpolation point index on the top of the master cell
*/
void ImposeFirstOrderNonhydroDirichletBoundaryCondition(double *dest, int StartPoint, int Np, \
	int Np2d, int NonzeroPerColumn, double *mass3d, double *mass2d, double *J, \
	double *J2d, double *Eid){
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
	AssembleContributionIntoRowAndColumn(TempContribution, Mass2d, Eid, Eid, Np, Np2d, -1);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);
	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
	free(Mass3d);
	free(InvMass3d);
	free(Mass2d);
	free(TempContribution);
	free(Contribution);
}

/*The following function is used to assemble terms corresponding to
$- \int_{{\epsilon}_i} [p_h] \cdot \boldsymbol{\left\{v\right\}}d\boldsymbol{x}$
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForFirstOrderTerm
and ImposeNonhydroDirichletBoundaryCondition.
The rest are as follows:
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
UpEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
*/
void GetLocalToUpElementContributionForFirstOrderTerm(double *dest, int StartPoint, int Np, int Np2d, int NonzeroPerColumn, \
	double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *UpEid){
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
	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, Mass2d, -0.5, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, UpEid, LocalEid, Np, Np2d, 1);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvUpMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);
	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
	free(UpMass3d);
	free(InvUpMass3d);
	free(Mass2d);
	free(TempContribution);
	free(Contribution);
	free(TempMass2d);
}

/*The following funciton is used to assemble terms corresponding to
$- \int_{{\epsilon}_i} [p_h] \cdot \boldsymbol{\left\{v\right\}}d\boldsymbol{x}$
Explanation of most variables can be found in function GetLocalVolumuIntegralTermForFirstOrderTerm
and ImposeNonhydroDirichletBoundaryCondition.
The rest are as follows:
LocalEid: A pointer to the start of the interpolation point index on the local face of the master cell
DownEid: A pointer to the start of the interpolation point index on the adjacent face of the master cell
*/
void GetLocalToDownElementContributionForFirstOrderTerm(double *dest, int StartPoint, int Np, int Np2d, int NonzeroPerColumn, \
	double *mass3d, double *mass2d, double *J, double *J2d, double *LocalEid, double *DownEid){
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
	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, Mass2d, 0.5, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, DownEid, LocalEid, Np, Np2d, 1);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvDownMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);
	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
	free(DownMass3d);
	free(InvDownMass3d);
	free(Mass2d);
	free(TempContribution);
	free(Contribution);
	free(TempMass2d);
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

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < Ele2d; e++){
		int StartPoint;

		for (int L = 0; L < Nlayer; L++){
			//Index of the studied element
			int LocalEle = e*Nlayer + L + 1;
			//Index of the element that located upside of the studied element 
			int UpEle = (int)EToE[e*Nlayer*Nface + L*Nface + Nface - 1];
			//Index of the element that located downside of the studied element
			int DownEle = (int)EToE[e*Nlayer*Nface + L*Nface + Nface - 2];

			GetLocalVolumuIntegralTermForFirstOrderTerm(sr, jcs[e*Nlayer*Np + L*Np], \
				Np, jcs[e*Nlayer*Np + L*Np + 1] - jcs[e*Nlayer*Np + L*Np], M3d, \
				Dt, tz + e*Nlayer*Np + L*Np, J + e*Nlayer*Np + L*Np, L);

			if (UpEle == DownEle){ // we have only one layer in vertical direction, since L=0, we omit them in startpoint
				StartPoint = (int)jcs[e * 1 * Np];
				ImposeFirstOrderNonhydroDirichletBoundaryCondition(sr, StartPoint, Np, \
					Np2d, jcs[e * 1 * Np + 1] - jcs[e * 1 * Np], M3d, M2d, J + e * 1 * Np, \
					J2d + e*Np2d, UpEidM);
			}
			else{//two or more layers included in vertical direction
				if (LocalEle == UpEle){//This is the top most cell, and L=0
					StartPoint = (int)jcs[e*Nlayer*Np + L*Np] + 0;

					ImposeFirstOrderNonhydroDirichletBoundaryCondition(sr, StartPoint, Np, \
						Np2d, jcs[e * Nlayer * Np + 1] - jcs[e * Nlayer * Np], M3d, M2d, J + e * Nlayer * Np, \
						J2d + e*Np2d, UpEidM);

					GetLocalToDownElementContributionForFirstOrderTerm(sr, StartPoint + Np, Np, Np2d, \
						jcs[e * Nlayer * Np + 1] - jcs[e * Nlayer * Np], \
						M3d, M2d, J + e * Nlayer * Np + Np, J2d + e*Np2d, BotEidM, UpEidM);
				}
				else if (LocalEle == DownEle){// This is the bottom most cell.
					/*For this situation, this is the bottom most cell, since Newmann boundary
					condition is added, for the first order derivative about non-hydrostatic
					pressure, we do nothing here*/
					StartPoint = (int)jcs[e*Nlayer*Np + L*Np + 1] - 2 * Np;

					GetLocalToUpElementContributionForFirstOrderTerm(sr, StartPoint, Np, Np2d, \
						jcs[e*Nlayer*Np + L*Np + 1] - jcs[e*Nlayer*Np + L*Np], \
						M3d, M2d, J + (e*Nlayer + L - 1)*Np, J2d + e*Np2d, UpEidM, BotEidM);

				}
				else{
					/*
					Element locates in the middle of the studied column, and there are elements located both upside and downside
					*/
					StartPoint = (int)jcs[e*Nlayer*Np + L*Np];

					GetLocalToUpElementContributionForFirstOrderTerm(sr, StartPoint, Np, Np2d, \
						jcs[e*Nlayer*Np + L*Np + 1] - jcs[e*Nlayer*Np + L*Np], \
						M3d, M2d, J + (e*Nlayer + L - 1)*Np, J2d + e*Np2d, UpEidM, BotEidM);

					StartPoint = (int)jcs[e*Nlayer*Np + L*Np] + 2 * Np;

					GetLocalToDownElementContributionForFirstOrderTerm(sr, StartPoint, Np, Np2d, \
						jcs[e * Nlayer * Np + 1] - jcs[e * Nlayer * Np], \
						M3d, M2d, J + e * Nlayer * Np + Np, J2d + e*Np2d, BotEidM, UpEidM);
				}
			}
		}
	}

	free(UpEidM);
	free(BotEidM);
}