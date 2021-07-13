#include "NdgMath.h"
//#include <omp.h>

void Add(double *dest, double *sourcea, double *sourceb, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] + sourceb[i];
}

void AddByConstant(double *dest, double *sourcea, double ConstData, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] + ConstData;
}

/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the column index, and has been checked*/
void AssembleContributionIntoColumn(double *dest, double *source, double *column, int Np3d, int Np2d)
{
	for (int colI = 0; colI < Np2d; colI++){
		for (int RowI = 0; RowI < Np3d; RowI++)
			dest[((int)column[colI] - 1)*Np3d + RowI] += source[colI*Np3d + RowI];
	}
}
/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the row index, and has been checked*/
void AssembleContributionIntoRow(double *dest, double *source, double *Row, int Np3d, int Np2d)
{
	for (int colI = 0; colI < Np3d; colI++){
		for (int RowI = 0; RowI < Np2d; RowI++)
			dest[colI*Np3d + (int)Row[RowI] - 1] += source[colI*Np2d + RowI];
	}
}
/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the row index and the column index, and has been checked*/
void AssembleContributionIntoRowAndColumn(double *dest, double *source, double *Row, double *column, int Np3d, int Np2d, int Flag)
{
	for (int colI = 0; colI < Np2d; colI++)
	{
		for (int RowI = 0; RowI < Np2d; RowI++)
			dest[((int)column[colI] - 1)*Np3d + (int)Row[RowI] - 1] += Flag * source[colI*Np2d + RowI];
	}
}

void AssembleDataIntoPoint(double *dest, double *source, double *PIndex, int Size){
	for (int p = 0; p < Size; p++){
		dest[(int)PIndex[p] - 1] += source[p];
	}
}

/*Note: this function is used to assemble the element mass matrix and the physical diff matrix, and has been verified.
Left multiply the matrix source with a diagonal matrix composed of element contained in coe. 
*/
void DiagMultiply(double *dest, const double *source, const double *coe, int Np)
{
	for (int colI = 0; colI < Np; colI++){
		for (int RowI = 0; RowI < Np; RowI++)
			dest[colI*Np + RowI] = coe[RowI] * source[colI*Np + RowI];
	}
}

void DiagRightMultiply(double *dest, const double *source, const double *coe, int Np)
{
	for (int colI = 0; colI < Np; colI++){
		for (int RowI = 0; RowI < Np;RowI++)
			dest[colI*Np + RowI] = coe[colI] * source[colI*Np + RowI];
	}
}

void DotCriticalDivide(double *dest, double *source, double *criticalValue, double *Depth, int size){
	for (int i = 0; i < size; i++){
		if (Depth[i] >= *criticalValue)
			dest[i] = source[i] / Depth[i];
		else
			dest[i] = 0;
	}
}

void DotDivide(double *dest, double *source, double *Coefficient, int size){
	for (int i = 0; i < size; i++)
		dest[i] = source[i] / Coefficient[i];
}

void DotDivideByConstant(double *dest, double *Source, double Coefficient, int Np){
	for (int i = 0; i < Np; i++)
		dest[i] = Source[i] / Coefficient;
}

void DotProduct(double *dest, double *sourcea, double *sourceb, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] * sourceb[i];
}


void FetchInnerEdgeFacialValue(double *fm, double *fp, double *source, \
	double *FToE, double *FToN1, double *FToN2, int Np, int Nfp){
	int ind1 = ((int)FToE[0] - 1)*Np;
	int ind2 = ((int)FToE[1] - 1)*Np;
	for (int i = 0; i < Nfp; i++){
		fm[i] = source[ind1 + (int)FToN1[i] - 1];
		fp[i] = source[ind2 + (int)FToN2[i] - 1];
	}
}

void FetchBoundaryEdgeFacialValue(double *fm, double *source, \
	double *FToE, double *FToN1, int Np, int Nfp){
	int ind1 = ((int)FToE[0] - 1)*Np;
	for (int i = 0; i < Nfp; i++){
		fm[i] = source[ind1 + (int)FToN1[i] - 1];
	}
}

void Flip(double *dest, int size){
	double Temp;
	int i, j;
	i = 0;
	j = size - 1;
	while (i < j){
		Temp = dest[i];
		dest[i] = dest[j];
		dest[j] = Temp;
		++i;
		--j;
	}
}


void GetFacialFluxTerm2d(double *dest, double *hu, double *hv, double *nx, double *ny, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = hu[i] * nx[i] + hv[i] * ny[i];
}

void GetScalarFluxTerm2d(double *dest, double *C, double *Vector, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = C[i] * Vector[i];
}

void GetMeshAverageValue(double *dest, double *LAV, char *transA, char *transB, ptrdiff_t *ROPA, ptrdiff_t *COPB, ptrdiff_t *COPA, double *Alpha, double *A, \
	ptrdiff_t *LDA, double *fphys, double *Jacobian, ptrdiff_t *LDB, double *Beta, ptrdiff_t *LDC, double *wq){

	GetMeshIntegralValue(dest, transA, transB, ROPA, COPB, COPA, Alpha, A, LDA, fphys, Jacobian, LDB, Beta, LDC, wq);

	(*dest) = (*dest) / (*LAV);
}

void GetMeshIntegralValue(double *dest, char *transA, char *transB, ptrdiff_t *ROPA, ptrdiff_t *COPB, ptrdiff_t *COPA, double *Alpha, double *A, \
	ptrdiff_t *LDA, double *fphys, double *Jacobian, ptrdiff_t *LDB, double *Beta, ptrdiff_t *LDC, double *wq){

	double Jq[(int)(*LDC)], fq[(int)(*LDC)];

	// map the node values fvar to quadrature nodes by
	// \f$ fq = Vq * fvar \f$
	dgemm(transA, transB, ROPA, COPB, COPA, Alpha, A,
		LDA, fphys, LDB, Beta, fq, LDC);
	dgemm(transA, transB, ROPA, COPB, COPA, Alpha, A,
		LDA, Jacobian, LDB, Beta, Jq, LDC);

	for (int n = 0; n < (int)(*LDC); n++) {
		(*dest) += wq[n] * Jq[n] * fq[n];
	}
}

void GetVolumnIntegral1d(double *dest, ptrdiff_t *RowOPA, ptrdiff_t *ColOPB, ptrdiff_t *ColOPA, double *alpha, \
	double *Dt, ptrdiff_t *LDA, double *B, ptrdiff_t *LDB, double *Beta, ptrdiff_t *LDC, double *tz){

	dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dt, LDA, B, LDB, Beta, dest, LDC);

	DotProduct(dest, dest, tz, (int)(*LDC));


}

void GetVolumnIntegral2d(double *dest, double *tempdest, ptrdiff_t *RowOPA, ptrdiff_t *ColOPB, ptrdiff_t *ColOPA, double *alpha, \
	double *Dr, double *Ds, ptrdiff_t *LDA, double *B, ptrdiff_t *LDB, double *Beta, ptrdiff_t *LDC, double *rx, double *sx){
	/*$Dr*u(v,\theta)$*/
	dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dr, LDA, B, LDB, Beta, dest, LDC);
	/*$Ds*u(v,\theta)$*/
	dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Ds, LDA, B, LDB, Beta, tempdest, LDC);
	/*$rx\cdot Dr*u(v,\theta)$*/
	DotProduct(dest, dest, rx, (int)(*LDC));
	/*$sx\cdot Ds*u(v,\theta)$*/
	DotProduct(tempdest, tempdest, sx, (int)(*LDC));
	/*rx\cdot Dr*u(v,\theta) + sx\cdot Ds*u(v,\theta)*/
	Add(dest, dest, tempdest, (int)(*LDC));
}

void GetVolumnIntegral3d(double *dest, double *tempdest, ptrdiff_t *RowOPA, ptrdiff_t *ColOPB, ptrdiff_t *ColOPA, double *alpha, \
	double *Dr, double *Ds, double *Dt, double *E, double *G, double *H, ptrdiff_t *LDA, ptrdiff_t *LDB, double *Beta, ptrdiff_t *LDC, \
	double *rx, double *sx, double *ry, double *sy, double *tz, int Nvar, int Np, int K)
{
	for (int n = 0; n < Nvar; n++){
		/*$Dr*E$*/
		dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dr, LDA, E + n*Np*K, LDB, Beta, dest + n*Np*K, LDC);
		/*$rx\cdot Dr*E$*/
		DotProduct(dest + n*Np*K, dest + n*Np*K, rx, (int)(*LDC));
		/*$Ds*E$*/
		dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Ds, LDA, E + n*Np*K, LDB, Beta, tempdest, LDC);
		/*$sx\cdot Ds*E$*/
		DotProduct(tempdest, tempdest, sx, (int)(*LDC));
		/*rx\cdot Dr*E + sx\cdot Ds*E*/
		Add(dest + n*Np*K, dest + n*Np*K, tempdest, (int)(*LDC));
		/*$Dr*G$*/
		dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dr, LDA, G + n*Np*K, LDB, Beta, tempdest, LDC);
		/*$ry\cdot Dr*G$*/
		DotProduct(tempdest, tempdest, ry, (int)(*LDC));
		/*rx\cdot Dr*E + sx\cdot Ds*E + ry\cdot Dr*G*/
		Add(dest + n*Np*K, dest + n*Np*K, tempdest, (int)(*LDC));
		/*$Ds*G$*/
		dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Ds, LDA, G + n*Np*K, LDB, Beta, tempdest, LDC);
		/*$sy\cdot Ds*G$*/
		DotProduct(tempdest, tempdest, sy, (int)(*LDC));
		/*rx\cdot Dr*E + sx\cdot Ds*E + ry\cdot Dr*G + sy\cdot Ds*G */
		Add(dest + n*Np*K, dest + n*Np*K, tempdest, (int)(*LDC));
		/*$Dt*H$*/
		dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dt, LDA, H + n*Np*K, LDB, Beta, tempdest, LDC);
		/*$tz\cdot Dt*H$*/
		DotProduct(tempdest, tempdest, tz, (int)(*LDC));
		/*rx\cdot Dr*E + sx\cdot Ds*E + ry\cdot Dr*G + sy\cdot Ds*G + tz\cdot Dt*H*/
		Add(dest + n*Np*K, dest + n*Np*K, tempdest, (int)(*LDC));
	}
}

/*The following function is used to lump the mass matrix, i.e. sum each row of the mass matrix and put the sum at the diagonal position*/
void MassLumping(double *dest, double *source, int Np){
	memset(dest, 0, Np*Np*sizeof(double));
	double Tempdata;
	for (int row = 0; row < Np; row++){
		Tempdata = 0;
		for (int col = 0; col < Np; col++){
			Tempdata += source[row + col*Np];
		}
		dest[row + row*Np] = Tempdata;
	}
}

/*This function is used to get the inverse matrix of the given matrix, and has been verified*/
void MatrixInverse(double *dest, ptrdiff_t Np)
{
	ptrdiff_t info;
	ptrdiff_t *IPIV = malloc(Np*sizeof(ptrdiff_t));
	ptrdiff_t LWORK = Np * Np;
	double *WORK = malloc(LWORK*sizeof(double));
	dgetrf(&Np, &Np, dest, &Np, IPIV, &info);
	dgetri(&Np, dest, &Np, IPIV, WORK, &LWORK, &info);
	free(IPIV);
	free(WORK);
}

/* This function invokes dgemm implemented by blas to calculate C = alpha*A*B + beta*C */
void MatrixMultiply(char* TRANSA, char* TRANSB, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, double ALPHA, double *A,
	ptrdiff_t LDA, double *B, ptrdiff_t LDB, double BETA, double *C, ptrdiff_t LDC)	{

	dgemm(TRANSA, TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);

}

void Minus(double *dest, double *sourcea, double *sourceb, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] - sourceb[i];
}

void MultiEdgeContributionByLiftOperator(double *SrcAndDest, double *TempSource, ptrdiff_t *RowOPA, ptrdiff_t *ColOPB, ptrdiff_t *ColOPA, double *Alpha, \
	double *A, ptrdiff_t *LDA, ptrdiff_t *LDB, double *Beta, ptrdiff_t *LDC, double *J, int Np){

	dgemm("N", "N", RowOPA, ColOPB, ColOPA, Alpha, A, LDA, SrcAndDest, LDB, Beta, TempSource, \
		LDC);
	DotDivide(SrcAndDest, TempSource, J, Np);
}

void MultiplyByConstant(double *dest, double *Source, double Coefficient, int Np){
	for (int i = 0; i < Np; i++)
		dest[i] = Source[i] * Coefficient;
}

void NdgExtend2dField(double *dest, double *source, int Np2d, int Index, int Np3d, int NLayer, int Nz){
	//直接把第一层循环中的k提出来，放置到调用层，利用openmp并行，传入的K2d替换成序号k
	//for (int k = 0; k < K2d; k++){
	for (int Layer = 0; Layer < NLayer; Layer++){
		for (int N = 0; N < Nz + 1; N++){
			for (int i = 0; i < Np2d; i++){
				dest[Index*NLayer*Np3d + Layer*Np3d + N*Np2d + i] = \
					source[Index*Np2d + i];
			}
		}
	}
}

void RepmatValue(double *dest, double *source, int size, int Npz){
	for (int i = 0; i < Npz; i++){
		for (int j = 0; j < size; j++){
			dest[i*size + j] = source[j];
		}
	}
}
/*Note: source[i] equals to zero is not excluded*/
void ReverseValue(double *dest, double *source, int size){
	for (int i = 0; i < size; i++)
		dest[i] = 1.0 / source[i];
}

void Sort(double *dest, int Num){
	// Adopt bubble method to sort the given array.
	int isSorted;
	double temp;
	for (int i = 0; i<Num - 1; i++){
		//Suppose the rest data are already sorted
		isSorted = 1;
		for (int j = 0; j<Num - 1 - i; j++){
			if (dest[j] > dest[j + 1]){
				temp = dest[j];
				dest[j] = dest[j + 1];
				dest[j + 1] = temp;
				isSorted = 0;  //Once swapped, it indicates the data is not sorted
			}
		}
		if (isSorted) break; //if no swap needed, it indicates the data is already sorted.
	}
}

void StrongFormBoundaryEdgeRHS(int edgeIndex, double *FToE, double *FToF, int Np, int K, \
	int Nfp, double *FToN1, double *fluxM_, double *fluxS_, double *Js, double *Mb, double *rhs){
	const int e1 = (int)FToE[2 * edgeIndex] - 1;
	const int f1 = (int)FToF[2 * edgeIndex] - 1;
	const int ind1 = e1 * Np - 1 + f1 * Np * K;
	const int ind = edgeIndex * Nfp;
	double *rhsM = malloc(Nfp*sizeof(double));
	memset(rhsM, 0, Nfp*sizeof(double));
	for (int n = 0; n < Nfp; n++) {
		const int sk = n + ind;
		double dfM = fluxM_[sk] - fluxS_[sk];
		double j = Js[sk];
		double *mb = Mb + n * Nfp;
		for (int m = 0; m < Nfp; m++) {
			rhsM[m] += mb[m] * j * dfM;
		}
	}
	for (int n = 0; n < Nfp; n++) {
		const int sk = n + ind;
		const int m1 = (int)FToN1[sk] + ind1;
		rhs[m1] += rhsM[n];
	}
	free(rhsM);
}

void StrongFormInnerEdgeRHS(int edgeIndex, double *FToE, double *FToF, int Np, int K,\
	int Nfp, double *FToN1, double *FToN2, double *fluxM_, double *fluxP_,\
	                  double *fluxS_, double *Js, double *Mb, double *rhs){
	const int e1 = (int)FToE[2 * edgeIndex] - 1;
	const int e2 = (int)FToE[2 * edgeIndex + 1] - 1;
    const int f1 = (int)FToF[2 * edgeIndex] - 1;
    const int f2 = (int)FToF[2 * edgeIndex + 1] - 1;
	const int ind1 = e1 * Np - 1 + f1 * Np * K;
	const int ind2 = e2 * Np - 1 + f2 * Np * K;
	const int ind = edgeIndex * Nfp;
	double *rhsM = malloc(Nfp*sizeof(double));
	double *rhsP = malloc(Nfp*sizeof(double));
	memset(rhsM, 0, Nfp*sizeof(double));
	memset(rhsP, 0, Nfp*sizeof(double));
	for (int n = 0; n < Nfp; n++) {
		const int sk = n + ind;
		double dfM = fluxM_[sk] - fluxS_[sk];
		double dfP = fluxP_[sk] - fluxS_[sk];
		double j = Js[sk];
		double *mb = Mb + n * Nfp;
		for (int m = 0; m < Nfp; m++) {
			rhsM[m] += mb[m] * j * dfM;
			rhsP[m] -= mb[m] * j * dfP;
		}
	}
	for (int n = 0; n < Nfp; n++) {
	const int sk = n + ind;
	const int m1 = (int)FToN1[sk] + ind1;
	const int m2 = (int)FToN2[sk] + ind2;
	rhs[m1] += rhsM[n];
	rhs[m2] += rhsP[n];
	}
	free(rhsM);
	free(rhsP);
}