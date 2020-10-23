#include "NdgMath.h"
//#include <omp.h>
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

void RepmatValue(double *dest, double *source, int Layer){
	for (int i = 0; i < Layer; i++)
		dest[i] = *source;
}

void DotProduct(double *dest, double *sourcea, double *sourceb, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] * sourceb[i];
}

void DotDivide(double *dest, double *source, double *Coefficient, int size){
	for (int i = 0; i < size; i++)
		dest[i] = source[i] / Coefficient[i];
}

void DotCriticalDivide(double *dest, double *source, double *criticalValue, double *Depth, int size){
	for (int i = 0; i < size; i++){
		if (Depth[i] >= *criticalValue)
			dest[i] = source[i] / Depth[i];
		else
			dest[i] = 0;
	}
}

void DotDivideByConstant(double *dest, double *Source, double Coefficient, int Np){
	for (int i = 0; i < Np; i++)
		dest[i] = Source[i] / Coefficient;
}


void Add(double *dest, double *sourcea, double *sourceb, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] + sourceb[i];
}

void Minus(double *dest, double *sourcea, double *sourceb, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] - sourceb[i];
}


void StrongFormInnerEdgeRHS(int edgeIndex, double *FToE, int Np,\
	int Nfp, double *FToN1, double *FToN2, double *fluxM_, double *fluxP_,\
	                  double *fluxS_, double *Js, double *Mb, double *rhs){
	const int e1 = (int)FToE[2 * edgeIndex] - 1;
	const int e2 = (int)FToE[2 * edgeIndex + 1] - 1;
	const int ind1 = e1 * Np - 1;
	const int ind2 = e2 * Np - 1;
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

void StrongFormBoundaryEdgeRHS(int edgeIndex, double *FToE, int Np, int Nfp, \
	double *FToN1, double *fluxM_, double *fluxS_, double *Js, double *Mb, double *rhs){
	const int e1 = (int)FToE[2 * edgeIndex] - 1;
	const int ind1 = e1 * Np - 1;
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

void NdgExtend2dField(double *dest, double *source, int Np2d, int Index, int Np3d, int K3d, int NLayer, int Nz){
	//ֱ�Ӱѵ�һ��ѭ���е�k����������õ����ò㣬����openmp���У������K2d�滻�����k
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

void GetFacialFluxTerm2d(double *dest, double *hu, double *hv, double *nx, double *ny, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = hu[i] * nx[i] + hv[i]*ny[i];
}

void MultiEdgeContributionByLiftOperator(double *SrcAndDest, double *TempSource, ptrdiff_t *RowOPA, ptrdiff_t *ColOPB, ptrdiff_t *ColOPA, double *Alpha, \
	double *A, ptrdiff_t *LDA, ptrdiff_t *LDB, double *Beta, ptrdiff_t *LDC, double *J, int Np){

	dgemm("N", "N", RowOPA, ColOPB, ColOPA, Alpha, A, LDA, SrcAndDest, LDB, Beta, TempSource, \
		LDC);
	DotDivide(SrcAndDest, TempSource, J, Np);
}