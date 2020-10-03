#include "NdgMath.h"
//#include <omp.h>
void FetchInnerEdgeFacialValue(double *fm, double *fp, double *source, \
	double *FToE, double *FToN1, double *FToN2, int Np, int Nfp){
	int ind1 = ((int)FToE[0] - 1)*Np;
	int ind2 = ((int)FToE[1] - 1)*Np;
	for (int i = 0; i < Nfp; i++){
		fm[i] = source[ind1 + (int)FToN1[i]];
		fp[i] = source[ind2 + (int)FToN2[i]];
	}
}

void FetchBoundaryEdgeFacialValue(double *fm, double *source, \
	double *FToE, double *FToN1, int Np, int Nfp){
	int ind1 = ((int)FToE[0] - 1)*Np;
	for (int i = 0; i < Nfp; i++){
		fm[i] = source[ind1 + (int)FToN1[i]];
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
#pragma omp atomic  
	rhs[m1] += rhsM[n];
#pragma omp atomic  
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
#pragma omp atomic  
		rhs[m1] += rhsM[n];
	}
	free(rhsM);
}
