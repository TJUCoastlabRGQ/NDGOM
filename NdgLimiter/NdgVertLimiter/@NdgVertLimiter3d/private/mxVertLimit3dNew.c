#include "../../../../NdgMath/NdgMath.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#define INF 10.0e9

void GetMinimumSlope(double *, double *, double *, int);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *fphys = mxGetPr(prhs[0]);
	double *varFieldIndex = mxGetPr(prhs[1]);
	int Nvar = (int)mxGetNumberOfElements(prhs[1]);
	int Np = (int)mxGetScalar(prhs[2]);
	int K = (int)mxGetScalar(prhs[3]);
	double *J = mxGetPr(prhs[4]);
	double *wq = mxGetPr(prhs[5]);
	double *Vq = mxGetPr(prhs[6]);
	int RVq = (int)mxGetM(prhs[6]);
	int CVq = (int)mxGetN(prhs[6]);
	double *fext = mxGetPr(prhs[7]);
	double gra = mxGetScalar(prhs[8]);
	double hcrit = mxGetScalar(prhs[9]);
	signed char *BEftype = (signed char *)mxGetData(prhs[10]);
	signed char *BBEftype = (signed char *)mxGetData(prhs[11]);
	signed char *SBEftype = (signed char *)mxGetData(prhs[12]);

	/*
	13 14 15
	BoundaryEdge, BottomBoundaryEdge, SurfaceBoundaryEdge
	to be read as required
	*/
	/*Number of vertex*/
	int Nv = (int)mxGetScalar(prhs[16]);
	/*Number of elements connected to each vertex*/
	double *Nvc = mxGetPr(prhs[17]);
	/*Index of elements that connect to each studied vertex*/
	double *VToK = mxGetPr(prhs[18]);
	/*Maximum number of element connected to a vertex*/
	int maxNk = (int)mxGetM(prhs[18]);
	/*Node index on each face of the studied 2d master cell*/
	double *Fmask2d = mxGetPr(prhs[19]);
	int Nfp2d = (int)mxGetM(prhs[19]);
	int Np2d = (int)mxGetScalar(prhs[20]);
	int Nz = (int)mxGetScalar(prhs[21]);
	/*Number of vertex of the studied 2d cell*/
	int Nv2d = (int)mxGetScalar(prhs[22]);
	double *LAV = mxGetPr(prhs[23]);
	double *EToV = mxGetPr(prhs[24]);
	/*Number of vertex of the studied 3d cell*/
	int Nv3d = (int)mxGetM(prhs[24]);
	//double *Ave = mxGetPr(prhs[25]);

	size_t NdimOut = 3;
	mwSize dimOut[3] = { Np, K, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *flimit = mxGetPr(plhs[0]);
	
	char *transA = "N";
	char *transB = "N";
	ptrdiff_t ROPA = (ptrdiff_t)RVq;
	ptrdiff_t COPA = (ptrdiff_t)CVq;
	ptrdiff_t COPB = 1;
	double Alpha = 1.0;
	double *A = Vq;
	ptrdiff_t LDA = (ptrdiff_t)RVq;
	ptrdiff_t LDB = (ptrdiff_t)RVq;
	ptrdiff_t LDC = LDA;
	double Beta = 0.0;

	/*Calculate the average value first*/
	double *Ave = malloc(K*Nvar*sizeof(double));
	memset(Ave, 0, K*Nvar*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int var = 0; var < Nvar; var++){
			GetMeshAverageValue(Ave + var*K + k, LAV + k, transA, transB, &ROPA, &COPB, &COPA, &Alpha, A,\
				&LDA, fphys + ((int)varFieldIndex[var]-1)*Np*K + k*Np, J + k*Np, &LDB, &Beta, &LDC, wq);
		}
	}

	/*Calculate the fvmin and fvmax without considering boundary condition*/
	double *fmin = malloc(Nv*Nvar*sizeof(double));
	double *fmax = malloc(Nv*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int n = 0; n < Nv; n++){
		for (int var = 0; var < Nvar; var++){
			int Nk = (int)Nvc[n]; // number of cells connecting to vertex n
			fmax[var*Nv+n] = -INF;
			fmin[var*Nv+n] = INF;
			for (int k = 0; k<Nk; k++) {
				int cellId = (int)VToK[n*maxNk + k] - 1;
				double temp = Ave[var*K + cellId];
				fmax[var*Nv + n] = max(fmax[var*Nv + n], temp);
				fmin[var*Nv + n] = min(fmin[var*Nv + n], temp);
			}
		}
	}

/*Boundary condition to be considered here next*/


	/*
	Boundary Edge part to be considered here, this part is used to alter fmax and fmin
	*/

	/*Limit the physical value according to fmax and fmin*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int var = 0; var < Nvar; var++){
			double *Data = fphys + ((int)varFieldIndex[var] - 1)*Np*K;
			/*Copy the original data into flimit, if the cell is not troubled, just output it as original*/
			for (int i = 0; i < Np; i++){
				flimit[var*Np*K + k*Np + i] = Data[k*Np + i];
			}

			/*Identify whether the studied cell needs to be corrected or not*/
			int flag = 0;
			/*study the bottomost face and uppermost face, respectively*/
			for (int L = 0; L < 2; L++){
				for (int i = 0; i < Nv2d; i++){
					/*Global index of the studied vertex of the studied cell*/
					int nodeId = k * Np + L*Np2d*Nz + (int)Fmask2d[i * Nfp2d] - 1;
					int vertId = (int)EToV[k * Nv3d + L*Nv2d + i] - 1;
					/*If value at the studied vertex is greater than the maxmium or smaller than the minimum value*/
					if (Data[nodeId] > fmax[var*Nv + vertId] || Data[nodeId] < fmin[var*Nv + vertId]){
						flag = 1;
					}
				}
			}
			/*if the cell is flagged as a troubled one*/
			if (flag){
				double *LambdaMax = malloc(Nv3d*sizeof(double));
				double *LambdaMin = malloc(Nv3d*sizeof(double));
				for (int L = 0; L < 2; L++){
					for (int i = 0; i < Nv2d; i++){
						/*Global index of the studied vertex of the studied cell*/
						int nodeId = k * Np + L*Np2d*Nz + (int)Fmask2d[i * Nfp2d] - 1;
						int vertId = (int)EToV[k * Nv3d + L*Nv2d + i] - 1;
						/*Calculate the maxmum slope parameter to satisfy the required maxmum/minimum requirement at each node*/
						if (Data[nodeId] > fmax[var*Nv + vertId]){
							LambdaMax[L*Nv2d + i] = (fmax[var*Nv + vertId] - Ave[var*K + k]) / (Data[nodeId] - Ave[var*K + k]);
							LambdaMin[L*Nv2d + i] = 1.0;
						} 
						else if (Data[nodeId] < fmin[var*Nv + vertId]){
							LambdaMin[L*Nv2d + i] = (fmin[var*Nv + vertId] - Ave[var*K + k]) / (Data[nodeId] - Ave[var*K + k]);
							LambdaMax[L*Nv2d + i] = 1.0;
						}
						else{
							LambdaMin[L*Nv2d + i] = 1.0;
							LambdaMax[L*Nv2d + i] = 1.0;
						}
					}
				}
				double Lambda;
				/*Get the mininum Lambda*/
				GetMinimumSlope(&Lambda, LambdaMax, LambdaMin, Nv3d);
				/*Here we only alter the vertex value only, if high order case is considered, we need to change this part*/
				for (int L = 0; L < 2; L++){
					for (int i = 0; i < Nv2d; i++){
						/*Global index of the studied vertex of the studied cell*/
						int nodeId = k * Np + L*Np2d*Nz + (int)Fmask2d[i * Nfp2d] - 1;
						flimit[var*Np*K + nodeId] = Lambda * flimit[var*Np*K + nodeId] + (1 - Lambda)*Ave[var*K + k];
					}
				}
				free(LambdaMax);
				free(LambdaMin);
			}
		}
	}
	free(Ave);
	free(fmin);
	free(fmax);
}

/*This function is used to get the minimum slope parameter*/
void GetMinimumSlope(double *dest, double *LamMax, double *LamMin, int Nv){
	*dest = 100;
	for (int i = 0; i < Nv; i++){
		*dest = min(*dest, LamMax[i]);
		*dest = min(*dest, LamMin[i]);
	}
	
}