#include "../../../../NdgMath/NdgMath.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#define INF 10.0e9

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
	int Nv = (int)mxGetScalar(prhs[16]);
	double *Nvc = mxGetPr(prhs[17]);
	double *VToK = mxGetPr(prhs[18]);
	int maxNk = (int)mxGetM(prhs[18]);
	double *Fmask2d = mxGetPr(prhs[19]);
	int Nfp2d = (int)mxGetM(prhs[19]);
	int Np2d = (int)mxGetScalar(prhs[20]);
	int Nz = (int)mxGetScalar(prhs[21]);
	int Nv2d = (int)mxGetScalar(prhs[22]);
	
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

	size_t NdimOut = 3;
	mwSize dimOut[3] = { Np, K, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *flimit = mxGetPr(plhs[0]);

	/*Calculate the average alue first*/
	double *Ave = malloc(K*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int var = 0; var < Nvar; var++){
			GetMeshAverageValue(var*K + k, LAV + k, transA, transB, &ROPA, &COPB, &COPA, &Alpha, A,\
				&LDA, fphys + (varFieldIndex[var]-1)*Np*K + k*Np, Jacobian + k*Np, &LDB, &Beta, &LDC, wq);
		}
	}

	/*Calculate the fvmin and fvmax without considering boundary condition*/
	double *fmin = malloc(Nv*Nvar*sizeof(double));
	double *fmax = malloc(Nv*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < Nv; i++){
		for (int var = 0; var < Nvar; var++){
			int Nk = (int)Nvc[n]; // number of cells connecting to vertex n
			fmax[var*K+n] = -INF;
			fmin[var*K+n] = INF;
			for (int k = 0; k<Nk; k++) {
				int cellId = (int)VToK[n*maxNk + k] - 1;
				double temp = Ave[var*K + cellId];
				fmax[var*K + n] = max(fmax[var*K + n], temp);
				fmin[var*K + n] = min(fmin[var*K + n], temp);
			}
		}
	}

/*Boundary condition to be considered here next*/


	/*
	*
	Boundary Edge part to be considered here, this part is used to alter fmax and fmin
	*
	*
	*/

	/*Limit the physical value according to fmax and fmin*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int var = 0; var < Nvar; var++){
			double *Data = fphys + (varFieldIndex[var] - 1)*Np*K;

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
					int vertId = (int)EToV[k * Nv + L*Nv2d + i] - 1;
					/*If value at the studied vertex is greater than the maxmium or smaller than the minimum value*/
					if (Data[nodeId] > fmax[vertId] || Data[nodeId] < fmin[vertId]){
						flag = 1;
					}
				}
			}
			/*if the cell is flagged as a troubled one*/
			if (flag){
				double LambdaMax[Nv];
				double LambdaMin[Nv];
				for (int L = 0; L < 2; L++){
					for (int i = 0; i < Nv2d; i++){
						/*Global index of the studied vertex of the studied cell*/
						int nodeId = k * Np + L*Np2d*Nz + (int)Fmask2d[i * Nfp2d] - 1;
						int vertId = (int)EToV[k * Nv + L*Nv2d + i] - 1;
						/*Calculate the maxmum slope parameter to satisfy the required maxmum/minimum requirement at each node*/
						if (Data[nodeId] > fmax[vertId]){
							LambdaMax[L*Nv2d + i] = (fmax[vertId] - Ave[var*K + k]) / (Data[nodeId] - Ave[var*K + k]);
							LambdaMin[L*Nv2d + i] = 1.0;
						} 
						else {
							LambdaMin[L*Nv2d + i] = (fmin[vertId] - Ave[var*K + k]) / (Data[nodeId] - Ave[var*K + k]);
							LambdaMax[L*Nv2d + i] = 1.0;
						}
					}
				}
				double Lambda;
				/*Get the mininum Lambda*/
				GetMinimumSlope(&Lambda, LambdaMax, LambdaMin, Nv);
				/*Here we only alter the vertex value only, if high order case is considered, we need to change this part*/
				for (int L = 0; L < 2; L++){
					for (int i = 0; i < Nv2d; i++){
						/*Global index of the studied vertex of the studied cell*/
						int nodeId = k * Np + L*Np2d*Nz + (int)Fmask2d[i * Nfp2d] - 1;
						flimit[var*Np*K + nodeId] = Lambda * flimit[var*Np*K + nodeId] + (1 - Lambda)*Ave[var*K + k]
					}
				}
			}
		}
	}
}

/*This function is used to get the minimum slope parameter*/
void GetMinimumSlope(double *dest, double *LamMax, double *LamMin, int Nv){
	*dest = 100;
	for (int i = 0; i < Nv; i++){
		*dest = min(*dest, LamMax[i]);
		*dest = min(*dest, LamMin[i]);
	}
	
}