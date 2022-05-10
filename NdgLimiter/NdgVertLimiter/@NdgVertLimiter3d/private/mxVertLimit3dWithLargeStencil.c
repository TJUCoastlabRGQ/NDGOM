#include "../../../../NdgMath/NdgMath.h"
#include "../../../../NdgMath/NdgSWE.h"
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#define INF 10.0e9


#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)


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

	/*Number of vertex*/
	int Nv = (int)mxGetScalar(prhs[7]);
	/*Number of elements connected to each vertex*/
	double *Nvc = mxGetPr(prhs[8]);
	/*Index of elements that connect to each studied vertex*/
	double *VToK = mxGetPr(prhs[9]);
	/*Maximum number of element connected to a vertex*/
	int maxNk = (int)mxGetM(prhs[9]);
	/*Node index on each face of the studied 2d master cell*/
	double *Fmask2d = mxGetPr(prhs[10]);
	int Nfp2d = (int)mxGetM(prhs[10]);
	int Np2d = (int)mxGetScalar(prhs[11]);
	int Nz = (int)mxGetScalar(prhs[12]);
	int Nh = (int)mxGetScalar(prhs[13]);
	/*Number of vertex of the studied 2d cell*/
	int Nv2d = (int)mxGetScalar(prhs[14]);
	double *LAV = mxGetPr(prhs[15]);
	double *EToV = mxGetPr(prhs[16]);
	/*Number of vertex of the studied 3d cell*/
	int Nv3d = (int)mxGetM(prhs[16]);

	int NLayer = (int)mxGetScalar(prhs[17]);
	int MNv2d = (int)mxGetScalar(prhs[18]);
	double *MNvc2d = mxGetPr(prhs[19]);
	double *InvVandVert = mxGetPr(prhs[20]);
	double *VandInterp = mxGetPr(prhs[21]);

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
	int k, var;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS) private(var) if(K>1000)
#endif
	for (k = 0; k < K; k++){
		for (var = 0; var < Nvar; var++){
			GetMeshAverageValue(Ave + var*K + k, LAV + k, transA, transB, &ROPA, &COPB, &COPA, &Alpha, A,\
				&LDA, fphys + ((int)varFieldIndex[var]-1)*Np*K + k*Np, J + k*Np, &LDB, &Beta, &LDC, wq);
		}
	}

	/*Calculate the fvmin and fvmax without considering boundary condition*/
	double *fmin = malloc(Nv*Nvar*sizeof(double));
	double *fmax = malloc(Nv*Nvar*sizeof(double));
	int n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS) private(var, k) if(Nv>1000)
#endif
	for (n = 0; n < Nv; n++){
		for (var = 0; var < Nvar; var++){
			int Nk = (int)Nvc[n]; // number of cells connecting to vertex n
			fmax[var*Nv+n] = -INF;
			fmin[var*Nv+n] = INF;
			for (k = 0; k<Nk; k++) {
				int cellId = (int)VToK[n*maxNk + k] - 1;
				double temp = Ave[var*K + cellId];
				fmax[var*Nv + n] = max(fmax[var*Nv + n], temp);
				fmin[var*Nv + n] = min(fmin[var*Nv + n], temp);
			}
		}
	}
	
	int i, L;
	/*Limit the physical value according to fmax and fmin*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS) private(var,i,L) if(K>1000)
#endif
	for (k = 0; k < K; k++) {
		for (var = 0; var < Nvar; var++){
			double *Data = fphys + ((int)varFieldIndex[var] - 1)*Np*K;
			/*Copy the original data into flimit, if the cell is not troubled, just output it as original*/
			for (i = 0; i < Np; i++){
				flimit[var*Np*K + k*Np + i] = Data[k*Np + i];
			}

			/*Identify whether the studied cell needs to be corrected or not*/
			int flag = 0;
			/*study the bottomost face and uppermost face, respectively*/
			for (L = 0; L < 2; L++){
				for (i = 0; i < Nv2d; i++){
					/*Global index of the studied vertex of the studied cell. Here, nodeId and vertId must refer to the same point*/
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
				for (L = 0; L < 2; L++){
					for (i = 0; i < Nv2d; i++){
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
				double *LVertValue = malloc(2 * Nv2d*sizeof(double));
				/*Here we alter the vertex value first, then we calculate the mode coefficient according to the vertex value, 
				* then we multiply the mode coefficient by the vandemonde matrix corresponding to the first order to regain the 
				* intepolation point date value
				*/
				for (L = 0; L < 2; L++){
					for (i = 0; i < Nv2d; i++){
						/*Global index of the studied vertex of the studied cell*/
						int nodeId = k * Np + L*Np2d*Nz + (int)Fmask2d[i * Nfp2d] - 1;
						LVertValue[L*Nv2d + i] = Lambda * flimit[var*Np*K + nodeId] + (1 - Lambda)*Ave[var*K + k];
						//flimit[var*Np*K + nodeId] = Lambda * flimit[var*Np*K + nodeId] + (1 - Lambda)*Ave[var*K + k];
					}
				}
				double alpha = 1.0, beta = 0.0;

				ptrdiff_t RowA, ColB, ColA;

				double *fmod = malloc(2 * Nv2d*sizeof(double));

				RowA = 2 * Nv2d, ColB = 1, ColA = 2 * Nv2d;

				MatrixMultiply("n", "n", RowA, ColB, ColA, alpha, InvVandVert, \
					RowA, LVertValue, RowA, beta, fmod, RowA);

				RowA = Np, ColB = 1, ColA = 2 * Nv2d;

				MatrixMultiply("n", "n", RowA, ColB, ColA, alpha, VandInterp, \
					RowA, fmod, RowA, beta, flimit + var*Np*K + k*Np, RowA);

				free(LVertValue);
				free(fmod);
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
	*dest = INF;
	int i;
	for (i = 0; i < Nv; i++){
		*dest = min(*dest, LamMax[i]);
		*dest = min(*dest, LamMin[i]);
	}
}
