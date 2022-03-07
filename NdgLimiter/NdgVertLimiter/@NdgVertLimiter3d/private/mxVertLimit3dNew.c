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
void GetBCInvolvedLimitScope(double *, double *, double *, int , double *, double *);
void GetSBEInvolvedLimitScope(double *, double *, double *, double *, int, int, int, int, double *, int , double *, int );

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
	double Hcrit = mxGetScalar(prhs[9]);
	signed char *ftype = (signed char *)mxGetData(prhs[10]);


	/*For boundary edge object*/
	mxArray *TempBENe = mxGetField(prhs[11], 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBENfp = mxGetField(prhs[11], 0, "Nfp");
	int BENfp = mxGetScalar(TempBENfp);
	mxArray *TempBEnx = mxGetField(prhs[11], 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(prhs[11], 0, "ny");
	double *BEny = mxGetPr(TempBEny);
	mxArray *TempBEFToE = mxGetField(prhs[11], 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToN1 = mxGetField(prhs[11], 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);
	mxArray *TempBEFToF = mxGetField(prhs[11], 0, "FToF");
	double *BEFToF = mxGetPr(TempBEFToF);
	mxArray *TempBEFToV = mxGetField(prhs[11], 0, "FToV");
	double *BEFToV = mxGetPr(TempBEFToV);

	/*For bottom boundary edge object*/
	mxArray *TempBotBENe = mxGetField(prhs[12], 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBEFToE = mxGetField(prhs[12], 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToN1 = mxGetField(prhs[12], 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);
	mxArray *TempBotBENfp = mxGetField(prhs[12], 0, "Nfp");
	int BotBENfp = mxGetScalar(TempBotBENfp);
	mxArray *TempBotBEFToV = mxGetField(prhs[12], 0, "FToV");
	double *BotBEFToV = mxGetPr(TempBotBEFToV);

	/*For surface boundary edge object*/
	mxArray *TempSurfBENe = mxGetField(prhs[13], 0, "Ne");
	int SurfBENe = (int)mxGetScalar(TempSurfBENe);
	mxArray *TempSurfBEFToE = mxGetField(prhs[13], 0, "FToE");
	double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
	mxArray *TempSurfBEFToN1 = mxGetField(prhs[13], 0, "FToN1");
	double *SurfBEFToN1 = mxGetPr(TempSurfBEFToN1);
	mxArray *TempSurfBENfp = mxGetField(prhs[13], 0, "Nfp");
	int SurfBENfp = mxGetScalar(TempSurfBENfp);
	mxArray *TempSurfBEFToV = mxGetField(prhs[13], 0, "FToV");
	double *SurfBEFToV = mxGetPr(TempSurfBEFToV);


	/*Number of vertex*/
	int Nv = (int)mxGetScalar(prhs[14]);
	/*Number of elements connected to each vertex*/
	double *Nvc = mxGetPr(prhs[15]);
	/*Index of elements that connect to each studied vertex*/
	double *VToK = mxGetPr(prhs[16]);
	/*Maximum number of element connected to a vertex*/
	int maxNk = (int)mxGetM(prhs[16]);
	/*Node index on each face of the studied 2d master cell*/
	double *Fmask2d = mxGetPr(prhs[17]);
	int Nfp2d = (int)mxGetM(prhs[17]);
	int Np2d = (int)mxGetScalar(prhs[18]);
	int Nz = (int)mxGetScalar(prhs[19]);
	int Nh = (int)mxGetScalar(prhs[20]);
	/*Number of vertex of the studied 2d cell*/
	int Nv2d = (int)mxGetScalar(prhs[21]);
	double *LAV = mxGetPr(prhs[22]);
	double *EToV = mxGetPr(prhs[23]);
	/*Number of vertex of the studied 3d cell*/
	int Nv3d = (int)mxGetM(prhs[23]);

	int NLayer = (int)mxGetScalar(prhs[24]);
	int MNv2d = (int)mxGetScalar(prhs[25]);
	double *MNvc2d = mxGetPr(prhs[26]);
	double *InvVandVert = mxGetPr(prhs[27]);
	double *VandInterp = mxGetPr(prhs[28]);

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
	/*Impose boundary condition first, then to calculate the local Riemann problem*/
	/*************************************************************************************************************************************/
	             /**************************************Boundary Edge Part*******************************************************/
	                            /********************************************************************/
	/*Allocate memory for fm and fp defined over boundary edges. Here, variables correspond to hu, hv, h, hT, hS,
	sediment and other passive transport material, vertical velocity omega is not included for boundary edges*/

	double *huM = NULL, *hvM = NULL, *hM = NULL;

	double *hu = fphys, *hv = fphys + Np*K, \
		*h = fphys + 3 * Np*K, *z = fphys + 5 * Np*K;
	/*
	double *fm = malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
	huM = fm, hvM = fm + BENfp*BENe, hM = fm + 2 * BENfp*BENe;
	double *fp = malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
	double *zM = malloc(BENfp*BENe*sizeof(double));
	double *zP = malloc(BENfp*BENe*sizeof(double));
	//Since we only need the vertex value, so we choose to use the vertex value only
	double *fRiemann = malloc(4*BENe*Nvar*sizeof(double));

	// Fetch variable fm and fp first, then impose boundary condition and conduct hydrostatic reconstruction.
	//Finally, calculate the local Riemann problem to get variable at the interface
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
		FetchBoundaryEdgeFacialValue(huM + face*BENfp, hu, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*BENfp, hv, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(zM + face*BENfp, z, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		// The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
		// here 1 stands for the memory occupied by water depth h
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(fm + (field + 1)*BENe*BENfp + face*BENfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
		}
		ImposeBoundaryCondition(&gra, type, BEnx + face*BENfp, BEny + face*BENfp, fm + face*BENfp, fp + face*BENfp, \
			zM + face*BENfp, zP + face*BENfp, fext + face*BENfp, BENfp, Nvar + 1, BENe);

		EvaluateHydroStaticReconstructValue(Hcrit, fm + face*BENfp, fp + face*BENfp, zM + face*BENfp, zP + face*BENfp, BENfp, Nvar + 1, BENe);
		EvaluateVerticalFaceRiemannProblem(fRiemann + face*4, fm + face*BENfp, fp + face*BENfp, \
			BEnx + face*BENfp, BEny + face*BENfp, &gra, Hcrit, BENe, BENfp, Nvar, Nh, Nz);
	}

	
	// Boundary Edge part to be considered here, this part is used to alter fmax and fmin
	

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int n = 0; n < Nvar; n++){
		GetBCInvolvedLimitScope(fmax + n*Nv, fmin + n*Nv, BEFToF, BENe, BEFToV, fRiemann + 4 * BENe * n);
	}
	free(fm);
	free(fp);
	free(zM);
	free(zP);
	free(fRiemann);
	*/

	/*
	double *Surffm = malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
	huM = Surffm;
	hvM = Surffm + SurfBENfp*SurfBENe;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		FetchBoundaryEdgeFacialValue(huM + face*SurfBENfp, hu, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*SurfBENfp, hv, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(Surffm + field*SurfBENe*SurfBENfp + face*SurfBENfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				SurfBEFToE + 2 * face, SurfBEFToN1 + SurfBENfp*face, Np, SurfBENfp);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int n = 0; n < Nvar; n++){
		GetSBEInvolvedLimitScope(fmax + n*Nv, fmin + n*Nv, SurfBEFToV, Surffm, SurfBENe, Np2d, Nfp2d, Nv2d, Fmask2d, MNv2d, MNvc2d, 0);
	}

	free(Surffm);
*/
	double *BotBEfm = malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
	huM = BotBEfm;
	hvM = BotBEfm + BotBENfp*BotBENe;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		FetchBoundaryEdgeFacialValue(huM + face*BotBENfp, hu, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*BotBENfp, hv, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		for (int field = 2; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(BotBEfm + field*BotBENe*BotBENfp + face*BotBENfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				BotBEFToE + 2 * face, BotBEFToN1 + BotBENfp*face, Np, BotBENfp);
		}
	}
	//int Nv, double *Nvc2d, int NvB
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int n = 0; n < Nvar; n++){
		GetSBEInvolvedLimitScope(fmax + n*Nv, fmin + n*Nv, BotBEFToV, BotBEfm, BotBENe, Np2d, Nfp2d, Nv2d, Fmask2d, MNv2d, MNvc2d, NLayer*MNv2d);
	}

	free(BotBEfm);
	
	
/******************************************************************Boundary edge part finished************************************************************************************/

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
				double *LVertValue = malloc(2 * Nv2d*sizeof(double));
				/*Here we alter the vertex value first, then we calculate the mode coefficient according to the vertex value, 
				* then we multiply the mode coefficient by the vandemonde matrix corresponding to the first order to regain the 
				* intepolation point date value
				*/
				for (int L = 0; L < 2; L++){
					for (int i = 0; i < Nv2d; i++){
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

	/*Calculate the average value first*/
/*	double *AfterAve = malloc(K*Nvar*sizeof(double));
	memset(AfterAve, 0, K*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int var = 0; var < Nvar; var++){
			GetMeshAverageValue(AfterAve + var*K + k, LAV + k, transA, transB, &ROPA, &COPB, &COPA, &Alpha, A, \
				&LDA, flimit + var*Np*K + k*Np, J + k*Np, &LDB, &Beta, &LDC, wq);
			printf("The difference is:%f\n", pow(10.0,12.0)*(*(AfterAve + var*K + k) - *(Ave + var*K + k)) + 1.0);
		}
	}

	free(AfterAve);
	*/

	free(Ave);
	free(fmin);
	free(fmax);
}

/*This function is used to get the minimum slope parameter*/
void GetMinimumSlope(double *dest, double *LamMax, double *LamMin, int Nv){
	*dest = INF;
	for (int i = 0; i < Nv; i++){
		*dest = min(*dest, LamMax[i]);
		*dest = min(*dest, LamMin[i]);
	}
	
}

/*
* Purpose: This function is used to get the minimum and maximum physical value at the vertex for boundary edge
*
* Input:
*      double[Nv x Ne] fmax the maximum value defined at the vertex according to the average value of cells adjacent to the studied vertex
*      double[Nv x Ne] fmin the minimum value defined at the vertex according to the average value of cells adjacent to the studied vertex
*      double[2 x Ne] FToF the face to face relationship at the boundary
* 	   int Ne number of face located on boundary edges.
*      double[4 x Ne] FToV the global vertex index for boundary edge
*      double[4 x Ne] fRiemann the Riemann solution of variables for boundary edge
* Output:
* 		double[Nv x Ne] fmax the maximum value defined at the vertex with boundary condition considered
*       double[Nv x Ne] fmin the minimum value defined at the vertex with boundary condition considered
*/

void GetBCInvolvedLimitScope(double *fmax, double *fmin, double *FToF, int Ne, double *FToV, double *fRiemann){
	int v;
	for (int i = 0; i < Ne; i++){
		int face = (int)FToF[i];
		if (face >= 3){
			v = (int)FToV[4 * i]-1;
			fmax[v] = max(fmax[v], fRiemann[4 * i + 1]);
			fmin[v] = min(fmin[v], fRiemann[4 * i + 1]);
			v = (int)FToV[4 * i + 1]-1;
			fmax[v] = max(fmax[v], fRiemann[4 * i]);
			fmin[v] = min(fmin[v], fRiemann[4 * i]);
			v = (int)FToV[4 * i + 2]-1;
			fmax[v] = max(fmax[v], fRiemann[4 * i + 2]);
			fmin[v] = min(fmin[v], fRiemann[4 * i + 2]);
			v = (int)FToV[4 * i + 3]-1;
			fmax[v] = max(fmax[v], fRiemann[4 * i + 3]);
			fmin[v] = min(fmin[v], fRiemann[4 * i + 3]);
		}
		else{
			v = (int)FToV[4 * i]-1;
			fmax[v] = max(fmax[v], fRiemann[4 * i]);
			fmin[v] = min(fmin[v], fRiemann[4 * i]);
			v = (int)FToV[4 * i + 1]-1;
			fmax[v] = max(fmax[v], fRiemann[4 * i + 1]);
			fmin[v] = min(fmin[v], fRiemann[4 * i + 1]);
			v = (int)FToV[4 * i + 2]-1;
			fmax[v] = max(fmax[v], fRiemann[4 * i + 3]);
			fmin[v] = min(fmin[v], fRiemann[4 * i + 3]);
			v = (int)FToV[4 * i + 3]-1;
			fmax[v] = max(fmax[v], fRiemann[4 * i + 2]);
			fmin[v] = min(fmin[v], fRiemann[4 * i + 2]);
		}
	}
}


/*
* Purpose: This function is used to get the minimum and maximum physical value at the surface and bottom boundary vertex
*
* Input:
*      double[Nv x Ne] fmax the maximum value defined at the vertex according to the average value of cells adjacent to the studied vertex
*      double[Nv x Ne] fmin the minimum value defined at the vertex according to the average value of cells adjacent to the studied vertex
*      double[Nv2d x BENe]  BEFToV the global vertex index for boundary edge
* 	   double[Nfp2d x BENe]  fm the local physical value defined over the boundary edge
* 	   int BENe number of elements located on boundary edges. Note: boundary edges here means surface boundary edge and bottom boundary edge
*      int Np2d number of interpolation points for boundary edge element
*      int Nfp2d number of interpolation points located on each edge of the boundary edge element
*      int Nv2d number of vertex for each boundary edge element
*      double[Nfp2d x Nv2d] Fmask index of the interpolation points on edge of the two dimensional master cell
* Output:
* 		double[Nv x Ne] fmax the maximum value defined at the vertex with boundary condition considered
*       double[Nv x Ne] fmin the minimum value defined at the vertex with boundary condition considered
*/

void GetSBEInvolvedLimitScope(double *fmax, double *fmin, double *BEFToV, double *fm, int BENe, int Np2d, int Nfp2d, int Nv2d, double *Fmask2d, int MNv2d, double *Nvc2d, int NvB){
	int v, Index;
	double Data;
	double *Ave = malloc(MNv2d*sizeof(double));
	memset(Ave, 0, MNv2d*sizeof(double));
	//How to calculate the average value
	if (Nv2d == 4){
		for (int i = 0; i < BENe; i++){
			for (int j = 0; j < 2; j++){
				v = (int)BEFToV[i*Nv2d + j] - 1 - NvB;
				Index = (int)Fmask2d[j * Nfp2d] - 1;
				Data = fm[i*Np2d + Index];
				Ave[v] += Data;
			}
			/*The left two vertex, because their order in FToV and Fmask is not the same*/
			v = (int)BEFToV[i*Nv2d + 2] - 1 - NvB;
			Index = (int)Fmask2d[3 * Nfp2d] - 1;
			Data = fm[i*Np2d + Index];
			Ave[v] += Data;

			v = (int)BEFToV[i*Nv2d + 3] - 1 - NvB;
			Index = (int)Fmask2d[2 * Nfp2d] - 1;
			Data = fm[i*Np2d + Index];
			Ave[v] += Data;
		}
	}
	else if (Nv2d == 3){
		for (int i = 0; i < BENe; i++){
			for (int j = 0; j < Nv2d; j++){
				v = (int)BEFToV[i*Nv2d + j] - 1 - NvB;
				Index = (int)Fmask2d[j * Nfp2d] - 1;
				Data = fm[i*Np2d + Index];
				Ave[v] += Data;
			}
		}
	}

	for (int i = 0; i < MNv2d; i++){
		Ave[i] = Ave[i] / Nvc2d[i];
	}

	for (int i = 0; i < BENe; i++){
		for (int j = 0; j < Nv2d; j++){
			v = (int)BEFToV[i*Nv2d + j] - 1;
			fmax[v] = max(fmax[v], Ave[v - NvB]);
			fmin[v] = min(fmin[v], Ave[v - NvB]);
		}
	}
	free(Ave);
}
