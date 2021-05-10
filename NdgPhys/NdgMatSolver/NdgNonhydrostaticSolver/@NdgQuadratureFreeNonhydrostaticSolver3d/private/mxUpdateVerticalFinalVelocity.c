#include "SWENonhydrostatic3d.h"

void FetchBoundaryData(double *dest, double *source, double *destIndex, double *sourceIndex, int size)
{
	for (int i = 0; i < size; i++)
		dest[(int)destIndex[i] - 1] = source[(int)sourceIndex[i] - 1];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	double *fphys = mxGetPr(prhs[0]);
	double *varIndex = mxGetPr(prhs[1]);
	double *Hcrit = mxGetPr(prhs[2]);
	double *PSPX = mxGetPr(prhs[3]);
	double *PSPY = mxGetPr(prhs[4]);
	double *PUPX = mxGetPr(prhs[5]);
	double *PVPY = mxGetPr(prhs[6]);
	double *PUPS = mxGetPr(prhs[7]);
	double *PVPS = mxGetPr(prhs[8]);
	double *RHSCoeMatrix = mxGetPr(prhs[9]);
	double *VertCoeMatrix = mxGetPr(prhs[10]);
	const mxArray *mesh = prhs[11];
	const mxArray *cell = prhs[12];
	const mxArray *BottomBoundaryEdge = prhs[13];
	double *PUVPXY = mxGetPr(prhs[14]);

	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);

	mxArray *TempNfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
	int Np2d = (int)mxGetScalar(TempNfp);

	mxArray *TempNe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int K2d = (int)mxGetScalar(TempNe);

	mxArray *TempNLayer = mxGetField(mesh, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);

	double *BotEidM = malloc(Np2d*sizeof(double));
	double *UpEidM = malloc(Np2d*sizeof(double));

	mxArray *TempFmask = mxGetField(cell, 0, "Fmask");
	int maxNfp = (int)mxGetM(TempFmask);
	int Nface = (int)mxGetN(TempFmask);
	double *Fmask = mxGetPr(TempFmask);

	for (int i = 0; i < Np2d; i++){
		BotEidM[i] = Fmask[(Nface-2)*maxNfp+i];
		UpEidM[i] = Fmask[(Nface-1)*maxNfp+i];
	}

	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *VerticalVelocity = mxGetPr(plhs[0]);
	double *VSTempVerticalVelocity = malloc(Np*K*sizeof(double));
	double *RHS = malloc(Np*K*sizeof(double));
	memset(RHS, 0, Np*K*sizeof(double));

	double *BottomVerticalVelocity = malloc(Np2d*K2d*sizeof(double));
	double *TempBottomVertVelocityX = malloc(Np2d*K2d*sizeof(double));
	double *TempBottomVertVelocityY = malloc(Np2d*K2d*sizeof(double));
	double *ubot = malloc(Np2d * K2d*sizeof(double));
	double *vbot = malloc(Np2d*K2d*sizeof(double));
	double *Hbot = malloc(Np2d*K2d*sizeof(double));
	double *zx = malloc(Np2d*K2d*sizeof(double));
	double *zy = malloc(Np2d*K2d*sizeof(double));

	double *Hu = fphys + ((int)varIndex[0] - 1)*Np*K;
	double *Hv = fphys + ((int)varIndex[1] - 1)*Np*K;
	double *H = fphys + ((int)varIndex[3] - 1)*Np*K;
	double *zx3d = fphys + ((int)varIndex[4] - 1)*Np*K;
	double *zy3d = fphys + ((int)varIndex[5] - 1)*Np*K;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int i = 0; i < Np; i++){
//			RHS[k*Np + i] = -1.0 * H[k*Np + i] * \
				(PUPX[k*Np + i] + PVPY[k*Np + i] + PSPX[k*Np + i] * PUPS[k*Np + i] + \
				PSPY[k*Np + i] * PVPS[k*Np + i]);
			RHS[k*Np + i] = -1.0 * H[k*Np + i] * \
				(PUVPXY[k*Np + i] + PSPX[k*Np + i] * PUPS[k*Np + i] + \
				PSPY[k*Np + i] * PVPS[k*Np + i]);
		}
	}

	mxArray *TempBotBENfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
	int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBEFToE = mxGetField(BottomBoundaryEdge, 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToN1 = mxGetField(BottomBoundaryEdge, 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		FetchBoundaryEdgeFacialValue(ubot + face*BotBENfp, Hu, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

		FetchBoundaryEdgeFacialValue(vbot + face*BotBENfp, Hv, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

		FetchBoundaryEdgeFacialValue(Hbot + face*BotBENfp, H, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

		FetchBoundaryEdgeFacialValue(zx + face*BotBENfp, zx3d, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

		FetchBoundaryEdgeFacialValue(zy + face*BotBENfp, zy3d, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

		//For variable u, v, and w
		DotCriticalDivide(ubot + face*BotBENfp, ubot + face*BotBENfp, Hcrit, Hbot + face*BotBENfp, BotBENfp);
		DotCriticalDivide(vbot + face*BotBENfp, vbot + face*BotBENfp, Hcrit, Hbot + face*BotBENfp, BotBENfp);

		DotProduct(TempBottomVertVelocityX + face*BotBENfp, ubot + face*BotBENfp, zx + face*BotBENfp, BotBENfp);

		DotProduct(TempBottomVertVelocityY + face*BotBENfp, vbot + face*BotBENfp, zy + face*BotBENfp, BotBENfp);

		/*$w_{bot} = \frac{\partial z}{\partial t} + u_d\frac{\partial z}{\partial x} + v_d\frac{\partial z}{\partial y}$*/
		Add(BottomVerticalVelocity + face*BotBENfp, TempBottomVertVelocityX + face*BotBENfp, TempBottomVertVelocityY + face*BotBENfp, BotBENfp);
	}

	double *EleBotVertVelocity = malloc(Np*K2d*sizeof(double));

	memset(EleBotVertVelocity, 0, Np*K2d*sizeof(double));

	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){

		AssembleDataIntoPoint(EleBotVertVelocity + k*Np, BottomVerticalVelocity + k*Np2d, BotEidM, Np2d);

		dgemm("N", "N", &np, &oneI, &np, &one, RHSCoeMatrix + k*NLayer*Np*Np + (NLayer - 1)*Np*Np, \
			&np, RHS + k*NLayer*Np + (NLayer - 1)*Np, &np, &zero, VerticalVelocity + k*NLayer*Np + (NLayer - 1)*Np, \
			&np);
		dgemm("N", "N", &np, &oneI, &np, &one, VertCoeMatrix + k*NLayer*Np*Np + (NLayer - 1)*Np*Np, \
			&np, EleBotVertVelocity + k*Np, &np, &zero, VSTempVerticalVelocity + k*NLayer*Np + (NLayer - 1)*Np, \
			&np);
		Add(VerticalVelocity + k*NLayer*Np + (NLayer - 1)*Np, VerticalVelocity + k*NLayer*Np + (NLayer - 1)*Np, \
			VSTempVerticalVelocity + k*NLayer*Np + (NLayer - 1)*Np, Np);

		for (int L = 1; L < NLayer; L++){
			FetchBoundaryData(EleBotVertVelocity + k*Np, VerticalVelocity + k*NLayer*Np + (NLayer - L)*Np, BotEidM, UpEidM, Np2d);
			dgemm("N", "N", &np, &oneI, &np, &one, RHSCoeMatrix + k*NLayer*Np*Np + (NLayer - L - 1)*Np*Np, \
				&np, RHS + k*NLayer*Np + (NLayer - L - 1)*Np, &np, &zero, VerticalVelocity + k*NLayer*Np + (NLayer - L - 1)*Np, \
				&np);
			dgemm("N", "N", &np, &oneI, &np, &one, VertCoeMatrix + k*NLayer*Np*Np + (NLayer - L - 1)*Np*Np, \
				&np, EleBotVertVelocity + k*Np, &np, &zero, VSTempVerticalVelocity + k*NLayer*Np + (NLayer - L - 1)*Np, \
				&np);
			Add(VerticalVelocity + k*NLayer*Np + (NLayer - L - 1)*Np, VerticalVelocity + k*NLayer*Np + (NLayer - L - 1)*Np, \
				VSTempVerticalVelocity + k*NLayer*Np + (NLayer - L - 1)*Np, Np);

		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int p = 0; p < Np; p++){
			VerticalVelocity[k*Np + p] = H[k*Np + p] * VerticalVelocity[k*Np + p];
		}
	}

	free(BotEidM);
	free(UpEidM);
	free(VSTempVerticalVelocity);
	free(RHS);
	free(BottomVerticalVelocity);
	free(TempBottomVertVelocityX);
	free(TempBottomVertVelocityY);
	free(ubot);
	free(vbot);
	free(Hbot);
	free(zx);
	free(zy);
	free(EleBotVertVelocity);
}