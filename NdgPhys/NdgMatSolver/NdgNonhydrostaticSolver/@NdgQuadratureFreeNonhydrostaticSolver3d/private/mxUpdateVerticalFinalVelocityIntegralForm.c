
#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE.h"
#include "../../../../../NdgMath/NdgSWE3D.h"
#include "../../../../../NdgMath/NdgMemory.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *fphys = mxGetPr(prhs[0]);
	double *varIndex = mxGetPr(prhs[1]);
	double *Hcrit = mxGetPr(prhs[2]);
	double *PSPX = mxGetPr(prhs[3]);
	double *PSPY = mxGetPr(prhs[4]);
	double *PUPX = mxGetPr(prhs[5]);
	double *PVPY = mxGetPr(prhs[6]);
	double *PUPS = mxGetPr(prhs[7]);
	double *PVPS = mxGetPr(prhs[8]);
	const mxArray *mesh = prhs[9];
	const mxArray *cell = prhs[10];
	const mxArray *BottomBoundaryEdge = prhs[11];

	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);
	mxArray *TempJz = mxGetField(mesh, 0, "Jz");
	double *Jz = mxGetPr(TempJz);
	mxArray *TempNLayer = mxGetField(mesh, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);

	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempNpz = mxGetField(cell, 0, "Npz");
	int Npz = (int)mxGetScalar(TempNpz);
	mxArray *TempVint = mxGetField(cell, 0, "Vint");
	double *Vint = mxGetPr(TempVint);
	mxArray *TempV3d = mxGetField(cell, 0, "V");
	double *V3d = mxGetPr(TempV3d);
	double *InvV3d = malloc(Np*Np * sizeof(double));
	memcpy(InvV3d, V3d, Np*Np * sizeof(double));
	MatrixInverse(InvV3d, (ptrdiff_t)Np);
	mxArray *TempNz = mxGetField(cell, 0, "Nz");
	int Nz = (int)mxGetScalar(TempNz);
	mxArray *TempFmask = mxGetField(cell, 0, "Fmask");
	int maxNfp = (int)mxGetM(TempFmask);
	int Nface = (int)mxGetN(TempFmask);
	double *Fmask = mxGetPr(TempFmask);

	mxArray *TempBotBENfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
	int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	int Np2d = BotBENfp;
	mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBEFToE = mxGetField(BottomBoundaryEdge, 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToN1 = mxGetField(BottomBoundaryEdge, 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);
	mxArray *TempNe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int K2d = (int)mxGetScalar(TempNe);

	double *BotEidM = malloc(Np2d * sizeof(double));
	double *UpEidM = malloc(Np2d * sizeof(double));

	for (int i = 0; i < Np2d; i++) {
		BotEidM[i] = Fmask[(Nface - 2)*maxNfp + i];
		UpEidM[i] = Fmask[(Nface - 1)*maxNfp + i];
	}

	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *VerticalVelocity = mxGetPr(plhs[0]);

	double *ubot = malloc(Np2d * K2d * sizeof(double));
	double *vbot = malloc(Np2d*K2d * sizeof(double));
	double *wbot = malloc(Np2d*K2d * sizeof(double));
	double *Hbot = malloc(Np2d*K2d * sizeof(double));
	double *zx = malloc(Np2d*K2d * sizeof(double));
	double *zy = malloc(Np2d*K2d * sizeof(double));

	double *Hu = fphys + ((int)varIndex[0] - 1)*Np*K;
	double *Hv = fphys + ((int)varIndex[1] - 1)*Np*K;
	double *H = fphys + ((int)varIndex[3] - 1)*Np*K;
	double *zx3d = fphys + ((int)varIndex[4] - 1)*Np*K;
	double *zy3d = fphys + ((int)varIndex[5] - 1)*Np*K;
	double *RHS = malloc(Np*K * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int i = 0; i < Np; i++) {
			RHS[k*Np + i] = -1.0 * H[k*Np + i] * \
				(PUPX[k*Np + i] + PVPY[k*Np + i] + PSPX[k*Np + i] * PUPS[k*Np + i] + \
					PSPY[k*Np + i] * PVPS[k*Np + i]);
		}
	}

	double *fmod = malloc(K2d*Np * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		VerticalIntegralFromBottom(VerticalVelocity + k*NLayer*Np, RHS + k*NLayer*Np, Jz + k*NLayer*Np, fmod + k*Np, NLayer, (ptrdiff_t)Np, InvV3d, Np2d, Npz, Vint);
	}

	double *TempBottomVertVelocityX = malloc(Np2d*K2d * sizeof(double));
	double *TempBottomVertVelocityY = malloc(Np2d*K2d * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++) {
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
		Add(wbot + face*BotBENfp, TempBottomVertVelocityX + face*BotBENfp, TempBottomVertVelocityY + face*BotBENfp, BotBENfp);
	}
	double *wbot3d = malloc(Np*K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < BotBENe; k++) {
		NdgExtend2dField(wbot3d, wbot, BotBENfp, k, Np, NLayer, Nz);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int p = 0; p < Np; p++) {
			VerticalVelocity[k*Np + p] += wbot3d[k*Np+p];
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int p = 0; p < Np; p++) {
			VerticalVelocity[k*Np + p] = VerticalVelocity[k*Np + p] * H[k*Np + p];
		}
	}

	free(InvV3d);
	free(BotEidM);
	free(UpEidM);
	free(RHS);
	free(ubot);
	free(vbot);
	free(wbot);
	free(Hbot);
	free(zx);
	free(zy);
	free(fmod);
	free(TempBottomVertVelocityX);
	free(TempBottomVertVelocityY);
	free(wbot3d);
}