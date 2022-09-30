#include "SWENonhydrostatic3d.h"

extern double *NonhydroUbot, *NonhydroVbot, *NonhydroHbot, \
*NonhydroZxbot, *NonhydroZybot;

extern char *NHVertVelocityInitialized;

void MyExit()
{
	if (!strcmp("True", NHVertVelocityInitialized)){
		SWENonhydroVertVelocityMemoryDeAllocation();
		NHVertVelocityInitialized = "False";
	}
	return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	const mxArray *cell = prhs[0];
	const mxArray *BottomBoundaryEdge = prhs[1];
	double *fphys = mxGetPr(prhs[2]);
	double *varIndex = mxGetPr(prhs[3]);
	const mxArray *mesh = prhs[4];
	double Hcrit = mxGetScalar(prhs[5]);

	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);

	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBENfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
	int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	mxArray *TempBotBEFToE = mxGetField(BottomBoundaryEdge, 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToN1 = mxGetField(BottomBoundaryEdge, 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);

	plhs[0] = mxCreateDoubleMatrix(BotBENfp, BotBENe, mxREAL);
	double *W = mxGetPr(plhs[0]);

	if (!strcmp("False", NHVertVelocityInitialized)){
		SWENonhydroVertVelocityMemoryAllocation(BotBENfp, BotBENe);
	}

	double *Hu = fphys + ((int)varIndex[0] - 1)*Np*K;
	double *Hv = fphys + ((int)varIndex[1] - 1)*Np*K;
	double *H = fphys + ((int)varIndex[3] - 1)*Np*K;
	double *Zx = fphys + ((int)varIndex[4] - 1)*Np*K;
	double *Zy = fphys + ((int)varIndex[5] - 1)*Np*K;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		/*The date is stored first in U is Hu*/
		FetchBoundaryEdgeFacialValue(NonhydroUbot + face*BotBENfp, Hu, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*The date is stored first in V is Hv*/
		FetchBoundaryEdgeFacialValue(NonhydroVbot + face*BotBENfp, Hv, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(NonhydroHbot + face*BotBENfp, H, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*Dot divide Hu by H to get u*/
		DotCriticalDivide(NonhydroUbot + face*BotBENfp, NonhydroUbot + face*BotBENfp, &Hcrit, NonhydroHbot + face*BotBENfp, BotBENfp);
		/*Dot divide Hv by H to get v*/
		DotCriticalDivide(NonhydroVbot + face*BotBENfp, NonhydroVbot + face*BotBENfp, &Hcrit, NonhydroHbot + face*BotBENfp, BotBENfp);
		/*$Tempz_x=\frac{\partial z}{\partial x}$*/
		FetchBoundaryEdgeFacialValue(NonhydroZxbot + face*BotBENfp, Zx, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*$u_dz_x$, at present, the data is stored in TempZx2d*/
		DotProduct(NonhydroZxbot + face*BotBENfp, NonhydroUbot + face*BotBENfp, NonhydroZxbot + face*BotBENfp, BotBENfp);
		FetchBoundaryEdgeFacialValue(NonhydroZybot + face*BotBENfp, Zy, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		/*$v_dz_y$, at present, the data is stored in TempZy2d*/
		DotProduct(NonhydroZybot + face*BotBENfp, NonhydroVbot + face*BotBENfp, NonhydroZybot + face*BotBENfp, BotBENfp);
		/*$w_{bot} = \frac{\partial z}{\partial t} + u_d\frac{\partial z}{\partial x} + v_d\frac{\partial z}{\partial y}$*/
		Add(W + face*BotBENfp, NonhydroZxbot + face*BotBENfp, NonhydroZybot + face*BotBENfp, BotBENfp);
	}

}