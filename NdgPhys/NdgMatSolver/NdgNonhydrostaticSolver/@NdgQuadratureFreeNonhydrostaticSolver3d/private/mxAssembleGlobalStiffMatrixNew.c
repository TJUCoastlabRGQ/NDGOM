#include "SWENonhydrostatic3d.h"
#include "stdio.h"

extern double *K33, *InvSquaHeight, *InnerEdgeTau, *BottomEdgeTau, *SurfaceEdgeTau;
extern char *GlobalStiffMatrixInitialized;

void MyExit()
{
	if (!strcmp("True", GlobalStiffMatrixInitialized)){
		GlobalStiffMatrixMemoryDeAllocation();
		GlobalStiffMatrixInitialized = "False";
	}
	return;
}

void GetInverseSquareHeight(double *, double *, double, int); 

void GetLocalVolumnIntegralTerm(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, double *M3d, \
	double *Dr, double *Ds, double *Dt, double *rx, double *sx, double *ry, double *sy, double *tz, double *J, double *K13, \
	double *K23, double *K33, int LocalEle);

void GetLocalToDownFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int DownEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *DownEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33, double Tau);

void GetLocalToUpFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int UpEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *UpEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33, double Tau);


void GetLocalDownFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33, double Tau);

void GetLocalUpFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33, double Tau);

void ImposeDirichletBoundaryCondition(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, \
	double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, \
	double *K23, double *K33, double Tau, int LocalEle);

void GetLocalToAdjacentFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *M2d, double *J2d, double *LocalEid, double *AdjEid, double *Dt, double *LocalTz, \
	double *AdjTz, double *nx, double *ny, double *LocalK13, double *LocalK23, \
	double *AdjK13, double *AdjK23, double Tau, int LocalEle, int AdjEle);

void GetLocalFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *M2d, double *J2d, double *LocalEid, double *Dt, double *tz, double *nx, \
	double *ny, double *K13, double *K23, double Tau, int LocalEle);
/*This function is used to assemble the global stiff matrix used in the 
three-dimensional non-hydrostatic solver, the input parameters are organized as follows:
PNPS: $\frac{\partial q}{\partial \sigma}$, indexed as 0;
SPNPX: $\frac{\partial^2 q}{\partial x^2}$, indexed as 1;
SPNPY: $\frac{\partial^2 q}{\partial y^2}$, indexed as 2;
SPNPS: $\frac{\partial^2 q}{\partial \sigma^2}$, indexed as 3;
MSPNPX: $\frac{\partial}{\partial x}\left (\frac{\partial q}{\partial \sigma}\right )$, indexed as 4;
MSPNPY: $\frac{\partial}{\partial y}\left (\frac{\partial q}{\partial \sigma}\right )$, indexed as 5;
PSPX: $\frac{\partial \sigma}{\partial x}$, indexed as 6;
PSPY: $\frac{\partial \sigma}{\partial y}$, indexed as 7;
SPSPX: $\frac{\partial^2 \sigma}{\partial x^2}$, indexed as 8;
SPSPY: $\frac{\partial^2 \sigma}{\partial y^2}$, indexed as 9;
SQPSPX: $\left (\frac{\partial \sigma}{\partial x}\right )^2$, indexed as 10;
SQPSPY: $\left (\frac{\partial \sigma}{\partial y}\right )^2$, indexed as 11;
Height: The water depth, indexed as 12;
Hcrit: The critical water depth, indexed as 13;
BoundaryEdge2d: The two-dimensional boundary edge, indexed as 14;
mesh: The three-dimensional mesh object, indexed as 15;
cell: The three-dimensional master cell, indexed as 16;
BoundaryEdge: The three-dimensional mesh object, indexed as 17;
ftype2d: The two-dimensional boundary edge type, indexed as 18;
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexAtExit(&MyExit);
	double *SPNPX = mxGetPr(prhs[0]);
	mwIndex *jcSPNPX = mxGetJc(prhs[0]);
	mwIndex *irSPNPX = mxGetIr(prhs[0]);

	double *SPNPY = mxGetPr(prhs[1]);
	mwIndex *jcSPNPY = mxGetJc(prhs[1]);
	mwIndex *irSPNPY = mxGetIr(prhs[1]);

	/*$\frac{\partial \sigma}{\partial x^*}$*/
	double *K13 = mxGetPr(prhs[2]);
	/*$\frac{\partial \sigma}{\partial y^*}$*/
	double *K23 = mxGetPr(prhs[3]);

	double *SQPSPX = mxGetPr(prhs[4]);
	double *SQPSPY = mxGetPr(prhs[5]);

	double Hcrit = mxGetScalar(prhs[6]);
	double *Height = mxGetPr(prhs[7]);

	const mxArray *mesh = prhs[8];

	mxArray *Temprx = mxGetField(mesh, 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(mesh, 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(mesh, 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *TempLAV = mxGetField(mesh, 0, "LAV");
	double *LAV = mxGetPr(TempLAV);
	mxArray *Tempsy = mxGetField(mesh, 0, "sy");
	double *sy = mxGetPr(Tempsy);
	mxArray *Temptz = mxGetField(mesh, 0, "tz");
	double *tz = mxGetPr(Temptz);
	mxArray *TempJ = mxGetField(mesh, 0, "J");
	double *J = mxGetPr(TempJ);
	mxArray *TempEToE = mxGetField(mesh, 0, "EToE");
	double *EToE = mxGetPr(TempEToE);
	mxArray *TempNlayer = mxGetField(mesh, 0, "Nz");
	int Nlayer = (int)mxGetScalar(TempNlayer);
	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);


	const mxArray *cell = prhs[9];
	mxArray *TempFmask = mxGetField(cell, 0, "Fmask");
	double *Fmask = mxGetPr(TempFmask);
	int maxNfp = (int)mxGetM(TempFmask);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	int Nface2d = Nface - 2;
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempN = mxGetField(cell, 0, "N");
	int N = (int)mxGetScalar(TempN);
	mxArray *TempNz = mxGetField(cell, 0, "Nz");
	int Nz = (int)mxGetScalar(TempNz);
	mxArray *TempMass3d = mxGetField(cell, 0, "M");
	double  *M3d = mxGetPr(TempMass3d);
	mxArray *TempDr = mxGetField(cell, 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(cell, 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempNfp = mxGetField(cell, 0, "Nfp");
	double *Nfp = mxGetPr(TempNfp);
	int Np2d = (int)Nfp[Nface - 1];
	double *UpEidM = malloc(Np2d*sizeof(double));
	double *BotEidM = malloc(Np2d*sizeof(double));
	for (int i = 0; i < Np2d; i++){
		UpEidM[i] = Fmask[(Nface - 1)*maxNfp + i];
		BotEidM[i] = Fmask[(Nface - 2)*maxNfp + i];
	}

	const mxArray *InnerEdge = prhs[10];
	mxArray *TempIENe = mxGetField(InnerEdge, 0, "Ne");
	int IENe = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp = mxGetField(InnerEdge, 0, "Nfp");
	int IENfp = (int)mxGetScalar(TempIENfp);
	mxArray *TempIEMb = mxGetField(InnerEdge, 0, "M");
	double *IEMb = mxGetPr(TempIEMb);
	mxArray *TempIEJs = mxGetField(InnerEdge, 0, "Js");
	double *IEJs = mxGetPr(TempIEJs);
	mxArray *TempIEnx = mxGetField(InnerEdge, 0, "nx");
	double *IEnx = mxGetPr(TempIEnx);
	mxArray *TempIEny = mxGetField(InnerEdge, 0, "ny");
	double *IEny = mxGetPr(TempIEny);
	mxArray *TempIELAV = mxGetField(InnerEdge, 0, "LAV");
	double *IELAV = mxGetPr(TempIELAV);
	mxArray *TempIEFToE = mxGetField(InnerEdge, 0, "FToE");
	double *IEFToE = mxGetPr(TempIEFToE);
	mxArray *TempIEFToF = mxGetField(InnerEdge, 0, "FToF");
	double *IEFToF = mxGetPr(TempIEFToF);
	mxArray *TempIEFToN1 = mxGetField(InnerEdge, 0, "FToN1");
	double *IEFToN1 = mxGetPr(TempIEFToN1);
	mxArray *TempIEFToN2 = mxGetField(InnerEdge, 0, "FToN2");
	double *IEFToN2 = mxGetPr(TempIEFToN2);

	const mxArray *BottomEdge = prhs[11];
	mxArray *TempBotENe = mxGetField(BottomEdge, 0, "Ne");
	int BotENe = (int)mxGetScalar(TempBotENe);
	mxArray *TempBotENfp = mxGetField(BottomEdge, 0, "Nfp");
	int BotENfp = (int)mxGetScalar(TempBotENfp);
	mxArray *TempBotEFToE = mxGetField(BottomEdge, 0, "FToE");
	double *BotEFToE = mxGetPr(TempBotEFToE);
	mxArray *TempBotEFToN1 = mxGetField(BottomEdge, 0, "FToN1");
	double *BotEFToN1 = mxGetPr(TempBotEFToN1);
	mxArray *TempBotEFToN2 = mxGetField(BottomEdge, 0, "FToN2");
	double *BotEFToN2 = mxGetPr(TempBotEFToN2);
	mxArray *TempBotELAV = mxGetField(BottomEdge, 0, "LAV");
	double *BotELAV = mxGetPr(TempBotELAV);

	const mxArray *SurfaceBoundaryEdge = prhs[12];
	mxArray *TempSurfBENe = mxGetField(SurfaceBoundaryEdge, 0, "Ne");
	int SurfBENe = (int)mxGetScalar(TempSurfBENe);
	mxArray *TempSurfBENfp = mxGetField(SurfaceBoundaryEdge, 0, "Nfp");
	int SurfBENfp = (int)mxGetScalar(TempSurfBENfp);
	mxArray *TempSurfBEFToE = mxGetField(SurfaceBoundaryEdge, 0, "FToE");
	double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
	mxArray *TempSurfBEFToN1 = mxGetField(SurfaceBoundaryEdge, 0, "FToN1");
	double *SurfBEFToN1 = mxGetPr(TempSurfBEFToN1);
	mxArray *TempSurfBEFToN2 = mxGetField(SurfaceBoundaryEdge, 0, "FToN2");
	double *SurfBEFToN2 = mxGetPr(TempSurfBEFToN2);
	mxArray *TempSurfBELAV = mxGetField(SurfaceBoundaryEdge, 0, "LAV");
	double *SurfBELAV = mxGetPr(TempSurfBELAV);


	double *M2d = mxGetPr(prhs[13]);
	double *J2d = mxGetPr(prhs[14]);
	int K2d = (int)mxGetScalar(prhs[15]);

	char* BoundaryType;
	BoundaryType = mxArrayToString(prhs[16]);

	if (!strcmp("False", GlobalStiffMatrixInitialized)){
		GlobalStiffMatrixMemoryAllocation(Np, K, IENe, BotENe, SurfBENe);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		GetInverseSquareHeight(InvSquaHeight + Np*k, Height + Np*k, Hcrit, Np);
		for (int p = 0; p < Np; p++)
			K33[k*Np + p] = SQPSPX[k*Np + p] + SQPSPY[k*Np + p] + \
			InvSquaHeight[k*Np + p];
	}

	int row, col;
	double *sr;
	mwIndex *irs, *jcs;

	row = (int)mxGetM(prhs[0]);
	col = (int)mxGetN(prhs[0]);

	plhs[0] = mxCreateSparse(row, col, jcSPNPX[col], mxREAL);
	sr = mxGetPr(plhs[0]);
	irs = mxGetIr(plhs[0]);
	jcs = mxGetJc(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int i = 0; i < col + 1; i++){
		jcs[i] = jcSPNPX[i];
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int i = 0; i < (int)jcSPNPX[col]; i++){
		irs[i] = irSPNPX[i];
		sr[i] = SPNPX[i] + SPNPY[i];
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int i = 0; i < IENe; i++){
		CalculatePenaltyParameter(InnerEdgeTau, IEFToE, IEFToN1, IEFToN2, Np, IENfp, \
			i, K13, K23, K33, IELAV, LAV, max(N, Nz), Nface);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int i = 0; i < BotENe; i++){
		CalculatePenaltyParameter(BottomEdgeTau, BotEFToE, BotEFToN1, BotEFToN2, Np, BotENfp, \
			i, K13, K23, K33, BotELAV, LAV, max(N, Nz), Nface);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int i = 0; i < SurfBENe; i++){
		CalculatePenaltyParameter(SurfaceEdgeTau, SurfBEFToE, SurfBEFToN1, SurfBEFToN2, Np, SurfBENfp, \
			i, K13, K23, K33, SurfBELAV, LAV, max(N, Nz), Nface);
	}
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int ele = 0; ele < K; ele++){
		int EleNumber = 0;
		double *TempEToE = malloc((Nface + 1)*sizeof(double));
		FindUniqueElementAndSortOrder(TempEToE, EToE + ele*Nface, &EleNumber, Nface, ele + 1);

		for (int i = 0; i < EleNumber; i++){
			if (ele + 1 == (int)TempEToE[i]){
				LocalStartPoint = jcs[ele*Np] + i*Np;
				GetLocalVolumnIntegralTerm(sr, irs, jcs, \
					Np, M3d,Dr, Ds, Dt, rx + ele*Np, sx + ele*Np, ry + ele*Np, \
					sy + ele*Np, tz + ele*Np, J + ele*Np, K13 + ele*Np, \
					K23 + ele*Np, K33 + ele*Np, ele + 1);
				break;
			}
		}

		double *Facialnx = malloc(Nface2d*IENfp*sizeof(double));
		double *Facialny = malloc(Nface2d*IENfp*sizeof(double));
		int *GlobalFace = malloc(Nface2d*sizeof(int));
		int *AdjEle = malloc(Nface2d*sizeof(int));
		int *ReverseFlag = malloc(Nface2d*sizeof(int));
		int InternalFace = 0;

		FindFaceAndDirectionVector(Facialnx, GlobalFace, AdjEle, \
			&InternalFace, ReverseFlag, IENfp, ele + 1, IEFToE, IEFToF, IEnx, IENe, Nface2d);

		FindFaceAndDirectionVector(Facialny, GlobalFace, AdjEle, \
			&InternalFace, ReverseFlag, IENfp, ele + 1, IEFToE, IEFToF, IEny, IENe, Nface2d);

		double *LocalEidM = malloc(IENfp*sizeof(double));

		double *AdjEidM = malloc(IENfp*sizeof(double));

		for (int i = 0; i < InternalFace; i++){
			if (ReverseFlag[i] == 0){
				for (int p = 0; p < IENfp; p++){
					LocalEidM[p] = IEFToN1[GlobalFace[i] * IENfp + p];
					AdjEidM[p] = IEFToN2[GlobalFace[i] * IENfp + p];
				}
			}
			else{
				for (int p = 0; p < IENfp; p++){
					LocalEidM[p] = IEFToN2[GlobalFace[i] * IENfp + p];
					AdjEidM[p] = IEFToN1[GlobalFace[i] * IENfp + p];
				}
			}

			GetLocalFacialContributionInHorizontalDirection(sr, irs, jcs, Np, IENfp, \
				IEMb, IEJs + GlobalFace[i] * IENfp, LocalEidM, Dt, tz + ele*Np, Facialnx + i * IENfp, \
				Facialny + i*IENfp, K13 + ele*Np, K23 + ele*Np, *(InnerEdgeTau + GlobalFace[i]), ele + 1);

			for (int j = 0; j < EleNumber; j++){
				if ((int)TempEToE[j] == AdjEle[i]){
					GetLocalToAdjacentFacialContributionInHorizontalDirection(sr, irs, jcs, Np, IENfp,\
						IEMb, IEJs + GlobalFace[i] * IENfp, LocalEidM, AdjEidM, Dt, tz + ele*Np, \
						tz + (int)(TempEToE[j] - 1)*Np, Facialnx + i*IENfp, Facialny + i*IENfp, K13 + ele*Np, K23 + ele*Np, \
						K13 + (int)(TempEToE[j] - 1)*Np, K23 + (int)(TempEToE[j] - 1)*Np, *(InnerEdgeTau + GlobalFace[i]), ele + 1, AdjEle[i]);

				}
			}
		}
		free(TempEToE);
		free(Facialnx);
		free(Facialny);
		free(GlobalFace);
		free(AdjEle);
		free(ReverseFlag);
		free(LocalEidM);
		free(AdjEidM);
	}
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int ele = 0; ele < K2d; ele++){
		int StartPoint;
		for (int L = 0; L < Nlayer; L++){
			//Index of the studied element
			int LocalEle = ele*Nlayer + L + 1;
			//Index of the element that located upside of the studied element 
			int UpEle = (int)EToE[ele*Nlayer*Nface + L*Nface + Nface - 1];
			//Index of the element that located downside of the studied element
			int DownEle = (int)EToE[ele*Nlayer*Nface + L*Nface + Nface - 2];

			int EleNumber = 0;

			double *TempEToE = malloc((Nface + 1)*sizeof(double));

			FindUniqueElementAndSortOrder(TempEToE, EToE + (LocalEle - 1)*Nface, &EleNumber, Nface, LocalEle);

			if (UpEle == DownEle){//Only one layer in the vertical direction, impose the Dirichlet boundary condition
				for (int i = 0; i < EleNumber; i++){
					if (LocalEle == (int)TempEToE[i]){
						if (!strcmp(BoundaryType, "Dirichlet")){
							ImposeDirichletBoundaryCondition(sr, irs, jcs, Np, Np2d, M2d, \
								J2d + ele*Np2d, UpEidM, rx + (LocalEle - 1)*Np, \
								sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
								tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np,\
								K33 + (LocalEle - 1)*Np, *(SurfaceEdgeTau + ele), LocalEle);
						}
						break;
					}
				}
			}
			else{//two or more layers included in vertical direction
				if (LocalEle == UpEle){//This is the top most cell, and L=0
					for (int i = 0; i < EleNumber; i++){
						if (LocalEle == (int)TempEToE[i]){
							StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
							if (!strcmp(BoundaryType, "Dirichlet")){
								ImposeDirichletBoundaryCondition(sr, irs, jcs, Np, Np2d, M2d, \
									J2d + ele*Np2d, UpEidM, rx + (LocalEle - 1)*Np, \
									sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
									tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, \
									K33 + (LocalEle - 1)*Np, *(SurfaceEdgeTau + ele), LocalEle);
							}
							break;
						}
					}
					
					GetLocalDownFacialContribution(sr, irs, jcs, LocalEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, BotEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np,\
						K33 + (LocalEle - 1)*Np, *(BottomEdgeTau + L*K2d + ele));
					
					GetLocalToDownFacialContribution(sr, irs, jcs, LocalEle, DownEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, BotEidM, UpEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (DownEle - 1)*Np, sx + (DownEle - 1)*Np, \
						ry + (DownEle - 1)*Np, sy + (DownEle - 1)*Np, tz + (DownEle - 1)*Np, Dr, Ds, Dt,\
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np,\
						K13 + (DownEle - 1)*Np, K23 + (DownEle - 1)*Np, K33 + (DownEle - 1)*Np, *(BottomEdgeTau + L*K2d + ele));
						

				}
				else if (LocalEle == DownEle){// This is the bottom most cell.
					

					GetLocalUpFacialContribution(sr, irs, jcs, LocalEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np,\
						K33 + (LocalEle - 1)*Np, *(BottomEdgeTau + (L-1)*K2d + ele));
					
					GetLocalToUpFacialContribution(sr, irs, jcs, LocalEle, UpEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, BotEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (UpEle - 1)*Np, sx + (UpEle - 1)*Np, \
						ry + (UpEle - 1)*Np, sy + (UpEle - 1)*Np, tz + (UpEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (UpEle - 1)*Np, K23 + (UpEle - 1)*Np, K33 + (UpEle - 1)*Np, *(BottomEdgeTau + (L - 1)*K2d + ele));
						
				}
				else{
					
					GetLocalUpFacialContribution(sr, irs, jcs, LocalEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np,\
						K33 + (LocalEle - 1)*Np, *(BottomEdgeTau + (L - 1)*K2d + ele));
					
					GetLocalToUpFacialContribution(sr, irs, jcs, LocalEle, UpEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, UpEidM, BotEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (UpEle - 1)*Np, sx + (UpEle - 1)*Np, \
						ry + (UpEle - 1)*Np, sy + (UpEle - 1)*Np, tz + (UpEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (UpEle - 1)*Np, K23 + (UpEle - 1)*Np, K33 + (UpEle - 1)*Np, *(BottomEdgeTau + (L - 1)*K2d + ele));
                     
					GetLocalDownFacialContribution(sr, irs, jcs, LocalEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, BotEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, \
						K33 + (LocalEle - 1)*Np, *(BottomEdgeTau + L*K2d + ele));
                      
					GetLocalToDownFacialContribution(sr, irs, jcs, LocalEle, DownEle, Np, Np2d, M2d, \
						J2d + ele*Np2d, BotEidM, UpEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (DownEle - 1)*Np, sx + (DownEle - 1)*Np, \
						ry + (DownEle - 1)*Np, sy + (DownEle - 1)*Np, tz + (DownEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (DownEle - 1)*Np, K23 + (DownEle - 1)*Np, K33 + (DownEle - 1)*Np, *(BottomEdgeTau + L*K2d + ele));
						

				}
			}
			free(TempEToE);
		}
	}
	free(UpEidM);
	free(BotEidM);
}

void GetLocalVolumnIntegralTerm(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, double *M3d, \
	double *Dr, double *Ds, double *Dt, double *rx, double *sx, double *ry, double *sy, double *tz, double *J, double *K13, \
	double *K23, double *K33, int LocalEle){
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));
	double *VertDiffMatrix = malloc(Np*Np*sizeof(double));
	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, M3d, J, Np);
	double *LumpedMass3d = malloc(Np*Np*sizeof(double));
	MassLumping(LumpedMass3d, Mass3d, Np);
	double *TempContributionBuff = malloc(Np*Np*sizeof(double));
	double *ContributionBuff = malloc(Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	memset(Contribution, 0, Np*Np*sizeof(double));

	/*For term $$k_{13}\frac{\partial v}{\partial \sigma}\frac{\partial p_h}{\partial x}$$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);
	DiagMultiply(VertDiffMatrix, VertDiffMatrix, K13, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, VertDiffMatrix,
		(ptrdiff_t)Np, LumpedMass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{23}\frac{\partial v}{\partial \sigma}\frac{\partial p_h}{\partial y}$$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);
	DiagMultiply(VertDiffMatrix, VertDiffMatrix, K23, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, VertDiffMatrix,
		(ptrdiff_t)Np, LumpedMass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{31}\frac{\partial v}{\partial x}\frac{\partial p_h}{\partial \sigma}$$.Here  $k_{31}=k_{13}$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, LumpedMass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VertDiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{32}\frac{\partial v}{\partial y}\frac{\partial p_h}{\partial \sigma}$$. Here, $k_{32}=k_{23}$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, LumpedMass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VertDiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{33}\frac{\partial v}{\partial \sigma}\frac{\partial p_h}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, LumpedMass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VertDiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	AssembleVolumnContributionIntoSparseMatrix(dest, Ir, Jc, Np, Contribution, LocalEle);

	free(DiffMatrix);
	free(TempDiffMatrix);
	free(VertDiffMatrix);
	free(Mass3d);
	free(LumpedMass3d);
	free(TempContributionBuff);
	free(ContributionBuff);
	free(Contribution);
}

void GetLocalToDownFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int DownEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *DownEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33, double Tau){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p}{\partial x}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial p}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p}{\partial y}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial p}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Adjtz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Localtz, Np);
	/*For term $k_{33}\frac{\partial p}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, DownEid, LocalEid, Np, Np2d, 1.0);

	double *SortedDownEid = malloc(Np2d*sizeof(double));
	memcpy(SortedDownEid, DownEid, Np2d);
	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Np2d);
	Sort(SortedDownEid, Np2d);
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedLocalEid, SortedDownEid, Np, Np2d, TempContribution, LocalEle, DownEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempMass2d);
	free(SortedDownEid);
	free(SortedLocalEid);
}

void GetLocalToUpFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int UpEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *UpEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33, double Tau){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p}{\partial x}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial p}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p}{\partial y}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial p}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1 * 0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Adjtz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Localtz, Np);
	/*For term $k_{33}\frac{\partial p}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1 * 0.5, Np*Np2d);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, UpEid, LocalEid, Np, Np2d, 1.0);

	double *SortedUpEid = malloc(Np2d*sizeof(double));
	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedUpEid, UpEid, Np2d);
	memcpy(SortedLocalEid, LocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedLocalEid, SortedUpEid, Np, Np2d, TempContribution, LocalEle, UpEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempMass2d);
	free(SortedUpEid);
	free(SortedLocalEid);
}


void GetLocalDownFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33, double Tau){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Np2d*sizeof(double));
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedLocalEid, SortedLocalEid, Np, Np2d, TempContribution, LocalEle, LocalEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempMass2d);
	free(SortedLocalEid);
}

void GetLocalUpFacialContribution(double *dest, mwIndex *irs, mwIndex *jcs, int LocalEle, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33, double Tau){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *EdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, EdgeContribution, (ptrdiff_t)Np2d);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Np2d);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	double *SortedLocalEid = malloc(Np2d*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Np2d*sizeof(double));
	Sort(SortedLocalEid, Np2d);

	AssembleFacialContributionIntoSparseMatrix(dest, irs, jcs, SortedLocalEid, SortedLocalEid, Np, Np2d, TempContribution, LocalEle, LocalEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(TempMass2d);
	free(SortedLocalEid);
}

void ImposeDirichletBoundaryCondition(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Np2d, \
	double *M2d, double *J2d, double *LocalEid, double *rx, double *sx,\
	double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, \
	double *K23, double *K33, double Tau, int LocalEle){

	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Np2d*sizeof(double));

	double *BoundaryEdgeContribution = malloc(Np*Np2d*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Np2d, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np2d, (ptrdiff_t)Np2d, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Np2d, FacialMass2d, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np2d, (ptrdiff_t)Np, (ptrdiff_t)Np2d, 1.0, FacialMass2d,
		(ptrdiff_t)Np2d, FacialDiffMatrix, (ptrdiff_t)Np2d, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np2d);

	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Np2d);

	double *TempMass2d = malloc(Np2d*Np2d*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Np2d*Np2d);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Np2d, -1.0);

	double *SortedLocalEid = malloc(Nfp*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Nfp*sizeof(double));
	Sort(SortedLocalEid, Nfp);

	AssembleFacialContributionIntoSparseMatrix(dest, Ir, Jc, SortedLocalEid, SortedLocalEid, Np, Nfp, TempContribution, LocalEle, LocalEle);

	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(BoundaryEdgeContribution);
	free(TempMass2d);
	free(SortedLocalEid);
}


void GetLocalToAdjacentFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *M2d, double *J2d, double *LocalEid, double *AdjEid, double *Dt, double *LocalTz, \
	double *AdjTz, double *nx, double *ny, double *LocalK13, double *LocalK23, \
	double *AdjK13, double *AdjK23, double Tau, int LocalEle, int AdjEle){
	// \D5\E2\C0\EF\B5\C4nx \BA\CDny \D3ɱ\BE\B5ص\A5Ԫָ\B3\F6
	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *AdjDiff = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *InnerEdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{13}\frac{\partial v}{\partial \sigma}\right\}[p_h]_x $$*/
	DiagMultiply(AdjDiff, Dt, AdjTz, Np);
	/*For term $k_{13}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(AdjDiff, AdjDiff, AdjK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, AdjDiff, AdjEid, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*nx[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial v}{\partial \sigma}\right\}[p_h]_y $$*/
	DiagMultiply(AdjDiff, Dt, AdjTz, Np);
	/*For term $k_{23}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(AdjDiff, AdjDiff, AdjK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, AdjDiff, AdjEid, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*ny[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$k_{13}\frac{\partial p_h}{\partial \sigma} [v]_x$$*/
	double *LocalDiff = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiff, Dt, LocalTz, Np);
	/*For term $k_{13}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5*nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEid, Np, Nfp);

	/*For term $$k_{23}\frac{\partial p_h}{\partial \sigma} [v]_y$$*/
	DiagMultiply(LocalDiff, Dt, LocalTz, Np);
	/*For term $k_{23}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5*ny[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEid, Np, Nfp);

	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Nfp*Nfp);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, AdjEid, LocalEid, Np, Nfp, 1.0);

	double *SortedAdjEid = malloc(Nfp*sizeof(double));
	memcpy(SortedAdjEid, AdjEid, Nfp);
	double *SortedLocalEid = malloc(Nfp*sizeof(double));
	memcpy(SortedLocalEid, LocalEid, Nfp);
	Sort(SortedAdjEid, Nfp);
	Sort(SortedLocalEid, Nfp);

	AssembleFacialContributionIntoSparseMatrix(dest, Ir, Jc, SortedLocalEid, SortedAdjEid, Np, Nfp, TempContribution, LocalEle, AdjEle);

	free(FacialMass2d);
	free(TempContribution);
	free(AdjDiff);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(LocalDiff);
	free(TempMass2d);
	free(SortedLocalEid);
	free(SortedAdjEid);
}

void GetLocalFacialContributionInHorizontalDirection(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, int Nfp, \
	double *M2d, double *J2d, double *LocalEid, double *Dt, double *tz, double *nx, \
	double *ny, double *K13, double *K23, double Tau, int LocalEle){

	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));

	double *LocalDiff = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *InnerEdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{13}\frac{\partial v}{\partial \sigma}\right\}[p_h]_x $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{13}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*nx[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial v}{\partial \sigma}\right\}[p_h]_y $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{23}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*ny[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{13}\frac{\partial p_h}{\partial \sigma}\right\}[v]_x $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{13}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial p_h}{\partial \sigma}\right\}[v]_y $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{23}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*ny[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEid, Np, Nfp);

	/*The follow part is for term $ - \int_{ \epsilon_i }\tau^k[p_h][s]d\boldsymbol{ x }$*/
	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));
	MultiplyByConstant(TempMass2d, FacialMass2d, Tau, Nfp*Nfp);
	AssembleContributionIntoRowAndColumn(TempContribution, TempMass2d, LocalEid, LocalEid, Np, Nfp, -1.0);

	double *SortedLocalEidM = malloc(Nfp*sizeof(double));
	memcpy(SortedLocalEidM, LocalEid, Nfp*sizeof(double));
	Sort(SortedLocalEidM, Nfp);

	AssembleFacialContributionIntoSparseMatrix(dest, Ir, Jc, SortedLocalEidM, SortedLocalEidM, Np, Nfp, TempContribution, LocalEle, LocalEle);
	free(FacialMass2d);
	free(TempContribution);
	free(LocalDiff);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(TempMass2d);
	free(SortedLocalEidM);
}

void GetInverseSquareHeight(double *dest, double *source, double Hcrit, int Np){
	for (int i = 0; i < Np; i++){
		if (source[i] >= Hcrit){
			dest[i] = 1.0 / source[i] / source[i];
		}
		else{
			dest[i] = 0;
		}
	}
}
