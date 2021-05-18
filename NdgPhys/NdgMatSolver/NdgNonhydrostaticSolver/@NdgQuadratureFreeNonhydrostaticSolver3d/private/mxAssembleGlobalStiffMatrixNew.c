#include "SWENonhydrostatic3d.h"

void GetInverseSquareHeight(double *, double *, double, int);

void GetLocalToDownFacialContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int DownEle, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, double *J, double *J2d, double *LocalEid, double *DownEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33);

void GetLocalToUpFacialContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int UpEle, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, double *J, double *J2d, double *LocalEid, double *UpEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33);

void GetLocalDownFacialContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, double *J, double *J2d, double *Eid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33);

void GetLocalUpFacialContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, double *J, double *J2d, double *Eid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33);

void GetLocalFacialContributionInHorizontalDirection(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *M3d, double *M2d, double *J, double *J2d, double *EidM, double *Dt, double *tz, double *nx, \
	double *ny, double *K13, double *K23);

void GetLocalToAdjacentFacialContributionInHorizontalDirection(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *M3d, double *M2d, double *AdjacentJ, double *J2d, double *LocalEid, double *AdjEid, double *Dt, double *LocalTz, \
	double *AdjTz, double *nx, double *ny, double *LocalK13, double *LocalK23, \
	double *AdjK13, double *AdjK23);

void ImposeDirichletBoundaryCondition(double *dest, int StartPoint, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, \
	double *J, double *J2d, double *LocalEid, double *rx, double *sx, \
	double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33);

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

	mxArray *mesh = prhs[8];

	mxArray *Temprx = mxGetField(mesh, 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(mesh, 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(mesh, 0, "ry");
	double *ry = mxGetPr(Tempry);
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


	mxArray *cell = prhs[9];
	mxArray *TempFmask = mxGetField(cell, 0, "Fmask");
	double *Fmask = mxGetPr(TempFmask);
	int maxNfp = (int)mxGetM(TempFmask);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempMass3d = mxGetField(cell, 0, "M");
	double  *Mass3d = mxGetPr(TempMass3d);
	mxArray *TempDr = mxGetField(cell, 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(cell, 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempNfp = mxGetField(cell, 0, "Nfp");
	double *Nfp = mxGetPr(TempNfp);
	int HorNfp = (int)Nfp[0];
	int VertNfp = (int)Nfp[Nface - 1];
	double *UpEidM = malloc(VertNfp*sizeof(double));
	double *BotEidM = malloc(VertNfp*sizeof(double));
	for (int i = 0; i < VertNfp; i++){
		UpEidM[i] = Fmask[(Nface - 1)*maxNfp + i];
		BotEidM[i] = Fmask[(Nface - 2)*maxNfp + i];
	}

	mxArray *InnerEdge = prhs[10];
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

	double *K33 = malloc(Np*K*sizeof(double));

	double *InvSquaHeight = malloc(Np*K*sizeof(double));

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
	for (int ele = 0; ele < Ele3d; ele++){
		int EleNumber = 0;
		int StartPoint;
		int LocalStartPoint;
		double *TempEToE = malloc((Nface + 1)*sizeof(double));
		FindUniqueElementAndSortOrder(TempEToE, EToE + ele*Nface, &EleNumber, Nface, ele + 1);

		for (int i = 0; i < EleNumber; i++){
			if (ele + 1 == (int)TempEToE[i]){
				LocalStartPoint = jcs[ele*Np] + i*Np;
				GetLocalVolumnIntegralTerm(sr, LocalStartPoint, \
					Np, jcs[ele*Np + 1] - jcs[ele*Np], M3d, \
					Dr, Ds, Dt, rx + ele*Np, sx + ele*Np, ry + ele*Np, \
					sy + ele*Np, tz + ele*Np, J + ele*Np, K13 + ele*Np, \
					K23 + ele*Np, K33 + ele*Np);
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
			&InternalFace, ReverseFlag, IENfp, ele + 1, FToE, FToF, IEnx, IENe, Nface2d);

		FindFaceAndDirectionVector(Facialny, GlobalFace, AdjEle, \
			&InternalFace, ReverseFlag, IENfp, ele + 1, FToE, FToF, IEny, IENe, Nface2d);

		double *LocalEidM = malloc(IENfp*sizeof(double));

		double *AdjEidM = malloc(IENfp*sizeof(double));

		for (int i = 0; i < InternalFace; i++){
			if (ReverseFlag[i] == 0){
				for (int p = 0; p < IENfp; p++){
					LocalEidM[p] = FToN1[GlobalFace[i] * IENfp + p];
					AdjEidM[p] = FToN2[GlobalFace[i] * IENfp + p];
				}
			}
			else{
				for (int p = 0; p < IENfp; p++){
					LocalEidM[p] = FToN2[GlobalFace[i] * IENfp + p];
					AdjEidM[p] = FToN1[GlobalFace[i] * IENfp + p];
				}
			}

			GetLocalFacialContributionInHorizontalDirection(sr, LocalStartPoint, Np, IENfp, jcs[ele*Np + 1] - jcs[ele*Np], \
				Mass3d, IEMb, J + ele*Np, Js + GlobalFace[i] * IENfp, LocalEidM, Dt, tz + ele*Np, Facialnx + i * IENfp, \
				Facialny + i*IENfp, K13 + ele*Np, K23 + ele*Np);

			for (int j = 0; j < EleNumber; j++){
				if ((int)TempEToE[j] == AdjEle[i]){
					StartPoint = jcs[ele*Np] + j*Np;
					GetLocalToAdjacentFacialContributionInHorizontalDirection(sr, StartPoint, Np, IENfp, jcs[ele*Np + 1] - jcs[ele*Np], \
						Mass3d, IEMb, J + (int)(TempEToE[j] - 1)*Np, Js + GlobalFace[i] * IENfp, LocalEidM, AdjEidM, Dt, tz + ele*Np, \
						tz + (int)(TempEToE[j] - 1)*Np, Facialnx + i*IENfp, Facialny + i*IENfp, K13 + ele*Np, K23 + ele*Np, \
						K13 + (int)(TempEToE[j] - 1)*Np, K23 + (int)(TempEToE[j] - 1)*Np);
				}
			}
		}
		free(TempEToE);
		free(FacialVector);
		free(GlobalFace);
		free(AdjEle);
		free(ReverseFlag);
		free(LocalEidM);
		free(AdjEidM);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int ele = 0; ele < Ele2d; ele++){
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
						StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
						ImposeDirichletBoundaryCondition(sr, StartPoint, Np, VertNfp, \
							jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
							J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, rx + (LocalEle - 1)*Np, \
							sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, , sy + (LocalEle - 1)*Np, \
							tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np);
						break;
					}
				}
			}
			else{//two or more layers included in vertical direction
				if (LocalEle == UpEle){//This is the top most cell, and L=0
					for (int i = 0; i < EleNumber; i++){
						if (LocalEle == (int)TempEToE[i]){
							StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
							ImposeDirichletBoundaryCondition(sr, StartPoint, Np, VertNfp, \
								jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
								J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, rx + (LocalEle - 1)*Np, \
								sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
								tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np);
							break;
						}
					}

					GetLocalDownFacialContribution(sr, jcs, TempEToE, EleNumber, LocalEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, BotEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np);

					GetLocalToDownFacialContribution(sr, jcs, TempEToE, EleNumber, LocalEle, DownEle, Np, VertNfp,\
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (DownEle - 1)*Np, J2d + ele*VertNfp, BotEidM, UpEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np,\
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (DownEle - 1)*Np, sx + (DownEle - 1)*Np, \
						ry + (DownEle - 1)*Np, sy + (DownEle - 1)*Np, tz + (DownEle - 1)*Np, Dr, Ds, Dt,\
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np,\
						K13 + (DownEle - 1)*Np, K23 + (DownEle - 1)*Np, K33 + (DownEle - 1)*Np);

				}
				else if (LocalEle == DownEle){// This is the bottom most cell.

					GetLocalUpFacialContribution(sr, jcs, TempEToE, EleNumber, LocalEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np);

					GetLocalToUpFacialContribution(sr, jcs, TempEToE, EleNumber, LocalEle, UpEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (DownEle - 1)*Np, J2d + ele*VertNfp, BotEidM, UpEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (UpEle - 1)*Np, sx + (UpEle - 1)*Np, \
						ry + (UpEle - 1)*Np, sy + (UpEle - 1)*Np, tz + (UpEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (UpEle - 1)*Np, K23 + (UpEle - 1)*Np, K33 + (UpEle - 1)*Np);
				}
				else{
					/*Upper surface*/
					GetLocalUpFacialContribution(sr, jcs, TempEToE, EleNumber, LocalEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, UpEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np);

					GetLocalToUpFacialContribution(sr, jcs, TempEToE, EleNumber, LocalEle, UpEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (DownEle - 1)*Np, J2d + ele*VertNfp, BotEidM, UpEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (UpEle - 1)*Np, sx + (UpEle - 1)*Np, \
						ry + (UpEle - 1)*Np, sy + (UpEle - 1)*Np, tz + (UpEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (UpEle - 1)*Np, K23 + (UpEle - 1)*Np, K33 + (UpEle - 1)*Np);

					/*Bottom surface*/
					GetLocalDownFacialContribution(sr, jcs, TempEToE, EleNumber, LocalEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (LocalEle - 1)*Np, J2d + ele*VertNfp, BotEidM, rx + (LocalEle - 1)*Np, \
						sx + (LocalEle - 1)*Np, ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, \
						tz + (LocalEle - 1)*Np, Dr, Ds, Dt, K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np);

					GetLocalToDownFacialContribution(sr, jcs, TempEToE, EleNumber, LocalEle, DownEle, Np, VertNfp, \
						jcs[(LocalEle - 1)*Np + 1] - jcs[(LocalEle - 1)*Np], M3d, M2d, \
						J + (DownEle - 1)*Np, J2d + ele*VertNfp, BotEidM, UpEidM, rx + (LocalEle - 1)*Np, sx + (LocalEle - 1)*Np, \
						ry + (LocalEle - 1)*Np, sy + (LocalEle - 1)*Np, tz + (LocalEle - 1)*Np, \
						rx + (DownEle - 1)*Np, sx + (DownEle - 1)*Np, \
						ry + (DownEle - 1)*Np, sy + (DownEle - 1)*Np, tz + (DownEle - 1)*Np, Dr, Ds, Dt, \
						K13 + (LocalEle - 1)*Np, K23 + (LocalEle - 1)*Np, K33 + (LocalEle - 1)*Np, \
						K13 + (DownEle - 1)*Np, K23 + (DownEle - 1)*Np, K33 + (DownEle - 1)*Np);

				}
			}
		}
	}
	free(UpEidM);
	free(BotEidM);
	free(InvSquaHeight);
	free(K33);
}

void GetLocalToDownFacialContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int DownEle, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, double *J, double *J2d, double *LocalEid, double *DownEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvDownMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvDownMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvDownMass3d, (ptrdiff_t)Np);
	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{31}\frac{\partial p}{\partial x}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial p}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial p}{\partial y}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial p}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Adjtz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, DownEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial p}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Localtz, Np);
	/*For term $k_{33}\frac{\partial p}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, DownEid, Np, Nfp);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvDownMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	int StartPoint;

	for (int i = 0; i < EleNumber; i++){
		if ((int)TempEToE[i] == DownEle){
			StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}

	free(Mass3d);
	free(InvDownMass3d);
	free(FacialMass2d);
	free(TempContribution);
	free(Contribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}

void GetLocalToUpFacialContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int UpEle, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, double *J, double *J2d, double *LocalEid, double *UpEid, double *Localrx, \
	double *Localsx, double *Localry, double *Localsy, double *Localtz, double *Adjrx, double *Adjsx, double *Adjry, double *Adjsy, \
	double *Adjtz, double *Dr, double *Ds, double *Dt, double *LocalK13, double *LocalK23, double *LocalK33, double *AdjK13, \
	double *AdjK23, double *AdjK33){

	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, mass3d, J, Np);
	double *InvUpMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvUpMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvUpMass3d, (ptrdiff_t)Np);
	double *FacialMass2d = malloc(Np2d*Np2d*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Np2d);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{31}\frac{\partial p}{\partial x}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localrx, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial p}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Adjry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Adjsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial p}{\partial y}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, Localry, Np);
	DiagMultiply(TempDiffMatrix, Ds, Localsy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial p}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1 * 0.5, Np*Nfp);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Adjtz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, AdjK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, UpEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial p}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, Localtz, Np);
	/*For term $k_{33}\frac{\partial p}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, DiffMatrix, LocalK33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1 * 0.5, Np*Nfp);
	AssembleContributionIntoRow(TempContribution, EdgeContribution, UpEid, Np, Nfp);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvUpMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	int StartPoint;

	for (int i = 0; i < EleNumber; i++){
		if ((int)TempEToE[i] == UpEle){
			StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}

	free(Mass3d);
	free(InvUpMass3d);
	free(FacialMass2d);
	free(TempContribution);
	free(Contribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
}



void GetLocalDownFacialContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, double *J, double *J2d, double *Eid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33){

	double *LocalMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalMass3d, M3d, J, Np);
	double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, -1*0.5, Np*Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvAdjMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	int StartPoint;

	for (int i = 0; i < EleNumber; i++){
		if ((int)TempEToE[i] == LocalEle){
			StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}

	free(LocalMass3d);
	free(InvLocalMass3d);
	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);

}


void GetLocalUpFacialContribution(double *dest, mwIndex *jcs, double *TempEToE, int EleNumber, int LocalEle, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, double *J, double *J2d, double *Eid, double *rx, double *sx, double *ry, \
	double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33){

	double *LocalMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalMass3d, M3d, J, Np);
	double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *EdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoColumn(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, EdgeContribution, (ptrdiff_t)Np);

	MultiplyByConstant(EdgeContribution, EdgeContribution, 0.5, Np*Nfp);

	AssembleContributionIntoRow(TempContribution, EdgeContribution, LocalEid, Np, Nfp);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvAdjMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	int StartPoint;

	for (int i = 0; i < EleNumber; i++){
		if ((int)TempEToE[i] == LocalEle){
			StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}

	free(LocalMass3d);
	free(InvLocalMass3d);
	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(EdgeContribution);
	free(Contribution);

}

void ImposeDirichletBoundaryCondition(double *dest, int StartPoint, int Np, int Nfp, \
	int NonzeroPerColumn, double *M3d, double *M2d, \
	double *J, double *J2d, double *LocalEid, double *rx, double *sx,\
	double *ry, double *sy, double *tz, double *Dr, double *Ds, double *Dt, double *K13, double *K23, double *K33){

	double *LocalMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalMass3d, M3d, J, Np);
	double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	double *DiffMatrix = malloc(Np*Np*sizeof(double));

	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *BoundaryEdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{31}\frac{\partial v}{\partial x}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{31}\frac{\partial v}{\partial x}$. Here, $k_{31} = k_{13}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{31}\frac{\partial p_h}{\partial x}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial v}{\partial y}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	/*For term $k_{32}\frac{\partial v}{\partial y}$. Here, $k_{32} = k_{23}$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{32}\frac{\partial p_h}{\partial y}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial v}{\partial \sigma}\right\}[p_h]_{\sigma} $$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	/*For term $k_{33}\frac{\partial v}{\partial \sigma}$.$*/
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, DiffMatrix, LocalEid, Nfp, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoColumn(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Nfp);

	/*For term $$ \left \{k_{33}\frac{\partial p_h}{\partial \sigma}\right\}[v]_{\sigma} $$*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialMass2d,
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, BoundaryEdgeContribution, (ptrdiff_t)Np);
	AssembleContributionIntoRow(TempContribution, BoundaryEdgeContribution, LocalEid, Np, Nfp);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvAdjMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, VolumnContribution, (ptrdiff_t)Np);

	int StartPoint;

	for (int i = 0; i < EleNumber; i++){
		if ((int)TempEToE[i] == LocalEle){
			StartPoint = jcs[(LocalEle - 1)*Np] + i*Np;
			AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
			break;
		}
	}

	free(LocalMass3d);
	free(InvLocalMass3d);
	free(FacialMass2d);
	free(TempContribution);
	free(DiffMatrix);
	free(TempDiffMatrix);
	free(FacialDiffMatrix);
	free(BoundaryEdgeContribution);
	free(Contribution);
}



void GetLocalToAdjacentFacialContributionInHorizontalDirection(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *M3d, double *M2d, double *AdjacentJ, double *J2d, double *LocalEid, double *AdjEid, double *Dt, double *LocalTz, \
	double *AdjTz, double *nx, double *ny, double *LocalK13, double *LocalK23, \
	double *AdjK13, double *AdjK23){
	// 这里的nx 和ny 由本地单元指出
	double *AdjMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(AdjMass3d, M3d, AdjacentJ, Np);
	double *InvAdjMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvAdjMass3d, AdjMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvAdjMass3d, (ptrdiff_t)Np);
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

	double *LocalDiff = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalDiff, Dt, LocalTz, Np);
	/*For term $k_{13}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, LocalK13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5*nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEid, Np, Nfp);

	DiagMultiply(LocalDiff, Dt, LocalTz, Np);
	/*For term $k_{23}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, LocalK23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, LocalEid, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, -0.5*ny[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, AdjEid, Np, Nfp);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvAdjMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, VolumnContribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, VolumnContribution, NonzeroPerColumn, Np);

	free(AdjMass3d);
	free(InvAdjMass3d);
	free(FacialMass2d);
	free(TempContribution);
	free(AdjDiff);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(LocalDiff);
}

void GetLocalFacialContributionInHorizontalDirection(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *M3d, double *M2d, double *J, double *J2d, double *EidM, double *Dt, double *tz, double *nx, \
	double *ny, double *K13, double *K23){

	double *LocalMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(LocalMass3d, M3d, J, Np);
	double *InvLocalMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvLocalMass3d, LocalMass3d, Np*Np*sizeof(double));
	MatrixInverse(InvLocalMass3d, (ptrdiff_t)Np);
	double *FacialMass2d = malloc(Nfp*Nfp*sizeof(double));
	DiagMultiply(FacialMass2d, M2d, J2d, Nfp);
	double *TempContribution = malloc(Np*Np*sizeof(double));
	memset(TempContribution, 0, Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));

	double *LocalDiff = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *InnerEdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{13}\frac{\partial v}{\partial \sigma}\right\}[p_h]_x $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{13}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, EidM, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*nx[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial v}{\partial \sigma}\right\}[p_h]_y $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{23}\frac{\partial v}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, EidM, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*ny[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*For term $$ \left \{k_{13}\frac{\partial p_h}{\partial \sigma}\right\}[v]_x $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{13}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, EidM, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial p_h}{\partial \sigma}\right\}[v]_y $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	/*For term $k_{23}\frac{\partial p_h}{\partial \sigma}$*/
	DiagMultiply(LocalDiff, LocalDiff, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, EidM, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*ny[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);
	free(LocalMass3d);
	free(InvLocalMass3d);
	free(FacialMass2d);
	free(TempContribution);
	free(LocalDiff);
	free(FacialDiffMatrix);
	free(InnerEdgeContribution);
	free(Contribution);
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