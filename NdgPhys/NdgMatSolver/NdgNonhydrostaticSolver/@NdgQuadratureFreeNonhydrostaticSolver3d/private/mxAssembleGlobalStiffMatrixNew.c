#include "SWENonhydrostatic3d.h"

void ImposeNonhydroBoundaryCondition(double *, mwIndex *, mwIndex *, double *, double *, \
	const mxArray *, const mxArray *, const mxArray *, const mxArray *, const mxArray *);


void ImposeNewmannBoundaryCondition(double *, mwIndex *, mwIndex *, int , \
	double *, double *, int , int , double *, double *, double *, double *, \
	double *, int , double *, double *, double *, const mxArray *, const mxArray *);

void SumInColumn(double *, double *, int);

int GetGlobalFace(int , int , double *, double *, int);

void GetFaceTypeAndFaceOrder(int *, int *, int *, double *, double *, signed char *, int);

void GetInverseSquareHeight(double *, double *, double, int);

int bin(int, int, int, mwIndex *);
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
					Dr, Ds, Dt, rx + ele*Np, sx + ele*Np, ry + ele*Np,\
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
				Mass3d, IEMb, J + ele*Np, Js + GlobalFace[i] * IENfp, LocalEidM, Dt, tz + ele*Np, Facialnx + i * IENfp,\
				Facialny + i*IENfp, K13 + ele*Np, K23 + ele*Np);

			for (int j = 0; j < EleNumber; j++){
				if ((int)TempEToE[j] == AdjEle[i]){
					StartPoint = jcs[ele*Np] + j*Np;
					GetLocalToAdjacentFacialContributionInHorizontalDirection(sr, StartPoint, Np, IENfp, jcs[ele*Np + 1] - jcs[ele*Np], \
						Mass3d, IEMb, J + (int)(TempEToE[j] - 1)*Np, Js + GlobalFace[i] * IENfp, LocalEidM, AdjEidM, Dt, tz+ele*Np, \
						tz + (int)(TempEToE[j] - 1)*Np, Facialnx + i*IENfp, Facialny + i*IENfp, K13 + ele*Np, K23 + ele*Np,\
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
































	double *SPNPS = mxGetPr(prhs[3]);
	mwIndex *jcSPNPS = mxGetJc(prhs[3]);
	mwIndex *irSPNPS = mxGetIr(prhs[3]);

	double *MSPNPX = mxGetPr(prhs[4]);
	mwIndex *jcMSPNPX = mxGetJc(prhs[4]);
	mwIndex *irMSPNPX = mxGetIr(prhs[4]);

	double *MSPNPY = mxGetPr(prhs[5]);
	mwIndex *jcMSPNPY = mxGetJc(prhs[5]);
	mwIndex *irMSPNPY = mxGetIr(prhs[5]);

	int row, col;

	row = (int)mxGetM(prhs[5]);
	col = (int)mxGetN(prhs[5]);

	double *sr;
	mwIndex *irs, *jcs;

	double *SPSPX = mxGetPr(prhs[8]);
	double *SPSPY = mxGetPr(prhs[9]);


	double *Height = mxGetPr(prhs[12]);
	int Np = (int)mxGetM(prhs[12]);
	int  K = (int)mxGetN(prhs[12]);

	double *InvSquaHeight = malloc(Np*K*sizeof(double));

	double Hcrit = mxGetScalar(prhs[13]);

	plhs[0] = mxCreateSparse(row, col, jcMSPNPY[col], mxREAL);
	sr = mxGetPr(plhs[0]);
	irs = mxGetIr(plhs[0]);
	jcs = mxGetJc(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		GetInverseSquareHeight(InvSquaHeight + Np*k, Height + Np*k, Hcrit, Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < col; i++)
	{
		int Num = (int)jcMSPNPY[i + 1] - (int)jcMSPNPY[i];
		double *temprhsu = sr + (int)jcMSPNPY[i];
		memset(temprhsu, 0, Num*sizeof(double));

		/*How many non-zero Np included in column i*/
		int NumSPNPX = (int)jcSPNPX[i + 1] - (int)jcSPNPX[i];
		for (int j = 0; j < NumSPNPX; j++){
			int rowIndex = (int)irSPNPX[j + (int)jcSPNPX[i]];
			int Order = bin(rowIndex, 0, Num - 1, irMSPNPY + (int)jcMSPNPY[i]);
			/*The position of the studied NP element in the input sparse matrix corresponding to non-hydrostatic pressure*/
			int Index = (int)jcSPNPX[i] + j;
			temprhsu[Order] = SPNPX[Index];
		}

		int NumSPNPY = (int)jcSPNPY[i + 1] - (int)jcSPNPY[i];
		for (int j = 0; j < NumSPNPY; j++){
			/*The position of the studied NP element in the input sparse matrix corresponding to non-hydrostatic pressure*/
			/*Index of the row that contains the studied non-zero Np in the studied column*/
			int rowIndex = (int)irSPNPY[j + (int)jcSPNPY[i]];
			int Order = bin(rowIndex, 0, Num - 1, irMSPNPY + (int)jcMSPNPY[i]);
			int Index = (int)jcSPNPY[i] + j;
			temprhsu[Order] = temprhsu[Order] + SPNPY[Index];
		}

		int NumSPNPS = (int)jcSPNPS[i + 1] - (int)jcSPNPS[i];
		for (int j = 0; j < NumSPNPS; j++){
			int rowIndex = (int)irSPNPS[j + (int)jcSPNPS[i]];
			int Order = bin(rowIndex, 0, Num - 1, irMSPNPY + (int)jcMSPNPY[i]);
			int Index = (int)jcSPNPS[i] + j;
//			temprhsu[Order] = temprhsu[Order] + InvSquaHeight[rowIndex]*SPNPS[Index];
			temprhsu[Order] = temprhsu[Order] + ( SQPSPX[rowIndex] + SQPSPY[rowIndex] + \
				InvSquaHeight[rowIndex] )*SPNPS[Index]; 
		}

		int NumMSPNPY = (int)jcMSPNPY[i + 1] - (int)jcMSPNPY[i];

		for (int j = 0; j < NumMSPNPY; j++){
			int rowIndex = (int)irMSPNPY[j + (int)jcMSPNPY[i]];
			int Order = bin(rowIndex, 0, Num - 1, irMSPNPY + (int)jcMSPNPY[i]);
			int Index = (int)jcMSPNPY[i] + j;
//			temprhsu[Order] = temprhsu[Order] + 2 * (PSPX[rowIndex] * MSPNPX[Index] + \
				PSPY[rowIndex] * MSPNPY[Index]);
			temprhsu[Order] = temprhsu[Order] + (PSPX[rowIndex] * MSPNPX[Index] + \
				PSPY[rowIndex] * MSPNPY[Index]);
		}

		int NumPNPS = (int)jcPNPS[i + 1] - (int)jcPNPS[i];

		for (int j = 0; j < NumPNPS; j++){
			int rowIndex = (int)irPNPS[j + (int)jcPNPS[i]];
			int Order = bin(rowIndex, 0, Num - 1, irMSPNPY + (int)jcMSPNPY[i]);
			int Index = (int)jcPNPS[i] + j;
			temprhsu[Order] = temprhsu[Order] + (SPSPX[rowIndex] + \
				SPSPY[rowIndex]) * PNPS[Index];
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int i = 0; i < col + 1; i++){
		jcs[i] = jcMSPNPY[i];
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int i = 0; i < (int)jcMSPNPY[col]; i++){
		irs[i] = irMSPNPY[i];
	}


	ImposeNonhydroBoundaryCondition(sr, irs, jcs, PSPX, PSPY, prhs[14], prhs[15], prhs[16], prhs[17], prhs[18]);

	free(InvSquaHeight);
	free(K33);
}

void GetLocalToAdjacentFacialContributionInHorizontalDirection(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *M3d, double *M2d, double *AdjacentJ, double *Js, double *LocalEid, double *AdjEid, double *Dt, double *LocalTz, \
	double *AdjTz, double *nx, double *ny, double *LocalK13, double *LocalK23, \
	double *AdjK13, double *AdjK23){\

	// 这里的nx 和ny 由本地单元指出

}

void FetchFacialValue(double *dest, double *Eid, double *Source, int Nfp){
	for (int i = 0; i < Nfp; i++){
		dest[i] = Source[(int)Eid[i] - 1];
	}
}

void GetLocalFacialContributionInHorizontalDirection(double *dest, int StartPoint, int Np, int Nfp, int NonzeroPerColumn, \
	double *M3d, double *M2d, double *J, double *Js, double *EidM, double *Dt, double *tz, double *nx, \
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

	double *LocalDiff = malloc(Np*Np*sizeof(double));

	double *FacialDiffMatrix = malloc(Np*Nfp*sizeof(double));

	double *InnerEdgeContribution = malloc(Np*Nfp*sizeof(double));
	/*For term $$ \left \{k_{13}\frac{\partial v}{\partial \sigma}\right\}[p_h]_x $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	DiagMultiply(LocalDiff, LocalDiff, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, EidM, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*nx[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial v}{\partial \sigma}\right\}[p_h]_y $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	DiagMultiply(LocalDiff, LocalDiff, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, EidM, Nfp, Np);
	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Nfp, (ptrdiff_t)Nfp, 1.0, FacialDiffMatrix,
		(ptrdiff_t)Nfp, FacialMass2d, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Np);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*ny[0], Np*Nfp);
	AssembleContributionIntoColumn(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*For term $$ \left \{k_{13}\frac{\partial p_h}{\partial \sigma}\right\}[v]_x $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	DiagMultiply(LocalDiff, LocalDiff, K13, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, EidM, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*nx[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*For term $$ \left \{k_{23}\frac{\partial p_h}{\partial \sigma}\right\}[v]_y $$*/
	DiagMultiply(LocalDiff, Dt, tz, Np);
	DiagMultiply(LocalDiff, LocalDiff, K23, Np);
	AssembleFacialDiffMatrix(FacialDiffMatrix, LocalDiff, EidM, Nfp, Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Nfp, (ptrdiff_t)Np, (ptrdiff_t)Nfp, 1.0, FacialMass2d, \
		(ptrdiff_t)Nfp, FacialDiffMatrix, (ptrdiff_t)Nfp, 0.0, InnerEdgeContribution, (ptrdiff_t)Nfp);
	MultiplyByConstant(InnerEdgeContribution, InnerEdgeContribution, 0.5*ny[0], Np*Nfp);
	AssembleContributionIntoRow(TempContribution, InnerEdgeContribution, LocalEidM, Np, Nfp);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvLocalMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, VolumnContribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, VolumnContribution, NonzeroPerColumn, Np);

}


void GetLocalVolumnIntegralTerm(double *dest, int StartPoint, int Np, int NonzeroPerColumn, double *M3d, \
	double *Dr, double *Ds, double *Dt, double *rx, double *sx, double *ry, double *sy, double *tz, double *J, double *K13, \
	double *K23, double *K33){
	double *DiffMatrix = malloc(Np*Np*sizeof(double));
	double *TempDiffMatrix = malloc(Np*Np*sizeof(double));
	double *VertDiffMatrix = malloc(Np*Np*sizeof(double));
	double *Mass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(Mass3d, M3d, J, Np);
	double *InvMass3d = malloc(Np*Np*sizeof(double));
	memcpy(InvMass3d, Mass3d, Np*Np*sizeof(double));
	MatrixInverse(InvMass3d, (ptrdiff_t)Np);
	double *TempContributionBuff = malloc(Np*Np*sizeof(double));
	double *ContributionBuff = malloc(Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	memset(Contribution, 0, Np*Np*sizeof(double));

	double *VolumnContribution = malloc(Np*Np*sizeof(double));

	/*For term $$k_{13}\frac{\partial v}{\partial \sigma}\frac{\partial p_h}{\partial x}$$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);
	DiagMultiply(VertDiffMatrix, VertDiffMatrix, K13, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, VertDiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
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
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, DiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{31}\frac{\partial v}{\partial x}\frac{\partial p_h}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, Dr, rx, Np);
	DiagMultiply(TempDiffMatrix, Ds, sx, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K31, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VertDiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{32}\frac{\partial v}{\partial y}\frac{\partial p_h}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, Dr, ry, Np);
	DiagMultiply(TempDiffMatrix, Ds, sy, Np);
	Add(DiffMatrix, DiffMatrix, TempDiffMatrix, Np*Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K32, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VertDiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*For term $$k_{33}\frac{\partial v}{\partial \sigma}\frac{\partial p_h}{\partial \sigma}$$*/
	DiagMultiply(DiffMatrix, Dt, tz, Np);
	DiagMultiply(DiffMatrix, DiffMatrix, K33, Np);

	DiagMultiply(VertDiffMatrix, Dt, tz, Np);

	MatrixMultiply("T", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, DiffMatrix,
		(ptrdiff_t)Np, Mass3d, (ptrdiff_t)Np, 0.0, TempContributionBuff, (ptrdiff_t)Np);
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, TempContributionBuff,
		(ptrdiff_t)Np, VertDiffMatrix, (ptrdiff_t)Np, 0.0, ContributionBuff, (ptrdiff_t)Np);
	MultiplyByConstant(ContributionBuff, ContributionBuff, -1.0, Np*Np);

	Add(Contribution, Contribution, ContributionBuff, Np*Np);

	/*Multiply by inverse mass matrix*/
	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvMass3d,
		(ptrdiff_t)Np, Contribution, (ptrdiff_t)Np, 0.0, VolumnContribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, VolumnContribution, NonzeroPerColumn, Np);

	free(DiffMatrix);
	free(TempDiffMatrix);
	free(VertDiffMatrix);
	free(Mass3d);
	free(InvMass3d);
	free(TempContributionBuff);
	free(ContributionBuff);
	free(Contribution);
	free(VolumnContribution);
}


void ImposeNonhydroBoundaryCondition(double *dest, mwIndex *irs, mwIndex *jcs, double *PSPX, double *PSPY, \
	const mxArray *BoundaryEdge2d, const mxArray *mesh, const mxArray *cell, const mxArray *BoundaryEdge, const mxArray *Arrayftype2d){


	mxArray *TempNlayer = mxGetField(mesh, 0, "Nz");
	int Nlayer = (int)mxGetScalar(TempNlayer);
	mxArray *TempEToE = mxGetField(mesh, 0, "EToE");
	double *EToE = mxGetPr(TempEToE);
	mxArray *TempJ = mxGetField(mesh, 0, "J");
	double *J = mxGetPr(TempJ);

	mxArray *TempJs = mxGetField(BoundaryEdge, 0, "Js");
	double *Js = mxGetPr(TempJs);
	int Nfp = (int)mxGetM(TempJs);
	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempFToF = mxGetField(BoundaryEdge, 0, "FToF");
	double *FToF = mxGetPr(TempFToF);
	mxArray *TempFToE = mxGetField(BoundaryEdge, 0, "FToE");
	double *FToE = mxGetPr(TempFToE);
	mxArray *Tempnx = mxGetField(BoundaryEdge, 0, "nx");
	double *nx = mxGetPr(Tempnx);
	mxArray *Tempny = mxGetField(BoundaryEdge, 0, "ny");
	double *ny = mxGetPr(Tempny);
	mxArray *TempLMass2d = mxGetField(BoundaryEdge, 0, "M");
	double  *LMass2d = mxGetPr(TempLMass2d);

	mxArray *TempFToF2d = mxGetField(BoundaryEdge2d, 0, "FToF");
	double *FToF2d = mxGetPr(TempFToF2d);
	mxArray *TempFToE2d = mxGetField(BoundaryEdge2d, 0, "FToE");
	double *FToE2d = mxGetPr(TempFToE2d);
	/*
	mxArray *Tempftype2d = mxGetField(BoundaryEdge2d, 0, "ftype");	
    signed char *ftype2d = (signed char *)mxGetData(Tempftype2d);
	*/
	signed char *ftype2d = (signed char *)mxGetData(Arrayftype2d);
	mxArray *TempBENe2d = mxGetField(BoundaryEdge2d, 0, "Ne");
	int BENe2d = (int)mxGetScalar(TempBENe2d);

	mxArray *TempFmask = mxGetField(cell, 0, "Fmask");
	double *Fmask = mxGetPr(TempFmask);
	int maxNfp = (int)mxGetM(TempFmask);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempMass3d = mxGetField(cell, 0, "M");
	double  *Mass3d = mxGetPr(TempMass3d);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int edge = 0; edge < BENe2d; edge++){
		int Flag = 0;
		int face2d = 0;
		int ele2d = 0;
		int Nface2d = Nface - 2;

		double *FpIndex = malloc(Nfp*sizeof(double));

		GetFaceTypeAndFaceOrder(&Flag, &face2d, &ele2d, FToF2d, FToE2d, ftype2d, edge);

		for (int p = 0; p < Nfp; p++){
			FpIndex[p] = Fmask[maxNfp*(face2d - 1) + p];
		}

		if(Flag ==1 ){ //Impose Newmann boundary condition here.

			double *TempEToE = NULL, *TempJ = NULL, *TempJs = NULL;

			int GlobalFace, LocalEle;

			for (int L = 0; L < Nlayer; L++){

				LocalEle = (ele2d - 1)*Nlayer + L + 1;

				TempEToE = EToE + (ele2d - 1)*Nlayer*Nface + L*Nface;

				TempJ = J + (ele2d - 1)*Nlayer*Np + L*Np;

				GlobalFace = GetGlobalFace(face2d, BENe, FToE, FToF, LocalEle);

				TempJs = Js + GlobalFace * Nfp;

				ImposeNewmannBoundaryCondition(dest, irs, jcs, LocalEle, \
					PSPX, PSPY, Np, Nfp, Mass3d, TempJ, TempJs, LMass2d, \
					TempEToE, Nface, FpIndex, nx + GlobalFace * Nfp, ny + GlobalFace * Nfp, \
					mesh, cell);
			}		
		}
		free(FpIndex);
	}
}

void ImposeNewmannBoundaryCondition(double *dest, mwIndex *Irs, mwIndex *Jcs, int LocalEle,\
	double *PsPx, double *PsPy, int Np, int Nfp, double *M3d, double *J, double *Js, double *M2d, \
	double *EToE, int Nface, double *FpIndex, double *nx, double *ny, const mxArray *mesh, const mxArray *cell){

	double *TempPsPx = malloc(Nfp*sizeof(double));
	double *TempPsPy = malloc(Nfp*sizeof(double));
	double *TempPsP  = malloc(Nfp*sizeof(double));
	/*$\frac{\partial \sigma}{\partial x}$*/
	FetchFacialData(TempPsPx, PsPx + (LocalEle - 1)*Np, FpIndex, Nfp);
	/*$\frac{\partial \sigma}{\partial y}$*/
	FetchFacialData(TempPsPy, PsPy + (LocalEle - 1)*Np, FpIndex, Nfp);
	/*$\frac{\partial \sigma}{\partial x}n_x$*/
	DotProduct(TempPsPx, TempPsPx, nx, Nfp);
	/*$\frac{\partial \sigma}{\partial y}n_y$*/
	DotProduct(TempPsPy, TempPsPy, ny, Nfp);
	/*$-\frac{\partial \sigma}{\partial x}n_x$*/
	DotDivideByConstant(TempPsPx, TempPsPx, -1.0, Nfp);
	/*$-\frac{\partial \sigma}{\partial y}n_y$*/
	DotDivideByConstant(TempPsPy, TempPsPy, -1.0, Nfp);

	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempTz = mxGetField(mesh, 0, "tz");
	double  *tz = mxGetPr(TempTz);

	/*Withdraw the data in $\frac{\partial p}{\partial \sigma}$, and store them in TempPNPS*/
	double *Dz = malloc(Np*Np*sizeof(double));
	DiagMultiply(Dz, Dt, tz + (LocalEle - 1)*Np, Np);
	double *TempFacialData = malloc(Nfp*sizeof(double));
	double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));
	double *TempContribution = malloc(Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	double *InvEleMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(InvEleMass3d, M3d, J, Np);
	MatrixInverse(InvEleMass3d, (ptrdiff_t)Np);
	DiagMultiply(EleMass2d, M2d, Js, Nfp);

	int UniNum = 0, StartPoint;
	
	/*Find the exact place where to fill in the data, and the place is stored in StartPoint*/
	double *TempEToE = malloc((Nface + 1)*sizeof(double));

	FindUniqueElementAndSortOrder(TempEToE, EToE, &UniNum, Nface, LocalEle);

	int NonzeroPerColumn = Jcs[(LocalEle - 1)*Np + 1] - Jcs[(LocalEle - 1)*Np];

	for (int j = 0; j < UniNum; j++){
		if ((int)TempEToE[j] == LocalEle)
			StartPoint = (int)Jcs[(LocalEle - 1)*Np] + j*Np;
	}
	double *ContributionPerPoint = malloc(Np*sizeof(double));
	double *TempMass2d = malloc(Nfp*Nfp*sizeof(double));

	double column = 1.0;

	memset(TempContribution, 0, Np*Np*sizeof(double));
	/*Begin to assemble the local contribution in to the Newmann boundary condition*/
	for (int j = 0; j < Np; j++){
		memset(ContributionPerPoint, 0, Np*sizeof(double));
		
		FetchFacialData(TempFacialData, Dz + j*Np, FpIndex, Nfp);
		/*$-\frac{\partial \sigma}{\partial x}n_x\frac{\partial p}{\partial \sigma}$*/
		DotProduct(TempPsPx, TempPsPx, TempFacialData, Nfp);
		/*$-\frac{\partial \sigma}{\partial y}n_y\frac{\partial p}{\partial \sigma}$*/
		DotProduct(TempPsPy, TempPsPy, TempFacialData, Nfp);
		/*$-\frac{\partial \sigma}{\partial x}n_x\frac{\partial p}{\partial \sigma}-\frac{\partial \sigma}{\partial y}n_y\frac{\partial p}{\partial \sigma}$*/
		Add(TempPsP, TempPsPx, TempPsPy, Nfp);

		DiagMultiply(TempMass2d, EleMass2d, TempPsP, Nfp);

		SumInColumn(ContributionPerPoint, TempMass2d, Nfp);

		AssembleContributionIntoColumn(TempContribution + j*Np, ContributionPerPoint, &column, Np, 1);

	}

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvEleMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(TempPsPx);
	free(TempPsPy);
	free(TempPsP);
	free(Dz);
	free(TempFacialData);
	free(EleMass2d);
	free(TempMass2d);
	free(TempContribution);
	free(Contribution);
	free(InvEleMass3d);
	free(TempEToE);
	free(ContributionPerPoint);
}

void SumInColumn(double *dest, double *Source, int Np){
	for (int col = 0; col < Np; col++){
		for (int Row = 0; Row < Np; Row++){
			dest[col] += Source[col*Np + Row];
		}
	}
}

int GetGlobalFace(int face, int Ne, double *FToE, double *FToF, int LocalElement){
	for (int i = 0; i < Ne; i++){
		if ((int)FToF[2 * i] == face && (int)FToE[2 * i] == LocalElement){
			return i;
			break;
		}
	}
	return -1; //failed
}

void GetFaceTypeAndFaceOrder(int *flag, int *face, int *ele, double *FToF, double *FToE, signed char *ftype, int edge){
	(*ele) = (int)FToE[2 * edge];
	(*face) = (int)FToF[2 * edge];
	/*At present, only Newmann bounary is considered*/
	if ((NdgEdgeType)ftype[edge] == NdgEdgeSlipWall){
		(*flag) = 1;
	}
	else{
		(*flag) = 0;
	}
}

/*The following function is used to find the order of the input parameter aim in matrix a*/
int bin(int aim, int low, int high, mwIndex *a)
{
	int mid;
	while (low <= high)
	{
		mid = (low + high) / 2;
		if ((int)a[mid] == aim)return mid;
		else if ((int)a[mid]<aim)low = mid + 1;
		else high = mid - 1;
	}
	return -1;//Not Find
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