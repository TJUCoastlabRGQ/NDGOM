#include "SWENonhydrostatic3d.h"

void ImposeNonhydroBoundaryCondition(double *, mwIndex *, mwIndex *, double *, double *, double *, \
	mwIndex *, mwIndex *, const mxArray *, const mxArray *, \
	const mxArray *, const mxArray *);


void ImposeNewmannBoundaryCondition(double *, mwIndex *, mwIndex *, int, int, double *, mwIndex *, \
	mwIndex *, double *, double *, int, int, double *, double *, double *, double *, \
	double *, int, double *, double *, double *);

void SumInColumn(double *, double *, int);

int GetGlobalFace(int , int , double *, double *, int);

void GetFaceTypeAndFaceOrder(int *, int *, int *, double *, double *, signed char *, int);

void GetInverseSquareHeight(double *, double *, double, int);

int bin(int, int, int, mwIndex *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *PNPS = mxGetPr(prhs[0]);
	mwIndex *jcPNPS = mxGetJc(prhs[0]);
	mwIndex *irPNPS = mxGetIr(prhs[0]);

	double *SPNPX = mxGetPr(prhs[1]);
	mwIndex *jcSPNPX = mxGetJc(prhs[1]);
	mwIndex *irSPNPX = mxGetIr(prhs[1]);

	double *SPNPY = mxGetPr(prhs[2]);
	mwIndex *jcSPNPY = mxGetJc(prhs[2]);
	mwIndex *irSPNPY = mxGetIr(prhs[2]);

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

	double *PSPX = mxGetPr(prhs[6]);
	double *PSPY = mxGetPr(prhs[7]);

	double *SPSPX = mxGetPr(prhs[8]);
	double *SPSPY = mxGetPr(prhs[9]);
	double *SQPSPX = mxGetPr(prhs[10]);
	double *SQPSPY = mxGetPr(prhs[11]);

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
			temprhsu[Order] = temprhsu[Order] + (pow(SQPSPX[rowIndex], 2) + pow(SQPSPY[rowIndex], 2) + \
				InvSquaHeight[rowIndex])*SPNPS[Index]; 
		}

		int NumMSPNPY = (int)jcMSPNPY[i + 1] - (int)jcMSPNPY[i];

		for (int j = 0; j < NumMSPNPY; j++){
			int rowIndex = (int)irMSPNPY[j + (int)jcMSPNPY[i]];
			int Order = bin(rowIndex, 0, Num - 1, irMSPNPY + (int)jcMSPNPY[i]);
			int Index = (int)jcMSPNPY[i] + j;
			temprhsu[Order] = temprhsu[Order] + 2 * (PSPX[rowIndex] * MSPNPX[Index] + \
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


	ImposeNonhydroBoundaryCondition(sr, irs, jcs, PSPX, PSPY, PNPS, irPNPS, jcPNPS, \
	     prhs[16], prhs[17], prhs[18], prhs[19]);

	free(InvSquaHeight);
}


void ImposeNonhydroBoundaryCondition(double *dest, mwIndex *irs, mwIndex *jcs, double *PSPX, double *PSPY, double *PNPS, \
	mwIndex *irPNS, mwIndex *jcPNS,  const mxArray *BoundaryEdge2d, const mxArray *mesh, \
	const mxArray *cell, const mxArray *BoundaryEdge){


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
	mxArray *Tempftype2d = mxGetField(BoundaryEdge2d, 0, "ftype");
	signed char *ftype2d = mxGetPr(Tempftype2d);
	mxArray *TempBENe2d = mxGetField(BoundaryEdge2d, 0, "Ne");
	double *BENe2d = mxGetPr(TempBENe2d);

	mxArray *TempFmask = mxGetField(cell, 0, "Fmask");
	double *Fmask = mxGetPr(TempFmask);
	int maxNfp = (int)mxGetM(TempFmask);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetM(TempNface);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempMass3d = mxGetField(cell, 0, "M");
	double  *Mass3d = mxGetPr(TempMass3d);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int edge = 0; edge < (int)BENe2d; edge++){
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


			int GlobalFace;

			for (int L = 0; L < Nlayer; L++){
				//Index of the studied element
				int LocalEle = (ele2d - 1)*Nlayer + L + 1;
				//Index of the element that located upside of the studied element 
				int UpEle = (int)EToE[(ele2d - 1)*Nlayer*Nface + L*Nface + Nface - 1];
				//Index of the element that located downside of the studied element
				int DownEle = (int)EToE[(ele2d - 1)*Nlayer*Nface + L*Nface + Nface - 2];

				TempEToE = EToE + (ele2d - 1)*Nlayer*Nface + L*Nface;

				TempJ = J + (ele2d - 1)*Nlayer*Np + L*Np;

				GlobalFace = GetGlobalFace(face2d, BENe, FToE, FToF, LocalEle);

				TempJs = Js + GlobalFace * Nfp;

				ImposeNewmannBoundaryCondition(dest, irs, jcs, LocalEle, LocalEle, PNPS, irPNS, jcPNS, \
					PSPX, PSPY, Np, Nfp, Mass3d, TempJ, TempJs, LMass2d, \
					TempEToE, Nface, FpIndex, nx + GlobalFace * Nfp, ny + GlobalFace * Nfp);

				if (UpEle == DownEle){ // we have only one layer in vertical direction, since L=0
					// Doing nothing
				}
				else{//two or more layers included in vertical direction
					if (LocalEle == UpEle){//This is the top most cell, and L=0

						TempJ = J + (ele2d - 1)*Nlayer*Np + (L+1)*Np;

						GlobalFace = GetGlobalFace(face2d, BENe, FToE, FToF, DownEle);

						TempJs = Js + GlobalFace * Nfp;
						/*Consider the impact of local element have on the element downside*/
						ImposeNewmannBoundaryCondition(dest, irs, jcs, LocalEle, DownEle, PNPS, irPNS, jcPNS, \
							PSPX, PSPY, Np, Nfp, Mass3d, TempJ, TempJs, LMass2d, \
							TempEToE, Nface, FpIndex, nx + GlobalFace * Nfp, ny + GlobalFace * Nfp);
					}
					else if (LocalEle == DownEle){// This is the bottom most cell.
						/*For this situation, this is the bottom most cell, since Newmann boundary
						condition is added, for the first order derivative about non-hydrostatic
						pressure, we do nothing here*/
						TempJ = J + (ele2d - 1)*Nlayer*Np + (L - 1)*Np;

						GlobalFace = GetGlobalFace(face2d, BENe, FToE, FToF, UpEle);

						TempJs = Js + GlobalFace * Nfp;
						/*Consider the impact of local element have on the element upside*/
						ImposeNewmannBoundaryCondition(dest, irs, jcs, LocalEle, UpEle, PNPS, irPNS, jcPNS, \
							PSPX, PSPY, Np, Nfp, Mass3d, TempJ, TempJs, LMass2d, \
							TempEToE, Nface, FpIndex, nx + GlobalFace * Nfp, ny + GlobalFace * Nfp);
					}
					else{
						/*Consider the influence of local element have on the adjacent cell upside*/
						TempJ = J + (ele2d - 1)*Nlayer*Np + (L - 1)*Np;

						GlobalFace = GetGlobalFace(face2d, BENe, FToE, FToF, UpEle);

						TempJs = Js + GlobalFace * Nfp;
						/*Consider the impact of local element have on the element upside*/
						ImposeNewmannBoundaryCondition(dest, irs, jcs, LocalEle, UpEle, PNPS, irPNS, jcPNS, \
							PSPX, PSPY, Np, Nfp, Mass3d, TempJ, TempJs, LMass2d, \
							TempEToE, Nface, FpIndex, nx + GlobalFace * Nfp, ny + GlobalFace * Nfp);

						/*Consider the influence of local element have on the adjacent cell downside*/
						TempJ = J + (ele2d - 1)*Nlayer*Np + (L + 1)*Np;

						GlobalFace = GetGlobalFace(face2d, BENe, FToE, FToF, DownEle);

						TempJs = Js + GlobalFace * Nfp;
						/*Consider the impact of local element have on the element downside*/
						ImposeNewmannBoundaryCondition(dest, irs, jcs, LocalEle, DownEle, PNPS, irPNS, jcPNS, \
							PSPX, PSPY, Np, Nfp, Mass3d, TempJ, TempJs, LMass2d, \
							TempEToE, Nface, FpIndex, nx + GlobalFace * Nfp, ny + GlobalFace * Nfp);
					}
				}
			}		
		}
		free(FpIndex);
	}
}

void ImposeNewmannBoundaryCondition(double *dest, mwIndex *Irs, mwIndex *Jcs, int LocalEle, int AdjOrLocalEle, double *PsSource, mwIndex *irPNS, \
	mwIndex *jcPNS, double *PsPx, double *PsPy, int Np, int Nfp, double *M3d, double *J, double *Js, double *M2d, \
	double *EToE, int Nface, double *FpIndex, double *nx, double *ny){

	double *TempPsPx = malloc(Nfp*sizeof(double));
	double *TempPsPy = malloc(Nfp*sizeof(double));
	double *TempPsP  = malloc(Nfp*sizeof(double));
	/*$\frac{\partial \sigma}{\partial x}$*/
	FetchFacialData(TempPsPx, PsPx + (AdjOrLocalEle - 1)*Np, FpIndex, Nfp);  
	/*$\frac{\partial \sigma}{\partial y}$*/
	FetchFacialData(TempPsPy, PsPy + (AdjOrLocalEle - 1)*Np, FpIndex, Nfp); 
	/*$\frac{\partial \sigma}{\partial x}n_x$*/
	DotProduct(TempPsPx, TempPsPx, nx, Nfp);
	/*$\frac{\partial \sigma}{\partial y}n_y$*/
	DotProduct(TempPsPy, TempPsPy, ny, Nfp);
	/*$-\frac{\partial \sigma}{\partial x}n_x$*/
	DotDivideByConstant(TempPsPx, TempPsPx, -1.0, Nfp);
	/*$-\frac{\partial \sigma}{\partial y}n_y$*/
	DotDivideByConstant(TempPsPy, TempPsPy, -1.0, Nfp);

	/*Withdraw the data in $\frac{\partial p}{\partial \sigma}$, and store them in TempPNPS*/
	double *TempPNPS = malloc(Np*Np*sizeof(double));
	double *TempEToEInVert = malloc(3 * sizeof(double));
	int UniNum = 0;
	int StartPoint;
	double *EToEInVert = malloc(2 * sizeof(double));
	EToEInVert[0] = EToE[Nface - 2];
	EToEInVert[1] = EToE[Nface - 1];
	FindUniqueElementAndSortOrder(TempEToEInVert, EToEInVert, &UniNum, 2, LocalEle);
	for (int j = 0; j < UniNum; j++){
		if ((int)TempEToEInVert[j] == AdjOrLocalEle)//DownEle
			StartPoint = (int)jcPNS[(LocalEle - 1)*Np] + j*Np;
	}
	int NonzeroPerColumn = jcPNS[(LocalEle - 1)*Np + 1] - jcPNS[(LocalEle - 1)*Np];
	FetchDataInSparseMatrix(TempPNPS, PsSource + StartPoint, NonzeroPerColumn, Np);
	

	double *TempFacialData = malloc(Nfp*sizeof(double));
	double *EleMass2d = malloc(Nfp*Nfp*sizeof(double));
	double *TempContribution = malloc(Np*Np*sizeof(double));
	double *Contribution = malloc(Np*Np*sizeof(double));
	double *InvEleMass3d = malloc(Np*Np*sizeof(double));
	DiagMultiply(InvEleMass3d, M3d, J, Np);//J down
	MatrixInverse(InvEleMass3d, (ptrdiff_t)Np);
	DiagMultiply(EleMass2d, M2d, Js, Nfp); //Js down
	
	/*Find the exact place where to fill in the data, and the place is stored in StartPoint*/
	double *TempEToE = malloc((Nface + 1)*sizeof(double));
	FindUniqueElementAndSortOrder(TempEToE, EToE, &UniNum, Nface, LocalEle);
	NonzeroPerColumn = Jcs[(LocalEle - 1)*Np + 1] - Jcs[(LocalEle - 1)*Np];

	for (int j = 0; j < UniNum; j++){
		if ((int)TempEToE[j] == AdjOrLocalEle)
			StartPoint = (int)Jcs[(LocalEle - 1)*Np] + j*Np;
	}
	double *ContributionPerPoint = malloc(Np*sizeof(double));

	double column = 1.0;

	memset(TempContribution, 0, Np*Np*sizeof(double));
	/*Begin to assemble the local contribution in to the Newmann boundary condition*/
	for (int j = 0; j < Np; j++){
		memset(ContributionPerPoint, 0, Np*sizeof(double));
		
		FetchFacialData(TempFacialData, TempPNPS + j*Np, FpIndex, Nfp);
		/*$-\frac{\partial \sigma}{\partial x}n_x\frac{\partial p}{\partial \sigma}$*/
		DotProduct(TempPsPx, TempPsPx, TempFacialData, Nfp);
		/*$-\frac{\partial \sigma}{\partial y}n_y\frac{\partial p}{\partial \sigma}$*/
		DotProduct(TempPsPy, TempPsPy, TempFacialData, Nfp);
		/*$-\frac{\partial \sigma}{\partial x}n_x\frac{\partial p}{\partial \sigma}-\frac{\partial \sigma}{\partial y}n_y\frac{\partial p}{\partial \sigma}$*/
		Add(TempPsP, TempPsPx, TempPsPy, Nfp);

		DiagMultiply(EleMass2d, EleMass2d, TempPsP, Nfp);

		SumInColumn(ContributionPerPoint, EleMass2d, Nfp);

		AssembleContributionIntoColumn(TempContribution + j*Np, EleMass2d, &column, Np, 1);

	}

	MatrixMultiply("N", "N", (ptrdiff_t)Np, (ptrdiff_t)Np, (ptrdiff_t)Np, 1.0, InvEleMass3d,
		(ptrdiff_t)Np, TempContribution, (ptrdiff_t)Np, 0.0, Contribution, (ptrdiff_t)Np);

	AssembleContributionIntoSparseMatrix(dest + StartPoint, Contribution, NonzeroPerColumn, Np);

	free(TempPsPx);
	free(TempPsPy);
	free(TempPsP);
	free(TempPNPS);
	free(TempEToEInVert);
	free(EToEInVert);
	free(TempFacialData);
	free(EleMass2d);
	free(TempContribution);
	free(Contribution);
	free(InvEleMass3d);
	free(TempEToE);
	free(ContributionPerPoint);

}

void SumInColumn(double *dest, double *Source, int Np){
	for (int col = 0; col < Np; col++){
		for (int Row = 0; Row < Np; Row++){
			dest[col] = dest[col*Np + Row];
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
			dest[i] = 1 / source[i] / source[i];
		}
		else{
			dest[i] = 0;
		}
	}
}