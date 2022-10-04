#include "SWENonhydrostatic3d.h"

int bin(int , int , int , mwIndex *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int Col = (int)mxGetN(prhs[0]);
	int Row = (int)mxGetM(prhs[0]);
	mwIndex *Ir = mxGetIr(prhs[0]);
	mwIndex *Jc = mxGetJc(prhs[0]);
	double *InputData = mxGetPr(prhs[0]);
	
	plhs[0] = mxCreateSparse(Row, Col, Jc[Col], mxREAL);
	double *OutputData = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);
	mwIndex *jcs = mxGetJc(plhs[0]);

	memcpy(irs, Ir, Jc[Col] * sizeof(mwIndex));
	memcpy(jcs, Jc, (Col + 1)*sizeof(mwIndex));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif  
	for (int Lcol = 0; Lcol < Col; Lcol++){
		double LTempdate, UTempdate;
		int LStartPoint = (int)Jc[Lcol], UStartPoint;
		int LRow, URow;
		for (int row = 0; row < (int)(Jc[Lcol + 1]) - LStartPoint; row++){
			/*We only consider the lower triangle part*/
			if ((int)(Ir[LStartPoint + row]) > Lcol){
				LTempdate = InputData[LStartPoint + row];
				LRow = (int)(Ir[LStartPoint + row]);
				UStartPoint = (int)Jc[LRow];
				int P = bin(Lcol, 0, (int)Jc[LRow + 1] - UStartPoint, Ir + UStartPoint);
				UTempdate = InputData[UStartPoint + P];
				/*Multiply the data by -1 to form the positive definite matrix*/
				/*The lower triangle part*/
				OutputData[LStartPoint + row] = -1.0 * (LTempdate + UTempdate) / 2.0;
				/*The upper triangle part*/
				OutputData[UStartPoint + P] = -1.0 * (LTempdate + UTempdate) / 2.0;
			}
			else if ((int)(Ir[LStartPoint + row]) == Lcol){
				/*For date on the diagonal, multiply the date by -1.0*/
				OutputData[LStartPoint + row] = -1.0 * InputData[LStartPoint + row];
			}
				
		}
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