#include "SWENonhydrostatic3d.h"
/*
This function is used to assemble the following terms into the 
global stiff matrix to get the final global stiff matrix. These
terms are given as:
$$\frac{1}{D}\frac{\partial D}{\partial x^*}\frac{\partial p}{\partial x} + 
\frac{1}{D}\frac{\partial D}{\partial y^*}\frac{\partial p}{\partial y} -
\left(\frac{1}{D}\frac{\partial z_b}{\partial x}+\frac{1+\sigma}{D}\frac{\partial D}{\partial x^*}\right)\frac{1}{D}\frac{\partial D}{\partial x^*}\frac{\partial p}{\partial \sigma} - 
\left(\frac{1}{D}\frac{\partial z_b}{\partial y}+\frac{1+\sigma}{D}\frac{\partial D}{\partial y^*}\right)\frac{1}{D}\frac{\partial D}{\partial y^*}\frac{\partial p}{\partial \sigma}
$$
In this part the term correspoing to the above terms are multiplied by the mass matrix.
*/

void GetInverseHeight(double *, double *, double, int);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int Np = (int)mxGetScalar(prhs[0]);
	int K = (int)mxGetScalar(prhs[1]);
	double Hcrit = mxGetScalar(prhs[2]);
	double *EToE = mxGetPr(prhs[3]);
	int Nface = (int)mxGetScalar(prhs[4]);
	double *Mass3d = mxGetPr(prhs[5]);
	double *J = mxGetPr(prhs[6]);

	double *StiffMatrix = mxGetPr(prhs[7]);
	mwIndex *IrStiffMatrix = mxGetIr(prhs[7]);
	mwIndex *JcStiffMatrix = mxGetJc(prhs[7]);
	int row = (int)mxGetM(prhs[7]);
	int col = (int)mxGetN(prhs[7]);

	double *PNPX = mxGetPr(prhs[8]);
	mwIndex *IrPNPX = mxGetIr(prhs[8]);
	mwIndex *JcPNPX = mxGetJc(prhs[8]);

	double *PNPY = mxGetPr(prhs[9]);
	mwIndex *IrPNPY = mxGetIr(prhs[9]);
	mwIndex *JcPNPY = mxGetJc(prhs[9]);

	double *PNPS = mxGetPr(prhs[10]);
	mwIndex *IrPNPS = mxGetIr(prhs[10]);
	mwIndex *JcPNPS = mxGetJc(prhs[10]);

	double *Height = mxGetPr(prhs[11]);
	double *Hx = mxGetPr(prhs[12]);
	double *Hy = mxGetPr(prhs[13]);
	double *Zx = mxGetPr(prhs[14]);
	double *Zy = mxGetPr(prhs[15]);

	double *Z = mxGetPr(prhs[16]);

	double *UniEleNumber = mxGetPr(prhs[17]);
	double *UniEle = mxGetPr(prhs[18]);

	double *InverseHeight = malloc(Np*K * sizeof(double));

	plhs[0] = mxCreateSparse(row, col, JcStiffMatrix[col], mxREAL);
	double *sr = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);
	mwIndex *jcs = mxGetJc(plhs[0]);
	/*Copy the input stiff matrix into the output one*/
	memcpy(sr, StiffMatrix, JcStiffMatrix[col] * sizeof(double));
	memcpy(irs, IrStiffMatrix, JcStiffMatrix[col] * sizeof(mwIndex));
	memcpy(jcs, JcStiffMatrix, (col+1) * sizeof(mwIndex));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		GetInverseHeight(InverseHeight + Np*k, Height + Np*k, Hcrit, Np);
	}

	ptrdiff_t One = 1;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < K; e++) {
		int OutNum = (int)UniEleNumber[e];
		double *UniElement = UniEle + e * (Nface + 1);
		//FindUniqueElementAndSortOrder(UniElement, EToE + e*Nface, &OutNum, Nface, e + 1);
		double *EleMass3d = malloc(Np*Np*OutNum * sizeof(double));
		for (int i = 0; i < OutNum; i++) {
			DiagMultiply(EleMass3d+i*Np*Np,Mass3d,J+((int)UniElement[i]-1)*Np,Np);
		}
		double *TempStiffMatrix = malloc(Np * sizeof(double));
		/*The following space is used to store the coefficient in direction x*/
		double *Coex = malloc(Np * sizeof(double));
		/*The following space is used to store the coefficient in direction y*/
		double *Coey = malloc(Np * sizeof(double));
		double *TempCoex = malloc(Np * sizeof(double));
		double *TempCoey = malloc(Np * sizeof(double));

		for (int p = 0; p < Np; p++) {
			int StartPoint = JcPNPX[e*Np + p];
			int Nonzero = (int)JcPNPX[e*Np + p + 1] - StartPoint;
			/*How many elements are contained in the current row, we note that Np points are influenced for each cell for this special case*/
			/*This can be seen from the initialization stage. We note that less than Np points are actually influenced for facial integral part,
			To faciliate the assemble of the sparse matrix, space for the first order term is enlarged.*/
			int Ele = Nonzero / Np;
			for (int i = 0; i < Ele; i++) {
				/*The index of the influenced element*/
				int EI = (int)IrPNPX[StartPoint + i*Np] / Np + 1;
				/*The index of the studied element in UniElement, this is used when call the mass matrix*/
				int UEI;
				for (int j = 0; j < OutNum; j++) {
					if (EI == (int)UniElement[j]) {
						UEI = j;
						break;
					}
				}
				/*$\frac{1}{D}\frac{\partial D}{\partial x^*}$*/
				DotProduct(Coex, InverseHeight + (EI - 1)*Np, Hx + (EI - 1)*Np, Np);
				/*$\frac{1}{D}\frac{\partial D}{\partial y^*}$*/
				DotProduct(Coey, InverseHeight + (EI - 1)*Np, Hy + (EI - 1)*Np, Np);
				/*$\frac{1}{D}\frac{\partial D}{\partial x^*}\frac{\partial p}{\partial x}$*/
				DotProduct(Coex, Coex, PNPX + StartPoint + i*Np, Np);
				/*$\frac{1}{D}\frac{\partial D}{\partial y^*}\frac{\partial p}{\partial y}$*/
				DotProduct(Coey, Coey, PNPY + StartPoint + i*Np, Np);
				Add(Coex, Coex, Coey, Np);
				MatrixMultiply("N", "N", (ptrdiff_t)Np, One, (ptrdiff_t)Np, 1.0, EleMass3d + UEI*Np*Np, (ptrdiff_t)Np, Coex, \
					(ptrdiff_t)Np, 0.0, TempStiffMatrix, Np);
				int StiffStartPoint = jcs[e*Np + p];
				for (int rp = 0; rp < jcs[e*Np + p + 1] - StiffStartPoint; rp++) {
					if (irs[StiffStartPoint + rp] == (EI - 1)*Np) {
						for (int i = 0; i < Np; i++) {
							sr[StiffStartPoint + rp + i] += TempStiffMatrix[i];
						}
						break;
					}
				}
			}
		}

		for (int p = 0; p < Np; p++) {
			int StartPoint = JcPNPS[e*Np + p];
			int Nonzero = JcPNPS[e*Np + p + 1] - StartPoint;
			int Ele = Nonzero / Np;
			for (int i = 0; i < Ele; i++) {
				/*The index of the influenced element*/
				int EI = (int)IrPNPS[StartPoint + i*Np] / Np + 1;
				/*The index of the studied element in UniElement, this is used when call the mass matrix*/
				int UEI;
				for (int j = 0; j < OutNum; j++) {
					if (EI == (int)UniElement[j]) {
						UEI = j;
						break;
					}
				}
				/*$\frac{1}{D}\frac{\partial z_b}{\partial x}$*/
				DotProduct(TempCoex, InverseHeight + (EI - 1)*Np, Zx + (EI - 1)*Np, Np);
				/*$1+\sigma$*/
				AddByConstant(Coex, Z + (EI - 1)*Np, 1.0, Np);
				/*$\frac{1+\sigma}{D}$*/
				DotProduct(Coex, Coex, InverseHeight + (EI - 1)*Np, Np);
				/*$\frac{1+\sigma}{D}\frac{\partial D}{\partial x}$*/
				DotProduct(Coex, Coex, Hx + (EI - 1)*Np, Np);
				/*$\frac{1}{D}\frac{\partial z_b}{\partial x} + \frac{1+\sigma}{D}\frac{\partial D}{\partial x}$*/
				Add(Coex, Coex, TempCoex, Np);
				/*$\left (\frac{1}{D}\frac{\partial z_b}{\partial x} + \frac{1+\sigma}{D}\frac{\partial D}{\partial x}\right)\frac{1}{D}$*/
				DotProduct(Coex, Coex, InverseHeight + (EI - 1)*Np, Np);
				/*$\left (\frac{1}{D}\frac{\partial z_b}{\partial x} + \frac{1+\sigma}{D}\frac{\partial D}{\partial x}\right)\frac{1}{D}\frac{\partial D}{\partial x}$*/
				DotProduct(Coex, Coex, Hx + (EI - 1)*Np, Np);

				/*$\frac{1}{D}\frac{\partial z_b}{\partial y}$*/
				DotProduct(TempCoey, InverseHeight + (EI - 1)*Np, Zy + (EI - 1)*Np, Np);
				/*$1+\sigma$*/
				AddByConstant(Coey, Z + (EI - 1)*Np, 1.0, Np);
				/*$\frac{1+\sigma}{D}$*/
				DotProduct(Coey, Coey, InverseHeight + (EI - 1)*Np, Np);
				/*$\frac{1+\sigma}{D}\frac{\partial D}{\partial y}$*/
				DotProduct(Coey, Coey, Hy + (EI - 1)*Np, Np);
				/*$\frac{1}{D}\frac{\partial z_b}{\partial y} + \frac{1+\sigma}{D}\frac{\partial D}{\partial y}$*/
				Add(Coey, Coey, TempCoey, Np);
				/*$\left (\frac{1}{D}\frac{\partial z_b}{\partial y} + \frac{1+\sigma}{D}\frac{\partial D}{\partial y}\right)\frac{1}{D}$*/
				DotProduct(Coey, Coey, InverseHeight + (EI - 1)*Np, Np);
				/*$\left (\frac{1}{D}\frac{\partial z_b}{\partial y} + \frac{1+\sigma}{D}\frac{\partial D}{\partial y}\right)\frac{1}{D}\frac{\partial D}{\partial y}$*/
				DotProduct(Coey, Coey, Hy + (EI - 1)*Np, Np);

				Add(Coex, Coex, Coey, Np);

				DotProduct(Coex, Coex, PNPS + StartPoint + i*Np, Np);

				MatrixMultiply("N", "N", (ptrdiff_t)Np, One, (ptrdiff_t)Np, 1.0, EleMass3d + UEI*Np*Np, (ptrdiff_t)Np, Coex, \
					(ptrdiff_t)Np, 0.0, TempStiffMatrix, Np);
				int StiffStartPoint = jcs[e*Np + p];
				for (int rp = 0; rp < jcs[e*Np + p + 1] - StiffStartPoint; rp++) {
					if (irs[StiffStartPoint + rp] == (EI - 1)*Np) {
						for (int i = 0; i < Np; i++) {
							sr[StiffStartPoint + rp + i] -= TempStiffMatrix[i];
						}
						break;
					}
				}
			}
		}
		free(Coex);
		free(Coey);
		free(TempCoex);
		free(TempCoey);
		free(TempStiffMatrix);
		free(EleMass3d);
//		free(UniElement);
	}
	free(InverseHeight);
}

void GetInverseHeight(double *dest, double *source, double Hcritical, int size) {
	for (int i = 0; i < size; i++) {
		if (source[i] >= Hcritical) {
			dest[i] = 1.0 / source[i];
		}
	}
}