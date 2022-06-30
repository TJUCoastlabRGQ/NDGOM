
#include  "mxImplicitVerticalEddyViscosity.h"

int *Ir = NULL, *Jc = NULL, *QuatrIr = NULL, *QuatrJc = NULL;

char *ImVertEddyInitialized = "False";

int NNZ, QuatrNNZ;

/*Function checked*/
void AssembleContributionIntoSparseMatrix(double *dest, double *src, int NonzeroNum, int Np) {
	double *Tempdest = dest;
	for (int col = 0; col < Np; col++) {
		Tempdest = dest + col * NonzeroNum;
		for (int row = 0; row < Np; row++) {
			Tempdest[row] += src[col*Np + row];
		}
	}
}

/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the column index, and has been checked*/
void AssembleContributionIntoColumn(double *dest, double *source, double *column, int Np3d, int Np2d)
{
	for (int colI = 0; colI < Np2d; colI++) {
		for (int RowI = 0; RowI < Np3d; RowI++)
			dest[((int)column[colI] - 1)*Np3d + RowI] += source[colI*Np3d + RowI];
	}
}
/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the row index, and has been checked*/
void AssembleContributionIntoRow(double *dest, double *source, double *Row, int Np3d, int Np2d)
{
	for (int colI = 0; colI < Np3d; colI++) {
		for (int RowI = 0; RowI < Np2d; RowI++)
			dest[colI*Np3d + (int)Row[RowI] - 1] += source[colI*Np2d + RowI];
	}
}
/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the row index and the column index, and has been checked*/
void AssembleContributionIntoRowAndColumn(double *dest, double *source, double *Row, double *column, int Np3d, int Np2d, int Flag)
{
	for (int colI = 0; colI < Np2d; colI++)
	{
		for (int RowI = 0; RowI < Np2d; RowI++)
			dest[((int)column[colI] - 1)*Np3d + (int)Row[RowI] - 1] += Flag * source[colI*Np2d + RowI];
	}
}

void AssembleDataIntoPoint(double *dest, double *source, double *PIndex, int Size) {
	for (int p = 0; p < Size; p++) {
		dest[(int)PIndex[p] - 1] += source[p];
	}
}

/*Note: this function is used to assemble the element mass matrix and the physical diff matrix, and has been verified.
Left multiply the matrix source with a diagonal matrix composed of element contained in coe.
*/
void DiagMultiply(double *dest, const double *source, const double *coe, int Np)
{
	for (int colI = 0; colI < Np; colI++) {
		for (int RowI = 0; RowI < Np; RowI++)
			dest[colI*Np + RowI] = coe[RowI] * source[colI*Np + RowI];
	}
}

void DiagLeftMultiplyUnsymmetric(double *dest, const double *source, const double *coe, int Row, int Col)
{
	for (int colI = 0; colI < Col; colI++) {
		for (int RowI = 0; RowI < Row; RowI++)
			dest[colI*Row + RowI] = coe[RowI] * source[colI*Row + RowI];
	}
}

void DiagRightMultiply(double *dest, const double *source, const double *coe, int Np)
{
	for (int colI = 0; colI < Np; colI++) {
		for (int RowI = 0; RowI < Np; RowI++)
			dest[colI*Np + RowI] = coe[colI] * source[colI*Np + RowI];
	}
}

void DiagRightMultiplyUnsymmetric(double *dest, const double *source, const double *coe, int Row, int Col)
{
	for (int colI = 0; colI < Col; colI++) {
		for (int RowI = 0; RowI < Row; RowI++)
			dest[colI*Row + RowI] = coe[colI] * source[colI*Row + RowI];
	}
}

void ImMatrixInverse(double *dest, lapack_int Np)
{
	lapack_int *IPIV = mkl_malloc(Np * sizeof(lapack_int), 64);

	LAPACKE_dgetrf(LAPACK_COL_MAJOR, Np, Np, dest, Np, IPIV);

	LAPACKE_dgetri(LAPACK_COL_MAJOR, Np, dest, Np, IPIV);

	mkl_free(IPIV);
}

void MultiplyByConstant(double *dest, double *Source, double Coefficient, int Np) {
	for (int i = 0; i < Np; i++)
		dest[i] = Source[i] * Coefficient;
}

void ImEddyVisInVertAllocation(int Np, int Nlayer) {

	NNZ = Np * Np * 3 * Nlayer - 2 * Np * Np;

	QuatrNNZ = 2 * NNZ + 2 * Np*Nlayer;

	Ir = malloc(NNZ * sizeof(int));

	QuatrIr = malloc(QuatrNNZ * sizeof(int));

	Jc = malloc((Nlayer * Np + 1)*sizeof(int));
	// The size has doubled
	QuatrJc = malloc((2*Nlayer * Np + 1) * sizeof(int));

	while (Ir == NULL) {
		Ir = malloc(NNZ * sizeof(int));
	}
	while (QuatrIr == NULL)
	{
		QuatrIr = malloc(QuatrNNZ * sizeof(int));
	}
	while (Jc == NULL) {
		Jc = malloc((Nlayer*Np+1) * sizeof(int));
	}
	while (QuatrJc == NULL) {
		QuatrJc = malloc((2*Nlayer * Np + 1) * sizeof(int));
	}

	int *SingleColumn = malloc(Np*Nlayer * sizeof(int));

	int *SingleRow;

	int SingleNonzero;

	if (Nlayer == 1) {
		SingleNonzero = Np*Np;
		// The row index of the non-zero elements
		SingleRow = malloc(SingleNonzero * sizeof(int));
		for (int ic = 0; ic < Np; ic++) {
			// Number of total non-zero elements
			SingleColumn[ic] = (ic + 1)*Np;
			for (int jr = 0; jr < Np; jr++) {
				SingleRow[ic*Np + jr] = jr;
			}
		}
	}
	else if (Nlayer == 2) {
		SingleNonzero = Np*Np * 2 * 2;
		SingleRow = malloc(SingleNonzero * sizeof(int));
		for (int ic = 0; ic < 2 * Np; ic++) {
			SingleColumn[ic] = (ic + 1)*Np * 2;
			for (int jr = 0; jr < 2 * Np; jr++) {
				SingleRow[ic*Np * 2 + jr] = jr;
			}
		}
	}
	else {
		SingleNonzero = Np*Np * 2 * 2 + (Nlayer - 2)*Np*Np * 3;
		SingleRow = malloc(SingleNonzero * sizeof(int));
		int NSumR = 0;
		for (int Layer = 0; Layer < 1; Layer++) {
			for (int ic = Layer*Np; ic < (Layer + 1)*Np; ic++) {
				SingleColumn[ic] = (ic + 1) * 2 * Np;
				for (int jr = 0; jr < 2 * Np; jr++) {
					SingleRow[NSumR + ic*Np * 2 + jr] = jr;
				}
			}
		}
		NSumR = Np * Np * 2;
		for (int Layer = 1; Layer < Nlayer - 1; Layer++) {
			for (int ic = Layer*Np; ic < (Layer + 1)*Np; ic++) {
				SingleColumn[ic] = SingleColumn[ic - 1] + 3 * Np;
				for (int jr = 0; jr < 3 * Np; jr++) {
					SingleRow[NSumR + (ic - Np) * 3 * Np + jr] = (Layer - 1)*Np + jr;
				}
			}
		}
		NSumR = Np*Np * 2 + (Nlayer - 2)*Np * 3 * Np;
		for (int Layer = Nlayer - 1; Layer < Nlayer; Layer++) {
			for (int ic = Layer*Np; ic < (Layer + 1)*Np; ic++) {
				SingleColumn[ic] = SingleColumn[ic - 1] + 2 * Np;
				for (int jr = 0; jr < 2 * Np; jr++) {
					SingleRow[NSumR + (ic - (Nlayer - 1)*Np) * 2 * Np + jr] = (Layer - 1)*Np + jr;   //²âÊÔ
				}
			}
		}
	}

	/*Note data are indexed begin with one in pardiso, so we have to add Ir and Jc*/
	for (int j = 0; j < SingleNonzero; j++) {
		Ir[j] = SingleRow[j] + 1;
	}
	Jc[0] = 1;
	for (int j = 0; j < Nlayer*Np; j++) {
		Jc[1 + j] = SingleColumn[j] + 1;
	}

	for (int p = 0; p < Nlayer*Np; p++) {
		if (p == 0) {
			for (int p1 = 0; p1 < SingleColumn[0]; p1++) {
				QuatrIr[p1] = SingleRow[p1] + 1;
			}
			QuatrIr[SingleColumn[0]] = Nlayer*Np + 1;
		}
		else {
			for (int p1 = 0; p1 < SingleColumn[p] - SingleColumn[p - 1]; p1++) {
				QuatrIr[SingleColumn[p - 1] + p + p1] = SingleRow[SingleColumn[p - 1] + p1]+1;
			}
			QuatrIr[SingleColumn[p] + p] = Nlayer*Np + p + 1;
		}
		
	}
	for (int p = Nlayer*Np; p < 2 * Nlayer*Np; p++) {
		if (p == Nlayer*Np) {
			QuatrIr[SingleNonzero + Nlayer*Np + p - Nlayer*Np] = p - Nlayer*Np + 1;
			for (int p1 = 0; p1 < SingleColumn[p- Nlayer*Np]; p1++) {
				//+1 for the upper diagnoal element
				QuatrIr[SingleNonzero + Nlayer*Np + p1+1] = SingleRow[p1] + Nlayer*Np + 1;
			}
		}
		else {
			QuatrIr[SingleNonzero + Nlayer*Np + SingleColumn[p-Nlayer*Np-1] + p - Nlayer*Np] = p - Nlayer*Np + 1;
			for (int p1 = 0; p1 < SingleColumn[p- Nlayer*Np] - SingleColumn[p - Nlayer*Np - 1]; p1++) {
				QuatrIr[SingleNonzero + Nlayer*Np + SingleColumn[p - Nlayer*Np-1] + p - Nlayer*Np + p1 + 1] = \
					SingleRow[SingleColumn[p - 1 - Nlayer*Np] + p1]+ Nlayer*Np + 1;
			}
		}
	}
	QuatrJc[0] = 1;
	for (int j = 0; j < Nlayer*Np; j++) {
		QuatrJc[1 + j] = SingleColumn[j] + 1 + j + 1;
	}
	for (int j = Nlayer*Np; j < 2*Nlayer*Np; j++) {
		QuatrJc[1 + j] = SingleColumn[Nlayer*Np-1] + Nlayer*Np + SingleColumn[j- Nlayer*Np] + 1 + j - Nlayer*Np + 1;
	}

   /*Set the thread for pardiso part to be one, set this part for once*/
	int thread = 1;
	mkl_set_num_threads(thread);
	ImVertEddyInitialized = "True";
	free(SingleColumn);
	free(SingleRow);
}



void ImEddyVisInVertDeAllocation() 
{
	free(Ir), Ir = NULL;
	free(QuatrIr), QuatrIr = NULL;
	free(Jc), Jc = NULL;
	free(QuatrJc), QuatrJc = NULL;
	ImVertEddyInitialized = "False";
}



// n is the leading dimension of the stiff matrix
void SparseEquationSolve(double *dest, MKL_INT n, double *StiffMatrix, double *RHS, int Nlayer, int Np, \
	int Index, int Nvar, int K3d, int Nonzero, int Column, int *RowIndex, int *ColumnNonzero) {

	if (Nvar >= 1) {
		int *TempIr = malloc(Nonzero * sizeof(int));
		memcpy(TempIr, RowIndex, Nonzero * sizeof(int));

		int *TempJc = malloc((Column + 1) * sizeof(int));
		memcpy(TempJc, ColumnNonzero, (Column + 1) * sizeof(int));

		MKL_INT iparm[64];
		void *pt[64];
		MKL_INT i;

		for (i = 0; i < 64; i++)
		{
			iparm[i] = 0;
		}
		iparm[0] = 1;         /* No solver default */
		iparm[1] = 2;         /* Fill-in reordering from METIS */
		iparm[3] = 0;         /* No iterative-direct algorithm */
		iparm[4] = 0;         /* No user fill-in reducing permutation */
		iparm[5] = 0;         /* Write solution into x */
		iparm[6] = 0;         /* Not in use */
		iparm[7] = 2;         /* Max numbers of iterative refinement steps */
		iparm[8] = 0;         /* Not in use */
		iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
		iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
		iparm[11] = 0;        /* Conjugate/transpose solve */
		iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
		iparm[13] = 0;        /* Output: Number of perturbed pivots */
		iparm[14] = 0;        /* Not in use */
		iparm[15] = 0;        /* Not in use */
		iparm[16] = 0;        /* Not in use */
		iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
		iparm[18] = -1;       /* Output: Mflops for LU factorization */
		iparm[19] = 0;        /* Output: Numbers of CG Iterations */

		for (i = 0; i < 64; i++)
		{
			pt[i] = 0;
		}

		MKL_INT mtype = 11;       /* Real unsymmetric matrix */

		MKL_INT maxfct, mnum, error, msglvl;
		MKL_INT Tempmaxfct, Tempmnum, Temperror, Tempmsglvl, Tempidum;
		double Tempddum;
		maxfct = 1;           /* Maximum number of numerical factorizations. */
		mnum = 1;             /* Which factorization to use. */
		msglvl = 0;           /* Print statistical information  */
		error = 0;            /* Initialize error flag */

		double ddum;          /* Double dummy */

		MKL_INT idum;         /* Integer dummy. */

		PardisoReOrderAndSymFactorize(iparm, pt, n, StiffMatrix, TempJc, TempIr, \
			&maxfct, &mnum, &msglvl, &error, &mtype, &ddum, &idum);

		for (int var = 0; var < Nvar; var++) {

			Tempmaxfct = maxfct, Tempmnum = mnum, \
				Temperror = error, Tempmsglvl = msglvl, Tempidum = idum, Tempddum = ddum;

			PardisoFactorize(iparm, pt, n, StiffMatrix + var*Nonzero, TempJc, TempIr, \
				&Tempmaxfct, &Tempmnum, &Tempmsglvl, &Temperror, &mtype, &Tempddum, &Tempidum);

			PardisoSolve(dest + var * Np*K3d + Index*Nlayer*Np, StiffMatrix + var*Nonzero, \
				RHS + var * Np*K3d + Index*Nlayer*Np, n, iparm, TempJc, TempIr, \
				&Tempmaxfct, &Tempmnum, &Tempmsglvl, &Temperror, &mtype, pt, &Tempddum, &Tempidum);
		}

		PardisoFree(iparm, pt, n, StiffMatrix, TempJc, TempIr, \
			&Tempmaxfct, &Tempmnum, &Tempmsglvl, &Temperror, &mtype, &Tempddum, &Tempidum);

		free(TempIr);
		free(TempJc);
	}
}

void PardisoFree(MKL_INT *iparm, void *pt, MKL_INT n, double *a, MKL_INT *ia, MKL_INT *ja, \
	MKL_INT *maxfct, MKL_INT *mnum, MKL_INT *msglvl, MKL_INT *error, MKL_INT *mtype, double *ddum, MKL_INT *idum) {
	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	MKL_INT NRHS = 1;
	MKL_INT phase = -1;           /* Release internal memory. */
	PARDISO(pt, maxfct, mnum, mtype, &phase,
		&n, ddum, ia, ja, idum, &NRHS,
		iparm, msglvl, ddum, ddum, error);
}

void PardisoReOrderAndSymFactorize(MKL_INT *iparm, void *pt, MKL_INT n, double *a, MKL_INT *ia, MKL_INT *ja, \
	MKL_INT *maxfct, MKL_INT *mnum, MKL_INT *msglvl, MKL_INT *error, MKL_INT *mtype, double *ddum, MKL_INT *idum) {

	MKL_INT NRHS = 1;

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	MKL_INT phase;

	phase = 11;
	PARDISO(pt, maxfct, mnum, mtype, &phase,
		&n, a, ia, ja, idum, &NRHS, iparm, msglvl, ddum, ddum, error);

}

// Here ia corresponds to Jc, ja corresponds to Ir
void PardisoFactorize(MKL_INT *iparm, void *pt, MKL_INT n, double *a, MKL_INT *ia, MKL_INT *ja, \
	MKL_INT *maxfct, MKL_INT *mnum, MKL_INT *msglvl, MKL_INT *error, MKL_INT *mtype, double *ddum, MKL_INT *idum) {

	MKL_INT NRHS = 1;

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	MKL_INT phase;

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, maxfct, mnum, mtype, &phase,
		&n, a, ia, ja, idum, &NRHS, iparm, msglvl, ddum, ddum, error);
}

void PardisoSolve(double *dest, double *StiffMatrix, double *RHS, MKL_INT n, int *iparm, \
	MKL_INT *ia, MKL_INT *ja, MKL_INT *maxfct, MKL_INT *mnum, MKL_INT *msglvl, \
	MKL_INT *error, MKL_INT *mtype, void *pt, double *ddum, MKL_INT *idum)
{
	/* -------------------------------------------------------------------- */
	/* .. Solution phase. */
	/* -------------------------------------------------------------------- */
	MKL_INT phase = 33;

	MKL_INT NRHS = 1;

	// Transpose solve is used for systems in CSC format
	iparm[11] = 2;

	PARDISO(pt, maxfct, mnum, mtype, &phase,
		&n, StiffMatrix, ia, ja, idum, &NRHS, iparm, msglvl, RHS, dest, error);
}

void SparseMatrixMultiply(double *dest, double *OPA, double *B, MKL_INT ROWA, MKL_INT *ia, MKL_INT *ja) {

	const double alpha = 1.0;

	const double beta = 0.0;

	const sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;

	sparse_matrix_t A;

	mkl_sparse_d_create_csc(&A, SPARSE_INDEX_BASE_ONE, ROWA, ROWA, ia, ia + 1, ja, OPA);

	struct matrix_descr descr;

	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	//	descr.mode = SPARSE_FILL_MODE_LOWER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	const sparse_layout_t layout = SPARSE_LAYOUT_COLUMN_MAJOR;
	/*Number of columns of matrix C.*/
	const MKL_INT columns = 1;

	const MKL_INT ldc = ROWA;

	mkl_sparse_d_mm(operation, alpha, A, descr, layout, B, columns, ROWA, beta, dest, ldc);

	mkl_sparse_destroy(A);
}