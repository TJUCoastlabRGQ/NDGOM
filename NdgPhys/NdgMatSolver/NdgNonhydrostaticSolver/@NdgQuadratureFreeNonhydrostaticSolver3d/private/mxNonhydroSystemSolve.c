#include <mex.h>

#include <string.h>

#include "mkl_pardiso.h"

#include "mkl_types.h"

#include "mkl_spblas.h"

#include "mkl.h"

int *NonhydroIr = NULL, *NonhydroJc = NULL;

double *NonhydroStiffMatrix = NULL;

MKL_INT iparm[64];

void *pt[64];

MKL_INT maxfct, mnum, error, msglvl;

MKL_INT Tempmaxfct, Tempmnum, Temperror, Tempmsglvl, Tempidum;

double Tempddum;

double ddum;          /* Double dummy */

MKL_INT idum;         /* Integer dummy. */

MKL_INT mtype;

int LeadDim;

char *NonhydroSystemInitialized = "False";

void NonhydroMemoryDeAllocation();

void NonhydroMemoryAllocationAndInitialize(mwIndex *, mwIndex *, int , int , int , double *);

void PardisoFactorize(MKL_INT *, void *, MKL_INT, double *, \
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, double *, MKL_INT *);

void PardisoReOrderAndSymFactorize(MKL_INT *, void *, MKL_INT, double *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, double *, MKL_INT *);

void PardisoFree(MKL_INT *, void *, MKL_INT, double *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, double *, MKL_INT *);

void PardisoSolve(double *, double *, double *, MKL_INT, int *, \
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, void *, double *, MKL_INT *);

void MyExit()
{
    if ( !strcmp("True", NonhydroSystemInitialized) ){
        NonhydroMemoryDeAllocation();
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexAtExit(&MyExit);

	double *StiffMatrix = mxGetPr(prhs[0]);
	mwIndex *TempIr = mxGetIr(prhs[0]);
	mwIndex *TempJc = mxGetJc(prhs[0]);
	int NonZero = (int)mxGetScalar(prhs[1]);
	int Np = (int)mxGetScalar(prhs[2]);
	int K = (int)mxGetScalar(prhs[3]);
	double *RHS = mxGetPr(prhs[4]);

	plhs[0] = mxCreateDoubleMatrix(Np*K, 1, mxREAL);
	double *NonhydroPressure = mxGetPr(plhs[0]);

	/*If not initialized, initialize first*/
	if (!strcmp("False", NonhydroSystemInitialized))
	{
		NonhydroMemoryAllocationAndInitialize(TempIr, TempJc, NonZero, Np, K, StiffMatrix);
	}
	// Copy the data from the input stiff matrix to the space allocated in this function
	memcpy(NonhydroStiffMatrix, StiffMatrix, NonZero * sizeof(double));

	PardisoSolve(NonhydroPressure, NonhydroStiffMatrix, \
		RHS, Np*K, iparm, NonhydroJc, NonhydroIr, \
		&Tempmaxfct, &Tempmnum, &Tempmsglvl, &Temperror, &mtype, pt, &Tempddum, &Tempidum);

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


void NonhydroMemoryAllocationAndInitialize(mwIndex *Ir, mwIndex *Jc, int NonZero, int Np, int K, double *StiffMatrix)
{
	/*Set the thread for pardiso part, set this part for once*/
	int thread = 6;
	mkl_set_num_threads(thread);

	int index;
	NonhydroIr = malloc(NonZero * sizeof(int));
	for (index = 0; index < NonZero; index++) 
	{
		NonhydroIr[index] = Ir[index] + 1;
	}
	NonhydroJc = malloc((Np*K + 1) * sizeof(int));
	for (index = 0; index < Np*K + 1; index++) 
	{
		NonhydroJc[index] = Jc[index] + 1;
	}

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

	mtype = 11;       /* Real unsymmetric matrix */

	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;             /* Which factorization to use. */
	msglvl = 0;           /* Print statistical information  */
	error = 0;            /* Initialize error flag */

	LeadDim = Np * K;

	NonhydroStiffMatrix = malloc(NonZero * sizeof(double));

	memcpy(NonhydroStiffMatrix, StiffMatrix, NonZero * sizeof(double));

	PardisoReOrderAndSymFactorize(iparm, pt, LeadDim, NonhydroStiffMatrix, NonhydroJc, NonhydroIr, \
		&maxfct, &mnum, &msglvl, &error, &mtype, &ddum, &idum);

	Tempmaxfct = maxfct, Tempmnum = mnum, \
		Temperror = error, Tempmsglvl = msglvl, Tempidum = idum, Tempddum = ddum;

	PardisoFactorize(iparm, pt, LeadDim, NonhydroStiffMatrix, NonhydroJc, NonhydroIr, \
		&Tempmaxfct, &Tempmnum, &Tempmsglvl, &Temperror, &mtype, &Tempddum, &Tempidum);

	NonhydroSystemInitialized = "True";

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

void NonhydroMemoryDeAllocation()
{
	PardisoFree(iparm, pt, LeadDim, NonhydroStiffMatrix, NonhydroJc, NonhydroIr, \
		&Tempmaxfct, &Tempmnum, &Tempmsglvl, &Temperror, &mtype, &Tempddum, &Tempidum);

	free(NonhydroIr); NonhydroIr = NULL;
	free(NonhydroJc); NonhydroJc = NULL;
	free(NonhydroStiffMatrix); NonhydroStiffMatrix = NULL;
	NonhydroSystemInitialized = "False";
}