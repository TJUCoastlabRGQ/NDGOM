#include <stdio.h>

#include <stdlib.h>

#include <string.h>

#include <mex.h>

#include <math.h>

#include "mkl_pardiso.h"

#include "mkl_types.h"

#include "mkl_spblas.h"

#include "mkl.h"

void ImEddyVisInVertAllocation( int, int);

void ImEddyVisInVertDeAllocation();

void PardisoFactorize(MKL_INT *, void *, MKL_INT , double *,\
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, double *, MKL_INT *);

void PardisoReOrderAndSymFactorize(MKL_INT *, void *, MKL_INT , double *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, double *, MKL_INT *);

void PardisoFree(MKL_INT *, void *, MKL_INT , double *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, double *, MKL_INT *);

void SparseEquationSolve(double *, MKL_INT , double *, double *, int ,\
	int , int , int , int , int , int , int *, int *);

void PardisoSolve(double *, double *, double *, MKL_INT , int *, \
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, void *, double *, MKL_INT *);

void AssembleContributionIntoSparseMatrix(double *, double *, int , int );

void SparseMatrixMultiply(double *, double *, double *, MKL_INT, MKL_INT *, MKL_INT *);

void MultiplyByConstant(double *dest, double *Source, double Coefficient, int Np);

/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the column index, and has been checked*/
void AssembleContributionIntoColumn(double *dest, double *source, double *column, int Np3d, int Np2d);
/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the row index, and has been checked*/
void AssembleContributionIntoRow(double *dest, double *source, double *Row, int Np3d, int Np2d);
/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the row index and the column index, and has been checked*/
void AssembleContributionIntoRowAndColumn(double *dest, double *source, double *Row, double *column, int Np3d, int Np2d, int Flag);

void AssembleDataIntoPoint(double *dest, double *source, double *PIndex, int Size);

void DiagMultiply(double *dest, const double *source, const double *coe, int Np);

void DiagLeftMultiplyUnsymmetric(double *dest, const double *source, const double *coe, int Row, int Col);

void DiagRightMultiply(double *dest, const double *source, const double *coe, int Np);

void DiagRightMultiplyUnsymmetric(double *dest, const double *source, const double *coe, int Row, int Col);

void ImMatrixInverse(double *, lapack_int );