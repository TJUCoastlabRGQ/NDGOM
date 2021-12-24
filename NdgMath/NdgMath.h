#ifndef _NdgMath_H
#define _NdgMath_H

#ifdef _OPENMP
#include <omp.h>
#endif


#include <mex.h>
#include "matrix.h"
#include "lapack.h"
#include <string.h>
#include <blas.h>
#include <stdlib.h>
#include <math.h>

#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#if !defined(_WIN32)
#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)
#endif

/*Note: this function is used to assemble the element mass matrix and the physical diff matrix, and has been verified*/
void DiagMultiply(double *, const double *, const double *, int );

void DiagRightMultiply(double *, const double *, const double *, int );

void Add(double *, double *, double *, int );

void AddByConstant(double *, double *, double, int);

void AssembleDataIntoPoint(double *, double *, double *, int);

void Minus(double *, double *, double *, int);

void MatrixMultiply(char*, char*, ptrdiff_t, ptrdiff_t, ptrdiff_t, double, double *, ptrdiff_t, \
	double *, ptrdiff_t, double, double *, ptrdiff_t);

void MassLumping(double *, double *, int);

void MatrixInverse(double *, ptrdiff_t );

void DotProduct(double *, double *, double *, int );

void DotDivide(double *, double *, double *, int );

void DotDivideByConstant(double *, double *, double , int);

void MultiplyByConstant(double *, double *, double, int);

void DotCriticalDivide(double *, double *, double *, double *, int);

int Sign(double *);

void Sort(double *, int);

void StrongFormInnerEdgeRHS(int , double *, double *, int , int ,int , double *, double *, double *, \
	double *,double *, double *, double *, double *);

void StrongFormBoundaryEdgeRHS(int , double *, double *, int , int ,  int , double *, double *, double *, \
	double *, double *, double *);

void FetchInnerEdgeFacialValue(double *, double *, double *, double *, double *, double *, int , int );

void FetchBoundaryEdgeFacialValue(double *, double *, double *, double *, int , int );

void Flip(double *, int);

void ReverseValue(double *, double *, int );

void RepmatValue(double *, double *, int , int );

void AssembleContributionIntoColumn(double *, double *, double *, int, int);

void AssembleContributionIntoRow(double *, double *, double *, int, int);

void AssembleContributionIntoRowAndColumn(double *, double *, double *, double *, int, int, int);

void NdgExtend2dField(double *, double *, int, int, int, int, int);

void GetVolumnIntegral1d(double *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, \
	double *, ptrdiff_t *, double *, ptrdiff_t *, double *, ptrdiff_t *, double *);

void GetVolumnIntegral2d(double *, double *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, double *, \
	double *, ptrdiff_t *, double *, ptrdiff_t *, double *, ptrdiff_t *, double *, double *);

void GetVolumnIntegral3d(double *, double *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, \
	double *, double *, double *, double *, double *, double *, ptrdiff_t *, ptrdiff_t *, double *, ptrdiff_t *, \
	double *, double *, double *, double *, double *, int , int , int );

void GetFacialFluxTerm2d(double *, double *, double *, double *, double *, int);

void GetScalarFluxTerm2d(double *, double *, double *, int);

void MultiEdgeContributionByLiftOperator(double *, double *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, \
	double *, ptrdiff_t *, ptrdiff_t *, double *, ptrdiff_t *, double *, int);

void GetMeshIntegralValue(double *, char *, char *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, double *, \
	ptrdiff_t *, double *, double *, ptrdiff_t *, double *, ptrdiff_t *, double *);

void GetMeshAverageValue(double *, double *, char *, char *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, double *, \
	ptrdiff_t *, double *, double *, ptrdiff_t *, double *, ptrdiff_t *, double *);

#endif