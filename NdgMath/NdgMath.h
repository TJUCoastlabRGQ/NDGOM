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

void Add(double *, double *, double *, int );

void Minus(double *, double *, double *, int);

void DotProduct(double *, double *, double *, int );

void DotDivide(double *, double *, double *, int );

void DotDivideByConstant(double *, double *, double , int);

void DotCriticalDivide(double *, double *, double *, double *, int);

void StrongFormInnerEdgeRHS(int , double *, double *, int , int ,int , double *, double *, double *, double *,double *, double *, double *, double *);

void StrongFormBoundaryEdgeRHS(int , double *, double *, int , int ,  int , double *, double *, double *, double *, double *, double *);

void FetchInnerEdgeFacialValue(double *, double *, double *, double *, double *, double *, int , int );

void FetchBoundaryEdgeFacialValue(double *, double *, double *, double *, int , int );

void RepmatValue(double *, double *, int);

void AssembleContributionIntoColumn(double *, double *, double *, int, int);

void AssembleContributionIntoRow(double *, double *, double *, int, int);

void AssembleContributionIntoRowAndColumn(double *, double *, double *, double *, int, int, int);

void NdgExtend2dField(double *, double *, int, int, int, int, int);

void GetVolumnIntegral2d(double *, double *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, double *, \
	double *, ptrdiff_t *, double *, ptrdiff_t *, double *, ptrdiff_t *, double *, double *);

void GetVolumnIntegral3d(double *, double *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, \
	double *, double *, double *, double *, double *, double *, ptrdiff_t *, ptrdiff_t *, double *, ptrdiff_t *, \
	double *, double *, double *, double *, double *, int , int , int );

void GetFacialFluxTerm2d(double *, double *, double *, double *, double *, int);

void MultiEdgeContributionByLiftOperator(double *, double *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *, \
	double *, ptrdiff_t *, ptrdiff_t *, double *, ptrdiff_t *, double *, int);

#endif