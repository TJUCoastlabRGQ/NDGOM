#ifndef _NdgMath_H
#define _NdgMath_H

#include <mex.h>
#include "lapack.h"
#include <string.h>
#include <blas.h>
#include <stdlib.h>

void Add(double *, double *, double *, int );

void Minus(double *, double *, double *, int);

void DotProduct(double *, double *, double *, int );

void DotDivide(double *, double *, double *, int );

void DotDivideByConstant(double *, double *, double , int);

void DotCriticalDivide(double *, double *, double *, double *, int);

void StrongFormInnerEdgeRHS(int, double *, int, int, double *, double *, double *, double *, double *, double *, double *, double *);

void StrongFormBoundaryEdgeRHS(int , double *, int , int ,double *, double *, double *, double *, double *, double *);

void FetchInnerEdgeFacialValue(double *, double *, double *, double *, double *, double *, int , int );

void FetchBoundaryEdgeFacialValue(double *, double *, double *, double *, int , int );

void RepmatValue(double *, double *, int);


#endif