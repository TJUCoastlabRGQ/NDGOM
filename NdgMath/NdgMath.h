#ifndef _NdgMath_H
#define _NdgMath_H

#include <mex.h>
#include "lapack.h"
#include <string.h>
#include <blas.h>
#include <stdlib.h>

void DotProduct(double *, double *, int );

void DotDivide(double *, double *, double *, int );

void StrongFormInnerEdgeRHS(int, double *, int, int, double *, double *, double *, double *, double *, double *, double *, double *);

void StrongFormBoundaryEdgeRHS(int , double *, int , int ,double *, double *, double *, double *, double *, double *);

#endif