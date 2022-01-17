#include <stdio.h>

#include <stdlib.h>

#include <string.h>

#include <mex.h>

#include <math.h>

#include "mkl_pardiso.h"

#include "mkl_types.h"

#include "mkl_spblas.h"


int *Ir = NULL, *Jc = NULL;

char *VertEddyInitialized = "False";

int NNZ;

void ImEddyVisInVertAllocation( int, int, int);

void ImEddyVisInVertDeAllocation();

void PardisoFactorize(MKL_INT *, void *, MKL_INT , double *,\
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, double *, MKL_INT *);

void EquationSolve(double *, MKL_INT , double *, double *);

void PardisoSolve(double *, double *, double *, MKL_INT , int *, \
	MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, \
	MKL_INT *, MKL_INT *, void *, double *, MKL_INT *);