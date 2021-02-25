#ifndef _SWENONHYDR3D
#define _SWENONHYDR3D

#include "mex.h"

#include <stdio.h>

#include <stdlib.h>

#include <string.h>

#include "../../../../../NdgMath/NdgMath.h"

#include "../../../../../NdgMath/NdgSWE.h"

#include "../../../../../NdgMath/NdgSWE3D.h"

#include "../../../../../NdgMath/NdgMemory.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void AssembleContributionIntoSparseMatrix(double *dest, double *src, int NonzeroNum, int Np);

void AssembleFacialDiffMatrix(double *, double *, double *, int, int);

void EvaluateNonhydroVerticalFaceSurfFlux(double *, double *, double *, int );

void EvaluateNonhydroVerticalFaceNumFlux_Central(double *, double *, double *, double *, int );

void GetSparsePatternInVerticalDirection(mwIndex *, mwIndex *, int, int, int);

void GetSparsePatternInHorizontalDirection(mwIndex *, mwIndex *, double *, int, int, int, int);

void FetchFacialData(double *, double *, double *, int );

void FetchDataInSparseMatrix(double *, double *, int , int );

void FindUniqueElementAndSortOrder(double *, double *, int *, int, int );

void FindFaceAndDirectionVector(double *, int *, int *, int *, int, int, int, double *, double *, double *, int);

void FindFaceAndDirectionVectorAtBoundary(double *, int *, int *, int, int, double *, double *, int, double *, int);

#endif