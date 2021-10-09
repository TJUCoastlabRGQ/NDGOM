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

void AssembleFacialContributionIntoSparseMatrix(double *, mwIndex *, mwIndex *, double *, double *, int, int, double *, int, int);

void AssembleVolumnContributionIntoSparseMatrix(double *, mwIndex *, mwIndex *, int, double *, int);

void AssembleFacialDiffMatrix(double *, double *, double *, int, int);

void AssembleFacialContributionForFirstOrderTermIntoSparseMatrix(double *, mwIndex *, mwIndex *, double *, int, int, double *, int, int);

void EvaluateNonhydroVerticalFaceSurfFlux(double *, double *, double *, int );

void EvaluateNonhydroVerticalFaceNumFlux_Central(double *, double *, double *, double *, int );

void GetSparsePattern(mwIndex *, mwIndex *, double *, double *, double *, double *, double *, \
	double *, double *, int, int, int, int, int, int , int );

void GetSparsePatternForHorizontalFirstOrderTerm(mwIndex *, mwIndex *, double *, double *, double *, double *, \
	int , int , int , int , int , int , double *);

void GetSparsePatternForVerticalFirstOrderTerm(mwIndex *, mwIndex *, double *, double *, double *, double *, \
	int , int , int , int , int , int , double *);

void FetchFacialData(double *, double *, double *, int );

void FetchDataInSparseMatrix(double *, double *, int , int );

void FindUniqueElementAndSortOrder(double *, double *, double *, double, int );

void FindFaceAndDirectionVector(double *, double *, double *, double *, double *, int , int , double *, double *, double *, int , int );

void FindGlobalBottomEdgeFace(int *, double *, int , int , int );

void FindFaceAndDirectionVectorAtBoundary(double *, int *, int *, int, int, double *, double *, int, double *, int);

void FindFaceAndFacialPoint(double *, double *, int, double *, double *, double *, int, int, int);

void CalculatePenaltyParameter(double *, double *, double *, double *, int , int , \
	int , double *, double *, double *, double *, double *, int , int );

#endif