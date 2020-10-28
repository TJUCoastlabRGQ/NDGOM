#ifndef _NDGSWE_H
#define _NDGSWE_H

#include <mex.h>
#include <math.h>

#define TOLERR 1e-6

typedef enum {
	NdgEdgeInner = 0,
	NdgEdgeGaussEdge = 1,
	NdgEdgeSlipWall = 2,
	NdgEdgeNonSlipWall = 3,
	NdgEdgeZeroGrad = 4,
	NdgEdgeClamped = 5,
	NdgEdgeClampedDepth = 6,
	NdgEdgeClampedVel = 7,
	NdgEdgeFlather = 8,
	NdgEdgeNonLinearFlather = 9,
	NdgEdgeNonLinearFlatherFlow = 10,
	NdgEdgeNonReflectingFlux = 11
} NdgEdgeType;

void ImposeBoundaryCondition(double *, NdgEdgeType , double *, double *, double *, double *, double *, double *, double *,\
	double *, double *, double *, double *, double *, double *);

void EvaluateHydroStaticReconstructValue(double hmin, double *, double *, double *, double *, double *, double *, double *,\
	double *);

void EvaluateFlowRateByDeptheThreshold(double hmin, double *hM, double *huM, double *hvM, double *um, double *vm);

void GetPCENumericalFluxTerm(double *, double *, double *, double *, double *, double *, double *, double *, double *, int, \
	double, double);

#endif