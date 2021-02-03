#ifndef _NDGSWE_H
#define _NDGSWE_H

#include <mex.h>
#include <math.h>

#define TOLERR 1e-6

#define EPS 1e-6

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

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

void ImposeBoundaryCondition(double *, NdgEdgeType , double *, double *, double *, double *, \
	double *, double *, double *, int , int , int );

void EvaluateHydroStaticReconstructValue(double , double *, double *, double *, double *, int , int , int );

//void EvaluateFlowRateByDeptheThreshold(double hmin, double *hM, double *huM, double *hvM, double *um, double *vm);

void GetPCENumericalFluxTerm_HLLC_LU(double *, double *, double *, double *, double *, double *, double , int , int );

void GetPCENumericalFluxTerm_HLL(double *, double *, double *, double *, double *, double *, double, int, int);

void GetPCENumericalFluxTerm_HLLC_LAI(double *, double *, double *, double *, double *, double *, double, int, int);

void EvaluateVerticalFaceNumFlux_HLLC_LAI(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne);

void EvaluateVerticalFaceNumFlux_HLL(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne);

void EvaluatePhysicalVariableByDepthThreshold(double , double *, double *, double *);

/** Rotate flux to outward normal direction */
void RotateFluxToNormal2d(double *, double *, double *, double *, double *, double *);

void EvaluateVerticalFaceRiemannProblem(double *, double *, double *, double *, double *, \
	double *, double, int, int, int, int, int);

#endif