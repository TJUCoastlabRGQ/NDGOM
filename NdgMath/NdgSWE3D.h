#ifndef _NDGSWE3D_H
#define _NDGSWE3D_H

#include "mex.h"

void EvaluateVerticalFaceSurfFlux(double *dest, double *fm, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne);

void EvaluateHorizontalFaceSurfFlux(double *flux, double *fm, double *nz, double Hcrit, int Nfp, int Nvar, int Ne);

void EvaluatePhysicalVariableByDepthThreshold(double hmin, double *h, double *variable, double *outPut);

void EvaluateHorizontalFaceNumFlux(double *FluxS, double *fm, double *fp, double *nz, double Hcrit, int Nfp, int Nvar, int Ne);

/** Rotate flux to outward normal direction */
void RotateFluxToNormal2d(double *hu, double *hv, double *nx, double *ny, double *qn, double *qv);

void EvaluateFluxTerm2d(double hmin, double *gra, double *h, double *hu, double *hv, double *E);

void EvaluatePrebalanceVolumeTerm(double *, double *, double *, double *, double *, int , double *, int , int , double );

void GetModCoefficient(double *, double *, double *, int , int );

void GetIntegralValue(double *, int, double *, double *);

void VerticalColumnIntegralField3d(double *, int, double *, double *, double *, double *, double *, double *, double *, int, int);

void VerticalFaceColumnIntegral(double *, double *, double *, double *, ptrdiff_t , double *, int , double *, ptrdiff_t , int );

void VerticalIntegralFromBottom(double *, double *, double *, double *, int, ptrdiff_t, double *, int, int, double *);

void VerticalIntegralFromSurface(double *, double *, double *, double *, int, ptrdiff_t, double *, int, int, double *);

void VerticalRepmatFacialValue(double *, double *, int , int , int , int , int );

#endif