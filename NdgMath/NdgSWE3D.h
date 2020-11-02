#ifndef _NDGSWE3D_H
#define _NDGSWE3D_H

#define EPS 1e-6

void EvaluateVerticalFaceSurfFlux(double *dest, double *fm, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne);

void EvaluatePhysicalVariableByDepthThreshold(double hmin, double *h, double *variable, double *outPut);

void EvaluateVerticalFaceNumFlux(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne);

/** Rotate flux to outward normal direction */
void RotateFluxToNormal2d(double *hu, double *hv, double *nx, double *ny, double *qn, double *qv);

void EvaluateFluxTerm2d(double hmin, double *gra, double *h, double *hu, double *hv, double *E);

#endif