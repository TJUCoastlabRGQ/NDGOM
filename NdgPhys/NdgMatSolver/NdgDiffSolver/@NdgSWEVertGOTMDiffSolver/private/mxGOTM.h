#ifndef MXGOTM_H
#define MXGOTM_H

#include "mex.h"
/*number of elements and interpolation points for the two-dimensional mesh*/
extern int K2d, Np2d;
/*number of elements and interpolation points for the three-dimensional mesh*/
extern int K3d, Np3d;
/*layers of the three dimensional shallow water model*/
extern long long int nlev;
/*the critical water depth, when the depth is less than this value, velosities are set to be zero*/
extern double hcrit;
/*the physical roughness length scale*/
extern double h0b;
/*final time of the */
extern double finalTime;

extern double *VCV, *tkeGOTM, *epsGOTM, *LGOTM, *nuhGOTM, *numGOTM, *layerHeight, *huCentralDate, \
*hvCentralDate, *huVerticalLine, *hvVerticalLine, *shearFrequencyDate, *buoyanceFrequencyDate, *BottomFrictionLength, \
*BottomFrictionVelocity, *SurfaceFrictionLength, *SurfaceFrictionVelocity, *eddyViscosityDate;

/*This function is used to interpolate the physical value from interpolation point to central point in vertical direction*/
void InterpolationToCentralPoint(double *fphys, double *dest, ptrdiff_t *Np2d, ptrdiff_t *K3d, ptrdiff_t *Np3d);
/*This function is used to map the date located at the central point to the vertical line that the GOTM adapted*/
void mapCentralPointDateToVerticalDate(double *CentralPointData, double *VerticalLineTempdata);
/*This function is used to calculate the shear production term*/
void CalculateShearFrequencyDate(double *H2d);

void CalculateBuoyanceFrequencyDate();
// CalculateBuoyanceProductionTerm to be added 
//Here, z0b is the bottom roughness, utaub is the friction velocity, z0s is the surface roughness
void CalculateLengthScaleAndShearVelocity(double *H2d, double *DragCoefficient, double *Taux, double *Tauy );
/*This is the interface to do_turbulence function in GOTM*/
//void DGDoTurbulence(double *dt, double *H2d, double *ShearProductionDate, double *buoyanceProductionDate, double *utaus, double *utaub, double *z0s, double *z0b);
void DGDoTurbulence(double *TimeStep, double *H2d, double *Grass);
/*This function is used to map the date calculated by GOTM to the output matrix*/
void mapVedgeDateToDof(double *VedgeDate, double *exportDofContainer);
/*This function is used to initialize the GOTM module*/
 void InitTurbulenceModelGOTM(long long int *, char *, long long int);
/*This function is used to calculate the water depth for each layer*/
 void CalculateWaterDepth(double *source);

 void getGotmDate(int);

 void setGotmDate(int);

/*The following is the GOTM part*/

void __stdcall TURBULENCE_mp_INIT_TURBULENCE(long long int *, char *, long long int *, long long int );

void __stdcall MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(long long int *);

void __stdcall TURBULENCE_mp_DO_TURBULENCE(long long int *, double *, double *, double *,
	double *, double *, double *, double *, double *, double *, double *);

void __stdcall MTRIDIAGONAL_mp_CLEAN_TRIDIAGONAL();

void __stdcall TURBULENCE_mp_CLEAN_TURBULENCE();

double * TURBULENCE_mp_TKE;
double * TURBULENCE_mp_EPS;
double * TURBULENCE_mp_NUM;
double * TURBULENCE_mp_NUH;
double * TURBULENCE_mp_L;

#endif