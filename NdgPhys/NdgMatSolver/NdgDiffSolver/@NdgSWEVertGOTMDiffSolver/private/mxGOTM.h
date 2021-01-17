#ifndef MXGOTM_H
#define MXGOTM_H

#include "mex.h"


extern double *tkeGOTM, *epsGOTM, *LGOTM, *nuhGOTM, *numGOTM, *layerHeight, *huCentralDate, \
*hvCentralDate, *huVerticalLine, *hvVerticalLine, *shearFrequencyDate, *buoyanceFrequencyDate, *BottomFrictionLength, \
*BottomFrictionVelocity, *SurfaceFrictionLength, *SurfaceFrictionVelocity, *eddyViscosityDate;

extern signed char *GOTMInitialized;

/*This function is used to interpolate the physical value from interpolation point to central point in vertical direction*/
void InterpolationToCentralPoint(double *, double *, ptrdiff_t *, ptrdiff_t *, ptrdiff_t *, double *);
/*This function is used to map the date located at the central point to the vertical line that the GOTM adapted*/
void mapCentralPointDateToVerticalDate(double *, double *, int , long long int , int );
/*This function is used to calculate the shear production term*/
void CalculateShearFrequencyDate(double *, int , int , double , long long int );

void CalculateBuoyanceFrequencyDate(int , int , long long int );
// CalculateBuoyanceProductionTerm to be added 
//Here, z0b is the bottom roughness, utaub is the friction velocity, z0s is the surface roughness
void CalculateLengthScaleAndShearVelocity(double *, double, double *, double *, double *, int, int, int, double, long long int);/*This is the interface to do_turbulence function in GOTM*/
//void DGDoTurbulence(double *dt, double *H2d, double *ShearProductionDate, double *buoyanceProductionDate, double *utaus, double *utaub, double *z0s, double *z0b);
void DGDoTurbulence(double *, double *, double , double *, int, int, long long int);
/*This function is used to map the date calculated by GOTM to the output matrix*/
void mapVedgeDateToDof(double *, double *, int, int, int, long long int);
/*This function is used to initialize the GOTM module*/
void InitTurbulenceModelGOTM(long long int *, char *, long long int, long long int, int, int);
/*This function is used to calculate the water depth for each layer*/
void CalculateWaterDepth(double *, int, int, double, long long int);

void getGotmDate(int, long long int);

void setGotmDate(int, long long int);

/*The following is the GOTM part*/
/*
void __stdcall __turbulence_MOD_init_turbulence(long long int *, char *, long long int *, long long int );

void __stdcall __mtridiagonal_MOD_init_tridiagonal(long long int *);

void __stdcall __turbulence_MOD_do_turbulence(long long int *, double *, double *, double *,
	double *, double *, double *, double *, double *, double *, double *);
__turbulence_MOD_do_turbulence
void __stdcall __mtridiagonal_MOD_clean_tridiagonal();

void __stdcall __turbulence_MOD_clean_turbulence();
*/
void __turbulence_MOD_init_turbulence(long long int *, char *, long long int *, long long int );

void __mtridiagonal_MOD_init_tridiagonal(long long int *);

void __turbulence_MOD_do_turbulence(long long int *, double *, double *, double *,
	double *, double *, double *, double *, double *, double *, double *);
void 

__mtridiagonal_MOD_clean_tridiagonal();

void __turbulence_MOD_clean_turbulence();

double * __turbulence_MOD_tke;
double * __turbulence_MOD_eps;
double * __turbulence_MOD_num;
double * __turbulence_MOD_nuh;
double * __turbulence_MOD_l;

#endif
