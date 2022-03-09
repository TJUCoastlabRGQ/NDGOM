#ifndef MXGOTM_H
#define MXGOTM_H

#include "mex.h"

#if !(defined(_WIN32) || defined(_WIN64))
#define TURBULENCE_mp_INIT_TURBULENCE __turbulence_MOD_init_turbulence
#define TURBULENCE_mp_DO_TURBULENCE __turbulence_MOD_do_turbulence
#define TURBULENCE_mp_CLEAN_TURBULENCE __turbulence_MOD_clean_turbulence
#define MTRIDIAGONAL_mp_INIT_TRIDIAGONAL __mtridiagonal_MOD_init_tridiagonal
#define MTRIDIAGONAL_mp_CLEAN_TRIDIAGONAL __mtridiagonal_MOD_clean_tridiagonal
#define TURBULENCE_mp_TKE __turbulence_MOD_tke
#define TURBULENCE_mp_EPS __turbulence_MOD_eps
#define TURBULENCE_mp_NUM __turbulence_MOD_num
#define TURBULENCE_mp_NUH __turbulence_MOD_nuh
#define TURBULENCE_mp_L __turbulence_MOD_l
#endif

#if !defined(_WIN32)
#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)
#endif

extern double *tkeGOTM, *epsGOTM, *LGOTM, *nuhGOTM, *numGOTM, *layerHeight, *huCentralDate, *huCentralDateNew,\
*hvCentralDate, *hvCentralDateNew, *huVerticalLine, *huVerticalLineNew, *hvVerticalLine, *hvVerticalLineNew, \
*shearFrequencyDate, *buoyanceFrequencyDate, *BottomFrictionLength, *BottomFrictionVelocity, \
*SurfaceFrictionLength, *SurfaceFrictionVelocity, *eddyViscosityDate, *rhoCentralDate, \
*rhoVerticalLine, *hcenter;

//*eddyDiffusionDate, *eddyTKEDate, *eddyLengthDate, *eddyEPSDate,

extern char *GOTMInitialized;

/*This function is used to interpolate the physical value from interpolation point to central point in vertical direction*/
void InterpolationToCentralPoint(double *fphys, double *dest, int K2d, int Np2d, int Np3d, \
	int nlayer, double *J2d, double *wq2d, double *Vq2d, ptrdiff_t RVq2d, ptrdiff_t Cvq2d, double *LAV);/*This function is used to map the date located at the central point to the vertical line that the GOTM adapted*/

void mapCentralPointDateToVerticalDate(double *centralDate, double *verticalLineDate, int K2d, \
	int nlev);
/*This function is used to calculate the shear production term*/
void CalculateShearFrequencyDate(int K2d, double hcrit, int nlev);

//void CalculateBuoyanceFrequencyDate(int , int , long long int );
void CalculateBuoyanceFrequencyDate(int Np2d, int K2d, double hcrit, int nlev, \
	double gra, double rho0);
// CalculateBuoyanceProductionTerm to be added 
//Here, z0b is the bottom roughness, utaub is the friction velocity, z0s is the surface roughness
void CalculateLengthScaleAndShearVelocity(double z0b, double z0s, double hcrit, double *DragCoefficient, \
	double *Taux, double *Tauy, int Np2d, int K2d, int nlev);

void DGDoTurbulence(double *TimeStep, double hcrit, double *Grass, int K2d, long long int nlev);
/*This function is used to map the date calculated by GOTM to the output matrix*/
void mapVedgeDateToDof(double *SourceDate, double *DestinationDate, int Np2d, int K2d, int Np3d, int nlev);
/*This function is used to initialize the GOTM module*/
void InitTurbulenceModelGOTM(long long int *NameList, char * buf, long long int buflen, \
	long long int nlev, int K2d);
/*This function is used to calculate the water depth for each layer*/
void CalculateWaterDepth(int K2d, double hcrit, int nlev);

void getGotmDate(int, int);

void setGotmDate(int, int);

/*The following is the GOTM part*/

void TURBULENCE_mp_INIT_TURBULENCE(long long int *, char *, long long int *, long long int );

void MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(long long int *);

void TURBULENCE_mp_DO_TURBULENCE(long long int *, double *, double *, double *,
	double *, double *, double *, double *, double *, double *, double *);

void MTRIDIAGONAL_mp_CLEAN_TRIDIAGONAL();

void TURBULENCE_mp_CLEAN_TURBULENCE();

void GetElementCentralData(double *dest, double *source, double *Jacobian, \
	double *wq, double *Vq, ptrdiff_t RVq, ptrdiff_t CVq, double *LAV);

double * TURBULENCE_mp_TKE;
double * TURBULENCE_mp_EPS;
double * TURBULENCE_mp_NUM;
double * TURBULENCE_mp_NUH;
double * TURBULENCE_mp_L;

#endif
