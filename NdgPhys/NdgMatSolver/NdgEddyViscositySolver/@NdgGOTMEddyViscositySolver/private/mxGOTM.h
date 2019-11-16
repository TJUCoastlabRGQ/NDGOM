#ifndef __mxGOTM_H__
#define __mxGOTM_H__
#include "mex.h"
/*number of elements and interpolation points for the two-dimensional mesh*/
int K2d, Np2d;
/*number of elements and interpolation points for the three-dimensional mesh*/
int K3d, Np3d;
/*layers of the three dimensional shallow water model*/
long long int nlev;
/*the critical water depth, when the depth is less than this value, velosities are set to be zero*/
double hcrit;
/*the physical roughness length scale*/
double h0b;
/*the Von kamma constant*/
double kappa = 0.4;
/*the limit value for the bottom roughness*/
double z0s_min = 0.02;
/*this value is used for computing bottom roughness*/
double avmolu = 1.3e-6;
/*final time of the */
double finalTime;
/*The iteration times to calculate the friction velocity, this is set to be 2 by default*/
int MaxItration = 2;
/*Whether GOTM is initialized or not, this variable is false by default*/
char *Initialized = 'False';
/*This function is used to interpolate the physical value from interpolation point to central point in vertical direction*/
void InterpolateToCentralPoint(double *InterpolationMatrix, double *FphysField, double *destination);
/*This function is used to map the date located at the central point to the vertical line that the GOTM adapted*/
void MapCentralPointDataToVerticalLine(double *CentralPointData, double *VerticalLineTempdata);
/*This function is used to calculate the shear production term*/
void CalculateShearProductioinTerm(double *H2d, double *UVerticalLineTempdata, double *VVerticalLineTempdata, double *ShearProductionData);
// CalculateBuoyanceProductionTerm to be added 
//Here, z0b is the bottom roughness, utaub is the friction velocity, z0s is the surface roughness
void CalculateLengthScaleAndShearVelocity(double *H2d, double* z0b, double *utaub, double *z0s, double *utaus); 
/*This is the interface to do_turbulence function in GOTM*/
void DGDoTurbulence(double *H2d, double *ShearProductionDate, double *buoyanceProductionDate, double *utaus, double *utaub, double *z0s, double *z0b);
/*This function is used to map the date calculated by GOTM to the output matrix*/
void mapVedgeDateToDof(double *VedgeDate, double *exportDofContainer);
/*This function is used to initialize the GOTM module*/
void InitTurbulenceModelGOTM(long long int *, char * , long long int *, long long int *)

/*The following is the GOTM part*/

void __stdcall TURBULENCE_mp_INIT_TURBULENCE(long long int *, char *, long long int *, long long int);

//The dealocation part need to be considered

void __stdcall MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(long long int *);

double * TURBULENCE_mp_NUM;

#endif