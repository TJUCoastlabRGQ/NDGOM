#include "..\..\..\..\..\NdgMath\NdgMath.h"

typedef enum {
	One = 1,
	Two = 2,
	Three = 3
} NdgMeshType;


void GetFacialFluxTerm(double *, double *, double *, int);

void GetCentralNumFlux(double *, double *, double *, double *, int);

void GetIPNumFlux(double *, double *, double *, double *, double *, double *, double *, int , double *, double );

void GetIPBoundaryNumFlux(double *, double *, double *, double *, double *, double *, int, double *, double);

void GetIEContributionToAuxialaryVariable(double *, int, int, int, int, double *, double *, double *, \
	double *, double *, int, double *, double *, double *, double *, double *, double *);

void GetBEContributionToAuxialaryVariable(double *, int , int , int , int , double *, double *, double *, \
	double *, int , double *, double *, double *, double *, double *);


void GetIEContributionToRHS(double *, double *, double *, int , double *, double *, double *, \
	double *, int , int , int , int , double *, double *, double *, double *, double *, double *, double , \
	double *, double *, double *, double *, double *, double *, double *);


void GetBEContributionToRHS(double *, double *, int , double *, double *, double *, \
	int , int , int , int , double *, double *, double *, double *, double *, double *, double , \
	double *, double *, double *, double *, double *);