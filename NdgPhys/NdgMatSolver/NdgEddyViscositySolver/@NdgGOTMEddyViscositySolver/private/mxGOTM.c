#include "mxGOTM.h"

/*the Von kamma constant*/
double kappa = 0.4;
/*the limit value for the bottom roughness*/
double z0s_min = 0.02;
/*this value is used for computing bottom roughness*/
double avmolu = 1.3e-6;
/*The iteration times to calculate the friction velocity, this is set to be 2 by default*/
int NMaxItration = 2;


void InitTurbulenceModelGOTM(long long int *NameList, char * buf, long long int *NumOfVertLayer, long long int buflen){
	TURBULENCE_mp_INIT_TURBULENCE(NameList, buf, NumOfVertLayer, buflen);
	MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(NumOfVertLayer);
}