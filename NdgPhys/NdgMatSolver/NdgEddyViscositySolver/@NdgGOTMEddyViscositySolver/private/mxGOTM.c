#include "mxGOTM.h"

Void InitTurbulenceModelGOTM(long long int *NameList, char * buf, long long int *NumOfVertLayer, long long int *buflen){
	TURBULENCE_mp_INIT_TURBULENCE(Namelist, buf, NumOfVertLayer, *buflen);
	MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(NumOfVertLayer);
}