#include "mxGOTM.h"
#include "blas.h"

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

void InterpolationToCentralPoint(double *VCV, double *fphys, double *dest){
	char *chn = "N";
	double alpha = 1;
	double beta = 0;
	dgemm(chn, chn, &Np2d, &K3d, &K3d, &alpha, VCV, &Np2d, fphys, &Np3d, &beta, dest, &Np2d);
}

void mapCentralPointDateToVerticalDate(double *centralDate, double *verticalLineDate){
	//This has been verified by tests
	for (int k = 0; k < K2d; k++){
		for (int L = 1; L < nlev + 1; L++){
			for (int p = 0; p < Np2d; p++){
				verticalLineData(k*(nlev + 1)*Np2d + p*(nlev + 1) + L) = \
					centralDate(k*nlev*Np2d + (nlev - L)*Np2d + p);
			}
		}
	}
}

void CalculateShearProductionDate(double *H2d, double *huVertial, double *hvVertical, double *shearProductionDate){

}