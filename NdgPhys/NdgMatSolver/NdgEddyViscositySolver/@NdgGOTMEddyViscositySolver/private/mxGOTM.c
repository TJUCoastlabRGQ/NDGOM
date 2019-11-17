#include "mxGOTM.h"
#include "blas.h"
#include <math.h>

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
				verticalLineData[k*(nlev + 1)*Np2d + p*(nlev + 1) + L] = \
					centralDate[k*nlev*Np2d + (nlev - L)*Np2d + p];
			}
		}
	}
}

void CalculateWaterDepth(double *H2d, double *TempLayerHeight){
   /*if the water depth is larger than the threshold value, the layer height is calculated by the ratio of water depth to number of lyers.
   *For the current version, we only consider the equalspace division in vertical. When the water depth is less than the threshold, the water 
   *depth is set to be zero
   */
	for (int p = 0; p < Np2d * K2d; p++){
		if (H2d[p] >= hcrit){
			for (int L = 1; L < nlev + 1; L++){
				TempLayerHeight[p * (nlev + 1) + L] = H2d[p] / nelv;
			}
		}
	}

}


void CalculateShearFrequencyDate(double *H2d, double *LayerHeight, double *huVertial, double *hvVertical, double *shearProductionDest){
	//SS = $(\frac{\partial u}{\partial x})^2+(\frac{\partial v}{\partial y})^2$
	for (int p = 0; p < Np2d*K2d; p++){
		if (H2d(p) >= hcrit){
			for (int L = 1; L < nlev; L++){
				shearProductionDest[p*(nlev + 1) + L] = ((huVertial[p*(nlev + 1) + L + 1] - huVertial[p*(nlev + 1) + L]) / H2d[p] / (0.5*(LayerHeight[L + 1] + LayerHeight[L]))) ^ 2 + \
					((hvVertial[p*(nlev + 1) + L + 1] - hvVertial[p*(nlev + 1) + L)] / H2d[p] / (0.5*(LayerHeight[L + 1] + LayerHeight[L]))) ^ 2;
			}
			//For each vertical segment, we have SS(0) = SS(1), SS(nlev) = SS(nlev - 1)
			shearProductionDest[p*(nlev + 1)] = shearProductionDest[p*(nlev + 1) + 1];
			shearProductionDest[p*(nlev + 1) + nlev] = shearProductionDest[p*(nlev + 1) + nlev - 1];
		}
	}

}

void CalculateLengthScaleAndShearVelocity(double *H2d, double* LayerHeight, double *huvl, double *hvvl, double* z0b, double *utaub, double *z0s, double *utaus){
	/*for surface friction length, another way is the charnock method*/
	for (int p = 0; p < Np2d * K2d; p++){
		z0s[p] = z0s_min;
		utaus[p] = 0;
		double Tempz0b = 0;
		if (H2d[p] >= hcrit){
			for (int itera = 0; itera < NMaxItration; itera++){
				Tempz0b = 0.1*avmolu / max(avmolu, utaub[p]) + 0.03*h0b;
				double rr = kappa / log((Tempz0b + LayerHeight[p*(nlev + 1) + 1] / 2) / Tempz0b) + 0.03*h0b;
				utaub[p] = rr * sqrt((huvl[p*(nlev + 1) + 1] / H2d[p]) ^ 2 + (hvvl[p*(nlev + 1) + 1] / H2d[p]) ^ 2);
			}
			z0b[p] = Tempz0b;
		}
	}
}

 void DGDoTurbulence(double *TimeStep, double *H2d, double *utaus, double *utaub, double *z0s, double *z0b, double *LayerHeight, double *NN, double *SS, double *Grass, double *EddyViscosity){
	 //For the current version, grass is not considered
	 for (int p = 0; p < Np2d * K2d; p++){
		 if (H2d[p] >= hcrit){
			 TURBULENCE_mp_DO_TURBULENCE(&nlev, TimeStep, H2d + p, utaus + p, utaub + p, z0s + p, z0b + p, LayerHeight + p*( nlev + 1),\
				 NN + p*(nlev + 1), SS + p*(nlev + 1), NULL);
			 for (int L = 0; L < nlev + 1; L++){
				 EddyViscosity[p*(nlev + 1) + L] = TURBULENCE_mp_NUM[L];
			 }
		 }
	 }

}

 void mapVedgeDateToDof(double *SourceDate, double *DestinationDate){


 }