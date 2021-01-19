#include "mxGOTM.h"
#include "blas.h"
#include <math.h>
#include <omp.h>
#include "../../../../../NdgMath/NdgMemory.h"

/*the Von kamma constant*/
double kappa = 0.41;
/*the limit value for the surface roughness*/
double z0s_min = 0.02;
//double z0s_min = 0;
/*this value is used for computing bottom roughness*/
double avmolu = 1.3e-6;

void getGotmDate(int index, long long int nlev){
	for (int i = 0; i < nlev + 1; i++){
		tkeGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_TKE[i];
		epsGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_EPS[i];
		LGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_L[i];
		nuhGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_NUH[i];
		numGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_NUM[i];
	}
}

void setGotmDate(int index, long long int nlev){
	for (int i = 0; i < nlev + 1; i++){
		TURBULENCE_mp_TKE[i] = tkeGOTM[i + index*(nlev + 1)];
		TURBULENCE_mp_EPS[i] = epsGOTM[i + index*(nlev + 1)];
		TURBULENCE_mp_L[i] = LGOTM[i + index*(nlev + 1)];
		TURBULENCE_mp_NUH[i] = nuhGOTM[i + index*(nlev + 1)];
		TURBULENCE_mp_NUM[i] = numGOTM[i + index*(nlev + 1)];
	}
}


void InitTurbulenceModelGOTM(long long int *NameList, char * buf, long long int buflen, long long int nlev, int Np2d, int K2d){

	TURBULENCE_mp_INIT_TURBULENCE(NameList, buf, &nlev, buflen);

	MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(&nlev);

	for (int i = 0; i < Np2d*K2d; i++){
		getGotmDate(i, nlev);
	}
}

void InterpolationToCentralPoint(double *fphys, double *dest, ptrdiff_t *Np2d, ptrdiff_t *K3d, ptrdiff_t *Np3d, double *VCV){
	char *chn = "N";
	double alpha = 1;
	double beta = 0;
	ptrdiff_t Col = 1;
    //ptrdiff_t TempNp2d = (ptrdiff_t)Np2d, TempK3d = (ptrdiff_t)K3d, TempNp3d = (ptrdiff_t)Np3d;
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int i = 0; i < (int)(*K3d); i++){
		dgemm(chn, chn, Np2d, &Col, Np3d, &alpha, VCV, Np2d, fphys + i*(int)(*Np3d), Np3d, &beta, dest + i*(int)(*Np2d), Np2d);
	}
}

void mapCentralPointDateToVerticalDate(double *centralDate, double *verticalLineDate, int K2d, long long int nlev, int Np2d){
	//This has been verified by tests
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K2d; k++){
		for (int L = 1; L < nlev + 1; L++){
			for (int p = 0; p < Np2d; p++){
				verticalLineDate[k*(nlev + 1)*Np2d + p*(nlev + 1) + L] = \
					centralDate[k*nlev*Np2d + (nlev - L)*Np2d + p];
			}
		}
	}
}

void CalculateWaterDepth(double *H2d, int Np2d, int K2d, double hcrit, long long int nlev){
   /*if the water depth is larger than the threshold value, the layer height is calculated by the ratio of water depth to number of lyers.
   *For the current version, we only consider the equalspace division in vertical. When the water depth is less than the threshold, the water 
   *depth is set to be zero
   */
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int p = 0; p < Np2d * K2d; p++){
		if (H2d[p] >= hcrit){
			for (int L = 1; L < nlev + 1; L++){
				layerHeight[p * (nlev + 1) + L] = H2d[p] / nlev;
			}
		}
	}

}


void CalculateShearFrequencyDate(double *H2d, int Np2d, int K2d, double hcrit, long long int nlev){
	//SS = $(\frac{\partial u}{\partial x})^2+(\frac{\partial v}{\partial y})^2$
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int p = 0; p < Np2d*K2d; p++){
		if (H2d[p] >= hcrit){
			for (int L = 1; L < nlev; L++){
				shearFrequencyDate[p*(nlev + 1) + L] = pow((huVerticalLine[p*(nlev + 1) + L + 1] - huVerticalLine[p*(nlev + 1) + L]) / H2d[p] / (0.5*(layerHeight[p*(nlev + 1) + L + 1] + layerHeight[p*(nlev + 1) + L])), 2) \
					+ pow((hvVerticalLine[p*(nlev + 1) + L + 1] - hvVerticalLine[p*(nlev + 1) + L]) / H2d[p] / (0.5*(layerHeight[p*(nlev + 1) + L + 1] + layerHeight[p*(nlev + 1) + L])), 2);
			}
			//For each vertical segment, we have SS(0) = SS(1), SS(nlev) = SS(nlev - 1)
			shearFrequencyDate[p*(nlev + 1)] = shearFrequencyDate[p*(nlev + 1) + 1];
			shearFrequencyDate[p*(nlev + 1) + nlev] = shearFrequencyDate[p*(nlev + 1) + nlev - 1];
		}
	}

}

void CalculateLengthScaleAndShearVelocity(double *H2d, double hcrit, double *DragCoefficient, double *Taux, double *Tauy, int Np2d, int K2d, int NMaxItration, double h0b, long long int nlev){
	/*for surface friction length, another way is the charnock method*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int p = 0; p < Np2d * K2d; p++){
		SurfaceFrictionLength[p] = z0s_min;
		SurfaceFrictionVelocity[p] = 0;
        DragCoefficient[p] = 0;
		double Tempz0b = 0;
        double rr = 0;
		if (H2d[p] >= hcrit){
			for (int itera = 0; itera < NMaxItration; itera++){
				Tempz0b = 0.1*avmolu / max(avmolu, BottomFrictionVelocity[p]) + 0.03*h0b;
				Tempz0b = 0.0015;
				/*$rr=\frac{\kappa}{log(\frac{z0b+\frac{h(1)}{2}}{z0b})$, where z0b is the bottom roughness length and h(1) is the height of the first layer of cell,
				*detail of this part can refer to the friction.F90 file in GOTM
				*/
				rr = kappa / log((Tempz0b + layerHeight[p*(nlev + 1) + 1] / 2) / Tempz0b);
				BottomFrictionVelocity[p] = rr * sqrt(pow((huVerticalLine[p*(nlev + 1) + 1] / H2d[p]), 2) + pow((hvVerticalLine[p*(nlev + 1) + 1] / H2d[p]), 2));
			}
			BottomFrictionLength[p] = Tempz0b;
            /*Note that, the surface-stress here must be devided by rho_0*/
            SurfaceFrictionVelocity[p] = pow(pow(Taux[p],2) + pow(Tauy[p],2),1/4);
            DragCoefficient[p] = pow(rr,2);
            /*According to GOTM, DragCoefficient[p] = pow(rr,2)/layerHeight[p*(nlev + 1) + 1] should be adopted, see uequation.F90*/
		}
	}
}

 void DGDoTurbulence(double *TimeStep, double *H2d, double hcrit, double *Grass, int Np2d, int K2d, long long int nlev){
	 //For the current version, grass is not considered
	 for (int p = 0; p < Np2d * K2d; p++){
		 if (H2d[p] >= hcrit){
			 setGotmDate(p, nlev);
			 TURBULENCE_mp_DO_TURBULENCE(&nlev, TimeStep, H2d + p, SurfaceFrictionVelocity + p, BottomFrictionVelocity + p, SurfaceFrictionLength + p, \
				 BottomFrictionLength + p, layerHeight + p*(nlev + 1), buoyanceFrequencyDate + p*(nlev + 1), shearFrequencyDate + p*(nlev + 1), Grass);
			 getGotmDate(p, nlev);
			 for (int L = 0; L < nlev + 1; L++){
				 eddyViscosityDate[p*(nlev + 1) + L] = TURBULENCE_mp_NUM[L];
			 }
		 }
	 }

}

 void mapVedgeDateToDof(double *SourceDate, double *DestinationDate, int Np2d, int K2d, int Np3d, long long int nlev){
#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	 for (int k = 0; k < K2d; k++){
		 for (int p = 0; p < Np2d; p++){
			 DestinationDate[k*nlev*Np3d + (nlev - 1)*Np3d + p] = SourceDate[k*Np2d*(nlev + 1) + p*(nlev+1)];//the down face of the bottommost cell for each column
			 DestinationDate[k*nlev*Np3d + p + Np2d] = SourceDate[k*Np2d*(nlev + 1) + p*(nlev+1) + nlev];//the upper face of the topmost cell for each column
			 for (int L = 1; L < nlev; L++){
				 DestinationDate[k*nlev*Np3d + (nlev - L)*Np3d + p + Np2d] = SourceDate[k*Np2d*(nlev + 1) + p*(nlev+1) + L];  //The top layer of the down cell
				 DestinationDate[k*nlev*Np3d + (nlev - L-1)*Np3d + p] = SourceDate[k*Np2d*(nlev + 1) + p*(nlev+1) + L];//The bottom layer of the up cell
			 }
		 }
	 }
 }

 void CalculateBuoyanceFrequencyDate(int Np2d, int K2d, long long int nlev){

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	 for (int i = 0; i < Np2d*K2d*(nlev + 1); i++)
		 buoyanceFrequencyDate[i] = 0;
 }