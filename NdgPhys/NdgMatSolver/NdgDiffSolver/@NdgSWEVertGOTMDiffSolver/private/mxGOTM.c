#include "mxGOTM.h"
#include "blas.h"
#include <math.h>
#include <omp.h>
#include "../../../../../NdgMath/NdgMemory.h"

/*the Von kamma constant*/
double kappa = 0.41;

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


void InitTurbulenceModelGOTM(long long int *NameList, char * buf, long long int buflen,\
	long long int nlev, int Np2d, int K2d){

	TURBULENCE_mp_INIT_TURBULENCE(NameList, buf, &nlev, buflen);

	MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(&nlev);

	for (int i = 0; i < Np2d*K2d; i++){
		getGotmDate(i, nlev);
	}
}

void InterpolationToCentralPoint(double *fphys, double *dest, ptrdiff_t *Np2d, ptrdiff_t *K3d,\
	ptrdiff_t *Np3d, double *VCV){
	char *chn = "N";
	double alpha = 1;
	double beta = 0;
	ptrdiff_t Col = 1;
    //ptrdiff_t TempNp2d = (ptrdiff_t)Np2d, TempK3d = (ptrdiff_t)K3d, TempNp3d = (ptrdiff_t)Np3d;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < (int)(*K3d); i++){
		dgemm(chn, chn, Np2d, &Col, Np3d, &alpha, VCV, Np2d, fphys + i*(int)(*Np3d), Np3d, &beta, dest + i*(int)(*Np2d), Np2d);
	}
}

void mapCentralPointDateToVerticalDate(double *centralDate, double *verticalLineDate, int K2d, \
	long long int nlev, int Np2d){
	//This has been verified by tests
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
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
#pragma omp parallel for num_threads(DG_THREADS)
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
#pragma omp parallel for num_threads(DG_THREADS)
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

void CalculateLengthScaleAndShearVelocity(double z0b, double z0s, double *H2d, double hcrit, double *DragCoefficient, \
	double *Taux, double *Tauy, int Np2d, int K2d, long long int nlev) {
	/*for surface friction length, another way is the charnock method*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int p = 0; p < Np2d * K2d; p++) {
		SurfaceFrictionLength[p] = z0s;
		SurfaceFrictionVelocity[p] = 0;
		DragCoefficient[p] = 0;
		double rr = 0;
		if (H2d[p] >= hcrit) {
			rr = kappa / log((z0b + layerHeight[p*(nlev + 1) + 1] / 2) / z0b);
			BottomFrictionVelocity[p] = rr * sqrt(pow((huVerticalLine[p*(nlev + 1) + 1] / H2d[p]), 2) + pow((hvVerticalLine[p*(nlev + 1) + 1] / H2d[p]), 2));
			BottomFrictionLength[p] = z0b;
			/*Note that, the surface-stress here must be devided by rho_0, the following calculation is according to friction.F90 in gotm*/
			/*Formula by GOTM*/
			//SurfaceFrictionVelocity[p] = pow(pow(Taux[p] / rho0, 2) + pow(Tauy[p] / rho0, 2), 1 / 4);
			/*Formula by SLIM, we note that Taux and Tauy here have already devided by rho0*/
			SurfaceFrictionVelocity[p] = sqrt(hypot(Taux[p], Tauy[p]));
			/*Formula by SLIM*/
			//SurfaceFrictionVelocity[p] = pow(pow(pow(Taux[p], 2) + pow(Tauy[p], 2),1/2)/rho0, 1 / 2);
			/*Formula by SLIM*/
			DragCoefficient[p] = pow(rr,2);   //Cd = rr^2;
			/*Formula by GOTM, see friction.F90(Line 122) and uequation.F90(Line 182)*/
			//DragCoefficient[p] = pow(rr, 2) / layerHeight[p*(nlev + 1) + 1]
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
			//	 eddyDiffusionDate[p*(nlev + 1) + L] = TURBULENCE_mp_NUH[L];
				 eddyTKEDate[p*(nlev + 1) + L] = TURBULENCE_mp_TKE[L];
				 eddyLengthDate[p*(nlev + 1) + L] = TURBULENCE_mp_L[L];
				 eddyEPSDate[p*(nlev + 1) + L] = TURBULENCE_mp_EPS[L];
			 }
		 }
	 }

}

 void mapVedgeDateToDof(double *SourceDate, double *DestinationDate, int Np2d, int K2d, int Np3d, long long int nlev){
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
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

 void CalculateBuoyanceFrequencyDate(double *H2d, int Np2d, int K2d, double hcrit, long long int nlev, \
	 double gra, double rho0){

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	 for (int p = 0; p < Np2d*K2d; p++) {
		 if (H2d[p] >= hcrit) {
			 for (int L = 1; L < nlev; L++) {
				 shearFrequencyDate[p*(nlev + 1) + L] = -1 * gra/rho0*(rhoVerticalLine[p*(nlev + 1) + L + 1] - rhoVerticalLine[p*(nlev + 1) + L]) / (0.5*(layerHeight[p*(nlev + 1) + L + 1] + layerHeight[p*(nlev + 1) + L]));
				 shearFrequencyDate[p*(nlev + 1) + L] = max(shearFrequencyDate[p*(nlev + 1) + L], 0.0);
			 }
			 //For each vertical segment, we have NN(0) = NN(1), NN(nlev) = NN(nlev - 1)
			 buoyanceFrequencyDate[p*(nlev + 1)] = buoyanceFrequencyDate[p*(nlev + 1) + 1];
			 buoyanceFrequencyDate[p*(nlev + 1) + nlev] = buoyanceFrequencyDate[p*(nlev + 1) + nlev - 1];
		 }
		 else
		 {
			 for (int L = 0; L < nlev; L++) {
				 shearFrequencyDate[p*(nlev + 1) + L] = 0;
			 }
		 }
	 }
 }
