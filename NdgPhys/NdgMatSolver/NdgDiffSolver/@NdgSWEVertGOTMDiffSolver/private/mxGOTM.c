#include "mxGOTM.h"
#include "blas.h"
#include <math.h>
#include <omp.h>
#include "../../../../../NdgMath/NdgMemory.h"
#include "../../../../../NdgMath/NdgMath.h"

/*the Von kamma constant*/
double kappa = 0.40;

void getGotmDate(int index, int nlev){
	for (int i = 0; i < nlev + 1; i++){
		tkeGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_TKE[i];
		epsGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_EPS[i];
		LGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_L[i];
		nuhGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_NUH[i];
		numGOTM[index*(nlev + 1) + i] = TURBULENCE_mp_NUM[i];
	}
}

void setGotmDate(int index, int nlev){
	for (int i = 0; i < nlev + 1; i++){
		TURBULENCE_mp_TKE[i] = tkeGOTM[i + index*(nlev + 1)];
		TURBULENCE_mp_EPS[i] = epsGOTM[i + index*(nlev + 1)];
		TURBULENCE_mp_L[i] = LGOTM[i + index*(nlev + 1)];
		TURBULENCE_mp_NUH[i] = nuhGOTM[i + index*(nlev + 1)];
		TURBULENCE_mp_NUM[i] = numGOTM[i + index*(nlev + 1)];
	}
}


void InitTurbulenceModelGOTM(long long int *NameList, char * buf, long long int buflen,\
	long long int nlev, int K2d){

	TURBULENCE_mp_INIT_TURBULENCE(NameList, buf, &nlev, buflen);

	MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(&nlev);

	for (int i = 0; i < K2d; i++){
		getGotmDate(i, (int)nlev);
	}
}

void InterpolationToCentralPoint(double *fphys, double *dest, int K2d, int Np3d, \
	int nlayer, double *J, double *wq, double *Vq, ptrdiff_t RVq, ptrdiff_t Cvq, \
	double *LAV) {
	char *transA = "N";
	char *transB = "N";
	double *A = Vq;
	ptrdiff_t ROPA = RVq;
	ptrdiff_t COPA = Cvq;
	ptrdiff_t COPB = 1;
	double Alpha = 1.0;
	double Beta = 0.0;
	ptrdiff_t LDA = RVq;
	ptrdiff_t LDB = RVq;
	ptrdiff_t LDC = LDA;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++) {
		for (int L = 0; L < nlayer; L++) {
//			GetMeshAverageValue(dest + i*nlayer + L, LAV + i*nlayer + L, transA, transB,\
				&ROPA, &COPB, &COPA, &Alpha, A, &LDA, fphys + i*nlayer*Np3d + L*Np3d,\
				J + i*nlayer*Np3d + L*Np3d, &LDB, &Beta, &LDC, wq);
			dest[i*nlayer + L] = 0.0;
			GetMeshAverageValue(dest + i*nlayer + L, LAV + i*nlayer + nlayer-1-L, transA, transB, \
				&ROPA, &COPB, &COPA, &Alpha, A, &LDA, fphys + i*nlayer*Np3d + (nlayer-1-L)*Np3d, \
				J + i*nlayer*Np3d + (nlayer - 1 - L)*Np3d, &LDB, &Beta, &LDC, wq);
		}
	}
}

void InterpolationToCentralPointO(double *fphys, double *dest, int K2d, int Np2d, int Np3d, \
	int nlayer, double *J2d, double *wq2d, double *Vq2d, ptrdiff_t RVq2d, ptrdiff_t Cvq2d,\
    double *LAV) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++) {
		for (int L = 0; L < nlayer; L++) {
			double BottomAve = 0.0;
			double SurfaceAve = 0.0;
			//The average data in dest is arranged from bottom to surface
			GetElementCentralData(&BottomAve, fphys + i*nlayer*Np3d + (nlayer - L - 1)*Np3d, J2d + i*Np2d, wq2d, Vq2d, RVq2d, Cvq2d, LAV + i);
			GetElementCentralData(&SurfaceAve, fphys + i*nlayer*Np3d + (nlayer - L - 1)*Np3d + Np3d - Np2d, J2d + i*Np2d, wq2d, Vq2d, RVq2d, Cvq2d, LAV+i);
			dest[i*nlayer + L] = (BottomAve + SurfaceAve) / 2.0;
		}
	}
}

void mapCentralPointDateToVerticalDate(double *centralDate, double *verticalLineDate, int K2d, \
	int nlev){
	//This has been verified by tests
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		for (int L = 0; L < nlev; L++){
			/*The first data for each vertical line is left undefined*/
			verticalLineDate[k*(nlev + 1) + L + 1] = \
				centralDate[k*nlev + L];
		}
	}
}

void CalculateWaterDepth(int K2d, double hcrit, int nlev){
   /*if the water depth is larger than the threshold value, the layer height is calculated by the ratio of water depth to number of lyers.
   *For the current version, we only consider the equalspace division in vertical. When the water depth is less than the threshold, the water 
   *depth is set to be zero
   */
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++){
		if (hcenter[i] >= hcrit){
			for (int L = 1; L < nlev + 1; L++){ 
				/*The first value of layer height for each vertical line is left undefined*/
				layerHeight[i * (nlev + 1) + L] = hcenter[i] / nlev;
			}
		}
	}

}

void CalculateShearFrequencyDate(int K2d, double hcrit, int nlev){
	//SS = $(\frac{\partial u}{\partial x})^2+(\frac{\partial v}{\partial y})^2$
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++){
		// L stands for lower, U stands for upper, Old stands for the data at previous step(without vertical eddy considered)
		// New stands for the data at the updated step(with vertical eddy considered)
		double uLOld, uUOld, vLOld, vUOld, uLNew, uUNew, vLNew, vUNew;
		if (hcenter[i] >= hcrit){
			for (int L = 1; L < nlev; L++){
				uLOld = huVerticalLine[i*(nlev + 1) + L];
				uUOld = huVerticalLine[i*(nlev + 1) + L + 1];
				uLNew = huVerticalLineNew[i*(nlev + 1) + L];
				uUNew = huVerticalLineNew[i*(nlev + 1) + L + 1];
				vLOld = hvVerticalLine[i*(nlev + 1) + L];
				vUOld = hvVerticalLine[i*(nlev + 1) + L + 1];
				vLNew = hvVerticalLineNew[i*(nlev + 1) + L];
				vUNew = hvVerticalLineNew[i*(nlev + 1) + L + 1];
				/*$SS = \left(\frac{\partial u}{\partial z}\right)^2 + \left(\frac{\partial v}{\partial z}\right)^2 = \\
				\frac{(\hat u_{j+1} - \hat u_j)\times 0.5\times (\hat u_{j+1} + u_{j+1} - \hat u_j - u_j)}{(z_{j+1}-z_j)^2} + 
				\frac{(\hat v_{j+1} - \hat v_j)\times 0.5\times (\hat v_{j+1} + v_{j+1} - \hat v_j - v_j)}{(z_{j+1}-z_j)^2}$,
				here $\hat u_{j+1}$ stands for uUNew, $\hat u_{j}$ uLNew, $u_{j+1}$ uUOld, $u_{j}$ uLOld. Also,  
				$\hat v_{j+1}$ stands for vUNew, $\hat v_{j}$ vLNew, $v_{j+1}$ vUOld, $v_{j}$ vLOld.
				*/
				shearFrequencyDate[i*(nlev + 1) + L] = max( (uUNew - uLNew)*0.5*(uUNew + uUOld - uLNew - uLOld)/ pow(0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]),2.0) + \
					(vUNew - vLNew)*0.5*(vUNew + vUOld - vLNew - vLOld) / pow(0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]), 2.0), 0);
			}
			//For each vertical segment, we have SS(0) = SS(1), SS(nlev) = SS(nlev - 1)
			shearFrequencyDate[i*(nlev + 1)] = shearFrequencyDate[i*(nlev + 1) + 1];
			shearFrequencyDate[i*(nlev + 1) + nlev] = shearFrequencyDate[i*(nlev + 1) + nlev - 1];
		}
		else {
			for (int L = 0; L < nlev + 1; L++) {
				shearFrequencyDate[i*(nlev + 1) + L] = 0;
			}
		}
	}

}

void CalculateLengthScaleAndShearVelocity(double z0b, double z0s, double hcrit, double *DragCoefficient, \
	double *Taux, double *Tauy, int Np2d, int K2d, int nlev) {
	/*for surface friction length, another way is the charnock method*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++) {
		SurfaceFrictionLength[i] = z0s;
		//SurfaceFrictionLength[i] = 0.0030000000260770321;
		SurfaceFrictionVelocity[i] = 0;
		for (int p = 0; p < Np2d; p++) {
			DragCoefficient[i*Np2d + p] = 0;
		}		
		double rr = 0;
		if (hcenter[i] >= hcrit) {
			rr = kappa / log((z0b + layerHeight[i*(nlev + 1) + 1] / 2) / z0b);
			BottomFrictionVelocity[i] = rr * sqrt(pow((huVerticalLine[i*(nlev + 1) + 1] / hcenter[i]), 2) + pow((hvVerticalLine[i*(nlev + 1) + 1] / hcenter[i]), 2));
			BottomFrictionLength[i] = z0b;
			//BottomFrictionLength[i] = 0.0030000000260770321;
			/*Formula by SLIM, we note that Taux and Tauy here have already devided by rho0*/
			SurfaceFrictionVelocity[i] = sqrt(hypot(Taux[i], Tauy[i]));
			//SurfaceFrictionVelocity[i] = 0.031622776601683791;
			/*Formula by SLIM*/
			// DragCoefficient[i] = pow(rr,2);   //Cd = rr^2;
			for (int p = 0; p < Np2d; p++) {
				DragCoefficient[i*Np2d + p] = pow(rr, 2);
			}
			/*Formula by GOTM, see friction.F90(Line 122) and uequation.F90(Line 182)*/
			//DragCoefficient[p] = pow(rr, 2) / layerHeight[p*(nlev + 1) + 1]
		}
	}
}

 void DGDoTurbulence(double *TimeStep, double hcrit, double *Grass, int K2d, long long int nlev){
	 //For the current version, grass is not considered
	 for (int i = 0; i < K2d; i++){
		 if (hcenter[i] >= hcrit){
			 setGotmDate(i, (int)nlev);
			 TURBULENCE_mp_DO_TURBULENCE(&nlev, TimeStep, hcenter + i, SurfaceFrictionVelocity + i, BottomFrictionVelocity + i, SurfaceFrictionLength + i, \
				 BottomFrictionLength + i, layerHeight + i*(nlev + 1), buoyanceFrequencyDate + i*(nlev + 1), shearFrequencyDate + i*(nlev + 1), Grass);
			 getGotmDate(i, (int)nlev);
		 }
	 }

}

 void mapVedgeDateToDof(double *SourceDate, double *DestinationDate, int Np2d, int K2d, int Np3d, int nlev){
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	 for (int k = 0; k < K2d; k++){
		 for (int L = 0; L < nlev; L++) {
			 for (int p = 0; p < Np2d; p++) {
				 DestinationDate[k*nlev*Np3d + (nlev - L - 1)*Np3d + p] = SourceDate[k*(nlev + 1) + L];//The bottom layer of the up cell
				 DestinationDate[k*nlev*Np3d + (nlev - L - 1)*Np3d + Np3d - Np2d + p] = SourceDate[k*(nlev + 1) + L + 1];  //The top layer of the down cell
			 }
		 }
	 }
 }

 void CalculateBuoyanceFrequencyDate( int Np2d, int K2d, double hcrit, int nlev, \
	 double gra, double rho0){

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	 for (int i = 0; i < K2d; i++) {
		 if (hcenter[i] >= hcrit) {
			 for (int L = 1; L < nlev; L++) {
				 buoyanceFrequencyDate[i*(nlev + 1) + L] = -1 * gra/rho0*(rhoVerticalLine[i*(nlev + 1) + L + 1] - rhoVerticalLine[i*(nlev + 1) + L]) / (0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]));
			 }
			 //For each vertical segment, we have NN(0) = NN(1), NN(nlev) = NN(nlev - 1)
			 buoyanceFrequencyDate[i*(nlev + 1)] = buoyanceFrequencyDate[i*(nlev + 1) + 1];
			 buoyanceFrequencyDate[i*(nlev + 1) + nlev] = buoyanceFrequencyDate[i*(nlev + 1) + nlev - 1];
		 }
		 else
		 {
			 for (int L = 0; L < nlev + 1; L++) {
				 buoyanceFrequencyDate[i*(nlev + 1) + L] = 0;
			 }
		 }
	 }
 }

 void GetElementCentralData(double *dest, double *source, double *Jacobian, \
	 double *wq, double *Vq, ptrdiff_t RVq, ptrdiff_t CVq, double *LAV) {
	 char *transA = "N";
	 char *transB = "N";
	 ptrdiff_t COPB = 1;
	 double Alpha = 1.0;
	 double Beta = 0.0;
	 ptrdiff_t LDA = RVq;
	 ptrdiff_t LDB = RVq;
	 ptrdiff_t LDC = LDA;
	 GetMeshAverageValue(dest, LAV, transA, transB, &RVq, &COPB, &CVq, &Alpha, Vq, \
		 &LDA, source, Jacobian, &LDB, &Beta, &LDC, wq);
 }
