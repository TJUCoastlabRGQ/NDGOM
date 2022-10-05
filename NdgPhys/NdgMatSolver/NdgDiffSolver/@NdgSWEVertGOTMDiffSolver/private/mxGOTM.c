#include "mxGOTM.h"
#include "blas.h"
#include <math.h>
#include <omp.h>
#include "../../../../../NdgMath/NdgMemory.h"
#include "../../../../../NdgMath/NdgMath.h"
#include <stdio.h>
#include "../../../../../Application/SWE/SWE3d/@SWEBaroclinic3d/private/eqstate.h"

/*the Von kamma constant*/
double kappa = 0.41587;

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
	long long int nlev, int K2d, int Np2d){

	TURBULENCE_mp_INIT_TURBULENCE(NameList, buf, &nlev, buflen);

	MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(&nlev);

	for (int i = 0; i < Np2d*K2d; i++){
		getGotmDate(i, (int)nlev);
	}
}

void InterpolationToCentralPointO(double *fphys, double *dest, int K2d, int Np3d, \
	int nlayer, double *J, double *wq, double *Vq, ptrdiff_t RVq, ptrdiff_t Cvq, \
	double *LAV3d) {
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
//			GetMeshAverageValue(dest + i*nlayer + L, LAV3d + i*nlayer + L, transA, transB,\
				&ROPA, &COPB, &COPA, &Alpha, A, &LDA, fphys + i*nlayer*Np3d + L*Np3d,\
				J + i*nlayer*Np3d + L*Np3d, &LDB, &Beta, &LDC, wq);
			dest[i*nlayer + L] = 0.0;
			GetMeshAverageValue(dest + i*nlayer + L, LAV3d + i*nlayer + nlayer-1-L, transA, transB, \
				&ROPA, &COPB, &COPA, &Alpha, A, &LDA, fphys + i*nlayer*Np3d + (nlayer-1-L)*Np3d, \
				J + i*nlayer*Np3d + (nlayer - 1 - L)*Np3d, &LDB, &Beta, &LDC, wq);
		}
	}
}

void InterpolationToCentralPoint(double *fphys, double *dest, ptrdiff_t *Np2d, ptrdiff_t *K3d, ptrdiff_t *Np3d, double *VCV) {

	char *chn = "N";
	double alpha = 1;
	double beta = 0;
	ptrdiff_t Col = 1;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < (int)(*K3d); i++) {
		dgemm(chn, chn, Np2d, &Col, Np3d, &alpha, VCV, Np2d, fphys + i*(int)(*Np3d), Np3d, &beta, dest + i*(int)(*Np2d), Np2d);
	}
	/*
	for (int i = 0; i < K2d; i++) {
		for (int L = 0; L < nlayer; L++) {
			double BottomAve = 0.0;
			double SurfaceAve = 0.0;
			//The average data in dest is arranged from bottom to surface
			GetElementCentralData(&BottomAve, fphys + i*nlayer*Np3d + (nlayer - L - 1)*Np3d, J2d + i*Np2d, wq2d, Vq2d, RVq2d, Cvq2d, LAV2d + i);
			GetElementCentralData(&SurfaceAve, fphys + i*nlayer*Np3d + (nlayer - L)*Np3d - Np2d, J2d + i*Np2d, wq2d, Vq2d, RVq2d, Cvq2d, LAV2d +i);
			dest[i*nlayer + L] = (BottomAve + SurfaceAve) / 2.0;
		}
	}
	*/
}


void mapCentralPointDateToVerticalDate(double *centralDate, double *verticalLineDate, int K2d, \
	int nlev, int Np2d){
	//This has been verified by tests
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		for (int L = 1; L < nlev + 1; L++) {
			for (int p = 0; p < Np2d; p++) {
				verticalLineDate[k*(nlev + 1)*Np2d + p*(nlev + 1) + L] = \
					centralDate[k*nlev*Np2d + (nlev - L)*Np2d + p];
			}
		}
	}
	/*
	for (int k = 0; k < K2d; k++){
		for (int L = 0; L < nlev; L++){
			//The first data for each vertical line is left undefined
			verticalLineDate[k*(nlev + 1) + L + 1] = \
				centralDate[k*nlev + L];
		}
	}
	*/
}

void CalculateWaterDepth(double *H2d, int Np2d, int K2d, double hcrit, long long int nlev) {
   /*if the water depth is larger than the threshold value, the layer height is calculated by the ratio of water depth to number of lyers.
   *For the current version, we only consider the equalspace division in vertical. When the water depth is less than the threshold, the water 
   *depth is set to be zero
   */
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int p = 0; p < Np2d * K2d; p++) {
		if (H2d[p] >= hcrit) {
			for (int L = 1; L < nlev + 1; L++) {
				layerHeight[p * (nlev + 1) + L] = H2d[p] / nlev;
			}
		}
	}
}

void CalculateShearFrequencyDate(double *H2d, int Np2d, int K2d, double hcrit, int nlev){
	//SS = $(\frac{\partial u}{\partial x})^2+(\frac{\partial v}{\partial y})^2$
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d*Np2d; i++){
		// L stands for lower, U stands for upper, Old stands for the data at previous step(without vertical eddy considered)
		// New stands for the data at the updated step(with vertical eddy considered)
		double uLOld, uUOld, vLOld, vUOld, uLNew, uUNew, vLNew, vUNew;
		if (H2d[i] >= hcrit){
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
				//shearFrequencyDate[i*(nlev + 1) + L] = max( (uUNew - uLNew)*0.5*(uUNew + uUOld - uLNew - uLOld)/ pow(0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]),2.0) + \
					(vUNew - vLNew)*0.5*(vUNew + vUOld - vLNew - vLOld) / pow(0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]), 2.0), 0);
				shearFrequencyDate[i*(nlev + 1) + L] = max( (uUNew - uLNew)*0.5*(uUNew + uUOld - uLNew - uLOld) / pow(0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]), 2.0) + \
					(vUNew - vLNew)*0.5*(vUNew + vUOld - vLNew - vLOld) / pow(0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]), 2.0), 0.0 );
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

void CalculateLengthScaleAndShearVelocity(double *h2d, double z0b, double z0s, double hcrit, double *DragCoefficient, \
	double *Taux, double *Tauy, int Np2d, int K2d, int nlev) {
	/*for surface friction length, another way is the charnock method*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d*Np2d; i++) {	
		double rr = 0;
		if (h2d[i] >= hcrit) {
			rr = kappa / log((z0b + layerHeight[i*(nlev + 1) + 1] / 2) / z0b);
			BottomFrictionVelocity[i] = rr * sqrt(pow((huVerticalLine[i*(nlev + 1) + 1] / h2d[i]), 2.0) + pow((hvVerticalLine[i*(nlev + 1) + 1] / h2d[i]), 2.0));
			BottomFrictionLength[i] = z0b;
			/*We note that Taux and Tauy here have already devided by rho0*/
			SurfaceFrictionVelocity[i] = sqrt(hypot(Taux[i], Tauy[i]));
			SurfaceFrictionLength[i] = z0s;
			DragCoefficient[i] = pow(rr, 2.0);
			/*
			for (int p = 0; p < Np2d; p++) {
				DragCoefficient[i*Np2d + p] = pow(rr, 2.0);
			}
			*/
		}
		else {
			BottomFrictionVelocity[i] = 0;
			BottomFrictionLength[i] = 0;
			SurfaceFrictionVelocity[i] = 0;
			SurfaceFrictionLength[i] = 0;
			for (int p = 0; p < Np2d; p++) {
				DragCoefficient[i*Np2d + p] = 0;
			}

		}
	}
}

 void DGDoTurbulence(double *h2d, double *TimeStep, double hcrit, double *Grass, int Np2d, int K2d, long long int nlev){
	 //For the current version, grass is not considered
	 for (int i = 0; i < K2d*Np2d; i++){
		 if (h2d[i] >= hcrit){
			 setGotmDate(i, (int)nlev);
			 TURBULENCE_mp_DO_TURBULENCE(&nlev, TimeStep, h2d + i, SurfaceFrictionVelocity + i, BottomFrictionVelocity + i, SurfaceFrictionLength + i, \
				 BottomFrictionLength + i, layerHeight + i*(nlev + 1), buoyanceFrequencyDate + i*(nlev + 1), shearFrequencyDate + i*(nlev + 1), Grass);
			 getGotmDate(i, (int)nlev);
		 }
	 }
}

 void mapVedgeDateToDof(double *SourceDate, double *DestinationDate, int Np2d, int K2d, int Np3d, long long int nlev) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	 for (int k = 0; k < K2d; k++) {
		 for (int p = 0; p < Np2d; p++) {
			 DestinationDate[k*nlev*Np3d + (nlev - 1)*Np3d + p] = SourceDate[k*Np2d*(nlev + 1) + p*(nlev + 1)];//the down face of the bottommost cell for each column
			 DestinationDate[k*nlev*Np3d + p + Np2d] = SourceDate[k*Np2d*(nlev + 1) + p*(nlev + 1) + nlev];//the upper face of the topmost cell for each column
			 for (int L = 1; L < nlev; L++) {
				 DestinationDate[k*nlev*Np3d + (nlev - L)*Np3d + p + Np2d] = SourceDate[k*Np2d*(nlev + 1) + p*(nlev + 1) + L];  //The top layer of the down cell
				 DestinationDate[k*nlev*Np3d + (nlev - L - 1)*Np3d + p] = SourceDate[k*Np2d*(nlev + 1) + p*(nlev + 1) + L];//The bottom layer of the up cell
			 }
		 }
	 }
 }

 void mapDofDateToVedge(double *SourceDate, double *DestinationDate, int K2d, int Np2d, int nlev) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	 for (int ele = 0; ele < K2d; ele++) {
		 for (int p = 0; p < Np2d; p++) {
			 //The first point on the bottom
			 DestinationDate[ele*Np2d*(nlev+1) + p*(nlev + 1)] = SourceDate[ele*Np2d*2*nlev + (nlev-1)*2*Np2d+p];
			 for (int L = 1; L < nlev; L++) {
				 DestinationDate[ele*Np2d*(nlev + 1) + p*(nlev + 1) + L] = 0.5 * SourceDate[ele*Np2d * 2 * nlev + (nlev - L) * 2 * Np2d + Np2d + p] + \
					 0.5*SourceDate[ele*Np2d * 2 * nlev + (nlev - L - 1) * 2 * Np2d + p];
			 }
			 //The point on the surface
			 DestinationDate[ele*Np2d*(nlev + 1) + p*(nlev + 1)+nlev] = SourceDate[ele*Np2d * 2 * nlev + Np2d + p];
		 }
	 }
 }


 void CalculateBuoyanceFrequencyDate(double *hT, double *hS, double hcrit, int K2d, \
	 int Np2d, int Np3d, int nlev, double gra, double rho0, double *J2d, double *wq2d,\
	 double *Vq2d, ptrdiff_t RVq2d, ptrdiff_t Cvq2d, double *LAV, char *type, double T0,\
	 double S0, double alphaT, double betaS){
	 memset(buoyanceFrequencyDate, 0.0, K2d*Np2d*(nlev + 1) * sizeof(double));
/*
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	 for (int i = 0; i < K2d; i++) {
		 //Values at the bottom surface and the upper surface at the same location
		 double TBottomAve, TUpAve, SBottomAve, SUpAve;
		 //The average T and S at the interface
		 double TSurfaceAve, SSurfaceAve;
		 //Density calculated according to the upper elemental information and lower elemental information
		 double rhoUp, rhoDown;
		 //Contribution to buoyance due to T and S
		 double NNT, NNS;

		 if (hcenter[i] >= hcrit) {

			 for (int L = 1; L < nlev; L++) {

				 GetElementCentralData(&TBottomAve, hT + i*nlev*Np3d + (nlev - L + 1)*Np3d - Np2d, \
					 J2d + i*Np2d, wq2d, Vq2d, RVq2d, Cvq2d, LAV + i);

				 GetElementCentralData(&TUpAve, hT + i*nlev*Np3d + (nlev - L - 1)*Np3d, \
					 J2d + i*Np2d, wq2d, Vq2d, RVq2d, Cvq2d, LAV + i);

				 TSurfaceAve = max((TBottomAve + TUpAve) / 2.0, 0.0);

				 if (!strcmp(type, "Jackett05")) {
					 EosByFeistel(&rhoUp, TSurfaceAve / hcenter[i], max(hSVerticalLine[i*(nlev + 1) + L + 1] / hcenter[i], 0.0));

					 EosByFeistel(&rhoDown, TSurfaceAve / hcenter[i], max(hSVerticalLine[i*(nlev + 1) + L] / hcenter[i], 0.0));
				 }
				 else if (!strcmp(type, "UNESCO83")) {
					 EosByUNESCO(&rhoUp, TSurfaceAve / hcenter[i], max(hSVerticalLine[i*(nlev + 1) + L + 1] / hcenter[i], 0.0));

					 EosByUNESCO(&rhoDown, TSurfaceAve / hcenter[i], max(hSVerticalLine[i*(nlev + 1) + L] / hcenter[i], 0.0));
				 }
				 else if(!strcmp(type, "Linear")) {
					 EosByLinear(&rhoUp, TSurfaceAve / hcenter[i], max(hSVerticalLine[i*(nlev + 1) + L + 1] / hcenter[i], 0.0), rho0, T0, S0, alphaT, betaS);

					 EosByLinear(&rhoDown, TSurfaceAve / hcenter[i], max(hSVerticalLine[i*(nlev + 1) + L] / hcenter[i], 0.0), rho0, T0, S0, alphaT, betaS);
				 }
				 else {
					 printf("Equation of state(EOS) needs to be pointed for this part!\n");
					 break;
				 }
				 
				 NNT = -1.0*gra / rho0*(rhoUp - rhoDown) / (0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]));

				 GetElementCentralData(&SBottomAve, hS + i*nlev*Np3d + (nlev - L + 1)*Np3d - Np2d, \
					 J2d + i*Np2d, wq2d, Vq2d, RVq2d, Cvq2d, LAV + i);

				 GetElementCentralData(&SUpAve, hS + i*nlev*Np3d + (nlev - L - 1)*Np3d, \
					 J2d + i*Np2d, wq2d, Vq2d, RVq2d, Cvq2d, LAV + i);

				 SSurfaceAve = max((SBottomAve + SUpAve) / 2.0, 0.0);

				 if (!strcmp(type, "Jackett05")) {
					 EosByFeistel(&rhoUp, max(hTVerticalLine[i*(nlev + 1) + L + 1] / hcenter[i], 0.0), SSurfaceAve / hcenter[i]);

					 EosByFeistel(&rhoDown, max(hTVerticalLine[i*(nlev + 1) + L] / hcenter[i], 0.0), SSurfaceAve / hcenter[i]);
				 }
				 else if (!strcmp(type, "UNESCO83")) {
					 EosByUNESCO(&rhoUp, max(hTVerticalLine[i*(nlev + 1) + L + 1] / hcenter[i], 0.0), SSurfaceAve / hcenter[i]);

					 EosByUNESCO(&rhoDown, max(hTVerticalLine[i*(nlev + 1) + L] / hcenter[i], 0.0), SSurfaceAve / hcenter[i]);
				 }
				 else if (!strcmp(type, "Linear")) {
					 EosByLinear(&rhoUp, max(hTVerticalLine[i*(nlev + 1) + L + 1] / hcenter[i], 0.0), SSurfaceAve / hcenter[i], rho0, T0, S0, alphaT, betaS);

					 EosByLinear(&rhoDown, max(hTVerticalLine[i*(nlev + 1) + L] / hcenter[i], 0.0), SSurfaceAve / hcenter[i], rho0, T0, S0, alphaT, betaS);
				 }
				 else {
					 printf("Equation of state(EOS) needs to be pointed for this part!\n");
					 break;
				 }

				 NNS = -1.0*gra / rho0*(rhoUp - rhoDown) / (0.5*(layerHeight[i*(nlev + 1) + L + 1] + layerHeight[i*(nlev + 1) + L]));

				 buoyanceFrequencyDate[i*(nlev + 1) + L] = NNT + NNS;

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
	 */
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
	 (*dest) = 0.0;
	 GetMeshAverageValue(dest, LAV, transA, transB, &RVq, &COPB, &CVq, &Alpha, Vq, \
		 &LDA, source, Jacobian, &LDB, &Beta, &LDC, wq);
 }
