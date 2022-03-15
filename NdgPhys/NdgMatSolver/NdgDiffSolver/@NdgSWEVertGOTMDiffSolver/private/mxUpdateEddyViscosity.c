#include "mxGOTM.h"
#include <string.h>
#include <math.h>
#include "../../../../../NdgMath/NdgMemory.h"
#include "../../../../../NdgMath/NdgMath.h"

/*
* Purpose: This function is used to update the eddy viscosity according to the flow field
*
* Input:
*      double[Np2d x Np3d] VCV, the interpolation matrix used to calculate the value at the central point
* 	   double[Np2d x K2d] WaterDepth, water depth field
*      double[Np3d x K3d] hu, the three dimensional horizontal momentum field in x direction
* 	   double[Np3d x K3d] hv, the three dimensional horizontal momentum field in y direction
* 	   double[1] dt, time step, the time step used by the three-dimensional hydraulic model
*      double[1] time, time point, the time that the simulation has arrived at
* Output:
*      double[Np3d x K3d], nu, the updated eddy viscosity
*/

void MyExit()
{
	if (!strcmp("True", GOTMInitialized)){
		GotmSolverMemoryDeAllocation();
		TURBULENCE_mp_CLEAN_TURBULENCE();
		MTRIDIAGONAL_mp_CLEAN_TRIDIAGONAL();
	}
	return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexAtExit(&MyExit);
	int Np2d = (int)mxGetScalar(prhs[0]);
	int K2d = (int)mxGetScalar(prhs[1]);
	int Np3d = (int)mxGetScalar(prhs[2]);
	int K3d = (int)mxGetScalar(prhs[3]);
	long long int nlev = (long long int)mxGetScalar(prhs[4]);
	double hcrit = mxGetScalar(prhs[5]);
	double *VCV = mxGetPr(prhs[6]);
	double* h = mxGetPr(prhs[7]);
	double* hu = mxGetPr(prhs[8]);
	double* hv = mxGetPr(prhs[9]);
	double dt = mxGetScalar(prhs[11]);
	double* WindTaux = mxGetPr(prhs[12]);
	double* WindTauy = mxGetPr(prhs[13]);
	double* rho = mxGetPr(prhs[14]);
	double z0s = mxGetScalar(prhs[15]);
	double z0b = mxGetScalar(prhs[16]);
	double gra = mxGetScalar(prhs[17]);
	double rho0 = mxGetScalar(prhs[18]);
	int Interface = (int)nlev + 1;
	double *J2d = mxGetPr(prhs[19]);
	double *wq2d = mxGetPr(prhs[20]);
	double *Vq2d = mxGetPr(prhs[21]);
	int RVq2d = (int)mxGetM(prhs[21]);
	int CVq2d = (int)mxGetN(prhs[21]);
	double *LAV2d = mxGetPr(prhs[22]);
	double *huNew = mxGetPr(prhs[23]);
	double *hvNew = mxGetPr(prhs[24]);

	double *J3d = mxGetPr(prhs[25]);
	double *wq3d = mxGetPr(prhs[26]);
	double *Vq3d = mxGetPr(prhs[27]);
	int RVq3d = (int)mxGetM(prhs[27]);
	int CVq3d = (int)mxGetN(prhs[27]);
	double *LAV3d = mxGetPr(prhs[28]);

//	int NMaxItration = 3;

	if (!strcmp("False", GOTMInitialized)){
		GotmSolverMemoryAllocation(Interface, Np2d, K3d, K2d);
		/* Get number of characters in the input string.  Allocate enough
		memory to hold the converted string. */
		long long int buflen = (long long int)mxGetN(prhs[10]) + 1;
		char *buf = malloc(buflen);
		/* Copy the string data into buf. */
		mxGetString(prhs[10], buf, (mwSize)buflen);
		long long int _nNamelist = 2;
		InitTurbulenceModelGOTM(&_nNamelist, buf, buflen, nlev, K2d);
		free(buf);
	}
		plhs[0] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
		plhs[2] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		//plhs[4] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		//plhs[5] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		//plhs[6] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
//		plhs[4] = mxCreateDoubleMatrix((int)nlev + 1, K2d, mxREAL);
//		double *OutShearFrequencyData = mxGetPr(plhs[4]);
//		plhs[4] = mxCreateDoubleMatrix(K2d, 1, mxREAL);
//		double *Hcenter = mxGetPr(plhs[4]);
//		plhs[5] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
//		double *PtrOutLength = mxGetPr(plhs[5]);

		double *PtrOutEddyViscosity = mxGetPr(plhs[0]);
        double *PtrOutDragCoefficient = mxGetPr(plhs[1]);
		/*DiffusionParameter for T*/
//		double *PtrOutDiffusionCoeForT = mxGetPr(plhs[2]);
        /*DiffusionParameter for S*/
//		double *PtrOutDiffusionCoeForS = mxGetPr(plhs[3]);
		/*Turbulence kenetic energy*/
		double *PtrOutTKE = mxGetPr(plhs[2]);
		/*Turbulence kenetic length*/
		//double *PtrOutLength = mxGetPr(plhs[3]);
		/*Turbulence kinetic energy dissipation rate*/
		double *PtrOutEPS = mxGetPr(plhs[3]);

		ptrdiff_t TempNp2d = (ptrdiff_t)Np2d;
		ptrdiff_t TempNp3d = (ptrdiff_t)Np3d;
		ptrdiff_t TempK3d = (ptrdiff_t)K3d;

		double *huCentralDateO = malloc(K3d * sizeof(double));
		double *huCentralDateNewO = malloc(K3d * sizeof(double));
		// Original version
		plhs[4] = mxCreateDoubleMatrix(1, K3d, mxREAL);
		double *UcentralOutO = mxGetPr(plhs[4]);
		// New Version
		plhs[5] = mxCreateDoubleMatrix(1, K3d, mxREAL);
		double *UcentralOut = mxGetPr(plhs[5]);

		/*Calculate the water depth at cell center first*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int i = 0; i < K2d; i++) {
			GetElementCentralData(hcenter + i, h + i*Np2d, J2d + i*Np2d, wq2d, Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d + i);
		}

//		memcpy(Hcenter, hcenter, K2d*sizeof(double));
		InterpolationToCentralPoint(hu, huCentralDate, K2d, Np3d, (int)nlev, J3d, wq3d, Vq3d, (ptrdiff_t)RVq3d, (ptrdiff_t)CVq3d, LAV3d);
		InterpolationToCentralPoint(huNew, huCentralDateNew, K2d, Np3d, (int)nlev, J3d, wq3d, Vq3d, (ptrdiff_t)RVq3d, (ptrdiff_t)CVq3d, LAV3d);

//		InterpolationToCentralPointO(hu, huCentralDateO, K2d, Np2d, Np3d, (int)nlev, J2d, wq2d, \
			Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d);

//		InterpolationToCentralPointO(huNew, huCentralDateNewO, K2d, Np2d, Np3d, (int)nlev, J2d, wq2d, \
			Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d);
/*
		for (int k = 0; k < K3d; k++) {
			printf("For element %d\n",k);
			printf("For huCentralDate, the difference is:%f\n",huCentralDate[k] - huCentralDateO[k]);
			printf("For huCentralDateNew, the difference is:%f\n", huCentralDateNew[k] - huCentralDateNewO[k]);
		}
		*/
		// Original version
		memcpy(UcentralOutO, huCentralDateNewO, K3d * sizeof(double));
		//New version, implemented in matlab
		memcpy(UcentralOut, huCentralDateNew, K3d * sizeof(double));

//		memcpy(hucenterOutput, huCentralDate, K3d*sizeof(double));

		InterpolationToCentralPoint(hv, hvCentralDate, K2d, Np3d, (int)nlev, J3d, wq3d, Vq3d, (ptrdiff_t)RVq3d, (ptrdiff_t)CVq3d, LAV3d);
		InterpolationToCentralPoint(hvNew, hvCentralDateNew, K2d, Np3d, (int)nlev, J3d, wq3d, Vq3d, (ptrdiff_t)RVq3d, (ptrdiff_t)CVq3d, LAV3d);

		/*The gradient about rho in vertical direction is calculated according to rho directly, not T and S.
		Details about the latter manner can be found in Tuomas and Vincent(2012, Ocean modelling)
		*/
		InterpolationToCentralPoint(rho, rhoCentralDate, K2d, Np3d, (int)nlev, J3d, wq3d, Vq3d, (ptrdiff_t)RVq3d, (ptrdiff_t)CVq3d, LAV3d);
		//Tc to be continued
		//Sc to be continued
		mapCentralPointDateToVerticalDate(huCentralDate, huVerticalLine, K2d, (int)nlev);
		mapCentralPointDateToVerticalDate(huCentralDateNew, huVerticalLineNew, K2d, (int)nlev);
//		memcpy(huVerticalLineData, huVerticalLine, K2d*((int)nlev + 1) * sizeof(double));

		mapCentralPointDateToVerticalDate(hvCentralDate, hvVerticalLine, K2d, (int)nlev);
		mapCentralPointDateToVerticalDate(hvCentralDateNew, hvVerticalLineNew, K2d, (int)nlev);
		/*The gradient about rho in vertical direction is calculated according to rho directly, not T and S*/
		mapCentralPointDateToVerticalDate(rhoCentralDate, rhoVerticalLine, K2d, (int)nlev);
		//Tvl to be continued
		//Svl to be continued
		CalculateWaterDepth(K2d, hcrit, (int)nlev);

		CalculateShearFrequencyDate(K2d, hcrit, (int)nlev);
		//Checked
//		memcpy(OutShearFrequencyData, shearFrequencyDate, K2d*((int)nlev + 1) * sizeof(double));

		CalculateBuoyanceFrequencyDate( Np2d, K2d, hcrit, (int)nlev, gra, rho0);

		CalculateLengthScaleAndShearVelocity(z0b, z0s, hcrit, PtrOutDragCoefficient, WindTaux, WindTauy, Np2d, K2d, (int)nlev);

		DGDoTurbulence(&dt, hcrit, NULL, K2d, nlev);

		mapVedgeDateToDof(numGOTM, PtrOutEddyViscosity, Np2d, K2d, Np3d, (int)nlev);
//		mapVedgeDateToDof(nuhGOTM, PtrOutDiffusionCoeForT, Np2d, K2d, Np3d, nlev);

//		mapVedgeDateToDof(numGOTM, PtrOutDiffusionCoeForS, Np2d, K2d, Np3d, nlev);

		mapVedgeDateToDof(tkeGOTM, PtrOutTKE, Np2d, K2d, Np3d, (int)nlev);

//		mapVedgeDateToDof(LGOTM, PtrOutLength, Np2d, K2d, Np3d, nlev);

		mapVedgeDateToDof(epsGOTM, PtrOutEPS, Np2d, K2d, Np3d, (int)nlev);

		free(huCentralDateO);

		free(huCentralDateNewO);

}