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
	double z0s = mxGetScalar(prhs[14]);
	double z0b = mxGetScalar(prhs[15]);
	double gra = mxGetScalar(prhs[16]);
	double rho0 = mxGetScalar(prhs[17]);
	int Interface = (int)nlev + 1;
	double *J2d = mxGetPr(prhs[18]);
	double *wq2d = mxGetPr(prhs[19]);
	double *Vq2d = mxGetPr(prhs[20]);
	int RVq2d = (int)mxGetM(prhs[20]);
	int CVq2d = (int)mxGetN(prhs[20]);
	double *LAV2d = mxGetPr(prhs[21]);
	double *huNew = mxGetPr(prhs[22]);
	double *hvNew = mxGetPr(prhs[23]);
	double *J3d = mxGetPr(prhs[24]);
	double *wq3d = mxGetPr(prhs[25]);
	double *Vq3d = mxGetPr(prhs[26]);
	int RVq3d = (int)mxGetM(prhs[26]);
	int CVq3d = (int)mxGetN(prhs[26]);
	double *LAV3d = mxGetPr(prhs[27]);
	double *hT = mxGetPr(prhs[28]);
	double *hS = mxGetPr(prhs[29]);
	double T0 = mxGetScalar(prhs[30]);
	double S0 = mxGetScalar(prhs[31]);
	double alphaT = mxGetScalar(prhs[32]);
	double betaS = mxGetScalar(prhs[33]);
	char* EosType;
	EosType = mxArrayToString(prhs[34]);

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
		plhs[4] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		double *PtrOutEddyViscosity = mxGetPr(plhs[0]);
        double *PtrOutDragCoefficient = mxGetPr(plhs[1]);
		/*Turbulence kenetic energy*/
		double *PtrOutTKE = mxGetPr(plhs[2]);
		/*Turbulence kinetic energy dissipation rate*/
		double *PtrOutEPS = mxGetPr(plhs[3]);
		/*DiffusionParameter for T and S*/
		double *PtrOutDiffusionCoeForT = mxGetPr(plhs[4]);

		ptrdiff_t TempNp2d = (ptrdiff_t)Np2d;

		ptrdiff_t TempNp3d = (ptrdiff_t)Np3d;

		ptrdiff_t TempK3d = (ptrdiff_t)K3d;

		/*Calculate the water depth at cell center first*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int i = 0; i < K2d; i++) {
			GetElementCentralData(hcenter + i, h + i*Np2d, J2d + i*Np2d, wq2d, Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d + i);
		}

		InterpolationToCentralPoint(hu, huCentralDate, K2d, Np2d, Np3d, (int)nlev, J2d, wq2d, Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d);
		
		InterpolationToCentralPoint(huNew, huCentralDateNew, K2d, Np2d, Np3d, (int)nlev, J2d, wq2d, Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d);

		InterpolationToCentralPoint(hv, hvCentralDate, K2d, Np2d, Np3d, (int)nlev, J2d, wq2d, Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d);
		
		InterpolationToCentralPoint(hvNew, hvCentralDateNew, K2d, Np2d, Np3d, (int)nlev, J2d, wq2d, Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d);

		InterpolationToCentralPoint(hT, hTCentralData, K2d, Np2d, Np3d, (int)nlev, J2d, wq2d, Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d);

		InterpolationToCentralPoint(hS, hSCentralData, K2d, Np2d, Np3d, (int)nlev, J2d, wq2d, Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d);

		mapCentralPointDateToVerticalDate(huCentralDate, huVerticalLine, K2d, (int)nlev);

		mapCentralPointDateToVerticalDate(huCentralDateNew, huVerticalLineNew, K2d, (int)nlev);

		mapCentralPointDateToVerticalDate(hvCentralDate, hvVerticalLine, K2d, (int)nlev);

		mapCentralPointDateToVerticalDate(hvCentralDateNew, hvVerticalLineNew, K2d, (int)nlev);

		/*The gradient about rho in vertical direction is calculated according to T and S*/
		mapCentralPointDateToVerticalDate(hTCentralData, hTVerticalLine, K2d, (int)nlev);

		mapCentralPointDateToVerticalDate(hSCentralData, hSVerticalLine, K2d, (int)nlev);

		CalculateWaterDepth(K2d, hcrit, (int)nlev);

		CalculateShearFrequencyDate(K2d, hcrit, (int)nlev);

		CalculateBuoyanceFrequencyDate(hT, hS, hcrit, K2d,  Np2d, Np3d, (int)nlev, gra, rho0, J2d, wq2d, \
			Vq2d, (ptrdiff_t)RVq2d, (ptrdiff_t)CVq2d, LAV2d, EosType, T0, S0, alphaT, betaS);

		CalculateLengthScaleAndShearVelocity(z0b, z0s, hcrit, PtrOutDragCoefficient, WindTaux, WindTauy, Np2d, K2d, (int)nlev);

		DGDoTurbulence(&dt, hcrit, NULL, K2d, nlev);

		mapVedgeDateToDof(numGOTM, PtrOutEddyViscosity, Np2d, K2d, Np3d, (int)nlev);

		mapVedgeDateToDof(nuhGOTM, PtrOutDiffusionCoeForT, Np2d, K2d, Np3d, (int)nlev);

		mapVedgeDateToDof(tkeGOTM, PtrOutTKE, Np2d, K2d, Np3d, (int)nlev);

//		mapVedgeDateToDof(LGOTM, PtrOutLength, Np2d, K2d, Np3d, nlev);

		mapVedgeDateToDof(epsGOTM, PtrOutEPS, Np2d, K2d, Np3d, (int)nlev);

}