#include "mxGOTM.h"
#include <string.h>
#include <math.h>
#include "../../../../../NdgMath/NdgMemory.h"

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
	int Num2d = Np2d*K2d;
	int Interface = (int)nlev + 1;
//	int NMaxItration = 3;

	if (!strcmp("False", GOTMInitialized)){
		GotmSolverMemoryAllocation(Num2d, Interface, Np2d, K3d);
		/* Get number of characters in the input string.  Allocate enough
		memory to hold the converted string. */
		long long int buflen = (long long int)mxGetN(prhs[10]) + 1;
		char *buf = malloc(buflen);
		/* Copy the string data into buf. */
		mxGetString(prhs[10], buf, (mwSize)buflen);
		long long int _nNamelist = 2;
		InitTurbulenceModelGOTM(&_nNamelist, buf, buflen, nlev, Np2d, K2d);
		free(buf);
	}
		plhs[0] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
		plhs[2] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		//plhs[4] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		//plhs[5] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		//plhs[6] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);

		//mapVedgeDateToDof(eddyTKEDate, PtrOutTKE, Np2d, K2d, Np3d, nlev);

		//mapVedgeDateToDof(eddyLengthDate, PtrOutLength, Np2d, K2d, Np3d, nlev);

		//mapVedgeDateToDof(eddyEPSDate, PtrOutEPS, Np2d, K2d, Np3d, nlev);
		
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

		InterpolationToCentralPoint(hu, huCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV);

		InterpolationToCentralPoint(hv, hvCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV);
		/*The gradient about rho in vertical direction is calculated according to rho directly, not T and S.
		Details about the latter manner can be found in Tuomas and Vincent(2012, Ocean modelling)
		*/
		InterpolationToCentralPoint(rho, rhoCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV);

		//Tc to be continued
		//Sc to be continued
		mapCentralPointDateToVerticalDate(huCentralDate, huVerticalLine, K2d, nlev, Np2d);

		mapCentralPointDateToVerticalDate(hvCentralDate, hvVerticalLine, K2d, nlev, Np2d);
		/*The gradient about rho in vertical direction is calculated according to rho directly, not T and S*/
		mapCentralPointDateToVerticalDate(rhoCentralDate, rhoVerticalLine, K2d, nlev, Np2d);
		//Tvl to be continued
		//Svl to be continued
		CalculateWaterDepth(h, Np2d, K2d, hcrit, nlev);

		CalculateShearFrequencyDate(h, Np2d, K2d, hcrit, nlev);

		CalculateBuoyanceFrequencyDate(h, Np2d, K2d, hcrit, nlev, gra, rho0);

		CalculateLengthScaleAndShearVelocity(z0b, z0s, h, hcrit, PtrOutDragCoefficient, WindTaux, WindTauy, Np2d, K2d, nlev);

		DGDoTurbulence(&dt, h, hcrit, NULL, Np2d, K2d, nlev);

		mapVedgeDateToDof(eddyViscosityDate, PtrOutEddyViscosity, Np2d, K2d, Np3d, nlev);
/* If the following parts are to be exported, the corresponding parts in DGDoTurbulence in file mxGOTM.c need to be activated*/
//		mapVedgeDateToDof(eddyDiffusionDate, PtrOutDiffusionCoeForT, Np2d, K2d, Np3d, nlev);

//		mapVedgeDateToDof(eddyDiffusionDate, PtrOutDiffusionCoeForS, Np2d, K2d, Np3d, nlev);

		mapVedgeDateToDof(eddyTKEDate, PtrOutTKE, Np2d, K2d, Np3d, nlev);

//		mapVedgeDateToDof(eddyLengthDate, PtrOutLength, Np2d, K2d, Np3d, nlev);

		mapVedgeDateToDof(eddyEPSDate, PtrOutEPS, Np2d, K2d, Np3d, nlev);

}