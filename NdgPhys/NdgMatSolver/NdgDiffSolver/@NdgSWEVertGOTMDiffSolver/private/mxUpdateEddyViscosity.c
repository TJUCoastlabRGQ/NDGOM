#include "mxGOTM.h"
#include <string.h>
#include <math.h>
#include "../../../../../NdgMath/NdgMemory.h"
//#include <blas.h>

//#define NLHS 1
//#define NRHS 6

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
	double h0b = mxGetScalar(prhs[14]);
	double *VCV = mxGetPr(prhs[7]);
	int Num2d = Np2d*K2d;
	int Interface = (int)nlev + 1;
	int NMaxItration = 3;

	if (!strcmp("False", GOTMInitialized)){
		GotmSolverMemoryAllocation(Num2d, Interface, Np2d, K3d);
		/* Get number of characters in the input string.  Allocate enough
		memory to hold the converted string. */
		long long int buflen = (long long int)mxGetN(prhs[11]) + 1;
		char *buf = malloc(buflen);
		/* Copy the string data into buf. */
		mxGetString(prhs[11], buf, (mwSize)buflen);
		long long int _nNamelist = 2;
		InitTurbulenceModelGOTM(&_nNamelist, buf, buflen, nlev, Np2d, K2d);
//		TURBULENCE_mp_CLEAN_TURBULENCE();
//		TURBULENCE_mp_INIT_TURBULENCE(&_nNamelist, buf, &nlev, buflen);
		free(buf);
	}
	
		double* h = mxGetPr(prhs[8]);
		double* hu = mxGetPr(prhs[9]);
		double* hv = mxGetPr(prhs[10]);
		double dt = (double)mxGetScalar(prhs[12]);
        double* WindTaux = mxGetPr(prhs[15]);
        double* WindTauy = mxGetPr(prhs[16]);
		
		plhs[0] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
		
		double *PtrOutEddyViscosity = mxGetPr(plhs[0]);
        double *PteOutDragCoefficient = mxGetPr(plhs[1]);
		// TDiffusionParameter to be continued
		// SDiffusionParameter to be continued

		ptrdiff_t TempNp2d = (ptrdiff_t)Np2d;
		ptrdiff_t TempNp3d = (ptrdiff_t)Np3d;
		ptrdiff_t TempK3d = (ptrdiff_t)K3d;

		InterpolationToCentralPoint(hu, huCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV);

		InterpolationToCentralPoint(hv, hvCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV);
		//Tc to be continued
		//Sc to be continued
		mapCentralPointDateToVerticalDate(huCentralDate, huVerticalLine, K2d, nlev, Np2d);
		mapCentralPointDateToVerticalDate(hvCentralDate, hvVerticalLine, K2d, nlev, Np2d);
		//Tvl to be continued
		//Svl to be continued
		CalculateWaterDepth(h, Np2d, K2d, hcrit, nlev);

		CalculateShearFrequencyDate(h, Np2d, K2d, hcrit, nlev);

		CalculateBuoyanceFrequencyDate(Np2d, K2d, nlev);// At present, all are set to zero

		CalculateLengthScaleAndShearVelocity(h, hcrit, PteOutDragCoefficient, WindTaux, WindTauy, Np2d, K2d, NMaxItration, h0b, nlev);

		DGDoTurbulence(&dt, h, hcrit, NULL, Np2d, K2d, nlev);

		mapVedgeDateToDof(eddyViscosityDate, PtrOutEddyViscosity, Np2d, K2d, Np3d, nlev);

}