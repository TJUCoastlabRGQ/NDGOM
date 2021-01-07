#include "mxGOTM.h"
#include<string.h>
#include<math.h>
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
long long int nlev;
double hcrit, finalTime, h0b;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double time = mxGetScalar(prhs[13]);
	int Np2d = (int)mxGetScalar(prhs[0]);
	int K2d = (int)mxGetScalar(prhs[1]);
	int Np3d = (int)mxGetScalar(prhs[2]);
	int K3d = (int)mxGetScalar(prhs[3]);
	long long int nlev = (long long int)mxGetScalar(prhs[4]);
	double hcrit = mxGetScalar(prhs[5]);
	double finalTime = mxGetScalar(prhs[6]);
	double h0b = mxGetScalar(prhs[14]);
	double *VCV = mxGetPr(prhs[7]);
	int Num2d = Np2d*K2d;
	int Interface = (int)nlev + 1;

	if (!strcmp("False", GOTMInitialized)){
		GotmSolverMemoryAllocation(int Num2d, int Interface, int Np2d, int K2d, int K3d);
		/* Get number of characters in the input string.  Allocate enough
		memory to hold the converted string. */
		buflen = (long long int)mxGetN(prhs[11]) + 1;
		buf = mxMalloc(buflen);
		/* Copy the string data into buf. */
		status = mxGetString(prhs[11], buf, (mwSize)buflen);
		long long int _nNamelist = 2;
		InitTurbulenceModelGOTM(&_nNamelist, buf, buflen - 1);
	}
		/*Memory allocation part*/

		char *buf;
		long long int buflen;
		int status;

	/*Simulation is continuing*/
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

		InterpolationToCentralPoint(hu, huCentralDate, &TempNp2d, &TempK3d, &TempNp3d);
		InterpolationToCentralPoint(hv, hvCentralDate, &TempNp2d, &TempK3d, &TempNp3d);
		//Tc to be continued
		//Sc to be continued
		mapCentralPointDateToVerticalDate(huCentralDate, huVerticalLine);
		mapCentralPointDateToVerticalDate(hvCentralDate, hvVerticalLine);
		//Tvl to be continued
		//Svl to be continued
		CalculateWaterDepth(h);

		CalculateShearFrequencyDate(h);

		CalculateBuoyanceFrequencyDate();// At present, all are set to zero

		CalculateLengthScaleAndShearVelocity(h, PteOutDragCoefficient, WindTaux, WindTauy);

		DGDoTurbulence(&dt, h, NULL);

		mapVedgeDateToDof(eddyViscosityDate, PtrOutEddyViscosity);


}