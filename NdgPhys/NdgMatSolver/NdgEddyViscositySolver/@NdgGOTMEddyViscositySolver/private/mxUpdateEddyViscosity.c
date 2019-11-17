#include "mxGOTM.h"
//#include <blas.h>

#define NLHS 1
#define NRHS 6

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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* check input & output */
	if (nrhs != NRHS) {
		mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
		mexPrintf("%d inputs required.\n", NRHS);
	}

	if (nlhs != NLHS) {
		mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
		mexPrintf("%d inputs required.\n", NLHS);
	}
	double* VCV = mxGetPr(prhs[0]);
	double* h = mxGetPr(prhs[1]);
	double* hu = mxGetPr(prhs[2]);
	double* hv = mxGetPr(prhs[3]);
	double dt = (double)mxGetScalar(prhs[4]);
	double time = (double)mxGetScalar(prhs[5]);
	if (time == finalTime){
		mexPrintf("This is the end of the simulation, and we are ready to dealocation the memory");
		// to be continued
	else//Update the eddy viscosity
		mxArray *huCentralDate = mxCreateDoubleMatrix(Np2d, K3d, mxREAL);
	    double *huC = mxGetPr(uCentralDate);
	    mxArray *hvCentralDate = mxCreateDoubleMatrix(Np2d, K3d, mxREAL);
		double *hvC = mxGetPr(vCentralDate);
		//S and T to be continued
		mxArray *huVertialLine = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		double *huvl = mxGetPr(uVertialLine);
		mxArray *hvVertialLine = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		double *hvvl = mxGetPr(vVertialLine);
		//S and T to be continued
		mxArray *shearFrequencyDate = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		double *shearFrequency = mxGetPr(shearFrequencyDate);
		//Buoyance frequency to be continued
		mxArray *eddyViscosityDate = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		double *eddyViscosity = mxGetPr(eddyViscosityDate);
		// TDiffusionParameter to be continued
		// SDiffusionParameter to be continued

		InterpolationToCentralPoint(VCV, hu, huC);
		InterpolationToCentralPoint(VCV, hv, hvC);
		//Tc to be continued
		//Sc to be continued
		mapCentralPointDateToVerticalDate(huC, huvl);
		mapCentralPointDateToVerticalDate(hvC, hvvl);
		//Tuvl to be continued
		//Suvl to be continued



	}
	mxDestroyArray(uCentralDate);
	mxDestroyArray(vCentralDate);
	mxDestroyArray(uVertialLine);
	mxDestroyArray(vVertialLine);
	mxDestroyArray(shearFrequencyDate);
//to do destroy the eddy viscosity array after output

}
