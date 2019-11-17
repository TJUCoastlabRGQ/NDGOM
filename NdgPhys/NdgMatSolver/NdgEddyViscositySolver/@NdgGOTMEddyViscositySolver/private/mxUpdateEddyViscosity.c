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
		mxArray *layerHeight = mxCreateDoubleMatrix(nlev + 1, Np2d, mxREAL);
	    double *PtrLayerHeight = mxGetPr(layerHeight);
		mxArray *huCentralDate = mxCreateDoubleMatrix(Np2d, K3d, mxREAL);
		double *PtrhuCentralDate = mxGetPr(huCentralDate);
	    mxArray *hvCentralDate = mxCreateDoubleMatrix(Np2d, K3d, mxREAL);
		double *PtrhvCentralDate = mxGetPr(hvCentralDate);
		//S and T to be continued
		mxArray *huVertialLine = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		double *PtrhuVertialLine = mxGetPr(huVertialLine);
		mxArray *hvVertialLine = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		double *PtrhvVertialLine = mxGetPr(hvVertialLine);
		//S and T to be continued
		mxArray *shearFrequencyDate = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		double *PtrshearFrequency = mxGetPr(shearFrequencyDate);
		mxArray *buoyanceFrequencyDate = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		mxArray *PtrbuoyanceFrequencyDate = mxGetPr(buoyanceFrequencyDate);
		//Buoyance frequency to be continued
		mxArray *BottomFrictionLength = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
		double *PtrBottomFritionLength = mxGetPr(BottoFrictionLength);

		mxArray *BottomFrictionVelocity = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
		double *PtrBottomFritionVelocity = mxGetPr(BottoFrictionVelocity);

		mxArray *SurfaceFrictionVelocity = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
		double *PtrSurfaceFritionVelocity = mxGetPr(SurfaceFrictionVelocity);

		mxArray *SurfaceFrictionLength = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
		double *PtrSurfaceFrictionLength = mxGetPr(Np2d, K2d, mxREAL);


		mxArray *eddyViscosityDate = mxCreateDoubleMatrix(nlev + 1, Np2d*K2d, mxREAL);
		double *PtreddyViscosity = mxGetPr(eddyViscosityDate);

		plhs[0] = mxCreateDoubleMatrix(Np3d, K3d);
		double *PtrOutEddyViscosity = mxGetPr(plhs[0]);
		// TDiffusionParameter to be continued
		// SDiffusionParameter to be continued

		InterpolationToCentralPoint(VCV, hu, PtrhuCentralDate);
		InterpolationToCentralPoint(VCV, hv, PtrhvCentralDate);
		//Tc to be continued
		//Sc to be continued
		mapCentralPointDateToVerticalDate(PtrhuCentralDate, PtrhuVertialLine);
		mapCentralPointDateToVerticalDate(PtrhvCentralDate, PtrhvVertialLine);
		//Tvl to be continued
		//Svl to be continued
		CalculateWaterDepth(double * h, double *PtrLayerHeight);

		CalculateShearFrequencyDate(double *h, double *PtrLayerHeight, double *PtrhuVertialLine, double *PtrhvVertialLine, double *PtrshearFrequency);

		CalculateLengthScaleAndShearVelocity(double *h, double* PtrLayerHeight, double *PtrhuVertialLine, double *PtrhvVertialLine, double* PtrBottomFritionLength, \
                 double *PtrBottomFritionVelocity, double *PtrSurfaceFrictionLength, double *PtrSurfaceFritionVelocity);

		DGDoTurbulence(&dt, h, PtrSurfaceFritionVelocity, PtrBottomFritionVelocity, PtrSurfaceFrictionLength, PtrBottomFritionLength, \
			PtrLayerHeight, PtrbuoyanceFrequencyDate, PtrshearFrequency, NULL, PtreddyViscosity);

		mapVedgeDateToDof(PtreddyViscosity, PtrOutEddyViscosity);
	}
	mxDestroyArray(layerHeight);
	mxDestroyArray(uCentralDate);
	mxDestroyArray(vCentralDate);
	mxDestroyArray(uVertialLine);
	mxDestroyArray(vVertialLine);
	mxDestroyArray(shearFrequencyDate);
	mxDestroyArray(BottomFrictionLength);
	mxDestroyArray(BottomFrictionVelocity);
	mxDestroyArray(SurfaceFrictionLength);
	mxDestroyArray(SurfaceFrictionVelocity);
	
//to do destroy the eddy viscosity array after output

}
