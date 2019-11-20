#include "mxGOTM.h"
#include<string.h>
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

char *Initialized = "False";
int Np2d, K2d, Np3d, K3d;
long long int nlev;
double hcrit, finalTime, h0b;
double *VCV = NULL;
double *tkeGOTM = NULL;
double *epsGOTM = NULL;
double *LGOTM = NULL;
double *nuhGOTM = NULL;
double *numGOTM = NULL;
double *layerHeight = NULL;
double *huCentralDate = NULL;
double *hvCentralDate = NULL;
double *huVerticalLine = NULL;
double *hvVerticalLine = NULL;
double *shearFrequencyDate = NULL;
double *buoyanceFrequencyDate = NULL;
double *BottomFrictionLength = NULL;
double *BottomFrictionVelocity = NULL;
double *SurfaceFrictionLength = NULL;
double *SurfaceFrictionVelocity = NULL;
double *eddyViscosityDate = NULL;

void memoryCheck(double *Ptr) {
	if (Ptr == NULL){
		mexPrintf("Memory allocation failed");
		exit(1);
	}
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* check input & output */
	//if (nrhs != NRHS) {
	//	mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
	//	mexPrintf("%d inputs required.\n", NRHS);
	//}

	//if (nlhs != NLHS) {
	//	mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
	//	mexPrintf("%d inputs required.\n", NLHS);
	//}

	if (strcmp("True", Initialized)){
		mexPrintf("The turbulence model has not been initialized yet, and we are ready to do it\n");
		/*Memory allocation part*/
		Np2d = (int)mxGetScalar(prhs[0]);
		K2d = (int)mxGetScalar(prhs[1]);
		Np3d = (int)mxGetScalar(prhs[2]);
		K3d = (int)mxGetScalar(prhs[3]);
		nlev = (long long int)mxGetScalar(prhs[4]);
		hcrit = mxGetScalar(prhs[5]);
		finalTime = mxGetScalar(prhs[6]);
		h0b = mxGetScalar(prhs[14]);
		VCV = (double*)malloc(sizeof(double)*(int)mxGetNumberOfElements(prhs[7]));
		double *PtrVCV = mxGetPr(prhs[7]);
		for (int i = 0; i < (int)mxGetNumberOfElements(prhs[7]); i++)
			VCV[i] = PtrVCV[i];
		memoryCheck(VCV);
		int Num2d = Np2d*K2d;
		int Interface = (int)nlev + 1;
		tkeGOTM = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(tkeGOTM);
		epsGOTM = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(epsGOTM);
		LGOTM = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(LGOTM);
		nuhGOTM = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(nuhGOTM);
		numGOTM = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(numGOTM);
		layerHeight = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(layerHeight);
		huCentralDate = (double*)malloc(sizeof(double)*(Np2d*K3d)); memoryCheck(huCentralDate);
		hvCentralDate = (double*)malloc(sizeof(double)*(Np2d*K3d)); memoryCheck(hvCentralDate);
		huVerticalLine = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(huVerticalLine);
		hvVerticalLine = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(hvVerticalLine);
		shearFrequencyDate = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(shearFrequencyDate);
		buoyanceFrequencyDate = (double*)malloc(sizeof(double)*(Num2d*Interface)); memoryCheck(buoyanceFrequencyDate);
		BottomFrictionLength = (double*)malloc(sizeof(double)*Num2d); memoryCheck(BottomFrictionLength);
		BottomFrictionVelocity = (double*)malloc(sizeof(double)*Num2d); memoryCheck(BottomFrictionVelocity);
		SurfaceFrictionLength = (double*)malloc(sizeof(double)*Num2d); memoryCheck(SurfaceFrictionLength);
		SurfaceFrictionVelocity = (double*)malloc(sizeof(double)*Num2d); memoryCheck(SurfaceFrictionVelocity);
		eddyViscosityDate = (double*)malloc(sizeof(double)*(Num2d * Interface)); memoryCheck(eddyViscosityDate);

		char *buf;
		long long int buflen;
		int status;
		if (!mxIsChar(prhs[11]) || (mxGetM(prhs[11]) != 1))  {
			mexErrMsgIdAndTxt("MATLAB:mxmalloc:invalidInput",
				"Input argument must be the name of the file that contain the information of the turbulence model.");
		}
		/* Get number of characters in the input string.  Allocate enough
		memory to hold the converted string. */
		buflen = (long long int)mxGetN(prhs[11]) + 1;
		buf = mxMalloc(buflen);
		/* Copy the string data into buf. */
		status = mxGetString(prhs[11], buf, (mwSize)buflen);
		long long int _nNamelist = 2;
		InitTurbulenceModelGOTM(&_nNamelist, buf, buflen - 1);
		Initialized = "True";
	}
	else if (finalTime = mxGetScalar(prhs[13])){
		free(VCV);
		free(nuhGOTM); free(numGOTM); free(tkeGOTM); free(epsGOTM); free(LGOTM);
		free(layerHeight);
		free(huCentralDate); free(hvCentralDate);
		free(huVerticalLine); free(hvVerticalLine);
		free(shearFrequencyDate); free(buoyanceFrequencyDate);
		free(BottomFrictionLength); free(BottomFrictionVelocity);
		free(SurfaceFrictionLength);free(SurfaceFrictionVelocity);
		free(eddyViscosityDate);
	}
	else
	{
		double* h = mxGetPr(prhs[8]);
		double* hu = mxGetPr(prhs[9]);
		double* hv = mxGetPr(prhs[10]);
		double dt = (double)mxGetScalar(prhs[12]);
		plhs[0] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
		double *PtrOutEddyViscosity = mxGetPr(plhs[0]);
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

		CalculateLengthScaleAndShearVelocity(h);

		DGDoTurbulence(&dt, h, NULL);

		mapVedgeDateToDof(eddyViscosityDate, PtrOutEddyViscosity);
	}
}