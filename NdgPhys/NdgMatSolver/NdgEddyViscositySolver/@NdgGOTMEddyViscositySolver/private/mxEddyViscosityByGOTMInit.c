#include "mxGOTM.h"

#define NLHS 0
#define NRHS 9

//mxEddyViscosityByGOTMInit(mesh2d(1).K, mesh2d(1).cell.Np, mesh3d(1).K, mesh3d(1).cell.Np, ...
//	mesh3d(1).Nz, hcrit, obj.h0b, physClass.finalTime);



/*
* Purpose: This function is used to initialize all the global variables contained in file GOTM.h, and also initialize GOTM
*
* Input:
*      int[1] K2d element number of the two-dimensional mesh
* 	   int[1] Np2d number of intepolation points for the two-dimensional standard cell
*      int[1] K3d element number of the three-dimensional mesh
* 	   int[1] Np3d number of intepolation points for the three-dimensional standard cell
* 	   long long int[1] nlev number of vertical layers
*      double[1] hcrit the threshold value of water depth for wet-dry cell
*      double[1] h0b the physical bottom roughness length 
*      double[1] finalTime the last time of the whole simulation
*      char*     GOTMNameList the file containing the variable definition of the turbulence model
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

	int K2d = (int)mxGetScalar(prhs[0]);
	int Np2d = (int)mxGetScalar(prhs[1]);
	int K3d = (int)mxGetScalar(prhs[2]);
	int Np3d = (int)mxGetScalar(prhs[3]);
	long long int nlev = (long long int)mxGetScalar(prhs[4]);
	double hcrit = (double)mxGetScalar(prhs[5]);
	double h0b = (double)mxGetScalar(prhs[6]);
	double finalTime = (double)mxGetScalar(prhs[7]);
	

	char *buf;
	long long int buflen;
	int status;

	if (!mxIsChar(prhs[8]) || (mxGetM(prhs[8]) != 1))  {
		mexErrMsgIdAndTxt("MATLAB:mxmalloc:invalidInput",
			"Input argument must be the name of the file that contain the information of the turbulence model.");
	}

	/* Get number of characters in the input string.  Allocate enough
	memory to hold the converted string. */

	buflen = (long long int)mxGetN(prhs[8]) + 1;
	buf = mxMalloc(buflen);

	/* Copy the string data into buf. */
	status = mxGetString(prhs[8], buf, (mwSize)buflen);

	long long int  _nNamelist = 2;

	if (!strcmp("True", Initialized)){
		mexPrintf("The turbulence model has been initialized already\n");
	}
	else{
		mexPrintf("The turbulence model is not initialized, and we are ready to initialize it\n");
		InitTurbulenceModelGOTM(&_nNamelist, buf, &nlev, buflen);
		Initialized = "True";
	}


}