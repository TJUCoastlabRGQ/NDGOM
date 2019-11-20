#include "mex.h"
#include <string.h>
#include "Date.h"

#define NLHS 0
#define NRHS 2



/*
* Purpose: This function is used to initialize all the global variables contained in file GOTM.h, and also initialize GOTM
*
* Input:
*      int[1] K2d element number of the two-dimensional mesh
* 	   int[1] Np2d number of intepolation points for the two-dimensional standNard cell
*      int[1] K3d element number of the three-dimensional mesh
* 	   int[1] Np3d number of intepolation points for the three-dimensional standard cell
* 	   long long int[1] nlev number of vertical layers
*      double[1] hcrit the threshold value of water depth for wet-dry cell
*      double[1] h0b the physical bottom roughness length 
*      double[1] finalTime the last time of the whole simulation
*      char*     GOTMNameList the file containing the variable definition of the turbulence model
*/
char *Initialized = "False";
 int Np2d;
 int K2d;
 int *TestDate;
//double Date[K2d][Np2d];
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
	if(strcmp("True", Initialized)){
		mexPrintf("The turbulence model is not initialized, and we are ready to initialize it\n");
	    Np2d = (int)mxGetScalar(prhs[0]);
		K2d = (int)mxGetScalar(prhs[1]);
		//int *TestDate;
		TestDate = (int *)malloc(sizeof(int)*Np2d*K2d);
        for(int i = 0; i<K2d*Np2d;i++){
			TestDate[i] = i;
        }
		PrintTest(TestDate);
		Initialized = "True";
	}
    else{
        mexPrintf("The turbulence model has been initialized already\n");
		PrintTest(TestDate);
        free(TestDate);
        Initialized = "False";
//        for(int i = 0; i<K2d*Np2d;i++)
//			mexPrintf("The %dth date is %f\n", i, TestDate[i]);
    }
}

void PrintTest(int *data){
	for (int i = 0; i<K2d*Np2d; i++)
		mexPrintf("The %dth date is %d\n", i, data[i]);
}