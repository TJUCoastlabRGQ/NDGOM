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
extern int Np2d;
extern int K2d;
double Date[K2d][Np2d];
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
        for(int i = 0; i<K2d;i++){
			for (int j = 0; j < Np2d;j++)
                Date[i][j] = i*K2d+j;
        }
		Initialized = "True";
	}
    else{
        mexPrintf("The turbulence model has been initialized already\n");
        for(int i = 0; i<K2d*Np2d;i++)
            mexPrintf("The %dth date is %f",i,*(Date[0]+j));
    }
}