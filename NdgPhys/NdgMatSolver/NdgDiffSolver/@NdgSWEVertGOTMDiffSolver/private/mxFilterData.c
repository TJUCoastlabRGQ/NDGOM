#include <string.h>
#include <math.h>
#include "mex.h"
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
	int Np2d = (int)mxGetScalar(prhs[0]);
	int K2d = (int)mxGetScalar(prhs[1]);
	int nlev = (int)mxGetScalar(prhs[2]);
    int Nz = (int)mxGetScalar(prhs[3]);
    double *fphys = mxGetPr(prhs[4]);
    mwSize DimOfRHS = mxGetNumberOfDimensions(prhs[4]);
    int Nvar;
    const size_t *PRHS;
    PRHS = mxGetDimensions(prhs[4]);
    if (DimOfRHS == 2)
        Nvar = 1;
    else{
        Nvar = (int)PRHS[2];
    }
    int Np = (int)PRHS[0];
    int K3d = (int)PRHS[1];
    const size_t NdimOut = 3;
    const mwSize dimOut[3] = {Np,K3d,Nvar};
    plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);    
    double *Outfphys = mxGetPr(plhs[0]);
    int i;
    		
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (i = 0; i < K2d; i++) {
		double TempData, *Tempfphys, *Finalfphys;
		int p, var, L;
		for (var = 0; var<Nvar; var++) {
			Tempfphys = fphys + var*Np*K3d + i*nlev*Np;
			Finalfphys = Outfphys + var*Np*K3d + i*nlev*Np;
			if (nlev >= 2) {
				//For the first layer
				//The top face
				for (p = Np - Np2d; p < Np; p++) {
					Finalfphys[p] = Tempfphys[p];
				}
				//The middle part
				for (p = Np2d; p < Nz*Np2d; p++) {
					Finalfphys[p] = Tempfphys[p];
				}
				//The bottom face
				for (p = 0; p < Np2d; p++) {
					TempData = 0.5 * Tempfphys[p] + \
						0.5*Tempfphys[Np + Np - Np2d + p];
					Finalfphys[p] = TempData;
					Finalfphys[Np + Np - Np2d + p] = TempData;
				}
				//For the middle layer
				for (L = 1; L < nlev - 1; L++) {
					for (p = 0; p < Np2d; p++) {
						//The bottom surface data
						TempData = 0.5 * Tempfphys[L*Np + p] + \
							0.5*Tempfphys[L*Np + Np + Np - Np2d + p];
						Finalfphys[L*Np + p] = TempData;
						Finalfphys[L*Np + Np + Np - Np2d + p] = TempData;
					}
					for (p = Np2d; p < Nz*Np2d; p++) {
						//The inner data
						Finalfphys[L*Np + p] = Tempfphys[L*Np + p];
					}
				}
				//For the bottom most layer
				for (L = nlev - 1; L < nlev; L++) {
					for (p = 0; p < Nz*Np2d; p++) {
						Finalfphys[L*Np + p] = Tempfphys[L*Np + p];
					}
				}
			}
			else {//nlev = 1
				for (p = 0; p<Np; p++) {
					Finalfphys[p] = Tempfphys[p];
				}
			}
		}
	}
}