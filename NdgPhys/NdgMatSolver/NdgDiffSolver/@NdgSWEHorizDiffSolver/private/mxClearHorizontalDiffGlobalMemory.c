#include <mex.h>
#include "..\..\..\..\..\NdgMath\NdgMemory.h"
#include <string.h>

extern char *HorDiffInitialized;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	if (!strcmp("True", HorDiffInitialized)){
		HorizDiffMemoryDeAllocation();
		HorDiffInitialized = "False";
	}
	return;
}