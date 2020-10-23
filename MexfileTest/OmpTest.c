#include "mex.h"
#include <stdio.h>
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	printf("max threads = %d\n", omp_get_max_threads());
#pragma omp parallel
	{
		printf("ID = %d\n", omp_get_thread_num());
		printf("nThreads = %d\n", omp_get_num_threads());
	}
	printf("End\n");
	return;

}