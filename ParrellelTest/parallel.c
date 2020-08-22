#include "mex.h"

#include <omp.h>


void SumForNumber()
{
    int a = 0;
    for (int i = 0; i<100000000; i++)
        a++;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#pragma omp parallel for num_threads(12)
    for (int i = 0; i<100; i++){
        SumForNumber();
    }
        
  return;
}