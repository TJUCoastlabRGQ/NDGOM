
#include "mex.h"
#include <math.h>
#include "blas.h"
#include "stdio.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 7
#define NLHS 1

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

/*
 * Purpose: This function is used to limit the physical value to satisfy the maxmium principle
 * This function is programmed according to ( Philippe Delandmeter, 2017 ).
 *
 * Input:
 *      double[Np x K] fphys the physical field to be limited.
 * 		double[1 x K]  avar the average value of each cell for the physical field
 * 		double[Nface x K]  EToE the topological relation of the studied mesh
 * 		double[Npz x Npz]  V1d the inversed one-dimensional Vandmonde matrix, used to calculate the average value in vertical direction for each line
 *      int[1] Npz number of interpolation points in vertical direction
 * 		int[1] Nph number of interpolation points in horizontal direction
 *      double[Npz x Npz]  OV1d the one-dimensional Vandmonde matrix, the first value of this matrix is used to calculate the average value
 * Output:
 * 		double[Np x K] limfphys the limited physical field.
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* check input & output */
    if (nrhs != NRHS) {
        mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NRHS);
    }
    
    if (nlhs != NLHS) {
        mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NLHS);
    }
    
    /* get inputs */
    double *fphys = mxGetPr(prhs[0]);
    double *avar = mxGetPr(prhs[1]);
    double *EToE = mxGetPr(prhs[2]);
    /*This is the inversed one-dimensional Vandmande matrix*/
    double *V1d = mxGetPr(prhs[3]);
    const int Npz = (int)mxGetScalar(prhs[4]);
    const int Nph = (int)mxGetScalar(prhs[5]);
    /*This is the one-dimensional VandMande matrix*/
    double *OV1d = mxGetPr(prhs[6]);
    size_t Np = mxGetM(prhs[0]);
    size_t K = mxGetN(prhs[0]);
    
    size_t Nface = mxGetM(prhs[2]); // number of faces of the computation cell
    
    /* allocate output array */
    plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
    double *limfphys = mxGetPr(plhs[0]);
    
    const double one = 1;
    const double zero = 0;
    ptrdiff_t one_ptrdiff = 1;
    ptrdiff_t Np_ptrdiff = Npz;
    ptrdiff_t Nq_ptrdiff = Npz;
    char *tran = "N";
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    
    for (int k = 0; k < K; k++) {
        /*This part is used to determine the maxmum and minimum allowable value of the studied cell with index k */
        double amax = -1*pow(10,10), amin = pow(10,10);
        for (int f = 0; f < Nface; f++){
            /*Only the adjacent cell is considered, and the average value of the studied cell is not included*/
            if ((int)EToE[k*Nface + f]-1 != k){
                amax = max(amax, avar[(mwIndex)EToE[k*Nface + f]-1]);
                amin = min(amin, avar[(mwIndex)EToE[k*Nface + f]-1]);
            }
            else
                continue;
        }
        /*This part is used to decide whether the studied cell is problematic following Cockburn and Shu*/
        int flag = 0;
        for (int i = 0; i < Np; i++){
            
            if (fphys[k*Np + i] > amax || fphys[k*Np + i] < amin){
                flag = 1;
                break;
            }
        }
        
        if (flag == 0) //the studied cell satisfy the maximum principle
            continue;
        else{// the studied cell is problematic
            // Calculate the average value of the vertical line
            double *avL = (double*)malloc(sizeof(double) * Nph), \
                    *tempValue = (double*)malloc(sizeof(double) * Npz), \
                    *fmod = (double*)malloc(sizeof(double)*Npz);
            for (int i = 0; i < Nph; i++){
                for (int j = 0; j < Npz; j++){
                    //Fetch the original value in each line and store them in tempValue, from top to down
                    *(tempValue+j) = fphys[k*Np + Nph*(Npz - 1) + i - j*Nph];
					//*(tempValue + j) = fphys[k*Np + i - j*Nph];
                }
                // Calculate the corresponding mode coefficients for each vertical line
                dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, V1d,
                        &Nq_ptrdiff, tempValue, &Np_ptrdiff, &zero, fmod, &Nq_ptrdiff);
                // the final average data for each vertical line
                *( avL + i ) = fmod[0] * OV1d[0];
                //Next limit the corresponding date
                /*Calculate the slope parameter first*/
                double Lambda;
                if (fphys[k*Np + Nph * (Npz - 1) + i] > amax)
                    Lambda = (amax - avar[k]) / ( fphys[k*Np + Nph*(Npz - 1) + i] - avar[k] + pow(10,-10) );
                else if (fphys[k*Np + Nph * (Npz - 1) + i] < amin)
                    Lambda = (avar[k] - amin) / (avar[k] - fphys[k*Np + Nph*(Npz - 1) + i] + pow(10, -10));
                else
                    Lambda = 1;
                /*Limit the value in vertical direction*/
                for (int j = 0; j < Npz; j++){
                    fphys[k*Np + Nph*(Npz - 1) + i - j*Nph] = Lambda * fphys[k*Np + Nph*(Npz - 1) + i - j*Nph]\
                            + (1-Lambda) * avL[i];
                }
            }
            free(avL);
            free(tempValue);
            free(fmod);
        }
        /*Next to limit the whole computation cell follows Cockburn and Shu*/
        /*This part is used to decide whether the studied cell is problematic following Cockburn and Shu*/
        flag = 0;
        for (int i = 0; i < Np; i++){
            if (fphys[k*Np + i] > amax || fphys[k*Np + i] < amin){
                flag = 1;
                break;
            }
        }
        if (flag == 0)
            continue;
        else{
            double Lambda = 1;
            for (int i=0;i<Np;i++){
                if ( fphys[k*Np + i] > amax )
                    Lambda = min( Lambda, ( amax - avar[k] )/( fphys[k*Np + i] - avar[k] + pow(10, -10)) );
                else if(fphys[k*Np + i] < amin)
                    Lambda = min( Lambda, ( avar[k] - amin )/( avar[k] - fphys[k*Np + i] + pow(10, -10)) );
            }
            Lambda = max( 0, Lambda );
            for(int i=0; i<Np;i++){
                fphys[k*Np + i] = Lambda * fphys[k*Np + i] + (1 - Lambda)*avar[k];
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    /*Copy the limited value to the output*/
    for(int i=0; i<Np*K; i++)
        limfphys[i] = fphys[i];
}
