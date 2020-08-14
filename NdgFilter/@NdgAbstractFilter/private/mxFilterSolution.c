#include "mex.h"
#include <math.h>
#include "blas.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 9
#define NLHS 1

#define EPS pow(10,-10)

/*
 * Purpose: This function is used to filter the physical value to eliminate unphysical oscillation
 * This function is programmed according to ( Philippe Delandmeter, 2017 ), with the improvement that we have extend the origin version to higher order case.
 *
 * Input:
 *       double[Np x K] fphys the physical field to be limited.
 * 		double[1 x K]  avar the average value of each cell for the physical field
 * 		double[Nface x K]  EToE the topological relation of the studied mesh
 * 		double[Npz x Npz]  V1d the inversed one-dimensional Vandmonde matrix, used to calculate the average value in vertical direction for each line
 *       int[1] Npz number of interpolation points in vertical direction
 * 		int[1] Nph number of interpolation points in horizontal direction
 *       double[Npz x Npz]  OV1d the one-dimensional Vandmonde matrix, the first value of this matrix is used to calculate the average value
 * Output:
 * 		double[Np x K] limfphys the limited physical field.
 */
//CalculateSquareValue( wq, J, Vq, squafphy, sumSqua, Nq, Np)
void CalculateSquareValue( double *w, double *J, double *Vqua, double *sfphy, double *sum, size_t Nq, size_t Np){
    
    const double one = 1;
    const double zero = 0;
    ptrdiff_t one_ptrdiff = 1;
    ptrdiff_t Np_ptrdiff = Np;
    ptrdiff_t Nq_ptrdiff = Nq;
    char *tran = "N";
    //double Jq[Nq], fq[Nq];
	double *Jq = (double*)malloc(sizeof(double)*Nq), \
		*fq = (double*)malloc(sizeof(double)*Nq);
    
    // map the node values fvar to quadrature nodes by
    // \f$ fq = Vq * fvar \f$
    dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, Vqua,
            &Nq_ptrdiff, sfphy, &Np_ptrdiff, &zero, fq, &Nq_ptrdiff);
    dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, Vqua,
            &Nq_ptrdiff, J, &Np_ptrdiff, &zero, Jq, &Nq_ptrdiff);
    
    for (int n = 0; n < Nq; n++) {
        *sum += w[n] * Jq[n] * fq[n];
    }
	free(Jq);
	free(fq);
}

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
    double *VF = mxGetPr(prhs[0]);
    double SR = mxGetScalar(prhs[1]);
    double *FS = mxGetPr(prhs[2]);
    /*This is the inversed one-dimensional Vandmande matrix*/
    double *wq = mxGetPr(prhs[3]);
    double *J = mxGetPr(prhs[4]);
    double *Vq = mxGetPr(prhs[5]);
    /*This is the one-dimensional VandMande matrix*/
    double *fphys = mxGetPr(prhs[6]);
    double *IV = mxGetPr(prhs[7]);
    double *V = mxGetPr(prhs[8]);
    
    size_t Np = mxGetM(prhs[6]);
    size_t K = mxGetN(prhs[6]);
    size_t Nq = mxGetM(prhs[5]); // number of quadrature points
    
    /* allocate output array */
    plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
    double *limfphys = mxGetPr(plhs[0]);
    
    const double one = 1;
    const double zero = 0;
    ptrdiff_t one_ptrdiff = 1;
    ptrdiff_t Np_ptrdiff = Np; ///////Need to be checked again
    ptrdiff_t Nq_ptrdiff = Np;
    char *tran = "N";
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for(int k=0; k<K; k++){
        //Calculate mode coefficient first
        double *fmod = (double*)malloc(sizeof(double)*Np),\
                *trunfphy = (double*)malloc(sizeof(double)*Np),\
                *squafphy = (double*)malloc(sizeof(double)*Np),\
                *trunsquafphy = (double*)malloc(sizeof(double)*Np);
        dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, IV,
                &Nq_ptrdiff, fphys + k * Np, &Np_ptrdiff, &zero, fmod, &Nq_ptrdiff);
        //Calculate the truncted solution
        dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, VF,
                &Nq_ptrdiff, fmod, &Np_ptrdiff, &zero, trunfphy, &Nq_ptrdiff);
        // Calculate the square of the solution
        double sumSqua = 0;
        for(int i=0; i<Np; i++)
            squafphy[i] = pow(*(fphys + k * Np + i), 2);
        CalculateSquareValue( wq, J, Vq, squafphy, &sumSqua, Nq, Np);
        
        // Calculate the square of the solution minus the truncted solution
        double trunctedSumSqua = 0;
        for(int i=0; i<Np; i++)
            trunsquafphy[i] = pow(*(fphys + k * Np + i) - *( trunfphy + i ), 2);
        CalculateSquareValue( wq, J, Vq, trunsquafphy, &trunctedSumSqua, Nq, Np);
        
        //Calculate SOmega
        double SOmega = 0;
        SOmega = trunctedSumSqua/( sumSqua + EPS );
        
        //Calculate SI
        double SI;
        SI = log10( SOmega + EPS );
        // Decide whether to filter the studied cell
        while(SI > SR)
        {

			sumSqua = 0;

			trunctedSumSqua = 0;
            // filter the solution
            for(int i = 0;i<Np;i++)
                fmod[i] = fmod[i]* FS[i];
            //Calculate the solution after filter
            dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, V,
                    &Nq_ptrdiff, fmod, &Np_ptrdiff, &zero, fphys + k * Np, &Nq_ptrdiff);
            
            //Calculate the truncted solution
            dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, VF,
                    &Nq_ptrdiff, fmod, &Np_ptrdiff, &zero, trunfphy, &Nq_ptrdiff);
            // Calculate the square of the solution
            for(int i=0; i<Np; i++)
                squafphy[i] = pow(*(fphys + k * Np + i), 2);
            CalculateSquareValue( wq, J, Vq, squafphy, &sumSqua, Nq, Np);
            
            // Calculate the square of the solution minus the truncted solution
            for(int i=0; i<Np; i++)
                trunsquafphy[i] = pow(*(fphys + k * Np + i) - *( trunfphy + i ), 2);
            CalculateSquareValue( wq, J, Vq, trunsquafphy, &trunctedSumSqua, Nq, Np);
            
            //Calculate SOmega
            SOmega = trunctedSumSqua/( sumSqua + EPS );
            
            //Calculate SI
            SI = log10( SOmega + EPS );
            
        }
        free(fmod);
        free(trunfphy);
        free(squafphy);
        free(trunsquafphy);
    }
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    /*Copy the limited value to the output*/
    for(int i=0; i<Np*K; i++)
        limfphys[i] = fphys[i];
}