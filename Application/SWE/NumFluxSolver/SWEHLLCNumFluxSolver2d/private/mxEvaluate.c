#include <math.h>
#include "mex.h"
#include "../../SWEAbstractNumFluxSolver2d/private/SWENumFlux2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Calculation of the HLL numerical flux
 * @param [in] hmin - threadhold of the water depth
 * @param [in] gra - gravity accelerated
 * @param [in] hL, qnL, qvL - water depth on local element node
 * @param [in] hR, qnR, qvR - water depth on adjacent element node
 * @param [out] Fhn, Fqxn, Fqyn - HLL numerical flux;
 */
void evaluateHLLCFormula(double hmin,
                        double gra,
                        double hL,
                        double qnL,
                        double qvL,
                        double hR,
                        double qnR,
                        double qvR,
                        double* Fhn,
                        double* Fqxn,
                        double* Fqyn) {
  double eps = 1e-6;
  double uL, vL;
  double uR, vR;
  double AR = 0, AL = 0;
  double hstar=0, ustar=0;
  double PQL = 0, PQR = 0;
  double SL, SR;
  double DENOM;
  double SSTAR;
  int SPY = 0;
  double QL[3], QR[3];
  double QSTARL[3], QSTARR[3];
  double POND = 0;
  double FL[3], FR[3];
  double FSTARL[3], FSTARR[3];
  
  evaluateFlowRateByDeptheThreshold(hmin, hR, qnR, qvR, &uR, &vR);
  evaluateFlowRateByDeptheThreshold(hmin, hL, qnL, qvL, &uL, &vL);
  
  if (!(hL < eps && hR < eps))
  {
	  /*CELERITIES*/
	  AL = sqrt(gra*hL), AR = sqrt(gra*hR);
	  
	  /*STAR VARIABLES*/
	  hstar = 0.5*(hL + hR) - 0.25*(uR - uL)*(hL + hR) / (AL + AR);
	  ustar = 0.5*(uL + uR) - (hR - hL)*(AL + AR) / (hL + hR);
	  
	  /*IT WILL DEPEND IF WE ARE IN PRESENCE OF SHOCK OR RAREFACTION WAVE*/
	  if (hstar < hL)
	  {
		  PQL = 1;
	  }else if (hL>eps){
		  PQL = sqrt(0.5*(hstar + hL)*hstar / pow(hL,2));
	  }

	  if (hstar < hR)
	  {
		  PQR = 1;
	  }else if (hR>eps){
		  PQR = sqrt(0.5*(hstar + hR)*hstar / pow(hR,2));
	  }
  }
  /*COMPUTE SL, SR AND SSTAR  (WE CONSIDER DRY CASES)*/
  
  if (hL > eps && hR > eps){
	  SL = uL - AL*PQL;
	  SR = uR + AR*PQR;
  }else if (hL > eps && hR <= eps){
	  SL = uL - AL;
	  SR = uL + 2 * AL;
  }else{
	  SL = uR - 2 * AR;
	  SR = uR + AR;
  }


  DENOM = hR*(uR - SR) - hL*(uL - SL);

  if (abs(DENOM) < eps){
	  SSTAR = ustar;
  }else{
	  SSTAR = (SL*hR*(uR - SR) - SR*hL*(uL - SL)) / DENOM;
  }

  /*COMPUTE QL AND QR*/
  QL[0] = hL;
  QL[1] = hL*uL;
  QL[2] = hL*vL;

  QR[0] = hR;
  QR[1] = hR*uR;
  QR[2] = hR*vR;

  /*COMPUTE QSTARL AND QSTARR*/
  if (abs(SL - SSTAR) > eps){
	  POND = hL*(SL-uL) / (SL-SSTAR);
  }

  QSTARL[0] = POND;
  QSTARL[1] = POND*SSTAR;
  QSTARL[2] = POND*vL;

  POND = 0;

  if (abs(SR - SSTAR) > eps){
	  POND = hR*(SR - uR) / (SR - SSTAR);
  }

  QSTARR[0] = POND;
  QSTARR[1] = POND*SSTAR;
  QSTARR[2] = POND*vR;

  /*Compute FL AND FR*/

  FL[0] = hL*uL;
  FL[1] = hL*pow(uL,2) + 0.5*gra*pow(hL,2);
  FL[2] = hL*uL*vL;

  FR[0] = hR*uR;
  FR[1] = hR*pow(uR,2) + 0.5*gra*pow(hR,2);
  FR[2] = hR*uR*vR;

  /*Compute FSTARL AND FSTARR*/
 
  FSTARL[0] = FL[0] + SL*(QSTARL[0] - QL[0]);
  FSTARL[1] = FL[1] + SL*(QSTARL[1] - QL[1]);
  FSTARL[2] = FL[2] + SL*(QSTARL[2] - QL[2]);

  FSTARR[0] = FR[0] + SR*(QSTARR[0] - QR[0]);
  FSTARR[1] = FR[1] + SR*(QSTARR[1] - QR[1]);
  FSTARR[2] = FR[2] + SR*(QSTARR[2] - QR[2]);

  /*AND FINALLY THE HLLC FLUX (BEFORE ROTATION)*/
  if (0 < SL){
	  *Fhn = FL[0];
	  *Fqxn = FL[1];
	  *Fqyn = FL[2];
	  SPY = 1;
  }else if (0 < SSTAR && 0 > SL){
	  *Fhn = FSTARL[0];
	  *Fqxn = FSTARL[1];
	  *Fqyn = FSTARL[2];
	  SPY = 1;
  }else if (0 > SSTAR && 0 < SR){
	  *Fhn = FSTARR[0];
	  *Fqxn = FSTARR[1];
	  *Fqyn = FSTARR[2];
	  SPY = 1;
  }else{
	  *Fhn = FR[0];
	  *Fqxn = FR[1];
	  *Fqyn = FR[2];
	  SPY = 1;
  }
  if (SPY == 0){
	  printf("Error occured when calculating the HLLC flux! check please!");
  }
  return;
}

#define NRHS 6
#define NLHS 1

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  /* check input & output */
  if (nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NLHS);
  }
  double hcrit = mxGetScalar(prhs[0]);
  double gra = mxGetScalar(prhs[1]);
  double* nx = mxGetPr(prhs[2]);
  double* ny = mxGetPr(prhs[3]);
  double* fm = mxGetPr(prhs[4]);
  double* fp = mxGetPr(prhs[5]);

  const mwSize* dims = mxGetDimensions(prhs[4]);
  const size_t TNfp = dims[0];
  const size_t K = dims[1];

  const size_t NdimOut = 3;
  const mwSize dimOut[3] = {TNfp, K, 3};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double* Fh = mxGetPr(plhs[0]);
  double* Fqx = Fh + TNfp * K;
  double* Fqy = Fh + 2 * TNfp * K;

  double* hL = fm;
  double* hum = fm + K * TNfp;
  double* hvm = fm + 2 * K * TNfp;

  double* hR = fp;
  double* hup = fp + K * TNfp;
  double* hvp = fp + 2 * K * TNfp;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;
      double fm[3] = {hL[sk], hum[sk], hvm[sk]};
      double fp[3] = {hR[sk], hup[sk], hvp[sk]};
      const double nx_ = nx[sk];
      const double ny_ = ny[sk];
	  double Fhns, Fqns, Fqvs;
      const double qnL = fm[1] * nx_ + fm[2] * ny_;
      const double qvL = -fm[1] * ny_ + fm[2] * nx_;
      const double qnR = fp[1] * nx_ + fp[2] * ny_;
      const double qvR = -fp[1] * ny_ + fp[2] * nx_;	  
	  evaluateHLLCFormula(hcrit, gra, fm[0],qnL, qvL, fp[0],qnR ,qvR ,&Fhns, &Fqns, &Fqvs);
	  Fh[sk] = Fhns;
      Fqx[sk] = (Fqns * nx_ - Fqvs * ny_);
      Fqy[sk] = (Fqns * ny_ + Fqvs * nx_);

#if 0
      mexPrintf("k=%d, n=%d, hL=%f, hR=%f, hum=%f, hup=%f, hvm=%f, hvp=%f\n", k,
                n, fm[0], fp[0], qnL, qnR, qvL, qvR);
      mexPrintf("k=%d, n=%d, Fh=%f, Fqxn=%f, Fqyn=%f\n", k, n, Fhns, Fqns,
                Fqvs);
      mexPrintf("k=%d, n=%d, Fh=%f, Fqx=%f, Fqy=%f\n", k, n, Fh[sk], Fqx[sk],
                Fqy[sk]);
#endif
    }
  }

  return;
}
