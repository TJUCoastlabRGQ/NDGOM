#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE.h"
#include <omp.h>

void EvaluateNumFlux_HLLC_LAI(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne);


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	double hcrit = mxGetScalar(prhs[0]);
	double gra = mxGetScalar(prhs[1]);
	double* nx = mxGetPr(prhs[2]);
	double* ny = mxGetPr(prhs[3]);
	double* fm = mxGetPr(prhs[4]);
	double* fp = mxGetPr(prhs[5]);
	int Nvar = (int)mxGetScalar(prhs[6]);

	const mwSize* dims = mxGetDimensions(prhs[4]);
	const size_t TNfp = dims[0];
	const size_t Ne = dims[1];

	const size_t NdimOut = 3;
	const mwSize dimOut[3] = { TNfp, Ne, 3 };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

	double* Fh = mxGetPr(plhs[0]);
	double* Fqx = Fh + TNfp * Ne;
	double* Fqy = Fh + 2 * TNfp * Ne;

	double* hL = fm;
	double* hum = fm + Ne * TNfp;
	double* hvm = fm + 2 * Ne * TNfp;

	double* hR = fp;
	double* hup = fp + Ne * TNfp;
	double* hvp = fp + 2 * Ne * TNfp;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < Ne; face++){

		EvaluateNumFlux_HLLC_LAI(Fh + face*TNfp, fm + face*TNfp, fp + face*TNfp, \
			nx + face*TNfp, ny + face*TNfp, &gra, hcrit, TNfp, Nvar, Ne);
	}

	return;
}

/*HLLC numerical flux is adopted in vertical face, */
void EvaluateNumFlux_HLLC_LAI(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne){
	/*The HLLC numerical flux presented in LAI(2012, Modeling one- and two-dimensional shallow water flows with discontinuous Galerkin method) is adopted*/
	double HR, HL, UR, UL, VR, VL;
	double SL, SM, SR;
	double USTAR, HSTAR;
	double *hm = fm, *hum = fm + Nfp*Ne, *hvm = fm + 2 * Nfp*Ne;
	double *hp = fp, *hup = fp + Nfp*Ne, *hvp = fp + 2 * Nfp*Ne;

	double *QL = malloc(Nvar*sizeof(double)), *QR = malloc(Nvar*sizeof(double));
	double *FL = malloc(Nvar*sizeof(double)), *FR = malloc(Nvar*sizeof(double));
	double *FSTARL = malloc(Nvar*sizeof(double)), *FSTARR = malloc(Nvar*sizeof(double));
	double *VARIABLEL = malloc(Nvar*sizeof(double)), *VARIABLER = malloc(Nvar*sizeof(double));
	double *QSTARL = malloc(Nvar*sizeof(double)), *QSTARR = malloc(Nvar*sizeof(double));
	for (int i = 0; i < Nfp; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, HUt, H theta*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);
		for (int n = 3; n < Nvar; n++)
		{
			QL[n] = *(fm + n*Nfp*Ne + i);
			QR[n] = *(fp + n*Nfp*Ne + i);
		}

		int SPY = 0;
		/*Compute the original variable, u, v and theta, in normal direction*/
		/*Water depth h comes first*/
		VARIABLEL[0] = QL[0];
		VARIABLER[0] = QR[0];
		for (int n = 1; n < Nvar; n++){
			/*In this part, we calculate the original variable. If the water depth is smaller than the critical value,
			we directly set the value of the original variable to be zero, and the original value is further used to calculate
			the flux term in left and right side. The question here is whether the water depth should be set to zero if it is
			smaller than the critical value, since it is used when calcualte the flux term. In the current version, we leave them
			unchanged*/
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLEL, QL + n, VARIABLEL + n);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLER, QR + n, VARIABLER + n);
		}
		/*Assign water depth, u and v in both sides to HL(R), UL(R), VL(R)*/
		HL = VARIABLEL[0], HR = VARIABLER[0];
		UL = VARIABLEL[1], UR = VARIABLER[1];
		VL = VARIABLEL[2], VR = VARIABLER[2];
		/*Here, we use the critical value to categorize different condition, I wonder whether this value can be
		directly set to zero, if this value is set to zero, then they calculation in the further part can be looked
		to be very smooth*/
		if (!(HL < Hcrit && HR < Hcrit))
		{

			/*Compute FL AND FR*/
			FL[0] = HL*UL;
			FL[1] = HL*pow(UL, 2) + 0.5*(*gra)*pow(HL, 2);
			FL[2] = HL*UL*VL;
			for (int n = 3; n < Nvar; n++){
				FL[n] = FL[0] * VARIABLEL[n];
			}

			FR[0] = HR*UR;
			FR[1] = HR*pow(UR, 2) + 0.5*(*gra)*pow(HR, 2);
			FR[2] = HR*UR*VR;
			for (int n = 3; n < Nvar; n++){
				FR[n] = FR[0] * VARIABLER[n];
			}

			if ((HL>Hcrit) & (HR > Hcrit)){
				USTAR = 0.5 * (UL + UR) + sqrt(*gra*HL) - sqrt(*gra*HR);
				HSTAR = 1.0 / (*gra)*pow(0.5*(sqrt(*gra*HL) + sqrt(*gra*HR)) + 0.25*(UL - UR), 2);
				SL = min(UL - sqrt(*gra * HL), USTAR - sqrt(*gra*HSTAR));
				SR = max(UR + sqrt(*gra * HR), USTAR + sqrt(*gra*HSTAR));
				SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
				/*Compute QSTARL AND QSTARR*/
				QSTARL[0] = HL*((SL - UL) / (SL - SM)) * 1;
				QSTARL[1] = HL*((SL - UL) / (SL - SM)) * SM;
				QSTARR[0] = HR*((SR - UR) / (SR - SM)) * 1;
				QSTARR[1] = HR*((SR - UR) / (SR - SM)) * SM;
				for (int i = 2; i<Nvar; i++){
					QSTARL[i] = HL*((SL - UL) / (SL - SM)) * VARIABLEL[i];
					QSTARR[i] = HR*((SR - UR) / (SR - SM)) * VARIABLER[i];
				}
				/*Compute FSTARL and FSTARR*/
				for (int n = 0; n < Nvar; n++){
					FSTARL[n] = FL[n] + SL * (QSTARL[n] - QL[n]);
					FSTARR[n] = FR[n] + SR * (QSTARR[n] - QR[n]);
				}
				/*AND FINALLY THE HLLC FLUX (BEFORE ROTATION)*/
				if (SL >= 0){
					dest[i] = FL[0];
					dest[i + Nfp*Ne] = (FL[1] * nx[i] - FL[2] * ny[i]);
					dest[i + 2 * Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * Nfp*Ne] = FL[n];
					}
					SPY = 1;
				}
				else if (SM >= 0 && SL < 0){
					dest[i] = FSTARL[0];
					dest[i + Nfp*Ne] = (FSTARL[1] * nx[i] - FSTARL[2] * ny[i]);
					dest[i + 2 * Nfp*Ne] = (FSTARL[1] * ny[i] + FSTARL[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * Nfp*Ne] = FSTARL[n];
					}
					SPY = 1;
				}
				else if (SM < 0 && SR > 0){
					dest[i] = FSTARR[0];
					dest[i + Nfp*Ne] = (FSTARR[1] * nx[i] - FSTARR[2] * ny[i]);
					dest[i + 2 * Nfp*Ne] = (FSTARR[1] * ny[i] + FSTARR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * Nfp*Ne] = FSTARR[n];
					}
					SPY = 1;
				}
				else if (0 >= SR){
					dest[i] = FR[0];
					dest[i + Nfp*Ne] = (FR[1] * nx[i] - FR[2] * ny[i]);
					dest[i + 2 * Nfp*Ne] = (FR[1] * ny[i] + FR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * Nfp*Ne] = FR[n];
					}
					SPY = 1;
				}

			}
			else if (HL>Hcrit && HR <= Hcrit){
				SL = UL - sqrt(*gra*HL);
				SR = UL + 2 * sqrt(*gra*HL);
				//是否将SM的计算移到这里来
				//				SM = SR;
				if (SL >= 0){
					dest[i] = FL[0];
					dest[i + Nfp*Ne] = (FL[1] * nx[i] - FL[2] * ny[i]);
					dest[i + 2 * Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * Nfp*Ne] = FL[n];
					}
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					double Fqxn = (SR * FL[1] - SL * FR[1] + SL * SR * (HR*UR - HL*UL)) / (SR - SL);
					double Fqyn = (SR * FL[2] - SL * FR[2] + SL * SR * (HR*VR - HL*VL)) / (SR - SL);
					dest[i] = (SR * FL[0] - SL * FR[0] + SL * SR * (HR - HL)) / (SR - SL);
					dest[i + Nfp*Ne] = (Fqxn * nx[i] - Fqyn * ny[i]);
					dest[i + 2 * Nfp*Ne] = (Fqxn * ny[i] + Fqyn * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n*Nfp*Ne] = (SR * FL[n] - SL * FR[n] + SL * SR * (QR[n] - QL[n])) / (SR - SL);
					}
					SPY = 1;
				}
				else if (SR <= 0){
					dest[i] = FR[0];
					dest[i + Nfp*Ne] = (FR[1] * nx[i] - FR[2] * ny[i]);
					dest[i + 2 * Nfp*Ne] = (FR[1] * ny[i] + FR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * Nfp*Ne] = FR[n];
					}
					SPY = 1;
				}
			}
			else if (HL <= Hcrit && HR > Hcrit){
				SL = UR - 2 * sqrt(*gra*HR);
				SR = UR + sqrt(*gra*HR);
				//				SM = SL;
				if (SL >= 0){
					dest[i] = FL[0];
					dest[i + Nfp*Ne] = (FL[1] * nx[i] - FL[2] * ny[i]);
					dest[i + 2 * Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * Nfp*Ne] = FL[n];
					}
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					double Fqxn = (SR * FL[1] - SL * FR[1] + SL * SR * (HR*UR - HL*UL)) / (SR - SL);
					double Fqyn = (SR * FL[2] - SL * FR[2] + SL * SR * (HR*VR - HL*VL)) / (SR - SL);
					dest[i] = (SR * FL[0] - SL * FR[0] + SL * SR * (HR - HL)) / (SR - SL);
					dest[i + Nfp*Ne] = (Fqxn * nx[i] - Fqyn * ny[i]);
					dest[i + 2 * Nfp*Ne] = (Fqxn * ny[i] + Fqyn * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n*Nfp*Ne] = (SR * FL[n] - SL * FR[n] + SL * SR * (QR[n] - QL[n])) / (SR - SL);
					}
					SPY = 1;
				}
				else if (SR <= 0){
					dest[i] = FR[0];
					dest[i + Nfp*Ne] = (FR[1] * nx[i] - FR[2] * ny[i]);
					dest[i + 2 * Nfp*Ne] = (FR[1] * ny[i] + FR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * Nfp*Ne] = FR[n];
					}
					SPY = 1;
				}
			}
			//SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
			if (SPY == 0){
				printf("Error occured when calculating the HLLC flux! check please!");
			}
		}
	}
	free(QL), free(QR);
	free(QSTARL), free(QSTARR);
	free(FL), free(FR);
	free(FSTARL), free(FSTARR);
	free(VARIABLEL), free(VARIABLER);
}
