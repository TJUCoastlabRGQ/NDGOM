#include "NdgSWE3D.h"
#include "NdgSWE.h"

void EvaluateVerticalFaceSurfFlux(double *dest, double *fm, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne){
	double u, v, theta;
	double *hu = fm, *hv = fm + Nfp*Ne, *h = fm + 2 * Nfp*Ne;
	for (int i = 0; i < Nfp; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hu + i, &u);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hv + i, &v);
		dest[i] = (hu[i] * u + 1 / 2 * *gra*h[i] * h[i])*nx[i] + hu[i] * v*ny[i];
		dest[i + Ne*Nfp] = hu[i] * v*nx[i] + (hv[i] * v + 1 / 2 * *gra*h[i] * h[i])*ny[i];
		for (int field = 3; field < Nvar+1; field++){
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, fm + field*Ne*Nfp + i, &theta);
			dest[i + (field-1)*Ne*Nfp] = hu[i] * theta * nx[i] + hv[i] * theta * ny[i];
		}
	}
}

void EvaluateHorizontalFaceSurfFlux(double *flux, double *fm, double *nz, double Hcrit, int Nfp, int Nvar, int Ne){
	double *hum = fm, *hvm = fm + Nfp*Ne, *omega = fm + 2 * Nfp*Ne, *H = fm + 3 * Nfp*Ne;
	double u, v, variable, theta;
	for (int i = 0; i < Nfp; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, H + i, hum + i, &u);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, H + i, hvm + i, &v);
		*(flux + i) = u*omega[i]*nz[i];
		*(flux + Ne*Nfp + i) = v*omega[i]*nz[i];
		for (int n = 2; n < Nvar; n++){
			/*Here 2*Ne*Nfp stands for the memory occupied by h and omega*/
			variable = *(fm + 2 * Ne*Nfp + n*Ne*Nfp + i);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, H + i, &variable, &theta);
			*(flux + n*Ne*Nfp + i) = theta*omega[i]*nz[i];
		}
	}
}

/*HLLC numerical flux is adopted in vertical face*/
void EvaluateVerticalFaceNumFlux(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne){
	double HR, HL, UR, UL, VR, VL;
	double SL, SSTAR, SR;
	double *hum = fm, *hvm = fm + Nfp*Ne, *hm = fm + 2 * Nfp*Ne;
	double *hup = fp, *hvp = fp + Nfp*Ne, *hp = fp + 2 * Nfp*Ne;

	double *QL = malloc((Nvar + 1)*sizeof(double)), *QR = malloc((Nvar + 1)*sizeof(double));
	double *QSTARL = malloc((Nvar + 1)*sizeof(double)), *QSTARR = malloc((Nvar + 1)*sizeof(double));
	double *FL = malloc((Nvar + 1)*sizeof(double)), *FR = malloc((Nvar + 1)*sizeof(double));
	double *FSTARL = malloc((Nvar + 1)*sizeof(double)), *FSTARR = malloc((Nvar + 1)*sizeof(double));
	double *FLX = malloc((Nvar + 1)*sizeof(double));
	double *VARIABLEL = malloc((Nvar + 1)*sizeof(double)), *VARIABLER = malloc((Nvar + 1)*sizeof(double));
	for (int i = 0; i < Nfp; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, HUt, H theta*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);
		for (int n = 3; n < Nvar + 1; n++)
		{
			QL[n] = *(fm + n*Nfp*Ne + i);
			QR[n] = *(fp + n*Nfp*Ne + i);
		}

		double AR = 0, AL = 0;   //checked
		double HSTAR = 0, USTAR = 0;
		double PQL = 0, PQR = 0;
		double DENOM;
		int SPY = 0;
		double POND = 0;
		/*Compute the original variable, u, v and theta, in normal direction*/
		/*Water depth h comes first*/
		VARIABLEL[0] = QL[0];
		VARIABLER[0] = QR[0];
		for (int n = 1; n < Nvar + 1; n++){
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLEL, QL + n, VARIABLEL+n);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLER, QR + n, VARIABLER+n);
		}
		/*Assign water depth, u and v in both sides to HL(R), UL(R), VL(R)*/
		HL = VARIABLEL[0], HR = VARIABLER[0];
		UL = VARIABLEL[1], UR = VARIABLER[1];
		VL = VARIABLEL[2], VR = VARIABLER[2];

		if (!(HL < Hcrit && HR < Hcrit))
		{
			/*CELERITIES*/
			AL = sqrt(*gra*HL), AR = sqrt(*gra*HR);
			/*STAR VARIABLES*/
			HSTAR = 0.5*(HL + HR) - 0.25*(UR - UL)*(HL + HR) / (AL + AR);
			USTAR = 0.5*(UL + UR) - (HR - HL)*(AL + AR) / (HL + HR);
			/*IT WILL DEPEND IF WE ARE IN PRESENCE OF SHOCK OR RAREFACTION WAVE*/
			if (HSTAR < HL)
			{
				PQL = 1;
			}
			else if (HL>Hcrit){
				PQL = sqrt(0.5*(HSTAR + HL)*HSTAR / pow(HL, 2));
			}

			if (HSTAR < HR)
			{
				PQR = 1;
			}
			else if (HR>Hcrit){
				PQR = sqrt(0.5*(HSTAR + HR)*HSTAR / pow(HR, 2));
			}
		}

		if (HL > Hcrit && HR > Hcrit){
			SL = UL - AL*PQL;
			SR = UR + AR*PQR;
		}
		else if (HL > Hcrit && HR <= Hcrit){
			SL = UL - AL;
			SR = UL + 2 * AL;
		}
		else{
			SL = UR - 2 * AR;
			SR = UR + AR;    
		}

		DENOM = HR*(UR - SR) - HL*(UL - SL);

		if (abs(DENOM) < Hcrit){
			SSTAR = USTAR;
		}
		else{
			SSTAR = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / DENOM;
		}
		/*COMPUTE QSTARL AND QSTARR*/
		if (abs(SL - SSTAR) > EPS){
			/*POND is initialized as zero*/
			POND = HL*(SL - UL) / (SL - SSTAR);
		}

		QSTARL[0] = POND;
		QSTARL[1] = POND*SSTAR;
		QSTARL[2] = POND*VL;
		for (int n = 3; n < Nvar+1; n++)
		{
			QSTARL[n] = POND*VARIABLEL[n];
		}

		POND = 0;
		if (abs(SR - SSTAR) > EPS){
			/*POND is initialized as zero*/
			POND = HR*(SR - UR) / (SR - SSTAR);
		}

		QSTARR[0] = POND;
		QSTARR[1] = POND*SSTAR;
		QSTARR[2] = POND*VR;
		for (int n = 3; n < Nvar+1; n++)
		{
			QSTARR[n] = POND*VARIABLER[n];
		}

		/*Compute FL AND FR*/
		FL[0] = HL*UL;
		FL[1] = HL*pow(UL, 2) + 0.5*(*gra)*pow(HL, 2);
		FL[2] = FL[0]*VL;
		for (int n = 3; n < Nvar + 1; n++){
			FL[n] = FL[0] * VARIABLEL[n];
		}

		FR[0] = HR*UR;
		FR[1] = HR*pow(UR, 2) + 0.5*(*gra)*pow(HR, 2);
		FR[2] = HR*UR*VR;
		for (int n = 3; n < Nvar + 1; n++){
			FR[n] = FR[0] * VARIABLER[n];
		}
		/*Compute FSTARL AND FSTARR*/
		FSTARL[0] = FL[0] + SL*(QSTARL[0] - QL[0]);
		FSTARL[1] = FL[1] + SL*(QSTARL[1] - QL[1]);
		FSTARL[2] = FL[2] + SL*(QSTARL[2] - QL[2]);
		for (int n = 3; n < Nvar + 1; n++){
			FSTARL[n] = FL[n] + SL*(QSTARL[n] - QL[n]);
		}

		FSTARR[0] = FR[0] + SR*(QSTARR[0] - QR[0]);
		FSTARR[1] = FR[1] + SR*(QSTARR[1] - QR[1]);
		FSTARR[2] = FR[2] + SR*(QSTARR[2] - QR[2]);
		for (int n = 3; n < Nvar + 1; n++){
			FSTARR[n] = FR[n] + SR*(QSTARR[n] - QR[n]);
		}
		/*AND FINALLY THE HLLC FLUX (BEFORE ROTATION)*/
		if (0 < SL){
			//*Fhn = FL[0];
			dest[i] = (FL[1] * nx[i] - FL[2] * ny[i]);
			dest[i + Nfp*Ne] = (FL[1] * ny[i] + FL[2] * nx[i]);
			for (int n = 3; n < Nvar+1; n++){
				dest[i + (n-1) * Nfp*Ne] = FL[n];
			}
			SPY = 1;
		}
		else if (0 < SSTAR && 0 > SL){
			dest[i] = (FSTARL[1] * nx[i] - FSTARL[2] * ny[i]);
			dest[i + Nfp*Ne] = (FSTARL[1] * ny[i] + FSTARL[2] * nx[i]);
			for (int n = 3; n < Nvar+1; n++){
				dest[i + (n - 1) * Nfp*Ne] = FSTARL[n];
			}
			SPY = 1;
		}
		else if (0 > SSTAR && 0 < SR){
			dest[i] = (FSTARR[1] * nx[i] - FSTARR[2] * ny[i]);
			dest[i + Nfp*Ne] = (FSTARR[1] * ny[i] + FSTARR[2] * nx[i]);
			for (int n = 3; n < Nvar + 1; n++){
				dest[i + (n - 1) * Nfp*Ne] = FSTARR[n];
			}
			SPY = 1;
		}
		else{
			dest[i] = (FR[1] * nx[i] - FR[2] * ny[i]);
			dest[i + Nfp*Ne] = (FR[1] * ny[i] + FR[2] * nx[i]);
			for (int n = 3; n < Nvar + 1; n++){
				dest[i + (n - 1) * Nfp*Ne] = FR[n];
			}
			SPY = 1;
		}
		if (SPY == 0){
			printf("Error occured when calculating the HLLC flux! check please!");
		}
	}
	free(QL), free(QR);
	free(QSTARL), free(QSTARR);
	free(FL), free(FR);
	free(FSTARL), free(FSTARR);
	free(FLX);
	free(VARIABLEL), free(VARIABLER);
}

void EvaluateHorizontalFaceNumFlux(double *FluxS, double *fm, double *fp, double *nz, double Hcrit, int Nfp, int Nvar, int Ne){
	double *hum = fm, *hvm = fm + Nfp*Ne, *omegam = fm + 2 * Nfp*Ne, *Hm = fm + 3 * Nfp*Ne;
	double *hup = fp, *hvp = fp + Nfp*Ne, *omegap = fp + 2 * Nfp*Ne, *Hp = fp + 3 * Nfp*Ne;
	double um, vm, up, vp, variablem, variablep, thetam, thetap;
	for (int i = 0; i < Nfp; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, hum + i, &um);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, hvm + i, &vm);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, hup + i, &up);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, hvp + i, &vp);
		*(FluxS + i) = (um*omegam[i] + up*omegap[i])/2 * nz[i];
		*(FluxS + Ne*Nfp + i) = (vm*omegam[i] + vp*omegap[i]) / 2 * nz[i];
		for (int n = 2; n < Nvar; n++){
			/*Here 2*Ne*Nfp stands for the memory occupied by h and omega*/
			variablem = *(fm + 2 * Ne*Nfp + n*Ne*Nfp + i);
			variablep = *(fp + 2 * Ne*Nfp + n*Ne*Nfp + i);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, &variablem, &thetam);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, &variablep, &thetap);
			*(FluxS + n*Ne*Nfp + i) = (thetam*omegam[i] + thetap*omegap[i]) / 2 * nz[i];
		}
	}

}

void EvaluatePrebalanceVolumeTerm(double *Edest, double *Gdest, double *Hdest, double *fphys, \
	double *varIndex, int Nvar, double *gra, int Np, int K, double Hcrit)
{
	double *hu = fphys, *hv = fphys + Np*K, *omega = fphys + 2 * Np*K;
	double *h = fphys + 3 * Np*K, *z = fphys + 5 * Np*K;
	double u, v, variable, theta;
	for (int i = 0; i < Np; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hu + i, &u);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hv + i, &v);
		/*NOTE THAT: THE WET AND DRY NEED TO BE TACKLED HERE*/
		Edest[i] = hu[i] * u + 1 / 2 * (*gra)*(pow(h[i], 2) - pow(z[i], 2));
		Gdest[i] = hu[i] * v;
		Hdest[i] = u*omega[i];
		Edest[i + Np*K] = hu[i] * v;
		Gdest[i + Np*K] = hv[i] * v + 1 / 2 * (*gra)*(pow(h[i], 2) - pow(z[i], 2));
		Hdest[i + Np*K] = v*omega[i];
		/*If temperature, salt, sediment or other passive transport material is included, the following part is used to calculate
		the volume term corresponding to these terms*/
		for (int n = 2; n < Nvar; n++){
			variable = *(fphys + ((int)varIndex[n] - 1)*Np*K + i);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, &variable, &theta);
			Edest[i + n*Np*K] = hu[i] * theta;
			Gdest[i + n*Np*K] = hv[i] * theta;
			Hdest[i + n*Np*K] = omega[i] * theta;
		}
	}
}

/** Evaluate flux term in surface integration */
/*
void EvaluateFluxTerm2d( double hmin, ///< water threshold
	double *gra,  ///< gravity acceleration
	double *h,    ///< water depth
	double *hu,   ///< flux variable
	double *hv,   ///< flux variable
	double *E        ///< surface integral flux term
	) {
	double u, v;
	if (h > hmin) {
		u = hu / h;
		v = hv / h;
	}
	else {
		u = 0.0;
		v = 0.0;
	}
	double huv = h * u * v;
    double h2 = h * h;
	E[0] = hu;
	E[1] = h * u * u + 0.5 * gra * h2;
	return;
}
*/