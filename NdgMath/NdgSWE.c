#include "NdgSWE.h"

/*This function is used to impose the boundary condition for the pure hydrulic problem*/
void ImposeBoundaryCondition(double gra, NdgEdgeType type, double *nx, double *ny, double *fm, double *fp, \
	double *zM, double *zP, double *fext, int Nfp, int Nvar, int Ne){
	/*All other fields except H, Hu, Hv are set to be the inner value at the boundaries unless they are prescribed out of the program */
	// assign the local node values
	double *huM = fm, *hvM = fm + Nfp*Ne, *hM = fm + 2 * Nfp*Ne;
	double *huP = fp, *hvP = fp + Nfp*Ne, *hP = fp + 2 * Nfp*Ne;
	double *huE = fext, *hvE = fext + Nfp*Ne, *hE = fext + 2 * Nfp*Ne;
	for (int i = 0; i < Nfp; i++){
		zP[i] = zM[i];
	}

	// get next node values
	if (type == NdgEdgeInner) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hM[i];
			huP[i] = huM[i];
			hvP[i] = hvM[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeZeroGrad) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hM[i];
			huP[i] = huM[i];
			hvP[i] = hvM[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeClamped) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hE[i];
			huP[i] = huE[i];
			hvP[i] = hvE[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fext[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeClampedDepth) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hE[i];
			huP[i] = huM[i];
			hvP[i] = hvM[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeClampedVel) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hM[i];
			huP[i] = huE[i];
			hvP[i] = hvE[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeSlipWall) {
		for (int i = 0; i < Nfp; i++){
			const double qxM = huM[i];
			const double qyM = hvM[i];
			double qnM = qxM * nx[i] + qyM * ny[i];   // outward normal flux
			double qvM = -qxM * ny[i] + qyM * nx[i];  // outward tangential flux
			// adjacent value
			hP[i] = hM[i];
			huP[i] = (-qnM) * nx[i] - qvM * ny[i];
			hvP[i] = (-qnM) * ny[i] + qvM * nx[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeNonSlipWall) {
		for (int i = 0; i < Nfp; i++){
			hP[i] = hM[i];
			huP[i] = -huM[i];
			hvP[i] = -hvM[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeFlather) {
		for (int i = 0; i < Nfp; i++){
			const double uE = huE[i] / hE[i];
			const double vE = hvE[i] / hE[i];
			const double unE = uE * nx[i] + vE * ny[i];   // outward normal flux
			const double uvE = -uE * ny[i] + vE * nx[i];  // tangential flux
			const double un = unE - sqrt(-*gra / zP[i]) * (hE[i] - hM[i]);
			const double uv = uvE;
			hP[i] = hM[i];
			huP[i] = (un * nx[i] - uv * ny[i]) * hM[i];
			hvP[i] = (un * ny[i] + uv * nx[i]) * hM[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeNonLinearFlather) {
		for (int i = 0; i < Nfp; i++){
			const double uE = huE[i] / hE[i];
			const double vE = hvE[i] / hE[i];
			const double unE = uE * nx[i] + vE * ny[i];  // outward normal flux
			const double un = huM[i] / hM[i] * nx[i] + hvM[i] / hM[i] * ny[i];
			const double temp = 0.5 * (un - unE) + sqrt(*gra * hE[i]);
			hP[i] = temp * temp / *gra;
			huP[i] = huE[i];
			hvP[i] = hvE[i];
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else if (type == NdgEdgeNonLinearFlatherFlow) {
		for (int i = 0; i < Nfp; i++){
			const double uE = huE[i] / hE[i];
			const double vE = hvE[i] / hE[i];
			const double unE = uE * nx[i] + vE * ny[i];   // outward normal flux
			const double uvE = -uE * ny[i] + vE * nx[i];  // tangential flux
			const double un = unE - 2 * sqrt(*gra * hE[i]) + 2 * sqrt(*gra * hM[i]);
			const double uv = uvE;
			hP[i] = hM[i];
			huP[i] = (un * nx[i] - uv * ny[i]) * hM[i];
			hvP[i] = (un * ny[i] + uv * nx[i]) * hM[i];
			//   } else if (type == NdgEdgeNonReflectingFlux) {
			//     const double unM = surf->huM / surf->hM * nx + surf->hvM / surf->hM *
			//     ny; const double RLP = unM + 2 * sqrt(gra * surf->hM); const double
			//     hs = fext->h[idM]; surf->hP = hs; const double RRM = RLP - 4 *
			//     sqrt(gra * hs); const double un = 0.5 * (RRM + RLP); const double uvM
			//     = -surf->huM / surf->hM * ny + surf->hvM / surf->hM * nx; surf->huP =
			//     (un * nx - uvM * ny) * surf->hM; surf->hvP = (un * ny + uvM * nx) *
			//     surf->hM;
			for (n = 3; n < Nvar; n++){
				fp[n*Nfp*Ne + i] = fm[n*Nfp*Ne + i];
			}
		}
	}
	else {
		printf("Matlab:%s:Unknown boundary type: %d\n", __FILE__, type);
	}
	return;
}

void EvaluateHydroStaticReconstructValue(double hmin, double *fm, double *fp, double *zM, double *zP, int Nfp, int Nvar, int Ne)
{
	/*Passive transport substances, such as temperature and salt are not reconstructed*/
	double um, vm, up, vp, etaM, etaP, zstar;
	double *huM = fm, *huP = fp;
	double *hvM = fm + Nfp*Ne, *hvP = fp + Nfp*Ne;
	double *hM = fm + 2 * Nfp*Ne, *hP = fp + 2 * Nfp*Ne;
	for (int i = 0; i < Nfp; i++){
		zstar = max(zM[i], zP[i]);
		EvaluatePhysicalVariableByDepthThreshold(hmin, hM + i, huM + i, &um);
		EvaluatePhysicalVariableByDepthThreshold(hmin, hM + i, hvM + i, &vm);
		EvaluatePhysicalVariableByDepthThreshold(hmin, hP + i, huP + i, &up);
		EvaluatePhysicalVariableByDepthThreshold(hmin, hP + i, hvP + i, &vp);
		etaM = hM[i] + zM[i];
		etaP = hP[i] + zP[i];
		zstar = min(etaM, zstar);
		huM[i] = hM[i] * um;
		hvM[i] = hM[i] * vm;
		huP[i] = hP[i] * up;
		hvP[i] = hP[i] * vp;
		for (n = 3; n < Nvar; n++){
			double variable;
			/*Passive transport substances, such as temperature and salt are reconstructed following the same strategy for hu and hv*/
			EvaluatePhysicalVariableByDepthThreshold(hmin, hM + i, fm + n*Nfp*Ne + i, &variable);
			fm[n*Nfp*Ne + i] = hM[i] * variable;
			EvaluatePhysicalVariableByDepthThreshold(hmin, hP + i, fp + n*Nfp*Ne + i, &variable);
			fp[n*Nfp*Ne + i] = hP[i] * variable;
		}
		hM[i] = etaM - zstar;
		hP[i] = max(0, etaP - zstar) - max(0, zP[i] - zstar);
		zM[i] = zstar;
		zP[i] = zstar;
	}
}
/*
void EvaluateFlowRateByDeptheThreshold(double hmin, double *h, double *hu, double *hv, double *um, double *vm)
{
	if (*h > hmin) {
		//     const double sqrt2 = 1.414213562373095;
		//     double h4 = pow(h, 4);
		//     *u = sqrt2 * h * hu / sqrt( h4 + max( hcrit, h4 ) );
		//     *v = sqrt2 * h * hv / sqrt( h4 + max( hcrit, h4 ) );
		*um = *hu / *h;
		*vm = *hv / *h;
	}
	else {
		*um = 0.0;
		*vm = 0.0;
	}
}
*/

/*This function is used to calcualte the numerical flux in the primitive continuity equation(PCE) */
void GetPCENumericalFluxTerm(double *dest, double *fm, double *fp, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Ne){
	/*HU ranks first, HV second and H third*/
	double HR, HL, UR, UL, VR, VL;
	double SL, SSTAR, SR;
	double *hum = fm, hvm = fm + Nfp*Ne, hm = fm + 2 * Nfp*Ne;
	double *hup = fp, hvp = fp + Nfp*Ne, hp = fp + 2 * Nfp*Ne;

	double *QL = malloc(3 * sizeof(double)), *QR = malloc(3 * sizeof(double));
	double *VARIABLEL = malloc(3 * sizeof(double)), *VARIABLER = malloc(3 * sizeof(double));
	double QSTARL, QSTARR;
	double FL, FR;
	double FSTARL, FSTARR;
	double FLX;
	for (int i = 0; i < Nfp; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, and HUt*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);

		double AR = 0, AL = 0;   //checked
		double HSTAR = 0, USTAR = 0;
		double PQL = 0, PQR = 0;
		double DENOM;
		double SSTAR;
		int SPY = 0;
		double POND = 0;
		/*Compute the original variable, u, v and theta, in normal direction*/
		/*Water depth h comes first*/
		VARIABLEL[0] = QL[0];
		VARIABLER[0] = QR[0];
		for (int n = 1; n < 3; n++){
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLEL, QL + n, VARIABLEL + n);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, VARIABLER, QR + n, VARIABLER + n);
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
			else if (HR>Hcirt){
				PQR = sqrt(0.5*(HSTAR + HR)*HSTAR / pow(HR, 2));
			}
		}

		if (HL > Hcirt && HR > Hcrit){
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

		QSTARL = POND;

		POND = 0;
		if (abs(SR - SSTAR) > EPS){
			/*POND is initialized as zero*/
			POND = hR*(SR - uR) / (SR - SSTAR);
		}

		QSTARR = POND;

		/*Compute FL AND FR*/
		FL = HL*UL;

		FR = HR*UR;
		/*Compute FSTARL AND FSTARR*/
		FSTARL = FL + SL*(QSTARL - HL);

		FSTARR = FR + SR*(QSTARR - HR);

		/*AND FINALLY THE HLLC FLUX (BEFORE ROTATION)*/
		if (0 < SL){
			dest[i] = FL;
			SPY = 1;
		}
		else if (0 < SSTAR && 0 > SL){
			dest[i] = FSTARL;
			SPY = 1;
		}
		else if (0 > SSTAR && 0 < SR){
			dest[i] = FSTARR;
			SPY = 1;
		}
		else{
			dest[i] = FR
				SPY = 1;
		}
		if (SPY == 0){
			printf("Error occured when calculating the HLLC flux! check please!");
		}
	}
	free(QL), free(QR);
	free(VARIABLEL), free(VARIABLER);
}

void EvaluatePhysicalVariableByDepthThreshold(double hmin, double *h, double *variable, double *outPut){
	if (*h > hmin){
		*outPut = *variable / *h;
	}
	else{
		*outPut = 0;
	}
}