#include "NdgSWE.h"

//This is used to impose the boundary condition for the pure hydrulic problem
void ImposeBoundaryCondition(double *gra, NdgEdgeType type, double *nx, double *ny, double *hM, double *huM, double *hvM, \
	double *zM, double *hE, double *huE, double *hvE, double *hP, double *huP, double *hvP, double *zP){
	// assign the local node values
	*zP = *zM;

	// get next node values
	if (type == NdgEdgeInner) {
		*hP = *hM;
		*huP = *huM;
		*hvP = *hvM;
	}
	else if (type == NdgEdgeZeroGrad) {
		*hP = *hM;
		*huP = *huM;
		*hvP = *hvM;
	}
	else if (type == NdgEdgeClamped) {
		*hP = *hE;
		*huP = *huE;
		*hvP = *hvE;
	}
	else if (type == NdgEdgeClampedDepth) {
		*hP = *hE;
		*huP = *huM;
		*hvP = *hvM;
	}
	else if (type == NdgEdgeClampedVel) {
		*hP = *hM;
		*huP = *huE;
		*hvP = *hvE;
	}
	else if (type == NdgEdgeSlipWall) {
		const double qxM = *huM;
		const double qyM = *hvM;
		double qnM = qxM * (*nx) + qyM * (*ny);   // outward normal flux
		double qvM = -qxM * (*ny) + qyM * (*nx);  // outward tangential flux
		// adjacent value
		*hP = *hM;
		*huP = (-qnM) * (*nx) - qvM * (*ny);
		*hvP = (-qnM) * (*ny) + qvM * (*nx);
	}
	else if (type == NdgEdgeNonSlipWall) {
		*hP = *hM;
		*huP = -*huM;
		*hvP = -*hvM;
	}
	else if (type == NdgEdgeFlather) {
		const double uE = *huE / *hE;
		const double vE = *hvE / *hE;
		const double unE = uE * (*nx) + vE * (*ny);   // outward normal flux
		const double uvE = -uE * (*ny) + vE * (*nx);  // tangential flux
		const double un = unE - sqrt(-*gra / *zP) * (*hE - *hM);
		const double uv = uvE;
		*hP = *hM;
		*huP = (un * (*nx) - uv * (*ny)) * (*hM);
		*hvP = (un * (*ny) + uv * (*nx)) * (*hM);
	}
	else if (type == NdgEdgeNonLinearFlather) {
		const double uE = *huE / *hE;
		const double vE = *hvE / *hE;
		const double unE = uE * (*nx) + vE * (*ny);  // outward normal flux
		const double un = *huM / *hM * (*nx) + *hvM / *hM * (*ny);
		const double temp = 0.5 * (un - unE) + sqrt(*gra * (*hE));
		*hP = temp * temp / *gra;
		*huP = *huE;
		*hvP = *hvE;
	}
	else if (type == NdgEdgeNonLinearFlatherFlow) {
		const double uE = *huE / *hE;
		const double vE = *hvE / *hE;
		const double unE = uE * (*nx) + vE * (*ny);   // outward normal flux
		const double uvE = -uE * (*ny) + vE * (*nx);  // tangential flux
		const double un = unE - 2 * sqrt(*gra * (*hE)) + 2 * sqrt(*gra * (*hM));
		const double uv = uvE;
		*hP = *hM;
		*huP = (un * (*nx) - uv * (*ny)) * (*hM);
		*hvP = (un * (*ny) + uv * (*nx)) * (*hM);
		//   } else if (type == NdgEdgeNonReflectingFlux) {
		//     const double unM = surf->huM / surf->hM * nx + surf->hvM / surf->hM *
		//     ny; const double RLP = unM + 2 * sqrt(gra * surf->hM); const double
		//     hs = fext->h[idM]; surf->hP = hs; const double RRM = RLP - 4 *
		//     sqrt(gra * hs); const double un = 0.5 * (RRM + RLP); const double uvM
		//     = -surf->huM / surf->hM * ny + surf->hvM / surf->hM * nx; surf->huP =
		//     (un * nx - uvM * ny) * surf->hM; surf->hvP = (un * ny + uvM * nx) *
		//     surf->hM;
	}
	else {
		printf("Matlab:%s:Unknown boundary type: %d\n", __FILE__, type);
	}
	return;
}

void EvaluateHydroStaticReconstructValue(double hmin, double *hM, double *huM, double *hvM, double *zM, double *hP, \
double *huP, double *hvP, double *zP){
	double zstar = max(*zM, *zP);
	double um, vm, up, vp;
	EvaluateFlowRateByDeptheThreshold(hmin, hM, huM, hvM, &um, &vm);
	EvaluateFlowRateByDeptheThreshold(hmin, hP, huP, hvP, &up, &vp);
	const double etaM = *hM + *zM;
	const double etaP = *hP + *zP;
	zstar = min(etaM, zstar);  // z* = min( \eta^-, z* )
	*huM = *hM * um;
	*hvM = *hM * vm;
	*huP = *hP * up;
	*hvP = *hP * vp;
	*hM = etaM - zstar;
	*hP = max(0, etaP - zstar) - max(0, *zP - zstar);
	*zM = zstar;
	*zP = zstar;
}

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