#include "NdgSWE3D.h"
#include "NdgSWE.h"

void EvaluateVerticalFaceSurfFlux(double *dest, double *fm, double *nx, double *ny, double *gra, double Hcrit, int Nfp, int Nvar, int Ne){
	double u, v, theta;
	double *hu = fm, *hv = fm + Nfp*Ne, *h = fm + 2 * Nfp*Ne;
	for (int i = 0; i < Nfp; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hu + i, &u);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i, hv + i, &v);
		dest[i] = (hu[i] * u + 0.5 * (*gra)*pow(h[i],2))*nx[i] + hu[i] * v*ny[i];
		dest[i + Ne*Nfp] = hu[i] * v*nx[i] + (hv[i] * v + 0.5 * (*gra) * pow(h[i], 2))*ny[i];
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

void EvaluateHorizontalFaceNumFlux(double *FluxS, double *fm, double *fp, double *nz, double Hcrit, int Nfp, int Nvar, int Ne){
	double *hum = fm, *hvm = fm + Nfp*Ne, *omegam = fm + 2 * Nfp*Ne, *Hm = fm + 3 * Nfp*Ne;
	double *hup = fp, *hvp = fp + Nfp*Ne, *omegap = fp + 2 * Nfp*Ne, *Hp = fp + 3 * Nfp*Ne;
	double um, vm, up, vp, variablem, variablep, thetam, thetap;
	for (int i = 0; i < Nfp; i++){
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, hum + i, &um);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, hvm + i, &vm);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, hup + i, &up);
		EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, hvp + i, &vp);
		*(FluxS + i) = 0.5*(um*omegam[i] + up*omegap[i]) * nz[i];
		*(FluxS + Ne*Nfp + i) = 0.5*(vm*omegam[i] + vp*omegap[i]) * nz[i];
		for (int n = 2; n < Nvar; n++){
			/*Here 2*Ne*Nfp stands for the memory occupied by h and omega*/
			variablem = *(fm + 2 * Ne*Nfp + n*Ne*Nfp + i);
			variablep = *(fp + 2 * Ne*Nfp + n*Ne*Nfp + i);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hm + i, &variablem, &thetam);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, Hp + i, &variablep, &thetap);
			*(FluxS + n*Ne*Nfp + i) = 0.5*(thetam*omegam[i] + thetap*omegap[i]) * nz[i];
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
		Edest[i] = hu[i] * u + 0.5 * (*gra)*(pow(h[i], 2) - pow(z[i], 2));
		Gdest[i] = hu[i] * v;
		Hdest[i] = u*omega[i];
		Edest[i + Np*K] = hu[i] * v;
		Gdest[i + Np*K] = hv[i] * v + 0.5 * (*gra)*(pow(h[i], 2) - pow(z[i], 2));
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