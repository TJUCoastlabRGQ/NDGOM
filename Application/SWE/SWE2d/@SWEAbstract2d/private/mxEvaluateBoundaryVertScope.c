#include "../../../../../NdgMath/NdgMath.h"
#include "../../../../../NdgMath/NdgSWE.h"

void EvaluateFaceRiemannProblem(double *, double *, double *, double *, double *, \
	double *, double , int , int , int , int );

void GetBCInvolvedLimitScope(double *, double *, double *, int , double *, double *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/*Allocate memory for fm and fp defined over boundary edges. Here, variables correspond to hu, hv, h, hT, hS,
	sediment and other passive transport material, vertical velocity omega is not included for boundary edges*/
	double *fphys = mxGetPr(prhs[0]);
	const mwSize *dims = mxGetDimensions(prhs[0]);
	const size_t Np = dims[0];
	const size_t K = dims[1];
	int BENfp = (int)mxGetScalar(prhs[1]);
	int BENe = (int)mxGetScalar(prhs[2]);
	int Nvar = (int)mxGetScalar(prhs[3]);
	double *BEFToE = mxGetPr(prhs[4]);
	double *BEFToN1 = mxGetPr(prhs[5]);
	double *BEnx = mxGetPr(prhs[6]);
	double *BEny = mxGetPr(prhs[7]);
	double *varFieldIndex = mxGetPr(prhs[8]);
	signed char *ftype = (signed char *)mxGetData(prhs[9]);
	//For the input exterior value, h comes first, then hu and finally hv, we need to consider this when impose boundary condition
	double *Tempfext = mxGetPr(prhs[10]);
	double *fext = malloc(BENfp*BENe*Nvar*sizeof(double));
	//hu first
	memcpy(fext, Tempfext + BENfp*BENe, BENfp*BENe*sizeof(double));
	//hv next
	memcpy(fext + BENfp*BENe, Tempfext + 2 * BENfp*BENe, BENfp*BENe*sizeof(double));
	//h next
	memcpy(fext + 2*BENfp*BENe, Tempfext, BENfp*BENe*sizeof(double));
	for (int i = 3; i < Nvar; i++){
		memcpy(fext + i * BENfp*BENe, Tempfext + i * BENfp*BENe, BENfp*BENe*sizeof(double));
	}

	double *fvmax = mxGetPr(prhs[11]);
	double *fvmin = mxGetPr(prhs[12]);
	double Hcrit = mxGetScalar(prhs[13]);
	double gra = mxGetScalar(prhs[14]);
	int Nh = (int)mxGetScalar(prhs[15]);
	int Nv = (int)mxGetScalar(prhs[16]);
	double *BEFToF = mxGetPr(prhs[17]);
	double *BEFToV = mxGetPr(prhs[18]);

	size_t NdimOut = 3;
	mwSize dimOut[3] = { Nv, 1, Nvar };
	plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *fmax = mxGetPr(plhs[0]);
	double *fmin = mxGetPr(plhs[1]);
	memcpy(fmax, fvmax, Nv*Nvar*sizeof(double));
	memcpy(fmin, fvmin, Nv*Nvar*sizeof(double));

	double *h = fphys, *hu = fphys + Np*K, \
		*hv = fphys + 2 * Np*K, *z = fphys + 3 * Np*K;

	double *fm = malloc(BENfp*BENe*Nvar*sizeof(double));
	double *huM = fm, *hvM = fm + BENfp*BENe, *hM = fm + 2 * BENfp*BENe;
	double *fp = malloc(BENfp*BENe*Nvar*sizeof(double));
	double *zM = malloc(BENfp*BENe*sizeof(double));
	double *zP = malloc(BENfp*BENe*sizeof(double));
	/*Since we only need the vertex value, so we choose to use the vertex value only*/
	double *fRiemann = malloc(2 * BENe*Nvar*sizeof(double));

	/*Fetch variable fm and fp first, then impose boundary condition and conduct hydrostatic reconstruction.
	Finally, calculate the local Riemann problem to get variable at the interface*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
		FetchBoundaryEdgeFacialValue(huM + face*BENfp, hu, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(hvM + face*BENfp, hv, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(zM + face*BENfp, z, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
		here 1 stands for the memory occupied by water depth h*/
		for (int field = 3; field < Nvar; field++){
			FetchBoundaryEdgeFacialValue(fm + field*BENe*BENfp + face*BENfp, \
				fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
				BEFToE + 2 * face, BEFToN1 + BENfp*face, Np, BENfp);
		}
		ImposeBoundaryCondition(&gra, type, BEnx + face*BENfp, BEny + face*BENfp, fm + face*BENfp, fp + face*BENfp, \
			zM + face*BENfp, zP + face*BENfp, fext + face*BENfp, BENfp, Nvar, BENe, varFieldIndex);

		EvaluateHydroStaticReconstructValue(Hcrit, fm + face*BENfp, fp + face*BENfp, zM + face*BENfp, zP + face*BENfp, BENfp, Nvar, BENe);
		EvaluateFaceRiemannProblem(fRiemann + face * 2, fm + face*BENfp, fp + face*BENfp, \
			BEnx + face*BENfp, BEny + face*BENfp, &gra, Hcrit, BENe, BENfp, Nvar, Nh);
	}

	/*
	Boundary Edge part to be considered here, this part is used to alter fmax and fmin
	*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int n = 0; n < Nvar; n++){
		GetBCInvolvedLimitScope(fmax + n*Nv, fmin + n*Nv, BEFToF, BENe, BEFToV, fRiemann + 2 * BENe * n);
	}
	free(fm);
	free(fp);
	free(zM);
	free(zP);
	free(fRiemann);
	free(fext);
}

void GetBCInvolvedLimitScope(double *fmax, double *fmin, double *FToF, int Ne, double *FToV, double *fRiemann){
	int v;
	for (int i = 0; i < Ne; i++){
			v = (int)FToV[2 * i] - 1;
			fmax[v] = max(fmax[v], fRiemann[2 * i]);
			fmin[v] = min(fmin[v], fRiemann[2 * i]);

			v = (int)FToV[2 * i + 1] - 1;
			fmax[v] = max(fmax[v], fRiemann[2 * i+1]);
			fmin[v] = min(fmin[v], fRiemann[2 * i+1]);
	}
}

void EvaluateFaceRiemannProblem(double *dest, double *Tempfm, double *Tempfp, double *Tempnx, double *Tempny, \
	double *gra, double Hcrit, int Ne, int Nfp, int Nvar, int Nh){
	/*The HLLC numerical flux presented in LAI(2012, Modeling one- and two-dimensional shallow water flows with discontinuous Galerkin method) is adopted*/
	double HR, HL, UR, UL, VR, VL;
	double SL, SM, SR;
	double USTAR, HSTAR;
	/*This part is used to extract the vertex value*/
	double *fm = malloc(2 * Nvar *sizeof(double)), *fp = malloc(2 * Nvar *sizeof(double));
	double *hum = fm, *hvm = fm + 2, *hup = fp, *hvp = fp + 2, *hm = fm + 2 * 2, *hp = fp + 2 * 2;
	double *nx = malloc(2 * sizeof(double)), *ny = malloc(2 * sizeof(double));
	nx[0] = Tempnx[0], nx[1] = Tempnx[Nh];
	ny[0] = Tempny[0], ny[1] = Tempny[Nh];
	for (int n = 0; n < Nvar; n++){
			fm[n * 2] = Tempfm[n*Nfp*Ne];
			fp[n * 2] = Tempfp[n*Nfp*Ne];
			fm[n * 2 + 1] = Tempfm[n*Nfp*Ne + Nh];
			fp[n * 2 + 1] = Tempfp[n*Nfp*Ne + Nh];
	}

	double *QL = malloc(Nvar*sizeof(double)), *QR = malloc(Nvar*sizeof(double));
	double *VARIABLEL = malloc(Nvar*sizeof(double)), *VARIABLER = malloc(Nvar*sizeof(double));
	double *QSTARL = malloc(Nvar*sizeof(double)), *QSTARR = malloc(Nvar*sizeof(double));
	for (int i = 0; i < 2; i++){
		/*Rotate variable to normal and tangential direction*/
		/*Here Q stands for H, HUn, HUt, H theta*/
		QL[0] = *(hm + i), QR[0] = *(hp + i);
		RotateFluxToNormal2d(hum + i, hvm + i, nx + i, ny + i, QL + 1, QL + 2);
		RotateFluxToNormal2d(hup + i, hvp + i, nx + i, ny + i, QR + 1, QR + 2);
		for (int n = 3; n < Nvar; n++)
		{
			QL[n] = *(fm + 2 * n + i);
			QR[n] = *(fp + 2 * n + i);
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
				/*AND FINALLY THE HLLC FLUX (AFTER ROTATION)*/
				if (SL >= 0){
					//*Fhn = FL[0];
					dest[i] = QL[0];
					dest[i + 2 * Ne] = (QL[1] * nx[i] - QL[2] * ny[i]);
					dest[i + 4 * Ne] = (QL[1] * ny[i] + QL[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QL[n];
					}
					SPY = 1;
				}
				else if (SM >= 0 && SL < 0){
					dest[i] = QSTARL[0];
					dest[i + 2 * Ne] = (QSTARL[1] * nx[i] - QSTARL[2] * ny[i]);
					dest[i + 4 * Ne] = (QSTARL[1] * ny[i] + QSTARL[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QSTARL[n];
					}
					SPY = 1;
				}
				else if (SM < 0 && SR > 0){
					dest[i] = QSTARR[0];
					dest[i + 2 * Ne] = (QSTARR[1] * nx[i] - QSTARR[2] * ny[i]);
					dest[i + 4 * Ne] = (QSTARR[1] * ny[i] + QSTARR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QSTARR[n];
					}
					SPY = 1;
				}
				else if (0 >= SR){
					dest[i] = QR[0];
					dest[i + 2 * Ne] = (QR[1] * nx[i] - QR[2] * ny[i]);
					dest[i + 4 * Ne] = (QR[1] * ny[i] + QR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QR[n];
					}
					SPY = 1;
				}

			}
			else if (HL>Hcrit && HR <= Hcrit){
				/*For this situation, the three wave structure has degenerated to the two wave structure.
				*For this special case, i.e. HR <= Hcrit, SM = SR.
				*/
				SL = UL - sqrt(*gra*HL);
				SR = UL + 2 * sqrt(*gra*HL);
				if (SL >= 0){
					dest[i] = QL[0];
					dest[i + 2 * Ne] = (QL[1] * nx[i] - QL[2] * ny[i]);
					dest[i + 4 * Ne] = (QL[1] * ny[i] + QL[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QL[n];
					}
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					QSTARL[0] = HL *(SL - UL) / (SL - SR);
					QSTARL[1] = HL *(SL - UL) / (SL - SR)*SR;
					QSTARL[2] = HL *(SL - UL) / (SL - SR)*VL;
					for (int i = 3; i<Nvar; i++){
						QSTARL[i] = HL*((SL - UL) / (SL - SR)) * VARIABLEL[i];
					}
					dest[i] = QSTARL[0];
					dest[i + 2 * Ne] = QSTARL[1] * nx[i] - QSTARL[2] * ny[i];
					dest[i + 4 * Ne] = QSTARL[1] * ny[i] + QSTARL[2] * nx[i];
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QSTARL[n];
					}
					SPY = 1;
				}
				else if (SR <= 0){
					/*Actually, I wonder whether this situation would occur*/
					dest[i] = QR[0];
					dest[i + 2 * Ne] = (QR[1] * nx[i] - QR[2] * ny[i]);
					dest[i + 4 * Ne] = (QR[1] * ny[i] + QR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QR[n];
					}
					SPY = 1;
				}
			}
			else if (HL <= Hcrit && HR > Hcrit){
				/*
				* For this special case, i.e. HL <= Hcrit, SM = SL.
				*/
				SL = UR - 2 * sqrt(*gra*HR);
				SR = UR + sqrt(*gra*HR);
				if (SL >= 0){
					dest[i] = QL[0];
					dest[i + 2 * Ne] = (QL[1] * nx[i] - QL[2] * ny[i]);
					dest[i + 4 * Ne] = (QL[1] * ny[i] + QL[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QL[n];
					}
					SPY = 1;
				}
				else if (SL < 0 && 0 < SR){
					QSTARR[0] = HR *(SR - UR) / (SR - SL);
					QSTARR[1] = HR *(SR - UR) / (SR - SL)*SL;
					QSTARR[2] = HR *(SR - UR) / (SR - SL)*VR;
					for (int i = 3; i<Nvar; i++){
						QSTARR[i] = HR*((SR - UR) / (SR - SL)) * VARIABLER[i];
					}
					dest[i] = QSTARR[0];
					dest[i + 2 * Ne] = (QSTARR[1] * nx[i] - QSTARR[2] * ny[i]);
					dest[i + 4 * Ne] = (QSTARR[1] * ny[i] + QSTARR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QSTARR[n];
					}
					SPY = 1;
				}
				else if (SR <= 0){
					dest[i] = QR[0];
					dest[i + 2 * Ne] = (QR[1] * nx[i] - QR[2] * ny[i]);
					dest[i + 4 * Ne] = (QR[1] * ny[i] + QR[2] * nx[i]);
					for (int n = 3; n < Nvar; n++){
						dest[i + n * 2 * Ne] = QR[n];
					}
					SPY = 1;
				}
			}
			//SM = (SL*HR*(UR - SR) - SR*HL*(UL - SL)) / (HR*(UR - SR) - HL*(UL - SL));
			if (SPY == 0){
				printf("Error occured when calculating the HLLC Riemann problem! check please!");
			}
		}
	}
	free(QL), free(QR);
	free(QSTARL), free(QSTARR);
	free(VARIABLEL), free(VARIABLER);
	free(fm), free(fp);
	free(nx), free(ny);
}
