#include "SWENonhydrostatic3d.h"

/*This function is used to update the final velocity following:
$\left (hu\right )^{n+1}=\left (hu\right )^*-\frac{h\Delta t}{\rho_0}\left(\frac{\partial q}{\partial x}+\frac{\partial \sigma}{\partial x}\frac{\partial q}{\partial \sigma}\right)$
$\left (hv\right )^{n+1}=\left (hv\right )^*-\frac{h\Delta t}{\rho_0}\left(\frac{\partial q}{\partial y}+\frac{\partial \sigma}{\partial y}\frac{\partial q}{\partial \sigma}\right)$
$\left (hw\right )^{n+1}=\left (hw\right )^*-\frac{\Delta t}{\rho_0}\frac{\partial q}{\partial \sigma}$
The input parameters are organized as follows:
NonhydroPressure: The nonhydrostatic pressure, indexed as 0;
fphys: The three-dimensional physical field, indexed as 1;
varIndex: The variable index, indexed as 2;
rho: The water density, indexed as 3;
dt: The time step, indexed as 4;
PSPX: $\frac{\partial \sigma}{\partial x}$, indexed as 5;
PSPY: $\frac{\partial \sigma}{\partial y}$, indexed as 6;
mesh: The three-dimensional mesh object, indexed as 7;
cell: The three-dimensional master cell, indexed as 8;
InnerEdge: The three-dimensional inner edge object, indexed as 9;
BoundaryEdge: The three-dimensional boundary edge object, indexed as 10;
BottomEdge: The three-dimensional bottom edge object, indexed as 11;
BottomBoundaryEdge: The three-dimensional bottom boundary edge object, indexed as 12;
SurfaceBoundaryEdge: The three-dimensional surface boundary edge object, indexed as 13;
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *NonhydroPressure = mxGetPr(prhs[0]);
	int Np = (int)mxGetM(prhs[0]);
	int K = (int)mxGetN(prhs[0]);

	double *fphys = mxGetPr(prhs[1]);
	double *varIndex = mxGetPr(prhs[2]);
	double rho = mxGetScalar(prhs[3]);
	double dt = mxGetScalar(prhs[4]);

	double *PSPX = mxGetPr(prhs[5]);
	double *PSPY = mxGetPr(prhs[6]);

	double *hu = fphys + ((int)varIndex[0] - 1)*Np*K;
	double *hv = fphys + ((int)varIndex[1] - 1)*Np*K;
	double *hw = fphys + ((int)varIndex[2] - 1)*Np*K;
	double  *h = fphys + ((int)varIndex[3] - 1)*Np*K;

	const size_t ndimOut = 3;
	const mwSize dimOut[3] = { Np, K, 3 };
	plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutVariable = mxGetPr(plhs[0]);
	memcpy(OutVariable, hu, Np*K*sizeof(double));
	memcpy(OutVariable + Np*K, hv, Np*K*sizeof(double));
	memcpy(OutVariable + 2*Np*K, hw, Np*K*sizeof(double));

	mxArray *Temprx = mxGetField(prhs[7], 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(prhs[7], 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(prhs[7], 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *Tempsy = mxGetField(prhs[7], 0, "sy");
	double *sy = mxGetPr(Tempsy);
	mxArray *Temptz = mxGetField(prhs[7], 0, "tz");
	double *tz = mxGetPr(Temptz);
	mxArray *TempJ = mxGetField(prhs[7], 0, "J");
	double *J = mxGetPr(TempJ);

	mxArray *TempDr = mxGetField(prhs[8], 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(prhs[8], 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(prhs[8], 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempNface = mxGetField(prhs[8], 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);// Substract the bottom face and surface face from the total face
	mxArray *TempInvM = mxGetField(prhs[8], 0, "invM");
	double *invM = mxGetPr(TempInvM);

	mxArray *TempIENe = mxGetField(prhs[9], 0, "Ne");
	int IENe = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp = mxGetField(prhs[9], 0, "Nfp");
	int IENfp = (int)mxGetScalar(TempIENfp);
	mxArray *TempIEMb = mxGetField(prhs[9], 0, "M");
	double *IEMb = mxGetPr(TempIEMb);
	mxArray *TempIEJs = mxGetField(prhs[9], 0, "Js");
	double *IEJs = mxGetPr(TempIEJs);
	mxArray *TempIEnx = mxGetField(prhs[9], 0, "nx");
	double *IEnx = mxGetPr(TempIEnx);
	mxArray *TempIEny = mxGetField(prhs[9], 0, "ny");
	double *IEny = mxGetPr(TempIEny);
	mxArray *TempIELAV = mxGetField(prhs[9], 0, "LAV");
	double *IELAV = mxGetPr(TempIELAV);
	mxArray *TempIEFToE = mxGetField(prhs[9], 0, "FToE");
	double *IEFToE = mxGetPr(TempIEFToE);
	mxArray *TempIEFToF = mxGetField(prhs[9], 0, "FToF");
	double *IEFToF = mxGetPr(TempIEFToF);
	mxArray *TempIEFToN1 = mxGetField(prhs[9], 0, "FToN1");
	double *IEFToN1 = mxGetPr(TempIEFToN1);
	mxArray *TempIEFToN2 = mxGetField(prhs[9], 0, "FToN2");
	double *IEFToN2 = mxGetPr(TempIEFToN2);

	mxArray *TempBENe = mxGetField(prhs[10], 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBENfp = mxGetField(prhs[10], 0, "Nfp");
	int BENfp = mxGetScalar(TempBENfp);
	mxArray *TempBEMb = mxGetField(prhs[10], 0, "M");
	double *BEMb = mxGetPr(TempBEMb);
	mxArray *TempBEJs = mxGetField(prhs[10], 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	mxArray *TempBEnx = mxGetField(prhs[10], 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(prhs[10], 0, "ny");
	double *BEny = mxGetPr(TempBEny);
	mxArray *TempBELAV = mxGetField(prhs[10], 0, "LAV");
	double *BELAV = mxGetPr(TempBELAV);
	mxArray *TempBEFToE = mxGetField(prhs[10], 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToF = mxGetField(prhs[10], 0, "FToF");
	double *BEFToF = mxGetPr(TempBEFToF);
	mxArray *TempBEFToN1 = mxGetField(prhs[10], 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);
	mxArray *Tempftype = mxGetField(prhs[10], 0, "ftype");
	signed char *ftype = (signed char *)mxGetData(TempBEFToN1);

	/*For bottom edge object*/
	mxArray *TempBotENe = mxGetField(prhs[11], 0, "Ne");
	int BotENe = (int)mxGetScalar(TempBotENe);
	mxArray *TempBotENfp = mxGetField(prhs[11], 0, "Nfp");
	int BotENfp = (int)mxGetScalar(TempBotENfp);
	mxArray *TempBotEMb = mxGetField(prhs[11], 0, "M");
	double *BotEMb = mxGetPr(TempBotEMb);
	mxArray *TempBotEJs = mxGetField(prhs[11], 0, "Js");
	double *BotEJs = mxGetPr(TempBotEJs);
	mxArray *TempBotEnz = mxGetField(prhs[11], 0, "nz");
	double *BotEnz = mxGetPr(TempBotEnz);
	mxArray *TempBotEFToE = mxGetField(prhs[11], 0, "FToE");
	double *BotEFToE = mxGetPr(TempBotEFToE);
	mxArray *TempBotEFToF = mxGetField(prhs[11], 0, "FToF");
	double *BotEFToF = mxGetPr(TempBotEFToF);
	mxArray *TempBotEFToN1 = mxGetField(prhs[11], 0, "FToN1");
	double *BotEFToN1 = mxGetPr(TempBotEFToN1);
	mxArray *TempBotEFToN2 = mxGetField(prhs[11], 0, "FToN2");
	double *BotEFToN2 = mxGetPr(TempBotEFToN2);

	mxArray *TempBotBENe = mxGetField(prhs[12], 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBENfp = mxGetField(prhs[12], 0, "Nfp");
	int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	mxArray *TempBotBEMb = mxGetField(prhs[12], 0, "M");
	double *BotBEMb = mxGetPr(TempBotBEMb);
	mxArray *TempBotBEJs = mxGetField(prhs[12], 0, "Js");
	double *BotBEJs = mxGetPr(TempBotBEJs);
	mxArray *TempBotBEnz = mxGetField(prhs[12], 0, "nz");
	double *BotBEnz = mxGetPr(TempBotBEnz);
	mxArray *TempBotBEFToE = mxGetField(prhs[12], 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToF = mxGetField(prhs[12], 0, "FToF");
	double *BotBEFToF = mxGetPr(TempBotBEFToF);
	mxArray *TempBotBEFToN1 = mxGetField(prhs[12], 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);

	/*For surface boundary edge object*/
	mxArray *TempSurfBENe = mxGetField(prhs[13], 0, "Ne");
	int SurfBENe = (int)mxGetScalar(TempSurfBENe);
	mxArray *TempSurfBENfp = mxGetField(prhs[13], 0, "Nfp");
	int SurfBENfp = (int)mxGetScalar(TempSurfBENfp);
	mxArray *TempSurfBEMb = mxGetField(prhs[13], 0, "M");
	double *SurfBEMb = mxGetPr(TempSurfBEMb);
	mxArray *TempSurfBEJs = mxGetField(prhs[13], 0, "Js");
	double *SurfBEJs = mxGetPr(TempSurfBEJs);
	mxArray *TempSurfBEnz = mxGetField(prhs[13], 0, "nz");
	double *SurfBEnz = mxGetPr(TempSurfBEnz);
	mxArray *TempSurfBEFToE = mxGetField(prhs[13], 0, "FToE");
	double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
	mxArray *TempSurfBEFToF = mxGetField(prhs[13], 0, "FToF");
	double *SurfBEFToF = mxGetPr(TempSurfBEFToF);
	mxArray *TempSurfBEFToN1 = mxGetField(prhs[13], 0, "FToN1");
	double *SurfBEFToN1 = mxGetPr(TempSurfBEFToN1);

	double *PNPX = malloc(Np*K*sizeof(double));
	double *PNPY = malloc(Np*K*sizeof(double));
	double *PNPS = malloc(Np*K*sizeof(double));

	double *qM = malloc(IENe*IENfp*sizeof(double));
	double *qP = malloc(IENe*IENfp*sizeof(double));
	double *QIEfluxSx = malloc(IENe*IENfp*sizeof(double));
	double *QIEfluxSy = malloc(IENe*IENfp*sizeof(double));
	double *NonhydroERHS = malloc(Np*K * 2 * Nface*sizeof(double));
	double *QIEfluxMx = malloc(IENe*IENfp*sizeof(double));
	double *QIEfluxMy = malloc(IENe*IENfp*sizeof(double));
	double *QIEfluxPx = malloc(IENe*IENfp*sizeof(double));
	double *QIEfluxPy = malloc(IENe*IENfp*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		/*Fetch variable IEfm and IEfp first*/
		FetchInnerEdgeFacialValue(qM + face*IENfp, qP + face*IENfp, NonhydroPressure, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(QIEfluxMx + face*IENfp, qM + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(QIEfluxPx + face*IENfp, qP + face*IENfp, IEnx + face*IENfp, IENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(QIEfluxMy + face*IENfp, qM + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(QIEfluxPy + face*IENfp, qP + face*IENfp, IEny + face*IENfp, IENfp);

		EvaluateNonhydroVerticalFaceNumFlux_Central(QIEfluxSx + face*IENfp, qM + face*IENfp, qP + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(QIEfluxSy + face*IENfp, qM + face*IENfp, qP + face*IENfp, IEny + face*IENfp, IENfp);
	}

	memset(NonhydroERHS, 0, Np*K * 2 * Nface*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, QIEfluxMx, \
			QIEfluxPx, QIEfluxSx, IEJs, IEMb, NonhydroERHS);
		StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, QIEfluxMy, \
			QIEfluxPy, QIEfluxSy, IEJs, IEMb, NonhydroERHS + Np*K*Nface);
	}

	double *qbM = malloc(BENe*BENfp*sizeof(double));
	double *QBEfluxMx = malloc(BENe*BENfp*sizeof(double));
	double *QBEfluxMy = malloc(BENe*BENfp*sizeof(double));
	double *QBEfluxSx = malloc(BENe*BENfp*sizeof(double));
	double *QBEfluxSy = malloc(BENe*BENfp*sizeof(double));
	double *NonhydroTempFacialIntegral = malloc(Np*K*sizeof(double));
	double *TempVolumeIntegral = malloc(Np*K*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition, unused at present
		FetchBoundaryEdgeFacialValue(qbM + face*BENfp, NonhydroPressure, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(QBEfluxMx + face*BENfp, qbM + face*BENfp, BEnx + face*BENfp, BENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(QBEfluxMy + face*BENfp, qbM + face*BENfp, BEny + face*BENfp, BENfp);
		if (type == NdgEdgeSlipWall){
			DotProduct(QBEfluxSx + face*BENfp, qbM + face*SurfBENfp, BEnx + face*BENfp, BENfp);
			DotProduct(QBEfluxSy + face*BENfp, qbM + face*SurfBENfp, BEny + face*BENfp, BENfp);
		}
		else{//impose zero boundary condition at other boundary
			for (int i = 0; i < BENfp; i++){
				QBEfluxSx[face*BENfp + i] = 0;
				QBEfluxSy[face*BENfp + i] = 0;
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		StrongFormBoundaryEdgeRHS(face, BEFToE, BEFToF, Np, K, BENfp, BEFToN1, QBEfluxMx, QBEfluxSx, BEJs, BEMb, NonhydroERHS);
		StrongFormBoundaryEdgeRHS(face, BEFToE, BEFToF, Np, K, BENfp, BEFToN1, QBEfluxMy, QBEfluxSy, BEJs, BEMb, NonhydroERHS + Np*K*Nface);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field<2; field++){
			for (int face = 1; face<Nface; face++){
				Add(NonhydroERHS + field*Np*K*Nface + k*Np, NonhydroERHS + field*Np*K*Nface + k*Np, NonhydroERHS + field*Np*K*Nface + face*Np*K + k*Np, Np);
			}
		}
	}

	ptrdiff_t np = Np;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < 2; field++){

			MultiEdgeContributionByLiftOperator(NonhydroERHS + field*Np*K*Nface + k*Np, NonhydroTempFacialIntegral + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
		GetVolumnIntegral2d(PNPX + k*Np, TempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, NonhydroPressure + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);
		/*$\bold{r_y}\cdot (Dr*hv)+\bold{s_y}\cdot (Ds*hv)$*/
		GetVolumnIntegral2d(PNPY + k*Np, TempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, NonhydroPressure + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){

		Minus(PNPX + k*Np, PNPX + k*Np, NonhydroERHS + k*Np, Np);

		Minus(PNPY + k*Np, PNPY + k*Np, NonhydroERHS + Np*K*Nface + k*Np, Np);

	}

/*********************************************************Partial derivative in vertical direction*************************************************************************/

	double *qBotEM = malloc(BotENe*BotENfp*sizeof(double));
	double *qBotEP = malloc(BotENe*BotENfp*sizeof(double));
	double *QBotEfluxS = malloc(BotENe*BotENfp*sizeof(double));
	double *QBotEfluxM = malloc(BotENe*BotENfp*sizeof(double));
	double *QBotEfluxP = malloc(BotENe*BotENfp*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++){
		/*Fetch variable IEfm and IEfp first*/
		FetchInnerEdgeFacialValue(qBotEM + face*BotENfp, qBotEP + face*BotENfp, NonhydroPressure, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(QBotEfluxM + face*BotENfp, qBotEM + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(QBotEfluxP + face*BotENfp, qBotEP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);

		EvaluateNonhydroVerticalFaceNumFlux_Central(QBotEfluxS + face*BotENfp, qBotEM + face*BotENfp, qBotEP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
	}

	memset(NonhydroERHS, 0, Np*K * Nface *sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		StrongFormInnerEdgeRHS(face, BotEFToE, BotEFToF, Np, K, BotENfp, BotEFToN1, BotEFToN2, QBotEfluxM, \
			QBotEfluxP, QBotEfluxS, BotEJs, BotEMb, NonhydroERHS);
	}

	double *qBotBEM = malloc(BotBENe*BotBENfp*sizeof(double));
	double *QBotBEfluxM = malloc(BotBENe*BotBENfp*sizeof(double));
	double *QBotBEfluxS = malloc(BotBENe*BotBENfp*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		FetchBoundaryEdgeFacialValue(qBotBEM + face*BotBENfp, NonhydroPressure, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(QBotBEfluxM + face*BotBENfp, qBotBEM + face*BotBENfp, BotBEnz + face*BotBENfp, BotBENfp);
		DotProduct(QBotBEfluxS + face*BotBENfp, qBotBEM + face*BotBENfp, BotBEnz + face*BotBENfp, BotBENfp);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		StrongFormBoundaryEdgeRHS(face, BotBEFToE, BotBEFToF, Np, K, BotBENfp, BotBEFToN1, QBotBEfluxM, QBotBEfluxS, BotBEJs, BotBEMb, NonhydroERHS);
	}

	double *qSurfBEM = malloc(SurfBENe*SurfBENfp*sizeof(double));
	double *QSurfBEfluxM = malloc(SurfBENe*SurfBENfp*sizeof(double));
	double *QSurfBEfluxS = malloc(SurfBENe*SurfBENfp*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		FetchBoundaryEdgeFacialValue(qSurfBEM + face*SurfBENfp, NonhydroPressure, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(QSurfBEfluxM + face*SurfBENfp, qSurfBEM + face*SurfBENfp, SurfBEnz + face*SurfBENfp, SurfBENfp);
		/*At dirichelet boundary, the numerical flux is given as the Dirichlet boundary*/
		for (int p = 0; p < SurfBENfp; p++){
			QSurfBEfluxS[face*SurfBENfp + p] = 0;
		}
//		DotProduct(QSurfBEfluxS + face*SurfBENfp, qSurfBEM + face*SurfBENfp, SurfBEnz + face*SurfBENfp, SurfBENfp);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		StrongFormBoundaryEdgeRHS(face, SurfBEFToE, SurfBEFToF, Np, K, SurfBENfp, SurfBEFToN1, QSurfBEfluxM, QSurfBEfluxS, SurfBEJs, SurfBEMb, NonhydroERHS);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field<1; field++){
			for (int face = 1; face<Nface; face++){
				Add(NonhydroERHS + field*Np*K*Nface + k*Np, NonhydroERHS + field*Np*K*Nface + k*Np, NonhydroERHS + field*Np*K*Nface + face*Np*K + k*Np, Np);
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int field = 0; field < 1; field++){

			MultiEdgeContributionByLiftOperator(NonhydroERHS + field*Np*K*Nface + k*Np, NonhydroTempFacialIntegral + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
	for (int k = 0; k < K; k++){
		/*$\bold{t_z}\cdot (Dt*u)$*/
		GetVolumnIntegral1d(PNPS + k*Np, &np, &oneI, &np, &one, \
			Dt, &np, NonhydroPressure + k*Np, &np, &zero, &np, tz + k*Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){

		Minus(PNPS + k*Np, PNPS + k*Np, NonhydroERHS + k*Np, Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int p = 0; p < Np; p++){
			hu[k*Np + p] = hu[k*Np + p] - h[k*Np + p] * dt / rho*(PNPX[k*Np + p] + PSPX[k*Np + p] * PNPS[k*Np + p]);
			hv[k*Np + p] = hv[k*Np + p] - h[k*Np + p] * dt / rho*(PNPY[k*Np + p] + PSPY[k*Np + p] * PNPS[k*Np + p]);
			hw[k*Np + p] = hw[k*Np + p] - dt / rho*PNPS[k*Np + p];
		}
	}


	free(PNPX);
	free(PNPY);
	free(PNPS);

	free(qM);
	free(qP);
	free(QIEfluxSx);
	free(QIEfluxSy);
	free(NonhydroERHS);
	free(QIEfluxMx);
	free(QIEfluxMy);
	free(QIEfluxPx);
	free(QIEfluxPy);

	free(qbM);
	free(QBEfluxMx);
	free(QBEfluxMy);
	free(QBEfluxSx);
	free(QBEfluxSy);
	free(NonhydroTempFacialIntegral);
	free(TempVolumeIntegral);

	free(qBotEM);
	free(qBotEP);
	free(QBotEfluxS);
	free(QBotEfluxM);
	free(QBotEfluxP);

	free(qBotBEM);
	free(QBotBEfluxM);
	free(QBotBEfluxS);

	free(qSurfBEM);
	free(QSurfBEfluxM);
	free(QSurfBEfluxS);

}