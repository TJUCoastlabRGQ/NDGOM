#include "SWENonhydrostatic3d.h"

extern double  *NonhydroVariable, \
*PSPX, *PSPY, *TempPSPX, *TempPSPY, *NonhydroHU2d, *NonhydroHV2d,\
*TempNonhydroHU2d, *TempNonhydroHV2d, *TempNonhydrofield3d, *Nonhydrofmod,\
*NonhydroIEfm, *NonhydroIEfp, *NonhydroIEFluxM, *NonhydroIEFluxP, *NonhydroIEFluxS,\
*NonhydroERHS, *NonhydroBEfm, *NonhydroBEfp, *NonhydroBEFluxM, *NonhydroBEFluxS,\
*NonhydroTempFacialIntegral, *NonhydroTempVolumeIntegral, *NonhydrozM, *NonhydrozP,\
*CoePS, *NonhydroBotEfm, *NonhydroBotEfp, *NonhydroBotEFluxM, *NonhydroBotEFluxP,\
*NonhydroBotEFluxS, *NonhydroBotBEfm, *NonhydroBotBEFluxM, *NonhydroBotBEFluxS,\
*NonhydroSurfBEfm, *NonhydroSurfBEFluxM, *NonhydroSurfBEFluxS, *NonhydroRHS2d,\
*Weta, *Wbot, *ueta, *ubot, *veta, *vbot, *etax, *etay, *TempZx2d, *TempZy2d, *NonhydroIEfm2d,\
*NonhydroIEfp2d, *NonhydroIEFluxM2d, *NonhydroIEFluxP2d, *NonhydroIEFluxS2d, *NonhydroERHS2d,\
*NonhydroPCEVolumeIntegralX, *NonhydroPCETempVolumeIntegralX, *NonhydroPCEVolumeIntegralY,\
*NonhydroPCETempVolumeIntegralY, *NonhydroBEfm2d, *NonhydroBEzM2d, *NonhydroBEfp2d,\
*NonhydroBEzP2d, *NonhydroBEFluxS2d, *NonhydroBEFluxM2d, *NonhydroPCETempFacialIntegral,\
*NonhydroIEfmod, *NonhydroBEfmod;

//extern double *PHPX, *PHPY;

extern char *SWENonhydro3dInitialized;

void MyExit()
{
	if (!strcmp("True", SWENonhydro3dInitialized)){
		SWENonhydro3dMemoryDeAllocation();
		SWENonhydro3dInitialized = "False";
	}
	return;
}

void GetFirstOrderPartialDerivativeInHorizontalDirection( double *, double *, double *, double *, double *, \
	double *, const mxArray *, const mxArray *, const mxArray *, \
	const mxArray *, double *, double *, double *, double *, double *, double *, double, double *, double, signed char *);

void GetFirstOrderPartialDerivativeInVerticalDirection(double *, double *, double *, const mxArray *, const mxArray *, \
	const mxArray *, const mxArray *, const mxArray *, double *, double *, double *, double *, \
	double *);

void GetVerticalVelocityAtSurfaceAndBottom(double *, double *, const mxArray *, const mxArray *, \
	const mxArray *, const mxArray *, const mxArray *, const mxArray *, const mxArray *, \
	const mxArray *, const mxArray *, signed char *, signed char *, double *, double *, double *, \
	double *, double *, double *, double *, double *, double *, double *, double *, double *, double , \
	double , double *, double *, double *, double *, double *, double *);

void DotCriticalDevideByLocalValue(double *, double *, int , double );

void DotCriticalDevideByAveragedValue(double *, double *, double *, int , double );


/*The following function is used to assemble the final partial derivative of $\sigma$ with respect to $x$ or
$y$. They are explicitly given as $\frac{1}{H}\frac{\partial h}{\partial x}-\frac{\left (1+\sigma\right )}
{H}\frac{\partial H}{\partial x}$and $\frac{1}{H}\frac{\partial h}{\partial y}-\frac{\left (1+\sigma\right )}
{H}\frac{\partial H}{\partial y}$,with $h$ the still water depth such that $\frac{\partial h}{\partial x}$ 
is equal to $-\frac{\partial z}{\partial x}$
*/

void AssembleFirstOrderPartialDerivativeInSigma(double *dest, double *Tempdest, double *phpd, double *hcrit, double *zs, \
	double *h, int Np, double *Coe, double *z){
	/*$\frac{1}{H}\frac{\partial z}{\partial x}$*/
	DotCriticalDivide(Tempdest, zs, hcrit, h, Np);
	/*$-\frac{1}{H}\frac{\partial z}{\partial x}$*/
	DotDivideByConstant(Tempdest, Tempdest, -1.0, Np);
	/*$1+\sigma$*/
	AddByConstant(Coe, z, 1, Np);
	/*$\frac{1+\sigma}{H}$*/
	DotCriticalDivide(Coe, Coe, hcrit, h, Np);
	/*$-\frac{1+\sigma}{H}$*/
	DotDivideByConstant(Coe, Coe, -1.0, Np);
	/*$-\frac{1+\sigma}{H}\frac{\partial H}{\partial x}$*/
	DotProduct(dest, phpd, Coe, Np);
	/*$\frac{ 1 }{H}\frac{ \partial h }{\partial x}-\frac{ \left(1 + \sigma\right) }
	{H}\frac{ \partial H }{\partial x}$*/
	Add(dest, dest, Tempdest, Np);
}

/*The following function is used to calculate $\left(\frac{\partial \sigma}{\partial x}\right )^2$ and
$\left(\frac{\partial \sigma}{\partial y}\right )^2$
*/
void AssembleSecondOrderPartialDerivativeInSigma(double *dest, double *source, int Np){
	for (int i = 0; i < Np; i++){
		dest[i] = source[i] * source[i];
	}
}

/*This function is used to calculate the variable used in the non-hydrostatic solver, variables calculated
in this part include $\frac{\partial u}{\partial x}$, $\frac{\partial v}{\partial y}$, 
$\frac{\partial \sigma}{\partial x}$, $\frac{\partial \sigma}{\partial y}$,
$\left (\frac{\partial \sigma}{\partial x}\right )^2$, $\left (\frac{\partial \sigma}{\partial y}\right )^2$
$\frac{\partial^2 \sigma}{\partial x^2}$, $\frac{\partial^2 \sigma}{\partial y^2}$
The input variables are as follows:
hcirt: the critical water depth, the input order is 0.
mesh: the three-dimensional mesh object, the input order is 1.
cell: the three-dimensional master cell, the input order is 2.
InnerEdge: the three-dimensional inner edge object, the input order is 3.
BoundaryEdge: the three-dimensional boundary edge object, the input order is 4.
BottomEdge: the three-dimensional bottom edge object, the input order is 5.
BottomBoundaryEdge: the three-dimensional bottom boundary edge object, the input order is 6.
SurfaceBoundaryEdge: the three-dimensional surface boundary edge object, the input order is 7.
fphys: the three-dimensional physical field, the input order is 8.
varIndex: the variable index, the input order is 9.
ftype3d: the three-dimensional boundary edge type, the input order is 10.
gra: the accelaration term due to gravity, the input order is 11.
fext: the three-dimensional exterior data at the boundary edge, the input order is 12.
h2d: the two-dimensional water depth field, the input order is 13.
z2d: the two-dimensional bottom elevation field, the input order is 14.
fext2d: the two-dimensional exterior data at the boundary edge, the input order is 15.
mesh2d: the two-dimensional mesh object, the input order is 16.
InnerEdge2d: the two-dimensional inner edge object, the input order is 17.
BoundaryEdge2d: the two-dimensional boundary edge object, the input order is 18.
cell2d: the two-dimensional master cell, the input order is 19.
ftype2d: the two-dimensional boundary edge type, the input order is 20.
*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	mexAtExit(&MyExit);
	double *Hcrit = mxGetPr(prhs[0]);

	const mxArray *mesh = prhs[1];
	const mxArray *cell = prhs[2];
	const mxArray *InnerEdge = prhs[3];
	const mxArray *BoundaryEdge = prhs[4];
	const mxArray *BottomEdge = prhs[5];
	const mxArray *BottomBoundaryEdge = prhs[6];
	const mxArray *SurfaceBoundaryEdge = prhs[7];

	double *fphys = mxGetPr(prhs[8]);
	const size_t *PRHS;
	PRHS = mxGetDimensions(prhs[8]);
	int Np = (int)PRHS[0];
	int K = (int)PRHS[1];

	mxArray *Tempz = mxGetField(mesh, 0, "z");
	double *z = mxGetPr(Tempz);
	mxArray *TempNLayer = mxGetField(mesh, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);
	mxArray *TempJz = mxGetField(mesh, 0, "Jz");
	double *Jz = mxGetPr(TempJz);

	mxArray *TempV3d = mxGetField(cell, 0, "V");
	double *V3d = mxGetPr(TempV3d);
	double *InvV3d = malloc(Np*Np*sizeof(double));
	memcpy(InvV3d, V3d, Np*Np*sizeof(double));
	MatrixInverse(InvV3d, (ptrdiff_t)Np);

	//Nvar = (int)PRHS[2];
	double *varIndex = mxGetPr(prhs[9]);
	int var = (int)mxGetNumberOfElements(prhs[9]);
	double *hu = fphys + ((int)varIndex[0] - 1)*Np*K;
	double *hv = fphys + ((int)varIndex[1] - 1)*Np*K;
	double *hw = fphys + ((int)varIndex[2] - 1)*Np*K;
	double  *h = fphys + ((int)varIndex[3] - 1)*Np*K;
	double *zx = fphys + ((int)varIndex[4] - 1)*Np*K;
	double *zy = fphys + ((int)varIndex[5] - 1)*Np*K;
	double *zbot = fphys + ((int)varIndex[6] - 1)*Np*K;

	signed char *ftype = (signed char *)mxGetData(prhs[10]);
	double gra = mxGetScalar(prhs[11]);
	double *fext = mxGetPr(prhs[12]);
	double *h2d = mxGetPr(prhs[13]);
	double *z2d = mxGetPr(prhs[14]);
	double *fext2d = mxGetPr(prhs[15]);
	signed char *ftype2d = (signed char *)mxGetData(prhs[20]);

	const mxArray *mesh2d = prhs[16];
	mxArray *TempK2d = mxGetField(mesh2d, 0, "K");
	int K2d = (int)mxGetScalar(TempK2d);

	const mxArray *cell2d = prhs[19];
	mxArray *TempNp2d = mxGetField(cell2d, 0, "Np");
	int Np2d = (int)mxGetScalar(TempNp2d);
	mxArray *TempV2d = mxGetField(cell2d, 0, "V");
	double *V2d = mxGetPr(TempV2d);
	mxArray *TempNface2d = mxGetField(cell2d, 0, "Nface");
	int Nface2d = (int)mxGetScalar(TempNface2d);

	const mxArray *InnerEdge2d = prhs[17];
	const mxArray *BoundaryEdge2d = prhs[18];

	double *pept = mxGetPr(prhs[21]);

	if (!strcmp("False", SWENonhydro3dInitialized)){
	   mxArray *TempNface = mxGetField(cell, 0, "Nface");
	   int Nface3d = (int)mxGetScalar(TempNface);
	   mxArray *TempIENfp = mxGetField(InnerEdge, 0, "Nfp");
	   int IENfp = (int)mxGetScalar(TempIENfp);
	   mxArray *TempIENe = mxGetField(InnerEdge, 0, "Ne");
	   int IENe = (int)mxGetScalar(TempIENe);

	   mxArray *TempBENfp = mxGetField(BoundaryEdge, 0, "Nfp");
	   int BENfp = (int)mxGetScalar(TempBENfp);
	   mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	   int BENe = (int)mxGetScalar(TempBENe);

	   mxArray *TempBotENfp = mxGetField(BottomEdge, 0, "Nfp");
	   int BotENfp = (int)mxGetScalar(TempBotENfp);
	   mxArray *TempBotENe = mxGetField(BottomEdge, 0, "Ne");
	   int BotENe = (int)mxGetScalar(TempBotENe);

	   mxArray *TempBotBENfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
	   int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	   mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	   int BotBENe = (int)mxGetScalar(TempBotBENe);

	   mxArray *TempSurfBENfp = mxGetField(SurfaceBoundaryEdge, 0, "Nfp");
	   int SurfBENfp = (int)mxGetScalar(TempSurfBENfp);
	   mxArray *TempSurfBENe = mxGetField(SurfaceBoundaryEdge, 0, "Ne");
	   int SurfBENe = (int)mxGetScalar(TempSurfBENe);

	   mxArray *TempIENfp2d = mxGetField(InnerEdge2d, 0, "Nfp");
	   int IENfp2d = (int)mxGetScalar(TempIENfp2d);
	   mxArray *TempIENe2d = mxGetField(InnerEdge2d, 0, "Ne");
	   int IENe2d = (int)mxGetScalar(TempIENe2d);

	   mxArray *TempBENfp2d = mxGetField(BoundaryEdge2d, 0, "Nfp");
	   int BENfp2d = (int)mxGetScalar(TempBENfp2d);
	   mxArray *TempBENe2d = mxGetField(BoundaryEdge2d, 0, "Ne");
	   int BENe2d = (int)mxGetScalar(TempBENe2d);

		SWENonhydro3dMemoryAllocation(Np, K, IENfp, IENe, Nface3d, \
			 BENfp, BENe, BotENfp, BotENe, BotBENfp, BotBENe, SurfBENfp, \
			 SurfBENe, Np2d, K2d, IENfp2d, IENe2d, Nface2d, BENfp2d, \
			 BENe2d);
	}

	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PSPX = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PSPY = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *SQPSPX = mxGetPr(plhs[2]);
	plhs[3] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *SQPSPY = mxGetPr(plhs[3]);
	plhs[4] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PUPX = mxGetPr(plhs[4]);
	plhs[5] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PUPY = mxGetPr(plhs[5]);
	plhs[6] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PVPX = mxGetPr(plhs[6]);
	plhs[7] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PVPY = mxGetPr(plhs[7]);
	plhs[8] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PUPS = mxGetPr(plhs[8]);
	plhs[9] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PVPS = mxGetPr(plhs[9]);
	plhs[10] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PWPS = mxGetPr(plhs[10]);
	plhs[11] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
	double *Wnew = mxGetPr(plhs[11]);
	/*The horizontal velocity at the bottom boundary*/
	plhs[12] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
	double *Unew = mxGetPr(plhs[12]);
	/*The horizontal velocity at the bottom boundary*/
	plhs[13] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
	double *Vnew = mxGetPr(plhs[13]);

	plhs[14] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PHPX = mxGetPr(plhs[14]);

	plhs[15] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *PHPY = mxGetPr(plhs[15]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field < 3; field++){
			//For variable u, v, and w
			DotCriticalDivide(NonhydroVariable + field*Np*K + k*Np, \
				fphys + (int)(varIndex[field] - 1)*Np*K + k*Np, Hcrit, \
				h + k*Np, Np);
		}
	}

	double *u = NonhydroVariable;
	double *v = NonhydroVariable + Np*K;
	double *w = NonhydroVariable + 2*Np*K;


	GetFirstOrderPartialDerivativeInHorizontalDirection(PHPX, PHPY, PUPX, PUPY, PVPX, PVPY, mesh, cell, InnerEdge, BoundaryEdge, h, \
		u, v, hu, hv, zbot, gra, fext, *Hcrit, ftype);


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){

		AssembleFirstOrderPartialDerivativeInSigma(PSPX + k*Np, TempPSPX + k*Np, PHPX + k*Np, Hcrit, zx + k*Np, \
			h + k*Np, Np, CoePS + k*Np, z + k*Np);

		AssembleFirstOrderPartialDerivativeInSigma(PSPY + k*Np, TempPSPY + k*Np, PHPY + k*Np, Hcrit, zy + k*Np, \
			h + k*Np, Np, CoePS + k*Np, z + k*Np);

		AssembleSecondOrderPartialDerivativeInSigma(SQPSPX + k*Np, PSPX + k*Np, Np);

		AssembleSecondOrderPartialDerivativeInSigma(SQPSPY + k*Np, PSPY + k*Np, Np);

	}

	memset(TempNonhydroHU2d, 0, Np2d*K2d*sizeof(double));
	memset(TempNonhydroHV2d, 0, Np2d*K2d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++){
		VerticalColumnIntegralField3d(NonhydroHU2d + i*Np2d, Np2d, V2d, TempNonhydroHU2d + i*Np2d, \
			TempNonhydrofield3d + i*Np*NLayer, hu + i*Np*NLayer, Jz + i*Np*NLayer, \
			Nonhydrofmod + i*Np*NLayer, InvV3d, Np, NLayer);

		VerticalColumnIntegralField3d(NonhydroHV2d + i*Np2d, Np2d, V2d, TempNonhydroHV2d + i*Np2d, \
			TempNonhydrofield3d + i*Np*NLayer, hv + i*Np*NLayer, Jz + i*Np*NLayer, \
			Nonhydrofmod + i*Np*NLayer, InvV3d, Np, NLayer);
	}

	GetVerticalVelocityAtSurfaceAndBottom(Weta, Wbot, cell, SurfaceBoundaryEdge, BottomBoundaryEdge, InnerEdge, BoundaryEdge, mesh2d, InnerEdge2d, BoundaryEdge2d, \
		cell2d, ftype2d, ftype, PHPX, PHPY, NonhydroHU2d, hu, NonhydroHV2d, hv, z2d, zbot, h2d, h, fext2d, fext, gra, *Hcrit, zx, zy, u, v, w, pept);

	memcpy(Wnew, Wbot, Np2d*K2d*sizeof(double));
	//memcpy(Wnew, Weta, Np2d*K2d * sizeof(double));

	memcpy(Unew, ubot, Np2d*K2d * sizeof(double));

	memcpy(Vnew, vbot, Np2d*K2d * sizeof(double));

	GetFirstOrderPartialDerivativeInVerticalDirection(PUPS, PVPS, PWPS, mesh, cell, BottomEdge, BottomBoundaryEdge, SurfaceBoundaryEdge, u, v, w, Weta, Wbot);

	free(InvV3d);

}

void DotCriticalDevideByLocalValue(double *dest, double *h, int Nfp, double Hcrit){
	for (int i = 0; i < Nfp; i++){
		if (h[i]>=Hcrit){
			dest[i] = dest[i] / h[i];
		}
		else{
			dest[i] = 0.0;
		}
	}
}

void DotCriticalDevideByAveragedValue(double *dest, double *hM, double *hP, int Nfp, double Hcrit){
	for (int i = 0; i < Nfp; i++){
		if (hM[i] >= Hcrit && hP[i] >= Hcrit){
			dest[i] = dest[i] / ((hM[i] + hP[i]) / 2);
		}
	}
}

void GetFirstOrderPartialDerivativeInVerticalDirection(double *PupsDest, double *PvpsDest, double *PwpsDest, const mxArray *mesh, const mxArray *cell, \
	const mxArray *BottomEdge, const mxArray *BottomBoundaryEdge, const mxArray *SurfaceBoundaryEdge, double *u, double *v, double *w, double *ws, \
	double *wb){
	mxArray *Temprx = mxGetField(mesh, 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(mesh, 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(mesh, 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *Tempsy = mxGetField(mesh, 0, "sy");
	double *sy = mxGetPr(Tempsy);
	mxArray *Temptz = mxGetField(mesh, 0, "tz");
	double *tz = mxGetPr(Temptz);
	mxArray *TempJ = mxGetField(mesh, 0, "J");
	double *J = mxGetPr(TempJ);
	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	mxArray *TempDr = mxGetField(cell, 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(cell, 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempInvM = mxGetField(cell, 0, "invM");
	double *invM = mxGetPr(TempInvM);

	/*For bottom edge object*/
	mxArray *TempBotENe = mxGetField(BottomEdge, 0, "Ne");
	int BotENe = (int)mxGetScalar(TempBotENe);
	mxArray *TempBotENfp = mxGetField(BottomEdge, 0, "Nfp");
	int BotENfp = (int)mxGetScalar(TempBotENfp);
	mxArray *TempBotEMb = mxGetField(BottomEdge, 0, "M");
	double *BotEMb = mxGetPr(TempBotEMb);
	mxArray *TempBotEJs = mxGetField(BottomEdge, 0, "Js");
	double *BotEJs = mxGetPr(TempBotEJs);
	mxArray *TempBotEnz = mxGetField(BottomEdge, 0, "nz");
	double *BotEnz = mxGetPr(TempBotEnz);
	mxArray *TempBotEFToE = mxGetField(BottomEdge, 0, "FToE");
	double *BotEFToE = mxGetPr(TempBotEFToE);
	mxArray *TempBotEFToF = mxGetField(BottomEdge, 0, "FToF");
	double *BotEFToF = mxGetPr(TempBotEFToF);
	mxArray *TempBotEFToN1 = mxGetField(BottomEdge, 0, "FToN1");
	double *BotEFToN1 = mxGetPr(TempBotEFToN1);
	mxArray *TempBotEFToN2 = mxGetField(BottomEdge, 0, "FToN2");
	double *BotEFToN2 = mxGetPr(TempBotEFToN2);

	/*For bottom boundary edge object*/
	mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBENfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
	int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	mxArray *TempBotBEMb = mxGetField(BottomBoundaryEdge, 0, "M");
	double *BotBEMb = mxGetPr(TempBotBEMb);
	mxArray *TempBotBEJs = mxGetField(BottomBoundaryEdge, 0, "Js");
	double *BotBEJs = mxGetPr(TempBotBEJs);
	mxArray *TempBotBEnz = mxGetField(BottomBoundaryEdge, 0, "nz");
	double *BotBEnz = mxGetPr(TempBotBEnz);
	mxArray *TempBotBEFToE = mxGetField(BottomBoundaryEdge, 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToF = mxGetField(BottomBoundaryEdge, 0, "FToF");
	double *BotBEFToF = mxGetPr(TempBotBEFToF);
	mxArray *TempBotBEFToN1 = mxGetField(BottomBoundaryEdge, 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);

	/*For surface boundary edge object*/
	mxArray *TempSurfBENe = mxGetField(SurfaceBoundaryEdge, 0, "Ne");
	int SurfBENe = (int)mxGetScalar(TempSurfBENe);
	mxArray *TempSurfBENfp = mxGetField(SurfaceBoundaryEdge, 0, "Nfp");
	int SurfBENfp = (int)mxGetScalar(TempSurfBENfp);
	mxArray *TempSurfBEMb = mxGetField(SurfaceBoundaryEdge, 0, "M");
	double *SurfBEMb = mxGetPr(TempSurfBEMb);
	mxArray *TempSurfBEJs = mxGetField(SurfaceBoundaryEdge, 0, "Js");
	double *SurfBEJs = mxGetPr(TempSurfBEJs);
	mxArray *TempSurfBEnz = mxGetField(SurfaceBoundaryEdge, 0, "nz");
	double *SurfBEnz = mxGetPr(TempSurfBEnz);
	mxArray *TempSurfBEFToE = mxGetField(SurfaceBoundaryEdge, 0, "FToE");
	double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
	mxArray *TempSurfBEFToF = mxGetField(SurfaceBoundaryEdge, 0, "FToF");
	double *SurfBEFToF = mxGetPr(TempSurfBEFToF);
	mxArray *TempSurfBEFToN1 = mxGetField(SurfaceBoundaryEdge, 0, "FToN1");
	double *SurfBEFToN1 = mxGetPr(TempSurfBEFToN1);


	double *uM = NonhydroBotEfm, *vM = NonhydroBotEfm + BotENfp*BotENe, *wM = NonhydroBotEfm + 2 * BotENfp*BotENe;
	double *uP = NonhydroBotEfp, *vP = NonhydroBotEfp + BotENfp*BotENe, *wP = NonhydroBotEfp + 2 * BotENfp*BotENe;
	/*Allocate memory for BotEFluxM, BotEFluxP and BotEFluxS, and calculate these flux term*/


	double *NonhydroUBotEFluxM = NonhydroBotEFluxM, *NonhydroVBotEFluxM = NonhydroBotEFluxM + BotENfp*BotENe, \
		*NonhydroWBotEFluxM = NonhydroBotEFluxM + 2*BotENfp*BotENe;
	double *NonhydroUBotEFluxP = NonhydroBotEFluxP, *NonhydroVBotEFluxP = NonhydroBotEFluxP + BotENfp*BotENe, \
		*NonhydroWBotEFluxP = NonhydroBotEFluxP + 2 * BotENfp*BotENe;
	double *UBotEFluxS = NonhydroBotEFluxS, *VBotEFluxS = NonhydroBotEFluxS + BotENfp*BotENe, \
		*WBotEFluxS = NonhydroBotEFluxS + 2 * BotENfp*BotENe;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++){
		/*Fetch variable BotEfm and BotEfp first*/
		FetchInnerEdgeFacialValue(uM + face*BotENfp, uP + face*BotENfp, u, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		FetchInnerEdgeFacialValue(vM + face*BotENfp, vP + face*BotENfp, v, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);
		FetchInnerEdgeFacialValue(wM + face*BotENfp, wP + face*BotENfp, w, BotEFToE + 2 * face, \
			BotEFToN1 + BotENfp*face, BotEFToN2 + BotENfp*face, Np, BotENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroUBotEFluxM + face*BotENfp, uM + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroUBotEFluxP + face*BotENfp, uP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroVBotEFluxM + face*BotENfp, vM + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroVBotEFluxP + face*BotENfp, vP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroWBotEFluxM + face*BotENfp, wM + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroWBotEFluxP + face*BotENfp, wP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);

		EvaluateNonhydroVerticalFaceNumFlux_Central(UBotEFluxS + face*BotENfp, uM + face*BotENfp, uP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(VBotEFluxS + face*BotENfp, vM + face*BotENfp, vP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(WBotEFluxS + face*BotENfp, wM + face*BotENfp, wP + face*BotENfp, BotEnz + face*BotENfp, BotENfp);
	}
	/*This part is also used by partial derivative in horizontal direction*/

	memset(NonhydroERHS, 0, Np*K * 3 * Nface*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++){
		for (int field = 0; field < 3; field++){
			StrongFormInnerEdgeRHS(face, BotEFToE, BotEFToF, Np, K, BotENfp, BotEFToN1, BotEFToN2, NonhydroBotEFluxM + field*BotENe*BotENfp, \
				NonhydroBotEFluxP + field*BotENe*BotENfp, NonhydroBotEFluxS + field*BotENe*BotENfp, BotEJs, BotEMb, NonhydroERHS + field*Np*K*Nface);
		}
	}

	uM = NonhydroBotBEfm;
	vM = NonhydroBotBEfm + BotBENfp*BotBENe;
	wM = NonhydroBotBEfm + 2 * BotBENfp*BotBENe;


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		FetchBoundaryEdgeFacialValue(uM + face*BotBENfp, u, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(vM + face*BotBENfp, v, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
		FetchBoundaryEdgeFacialValue(wM + face*BotBENfp, w, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroBotBEFluxM + face*BotBENfp, uM + face*BotBENfp, BotBEnz + face*BotBENfp, BotBENfp);
		DotProduct(NonhydroBotBEFluxS + face*BotBENfp, uM + face*BotBENfp, BotBEnz + face*BotBENfp, BotBENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroBotBEFluxM + BotBENe*BotBENfp + face*BotBENfp, vM + face*BotBENfp, BotBEnz + face*BotBENfp, BotBENfp);
		DotProduct(NonhydroBotBEFluxS + BotBENe*BotBENfp + face*BotBENfp, vM + face*BotBENfp, BotBEnz + face*BotBENfp, BotBENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroBotBEFluxM + 2*BotBENe*BotBENfp + face*BotBENfp, wM + face*BotBENfp, BotBEnz + face*BotBENfp, BotBENfp);
		/*Here, the vertical velocity at the bottom boundary is imposed as the numerical flux*/
		DotProduct(NonhydroBotBEFluxS + 2 * BotBENe*BotBENfp + face*BotBENfp, wb + face*BotBENfp, BotBEnz + face*BotBENfp, BotBENfp);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		for (int field = 0; field < 3; field++){
			StrongFormBoundaryEdgeRHS(face, BotBEFToE, BotBEFToF, Np, K, BotBENfp, BotBEFToN1, NonhydroBotBEFluxM + field*BotBENe*BotBENfp, NonhydroBotBEFluxS + field*BotBENe*BotBENfp, BotBEJs, BotBEMb, NonhydroERHS + field*Np*K*Nface);
		}
	}

	uM = NonhydroSurfBEfm;
	vM = NonhydroSurfBEfm + SurfBENfp*SurfBENe;
	wM = NonhydroSurfBEfm + 2 * SurfBENfp*SurfBENe;



#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		FetchBoundaryEdgeFacialValue(uM + face*SurfBENfp, u, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(vM + face*SurfBENfp, v, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
		FetchBoundaryEdgeFacialValue(wM + face*SurfBENfp, w, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroSurfBEFluxM + face*SurfBENfp, uM + face*SurfBENfp, SurfBEnz + face*SurfBENfp, SurfBENfp);
		DotProduct(NonhydroSurfBEFluxS + face*SurfBENfp, uM + face*SurfBENfp, SurfBEnz + face*SurfBENfp, SurfBENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroSurfBEFluxM + SurfBENe*SurfBENfp + face*SurfBENfp, vM + face*SurfBENfp, SurfBEnz + face*SurfBENfp, SurfBENfp);
		DotProduct(NonhydroSurfBEFluxS + SurfBENe*SurfBENfp + face*SurfBENfp, vM + face*SurfBENfp, SurfBEnz + face*SurfBENfp, SurfBENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(NonhydroSurfBEFluxM + 2 * SurfBENe*SurfBENfp + face*SurfBENfp, wM + face*SurfBENfp, SurfBEnz + face*SurfBENfp, SurfBENfp);
		/*Here, the vertical velocity at the bottom boundary is imposed as the numerical flux*/
		DotProduct(NonhydroSurfBEFluxS + 2 * SurfBENe*SurfBENfp + face*SurfBENfp, ws + face*SurfBENfp, SurfBEnz + face*SurfBENfp, SurfBENfp);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		for (int field = 0; field < 3; field++){
			StrongFormBoundaryEdgeRHS(face, SurfBEFToE, SurfBEFToF, Np, K, SurfBENfp, SurfBEFToN1, NonhydroSurfBEFluxM + field*SurfBENe*SurfBENfp, NonhydroSurfBEFluxS + field*SurfBENe*SurfBENfp, SurfBEJs, SurfBEMb, NonhydroERHS + field*Np*K*Nface);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field<3; field++){
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
		for (int field = 0; field < 3; field++){

			MultiEdgeContributionByLiftOperator(NonhydroERHS + field*Np*K*Nface + k*Np, NonhydroTempFacialIntegral + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		/*$\bold{t_z}\cdot (Dt*u)$*/
		GetVolumnIntegral1d(PupsDest + k*Np, &np, &oneI, &np, &one, \
			Dt, &np, u + k*Np, &np, &zero, &np, tz + k*Np);
		/*$\bold{t_z}\cdot (Dt*v)$*/
		GetVolumnIntegral1d(PvpsDest + k*Np, &np, &oneI, &np, &one, \
			Dt, &np, v + k*Np, &np, &zero, &np, tz + k*Np);
		/*$\bold{t_z}\cdot (Dt*w)$*/
		GetVolumnIntegral1d(PwpsDest + k*Np, &np, &oneI, &np, &one, \
			Dt, &np, w + k*Np, &np, &zero, &np, tz + k*Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){

		Minus(PupsDest + k*Np, PupsDest + k*Np, NonhydroERHS + k*Np, Np);

		Minus(PvpsDest + k*Np, PvpsDest + k*Np, NonhydroERHS + Np*K*Nface + k*Np, Np);

		Minus(PwpsDest + k*Np, PwpsDest + k*Np, NonhydroERHS + 2 * Np*K*Nface + k*Np, Np);

	}
}

void GetVerticalVelocityAtSurfaceAndBottom(double *Wetadest, double *Wbotdest, const mxArray *cell, const mxArray *SurfaceBoundaryEdge, \
	const mxArray *BottomBoundaryEdge, const mxArray *InnerEdge, const mxArray *BoundaryEdge, const mxArray *mesh2d, const mxArray *InnerEdge2d,\
	const mxArray *BoundaryEdge2d, const mxArray *cell2d, signed char *ftype2d, signed char *ftype3d, double *PHPX, double *PHPY, double *hu2d, \
	double *hu, double *hv2d, double *hv, double *z2d, double *z3d, double *h2d, double *h3d, double *fext2d, double *fext3d, double gra, \
	double Hcrit, double *zx3d, double *zy3d, double *u, double *v, double *w, double *pept){
	/******************************For the following version, the flux at the boundary is calculated by depth-integration of it's three dimension value**************************************/
	/*
	mxArray *Temprx2d = mxGetField(mesh2d, 0, "rx");
	double *rx2d = mxGetPr(Temprx2d);
	int Np2d = (int)mxGetM(Temprx2d);
	int K2d = (int)mxGetN(Temprx2d);
	mxArray *Tempsx2d = mxGetField(mesh2d, 0, "sx");
	double *sx2d = mxGetPr(Tempsx2d);
	mxArray *Tempry2d = mxGetField(mesh2d, 0, "ry");
	double *ry2d = mxGetPr(Tempry2d);
	mxArray *Tempsy2d = mxGetField(mesh2d, 0, "sy");
	double *sy2d = mxGetPr(Tempsy2d);
	mxArray *TempJ2d = mxGetField(mesh2d, 0, "J");
	double *J2d = mxGetPr(TempJ2d);

	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np3d = (int)mxGetScalar(TempNp);

	//Properties contained in two dimensional inner edge
	mxArray *TempIENe2d = mxGetField(InnerEdge2d, 0, "Ne");
	int IENe2d = (int)mxGetScalar(TempIENe2d);
	mxArray *TempIENfp2d = mxGetField(InnerEdge2d, 0, "Nfp");
	int IENfp2d = (int)mxGetScalar(TempIENfp2d);
	mxArray *TempIEFToE2d = mxGetField(InnerEdge2d, 0, "FToE");
	double *IEFToE2d = mxGetPr(TempIEFToE2d);
	mxArray *TempIEFToF2d = mxGetField(InnerEdge2d, 0, "FToF");
	double *IEFToF2d = mxGetPr(TempIEFToF2d);
	mxArray *TempIEFToN12d = mxGetField(InnerEdge2d, 0, "FToN1");
	double *IEFToN12d = mxGetPr(TempIEFToN12d);
	mxArray *TempIEFToN22d = mxGetField(InnerEdge2d, 0, "FToN2");
	double *IEFToN22d = mxGetPr(TempIEFToN22d);
	mxArray *TempIEnx2d = mxGetField(InnerEdge2d, 0, "nx");
	double *IEnx2d = mxGetPr(TempIEnx2d);
	mxArray *TempIEny2d = mxGetField(InnerEdge2d, 0, "ny");
	double *IEny2d = mxGetPr(TempIEny2d);
	mxArray *TempIEJs2d = mxGetField(InnerEdge2d, 0, "Js");
	double *IEJs2d = mxGetPr(TempIEJs2d);
	mxArray *TempIEMb2d = mxGetField(InnerEdge2d, 0, "M");
	double *IEMb2d = mxGetPr(TempIEMb2d);

	//Properties contained in three dimensional inner edge
	mxArray *TempIENe = mxGetField(InnerEdge, 0, "Ne");
	int IENe3d = (int)mxGetScalar(TempIENe);
	mxArray *TempNLayer = mxGetField(InnerEdge, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);
	mxArray *TempIENfp3d = mxGetField(InnerEdge, 0, "Nfp");
	int IENfp3d = (int)mxGetScalar(TempIENfp3d);
	mxArray *TempIEFToE3d = mxGetField(InnerEdge, 0, "FToE");
	double *IEFToE3d = mxGetPr(TempIEFToE3d);
	mxArray *TempIEFToF3d = mxGetField(InnerEdge, 0, "FToF");
	double *IEFToF3d = mxGetPr(TempIEFToF3d);
	mxArray *TempIEFToN13d = mxGetField(InnerEdge, 0, "FToN1");
	double *IEFToN13d = mxGetPr(TempIEFToN13d);
	mxArray *TempIEFToN23d = mxGetField(InnerEdge, 0, "FToN2");
	double *IEFToN23d = mxGetPr(TempIEFToN23d);
	mxArray *TempIEnx3d = mxGetField(InnerEdge, 0, "nx");
	double *IEnx3d = mxGetPr(TempIEnx3d);
	mxArray *TempIEny3d = mxGetField(InnerEdge, 0, "ny");
	double *IEny3d = mxGetPr(TempIEny3d);
	mxArray *TempIEJz3d = mxGetField(InnerEdge, 0, "Jz");
	double *IEJz3d = mxGetPr(TempIEJz3d);
	mxArray *TempV1d = mxGetField(InnerEdge, 0, "V1d");
	double *V1d = mxGetPr(TempV1d);
	mxArray *TempV2d = mxGetField(InnerEdge, 0, "V2d");
	double *V2d = mxGetPr(TempV2d);
	double *InvV2d = malloc(IENfp3d * IENfp3d * sizeof(double));
	memcpy(InvV2d, V2d, IENfp3d * IENfp3d * sizeof(double));
	MatrixInverse(InvV2d, (ptrdiff_t)IENfp3d);

	//Properties contained in two dimensional boundary edge
	mxArray *TempBENe2d = mxGetField(BoundaryEdge2d, 0, "Ne");
	int BENe2d = (int)mxGetScalar(TempBENe2d);
	mxArray *TempBENfp2d = mxGetField(BoundaryEdge2d, 0, "Nfp");
	int BENfp2d = (int)mxGetScalar(TempBENfp2d);
	mxArray *TempBEFToE2d = mxGetField(BoundaryEdge2d, 0, "FToE");
	double *BEFToE2d = mxGetPr(TempBEFToE2d);
	mxArray *TempBEFToF2d = mxGetField(BoundaryEdge2d, 0, "FToF");
	double *BEFToF2d = mxGetPr(TempBEFToF2d);
	mxArray *TempBEFToN12d = mxGetField(BoundaryEdge2d, 0, "FToN1");
	double *BEFToN12d = mxGetPr(TempBEFToN12d);
	mxArray *TempBEnx2d = mxGetField(BoundaryEdge2d, 0, "nx");
	double *BEnx2d = mxGetPr(TempBEnx2d);
	mxArray *TempBEny2d = mxGetField(BoundaryEdge2d, 0, "ny");
	double *BEny2d = mxGetPr(TempBEny2d);
	mxArray *TempBEJs2d = mxGetField(BoundaryEdge2d, 0, "Js");
	double *BEJs2d = mxGetPr(TempBEJs2d);
	mxArray *TempBEMb2d = mxGetField(BoundaryEdge2d, 0, "M");
	double *BEMb2d = mxGetPr(TempBEMb2d);

	//Properties contained in three dimensional inner edge
	mxArray *TempBENe3d = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe3d = (int)mxGetScalar(TempBENe3d);
	mxArray *TempBENfp3d = mxGetField(BoundaryEdge, 0, "Nfp");
	int BENfp3d = (int)mxGetScalar(TempBENfp3d);
	mxArray *TempBEFToE3d = mxGetField(BoundaryEdge, 0, "FToE");
	double *BEFToE3d = mxGetPr(TempBEFToE3d);
	mxArray *TempBEFToF3d = mxGetField(BoundaryEdge, 0, "FToF");
	double *BEFToF3d = mxGetPr(TempBEFToF3d);
	mxArray *TempBEFToN13d = mxGetField(BoundaryEdge, 0, "FToN1");
	double *BEFToN13d = mxGetPr(TempBEFToN13d);
	mxArray *TempBEnx3d = mxGetField(BoundaryEdge, 0, "nx");
	double *BEnx3d = mxGetPr(TempBEnx3d);
	mxArray *TempBEny3d = mxGetField(BoundaryEdge, 0, "ny");
	double *BEny3d = mxGetPr(TempBEny3d);
	mxArray *TempBEJz3d = mxGetField(BoundaryEdge, 0, "Jz");
	double *BEJz3d = mxGetPr(TempBEJz3d);

	//Data contained in two-dimensional standard cell
	mxArray *TempDr2d = mxGetField(cell2d, 0, "Dr");
	double *Dr2d = mxGetPr(TempDr2d);
	mxArray *TempDs2d = mxGetField(cell2d, 0, "Ds");
	double *Ds2d = mxGetPr(TempDs2d);
	mxArray *TempNface = mxGetField(cell2d, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempinvM2d = mxGetField(cell2d, 0, "invM");
	double *invM2d = mxGetPr(TempinvM2d);

	//For bottom boundary edge object
	mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
	int BotBENe = (int)mxGetScalar(TempBotBENe);
	mxArray *TempBotBENfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
	int BotBENfp = (int)mxGetScalar(TempBotBENfp);
	mxArray *TempBotBEFToE = mxGetField(BottomBoundaryEdge, 0, "FToE");
	double *BotBEFToE = mxGetPr(TempBotBEFToE);
	mxArray *TempBotBEFToN1 = mxGetField(BottomBoundaryEdge, 0, "FToN1");
	double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);
	//For surface boundary edge object
	mxArray *TempSurfBENe = mxGetField(SurfaceBoundaryEdge, 0, "Ne");
	int SurfBENe = (int)mxGetScalar(TempSurfBENe);
	mxArray *TempSurfBENfp = mxGetField(SurfaceBoundaryEdge, 0, "Nfp");
	int SurfBENfp = (int)mxGetScalar(TempSurfBENfp);
	mxArray *TempSurfBEFToE = mxGetField(SurfaceBoundaryEdge, 0, "FToE");
	double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
	mxArray *TempSurfBEFToN1 = mxGetField(SurfaceBoundaryEdge, 0, "FToN1");
	double *SurfBEFToN1 = mxGetPr(TempSurfBEFToN1);

	double *IEhuM2d = NonhydroIEfm2d, *IEhvM2d = NonhydroIEfm2d + IENfp2d * IENe2d, \
		*IEhM2d = NonhydroIEfm2d + 2 * IENfp2d*IENe2d;

	double *IEhuP2d = NonhydroIEfp2d, *IEhvP2d = NonhydroIEfp2d + IENfp2d * IENe2d, \
		*IEhP2d = NonhydroIEfp2d + 2 * IENfp2d*IENe2d;

	memset(NonhydroIEFluxM2d, 0, IENfp2d*IENe2d*sizeof(double));

	memset(NonhydroIEFluxP2d, 0, IENfp2d*IENe2d*sizeof(double));

	memset(NonhydroIEFluxS2d, 0, IENfp2d*IENe2d*sizeof(double));

	memset(NonhydroERHS2d, 0, Np2d*K2d*Nface*sizeof(double));

	ptrdiff_t np = Np2d;
	ptrdiff_t oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		//$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$
		GetVolumnIntegral2d(NonhydroPCEVolumeIntegralX + k*Np2d, NonhydroPCETempVolumeIntegralX + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hu2d + k*Np2d, &np, &zero, &np, rx2d + k*Np2d, sx2d + k*Np2d);
		//$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$
		GetVolumnIntegral2d(NonhydroPCEVolumeIntegralY + k*Np2d, NonhydroPCETempVolumeIntegralY + k*Np2d, &np, &oneI, &np, &one, \
			Dr2d, Ds2d, &np, hv2d + k*Np2d, &np, &zero, &np, ry2d + k*Np2d, sy2d + k*Np2d);

		Add(NonhydroRHS2d + k*Np2d, NonhydroPCEVolumeIntegralX + k*Np2d, NonhydroPCEVolumeIntegralY + k*Np2d, Np2d);
	}
	//Two dimensional inner edge flux part

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		FetchInnerEdgeFacialValue(IEhM2d + e*IENfp2d, IEhP2d + e*IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhuM2d + e*IENfp2d, IEhuP2d + e*IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		FetchInnerEdgeFacialValue(IEhvM2d + e*IENfp2d, IEhvP2d + e*IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
		GetFacialFluxTerm2d(NonhydroIEFluxM2d + e*IENfp2d, IEhuM2d + e*IENfp2d, IEhvM2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
		GetFacialFluxTerm2d(NonhydroIEFluxP2d + e*IENfp2d, IEhuP2d + e*IENfp2d, IEhvP2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
//		GetPCENumericalFluxTerm_HLLC_LAI(NonhydroIEFluxS2d + e*IENfp2d, NonhydroIEfm2d + e*IENfp2d, NonhydroIEfp2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, &gra, Hcrit, IENfp2d, IENe2d);
	}
	double *IEhuM3d = NonhydroIEfm;
	double *IEhuP3d = NonhydroIEfp;
	double *IEhvM3d = NonhydroIEfm + IENe3d * IENfp3d;
	double *IEhvP3d = NonhydroIEfp + IENe3d * IENfp3d;
	double *IEhM3d = NonhydroIEfm + 2 * IENe3d * IENfp3d;
	double *IEhP3d = NonhydroIEfp + 2 * IENe3d * IENfp3d;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		FetchInnerEdgeFacialValue(IEhM3d + e*IENfp3d, IEhP3d + e*IENfp3d, h3d, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		FetchInnerEdgeFacialValue(IEhuM3d + e*IENfp3d, IEhuP3d + e*IENfp3d, hu, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		FetchInnerEdgeFacialValue(IEhvM3d + e*IENfp3d, IEhvP3d + e*IENfp3d, hv, IEFToE3d + 2 * e, IEFToN13d + e*IENfp3d, IEFToN23d + e*IENfp3d, Np3d, IENfp3d);
		GetPCENumericalFluxTerm_HLLC_LAI(NonhydroIEFluxS + e*IENfp3d, NonhydroIEfm + e*IENfp3d, NonhydroIEfp + e*IENfp3d, IEnx3d + e*IENfp3d, IEny3d + e*IENfp3d, &gra, Hcrit, IENfp3d, IENe3d);
	}

	memset(NonhydroIEfmod, 0, IENe2d*IENfp3d*sizeof(double));

	//void VerticalFaceColumnIntegral(double *dest, double *source, double *fmod, double *InvV2d, int Nfp2d, double *Jz, int Nlayer, double *V1d, int LNfp2d, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		VerticalFaceColumnIntegral(NonhydroIEFluxS2d + e*IENfp2d, NonhydroIEFluxS + e*NLayer*IENfp3d, NonhydroIEfmod + e*IENfp3d, InvV2d, (ptrdiff_t)IENfp3d, IEJz3d + e*NLayer*IENfp3d, NLayer, V1d, (ptrdiff_t)IENfp2d, (int)(*(IEFToF3d + e*NLayer * 2)));
	}

	double *BEhuM3d = NonhydroBEfm;
	double *BEhvM3d = NonhydroBEfm + BENe3d * BENfp3d;
	double *BEhM3d = NonhydroBEfm + 2 * BENe3d * BENfp3d;

	int Nfield = 2;
	//The following void pointer is added on 08/25/2021 to accomadate the usage of function ImposeBoundaryCondition
	double *varFieldIndex = NULL;

	//fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype3d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhuM3d + e*BENfp3d, hu, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(BEhvM3d + e*BENfp3d, hv, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(BEhM3d + e*BENfp3d, h3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);
		FetchBoundaryEdgeFacialValue(NonhydrozM + e*BENfp3d, z3d, BEFToE3d + 2 * e, BEFToN13d + e*BENfp3d, Np3d, BENfp3d);

		ImposeBoundaryCondition(&gra, type, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, NonhydroBEfm + e*BENfp3d, NonhydroBEfp + e*BENfp3d, \
			NonhydrozM + e*BENfp3d, NonhydrozP + e*BENfp3d, fext3d + e*BENfp3d, BENfp3d, Nfield, BENe3d, varFieldIndex);
		EvaluateHydroStaticReconstructValue(Hcrit, NonhydroBEfm + e*BENfp3d, NonhydroBEfp + e*BENfp3d, NonhydrozM + e*BENfp3d, NonhydrozP + e*BENfp3d, BENfp3d, Nfield, BENe3d);
		GetPCENumericalFluxTerm_HLLC_LAI(NonhydroBEFluxS + e*BENfp3d, NonhydroBEfm + e*BENfp3d, NonhydroBEfp + e*BENfp3d, BEnx3d + e*BENfp3d, BEny3d + e*BENfp3d, &gra, Hcrit, BENfp3d, BENe3d);
	}

	memset(NonhydroBEfmod, 0, BENe2d*BENfp3d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		VerticalFaceColumnIntegral(NonhydroBEFluxS2d + e*BENfp2d, NonhydroBEFluxS + e*NLayer*BENfp3d, NonhydroBEfmod + e*BENfp3d, InvV2d, (ptrdiff_t)BENfp3d, BEJz3d + e*NLayer*BENfp3d, NLayer, V1d, (ptrdiff_t)BENfp2d, (int)(*(BEFToF3d + e*NLayer * 2)));
	}

	double *BEhuM2d = NonhydroBEfm2d, *BEhvM2d = NonhydroBEfm2d + BENe2d * BENfp2d, \
		*BEhM2d = NonhydroBEfm2d + 2 * BENe2d * BENfp2d;

	memset(NonhydroBEFluxM2d, 0, BENe2d*BENfp2d*sizeof(double));
	//fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
		FetchBoundaryEdgeFacialValue(BEhM2d + e*BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhuM2d + e*BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(BEhvM2d + e*BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		FetchBoundaryEdgeFacialValue(NonhydroBEzM2d + e*BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
		ImposeBoundaryCondition(&gra, type, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, NonhydroBEfm2d + e*BENfp2d, NonhydroBEfp2d + e*BENfp2d, \
			NonhydroBEzM2d + e*BENfp2d, NonhydroBEzP2d + e*BENfp2d, fext2d + e*BENfp2d, BENfp2d, Nfield, BENe2d, varFieldIndex);
		EvaluateHydroStaticReconstructValue(Hcrit, NonhydroBEfm2d + e*BENfp2d, NonhydroBEfp2d + e*BENfp2d, NonhydroBEzM2d + e*BENfp2d, NonhydroBEzP2d + e*BENfp2d, BENfp2d, Nfield, BENe2d);
		GetFacialFluxTerm2d(NonhydroBEFluxM2d + e*BENfp2d, BEhuM2d + e*BENfp2d, BEhvM2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BENfp2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, NonhydroIEFluxM2d, NonhydroIEFluxP2d, NonhydroIEFluxS2d, IEJs2d, IEMb2d, NonhydroERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif    
	for (int e = 0; e < BENe2d; e++){
		StrongFormBoundaryEdgeRHS(e, BEFToE2d, BEFToF2d, Np2d, K2d, BENfp2d, BEFToN12d, NonhydroBEFluxM2d, NonhydroBEFluxS2d, BEJs2d, BEMb2d, NonhydroERHS2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		for (int face = 1; face<Nface; face++){
			Add(NonhydroERHS2d + k*Np2d, NonhydroERHS2d + k*Np2d, NonhydroERHS2d + face*Np2d*K2d + k*Np2d, Np2d);
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		MultiEdgeContributionByLiftOperator(NonhydroERHS2d + k*Np2d, NonhydroPCETempFacialIntegral + k*Np2d, &np, &oneI, &np, \
			&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
	}

	//Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		Minus(NonhydroRHS2d + k*Np2d, NonhydroERHS2d + k*Np2d, NonhydroRHS2d + k*Np2d, Np2d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		FetchBoundaryEdgeFacialValue(ueta + face*SurfBENfp, u, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np3d, SurfBENfp);
		FetchBoundaryEdgeFacialValue(veta + face*SurfBENfp, v, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np3d, SurfBENfp);
		//$\eta_x=\frac{\partial H}{\partial x}$, here $\eta_x$ is a buff for $H_x$
		FetchBoundaryEdgeFacialValue(etax + face*SurfBENfp, PHPX, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np3d, SurfBENfp);
		//$Tempz_x=\frac{\partial z}{\partial x}$
		FetchBoundaryEdgeFacialValue(TempZx2d + face*SurfBENfp, zx3d, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np3d, SurfBENfp);
		//$\eta_x=\H_x + Tempz_x$, since $\eta = H+z$
		Add(etax + face*SurfBENfp, etax + face*SurfBENfp, TempZx2d + face*SurfBENfp, SurfBENfp);
		//$u_{\eta}\eta_x$, at present, the data is stored in \eta_x
		DotProduct(etax + face*SurfBENfp, ueta + face*SurfBENfp, etax + face*SurfBENfp, SurfBENfp);

		FetchBoundaryEdgeFacialValue(etay + face*SurfBENfp, PHPY, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np3d, SurfBENfp);
		FetchBoundaryEdgeFacialValue(TempZy2d + face*SurfBENfp, zy3d, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np3d, SurfBENfp);
		Add(etay + face*SurfBENfp, etay + face*SurfBENfp, TempZy2d + face*SurfBENfp, SurfBENfp);
		DotProduct(etay + face*SurfBENfp, veta + face*SurfBENfp, etay + face*SurfBENfp, SurfBENfp);
		//$w_{\eta} = \frac{\partial \eta}{\partial t} + u_{\eta}\frac{\partial \eta}{\partial x} + v_{\eta}\frac{\partial \eta}{\partial y}$
		Add(Wetadest + face*SurfBENfp, NonhydroRHS2d + face*SurfBENfp, etax + face*SurfBENfp, SurfBENfp);
		Add(Wetadest + face*SurfBENfp, Wetadest + face*SurfBENfp, etay + face*SurfBENfp, SurfBENfp);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		FetchBoundaryEdgeFacialValue(ubot + face*BotBENfp, u, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np3d, BotBENfp);

		FetchBoundaryEdgeFacialValue(vbot + face*BotBENfp, v, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np3d, BotBENfp);
		//$Tempz_x=\frac{\partial z}{\partial x}$
		FetchBoundaryEdgeFacialValue(TempZx2d + face*BotBENfp, zx3d, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np3d, BotBENfp);
		//$u_dz_x$, at present, the data is stored in TempZx2d
		DotProduct(TempZx2d + face*BotBENfp, ubot + face*BotBENfp, TempZx2d + face*BotBENfp, BotBENfp);

		FetchBoundaryEdgeFacialValue(TempZy2d + face*BotBENfp, zy3d, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np3d, BotBENfp);
		//$v_dz_y$, at present, the data is stored in TempZy2d
		DotProduct(TempZy2d + face*BotBENfp, vbot + face*BotBENfp, TempZy2d + face*BotBENfp, BotBENfp);

		//$w_{bot} = \frac{\partial z}{\partial t} + u_d\frac{\partial z}{\partial x} + v_d\frac{\partial z}{\partial y}$
		Add(Wbotdest + face*SurfBENfp, TempZx2d + face*BotBENfp, TempZy2d + face*BotBENfp, BotBENfp);
	}

	free(InvV2d);
	*/
	/******************************For the following version, the flux at the boundary is calculated by the Riemann solver**************************************/
/*
mxArray *Temprx2d = mxGetField(mesh2d, 0, "rx");
double *rx2d = mxGetPr(Temprx2d);
int Np2d = (int)mxGetM(Temprx2d);
int K2d = (int)mxGetN(Temprx2d);
mxArray *Tempsx2d = mxGetField(mesh2d, 0, "sx");
double *sx2d = mxGetPr(Tempsx2d);
mxArray *Tempry2d = mxGetField(mesh2d, 0, "ry");
double *ry2d = mxGetPr(Tempry2d);
mxArray *Tempsy2d = mxGetField(mesh2d, 0, "sy");
double *sy2d = mxGetPr(Tempsy2d);
mxArray *TempJ2d = mxGetField(mesh2d, 0, "J");
double *J2d = mxGetPr(TempJ2d);
*/
mxArray *TempNp = mxGetField(cell, 0, "Np");
int Np = (int)mxGetScalar(TempNp);
/*
//Properties contained in two dimensional inner edge
mxArray *TempIENe2d = mxGetField(InnerEdge2d, 0, "Ne");
int IENe2d = (int)mxGetScalar(TempIENe2d);
mxArray *TempIENfp2d = mxGetField(InnerEdge2d, 0, "Nfp");
int IENfp2d = (int)mxGetScalar(TempIENfp2d);
mxArray *TempIEFToE2d = mxGetField(InnerEdge2d, 0, "FToE");
double *IEFToE2d = mxGetPr(TempIEFToE2d);
mxArray *TempIEFToF2d = mxGetField(InnerEdge2d, 0, "FToF");
double *IEFToF2d = mxGetPr(TempIEFToF2d);
mxArray *TempIEFToN12d = mxGetField(InnerEdge2d, 0, "FToN1");
double *IEFToN12d = mxGetPr(TempIEFToN12d);
mxArray *TempIEFToN22d = mxGetField(InnerEdge2d, 0, "FToN2");
double *IEFToN22d = mxGetPr(TempIEFToN22d);
mxArray *TempIEnx2d = mxGetField(InnerEdge2d, 0, "nx");
double *IEnx2d = mxGetPr(TempIEnx2d);
mxArray *TempIEny2d = mxGetField(InnerEdge2d, 0, "ny");
double *IEny2d = mxGetPr(TempIEny2d);
mxArray *TempIEJs2d = mxGetField(InnerEdge2d, 0, "Js");
double *IEJs2d = mxGetPr(TempIEJs2d);
mxArray *TempIEMb2d = mxGetField(InnerEdge2d, 0, "M");
double *IEMb2d = mxGetPr(TempIEMb2d);

//Properties contained in two dimensional boundary edge
mxArray *TempBENe2d = mxGetField(BoundaryEdge2d, 0, "Ne");
int BENe2d = (int)mxGetScalar(TempBENe2d);
mxArray *TempBENfp2d = mxGetField(BoundaryEdge2d, 0, "Nfp");
int BENfp2d = (int)mxGetScalar(TempBENfp2d);
mxArray *TempBEFToE2d = mxGetField(BoundaryEdge2d, 0, "FToE");
double *BEFToE2d = mxGetPr(TempBEFToE2d);
mxArray *TempBEFToF2d = mxGetField(BoundaryEdge2d, 0, "FToF");
double *BEFToF2d = mxGetPr(TempBEFToF2d);
mxArray *TempBEFToN12d = mxGetField(BoundaryEdge2d, 0, "FToN1");
double *BEFToN12d = mxGetPr(TempBEFToN12d);
mxArray *TempBEnx2d = mxGetField(BoundaryEdge2d, 0, "nx");
double *BEnx2d = mxGetPr(TempBEnx2d);
mxArray *TempBEny2d = mxGetField(BoundaryEdge2d, 0, "ny");
double *BEny2d = mxGetPr(TempBEny2d);
mxArray *TempBEJs2d = mxGetField(BoundaryEdge2d, 0, "Js");
double *BEJs2d = mxGetPr(TempBEJs2d);
mxArray *TempBEMb2d = mxGetField(BoundaryEdge2d, 0, "M");
double *BEMb2d = mxGetPr(TempBEMb2d);

//Data contained in two-dimensional standard cell
mxArray *TempDr2d = mxGetField(cell2d, 0, "Dr");
double *Dr2d = mxGetPr(TempDr2d);
mxArray *TempDs2d = mxGetField(cell2d, 0, "Ds");
double *Ds2d = mxGetPr(TempDs2d);
mxArray *TempNface = mxGetField(cell2d, 0, "Nface");
int Nface = (int)mxGetScalar(TempNface);
mxArray *TempinvM2d = mxGetField(cell2d, 0, "invM");
double *invM2d = mxGetPr(TempinvM2d);
*/
//For bottom boundary edge object
mxArray *TempBotBENe = mxGetField(BottomBoundaryEdge, 0, "Ne");
int BotBENe = (int)mxGetScalar(TempBotBENe);
mxArray *TempBotBENfp = mxGetField(BottomBoundaryEdge, 0, "Nfp");
int BotBENfp = (int)mxGetScalar(TempBotBENfp);
mxArray *TempBotBEFToE = mxGetField(BottomBoundaryEdge, 0, "FToE");
double *BotBEFToE = mxGetPr(TempBotBEFToE);
mxArray *TempBotBEFToN1 = mxGetField(BottomBoundaryEdge, 0, "FToN1");
double *BotBEFToN1 = mxGetPr(TempBotBEFToN1);
//For surface boundary edge object
mxArray *TempSurfBENe = mxGetField(SurfaceBoundaryEdge, 0, "Ne");
int SurfBENe = (int)mxGetScalar(TempSurfBENe);
mxArray *TempSurfBENfp = mxGetField(SurfaceBoundaryEdge, 0, "Nfp");
int SurfBENfp = (int)mxGetScalar(TempSurfBENfp);
mxArray *TempSurfBEFToE = mxGetField(SurfaceBoundaryEdge, 0, "FToE");
double *SurfBEFToE = mxGetPr(TempSurfBEFToE);
mxArray *TempSurfBEFToN1 = mxGetField(SurfaceBoundaryEdge, 0, "FToN1");
double *SurfBEFToN1 = mxGetPr(TempSurfBEFToN1);
/*
double *IEhuM2d = NonhydroIEfm2d, *IEhvM2d = NonhydroIEfm2d + IENfp2d * IENe2d, \
*IEhM2d = NonhydroIEfm2d + 2 * IENfp2d*IENe2d;

double *IEhuP2d = NonhydroIEfp2d, *IEhvP2d = NonhydroIEfp2d + IENfp2d * IENe2d, \
*IEhP2d = NonhydroIEfp2d + 2 * IENfp2d*IENe2d;

memset(NonhydroIEFluxM2d, 0, IENfp2d*IENe2d * sizeof(double));

memset(NonhydroIEFluxP2d, 0, IENfp2d*IENe2d * sizeof(double));

memset(NonhydroIEFluxS2d, 0, IENfp2d*IENe2d * sizeof(double));

memset(NonhydroERHS2d, 0, Np2d*K2d*Nface * sizeof(double));

ptrdiff_t np = Np2d;
ptrdiff_t oneI = 1;
double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int k = 0; k < K2d; k++) {
	//$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$
	GetVolumnIntegral2d(NonhydroPCEVolumeIntegralX + k*Np2d, NonhydroPCETempVolumeIntegralX + k*Np2d, &np, &oneI, &np, &one, \
		Dr2d, Ds2d, &np, hu2d + k*Np2d, &np, &zero, &np, rx2d + k*Np2d, sx2d + k*Np2d);
	//$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$
	GetVolumnIntegral2d(NonhydroPCEVolumeIntegralY + k*Np2d, NonhydroPCETempVolumeIntegralY + k*Np2d, &np, &oneI, &np, &one, \
		Dr2d, Ds2d, &np, hv2d + k*Np2d, &np, &zero, &np, ry2d + k*Np2d, sy2d + k*Np2d);

	Add(NonhydroRHS2d + k*Np2d, NonhydroPCEVolumeIntegralX + k*Np2d, NonhydroPCEVolumeIntegralY + k*Np2d, Np2d);
}
//Two dimensional inner edge flux part

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int e = 0; e < IENe2d; e++) {
	FetchInnerEdgeFacialValue(IEhM2d + e*IENfp2d, IEhP2d + e*IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
	FetchInnerEdgeFacialValue(IEhuM2d + e*IENfp2d, IEhuP2d + e*IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
	FetchInnerEdgeFacialValue(IEhvM2d + e*IENfp2d, IEhvP2d + e*IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e*IENfp2d, IEFToN22d + e*IENfp2d, Np2d, IENfp2d);
	GetFacialFluxTerm2d(NonhydroIEFluxM2d + e*IENfp2d, IEhuM2d + e*IENfp2d, IEhvM2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
	GetFacialFluxTerm2d(NonhydroIEFluxP2d + e*IENfp2d, IEhuP2d + e*IENfp2d, IEhvP2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, IENfp2d);
	GetPCENumericalFluxTerm_HLLC_LAI(NonhydroIEFluxS2d + e*IENfp2d, NonhydroIEfm2d + e*IENfp2d, NonhydroIEfp2d + e*IENfp2d, IEnx2d + e*IENfp2d, IEny2d + e*IENfp2d, &gra, Hcrit, IENfp2d, IENe2d);
}

double *BEhuM2d = NonhydroBEfm2d, *BEhvM2d = NonhydroBEfm2d + BENe2d * BENfp2d, \
*BEhM2d = NonhydroBEfm2d + 2 * BENe2d * BENfp2d;

memset(NonhydroBEFluxS2d, 0, BENe2d*BENfp2d * sizeof(double));

memset(NonhydroBEFluxM2d, 0, BENe2d*BENfp2d * sizeof(double));
int Nfield = 2;
double *varFieldIndex = NULL;
//fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int e = 0; e < BENe2d; e++) {
	NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
	FetchBoundaryEdgeFacialValue(BEhM2d + e*BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
	FetchBoundaryEdgeFacialValue(BEhuM2d + e*BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
	FetchBoundaryEdgeFacialValue(BEhvM2d + e*BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
	FetchBoundaryEdgeFacialValue(NonhydroBEzM2d + e*BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e*BENfp2d, Np2d, BENfp2d);
	ImposeBoundaryCondition(&gra, type, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, NonhydroBEfm2d + e*BENfp2d, NonhydroBEfp2d + e*BENfp2d, \
		NonhydroBEzM2d + e*BENfp2d, NonhydroBEzP2d + e*BENfp2d, fext2d + e*BENfp2d, BENfp2d, Nfield, BENe2d, varFieldIndex);
	EvaluateHydroStaticReconstructValue(Hcrit, NonhydroBEfm2d + e*BENfp2d, NonhydroBEfp2d + e*BENfp2d, NonhydroBEzM2d + e*BENfp2d, NonhydroBEzP2d + e*BENfp2d, BENfp2d, Nfield, BENe2d);
	GetFacialFluxTerm2d(NonhydroBEFluxM2d + e*BENfp2d, BEhuM2d + e*BENfp2d, BEhvM2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, BENfp2d);
	GetPCENumericalFluxTerm_HLLC_LAI(NonhydroBEFluxS2d + e*BENfp2d, NonhydroBEfm2d + e*BENfp2d, NonhydroBEfp2d + e*BENfp2d, BEnx2d + e*BENfp2d, BEny2d + e*BENfp2d, &gra, Hcrit, BENfp2d, BENe2d);
}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int e = 0; e < IENe2d; e++) {
	StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, NonhydroIEFluxM2d, NonhydroIEFluxP2d, NonhydroIEFluxS2d, IEJs2d, IEMb2d, NonhydroERHS2d);
}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif    
for (int e = 0; e < BENe2d; e++) {
	StrongFormBoundaryEdgeRHS(e, BEFToE2d, BEFToF2d, Np2d, K2d, BENfp2d, BEFToN12d, NonhydroBEFluxM2d, NonhydroBEFluxS2d, BEJs2d, BEMb2d, NonhydroERHS2d);
}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int k = 0; k < K2d; k++) {
	for (int face = 1; face<Nface; face++) {
		Add(NonhydroERHS2d + k*Np2d, NonhydroERHS2d + k*Np2d, NonhydroERHS2d + face*Np2d*K2d + k*Np2d, Np2d);
	}
}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int k = 0; k < K2d; k++) {
	MultiEdgeContributionByLiftOperator(NonhydroERHS2d + k*Np2d, NonhydroPCETempFacialIntegral + k*Np2d, &np, &oneI, &np, \
		&one, invM2d, &np, &np, &zero, &np, J2d + k*Np2d, Np2d);
}

//Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int k = 0; k < K2d; k++) {
	Minus(NonhydroRHS2d + k*Np2d, NonhydroERHS2d + k*Np2d, NonhydroRHS2d + k*Np2d, Np2d);
}
*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int face = 0; face < SurfBENe; face++) {
	FetchBoundaryEdgeFacialValue(ueta + face*SurfBENfp, u, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
	FetchBoundaryEdgeFacialValue(veta + face*SurfBENfp, v, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
	//$\eta_x=\frac{\partial H}{\partial x}$, here $\eta_x$ is a buff for $H_x$
	FetchBoundaryEdgeFacialValue(etax + face*SurfBENfp, PHPX, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
	//$Tempz_x=\frac{\partial z}{\partial x}$
	FetchBoundaryEdgeFacialValue(TempZx2d + face*SurfBENfp, zx3d, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
	//$\eta_x=\H_x + Tempz_x$, since $\eta = H+z$
	Add(etax + face*SurfBENfp, etax + face*SurfBENfp, TempZx2d + face*SurfBENfp, SurfBENfp);
	//$u_{\eta}\eta_x$, at present, the data is stored in \eta_x
	DotProduct(etax + face*SurfBENfp, ueta + face*SurfBENfp, etax + face*SurfBENfp, SurfBENfp);

	FetchBoundaryEdgeFacialValue(etay + face*SurfBENfp, PHPY, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
	FetchBoundaryEdgeFacialValue(TempZy2d + face*SurfBENfp, zy3d, SurfBEFToE + 2 * face, SurfBEFToN1 + face*SurfBENfp, Np, SurfBENfp);
	Add(etay + face*SurfBENfp, etay + face*SurfBENfp, TempZy2d + face*SurfBENfp, SurfBENfp);
	DotProduct(etay + face*SurfBENfp, veta + face*SurfBENfp, etay + face*SurfBENfp, SurfBENfp);
	//$w_{\eta} = \frac{\partial \eta}{\partial t} + u_{\eta}\frac{\partial \eta}{\partial x} + v_{\eta}\frac{\partial \eta}{\partial y}$
	//Add(Wetadest + face*SurfBENfp, NonhydroRHS2d + face*SurfBENfp, etax + face*SurfBENfp, SurfBENfp);
	Add(Wetadest + face*SurfBENfp, pept + face*SurfBENfp, etax + face*SurfBENfp, SurfBENfp);
	Add(Wetadest + face*SurfBENfp, Wetadest + face*SurfBENfp, etay + face*SurfBENfp, SurfBENfp);
}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
for (int face = 0; face < BotBENe; face++) {
	FetchBoundaryEdgeFacialValue(ubot + face*BotBENfp, u, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);

	FetchBoundaryEdgeFacialValue(vbot + face*BotBENfp, v, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
	//$Tempz_x=\frac{\partial z}{\partial x}$
	FetchBoundaryEdgeFacialValue(TempZx2d + face*BotBENfp, zx3d, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
	//$u_dz_x$, at present, the data is stored in TempZx2d
	DotProduct(TempZx2d + face*BotBENfp, ubot + face*BotBENfp, TempZx2d + face*BotBENfp, BotBENfp);

	FetchBoundaryEdgeFacialValue(TempZy2d + face*BotBENfp, zy3d, BotBEFToE + 2 * face, BotBEFToN1 + face*BotBENfp, Np, BotBENfp);
	//$v_dz_y$, at present, the data is stored in TempZy2d
	DotProduct(TempZy2d + face*BotBENfp, vbot + face*BotBENfp, TempZy2d + face*BotBENfp, BotBENfp);

	//$w_{bot} = \frac{\partial z}{\partial t} + u_d\frac{\partial z}{\partial x} + v_d\frac{\partial z}{\partial y}$
	Add(Wbotdest + face*SurfBENfp, TempZx2d + face*BotBENfp, TempZy2d + face*BotBENfp, BotBENfp);
}

}

void GetFirstOrderPartialDerivativeInHorizontalDirection( double *PHPX, double *PHPY, double *PUPX, \
	double *PUPY, double *PVPX, double *PVPY, const mxArray *mesh, const mxArray *cell, const mxArray *InnerEdge, \
	const mxArray *BoundaryEdge, double *h, double *u, double *v, double *hu, double *hv, double *z, \
	double gra, double *fext, double Hcrit, signed char *ftype){

	mxArray *Temprx = mxGetField(mesh, 0, "rx");
	double *rx = mxGetPr(Temprx);
	mxArray *Tempsx = mxGetField(mesh, 0, "sx");
	double *sx = mxGetPr(Tempsx);
	mxArray *Tempry = mxGetField(mesh, 0, "ry");
	double *ry = mxGetPr(Tempry);
	mxArray *Tempsy = mxGetField(mesh, 0, "sy");
	double *sy = mxGetPr(Tempsy);
	mxArray *Temptz = mxGetField(mesh, 0, "tz");
	double *tz = mxGetPr(Temptz);
	mxArray *TempJ = mxGetField(mesh, 0, "J");
	double *J = mxGetPr(TempJ);
	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	mxArray *TempDr = mxGetField(cell, 0, "Dr");
	double *Dr = mxGetPr(TempDr);
	mxArray *TempDs = mxGetField(cell, 0, "Ds");
	double *Ds = mxGetPr(TempDs);
	mxArray *TempDt = mxGetField(cell, 0, "Dt");
	double *Dt = mxGetPr(TempDt);
	mxArray *TempNface = mxGetField(cell, 0, "Nface");
	int Nface = (int)mxGetScalar(TempNface);
	mxArray *TempInvM = mxGetField(cell, 0, "invM");
	double *invM = mxGetPr(TempInvM);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);

	mxArray *TempIENe = mxGetField(InnerEdge, 0, "Ne");
	int IENe = (int)mxGetScalar(TempIENe);
	mxArray *TempIENfp = mxGetField(InnerEdge, 0, "Nfp");
	int IENfp = (int)mxGetScalar(TempIENfp);
	mxArray *TempIEMb = mxGetField(InnerEdge, 0, "M");
	double *IEMb = mxGetPr(TempIEMb);
	mxArray *TempIEJs = mxGetField(InnerEdge, 0, "Js");
	double *IEJs = mxGetPr(TempIEJs);
	mxArray *TempIEnx = mxGetField(InnerEdge, 0, "nx");
	double *IEnx = mxGetPr(TempIEnx);
	mxArray *TempIEny = mxGetField(InnerEdge, 0, "ny");
	double *IEny = mxGetPr(TempIEny);
	mxArray *TempIELAV = mxGetField(InnerEdge, 0, "LAV");
	double *IELAV = mxGetPr(TempIELAV);
	mxArray *TempIEFToE = mxGetField(InnerEdge, 0, "FToE");
	double *IEFToE = mxGetPr(TempIEFToE);
	mxArray *TempIEFToF = mxGetField(InnerEdge, 0, "FToF");
	double *IEFToF = mxGetPr(TempIEFToF);
	mxArray *TempIEFToN1 = mxGetField(InnerEdge, 0, "FToN1");
	double *IEFToN1 = mxGetPr(TempIEFToN1);
	mxArray *TempIEFToN2 = mxGetField(InnerEdge, 0, "FToN2");
	double *IEFToN2 = mxGetPr(TempIEFToN2);

	mxArray *TempBENe = mxGetField(BoundaryEdge, 0, "Ne");
	int BENe = (int)mxGetScalar(TempBENe);
	mxArray *TempBENfp = mxGetField(BoundaryEdge, 0, "Nfp");
	int BENfp = mxGetScalar(TempBENfp);
	mxArray *TempBEMb = mxGetField(BoundaryEdge, 0, "M");
	double *BEMb = mxGetPr(TempBEMb);
	mxArray *TempBEJs = mxGetField(BoundaryEdge, 0, "Js");
	double *BEJs = mxGetPr(TempBEJs);
	mxArray *TempBEnx = mxGetField(BoundaryEdge, 0, "nx");
	double *BEnx = mxGetPr(TempBEnx);
	mxArray *TempBEny = mxGetField(BoundaryEdge, 0, "ny");
	double *BEny = mxGetPr(TempBEny);
	mxArray *TempBELAV = mxGetField(BoundaryEdge, 0, "LAV");
	double *BELAV = mxGetPr(TempBELAV);
	mxArray *TempBEFToE = mxGetField(BoundaryEdge, 0, "FToE");
	double *BEFToE = mxGetPr(TempBEFToE);
	mxArray *TempBEFToF = mxGetField(BoundaryEdge, 0, "FToF");
	double *BEFToF = mxGetPr(TempBEFToF);
	mxArray *TempBEFToN1 = mxGetField(BoundaryEdge, 0, "FToN1");
	double *BEFToN1 = mxGetPr(TempBEFToN1);



	double *uM = NonhydroIEfm, *vM = NonhydroIEfm + IENe*IENfp, \
		*hM = NonhydroIEfm + 2 * IENe*IENfp;
	double *uP = NonhydroIEfp, *vP = NonhydroIEfp + IENe*IENfp, \
		*hP = NonhydroIEfp + 2 * IENe*IENfp;

	double *UIEfluxM = NonhydroIEFluxM, *VIEfluxM = NonhydroIEFluxM + IENe*IENfp, \
		*HIEfluxMx = NonhydroIEFluxM + 2 * IENe*IENfp, *HIEfluxMy = NonhydroIEFluxM + 3 * IENe*IENfp, \
		*UIEfluxMy = NonhydroIEFluxM + 4 * IENe*IENfp, *VIEfluxMx = NonhydroIEFluxM + 5 * IENe*IENfp;

	double *UIEfluxP = NonhydroIEFluxP, *VIEfluxP = NonhydroIEFluxP + IENe*IENfp, \
		*HIEfluxPx = NonhydroIEFluxP + 2 * IENe*IENfp, *HIEfluxPy = NonhydroIEFluxP + 3 * IENe*IENfp, \
		*UIEfluxPy = NonhydroIEFluxP + 4 * IENe*IENfp, *VIEfluxPx = NonhydroIEFluxP + 5 * IENe*IENfp;

	double *UIEfluxS = NonhydroIEFluxS, *VIEfluxS = NonhydroIEFluxS + IENe*IENfp, \
		*HIEfluxSx = NonhydroIEFluxS + 2 * IENe*IENfp, *HIEfluxSy = NonhydroIEFluxS + 3 * IENe*IENfp, \
		*UIEfluxSy = NonhydroIEFluxS + 4 * IENe*IENfp, *VIEfluxSx = NonhydroIEFluxS + 5 * IENe*IENfp;
    /*The following null pointer is added on 08/25/2021 to accomdate the usage of function ImposeBoundaryCondition,
     defined in SWE.c. For this part, this parameter is actually useless*/
    double *varFieldIndex = NULL;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		/*Fetch variable IEfm and IEfp first*/
		FetchInnerEdgeFacialValue(hM + face*IENfp, hP + face*IENfp, h, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		FetchInnerEdgeFacialValue(uM + face*IENfp, uP + face*IENfp, u, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);
		FetchInnerEdgeFacialValue(vM + face*IENfp, vP + face*IENfp, v, IEFToE + 2 * face, \
			IEFToN1 + IENfp*face, IEFToN2 + IENfp*face, Np, IENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(UIEfluxM + face*IENfp, uM + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(UIEfluxP + face*IENfp, uP + face*IENfp, IEnx + face*IENfp, IENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(UIEfluxMy + face*IENfp, uM + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(UIEfluxPy + face*IENfp, uP + face*IENfp, IEny + face*IENfp, IENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(VIEfluxM + face*IENfp, vM + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(VIEfluxP + face*IENfp, vP + face*IENfp, IEny + face*IENfp, IENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(VIEfluxMx + face*IENfp, vM + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(VIEfluxPx + face*IENfp, vP + face*IENfp, IEnx + face*IENfp, IENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(HIEfluxMx + face*IENfp, hM + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(HIEfluxPx + face*IENfp, hP + face*IENfp, IEnx + face*IENfp, IENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(HIEfluxMy + face*IENfp, hM + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceSurfFlux(HIEfluxPy + face*IENfp, hP + face*IENfp, IEny + face*IENfp, IENfp);

		EvaluateNonhydroVerticalFaceNumFlux_Central(UIEfluxS + face*IENfp, uM + face*IENfp, uP + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(VIEfluxS + face*IENfp, vM + face*IENfp, vP + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(UIEfluxSy + face*IENfp, uM + face*IENfp, uP + face*IENfp, IEny + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(VIEfluxSx + face*IENfp, vM + face*IENfp, vP + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(HIEfluxSx + face*IENfp, hM + face*IENfp, hP + face*IENfp, IEnx + face*IENfp, IENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(HIEfluxSy + face*IENfp, hM + face*IENfp, hP + face*IENfp, IEny + face*IENfp, IENfp);
	}

	memset(NonhydroERHS, 0, Np * K * 6 * (Nface - 2) * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		for (int field = 0; field < 6; field++){
			StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, NonhydroIEFluxM + field*IENe*IENfp, \
				NonhydroIEFluxP + field*IENe*IENfp, NonhydroIEFluxS + field*IENe*IENfp, IEJs, IEMb, NonhydroERHS + field*Np*K*(Nface - 2));
		}
	}

	/**************************************Boundary Edge Part*******************************************************/
	/********************************************************************/

	/*Allocate memory for fm and fp defined over boundary edges. Here, variables correspond to hu, hv, hw, h */
	uM = NonhydroBEfm, vM = NonhydroBEfm + BENfp*BENe, hM = NonhydroBEfm + 2 * BENfp*BENe;
    uP = NonhydroBEfm, vP = NonhydroBEfm + BENfp*BENe, hP = NonhydroBEfm + 2 * BENfp*BENe;
//	uP = NonhydroBEfp, vP = NonhydroBEfp + BENfp*BENe, hP = NonhydroBEfp + 2 * BENfp*BENe;

//	double *Hext = fext + 2 * BENfp*BENe;

	double *UBEfluxM = NonhydroBEFluxM, *VBEfluxM = NonhydroBEFluxM + BENe*BENfp, \
		*HBEfluxMx = NonhydroBEFluxM + 2 * BENe*BENfp, *HBEfluxMy = NonhydroBEFluxM + 3 * BENe*BENfp, \
		*UBEfluxMy = NonhydroBEFluxM + 4 * BENe*BENfp, *VBEfluxMx = NonhydroBEFluxM + 5 * BENe*BENfp;

	double *UBEfluxS = NonhydroBEFluxS, *VBEfluxS = NonhydroBEFluxS + BENe*BENfp, \
		*HBEfluxSx = NonhydroBEFluxS + 2 * BENe*BENfp, *HBEfluxSy = NonhydroBEFluxS + 3 * BENe*BENfp, \
		*UBEfluxSy = NonhydroBEFluxS + 4 * BENe*BENfp, *VBEfluxSx = NonhydroBEFluxS + 5 * BENe*BENfp;

	/*Fetch variable BEfm and BEfp first, then impose boundary condition and conduct hydrostatic reconstruction.
	Finally, calculate local flux term, adjacent flux term and numerical flux term*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
		FetchBoundaryEdgeFacialValue(uM + face*BENfp, hu, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(vM + face*BENfp, hv, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(hM + face*BENfp, h, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);
		FetchBoundaryEdgeFacialValue(NonhydrozM + face*BENfp, z, BEFToE + 2 * face, BEFToN1 + face*BENfp, Np, BENfp);

		ImposeBoundaryCondition(&gra, type, BEnx + face*BENfp, BEny + face*BENfp, NonhydroBEfm + face*BENfp, NonhydroBEfp + face*BENfp, \
			NonhydrozM + face*BENfp, NonhydrozP + face*BENfp, fext + face*BENfp, BENfp, 2, BENe, varFieldIndex);

		EvaluateHydroStaticReconstructValue(Hcrit, NonhydroBEfm + face*BENfp, NonhydroBEfp + face*BENfp, NonhydrozM + face*BENfp, NonhydrozP + face*BENfp, BENfp, 2, BENe);

		
		if (type == NdgEdgeSlipWall) {
			for (int i = 0; i < BENfp; i++) {
				hP[i + face*BENfp] = hM[i + face*BENfp];
				uP[i + face*BENfp] = uM[i + face*BENfp];
				vP[i + face*BENfp] = vM[i + face*BENfp];
			}
		}
	
        if (type == NdgEdgeClampedVel || type == NdgEdgeClampedDepth) {
			for (int i = 0; i < BENfp; i++){
				hP[i + face*BENfp] = hM[i + face*BENfp];
				uP[i + face*BENfp] = uM[i + face*BENfp];
				vP[i + face*BENfp] = vM[i + face*BENfp];
			}
		}
		
		DotCriticalDivide(uM + face*BENfp, uM + face*BENfp, &Hcrit, hM + face*BENfp, BENfp);
		DotCriticalDivide(uP + face*BENfp, uP + face*BENfp, &Hcrit, hP + face*BENfp, BENfp);
//      MultiplyByConstant(uP + face*BENfp, uP + face*BENfp, 0.95, BENfp);
//		DotCriticalDivide(uP + face*BENfp, uP + face*BENfp, &Hcrit, Hext + face*BENfp, BENfp);
		DotCriticalDivide(vM + face*BENfp, vM + face*BENfp, &Hcrit, hM + face*BENfp, BENfp);
		DotCriticalDivide(vP + face*BENfp, vP + face*BENfp, &Hcrit, hP + face*BENfp, BENfp);
//      MultiplyByConstant(vP + face*BENfp, vP + face*BENfp, 0.95, BENfp);
//		DotCriticalDivide(vP + face*BENfp, vP + face*BENfp, &Hcrit, Hext + face*BENfp, BENfp);
/*		
       if (type == NdgEdgeClampedVel) {
			MultiplyByConstant(uP + face*BENfp, uP + face*BENfp, 1.05, BENfp);
			MultiplyByConstant(vP + face*BENfp, vP + face*BENfp, 1.05, BENfp);
		}
	*/	
		EvaluateNonhydroVerticalFaceSurfFlux(UBEfluxM + face*BENfp, uM + face*BENfp, BEnx + face*BENfp, BENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(UBEfluxMy + face*BENfp, uM + face*BENfp, BEny + face*BENfp, BENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(VBEfluxM + face*BENfp, vM + face*BENfp, BEny + face*BENfp, BENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(VBEfluxMx + face*BENfp, vM + face*BENfp, BEnx + face*BENfp, BENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(HBEfluxMx + face*BENfp, hM + face*BENfp, BEnx + face*BENfp, BENfp);

		EvaluateNonhydroVerticalFaceSurfFlux(HBEfluxMy + face*BENfp, hM + face*BENfp, BEny + face*BENfp, BENfp);

		EvaluateNonhydroVerticalFaceNumFlux_Central(UBEfluxS + face*BENfp, uM + face*BENfp, uP + face*BENfp, BEnx + face*BENfp, BENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(VBEfluxS + face*BENfp, vM + face*BENfp, vP + face*BENfp, BEny + face*BENfp, BENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(UBEfluxSy + face*BENfp, uM + face*BENfp, uP + face*BENfp, BEny + face*BENfp, BENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(VBEfluxSx + face*BENfp, vM + face*BENfp, vP + face*BENfp, BEnx + face*BENfp, BENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(HBEfluxSx + face*BENfp, hM + face*BENfp, hP + face*BENfp, BEnx + face*BENfp, BENfp);
		EvaluateNonhydroVerticalFaceNumFlux_Central(HBEfluxSy + face*BENfp, hM + face*BENfp, hP + face*BENfp, BEny + face*BENfp, BENfp);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		for (int field = 0; field < 6; field++){
			StrongFormBoundaryEdgeRHS(face, BEFToE, BEFToF, Np, K, BENfp, BEFToN1, NonhydroBEFluxM + field*BENe*BENfp, NonhydroBEFluxS + field*BENe*BENfp, BEJs, BEMb, NonhydroERHS + field*Np*K*(Nface-2));
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		for (int field = 0; field<6; field++){
			for (int face = 1; face<Nface - 2; face++){
				Add(NonhydroERHS + field*Np*K*(Nface - 2) + k*Np, NonhydroERHS + field*Np*K*(Nface - 2) + k*Np, NonhydroERHS + field*Np*K*(Nface - 2) + face*Np*K + k*Np, Np);
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
		for (int field = 0; field < 6; field++){
			MultiEdgeContributionByLiftOperator(NonhydroERHS + field*Np*K*(Nface - 2) + k*Np, NonhydroTempFacialIntegral + k*Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k*Np, Np);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
		GetVolumnIntegral2d(PUPX + k*Np, NonhydroTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, u + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);

		GetVolumnIntegral2d(PUPY + k*Np, NonhydroTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, u + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		/*$\bold{r_y}\cdot (Dr*hv)+\bold{s_y}\cdot (Ds*hv)$*/
		GetVolumnIntegral2d(PVPY + k*Np, NonhydroTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, v + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);

		GetVolumnIntegral2d(PVPX + k*Np, NonhydroTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, v + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);

		GetVolumnIntegral2d(PHPX + k*Np, NonhydroTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, h + k*Np, &np, &zero, &np, rx + k*Np, sx + k*Np);

		GetVolumnIntegral2d(PHPY + k*Np, NonhydroTempVolumeIntegral + k*Np, &np, &oneI, &np, &one, \
			Dr, Ds, &np, h + k*Np, &np, &zero, &np, ry + k*Np, sy + k*Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){

		Minus(PUPX + k*Np, PUPX + k*Np, NonhydroERHS + k*Np, Np);

		Minus(PVPY + k*Np, PVPY + k*Np, NonhydroERHS + Np*K*(Nface-2) + k*Np, Np);

		Minus(PHPX + k*Np, PHPX + k*Np, NonhydroERHS + 2 * Np*K*(Nface - 2) + k*Np, Np);

		Minus(PHPY + k*Np, PHPY + k*Np, NonhydroERHS + 3 * Np*K*(Nface - 2) + k*Np, Np);

		Minus(PUPY + k*Np, PUPY + k*Np, NonhydroERHS + 4 * Np*K*(Nface - 2) + k*Np, Np);

		Minus(PVPX + k*Np, PVPX + k*Np, NonhydroERHS + 5 * Np*K*(Nface - 2) + k*Np, Np);

	}

}