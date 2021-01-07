# include "NdgMemory.h"

/*This is for GOTM part*/
double *VCV = NULL, *tkeGOTM = NULL, *epsGOTM = NULL, *LGOTM = NULL, *nuhGOTM = NULL,\
*numGOTM = NULL, *layerHeight = NULL, *huCentralDate = NULL, *hvCentralDate = NULL,\
*huVerticalLine = NULL, *hvVerticalLine = NULL, *shearFrequencyDate = NULL, *buoyanceFrequencyDate = NULL,\
*BottomFrictionLength = NULL, *BottomFrictionVelocity = NULL, *SurfaceFrictionLength = NULL,\
*SurfaceFrictionVelocity = NULL, *eddyViscosityDate = NULL;

signed char *GOTMInitialized = "False";

void GotmSolverMemoryAllocation(int Num2d, int Interface, int Np2d, int K2d, int K3d){
	tkeGOTM = malloc(sizeof(double)*(Num2d*Interface)); 
	epsGOTM = malloc(sizeof(double)*(Num2d*Interface)); 
	LGOTM = malloc(sizeof(double)*(Num2d*Interface));
	nuhGOTM = malloc(sizeof(double)*(Num2d*Interface)); 
	numGOTM = malloc(sizeof(double)*(Num2d*Interface));
	layerHeight = malloc(sizeof(double)*(Num2d*Interface));
	huCentralDate = malloc(sizeof(double)*(Np2d*K3d));
	hvCentralDate = malloc(sizeof(double)*(Np2d*K3d));
	huVerticalLine = malloc(sizeof(double)*(Num2d*Interface));
	hvVerticalLine = malloc(sizeof(double)*(Num2d*Interface));
	shearFrequencyDate = malloc(sizeof(double)*(Num2d*Interface));
	buoyanceFrequencyDate = malloc(sizeof(double)*(Num2d*Interface));
	BottomFrictionLength = malloc(sizeof(double)*Num2d);
	BottomFrictionVelocity = malloc(sizeof(double)*Num2d);
	SurfaceFrictionLength = malloc(sizeof(double)*Num2d);
	SurfaceFrictionVelocity = malloc(sizeof(double)*Num2d);
	eddyViscosityDate = malloc(sizeof(double)*(Num2d * Interface));
	GOTMInitialized = "True";
}

void GotmSolverMemoryDeAllocation(){
	free(nuhGOTM); nuhGOTM = NULL;
	free(numGOTM); numGOTM = NULL;
	free(tkeGOTM); tkeGOTM = NULL;
	free(epsGOTM); epsGOTM = NULL;
	free(LGOTM); LGOTM = NULL;
	free(layerHeight); layerHeight = NULL;
	free(huCentralDate); huCentralDate = NULL;
	free(hvCentralDate); hvCentralDate = NULL;
	free(huVerticalLine); huVerticalLine = NULL;
	free(hvVerticalLine); hvVerticalLine = NULL;
	free(shearFrequencyDate); shearFrequencyDate = NULL;
	free(buoyanceFrequencyDate); buoyanceFrequencyDate = NULL;
	free(BottomFrictionLength); BottomFrictionLength = NULL;
	free(BottomFrictionVelocity); BottomFrictionVelocity = NULL;
	free(SurfaceFrictionLength); SurfaceFrictionLength = NULL;
	free(SurfaceFrictionVelocity); SurfaceFrictionVelocity = NULL;
	free(eddyViscosityDate); eddyViscosityDate = NULL;
	GOTMInitialized = "False";
}


/*This is for vertical velocity solver part*/
double *VSrhs2d = NULL, *VSIEfm2d = NULL, *VSIEfp2d = NULL, *VSIEFluxM2d = NULL, \
*VSIEFluxP2d = NULL, *VSIEFluxS2d = NULL, *VSERHS2d = NULL, *VSVolumeIntegralX = NULL, \
*VSTempVolumeIntegralX = NULL, *VSVolumeIntegralY = NULL, *VSTempVolumeIntegralY = NULL, \
*VSBEfm2d = NULL, *VSBEzM2d = NULL, *VSBEfp2d = NULL, *VSBEzP2d = NULL, *VSBEFluxS2d = NULL, \
*VSBEFluxM2d = NULL, *VSTempFacialIntegral = NULL, *VSfield2d = NULL, *VSrhs3d = NULL, \
*VSIEfm3d = NULL, *VSIEfp3d = NULL, *VSIEFluxM3d = NULL, *VSIEFluxP3d = NULL, *VSIEFluxS3d = NULL, \
*VSERHS3d = NULL, *VSVolumeIntegralX3d = NULL, *VSTempVolumeIntegralX3d = NULL, \
*VSVolumeIntegralY3d = NULL, *VSTempVolumeIntegralY3d = NULL, *VSBEfm3d = NULL, \
*VSBEzM3d = NULL, *VSBEfp3d = NULL, *VSBEzP3d = NULL, *VSBEFluxS3d = NULL, *VSBEFluxM3d = NULL, \
*VSTempFacialIntegral3d = NULL, *VSTempVerticalVelocity = NULL, *VSBotVertVelocity = NULL;

signed char *VertVelocityInitialized = "False";

void VertVelocitySolverMemoryAllocation(int Np2d, int K2d, int IENfp2d, int IENe2d, int Nface, int BENe2d, int BENfp2d, int Np3d,\
	int K3d, int IENfp3d, int IENe3d, int BENe3d, int BENfp3d){
	VSrhs2d = malloc(Np2d*K2d*sizeof(double));
	VSIEfm2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	VSIEfp2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	VSIEFluxM2d = malloc(IENfp2d*IENe2d*sizeof(double));
	VSIEFluxP2d = malloc(IENfp2d*IENe2d*sizeof(double));
	VSIEFluxS2d = malloc(IENfp2d*IENe2d*sizeof(double));
	VSERHS2d = malloc(Np2d*K2d*Nface*sizeof(double));
	VSVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	VSTempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	VSVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	VSTempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	VSBEfm2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	VSBEzM2d = malloc(BENe2d * BENfp2d*sizeof(double));
	VSBEfp2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	VSBEzP2d = malloc(BENe2d * BENfp2d*sizeof(double));
	VSBEFluxS2d = malloc(BENe2d*BENfp2d*sizeof(double));
	VSBEFluxM2d = malloc(BENe2d*BENfp2d*sizeof(double));
	VSTempFacialIntegral = malloc(Np2d*K2d*sizeof(double));
	VSfield2d = malloc(Np3d*K3d*sizeof(double));
	VSrhs3d = malloc(Np3d*K3d*sizeof(double));
	VSIEfm3d = malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	VSIEfp3d = malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	VSIEFluxM3d = malloc(IENfp3d*IENe3d*sizeof(double));
	VSIEFluxP3d = malloc(IENfp3d*IENe3d*sizeof(double));
	VSIEFluxS3d = malloc(IENfp3d*IENe3d*sizeof(double));
	VSERHS3d = malloc(Np3d*K3d*Nface*sizeof(double));
	VSVolumeIntegralX3d = malloc(Np3d*K3d*sizeof(double));
	VSTempVolumeIntegralX3d = malloc(Np3d*K3d*sizeof(double));
	VSVolumeIntegralY3d = malloc(Np3d*K3d*sizeof(double));
	VSTempVolumeIntegralY3d = malloc(Np3d*K3d*sizeof(double));
	VSBEfm3d = malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	VSBEzM3d = malloc(BENe3d * BENfp3d * sizeof(double));
	VSBEfp3d = malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	VSBEzP3d = malloc(BENe3d * BENfp3d * sizeof(double));
	VSBEFluxS3d = malloc(BENe3d*BENfp3d*sizeof(double));
	VSBEFluxM3d = malloc(BENe3d*BENfp3d*sizeof(double));
	VSTempVerticalVelocity = malloc(Np3d*K3d*sizeof(double));
	VSBotVertVelocity = malloc(Np3d*K2d*sizeof(double));
	VertVelocityInitialized = "True";
}

void VertVelocitySolverMemoryDeAllocation(){
	free(VSrhs2d), VSrhs2d = NULL;
	free(VSIEfm2d), VSIEfm2d = NULL; 
	free(VSIEfp2d), VSIEfp2d = NULL;
	free(VSIEFluxM2d), VSIEFluxM2d = NULL;
	free(VSIEFluxP2d), VSIEFluxP2d = NULL;
	free(VSIEFluxS2d), VSIEFluxS2d = NULL;
	free(VSERHS2d), VSERHS2d = NULL;
	free(VSVolumeIntegralX), VSVolumeIntegralX = NULL;
	free(VSTempVolumeIntegralX), VSTempVolumeIntegralX = NULL;
	free(VSVolumeIntegralY), VSVolumeIntegralY = NULL;
	free(VSTempVolumeIntegralY), VSTempVolumeIntegralY = NULL;
	free(VSBEfm2d), VSBEfm2d = NULL;
	free(VSBEzM2d), VSBEzM2d = NULL;
	free(VSBEfp2d), VSBEfp2d = NULL;
	free(VSBEzP2d), VSBEzP2d = NULL;
	free(VSBEFluxS2d), VSBEFluxS2d = NULL;
	free(VSBEFluxM2d), VSBEFluxM2d = NULL;
	free(VSTempFacialIntegral), VSTempFacialIntegral = NULL;
	free(VSfield2d), VSfield2d = NULL;
	free(VSrhs3d), VSrhs3d = NULL;
	free(VSIEfm3d), VSIEfm3d = NULL;
	free(VSIEfp3d), VSIEfp3d = NULL;
	free(VSIEFluxM3d), VSIEFluxM3d = NULL;
	free(VSIEFluxP3d), VSIEFluxP3d = NULL;
	free(VSIEFluxS3d), VSIEFluxS3d = NULL;
	free(VSERHS3d), VSERHS3d = NULL;
	free(VSVolumeIntegralX3d), VSVolumeIntegralX3d = NULL;
	free(VSTempVolumeIntegralX3d), VSTempVolumeIntegralX3d = NULL;
	free(VSVolumeIntegralY3d), VSVolumeIntegralY3d = NULL;
	free(VSTempVolumeIntegralY3d), VSTempVolumeIntegralY3d = NULL;
	free(VSBEfm3d), VSBEfm3d = NULL;
	free(VSBEzM3d), VSBEzM3d = NULL;
	free(VSBEfp3d), VSBEfp3d = NULL;
	free(VSBEzP3d), VSBEzP3d = NULL;
	free(VSBEFluxS3d), VSBEFluxS3d = NULL;
	free(VSBEFluxM3d), VSBEFluxM3d = NULL;
	free(VSTempVerticalVelocity), VSTempVerticalVelocity = NULL;
	free(VSBotVertVelocity), VSBotVertVelocity = NULL;
	VertVelocityInitialized = "False";
}

/*This is for PCE Solver part*/
double *IEfm2d = NULL, *IEfp2d = NULL, *IEFluxM2d = NULL, *IEFluxP2d = NULL, *IEFluxS2d = NULL, \
*ERHS2d = NULL, *PCEVolumeIntegralX = NULL, *PCETempVolumeIntegralX = NULL, *PCEVolumeIntegralY = NULL, \
*PCETempVolumeIntegralY = NULL, *BEfm2d = NULL, *BEzM2d = NULL, *BEfp2d = NULL, *BEzP2d = NULL, \
*BEFluxS2d = NULL, *BEFluxM2d = NULL, *PCETempFacialIntegral = NULL;

signed char *PCEInitialized = "False";

void PCEMemoryAllocation(int IENfp2d, int IENe2d, int Np2d, int K2d, int Nface, int BENe2d, int BENfp2d){
	IEfm2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	IEfp2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	IEFluxM2d = malloc(IENfp2d*IENe2d*sizeof(double));
	IEFluxP2d = malloc(IENfp2d*IENe2d*sizeof(double));
	IEFluxS2d = malloc(IENfp2d*IENe2d*sizeof(double));
	ERHS2d = malloc(Np2d*K2d*Nface*sizeof(double));
	PCEVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	PCETempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	PCEVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	PCETempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	BEfm2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	BEzM2d = malloc(BENe2d * BENfp2d*sizeof(double));
	BEfp2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	BEzP2d = malloc(BENe2d * BENfp2d*sizeof(double));
	BEFluxS2d = malloc(BENe2d*BENfp2d*sizeof(double));
	BEFluxM2d = malloc(BENe2d*BENfp2d*sizeof(double));
	PCETempFacialIntegral = malloc(Np2d*K2d*sizeof(double));
	PCEInitialized = "True";
}

void PCEMemoryDeAllocation(){
	free(IEfm2d), IEfm2d = NULL;
	free(IEfp2d), IEfp2d = NULL;
	free(IEFluxM2d), IEFluxM2d = NULL;
	free(IEFluxP2d), IEFluxP2d = NULL;
	free(IEFluxS2d), IEFluxS2d = NULL;
	free(ERHS2d), ERHS2d = NULL;
	free(PCEVolumeIntegralX), PCEVolumeIntegralX = NULL;
	free(PCETempVolumeIntegralX), PCETempVolumeIntegralX = NULL;
	free(PCEVolumeIntegralY), PCEVolumeIntegralY = NULL;
	free(PCETempVolumeIntegralY), PCETempVolumeIntegralY = NULL;
	free(BEfm2d), BEfm2d = NULL;
	free(BEzM2d), BEzM2d = NULL;
	free(BEfp2d), BEfp2d = NULL;
	free(BEzP2d), BEzP2d = NULL;
	free(BEFluxS2d), BEFluxS2d = NULL;
	free(BEFluxM2d), BEFluxM2d = NULL;
	free(PCETempFacialIntegral), PCETempFacialIntegral = NULL;
	PCEInitialized = "False";
}


/*This is for advection memory part*/
double *TempFacialIntegral = NULL, *IEfm = NULL, *IEfp = NULL, *IEFluxM = NULL, *IEFluxP = NULL, \
*IEFluxS = NULL, *ERHS = NULL, *BEfm = NULL, *BEfp = NULL, *AdvzM = NULL, *AdvzP = NULL, *BEFluxM = NULL, \
*BEFluxS = NULL, *BotEfm = NULL, *BotEfp = NULL, *BotEFluxM = NULL, *BotEFluxP = NULL, *BotEFluxS = NULL, \
*BotBEfm = NULL, *BotBEFluxM = NULL, *BotBEFluxS = NULL, *SurfBEfm = NULL, *SurfBEFluxM = NULL, \
*SurfBEFluxS = NULL, *E = NULL, *G = NULL, *H = NULL, *TempVolumeIntegral = NULL;

signed char *AdvInitialized = "False";

void AdvMemoryAllocation(int Np, int K, int Nvar, int IENfp, int IENe, int Nface, int BENfp, int BENe, int BotENfp,\
                        int BotENe, int BotBENfp, int BotBENe, int SurfBENfp, int SurfBENe){
	TempFacialIntegral = malloc(Np*K*Nvar*sizeof(double));
	IEfm = malloc(IENfp*IENe*(Nvar + 1)*sizeof(double));
	IEfp = malloc(IENfp*IENe*(Nvar + 1)*sizeof(double));
	IEFluxM = malloc(IENfp*IENe*Nvar*sizeof(double));
	IEFluxP = malloc(IENfp*IENe*Nvar*sizeof(double));
	IEFluxS = malloc(IENfp*IENe*Nvar*sizeof(double));
	ERHS = malloc(Np*K*Nvar*Nface*sizeof(double));
	BEfm = malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
	BEfp = malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
	AdvzM = malloc(BENfp*BENe*sizeof(double));
	AdvzP = malloc(BENfp*BENe*sizeof(double));
	BEFluxM = malloc(BENfp*BENe*Nvar*sizeof(double));
	BEFluxS = malloc(BENfp*BENe*Nvar*sizeof(double));
	BotEfm = malloc(BotENfp*BotENe*(Nvar + 2)*sizeof(double));
	BotEfp = malloc(BotENfp*BotENe*(Nvar + 2)*sizeof(double));
	BotEFluxM = malloc(BotENfp*BotENe*Nvar*sizeof(double));
	BotEFluxP = malloc(BotENfp*BotENe*Nvar*sizeof(double));
	BotEFluxS = malloc(BotENfp*BotENe*Nvar*sizeof(double));
	BotBEfm = malloc(BotBENfp*BotBENe*(Nvar + 2)*sizeof(double));
	BotBEFluxM = malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
	BotBEFluxS = malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
	SurfBEfm = malloc(SurfBENfp*SurfBENe*(Nvar + 2)*sizeof(double));
	SurfBEFluxM = malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
	SurfBEFluxS = malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
	/*Allocate memory for E, G and H, and calculate these volume flux term*/
	E = malloc(Np*K*Nvar*sizeof(double));
	G = malloc(Np*K*Nvar*sizeof(double));
	H = malloc(Np*K*Nvar*sizeof(double));
	TempVolumeIntegral = malloc(Np*K*sizeof(double));
	AdvInitialized = "True";
}

void AdvMemoryDeAllocation(){
	free(TempFacialIntegral), TempFacialIntegral = NULL;
	free(IEfm), IEfm = NULL;
	free(IEfp), IEfp = NULL;
	free(IEFluxM), IEFluxM = NULL;
	free(IEFluxP), IEFluxP = NULL;
	free(IEFluxS), IEFluxS = NULL;
	free(ERHS), ERHS = NULL;
	free(BEfm), BEfm = NULL;
	free(BEfp), BEfp = NULL;
	free(AdvzM), AdvzM = NULL;
	free(AdvzP), AdvzP = NULL;
	free(BEFluxM), BEFluxM = NULL;
	free(BEFluxS), BEFluxS = NULL;
	free(BotEfm), BotEfm = NULL;
	free(BotEfp), BotEfp = NULL;
	free(BotEFluxM), BotEFluxM = NULL;
	free(BotEFluxP), BotEFluxP = NULL;
	free(BotEFluxS), BotEFluxS = NULL;
	free(BotBEfm), BotBEfm = NULL;
	free(BotBEFluxM), BotBEFluxM = NULL;
	free(BotBEFluxS), BotBEFluxS = NULL;
	free(SurfBEfm), BotBEFluxS = NULL;
	free(SurfBEFluxM), SurfBEFluxM = NULL;
	free(SurfBEFluxS), SurfBEFluxS = NULL;
	free(E), E = NULL;
	free(G), G = NULL;
	free(H), H = NULL;
	free(TempVolumeIntegral), TempVolumeIntegral = NULL;
	AdvInitialized = "False";
}


/*This is for horizontal diffusion memory part*/
double *HorDiffnv = NULL, *HorDiffvariable = NULL, *HorDiffBEfp = NULL, *HorDiffzM = NULL, \
*HorDiffzP = NULL, *HorDiffTempBEfp = NULL, *HorDiffTempBEfm = NULL, *HorDiffAVx = NULL, \
*HorDiffAVy = NULL, *HorDiffVx = NULL, *HorDiffTempVx = NULL, *HorDiffVy = NULL, *HorDiffTempVy = NULL, \
*HorDiffIEfm = NULL, *HorDiffAVIEfm = NULL, *HorDiffIEfp = NULL, *HorDiffAVIEfp = NULL, *HorDiffIEFluxM = NULL, \
*HorDiffIEFluxP = NULL, *HorDiffBEfm = NULL, *HorDiffIEFluxS = NULL, *HorDiffAVBEfm = NULL, *HorDiffBEFluxM = NULL, \
*HorDiffBEFluxS = NULL, *HorDiffERHSX = NULL, *HorDiffERHSY = NULL, *HorDiffLocalPrimitiveDiffTermX = NULL, \
*HorDiffLocalPrimitiveDiffTermY = NULL, *HorDiffLPDTIEfm = NULL, *HorDiffLPDTIEfp = NULL, *HorDiffLPDTBEfm = NULL, \
*HorDiffTempFacialIntegralX = NULL, *HorDiffTempFacialIntegralY = NULL, *HorDiffInnerEdgeTau = NULL, \
*HorDiffBoundaryEdgeTau = NULL, *HorDiffIEnvfm = NULL, *HorDiffIEnvfp = NULL, *HorDiffBEnvfm = NULL;

char *HorDiffInitialized = "False";

void HorizDiffMemoryAllocation(NdgMeshType type, int Np, int K, int Nvar, int tempNface, int BENfp, int BENe, int IENfp, int IENe){
	int Nfield;
	int Nface;
	if (type == Two){
		/*For 2d shallow water problem, no horizontal diffusion terms are included in the governing equation for water depth $H$*/
		Nfield = Nvar - 1;
		/*For 2d shallow water problem, the face number is equal to TempNface, since there is no surface edge and bottom edge*/
		Nface = tempNface;
	}
	else if (type == Three){
		Nfield = Nvar;
		/*For 3d shallow water problem, the face number is equal to TempNface - 2, since there
		* is surface edge and bottom edge is not considered for horizontal diffusion term*/
		Nface = tempNface - 2;
	}
	HorDiffzM = malloc(BENfp*BENe*sizeof(double));

	HorDiffzP = malloc(BENfp*BENe*sizeof(double));

	HorDiffnv = malloc(Np*K*sizeof(double));

	HorDiffTempBEfp = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));

	HorDiffTempBEfm = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
	/*Allocate memory for the original variable over boundary edge*/
	HorDiffBEfp = malloc(BENfp*BENe*Nfield*sizeof(double));
	/*Allocate memory for the local face value at boundary edge*/
	HorDiffBEfm = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the original variable $u,v$ and $\theta$*/
	HorDiffvariable = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for auxiallary variable $q_x=\frac{\partial u(v,\theta)}{\partial x}$*/
	HorDiffAVx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for auxiallary variable $q_y=\frac{\partial u(v,\theta)}{\partial y}$*/
	HorDiffAVy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for part of the volumn integral of auxiallary variable $q_{x1}=\bold{r_x}\cdot (D_r*u(v,\theta))$,
	and for part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffVx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the rest part of the volumn integral of auxiallary variable $q_{x2}=\bold{r_x}\cdot (D_r*u(v,\theta))$,
	and for the rest part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffTempVx = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for part of the volumn integral of auxiallary variable $q_{y1}=\bold{r_y}\cdot (D_r*u(v,\theta))$,
	and for part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffVy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the rest part of the volumn integral of auxiallary variable $q_{y2}=\bold{s_y}\cdot (D_s*u(v,\theta))$,
	and for the rest part of the volumn integral part of the second order operator $\frac{\partial q}{\partial y}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$*/
	HorDiffTempVy = malloc(Np*K*Nfield*sizeof(double));
	/*Allocate memory for the local face value at inner edge*/
	HorDiffIEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at inner edge*/
	HorDiffAVIEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent face value at inner edge*/
	HorDiffIEfp = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent face value of the auxialary variable at inner edge*/
	HorDiffAVIEfp = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local flux term at inner edge*/
	HorDiffIEFluxM = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent flux term at inner edge*/
	HorDiffIEFluxP = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at inner edge*/
	HorDiffIEFluxS = malloc(IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at boundary edge*/
	HorDiffAVBEfm = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the local flux term at boundary edge*/
	HorDiffBEFluxM = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at boundary edge*/
	HorDiffBEFluxS = malloc(BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for interior edge contribution to right hand side in x direction*/
	HorDiffERHSX = malloc(Np*K*Nfield*Nface*sizeof(double));
	/*Allocate memory for interior edge contribution to right hand side in y direction*/
	HorDiffERHSY = malloc(Np*K*Nfield*Nface*sizeof(double));

	/*Allocate memory for local primitive diffusion term in both x and y direction, $\nv\Nabla u(v,\theta)$*/
	/*This part is used because we need to extract the trace of the local derivative operator to compute the numerical flux*/
	HorDiffLocalPrimitiveDiffTermX = malloc(Np*K*Nfield*sizeof(double));
	HorDiffLocalPrimitiveDiffTermY = malloc(Np*K*Nfield*sizeof(double));
	HorDiffLPDTIEfm = malloc(IENfp*IENe*Nfield*sizeof(double));
	HorDiffLPDTIEfp = malloc(IENfp*IENe*Nfield*sizeof(double));
	HorDiffLPDTBEfm = malloc(BENfp*BENe*Nfield*sizeof(double));

	HorDiffTempFacialIntegralX = malloc(Np*K*Nfield*sizeof(double));
	HorDiffTempFacialIntegralY = malloc(Np*K*Nfield*sizeof(double));

	HorDiffInnerEdgeTau = malloc(IENe*IENfp*sizeof(double));
	HorDiffBoundaryEdgeTau = malloc(BENe*BENfp*sizeof(double));
	HorDiffIEnvfm = malloc(IENe*IENfp*sizeof(double));
	HorDiffIEnvfp = malloc(IENe*IENfp*sizeof(double));
	HorDiffBEnvfm = malloc(BENe*BENfp*sizeof(double));

	HorDiffInitialized = "True";


}

void HorizDiffMemoryDeAllocation()
{
	free(HorDiffnv), HorDiffnv = NULL;
	free(HorDiffvariable), HorDiffvariable = NULL;
	free(HorDiffBEfp), HorDiffBEfp = NULL;
	free(HorDiffzM), HorDiffzM = NULL;
	free(HorDiffzP), HorDiffzP = NULL;
	free(HorDiffTempBEfp), HorDiffTempBEfp = NULL;
	free(HorDiffTempBEfm), HorDiffTempBEfm = NULL;
	free(HorDiffAVx), HorDiffAVx = NULL;
	free(HorDiffAVy), HorDiffAVy = NULL;
	free(HorDiffVx), HorDiffVx = NULL;
	free(HorDiffTempVx), HorDiffVx = NULL;
	free(HorDiffVy), HorDiffVy = NULL;
	free(HorDiffTempVy), HorDiffTempVy = NULL;
	free(HorDiffIEfm), HorDiffIEfm = NULL;
	free(HorDiffAVIEfm), HorDiffAVIEfm = NULL;
	free(HorDiffIEfp), HorDiffIEfp = NULL;
	free(HorDiffAVIEfp), HorDiffAVIEfp = NULL;
	free(HorDiffIEFluxM), HorDiffIEFluxM = NULL;
	free(HorDiffIEFluxP), HorDiffIEFluxP = NULL;
	free(HorDiffBEfm), HorDiffBEfm = NULL;
	free(HorDiffIEFluxS), HorDiffIEFluxS = NULL;
	free(HorDiffAVBEfm), HorDiffAVBEfm = NULL;
	free(HorDiffBEFluxM), HorDiffBEFluxM = NULL;
	free(HorDiffBEFluxS), HorDiffBEFluxS = NULL;
	free(HorDiffERHSX), HorDiffERHSX = NULL;
	free(HorDiffERHSY), HorDiffERHSY = NULL;
	free(HorDiffLocalPrimitiveDiffTermX), HorDiffLocalPrimitiveDiffTermX = NULL;
	free(HorDiffLocalPrimitiveDiffTermY), HorDiffLocalPrimitiveDiffTermY = NULL;
	free(HorDiffLPDTIEfm), HorDiffLPDTIEfm = NULL;
	free(HorDiffLPDTIEfp), HorDiffLPDTIEfp = NULL;
	free(HorDiffLPDTBEfm), HorDiffLPDTBEfm = NULL;
	free(HorDiffTempFacialIntegralX), HorDiffTempFacialIntegralX = NULL;
	free(HorDiffTempFacialIntegralY), HorDiffTempFacialIntegralY = NULL;
	free(HorDiffInnerEdgeTau), HorDiffInnerEdgeTau = NULL;
	free(HorDiffBoundaryEdgeTau), HorDiffBoundaryEdgeTau = NULL;
	free(HorDiffIEnvfm), HorDiffIEnvfm = NULL;
	free(HorDiffIEnvfp), HorDiffIEnvfp = NULL;
	free(HorDiffBEnvfm), HorDiffBEnvfm = NULL;

	HorDiffInitialized = "False";
}