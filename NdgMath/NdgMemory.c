# include "NdgMemory.h"

void MemoryAllocationCheck(double *dest, int size){
    while(dest == NULL){
        dest = malloc(size); 
    }
}

/*This is for three-dimensional non-hydrostatic model*/

char *SWENonhydro3dInitialized = "False";

double *NonhydroHU2d = NULL, *NonhydroHV2d = NULL, \
*TempNonhydroHU2d = NULL, *TempNonhydroHV2d = NULL, *TempNonhydrofield3d = NULL, *Nonhydrofmod = NULL, \
*NonhydroVariable = NULL, *PUPX = NULL, *PVPY = NULL, *PHPX = NULL, *PHPY = NULL, *PSPX = NULL, \
*PSPY = NULL, *TempPSPX = NULL, *TempPSPY = NULL, *NonhydroIEfm = NULL, *NonhydroIEfp = NULL, \
*NonhydroIEFluxM = NULL, *NonhydroIEFluxP = NULL, *NonhydroIEFluxS = NULL, *NonhydroERHS = NULL, \
*NonhydroBEfm = NULL, *NonhydroBEfp = NULL, *NonhydroBEFluxM = NULL, *NonhydroBEFluxS = NULL, \
*NonhydroTempFacialIntegral = NULL, *NonhydroTempVolumeIntegral = NULL, *NonhydrozM = NULL, *NonhydrozP = NULL, \
*CoePS = NULL, *NonhydroBotEfm = NULL, *NonhydroBotEfp = NULL, \
*NonhydroBotEFluxM = NULL, *NonhydroBotEFluxP = NULL, *NonhydroBotEFluxS = NULL, \
*NonhydroBotBEfm = NULL, *NonhydroBotBEFluxM = NULL, *NonhydroBotBEFluxS = NULL, *NonhydroSurfBEfm = NULL, \
*NonhydroSurfBEFluxM = NULL, *NonhydroSurfBEFluxS = NULL, *NonhydroRHS2d = NULL, \
*Weta = NULL, *Wbot = NULL, *ueta = NULL, *ubot = NULL, *veta = NULL, *vbot = NULL, \
*etax = NULL, *etay = NULL, *TempZx2d = NULL, *TempZy2d = NULL, *NonhydroIEfm2d = NULL, \
*NonhydroIEfp2d = NULL, *NonhydroIEFluxM2d = NULL, *NonhydroIEFluxP2d = NULL, \
*NonhydroIEFluxS2d = NULL, *NonhydroERHS2d = NULL, *NonhydroPCEVolumeIntegralX = NULL, \
*NonhydroPCETempVolumeIntegralX = NULL, *NonhydroPCEVolumeIntegralY = NULL, \
*NonhydroPCETempVolumeIntegralY = NULL, *NonhydroBEfm2d = NULL, *NonhydroBEzM2d = NULL, \
*NonhydroBEfp2d = NULL, *NonhydroBEzP2d = NULL, *NonhydroBEFluxS2d = NULL, \
*NonhydroBEFluxM2d = NULL, *NonhydroPCETempFacialIntegral = NULL, *InvSquaHeight = NULL, \
*InvSHeight = NULL;//(Here, InvSquaHeight is used in mxAssembleGlobalStiffMatrix, InvSHeight is used in mxAssembleNonhydroRHS)

void SWENonhydro3dMemoryAllocation(int Np3d, int K3d, int IENfp, int IENe, int Nface3d,\
	int BENfp, int BENe, int BotENfp, int BotENe, int BotBENfp, int BotBENe, int SurfBENfp,\
	int SurfBENe, int Np2d, int K2d, int IENfp2d, int IENe2d, int Nface2d, int BENfp2d,\
	int BENe2d ){
	/*Water depth $H$ is not included*/
	/*In the main function*/
	NonhydroVariable = malloc(Np3d*K3d*3*sizeof(double));
	MemoryAllocationCheck(NonhydroVariable, Np3d*K3d * 3 * sizeof(double));
	PUPX = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(PUPX, Np3d*K3d*sizeof(double));
	PVPY = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(PVPY, Np3d*K3d*sizeof(double));
	PHPX = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(PHPX, Np3d*K3d*sizeof(double));
	PHPY = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(PHPY, Np3d*K3d*sizeof(double));
	PSPX = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(PSPX, Np3d*K3d*sizeof(double));
	PSPY = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(PSPY, Np3d*K3d*sizeof(double));
	TempPSPX = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(TempPSPX, Np3d*K3d*sizeof(double));
	TempPSPY = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(TempPSPY, Np3d*K3d*sizeof(double));

    /*For depth averaged velocity part*/
	NonhydroHU2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(NonhydroHU2d, Np2d*K2d*sizeof(double));
	NonhydroHV2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(NonhydroHV2d, Np2d*K2d*sizeof(double));
	TempNonhydroHU2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(TempNonhydroHU2d, Np2d*K2d*sizeof(double));
	TempNonhydroHV2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(TempNonhydroHV2d, Np2d*K2d*sizeof(double));
	TempNonhydrofield3d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(TempNonhydrofield3d, Np3d*K3d*sizeof(double));
	Nonhydrofmod = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(Nonhydrofmod, Np3d*K3d*sizeof(double));

	/*GetSecondOrderPartialDerivativeInHorizontalDirection*/
	NonhydroIEfm = malloc(IENe*IENfp * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEfm, IENe*IENfp * 3 * sizeof(double));
    NonhydroIEfp = malloc(IENe*IENfp * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEfp, IENe*IENfp * 3 * sizeof(double));
	NonhydroIEFluxM = malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxM, IENe*IENfp * 4 * sizeof(double));
	NonhydroIEFluxP = malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxP, IENe*IENfp * 4 * sizeof(double));
	NonhydroIEFluxS = malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxS, IENe*IENfp * 4 * sizeof(double));
	NonhydroERHS = malloc(4 * Np3d*K3d*Nface3d*sizeof(double));
	MemoryAllocationCheck(NonhydroERHS, 4 * Np3d*K3d*Nface3d*sizeof(double));

	NonhydroBEfm = malloc(BENe*BENfp * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEfm, BENe*BENfp * 3 * sizeof(double));
	NonhydroBEfp = malloc(BENe*BENfp * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEfp, BENe*BENfp * 3 * sizeof(double));
	NonhydroBEFluxM = malloc(BENe*BENfp * 4 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEFluxM, BENe*BENfp * 4 * sizeof(double));
	NonhydroBEFluxS = malloc(BENe*BENfp * 4 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEFluxS, BENe*BENfp * 4 * sizeof(double));
	NonhydroTempFacialIntegral = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(NonhydroTempFacialIntegral, Np3d*K3d*sizeof(double));
	NonhydroTempVolumeIntegral = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(NonhydroTempVolumeIntegral, Np3d*K3d*sizeof(double));

	/*GetFirstOrderPartialDerivativeInHorizontalDirection*/
	NonhydrozM = malloc(BENfp*BENe * sizeof(double));
	MemoryAllocationCheck(NonhydrozM, BENfp*BENe * sizeof(double));
	NonhydrozP = malloc(BENfp*BENe * sizeof(double));
	MemoryAllocationCheck(NonhydrozP, BENfp*BENe * sizeof(double));
	CoePS = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(CoePS, Np3d*K3d*sizeof(double));

	/*GetFirstOrderPartialDerivativeInVerticalDirection*/
	NonhydroBotEfm = malloc(BotENfp*BotENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBotEfm, BotENfp*BotENe * 3 * sizeof(double));
	NonhydroBotEfp = malloc(BotENfp*BotENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBotEfp, BotENfp*BotENe * 3 * sizeof(double));
	NonhydroBotEFluxM = malloc(BotENfp*BotENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBotEFluxM, BotENfp*BotENe * 3 * sizeof(double));
	NonhydroBotEFluxP = malloc(BotENfp*BotENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBotEFluxP, BotENfp*BotENe * 3 * sizeof(double));
	NonhydroBotEFluxS = malloc(BotENfp*BotENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBotEFluxS, BotENfp*BotENe * 3 * sizeof(double));
	NonhydroBotBEfm = malloc(BotBENfp*BotBENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBotBEfm, BotBENfp*BotBENe * 3 * sizeof(double));
	NonhydroBotBEFluxM = malloc(BotBENfp*BotBENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBotBEFluxM, BotBENfp*BotBENe * 3 * sizeof(double));
	NonhydroBotBEFluxS = malloc(BotBENfp*BotBENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBotBEFluxS, BotBENfp*BotBENe * 3 * sizeof(double));
	NonhydroSurfBEfm = malloc(SurfBENfp*SurfBENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroSurfBEfm, SurfBENfp*SurfBENe * 3 * sizeof(double));
	NonhydroSurfBEFluxM = malloc(SurfBENfp*SurfBENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroSurfBEFluxM, SurfBENfp*SurfBENe * 3 * sizeof(double));
	NonhydroSurfBEFluxS = malloc(SurfBENfp*SurfBENe * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroSurfBEFluxS, SurfBENfp*SurfBENe * 3 * sizeof(double));

	
	/*For surface vertical velocity and bottom vertical velocity part*/
	NonhydroRHS2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(NonhydroRHS2d, Np2d*K2d*sizeof(double));
	Weta = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(Weta, Np2d*K2d*sizeof(double));
	Wbot = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(Wbot, Np2d*K2d*sizeof(double));
	ueta = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ueta, Np2d*K2d*sizeof(double));
	ubot = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ubot, Np2d*K2d*sizeof(double));
	veta = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(veta, Np2d*K2d*sizeof(double));
	vbot = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(vbot, Np2d*K2d*sizeof(double));
	etax = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(etax, Np2d*K2d*sizeof(double));
	etay = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(etay, Np2d*K2d*sizeof(double));
	TempZx2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(TempZx2d, Np2d*K2d*sizeof(double));
	TempZy2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(TempZy2d, Np2d*K2d*sizeof(double));
	
	NonhydroIEfm2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEfm2d, IENfp2d*IENe2d * 3 * sizeof(double));
	NonhydroIEfp2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEfp2d, IENfp2d*IENe2d * 3 * sizeof(double));
	NonhydroIEFluxM2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxM2d, IENfp2d*IENe2d*sizeof(double));
	NonhydroIEFluxP2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxP2d, IENfp2d*IENe2d*sizeof(double));
	NonhydroIEFluxS2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxS2d, IENfp2d*IENe2d*sizeof(double));
	NonhydroERHS2d = malloc(Np2d*K2d*Nface2d*sizeof(double));
	MemoryAllocationCheck(NonhydroERHS2d, Np2d*K2d*Nface2d*sizeof(double));
	NonhydroPCEVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(NonhydroPCEVolumeIntegralX, Np2d*K2d*sizeof(double));
	NonhydroPCETempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(NonhydroPCETempVolumeIntegralX, Np2d*K2d*sizeof(double));
	NonhydroPCEVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(NonhydroPCEVolumeIntegralY, Np2d*K2d*sizeof(double));
	NonhydroPCETempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(NonhydroPCETempVolumeIntegralY, Np2d*K2d*sizeof(double));
	NonhydroBEfm2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEfm2d, BENe2d * BENfp2d * 3 * sizeof(double));
	NonhydroBEzM2d = malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(NonhydroBEzM2d, BENe2d * BENfp2d*sizeof(double));
	NonhydroBEfp2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEfp2d, BENe2d * BENfp2d * 3 * sizeof(double));
	NonhydroBEzP2d = malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(NonhydroBEzP2d, BENe2d * BENfp2d*sizeof(double));
	NonhydroBEFluxS2d = malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(NonhydroBEFluxS2d, BENe2d*BENfp2d*sizeof(double));
	NonhydroBEFluxM2d = malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(NonhydroBEFluxM2d, BENe2d*BENfp2d*sizeof(double));
	NonhydroPCETempFacialIntegral = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(NonhydroPCETempFacialIntegral, Np2d*K2d*sizeof(double));

	InvSquaHeight = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(InvSquaHeight, Np3d*K3d*sizeof(double));

	InvSHeight = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(InvSHeight, Np3d*K3d*sizeof(double));

	SWENonhydro3dInitialized = "True";
}

void SWENonhydro3dMemoryDeAllocation(){
	free(NonhydroVariable), NonhydroVariable = NULL;
	free(NonhydroERHS), NonhydroERHS = NULL;
	free(NonhydroHU2d), NonhydroHU2d = NULL;
	free(NonhydroHV2d), NonhydroHV2d = NULL; 
	free(TempNonhydroHU2d), TempNonhydroHU2d = NULL;
	free(TempNonhydroHV2d), TempNonhydroHV2d = NULL;
	free(TempNonhydrofield3d), TempNonhydrofield3d = NULL;
	free(Nonhydrofmod), Nonhydrofmod = NULL; 
	free(NonhydroVariable), NonhydroVariable = NULL;
	free(PUPX), PUPX = NULL;
	free(PVPY), PVPY = NULL;
	free(PHPX), PHPX = NULL;
	free(PHPY), PHPY = NULL;
	free(PSPX), PSPX = NULL; 
	free(PSPY), PSPY = NULL;
	free(TempPSPX), TempPSPX = NULL;
	free(TempPSPY), TempPSPY = NULL;
	free(NonhydroIEfm), NonhydroIEfm = NULL;
	free(NonhydroIEfp), NonhydroIEfp = NULL; 
	free(NonhydroIEFluxM), NonhydroIEFluxM = NULL;
	free(NonhydroIEFluxP), NonhydroIEFluxP = NULL;
	free(NonhydroIEFluxS), NonhydroIEFluxS = NULL;
	free(NonhydroERHS), NonhydroERHS = NULL; 
	free(NonhydroBEfm), NonhydroBEfm = NULL;
	free(NonhydroBEfp), NonhydroBEfp = NULL;
	free(NonhydroBEFluxM), NonhydroBEFluxM = NULL;
	free(NonhydroBEFluxS), NonhydroBEFluxS = NULL; 
	free(NonhydroTempFacialIntegral), NonhydroTempFacialIntegral = NULL;
	free(NonhydroTempVolumeIntegral), NonhydroTempVolumeIntegral = NULL;
	free(NonhydrozM), NonhydrozM = NULL;
	free(NonhydrozP), NonhydrozP = NULL; 
	free(TempPSPX), TempPSPX = NULL;
	free(TempPSPY), TempPSPY = NULL;
	free(CoePS), CoePS = NULL;
	free(NonhydroBotEfm), NonhydroBotEfm = NULL;
	free(NonhydroBotEfp), NonhydroBotEfp = NULL; 
	free(NonhydroBotEFluxM), NonhydroBotEFluxM = NULL;
	free(NonhydroBotEFluxP), NonhydroBotEFluxP = NULL;
	free(NonhydroBotEFluxS), NonhydroBotEFluxS = NULL; 
	free(NonhydroBotBEFluxM), NonhydroBotBEFluxM = NULL;
	free(NonhydroBotBEFluxS), NonhydroBotBEFluxS = NULL;
	free(NonhydroSurfBEfm), NonhydroSurfBEfm = NULL; 
	free(NonhydroSurfBEFluxM), NonhydroSurfBEFluxM = NULL;
	free(NonhydroSurfBEFluxS), NonhydroSurfBEFluxS = NULL;
	free(NonhydroRHS2d), NonhydroRHS2d = NULL; 
	free(Weta), Weta = NULL;
	free(Wbot), Wbot = NULL;
	free(ueta), ueta = NULL;
	free(ubot), ubot = NULL;
	free(veta), veta = NULL;
	free(vbot), vbot = NULL; 
	free(etax), etax = NULL;
	free(etay), etay = NULL;
	free(TempZx2d), TempZx2d = NULL;
	free(TempZy2d), TempZy2d = NULL;
	free(NonhydroIEfm2d), NonhydroIEfm2d = NULL; 
	free(NonhydroIEfp2d), NonhydroIEfp2d = NULL;
	free(NonhydroIEFluxM2d), NonhydroIEFluxM2d = NULL;
	free(NonhydroIEFluxP2d), NonhydroIEFluxP2d = NULL; 
	free(NonhydroIEFluxS2d), NonhydroIEFluxS2d = NULL;
	free(NonhydroERHS2d), NonhydroERHS2d = NULL;
	free(NonhydroPCEVolumeIntegralX), NonhydroPCEVolumeIntegralX = NULL; 
	free(NonhydroPCETempVolumeIntegralX), NonhydroPCETempVolumeIntegralX = NULL;
	free(NonhydroPCEVolumeIntegralY), NonhydroPCEVolumeIntegralY = NULL; 
	free(NonhydroPCETempVolumeIntegralY), NonhydroPCETempVolumeIntegralY = NULL;
	free(NonhydroBEfm2d), NonhydroBEfm2d = NULL;
	free(NonhydroBEzM2d), NonhydroBEzM2d = NULL; 
	free(NonhydroBEfp2d), NonhydroBEfp2d = NULL;
	free(NonhydroBEzP2d), NonhydroBEzP2d = NULL;
	free(NonhydroBEFluxS2d), NonhydroBEFluxS2d = NULL; 
	free(NonhydroBEFluxM2d), NonhydroBEFluxM2d = NULL;
	free(NonhydroPCETempFacialIntegral), NonhydroPCETempFacialIntegral = NULL;

	SWENonhydro3dInitialized = "False";
}

/*This is for vertical diffusion part*/
double *Tau = NULL;
char *VertDiffInitialized = "False";

void VertDiffMemoryAllocation(const int Np2d, int K2d, const int Nz){
	Tau = malloc(sizeof(double)*(Np2d*K2d*(Nz+1)));
    MemoryAllocationCheck(Tau, sizeof(double)*(Np2d*K2d*(Nz+1)));
	memset(Tau, 0, Np2d*K2d*(Nz + 1)*sizeof(double));
	VertDiffInitialized = "True";
}

void VertDiffMemoryDeAllocation(){
	free(Tau);
	Tau = NULL;
	VertDiffInitialized = "False";
}

/*This is for GOTM part*/
double *tkeGOTM = NULL, *epsGOTM = NULL, *LGOTM = NULL, *nuhGOTM = NULL,\
*numGOTM = NULL, *layerHeight = NULL, *huCentralDate = NULL, *hvCentralDate = NULL,\
*huVerticalLine = NULL, *hvVerticalLine = NULL, *shearFrequencyDate = NULL, *buoyanceFrequencyDate = NULL,\
*BottomFrictionLength = NULL, *BottomFrictionVelocity = NULL, *SurfaceFrictionLength = NULL,\
*SurfaceFrictionVelocity = NULL, *eddyViscosityDate = NULL;

char *GOTMInitialized = "False";

void GotmSolverMemoryAllocation(int Num2d, int Interface, int Np2d, int K3d){
	tkeGOTM = malloc(sizeof(double)*(Num2d*Interface)); 
    MemoryAllocationCheck(tkeGOTM, sizeof(double)*(Num2d*Interface));
	epsGOTM = malloc(sizeof(double)*(Num2d*Interface)); 
    MemoryAllocationCheck(epsGOTM, sizeof(double)*(Num2d*Interface));
	LGOTM = malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(LGOTM, sizeof(double)*(Num2d*Interface));
	nuhGOTM = malloc(sizeof(double)*(Num2d*Interface)); 
    MemoryAllocationCheck(nuhGOTM, sizeof(double)*(Num2d*Interface));
	numGOTM = malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(numGOTM, sizeof(double)*(Num2d*Interface));
	layerHeight = malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(layerHeight, sizeof(double)*(Num2d*Interface));
	huCentralDate = malloc(sizeof(double)*(Np2d*K3d));
    MemoryAllocationCheck(huCentralDate, sizeof(double)*(Np2d*K3d));
	hvCentralDate = malloc(sizeof(double)*(Np2d*K3d));
    MemoryAllocationCheck(hvCentralDate, sizeof(double)*(Np2d*K3d));
	huVerticalLine = malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(huVerticalLine, sizeof(double)*(Num2d*Interface));
	hvVerticalLine = malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(hvVerticalLine, sizeof(double)*(Num2d*Interface));
	shearFrequencyDate = malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(shearFrequencyDate, sizeof(double)*(Num2d*Interface));
	buoyanceFrequencyDate = malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(buoyanceFrequencyDate, sizeof(double)*(Num2d*Interface));
	BottomFrictionLength = malloc(sizeof(double)*Num2d);
    MemoryAllocationCheck(BottomFrictionLength, sizeof(double)*Num2d);
	BottomFrictionVelocity = malloc(sizeof(double)*Num2d);
    MemoryAllocationCheck(BottomFrictionVelocity, sizeof(double)*Num2d);
	SurfaceFrictionLength = malloc(sizeof(double)*Num2d);
    MemoryAllocationCheck(SurfaceFrictionLength, sizeof(double)*Num2d);
	SurfaceFrictionVelocity = malloc(sizeof(double)*Num2d);
    MemoryAllocationCheck(SurfaceFrictionVelocity, sizeof(double)*Num2d);
	eddyViscosityDate = malloc(sizeof(double)*(Num2d * Interface));
    MemoryAllocationCheck(eddyViscosityDate, sizeof(double)*(Num2d * Interface));
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

char *VertVelocityInitialized = "False";

void VertVelocitySolverMemoryAllocation(int Np2d, int K2d, int IENfp2d, int IENe2d, int Nface, int BENe2d, int BENfp2d, int Np3d,\
	int K3d, int IENfp3d, int IENe3d, int BENe3d, int BENfp3d){
	VSrhs2d = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(VSrhs2d,Np2d*K2d*sizeof(double));
	VSIEfm2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
    MemoryAllocationCheck(VSIEfm2d,IENfp2d*IENe2d * 3 * sizeof(double));
	VSIEfp2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
    MemoryAllocationCheck(VSIEfp2d,IENfp2d*IENe2d * 3 * sizeof(double));
	VSIEFluxM2d = malloc(IENfp2d*IENe2d*sizeof(double));
    MemoryAllocationCheck(VSIEFluxM2d,IENfp2d*IENe2d*sizeof(double));
	VSIEFluxP2d = malloc(IENfp2d*IENe2d*sizeof(double));
    MemoryAllocationCheck(VSIEFluxP2d, IENfp2d*IENe2d*sizeof(double));
	VSIEFluxS2d = malloc(IENfp2d*IENe2d*sizeof(double));
    MemoryAllocationCheck(VSIEFluxS2d, IENfp2d*IENe2d*sizeof(double));
	VSERHS2d = malloc(Np2d*K2d*Nface*sizeof(double));
    MemoryAllocationCheck(VSERHS2d, Np2d*K2d*Nface*sizeof(double));
	VSVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(VSVolumeIntegralX,Np2d*K2d*sizeof(double));
	VSTempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(VSTempVolumeIntegralX, Np2d*K2d*sizeof(double));
	VSVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(VSVolumeIntegralY,Np2d*K2d*sizeof(double));
	VSTempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(VSTempVolumeIntegralY,Np2d*K2d*sizeof(double));
	VSBEfm2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
    MemoryAllocationCheck(VSBEfm2d,BENe2d * BENfp2d * 3 * sizeof(double));
	VSBEzM2d = malloc(BENe2d * BENfp2d*sizeof(double));
    MemoryAllocationCheck(VSBEzM2d,BENe2d * BENfp2d*sizeof(double));
	VSBEfp2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
    MemoryAllocationCheck(VSBEfp2d,BENe2d * BENfp2d * 3 * sizeof(double));
	VSBEzP2d = malloc(BENe2d * BENfp2d*sizeof(double));
    MemoryAllocationCheck(VSBEzP2d,BENe2d * BENfp2d*sizeof(double));
	VSBEFluxS2d = malloc(BENe2d*BENfp2d*sizeof(double));
    MemoryAllocationCheck(VSBEFluxS2d,BENe2d*BENfp2d*sizeof(double));
	VSBEFluxM2d = malloc(BENe2d*BENfp2d*sizeof(double));
    MemoryAllocationCheck(VSBEFluxM2d,BENe2d*BENfp2d*sizeof(double));
	VSTempFacialIntegral = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(VSTempFacialIntegral,Np2d*K2d*sizeof(double));
	VSfield2d = malloc(Np3d*K3d*sizeof(double));
    MemoryAllocationCheck(VSfield2d,Np3d*K3d*sizeof(double));
	VSrhs3d = malloc(Np3d*K3d*sizeof(double));
    MemoryAllocationCheck(VSrhs3d,Np3d*K3d*sizeof(double));
	VSIEfm3d = malloc(IENfp3d*IENe3d * 3 * sizeof(double));
    MemoryAllocationCheck(VSIEfm3d,IENfp3d*IENe3d * 3 * sizeof(double));
	VSIEfp3d = malloc(IENfp3d*IENe3d * 3 * sizeof(double));
    MemoryAllocationCheck(VSIEfp3d,IENfp3d*IENe3d * 3 * sizeof(double));
	VSIEFluxM3d = malloc(IENfp3d*IENe3d*sizeof(double));
    MemoryAllocationCheck(VSIEFluxM3d,IENfp3d*IENe3d*sizeof(double));
	VSIEFluxP3d = malloc(IENfp3d*IENe3d*sizeof(double));
    MemoryAllocationCheck(VSIEFluxP3d,IENfp3d*IENe3d*sizeof(double));
	VSIEFluxS3d = malloc(IENfp3d*IENe3d*sizeof(double));
    MemoryAllocationCheck(VSIEFluxS3d,IENfp3d*IENe3d*sizeof(double));
	VSERHS3d = malloc(Np3d*K3d*Nface*sizeof(double));
    MemoryAllocationCheck(VSERHS3d,Np3d*K3d*Nface*sizeof(double));
	VSVolumeIntegralX3d = malloc(Np3d*K3d*sizeof(double));
    MemoryAllocationCheck(VSVolumeIntegralX3d,Np3d*K3d*sizeof(double));
	VSTempVolumeIntegralX3d = malloc(Np3d*K3d*sizeof(double));
    MemoryAllocationCheck(VSTempVolumeIntegralX3d,Np3d*K3d*sizeof(double));
	VSVolumeIntegralY3d = malloc(Np3d*K3d*sizeof(double));
    MemoryAllocationCheck(VSVolumeIntegralY3d,Np3d*K3d*sizeof(double));
	VSTempVolumeIntegralY3d = malloc(Np3d*K3d*sizeof(double));
    MemoryAllocationCheck(VSTempVolumeIntegralY3d,Np3d*K3d*sizeof(double));
	VSBEfm3d = malloc(BENe3d * BENfp3d * 3 * sizeof(double));
    MemoryAllocationCheck(VSBEfm3d,BENe3d * BENfp3d * 3 * sizeof(double));
	VSBEzM3d = malloc(BENe3d * BENfp3d * sizeof(double));
    MemoryAllocationCheck(VSBEzM3d,BENe3d * BENfp3d * sizeof(double));
	VSBEfp3d = malloc(BENe3d * BENfp3d * 3 * sizeof(double));
    MemoryAllocationCheck(VSBEfp3d,BENe3d * BENfp3d * 3 * sizeof(double));
	VSBEzP3d = malloc(BENe3d * BENfp3d * sizeof(double));
    MemoryAllocationCheck(VSBEzP3d,BENe3d * BENfp3d * sizeof(double));
	VSBEFluxS3d = malloc(BENe3d*BENfp3d*sizeof(double));
    MemoryAllocationCheck(VSBEFluxS3d,BENe3d*BENfp3d*sizeof(double));
	VSBEFluxM3d = malloc(BENe3d*BENfp3d*sizeof(double));
    MemoryAllocationCheck(VSBEFluxM3d,BENe3d*BENfp3d*sizeof(double));
	VSTempVerticalVelocity = malloc(Np3d*K3d*sizeof(double));
    MemoryAllocationCheck(VSTempVerticalVelocity,Np3d*K3d*sizeof(double));
	VSBotVertVelocity = malloc(Np3d*K2d*sizeof(double));
    MemoryAllocationCheck(VSBotVertVelocity,Np3d*K2d*sizeof(double));
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

char *PCEInitialized = "False";

void PCEMemoryAllocation(int IENfp2d, int IENe2d, int Np2d, int K2d, int Nface, int BENe2d, int BENfp2d){
	IEfm2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
    MemoryAllocationCheck(IEfm2d,IENfp2d*IENe2d * 3 * sizeof(double));
	IEfp2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
    MemoryAllocationCheck(IEfp2d,IENfp2d*IENe2d * 3 * sizeof(double));
	IEFluxM2d = malloc(IENfp2d*IENe2d*sizeof(double));
    MemoryAllocationCheck(IEFluxM2d,IENfp2d*IENe2d*sizeof(double));
	IEFluxP2d = malloc(IENfp2d*IENe2d*sizeof(double));
    MemoryAllocationCheck(IEFluxP2d,IENfp2d*IENe2d*sizeof(double));
	IEFluxS2d = malloc(IENfp2d*IENe2d*sizeof(double));
    MemoryAllocationCheck(IEFluxS2d,IENfp2d*IENe2d*sizeof(double));
	ERHS2d = malloc(Np2d*K2d*Nface*sizeof(double));
    MemoryAllocationCheck(ERHS2d,Np2d*K2d*Nface*sizeof(double));
	PCEVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(PCEVolumeIntegralX,Np2d*K2d*sizeof(double));
	PCETempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(PCETempVolumeIntegralX,Np2d*K2d*sizeof(double));
	PCEVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(PCEVolumeIntegralY,Np2d*K2d*sizeof(double));
	PCETempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(PCETempVolumeIntegralY,Np2d*K2d*sizeof(double));
	BEfm2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
    MemoryAllocationCheck(BEfm2d,BENe2d * BENfp2d * 3 * sizeof(double));
	BEzM2d = malloc(BENe2d * BENfp2d*sizeof(double));
    MemoryAllocationCheck(BEzM2d,BENe2d * BENfp2d*sizeof(double));
	BEfp2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
    MemoryAllocationCheck(BEfp2d,BENe2d * BENfp2d * 3 * sizeof(double));
	BEzP2d = malloc(BENe2d * BENfp2d*sizeof(double));
    MemoryAllocationCheck(BEzP2d,BENe2d * BENfp2d*sizeof(double));
	BEFluxS2d = malloc(BENe2d*BENfp2d*sizeof(double));
    MemoryAllocationCheck(BEFluxS2d,BENe2d*BENfp2d*sizeof(double));
	BEFluxM2d = malloc(BENe2d*BENfp2d*sizeof(double));
    MemoryAllocationCheck(BEFluxM2d,BENe2d*BENfp2d*sizeof(double));
	PCETempFacialIntegral = malloc(Np2d*K2d*sizeof(double));
    MemoryAllocationCheck(PCETempFacialIntegral,Np2d*K2d*sizeof(double));
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

char *AdvInitialized = "False";

void AdvMemoryAllocation(int Np, int K, int Nvar, int IENfp, int IENe, int Nface, int BENfp, int BENe, int BotENfp,\
                        int BotENe, int BotBENfp, int BotBENe, int SurfBENfp, int SurfBENe){
	TempFacialIntegral = malloc(Np*K*Nvar*sizeof(double));
    MemoryAllocationCheck(TempFacialIntegral,Np*K*Nvar*sizeof(double));
	IEfm = malloc(IENfp*IENe*(Nvar + 1)*sizeof(double));
    MemoryAllocationCheck(IEfm,IENfp*IENe*(Nvar + 1)*sizeof(double));
	IEfp = malloc(IENfp*IENe*(Nvar + 1)*sizeof(double));
    MemoryAllocationCheck(IEfp,IENfp*IENe*(Nvar + 1)*sizeof(double));
	IEFluxM = malloc(IENfp*IENe*Nvar*sizeof(double));
    MemoryAllocationCheck(IEFluxM,IENfp*IENe*Nvar*sizeof(double));
	IEFluxP = malloc(IENfp*IENe*Nvar*sizeof(double));
    MemoryAllocationCheck(IEFluxP,IENfp*IENe*Nvar*sizeof(double));
	IEFluxS = malloc(IENfp*IENe*Nvar*sizeof(double));
    MemoryAllocationCheck(IEFluxS,IENfp*IENe*Nvar*sizeof(double));
	ERHS = malloc(Np*K*Nvar*Nface*sizeof(double));
    MemoryAllocationCheck(ERHS,Np*K*Nvar*Nface*sizeof(double));
	BEfm = malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
    MemoryAllocationCheck(BEfm,BENfp*BENe*(Nvar + 1)*sizeof(double));
	BEfp = malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
    MemoryAllocationCheck(BEfp,BENfp*BENe*(Nvar + 1)*sizeof(double));
	AdvzM = malloc(BENfp*BENe*sizeof(double));
    MemoryAllocationCheck(AdvzM,BENfp*BENe*sizeof(double));
	AdvzP = malloc(BENfp*BENe*sizeof(double));
    MemoryAllocationCheck(AdvzP,BENfp*BENe*sizeof(double));
	BEFluxM = malloc(BENfp*BENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BEFluxM,BENfp*BENe*Nvar*sizeof(double));
	BEFluxS = malloc(BENfp*BENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BEFluxS,BENfp*BENe*Nvar*sizeof(double));
	BotEfm = malloc(BotENfp*BotENe*(Nvar + 2)*sizeof(double));
    MemoryAllocationCheck(BotEfm,BotENfp*BotENe*(Nvar + 2)*sizeof(double));
	BotEfp = malloc(BotENfp*BotENe*(Nvar + 2)*sizeof(double));
    MemoryAllocationCheck(BotEfp,BotENfp*BotENe*(Nvar + 2)*sizeof(double));
	BotEFluxM = malloc(BotENfp*BotENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotEFluxM,BotENfp*BotENe*Nvar*sizeof(double));
	BotEFluxP = malloc(BotENfp*BotENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotEFluxP,BotENfp*BotENe*Nvar*sizeof(double));
	BotEFluxS = malloc(BotENfp*BotENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotEFluxS,BotENfp*BotENe*Nvar*sizeof(double));
	BotBEfm = malloc(BotBENfp*BotBENe*(Nvar + 2)*sizeof(double));
    MemoryAllocationCheck(BotBEfm,BotBENfp*BotBENe*(Nvar + 2)*sizeof(double));
	BotBEFluxM = malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotBEFluxM,BotBENfp*BotBENe*Nvar*sizeof(double));
	BotBEFluxS = malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotBEFluxS,BotBENfp*BotBENe*Nvar*sizeof(double));
	SurfBEfm = malloc(SurfBENfp*SurfBENe*(Nvar + 2)*sizeof(double));
    MemoryAllocationCheck(SurfBEfm,SurfBENfp*SurfBENe*(Nvar + 2)*sizeof(double));
	SurfBEFluxM = malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
    MemoryAllocationCheck(SurfBEFluxM,SurfBENfp*SurfBENe*Nvar*sizeof(double));
	SurfBEFluxS = malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
    MemoryAllocationCheck(SurfBEFluxS,SurfBENfp*SurfBENe*Nvar*sizeof(double));
	/*Allocate memory for E, G and H, and calculate these volume flux term*/
	E = malloc(Np*K*Nvar*sizeof(double));
    MemoryAllocationCheck(E,Np*K*Nvar*sizeof(double));
	G = malloc(Np*K*Nvar*sizeof(double));
    MemoryAllocationCheck(G,Np*K*Nvar*sizeof(double));
	H = malloc(Np*K*Nvar*sizeof(double));
    MemoryAllocationCheck(H,Np*K*Nvar*sizeof(double));
	TempVolumeIntegral = malloc(Np*K*sizeof(double));
    MemoryAllocationCheck(TempVolumeIntegral,Np*K*sizeof(double));
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
    MemoryAllocationCheck(HorDiffzM,BENfp*BENe*sizeof(double));
	HorDiffzP = malloc(BENfp*BENe*sizeof(double));
    MemoryAllocationCheck(HorDiffzP,BENfp*BENe*sizeof(double));
	HorDiffnv = malloc(Np*K*sizeof(double));
    MemoryAllocationCheck(HorDiffnv,Np*K*sizeof(double));
	HorDiffTempBEfp = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
    MemoryAllocationCheck(HorDiffTempBEfp,BENe*BENfp*(Nfield + 1)*sizeof(double));
	HorDiffTempBEfm = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
    MemoryAllocationCheck(HorDiffTempBEfm,BENe*BENfp*(Nfield + 1)*sizeof(double));
	/*Allocate memory for the original variable over boundary edge*/
	HorDiffBEfp = malloc(BENfp*BENe*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffBEfp,BENfp*BENe*Nfield*sizeof(double));
	/*Allocate memory for the local face value at boundary edge*/
	HorDiffBEfm = malloc(BENe*BENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffBEfm,BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the original variable $u,v$ and $\theta$*/
	HorDiffvariable = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffvariable,Np*K*Nfield*sizeof(double));
	/*Allocate memory for auxiallary variable $q_x=\frac{\partial u(v,\theta)}{\partial x}$*/
	HorDiffAVx = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVx,Np*K*Nfield*sizeof(double));
	/*Allocate memory for auxiallary variable $q_y=\frac{\partial u(v,\theta)}{\partial y}$*/
	HorDiffAVy = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVy,Np*K*Nfield*sizeof(double));
	/*Allocate memory for part of the volumn integral of auxiallary variable $q_{x1}=\bold{r_x}\cdot (D_r*u(v,\theta))$,
	and for part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffVx = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffVx,Np*K*Nfield*sizeof(double));
	/*Allocate memory for the rest part of the volumn integral of auxiallary variable $q_{x2}=\bold{r_x}\cdot (D_r*u(v,\theta))$,
	and for the rest part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffTempVx = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffTempVx,Np*K*Nfield*sizeof(double));
	/*Allocate memory for part of the volumn integral of auxiallary variable $q_{y1}=\bold{r_y}\cdot (D_r*u(v,\theta))$,
	and for part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffVy = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffVy,Np*K*Nfield*sizeof(double));
	/*Allocate memory for the rest part of the volumn integral of auxiallary variable $q_{y2}=\bold{s_y}\cdot (D_s*u(v,\theta))$,
	and for the rest part of the volumn integral part of the second order operator $\frac{\partial q}{\partial y}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$*/
	HorDiffTempVy = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffTempVy,Np*K*Nfield*sizeof(double));
	/*Allocate memory for the local face value at inner edge*/
	HorDiffIEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEfm,IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at inner edge*/
	HorDiffAVIEfm = malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVIEfm,IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent face value at inner edge*/
	HorDiffIEfp = malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEfp,IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent face value of the auxialary variable at inner edge*/
	HorDiffAVIEfp = malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVIEfp,IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local flux term at inner edge*/
	HorDiffIEFluxM = malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEFluxM,IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the adjacent flux term at inner edge*/
	HorDiffIEFluxP = malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEFluxP,IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at inner edge*/
	HorDiffIEFluxS = malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEFluxS,IENe*IENfp*Nfield*sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at boundary edge*/
	HorDiffAVBEfm = malloc(BENe*BENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVBEfm,BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the local flux term at boundary edge*/
	HorDiffBEFluxM = malloc(BENe*BENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffBEFluxM,BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for the numerical flux term at boundary edge*/
	HorDiffBEFluxS = malloc(BENe*BENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffBEFluxS,BENe*BENfp*Nfield*sizeof(double));
	/*Allocate memory for interior edge contribution to right hand side in x direction*/
	HorDiffERHSX = malloc(Np*K*Nfield*Nface*sizeof(double));
    MemoryAllocationCheck(HorDiffERHSX,Np*K*Nfield*Nface*sizeof(double));
	/*Allocate memory for interior edge contribution to right hand side in y direction*/
	HorDiffERHSY = malloc(Np*K*Nfield*Nface*sizeof(double));
    MemoryAllocationCheck(HorDiffERHSY,Np*K*Nfield*Nface*sizeof(double));

	/*Allocate memory for local primitive diffusion term in both x and y direction, $\nu\nabla u(v,\theta)$*/
	/*This part is used because we need to extract the trace of the local derivative operator to compute the numerical flux*/
	HorDiffLocalPrimitiveDiffTermX = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLocalPrimitiveDiffTermX,Np*K*Nfield*sizeof(double));
	HorDiffLocalPrimitiveDiffTermY = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLocalPrimitiveDiffTermY,Np*K*Nfield*sizeof(double));
	HorDiffLPDTIEfm = malloc(IENfp*IENe*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLPDTIEfm,IENfp*IENe*Nfield*sizeof(double));
	HorDiffLPDTIEfp = malloc(IENfp*IENe*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLPDTIEfp,IENfp*IENe*Nfield*sizeof(double));
	HorDiffLPDTBEfm = malloc(BENfp*BENe*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLPDTBEfm,BENfp*BENe*Nfield*sizeof(double));
	HorDiffTempFacialIntegralX = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffTempFacialIntegralX,Np*K*Nfield*sizeof(double));
	HorDiffTempFacialIntegralY = malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffTempFacialIntegralY,Np*K*Nfield*sizeof(double));
	HorDiffInnerEdgeTau = malloc(IENe*IENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffInnerEdgeTau,IENe*IENfp*sizeof(double));
	HorDiffBoundaryEdgeTau = malloc(BENe*BENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffBoundaryEdgeTau,BENe*BENfp*sizeof(double));
	HorDiffIEnvfm = malloc(IENe*IENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffIEnvfm,IENe*IENfp*sizeof(double));
	HorDiffIEnvfp = malloc(IENe*IENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffIEnvfp,IENe*IENfp*sizeof(double));
	HorDiffBEnvfm = malloc(BENe*BENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffBEnvfm,BENe*BENfp*sizeof(double));
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