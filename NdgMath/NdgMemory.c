# include "NdgMemory.h"

void MemoryAllocationCheck(double *dest, int size){
    while(dest == NULL){
        dest = malloc(size); 
    }
}

/*The following part is for calculation of density, and is called from mxCalculateDensityField.c*/
double *BaroclinicT = NULL, *BaroclinicS = NULL, *BaroclinicDTS = NULL;

char *BaroDensityInitialized = "False";

void BaroDensityMemoryAllocation( int Np, int K) {
	BaroclinicT = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicT, Np*K * sizeof(double));
	BaroclinicS = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicS, Np*K * sizeof(double));
	BaroclinicDTS = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicDTS, Np*K * sizeof(double));
	BaroDensityInitialized = "True";
}

void BaroDensityMemoryDeAllocation(int Np, int K) {
	free(BaroclinicT), BaroclinicT = NULL;
	free(BaroclinicS), BaroclinicS = NULL;
	free(BaroclinicDTS), BaroclinicDTS = NULL;
	BaroDensityInitialized = "False";
}

/*The following part is for calculation of baroclinic term, and is called from mxCalculateBaroclinicTerm.c*/

double *BaroclinicPDPX = NULL, *BaroclinicPDPY = NULL, *BaroclinicPRHOPX = NULL, *BaroclinicPRHOPY = NULL,\
*BaroclinicPRHOPS = NULL, *BaroclinicInXPartOne = NULL, *BaroclinicInXPartTwo = NULL, *BaroclinicInYPartOne = NULL,\
*BaroclinicInYPartTwo = NULL, *BaroclinicInXTempRHS = NULL, *BaroclinicInYTempRHS = NULL, *Baroclinicfmod = NULL,\
*BaroclinicBotEfm = NULL, *BaroclinicBotEfp = NULL, *BaroclinicBotEFluxM = NULL, *BaroclinicBotEFluxP = NULL, \
*BaroclinicBotEFluxS = NULL, *BaroclinicIEfm = NULL, *BaroclinicIEfp = NULL, *BaroclinicIEfluxM = NULL, \
*BaroclinicIEfluxP = NULL, *BaroclinicIEfluxS = NULL, *BaroclinicERHS = NULL, *BaroclinicTempFacialIntegral = NULL,\
*BaroclinicTempVolumeIntegral = NULL, *BaroclinicBEfm = NULL, *BaroclinicBEfp = NULL, *BaroclinicBEfluxM = NULL, \
*BaroclinicBEfluxS = NULL;

char *BaroclinicPartInitialized = "False";

void BaroclinicPartMemoryAllocation(int Np, int K, int K2d, int BotENe, int BotENfp, int IENe, \
	int IENfp, int Nface,int BENe, int BENfp) {
	BaroclinicPDPX = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPDPX, Np*K * sizeof(double));
	BaroclinicPDPY = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPDPY, Np*K * sizeof(double));
	BaroclinicPRHOPX = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPRHOPX, Np*K * sizeof(double));
	BaroclinicPRHOPY = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPRHOPY, Np*K * sizeof(double));
	BaroclinicPRHOPS = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPRHOPS, Np*K * sizeof(double));
	BaroclinicInXPartOne = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInXPartOne, Np*K * sizeof(double));
	BaroclinicInXPartTwo = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInXPartTwo, Np*K * sizeof(double));
	BaroclinicInYPartOne = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInYPartOne, Np*K * sizeof(double));
	BaroclinicInYPartTwo = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInYPartTwo, Np*K * sizeof(double));
	BaroclinicInXTempRHS = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInXTempRHS, Np*K * sizeof(double));
	BaroclinicInYTempRHS = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInYTempRHS, Np*K * sizeof(double));
	Baroclinicfmod = malloc(Np*K2d * sizeof(double));
	MemoryAllocationCheck(Baroclinicfmod, Np*K2d * sizeof(double));
	BaroclinicBotEfm = malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEfm, Np*K * sizeof(double));
	BaroclinicBotEfp = malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEfp, BotENe*BotENfp * sizeof(double));
	BaroclinicBotEFluxM = malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEFluxM, BotENe*BotENfp * sizeof(double));
	BaroclinicBotEFluxP = malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEFluxP, BotENe*BotENfp * sizeof(double));
	BaroclinicBotEFluxS = malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEFluxS, BotENe*BotENfp * sizeof(double));
	BaroclinicTempFacialIntegral = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicTempFacialIntegral, Np*K * sizeof(double));
	BaroclinicIEfm = malloc(IENe*IENfp * 2 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfm, IENe*IENfp * 2 * sizeof(double));
	BaroclinicIEfp = malloc(IENe*IENfp * 2 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfp, IENe*IENfp * 2 * sizeof(double));
	BaroclinicIEfluxM = malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfluxM, IENe*IENfp * 4 * sizeof(double));
	BaroclinicIEfluxP = malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfluxP, IENe*IENfp * 4 * sizeof(double));
	BaroclinicIEfluxS = malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfluxS, IENe*IENfp * 4 * sizeof(double));
	BaroclinicERHS = malloc(4 * Np*K*(Nface - 2) * sizeof(double));
	MemoryAllocationCheck(BaroclinicERHS, 4 * Np*K*(Nface - 2) * sizeof(double));
	BaroclinicTempVolumeIntegral = malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicTempVolumeIntegral, Np*K * sizeof(double));
	BaroclinicBEfm = malloc(BENe*BENfp * 2 * sizeof(double));
	MemoryAllocationCheck(BaroclinicBEfm, BENe*BENfp * 2 * sizeof(double));
	BaroclinicBEfp = malloc(BENe*BENfp * 3 * sizeof(double));
	MemoryAllocationCheck(BaroclinicBEfp, BENe*BENfp * 3 * sizeof(double));
	BaroclinicBEfluxM = malloc(BENe*BENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicBEfluxM, BENe*BENfp * 4 * sizeof(double));
	BaroclinicBEfluxS = malloc(BENe*BENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicBEfluxS, BENe*BENfp * 4 * sizeof(double));
	BaroclinicPartInitialized = "True";
}

void BaroclinicPartMemoryDeAllocation() {
	free(BaroclinicPDPX), BaroclinicPDPX = NULL;
	free(BaroclinicPDPY), BaroclinicPDPY = NULL;
	free(BaroclinicPRHOPX), BaroclinicPRHOPX = NULL;
	free(BaroclinicPRHOPY), BaroclinicPRHOPY = NULL;
	free(BaroclinicPRHOPS), BaroclinicPRHOPS = NULL;
	free(BaroclinicInXPartOne), BaroclinicInXPartOne = NULL;
	free(BaroclinicInXPartTwo), BaroclinicInXPartTwo = NULL;
	free(BaroclinicInYPartOne), BaroclinicInYPartOne = NULL;
	free(BaroclinicInYPartTwo), BaroclinicInYPartTwo = NULL;
	free(BaroclinicInXTempRHS), BaroclinicInXTempRHS = NULL;
	free(BaroclinicInYTempRHS), BaroclinicInYTempRHS = NULL;
	free(Baroclinicfmod), Baroclinicfmod = NULL;
	free(BaroclinicBotEfm), BaroclinicBotEfm = NULL;
	free(BaroclinicBotEfp), BaroclinicBotEfp = NULL;
	free(BaroclinicBotEFluxM), BaroclinicBotEFluxM = NULL;
	free(BaroclinicBotEFluxP), BaroclinicBotEFluxP = NULL;
	free(BaroclinicBotEFluxS), BaroclinicBotEFluxS = NULL;
	free(BaroclinicTempFacialIntegral), BaroclinicTempFacialIntegral = NULL;
	free(BaroclinicIEfm), BaroclinicIEfm = NULL;
	free(BaroclinicIEfp), BaroclinicIEfp = NULL;
	free(BaroclinicIEfluxM), BaroclinicIEfluxM = NULL;
	free(BaroclinicIEfluxP), BaroclinicIEfluxP = NULL;
	free(BaroclinicIEfluxS), BaroclinicIEfluxS = NULL;
	free(BaroclinicERHS), BaroclinicERHS = NULL;
	free(BaroclinicTempFacialIntegral), BaroclinicTempFacialIntegral = NULL;
	free(BaroclinicTempVolumeIntegral), BaroclinicTempVolumeIntegral = NULL;
	free(BaroclinicBEfm), BaroclinicBEfm = NULL;
	free(BaroclinicBEfp), BaroclinicBEfp = NULL;
	free(BaroclinicBEfluxM), BaroclinicBEfluxM = NULL;
	free(BaroclinicBEfluxS), BaroclinicBEfluxS = NULL;
	BaroclinicPartInitialized = "False";
}

double *TimeIntervalu = NULL, *TimeIntervalv = NULL, *TimeIntervalw = NULL, *TimeIntervalt = NULL;
char *TimeIntervalInitialized = "False";

void SWENH3dTimeIntervalMemoryAllocation(int Np, int K3d, int K2d){
	TimeIntervalu = malloc(Np*K3d*sizeof(double));
	MemoryAllocationCheck(TimeIntervalu, Np*K3d*sizeof(double));
	TimeIntervalv = malloc(Np*K3d*sizeof(double));
	MemoryAllocationCheck(TimeIntervalv, Np*K3d*sizeof(double));
	TimeIntervalw = malloc(Np*K3d*sizeof(double));
	MemoryAllocationCheck(TimeIntervalw, Np*K3d*sizeof(double));
	TimeIntervalt = malloc(K2d*sizeof(double));
	MemoryAllocationCheck(TimeIntervalw, K2d*sizeof(double));

	TimeIntervalInitialized = "True";
}

void SWENH3dTimeIntervalMemoryDeAllocation(){
	free(TimeIntervalu); TimeIntervalu = NULL;
	free(TimeIntervalv); TimeIntervalv = NULL;
	free(TimeIntervalw); TimeIntervalw = NULL;
	free(TimeIntervalt); TimeIntervalt = NULL;
	TimeIntervalInitialized = "False";
}

/*This is for three-dimensional non-hydrostatic model*/

char *SWENonhydro3dInitialized = "False";

double *NonhydroHU2d = NULL, *NonhydroHV2d = NULL, \
*TempNonhydroHU2d = NULL, *TempNonhydroHV2d = NULL, *TempNonhydrofield3d = NULL, *Nonhydrofmod = NULL, \
*NonhydroVariable = NULL, *PSPX = NULL, \
*PSPY = NULL, *TempPSPX = NULL,*TempPSPY = NULL, *NonhydroIEfm = NULL, *NonhydroIEfp = NULL, \
*NonhydroIEFluxM = NULL, *NonhydroIEFluxP = NULL, *NonhydroIEFluxS = NULL, *NonhydroERHS = NULL, \
*NonhydroBEfm = NULL, *NonhydroBEfp = NULL, *NonhydroBEFluxM = NULL, *NonhydroBEFluxS = NULL, \
*NonhydroTempFacialIntegral = NULL, *NonhydroTempVolumeIntegral = NULL,*NonhydrozM = NULL, *NonhydrozP = NULL, \
*CoePS = NULL,*NonhydroBotEfm = NULL,*NonhydroBotEfp = NULL, *NonhydroBotEFluxM = NULL,*NonhydroBotEFluxP = NULL,*NonhydroBotEFluxS = NULL, \
*NonhydroBotBEfm = NULL,*NonhydroBotBEFluxM = NULL,*NonhydroBotBEFluxS = NULL,*NonhydroSurfBEfm = NULL, \
*NonhydroSurfBEFluxM = NULL,*NonhydroSurfBEFluxS = NULL,*NonhydroRHS2d = NULL, \
*Weta = NULL,*Wbot = NULL,*ueta = NULL,*ubot = NULL,*veta = NULL,*vbot = NULL, \
*etax = NULL,*etay = NULL, *TempZx2d = NULL,*TempZy2d = NULL,*NonhydroIEfm2d = NULL, \
*NonhydroIEfp2d = NULL,*NonhydroIEFluxM2d = NULL,*NonhydroIEFluxP2d = NULL, \
*NonhydroIEFluxS2d = NULL,*NonhydroERHS2d = NULL,*NonhydroPCEVolumeIntegralX = NULL, \
*NonhydroPCETempVolumeIntegralX = NULL,*NonhydroPCEVolumeIntegralY = NULL, \
*NonhydroPCETempVolumeIntegralY = NULL, *NonhydroBEfm2d = NULL, *NonhydroBEzM2d = NULL, \
*NonhydroBEfp2d = NULL,*NonhydroBEzP2d = NULL,*NonhydroBEFluxS2d = NULL, \
*NonhydroBEFluxM2d = NULL,*NonhydroPCETempFacialIntegral = NULL, *NonhydroIEfmod = NULL, \
*NonhydroBEfmod = NULL;
//double *PHPX = NULL, *PHPY = NULL;

/*The following space is allocated in file mxCalculatePartialDerivative.c and file mxCalculatePartialDerivativeUpdated.c.
* We note that, ONLY ONE function can be called any time, Or the program would crash
*/
void SWENonhydro3dMemoryAllocation(int Np3d, int K3d, int IENfp, int IENe, int Nface3d,\
	int BENfp, int BENe, int BotENfp, int BotENe, int BotBENfp, int BotBENe, int SurfBENfp,\
	int SurfBENe, int Np2d, int K2d, int IENfp2d, int IENe2d, int Nface2d, int BENfp2d,\
	int BENe2d ){
	/*Water depth $H$ is not included*/
	/*In the main function*/
	NonhydroVariable = malloc(Np3d*K3d*3*sizeof(double));
	MemoryAllocationCheck(NonhydroVariable, Np3d*K3d * 3 * sizeof(double));
	/*
	PHPX = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(PHPX, Np3d*K3d*sizeof(double));
	PHPY = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(PHPY, Np3d*K3d*sizeof(double));
	*/
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
	NonhydroIEFluxM = malloc(IENe*IENfp * 6 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxM, IENe*IENfp * 6 * sizeof(double));
	NonhydroIEFluxP = malloc(IENe*IENfp * 6 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxP, IENe*IENfp * 6 * sizeof(double));
	NonhydroIEFluxS = malloc(IENe*IENfp * 6 * sizeof(double));
	MemoryAllocationCheck(NonhydroIEFluxS, IENe*IENfp * 6 * sizeof(double));
	/*Space allocated for the right hand correspoinding to PUPX PUPY PVPX PVPY PHPX PHPY,
	This space is bigger than that for PUPS PVPS and PWPS
	*/ 
	NonhydroERHS = malloc(6 * Np3d*K3d*(Nface3d-2)*sizeof(double));
	MemoryAllocationCheck(NonhydroERHS, 6 * Np3d*K3d*(Nface3d - 2) *sizeof(double));

	NonhydroBEfm = malloc(BENe*BENfp * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEfm, BENe*BENfp * 3 * sizeof(double));
	NonhydroBEfp = malloc(BENe*BENfp * 3 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEfp, BENe*BENfp * 3 * sizeof(double));
	NonhydroBEFluxM = malloc(BENe*BENfp * 6 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEFluxM, BENe*BENfp * 6 * sizeof(double));
	NonhydroBEFluxS = malloc(BENe*BENfp * 6 * sizeof(double));
	MemoryAllocationCheck(NonhydroBEFluxS, BENe*BENfp * 6 * sizeof(double));
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

	/*The following memory is allocated in the updated version of function */
	NonhydroIEfmod = malloc(IENe2d*IENfp*sizeof(double));
	MemoryAllocationCheck(NonhydroIEfmod, IENe2d*IENfp*sizeof(double));
	NonhydroBEfmod = malloc(BENe2d*BENfp*sizeof(double));
	MemoryAllocationCheck(NonhydroBEfmod, BENe2d*BENfp*sizeof(double));

	SWENonhydro3dInitialized = "True";
}

void SWENonhydro3dMemoryDeAllocation(){
	free(NonhydroVariable), NonhydroVariable = NULL;
	/*
	free(PHPX), PHPX = NULL;
	free(PHPY), PHPY = NULL;
	*/
	free(PSPX), PSPX = NULL;
	free(PSPY), PSPY = NULL;
	free(TempPSPX), TempPSPX = NULL;
	free(TempPSPY), TempPSPY = NULL;
	free(NonhydroHU2d), NonhydroHU2d = NULL;
	free(NonhydroHV2d), NonhydroHV2d = NULL;
	free(TempNonhydroHU2d), TempNonhydroHU2d = NULL;
	free(TempNonhydroHV2d), TempNonhydroHV2d = NULL;
	free(TempNonhydrofield3d), TempNonhydrofield3d = NULL;
	free(Nonhydrofmod), Nonhydrofmod = NULL;
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
	free(CoePS), CoePS = NULL;
	free(NonhydroBotEfm), NonhydroBotEfm = NULL;
	free(NonhydroBotEfp), NonhydroBotEfp = NULL;
	free(NonhydroBotEFluxM), NonhydroBotEFluxM = NULL;
	free(NonhydroBotEFluxP), NonhydroBotEFluxP = NULL;
	free(NonhydroBotEFluxS), NonhydroBotEFluxS = NULL;
	free(NonhydroBotBEfm); NonhydroBotBEfm = NULL;
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
	free(NonhydroIEfmod), NonhydroIEfmod = NULL;
	free(NonhydroBEfmod), NonhydroBEfmod = NULL;
	SWENonhydro3dInitialized = "False";
}

/*The following global space is used in file mxCalculateBottomVerticalVelocity.c*/
double *NonhydroUbot = NULL, *NonhydroVbot = NULL, *NonhydroHbot = NULL, \
*NonhydroZxbot = NULL, *NonhydroZybot = NULL;

char *NHVertVelocityInitialized = "False";
void SWENonhydroVertVelocityMemoryAllocation( int Nfp, int Ne){
	NonhydroUbot = malloc(Nfp*Ne*sizeof(double));
	MemoryAllocationCheck(NonhydroUbot, Nfp*Ne*sizeof(double));
	NonhydroVbot = malloc(Nfp*Ne*sizeof(double));
	MemoryAllocationCheck(NonhydroVbot, Nfp*Ne*sizeof(double));
	NonhydroHbot = malloc(Nfp*Ne*sizeof(double));
	MemoryAllocationCheck(NonhydroHbot, Nfp*Ne*sizeof(double));
	NonhydroZxbot = malloc(Nfp*Ne*sizeof(double));
	MemoryAllocationCheck(NonhydroZxbot, Nfp*Ne*sizeof(double));
	NonhydroZybot = malloc(Nfp*Ne*sizeof(double));
	MemoryAllocationCheck(NonhydroZybot, Nfp*Ne*sizeof(double));
	NHVertVelocityInitialized = "True";
}

void SWENonhydroVertVelocityMemoryDeAllocation(){
	free(NonhydroUbot), NonhydroUbot = NULL;
	free(NonhydroVbot), NonhydroVbot = NULL;
	free(NonhydroHbot), NonhydroHbot = NULL;
	free(NonhydroZxbot), NonhydroZxbot = NULL;
	free(NonhydroZybot), NonhydroZybot = NULL;
	NHVertVelocityInitialized = "False";
}


/*The following global space is used in file mxAssembleGlobalStiffMatrixNew.c*/

double *K33 = NULL, *InvSquaHeight = NULL, *InnerEdgeTau = NULL, *BottomEdgeTau = NULL, *SurfaceEdgeTau = NULL;
char *GlobalStiffMatrixInitialized = "False";

void GlobalStiffMatrixMemoryAllocation(int Np3d, int K3d, int IENe, int BotENe, int SurfBENe){
	K33 = malloc(sizeof(double)*Np3d*K3d);
	MemoryAllocationCheck(K33, sizeof(double)*Np3d*K3d);
	InvSquaHeight = malloc(sizeof(double)*Np3d*K3d);
	MemoryAllocationCheck(InvSquaHeight, sizeof(double)*Np3d*K3d);
	InnerEdgeTau = malloc(IENe*sizeof(double));
	MemoryAllocationCheck(InnerEdgeTau, sizeof(double)*IENe);
	BottomEdgeTau = malloc(BotENe*sizeof(double));
	MemoryAllocationCheck(BottomEdgeTau, BotENe*sizeof(double));
	SurfaceEdgeTau = malloc(SurfBENe*sizeof(double));
	MemoryAllocationCheck(SurfaceEdgeTau, SurfBENe*sizeof(double));
	GlobalStiffMatrixInitialized = "True";
}

void GlobalStiffMatrixMemoryDeAllocation(){
	free(K33); K33 = NULL;
	free(InvSquaHeight); InvSquaHeight = NULL;
	free(InnerEdgeTau); InnerEdgeTau = NULL;
	free(BottomEdgeTau); BottomEdgeTau = NULL;
	free(SurfaceEdgeTau); SurfaceEdgeTau = NULL;
	GlobalStiffMatrixInitialized = "False";
}

/*The following global space is used in file mxImposeBoundaryCondition.c*/

double *ImposeBCsInvSquaHeight = NULL, *ImposeBCsK33 = NULL, *BETau = NULL, *ImposeBCsNewmannData = NULL, *ImposeBCsWx = NULL, *ImposeBCsWy = NULL, \
*ImposeBCsWxRHS2d = NULL, *ImposeBCsWyRHS2d = NULL, *ImposeBCsWIEFluxMx2d = NULL, *ImposeBCsWIEFluxMy2d = NULL, *ImposeBCsWIEFluxPx2d = NULL, *ImposeBCsWIEFluxPy2d = NULL, \
*ImposeBCsWIEFluxSx2d = NULL, *ImposeBCsWIEFluxSy2d = NULL, *ImposeBCsVolumeIntegralX = NULL, *ImposeBCsTempVolumeIntegralX = NULL, *ImposeBCsVolumeIntegralY = NULL, \
*ImposeBCsTempVolumeIntegralY = NULL, *ImposeBCsIEfm = NULL, *ImposeBCsIEfp = NULL, *ImposeBCsERHSx = NULL, *ImposeBCsERHSy = NULL, *ImposeBCsTempFacialIntegral = NULL, \
*ImposeBCsBotBEU = NULL, *ImposeBCsBotBEV = NULL, *ImposeBCsBotBEH = NULL, *ImposeBCsBotBECoe = NULL, *ImposeBCsBotBEPWPS = NULL, \
*ImposeBCsWs = NULL, *ImposeBCsPUPX = NULL, *ImposeBCsPUPY = NULL, *ImposeBCsPUPS = NULL, *ImposeBCsPVPX = NULL, *ImposeBCsPVPY = NULL, \
*ImposeBCsPVPS = NULL, *ImposeBCsTempNewmannData = NULL, *ImposeBCsPWPH = NULL, *ImposeBCsPHPX = NULL, *ImposeBCsPHPY = NULL, \
*ImposeBCsPzPx = NULL, *ImposeBCsPzPy = NULL, *ImposeBCsTempDataInX = NULL, *ImposeBCsTempDataInY = NULL, *ImposeBCsTempDataInZ = NULL, \
*ImposeBCsInverseH = NULL, *ImposeBCsPNPS = NULL, *ImposeBCsPNPX = NULL, *ImposeBCsPNPY = NULL, *ImposeBCsGPEtaPX = NULL, \
*ImposeBCsDPNPSX = NULL, *ImposeBCsGPEtaPY = NULL, *ImposeBCsDPNPSY = NULL, *ImposeBCsInvHeight = NULL;
char *ImposeBoundaryInitialized = "False";

void SWENH3dImposeBoundaryMemoryAllocation(int Np, int K, int BENe, int Np2d, int K2d, int Nface2d, int IENfp2d, int IENe2d, int BotBENe, int BotBENfp){
	ImposeBCsInvSquaHeight = malloc(sizeof(double)*Np*K);
	MemoryAllocationCheck(ImposeBCsInvSquaHeight, sizeof(double)*Np*K);
	ImposeBCsK33 = malloc(sizeof(double)*Np*K);
	MemoryAllocationCheck(ImposeBCsK33, sizeof(double)*Np*K);
	BETau = malloc(BENe*sizeof(double));
	MemoryAllocationCheck(BETau, sizeof(double)*BENe);
	ImposeBCsNewmannData = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsNewmannData, Np2d*K2d*sizeof(double));
	ImposeBCsWx = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWx, Np2d*K2d*sizeof(double));
	ImposeBCsWy = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWy, Np2d*K2d*sizeof(double));
	ImposeBCsWxRHS2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWxRHS2d, Np2d*K2d*sizeof(double));
	ImposeBCsWyRHS2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWyRHS2d, Np2d*K2d*sizeof(double));
	ImposeBCsWIEFluxMx2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWIEFluxMx2d, IENfp2d*IENe2d*sizeof(double));
	ImposeBCsWIEFluxMy2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWIEFluxMy2d, IENfp2d*IENe2d*sizeof(double));
	ImposeBCsWIEFluxPx2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWIEFluxPx2d, IENfp2d*IENe2d*sizeof(double));
	ImposeBCsWIEFluxPy2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWIEFluxPy2d, IENfp2d*IENe2d*sizeof(double));
	ImposeBCsWIEFluxSx2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWIEFluxSx2d, IENfp2d*IENe2d*sizeof(double));
	ImposeBCsWIEFluxSy2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWIEFluxSy2d, IENfp2d*IENe2d*sizeof(double));
	ImposeBCsVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsVolumeIntegralX, Np2d*K2d*sizeof(double));
	ImposeBCsTempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsTempVolumeIntegralX, Np2d*K2d*sizeof(double));
	ImposeBCsVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsVolumeIntegralY, Np2d*K2d*sizeof(double));
	ImposeBCsTempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsTempVolumeIntegralY, Np2d*K2d*sizeof(double));
	ImposeBCsIEfm = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsIEfm, IENfp2d*IENe2d*sizeof(double));
	ImposeBCsIEfp = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsIEfp, IENfp2d*IENe2d*sizeof(double));
	ImposeBCsERHSx = malloc(Np2d*K2d*Nface2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsERHSx, Np2d*K2d*Nface2d*sizeof(double));
	ImposeBCsERHSy = malloc(Np2d*K2d*Nface2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsERHSy, Np2d*K2d*Nface2d*sizeof(double));
	ImposeBCsTempFacialIntegral = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(ImposeBCsTempFacialIntegral, Np2d*K2d*sizeof(double));
	ImposeBCsBotBEU = malloc(BotBENe*BotBENfp*sizeof(double));
	MemoryAllocationCheck(ImposeBCsBotBEU, BotBENe*BotBENfp*sizeof(double));
	ImposeBCsBotBEV = malloc(BotBENe*BotBENfp*sizeof(double));
	MemoryAllocationCheck(ImposeBCsBotBEV, BotBENe*BotBENfp*sizeof(double));
	ImposeBCsBotBEH = malloc(BotBENe*BotBENfp*sizeof(double));
	MemoryAllocationCheck(ImposeBCsBotBEH, BotBENe*BotBENfp*sizeof(double));
	ImposeBCsBotBECoe = malloc(BotBENe*BotBENfp*sizeof(double));
	MemoryAllocationCheck(ImposeBCsBotBECoe, BotBENe*BotBENfp*sizeof(double));
	ImposeBCsBotBEPWPS = malloc(BotBENe*BotBENfp*sizeof(double));
	MemoryAllocationCheck(ImposeBCsBotBEPWPS, BotBENe*BotBENfp*sizeof(double));
	ImposeBCsWs = malloc(BotBENe*BotBENfp*sizeof(double));
	MemoryAllocationCheck(ImposeBCsWs, BotBENe*BotBENfp*sizeof(double));
	ImposeBCsPUPX = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPUPX, BotBENe*BotBENfp * sizeof(double));
	ImposeBCsPUPY = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPUPY, BotBENe*BotBENfp * sizeof(double));
	ImposeBCsPUPS = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPUPS, BotBENe*BotBENfp * sizeof(double));
	ImposeBCsPVPX = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPVPX, BotBENe*BotBENfp * sizeof(double));
	ImposeBCsPVPY = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPVPY, BotBENe*BotBENfp * sizeof(double));
	ImposeBCsPVPS = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPVPS, BotBENe*BotBENfp * sizeof(double));
	ImposeBCsTempNewmannData = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsTempNewmannData, BotBENe*BotBENfp * sizeof(double));
	/*$\frac{\partial w}{\partial H}$*/
	ImposeBCsPWPH = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPWPH, BotBENe*BotBENfp * sizeof(double));
	/*$\frac{\partial H}{\partial x}$*/
	ImposeBCsPHPX = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPHPX, BotBENe*BotBENfp * sizeof(double));
	/*$\frac{\partial H}{\partial y}$*/
	ImposeBCsPHPY = malloc(BotBENe*BotBENfp * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPHPY, BotBENe*BotBENfp * sizeof(double));
	ImposeBCsPzPx = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPzPx, Np2d*K2d * sizeof(double));
	ImposeBCsPzPy = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPzPy, Np2d*K2d * sizeof(double));
	ImposeBCsTempDataInX = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsTempDataInX, Np2d*K2d * sizeof(double));
	ImposeBCsTempDataInY = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsTempDataInY, Np2d*K2d * sizeof(double));
	ImposeBCsTempDataInZ = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsTempDataInZ, Np2d*K2d * sizeof(double));
	ImposeBCsInverseH = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsInverseH, Np2d*K2d * sizeof(double));
	ImposeBCsPNPS = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPNPS, Np2d*K2d * sizeof(double));
	ImposeBCsPNPX = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPNPX, Np2d*K2d * sizeof(double));
	ImposeBCsPNPY = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsPNPY, Np2d*K2d * sizeof(double));
	ImposeBCsGPEtaPX = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsGPEtaPX, Np2d*K2d * sizeof(double));
	ImposeBCsDPNPSX = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsDPNPSX, Np2d*K2d * sizeof(double));
	ImposeBCsGPEtaPY = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsGPEtaPY, Np2d*K2d * sizeof(double));
	ImposeBCsDPNPSY = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(ImposeBCsDPNPSY, Np2d*K2d * sizeof(double));
	ImposeBCsInvHeight = malloc(Np*	K * sizeof(double));
	MemoryAllocationCheck(ImposeBCsInvHeight, Np*K * sizeof(double));
	ImposeBoundaryInitialized = "True";
}

void SWENH3dImposeBoundaryMemoryDeAllocation(){
	free(ImposeBCsInvSquaHeight); ImposeBCsInvSquaHeight = NULL;
	free(ImposeBCsK33); ImposeBCsK33 = NULL;
	free(BETau); BETau = NULL;
	free(ImposeBCsNewmannData), ImposeBCsNewmannData = NULL;
	free(ImposeBCsWx), ImposeBCsWx = NULL;
	free(ImposeBCsWy), ImposeBCsWy = NULL;
	free(ImposeBCsWxRHS2d), ImposeBCsWxRHS2d = NULL;
	free(ImposeBCsWyRHS2d), ImposeBCsWyRHS2d = NULL;
	free(ImposeBCsWIEFluxMx2d), ImposeBCsWIEFluxMx2d = NULL;
	free(ImposeBCsWIEFluxMy2d), ImposeBCsWIEFluxMy2d = NULL;
	free(ImposeBCsWIEFluxPx2d), ImposeBCsWIEFluxPx2d = NULL;
	free(ImposeBCsWIEFluxPy2d), ImposeBCsWIEFluxPy2d = NULL;
	free(ImposeBCsWIEFluxSx2d), ImposeBCsWIEFluxSx2d = NULL;
	free(ImposeBCsWIEFluxSy2d), ImposeBCsWIEFluxSy2d = NULL;
	free(ImposeBCsVolumeIntegralX), ImposeBCsVolumeIntegralX = NULL;
	free(ImposeBCsTempVolumeIntegralX), ImposeBCsTempVolumeIntegralX = NULL;
	free(ImposeBCsVolumeIntegralY), ImposeBCsVolumeIntegralY = NULL;
	free(ImposeBCsTempVolumeIntegralY), ImposeBCsTempVolumeIntegralY = NULL;
	free(ImposeBCsIEfm), ImposeBCsIEfm = NULL;
	free(ImposeBCsIEfp), ImposeBCsIEfp = NULL;
	free(ImposeBCsERHSx), ImposeBCsERHSx = NULL;
	free(ImposeBCsERHSy), ImposeBCsERHSy = NULL;
	free(ImposeBCsTempFacialIntegral), ImposeBCsTempFacialIntegral = NULL;
	free(ImposeBCsBotBEU), ImposeBCsBotBEU = NULL;
	free(ImposeBCsBotBEV), ImposeBCsBotBEV = NULL;
	free(ImposeBCsBotBEH), ImposeBCsBotBEH = NULL;
	free(ImposeBCsBotBECoe), ImposeBCsBotBECoe = NULL;
	free(ImposeBCsBotBEPWPS), ImposeBCsBotBEPWPS = NULL;
	free(ImposeBCsWs), ImposeBCsWs = NULL;
	free(ImposeBCsPUPX), ImposeBCsPUPX = NULL;
	free(ImposeBCsPUPY), ImposeBCsPUPY = NULL;
	free(ImposeBCsPUPS), ImposeBCsPUPS = NULL;
	free(ImposeBCsPVPX), ImposeBCsPVPX = NULL;
	free(ImposeBCsPVPY), ImposeBCsPVPY = NULL;
	free(ImposeBCsPVPS), ImposeBCsPVPS = NULL;
	free(ImposeBCsTempNewmannData), ImposeBCsTempNewmannData = NULL;
	free(ImposeBCsPWPH), ImposeBCsPWPH = NULL;
	free(ImposeBCsPHPX), ImposeBCsPHPX = NULL;
	free(ImposeBCsPHPY), ImposeBCsPHPY = NULL;
	free(ImposeBCsPzPx), ImposeBCsPzPx = NULL;
	free(ImposeBCsPzPy), ImposeBCsPzPy = NULL;
	free(ImposeBCsTempDataInX), ImposeBCsTempDataInX = NULL;
	free(ImposeBCsTempDataInY), ImposeBCsTempDataInY = NULL;
	free(ImposeBCsTempDataInZ), ImposeBCsTempDataInZ = NULL;
	free(ImposeBCsInverseH), ImposeBCsInverseH = NULL;
	free(ImposeBCsPNPS), ImposeBCsPNPS = NULL;
	free(ImposeBCsPNPX), ImposeBCsPNPX = NULL;
	free(ImposeBCsPNPY), ImposeBCsPNPY = NULL;
	free(ImposeBCsGPEtaPX), ImposeBCsGPEtaPX = NULL;
	free(ImposeBCsDPNPSX), ImposeBCsDPNPSX = NULL;
	free(ImposeBCsGPEtaPY), ImposeBCsGPEtaPY = NULL;
	free(ImposeBCsDPNPSY), ImposeBCsDPNPSY = NULL;
	free(ImposeBCsInvHeight), ImposeBCsInvHeight = NULL;
	ImposeBoundaryInitialized = "False";
}


/*This is for vertical diffusion part*/
double *Tau = NULL, *u2d = NULL, *v2d = NULL;
char *VertDiffInitialized = "False";

void VertDiffMemoryAllocation(const int Np2d, int K2d, const int Nz){
	Tau = malloc(sizeof(double)*(Np2d*K2d*(Nz+1)));
    MemoryAllocationCheck(Tau, sizeof(double)*(Np2d*K2d*(Nz+1)));
	memset(Tau, 0, Np2d*K2d*(Nz + 1)*sizeof(double));
	u2d = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(u2d, sizeof(double)*(Np2d*K2d));
	v2d = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(v2d, sizeof(double)*(Np2d*K2d));
	VertDiffInitialized = "True";
}

void VertDiffMemoryDeAllocation(){
	free(Tau); Tau = NULL;
	free(u2d); u2d = NULL;
	free(v2d); v2d = NULL;
	VertDiffInitialized = "False";
}

/*This is for sparse vertical diffusion part*/
double *ImTau = NULL, *Imu2d = NULL, *Imv2d = NULL, *GlobalSystemRHS = NULL;
char *ImVertDiffInitialized = "False";

void ImVertDiffMemoryAllocation(const int Np2d, int K2d, const int Nz, int Np3d, int Nvar) {
	ImTau = malloc(sizeof(double)*(Np2d*K2d*(Nz + 1)));
	MemoryAllocationCheck(ImTau, sizeof(double)*(Np2d*K2d*(Nz + 1)));
	memset(ImTau, 0, Np2d*K2d*(Nz + 1) * sizeof(double));
	Imu2d = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(Imu2d, sizeof(double)*(Np2d*K2d));
	Imv2d = malloc(Np2d*K2d * sizeof(double));
	MemoryAllocationCheck(Imv2d, sizeof(double)*(Np2d*K2d));
	GlobalSystemRHS = malloc(Np3d*K2d * Nz * Nvar * sizeof(double));
	MemoryAllocationCheck(GlobalSystemRHS, Np3d*K2d * Nz * Nvar * sizeof(double));
	ImVertDiffInitialized = "True";
}

void ImVertDiffMemoryDeAllocation() {
	free(ImTau); ImTau = NULL;
	free(Imu2d); Imu2d = NULL;
	free(Imv2d); Imv2d = NULL;
	free(GlobalSystemRHS); GlobalSystemRHS = NULL;
	ImVertDiffInitialized = "False";
}

/*This is for GOTM part*/
double *tkeGOTM = NULL, *epsGOTM = NULL, *LGOTM = NULL, *nuhGOTM = NULL, \
*numGOTM = NULL, *layerHeight = NULL, *huCentralDate = NULL, *hvCentralDate = NULL, \
*huVerticalLine = NULL, *hvVerticalLine = NULL, *shearFrequencyDate = NULL, *buoyanceFrequencyDate = NULL, \
*BottomFrictionLength = NULL, *BottomFrictionVelocity = NULL, *SurfaceFrictionLength = NULL, \
*SurfaceFrictionVelocity = NULL, *eddyViscosityDate = NULL, *rhoCentralDate = NULL, *rhoVerticalLine = NULL, \
*hcenter = NULL, *huCentralDateNew = NULL, *hvCentralDateNew = NULL, *huVerticalLineNew = NULL, \
*hvVerticalLineNew = NULL, *hTCentralData = NULL, *hSCentralData = NULL, \
*hTVerticalLine = NULL, *hSVerticalLine = NULL;
//*eddyDiffusionDate = NULL, *eddyTKEDate = NULL, *eddyLengthDate = NULL, *eddyEPSDate = NULL, 

char *GOTMInitialized = "False";

void GotmSolverMemoryAllocation( int Interface, int Np2d, int K3d, int K2d ){
	tkeGOTM = malloc(sizeof(double)*(K2d*Interface)); 
    MemoryAllocationCheck(tkeGOTM, sizeof(double)*(K2d*Interface));
	epsGOTM = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(epsGOTM, sizeof(double)*(K2d*Interface));
	LGOTM = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(LGOTM, sizeof(double)*(K2d*Interface));
	nuhGOTM = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(nuhGOTM, sizeof(double)*(K2d*Interface));
	numGOTM = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(numGOTM, sizeof(double)*(K2d*Interface));
	layerHeight = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(layerHeight, sizeof(double)*(K2d*Interface));
	huCentralDate = malloc(sizeof(double)*(K3d));
    MemoryAllocationCheck(huCentralDate, sizeof(double)*(K3d));
	huCentralDateNew = malloc(sizeof(double)*(K3d));
	MemoryAllocationCheck(huCentralDateNew, sizeof(double)*(K3d));
	hvCentralDate = malloc(sizeof(double)*(K3d));
    MemoryAllocationCheck(hvCentralDate, sizeof(double)*(K3d));
	hvCentralDateNew = malloc(sizeof(double)*(K3d));
	MemoryAllocationCheck(hvCentralDateNew, sizeof(double)*(K3d));
	//The space occupied by the first data in huVerticalLine is not used
	huVerticalLine = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(huVerticalLine, sizeof(double)*(K2d*Interface));
	huVerticalLineNew = malloc(sizeof(double)*(K2d*Interface));
	MemoryAllocationCheck(huVerticalLineNew, sizeof(double)*(K2d*Interface));
	hvVerticalLine = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(hvVerticalLine, sizeof(double)*(K2d*Interface));
	hvVerticalLineNew = malloc(sizeof(double)*(K2d*Interface));
	MemoryAllocationCheck(hvVerticalLineNew, sizeof(double)*(K2d*Interface));
	shearFrequencyDate = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(shearFrequencyDate, sizeof(double)*(K2d*Interface));
	buoyanceFrequencyDate = malloc(sizeof(double)*(K2d*Interface));
    MemoryAllocationCheck(buoyanceFrequencyDate, sizeof(double)*(K2d*Interface));
	BottomFrictionLength = malloc(sizeof(double)*K2d);
    MemoryAllocationCheck(BottomFrictionLength, sizeof(double)*K2d);
	BottomFrictionVelocity = malloc(sizeof(double)*K2d);
    MemoryAllocationCheck(BottomFrictionVelocity, sizeof(double)*K2d);
	SurfaceFrictionLength = malloc(sizeof(double)*K2d);
    MemoryAllocationCheck(SurfaceFrictionLength, sizeof(double)*K2d);
	SurfaceFrictionVelocity = malloc(sizeof(double)*K2d);
    MemoryAllocationCheck(SurfaceFrictionVelocity, sizeof(double)*K2d);
	eddyViscosityDate = malloc(sizeof(double)*(K2d * Interface));
    MemoryAllocationCheck(eddyViscosityDate, sizeof(double)*(K2d * Interface));
	rhoCentralDate = malloc(sizeof(double)*(K3d));
	MemoryAllocationCheck(rhoCentralDate, sizeof(double)*(K3d));
	rhoVerticalLine = malloc(sizeof(double)*(K2d*Interface));
	MemoryAllocationCheck(rhoVerticalLine, sizeof(double)*(K2d*Interface));

	/*
	eddyDiffusionDate = malloc(sizeof(double)*(K2d * Interface));
	MemoryAllocationCheck(eddyDiffusionDate, sizeof(double)*(K2d * Interface));
	eddyTKEDate = malloc(sizeof(double)*(K2d * Interface));
	MemoryAllocationCheck(eddyTKEDate, sizeof(double)*(K2d * Interface));
	eddyLengthDate = malloc(sizeof(double)*(K2d * Interface));
	MemoryAllocationCheck(eddyLengthDate, sizeof(double)*(K2d * Interface));
	eddyEPSDate = malloc(sizeof(double)*(K2d * Interface));
	MemoryAllocationCheck(eddyEPSDate, sizeof(double)*(K2d * Interface));
	*/
	hcenter = malloc(K2d * sizeof(double));
	MemoryAllocationCheck(hcenter, sizeof(double)*K2d);
	hTCentralData = malloc(K3d * sizeof(double));
	hSCentralData = malloc(K3d * sizeof(double));
	hTVerticalLine = malloc(sizeof(double)*(K2d*Interface));
	hSVerticalLine = malloc(sizeof(double)*(K2d*Interface));
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
	free(huCentralDateNew); huCentralDateNew = NULL;
	free(hvCentralDate); hvCentralDate = NULL;
	free(hvCentralDateNew); hvCentralDateNew = NULL;
	free(huVerticalLine); huVerticalLine = NULL;
	free(huVerticalLineNew); huVerticalLineNew = NULL;
	free(hvVerticalLine); hvVerticalLine = NULL;
	free(hvVerticalLineNew); hvVerticalLineNew = NULL;
	free(shearFrequencyDate); shearFrequencyDate = NULL;
	free(buoyanceFrequencyDate); buoyanceFrequencyDate = NULL;
	free(BottomFrictionLength); BottomFrictionLength = NULL;
	free(BottomFrictionVelocity); BottomFrictionVelocity = NULL;
	free(SurfaceFrictionLength); SurfaceFrictionLength = NULL;
	free(SurfaceFrictionVelocity); SurfaceFrictionVelocity = NULL;
	free(eddyViscosityDate); eddyViscosityDate = NULL;
	free(rhoCentralDate); rhoCentralDate = NULL;
	free(rhoVerticalLine); rhoVerticalLine = NULL;

	/*
	free(eddyDiffusionDate); eddyDiffusionDate = NULL;
	free(eddyTKEDate); eddyTKEDate = NULL;
	free(eddyLengthDate); eddyLengthDate = NULL;
	free(eddyEPSDate); eddyEPSDate = NULL;
	*/
	free(hcenter); hcenter = NULL;
	free(hTCentralData); hTCentralData = NULL;
	free(hSCentralData); hSCentralData = NULL;
	free(hTVerticalLine); hTVerticalLine = NULL;
	free(hSVerticalLine); hSVerticalLine = NULL;
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
	VSTempFacialIntegral3d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(VSTempFacialIntegral3d, Np3d*K3d*sizeof(double));
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
	free(VSTempFacialIntegral3d), VSTempFacialIntegral3d = NULL;
	VertVelocityInitialized = "False";
}

/*This is for updated vertical velocity solver part*/
double *UpdatedVSrhs2d = NULL, *UpdatedVSIEfm2d = NULL, *UpdatedVSIEfp2d = NULL, *UpdatedVSIEFluxM2d = NULL, \
*UpdatedVSIEFluxP2d = NULL, *UpdatedVSIEFluxS2d = NULL, *UpdatedVSERHS2d = NULL, *UpdatedVSVolumeIntegralX = NULL, \
*UpdatedVSTempVolumeIntegralX = NULL, *UpdatedVSVolumeIntegralY = NULL, *UpdatedVSTempVolumeIntegralY = NULL, \
*UpdatedVSBEfm2d = NULL, *UpdatedVSBEzM2d = NULL, *UpdatedVSBEfp2d = NULL, *UpdatedVSBEzP2d = NULL, *UpdatedVSBEFluxS2d = NULL, \
*UpdatedVSBEFluxM2d = NULL, *UpdatedVSTempFacialIntegral = NULL, *UpdatedVSfield2d = NULL, *UpdatedVSrhs3d = NULL, \
*UpdatedVSIEfm3d = NULL, *UpdatedVSIEfp3d = NULL, *UpdatedVSIEFluxM3d = NULL, *UpdatedVSIEFluxP3d = NULL, *UpdatedVSIEFluxS3d = NULL, \
*UpdatedVSERHS3d = NULL, *UpdatedVSVolumeIntegralX3d = NULL, *UpdatedVSTempVolumeIntegralX3d = NULL, \
*UpdatedVSVolumeIntegralY3d = NULL, *UpdatedVSTempVolumeIntegralY3d = NULL, *UpdatedVSBEfm3d = NULL, \
*UpdatedVSBEzM3d = NULL, *UpdatedVSBEfp3d = NULL, *UpdatedVSBEzP3d = NULL, *UpdatedVSBEFluxS3d = NULL, *UpdatedVSBEFluxM3d = NULL, \
*UpdatedVSTempFacialIntegral3d = NULL, *UpdatedVSIEfmod = NULL, *UpdatedVSBEfmod = NULL, *Updatedfmod = NULL;

char *UpdatedVertVelocityInitialized = "False";

void UpdatedVertVelocitySolverMemoryAllocation(int Np2d, int K2d, int IENfp2d, int IENe2d, int Nface, int BENe2d, int BENfp2d, int Np3d, \
	int K3d, int IENfp3d, int IENe3d, int BENe3d, int BENfp3d){
	UpdatedVSrhs2d = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSrhs2d, Np2d*K2d*sizeof(double));
	UpdatedVSIEfm2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfm2d, IENfp2d*IENe2d * 3 * sizeof(double));
	UpdatedVSIEfp2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfp2d, IENfp2d*IENe2d * 3 * sizeof(double));
	UpdatedVSIEFluxM2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxM2d, IENfp2d*IENe2d*sizeof(double));
	UpdatedVSIEFluxP2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxP2d, IENfp2d*IENe2d*sizeof(double));
	UpdatedVSIEFluxS2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxS2d, IENfp2d*IENe2d*sizeof(double));
	UpdatedVSERHS2d = malloc(Np2d*K2d*Nface*sizeof(double));
	MemoryAllocationCheck(UpdatedVSERHS2d, Np2d*K2d*Nface*sizeof(double));
	UpdatedVSVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSVolumeIntegralX, Np2d*K2d*sizeof(double));
	UpdatedVSTempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempVolumeIntegralX, Np2d*K2d*sizeof(double));
	UpdatedVSVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSVolumeIntegralY, Np2d*K2d*sizeof(double));
	UpdatedVSTempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempVolumeIntegralY, Np2d*K2d*sizeof(double));
	UpdatedVSBEfm2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfm2d, BENe2d * BENfp2d * 3 * sizeof(double));
	UpdatedVSBEzM2d = malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEzM2d, BENe2d * BENfp2d*sizeof(double));
	UpdatedVSBEfp2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfp2d, BENe2d * BENfp2d * 3 * sizeof(double));
	UpdatedVSBEzP2d = malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEzP2d, BENe2d * BENfp2d*sizeof(double));
	UpdatedVSBEFluxS2d = malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxS2d, BENe2d*BENfp2d*sizeof(double));
	UpdatedVSBEFluxM2d = malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxM2d, BENe2d*BENfp2d*sizeof(double));
	UpdatedVSTempFacialIntegral = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempFacialIntegral, Np2d*K2d*sizeof(double));
	UpdatedVSfield2d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSfield2d, Np3d*K3d*sizeof(double));
	UpdatedVSrhs3d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSrhs3d, Np3d*K3d*sizeof(double));
	UpdatedVSIEfm3d = malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfm3d, IENfp3d*IENe3d * 3 * sizeof(double));
	UpdatedVSIEfp3d = malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfp3d, IENfp3d*IENe3d * 3 * sizeof(double));
	UpdatedVSIEFluxM3d = malloc(IENfp3d*IENe3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxM3d, IENfp3d*IENe3d*sizeof(double));
	UpdatedVSIEFluxP3d = malloc(IENfp3d*IENe3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxP3d, IENfp3d*IENe3d*sizeof(double));
	UpdatedVSIEFluxS3d = malloc(IENfp3d*IENe3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxS3d, IENfp3d*IENe3d*sizeof(double));
	UpdatedVSERHS3d = malloc(Np3d*K3d*Nface*sizeof(double));
	MemoryAllocationCheck(UpdatedVSERHS3d, Np3d*K3d*Nface*sizeof(double));
	UpdatedVSVolumeIntegralX3d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSVolumeIntegralX3d, Np3d*K3d*sizeof(double));
	UpdatedVSTempVolumeIntegralX3d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempVolumeIntegralX3d, Np3d*K3d*sizeof(double));
	UpdatedVSVolumeIntegralY3d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSVolumeIntegralY3d, Np3d*K3d*sizeof(double));
	UpdatedVSTempVolumeIntegralY3d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempVolumeIntegralY3d, Np3d*K3d*sizeof(double));
	UpdatedVSBEfm3d = malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfm3d, BENe3d * BENfp3d * 3 * sizeof(double));
	UpdatedVSBEzM3d = malloc(BENe3d * BENfp3d * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEzM3d, BENe3d * BENfp3d * sizeof(double));
	UpdatedVSBEfp3d = malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfp3d, BENe3d * BENfp3d * 3 * sizeof(double));
	UpdatedVSBEzP3d = malloc(BENe3d * BENfp3d * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEzP3d, BENe3d * BENfp3d * sizeof(double));
	UpdatedVSBEFluxS3d = malloc(BENe3d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxS3d, BENe3d*BENfp3d*sizeof(double));
	UpdatedVSBEFluxM3d = malloc(BENe3d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxM3d, BENe3d*BENfp3d*sizeof(double));
	UpdatedVSTempFacialIntegral3d = malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempFacialIntegral3d, Np3d*K3d*sizeof(double));
	UpdatedVSBEFluxM3d = malloc(BENe3d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxM3d, BENe3d*BENfp3d*sizeof(double));
	UpdatedVSIEfmod = malloc(IENe2d*IENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfmod, IENe2d*IENfp3d*sizeof(double));
	UpdatedVSBEfmod = malloc(BENe2d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfmod, BENe2d*BENfp3d*sizeof(double));
	Updatedfmod = malloc(K2d*Np3d*sizeof(double));
	MemoryAllocationCheck(Updatedfmod, K2d*Np3d*sizeof(double));
	UpdatedVertVelocityInitialized = "True";
}

void UpdatedVertVelocitySolverMemoryDeAllocation(){
	free(UpdatedVSrhs2d), UpdatedVSrhs2d = NULL;
	free(UpdatedVSIEfm2d), UpdatedVSIEfm2d = NULL;
	free(UpdatedVSIEfp2d), UpdatedVSIEfp2d = NULL;
	free(UpdatedVSIEFluxM2d), UpdatedVSIEFluxM2d = NULL;
	free(UpdatedVSIEFluxP2d), UpdatedVSIEFluxP2d = NULL;
	free(UpdatedVSIEFluxS2d), UpdatedVSIEFluxS2d = NULL;
	free(UpdatedVSERHS2d), UpdatedVSERHS2d = NULL;
	free(UpdatedVSVolumeIntegralX), UpdatedVSVolumeIntegralX = NULL;
	free(UpdatedVSTempVolumeIntegralX), UpdatedVSTempVolumeIntegralX = NULL;
	free(UpdatedVSVolumeIntegralY), UpdatedVSVolumeIntegralY = NULL;
	free(UpdatedVSTempVolumeIntegralY), UpdatedVSTempVolumeIntegralY = NULL;
	free(UpdatedVSBEfm2d), UpdatedVSBEfm2d = NULL;
	free(UpdatedVSBEzM2d), UpdatedVSBEzM2d = NULL;
	free(UpdatedVSBEfp2d), UpdatedVSBEfp2d = NULL;
	free(UpdatedVSBEzP2d), UpdatedVSBEzP2d = NULL;
	free(UpdatedVSBEFluxS2d), UpdatedVSBEFluxS2d = NULL;
	free(UpdatedVSBEFluxM2d), UpdatedVSBEFluxM2d = NULL;
	free(UpdatedVSTempFacialIntegral), UpdatedVSTempFacialIntegral = NULL;
	free(UpdatedVSfield2d), UpdatedVSfield2d = NULL;
	free(UpdatedVSrhs3d), UpdatedVSrhs3d = NULL;
	free(UpdatedVSIEfm3d), UpdatedVSIEfm3d = NULL;
	free(UpdatedVSIEfp3d), UpdatedVSIEfp3d = NULL;
	free(UpdatedVSIEFluxM3d), UpdatedVSIEFluxM3d = NULL;
	free(UpdatedVSIEFluxP3d), UpdatedVSIEFluxP3d = NULL;
	free(UpdatedVSIEFluxS3d), UpdatedVSIEFluxS3d = NULL;
	free(UpdatedVSERHS3d), UpdatedVSERHS3d = NULL;
	free(UpdatedVSVolumeIntegralX3d), UpdatedVSVolumeIntegralX3d = NULL;
	free(UpdatedVSTempVolumeIntegralX3d), UpdatedVSTempVolumeIntegralX3d = NULL;
	free(UpdatedVSVolumeIntegralY3d), UpdatedVSVolumeIntegralY3d = NULL;
	free(UpdatedVSTempVolumeIntegralY3d), UpdatedVSTempVolumeIntegralY3d = NULL;
	free(UpdatedVSBEfm3d), UpdatedVSBEfm3d = NULL;
	free(UpdatedVSBEzM3d), UpdatedVSBEzM3d = NULL;
	free(UpdatedVSBEfp3d), UpdatedVSBEfp3d = NULL;
	free(UpdatedVSBEzP3d), UpdatedVSBEzP3d = NULL;
	free(UpdatedVSBEFluxS3d), UpdatedVSBEFluxS3d = NULL;
	free(UpdatedVSBEFluxM3d), UpdatedVSBEFluxM3d = NULL;
	free(UpdatedVSTempFacialIntegral3d), UpdatedVSTempFacialIntegral3d = NULL;
	free(UpdatedVSIEfmod); UpdatedVSIEfmod = NULL;
	free(UpdatedVSBEfmod); UpdatedVSBEfmod = NULL;
	free(Updatedfmod); Updatedfmod = NULL;
	UpdatedVertVelocityInitialized = "False";
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

/*This is for Updated PCE Solver part*/
double *PCEUpdatedIEfm2d = NULL, *PCEUpdatedIEfp2d = NULL, *PCEUpdatedIEFluxM2d = NULL, *PCEUpdatedIEFluxP2d = NULL, *PCEUpdatedIEFluxS2d = NULL, \
*PCEUpdatedERHS2d = NULL, *PCEUpdatedVolumeIntegralX = NULL, *PCEUpdatedTempVolumeIntegralX = NULL, *PCEUpdatedVolumeIntegralY = NULL, \
*PCEUpdatedTempVolumeIntegralY = NULL, *PCEUpdatedBEfm2d = NULL, *PCEUpdatedBEzM2d = NULL, *PCEUpdatedBEfp2d = NULL, *PCEUpdatedBEzP2d = NULL, \
*PCEUpdatedBEFluxS2d = NULL, *PCEUpdatedBEFluxM2d = NULL, *PCEUpdatedPCETempFacialIntegral = NULL, *PCEUpdatedIEfmod = NULL, *PCEUpdatedBEfmod = NULL, \
*PCEUpdatedIEfm3d = NULL, *PCEUpdatedIEfp3d = NULL, *PCEUpdatedIEFluxS3d = NULL, *PCEUpdatedBEfm3d = NULL, *PCEUpdatedBEfp3d = NULL, *PCEUpdatedBEFluxS3d = NULL, \
*PCEUpdatedBEzM3d = NULL, *PCEUpdatedBEzP3d = NULL;

char *PCEUpdatedInitialized = "False";

void PCEUpdatedMemoryAllocation(int IENfp2d, int IENe2d, int Np2d, int K2d, int Nface, int BENe2d, int BENfp2d, \
	int IENfp3d, int IENe3d, int BENe3d, int BENfp3d){
	PCEUpdatedIEfm2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfm2d, IENfp2d*IENe2d * 3 * sizeof(double));
	PCEUpdatedIEfp2d = malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfp2d, IENfp2d*IENe2d * 3 * sizeof(double));
	PCEUpdatedIEFluxM2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEFluxM2d, IENfp2d*IENe2d*sizeof(double));
	PCEUpdatedIEFluxP2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEFluxP2d, IENfp2d*IENe2d*sizeof(double));
	PCEUpdatedIEFluxS2d = malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEFluxS2d, IENfp2d*IENe2d*sizeof(double));
	PCEUpdatedERHS2d = malloc(Np2d*K2d*Nface*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedERHS2d, Np2d*K2d*Nface*sizeof(double));
	PCEUpdatedVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedVolumeIntegralX, Np2d*K2d*sizeof(double));
	PCEUpdatedTempVolumeIntegralX = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedTempVolumeIntegralX, Np2d*K2d*sizeof(double));
	PCEUpdatedVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedVolumeIntegralY, Np2d*K2d*sizeof(double));
	PCEUpdatedTempVolumeIntegralY = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedTempVolumeIntegralY, Np2d*K2d*sizeof(double));
	PCEUpdatedBEfm2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfm2d, BENe2d * BENfp2d * 3 * sizeof(double));
	PCEUpdatedBEzM2d = malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEzM2d, BENe2d * BENfp2d*sizeof(double));
	PCEUpdatedBEfp2d = malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfp2d, BENe2d * BENfp2d * 3 * sizeof(double));
	PCEUpdatedBEzP2d = malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEzP2d, BENe2d * BENfp2d*sizeof(double));
	PCEUpdatedBEFluxS2d = malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEFluxS2d, BENe2d*BENfp2d*sizeof(double));
	PCEUpdatedBEFluxM2d = malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEFluxM2d, BENe2d*BENfp2d*sizeof(double));
	PCEUpdatedPCETempFacialIntegral = malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedPCETempFacialIntegral, Np2d*K2d*sizeof(double));

	PCEUpdatedIEfmod = malloc(IENe2d*IENfp3d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfmod, IENe2d*IENfp3d*sizeof(double));
	PCEUpdatedBEfmod = malloc(BENe2d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfmod, BENe2d*BENfp3d*sizeof(double));

	PCEUpdatedIEfm3d = malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfm3d, IENfp3d*IENe3d * 3 * sizeof(double));

	PCEUpdatedIEfp3d = malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfp3d, IENfp3d*IENe3d * 3 * sizeof(double));

	PCEUpdatedIEFluxS3d = malloc(IENfp3d*IENe3d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEFluxS3d, IENfp3d*IENe3d*sizeof(double));

	PCEUpdatedBEfm3d = malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfm3d, BENe3d * BENfp3d * 3 * sizeof(double));

	PCEUpdatedBEfp3d = malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfp3d, BENe3d * BENfp3d * 3 * sizeof(double));

	PCEUpdatedBEzM3d = malloc(BENe3d * BENfp3d * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEzM3d, BENe3d * BENfp3d * sizeof(double));

	PCEUpdatedBEzP3d = malloc(BENe3d * BENfp3d * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEzP3d, BENe3d * BENfp3d * sizeof(double));

	PCEUpdatedBEFluxS3d = malloc(BENe3d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEFluxS3d, BENe3d*BENfp3d*sizeof(double));

	PCEUpdatedInitialized = "True";
}

void PCEUpdatedMemoryDeAllocation(){
	free(PCEUpdatedIEfm2d), PCEUpdatedIEfm2d = NULL;
	free(PCEUpdatedIEfp2d), PCEUpdatedIEfp2d = NULL;
	free(PCEUpdatedIEFluxM2d), PCEUpdatedIEFluxM2d = NULL;
	free(PCEUpdatedIEFluxP2d), PCEUpdatedIEFluxP2d = NULL;
	free(PCEUpdatedIEFluxS2d), PCEUpdatedIEFluxS2d = NULL;
	free(PCEUpdatedERHS2d), PCEUpdatedERHS2d = NULL;
	free(PCEUpdatedVolumeIntegralX), PCEUpdatedVolumeIntegralX = NULL;
	free(PCEUpdatedTempVolumeIntegralX), PCEUpdatedTempVolumeIntegralX = NULL;
	free(PCEUpdatedVolumeIntegralY), PCEUpdatedVolumeIntegralY = NULL;
	free(PCEUpdatedTempVolumeIntegralY), PCEUpdatedTempVolumeIntegralY = NULL;
	free(PCEUpdatedBEfm2d), PCEUpdatedBEfm2d = NULL;
	free(PCEUpdatedBEzM2d), PCEUpdatedBEzM2d = NULL;
	free(PCEUpdatedBEfp2d), PCEUpdatedBEfp2d = NULL;
	free(PCEUpdatedBEzP2d), PCEUpdatedBEzP2d = NULL;
	free(PCEUpdatedBEFluxS2d), PCEUpdatedBEFluxS2d = NULL;
	free(PCEUpdatedBEFluxM2d), PCEUpdatedBEFluxM2d = NULL;
	free(PCEUpdatedPCETempFacialIntegral), PCEUpdatedPCETempFacialIntegral = NULL;
	free(PCEUpdatedIEfmod), PCEUpdatedIEfmod = NULL;
	free(PCEUpdatedBEfmod), PCEUpdatedBEfmod = NULL;
	free(PCEUpdatedIEfm3d), PCEUpdatedIEfm3d = NULL;
	free(PCEUpdatedIEfp3d), PCEUpdatedIEfp3d = NULL;
	free(PCEUpdatedBEfm3d), PCEUpdatedBEfm3d = NULL;
	free(PCEUpdatedBEfp3d), PCEUpdatedBEfp3d = NULL;
	free(PCEUpdatedBEzM3d), PCEUpdatedBEzM3d = NULL;
	free(PCEUpdatedBEzP3d), PCEUpdatedBEzP3d = NULL;
	free(PCEUpdatedBEFluxS3d), PCEUpdatedBEFluxS3d = NULL;
	PCEUpdatedInitialized = "False";
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