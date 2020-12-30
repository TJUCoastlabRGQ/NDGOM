# include "NdgMemory.h"
#include <stdlib.h>

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
	HorDiffzM = malloc(BENfp*BENe*sizeof(double));
	HorDiffzP = malloc(BENfp*BENe*sizeof(double));
	if (type == Two){
		/*For 2d shallow water problem, no horizontal diffusion terms are included in the governing equation for water depth $H$*/
		Nfield = Nvar - 1;
		/*For 2d shallow water problem, the face number is equal to TempNface, since there is no surface edge and bottom edge*/
		Nface = tempNface;
		/*Allocate memory for the original variable $u,v$ and $\theta$*/
		HorDiffvariable = malloc(Np*K*Nfield*sizeof(double));
		/*Allocate memory for the original variable over boundary edge*/
		HorDiffBEfp = malloc(BENfp*BENe*Nfield*sizeof(double));
		/*Allocate memory for the local face value at boundary edge*/
		HorDiffTempBEfp = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
		HorDiffTempBEfm = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
	}
	else if (type == Three){
		Nfield = Nvar;
		/*For 3d shallow water problem, the face number is equal to TempNface - 2, since there
		* is surface edge and bottom edge is not considered for horizontal diffusion term*/
		Nface = tempNface - 2;
		/*Allocate memory for the original variable $u,v$ and $\theta$*/
		HorDiffvariable = malloc(Np*K*Nfield*sizeof(double));

		HorDiffTempBEfp = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
		HorDiffTempBEfm = malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
		/*Allocate memory for the original variable over boundary edge*/
		HorDiffBEfp = malloc(BENfp*BENe*Nfield*sizeof(double));
	}
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
	/*Allocate memory for the local face value at boundary edge*/
	HorDiffBEfm = malloc(BENe*BENfp*Nfield*sizeof(double));
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
	free(HorDiffnv);
	free(HorDiffvariable);
	free(HorDiffBEfp); 
	free(HorDiffzM);
	free(HorDiffzP);
	free(HorDiffTempBEfp);
	free(HorDiffTempBEfm);
	free(HorDiffAVx);
	free(HorDiffAVy);
	free(HorDiffVx);
	free(HorDiffTempVx);
	free(HorDiffVy);
	free(HorDiffTempVy);
	free(HorDiffIEfm);
	free(HorDiffAVIEfm);
	free(HorDiffIEfp);
	free(HorDiffAVIEfp);
	free(HorDiffIEFluxM);
	free(HorDiffIEFluxP);
	free(HorDiffBEfm);
	free(HorDiffIEFluxS);
	free(HorDiffAVBEfm);
	free(HorDiffBEFluxM);
	free(HorDiffBEFluxS);
	free(HorDiffERHSX);
	free(HorDiffERHSY);
	free(HorDiffLocalPrimitiveDiffTermX);
	free(HorDiffLocalPrimitiveDiffTermY);
	free(HorDiffLPDTIEfm);
	free(HorDiffLPDTIEfp);
	free(HorDiffLPDTBEfm);
	free(HorDiffTempFacialIntegralX);
	free(HorDiffTempFacialIntegralY);
	free(HorDiffInnerEdgeTau);
	free(HorDiffBoundaryEdgeTau);
	free(HorDiffIEnvfm);
	free(HorDiffIEnvfp);
	free(HorDiffBEnvfm);

	HorDiffInitialized = "True";
}