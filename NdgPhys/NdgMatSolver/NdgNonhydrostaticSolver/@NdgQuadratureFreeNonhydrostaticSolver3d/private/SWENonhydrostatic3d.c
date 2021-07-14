#include "SWENonhydrostatic3d.h"

/*Function checked*/
void AssembleFacialDiffMatrix(double *dest, double *source, double *Eid, int Np2d, int Np3d)
{
	int Index = 0;
	for (int colI = 0; colI < Np3d; colI++){
		for (int RowI = 0; RowI < Np2d; RowI++)
		{
			dest[Index] = source[colI*Np3d + (int)Eid[RowI] - 1];
			Index++;
		}
	}
}

/*We note here that, when the local faicl contribution is assembled to the sparse matrix, AdjEid here refers to the LocalEid, 
* and AdjEle is the local element.
*/
void AssembleFacialContributionIntoSparseMatrix(double *dest, mwIndex *Ir, mwIndex *Jc, double *LocalEid, double *AdjEid, int Np, int Nfp, double *Tempdest, int LocalEle, int AdjEle){
	for (int i = 0; i < Nfp; i++){
		int Sp = (int)Jc[(LocalEle - 1)*Np + (int)LocalEid[i] - 1];
		for (int j = 0; j < (int)Jc[(LocalEle - 1)*Np + (int)LocalEid[i]] - Sp; j++){
			if (Ir[Sp + j] == (mwIndex)((AdjEle - 1)*Np + AdjEid[0] - 1)){
				for (int p = 0; p < Nfp; p++){
					dest[Sp + j + p] += Tempdest[((int)LocalEid[i] - 1)*Np + (int)AdjEid[p] - 1];
				}
			}
			continue;
		}
	}
}

void AssembleVolumnContributionIntoSparseMatrix(double *dest, mwIndex *Ir, mwIndex *Jc, int Np, double *Tempdest, int LocalEle){
	for (int i = 0; i < Np; i++){
		int Sp = (int)Jc[(LocalEle - 1)*Np + i];
		for (int j = 0; j < (int)Jc[(LocalEle - 1)*Np + i + 1] - Sp; j++){
			if (Ir[Sp + j] == (mwIndex)((LocalEle - 1)*Np)){
				for (int p = 0; p < Np; p++){
					dest[Sp + j + p] += Tempdest[i*Np + p];
				}
			}
			continue;
		}
	}
}

void CalculatePenaltyParameter(double *Tau, double *FToE, double *FToN1, double *FToN2, int Np, int Nfp,\
	int Index, double *K13, double *K23, double *K33, double *ELAV, double *MLAV, int P, int Nface){
	
	int LocalEle = (int)FToE[2*Index];
	int AdjEle = (int)FToE[2 * Index + 1];
	double *LocalK13 = K13 + (LocalEle - 1)*Np;
	double *LocalK23 = K23 + (LocalEle - 1)*Np;
	double *LocalK33 = K33 + (LocalEle - 1)*Np;
	double *AdjK13 = K13 + (AdjEle - 1)*Np;
	double *AdjK23 = K23 + (AdjEle - 1)*Np;
	double *AdjK33 = K33 + (AdjEle - 1)*Np;
	double MaximumK = 1.0;
	for (int i = 0; i < Nfp; i++){
		MaximumK = max(MaximumK, abs(LocalK13[(int)FToN1[Index*Nfp + i]-1]));
		MaximumK = max(MaximumK, abs(LocalK23[(int)FToN1[Index*Nfp + i]-1]));
		MaximumK = max(MaximumK, LocalK33[(int)FToN1[Index*Nfp + i]-1]);
		MaximumK = max(MaximumK, abs(AdjK13[(int)FToN2[Index*Nfp + i]-1]));
		MaximumK = max(MaximumK, abs(AdjK23[(int)FToN2[Index*Nfp + i]-1]));
		MaximumK = max(MaximumK, AdjK33[(int)FToN2[Index*Nfp + i]-1]);
	}
	*(Tau + Index) = 10 * 10 * 10 * MaximumK*(P + 1.0)*(P + 3.0) / 3.0*Nface / 2.0*max(ELAV[Index] / MLAV[LocalEle - 1], ELAV[Index] / MLAV[AdjEle - 1]);
}

void EvaluateNonhydroVerticalFaceSurfFlux(double *dest, double *fm, double *n, int Nfp){
	for (int i = 0; i < Nfp; i++){
		dest[i] = fm[i] * n[i];
	}
}

void EvaluateNonhydroVerticalFaceNumFlux_Central(double *dest, double *fm, double *fp, \
	double *n, int Nfp){
	for (int i = 0; i < Nfp; i++){
		dest[i] = (fm[i] + fp[i]) * n[i] * 0.5;
	}
}

/*Function checked*/

void FetchDataInSparseMatrix(double *dest, double *src, int NonzeroNum, int Np){
	double *Tempdest = dest;
	for (int col = 0; col < Np; col++){
		Tempdest = dest + col*Np;
		for (int row = 0; row < Np; row++){
			Tempdest[row] = src[col*NonzeroNum + row];
		}
	}

}

void FetchFacialData(double *dest, double *src, double *FpIndex, int Nfp){

	for (int i = 0; i < Nfp; i++){
		dest[i] = src[(int)FpIndex[i] - 1 ];
	}

}

void FindUniqueElementAndSortOrder(double *dest, double *Src, int *OutNum, int InNum, int LocalElement){
	//First, take the unique value in Src and store them in dest.
	(*OutNum)= 1;
	dest[(*OutNum) - 1] = LocalElement;
	for (int i = 0; i < InNum; i++){
		int Flag;
		int CurrentCell = (int)Src[i];
		for (int j = 0; j < (*OutNum); j++){
			if (CurrentCell == dest[j]){
				Flag = 0;
				break;
			}
			else{
				Flag = 1;
			}
		}
		if (Flag){
			(*OutNum)++;
			dest[(*OutNum) - 1] = CurrentCell;
			}
	}

	Sort(dest, *OutNum);
}

void FindFaceAndDirectionVector(double *FacialVector, int *GlobalFace, \
	int *AdjEle, int *InternalFace, int *Flag, int Nfp, int LocalEle, double *FToE, \
	double *FToF, double *Vector, int IENe, int Nface2d){
	*InternalFace = 0;
	for (int i = 0; i < IENe; i++){
		if ((int)FToE[2 * i] == LocalEle){
			AdjEle[*InternalFace] = (int)FToE[2 * i + 1];
			for (int p = 0; p < Nfp; p++){
				FacialVector[(*InternalFace)*Nfp + p] = Vector[i*Nfp + p];
			}
			GlobalFace[*InternalFace] = i;
			Flag[*InternalFace] = 0;
			(*InternalFace)++;
		}
		if ((int)FToE[2 * i + 1] == LocalEle){
			AdjEle[*InternalFace] = (int)FToE[2 * i];
			for (int p = 0; p < Nfp; p++){
				FacialVector[(*InternalFace)*Nfp + p] = -1 * Vector[i*Nfp + p];
			}
			GlobalFace[*InternalFace] = i;
			Flag[*InternalFace] = 1;
			(*InternalFace)++;
		}
		/*If internal face number equals to Nface2d, we have find all the internal face, just skip the loop*/
		if ((*InternalFace) == Nface2d) break;
	}
}

void FindFaceAndDirectionVectorAtBoundary(double *FacialVector, int *GlobalFace, int *LocalFace, \
	int Nfp, int LocalEle, double *FToE, double *FToF, int face, double *Vector, int BENe){

	for (int f = 0; f < BENe; f++){
		int TempLocalEle = (int)FToE[2 * f];
		int TempLocalFace = (int)FToF[2 * f];
		if (LocalEle == TempLocalEle && TempLocalFace == face){
			(*GlobalFace) = f;
			(*LocalFace) = (int)FToF[2 * f];
			for (int p = 0; p < Nfp; p++){
				FacialVector[p] = Vector[f*Nfp + p];
			}
			break;
		}
	}
}

void FindGlobalBottomEdgeFace(int *GlobalFace, double *FToE, int LocalEle, int AdjacentEle, int Ne){
	for (int f = 0; f < Ne; f++){
		int TempLocalEle = (int)FToE[2 * f];
		int TempAdjEle = (int)FToE[2 * f + 1];
		if (LocalEle == TempLocalEle || LocalEle == TempAdjEle){
			if (AdjacentEle == TempLocalEle || AdjacentEle == TempAdjEle){
				if (AdjacentEle == TempAdjEle && LocalEle == TempLocalEle){
					(*GlobalFace) = f;
					break;
				}
				else if (AdjacentEle == TempLocalEle && LocalEle == TempAdjEle){
					(*GlobalFace) = f;
					break;
				}
				else{
					printf("Problems occured when finding the topological relation for the three dimensional nonhydrostatic model, check again!\n");
					exit(0);
				}
			}
		}
	}
}

/*We note here that Nface2d in the function call can be Nface, when we study the mixed second order derivative in horizontal direction*/
void GetSparsePatternInHorizontalDirection(mwIndex *TempIr, mwIndex *TempJc, double *EToE, double *IEFToE, double *IEFToN1, double *IEFToN2, double *BotEFToE, \
	double *BotEFToN1, double *BotEFToN2, int Nface, int IENfp, int BotENfp, int Np, int Ele3d, int IENe, int BotENe){

	double *TempEToE = malloc((Nface+1)*Ele3d*sizeof(double));
	int *UniNum = malloc(Ele3d*sizeof(int));

	//First we need to know how many elements are adjacent to a studied element.
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < Ele3d; e++){
		FindUniqueElementAndSortOrder(TempEToE + e*(Nface + 1), EToE + e*Nface, UniNum + e, Nface, e+1);
	}
	/*
	* next we'll move on the the construction of the index of the sparse matrix,
	* we first need to know how many nonzeros are contained in each column, and
	* this part can't be parallised since the number of nonzeros of a studied
	* column is influenced by that of the columns come before the studied one.
	*/
	/*The row index, with (Nface+1)*Np the possible maximum influence scope and Np the number of point for each cell*/
	mwIndex *Ir = malloc((Nface + 1)*Np*Np*sizeof(mwIndex));
	/*Number of point influenced by each point*/
	int *NumRowPerPoint = malloc(Np*sizeof(int));
	/*The position the row index of each point has been inserted*/
	int *CurrentPosition = malloc(Np*sizeof(int));

	double *LocalEidM = NULL, *AdjacentEidM = NULL;

	for (int i = 0; i < Ele3d; i++){
		memset(Ir, 0, (Nface + 1)*Np*Np*sizeof(mwIndex));
		memset(NumRowPerPoint, 0, Np*sizeof(int));
		memset(CurrentPosition, 0, Np*sizeof(int));
		for (int e = 0; e < *(UniNum + i);e++){
			int TempEle = (int)TempEToE[i*(Nface + 1) + e];
			/*The cell studied is located upside or downside of cell i*/
			if (TempEle != (i + 1) && (TempEle == EToE[i*Nface + Nface - 2] || TempEle == EToE[i*Nface + Nface - 1])){
				LocalEidM = malloc(BotENfp*sizeof(double));
				AdjacentEidM = malloc(BotENfp*sizeof(double));
				FindFaceAndFacialPoint(LocalEidM, AdjacentEidM, BotENfp, BotEFToE, BotEFToN1, BotEFToN2, BotENe, i + 1, TempEle);
				for (int p = 0; p < BotENfp; p++){
					/*The current point influence the point indexed by AdjacentEidM*/
					NumRowPerPoint[(int)LocalEidM[p] - 1] += BotENfp;
					/*The row index of the point been influenced*/
					for (int j = 0; j < BotENfp; j++){
						Ir[(int)(LocalEidM[p] - 1)*Np*(Nface + 1) + CurrentPosition[(int)(LocalEidM[p] - 1)] + j] = (mwIndex)(AdjacentEidM[j] - 1 + (TempEle - 1)*Np);
					}
				}
				for (int p = 0; p < BotENfp; p++){
					/*Increase the insert position of each point by BotENfp*/
					CurrentPosition[(int)(LocalEidM[p] - 1)] += BotENfp;
				}
				free(LocalEidM), LocalEidM = NULL;
				free(AdjacentEidM), AdjacentEidM = NULL;
			}
			/*The cell studied is located adjacent to cell in horizontal direction*/
			else if (TempEle != (i + 1) && (TempEle != EToE[i*Nface + Nface - 2] && TempEle != EToE[i*Nface + Nface - 1])){
				LocalEidM = malloc(IENfp*sizeof(double));
				AdjacentEidM = malloc(IENfp*sizeof(double));
				FindFaceAndFacialPoint(LocalEidM, AdjacentEidM, IENfp, IEFToE, IEFToN1, IEFToN2, IENe, i+1, TempEle);
				for (int p = 0; p < IENfp; p++){
					NumRowPerPoint[(int)LocalEidM[p] - 1] += IENfp;
					for (int j = 0; j < IENfp; j++){
						Ir[(int)(LocalEidM[p] - 1)*Np*(Nface + 1) + CurrentPosition[(int)(LocalEidM[p] - 1)] + j] = (mwIndex)(AdjacentEidM[j] - 1 + (TempEle - 1)*Np);
					}
				}
				for (int p = 0; p < IENfp; p++){
					CurrentPosition[(int)(LocalEidM[p] - 1)] += IENfp;
				}

				free(LocalEidM), LocalEidM = NULL;
				free(AdjacentEidM), AdjacentEidM = NULL;
			}
			/*The cell studied is cell i itself*/
			else if (TempEle == (i + 1)){
				for (int p = 0; p < Np; p++){
					NumRowPerPoint[p] += Np;
					for (int j = 0; j < Np; j++){
						Ir[p*Np*(Nface + 1) + CurrentPosition[p] + j] = (mwIndex)(j + (TempEle - 1)*Np);
					}
				}
				for (int p = 0; p < Np; p++){
					CurrentPosition[p] += Np;
				}
			}
		}
		for (int p = 0; p < Np; p++){
			TempJc[i*Np + p + 1] = TempJc[i*Np + p] + NumRowPerPoint[p];
		}
		for (int col = 0; col < Np; col++){
			for (int p = 0; p < NumRowPerPoint[col]; p++){
				/*The row index in the sparse matrix*/
				TempIr[TempJc[i*Np + col] + p] = Ir[col * (Nface + 1) * Np + p];
			}
		}

	}
	free(TempEToE);
	free(UniNum);
	free(Ir);
	free(NumRowPerPoint);
	free(CurrentPosition);
}

void FindFaceAndFacialPoint(double *Localdest, double *Adjacentdest, int Nfp, double *FToE, double *FToN1, double *FToN2, int Ne, int LocalEle, int AdjEle){
	for (int i = 0; i < Ne; i++){
		if ((int)FToE[2 * i] == LocalEle && (int)FToE[2 * i + 1] == AdjEle){
			for (int p = 0; p < Nfp; p++){
				Localdest[p] = FToN1[i*Nfp + p];
				Adjacentdest[p] = FToN2[i*Nfp + p];
			}
			Sort(Localdest, Nfp);
			Sort(Adjacentdest, Nfp);
			break;
		}
		else if ((int)FToE[2 * i + 1] == LocalEle && (int)FToE[2 * i] == AdjEle){
			for (int p = 0; p < Nfp; p++){
				Localdest[p] = FToN2[i*Nfp + p];
				Adjacentdest[p] = FToN1[i*Nfp + p];
			}
			Sort(Localdest, Nfp);
			Sort(Adjacentdest, Nfp);
			break;
		}
	}
}
