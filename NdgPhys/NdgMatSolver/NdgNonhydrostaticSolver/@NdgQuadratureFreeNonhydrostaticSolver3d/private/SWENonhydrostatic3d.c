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

/*Function checked*/
void AssembleContributionIntoSparseMatrix(double *dest, double *src, int NonzeroNum, int Np){
	double *Tempdest = dest;
	for (int col = 0; col < Np; col++){
		Tempdest = dest + col * NonzeroNum;
		for (int row = 0; row < Np; row++){
			Tempdest[row] += src[col*Np + row];
		}
	}
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
	// Adopt bubble method to sort the given array.
	int isSorted;
	double temp;
	for (int i = 0; i<(*OutNum) - 1; i++){
		//Suppose the rest data are already sorted
		isSorted = 1;
		for (int j = 0; j<(*OutNum) - 1 - i; j++){
			if (dest[j] > dest[j + 1]){
				temp = dest[j];
				dest[j] = dest[j + 1];
				dest[j + 1] = temp;
				isSorted = 0;  //Once swapped, it indicates the data is not sorted
			}
		}
		if (isSorted) break; //if no swap needed, it indicates the data is already sorted.
	}
}

void FindFaceAndDirectionVector(double *FacialVector, int *GlobalFace, int *LocalFace, int *AdjacentFace, \
	int Nfp, int LocalEle, int AdjacentEle, double *FToE, double *FToF, double *Vector, int IENe){

	for (int f = 0; f < IENe; f++){
		int TempLocalEle = (int)FToE[2 * f];
		int TempAdjEle = (int)FToE[2 * f + 1];
		if (LocalEle == TempLocalEle || LocalEle == TempAdjEle){
			if (AdjacentEle == TempLocalEle || AdjacentEle == TempAdjEle){
				if (AdjacentEle == TempAdjEle && LocalEle == TempLocalEle){
					(*GlobalFace) = f;
					(*LocalFace) = (int)FToF[2 * f];
					(*AdjacentFace) = (int)FToF[2 * f + 1];
					for (int p = 0; p < Nfp; p++){
						FacialVector[p] = Vector[f*Nfp + p];
					}
					break;
				}
				else if (AdjacentEle == TempLocalEle && LocalEle == TempAdjEle){
					(*GlobalFace) = f;
					(*LocalFace) = (int)FToF[2 * f + 1];
					(*AdjacentFace) = (int)FToF[2 * f];
					for (int p = 0; p < Nfp; p++){
						FacialVector[p] = -1 * Vector[f*Nfp + p];
					}
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

void GetSparsePatternInVerticalDirection(mwIndex *TempIr, mwIndex *TempJc, int Np, int Nlayer, int Ele2d){
	int *SingleColumn = malloc(Np*Nlayer*sizeof(int));
	int *SingleRow;
	int SingleNonzero;
	if (Nlayer == 1){
		SingleNonzero = Np*Np;
		SingleRow = malloc(SingleNonzero*sizeof(int));
		for (int ic = 0; ic < Np; ic++){
			SingleColumn[ic] = (ic + 1)*Np;
			for (int jr = 0; jr < Np; jr++){
				SingleRow[ic*Np + jr] = jr;
			}
		}
	}
	else if (Nlayer == 2){
		SingleNonzero = Np*Np * 2 * 2;
		SingleRow = malloc(SingleNonzero*sizeof(int));
		for (int ic = 0; ic < 2 * Np; ic++){
			SingleColumn[ic] = (ic + 1)*Np * 2;
			for (int jr = 0; jr < 2 * Np; jr++){
				SingleRow[ic*Np * 2 + jr] = jr;
			}
		}
	}
	else{
		SingleNonzero = Np*Np * 2 * 2 + (Nlayer - 2)*Np*Np * 3;
		SingleRow = malloc(SingleNonzero*sizeof(int));
		int NSumR = 0;
		for (int Layer = 0; Layer < 1; Layer++){
			for (int ic = Layer*Np; ic < (Layer + 1)*Np; ic++){
				SingleColumn[ic] = (ic + 1) * 2 * Np;
				for (int jr = 0; jr < 2 * Np; jr++){
					SingleRow[NSumR + ic*Np * 2 + jr] = jr;
				}
			}
		}
		NSumR = Np * Np * 2;
		for (int Layer = 1; Layer < Nlayer - 1; Layer++){
			for (int ic = Layer*Np; ic < (Layer + 1)*Np; ic++){
				SingleColumn[ic] = SingleColumn[ic - 1] + 3 * Np;
				for (int jr = 0; jr < 3 * Np; jr++){
					SingleRow[NSumR + (ic - Np) * 3 * Np + jr] = (Layer - 1)*Np + jr;
				}
			}
		}
		NSumR = Np*Np * 2 + (Nlayer - 2)*Np * 3 * Np;
		for (int Layer = Nlayer - 1; Layer < Nlayer; Layer++){
			for (int ic = Layer*Np; ic < (Layer + 1)*Np; ic++){
				SingleColumn[ic] = SingleColumn[ic - 1] + 2 * Np;
				for (int jr = 0; jr < 2 * Np; jr++){
					SingleRow[NSumR + (ic - (Nlayer - 1)*Np) * 2 * Np + jr] = (Layer - 1)*Np + jr;   //测试
				}
			}
		}
	}

	//接下来对稀疏矩阵的索引坐标进行复制。
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < Ele2d; i++){

		for (int j = 0; j < SingleNonzero; j++){
			TempIr[i*SingleNonzero + j] = SingleRow[j] + i*Nlayer*Np;
		}
		for (int j = 0; j < Nlayer*Np; j++){
			TempJc[1 + i*Nlayer*Np + j] = i * SingleColumn[Nlayer*Np - 1] + SingleColumn[j];
		}
	}

	free(SingleColumn);
	free(SingleRow);
}


/*We note here that Nface2d in the function call can be Nface, when we study the mixed second order derivative in horizontal direction*/
void GetSparsePatternInHorizontalDirection(mwIndex *TempIr, mwIndex *TempJc, double *EToE, int Nface2d, int Nface, int Ele3d, int Np){
	double *TempEToE = malloc((Nface2d+1)*Ele3d*sizeof(double));
	int *UniNum = malloc(Ele3d*sizeof(int));

	//First we need to know how many elements are adjacent to a studied element.
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < Ele3d; e++){
		FindUniqueElementAndSortOrder(TempEToE + e*(Nface2d + 1), EToE + e*Nface, UniNum + e, Nface2d, e+1);
	}
	/*
	* next we'll move on the the construction of the index of the sparse matrix,
	* we first need to know how many nonzeros are contained in each column, and
	* this part can't be parallised since the number of nonzeros of a studied
	* column is influenced by that of the columns come before the studied one.
	*/
	
	for (int i = 0; i < Ele3d; i++){
		int NumRow = *(UniNum + i)*Np;
		for (int j = 0; j < Np; j++){
			TempJc[i*Np + j + 1] = TempJc[i*Np + j] + NumRow;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < Ele3d; i++){
		for (int e = 0; e < *(UniNum + i); e++){
			int ele = (int)TempEToE[i*(Nface2d+1) + e];
			for (int col = 0; col < Np; col++){
				for (int row = 0; row < Np; row++){
					TempIr[(int)TempJc[i*Np + col] + e*Np + row] = row + (ele - 1)*Np;
				}
			}
		}
	}
	free(TempEToE);
	free(UniNum);
}