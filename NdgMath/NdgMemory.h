#ifndef _NdgMemory_H
#define _NdgMemory_H
#include <stdlib.h>
#include <string.h>
typedef enum {
	One = 1,
	Two = 2,
	Three = 3
} NdgMeshType;

void BaroclinicPartMemoryAllocation(int , int , int , int , int , int , int , int, int, int);

void BaroclinicPartMemoryDeAllocation();

void BaroDensityMemoryAllocation(int, int);

void BaroDensityMemoryDeAllocation();

void SWENH3dTimeIntervalMemoryAllocation(int, int, int);

void SWENH3dTimeIntervalMemoryDeAllocation();

void SWENH3dImposeBoundaryMemoryAllocation(int, int, int, int, int, int, int, int, int, int);

void SWENH3dImposeBoundaryMemoryDeAllocation();

void SWENonhydro3dMemoryAllocation(int , int , int , int , int , int , int , int , int , int , \
	int , int , int , int , int , int , int , int , int , int );

void SWENonhydro3dMemoryDeAllocation();

void SWENonhydroVertVelocityMemoryAllocation(int, int);

void SWENonhydroVertVelocityMemoryDeAllocation();

void GlobalStiffMatrixMemoryAllocation(int, int, int, int, int, int, int, int, int, int);

void GlobalStiffMatrixMemoryDeAllocation();

void MemoryAllocationCheck(void *, int );

void VertDiffMemoryAllocation(const int, int, const int);

void VertDiffMemoryDeAllocation();

void GotmSolverMemoryAllocation(int, int, int, int);

void GotmSolverMemoryDeAllocation();

void HorizDiffMemoryAllocation(NdgMeshType , int , int , int , int , int , int , int , int );

void HorizDiffMemoryDeAllocation();

void AdvMemoryAllocation(int, int, int, int, int, int, int, int, int, int, int, int, int, int);

void AdvMemoryDeAllocation();

void PCEMemoryAllocation(int ,int ,int ,int ,int ,int ,int );

void PCEMemoryDeAllocation();

void PCEUpdatedMemoryAllocation(int, int, int, int, int, int, int, \
	int, int, int, int);

void PCEUpdatedMemoryDeAllocation();

void VertVelocitySolverMemoryAllocation(int, int, int, int, int, int, int, int, \
	int, int, int, int, int);

void UpdatedVertVelocitySolverMemoryAllocation(int, int, int, int, int, int, int, int, \
	int, int, int, int, int);

void UpdatedVertVelocitySolverMemoryDeAllocation();

void VertVelocitySolverMemoryDeAllocation();

void ImVertDiffMemoryAllocation(const int , int , const int , int , int );

void ImVertDiffMemoryDeAllocation();

#endif
