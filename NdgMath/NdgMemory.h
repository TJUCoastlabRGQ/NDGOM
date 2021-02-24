#ifndef _NdgMemory_H
#define _NdgMemory_H
#include <stdlib.h>
#include <string.h>
typedef enum {
	One = 1,
	Two = 2,
	Three = 3
} NdgMeshType;

void SWENonhydro3dMemoryAllocation(int , int , int , int , int , int , int , int , int , int , \
	int , int , int , int , int , int , int , int , int , int );

void SWENonhydro3dMemoryDeAllocation();

void MemoryAllocationCheck(double *, int );

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

void VertVelocitySolverMemoryAllocation(int, int, int, int, int, int, int, int, \
	int, int, int, int, int);

void VertVelocitySolverMemoryDeAllocation();

#endif
