#ifndef _NdgMemory_H
#define _NdgMemory_H
#include <stdlib.h>

typedef enum {
	One = 1,
	Two = 2,
	Three = 3
} NdgMeshType;

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
