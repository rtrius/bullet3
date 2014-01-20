/*
BulletFluids 
Copyright (c) 2013 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef B3_FLUID_HASH_GRID_OPENCL_H
#define B3_FLUID_HASH_GRID_OPENCL_H

#include "Bullet3Common/b3Scalar.h"

#include "Bullet3OpenCL/ParallelPrimitives/b3RadixSort32CL.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3PrefixScanCL.h"

#include "Bullet3Fluids/Sph/b3FluidSph.h"
class b3FluidSphOpenCL;


#define B3_FLUID_HASH_GRID_COORD_RANGE 64
#define B3_FLUID_HASH_GRID_NUM_CELLS B3_FLUID_HASH_GRID_COORD_RANGE * B3_FLUID_HASH_GRID_COORD_RANGE * B3_FLUID_HASH_GRID_COORD_RANGE

class b3FluidHashGridOpenCL
{
public:	
	b3OpenCLArray<b3FluidGridIterator> m_cellContents;
	
	b3FluidHashGridOpenCL(cl_context context, cl_command_queue queue);
};

class b3FluidHashGridOpenCLProgram
{
	cl_program m_hashGridProgram;
	cl_kernel m_generateValueIndexPairsModuloKernel;
	cl_kernel m_rearrangeParticleArraysKernel;
	
	cl_kernel m_resetGridCellsModuloKernel;
	cl_kernel m_detectIndexRangesModuloKernel;

	b3OpenCLArray<b3Vector3> m_tempBufferCL;		//Used to rearrange fluid particle arrays(position, velocity, etc.)
	//b3AlignedObjectArray<b3Vector3> m_tempBufferVector;
	//b3AlignedObjectArray<void*> m_tempBufferVoid;
	
	b3RadixSort32CL m_radixSorter;
	b3OpenCLArray<b3SortData> m_valueIndexPairs;
	//b3AlignedObjectArray<b3SortData> m_valueIndexPairsHost;
	
public:
	b3FluidHashGridOpenCLProgram(cl_context context, cl_device_id device, cl_command_queue queue);
	~b3FluidHashGridOpenCLProgram();
	
	void insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
								 b3FluidSph* fluid, b3FluidSphOpenCL* fluidData, b3FluidHashGridOpenCL* gridData);
	
	
private:
	void generateValueIndexPairsModulo(cl_command_queue commandQueue, int numFluidParticles, b3Scalar cellSize, cl_mem fluidPositionsBuffer);
	void rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, cl_mem fluidBuffer);
};

#endif