/*
BulletFluids 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef B3_FLUID_SORTING_GRID_OPENCL_H
#define B3_FLUID_SORTING_GRID_OPENCL_H

#include "Bullet3Common/b3Scalar.h"

#include "Bullet3OpenCL/ParallelPrimitives/b3RadixSort32CL.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3PrefixScanCL.h"

#include "Bullet3Fluids/Sph/b3FluidSphTypedData.h"
#include "Bullet3Fluids/Sph/b3FluidSph.h"
#include "Bullet3Fluids/Sph/b3FluidSortingGrid.h"

class b3FluidSortingGrid;
class b3FluidSphOpenCL;

///Manages OpenCL buffers corresponding to a b3FluidSortingGrid.
class b3FluidSortingGridOpenCL : public b3FluidSphTypedData
{
public:	
	b3OpenCLArray<int> m_numActiveCells;
	
	b3OpenCLArray<b3FluidGridCombinedPos> m_activeCells;
	b3OpenCLArray<b3FluidGridIterator> m_cellContents;
	b3OpenCLArray<b3FluidSortingGrid::FoundCellsGpu> m_foundCells;
	
	//Size == m_numActiveCells
		//b3SortData.m_key == grid cell value(value to sort by), converted to y/z-axis orientation
		//b3SortData.m_value == x-axis oriented grid cell index
		b3OpenCLArray<b3SortData> m_yOrientedPairs;
		b3OpenCLArray<b3SortData> m_zOrientedPairs;
		
		b3OpenCLArray<int> m_yIndex;	//Parallel array to m_activeCells
		b3OpenCLArray<int> m_zIndex;	//Parallel array to m_activeCells
	
	b3FluidSortingGridOpenCL(cl_context context, cl_command_queue queue);
		
	virtual b3FluidSphDataType getType() const { return FSDT_b3FluidSortingGridOpenCL; }
	
	void writeToOpenCL(cl_command_queue queue, b3FluidSortingGrid& sortingGrid);
	void readFromOpenCL(cl_command_queue queue, b3FluidSortingGrid& sortingGrid);
	
	int getNumActiveCells() const;
};

///Parallelized implementation of b3FluidSortingGridOpenCLProgram::generateUniques_serial() for OpenCL.
class b3FluidSortingGridOpenCLProgram_GenerateUniques
{
	cl_program m_sortingGridProgram;
	cl_kernel m_markUniquesKernel;
	cl_kernel m_storeUniquesAndIndexRangesKernel;
	
	cl_kernel m_convertCellValuesAndLoadCellIndexKernel;
	cl_kernel m_writebackReorientedCellIndiciesKernel;
	
	b3PrefixScanCL m_prefixScanner;
	
	b3OpenCLArray<unsigned int> m_tempInts;
	b3OpenCLArray<unsigned int> m_scanResults;
	
public:
	b3FluidSortingGridOpenCLProgram_GenerateUniques(cl_context context, cl_device_id device, cl_command_queue queue);
	~b3FluidSortingGridOpenCLProgram_GenerateUniques();
	
	void generateUniques(cl_command_queue commandQueue, const b3OpenCLArray<b3SortData>& valueIndexPairs,
							b3FluidSortingGridOpenCL* gridData, int numFluidParticles, b3RadixSort32CL& radixSorter);
};

///Implements b3FluidSortingGrid::insertParticles() for OpenCL.
///@remarks
///#define B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT is not supported.
class b3FluidSortingGridOpenCLProgram
{
	cl_program m_sortingGridProgram;
	cl_kernel m_generateValueIndexPairsKernel;
	cl_kernel m_rearrangeParticleArraysKernel;
	cl_kernel m_generateUniquesKernel;

	b3OpenCLArray<b3Vector3> m_tempBufferCL;		//Used to rearrange fluid particle arrays(position, velocity, etc.)
	b3AlignedObjectArray<b3Vector3> m_tempBufferVector;
	b3AlignedObjectArray<void*> m_tempBufferVoid;
	
	b3RadixSort32CL m_radixSorter;
	b3OpenCLArray<b3SortData> m_valueIndexPairs;
	b3AlignedObjectArray<b3SortData> m_valueIndexPairsHost;
	
	b3FluidSortingGridOpenCLProgram_GenerateUniques m_generateUniquesProgram;
	
public:
	b3FluidSortingGridOpenCLProgram(cl_context context, cl_device_id device, cl_command_queue queue);
	~b3FluidSortingGridOpenCLProgram();
	
	void insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
								 b3FluidSph* fluid, b3FluidSphOpenCL* fluidData, b3FluidSortingGridOpenCL* gridData);
	
	//This can only be called after insertParticlesIntoGrid() is for the current fluid
	//and before insertParticlesIntoGrid() is called for the next fluid
	//void rearrangeParticlesOnHost(b3FluidSph* fluid);
	
private:
	void generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, b3Scalar cellSize, cl_mem fluidPositionsBuffer);
	void rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, cl_mem fluidBuffer);
	void generateUniques_serial(cl_command_queue commandQueue, int numFluidParticles, b3FluidSortingGridOpenCL* gridData);
};

#endif