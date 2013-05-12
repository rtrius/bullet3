/*
Bullet-FLUIDS 
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
#include "b3FluidSortingGridOpenCL.h"

#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"

#include "BulletFluids/Sph/b3FluidSortingGrid.h"

#include "b3FluidSphOpenCL.h"
#include "fluidSphCL.h"


#define SWAP_REARRANGED_ARRAY
#ifdef SWAP_REARRANGED_ARRAY
#include <cstring>	//memcpy()
#endif

// /////////////////////////////////////////////////////////////////////////////
//class b3FluidSortingGridOpenCL
// /////////////////////////////////////////////////////////////////////////////
b3FluidSortingGridOpenCL::b3FluidSortingGridOpenCL(cl_context context, cl_command_queue queue) 
: m_numActiveCells(context, queue), m_activeCells(context, queue), m_foundCells(context, queue), m_cellContents(context, queue) 
{
	m_numActiveCells.resize(1);
}

void b3FluidSortingGridOpenCL::writeToOpenCL(cl_command_queue queue, b3FluidSortingGrid& sortingGrid)
{
	int numActiveCells = sortingGrid.internalGetActiveCells().size();
	
	m_numActiveCells.copyFromHostPointer(&numActiveCells, 1, 0, false);
	
	m_activeCells.copyFromHost( sortingGrid.internalGetActiveCells(), false );
	m_cellContents.copyFromHost( sortingGrid.internalGetCellContents(), false );
	
	clFinish(queue);
}
void b3FluidSortingGridOpenCL::readFromOpenCL(cl_command_queue queue, b3FluidSortingGrid& sortingGrid)
{
	m_activeCells.copyToHost( sortingGrid.internalGetActiveCells(), false );
	m_cellContents.copyToHost( sortingGrid.internalGetCellContents(), false );
	
	clFinish(queue);
}

int b3FluidSortingGridOpenCL::getNumActiveCells() const
{
	int numActiveCells;
	m_numActiveCells.copyToHostPointer(&numActiveCells, 1, 0, true);
	
	return numActiveCells;
}


// /////////////////////////////////////////////////////////////////////////////
//class b3FluidSortingGridOpenCLProgram_GenerateUniques
// /////////////////////////////////////////////////////////////////////////////
b3FluidSortingGridOpenCLProgram_GenerateUniques::b3FluidSortingGridOpenCLProgram_GenerateUniques
(
	cl_context context, 
	cl_device_id device,
	cl_command_queue queue
)
: m_prefixScanner(context, device, queue), m_tempInts(context, queue), m_scanResults(context, queue)
{
	const char CL_SORTING_GRID_PROGRAM_PATH[] = "src/BulletFluidsOpenCL/fluidSph.cl";
	
	const char* kernelSource = fluidSphCL;	//fluidSphCL.h
	cl_int error;
	char* additionalMacros = 0;
	sortingGrid_program = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, 
																	additionalMacros, CL_SORTING_GRID_PROGRAM_PATH);
	b3Assert(sortingGrid_program);
	
	kernel_markUniques = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "markUniques", &error, sortingGrid_program, additionalMacros );
	b3Assert(kernel_markUniques);
	kernel_storeUniques = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "storeUniques", &error, sortingGrid_program, additionalMacros );
	b3Assert(kernel_storeUniques);
	kernel_setZero = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "setZero", &error, sortingGrid_program, additionalMacros );
	b3Assert(kernel_setZero);
	kernel_countUniques = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "countUniques", &error, sortingGrid_program, additionalMacros );
	b3Assert(kernel_countUniques);
	kernel_generateIndexRanges = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "generateIndexRanges", &error, sortingGrid_program, additionalMacros );
	b3Assert(kernel_generateIndexRanges);
}
b3FluidSortingGridOpenCLProgram_GenerateUniques::~b3FluidSortingGridOpenCLProgram_GenerateUniques()
{
	clReleaseKernel(kernel_markUniques);
	clReleaseKernel(kernel_storeUniques);
	clReleaseKernel(kernel_setZero);
	clReleaseKernel(kernel_countUniques);
	clReleaseKernel(kernel_generateIndexRanges);
	clReleaseProgram(sortingGrid_program);
}

void b3FluidSortingGridOpenCLProgram_GenerateUniques::generateUniques(cl_command_queue commandQueue,
																	const b3OpenCLArray<b3SortData>& valueIndexPairs,
																	b3FluidSortingGridOpenCL* gridData, int numFluidParticles)
{
	if( m_tempInts.size() < numFluidParticles ) m_tempInts.resize(numFluidParticles, false);
	if( m_scanResults.size() < numFluidParticles ) m_scanResults.resize(numFluidParticles, false);
	
	b3OpenCLArray<int>& out_numActiveCells = gridData->m_numActiveCells;
	b3OpenCLArray<b3FluidGridCombinedPos>& out_sortGridValues = gridData->m_activeCells;
	b3OpenCLArray<b3FluidGridIterator>& out_iterators = gridData->m_cellContents;
	
	unsigned int numUniques = 0;
	
	//Detect unique values
	{
		//If the element to the right is different(or out of bounds), set 1; set 0 otherwise
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( valueIndexPairs.getBufferCL() ),
				b3BufferInfoCL( m_tempInts.getBufferCL() )
			};
	
			b3LauncherCL launcher(commandQueue, kernel_markUniques);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
		
		//
		{
			m_prefixScanner.execute(m_tempInts, m_scanResults, numFluidParticles, &numUniques);		//Exclusive scan
			++numUniques;	//Prefix scanner returns last index if the array is filled with 1 and 0; add 1 to get size
			
			int numActiveCells = static_cast<int>(numUniques);
			out_numActiveCells.copyFromHostPointer(&numActiveCells, 1, 0, true);
			
			out_sortGridValues.resize(numActiveCells, false);
			out_iterators.resize(numActiveCells, false);
		}
		
		//Use scan results to store unique b3FluidGridCombinedPos
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( valueIndexPairs.getBufferCL() ),
				b3BufferInfoCL( m_tempInts.getBufferCL() ),
				b3BufferInfoCL( m_scanResults.getBufferCL() ),
				b3BufferInfoCL( out_sortGridValues.getBufferCL() )
			};
			
			b3LauncherCL launcher(commandQueue, kernel_storeUniques);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
	}
	
	//Count number of each unique value and use count to generate index ranges
	{
		//Set m_tempInts to 0; it will store the number of particles corresponding to each b3FluidGridCombinedPos(num particles in each cell)
		{
			b3BufferInfoCL bufferInfo[] = { b3BufferInfoCL( m_tempInts.getBufferCL() ) };
		
			b3LauncherCL launcher(commandQueue, kernel_setZero);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numUniques);
			
			launcher.launch1D(numUniques);
		}
		
		//Use binary search and atomic increment to count the number of particles in each cell
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( valueIndexPairs.getBufferCL() ),
				b3BufferInfoCL( out_sortGridValues.getBufferCL() ),
				b3BufferInfoCL( m_tempInts.getBufferCL() )
			};
			
			b3LauncherCL launcher(commandQueue, kernel_countUniques);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numUniques);
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
		
		//Exclusive scan
		m_prefixScanner.execute(m_tempInts, m_scanResults, numUniques, 0);	
		
		//Use scan results to generate index ranges(check element to the right)
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( m_scanResults.getBufferCL() ),
				b3BufferInfoCL( out_iterators.getBufferCL() )
			};
			
			b3LauncherCL launcher(commandQueue, kernel_generateIndexRanges);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst( static_cast<int>(numUniques) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numUniques);
		}
	}
}
	
// /////////////////////////////////////////////////////////////////////////////
//class b3FluidSortingGridOpenCLProgram
// /////////////////////////////////////////////////////////////////////////////
b3FluidSortingGridOpenCLProgram::b3FluidSortingGridOpenCLProgram(cl_context context, cl_command_queue queue, cl_device_id device)
: m_tempBufferCL(context, queue), m_radixSorter(context, device, queue), 
	m_valueIndexPairs(context, queue), m_generateUniquesProgram(context, device, queue)
{
	const char CL_SORTING_GRID_PROGRAM_PATH[] = "src/BulletFluidsOpenCL/fluidSph.cl";
	
	const char* kernelSource = fluidSphCL;	//fluidSphCL.h
	cl_int error;
	char* additionalMacros = 0;
	sortingGrid_program = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, 
																	additionalMacros, CL_SORTING_GRID_PROGRAM_PATH);
	b3Assert(sortingGrid_program);

	kernel_generateValueIndexPairs = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "generateValueIndexPairs", &error, sortingGrid_program, additionalMacros );
	b3Assert(kernel_generateValueIndexPairs);
	kernel_rearrangeParticleArrays = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "rearrangeParticleArrays", &error, sortingGrid_program, additionalMacros );
	b3Assert(kernel_rearrangeParticleArrays);
	kernel_generateUniques = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "generateUniques", &error, sortingGrid_program, additionalMacros );
	b3Assert(kernel_generateUniques);
}
b3FluidSortingGridOpenCLProgram::~b3FluidSortingGridOpenCLProgram()
{
	clReleaseKernel(kernel_generateValueIndexPairs);
	clReleaseKernel(kernel_rearrangeParticleArrays);
	clReleaseKernel(kernel_generateUniques);
	clReleaseProgram(sortingGrid_program);
}

template<typename T>
void rearrangeToMatchSortedValues2(const b3AlignedObjectArray<b3SortData>& sortedValues, 
									b3AlignedObjectArray<T>& temp, b3AlignedObjectArray<T>& out_rearranged)
{
	temp.resize( out_rearranged.size() );
	
	for(int i = 0; i < out_rearranged.size(); ++i)
	{
		int oldIndex = sortedValues[i].m_value;
		int newIndex = i;
			
		temp[newIndex] = out_rearranged[oldIndex];
	}

#ifndef SWAP_REARRANGED_ARRAY
	out_rearranged = temp;
#else
	const int SIZEOF_ARRAY = sizeof(b3AlignedObjectArray<T>);
	
	char swap[SIZEOF_ARRAY];
		
	memcpy(swap, &temp, SIZEOF_ARRAY);
	memcpy(&temp, &out_rearranged, SIZEOF_ARRAY);
	memcpy(&out_rearranged, swap, SIZEOF_ARRAY);
#endif
}
void b3FluidSortingGridOpenCLProgram::insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
															  b3FluidSph* fluid, b3FluidSphOpenCL* fluidData, b3FluidSortingGridOpenCL* gridData)
{
#ifdef B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT
	b3Assert(0);	//Not implemented
#endif

	B3_PROFILE("b3FluidSortingGridOpenCLProgram::insertParticlesIntoGrid()");
	
	int numFluidParticles = fluid->numParticles();
	m_tempBufferCL.resize(numFluidParticles);
	m_valueIndexPairs.resize(numFluidParticles);
	
	//Cannot check number of nonempty grid cells before generateUniques();
	//temporarily resize m_activeCells and m_cellContents
	//to handle the case where each particle occupies a different grid cell.
	gridData->m_activeCells.resize(numFluidParticles);
	gridData->m_cellContents.resize(numFluidParticles);
	
	//
	{
		B3_PROFILE("generateValueIndexPairs()");
		generateValueIndexPairs( commandQueue, numFluidParticles, fluid->getGrid().getCellSize(), fluidData->m_pos.getBufferCL() );
		
		clFinish(commandQueue);
	}
	
	//Note that b3RadixSort32CL uses b3SortData, while b3FluidSortingGrid uses b3FluidGridValueIndexPair.
	//b3SortData.m_key == b3FluidGridValueIndexPair.m_value (value to sort by)
	//b3SortData.m_value == b3FluidGridValueIndexPair.m_index (fluid particle index)
	{
		B3_PROFILE("radix sort");
		m_radixSorter.execute(m_valueIndexPairs, 32);
		
		clFinish(commandQueue);
	}
	
	
	//
	{
		B3_PROFILE("rearrange device");
		
		rearrangeParticleArrays( commandQueue, numFluidParticles, fluidData->m_pos.getBufferCL() );
		fluidData->m_pos.copyFromOpenCLArray(m_tempBufferCL);
		
		rearrangeParticleArrays( commandQueue, numFluidParticles, fluidData->m_vel.getBufferCL() );
		fluidData->m_vel.copyFromOpenCLArray(m_tempBufferCL);
		
		rearrangeParticleArrays( commandQueue, numFluidParticles, fluidData->m_vel_eval.getBufferCL() );
		fluidData->m_vel_eval.copyFromOpenCLArray(m_tempBufferCL);
		
		rearrangeParticleArrays( commandQueue, numFluidParticles, fluidData->m_accumulatedForce.getBufferCL() );
		fluidData->m_accumulatedForce.copyFromOpenCLArray(m_tempBufferCL);
		
		clFinish(commandQueue);
	}
	
	//
	const bool USE_PARALLEL_GENERATE_UNIQUES = true;
	if(USE_PARALLEL_GENERATE_UNIQUES)
	{
		B3_PROFILE("generateUniques_parallel()");
		
		m_generateUniquesProgram.generateUniques(commandQueue, m_valueIndexPairs, gridData, numFluidParticles);
		clFinish(commandQueue);
	}
	else
	{
		B3_PROFILE("generateUniques_serial()");
		generateUniques_serial(commandQueue, numFluidParticles, gridData);
		
		int numActiveCells = gridData->getNumActiveCells();
		gridData->m_activeCells.resize(numActiveCells);
		gridData->m_cellContents.resize(numActiveCells);
		clFinish(commandQueue);
	}
	
	{
		B3_PROFILE("copy valueIndexPairs to host");
		m_valueIndexPairs.copyToHost(m_valueIndexPairsHost, true);
	}
	
	clFinish(commandQueue);
}
void b3FluidSortingGridOpenCLProgram::rearrangeParticlesOnHost(b3FluidSph* fluid)
{
	B3_PROFILE("rearrange host");
		
	b3FluidParticles& particles = fluid->internalGetParticles();
	rearrangeToMatchSortedValues2(m_valueIndexPairsHost, m_tempBufferVector, particles.m_pos);
	rearrangeToMatchSortedValues2(m_valueIndexPairsHost, m_tempBufferVector, particles.m_vel);
	rearrangeToMatchSortedValues2(m_valueIndexPairsHost, m_tempBufferVector, particles.m_vel_eval);
	rearrangeToMatchSortedValues2(m_valueIndexPairsHost, m_tempBufferVector, particles.m_accumulatedForce);
	rearrangeToMatchSortedValues2(m_valueIndexPairsHost, m_tempBufferVoid, particles.m_userPointer);
}

void b3FluidSortingGridOpenCLProgram::generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, 
															  b3Scalar cellSize, cl_mem fluidPositionsBuffer)
{
	b3BufferInfoCL bufferInfo[] = 
	{
		b3BufferInfoCL( fluidPositionsBuffer ),
		b3BufferInfoCL( m_valueIndexPairs.getBufferCL() )
	};
	
	b3LauncherCL launcher(commandQueue, kernel_generateValueIndexPairs);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(cellSize);
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
}
void b3FluidSortingGridOpenCLProgram::rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, cl_mem fluidBuffer)
{
	b3BufferInfoCL bufferInfo[] = 
	{
		b3BufferInfoCL( m_valueIndexPairs.getBufferCL() ),
		b3BufferInfoCL( fluidBuffer ),
		b3BufferInfoCL( m_tempBufferCL.getBufferCL() )
	};
	
	b3LauncherCL launcher(commandQueue, kernel_rearrangeParticleArrays);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
}

void b3FluidSortingGridOpenCLProgram::generateUniques_serial(cl_command_queue commandQueue, int numFluidParticles, b3FluidSortingGridOpenCL* gridData)
{
	b3BufferInfoCL bufferInfo[] = 
	{
		b3BufferInfoCL( m_valueIndexPairs.getBufferCL() ),
		b3BufferInfoCL( gridData->m_activeCells.getBufferCL() ),
		b3BufferInfoCL( gridData->m_cellContents.getBufferCL() ),
		b3BufferInfoCL( gridData->m_numActiveCells.getBufferCL() )
	};
	
	b3LauncherCL launcher(commandQueue, kernel_generateUniques);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(1);
}
