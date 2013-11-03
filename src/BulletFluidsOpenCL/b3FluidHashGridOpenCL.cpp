/*
Bullet-FLUIDS 
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
#include "b3FluidHashGridOpenCL.h"

#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"

#include "BulletFluids/Sph/b3FluidSortingGrid.h"

#include "b3FluidSphOpenCL.h"
#include "fluidSphCL2.h"

// /////////////////////////////////////////////////////////////////////////////
//class b3FluidHashGridOpenCL
// /////////////////////////////////////////////////////////////////////////////
b3FluidHashGridOpenCL::b3FluidHashGridOpenCL(cl_context context, cl_command_queue queue)
: m_cellContents(context, queue)
{
	m_cellContents.resize(B3_FLUID_HASH_GRID_NUM_CELLS);
}


	
// /////////////////////////////////////////////////////////////////////////////
//class b3FluidHashGridOpenCLProgram
// /////////////////////////////////////////////////////////////////////////////
b3FluidHashGridOpenCLProgram::b3FluidHashGridOpenCLProgram(cl_context context, cl_device_id device, cl_command_queue queue)
: m_tempBufferCL(context, queue), m_radixSorter(context, device, queue), m_valueIndexPairs(context, queue)
{
	const char CL_HASH_GRID_PROGRAM_PATH[] = "src/BulletFluidsOpenCL/fluidSph2.cl";
	
	const char* kernelSource = fluidSphCL2;	//fluidSphCL2.h
	cl_int error;
	char* additionalMacros = 0;
	m_hashGridProgram = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, 
																	additionalMacros, CL_HASH_GRID_PROGRAM_PATH);
	b3Assert(m_hashGridProgram);

	m_generateValueIndexPairsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "generateValueIndexPairs", &error, m_hashGridProgram, additionalMacros );
	b3Assert(m_generateValueIndexPairsKernel);
	m_rearrangeParticleArraysKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "rearrangeParticleArrays", &error, m_hashGridProgram, additionalMacros );
	b3Assert(m_rearrangeParticleArraysKernel);
	
	m_resetGridCellsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "resetGridCells", &error, m_hashGridProgram, additionalMacros );
	b3Assert(m_resetGridCellsKernel);
	m_detectIndexRangesKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "detectIndexRanges", &error, m_hashGridProgram, additionalMacros );
	b3Assert(m_detectIndexRangesKernel);
}
b3FluidHashGridOpenCLProgram::~b3FluidHashGridOpenCLProgram()
{
	clReleaseKernel(m_generateValueIndexPairsKernel);
	clReleaseKernel(m_rearrangeParticleArraysKernel);
	
	clReleaseKernel(m_resetGridCellsKernel);
	clReleaseKernel(m_detectIndexRangesKernel);
	
	clReleaseProgram(m_hashGridProgram);
}

void b3FluidHashGridOpenCLProgram::insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
															  b3FluidSph* fluid, b3FluidSphOpenCL* fluidData, b3FluidHashGridOpenCL* gridData)
{
	B3_PROFILE("b3FluidHashGridOpenCLProgram::insertParticlesIntoGrid()");
	
	int numFluidParticles = fluid->numParticles();
	m_tempBufferCL.resize(numFluidParticles);
	m_valueIndexPairs.resize(numFluidParticles);
	
	b3Scalar gridCellSize = fluid->getGrid().getCellSize();
	
	//
	{
		B3_PROFILE("generateValueIndexPairs()");
		generateValueIndexPairs( commandQueue, numFluidParticles, gridCellSize, fluidData->m_pos.getBufferCL() );
		
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
	
	//Reset grid cells
	{
		b3BufferInfoCL bufferInfo[] = 
		{
			b3BufferInfoCL( gridData->m_cellContents.getBufferCL() )
		};
		
		b3LauncherCL launcher(commandQueue, m_resetGridCellsKernel);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
		launcher.setConst(B3_FLUID_HASH_GRID_NUM_CELLS);
		
		launcher.launch1D(B3_FLUID_HASH_GRID_NUM_CELLS);
	}
	
	//Detect index ranges for each cell
	{
		b3BufferInfoCL bufferInfo[] = 
		{
			b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
			b3BufferInfoCL( m_valueIndexPairs.getBufferCL() ),
			b3BufferInfoCL( gridData->m_cellContents.getBufferCL() )
		};
		
		b3LauncherCL launcher(commandQueue, m_detectIndexRangesKernel);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
		launcher.setConst(gridCellSize);
		launcher.setConst(numFluidParticles);
		
		launcher.launch1D(numFluidParticles);
	}
	
	
}
void b3FluidHashGridOpenCLProgram::generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, 
															  b3Scalar cellSize, cl_mem fluidPositionsBuffer)
{
	b3BufferInfoCL bufferInfo[] = 
	{
		b3BufferInfoCL( fluidPositionsBuffer ),
		b3BufferInfoCL( m_valueIndexPairs.getBufferCL() )
	};
	
	b3LauncherCL launcher(commandQueue, m_generateValueIndexPairsKernel);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(cellSize);
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
}
void b3FluidHashGridOpenCLProgram::rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, cl_mem fluidBuffer)
{
	b3BufferInfoCL bufferInfo[] = 
	{
		b3BufferInfoCL( m_valueIndexPairs.getBufferCL() ),
		b3BufferInfoCL( fluidBuffer ),
		b3BufferInfoCL( m_tempBufferCL.getBufferCL() )
	};
	
	b3LauncherCL launcher(commandQueue, m_rearrangeParticleArraysKernel);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
}
