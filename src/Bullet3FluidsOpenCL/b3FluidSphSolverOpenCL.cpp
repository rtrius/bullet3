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
#include "b3FluidSphSolverOpenCL.h"

#include "Bullet3Common/b3Logging.h"		//B3_PROFILE(name) macro

#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"

#include "Bullet3Fluids/Sph/b3FluidSphParameters.h"

#include "fluidSphCL.h"

b3FluidSphSolverOpenCL::b3FluidSphSolverOpenCL(cl_context context, cl_device_id device, cl_command_queue queue)
: m_sortingGridProgram(context, device, queue), m_hashGridProgram(context, device, queue), m_fluidRigidInteractor(context, device, queue)
{
	m_context = context;
	m_commandQueue = queue;

	//
	const char CL_PROGRAM_PATH[] = "src/Bullet3FluidsOpenCL/fluidSph.cl";
	
	const char* kernelSource = fluidSphCL;	//fluidSphCL.h
	cl_int error;
	char* additionalMacros = 0;
	m_fluidsProgram = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, 
																	additionalMacros, CL_PROGRAM_PATH);
	b3Assert(m_fluidsProgram);
	
	m_findNeighborCellsPerCellKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "findNeighborCellsPerCell", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_findNeighborCellsPerCellKernel);
	m_findGridCellIndexPerParticleKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "findGridCellIndexPerParticle", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_findGridCellIndexPerParticleKernel);
	m_sphComputePressureKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "sphComputePressure", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_sphComputePressureKernel);
	m_sphComputeForceKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "sphComputeForce", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_sphComputeForceKernel);
	m_applyForcesKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "applyForces", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_applyForcesKernel);
	m_collideAabbImpulseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "collideAabbImpulse", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_collideAabbImpulseKernel);
	m_integratePositionKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "integratePositions", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_integratePositionKernel);
	
	m_sphComputePressureModuloKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "sphComputePressureModulo", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_sphComputePressureModuloKernel);
	m_sphComputeForceModuloKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "sphComputeForceModulo", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_sphComputeForceModuloKernel);
}

b3FluidSphSolverOpenCL::~b3FluidSphSolverOpenCL()
{
	clReleaseKernel(m_sphComputePressureModuloKernel);
	clReleaseKernel(m_sphComputeForceModuloKernel);
	
	clReleaseKernel(m_findNeighborCellsPerCellKernel);
	clReleaseKernel(m_findGridCellIndexPerParticleKernel);
	clReleaseKernel(m_sphComputePressureKernel);
	clReleaseKernel(m_sphComputeForceKernel);
	clReleaseKernel(m_applyForcesKernel);
	clReleaseKernel(m_collideAabbImpulseKernel);
	clReleaseKernel(m_integratePositionKernel);
	
	clReleaseProgram(m_fluidsProgram);
}

void b3FluidSphSolverOpenCL::stepSimulation(b3FluidSph* fluid, RigidBodyGpuData& rbData)
{	
	B3_PROFILE("b3FluidSphSolverOpenCL::stepSimulation()");
	
#ifdef B3_USE_DOUBLE_PRECISION
	b3Assert(0 && "B3_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
	return;
#endif	

	if( !fluid->numParticles() ) 
	{
		//Update AABB
		b3FluidSortingGrid& grid = fluid->internalGetGrid();
		b3Vector3& pointAabbMin = grid.internalGetPointAabbMin();
		b3Vector3& pointAabbMax = grid.internalGetPointAabbMax();
		
		pointAabbMin.setValue(0,0,0);
		pointAabbMax.setValue(0,0,0);
		
		return;
	}
	
//B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT is not supported when using OpenCL grid update.
#ifdef B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT
	const bool UPDATE_GRID_ON_GPU = false;
	b3Assert(0);	//	current implementation requires GPU grid update
#else
	const bool UPDATE_GRID_ON_GPU = true;
#endif
	
	if(!UPDATE_GRID_ON_GPU) fluid->insertParticlesIntoGrid();
	
	
	//Allocate solver data b3FluidSphOpenCL and b3FluidSortingGridOpenCL
	b3FluidSphOpenCL* fluidDataCL = static_cast<b3FluidSphOpenCL*>( fluid->getSolverDataGpu() );
	b3FluidSortingGridOpenCL* gridDataCL = static_cast<b3FluidSortingGridOpenCL*>( fluid->getGridDataGpu() );
	b3FluidHashGridOpenCL* hashGridDataCL = static_cast<b3FluidHashGridOpenCL*>( fluid->getGridDataGpu() );
	{
		if(fluidDataCL && fluidDataCL->getType() != FSDT_b3FluidSphOpenCL)
		{
			b3AlignedFree(fluidDataCL);
			fluidDataCL = 0;
		}
		if(!fluidDataCL)
		{
			void* ptr = b3AlignedAlloc( sizeof(b3FluidSphOpenCL), 16 );
			fluidDataCL = new(ptr) b3FluidSphOpenCL(m_context, m_commandQueue);
			fluid->setSolverDataGpu(fluidDataCL);
		}
		
		if(!USE_HASH_GRID)
		{
			if(gridDataCL && gridDataCL->getType() != FSDT_b3FluidSortingGridOpenCL)
			{
				b3AlignedFree(gridDataCL);
				gridDataCL = 0;
			}
			if(!gridDataCL)
			{
				void* ptr = b3AlignedAlloc( sizeof(b3FluidSortingGridOpenCL), 16 );
				gridDataCL = new(ptr) b3FluidSortingGridOpenCL(m_context, m_commandQueue);
				fluid->setGridDataGpu(gridDataCL);
			}
		}
		else
		{
			if(hashGridDataCL && hashGridDataCL->getType() != FSDT_b3FluidHashGridOpenCL)
			{
				b3AlignedFree(hashGridDataCL);
				gridDataCL = 0;
			}
			if(!hashGridDataCL)
			{
				void* ptr = b3AlignedAlloc( sizeof(b3FluidHashGridOpenCL), 16 );
				hashGridDataCL = new(ptr) b3FluidHashGridOpenCL(m_context, m_commandQueue);
				fluid->setGridDataGpu(hashGridDataCL);
			}
		}
	}
	
	//Write data from CPU to OpenCL
	{
		B3_PROFILE("writeToOpenCL");
		
		const b3FluidSphParameters& FP = fluid->getParameters();
			
		if(!USE_HASH_GRID && !UPDATE_GRID_ON_GPU) gridDataCL->writeToOpenCL( m_commandQueue, fluid->internalGetGrid() );
		fluidDataCL->writeToOpenCL(m_commandQueue, fluid);
	}
	
	//
	{
		B3_PROFILE("calculate sph force");
		
		int numFluidParticles = fluid->numParticles();
			
		if(!USE_HASH_GRID)
		{
			if(UPDATE_GRID_ON_GPU)
				m_sortingGridProgram.insertParticlesIntoGrid(m_context, m_commandQueue, fluid, fluidDataCL, gridDataCL);
			
			int numActiveCells = gridDataCL->getNumActiveCells();
			
			findNeighborCells( numActiveCells, numFluidParticles, gridDataCL, fluidDataCL);
			sphComputePressure( numFluidParticles, gridDataCL, fluidDataCL, fluid->getGrid().getCellSize() );
			sphComputeForce( numFluidParticles, gridDataCL, fluidDataCL, fluid->getGrid().getCellSize() );
		}
		else
		{
			m_hashGridProgram.insertParticlesIntoGrid(m_context, m_commandQueue, fluid, fluidDataCL, hashGridDataCL);
		
			sphComputePressureModulo( numFluidParticles, hashGridDataCL, fluidDataCL, fluid->getGrid().getCellSize() );
			sphComputeForceModulo( numFluidParticles, hashGridDataCL, fluidDataCL, fluid->getGrid().getCellSize() );
		}
		
		//The previous vel(velocity at t-1/2) is needed to update vel_eval for leapfrog integration
		//Since vel_eval(velocity at t) is used only for SPH force computation,
		//it is possible to store the previous velocity in m_velocityEval
		fluidDataCL->m_velocityEval.copyFromOpenCLArray(fluidDataCL->m_velocity);
		
		clFinish(m_commandQueue);
	}
	
	const bool GPU_INTEGRATE = true;
	if(GPU_INTEGRATE)
	{
		B3_PROFILE("apply boundary impulses, integrate");
		
		int numFluidParticles = fluid->numParticles();
	
		{
			b3BufferInfoCL bufferInfo[] = 
			{ 
				b3BufferInfoCL( fluidDataCL->m_parameters.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_accumulatedForce.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_sph_force.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocity.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocityEval.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_applyForcesKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
		
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( fluidDataCL->m_parameters.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_position.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocity.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocityEval.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_collideAabbImpulseKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
		
		if(!USE_HASH_GRID) m_fluidRigidInteractor.interact(fluidDataCL, gridDataCL, 0, rbData);
		else m_fluidRigidInteractor.interact(fluidDataCL, 0, hashGridDataCL, rbData);
		
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( fluidDataCL->m_parameters.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_position.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocity.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocityEval.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_integratePositionKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
		
		
		clFinish(m_commandQueue);
	}
	
	//Read data from OpenCL to CPU
	{
		B3_PROFILE("readFromOpenCL");
		
		{
			if(!USE_HASH_GRID && UPDATE_GRID_ON_GPU) gridDataCL->readFromOpenCL( m_commandQueue, fluid->internalGetGrid() );
		
			if( m_tempSphForce.size() < fluid->numParticles() ) m_tempSphForce.resize( fluid->numParticles() );
			fluidDataCL->readFromOpenCL(m_commandQueue, fluid);
			
			if(!GPU_INTEGRATE)
			{
				b3Assert(0);
			
				b3FluidParticles& particles = fluid->internalGetParticles();
				particles.m_velocityEval = particles.m_velocity;
			
				applySphForce(fluid, m_tempSphForce);
				
				b3FluidSphSolver::applyForcesSingleFluid(fluid);
				b3FluidSphSolver::applyAabbImpulsesSingleFluid(fluid);
				
				b3FluidSphSolver::integratePositionsSingleFluid( fluid->getParameters(), particles );
			}
		}
	}
}


void b3FluidSphSolverOpenCL::findNeighborCells(int numActiveGridCells, int numFluidParticles, 
												b3FluidSortingGridOpenCL* gridData, b3FluidSphOpenCL* fluidData)
{
	B3_PROFILE("findNeighborCells");
	
	//Perform 9 binary searches per cell, to locate the 27 neighbor cells
	{
		gridData->m_foundCells.resize(numActiveGridCells);
		
		b3BufferInfoCL bufferInfo[] = 
		{ 
			b3BufferInfoCL( gridData->m_numActiveCells.getBufferCL() ),
			b3BufferInfoCL( gridData->m_activeCells.getBufferCL() ),
			b3BufferInfoCL( gridData->m_cellContents.getBufferCL() ),
			b3BufferInfoCL( gridData->m_foundCells.getBufferCL() )
		};
		
		b3LauncherCL launcher(m_commandQueue, m_findNeighborCellsPerCellKernel);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
		
		launcher.launch1D(numActiveGridCells);
	}
	
	//For each particle, locate the grid cell that they are contained in so that 
	//they can use the results from m_findNeighborCellsPerCellKernel, executed above
	{
		b3BufferInfoCL bufferInfo[] = 
		{
			b3BufferInfoCL( gridData->m_numActiveCells.getBufferCL() ),
			b3BufferInfoCL( gridData->m_cellContents.getBufferCL() ),
			b3BufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
		};
		
		b3LauncherCL launcher(m_commandQueue, m_findGridCellIndexPerParticleKernel);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
		
		launcher.launch1D(numActiveGridCells);
	}
	
	clFinish(m_commandQueue);
}
void b3FluidSphSolverOpenCL::sphComputePressure(int numFluidParticles, b3FluidSortingGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize) 
{
	B3_PROFILE("sphComputePressure");
	
	b3BufferInfoCL bufferInfo[] = 
	{ 
		b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_density.getBufferCL() ),
		b3BufferInfoCL( gridData->m_foundCells.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
	};
	
	b3LauncherCL launcher(m_commandQueue, m_sphComputePressureKernel);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	clFinish(m_commandQueue);
}
void b3FluidSphSolverOpenCL::sphComputeForce(int numFluidParticles, b3FluidSortingGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize) 
{
	B3_PROFILE("sphComputeForce");
	
	b3BufferInfoCL bufferInfo[] = 
	{ 
		b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_velocityEval.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_sph_force.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_density.getBufferCL() ),
		b3BufferInfoCL( gridData->m_foundCells.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
	};
	
	b3LauncherCL launcher(m_commandQueue, m_sphComputeForceKernel);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	clFinish(m_commandQueue);
}


void b3FluidSphSolverOpenCL::sphComputePressureModulo(int numFluidParticles, b3FluidHashGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize) 
{
	B3_PROFILE("sphComputePressureModulo");
	
	b3BufferInfoCL bufferInfo[] = 
	{ 
		b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_density.getBufferCL() ),
		b3BufferInfoCL( gridData->m_cellContents.getBufferCL() )
	};
	
	b3LauncherCL launcher(m_commandQueue, m_sphComputePressureModuloKernel);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(cellSize);
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	clFinish(m_commandQueue);
}
void b3FluidSphSolverOpenCL::sphComputeForceModulo(int numFluidParticles, b3FluidHashGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize) 
{
	B3_PROFILE("sphComputeForceModulo");
	
	b3BufferInfoCL bufferInfo[] = 
	{ 
		b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_velocityEval.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_sph_force.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_density.getBufferCL() ),
		b3BufferInfoCL( gridData->m_cellContents.getBufferCL() )
	};
	
	b3LauncherCL launcher(m_commandQueue, m_sphComputeForceModuloKernel);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(cellSize);
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	clFinish(m_commandQueue);
}

