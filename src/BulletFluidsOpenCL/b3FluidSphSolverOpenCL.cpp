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
#include "b3FluidSphSolverOpenCL.h"

#include "Bullet3Common/b3Quickprof.h"		//B3_PROFILE(name) macro

#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"

#include "BulletFluids/Sph/b3FluidSphParameters.h"

#include "fluidSphCL.h"

b3FluidSphSolverOpenCL::b3FluidSphSolverOpenCL(cl_context context, cl_command_queue queue, cl_device_id device)
: m_globalFluidParams(context, queue), m_sortingGridProgram(context, queue, device)
{
	m_context = context;
	m_commandQueue = queue;

	//
	const char CL_PROGRAM_PATH[] = "src/BulletFluidsOpenCL/fluidSph.cl";
	
	const char* kernelSource = fluidSphCL;	//fluidSphCL.h
	cl_int error;
	char* additionalMacros = 0;
	m_fluidsProgram = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, 
																	additionalMacros, CL_PROGRAM_PATH);
	b3Assert(m_fluidsProgram);
	
	m_kernel_findNeighborCellsPerCell = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "findNeighborCellsPerCell", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_kernel_findNeighborCellsPerCell);
	m_kernel_findGridCellIndexPerParticle = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "findGridCellIndexPerParticle", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_kernel_findGridCellIndexPerParticle);
	m_kernel_sphComputePressure = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "sphComputePressure", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_kernel_sphComputePressure);
	m_kernel_sphComputeForce = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "sphComputeForce", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_kernel_sphComputeForce);
	m_applyForcesKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "applyForces", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_applyForcesKernel);
	m_collideAabbImpulseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "collideAabbImpulse", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_collideAabbImpulseKernel);
	m_integratePositionKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "integratePositions", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_integratePositionKernel);
}

b3FluidSphSolverOpenCL::~b3FluidSphSolverOpenCL()
{
	clReleaseKernel(m_kernel_findNeighborCellsPerCell);
	clReleaseKernel(m_kernel_findGridCellIndexPerParticle);
	clReleaseKernel(m_kernel_sphComputePressure);
	clReleaseKernel(m_kernel_sphComputeForce);
	clReleaseKernel(m_applyForcesKernel);
	clReleaseKernel(m_collideAabbImpulseKernel);
	clReleaseKernel(m_integratePositionKernel);
	clReleaseProgram(m_fluidsProgram);
}

///BULLET_2_TO_3_PLACEHOLDER
inline void resolveAabbCollision_impulse(const b3FluidSphParametersGlobal& FG, const b3FluidSphParametersLocal& FL, const b3Vector3& velocity, 
										const b3Vector3& normal, b3Scalar distance, b3Vector3& out_impulse)
{
	if( distance < b3Scalar(0.0) )	//Negative distance indicates penetration
	{
		b3Scalar penetratingMagnitude = velocity.dot(-normal);
		if( penetratingMagnitude < b3Scalar(0.0) ) penetratingMagnitude = b3Scalar(0.0);
		
		b3Vector3 penetratingVelocity = -normal * penetratingMagnitude;
		b3Vector3 tangentialVelocity = velocity - penetratingVelocity;
		
		penetratingVelocity *= b3Scalar(1.0) + FL.m_boundaryRestitution;
		
		b3Scalar positionError = (-distance) * (FG.m_simulationScale/FG.m_timeStep) * FL.m_boundaryErp;
		
		out_impulse += -( penetratingVelocity + (-normal*positionError) + tangentialVelocity * FL.m_boundaryFriction );
	}
}
void accumulateBoundaryImpulse(const b3FluidSphParametersGlobal& FG, b3Scalar simScaleParticleRadius,
								const b3FluidSphParametersLocal& FL, b3FluidParticles& particles, int particleIndex,
								b3Vector3& out_accumulatedImpulse)
{
	int i = particleIndex;
	
	const b3Scalar radius = simScaleParticleRadius;
	const b3Scalar simScale = FG.m_simulationScale;
	
	const b3Vector3& boundaryMin = FL.m_aabbBoundaryMin;
	const b3Vector3& boundaryMax = FL.m_aabbBoundaryMax;
	
	const b3Vector3& pos = particles.m_pos[i];
	const b3Vector3& vel = particles.m_vel[i];
	
	b3Vector3& impulse = out_accumulatedImpulse;
	resolveAabbCollision_impulse( FG, FL, vel, b3Vector3( 1.0, 0.0, 0.0), ( pos.getX() - boundaryMin.getX() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, b3Vector3(-1.0, 0.0, 0.0), ( boundaryMax.getX() - pos.getX() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, b3Vector3(0.0,  1.0, 0.0), ( pos.getY() - boundaryMin.getY() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, b3Vector3(0.0, -1.0, 0.0), ( boundaryMax.getY() - pos.getY() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, b3Vector3(0.0, 0.0,  1.0), ( pos.getZ() - boundaryMin.getZ() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, b3Vector3(0.0, 0.0, -1.0), ( boundaryMax.getZ() - pos.getZ() )*simScale - radius, impulse );
}
void applyAabbImpulsesSingleFluid(const b3FluidSphParametersGlobal& FG, b3FluidSph* fluid)
{
	B3_PROFILE("applyAabbImpulsesSingleFluid()");
	
	const b3FluidSphParametersLocal& FL = fluid->getLocalParameters();
	b3FluidParticles& particles = fluid->internalGetParticles();
	
	const b3Scalar simScaleParticleRadius = FL.m_particleRadius * FG.m_simulationScale;
	
	for(int i = 0; i < particles.size(); ++i)
	{
		b3Vector3 aabbImpulse(0, 0, 0);
		accumulateBoundaryImpulse(FG, simScaleParticleRadius, FL, particles, i, aabbImpulse);
		
		b3Vector3& vel = particles.m_vel[i];
		b3Vector3& vel_eval = particles.m_vel_eval[i];
		
		//Leapfrog integration
		b3Vector3 velNext = vel + aabbImpulse;
		vel_eval = (vel + velNext) * b3Scalar(0.5);
		vel = velNext;
	}
}
///BULLET_2_TO_3_PLACEHOLDER

void b3FluidSphSolverOpenCL::updateGridAndCalculateSphForces(const b3FluidSphParametersGlobal& FG, b3FluidSph** fluids, int numFluids)
{	
	B3_PROFILE("b3FluidSphSolverOpenCL::updateGridAndCalculateSphForces()");
	
#ifdef B3_USE_DOUBLE_PRECISION
	b3Assert(0 && "B3_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
	return;
#endif	

	b3AlignedObjectArray<b3FluidSph*> validFluids;
	for(int i = 0; i < numFluids; ++i) 
	{
		if( fluids[i]->numParticles() ) 
		{
			validFluids.push_back( fluids[i] );
		}
		else
		{
			//Update AABB
			b3FluidSortingGrid& grid = fluids[i]->internalGetGrid();
			b3Vector3& pointAabbMin = grid.internalGetPointAabbMin();
			b3Vector3& pointAabbMax = grid.internalGetPointAabbMax();
			
			pointAabbMin.setValue(0,0,0);
			pointAabbMax.setValue(0,0,0);
		}
	}

	int numValidFluids = validFluids.size();
	
//B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT is not supported when using OpenCL grid update.
#ifdef B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT
	const bool UPDATE_GRID_ON_GPU = false;
	b3Assert(0);	//	current implementation requires GPU grid update
#else
	const bool UPDATE_GRID_ON_GPU = true;
#endif
	
	if(!UPDATE_GRID_ON_GPU)
		for(int i = 0; i < numValidFluids; ++i) validFluids[i]->insertParticlesIntoGrid();
	
	//Write data from CPU to OpenCL
	m_globalFluidParams.resize(1);
	m_globalFluidParams.copyFromHostPointer(&FG, 1, 0, false);
	clFinish(m_commandQueue);

	//resize m_gridData to numValidFluids
	{
		while(m_gridData.size() > numValidFluids)
		{
			b3FluidSortingGridOpenCL* lastElement = m_gridData[m_gridData.size() - 1];
			lastElement->~b3FluidSortingGridOpenCL();
			b3AlignedFree(lastElement);
			
			m_gridData.pop_back();
		}
		while(m_gridData.size() < numValidFluids)
		{
			void* ptr = b3AlignedAlloc( sizeof(b3FluidSortingGridOpenCL), 16 );
			b3FluidSortingGridOpenCL* newElement = new(ptr) b3FluidSortingGridOpenCL(m_context, m_commandQueue);
			
			m_gridData.push_back(newElement);
		}
	}
	//resize m_fluidData to numValidFluids
	{
		while(m_fluidData.size() > numValidFluids)
		{
			b3FluidSphOpenCL* lastElement = m_fluidData[m_fluidData.size() - 1];
			lastElement->~b3FluidSphOpenCL();
			b3AlignedFree(lastElement);
			
			m_fluidData.pop_back();
		}
		while(m_fluidData.size() < numValidFluids)
		{
			void* ptr = b3AlignedAlloc( sizeof(b3FluidSphOpenCL), 16 );
			b3FluidSphOpenCL* newElement = new(ptr) b3FluidSphOpenCL(m_context, m_commandQueue);
			
			m_fluidData.push_back(newElement);
		}
	}
	
	{
		B3_PROFILE("writeToOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			const b3FluidSphParametersLocal& FL = validFluids[i]->getLocalParameters();
			
			if(!UPDATE_GRID_ON_GPU) m_gridData[i]->writeToOpenCL( m_commandQueue, validFluids[i]->internalGetGrid() );
			m_fluidData[i]->writeToOpenCL( m_commandQueue, FL, validFluids[i]->internalGetParticles() );
		}
	}
	
	//
	{
		B3_PROFILE("calculate sph force");
	
		for(int i = 0; i < numValidFluids; ++i)
		{
			b3FluidSph* fluid = validFluids[i];
			b3FluidSortingGridOpenCL* gridData = m_gridData[i];
			b3FluidSphOpenCL* fluidData = m_fluidData[i];
			
			if(UPDATE_GRID_ON_GPU)
				m_sortingGridProgram.insertParticlesIntoGrid(m_context, m_commandQueue, fluid, m_fluidData[i], m_gridData[i]);
			
			int numActiveCells = gridData->getNumActiveCells();
			int numFluidParticles = fluid->numParticles();
			
			findNeighborCells( numActiveCells, numFluidParticles, gridData, fluidData);
			sphComputePressure( numFluidParticles, gridData, fluidData, fluid->getGrid().getCellSize() );
			sphComputeForce( numFluidParticles, gridData, fluidData, fluid->getGrid().getCellSize() );
			
			clFlush(m_commandQueue);
			
			//This branch is executed on CPU, while findNeighborCells() and sphComputePressure/Force() is simultaneously executed on GPU
			if(0 && UPDATE_GRID_ON_GPU)
			{
				B3_PROFILE("simultaneous rearrange/AABB update");
			
				m_sortingGridProgram.rearrangeParticlesOnHost(fluid);
				
				{
					B3_PROFILE("Update AABB - CPU");
					b3Vector3 aabbMin(B3_LARGE_FLOAT, B3_LARGE_FLOAT, B3_LARGE_FLOAT);
					b3Vector3 aabbMax(-B3_LARGE_FLOAT, -B3_LARGE_FLOAT, -B3_LARGE_FLOAT);
					for(int i = 0; i < fluid->numParticles(); ++i) 
					{
						const b3Vector3& position = fluid->getPosition(i);
							
						aabbMin.setMin(position);
						aabbMax.setMax(position);
					}
					
					b3FluidSortingGrid& grid = fluid->internalGetGrid();
					b3Vector3& pointAabbMin = grid.internalGetPointAabbMin();
					b3Vector3& pointAabbMax = grid.internalGetPointAabbMax();
					
					pointAabbMin = aabbMin;
					pointAabbMax = aabbMax;
				}
			}
			
			clFinish(m_commandQueue);
		}
	}
	
	const bool GPU_INTEGRATE = true;
	if(GPU_INTEGRATE)
	{
		B3_PROFILE("apply boundary impulses, integrate");
		
		for(int i = 0; i < numValidFluids; ++i)
		{
			b3FluidSph* fluid = validFluids[i];
			b3FluidSphOpenCL* fluidData = m_fluidData[i];
			int numFluidParticles = fluid->numParticles();
		
			{
				b3BufferInfoCL bufferInfo[] = 
				{ 
					b3BufferInfoCL( m_globalFluidParams.getBufferCL() ), 
					b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_accumulatedForce.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_sph_force.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_vel.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_vel_eval.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_applyForcesKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidParticles);
				
				launcher.launch1D(numFluidParticles);
			}
			
			{
				b3BufferInfoCL bufferInfo[] = 
				{ 
					b3BufferInfoCL( m_globalFluidParams.getBufferCL() ), 
					b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_vel.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_vel_eval.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_collideAabbImpulseKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidParticles);
				
				launcher.launch1D(numFluidParticles);
			}
			
			{
				b3BufferInfoCL bufferInfo[] = 
				{ 
					b3BufferInfoCL( m_globalFluidParams.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_vel.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_integratePositionKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidParticles);
				
				launcher.launch1D(numFluidParticles);
			}
		}
		
		clFinish(m_commandQueue);
	}
	
	//Read data from OpenCL to CPU
	{
		B3_PROFILE("readFromOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			if(UPDATE_GRID_ON_GPU) m_gridData[i]->readFromOpenCL( m_commandQueue, validFluids[i]->internalGetGrid() );
		
			if( m_tempSphForce.size() < validFluids[i]->numParticles() ) m_tempSphForce.resize( validFluids[i]->numParticles() );
			m_fluidData[i]->readFromOpenCL( m_commandQueue, validFluids[i]->internalGetParticles(), m_tempSphForce);
			
			if(!GPU_INTEGRATE)
			{
				applySphForce(FG, validFluids[i], m_tempSphForce);
				
				b3FluidSphSolver::applyForcesSingleFluid(FG, validFluids[i]);
				applyAabbImpulsesSingleFluid(FG, validFluids[i]);
				
				b3FluidSphSolver::integratePositionsSingleFluid( FG, validFluids[i]->internalGetParticles() );
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
		
		b3LauncherCL launcher(m_commandQueue, m_kernel_findNeighborCellsPerCell);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
		
		launcher.launch1D(numActiveGridCells);
	}
	
	//For each particle, locate the grid cell that they are contained in so that 
	//they can use the results from m_kernel_findNeighborCellsPerCell, executed above
	{
		b3BufferInfoCL bufferInfo[] = 
		{
			b3BufferInfoCL( gridData->m_numActiveCells.getBufferCL() ),
			b3BufferInfoCL( gridData->m_cellContents.getBufferCL() ),
			b3BufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
		};
		
		b3LauncherCL launcher(m_commandQueue, m_kernel_findGridCellIndexPerParticle);
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
		b3BufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_density.getBufferCL() ),
		b3BufferInfoCL( gridData->m_foundCells.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
	};
	
	b3LauncherCL launcher(m_commandQueue, m_kernel_sphComputePressure);
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
		b3BufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_vel_eval.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_sph_force.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_density.getBufferCL() ),
		b3BufferInfoCL( gridData->m_foundCells.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
	};
	
	b3LauncherCL launcher(m_commandQueue, m_kernel_sphComputeForce);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	clFinish(m_commandQueue);
}
