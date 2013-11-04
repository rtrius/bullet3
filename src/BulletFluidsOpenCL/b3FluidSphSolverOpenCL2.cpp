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
#include "b3FluidSphSolverOpenCL2.h"

#include "Bullet3Common/b3Quickprof.h"		//B3_PROFILE(name) macro

#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"

#include "BulletFluids/Sph/b3FluidSphParameters.h"

#include "fluidSphCL2.h"

b3FluidSphSolverOpenCL2::b3FluidSphSolverOpenCL2(cl_context context, cl_device_id device, cl_command_queue queue)
: m_globalFluidParams(context, queue), m_hashGridProgram(context, device, queue), m_fluidRigidInteractor(context, device, queue)
{
	m_context = context;
	m_commandQueue = queue;

	//
	const char CL_PROGRAM_PATH[] = "src/BulletFluidsOpenCL/fluidSph2.cl";
	const char* kernelSource = fluidSphCL2;	//fluidSphCL2.h
	cl_int error;
	char* additionalMacros = 0;
	m_fluidsProgram = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, 
																	additionalMacros, CL_PROGRAM_PATH);
	b3Assert(m_fluidsProgram);

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
}

b3FluidSphSolverOpenCL2::~b3FluidSphSolverOpenCL2()
{
	clReleaseKernel(m_sphComputePressureKernel);
	clReleaseKernel(m_sphComputeForceKernel);
	
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

void b3FluidSphSolverOpenCL2::stepSimulation(const b3FluidSphParametersGlobal& FG, b3FluidSph** fluids, int numFluids, RigidBodyGpuData& rbData)
{	
	B3_PROFILE("b3FluidSphSolverOpenCL2::stepSimulation()");
	
#ifdef B3_USE_DOUBLE_PRECISION
	b3Assert(0 && "B3_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
	return;
#endif	

	b3AlignedObjectArray<b3FluidSph*> validFluids;
	for(int i = 0; i < numFluids; ++i) 
	{
		fluids[i]->setFluidDataCL(0);
		fluids[i]->setGridDataCL(0);
		
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
	
	//if(!UPDATE_GRID_ON_GPU)
	//	for(int i = 0; i < numValidFluids; ++i) validFluids[i]->insertParticlesIntoGrid();
	
	//Write data from CPU to OpenCL
	m_globalFluidParams.resize(1);
	m_globalFluidParams.copyFromHostPointer(&FG, 1, 0, false);
	clFinish(m_commandQueue);
	
	//resize m_gridData to numValidFluids
	{
		while(m_gridData.size() > numValidFluids)
		{
			b3FluidHashGridOpenCL* lastElement = m_gridData[m_gridData.size() - 1];
			lastElement->~b3FluidHashGridOpenCL();
			b3AlignedFree(lastElement);
			
			m_gridData.pop_back();
		}
		while(m_gridData.size() < numValidFluids)
		{
			void* ptr = b3AlignedAlloc( sizeof(b3FluidHashGridOpenCL), 16 );
			b3FluidHashGridOpenCL* newElement = new(ptr) b3FluidHashGridOpenCL(m_context, m_commandQueue);
			
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
			
			//if(!UPDATE_GRID_ON_GPU) m_gridData[i]->writeToOpenCL( m_commandQueue, validFluids[i]->internalGetGrid() );
			m_fluidData[i]->writeToOpenCL( m_commandQueue, FL, validFluids[i]->internalGetParticles() );
			
			validFluids[i]->setFluidDataCL( static_cast<void*>(m_fluidData[i]) );
			validFluids[i]->setGridDataCL( static_cast<void*>(m_gridData[i]) );
		}
	}
	
	//
	{
		B3_PROFILE("calculate sph force");
	
		for(int i = 0; i < numValidFluids; ++i)
		{
			b3FluidSph* fluid = validFluids[i];
			b3FluidHashGridOpenCL* gridData = m_gridData[i];
			b3FluidSphOpenCL* fluidData = m_fluidData[i];
			
			//if(UPDATE_GRID_ON_GPU)
				m_hashGridProgram.insertParticlesIntoGrid(m_context, m_commandQueue, fluid, fluidData, gridData);
			
			int numFluidParticles = fluid->numParticles();
			sphComputePressure( numFluidParticles, gridData, fluidData, fluid->getGrid().getCellSize() );
			sphComputeForce( numFluidParticles, gridData, fluidData, fluid->getGrid().getCellSize() );
			
			//The previous vel(velocity at t-1/2) is needed to update vel_eval for leapfrog integration
			//Since vel_eval(velocity at t) is used only for SPH force computation,
			//it is possible to store the previous velocity in m_vel_eval
			fluidData->m_vel_eval.copyFromOpenCLArray(fluidData->m_vel);
			
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
			b3FluidHashGridOpenCL* gridData = m_gridData[i];
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
			
			m_fluidRigidInteractor.interact(m_globalFluidParams, fluidData, 0, gridData, rbData);
			
			{
				b3BufferInfoCL bufferInfo[] = 
				{
					b3BufferInfoCL( m_globalFluidParams.getBufferCL() ), 
					b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_vel.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_vel_eval.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_integratePositionKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidParticles);
				
				launcher.launch1D(numFluidParticles);
			}
		}
		
		clFinish(m_commandQueue);
	}
	
	//Send data from OpenCL to CPU
	{
		B3_PROFILE("readFromOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			//if(UPDATE_GRID_ON_GPU) m_gridData[i]->readFromOpenCL( m_commandQueue, validFluids[i]->internalGetGrid() );
		
			if( m_tempSphForce.size() < validFluids[i]->numParticles() ) m_tempSphForce.resize( validFluids[i]->numParticles() );
			m_fluidData[i]->readFromOpenCL( m_commandQueue, validFluids[i]->internalGetParticles(), m_tempSphForce);
			
			if(!GPU_INTEGRATE)
			{
				b3FluidParticles& particles = validFluids[i]->internalGetParticles();
				particles.m_vel_eval = particles.m_vel;
			
				applySphForce(FG, validFluids[i], m_tempSphForce);
				
				b3FluidSphSolver::applyForcesSingleFluid(FG, validFluids[i]);
				applyAabbImpulsesSingleFluid(FG, validFluids[i]);
				
				b3FluidSphSolver::integratePositionsSingleFluid( FG, validFluids[i]->getLocalParameters(), particles );
			}
		}
	}
}

void b3FluidSphSolverOpenCL2::sphComputePressure(int numFluidParticles, b3FluidHashGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize) 
{
	B3_PROFILE("sphComputePressure");
	
	b3BufferInfoCL bufferInfo[] = 
	{ 
		b3BufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
		b3BufferInfoCL( fluidData->m_density.getBufferCL() ),
		b3BufferInfoCL( gridData->m_cellContents.getBufferCL() )
	};
	
	b3LauncherCL launcher(m_commandQueue, m_sphComputePressureKernel);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(cellSize);
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	clFinish(m_commandQueue);
}
void b3FluidSphSolverOpenCL2::sphComputeForce(int numFluidParticles, b3FluidHashGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize) 
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
		b3BufferInfoCL( gridData->m_cellContents.getBufferCL() )
	};
	
	b3LauncherCL launcher(m_commandQueue, m_sphComputeForceKernel);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
	launcher.setConst(cellSize);
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	clFinish(m_commandQueue);
}
