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
#include "b3FluidSphOpenCL.h"

#include "Bullet3Fluids/Sph/b3FluidSph.h"

b3FluidSphOpenCL::b3FluidSphOpenCL(cl_context context, cl_command_queue queue) :
	m_parameters(context, queue),
	m_position(context, queue),
	m_velocity(context, queue),
	m_velocityEval(context, queue),
	m_accumulatedForce(context, queue),
	m_sph_force(context, queue),
	m_density(context, queue),
	m_cellIndex(context, queue)
{
	m_parameters.resize(1);
}
	
	
void b3FluidSphOpenCL::writeToOpenCL(cl_command_queue queue, b3FluidSph* fluid)
{
	const b3FluidSphSyncronizationFlags& syncFlags = fluid->getGpuSyncFlags();
	const b3FluidParticles& particles = fluid->getParticles();
	const b3FluidSphParameters& FP = fluid->getParameters();
	m_parameters.copyFromHostPointer(&FP, 1, 0, false);
	
	if( fluid->needsWriteStateToGpu() )
	{
		fluid->shouldWriteStateToGpu(false);
		
		int numParticles = particles.size();
		resize(numParticles);
	
		m_position.copyFromHost(particles.m_position, false);
		m_velocity.copyFromHost(particles.m_velocity, false);
		m_velocityEval.copyFromHost(particles.m_velocityEval, false);
		m_accumulatedForce.copyFromHost(particles.m_accumulatedForce, false);
	}
	else
	{
		if(syncFlags.m_writeForces) m_accumulatedForce.copyFromHost(particles.m_accumulatedForce, false);
	}
	
	clFinish(queue);
}
	
void b3FluidSphOpenCL::readFromOpenCL(cl_command_queue queue, b3FluidSph* fluid)
{
	const b3FluidSphSyncronizationFlags& syncFlags = fluid->getGpuSyncFlags();
	
	b3FluidParticles& particles = fluid->internalGetParticles();

	if(syncFlags.m_syncPosition) m_position.copyToHost(particles.m_position, false);
	if(syncFlags.m_syncVelocity) m_velocity.copyToHost(particles.m_velocity, false);
	if(syncFlags.m_syncVelocityEval) m_velocityEval.copyToHost(particles.m_velocityEval, false);
	
	//m_sph_force.copyToHost(sphForce, false);
	
	clFinish(queue);
}

void b3FluidSphOpenCL::resize(int size)
{
	const bool copyOldContents = true;
	m_position.resize(size, copyOldContents);
	m_velocity.resize(size, copyOldContents);
	m_velocityEval.resize(size, copyOldContents);
	m_accumulatedForce.resize(size, copyOldContents);
	
	m_sph_force.resize(size, copyOldContents);
	m_density.resize(size, copyOldContents);
	m_cellIndex.resize(size, copyOldContents);
	
}
	