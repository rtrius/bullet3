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

#include "Bullet3Fluids/Sph/b3FluidSphParameters.h"
#include "Bullet3Fluids/Sph/b3FluidParticles.h"

void b3FluidSphOpenCL::writeToOpenCL(cl_command_queue queue, const b3FluidSphParametersLocal& FL, b3FluidParticles& particles)
{
	if(m_initialized) return;

	m_localParameters.resize(1);
	m_localParameters.copyFromHostPointer(&FL, 1, 0, false);
	
	int numParticles = particles.size();
	m_pos.resize(numParticles);
	m_vel.resize(numParticles);
	m_vel_eval.resize(numParticles);
	m_accumulatedForce.resize(numParticles);
	m_sph_force.resize(numParticles);
	m_density.resize(numParticles);
	m_cellIndex.resize(numParticles);
	
	m_pos.copyFromHost(particles.m_pos, false);
	m_vel.copyFromHost(particles.m_vel, false);
	m_vel_eval.copyFromHost(particles.m_vel_eval, false);
	m_accumulatedForce.copyFromHost(particles.m_accumulatedForce, false);
	
	clFinish(queue);
	m_initialized = true;
}
	
void b3FluidSphOpenCL::readFromOpenCL(cl_command_queue queue, b3FluidParticles& particles, b3AlignedObjectArray<b3Vector3>& sphForce)
{
	m_pos.copyToHost(particles.m_pos, false);
	/*m_vel.copyToHost(particles.m_vel, false);
	m_vel_eval.copyToHost(particles.m_vel_eval, false);
	m_accumulatedForce.copyToHost(particles.m_accumulatedForce, false);
	m_sph_force.copyToHost(sphForce, false);*/
	clFinish(queue);
}
	