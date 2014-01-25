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
//Portions of this file based on FLUIDS v.2 - SPH Fluid Simulator for CPU and GPU
//Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

#include "b3FluidSph.h"

#include "Bullet3Common/b3Logging.h"		//B3_PROFILE(name) macro
#include "Bullet3Common/b3Random.h"			//b3rand(), B3_RAND_MAX
#include "Bullet3Geometry/b3AabbUtil.h"		//b3TestPointAgainstAabb2()

#include "b3FluidSortingGrid.h"

b3FluidSph::b3FluidSph(b3FluidSphSolver* solver, int maxNumParticles)
{
	m_solver = solver;

	m_solverData = 0;
	
	m_solverDataGpu = 0;
	m_gridDataGpu = 0;
	
	m_copyParticleDataToGpu = false;	//Use incremental update by default
	
	setMaxParticles(maxNumParticles);
	m_grid.setCellSize(m_parameters.m_simulationScale, m_parameters.m_sphSmoothRadius);
}
b3FluidSph::~b3FluidSph()
{
}

void b3FluidSph::setMaxParticles(int maxNumParticles)
{
	m_copyParticleDataToGpu = true;
	
	if( maxNumParticles < m_particles.size() )m_particles.resize(maxNumParticles);
	m_particles.setMaxParticles(maxNumParticles);
}

void b3FluidSph::removeAllParticles()
{
	m_copyParticleDataToGpu = true;
	
	m_particles.resize(0);
	m_updates.clear();
	m_grid.clear();
}

void b3FluidSph::applyUpdates()
{
	m_copyParticleDataToGpu = true;

	//Create particles
	{
		int numCreatedParticles = m_updates.m_addedParticlePositions.size();
		printf("numCreatedParticles : %d \n", numCreatedParticles);
		for(int i = 0; i < numCreatedParticles; ++i)
		{
			int particleIndex = m_particles.addParticle(m_updates.m_addedParticlePositions[i]);
			if( particleIndex != numParticles() ) 
			{
				m_particles.m_velocity[particleIndex] = m_updates.m_addedParticleVelocities[i];
				m_particles.m_velocityEval[particleIndex] = m_updates.m_addedParticleVelocities[i];
			}
		}
	}
	
	//Set position and velocity
	{
		for(int i = 0; i < m_updates.m_updatedPositions.size(); ++i)
		{
			int particleIndex = m_updates.m_updatedPositionsIndex[i];
			m_particles.m_position[particleIndex] = m_updates.m_updatedPositions[i];
		}
		
		for(int i = 0; i < m_updates.m_updatedVelocities.size(); ++i)
		{
			int particleIndex = m_updates.m_updatedVelocitiesIndex[i];
			m_particles.m_velocity[particleIndex] = m_updates.m_updatedVelocities[i];
			m_particles.m_velocityEval[particleIndex] = m_updates.m_updatedVelocities[i];
		}
	}
	
	//Remove marked particles
	int numRemovedParticles = m_updates.m_removedParticleIndices.size();
	if(numRemovedParticles)
	{
		int numParticlesBeforeRemove = numParticles();
		int numParticlesAfterRemove = numParticlesBeforeRemove - numRemovedParticles;
		
		//Load indices into updates.m_removeSwapSourceCpu and updates.m_removeSwapTargetCpu
		m_updates.prepareToRemoveParticles(numParticlesBeforeRemove);
		
		for(int i = 0; i < numRemovedParticles; ++i)
		{
			int source = m_updates.m_removeSwapSourceCpu[i];	//Non-removed particle index
			int target = m_updates.m_removeSwapTargetCpu[i];	//Removed particle index
			
			m_particles.m_position[target] = m_particles.m_position[source];
			m_particles.m_velocity[target] = m_particles.m_velocity[source];
			m_particles.m_velocityEval[target] = m_particles.m_velocityEval[source];
			//m_particles.m_userPointer[target] = m_particles.m_userPointer[source];
			
			m_particles.m_accumulatedForce[target] = m_particles.m_accumulatedForce[source];
		}
		
		m_particles.resize(numParticlesAfterRemove);
	}
	
	m_updates.clear();
}

void b3FluidSph::insertParticlesIntoGrid()
{	
	m_copyParticleDataToGpu = true;
	
	B3_PROFILE("b3FluidSph::insertParticlesIntoGrid()");
	
	m_grid.setCellSize(m_parameters.m_simulationScale, m_parameters.m_sphSmoothRadius);
	
	//
	m_grid.clear();
	m_grid.insertParticles(m_particles);
}


// /////////////////////////////////////////////////////////////////////////////
// struct b3FluidEmitter
// /////////////////////////////////////////////////////////////////////////////
void b3FluidEmitter::emit()
{
	b3Assert(m_fluid);
	
	//m_particleIndicies.resize(0);
	
	if(!m_active) return;
	
	b3Quaternion rigidRotation = b3Quaternion::getIdentity();
	//if(m_attachTo) m_attachTo->getWorldTransform().getBasis().getRotation(rigidRotation);
	
	b3Vector3 velocity = b3QuatRotate(rigidRotation * m_rotation, m_direction) * m_speed;
	
	for(int i = 0; i < m_positions.size(); ++i)
	{
		b3Vector3 position = b3QuatRotate(m_rotation, m_positions[i]) + m_center;
		//if(m_attachTo) position = b3QuatRotate(rigidRotation, position) + m_attachTo->getWorldTransform().getOrigin();
		
		int index = m_fluid->addParticleCached(position, velocity);
		//if( index != m_fluid->numParticles() ) m_particleIndicies.push_back(index);
		/*
		if( index != m_fluid->numParticles() ) 
		{
			m_fluid->setVelocity(index, velocity);
			m_particleIndicies.push_back(index);
		}
		else if(m_useRandomIfAllParticlesAllocated)
		{
			index = ( m_fluid->numParticles() - 1 ) * GEN_rand() / GEN_RAND_MAX;		//Random index
		
			m_fluid->setPosition(index, position);
			m_fluid->setVelocity(index, velocity);
			m_particleIndicies.push_back(index);
		}*/
	}
}
void b3FluidEmitter::addVolume(b3FluidSph* fluid, const b3Vector3& min, const b3Vector3& max, b3Scalar spacing, int maxCreated)
{
	int numAdded = 0;

	for(b3Scalar z = max.getZ(); z >= min.getZ(); z -= spacing) 
		for(b3Scalar y = min.getY(); y <= max.getY(); y += spacing) 
			for(b3Scalar x = min.getX(); x <= max.getX(); x += spacing) 
			{
				fluid->addParticleCached( b3MakeVector3(x,y,z) );
				++numAdded;
				
				if(numAdded >= maxCreated) return;
			}
}
