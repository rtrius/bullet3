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
	
	m_copyParticleDataToGpu = true;
	
	setMaxParticles(maxNumParticles);
	m_grid.setCellSize(m_parameters.m_simulationScale, m_parameters.m_sphSmoothRadius);
}
b3FluidSph::~b3FluidSph()
{
}

void b3FluidSph::setMaxParticles(int maxNumParticles)
{
	if( maxNumParticles < m_particles.size() )m_particles.resize(maxNumParticles);
	m_particles.setMaxParticles(maxNumParticles);
}

void b3FluidSph::removeAllParticles()
{
	m_particles.resize(0);
	
	m_updates.clear();
	
	m_grid.clear();
}

//Assumes that out_unique is already sorted.
//Removes duplicates; rearranges array such that all unique values are in the range [0, uniqueSize).
void makeUniqueInt(b3AlignedObjectArray<int>& out_unique)
{
	int uniqueSize = 0;
	if( out_unique.size() ) 
	{
		uniqueSize = 1;
		for(int i = 1; i < out_unique.size(); ++i)
		{
			if( out_unique[i] != out_unique[i-1] )
			{
				out_unique[uniqueSize] = out_unique[i];
				++uniqueSize;
			}
		}
	}
	
	out_unique.resize(uniqueSize);
}
struct AscendingSortPredicate { inline bool operator() (const int& a, const int& b) const { return (a < b); } };
void b3FluidSph::applyUpdates()
{
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
	{
		//makeUnique() assumes that the array is sorted
		//m_updates.m_removedParticleIndices.heapSort( AscendingSortPredicate() );
		m_updates.m_removedParticleIndices.quickSort( AscendingSortPredicate() );
		
		//Remove duplicate indicies
		makeUniqueInt(m_updates.m_removedParticleIndices);
		
		//Since removing elements from the array invalidates(higher) indicies,
		//elements should be removed in descending order.
		for(int i = m_updates.m_removedParticleIndices.size() - 1; i >= 0; --i) 
			m_particles.removeParticle( m_updates.m_removedParticleIndices[i] );
	}
	
	m_updates.clear();
}

void b3FluidSph::insertParticlesIntoGrid()
{	
	B3_PROFILE("b3FluidSph::insertParticlesIntoGrid()");
	
	m_grid.setCellSize(m_parameters.m_simulationScale, m_parameters.m_sphSmoothRadius);
	
	//
	m_grid.clear();
	m_grid.insertParticles(m_particles);
}


// /////////////////////////////////////////////////////////////////////////////
// struct b3FluidEmitter
// /////////////////////////////////////////////////////////////////////////////
void b3FluidEmitter::emit(b3FluidSph* fluid, int numParticles, b3Scalar spacing)
{
	/*
	int x = static_cast<int>( b3Sqrt(static_cast<b3Scalar>(numParticles)) );
	
	for(int i = 0; i < numParticles; i++) 
	{
		b3Scalar ang_rand = ( static_cast<b3Scalar>(b3rand()*b3Scalar(2.0)/B3_RAND_MAX) - b3Scalar(1.0) ) * m_yawSpread;
		b3Scalar tilt_rand = ( static_cast<b3Scalar>(b3rand()*b3Scalar(2.0)/B3_RAND_MAX) - b3Scalar(1.0) ) * m_pitchSpread;
		
		b3Vector3 dir = b3MakeVector3( 	b3Cos((m_yaw + ang_rand) * B3_RADS_PER_DEG) * b3Sin((m_pitch + tilt_rand) * B3_RADS_PER_DEG) * m_velocity,
										b3Cos((m_pitch + tilt_rand) * B3_RADS_PER_DEG) * m_velocity,
										b3Sin((m_yaw + ang_rand) * B3_RADS_PER_DEG) * b3Sin((m_pitch + tilt_rand) * B3_RADS_PER_DEG) * m_velocity );
		
		b3Vector3 position = b3MakeVector3( spacing*(i/x), spacing*(i%x), 0 );
		position += m_position;
		
		int index = fluid->addParticle(position);
		
		if( index != fluid->numParticles() ) fluid->setVelocity(index, dir);
		else if(m_useRandomIfAllParticlesAllocated)
		{
			index = ( fluid->numParticles() - 1 ) * b3rand() / B3_RAND_MAX;		//Random index
		
			fluid->setPosition(index, position);
			fluid->setVelocity(index, dir);
		}
	}
	*/
}
void b3FluidEmitter::addVolume(b3FluidSph* fluid, const b3Vector3& min, const b3Vector3& max, b3Scalar spacing)
{
	for(b3Scalar z = max.getZ(); z >= min.getZ(); z -= spacing) 
		for(b3Scalar y = min.getY(); y <= max.getY(); y += spacing) 
			for(b3Scalar x = min.getX(); x <= max.getX(); x += spacing) 
			{
				fluid->addParticleCached( b3MakeVector3(x,y,z) );
			}
}
