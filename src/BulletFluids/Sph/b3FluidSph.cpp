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
//Portions of this file based on FLUIDS v.2 - SPH Fluid Simulator for CPU and GPU
//Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

#include "b3FluidSph.h"

#include "Bullet3Common/b3Logging.h"		//B3_PROFILE(name) macro
#include "Bullet3Common/b3Random.h"			//b3rand(), B3_RAND_MAX
#include "Bullet3Geometry/b3AabbUtil.h"		//b3TestPointAgainstAabb2()

#include "b3FluidSortingGrid.h"
//#include "b3FluidSphCollisionShape.h"

b3FluidSph::b3FluidSph(const b3FluidSphParametersGlobal& FG, int maxNumParticles)
{
	m_overrideSolver = 0;
	m_overrideParameters = 0;

	m_fluidDataCL = 0;
	m_gridDataCL = 0;
	
	setMaxParticles(maxNumParticles);
	setGridCellSize(FG);
	
	//b3CollisionObject
	///BULLET_2_TO_3_PLACEHOLDER
	/*
	{
		m_worldTransform.setIdentity();
		m_internalType = CO_USER_TYPE;	// replace later with CO_FLUID_SPH
		
		void* ptr = b3AlignedAlloc( sizeof(b3FluidSphCollisionShape), 16 );
		m_collisionShape = new(ptr) b3FluidSphCollisionShape(this);
		m_collisionShape->setMargin( b3Scalar(0.25) );	//Arbitrary value
		
		m_rootCollisionShape = m_collisionShape;
	}
	*/
}
b3FluidSph::~b3FluidSph()
{
	//b3CollisionObject
	///BULLET_2_TO_3_PLACEHOLDER
	/*
	{
		m_collisionShape->~b3CollisionShape();
		b3AlignedFree(m_collisionShape);
	}
	*/
}

void b3FluidSph::setGridCellSize(const b3FluidSphParametersGlobal& FG)
{
	m_grid.setCellSize(FG.m_simulationScale, FG.m_sphSmoothRadius);
}

void b3FluidSph::setMaxParticles(int maxNumParticles)
{
	if( maxNumParticles < m_particles.size() )m_particles.resize(maxNumParticles);
	m_particles.setMaxParticles(maxNumParticles);
}

void b3FluidSph::removeAllParticles()
{
	m_particles.resize(0);
	
	m_removedFluidIndicies.resize(0);
	
	m_grid.clear();
}

b3Scalar b3FluidSph::getValue(b3Scalar x, b3Scalar y, b3Scalar z) const
{
	const b3Scalar worldSphRadius = m_grid.getCellSize();	//Grid cell size == sph interaction radius, at world scale
	const b3Scalar R2 = worldSphRadius * worldSphRadius;
	
	b3FluidSortingGrid::FoundCells foundCells;
	m_grid.findCells( b3MakeVector3(x,y,z), foundCells );
		
	b3Scalar sum = 0.0;
	for(int cell = 0; cell < b3FluidSortingGrid::NUM_FOUND_CELLS; cell++) 
	{
		b3FluidGridIterator& FI = foundCells.m_iterators[cell];
			
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			const b3Vector3& position = m_particles.m_pos[n];
			b3Scalar dx = x - position.getX();
			b3Scalar dy = y - position.getY();
			b3Scalar dz = z - position.getZ();
			b3Scalar distanceSquared = dx*dx + dy*dy + dz*dz;
				
			if(distanceSquared < R2) sum += R2 / distanceSquared;
		}
	}
	
	return sum;
}	
b3Vector3 b3FluidSph::getGradient(b3Scalar x, b3Scalar y, b3Scalar z) const
{
	const b3Scalar worldSphRadius = m_grid.getCellSize();	//Grid cell size == sph interaction radius, at world scale
	const b3Scalar R2 = worldSphRadius*worldSphRadius;
	
	b3FluidSortingGrid::FoundCells foundCells;
	m_grid.findCells( b3MakeVector3(x,y,z), foundCells );
	
	b3Vector3 normal = b3MakeVector3(0,0,0);
	for(int cell = 0; cell < b3FluidSortingGrid::NUM_FOUND_CELLS; cell++)
	{
		b3FluidGridIterator& FI = foundCells.m_iterators[cell];
			
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			const b3Vector3& position = m_particles.m_pos[n];
			b3Scalar dx = x - position.getX();
			b3Scalar dy = y - position.getY();
			b3Scalar dz = z - position.getZ();
			b3Scalar distanceSquared = dx*dx + dy*dy + dz*dz;
			
			if( b3Scalar(0.0) < distanceSquared && distanceSquared < R2 ) 
			{
				distanceSquared = b3Scalar(2.0)*R2 / (distanceSquared*distanceSquared);
				
				b3Vector3 particleNorm = b3MakeVector3(dx * distanceSquared, dy * distanceSquared, dz * distanceSquared);
				normal += particleNorm;
			}
		}
	}
	
	normal.normalize();
	return normal;
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
void b3FluidSph::removeMarkedParticles()
{
	//makeUnique() assumes that the array is sorted
	//m_removedFluidIndicies.heapSort( AscendingSortPredicate() );
	m_removedFluidIndicies.quickSort( AscendingSortPredicate() );
	
	//Remove duplicate indicies
	makeUniqueInt(m_removedFluidIndicies);
	
	//Since removing elements from the array invalidates(higher) indicies,
	//elements should be removed in descending order.
	for(int i = m_removedFluidIndicies.size() - 1; i >= 0; --i) m_particles.removeParticle( m_removedFluidIndicies[i] );
	
	m_removedFluidIndicies.resize(0);
}

void b3FluidSph::insertParticlesIntoGrid()
{	
	B3_PROFILE("b3FluidSph::insertParticlesIntoGrid()");
	
	//
	m_grid.clear();
	m_grid.insertParticles(m_particles);
}


// /////////////////////////////////////////////////////////////////////////////
// struct b3FluidEmitter
// /////////////////////////////////////////////////////////////////////////////
void b3FluidEmitter::emit(b3FluidSph* fluid, int numParticles, b3Scalar spacing)
{
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
}
void b3FluidEmitter::addVolume(b3FluidSph* fluid, const b3Vector3& min, const b3Vector3& max, b3Scalar spacing)
{
	for(b3Scalar z = max.getZ(); z >= min.getZ(); z -= spacing) 
		for(b3Scalar y = min.getY(); y <= max.getY(); y += spacing) 
			for(b3Scalar x = min.getX(); x <= max.getX(); x += spacing) 
			{
				fluid->addParticle( b3MakeVector3(x,y,z) );
			}
}


// /////////////////////////////////////////////////////////////////////////////
// struct b3FluidAbsorber
// /////////////////////////////////////////////////////////////////////////////
struct b3FluidAbsorberCallback : public b3FluidSortingGrid::AabbCallback
{
	b3FluidSph* m_fluidSph;

	const b3Vector3& m_min;
	const b3Vector3& m_max;

	b3FluidAbsorberCallback(b3FluidSph* fluidSph, const b3Vector3& min, const b3Vector3& max) 
	: m_fluidSph(fluidSph), m_min(min), m_max(max) {}
	
	virtual bool processParticles(const b3FluidGridIterator FI, const b3Vector3& aabbMin, const b3Vector3& aabbMax)
	{
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			if( b3TestPointAgainstAabb2( m_min, m_max, m_fluidSph->getPosition(n) ) ) m_fluidSph->markParticleForRemoval(n);
		}
	
		return true;
	}
};
void b3FluidAbsorber::absorb(b3FluidSph* fluid)
{
	const b3FluidSortingGrid& grid = fluid->getGrid();
	
	b3FluidAbsorberCallback absorber(fluid, m_min, m_max);
	grid.forEachGridCell(m_min, m_max, absorber);
}
