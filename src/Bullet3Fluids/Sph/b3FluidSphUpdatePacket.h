#ifndef B3_FLUID_SPH_UPDATE_PACKET
#define B3_FLUID_SPH_UPDATE_PACKET

#include "Bullet3Common/b3AlignedObjectArray.h"
#include "Bullet3Common/b3Vector3.h"

///@brief Allows incremental update of particle data on the GPU.
///@remarks
///Note that the indices of particles change from frame to frame - make sure that the indices are correct before calling any functions.
///The update process occurs as follows:
/// - First, added particles are created.
/// - Next, the position and velocity is set.
/// - Last, particles are removed.
struct b3FluidSphUpdatePacket
{
	///The lengths of each pair of arrays is equal, for example m_updatedPositionsIndex.size() == m_updatedPositions.size()
	b3AlignedObjectArray<int> m_updatedPositionsIndex;
	b3AlignedObjectArray<b3Vector3> m_updatedPositions;
	
	b3AlignedObjectArray<int> m_updatedVelocitiesIndex;
	b3AlignedObjectArray<b3Vector3> m_updatedVelocities;
	
	b3AlignedObjectArray<b3Vector3> m_addedParticlePositions;
	b3AlignedObjectArray<b3Vector3> m_addedParticleVelocities;
	
	//
	b3AlignedObjectArray<int> m_removedParticleIndices;
	
	enum UpdateType
	{
		UT_Position,	///<Set the position of a particle
		UT_Velocity,	///<Set the velocity of a particle
		//UT_Force,		///<Set the force on a particle(total force, does not accumulate with previously applied force)
	};

	void clear()
	{
		m_updatedPositionsIndex.resize(0);
		m_updatedPositions.resize(0);
		
		m_updatedVelocitiesIndex.resize(0);
		m_updatedVelocities.resize(0);
		
		m_addedParticlePositions.resize(0);
		m_addedParticleVelocities.resize(0);
		
		m_removedParticleIndices.resize(0);
	}
	
	///The returned index is only valid on the frame that the particle was created
	int addParticle(int maxParticles, int numParticles, const b3Vector3& position, const b3Vector3& velocity) 
	{
		int newIndex = numParticles + m_addedParticlePositions.size() + 1;
		if(newIndex >= maxParticles) return numParticles;
	
		m_addedParticlePositions.push_back(position);
		m_addedParticleVelocities.push_back(velocity);
		
		return newIndex;
	}
	
	///Do not add multiple updates with the same index; threads will collide resuting in unpredictable position/velocity
	void addUpdate(b3FluidSphUpdatePacket::UpdateType type, int index, const b3Vector3& positionVelocityOrForce)
	{
		switch(type)
		{
			case b3FluidSphUpdatePacket::UT_Position:
				m_updatedPositionsIndex.push_back(index);
				m_updatedPositions.push_back(positionVelocityOrForce);
				break;
			case b3FluidSphUpdatePacket::UT_Velocity:
				m_updatedVelocitiesIndex.push_back(index);
				m_updatedVelocities.push_back(positionVelocityOrForce);
				break;
				
			default:
				b3Assert(0);
				break;
		}
	}
	
	///Duplicate indices are ignored, so a particle may be marked twice without any issues
	void markParticleForRemoval(int index) { m_removedParticleIndices.push_back(index); }
};

#endif
