#ifndef B3_FLUID_SPH_UPDATE_PACKET_H
#define B3_FLUID_SPH_UPDATE_PACKET_H

#include "Bullet3Common/b3AlignedObjectArray.h"
#include "Bullet3Common/b3Vector3.h"

///@brief Allows incremental update of particle data on the GPU.
///@remarks
///Note that the indices of particles change from frame to frame - make sure that the indices are correct before calling any functions.
///The update process occurs as follows:
/// - First, added particles are created.
/// - Next, the position and velocity is set.
/// - Last, particles are removed.
///@par 
///The update data still has to be transferred to the GPU every frame, so creating
///a large number of updates will decrease performance.
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
	
	//
	b3AlignedObjectArray<int> m_removeSwapSourceCpu;
	b3AlignedObjectArray<int> m_removeSwapTargetCpu;

	void clear()
	{
		m_updatedPositionsIndex.resize(0);
		m_updatedPositions.resize(0);
		
		m_updatedVelocitiesIndex.resize(0);
		m_updatedVelocities.resize(0);
		
		m_addedParticlePositions.resize(0);
		m_addedParticleVelocities.resize(0);
		
		m_removedParticleIndices.resize(0);
		
		m_removeSwapSourceCpu.resize(0);
		m_removeSwapTargetCpu.resize(0);
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
	
	///Do not add multiple updates with the same index; threads will collide resulting in unpredictable position
	void setPosition(int index, const b3Vector3& position)
	{
		m_updatedPositionsIndex.push_back(index);
		m_updatedPositions.push_back(position);
	}
	
	///Do not add multiple updates with the same index; threads will collide resulting in unpredictable velocity
	void setVelocity(int index, const b3Vector3& velocity)
	{
		m_updatedVelocitiesIndex.push_back(index);
		m_updatedVelocities.push_back(velocity);
	}
	
	///Duplicate indices are ignored, so a particle may be marked twice without any issues
	void markParticleForRemoval(int index) { m_removedParticleIndices.push_back(index); }
	
	//Load m_removeSwapSourceCpu and m_removeSwapTargetCpu
	//To remove particles, copy particle attributes(position, velocity, etc.)
	//from indicies in m_removeSwapSourceCpu to indicies in m_removeSwapTargetCpu,
	//then shrink the arrays by the number of particles removed.
	void prepareToRemoveParticles(int numParticlesBeforeRemove)
	{
		//Process for removing particles:
		//For instance there can be the array of particles:
		// N R N R N N R N N R
		//Where N is an index that is not marked to be removed, and R is marked to be removed.
		//In this case there are 10 particles, with 6 marked N and 4 marked R.
		//
		//Split the array at the last index after particles are removed(6 particles remain, so split occurs at 6th particle).
		// N R N R N N / R N N R
		//The main point is to notice that when the array is divided at the number of particles remaining,
		//the number of particles marked R to the left of the split is always the same as the number of particles marked N to the right.
		//
		//As a result, we can swap the particles marked N on the right with the particles marked R on the left.
		//
		//       |-Swap 3, 8-| 
		//       |           |
		// N R N R N N / R N N R        <-- 2 R on left, 2 N on right
		// 0 1 2 3 4 5   6 7 8 9
		//   |             |
		//   |- Swap 1, 7 -|
		//
		// N N N N N N / R R R R
		//and finally resize/truncate the array to eliminate the removed particles.
		//
		//In order to ensure that the result is deterministic, the array of marked particles
		//is sorted in ascending order and the source and target arrays for the swap are also 
		//sorted in ascending order.
		//
		//In the above example, 
		//    m_removeSwapSourceCpu == [7, 8], and
		//    m_removeSwapTargetCpu == [1, 3]
		//note ascending order of both arrays.
		
		int numRemovedParticles = m_removedParticleIndices.size();
		int numParticlesAfterRemove = numParticlesBeforeRemove - numRemovedParticles;
		
		//makeUniqueInt() assumes that the array is sorted
		m_removedParticleIndices.quickSort( b3FluidSphUpdatePacket::AscendingSortPredicate() );
	
		//Remove duplicate indicies
		b3FluidSphUpdatePacket::makeUniqueInt(m_removedParticleIndices);
		
		//Find indices to the right of the split that are marked N(that is, not marked for remove) and use them as the swap source.
		{
			m_removeSwapSourceCpu.resize(0);
		
			int removedIndex = 0;
			for(int i = numParticlesAfterRemove; i < numParticlesBeforeRemove; ++i) 
			{
				while( i > m_removedParticleIndices[removedIndex] && removedIndex < m_removedParticleIndices.size() ) ++removedIndex;
			
				if(i != m_removedParticleIndices[removedIndex]) m_removeSwapSourceCpu.push_back(i);
			}
		}
		
		//Find indices to the left of the split that are marked R(remove), and use them as the swap target.
		m_removeSwapTargetCpu.resize(0);
		for(int i = 0; i < numRemovedParticles; ++i)
		{
			if(m_removedParticleIndices[i] < numParticlesAfterRemove) m_removeSwapTargetCpu.push_back(m_removedParticleIndices[i]);
			else break;
		}
	}
	
	//Assumes that out_unique is already sorted.
	//Removes duplicates; rearranges array such that all unique values are in the range [0, uniqueSize).
	static void makeUniqueInt(b3AlignedObjectArray<int>& out_unique)
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
};

#endif
