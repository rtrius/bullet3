

#include "Bullet3Collision/NarrowPhaseCollision/shared/b3UpdateAabbs.h"


__kernel void initializeGpuAabbsFull(const int numActiveNodes, 
									__global int* activeNodeIndices, 
									__global b3RigidBodyData_t* gBodies,
									__global b3Collidable_t* collidables,
									__global b3Aabb_t* plocalShapeAABB,
									__global b3Aabb_t* pAABB)
{
	int nodeID = get_global_id(0);
	if( nodeID < numActiveNodes )
	{
		b3ComputeWorldAabb(activeNodeIndices[nodeID], nodeID, gBodies, collidables, plocalShapeAABB, pAABB);
	}
}

__kernel void loadLargeAndSmallAabbIndices(__global const int* activeNodeIndices, 
										__global const b3RigidBodyData_t* rigidBodyData,
										
										__global int* numLargeAabbs, 
										__global int* numSmallAabbs, 
										__global int* largeAabbIndices, 
										__global int* smallAabbIndices, 
										int maxLargeAabbs, int maxSmallAabbs, int numActiveNodes)
{
	int activeIndex = get_global_id(0);
	if(activeIndex >= numActiveNodes) return;

	int rigidBodyIndex = activeNodeIndices[activeIndex];

	int isLargeAabb = (rigidBodyData[rigidBodyIndex].m_invMass == 0.f);
	int maxAabbs = (isLargeAabb) ? maxLargeAabbs : maxSmallAabbs;
	__global int* numAabbs = (isLargeAabb) ? numLargeAabbs : numSmallAabbs;
	__global int* aabbIndices = (isLargeAabb) ?  largeAabbIndices : smallAabbIndices;
	
	int aabbIndex = atomic_inc(numAabbs);
	if(aabbIndex < maxAabbs) aabbIndices[aabbIndex] = activeIndex;
}

__kernel void clearOverlappingPairsKernel(  __global int4* pairs, int numPairs)
{
	int pairId = get_global_id(0);
	if( pairId< numPairs )
	{
		pairs[pairId].z = 0xffffffff;
	}
}