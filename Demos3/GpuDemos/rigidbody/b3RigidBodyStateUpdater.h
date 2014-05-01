#ifndef B3_RIGID_BODY_STATE_UPDATER
#define B3_RIGID_BODY_STATE_UPDATER

class b3RigidBodyStateUpdater
{
	b3AlignedObjectArray<b3RigidBodyData> m_addedRigid;

	b3AlignedObjectArray<int> m_removedIndices;

	//Used only during applyUpdatesCpu()
	b3AlignedObjectArray<int> m_tempNewRigidIndices;
	b3AlignedObjectArray<b3InertiaData> m_tempAddedInertia;
	
public:
	b3RigidBodyStateUpdater() {}
	virtual ~b3RigidBodyStateUpdater() {}
	
	void addRigidBody(int collidableIndex, const b3Vector3& position, const b3Quaternion& orientation, b3Scalar mass,
						const b3Vector3& linearVelocity = b3MakeVector3(0,0,0), const b3Vector3& angularVelocity = b3MakeVector3(0,0,0) )
	{
		b3RigidBodyData rigid;
		
		rigid.m_pos = position;
		rigid.m_quat = orientation;
		rigid.m_angVel = linearVelocity;
		rigid.m_linVel = angularVelocity;
		
		rigid.m_collidableIdx = collidableIndex;
		rigid.m_invMass = (mass) ? b3Scalar(1.0) / mass : b3Scalar(0.0);
		rigid.m_restituitionCoeff = b3Scalar(0.0);
		rigid.m_frictionCoeff = b3Scalar(1.0);		//Same restitution/friction defaults in b3GpuNarrowPhase::registerRigidBody()
		
		m_addedRigid.push_back(rigid);
	}
	
	///@param rigidBodyIndex Must be an active rigid body index; the same index cannot be used twice.
	void markRigidBodyForRemoval(int rigidBodyIndex) { m_removedIndices.push_back(rigidBodyIndex); }
	
	void applyUpdatesCpu(b3GpuRigidBodyState* rigidState,
						b3GpuNarrowPhaseInternalData* narrowphaseInternalData)
	{
		b3AlignedObjectArray<b3RigidBodyData>* rigidBodiesCpu = narrowphaseInternalData->m_bodyBufferCPU;
		b3AlignedObjectArray<b3InertiaData>* inertiaCpu = narrowphaseInternalData->m_inertiaBufferCPU;
		b3AlignedObjectArray<int>* availableIndicesCpu = rigidState->m_availableRigidIndicesCPU;
		b3AlignedObjectArray<int>* usedIndicesCpu = rigidState->m_usedRigidIndicesCPU;
		
		b3OpenCLArray<b3RigidBodyData>* rigidBodiesGpu = narrowphaseInternalData->m_bodyBufferGPU;
		b3OpenCLArray<b3InertiaData>* inertiaGpu = narrowphaseInternalData->m_inertiaBufferGPU;
		b3OpenCLArray<int>* availableIndicesGpu = rigidState->m_availableRigidIndicesGPU;
		b3OpenCLArray<int>* usedIndicesGpu = rigidState->m_usedRigidIndicesGPU;	
		
		//Assuming that the GPU data is more current than CPU
		{
			B3_PROFILE("Transfer rigid body data to CPU");
			rigidBodiesGpu->copyToHost(*rigidBodiesCpu);
			inertiaGpu->copyToHost(*inertiaCpu);
			availableIndicesGpu->copyToHost(*availableIndicesCpu);
			usedIndicesGpu->copyToHost(*usedIndicesCpu);
		}
		
		//
		int numRemovedRigids = m_removedIndices.size();
		if(numRemovedRigids)
		{
			B3_PROFILE("Remove rigids from active array");
			
			//Update rigid body count
			{
				int numRemovedLargeRigids = 0;
				int numRemovedSmallRigids = 0;
				
				for(int i = 0; i < numRemovedRigids; ++i)
				{
					int removedRigidIndex = m_removedIndices[i];
					const b3RigidBodyData& removedRigid = (*rigidBodiesCpu)[removedRigidIndex];
					
					if( removedRigid.m_invMass == b3Scalar(0.0) ) ++numRemovedLargeRigids;
					else ++numRemovedSmallRigids;
				}
				
				rigidState->m_numLargeRigidBodies -= numRemovedLargeRigids;
				rigidState->m_numSmallRigidBodies -= numRemovedSmallRigids;
				narrowphaseInternalData->m_numAcceleratedRigidBodies -= numRemovedRigids;
			}
			
			//Update used indices
			{
				for(int i = 0; i < numRemovedRigids; ++i)
				{
					int indexToRemove = m_removedIndices[i];
					
					//Slow linear search results in worst case O(N^2)
					usedIndicesCpu->remove(indexToRemove);
					availableIndicesCpu->push_back(indexToRemove);
				}
				
				m_removedIndices.resize(0);
			}
		}
		
		
		//Add new rigids
		int numAddedRigids = m_addedRigid.size();
		if(numAddedRigids)
		{
			//Compute inertia tensor
			{
				m_tempAddedInertia.resize(numAddedRigids);
				
				for(int i = 0; i < numAddedRigids; ++i)
				{
					const b3RigidBodyData& rigid = m_addedRigid[i];
					
					b3InertiaData inertia;
					if( rigid.m_invMass == b3Scalar(0.0) )
					{
						inertia.m_initInvInertia.setValue(0,0,0, 0,0,0, 0,0,0);
						inertia.m_invInertiaWorld.setValue(0,0,0, 0,0,0, 0,0,0);
					}
					else
					{
						int collidableIndex = rigid.m_collidableIdx;
				
						//Approximate using inertia of a box(AABB)
						const b3SapAabb& aabb = (*narrowphaseInternalData->m_localShapeAABBCPU)[collidableIndex];
						b3Vector3 edge = aabb.m_maxVec - aabb.m_minVec;
						
						//For consistency with b3GpuNarrowPhase::registerRigidBody(), artifically increase the inertia
						edge *= b3Scalar(2.0);
						
						b3Vector3 localInertia = b3MakeVector3(edge.y * edge.y + edge.z * edge.z,
																edge.x * edge.x + edge.z * edge.z,
																edge.x * edge.x + edge.y * edge.y);
																
						b3Scalar mass = b3Scalar(1.0) / rigid.m_invMass;
						localInertia *= mass / b3Scalar(12.0);
						
						b3Vector3 invLocalInertia;
						invLocalInertia[0] = 1.f/localInertia[0];
						invLocalInertia[1] = 1.f/localInertia[1];
						invLocalInertia[2] = 1.f/localInertia[2];
						invLocalInertia[3] = 0.f;
						
						inertia.m_initInvInertia.setValue(
							invLocalInertia[0],		0,						0,
							0,						invLocalInertia[1],		0,
							0,						0,						invLocalInertia[2]);

						b3Matrix3x3 m (rigid.m_quat);

						inertia.m_invInertiaWorld = m.scaled(invLocalInertia) * m.transpose();
					}
					
					m_tempAddedInertia[i] = inertia;
				}
			}
			
			//Allocate new rigid indices
			{
				m_tempNewRigidIndices.resize(0);
				
				if( availableIndicesCpu->size() >= numAddedRigids )
				{
					//There are enough gaps in narrowphaseInternalData->m_bodyBufferCPU/GPU
					for(int i = 0; i < numAddedRigids; ++i)
					{
						m_tempNewRigidIndices.push_back( (*availableIndicesCpu)[ availableIndicesCpu->size() - 1 ] );
						availableIndicesCpu->pop_back();
					}
				}
				else
				{
					//Not enough empty rigids, need to resize the arrays(allocate completely new indices)
					m_tempNewRigidIndices = (*availableIndicesCpu);
					availableIndicesCpu->resize(0);
					
					int numNewRigidsInNewArraySlots = numAddedRigids - availableIndicesCpu->size();
					
					int newArraySize = narrowphaseInternalData->m_numAcceleratedRigidBodies + numNewRigidsInNewArraySlots;
					if( newArraySize > rigidBodiesCpu->size() )
					{
						rigidBodiesCpu->resize(newArraySize);
						inertiaCpu->resize(newArraySize);
					}
					
					int firstNewIndex = narrowphaseInternalData->m_numAcceleratedRigidBodies;
					for(int i = 0; i < numNewRigidsInNewArraySlots; ++i) m_tempNewRigidIndices.push_back(firstNewIndex + i);
				}
			}
			
			//Update arrays and rigid body count
			{
				int numNewLargeRigids = 0;
				int numNewSmallRigids = 0;
				
				for(int i = 0; i < numAddedRigids; ++i)
				{
					int rigidIndex = m_tempNewRigidIndices[i];
					(*rigidBodiesCpu)[rigidIndex] = m_addedRigid[i];
					(*inertiaCpu)[rigidIndex] = m_tempAddedInertia[i];
					
					if( m_addedRigid[i].m_invMass == b3Scalar(0.0) ) ++numNewLargeRigids;
					else ++numNewSmallRigids;
					
					usedIndicesCpu->push_back(rigidIndex);
				}
				
				rigidState->m_numLargeRigidBodies += numNewLargeRigids;
				rigidState->m_numSmallRigidBodies += numNewSmallRigids;
				narrowphaseInternalData->m_numAcceleratedRigidBodies += numAddedRigids;
			}
			
			m_addedRigid.resize(0);
		}
		
		//
		{
			B3_PROFILE("Transfer rigid body data to GPU");
			rigidBodiesGpu->copyFromHost(*rigidBodiesCpu);
			inertiaGpu->copyFromHost(*inertiaCpu);
			availableIndicesGpu->copyFromHost(*availableIndicesCpu);
			usedIndicesGpu->copyFromHost(*usedIndicesCpu);
		}
	}
};

#endif
