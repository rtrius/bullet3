#ifndef B3_GPU_NARROWPHASE_H
#define B3_GPU_NARROWPHASE_H

#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLInclude.h"
#include "Bullet3Common/b3AlignedObjectArray.h"
#include "Bullet3Common/b3Vector3.h"

#include "Bullet3OpenCL/B3State/b3StateAabbs.h"
#include "Bullet3OpenCL/B3State/b3StateOverlappingPairs.h"
#include "Bullet3OpenCL/B3State/b3StateRigidBodies.h"
#include "Bullet3OpenCL/B3State/b3StateRigidCollidables.h"
#include "Bullet3OpenCL/B3State/b3StateRigidContacts.h"

class b3GpuNarrowPhase
{
protected:

	struct b3GpuNarrowPhaseInternalData*	m_data;
	int m_acceleratedCompanionShapeIndex;
	int m_maxTriConvexPairCapacity;

	cl_context m_context;
	cl_device_id m_device;
	cl_command_queue m_queue;


public:
	b3GpuNarrowPhase(cl_context vtx, cl_device_id dev, cl_command_queue q, const struct b3Config& config);

	virtual ~b3GpuNarrowPhase(void);

	//virtual void computeContacts(cl_mem broadphasePairs, int numBroadphasePairs, cl_mem aabbsWorldSpace, int numObjects);
	virtual void computeContacts(b3StateAabbs& worldSpaceAabbs, b3StateOverlappingPairs& pairs, 
								b3StateRigidBodies& rigids, b3StateRigidCollidables& collidables,
								b3StateRigidContacts& contacts);


};

#endif //B3_GPU_NARROWPHASE_H

