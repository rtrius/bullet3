#ifndef B3_GPU_NARROWPHASE_H
#define B3_GPU_NARROWPHASE_H

#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLInclude.h"
#include "Bullet3Common/b3AlignedObjectArray.h"
#include "Bullet3Common/b3Vector3.h"

class b3GpuNarrowPhase
{
protected:

	struct b3GpuNarrowPhaseInternalData*	m_data;
	int m_acceleratedCompanionShapeIndex;
	int	m_static0Index;
	int m_maxTriConvexPairCapacity;

	cl_context m_context;
	cl_device_id m_device;
	cl_command_queue m_queue;


public:

	


	b3GpuNarrowPhase(cl_context vtx, cl_device_id dev, cl_command_queue q, const struct b3Config& config);

	virtual ~b3GpuNarrowPhase(void);

	void	writeAllBodiesToGpu();
	void  reset();
	void	readbackAllBodiesToCpu();

	virtual void computeContacts(cl_mem broadphasePairs, int numBroadphasePairs, cl_mem aabbsWorldSpace, int numObjects);
	

	int getStatic0Index() const
	{
		return m_static0Index;
	}
	const b3GpuNarrowPhaseInternalData*	getInternalData() const
	{
			return m_data;
	}

	b3GpuNarrowPhaseInternalData*	getInternalData()
	{
			return m_data;
	}

};

#endif //B3_GPU_NARROWPHASE_H

