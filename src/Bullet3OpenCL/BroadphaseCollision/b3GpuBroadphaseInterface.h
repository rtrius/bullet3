
#ifndef B3_GPU_BROADPHASE_INTERFACE_H
#define B3_GPU_BROADPHASE_INTERFACE_H

#include "Bullet3OpenCL/B3State/b3StateAabbs.h"
#include "Bullet3OpenCL/B3State/b3StateOverlappingPairs.h"

class b3GpuBroadphaseInterface
{
public:

	typedef class b3GpuBroadphaseInterface* (CreateFunc)(cl_context ctx,cl_device_id device, cl_command_queue q);

	virtual ~b3GpuBroadphaseInterface(){}

	virtual void computeOverlappingPairs(b3StateAabbs& input, b3StateOverlappingPairs& output, int maxPairs) = 0;

};

#endif //B3_GPU_BROADPHASE_INTERFACE_H
