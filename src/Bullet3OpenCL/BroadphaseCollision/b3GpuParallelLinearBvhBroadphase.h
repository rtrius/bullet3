/*
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef B3_GPU_PARALLEL_LINEAR_BVH_BROADPHASE_H
#define B3_GPU_PARALLEL_LINEAR_BVH_BROADPHASE_H

#include "Bullet3OpenCL/BroadphaseCollision/b3GpuBroadphaseInterface.h"

#include "b3GpuParallelLinearBvh.h"

class b3GpuParallelLinearBvhBroadphase : public b3GpuBroadphaseInterface
{
	b3GpuParallelLinearBvh m_plbvh;

public:
	b3GpuParallelLinearBvhBroadphase(cl_context context, cl_device_id device, cl_command_queue queue);
	virtual ~b3GpuParallelLinearBvhBroadphase() {}

	virtual void computeOverlappingPairs(b3StateAabbs& input, b3StateOverlappingPairs& output, int maxPairs);

	static b3GpuBroadphaseInterface* CreateFunc(cl_context context, cl_device_id device, cl_command_queue queue)
	{
		return new b3GpuParallelLinearBvhBroadphase(context, device, queue);
	}
};

#endif
