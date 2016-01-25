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

#include "b3GpuParallelLinearBvhBroadphase.h"


b3GpuParallelLinearBvhBroadphase::b3GpuParallelLinearBvhBroadphase(cl_context context, cl_device_id device, cl_command_queue queue) :
m_plbvh(context, device, queue)
{
}

void b3GpuParallelLinearBvhBroadphase::computeOverlappingPairs(b3StateAabbs& input, b3StateOverlappingPairs& output, int maxPairs)
{
	//Reconstruct BVH
	m_plbvh.build(input.m_aabbsGpu, input.m_smallAabbsMappingGpu, input.m_largeAabbsMappingGpu);

	//
	output.m_overlappingPairsGpu.resize(maxPairs);
	m_plbvh.calculateOverlappingPairs(output.m_overlappingPairsGpu);
}
