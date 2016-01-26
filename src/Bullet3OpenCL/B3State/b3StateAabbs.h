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
#ifndef B3_STATE_AABBS
#define B3_STATE_AABBS

#include "Bullet3Common/shared/b3Int4.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3SapAabb.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"

struct b3StateAabbs
{
	b3OpenCLArray<b3SapAabb> m_aabbsGpu;
	b3OpenCLArray<int> m_smallAabbsMappingGpu;
	b3OpenCLArray<int> m_largeAabbsMappingGpu;
	
	b3AlignedObjectArray<b3SapAabb> m_aabbsCpu;
	b3AlignedObjectArray<int> m_smallAabbsMappingCpu;
	b3AlignedObjectArray<int> m_largeAabbsMappingCpu;

	b3StateAabbs(cl_context context, cl_command_queue queue);
	virtual ~b3StateAabbs() {}
	
	virtual void createProxy(const b3Vector3& aabbMin, const b3Vector3& aabbMax, int userPtr, 
							short int collisionFilterGroup, short int collisionFilterMask);
	virtual void createLargeProxy(const b3Vector3& aabbMin, const b3Vector3& aabbMax, int userPtr, 
							short int collisionFilterGroup, short int collisionFilterMask);
	
	//call writeAabbsToGpu after done making all changes (createProxy etc)
	virtual void writeAabbsToGpu();
};

#endif
