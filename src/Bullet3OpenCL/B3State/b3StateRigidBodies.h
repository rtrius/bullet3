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

#ifndef B3_STATE_RIGID_BODIES
#define B3_STATE_RIGID_BODIES

#include "Bullet3Collision/NarrowPhaseCollision/shared/b3RigidBodyData.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"

struct b3StateRigidBodies
{
	int m_maxRigidBodies;
	int m_numRigidBodies;

	b3OpenCLArray<b3RigidBodyData> m_bodyBufferGPU;
	b3OpenCLArray<b3InertiaData> m_inertiaBufferGPU;

	b3AlignedObjectArray<b3RigidBodyData> m_bodyBufferCPU;
	b3AlignedObjectArray<b3InertiaData> m_inertiaBufferCPU;


	b3StateRigidBodies(cl_context context, cl_command_queue queue, int maxRigidBodies);
	virtual ~b3StateRigidBodies();

	int registerRigidBody(int collidableIndex, float mass, const float* position, const float* orientation, const float* aabbMin, const float* aabbMax, bool writeToGpu);
	void setObjectTransform(const float* position, const float* orientation, int bodyIndex);

	void writeAllBodiesToGpu();
	void reset() { m_numRigidBodies = 0; }
	void readbackAllBodiesToCpu();
	bool getObjectTransformFromCpu(float* position, float* orientation, int bodyIndex) const;

	void setObjectTransformCpu(float* position, float* orientation, int bodyIndex);
	void setObjectVelocityCpu(float* linVel, float* angVel, int bodyIndex);
};

#endif
