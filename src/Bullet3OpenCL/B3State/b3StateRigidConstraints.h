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

#ifndef B3_STATE_RIGID_CONSTRAINTS
#define B3_STATE_RIGID_CONSTRAINTS

#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3OpenCL/RigidBody/b3GpuGenericConstraint.h"

struct b3StateRigidConstraints
{
	int	m_constraintUid;

	b3OpenCLArray<b3GpuGenericConstraint> m_constraintsGpu;
	b3AlignedObjectArray<b3GpuGenericConstraint> m_constraintsCpu;

	b3StateRigidConstraints(cl_context context, cl_command_queue queue);
	virtual ~b3StateRigidConstraints();

	int createPoint2PointConstraint(int bodyA, int bodyB, const float* pivotInA, const float* pivotInB, float breakingThreshold);
	int createFixedConstraint(int bodyA, int bodyB, const float* pivotInA, const float* pivotInB, const float* relTargetAB, float breakingThreshold);

	void removeConstraintByUid(int uid);

	void reset();

	void writeConstraintsToGpu();
	void copyConstraintsToHost();
};

#endif

