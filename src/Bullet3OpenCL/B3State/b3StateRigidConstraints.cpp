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
#include "b3StateRigidConstraints.h"

b3StateRigidConstraints::b3StateRigidConstraints(cl_context context, cl_command_queue queue) :
	m_constraintUid(0),
	m_constraintsGpu(context, queue)
{

}

b3StateRigidConstraints::~b3StateRigidConstraints()
{

}


int b3StateRigidConstraints::createPoint2PointConstraint(int bodyA, int bodyB, const float* pivotInA, const float* pivotInB, float breakingThreshold)
{
	b3GpuGenericConstraint c;
	c.m_uid = m_constraintUid++;
	c.m_flags = B3_CONSTRAINT_FLAG_ENABLED;
	c.m_rbA = bodyA;
	c.m_rbB = bodyB;
	c.m_pivotInA.setValue(pivotInA[0], pivotInA[1], pivotInA[2]);
	c.m_pivotInB.setValue(pivotInB[0], pivotInB[1], pivotInB[2]);
	c.m_breakingImpulseThreshold = breakingThreshold;
	c.m_constraintType = B3_GPU_POINT2POINT_CONSTRAINT_TYPE;

	m_constraintsCpu.push_back(c);
	return c.m_uid;
}
int b3StateRigidConstraints::createFixedConstraint(int bodyA, int bodyB, const float* pivotInA, const float* pivotInB, const float* relTargetAB, float breakingThreshold)
{
	b3GpuGenericConstraint c;
	c.m_uid = m_constraintUid++;
	c.m_flags = B3_CONSTRAINT_FLAG_ENABLED;
	c.m_rbA = bodyA;
	c.m_rbB = bodyB;
	c.m_pivotInA.setValue(pivotInA[0], pivotInA[1], pivotInA[2]);
	c.m_pivotInB.setValue(pivotInB[0], pivotInB[1], pivotInB[2]);
	c.m_relTargetAB.setValue(relTargetAB[0], relTargetAB[1], relTargetAB[2], relTargetAB[3]);
	c.m_breakingImpulseThreshold = breakingThreshold;
	c.m_constraintType = B3_GPU_FIXED_CONSTRAINT_TYPE;

	m_constraintsCpu.push_back(c);
	return c.m_uid;
}



void  b3StateRigidConstraints::removeConstraintByUid(int uid)
{
	//slow linear search
	m_constraintsGpu.copyToHost(m_constraintsCpu);

	//remove
	for (int i = 0; i<m_constraintsCpu.size(); i++)
	{
		if (m_constraintsCpu[i].m_uid == uid)
		{
			m_constraintsCpu.swap(i, m_constraintsCpu.size() - 1);
			m_constraintsCpu.pop_back();

			break;
		}
	}

	if (m_constraintsCpu.size())
	{
		m_constraintsGpu.copyFromHost(m_constraintsCpu);
	}
	else
	{
		m_constraintsGpu.resize(0);
	}

}


void b3StateRigidConstraints::reset()
{
	m_constraintsGpu.resize(0);
	m_constraintsCpu.resize(0);
}


void b3StateRigidConstraints::copyConstraintsToHost()
{
	m_constraintsGpu.copyToHost(m_constraintsCpu);
}


void b3StateRigidConstraints::writeConstraintsToGpu()
{
	m_constraintsGpu.copyFromHost(m_constraintsCpu);
}