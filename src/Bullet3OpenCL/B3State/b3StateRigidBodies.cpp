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
#include "b3StateRigidBodies.h"


b3StateRigidBodies::b3StateRigidBodies(cl_context context, cl_command_queue queue, int maxRigidBodies) :
	m_maxRigidBodies(maxRigidBodies),
	m_numRigidBodies(0),
	m_bodyBufferGPU(context, queue), 
	m_inertiaBufferGPU(context, queue)
{
	m_bodyBufferGPU.resize(maxRigidBodies);
	m_inertiaBufferGPU.resize(maxRigidBodies);
	m_bodyBufferCPU.resize(maxRigidBodies);
	m_inertiaBufferCPU.resize(maxRigidBodies);
}
b3StateRigidBodies::~b3StateRigidBodies()
{
}

int b3StateRigidBodies::registerRigidBody(int collidableIndex, float mass, const float* position, const float* orientation, const float* aabbMinPtr, const float* aabbMaxPtr, bool writeToGpu)
{
	b3Vector3 aabbMin = b3MakeVector3(aabbMinPtr[0], aabbMinPtr[1], aabbMinPtr[2]);
	b3Vector3 aabbMax = b3MakeVector3(aabbMaxPtr[0], aabbMaxPtr[1], aabbMaxPtr[2]);


	if (m_numRigidBodies >= m_maxRigidBodies)
	{
		b3Error("registerRigidBody: exceeding the number of rigid bodies, %d > %d \n", m_numRigidBodies, m_maxRigidBodies);
		return -1;
	}

	m_bodyBufferCPU.resize(m_numRigidBodies + 1);

	b3RigidBodyData& body = m_bodyBufferCPU[m_numRigidBodies];

	float friction = 1.f;
	float restitution = 0.f;

	body.m_frictionCoeff = friction;
	body.m_restituitionCoeff = restitution;
	body.m_angVel = b3MakeVector3(0, 0, 0);
	body.m_linVel = b3MakeVector3(0, 0, 0);
	body.m_pos = b3MakeVector3(position[0], position[1], position[2]);
	body.m_quat.setValue(orientation[0], orientation[1], orientation[2], orientation[3]);
	body.m_collidableIdx = collidableIndex;

	body.m_invMass = mass ? 1.f / mass : 0.f;

	if (writeToGpu)
	{
		m_bodyBufferGPU.copyFromHostPointer(&body, 1, m_numRigidBodies);
	}

	b3InertiaData& shapeInfo = m_inertiaBufferCPU[m_numRigidBodies];

	if (mass == 0.f)
	{
		shapeInfo.m_initInvInertia.setValue(0, 0, 0, 0, 0, 0, 0, 0, 0);
		shapeInfo.m_invInertiaWorld.setValue(0, 0, 0, 0, 0, 0, 0, 0, 0);
	}
	else
	{
		b3Assert(body.m_collidableIdx >= 0);

		//approximate using the aabb of the shape

		b3Vector3 halfExtents = (aabbMax - aabbMin);//*0.5f;//fake larger inertia makes demos more stable ;-)

		b3Vector3 localInertia;

		float lx = 2.f*halfExtents[0];
		float ly = 2.f*halfExtents[1];
		float lz = 2.f*halfExtents[2];

		localInertia.setValue((mass / 12.0f) * (ly*ly + lz*lz),
			(mass / 12.0f) * (lx*lx + lz*lz),
			(mass / 12.0f) * (lx*lx + ly*ly));

		b3Vector3 invLocalInertia;
		invLocalInertia[0] = 1.f / localInertia[0];
		invLocalInertia[1] = 1.f / localInertia[1];
		invLocalInertia[2] = 1.f / localInertia[2];
		invLocalInertia[3] = 0.f;

		shapeInfo.m_initInvInertia.setValue(
			invLocalInertia[0], 0, 0,
			0, invLocalInertia[1], 0,
			0, 0, invLocalInertia[2]);

		b3Matrix3x3 m(body.m_quat);

		shapeInfo.m_invInertiaWorld = m.scaled(invLocalInertia) * m.transpose();
	}

	if (writeToGpu)
		m_inertiaBufferGPU.copyFromHostPointer(&shapeInfo, 1, m_numRigidBodies);


	return m_numRigidBodies++;
}

void	b3StateRigidBodies::writeAllBodiesToGpu()
{
	m_bodyBufferGPU.resize(m_numRigidBodies);
	m_inertiaBufferGPU.resize(m_numRigidBodies);

	if (m_numRigidBodies)
	{
		m_bodyBufferGPU.copyFromHostPointer(&m_bodyBufferCPU[0], m_numRigidBodies);
		m_inertiaBufferGPU.copyFromHostPointer(&m_inertiaBufferCPU[0], m_numRigidBodies);
	}
}

void	b3StateRigidBodies::readbackAllBodiesToCpu()
{
	m_bodyBufferGPU.copyToHostPointer(&m_bodyBufferCPU[0], m_numRigidBodies);
}

void b3StateRigidBodies::setObjectTransformCpu(float* position, float* orientation, int bodyIndex)
{
	if (bodyIndex >= 0 && bodyIndex<m_bodyBufferCPU.size())
	{
		m_bodyBufferCPU[bodyIndex].m_pos = b3MakeVector3(position[0], position[1], position[2]);
		m_bodyBufferCPU[bodyIndex].m_quat.setValue(orientation[0], orientation[1], orientation[2], orientation[3]);
	}
	else
	{
		b3Warning("setObjectVelocityCpu out of range.\n");
	}
}
void b3StateRigidBodies::setObjectVelocityCpu(float* linVel, float* angVel, int bodyIndex)
{
	if (bodyIndex >= 0 && bodyIndex<m_bodyBufferCPU.size())
	{
		m_bodyBufferCPU[bodyIndex].m_linVel = b3MakeVector3(linVel[0], linVel[1], linVel[2]);
		m_bodyBufferCPU[bodyIndex].m_angVel = b3MakeVector3(angVel[0], angVel[1], angVel[2]);
	}
	else
	{
		b3Warning("setObjectVelocityCpu out of range.\n");
	}
}

bool b3StateRigidBodies::getObjectTransformFromCpu(float* position, float* orientation, int bodyIndex) const
{
	if (bodyIndex >= 0 && bodyIndex<m_bodyBufferCPU.size())
	{
		position[0] = m_bodyBufferCPU[bodyIndex].m_pos.x;
		position[1] = m_bodyBufferCPU[bodyIndex].m_pos.y;
		position[2] = m_bodyBufferCPU[bodyIndex].m_pos.z;
		position[3] = 1.f;

		orientation[0] = m_bodyBufferCPU[bodyIndex].m_quat.x;
		orientation[1] = m_bodyBufferCPU[bodyIndex].m_quat.y;
		orientation[2] = m_bodyBufferCPU[bodyIndex].m_quat.z;
		orientation[3] = m_bodyBufferCPU[bodyIndex].m_quat.w;
		return true;
	}

	b3Warning("getObjectTransformFromCpu out of range.\n");
	return false;
}

