/*
Copyright (c) 2013 Advanced Micro Devices, Inc.  

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
//Originally written by Erwin Coumans

#ifndef B3_GPU_RIGIDBODY_PIPELINE_H
#define B3_GPU_RIGIDBODY_PIPELINE_H

#include "Bullet3OpenCL/Initialize/b3OpenCLInclude.h"
#include "Bullet3Collision/NarrowPhaseCollision/b3Config.h"

#include "Bullet3Common/b3AlignedObjectArray.h"
#include "Bullet3Collision/NarrowPhaseCollision/b3RaycastInfo.h"


#include "Bullet3OpenCL/B3State/b3StateAabbs.h"
#include "Bullet3OpenCL/B3State/b3StateOverlappingPairs.h"
#include "Bullet3OpenCL/B3State/b3StateRigidCollidables.h"
#include "Bullet3OpenCL/B3State/b3StateRigidBodies.h"
#include "Bullet3OpenCL/B3State/b3StateRigidContacts.h"


class b3GpuRigidBodyPipeline
{
protected:
	struct b3GpuRigidBodyPipelineInternalData*	m_data;

public:
	b3StateAabbs m_aabbs;
	b3StateOverlappingPairs m_pairs;
	b3StateRigidCollidables m_collidables;
	b3StateRigidBodies m_rigidBodies;
	b3StateRigidContacts m_contacts;

	b3GpuRigidBodyPipeline(cl_context ctx,cl_device_id device, cl_command_queue  q , class b3GpuNarrowPhase* narrowphase, class b3GpuBroadphaseInterface* broadphase, const b3Config& config);
	virtual ~b3GpuRigidBodyPipeline();

	void	stepSimulation(float deltaTime);

	//void	castRays(const b3AlignedObjectArray<b3RayInfo>& rays,	b3AlignedObjectArray<b3RayHit>& hitResults);


};

#endif //B3_GPU_RIGIDBODY_PIPELINE_H