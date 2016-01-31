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

#include "b3GpuRigidBodyPipeline.h"
#include "b3GpuRigidBodyPipelineInternalData.h"
#include "kernels/integrateKernel.h"
#include "kernels/updateAabbsKernel.h"

#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"
#include "b3GpuNarrowPhase.h"
#include "Bullet3Geometry/b3AabbUtil.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3SapAabb.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3GpuBroadphaseInterface.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3Dynamics/ConstraintSolver/b3PgsJacobiSolver.h"
#include "Bullet3OpenCL/RigidBody/b3GpuPgsContactSolver.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3UpdateAabbs.h"
#include "Bullet3Collision/BroadPhaseCollision/b3DynamicBvhBroadphase.h"


#define B3_RIGIDBODY_INTEGRATE_PATH "src/Bullet3OpenCL/RigidBody/kernels/integrateKernel.cl"
#define B3_RIGIDBODY_UPDATEAABB_PATH "src/Bullet3OpenCL/RigidBody/kernels/updateAabbsKernel.cl"


#include "b3GpuJacobiContactSolver.h"

#include "Bullet3Collision/NarrowPhaseCollision/shared/b3RigidBodyData.h"
#include "Bullet3Collision/NarrowPhaseCollision/b3Contact4.h"

#include "Bullet3Collision/NarrowPhaseCollision/b3Config.h"

#include "Bullet3OpenCL/BroadphaseCollision/b3GpuParallelLinearBvhBroadphase.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhase.h"
#include "Bullet3OpenCL/RigidBody/b3GpuPgsConstraintSolver.h"
#include "Bullet3OpenCL/Raycast/b3GpuRaycast.h"

	
#include "Bullet3Dynamics/shared/b3IntegrateTransforms.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhaseInternalData.h"

#define MAX_RIGID_BODIES 32*1024			//todo: move 
#define MAX_CONTACTS 16*MAX_RIGID_BODIES		//todo: move

b3GpuRigidBodyPipeline::b3GpuRigidBodyPipeline(cl_context ctx, cl_device_id device, cl_command_queue q) :
m_aabbs(ctx, q),
m_pairs(ctx, q),
m_collidables(ctx, q),
m_rigidBodies(ctx, q, MAX_RIGID_BODIES),
m_constraints(ctx, q),
m_contacts(ctx, q, MAX_CONTACTS)
{
	m_data = new b3GpuRigidBodyPipelineInternalData;
	m_data->m_config = b3Config();
	m_data->m_context = ctx;
	m_data->m_device = device;
	m_data->m_queue = q;

	m_data->m_broadphase = new b3GpuParallelLinearBvhBroadphase(ctx, device, q);
	m_data->m_narrowphase = new b3GpuNarrowPhase(ctx, device, q, m_data->m_config);
	//m_data->m_jointSolver = new b3GpuPgsConstraintSolver(ctx, device, q, true);
	m_data->m_contactSolver = new b3GpuPgsContactSolver(ctx, device, q, m_data->m_config.m_maxBroadphasePairs);
	//m_data->m_raycaster = new b3GpuRaycast(ctx,device,q);

	m_data->m_gravity.setValue(0.f,-9.8f,0.f);

	cl_int errNum=0;

	{
		cl_program prog = b3OpenCLUtils::compileCLProgramFromString(m_data->m_context,m_data->m_device,integrateKernelCL,&errNum,"",B3_RIGIDBODY_INTEGRATE_PATH);
		b3Assert(errNum==CL_SUCCESS);
		m_data->m_integrateTransformsKernel = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,integrateKernelCL, "integrateTransformsKernel",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);
		clReleaseProgram(prog);
	}
	{
		cl_program prog = b3OpenCLUtils::compileCLProgramFromString(m_data->m_context,m_data->m_device,updateAabbsKernelCL,&errNum,"",B3_RIGIDBODY_UPDATEAABB_PATH);
		b3Assert(errNum==CL_SUCCESS);
		m_data->m_updateAabbsKernel = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,updateAabbsKernelCL, "initializeGpuAabbsFull",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);


		m_data->m_clearOverlappingPairsKernel = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,updateAabbsKernelCL, "clearOverlappingPairsKernel",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);

		clReleaseProgram(prog);
	}


}

b3GpuRigidBodyPipeline::~b3GpuRigidBodyPipeline()
{
	if (m_data->m_integrateTransformsKernel)
		clReleaseKernel(m_data->m_integrateTransformsKernel);
	
	if (m_data->m_updateAabbsKernel)
		clReleaseKernel(m_data->m_updateAabbsKernel);
	
	if (m_data->m_clearOverlappingPairsKernel)
		clReleaseKernel(m_data->m_clearOverlappingPairsKernel);
	
	//delete m_data->m_jointSolver;
	delete m_data->m_contactSolver;
	//delete m_data->m_raycaster;
	
	delete m_data;
}




void	b3GpuRigidBodyPipeline::stepSimulation(float deltaTime)
{
	int numBodies = m_rigidBodies.m_numRigidBodies;

	//update worldspace AABBs from local AABB/worldtransform
	{
		B3_PROFILE("setupGpuAabbs");

		if (!numBodies) return;

		{
			cl_mem bodies = m_rigidBodies.m_bodyBufferGPU.getBufferCL();
			cl_mem collidables = m_collidables.m_collidablesGpu.getBufferCL();
			cl_mem localAabbs = m_collidables.m_localShapeAabbGpu.getBufferCL();
			cl_mem worldAabbs = m_aabbs.m_aabbsGpu.getBufferCL();

			b3LauncherCL launcher(m_data->m_queue, m_data->m_updateAabbsKernel, "m_updateAabbsKernel");
			launcher.setConst(numBodies);
			launcher.setBuffer(bodies);
			launcher.setBuffer(collidables);
			launcher.setBuffer(localAabbs);
			launcher.setBuffer(worldAabbs);
			launcher.launch1D(numBodies);
		}
	}

	int numPairs = 0;

	//compute overlapping pairs
	{
		m_data->m_broadphase->computeOverlappingPairs(m_aabbs, m_pairs, m_data->m_config.m_maxBroadphasePairs);
		numPairs = m_pairs.m_overlappingPairsGpu.size();
	}

	//compute contact points
	int numContacts = 0;
	cl_mem contacts = 0;

	printf("numBodies:%d \n", numBodies);
	printf("numPairs:%d \n", numPairs);
	if (numPairs)
	{
		//mark the contacts for each pair as 'unused'
		if (numPairs)
		{
			cl_mem pairs = m_pairs.m_overlappingPairsGpu.getBufferCL();

			b3LauncherCL launcher(m_data->m_queue,m_data->m_clearOverlappingPairsKernel,"clearOverlappingPairsKernel");
			launcher.setBuffer(pairs);
			launcher.setConst(numPairs);
			launcher.launch1D(numPairs);
		}

		m_data->m_narrowphase->computeContacts(m_aabbs, m_pairs, m_rigidBodies, m_collidables, m_contacts);
		numContacts = m_contacts.m_pContactBuffersGPU[m_contacts.m_currentContactBuffer]->size();
		contacts = m_contacts.m_pContactBuffersGPU[m_contacts.m_currentContactBuffer]->getBufferCL();
	}
	printf("numContacts:%d \n", numContacts);

	if (numContacts)
	{
		int static0Index = 0;	//todo: check usage
		b3JacobiSolverInfo solverInfo;
		m_data->m_contactSolver->solveContacts(numBodies, m_rigidBodies.m_bodyBufferGPU.getBufferCL(), m_rigidBodies.m_inertiaBufferGPU.getBufferCL(),
												numContacts, contacts, m_data->m_config, static0Index);
	}
		
	{
		//integrate
		float angularDamp = 0.99f;

		{
			b3LauncherCL launcher(m_data->m_queue, m_data->m_integrateTransformsKernel, "m_integrateTransformsKernel");
			launcher.setBuffer( m_rigidBodies.m_bodyBufferGPU.getBufferCL() );

			launcher.setConst(numBodies);
			launcher.setConst(deltaTime);
			launcher.setConst(angularDamp);
			launcher.setConst(m_data->m_gravity);
			launcher.launch1D(numBodies);
		}
	}
}








#ifdef GPU_API_REDESIGN
void	b3GpuRigidBodyPipeline::castRays(const b3AlignedObjectArray<b3RayInfo>& rays,	b3AlignedObjectArray<b3RayHit>& hitResults)
{
	this->m_data->m_raycaster->castRays(rays,hitResults,
		getNumBodies(),this->m_data->m_narrowphase->getBodiesCpu(),
		m_data->m_narrowphase->getNumCollidablesGpu(), m_data->m_narrowphase->getCollidablesCpu(),
		m_data->m_narrowphase->getInternalData(), m_data->m_broadphaseSap);

}
#endif
