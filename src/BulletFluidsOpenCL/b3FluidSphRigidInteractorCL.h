#ifndef B3_FLUID_SPH_RIGID_INTERACTOR_CL_H
#define B3_FLUID_SPH_RIGID_INTERACTOR_CL_H

#include "Bullet3Collision/NarrowPhaseCollision/b3RigidBodyCL.h"

#include "Bullet3OpenCL/NarrowPhaseCollision/b3Collidable.h"
#include "Bullet3OpenCL/NarrowphaseCollision/b3ConvexPolyhedronCL.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3SapAabb.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3GpuSapBroadphase.h"

#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"

#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhase.h"

#include "b3FluidSphOpenCL.h"
#include "b3FluidSortingGridOpenCL.h"

#include "fluidSphRigidCL.h"

//Data for all GPU rigid bodies
struct RigidBodyGpuData
{
	//Rigid bodies
	int m_numRigidBodies;
	int m_numRigidBodyInertias;
	int m_numWorldSpaceAabbs;
	cl_mem m_rigidBodies;			//b3RigidBody
	cl_mem m_rigidBodyInertias;		//b3InertiaCL
	cl_mem m_worldSpaceAabbs;		//b3SapAabb
	
	//Shapes
	int m_numCollidables;
	cl_mem m_collidables;			//b3Collidable
	
	//Shape data
	int m_numConvexPolyhedra;
	int m_numFaces;
	int m_numConvexIndicies;
	int m_numConvexVertices;
	cl_mem m_convexPolyhedra;		//b3ConvexPolyhedronCL
	cl_mem m_faces;					//btGpuFace
	cl_mem m_convexIndices; 		//int
	cl_mem m_convexVertices;		//b3Vector3
	
	void load(b3GpuSapBroadphase* broadphase, b3GpuNarrowPhase* narrowPhase)
	{
		m_numRigidBodies = narrowPhase->getNumBodiesGpu();
		m_numRigidBodyInertias = narrowPhase->getNumBodyInertiasGpu();
		m_numWorldSpaceAabbs = broadphase->getNumAabbWS();
		m_rigidBodies = narrowPhase->getBodiesGpu();
		m_rigidBodyInertias = narrowPhase->getBodyInertiasGpu();
		m_worldSpaceAabbs = broadphase->getAabbBufferWS();		//if useDbvt == true in b3GpuRigidBodyPipeline.cpp this is incorrect
		b3Assert(m_numRigidBodies == m_numRigidBodyInertias);
		b3Assert(m_numRigidBodies == m_numWorldSpaceAabbs);
		
		m_numCollidables = narrowPhase->getNumCollidablesGpu();
		m_collidables = narrowPhase->getCollidablesGpu();
		
		m_numConvexPolyhedra = narrowPhase->getNumConvexPolyhedraGpu();
		m_numFaces = narrowPhase->getNumFacesGpu();
		m_numConvexIndicies = narrowPhase->getNumConvexIndiciesGpu();
		m_numConvexVertices = narrowPhase->getNumConvexVerticesGpu();
		m_convexPolyhedra = narrowPhase->getConvexPolyhedraGpu();
		m_faces = narrowPhase->getFacesGpu();
		m_convexIndices = narrowPhase->getConvexIndiciesGpu();
		m_convexVertices = narrowPhase->getConvexVerticesGpu();
		
		if(0)
		{
			printf("RigidBodyGpuData::load()\n");
			printf("m_numRigidBodies: %d \n", m_numRigidBodies);
			printf("m_numRigidBodyInertias: %d \n", m_numRigidBodyInertias);
			printf("m_numWorldSpaceAabbs: %d \n", m_numWorldSpaceAabbs);
			printf("m_numCollidables: %d \n", m_numCollidables);
			printf("m_numConvexPolyhedra: %d \n", m_numConvexPolyhedra);
			printf("m_numFaces: %d \n", m_numFaces);
			printf("m_numConvexIndicies: %d \n", m_numConvexIndicies);
			printf("m_numConvexVertices: %d \n", m_numConvexVertices);
			printf("\n");
		}
	}
};

#define MAX_FLUID_RIGID_PAIRS 32
#define MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE 4
#define MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID 256

///Contains the indicies of rigid bodies whose AABB intersects with that of a single fluid particle
struct FluidRigidPairs
{
	int m_numIndicies;
	int m_rigidIndicies[MAX_FLUID_RIGID_PAIRS];
};

struct FluidRigidContacts
{
	int m_numContacts;
	int m_rigidIndicies[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];
	b3Scalar m_distances[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];
	b3Vector3 m_pointsOnRigid[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];
	b3Vector3 m_normalsOnRigid[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];
};

struct RigidFluidContacts
{
	int m_numContacts;
	int m_fluidIndicies[MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID];
	int m_contactIndicies[MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID];
};

///Handles collision detection and response between SPH particles and (GPU) rigid bodies
class b3FluidSphRigidInteractorCL
{
	cl_context m_context;
	cl_command_queue m_commandQueue;
	
	b3OpenCLArray<FluidRigidPairs> m_pairs;
	b3OpenCLArray<FluidRigidContacts> m_fluidRigidContacts;
	b3OpenCLArray<RigidFluidContacts> m_rigidFluidContacts;
	b3OpenCLArray<b3Vector3> m_fluidVelocities;
	
	cl_program m_fluidRigidProgram;
	
	cl_kernel m_clearFluidRigidPairsAndContactsKernel;
	cl_kernel m_fluidRigidBroadphaseKernel;
	cl_kernel m_fluidRigidNarrowphaseKernel;
	cl_kernel m_resolveFluidRigidCollisionsKernel;
	
	cl_kernel m_clearRigidFluidContactsKernel;
	cl_kernel m_mapRigidFluidContactsKernel;
	cl_kernel m_resolveRigidFluidCollisionsKernel;
	
public:
	b3FluidSphRigidInteractorCL(cl_context context, cl_device_id device, cl_command_queue queue) 
	: m_pairs(context, queue), m_fluidRigidContacts(context, queue), m_rigidFluidContacts(context, queue), m_fluidVelocities(context, queue)
	{
		m_context = context;
		m_commandQueue = queue;
	
		const char CL_PROGRAM_PATH[] = "src/BulletFluidsOpenCL/fluidSphRigid.cl";
		
		const char* kernelSource = fluidSphRigidCL;	//fluidSphRigidCL.h
		cl_int error;
		char* additionalMacros = 0;
		m_fluidRigidProgram = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, additionalMacros, CL_PROGRAM_PATH);
		b3Assert(m_fluidRigidProgram);
		
		m_clearFluidRigidPairsAndContactsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "clearFluidRigidPairsAndContacts", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_clearFluidRigidPairsAndContactsKernel);
		m_fluidRigidBroadphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidRigidBroadphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidRigidBroadphaseKernel);
		m_fluidRigidNarrowphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidRigidNarrowphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidRigidNarrowphaseKernel);
		m_resolveFluidRigidCollisionsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "resolveFluidRigidCollisions", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_resolveFluidRigidCollisionsKernel);
		
		m_clearRigidFluidContactsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "clearRigidFluidContacts", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_clearRigidFluidContactsKernel);
		m_mapRigidFluidContactsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "mapRigidFluidContacts", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_mapRigidFluidContactsKernel);
		m_resolveRigidFluidCollisionsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "resolveRigidFluidCollisions", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_resolveRigidFluidCollisionsKernel);
	}
	
	virtual ~b3FluidSphRigidInteractorCL()
	{
		clReleaseKernel(m_clearFluidRigidPairsAndContactsKernel);
		clReleaseKernel(m_fluidRigidBroadphaseKernel);
		clReleaseKernel(m_fluidRigidNarrowphaseKernel);
		clReleaseKernel(m_resolveFluidRigidCollisionsKernel);
		
		clReleaseKernel(m_clearRigidFluidContactsKernel);
		clReleaseKernel(m_mapRigidFluidContactsKernel);
		clReleaseKernel(m_resolveRigidFluidCollisionsKernel);
		
		clReleaseProgram(m_fluidRigidProgram);
	}
	

	void interact(const b3OpenCLArray<b3FluidSphParametersGlobal>& globalFluidParams, b3FluidSphOpenCL* fluidData, 
					b3FluidSortingGridOpenCL* gridData, RigidBodyGpuData& rigidBodyData)
	{
		clFinish(m_commandQueue);
	
		B3_PROFILE("b3FluidSphRigidInteractorCL::interact()");
	
		int numRigidBodies = rigidBodyData.m_numRigidBodies;
		int numFluidParticles = fluidData->m_pos.size();
		int numGridCells = gridData->getNumActiveCells();
		
		if(!numRigidBodies || !numFluidParticles || !numGridCells) return;
		
		
		if(m_pairs.size() < numFluidParticles) m_pairs.resize(numFluidParticles);
		if(m_fluidRigidContacts.size() < numFluidParticles) m_fluidRigidContacts.resize(numFluidParticles);
		if(m_rigidFluidContacts.size() < numRigidBodies) m_rigidFluidContacts.resize(numRigidBodies);
		
		//Clear broadphase pairs and contacts
		{
			b3BufferInfoCL bufferInfo[] = 
			{ 
				b3BufferInfoCL( m_pairs.getBufferCL() ),
				b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_clearFluidRigidPairsAndContactsKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
	
		//Broadphase
		{
			//Quantize the rigid AABB into fluid grid coordinates,
			//then test each particle in each intersecting grid cell;
			//if the rigid AABB intersects with the particle AABB, add it
			//to a per-particle array of broadphase pairs.
			
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( globalFluidParams.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
				
				b3BufferInfoCL( gridData->m_activeCells.getBufferCL() ),
				b3BufferInfoCL( gridData->m_cellContents.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_worldSpaceAabbs ),
				
				b3BufferInfoCL( m_pairs.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_fluidRigidBroadphaseKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numGridCells);
			launcher.setConst(numRigidBodies);
			
			launcher.launch1D(numRigidBodies);
		}
		
		//Narrowphase
		{
			//Use the Separating Axis Test for each fluid-rigid pair
			
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( globalFluidParams.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
				b3BufferInfoCL( m_pairs.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_collidables ),
				b3BufferInfoCL( rigidBodyData.m_convexPolyhedra ),
				b3BufferInfoCL( rigidBodyData.m_faces ),
				b3BufferInfoCL( rigidBodyData.m_convexIndices ),
				b3BufferInfoCL( rigidBodyData.m_convexVertices ),
				
				b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_fluidRigidNarrowphaseKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
		
		//m_resolveFluidRigidCollisionsKernel, executed below, overwrites fluid particle velocities(m_vel);
		//the pre-collision particle velocities are needed to calculate impulses for the dynamic rigid bodies
		m_fluidVelocities.copyFromOpenCLArray(fluidData->m_vel);
		
		//Resolve Collisions - apply impulses to fluid particles
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( globalFluidParams.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_rigidBodyInertias ),
				b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() ),
				
				b3BufferInfoCL( fluidData->m_vel.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_vel_eval.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_resolveFluidRigidCollisionsKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
		}
		
		
		//Map fluid contacts to rigid bodies
		{
			//Clear rigid side contacts
			{
				b3BufferInfoCL bufferInfo[] = 
				{ 
					b3BufferInfoCL( m_rigidFluidContacts.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_clearRigidFluidContactsKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numRigidBodies);
				
				launcher.launch1D(numRigidBodies);
			}
			
			//Map fluid to rigid
			{
				b3BufferInfoCL bufferInfo[] = 
				{ 
					b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() ),
					b3BufferInfoCL( m_rigidFluidContacts.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_mapRigidFluidContactsKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidParticles);
				
				launcher.launch1D(numFluidParticles);
			}
		}
		
		//Resolve Collisions - apply impulses to rigid bodies
		{
			b3BufferInfoCL bufferInfo[] = 
			{ 
				b3BufferInfoCL( globalFluidParams.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_rigidBodyInertias ),
				b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() ),
				b3BufferInfoCL( m_rigidFluidContacts.getBufferCL() ),
				
				b3BufferInfoCL( m_fluidVelocities.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_resolveRigidFluidCollisionsKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numRigidBodies);
			
			launcher.launch1D(numRigidBodies);
		}
		
		clFinish(m_commandQueue);
	}
	
};

#endif