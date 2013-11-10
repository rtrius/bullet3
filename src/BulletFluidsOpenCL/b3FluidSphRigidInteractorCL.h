#ifndef B3_FLUID_SPH_RIGID_INTERACTOR_CL_H
#define B3_FLUID_SPH_RIGID_INTERACTOR_CL_H

#include "Bullet3Collision/NarrowPhaseCollision/b3RigidBodyCL.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"

#include "Bullet3OpenCL/NarrowphaseCollision/b3ConvexPolyhedronCL.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3SapAabb.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3GpuSapBroadphase.h"

#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"

#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhase.h"

#include "b3FluidSphOpenCL.h"
#include "b3FluidSortingGridOpenCL.h"
#include "b3FluidHashGridOpenCL.h"

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
	
	//Concave BVH data
	int m_numBvhInfo;
	int m_numBvhSubtreeInfo;
	int m_numBvhNodes;
	cl_mem m_bvhInfo;				//b3BvhInfo
	cl_mem m_bvhSubtreeInfo;		//b3BvhSubtreeInfo
	cl_mem m_bvhNodes;				//b3QuantizedBvhNode
	
	//Compound shape data
	int m_numGpuChildShapes;
	cl_mem m_gpuChildShapes; 		//b3GpuChildShape
	
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
		
		m_numBvhInfo = narrowPhase->getNumBvhInfoGpu();
		m_numBvhSubtreeInfo = narrowPhase->getNumBvhSubtreeInfoGpu();
		m_numBvhNodes = narrowPhase->getNumBvhNodesGpu();
		m_bvhInfo = narrowPhase->getBvhInfoGpu();
		m_bvhSubtreeInfo = narrowPhase->getBvhSubtreeInfoGpu();
		m_bvhNodes = narrowPhase->getBvhNodesGpu();
		
		m_numGpuChildShapes = narrowPhase->getNumGpuChildShapesGpu();
		m_gpuChildShapes = narrowPhase->getGpuChildShapesGpu();
	
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
#define MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE 16
#define MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID 256

///Contains the indicies of rigid bodies whose AABB intersects with that of a single fluid particle
struct FluidRigidPairs
{
	int m_numIndicies;
	int m_rigidIndicies[MAX_FLUID_RIGID_PAIRS];
	int m_rigidSubIndicies[MAX_FLUID_RIGID_PAIRS];	//Only used if the rigid has triangle mesh or compound shape
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
	///Maximum number of rigid bodies with large AABBs(AABBs larger than the fluid's grid) that are considered for collision
	static const int MAX_LARGE_AABB_RIGIDS = 32;	
	
	cl_context m_context;
	cl_command_queue m_commandQueue;
	
	b3OpenCLArray<FluidRigidPairs> m_pairs;
	b3OpenCLArray<FluidRigidPairs> m_midphasePairs;
	b3OpenCLArray<FluidRigidContacts> m_fluidRigidContacts;
	b3OpenCLArray<RigidFluidContacts> m_rigidFluidContacts;
	b3OpenCLArray<b3Vector3> m_fluidVelocities;
	
	//Large AABB is defined here as exceeding the extent of the fluid's grid;
	//it is not related to the b3GpuSapBroadphase
	b3OpenCLArray<int> m_numLargeAabbRigid;
	b3OpenCLArray<int> m_largeAabbRigidIndicies;
	
	cl_program m_fluidRigidProgram;
	
	cl_kernel m_clearFluidRigidPairsAndContactsKernel;
	cl_kernel m_detectLargeAabbRigidsKernel;
	cl_kernel m_fluidLargeRigidBroadphaseKernel;
	cl_kernel m_fluidSmallRigidBroadphaseKernel;
	cl_kernel m_fluidSmallRigidBroadphaseModuloKernel;
	cl_kernel m_fluidRigidMidphaseKernel;
	cl_kernel m_fluidRigidNarrowphaseKernel;
	cl_kernel m_resolveFluidRigidCollisionsKernel;
	
	cl_kernel m_clearRigidFluidContactsKernel;
	cl_kernel m_mapRigidFluidContactsKernel;
	cl_kernel m_resolveRigidFluidCollisionsKernel;
	
public:
	b3FluidSphRigidInteractorCL(cl_context context, cl_device_id device, cl_command_queue queue) 
	:	m_pairs(context, queue), m_midphasePairs(context, queue), m_fluidRigidContacts(context, queue), 
		m_rigidFluidContacts(context, queue), m_fluidVelocities(context, queue),
		m_numLargeAabbRigid(context, queue), m_largeAabbRigidIndicies(context, queue)
	{
		m_context = context;
		m_commandQueue = queue;
	
		//
		m_numLargeAabbRigid.resize(1);
		m_largeAabbRigidIndicies.resize(MAX_LARGE_AABB_RIGIDS);
	
		//
		const char CL_PROGRAM_PATH[] = "src/BulletFluidsOpenCL/fluidSphRigid.cl";
		
		const char* kernelSource = fluidSphRigidCL;	//fluidSphRigidCL.h
		cl_int error;
		char* additionalMacros = 0;
		m_fluidRigidProgram = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, additionalMacros, CL_PROGRAM_PATH);
		b3Assert(m_fluidRigidProgram);
		
		m_clearFluidRigidPairsAndContactsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "clearFluidRigidPairsAndContacts", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_clearFluidRigidPairsAndContactsKernel);
		m_detectLargeAabbRigidsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "detectLargeAabbRigids", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_detectLargeAabbRigidsKernel);
		m_fluidLargeRigidBroadphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidLargeRigidBroadphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidLargeRigidBroadphaseKernel);
		m_fluidSmallRigidBroadphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidSmallRigidBroadphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidSmallRigidBroadphaseKernel);
		m_fluidSmallRigidBroadphaseModuloKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidSmallRigidBroadphaseModulo", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidSmallRigidBroadphaseModuloKernel);
		m_fluidRigidMidphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidRigidMidphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidRigidMidphaseKernel);
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
		clReleaseKernel(m_detectLargeAabbRigidsKernel);
		clReleaseKernel(m_fluidLargeRigidBroadphaseKernel);
		clReleaseKernel(m_fluidSmallRigidBroadphaseKernel);
		clReleaseKernel(m_fluidSmallRigidBroadphaseModuloKernel);
		clReleaseKernel(m_fluidRigidMidphaseKernel);
		clReleaseKernel(m_fluidRigidNarrowphaseKernel);
		clReleaseKernel(m_resolveFluidRigidCollisionsKernel);
		
		clReleaseKernel(m_clearRigidFluidContactsKernel);
		clReleaseKernel(m_mapRigidFluidContactsKernel);
		clReleaseKernel(m_resolveRigidFluidCollisionsKernel);
		
		clReleaseProgram(m_fluidRigidProgram);
	}
	
	//Either gridData or moduloGridData must be nonzero, but not both; it is used to determine which grid type is being used
	void interact(const b3OpenCLArray<b3FluidSphParametersGlobal>& globalFluidParams, b3FluidSphOpenCL* fluidData, 
					b3FluidSortingGridOpenCL* gridData, b3FluidHashGridOpenCL* moduloGridData, RigidBodyGpuData& rigidBodyData)
	{
		b3Assert(gridData || moduloGridData);
		b3Assert( !(gridData && moduloGridData));
	
		clFinish(m_commandQueue);
	
		B3_PROFILE("b3FluidSphRigidInteractorCL::interact()");
	
		int numRigidBodies = rigidBodyData.m_numRigidBodies;
		int numFluidParticles = fluidData->m_pos.size();
		int numGridCells = (gridData) ? gridData->getNumActiveCells() : B3_FLUID_HASH_GRID_NUM_CELLS;
		
		if(!numRigidBodies || !numFluidParticles || !numGridCells) return;
		
		//
		if(m_pairs.size() < numFluidParticles) m_pairs.resize(numFluidParticles);
		if(m_midphasePairs.size() < numFluidParticles) m_midphasePairs.resize(numFluidParticles);
		if(m_fluidRigidContacts.size() < numFluidParticles) m_fluidRigidContacts.resize(numFluidParticles);
		if(m_rigidFluidContacts.size() < numRigidBodies) m_rigidFluidContacts.resize(numRigidBodies);
		
		//Clear broadphase pairs and contacts
		{
			B3_PROFILE("m_clearFluidRigidPairsAndContactsKernel");
		
			b3BufferInfoCL bufferInfo[] = 
			{ 
				b3BufferInfoCL( m_pairs.getBufferCL() ),
				b3BufferInfoCL( m_midphasePairs.getBufferCL() ),
				b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_clearFluidRigidPairsAndContactsKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
			clFinish(m_commandQueue);
		}
	
		//Broadphase - large AABB rigids
		{
			
			const int reset = 0;
			m_numLargeAabbRigid.copyFromHostPointer(&reset, 1);
			
			//Detect rigids with large AABBs
			{
				B3_PROFILE("m_detectLargeAabbRigidsKernel");
				
				b3BufferInfoCL bufferInfo[] = 
				{
					b3BufferInfoCL( globalFluidParams.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
					
					b3BufferInfoCL( rigidBodyData.m_worldSpaceAabbs ),
					b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
					b3BufferInfoCL( rigidBodyData.m_collidables ),
					
					b3BufferInfoCL( m_numLargeAabbRigid.getBufferCL() ),
					b3BufferInfoCL( m_largeAabbRigidIndicies.getBufferCL() ),
				};
				
				b3LauncherCL launcher(m_commandQueue, m_detectLargeAabbRigidsKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(MAX_LARGE_AABB_RIGIDS);
				launcher.setConst(numRigidBodies);
				
				launcher.launch1D(numRigidBodies);
				clFinish(m_commandQueue);
			}
			
			//Intersection test with fluid AABBs
			{
				B3_PROFILE("m_fluidLargeRigidBroadphaseKernel");
				
				b3BufferInfoCL bufferInfo[] = 
				{
					b3BufferInfoCL( globalFluidParams.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
				
					b3BufferInfoCL( rigidBodyData.m_worldSpaceAabbs ),
					b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
					b3BufferInfoCL( rigidBodyData.m_collidables ),
					
					b3BufferInfoCL( m_numLargeAabbRigid.getBufferCL() ),
					b3BufferInfoCL( m_largeAabbRigidIndicies.getBufferCL() ),
					
					b3BufferInfoCL( m_pairs.getBufferCL() ),
					b3BufferInfoCL( m_midphasePairs.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_fluidLargeRigidBroadphaseKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(MAX_LARGE_AABB_RIGIDS);
				launcher.setConst(numFluidParticles);
				
				launcher.launch1D(numFluidParticles);
				clFinish(m_commandQueue);
			}
		}
		
		//Broadphase - small AABB rigids
		if(gridData)
		{
			B3_PROFILE("m_fluidSmallRigidBroadphaseKernel");
				
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
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_collidables ),
				
				b3BufferInfoCL( m_pairs.getBufferCL() ),
				b3BufferInfoCL( m_midphasePairs.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_fluidSmallRigidBroadphaseKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numGridCells);
			launcher.setConst(numRigidBodies);
			
			launcher.launch1D(numRigidBodies);
			clFinish(m_commandQueue);
		}
		else if(moduloGridData)
		{
			B3_PROFILE("m_fluidSmallRigidBroadphaseModuloKernel");
			
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( globalFluidParams.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
				
				b3BufferInfoCL( moduloGridData->m_cellContents.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_worldSpaceAabbs ),
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_collidables ),
				
				b3BufferInfoCL( m_pairs.getBufferCL() ),
				b3BufferInfoCL( m_midphasePairs.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_fluidSmallRigidBroadphaseModuloKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numRigidBodies);
			
			launcher.launch1D(numRigidBodies);
			clFinish(m_commandQueue);
		}
		else
		{
			b3Assert(0);	//No grid data
		}
	
		//Midphase - for triangle mesh and compound shapes
		{
			B3_PROFILE("m_fluidRigidMidphaseKernel");
		
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_pos.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_collidables ),
				
				b3BufferInfoCL( rigidBodyData.m_bvhInfo ),
				b3BufferInfoCL( rigidBodyData.m_bvhSubtreeInfo ),
				b3BufferInfoCL( rigidBodyData.m_bvhNodes ),
				
				b3BufferInfoCL( m_midphasePairs.getBufferCL() ),
				b3BufferInfoCL( m_pairs.getBufferCL() ),
			};
			
			b3LauncherCL launcher(m_commandQueue, m_fluidRigidMidphaseKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
			clFinish(m_commandQueue);
		}
		
		//Narrowphase
		{
			B3_PROFILE("m_fluidRigidNarrowphaseKernel");
			
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
				b3BufferInfoCL( rigidBodyData.m_gpuChildShapes ),
				
				b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_fluidRigidNarrowphaseKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
			clFinish(m_commandQueue);
		}
		
		//m_resolveFluidRigidCollisionsKernel, executed below, overwrites fluid particle velocities(m_vel);
		//the pre-collision particle velocities are needed to calculate impulses for the dynamic rigid bodies
		m_fluidVelocities.copyFromOpenCLArray(fluidData->m_vel);
		
		//Resolve Collisions - apply impulses to fluid particles
		{
			B3_PROFILE("m_resolveFluidRigidCollisionsKernel");
		
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
			clFinish(m_commandQueue);
		}
		
		
		//Map fluid contacts to rigid bodies
		//Since applying impulses simultaneously(in 1 kernel) to both fluid and rigid body would require
		//syncronization between threads, we instead run 2 kernels - the first iterates through all contacts with
		//1 thread per fluid particle, and the second does the same with 1 thread per rigid body.
		{
			//Clear rigid side contacts
			{
				B3_PROFILE("m_clearRigidFluidContactsKernel");
			
				b3BufferInfoCL bufferInfo[] = 
				{ 
					b3BufferInfoCL( m_rigidFluidContacts.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_clearRigidFluidContactsKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numRigidBodies);
				
				launcher.launch1D(numRigidBodies);
				clFinish(m_commandQueue);
			}
			
			//Map fluid to rigid
			{
				B3_PROFILE("m_mapRigidFluidContactsKernel");
				
				b3BufferInfoCL bufferInfo[] = 
				{ 
					b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() ),
					b3BufferInfoCL( m_rigidFluidContacts.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_mapRigidFluidContactsKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidParticles);
				
				launcher.launch1D(numFluidParticles);
				clFinish(m_commandQueue);
			}
		}
		
		//Resolve Collisions - apply impulses to rigid bodies
		{
			B3_PROFILE("m_resolveRigidFluidCollisionsKernel");
				
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
			clFinish(m_commandQueue);
		}
		
		clFinish(m_commandQueue);
	}
	
};

#endif
