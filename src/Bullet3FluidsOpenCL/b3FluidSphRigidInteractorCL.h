#ifndef B3_FLUID_SPH_RIGID_INTERACTOR_CL_H
#define B3_FLUID_SPH_RIGID_INTERACTOR_CL_H

#include "Bullet3Collision/NarrowPhaseCollision/b3RigidBodyCL.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"

#include "Bullet3OpenCL/NarrowphaseCollision/b3ConvexPolyhedronCL.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3SapAabb.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3GpuSapBroadphase.h"

#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3FillCL.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3RadixSort32CL.h"

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
			printf( "narrowPhase->getNumRigidBodies(): %d \n", narrowPhase->getNumRigidBodies() );
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
		b3Assert(m_numRigidBodies == m_numRigidBodyInertias);
		b3Assert(m_numRigidBodies == m_numWorldSpaceAabbs);
	}
};



typedef struct
{
	int m_fluidParticleIndex;
	int m_rigidIndex;
	int m_rigidSubIndex;		//Only used if the rigid has triangle mesh or compound shape
	int m_rigidShapeType;		//For sorting
} FluidRigidPair;

typedef struct
{
	b3Scalar m_distance;
	b3Vector3 m_pointOnRigid;
	b3Vector3 m_normalOnRigid;
} FluidRigidContact;

///Handles collision detection and response between SPH particles and (GPU) rigid bodies
class b3FluidSphRigidInteractorCL
{
	//If either of these values is exceeded, the interactor silently discards additional contacts/broadphase pairs
	static const int MAX_FLUID_RIGID_PAIRS = 131072*4; 	//Also the max number of fluid-rigid contacts
	static const int MAX_MIDPHASE_PAIRS = 131072*1;

	///Maximum number of rigid bodies with large AABBs(AABBs larger than the fluid's grid) that are considered for collision
	static const int MAX_LARGE_AABB_RIGIDS = 32;	
	
	cl_context m_context;
	cl_command_queue m_commandQueue;
	
	b3FillCL m_fill;
	b3RadixSort32CL m_radixSorter;
	
	//1 element per pair of fluid particle and (convex hull, triangle, or sphere shape) rigid with intersecting AABB
	b3OpenCLArray<int> m_numPairs;
	b3OpenCLArray<FluidRigidPair> m_pairs;
	b3OpenCLArray<FluidRigidContact> m_fluidRigidContacts;
	b3OpenCLArray<FluidRigidPair> m_pairsTemp;
	b3OpenCLArray<FluidRigidContact> m_fluidRigidContactsTemp;
	b3OpenCLArray<b3SortData> m_sortByFluidIndexData;
	b3OpenCLArray<b3SortData> m_fluidToRigidMap;		//m_key == rigid body index, m_value == contact index
	
	//1 element per rigid body
	b3OpenCLArray<int> m_firstContactIndexPerRigid;	//Contains indices of m_fluidToRigidMap, which is where the per rigid contact indices are 
	b3OpenCLArray<int> m_numContactsPerRigid;
	
	//1 element per pair of fluid particle and (trimesh or compound shape) rigid with intersecting AABB
	b3OpenCLArray<int> m_numMidphasePairs;
	b3OpenCLArray<FluidRigidPair> m_midphasePairs;
	
	//1 element per fluid particle
	b3OpenCLArray<int> m_firstContactIndexPerParticle;
	b3OpenCLArray<int> m_numContactsPerParticle;
	b3OpenCLArray<b3Vector3> m_fluidVelocities;
	
	//Large AABB is defined here as exceeding the extent of the fluid's grid;
	//it is not related to the b3GpuSapBroadphase
	b3OpenCLArray<int> m_numLargeAabbRigid;
	b3OpenCLArray<int> m_largeAabbRigidIndicies;
	
	cl_program m_fluidRigidProgram;
	
	cl_kernel m_detectLargeAabbRigidsKernel;
	cl_kernel m_fluidLargeRigidBroadphaseKernel;
	cl_kernel m_fluidSmallRigidBroadphaseKernel;
	//cl_kernel m_fluidSmallRigidBroadphaseModuloKernel;
	cl_kernel m_fluidRigidMidphaseKernel;
	cl_kernel m_fluidRigidNarrowphaseKernel;
	cl_kernel m_loadSortDataKernel;
	cl_kernel m_rearrangePairsAndContactsKernel;
	cl_kernel m_findPerParticleContactRangeKernel;
	cl_kernel m_resolveFluidRigidCollisionsKernel;
	
	cl_kernel m_mapFluidToRigidContactsKernel;
	cl_kernel m_findPerRigidContactRangeKernel;
	cl_kernel m_resolveRigidFluidCollisionsKernel;
	
public:
	b3FluidSphRigidInteractorCL(cl_context context, cl_device_id device, cl_command_queue queue) 
	:	m_fill(context, device, queue), 
		m_radixSorter(context, device, queue),
		
		m_numPairs(context, queue), 
		m_pairs(context, queue), 
		m_fluidRigidContacts(context, queue),
		m_pairsTemp(context, queue), 
		m_fluidRigidContactsTemp(context, queue), 
		m_sortByFluidIndexData(context, queue), 
		m_fluidToRigidMap(context, queue), 
		
		m_firstContactIndexPerRigid(context, queue), 
		m_numContactsPerRigid(context, queue),
		
		m_numMidphasePairs(context, queue), 
		m_midphasePairs(context, queue),
		
		m_fluidVelocities(context, queue), 
		m_firstContactIndexPerParticle(context, queue), 
		m_numContactsPerParticle(context, queue),
		
		m_numLargeAabbRigid(context, queue), 
		m_largeAabbRigidIndicies(context, queue)
	{
		m_context = context;
		m_commandQueue = queue;
	
		//
		{
			m_numPairs.resize(1);
			m_pairs.resize(MAX_FLUID_RIGID_PAIRS);
			m_fluidRigidContacts.resize(MAX_FLUID_RIGID_PAIRS);
			m_pairsTemp.resize(MAX_FLUID_RIGID_PAIRS);
			m_fluidRigidContactsTemp.resize(MAX_FLUID_RIGID_PAIRS);
			m_sortByFluidIndexData.resize(MAX_FLUID_RIGID_PAIRS);
			m_fluidToRigidMap.resize(MAX_FLUID_RIGID_PAIRS);
			
			m_numMidphasePairs.resize(1);
			m_midphasePairs.resize(MAX_MIDPHASE_PAIRS);
			
			m_numLargeAabbRigid.resize(1);
			m_largeAabbRigidIndicies.resize(MAX_LARGE_AABB_RIGIDS);
		}
		
		//
		const char CL_PROGRAM_PATH[] = "src/Bullet3FluidsOpenCL/fluidSphRigid.cl";
		
		const char* kernelSource = fluidSphRigidCL;	//fluidSphRigidCL.h
		cl_int error;
		char* additionalMacros = 0;
		m_fluidRigidProgram = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, additionalMacros, CL_PROGRAM_PATH);
		b3Assert(m_fluidRigidProgram);
		
		m_detectLargeAabbRigidsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "detectLargeAabbRigids", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_detectLargeAabbRigidsKernel);
		m_fluidLargeRigidBroadphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidLargeRigidBroadphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidLargeRigidBroadphaseKernel);
		m_fluidSmallRigidBroadphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidSmallRigidBroadphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidSmallRigidBroadphaseKernel);
	//	m_fluidSmallRigidBroadphaseModuloKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidSmallRigidBroadphaseModulo", &error, m_fluidRigidProgram, additionalMacros );
	//	b3Assert(m_fluidSmallRigidBroadphaseModuloKernel);
		m_fluidRigidMidphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidRigidMidphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidRigidMidphaseKernel);
		m_fluidRigidNarrowphaseKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "fluidRigidNarrowphase", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_fluidRigidNarrowphaseKernel);
		m_loadSortDataKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "loadSortData", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_loadSortDataKernel);
		m_rearrangePairsAndContactsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "rearrangePairsAndContacts", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_rearrangePairsAndContactsKernel);
		m_findPerParticleContactRangeKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "findPerParticleContactRange", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_findPerParticleContactRangeKernel);
		m_resolveFluidRigidCollisionsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "resolveFluidRigidCollisions", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_resolveFluidRigidCollisionsKernel);
		
		m_mapFluidToRigidContactsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "mapFluidToRigidContacts", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_mapFluidToRigidContactsKernel);
		m_findPerRigidContactRangeKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "findPerRigidContactRange", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_findPerRigidContactRangeKernel);
		m_resolveRigidFluidCollisionsKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "resolveRigidFluidCollisions", &error, m_fluidRigidProgram, additionalMacros );
		b3Assert(m_resolveRigidFluidCollisionsKernel);
	}
	
	virtual ~b3FluidSphRigidInteractorCL()
	{
		clReleaseKernel(m_detectLargeAabbRigidsKernel);
		clReleaseKernel(m_fluidLargeRigidBroadphaseKernel);
		clReleaseKernel(m_fluidSmallRigidBroadphaseKernel);
		//clReleaseKernel(m_fluidSmallRigidBroadphaseModuloKernel);
		clReleaseKernel(m_fluidRigidMidphaseKernel);
		clReleaseKernel(m_fluidRigidNarrowphaseKernel);
		clReleaseKernel(m_loadSortDataKernel);
		clReleaseKernel(m_rearrangePairsAndContactsKernel);
		clReleaseKernel(m_findPerParticleContactRangeKernel);
		clReleaseKernel(m_resolveFluidRigidCollisionsKernel);
		
		clReleaseKernel(m_mapFluidToRigidContactsKernel);
		clReleaseKernel(m_findPerRigidContactRangeKernel);
		clReleaseKernel(m_resolveRigidFluidCollisionsKernel);
		
		clReleaseProgram(m_fluidRigidProgram);
	}
	
	//Either gridData or moduloGridData must be nonzero, but not both; it is used to determine which grid type is being used
	void interact(b3FluidSphOpenCL* fluidData,  b3FluidSortingGridOpenCL* gridData, 
					b3FluidHashGridOpenCL* moduloGridData, RigidBodyGpuData& rigidBodyData)
	{
		b3Assert(gridData || moduloGridData);
		b3Assert( !(gridData && moduloGridData));
	
		clFinish(m_commandQueue);
	
		B3_PROFILE("b3FluidSphRigidInteractorCL::interact()");
	
		int numRigidBodies = rigidBodyData.m_numRigidBodies;
		int numFluidParticles = fluidData->m_position.size();
		int numGridCells = (gridData) ? gridData->getNumActiveCells() : B3_FLUID_HASH_GRID_NUM_CELLS;
		
		if(!numRigidBodies || !numFluidParticles || !numGridCells) return;
		
		if( m_firstContactIndexPerParticle.size() < numFluidParticles ) m_firstContactIndexPerParticle.resize(numFluidParticles);
		if( m_numContactsPerParticle.size() < numFluidParticles ) m_numContactsPerParticle.resize(numFluidParticles);
		
		if( m_firstContactIndexPerRigid.size() < numRigidBodies ) m_firstContactIndexPerRigid.resize(numRigidBodies);
		if( m_numContactsPerRigid.size() < numRigidBodies ) m_numContactsPerRigid.resize(numRigidBodies);
		
		//Set number of pairs, contacts, and midphase pairs to 0 (numPairs == numContacts)
		{
			const int reset = 0;
			m_numPairs.copyFromHostPointer(&reset, 1);
			m_numMidphasePairs.copyFromHostPointer(&reset, 1);
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
					b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
					
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
				
				//
				int numLargeRigidAabbs = -1;
				m_numLargeAabbRigid.copyToHostPointer(&numLargeRigidAabbs, 1);
				clFinish(m_commandQueue);
				
				if(numLargeRigidAabbs > MAX_LARGE_AABB_RIGIDS) m_numLargeAabbRigid.copyFromHostPointer(&MAX_LARGE_AABB_RIGIDS, 1);
				clFinish(m_commandQueue);
			}
		
			//Large Rigid AABB intersection test with fluid AABBs
			{
				B3_PROFILE("m_fluidLargeRigidBroadphaseKernel");
				
				b3BufferInfoCL bufferInfo[] = 
				{
					b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
					b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
				
					b3BufferInfoCL( rigidBodyData.m_worldSpaceAabbs ),
					b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
					b3BufferInfoCL( rigidBodyData.m_collidables ),
					
					b3BufferInfoCL( m_numLargeAabbRigid.getBufferCL() ),
					b3BufferInfoCL( m_largeAabbRigidIndicies.getBufferCL() ),
					
					b3BufferInfoCL( m_numPairs.getBufferCL() ),
					b3BufferInfoCL( m_numMidphasePairs.getBufferCL() ),
					b3BufferInfoCL( m_pairs.getBufferCL() ),
					b3BufferInfoCL( m_midphasePairs.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_fluidLargeRigidBroadphaseKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(MAX_LARGE_AABB_RIGIDS);
				launcher.setConst(MAX_FLUID_RIGID_PAIRS);
				launcher.setConst(MAX_MIDPHASE_PAIRS);
				launcher.setConst(numFluidParticles);
				
				launcher.launch1D(numFluidParticles);
				clFinish(m_commandQueue);
				
				//
				int numPairs = -1;
				int numMidphasePairs = -1;
				m_numPairs.copyToHostPointer(&numPairs, 1);
				m_numMidphasePairs.copyToHostPointer(&numMidphasePairs, 1);
				clFinish(m_commandQueue);
				
				if(numPairs > MAX_FLUID_RIGID_PAIRS) m_numPairs.copyFromHostPointer(&MAX_FLUID_RIGID_PAIRS, 1);
				if(numMidphasePairs > MAX_MIDPHASE_PAIRS) m_numMidphasePairs.copyFromHostPointer(&MAX_MIDPHASE_PAIRS, 1);
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
				b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
				
				b3BufferInfoCL( gridData->m_activeCells.getBufferCL() ),
				b3BufferInfoCL( gridData->m_cellContents.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_worldSpaceAabbs ),
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_collidables ),
				
				b3BufferInfoCL( m_numPairs.getBufferCL() ),
				b3BufferInfoCL( m_numMidphasePairs.getBufferCL() ),
				b3BufferInfoCL( m_pairs.getBufferCL() ),
				b3BufferInfoCL( m_midphasePairs.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_fluidSmallRigidBroadphaseKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(MAX_FLUID_RIGID_PAIRS);
			launcher.setConst(MAX_MIDPHASE_PAIRS);
			launcher.setConst(numGridCells);
			launcher.setConst(numRigidBodies);
			
			launcher.launch1D(numRigidBodies);
			clFinish(m_commandQueue);
				
			//
			int numPairs = -1;
			int numMidphasePairs = -1;
			m_numPairs.copyToHostPointer(&numPairs, 1);
			m_numMidphasePairs.copyToHostPointer(&numMidphasePairs, 1);
			clFinish(m_commandQueue);
			
			if(1) printf("numMidphasePairs/max: %d / %d\n", numMidphasePairs, MAX_MIDPHASE_PAIRS);
			
			if(numPairs > MAX_FLUID_RIGID_PAIRS) m_numPairs.copyFromHostPointer(&MAX_FLUID_RIGID_PAIRS, 1);
			if(numMidphasePairs > MAX_MIDPHASE_PAIRS) m_numMidphasePairs.copyFromHostPointer(&MAX_MIDPHASE_PAIRS, 1);
			clFinish(m_commandQueue);
		}
		else if(moduloGridData)
		{
			b3Assert(0);
		
			/*
			B3_PROFILE("m_fluidSmallRigidBroadphaseModuloKernel");
			
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
				
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
			*/
		}
		else
		{
			b3Assert(0);	//No grid data
		}
		
		//Midphase - for triangle mesh and compound shapes
		{
			B3_PROFILE("m_fluidRigidMidphaseKernel");
			
			int numMidphasePairs = -1;
			m_numMidphasePairs.copyToHostPointer(&numMidphasePairs, 1);
			clFinish(m_commandQueue);
		
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_collidables ),
				
				b3BufferInfoCL( rigidBodyData.m_bvhInfo ),
				b3BufferInfoCL( rigidBodyData.m_bvhSubtreeInfo ),
				b3BufferInfoCL( rigidBodyData.m_bvhNodes ),
				
				b3BufferInfoCL( m_midphasePairs.getBufferCL() ),
				
				b3BufferInfoCL( m_numPairs.getBufferCL() ),
				b3BufferInfoCL( m_pairs.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_fluidRigidMidphaseKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(MAX_FLUID_RIGID_PAIRS);
			launcher.setConst(numMidphasePairs);
			
			launcher.launch1D(numMidphasePairs);
			clFinish(m_commandQueue);
			
			//
			int numPairs = -1;
			m_numPairs.copyToHostPointer(&numPairs, 1);
			clFinish(m_commandQueue);
			
			if(numPairs > MAX_FLUID_RIGID_PAIRS) m_numPairs.copyFromHostPointer(&MAX_FLUID_RIGID_PAIRS, 1);
			clFinish(m_commandQueue);
		}
		
		//Optimization - sort the fluid-rigid pairs by the rigid body shape type to reduce divergence (optional)
		{
		}
		
		//Narrowphase
		{
			B3_PROFILE("m_fluidRigidNarrowphaseKernel");
			
			int numFluidRigidPairs = -1;
			m_numPairs.copyToHostPointer(&numFluidRigidPairs, 1);
			clFinish(m_commandQueue);
			
			if(1) printf("numFluidRigidPairs/max: %d / %d\n", numFluidRigidPairs, MAX_FLUID_RIGID_PAIRS);
			
			//Perform sphere-(convex hull, triangle, plane, or sphere) collision for each pair 
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
				b3BufferInfoCL( fluidData->m_position.getBufferCL() ),
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
			launcher.setConst(numFluidRigidPairs);
			
			launcher.launch1D(numFluidRigidPairs);
			clFinish(m_commandQueue);
		}
		
		//m_resolveFluidRigidCollisionsKernel, executed below, overwrites fluid particle velocities(m_velocity);
		//the pre-collision particle velocities are needed to calculate impulses for the dynamic rigid bodies
		m_fluidVelocities.copyFromOpenCLArray(fluidData->m_velocity);
		
		//Sort pairs and contacts by fluid particle index
		{
			B3_PROFILE("sort fluid-rigid pairs and contacts");
		
			//	duplicated (see above)
			int numFluidRigidPairs = -1;
			m_numPairs.copyToHostPointer(&numFluidRigidPairs, 1);
			clFinish(m_commandQueue);
				
			//Load each entry in m_sortByFluidIndexData with:
			//b3SortData.m_key == fluid particle index (value to sort by)
			//b3SortData.m_value == fluid-rigid pair index
			{
				B3_PROFILE("m_loadSortDataKernel");
				
				m_sortByFluidIndexData.resize(numFluidRigidPairs);
				
				b3BufferInfoCL bufferInfo[] = 
				{
					b3BufferInfoCL( m_pairs.getBufferCL() ),
					b3BufferInfoCL( m_sortByFluidIndexData.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_loadSortDataKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidRigidPairs);
				
				launcher.launch1D(numFluidRigidPairs);
				clFinish(m_commandQueue);
			}
			
			{
				B3_PROFILE("radix sort m_sortByFluidIndexData");
				m_radixSorter.execute(m_sortByFluidIndexData, 32);
				
				clFinish(m_commandQueue);
			}
			
			{
				B3_PROFILE("duplicate m_pairs and m_fluidRigidContacts for rearrange");
				
				//m_pairs.copyToCL( m_pairsTemp.getBufferCL(), numFluidRigidPairs );
				//m_fluidRigidContacts.copyToCL( m_fluidRigidContactsTemp.getBufferCL(), numFluidRigidPairs );
			
				m_pairsTemp.copyFromOpenCLArray(m_pairs);
				m_fluidRigidContactsTemp.copyFromOpenCLArray(m_fluidRigidContacts);
				clFinish(m_commandQueue);
			}
			
			//Rearrange m_pairs and m_fluidRigidContacts
			{
				B3_PROFILE("m_rearrangePairsAndContactsKernel");
			
				b3BufferInfoCL bufferInfo[] = 
				{
					b3BufferInfoCL( m_sortByFluidIndexData.getBufferCL() ),
					
					b3BufferInfoCL( m_pairsTemp.getBufferCL() ),
					b3BufferInfoCL( m_fluidRigidContactsTemp.getBufferCL() ),
					
					b3BufferInfoCL( m_pairs.getBufferCL() ),
					b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_rearrangePairsAndContactsKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidRigidPairs);
				
				launcher.launch1D(numFluidRigidPairs);
				clFinish(m_commandQueue);
			}
		}
		
		//Use atomic_min, atomic_inc to find the range of contact indices for each fluid particle
		{
			B3_PROFILE("Find contact index ranges");
		
			//	duplicated (see above x2)
			int numFluidRigidPairs = -1;
			m_numPairs.copyToHostPointer(&numFluidRigidPairs, 1);
			clFinish(m_commandQueue);
				
			{
				const int firstContactIndex = numFluidRigidPairs;
				const int numContacts = 0;
				m_fill.execute(m_firstContactIndexPerParticle, firstContactIndex, numFluidParticles, 0);
				m_fill.execute(m_numContactsPerParticle, numContacts, numFluidParticles, 0);
				clFinish(m_commandQueue);
			}
			
			{
				b3BufferInfoCL bufferInfo[] = 
				{
					b3BufferInfoCL( m_pairs.getBufferCL() ),
					b3BufferInfoCL( m_firstContactIndexPerParticle.getBufferCL() ),
					b3BufferInfoCL( m_numContactsPerParticle.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_findPerParticleContactRangeKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numFluidRigidPairs);
				
				launcher.launch1D(numFluidRigidPairs);
				clFinish(m_commandQueue);
			}
		}
		
		//Resolve Collisions - apply impulses to fluid particles
		{
			B3_PROFILE("m_resolveFluidRigidCollisionsKernel");
		
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
				
				b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
				b3BufferInfoCL( rigidBodyData.m_rigidBodyInertias ),
				
				b3BufferInfoCL( m_pairs.getBufferCL() ),
				b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() ),
				b3BufferInfoCL( m_firstContactIndexPerParticle.getBufferCL() ),
				b3BufferInfoCL( m_numContactsPerParticle.getBufferCL() ),
					
				b3BufferInfoCL( fluidData->m_velocity.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_resolveFluidRigidCollisionsKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numFluidParticles);
			
			launcher.launch1D(numFluidParticles);
			clFinish(m_commandQueue);
		}
		
		const bool APPLY_IMPULSES_TO_RIGID_BODIES = true;
		if(APPLY_IMPULSES_TO_RIGID_BODIES)
		{
			//	duplicated (see above x3)
			int numFluidRigidPairs = -1;
			m_numPairs.copyToHostPointer(&numFluidRigidPairs, 1);
			clFinish(m_commandQueue);
			
			//Since applying impulses simultaneously(in 1 kernel) to both fluid and rigid body would require
			//syncronization between threads, we instead run 2 kernels - the first iterates through all contacts with
			//1 thread per fluid particle, and the second does the same with 1 thread per rigid body.
			
			//Clear rigid side contacts
			{
				const int firstContactIndex = numFluidRigidPairs;
				const int numContacts = 0;
				m_fill.execute(m_firstContactIndexPerRigid, firstContactIndex, numRigidBodies, 0);
				m_fill.execute(m_numContactsPerRigid, numContacts, numRigidBodies, 0);
				clFinish(m_commandQueue);
			}
			
			//Map contacts, which are sorted by fluid particle index, to each rigid body.
			//Rather than duplicating the contact data and sorting that by the rigid body index,
			//create an array of {rigid index, contact index} pairs and sort that instead.
			{
				//Load each entry in m_fluidToRigidMap with:
				//b3SortData.m_key == rigid body index (value to sort by)
				//b3SortData.m_value == contact index (same as pair index)
				{
					B3_PROFILE("m_mapFluidToRigidContactsKernel");
					
					m_fluidToRigidMap.resize(numFluidRigidPairs);
					
					b3BufferInfoCL bufferInfo[] = 
					{
						b3BufferInfoCL( m_pairs.getBufferCL() ),
						b3BufferInfoCL( m_fluidToRigidMap.getBufferCL() )
					};
					
					b3LauncherCL launcher(m_commandQueue, m_mapFluidToRigidContactsKernel);
					launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
					launcher.setConst(numFluidRigidPairs);
					
					launcher.launch1D(numFluidRigidPairs);
					clFinish(m_commandQueue);
				}
				
				{
					B3_PROFILE("radix sort m_fluidToRigidMap");
					m_radixSorter.execute(m_fluidToRigidMap, 32);
					
					clFinish(m_commandQueue);
				}
				
				//Use atomic_min and atomic_inc to find the range of indices in m_fluidToRigidMap for each rigid body
				{
					b3BufferInfoCL bufferInfo[] = 
					{
						b3BufferInfoCL( m_fluidToRigidMap.getBufferCL() ),
						b3BufferInfoCL( m_firstContactIndexPerRigid.getBufferCL() ),
						b3BufferInfoCL( m_numContactsPerRigid.getBufferCL() )
					};
					
					b3LauncherCL launcher(m_commandQueue, m_findPerRigidContactRangeKernel);
					launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
					launcher.setConst(numFluidRigidPairs);
					
					launcher.launch1D(numFluidRigidPairs);
					clFinish(m_commandQueue);
				}
			}
			
			//Resolve Collisions - apply impulses to rigid bodies
			{
				B3_PROFILE("m_resolveRigidFluidCollisionsKernel");
					
				b3BufferInfoCL bufferInfo[] = 
				{ 
					b3BufferInfoCL( fluidData->m_parameters.getBufferCL() ),
					
					b3BufferInfoCL( rigidBodyData.m_rigidBodies ),
					b3BufferInfoCL( rigidBodyData.m_rigidBodyInertias ),
					
					b3BufferInfoCL( m_pairs.getBufferCL() ),
					b3BufferInfoCL( m_fluidRigidContacts.getBufferCL() ),
					b3BufferInfoCL( m_fluidToRigidMap.getBufferCL() ),
					b3BufferInfoCL( m_firstContactIndexPerRigid.getBufferCL() ),
					b3BufferInfoCL( m_numContactsPerRigid.getBufferCL() ),
					
					b3BufferInfoCL( m_fluidVelocities.getBufferCL() )
				};
				
				b3LauncherCL launcher(m_commandQueue, m_resolveRigidFluidCollisionsKernel);
				launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(numRigidBodies);
				
				launcher.launch1D(numRigidBodies);
				clFinish(m_commandQueue);
			}
		}
		
		clFinish(m_commandQueue);
	}
	
};

#endif
