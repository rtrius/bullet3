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

//#include "fluidSphRigidCL.h"

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
	
	void load(b3GpuBroadphaseInterface* broadphase, b3GpuNarrowPhase* narrowPhase)
	{
		m_numRigidBodies = narrowPhase->getNumBodiesGpu();
		m_numRigidBodyInertias = narrowPhase->getNumBodyInertiasGpu();
		
		m_numWorldSpaceAabbs = broadphase->getAllAabbsGPU().size();
		m_rigidBodies = narrowPhase->getBodiesGpu();
		m_rigidBodyInertias = narrowPhase->getBodyInertiasGpu();
		m_worldSpaceAabbs = broadphase->getAllAabbsGPU().getBufferCL();		//if useDbvt == true in b3GpuRigidBodyPipeline.cpp this is incorrect
		
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
	static const int MAX_FLUID_RIGID_PAIRS; 	//Also the max number of fluid-rigid contacts
	static const int MAX_MIDPHASE_PAIRS;

	///Maximum number of rigid bodies with large AABBs(AABBs larger than the fluid's grid) that are considered for collision
	static const int MAX_LARGE_AABB_RIGIDS;	
	
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
	cl_kernel m_fluidSmallRigidBroadphaseModuloKernel;
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
	b3FluidSphRigidInteractorCL(cl_context context, cl_device_id device, cl_command_queue queue);
	
	virtual ~b3FluidSphRigidInteractorCL();
	
	//Either gridData or moduloGridData must be nonzero, but not both; it is used to determine which grid type is being used
	void interact(b3FluidSphOpenCL* fluidData,  b3FluidSortingGridOpenCL* gridData, 
					b3FluidHashGridOpenCL* moduloGridData, RigidBodyGpuData& rigidBodyData);
	
};

#endif
