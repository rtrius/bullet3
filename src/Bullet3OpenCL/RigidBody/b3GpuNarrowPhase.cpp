#include "b3GpuNarrowPhase.h"


#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3ConvexPolyhedronData.h"
#include "Bullet3OpenCL/NarrowphaseCollision/b3ConvexHullContact.h"
#include <string.h>
#include "Bullet3Collision/NarrowPhaseCollision/b3Config.h"
#include "Bullet3OpenCL/NarrowphaseCollision/b3BvhInfo.h"

#include "b3GpuNarrowPhaseInternalData.h"
#include "Bullet3OpenCL/NarrowphaseCollision/b3QuantizedBvh.h"

b3GpuNarrowPhase::b3GpuNarrowPhase(cl_context ctx, cl_device_id device, cl_command_queue queue, const b3Config& config)
	:m_data(0), m_maxTriConvexPairCapacity(256 * 1024),
m_context(ctx),
m_device(device),
m_queue(queue)
{
	m_data = new b3GpuNarrowPhaseInternalData();
	memset(m_data,0,sizeof(b3GpuNarrowPhaseInternalData));
    
	m_data->m_config = config;
	
	m_data->m_gpuSatCollision = new GpuSatCollision(ctx,device,queue);
	
	m_data->m_triangleConvexPairs = new b3OpenCLArray<b3Int4>(m_context,m_queue, config.m_maxTriConvexPairCapacity);
	m_data->m_worldVertsB1GPU = new b3OpenCLArray<b3Vector3>(ctx,queue,config.m_maxConvexBodies*config.m_maxVerticesPerFace);
    m_data->m_clippingFacesOutGPU = new  b3OpenCLArray<b3Int4>(ctx,queue,config.m_maxConvexBodies);
    m_data->m_worldNormalsAGPU = new  b3OpenCLArray<b3Vector3>(ctx,queue,config.m_maxConvexBodies);
	m_data->m_worldVertsA1GPU = new b3OpenCLArray<b3Vector3>(ctx,queue,config.m_maxConvexBodies*config.m_maxVerticesPerFace);
    m_data->m_worldVertsB2GPU = new  b3OpenCLArray<b3Vector3>(ctx,queue,config.m_maxConvexBodies*config.m_maxVerticesPerFace);
}

b3GpuNarrowPhase::~b3GpuNarrowPhase()
{
	delete m_data->m_gpuSatCollision;
	
	delete m_data->m_triangleConvexPairs;
	delete m_data->m_worldVertsB1GPU;
    delete m_data->m_clippingFacesOutGPU;
    delete m_data->m_worldNormalsAGPU;
	delete m_data->m_worldVertsA1GPU;
    delete m_data->m_worldVertsB2GPU;
}

void b3GpuNarrowPhase::computeContacts(b3StateAabbs& worldSpaceAabbs, b3StateOverlappingPairs& pairs,
										b3StateRigidBodies& rigids, b3StateRigidCollidables& collidables,
										b3StateRigidContacts& contacts)
{
	//swap buffer
	contacts.m_currentContactBuffer = 1 - contacts.m_currentContactBuffer;

	int maxContactCapacity = m_data->m_config.m_maxContactCapacity;
	int compoundPairCapacity = m_data->m_config.m_compoundPairCapacity;
	int maxTriConvexPairCapacity = m_data->m_config.m_maxTriConvexPairCapacity;
	int nContactOut = 0;
	int numTriConvexPairsOut = 0;
	

	int numObjects = rigids.m_numRigidBodies;

	b3OpenCLArray<b3Aabb> clAabbArrayWorldSpace(m_context,m_queue);
	clAabbArrayWorldSpace.setFromOpenCLBuffer( worldSpaceAabbs.m_aabbsGpu.getBufferCL(), numObjects );

	b3OpenCLArray<b3Aabb> clAabbArrayLocalSpace(m_context,m_queue);
	clAabbArrayLocalSpace.setFromOpenCLBuffer( collidables.m_localShapeAabbGpu.getBufferCL(), numObjects );

	m_data->m_gpuSatCollision->computeConvexConvexContactsGPUSAT(
		&pairs.m_overlappingPairsGpu, pairs.m_overlappingPairsGpu.size(),
		&rigids.m_bodyBufferGPU,
		contacts.m_pContactBuffersGPU[contacts.m_currentContactBuffer],
		nContactOut,
		contacts.m_pContactBuffersGPU[1 - contacts.m_currentContactBuffer],
		maxContactCapacity,
		compoundPairCapacity,
		collidables.m_convexPolyhedraGpu,
		collidables.m_convexVerticesGpu,
		collidables.m_uniqueEdgesGpu,
		collidables.m_convexFacesGpu,
		collidables.m_convexIndicesGpu,
		collidables.m_collidablesGpu,
		collidables.m_childShapesGpu,
		clAabbArrayWorldSpace,
		clAabbArrayLocalSpace,
		*m_data->m_worldVertsB1GPU,
		*m_data->m_clippingFacesOutGPU,
		*m_data->m_worldNormalsAGPU,
		*m_data->m_worldVertsA1GPU,
		*m_data->m_worldVertsB2GPU,
		collidables.m_bvhData,
		&collidables.m_treeNodesGpu,
		&collidables.m_subTreesGpu,
		&collidables.m_bvhInfoGpu,
		numObjects,
		maxTriConvexPairCapacity,
		*m_data->m_triangleConvexPairs,
		numTriConvexPairsOut
		);
}



