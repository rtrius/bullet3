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
	:m_data(0), m_static0Index(-1), m_maxTriConvexPairCapacity(256 * 1024),
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







void b3GpuNarrowPhase::computeContacts(cl_mem broadphasePairs, int numBroadphasePairs, cl_mem aabbsWorldSpace, int numObjects)
{
#ifdef GPU_API_REDESIGN
	cl_mem aabbsLocalSpace = m_data->m_localShapeAABBGPU->getBufferCL();

	int nContactOut = 0;

	//swap buffer
	m_data->m_currentContactBuffer=1-m_data->m_currentContactBuffer;

	int curSize = m_data->m_pBufContactBuffersGPU[m_data->m_currentContactBuffer]->size();

	int maxTriConvexPairCapacity = m_data->m_config.m_maxTriConvexPairCapacity;
	int numTriConvexPairsOut=0;
	
	b3OpenCLArray<b3Int4> broadphasePairsGPU(m_context,m_queue);
	broadphasePairsGPU.setFromOpenCLBuffer(broadphasePairs,numBroadphasePairs);

	


	b3OpenCLArray<b3Aabb> clAabbArrayWorldSpace(this->m_context,this->m_queue);
	clAabbArrayWorldSpace.setFromOpenCLBuffer(aabbsWorldSpace,numObjects);

	b3OpenCLArray<b3Aabb> clAabbArrayLocalSpace(this->m_context,this->m_queue);
	clAabbArrayLocalSpace.setFromOpenCLBuffer(aabbsLocalSpace,numObjects);

	m_data->m_gpuSatCollision->computeConvexConvexContactsGPUSAT(
		&broadphasePairsGPU, numBroadphasePairs,
		m_data->m_bodyBufferGPU,
		m_data->m_pBufContactBuffersGPU[m_data->m_currentContactBuffer],
		nContactOut,
		m_data->m_pBufContactBuffersGPU[1-m_data->m_currentContactBuffer],
		m_data->m_config.m_maxContactCapacity,
		m_data->m_config.m_compoundPairCapacity,
		*m_data->m_convexPolyhedraGPU,
		*m_data->m_convexVerticesGPU,
		*m_data->m_uniqueEdgesGPU,
		*m_data->m_convexFacesGPU,
		*m_data->m_convexIndicesGPU,
		*m_data->m_collidablesGPU,
		*m_data->m_gpuChildShapes,
		clAabbArrayWorldSpace,
		clAabbArrayLocalSpace,
		*m_data->m_worldVertsB1GPU,
		*m_data->m_clippingFacesOutGPU,
		*m_data->m_worldNormalsAGPU,
		*m_data->m_worldVertsA1GPU,
		*m_data->m_worldVertsB2GPU,
		m_data->m_bvhData,
		m_data->m_treeNodesGPU,
		m_data->m_subTreesGPU,
		m_data->m_bvhInfoGPU,
		numObjects,
		maxTriConvexPairCapacity,
		*m_data->m_triangleConvexPairs,
		numTriConvexPairsOut
		);

#endif
}



