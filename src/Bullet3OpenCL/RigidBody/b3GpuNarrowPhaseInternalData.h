
#ifndef B3_GPU_NARROWPHASE_INTERNAL_DATA_H
#define B3_GPU_NARROWPHASE_INTERNAL_DATA_H

#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3Collision/NarrowPhaseCollision/b3Config.h"

#include "Bullet3OpenCL/Initialize/b3OpenCLInclude.h"
#include "Bullet3Common/b3AlignedObjectArray.h"
#include "Bullet3Common/b3Vector3.h"


#include "Bullet3Common/shared/b3Int4.h"
#include "Bullet3Common/shared/b3Int2.h"


class b3ConvexUtility;

struct b3GpuNarrowPhaseInternalData
{
	//GPU_API_REDESIGN
	//temp? buffers for narrowphase
    b3OpenCLArray<b3Vector3>* m_worldVertsB1GPU;
    b3OpenCLArray<b3Int4>* m_clippingFacesOutGPU;
    b3OpenCLArray<b3Vector3>* m_worldNormalsAGPU;
    b3OpenCLArray<b3Vector3>* m_worldVertsA1GPU;
	b3OpenCLArray<b3Vector3>* m_worldVertsB2GPU;
	b3OpenCLArray<b3Int4>*			m_triangleConvexPairs;
    
	struct GpuSatCollision*	m_gpuSatCollision;
	    
  
	b3Config	m_config;
    

};

#endif //B3_GPU_NARROWPHASE_INTERNAL_DATA_H
