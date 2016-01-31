
#include "b3GpuRaycast.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3RigidBodyData.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3ConvexPolyhedronData.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhaseInternalData.h"


#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3FillCL.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3RadixSort32CL.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3GpuBroadphaseInterface.h"
#include "Bullet3OpenCL/BroadphaseCollision/b3GpuParallelLinearBvh.h"

#include "Bullet3OpenCL/Raycast/kernels/rayCastKernels.h"


#define B3_RAYCAST_PATH "src/Bullet3OpenCL/Raycast/kernels/rayCastKernels.cl"



struct b3GpuRaycastInternalData
{
	cl_context m_context;
	cl_device_id m_device;
	cl_command_queue  m_q;
	cl_kernel m_raytraceKernel;
	cl_kernel m_raytracePairsKernel;
	cl_kernel m_findRayRigidPairIndexRanges;
	cl_kernel m_findFirstHitPerRayKernel;
	
	cl_kernel m_plbvhRayTraverseFirstHitKernel;
	cl_kernel m_plbvhLargeAabbRayTestFirstHitKernel;
	
	b3GpuParallelLinearBvh* m_plbvh;
	b3RadixSort32CL* m_radixSorter;
	b3FillCL* m_fill;
	
	//1 element per ray
	b3OpenCLArray<b3RayInfo>* m_gpuRays;
	b3OpenCLArray<b3RayHit>* m_gpuHitResults;
	b3OpenCLArray<int>* m_firstRayRigidPairIndexPerRay;
	b3OpenCLArray<int>* m_numRayRigidPairsPerRay;
	
	//1 element per (ray index, rigid index) pair, where the ray intersects with the rigid's AABB
	b3OpenCLArray<int>* m_gpuNumRayRigidPairs;
	b3OpenCLArray<b3Int2>* m_gpuRayRigidPairs;					//x == ray index, y == rigid index
	b3OpenCLArray<b3Vector3>* m_normalAndHitFractionPerPair;	//normal == x,y,z, hit fraction == w
};

b3GpuRaycast::b3GpuRaycast(cl_context ctx,cl_device_id device, cl_command_queue  q)
{
	m_data = new b3GpuRaycastInternalData;
	m_data->m_context = ctx;
	m_data->m_device = device;
	m_data->m_q = q;
	m_data->m_raytraceKernel = 0;
	m_data->m_raytracePairsKernel = 0;
	m_data->m_findRayRigidPairIndexRanges = 0;
	m_data->m_findFirstHitPerRayKernel = 0;
	
	m_data->m_plbvhRayTraverseFirstHitKernel = 0;
	m_data->m_plbvhLargeAabbRayTestFirstHitKernel = 0;

	m_data->m_plbvh = new b3GpuParallelLinearBvh(ctx, device, q);
	m_data->m_radixSorter = new b3RadixSort32CL(ctx, device, q);
	m_data->m_fill = new b3FillCL(ctx, device, q);
	
	m_data->m_gpuRays = new b3OpenCLArray<b3RayInfo>(ctx, q);
	m_data->m_gpuHitResults = new b3OpenCLArray<b3RayHit>(ctx, q);
	m_data->m_firstRayRigidPairIndexPerRay = new b3OpenCLArray<int>(ctx, q);
	m_data->m_numRayRigidPairsPerRay = new b3OpenCLArray<int>(ctx, q);
	
	m_data->m_gpuNumRayRigidPairs = new b3OpenCLArray<int>(ctx, q);
	m_data->m_gpuRayRigidPairs = new b3OpenCLArray<b3Int2>(ctx, q);
	m_data->m_normalAndHitFractionPerPair = new b3OpenCLArray<b3Vector3>(ctx, q);

	{
		cl_int errNum=0;
		cl_program prog = b3OpenCLUtils::compileCLProgramFromString(m_data->m_context,m_data->m_device,rayCastKernelCL,&errNum,"",B3_RAYCAST_PATH);
		b3Assert(errNum==CL_SUCCESS);
		m_data->m_raytraceKernel = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,rayCastKernelCL, "rayCastKernel",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);
		m_data->m_raytracePairsKernel = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,rayCastKernelCL, "rayCastPairsKernel",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);
		m_data->m_findRayRigidPairIndexRanges = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,rayCastKernelCL, "findRayRigidPairIndexRanges",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);
		m_data->m_findFirstHitPerRayKernel = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,rayCastKernelCL, "findFirstHitPerRay",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);
		
		m_data->m_plbvhRayTraverseFirstHitKernel = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,rayCastKernelCL, "plbvhRayTraverseFirstHit",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);
		m_data->m_plbvhLargeAabbRayTestFirstHitKernel = b3OpenCLUtils::compileCLKernelFromString(m_data->m_context, m_data->m_device,rayCastKernelCL, "plbvhLargeAabbRayTestFirstHit",&errNum,prog);
		b3Assert(errNum==CL_SUCCESS);
		
		clReleaseProgram(prog);
	}

	//
	m_maxRayRigidPairs = 16384;	//Arbitrary value
}

b3GpuRaycast::~b3GpuRaycast()
{
	clReleaseKernel(m_data->m_raytraceKernel);
	clReleaseKernel(m_data->m_raytracePairsKernel);
	clReleaseKernel(m_data->m_findRayRigidPairIndexRanges);
	clReleaseKernel(m_data->m_findFirstHitPerRayKernel);
	
	clReleaseKernel(m_data->m_plbvhRayTraverseFirstHitKernel);
	clReleaseKernel(m_data->m_plbvhLargeAabbRayTestFirstHitKernel);
	
	delete m_data->m_plbvh;
	delete m_data->m_radixSorter;
	delete m_data->m_fill;
	
	delete m_data->m_gpuRays;
	delete m_data->m_gpuHitResults;
	delete m_data->m_firstRayRigidPairIndexPerRay;
	delete m_data->m_numRayRigidPairsPerRay;
	
	delete m_data->m_gpuNumRayRigidPairs;
	delete m_data->m_gpuRayRigidPairs;
	delete m_data->m_normalAndHitFractionPerPair;
	
	delete m_data;
}

void b3GpuRaycast::castRays(const b3AlignedObjectArray<b3RayInfo>& rays, b3AlignedObjectArray<b3RayHit>& hitResults,
	int numBodies, b3OpenCLArray<b3RigidBodyData>& rigidBodies, b3StateRigidCollidables& collidables, b3StateAabbs& aabbs)
{
	if(0)
	{
		castRaysUsingPairs(rays, hitResults, numBodies, rigidBodies, collidables, aabbs);
		return;
	}

	B3_PROFILE("b3GpuRaycast::castRays()");
	
	{
		B3_PROFILE("raycast copyFromHost");
		m_data->m_gpuRays->copyFromHost(rays);
		m_data->m_gpuHitResults->copyFromHost(hitResults);
	}
	
	//
	const bool USE_BRUTE_FORCE_RAYCAST = false;
	if(USE_BRUTE_FORCE_RAYCAST)
	{
		B3_PROFILE("raycast launch1D");

		b3LauncherCL launcher(m_data->m_q,m_data->m_raytraceKernel,"m_raytraceKernel");
		int numRays = rays.size();
		launcher.setConst(numRays);

		launcher.setBuffer(m_data->m_gpuRays->getBufferCL());
		launcher.setBuffer(m_data->m_gpuHitResults->getBufferCL());

		launcher.setConst(numBodies);
		launcher.setBuffer(rigidBodies.getBufferCL());
		launcher.setBuffer(collidables.m_collidablesGpu.getBufferCL());
		launcher.setBuffer(collidables.m_convexFacesGpu.getBufferCL());
		launcher.setBuffer(collidables.m_convexPolyhedraGpu.getBufferCL());
		
		launcher.launch1D(numRays);
		clFinish(m_data->m_q);
	}
	else
	{
		m_data->m_plbvh->build(aabbs.m_aabbsGpu, aabbs.m_smallAabbsMappingGpu, aabbs.m_largeAabbsMappingGpu);
		
		//Simultaneously traverse BVH and find first/closest hit(computes ray-rigid intersection and ray-AABB intersection)
		{
			B3_PROFILE("m_plbvhRayTraverseFirstHitKernel");
			
			int numRays = hitResults.size();
		
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( m_data->m_plbvh->getLeafNodeAabbs().getBufferCL() ),
				
				b3BufferInfoCL( m_data->m_plbvh->getRootNodeIndex().getBufferCL() ),
				b3BufferInfoCL( m_data->m_plbvh->getInternalNodeChildNodes().getBufferCL() ),
				b3BufferInfoCL( m_data->m_plbvh->getInternalNodeAabbs().getBufferCL() ),
				b3BufferInfoCL( m_data->m_plbvh->getMortonCodesAndAabbIndices().getBufferCL() ),
		
				b3BufferInfoCL( rigidBodies.getBufferCL() ),
				b3BufferInfoCL( collidables.m_collidablesGpu.getBufferCL() ),
				b3BufferInfoCL( collidables.m_convexFacesGpu.getBufferCL() ),
				b3BufferInfoCL( collidables.m_convexPolyhedraGpu.getBufferCL() ),
				
				b3BufferInfoCL( m_data->m_gpuRays->getBufferCL() ),
				b3BufferInfoCL( m_data->m_gpuHitResults->getBufferCL() )
			};
			
			const int UNUSED = -1;
			b3LauncherCL launcher(m_data->m_q, m_data->m_plbvhRayTraverseFirstHitKernel, "m_plbvhRayTraverseFirstHitKernel");
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(UNUSED);
			launcher.setConst(numRays);
			
			launcher.launch1D(numRays);
			clFinish(m_data->m_q);
		}
		
		//Compute closest hit for large AABBs(only small AABBs are in the BVH)
		{
			B3_PROFILE("m_plbvhLargeAabbRayTestFirstHitKernel");
			
			b3OpenCLArray<b3SapAabb>& largeAabbs = m_data->m_plbvh->getLargeAabbs();
			
			int numRays = hitResults.size();
			int numLargeAabbs = largeAabbs.size();
		
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( largeAabbs.getBufferCL() ),

				b3BufferInfoCL(rigidBodies.getBufferCL()),
				b3BufferInfoCL(collidables.m_collidablesGpu.getBufferCL()),
				b3BufferInfoCL(collidables.m_convexFacesGpu.getBufferCL()),
				b3BufferInfoCL(collidables.m_convexPolyhedraGpu.getBufferCL()),
				
				b3BufferInfoCL( m_data->m_gpuRays->getBufferCL() ),
				b3BufferInfoCL( m_data->m_gpuHitResults->getBufferCL() )
			};
			
			b3LauncherCL launcher(m_data->m_q, m_data->m_plbvhLargeAabbRayTestFirstHitKernel, "m_plbvhLargeAabbRayTestFirstHitKernel");
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numLargeAabbs);
			launcher.setConst(numRays);
			
			launcher.launch1D(numRays);
			clFinish(m_data->m_q);
		}
	}
	
	//
	{
		B3_PROFILE("raycast copyToHost");
		m_data->m_gpuHitResults->copyToHost(hitResults);
	}
}

void b3GpuRaycast::castRaysUsingPairs(const b3AlignedObjectArray<b3RayInfo>& rays, b3AlignedObjectArray<b3RayHit>& hitResults,
	int numBodies, b3OpenCLArray<b3RigidBodyData>& rigidBodies, b3StateRigidCollidables& collidables, b3StateAabbs& aabbs)
{
	B3_PROFILE("b3GpuRaycast::castRaysUsingPairs()");

	{
		B3_PROFILE("raycast copyFromHost");
		m_data->m_gpuRays->copyFromHost(rays);
		m_data->m_gpuHitResults->copyFromHost(hitResults);
		
	}
	
	{
		int numRays = hitResults.size();
		{
			m_data->m_firstRayRigidPairIndexPerRay->resize(numRays);
			m_data->m_numRayRigidPairsPerRay->resize(numRays);
			
			m_data->m_gpuNumRayRigidPairs->resize(1);
			
			m_data->m_gpuRayRigidPairs->resize(m_maxRayRigidPairs);
			m_data->m_normalAndHitFractionPerPair->resize(m_maxRayRigidPairs);
		}

		m_data->m_plbvh->build(aabbs.m_aabbsGpu, aabbs.m_smallAabbsMappingGpu, aabbs.m_largeAabbsMappingGpu);

		m_data->m_plbvh->testRaysAgainstBvhAabbs(*m_data->m_gpuRays, *m_data->m_gpuNumRayRigidPairs, *m_data->m_gpuRayRigidPairs);
		
		int numRayRigidPairs = -1;
		m_data->m_gpuNumRayRigidPairs->copyToHostPointer(&numRayRigidPairs, 1);
		if( numRayRigidPairs > m_data->m_gpuRayRigidPairs->size() )
		{
			numRayRigidPairs = m_data->m_gpuRayRigidPairs->size();
			m_data->m_gpuNumRayRigidPairs->copyFromHostPointer(&numRayRigidPairs, 1);
		}
		
		m_data->m_gpuRayRigidPairs->resize(numRayRigidPairs);	//Radix sort needs b3OpenCLArray::size() to be correct
		
		//Sort ray-rigid pairs by ray index
		{
			B3_PROFILE("sort ray-rigid pairs");
			m_data->m_radixSorter->execute( *reinterpret_cast< b3OpenCLArray<b3SortData>* >(m_data->m_gpuRayRigidPairs) );
			clFinish(m_data->m_q);
		}
		
		//Determine number of ray-rigid pairs associated with each ray, and first/lowest ray-rigid pair index
		{
			B3_PROFILE("detect ray-rigid pair index ranges");
			
			{
				B3_PROFILE("reset ray-rigid pair index ranges");
				
				m_data->m_fill->execute(*m_data->m_firstRayRigidPairIndexPerRay, numRayRigidPairs, numRays);	//atomic_min used to find first index
				m_data->m_fill->execute(*m_data->m_numRayRigidPairsPerRay, 0, numRays);
				clFinish(m_data->m_q);
			}
			
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( m_data->m_gpuRayRigidPairs->getBufferCL() ),
				
				b3BufferInfoCL( m_data->m_firstRayRigidPairIndexPerRay->getBufferCL() ),
				b3BufferInfoCL( m_data->m_numRayRigidPairsPerRay->getBufferCL() )
			};
			
			b3LauncherCL launcher(m_data->m_q, m_data->m_findRayRigidPairIndexRanges, "m_findRayRigidPairIndexRanges");
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numRayRigidPairs);
			
			launcher.launch1D(numRayRigidPairs);
			clFinish(m_data->m_q);
		}
		
		//Perform narrowphase/intersection test for all ray-rigid pairs
		{
			B3_PROFILE("ray-rigid intersection");
			
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( m_data->m_gpuRays->getBufferCL() ),
				b3BufferInfoCL( m_data->m_gpuHitResults->getBufferCL() ),

				b3BufferInfoCL(rigidBodies.getBufferCL()),
				b3BufferInfoCL(collidables.m_collidablesGpu.getBufferCL()),
				b3BufferInfoCL(collidables.m_convexFacesGpu.getBufferCL()),
				b3BufferInfoCL(collidables.m_convexPolyhedraGpu.getBufferCL()),
				
				b3BufferInfoCL( m_data->m_gpuRayRigidPairs->getBufferCL() ),
				b3BufferInfoCL( m_data->m_normalAndHitFractionPerPair->getBufferCL() )
			};
			
			b3LauncherCL launcher(m_data->m_q, m_data->m_raytracePairsKernel, "m_raytracePairsKernel");
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numRayRigidPairs);
			
			launcher.launch1D(numRayRigidPairs);
			clFinish(m_data->m_q);
		}
		
		//For each ray, find the ray-rigid pair with the nearest(smallest) hit fraction
		//and copy the result from m_data->m_normalAndHitFractionPerPair into m_data->m_gpuHitResults
		{
			B3_PROFILE("find first hit");
			
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( m_data->m_gpuRays->getBufferCL() ),
				b3BufferInfoCL( m_data->m_gpuRayRigidPairs->getBufferCL() ),
				b3BufferInfoCL( m_data->m_firstRayRigidPairIndexPerRay->getBufferCL() ),
				b3BufferInfoCL( m_data->m_numRayRigidPairsPerRay->getBufferCL() ),
				b3BufferInfoCL( m_data->m_normalAndHitFractionPerPair->getBufferCL() ),
				b3BufferInfoCL( m_data->m_gpuHitResults->getBufferCL() )
			};
			
			b3LauncherCL launcher(m_data->m_q, m_data->m_findFirstHitPerRayKernel, "m_findFirstHitPerRayKernel");
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numRays);
			
			launcher.launch1D(numRays);
			clFinish(m_data->m_q);
		}
	}
	
	//
	{
		B3_PROFILE("raycast copyToHost");
		m_data->m_gpuHitResults->copyToHost(hitResults);
	}
}



