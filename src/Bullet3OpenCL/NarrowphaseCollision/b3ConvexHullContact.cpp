/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2011 Advanced Micro Devices, Inc.  http://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

bool findSeparatingAxisOnGpu = true;
bool splitSearchSepAxisConcave = false;
bool splitSearchSepAxisConvex = true;
bool useMprGpu = false;//use mpr for edge-edge  (+contact point) or sat. Needs testing on main OpenCL platforms, before enabling...
bool bvhTraversalKernelGPU = true;
bool findConcaveSeparatingAxisKernelGPU = true;
bool reduceConcaveContactsOnGPU = true;//false;
bool reduceConvexContactsOnGPU = true;//false;
bool findConvexClippingFacesGPU = true;


static int myframecount=0;///for testing

///This file was written by Erwin Coumans
///Separating axis rest based on work from Pierre Terdiman, see
///And contact clipping based on work from Simon Hobbs

//#define B3_DEBUG_SAT_FACE

//#define CHECK_ON_HOST

#ifdef CHECK_ON_HOST
//#define PERSISTENT_CONTACTS_HOST
#endif

int b3g_actualSATPairTests=0;

#include "b3ConvexHullContact.h"
#include <string.h>//memcpy
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3ConvexPolyhedronData.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3MprPenetration.h"

#include "Bullet3OpenCL/NarrowphaseCollision/b3ContactCache.h"
#include "Bullet3Geometry/b3AabbUtil.h"

typedef b3AlignedObjectArray<b3Vector3> b3VertexArray;


#include <float.h> //for FLT_MAX
#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
//#include "AdlQuaternion.h"

#include "kernels/satKernels.h"
#include "kernels/mprKernels.h"

#include "kernels/satConcaveKernels.h"

#include "kernels/satClipHullContacts.h"
#include "kernels/bvhTraversal.h"
#include "kernels/primitiveContacts.h"


#include "Bullet3Geometry/b3AabbUtil.h"

#define BT_NARROWPHASE_SAT_PATH "src/Bullet3OpenCL/NarrowphaseCollision/kernels/sat.cl"
#define BT_NARROWPHASE_SAT_CONCAVE_PATH "src/Bullet3OpenCL/NarrowphaseCollision/kernels/satConcave.cl"

#define BT_NARROWPHASE_MPR_PATH "src/Bullet3OpenCL/NarrowphaseCollision/kernels/mpr.cl"


#define BT_NARROWPHASE_CLIPHULL_PATH "src/Bullet3OpenCL/NarrowphaseCollision/kernels/satClipHullContacts.cl"
#define BT_NARROWPHASE_BVH_TRAVERSAL_PATH "src/Bullet3OpenCL/NarrowphaseCollision/kernels/bvhTraversal.cl"
#define BT_NARROWPHASE_PRIMITIVE_CONTACT_PATH "src/Bullet3OpenCL/NarrowphaseCollision/kernels/primitiveContacts.cl"


#ifndef __global
#define __global
#endif

#ifndef __kernel
#define __kernel
#endif


#include "Bullet3Collision/NarrowPhaseCollision/shared/b3BvhTraversal.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3FindConcaveSatAxis.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3ClipFaces.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3NewContactReduction.h"



#define dot3F4 b3Dot

GpuSatCollision::GpuSatCollision(cl_context ctx,cl_device_id device, cl_command_queue  q )
:m_context(ctx),
m_device(device),
m_queue(q),
m_findSeparatingAxisKernel(0),
m_findSeparatingAxisVertexFaceKernel(0),
m_findSeparatingAxisEdgeEdgeKernel(0),
m_totalContactsOut(m_context, m_queue),
m_sepNormals(m_context, m_queue),
m_hasSeparatingNormals(m_context, m_queue),
m_concaveSepNormals(m_context, m_queue),
m_concaveHasSeparatingNormals(m_context,m_queue),
m_numConcavePairsOut(m_context, m_queue),
m_gpuCompoundPairs(m_context, m_queue),
m_gpuCompoundSepNormals(m_context, m_queue),
m_gpuHasCompoundSepNormals(m_context, m_queue),
m_numCompoundPairsOut(m_context, m_queue),
m_dmins(m_context,m_queue),
m_unitSphereDirections(m_context,m_queue)
{
	m_totalContactsOut.push_back(0);
	
	cl_int errNum=0;

	if (1)
	{
		const char* mprSrc = mprKernelsCL;
		
		const char* srcConcave = satConcaveKernelsCL;
		char flags[1024]={0};
//#ifdef CL_PLATFORM_INTEL
//		sprintf(flags,"-g -s \"%s\"","C:/develop/bullet3_experiments2/opencl/gpu_narrowphase/kernels/sat.cl");
//#endif
		m_mprPenetrationKernel  = 0;
		m_findSeparatingAxisUnitSphereKernel = 0;

		if (useMprGpu)
		{
			cl_program mprProg = b3OpenCLUtils::compileCLProgramFromString(m_context,m_device,mprSrc,&errNum,flags,BT_NARROWPHASE_MPR_PATH);
			b3Assert(errNum==CL_SUCCESS);
		
			m_mprPenetrationKernel  = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,mprSrc, "mprPenetrationKernel",&errNum,mprProg );
			b3Assert(m_mprPenetrationKernel);
			b3Assert(errNum==CL_SUCCESS);

			m_findSeparatingAxisUnitSphereKernel =  b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,mprSrc, "findSeparatingAxisUnitSphereKernel",&errNum,mprProg );
			b3Assert(m_findSeparatingAxisUnitSphereKernel);
            b3Assert(errNum==CL_SUCCESS);


			int numDirections = sizeof(unitSphere162)/sizeof(b3Vector3);
			m_unitSphereDirections.resize(numDirections);
			m_unitSphereDirections.copyFromHostPointer(unitSphere162,numDirections,0,true);


		}


		cl_program satProg = b3OpenCLUtils::compileCLProgramFromString(m_context,m_device,satKernelsCL,&errNum,flags,BT_NARROWPHASE_SAT_PATH);
		b3Assert(errNum==CL_SUCCESS);

		cl_program satConcaveProg = b3OpenCLUtils::compileCLProgramFromString(m_context,m_device,srcConcave,&errNum,flags,BT_NARROWPHASE_SAT_CONCAVE_PATH);
		b3Assert(errNum==CL_SUCCESS);

		m_findSeparatingAxisKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,satKernelsCL, "findSeparatingAxisKernel",&errNum,satProg );
		b3Assert(m_findSeparatingAxisKernel);
		b3Assert(errNum==CL_SUCCESS);


		m_findSeparatingAxisVertexFaceKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,satKernelsCL, "findSeparatingAxisVertexFaceKernel",&errNum,satProg );
		b3Assert(m_findSeparatingAxisVertexFaceKernel);

		m_findSeparatingAxisEdgeEdgeKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,satKernelsCL, "findSeparatingAxisEdgeEdgeKernel",&errNum,satProg );
		b3Assert(m_findSeparatingAxisVertexFaceKernel);


		m_findConcaveSeparatingAxisKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,satKernelsCL, "findConcaveSeparatingAxisKernel",&errNum,satProg );
		b3Assert(m_findConcaveSeparatingAxisKernel);
		b3Assert(errNum==CL_SUCCESS);
        
        m_findConcaveSeparatingAxisVertexFaceKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcConcave, "findConcaveSeparatingAxisVertexFaceKernel",&errNum,satConcaveProg );
		b3Assert(m_findConcaveSeparatingAxisVertexFaceKernel);
		b3Assert(errNum==CL_SUCCESS);
        
        m_findConcaveSeparatingAxisEdgeEdgeKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcConcave, "findConcaveSeparatingAxisEdgeEdgeKernel",&errNum,satConcaveProg );
		b3Assert(m_findConcaveSeparatingAxisEdgeEdgeKernel);
		b3Assert(errNum==CL_SUCCESS);
        
     
        
		
		m_findCompoundPairsKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,satKernelsCL, "findCompoundPairsKernel",&errNum,satProg );
		b3Assert(m_findCompoundPairsKernel);
		b3Assert(errNum==CL_SUCCESS);
		m_processCompoundPairsKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,satKernelsCL, "processCompoundPairsKernel",&errNum,satProg );
		b3Assert(m_processCompoundPairsKernel);
		b3Assert(errNum==CL_SUCCESS);
	}

	if (1)
	{
		const char* srcClip = satClipKernelsCL;

		char flags[1024]={0};
//#ifdef CL_PLATFORM_INTEL
//		sprintf(flags,"-g -s \"%s\"","C:/develop/bullet3_experiments2/opencl/gpu_narrowphase/kernels/satClipHullContacts.cl");
//#endif

		cl_program satClipContactsProg = b3OpenCLUtils::compileCLProgramFromString(m_context,m_device,srcClip,&errNum,flags,BT_NARROWPHASE_CLIPHULL_PATH);
		b3Assert(errNum==CL_SUCCESS);

		m_clipHullHullKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcClip, "clipHullHullKernel",&errNum,satClipContactsProg);
		b3Assert(errNum==CL_SUCCESS);

		m_clipCompoundsHullHullKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcClip, "clipCompoundsHullHullKernel",&errNum,satClipContactsProg);
		b3Assert(errNum==CL_SUCCESS);
		

        m_findClippingFacesKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcClip, "findClippingFacesKernel",&errNum,satClipContactsProg);
		b3Assert(errNum==CL_SUCCESS);

        m_clipFacesAndFindContacts = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcClip, "clipFacesAndFindContactsKernel",&errNum,satClipContactsProg);
		b3Assert(errNum==CL_SUCCESS);        

		m_clipHullHullConcaveConvexKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcClip, "clipHullHullConcaveConvexKernel",&errNum,satClipContactsProg);
		b3Assert(errNum==CL_SUCCESS);

//		m_extractManifoldAndAddContactKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcClip, "extractManifoldAndAddContactKernel",&errNum,satClipContactsProg);
	//	b3Assert(errNum==CL_SUCCESS);

        m_newContactReductionKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcClip,
                            "newContactReductionKernel",&errNum,satClipContactsProg);
		b3Assert(errNum==CL_SUCCESS);
	}
   else
	{
		m_clipHullHullKernel=0;
		m_clipCompoundsHullHullKernel = 0;
        m_findClippingFacesKernel = 0;
        m_newContactReductionKernel=0;
        m_clipFacesAndFindContacts = 0;
		m_clipHullHullConcaveConvexKernel = 0;
//		m_extractManifoldAndAddContactKernel = 0;
	}

	 if (1)
	{
		const char* srcBvh = bvhTraversalKernelCL;
		cl_program bvhTraversalProg = b3OpenCLUtils::compileCLProgramFromString(m_context,m_device,srcBvh,&errNum,"",BT_NARROWPHASE_BVH_TRAVERSAL_PATH);
		b3Assert(errNum==CL_SUCCESS);

		m_bvhTraversalKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,srcBvh, "bvhTraversalKernel",&errNum,bvhTraversalProg,"");
		b3Assert(errNum==CL_SUCCESS);

	}
        
	 {
		 const char* primitiveContactsSrc = primitiveContactsKernelsCL;
		cl_program primitiveContactsProg = b3OpenCLUtils::compileCLProgramFromString(m_context,m_device,primitiveContactsSrc,&errNum,"",BT_NARROWPHASE_PRIMITIVE_CONTACT_PATH);
		b3Assert(errNum==CL_SUCCESS);

		m_primitiveContactsKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,primitiveContactsSrc, "primitiveContactsKernel",&errNum,primitiveContactsProg,"");
		b3Assert(errNum==CL_SUCCESS);

		m_findConcaveSphereContactsKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,primitiveContactsSrc, "findConcaveSphereContactsKernel",&errNum,primitiveContactsProg );
		b3Assert(errNum==CL_SUCCESS);
		b3Assert(m_findConcaveSphereContactsKernel);

		m_processCompoundPairsPrimitivesKernel = b3OpenCLUtils::compileCLKernelFromString(m_context, m_device,primitiveContactsSrc, "processCompoundPairsPrimitivesKernel",&errNum,primitiveContactsProg,"");
		b3Assert(errNum==CL_SUCCESS);
		b3Assert(m_processCompoundPairsPrimitivesKernel);
		 
	 }
	

}

GpuSatCollision::~GpuSatCollision()
{
	
	if (m_findSeparatingAxisVertexFaceKernel)
		clReleaseKernel(m_findSeparatingAxisVertexFaceKernel);

	if (m_findSeparatingAxisEdgeEdgeKernel)
		clReleaseKernel(m_findSeparatingAxisEdgeEdgeKernel);

	if (m_findSeparatingAxisUnitSphereKernel)
		clReleaseKernel(m_findSeparatingAxisUnitSphereKernel);

	if (m_mprPenetrationKernel)
		clReleaseKernel(m_mprPenetrationKernel);


	if (m_findSeparatingAxisKernel)
		clReleaseKernel(m_findSeparatingAxisKernel);

    if (m_findConcaveSeparatingAxisVertexFaceKernel)
        clReleaseKernel(m_findConcaveSeparatingAxisVertexFaceKernel);

    
    if (m_findConcaveSeparatingAxisEdgeEdgeKernel)
        clReleaseKernel(m_findConcaveSeparatingAxisEdgeEdgeKernel);
    
	if (m_findConcaveSeparatingAxisKernel)
		clReleaseKernel(m_findConcaveSeparatingAxisKernel);

	if (m_findCompoundPairsKernel)
		clReleaseKernel(m_findCompoundPairsKernel);

	if (m_processCompoundPairsKernel)
		clReleaseKernel(m_processCompoundPairsKernel);
    
    if (m_findClippingFacesKernel)
        clReleaseKernel(m_findClippingFacesKernel);
   
    if (m_clipFacesAndFindContacts)
        clReleaseKernel(m_clipFacesAndFindContacts);
    if (m_newContactReductionKernel)
        clReleaseKernel(m_newContactReductionKernel);
	if (m_primitiveContactsKernel)
		clReleaseKernel(m_primitiveContactsKernel);
    
	if (m_findConcaveSphereContactsKernel)
		clReleaseKernel(m_findConcaveSphereContactsKernel);

	if (m_processCompoundPairsPrimitivesKernel)
		clReleaseKernel(m_processCompoundPairsPrimitivesKernel);

	if (m_clipHullHullKernel)
		clReleaseKernel(m_clipHullHullKernel);
	if (m_clipCompoundsHullHullKernel)
		clReleaseKernel(m_clipCompoundsHullHullKernel);

	if (m_clipHullHullConcaveConvexKernel)
		clReleaseKernel(m_clipHullHullConcaveConvexKernel);
//	if (m_extractManifoldAndAddContactKernel)
	//	clReleaseKernel(m_extractManifoldAndAddContactKernel);

	if (m_bvhTraversalKernel)
		clReleaseKernel(m_bvhTraversalKernel);

}







int numAabbChecks = 0;
int maxNumAabbChecks = 0;
int maxDepth = 0;



																
																
void GpuSatCollision::computeConvexConvexContactsGPUSAT( b3OpenCLArray<b3Int4>* pairs, int nPairs,
			const b3OpenCLArray<b3RigidBodyData>* bodyBuf,
			b3OpenCLArray<b3Contact4>* contactOut, int& nContacts,
			const b3OpenCLArray<b3Contact4>* oldContacts,
			int maxContactCapacity,
			int compoundPairCapacity,
			const b3OpenCLArray<b3ConvexPolyhedronData>& convexData,
			const b3OpenCLArray<b3Vector3>& gpuVertices,
			const b3OpenCLArray<b3Vector3>& gpuUniqueEdges,
			const b3OpenCLArray<b3GpuFace>& gpuFaces,
			const b3OpenCLArray<int>& gpuIndices,
			const b3OpenCLArray<b3Collidable>& gpuCollidables,
			const b3OpenCLArray<b3GpuChildShape>& gpuChildShapes,

			const b3OpenCLArray<b3Aabb>& clAabbsWorldSpace,
			const b3OpenCLArray<b3Aabb>& clAabbsLocalSpace,

            b3OpenCLArray<b3Vector3>& worldVertsB1GPU,
            b3OpenCLArray<b3Int4>& clippingFacesOutGPU,
            b3OpenCLArray<b3Vector3>& worldNormalsAGPU,
            b3OpenCLArray<b3Vector3>& worldVertsA1GPU,
            b3OpenCLArray<b3Vector3>& worldVertsB2GPU,    
			b3AlignedObjectArray<class b3OptimizedBvh*>& bvhDataUnused,
			b3OpenCLArray<b3QuantizedBvhNode>*	treeNodesGPU,
			b3OpenCLArray<b3BvhSubtreeInfo>*	subTreesGPU,
			b3OpenCLArray<b3BvhInfo>*	bvhInfo,

			int numObjects,
			int maxTriConvexPairCapacity,
			b3OpenCLArray<b3Int4>& triangleConvexPairsOut,
			int& numTriConvexPairsOut
			)
{
	myframecount++;

	if (!nPairs)
		return;

	{
		if (nPairs)
		{
			m_totalContactsOut.copyFromHostPointer(&nContacts,1,0,true);

			B3_PROFILE("primitiveContactsKernel");
			b3BufferInfoCL bInfo[] = {
				b3BufferInfoCL( pairs->getBufferCL(), true ), 
				b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
				b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
				b3BufferInfoCL( convexData.getBufferCL(),true),
				b3BufferInfoCL( gpuVertices.getBufferCL(),true),
				b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
				b3BufferInfoCL( gpuFaces.getBufferCL(),true),
				b3BufferInfoCL( gpuIndices.getBufferCL(),true),
				b3BufferInfoCL( contactOut->getBufferCL()),
				b3BufferInfoCL( m_totalContactsOut.getBufferCL())	
			};
			
			b3LauncherCL launcher(m_queue, m_primitiveContactsKernel,"m_primitiveContactsKernel");
			launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst( nPairs  );
			launcher.setConst(maxContactCapacity);
			int num = nPairs;
			launcher.launch1D( num);
			clFinish(m_queue);
		
			nContacts = m_totalContactsOut.at(0);
			contactOut->resize(nContacts);
		}
	}


	
	B3_PROFILE("computeConvexConvexContactsGPUSAT");
   // printf("nContacts = %d\n",nContacts);
    
	
	m_sepNormals.resize(nPairs);
	m_hasSeparatingNormals.resize(nPairs);
	
	int concaveCapacity=maxTriConvexPairCapacity;
	m_concaveSepNormals.resize(concaveCapacity);
	m_concaveHasSeparatingNormals.resize(concaveCapacity);
	m_numConcavePairsOut.resize(0);
	m_numConcavePairsOut.push_back(0);

	
	m_gpuCompoundPairs.resize(compoundPairCapacity);

	m_gpuCompoundSepNormals.resize(compoundPairCapacity);
	
	
	m_gpuHasCompoundSepNormals.resize(compoundPairCapacity);
	
	m_numCompoundPairsOut.resize(0);
	m_numCompoundPairsOut.push_back(0);

	int numCompoundPairs = 0;

	int numConcavePairs =0;

	{
		clFinish(m_queue);
		if (findSeparatingAxisOnGpu)
		{
			m_dmins.resize(nPairs);
			if (splitSearchSepAxisConvex)
			{
					

				if (useMprGpu)
				{
					nContacts = m_totalContactsOut.at(0);
					{
						B3_PROFILE("mprPenetrationKernel");
						b3BufferInfoCL bInfo[] = { 
							b3BufferInfoCL( pairs->getBufferCL(), true ), 
							b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
							b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
							b3BufferInfoCL( convexData.getBufferCL(),true),
							b3BufferInfoCL( gpuVertices.getBufferCL(),true),
							b3BufferInfoCL( m_sepNormals.getBufferCL()),
							b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL()),
							b3BufferInfoCL( contactOut->getBufferCL()),
							b3BufferInfoCL( m_totalContactsOut.getBufferCL())
						};

						b3LauncherCL launcher(m_queue, m_mprPenetrationKernel,"mprPenetrationKernel");
						launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );

						launcher.setConst(maxContactCapacity);
						launcher.setConst( nPairs  );

						int num = nPairs;
						launcher.launch1D( num);
						clFinish(m_queue);
						/*
						b3AlignedObjectArray<int>hostHasSepAxis;
						m_hasSeparatingNormals.copyToHost(hostHasSepAxis);
						b3AlignedObjectArray<b3Vector3>hostSepAxis;
						m_sepNormals.copyToHost(hostSepAxis);
						*/
						nContacts = m_totalContactsOut.at(0);
						contactOut->resize(nContacts);
					//	printf("nContacts (after mprPenetrationKernel) = %d\n",nContacts);
						if (nContacts>maxContactCapacity)
						{
                
							b3Error("Error: contacts exceeds capacity (%d/%d)\n", nContacts, maxContactCapacity);
							nContacts = maxContactCapacity;
						}

					}
				}
				
				if (1)
				{

					if (1)
					{
					{
						B3_PROFILE("findSeparatingAxisVertexFaceKernel");
						b3BufferInfoCL bInfo[] = { 
							b3BufferInfoCL( pairs->getBufferCL(), true ), 
							b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
							b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
							b3BufferInfoCL( convexData.getBufferCL(),true),
							b3BufferInfoCL( gpuVertices.getBufferCL(),true),
							b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
							b3BufferInfoCL( gpuFaces.getBufferCL(),true),
							b3BufferInfoCL( gpuIndices.getBufferCL(),true),
							b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
							b3BufferInfoCL( m_sepNormals.getBufferCL()),
							b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL()),
							b3BufferInfoCL( m_dmins.getBufferCL())
						};

						b3LauncherCL launcher(m_queue, m_findSeparatingAxisVertexFaceKernel,"findSeparatingAxisVertexFaceKernel");
						launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
						launcher.setConst( nPairs  );

						int num = nPairs;
						launcher.launch1D( num);
						clFinish(m_queue);
					}


					int numDirections = sizeof(unitSphere162)/sizeof(b3Vector3);
					
					{
						B3_PROFILE("findSeparatingAxisEdgeEdgeKernel");
						b3BufferInfoCL bInfo[] = { 
							b3BufferInfoCL( pairs->getBufferCL(), true ), 
							b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
							b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
							b3BufferInfoCL( convexData.getBufferCL(),true),
							b3BufferInfoCL( gpuVertices.getBufferCL(),true),
							b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
							b3BufferInfoCL( gpuFaces.getBufferCL(),true),
							b3BufferInfoCL( gpuIndices.getBufferCL(),true),
							b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
							b3BufferInfoCL( m_sepNormals.getBufferCL()),
							b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL()),
							b3BufferInfoCL( m_dmins.getBufferCL()),
							b3BufferInfoCL( m_unitSphereDirections.getBufferCL(),true)

						};

						b3LauncherCL launcher(m_queue, m_findSeparatingAxisEdgeEdgeKernel,"findSeparatingAxisEdgeEdgeKernel");
						launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
						launcher.setConst( numDirections);
						launcher.setConst( nPairs  );
						int num = nPairs;
						launcher.launch1D( num);
						clFinish(m_queue);

					}
					}
					if (useMprGpu)
					{
						B3_PROFILE("findSeparatingAxisUnitSphereKernel");
						b3BufferInfoCL bInfo[] = { 
								b3BufferInfoCL( pairs->getBufferCL(), true ), 
								b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
								b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
								b3BufferInfoCL( convexData.getBufferCL(),true),
								b3BufferInfoCL( gpuVertices.getBufferCL(),true),
								b3BufferInfoCL( m_unitSphereDirections.getBufferCL(),true),
								b3BufferInfoCL( m_sepNormals.getBufferCL()),
								b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL()),
								b3BufferInfoCL( m_dmins.getBufferCL())
						};

						b3LauncherCL launcher(m_queue, m_findSeparatingAxisUnitSphereKernel,"findSeparatingAxisUnitSphereKernel");
						launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
						int numDirections = sizeof(unitSphere162)/sizeof(b3Vector3);
						launcher.setConst( numDirections);

						launcher.setConst( nPairs  );
                                                
						int num = nPairs;
						launcher.launch1D( num);
						clFinish(m_queue);
					}
			}
				

			} else
			{
				B3_PROFILE("findSeparatingAxisKernel");
				b3BufferInfoCL bInfo[] = { 
					b3BufferInfoCL( pairs->getBufferCL(), true ), 
					b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
					b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
					b3BufferInfoCL( convexData.getBufferCL(),true),
					b3BufferInfoCL( gpuVertices.getBufferCL(),true),
					b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
					b3BufferInfoCL( gpuFaces.getBufferCL(),true),
					b3BufferInfoCL( gpuIndices.getBufferCL(),true),
					b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
					b3BufferInfoCL( m_sepNormals.getBufferCL()),
					b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL())
				};

				b3LauncherCL launcher(m_queue, m_findSeparatingAxisKernel,"m_findSeparatingAxisKernel");
				launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst( nPairs  );

				int num = nPairs;
				launcher.launch1D( num);
				clFinish(m_queue);
			}
			
			
		}
        
        
        numCompoundPairs = m_numCompoundPairsOut.at(0);
        bool useGpuFindCompoundPairs=true;
        if (useGpuFindCompoundPairs)
        {
            B3_PROFILE("findCompoundPairsKernel");
            b3BufferInfoCL bInfo[] = 
            { 
                b3BufferInfoCL( pairs->getBufferCL(), true ), 
                b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
                b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
                b3BufferInfoCL( convexData.getBufferCL(),true),
                b3BufferInfoCL( gpuVertices.getBufferCL(),true),
                b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
                b3BufferInfoCL( gpuFaces.getBufferCL(),true),
                b3BufferInfoCL( gpuIndices.getBufferCL(),true),
                b3BufferInfoCL( clAabbsLocalSpace.getBufferCL(),true),
                b3BufferInfoCL( gpuChildShapes.getBufferCL(),true),
                b3BufferInfoCL( m_gpuCompoundPairs.getBufferCL()),
                b3BufferInfoCL( m_numCompoundPairsOut.getBufferCL()),
                b3BufferInfoCL(subTreesGPU->getBufferCL()),
                b3BufferInfoCL(treeNodesGPU->getBufferCL()),
                b3BufferInfoCL(bvhInfo->getBufferCL())
            };

            b3LauncherCL launcher(m_queue, m_findCompoundPairsKernel,"m_findCompoundPairsKernel");
            launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
            launcher.setConst( nPairs  );
            launcher.setConst( compoundPairCapacity);

            int num = nPairs;
            launcher.launch1D( num);
            clFinish(m_queue);

            numCompoundPairs = m_numCompoundPairsOut.at(0);
            //printf("numCompoundPairs =%d\n",numCompoundPairs );
            if (numCompoundPairs)
            {
                //printf("numCompoundPairs=%d\n",numCompoundPairs);
            }
            

        }
		if (numCompoundPairs)
		{
			printf("numCompoundPairs=%d\n",numCompoundPairs);
		}

        if (numCompoundPairs > compoundPairCapacity)
        {
            b3Error("Exceeded compound pair capacity (%d/%d)\n", numCompoundPairs,  compoundPairCapacity);
            numCompoundPairs = compoundPairCapacity;
        }

        

        m_gpuCompoundPairs.resize(numCompoundPairs);
        m_gpuHasCompoundSepNormals.resize(numCompoundPairs);
        m_gpuCompoundSepNormals.resize(numCompoundPairs);
        

        if (numCompoundPairs)
        {
            B3_PROFILE("processCompoundPairsPrimitivesKernel");
            b3BufferInfoCL bInfo[] = 
            { 
                b3BufferInfoCL( m_gpuCompoundPairs.getBufferCL(), true ), 
                b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
                b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
                b3BufferInfoCL( convexData.getBufferCL(),true),
                b3BufferInfoCL( gpuVertices.getBufferCL(),true),
                b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
                b3BufferInfoCL( gpuFaces.getBufferCL(),true),
                b3BufferInfoCL( gpuIndices.getBufferCL(),true),
                b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
                b3BufferInfoCL( gpuChildShapes.getBufferCL(),true),
                b3BufferInfoCL( contactOut->getBufferCL()),
                b3BufferInfoCL( m_totalContactsOut.getBufferCL())	
            };

            b3LauncherCL launcher(m_queue, m_processCompoundPairsPrimitivesKernel,"m_processCompoundPairsPrimitivesKernel");
            launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
            launcher.setConst( numCompoundPairs  );
            launcher.setConst(maxContactCapacity);

            int num = numCompoundPairs;
            launcher.launch1D( num);
            clFinish(m_queue);
            nContacts = m_totalContactsOut.at(0);
            //printf("nContacts (after processCompoundPairsPrimitivesKernel) = %d\n",nContacts);
            if (nContacts>maxContactCapacity)
            {
                
                b3Error("Error: contacts exceeds capacity (%d/%d)\n", nContacts, maxContactCapacity);
                nContacts = maxContactCapacity;
            }
        }
        

        if (numCompoundPairs)
        {
            B3_PROFILE("processCompoundPairsKernel");
            b3BufferInfoCL bInfo[] = 
            { 
                b3BufferInfoCL( m_gpuCompoundPairs.getBufferCL(), true ), 
                b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
                b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
                b3BufferInfoCL( convexData.getBufferCL(),true),
                b3BufferInfoCL( gpuVertices.getBufferCL(),true),
                b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
                b3BufferInfoCL( gpuFaces.getBufferCL(),true),
                b3BufferInfoCL( gpuIndices.getBufferCL(),true),
                b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
                b3BufferInfoCL( gpuChildShapes.getBufferCL(),true),
                b3BufferInfoCL( m_gpuCompoundSepNormals.getBufferCL()),
                b3BufferInfoCL( m_gpuHasCompoundSepNormals.getBufferCL())
            };

            b3LauncherCL launcher(m_queue, m_processCompoundPairsKernel,"m_processCompoundPairsKernel");
            launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
            launcher.setConst( numCompoundPairs  );

            int num = numCompoundPairs;
            launcher.launch1D( num);
            clFinish(m_queue);
        
        }


        //printf("numConcave  = %d\n",numConcave);

    

//		printf("hostNormals.size()=%d\n",hostNormals.size());
		//int numPairs = pairCount.at(0);
		
		
		
	}
	int vertexFaceCapacity = 64;


		
	{
		//now perform the tree query on GPU
			
					
				
				
		if (treeNodesGPU->size() && treeNodesGPU->size())
		{
			if (bvhTraversalKernelGPU)
			{
						
				B3_PROFILE("m_bvhTraversalKernel");
						
						
				numConcavePairs = m_numConcavePairsOut.at(0);
						
				b3LauncherCL launcher(m_queue, m_bvhTraversalKernel,"m_bvhTraversalKernel");
				launcher.setBuffer( pairs->getBufferCL());
				launcher.setBuffer(  bodyBuf->getBufferCL());
				launcher.setBuffer( gpuCollidables.getBufferCL());
				launcher.setBuffer( clAabbsWorldSpace.getBufferCL());
				launcher.setBuffer( triangleConvexPairsOut.getBufferCL());
				launcher.setBuffer( m_numConcavePairsOut.getBufferCL());
				launcher.setBuffer( subTreesGPU->getBufferCL());
				launcher.setBuffer( treeNodesGPU->getBufferCL());
				launcher.setBuffer( bvhInfo->getBufferCL());
						
				launcher.setConst( nPairs  );
				launcher.setConst( maxTriConvexPairCapacity);
				int num = nPairs;
				launcher.launch1D( num);
				clFinish(m_queue);
				numConcavePairs = m_numConcavePairsOut.at(0);
			} 

				//printf("numConcavePairs=%d (max = %d\n",numConcavePairs,maxTriConvexPairCapacity);
						
			if (numConcavePairs > maxTriConvexPairCapacity)
			{
				static int exceeded_maxTriConvexPairCapacity_count = 0;
				b3Error("Exceeded the maxTriConvexPairCapacity (found %d but max is %d, it happened %d times)\n",
					numConcavePairs,maxTriConvexPairCapacity,exceeded_maxTriConvexPairCapacity_count++);
				numConcavePairs = maxTriConvexPairCapacity;
			}
			triangleConvexPairsOut.resize(numConcavePairs);
	
			if (numConcavePairs)
			{
				clippingFacesOutGPU.resize(numConcavePairs);
				worldNormalsAGPU.resize(numConcavePairs);
				worldVertsA1GPU.resize(vertexFaceCapacity*(numConcavePairs));
				worldVertsB1GPU.resize(vertexFaceCapacity*(numConcavePairs));


				if (findConcaveSeparatingAxisKernelGPU)
				{


					//now perform a SAT test for each triangle-convex element (stored in triangleConvexPairsOut)
                    if (splitSearchSepAxisConcave)
                    {
                        //printf("numConcavePairs = %d\n",numConcavePairs);
                        m_dmins.resize(numConcavePairs);
                        {
                            B3_PROFILE("findConcaveSeparatingAxisVertexFaceKernel");
                            b3BufferInfoCL bInfo[] = {
                                b3BufferInfoCL( triangleConvexPairsOut.getBufferCL() ),
                                b3BufferInfoCL( bodyBuf->getBufferCL(),true),
                                b3BufferInfoCL( gpuCollidables.getBufferCL(),true),
                                b3BufferInfoCL( convexData.getBufferCL(),true),
                                b3BufferInfoCL( gpuVertices.getBufferCL(),true),
                                b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
                                b3BufferInfoCL( gpuFaces.getBufferCL(),true),
                                b3BufferInfoCL( gpuIndices.getBufferCL(),true),
                                b3BufferInfoCL( gpuChildShapes.getBufferCL(),true),
                                b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
                                b3BufferInfoCL( m_concaveSepNormals.getBufferCL()),
                                b3BufferInfoCL( m_concaveHasSeparatingNormals.getBufferCL()),
                                b3BufferInfoCL( clippingFacesOutGPU.getBufferCL()),
                                b3BufferInfoCL( worldVertsA1GPU.getBufferCL()),
                                b3BufferInfoCL(worldNormalsAGPU.getBufferCL()),
                                b3BufferInfoCL(worldVertsB1GPU.getBufferCL()),
                                b3BufferInfoCL(m_dmins.getBufferCL())
                            };
                            
                            b3LauncherCL launcher(m_queue, m_findConcaveSeparatingAxisVertexFaceKernel,"m_findConcaveSeparatingAxisVertexFaceKernel");
                            launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
                            launcher.setConst(vertexFaceCapacity);
                            launcher.setConst( numConcavePairs  );
                            
                            int num = numConcavePairs;
                            launcher.launch1D( num);
                            clFinish(m_queue);

                            
                        }
//                        numConcavePairs = 0;
                        if (1)
                        {
                            B3_PROFILE("findConcaveSeparatingAxisEdgeEdgeKernel");
                            b3BufferInfoCL bInfo[] = {
                                b3BufferInfoCL( triangleConvexPairsOut.getBufferCL() ),
                                b3BufferInfoCL( bodyBuf->getBufferCL(),true),
                                b3BufferInfoCL( gpuCollidables.getBufferCL(),true),
                                b3BufferInfoCL( convexData.getBufferCL(),true),
                                b3BufferInfoCL( gpuVertices.getBufferCL(),true),
                                b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
                                b3BufferInfoCL( gpuFaces.getBufferCL(),true),
                                b3BufferInfoCL( gpuIndices.getBufferCL(),true),
                                b3BufferInfoCL( gpuChildShapes.getBufferCL(),true),
                                b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
                                b3BufferInfoCL( m_concaveSepNormals.getBufferCL()),
                                b3BufferInfoCL( m_concaveHasSeparatingNormals.getBufferCL()),
                                b3BufferInfoCL( clippingFacesOutGPU.getBufferCL()),
                                b3BufferInfoCL( worldVertsA1GPU.getBufferCL()),
                                b3BufferInfoCL(worldNormalsAGPU.getBufferCL()),
                                b3BufferInfoCL(worldVertsB1GPU.getBufferCL()),
                                b3BufferInfoCL(m_dmins.getBufferCL())
                            };
                            
                            b3LauncherCL launcher(m_queue, m_findConcaveSeparatingAxisEdgeEdgeKernel,"m_findConcaveSeparatingAxisEdgeEdgeKernel");
                            launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
                            launcher.setConst(vertexFaceCapacity);
                            launcher.setConst( numConcavePairs  );
                            
                            int num = numConcavePairs;
                            launcher.launch1D( num);
                            clFinish(m_queue);
                        }
                      
                        
                        // numConcavePairs = 0;
                        
                        
                        
                        
                        
                        
                    } else
                    {
                        B3_PROFILE("findConcaveSeparatingAxisKernel");
                        b3BufferInfoCL bInfo[] = { 
                            b3BufferInfoCL( triangleConvexPairsOut.getBufferCL() ), 
                            b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
                            b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
                            b3BufferInfoCL( convexData.getBufferCL(),true),
                            b3BufferInfoCL( gpuVertices.getBufferCL(),true),
                            b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
                            b3BufferInfoCL( gpuFaces.getBufferCL(),true),
                            b3BufferInfoCL( gpuIndices.getBufferCL(),true),
                            b3BufferInfoCL( gpuChildShapes.getBufferCL(),true),
                            b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
                            b3BufferInfoCL( m_concaveSepNormals.getBufferCL()),
                            b3BufferInfoCL( m_concaveHasSeparatingNormals.getBufferCL()),
                            b3BufferInfoCL( clippingFacesOutGPU.getBufferCL()),
                            b3BufferInfoCL( worldVertsA1GPU.getBufferCL()),
                            b3BufferInfoCL(worldNormalsAGPU.getBufferCL()),
                            b3BufferInfoCL(worldVertsB1GPU.getBufferCL())
                        };

                        b3LauncherCL launcher(m_queue, m_findConcaveSeparatingAxisKernel,"m_findConcaveSeparatingAxisKernel");
                        launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
                        launcher.setConst(vertexFaceCapacity);
                        launcher.setConst( numConcavePairs  );

                        int num = numConcavePairs;
                        launcher.launch1D( num);
                        clFinish(m_queue);
                    }
                    
                    
				}

			}
		}
		
		
	}

	if (numConcavePairs)
	{
			if (numConcavePairs)
		{
			B3_PROFILE("findConcaveSphereContactsKernel");
				nContacts = m_totalContactsOut.at(0);
//				printf("nContacts1 = %d\n",nContacts);
			b3BufferInfoCL bInfo[] = { 
				b3BufferInfoCL( triangleConvexPairsOut.getBufferCL() ), 
				b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
				b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
				b3BufferInfoCL( convexData.getBufferCL(),true),
				b3BufferInfoCL( gpuVertices.getBufferCL(),true),
				b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
				b3BufferInfoCL( gpuFaces.getBufferCL(),true),
				b3BufferInfoCL( gpuIndices.getBufferCL(),true),
				b3BufferInfoCL( clAabbsWorldSpace.getBufferCL(),true),
				b3BufferInfoCL( contactOut->getBufferCL()),
				b3BufferInfoCL( m_totalContactsOut.getBufferCL())
			};

			b3LauncherCL launcher(m_queue, m_findConcaveSphereContactsKernel,"m_findConcaveSphereContactsKernel");
			launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );

			launcher.setConst( numConcavePairs  );
			launcher.setConst(maxContactCapacity);

			int num = numConcavePairs;
			launcher.launch1D( num);
			clFinish(m_queue);
			nContacts = m_totalContactsOut.at(0);
			
			if (nContacts >= maxContactCapacity)
			{
				b3Error("Error: contacts exceeds capacity (%d/%d)\n", nContacts, maxContactCapacity);
				nContacts = maxContactCapacity;
			}
		}
		
	}



#ifdef __APPLE__
	bool contactClippingOnGpu = true;
#else
	bool contactClippingOnGpu = true;
#endif

	if (contactClippingOnGpu)
	{
		m_totalContactsOut.copyFromHostPointer(&nContacts,1,0,true);
//		printf("nContacts3 = %d\n",nContacts);


		//B3_PROFILE("clipHullHullKernel");

		bool breakupConcaveConvexKernel = true;

#ifdef __APPLE__
		//actually, some Apple OpenCL platform/device combinations work fine...
		breakupConcaveConvexKernel = true;
#endif
		//concave-convex contact clipping
		if (numConcavePairs)
		{
			//			printf("numConcavePairs = %d\n", numConcavePairs);
			//		nContacts = m_totalContactsOut.at(0);
			//	printf("nContacts before = %d\n", nContacts);

			if (breakupConcaveConvexKernel)
			{

				worldVertsB2GPU.resize(vertexFaceCapacity*numConcavePairs);


				
				{

					if (1)
					{



						B3_PROFILE("clipFacesAndFindContacts");
						//nContacts = m_totalContactsOut.at(0);
						//int h = m_hasSeparatingNormals.at(0);
						//int4 p = clippingFacesOutGPU.at(0);
						b3BufferInfoCL bInfo[] = {
							b3BufferInfoCL( m_concaveSepNormals.getBufferCL()),
							b3BufferInfoCL( m_concaveHasSeparatingNormals.getBufferCL()),
							b3BufferInfoCL( clippingFacesOutGPU.getBufferCL()),
							b3BufferInfoCL( worldVertsA1GPU.getBufferCL()),
							b3BufferInfoCL( worldNormalsAGPU.getBufferCL()),
							b3BufferInfoCL( worldVertsB1GPU.getBufferCL()),
							b3BufferInfoCL( worldVertsB2GPU.getBufferCL())
						};
						b3LauncherCL launcher(m_queue, m_clipFacesAndFindContacts,"m_clipFacesAndFindContacts");
						launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
						launcher.setConst(vertexFaceCapacity);

						launcher.setConst( numConcavePairs  );
						int debugMode = 0;
						launcher.setConst( debugMode);
						int num = numConcavePairs;
						launcher.launch1D( num);
						clFinish(m_queue);
					}
				}
				//contactReduction
				{
					int newContactCapacity=nContacts+numConcavePairs; 
					contactOut->reserve(newContactCapacity);
					if (reduceConcaveContactsOnGPU)
					{
//						printf("newReservation = %d\n",newReservation);
						{
							B3_PROFILE("newContactReductionKernel");
							b3BufferInfoCL bInfo[] =
							{
								b3BufferInfoCL( triangleConvexPairsOut.getBufferCL(), true ),
								b3BufferInfoCL( bodyBuf->getBufferCL(),true),
								b3BufferInfoCL( m_concaveSepNormals.getBufferCL()),
								b3BufferInfoCL( m_concaveHasSeparatingNormals.getBufferCL()),
								b3BufferInfoCL( contactOut->getBufferCL()),
								b3BufferInfoCL( clippingFacesOutGPU.getBufferCL()),
								b3BufferInfoCL( worldVertsB2GPU.getBufferCL()),
								b3BufferInfoCL( m_totalContactsOut.getBufferCL())
							};

							b3LauncherCL launcher(m_queue, m_newContactReductionKernel,"m_newContactReductionKernel");
							launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
							launcher.setConst(vertexFaceCapacity);
							launcher.setConst(newContactCapacity);
							launcher.setConst( numConcavePairs  );
							int num = numConcavePairs;

							launcher.launch1D( num);
						}
						nContacts = m_totalContactsOut.at(0);
						contactOut->resize(nContacts);

						//printf("contactOut4 (after newContactReductionKernel) = %d\n",nContacts);
					}

				}
				//re-use?


			} else
			{
				B3_PROFILE("clipHullHullConcaveConvexKernel");
				nContacts = m_totalContactsOut.at(0);
				int newContactCapacity = contactOut->capacity();

				//printf("contactOut5 = %d\n",nContacts);
				b3BufferInfoCL bInfo[] = { 
					b3BufferInfoCL( triangleConvexPairsOut.getBufferCL(), true ), 
					b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
					b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
					b3BufferInfoCL( convexData.getBufferCL(),true),
					b3BufferInfoCL( gpuVertices.getBufferCL(),true),
					b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
					b3BufferInfoCL( gpuFaces.getBufferCL(),true),
					b3BufferInfoCL( gpuIndices.getBufferCL(),true),
					b3BufferInfoCL( gpuChildShapes.getBufferCL(),true),
					b3BufferInfoCL( m_concaveSepNormals.getBufferCL()),
					b3BufferInfoCL( contactOut->getBufferCL()),
					b3BufferInfoCL( m_totalContactsOut.getBufferCL())	
				};
				b3LauncherCL launcher(m_queue, m_clipHullHullConcaveConvexKernel,"m_clipHullHullConcaveConvexKernel");
				launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst(newContactCapacity);
				launcher.setConst( numConcavePairs  );
				int num = numConcavePairs;
				launcher.launch1D( num);
				clFinish(m_queue);
				nContacts = m_totalContactsOut.at(0);
				contactOut->resize(nContacts);
				//printf("contactOut6 = %d\n",nContacts);
				b3AlignedObjectArray<b3Contact4> cpuContacts;
				contactOut->copyToHost(cpuContacts);
			}
			//			printf("nContacts after = %d\n", nContacts);
		}//numConcavePairs



		//convex-convex contact clipping
		
		bool breakupKernel = false;

#ifdef __APPLE__
		breakupKernel = true;
#endif

#ifdef CHECK_ON_HOST
	bool computeConvexConvex = false;
#else
	bool computeConvexConvex = true;
#endif//CHECK_ON_HOST
		if (computeConvexConvex)
		{
			B3_PROFILE("clipHullHullKernel");
		if (breakupKernel)
		{




			worldVertsB1GPU.resize(vertexFaceCapacity*nPairs);
			clippingFacesOutGPU.resize(nPairs);
			worldNormalsAGPU.resize(nPairs);
			worldVertsA1GPU.resize(vertexFaceCapacity*nPairs);
			worldVertsB2GPU.resize(vertexFaceCapacity*nPairs);

			if (findConvexClippingFacesGPU)
			{
				B3_PROFILE("findClippingFacesKernel");
				b3BufferInfoCL bInfo[] = {
					b3BufferInfoCL( pairs->getBufferCL(), true ),
					b3BufferInfoCL( bodyBuf->getBufferCL(),true),
					b3BufferInfoCL( gpuCollidables.getBufferCL(),true),
					b3BufferInfoCL( convexData.getBufferCL(),true),
					b3BufferInfoCL( gpuVertices.getBufferCL(),true),
					b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
					b3BufferInfoCL( gpuFaces.getBufferCL(),true), 
					b3BufferInfoCL( gpuIndices.getBufferCL(),true),
					b3BufferInfoCL( m_sepNormals.getBufferCL()),
					b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL()),
					b3BufferInfoCL( clippingFacesOutGPU.getBufferCL()),
					b3BufferInfoCL( worldVertsA1GPU.getBufferCL()),
					b3BufferInfoCL( worldNormalsAGPU.getBufferCL()),
					b3BufferInfoCL( worldVertsB1GPU.getBufferCL())
				};

				b3LauncherCL launcher(m_queue, m_findClippingFacesKernel,"m_findClippingFacesKernel");
				launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst( vertexFaceCapacity);
				launcher.setConst( nPairs  );
				int num = nPairs;
				launcher.launch1D( num);
				clFinish(m_queue);

			} else
			{
				
				float minDist = -1e30f;
				float maxDist = 0.02f;

				b3AlignedObjectArray<b3ConvexPolyhedronData> hostConvexData;
				convexData.copyToHost(hostConvexData);
				b3AlignedObjectArray<b3Collidable> hostCollidables;
				gpuCollidables.copyToHost(hostCollidables);

				b3AlignedObjectArray<int> hostHasSepNormals;
				m_hasSeparatingNormals.copyToHost(hostHasSepNormals);
				b3AlignedObjectArray<b3Vector3> cpuSepNormals;
				m_sepNormals.copyToHost(cpuSepNormals);

				b3AlignedObjectArray<b3Int4> hostPairs;
				pairs->copyToHost(hostPairs);
				b3AlignedObjectArray<b3RigidBodyData> hostBodyBuf;
				bodyBuf->copyToHost(hostBodyBuf);


				//worldVertsB1GPU.resize(vertexFaceCapacity*nPairs);
				b3AlignedObjectArray<b3Vector3> worldVertsB1CPU;
				worldVertsB1GPU.copyToHost(worldVertsB1CPU);

				b3AlignedObjectArray<b3Int4> clippingFacesOutCPU;
				clippingFacesOutGPU.copyToHost(clippingFacesOutCPU);

				b3AlignedObjectArray<b3Vector3> worldNormalsACPU;
				worldNormalsACPU.resize(nPairs);

				b3AlignedObjectArray<b3Vector3> worldVertsA1CPU;
				worldVertsA1CPU.resize(worldVertsA1GPU.size());
			
			
				b3AlignedObjectArray<b3Vector3> hostVertices;
				gpuVertices.copyToHost(hostVertices);
				b3AlignedObjectArray<b3GpuFace> hostFaces;
				gpuFaces.copyToHost(hostFaces);
				b3AlignedObjectArray<int> hostIndices;
				gpuIndices.copyToHost(hostIndices);
				

				for (int i=0;i<nPairs;i++)
				{

					int bodyIndexA = hostPairs[i].x;
					int bodyIndexB = hostPairs[i].y;
			
					int collidableIndexA = hostBodyBuf[bodyIndexA].m_collidableIdx;
					int collidableIndexB = hostBodyBuf[bodyIndexB].m_collidableIdx;
			
					int shapeIndexA = hostCollidables[collidableIndexA].m_shapeIndex;
					int shapeIndexB = hostCollidables[collidableIndexB].m_shapeIndex;
			

					if (hostHasSepNormals[i])
					{
						b3FindClippingFaces(cpuSepNormals[i],
							&hostConvexData[shapeIndexA],
							&hostConvexData[shapeIndexB],
							hostBodyBuf[bodyIndexA].m_pos,hostBodyBuf[bodyIndexA].m_quat,
							hostBodyBuf[bodyIndexB].m_pos,hostBodyBuf[bodyIndexB].m_quat,
							&worldVertsA1CPU.at(0),&worldNormalsACPU.at(0),
							&worldVertsB1CPU.at(0),
							vertexFaceCapacity,minDist,maxDist,
							&hostVertices.at(0),&hostFaces.at(0),
							&hostIndices.at(0),
							&hostVertices.at(0),&hostFaces.at(0),
							&hostIndices.at(0),&clippingFacesOutCPU.at(0),i);
					}
				}

				clippingFacesOutGPU.copyFromHost(clippingFacesOutCPU);
				worldVertsA1GPU.copyFromHost(worldVertsA1CPU);
				worldNormalsAGPU.copyFromHost(worldNormalsACPU);
				worldVertsB1GPU.copyFromHost(worldVertsB1CPU);

			}





			///clip face B against face A, reduce contacts and append them to a global contact array
			if (1)
			{
			
				{
					B3_PROFILE("clipFacesAndFindContacts");
					//nContacts = m_totalContactsOut.at(0);
					//int h = m_hasSeparatingNormals.at(0);
					//int4 p = clippingFacesOutGPU.at(0);
					b3BufferInfoCL bInfo[] = {
						b3BufferInfoCL( m_sepNormals.getBufferCL()),
						b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL()),
						b3BufferInfoCL( clippingFacesOutGPU.getBufferCL()),
						b3BufferInfoCL( worldVertsA1GPU.getBufferCL()),
						b3BufferInfoCL( worldNormalsAGPU.getBufferCL()),
						b3BufferInfoCL( worldVertsB1GPU.getBufferCL()),
						b3BufferInfoCL( worldVertsB2GPU.getBufferCL())
					};

					b3LauncherCL launcher(m_queue, m_clipFacesAndFindContacts,"m_clipFacesAndFindContacts");
					launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
					launcher.setConst(vertexFaceCapacity);

					launcher.setConst( nPairs  );
					int debugMode = 0;
					launcher.setConst( debugMode);
					int num = nPairs;
					launcher.launch1D( num);
					clFinish(m_queue);
				} 

				{
					nContacts = m_totalContactsOut.at(0);
					//printf("nContacts = %d\n",nContacts);

					int newContactCapacity = nContacts+nPairs;
					contactOut->reserve(newContactCapacity);

					if (reduceConvexContactsOnGPU)
					{
						{
							B3_PROFILE("newContactReductionKernel");
							b3BufferInfoCL bInfo[] =
							{
								b3BufferInfoCL( pairs->getBufferCL(), true ),
								b3BufferInfoCL( bodyBuf->getBufferCL(),true),
								b3BufferInfoCL( m_sepNormals.getBufferCL()),
								b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL()),
								b3BufferInfoCL( contactOut->getBufferCL()),
								b3BufferInfoCL( clippingFacesOutGPU.getBufferCL()),
								b3BufferInfoCL( worldVertsB2GPU.getBufferCL()),
								b3BufferInfoCL( m_totalContactsOut.getBufferCL())
							};

							b3LauncherCL launcher(m_queue, m_newContactReductionKernel,"m_newContactReductionKernel");
							launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
							launcher.setConst(vertexFaceCapacity);
							launcher.setConst(newContactCapacity);
							launcher.setConst( nPairs  );
							int num = nPairs;

							launcher.launch1D( num);
						}
						nContacts = m_totalContactsOut.at(0);
						contactOut->resize(nContacts);
					} else
					{

						volatile int nGlobalContactsOut = nContacts;
						b3AlignedObjectArray<b3Int4> hostPairs;
						pairs->copyToHost(hostPairs);
						b3AlignedObjectArray<b3RigidBodyData> hostBodyBuf;
						bodyBuf->copyToHost(hostBodyBuf);
						b3AlignedObjectArray<b3Vector3> hostSepNormals;
						m_sepNormals.copyToHost(hostSepNormals);
						b3AlignedObjectArray<int> hostHasSepAxis;
						m_hasSeparatingNormals.copyToHost(hostHasSepAxis);
						b3AlignedObjectArray<b3Contact4> hostContactsOut;
						contactOut->copyToHost(hostContactsOut);
						hostContactsOut.resize(newContactCapacity);

						b3AlignedObjectArray<b3Int4> hostClippingFaces;
						clippingFacesOutGPU.copyToHost(hostClippingFaces);
						b3AlignedObjectArray<b3Vector3> worldVertsB2CPU;
						worldVertsB2GPU.copyToHost(worldVertsB2CPU);

						for (int i=0;i<nPairs;i++)
						{
							b3NewContactReductionKernel(&hostPairs.at(0),
								&hostBodyBuf.at(0),
								&hostSepNormals.at(0),
								&hostHasSepAxis.at(0),
								&hostContactsOut.at(0),
								&hostClippingFaces.at(0),
								&worldVertsB2CPU.at(0),
								&nGlobalContactsOut,
								vertexFaceCapacity,
								newContactCapacity,
								nPairs,
								i);
						}

						nContacts = nGlobalContactsOut;
						m_totalContactsOut.copyFromHostPointer(&nContacts,1,0,true);
						hostContactsOut.resize(nContacts);
						//printf("contactOut4 (after newContactReductionKernel) = %d\n",nContacts);
						contactOut->copyFromHost(hostContactsOut);
					}
					//                    b3Contact4 pt = contactOut->at(0);
					//                  printf("nContacts = %d\n",nContacts);
				}
			}
		}            
		else//breakupKernel
		{

			if (nPairs)
			{
				b3BufferInfoCL bInfo[] = {
					b3BufferInfoCL( pairs->getBufferCL(), true ), 
					b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
					b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
					b3BufferInfoCL( convexData.getBufferCL(),true),
					b3BufferInfoCL( gpuVertices.getBufferCL(),true),
					b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
					b3BufferInfoCL( gpuFaces.getBufferCL(),true),
					b3BufferInfoCL( gpuIndices.getBufferCL(),true),
					b3BufferInfoCL( m_sepNormals.getBufferCL()),
					b3BufferInfoCL( m_hasSeparatingNormals.getBufferCL()),
					b3BufferInfoCL( contactOut->getBufferCL()),
					b3BufferInfoCL( m_totalContactsOut.getBufferCL())	
				};
				b3LauncherCL launcher(m_queue, m_clipHullHullKernel,"m_clipHullHullKernel");
				launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
				launcher.setConst( nPairs  );
				launcher.setConst(maxContactCapacity);

				int num = nPairs;
				launcher.launch1D( num);
				clFinish(m_queue);

				nContacts = m_totalContactsOut.at(0);
				if (nContacts >= maxContactCapacity)
				{
					b3Error("Exceeded contact capacity (%d/%d)\n",nContacts,maxContactCapacity);
					nContacts = maxContactCapacity;
				}
				contactOut->resize(nContacts);
			}
		}


		int nCompoundsPairs = m_gpuCompoundPairs.size();

		if (nCompoundsPairs)
		{
			b3BufferInfoCL bInfo[] = {
				b3BufferInfoCL( m_gpuCompoundPairs.getBufferCL(), true ), 
				b3BufferInfoCL( bodyBuf->getBufferCL(),true), 
				b3BufferInfoCL( gpuCollidables.getBufferCL(),true), 
				b3BufferInfoCL( convexData.getBufferCL(),true),
				b3BufferInfoCL( gpuVertices.getBufferCL(),true),
				b3BufferInfoCL( gpuUniqueEdges.getBufferCL(),true),
				b3BufferInfoCL( gpuFaces.getBufferCL(),true),
				b3BufferInfoCL( gpuIndices.getBufferCL(),true),
				b3BufferInfoCL( gpuChildShapes.getBufferCL(),true),
				b3BufferInfoCL( m_gpuCompoundSepNormals.getBufferCL(),true),
				b3BufferInfoCL( m_gpuHasCompoundSepNormals.getBufferCL(),true),
				b3BufferInfoCL( contactOut->getBufferCL()),
				b3BufferInfoCL( m_totalContactsOut.getBufferCL())	
			};
			b3LauncherCL launcher(m_queue, m_clipCompoundsHullHullKernel,"m_clipCompoundsHullHullKernel");
			launcher.setBuffers( bInfo, sizeof(bInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst( nCompoundsPairs  );
			launcher.setConst(maxContactCapacity);

			int num = nCompoundsPairs;
			launcher.launch1D( num);
			clFinish(m_queue);

			nContacts = m_totalContactsOut.at(0);
			if (nContacts>maxContactCapacity)
			{

				b3Error("Error: contacts exceeds capacity (%d/%d)\n", nContacts, maxContactCapacity);
				nContacts = maxContactCapacity;
			}
			contactOut->resize(nContacts);
		}//if nCompoundsPairs
		}
	}//contactClippingOnGpu

}
