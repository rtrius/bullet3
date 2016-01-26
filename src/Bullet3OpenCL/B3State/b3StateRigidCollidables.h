/*
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef B3_STATE_RIGID_COLLIDABLES
#define B3_STATE_RIGID_COLLIDABLES

#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3ConvexPolyhedronData.h"
#include "Bullet3Collision/NarrowPhaseCollision/b3ConvexUtility.h"

#include "Bullet3OpenCL/BroadphaseCollision/b3SapAabb.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3OpenCL/NarrowphaseCollision/b3QuantizedBvh.h"
#include "Bullet3OpenCL/NarrowphaseCollision/b3BvhInfo.h"

#define B3_INVALID_COLLIDABLE_INDEX -1
#define B3_INVALID_SHAPE_INDEX -1

struct b3StateCollidablesLimits
{
	int m_maxConvexShapes;
	int	m_maxCompoundChildShapes;
	
	b3StateCollidablesLimits() :
		m_maxConvexShapes(32 * 1024),
		m_maxCompoundChildShapes(8192) {}
};

struct b3StateRigidCollidables
{
	b3StateCollidablesLimits m_limits;

	int m_numCollidables;

	//
	b3AlignedObjectArray<b3Collidable> m_collidablesCpu;
	b3AlignedObjectArray<b3SapAabb> m_localShapeAabbCpu;	//m_localShapeAabbCpu has same number of elements as m_collidablesCpu

	b3OpenCLArray<b3Collidable> m_collidablesGpu;
	b3OpenCLArray<b3SapAabb> m_localShapeAabbGpu;
	
	//Data used by convex, compound, and trimesh shapes
	b3AlignedObjectArray<b3ConvexPolyhedronData> m_convexPolyhedra;		//Each element references a set of faces, edges and vertices
	b3AlignedObjectArray<b3GpuFace> m_convexFaces;						//Each element references a set of indices
	b3AlignedObjectArray<b3Vector3> m_uniqueEdges;
	b3AlignedObjectArray<b3Vector3> m_convexVertices;
	b3AlignedObjectArray<int> m_convexIndices;

	b3OpenCLArray<b3ConvexPolyhedronData> m_convexPolyhedraGpu;
	b3OpenCLArray<b3GpuFace> m_convexFacesGpu;
	b3OpenCLArray<b3Vector3> m_uniqueEdgesGpu;
	b3OpenCLArray<b3Vector3> m_convexVerticesGpu;
	b3OpenCLArray<int> m_convexIndicesGpu;

	//Data for compound shapes and trimesh shapes
	b3AlignedObjectArray<b3BvhInfo> m_bvhInfoCpu;
	b3AlignedObjectArray<b3QuantizedBvhNode> m_treeNodesCpu;
	b3AlignedObjectArray<b3BvhSubtreeInfo> m_subTreesCpu;

	b3OpenCLArray<b3BvhInfo> m_bvhInfoGpu;
	b3OpenCLArray<b3QuantizedBvhNode> m_treeNodesGpu;
	b3OpenCLArray<b3BvhSubtreeInfo>	m_subTreesGpu;

	//Data for compound shapes only
	b3AlignedObjectArray<b3GpuChildShape> m_childShapes;
	b3OpenCLArray<b3GpuChildShape> m_childShapesGpu;

	//CPU only
	b3AlignedObjectArray<b3ConvexUtility*> m_convexData;
	b3AlignedObjectArray<class b3OptimizedBvh*> m_bvhData;
	b3AlignedObjectArray<class b3TriangleIndexVertexArray*> m_meshInterfaces;

	//
	b3StateRigidCollidables(cl_context context, cl_command_queue queue);
	virtual~b3StateRigidCollidables();

	int registerSphereShape(float radius);
	int registerPlaneShape(const b3Vector3& planeNormal, float planeConstant);

	int	registerConvexHullShape(const float* vertices, int strideInBytes, int numVertices, const float* scaling);
	int	registerConvexHullShape(b3ConvexUtility* utilPtr);

	int registerCompoundShape(b3AlignedObjectArray<b3GpuChildShape>* childShapes);

	int	registerConcaveMesh(b3AlignedObjectArray<b3Vector3>* vertices, b3AlignedObjectArray<int>* indices, const float* scaling);
	
	//
	void writeCollidablesToGpu();
	void reset();


protected:
	void resetCpuOnlyArrays();

	int allocateCollidable();
	int registerFace(const b3Vector3& faceNormal, float faceConstant);
	int registerConvexHullShapeInternal(class b3ConvexUtility* convexPtr, b3Collidable& col);
	int registerConcaveMeshShape(b3AlignedObjectArray<b3Vector3>* vertices, b3AlignedObjectArray<int>* indices, b3Collidable& col, const float* scaling);
};


#endif
