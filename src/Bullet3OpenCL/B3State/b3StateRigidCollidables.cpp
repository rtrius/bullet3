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
#include "b3StateRigidCollidables.h"

#include "Bullet3Geometry/b3AabbUtil.h"

#include "Bullet3OpenCL/NarrowphaseCollision/b3OptimizedBvh.h"
#include "Bullet3OpenCL/NarrowphaseCollision/b3TriangleIndexVertexArray.h"

b3StateRigidCollidables::b3StateRigidCollidables(cl_context context, cl_command_queue queue) :
	m_collidablesGpu(context, queue),
	m_localShapeAabbGpu(context, queue),
	m_convexPolyhedraGpu(context, queue),
	m_convexFacesGpu(context, queue),
	m_uniqueEdgesGpu(context, queue),
	m_convexVerticesGpu(context, queue),
	m_convexIndicesGpu(context, queue),
	m_subTreesGpu(context, queue),
	m_treeNodesGpu(context, queue),
	m_bvhInfoGpu(context, queue),
	m_childShapesGpu(context, queue)
{
	m_numCollidables = 0;
}
b3StateRigidCollidables::~b3StateRigidCollidables()
{
	resetCpuOnlyArrays();
}

int b3StateRigidCollidables::registerSphereShape(float radius)
{
	int collidableIndex = allocateCollidable();
	if (collidableIndex == B3_INVALID_COLLIDABLE_INDEX) return collidableIndex;

	b3Collidable& col = m_collidablesCpu[collidableIndex];
	col.m_shapeType = SHAPE_SPHERE;
	col.m_shapeIndex = B3_INVALID_SHAPE_INDEX;
	col.m_radius = radius;

	{
		b3SapAabb aabb;
		b3Vector3 myAabbMin = b3MakeVector3(-radius, -radius, -radius);
		b3Vector3 myAabbMax = b3MakeVector3(radius, radius, radius);

		aabb.m_min[0] = myAabbMin[0];
		aabb.m_min[1] = myAabbMin[1];
		aabb.m_min[2] = myAabbMin[2];
		aabb.m_minIndices[3] = 0;

		aabb.m_max[0] = myAabbMax[0];
		aabb.m_max[1] = myAabbMax[1];
		aabb.m_max[2] = myAabbMax[2];
		aabb.m_signedMaxIndices[3] = 0;

		m_localShapeAabbCpu.push_back(aabb);
	}

	return collidableIndex;
}



int	b3StateRigidCollidables::registerPlaneShape(const b3Vector3& planeNormal, float planeConstant)
{
	int collidableIndex = allocateCollidable();
	if (collidableIndex == B3_INVALID_COLLIDABLE_INDEX) return collidableIndex;


	b3Collidable& col = m_collidablesCpu[collidableIndex];
	col.m_shapeType = SHAPE_PLANE;
	col.m_shapeIndex = registerFace(planeNormal, planeConstant);
	col.m_radius = planeConstant;

	if (col.m_shapeIndex >= 0)
	{
		b3SapAabb aabb;
		aabb.m_min[0] = -1e30f;
		aabb.m_min[1] = -1e30f;
		aabb.m_min[2] = -1e30f;
		aabb.m_minIndices[3] = 0;

		aabb.m_max[0] = 1e30f;
		aabb.m_max[1] = 1e30f;
		aabb.m_max[2] = 1e30f;
		aabb.m_signedMaxIndices[3] = 0;

		m_localShapeAabbCpu.push_back(aabb);
	}

	return collidableIndex;
}



int b3StateRigidCollidables::registerConvexHullShape(const float* vertices, int strideInBytes, int numVertices, const float* scaling)
{
	b3AlignedObjectArray<b3Vector3> verts;

	unsigned char* vts = (unsigned char*)vertices;
	for (int i = 0; i<numVertices; i++)
	{
		float* vertex = (float*)&vts[i*strideInBytes];
		verts.push_back(b3MakeVector3(vertex[0] * scaling[0], vertex[1] * scaling[1], vertex[2] * scaling[2]));
	}

	b3ConvexUtility* utilPtr = new b3ConvexUtility();
	bool merge = true;
	if (numVertices)
	{
		utilPtr->initializePolyhedralFeatures(&verts[0], verts.size(), merge);
	}

	int collidableIndex = registerConvexHullShape(utilPtr);
	delete utilPtr;
	return collidableIndex;
}

int b3StateRigidCollidables::registerConvexHullShape(b3ConvexUtility* utilPtr)
{
	int collidableIndex = allocateCollidable();
	if (collidableIndex == B3_INVALID_COLLIDABLE_INDEX) return collidableIndex;

	b3Collidable& col = m_collidablesCpu[collidableIndex];
	col.m_shapeType = SHAPE_CONVEX_HULL;

	{
		b3Vector3 localCenter = b3MakeVector3(0, 0, 0);
		for (int i = 0; i<utilPtr->m_vertices.size(); i++)
			localCenter += utilPtr->m_vertices[i];
		localCenter *= (1.f / utilPtr->m_vertices.size());
		utilPtr->m_localCenter = localCenter;

		col.m_shapeIndex = registerConvexHullShapeInternal(utilPtr, col);
	}

	if (col.m_shapeIndex >= 0)
	{
		b3Vector3 myAabbMin = b3MakeVector3(1e30f, 1e30f, 1e30f);
		b3Vector3 myAabbMax = b3MakeVector3(-1e30f, -1e30f, -1e30f);

		for (int i = 0; i<utilPtr->m_vertices.size(); i++)
		{
			myAabbMin.setMin(utilPtr->m_vertices[i]);
			myAabbMax.setMax(utilPtr->m_vertices[i]);
		}

		b3SapAabb aabb;
		aabb.m_min[0] = myAabbMin[0];
		aabb.m_min[1] = myAabbMin[1];
		aabb.m_min[2] = myAabbMin[2];
		aabb.m_minIndices[3] = 0;

		aabb.m_max[0] = myAabbMax[0];
		aabb.m_max[1] = myAabbMax[1];
		aabb.m_max[2] = myAabbMax[2];
		aabb.m_signedMaxIndices[3] = 0;

		m_localShapeAabbCpu.push_back(aabb);
	}

	return collidableIndex;
}

int	b3StateRigidCollidables::registerCompoundShape(b3AlignedObjectArray<b3GpuChildShape>* childShapes)
{
	int collidableIndex = allocateCollidable();
	if (collidableIndex == B3_INVALID_COLLIDABLE_INDEX) return collidableIndex;

	b3Collidable& col = m_collidablesCpu[collidableIndex];
	col.m_shapeType = SHAPE_COMPOUND_OF_CONVEX_HULLS;
	col.m_shapeIndex = m_childShapes.size();
	col.m_compoundBvhIndex = m_bvhInfoCpu.size();
	col.m_numChildShapes = childShapes->size();

	{
		b3Assert(col.m_shapeIndex + childShapes->size()< m_limits.m_maxCompoundChildShapes);
		for (int i = 0; i<childShapes->size(); i++) m_childShapes.push_back( childShapes->at(i) );
	}
	
	b3SapAabb aabbLocalSpace;
	b3Vector3 myAabbMin = b3MakeVector3(1e30f, 1e30f, 1e30f);
	b3Vector3 myAabbMax = b3MakeVector3(-1e30f, -1e30f, -1e30f);

	b3AlignedObjectArray<b3Aabb> childLocalAabbs;
	childLocalAabbs.resize(childShapes->size());

	//compute local AABB of the compound of all children
	for (int i = 0; i<childShapes->size(); i++)
	{
		int childColIndex = childShapes->at(i).m_shapeIndex;
		b3Collidable& childCol = m_collidablesCpu[childColIndex];
		b3SapAabb aabbLoc = m_localShapeAabbCpu[childColIndex];

		b3Vector3 childLocalAabbMin = b3MakeVector3(aabbLoc.m_min[0], aabbLoc.m_min[1], aabbLoc.m_min[2]);
		b3Vector3 childLocalAabbMax = b3MakeVector3(aabbLoc.m_max[0], aabbLoc.m_max[1], aabbLoc.m_max[2]);
		b3Vector3 aMin, aMax;
		b3Scalar margin(0.f);
		b3Transform childTr;
		childTr.setIdentity();

		childTr.setOrigin(childShapes->at(i).m_childPosition);
		childTr.setRotation(b3Quaternion(childShapes->at(i).m_childOrientation));
		b3TransformAabb(childLocalAabbMin, childLocalAabbMax, margin, childTr, aMin, aMax);
		myAabbMin.setMin(aMin);
		myAabbMax.setMax(aMax);
		childLocalAabbs[i].m_min[0] = aMin[0];
		childLocalAabbs[i].m_min[1] = aMin[1];
		childLocalAabbs[i].m_min[2] = aMin[2];
		childLocalAabbs[i].m_min[3] = 0;
		childLocalAabbs[i].m_max[0] = aMax[0];
		childLocalAabbs[i].m_max[1] = aMax[1];
		childLocalAabbs[i].m_max[2] = aMax[2];
		childLocalAabbs[i].m_max[3] = 0;
	}

	aabbLocalSpace.m_min[0] = myAabbMin[0];
	aabbLocalSpace.m_min[1] = myAabbMin[1];
	aabbLocalSpace.m_min[2] = myAabbMin[2];
	aabbLocalSpace.m_minIndices[3] = 0;

	aabbLocalSpace.m_max[0] = myAabbMax[0];
	aabbLocalSpace.m_max[1] = myAabbMax[1];
	aabbLocalSpace.m_max[2] = myAabbMax[2];
	aabbLocalSpace.m_signedMaxIndices[3] = 0;

	m_localShapeAabbCpu.push_back(aabbLocalSpace);


	b3QuantizedBvh* bvh = new b3QuantizedBvh;
	bvh->setQuantizationValues(myAabbMin, myAabbMax);
	QuantizedNodeArray&	nodes = bvh->getLeafNodeArray();
	int numNodes = childShapes->size();

	for (int i = 0; i<numNodes; i++)
	{
		b3QuantizedBvhNode node;
		b3Vector3 aabbMin, aabbMax;
		aabbMin = (b3Vector3&)childLocalAabbs[i].m_min;
		aabbMax = (b3Vector3&)childLocalAabbs[i].m_max;

		bvh->quantize(&node.m_quantizedAabbMin[0], aabbMin, 0);
		bvh->quantize(&node.m_quantizedAabbMax[0], aabbMax, 1);
		int partId = 0;
		node.m_escapeIndexOrTriangleIndex = (partId << (31 - MAX_NUM_PARTS_IN_BITS)) | i;
		nodes.push_back(node);
	}
	bvh->buildInternal();	///buildInternal is expert use only: assumes that setQuantizationValues and LeafNodeArray are initialized

	int numSubTrees = bvh->getSubtreeInfoArray().size();

	b3BvhInfo bvhInfo;

	bvhInfo.m_aabbMin = bvh->m_bvhAabbMin;
	bvhInfo.m_aabbMax = bvh->m_bvhAabbMax;
	bvhInfo.m_quantization = bvh->m_bvhQuantization;
	bvhInfo.m_numNodes = numNodes;
	bvhInfo.m_numSubTrees = numSubTrees;
	bvhInfo.m_nodeOffset = m_treeNodesCpu.size();
	bvhInfo.m_subTreeOffset = m_subTreesCpu.size();

	int numNewNodes = bvh->getQuantizedNodeArray().size();

	for (int i = 0; i<numNewNodes - 1; i++)
	{
		if (bvh->getQuantizedNodeArray()[i].isLeafNode())
		{
			int orgIndex = bvh->getQuantizedNodeArray()[i].getTriangleIndex();

			b3Vector3 nodeMinVec = bvh->unQuantize(bvh->getQuantizedNodeArray()[i].m_quantizedAabbMin);
			b3Vector3 nodeMaxVec = bvh->unQuantize(bvh->getQuantizedNodeArray()[i].m_quantizedAabbMax);

			for (int c = 0; c<3; c++)
			{
				if (childLocalAabbs[orgIndex].m_min[c] < nodeMinVec[c])
				{
					printf("min org (%f) and new (%f) ? at i:%d,c:%d\n", childLocalAabbs[i].m_min[c], nodeMinVec[c], i, c);
				}
				if (childLocalAabbs[orgIndex].m_max[c] > nodeMaxVec[c])
				{
					printf("max org (%f) and new (%f) ? at i:%d,c:%d\n", childLocalAabbs[i].m_max[c], nodeMaxVec[c], i, c);
				}

			}
		}

	}

	m_bvhInfoCpu.push_back(bvhInfo);

	int numNewSubtrees = bvh->getSubtreeInfoArray().size();
	m_subTreesCpu.reserve(m_subTreesCpu.size() + numNewSubtrees);
	for (int i = 0; i<numNewSubtrees; i++)
	{
		m_subTreesCpu.push_back(bvh->getSubtreeInfoArray()[i]);
	}
	int numNewTreeNodes = bvh->getQuantizedNodeArray().size();

	for (int i = 0; i<numNewTreeNodes; i++)
	{
		m_treeNodesCpu.push_back(bvh->getQuantizedNodeArray()[i]);
	}

	return collidableIndex;
}


int b3StateRigidCollidables::registerConcaveMesh(b3AlignedObjectArray<b3Vector3>* vertices, b3AlignedObjectArray<int>* indices, const float* scaling1)
{
	int collidableIndex = allocateCollidable();
	if (collidableIndex == B3_INVALID_COLLIDABLE_INDEX) return collidableIndex;

	b3Vector3 scaling = b3MakeVector3(scaling1[0], scaling1[1], scaling1[2]);
	b3Collidable& col = m_collidablesCpu[collidableIndex];

	col.m_shapeType = SHAPE_CONCAVE_TRIMESH;
	col.m_shapeIndex = registerConcaveMeshShape(vertices, indices, col, scaling);
	col.m_bvhIndex = m_bvhInfoCpu.size();


	b3SapAabb aabb;
	b3Vector3 myAabbMin = b3MakeVector3(1e30f, 1e30f, 1e30f);
	b3Vector3 myAabbMax = b3MakeVector3(-1e30f, -1e30f, -1e30f);

	for (int i = 0; i<vertices->size(); i++)
	{
		b3Vector3 vtx(vertices->at(i)*scaling);
		myAabbMin.setMin(vtx);
		myAabbMax.setMax(vtx);
	}
	aabb.m_min[0] = myAabbMin[0];
	aabb.m_min[1] = myAabbMin[1];
	aabb.m_min[2] = myAabbMin[2];
	aabb.m_minIndices[3] = 0;

	aabb.m_max[0] = myAabbMax[0];
	aabb.m_max[1] = myAabbMax[1];
	aabb.m_max[2] = myAabbMax[2];
	aabb.m_signedMaxIndices[3] = 0;

	m_localShapeAabbCpu.push_back(aabb);

	b3OptimizedBvh* bvh = new b3OptimizedBvh();

	bool useQuantizedAabbCompression = true;
	b3TriangleIndexVertexArray* meshInterface = new b3TriangleIndexVertexArray();
	m_meshInterfaces.push_back(meshInterface);
	b3IndexedMesh mesh;
	mesh.m_numTriangles = indices->size() / 3;
	mesh.m_numVertices = vertices->size();
	mesh.m_vertexBase = (const unsigned char *)&vertices->at(0).x;
	mesh.m_vertexStride = sizeof(b3Vector3);
	mesh.m_triangleIndexStride = 3 * sizeof(int);// or sizeof(int)
	mesh.m_triangleIndexBase = (const unsigned char *)&indices->at(0);

	meshInterface->addIndexedMesh(mesh);
	bvh->build(meshInterface, useQuantizedAabbCompression, (b3Vector3&)aabb.m_min, (b3Vector3&)aabb.m_max);
	m_bvhData.push_back(bvh);
	int numNodes = bvh->getQuantizedNodeArray().size();
	int numSubTrees = bvh->getSubtreeInfoArray().size();

	b3BvhInfo bvhInfo;

	bvhInfo.m_aabbMin = bvh->m_bvhAabbMin;
	bvhInfo.m_aabbMax = bvh->m_bvhAabbMax;
	bvhInfo.m_quantization = bvh->m_bvhQuantization;
	bvhInfo.m_numNodes = numNodes;
	bvhInfo.m_numSubTrees = numSubTrees;
	bvhInfo.m_nodeOffset = m_treeNodesCpu.size();
	bvhInfo.m_subTreeOffset = m_subTreesCpu.size();

	m_bvhInfoCpu.push_back(bvhInfo);


	int numNewSubtrees = bvh->getSubtreeInfoArray().size();
	m_subTreesCpu.reserve(m_subTreesCpu.size() + numNewSubtrees);
	for (int i = 0; i<numNewSubtrees; i++)
	{
		m_subTreesCpu.push_back(bvh->getSubtreeInfoArray()[i]);
	}
	int numNewTreeNodes = bvh->getQuantizedNodeArray().size();

	for (int i = 0; i<numNewTreeNodes; i++)
	{
		m_treeNodesCpu.push_back(bvh->getQuantizedNodeArray()[i]);
	}

	return collidableIndex;
}


void	b3StateRigidCollidables::writeCollidablesToGpu()
{
	if (m_collidablesCpu.size()) m_collidablesGpu.copyFromHost(m_collidablesCpu);
	if (m_localShapeAabbCpu.size()) m_localShapeAabbGpu.copyFromHost(m_localShapeAabbCpu);

	m_convexPolyhedraGpu.copyFromHost(m_convexPolyhedra);
	m_convexFacesGpu.copyFromHost(m_convexFaces);
	m_uniqueEdgesGpu.copyFromHost(m_uniqueEdges);
	m_convexVerticesGpu.copyFromHost(m_convexVertices);
	m_convexIndicesGpu.copyFromHost(m_convexIndices);

	m_bvhInfoGpu.copyFromHost(m_bvhInfoCpu);
	m_treeNodesGpu.copyFromHost(m_treeNodesCpu);
	m_subTreesGpu.copyFromHost(m_subTreesCpu);

	m_childShapesGpu.copyFromHost(m_childShapes);
}


void b3StateRigidCollidables::reset()
{
	m_numCollidables = 0;

	m_collidablesCpu.resize(0);
	m_localShapeAabbCpu.resize(0);

	m_convexPolyhedra.resize(0);
	m_convexFaces.resize(0);
	m_uniqueEdges.resize(0);
	m_convexVertices.resize(0);
	m_convexIndices.resize(0);

	m_bvhInfoCpu.resize(0);
	m_treeNodesCpu.resize(0);
	m_subTreesCpu.resize(0);

	m_childShapes.resize(0);

	resetCpuOnlyArrays();
}

void b3StateRigidCollidables::resetCpuOnlyArrays()
{
	//m_convexData[] pointers are not allocated by b3StateRigidCollidables so do not delete it here
	for (int i = 0; i<m_bvhData.size(); i++) delete m_bvhData[i];
	for (int i = 0; i<m_meshInterfaces.size(); i++) delete m_meshInterfaces[i];

	m_convexData.resize(0);
	m_bvhData.resize(0);
	m_meshInterfaces.resize(0);
}


int	b3StateRigidCollidables::allocateCollidable()
{
	int curSize = m_collidablesCpu.size();
	if (curSize<m_limits.m_maxConvexShapes)
	{
		m_collidablesCpu.expand();
		return curSize;
	}
	else
	{
		b3Error("allocateCollidable out-of-range %d\n", m_limits.m_maxConvexShapes);
	}
	return B3_INVALID_COLLIDABLE_INDEX;

}

int b3StateRigidCollidables::registerFace(const b3Vector3& faceNormal, float faceConstant)
{
	int faceOffset = m_convexFaces.size();
	b3GpuFace& face = m_convexFaces.expand();
	face.m_plane = b3MakeVector3(faceNormal.x, faceNormal.y, faceNormal.z, faceConstant);
	return faceOffset;
}

int b3StateRigidCollidables::registerConvexHullShapeInternal(b3ConvexUtility* convexPtr, b3Collidable& col)
{

	m_convexData.resize(m_numCollidables + 1);
	m_convexPolyhedra.resize(m_numCollidables + 1);


	b3ConvexPolyhedronData& convex = m_convexPolyhedra.at(m_convexPolyhedra.size() - 1);
	convex.mC = convexPtr->mC;
	convex.mE = convexPtr->mE;
	convex.m_extents = convexPtr->m_extents;
	convex.m_localCenter = convexPtr->m_localCenter;
	convex.m_radius = convexPtr->m_radius;

	convex.m_numUniqueEdges = convexPtr->m_uniqueEdges.size();
	int edgeOffset = m_uniqueEdges.size();
	convex.m_uniqueEdgesOffset = edgeOffset;

	m_uniqueEdges.resize(edgeOffset + convex.m_numUniqueEdges);

	//convex data here
	for (int i = 0; i<convexPtr->m_uniqueEdges.size(); i++)
	{
		m_uniqueEdges[edgeOffset + i] = convexPtr->m_uniqueEdges[i];
	}

	int faceOffset = m_convexFaces.size();
	convex.m_faceOffset = faceOffset;
	convex.m_numFaces = convexPtr->m_faces.size();

	m_convexFaces.resize(faceOffset + convex.m_numFaces);


	for (int i = 0; i<convexPtr->m_faces.size(); i++)
	{
		m_convexFaces[convex.m_faceOffset + i].m_plane = b3MakeVector3(convexPtr->m_faces[i].m_plane[0],
			convexPtr->m_faces[i].m_plane[1],
			convexPtr->m_faces[i].m_plane[2],
			convexPtr->m_faces[i].m_plane[3]);


		int indexOffset = m_convexIndices.size();
		int numIndices = convexPtr->m_faces[i].m_indices.size();
		m_convexFaces[convex.m_faceOffset + i].m_numIndices = numIndices;
		m_convexFaces[convex.m_faceOffset + i].m_indexOffset = indexOffset;
		m_convexIndices.resize(indexOffset + numIndices);
		for (int p = 0; p<numIndices; p++)
		{
			m_convexIndices[indexOffset + p] = convexPtr->m_faces[i].m_indices[p];
		}
	}

	convex.m_numVertices = convexPtr->m_vertices.size();
	int vertexOffset = m_convexVertices.size();
	convex.m_vertexOffset = vertexOffset;

	m_convexVertices.resize(vertexOffset + convex.m_numVertices);
	for (int i = 0; i<convexPtr->m_vertices.size(); i++)
	{
		m_convexVertices[vertexOffset + i] = convexPtr->m_vertices[i];
	}

	m_convexData[m_numCollidables] = convexPtr;



	return m_numCollidables++;
}


int b3StateRigidCollidables::registerConcaveMeshShape(b3AlignedObjectArray<b3Vector3>* vertices, b3AlignedObjectArray<int>* indices, b3Collidable& col, const float* scaling1)
{
	b3Vector3 scaling = b3MakeVector3(scaling1[0], scaling1[1], scaling1[2]);

	m_convexData.resize(m_numCollidables + 1);
	m_convexPolyhedra.resize(m_numCollidables + 1);


	b3ConvexPolyhedronData& convex = m_convexPolyhedra.at(m_convexPolyhedra.size() - 1);
	convex.mC = b3MakeVector3(0, 0, 0);
	convex.mE = b3MakeVector3(0, 0, 0);
	convex.m_extents = b3MakeVector3(0, 0, 0);
	convex.m_localCenter = b3MakeVector3(0, 0, 0);
	convex.m_radius = 0.f;

	convex.m_numUniqueEdges = 0;
	int edgeOffset = m_uniqueEdges.size();
	convex.m_uniqueEdgesOffset = edgeOffset;

	int faceOffset = m_convexFaces.size();
	convex.m_faceOffset = faceOffset;

	convex.m_numFaces = indices->size() / 3;
	m_convexFaces.resize(faceOffset + convex.m_numFaces);
	m_convexIndices.reserve(convex.m_numFaces * 3);
	for (int i = 0; i<convex.m_numFaces; i++)
	{
		b3Vector3 vert0(vertices->at(indices->at(i * 3))*scaling);
		b3Vector3 vert1(vertices->at(indices->at(i * 3 + 1))*scaling);
		b3Vector3 vert2(vertices->at(indices->at(i * 3 + 2))*scaling);

		b3Vector3 normal = ((vert1 - vert0).cross(vert2 - vert0)).normalize();
		b3Scalar c = -(normal.dot(vert0));

		m_convexFaces[convex.m_faceOffset + i].m_plane = b3MakeVector4(normal.x, normal.y, normal.z, c);
		int indexOffset = m_convexIndices.size();
		int numIndices = 3;
		m_convexFaces[convex.m_faceOffset + i].m_numIndices = numIndices;
		m_convexFaces[convex.m_faceOffset + i].m_indexOffset = indexOffset;
		m_convexIndices.resize(indexOffset + numIndices);
		for (int p = 0; p<numIndices; p++)
		{
			int vi = indices->at(i * 3 + p);
			m_convexIndices[indexOffset + p] = vi;
		}
	}

	convex.m_numVertices = vertices->size();
	int vertexOffset = m_convexVertices.size();
	convex.m_vertexOffset = vertexOffset;
	m_convexVertices.resize(vertexOffset + convex.m_numVertices);
	for (int i = 0; i<vertices->size(); i++)
	{
		m_convexVertices[vertexOffset + i] = vertices->at(i)*scaling;
	}

	m_convexData[m_numCollidables] = 0;


	return m_numCollidables++;
}
