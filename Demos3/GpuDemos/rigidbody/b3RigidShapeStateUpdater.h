#ifndef B3_RIGID_SHAPE_STATE_UPDATER
#define B3_RIGID_SHAPE_STATE_UPDATER

#include "Bullet3Collision/NarrowPhaseCollision/b3ConvexUtility.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3ConvexPolyhedronData.h"

///Updates rigid body shapes; will cause crashes if shapes other than convex hulls are used
class b3RigidShapeStateUpdater
{
	b3AlignedObjectArray<b3ConvexUtility*> m_addedShapes;

	b3AlignedObjectArray<int> m_removedCollidableIndices;
	
	//Used only during applyUpdatesCpu()
	b3AlignedObjectArray<int> m_tempNewCollidableIndices;
	b3AlignedObjectArray<b3SapAabb> m_tempNewLocalAabbs;
	b3AlignedObjectArray<b3Collidable> m_tempNewCollidables;
	b3AlignedObjectArray<b3GpuFace> m_tempNewFaces;
	
	b3AlignedObjectArray<b3ConvexPolyhedronData> m_tempConvexPolyhedra;
	b3AlignedObjectArray<b3GpuFace> m_tempConvexFaces;
	b3AlignedObjectArray<b3Vector3> m_tempUniqueEdges;
	b3AlignedObjectArray<b3Vector3> m_tempConvexVertices;
	b3AlignedObjectArray<int> m_tempConvexIndices;
	
public:
	b3RigidShapeStateUpdater() {}
	virtual ~b3RigidShapeStateUpdater() {}
	
	///Returns a temporary index that is converted into a collidable index with applyUpdatesCpu() is called.
	int addShape(const b3AlignedObjectArray<b3Vector3>& convexHullVertices)
	{
		b3Assert( convexHullVertices.size() );
	
		b3ConvexUtility* utilPtr = new b3ConvexUtility();
		utilPtr->initializePolyhedralFeatures( &convexHullVertices[0], convexHullVertices.size() );

		int tempShapeIndex = m_addedShapes.size();
		m_addedShapes.push_back(utilPtr);	//delete in applyUpdatesCpu()
		
		return tempShapeIndex;
	}
	
	///@param collidableIndex Must be an active shape index; the same index cannot be used twice.
	void markShapeForRemoval(int collidableIndex) { m_removedCollidableIndices.push_back(collidableIndex); }
	
	///@param out_tempToCollidableIndexMap If pointer is nonzero, this array is loaded with the indices of created collidables.
	///addShape() returns an index to this array; for instance, if addShape() returns 5 then 
	///(*out_tempToCollidableIndexMap)[5] will contain the collidable index for that shape.
	void applyUpdatesCpu(b3GpuRigidShapeState* shapeState,
						b3GpuNarrowPhaseInternalData* narrowphaseInternalData,
						b3AlignedObjectArray<int>* out_tempToCollidableIndexMap = 0)
	{
		b3AlignedObjectArray<int>* availableIndicesCpu = shapeState->m_availableShapeIndicesCPU;
		b3AlignedObjectArray<int>* usedIndicesCpu = shapeState->m_usedShapeIndicesCPU;
		
		b3OpenCLArray<int>* availableIndicesGpu = shapeState->m_availableShapeIndicesGPU;
		b3OpenCLArray<int>* usedIndicesGpu = shapeState->m_usedShapeIndicesGPU;
		
		//Parallel Arrays; may have gaps/unused elements
		b3AlignedObjectArray<b3SapAabb>* localShapeAABBCpu = narrowphaseInternalData->m_localShapeAABBCPU;
		b3AlignedObjectArray<b3Collidable>& collidablesCpu = narrowphaseInternalData->m_collidablesCPU;
		
		//Contiguous Arrays
		b3AlignedObjectArray<b3ConvexPolyhedronData>& convexPolyhedraCpu = narrowphaseInternalData->m_convexPolyhedra;
		b3AlignedObjectArray<b3GpuFace>& convexFacesCpu = narrowphaseInternalData->m_convexFaces;
		b3AlignedObjectArray<b3Vector3>& uniqueEdgesCpu = narrowphaseInternalData->m_uniqueEdges;
		b3AlignedObjectArray<b3Vector3>& convexVerticesCpu = narrowphaseInternalData->m_convexVertices;
		b3AlignedObjectArray<int>& convexIndicesCpu = narrowphaseInternalData->m_convexIndices;
		
		//
		b3OpenCLArray<b3SapAabb>* localShapeAABBGpu = narrowphaseInternalData->m_localShapeAABBGPU;
		b3OpenCLArray<b3Collidable>* collidablesGpu = narrowphaseInternalData->m_collidablesGPU;
		
		b3OpenCLArray<b3ConvexPolyhedronData>* convexPolyhedraGpu = narrowphaseInternalData->m_convexPolyhedraGPU;
		b3OpenCLArray<b3GpuFace>* convexFacesGpu = narrowphaseInternalData->m_convexFacesGPU;
		b3OpenCLArray<b3Vector3>* uniqueEdgesGpu = narrowphaseInternalData->m_uniqueEdgesGPU;
		b3OpenCLArray<b3Vector3>* convexVerticesGpu = narrowphaseInternalData->m_convexVerticesGPU;
		b3OpenCLArray<int>* convexIndicesGpu = narrowphaseInternalData->m_convexIndicesGPU;
		
		//Not used by convex hull
		//b3OpenCLArray<b3GpuChildShape>* childShapesGpu = narrowphaseInternalData->m_gpuChildShapes;
		//b3OpenCLArray<b3BvhInfo>* bvhInfoGpu = narrowphaseInternalData->m_bvhInfoGPU;
		//b3OpenCLArray<b3QuantizedBvhNode>* treeNodesGpu = narrowphaseInternalData->m_treeNodesGPU;
		//b3OpenCLArray<b3BvhSubtreeInfo>* subTreesGpu = narrowphaseInternalData->m_subTreesGPU;
		
		//Not used
		//b3AlignedObjectArray<b3ConvexUtility*>* convexData = narrowphaseInternalData->m_convexData;
		
		
		//Assuming that the GPU data is more current than CPU
		{
			B3_PROFILE("Transfer shape data to CPU");
			
			availableIndicesGpu->copyToHost(*availableIndicesCpu);
			usedIndicesGpu->copyToHost(*usedIndicesCpu);
			
			localShapeAABBGpu->copyToHost(*localShapeAABBCpu);
			collidablesGpu->copyToHost(collidablesCpu);
			
			convexPolyhedraGpu->copyToHost(convexPolyhedraCpu);
			convexFacesGpu->copyToHost(convexFacesCpu);
			uniqueEdgesGpu->copyToHost(uniqueEdgesCpu);
			convexVerticesGpu->copyToHost(convexVerticesCpu);
			convexIndicesGpu->copyToHost(convexIndicesCpu);
		}
		
		//
		int numRemovedShapes = m_removedCollidableIndices.size();
		if(numRemovedShapes)
		{
			B3_PROFILE("Remove convex hull shapes");
			
			//Update available/used b3Collidable indices
			{
				for(int i = 0; i < numRemovedShapes; ++i)
				{
					int removedCollidableIndex = m_removedCollidableIndices[i];
					
					//Slow linear search results in worst case O(N^2)
					usedIndicesCpu->remove(removedCollidableIndex);
					availableIndicesCpu->push_back(removedCollidableIndex);
				}
				
				m_removedCollidableIndices.resize(0);
			}
			
			//Regenerate contiguous arrays with remaining collidables
			{
				m_tempConvexPolyhedra.resize(0);
				m_tempConvexFaces.resize(0);
				m_tempUniqueEdges.resize(0);
				m_tempConvexVertices.resize(0);
				m_tempConvexIndices.resize(0);
			 
				for(int i = 0; i < usedIndicesCpu->size(); ++i)
				{
					int collidableIndex = (*usedIndicesCpu)[i];
					b3Collidable& collidable = collidablesCpu[collidableIndex];
					b3Assert(collidable.m_shapeType == SHAPE_CONVEX_HULL);
				
					const b3ConvexPolyhedronData& polyhedra = convexPolyhedraCpu[collidable.m_shapeIndex];
					
					//Copy b3ConvexPolyhedronData and update offsets
					b3ConvexPolyhedronData updatedPolyhedra = polyhedra;
					updatedPolyhedra.m_faceOffset = m_tempConvexFaces.size();
					updatedPolyhedra.m_uniqueEdgesOffset = m_tempUniqueEdges.size();
					updatedPolyhedra.m_vertexOffset = m_tempConvexVertices.size();
					
					collidable.m_shapeIndex = m_tempConvexPolyhedra.size();
					m_tempConvexPolyhedra.push_back(updatedPolyhedra);
					
					//Copy faces, uniqueEdges, vertices
					{
						for(int n = 0; n < polyhedra.m_numFaces; ++n)
						{
							const b3GpuFace& face = convexFacesCpu[polyhedra.m_faceOffset + n];
							
							//
							b3GpuFace updatedFace = face;
							updatedFace.m_indexOffset = m_tempConvexIndices.size();
							m_tempConvexFaces.push_back(updatedFace);
							
							//Copy convex hull indices
							for(int j = 0; j < face.m_numIndices; ++j)
								m_tempConvexIndices.push_back( convexIndicesCpu[face.m_indexOffset + j] );
						}
						
						for(int n = 0; n < polyhedra.m_numUniqueEdges; ++n)
							m_tempUniqueEdges.push_back( uniqueEdgesCpu[polyhedra.m_uniqueEdgesOffset + n] );
						for(int n = 0; n < polyhedra.m_numVertices; ++n)
							m_tempConvexVertices.push_back( convexVerticesCpu[polyhedra.m_vertexOffset + n] );
					}
				}
				
				convexPolyhedraCpu = m_tempConvexPolyhedra;
				convexFacesCpu = m_tempConvexFaces;
				uniqueEdgesCpu = m_tempUniqueEdges;
				convexVerticesCpu = m_tempConvexVertices;
				convexIndicesCpu = m_tempConvexIndices;
			}
			
			narrowphaseInternalData->m_numAcceleratedShapes -= numRemovedShapes;
		}
		
		//
		int numAddedShapes = m_addedShapes.size();
		if(numAddedShapes)
		{
			B3_PROFILE("Create convex hull and add collidable");
			
			//Append new shape data to contiguous arrays
			{
				m_tempNewLocalAabbs.resize(numAddedShapes);
				m_tempNewCollidables.resize(numAddedShapes);
				
				for(int i = 0; i < numAddedShapes; ++i)
				{
					b3ConvexUtility* utilPtr = m_addedShapes[i];
					m_tempNewFaces.resize(0);
					
					//Update arrays pointed to by b3GpuFace
					{
						for(int j = 0; j < utilPtr->m_faces.size(); ++j)
						{
							const b3MyFace& face = utilPtr->m_faces[j];
						
							b3GpuFace gpuFace;
							gpuFace.m_plane = b3MakeVector3(face.m_plane[0], face.m_plane[1], face.m_plane[2], face.m_plane[3]);
							gpuFace.m_indexOffset = convexIndicesCpu.size();
							gpuFace.m_numIndices = face.m_indices.size();
							m_tempNewFaces.push_back(gpuFace);
						
							for(int k = 0; k < face.m_indices.size(); ++k) convexIndicesCpu.push_back(face.m_indices[k]);
						}
					}
					
					//Update arrays pointed to by b3ConvexPolyhedronData
					int faceOffset = convexFacesCpu.size();
					int uniqueEdgeOffset = uniqueEdgesCpu.size();
					int vertexOffset = convexVerticesCpu.size();
					{
						for(int n = 0; n < utilPtr->m_faces.size(); ++n) convexFacesCpu.push_back(m_tempNewFaces[n]);
						for(int n = 0; n < utilPtr->m_uniqueEdges.size(); ++n) uniqueEdgesCpu.push_back(utilPtr->m_uniqueEdges[n]);
						for(int n = 0; n < utilPtr->m_vertices.size(); ++n) convexVerticesCpu.push_back(utilPtr->m_vertices[n]);
					}
					
					//
					int newPolyhedronIndex = convexPolyhedraCpu.size();
					{
						b3ConvexPolyhedronData convex;
						convex.mC = utilPtr->mC;
						convex.mE = utilPtr->mE;
						convex.m_extents = utilPtr->m_extents;
						convex.m_localCenter = utilPtr->m_localCenter;
						convex.m_radius = utilPtr->m_radius;
						
						convex.m_numFaces = utilPtr->m_faces.size();
						convex.m_numUniqueEdges = utilPtr->m_uniqueEdges.size();
						convex.m_numVertices = utilPtr->m_vertices.size();
						convex.m_faceOffset = faceOffset;
						convex.m_uniqueEdgesOffset = uniqueEdgeOffset;
						convex.m_vertexOffset = vertexOffset;
						
						convexPolyhedraCpu.push_back(convex);
					}
					
					//
					{
						b3Collidable collidable;
						collidable.m_shapeType = SHAPE_CONVEX_HULL;
						collidable.m_shapeIndex = newPolyhedronIndex;
						m_tempNewCollidables[i] = collidable;
					}
					
					//
					{
						b3SapAabb localAabb;
						localAabb.m_minVec = b3MakeVector3(B3_LARGE_FLOAT, B3_LARGE_FLOAT, B3_LARGE_FLOAT);
						localAabb.m_maxVec = b3MakeVector3(-B3_LARGE_FLOAT, -B3_LARGE_FLOAT, -B3_LARGE_FLOAT);

						for(int n = 0; n < utilPtr->m_vertices.size(); n++)
						{
							localAabb.m_minVec.setMin(utilPtr->m_vertices[n]);
							localAabb.m_maxVec.setMax(utilPtr->m_vertices[n]);
						}
						localAabb.m_minIndices[3] = 0;
						localAabb.m_signedMaxIndices[3] = 0;
						m_tempNewLocalAabbs[i] = localAabb;
					}
				}
				
				for(int i = 0; i < m_addedShapes.size(); ++i) delete m_addedShapes[i];
				m_addedShapes.resize(0);
			}
		
			//Allocate indices for new collidables
			{
				m_tempNewCollidableIndices.resize(0);
				
				if( availableIndicesCpu->size() >= numAddedShapes )
				{
					for(int i = 0; i < numAddedShapes; ++i)
					{
						m_tempNewCollidableIndices.push_back( (*availableIndicesCpu)[ availableIndicesCpu->size() - 1 ] );
						availableIndicesCpu->pop_back();
					}
				}
				else
				{
					m_tempNewCollidableIndices = (*availableIndicesCpu);
					availableIndicesCpu->resize(0);
					
					int numNewIndices = numAddedShapes - availableIndicesCpu->size();
					
					int newArraySize = narrowphaseInternalData->m_numAcceleratedShapes + numNewIndices;
					if( newArraySize > collidablesCpu.size() )
					{
						localShapeAABBCpu->resize(newArraySize);
						collidablesCpu.resize(newArraySize);
					}
					
					int firstNewIndex = narrowphaseInternalData->m_numAcceleratedShapes;
					for(int i = 0; i < numNewIndices; ++i) m_tempNewCollidableIndices.push_back(firstNewIndex + i);
				}
			}
			
			if(out_tempToCollidableIndexMap) out_tempToCollidableIndexMap->resize(0);
			
			//Copy data to parallel arrays
			for(int i = 0; i < numAddedShapes; ++i)
			{
				int newCollidableIndex = m_tempNewCollidableIndices[i];
				(*localShapeAABBCpu)[newCollidableIndex] = m_tempNewLocalAabbs[i];
				collidablesCpu[newCollidableIndex] = m_tempNewCollidables[i];
				
				usedIndicesCpu->push_back(newCollidableIndex);
				if(out_tempToCollidableIndexMap) out_tempToCollidableIndexMap->push_back(newCollidableIndex);
			}
			
			narrowphaseInternalData->m_numAcceleratedShapes += numAddedShapes;
		}
		
		{
			B3_PROFILE("Transfer shape data to GPU");
			
			availableIndicesGpu->copyFromHost(*availableIndicesCpu);
			usedIndicesGpu->copyFromHost(*usedIndicesCpu);
			
			localShapeAABBGpu->copyFromHost(*localShapeAABBCpu);
			collidablesGpu->copyFromHost(collidablesCpu);
			
			convexPolyhedraGpu->copyFromHost(convexPolyhedraCpu);
			convexFacesGpu->copyFromHost(convexFacesCpu);
			uniqueEdgesGpu->copyFromHost(uniqueEdgesCpu);
			convexVerticesGpu->copyFromHost(convexVerticesCpu);
			convexIndicesGpu->copyFromHost(convexIndicesCpu);
		}
	}
};


#endif
