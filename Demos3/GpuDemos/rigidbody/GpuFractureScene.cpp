#include "GpuFractureScene.h"

#include "OpenGLWindow/GLInstancingRenderer.h"

#include "Bullet3Common/b3Random.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"
#include "Bullet3Geometry/b3AabbUtil.h"
#include "Bullet3Geometry/b3ConvexHullComputer.h"

#include "Bullet3OpenCL/RigidBody/b3GpuRigidBodyPipeline.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhase.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhaseInternalData.h"
#include "Bullet3OpenCL/Raycast/b3GpuRaycast.h"

#include "../GpuDemoInternalData.h"
#include "GpuRigidBodyDemoInternalData.h"
#include "b3RigidShapeStateUpdater.h"
#include "b3RigidBodyStateUpdater.h"

b3RigidShapeStateUpdater shapeUpdater;
b3RigidBodyStateUpdater rigidUpdater;

bool needsResetRenderState = true;		//Set true whenever a rigid body or shape is added/removed

void GpuFractureScene::setupScene(const ConstructionInfo& ci)
{
	needsResetRenderState = true;

	ConstructionInfo ci2 = ci;
	ci2.arraySizeX = 10;
	ci2.arraySizeY = 10;
	ci2.arraySizeZ = 10;
	
	m_primRenderer = ci.m_primRenderer;
	m_raycaster = new b3GpuRaycast(m_clData->m_clContext,m_clData->m_clDevice,m_clData->m_clQueue);

	GpuConvexScene::createStaticEnvironment(ci2);
	m_data->m_rigidBodyPipeline->writeAllInstancesToGpu();
	m_data->m_np->writeAllBodiesToGpu();
	
	createDynamicsObjects(ci2);


	float camPos[4]={0,0,0,0};
	m_instancingRenderer->setCameraTargetPosition(camPos);
	m_instancingRenderer->setCameraDistance(150);
	m_instancingRenderer->setCameraYaw(30);
	m_instancingRenderer->setCameraPitch(225);
	
	m_instancingRenderer->updateCamera();
}
int GpuFractureScene::createDynamicsObjects(const ConstructionInfo& ci)
{
	//Create a single dynamic rigid body box; each dynamic rigid body must have a unique shape for this demo
	{
		b3GpuNarrowPhaseInternalData* npInternalData = m_data->m_np->getInternalData();
		b3GpuRigidBodyState* rigidState = m_data->m_np->getRigidBodyState();
		b3GpuRigidShapeState* shapeState = m_data->m_np->getRigidShapeState();
	
		b3AlignedObjectArray<b3Vector3> vertices;
		vertices.push_back( b3MakeVector3(1, 1, 1) );
		vertices.push_back( b3MakeVector3(1, 1, -1) );
		vertices.push_back( b3MakeVector3(1, -1, 1) );
		vertices.push_back( b3MakeVector3(1, -1, -1) );
		
		vertices.push_back( b3MakeVector3(-1, 1, 1) );
		vertices.push_back( b3MakeVector3(-1, 1, -1) );
		vertices.push_back( b3MakeVector3(-1, -1, 1) );
		vertices.push_back( b3MakeVector3(-1, -1, -1) );
		
		const b3Scalar SCALE(50.0);
		for(int i = 0; i < vertices.size(); ++i) vertices[i] *= SCALE;
		
		shapeUpdater.addShape(vertices);
		
		b3AlignedObjectArray<int> collidableIndices;
		shapeUpdater.applyUpdatesCpu(shapeState, npInternalData, &collidableIndices);
	
		const b3Scalar MASS(10.0);
		rigidUpdater.addRigidBody(collidableIndices[0], b3MakeVector3(0, SCALE, 0), b3Quaternion(0,0,0,1), MASS);
		rigidUpdater.applyUpdatesCpu(rigidState, npInternalData);
		
		return 1;
	}
}


void GpuFractureScene::destroyScene()
{
	delete m_raycaster;
	m_raycaster = 0;
}

void GpuFractureScene::renderScene()
{
	const int INVALID_RENDER_SHAPE = -1;
	const int INVALID_RENDERER_INSTANCE = -1;

	const b3GpuRigidBodyState* rigidState = m_data->m_np->getRigidBodyState();
	const b3GpuNarrowPhaseInternalData* npInternalData = m_data->m_np->getInternalData();
	
	static b3AlignedObjectArray<int> m_rigidIndexToRenderInstanceIndexMap;

	//	compound/non-convex(trimesh)/sphere/plane support not implemented
	if(needsResetRenderState)
	{
		B3_PROFILE("Reset rendersystem state");
		
		needsResetRenderState = false;
		
		static b3AlignedObjectArray<int> m_usedRigidIndicesCpu;
		
		static b3AlignedObjectArray<b3RigidBodyData> m_rigidBodiesCpu;
		
		static b3AlignedObjectArray<b3Collidable> m_collidableCpu;
		static b3AlignedObjectArray<b3ConvexPolyhedronData> m_convexPolyhedraCpu;
		static b3AlignedObjectArray<b3GpuFace> m_convexFacesCpu;
		static b3AlignedObjectArray<b3Vector3> m_uniqueEdgesCpu;
		static b3AlignedObjectArray<b3Vector3> m_convexVerticesCpu;
		static b3AlignedObjectArray<int> m_convexIndicesCpu;
		
		{
			B3_PROFILE("Transfer rigid bodies/shapes to CPU");
			
			rigidState->m_usedRigidIndicesGPU->copyToHost(m_usedRigidIndicesCpu);
			
			npInternalData->m_bodyBufferGPU->copyToHost(m_rigidBodiesCpu);
			
			npInternalData->m_collidablesGPU->copyToHost(m_collidableCpu);
			npInternalData->m_convexPolyhedraGPU->copyToHost(m_convexPolyhedraCpu);
			npInternalData->m_convexFacesGPU->copyToHost(m_convexFacesCpu);
			npInternalData->m_uniqueEdgesGPU->copyToHost(m_uniqueEdgesCpu);
			npInternalData->m_convexVerticesGPU->copyToHost(m_convexVerticesCpu);
			npInternalData->m_convexIndicesGPU->copyToHost(m_convexIndicesCpu);
		}
		
		//
		int numRigidBodies = m_usedRigidIndicesCpu.size();
		{
			B3_PROFILE("remove all render instances");
			m_instancingRenderer->removeAllInstances();
			
			m_rigidIndexToRenderInstanceIndexMap.resize(numRigidBodies);
			for(int i = 0; i < numRigidBodies; ++i) m_rigidIndexToRenderInstanceIndexMap[i] = INVALID_RENDERER_INSTANCE;
		}
		
		//Re-add convex shapes
		int numCollidables = m_collidableCpu.size();
		static b3AlignedObjectArray<int> collidableToRenderShapeMapping;
		collidableToRenderShapeMapping.resize(numCollidables);
		for(int i = 0; i < numCollidables; ++i) collidableToRenderShapeMapping[i] = INVALID_RENDER_SHAPE;
		
		//Convert each convex shape into the renderer's format
		for(int i = 0; i < numCollidables; ++i)
		{
			static b3AlignedObjectArray<float> vertices;
			static b3AlignedObjectArray<int> triangleIndices;
			vertices.resize(0);
			triangleIndices.resize(0);
		
			const b3Collidable& collidable = m_collidableCpu[i];
			if(collidable.m_shapeType != SHAPE_CONVEX_HULL) continue;
			//else b3Assert(0);
		
			const b3ConvexPolyhedronData& convexHull = m_convexPolyhedraCpu[collidable.m_shapeIndex];
			int numVertices = convexHull.m_numVertices;
			
			b3Vector3 center = b3MakeVector3(0, 0, 0);
			{
				for(int n = 0; n < numVertices; ++n) center += m_convexVerticesCpu[ convexHull.m_vertexOffset + n ];
				center /= static_cast<float>(numVertices);
			}
			
			//Convert vertex format
			//For rendering, "vertices must be in the format x,y,z,w(unused), nx,ny,nz, u,v"
			for(int n = 0; n < numVertices; ++n)
			{
				const b3Vector3& vertex = m_convexVerticesCpu[ convexHull.m_vertexOffset + n ];
				vertices.push_back(vertex.x);
				vertices.push_back(vertex.y);
				vertices.push_back(vertex.z);
				vertices.push_back(0.0f);
				
				b3Vector3 normal = (vertex - center).normalized(); 
				vertices.push_back(normal.x);
				vertices.push_back(normal.y);
				vertices.push_back(normal.z);
				
				vertices.push_back(0.5f);
				vertices.push_back(0.5f);
			}
			
			//Convert each polygon in the convex hull into triangles
			for(int n = 0; n < convexHull.m_numFaces; ++n)
			{
				//Since it is a convex hull, the face is also a convex polygon
				const b3GpuFace& face = m_convexFacesCpu[ convexHull.m_faceOffset + n ];
			
				//Assuming that the indices specify vertices in clockwise or counter-clockwise order.
				//That is, that index[0] and index[1] form an edge, index[1] and index[2] forms the next edge, and so on.
				//In such a case, the polygon can be converted into triangles using
				//indices (0,1,2), (0,2,3), (0,3,4), and so on.
				for(int j = 2; j < face.m_numIndices; ++j)
				{
					triangleIndices.push_back( m_convexIndicesCpu[face.m_indexOffset + 0] );
					triangleIndices.push_back( m_convexIndicesCpu[face.m_indexOffset + j - 1] );
					triangleIndices.push_back( m_convexIndicesCpu[face.m_indexOffset + j] );
				}
			}
			
			int renderShapeIndex = m_instancingRenderer->registerShape( reinterpret_cast<const float*>(&vertices[0]), numVertices, 
																		&triangleIndices[0], triangleIndices.size() );
			collidableToRenderShapeMapping[i] = renderShapeIndex;
			
			
			//Re-add rigid body instances
			//Current renderer implementation requires that 
			//all instances of a shape are added right after registerShape()
			{
				B3_PROFILE("Create rigid body instances");
			
				b3Vector4 COLORS[4] =
				{
					b3MakeVector4(0,1,1,1),
					b3MakeVector4(1,0,0,1),
					b3MakeVector4(0,1,0,1),
					b3MakeVector4(1,1,0,1),
				};
				
				const float SCALING[3] = { 1.0f, 1.0f, 1.0f };
				
				int collidableIndex = i;
				for(int n = 0; n < numRigidBodies; ++n) 
				{
					const b3RigidBodyData& rigid = m_rigidBodiesCpu[ m_usedRigidIndicesCpu[n] ];
				
					//Render instances must be added in order, from lowest shape index to highest
					int renderShapeIndex = collidableToRenderShapeMapping[rigid.m_collidableIdx];
					if(renderShapeIndex != INVALID_RENDER_SHAPE && renderShapeIndex == collidableIndex)
					{
						int renderInstanceIndex = m_instancingRenderer->registerGraphicsInstance(renderShapeIndex, reinterpret_cast<const float*>(&rigid.m_pos), 
																	reinterpret_cast<const float*>(&rigid.m_quat), COLORS[n % 4], SCALING);
						m_rigidIndexToRenderInstanceIndexMap[ m_usedRigidIndicesCpu[n] ] = renderInstanceIndex;
					}
				}
			}
		}
		
	}
	else
	{
		static b3AlignedObjectArray<b3RigidBodyData> rigidBodies;
		npInternalData->m_bodyBufferGPU->copyToHost(rigidBodies);
		
		int numRigidBodies = rigidState->m_usedRigidIndicesCPU->size();
		for(int i = 0; i < numRigidBodies; ++i)
		{
			int renderInstanceIndex = m_rigidIndexToRenderInstanceIndexMap[i];
			if(renderInstanceIndex != INVALID_RENDERER_INSTANCE)
			{
				m_instancingRenderer->writeSingleInstanceTransformToCPU(&rigidBodies[i].m_pos.x, &rigidBodies[i].m_quat.x, renderInstanceIndex);
			}
		}
	}
	
	{
		B3_PROFILE("write transforms renderer");
		m_instancingRenderer->writeTransforms();
	}
	
	GpuBoxPlaneScene::renderScene();
}


#include <set>

//Voronoi fracture and shatter code copyright (c) 2011 Alain Ducharme
//  - Note that demo's visual cracks between voronoi shards are NOT present in the internally generated voronoi mesh!
struct VoronoiFracture
{
	b3AlignedObjectArray<int> m_convexHullVerticesOffset;
	b3AlignedObjectArray<int> m_convexHullVerticesCount;
	b3AlignedObjectArray<b3Vector3> m_convexHullVertices;
	
	b3AlignedObjectArray<b3Vector3> m_rigidPosition;
	b3AlignedObjectArray<b3Scalar> m_rigidMass;

	static void getVerticesInsidePlanes(const b3AlignedObjectArray<b3Vector3>& planes,
										b3AlignedObjectArray<b3Vector3>& verticesOut,
										std::set<int>& planeIndicesOut)
	{
		const b3Scalar TOLERANCE(0.0001);	//Distance at which a point and halfspace(plane) are considered to be not intersecting
	
		// Based on b3GeometryUtil.cpp (Gino van den Bergen / Erwin Coumans)
		verticesOut.resize(0);
		planeIndicesOut.clear();
		const int numPlanes = planes.size();
		int i, j, k, l;
		for (i=0;i<numPlanes;i++)
		{
			const b3Vector3& N1 = planes[i];
			for (j=i+1;j<numPlanes;j++)
			{
				const b3Vector3& N2 = planes[j];
				b3Vector3 n1n2 = N1.cross(N2);
				if (n1n2.length2() > b3Scalar(0.0001))
				{
					for (k=j+1;k<numPlanes;k++)
					{
						const b3Vector3& N3 = planes[k];
						b3Vector3 n2n3 = N2.cross(N3);
						b3Vector3 n3n1 = N3.cross(N1);
						if ((n2n3.length2() > b3Scalar(0.0001)) && (n3n1.length2() > b3Scalar(0.0001) ))
						{
							b3Scalar quotient = (N1.dot(n2n3));
							if (b3Fabs(quotient) > b3Scalar(0.0001))
							{
								b3Vector3 potentialVertex = (n2n3 * N1[3] + n3n1 * N2[3] + n1n2 * N3[3]) * (b3Scalar(-1.) / quotient);
								for (l=0; l<numPlanes; l++)
								{
									const b3Vector3& NP = planes[l];
									if (b3Scalar(NP.dot(potentialVertex))+b3Scalar(NP[3]) > TOLERANCE)
										break;
								}
								if (l == numPlanes)
								{
									// vertex (three plane intersection) inside all planes
									verticesOut.push_back(potentialVertex);
									planeIndicesOut.insert(i);
									planeIndicesOut.insert(j);
									planeIndicesOut.insert(k);
								}
							}
						}
					}
				}
			}
		}
	}
	
	struct DistanceSortPredicate
	{
		b3Vector3 m_point;
	
		DistanceSortPredicate(const b3Vector3& point) : m_point(point) {}
	
		bool operator()(const b3Vector3& p1, const b3Vector3& p2) const
		{
			float v1 = (p1-m_point).length2();
			float v2 = (p2-m_point).length2();
			bool result0 = v1 < v2;
			//bool result1 = ((b3Scalar)(p1-m_point).length2()) < ((b3Scalar)(p2-m_point).length2());
			//apparently result0 is not always result1, because extended precision used in registered is different from precision when values are stored in memory
			return result0;
		}
	};
	
	// points define voronoi cells in world space (avoid duplicates)
	// verts = source (convex hull) mesh vertices in local space
	// bbq & bbt = source (convex hull) mesh quaternion rotation and translation
	// matDensity = Material density for voronoi shard mass calculation
	void voronoiConvexHullShatter(const b3AlignedObjectArray<b3Vector3>& points,
									const b3AlignedObjectArray<b3Vector3>& verts,
									const b3Quaternion& bbq, const b3Vector3& bbt, b3Scalar matDensity)
	{
		m_convexHullVerticesOffset.resize(0);
		m_convexHullVerticesCount.resize(0);
		m_convexHullVertices.resize(0);
		m_rigidPosition.resize(0);
		m_rigidMass.resize(0);
	
		b3ConvexHullComputer convexHC;
		b3AlignedObjectArray<b3Vector3> vertices, chverts;
		b3AlignedObjectArray<b3Vector3> sortedVoronoiPoints;
		b3AlignedObjectArray<b3Vector3> planes, convexPlanes;
		std::set<int> planeIndices;
		
		sortedVoronoiPoints.copyFromArray(points);
		
		// Convert verts to world space and get convexPlanes
		int numverts = verts.size();
		chverts.resize(verts.size());
		for (int i=0; i < numverts; i++) chverts[i] = b3QuatRotate(bbq, verts[i]) + bbt;
		
		//b3GeometryUtil::getPlaneEquationsFromVertices(chverts, convexPlanes);
		// Using convexHullComputer faster than getPlaneEquationsFromVertices for large meshes...
		convexHC.compute(&chverts[0].getX(), sizeof(b3Vector3), numverts, 0.0, 0.0);
		int numFaces = convexHC.faces.size();
		
		for (int i=0; i < numFaces; i++)
		{
			const b3ConvexHullComputer::Edge* edge = &convexHC.edges[convexHC.faces[i]];
			int v0 = edge->getSourceVertex();
			int v1 = edge->getTargetVertex();
			edge = edge->getNextEdgeOfFace();
			int v2 = edge->getTargetVertex();
			b3Vector3 plane = (convexHC.vertices[v1]-convexHC.vertices[v0]).cross(convexHC.vertices[v2]-convexHC.vertices[v0]).normalize();
			plane[3] = -plane.dot(convexHC.vertices[v0]);
			convexPlanes.push_back(plane);
		}
		

		int numpoints = points.size();
		
		//	construct voronoi shards
		int cellnum = 0;
		for (int i=0; i < numpoints ;i++)
		{
			b3Vector3 curVoronoiPoint = points[i];
			planes.copyFromArray(convexPlanes);
			
			int numconvexPlanes = convexPlanes.size();
			for (int j=0; j < numconvexPlanes ;j++) planes[j][3] += planes[j].dot(curVoronoiPoint);	//	point-plane distance
			
			b3Scalar maxDistance = B3_INFINITY;
			
			//Sort points by distance to curVoronoiPoint(closest to farthest) 
			sortedVoronoiPoints.heapSort( DistanceSortPredicate(curVoronoiPoint) );
			
			//	generate a set of vertices for a voronoi region (alternatively, get the set of planes that compose a voronoi region)
			for (int j=1; j < numpoints; j++)
			{
				//	the algorithm can be thought of as clipping away volume(removing and adding planes and vertices)
				//	from the original convex hull, each voronoi point generates a plane
				b3Vector3 normal = sortedVoronoiPoints[j] - curVoronoiPoint;
				b3Scalar nlength = normal.length();
				
				//	since sortedVoronoiPoints is sorted by distance, first element out of range
				//	means that all following are also out of range
				if (nlength > maxDistance) break;
				
				//	create plane(pointing outwards i.e. negative == penetrating), and calculate distance from curVoronoiPoint to point in plane
				b3Vector3 plane = normal.normalized();
				plane[3] = -nlength / b3Scalar(2.);
				planes.push_back(plane);
				
				//	(clip away volume not enclosed by all planes)
				//	clear vertices/planeIndices
				//	get all vertices resulting from all plane triplet intersections
				//	of those vertices, keep only those enclosed by all planes
				//	load planeIndices with indices of planes that contribute to the remaining vertices
				getVerticesInsidePlanes(planes, vertices, planeIndices);
				if (vertices.size() == 0) break;
				
				//	load planes with planes marked in planeIndices ( planeIndices.size() <= planes.size() always )
				//	(remove planes that do not contribute to a vertex)
				int numplaneIndices = planeIndices.size();
				if (numplaneIndices != planes.size())
				{
					std::set<int>::iterator planeIndicesIter = planeIndices.begin();
					for (int k=0; k < numplaneIndices; k++)
					{
						if (k != *planeIndicesIter)
							planes[k] = planes[*planeIndicesIter];
						planeIndicesIter++;
					}
					planes.resize(numplaneIndices);
				}
				
				//	load maxdistance with largest distance vertex x2
				//	(x2 distance since we do not test against the plane but the voronoi point)
				maxDistance = vertices[0].length();
				for (int k=1; k < vertices.size(); k++)
				{
					b3Scalar distance = vertices[k].length();
					if (maxDistance < distance) maxDistance = distance;
				}
				maxDistance *= b3Scalar(2.);
			}
			
			if (vertices.size() == 0) continue;
			
			// Clean-up voronoi convex shard vertices and generate edges & faces
			convexHC.compute(&vertices[0].getX(), sizeof(b3Vector3), vertices.size(),0.0,0.0);

			// At this point we have a complete 3D voronoi shard mesh contained in convexHC

			// Calculate volume and center of mass (Stan Melax volume integration)
			numFaces = convexHC.faces.size();
			b3Scalar volume = b3Scalar(0.);
			b3Vector3 com = b3MakeVector3(0, 0, 0);
			for (int j=0; j < numFaces; j++)
			{
				const b3ConvexHullComputer::Edge* edge = &convexHC.edges[convexHC.faces[j]];
				int v0 = edge->getSourceVertex();
				int v1 = edge->getTargetVertex();
				edge = edge->getNextEdgeOfFace();
				int v2 = edge->getTargetVertex();
				while (v2 != v0)
				{
					// Counter-clockwise triangulated voronoi shard mesh faces (v0-v1-v2) and edges here...
					b3Scalar vol = convexHC.vertices[v0].triple(convexHC.vertices[v1], convexHC.vertices[v2]);
					volume += vol;
					com += vol * (convexHC.vertices[v0] + convexHC.vertices[v1] + convexHC.vertices[v2]);
					edge = edge->getNextEdgeOfFace();
					v1 = v2;
					v2 = edge->getTargetVertex();
				}
			}
			com /= volume * b3Scalar(4.);
			volume /= b3Scalar(6.);

			// Shift all vertices relative to center of mass
			int numVerts = convexHC.vertices.size();
			for (int j=0; j < numVerts; j++) convexHC.vertices[j] -= com;
			
			// Create rigid body shards
			b3Vector3 shardPosition = curVoronoiPoint + com; // Shard's adjusted location
			b3Scalar shardMass(volume * matDensity);
			
			m_convexHullVerticesOffset.push_back( m_convexHullVertices.size() );
			m_convexHullVerticesCount.push_back(numVerts);
			for (int j = 0; j < numVerts; ++j) m_convexHullVertices.push_back(convexHC.vertices[j]);
			m_rigidPosition.push_back(shardPosition);
			m_rigidMass.push_back(shardMass);	
			
			cellnum ++;
		}
		printf("Generated %d voronoi b3RigidBody shards\n", cellnum);
	}
};
	
void GpuFractureScene::clientMoveAndDisplay()
{
	{
		B3_PROFILE("stepSimulation");
		m_data->m_rigidBodyPipeline->stepSimulation(1./60.f);
	}
}


bool GpuFractureScene::mouseMoveCallback(float x,float y)
{
	return false;
}

bool GpuFractureScene::mouseButtonCallback(int button, int state, float x, float y)
{
	//Button 0 == Left Mouse, State 1 == Pressed
	if ( button == 0 && state == 1 && (m_data->m_altPressed == 0 && m_data->m_controlPressed == 0) )
	{
		b3AlignedObjectArray<b3RayInfo> rays;
		b3AlignedObjectArray<b3RayHit> hitResults;
		
		b3Vector3 camPos;
		m_instancingRenderer->getCameraPosition(camPos);

		b3RayInfo ray;
		ray.m_from = camPos;
		ray.m_to = getRayTo(x,y);
		rays.push_back(ray);
		
		b3RayHit hit;
		hit.m_hitFraction = 1.f;
		hitResults.push_back(hit);
		
		m_data->m_rigidBodyPipeline->castRays(rays, hitResults);
		if ( hitResults[0].m_hitFraction < b3Scalar(1.0) )
		{
			int rigidBodyIndex = hitResults[0].m_hitBody;
			b3RigidBodyData rigidBody = m_data->m_np->getInternalData()->m_bodyBufferGPU->at(rigidBodyIndex);
			
			if ( rigidBody.m_invMass != b3Scalar(0.0) )
			{
				b3GpuNarrowPhaseInternalData* npInternalData = m_data->m_np->getInternalData();
				b3GpuRigidBodyState* rigidState = m_data->m_np->getRigidBodyState();
				b3GpuRigidShapeState* shapeState = m_data->m_np->getRigidShapeState();
					
				rigidUpdater.markRigidBodyForRemoval(rigidBodyIndex);
				
				//Create multiple new rigids using voronoi fracture
				{
					static b3AlignedObjectArray<b3Vector3> convexHullVertices;
					static b3AlignedObjectArray<b3Vector3> voronoiPoints;
					static b3AlignedObjectArray<b3Vector3> tempVertices;
					convexHullVertices.resize(0);
					voronoiPoints.resize(0);
					tempVertices.resize(0);
					
					b3Collidable collidable = npInternalData->m_collidablesGPU->at(rigidBody.m_collidableIdx);
					b3Assert(collidable.m_shapeType == SHAPE_CONVEX_HULL);
					
					b3ConvexPolyhedronData convexData = npInternalData->m_convexPolyhedraGPU->at(collidable.m_shapeIndex);
					
					//Assuming that vertex data on CPU and GPU is the same(syncronized)
					for(int i = 0; i < convexData.m_numVertices; ++i) 
						convexHullVertices.push_back( npInternalData->m_convexVertices[convexData.m_vertexOffset + i] );
						
					shapeUpdater.markShapeForRemoval(rigidBody.m_collidableIdx);
					
					//Generate voronoi points in world space
					const int NUM_VORONOI_POINTS = 32;
					for(int i = 0; i < NUM_VORONOI_POINTS; ++i) 
					{
						//Assuming that the local AABB on CPU and GPU is the same(syncronized)
						const b3SapAabb& localAabb = npInternalData->m_localShapeAABBCPU->at(rigidBody.m_collidableIdx);
						const b3Vector3& min = localAabb.m_minVec;
						const b3Vector3& max = localAabb.m_maxVec;
						
						b3Vector3 randomInRigidSpace = b3MakeVector3( b3RandRange(min.x, max.x), b3RandRange(min.y, max.y), b3RandRange(min.z, max.z) );
						b3Vector3 randomInWorldSpace = b3QuatRotate(rigidBody.m_quat, randomInRigidSpace) + rigidBody.m_pos;
						
						voronoiPoints.push_back(randomInWorldSpace);
					}
					
					//Compute voronoi fracture; use voronoi points to generate convex hulls
					static VoronoiFracture fracturer;
					
					const b3Scalar DENSITY(1.0);
					fracturer.voronoiConvexHullShatter(voronoiPoints, convexHullVertices, rigidBody.m_quat, rigidBody.m_pos, DENSITY);
					
					//Create shapes
					int numCreatedShapes = fracturer.m_convexHullVerticesCount.size();
					for(int i = 0; i < numCreatedShapes; ++i)
					{
						tempVertices.resize(0);
					
						int offset = fracturer.m_convexHullVerticesOffset[i];
						int numVertices = fracturer.m_convexHullVerticesCount[i];
						for(int n = 0; n < numVertices; ++n) tempVertices.push_back( fracturer.m_convexHullVertices[offset + n] );
						shapeUpdater.addShape(tempVertices);
					}
					
					//Apply shape updates
					b3AlignedObjectArray<int> collidableIndices;
					shapeUpdater.applyUpdatesCpu(shapeState, npInternalData, &collidableIndices);
					
					//Create rigid bodies
					int numCreatedRigids = fracturer.m_rigidPosition.size();
					for(int i = 0; i < numCreatedRigids; ++i)
					{
						b3Vector3 relativePosition = fracturer.m_rigidPosition[i] - rigidBody.m_pos;
						b3Vector3 linearVelocity = rigidBody.m_linVel + rigidBody.m_angVel.cross(relativePosition);
						
						b3Scalar mass = 1.0; //fracturer.m_rigidMass[i];
						rigidUpdater.addRigidBody(collidableIndices[i], fracturer.m_rigidPosition[i],
													b3Quaternion(0,0,0,1), mass, linearVelocity, rigidBody.m_angVel);
					}
				}
			
				rigidUpdater.applyUpdatesCpu(rigidState, npInternalData);
				printf("m_numAcceleratedRigidBodies: %d\n", npInternalData->m_numAcceleratedRigidBodies);
				
				needsResetRenderState = true;
				return true;
			}
		}
	}
	
	return false;
}
