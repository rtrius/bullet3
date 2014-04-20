#include "GpuFractureScene.h"

#include "OpenGLWindow/GLInstancingRenderer.h"

#include "Bullet3Common/b3Random.h"
#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"
#include "Bullet3Geometry/b3AabbUtil.h"

#include "Bullet3OpenCL/RigidBody/b3GpuRigidBodyPipeline.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhase.h"
#include "Bullet3OpenCL/RigidBody/b3GpuNarrowPhaseInternalData.h"
#include "Bullet3OpenCL/Raycast/b3GpuRaycast.h"

#include "../GpuDemoInternalData.h"
#include "GpuRigidBodyDemoInternalData.h"
#include "b3RigidShapeStateUpdater.h"
#include "b3RigidBodyStateUpdater.h"

void GpuFractureScene::setupScene(const ConstructionInfo& ci)
{
	ConstructionInfo ci2 = ci;
	ci2.arraySizeX = 10;
	ci2.arraySizeY = 10;
	ci2.arraySizeZ = 10;
	
	m_primRenderer = ci.m_primRenderer;
	m_raycaster = new b3GpuRaycast(m_clData->m_clContext,m_clData->m_clDevice,m_clData->m_clQueue);

	GpuConvexScene::createStaticEnvironment(ci2);
	createDynamicsObjects(ci2);

	m_data->m_rigidBodyPipeline->writeAllInstancesToGpu();

	float camPos[4]={0,0,0,0};
	m_instancingRenderer->setCameraTargetPosition(camPos);
	m_instancingRenderer->setCameraDistance(150);
	m_instancingRenderer->setCameraYaw(30);
	m_instancingRenderer->setCameraPitch(225);
	
	m_instancingRenderer->updateCamera();
}
int GpuFractureScene::createDynamicsObjects(const ConstructionInfo& ci)
{
	return GpuBoxPlaneScene::createDynamicsObjects(ci);
}


void GpuFractureScene::destroyScene()
{
	delete m_raycaster;
	m_raycaster = 0;
}


void GpuFractureScene::renderScene()
{
	//	need to add compound/non-convex/sphere/plane support
	const bool RESET_RENDER_STATE_EVERY_FRAME = 1;
	if(RESET_RENDER_STATE_EVERY_FRAME)
	{
		B3_PROFILE("Reset rendersystem state");
		
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
			
			const b3GpuRigidBodyState* rigidState = m_data->m_np->getRigidBodyState();
			rigidState->m_usedRigidIndicesGPU->copyToHost(m_usedRigidIndicesCpu);
			
			const b3GpuNarrowPhaseInternalData* npInternalData = m_data->m_np->getInternalData();
			npInternalData->m_bodyBufferGPU->copyToHost(m_rigidBodiesCpu);
			
			npInternalData->m_collidablesGPU->copyToHost(m_collidableCpu);
			npInternalData->m_convexPolyhedraGPU->copyToHost(m_convexPolyhedraCpu);
			npInternalData->m_convexFacesGPU->copyToHost(m_convexFacesCpu);
			npInternalData->m_uniqueEdgesGPU->copyToHost(m_uniqueEdgesCpu);
			npInternalData->m_convexVerticesGPU->copyToHost(m_convexVerticesCpu);
			npInternalData->m_convexIndicesGPU->copyToHost(m_convexIndicesCpu);
		}
		
		//
		m_instancingRenderer->removeAllInstances();
		
		//Re-add convex shapes
		const int INVALID_RENDER_SHAPE = -1;
		
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
				const float COLOR[4] = { 0.0f, 0.7f, 1.0f, 1.0f };
				const float SCALING[3] = { 1.0f, 1.0f, 1.0f };
				
				int numRigidBodies = m_usedRigidIndicesCpu.size();
				
				int collidableIndex = i;
				for(int n = 0; n < numRigidBodies; ++n) 
				{
					const b3RigidBodyData& rigid = m_rigidBodiesCpu[ m_usedRigidIndicesCpu[n] ];
				
					//Render instances must be added in order, from lowest shape index to highest
					int renderShapeIndex = collidableToRenderShapeMapping[rigid.m_collidableIdx];
					if(renderShapeIndex != INVALID_RENDER_SHAPE && renderShapeIndex == collidableIndex)
						m_instancingRenderer->registerGraphicsInstance(renderShapeIndex, reinterpret_cast<const float*>(&rigid.m_pos), 
																	reinterpret_cast<const float*>(&rigid.m_quat), COLOR, SCALING);
				}
			}
		}
		
		m_instancingRenderer->writeTransforms();
	}
	
	GpuBoxPlaneScene::renderScene();
}

b3RigidShapeStateUpdater shapeUpdater;
b3RigidBodyStateUpdater rigidUpdater;
	
void GpuFractureScene::clientMoveAndDisplay()
{
	b3GpuNarrowPhaseInternalData* npInternalData = m_data->m_np->getInternalData();
	b3GpuRigidBodyState* rigidState = m_data->m_np->getRigidBodyState();
	b3GpuRigidShapeState* shapeState = m_data->m_np->getRigidShapeState();
	
	const bool TEST_RIGID_REMOVE = false;
	if(TEST_RIGID_REMOVE)
	{
		static b3AlignedObjectArray<b3RigidBodyData> m_rigidBodiesCpu;
		static b3AlignedObjectArray<int> m_usedRigidIndices;
		
		{
			B3_PROFILE("Transfer rigid data to CPU");
			npInternalData->m_bodyBufferGPU->copyToHost(m_rigidBodiesCpu);
			rigidState->m_usedRigidIndicesGPU->copyToHost(m_usedRigidIndices);
		}
		
		b3SapAabb removeAabb;
		removeAabb.m_minVec = b3MakeVector3(-10.f, 0.f, -10.f);
		removeAabb.m_maxVec = b3MakeVector3(10.f, 2.f, 10.f);
		
		for(int i = 0; i < m_usedRigidIndices.size(); ++i)
		{
			int activeIndex = m_usedRigidIndices[i];
		
			if( b3TestPointAgainstAabb2(removeAabb.m_minVec, removeAabb.m_maxVec, m_rigidBodiesCpu[activeIndex].m_pos) ) 
				rigidUpdater.markRigidBodyForRemoval(activeIndex);
		}
	}
	
	const bool TEST_RIGID_ADD = false;
	if(TEST_RIGID_ADD)
	{
		static int counter = 0;
		counter++;
		
		if(0 < counter && counter < 2000 && counter % 2)
		{
			const int COLLIDABLE_INDEX = 1;
			const b3Scalar MASS(1.0);
			rigidUpdater.addRigidBody(COLLIDABLE_INDEX, b3MakeVector3(5, 60, 5), b3Quaternion(0,0,0,1), MASS);
		}
	}
	
	const bool TEST_RIGID_ADD_AND_SHAPE_ADD = true;
	if(TEST_RIGID_ADD_AND_SHAPE_ADD)
	{
		static int counter = 0;
		counter++;
		
		if(0 < counter && counter < 2000 && counter % 2)
		{
			b3AlignedObjectArray<b3Vector3> vertices;
			vertices.push_back( b3MakeVector3(1, 1, 1) );
			vertices.push_back( b3MakeVector3(1, 1, -1) );
			vertices.push_back( b3MakeVector3(1, -1, 1) );
			vertices.push_back( b3MakeVector3(1, -1, -1) );
			
			vertices.push_back( b3MakeVector3(-1, 1, 1) );
			vertices.push_back( b3MakeVector3(-1, 1, -1) );
			vertices.push_back( b3MakeVector3(-1, -1, 1) );
			vertices.push_back( b3MakeVector3(-1, -1, -1) );
			
			for(int i = 0; i < vertices.size(); ++i) 
			{
				const b3Scalar MIN(-0.5);
				const b3Scalar MAX(0.5);
				vertices[i] += b3MakeVector3( b3RandRange(MIN, MAX), b3RandRange(MIN, MAX), b3RandRange(MIN, MAX) );
			}
			
			shapeUpdater.addShape(vertices);
			
			b3AlignedObjectArray<int> collidableIndices;
			shapeUpdater.applyUpdatesCpu(shapeState, npInternalData, &collidableIndices);
			
			const b3Scalar MASS(1.0);
			rigidUpdater.addRigidBody(collidableIndices[0], b3MakeVector3(5, 60, 5), b3Quaternion(0,0,0,1), MASS);
		}
	}
	
	rigidUpdater.applyUpdatesCpu(rigidState, npInternalData);
	
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
				rigidUpdater.markRigidBodyForRemoval(rigidBodyIndex);
			
				return true;
			}
		}
	}
	
	return false;
}
