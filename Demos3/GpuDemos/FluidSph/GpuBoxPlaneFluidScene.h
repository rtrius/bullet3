#ifndef GPU_BOX_PLANE_FLUID_SCENE_H
#define GPU_BOX_PLANE_FLUID_SCENE_H

#include "../rigidbody/GpuConvexScene.h"

#include "OpenGLWindow/GLInstancingRenderer.h"
#include "OpenGLWindow/ShapeData.h"
#include "Bullet3Common/b3Quaternion.h"
#include "../GpuDemoInternalData.h"
#include "../rigidbody/GpuRigidBodyDemoInternalData.h"

#include "ScreenSpaceFluidRendererGL.h"
#include "BulletFluids/Sph/b3FluidSph.h"
#include "BulletFluidsOpenCL/b3FluidSphSolverOpenCL.h"
#include "BulletFluidsOpenCL/b3FluidSphOpenCL.h"


#include "../rigidbody/ConcaveScene.h"
class GpuBoxPlaneFluidScene : public ConcaveScene //public GpuBoxPlaneScene
{
public:

	b3FluidSphSolverOpenCL* m_solver;
	b3FluidSphParametersGlobal m_globalParameters;
	b3FluidSph* m_sphFluid;

	ScreenSpaceFluidRendererGL* m_fluidRenderer;
	
	GpuBoxPlaneFluidScene()
	{
		m_sphFluid = new b3FluidSph(m_globalParameters, 131072);
		
		b3FluidSphParametersLocal FL = m_sphFluid->getLocalParameters();
		
		b3Scalar EXTENT(400.0);
		FL.m_aabbBoundaryMin = b3Vector3(-EXTENT, -EXTENT, -EXTENT);
		FL.m_aabbBoundaryMax = b3Vector3(EXTENT, EXTENT*b3Scalar(2.0), EXTENT);
		FL.m_enableAabbBoundary = 1;
		
		FL.m_particleMass = b3Scalar(0.01);
		//FL.m_boundaryErp = b3Scalar(0.9);
		FL.m_particleRadius = b3Scalar(4.0);
		
		m_sphFluid->setLocalParameters(FL);
		
		m_solver = 0;
		m_fluidRenderer = 0;
	}
	virtual ~GpuBoxPlaneFluidScene()
	{ 
		delete m_sphFluid; 
		
		if(m_solver) delete m_solver;
		if(m_fluidRenderer) delete m_fluidRenderer;
	}
	virtual const char* getName()
	{
		return "FluidTest";
	}

	static GpuDemo* MyCreateFunc()
	{
		GpuDemo* demo = new GpuBoxPlaneFluidScene;
		return demo;
	}
	
	virtual void setupScene(const ConstructionInfo& ci)
	{
		b3Assert(m_window);
		int width, height;
		m_window->getRenderingResolution(width, height);
		m_fluidRenderer = new ScreenSpaceFluidRendererGL(width, height);
	
		//GpuBoxPlaneScene::setupScene(ci);
		ConcaveScene::setupScene(ci);
		m_solver = new b3FluidSphSolverOpenCL(m_clData->m_clContext, m_clData->m_clDevice, m_clData->m_clQueue);
	}
	
	//virtual int	createDynamicsObjects(const ConstructionInfo& ci)
	virtual void	createDynamicObjects(const ConstructionInfo& ci)
	{
		//int numObjects = GpuBoxPlaneScene::createDynamicsObjects(ci);
		ConcaveScene::createDynamicObjects(ci);
		
		b3Vector3 OFFSET(125, 0, 0);
		b3Scalar EXTENT(45.0);
		b3Vector3 MIN(-EXTENT, b3Scalar(0.0), -EXTENT);
		b3Vector3 MAX(EXTENT, EXTENT, EXTENT);
		b3FluidEmitter::addVolume( m_sphFluid, MIN + OFFSET, MAX + OFFSET, b3Scalar(1.3) );
		
		//return numObjects;
	}
	
	virtual void renderScene()
	{
		//GpuBoxPlaneScene::renderScene();
		ConcaveScene::renderScene();
		/*
		float r = 0.5f;
		float g = 0.8f;
		float b = 1.0f;
		*/
		float r = 0.7f;
		float g = 0.7f;
		float b = 0.7f;
				
		//Beer's law constants
		//Controls the darkening of the fluid's color based on its thickness
		//For a constant k, (k > 1) == darkens faster; (k < 1) == darkens slower; (k == 0) == disable
		/*
		float absorptionR = 0.5;	
		float absorptionG = 0.5;
		float absorptionB = 0.5;
		*/
		float absorptionR = 0;	
		float absorptionG = 0;
		float absorptionB = 0;
		
		{
			B3_PROFILE("render fluid");
		
			const bool USE_MAPPED_BUFFER = false;
			if(USE_MAPPED_BUFFER)
			{
				b3Assert( sizeof(b3Vector3) == 16 );
				
				int numParticles = m_sphFluid->numParticles();
				int targetBufferSize = sizeof(b3Vector3) * numParticles;
				//int currentBufferSize = 0;
				//glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &currentBufferSize);
			
				//Resize VBO
				GLuint particlePositionVbo = m_fluidRenderer->getPositionVertexBuffer();
				
				glBindBuffer(GL_ARRAY_BUFFER, particlePositionVbo);
				//if(currentBufferSize < targetBufferSize) glBufferData(GL_ARRAY_BUFFER, targetBufferSize, 0, GL_DYNAMIC_DRAW);
				glBufferData(GL_ARRAY_BUFFER, targetBufferSize, 0, GL_DYNAMIC_DRAW);
				b3Vector3* glBuffer = static_cast<b3Vector3*>( glMapBufferRange(GL_ARRAY_BUFFER, 0, targetBufferSize, GL_MAP_WRITE_BIT) );
				
				//Copy particle positions from CL to GL buffer
				const void* sphOpenClObject = m_sphFluid->getFluidDataCL();
				const b3FluidSphOpenCL* sphDataCL = static_cast<const b3FluidSphOpenCL*>(sphOpenClObject);
				if(sphDataCL) sphDataCL->m_pos.copyToHostPointer( static_cast<b3Vector3*>(glBuffer), numParticles, 0 );
				else 
				{
					printf("sphDataCL: %d \n", sphDataCL);
					b3Assert(0);
				}
				
				glUnmapBuffer(GL_ARRAY_BUFFER);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
			}
			
			
			//m_fluidRenderer->render(m_sphFluid->getParticles().m_pos, 1.75f, r, g, b, absorptionR, absorptionG, absorptionB, !USE_MAPPED_BUFFER);
			m_fluidRenderer->render(m_sphFluid->getParticles().m_pos, 1.2f, r, g, b, absorptionR, absorptionG, absorptionB, !USE_MAPPED_BUFFER);
		}
	}
	
	virtual void clientMoveAndDisplay()
	{
		RigidBodyGpuData rbData;
		rbData.load(m_data->m_bp, m_data->m_np);
		
		m_solver->stepSimulation(m_globalParameters, &m_sphFluid, 1, rbData);
		
		static int counter = 0;
		if(++counter >= 100)
		{
			counter = 0;
			printf("m_sphFluid->numParticles(): %d \n", m_sphFluid->numParticles());
		}
		
		//GpuBoxPlaneScene::clientMoveAndDisplay();
		ConcaveScene::clientMoveAndDisplay();
	}
	
	virtual void resize(int width, int height)
	{
		if(m_fluidRenderer)
		{
			m_fluidRenderer->setWindowResolution(width, height);
			m_fluidRenderer->setRenderingResolution(width, height);
			//m_fluidRenderer->setRenderingResolution(width/2, height/2);
		}
	}

};

#endif