#ifndef GPU_BOX_PLANE_FLUID_SCENE_H
#define GPU_BOX_PLANE_FLUID_SCENE_H

#include "../rigidbody/GpuConvexScene.h"

#include "OpenGLWindow/GLInstancingRenderer.h"
#include "OpenGLWindow/ShapeData.h"
#include "Bullet3Common/b3Quaternion.h"
#include "../GpuDemoInternalData.h"

#include "ScreenSpaceFluidRendererGL.h"
#include "BulletFluids/Sph/b3FluidSph.h"
#include "BulletFluidsOpenCL/b3FluidSphSolverOpenCL.h"
class GpuBoxPlaneFluidScene : public GpuBoxPlaneScene
{
public:

	b3FluidSphSolverOpenCL* m_solver;
	b3FluidSphParametersGlobal m_globalParameters;
	b3FluidSph *m_sphFluid;

	GpuBoxPlaneFluidScene()
	{
		m_sphFluid = new b3FluidSph(m_globalParameters, 65536);
		
		b3FluidSphParametersLocal FL = m_sphFluid->getLocalParameters();
		
		b3Scalar EXTENT(35.0);
		FL.m_aabbBoundaryMin = b3Vector3(-EXTENT, 0, -EXTENT);
		FL.m_aabbBoundaryMax = b3Vector3(EXTENT, EXTENT*b3Scalar(2.0), EXTENT);
		FL.m_enableAabbBoundary = 1;
		
		m_sphFluid->setLocalParameters(FL);
		
		m_solver = 0;
	}
	virtual ~GpuBoxPlaneFluidScene()
	{ 
		delete m_sphFluid; 
		if(m_solver) delete m_solver;
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
		GpuBoxPlaneScene::setupScene(ci);
		m_solver = new b3FluidSphSolverOpenCL(m_clData->m_clContext, m_clData->m_clQueue, m_clData->m_clDevice);
	}
	
	virtual int	createDynamicsObjects(const ConstructionInfo& ci)
	{
		int numObjects = GpuBoxPlaneScene::createDynamicsObjects(ci);
		
		b3Scalar EXTENT(30.0);
		b3Vector3 MIN(-EXTENT, EXTENT, -EXTENT);
		b3Vector3 MAX(EXTENT, EXTENT*b3Scalar(2.0), EXTENT);
		b3FluidEmitter::addVolume( m_sphFluid, MIN, MAX, b3Scalar(1.0) );
		
		return numObjects;
	}
	
	virtual void renderScene()
	{
		GpuBoxPlaneScene::renderScene();
		
		static ScreenSpaceFluidRendererGL* fluidRenderer = 0;
		if(!fluidRenderer) fluidRenderer = new ScreenSpaceFluidRendererGL(1024, 768);
		//fluidRenderer->setRenderingResolution(512, 384);
		
		float r = 0.5f;
		float g = 0.8f;
		float b = 1.0f;
				
		//Beer's law constants
		//Controls the darkening of the fluid's color based on its thickness
		//For a constant k, (k > 1) == darkens faster; (k < 1) == darkens slower; (k == 0) == disable
		float absorptionR = 2.5;	
		float absorptionG = 1.0;
		float absorptionB = 0.5;
		
		fluidRenderer->render(m_sphFluid->getParticles().m_pos, 1.2f, r, g, b, absorptionR, absorptionG, absorptionB);
	}
	
	virtual void clientMoveAndDisplay()
	{
		m_solver->updateGridAndCalculateSphForces(m_globalParameters, &m_sphFluid, 1);
		
		static int counter = 0;
		if(++counter >= 100)
		{
			counter = 0;
			printf("m_sphFluid->numParticles(): %d \n", m_sphFluid->numParticles());
		}
		//printf("m_sphFluid->getParticles().m_vel[0].y: %f \n", m_sphFluid->getParticles().m_vel[0].y);
		
		GpuBoxPlaneScene::clientMoveAndDisplay();
	}

};

#endif