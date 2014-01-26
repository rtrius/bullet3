#ifndef GPU_BOX_PLANE_FLUID_SCENE_H
#define GPU_BOX_PLANE_FLUID_SCENE_H

#include "../rigidbody/GpuConvexScene.h"

#include "OpenGLWindow/GLInstancingRenderer.h"
#include "OpenGLWindow/ShapeData.h"
#include "Bullet3Common/b3Quaternion.h"
#include "../GpuDemoInternalData.h"
#include "../rigidbody/GpuRigidBodyDemoInternalData.h"

#include "ScreenSpaceFluidRendererGL.h"
#include "Bullet3Fluids/Sph/b3FluidSph.h"
#include "Bullet3FluidsOpenCL/b3FluidSphSolverOpenCL.h"
#include "Bullet3FluidsOpenCL/b3FluidSphOpenCL.h"


#include "../rigidbody/ConcaveScene.h"
#include "../rigidbody/GpuCompoundScene.h"

//#define BASE_DEMO_CLASS ConcaveScene
#define BASE_DEMO_CLASS GpuBoxPlaneScene
//#define BASE_DEMO_CLASS GpuCompoundPlaneScene
const bool CONCAVE_SCENE = false;

void b3Matrix4x4Mul16(const float aIn[16], const float bIn[16], float result[16])
{
	for (int j=0;j<4;j++)
		for (int i=0;i<4;i++)
			result[j*4+i] = aIn[0*4+i] * bIn[j*4+0] + aIn[1*4+i] * bIn[j*4+1] + aIn[2*4+i] * bIn[j*4+2] + aIn[3*4+i] * bIn[j*4+3];
}

class GpuBoxPlaneFluidScene : public BASE_DEMO_CLASS
{
public:
	static const int MAX_FLUID_PARTICLES = 131072;
	static const int INITIAL_NUM_FLUID_PARTICLES = 131072;

	b3FluidSphSolverOpenCL* m_solver;

	b3FluidSph* m_sphFluid;
	
	GpuBoxPlaneFluidScene()
	{
		m_sphFluid = new b3FluidSph(0, MAX_FLUID_PARTICLES);
		
		b3FluidSphParameters FP = m_sphFluid->getParameters();
		
		b3Scalar EXTENT(100.0);
		if(CONCAVE_SCENE) EXTENT = b3Scalar(400.0);
		
		FP.m_aabbBoundaryMin = b3MakeVector3(-EXTENT, -EXTENT, -EXTENT);
		FP.m_aabbBoundaryMax = b3MakeVector3(EXTENT, EXTENT*b3Scalar(2.0), EXTENT);
		FP.m_enableAabbBoundary = 1;
		
		FP.m_particleMass = b3Scalar(0.005);
		
		if(CONCAVE_SCENE) 
		{
			FP.m_boundaryErp = b3Scalar(0.05);
			FP.m_particleRadius = b3Scalar(4.0);
		}
		
		m_sphFluid->setParameters(FP);
		
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
		BASE_DEMO_CLASS::setupScene(ci);
		
		m_solver = new b3FluidSphSolverOpenCL(m_clData->m_clContext, m_clData->m_clDevice, m_clData->m_clQueue);
		
		{
			b3Vector3 OFFSET = b3MakeVector3(0, 0, 0);
			b3Scalar EXTENT(90.0);
			
			if(CONCAVE_SCENE)
			{
				OFFSET = b3MakeVector3(125, 0, 0);
				EXTENT = b3Scalar(45.0);
			}
			
			b3Vector3 MIN = b3MakeVector3(-EXTENT, b3Scalar(0.0), -EXTENT);
			b3Vector3 MAX = b3MakeVector3(EXTENT, EXTENT, EXTENT);
			b3FluidEmitter::addVolume( m_sphFluid, MIN + OFFSET, MAX + OFFSET, b3Scalar(1.3), INITIAL_NUM_FLUID_PARTICLES );
		}
	}
	
	
	virtual void renderScene()
	{
		BASE_DEMO_CLASS::renderScene();
		
		B3_PROFILE("render fluid");
		
		const float SPHERE_SIZE(1.25);
		const float COLOR[4] = {0.5f, 0.8f, 1.0f, 1.0f};
		
		const bool USE_BULLET3_RENDERER = 1;
		if(USE_BULLET3_RENDERER)
		{
			int numParticles = m_sphFluid->numParticles();
			if(numParticles)
			{
				const float* positions = reinterpret_cast<const float*>(&m_sphFluid->getParticles().m_position[0]);
				m_instancingRenderer->drawSpheres(positions, COLOR, numParticles, sizeof(float)*4, SPHERE_SIZE);
			}
		}
		else
		{
			b3Assert(m_window);
			int width, height;
			m_window->getRenderingResolution(width, height);
		
			static ScreenSpaceFluidRendererGL* fluidRenderer = 0;
			if(!fluidRenderer) fluidRenderer = new ScreenSpaceFluidRendererGL(width, height);
			
			int rendererWidth, rendererHeight;
			fluidRenderer->getWindowResolution(rendererWidth, rendererHeight);
			if(width != rendererWidth || height != rendererHeight)
			{
				fluidRenderer->setWindowResolution(width, height);
				fluidRenderer->setRenderingResolution(width, height);
			}
			
			//Beer's law constants - controls the darkening of the fluid's color based on its thickness
			//For a constant k, (k > 1) == darkens faster; (k < 1) == darkens slower; (k == 0) == disable
			float absorptionR = 0.5, absorptionG = 0.5, absorptionB = 0.5;
			
			const float* projectionMatrix = m_instancingRenderer->getProjectionMatrix();
			const float* modelviewMatrix = m_instancingRenderer->getModelviewMatrix();
			float modelviewProjectionMatrix[16];
			b3Matrix4x4Mul16(projectionMatrix, modelviewMatrix, modelviewProjectionMatrix);
			
			fluidRenderer->render(projectionMatrix, modelviewMatrix, modelviewProjectionMatrix, m_sphFluid->getParticles().m_position, 
									SPHERE_SIZE, COLOR[0], COLOR[1], COLOR[2], absorptionR, absorptionG, absorptionB, true);
		}
	}
	
	virtual void clientMoveAndDisplay()
	{
		RigidBodyGpuData rbData;
		rbData.load(m_data->m_bp, m_data->m_np);
		
		m_solver->stepSimulation(m_sphFluid, rbData);
		
		static int counter = 0;
		if(++counter >= 100)
		{
			counter = 0;
			printf("m_sphFluid->numParticles(): %d \n", m_sphFluid->numParticles());
		}
		
		const bool TEST_PARTICLE_ADD = true;
		if(TEST_PARTICLE_ADD)
		{
			b3FluidEmitter emitter;
			emitter.m_fluid = m_sphFluid;
			
			emitter.m_center = b3MakeVector3(20, 50, 0);
			emitter.m_rotation = b3Quaternion(0, 75, 0);
			emitter.m_speed = 2.5;
			
			b3Scalar spacing = m_sphFluid->getEmitterSpacing() * 1.5f;
			emitter.m_positions.push_back( b3MakeVector3(0, 0, 0) );
			
			emitter.m_positions.push_back( b3MakeVector3(spacing, 0, 0) );
			emitter.m_positions.push_back( b3MakeVector3(-spacing, 0, 0) );
			emitter.m_positions.push_back( b3MakeVector3(0, spacing, 0) );
			emitter.m_positions.push_back( b3MakeVector3(0, -spacing, 0) );
			
			emitter.m_positions.push_back( b3MakeVector3(spacing, spacing, 0) );
			emitter.m_positions.push_back( b3MakeVector3(spacing, -spacing, 0) );
			emitter.m_positions.push_back( b3MakeVector3(-spacing, spacing, 0) );
			emitter.m_positions.push_back( b3MakeVector3(-spacing, -spacing, 0) );
			
			emitter.m_positions.push_back( b3MakeVector3(spacing*2.0, 0, 0) );
			emitter.m_positions.push_back( b3MakeVector3(-spacing*2.0, 0, 0) );
			emitter.m_positions.push_back( b3MakeVector3(0, spacing*2.0, 0) );
			emitter.m_positions.push_back( b3MakeVector3(0, -spacing*2.0, 0) );

			emitter.emit();
		}
		
		const bool TEST_PARTICLE_REMOVE = true;
		if(TEST_PARTICLE_REMOVE)
		{
			const b3Scalar EXTENT(10.0);
			b3Vector3 min = b3MakeVector3(-EXTENT, 0.0, -EXTENT);
			b3Vector3 max = b3MakeVector3(EXTENT, 60.0, EXTENT);
			
			for(int i = 0; i < m_sphFluid->numParticles(); ++i)
			{
				const b3Vector3& position = m_sphFluid->getParticles().m_position[i];
			
				if( min.x < position.x && position.x < max.x
				 && min.y < position.y && position.y < max.y
				 && min.z < position.z && position.z < max.z ) m_sphFluid->markParticleForRemoval(i);
			}
		}
		
		BASE_DEMO_CLASS::clientMoveAndDisplay();
	}

};

#endif