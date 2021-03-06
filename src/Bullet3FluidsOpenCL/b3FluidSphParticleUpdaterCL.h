
#ifndef B3_FLUID_SPH_PARTICLE_UPDATER_CL_H
#define B3_FLUID_SPH_PARTICLE_UPDATER_CL_H

#include "Bullet3Common/b3Vector3.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3FillCL.h"

class b3FluidSph;
class b3FluidSphOpenCL;

///Performs incremential particle add, update, and remove on the GPU
class b3FluidSphParticleUpdaterCL
{
protected:
	cl_context m_context;
	cl_command_queue m_commandQueue;
	
	cl_program m_fluidsProgram;
	cl_kernel m_applyParticleUpdatesKernel;
	cl_kernel m_swapRemovedParticlesKernel;
	
	b3FillCL m_fill;
	
	//
	b3OpenCLArray<b3Vector3> m_createdPosition;
	b3OpenCLArray<b3Vector3> m_createdVelocity;
	
	//
	b3OpenCLArray<int> m_updatedPositionIndices;
	b3OpenCLArray<b3Vector3> m_updatedPosition;
	b3OpenCLArray<int> m_updatedVelocityIndices;
	b3OpenCLArray<b3Vector3> m_updatedVelocity;
	
	//
	b3OpenCLArray<int> m_removeSwapSource;
	b3OpenCLArray<int> m_removeSwapTarget;
	
public:
	b3FluidSphParticleUpdaterCL(cl_context context, cl_device_id device, cl_command_queue queue);
	virtual ~b3FluidSphParticleUpdaterCL();

	//GPU version of b3FluidSph::applyUpdates()
	void createParticlesApplyUpdatesAndRemoveParticles(b3FluidSph* fluid, b3FluidSphOpenCL* fluidDataCL);
};

#endif
