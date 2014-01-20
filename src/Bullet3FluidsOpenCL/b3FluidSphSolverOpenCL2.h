/*
BulletFluids 
Copyright (c) 2013 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef B3_FLUID_SPH_SOLVER_OPENCL2_H
#define B3_FLUID_SPH_SOLVER_OPENCL2_H

#include "Bullet3Common/b3AlignedObjectArray.h"

#include "Bullet3Fluids/Sph/b3FluidSphSolver.h"

#include "Bullet3FluidsOpenCL/b3FluidSphRigidInteractorCL.h"

#include "b3FluidSphOpenCL.h"
#include "b3FluidHashGridOpenCL.h"

struct b3FluidSphParametersGlobal;
class b3FluidSph;

///Uses an infinite hashed grid with collisions, as opposed to a statically sized(1024^3) grid with no collisions as in b3FluidSphSolverOpenCL
class b3FluidSphSolverOpenCL2 : public b3FluidSphSolver
{
	cl_context m_context;
	cl_command_queue m_commandQueue;
	
	cl_program m_fluidsProgram;
	
	cl_kernel m_sphComputePressureModuloKernel;
	cl_kernel m_sphComputeForceModuloKernel;
	
	cl_kernel m_applyForcesKernel;
	cl_kernel m_collideAabbImpulseKernel;
	cl_kernel m_integratePositionKernel;
	
	b3OpenCLArray<b3FluidSphParametersGlobal> m_globalFluidParams;
	
	b3AlignedObjectArray<b3FluidSphOpenCL*> m_fluidData;
	b3AlignedObjectArray<b3FluidHashGridOpenCL*> m_gridData;
	
	b3FluidHashGridOpenCLProgram m_hashGridProgram;
	b3FluidSphRigidInteractorCL m_fluidRigidInteractor;
	
	b3AlignedObjectArray<b3Vector3> m_tempSphForce;
	
public:	
	b3FluidSphSolverOpenCL2(cl_context context, cl_device_id device, cl_command_queue queue);
	virtual ~b3FluidSphSolverOpenCL2();
	
	//	remove/rename
	virtual void updateGridAndCalculateSphForces(const b3FluidSphParametersGlobal& FG, b3FluidSph** fluids, int numFluids) { b3Assert(0); }
	
	virtual void stepSimulation(const b3FluidSphParametersGlobal& FG, b3FluidSph** fluids, int numFluids, RigidBodyGpuData& rbData);
	
private:
	void sphComputePressureModulo(int numFluidParticles, b3FluidHashGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize);
	void sphComputeForceModulo(int numFluidParticles, b3FluidHashGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize);
};

#endif


