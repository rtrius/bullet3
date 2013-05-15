/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef B3_FLUID_SPH_SOLVER_OPENCL_H
#define B3_FLUID_SPH_SOLVER_OPENCL_H

#include "Bullet3Common/b3AlignedObjectArray.h"

#include "BulletFluids/Sph/b3FluidSphSolver.h"

#include "BulletFluidsOpenCL/b3FluidSphRigidInteractorCL.h"

#include "b3FluidSphOpenCL.h"
#include "b3FluidSortingGridOpenCL.h"

struct b3FluidSphParametersGlobal;
class b3FluidSph;

///@brief Solver that uses the GPU to accelerate SPH force calculation.
///@remarks
///Does not implement fluid-fluid interactions.
class b3FluidSphSolverOpenCL : public b3FluidSphSolver
{
	cl_context m_context;
	cl_command_queue m_commandQueue;
	
	cl_program m_fluidsProgram;
	cl_kernel m_kernel_findNeighborCellsPerCell;
	cl_kernel m_kernel_findGridCellIndexPerParticle;
	cl_kernel m_kernel_sphComputePressure;
	cl_kernel m_kernel_sphComputeForce;

	cl_kernel m_applyForcesKernel;
	cl_kernel m_collideAabbImpulseKernel;
	cl_kernel m_integratePositionKernel;
	
	b3OpenCLArray<b3FluidSphParametersGlobal> m_globalFluidParams;
	
	b3AlignedObjectArray<b3FluidSphOpenCL*> m_fluidData;
	b3AlignedObjectArray<b3FluidSortingGridOpenCL*> m_gridData;
	
	b3FluidSortingGridOpenCLProgram m_sortingGridProgram;
	b3FluidSphRigidInteractorCL m_fluidRigidInteractor;
	
	b3AlignedObjectArray<b3Vector3> m_tempSphForce;
	
public:	
	b3FluidSphSolverOpenCL(cl_context context, cl_device_id device, cl_command_queue queue);
	virtual ~b3FluidSphSolverOpenCL();
	
	//	remove/rename
	virtual void updateGridAndCalculateSphForces(const b3FluidSphParametersGlobal& FG, b3FluidSph** fluids, int numFluids) { b3Assert(0); }
	
	virtual void stepSimulation(const b3FluidSphParametersGlobal& FG, b3FluidSph** fluids, int numFluids, RigidBodyGpuData& rbData);
	
private:
	void findNeighborCells(int numActiveGridCells, int numFluidParticles, b3FluidSortingGridOpenCL* gridData, b3FluidSphOpenCL* fluidData);
	void sphComputePressure(int numFluidParticles, b3FluidSortingGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize);
	void sphComputeForce(int numFluidParticles, b3FluidSortingGridOpenCL* gridData, b3FluidSphOpenCL* fluidData, b3Scalar cellSize);
};

#endif


