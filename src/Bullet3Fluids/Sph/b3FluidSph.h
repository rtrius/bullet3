/*
BulletFluids 
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
//Portions of this file based on FLUIDS v.2 - SPH Fluid Simulator for CPU and GPU
//Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

#ifndef B3_FLUID_SPH_H
#define B3_FLUID_SPH_H

#include "Bullet3Common/b3Quaternion.h"

#include "b3FluidParticles.h"
#include "b3FluidSphParameters.h"
#include "b3FluidSortingGrid.h"
#include "b3FluidSphUpdatePacket.h"

///@brief Specifies b3FluidSph data copied between the CPU and GPU every frame.
///@remarks
/// - Members beginning with m_sync* determine a copy from GPU to CPU at the end of a frame
/// - Members beginning with m_write* determine a copy from CPU to GPU at the beginning of a frame
///@par
///The transfer of particle data from CPU to GPU and back is a fairly expensive operation,
///so as many flags as possible should be set to false.
struct b3FluidSphSyncronizationFlags
{
	bool m_syncPosition;			///<b3FluidParticles.m_position ( b3FluidSph::getParticles() )
	bool m_syncVelocity;			///<b3FluidParticles.m_velocity ( b3FluidSph::getParticles() )
	bool m_syncVelocityEval;		///<b3FluidParticles.m_velocityEval ( b3FluidSph::getParticles() )
	//bool m_syncUserPointer;			///<b3FluidParticles.m_userPointer ( b3FluidSph::getParticles() )
	
	//bool m_syncGridState;			///<If true, the state of the b3FluidSortingGrid is copied back to CPU every frame
	//boom m_syncRigidContacts;		///<If true, data from collision with rigid bodies is copied back to CPU every frame
	
	bool m_writeForces;			///<If false, applied forces are not written to GPU( b3FluidSph::applyForce() will have no effect )
	
	b3FluidSphSyncronizationFlags()
	{
		m_syncPosition = true;
		m_syncVelocity = false;
		m_syncVelocityEval = false;
		//m_syncUserPointer = false;
		
		//m_syncGridState = false;
		//m_syncRigidContacts = false;
		
		m_writeForces = false;
	}
};

class b3FluidSphSolver;
class b3FluidSphTypedData;

///@brief Main class for particle fluids.
///@remarks b3FluidSph is basically a manager class for several components used to simulate a SPH fluid.
///These components are:
/// - b3FluidSphParameters, which defines the material of the fluid
/// - b3FluidParticles, which manages the per-particle properities such as position, velocity, and accumulated force
/// - b3FluidSortingGrid, a uniform grid broadphase for accelerating the SPH force calculation
/// - b3FluidSphUpdatePacket, which is used to incrementally add particles, remove particles, and set per-particle properities 
/// - b3FluidSphSyncronizationFlags, which is used to specify data copied back to and from the CPU every frame(can be ignored if using a CPU solver)
///@remarks Each b3FluidSph cannot interact with another b3FluidSph.
///By defining another b3FluidSphSolver, it is also possible to use this as a generic particle system.
class b3FluidSph
{
protected:
	b3FluidSphParameters m_parameters;
	b3FluidParticles m_particles;
	b3FluidSortingGrid m_grid;
	
	b3FluidSphUpdatePacket m_updates;
	b3FluidSphSyncronizationFlags m_gpuSyncFlags;
	
	b3FluidSphSolver* m_solver;
	
	//This pointer is assigned and managed by a CPU b3FluidSphSolver, if it is being used
	b3FluidSphTypedData* m_solverData;
	
	//These pointers are assigned and managed by an OpenCL solver, if one is being used 
	b3FluidSphTypedData* m_solverDataGpu;	//e.g. b3FluidSphOpenCL
	b3FluidSphTypedData* m_gridDataGpu;		//e.g. b3FluidSortingGridOpenCL
	
	bool m_copyParticleDataToGpu;	
	
public:
	b3FluidSph(b3FluidSphSolver* solver, int maxNumParticles);
	virtual ~b3FluidSph();
	
	int	numParticles() const { return m_particles.size(); }
	int getMaxParticles() const { return m_particles.getMaxParticles(); }
	
	///Contains fluid material definition(e.g. viscosity) and rigid body interaction parameters 
	const b3FluidSphParameters& getParameters() const { return m_parameters; }
	b3FluidSphParameters& getParameters() { return m_parameters; }
	void setParameters(const b3FluidSphParameters& FP) { m_parameters = FP; }
	
	///If using a GPU(OpenCL) solver, the sync flags determine what data is copied back to the CPU every frame 
	const b3FluidSphSyncronizationFlags& getGpuSyncFlags() const { return m_gpuSyncFlags; }
	b3FluidSphSyncronizationFlags& getGpuSyncFlags() { return m_gpuSyncFlags; }
	void setGpuSyncFlags(const b3FluidSphSyncronizationFlags& gpuSyncFlags) { m_gpuSyncFlags = gpuSyncFlags; }
	
	///The solver determines how the particles interact with each other
	void setSolver(b3FluidSphSolver* solver) { m_solver = solver; }
	b3FluidSphSolver* getSolver() const { return m_solver; }
	
	b3Scalar getEmitterSpacing() const { return m_parameters.m_particleDist / m_parameters.m_simulationScale; }
	
	///Use to access particle positions and velocities; if using a GPU solver make sure that getGpuSyncFlags() is set correctly.
	const b3FluidParticles& getParticles() const { return m_particles; }
	const b3FluidSortingGrid& getGrid() const { return m_grid; }
	const b3FluidSphUpdatePacket& getUpdates() const { return m_updates; }
	
	///Does not distinguish whether the data on CPU or GPU is current; use with caution. 
	b3FluidParticles& internalGetParticles() { return m_particles; }
	b3FluidSortingGrid& internalGetGrid() { return m_grid; }
	b3FluidSphUpdatePacket& internalGetUpdates() { return m_updates; }
	
	///If true, CPU position, velocity, accumulated force is copied to GPU on next simulation step.
	///This is automatically set to false after the data is copied.
	bool needsWriteStateToGpu() const { return m_copyParticleDataToGpu; }
	void shouldWriteStateToGpu(bool transferCpuToGpuNextFrame) { m_copyParticleDataToGpu = transferCpuToGpuNextFrame; }
	
	//These functions assume that the data on CPU is current; if a GPU solver is being used they should only be called during initialization
	//and even then the incremental update functions such as addParticleCached() should be used.
	///@name Calling any of these functions will overwrite the GPU data with the CPU data on the next frame( shouldWriteStateToGpu(true) ).
	///@{
		void setMaxParticles(int maxNumParticles);	///<Removes particles if( maxNumParticles < numParticles() ).
		void removeAllParticles();		///<Also clears any updates or added particles.
		
		void applyUpdates(); 			///<Automatically called during b3FluidRigidDynamicsWorld::stepSimulation(); invalidates grid.
		void insertParticlesIntoGrid(); ///<Automatically called during b3FluidRigidDynamicsWorld::stepSimulation(); updates the grid.
	///@}

	//If using a GPU solver, adding and removing particles will cause the CPU particle arrays to be resized, 
	//but the data is not copied over unless getGpuSyncFlags() is set.
	///@name Cached particle update. Note that these changes are not applied until applyUpdates(), or (for GPU solvers) stepSimulation() is called.
	///@{
		///Returns a particle index; creates a new particle if numParticles() < getMaxParticles(), returns numParticles() otherwise.
		///The particle indicies change during each internal simulation step, so the returned index should be used only for initialization.
		int addParticleCached( const b3Vector3& position, const b3Vector3& velocity = b3MakeVector3(0,0,0) ) 
		{
			return m_updates.addParticle( getMaxParticles(), numParticles(), position, velocity); 
		}
		
		///Do not use the same index twice; will result in thread collisions.
		///Avoid placing particles at the same position; particles with same position and velocity,
		///as they will experience identical SPH forces and not seperate.
		void setPositionCached(int index, const b3Vector3& position) { m_updates.setPosition(index, position); }
		
		///Do not use the same index twice; will result in thread collisions.
		///Velocity is at simulation scale; make sure to multiply by getParameters().m_simulationScale before applying if at world scale.
		void setVelocityCached(int index, const b3Vector3& velocity) { m_updates.setVelocity(index, velocity); }
		
		///Do not use the same index twice; will result in thread collisions.
		//void setParticleUserPointer(int index, void* userPointer) { m_particles.m_userPointer[index] = userPointer; }
		
		///Duplicate indicies are ignored, so a particle may be marked twice without any issues.
		void markParticleForRemoval(int index) { m_updates.markParticleForRemoval(index); }
		
		///Accumulates a simulation scale force that is applied, and then set to 0 during b3FluidRigidDynamicsWorld::stepSimulation().
		void applyForce(int index, const b3Vector3& force) { m_particles.m_accumulatedForce[index] += force; }
	///@}
	
	///@name Solver specific data managed by the solver. Can be accessed but should not be modified.
	///@{
		void setSolverDataCpu(b3FluidSphTypedData* solverData) { m_solverData = solverData; }
		b3FluidSphTypedData* getSolverDataCpu() const { return m_solverData; }
		
		void setSolverDataGpu(b3FluidSphTypedData* sphData) { m_solverDataGpu = sphData; }
		b3FluidSphTypedData* getSolverDataGpu() const { return m_solverDataGpu; }
		void setGridDataGpu(b3FluidSphTypedData* gridData) { m_gridDataGpu = gridData; }
		b3FluidSphTypedData* getGridDataGpu() const { return m_gridDataGpu; }
	///@}
};

///@brief Adds particles to a b3FluidSph.
class b3FluidEmitter
{
public:
	b3FluidSph* m_fluid;
	
	///If this is nonzero, m_position is also offset by the btCollisionObjects' position and orientation
	///must be either a btCollisionObject or btRigidBody
	//btCollisionObject* m_attachTo;	
	
	///Each position in this array corresponds to a particle created when emit() is called
	b3AlignedObjectArray<b3Vector3> m_positions;
	
	///This transform is applied to m_positions and m_direction; it defines the world transform of the emitter
	b3Vector3 m_center;
	b3Quaternion m_rotation;
	
	///Must have length() == 1.0; by default points at -z; it is generally better to change m_rotation instead of this
	b3Vector3 m_direction;
	b3Scalar m_speed;
	
	///Emit particles if true
	bool m_active;
	
	///If the maximum number of particles in a btFluidSph is reached,
	///reset and use a randomly selected existing particle
	//bool m_useRandomIfAllParticlesAllocated;
	
	///Contains the indicies of particles emitted on the last call to emit()
	//btAlignedObjectArray<int> m_particleIndicies;
	
	b3FluidEmitter()
	{
		m_fluid = 0;
		//m_attachTo = 0;
		
		m_center = b3MakeVector3(0, 0, 0);
		m_rotation = b3Quaternion(0, 0, 0, 1);
		m_direction = b3MakeVector3(0, 0, -1);
		m_speed = b3Scalar(1.0);
		
		m_active = true;
		//m_useRandomIfAllParticlesAllocated = true;
	}
	
	///This does not need to be called if the emitter is attached using btFluidRigidDynamicsWorld::addSphEmitter()
	void emit();


	static void addVolume(b3FluidSph* fluid, const b3Vector3& min, const b3Vector3& max, b3Scalar spacing, int maxCreated);
};

#endif


