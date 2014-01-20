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

//#include "BulletCollision/CollisionDispatch/b3CollisionObject.h"

#include "b3FluidParticles.h"
#include "b3FluidSphParameters.h"
#include "b3FluidSortingGrid.h"

///BULLET_2_TO_3_PLACEHOLDER
#define CO_USER_TYPE 1111
struct b3CollisionObject 
{ 
	int getInternalType() const { return 0; }
};
struct b3CollisionShape 
{ 
};
///

///Describes a single contact between a b3FluidSph particle and a b3CollisionObject or b3RigidBody.
struct b3FluidSphRigidContact
{
	int m_fluidParticleIndex;
	
	b3Vector3 m_normalOnObject;
	b3Vector3 m_hitPointWorldOnObject;
	b3Scalar m_distance;
};

///Contains all b3FluidSphRigidContact between a b3FluidSph and a b3CollisionObject.
struct b3FluidSphRigidContactGroup
{
	const b3CollisionObject* m_object;
	b3AlignedObjectArray<b3FluidSphRigidContact> m_contacts;
	
	void addContact(const b3FluidSphRigidContact &contact) { m_contacts.push_back(contact); }
	int numContacts() const { return m_contacts.size(); }
};

class b3FluidSphSolver;

///@brief Main fluid class. Coordinates a set of b3FluidParticles with material definition and grid broadphase.
class b3FluidSph //: public b3CollisionObject
{
protected:
	b3FluidSphParametersLocal	m_localParameters;
	
	b3FluidSortingGrid		m_grid;
	
	b3FluidParticles 		m_particles;
	
	b3AlignedObjectArray<int> m_removedFluidIndicies;

	b3AlignedObjectArray<const b3CollisionObject*> m_intersectingRigidAabb;	///<Contains b3CollisionObject/b3RigidBody(not b3Softbody)
	b3AlignedObjectArray<b3FluidSphRigidContactGroup> m_rigidContacts;
	
	//If either override is set, the fluid is passed separately to the solver(b3FluidSph-b3FluidSph interaction is disabled)
	b3FluidSphSolver* m_overrideSolver;
	b3FluidSphParametersGlobal* m_overrideParameters;
	
	
	//These pointers are assigned and managed by an OpenCL solver, if one is being used 
	void* m_fluidDataCL;	//b3FluidSphOpenCL
	void* m_gridDataCL;		//b3FluidSortingGridOpenCL
	
public:
	///@param FG Reference returned by b3FluidRigidDynamicsWorld::getGlobalParameters().
	b3FluidSph(const b3FluidSphParametersGlobal& FG, int maxNumParticles);
	virtual ~b3FluidSph();
	
	int	numParticles() const { return m_particles.size(); }
	int getMaxParticles() const { return m_particles.getMaxParticles(); }
	void setMaxParticles(int maxNumParticles);	///<Removes particles if( maxNumParticles < numParticles() ).
	
	///Returns a particle index; creates a new particle if numParticles() < getMaxParticles(), returns numParticles() otherwise.
	int addParticle(const b3Vector3& position) { return m_particles.addParticle(position); }
	
	///Duplicate indicies are ignored, so a particle may be marked twice without any issues.
	void markParticleForRemoval(int index) { m_removedFluidIndicies.push_back(index); }
	
	void removeAllParticles();
	void removeMarkedParticles();	///<Automatically called during b3FluidRigidDynamicsWorld::stepSimulation(); invalidates grid.
	void insertParticlesIntoGrid(); ///<Automatically called during b3FluidRigidDynamicsWorld::stepSimulation(); updates the grid.
	
	///Avoid placing particles at the same position; particles with same position and velocity will experience identical SPH forces and not seperate.
	void setPosition(int index, const b3Vector3& position) { m_particles.m_pos[index] = position; }
	
	///Sets both velocities; getVelocity() and getEvalVelocity().
	void setVelocity(int index, const b3Vector3& velocity) 
	{
		m_particles.m_vel[index] = velocity;
		m_particles.m_vel_eval[index] = velocity;
	}
	
	///Accumulates a simulation scale force that is applied, and then set to 0 during b3FluidRigidDynamicsWorld::stepSimulation().
	void applyForce(int index, const b3Vector3& force) { m_particles.m_accumulatedForce[index] += force; }
	
	const b3Vector3& getPosition(int index) const { return m_particles.m_pos[index]; }
	const b3Vector3& getVelocity(int index) const { return m_particles.m_vel[index]; }			///<Returns the 'current+(1/2)*timestep' velocity.
	const b3Vector3& getEvalVelocity(int index) const { return m_particles.m_vel_eval[index]; } ///<Returns the current velocity.
	
	void setParticleUserPointer(int index, void* userPointer) { m_particles.m_userPointer[index] = userPointer; }
	void* getParticleUserPointer(int index) const { return m_particles.m_userPointer[index]; }
	//
	const b3FluidSortingGrid& getGrid() const { return m_grid; }
	
	///@param FG Reference returned by b3FluidRigidDynamicsWorld::getGlobalParameters().
	void setGridCellSize(const b3FluidSphParametersGlobal& FG);
	
	//Parameters
	const b3FluidSphParametersLocal& getLocalParameters() const { return m_localParameters; }
	b3FluidSphParametersLocal& getLocalParameters() { return m_localParameters; }
	void setLocalParameters(const b3FluidSphParametersLocal& FP) { m_localParameters = FP; }
	b3Scalar getEmitterSpacing(const b3FluidSphParametersGlobal& FG) const { return m_localParameters.m_particleDist / FG.m_simulationScale; }
	
	///If solver is not 0, then it will be used instead of the solver specified by b3FluidRigidDynamicsWorld::getFluidSolver()
	void setOverrideSolver(b3FluidSphSolver* solver) { m_overrideSolver = solver; }
	b3FluidSphSolver* getOverrideSolver() const { return m_overrideSolver; }
	
	///If parameters is not 0, then it will be used instead of the parameters specified by b3FluidRigidDynamicsWorld::getGlobalParameters()
	void setOverrideParameters(b3FluidSphParametersGlobal* parameters) { m_overrideParameters = parameters; }
	b3FluidSphParametersGlobal* getOverrideParameters() const { return m_overrideParameters; }
	
	//Metablobs	
	b3Scalar getValue(b3Scalar x, b3Scalar y, b3Scalar z) const;
	b3Vector3 getGradient(b3Scalar x, b3Scalar y, b3Scalar z) const;

	const b3FluidParticles& getParticles() const { return m_particles; }
	b3FluidParticles& internalGetParticles() { return m_particles; }
	b3FluidSortingGrid& internalGetGrid() { return m_grid; }
	
	//b3FluidSph-Rigid collisions
	void internalClearRigidContacts()
	{
		m_intersectingRigidAabb.clear();
		m_rigidContacts.clear();
	}
	b3AlignedObjectArray<const b3CollisionObject*>& internalGetIntersectingRigidAabbs() { return m_intersectingRigidAabb; }
	b3AlignedObjectArray<b3FluidSphRigidContactGroup>& internalGetRigidContacts() { return m_rigidContacts; }
	
	const b3AlignedObjectArray<const b3CollisionObject*>& getIntersectingRigidAabbs() const { return m_intersectingRigidAabb; }
	const b3AlignedObjectArray<b3FluidSphRigidContactGroup>& getRigidContacts() const { return m_rigidContacts; }
	
	//b3CollisionObject
	virtual void setCollisionShape(b3CollisionShape *collisionShape) { b3Assert(0); }
	
	virtual void getAabb(b3Vector3& aabbMin, b3Vector3& aabbMax) const
	{
		m_grid.getPointAabb(aabbMin, aabbMax);
		
		b3Scalar radius = m_localParameters.m_particleRadius;
		b3Vector3 extent = b3MakeVector3(radius, radius, radius);
		
		aabbMin -= extent;
		aabbMax += extent;
	}

	static const b3FluidSph* upcast(const b3CollisionObject* colObj)
	{
		// replace later with CO_FLUID_SPH
		return (colObj->getInternalType() == CO_USER_TYPE) ? (const b3FluidSph*)colObj : 0;
	}
	static b3FluidSph* upcast(b3CollisionObject* colObj)
	{
		// replace later with CO_FLUID_SPH
		return (colObj->getInternalType() == CO_USER_TYPE) ? (b3FluidSph*)colObj : 0;
	}
	
	void setFluidDataCL(void* sphData) { m_fluidDataCL = sphData; }
	const void* getFluidDataCL() const { return m_fluidDataCL; }
	void setGridDataCL(void* gridData) { m_gridDataCL = gridData; }
	const void* getGridDataCL() const { return m_gridDataCL; }
};

///@brief Adds particles to a b3FluidSph.
struct b3FluidEmitter
{
	b3Vector3 m_position;

	b3Scalar m_velocity;
	
	b3Scalar m_yaw;
	b3Scalar m_pitch;
	
	b3Scalar m_yawSpread;
	b3Scalar m_pitchSpread;
	
	bool m_useRandomIfAllParticlesAllocated;
	
	b3FluidEmitter()
	{
		m_position.setValue(0,0,0);
		m_yaw = b3Scalar(0);
		m_pitch = b3Scalar(0); 
		m_velocity = b3Scalar(0);
		m_yawSpread = b3Scalar(0);
		m_pitchSpread = b3Scalar(0);
		m_useRandomIfAllParticlesAllocated = true;
	}
	
	void emit(b3FluidSph* fluid, int numParticles, b3Scalar spacing);

	static void addVolume(b3FluidSph* fluid, const b3Vector3& min, const b3Vector3& max, b3Scalar spacing);
};

///@brief Marks particles from a b3FluidSph for removal; see b3FluidSph::removeMarkedParticles().
struct b3FluidAbsorber
{
	b3Vector3 m_min;
	b3Vector3 m_max;
	
	//int m_maxParticlesRemoved;
	//	add velocity limit / max particles removed, etc.?
	
	b3FluidAbsorber()
	{
		m_min.setValue(0,0,0);
		m_max.setValue(0,0,0);
	}
	
	void absorb(b3FluidSph* fluid);
};

#endif


